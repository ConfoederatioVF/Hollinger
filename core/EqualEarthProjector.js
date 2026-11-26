"use strict";

/**
 * Analytic Equal Earth projection (Savric, Patterson & Jenny, 2018) with precomputed
 * inverse pixel lookup tables for rapid 4K equirectangular-to-Equal-Earth rendering.
 */
class EqualEarthProjector {
  /**
   * @param {Object} arg0_canvas_config
   * @param {Object} [arg1_src_grid]
   */
  constructor (arg0_canvas_config, arg1_src_grid) {
    //Convert from parameters
    let canvas_cfg = arg0_canvas_config;
    let src_grid = arg1_src_grid ? arg1_src_grid : { height: 2160, width: 4320 };

    //Declare local instance variables
    let A1 = 1.340264;
    let A2 = -0.081106;
    let A3 = 0.000893;
    let A4 = 0.003796;
    let M = Math.sqrt(3.0)/2.0;

    //Compute exact unit-sphere Equal Earth bounding half-extents:
    //  x_max occurs at phi = 0, lambda = pi -> x_max = (2*sqrt(3)*pi) / (3*A1)
    //  y_max occurs at phi = pi/2 -> theta_max = asin(M) = pi/3
    let theta_max = Math.asin(M);
    let t2 = theta_max*theta_max;
    let t6 = t2*t2*t2;
    let x_max = (2.0*Math.sqrt(3.0)*Math.PI)/(3.0*A1);
    let y_max = theta_max*(A1 + A2*t2 + t6*(A3 + A4*t2));

    let avail_h = canvas_cfg.height - canvas_cfg.margin_top - canvas_cfg.margin_bottom;
    let avail_w = canvas_cfg.width - canvas_cfg.margin_left - canvas_cfg.margin_right;
    let scale = Math.min(avail_w/(2.0*x_max), avail_h/(2.0*y_max));
    let cx = canvas_cfg.margin_left + avail_w*0.5;
    let cy = canvas_cfg.margin_top + avail_h*0.5;

    this.bilinear = Boolean(canvas_cfg.bilinear_sampling);
    this.central_meridian_rad = ((canvas_cfg.central_meridian_deg || 0.0)*Math.PI)/180.0;
    this.height = canvas_cfg.height;
    this.map_top_px = cy - y_max*scale;
    this.src_height = src_grid.height;
    this.src_width = src_grid.width;
    this.width = canvas_cfg.width;

    //Precompute inverse LUT for every output pixel
    let total_px = this.width*this.height;
    this.mask = new Uint8Array(total_px);
    this.src_col = new Float32Array(total_px);
    this.src_row = new Float32Array(total_px);

    for (let py = 0; py < this.height; py++) {
      let y_proj = (cy - (py + 0.5))/scale;
      let row_offset = py*this.width;

      if (Math.abs(y_proj) > y_max) continue;

      //Newton-Raphson solve for theta (depends only on py!)
      let theta = y_proj/A1;
      for (let i = 0; i < 12; i++) {
        let th2 = theta*theta;
        let th6 = th2*th2*th2;
        let fn = theta*(A1 + A2*th2 + th6*(A3 + A4*th2)) - y_proj;
        let fd = A1 + 3.0*A2*th2 + th6*(7.0*A3 + 9.0*A4*th2);
        let delta = fn/fd;
        theta -= delta;
        if (Math.abs(delta) < 1e-12) break;
      }

      let th2 = theta*theta;
      let th6 = th2*th2*th2;
      let denom = M*(A1 + 3.0*A2*th2 + th6*(7.0*A3 + 9.0*A4*th2));
      let cos_theta = Math.cos(theta);
      let sin_phi = Math.min(1.0, Math.max(-1.0, Math.sin(theta)/M));
      let phi = Math.asin(sin_phi);
      let lat_deg = (phi*180.0)/Math.PI;
      let src_r = ((90.0 - lat_deg)/180.0)*this.src_height - 0.5;
      let clamped_r = Math.min(this.src_height - 1.0001, Math.max(0.0, src_r));

      for (let px = 0; px < this.width; px++) {
        let x_proj = ((px + 0.5) - cx)/scale;
        let lambda = (x_proj*denom)/cos_theta;

        if (Math.abs(lambda) > Math.PI) continue;

        let lon_rad = lambda + this.central_meridian_rad;
        if (lon_rad > Math.PI) lon_rad -= 2.0*Math.PI;
        if (lon_rad < -Math.PI) lon_rad += 2.0*Math.PI;

        let lon_deg = (lon_rad*180.0)/Math.PI;
        let src_c = ((lon_deg + 180.0)/360.0)*this.src_width - 0.5;
        if (src_c < 0.0) src_c += this.src_width;
        if (src_c >= this.src_width) src_c -= this.src_width;

        let idx = row_offset + px;
        this.mask[idx] = 1;
        this.src_row[idx] = clamped_r;
        this.src_col[idx] = src_c;
      }
    }

    //Precompute 1-2px utilitarian globe border ring outside the oval
    this.border = new Uint8Array(total_px);
    for (let py = 1; py < this.height - 1; py++) {
      let row_offset = py*this.width;
      for (let px = 1; px < this.width - 1; px++) {
        let idx = row_offset + px;
        if (this.mask[idx] === 0) {
          if (
            this.mask[idx - 1] === 1 ||
            this.mask[idx + 1] === 1 ||
            this.mask[idx - this.width] === 1 ||
            this.mask[idx + this.width] === 1
          )
            this.border[idx] = 1;
        }
      }
    }

    //Precompute 4x4 subpixel anti-aliased edge alpha and nearest-valid-source LUT fallback for borderless rendering
    this.edge_alpha = new Uint8Array(total_px);
    for (let i = 0; i < total_px; i++) {
      if (this.mask[i] === 1) this.edge_alpha[i] = 255;
    }

    let sub_offsets = [0.125, 0.375, 0.625, 0.875];
    for (let py = 1; py < this.height - 1; py++) {
      let row_offset = py*this.width;
      for (let px = 1; px < this.width - 1; px++) {
        let idx = row_offset + px;
        let m = this.mask[idx];
        let is_boundary = (
          this.mask[idx - 1] !== m ||
          this.mask[idx + 1] !== m ||
          this.mask[idx - this.width] !== m ||
          this.mask[idx + this.width] !== m
        );
        if (!is_boundary) continue;

        let inside_count = 0;
        for (let sy_i = 0; sy_i < 4; sy_i++) {
          let y_proj = (cy - (py + sub_offsets[sy_i]))/scale;
          if (Math.abs(y_proj) > y_max) continue;

          let theta = y_proj/A1;
          for (let it = 0; it < 12; it++) {
            let th2 = theta*theta;
            let th6 = th2*th2*th2;
            let fn = theta*(A1 + A2*th2 + th6*(A3 + A4*th2)) - y_proj;
            let fd = A1 + 3.0*A2*th2 + th6*(7.0*A3 + 9.0*A4*th2);
            let delta = fn/fd;
            theta -= delta;
            if (Math.abs(delta) < 1e-12) break;
          }

          let th2 = theta*theta;
          let th6 = th2*th2*th2;
          let denom = M*(A1 + 3.0*A2*th2 + th6*(7.0*A3 + 9.0*A4*th2));
          let cos_theta = Math.cos(theta);

          for (let sx_i = 0; sx_i < 4; sx_i++) {
            let x_proj = ((px + sub_offsets[sx_i]) - cx)/scale;
            let lambda = (x_proj*denom)/cos_theta;
            if (Math.abs(lambda) <= Math.PI) inside_count++;
          }
        }

        this.edge_alpha[idx] = Math.round((inside_count/16.0)*255.0);

        if (m === 0 && this.edge_alpha[idx] > 0) {
          let neighbours = [idx - 1, idx + 1, idx - this.width, idx + this.width];
          for (let n_i = 0; n_i < neighbours.length; n_i++) {
            let n_idx = neighbours[n_i];
            if (this.mask[n_idx] === 1) {
              this.src_row[idx] = this.src_row[n_idx];
              this.src_col[idx] = this.src_col[n_idx];
              break;
            }
          }
        }
      }
    }
  }

  /**
   * Projects a 4320x2160 signed float32 km2 grid onto a 3840x2160 RGB24 buffer using the QBAM colour palette.
   *
   * @param {Float32Array} arg0_land_km2 - Signed 4320x2160 km2 grid (>0 land, <0 glacier, 0 sea).
   * @param {Float64Array|Float32Array} arg1_cell_area_by_row
   * @param {Object} arg2_palette_rgb - { bg, border, glacier, land, sea, land_blend, glacier_blend }
   * @param {Uint8Array} [arg3_out_rgb] - Reusable output RGB24 buffer (width*height*3).
   *
   * @returns {Uint8Array}
   */
  projectToRgb (arg0_land_km2, arg1_cell_area_by_row, arg2_palette_rgb, arg3_out_rgb) {
    //Convert from parameters
    let land_km2 = arg0_land_km2;
    let cell_area_by_row = arg1_cell_area_by_row;
    let pal = arg2_palette_rgb;
    let out_rgb = arg3_out_rgb ? arg3_out_rgb : new Uint8Array(this.width*this.height*3);

    //Declare local instance variables
    let bg_r = pal.bg[0], bg_g = pal.bg[1], bg_b = pal.bg[2];
    let bd_r = pal.border ? pal.border[0] : bg_r;
    let bd_g = pal.border ? pal.border[1] : bg_g;
    let bd_b = pal.border ? pal.border[2] : bg_b;
    let border = this.border;
    let gl_blend = Boolean(pal.glacier_blend);
    let gl_r = pal.glacier[0], gl_g = pal.glacier[1], gl_b = pal.glacier[2];
    let has_border = Boolean(pal.border);
    let ln_blend = Boolean(pal.land_blend);
    let ln_r = pal.land[0], ln_g = pal.land[1], ln_b = pal.land[2];
    let mask = this.mask;
    let sea_r = pal.sea[0], sea_g = pal.sea[1], sea_b = pal.sea[2];
    let src_col = this.src_col;
    let src_row = this.src_row;
    let src_w = this.src_width;
    let total_px = this.width*this.height;

    //Pre-normalise source grid into RGB24 once per frame for crisp, sign-safe bilinear projection!
    //(Sign-safe: avoids blending positive land and negative glacier through zero!)
    let src_total = this.src_width*this.src_height;
    if (!this._src_rgb_cache || this._src_rgb_cache.length !== src_total*3)
      this._src_rgb_cache = new Uint8Array(src_total*3);
    let src_rgb = this._src_rgb_cache;

    for (let y = 0; y < this.src_height; y++) {
      let inv_area = 1.0/cell_area_by_row[y];
      let row_off = y*src_w;
      let y_up = (y > 0) ? (y - 1)*src_w : row_off;
      let y_dn = (y < this.src_height - 1) ? (y + 1)*src_w : row_off;

      for (let x = 0; x < src_w; x++) {
        let idx = row_off + x;
        let val = land_km2[idx];
        let o3 = idx*3;

        if (val > 0.0) {
          let f = ln_blend ? Math.min(1.0, val*inv_area) : 1.0;
          if (f >= 0.995) {
            src_rgb[o3] = ln_r;
            src_rgb[o3 + 1] = ln_g;
            src_rgb[o3 + 2] = ln_b;
          } else {
            //Determine whether fractional complement is glacier (inland ice margin / mountain glacier) or water
            let x_l = (x === 0) ? (src_w - 1) : (x - 1);
            let x_r = (x === src_w - 1) ? 0 : (x + 1);
            let n0 = land_km2[row_off + x_l];
            let n1 = land_km2[row_off + x_r];
            let n2 = land_km2[y_up + x];
            let n3 = land_km2[y_dn + x];
            let adj_sea = (n0 === 0.0 || n1 === 0.0 || n2 === 0.0 || n3 === 0.0);
            let adj_ice = (n0 < 0.0 || n1 < 0.0 || n2 < 0.0 || n3 < 0.0);

            if (adj_ice && !adj_sea) {
              //Inland land-glacier margin: blend between glacier greyish-white and land lime (never blue sea!)
              src_rgb[o3] = (gl_r + (ln_r - gl_r)*f + 0.5) | 0;
              src_rgb[o3 + 1] = (gl_g + (ln_g - gl_g)*f + 0.5) | 0;
              src_rgb[o3 + 2] = (gl_b + (ln_b - gl_b)*f + 0.5) | 0;
            } else {
              let f_vis = Math.sqrt(f);
              src_rgb[o3] = (sea_r + (ln_r - sea_r)*f_vis + 0.5) | 0;
              src_rgb[o3 + 1] = (sea_g + (ln_g - sea_g)*f_vis + 0.5) | 0;
              src_rgb[o3 + 2] = (sea_b + (ln_b - sea_b)*f_vis + 0.5) | 0;
            }
          }
        } else if (val < 0.0) {
          let f = gl_blend ? Math.min(1.0, (-val)*inv_area) : 1.0;
          if (f >= 0.995) {
            src_rgb[o3] = gl_r;
            src_rgb[o3 + 1] = gl_g;
            src_rgb[o3 + 2] = gl_b;
          } else {
            let x_l = (x === 0) ? (src_w - 1) : (x - 1);
            let x_r = (x === src_w - 1) ? 0 : (x + 1);
            let n0 = land_km2[row_off + x_l];
            let n1 = land_km2[row_off + x_r];
            let n2 = land_km2[y_up + x];
            let n3 = land_km2[y_dn + x];
            let adj_sea = (n0 === 0.0 || n1 === 0.0 || n2 === 0.0 || n3 === 0.0);
            let adj_land = (n0 > 0.0 || n1 > 0.0 || n2 > 0.0 || n3 > 0.0);

            if (adj_land && !adj_sea) {
              src_rgb[o3] = (ln_r + (gl_r - ln_r)*f + 0.5) | 0;
              src_rgb[o3 + 1] = (ln_g + (gl_g - ln_g)*f + 0.5) | 0;
              src_rgb[o3 + 2] = (ln_b + (gl_b - ln_b)*f + 0.5) | 0;
            } else {
              let f_vis = Math.sqrt(f);
              src_rgb[o3] = (sea_r + (gl_r - sea_r)*f_vis + 0.5) | 0;
              src_rgb[o3 + 1] = (sea_g + (gl_g - sea_g)*f_vis + 0.5) | 0;
              src_rgb[o3 + 2] = (sea_b + (gl_b - sea_b)*f_vis + 0.5) | 0;
            }
          }
        } else {
          src_rgb[o3] = sea_r;
          src_rgb[o3 + 1] = sea_g;
          src_rgb[o3 + 2] = sea_b;
        }
      }
    }

    if (this.bilinear) {
      for (let i = 0; i < total_px; i++) {
        let o3 = i*3;
        if (mask[i] === 0) {
          if (has_border && border[i] === 1) {
            out_rgb[o3] = bd_r;
            out_rgb[o3 + 1] = bd_g;
            out_rgb[o3 + 2] = bd_b;
          } else {
            out_rgb[o3] = bg_r;
            out_rgb[o3 + 1] = bg_g;
            out_rgb[o3 + 2] = bg_b;
          }
          continue;
        }

        let r_f = src_row[i];
        let c_f = src_col[i];
        let r0 = r_f | 0;
        let c0 = c_f | 0;
        let r1 = (r0 + 1 < this.src_height) ? (r0 + 1) : r0;
        let c1 = (c0 + 1 < src_w) ? (c0 + 1) : 0;
        let dr = r_f - r0;
        let dc = c_f - c0;
        let w00 = (1.0 - dr)*(1.0 - dc);
        let w01 = (1.0 - dr)*dc;
        let w10 = dr*(1.0 - dc);
        let w11 = dr*dc;

        let p00 = (r0*src_w + c0)*3;
        let p01 = (r0*src_w + c1)*3;
        let p10 = (r1*src_w + c0)*3;
        let p11 = (r1*src_w + c1)*3;

        out_rgb[o3] = (src_rgb[p00]*w00 + src_rgb[p01]*w01 + src_rgb[p10]*w10 + src_rgb[p11]*w11 + 0.5) | 0;
        out_rgb[o3 + 1] = (src_rgb[p00 + 1]*w00 + src_rgb[p01 + 1]*w01 + src_rgb[p10 + 1]*w10 + src_rgb[p11 + 1]*w11 + 0.5) | 0;
        out_rgb[o3 + 2] = (src_rgb[p00 + 2]*w00 + src_rgb[p01 + 2]*w01 + src_rgb[p10 + 2]*w10 + src_rgb[p11 + 2]*w11 + 0.5) | 0;
      }
    } else {
      for (let i = 0; i < total_px; i++) {
        let o3 = i*3;
        if (mask[i] === 0) {
          if (has_border && border[i] === 1) {
            out_rgb[o3] = bd_r;
            out_rgb[o3 + 1] = bd_g;
            out_rgb[o3 + 2] = bd_b;
          } else {
            out_rgb[o3] = bg_r;
            out_rgb[o3 + 1] = bg_g;
            out_rgb[o3 + 2] = bg_b;
          }
          continue;
        }

        let r0 = Math.round(src_row[i]);
        let c0 = Math.round(src_col[i]) % src_w;
        let p00 = (r0*src_w + c0)*3;
        out_rgb[o3] = src_rgb[p00];
        out_rgb[o3 + 1] = src_rgb[p00 + 1];
        out_rgb[o3 + 2] = src_rgb[p00 + 2];
      }
    }

    //Return statement
    return out_rgb;
  }


  /**
   * Projects a categorical integer grid onto an RGB frame using a flat LUT, optionally applying
   * class-preserving bilinear coastal antialiasing from physical land/glaciated fractions
   * with borderless 4x4 subpixel anti-aliased globe edges.
   *
   * @param {Uint8Array|Uint16Array|Uint32Array} arg0_class_grid - Source 2D class ID grid (row-major).
   * @param {Uint8Array} arg1_lut_rgb - Flat RGB lookup table (length >= (max_id + 1) * 3).
   * @param {Array<number>} [arg2_bg_rgb] - Canvas background RGB [r,g,b].
   * @param {Object} [arg3_options] - Optional { land_km2, cell_area_by_row, ice_class_id }.
   * @returns {Uint8Array} Projected RGB frame buffer (width * height * 3).
   */
  projectCategoricalToRgb (arg0_class_grid, arg1_lut_rgb, arg2_bg_rgb, arg3_options) {
    //Convert from parameters
    let class_grid = arg0_class_grid;
    let lut_rgb = arg1_lut_rgb;
    let bg = arg2_bg_rgb || [0, 0, 0];
    let options = arg3_options ? arg3_options : {};

    //Declare local instance variables
    let bg_r = bg[0];
    let bg_g = bg[1];
    let bg_b = bg[2];
    let cell_area_by_row = options.cell_area_by_row || null;
    let edge_alpha = this.edge_alpha;
    let land_km2 = options.land_km2 || null;
    let lut_len = lut_rgb.length;
    let n_pix = this.width*this.height;
    let out_rgb = new Uint8Array(n_pix*3);
    let src_col = this.src_col;
    let src_h = this.src_height;
    let src_row = this.src_row;
    let src_w = this.src_width;

    if (land_km2 && cell_area_by_row) {
      let ice_id = (options.ice_class_id !== undefined) ? options.ice_class_id : 0;
      let i3 = ice_id*3;
      let ice_r = (i3 + 2 < lut_len) ? lut_rgb[i3] : 255;
      let ice_g = (i3 + 2 < lut_len) ? lut_rgb[i3 + 1] : 255;
      let ice_b = (i3 + 2 < lut_len) ? lut_rgb[i3 + 2] : 255;
      let sea_r = (lut_len >= 3) ? lut_rgb[0] : bg_r;
      let sea_g = (lut_len >= 3) ? lut_rgb[1] : bg_g;
      let sea_b = (lut_len >= 3) ? lut_rgb[2] : bg_b;
      let src_total = src_w*src_h;

      if (!this._cat_rgb_cache || this._cat_rgb_cache.length !== src_total*3)
        this._cat_rgb_cache = new Uint8Array(src_total*3);
      if (!this._cat_alpha_cache || this._cat_alpha_cache.length !== src_total)
        this._cat_alpha_cache = new Float32Array(src_total);

      let src_rgb = this._cat_rgb_cache;
      let src_alpha = this._cat_alpha_cache;

      for (let y = 0; y < src_h; y++) {
        let inv_area = 1.0/cell_area_by_row[y];
        let row_off = y*src_w;
        let y_up = (y > 0) ? (y - 1)*src_w : row_off;
        let y_dn = (y < src_h - 1) ? (y + 1)*src_w : row_off;

        for (let x = 0; x < src_w; x++) {
          let idx = row_off + x;
          let cid = class_grid[idx];
          let o3 = idx*3;

          if (cid === 0) {
            src_alpha[idx] = 0.0;
            src_rgb[o3] = sea_r;
            src_rgb[o3 + 1] = sea_g;
            src_rgb[o3 + 2] = sea_b;
            continue;
          }

          let c3 = cid*3;
          let cr = (c3 + 2 < lut_len) ? lut_rgb[c3] : bg_r;
          let cg = (c3 + 2 < lut_len) ? lut_rgb[c3 + 1] : bg_g;
          let cb = (c3 + 2 < lut_len) ? lut_rgb[c3 + 2] : bg_b;
          let val = land_km2[idx];
          let f = Math.min(1.0, Math.abs(val)*inv_area);

          if (f >= 0.995) {
            src_alpha[idx] = 1.0;
            src_rgb[o3] = cr;
            src_rgb[o3 + 1] = cg;
            src_rgb[o3 + 2] = cb;
          } else {
            let x_l = (x === 0) ? (src_w - 1) : (x - 1);
            let x_r = (x === src_w - 1) ? 0 : (x + 1);
            let n0 = land_km2[row_off + x_l];
            let n1 = land_km2[row_off + x_r];
            let n2 = land_km2[y_up + x];
            let n3 = land_km2[y_dn + x];
            let adj_sea = (n0 === 0.0 || n1 === 0.0 || n2 === 0.0 || n3 === 0.0);
            let adj_ice = (n0 < 0.0 || n1 < 0.0 || n2 < 0.0 || n3 < 0.0);
            let adj_land = (n0 > 0.0 || n1 > 0.0 || n2 > 0.0 || n3 > 0.0);

            if (val > 0.0 && adj_ice && !adj_sea && ice_id > 0) {
              src_alpha[idx] = 1.0;
              src_rgb[o3] = (ice_r + (cr - ice_r)*f + 0.5) | 0;
              src_rgb[o3 + 1] = (ice_g + (cg - ice_g)*f + 0.5) | 0;
              src_rgb[o3 + 2] = (ice_b + (cb - ice_b)*f + 0.5) | 0;
            } else if (val < 0.0 && adj_land && !adj_sea) {
              src_alpha[idx] = 1.0;
              src_rgb[o3] = cr;
              src_rgb[o3 + 1] = cg;
              src_rgb[o3 + 2] = cb;
            } else {
              src_alpha[idx] = f;
              src_rgb[o3] = cr;
              src_rgb[o3 + 1] = cg;
              src_rgb[o3 + 2] = cb;
            }
          }
        }
      }

      for (let i = 0; i < n_pix; i++) {
        let o3 = i*3;
        let a = edge_alpha[i];

        if (a === 0) {
          out_rgb[o3] = bg_r;
          out_rgb[o3 + 1] = bg_g;
          out_rgb[o3 + 2] = bg_b;
          continue;
        }

        let r_f = src_row[i];
        let c_f = src_col[i];
        let r0 = r_f | 0;
        let c0 = c_f | 0;
        let r1 = (r0 + 1 < src_h) ? (r0 + 1) : r0;
        let c1 = (c0 + 1 < src_w) ? (c0 + 1) : 0;
        let dr = r_f - r0;
        let dc = c_f - c0;
        let w00 = (1.0 - dr)*(1.0 - dc);
        let w01 = (1.0 - dr)*dc;
        let w10 = dr*(1.0 - dc);
        let w11 = dr*dc;

        let i00 = r0*src_w + c0;
        let i01 = r0*src_w + c1;
        let i10 = r1*src_w + c0;
        let i11 = r1*src_w + c1;
        let a00 = src_alpha[i00];
        let a01 = src_alpha[i01];
        let a10 = src_alpha[i10];
        let a11 = src_alpha[i11];
        let alpha = w00*a00 + w01*a01 + w10*a10 + w11*a11;

        let pr, pg, pb;
        if (alpha <= 0.001) {
          pr = sea_r;
          pg = sea_g;
          pb = sea_b;
        } else if (alpha >= 0.999) {
          let rn = Math.min(src_h - 1, Math.max(0, Math.round(r_f)));
          let cn = Math.round(c_f)%src_w;
          let n3 = (rn*src_w + cn)*3;
          pr = src_rgb[n3];
          pg = src_rgb[n3 + 1];
          pb = src_rgb[n3 + 2];
        } else {
          let best_i = i00;
          let best_w = w00*a00;
          if (w01*a01 > best_w) { best_w = w01*a01; best_i = i01; }
          if (w10*a10 > best_w) { best_w = w10*a10; best_i = i10; }
          if (w11*a11 > best_w) { best_i = i11; }
          let b3 = best_i*3;
          let cr = src_rgb[b3];
          let cg = src_rgb[b3 + 1];
          let cb = src_rgb[b3 + 2];
          let alpha_vis = Math.sqrt(alpha);
          pr = (sea_r + (cr - sea_r)*alpha_vis + 0.5) | 0;
          pg = (sea_g + (cg - sea_g)*alpha_vis + 0.5) | 0;
          pb = (sea_b + (cb - sea_b)*alpha_vis + 0.5) | 0;
        }

        if (a === 255) {
          out_rgb[o3] = pr;
          out_rgb[o3 + 1] = pg;
          out_rgb[o3 + 2] = pb;
        } else {
          let inv = 255 - a;
          out_rgb[o3] = (pr*a + bg_r*inv + 128) >> 8;
          out_rgb[o3 + 1] = (pg*a + bg_g*inv + 128) >> 8;
          out_rgb[o3 + 2] = (pb*a + bg_b*inv + 128) >> 8;
        }
      }

      //Return statement
      return out_rgb;
    }

    for (let i = 0; i < n_pix; i++) {
      let o3 = i*3;
      let a = edge_alpha[i];

      if (a === 0) {
        out_rgb[o3] = bg_r;
        out_rgb[o3 + 1] = bg_g;
        out_rgb[o3 + 2] = bg_b;
        continue;
      }

      let r = Math.min(src_h - 1, Math.max(0, Math.round(src_row[i])));
      let c = Math.round(src_col[i])%src_w;
      let class_id = class_grid[r*src_w + c];
      let l3 = class_id*3;
      let cr = (l3 + 2 < lut_len) ? lut_rgb[l3] : bg_r;
      let cg = (l3 + 2 < lut_len) ? lut_rgb[l3 + 1] : bg_g;
      let cb = (l3 + 2 < lut_len) ? lut_rgb[l3 + 2] : bg_b;

      if (a === 255) {
        out_rgb[o3] = cr;
        out_rgb[o3 + 1] = cg;
        out_rgb[o3 + 2] = cb;
      } else {
        let inv = 255 - a;
        out_rgb[o3] = (cr*a + bg_r*inv + 128) >> 8;
        out_rgb[o3 + 1] = (cg*a + bg_g*inv + 128) >> 8;
        out_rgb[o3 + 2] = (cb*a + bg_b*inv + 128) >> 8;
      }
    }

    //Return statement
    return out_rgb;
  }
}

module.exports = EqualEarthProjector;
