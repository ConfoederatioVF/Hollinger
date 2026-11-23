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
            let isolated_atoll = (n0 === 0.0 && n1 === 0.0 && n2 === 0.0 && n3 === 0.0);

            if (adj_ice && !adj_sea) {
              //Inland land-glacier margin: blend between glacier greyish-white and land lime (never blue sea!)
              src_rgb[o3] = (gl_r + (ln_r - gl_r)*f + 0.5) | 0;
              src_rgb[o3 + 1] = (gl_g + (ln_g - gl_g)*f + 0.5) | 0;
              src_rgb[o3 + 2] = (gl_b + (ln_b - gl_b)*f + 0.5) | 0;
            } else {
              let f_vis = f;
              if (!isolated_atoll) {
                if (f <= 0.42) {
                  f_vis = 0.0;
                } else if (f >= 0.62) {
                  f_vis = 1.0;
                } else {
                  let t = (f - 0.42)/0.20;
                  f_vis = t*t*(3.0 - 2.0*t);
                }
              }
              src_rgb[o3] = (sea_r + (ln_r - sea_r)*f_vis + 0.5) | 0;
              src_rgb[o3 + 1] = (sea_g + (ln_g - sea_g)*f_vis + 0.5) | 0;
              src_rgb[o3 + 2] = (sea_b + (ln_b - sea_b)*f_vis + 0.5) | 0;
            }
          }
        } else if (val < 0.0) {
          let f = gl_blend ? Math.min(1.0, (-val)*inv_area) : 1.0;
          src_rgb[o3] = (sea_r + (gl_r - sea_r)*f + 0.5) | 0;
          src_rgb[o3 + 1] = (sea_g + (gl_g - sea_g)*f + 0.5) | 0;
          src_rgb[o3 + 2] = (sea_b + (gl_b - sea_b)*f + 0.5) | 0;
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
}

module.exports = EqualEarthProjector;
