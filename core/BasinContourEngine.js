"use strict";

let fs = require("node:fs");
let path = require("node:path");

/**
 * 30-arcsecond (~900 m) ETOPO 2022 topographic contour & connected-component basin flood engine.
 * Drives physically accurate shorelines for:
 *  - Closed/sill-controlled basins (Caspian Sea, Black Sea Neoeuxinian, Aral Sea, Lake Chad / Mega-Lake Chad, Lake Urmia)
 *  - Ocean-connected Holocene prograding deltas (Mesopotamia/Persian Gulf, China Bohai/Yellow/Yangtze/PRD, Nile, Po)
 */
class BasinContourEngine {
  /**
   * @param {string} arg0_intermediate_dir
   */
  constructor (arg0_intermediate_dir) {
    //Convert from parameters
    let intermediate_dir = arg0_intermediate_dir;

    //Declare local instance variables
    let meta = JSON.parse(fs.readFileSync(path.join(intermediate_dir, "etopo_meta.json"), "utf8"));
    let buf = fs.readFileSync(path.join(intermediate_dir, "etopo_bed_halfm.i16"));

    this.bed_i16 = new Int16Array(buf.buffer, buf.byteOffset, buf.byteLength/2);
    this.factor_x = meta.factor_x;
    this.factor_y = meta.factor_y;
    this.scale_m = meta.scale_m_per_unit;
    this.src_h = meta.src_h;
    this.src_w = meta.src_w;
    this.sub_per_cell = meta.factor_x*meta.factor_y;
  }

  /**
   * Evaluates a piecewise-linear [year, value] curve at target_year.
   *
   * @param {Array<Array<number>>} arg0_curve
   * @param {number} arg1_year
   *
   * @returns {number}
   */
  static evalCurve (arg0_curve, arg1_year) {
    //Convert from parameters
    let curve = arg0_curve;
    let year = arg1_year;

    //Guard clauses
    if (!curve || curve.length === 0) return 0.0;
    if (year <= curve[0][0]) return curve[0][1];
    if (year >= curve[curve.length - 1][0]) return curve[curve.length - 1][1];

    for (let i = 0; i < curve.length - 1; i++) {
      let y0 = curve[i][0];
      let y1 = curve[i + 1][0];
      if (year >= y0 && year <= y1) {
        if (y1 === y0) return curve[i + 1][1];
        let t = (year - y0)/(y1 - y0);
        return curve[i][1] + t*(curve[i + 1][1] - curve[i][1]);
      }
    }

    //Return statement
    return curve[curve.length - 1][1];
  }

  /**
   * Performs 30" connected-component flood-fill inside bbox from the lowest point near seed
   * for all sub-cells with elevation <= water_level_m, returning sparse 5' [pixel_idx, water_frac] pairs.
   *
   * @param {Object} arg0_bbox - { lat_min, lat_max, lon_min, lon_max }
   * @param {Object} arg1_seed - { lat, lon }
   * @param {number} arg2_water_level_m
   * @param {Object} [arg3_options]
   *  @param {boolean} [arg3_options.seed_edge_sea=false] - Seed all sub-cells <= 0 on bbox sea edge / below sea
   *  @param {boolean} [arg3_options.clean_stipples=false] - Morphological closing of tiny interior speckles
   *  @param {boolean} [arg3_options.smooth_substep=false] - Interpolate spatially between 0.5m DEM contours via dual-BFS shoreline distance
   *
   * @returns {Array<Array<number>>} [[pixel_idx_5m, water_frac], ...]
   */
  floodSeededContour (arg0_bbox, arg1_seed, arg2_water_level_m, arg3_options) {
    //Convert from parameters
    let bbox = arg0_bbox;
    let seed = arg1_seed;
    let water_level_m = arg2_water_level_m;
    let options = arg3_options ? arg3_options : {};

    //Declare local instance variables
    let bed = this.bed_i16;
    let fx = this.factor_x;
    let fy = this.factor_y;
    let r0 = Math.max(0, Math.floor(((90.0 - bbox.lat_max)/180.0)*2160));
    let r1 = Math.min(2160, Math.ceil(((90.0 - bbox.lat_min)/180.0)*2160));
    let c0 = Math.max(0, Math.floor(((bbox.lon_min + 180.0)/360.0)*4320));
    let c1 = Math.min(4320, Math.ceil(((bbox.lon_max + 180.0)/360.0)*4320));
    let sy0 = r0*fy;
    let sy1 = r1*fy;
    let sx0 = c0*fx;
    let sx1 = c1*fx;
    let h = sy1 - sy0;
    let w = sx1 - sx0;
    let lv_exact = water_level_m/this.scale_m;
    let lv_lo = Math.floor(lv_exact);
    let frac_sub = lv_exact - lv_lo;
    let do_substep = Boolean(options.smooth_substep) && frac_sub > 0.02 && frac_sub < 0.98;
    let lv_hm = do_substep ? (lv_lo + 1) : Math.round(lv_exact);
    let src_w = this.src_w;

    let vis = new Uint8Array(w*h);
    let q = new Int32Array(w*h);
    let qh = 0;
    let qt = 0;

    if (options.seed_edge_sea) {
      //Seed every modern marine sub-cell (z <= 0) inside bbox that is reachable from the sea seed
      let sr = Math.min(h - 1, Math.max(0, Math.floor(((90.0 - seed.lat)/180.0)*this.src_h) - sy0));
      let sc = Math.min(w - 1, Math.max(0, Math.floor(((seed.lon + 180.0)/360.0)*this.src_w) - sx0));

      //Locate lowest sub-cell within 25 sub-pixels (~22 km) of sea_seed
      let best_z = 32767;
      let best_r = sr;
      let best_c = sc;
      for (let dy = -25; dy <= 25; dy++) {
        for (let dx = -25; dx <= 25; dx++) {
          let rr = sr + dy;
          let cc = sc + dx;
          if (rr < 0 || rr >= h || cc < 0 || cc >= w) continue;
          let z = bed[(sy0 + rr)*src_w + (sx0 + cc)];
          if (z < best_z) {
            best_z = z;
            best_r = rr;
            best_c = cc;
          }
        }
      }

      if (best_z <= lv_hm) {
        let s_idx = best_r*w + best_c;
        vis[s_idx] = 1;
        q[qt++] = s_idx;
      }
    } else {
      let sr = Math.min(h - 1, Math.max(0, Math.floor(((90.0 - seed.lat)/180.0)*this.src_h) - sy0));
      let sc = Math.min(w - 1, Math.max(0, Math.floor(((seed.lon + 180.0)/360.0)*this.src_w) - sx0));
      let best_z = 32767;
      let best_r = sr;
      let best_c = sc;

      for (let dy = -20; dy <= 20; dy++) {
        for (let dx = -20; dx <= 20; dx++) {
          let rr = sr + dy;
          let cc = sc + dx;
          if (rr < 0 || rr >= h || cc < 0 || cc >= w) continue;
          let z = bed[(sy0 + rr)*src_w + (sx0 + cc)];
          if (z < best_z) {
            best_z = z;
            best_r = rr;
            best_c = cc;
          }
        }
      }

      if (best_z <= lv_hm) {
        let s_idx = best_r*w + best_c;
        vis[s_idx] = 1;
        q[qt++] = s_idx;
      }
    }

    //Iterative 4-neighbour BFS at 30" (~900 m) resolution
    for (; qh < qt; ) {
      let cur = q[qh++];
      let cy = (cur/w) | 0;
      let cx = cur - cy*w;
      let g_row = (sy0 + cy)*src_w + sx0;

      if (cx > 0) {
        let ni = cur - 1;
        if (vis[ni] === 0 && bed[g_row + cx - 1] <= lv_hm) {
          vis[ni] = 1;
          q[qt++] = ni;
        }
      }
      if (cx < w - 1) {
        let ni = cur + 1;
        if (vis[ni] === 0 && bed[g_row + cx + 1] <= lv_hm) {
          vis[ni] = 1;
          q[qt++] = ni;
        }
      }
      if (cy > 0) {
        let ni = cur - w;
        if (vis[ni] === 0 && bed[g_row - src_w + cx] <= lv_hm) {
          vis[ni] = 1;
          q[qt++] = ni;
        }
      }
      if (cy < h - 1) {
        let ni = cur + w;
        if (vis[ni] === 0 && bed[g_row + src_w + cx] <= lv_hm) {
          vis[ni] = 1;
          q[qt++] = ni;
        }
      }
    }

    //Sub-0.5m spatial shoreline interpolation between lv_lo and lv_lo + 1 via dual-BFS sea/rim distance
    if (do_substep && qt > 0) {
      let core = new Uint8Array(w*h);
      let cq = new Int32Array(w*h);
      let cqh = 0;
      let cqt = 0;
      let s0 = q[0];
      let sr0 = (s0/w) | 0;
      let sc0 = s0 - sr0*w;
      if (bed[(sy0 + sr0)*src_w + (sx0 + sc0)] <= lv_lo) {
        core[s0] = 1;
        cq[cqt++] = s0;
      }

      for (; cqh < cqt; ) {
        let cur = cq[cqh++];
        let cy = (cur/w) | 0;
        let cx = cur - cy*w;
        let g_row = (sy0 + cy)*src_w + sx0;

        if (cx > 0) {
          let ni = cur - 1;
          if (core[ni] === 0 && vis[ni] === 1 && bed[g_row + cx - 1] <= lv_lo) {
            core[ni] = 1;
            cq[cqt++] = ni;
          }
        }
        if (cx < w - 1) {
          let ni = cur + 1;
          if (core[ni] === 0 && vis[ni] === 1 && bed[g_row + cx + 1] <= lv_lo) {
            core[ni] = 1;
            cq[cqt++] = ni;
          }
        }
        if (cy > 0) {
          let ni = cur - w;
          if (core[ni] === 0 && vis[ni] === 1 && bed[g_row - src_w + cx] <= lv_lo) {
            core[ni] = 1;
            cq[cqt++] = ni;
          }
        }
        if (cy < h - 1) {
          let ni = cur + w;
          if (core[ni] === 0 && vis[ni] === 1 && bed[g_row + src_w + cx] <= lv_lo) {
            core[ni] = 1;
            cq[cqt++] = ni;
          }
        }
      }

      let d_sea = new Int32Array(w*h);
      let d_rim = new Int32Array(w*h);
      let sq = new Int32Array(w*h);
      let rq = new Int32Array(w*h);
      let sqh = 0;
      let sqt = 0;
      let rqh = 0;
      let rqt = 0;

      for (let idx = 0; idx < w*h; idx++) {
        if (vis[idx] === 0 || core[idx] === 1) continue;
        let cy = (idx/w) | 0;
        let cx = idx - cy*w;
        let near_core = ((cx > 0 && core[idx - 1] === 1) ||
          (cx < w - 1 && core[idx + 1] === 1) ||
          (cy > 0 && core[idx - w] === 1) ||
          (cy < h - 1 && core[idx + w] === 1));
        if (near_core) {
          d_sea[idx] = 1;
          sq[sqt++] = idx;
        }
        let near_dry = ((cx > 0 && vis[idx - 1] === 0) ||
          (cx < w - 1 && vis[idx + 1] === 0) ||
          (cy > 0 && vis[idx - w] === 0) ||
          (cy < h - 1 && vis[idx + w] === 0));
        if (near_dry) {
          d_rim[idx] = 1;
          rq[rqt++] = idx;
        }
      }

      for (; sqh < sqt; ) {
        let cur = sq[sqh++];
        let cy = (cur/w) | 0;
        let cx = cur - cy*w;
        let nd = d_sea[cur] + 1;
        if (cx > 0 && vis[cur - 1] === 1 && core[cur - 1] === 0 && d_sea[cur - 1] === 0) {
          d_sea[cur - 1] = nd;
          sq[sqt++] = cur - 1;
        }
        if (cx < w - 1 && vis[cur + 1] === 1 && core[cur + 1] === 0 && d_sea[cur + 1] === 0) {
          d_sea[cur + 1] = nd;
          sq[sqt++] = cur + 1;
        }
        if (cy > 0 && vis[cur - w] === 1 && core[cur - w] === 0 && d_sea[cur - w] === 0) {
          d_sea[cur - w] = nd;
          sq[sqt++] = cur - w;
        }
        if (cy < h - 1 && vis[cur + w] === 1 && core[cur + w] === 0 && d_sea[cur + w] === 0) {
          d_sea[cur + w] = nd;
          sq[sqt++] = cur + w;
        }
      }

      for (; rqh < rqt; ) {
        let cur = rq[rqh++];
        let cy = (cur/w) | 0;
        let cx = cur - cy*w;
        let nd = d_rim[cur] + 1;
        if (cx > 0 && vis[cur - 1] === 1 && core[cur - 1] === 0 && d_rim[cur - 1] === 0) {
          d_rim[cur - 1] = nd;
          rq[rqt++] = cur - 1;
        }
        if (cx < w - 1 && vis[cur + 1] === 1 && core[cur + 1] === 0 && d_rim[cur + 1] === 0) {
          d_rim[cur + 1] = nd;
          rq[rqt++] = cur + 1;
        }
        if (cy > 0 && vis[cur - w] === 1 && core[cur - w] === 0 && d_rim[cur - w] === 0) {
          d_rim[cur - w] = nd;
          rq[rqt++] = cur - w;
        }
        if (cy < h - 1 && vis[cur + w] === 1 && core[cur + w] === 0 && d_rim[cur + w] === 0) {
          d_rim[cur + w] = nd;
          rq[rqt++] = cur + w;
        }
      }

      for (let idx = 0; idx < w*h; idx++) {
        if (vis[idx] === 0 || core[idx] === 1) continue;
        let ds = d_sea[idx];
        let dr = d_rim[idx];
        if (ds === 0) {
          vis[idx] = 0;
        } else {
          let eff_dr = (dr > 0) ? dr : 9999;
          let u = (ds - 0.5)/(ds + eff_dr - 1.0);
          if (u > frac_sub) vis[idx] = 0;
        }
      }
    }

    //Optional morphological closing (1-pass fill of tiny interior DEM speckles) for smooth alluvial delta plains
    if (options.clean_stipples) {
      let tmp = new Uint8Array(vis);
      for (let cy = 1; cy < h - 1; cy++) {
        let row = cy*w;
        for (let cx = 1; cx < w - 1; cx++) {
          let i = row + cx;
          if (vis[i] === 1) continue;
          let nbr = vis[i - 1] + vis[i + 1] + vis[i - w] + vis[i + w] +
            vis[i - w - 1] + vis[i - w + 1] + vis[i + w - 1] + vis[i + w + 1];
          if (nbr >= 5) tmp[i] = 1;
        }
      }
      vis = tmp;
    }

    //Aggregate 10x10 30" sub-cells into 5' sparse [pixel_idx, water_frac] pairs
    let pairs = [];
    let out_h = r1 - r0;
    let out_w = c1 - c0;
    let inv_sub = 1.0/this.sub_per_cell;

    for (let ly = 0; ly < out_h; ly++) {
      let gy = r0 + ly;
      for (let lx = 0; lx < out_w; lx++) {
        let hits = 0;
        let sub_r0 = ly*fy;
        let sub_c0 = lx*fx;

        for (let dy = 0; dy < fy; dy++) {
          let row_off = (sub_r0 + dy)*w + sub_c0;
          for (let dx = 0; dx < fx; dx++)
            if (vis[row_off + dx] === 1) hits++;
        }

        if (hits > 0) {
          let gx = c0 + lx;
          pairs.push([gy*4320 + gx, Math.round(hits*inv_sub*10000)/10000]);
        }
      }
    }

    //Return statement
    return pairs;
  }

  /**
   * Floods an ice-dammed proglacial lake on the 30" DEM against an isostatically tilted water plane.
   * A sub-cell is water when it is 4-connected to a seed, its bedrock lies above sea level and below the tilted plane, and its 5' parent cell is not glaciated.
   *
   * @param {Object} arg0_spec - { bbox, seeds: [{ lat, lon }], tilt: { ref_lat, ref_lon, azimuth_deg } }
   * @param {number} arg1_level_ref_m - Water plane elevation at the tilt reference point (m).
   * @param {number} arg2_gradient_m_per_km - Water plane rise per km along the tilt azimuth.
   * @param {Float32Array} arg3_ice_5m - 4320x2160 glaciated fraction for the target year.
   *
   * @returns {Array<Array<number>>} [[pixel_idx_5m, water_frac], ...]
   */
  floodTiltedProglacial (arg0_spec, arg1_level_ref_m, arg2_gradient_m_per_km, arg3_ice_5m) {
    //Convert from parameters
    let spec = arg0_spec;
    let level_ref_m = arg1_level_ref_m;
    let gradient = arg2_gradient_m_per_km;
    let ice_5m = arg3_ice_5m;

    //Declare local instance variables
    let az = (spec.tilt.azimuth_deg*Math.PI)/180.0;
    let bbox = spec.bbox;
    let bed = this.bed_i16;
    let fx = this.factor_x;
    let fy = this.factor_y;
    let r0 = Math.max(0, Math.floor(((90.0 - bbox.lat_max)/180.0)*2160));
    let r1 = Math.min(2160, Math.ceil(((90.0 - bbox.lat_min)/180.0)*2160));
    let c0 = Math.max(0, Math.floor(((bbox.lon_min + 180.0)/360.0)*4320));
    let c1 = Math.min(4320, Math.ceil(((bbox.lon_max + 180.0)/360.0)*4320));
    let sy0 = r0*fy;
    let sx0 = c0*fx;
    let h = (r1 - r0)*fy;
    let w = (c1 - c0)*fx;
    let src_w = this.src_w;
    let scale_m = this.scale_m;
    let sub_deg = 180.0/this.src_h;
    let km_per_deg = 111.2;

    //Precompute plane offsets: level(row, col) = row_level[row] + col_slope[row]*col
    let row_level = new Float64Array(h);
    let col_slope = new Float64Array(h);
    for (let y = 0; y < h; y++) {
      let lat = 90.0 - (sy0 + y + 0.5)*sub_deg;
      let lon_left = -180.0 + (sx0 + 0.5)*sub_deg;
      let kx = km_per_deg*Math.cos((lat*Math.PI)/180.0)*Math.sin(az);
      let d_km = (lat - spec.tilt.ref_lat)*km_per_deg*Math.cos(az) + (lon_left - spec.tilt.ref_lon)*kx;
      row_level[y] = level_ref_m + gradient*d_km;
      col_slope[y] = gradient*kx*sub_deg;
    }

    let passable = (arg0_y, arg1_x) => {
      let z_m = bed[(sy0 + arg0_y)*src_w + sx0 + arg1_x]*scale_m;
      if (z_m <= 0.0 || z_m > row_level[arg0_y] + col_slope[arg0_y]*arg1_x) return false;
      let p_idx = (r0 + ((arg0_y/fy) | 0))*4320 + c0 + ((arg1_x/fx) | 0);

      //Return statement
      return ice_5m[p_idx] < 0.5;
    };

    let vis = new Uint8Array(w*h);
    let q = new Int32Array(w*h);
    let qh = 0;
    let qt = 0;

    //Seed from the lowest passable sub-cell near each seed point
    for (let i = 0; i < spec.seeds.length; i++) {
      let sr = Math.floor(((90.0 - spec.seeds[i].lat)/180.0)*this.src_h) - sy0;
      let sc = Math.floor(((spec.seeds[i].lon + 180.0)/360.0)*this.src_w) - sx0;
      let best = -1;
      let best_z = 32767;

      for (let dy = -20; dy <= 20; dy++)
        for (let dx = -20; dx <= 20; dx++) {
          let rr = sr + dy;
          let cc = sc + dx;
          if (rr < 0 || rr >= h || cc < 0 || cc >= w) continue;
          let z = bed[(sy0 + rr)*src_w + sx0 + cc];
          if (z < best_z && passable(rr, cc)) {
            best_z = z;
            best = rr*w + cc;
          }
        }

      if (best >= 0 && vis[best] === 0) {
        vis[best] = 1;
        q[qt++] = best;
      }
    }

    //Iterative 4-neighbour BFS at 30" resolution
    for (; qh < qt; ) {
      let cur = q[qh++];
      let cy = (cur/w) | 0;
      let cx = cur - cy*w;
      let nbrs = [[cy, cx - 1], [cy, cx + 1], [cy - 1, cx], [cy + 1, cx]];

      for (let k = 0; k < 4; k++) {
        let ny = nbrs[k][0];
        let nx = nbrs[k][1];
        if (ny < 0 || ny >= h || nx < 0 || nx >= w) continue;
        let ni = ny*w + nx;
        if (vis[ni] === 0 && passable(ny, nx)) {
          vis[ni] = 1;
          q[qt++] = ni;
        }
      }
    }

    //Two-pass majority closing to remove single sub-cell DEM speckle along the shoreline
    for (let pass = 0; pass < 2; pass++) {
      let tmp = new Uint8Array(vis);
      for (let cy = 1; cy < h - 1; cy++) {
        let row = cy*w;
        for (let cx = 1; cx < w - 1; cx++) {
          let i = row + cx;
          let nbr = vis[i - 1] + vis[i + 1] + vis[i - w] + vis[i + w] +
            vis[i - w - 1] + vis[i - w + 1] + vis[i + w - 1] + vis[i + w + 1];
          if (vis[i] === 0 && nbr >= 5) tmp[i] = 1;
          if (vis[i] === 1 && nbr <= 2) tmp[i] = 0;
        }
      }
      vis = tmp;
    }

    //Aggregate 30" sub-cells into 5' sparse [pixel_idx, water_frac] pairs, never over glaciated parents
    let pairs = [];
    let inv_sub = 1.0/this.sub_per_cell;
    for (let ly = 0; ly < r1 - r0; ly++) {
      for (let lx = 0; lx < c1 - c0; lx++) {
        let p_idx = (r0 + ly)*4320 + c0 + lx;
        if (ice_5m[p_idx] >= 0.5) continue;
        let hits = 0;
        for (let dy = 0; dy < fy; dy++) {
          let row_off = (ly*fy + dy)*w + lx*fx;
          for (let dx = 0; dx < fx; dx++)
            if (vis[row_off + dx] === 1) hits++;
        }
        if (hits > 0)
          pairs.push([p_idx, Math.round(hits*inv_sub*10000)/10000]);
      }
    }

    //Return statement
    return pairs;
  }


  /**
   * Computes a fast 2-pass Chamfer Euclidean distance transform from all 1-pixels in a binary grid.
   *
   * @param {Uint8Array} arg0_mask
   * @param {number} arg1_w
   * @param {number} arg2_h
   *
   * @returns {Float32Array}
   */
  static computeChamferDistance (arg0_mask, arg1_w, arg2_h) {
    //Convert from parameters
    let mask = arg0_mask;
    let w = arg1_w;
    let h = arg2_h;

    //Declare local instance variables
    let d_diag = 1.41421356;
    let d_orth = 1.0;
    let dist = new Float32Array(w*h);
    let inf = 1e6;

    for (let i = 0; i < dist.length; i++)
      dist[i] = (mask[i] === 1) ? 0.0 : inf;

    //Forward pass (top-left -> bottom-right)
    for (let y = 0; y < h; y++) {
      let row = y*w;
      for (let x = 0; x < w; x++) {
        let i = row + x;
        let v = dist[i];
        if (v === 0.0) continue;
        if (x > 0 && dist[i - 1] + d_orth < v) v = dist[i - 1] + d_orth;
        if (y > 0) {
          if (dist[i - w] + d_orth < v) v = dist[i - w] + d_orth;
          if (x > 0 && dist[i - w - 1] + d_diag < v) v = dist[i - w - 1] + d_diag;
          if (x < w - 1 && dist[i - w + 1] + d_diag < v) v = dist[i - w + 1] + d_diag;
        }
        dist[i] = v;
      }
    }

    //Backward pass (bottom-right -> top-left)
    for (let y = h - 1; y >= 0; y--) {
      let row = y*w;
      for (let x = w - 1; x >= 0; x--) {
        let i = row + x;
        let v = dist[i];
        if (v === 0.0) continue;
        if (x < w - 1 && dist[i + 1] + d_orth < v) v = dist[i + 1] + d_orth;
        if (y < h - 1) {
          if (dist[i + w] + d_orth < v) v = dist[i + w] + d_orth;
          if (x > 0 && dist[i + w - 1] + d_diag < v) v = dist[i + w - 1] + d_diag;
          if (x < w - 1 && dist[i + w + 1] + d_diag < v) v = dist[i + w + 1] + d_diag;
        }
        dist[i] = v;
      }
    }

    //Return statement
    return dist;
  }

  /**
   * Computes a 2D Signed Distance Field (+ inside binary mask, - outside binary mask).
   *
   * @param {Uint8Array} arg0_mask
   * @param {number} arg1_w
   * @param {number} arg2_h
   *
   * @returns {Float32Array}
   */
  static computeSignedDistanceField (arg0_mask, arg1_w, arg2_h) {
    //Convert from parameters
    let mask = arg0_mask;
    let w = arg1_w;
    let h = arg2_h;

    //Declare local instance variables
    let inv = new Uint8Array(w*h);
    for (let i = 0; i < mask.length; i++)
      inv[i] = (mask[i] === 1) ? 0 : 1;

    let d_out = this.computeChamferDistance(mask, w, h);
    let d_in = this.computeChamferDistance(inv, w, h);
    let sdf = new Float32Array(w*h);

    for (let i = 0; i < sdf.length; i++)
      sdf[i] = (mask[i] === 1) ? d_in[i] : -d_out[i];

    //Return statement
    return sdf;
  }

  /**
   * Rasterises a GeoJSON Polygon/MultiPolygon onto a local 30" sub-grid inside bbox.
   *
   * @param {Object} arg0_geometry
   * @param {number} arg1_sy0
   * @param {number} arg2_sx0
   * @param {number} arg3_w
   * @param {number} arg4_h
   *
   * @returns {Uint8Array}
   */
  rasteriseGeometryOn30s (arg0_geometry, arg1_sy0, arg2_sx0, arg3_w, arg4_h) {
    //Convert from parameters
    let geometry = arg0_geometry;
    let sy0 = arg1_sy0;
    let sx0 = arg2_sx0;
    let w = arg3_w;
    let h = arg4_h;

    //Declare local instance variables
    let mask = new Uint8Array(w*h);

    //Guard clause
    if (!geometry || !geometry.coordinates) return mask;

    let polys = (geometry.type === "Polygon") ? [geometry.coordinates] : (geometry.type === "MultiPolygon" ? geometry.coordinates : []);
    let scale_x = this.src_w/360.0;
    let scale_y = this.src_h/180.0;
    let xs_list = [];

    for (let p = 0; p < polys.length; p++) {
      let rings = polys[p];
      let edge_x0 = [], edge_y0 = [], edge_x1 = [], edge_y1 = [];
      let p_min_y = h;
      let p_max_y = -1;

      for (let r = 0; r < rings.length; r++) {
        let ring = rings[r];
        if (!ring || ring.length < 3) continue;
        for (let k = 0; k < ring.length - 1; k++) {
          let y0 = (90.0 - ring[k][1])*scale_y - sy0;
          let y1 = (90.0 - ring[k + 1][1])*scale_y - sy0;
          if (y0 === y1) continue;
          let x0 = (ring[k][0] + 180.0)*scale_x - sx0;
          let x1 = (ring[k + 1][0] + 180.0)*scale_x - sx0;
          if (y0 < p_min_y) p_min_y = y0;
          if (y1 < p_min_y) p_min_y = y1;
          if (y0 > p_max_y) p_max_y = y0;
          if (y1 > p_max_y) p_max_y = y1;
          edge_x0.push(x0);
          edge_y0.push(y0);
          edge_x1.push(x1);
          edge_y1.push(y1);
        }
      }

      let n_edges = edge_x0.length;
      if (n_edges === 0) continue;

      let y_start = Math.max(0, Math.floor(p_min_y));
      let y_end = Math.min(h - 1, Math.ceil(p_max_y));

      for (let y = y_start; y <= y_end; y++) {
        let scan_y = y + 0.5;
        xs_list.length = 0;
        for (let i = 0; i < n_edges; i++) {
          let ey0 = edge_y0[i];
          let ey1 = edge_y1[i];
          if ((ey0 <= scan_y && ey1 > scan_y) || (ey1 <= scan_y && ey0 > scan_y)) {
            let t = (scan_y - ey0)/(ey1 - ey0);
            xs_list.push(edge_x0[i] + t*(edge_x1[i] - edge_x0[i]));
          }
        }
        if (xs_list.length < 2) continue;
        xs_list.sort((a, b) => a - b);
        let row_off = y*w;
        for (let i = 0; i < xs_list.length - 1; i += 2) {
          let xs = Math.max(0, Math.ceil(xs_list[i] - 0.5));
          let xe = Math.min(w - 1, Math.floor(xs_list[i + 1] - 0.5));
          for (let x = xs; x <= xe; x++)
            mask[row_off + x] = 1;
        }
      }
    }

    //Return statement
    return mask;
  }

  /**
   * Builds crisp, zero-ghost, zero-stipple migrating paleoshoreline water-addition masks across all
   * target years using 30" Signed Distance Field (SDF) interpolation between historically calibrated
   * keyframe paleocoastline polygons while preserving bedrock hills above max_alluvial_elev_m.
   *
   * @param {Object} arg0_spec - { bbox, max_alluvial_elev_m, min_year, max_year, preserve_hills, keyframes: [{ year, geometry }] }
   * @param {Array<number>} arg1_years
   *
   * @returns {Object} { bbox, level_table, year_to_level_id }
   */
  buildSdfDeltaSeries (arg0_spec, arg1_years) {
    //Convert from parameters
    let spec = arg0_spec;
    let years = arg1_years;

    //Declare local instance variables
    let bbox = spec.bbox;
    let bed = this.bed_i16;
    let fx = this.factor_x;
    let fy = this.factor_y;
    let inv_sub = 1.0/this.sub_per_cell;
    let max_elev_hm = Math.round(((spec.max_alluvial_elev_m !== undefined) ? spec.max_alluvial_elev_m : 14.0)/this.scale_m);
    let preserve_hills = (spec.preserve_hills !== false);
    let r0 = Math.max(0, Math.floor(((90.0 - bbox.lat_max)/180.0)*2160));
    let r1 = Math.min(2160, Math.ceil(((90.0 - bbox.lat_min)/180.0)*2160));
    let c0 = Math.max(0, Math.floor(((bbox.lon_min + 180.0)/360.0)*4320));
    let c1 = Math.min(4320, Math.ceil(((bbox.lon_max + 180.0)/360.0)*4320));
    let out_h = r1 - r0;
    let out_w = c1 - c0;
    let sy0 = r0*fy;
    let sx0 = c0*fx;
    let h = out_h*fy;
    let w = out_w*fx;
    let src_w = this.src_w;

    //Build 30" ocean-connected modern sea mask & alluvial allow mask (filtering out genuine rocky hills > max_alluvial_elev_m)
    let allow = new Uint8Array(w*h);
    let modern_sea = new Uint8Array(w*h);
    let sq = new Int32Array(w*h);
    let sqh = 0;
    let sqt = 0;

    for (let cy = 0; cy < h; cy++) {
      let g_row = (sy0 + cy)*src_w + sx0;
      let l_row = cy*w;
      for (let cx = 0; cx < w; cx++) {
        let z = bed[g_row + cx];
        allow[l_row + cx] = (!preserve_hills || z <= max_elev_hm) ? 1 : 0;
        //Seed ocean-connected modern sea from deep water (z <= -6 -> -3m) or bbox boundary sea cells (z <= 0)
        let on_edge = (cy === 0 || cy === h - 1 || cx === 0 || cx === w - 1);
        if ((on_edge && z <= 0) || z <= -10) {
          let li = l_row + cx;
          modern_sea[li] = 1;
          sq[sqt++] = li;
        }
      }
    }

    for (; sqh < sqt; ) {
      let cur = sq[sqh++];
      let cy = (cur/w) | 0;
      let cx = cur - cy*w;
      let g_row = (sy0 + cy)*src_w + sx0;

      if (cx > 0 && modern_sea[cur - 1] === 0 && bed[g_row + cx - 1] <= 0) {
        modern_sea[cur - 1] = 1; sq[sqt++] = cur - 1;
      }
      if (cx < w - 1 && modern_sea[cur + 1] === 0 && bed[g_row + cx + 1] <= 0) {
        modern_sea[cur + 1] = 1; sq[sqt++] = cur + 1;
      }
      if (cy > 0 && modern_sea[cur - w] === 0 && bed[g_row - src_w + cx] <= 0) {
        modern_sea[cur - w] = 1; sq[sqt++] = cur - w;
      }
      if (cy < h - 1 && modern_sea[cur + w] === 0 && bed[g_row + src_w + cx] <= 0) {
        modern_sea[cur + w] = 1; sq[sqt++] = cur + w;
      }
    }

    //Clean single-pixel hill noise inside allow mask so bedrock islands are solid and alluvial plains have zero stipples
    if (preserve_hills) {
      let clean_allow = new Uint8Array(allow);
      for (let cy = 1; cy < h - 1; cy++) {
        let row = cy*w;
        for (let cx = 1; cx < w - 1; cx++) {
          let i = row + cx;
          let s = allow[i - 1] + allow[i + 1] + allow[i - w] + allow[i + w];
          if (allow[i] === 0 && s >= 3) clean_allow[i] = 1;
          else if (allow[i] === 1 && s <= 1) clean_allow[i] = 0;
        }
      }
      allow = clean_allow;
    }

    //Precompute 30" SDF for each keyframe polygon (and append exact modern_sea SDF at end_year if modern anchor requested)
    let kf_sorted = spec.keyframes.slice().sort((a, b) => a.year - b.year);
    let sdf_frames = [];
    for (let i = 0; i < kf_sorted.length; i++) {
      let kf = kf_sorted[i];
      let mask = kf.use_modern_sea
        ? modern_sea
        : this.rasteriseGeometryOn30s(kf.geometry, sy0, sx0, w, h);
      let sdf = BasinContourEngine.computeSignedDistanceField(mask, w, h);
      sdf_frames.push({ sdf: sdf, year: kf.year });
    }

    let level_table = [];
    let year_to_level_id = {};
    let min_year = (spec.min_year !== undefined) ? spec.min_year : -5500;
    let max_year = (spec.max_year !== undefined) ? spec.max_year : 2024;
    let cache = new Map();

    for (let yi = 0; yi < years.length; yi++) {
      let yr = years[yi];
      if (yr < min_year || yr > max_year) continue;

      let f0 = sdf_frames[0];
      let f1 = sdf_frames[0];
      if (yr <= sdf_frames[0].year) {
        f0 = sdf_frames[0]; f1 = sdf_frames[0];
      } else if (yr >= sdf_frames[sdf_frames.length - 1].year) {
        f0 = sdf_frames[sdf_frames.length - 1]; f1 = f0;
      } else {
        for (let k = 0; k < sdf_frames.length - 1; k++) {
          if (yr >= sdf_frames[k].year && yr <= sdf_frames[k + 1].year) {
            f0 = sdf_frames[k];
            f1 = sdf_frames[k + 1];
            break;
          }
        }
      }

      let w1_raw = (f0.year === f1.year) ? 0.0 : (yr - f0.year)/(f1.year - f0.year);
      //Quantise blend weight to 2% steps for deduplication without visual stepping
      let w1 = Math.round(w1_raw*50.0)/50.0;
      let w0 = 1.0 - w1;
      let key = `${f0.year}_${f1.year}_${w1.toFixed(2)}`;

      if (!cache.has(key)) {
        let s0 = f0.sdf;
        let s1 = f1.sdf;
        let f_new_grid = new Float32Array(out_h*out_w);
        let f_tot_grid = new Float32Array(out_h*out_w);
        let pairs = [];

        for (let ly = 0; ly < out_h; ly++) {
          let sub_r0 = ly*fy;
          for (let lx = 0; lx < out_w; lx++) {
            let sub_c0 = lx*fx;
            let hits_new = 0;
            let hits_tot = 0;

            for (let dy = 0; dy < fy; dy++) {
              let row_off = (sub_r0 + dy)*w + sub_c0;
              for (let dx = 0; dx < fx; dx++) {
                let si = row_off + dx;
                if (allow[si] === 0) continue;
                let d_val = w0*s0[si] + w1*s1[si];
                if (d_val > 0.0) {
                  hits_tot++;
                  if (modern_sea[si] === 0) hits_new++;
                }
              }
            }

            let li = ly*out_w + lx;
            f_new_grid[li] = hits_new*inv_sub;
            f_tot_grid[li] = hits_tot*inv_sub;
          }
        }

        for (let ly = 0; ly < out_h; ly++) {
          let gy = r0 + ly;
          let y_a = Math.max(0, ly - 3);
          let y_b = Math.min(out_h - 1, ly + 3);
          for (let lx = 0; lx < out_w; lx++) {
            let li = ly*out_w + lx;
            let f_new = f_new_grid[li];
            let f_tot = f_tot_grid[li];
            let max_nb_new = 0.0;

            if (f_tot >= 0.85) {
              let x_a = Math.max(0, lx - 3);
              let x_b = Math.min(out_w - 1, lx + 3);
              for (let ny = y_a; ny <= y_b; ny++) {
                let n_off = ny*out_w;
                for (let nx = x_a; nx <= x_b; nx++) {
                  let v = f_new_grid[n_off + nx];
                  if (v > max_nb_new) max_nb_new = v;
                }
              }
            }

            let loss_zero = (f_tot >= 0.85 && max_nb_new > 0.01) ? Math.min(1.0, max_nb_new/0.25) : 0.0;
            if (f_new > 0.01 || loss_zero > 0.01) {
              let gx = c0 + lx;
              pairs.push([gy*4320 + gx, Math.round(f_new*10000)/10000, Math.round(loss_zero*10000)/10000]);
            }
          }
        }

        cache.set(key, level_table.length);
        level_table.push({ key: key, pairs: pairs });
      }

      year_to_level_id[yr] = cache.get(key);
    }

    //Return statement
    return {
      bbox: bbox,
      level_table: level_table,
      mode: "proportional",
      year_to_level_id: year_to_level_id,
    };
  }
}

module.exports = BasinContourEngine;
