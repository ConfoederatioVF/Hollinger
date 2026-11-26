"use strict";

/**
 * Holocene estuarine and deltaic sediment back-stripping model for built-up Arctic river mouths.
 *
 * Modern ETOPO bedrock inside infilled estuaries (Ob, Taz, Yenisei, Khatanga) and Holocene deltas (Lena, Yana, Indigirka, Kolyma) records the post-transgression sediment surface, not the pre-Holocene incised valley floor. This model reconstructs the palaeo-surface as `bed - thickness*(1 - completion(t))`, so estuaries flood when sea level reaches the incised valley floor and deltas prograde as infill completes.
 */
class HoloceneInfillModel {
  /**
   * @param {Object} arg0_spec - holocene_infill section of hydro_geomorph_basins.json5.
   * @param {Int16Array} arg1_bed_i16 - 30" bedrock DEM in half-metre units.
   * @param {Object} arg2_meta - ETOPO metadata (src_w, src_h, factor_x, factor_y, scale_m_per_unit).
   */
  constructor (arg0_spec, arg1_bed_i16, arg2_meta) {
    //Convert from parameters
    let spec = arg0_spec ? arg0_spec : {};
    let bed = arg1_bed_i16;
    let meta = arg2_meta;

    //Declare local instance variables
    let all_domains = spec.domains || [];
    let fx = meta.factor_x;
    let fy = meta.factor_y;
    let scale_m = meta.scale_m_per_unit;
    let src_w = meta.src_w;
    let z_floor = (spec.z_floor_m !== undefined) ? spec.z_floor_m : -30.0;
    let z_full_lo = (spec.z_full_low_m !== undefined) ? spec.z_full_low_m : -18.0;
    let z_full_hi = (spec.z_full_high_m !== undefined) ? spec.z_full_high_m : 2.0;
    let z_top = (spec.z_top_m !== undefined) ? spec.z_top_m : 8.0;
    let edge_fade_km = (spec.edge_fade_km !== undefined) ? spec.edge_fade_km : 25.0;

    this.completion_curve = spec.completion_curve || [[-10000, 1.0]];
    this.domains = [];
    this.max_5m = new Float32Array(4320*2160);
    this.scale_m = scale_m;
    this.src_w = src_w;

    for (let di = 0; di < all_domains.length; di++) {
      let dm = all_domains[di];
      let sy0 = Math.max(0, Math.floor(((90.0 - dm.bbox.lat_max)/180.0)*meta.src_h));
      let sy1 = Math.min(meta.src_h, Math.ceil(((90.0 - dm.bbox.lat_min)/180.0)*meta.src_h));
      let sx0 = Math.max(0, Math.floor(((dm.bbox.lon_min + 180.0)/360.0)*src_w));
      let sx1 = Math.min(src_w, Math.ceil(((dm.bbox.lon_max + 180.0)/360.0)*src_w));
      let h = sy1 - sy0;
      let w = sx1 - sx0;
      let thick_dm = new Uint8Array(w*h);
      let vis = new Uint8Array(w*h);
      let q = new Int32Array(w*h);
      let qh = 0;
      let qt = 0;
      let sr = Math.floor(((90.0 - dm.seed.lat)/180.0)*meta.src_h) - sy0;
      let sc = Math.floor(((dm.seed.lon + 180.0)/360.0)*src_w) - sx0;

      //Guard clause: seed must lie inside its domain
      if (sr < 0 || sr >= h || sc < 0 || sc >= w) {
        console.error(`[HoloceneInfill] Seed outside domain: ${dm.id}`);
        continue;
      }

      let isEligible = (arg0_idx) => {
        let ly = (arg0_idx/w) | 0;
        let z = bed[(sy0 + ly)*src_w + sx0 + (arg0_idx - ly*w)]*scale_m;

        //Return statement
        return (z >= z_floor && z <= z_top);
      };

      //Locate nearest eligible sub-cell to the seed within ~10 km
      let best = -1;
      for (let d = 0; d <= 12 && best < 0; d++)
        for (let dy = -d; dy <= d && best < 0; dy++)
          for (let dx = -d; dx <= d && best < 0; dx++) {
            let rr = sr + dy;
            let cc = sc + dx;
            if (rr < 0 || rr >= h || cc < 0 || cc >= w) continue;
            if (isEligible(rr*w + cc)) best = rr*w + cc;
          }
      if (best < 0) {
        console.error(`[HoloceneInfill] No eligible seed cell: ${dm.id}`);
        continue;
      }
      vis[best] = 1;
      q[qt++] = best;

      //4-connected growth through the estuarine/deltaic elevation band inside the domain
      for (; qh < qt; ) {
        let cur = q[qh++];
        let cy = (cur/w) | 0;
        let cx = cur - cy*w;
        if (cx > 0 && vis[cur - 1] === 0 && isEligible(cur - 1)) {
          vis[cur - 1] = 1; q[qt++] = cur - 1;
        }
        if (cx < w - 1 && vis[cur + 1] === 0 && isEligible(cur + 1)) {
          vis[cur + 1] = 1; q[qt++] = cur + 1;
        }
        if (cy > 0 && vis[cur - w] === 0 && isEligible(cur - w)) {
          vis[cur - w] = 1; q[qt++] = cur - w;
        }
        if (cy < h - 1 && vis[cur + w] === 0 && isEligible(cur + w)) {
          vis[cur + w] = 1; q[qt++] = cur + w;
        }
      }

      //Thickness: full inside the estuary floor / delta plain band, tapering to the top and floor of the band and towards domain edges
      let km_y = 111.2/120.0;
      for (let ly = 0; ly < h; ly++) {
        let lat = 90.0 - (sy0 + ly + 0.5)/120.0;
        let km_x = km_y*Math.max(0.05, Math.cos((lat*Math.PI)/180.0));
        let d_edge_y = Math.min(ly, h - 1 - ly)*km_y;
        let row = (sy0 + ly)*src_w + sx0;

        for (let lx = 0; lx < w; lx++) {
          let li = ly*w + lx;
          if (vis[li] === 0) continue;
          let z = bed[row + lx]*scale_m;
          let fz = 1.0;
          if (z > z_full_hi) {
            fz = Math.max(0.0, (z_top - z)/(z_top - z_full_hi));
          } else if (z < z_full_lo) {
            fz = Math.max(0.0, (z - z_floor)/(z_full_lo - z_floor));
          }
          let d_edge = Math.min(d_edge_y, Math.min(lx, w - 1 - lx)*km_x);
          let fe = Math.min(1.0, d_edge/edge_fade_km);
          let s_m = dm.max_thickness_m*fz*fe;
          thick_dm[li] = Math.min(255, Math.round(s_m*10.0));

          let g_idx = Math.floor((sy0 + ly)/fy)*4320 + Math.floor((sx0 + lx)/fx);
          if (s_m > this.max_5m[g_idx]) this.max_5m[g_idx] = s_m;
        }
      }

      this.domains.push({ h: h, id: dm.id, sx0: sx0, sy0: sy0, thick_dm: thick_dm, w: w });
      console.log(`[HoloceneInfill]   ${dm.id}: ${qt} 30" sub-cells, max ${dm.max_thickness_m} m`);
    }
  }


  /**
   * Returns the fraction of Holocene infill already deposited at target_year (0 = pre-transgression valley, 1 = modern surface).
   *
   * @param {number} arg0_year
   *
   * @returns {number}
   */
  getCompletion (arg0_year) {
    //Convert from parameters
    let year = arg0_year;

    //Declare local instance variables
    let curve = this.completion_curve;

    //Guard clauses
    if (year <= curve[0][0]) return curve[0][1];
    if (year >= curve[curve.length - 1][0]) return curve[curve.length - 1][1];

    for (let i = 0; i < curve.length - 1; i++) {
      if (year >= curve[i][0] && year <= curve[i + 1][0]) {
        let t = (year - curve[i][0])/(curve[i + 1][0] - curve[i][0]);

        //Return statement
        return curve[i][1] + t*(curve[i + 1][1] - curve[i][1]);
      }
    }

    //Return statement
    return 1.0;
  }


  /**
   * Returns the missing infill thickness (m) of a 30" sub-cell at the current completion fraction.
   *
   * @param {number} arg0_sy - Global 30" row.
   * @param {number} arg1_sx - Global 30" column.
   * @param {number} arg2_missing_frac - 1 - completion(t).
   *
   * @returns {number}
   */
  getMissingThicknessM (arg0_sy, arg1_sx, arg2_missing_frac) {
    //Convert from parameters
    let sy = arg0_sy;
    let sx = arg1_sx;
    let missing_frac = arg2_missing_frac;

    //Declare local instance variables
    let s_max = 0;

    for (let i = 0; i < this.domains.length; i++) {
      let dm = this.domains[i];
      let ly = sy - dm.sy0;
      let lx = sx - dm.sx0;
      if (ly < 0 || ly >= dm.h || lx < 0 || lx >= dm.w) continue;
      let v = dm.thick_dm[ly*dm.w + lx];
      if (v > s_max) s_max = v;
    }

    //Return statement
    return s_max*0.1*missing_frac;
  }
}

module.exports = HoloceneInfillModel;
