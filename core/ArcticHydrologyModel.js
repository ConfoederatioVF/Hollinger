"use strict";

/**
 * Arctic Siberian surface hydrology simulation on the 1-arcminute ETOPO 2022 surface for a series of relative sea-level tiers:
 *  1. Priority-Flood D8 drainage and upstream area (Barnes et al., 2014).
 *  2. Channel water area from downstream hydraulic geometry `width_m = width_coeff*Q^width_exp` (Moody & Troutman, 2002).
 *  3. Floodplain inundation: Height Above Nearest Drainage (HAND, Nobre et al., 2011) against the overbank flood stage from Manning's equation at `flood_factor*Q`.
 *  4. Fill-and-spill depression lakes with a steady-state water balance: a depression fills to its spill level when catchment runoff exceeds net lake evaporation over its area, otherwise it holds the equilibrium area, filling its lowest cells first.
 */
class ArcticHydrologyModel {
  /**
   * Runs the hydrology simulation for every sea-level tier and aggregates water area into 5' cell fractions.
   *
   * @param {Object} arg0_spec - arctic_hydrology section of hydro_geomorph_basins.json5.
   * @param {BasinContourEngine} arg1_engine - Engine holding the 30" bedrock DEM.
   *
   * @returns {Object} { cells: [[pixel_idx, frac_tier0, ...]], rsl_levels_m, sectors }
   */
  static build (arg0_spec, arg1_engine) {
    //Convert from parameters
    let spec = arg0_spec;
    let engine = arg1_engine;

    //Guard clause
    if (!spec || !spec.watershed_bbox || !spec.bbox || !spec.rsl_levels_m)
      return { cells: [], rsl_levels_m: [], sectors: [] };

    //Declare local instance variables
    let a_chan = spec.channel_min_acc_km2;
    let bed = engine.bed_i16;
    let disp = spec.bbox;
    let dx8 = [-1, 0, 1, -1, 1, -1, 0, 1];
    let dy8 = [-1, -1, -1, 0, 0, 1, 1, 1];
    let e_inland = (spec.lake_net_evaporation_inland_m_per_yr !== undefined) ? spec.lake_net_evaporation_inland_m_per_yr : spec.lake_net_evaporation_m_per_yr;
    let e_net = spec.lake_net_evaporation_m_per_yr;
    let fade_deg = spec.edge_fade_deg;
    let flood_gain = Math.pow(spec.flood_factor, 0.6) - 1.0;
    let fp_lat_n = (spec.floodplain_lat_north_deg !== undefined) ? spec.floodplain_lat_north_deg : 66.5;
    let fp_lat_s = (spec.floodplain_lat_south_deg !== undefined) ? spec.floodplain_lat_south_deg : 62.5;
    let fp_z_full = (spec.floodplain_elev_full_m !== undefined) ? spec.floodplain_elev_full_m : 35.0;
    let fp_z_zero = (spec.floodplain_elev_zero_m !== undefined) ? spec.floodplain_elev_zero_m : 75.0;
    let lake_min_area = spec.lake_min_area_km2;
    let lake_min_depth = spec.lake_min_depth_m;
    let lk_lat_n = (spec.lake_lat_north_deg !== undefined) ? spec.lake_lat_north_deg : 66.5;
    let lk_lat_s = (spec.lake_lat_south_deg !== undefined) ? spec.lake_lat_south_deg : 62.5;
    let lk_z_full = (spec.lake_elev_full_m !== undefined) ? spec.lake_elev_full_m : 45.0;
    let lk_z_zero = (spec.lake_elev_zero_m !== undefined) ? spec.lake_elev_zero_m : 95.0;
    let manning_n = spec.manning_n;
    let rsl_levels = spec.rsl_levels_m.slice().sort((a, b) => a - b);
    let runoff_m = spec.runoff_m_per_yr;
    let scale_m = engine.scale_m;
    let src_w = engine.src_w;
    let stage_max = spec.flood_stage_max_m;
    let wb = spec.watershed_bbox;

    //1-arcminute surfaces: minimum (sea connectivity through incised channels) and mean (land surface for routing and inundation)
    let r0 = Math.max(0, Math.floor((90.0 - wb.lat_max)*60));
    let r1 = Math.min(10800, Math.ceil((90.0 - wb.lat_min)*60));
    let c0 = Math.max(0, Math.floor((wb.lon_min + 180.0)*60));
    let c1 = Math.min(21600, Math.ceil((wb.lon_max + 180.0)*60));
    let h = r1 - r0;
    let w = c1 - c0;
    let n = w*h;
    let z_mean = new Float32Array(n);
    let z_min = new Float32Array(n);

    for (let y = 0; y < h; y++) {
      let row_a = (r0 + y)*2*src_w;
      let row_b = row_a + src_w;
      for (let x = 0; x < w; x++) {
        let sc = (c0 + x)*2;
        let v0 = bed[row_a + sc];
        let v1 = bed[row_a + sc + 1];
        let v2 = bed[row_b + sc];
        let v3 = bed[row_b + sc + 1];
        z_min[y*w + x] = Math.min(v0, v1, v2, v3)*scale_m;
        z_mean[y*w + x] = (v0 + v1 + v2 + v3)*0.25*scale_m;
      }
    }

    let cos_row = new Float32Array(h);
    let km2_row = new Float32Array(h);
    let lat_row = new Float32Array(h);
    for (let y = 0; y < h; y++) {
      lat_row[y] = wb.lat_max - (y + 0.5)/60.0;
      cos_row[y] = Math.max(0.15, Math.cos((lat_row[y]*Math.PI)/180.0));
      km2_row[y] = cos_row[y]*1.853*1.853;
    }

    let dr0 = Math.max(0, Math.floor((90.0 - disp.lat_max)*12));
    let dr1 = Math.min(2160, Math.ceil((90.0 - disp.lat_min)*12));
    let dc0 = Math.max(0, Math.floor((disp.lon_min + 180.0)*12));
    let dc1 = Math.min(4320, Math.ceil((disp.lon_max + 180.0)*12));
    let out_h = dr1 - dr0;
    let out_w = dc1 - dc0;
    let cell5_km2 = new Float32Array(out_h);
    let fade5 = new Float32Array(out_h*out_w);
    for (let ly = 0; ly < out_h; ly++) {
      let lat = 90.0 - (dr0 + ly + 0.5)/12.0;
      cell5_km2[ly] = 9.266*9.266*Math.max(0.05, Math.cos((lat*Math.PI)/180.0));
      for (let lx = 0; lx < out_w; lx++) {
        let lon = -180.0 + (dc0 + lx + 0.5)/12.0;
        let d_edge = Math.min(lat - disp.lat_min, disp.lat_max - lat, lon - disp.lon_min, disp.lon_max - lon);
        fade5[ly*out_w + lx] = Math.min(1.0, Math.max(0.0, d_edge/fade_deg));
      }
    }

    //Reusable buffers across tiers
    let acc = new Float32Array(n);
    let drain_z = new Float32Array(n);
    let fill = new Float32Array(n);
    let h_idx = new Int32Array(n);
    let h_val = new Float32Array(n);
    let label = new Int32Array(n);
    let order = new Int32Array(n);
    let parent = new Int32Array(n);
    let sea = new Uint8Array(n);
    let sea_q = new Int32Array(n);
    let stage = new Float32Array(n);
    let vis = new Uint8Array(n);
    let water_km2 = new Float32Array(out_h*out_w);
    let cell_map = new Map();

    for (let ti = 0; ti < rsl_levels.length; ti++) {
      let rsl_m = rsl_levels[ti];
      label.fill(-1);
      parent.fill(-1);
      sea.fill(0);
      vis.fill(0);
      water_km2.fill(0.0);

      //1. Ocean connected from the northern rim and open domain edges below sea level
      let sqh = 0;
      let sqt = 0;
      for (let y = 0; y < h; y++)
        for (let x = 0; x < w; x++) {
          let i = y*w + x;
          let is_edge = (y === 0 || y === h - 1 || x === 0 || x === w - 1);
          if (is_edge && z_min[i] <= rsl_m) {
            sea[i] = 1;
            sea_q[sqt++] = i;
          }
        }
      for (; sqh < sqt; ) {
        let cur = sea_q[sqh++];
        let cy = (cur/w) | 0;
        let cx = cur - cy*w;
        if (cx > 0 && sea[cur - 1] === 0 && z_min[cur - 1] <= rsl_m) { sea[cur - 1] = 1; sea_q[sqt++] = cur - 1; }
        if (cx < w - 1 && sea[cur + 1] === 0 && z_min[cur + 1] <= rsl_m) { sea[cur + 1] = 1; sea_q[sqt++] = cur + 1; }
        if (cy > 0 && sea[cur - w] === 0 && z_min[cur - w] <= rsl_m) { sea[cur - w] = 1; sea_q[sqt++] = cur - w; }
        if (cy < h - 1 && sea[cur + w] === 0 && z_min[cur + w] <= rsl_m) { sea[cur + w] = 1; sea_q[sqt++] = cur + w; }
      }

      //2. Priority-Flood spill surface and D8 receivers
      let hs = 0;
      let pushHeap = (arg0_idx, arg1_val) => {
        let k = hs++;
        for (; k > 0; ) {
          let p = (k - 1) >> 1;
          if (h_val[p] <= arg1_val) break;
          h_idx[k] = h_idx[p];
          h_val[k] = h_val[p];
          k = p;
        }
        h_idx[k] = arg0_idx;
        h_val[k] = arg1_val;
      };
      let popHeap = () => {
        let top_idx = h_idx[0];
        let top_val = h_val[0];
        hs--;
        if (hs > 0) {
          let idx = h_idx[hs];
          let val = h_val[hs];
          let k = 0;
          for (;;) {
            let l = 2*k + 1;
            let r = l + 1;
            let s = k;
            if (l < hs && h_val[l] < val) s = l;
            if (r < hs && h_val[r] < ((s === k) ? val : h_val[l])) s = r;
            if (s === k) break;
            h_idx[k] = h_idx[s];
            h_val[k] = h_val[s];
            k = s;
          }
          h_idx[k] = idx;
          h_val[k] = val;
        }

        //Return statement
        return [top_idx, top_val];
      };

      for (let i = 0; i < n; i++) {
        if (sea[i] === 1) {
          vis[i] = 1;
          pushHeap(i, rsl_m);
        }
      }

      let ord_n = 0;
      for (; hs > 0; ) {
        let [cur, cz] = popHeap();
        order[ord_n++] = cur;
        fill[cur] = cz;
        let cy = (cur/w) | 0;
        let cx = cur - cy*w;
        for (let d = 0; d < 8; d++) {
          let ny = cy + dy8[d];
          let nx = cx + dx8[d];
          if (ny < 0 || ny >= h || nx < 0 || nx >= w) continue;
          let ni = ny*w + nx;
          if (vis[ni] === 1) continue;
          vis[ni] = 1;
          parent[ni] = cur;
          pushHeap(ni, Math.max(z_mean[ni], cz + 1.0e-5));
        }
      }

      for (let y = 0; y < h; y++)
        for (let x = 0; x < w; x++)
          acc[y*w + x] = (sea[y*w + x] === 1) ? 0.0 : km2_row[y];
      for (let k = ord_n - 1; k >= 0; k--) {
        let c = order[k];
        let p = parent[c];
        if (p >= 0 && sea[c] === 0) acc[p] += acc[c];
      }

      let addWater = (arg0_i, arg1_km2) => {
        let cy = (arg0_i/w) | 0;
        let gy = Math.floor((r0 + cy)/5);
        let gx = Math.floor((c0 + arg0_i - cy*w)/5);
        if (gy < dr0 || gy >= dr1 || gx < dc0 || gx >= dc1) return;
        water_km2[(gy - dr0)*out_w + (gx - dc0)] += arg1_km2;
      };

      //3. Channels (hydraulic-geometry width) and HAND floodplain inundation at the Manning overbank stage
      let water_class = new Uint8Array(n);
      for (let k = 0; k < ord_n; k++) {
        let i = order[k];
        if (sea[i] === 1) continue;
        let p = parent[i];
        let cy = (i/w) | 0;

        if (acc[i] >= a_chan) {
          let q_m3s = (acc[i]*1.0e6*runoff_m)/31557600.0;
          let width_m = spec.width_coeff*Math.pow(q_m3s, spec.width_exp);
          let j = i;
          for (let s = 0; s < spec.slope_reaches && parent[j] >= 0 && sea[parent[j]] === 0; s++) j = parent[j];
          let jy = (j/w) | 0;
          let dist_m = Math.max(1000.0, Math.hypot((jy - cy)*1853.0, (j - jy*w - (i - cy*w))*1853.0*cos_row[cy]));
          let slope = Math.max(spec.min_slope, (fill[i] - fill[j])/dist_m);
          let depth_m = Math.pow((q_m3s*manning_n)/(width_m*Math.sqrt(slope)), 0.6);
          stage[i] = Math.min(stage_max, depth_m*flood_gain);
          drain_z[i] = fill[i];
          water_class[i] = 1;

          let reach_km = 1.853;
          if (p >= 0) {
            let py = (p/w) | 0;
            reach_km = Math.hypot((py - cy)*1.853, (p - py*w - (i - cy*w))*1.853*cos_row[cy]);
          }
          addWater(i, Math.min(km2_row[cy], (width_m/1000.0)*reach_km));
        } else if (p >= 0 && sea[p] === 0 && (water_class[p] === 1 || water_class[p] === 2)) {
          drain_z[i] = drain_z[p];
          stage[i] = stage[p];
          let w_lat = this.smoothstep(fp_lat_s, fp_lat_n, lat_row[cy]);
          let w_z = 1.0 - this.smoothstep(fp_z_full, fp_z_zero, z_mean[i]);
          let stage_eff = stage[i]*w_lat*w_z;
          let is_pit = (fill[i] - z_mean[i] >= lake_min_depth);
          if (!is_pit && z_min[i] > 0.0 && stage_eff > 0.05 && fill[i] - drain_z[i] <= stage_eff) {
            water_class[i] = 2;
            addWater(i, km2_row[cy]);
          }
        }
      }

      //4. Fill-and-spill depression lakes with steady-state climate-graded water balance
      let n_lakes = 0;
      for (let i0 = 0; i0 < n; i0++) {
        if (sea[i0] === 1 || label[i0] !== -1 || fill[i0] - z_mean[i0] < lake_min_depth) continue;
        let area = 0.0;
        let catch_km2 = 0.0;
        let cells = [i0];
        let lat_sum = 0.0;
        let z_sum = 0.0;
        label[i0] = i0;
        for (let t = 0; t < cells.length; t++) {
          let cur = cells[t];
          let cy = (cur/w) | 0;
          let cx = cur - cy*w;
          area += km2_row[cy];
          if (acc[cur] > catch_km2) catch_km2 = acc[cur];
          lat_sum += lat_row[cy];
          z_sum += z_mean[cur];
          for (let d = 0; d < 8; d++) {
            let ny = cy + dy8[d];
            let nx = cx + dx8[d];
            if (ny < 0 || ny >= h || nx < 0 || nx >= w) continue;
            let ni = ny*w + nx;
            if (label[ni] === -1 && sea[ni] === 0 && fill[ni] - z_mean[ni] >= lake_min_depth) {
              label[ni] = i0;
              cells.push(ni);
            }
          }
        }
        if (area < lake_min_area) continue;

        let mean_lat = lat_sum/cells.length;
        let mean_z = z_sum/cells.length;
        let w_lat = this.smoothstep(lk_lat_s, lk_lat_n, mean_lat);
        let w_z = 1.0 - this.smoothstep(lk_z_full, lk_z_zero, mean_z);
        let w_lake = w_lat*w_z;
        if (w_lake <= 0.005) continue;

        let e_eff = e_net + (e_inland - e_net)*(1.0 - w_lat);
        let a_eq = ((catch_km2*runoff_m)/e_eff)*w_lake;
        let a_cap = Math.min(a_eq, area*w_lake);
        if (a_cap < area) cells.sort((a, b) => z_mean[a] - z_mean[b]);
        let filled = 0.0;
        for (let t = 0; t < cells.length; t++) {
          let c = cells[t];
          let cy = (c/w) | 0;
          if (filled >= a_cap) break;
          filled += km2_row[cy];
          if (water_class[c] === 0) addWater(c, km2_row[cy]);
          water_class[c] = 3;
        }
        n_lakes++;
      }

      //5. Store 5' water fractions for this tier
      let n_cells = 0;
      for (let ly = 0; ly < out_h; ly++) {
        for (let lx = 0; lx < out_w; lx++) {
          let l5 = ly*out_w + lx;
          let f = Math.min(1.0, water_km2[l5]/cell5_km2[ly])*fade5[l5];
          if (f < 0.005) continue;
          let p_idx = (dr0 + ly)*4320 + dc0 + lx;
          if (!cell_map.has(p_idx)) cell_map.set(p_idx, new Float32Array(rsl_levels.length));
          cell_map.get(p_idx)[ti] = f;
          n_cells++;
        }
      }
      console.log(`[ArcticHydrology]   RSL=${rsl_m} m -> ${n_cells} 5' water cells, ${n_lakes} lakes`);
    }

    let cells = [];
    for (let [p_idx, fr] of cell_map) {
      let row = [p_idx];
      for (let t = 0; t < fr.length; t++) row.push(Math.round(fr[t]*1000)/1000);
      cells.push(row);
    }

    //Return statement
    return {
      cells: cells,
      rsl_levels_m: rsl_levels,
      sectors: spec.sectors || [],
    };
  }


  /**
   * Cubic Hermite smoothstep in [0, 1] between edge0 and edge1.
   *
   * @param {number} arg0_edge0
   * @param {number} arg1_edge1
   * @param {number} arg2_val
   *
   * @returns {number}
   */
  static smoothstep (arg0_edge0, arg1_edge1, arg2_val) {
    //Convert from parameters
    let edge0 = arg0_edge0;
    let edge1 = arg1_edge1;
    let val = arg2_val;

    //Guard clause
    if (edge1 === edge0) return (val >= edge1) ? 1.0 : 0.0;

    //Declare local instance variables
    let t = Math.min(1.0, Math.max(0.0, (val - edge0)/(edge1 - edge0)));

    //Return statement
    return t*t*(3.0 - 2.0*t);
  }
}

module.exports = ArcticHydrologyModel;
