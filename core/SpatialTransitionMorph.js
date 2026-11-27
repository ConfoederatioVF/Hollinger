"use strict";

/**
 * Geodesic Signed-Distance-Field (SDF) spatial front-advection engine for keyframe transitions.
 * Eliminates translucent temporal cross-fade ghosts along retreating ice sheets, proglacial lakes,
 * and migrating Holocene shorelines by propagating phase boundaries smoothly across space.
 */
class SpatialTransitionMorph {
  /**
   * Maps a signed km2 cell value to its discrete physical phase: +1 land, -1 ice, 0 water.
   *
   * @param {number} arg0_val
   *
   * @returns {number}
   */
  static getPhase (arg0_val) {
    //Convert from parameters
    let val = arg0_val;

    //Return statement
    if (val > 0.0) return 1;
    if (val < 0.0) return -1;
    return 0;
  }

  /**
   * Computes O(M) 8-connected Chamfer distances within a compact subset of grid pixels using Dial's buckets.
   *
   * @param {Int32Array} arg0_group_pixels - Flat grid indices of pixels in this phase-transition group.
   * @param {Int32Array} arg1_grid_to_local - Reusable w*h lookup table mapping grid idx -> local idx (or -1).
   * @param {Int32Array} arg2_seed_dists - Initial distances (10/14 on boundary seeds, INF otherwise).
   * @param {number} arg3_w
   * @param {number} arg4_h
   *
   * @returns {Int32Array}
   */
  static computeGroupDialDistances (arg0_group_pixels, arg1_grid_to_local, arg2_seed_dists, arg3_w, arg4_h) {
    //Convert from parameters
    let group_pixels = arg0_group_pixels;
    let grid_to_local = arg1_grid_to_local;
    let seed_dists = arg2_seed_dists;
    let w = arg3_w;
    let h = arg4_h;

    //Declare local instance variables
    let inf = 0x3fffffff;
    let m = group_pixels.length;
    let dist = new Int32Array(seed_dists);
    let num_buckets = 16;
    let buckets = [];
    let max_d = 0;
    let min_d = inf;

    for (let b = 0; b < num_buckets; b++)
      buckets.push([]);

    for (let li = 0; li < m; li++) {
      let d = dist[li];
      if (d < inf) {
        buckets[d & 15].push(li);
        if (d < min_d) min_d = d;
        if (d > max_d) max_d = d;
      }
    }

    if (min_d === inf) return dist;

    let dx_arr = [-1, 1, 0, 0, -1, 1, -1, 1];
    let dy_arr = [0, 0, -1, 1, -1, -1, 1, 1];
    let step_cost = [10, 10, 10, 10, 14, 14, 14, 14];
    let cur_d = min_d;

    for (; cur_d <= max_d; ) {
      let b_idx = cur_d & 15;
      let bucket = buckets[b_idx];

      if (bucket.length === 0) {
        cur_d++;
        continue;
      }

      let u = bucket.pop();
      if (dist[u] < cur_d) continue;

      let g_idx = group_pixels[u];
      let gy = (g_idx/w) | 0;
      let gx = g_idx - gy*w;

      for (let k = 0; k < 8; k++) {
        let ny = gy + dy_arr[k];
        if (ny < 0 || ny >= h) continue;
        let nx = gx + dx_arr[k];
        if (nx < 0) nx += w;
        else if (nx >= w) nx -= w;

        let n_grid = ny*w + nx;
        let v = grid_to_local[n_grid];
        if (v < 0) continue;

        let nd = cur_d + step_cost[k];
        if (nd < dist[v]) {
          dist[v] = nd;
          buckets[nd & 15].push(v);
          if (nd > max_d) max_d = nd;
        }
      }
    }

    //Return statement
    return dist;
  }

  /**
   * Builds a spatial front-advection morph plan between two landarea keyframes.
   *
   * @param {Float32Array} arg0_left
   * @param {Float32Array} arg1_right
   * @param {number} arg2_w
   * @param {number} arg3_h
   * @param {Object} [arg4_scratch] - Optional reusable scratch { grid_to_local: Int32Array(w*h) }.
   *
   * @returns {Object}
   */
  static buildPlan (arg0_left, arg1_right, arg2_w, arg3_h, arg4_scratch) {
    //Convert from parameters
    let left = arg0_left;
    let right = arg1_right;
    let w = arg2_w;
    let h = arg3_h;
    let scratch = arg4_scratch ? arg4_scratch : {};

    //Declare local instance variables
    let n = w*h;
    if (!scratch.grid_to_local || scratch.grid_to_local.length !== n) {
      scratch.grid_to_local = new Int32Array(n);
      scratch.grid_to_local.fill(-1);
    }
    let grid_to_local = scratch.grid_to_local;
    let inf = 0x3fffffff;
    let groups = [[], [], [], []];
    let same_diff_list = [];

    for (let i = 0; i < n; i++) {
      let a = left[i];
      let b = right[i];
      if (a === b) continue;

      let p0 = this.getPhase(a);
      let p1 = this.getPhase(b);

      if (p0 === p1) {
        same_diff_list.push(i);
      } else if (p0 === -1 && p1 !== -1) {
        groups[0].push(i); //Ice retreat
      } else if (p0 !== -1 && p1 === -1) {
        groups[1].push(i); //Ice advance
      } else if (p0 === 1 && p1 === 0) {
        groups[2].push(i); //Flooding / transgression
      } else if (p0 === 0 && p1 === 1) {
        groups[3].push(i); //Drainage / regression / progradation
      }
    }

    let total_active = groups[0].length + groups[1].length + groups[2].length + groups[3].length;
    let active_half_w = new Float32Array(total_active);
    let active_idx = new Int32Array(total_active);
    let active_sharp = new Uint8Array(total_active);
    let active_tau = new Float32Array(total_active);
    let dx_arr = [-1, 1, 0, 0, -1, 1, -1, 1];
    let dy_arr = [0, 0, -1, 1, -1, -1, 1, 1];
    let step_cost = [10, 10, 10, 10, 14, 14, 14, 14];
    let write_ptr = 0;

    for (let g = 0; g < 4; g++) {
      let g_list = groups[g];
      let m = g_list.length;
      if (m === 0) continue;

      let group_pixels = new Int32Array(g_list);
      for (let li = 0; li < m; li++)
        grid_to_local[group_pixels[li]] = li;

      let seed_start = new Int32Array(m);
      let seed_end = new Int32Array(m);
      seed_start.fill(inf);
      seed_end.fill(inf);

      for (let li = 0; li < m; li++) {
        let g_idx = group_pixels[li];
        let gy = (g_idx/w) | 0;
        let gx = g_idx - gy*w;
        let best_s = inf;
        let best_e = inf;

        for (let k = 0; k < 8; k++) {
          let ny = gy + dy_arr[k];
          if (ny < 0 || ny >= h) continue;
          let nx = gx + dx_arr[k];
          if (nx < 0) nx += w;
          else if (nx >= w) nx -= w;

          let n_grid = ny*w + nx;
          if (grid_to_local[n_grid] >= 0) continue;

          let c = step_cost[k];
          let np0 = this.getPhase(left[n_grid]);
          let np1 = this.getPhase(right[n_grid]);
          let is_start = false;
          let is_end = false;

          if (g === 0) {
            is_start = (np0 !== -1);
            is_end = (np1 === -1);
          } else if (g === 1) {
            is_start = (np0 === -1);
            is_end = (np1 !== -1);
          } else if (g === 2) {
            is_start = (np0 === 0);
            is_end = (np1 === 1);
          } else {
            is_start = (np0 === 1);
            is_end = (np1 === 0);
          }

          if (is_start && c < best_s) best_s = c;
          if (is_end && c < best_e) best_e = c;
        }

        seed_start[li] = best_s;
        seed_end[li] = best_e;
      }

      let dist_start = this.computeGroupDialDistances(group_pixels, grid_to_local, seed_start, w, h);
      let dist_end = this.computeGroupDialDistances(group_pixels, grid_to_local, seed_end, w, h);

      //Fallback maxima per connected component where one or both endpoints are unreachable (isolated appearing/vanishing features)
      let comp_id = new Int32Array(m);
      comp_id.fill(-1);
      let comp_max_s = [];
      let comp_max_e = [];

      for (let li = 0; li < m; li++) {
        if (comp_id[li] >= 0) continue;
        if (dist_start[li] < inf && dist_end[li] < inf) continue;

        let cid = comp_max_s.length;
        let max_s = 0;
        let max_e = 0;
        let queue = [li];
        comp_id[li] = cid;

        for (let qh = 0; qh < queue.length; qh++) {
          let u = queue[qh];
          if (dist_start[u] < inf && dist_start[u] > max_s) max_s = dist_start[u];
          if (dist_end[u] < inf && dist_end[u] > max_e) max_e = dist_end[u];

          let g_idx = group_pixels[u];
          let gy = (g_idx/w) | 0;
          let gx = g_idx - gy*w;

          for (let k = 0; k < 8; k++) {
            let ny = gy + dy_arr[k];
            if (ny < 0 || ny >= h) continue;
            let nx = gx + dx_arr[k];
            if (nx < 0) nx += w;
            else if (nx >= w) nx -= w;

            let v = grid_to_local[ny*w + nx];
            if (v >= 0 && comp_id[v] < 0) {
              comp_id[v] = cid;
              queue.push(v);
            }
          }
        }

        comp_max_s.push(max_s);
        comp_max_e.push(max_e);
      }

      for (let li = 0; li < m; li++) {
        let idx = group_pixels[li];
        let ds = dist_start[li];
        let de = dist_end[li];
        let hw = 0.5;
        let sharp = 0;
        let tau = 0.5;

        if (ds < inf && de < inf) {
          let d_sum = ds + de;
          tau = ds/d_sum;
          hw = Math.min(0.5, 6.5/d_sum);
        } else if (ds < inf) {
          let denom = comp_max_s[comp_id[li]] + 10.0;
          tau = ds/denom;
          hw = Math.min(0.5, 6.5/denom);
        } else if (de < inf) {
          let denom = comp_max_e[comp_id[li]] + 10.0;
          tau = 1.0 - de/denom;
          hw = Math.min(0.5, 6.5/denom);
        }

        if (g === 0 && right[idx] > 0.0) sharp = 1;
        else if (g === 1 && left[idx] > 0.0) sharp = 1;

        active_idx[write_ptr] = idx;
        active_tau[write_ptr] = tau;
        active_half_w[write_ptr] = hw;
        active_sharp[write_ptr] = sharp;
        write_ptr++;
      }

      for (let li = 0; li < m; li++)
        grid_to_local[group_pixels[li]] = -1;
    }

    //Return statement
    return {
      active_half_w: active_half_w,
      active_idx: active_idx,
      active_sharp: active_sharp,
      active_tau: active_tau,
      height: h,
      same_diff_idx: new Int32Array(same_diff_list),
      width: w,
    };
  }

  /**
   * Interpolates a signed 4320x2160 landarea grid using precomputed spatial front advection.
   *
   * @param {Float32Array} arg0_left
   * @param {Float32Array} arg1_right
   * @param {Object} arg2_plan
   * @param {number} arg3_t
   * @param {Float32Array} arg4_out
   *
   * @returns {Float32Array}
   */
  static interpolateLandGrid (arg0_left, arg1_right, arg2_plan, arg3_t, arg4_out) {
    //Convert from parameters
    let left = arg0_left;
    let right = arg1_right;
    let plan = arg2_plan;
    let t = arg3_t;
    let out = arg4_out;

    //Guard clauses
    if (t <= 0.0) {
      out.set(left);
      return out;
    }
    if (t >= 1.0) {
      out.set(right);
      return out;
    }

    //Declare local instance variables
    let active_hw = plan.active_half_w;
    let active_idx = plan.active_idx;
    let active_sharp = plan.active_sharp;
    let active_tau = plan.active_tau;
    let same_idx = plan.same_diff_idx;
    let w0 = 1.0 - t;

    out.set(left);

    for (let i = 0; i < same_idx.length; i++) {
      let idx = same_idx[i];
      out[idx] = w0*left[idx] + t*right[idx];
    }

    for (let i = 0; i < active_idx.length; i++) {
      let idx = active_idx[i];
      let tau = active_tau[i];

      if (active_sharp[i] === 1) {
        out[idx] = (t < tau) ? left[idx] : right[idx];
        continue;
      }

      let hw = active_hw[i];
      if (t <= tau - hw) {
        out[idx] = left[idx];
      } else if (t >= tau + hw) {
        out[idx] = right[idx];
      } else {
        let s = (t - (tau - hw))/(2.0*hw);
        out[idx] = (1.0 - s)*left[idx] + s*right[idx];
      }
    }

    //Return statement
    return out;
  }

  /**
   * Projects a source-grid morph plan onto Equal Earth canvas coordinates once per keyframe interval.
   *
   * @param {Object} arg0_plan
   * @param {Object} arg1_projector
   * @param {Object} [arg2_scratch]
   *
   * @returns {Object}
   */
  static projectPlanToCanvas (arg0_plan, arg1_projector, arg2_scratch) {
    //Convert from parameters
    let plan = arg0_plan;
    let projector = arg1_projector;
    let scratch = arg2_scratch ? arg2_scratch : {};

    //Declare local instance variables
    let active_hw = plan.active_half_w;
    let active_idx = plan.active_idx;
    let active_tau = plan.active_tau;
    let src_h = plan.height;
    let src_w = plan.width;
    let n_src = src_w*src_h;

    if (!scratch.tau_map || scratch.tau_map.length !== n_src) {
      scratch.tau_map = new Float32Array(n_src);
      scratch.tau_map.fill(-1.0);
    }
    if (!scratch.hw_map || scratch.hw_map.length !== n_src) {
      scratch.hw_map = new Float32Array(n_src);
    }

    let hw_map = scratch.hw_map;
    let tau_map = scratch.tau_map;

    for (let i = 0; i < active_idx.length; i++) {
      let idx = active_idx[i];
      tau_map[idx] = active_tau[i];
      hw_map[idx] = active_hw[i];
    }

    let edge_alpha = projector.edge_alpha;
    let n_pix = projector.width*projector.height;
    let out_hw_list = [];
    let out_idx_list = [];
    let out_tau_list = [];
    let src_col = projector.src_col;
    let src_row = projector.src_row;

    for (let p = 0; p < n_pix; p++) {
      if (edge_alpha[p] === 0) continue;

      let r_f = src_row[p];
      let c_f = src_col[p];
      let r0 = r_f | 0;
      let c0 = c_f | 0;
      let r1 = (r0 + 1 < src_h) ? (r0 + 1) : r0;
      let c1 = (c0 + 1 < src_w) ? (c0 + 1) : 0;

      let i00 = r0*src_w + c0;
      let i01 = r0*src_w + c1;
      let i10 = r1*src_w + c0;
      let i11 = r1*src_w + c1;

      let t00 = tau_map[i00];
      let t01 = tau_map[i01];
      let t10 = tau_map[i10];
      let t11 = tau_map[i11];

      if (t00 < 0.0 && t01 < 0.0 && t10 < 0.0 && t11 < 0.0) continue;

      let dr = r_f - r0;
      let dc = c_f - c0;
      let w00 = (1.0 - dr)*(1.0 - dc);
      let w01 = (1.0 - dr)*dc;
      let w10 = dr*(1.0 - dc);
      let w11 = dr*dc;

      let w_sum = 0.0;
      let tau_sum = 0.0;
      let hw_sum = 0.0;
      let cnt = 0;
      let eq_tau = 0.0;
      let eq_hw = 0.0;

      if (t00 >= 0.0) { w_sum += w00; tau_sum += w00*t00; hw_sum += w00*hw_map[i00]; eq_tau += t00; eq_hw += hw_map[i00]; cnt++; }
      if (t01 >= 0.0) { w_sum += w01; tau_sum += w01*t01; hw_sum += w01*hw_map[i01]; eq_tau += t01; eq_hw += hw_map[i01]; cnt++; }
      if (t10 >= 0.0) { w_sum += w10; tau_sum += w10*t10; hw_sum += w10*hw_map[i10]; eq_tau += t10; eq_hw += hw_map[i10]; cnt++; }
      if (t11 >= 0.0) { w_sum += w11; tau_sum += w11*t11; hw_sum += w11*hw_map[i11]; eq_tau += t11; eq_hw += hw_map[i11]; cnt++; }

      let tau_out = (w_sum > 1e-6) ? (tau_sum/w_sum) : (eq_tau/cnt);
      let hw_out = (w_sum > 1e-6) ? (hw_sum/w_sum) : (eq_hw/cnt);

      out_idx_list.push(p);
      out_tau_list.push(tau_out);
      out_hw_list.push(hw_out);
    }

    for (let i = 0; i < active_idx.length; i++)
      tau_map[active_idx[i]] = -1.0;

    //Return statement
    return {
      out_hw: new Float32Array(out_hw_list),
      out_idx: new Int32Array(out_idx_list),
      out_tau: new Float32Array(out_tau_list),
    };
  }

  /**
   * Builds a dual-layer geodesic boundary-morph plan between two categorical keyframes
   * (Layer 1: physical land/ice/water phase advection; Layer 2: target-class front expansion and contraction).
   *
   * @param {Uint8Array|Uint16Array|Uint32Array} arg0_class0 - Start keyframe class ID grid.
   * @param {Uint8Array|Uint16Array|Uint32Array} arg1_class1 - End keyframe class ID grid.
   * @param {Float32Array} arg2_land0 - Start keyframe signed km2 land grid.
   * @param {Float32Array} arg3_land1 - End keyframe signed km2 land grid.
   * @param {number} arg4_w - Source grid width.
   * @param {number} arg5_h - Source grid height.
   * @param {Object} [arg6_scratch] - Optional reusable scratch buffers.
   *
   * @returns {Object} Categorical morph plan.
   */
  static buildCategoricalPlan (arg0_class0, arg1_class1, arg2_land0, arg3_land1, arg4_w, arg5_h, arg6_scratch) {
    //Convert from parameters
    let class0 = arg0_class0;
    let class1 = arg1_class1;
    let land0 = arg2_land0;
    let land1 = arg3_land1;
    let w = arg4_w;
    let h = arg5_h;
    let scratch = arg6_scratch ? arg6_scratch : {};

    //Declare local instance variables
    let land_plan = this.buildPlan(land0, land1, w, h, scratch);
    let n = w*h;
    if (!scratch.grid_to_local || scratch.grid_to_local.length !== n) {
      scratch.grid_to_local = new Int32Array(n);
      scratch.grid_to_local.fill(-1);
    }
    let grid_to_local = scratch.grid_to_local;
    let inf = 0x3fffffff;
    let same_class_coastal = [];
    let target_groups = new Map();
    let total_cat_active = 0;

    for (let i = 0; i < n; i++) {
      let p0 = this.getPhase(land0[i]);
      let p1 = this.getPhase(land1[i]);
      if (p0 !== p1) continue;

      let c0 = class0[i];
      let c1 = class1[i];
      if (c0 !== c1) {
        let list = target_groups.get(c1);
        if (!list) {
          list = [];
          target_groups.set(c1, list);
        }
        list.push(i);
        total_cat_active++;
      } else if (land0[i] !== land1[i]) {
        same_class_coastal.push(i);
      }
    }

    let n_land = land_plan.active_idx.length;
    let total_active = n_land + total_cat_active + same_class_coastal.length;
    let active_group_key = new Int32Array(total_active);
    let active_half_w = new Float32Array(total_active);
    let active_idx = new Int32Array(total_active);
    let active_tau = new Float32Array(total_active);
    let dx_arr = [-1, 1, 0, 0, -1, 1, -1, 1];
    let dy_arr = [0, 0, -1, 1, -1, -1, 1, 1];
    let step_cost = [10, 10, 10, 10, 14, 14, 14, 14];
    let write_ptr = 0;

    //Layer 1: Physical phase transitions (sharp land<->ice, smoothed water coastlines)
    for (let i = 0; i < n_land; i++) {
      let is_sharp = land_plan.active_sharp[i] === 1;
      active_idx[write_ptr] = land_plan.active_idx[i];
      active_tau[write_ptr] = land_plan.active_tau[i];
      active_half_w[write_ptr] = is_sharp ? 0.0 : land_plan.active_half_w[i];
      active_group_key[write_ptr] = is_sharp ? -1 : -2;
      write_ptr++;
    }

    //Layer 2: Same-phase categorical transitions grouped by target class k
    for (let [k, g_list] of target_groups.entries()) {
      let m = g_list.length;
      if (m === 0) continue;

      let group_pixels = new Int32Array(g_list);
      for (let li = 0; li < m; li++)
        grid_to_local[group_pixels[li]] = li;

      let seed_start = new Int32Array(m);
      let seed_end = new Int32Array(m);
      seed_start.fill(inf);
      seed_end.fill(inf);

      for (let li = 0; li < m; li++) {
        let g_idx = group_pixels[li];
        let gy = (g_idx/w) | 0;
        let gx = g_idx - gy*w;
        let best_s = inf;
        let best_e = inf;

        for (let d = 0; d < 8; d++) {
          let ny = gy + dy_arr[d];
          if (ny < 0 || ny >= h) continue;
          let nx = gx + dx_arr[d];
          if (nx < 0) nx += w;
          else if (nx >= w) nx -= w;

          let n_grid = ny*w + nx;
          if (grid_to_local[n_grid] >= 0) continue;

          let c = step_cost[d];
          let nc0 = class0[n_grid];
          let nc1 = class1[n_grid];

          if (nc0 === k && c < best_s) best_s = c;
          if (!(nc0 === 0 && nc1 === 0) && nc1 !== k && c < best_e) best_e = c;
        }

        seed_start[li] = best_s;
        seed_end[li] = best_e;
      }

      let dist_start = this.computeGroupDialDistances(group_pixels, grid_to_local, seed_start, w, h);
      let dist_end = this.computeGroupDialDistances(group_pixels, grid_to_local, seed_end, w, h);

      let comp_id = new Int32Array(m);
      comp_id.fill(-1);
      let comp_max_s = [];
      let comp_max_e = [];
      let has_isolated = false;

      for (let li = 0; li < m; li++) {
        if (comp_id[li] >= 0) continue;
        if (dist_start[li] < inf && dist_end[li] < inf) continue;

        let cid = comp_max_s.length;
        let max_s = 0;
        let max_e = 0;
        let queue = [li];
        comp_id[li] = cid;

        for (let qh = 0; qh < queue.length; qh++) {
          let u = queue[qh];
          if (dist_start[u] < inf && dist_start[u] > max_s) max_s = dist_start[u];
          if (dist_end[u] < inf && dist_end[u] > max_e) max_e = dist_end[u];

          let g_idx = group_pixels[u];
          let gy = (g_idx/w) | 0;
          let gx = g_idx - gy*w;

          for (let d = 0; d < 8; d++) {
            let ny = gy + dy_arr[d];
            if (ny < 0 || ny >= h) continue;
            let nx = gx + dx_arr[d];
            if (nx < 0) nx += w;
            else if (nx >= w) nx -= w;

            let v = grid_to_local[ny*w + nx];
            if (v >= 0 && comp_id[v] < 0) {
              comp_id[v] = cid;
              queue.push(v);
            }
          }
        }

        if (max_s === 0 && max_e === 0) has_isolated = true;
        comp_max_s.push(max_s);
        comp_max_e.push(max_e);
      }

      let comp_max_coast = null;
      let dist_coast = null;

      if (has_isolated) {
        let seed_coast = new Int32Array(m);
        seed_coast.fill(inf);

        for (let li = 0; li < m; li++) {
          if (dist_start[li] < inf || dist_end[li] < inf) continue;
          let g_idx = group_pixels[li];
          let gy = (g_idx/w) | 0;
          let gx = g_idx - gy*w;
          let best_c = inf;

          for (let d = 0; d < 8; d++) {
            let ny = gy + dy_arr[d];
            if (ny < 0 || ny >= h) continue;
            let nx = gx + dx_arr[d];
            if (nx < 0) nx += w;
            else if (nx >= w) nx -= w;

            let n_grid = ny*w + nx;
            if (class0[n_grid] === 0 && class1[n_grid] === 0 && step_cost[d] < best_c)
              best_c = step_cost[d];
          }
          seed_coast[li] = best_c;
        }

        dist_coast = this.computeGroupDialDistances(group_pixels, grid_to_local, seed_coast, w, h);
        comp_max_coast = new Int32Array(comp_max_s.length);

        for (let li = 0; li < m; li++) {
          let cid = comp_id[li];
          if (cid >= 0 && dist_coast[li] < inf && dist_coast[li] > comp_max_coast[cid])
            comp_max_coast[cid] = dist_coast[li];
        }
      }

      for (let li = 0; li < m; li++) {
        let idx = group_pixels[li];
        let gy = (idx/w) | 0;
        let gx = idx - gy*w;
        let ds = dist_start[li];
        let de = dist_end[li];
        let band_width_px = 2.0;
        let tau_geo = 0.5;

        if (ds < inf && de < inf) {
          tau_geo = ds/(ds + de);
          band_width_px = (ds + de)/10.0;
        } else if (ds < inf) {
          let denom = comp_max_s[comp_id[li]] + 10.0;
          tau_geo = ds/denom;
          band_width_px = denom/10.0;
        } else if (de < inf) {
          let denom = comp_max_e[comp_id[li]] + 10.0;
          tau_geo = 1.0 - de/denom;
          band_width_px = denom/10.0;
        } else if (dist_coast && dist_coast[li] < inf) {
          let denom = comp_max_coast[comp_id[li]] + 10.0;
          tau_geo = 1.0 - dist_coast[li]/denom;
          band_width_px = denom/10.0;
        }

        let lon = (gx/w)*6.283185307179586;
        let lat = (gy/h)*3.141592653589793;
        let raw_w = 0.55*Math.sin(lon*43.0 + Math.cos(lat*37.0)*1.5) + 0.45*Math.cos(lon*97.0 - lat*89.0 + 0.8);
        let s_wave = Math.min(1.0, Math.max(-1.0, raw_w*1.45));
        let amp = Math.min(0.42, 0.84/Math.max(2.0, band_width_px));
        let tau = Math.min(0.98, Math.max(0.02, tau_geo + amp*s_wave));

        active_idx[write_ptr] = idx;
        active_tau[write_ptr] = tau;
        active_half_w[write_ptr] = 0.0;
        active_group_key[write_ptr] = k;
        write_ptr++;
      }

      for (let li = 0; li < m; li++)
        grid_to_local[group_pixels[li]] = -1;
    }

    //Same-class fractional coastal adjustments (subtle coastal anti-aliased edge shift only)
    for (let i = 0; i < same_class_coastal.length; i++) {
      active_idx[write_ptr] = same_class_coastal[i];
      active_tau[write_ptr] = 0.5;
      active_half_w[write_ptr] = 0.5;
      active_group_key[write_ptr] = -10;
      write_ptr++;
    }

    //Return statement
    return {
      active_group_key: active_group_key,
      active_half_w: active_half_w,
      active_idx: active_idx,
      active_tau: active_tau,
      height: h,
      width: w,
    };
  }

  /**
   * Projects a dual-layer categorical morph plan onto Equal Earth canvas pixels without cross-boundary bleeding.
   *
   * @param {Object} arg0_plan - Output of buildCategoricalPlan.
   * @param {Object} arg1_projector - EqualEarthProjector instance.
   * @param {Uint8Array} arg2_rgb0 - Start keyframe projected RGB24 buffer.
   * @param {Uint8Array} arg3_rgb1 - End keyframe projected RGB24 buffer.
   * @param {Object} [arg4_scratch] - Optional reusable scratch buffers.
   *
   * @returns {Object} Canvas-space morph plan { out_hw, out_idx, out_tau }.
   */
  static projectCategoricalPlanToCanvas (arg0_plan, arg1_projector, arg2_rgb0, arg3_rgb1, arg4_scratch) {
    //Convert from parameters
    let plan = arg0_plan;
    let projector = arg1_projector;
    let rgb0 = arg2_rgb0;
    let rgb1 = arg3_rgb1;
    let scratch = arg4_scratch ? arg4_scratch : {};

    //Declare local instance variables
    let active_grp = plan.active_group_key;
    let active_hw = plan.active_half_w;
    let active_idx = plan.active_idx;
    let active_tau = plan.active_tau;
    let src_h = plan.height;
    let src_w = plan.width;
    let n_src = src_w*src_h;

    if (!scratch.tau_map || scratch.tau_map.length !== n_src) {
      scratch.tau_map = new Float32Array(n_src);
      scratch.tau_map.fill(-1.0);
    }
    if (!scratch.hw_map || scratch.hw_map.length !== n_src)
      scratch.hw_map = new Float32Array(n_src);
    if (!scratch.grp_map || scratch.grp_map.length !== n_src)
      scratch.grp_map = new Int32Array(n_src);

    let grp_map = scratch.grp_map;
    let hw_map = scratch.hw_map;
    let tau_map = scratch.tau_map;

    for (let i = 0; i < active_idx.length; i++) {
      let idx = active_idx[i];
      tau_map[idx] = active_tau[i];
      hw_map[idx] = active_hw[i];
      grp_map[idx] = active_grp[i];
    }

    let edge_alpha = projector.edge_alpha;
    let n_pix = projector.width*projector.height;
    let out_hw_list = [];
    let out_idx_list = [];
    let out_tau_list = [];
    let src_col = projector.src_col;
    let src_row = projector.src_row;

    for (let p = 0; p < n_pix; p++) {
      if (edge_alpha[p] === 0) continue;

      let p3 = p*3;
      if (rgb0[p3] === rgb1[p3] && rgb0[p3 + 1] === rgb1[p3 + 1] && rgb0[p3 + 2] === rgb1[p3 + 2])
        continue;

      let r_f = src_row[p];
      let c_f = src_col[p];
      let rn = Math.min(src_h - 1, Math.max(0, Math.round(r_f)));
      let cn = Math.round(c_f)%src_w;
      let i_near = rn*src_w + cn;

      let r0 = r_f | 0;
      let c0 = c_f | 0;
      let r1 = (r0 + 1 < src_h) ? (r0 + 1) : r0;
      let c1 = (c0 + 1 < src_w) ? (c0 + 1) : 0;

      let i00 = r0*src_w + c0;
      let i01 = r0*src_w + c1;
      let i10 = r1*src_w + c0;
      let i11 = r1*src_w + c1;

      let dr = r_f - r0;
      let dc = c_f - c0;
      let w00 = (1.0 - dr)*(1.0 - dc);
      let w01 = (1.0 - dr)*dc;
      let w10 = dr*(1.0 - dc);
      let w11 = dr*dc;

      let ref_grp = -999;
      if (tau_map[i_near] >= 0.0) {
        ref_grp = grp_map[i_near];
      } else {
        let best_w = -1.0;
        if (tau_map[i00] >= 0.0 && w00 > best_w) { best_w = w00; ref_grp = grp_map[i00]; }
        if (tau_map[i01] >= 0.0 && w01 > best_w) { best_w = w01; ref_grp = grp_map[i01]; }
        if (tau_map[i10] >= 0.0 && w10 > best_w) { best_w = w10; ref_grp = grp_map[i10]; }
        if (tau_map[i11] >= 0.0 && w11 > best_w) { ref_grp = grp_map[i11]; }
      }

      if (ref_grp === -999) {
        out_idx_list.push(p);
        out_tau_list.push(0.5);
        out_hw_list.push(0.0);
        continue;
      }

      let w_sum = 0.0;
      let tau_sum = 0.0;
      let hw_sum = 0.0;
      let eq_tau = 0.0;
      let eq_hw = 0.0;
      let cnt = 0;

      if (tau_map[i00] >= 0.0 && grp_map[i00] === ref_grp) { w_sum += w00; tau_sum += w00*tau_map[i00]; hw_sum += w00*hw_map[i00]; eq_tau += tau_map[i00]; eq_hw += hw_map[i00]; cnt++; }
      if (tau_map[i01] >= 0.0 && grp_map[i01] === ref_grp) { w_sum += w01; tau_sum += w01*tau_map[i01]; hw_sum += w01*hw_map[i01]; eq_tau += tau_map[i01]; eq_hw += hw_map[i01]; cnt++; }
      if (tau_map[i10] >= 0.0 && grp_map[i10] === ref_grp) { w_sum += w10; tau_sum += w10*tau_map[i10]; hw_sum += w10*hw_map[i10]; eq_tau += tau_map[i10]; eq_hw += hw_map[i10]; cnt++; }
      if (tau_map[i11] >= 0.0 && grp_map[i11] === ref_grp) { w_sum += w11; tau_sum += w11*tau_map[i11]; hw_sum += w11*hw_map[i11]; eq_tau += tau_map[i11]; eq_hw += hw_map[i11]; cnt++; }

      let allow_coast_blend = (
        tau_map[i_near] >= 0.0 &&
        grp_map[i_near] === -2
      );
      let tau_out = (w_sum > 1e-6) ? (tau_sum/w_sum) : (eq_tau/cnt);
      let hw_out = allow_coast_blend
        ? ((w_sum > 1e-6) ? (hw_sum/w_sum) : (eq_hw/cnt))
        : 0.0;

      out_idx_list.push(p);
      out_tau_list.push(tau_out);
      out_hw_list.push(hw_out);
    }

    for (let i = 0; i < active_idx.length; i++)
      tau_map[active_idx[i]] = -1.0;

    //Return statement
    return {
      out_hw: new Float32Array(out_hw_list),
      out_idx: new Int32Array(out_idx_list),
      out_tau: new Float32Array(out_tau_list),
    };
  }

  /**
   * Applies projected spatial front advection directly onto a categorical RGB24 frame buffer.
   *
   * @param {Uint8Array} arg0_rgb0
   * @param {Uint8Array} arg1_rgb1
   * @param {Object} arg2_canvas_plan
   * @param {number} arg3_t
   * @param {Uint8Array} arg4_out_rgb
   */
  static applyCanvasMorph (arg0_rgb0, arg1_rgb1, arg2_canvas_plan, arg3_t, arg4_out_rgb) {
    //Convert from parameters
    let rgb0 = arg0_rgb0;
    let rgb1 = arg1_rgb1;
    let cplan = arg2_canvas_plan;
    let t = arg3_t;
    let out_rgb = arg4_out_rgb;

    //Declare local instance variables
    let hw_arr = cplan.out_hw;
    let idx_arr = cplan.out_idx;
    let tau_arr = cplan.out_tau;

    for (let i = 0; i < idx_arr.length; i++) {
      let p3 = idx_arr[i]*3;
      let tau = tau_arr[i];
      let hw = hw_arr[i];

      if (hw <= 0.0) {
        let src = (t < tau) ? rgb0 : rgb1;
        out_rgb[p3] = src[p3];
        out_rgb[p3 + 1] = src[p3 + 1];
        out_rgb[p3 + 2] = src[p3 + 2];
      } else if (t <= tau - hw) {
        out_rgb[p3] = rgb0[p3];
        out_rgb[p3 + 1] = rgb0[p3 + 1];
        out_rgb[p3 + 2] = rgb0[p3 + 2];
      } else if (t >= tau + hw) {
        out_rgb[p3] = rgb1[p3];
        out_rgb[p3 + 1] = rgb1[p3 + 1];
        out_rgb[p3 + 2] = rgb1[p3 + 2];
      } else {
        let s256 = Math.round(((t - (tau - hw))/(2.0*hw))*256);
        let w0 = 256 - s256;
        out_rgb[p3] = (rgb0[p3]*w0 + rgb1[p3]*s256 + 128) >> 8;
        out_rgb[p3 + 1] = (rgb0[p3 + 1]*w0 + rgb1[p3 + 1]*s256 + 128) >> 8;
        out_rgb[p3 + 2] = (rgb0[p3 + 2]*w0 + rgb1[p3 + 2]*s256 + 128) >> 8;
      }
    }
  }
}

module.exports = SpatialTransitionMorph;
