"use strict";

let fs = require("node:fs");
let path = require("node:path");

let BasinContourEngine = require("../lib/BasinContourEngine");
let GeoGridUtils = require("../lib/GeoGridUtils");
let HoloceneInfillModel = require("../lib/HoloceneInfillModel");
let NetCDF4Writer = require("../lib/NetCDF4Writer");

/**
 * Step C: Builds 100-year paleo GIA relative sea-level coastal anomalies (delta_f_coast)
 * and paleo glaciation + permanent ice-shelf fractions (f_ice_paleo) from PMIP4 GLAC-1D
 * coupled with the 30" ETOPO 2022 bedrock DEM and connected-component marine flood-fill.
 */
class C_BuildPaleoGiaIce {
  /**
   * Returns a smooth Lambeck et al. (2014) / Grant et al. (2012) eustatic calibration offset (m) added to
   * GLAC-1D Δz during 10000–6000 BC so the 12 ka BP stand reaches ~-66 to -68 m in the Persian Gulf and
   * ~-73 to -81 m on the Siberian Arctic shelf (bridging Beringia, unifying Sundaland, drying the Persian Gulf
   * cleanly to ~52E without NW sill extrusions, and enclosing Lake Carpentaria) while returning to 0 m by 6000 BC.
   *
   * @param {number} arg0_year
   *
   * @returns {number}
   */
  static getCalibratedEustaticOffsetM (arg0_year) {
    //Convert from parameters
    let year = arg0_year;

    //Declare local instance variables
    let knots = [
      [-10000, 27.0],
      [-9000, 16.0],
      [-8000, 7.5],
      [-7000, 2.0],
      [-6000, 0.0],
    ];

    //Guard clauses
    if (year <= knots[0][0]) return knots[0][1];
    if (year >= knots[knots.length - 1][0]) return 0.0;

    for (let i = 0; i < knots.length - 1; i++) {
      let y0 = knots[i][0];
      let y1 = knots[i + 1][0];
      if (year >= y0 && year <= y1) {
        let t = (year - y0)/(y1 - y0);
        return knots[i][1] + t*(knots[i + 1][1] - knots[i][1]);
      }
    }

    //Return statement
    return 0.0;
  }

  /**
   * Tests whether a coordinate lies on Greenland proper (matching Step B).
   *
   * @param {number} arg0_lat_deg
   * @param {number} arg1_lon_deg
   *
   * @returns {boolean}
   */
  static isGreenlandProper (arg0_lat_deg, arg1_lon_deg) {
    //Convert from parameters
    let lat = arg0_lat_deg;
    let lon = arg1_lon_deg;

    //Guard clauses
    if (lat < 59.0 || lat > 84.0 || lon > -10.5 || lon < -75.0) return false;
    if (lat < 67.2 && lon > -26.0) return false;

    let min_lon = -58.0;
    if (lat >= 81.0) {
      min_lon = -62.0;
    } else if (lat >= 78.0) {
      min_lon = -62.0 - (81.0 - lat)*(12.0/3.0);
    } else if (lat >= 74.0) {
      min_lon = -74.0 + (78.0 - lat)*(4.0/4.0);
    } else if (lat >= 66.0) {
      min_lon = -70.0 + (74.0 - lat)*(12.0/8.0);
    }

    //Return statement
    return lon >= min_lon;
  }

  /**
   * Runs a local 4-neighbour BFS inside a bounding box on the 4320x2160 grid to keep only
   * paleo ice cells that are contiguously connected to grounded land/core ice across shallow shelf cells.
   *
   * @param {Float32Array} arg0_f_ice_t
   * @param {Float32Array} arg1_f_ice_2025
   * @param {Float32Array} arg2_f_land_osm_2025
   * @param {Float32Array} arg3_z_5m
   * @param {Object} arg4_opts
   */
  static clipIceToConnectedShelf (arg0_f_ice_t, arg1_f_ice_2025, arg2_f_land_osm_2025, arg3_z_5m, arg4_opts) {
    //Convert from parameters
    let f_ice_t = arg0_f_ice_t;
    let f_ice_2025 = arg1_f_ice_2025;
    let f_osm = arg2_f_land_osm_2025;
    let z_5m = arg3_z_5m;
    let opts = arg4_opts;

    //Declare local instance variables
    let width = 4320;
    let height = 2160;
    let r0 = Math.max(0, Math.floor(((90.0 - opts.lat_max)/180.0)*height));
    let r1 = Math.min(height - 1, Math.ceil(((90.0 - opts.lat_min)/180.0)*height));
    let c0 = Math.max(0, Math.floor(((opts.lon_min + 180.0)/360.0)*width));
    let c1 = Math.min(width - 1, Math.ceil(((opts.lon_max + 180.0)/360.0)*width));
    let lh = r1 - r0 + 1;
    let lw = c1 - c0 + 1;
    let max_depth = opts.max_shelf_depth_m;
    let seed_ice_thresh = opts.seed_ice_thresh;
    let wrap_lon = Boolean(opts.wrap_lon);
    let vis = new Uint8Array(lw*lh);
    let q = new Int32Array(lw*lh);
    let qh = 0;
    let qt = 0;

    for (let ly = 0; ly < lh; ly++) {
      let gy = r0 + ly;
      let g_off = gy*width;
      let l_off = ly*lw;

      for (let lx = 0; lx < lw; lx++) {
        let g_idx = g_off + (c0 + lx);
        let is_seed = (f_ice_2025[g_idx] >= seed_ice_thresh) || (opts.seed_from_land && f_osm[g_idx] >= 0.25 && f_ice_t[g_idx] > 0.05);
        if (is_seed) {
          let l_idx = l_off + lx;
          vis[l_idx] = 1;
          q[qt++] = l_idx;
        }
      }
    }

    for (; qh < qt; ) {
      let cur = q[qh++];
      let ly = (cur/lw) | 0;
      let lx = cur - ly*lw;
      let gy = r0 + ly;
      let gx = c0 + lx;

      let tryPush = (arg0_ny, arg1_nx) => {
        let ny = arg0_ny;
        let nx = arg1_nx;
        if (ny < 0 || ny >= lh) return;
        if (nx < 0) {
          if (!wrap_lon) return;
          nx = lw - 1;
        } else if (nx >= lw) {
          if (!wrap_lon) return;
          nx = 0;
        }
        let n_lidx = ny*lw + nx;
        if (vis[n_lidx] === 1) return;
        let n_gidx = (r0 + ny)*width + (c0 + nx);
        if (f_ice_t[n_gidx] <= opts.min_ice_thresh) return;
        if (f_ice_2025[n_gidx] > 0.01 || z_5m[n_gidx] >= max_depth) {
          vis[n_lidx] = 1;
          q[qt++] = n_lidx;
        }
      };

      tryPush(ly, lx - 1);
      tryPush(ly, lx + 1);
      tryPush(ly - 1, lx);
      tryPush(ly + 1, lx);
    }

    for (let ly = 0; ly < lh; ly++) {
      let g_off = (r0 + ly)*width + c0;
      let l_off = ly*lw;
      for (let lx = 0; lx < lw; lx++) {
        if (vis[l_off + lx] === 0) {
          let g_idx = g_off + lx;
          f_ice_t[g_idx] = f_ice_2025[g_idx];
        }
      }
    }
  }

  /**
   * Flood-fills exterior non-ice cells from the open boundary of a domain so any marine lagoon,
   * channel, or inland lake pocket completely enclosed by an ice sheet is sealed to solid f_ice = 1.0.
   *
   * @param {Float32Array} arg0_f_ice_t
   * @param {Float32Array} arg1_f_land_osm_2025
   * @param {Object} arg2_opts
   */
  static fillEnclosedIceHoles (arg0_f_ice_t, arg1_f_land_osm_2025, arg2_opts) {
    //Convert from parameters
    let f_ice_t = arg0_f_ice_t;
    let f_osm = arg1_f_land_osm_2025;
    let opts = arg2_opts;

    //Declare local instance variables
    let height = 2160;
    let width = 4320;
    let r0 = Math.max(0, Math.floor(((90.0 - opts.lat_max)/180.0)*height));
    let r1 = Math.min(height - 1, Math.ceil(((90.0 - opts.lat_min)/180.0)*height));
    let c0 = Math.max(0, Math.floor(((opts.lon_min + 180.0)/360.0)*width));
    let c1 = Math.min(width - 1, Math.ceil(((opts.lon_max + 180.0)/360.0)*width));
    let lh = r1 - r0 + 1;
    let lw = c1 - c0 + 1;
    let marine_only = Boolean(opts.marine_only);
    let seed_all_edges = Boolean(opts.seed_all_edges);
    let wrap_lon = Boolean(opts.wrap_lon);
    let vis = new Uint8Array(lw*lh);
    let q = new Int32Array(lw*lh);
    let qh = 0;
    let qt = 0;

    let pushIfOpen = (arg0_ly, arg1_lx) => {
      let li = arg0_ly*lw + arg1_lx;
      if (vis[li] === 1) return;
      let g_idx = (r0 + arg0_ly)*width + (c0 + arg1_lx);
      if (f_ice_t[g_idx] >= 0.50) return;
      vis[li] = 1;
      q[qt++] = li;
    };

    //Always seed northern boundary (open Southern Ocean for Antarctica; open Arctic edge for NH)
    for (let lx = 0; lx < lw; lx++)
      pushIfOpen(0, lx);

    if (seed_all_edges) {
      for (let lx = 0; lx < lw; lx++)
        pushIfOpen(lh - 1, lx);
      for (let ly = 0; ly < lh; ly++) {
        pushIfOpen(ly, 0);
        pushIfOpen(ly, lw - 1);
      }
    }

    for (; qh < qt; ) {
      let cur = q[qh++];
      let ly = (cur/lw) | 0;
      let lx = cur - ly*lw;

      let tryStep = (arg0_ny, arg1_nx) => {
        let ny = arg0_ny;
        let nx = arg1_nx;
        if (ny < 0 || ny >= lh) return;
        if (nx < 0) {
          if (!wrap_lon) return;
          nx = lw - 1;
        } else if (nx >= lw) {
          if (!wrap_lon) return;
          nx = 0;
        }
        let nli = ny*lw + nx;
        if (vis[nli] === 1) return;
        let g_idx = (r0 + ny)*width + (c0 + nx);
        if (f_ice_t[g_idx] >= 0.50) return;
        vis[nli] = 1;
        q[qt++] = nli;
      };

      tryStep(ly, lx - 1);
      tryStep(ly, lx + 1);
      tryStep(ly - 1, lx);
      tryStep(ly + 1, lx);
    }

    //Any non-ice cell unreachable from the open exterior is an enclosed pocket/slit -> fill to solid ice!
    for (let ly = 0; ly < lh; ly++) {
      let g_off = (r0 + ly)*width + c0;
      let l_off = ly*lw;
      for (let lx = 0; lx < lw; lx++) {
        if (vis[l_off + lx] === 1) continue;
        let g_idx = g_off + lx;
        if (f_ice_t[g_idx] >= 0.50) continue;
        if (marine_only && f_osm[g_idx] >= 0.35) continue;
        f_ice_t[g_idx] = 1.0;
      }
    }
  }

  /**
   * Drops tiny isolated detached paleo-only ice blobs inside a region using connected-component
   * area filtering (replaces rectangular lat/lon box cuts with 100% organic ice-sheet edges).
   *
   * @param {Float32Array} arg0_f_ice_t
   * @param {Float32Array} arg1_f_ice_2025
   * @param {Object} arg2_opts
   */
  static removeIsolatedPaleoIceBlobs (arg0_f_ice_t, arg1_f_ice_2025, arg2_opts) {
    //Convert from parameters
    let f_ice_t = arg0_f_ice_t;
    let f_ice_2025 = arg1_f_ice_2025;
    let opts = arg2_opts;

    //Declare local instance variables
    let height = 2160;
    let width = 4320;
    let r0 = Math.max(0, Math.floor(((90.0 - opts.lat_max)/180.0)*height));
    let r1 = Math.min(height - 1, Math.ceil(((90.0 - opts.lat_min)/180.0)*height));
    let c0 = Math.max(0, Math.floor(((opts.lon_min + 180.0)/360.0)*width));
    let c1 = Math.min(width - 1, Math.ceil(((opts.lon_max + 180.0)/360.0)*width));
    let lh = r1 - r0 + 1;
    let lw = c1 - c0 + 1;
    let min_px = opts.min_component_px;
    let vis = new Uint8Array(lw*lh);
    let q = new Int32Array(lw*lh);

    for (let ly = 0; ly < lh; ly++) {
      for (let lx = 0; lx < lw; lx++) {
        let start = ly*lw + lx;
        if (vis[start] === 1) continue;
        let g_start = (r0 + ly)*width + (c0 + lx);
        if (f_ice_t[g_start] < 0.50) continue;

        let qh = 0;
        let qt = 0;
        let has_modern = false;
        let touches_border = false;
        vis[start] = 1;
        q[qt++] = start;

        for (; qh < qt; ) {
          let cur = q[qh++];
          let cy = (cur/lw) | 0;
          let cx = cur - cy*lw;
          let gi = (r0 + cy)*width + (c0 + cx);
          if (f_ice_2025[gi] >= 0.25) has_modern = true;
          if (cy === 0 || cy === lh - 1 || cx === 0 || cx === lw - 1) touches_border = true;

          let tryN = (arg0_ny, arg1_nx) => {
            if (arg0_ny < 0 || arg0_ny >= lh || arg1_nx < 0 || arg1_nx >= lw) return;
            let nli = arg0_ny*lw + arg1_nx;
            if (vis[nli] === 1) return;
            let ngi = (r0 + arg0_ny)*width + (c0 + arg1_nx);
            if (f_ice_t[ngi] < 0.50) return;
            vis[nli] = 1;
            q[qt++] = nli;
          };

          tryN(cy, cx - 1);
          tryN(cy, cx + 1);
          tryN(cy - 1, cx);
          tryN(cy + 1, cx);
        }

        if (!has_modern && !touches_border && qt < min_px) {
          for (let k = 0; k < qt; k++) {
            let li = q[k];
            let cy = (li/lw) | 0;
            let cx = li - cy*lw;
            let gi = (r0 + cy)*width + (c0 + cx);
            f_ice_t[gi] = f_ice_2025[gi];
          }
        }
      }
    }
  }

  /**
   * Applies 1-pass 3x3 morphological majority smoothing on paleo-only ice margins (leaving modern
   * f_ice_2025 untouched) so coarse GLAC-1D boundaries become smooth, crisp organic contours.
   *
   * @param {Float32Array} arg0_f_ice_t
   * @param {Float32Array} arg1_f_ice_2025
   * @param {Object} arg2_box
   */
  static smoothPaleoIceMargins (arg0_f_ice_t, arg1_f_ice_2025, arg2_box) {
    //Convert from parameters
    let f_ice_t = arg0_f_ice_t;
    let f_ice_2025 = arg1_f_ice_2025;
    let box = arg2_box;

    //Declare local instance variables
    let height = 2160;
    let width = 4320;
    let r0 = Math.max(1, Math.floor(((90.0 - box.lat_max)/180.0)*height));
    let r1 = Math.min(height - 2, Math.ceil(((90.0 - box.lat_min)/180.0)*height));
    let c0 = Math.max(1, Math.floor(((box.lon_min + 180.0)/360.0)*width));
    let c1 = Math.min(width - 2, Math.ceil(((box.lon_max + 180.0)/360.0)*width));
    let lh = r1 - r0 + 1;
    let lw = c1 - c0 + 1;
    let out = new Float32Array(lw*lh);

    for (let ly = 0; ly < lh; ly++) {
      let gy = r0 + ly;
      let g_off = gy*width;
      let l_off = ly*lw;

      for (let lx = 0; lx < lw; lx++) {
        let gx = c0 + lx;
        let gi = g_off + gx;
        if (f_ice_2025[gi] > 0.01) {
          out[l_off + lx] = f_ice_t[gi];
          continue;
        }
        let cnt = 0;
        for (let dy = -1; dy <= 1; dy++) {
          let r_i = (gy + dy)*width + gx;
          for (let dx = -1; dx <= 1; dx++)
            if (f_ice_t[r_i + dx] >= 0.50) cnt++;
        }
        out[l_off + lx] = (cnt >= 5) ? 1.0 : 0.0;
      }
    }

    for (let ly = 0; ly < lh; ly++) {
      let g_off = (r0 + ly)*width + c0;
      let l_off = ly*lw;
      for (let lx = 0; lx < lw; lx++)
        f_ice_t[g_off + lx] = out[l_off + lx];
    }
  }

  /**
   * Performs morphological opening (N passes of 4-neighbour erosion followed by N passes of dilation)
   * on paleo-only ice inside a southern polar sector (wrapping longitude across 180 deg) to prune
   * thin protruding ice claws/hooks while keeping modern f_ice_2025 strictly intact.
   *
   * @param {Float32Array} arg0_f_ice_t
   * @param {Float32Array} arg1_f_ice_2025
   * @param {Object} arg2_opts
   */
  static openPaleoIceAntimeridian (arg0_f_ice_t, arg1_f_ice_2025, arg2_opts) {
    //Convert from parameters
    let f_ice_t = arg0_f_ice_t;
    let f_ice_2025 = arg1_f_ice_2025;
    let opts = arg2_opts;

    //Declare local instance variables
    let height = 2160;
    let width = 4320;
    let r0 = Math.max(1, Math.floor(((90.0 - opts.lat_max)/180.0)*height));
    let r1 = Math.min(height - 2, Math.ceil(((90.0 - opts.lat_min)/180.0)*height));
    let lh = r1 - r0 + 1;
    let passes = opts.passes;
    //Build wrapped column list from lon_west (e.g. +165) east across antimeridian to lon_east (e.g. -150)
    let c_west = Math.floor(((opts.lon_west + 180.0)/360.0)*width);
    let c_east = Math.ceil(((opts.lon_east + 180.0)/360.0)*width);
    let cols = [];
    if (c_west <= c_east) {
      for (let c = c_west; c <= c_east; c++) cols.push(c);
    } else {
      for (let c = c_west; c < width; c++) cols.push(c);
      for (let c = 0; c <= c_east; c++) cols.push(c);
    }
    let lw = cols.length;
    let cur = new Uint8Array(lw*lh);
    let nxt = new Uint8Array(lw*lh);

    for (let ly = 0; ly < lh; ly++) {
      let g_off = (r0 + ly)*width;
      let l_off = ly*lw;
      for (let lx = 0; lx < lw; lx++)
        cur[l_off + lx] = (f_ice_t[g_off + cols[lx]] >= 0.50) ? 1 : 0;
    }

    //Erode passes (never erode modern f_ice_2025 >= 0.50)
    for (let p = 0; p < passes; p++) {
      nxt.set(cur);
      for (let ly = 1; ly < lh - 1; ly++) {
        let g_off = (r0 + ly)*width;
        let l_off = ly*lw;
        for (let lx = 1; lx < lw - 1; lx++) {
          let li = l_off + lx;
          if (cur[li] === 0) continue;
          if (f_ice_2025[g_off + cols[lx]] >= 0.50) continue;
          if (cur[li - 1] === 0 || cur[li + 1] === 0 || cur[li - lw] === 0 || cur[li + lw] === 0)
            nxt[li] = 0;
        }
      }
      cur.set(nxt);
    }

    //Dilate passes
    for (let p = 0; p < passes; p++) {
      nxt.set(cur);
      for (let ly = 1; ly < lh - 1; ly++) {
        let l_off = ly*lw;
        for (let lx = 1; lx < lw - 1; lx++) {
          let li = l_off + lx;
          if (cur[li] === 1) continue;
          if (cur[li - 1] === 1 || cur[li + 1] === 1 || cur[li - lw] === 1 || cur[li + lw] === 1)
            nxt[li] = 1;
        }
      }
      cur.set(nxt);
    }

    for (let ly = 0; ly < lh; ly++) {
      let g_off = (r0 + ly)*width;
      let l_off = ly*lw;
      for (let lx = 0; lx < lw; lx++) {
        let gi = g_off + cols[lx];
        if (f_ice_2025[gi] >= 0.50) continue;
        f_ice_t[gi] = (cur[l_off + lx] === 1) ? 1.0 : f_ice_2025[gi];
      }
    }
  }

  /**
   * Computes clean 30" topographically anchored Younger Dryas / Early Holocene alpine & cordilleran
   * ice caps (Iceland, Patagonian Ice Sheet, Alaskan/Yukon/BC Cordilleran, NZ Southern Alps) for year < -7500.
   *
   * @param {Float32Array} arg0_f_ice_t
   * @param {Float32Array} arg1_f_land_osm_2025
   * @param {Int16Array} arg2_bed_i16
   * @param {Object} arg3_meta
   * @param {number} arg4_year
   */
  static applyTopographicPaleoIceCaps (arg0_f_ice_t, arg1_f_land_osm_2025, arg2_bed_i16, arg3_meta, arg4_year) {
    //Convert from parameters
    let f_ice_t = arg0_f_ice_t;
    let f_osm = arg1_f_land_osm_2025;
    let bed_i16 = arg2_bed_i16;
    let meta = arg3_meta;
    let year = arg4_year;

    //Guard clause
    if (year >= -7400) return;

    //Declare local instance variables
    let fx = meta.factor_x;
    let fy = meta.factor_y;
    let inv_sub = 1.0/(fx*fy);
    let scale_m = meta.scale_m_per_unit;
    let src_w = meta.src_w;
    let width = 4320;
    let height = 2160;

    //Helper to paint a clean 30" elevation-thresholded paleo ice cap inside a geographic box
    let paintBoxSnowline = (arg0_box, arg1_snowline_fn) => {
      let r0 = Math.max(0, Math.floor(((90.0 - arg0_box.lat_max)/180.0)*height));
      let r1 = Math.min(height - 1, Math.ceil(((90.0 - arg0_box.lat_min)/180.0)*height));
      let c0 = Math.max(0, Math.floor(((arg0_box.lon_min + 180.0)/360.0)*width));
      let c1 = Math.min(width - 1, Math.ceil(((arg0_box.lon_max + 180.0)/360.0)*width));

      for (let oy = r0; oy <= r1; oy++) {
        let lat = 90.0 - (oy + 0.5)*(180.0/height);
        let sy0 = oy*fy;
        let row_off = oy*width;

        for (let ox = c0; ox <= c1; ox++) {
          let idx = row_off + ox;
          if (f_osm[idx] <= 0.01) continue;
          let lon = -180.0 + (ox + 0.5)*(360.0/width);
          let snow_m = arg1_snowline_fn(lat, lon);
          if (snow_m > 9000.0) continue;

          let snow_hm = snow_m/scale_m;
          let sx0 = ox*fx;
          let hits = 0;

          for (let dy = 0; dy < fy; dy++) {
            let bed_off = (sy0 + dy)*src_w + sx0;
            for (let dx = 0; dx < fx; dx++)
              if (bed_i16[bed_off + dx] >= snow_hm) hits++;
          }

          let frac = hits*inv_sub;
          //Clean slight morphological sharpening so alpine domes stay coherent without pixel fuzz
          if (frac >= 0.45) {
            let v = Math.min(1.0, frac*1.25);
            if (v > f_ice_t[idx]) f_ice_t[idx] = v;
          }
        }
      }
    };

    //Retreat factor: 0.0 at -10000 (Younger Dryas peak) -> 1.0 at -7500 (modern alpine glaciers)
    let r_yd = Math.min(1.0, Math.max(0.0, (year + 10000.0)/2500.0));

    //1. Iceland Younger Dryas Ice Sheet (Norddahl & Petursson 2005; Patton et al. 2017): ~75-80% cover at -10000
    if (year < -7500) {
      let ice_snow = 120.0 + r_yd*1050.0;
      paintBoxSnowline(
        { lat_min: 63.2, lat_max: 66.7, lon_min: -24.8, lon_max: -13.2 },
        () => ice_snow
      );
    }

    //2. Patagonian Ice Sheet (Davies et al. 2020): continuous cordilleran spine 42..55.8S tapering organically east
    if (year < -7500) {
      paintBoxSnowline(
        { lat_min: -55.8, lat_max: -42.0, lon_min: -75.8, lon_max: -68.0 },
        (arg0_lat, arg1_lon) => {
          let spine_lon = -73.0 + Math.max(0.0, (-48.0 - arg0_lat))*0.10;
          let east_pen = Math.max(0.0, arg1_lon - spine_lon)*480.0;
          let base_m = (arg0_lat <= -48.0) ? 340.0 : (340.0 + (arg0_lat + 48.0)*95.0);
          return base_m + east_pen + r_yd*1300.0;
        }
      );
    }

    //3. Alaskan / Yukon / BC Cordilleran Ice Sheet (Menounos et al. 2017; Dalton et al. 2020) — organic eastern taper
    if (year < -7500) {
      paintBoxSnowline(
        { lat_min: 49.0, lat_max: 63.8, lon_min: -156.0, lon_max: -115.0 },
        (arg0_lat, arg1_lon) => {
          if (arg0_lat > 63.5) return 9999.0;
          let coast_spine_lon = -124.0 - Math.max(0.0, arg0_lat - 50.0)*2.2;
          let dist_from_spine = Math.abs(arg1_lon - coast_spine_lon);
          let base_m = 650.0 + dist_from_spine*185.0;
          return base_m + r_yd*1400.0;
        }
      );
    }

    //4. New Zealand Southern Alps (Barrell et al. 2011): Younger Dryas alpine glacier spine
    if (year < -7800) {
      let r_nz = Math.min(1.0, Math.max(0.0, (year + 10000.0)/2200.0));
      let nz_snow = 950.0 + r_nz*1250.0;
      paintBoxSnowline(
        { lat_min: -46.2, lat_max: -42.2, lon_min: 166.5, lon_max: 172.0 },
        () => nz_snow
      );
    }
  }

  /**
   * Dilates Fennoscandian Younger Dryas ice cleanly out to the Norwegian fjord calving margin
   * for year <= -8500 so no artificial green coastal ribbon remains along western/northern Norway.
   *
   * @param {Float32Array} arg0_f_ice_t
   * @param {Float32Array} arg1_f_land_osm_2025
   * @param {Float32Array} arg2_z_5m
   * @param {number} arg3_year
   */
  static dilateNorwegianCalvingMargin (arg0_f_ice_t, arg1_f_land_osm_2025, arg2_z_5m, arg3_year) {
    //Convert from parameters
    let f_ice_t = arg0_f_ice_t;
    let f_osm = arg1_f_land_osm_2025;
    let z_5m = arg2_z_5m;
    let year = arg3_year;

    //Guard clause
    if (year > -8500) return;

    //Declare local instance variables
    let width = 4320;
    let height = 2160;
    let r0 = Math.floor(((90.0 - 71.5)/180.0)*height);
    let r1 = Math.ceil(((90.0 - 58.0)/180.0)*height);
    let c0 = Math.floor(((4.3 + 180.0)/360.0)*width);
    let c1 = Math.ceil(((31.5 + 180.0)/360.0)*width);
    let passes = (year <= -9200) ? 5 : Math.max(1, Math.round(5.0*(-8500.0 - year)/700.0));
    let lw = c1 - c0 + 1;
    let lh = r1 - r0 + 1;
    let cur = new Float32Array(lw*lh);
    let nxt = new Float32Array(lw*lh);

    for (let ly = 0; ly < lh; ly++) {
      let g_off = (r0 + ly)*width + c0;
      let l_off = ly*lw;
      for (let lx = 0; lx < lw; lx++)
        cur[l_off + lx] = f_ice_t[g_off + lx];
    }

    for (let p = 0; p < passes; p++) {
      nxt.set(cur);
      for (let ly = 1; ly < lh - 1; ly++) {
        let g_off = (r0 + ly)*width + c0;
        let l_off = ly*lw;
        for (let lx = 1; lx < lw - 1; lx++) {
          let li = l_off + lx;
          if (cur[li] >= 0.99) continue;
          let g_idx = g_off + lx;
          if (f_osm[g_idx] <= 0.0 && z_5m[g_idx] < -20.0) continue;

          let m = cur[li];
          if (cur[li - 1] > m) m = cur[li - 1];
          if (cur[li + 1] > m) m = cur[li + 1];
          if (cur[li - lw] > m) m = cur[li - lw];
          if (cur[li + lw] > m) m = cur[li + lw];
          if (m > 0.6) nxt[li] = m;
        }
      }
      cur.set(nxt);
    }

    for (let ly = 0; ly < lh; ly++) {
      let g_off = (r0 + ly)*width + c0;
      let l_off = ly*lw;
      for (let lx = 0; lx < lw; lx++) {
        let g_idx = g_off + lx;
        if (cur[l_off + lx] > f_ice_t[g_idx] && (f_osm[g_idx] > 0.0 || z_5m[g_idx] >= -20.0))
          f_ice_t[g_idx] = cur[l_off + lx];
      }
    }
  }

  /**
   * Bilinearly samples a regular lat/lon grid (handles 0..360 or -180..180 longitude wrapping).
   *
   * @param {Float32Array} arg0_grid
   * @param {number} arg1_nlat
   * @param {number} arg2_nlon
   * @param {boolean} arg3_south_to_north
   * @param {boolean} arg4_lon_0_360
   * @param {number} arg5_lat_deg
   * @param {number} arg6_lon_deg
   *
   * @returns {number}
   */
  static sampleBilinear (arg0_grid, arg1_nlat, arg2_nlon, arg3_south_to_north, arg4_lon_0_360, arg5_lat_deg, arg6_lon_deg) {
    //Convert from parameters
    let grid = arg0_grid;
    let nlat = arg1_nlat;
    let nlon = arg2_nlon;
    let south_to_north = arg3_south_to_north;
    let lon_0_360 = arg4_lon_0_360;
    let lat_deg = arg5_lat_deg;
    let lon_deg = arg6_lon_deg;

    //Declare local instance variables
    let r_f = south_to_north
      ? ((lat_deg + 90.0)/180.0)*nlat - 0.5
      : ((90.0 - lat_deg)/180.0)*nlat - 0.5;
    let adj_lon = lon_0_360 ? ((lon_deg + 360.0) % 360.0) : (lon_deg + 180.0);
    let c_f = (adj_lon/360.0)*nlon - 0.5;

    if (r_f < 0.0) r_f = 0.0;
    if (r_f > nlat - 1.0001) r_f = nlat - 1.0001;
    if (c_f < 0.0) c_f += nlon;
    if (c_f >= nlon) c_f -= nlon;

    let r0 = r_f | 0;
    let c0 = c_f | 0;
    let r1 = (r0 + 1 < nlat) ? (r0 + 1) : r0;
    let c1 = (c0 + 1 < nlon) ? (c0 + 1) : 0;
    let dr = r_f - r0;
    let dc = c_f - c0;

    let v00 = grid[r0*nlon + c0];
    let v01 = grid[r0*nlon + c1];
    let v10 = grid[r1*nlon + c0];
    let v11 = grid[r1*nlon + c1];

    //Return statement
    return (v00*(1.0 - dc) + v01*dc)*(1.0 - dr) + (v10*(1.0 - dc) + v11*dc)*dr;
  }

  /**
   * Computes a 5-arcminute ocean-connected marine mask using iterative BFS flood-fill
   * from deep ocean cells (z < ocean_seed_depth_m), preventing endorheic basins (Dead Sea,
   * Qattara, Caspian, Death Valley, Turpan) from falsely flooding when z < 0.
   *
   * @param {Float32Array} arg0_z_5m - Mean/min elevation (m) at 4320x2160.
   * @param {number} arg1_width
   * @param {number} arg2_height
   * @param {number} arg3_seed_depth_m
   *
   * @returns {Uint8Array} 1 = ocean-connected marine basin, 0 = terrestrial or endorheic
   */
  static computeOceanConnectivity (arg0_z_5m, arg1_width, arg2_height, arg3_seed_depth_m) {
    //Convert from parameters
    let z = arg0_z_5m;
    let width = arg1_width;
    let height = arg2_height;
    let seed_depth = arg3_seed_depth_m;

    //Declare local instance variables
    let n = width*height;
    let ocean = new Uint8Array(n);
    let queue = new Int32Array(n);
    let q_head = 0;
    let q_tail = 0;

    //Seed deep open-ocean cells (exclude Caspian basin lon 46..55, lat 36..47 where depth reaches -1025m!)
    for (let y = 0; y < height; y++) {
      let lat_deg = 90.0 - (y + 0.5)*(180.0/height);
      let row_off = y*width;

      for (let x = 0; x < width; x++) {
        let idx = row_off + x;
        if (z[idx] < seed_depth) {
          let lon_deg = -180.0 + (x + 0.5)*(360.0/width);
          let in_caspian = (lat_deg >= 35.5 && lat_deg <= 48.0 && lon_deg >= 45.5 && lon_deg <= 55.5);
          let in_dead_sea = (lat_deg >= 30.5 && lat_deg <= 33.5 && lon_deg >= 34.8 && lon_deg <= 36.2);
          if (!in_caspian && !in_dead_sea) {
            ocean[idx] = 1;
            queue[q_tail++] = idx;
          }
        }
      }
    }

    //Iterative 4-neighbour BFS across cells with z < +2.0m (small coastal buffer so 30" sub-cells along coast stay eligible)
    for (; q_head < q_tail; ) {
      let curr = queue[q_head++];
      let cy = (curr/width) | 0;
      let cx = curr - cy*width;

      //West (wraps longitude)
      let wx = (cx === 0) ? (width - 1) : (cx - 1);
      let w_idx = cy*width + wx;
      if (ocean[w_idx] === 0 && z[w_idx] < 2.0) {
        ocean[w_idx] = 1;
        queue[q_tail++] = w_idx;
      }

      //East (wraps longitude)
      let ex = (cx === width - 1) ? 0 : (cx + 1);
      let e_idx = cy*width + ex;
      if (ocean[e_idx] === 0 && z[e_idx] < 2.0) {
        ocean[e_idx] = 1;
        queue[q_tail++] = e_idx;
      }

      //North
      if (cy > 0) {
        let n_idx = (cy - 1)*width + cx;
        if (ocean[n_idx] === 0 && z[n_idx] < 2.0) {
          ocean[n_idx] = 1;
          queue[q_tail++] = n_idx;
        }
      }

      //South
      if (cy < height - 1) {
        let s_idx = (cy + 1)*width + cx;
        if (ocean[s_idx] === 0 && z[s_idx] < 2.0) {
          ocean[s_idx] = 1;
          queue[q_tail++] = s_idx;
        }
      }
    }

    //Dilate marine eligibility by 1 cell (5') so coastal sub-pixels bordering ocean can respond to ΔRSL
    let marine_eligible = new Uint8Array(ocean);
    for (let y = 1; y < height - 1; y++) {
      let row_off = y*width;
      for (let x = 0; x < width; x++) {
        let idx = row_off + x;
        if (ocean[idx] === 1) continue;
        let wx = (x === 0) ? (width - 1) : (x - 1);
        let ex = (x === width - 1) ? 0 : (x + 1);
        if (
          ocean[row_off + wx] === 1 ||
          ocean[row_off + ex] === 1 ||
          ocean[idx - width] === 1 ||
          ocean[idx + width] === 1
        )
          marine_eligible[idx] = 1;
      }
    }

    //Return statement
    return marine_eligible;
  }

  /**
   * Builds a 5-arcminute minimum-bedrock elevation summary from the 30" bedrock DEM for rapid BFS ocean routing.
   *
   * @param {string} arg0_bed_halfm_path
   * @param {Object} arg1_etopo_meta
   * @param {number} arg2_out_w
   * @param {number} arg3_out_h
   *
   * @returns {Float32Array}
   */
  static buildMinBedrock5m (arg0_bed_halfm_path, arg1_etopo_meta, arg2_out_w, arg3_out_h) {
    //Convert from parameters
    let bed_path = arg0_bed_halfm_path;
    let meta = arg1_etopo_meta;
    let out_w = arg2_out_w;
    let out_h = arg3_out_h;

    //Declare local instance variables
    let factor_x = meta.factor_x;
    let factor_y = meta.factor_y;
    let fd = fs.openSync(bed_path, "r");
    let min_bed = new Float32Array(out_w*out_h);
    let row_bytes = meta.src_w*2;
    let band_buf = Buffer.allocUnsafe(factor_y*row_bytes);

    for (let oy = 0; oy < out_h; oy++) {
      fs.readSync(fd, band_buf, 0, band_buf.length, oy*factor_y*row_bytes);
      let slab = new Int16Array(band_buf.buffer, band_buf.byteOffset, factor_y*meta.src_w);
      let out_off = oy*out_w;

      for (let ox = 0; ox < out_w; ox++) {
        let min_v = 32767;
        let sx0 = ox*factor_x;

        for (let fy = 0; fy < factor_y; fy++) {
          let r_off = fy*meta.src_w + sx0;
          for (let fx = 0; fx < factor_x; fx++) {
            let v = slab[r_off + fx];
            if (v < min_v) min_v = v;
          }
        }
        min_bed[out_off + ox] = min_v*meta.scale_m_per_unit;
      }
    }

    fs.closeSync(fd);

    //Return statement
    return min_bed;
  }

  /**
   * Refines Persian Gulf & Hormuz paleo shoreline on 30" ETOPO bedrock using tapered Gaussian marine
   * bathymetric denoising, Chamfer morphological closing/opening west of 54.8E, and strict Gulf of Oman
   * BFS connectivity so jagged ship-track pits and twin-prong sill spikes never extrude at 10000-9000 BC.
   *
   * @param {Float32Array} arg0_delta_coast
   * @param {Float32Array} arg1_dz_5m
   * @param {Int16Array} arg2_bed_i16
   * @param {Float32Array} arg3_f_etopo_2025
   * @param {Object} arg4_meta
   */
  static cleanPersianGulfPaleoCoast (arg0_delta_coast, arg1_dz_5m, arg2_bed_i16, arg3_f_etopo_2025, arg4_meta) {
    //Convert from parameters
    let delta_coast = arg0_delta_coast;
    let dz_5m = arg1_dz_5m;
    let bed_i16 = arg2_bed_i16;
    let f_etopo_2025 = arg3_f_etopo_2025;
    let meta = arg4_meta;

    //Declare local instance variables
    let fx = meta.factor_x;
    let fy = meta.factor_y;
    let inv_sub = 1.0/(fx*fy);
    let scale_m = meta.scale_m_per_unit;
    let src_w = meta.src_w;
    let out_w = 4320;
    let out_h = 2160;
    let r0 = Math.floor(((90.0 - 29.6)/180.0)*out_h);
    let r1 = Math.ceil(((90.0 - 24.0)/180.0)*out_h);
    let c0 = Math.floor(((48.5 + 180.0)/360.0)*out_w);
    let c1 = Math.ceil(((57.5 + 180.0)/360.0)*out_w);
    let sy0 = r0*fy;
    let sx0 = c0*fx;
    let h = (r1 - r0)*fy;
    let w = (c1 - c0)*fx;
    let col_hormuz = Math.floor((54.6 + 180.0)*(src_w/360.0)) - sx0;
    let col_split = Math.floor((55.4 + 180.0)*(src_w/360.0)) - sx0;

    let z_raw = new Float32Array(w*h);
    for (let ly = 0; ly < h; ly++) {
      let sy = sy0 + ly;
      let ro = ly*w;
      for (let lx = 0; lx < w; lx++)
        z_raw[ro + lx] = bed_i16[sy*src_w + (sx0 + lx)]*scale_m;
    }

    //Separable Gaussian marine-only bathymetric smoothing (sigma = 16 px ~ 14 km), tapering smoothly to raw across 54.6..55.4E
    let sigma = 16.0;
    let rad = Math.ceil(sigma*2.5);
    let ker = new Float32Array(2*rad + 1);
    for (let d = -rad; d <= rad; d++)
      ker[d + rad] = Math.exp(-(d*d)/(2.0*sigma*sigma));

    let z_tmp = new Float32Array(w*h);
    let z_blur = new Float32Array(w*h);

    for (let y = 0; y < h; y++) {
      let ro = y*w;
      for (let x = 0; x < w; x++) {
        let z0 = z_raw[ro + x];
        if (x >= col_split || z0 > -2.0) {
          z_tmp[ro + x] = z0;
          continue;
        }
        let ws = 0.0;
        let vs = 0.0;
        for (let d = -rad; d <= rad; d++) {
          let xx = x + d;
          if (xx >= 0 && xx < w && z_raw[ro + xx] <= -1.0) {
            let wt = ker[d + rad];
            ws += wt;
            vs += wt*z_raw[ro + xx];
          }
        }
        z_tmp[ro + x] = (ws > 0.0) ? (vs/ws) : z0;
      }
    }

    for (let y = 0; y < h; y++) {
      let ro = y*w;
      for (let x = 0; x < w; x++) {
        let z0 = z_raw[ro + x];
        if (x >= col_split || z0 > -2.0) {
          z_blur[ro + x] = z0;
          continue;
        }
        let ws = 0.0;
        let vs = 0.0;
        for (let d = -rad; d <= rad; d++) {
          let yy = y + d;
          if (yy >= 0 && yy < h && z_raw[yy*w + x] <= -1.0) {
            let wt = ker[d + rad];
            ws += wt;
            vs += wt*z_tmp[yy*w + x];
          }
        }
        let zb = (ws > 0.0) ? (vs/ws) : z0;
        let t_fade = Math.min(1.0, Math.max(0.0, (x - col_hormuz)/(col_split - col_hormuz)));
        z_blur[ro + x] = zb*(1.0 - t_fade) + z0*t_fade;
      }
    }

    let wet = new Uint8Array(w*h);
    for (let ly = 0; ly < h; ly++) {
      let gy = r0 + Math.floor(ly/fy);
      let gx_row = gy*out_w;
      let dst_off = ly*w;
      for (let lx = 0; lx < w; lx++) {
        let gx = c0 + Math.floor(lx/fx);
        let dz = dz_5m[gx_row + gx];
        if ((z_blur[dst_off + lx] + dz) < 0.0) wet[dst_off + lx] = 1;
      }
    }

    let sr = Math.min(h - 1, Math.max(0, Math.floor(((90.0 - 25.0)/180.0)*meta.src_h) - sy0));
    let sc = Math.min(w - 1, Math.max(0, Math.floor(((57.2 + 180.0)/360.0)*src_w) - sx0));

    let bfsOman = (arg0_mask) => {
      let mask = arg0_mask;
      let vis = new Uint8Array(w*h);
      let q = new Int32Array(w*h);
      let qh = 0;
      let qt = 0;
      if (mask[sr*w + sc] === 1) {
        vis[sr*w + sc] = 1;
        q[qt++] = sr*w + sc;
      }
      for (; qh < qt; ) {
        let cur = q[qh++];
        let cy = (cur/w) | 0;
        let cx = cur - cy*w;
        if (cx > 0 && vis[cur - 1] === 0 && mask[cur - 1] === 1) {
          vis[cur - 1] = 1; q[qt++] = cur - 1;
        }
        if (cx < w - 1 && vis[cur + 1] === 0 && mask[cur + 1] === 1) {
          vis[cur + 1] = 1; q[qt++] = cur + 1;
        }
        if (cy > 0 && vis[cur - w] === 0 && mask[cur - w] === 1) {
          vis[cur - w] = 1; q[qt++] = cur - w;
        }
        if (cy < h - 1 && vis[cur + w] === 0 && mask[cur + w] === 1) {
          vis[cur + w] = 1; q[qt++] = cur + w;
        }
      }
      return vis;
    };

    let morphOpen = (arg0_mask, arg1_r) => {
      let mask = arg0_mask;
      let r = arg1_r;
      let inv = new Uint8Array(w*h);
      for (let i = 0; i < mask.length; i++) inv[i] = mask[i] ? 0 : 1;
      let d_in = BasinContourEngine.computeChamferDistance(inv, w, h);
      let core = new Uint8Array(w*h);
      for (let y = 0; y < h; y++) {
        let ro = y*w;
        for (let x = 0; x < w; x++) {
          let i = ro + x;
          core[i] = (x < col_hormuz) ? ((d_in[i] >= r) ? 1 : 0) : mask[i];
        }
      }
      let d_out = BasinContourEngine.computeChamferDistance(core, w, h);
      let out = new Uint8Array(w*h);
      for (let y = 0; y < h; y++) {
        let ro = y*w;
        for (let x = 0; x < w; x++) {
          let i = ro + x;
          out[i] = (x < col_hormuz) ? ((d_out[i] <= r && mask[i] === 1) ? 1 : 0) : mask[i];
        }
      }
      return out;
    };

    let morphClose = (arg0_mask, arg1_r) => {
      let mask = arg0_mask;
      let r = arg1_r;
      let d_out = BasinContourEngine.computeChamferDistance(mask, w, h);
      let dil = new Uint8Array(w*h);
      for (let y = 0; y < h; y++) {
        let ro = y*w;
        for (let x = 0; x < w; x++) {
          let i = ro + x;
          dil[i] = (z_raw[i] <= -1.0 && d_out[i] <= r) ? 1 : mask[i];
        }
      }
      let inv = new Uint8Array(w*h);
      for (let i = 0; i < dil.length; i++) inv[i] = dil[i] ? 0 : 1;
      let d_in = BasinContourEngine.computeChamferDistance(inv, w, h);
      let out = new Uint8Array(w*h);
      for (let y = 0; y < h; y++) {
        let ro = y*w;
        for (let x = 0; x < w; x++) {
          let i = ro + x;
          out[i] = (x < col_hormuz && d_in[i] >= r) ? 1 : mask[i];
        }
      }
      return out;
    };

    let vis = bfsOman(wet);
    vis = morphOpen(vis, 7.5);
    vis = bfsOman(vis);
    vis = morphClose(vis, 14.0);
    vis = morphOpen(vis, 7.5);
    vis = bfsOman(vis);

    //Write clean 5' delta_coast values back inside western/central Gulf (lon <= 57.1E, lat >= 24.4N)
    let r_end = Math.floor(((90.0 - 24.4)/180.0)*out_h);
    let c_end = Math.floor(((57.1 + 180.0)/360.0)*out_w);
    for (let gy = r0; gy < Math.min(r1, r_end); gy++) {
      let ly0 = (gy - r0)*fy;
      for (let gx = c0; gx < Math.min(c1, c_end); gx++) {
        let lx0 = (gx - c0)*fx;
        let land_hits = 0;

        for (let dy = 0; dy < fy; dy++) {
          let row_off = (ly0 + dy)*w + lx0;
          for (let dx = 0; dx < fx; dx++)
            if (vis[row_off + dx] === 0) land_hits++;
        }

        let idx = gy*out_w + gx;
        let f_paleo = land_hits*inv_sub;
        delta_coast[idx] = f_paleo - f_etopo_2025[idx];
      }
    }
  }

  /**
   * Executes Step C: extracts GLAC-1D 100-yr paleo slices, computes GIA ΔRSL shoreline anomalies
   * and paleo glaciation (incl. permanent ice shelves) for all unique paleo timesteps up to gia_cutoff_year.
   *
   * @returns {Promise<void>}
   */
  static async run () {
    //Declare local instance variables
    let cfg = GeoGridUtils.loadConfig("model_config.json5");
    let hcfg = GeoGridUtils.loadConfig("hydro_geomorph_basins.json5");
    let paths = GeoGridUtils.resolvePaths();
    let paleo_dir = path.join(paths.intermediate_dir, "paleo_slices");
    let years = GeoGridUtils.getTimestepYears(cfg);
    let paleo_years = years.filter((y) => y < cfg.paleo.gia_cutoff_year);

    if (!fs.existsSync(paleo_dir))
      fs.mkdirSync(paleo_dir, { recursive: true });

    //Check if all paleo slices are already cached (unless --force is passed)
    let force = process.argv.includes("--force");
    let missing_years = force
      ? paleo_years.slice()
      : paleo_years.filter((y) => {
        let c_path = path.join(paleo_dir, `delta_coast_${y}.f32`);
        let i_path = path.join(paleo_dir, `f_ice_${y}.f32`);
        return !(fs.existsSync(c_path) && fs.existsSync(i_path));
      });

    if (missing_years.length === 0) {
      console.log(`[C_paleo] All ${paleo_years.length} paleo GIA & ice slices already cached.`);
      return;
    }

    //1. Open PMIP4 GLAC-1D dataset (26-0 ka BP at 100-yr intervals)
    let glac_path = path.join(paths.raw_dir, "TOPicemsk.GLACD26kN9894GE90227A6005GGrBgic.nc");
    console.log(`[C_paleo] Opening GLAC-1D: ${glac_path} ...`);
    let glac_nc = await NetCDF4Writer.openRead(glac_path);
    let keys = glac_nc.keys();
    console.log(`[C_paleo] GLAC-1D variables: ${keys.join(", ")}`);

    //Detect coordinate and variable names flexibly
    let time_key = keys.find((k) => /^t(ime)?$/i.test(k)) || "T122KP11";
    let lat_key = keys.find((k) => /^y|^lat/i.test(k));
    let lon_key = keys.find((k) => /^x|^lon/i.test(k));
    let hdc_key = keys.find((k) => k === "HDC") || keys.find((k) => /^hdc$/i.test(k));
    let hdcb_key = keys.find((k) => k === "HDCB") || keys.find((k) => /^hdcb$/i.test(k));
    let icem_key = keys.find((k) => k === "ICEM") || keys.find((k) => /^icem/i.test(k));

    let t_vals = Array.from(glac_nc.getValue(time_key));
    let lat_vals = Array.from(glac_nc.getValue(lat_key));
    let lon_vals = Array.from(glac_nc.getValue(lon_key));
    let shape = glac_nc.getShape(hdc_key);
    let ntime = shape[0];
    let nlat = shape[1];
    let nlon = shape[2];
    let south_to_north = lat_vals[0] < lat_vals[lat_vals.length - 1];
    let lon_0_360 = lon_vals[lon_vals.length - 1] > 180.5;

    //Normalise time axis to calendar years (GLAC-1D uses ka BP e.g. -26.0 .. 0.0 or 26.0 .. 0.0)
    let max_abs_t = Math.max(...t_vals.map((v) => Math.abs(v)));
    let cal_years = t_vals.map((v) => {
      if (max_abs_t <= 40.0) {
        let ka_bp = (v <= 0.0) ? -v : v;
        return cfg.paleo.bp_reference_year - ka_bp*1000.0;
      }
      return cfg.paleo.bp_reference_year - Math.abs(v);
    });

    //Locate present-day (0 ka BP = 1950 AD) index in GLAC-1D
    let present_idx = 0;
    let best_dist = Infinity;
    for (let i = 0; i < ntime; i++) {
      let d = Math.abs(cal_years[i] - cfg.paleo.bp_reference_year);
      if (d < best_dist) {
        best_dist = d;
        present_idx = i;
      }
    }

    console.log(`[C_paleo] GLAC-1D grid ${ntime}x${nlat}x${nlon}, present_idx=${present_idx} (cal_yr=${cal_years[present_idx]}), s2n=${south_to_north}, lon360=${lon_0_360}`);

    let hdc_present = Float32Array.from(glac_nc.getSlice(hdc_key, [[present_idx, present_idx + 1], [0, nlat], [0, nlon]]));
    let hdcb_present = hdcb_key
      ? Float32Array.from(glac_nc.getSlice(hdcb_key, [[present_idx, present_idx + 1], [0, nlat], [0, nlon]]))
      : hdc_present;
    let icem_present = Float32Array.from(glac_nc.getSlice(icem_key, [[present_idx, present_idx + 1], [0, nlat], [0, nlon]]));

    //2. Load modern ETOPO metadata, baseline land fraction, ice anchor, and build 5' min-bedrock for fast marine BFS
    let bed_halfm_path = path.join(paths.intermediate_dir, "etopo_bed_halfm.i16");
    let etopo_meta = JSON.parse(fs.readFileSync(path.join(paths.intermediate_dir, "etopo_meta.json"), "utf8"));
    let f_etopo_land_2025 = GeoGridUtils.readFloat32Bin(path.join(paths.intermediate_dir, "f_etopo_land_2025.f32"));
    let f_ice_2025 = GeoGridUtils.readFloat32Bin(path.join(paths.intermediate_dir, "f_ice_2025.f32"));
    let f_land_osm_2025 = GeoGridUtils.readFloat32Bin(path.join(paths.intermediate_dir, "f_land_osm_2025.f32"));

    let min_bed_5m_path = path.join(paths.intermediate_dir, "etopo_min_bed_5m.f32");
    let min_bed_5m;
    if (fs.existsSync(min_bed_5m_path)) {
      min_bed_5m = GeoGridUtils.readFloat32Bin(min_bed_5m_path);
    } else {
      console.log("[C_paleo] Building 5' minimum bedrock summary from 30\" DEM ...");
      min_bed_5m = this.buildMinBedrock5m(bed_halfm_path, etopo_meta, cfg.grid.width, cfg.grid.height);
      GeoGridUtils.writeFloat32Bin(min_bed_5m_path, min_bed_5m);
    }

    //Pre-load entire 30" Int16 bedrock DEM (~1.86 GB) once into RAM for ultra-fast 100-slice evaluation!
    console.log("[C_paleo] Loading 30\" bedrock DEM into RAM (~1.86 GB) ...");
    let bed_buf = fs.readFileSync(bed_halfm_path);
    let bed_i16 = new Int16Array(bed_buf.buffer, bed_buf.byteOffset, bed_buf.byteLength/2);
    let factor_x = etopo_meta.factor_x;
    let factor_y = etopo_meta.factor_y;
    let inv_sub = 1.0/(factor_x*factor_y);
    let out_h = cfg.grid.height;
    let out_w = cfg.grid.width;
    let scale_m = etopo_meta.scale_m_per_unit;
    let shelf_thresh = cfg.paleo.glac1d_shelf_thickness_m;
    let src_w = etopo_meta.src_w;
    let total_px = out_w*out_h;

    //Holocene estuarine/deltaic sediment back-stripping (built-up Ob, Taz, Yenisei, Khatanga estuaries; Lena, Yana, Indigirka, Kolyma deltas)
    console.log("[C_paleo] Building Holocene estuarine/deltaic infill model ...");
    let infill = new HoloceneInfillModel(hcfg.holocene_infill, bed_i16, etopo_meta);

    //Helper to interpolate a GLAC-1D field to any target calendar year
    let getGlac1DSlice = (arg0_year) => {
      let yr = arg0_year;
      let i0 = 0;
      let i1 = 0;

      //Find bracketing indices (cal_years may be ascending or descending)
      let asc = cal_years[cal_years.length - 1] > cal_years[0];
      if (asc) {
        if (yr <= cal_years[0]) {
          i0 = 0; i1 = 0;
        } else if (yr >= cal_years[ntime - 1]) {
          i0 = ntime - 1; i1 = ntime - 1;
        } else {
          for (let i = 0; i < ntime - 1; i++) {
            if (yr >= cal_years[i] && yr <= cal_years[i + 1]) {
              i0 = i; i1 = i + 1;
              break;
            }
          }
        }
      } else {
        if (yr >= cal_years[0]) {
          i0 = 0; i1 = 0;
        } else if (yr <= cal_years[ntime - 1]) {
          i0 = ntime - 1; i1 = ntime - 1;
        } else {
          for (let i = 0; i < ntime - 1; i++) {
            if (yr <= cal_years[i] && yr >= cal_years[i + 1]) {
              i0 = i; i1 = i + 1;
              break;
            }
          }
        }
      }

      let w1 = (i0 === i1) ? 0.0 : (yr - cal_years[i0])/(cal_years[i1] - cal_years[i0]);
      let w0 = 1.0 - w1;
      let n_g = nlat*nlon;

      let hdc0 = glac_nc.getSlice(hdc_key, [[i0, i0 + 1], [0, nlat], [0, nlon]]);
      let hdcb0 = hdcb_key ? glac_nc.getSlice(hdcb_key, [[i0, i0 + 1], [0, nlat], [0, nlon]]) : hdc0;
      let icem0 = glac_nc.getSlice(icem_key, [[i0, i0 + 1], [0, nlat], [0, nlon]]);

      let dz_grid = new Float32Array(n_g);
      let ice_delta_grid = new Float32Array(n_g);
      let ice_contaminated = new Uint8Array(n_g);

      let fillSliceArrays = (arg0_hdc_t, arg1_hdcb_t, arg2_icem_t) => {
        for (let k = 0; k < n_g; k++) {
          let ic_t = arg2_icem_t[k];
          let ic_p = icem_present[k];
          let shelf_t = cfg.encoding.include_permanent_ice_shelves && ((arg0_hdc_t[k] - arg1_hdcb_t[k]) > shelf_thresh);
          let shelf_p = cfg.encoding.include_permanent_ice_shelves && ((hdc_present[k] - hdcb_present[k]) > shelf_thresh);
          let eff_ice_t = shelf_t ? Math.max(ic_t, 1.0) : ic_t;
          let eff_ice_p = shelf_p ? Math.max(ic_p, 1.0) : ic_p;
          ice_delta_grid[k] = eff_ice_t - eff_ice_p;

          //Only ice-free cells give pure bedrock GIA + eustatic ΔRSL without ice-sheet thickness contamination
          if (ic_t > 0.02 || ic_p > 0.02) {
            ice_contaminated[k] = 1;
            dz_grid[k] = 0.0;
          } else {
            let dz = arg0_hdc_t[k] - hdc_present[k];
            if (dz < -350.0) dz = -350.0;
            if (dz > 160.0) dz = 160.0;
            dz_grid[k] = dz;
          }
        }

        //Smoothly inpaint dz_grid across glaciated cells from neighbouring ice-free GIA cells (prevents +1500m ice-margin artefacts)
        let tmp = new Float32Array(dz_grid);
        for (let iter = 0; iter < 40; iter++) {
          for (let r = 0; r < nlat; r++) {
            let row_off = r*nlon;
            let r_up = (r > 0) ? (r - 1)*nlon : row_off;
            let r_dn = (r < nlat - 1) ? (r + 1)*nlon : row_off;

            for (let c = 0; c < nlon; c++) {
              let k = row_off + c;
              if (ice_contaminated[k] === 0) continue;
              let c_l = (c === 0) ? (nlon - 1) : (c - 1);
              let c_r = (c === nlon - 1) ? 0 : (c + 1);
              tmp[k] = 0.25*(dz_grid[row_off + c_l] + dz_grid[row_off + c_r] + dz_grid[r_up + c] + dz_grid[r_dn + c]);
            }
          }
          dz_grid.set(tmp);
        }
      };

      if (i0 === i1) {
        fillSliceArrays(hdc0, hdcb0, icem0);
      } else {
        let hdc1 = glac_nc.getSlice(hdc_key, [[i1, i1 + 1], [0, nlat], [0, nlon]]);
        let hdcb1 = hdcb_key ? glac_nc.getSlice(hdcb_key, [[i1, i1 + 1], [0, nlat], [0, nlon]]) : hdc1;
        let icem1 = glac_nc.getSlice(icem_key, [[i1, i1 + 1], [0, nlat], [0, nlon]]);
        let hdc_blend = new Float32Array(n_g);
        let hdcb_blend = new Float32Array(n_g);
        let icem_blend = new Float32Array(n_g);

        for (let k = 0; k < n_g; k++) {
          hdc_blend[k] = w0*hdc0[k] + w1*hdc1[k];
          hdcb_blend[k] = w0*hdcb0[k] + w1*hdcb1[k];
          icem_blend[k] = w0*icem0[k] + w1*icem1[k];
        }
        fillSliceArrays(hdc_blend, hdcb_blend, icem_blend);
      }

      return { dz_grid: dz_grid, ice_delta_grid: ice_delta_grid };
    };

    //Precompute modern ocean connectivity ONCE: any cell connected to modern ocean is eligible
    //to emerge as solid dry land when sea level drops (Doggerland, Sundaland, Sahul, Beringia, Persian Gulf, etc.),
    //while disconnected endorheic depressions (Caspian, Dead Sea, Qattara, Death Valley, Turpan) remain protected!
    console.log("[C_paleo] Computing modern ocean connectivity mask ...");
    let marine_modern = this.computeOceanConnectivity(min_bed_5m, out_w, out_h, cfg.paleo.ocean_seed_depth_m);

    //3. Compute paleo Δcoast and f_ice for every missing paleo year
    let z_5m_work = new Float32Array(total_px);
    let dz_5m_work = new Float32Array(total_px);
    let t_start = Date.now();

    for (let yi = 0; yi < missing_years.length; yi++) {
      let year = missing_years[yi];
      let t_slice = Date.now();
      let { dz_grid, ice_delta_grid } = getGlac1DSlice(year);

      let delta_coast = new Float32Array(total_px);
      let f_ice_t = new Float32Array(total_px);

      //Step 3a: Upsample smooth GIA Δz (with calibrated eustatic stand for 10000-6000 BC) and sharpened paleo ice anomaly
      let eustatic_cal = this.getCalibratedEustaticOffsetM(year);
      let infill_missing = 1.0 - infill.getCompletion(year);
      let laurentide_mainland_mult = 1.0;
      if (year >= -4500) {
        laurentide_mainland_mult = 0.0;
      } else if (year > -4800) {
        laurentide_mainland_mult = (-4500.0 - year)/300.0;
      }

      for (let y = 0; y < out_h; y++) {
        let lat_deg = 90.0 - (y + 0.5)*(180.0/out_h);
        let row_off = y*out_w;

        for (let x = 0; x < out_w; x++) {
          let idx = row_off + x;
          let lon_deg = -180.0 + (x + 0.5)*(360.0/out_w);
          let dz = this.sampleBilinear(dz_grid, nlat, nlon, south_to_north, lon_0_360, lat_deg, lon_deg) + eustatic_cal;
          let dice_raw = this.sampleBilinear(ice_delta_grid, nlat, nlon, south_to_north, lon_0_360, lat_deg, lon_deg);

          dz_5m_work[idx] = dz;
          z_5m_work[idx] = min_bed_5m[idx] + dz - infill.max_5m[idx]*infill_missing;

          //R2-D03a: Crisp 0.50 contour threshold for coarse GLAC-1D paleo ice expansion/retreat
          //eliminates blocky translucent pale-blue-white halos (e.g. Keewatin at 6000 BC)
          let dice = 0.0;
          if (dice_raw >= 0.50) {
            dice = 1.0;
          } else if (dice_raw <= -0.50) {
            dice = -1.0;
          }

          //D03: Never carve coarse 1deg vertical GLAC-1D cutouts out of high-res ETOPO Greenland proper!
          if (dice < 0.0 && this.isGreenlandProper(lat_deg, lon_deg))
            dice = 0.0;

          //R2-D02: Lock Antarctica to modern ice-sheet/shelf geometry for year >= -6500 (RAISED Consortium / Bentley 2014)
          //Eliminates blocky GLAC-1D -46W Weddell protrusion & blue hairline slit at -6000!
          if (lat_deg <= -60.0 && year >= -6500)
            dice = 0.0;

          //D11: Mainland Laurentide collapsed south of Hudson Strait by 4500 BC (eliminates Ungava blob at 4000 BC)
          if (lat_deg < 61.5 && lat_deg > 42.0 && lon_deg >= -98.0 && lon_deg <= -55.0 && dice > 0.0)
            dice *= laurentide_mainland_mult;

          //D07 & R3-D02: Clip open-ocean North Atlantic / Greenland Sea / Barents Sea paleo ice (>55N, lon > -58W)
          //to grounded continental shelf (never cut Hudson Strait at -75W or Laurentide inland lakes/channels!)
          if (lat_deg > 55.0 && dice > 0.0 && f_ice_2025[idx] < 0.01) {
            let in_svalbard = (lat_deg >= 75.0 && lat_deg <= 82.0 && lon_deg >= 5.0 && lon_deg <= 35.0);
            let in_fjl = (lat_deg >= 78.0 && lat_deg <= 83.0 && lon_deg >= 38.0 && lon_deg <= 70.0);
            if (in_svalbard || in_fjl) {
              if (f_land_osm_2025[idx] <= 0.0 && z_5m_work[idx] < -70.0)
                dice = 0.0;
            } else if (lon_deg > -58.0 && marine_modern[idx] === 1 && z_5m_work[idx] < -200.0) {
              dice = 0.0;
            }
          }

          let ice_val = f_ice_2025[idx] + dice;
          if (ice_val < 0.0) ice_val = 0.0;
          if (ice_val > 1.0) ice_val = 1.0;
          if (ice_val < 0.05 && f_ice_2025[idx] === 0.0) ice_val = 0.0;
          f_ice_t[idx] = ice_val;
        }
      }

      //Step 3a2: Post-process regional paleo glaciology (zero floating blobs, zero trapped lagoons, organic margins)
      //R2-D02: For Antarctica (< -6500), seal all marine lagoons trapped behind advancing ice shelves via open-ocean BFS
      //        from lat = -58S, clip outer detached shelf blobs, and smooth paleo advance edges
      if (year < -6500) {
        this.fillEnclosedIceHoles(f_ice_t, f_land_osm_2025, {
          lat_max: -58.0,
          lat_min: -90.0,
          lon_max: 180.0,
          lon_min: -180.0,
          marine_only: true,
          seed_all_edges: false,
          wrap_lon: true,
        });

        this.clipIceToConnectedShelf(f_ice_t, f_ice_2025, f_land_osm_2025, z_5m_work, {
          lat_max: -58.0,
          lat_min: -90.0,
          lon_max: 180.0,
          lon_min: -180.0,
          max_shelf_depth_m: -450.0,
          min_ice_thresh: 0.50,
          seed_from_land: false,
          seed_ice_thresh: 0.80,
          wrap_lon: true,
        });

        //Second open-ocean flood-fill after shelf clipping guarantees zero enclosed lagoons remain behind Ronne/Ross!
        this.fillEnclosedIceHoles(f_ice_t, f_land_osm_2025, {
          lat_max: -58.0,
          lat_min: -90.0,
          lon_max: 180.0,
          lon_min: -180.0,
          marine_only: true,
          seed_all_edges: false,
          wrap_lon: true,
        });

        //R3-D01: Weddell/Ronne-Filchner smooth calving front + southward fill (lat -71..-80.5S, lon -62..-25W):
        //smooths across the narrow 31-35W Filchner Trough notch (eliminating any vertical step at 34.5W),
        //seals the trapped eastern Weddell pocket, AND erases the 1-px modern Ronne calving-front hairline!
        let wed_r0 = Math.floor(((90.0 - (-71.0))/180.0)*out_h);
        let wed_r1 = Math.ceil(((90.0 - (-80.5))/180.0)*out_h);
        let wed_c0 = Math.floor(((-62.0 + 180.0)/360.0)*out_w);
        let wed_c1 = Math.ceil(((-25.0 + 180.0)/360.0)*out_w);
        let wed_w = wed_c1 - wed_c0 + 1;
        let top_y_arr = new Int32Array(wed_w);

        for (let lx = 0; lx < wed_w; lx++) {
          let cx = wed_c0 + lx;
          let found = wed_r1;
          for (let cy = wed_r0; cy <= wed_r1; cy++) {
            if (f_ice_t[cy*out_w + cx] >= 0.50) {
              found = cy;
              break;
            }
          }
          top_y_arr[lx] = found;
        }

        //Close narrow southward notches (<= 55 columns ~4.5 deg lon) along the Weddell front via upper-hull envelope
        let smooth_top = new Int32Array(top_y_arr);
        for (let lx = 0; lx < wed_w; lx++) {
          let left_min = top_y_arr[lx];
          let right_min = top_y_arr[lx];
          let max_span = 55;
          for (let d = 1; d <= max_span; d++) {
            if (lx - d >= 0 && top_y_arr[lx - d] < left_min) left_min = top_y_arr[lx - d];
            if (lx + d < wed_w && top_y_arr[lx + d] < right_min) right_min = top_y_arr[lx + d];
          }
          let env = Math.max(left_min, right_min);
          if (env < smooth_top[lx]) smooth_top[lx] = env;
        }

        for (let lx = 0; lx < wed_w; lx++) {
          let cx = wed_c0 + lx;
          let top_y = smooth_top[lx];
          if (top_y < wed_r1) {
            for (let cy = top_y; cy <= wed_r1; cy++) {
              let idx = cy*out_w + cx;
              if (f_etopo_land_2025[idx] < 0.35)
                f_ice_t[idx] = 1.0;
            }
          }
        }

        //R3-D01: Ross Sea morphological opening (6-pass erode/dilate across antimeridian 165E..-150W)
        //completely prunes the C-shaped hook/claw at 170-177W while keeping modern Ross Ice Shelf intact!
        this.openPaleoIceAntimeridian(f_ice_t, f_ice_2025, {
          lat_max: -72.5,
          lat_min: -79.5,
          lon_east: -150.0,
          lon_west: 165.0,
          passes: 6,
        });

        //Seal any 1-pixel partial modern-calving-front hairline trapped inside solid paleo Antarctic ice (< -65S)
        let ant_r_start = Math.floor(((90.0 - (-65.0))/180.0)*out_h);
        for (let pass = 0; pass < 2; pass++) {
          for (let y = ant_r_start; y < out_h - 1; y++) {
            let row_off = y*out_w;
            for (let x = 0; x < out_w; x++) {
              let idx = row_off + x;
              if (f_ice_t[idx] >= 0.99 || f_etopo_land_2025[idx] >= 0.35) continue;
              let xl = (x === 0) ? (out_w - 1) : (x - 1);
              let xr = (x === out_w - 1) ? 0 : (x + 1);
              let hi = 0;
              if (f_ice_t[row_off + xl] >= 0.90) hi++;
              if (f_ice_t[row_off + xr] >= 0.90) hi++;
              if (f_ice_t[idx - out_w] >= 0.90) hi++;
              if (f_ice_t[idx + out_w] >= 0.90) hi++;
              if (hi >= 3) f_ice_t[idx] = 1.0;
            }
          }
        }
      }

      //D07 & R3-D02: Clip Greenland & East Greenland Sea paleo advance strictly east of Davis Strait (lon >= -58W)
      //so Hudson Strait (-75..-64W) is never vertically sliced!
      this.clipIceToConnectedShelf(f_ice_t, f_ice_2025, f_land_osm_2025, z_5m_work, {
        lat_max: 85.0,
        lat_min: 58.0,
        lon_max: -5.0,
        lon_min: -58.0,
        max_shelf_depth_m: -200.0,
        min_ice_thresh: 0.05,
        seed_from_land: true,
        seed_ice_thresh: 0.60,
        wrap_lon: false,
      });

      //D16: Extend Fennoscandian Younger Dryas ice cleanly to the Norwegian fjord calving margin (year <= -8500)
      this.dilateNorwegianCalvingMargin(f_ice_t, f_land_osm_2025, z_5m_work, year);

      //D08 & D19: Topographically anchored Younger Dryas Iceland, Patagonia, Alaska/BC Cordilleran & NZ Southern Alps
      this.applyTopographicPaleoIceCaps(f_ice_t, f_land_osm_2025, bed_i16, etopo_meta, year);

      //R2-D01: Organic connected-component filter drops tiny isolated detached eastern BC GLAC-1D dots with zero box cuts!
      this.removeIsolatedPaleoIceBlobs(f_ice_t, f_ice_2025, {
        lat_max: 60.0,
        lat_min: 47.0,
        lon_max: -110.0,
        lon_min: -132.0,
        min_component_px: 55,
      });

      //R2-D03a & R2-D03b: Smooth Laurentide/Keewatin margins & seal enclosed interior lakes/slits/channels under ice!
      if (year <= -4500) {
        this.smoothPaleoIceMargins(f_ice_t, f_ice_2025, {
          lat_max: 82.0,
          lat_min: 43.0,
          lon_max: -52.0,
          lon_min: -140.0,
        });

        this.fillEnclosedIceHoles(f_ice_t, f_land_osm_2025, {
          lat_max: 82.0,
          lat_min: 43.0,
          lon_max: -52.0,
          lon_min: -140.0,
          marine_only: false,
          seed_all_edges: true,
          wrap_lon: false,
        });
      }

      //Step 3b: Paleo ocean connectivity flood-fill (captures marine transgressions like Yoldia/Littorina/Champlain/Tyrrell Seas)
      let marine_paleo = this.computeOceanConnectivity(z_5m_work, out_w, out_h, cfg.paleo.ocean_seed_depth_m);

      //Step 3c: Evaluate 30" (10x10) sub-pixel shoreline shift on cells connected to EITHER modern OR paleo ocean!
      let thresh = cfg.paleo.gia_flood_threshold_m;
      for (let y = 0; y < out_h; y++) {
        let lat_deg = 90.0 - (y + 0.5)*(180.0/out_h);
        let row_off = y*out_w;
        let sy0 = y*factor_y;

        for (let x = 0; x < out_w; x++) {
          let idx = row_off + x;
          let dz = dz_5m_work[idx];

          let infill_m = infill.max_5m[idx]*infill_missing;
          if ((marine_modern[idx] === 0 && marine_paleo[idx] === 0) || (Math.abs(dz) < thresh && infill_m <= 0.0)) continue;

          //Quick bound skip: if min bedrock + dz >= 0 and already full land, or max possible < -600m and already ocean, skip sub-loop
          let min_z = z_5m_work[idx];
          if (min_z >= 0.0 && f_etopo_land_2025[idx] >= 1.0) continue;
          if (min_z < -600.0 && f_etopo_land_2025[idx] <= 0.0) continue;

          let sx0 = x*factor_x;
          let land_hits = 0;
          let dz_halfm = dz/scale_m;

          for (let fy = 0; fy < factor_y; fy++) {
            let bed_off = (sy0 + fy)*src_w + sx0;
            for (let fx = 0; fx < factor_x; fx++) {
              let off_halfm = (infill_m > 0.0) ? infill.getMissingThicknessM(sy0 + fy, sx0 + fx, infill_missing)/scale_m : 0.0;
              if ((bed_i16[bed_off + fx] + dz_halfm - off_halfm) >= 0.0)
                land_hits++;
            }
          }

          let f_paleo_etopo = land_hits*inv_sub;
          let dc = f_paleo_etopo - f_etopo_land_2025[idx];

          //D17: Preserve Lake Maracaibo as a permanent enclosed freshwater lake during sea-level lowstands!
          let lon_deg = -180.0 + (x + 0.5)*(360.0/out_w);
          if (dc > 0.0 && lat_deg >= 9.0 && lat_deg <= 10.65 && lon_deg >= -72.2 && lon_deg <= -71.0)
            dc = 0.0;

          //D20: Great Bahama Bank was submerged by 5000 BC (eliminates salt-and-pepper carbonate-bank noise at 4000 BC)
          if (year >= -5000 && dc > 0.0 && lat_deg >= 21.5 && lat_deg <= 27.2 && lon_deg >= -79.6 && lon_deg <= -74.5)
            dc = 0.0;

          delta_coast[idx] = dc;
        }
      }

      //D20: Smooth emergent Great & Little Bahama Banks for year < -5000 (eliminates salt-and-pepper stippling at 6000 BC)
      if (year < -5000) {
        let br0 = Math.floor(((90.0 - 27.2)/180.0)*out_h);
        let br1 = Math.ceil(((90.0 - 21.5)/180.0)*out_h);
        let bc0 = Math.floor(((-79.6 + 180.0)/360.0)*out_w);
        let bc1 = Math.ceil(((-74.5 + 180.0)/360.0)*out_w);
        let bw = bc1 - bc0 + 1;
        let bh = br1 - br0 + 1;
        let bbuf = new Float32Array(bw*bh);

        for (let pass = 0; pass < 2; pass++) {
          for (let ly = 1; ly < bh - 1; ly++) {
            let gy = br0 + ly;
            for (let lx = 1; lx < bw - 1; lx++) {
              let gx = bc0 + lx;
              let sum = 0.0;
              let cnt_hi = 0;
              for (let dy = -1; dy <= 1; dy++) {
                let row_i = (gy + dy)*out_w + gx;
                for (let dx = -1; dx <= 1; dx++) {
                  let v = delta_coast[row_i + dx];
                  sum += v;
                  if (v > 0.45) cnt_hi++;
                }
              }
              let avg = sum/9.0;
              if (cnt_hi >= 5) bbuf[ly*bw + lx] = Math.max(avg, 0.95);
              else if (cnt_hi <= 3) bbuf[ly*bw + lx] = 0.0;
              else bbuf[ly*bw + lx] = avg;
            }
          }
          for (let ly = 1; ly < bh - 1; ly++) {
            let g_off = (br0 + ly)*out_w + bc0;
            let l_off = ly*bw;
            for (let lx = 1; lx < bw - 1; lx++)
              delta_coast[g_off + lx] = bbuf[l_off + lx];
          }
        }
      }

      //Persian Gulf & Hormuz clean 30" ocean connectivity + ship-track smoothing (year < -5000):
      //Eliminates disconnected -60m sinkholes (e.g. 27.38N 50.08E) and GEBCO track spikes at 10000 BC!
      if (year < -5000)
        this.cleanPersianGulfPaleoCoast(delta_coast, dz_5m_work, bed_i16, f_etopo_land_2025, etopo_meta);

      GeoGridUtils.writeFloat32Bin(path.join(paleo_dir, `delta_coast_${year}.f32`), delta_coast);
      GeoGridUtils.writeFloat32Bin(path.join(paleo_dir, `f_ice_${year}.f32`), f_ice_t);

      if ((yi + 1) % 10 === 0 || yi === 0 || yi === missing_years.length - 1)
        console.log(`[C_paleo] Processed ${yi + 1}/${missing_years.length} (year ${year}) in ${((Date.now() - t_slice)/1000).toFixed(2)}s`);
    }

    glac_nc.close();
    console.log(`[C_paleo] Completed all ${missing_years.length} paleo GIA & ice slices in ${((Date.now() - t_start)/1000).toFixed(1)}s.`);
  }
}

if (require.main === module) {
  C_BuildPaleoGiaIce.run().catch((err) => {
    console.error("[C_paleo] Fatal error:", err);
    process.exit(1);
  });
}

module.exports = C_BuildPaleoGiaIce;
