"use strict";

let fs = require("node:fs");
let path = require("node:path");

let BasinContourEngine = require("../../core/BasinContourEngine");
let GeoGridUtils = require("../../core/GeoGridUtils");
let NetCDF4Writer = require("../../core/NetCDF4Writer");

/**
 * Step E: Composites all 218 timesteps (-10000 to 2025) under strict signed precedence
 * (Land > 0, Glaciated + Permanent Floating Ice Shelves < 0, Sea/Lakes/Reservoirs = 0),
 * integrating:
 *  - Inclusive 2025 OSM coastlines + fixed GIA marine_modern || marine_paleo shelves
 *  - Official RCE Vos et al. (2021) Holocene Netherlands (9000 BC - 2000 AD) + modern IJsselmeer polders
 *  - 30" ETOPO Holocene delta progradation (Mesopotamia/Persian Gulf, China Bohai/Yellow/Yangtze/PRD, Nile, Po)
 *  - 30" ETOPO dynamic basins (Caspian, Black Sea Neoeuxinian, Aral South+North, Lake Chad / Mega-Chad, Urmia)
 *  - Scrubbed natural lakes + dated GDW v1.0 reservoirs (zero pre-dam ghost footprints) + Toshka Lakes (1998)
 */
class E_AssembleOutputs {
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
   * Checks whether a 5' pixel falls inside a lat/lon bounding box.
   *
   * @param {number} arg0_pixel_idx
   * @param {Object} arg1_bbox
   *
   * @returns {boolean}
   */
  static inBbox (arg0_pixel_idx, arg1_bbox) {
    //Convert from parameters
    let idx = arg0_pixel_idx;
    let b = arg1_bbox;

    //Declare local instance variables
    let r = (idx/4320) | 0;
    let c = idx - r*4320;
    let lat = 90.0 - (r + 0.5)*(180.0/2160.0);
    let lon = -180.0 + (c + 0.5)*(360.0/4320.0);

    //Return statement
    return (lat >= b.lat_min && lat <= b.lat_max && lon >= b.lon_min && lon <= b.lon_max);
  }

  /**
   * Applies simulated Arctic Siberian surface water: each 5' cell's water fraction is interpolated between relative sea-level tiers using the regional sea level, itself interpolated in longitude between sector sample points (modern years use 0 m).
   *
   * @param {Float32Array} arg0_f_land
   * @param {Float32Array} arg1_f_ice
   * @param {number} arg2_year
   * @param {Object} arg3_hyd - Parsed arctic_hydrology.json.
   * @param {number} arg4_wetness - Palaeo-hydrology weight for this year (0 = modern mapped water only).
   */
  static applyArcticHydrology (arg0_f_land, arg1_f_ice, arg2_year, arg3_hyd, arg4_wetness) {
    //Convert from parameters
    let f_land = arg0_f_land;
    let f_ice = arg1_f_ice;
    let year = arg2_year;
    let hyd = arg3_hyd;
    let wetness = arg4_wetness;

    //Guard clauses
    if (!hyd || !hyd.cells || hyd.cells.length === 0 || !hyd.sectors || hyd.sectors.length === 0) return;
    if (wetness <= 0.0) return;

    //Declare local instance variables
    let cells = hyd.cells;
    let levels = hyd.rsl_levels_m;
    let last = levels.length - 1;
    let rsls_raw = (hyd.rsl_by_year && hyd.rsl_by_year[year]) ? hyd.rsl_by_year[year] : null;
    let sectors = hyd.sectors.map((sec, si) => ({ lon: sec.sample_lon, rsl: rsls_raw ? rsls_raw[si] : 0.0 })).sort((a, b) => a.lon - b.lon);

    for (let k = 0; k < cells.length; k++) {
      let row = cells[k];
      let p_idx = row[0];
      if (f_land[p_idx] <= 0.0 || f_ice[p_idx] >= 0.50) continue;
      let lon = -180.0 + ((p_idx % 4320) + 0.5)*(360.0/4320.0);
      let rsl = sectors[0].rsl;

      if (lon >= sectors[sectors.length - 1].lon) {
        rsl = sectors[sectors.length - 1].rsl;
      } else if (lon > sectors[0].lon) {
        for (let si = 0; si < sectors.length - 1; si++) {
          if (lon >= sectors[si].lon && lon <= sectors[si + 1].lon) {
            let t = (lon - sectors[si].lon)/(sectors[si + 1].lon - sectors[si].lon);
            rsl = sectors[si].rsl + t*(sectors[si + 1].rsl - sectors[si].rsl);
            break;
          }
        }
      }

      let frac = row[1 + last];
      if (rsl <= levels[0]) {
        frac = row[1];
      } else if (rsl < levels[last]) {
        for (let ti = 0; ti < last; ti++) {
          if (rsl >= levels[ti] && rsl <= levels[ti + 1]) {
            let t = (rsl - levels[ti])/(levels[ti + 1] - levels[ti]);
            frac = row[1 + ti]*(1.0 - t) + row[2 + ti]*t;
            break;
          }
        }
      }

      frac *= wetness;
      if (frac <= 0.005) continue;
      let cap = 1.0 - frac;
      if (f_land[p_idx] > cap) f_land[p_idx] = Math.max(0.0, cap);
    }
  }


  /**
   * Builds a cleaned static natural-lake mask that removes:
   *  - Any pixels covered by non-Lake-Control GDW v1.0 anthropogenic reservoirs (fixes HydroLAKES Lake_type==1 leaks)
   *  - Static Caspian, Aral, Chad, and Urmia polygons (which are driven dynamically by 30" ETOPO contours!)
   *  - Extra artificial spillway lakes (e.g. Toshka Lakes 1998, Lake Nasser residuals)
   *
   * @param {Float32Array} arg0_f_natlake_raw
   * @param {Float32Array} arg1_f_gdw_union
   * @param {Object} arg2_hcfg
   *
   * @returns {{ cleaned: Float32Array, extra_timed: Array<Object> }}
   */
  static buildCleanedNatLakeMask (arg0_f_natlake_raw, arg1_f_gdw_union, arg2_hcfg) {
    //Convert from parameters
    let raw = arg0_f_natlake_raw;
    let gdw_union = arg1_f_gdw_union;
    let hcfg = arg2_hcfg;

    //Declare local instance variables
    let cleaned = new Float32Array(raw);
    let extra_timed = [];
    let total_px = raw.length;

    //1. Scrub all non-Lake-Control GDW anthropogenic reservoirs from static natural lakes
    if (gdw_union) {
      for (let i = 0; i < total_px; i++) {
        if (gdw_union[i] > 0.0 && cleaned[i] > 0.0) {
          let nv = cleaned[i] - gdw_union[i];
          cleaned[i] = (nv > 0.0) ? nv : 0.0;
        }
      }
    }

    //2. Scrub basins driven dynamically by 30" ETOPO contours (Caspian, Aral, Chad, Urmia)
    for (let bi = 0; bi < hcfg.seeded_basins.length; bi++) {
      let b = hcfg.seeded_basins[bi];
      if (!b.exclude_from_static_natlakes) continue;
      let r0 = Math.max(0, Math.floor(((90.0 - b.bbox.lat_max)/180.0)*2160));
      let r1 = Math.min(2160, Math.ceil(((90.0 - b.bbox.lat_min)/180.0)*2160));
      let c0 = Math.max(0, Math.floor(((b.bbox.lon_min + 180.0)/360.0)*4320));
      let c1 = Math.min(4320, Math.ceil(((b.bbox.lon_max + 180.0)/360.0)*4320));

      for (let r = r0; r < r1; r++) {
        let off = r*4320;
        for (let c = c0; c < c1; c++)
          cleaned[off + c] = 0.0;
      }
    }

    //3. Extract extra artificial water bodies (e.g. Toshka Lakes 1998, Lake Nasser residuals 1964) into dated overlays
    let extras = hcfg.extra_artificial_water_excludes || [];
    for (let ei = 0; ei < extras.length; ei++) {
      let ex = extras[ei];
      let r0 = Math.max(0, Math.floor(((90.0 - ex.bbox.lat_max)/180.0)*2160));
      let r1 = Math.min(2160, Math.ceil(((90.0 - ex.bbox.lat_min)/180.0)*2160));
      let c0 = Math.max(0, Math.floor(((ex.bbox.lon_min + 180.0)/360.0)*4320));
      let c1 = Math.min(4320, Math.ceil(((ex.bbox.lon_max + 180.0)/360.0)*4320));
      let pairs = [];

      for (let r = r0; r < r1; r++) {
        let off = r*4320;
        for (let c = c0; c < c1; c++) {
          let idx = off + c;
          if (cleaned[idx] > 0.0) {
            pairs.push([idx, cleaned[idx]]);
            cleaned[idx] = 0.0;
          }
        }
      }

      extra_timed.push({
        bbox: ex.bbox,
        id: ex.id,
        pairs: pairs,
        start_year: ex.start_year,
      });
    }

    //Return statement
    return { cleaned: cleaned, extra_timed: extra_timed };
  }

  /**
   * Applies official RCE / Vos et al. (2021) Holocene Netherlands maps inside the Dutch coastal domain.
   * For -9000..1850: interpolates Vos's 12 official paleogeographical epochs.
   * For 1851..2025: transitions seamlessly from Vos 1850 -> OSM 2025 using exact real shorelines
   * gated by historical polder completion dates (Wieringermeer 1930, Noordoostpolder 1942,
   * Oostelijk Flevoland 1957, Zuidelijk Flevoland 1968, Maasvlakte 1965-2013, Lauwersmeer 1969).
   *
   * @param {Float32Array} arg0_f_land
   * @param {number} arg1_year
   * @param {Object} arg2_nl_vos
   * @param {Object} arg3_state - Shared assembly state (f_etopo_land_2025, f_land_osm_2025, min_bed_5m).
   */
  static applyVosNetherlands (arg0_f_land, arg1_year, arg2_nl_vos, arg3_state) {
    //Convert from parameters
    let f_land = arg0_f_land;
    let year = arg1_year;
    let nl = arg2_nl_vos;
    let st = arg3_state;

    //Guard clause: earlier than -9000 uses full GIA Doggerland emergence
    if (!nl || !nl.epochs || nl.epochs.length === 0) return;
    if (year < -9000 || year >= 2015) return;

    let epochs = nl.epochs;
    let dom = nl.domain_frac;
    let lw = nl.local_w;
    let lh = nl.local_h;
    let r0 = nl.r0;
    let c0 = nl.c0;

    if (year <= 1850) {
      let e0 = epochs[0];
      let e1 = epochs[0];
      for (let i = 0; i < epochs.length - 1; i++) {
        if (year >= epochs[i].year && year <= epochs[i + 1].year) {
          e0 = epochs[i];
          e1 = epochs[i + 1];
          break;
        }
      }

      let w1 = (e0.year === e1.year) ? 0.0 : (year - e0.year)/(e1.year - e0.year);
      let arr0 = e0.f_land;
      let arr1 = e1.f_land;
      let gate = this.getVosSpatialGate(nl, st);
      let min_bed = st.min_bed_5m;
      let rsl = this.evalCurve(st.vos_rsl_curve_m, year);

      for (let ly = 0; ly < lh; ly++) {
        let gy = r0 + ly;
        let g_off = gy*4320 + c0;
        let l_off = ly*lw;

        for (let lx = 0; lx < lw; lx++) {
          let d = dom[l_off + lx];
          let g = (gate[l_off + lx] >= 0.5) ? 1.0 : 0.0;
          if (d <= 0.05 || g <= 0.0) continue;

          //Each cell switches once: flooding when sea level plus ~2 m tidal range reaches its bed, emergence (peat and delta infill) high ground first; the switch completes by the next epoch
          let a0 = arr0[l_off + lx];
          let a1 = arr1[l_off + lx];
          let vos_land = a0;
          if (Math.abs(a1 - a0) > 0.02) {
            let z = min_bed[g_off + lx];
            let switched = (w1 >= 0.999);
            if (a1 < a0) {
              switched = switched || (rsl + 2.0 >= z);
            } else {
              let z_norm = Math.min(1.0, Math.max(0.0, (z + 25.0)/30.0));
              switched = switched || (w1 >= (1.0 - z_norm)*0.9 + 0.05);
            }
            vos_land = switched ? a1 : a0;
          }

          let cur = f_land[g_off + lx];
          let blended = cur*(1.0 - d*g) + vos_land*g;
          f_land[g_off + lx] = Math.min(1.0, Math.max(0.0, blended));
        }
      }
      return;
    }

    //1851..2014: Exact pixel-level transition from Vos 1850 to OSM 2025 gated by real polder bboxes
    let ep1850 = epochs.find((e) => e.year === 1850);
    if (!ep1850) return;
    let arr1850 = ep1850.f_land;

    let ramp = (arg0_s, arg1_e) => {
      if (year <= arg0_s) return 0.0;
      if (year >= arg1_e) return 1.0;
      return (year - arg0_s)/(arg1_e - arg0_s);
    };

    let p_wieringer = ramp(1927, 1930);
    let p_nop = ramp(1936, 1942);
    let p_oflevo = ramp(1950, 1957);
    let p_zflevo = ramp(1959, 1968);
    let p_maasvlakte = ramp(1965, 2013);
    let p_general = ramp(1850, 1970);

    for (let ly = 0; ly < lh; ly++) {
      let gy = r0 + ly;
      let lat = 90.0 - (gy + 0.5)*(180.0/2160.0);
      let g_off = gy*4320 + c0;
      let l_off = ly*lw;

      for (let lx = 0; lx < lw; lx++) {
        let d = dom[l_off + lx];
        if (d <= 0.2) continue;

        let lon = -180.0 + (c0 + lx + 0.5)*(360.0/4320.0);
        let v1850 = arr1850[l_off + lx];
        let v2025 = f_land[g_off + lx];
        if (Math.abs(v2025 - v1850) < 0.02) continue;

        let prog = p_general;
        if (lat >= 52.74 && lat <= 52.95 && lon >= 4.86 && lon <= 5.18) {
          prog = p_wieringer;
        } else if (lat >= 52.62 && lat <= 52.87 && lon >= 5.52 && lon <= 6.02) {
          prog = p_nop;
        } else if (lat >= 52.43 && lat <= 52.62 && lon >= 5.36 && lon <= 5.96) {
          prog = p_oflevo;
        } else if (lat >= 52.23 && lat <= 52.49 && lon >= 5.08 && lon <= 5.65) {
          prog = p_zflevo;
        } else if (lat >= 51.88 && lat <= 52.04 && lon >= 3.90 && lon <= 4.20) {
          prog = p_maasvlakte;
        }

        let target = v1850 + (v2025 - v1850)*prog;
        f_land[g_off + lx] = Math.min(1.0, Math.max(0.0, v2025*(1.0 - d) + target*d));
      }
    }
  }

  /**
   * Builds (once) the spatial weight for Vos et al. overrides: Vos governs modern land cells, the GIA model governs the open shelf and anything outside the Vos map footprint. The switch is hard (weight thresholded at 0,5) so the two sources never cross-fade into translucent ribbons, and the switch line follows the modern coast rather than the map extent.
   *
   * @param {Object} arg0_nl_vos
   * @param {Object} arg1_state
   *
   * @returns {Float32Array}
   */
  static getVosSpatialGate (arg0_nl_vos, arg1_state) {
    //Convert from parameters
    let nl = arg0_nl_vos;
    let st = arg1_state;

    //Guard clause: cached
    if (nl._gate) return nl._gate;

    //Declare local instance variables
    let dom = nl.domain_frac;
    let gate = new Float32Array(nl.local_w*nl.local_h);
    let lh = nl.local_h;
    let lw = nl.local_w;
    let land_mask = new Uint8Array(lw*lh);
    let out_mask = new Uint8Array(lw*lh);

    for (let ly = 0; ly < lh; ly++)
      for (let lx = 0; lx < lw; lx++) {
        let l = ly*lw + lx;
        let g = (nl.r0 + ly)*4320 + nl.c0 + lx;
        if (st.f_etopo_land_2025[g] >= 0.5 || st.f_land_osm_2025[g] >= 0.5) land_mask[l] = 1;
        if (dom[l] <= 0.05 || ly === 0 || lx === 0 || ly === lh - 1 || lx === lw - 1) out_mask[l] = 1;
      }

    let d_coast = BasinContourEngine.computeChamferDistance(land_mask, lw, lh);
    let d_edge = BasinContourEngine.computeChamferDistance(out_mask, lw, lh);

    for (let l = 0; l < gate.length; l++) {
      //Exclude narrow barrier-island cells (fewer than 3 land neighbours) so they follow the GIA shelf rather than speckling it
      let ly = (l/lw) | 0;
      let lx = l - ly*lw;
      let nbrs = 0;
      if (lx > 0) nbrs += land_mask[l - 1];
      if (lx < lw - 1) nbrs += land_mask[l + 1];
      if (ly > 0) nbrs += land_mask[l - lw];
      if (ly < lh - 1) nbrs += land_mask[l + lw];
      let g_coast = (d_coast[l] <= 0.5 && nbrs >= 3) ? 1.0 : 0.0;
      let g_edge = Math.min(1.0, Math.max(0.0, d_edge[l]/3.0));
      gate[l] = g_coast*g_edge;
    }
    nl._gate = gate;

    //Return statement
    return gate;
  }


  /**
   * Composites the signed float32 km2 grid for a single target year.
   *
   * @param {number} arg0_year
   * @param {Object} arg1_state
   *
   * @returns {Float32Array}
   */
  static composeYearSlice (arg0_year, arg1_state) {
    //Convert from parameters
    let year = arg0_year;
    let st = arg1_state;

    //Declare local instance variables
    let cfg = st.cfg;
    let delta_coast = null;
    let f_ice = st.f_ice_2025;
    let f_land = new Float32Array(st.f_land_osm_2025);
    let total_px = cfg.grid.width*cfg.grid.height;

    //1. Apply paleo GIA coastal anomaly & paleo glaciation if year < gia_cutoff_year
    if (year < cfg.paleo.gia_cutoff_year) {
      let coast_path = path.join(st.paleo_dir, `delta_coast_${year}.f32`);
      let ice_path = path.join(st.paleo_dir, `f_ice_${year}.f32`);

      if (fs.existsSync(coast_path)) {
        delta_coast = GeoGridUtils.readFloat32Bin(coast_path);
        let f_etopo_2025 = st.f_etopo_land_2025;
        for (let i = 0; i < total_px; i++) {
          let dc = delta_coast[i];
          if (dc === 0.0) continue;

          let osm = f_land[i];
          let et = f_etopo_2025[i];
          if (dc > 0.0) {
            let denom = 1.0 - et;
            let gain = (denom > 1e-4) ? Math.min(1.0, dc/denom) : 1.0;
            f_land[i] = osm + (1.0 - osm)*gain;
          } else {
            let loss = (et > 1e-4) ? Math.min(1.0, (-dc)/et) : Math.min(1.0, -dc);
            f_land[i] = osm*(1.0 - loss);
          }
        }
      }

      if (fs.existsSync(ice_path))
        f_ice = GeoGridUtils.readFloat32Bin(ice_path);
    }

    //2. Apply Black Sea Neoeuxinian isolated phase override (-10000..-5000) so Sea of Azov & NW shelf dry accurately
    let bs = st.dem_basins.seeded.black_sea_neoeuxinian;
    if (bs && bs.year_to_level_id[year] !== undefined) {
      let r0 = Math.max(0, Math.floor(((90.0 - bs.bbox.lat_max)/180.0)*2160));
      let r1 = Math.min(2160, Math.ceil(((90.0 - bs.bbox.lat_min)/180.0)*2160));
      let c0 = Math.max(0, Math.floor(((bs.bbox.lon_min + 180.0)/360.0)*4320));
      let c1 = Math.min(4320, Math.ceil(((bs.bbox.lon_max + 180.0)/360.0)*4320));
      let lid = bs.year_to_level_id[year];
      let water_pairs = bs.level_table[lid].pairs;

      //Inside Black Sea basin during Neoeuxinian phase, set basin terrestrial then carve exact Neoeuxinian lake contour
      //Only override modern Black Sea / Sea of Azov cells (where f_etopo_land_2025 < 1.0)
      for (let r = r0; r < r1; r++) {
        let off = r*4320;
        for (let c = c0; c < c1; c++) {
          let idx = off + c;
          if (st.f_etopo_land_2025[idx] < 1.0)
            f_land[idx] = 1.0;
        }
      }
      for (let k = 0; k < water_pairs.length; k++) {
        let p_idx = water_pairs[k][0];
        let wf = water_pairs[k][1];
        f_land[p_idx] = Math.max(0.0, 1.0 - wf);
      }
    }

    //3. Apply 30" ETOPO & SDF-migrating Holocene alluvial delta progradation (ZERO ghost alpha-fade!)
    let delta_keys = Object.keys(st.dem_basins.deltas);
    for (let di = 0; di < delta_keys.length; di++) {
      let d = st.dem_basins.deltas[delta_keys[di]];
      let lid = d.year_to_level_id[year];
      if (lid === undefined) continue;
      let pairs = d.level_table[lid].pairs;
      let is_prop = (d.mode === "proportional");

      for (let k = 0; k < pairs.length; k++) {
        let pair = pairs[k];
        let p_idx = pair[0];
        if (is_prop && pair.length >= 3) {
          let f_new = pair[1];
          let loss_zero = pair[2];
          let et = st.f_etopo_land_2025[p_idx];
          if (et > 1e-4) {
            if (f_new > 0.0) {
              let loss = Math.min(1.0, f_new/et);
              f_land[p_idx] *= (1.0 - loss);
            }
          } else if (loss_zero > 0.0) {
            f_land[p_idx] *= (1.0 - loss_zero);
          }
        } else {
          let nv = f_land[p_idx] - pair[1];
          f_land[p_idx] = (nv > 0.0) ? nv : 0.0;
        }
      }
    }

    //4. Apply official RCE Vos et al. (2021) Atlas of the Holocene Netherlands (-9000..1926)
    this.applyVosNetherlands(f_land, year, st.nl_vos, st);

    //5. Apply modern/historic coastal reclamation (post-1927 NL Flevoland/IJsselmeer/Maasvlakte, Tokyo, Singapore, HK, China, Dubai, Bahrain, etc.)
    for (let i = 0; i < st.reclamation.length; i++) {
      let ev = st.reclamation[i];
      let progress = 1.0;
      if (year <= ev.start_year) {
        progress = 0.0;
      } else if (year < ev.end_year) {
        progress = (year - ev.start_year)/(ev.end_year - ev.start_year);
      }

      if (ev.action === "reclaimed_from_water" && progress < 1.0) {
        //Skip all nl_* polygons here because applyVosNetherlands now handles exact Vos 1850 -> OSM 2025 polder dates!
        if (ev.id.startsWith("nl_")) continue;

        let remove_factor = 1.0 - progress;
        for (let x = 0; x < ev.pairs.length; x++) {
          let p_idx = ev.pairs[x][0];
          //Scale subtraction by modern marine exposure (1 - gia_gain) so emergent Pleistocene seabed (e.g. Bahrain/Dubai at 10000 BC) is never carved into desert holes!
          let gia_scale = 1.0;
          if (delta_coast && delta_coast[p_idx] > 0.0) {
            let denom = 1.0 - st.f_etopo_land_2025[p_idx];
            let gain = (denom > 1e-4) ? Math.min(1.0, delta_coast[p_idx]/denom) : 1.0;
            gia_scale = Math.max(0.0, 1.0 - gain);
          }
          if (ev.pairs[x][1] >= 0.85 || st.f_etopo_land_2025[p_idx] < 0.05) {
            f_land[p_idx] *= Math.max(0.0, 1.0 - remove_factor*gia_scale);
          } else {
            let delta = ev.pairs[x][1]*remove_factor*gia_scale;
            let nv = f_land[p_idx] - delta;
            f_land[p_idx] = (nv > 0.0) ? nv : 0.0;
          }
        }
      } else if ((ev.action === "submerged_from_land" || ev.action === "canal_excavated") && progress < 1.0) {
        let restore_factor = 1.0 - progress;
        for (let x = 0; x < ev.pairs.length; x++) {
          let p_idx = ev.pairs[x][0];
          let delta = ev.pairs[x][1]*restore_factor;
          let nv = f_land[p_idx] + delta;
          f_land[p_idx] = (nv < 1.0) ? nv : 1.0;
        }
      }
    }

    //6. Subtract cleaned static natural lakes (with GDW reservoirs, Caspian/Aral/Chad/Urmia, and Toshka scrubbed!)
    //   R3-D03: On continental interior pixels where a natural lake covers >=50% of the 5' cell, snap to crisp open water (0.0)
    let f_natlake = st.f_natlake_cleaned;
    let f_osm_ref = st.f_land_osm_2025;
    for (let i = 0; i < total_px; i++) {
      let lk = f_natlake[i];
      if (lk > 0.0 && f_land[i] > 0.0) {
        if (lk >= 0.50 && f_osm_ref[i] >= 0.95) {
          f_land[i] = 0.0;
        } else {
          let nv = f_land[i] - lk;
          f_land[i] = (nv > 0.0) ? nv : 0.0;
        }
      }
    }

    //7. Apply extra dated artificial spillway lakes (e.g. Toshka Lakes 1998, Lake Nasser residuals 1964)
    for (let ei = 0; ei < st.extra_timed_lakes.length; ei++) {
      let ex = st.extra_timed_lakes[ei];
      if (year >= ex.start_year) {
        for (let k = 0; k < ex.pairs.length; k++) {
          let p_idx = ex.pairs[k][0];
          let nv = f_land[p_idx] - ex.pairs[k][1];
          f_land[p_idx] = (nv > 0.0) ? nv : 0.0;
        }
      }
    }

    //8. Apply 30" ETOPO dynamic closed basins (Caspian, Aral South+North, Lake Chad / Mega-Chad, Lake Urmia)
    let seeded_keys = Object.keys(st.dem_basins.seeded);
    for (let si = 0; si < seeded_keys.length; si++) {
      let b = st.dem_basins.seeded[seeded_keys[si]];
      if (b.is_marine_override) continue;
      let lid = b.year_to_level_id[year];
      if (lid === undefined) continue;
      let pairs = b.level_table[lid].pairs;
      for (let k = 0; k < pairs.length; k++) {
        let p_idx = pairs[k][0];
        let nv = f_land[p_idx] - pairs[k][1];
        f_land[p_idx] = (nv > 0.0) ? nv : 0.0;
      }
    }

    //8b. Apply ice-dammed proglacial lakes flooded on the 30" DEM (Lake Agassiz-Ojibway)
    let pg_keys = Object.keys(st.proglacial);
    for (let pi = 0; pi < pg_keys.length; pi++) {
      let pairs = st.proglacial[pg_keys[pi]].by_year[year];
      if (!pairs) continue;
      for (let k = 0; k < pairs.length; k++) {
        let p_idx = pairs[k][0];
        let nv = f_land[p_idx] - pairs[k][1];
        f_land[p_idx] = (nv > 0.0) ? nv : 0.0;
      }
    }

    //8c. Apply dynamic sill-controlled coastal lagoons (Garabogazköl: Caspian lowstands + 1980-1992 Soviet dam)
    let shore_band_lag = 0.03;
    for (let li = 0; li < st.lagoons.length; li++) {
      let lg = st.lagoons[li];
      let water_mult = this.evalCurve(lg.level_curve, year);
      if (water_mult <= 0.0) continue;

      for (let k = 0; k < lg.pairs.length; k++) {
        let p_idx = lg.pairs[k][0];
        let poly_frac = lg.pairs[k][1];
        let rank_norm = lg.pairs[k][2];
        let local_water = 0.0;

        if (water_mult >= 0.999 || rank_norm <= water_mult) {
          local_water = 1.0;
        } else if (rank_norm < water_mult + shore_band_lag) {
          local_water = 1.0 - (rank_norm - water_mult)/shore_band_lag;
        }

        let eff_sub = poly_frac*local_water;
        if (eff_sub >= 0.50 && st.f_land_osm_2025[p_idx] >= 0.95) {
          f_land[p_idx] = 0.0;
        } else if (eff_sub > 0.0) {
          let nv = f_land[p_idx] - eff_sub;
          f_land[p_idx] = (nv > 0.0) ? nv : 0.0;
        }
      }
    }

    //8d. Apply simulated Arctic Siberian surface water (channels, floodplains, fill-and-spill lakes) blended by local relative sea level (overlap-safe: water = max(existing water, simulated water))
    this.applyArcticHydrology(f_land, f_ice, year, st.arctic_hydrology, this.evalCurve(st.hydrology_wetness_curve, year));

    //8e. Persistent palaeo-estuary water (lower Ob, Taz, Yenisei): never drier than today inland of the modern estuary mouths
    let persist = (st.arctic_hydrology && st.arctic_hydrology.persistent_cells) ? st.arctic_hydrology.persistent_cells : [];
    for (let k = 0; k < persist.length; k++) {
      let p_idx = persist[k];
      if (f_land[p_idx] > st.f_land_osm_2025[p_idx]) f_land[p_idx] = st.f_land_osm_2025[p_idx];
    }

    //9. Apply dated GDW v1.0 reservoirs (strictly 0 subtraction before impound_year -> ZERO anachronistic ghosts!)
    let res = st.reservoirs;
    if (res && res.pixel_indices) {
      let p_arr = res.pixel_indices;
      let f_arr = res.fracs;
      let y_arr = res.impound_years;
      let suppressed_boxes = [];
      for (let ei = 0; ei < st.extra_timed_lakes.length; ei++) {
        if (year < st.extra_timed_lakes[ei].start_year && st.extra_timed_lakes[ei].bbox)
          suppressed_boxes.push(st.extra_timed_lakes[ei].bbox);
      }

      for (let i = 0; i < p_arr.length; i++) {
        if (year < y_arr[i]) continue;
        let p_idx = p_arr[i];
        if (suppressed_boxes.length > 0) {
          let skip = false;
          for (let bi = 0; bi < suppressed_boxes.length; bi++) {
            if (this.inBbox(p_idx, suppressed_boxes[bi])) {
              skip = true;
              break;
            }
          }
          if (skip) continue;
        }
        if (f_land[p_idx] > 0.0) {
          //R3-D03: Majority-water reservoir core pixels (>=45% impounded water on continental interior)
          //snap cleanly to open water so narrow-channel shield reservoirs (Manicouagan ring, Smallwood,
          //Caniapiscau, Robert-Bourassa) render as cohesive open-water basins after impoundment!
          if (f_arr[i] >= 0.45 && f_osm_ref[p_idx] >= 0.95) {
            f_land[p_idx] = 0.0;
          } else {
            let nv = f_land[p_idx] - f_arr[i];
            f_land[p_idx] = (nv > 0.0) ? nv : 0.0;
          }
        }
      }
    }

    //10. Apply remaining hypsometric proglacial/historic paleolakes (Lake Agassiz-Ojibway, Texcoco, Tulare, Lop Nur, Salton)
    let shore_band = 0.03;
    for (let i = 0; i < st.paleolakes.length; i++) {
      let lk = st.paleolakes[i];
      let water_mult = this.evalCurve(lk.curve, year);

      if (lk.action === "historical_water_body" && water_mult > 0.0) {
        for (let x = 0; x < lk.pairs.length; x++) {
          let p_idx = lk.pairs[x][0];
          let poly_frac = lk.pairs[x][1];
          let rank_norm = (lk.pairs[x].length >= 3) ? lk.pairs[x][2] : 0.0;
          let local_water = 0.0;

          if (water_mult >= 0.999 || rank_norm <= water_mult) {
            local_water = 1.0;
          } else if (rank_norm < water_mult + shore_band) {
            local_water = 1.0 - (rank_norm - water_mult)/shore_band;
          }

          if (local_water > 0.0) {
            let nv = f_land[p_idx] - poly_frac*local_water;
            f_land[p_idx] = (nv > 0.0) ? nv : 0.0;
          }
        }
      } else if (lk.action === "modern_lake_formed" && water_mult < 1.0) {
        let dry_mult = 1.0 - water_mult;
        for (let x = 0; x < lk.pairs.length; x++) {
          let p_idx = lk.pairs[x][0];
          let poly_frac = lk.pairs[x][1];
          let rank_norm = (lk.pairs[x].length >= 3) ? lk.pairs[x][2] : 0.0;
          let local_dry = 0.0;

          if (dry_mult >= 0.999 || rank_norm >= water_mult) {
            local_dry = 1.0;
          } else if (rank_norm > water_mult - shore_band) {
            local_dry = 1.0 - (water_mult - rank_norm)/shore_band;
          }

          if (local_dry > 0.0) {
            let nv = f_land[p_idx] + poly_frac*local_dry;
            f_land[p_idx] = (nv < 1.0) ? nv : 1.0;
          }
        }
      }
    }

    //10b. Apply explicit historic narrow land bridges & reef island chains (e.g. Adam's Bridge)
    for (let li = 0; li < st.land_bridges.length; li++) {
      let lb = st.land_bridges[li];
      let breach_t = 0.0;
      if (year >= lb.breach_end_year) {
        breach_t = 1.0;
      } else if (year > lb.breach_start_year) {
        breach_t = (year - lb.breach_start_year)/(lb.breach_end_year - lb.breach_start_year);
      }
      for (let ci = 0; ci < lb.pairs.length; ci++) {
        let p_idx = lb.pairs[ci][0];
        let min_pre = lb.pairs[ci][1];
        let min_post = lb.pairs[ci][2];
        let target_min = min_pre + (min_post - min_pre)*breach_t;
        if (target_min > f_land[p_idx])
          f_land[p_idx] = target_min;
      }
    }

    //11. Subtract glacial cover and apply strict signed precedence (Land > Glaciated + Permanent Ice Shelves > Sea)
    //    Includes D06b (zero false green strip along Antarctic floating ice-shelf calving fronts),
    //    D06c (zero false blue interior holes in Transantarctic/Ellsworth/Prince Charles Mountains),
    //    and D19 (majority-glaciated alpine peaks like NZ Southern Alps render as crisp glaciated terrain).
    let f_ice_eff = new Float32Array(f_ice);
    let f_etopo_2025 = st.f_etopo_land_2025;
    let f_osm_2025 = st.f_land_osm_2025;
    let width = cfg.grid.width;
    let height = cfg.grid.height;

    for (let y = 0; y < height; y++) {
      let lat_deg = 90.0 - (y + 0.5)*(180.0/height);
      let row_off = y*width;
      let is_ant_shelf_zone = (lat_deg <= -62.0);
      let is_ant_interior_zone = (lat_deg <= -65.0);

      for (let x = 0; x < width; x++) {
        let i = row_off + x;
        let ice = f_ice_eff[i];

        //D06b: Antarctic floating ice shelves & calving fronts (Ross, Ronne-Filchner, Amery, Larsen) -> never paint green land
        if (is_ant_shelf_zone && f_etopo_2025[i] < 0.50 && ice > 0.01) {
          f_land[i] = 0.0;
          continue;
        }
        //Also catch 1-px outer OSM shelf-edge pixels right along Ross & Ronne-Filchner fronts where ice > 0 or neighbours are shelf
        if (lat_deg <= -74.0 && f_etopo_2025[i] === 0.0 && f_land[i] > 0.0) {
          let lon_deg = -180.0 + (x + 0.5)*(360.0/width);
          let in_ross = (lat_deg <= -76.5 && (lon_deg >= 162.0 || lon_deg <= -150.0));
          let in_ronne = (lat_deg <= -74.0 && lon_deg >= -65.0 && lon_deg <= -25.0);
          if (in_ross || in_ronne) {
            if (ice < f_land[i]) f_ice_eff[i] = f_land[i];
            f_land[i] = 0.0;
            continue;
          }
        }

        //D06c: Antarctic continental interior (Transantarctic, Ellsworth, Prince Charles Mtns) -> never leave blue holes
        if (is_ant_interior_zone && f_osm_2025[i] >= 0.90) {
          if (ice >= 0.50) {
            f_land[i] = 0.0;
            if (f_osm_2025[i] > f_ice_eff[i]) f_ice_eff[i] = f_osm_2025[i];
          } else {
            f_land[i] = Math.max(f_land[i], f_osm_2025[i] - ice);
          }
          continue;
        }

        //R2-D03b: Grounded continental ice sheets (f_ice >= 0.50) strictly override inland lakes (f_natlake)
        if (ice >= 0.50 && f_land[i] <= 0.0 && f_osm_2025[i] >= 0.50) {
          f_ice_eff[i] = 1.0;
          continue;
        }

        if (ice > 0.0 && f_land[i] > 0.0) {
          //D19: Majority-glaciated high-alpine / polar ice cells (>=50% glacier; >=22% along narrow NZ Southern Alps crest)
          let lon_deg = -180.0 + (x + 0.5)*(360.0/width);
          let in_nz_alps = (lat_deg >= -46.0 && lat_deg <= -42.5 && lon_deg >= 167.0 && lon_deg <= 171.8);
          let gl_thresh = in_nz_alps ? 0.22 : 0.50;
          if (ice >= gl_thresh && ice >= f_land[i]*gl_thresh) {
            f_land[i] = 0.0;
            if (f_osm_2025[i] >= 0.95 && ice < 0.85) f_ice_eff[i] = Math.max(ice, 0.85);
          } else {
            let bare = f_land[i] - ice;
            f_land[i] = (bare > 0.0) ? bare : 0.0;
          }
        }
      }
    }

    let signed_km2 = GeoGridUtils.applySignedPrecedence(f_land, f_ice_eff, st.cell_area_by_row, {
      epsilon: cfg.encoding.min_fraction_epsilon,
      height: cfg.grid.height,
      width: cfg.grid.width,
    });

    //Return statement
    return signed_km2;
  }

  /**
   * Executes Step E: generates all 218 GeoPNGs (rasters/land_area_<year>.png) and landarea_Hollinger.nc.
   *
   * @returns {Promise<void>}
   */
  static async run () {
    //Declare local instance variables
    let cfg = GeoGridUtils.loadConfig("model_config.json5");
    let hcfg = GeoGridUtils.loadConfig("hydro_geomorph_basins.json5");
    let cell_area_by_row = GeoGridUtils.computeCellAreaByRow(cfg.grid);
    let coords = GeoGridUtils.getGridCoordinates(cfg.grid);
    let paths = GeoGridUtils.resolvePaths("landarea");
    let years = GeoGridUtils.getTimestepYears(cfg);
    let force = process.argv.includes("--force");

    if (!fs.existsSync(paths.rasters_dir))
      fs.mkdirSync(paths.rasters_dir, { recursive: true });

    let f_natlake_raw = GeoGridUtils.readFloat32Bin(path.join(paths.intermediate_dir, "f_natlake_2025.f32"));
    let gdw_union_path = path.join(paths.intermediate_dir, "f_gdw_reservoirs_union.f32");
    let f_gdw_union = fs.existsSync(gdw_union_path) ? GeoGridUtils.readFloat32Bin(gdw_union_path) : null;
    let { cleaned, extra_timed } = this.buildCleanedNatLakeMask(f_natlake_raw, f_gdw_union, hcfg);

    let lagoons_path = path.join(paths.intermediate_dir, "lagoon_basins_sparse.json");
    let hydrology_path = path.join(paths.intermediate_dir, "arctic_hydrology.json");
    let land_bridges = (hcfg.land_bridges || []).map((lb) => ({
      breach_end_year: lb.breach_end_year,
      breach_start_year: lb.breach_start_year,
      id: lb.id,
      pairs: lb.cells.map(([lat, lon, min_pre, min_post]) => {
        let gy = Math.min(2159, Math.max(0, Math.floor(((90.0 - lat)/180.0)*2160)));
        let gx = Math.min(4319, Math.max(0, Math.floor(((lon + 180.0)/360.0)*4320)));
        return [gy*4320 + gx, min_pre, min_post];
      }),
    }));

    let state = {
      cell_area_by_row: cell_area_by_row,
      cfg: cfg,
      dem_basins: JSON.parse(fs.readFileSync(path.join(paths.intermediate_dir, "dem_basins_by_year.json"), "utf8")),
      extra_timed_lakes: extra_timed,
      f_etopo_land_2025: GeoGridUtils.readFloat32Bin(path.join(paths.intermediate_dir, "f_etopo_land_2025.f32")),
      f_ice_2025: GeoGridUtils.readFloat32Bin(path.join(paths.intermediate_dir, "f_ice_2025.f32")),
      f_land_osm_2025: GeoGridUtils.readFloat32Bin(path.join(paths.intermediate_dir, "f_land_osm_2025.f32")),
      f_natlake_cleaned: cleaned,
      lagoons: fs.existsSync(lagoons_path) ? JSON.parse(fs.readFileSync(lagoons_path, "utf8")) : [],
      land_bridges: land_bridges,
      min_bed_5m: GeoGridUtils.readFloat32Bin(path.join(paths.intermediate_dir, "etopo_min_bed_5m.f32")),
      nl_vos: JSON.parse(fs.readFileSync(path.join(paths.intermediate_dir, "nl_vos_epochs.json"), "utf8")),
      paleo_dir: path.join(paths.intermediate_dir, "paleo_slices"),
      arctic_hydrology: fs.existsSync(hydrology_path) ? JSON.parse(fs.readFileSync(hydrology_path, "utf8")) : null,
      hydrology_wetness_curve: (hcfg.arctic_hydrology && hcfg.arctic_hydrology.palaeo_wetness_curve) ? hcfg.arctic_hydrology.palaeo_wetness_curve : [[-10000, 1.0]],
      paleolakes: JSON.parse(fs.readFileSync(path.join(paths.intermediate_dir, "paleolake_sparse.json"), "utf8")),
      proglacial: JSON.parse(fs.readFileSync(path.join(paths.intermediate_dir, "proglacial_lakes_by_year.json"), "utf8")),
      reclamation: JSON.parse(fs.readFileSync(path.join(paths.intermediate_dir, "reclamation_sparse.json"), "utf8")),
      reservoirs: JSON.parse(fs.readFileSync(path.join(paths.intermediate_dir, "reservoirs_sparse.json"), "utf8")),
      vos_rsl_curve_m: hcfg.rce_netherlands_vos.transition_rsl_curve_m,
    };

    let nc_path = paths.nc_path;
    let prefix = cfg.encoding.geopng_prefix;
    let t0 = Date.now();
    let write_queue = [];

    console.log(`[E_assemble] Writing ${years.length} signed float32 GeoPNGs & landarea.nc ...`);

    await NetCDF4Writer.writeHollingerNetCDF4(
      nc_path,
      {
        cell_area_km2: cell_area_by_row,
        chunk_shape: cfg.encoding.netcdf_chunk_shape,
        compression_level: cfg.encoding.netcdf_compression_level,
        lat: coords.lat,
        lon: coords.lon,
        years: years,
      },
      async (arg0_year, arg1_t_idx) => {
        let year = arg0_year;
        let t_idx = arg1_t_idx;
        let png_path = path.join(paths.rasters_dir, `${prefix}${year}.png`);
        let slice;

        if (!force && fs.existsSync(png_path)) {
          slice = GeoGridUtils.decodeFloat32GeoPNG(png_path).data;
        } else {
          slice = this.composeYearSlice(year, state);
          let p = GeoGridUtils.encodeFloat32GeoPNG(png_path, slice, {
            compression_level: 6,
            height: cfg.grid.height,
            width: cfg.grid.width,
          });
          write_queue.push(p);
          if (write_queue.length >= 6) {
            await Promise.all(write_queue);
            write_queue.length = 0;
          }
        }

        if ((t_idx + 1) % 20 === 0 || t_idx === 0 || t_idx === years.length - 1) {
          let land_km2 = 0.0;
          let ice_km2 = 0.0;
          for (let i = 0; i < slice.length; i++) {
            let v = slice[i];
            if (v > 0.0) land_km2 += v;
            else if (v < 0.0) ice_km2 -= v;
          }
          console.log(`[E_assemble]   [${t_idx + 1}/${years.length}] year=${year}: land=${(land_km2/1e6).toFixed(2)} M km2, glaciated=${(ice_km2/1e6).toFixed(2)} M km2 (${((Date.now() - t0)/1000).toFixed(1)}s)`);
        }

        return slice;
      }
    );

    if (write_queue.length > 0)
      await Promise.all(write_queue);

    let nc_mb = (fs.statSync(nc_path).size/1e6).toFixed(1);
    console.log(`[E_assemble] Completed landarea.nc (${nc_mb} MB) & ${years.length} GeoPNGs in ${((Date.now() - t0)/1000).toFixed(1)}s.`);
  }
}

if (require.main === module) {
  E_AssembleOutputs.run().catch((err) => {
    console.error("[E_assemble] Fatal error:", err);
    process.exit(1);
  });
}

module.exports = E_AssembleOutputs;
