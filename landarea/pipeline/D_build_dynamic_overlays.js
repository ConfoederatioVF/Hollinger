"use strict";

let fs = require("node:fs");
let path = require("node:path");

let shapefile = require("shapefile");

let ArcticHydrologyModel = require("../../core/ArcticHydrologyModel");
let BasinContourEngine = require("../../core/BasinContourEngine");
let C_BuildPaleoGiaIce = require("./C_build_paleo_gia_ice");
let GeoGridUtils = require("../../core/GeoGridUtils");
let NetCDF4Writer = require("../../core/NetCDF4Writer");
let ScanlineRasteriser = require("../../core/ScanlineRasteriser");

/**
 * Step D: Builds all dynamic terrestrial/hydrographic overlays:
 *  1. Global Dam Watch (GDW) v1.0 dated reservoirs (excluding regulated natural lakes) + union mask to scrub HydroLAKES leaks
 *  2. Official RCE / Vos et al. (2021) Atlas of the Holocene Netherlands (13 epochs, -9000..2000) + exact modern IJsselmeer polder masks
 *  3. 30" ETOPO topographic-contour closed basins (Caspian, Black Sea Neoeuxinian, Aral South+North, Lake Chad / Mega-Chad, Lake Urmia)
 *  4. 30" ETOPO Holocene prograding alluvial deltas (Mesopotamia/Persian Gulf, China Bohai/Yellow/Yangtze/PRD, Nile, Po)
 *  5. Remaining curated coastal reclamation & hypsometric proglacial/historic paleolakes (Lake Agassiz, Texcoco, Tulare, Lop Nur, Salton)
 */
class D_BuildDynamicOverlays {
  /**
   * Finds the reservoir polygon shapefile inside raw_sources/gdw_v10.
   *
   * @param {string} arg0_dir
   *
   * @returns {string|null}
   */
  static findReservoirShapefile (arg0_dir) {
    //Convert from parameters
    let dir = arg0_dir;

    //Guard clause
    if (!fs.existsSync(dir)) return null;

    //Declare local instance variables
    let all_shps = [];
    let walk = (arg0_current) => {
      let entries = fs.readdirSync(arg0_current, { withFileTypes: true });
      for (let i = 0; i < entries.length; i++) {
        let full = path.join(arg0_current, entries[i].name);
        if (entries[i].isDirectory()) {
          walk(full);
        } else if (entries[i].name.toLowerCase().endsWith(".shp")) {
          all_shps.push(full);
        }
      }
    };
    walk(dir);

    let res_shp = all_shps.find((p) => /reservoir/i.test(path.basename(p)));

    //Return statement
    return res_shp ? res_shp : (all_shps[0] || null);
  }

  /**
   * Rasterises a single GeoJSON geometry into compact sparse [pixel_idx, frac] pairs.
   *
   * @param {Object} arg0_geometry
   * @param {Float64Array} arg1_cell_area_by_row
   * @param {number} [arg2_supersample=20]
   *
   * @returns {Array<Array<number>>}
   */
  static rasteriseGeometrySparse (arg0_geometry, arg1_cell_area_by_row, arg2_supersample) {
    //Convert from parameters
    let geometry = arg0_geometry;
    let cell_area_by_row = arg1_cell_area_by_row;
    let supersample = (arg2_supersample !== undefined) ? arg2_supersample : 20;

    //Guard clause
    if (!geometry || !geometry.coordinates) return [];

    let max_lat = -90.0;
    let max_lon = -180.0;
    let min_lat = 90.0;
    let min_lon = 180.0;

    let scanRing = (arg0_ring) => {
      for (let i = 0; i < arg0_ring.length; i++) {
        let lon = arg0_ring[i][0];
        let lat = arg0_ring[i][1];
        if (lon < min_lon) min_lon = lon;
        if (lon > max_lon) max_lon = lon;
        if (lat < min_lat) min_lat = lat;
        if (lat > max_lat) max_lat = lat;
      }
    };

    let polys = (geometry.type === "Polygon") ? [geometry.coordinates] : (geometry.type === "MultiPolygon" ? geometry.coordinates : []);
    for (let i = 0; i < polys.length; i++)
      for (let x = 0; x < polys[i].length; x++)
        scanRing(polys[i][x]);

    if (min_lon >= max_lon || min_lat >= max_lat) return [];

    let height = 2160;
    let width = 4320;
    let r0 = Math.max(0, Math.floor(((90.0 - max_lat)/180.0)*height));
    let r1 = Math.min(height - 1, Math.ceil(((90.0 - min_lat)/180.0)*height));
    let c0 = Math.max(0, Math.floor(((min_lon + 180.0)/360.0)*width));
    let c1 = Math.min(width - 1, Math.ceil(((max_lon + 180.0)/360.0)*width));
    let local_h = r1 - r0 + 1;
    let local_w = c1 - c0 + 1;
    let sub_h = local_h*supersample;
    let sub_w = local_w*supersample;
    let local_mask = new Uint8Array(sub_h*sub_w);
    let scale_x = (width*supersample)/360.0;
    let scale_y = (height*supersample)/180.0;
    let origin_sx = c0*supersample;
    let origin_sy = r0*supersample;
    let total_hits = 0;

    for (let p = 0; p < polys.length; p++) {
      let rings = polys[p];
      let edge_x0 = [], edge_y0 = [], edge_x1 = [], edge_y1 = [];
      let p_min_sy = sub_h;
      let p_max_sy = -1;

      for (let r = 0; r < rings.length; r++) {
        let ring = rings[r];
        if (!ring || ring.length < 3) continue;
        for (let k = 0; k < ring.length - 1; k++) {
          let sy0 = (90.0 - ring[k][1])*scale_y - origin_sy;
          let sy1 = (90.0 - ring[k + 1][1])*scale_y - origin_sy;
          if (sy0 === sy1) continue;
          let sx0 = (ring[k][0] + 180.0)*scale_x - origin_sx;
          let sx1 = (ring[k + 1][0] + 180.0)*scale_x - origin_sx;
          if (sy0 < p_min_sy) p_min_sy = sy0;
          if (sy1 < p_min_sy) p_min_sy = sy1;
          if (sy0 > p_max_sy) p_max_sy = sy0;
          if (sy1 > p_max_sy) p_max_sy = sy1;
          edge_x0.push(sx0);
          edge_y0.push(sy0);
          edge_x1.push(sx1);
          edge_y1.push(sy1);
        }
      }

      let n_edges = edge_x0.length;
      if (n_edges === 0) continue;

      let sy_start = Math.max(0, Math.floor(p_min_sy));
      let sy_end = Math.min(sub_h - 1, Math.ceil(p_max_sy));
      let xs_list = [];

      for (let sy = sy_start; sy <= sy_end; sy++) {
        let scan_y = sy + 0.5;
        xs_list.length = 0;
        for (let i = 0; i < n_edges; i++) {
          let y0 = edge_y0[i];
          let y1 = edge_y1[i];
          if ((y0 <= scan_y && y1 > scan_y) || (y1 <= scan_y && y0 > scan_y)) {
            let t = (scan_y - y0)/(y1 - y0);
            xs_list.push(edge_x0[i] + t*(edge_x1[i] - edge_x0[i]));
          }
        }
        if (xs_list.length < 2) continue;
        xs_list.sort((a, b) => a - b);
        let row_off = sy*sub_w;
        for (let i = 0; i < xs_list.length - 1; i += 2) {
          let xs = Math.max(0, Math.ceil(xs_list[i] - 0.5));
          let xe = Math.min(sub_w - 1, Math.floor(xs_list[i + 1] - 0.5));
          for (let sx = xs; sx <= xe; sx++) {
            local_mask[row_off + sx] = 1;
            total_hits++;
          }
        }
      }
    }

    let pairs = [];
    let sub_per_cell = supersample*supersample;

    if (total_hits > 0) {
      for (let ly = 0; ly < local_h; ly++) {
        let gy = r0 + ly;
        for (let lx = 0; lx < local_w; lx++) {
          let cnt = 0;
          for (let fy = 0; fy < supersample; fy++) {
            let off = (ly*supersample + fy)*sub_w + lx*supersample;
            for (let fx = 0; fx < supersample; fx++)
              if (local_mask[off + fx] === 1) cnt++;
          }
          if (cnt > 0) {
            let gx = c0 + lx;
            pairs.push([gy*width + gx, cnt/sub_per_cell]);
          }
        }
      }
    } else if (polys.length > 0 && polys[0][0]) {
      let info = ScanlineRasteriser.computeRingAreaKm2(polys[0][0]);
      if (info.area_km2 >= 0.05) {
        let gy = Math.min(height - 1, Math.max(0, Math.floor(((90.0 - info.centroid_lat)/180.0)*height)));
        let gx = Math.min(width - 1, Math.max(0, Math.floor(((info.centroid_lon + 180.0)/360.0)*width)));
        pairs.push([gy*width + gx, Math.min(1.0, info.area_km2/cell_area_by_row[gy])]);
      }
    }

    //Return statement
    return pairs;
  }

  /**
   * Builds sparse reservoir overlay index from Global Dam Watch (GDW) v1.0, fixing string LAKE_CTRL ("Yes"/"Maybe")
   * and emitting a global reservoir union mask so HydroLAKES Lake_type==1 leaks are scrubbed cleanly.
   *
   * @param {Object} arg0_cfg
   * @param {Object} arg1_paths
   * @param {Float64Array} arg2_cell_area_by_row
   *
   * @returns {Promise<void>}
   */
  static async buildReservoirOverlay (arg0_cfg, arg1_paths, arg2_cell_area_by_row) {
    //Convert from parameters
    let cfg = arg0_cfg;
    let paths = arg1_paths;
    let cell_area_by_row = arg2_cell_area_by_row;

    //Declare local instance variables
    let force = process.argv.includes("--force");
    let out_path = path.join(paths.intermediate_dir, "reservoirs_sparse.json");
    let union_path = path.join(paths.intermediate_dir, "f_gdw_reservoirs_union.f32");

    if (!force && fs.existsSync(out_path) && fs.existsSync(union_path)) {
      console.log("[D_overlays] Cached: reservoirs_sparse.json & f_gdw_reservoirs_union.f32");
      return;
    }

    let shp_path = this.findReservoirShapefile(path.join(paths.raw_dir, "gdw_v10"));
    if (!shp_path) {
      console.error("[D_overlays] Could not locate GDW v1.0 reservoir shapefile.");
      return;
    }

    console.log(`[D_overlays] Streaming GDW v1.0 reservoirs: ${path.basename(shp_path)} ...`);
    let dbf_path = shp_path.replace(/\.shp$/i, ".dbf");
    let source = await shapefile.open(shp_path, fs.existsSync(dbf_path) ? dbf_path : undefined);
    let count = 0;
    let kept = 0;
    let min_area = cfg.hydrography.min_reservoir_area_km2;
    let fallback_year = cfg.hydrography.undated_reservoir_fallback_year;
    let skipped_natural_ctrl = 0;
    let t0 = Date.now();
    let union_mask = new Float32Array(4320*2160);

    let pixel_indices = [];
    let fracs = [];
    let impound_years = [];

    for (;;) {
      let result = await source.read();
      if (result.done) break;
      count++;

      let feat = result.value;
      let props = feat.properties ? feat.properties : {};

      //Exclude regulated natural lakes (GDW uses LAKE_CTRL = "Yes" | "Maybe", DAM_TYPE = "Lake Control Dam")
      //Allowlist true impounded mega-reservoirs mis-tagged LAKE_CTRL="Yes" in GDW (Manicouagan 1968, Notigi 1974)
      let res_name = String(props.RES_NAME || props.DAM_NAME || "");
      let force_keep_reservoir = /manicouagan|notigi/i.test(res_name);
      let res_type = String(props.RES_TYPE || props.DAM_TYPE || props.TYPE || "").toLowerCase();
      let ctrl_str = String(props.LAKE_CTRL || "").trim().toLowerCase();
      if (!force_keep_reservoir && (ctrl_str === "yes" || ctrl_str === "maybe" || ctrl_str === "1" || res_type.includes("lake control") || res_type.includes("natural"))) {
        skipped_natural_ctrl++;
        continue;
      }

      let area_km2 = Number(props.AREA_SKM || props.AREA_POLY || props.AREA_REP || 0);
      if (area_km2 > 0 && area_km2 < min_area) continue;
      if (area_km2 > 10000) {
        skipped_natural_ctrl++;
        continue;
      }

      let raw_year = Number(props.YEAR_DAM || props.YEAR || props.YEAR_SRC || -99);
      let impound_year = (raw_year >= 1800 && raw_year <= 2025) ? Math.round(raw_year) : fallback_year;

      //R3-D03: Inside large reservoirs (>=500 km2: Smallwood, Caniapiscau, Robert-Bourassa, Manicouagan, etc.),
      //strip sub-pixel micro-islet holes (<25 km2) while preserving all real islands >=25 km2 (incl. Rene-Levasseur 2020 km2!)
      let geom_use = feat.geometry;
      let res_ss = 10;
      if (area_km2 >= 500.0 && geom_use && geom_use.coordinates) {
        res_ss = 20;
        let filterPoly = (arg0_rings) => {
          if (!arg0_rings || arg0_rings.length <= 1) return arg0_rings;
          let kept_rings = [arg0_rings[0]];
          for (let r = 1; r < arg0_rings.length; r++) {
            if (ScanlineRasteriser.computeRingAreaKm2(arg0_rings[r]).area_km2 >= 25.0)
              kept_rings.push(arg0_rings[r]);
          }
          return kept_rings;
        };
        if (geom_use.type === "Polygon") {
          geom_use = { coordinates: filterPoly(geom_use.coordinates), type: "Polygon" };
        } else if (geom_use.type === "MultiPolygon") {
          geom_use = { coordinates: geom_use.coordinates.map(filterPoly), type: "MultiPolygon" };
        }
      }

      let pairs = this.rasteriseGeometrySparse(geom_use, cell_area_by_row, res_ss);
      if (pairs.length === 0) continue;

      kept++;
      for (let i = 0; i < pairs.length; i++) {
        let p_idx = pairs[i][0];
        let f = Math.round(pairs[i][1]*10000)/10000;
        pixel_indices.push(p_idx);
        fracs.push(f);
        impound_years.push(impound_year);
        if (f > union_mask[p_idx]) union_mask[p_idx] = f;
      }
    }

    fs.writeFileSync(out_path, JSON.stringify({
      count_reservoirs: kept,
      fracs: fracs,
      impound_years: impound_years,
      pixel_indices: pixel_indices,
    }));
    GeoGridUtils.writeFloat32Bin(union_path, union_mask);

    console.log(`[D_overlays] Saved reservoirs_sparse.json: ${kept} reservoirs (${skipped_natural_ctrl} regulated natural lakes skipped), ${pixel_indices.length} cell entries in ${((Date.now() - t0)/1000).toFixed(1)}s.`);
  }

  /**
   * Rasterises all 13 official RCE / Vos et al. (2021) Atlas of the Holocene Netherlands WFS GeoJSON epochs
   * into compact 5' Netherlands domain overrides + exact modern polder footprints.
   *
   * @param {Object} arg0_hcfg
   * @param {Object} arg1_paths
   */
  static buildVosNetherlands (arg0_hcfg, arg1_paths) {
    //Convert from parameters
    let hcfg = arg0_hcfg;
    let paths = arg1_paths;

    //Declare local instance variables
    let force = process.argv.includes("--force");
    let out_path = path.join(paths.intermediate_dir, "nl_vos_epochs.json");

    if (!force && fs.existsSync(out_path)) {
      console.log("[D_overlays] Cached: nl_vos_epochs.json");
      return;
    }

    let rce_dir = path.join(paths.raw_dir, "rce_netherlands");
    let vos = hcfg.rce_netherlands_vos;
    let water_set = new Set(vos.water_classes.map((s) => s.toLowerCase()));
    let early_water_set = new Set((vos.early_holocene_water_classes || vos.water_classes).map((s) => s.toLowerCase()));
    let early_max_year = (vos.early_holocene_max_year !== undefined) ? vos.early_holocene_max_year : -2500;
    let outside_set = new Set(vos.outside_classes.map((s) => s.toLowerCase()));

    //Fixed Netherlands bounding box on 4320x2160 grid (lat 50.6..53.8N, lon 3.1..7.4E)
    let r0 = Math.floor(((90.0 - 53.8)/180.0)*2160);
    let r1 = Math.ceil(((90.0 - 50.6)/180.0)*2160);
    let c0 = Math.floor(((3.1 + 180.0)/360.0)*4320);
    let c1 = Math.ceil(((7.4 + 180.0)/360.0)*4320);
    let local_h = r1 - r0;
    let local_w = c1 - c0;
    let s = 20;
    let sub_h = local_h*s;
    let sub_w = local_w*s;
    let scale_x = (4320*s)/360.0;
    let scale_y = (2160*s)/180.0;
    let origin_sx = c0*s;
    let origin_sy = r0*s;
    let sub_per_cell = s*s;

    let rasteriseMask = (arg0_features) => {
      let mask = new Uint8Array(sub_h*sub_w);
      let xs_list = [];

      for (let fi = 0; fi < arg0_features.length; fi++) {
        let geom = arg0_features[fi].geometry;
        if (!geom || !geom.coordinates) continue;
        let polys = (geom.type === "Polygon") ? [geom.coordinates] : (geom.type === "MultiPolygon" ? geom.coordinates : []);

        for (let p = 0; p < polys.length; p++) {
          let rings = polys[p];
          let edge_x0 = [], edge_y0 = [], edge_x1 = [], edge_y1 = [];
          let p_min_sy = sub_h;
          let p_max_sy = -1;

          for (let r = 0; r < rings.length; r++) {
            let ring = rings[r];
            if (!ring || ring.length < 3) continue;
            for (let k = 0; k < ring.length - 1; k++) {
              let sy0 = (90.0 - ring[k][1])*scale_y - origin_sy;
              let sy1 = (90.0 - ring[k + 1][1])*scale_y - origin_sy;
              if (sy0 === sy1) continue;
              let sx0 = (ring[k][0] + 180.0)*scale_x - origin_sx;
              let sx1 = (ring[k + 1][0] + 180.0)*scale_x - origin_sx;
              if (sy0 < p_min_sy) p_min_sy = sy0;
              if (sy1 < p_min_sy) p_min_sy = sy1;
              if (sy0 > p_max_sy) p_max_sy = sy0;
              if (sy1 > p_max_sy) p_max_sy = sy1;
              edge_x0.push(sx0);
              edge_y0.push(sy0);
              edge_x1.push(sx1);
              edge_y1.push(sy1);
            }
          }

          let n_edges = edge_x0.length;
          if (n_edges === 0) continue;

          let sy_start = Math.max(0, Math.floor(p_min_sy));
          let sy_end = Math.min(sub_h - 1, Math.ceil(p_max_sy));

          for (let sy = sy_start; sy <= sy_end; sy++) {
            let scan_y = sy + 0.5;
            xs_list.length = 0;
            for (let i = 0; i < n_edges; i++) {
              let y0 = edge_y0[i];
              let y1 = edge_y1[i];
              if ((y0 <= scan_y && y1 > scan_y) || (y1 <= scan_y && y0 > scan_y)) {
                let t = (scan_y - y0)/(y1 - y0);
                xs_list.push(edge_x0[i] + t*(edge_x1[i] - edge_x0[i]));
              }
            }
            if (xs_list.length < 2) continue;
            xs_list.sort((a, b) => a - b);
            let row_off = sy*sub_w;
            for (let i = 0; i < xs_list.length - 1; i += 2) {
              let xs = Math.max(0, Math.ceil(xs_list[i] - 0.5));
              let xe = Math.min(sub_w - 1, Math.floor(xs_list[i + 1] - 0.5));
              for (let sx = xs; sx <= xe; sx++)
                mask[row_off + sx] = 1;
            }
          }
        }
      }

      //Aggregate into local_h x local_w fractional Float32
      let frac = new Float32Array(local_h*local_w);
      for (let ly = 0; ly < local_h; ly++) {
        for (let lx = 0; lx < local_w; lx++) {
          let hits = 0;
          for (let fy = 0; fy < s; fy++) {
            let row_off = (ly*s + fy)*sub_w + lx*s;
            for (let fx = 0; fx < s; fx++)
              if (mask[row_off + fx] === 1) hits++;
          }
          frac[ly*local_w + lx] = Math.round((hits/sub_per_cell)*10000)/10000;
        }
      }
      return frac;
    };

    let epoch_grids = [];
    let domain_union = new Float32Array(local_h*local_w);
    let t0 = Date.now();

    for (let i = 0; i < vos.epochs.length; i++) {
      let ep = vos.epochs[i];
      let file_path = path.join(rce_dir, `vos_${ep.year}.geojson`);
      if (!fs.existsSync(file_path)) {
        console.error(`[D_overlays] Missing RCE Vos file: ${file_path}`);
        continue;
      }

      let gj = JSON.parse(fs.readFileSync(file_path, "utf8"));
      let land_feats = [];
      let nl_feats = [];
      let active_water_set = (ep.year <= early_max_year) ? early_water_set : water_set;

      for (let fi = 0; fi < gj.features.length; fi++) {
        let feat = gj.features[fi];
        let leg = String(feat.properties?.legenda || "").trim().toLowerCase();
        if (outside_set.has(leg)) continue;
        nl_feats.push(feat);
        if (!active_water_set.has(leg))
          land_feats.push(feat);
      }

      let f_land = rasteriseMask(land_feats);
      let f_dom = rasteriseMask(nl_feats);

      for (let k = 0; k < domain_union.length; k++)
        if (f_dom[k] > domain_union[k]) domain_union[k] = f_dom[k];

      epoch_grids.push({
        f_land: Array.from(f_land),
        year: ep.year,
      });
      console.log(`[D_overlays]   Vos NL ${ep.year}: ${land_feats.length}/${nl_feats.length} land features rasterised.`);
    }

    //Extract exact pixel indices inside the Netherlands domain (domain_union > 0.5)
    let out_obj = {
      c0: c0,
      c1: c1,
      domain_frac: Array.from(domain_union),
      epochs: epoch_grids,
      local_h: local_h,
      local_w: local_w,
      r0: r0,
      r1: r1,
    };

    fs.writeFileSync(out_path, JSON.stringify(out_obj));
    console.log(`[D_overlays] Built official RCE Vos et al. nl_vos_epochs.json (${epoch_grids.length} epochs, ${local_w}x${local_h} NL grid) in ${((Date.now() - t0)/1000).toFixed(1)}s.`);
  }

  /**
   * Precomputes exact 30" ETOPO topographic-contour masks for all seeded closed basins
   * (Caspian, Black Sea Neoeuxinian, Aral South/North, Lake Chad / Mega-Chad, Lake Urmia)
   * and Holocene alluvial prograding deltas (Mesopotamia/Persian Gulf, China Bohai/Yellow/Yangtze/PRD, Nile, Po)
   * across all 218 keyframes!
   *
   * @param {Object} arg0_cfg
   * @param {Object} arg1_hcfg
   * @param {Object} arg2_paths
   */
  static buildDemContourBasins (arg0_cfg, arg1_hcfg, arg2_paths) {
    //Convert from parameters
    let cfg = arg0_cfg;
    let hcfg = arg1_hcfg;
    let paths = arg2_paths;

    //Declare local instance variables
    let force = process.argv.includes("--force") || process.argv.includes("--force-basins");
    let out_path = path.join(paths.intermediate_dir, "dem_basins_by_year.json");

    if (!force && fs.existsSync(out_path)) {
      console.log("[D_overlays] Cached: dem_basins_by_year.json");
      return;
    }

    let engine = new BasinContourEngine(paths.intermediate_dir);
    let years = GeoGridUtils.getTimestepYears(cfg);
    let t0 = Date.now();
    let result = {
      deltas: {},
      seeded: {},
    };

    //1. Seeded closed / sill-controlled basins (deduplicate identical water levels for speed!)
    for (let i = 0; i < hcfg.seeded_basins.length; i++) {
      let b = hcfg.seeded_basins[i];
      let level_cache = new Map();
      let level_table = [];
      let year_to_level_id = {};

      for (let yi = 0; yi < years.length; yi++) {
        let y = years[yi];
        if (b.active_max_year !== undefined && y > b.active_max_year) continue;
        if (b.active_min_year !== undefined && y < b.active_min_year) continue;

        let z_raw = BasinContourEngine.evalCurve(b.level_curve_m, y);
        let z_q = Math.round(z_raw*2.0)/2.0; //Exact 0.5m resolution matching etopo_bed_halfm.i16!
        let key = z_q.toFixed(1);

        if (!level_cache.has(key)) {
          let pairs = engine.floodSeededContour(b.bbox, b.seed, z_q, {
            clean_stipples: Boolean(b.clean_stipples),
          });
          level_cache.set(key, level_table.length);
          level_table.push({ pairs: pairs, z_m: z_q });
        }
        year_to_level_id[y] = level_cache.get(key);
      }

      result.seeded[b.id] = {
        bbox: b.bbox,
        exclude_from_static_natlakes: Boolean(b.exclude_from_static_natlakes),
        is_marine_override: Boolean(b.is_marine_override),
        level_table: level_table,
        year_to_level_id: year_to_level_id,
      };
      console.log(`[D_overlays]   30" Seeded basin ${b.id}: ${level_table.length} unique elevation contours computed.`);
    }

    //2. Holocene prograding alluvial deltas (ocean-connected back-stripping on 30" DEM)
    //   Computes delta marine transgression mask for each unique backstrip depth (0.5m steps),
    //   subtracting the modern 0.0m ocean footprint so ONLY newly flooded historical delta plains are tagged!
    for (let i = 0; i < hcfg.delta_progradation_basins.length; i++) {
      let d = hcfg.delta_progradation_basins[i];
      let base_pairs = engine.floodSeededContour(d.bbox, d.sea_seed, 0.0, { seed_edge_sea: true });
      let base_map = new Map();
      for (let k = 0; k < base_pairs.length; k++)
        base_map.set(base_pairs[k][0], base_pairs[k][1]);

      let level_cache = new Map();
      let level_table = [];
      let year_to_level_id = {};

      for (let yi = 0; yi < years.length; yi++) {
        let y = years[yi];
        if (d.min_year !== undefined && y < d.min_year) continue;
        let bs_raw = BasinContourEngine.evalCurve(d.backstrip_curve_m, y);
        if (bs_raw < 0.05) continue;

        let bs_q = Math.round(bs_raw*10.0)/10.0;
        let key = bs_q.toFixed(1);

        if (!level_cache.has(key)) {
          let flooded_pairs = engine.floodSeededContour(d.bbox, d.sea_seed, bs_q, {
            clean_stipples: Boolean(d.clean_stipples),
            seed_edge_sea: true,
            smooth_substep: true,
          });
          //Compute proportional newly-flooded fraction (f_new) & total flooded fraction (f_tot) on bbox grid
          let r0 = Math.max(0, Math.floor(((90.0 - d.bbox.lat_max)/180.0)*2160));
          let r1 = Math.min(2160, Math.ceil(((90.0 - d.bbox.lat_min)/180.0)*2160));
          let c0 = Math.max(0, Math.floor(((d.bbox.lon_min + 180.0)/360.0)*4320));
          let c1 = Math.min(4320, Math.ceil(((d.bbox.lon_max + 180.0)/360.0)*4320));
          let bh = r1 - r0;
          let bw = c1 - c0;
          let f_new_grid = new Float32Array(bh*bw);
          let f_tot_grid = new Float32Array(bh*bw);

          for (let k = 0; k < flooded_pairs.length; k++) {
            let p_idx = flooded_pairs[k][0];
            let f_now = flooded_pairs[k][1];
            let gy = (p_idx/4320) | 0;
            let gx = p_idx - gy*4320;
            if (gy < r0 || gy >= r1 || gx < c0 || gx >= c1) continue;
            let li = (gy - r0)*bw + (gx - c0);
            let f_base = base_map.get(p_idx) || 0.0;
            f_tot_grid[li] = f_now;
            f_new_grid[li] = Math.max(0.0, f_now - f_base);
          }

          let added_water = [];
          for (let ly = 0; ly < bh; ly++) {
            let gy = r0 + ly;
            let y_a = Math.max(0, ly - 3);
            let y_b = Math.min(bh - 1, ly + 3);
            for (let lx = 0; lx < bw; lx++) {
              let li = ly*bw + lx;
              let f_new = f_new_grid[li];
              let f_tot = f_tot_grid[li];
              let max_nb_new = 0.0;

              if (f_tot >= 0.85) {
                let x_a = Math.max(0, lx - 3);
                let x_b = Math.min(bw - 1, lx + 3);
                for (let ny = y_a; ny <= y_b; ny++) {
                  let n_off = ny*bw;
                  for (let nx = x_a; nx <= x_b; nx++) {
                    let v = f_new_grid[n_off + nx];
                    if (v > max_nb_new) max_nb_new = v;
                  }
                }
              }

              let loss_zero = (f_tot >= 0.85 && max_nb_new > 0.01) ? Math.min(1.0, max_nb_new/0.25) : 0.0;
              if (f_new > 0.01 || loss_zero > 0.01) {
                added_water.push([gy*4320 + c0 + lx, Math.round(f_new*10000)/10000, Math.round(loss_zero*10000)/10000]);
              }
            }
          }
          level_cache.set(key, level_table.length);
          level_table.push({ backstrip_m: bs_q, pairs: added_water });
        }
        year_to_level_id[y] = level_cache.get(key);
      }

      result.deltas[d.id] = {
        bbox: d.bbox,
        level_table: level_table,
        mode: "proportional",
        year_to_level_id: year_to_level_id,
      };
      console.log(`[D_overlays]   30" Prograding delta ${d.id}: ${level_table.length} unique back-strip steps computed.`);
    }

    //3. Smooth 30" SDF-migrating paleoshoreline deltas (China Bohai/Jiangsu/Yangtze/PRD, Rann of Kutch,
    //   Chao Phraya/Bangkok Sea, Mekong, Ganges-Brahmaputra, Mississippi Birdfoot)
    let sdf_specs = hcfg.sdf_delta_progradation_basins || [];
    let sdf_gj_path = path.join(paths.curated_dir, "delta_paleoshorelines.geojson");
    if (sdf_specs.length > 0 && fs.existsSync(sdf_gj_path)) {
      let sdf_gj = JSON.parse(fs.readFileSync(sdf_gj_path, "utf8"));
      let by_delta = {};
      for (let i = 0; i < sdf_gj.features.length; i++) {
        let feat = sdf_gj.features[i];
        let did = feat.properties?.delta_id;
        if (!did) continue;
        if (!by_delta[did]) by_delta[did] = [];
        by_delta[did].push({
          geometry: feat.geometry,
          year: Number(feat.properties.year),
        });
      }

      for (let i = 0; i < sdf_specs.length; i++) {
        let sd = sdf_specs[i];
        let kfs = by_delta[sd.id] ? by_delta[sd.id].slice() : [];
        if (kfs.length === 0) {
          console.error(`[D_overlays] No paleoshoreline keyframes found for SDF delta: ${sd.id}`);
          continue;
        }
        if (sd.modern_anchor_year !== undefined)
          kfs.push({ use_modern_sea: true, year: sd.modern_anchor_year });

        result.deltas[sd.id] = engine.buildSdfDeltaSeries({
          bbox: sd.bbox,
          keyframes: kfs,
          max_alluvial_elev_m: sd.max_alluvial_elev_m,
          max_year: sd.max_year,
          min_year: sd.min_year,
          preserve_hills: sd.preserve_hills,
        }, years);
        console.log(`[D_overlays]   30" SDF migrating delta ${sd.id}: ${result.deltas[sd.id].level_table.length} unique shoreline steps computed.`);
      }
    }

    fs.writeFileSync(out_path, JSON.stringify(result));
    console.log(`[D_overlays] Built dem_basins_by_year.json in ${((Date.now() - t0)/1000).toFixed(1)}s.`);
  }

  /**
   * Builds remaining non-DEM curated coastal reclamation (excl. NL which is handled by RCE Vos + exact polders)
   * and remaining curated proglacial/historic paleolakes (Lake Agassiz-Ojibway, Texcoco, Tulare, Lop Nur, Salton).
   *
   * @param {Object} arg0_paths
   * @param {Float64Array} arg1_cell_area_by_row
   */
  static buildCuratedOverlays (arg0_paths, arg1_cell_area_by_row) {
    //Convert from parameters
    let paths = arg0_paths;
    let cell_area_by_row = arg1_cell_area_by_row;

    //1. Coastal land reclamation (non-Vos global regions + exact post-1927 NL modern IJsselmeer/Maasvlakte dates)
    let rec_registry = GeoGridUtils.loadConfig("reclamation_registry.json5");
    let rec_geojson = JSON.parse(fs.readFileSync(path.join(paths.curated_dir, "reclamation_polygons.geojson"), "utf8"));
    let rec_feat_map = {};
    for (let i = 0; i < rec_geojson.features.length; i++)
      rec_feat_map[rec_geojson.features[i].properties.id] = rec_geojson.features[i].geometry;

    //Skip pre-1927 Dutch events here because official RCE Vos et al. (2021) handles 9000 BC..1850 AD!
    let skip_nl_pre_1927 = new Set(["nl_haarlemmermeer", "nl_north_holland_polders", "nl_zuiderzee_formation"]);
    let rec_entries = [];
    for (let i = 0; i < rec_registry.events.length; i++) {
      let ev = rec_registry.events[i];
      if (skip_nl_pre_1927.has(ev.id)) continue;
      let geom = rec_feat_map[ev.id];
      if (!geom) continue;
      let pairs = this.rasteriseGeometrySparse(geom, cell_area_by_row, 20);
      rec_entries.push({
        action: ev.action,
        end_year: ev.end_year,
        id: ev.id,
        name: ev.name,
        pairs: pairs,
        start_year: ev.start_year,
      });
    }
    fs.writeFileSync(path.join(paths.intermediate_dir, "reclamation_sparse.json"), JSON.stringify(rec_entries));
    console.log(`[D_overlays] Built reclamation_sparse.json (${rec_entries.length} events).`);

    //2. Remaining curated paleolakes NOT superseded by 30" DEM contour basins (Agassiz, Texcoco, Tulare, Lop Nur, Salton)
    let superseded_lakes = new Set(["paleolake_mega_chad", "aral_sea_1960_extent", "aral_sea_west_deep", "lake_chad_1960", "lake_urmia_pre_2000"]);
    let lake_registry = GeoGridUtils.loadConfig("paleolake_registry.json5");
    let lake_geojson = JSON.parse(fs.readFileSync(path.join(paths.curated_dir, "paleolake_polygons.geojson"), "utf8"));
    let min_bed_path = path.join(paths.intermediate_dir, "etopo_min_bed_5m.f32");
    let min_bed_5m = fs.existsSync(min_bed_path) ? GeoGridUtils.readFloat32Bin(min_bed_path) : null;
    let lake_feat_map = {};
    for (let i = 0; i < lake_geojson.features.length; i++)
      lake_feat_map[lake_geojson.features[i].properties.id] = lake_geojson.features[i].geometry;

    let lake_entries = [];
    for (let i = 0; i < lake_registry.lakes.length; i++) {
      let lk = lake_registry.lakes[i];
      if (superseded_lakes.has(lk.id)) continue;
      let geom = lake_feat_map[lk.id];
      if (!geom) continue;
      let raw_pairs = this.rasteriseGeometrySparse(geom, cell_area_by_row, 20);
      let enriched = raw_pairs.map((pair) => ({
        elev: min_bed_5m ? min_bed_5m[pair[0]] : 0.0,
        frac: pair[1],
        idx: pair[0],
      }));
      enriched.sort((a, b) => a.elev - b.elev);
      let n_e = enriched.length;
      let ranked_pairs = enriched.map((item, rank) => [
        item.idx,
        Math.round(item.frac*10000)/10000,
        (n_e > 1) ? Math.round((rank/(n_e - 1))*10000)/10000 : 0.0,
      ]);

      lake_entries.push({
        action: lk.action,
        curve: lk.curve,
        id: lk.id,
        name: lk.name,
        pairs: ranked_pairs,
      });
    }
    fs.writeFileSync(path.join(paths.intermediate_dir, "paleolake_sparse.json"), JSON.stringify(lake_entries));
    console.log(`[D_overlays] Built hypsometric paleolake_sparse.json (${lake_entries.length} non-DEM paleolakes).`);
  }

  /**
   * Extracts dynamic lagoon footprints (e.g. Garabogazköl) from the HydroLAKES 2025 anchor and ranks cells inward from the outer shore using Chamfer distance so desiccation shrinks organically toward the lagoon core.
   *
   * @param {Object} arg0_hcfg
   * @param {Object} arg1_paths
   */
  static buildLagoonBasins (arg0_hcfg, arg1_paths) {
    //Convert from parameters
    let hcfg = arg0_hcfg;
    let paths = arg1_paths;

    //Declare local instance variables
    let force = process.argv.includes("--force") || process.argv.includes("--force-lagoons");
    let out_path = path.join(paths.intermediate_dir, "lagoon_basins_sparse.json");
    let specs = hcfg.lagoon_basins || [];

    if (!force && fs.existsSync(out_path)) {
      console.log("[D_overlays] Cached: lagoon_basins_sparse.json");
      return;
    }

    let f_natlake = GeoGridUtils.readFloat32Bin(path.join(paths.intermediate_dir, "f_natlake_2025.f32"));
    let f_osm = GeoGridUtils.readFloat32Bin(path.join(paths.intermediate_dir, "f_land_osm_2025.f32"));
    let out = [];

    for (let i = 0; i < specs.length; i++) {
      let sp = specs[i];
      let b = sp.bbox;
      let r0 = Math.max(0, Math.floor(((90.0 - b.lat_max)/180.0)*2160));
      let r1 = Math.min(2160, Math.ceil(((90.0 - b.lat_min)/180.0)*2160));
      let c0 = Math.max(0, Math.floor(((b.lon_min + 180.0)/360.0)*4320));
      let c1 = Math.min(4320, Math.ceil(((b.lon_max + 180.0)/360.0)*4320));
      let w = c1 - c0;
      let h = r1 - r0;
      let rim_mask = new Uint8Array(w*h);
      let in_lagoon = new Uint8Array(w*h);

      for (let ly = 0; ly < h; ly++) {
        let gy = r0 + ly;
        let lat = 90.0 - (gy + 0.5)/12.0;
        let g_off = gy*4320 + c0;
        let l_off = ly*w;
        for (let lx = 0; lx < w; lx++) {
          let lon = -180.0 + (c0 + lx + 0.5)/12.0;
          let g_idx = g_off + lx;
          //Separate Garabogazkol (OSM inland landmass + HydroLAKES water) from open Caspian Sea on the west spit
          let west_spit_cut = (lat < 41.25 && lon < 52.90);
          let member = (!west_spit_cut && f_osm[g_idx] >= 0.45 && f_natlake[g_idx] > 0.05);
          if (member) {
            in_lagoon[l_off + lx] = 1;
          } else {
            rim_mask[l_off + lx] = 1;
          }
        }
      }

      let dist = BasinContourEngine.computeChamferDistance(rim_mask, w, h);
      let max_d = 1.0;
      for (let k = 0; k < dist.length; k++)
        if (in_lagoon[k] === 1 && dist[k] > max_d) max_d = dist[k];

      let pairs = [];
      for (let ly = 0; ly < h; ly++) {
        let gy = r0 + ly;
        let g_off = gy*4320 + c0;
        let l_off = ly*w;
        for (let lx = 0; lx < w; lx++) {
          if (in_lagoon[l_off + lx] === 0) continue;
          let frac = f_natlake[g_off + lx];
          //rank_norm: 0 = deepest central lagoon, 1 = natural curved outer shoreline
          let rank_norm = 1.0 - Math.min(1.0, dist[l_off + lx]/max_d);
          pairs.push([
            g_off + lx,
            Math.round(frac*10000)/10000,
            Math.round(rank_norm*10000)/10000,
          ]);
        }
      }

      out.push({
        id: sp.id,
        level_curve: sp.level_curve,
        name: sp.name,
        pairs: pairs,
      });
      console.log(`[D_overlays]   Lagoon basin ${sp.id}: ${pairs.length} cells ranked from HydroLAKES.`);
    }

    fs.writeFileSync(out_path, JSON.stringify(out));
    console.log(`[D_overlays] Built lagoon_basins_sparse.json (${out.length} basins).`);
  }


  /**
   * Runs the Arctic Siberian hydrology simulation (channels, floodplains, fill-and-spill lakes) per relative sea-level tier and samples regional GLAC-1D + calibrated eustatic sea-level histories into arctic_hydrology.json.
   *
   * @param {Object} arg0_cfg
   * @param {Object} arg1_hcfg
   * @param {Object} arg2_paths
   *
   * @returns {Promise<void>}
   */
  static async buildPalaeoRivers (arg0_cfg, arg1_hcfg, arg2_paths) {
    //Convert from parameters
    let cfg = arg0_cfg;
    let hcfg = arg1_hcfg;
    let paths = arg2_paths;

    //Declare local instance variables
    let force = process.argv.includes("--force") || process.argv.includes("--force-palaeo-rivers");
    let out_path = path.join(paths.intermediate_dir, "arctic_hydrology.json");
    let spec = hcfg.arctic_hydrology;

    if (!force && fs.existsSync(out_path)) {
      console.log("[D_overlays] Cached: arctic_hydrology.json");
      return;
    }

    let engine = new BasinContourEngine(paths.intermediate_dir);
    let t0 = Date.now();
    let res = ArcticHydrologyModel.build(spec, engine);
    let sectors = res.sectors || [];

    //Sample GLAC-1D + calibrated eustatic RSL trajectory at each Siberian Arctic sector
    let glac_path = path.join(paths.raw_dir, "TOPicemsk.GLACD26kN9894GE90227A6005GGrBgic.nc");
    let glac_nc = await NetCDF4Writer.openRead(glac_path);
    let keys = glac_nc.keys();
    let time_key = keys.find((k) => /^t(ime)?$/i.test(k)) || "T122KP11";
    let lat_key = keys.find((k) => /^y|^lat/i.test(k));
    let lon_key = keys.find((k) => /^x|^lon/i.test(k));
    let hdc_key = keys.find((k) => k === "HDC") || keys.find((k) => /^hdc$/i.test(k));
    let t_vals = Array.from(glac_nc.getValue(time_key));
    let lat_vals = Array.from(glac_nc.getValue(lat_key));
    let lon_vals = Array.from(glac_nc.getValue(lon_key));
    let shape = glac_nc.getShape(hdc_key);
    let ntime = shape[0];
    let nlat = shape[1];
    let nlon = shape[2];
    let s2n = lat_vals[0] < lat_vals[lat_vals.length - 1];
    let lon360 = lon_vals[lon_vals.length - 1] > 180.5;
    let max_abs_t = Math.max(...t_vals.map((v) => Math.abs(v)));
    let cal_years = t_vals.map((v) => (max_abs_t <= 40.0) ? (cfg.paleo.bp_reference_year - Math.abs(v)*1000.0) : (cfg.paleo.bp_reference_year - Math.abs(v)));

    let present_idx = 0;
    let best_d = Infinity;
    for (let i = 0; i < ntime; i++) {
      let d = Math.abs(cal_years[i] - cfg.paleo.bp_reference_year);
      if (d < best_d) {
        best_d = d;
        present_idx = i;
      }
    }
    let hdc_present = Float32Array.from(glac_nc.getSlice(hdc_key, [[present_idx, present_idx + 1], [0, nlat], [0, nlon]]));

    let rsl_by_year = {};
    let asc = cal_years[cal_years.length - 1] > cal_years[0];
    let years = GeoGridUtils.getTimestepYears(cfg).filter((y) => y < cfg.paleo.gia_cutoff_year);

    for (let yi = 0; yi < years.length; yi++) {
      let yr = years[yi];
      let i0 = 0;
      let i1 = 0;

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
      let hdc0 = glac_nc.getSlice(hdc_key, [[i0, i0 + 1], [0, nlat], [0, nlon]]);
      let hdc1 = (i0 === i1) ? hdc0 : glac_nc.getSlice(hdc_key, [[i1, i1 + 1], [0, nlat], [0, nlon]]);
      let eust_cal = C_BuildPaleoGiaIce.getCalibratedEustaticOffsetM(yr);
      let row_rsls = [];

      for (let si = 0; si < sectors.length; si++) {
        let sec = sectors[si];
        let v0 = C_BuildPaleoGiaIce.sampleBilinear(hdc0, nlat, nlon, s2n, lon360, sec.sample_lat, sec.sample_lon);
        let v1 = (i0 === i1) ? v0 : C_BuildPaleoGiaIce.sampleBilinear(hdc1, nlat, nlon, s2n, lon360, sec.sample_lat, sec.sample_lon);
        let vp = C_BuildPaleoGiaIce.sampleBilinear(hdc_present, nlat, nlon, s2n, lon360, sec.sample_lat, sec.sample_lon);
        let rsl_m = -((w0*v0 + w1*v1 - vp) + eust_cal);
        let rsl_round = Math.round(rsl_m*10.0)/10.0;
        row_rsls.push(rsl_round);
      }

      rsl_by_year[yr] = row_rsls;
    }

    glac_nc.close();
    res.rsl_by_year = rsl_by_year;
    res.persistent_cells = this.listPersistentEstuaryCells(hcfg.holocene_infill, paths);
    fs.writeFileSync(out_path, JSON.stringify(res));
    console.log(`[D_overlays] Built arctic_hydrology.json (${res.cells.length} water cells, ${Object.keys(rsl_by_year).length} paleo years) in ${((Date.now() - t0)/1000).toFixed(1)}s.`);
  }


  /**
   * Lists modern estuary water cells (connected to each persistent_water domain seed through cells with OSM land fraction < 0.98) that remain water in every year as the palaeo river/lake.
   *
   * @param {Object} arg0_infill_spec - holocene_infill section of hydro_geomorph_basins.json5.
   * @param {Object} arg1_paths
   *
   * @returns {Array<number>} 5' pixel indices.
   */
  static listPersistentEstuaryCells (arg0_infill_spec, arg1_paths) {
    //Convert from parameters
    let spec = arg0_infill_spec ? arg0_infill_spec : {};
    let paths = arg1_paths;

    //Declare local instance variables
    let all_domains = (spec.domains || []).filter((dm) => dm.persistent_water);
    let f_osm = GeoGridUtils.readFloat32Bin(path.join(paths.intermediate_dir, "f_land_osm_2025.f32"));
    let out = new Set();

    for (let di = 0; di < all_domains.length; di++) {
      let dm = all_domains[di];
      let r0 = Math.floor((90.0 - dm.bbox.lat_max)*12);
      let r1 = Math.ceil((90.0 - dm.bbox.lat_min)*12);
      let c0 = Math.floor((dm.bbox.lon_min + 180.0)*12);
      let c1 = Math.ceil((dm.bbox.lon_max + 180.0)*12);
      let seed = Math.floor((90.0 - dm.seed.lat)*12)*4320 + Math.floor((dm.seed.lon + 180.0)*12);
      let queue = [seed];
      let seen = new Set([seed]);

      for (let t = 0; t < queue.length; t++) {
        let cur = queue[t];
        let cy = (cur/4320) | 0;
        let cx = cur - cy*4320;
        out.add(cur);
        let all_nb = [[cy - 1, cx], [cy + 1, cx], [cy, cx - 1], [cy, cx + 1]];
        for (let k = 0; k < 4; k++) {
          let ny = all_nb[k][0];
          let nx = all_nb[k][1];
          if (ny < r0 || ny >= r1 || nx < c0 || nx >= c1) continue;
          let ni = ny*4320 + nx;
          if (seen.has(ni) || f_osm[ni] >= 0.98) continue;
          seen.add(ni);
          queue.push(ni);
        }
      }
    }

    //Return statement
    return Array.from(out);
  }


  /**
   * Floods ice-dammed proglacial lakes (Lake Agassiz-Ojibway) on the 30" DEM against a tilted water plane,
   * bounded by each year's paleo ice front, writing proglacial_lakes_by_year.json.
   *
   * @param {Object} arg0_cfg
   * @param {Object} arg1_hcfg
   * @param {Object} arg2_paths
   */
  static buildProglacialLakes (arg0_cfg, arg1_hcfg, arg2_paths) {
    //Convert from parameters
    let cfg = arg0_cfg;
    let hcfg = arg1_hcfg;
    let paths = arg2_paths;

    //Declare local instance variables
    let force = process.argv.includes("--force") || process.argv.includes("--force-proglacial");
    let out_path = path.join(paths.intermediate_dir, "proglacial_lakes_by_year.json");
    let specs = hcfg.proglacial_lakes || [];

    if (!force && fs.existsSync(out_path)) {
      console.log("[D_overlays] Cached: proglacial_lakes_by_year.json");
      return;
    }

    let engine = new BasinContourEngine(paths.intermediate_dir);
    let paleo_dir = path.join(paths.intermediate_dir, "paleo_slices");
    let result = {};
    let t0 = Date.now();
    let years = GeoGridUtils.getTimestepYears(cfg);

    for (let i = 0; i < specs.length; i++) {
      let sp = specs[i];
      let by_year = {};
      let exclude_idx = this.listCellsInRing(sp.exclude_ring, sp.bbox);

      for (let yi = 0; yi < years.length; yi++) {
        let y = years[yi];
        if (y < sp.active_min_year || y > sp.active_max_year) continue;
        let ice_path = path.join(paleo_dir, `f_ice_${y}.f32`);
        if (!fs.existsSync(ice_path)) {
          console.error(`[D_overlays] Missing paleo ice slice: ${ice_path}`);
          continue;
        }

        let level = BasinContourEngine.evalCurve(sp.level_ref_curve_m, y);
        let gradient = BasinContourEngine.evalCurve(sp.gradient_curve_m_per_km, y);
        let block = GeoGridUtils.readFloat32Bin(ice_path);
        for (let k = 0; k < exclude_idx.length; k++)
          block[exclude_idx[k]] = 1.0;
        by_year[y] = engine.floodTiltedProglacial(sp, level, gradient, block);
      }

      result[sp.id] = { by_year: by_year };
      console.log(`[D_overlays]   30" proglacial lake ${sp.id}: ${Object.keys(by_year).length} keyframes flooded.`);
    }

    fs.writeFileSync(out_path, JSON.stringify(result));
    console.log(`[D_overlays] Built proglacial_lakes_by_year.json in ${((Date.now() - t0)/1000).toFixed(1)}s.`);
  }


  /**
   * Lists 5' pixel indices inside bbox whose centres fall within a [lon, lat] ring (even-odd rule).
   *
   * @param {Array<Array<number>>} [arg0_ring]
   * @param {Object} arg1_bbox - { lat_min, lat_max, lon_min, lon_max }
   *
   * @returns {Array<number>}
   */
  static listCellsInRing (arg0_ring, arg1_bbox) {
    //Convert from parameters
    let ring = arg0_ring ? arg0_ring : [];
    let bbox = arg1_bbox;

    //Declare local instance variables
    let out = [];
    let r0 = Math.max(0, Math.floor((90.0 - bbox.lat_max)*12));
    let r1 = Math.min(2160, Math.ceil((90.0 - bbox.lat_min)*12));
    let c0 = Math.max(0, Math.floor((bbox.lon_min + 180.0)*12));
    let c1 = Math.min(4320, Math.ceil((bbox.lon_max + 180.0)*12));

    //Guard clause
    if (ring.length < 4) return out;

    for (let r = r0; r < r1; r++) {
      let lat = 90.0 - (r + 0.5)/12;
      for (let c = c0; c < c1; c++) {
        let lon = -180.0 + (c + 0.5)/12;
        let inside = false;
        for (let k = 0, j = ring.length - 1; k < ring.length; j = k++) {
          let xi = ring[k][0], yi = ring[k][1], xj = ring[j][0], yj = ring[j][1];
          if (((yi > lat) !== (yj > lat)) && (lon < ((xj - xi)*(lat - yi))/(yj - yi) + xi))
            inside = !inside;
        }
        if (inside) out.push(r*4320 + c);
      }
    }

    //Return statement
    return out;
  }


  /**
   * Executes Step D.
   *
   * @returns {Promise<void>}
   */
  static async run () {
    //Declare local instance variables
    let cfg = GeoGridUtils.loadConfig("model_config.json5");
    let hcfg = GeoGridUtils.loadConfig("hydro_geomorph_basins.json5");
    let cell_area_by_row = GeoGridUtils.computeCellAreaByRow(cfg.grid);
    let paths = GeoGridUtils.resolvePaths("landarea");

    if (!fs.existsSync(paths.intermediate_dir))
      fs.mkdirSync(paths.intermediate_dir, { recursive: true });

    await this.buildReservoirOverlay(cfg, paths, cell_area_by_row);
    this.buildVosNetherlands(hcfg, paths);
    this.buildDemContourBasins(cfg, hcfg, paths);
    this.buildLagoonBasins(hcfg, paths);
    await this.buildPalaeoRivers(cfg, hcfg, paths);
    this.buildProglacialLakes(cfg, hcfg, paths);
    this.buildCuratedOverlays(paths, cell_area_by_row);
  }
}

if (require.main === module) {
  D_BuildDynamicOverlays.run().catch((err) => {
    console.error("[D_overlays] Fatal error:", err);
    process.exit(1);
  });
}

module.exports = D_BuildDynamicOverlays;
