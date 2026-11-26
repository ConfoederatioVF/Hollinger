"use strict";

let fs = require("node:fs");
let path = require("node:path");

let shapefile = require("shapefile");

let GeoGridUtils = require("../../core/GeoGridUtils");
let NetCDF4Writer = require("../../core/NetCDF4Writer");
let ScanlineRasteriser = require("../../core/ScanlineRasteriser");

/**
 * Step B: Constructs the inclusive 2025 AD modern anchor layers on the 4320x2160 WGS84 grid:
 *  - f_land_osm_2025.f32      (Exact sub-pixel OSM coastlines + preserved micro-islands/atolls)
 *  - f_natlake_2025.f32       (HydroLAKES v1.0 natural & regulated natural lakes, excl. reservoirs)
 *  - f_ice_2025.f32           (ETOPO 2022 polar ice sheets + permanent floating shelves + NE mountain glaciers)
 *  - f_etopo_land_2025.f32    (ETOPO 2022 30" baseline terrestrial fraction for ΔRSL differencing)
 *  - etopo_30s_bed_halfm.i16  (North-to-south 30" bedrock DEM in 0.5m units for paleo GIA shoreline modeling)
 */
class B_BuildModernAnchor {
  /**
   * Recursively locates the first .shp file inside a directory.
   *
   * @param {string} arg0_dir
   * @param {string} [arg1_substring]
   *
   * @returns {string|null}
   */
  static findShapefile (arg0_dir, arg1_substring) {
    //Convert from parameters
    let dir = arg0_dir;
    let substring = arg1_substring ? arg1_substring.toLowerCase() : "";

    //Guard clause
    if (!fs.existsSync(dir)) return null;

    //Declare local instance variables
    let entries = fs.readdirSync(dir, { withFileTypes: true });

    for (let i = 0; i < entries.length; i++) {
      let full = path.join(dir, entries[i].name);
      if (entries[i].isDirectory()) {
        let found = this.findShapefile(full, substring);
        if (found) return found;
      } else if (entries[i].name.toLowerCase().endsWith(".shp")) {
        if (!substring || entries[i].name.toLowerCase().includes(substring))
          return full;
      }
    }

    //Return statement
    return null;
  }

  /**
   * Tests whether a coordinate lies on Greenland proper (excluding Iceland and the Canadian Arctic
   * Archipelago west of Nares Strait / Baffin Bay / Davis Strait).
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
    //Exclude Iceland & Denmark Strait east of -26W south of 67.2N
    if (lat < 67.2 && lon > -26.0) return false;

    //Western maritime boundary separating Greenland from Ellesmere, Devon & Baffin Islands
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
   * Streams a shapefile through ScanlineRasteriser with optional feature filter callback.
   *
   * @param {string} arg0_shp_path
   * @param {Float64Array} arg1_cell_area_by_row
   * @param {Object} [arg2_options]
   *  @param {Function} [arg2_options.filter_fn]
   *  @param {string} [arg2_options.label]
   *  @param {number} [arg2_options.min_micro_area_km2=0.001]
   *  @param {number} [arg2_options.supersample=20]
   *
   * @returns {Promise<Float32Array>}
   */
  static async rasteriseShapefile (arg0_shp_path, arg1_cell_area_by_row, arg2_options) {
    //Convert from parameters
    let shp_path = arg0_shp_path;
    let cell_area_by_row = arg1_cell_area_by_row;
    let options = arg2_options ? arg2_options : {};

    //Initialise options
    if (!options.label) options.label = path.basename(shp_path);
    if (options.min_micro_area_km2 === undefined) options.min_micro_area_km2 = 0.001;
    if (options.supersample === undefined) options.supersample = 20;

    //Declare local instance variables
    let count = 0;
    let dbf_path = shp_path.replace(/\.shp$/i, ".dbf");
    let kept = 0;
    let rasteriser = new ScanlineRasteriser({
      min_micro_area_km2: options.min_micro_area_km2,
      supersample: options.supersample,
    });
    let t0 = Date.now();
    let source = await shapefile.open(shp_path, fs.existsSync(dbf_path) ? dbf_path : undefined);

    console.log(`[B_anchor] Streaming ${options.label} (S=${options.supersample}) ...`);

    for (;;) {
      let result = await source.read();
      if (result.done) break;
      count++;

      let feature = result.value;
      if (options.filter_fn && !options.filter_fn(feature)) continue;

      kept++;
      let min_hole_km2 = options.hole_filter_fn ? options.hole_filter_fn(feature) : 0.0;
      rasteriser.rasteriseGeometry(feature.geometry, { min_hole_area_km2: min_hole_km2 });

      if (count % 100000 === 0)
        console.log(`[B_anchor]   ${options.label}: ${count} features read (${kept} rasterised)`);
    }

    let frac = rasteriser.aggregateToFraction(cell_area_by_row);
    console.log(`[B_anchor] Finished ${options.label}: ${kept}/${count} features in ${((Date.now() - t0)/1000).toFixed(1)}s.`);

    //Return statement
    return frac;
  }

  /**
   * Processes ETOPO 2022 30" (or 60" fallback) bedrock and ice-surface NetCDF files to build:
   *  - f_ice_etopo_2025 (grounded polar ice sheets + permanent floating Antarctic ice shelves)
   *  - f_etopo_land_2025 (baseline elevation >= 0 fraction for paleo ΔRSL differencing)
   *  - etopo_bed_halfm.i16 (compact north-to-south bedrock DEM in 0.5m units)
   *
   * @param {Object} arg0_model_config
   * @param {Object} arg1_paths
   *
   * @returns {Promise<{ f_etopo_land: Float32Array, f_ice_etopo: Float32Array, meta: Object }>}
   */
  static async processEtopo2022 (arg0_model_config, arg1_paths) {
    //Convert from parameters
    let cfg = arg0_model_config;
    let paths = arg1_paths;

    //Declare local instance variables
    let bed_path = path.join(paths.raw_dir, "ETOPO_2022_v1_30s_N90W180_bed.nc");
    let surf_path = path.join(paths.raw_dir, "ETOPO_2022_v1_30s_N90W180_surface.nc");
    let ice_thresh = cfg.encoding.etopo_ice_thickness_threshold_m;
    let include_shelves = Boolean(cfg.encoding.include_permanent_ice_shelves);
    let out_h = cfg.grid.height;
    let out_w = cfg.grid.width;

    console.log("[B_anchor] Opening ETOPO 2022 bedrock & surface NetCDF-4 datasets ...");
    let bed_nc = await NetCDF4Writer.openRead(bed_path);
    let surf_nc = await NetCDF4Writer.openRead(surf_path);

    let bed_var = bed_nc.keys().includes("z") ? "z" : bed_nc.keys().find((k) => k !== "lat" && k !== "lon" && k !== "crs");
    let surf_var = surf_nc.keys().includes("z") ? "z" : surf_nc.keys().find((k) => k !== "lat" && k !== "lon" && k !== "crs");
    let lat_vals = bed_nc.getValue("lat");
    let shape = bed_nc.getShape(bed_var);
    let src_h = shape[0];
    let src_w = shape[1];
    let south_to_north = lat_vals[0] < lat_vals[lat_vals.length - 1];
    let factor_y = Math.round(src_h/out_h);
    let factor_x = Math.round(src_w/out_w);
    let sub_per_cell = factor_y*factor_x;

    console.log(`[B_anchor] ETOPO shape=${src_h}x${src_w} (${factor_y}x${factor_x} sub-cells/pixel), south_to_north=${south_to_north}`);

    let f_etopo_land = new Float32Array(out_w*out_h);
    let f_ice_etopo = new Float32Array(out_w*out_h);
    let bed_halfm_path = path.join(paths.intermediate_dir, "etopo_bed_halfm.i16");
    let bed_fd = fs.openSync(bed_halfm_path, "w");
    let band_out_rows = 30;
    let t0 = Date.now();

    for (let out_y0 = 0; out_y0 < out_h; out_y0 += band_out_rows) {
      let out_y1 = Math.min(out_h, out_y0 + band_out_rows);
      let n_out = out_y1 - out_y0;
      let sy_count = n_out*factor_y;

      //Map north-to-south output rows to source NetCDF row range
      let src_y0 = south_to_north ? (src_h - out_y1*factor_y) : (out_y0*factor_y);
      let src_y1 = src_y0 + sy_count;

      let bed_slab = bed_nc.getSlice(bed_var, [[src_y0, src_y1], [0, src_w]]);
      let surf_slab = surf_nc.getSlice(surf_var, [[src_y0, src_y1], [0, src_w]]);
      let bed_i16_ns = new Int16Array(sy_count*src_w);

      for (let local_oy = 0; local_oy < n_out; local_oy++) {
        let oy = out_y0 + local_oy;
        let out_row_off = oy*out_w;
        let lat_deg = 90.0 - (oy + 0.5)*(180.0/out_h);
        let is_polar_ice_zone = (lat_deg < -55.0 || lat_deg > 58.0);

        for (let fy = 0; fy < factor_y; fy++) {
          let ns_sy = local_oy*factor_y + fy;
          let slab_sy = south_to_north ? (sy_count - 1 - ns_sy) : ns_sy;
          let slab_row_off = slab_sy*src_w;
          let i16_row_off = ns_sy*src_w;

          for (let ox = 0; ox < out_w; ox++) {
            let sx_base = ox*factor_x;
            let land_hits = 0;
            let ice_hits = 0;

            for (let fx = 0; fx < factor_x; fx++) {
              let sx = sx_base + fx;
              let z_bed = bed_slab[slab_row_off + sx];
              let z_surf = surf_slab[slab_row_off + sx];

              let halfm = Math.round(z_bed*2.0);
              if (halfm < -32767) halfm = -32767;
              if (halfm > 32767) halfm = 32767;
              bed_i16_ns[i16_row_off + sx] = halfm;

              if (z_bed >= 0.0) land_hits++;

              let ice_thick = z_surf - z_bed;
              if (is_polar_ice_zone && ice_thick > ice_thresh && z_surf > 0.0) {
                //Grounded ice sheet (bed >= 0 or grounded below sea level) or permanent floating shelf (bed < 0, surf > 0)
                if (include_shelves || z_bed >= 0.0 || z_surf > 100.0)
                  ice_hits++;
              }
            }

            f_etopo_land[out_row_off + ox] += land_hits;
            f_ice_etopo[out_row_off + ox] += ice_hits;
          }
        }

        for (let ox = 0; ox < out_w; ox++) {
          f_etopo_land[out_row_off + ox] /= sub_per_cell;
          f_ice_etopo[out_row_off + ox] /= sub_per_cell;
        }
      }

      fs.writeSync(bed_fd, Buffer.from(bed_i16_ns.buffer, bed_i16_ns.byteOffset, bed_i16_ns.byteLength));

      if (out_y0 % 360 === 0)
        console.log(`[B_anchor]   ETOPO processed rows ${out_y0}..${out_y1}/${out_h} (${((Date.now() - t0)/1000).toFixed(1)}s)`);
    }

    fs.closeSync(bed_fd);
    bed_nc.close();
    surf_nc.close();

    let meta = {
      factor_x: factor_x,
      factor_y: factor_y,
      scale_m_per_unit: 0.5,
      src_h: src_h,
      src_w: src_w,
    };
    fs.writeFileSync(path.join(paths.intermediate_dir, "etopo_meta.json"), JSON.stringify(meta, null, 2));
    console.log(`[B_anchor] ETOPO 2022 processing completed in ${((Date.now() - t0)/1000).toFixed(1)}s.`);

    //Return statement
    return {
      f_etopo_land: f_etopo_land,
      f_ice_etopo: f_ice_etopo,
      meta: meta,
    };
  }

  /**
   * Executes Step B to build and cache all modern 2025 anchor arrays in intermediate/.
   *
   * @returns {Promise<void>}
   */
  static async run () {
    //Declare local instance variables
    let cfg = GeoGridUtils.loadConfig("model_config.json5");
    let cell_area_by_row = GeoGridUtils.computeCellAreaByRow(cfg.grid);
    let paths = GeoGridUtils.resolvePaths("landarea");
    let total_px = cfg.grid.width*cfg.grid.height;

    if (!fs.existsSync(paths.intermediate_dir))
      fs.mkdirSync(paths.intermediate_dir, { recursive: true });

    //1. Inclusive 2025 OSM Land Polygons (S=20 -> 15" + exact micro-island preservation)
    let osm_out_path = path.join(paths.intermediate_dir, "f_land_osm_2025.f32");
    if (!fs.existsSync(osm_out_path)) {
      let osm_shp = this.findShapefile(path.join(paths.raw_dir, "osm_land_polygons"));
      let f_osm = await this.rasteriseShapefile(osm_shp, cell_area_by_row, {
        label: "OSM Land Polygons",
        min_micro_area_km2: cfg.grid.min_island_area_km2,
        supersample: cfg.grid.vector_supersample,
      });
      GeoGridUtils.writeFloat32Bin(osm_out_path, f_osm);
    } else {
      console.log("[B_anchor] Cached: f_land_osm_2025.f32");
    }

    //2. HydroLAKES v1.0 Natural & Regulated Natural Lakes (Lake_type != 2 excludes impounded reservoirs)
    let force_lakes = process.argv.includes("--force") || process.argv.includes("--force-lakes");
    let lake_out_path = path.join(paths.intermediate_dir, "f_natlake_2025.f32");
    if (force_lakes || !fs.existsSync(lake_out_path)) {
      let lake_shp = this.findShapefile(path.join(paths.raw_dir, "hydrolakes_v10"));
      let f_lakes = await this.rasteriseShapefile(lake_shp, cell_area_by_row, {
        filter_fn: (arg0_feat) => {
          let props = arg0_feat.properties ? arg0_feat.properties : {};
          //Lake_type: 1 = Lake, 2 = Reservoir (handled temporally by GDW v1.0), 3 = Lake control (natural lake)
          if (Number(props.Lake_type) === 2) return false;
          if (props.Lake_area !== undefined && Number(props.Lake_area) < cfg.hydrography.min_lake_area_km2) return false;
          return true;
        },
        hole_filter_fn: (arg0_feat) => {
          let props = arg0_feat.properties ? arg0_feat.properties : {};
          //R3-D03: Inside large boreal lakes (>=500 km2, e.g. Lake Mistassini), drop sub-pixel <25 km2 micro-islets
          //while preserving all real islands >=25 km2 so narrow fjord-lake basins remain crisp open water at 5'
          return (Number(props.Lake_area || 0) >= 500.0) ? 25.0 : 0.0;
        },
        label: "HydroLAKES v1.0 (Natural Lakes)",
        min_micro_area_km2: cfg.hydrography.min_lake_area_km2,
        supersample: cfg.grid.vector_supersample,
      });
      GeoGridUtils.writeFloat32Bin(lake_out_path, f_lakes);
    } else {
      console.log("[B_anchor] Cached: f_natlake_2025.f32");
    }

    //3. ETOPO 2022 30" Bedrock DEM, ETOPO Baseline Land Fraction, and Polar Ice Sheets
    let etopo_land_path = path.join(paths.intermediate_dir, "f_etopo_land_2025.f32");
    let ice_etopo_path = path.join(paths.intermediate_dir, "f_ice_etopo_2025.f32");
    let f_ice_etopo;

    if (!fs.existsSync(etopo_land_path) || !fs.existsSync(ice_etopo_path)) {
      let etopo_res = await this.processEtopo2022(cfg, paths);
      GeoGridUtils.writeFloat32Bin(etopo_land_path, etopo_res.f_etopo_land);
      GeoGridUtils.writeFloat32Bin(ice_etopo_path, etopo_res.f_ice_etopo);
      f_ice_etopo = etopo_res.f_ice_etopo;
    } else {
      console.log("[B_anchor] Cached: f_etopo_land_2025.f32 & f_ice_etopo_2025.f32");
      f_ice_etopo = GeoGridUtils.readFloat32Bin(ice_etopo_path);
    }

    //4. Fuse Modern Ice: ETOPO 2022 Polar Sheets + Natural Earth Mountain Glaciers + Permanent Antarctic Ice Shelves
    let force_ice = process.argv.includes("--force") || process.argv.includes("--force-ice");
    let ice_out_path = path.join(paths.intermediate_dir, "f_ice_2025.f32");
    if (force_ice || !fs.existsSync(ice_out_path)) {
      let gl_shp = this.findShapefile(path.join(paths.raw_dir, "ne_10m_glaciated_areas"));
      let f_ne_gl = await this.rasteriseShapefile(gl_shp, cell_area_by_row, {
        label: "Natural Earth 10m Mountain Glaciers",
        supersample: cfg.grid.vector_supersample,
      });

      let f_ne_shelves = new Float32Array(total_px);
      if (cfg.encoding.include_permanent_ice_shelves) {
        let sh_shp = this.findShapefile(path.join(paths.raw_dir, "ne_10m_antarctic_ice_shelves"));
        f_ne_shelves = await this.rasteriseShapefile(sh_shp, cell_area_by_row, {
          label: "Natural Earth 10m Antarctic Permanent Ice Shelves",
          supersample: cfg.grid.vector_supersample,
        });
      }

      //In Greenland proper and Antarctica (where ETOPO 2022 30" already resolves BedMachine ice thickness at ~900m),
      //prefer ETOPO 30" + explicit permanent Antarctic shelves so coastal ice-free fjords are preserved sharply.
      //Everywhere outside Greenland proper & Antarctica (incl. Iceland, Ellesmere/Baffin, NZ Southern Alps, Alps,
      //Himalayas, Patagonia, Alaska, Svalbard), union max(ETOPO, NE_glaciers).
      let f_ice_2025 = new Float32Array(total_px);
      let height = cfg.grid.height;
      let width = cfg.grid.width;

      for (let y = 0; y < height; y++) {
        let lat_deg = 90.0 - (y + 0.5)*(180.0/height);
        let row_off = y*width;

        for (let x = 0; x < width; x++) {
          let idx = row_off + x;
          let lon_deg = -180.0 + (x + 0.5)*(360.0/width);
          let in_greenland = this.isGreenlandProper(lat_deg, lon_deg);
          let in_antarctica = (lat_deg <= -60.0);

          let v = f_ice_etopo[idx];
          if (in_antarctica) {
            if (f_ne_shelves[idx] > v) v = f_ne_shelves[idx];
          } else if (!in_greenland) {
            if (f_ne_gl[idx] > v) v = f_ne_gl[idx];
          }
          f_ice_2025[idx] = Math.min(1.0, Math.max(0.0, v));
        }
      }

      GeoGridUtils.writeFloat32Bin(ice_out_path, f_ice_2025);
    } else {
      console.log("[B_anchor] Cached: f_ice_2025.f32");
    }

    //Log 2025 anchor summary
    let f_osm_chk = GeoGridUtils.readFloat32Bin(osm_out_path);
    let f_lake_chk = GeoGridUtils.readFloat32Bin(lake_out_path);
    let f_ice_chk = GeoGridUtils.readFloat32Bin(ice_out_path);
    let sum_osm = 0.0, sum_lake = 0.0, sum_ice = 0.0;

    for (let y = 0; y < cfg.grid.height; y++) {
      let a = cell_area_by_row[y];
      let row_off = y*cfg.grid.width;
      for (let x = 0; x < cfg.grid.width; x++) {
        let idx = row_off + x;
        sum_osm += f_osm_chk[idx]*a;
        sum_lake += f_lake_chk[idx]*a;
        sum_ice += f_ice_chk[idx]*a;
      }
    }

    console.log(`[B_anchor] Summary 2025 Anchor -> OSM Land: ${(sum_osm/1e6).toFixed(2)} M km2 | Natural Lakes: ${(sum_lake/1e6).toFixed(2)} M km2 | Glaciers + Permanent Ice Shelves: ${(sum_ice/1e6).toFixed(2)} M km2`);
  }
}

if (require.main === module) {
  B_BuildModernAnchor.run().catch((err) => {
    console.error("[B_anchor] Fatal error:", err);
    process.exit(1);
  });
}

module.exports = B_BuildModernAnchor;
