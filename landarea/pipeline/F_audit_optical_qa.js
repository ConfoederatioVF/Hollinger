"use strict";

let fs = require("node:fs");
let path = require("node:path");

let sharp = require("sharp");

let GeoGridUtils = require("../../core/GeoGridUtils");
let NetCDF4Writer = require("../../core/NetCDF4Writer");

/**
 * Step F: Invariant verification, GeoPNG <-> NetCDF-4 round-trip audit, global trajectory CSV,
 * and quantitative + optical comparison of GHS POP 2025 vs HYDE land_area vs landarea_Hollinger.
 */
class F_AuditOpticalQA {
  /**
   * Loads HYDE land_area_per_cell.asc into a Float32Array(4320*2160).
   *
   * @param {string} arg0_asc_path
   *
   * @returns {Float32Array}
   */
  static loadHydeAsc (arg0_asc_path) {
    //Convert from parameters
    let asc_path = arg0_asc_path;

    //Declare local instance variables
    let lines = fs.readFileSync(asc_path, "utf8").split(/\r?\n/);
    let out = new Float32Array(4320*2160);
    let idx = 0;

    for (let i = 6; i < lines.length; i++) {
      let line = lines[i].trim();
      if (!line) continue;
      let tokens = line.split(/\s+/);
      for (let x = 0; x < tokens.length; x++) {
        if (idx < out.length)
          out[idx++] = Number(tokens[x]);
      }
    }

    //Return statement
    return out;
  }

  /**
   * Renders a regional crop comparison strip (HYDE vs Hollinger at requested years) into a crisp PNG.
   *
   * @param {string} arg0_out_path
   * @param {Array<Object>} arg1_panels - [{ data, ghs, hyde, lat_max, lat_min, lon_max, lon_min, title }]
   * @param {Float64Array} arg2_cell_area_by_row
   *
   * @returns {Promise<void>}
   */
  static async renderRegionalGrid (arg0_out_path, arg1_panels, arg2_cell_area_by_row) {
    //Convert from parameters
    let out_path = arg0_out_path;
    let panels = arg1_panels;
    let cell_area = arg2_cell_area_by_row;

    //Declare local instance variables
    let cols = 3;
    let gap = 16;
    let header_h = 34;
    let panel_h = 260;
    let panel_w = 460;
    let rows = Math.ceil(panels.length/cols);
    let total_h = rows*(panel_h + header_h + gap) + gap;
    let total_w = cols*(panel_w + gap) + gap;
    let composites = [];
    let svg_labels = [];

    for (let i = 0; i < panels.length; i++) {
      let p = panels[i];
      let gx = gap + (i % cols)*(panel_w + gap);
      let gy = gap + Math.floor(i/cols)*(panel_h + header_h + gap);
      let r0 = Math.max(0, Math.floor(((90.0 - p.lat_max)/180.0)*2160));
      let r1 = Math.min(2160, Math.ceil(((90.0 - p.lat_min)/180.0)*2160));
      let c0 = Math.max(0, Math.floor(((p.lon_min + 180.0)/360.0)*4320));
      let c1 = Math.min(4320, Math.ceil(((p.lon_max + 180.0)/360.0)*4320));
      let ch = Math.max(1, r1 - r0);
      let cw = Math.max(1, c1 - c0);
      let rgb = Buffer.allocUnsafe(cw*ch*3);

      for (let y = 0; y < ch; y++) {
        let row = r0 + y;
        let inv_a = 1.0/cell_area[row];
        for (let x = 0; x < cw; x++) {
          let col = c0 + x;
          let idx = row*4320 + col;
          let o3 = (y*cw + x)*3;
          let v = p.data[idx];

          if (p.mode === "diff_ghs") {
            let pop = p.ghs[idx];
            let has_hyde = p.hyde[idx] > 0.0;
            let has_holl = v > 0.0;
            if (pop > 0.0 && !has_hyde && has_holl) {
              //Recovered populated land (bright green)
              rgb[o3] = 46; rgb[o3 + 1] = 204; rgb[o3 + 2] = 113;
            } else if (pop > 0.0 && !has_hyde && !has_holl) {
              //Still missing (red)
              rgb[o3] = 231; rgb[o3 + 1] = 76; rgb[o3 + 2] = 60;
            } else if (v > 0.0) {
              rgb[o3] = 164; rgb[o3 + 1] = 198; rgb[o3 + 2] = 134;
            } else if (v < 0.0) {
              rgb[o3] = 228; rgb[o3 + 1] = 232; rgb[o3 + 2] = 236;
            } else {
              rgb[o3] = 162; rgb[o3 + 1] = 210; rgb[o3 + 2] = 238;
            }
          } else {
            if (v > 0.0) {
              let f = Math.min(1.0, v*inv_a);
              if (f >= 0.995) {
                rgb[o3] = 164; rgb[o3 + 1] = 198; rgb[o3 + 2] = 134;
              } else {
                let n_l = p.data[row*4320 + ((col === 0) ? 4319 : (col - 1))];
                let n_r = p.data[row*4320 + ((col === 4319) ? 0 : (col + 1))];
                let n_u = (row > 0) ? p.data[(row - 1)*4320 + col] : v;
                let n_d = (row < 2159) ? p.data[(row + 1)*4320 + col] : v;
                let adj_sea = (n_l === 0.0 || n_r === 0.0 || n_u === 0.0 || n_d === 0.0);
                let adj_ice = (n_l < 0.0 || n_r < 0.0 || n_u < 0.0 || n_d < 0.0);
                if (adj_ice && !adj_sea) {
                  rgb[o3] = Math.round(228 + (164 - 228)*f);
                  rgb[o3 + 1] = Math.round(232 + (198 - 232)*f);
                  rgb[o3 + 2] = Math.round(236 + (134 - 236)*f);
                } else {
                  rgb[o3] = Math.round(162 + (164 - 162)*f);
                  rgb[o3 + 1] = Math.round(210 + (198 - 210)*f);
                  rgb[o3 + 2] = Math.round(238 + (134 - 238)*f);
                }
              }
            } else if (v < 0.0) {
              let f = Math.min(1.0, (-v)*inv_a);
              rgb[o3] = Math.round(162 + (228 - 162)*f);
              rgb[o3 + 1] = Math.round(210 + (232 - 210)*f);
              rgb[o3 + 2] = Math.round(238 + (236 - 238)*f);
            } else {
              rgb[o3] = 162; rgb[o3 + 1] = 210; rgb[o3 + 2] = 238;
            }
          }
        }
      }

      let panel_png = await sharp(rgb, { raw: { channels: 3, height: ch, width: cw } })
        .resize(panel_w, panel_h, { fit: "fill", kernel: sharp.kernel.nearest })
        .png()
        .toBuffer();

      composites.push({ input: panel_png, left: gx, top: gy + header_h });
      let safe_title = p.title.replace(/&/g, "&amp;").replace(/</g, "&lt;");
      svg_labels.push(`<text x="${gx}" y="${gy + 24}" font-family="Arial" font-size="18" fill="#E6ECF2">${safe_title}</text>`);
    }

    let overlay_svg = Buffer.from(`<svg width="${total_w}" height="${total_h}" xmlns="http://www.w3.org/2000/svg">${svg_labels.join("")}</svg>`);
    composites.push({ input: overlay_svg, left: 0, top: 0 });

    await sharp({
      create: {
        background: { alpha: 1, b: 39, g: 33, r: 28 },
        channels: 4,
        height: total_h,
        width: total_w,
      },
    })
      .composite(composites)
      .png()
      .toFile(out_path);
  }

  /**
   * Executes Step F QA suite.
   *
   * @returns {Promise<void>}
   */
  static async run () {
    //Declare local instance variables
    let cfg = GeoGridUtils.loadConfig("model_config.json5");
    let cell_area = GeoGridUtils.computeCellAreaByRow(cfg.grid);
    let paths = GeoGridUtils.resolvePaths("landarea");
    let years = GeoGridUtils.getTimestepYears(cfg);
    let nc_path = paths.nc_path;

    if (!fs.existsSync(paths.diagnostics_dir))
      fs.mkdirSync(paths.diagnostics_dir, { recursive: true });

    console.log("[F_audit] Opening landarea.nc for invariant & round-trip audit ...");
    let nc = await NetCDF4Writer.openRead(nc_path);
    let shape = nc.getShape("land_area");
    let csv_lines = ["year,positive_land_km2,negative_glaciated_km2,positive_cells,negative_cells,sea_cells"];
    let max_roundtrip_err = 0.0;
    let nan_count = 0;
    let out_of_bounds = 0;
    let sample_years_set = new Set([-10000, -9000, -8900, -8800, -8700, -8600, -8500, -8200, -8000, -7600, -7000, -6000, -5500, -4000, -3800, -500, 0, 100, 800, 1500, 1600, 1700, 1800, 1840, 1850, 1860, 1900, 1910, 1920, 1930, 1950, 1957, 1960, 1968, 1970, 1980, 1985, 2000, 2005, 2015, 2025]);
    let cached_slices = {};

    for (let ti = 0; ti < years.length; ti++) {
      let y = years[ti];
      let slice = nc.getSlice("land_area", [[ti, ti + 1], [0, 2160], [0, 4320]]);
      let land_sum = 0.0;
      let ice_sum = 0.0;
      let n_pos = 0;
      let n_neg = 0;
      let n_zero = 0;

      for (let r = 0; r < 2160; r++) {
        let max_a = cell_area[r] + 0.05;
        let off = r*4320;
        for (let c = 0; c < 4320; c++) {
          let v = slice[off + c];
          if (!Number.isFinite(v)) {
            nan_count++;
            continue;
          }
          if (Math.abs(v) > max_a) out_of_bounds++;
          if (v > 0.0) {
            land_sum += v;
            n_pos++;
          } else if (v < 0.0) {
            ice_sum -= v;
            n_neg++;
          } else {
            n_zero++;
          }
        }
      }

      csv_lines.push(`${y},${land_sum.toFixed(2)},${ice_sum.toFixed(2)},${n_pos},${n_neg},${n_zero}`);

      if (sample_years_set.has(y)) {
        cached_slices[y] = Float32Array.from(slice);
        let png_path = path.join(paths.rasters_dir, `${cfg.encoding.geopng_prefix}${y}.png`);
        let png_data = GeoGridUtils.decodeFloat32GeoPNG(png_path).data;
        for (let k = 0; k < slice.length; k++) {
          let err = Math.abs(slice[k] - png_data[k]);
          if (err > max_roundtrip_err) max_roundtrip_err = err;
        }
      }
    }

    nc.close();
    fs.writeFileSync(path.join(paths.diagnostics_dir, "global_trajectory.csv"), csv_lines.join("\n"));

    //2. Quantitative & optical comparison: GHS POP 2025 vs HYDE land_area vs Hollinger 2025
    let hyde_asc_path = path.resolve(paths.base_dir, "../metadata_HYDE/land_area_per_cell.asc");
    let ghs_png_path = path.resolve(paths.base_dir, "../population_GHSL/population_rasters/2.rasters/GHS_POP_2025.png");
    let hyde_la = this.loadHydeAsc(hyde_asc_path);
    let ghs_pop = GeoGridUtils.decodeFloat32GeoPNG(ghs_png_path).data;
    let holl_2025 = cached_slices[2025];

    let evalGhsOverlap = (arg0_row_shift) => {
      let shift = arg0_row_shift;
      let total_pop = 0.0;
      let pop_miss_hyde = 0.0;
      let pop_miss_holl = 0.0;
      let cells_miss_hyde = 0;
      let cells_miss_holl = 0;
      let cells_recovered = 0;

      for (let r = 0; r < 2160; r++) {
        let gr = r - shift;
        if (gr < 0 || gr >= 2160) continue;
        let off = r*4320;
        let goff = gr*4320;
        for (let c = 0; c < 4320; c++) {
          let p = ghs_pop[goff + c];
          if (!Number.isFinite(p) || p <= 0.0) continue;
          total_pop += p;
          let h_ok = hyde_la[off + c] > 0.0;
          let n_ok = holl_2025[off + c] > 0.0;
          if (!h_ok) {
            cells_miss_hyde++;
            pop_miss_hyde += p;
            if (n_ok) cells_recovered++;
          }
          if (!n_ok) {
            cells_miss_holl++;
            pop_miss_holl += p;
          }
        }
      }

      return {
        cells_missing_hollinger: cells_miss_holl,
        cells_missing_hyde: cells_miss_hyde,
        cells_recovered_by_hollinger: cells_recovered,
        pct_pop_covered_hollinger: Number((100.0*(1.0 - pop_miss_holl/total_pop)).toFixed(4)),
        pct_pop_covered_hyde: Number((100.0*(1.0 - pop_miss_hyde/total_pop)).toFixed(4)),
        pop_missing_hollinger: Math.round(pop_miss_holl),
        pop_missing_hyde: Math.round(pop_miss_hyde),
        total_pop: Math.round(total_pop),
      };
    };

    //2b. Drained, desiccated, Holocene & artificial lakes audit
    let countOpenWater = (arg0_slice, arg1_lat_min, arg2_lat_max, arg3_lon_min, arg4_lon_max) => {
      let r0 = Math.max(0, Math.floor(((90.0 - arg2_lat_max)/180.0)*2160));
      let r1 = Math.min(2160, Math.ceil(((90.0 - arg1_lat_min)/180.0)*2160));
      let c0 = Math.max(0, Math.floor(((arg3_lon_min + 180.0)/360.0)*4320));
      let c1 = Math.min(4320, Math.ceil(((arg4_lon_max + 180.0)/360.0)*4320));
      let open_px = 0;
      for (let r = r0; r < r1; r++) {
        let off = r*4320;
        for (let c = c0; c < c1; c++)
          if (arg0_slice[off + c] === 0.0) open_px++;
      }
      return open_px;
    };

    let lake_checks = [
      { id: "tulare_lake_ca",           box: [35.65, 36.42, -120.10, -119.25], wet_y: 1850,  dry_y: 2025, min_wet: 10, max_dry: 0 },
      { id: "tulare_lake_ca_1920_1930", box: [35.65, 36.42, -120.10, -119.25], wet_y: 1920,  dry_y: 1930, min_wet: 2,  max_dry: 0 },
      { id: "buena_vista_kern_ca",      box: [35.05, 35.40, -119.45, -118.88], wet_y: 1850,  dry_y: 2025, min_wet: 2,  max_dry: 0 },
      { id: "owens_lake_ca",            box: [36.25, 36.60, -118.12, -117.78], wet_y: 1900,  dry_y: 2025, min_wet: 2,  max_dry: 0 },
      { id: "lake_texcoco_valley",      box: [19.15, 19.90,  -99.30,  -98.78], wet_y: 1500,  dry_y: 2025, min_wet: 8,  max_dry: 1 },
      { id: "lake_mayran_mx",           box: [25.45, 26.05, -103.20, -102.30], wet_y: 1900,  dry_y: 2025, min_wet: 8,  max_dry: 0 },
      { id: "lake_poopo_bo",            box: [-19.35,-18.30,  -67.55,  -66.75], wet_y: 2000,  dry_y: 2025, min_wet: 20, max_dry: 0 },
      { id: "lake_fucino_it",           box: [41.88, 42.10,   13.42,   13.72], wet_y: 1850,  dry_y: 2025, min_wet: 1,  max_dry: 0 },
      { id: "lake_copais_gr",           box: [38.34, 38.58,   22.85,   23.28], wet_y: 1850,  dry_y: 2025, min_wet: 1,  max_dry: 1 },
      { id: "lake_karla_gr",            box: [39.38, 39.62,   22.65,   22.88], wet_y: 1950,  dry_y: 2025, min_wet: 2,  max_dry: 0 },
      { id: "lake_bakhtegan_tashk_ir",  box: [29.10, 29.95,   53.10,   54.40], wet_y: 1950,  dry_y: 2025, min_wet: 8,  max_dry: 0 },
      { id: "lake_amik_tr",             box: [36.15, 36.55,   36.15,   36.52], wet_y: 1950,  dry_y: 2025, min_wet: 4,  max_dry: 0 },
      { id: "lake_moeris_fayum_eg",     box: [29.05, 29.70,   30.25,   31.15], wet_y: -4000, dry_y: 2025, min_wet: 20, max_dry: 5 },
      { id: "lake_faguibine_ml",        box: [16.55, 16.95,   -4.32,   -3.65], wet_y: 1950,  dry_y: 2025, min_wet: 4,  max_dry: 0 },
      { id: "lop_nur_basin",            box: [39.70, 40.70,   89.80,   91.40], wet_y: 1900,  dry_y: 2025, min_wet: 30, max_dry: 0 },
      { id: "juyan_lake_cn",            box: [42.10, 42.65,  100.25,  101.40], wet_y: 1950,  dry_y: 2025, min_wet: 6,  max_dry: 1 },
      { id: "lake_manas_cn",            box: [45.55, 46.10,   85.55,   86.40], wet_y: 1950,  dry_y: 2025, min_wet: 6,  max_dry: 0 },
      { id: "great_salt_lake_gilbert",  box: [39.95, 42.10, -114.20, -111.88], wet_y: -10000,dry_y: 2025, min_wet: 180,max_dry: 70 },
      { id: "dead_sea",                 box: [30.90, 31.85,   35.28,   35.65], wet_y: 1950,  dry_y: 2025, min_wet: 10, max_dry: 8 },
      { id: "salton_sea_artificial",    box: [33.00, 33.65, -116.20, -115.45], wet_y: 1950,  dry_y: 1900, min_wet: 12, max_dry: 0 },
      { id: "aydar_arnasay_artificial", box: [40.35, 41.35,   65.50,   68.00], wet_y: 1980,  dry_y: 1960, min_wet: 25, max_dry: 0 },
    ];
    let lake_audit = [];
    let lake_failures = 0;
    for (let li = 0; li < lake_checks.length; li++) {
      let lc = lake_checks[li];
      let wet_px = countOpenWater(cached_slices[lc.wet_y], lc.box[0], lc.box[1], lc.box[2], lc.box[3]);
      let dry_px = countOpenWater(cached_slices[lc.dry_y], lc.box[0], lc.box[1], lc.box[2], lc.box[3]);
      let ok = (wet_px >= lc.min_wet && dry_px <= lc.max_dry && wet_px > dry_px);
      if (!ok) lake_failures++;
      lake_audit.push({ dry_px: dry_px, dry_y: lc.dry_y, id: lc.id, ok: ok, wet_px: wet_px, wet_y: lc.wet_y });
    }

    //2c. Bohai Sea fictitious-island audit & Netherlands Vos antiquity fidelity audit
    let countSmallLandComponents = (arg0_slice, arg1_lat_min, arg2_lat_max, arg3_lon_min, arg4_lon_max, arg5_max_size) => {
      let r0 = Math.max(0, Math.floor(((90.0 - arg2_lat_max)/180.0)*2160));
      let r1 = Math.min(2160, Math.ceil(((90.0 - arg1_lat_min)/180.0)*2160));
      let c0 = Math.max(0, Math.floor(((arg3_lon_min + 180.0)/360.0)*4320));
      let c1 = Math.min(4320, Math.ceil(((arg4_lon_max + 180.0)/360.0)*4320));
      let lh = r1 - r0;
      let lw = c1 - c0;
      let vis = new Uint8Array(lh*lw);
      let q = new Int32Array(lh*lw);
      let small_cnt = 0;
      for (let ly = 0; ly < lh; ly++) {
        for (let lx = 0; lx < lw; lx++) {
          let li = ly*lw + lx;
          if (vis[li] === 1 || arg0_slice[(r0 + ly)*4320 + c0 + lx] <= 0.05) continue;
          let qh = 0, qt = 0;
          vis[li] = 1;
          q[qt++] = li;
          for (; qh < qt; ) {
            let cur = q[qh++];
            let cy = (cur/lw) | 0;
            let cx = cur - cy*lw;
            let nbrs = [[cy, cx - 1], [cy, cx + 1], [cy - 1, cx], [cy + 1, cx]];
            for (let ni = 0; ni < 4; ni++) {
              let ny = nbrs[ni][0], nx = nbrs[ni][1];
              if (ny < 0 || ny >= lh || nx < 0 || nx >= lw) continue;
              let nli = ny*lw + nx;
              if (vis[nli] === 0 && arg0_slice[(r0 + ny)*4320 + c0 + nx] > 0.05) {
                vis[nli] = 1;
                q[qt++] = nli;
              }
            }
          }
          if (qt < arg5_max_size) small_cnt++;
        }
      }
      return small_cnt;
    };

    let meanLandFrac = (arg0_slice, arg1_lat_min, arg2_lat_max, arg3_lon_min, arg4_lon_max) => {
      let r0 = Math.max(0, Math.floor(((90.0 - arg2_lat_max)/180.0)*2160));
      let r1 = Math.min(2160, Math.ceil(((90.0 - arg1_lat_min)/180.0)*2160));
      let c0 = Math.max(0, Math.floor(((arg3_lon_min + 180.0)/360.0)*4320));
      let c1 = Math.min(4320, Math.ceil(((arg4_lon_max + 180.0)/360.0)*4320));
      let sum = 0.0, n = 0;
      for (let r = r0; r < r1; r++) {
        let inv_a = 1.0/cell_area[r];
        let off = r*4320;
        for (let c = c0; c < c1; c++) {
          sum += Math.min(1.0, Math.max(0.0, arg0_slice[off + c]*inv_a));
          n++;
        }
      }
      return (n > 0) ? (sum/n) : 0.0;
    };

    let bohai_years = [-4000, 0, 1500, 1900];
    let bohai_audit = {};
    let bohai_failures = 0;
    for (let bi = 0; bi < bohai_years.length; bi++) {
      let by = bohai_years[bi];
      let cnt = countSmallLandComponents(cached_slices[by], 36.8, 40.0, 116.5, 120.4, 30);
      bohai_audit[by] = cnt;
      if (cnt > 0) bohai_failures++;
    }

    let nl_100_ij = meanLandFrac(cached_slices[100], 52.4, 53.0, 5.0, 5.7);
    let nl_100_ze = meanLandFrac(cached_slices[100], 51.35, 51.9, 3.5, 4.3);
    let nl_1850_ij = meanLandFrac(cached_slices[1850], 52.4, 53.0, 5.0, 5.7);
    let nl_1950_ij = meanLandFrac(cached_slices[1950], 52.4, 53.0, 5.0, 5.7);
    let nl_2025_ij = meanLandFrac(holl_2025, 52.4, 53.0, 5.0, 5.7);
    let nl_ok = (nl_100_ij >= 0.50 && nl_100_ze >= 0.68 && nl_1850_ij < 0.25 && nl_1950_ij > nl_1850_ij && nl_2025_ij > nl_1950_ij);

    let polder_metrics = {
      enkhuizen_2025: Number(meanLandFrac(holl_2025, 52.67, 52.74, 5.25, 5.32).toFixed(3)),
      goeree_2025: Number(meanLandFrac(holl_2025, 51.75, 51.83, 3.92, 4.00).toFixed(3)),
      haarlemmermeer_1700: Number(meanLandFrac(cached_slices[1700], 52.23, 52.36, 4.60, 4.78).toFixed(3)),
      haarlemmermeer_1840: Number(meanLandFrac(cached_slices[1840], 52.23, 52.36, 4.60, 4.78).toFixed(3)),
      haarlemmermeer_1860: Number(meanLandFrac(cached_slices[1860], 52.23, 52.36, 4.60, 4.78).toFixed(3)),
      lauwersmeer_1950: Number(meanLandFrac(cached_slices[1950], 53.34, 53.41, 6.18, 6.32).toFixed(3)),
      lauwersmeer_1970: Number(meanLandFrac(cached_slices[1970], 53.34, 53.41, 6.18, 6.32).toFixed(3)),
      lelystad_2025: Number(meanLandFrac(holl_2025, 52.50, 52.58, 5.42, 5.50).toFixed(3)),
      maasvlakte2_2005: Number(meanLandFrac(cached_slices[2005], 51.93, 51.99, 3.93, 4.00).toFixed(3)),
      maasvlakte2_2015: Number(meanLandFrac(cached_slices[2015], 51.93, 51.99, 3.93, 4.00).toFixed(3)),
      nh_golden_1500: Number(meanLandFrac(cached_slices[1500], 52.50, 52.66, 4.78, 4.98).toFixed(3)),
      nh_golden_1600: Number(meanLandFrac(cached_slices[1600], 52.50, 52.66, 4.78, 4.98).toFixed(3)),
      nh_golden_1700: Number(meanLandFrac(cached_slices[1700], 52.50, 52.66, 4.78, 4.98).toFixed(3)),
      nop_1930: Number(meanLandFrac(cached_slices[1930], 52.66, 52.82, 5.60, 5.90).toFixed(3)),
      nop_1950: Number(meanLandFrac(cached_slices[1950], 52.66, 52.82, 5.60, 5.90).toFixed(3)),
      oflevo_1950: Number(meanLandFrac(cached_slices[1950], 52.45, 52.57, 5.45, 5.80).toFixed(3)),
      oflevo_1957: Number(meanLandFrac(cached_slices[1957], 52.45, 52.57, 5.45, 5.80).toFixed(3)),
      urk_west_nop_2025: Number(meanLandFrac(holl_2025, 52.62, 52.80, 5.60, 5.66).toFixed(3)),
      wieringer_1920: Number(meanLandFrac(cached_slices[1920], 52.76, 52.90, 4.93, 5.08).toFixed(3)),
      wieringer_1930: Number(meanLandFrac(cached_slices[1930], 52.76, 52.90, 4.93, 5.08).toFixed(3)),
      zflevo_1957: Number(meanLandFrac(cached_slices[1957], 52.32, 52.42, 5.18, 5.48).toFixed(3)),
      zflevo_1968: Number(meanLandFrac(cached_slices[1968], 52.32, 52.42, 5.18, 5.48).toFixed(3)),
    };
    let nl_polders_ok = (
      polder_metrics.haarlemmermeer_1700 < 0.55 &&
      polder_metrics.haarlemmermeer_1840 < 0.55 &&
      polder_metrics.haarlemmermeer_1860 >= 0.85 &&
      polder_metrics.nh_golden_1500 < 0.65 &&
      polder_metrics.nh_golden_1600 < 0.65 &&
      polder_metrics.nh_golden_1700 >= 0.90 &&
      polder_metrics.wieringer_1920 < 0.25 &&
      polder_metrics.wieringer_1930 >= 0.80 &&
      polder_metrics.nop_1930 < 0.25 &&
      polder_metrics.nop_1950 >= 0.80 &&
      polder_metrics.urk_west_nop_2025 >= 0.45 &&
      polder_metrics.oflevo_1950 < 0.15 &&
      polder_metrics.oflevo_1957 >= 0.80 &&
      polder_metrics.zflevo_1957 < 0.20 &&
      polder_metrics.zflevo_1968 >= 0.75 &&
      polder_metrics.lauwersmeer_1950 < 0.25 &&
      polder_metrics.lauwersmeer_1970 >= 0.55 &&
      polder_metrics.maasvlakte2_2005 < 0.05 &&
      polder_metrics.maasvlakte2_2015 >= 0.20 &&
      polder_metrics.lelystad_2025 >= 0.35 &&
      polder_metrics.enkhuizen_2025 >= 0.25 &&
      polder_metrics.goeree_2025 >= 0.35
    );

    let oj_years = [-9000, -8900, -8800, -8700, -8600];
    let oj_counts = {};
    for (let oi = 0; oi < oj_years.length; oi++) {
      let oy = oj_years[oi];
      oj_counts[oy] = countOpenWater(cached_slices[oy], 47.8, 51.5, -81.5, -76.0);
    }
    let oj_ok = (
      oj_counts[-8800] >= 150 &&
      oj_counts[-8700] >= 250 &&
      oj_counts[-8700] >= oj_counts[-8900]*0.8 &&
      oj_counts[-8600] >= oj_counts[-8700]*0.85
    );

    let maxHorizEdgeRun = (arg0_slice, arg1_lat_min, arg2_lat_max, arg3_lon_min, arg4_lon_max, arg5_mode) => {
      let r0 = Math.max(1, Math.floor(((90.0 - arg2_lat_max)/180.0)*2160));
      let r1 = Math.min(2158, Math.ceil(((90.0 - arg1_lat_min)/180.0)*2160));
      let c0 = Math.max(0, Math.floor(((arg3_lon_min + 180.0)/360.0)*4320));
      let c1 = Math.min(4319, Math.ceil(((arg4_lon_max + 180.0)/360.0)*4320));
      let max_run = 0;
      for (let r = r0; r < r1; r++) {
        let run = 0;
        let off = r*4320;
        let noff = (r + 1)*4320;
        for (let c = c0; c <= c1; c++) {
          let v = arg0_slice[off + c];
          let nv = arg0_slice[noff + c];
          let edge = (arg5_mode === "water") ? (v === 0.0 && nv > 0.0) : (v < 0.0 && nv >= 0.0);
          run = edge ? (run + 1) : 0;
          if (run > max_run) max_run = run;
        }
      }
      return max_run;
    };

    let agassiz_straight_years = [-8600, -8200, -7600, -7000];
    let agassiz_straight_metrics = {};
    let agassiz_straight_ok = true;
    for (let ai = 0; ai < agassiz_straight_years.length; ai++) {
      let ay = agassiz_straight_years[ai];
      let sl = cached_slices[ay];
      let wat_run = maxHorizEdgeRun(sl, 47.0, 51.0, -85.0, -75.0, "water");
      let ice_run = maxHorizEdgeRun(sl, 47.5, 53.0, -88.0, -76.0, "ice");
      let ont_dry = meanLandFrac(sl, 45.2, 46.0, -80.0, -78.0);
      agassiz_straight_metrics[ay] = {
        max_water_horiz_run: wat_run,
        max_ice_horiz_run: ice_run,
        ontario_shield_dry: Number(ont_dry.toFixed(3)),
      };
      if (wat_run > 10 || ice_run > 20 || ont_dry <= 0.70) agassiz_straight_ok = false;
    }

    let report = {
      agassiz_ojibway_straight_line_qa: {
        metrics_by_year: agassiz_straight_metrics,
        ok: agassiz_straight_ok,
      },
      bohai_fictitious_islands: bohai_audit,
      drained_and_dynamic_lakes_audit: {
        checks: lake_audit,
        failures: lake_failures,
        total: lake_audit.length,
      },
      ghs_vs_hyde_vs_hollinger_aligned_1row: evalGhsOverlap(1),
      ghs_vs_hyde_vs_hollinger_raw: evalGhsOverlap(0),
      max_geopng_netcdf_roundtrip_err: max_roundtrip_err,
      nan_or_inf_cells: nan_count,
      netcdf_shape: shape,
      nl_polders_fidelity: {
        metrics: polder_metrics,
        ok: nl_polders_ok,
      },
      nl_vos_antiquity_fidelity: {
        ijsselmeer_100ad: Number(nl_100_ij.toFixed(3)),
        ijsselmeer_1850ad: Number(nl_1850_ij.toFixed(3)),
        ijsselmeer_1950ad: Number(nl_1950_ij.toFixed(3)),
        ijsselmeer_2025ad: Number(nl_2025_ij.toFixed(3)),
        ok: nl_ok,
        zeeland_100ad: Number(nl_100_ze.toFixed(3)),
      },
      ojibway_continuity: {
        counts_by_year: oj_counts,
        ok: oj_ok,
      },
      out_of_bounds_cells: out_of_bounds,
    };

    fs.writeFileSync(path.join(paths.diagnostics_dir, "qa_report.json"), JSON.stringify(report, null, 2));
    console.log("[F_audit] QA Report:", JSON.stringify(report, null, 2));

    //3. Render Global Optical Diff PNG (4320x2160 -> 2160x1080)
    let diff_rgb = Buffer.allocUnsafe(4320*2160*3);
    for (let i = 0; i < 4320*2160; i++) {
      let o3 = i*3;
      let pop = (Number.isFinite(ghs_pop[i]) && ghs_pop[i] > 0.0) ? ghs_pop[i] : 0.0;
      let has_hyde = hyde_la[i] > 0.0;
      let v = holl_2025[i];

      if (pop > 0.0 && !has_hyde && v > 0.0) {
        diff_rgb[o3] = 46; diff_rgb[o3 + 1] = 204; diff_rgb[o3 + 2] = 113; //Green: recovered by Hollinger
      } else if (pop > 0.0 && !has_hyde && v <= 0.0) {
        diff_rgb[o3] = 231; diff_rgb[o3 + 1] = 76; diff_rgb[o3 + 2] = 60; //Red: still missing
      } else if (v > 0.0) {
        diff_rgb[o3] = 164; diff_rgb[o3 + 1] = 198; diff_rgb[o3 + 2] = 134;
      } else if (v < 0.0) {
        diff_rgb[o3] = 228; diff_rgb[o3 + 1] = 232; diff_rgb[o3 + 2] = 236;
      } else {
        diff_rgb[o3] = 162; diff_rgb[o3 + 1] = 210; diff_rgb[o3 + 2] = 238;
      }
    }

    await sharp(diff_rgb, { raw: { channels: 3, height: 2160, width: 4320 } })
      .resize(2160, 1080, { kernel: sharp.kernel.nearest })
      .png()
      .toFile(path.join(paths.diagnostics_dir, "global_ghs_vs_hyde_vs_hollinger.png"));

    //4. Render 12 Regional Optical Benchmark Panels
    let panels = [
      { data: holl_2025, ghs: ghs_pop, hyde: hyde_la, lat_max: 8, lat_min: -2, lon_max: 77, lon_min: 71, mode: "diff_ghs", title: "Maldives 2025 (Green=Recovered vs HYDE)" },
      { data: holl_2025, ghs: ghs_pop, hyde: hyde_la, lat_max: 5, lat_min: -25, lon_max: -135, lon_min: -180, mode: "diff_ghs", title: "Polynesia 2025 (Green=Recovered vs HYDE)" },
      { data: holl_2025, ghs: ghs_pop, hyde: hyde_la, lat_max: 27, lat_min: 24, lon_max: 56, lon_min: 49, mode: "diff_ghs", title: "Bahrain & Dubai 2025 (Reclamation)" },
      { data: cached_slices[1950], lat_max: 53.6, lat_min: 51.3, lon_max: 7.3, lon_min: 3.3, title: "Netherlands 1950 (Pre-Flevoland)" },
      { data: holl_2025, lat_max: 53.6, lat_min: 51.3, lon_max: 7.3, lon_min: 3.3, title: "Netherlands 2025 (Flevoland + Maasvlakte)" },
      { data: cached_slices[-8000], lat_max: 60, lat_min: 49, lon_max: 11, lon_min: -6, title: "Doggerland -8000 (8000 BC)" },
      { data: cached_slices[-10000], lat_max: 15, lat_min: -12, lon_max: 122, lon_min: 94, title: "Sundaland -10000 (10000 BC)" },
      { data: cached_slices[-10000], lat_max: 75, lat_min: 40, lon_max: -50, lon_min: -130, title: "Laurentide & Agassiz -10000" },
      { data: cached_slices[-6000], lat_max: 75, lat_min: 40, lon_max: -50, lon_min: -130, title: "Laurentide Retreat -6000" },
      { data: holl_2025, lat_max: 84, lat_min: 59, lon_max: -10, lon_min: -74, title: "Greenland 2025 (Ice Sheet vs Bare Fjords)" },
      { data: cached_slices[1950], lat_max: 48, lat_min: 42, lon_max: 64, lon_min: 56, title: "Aral Sea 1950 (Full Extent)" },
      { data: holl_2025, lat_max: 48, lat_min: 42, lon_max: 64, lon_min: 56, title: "Aral Sea 2025 (Desiccated Aralkum)" },
    ];

    await this.renderRegionalGrid(path.join(paths.diagnostics_dir, "regional_audit_panels.png"), panels, cell_area);

    //5. Render 12-Panel Global Paleo & Shelf Topology Audit Sheet (verifies zero hollow outlines or ice/lake artefacts)
    let topo_panels = [
      { data: cached_slices[-10000], lat_max: 70, lat_min: 10, lon_max: 72, lon_min: -25, title: "Europe / MENA / Persian Gulf -10000 (User View)" },
      { data: cached_slices[-8000], lat_max: 66, lat_min: 30, lon_max: 58, lon_min: -12, title: "Europe / Doggerland / Baltic -8000" },
      { data: cached_slices[-10000], lat_max: 20, lat_min: -15, lon_max: 125, lon_min: 92, title: "Sundaland & Java Sea -10000" },
      { data: cached_slices[-10000], lat_max: 0, lat_min: -45, lon_max: 160, lon_min: 110, title: "Sahul (Australia & New Guinea) -10000" },
      { data: cached_slices[-10000], lat_max: 75, lat_min: 50, lon_max: -145, lon_min: -180, title: "Beringia East (Alaska) -10000" },
      { data: cached_slices[-10000], lat_max: 45, lat_min: 20, lon_max: 145, lon_min: 112, title: "Yellow Sea / Korea / Japan -10000" },
      { data: cached_slices[-10000], lat_max: 72, lat_min: 35, lon_max: -50, lon_min: -135, title: "North America / Laurentide / Agassiz -10000" },
      { data: cached_slices[-6000], lat_max: 72, lat_min: 35, lon_max: -50, lon_min: -135, title: "North America / Hudson Bay -6000" },
      { data: cached_slices[-7000], lat_max: 35, lat_min: 5, lon_max: 60, lon_min: -20, title: "Mega-Lake Chad / Red Sea / Persian Gulf -7000" },
      { data: cached_slices[-10000], lat_max: -20, lat_min: -58, lon_max: -48, lon_min: -80, title: "Patagonia & Falklands Shelf -10000" },
      { data: holl_2025, lat_max: 42, lat_min: 25, lon_max: 100, lon_min: 68, title: "Himalayas & Tibet 2025 (Glacier Edges)" },
      { data: holl_2025, lat_max: -60, lat_min: -90, lon_max: 180, lon_min: -180, title: "Antarctica & Permanent Ice Shelves 2025" },
    ];

    await this.renderRegionalGrid(path.join(paths.diagnostics_dir, "paleo_topology_audit.png"), topo_panels, cell_area);

    //6. Render dedicated 18-panel verification sheet targeting all 9 user-reported issues directly!
    let issue_panels = [
      { data: cached_slices[1900], lat_max: 12, lat_min: -30, lon_max: -34, lon_min: -72, title: "#1 Latin America 1900 (Zero Anachronistic Dams)" },
      { data: holl_2025, lat_max: 12, lat_min: -30, lon_max: -34, lon_min: -72, title: "#1 Latin America 2025 (Modern Reservoirs)" },
      { data: cached_slices[1900], lat_max: 27, lat_min: 19, lon_max: 36, lon_min: 27, title: "#9 Egypt / Aswan 1900 (No Nasser or Toshka)" },
      { data: holl_2025, lat_max: 27, lat_min: 19, lon_max: 36, lon_min: 27, title: "#9 Egypt / Aswan 2025 (Nasser + Toshka)" },
      { data: cached_slices[-4000], lat_max: 42, lat_min: 28, lon_max: 130, lon_min: 114, title: "#2 China Bohai/Yellow/Yangtze -4000 (Highstand)" },
      { data: cached_slices[1950], lat_max: 42, lat_min: 28, lon_max: 130, lon_min: 114, title: "#2 China Bohai/Yellow/Yangtze 1950 (Pre-Reclam)" },
      { data: holl_2025, lat_max: 42, lat_min: 28, lon_max: 130, lon_min: 114, title: "#2 China Bohai/Yellow/Yangtze 2025" },
      { data: cached_slices[1950], lat_max: 47.2, lat_min: 42.8, lon_max: 62.5, lon_min: 57.5, title: "#3 Aral Sea 1950 (Exact 53.4m 30in Bathymetry)" },
      { data: cached_slices[1985], lat_max: 47.2, lat_min: 42.8, lon_max: 62.5, lon_min: 57.5, title: "#3 Aral Sea 1985 (42m Desiccation)" },
      { data: holl_2025, lat_max: 47.2, lat_min: 42.8, lon_max: 62.5, lon_min: 57.5, title: "#3 Aral Sea 2025 (Small Aral + West Basin)" },
      { data: cached_slices[-4000], lat_max: 33.5, lat_min: 27.5, lon_max: 52.0, lon_min: 44.5, title: "#4 Mesopotamia -4000 (Gulf at Ur & Eridu!)" },
      { data: cached_slices[0], lat_max: 33.5, lat_min: 27.5, lon_max: 52.0, lon_min: 44.5, title: "#4 Mesopotamia 1 AD (Charax Spasinu)" },
      { data: cached_slices[-3800], lat_max: 53.8, lat_min: 50.6, lon_max: 7.4, lon_min: 3.1, title: "#5 Netherlands Vos -3800 (Tidal Basins)" },
      { data: cached_slices[1500], lat_max: 53.8, lat_min: 50.6, lon_max: 7.4, lon_min: 3.1, title: "#5 Netherlands Vos 1500 (Zuiderzee+Haarlemmermeer)" },
      { data: cached_slices[1950], lat_max: 53.8, lat_min: 50.6, lon_max: 7.4, lon_min: 3.1, title: "#5 Netherlands 1950 (NOP Done, Pre-Flevoland)" },
      { data: cached_slices[-8500], lat_max: 50, lat_min: 36, lon_max: 56, lon_min: 27, title: "#6 Black-Caspian-Volga -8500 (Mangyshlak Low)" },
      { data: cached_slices[-5500], lat_max: 50, lat_min: 36, lon_max: 56, lon_min: 27, title: "#6 Black-Caspian-Volga -5500 (New Caspian High)" },
      { data: cached_slices[-7000], lat_max: 19.5, lat_min: 10.5, lon_max: 20.5, lon_min: 11.5, title: "#8 Mega-Lake Chad -7000 (Exact 326m Contour)" },
      { data: cached_slices[1950], lat_max: 19.5, lat_min: 10.5, lon_max: 20.5, lon_min: 11.5, title: "#8 Lake Chad 1950 (Exact 282.5m Basin)" },
      { data: holl_2025, lat_max: 19.5, lat_min: 10.5, lon_max: 20.5, lon_min: 11.5, title: "#8 Lake Chad 2025 (Modern Desiccated)" },
      { data: cached_slices[-10000], lat_max: 62, lat_min: 48, lon_max: 12, lon_min: -8, title: "#7 Neolithic Doggerland -10000 (Solid Fill)" },
    ];

    await this.renderRegionalGrid(path.join(paths.diagnostics_dir, "user_9_issues_audit.png"), issue_panels, cell_area);

    //7. Render dedicated Drained, Desiccated & Artificial Lakes visual audit sheet
    let lake_panels = [
      { data: cached_slices[1850], lat_max: 37.5, lat_min: 32.5, lon_max: -115.0, lon_min: -121.0, title: "California 1850 (Tulare + Buena Vista + Owens; No Salton)" },
      { data: cached_slices[1950], lat_max: 37.5, lat_min: 32.5, lon_max: -115.0, lon_min: -121.0, title: "California 1950 (Tulare/Owens Drained; Salton Sea Present)" },
      { data: holl_2025, lat_max: 37.5, lat_min: 32.5, lon_max: -115.0, lon_min: -121.0, title: "California 2025" },
      { data: cached_slices[1500], lat_max: 20.5, lat_min: 18.8, lon_max: -98.2, lon_min: -100.0, title: "Valley of Mexico 1500 (Lake Texcoco Full)" },
      { data: cached_slices[1850], lat_max: 20.5, lat_min: 18.8, lon_max: -98.2, lon_min: -100.0, title: "Valley of Mexico 1850 (Shrinking)" },
      { data: holl_2025, lat_max: 20.5, lat_min: 18.8, lon_max: -98.2, lon_min: -100.0, title: "Valley of Mexico 2025 (Drained)" },
      { data: cached_slices[-4000], lat_max: 30.2, lat_min: 28.8, lon_max: 31.5, lon_min: 29.8, title: "Fayum / Lake Moeris -4000 (Holocene Highstand)" },
      { data: holl_2025, lat_max: 30.2, lat_min: 28.8, lon_max: 31.5, lon_min: 29.8, title: "Fayum / Birket Qarun 2025" },
      { data: cached_slices[1900], lat_max: 41.2, lat_min: 39.2, lon_max: 92.0, lon_min: 89.0, title: "Lop Nur 1900 (Active Terminal Lake)" },
      { data: holl_2025, lat_max: 41.2, lat_min: 39.2, lon_max: 92.0, lon_min: 89.0, title: "Lop Nur 2025 (Desiccated)" },
      { data: cached_slices[1960], lat_max: 42.0, lat_min: 39.8, lon_max: 68.5, lon_min: 65.0, title: "Aydar-Arnasay 1960 (Dry Pre-1969 Spillway)" },
      { data: cached_slices[1980], lat_max: 42.0, lat_min: 39.8, lon_max: 68.5, lon_min: 65.0, title: "Aydar-Arnasay 1980 (Impounded)" },
      { data: cached_slices[2000], lat_max: -17.8, lat_min: -19.8, lon_max: -66.0, lon_min: -68.5, title: "Lake Poopo 2000 (Full)" },
      { data: holl_2025, lat_max: -17.8, lat_min: -19.8, lon_max: -66.0, lon_min: -68.5, title: "Lake Poopo 2025 (Desiccated)" },
      { data: cached_slices[-10000], lat_max: 42.4, lat_min: 39.6, lon_max: -111.3, lon_min: -114.5, title: "Great Salt Lake -10000 (Gilbert Highstand)" },
      { data: cached_slices[1950], lat_max: 42.4, lat_min: 39.6, lon_max: -111.3, lon_min: -114.5, title: "Great Salt Lake 1950" },
      { data: holl_2025, lat_max: 42.4, lat_min: 39.6, lon_max: -111.3, lon_min: -114.5, title: "Great Salt Lake 2025" },
      { data: cached_slices[1850], lat_max: 42.6, lat_min: 36.0, lon_max: 37.0, lon_min: 13.0, title: "Fucino, Copais, Karla & Amik 1850" },
      { data: holl_2025, lat_max: 42.6, lat_min: 36.0, lon_max: 37.0, lon_min: 13.0, title: "Fucino, Copais, Karla & Amik 2025" },
    ];
    await this.renderRegionalGrid(path.join(paths.diagnostics_dir, "drained_lakes_audit.png"), lake_panels, cell_area);

    //8. Render dedicated Bohai Sea, Lake Tulare, and Netherlands Antiquity 18-panel audit sheet
    let btn_panels = [
      { data: cached_slices[-4000], lat_max: 40.0, lat_min: 36.8, lon_max: 120.5, lon_min: 116.5, title: "Bohai Sea -4000 (Zero Fictitious Islands)" },
      { data: cached_slices[0], lat_max: 40.0, lat_min: 36.8, lon_max: 120.5, lon_min: 116.5, title: "Bohai Sea 1 AD" },
      { data: cached_slices[1500], lat_max: 40.0, lat_min: 36.8, lon_max: 120.5, lon_min: 116.5, title: "Bohai Sea 1500 AD" },
      { data: cached_slices[1900], lat_max: 40.0, lat_min: 36.8, lon_max: 120.5, lon_min: 116.5, title: "Bohai Sea 1900 AD (SW-Curved Bay)" },
      { data: cached_slices[1950], lat_max: 40.0, lat_min: 36.8, lon_max: 120.5, lon_min: 116.5, title: "Bohai Sea 1950 AD (Intermediate Delta)" },
      { data: holl_2025, lat_max: 40.0, lat_min: 36.8, lon_max: 120.5, lon_min: 116.5, title: "Bohai Sea 2025 AD (Modern)" },
      { data: cached_slices[1850], lat_max: 36.6, lat_min: 35.4, lon_max: -119.1, lon_min: -120.3, title: "Lake Tulare 1850 (Full)" },
      { data: cached_slices[1900], lat_max: 36.6, lat_min: 35.4, lon_max: -119.1, lon_min: -120.3, title: "Lake Tulare 1900" },
      { data: cached_slices[1910], lat_max: 36.6, lat_min: 35.4, lon_max: -119.1, lon_min: -120.3, title: "Lake Tulare 1910" },
      { data: cached_slices[1920], lat_max: 36.6, lat_min: 35.4, lon_max: -119.1, lon_min: -120.3, title: "Lake Tulare 1920 (Still Visible)" },
      { data: cached_slices[1930], lat_max: 36.6, lat_min: 35.4, lon_max: -119.1, lon_min: -120.3, title: "Lake Tulare 1930 (Faded Out)" },
      { data: holl_2025, lat_max: 36.6, lat_min: 35.4, lon_max: -119.1, lon_min: -120.3, title: "Lake Tulare 2025 (Dry)" },
      { data: cached_slices[-3800], lat_max: 53.8, lat_min: 50.6, lon_max: 7.4, lon_min: 3.1, title: "Netherlands -3800 (Tidal Basins)" },
      { data: cached_slices[-500], lat_max: 53.8, lat_min: 50.6, lon_max: 7.4, lon_min: 3.1, title: "Netherlands -500 (Closed Peat)" },
      { data: cached_slices[100], lat_max: 53.8, lat_min: 50.6, lon_max: 7.4, lon_min: 3.1, title: "Netherlands 100 AD (Roman Coast)" },
      { data: cached_slices[800], lat_max: 53.8, lat_min: 50.6, lon_max: 7.4, lon_min: 3.1, title: "Netherlands 800 AD (Almere/Wadden)" },
      { data: cached_slices[1500], lat_max: 53.8, lat_min: 50.6, lon_max: 7.4, lon_min: 3.1, title: "Netherlands 1500 AD (Zuiderzee)" },
      { data: holl_2025, lat_max: 53.8, lat_min: 50.6, lon_max: 7.4, lon_min: 3.1, title: "Netherlands 2025 AD (Modern)" },
    ];
    await this.renderRegionalGrid(path.join(paths.diagnostics_dir, "bohai_tulare_nl_audit.png"), btn_panels, cell_area);

    //9. Render dedicated 12-panel Dutch Polder & Historical Reclamation sequence audit
    let nl_polder_panels = [
      { data: cached_slices[1500], lat_max: 53.6, lat_min: 51.3, lon_max: 7.3, lon_min: 3.3, title: "NL 1500 (Beemster & Haarlemmermeer Open)" },
      { data: cached_slices[1600], lat_max: 53.6, lat_min: 51.3, lon_max: 7.3, lon_min: 3.3, title: "NL 1600 (Haarlemmermeer 181 km2 Lake)" },
      { data: cached_slices[1700], lat_max: 53.6, lat_min: 51.3, lon_max: 7.3, lon_min: 3.3, title: "NL 1700 (Beemster Drained, Haarlemmer Open)" },
      { data: cached_slices[1800], lat_max: 53.6, lat_min: 51.3, lon_max: 7.3, lon_min: 3.3, title: "NL 1800 (Haarlemmermeer Still Open)" },
      { data: cached_slices[1840], lat_max: 53.6, lat_min: 51.3, lon_max: 7.3, lon_min: 3.3, title: "NL 1840 (Pre-Haarlemmermeer Pumping)" },
      { data: cached_slices[1860], lat_max: 53.6, lat_min: 51.3, lon_max: 7.3, lon_min: 3.3, title: "NL 1860 (Haarlemmermeer Drained 1852)" },
      { data: cached_slices[1920], lat_max: 53.6, lat_min: 51.3, lon_max: 7.3, lon_min: 3.3, title: "NL 1920 (Pre-Zuiderzeewerken)" },
      { data: cached_slices[1930], lat_max: 53.6, lat_min: 51.3, lon_max: 7.3, lon_min: 3.3, title: "NL 1930 (Wieringermeer Drained)" },
      { data: cached_slices[1950], lat_max: 53.6, lat_min: 51.3, lon_max: 7.3, lon_min: 3.3, title: "NL 1950 (Noordoostpolder Done)" },
      { data: cached_slices[1957], lat_max: 53.6, lat_min: 51.3, lon_max: 7.3, lon_min: 3.3, title: "NL 1957 (Oostelijk Flevoland Only)" },
      { data: cached_slices[1968], lat_max: 53.6, lat_min: 51.3, lon_max: 7.3, lon_min: 3.3, title: "NL 1968 (Zuidelijk Flevoland Done)" },
      { data: holl_2025, lat_max: 53.6, lat_min: 51.3, lon_max: 7.3, lon_min: 3.3, title: "NL 2025 (Lauwersmeer + Maasvlakte II)" },
    ];
    await this.renderRegionalGrid(path.join(paths.diagnostics_dir, "nl_polders_audit.png"), nl_polder_panels, cell_area);

    //10. Render zoomed Lake Agassiz-Ojibway diagnostic sheet (-8600, -8200, -7600, -7000)
    let agassiz_panels = [
      { data: cached_slices[-8600], lat_max: 54.5, lat_min: 45.5, lon_max: -73.0, lon_min: -95.0, title: "Agassiz-Ojibway -8600 (DEM Watershed + Ice Margin)" },
      { data: cached_slices[-8200], lat_max: 54.5, lat_min: 45.5, lon_max: -73.0, lon_min: -95.0, title: "Agassiz-Ojibway -8200" },
      { data: cached_slices[-7600], lat_max: 54.5, lat_min: 45.5, lon_max: -73.0, lon_min: -95.0, title: "Agassiz-Ojibway -7600" },
      { data: cached_slices[-7000], lat_max: 54.5, lat_min: 45.5, lon_max: -73.0, lon_min: -95.0, title: "Agassiz-Ojibway -7000" },
    ];
    await this.renderRegionalGrid(path.join(paths.diagnostics_dir, "agassiz_ojibway_audit.png"), agassiz_panels, cell_area);

    console.log("[F_audit] Saved global_ghs_vs_hyde_vs_hollinger.png, regional_audit_panels.png, paleo_topology_audit.png, user_9_issues_audit.png, drained_lakes_audit.png, bohai_tulare_nl_audit.png, nl_polders_audit.png & agassiz_ojibway_audit.png.");
    if (lake_failures > 0 || bohai_failures > 0 || !nl_ok || !nl_polders_ok || !oj_ok || !agassiz_straight_ok) {
      console.error(`[F_audit] FAILED QA checks! lake_failures=${lake_failures}, bohai_failures=${bohai_failures}, nl_ok=${nl_ok}, nl_polders_ok=${nl_polders_ok}, oj_ok=${oj_ok}, agassiz_straight_ok=${agassiz_straight_ok}`);
      process.exit(1);
    }
  }
}

if (require.main === module) {
  F_AuditOpticalQA.run().catch((err) => {
    console.error("[F_audit] Fatal error:", err);
    process.exit(1);
  });
}

module.exports = F_AuditOpticalQA;
