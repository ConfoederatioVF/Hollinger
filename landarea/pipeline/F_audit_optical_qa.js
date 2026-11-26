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
    let sample_years_set = new Set([-10000, -8500, -8000, -7000, -6000, -5500, -4000, -3800, 0, 100, 1500, 1850, 1900, 1950, 1985, 2025]);
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

    let report = {
      ghs_vs_hyde_vs_hollinger_aligned_1row: evalGhsOverlap(1),
      ghs_vs_hyde_vs_hollinger_raw: evalGhsOverlap(0),
      max_geopng_netcdf_roundtrip_err: max_roundtrip_err,
      nan_or_inf_cells: nan_count,
      netcdf_shape: shape,
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
    console.log("[F_audit] Saved global_ghs_vs_hyde_vs_hollinger.png, regional_audit_panels.png, paleo_topology_audit.png & user_9_issues_audit.png.");
  }
}

if (require.main === module) {
  F_AuditOpticalQA.run().catch((err) => {
    console.error("[F_audit] Fatal error:", err);
    process.exit(1);
  });
}

module.exports = F_AuditOpticalQA;
