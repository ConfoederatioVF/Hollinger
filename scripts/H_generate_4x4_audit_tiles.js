"use strict";

let fs = require("node:fs");
let path = require("node:path");

let sharp = require("sharp");

let GeoGridUtils = require("../lib/GeoGridUtils");

/**
 * Step H: Splits the global 4320x2160 Hollinger rasters into a 4x4 grid of 16 native-resolution
 * (1080x540, 90deg lon x 45deg lat) audit tiles across canonical historical epochs so the
 * adversarial Coastline Auditor can inspect every coastline on Earth at 1:1 5-arcminute clarity.
 */
class H_Generate4x4AuditTiles {
  /**
   * Renders a 4x4 tile set for a single year into diagnostics/audit_4x4/<year>/.
   *
   * @param {number} arg0_year
   * @param {Float64Array} arg1_cell_area
   * @param {Object} arg2_paths
   * @param {string} arg3_prefix
   *
   * @returns {Promise<Array<Object>>} Manifest of generated non-empty tiles
   */
  static async renderYearTiles (arg0_year, arg1_cell_area, arg2_paths, arg3_prefix) {
    //Convert from parameters
    let year = arg0_year;
    let cell_area = arg1_cell_area;
    let paths = arg2_paths;
    let prefix = arg3_prefix;

    //Declare local instance variables
    let png_path = path.join(paths.rasters_dir, `${prefix}${year}.png`);
    let out_dir = path.join(paths.diagnostics_dir, "audit_4x4", String(year));
    let data = GeoGridUtils.decodeFloat32GeoPNG(png_path).data;
    let tile_h = 540;
    let tile_w = 1080;
    let manifest = [];

    if (!fs.existsSync(out_dir))
      fs.mkdirSync(out_dir, { recursive: true });

    let display_year = (year === 0) ? "1 AD (0)" : ((year < 0) ? `${Math.abs(year)} BC` : `${year} AD`);

    for (let tr = 0; tr < 4; tr++) {
      let r0 = tr*tile_h;
      let lat_max = 90 - tr*45;
      let lat_min = lat_max - 45;

      for (let tc = 0; tc < 4; tc++) {
        let c0 = tc*tile_w;
        let lon_min = -180 + tc*90;
        let lon_max = lon_min + 90;
        let header_h = 30;
        let out_h = tile_h + header_h;
        let rgb = Buffer.allocUnsafe(tile_w*out_h*3);
        let land_px = 0;
        let ice_px = 0;

        //Fill dedicated 30px top header strip (#1C2127) so ZERO map pixels are ever obscured (fixes Great Lakes!)
        for (let hy = 0; hy < header_h; hy++) {
          let h_off = hy*tile_w*3;
          for (let hx = 0; hx < tile_w; hx++) {
            let o3 = h_off + hx*3;
            rgb[o3] = 28; rgb[o3 + 1] = 33; rgb[o3 + 2] = 39;
          }
        }

        for (let y = 0; y < tile_h; y++) {
          let gr = r0 + y;
          let lat_deg = 90.0 - (gr + 0.5)*(180.0/2160.0);
          let inv_a = 1.0/cell_area[gr];
          let g_off = gr*4320;
          let t_off = (y + header_h)*tile_w;

          for (let x = 0; x < tile_w; x++) {
            let gc = c0 + x;
            let v = data[g_off + gc];
            let o3 = (t_off + x)*3;

            if (v > 0.0) {
              land_px++;
              let f = Math.min(1.0, v*inv_a);
              if (f >= 0.995) {
                rgb[o3] = 164; rgb[o3 + 1] = 198; rgb[o3 + 2] = 134;
              } else {
                let n_l = data[g_off + ((gc === 0) ? 4319 : (gc - 1))];
                let n_r = data[g_off + ((gc === 4319) ? 0 : (gc + 1))];
                let n_u = (gr > 0) ? data[(gr - 1)*4320 + gc] : v;
                let n_d = (gr < 2159) ? data[(gr + 1)*4320 + gc] : v;
                let adj_sea = (n_l === 0.0 || n_r === 0.0 || n_u === 0.0 || n_d === 0.0);
                let adj_ice = (n_l < 0.0 || n_r < 0.0 || n_u < 0.0 || n_d < 0.0);
                let isolated_atoll = (n_l === 0.0 && n_r === 0.0 && n_u === 0.0 && n_d === 0.0);
                if (!adj_sea && adj_ice) {
                  rgb[o3] = Math.round(228 + (164 - 228)*f);
                  rgb[o3 + 1] = Math.round(232 + (198 - 232)*f);
                  rgb[o3 + 2] = Math.round(236 + (134 - 236)*f);
                } else {
                  //R2-D04: Crisp shoreline smoothstep [0.42..0.62] renders majority-water boreal lakes/reservoirs
                  //(Smallwood, Manicouagan ring, Caniapiscau, Robert-Bourassa, Mistassini) as sharp blue water
                  //while preserving isolated single-pixel oceanic micro-islands & atolls!
                  let f_vis = f;
                  if (!isolated_atoll) {
                    if (f <= 0.42) {
                      f_vis = 0.0;
                    } else if (f >= 0.62) {
                      f_vis = 1.0;
                    } else {
                      let t = (f - 0.42)/0.20;
                      f_vis = t*t*(3.0 - 2.0*t);
                    }
                  }
                  rgb[o3] = Math.round(162 + (164 - 162)*f_vis);
                  rgb[o3 + 1] = Math.round(210 + (198 - 210)*f_vis);
                  rgb[o3 + 2] = Math.round(238 + (134 - 238)*f_vis);
                }
              }
            } else if (v < 0.0) {
              ice_px++;
              let f = Math.min(1.0, (-v)*inv_a);
              rgb[o3] = Math.round(162 + (228 - 162)*f);
              rgb[o3 + 1] = Math.round(210 + (232 - 210)*f);
              rgb[o3 + 2] = Math.round(238 + (236 - 238)*f);
            } else {
              rgb[o3] = 162; rgb[o3 + 1] = 210; rgb[o3 + 2] = 238;
            }
          }
        }

        //Draw subtle 15-degree reference graticule lines + unobstructed top header label
        let svg_parts = [];
        for (let g = 1; g < 6; g++) {
          let gx = Math.round((g/6)*tile_w);
          let glon = lon_min + g*15;
          svg_parts.push(`<line x1="${gx}" y1="${header_h}" x2="${gx}" y2="${out_h}" stroke="#000000" stroke-opacity="0.14" stroke-width="1" />`);
          svg_parts.push(`<text x="${gx + 4}" y="${out_h - 6}" font-family="Arial" font-size="12" fill="#1C2127" fill-opacity="0.7">${glon}E</text>`);
        }
        for (let g = 1; g < 3; g++) {
          let gy = header_h + Math.round((g/3)*tile_h);
          let glat = lat_max - g*15;
          svg_parts.push(`<line x1="0" y1="${gy}" x2="${tile_w}" y2="${gy}" stroke="#000000" stroke-opacity="0.14" stroke-width="1" />`);
          svg_parts.push(`<text x="6" y="${gy - 4}" font-family="Arial" font-size="12" fill="#1C2127" fill-opacity="0.7">${glat}N</text>`);
        }

        let label = `Tile r${tr}_c${tc} | ${display_year} | lat ${lat_min}..${lat_max}, lon ${lon_min}..${lon_max}`;
        svg_parts.push(`<text x="12" y="21" font-family="Arial" font-size="15" fill="#F2F4F7">${label}</text>`);

        let svg = Buffer.from(`<svg width="${tile_w}" height="${out_h}" xmlns="http://www.w3.org/2000/svg">${svg_parts.join("")}</svg>`);
        let out_file = path.join(out_dir, `tile_r${tr}_c${tc}.png`);

        await sharp(rgb, { raw: { channels: 3, height: out_h, width: tile_w } })
          .composite([{ input: svg, left: 0, top: 0 }])
          .png({ compressionLevel: 6 })
          .toFile(out_file);

        manifest.push({
          file: out_file,
          ice_px: ice_px,
          land_px: land_px,
          lat_max: lat_max,
          lat_min: lat_min,
          lon_max: lon_max,
          lon_min: lon_min,
          tc: tc,
          tr: tr,
          year: year,
        });
      }
    }

    //Return statement
    return manifest;
  }

  /**
   * Executes Step H across key historical audit epochs.
   *
   * @returns {Promise<void>}
   */
  static async run () {
    //Declare local instance variables
    let cfg = GeoGridUtils.loadConfig("model_config.json5");
    let cell_area = GeoGridUtils.computeCellAreaByRow(cfg.grid);
    let paths = GeoGridUtils.resolvePaths();
    let years_arg = process.argv.find((a) => a.startsWith("--years="));
    let years = years_arg
      ? years_arg.slice(8).split(",").map((s) => parseInt(s.trim(), 10))
      : [-10000, -6000, -4000, 0, 1500, 1900, 2025];
    let all_manifest = [];

    for (let i = 0; i < years.length; i++) {
      let y = years[i];
      let m = await this.renderYearTiles(y, cell_area, paths, cfg.encoding.geopng_prefix);
      all_manifest.push(...m);
      console.log(`[H_audit4x4] Rendered 16 4x4 tiles for year ${y}.`);
    }

    let manifest_path = path.join(paths.diagnostics_dir, "audit_4x4", "manifest.json");
    fs.writeFileSync(manifest_path, JSON.stringify(all_manifest, null, 2));
    console.log(`[H_audit4x4] Completed ${all_manifest.length} tiles across ${years.length} epochs -> ${manifest_path}`);
  }
}

if (require.main === module) {
  H_Generate4x4AuditTiles.run().catch((err) => {
    console.error("[H_audit4x4] Fatal error:", err);
    process.exit(1);
  });
}

module.exports = H_Generate4x4AuditTiles;
