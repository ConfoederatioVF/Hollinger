"use strict";

let fs = require("node:fs");
let path = require("node:path");

let sharp = require("sharp");

let GeoGridUtils = require("../core/GeoGridUtils");

/**
 * Renders contact sheets of zoomed regional crops in the video palette for visual verification.
 *
 * Usage: node tests/regional_preview.js <name> <lat_min> <lat_max> <lon_min> <lon_max> <year0> <year1> [scale]
 */
class RegionalPreview {
  /**
   * @returns {Promise<void>}
   */
  static async run () {
    //Declare local instance variables
    let argv = process.argv.slice(2);
    let cfg = GeoGridUtils.loadConfig("model_config.json5");
    let rcfg = GeoGridUtils.loadConfig("render_config.json5");
    let cell_area = GeoGridUtils.computeCellAreaByRow(cfg.grid);
    let paths = GeoGridUtils.resolvePaths("landarea");
    let name = argv[0];
    let lat_min = parseFloat(argv[1]);
    let lat_max = parseFloat(argv[2]);
    let lon_min = parseFloat(argv[3]);
    let lon_max = parseFloat(argv[4]);
    let y0 = parseInt(argv[5], 10);
    let y1 = parseInt(argv[6], 10);
    let scale = argv[7] ? parseInt(argv[7], 10) : 3;
    let gl = GeoGridUtils.parseHexRgb(rcfg.palette.glacier_full);
    let ln = GeoGridUtils.parseHexRgb(rcfg.palette.land_full);
    let sea = GeoGridUtils.parseHexRgb(rcfg.palette.sea);
    let r0 = Math.floor((90 - lat_max)*12);
    let r1 = Math.ceil((90 - lat_min)*12);
    let c0 = Math.floor((lon_min + 180)*12);
    let c1 = Math.ceil((lon_max + 180)*12);
    let w = c1 - c0;
    let h = r1 - r0;
    let years = GeoGridUtils.getTimestepYears(cfg).filter((y) => y >= y0 && y <= y1);
    let out_dir = path.join(paths.diagnostics_dir, "regional_previews", name);
    let tiles = [];

    if (!fs.existsSync(out_dir))
      fs.mkdirSync(out_dir, { recursive: true });

    for (let i = 0; i < years.length; i++) {
      let data = GeoGridUtils.decodeFloat32GeoPNG(path.join(paths.rasters_dir, `${cfg.encoding.geopng_prefix}${years[i]}.png`)).data;
      let rgb = Buffer.alloc(w*h*3);

      for (let y = 0; y < h; y++)
        for (let x = 0; x < w; x++) {
          let v = data[(r0 + y)*4320 + c0 + x];
          let f = Math.min(1, Math.abs(v)/cell_area[r0 + y]);
          let col = (v > 0) ? ln : ((v < 0) ? gl : sea);
          let o = (y*w + x)*3;
          for (let k = 0; k < 3; k++)
            rgb[o + k] = (v === 0) ? sea[k] : Math.round(sea[k] + (col[k] - sea[k])*f);
        }

      let label = Buffer.from(`<svg width="${w*scale}" height="${h*scale}" xmlns="http://www.w3.org/2000/svg"><rect x="0" y="0" width="150" height="30" fill="#000"/><text x="6" y="22" font-family="Arial" font-size="20" fill="#fff">${years[i]}</text></svg>`);
      let tile = await sharp(rgb, { raw: { channels: 3, height: h, width: w } })
        .resize(w*scale, h*scale, { kernel: "nearest" })
        .composite([{ input: label, left: 0, top: 0 }])
        .png()
        .toBuffer();
      tiles.push(tile);
    }

    //Write contact sheets of up to 6 tiles each (3 columns)
    let per_sheet = 6;
    for (let s = 0; s*per_sheet < tiles.length; s++) {
      let chunk = tiles.slice(s*per_sheet, (s + 1)*per_sheet);
      let cols = Math.min(3, chunk.length);
      let rows = Math.ceil(chunk.length/3);
      let composites = chunk.map((t, k) => ({ input: t, left: (k % 3)*(w*scale + 4), top: Math.floor(k/3)*(h*scale + 4) }));
      let out = path.join(out_dir, `sheet_${String(s).padStart(2, "0")}_${years[s*per_sheet]}.png`);
      await sharp({ create: { background: "#000", channels: 3, height: rows*(h*scale + 4), width: cols*(w*scale + 4) } })
        .composite(composites)
        .png()
        .toFile(out);
      console.log(out);
    }
  }
}

RegionalPreview.run().catch((err) => {
  console.error(err);
  process.exit(1);
});
