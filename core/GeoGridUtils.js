"use strict";

let fs = require("node:fs");
let path = require("node:path");

let JSON5 = require("json5");
let pngjs = require("pngjs");
let sharp = require("sharp");

/**
 * Core WGS84 grid, signed float32 & uint32 GeoPNG I/O, timestep schedule, and scoped config utilities for landuse_Hollinger.
 */
class GeoGridUtils {
  /**
   * Applies strict signed precedence (Land > Glaciated > Sea) to produce float32 km2 grid values.
   *
   * @param {Float32Array} arg0_f_land - Ice-free terrestrial land fraction [0, 1].
   * @param {Float32Array} arg1_f_ice - Glaciated + permanent floating ice-shelf fraction [0, 1].
   * @param {Float64Array|Float32Array} arg2_cell_area_by_row - WGS84 cell area (km2) per row.
   * @param {Object} [arg3_options]
   *  @param {number} [arg3_options.epsilon=0.00025]
   *  @param {number} [arg3_options.height=2160]
   *  @param {number} [arg3_options.width=4320]
   *
   * @returns {Float32Array}
   */
  static applySignedPrecedence (arg0_f_land, arg1_f_ice, arg2_cell_area_by_row, arg3_options) {
    //Convert from parameters
    let f_land = arg0_f_land;
    let f_ice = arg1_f_ice;
    let cell_area_by_row = arg2_cell_area_by_row;
    let options = arg3_options ? arg3_options : {};

    //Initialise options
    if (options.epsilon === undefined) options.epsilon = 0.00025;
    if (options.height === undefined) options.height = 2160;
    if (options.width === undefined) options.width = 4320;

    //Declare local instance variables
    let eps = options.epsilon;
    let height = options.height;
    let out_array = new Float32Array(options.width*height);
    let width = options.width;

    for (let y = 0; y < height; y++) {
      let row_area = cell_area_by_row[y];
      let row_offset = y*width;

      for (let x = 0; x < width; x++) {
        let idx = row_offset + x;
        let land_frac = f_land[idx];

        if (land_frac > eps) {
          let clamped_land = Math.min(1.0, land_frac);
          out_array[idx] = clamped_land*row_area;
        } else {
          let ice_frac = f_ice[idx];
          if (ice_frac > eps) {
            let clamped_ice = Math.min(1.0, ice_frac);
            out_array[idx] = -clamped_ice*row_area;
          } else {
            out_array[idx] = 0.0;
          }
        }
      }
    }

    //Return statement
    return out_array;
  }

  /**
   * Computes exact WGS84 ellipsoidal cell area (km2) for each latitude row (north-to-south).
   *
   * @param {Object} [arg0_grid_config]
   *
   * @returns {Float64Array}
   */
  static computeCellAreaByRow (arg0_grid_config) {
    //Convert from parameters
    let grid_config = arg0_grid_config ? arg0_grid_config : {};

    //Initialise options
    let a = (grid_config.wgs84_semi_major_m !== undefined) ? grid_config.wgs84_semi_major_m : 6378137.0;
    let f = (grid_config.wgs84_flattening !== undefined) ? grid_config.wgs84_flattening : (1.0/298.257223563);
    let height = (grid_config.height !== undefined) ? grid_config.height : 2160;
    let width = (grid_config.width !== undefined) ? grid_config.width : 4320;

    //Declare local instance variables
    let area_by_row = new Float64Array(height);
    let b = a*(1.0 - f);
    let dlat_rad = Math.PI/height;
    let dlon_rad = (2.0*Math.PI)/width;
    let e2 = 1.0 - (b*b)/(a*a);
    let e = Math.sqrt(e2);
    let factor = (b*b*dlon_rad)/(2.0*1e6);

    //Exact authalic integral Q(phi) = sin(phi)/(1 - e^2 sin^2(phi)) + (1/(2e))*ln((1 + e sin(phi))/(1 - e sin(phi)))
    let computeQ = (arg0_phi) => {
      let s = Math.sin(arg0_phi);
      let es = e*s;
      return (s/(1.0 - es*es)) + (Math.log((1.0 + es)/(1.0 - es))/(2.0*e));
    };

    for (let y = 0; y < height; y++) {
      let phi_north = (Math.PI/2.0) - y*dlat_rad;
      let phi_south = (Math.PI/2.0) - (y + 1)*dlat_rad;
      area_by_row[y] = Math.abs(factor*(computeQ(phi_north) - computeQ(phi_south)));
    }

    //Return statement
    return area_by_row;
  }

  /**
   * Decodes a big-endian IEEE 754 float32 RGBA GeoPNG into a Float32Array.
   *
   * @param {string} arg0_file_path
   *
   * @returns {{ data: Float32Array, height: number, width: number }}
   */
  static decodeFloat32GeoPNG (arg0_file_path) {
    //Convert from parameters
    let file_path = arg0_file_path;

    //Declare local instance variables
    let png = pngjs.PNG.sync.read(fs.readFileSync(file_path));
    let raw_buf = png.data;
    let total_pixels = png.width*png.height;
    let out_data = new Float32Array(total_pixels);
    let view = new DataView(raw_buf.buffer, raw_buf.byteOffset, raw_buf.byteLength);

    for (let i = 0; i < total_pixels; i++)
      out_data[i] = view.getFloat32(i*4, false);

    //Return statement
    return {
      data: out_data,
      height: png.height,
      width: png.width,
    };
  }

  /**
   * Decodes a big-endian uint32 RGBA GeoPNG into a Uint32Array (matching SVEA GeoPNG int32 format).
   *
   * @param {string} arg0_file_path
   *
   * @returns {{ data: Uint32Array, height: number, width: number }}
   */
  static decodeUint32GeoPNG (arg0_file_path) {
    //Convert from parameters
    let file_path = arg0_file_path;

    //Declare local instance variables
    let png = pngjs.PNG.sync.read(fs.readFileSync(file_path));
    let raw_buf = png.data;
    let total_pixels = png.width*png.height;
    let out_data = new Uint32Array(total_pixels);
    let view = new DataView(raw_buf.buffer, raw_buf.byteOffset, raw_buf.byteLength);

    for (let i = 0; i < total_pixels; i++)
      out_data[i] = view.getUint32(i*4, false);

    //Return statement
    return {
      data: out_data,
      height: png.height,
      width: png.width,
    };
  }

  /**
   * Encodes a Float32Array as a big-endian IEEE 754 float32 RGBA GeoPNG using sharp.
   *
   * @param {string} arg0_file_path
   * @param {Float32Array} arg1_data
   * @param {Object} [arg2_options]
   *  @param {number} [arg2_options.compression_level=6]
   *  @param {number} [arg2_options.height=2160]
   *  @param {number} [arg2_options.width=4320]
   *
   * @returns {Promise<void>}
   */
  static async encodeFloat32GeoPNG (arg0_file_path, arg1_data, arg2_options) {
    //Convert from parameters
    let file_path = arg0_file_path;
    let data = arg1_data;
    let options = arg2_options ? arg2_options : {};

    //Initialise options
    if (options.compression_level === undefined) options.compression_level = 6;
    if (options.height === undefined) options.height = 2160;
    if (options.width === undefined) options.width = 4320;

    //Declare local instance variables
    let dir_path = path.dirname(file_path);
    let height = options.height;
    let total_pixels = options.width*height;
    let width = options.width;

    if (!fs.existsSync(dir_path))
      fs.mkdirSync(dir_path, { recursive: true });

    let rgba_buf = Buffer.allocUnsafe(total_pixels*4);
    for (let i = 0; i < total_pixels; i++)
      rgba_buf.writeFloatBE(data[i], i*4);

    await sharp(rgba_buf, {
      raw: {
        channels: 4,
        height: height,
        width: width,
      },
    })
      .png({
        compressionLevel: options.compression_level,
        adaptiveFiltering: false,
      })
      .toFile(file_path);
  }

  /**
   * Encodes a Uint32Array/Int32Array as a big-endian uint32 RGBA GeoPNG using sharp (matching SVEA GeoPNG int32 format).
   *
   * @param {string} arg0_file_path
   * @param {Uint32Array|Int32Array} arg1_data
   * @param {Object} [arg2_options]
   *  @param {number} [arg2_options.compression_level=6]
   *  @param {number} [arg2_options.height=2160]
   *  @param {number} [arg2_options.width=4320]
   *
   * @returns {Promise<void>}
   */
  static async encodeUint32GeoPNG (arg0_file_path, arg1_data, arg2_options) {
    //Convert from parameters
    let file_path = arg0_file_path;
    let data = arg1_data;
    let options = arg2_options ? arg2_options : {};

    //Initialise options
    if (options.compression_level === undefined) options.compression_level = 6;
    if (options.height === undefined) options.height = 2160;
    if (options.width === undefined) options.width = 4320;

    //Declare local instance variables
    let dir_path = path.dirname(file_path);
    let height = options.height;
    let total_pixels = options.width*height;
    let width = options.width;

    if (!fs.existsSync(dir_path))
      fs.mkdirSync(dir_path, { recursive: true });

    let rgba_buf = Buffer.allocUnsafe(total_pixels*4);
    for (let i = 0; i < total_pixels; i++)
      rgba_buf.writeUInt32BE(data[i] >>> 0, i*4);

    await sharp(rgba_buf, {
      raw: {
        channels: 4,
        height: height,
        width: width,
      },
    })
      .png({
        compressionLevel: options.compression_level,
        adaptiveFiltering: false,
      })
      .toFile(file_path);
  }

  /**
   * Builds the latitude and longitude centre arrays for the 5-arcminute WGS84 grid.
   *
   * @param {Object} [arg0_grid_config]
   *
   * @returns {{ lat: Float32Array, lon: Float32Array }}
   */
  static getGridCoordinates (arg0_grid_config) {
    //Convert from parameters
    let grid_config = arg0_grid_config ? arg0_grid_config : {};

    //Declare local instance variables
    let height = (grid_config.height !== undefined) ? grid_config.height : 2160;
    let lat = new Float32Array(height);
    let lat_step = 180.0/height;
    let width = (grid_config.width !== undefined) ? grid_config.width : 4320;
    let lon = new Float32Array(width);
    let lon_step = 360.0/width;

    for (let y = 0; y < height; y++)
      lat[y] = 90.0 - (y + 0.5)*lat_step;
    for (let x = 0; x < width; x++)
      lon[x] = -180.0 + (x + 0.5)*lon_step;

    //Return statement
    return { lat: lat, lon: lon };
  }

  /**
   * Builds the exact ordered array of 218 signed timestep years (-10000 to 2025).
   *
   * @param {Object} arg0_model_config
   *
   * @returns {Array<number>}
   */
  static getTimestepYears (arg0_model_config) {
    //Convert from parameters
    let model_config = arg0_model_config;

    //Declare local instance variables
    let t = model_config.timesteps;
    let years = [];

    for (let y = t.bce_start; y <= t.bce_end; y += t.bce_step)
      years.push(y);
    for (let y = t.ce_century_start; y <= t.ce_century_end; y += t.ce_century_step)
      years.push(y);
    for (let y = t.ce_decade_start; y <= t.ce_decade_end; y += t.ce_decade_step)
      years.push(y);
    for (let y = t.ce_annual_start; y <= t.ce_annual_end; y += t.ce_annual_step)
      years.push(y);

    //Return statement
    return years;
  }

  /**
   * Loads a scoped JSON5 configuration file from landuse_Hollinger/common/.
   *
   * @param {string} arg0_config_filename
   *
   * @returns {Object}
   */
  static loadConfig (arg0_config_filename) {
    //Convert from parameters
    let config_filename = arg0_config_filename;

    //Declare local instance variables
    let paths = this.resolvePaths();
    let full_path = path.join(paths.common_dir, config_filename);
    let raw_text = fs.readFileSync(full_path, "utf8");

    //Return statement
    return JSON5.parse(raw_text);
  }

  /**
   * Parses a hex colour string (#RRGGBB) into an [r, g, b] array.
   *
   * @param {string} arg0_hex
   *
   * @returns {Array<number>}
   */
  static parseHexRgb (arg0_hex) {
    //Convert from parameters
    let hex = String(arg0_hex).replace(/^#/, "").trim();

    //Declare local instance variables
    let num = parseInt(hex, 16);

    //Return statement
    return [(num >> 16) & 255, (num >> 8) & 255, num & 255];
  }

  /**
   * Reads a raw binary Float32Array file from disk.
   *
   * @param {string} arg0_file_path
   *
   * @returns {Float32Array}
   */
  static readFloat32Bin (arg0_file_path) {
    //Convert from parameters
    let file_path = arg0_file_path;

    //Declare local instance variables
    let buf = fs.readFileSync(file_path);
    let out = new Float32Array(buf.byteLength/4);
    new Uint8Array(out.buffer).set(buf);

    //Return statement
    return out;
  }

  /**
   * Resolves all canonical scoped directories inside histmap/1.data_raw/landuse_Hollinger/ for an indicator.
   *
   * @param {string} [arg0_indicator="landarea"]
   *
   * @returns {Object}
   */
  static resolvePaths (arg0_indicator) {
    //Convert from parameters
    let indicator = arg0_indicator ? arg0_indicator : "landarea";

    //Declare local instance variables
    let base_dir = path.resolve(__dirname, "..");
    let common_dir = path.join(base_dir, "common");
    let core_dir = path.join(base_dir, "core");
    let indicator_dir = path.join(base_dir, indicator);
    let tests_dir = path.join(base_dir, "tests");

    //Return statement
    return {
      base_dir: base_dir,
      common_dir: common_dir,
      config_dir: common_dir,
      core_dir: core_dir,
      curated_dir: path.join(indicator_dir, "1.raw_data", "curated"),
      diagnostics_dir: path.join(tests_dir, "diagnostics", indicator),
      exports_dir: path.join(base_dir, "exports"),
      indicator_dir: indicator_dir,
      intermediate_dir: path.join(indicator_dir, "2.intermediates"),
      lib_dir: core_dir,
      nc_path: path.join(indicator_dir, `${indicator}.nc`),
      pipeline_dir: path.join(indicator_dir, "pipeline"),
      rasters_dir: path.join(indicator_dir, "3.rasters"),
      raw_dir: path.join(indicator_dir, "1.raw_data"),
      renders_dir: path.join(base_dir, "exports", indicator),
      scripts_dir: path.join(indicator_dir, "pipeline"),
      tests_dir: tests_dir,
    };
  }

  /**
   * Writes a raw binary Float32Array file to disk.
   *
   * @param {string} arg0_file_path
   * @param {Float32Array} arg1_data
   */
  static writeFloat32Bin (arg0_file_path, arg1_data) {
    //Convert from parameters
    let file_path = arg0_file_path;
    let data = arg1_data;

    //Declare local instance variables
    let dir_path = path.dirname(file_path);

    if (!fs.existsSync(dir_path))
      fs.mkdirSync(dir_path, { recursive: true });

    fs.writeFileSync(file_path, Buffer.from(data.buffer, data.byteOffset, data.byteLength));
  }
}

module.exports = GeoGridUtils;
