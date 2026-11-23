"use strict";

let fs = require("node:fs");
let path = require("node:path");

let netcdfjs = require("netcdfjs");

/**
 * Pure Node.js NetCDF-3 / NetCDF-4 (HDF5) reader and CF-1.8 compressed NetCDF-4 writer powered by h5wasm and netcdfjs.
 */
class NetCDF4Writer {
  static _h5_module = null;
  static _mounted_dirs = new Set();

  /**
   * Lazily initialises the h5wasm/node WebAssembly module.
   *
   * @returns {Promise<Object>}
   */
  static async getH5 () {
    //Guard clause
    if (this._h5_module)
      return this._h5_module;

    //Declare local instance variables
    let h5wasm = await import("h5wasm/node");
    await h5wasm.ready;
    this._h5_module = h5wasm;

    //Return statement
    return h5wasm;
  }

  /**
   * Checks whether a file begins with the HDF5 magic signature (NetCDF-4) or classic CDF signature (NetCDF-3).
   *
   * @param {string} arg0_file_path
   *
   * @returns {string} "hdf5" | "netcdf3" | "unknown"
   */
  static detectFormat (arg0_file_path) {
    //Convert from parameters
    let file_path = arg0_file_path;

    //Declare local instance variables
    let buf = Buffer.alloc(8);
    let fd = fs.openSync(file_path, "r");
    fs.readSync(fd, buf, 0, 8, 0);
    fs.closeSync(fd);

    if (buf[0] === 0x89 && buf[1] === 0x48 && buf[2] === 0x44 && buf[3] === 0x46)
      return "hdf5";
    if (buf[0] === 0x43 && buf[1] === 0x44 && buf[2] === 0x46)
      return "netcdf3";

    //Return statement
    return "unknown";
  }

  /**
   * Opens a NetCDF-3 or NetCDF-4/HDF5 file and returns a uniform handle for reading variables and slices.
   *
   * @param {string} arg0_file_path
   *
   * @returns {Promise<Object>}
   */
  static async openRead (arg0_file_path) {
    //Convert from parameters
    let file_path = path.resolve(arg0_file_path);

    //Declare local instance variables
    let format = this.detectFormat(file_path);

    if (format === "hdf5") {
      let h5 = await this.getH5();
      let h5_file = new h5.File(file_path, "r");

      //Return statement
      return {
        close: () => h5_file.close(),
        format: "hdf5",
        getShape: (arg0_name) => h5_file.get(arg0_name).shape,
        getSlice: (arg0_name, arg1_ranges) => h5_file.get(arg0_name).slice(arg1_ranges),
        getValue: (arg0_name) => h5_file.get(arg0_name).value,
        keys: () => h5_file.keys(),
        raw: h5_file,
      };
    }

    if (format === "netcdf3") {
      let buf = fs.readFileSync(file_path);
      let reader = new netcdfjs.NetCDFReader(buf);
      let var_cache = new Map();
      let getVarCached = (arg0_name) => {
        if (!var_cache.has(arg0_name))
          var_cache.set(arg0_name, reader.getDataVariable(arg0_name));
        return var_cache.get(arg0_name);
      };

      //Return statement
      return {
        close: () => { var_cache.clear(); },
        format: "netcdf3",
        getShape: (arg0_name) => {
          let v = reader.variables.find((item) => item.name === arg0_name);
          if (!v) return [];
          return v.dimensions.map((dim_idx) => {
            let sz = reader.dimensions[dim_idx].size;
            return (sz === 0) ? reader.header.recordDimension.length : sz;
          });
        },
        getSlice: (arg0_name, arg1_ranges) => {
          let all_vals = getVarCached(arg0_name);
          let v = reader.variables.find((item) => item.name === arg0_name);
          let shape = v.dimensions.map((dim_idx) => {
            let sz = reader.dimensions[dim_idx].size;
            return (sz === 0) ? reader.header.recordDimension.length : sz;
          });
          if (shape.length === 3 && arg1_ranges.length === 3) {
            let t0 = arg1_ranges[0][0];
            //netcdfjs returns record variables as Array<Array|TypedArray> (one entry per record) or flat array
            if (Array.isArray(all_vals) && all_vals.length === shape[0] && (Array.isArray(all_vals[t0]) || ArrayBuffer.isView(all_vals[t0])))
              return Float32Array.from(all_vals[t0]);
            let slice_size = shape[1]*shape[2];
            let start = t0*slice_size;
            let out = new Float32Array(slice_size);
            for (let i = 0; i < slice_size; i++)
              out[i] = all_vals[start + i];
            return out;
          }
          return Float32Array.from(all_vals);
        },
        getValue: (arg0_name) => Float32Array.from(getVarCached(arg0_name)),
        keys: () => reader.variables.map((v) => v.name),
        raw: reader,
      };
    }

    console.error(`[NetCDF4Writer] Unsupported file format for: ${file_path}`);

    //Return statement
    return null;
  }

  /**
   * Writes the CF-1.8 compliant compressed NetCDF-4 file land_area_Hollinger.nc directly on disk via h5wasm/node.
   *
   * @param {string} arg0_output_path
   * @param {Object} arg1_spec
   *  @param {Float64Array|Float32Array} arg1_spec.cell_area_km2
   *  @param {Array<number>} arg1_spec.chunk_shape
   *  @param {number} arg1_spec.compression_level
   *  @param {Float32Array} arg1_spec.lat
   *  @param {Float32Array} arg1_spec.lon
   *  @param {Array<number>} arg1_spec.years
   * @param {Function} arg2_get_slice_fn - Async or sync callback(year, t_idx) -> Float32Array(2160*4320)
   *
   * @returns {Promise<void>}
   */
  static async writeHollingerNetCDF4 (arg0_output_path, arg1_spec, arg2_get_slice_fn) {
    //Convert from parameters
    let output_path = path.resolve(arg0_output_path);
    let spec = arg1_spec;
    let get_slice_fn = arg2_get_slice_fn;

    //Declare local instance variables
    let cell_area_f32 = Float32Array.from(spec.cell_area_km2);
    let chunk_shape = spec.chunk_shape ? spec.chunk_shape : [1, 216, 432];
    let comp_level = (spec.compression_level !== undefined) ? spec.compression_level : 4;
    let dir_path = path.dirname(output_path);
    let h5 = await this.getH5();
    let height = spec.lat.length;
    let num_times = spec.years.length;
    let width = spec.lon.length;
    let years_i32 = Int32Array.from(spec.years);

    if (!fs.existsSync(dir_path))
      fs.mkdirSync(dir_path, { recursive: true });
    if (fs.existsSync(output_path))
      fs.unlinkSync(output_path);

    let file = new h5.File(output_path, "w");

    //1. Global CF-1.8 metadata attributes
    file.create_attribute("Conventions", "CF-1.8");
    file.create_attribute("title", "land_area_Hollinger: Dynamic Global Ice-Free Land Area and Glaciation (10000 BC - 2025 AD)");
    file.create_attribute("institution", "Project 1509 - SVEA / Hollinger Geomorphology Pipeline");
    file.create_attribute("source", "OSM Land Polygons, ETOPO 2022 (Bedrock & Surface), PMIP4 GLAC-1D, HydroLAKES v1.0, Global Dam Watch v1.0, Natural Earth Glaciers & Antarctic Ice Shelves, Curated Historical Reclamation & Paleolake Registries");
    file.create_attribute("comment", "Signed float32 km2 per 5-arcminute WGS84 cell. Positive (>0) = ice-free terrestrial land; Negative (<0) = glaciated area (grounded ice sheets, mountain glaciers, and permanent year-round floating ice shelves); Zero (0) = sea, seasonal sea ice, lakes, and impounded reservoirs. Precedence: Land > Glaciated > Sea. Year 0 corresponds to 1 AD (matching HYDE).");
    file.create_attribute("history", `Created ${new Date().toISOString()} by land_area_Hollinger Node.js pipeline`);

    //2. Coordinate variables
    let time_ds = file.create_dataset({
      data: years_i32,
      dtype: "<i4",
      name: "time",
      shape: [num_times],
    });
    time_ds.create_attribute("standard_name", "time");
    time_ds.create_attribute("long_name", "Signed calendar year (negative = BC, 0 = 1 AD matching HYDE, positive = AD)");
    time_ds.create_attribute("units", "years");
    time_ds.create_attribute("axis", "T");

    let lat_ds = file.create_dataset({
      data: spec.lat,
      dtype: "<f4",
      name: "lat",
      shape: [height],
    });
    lat_ds.create_attribute("standard_name", "latitude");
    lat_ds.create_attribute("long_name", "Latitude of gridcell centre");
    lat_ds.create_attribute("units", "degrees_north");
    lat_ds.create_attribute("axis", "Y");

    let lon_ds = file.create_dataset({
      data: spec.lon,
      dtype: "<f4",
      name: "lon",
      shape: [width],
    });
    lon_ds.create_attribute("standard_name", "longitude");
    lon_ds.create_attribute("long_name", "Longitude of gridcell centre");
    lon_ds.create_attribute("units", "degrees_east");
    lon_ds.create_attribute("axis", "X");

    let area_ds = file.create_dataset({
      data: cell_area_f32,
      dtype: "<f4",
      name: "cell_area_km2",
      shape: [height],
    });
    area_ds.create_attribute("long_name", "Exact WGS84 ellipsoidal area of each 5-arcminute gridcell by latitude row");
    area_ds.create_attribute("units", "km2");

    //3. Primary 3D compressed float32 land_area dataset (chunked per timestep for streaming writes)
    let first_slice = await get_slice_fn(spec.years[0], 0);
    let land_ds = file.create_dataset({
      chunks: chunk_shape,
      compression: comp_level,
      data: first_slice,
      dtype: "<f4",
      maxshape: [num_times, height, width],
      name: "land_area",
      shape: [1, height, width],
    });
    if (num_times > 1)
      land_ds.resize([num_times, height, width]);
    land_ds.create_attribute("long_name", "Signed terrestrial land (>0) and glaciated (<0) area per gridcell; sea/water = 0");
    land_ds.create_attribute("units", "km2");
    land_ds.create_attribute("coordinates", "time lat lon");
    land_ds.create_attribute("precedence", "land > glaciated > sea");

    //4. Attach HDF5 / NetCDF-4 Dimension Scales so netCDF4-python, xarray, Panoply, and QGIS recognise named dimensions
    try {
      h5.Module.set_scale(file.file_id, "time", "time");
      h5.Module.set_scale(file.file_id, "lat", "lat");
      h5.Module.set_scale(file.file_id, "lon", "lon");
      h5.Module.attach_scale(file.file_id, "cell_area_km2", "lat", 0);
      h5.Module.attach_scale(file.file_id, "land_area", "time", 0);
      h5.Module.attach_scale(file.file_id, "land_area", "lat", 1);
      h5.Module.attach_scale(file.file_id, "land_area", "lon", 2);
    } catch (err) {
      console.error("[NetCDF4Writer] Warning: could not attach HDF5 dimension scales:", err.message);
    }

    //5. Stream remaining timestep slices directly into compressed HDF5 chunks
    for (let i = 1; i < num_times; i++) {
      let year = spec.years[i];
      let slice = await get_slice_fn(year, i);
      land_ds.write_slice([[i, i + 1], [0, height], [0, width]], slice);
    }

    file.close();
  }
}

module.exports = NetCDF4Writer;
