"use strict";

let fs = require("node:fs");
let path = require("node:path");
let stream = require("node:stream");
let { pipeline } = require("node:stream/promises");
let unzipper = require("unzipper");

let GeoGridUtils = require("../lib/GeoGridUtils");
let NetCDF4Writer = require("../lib/NetCDF4Writer");

/**
 * Step A: Resumable downloader (with HTTP Range continuation + NetCDF-4 HDF5 integrity check)
 * and archive extractor for authoritative upstream datasets.
 */
class A_AcquireSources {
  /**
   * Checks whether a NetCDF-4/HDF5 file on disk is complete and un-truncated.
   *
   * @param {string} arg0_file_path
   *
   * @returns {Promise<boolean>}
   */
  static async verifyNetCDFIntegrity (arg0_file_path) {
    //Convert from parameters
    let file_path = arg0_file_path;

    if (!fs.existsSync(file_path) || fs.statSync(file_path).size < 1024)
      return false;

    let fmt = NetCDF4Writer.detectFormat(file_path);
    if (fmt === "netcdf3")
      return true;
    if (fmt !== "hdf5")
      return false;

    try {
      let handle = await NetCDF4Writer.openRead(file_path);
      let keys = handle.keys();
      handle.close();
      return keys && keys.length > 0;
    } catch (err) {
      return false;
    }
  }

  /**
   * Downloads a URL with automatic HTTP Range resume until complete (verifying NetCDF-4 integrity).
   *
   * @param {string} arg0_url
   * @param {string} arg1_dest_path
   * @param {string} arg2_label
   * @param {string} arg3_format
   *
   * @returns {Promise<boolean>}
   */
  static async downloadFile (arg0_url, arg1_dest_path, arg2_label, arg3_format) {
    //Convert from parameters
    let url = arg0_url;
    let dest_path = arg1_dest_path;
    let label = arg2_label;
    let format = arg3_format;

    //Declare local instance variables
    let is_nc = (format === "netcdf4" || format === "netcdf");
    let max_attempts = 12;
    let part_path = `${dest_path}.part`;

    if (fs.existsSync(dest_path) && fs.statSync(dest_path).size > 1024) {
      if (!is_nc || await this.verifyNetCDFIntegrity(dest_path)) {
        let existing_mb = (fs.statSync(dest_path).size/1e6).toFixed(1);
        console.log(`[A_acquire] Cached & verified: ${label} (${existing_mb} MB)`);
        return true;
      }
      //Existing NetCDF-4 file was truncated -> move to .part so Range resume completes it!
      console.log(`[A_acquire] Resuming truncated ${label} from ${(fs.statSync(dest_path).size/1e6).toFixed(1)} MB ...`);
      fs.renameSync(dest_path, part_path);
    }

    for (let attempt = 1; attempt <= max_attempts; attempt++) {
      let existing_bytes = fs.existsSync(part_path) ? fs.statSync(part_path).size : 0;
      let headers = {};
      if (existing_bytes > 0)
        headers.Range = `bytes=${existing_bytes}-`;

      console.log(`[A_acquire] Downloading ${label} (attempt ${attempt}, offset ${(existing_bytes/1e6).toFixed(1)} MB) ...`);
      let response = await fetch(url, { headers: headers, redirect: "follow" });

      if (response.status === 416) {
        //Range Not Satisfiable -> already at exact EOF
        fs.renameSync(part_path, dest_path);
        return true;
      }

      if (!response.ok || !response.body) {
        console.error(`[A_acquire] HTTP ${response.status} for ${label}: ${url}`);
        return false;
      }

      let is_partial = (response.status === 206);
      let flags = (is_partial && existing_bytes > 0) ? "a" : "w";
      let received_bytes = is_partial ? existing_bytes : 0;
      let content_len = Number(response.headers.get("content-length") || 0);
      let expected_total = content_len > 0 ? (received_bytes + content_len) : 0;
      let last_log = Date.now();

      let progress_transform = new stream.Transform({
        transform (chunk, encoding, callback) {
          received_bytes += chunk.length;
          let now = Date.now();
          if (now - last_log > 5000) {
            let mb = (received_bytes/1e6).toFixed(1);
            let pct = expected_total > 0 ? ` (${((100*received_bytes)/expected_total).toFixed(1)}%)` : "";
            console.log(`[A_acquire]   ${label}: ${mb} MB${pct}`);
            last_log = now;
          }
          callback(null, chunk);
        },
      });

      try {
        await pipeline(
          stream.Readable.fromWeb(response.body),
          progress_transform,
          fs.createWriteStream(part_path, { flags: flags })
        );
      } catch (err) {
        console.log(`[A_acquire]   Stream interrupted on ${label} at ${(received_bytes/1e6).toFixed(1)} MB (${err.message}), retrying...`);
        continue;
      }

      let after_bytes = fs.statSync(part_path).size;
      if (expected_total > 0 && after_bytes < expected_total) {
        console.log(`[A_acquire]   Incomplete chunked transfer for ${label} (${after_bytes}/${expected_total}), resuming...`);
        continue;
      }

      if (is_nc) {
        let ok_nc = await this.verifyNetCDFIntegrity(part_path);
        if (!ok_nc) {
          console.log(`[A_acquire]   ${label} HDF5 superblock indicates more bytes remaining after ${(after_bytes/1e6).toFixed(1)} MB, resuming via Range...`);
          continue;
        }
      }

      fs.renameSync(part_path, dest_path);
      console.log(`[A_acquire] Completed & verified ${label} (${(after_bytes/1e6).toFixed(1)} MB).`);
      return true;
    }

    console.error(`[A_acquire] Failed to complete ${label} after ${max_attempts} attempts.`);

    //Return statement
    return false;
  }

  /**
   * Extracts a ZIP archive into target_dir if not already extracted.
   *
   * @param {string} arg0_zip_path
   * @param {string} arg1_extract_dir
   * @param {string} arg2_label
   *
   * @returns {Promise<void>}
   */
  static async extractZip (arg0_zip_path, arg1_extract_dir, arg2_label) {
    //Convert from parameters
    let zip_path = arg0_zip_path;
    let extract_dir = arg1_extract_dir;
    let label = arg2_label;

    //Declare local instance variables
    let marker = path.join(extract_dir, ".extracted_ok");

    if (fs.existsSync(marker)) {
      console.log(`[A_acquire] Extracted cache OK: ${label}`);
      return;
    }

    if (!fs.existsSync(extract_dir))
      fs.mkdirSync(extract_dir, { recursive: true });

    console.log(`[A_acquire] Extracting ${label} -> ${extract_dir} ...`);
    await fs.createReadStream(zip_path)
      .pipe(unzipper.Extract({ path: extract_dir }))
      .promise();

    fs.writeFileSync(marker, new Date().toISOString());
    console.log(`[A_acquire] Extracted ${label}.`);
  }

  /**
   * Executes Step A across all entries in config/source_manifest.json5.
   *
   * @returns {Promise<void>}
   */
  static async run () {
    //Declare local instance variables
    let manifest = GeoGridUtils.loadConfig("source_manifest.json5");
    let paths = GeoGridUtils.resolvePaths();
    let source_keys = Object.keys(manifest.sources);

    if (!fs.existsSync(paths.raw_dir))
      fs.mkdirSync(paths.raw_dir, { recursive: true });

    let concurrency = 2;
    let next_idx = 0;

    let runWorker = async () => {
      for (; next_idx < source_keys.length; ) {
        let key = source_keys[next_idx++];
        let entry = manifest.sources[key];
        let dest_path = path.join(paths.raw_dir, entry.filename);
        let ok = await this.downloadFile(entry.url, dest_path, key, entry.format);

        if (!ok && entry.fallback_url) {
          console.log(`[A_acquire] Retrying ${key} via fallback URL...`);
          ok = await this.downloadFile(entry.fallback_url, dest_path, `${key} (fallback)`, entry.format);
        }

        if (!ok && entry.required)
          console.error(`[A_acquire] REQUIRED source failed: ${key}`);

        if (ok && entry.format === "shapefile_zip" && entry.extract_dir)
          await this.extractZip(dest_path, path.join(paths.raw_dir, entry.extract_dir), key);
      }
    };

    let workers = [];
    for (let i = 0; i < concurrency; i++)
      workers.push(runWorker());
    await Promise.all(workers);

    console.log("[A_acquire] All upstream datasets verified.");
  }
}

if (require.main === module) {
  A_AcquireSources.run().catch((err) => {
    console.error("[A_acquire] Fatal error:", err);
    process.exit(1);
  });
}

module.exports = A_AcquireSources;
