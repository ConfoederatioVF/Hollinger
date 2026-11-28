/**
 * @fileoverview Parallel 2-pass libx264 compressor for all Hollinger 4K MP4 exports under exports/,
 * targeting <= 20 MB (decimal 20 000 000 bytes) with the _compressed suffix.
 */

"use strict";

let fs = require("node:fs");
let os = require("node:os");
let path = require("node:path");
let { spawn, spawnSync } = require("node:child_process");
let ffmpeg_path = require("ffmpeg-static");
let GeoGridUtils = require("./GeoGridUtils.js");


/**
 * Probes an MP4 file's exact duration in seconds using ffmpeg-static.
 *
 * @param {string} arg0_mp4_path - Input MP4 path.
 * @returns {number} Duration in seconds.
 */
function probeDurationSeconds (arg0_mp4_path) {
  //Convert from parameters
  let mp4_path = arg0_mp4_path;

  //Declare local instance variables
  let res = spawnSync(ffmpeg_path, ["-i", mp4_path], { encoding: "utf8" });
  let m = (res.stderr || "").match(/Duration:\s*(\d+):(\d+):(\d+(?:\.\d+)?)/);

  if (!m) throw new Error(`Could not probe duration for ${mp4_path}`);

  //Return statement
  return Number(m[1])*3600 + Number(m[2])*60 + Number(m[3]);
}


/**
 * Runs an ffmpeg command asynchronously and resolves on exit code 0.
 *
 * @param {Array<string>} arg0_args - ffmpeg CLI arguments.
 * @returns {Promise<void>}
 */
function runFfmpegAsync (arg0_args) {
  //Convert from parameters
  let ff_args = arg0_args;

  //Return statement
  return new Promise((resolve, reject) => {
    let proc = spawn(ffmpeg_path, ff_args, { stdio: ["ignore", "ignore", "pipe"] });
    let err = "";

    proc.stderr.on("data", (arg0_d) => { err += arg0_d.toString(); });
    proc.on("close", (arg0_code) => {
      if (arg0_code === 0) resolve();
      else reject(new Error(`ffmpeg exited ${arg0_code}: ${err.slice(-600)}`));
    });
  });
}


/**
 * Compresses a single MP4 via 2-pass libx264 with automatic retry if size exceeds 19.95 MB.
 *
 * @param {string} arg0_input_mp4 - Source MP4 path.
 * @param {object} arg1_comp_cfg - Compression settings from render_config.json5.
 * @param {string} arg2_pass_dir - Temporary passlog directory.
 * @param {number} arg3_threads - Thread count per ffmpeg instance.
 * @returns {Promise<{input: string, output: string, bytes: number, mb: string, mib: string}>} Result summary.
 */
async function compressSingleVideo (arg0_input_mp4, arg1_comp_cfg, arg2_pass_dir, arg3_threads) {
  //Convert from parameters
  let input_mp4 = arg0_input_mp4;
  let comp_cfg = arg1_comp_cfg;
  let pass_dir = arg2_pass_dir;
  let threads = arg3_threads || 6;

  //Declare local instance variables
  let duration_s = probeDurationSeconds(input_mp4);
  let max_bytes = Math.floor((comp_cfg.max_size_mb || 20.0)*1000000) - 50000;
  let null_sink = (process.platform === "win32") ? "NUL" : "/dev/null";
  let parsed = path.parse(input_mp4);
  let output_mp4 = path.join(parsed.dir, `${parsed.name}${comp_cfg.suffix || "_compressed"}${parsed.ext}`);
  let pass_log = path.join(pass_dir, parsed.name);
  let preset = comp_cfg.preset || "slow";
  let target_mb = comp_cfg.target_size_mb || 19.0;
  let tune = comp_cfg.tune || "animation";

  for (let attempt = 0; attempt < 3; attempt++) {
    let kbps = Math.max(150, Math.floor((target_mb*1000000*8)/(duration_s*1000)));
    let maxrate_k = Math.round(kbps*2.2);
    let bufsize_k = Math.round(kbps*4.4);

    console.log(`[Compress:${parsed.name}] Pass 1/2 (${kbps} kbps, target ${target_mb.toFixed(2)} MB)...`);

    await runFfmpegAsync([
      "-y",
      "-i", input_mp4,
      "-fps_mode", "passthrough",
      "-c:v", "libx264",
      "-pix_fmt", "yuv420p",
      "-preset", preset,
      "-tune", tune,
      "-b:v", `${kbps}k`,
      "-maxrate", `${maxrate_k}k`,
      "-bufsize", `${bufsize_k}k`,
      "-threads", String(threads),
      "-pass", "1",
      "-passlogfile", pass_log,
      "-an",
      "-f", "null",
      null_sink,
    ]);

    console.log(`[Compress:${parsed.name}] Pass 2/2 -> ${output_mp4}`);

    await runFfmpegAsync([
      "-y",
      "-i", input_mp4,
      "-fps_mode", "passthrough",
      "-c:v", "libx264",
      "-pix_fmt", "yuv420p",
      "-preset", preset,
      "-tune", tune,
      "-b:v", `${kbps}k`,
      "-maxrate", `${maxrate_k}k`,
      "-bufsize", `${bufsize_k}k`,
      "-threads", String(threads),
      "-pass", "2",
      "-passlogfile", pass_log,
      "-video_track_timescale", "60000",
      "-movflags", "+faststart",
      output_mp4,
    ]);

    let bytes = fs.statSync(output_mp4).size;
    let mb_dec = bytes/1000000;

    if (bytes <= max_bytes) {
      //Return statement
      return {
        input: input_mp4,
        output: output_mp4,
        bytes: bytes,
        mb: mb_dec.toFixed(2),
        mib: (bytes/(1024*1024)).toFixed(2),
      };
    }

    console.log(`[Compress:${parsed.name}] Output ${mb_dec.toFixed(2)} MB exceeds ${(max_bytes/1000000).toFixed(2)} MB; retrying...`);
    target_mb = target_mb*(18.6/mb_dec);
  }

  throw new Error(`Failed to compress ${input_mp4} below ${max_bytes} bytes after 3 attempts`);
}


/**
 * Compresses all 5 exported Hollinger MP4 videos in parallel.
 */
async function main () {
  //Declare local instance variables
  let rcfg = GeoGridUtils.loadConfig("render_config.json5");
  let anim = rcfg.animation;
  let cat_root = rcfg.categorical_animations;
  let comp_cfg = anim.compression || { target_size_mb: 19.0, max_size_mb: 20.0, suffix: "_compressed", preset: "slow", tune: "animation" };
  let exports_dir = GeoGridUtils.resolvePaths("landarea").exports_dir;
  let pass_dir = fs.mkdtempSync(path.join(os.tmpdir(), "hollinger_2pass_"));
  let target_map = {
    anthromes: path.join(exports_dir, cat_root.anthromes.output_filename),
    holdridge_lifezones: path.join(exports_dir, cat_root.holdridge_lifezones.output_filename),
    koeppen_geiger: path.join(exports_dir, cat_root.koeppen_geiger.output_filename),
    landarea: path.join(exports_dir, anim.output_filename),
    trewartha: path.join(exports_dir, cat_root.trewartha.output_filename),
  };
  let target_arg = process.argv.find((arg0_a) => arg0_a.startsWith("--target="));
  let targets = Object.values(target_map);

  if (target_arg) {
    let key = target_arg.slice(9).trim();
    if (!target_map[key]) throw new Error(`Unknown compression target '${key}'. Valid: ${Object.keys(target_map).join(", ")}`);
    targets = [target_map[key]];
  }

  try {
    let results = await Promise.all(
      targets.map((arg0_p) => compressSingleVideo(arg0_p, comp_cfg, pass_dir, 6))
    );

    for (let i = 0; i < results.length; i++) {
      let r = results[i];

      console.log(`[Compress:OK] ${r.output} | ${r.mb} MB (${r.mib} MiB, ${r.bytes} bytes)`);
    }
  } finally {
    fs.rmSync(pass_dir, { recursive: true, force: true });
  }
}


main().catch((arg0_err) => {
  console.error(arg0_err);
  process.exit(1);
});
