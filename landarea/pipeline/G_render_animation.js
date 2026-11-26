"use strict";

let { spawn } = require("node:child_process");
let fs = require("node:fs");
let path = require("node:path");
let { Worker, isMainThread, parentPort, workerData } = require("node:worker_threads");

let ffmpeg_path = require("ffmpeg-static");
let sharp = require("sharp");

let EqualEarthProjector = require("../../core/EqualEarthProjector");
let GeoGridUtils = require("../../core/GeoGridUtils");

/**
 * Step G: Utilitarian 4K Equal Earth QBAM Animation Renderer (12 026 frames, -10000 to 2025).
 *
 * - Projection: Equal Earth (Savric et al., 2018)
 * - Colourscheme: QBAM-inspired (desaturated lime land, greyish-white glaciers & permanent shelves, baby blue sea)
 * - Label: Minimal top-left Arial 48 pt (internal year 0 displays as "1 AD") + flat bottom-left colour key
 * - Framerate: 60 FPS up to 1 AD (-10000..-1), 30 FPS from 1 AD onwards (0..2025)
 *
 * Fast look-dev iteration:
 *   node scripts/G_render_animation.js --preview-years -10000,-6000,0,1950,2025
 */
class G_RenderAnimation {
  /**
   * Formats an internal integer year (-10000..2025) into the display video label.
   * Internal year 0 corresponds to 1 AD and displays as "1 AD".
   *
   * @param {number} arg0_year
   * @param {Object} arg1_label_cfg
   *
   * @returns {string}
   */
  static formatDisplayYear (arg0_year, arg1_label_cfg) {
    //Convert from parameters
    let year = arg0_year;
    let lcfg = arg1_label_cfg;

    if (year === 0)
      return lcfg.display_year_zero_as ? lcfg.display_year_zero_as : "1 AD";

    if (lcfg.display_mode === "signed")
      return String(year);

    //Return statement
    return (year < 0) ? `${Math.abs(year)} BC` : `${year} AD`;
  }

  /**
   * Builds the static colour key and title RGBA buffer once for fast alpha blitting onto every frame.
   *
   * @param {Object} arg0_rcfg
   * @param {number} arg1_map_top_px - Top edge of the projected map (px).
   *
   * @returns {Promise<{ buf: Buffer, height: number, width: number }>}
   */
  static async buildLegendOverlay (arg0_rcfg, arg1_map_top_px) {
    //Convert from parameters
    let rcfg = arg0_rcfg;
    let map_top_px = arg1_map_top_px;

    //Declare local instance variables
    let l = rcfg.legend;
    let w = rcfg.canvas.width;
    let h = rcfg.canvas.height;
    let font_px = Math.round((l.font_size_pt*rcfg.year_label.dpi)/72.0);
    let cursor_x = l.x_px;
    let elements = [];
    let swatch_y = h - l.margin_bottom_px - l.swatch_height_px;

    if (l.enabled) {
      for (let i = 0; i < l.items.length; i++) {
        let item = l.items[i];
        let hex = rcfg.palette[item.key];
        let safe_label = item.label.replace(/&/g, "&amp;").replace(/</g, "&lt;");

        elements.push(
          `<rect x="${cursor_x}" y="${swatch_y}" width="${l.swatch_width_px}" height="${l.swatch_height_px}" fill="${hex}" stroke="${l.border_color}" stroke-width="2" />`
        );
        let text_x = cursor_x + l.swatch_width_px + l.label_gap_px;
        let text_y = swatch_y + l.swatch_height_px - 3;
        elements.push(
          `<text x="${text_x}" y="${text_y}" font-family="${l.font_family}" font-size="${font_px}" fill="${l.text_color}">${safe_label}</text>`
        );
        cursor_x = text_x + Math.round(item.label.length*font_px*0.47) + l.item_gap_px;
      }
    }

    //Centred top title, vertically centred between the frame top and the top of the map
    if (rcfg.title && rcfg.title.enabled) {
      let tcfg = rcfg.title;
      let title_font_px = Math.round((tcfg.font_size_pt*rcfg.year_label.dpi)/72.0);
      let title_y = Math.round(map_top_px/2);
      let safe_title = tcfg.text.replace(/&/g, "&amp;").replace(/</g, "&lt;");

      elements.push(
        `<text x="${Math.round(w/2)}" y="${title_y}" text-anchor="middle" dominant-baseline="central" font-family="${tcfg.font_family}" font-size="${title_font_px}" fill="${tcfg.color}">${safe_title}</text>`
      );
    }

    let svg = Buffer.from(`<svg width="${w}" height="${h}" xmlns="http://www.w3.org/2000/svg">${elements.join("")}</svg>`);
    let { data } = await sharp(svg).raw().ensureAlpha().toBuffer({ resolveWithObject: true });

    //Return statement
    return { buf: data, height: h, width: w };
  }

  /**
   * Builds a small cropped top-left RGBA overlay for a single year label (Arial 48pt).
   *
   * @param {string} arg0_text
   * @param {Object} arg1_rcfg
   *
   * @returns {Promise<{ buf: Buffer, height: number, width: number }>}
   */
  static async buildYearLabelPatch (arg0_text, arg1_rcfg) {
    //Convert from parameters
    let text = arg0_text;
    let rcfg = arg1_rcfg;

    //Declare local instance variables
    let ycfg = rcfg.year_label;
    let font_px = Math.round((ycfg.font_size_pt*ycfg.dpi)/72.0);
    let patch_w = Math.round(font_px*6.5);
    let patch_h = Math.round(font_px*1.6);
    let safe_text = text.replace(/&/g, "&amp;").replace(/</g, "&lt;");
    let svg = Buffer.from(
      `<svg width="${patch_w}" height="${patch_h}" xmlns="http://www.w3.org/2000/svg">` +
      `<text x="0" y="${Math.round(patch_h/2)}" dominant-baseline="central" font-family="${ycfg.font_family}" font-size="${font_px}" fill="${ycfg.color}">${safe_text}</text>` +
      `</svg>`
    );

    let { data } = await sharp(svg).raw().ensureAlpha().toBuffer({ resolveWithObject: true });

    //Return statement
    return { buf: data, height: patch_h, width: patch_w };
  }

  /**
   * Alpha-blits an RGBA overlay patch directly onto an RGB24 frame buffer in-place.
   *
   * @param {Uint8Array} arg0_frame_rgb
   * @param {number} arg1_frame_w
   * @param {Buffer} arg2_patch_rgba
   * @param {number} arg3_patch_w
   * @param {number} arg4_patch_h
   * @param {number} arg5_dst_x
   * @param {number} arg6_dst_y
   */
  static blitRgbaOntoRgb (arg0_frame_rgb, arg1_frame_w, arg2_patch_rgba, arg3_patch_w, arg4_patch_h, arg5_dst_x, arg6_dst_y) {
    //Convert from parameters
    let frame = arg0_frame_rgb;
    let frame_w = arg1_frame_w;
    let patch = arg2_patch_rgba;
    let pw = arg3_patch_w;
    let ph = arg4_patch_h;
    let dx = arg5_dst_x;
    let dy = arg6_dst_y;

    for (let y = 0; y < ph; y++) {
      let dst_row = (dy + y)*frame_w;
      let src_row = y*pw;
      for (let x = 0; x < pw; x++) {
        let s4 = (src_row + x)*4;
        let a = patch[s4 + 3];
        if (a === 0) continue;
        let d3 = (dst_row + dx + x)*3;
        if (a === 255) {
          frame[d3] = patch[s4];
          frame[d3 + 1] = patch[s4 + 1];
          frame[d3 + 2] = patch[s4 + 2];
        } else {
          let af = a/255.0;
          let inv = 1.0 - af;
          frame[d3] = (frame[d3]*inv + patch[s4]*af + 0.5) | 0;
          frame[d3 + 1] = (frame[d3 + 1]*inv + patch[s4 + 1]*af + 0.5) | 0;
          frame[d3 + 2] = (frame[d3 + 2]*inv + patch[s4 + 2]*af + 0.5) | 0;
        }
      }
    }
  }

  /**
   * Sign-safe interpolation between two signed km2 keyframes (prevents glaciated->land crossing 0 sea).
   *
   * @param {Float32Array} arg0_left
   * @param {Float32Array} arg1_right
   * @param {number} arg2_t - Interpolation weight in [0, 1]
   * @param {Float32Array} arg3_out
   *
   * @returns {Float32Array}
   */
  static interpolateSignedSafe (arg0_left, arg1_right, arg2_t, arg3_out) {
    //Convert from parameters
    let left = arg0_left;
    let right = arg1_right;
    let t = arg2_t;
    let out = arg3_out;

    //Declare local instance variables
    let n = left.length;
    let w0 = 1.0 - t;

    if (t <= 0.0) {
      out.set(left);
      return out;
    }
    if (t >= 1.0) {
      out.set(right);
      return out;
    }

    for (let i = 0; i < n; i++) {
      let a = left[i];
      let b = right[i];

      if (a === b) {
        out[i] = a;
        continue;
      }

      //Separate positive land and negative glaciated magnitudes so transitions obey Land > Glaciated > Sea
      let la = (a > 0.0) ? a : 0.0;
      let lb = (b > 0.0) ? b : 0.0;
      let land_interp = w0*la + t*lb;

      if (land_interp > 0.01) {
        out[i] = land_interp;
      } else {
        let ia = (a < 0.0) ? -a : 0.0;
        let ib = (b < 0.0) ? -b : 0.0;
        let ice_interp = w0*ia + t*ib;
        out[i] = (ice_interp > 0.01) ? -ice_interp : 0.0;
      }
    }

    //Return statement
    return out;
  }

  /**
   * Parses palette RGB triplets from render_config.json5.
   *
   * @param {Object} arg0_rcfg
   *
   * @returns {Object}
   */
  static buildPaletteRgb (arg0_rcfg) {
    //Convert from parameters
    let p = arg0_rcfg.palette;

    //Return statement
    return {
      bg: GeoGridUtils.parseHexRgb(p.background),
      border: p.globe_border ? GeoGridUtils.parseHexRgb(p.globe_border) : null,
      glacier: GeoGridUtils.parseHexRgb(p.glacier_full),
      glacier_blend: p.glacier_coastal_blend,
      land: GeoGridUtils.parseHexRgb(p.land_full),
      land_blend: p.land_coastal_blend,
      sea: GeoGridUtils.parseHexRgb(p.sea),
    };
  }

  /**
   * Renders a contiguous range of annual frames [start_year..end_year] directly into an MP4 segment via ffmpeg.
   * Executed inside worker_threads for parallel 4K encoding.
   *
   * @param {Object} arg0_job
   *
   * @returns {Promise<string>} Output segment path
   */
  static async renderSegmentWorker (arg0_job) {
    //Convert from parameters
    let job = arg0_job;

    //Declare local instance variables
    let cfg = GeoGridUtils.loadConfig("model_config.json5");
    let rcfg = GeoGridUtils.loadConfig("render_config.json5");
    let cell_area = GeoGridUtils.computeCellAreaByRow(cfg.grid);
    let key_years = GeoGridUtils.getTimestepYears(cfg);
    let pal = this.buildPaletteRgb(rcfg);
    let paths = GeoGridUtils.resolvePaths("landarea");
    let prefix = cfg.encoding.geopng_prefix;
    let projector = new EqualEarthProjector(rcfg.canvas, cfg.grid);
    let legend = await this.buildLegendOverlay(rcfg, projector.map_top_px);
    let frame_rgb = new Uint8Array(rcfg.canvas.width*rcfg.canvas.height*3);
    let interp_buf = new Float32Array(cfg.grid.width*cfg.grid.height);

    //Spawn ffmpeg streaming pipe with unified 60000 timescale for lossless VFR concat
    let args = [
      "-y",
      "-f", "rawvideo",
      "-pix_fmt", "rgb24",
      "-s", `${rcfg.canvas.width}x${rcfg.canvas.height}`,
      "-r", String(job.fps),
      "-i", "pipe:0",
      "-c:v", rcfg.animation.video_codec,
      "-preset", rcfg.animation.preset,
      "-crf", String(rcfg.animation.crf),
      "-pix_fmt", rcfg.animation.pixel_format,
      "-threads", "3",
      "-video_track_timescale", "60000",
      "-movflags", "+faststart",
      job.out_path,
    ];

    let ff = spawn(ffmpeg_path, args, { stdio: ["pipe", "ignore", "pipe"] });
    let ff_err = "";
    ff.stderr.on("data", (chunk) => { ff_err += chunk.toString(); });

    let cache = new Map();
    let loadKeyframe = (arg0_y) => {
      if (cache.has(arg0_y)) return cache.get(arg0_y);
      let data = GeoGridUtils.decodeFloat32GeoPNG(path.join(paths.rasters_dir, `${prefix}${arg0_y}.png`)).data;
      if (cache.size >= 4) {
        let oldest = cache.keys().next().value;
        cache.delete(oldest);
      }
      cache.set(arg0_y, data);
      return data;
    };

    let k_ptr = 0;
    let step_years = job.step_years || 1;
    let total_frames = Math.floor((job.end_year - job.start_year)/step_years) + 1;
    let rendered = 0;

    for (let yr = job.start_year; yr <= job.end_year; yr += step_years) {
      for (; k_ptr < key_years.length - 2 && key_years[k_ptr + 1] < yr; )
        k_ptr++;

      let y0 = key_years[k_ptr];
      let y1 = key_years[Math.min(key_years.length - 1, k_ptr + 1)];
      let slice;

      if (yr === y0 || y0 === y1) {
        slice = loadKeyframe(y0);
      } else if (yr === y1) {
        slice = loadKeyframe(y1);
      } else {
        let left_arr = loadKeyframe(y0);
        let right_arr = loadKeyframe(y1);
        let t = (yr - y0)/(y1 - y0);
        slice = this.interpolateSignedSafe(left_arr, right_arr, t, interp_buf);
      }

      projector.projectToRgb(slice, cell_area, pal, frame_rgb);

      if (rcfg.legend.enabled)
        this.blitRgbaOntoRgb(frame_rgb, rcfg.canvas.width, legend.buf, legend.width, legend.height, 0, 0);

      if (rcfg.year_label.enabled) {
        let label = this.formatDisplayYear(yr, rcfg.year_label);
        let patch = await this.buildYearLabelPatch(label, rcfg);
        this.blitRgbaOntoRgb(frame_rgb, rcfg.canvas.width, patch.buf, patch.width, patch.height, rcfg.year_label.x_px, Math.round(projector.map_top_px/2 - patch.height/2));
      }

      let repeats = job.duplicate_frames ? 2 : 1;
      for (let r = 0; r < repeats; r++) {
        let ok = ff.stdin.write(Buffer.from(frame_rgb.buffer, frame_rgb.byteOffset, frame_rgb.byteLength));
        if (!ok)
          await new Promise((resolve) => ff.stdin.once("drain", resolve));
      }

      rendered++;
      if (rendered % 250 === 0 || rendered === total_frames)
        console.log(`[G_render] Segment ${job.id}: ${rendered}/${total_frames} frames (year ${yr}, ${job.fps} fps)`);
    }

    ff.stdin.end();
    let exit_code = await new Promise((resolve) => ff.on("close", resolve));
    if (exit_code !== 0)
      console.error(`[G_render] ffmpeg segment ${job.id} exited with code ${exit_code}:\n${ff_err.slice(-600)}`);

    //Return statement
    return job.out_path;
  }

  /**
   * Renders fast 4K PNG preview stills for specified years so styling in config/render_config.json5 can be iterated quickly.
   *
   * @param {Array<number>} arg0_years
   *
   * @returns {Promise<void>}
   */
  static async renderPreviews (arg0_years) {
    //Convert from parameters
    let preview_years = arg0_years;

    //Declare local instance variables
    let cfg = GeoGridUtils.loadConfig("model_config.json5");
    let rcfg = GeoGridUtils.loadConfig("render_config.json5");
    let cell_area = GeoGridUtils.computeCellAreaByRow(cfg.grid);
    let pal = this.buildPaletteRgb(rcfg);
    let paths = GeoGridUtils.resolvePaths("landarea");
    let preview_dir = path.join(paths.renders_dir, "previews");
    let prefix = cfg.encoding.geopng_prefix;
    let projector = new EqualEarthProjector(rcfg.canvas, cfg.grid);
    let legend = await this.buildLegendOverlay(rcfg, projector.map_top_px);
    let frame_rgb = new Uint8Array(rcfg.canvas.width*rcfg.canvas.height*3);

    if (!fs.existsSync(preview_dir))
      fs.mkdirSync(preview_dir, { recursive: true });

    for (let i = 0; i < preview_years.length; i++) {
      let yr = preview_years[i];
      let slice = GeoGridUtils.decodeFloat32GeoPNG(path.join(paths.rasters_dir, `${prefix}${yr}.png`)).data;
      projector.projectToRgb(slice, cell_area, pal, frame_rgb);

      if (rcfg.legend.enabled)
        this.blitRgbaOntoRgb(frame_rgb, rcfg.canvas.width, legend.buf, legend.width, legend.height, 0, 0);

      if (rcfg.year_label.enabled) {
        let label = this.formatDisplayYear(yr, rcfg.year_label);
        let patch = await this.buildYearLabelPatch(label, rcfg);
        this.blitRgbaOntoRgb(frame_rgb, rcfg.canvas.width, patch.buf, patch.width, patch.height, rcfg.year_label.x_px, Math.round(projector.map_top_px/2 - patch.height/2));
      }

      let out_png = path.join(preview_dir, `preview_${yr}.png`);
      await sharp(Buffer.from(frame_rgb.buffer, frame_rgb.byteOffset, frame_rgb.byteLength), {
        raw: { channels: 3, height: rcfg.canvas.height, width: rcfg.canvas.width },
      })
        .png({ compressionLevel: 6 })
        .toFile(out_png);

      console.log(`[G_render] Saved 4K preview: ${out_png}`);
    }
  }

  /**
   * Orchestrates parallel chunked rendering of all 12 026 annual frames (-10000..2025)
   * at 60 FPS up to 1 AD (-10000..-1) and 30 FPS from 1 AD onwards (0..2025).
   *
   * @returns {Promise<void>}
   */
  static async run () {
    //Check CLI preview mode
    let preview_idx = process.argv.indexOf("--preview-years");
    if (preview_idx !== -1 && process.argv[preview_idx + 1]) {
      let preview_years = process.argv[preview_idx + 1].split(",").map((s) => parseInt(s.trim(), 10));
      await this.renderPreviews(preview_years);
      return;
    }

    //Declare local instance variables
    let rcfg = GeoGridUtils.loadConfig("render_config.json5");
    let anim = rcfg.animation;
    let paths = GeoGridUtils.resolvePaths("landarea");
    let seg_dir = path.join(paths.renders_dir, "segments");
    let switch_yr = anim.fps_switch_year;
    let t0 = Date.now();

    if (!fs.existsSync(seg_dir))
      fs.mkdirSync(seg_dir, { recursive: true });

    //Partition BCE (-10000..-1 at 60 FPS) into 6 chunks and CE (0..2025 at 30 FPS) into 2 chunks (8 parallel workers)
    let jobs = [];
    let schedule = anim.schedule || [];
    let total_planned_frames = 0;

    for (let i = 0; i < schedule.length; i++) {
      let tier = schedule[i];
      let step = tier.step_years || 1;
      let tier_years = [];

      for (let yr = tier.start_year; yr <= tier.end_year; yr += step)
        tier_years.push(yr);

      total_planned_frames += tier_years.length;
      let chunks = tier.chunks || 1;
      let chunk_size = Math.ceil(tier_years.length/chunks);

      for (let c = 0; c < chunks; c++) {
        let slice = tier_years.slice(c*chunk_size, (c + 1)*chunk_size);
        if (slice.length === 0) continue;
        jobs.push({
          duplicate_frames: false,
          end_year: slice[slice.length - 1],
          fps: tier.fps,
          id: `${tier.id}_${c}`,
          out_path: path.join(seg_dir, `seg_${String(jobs.length).padStart(2, "0")}.mp4`),
          start_year: slice[0],
          step_years: step,
        });
      }
    }

    console.log(`[G_render] Launching ${jobs.length} parallel 4K Equal Earth QBAM workers for ${total_planned_frames} frames ...`);

    await Promise.all(jobs.map((job) => new Promise((resolve, reject) => {
      let w = new Worker(__filename, { workerData: job });
      w.on("message", resolve);
      w.on("error", reject);
      w.on("exit", (code) => {
        if (code !== 0) reject(new Error(`Worker ${job.id} exited with code ${code}`));
      });
    })));

    //Concatenate segments losslessly into final MP4
    let concat_list_path = path.join(seg_dir, "concat.txt");
    let concat_lines = jobs.map((j) => `file '${j.out_path.replace(/\\/g, "/")}'`);
    fs.writeFileSync(concat_list_path, concat_lines.join("\n"));

    let final_mp4 = path.join(paths.renders_dir, anim.output_filename);
    console.log(`[G_render] Concatenating ${jobs.length} segments -> ${final_mp4} ...`);

    let concat_proc = spawn(ffmpeg_path, [
      "-y",
      "-f", "concat",
      "-safe", "0",
      "-i", concat_list_path,
      "-c", "copy",
      "-movflags", "+faststart",
      final_mp4,
    ], { stdio: "inherit" });

    let code = await new Promise((resolve) => concat_proc.on("close", resolve));
    if (code !== 0) {
      console.error(`[G_render] Concat failed with exit code ${code}`);
      process.exit(code);
    }

    fs.rmSync(seg_dir, { recursive: true, force: true });

    //Also generate 5 standard preview PNGs for immediate visual inspection
    await this.renderPreviews([-10000, -6000, 0, 1950, 2025]);

    let mp4_mb = (fs.statSync(final_mp4).size/1e6).toFixed(1);
    console.log(`[G_render] Completed ${final_mp4} (${mp4_mb} MB) in ${((Date.now() - t0)/1000).toFixed(1)}s.`);
  }
}

if (!isMainThread) {
  G_RenderAnimation.renderSegmentWorker(workerData)
    .then((res) => parentPort.postMessage(res))
    .catch((err) => {
      console.error("[G_render Worker] Error:", err);
      process.exit(1);
    });
} else if (require.main === module) {
  G_RenderAnimation.run().catch((err) => {
    console.error("[G_render] Fatal error:", err);
    process.exit(1);
  });
}

module.exports = G_RenderAnimation;
