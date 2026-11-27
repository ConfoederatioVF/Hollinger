/**
 * @fileoverview Shared 4K Equal Earth categorical animation renderer for Hollinger indicators
 * (Trewartha, Köppen-Geiger, Holdridge Lifezones, Anthromes). Matches landarea's exact typography
 * (non-bold Arial, 72pt year label, 49.5pt title, #D8DEE6 legend text, sharp #4A5562 swatches),
 * borderless anti-aliased Equal Earth globe, 12 026-frame dual-FPS VFR schedule, and zero-occlusion
 * corner + bottombar legend packing.
 */

"use strict";

let fs = require("node:fs");
let path = require("node:path");
let { spawn } = require("node:child_process");
let { Worker, isMainThread, parentPort, workerData } = require("node:worker_threads");
let sharp = require("sharp");
let ffmpeg_path = require("ffmpeg-static");
let EqualEarthProjector = require("./EqualEarthProjector.js");
let GeoGridUtils = require("./GeoGridUtils.js");
let SpatialTransitionMorph = require("./SpatialTransitionMorph.js");

let text_width_cache = new Map();
let ICE_CLASS_IDS = {
  anthromes: 63,
  holdridge_lifezones: 1,
  koeppen_geiger: 30,
  trewartha: 17,
};


/**
 * Escapes XML special characters for SVG text nodes.
 *
 * @param {string} arg0_text - Raw text string.
 * @returns {string} Escaped XML text.
 */
function escapeXml (arg0_text) {
  //Convert from parameters
  let text = String(arg0_text);

  //Return statement
  return text
    .replace(/&/g, "&amp;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;")
    .replace(/"/g, "&quot;")
    .replace(/'/g, "&apos;");
}


/**
 * Formats integer simulation year for display (0 -> 1 AD, negative -> BC, positive -> AD).
 *
 * @param {number} arg0_year - Simulation year (-10000..2025).
 * @param {object} arg1_label_cfg - Year label configuration.
 * @returns {string} Formatted year label.
 */
function formatDisplayYear (arg0_year, arg1_label_cfg) {
  //Convert from parameters
  let year = arg0_year;
  let lcfg = arg1_label_cfg || {};

  //Return statement
  if (year === 0) return lcfg.display_year_zero_as ? lcfg.display_year_zero_as : "1 AD";
  if (lcfg.display_mode === "signed") return String(year);
  return (year < 0) ? `${Math.abs(year)} BC` : `${year} AD`;
}


/**
 * Measures exact rendered pixel width of a single SVG text string using sharp.
 *
 * @param {string} arg0_text - Text to measure.
 * @param {string} arg1_font_family - Font family name.
 * @param {number} arg2_font_px - Font size in pixels.
 * @returns {Promise<number>} Exact rendered width in pixels.
 */
async function measureTextWidthPx (arg0_text, arg1_font_family, arg2_font_px) {
  //Convert from parameters
  let text = String(arg0_text);
  let font_family = arg1_font_family;
  let font_px = arg2_font_px;

  //Declare local instance variables
  let cache_key = `${font_family}|${font_px}|${text}`;

  if (text_width_cache.has(cache_key)) return text_width_cache.get(cache_key);

  let pad_w = Math.max(128, Math.ceil(text.length*font_px*1.1) + 64);
  let pad_h = Math.ceil(font_px*1.8);
  let svg = Buffer.from(
    `<svg xmlns="http://www.w3.org/2000/svg" width="${pad_w}" height="${pad_h}">` +
    `<text x="16" y="${Math.round(font_px*1.2)}" fill="#FFFFFF" font-family="${escapeXml(font_family)}" font-size="${font_px}">${escapeXml(text)}</text>` +
    `</svg>`
  );
  let { data, info } = await sharp(svg).ensureAlpha().raw().toBuffer({ resolveWithObject: true });
  let min_x = info.width;
  let max_x = -1;

  for (let y = 0; y < info.height; y++) {
    let row_off = y*info.width;

    for (let x = 0; x < info.width; x++) {
      if (data[(row_off + x)*4 + 3] > 8) {
        if (x < min_x) min_x = x;
        if (x > max_x) max_x = x;
      }
    }
  }

  let measured_w = (max_x >= min_x) ? (max_x - min_x + 1) : 0;

  text_width_cache.set(cache_key, measured_w);

  //Return statement
  return measured_w;
}


/**
 * Loads categorical class metadata and builds an RGB lookup table.
 *
 * @param {object} arg0_target_cfg - Target entry from render_config.categorical_animations.
 * @param {object} arg1_cat_root_cfg - Root categorical_animations config.
 * @returns {{geopng_prefix: string, class_map: Map<number, object>, lut_rgb: Uint8Array}} Lookup bundle.
 */
function buildClassLookup (arg0_target_cfg, arg1_cat_root_cfg) {
  //Convert from parameters
  let target_cfg = arg0_target_cfg;
  let cat_root_cfg = arg1_cat_root_cfg;

  //Declare local instance variables
  let class_map = new Map();
  let ocean_rgb = cat_root_cfg.ocean_override_rgb || [16, 42, 67];
  let raw_cfg = GeoGridUtils.loadConfig(target_cfg.config_file);
  let sub_cfg = target_cfg.config_key ? raw_cfg[target_cfg.config_key] : raw_cfg;
  let classes = sub_cfg.classes || [];
  let geopng_prefix = (sub_cfg.encoding && sub_cfg.encoding.geopng_prefix)
    ? sub_cfg.encoding.geopng_prefix
    : (sub_cfg.geopng_prefix || `${target_cfg.preview_prefix}_`);
  let max_id = 0;

  let overrides = target_cfg.rgb_overrides || {};

  for (let i = 0; i < classes.length; i++) {
    let cls = classes[i];
    let rgb = (cls.id === 0)
      ? ocean_rgb
      : (overrides[String(cls.id)] ? overrides[String(cls.id)] : cls.rgb);

    if (cls.id > max_id) max_id = cls.id;
    class_map.set(cls.id, {
      id: cls.id,
      code: cls.code || String(cls.id),
      name: cls.name || String(cls.id),
      rgb: rgb,
    });
  }

  let lut_rgb = new Uint8Array((max_id + 1)*3);

  for (let i = 0; i <= max_id; i++) {
    let entry = class_map.get(i);

    if (entry) {
      lut_rgb[i*3] = entry.rgb[0];
      lut_rgb[i*3 + 1] = entry.rgb[1];
      lut_rgb[i*3 + 2] = entry.rgb[2];
    }
  }

  //Return statement
  return {
    geopng_prefix: geopng_prefix,
    class_map: class_map,
    lut_rgb: lut_rgb,
  };
}


/**
 * Renders the static title + corner-packed + bottom-bar categorical legend overlay as an RGBA buffer.
 *
 * @param {object} arg0_target_cfg - Target animation config.
 * @param {object} arg1_render_cfg - Global render_config.
 * @param {Map<number, object>} arg2_class_map - Class ID metadata map.
 * @param {number} arg3_map_top_px - Projected Equal Earth map top pixel Y.
 * @returns {Promise<{data: Buffer, width: number, height: number, svg_str: string, bboxes: Array<object>}>} RGBA overlay bundle.
 */
async function buildStaticOverlay (arg0_target_cfg, arg1_render_cfg, arg2_class_map, arg3_map_top_px) {
  //Convert from parameters
  let target_cfg = arg0_target_cfg;
  let render_cfg = arg1_render_cfg;
  let class_map = arg2_class_map;
  let map_top_px = arg3_map_top_px;

  //Declare local instance variables
  let bboxes = [];
  let c = render_cfg.canvas;
  let cat_root = render_cfg.categorical_animations;
  let corner_layout = target_cfg.corner_layout;
  let corner_names = ["top_left", "top_right", "bottom_left", "bottom_right"];
  let dpi = cat_root.dpi || 96;
  let font_family = cat_root.font_family || "Arial";
  let font_px = Math.round(((target_cfg.font_size_pt || 24)*dpi)/72.0);
  let label_gap = target_cfg.label_gap_px || 12;
  let svg_parts = [];
  let swatch_border = cat_root.swatch_border_color || "#4A5562";
  let swatch_border_w = cat_root.swatch_border_width_px || 2;
  let swatch_h = target_cfg.swatch_height_px || 26;
  let swatch_w = target_cfg.swatch_width_px || 44;
  let text_color = cat_root.text_color || "#D8DEE6";
  let title_cfg = render_cfg.title;
  let title_color = title_cfg.color || "#F2F4F7";
  let title_font = title_cfg.font_family || "Arial";
  let title_px = Math.round(((title_cfg.font_size_pt || 49.5)*dpi)/72.0);
  let title_y = Math.round(map_top_px/2);

  svg_parts.push(
    `<svg xmlns="http://www.w3.org/2000/svg" width="${c.width}" height="${c.height}">`,
    `<text x="${Math.round(c.width/2)}" y="${title_y}" fill="${title_color}" font-family="${escapeXml(title_font)}" font-size="${title_px}" text-anchor="middle" dominant-baseline="central">${escapeXml(target_cfg.title_text)}</text>`
  );

  //Render 4 Equal Earth corners
  for (let i = 0; i < corner_names.length; i++) {
    let corner_key = corner_names[i];
    let block = corner_layout[corner_key];

    if (!block) continue;

    let align_right = Boolean(block.align_right);
    let col_count = block.columns || 1;
    let col_gap = block.col_gap_px || 28;
    let cur_y = block.y;
    let line_h = block.line_height_px || target_cfg.line_height_px || 42;

    for (let x = 0; x < block.groups.length; x++) {
      let grp = block.groups[x];

      if (grp.header) {
        let hdr_w = await measureTextWidthPx(grp.header, font_family, font_px);

        if (align_right) {
          svg_parts.push(
            `<text x="${block.x}" y="${cur_y}" fill="${text_color}" font-family="${escapeXml(font_family)}" font-size="${font_px}" text-anchor="end">${escapeXml(grp.header)}</text>`
          );
          bboxes.push({ label: `${corner_key}:${grp.header}`, x0: block.x - hdr_w, y0: cur_y - font_px, x1: block.x, y1: cur_y + 4 });
        } else {
          svg_parts.push(
            `<text x="${block.x}" y="${cur_y}" fill="${text_color}" font-family="${escapeXml(font_family)}" font-size="${font_px}" text-anchor="start">${escapeXml(grp.header)}</text>`
          );
          bboxes.push({ label: `${corner_key}:${grp.header}`, x0: block.x, y0: cur_y - font_px, x1: block.x + hdr_w, y1: cur_y + 4 });
        }

        cur_y += line_h;
      }

      let items = grp.items || [];
      let rows_per_col = Math.ceil(items.length/col_count);
      let col_widths = [];

      for (let col_i = 0; col_i < col_count; col_i++) {
        let max_item_w = 0;

        for (let r_i = 0; r_i < rows_per_col; r_i++) {
          let item_idx = col_i*rows_per_col + r_i;

          if (item_idx >= items.length) break;

          let it = items[item_idx];
          let tw = await measureTextWidthPx(it.label, font_family, font_px);
          let total_w = swatch_w + label_gap + tw;

          if (total_w > max_item_w) max_item_w = total_w;
        }
        col_widths.push(max_item_w);
      }

      let col_offsets = [0];

      for (let col_i = 1; col_i < col_count; col_i++) {
        col_offsets.push(col_offsets[col_i - 1] + col_widths[col_i - 1] + col_gap);
      }

      for (let y = 0; y < items.length; y++) {
        let it = items[y];
        let cls = class_map.get(it.id);

        if (!cls) continue;

        let col_idx = Math.floor(y/rows_per_col);
        let row_idx = y%rows_per_col;
        let item_y = cur_y + row_idx*line_h;
        let swatch_y = item_y - Math.round(swatch_h*0.84);
        let fill_css = `rgb(${cls.rgb[0]},${cls.rgb[1]},${cls.rgb[2]})`;
        let tw = await measureTextWidthPx(it.label, font_family, font_px);

        if (align_right) {
          let item_x = block.x - col_offsets[col_idx];
          let sw_x = item_x - swatch_w;
          let txt_x = sw_x - label_gap;

          svg_parts.push(
            `<rect x="${sw_x}" y="${swatch_y}" width="${swatch_w}" height="${swatch_h}" fill="${fill_css}" stroke="${swatch_border}" stroke-width="${swatch_border_w}"/>`,
            `<text x="${txt_x}" y="${item_y}" fill="${text_color}" font-family="${escapeXml(font_family)}" font-size="${font_px}" text-anchor="end">${escapeXml(it.label)}</text>`
          );
          bboxes.push({ label: `${corner_key}:${it.label}`, x0: txt_x - tw, y0: swatch_y, x1: item_x, y1: Math.max(item_y + 4, swatch_y + swatch_h) });
        } else {
          let item_x = block.x + col_offsets[col_idx];
          let txt_x = item_x + swatch_w + label_gap;

          svg_parts.push(
            `<rect x="${item_x}" y="${swatch_y}" width="${swatch_w}" height="${swatch_h}" fill="${fill_css}" stroke="${swatch_border}" stroke-width="${swatch_border_w}"/>`,
            `<text x="${txt_x}" y="${item_y}" fill="${text_color}" font-family="${escapeXml(font_family)}" font-size="${font_px}" text-anchor="start">${escapeXml(it.label)}</text>`
          );
          bboxes.push({ label: `${corner_key}:${it.label}`, x0: item_x, y0: swatch_y, x1: txt_x + tw, y1: Math.max(item_y + 4, swatch_y + swatch_h) });
        }
      }

      cur_y += rows_per_col*line_h + Math.round(line_h*0.35);
    }
  }

  //Render horizontal bottom_bar legend rows
  if (corner_layout.bottom_bar && corner_layout.bottom_bar.rows) {
    let bb = corner_layout.bottom_bar;
    let group_gap = target_cfg.group_gap_px || 48;
    let item_gap = target_cfg.item_gap_px || 32;

    for (let r_i = 0; r_i < bb.rows.length; r_i++) {
      let row = bb.rows[r_i];
      let cursor_x = bb.x;
      let row_y = row.y;
      let swatch_y = row_y - Math.round(swatch_h*0.84);

      for (let g_i = 0; g_i < row.groups.length; g_i++) {
        let grp = row.groups[g_i];

        if (grp.header) {
          let hdr_w = await measureTextWidthPx(grp.header, font_family, font_px);

          svg_parts.push(
            `<text x="${cursor_x}" y="${row_y}" fill="${text_color}" font-family="${escapeXml(font_family)}" font-size="${font_px}" text-anchor="start">${escapeXml(grp.header)}</text>`
          );
          bboxes.push({ label: `bottom_bar:${grp.header}`, x0: cursor_x, y0: row_y - font_px, x1: cursor_x + hdr_w, y1: row_y + 4 });
          cursor_x += hdr_w + Math.round(item_gap*0.75);
        }

        for (let i_i = 0; i_i < grp.items.length; i_i++) {
          let it = grp.items[i_i];
          let cls = class_map.get(it.id);

          if (!cls) continue;

          let fill_css = `rgb(${cls.rgb[0]},${cls.rgb[1]},${cls.rgb[2]})`;
          let tw = await measureTextWidthPx(it.label, font_family, font_px);
          let txt_x = cursor_x + swatch_w + label_gap;

          svg_parts.push(
            `<rect x="${cursor_x}" y="${swatch_y}" width="${swatch_w}" height="${swatch_h}" fill="${fill_css}" stroke="${swatch_border}" stroke-width="${swatch_border_w}"/>`,
            `<text x="${txt_x}" y="${row_y}" fill="${text_color}" font-family="${escapeXml(font_family)}" font-size="${font_px}" text-anchor="start">${escapeXml(it.label)}</text>`
          );
          bboxes.push({ label: `bottom_bar:${it.label}`, x0: cursor_x, y0: swatch_y, x1: txt_x + tw, y1: Math.max(row_y + 4, swatch_y + swatch_h) });
          cursor_x = txt_x + tw + item_gap;
        }

        cursor_x += (group_gap - item_gap);
      }
    }
  }

  svg_parts.push("</svg>");

  let svg_str = svg_parts.join("\n");
  let { data, info } = await sharp(Buffer.from(svg_str))
    .ensureAlpha()
    .raw()
    .toBuffer({ resolveWithObject: true });

  //Return statement
  return {
    data: data,
    width: info.width,
    height: info.height,
    svg_str: svg_str,
    bboxes: bboxes,
  };
}


/**
 * Builds an RGBA patch buffer for the top-left year label matching landarea's exact styling.
 *
 * @param {number} arg0_year - Simulation year.
 * @param {object} arg1_render_cfg - Global render_config.
 * @param {number} arg2_map_top_px - Projected map top pixel Y.
 * @returns {Promise<{data: Buffer, width: number, height: number, x: number, y: number}>} Patch RGBA.
 */
async function buildYearLabelPatch (arg0_year, arg1_render_cfg, arg2_map_top_px) {
  //Convert from parameters
  let year = arg0_year;
  let render_cfg = arg1_render_cfg;
  let map_top_px = arg2_map_top_px;

  //Declare local instance variables
  let ycfg = render_cfg.year_label;
  let colour = ycfg.color || "#F2F4F7";
  let dpi = ycfg.dpi || 96;
  let font = ycfg.font_family || "Arial";
  let font_px = Math.round(((ycfg.font_size_pt || 72)*dpi)/72.0);
  let patch_h = Math.round(font_px*1.6);
  let patch_w = Math.round(font_px*6.5);
  let text = formatDisplayYear(year, ycfg);
  let svg = `<svg width="${patch_w}" height="${patch_h}" xmlns="http://www.w3.org/2000/svg">` +
    `<text x="0" y="${Math.round(patch_h/2)}" dominant-baseline="central" font-family="${escapeXml(font)}" font-size="${font_px}" fill="${colour}">${escapeXml(text)}</text>` +
    `</svg>`;
  let { data, info } = await sharp(Buffer.from(svg))
    .ensureAlpha()
    .raw()
    .toBuffer({ resolveWithObject: true });

  //Return statement
  return {
    data: data,
    width: info.width,
    height: info.height,
    x: ycfg.x_px || 56,
    y: Math.round(map_top_px/2 - info.height/2),
  };
}


/**
 * Alpha-blits an RGBA patch onto a full-frame RGB buffer in place.
 *
 * @param {Uint8Array} arg0_rgb - Destination RGB buffer.
 * @param {number} arg1_frame_w - Destination frame width.
 * @param {number} arg2_frame_h - Destination frame height.
 * @param {Buffer|Uint8Array} arg3_rgba - Source RGBA patch buffer.
 * @param {number} arg4_patch_w - Patch width.
 * @param {number} arg5_patch_h - Patch height.
 * @param {number} arg6_dst_x - Destination X offset.
 * @param {number} arg7_dst_y - Destination Y offset.
 */
function blitRgbaOntoRgb (arg0_rgb, arg1_frame_w, arg2_frame_h, arg3_rgba, arg4_patch_w, arg5_patch_h, arg6_dst_x, arg7_dst_y) {
  //Convert from parameters
  let rgb = arg0_rgb;
  let frame_w = arg1_frame_w;
  let frame_h = arg2_frame_h;
  let rgba = arg3_rgba;
  let patch_w = arg4_patch_w;
  let patch_h = arg5_patch_h;
  let dst_x = arg6_dst_x;
  let dst_y = arg7_dst_y;

  //Declare local instance variables
  let max_x = Math.min(patch_w, frame_w - dst_x);
  let max_y = Math.min(patch_h, frame_h - dst_y);

  for (let y = 0; y < max_y; y++) {
    let fy = dst_y + y;

    if (fy < 0) continue;

    for (let x = 0; x < max_x; x++) {
      let fx = dst_x + x;

      if (fx < 0) continue;

      let p4 = (y*patch_w + x)*4;
      let a = rgba[p4 + 3];

      if (a === 0) continue;

      let d3 = (fy*frame_w + fx)*3;

      if (a === 255) {
        rgb[d3] = rgba[p4];
        rgb[d3 + 1] = rgba[p4 + 1];
        rgb[d3 + 2] = rgba[p4 + 2];
      } else {
        let af = a/255.0;
        let inv = 1.0 - af;

        rgb[d3] = (rgb[d3]*inv + rgba[p4]*af + 0.5) | 0;
        rgb[d3 + 1] = (rgb[d3 + 1]*inv + rgba[p4 + 1]*af + 0.5) | 0;
        rgb[d3 + 2] = (rgb[d3 + 2]*inv + rgba[p4 + 2]*af + 0.5) | 0;
      }
    }
  }
}


/**
 * Builds a dilated globe occupancy mask (edge_alpha > 0) for clearance checking.
 *
 * @param {EqualEarthProjector} arg0_projector - Projector instance.
 * @param {number} arg1_clearance_px - Clearance radius in pixels.
 * @returns {Uint8Array} Binary halo mask (1 = inside globe or clearance halo).
 */
function buildDilatedGlobeMask (arg0_projector, arg1_clearance_px) {
  //Convert from parameters
  let projector = arg0_projector;
  let clearance_px = arg1_clearance_px;

  //Declare local instance variables
  let w = projector.width;
  let h = projector.height;
  let halo = new Uint8Array(w*h);
  let r_sq = clearance_px*clearance_px;

  for (let y = 0; y < h; y++) {
    for (let x = 0; x < w; x++) {
      let idx = y*w + x;

      if (projector.edge_alpha[idx] > 0) {
        let y0 = Math.max(0, y - clearance_px);
        let y1 = Math.min(h - 1, y + clearance_px);
        let x0 = Math.max(0, x - clearance_px);
        let x1 = Math.min(w - 1, x + clearance_px);

        for (let yy = y0; yy <= y1; yy++) {
          let dy = yy - y;

          for (let xx = x0; xx <= x1; xx++) {
            let dx = xx - x;

            if (dx*dx + dy*dy <= r_sq) halo[yy*w + xx] = 1;
          }
        }
      }
    }
  }

  //Return statement
  return halo;
}


/**
 * Validates typography hygiene, minimum font size, pairwise bbox collisions, and pixel-level occlusion.
 *
 * @param {object} arg0_target_cfg - Target config.
 * @param {object} arg1_cat_root - Root categorical_animations config.
 * @param {EqualEarthProjector} arg2_projector - Projector instance.
 * @param {Uint8Array} arg3_halo_mask - Dilated globe halo mask.
 * @param {{data: Buffer, width: number, height: number, svg_str: string, bboxes: Array<object>}} arg4_overlay - Static overlay.
 * @param {{data: Buffer, width: number, height: number, x: number, y: number}} arg5_year_patch - Worst-case year label patch.
 * @returns {object} Comprehensive audit report.
 */
function verifyLayoutAndOcclusion (arg0_target_cfg, arg1_cat_root, arg2_projector, arg3_halo_mask, arg4_overlay, arg5_year_patch) {
  //Convert from parameters
  let target_cfg = arg0_target_cfg;
  let cat_root = arg1_cat_root;
  let projector = arg2_projector;
  let halo_mask = arg3_halo_mask;
  let overlay = arg4_overlay;
  let year_patch = arg5_year_patch;

  //Declare local instance variables
  let bboxes = overlay.bboxes || [];
  let clearance_violations = 0;
  let edge_clip = 0;
  let first_violation = null;
  let globe_overlap = 0;
  let h = overlay.height;
  let item_overlap = 0;
  let label_overlap = 0;
  let margin_x = cat_root.edge_margin_x_px || 40;
  let margin_y = cat_root.edge_margin_y_px || 24;
  let min_pt = cat_root.min_font_size_pt || 24;
  let svg_issues = [];
  let w = overlay.width;

  if ((target_cfg.font_size_pt || 0) < min_pt) {
    svg_issues.push(`font_size_pt=${target_cfg.font_size_pt} < ${min_pt}`);
  }
  if (overlay.svg_str.includes("font-weight")) svg_issues.push("contains font-weight");
  if (overlay.svg_str.includes("\u2014")) svg_issues.push("contains em-dash");
  if (/\br[xy]\s*=/.test(overlay.svg_str)) svg_issues.push("contains rounded rect rx/ry");
  for (let i = 0; i < bboxes.length; i++) {
    if (/\s\/|\/\s/.test(bboxes[i].label)) {
      svg_issues.push(`space around slash in "${bboxes[i].label}"`);
      break;
    }
  }

  for (let i = 0; i < bboxes.length; i++) {
    let a = bboxes[i];

    for (let j = i + 1; j < bboxes.length; j++) {
      let b = bboxes[j];
      let pad = 4;

      if (
        a.x0 - pad < b.x1 &&
        a.x1 + pad > b.x0 &&
        a.y0 - pad < b.y1 &&
        a.y1 + pad > b.y0
      ) {
        item_overlap++;
        if (!first_violation) {
          first_violation = { type: "item_overlap", a: a.label, b: b.label, ax1: a.x1, bx0: b.x0 };
        }
      }
    }
  }

  for (let y = 0; y < h; y++) {
    for (let x = 0; x < w; x++) {
      let idx = y*w + x;
      let a = overlay.data[idx*4 + 3];

      if (a <= 12) continue;
      if (x < margin_x || x >= w - margin_x || y < margin_y || y >= h - margin_y) {
        edge_clip++;
        if (!first_violation) first_violation = { type: "edge_clip", x: x, y: y };
      }
      if (projector.edge_alpha[idx] > 0) {
        globe_overlap++;
        if (!first_violation) first_violation = { type: "globe_overlap", x: x, y: y };
      } else if (halo_mask[idx] === 1) {
        clearance_violations++;
        if (!first_violation) first_violation = { type: "clearance_violation", x: x, y: y };
      }
    }
  }

  for (let y = 0; y < year_patch.height; y++) {
    let fy = year_patch.y + y;

    if (fy < 0 || fy >= h) continue;

    for (let x = 0; x < year_patch.width; x++) {
      let fx = year_patch.x + x;

      if (fx < 0 || fx >= w) continue;

      let pa = year_patch.data[(y*year_patch.width + x)*4 + 3];

      if (pa <= 12) continue;

      let idx = fy*w + fx;

      if (overlay.data[idx*4 + 3] > 12) {
        label_overlap++;
        if (!first_violation) first_violation = { type: "label_overlap", x: fx, y: fy };
      }
      if (projector.edge_alpha[idx] > 0) {
        globe_overlap++;
        if (!first_violation) first_violation = { type: "year_globe_overlap", x: fx, y: fy };
      } else if (halo_mask[idx] === 1) {
        clearance_violations++;
        if (!first_violation) first_violation = { type: "year_clearance_violation", x: fx, y: fy };
      }
    }
  }

  let passed = (
    svg_issues.length === 0 &&
    item_overlap === 0 &&
    globe_overlap === 0 &&
    clearance_violations === 0 &&
    edge_clip === 0 &&
    label_overlap === 0
  );

  //Return statement
  return {
    passed: passed,
    svg_hygiene: svg_issues.length === 0 ? "PASS" : svg_issues.join("; "),
    item_overlap: item_overlap,
    globe_overlap: globe_overlap,
    clearance_violations: clearance_violations,
    edge_clip: edge_clip,
    label_overlap: label_overlap,
    first_violation: first_violation,
  };
}


/**
 * Renders preview PNG frames for a single target and verifies zero occlusion + typography hygiene.
 *
 * @param {string} arg0_target_key - Target key (e.g. "trewartha").
 * @param {Array<number>} arg1_years - Preview years to render.
 * @param {object} arg2_render_cfg - Global render_config.
 * @param {EqualEarthProjector} arg3_projector - Shared EqualEarthProjector.
 * @param {Uint8Array} arg4_halo_mask - Dilated globe halo mask.
 * @param {boolean} arg5_allow_violations - Whether to continue saving PNGs on audit failure.
 */
async function renderTargetPreviews (arg0_target_key, arg1_years, arg2_render_cfg, arg3_projector, arg4_halo_mask, arg5_allow_violations) {
  //Convert from parameters
  let target_key = arg0_target_key;
  let preview_years = arg1_years;
  let render_cfg = arg2_render_cfg;
  let projector = arg3_projector;
  let halo_mask = arg4_halo_mask;
  let allow_violations = Boolean(arg5_allow_violations);

  //Declare local instance variables
  let bg_rgb = GeoGridUtils.parseHexRgb(render_cfg.palette.background);
  let cat_root = render_cfg.categorical_animations;
  let target_cfg = cat_root[target_key];
  let paths = GeoGridUtils.resolvePaths(target_cfg.indicator);
  let lookup = buildClassLookup(target_cfg, cat_root);
  let overlay = await buildStaticOverlay(target_cfg, render_cfg, lookup.class_map, projector.map_top_px);
  let worst_year_patch = await buildYearLabelPatch(-10000, render_cfg, projector.map_top_px);
  let audit = verifyLayoutAndOcclusion(
    target_cfg,
    cat_root,
    projector,
    halo_mask,
    overlay,
    worst_year_patch
  );

  console.log(
    `[Preview Audit:${target_key}] globe_overlap=${audit.globe_overlap}, ` +
    `clearance_violations=${audit.clearance_violations}, edge_clip=${audit.edge_clip}, ` +
    `label_overlap=${audit.label_overlap}, item_overlap=${audit.item_overlap}, svg_hygiene=${audit.svg_hygiene}` +
    (audit.first_violation ? ` first=${JSON.stringify(audit.first_violation)}` : "")
  );

  let previews_dir = path.join(paths.renders_dir, "previews");

  fs.mkdirSync(previews_dir, { recursive: true });

  let model_cfg = GeoGridUtils.loadConfig("model_config.json5");
  let cell_area_by_row = GeoGridUtils.computeCellAreaByRow(model_cfg.grid);
  let land_paths = GeoGridUtils.resolvePaths("landarea");
  let land_prefix = model_cfg.encoding.geopng_prefix;
  let ice_class_id = (target_cfg.ice_class_id !== undefined) ? target_cfg.ice_class_id : (ICE_CLASS_IDS[target_key] || 0);
  let src_w = model_cfg.grid.width;
  let src_h = model_cfg.grid.height;
  let n_src = src_w*src_h;
  let morph_scratch = {
    grid_to_local: new Int32Array(n_src).fill(-1),
    grp_map: new Int32Array(n_src),
    hw_map: new Float32Array(n_src),
    tau_map: new Float32Array(n_src).fill(-1.0),
  };
  let timesteps = GeoGridUtils.getTimestepYears(model_cfg);
  let loadPreviewKeyframe = async (arg0_ts_year) => {
    let ts_year = arg0_ts_year;
    let png_path = path.join(paths.rasters_dir, `${lookup.geopng_prefix}${ts_year}.png`);
    let land_png_path = path.join(land_paths.rasters_dir, `${land_prefix}${ts_year}.png`);
    let decoded = await GeoGridUtils.decodeUint32GeoPNG(png_path);
    let land_decoded = GeoGridUtils.decodeFloat32GeoPNG(land_png_path);
    let rgb = projector.projectCategoricalToRgb(decoded.data, lookup.lut_rgb, bg_rgb, {
      cell_area_by_row: cell_area_by_row,
      ice_class_id: ice_class_id,
      land_km2: land_decoded.data,
    });

    return {
      classes: new Uint8Array(decoded.data),
      land: land_decoded.data,
      rgb: rgb,
    };
  };

  for (let i = 0; i < preview_years.length; i++) {
    let yr = preview_years[i];
    let ts_idx = 0;

    for (; ts_idx + 1 < timesteps.length && timesteps[ts_idx + 1] <= yr; ) {
      ts_idx++;
    }

    let y0 = timesteps[ts_idx];
    let y1 = (ts_idx + 1 < timesteps.length) ? timesteps[ts_idx + 1] : y0;
    let kf0 = await loadPreviewKeyframe(y0);
    let frame_rgb;

    if (y0 === y1 || yr <= y0) {
      frame_rgb = kf0.rgb;
    } else if (yr >= y1) {
      let kf1 = await loadPreviewKeyframe(y1);

      frame_rgb = kf1.rgb;
    } else {
      let kf1 = await loadPreviewKeyframe(y1);
      let src_plan = SpatialTransitionMorph.buildCategoricalPlan(
        kf0.classes,
        kf1.classes,
        kf0.land,
        kf1.land,
        src_w,
        src_h,
        morph_scratch
      );
      let canvas_morph = SpatialTransitionMorph.projectCategoricalPlanToCanvas(
        src_plan,
        projector,
        kf0.rgb,
        kf1.rgb,
        morph_scratch
      );
      let t = (yr - y0)/(y1 - y0);

      frame_rgb = new Uint8Array(kf0.rgb);
      SpatialTransitionMorph.applyCanvasMorph(kf0.rgb, kf1.rgb, canvas_morph, t, frame_rgb);
    }

    blitRgbaOntoRgb(
      frame_rgb,
      projector.width,
      projector.height,
      overlay.data,
      overlay.width,
      overlay.height,
      0,
      0
    );

    let yr_patch = await buildYearLabelPatch(yr, render_cfg, projector.map_top_px);

    blitRgbaOntoRgb(
      frame_rgb,
      projector.width,
      projector.height,
      yr_patch.data,
      yr_patch.width,
      yr_patch.height,
      yr_patch.x,
      yr_patch.y
    );

    let out_path = path.join(previews_dir, `${target_cfg.preview_prefix}_preview_${yr}.png`);

    await sharp(Buffer.from(frame_rgb), {
      raw: { width: projector.width, height: projector.height, channels: 3 },
    })
      .png({ compressionLevel: 6 })
      .toFile(out_path);

    console.log(`[Preview:${target_key}] Saved ${out_path}`);
  }

  if (!audit.passed && !allow_violations) {
    throw new Error(`Layout/occlusion audit failed for ${target_key}: ${JSON.stringify(audit)}`);
  }
}


/**
 * Worker thread entry: renders a contiguous segment of years to an H.264 MP4 chunk
 * matching landarea's exact 60000 timescale dual-FPS VFR stream.
 */
async function runSegmentWorker () {
  //Declare local instance variables
  let {
    target_key,
    segment_id,
    start_year,
    end_year,
    step_years,
    fps,
    duplicate_frames,
    output_mp4,
    render_cfg,
    timesteps,
  } = workerData;
  let step_y = step_years || 1;
  let anim = render_cfg.animation;
  let bg_rgb = GeoGridUtils.parseHexRgb(render_cfg.palette.background);
  let cat_root = render_cfg.categorical_animations;
  let target_cfg = cat_root[target_key];
  let paths = GeoGridUtils.resolvePaths(target_cfg.indicator);
  let lookup = buildClassLookup(target_cfg, cat_root);
  let canvas_cfg = Object.assign({}, render_cfg.canvas, { bilinear_sampling: false });
  let projector = new EqualEarthProjector(canvas_cfg);
  let overlay = await buildStaticOverlay(target_cfg, render_cfg, lookup.class_map, projector.map_top_px);
  let frame_w = projector.width;
  let frame_h = projector.height;
  let n_bytes = frame_w*frame_h*3;
  let blend_buf = new Uint8Array(n_bytes);
  let keyframe_cache = new Map();

  let model_cfg = GeoGridUtils.loadConfig("model_config.json5");
  let cell_area_by_row = GeoGridUtils.computeCellAreaByRow(model_cfg.grid);
  let land_paths = GeoGridUtils.resolvePaths("landarea");
  let land_prefix = model_cfg.encoding.geopng_prefix;
  let ice_class_id = (target_cfg.ice_class_id !== undefined) ? target_cfg.ice_class_id : (ICE_CLASS_IDS[target_key] || 0);
  let src_w = model_cfg.grid.width;
  let src_h = model_cfg.grid.height;
  let n_src = src_w*src_h;
  let morph_scratch = {
    grid_to_local: new Int32Array(n_src).fill(-1),
    grp_map: new Int32Array(n_src),
    hw_map: new Float32Array(n_src),
    tau_map: new Float32Array(n_src).fill(-1.0),
  };

  let getProjectedKeyframe = async (arg0_ts_year) => {
    let ts_year = arg0_ts_year;

    if (keyframe_cache.has(ts_year)) return keyframe_cache.get(ts_year);

    let png_path = path.join(paths.rasters_dir, `${lookup.geopng_prefix}${ts_year}.png`);
    let land_png_path = path.join(land_paths.rasters_dir, `${land_prefix}${ts_year}.png`);
    let decoded = await GeoGridUtils.decodeUint32GeoPNG(png_path);
    let land_decoded = GeoGridUtils.decodeFloat32GeoPNG(land_png_path);
    let rgb = projector.projectCategoricalToRgb(decoded.data, lookup.lut_rgb, bg_rgb, {
      cell_area_by_row: cell_area_by_row,
      ice_class_id: ice_class_id,
      land_km2: land_decoded.data,
    });

    blitRgbaOntoRgb(rgb, frame_w, frame_h, overlay.data, overlay.width, overlay.height, 0, 0);

    if (keyframe_cache.size >= 4) {
      let oldest_key = keyframe_cache.keys().next().value;

      keyframe_cache.delete(oldest_key);
    }
    let entry = { classes: new Uint8Array(decoded.data), land: land_decoded.data, rgb: rgb };
    keyframe_cache.set(ts_year, entry);
    return entry;
  };

  let ff_args = [
    "-y",
    "-f", "rawvideo",
    "-pix_fmt", "rgb24",
    "-s", `${frame_w}x${frame_h}`,
    "-r", String(fps),
    "-i", "pipe:0",
    "-c:v", anim.video_codec,
    "-preset", anim.preset,
    "-crf", String(anim.crf),
    "-pix_fmt", anim.pixel_format,
    "-threads", "3",
    "-video_track_timescale", "60000",
    "-movflags", "+faststart",
    output_mp4,
  ];
  let ff = spawn(ffmpeg_path, ff_args, { stdio: ["pipe", "ignore", "pipe"] });
  let ff_err = "";

  ff.stderr.on("data", (arg0_chunk) => {
    ff_err += arg0_chunk.toString();
  });

  let canvas_morph = null;
  let frame_count = 0;
  let morph_y0 = null;
  let morph_y1 = null;
  let ts_idx = 0;

  for (let yr = start_year; yr <= end_year; yr += step_y) {
    frame_count++;
    while (ts_idx + 1 < timesteps.length && timesteps[ts_idx + 1] <= yr) {
      ts_idx++;
    }

    let y0 = timesteps[ts_idx];
    let y1 = (ts_idx + 1 < timesteps.length) ? timesteps[ts_idx + 1] : y0;
    let kf0 = await getProjectedKeyframe(y0);
    let yr_patch = await buildYearLabelPatch(yr, render_cfg, projector.map_top_px);

    if (y0 === y1 || yr === y0) {
      blend_buf.set(kf0.rgb);
    } else {
      let t = (yr - y0)/(y1 - y0);
      let kf1 = await getProjectedKeyframe(y1);

      if (morph_y0 !== y0 || morph_y1 !== y1) {
        let src_plan = SpatialTransitionMorph.buildCategoricalPlan(
          kf0.classes,
          kf1.classes,
          kf0.land,
          kf1.land,
          src_w,
          src_h,
          morph_scratch
        );
        canvas_morph = SpatialTransitionMorph.projectCategoricalPlanToCanvas(
          src_plan,
          projector,
          kf0.rgb,
          kf1.rgb,
          morph_scratch
        );
        morph_y0 = y0;
        morph_y1 = y1;
      }

      blend_buf.set(kf0.rgb);
      SpatialTransitionMorph.applyCanvasMorph(kf0.rgb, kf1.rgb, canvas_morph, t, blend_buf);
    }

    blitRgbaOntoRgb(blend_buf, frame_w, frame_h, yr_patch.data, yr_patch.width, yr_patch.height, yr_patch.x, yr_patch.y);

    let repeats = duplicate_frames ? 2 : 1;

    for (let x = 0; x < repeats; x++) {
      let ok = ff.stdin.write(Buffer.from(blend_buf.buffer, blend_buf.byteOffset, blend_buf.byteLength));

      if (!ok) {
        await new Promise((resolve) => ff.stdin.once("drain", resolve));
      }
    }
  }

  ff.stdin.end();

  await new Promise((resolve, reject) => {
    ff.on("close", (arg0_code) => {
      if (arg0_code === 0) resolve();
      else reject(new Error(`ffmpeg segment ${segment_id} exited ${arg0_code}: ${ff_err.slice(-600)}`));
    });
  });

  parentPort.postMessage({ segment_id: segment_id, frames: frame_count });
}


/**
 * Encodes the full 12 026-frame 4K MP4 video for a single target using parallel worker segments
 * matching landarea's exact 6 BCE (60 FPS) + 2 CE (30 FPS) VFR concat schedule.
 *
 * @param {string} arg0_target_key - Target key.
 * @param {object} arg1_render_cfg - Global render_config.
 */
async function encodeTargetVideo (arg0_target_key, arg1_render_cfg) {
  //Convert from parameters
  let target_key = arg0_target_key;
  let render_cfg = arg1_render_cfg;

  //Declare local instance variables
  let anim = render_cfg.animation;
  let cat_root = render_cfg.categorical_animations;
  let target_cfg = cat_root[target_key];
  let paths = GeoGridUtils.resolvePaths(target_cfg.indicator);
  let timesteps = GeoGridUtils.getTimestepYears(GeoGridUtils.loadConfig("model_config.json5"));
  let seg_dir = path.join(paths.renders_dir, `segments_${target_key}`);
  let out_mp4 = path.join(paths.renders_dir, target_cfg.output_filename);
  let switch_yr = anim.fps_switch_year;
  let jobs = [];
  let t0 = Date.now();

  fs.mkdirSync(paths.renders_dir, { recursive: true });
  fs.mkdirSync(seg_dir, { recursive: true });

  //6-tier multi-FPS schedule across parallel workers
  let schedule = anim.schedule || [];
  let total_planned_frames = 0;

  for (let i = 0; i < schedule.length; i++) {
    let tier = schedule[i];
    let step = tier.step_years || 1;
    let tier_years = [];

    for (let yr = tier.start_year; yr <= tier.end_year; yr += step) {
      tier_years.push(yr);
    }

    total_planned_frames += tier_years.length;

    let chunks = tier.chunks || 1;
    let chunk_size = Math.ceil(tier_years.length/chunks);

    for (let c = 0; c < chunks; c++) {
      let slice = tier_years.slice(c*chunk_size, (c + 1)*chunk_size);

      if (slice.length === 0) continue;
      jobs.push({
        segment_id: `${tier.id}_${c}`,
        start_year: slice[0],
        end_year: slice[slice.length - 1],
        step_years: step,
        fps: tier.fps,
        duplicate_frames: false,
        output_mp4: path.join(seg_dir, `seg_${String(jobs.length).padStart(2, "0")}.mp4`),
      });
    }
  }

  console.log(`[Encode:${target_key}] Starting ${total_planned_frames} frames (${jobs.length} workers) -> ${out_mp4}`);

  await Promise.all(jobs.map((job) => new Promise((resolve, reject) => {
    let w = new Worker(__filename, {
      workerData: {
        target_key: target_key,
        segment_id: job.segment_id,
        start_year: job.start_year,
        end_year: job.end_year,
        step_years: job.step_years,
        fps: job.fps,
        duplicate_frames: job.duplicate_frames,
        output_mp4: job.output_mp4,
        render_cfg: render_cfg,
        timesteps: timesteps,
      },
    });

    w.on("message", (arg0_msg) => {
      console.log(`[Encode:${target_key}] Segment ${arg0_msg.segment_id} done (${arg0_msg.frames} frames).`);
    });
    w.on("error", reject);
    w.on("exit", (arg0_code) => {
      if (arg0_code === 0) resolve();
      else reject(new Error(`Worker ${job.segment_id} exited with code ${arg0_code}`));
    });
  })));

  let concat_list_path = path.join(seg_dir, "concat.txt");
  let concat_lines = jobs.map((j) => `file '${j.output_mp4.replace(/\\/g, "/")}'`);

  fs.writeFileSync(concat_list_path, concat_lines.join("\n"));

  await new Promise((resolve, reject) => {
    let ff = spawn(ffmpeg_path, [
      "-y",
      "-f", "concat",
      "-safe", "0",
      "-i", concat_list_path,
      "-c", "copy",
      "-movflags", "+faststart",
      out_mp4,
    ], { stdio: ["ignore", "ignore", "pipe"] });
    let err = "";

    ff.stderr.on("data", (arg0_d) => { err += arg0_d.toString(); });
    ff.on("close", (arg0_code) => {
      if (arg0_code === 0) resolve();
      else reject(new Error(`ffmpeg concat failed (${arg0_code}): ${err.slice(-500)}`));
    });
  });

  fs.rmSync(seg_dir, { recursive: true, force: true });

  let elapsed_s = ((Date.now() - t0)/1000).toFixed(1);
  let size_mb = (fs.statSync(out_mp4).size/(1024*1024)).toFixed(2);

  console.log(`[Encode:${target_key}] Complete: ${out_mp4} (${size_mb} MB, ${elapsed_s}s)`);
}


/**
 * Main CLI entry point.
 */
async function main () {
  //Declare local instance variables
  let args = process.argv.slice(2);
  let allow_violations = false;
  let preview_years = null;
  let target_arg = "all";

  for (let i = 0; i < args.length; i++) {
    if (args[i].startsWith("--target=")) target_arg = args[i].split("=")[1];
    else if (args[i] === "--allow-violations") allow_violations = true;
    else if (args[i] === "--preview-only") preview_years = [-10000, -6000, 0, 1950, 2025];
    else if (args[i].startsWith("--preview-years=")) {
      preview_years = args[i].split("=")[1].split(",").map(Number);
    } else if (args[i] === "--preview-years" && i + 1 < args.length) {
      preview_years = args[i + 1].split(",").map(Number);
      i++;
    }
  }

  let render_cfg = GeoGridUtils.loadConfig("render_config.json5");
  let all_targets = ["trewartha", "koeppen_geiger", "holdridge_lifezones", "anthromes"];
  let targets = (target_arg === "all") ? all_targets : [target_arg];

  if (preview_years) {
    let canvas_cfg = Object.assign({}, render_cfg.canvas, { bilinear_sampling: false });
    let projector = new EqualEarthProjector(canvas_cfg);
    let halo_mask = buildDilatedGlobeMask(
      projector,
      render_cfg.categorical_animations.border_clearance_px || 24
    );

    for (let i = 0; i < targets.length; i++) {
      await renderTargetPreviews(targets[i], preview_years, render_cfg, projector, halo_mask, allow_violations);
    }
  } else {
    for (let i = 0; i < targets.length; i++) {
      await encodeTargetVideo(targets[i], render_cfg);
    }
  }
}


if (isMainThread) {
  main().catch((arg0_err) => {
    console.error(arg0_err);
    process.exit(1);
  });
} else {
  runSegmentWorker().catch((arg0_err) => {
    console.error(arg0_err);
    process.exit(1);
  });
}
