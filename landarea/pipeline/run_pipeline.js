"use strict";

let A_AcquireSources = require("./A_acquire_sources");
let B_BuildModernAnchor = require("./B_build_modern_anchor");
let C_BuildPaleoGiaIce = require("./C_build_paleo_gia_ice");
let D_BuildDynamicOverlays = require("./D_build_dynamic_overlays");
let E_AssembleOutputs = require("./E_assemble_outputs");
let F_AuditOpticalQA = require("./F_audit_optical_qa");
let G_RenderAnimation = require("./G_render_animation");

/**
 * End-to-end orchestrator for the landarea_Hollinger Node.js pipeline.
 * Supports selective step execution via --only=A,B,C or --skip=G.
 */
class HollingerPipeline {
  /**
   * Executes pipeline steps A through G in strict dependency order.
   *
   * @returns {Promise<void>}
   */
  static async run () {
    //Declare local instance variables
    let only_arg = process.argv.find((a) => a.startsWith("--only="));
    let skip_arg = process.argv.find((a) => a.startsWith("--skip="));
    let only_set = only_arg ? new Set(only_arg.slice(7).toUpperCase().split(",")) : null;
    let skip_set = skip_arg ? new Set(skip_arg.slice(7).toUpperCase().split(",")) : new Set();
    let steps = [
      ["A", "Acquire Sources", () => A_AcquireSources.run()],
      ["B", "Build Modern 2025 Anchor", () => B_BuildModernAnchor.run()],
      ["C", "Build Paleo GIA & Deglaciation", () => C_BuildPaleoGiaIce.run()],
      ["D", "Build Dynamic Overlays (Reclamation, Reservoirs, Paleolakes)", () => D_BuildDynamicOverlays.run()],
      ["E", "Assemble NetCDF-4 & 218 GeoPNGs", () => E_AssembleOutputs.run()],
      ["F", "Run Optical & Quantitative QA (vs HYDE & GHS POP)", () => F_AuditOpticalQA.run()],
      ["G", "Render 4K Equal Earth QBAM Animation (12 026 Frames)", () => G_RenderAnimation.run()],
    ];
    let t0 = Date.now();

    for (let i = 0; i < steps.length; i++) {
      let [code, name, fn] = steps[i];
      if (only_set && !only_set.has(code)) continue;
      if (skip_set.has(code)) continue;

      let ts = Date.now();
      console.log(`\n========== STEP ${code}: ${name} ==========`);
      await fn();
      console.log(`========== STEP ${code} COMPLETE (${((Date.now() - ts)/1000).toFixed(1)}s) ==========\n`);
    }

    console.log(`[HollingerPipeline] All requested steps completed in ${((Date.now() - t0)/1000).toFixed(1)}s.`);
  }
}

if (require.main === module) {
  HollingerPipeline.run().catch((err) => {
    console.error("[HollingerPipeline] Fatal error:", err);
    process.exit(1);
  });
}

module.exports = HollingerPipeline;
