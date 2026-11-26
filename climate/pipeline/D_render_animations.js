/**
 * @fileoverview Pipeline Step D: Renders 4K Equal Earth export videos for Trewartha,
 * Köppen-Geiger, and Holdridge Life Zones.
 */

const path = require("path");
const { spawnSync } = require("child_process");


/**
 * Invokes CategoricalAnimationRenderer for all three climate targets, forwarding CLI flags.
 */
function run () {
  //Declare local instance variables
  let extra_args = process.argv.slice(2);
  let renderer_path = path.resolve(__dirname, "..", "..", "core", "CategoricalAnimationRenderer.js");
  let targets = ["trewartha", "koeppen_geiger", "holdridge_lifezones"];

  for (let i = 0; i < targets.length; i++) {
    let res = spawnSync(
      process.execPath,
      ["--max-old-space-size=16384", renderer_path, `--target=${targets[i]}`, ...extra_args],
      { stdio: "inherit" }
    );

    if (res.status !== 0) process.exit(res.status || 1);
  }
}


run();

