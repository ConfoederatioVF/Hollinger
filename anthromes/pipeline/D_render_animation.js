/**
 * @fileoverview Pipeline Step D: Renders 4K Equal Earth export video for Anthromes.
 */

const path = require("path");
const { spawnSync } = require("child_process");


/**
 * Invokes CategoricalAnimationRenderer for Anthromes, forwarding CLI flags.
 */
function run () {
  //Declare local instance variables
  let extra_args = process.argv.slice(2);
  let renderer_path = path.resolve(__dirname, "..", "..", "core", "CategoricalAnimationRenderer.js");
  let res = spawnSync(
    process.execPath,
    ["--max-old-space-size=16384", renderer_path, "--target=anthromes", ...extra_args],
    { stdio: "inherit" }
  );

  if (res.status !== 0) process.exit(res.status || 1);
}


run();

