"""End-to-end orchestrator for the anthromes indicator pipeline (--only= / --skip=)."""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import A_acquire_sources
import B_build_dynamic_biomes_and_masks
import C_classify_anthromes


STEPS = [
  ("A", "Acquire sources", A_acquire_sources.main),
  ("B", "Build dynamic biomes and masks", B_build_dynamic_biomes_and_masks.main),
  ("C", "Classify anthromes & write outputs", C_classify_anthromes.main),
]


def parse_args() -> argparse.Namespace:
  parser = argparse.ArgumentParser(description="Run the anthromes indicator pipeline.")
  parser.add_argument("--only", type=str, default="", help="Comma-separated step letters to run (e.g. B,C)")
  parser.add_argument("--skip", type=str, default="", help="Comma-separated step letters to skip (e.g. A)")
  return parser.parse_args()


def main() -> None:
  args = parse_args()
  only_set = {s.strip().upper() for s in args.only.split(",") if s.strip()}
  skip_set = {s.strip().upper() for s in args.skip.split(",") if s.strip()}

  t0 = time.time()
  for code, label, fn in STEPS:
    if only_set and code not in only_set:
      continue
    if code in skip_set:
      print(f"[SKIP] Step {code}: {label}")
      continue
    print(f"\n>>> Running Step {code}: {label}")
    fn()

  print(f"\nAnthromes pipeline finished in {time.time() - t0:.1f}s")


if __name__ == "__main__":
  main()
