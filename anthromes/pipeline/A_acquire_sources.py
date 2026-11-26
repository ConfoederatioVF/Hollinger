"""Step A: Acquire Anthromes 12K DGG v1 sources and verify local HYDE 3.5, landarea, and climate inputs."""

from __future__ import annotations

import shutil
import subprocess
import sys
import time
import zipfile
from pathlib import Path

import netCDF4 as nc
import requests

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "core"))
import hollinger_core as hc

if hasattr(sys.stdout, "reconfigure"):
  sys.stdout.reconfigure(line_buffering=True)


BROWSER_UA = (
  "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
  "AppleWebKit/537.36 (KHTML, like Gecko) "
  "Chrome/120.0.0.0 Safari/537.36"
)
BROWSER_HEADERS = {"User-Agent": BROWSER_UA}


def download_file_resumable(
  url: str,
  dest_path: Path,
  retries: int = 5,
  timeout: int = 120,
) -> Path:
  """Downloads a file with HTTP Range resume, browser User-Agent, and retry backoff."""
  dest_path.parent.mkdir(parents=True, exist_ok=True)
  part_path = dest_path.with_suffix(dest_path.suffix + ".part")

  if dest_path.exists() and dest_path.stat().st_size > 0:
    return dest_path

  curl_bin = shutil.which("curl.exe") or shutil.which("curl")
  if curl_bin is not None:
    cmd = [
      curl_bin,
      "-f",
      "-L",
      "-s",
      "-S",
      "-C",
      "-",
      "--retry",
      str(retries),
      "--connect-timeout",
      "30",
      "-A",
      BROWSER_UA,
      "-o",
      str(part_path),
      url,
    ]
    subprocess.run(cmd, check=True)
    if part_path.exists() and part_path.stat().st_size > 0:
      part_path.replace(dest_path)
      return dest_path

  for attempt in range(1, retries + 1):
    try:
      existing = part_path.stat().st_size if part_path.exists() else 0
      headers = dict(BROWSER_HEADERS)
      mode = "wb"
      if existing > 0:
        headers["Range"] = f"bytes={existing}-"
        mode = "ab"

      with requests.get(url, headers=headers, stream=True, timeout=timeout, allow_redirects=True) as resp:
        if resp.status_code == 416:
          if part_path.exists():
            part_path.unlink()
          existing = 0
          mode = "wb"
          resp.close()
          continue
        if resp.status_code == 200 and existing > 0:
          mode = "wb"
          existing = 0
        elif resp.status_code not in (200, 206):
          resp.raise_for_status()

        total_len = resp.headers.get("Content-Length")
        expected_total = existing + int(total_len) if total_len is not None else None
        with open(part_path, mode) as f:
          for chunk in resp.iter_content(chunk_size=1024 * 1024):
            if chunk:
              f.write(chunk)

        if expected_total is not None and part_path.stat().st_size < expected_total:
          raise IOError(f"Incomplete download ({part_path.stat().st_size} < {expected_total})")

      part_path.replace(dest_path)
      return dest_path
    except Exception as exc:
      print(f"  [WARN] Attempt {attempt}/{retries} failed for {dest_path.name}: {exc}")
      if attempt == retries:
        raise
      time.sleep(min(30, 2 ** attempt))

  return dest_path


def verify_netcdf_file(nc_path: Path, expected_var: str) -> None:
  """Verifies a NetCDF file exists, opens cleanly, and contains the expected variable."""
  if not nc_path.exists():
    raise FileNotFoundError(f"Missing required NetCDF: {nc_path}")
  with nc.Dataset(str(nc_path), "r") as ds:
    if expected_var not in ds.variables:
      raise KeyError(f"Variable '{expected_var}' not found in {nc_path}")
    v = ds.variables[expected_var]
    if v.ndim != 3 or v.shape[1:] != (2160, 4320):
      raise ValueError(f"Unexpected shape {v.shape} for '{expected_var}' in {nc_path}")
    _ = v[0, 0, 0]
    _ = v[-1, -1, -1]


def extract_supporting_grids(raw_zip_path: Path, out_dir: Path) -> None:
  """Extracts supporting_5m_grids.zip from raw-data.zip and unpacks target GeoTIFFs."""
  out_dir.mkdir(parents=True, exist_ok=True)
  target_files = [
    "potveg15.tif",
    "potvill20.tif",
    "maxln_cr.tif",
    "simple_regions.tif",
    "iso_cr.tif",
  ]
  if all((out_dir / fname).exists() and (out_dir / fname).stat().st_size > 0 for fname in target_files):
    print(f"  [OK] Supporting 5m grids already extracted in {out_dir}")
    return

  inner_zip_path = raw_zip_path.parent / "supporting_5m_grids.zip"
  if not inner_zip_path.exists() or inner_zip_path.stat().st_size == 0:
    with zipfile.ZipFile(raw_zip_path, "r") as zf:
      matches = [m for m in zf.namelist() if m.endswith("supporting_5m_grids.zip")]
      if not matches:
        raise FileNotFoundError(f"supporting_5m_grids.zip not found inside {raw_zip_path}")
      print(f"  Extracting {matches[0]} from {raw_zip_path.name}...")
      with zf.open(matches[0]) as src, open(inner_zip_path, "wb") as dst:
        while True:
          buf = src.read(4 * 1024 * 1024)
          if not buf:
            break
          dst.write(buf)

  with zipfile.ZipFile(inner_zip_path, "r") as zf:
    names = zf.namelist()
    for fname in target_files:
      matches = [m for m in names if m.endswith("/" + fname) or m == fname]
      if not matches:
        raise FileNotFoundError(f"{fname} not found inside {inner_zip_path}")
      dest = out_dir / fname
      if not dest.exists() or dest.stat().st_size == 0:
        print(f"  Extracting {fname} -> {dest}")
        with zf.open(matches[0]) as src, open(dest, "wb") as dst:
          while True:
            buf = src.read(4 * 1024 * 1024)
            if not buf:
              break
            dst.write(buf)


def main() -> None:
  cfg = hc.load_json5("anthromes_config.json5")
  paths = hc.resolve_paths("anthromes")
  raw_dir = paths["raw_dir"]
  raw_dir.mkdir(parents=True, exist_ok=True)
  base_dir = paths["base_dir"]
  src_cfg = cfg["sources"]

  # 1. Download Anthromes 12K DGG v1 sources
  downloads = [
    (src_cfg["anthromes12k_raw_data_url"], raw_dir / "raw-data.zip"),
    (src_cfg["anthromes12k_flowchart_url"], raw_dir / "anthromes_12K_DGG_v1_flowchart.pdf"),
    (src_cfg["anthromes12k_readme_url"], raw_dir / "README.txt"),
  ]
  for url, dest in downloads:
    print(f"Checking/downloading {dest.name} from {url}...")
    download_file_resumable(url, dest)
    print(f"  [OK] {dest.name} ({dest.stat().st_size:,} bytes)")

  # 2. Extract supporting_5m_grids.zip from raw-data.zip
  grids_dir = raw_dir / "supporting_5m_grids"
  extract_supporting_grids(raw_dir / "raw-data.zip", grids_dir)

  # 3. Verify sibling HYDE 3.5 NetCDFs
  hyde_dir = (base_dir / src_cfg["hyde35_rel_dir"]).resolve()
  hyde_vars = [
    ("cropland.nc", "cropland"),
    ("grazing_land.nc", "grazing_land"),
    ("urban_area.nc", "urban_area"),
    ("total_rice.nc", "total_rice"),
    ("total_irrigated.nc", "total_irrigated"),
    ("population.nc", "population"),
  ]
  for fname, var in hyde_vars:
    verify_netcdf_file(hyde_dir / fname, var)
  print(f"  [OK] Verified {len(hyde_vars)} HYDE 3.5 NetCDFs in {hyde_dir}")

  # 4. Verify HYDE land_area_per_cell.asc
  hyde_meta_dir = (base_dir / src_cfg["hyde_meta_rel_dir"]).resolve()
  asc_candidates = [
    hyde_meta_dir / "land_area_per_cell.asc",
    hyde_meta_dir.parent / "land_area_per_cell.asc",
  ]
  asc_path = next((p for p in asc_candidates if p.exists()), None)
  if asc_path is None:
    raise FileNotFoundError(f"land_area_per_cell.asc not found in {[str(p) for p in asc_candidates]}")
  print(f"  [OK] Verified HYDE land_area_per_cell.asc at {asc_path}")

  # 5. Verify landarea/landarea.nc and climate/holdridge_lifezones.nc
  verify_netcdf_file(base_dir / "landarea" / "landarea.nc", "land_area")
  verify_netcdf_file(base_dir / "climate" / "holdridge_lifezones.nc", "holdridge_lifezones")
  print("  [OK] Verified landarea/landarea.nc and climate/holdridge_lifezones.nc")


if __name__ == "__main__":
  main()
