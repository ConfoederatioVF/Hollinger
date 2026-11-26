"""Step A: Resumable, parallel source acquisition for the climate indicator."""

from __future__ import annotations

import concurrent.futures
import gzip
import os
import re
import shutil
import sys
import time
import zipfile
from pathlib import Path

import netCDF4 as nc
import requests
import tifffile

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "core"))
import hollinger_core as hc


def verify_netcdf_file(nc_path: Path, expected_var: str) -> bool:
  """Verifies a NetCDF file can be opened and has the expected variable."""
  if not nc_path.exists() or nc_path.stat().st_size < 1024:
    return False
  try:
    with nc.Dataset(str(nc_path), "r") as ds:
      if expected_var not in ds.variables:
        return False
      v = ds.variables[expected_var]
      if v.ndim < 2 or v.shape[0] == 0:
        return False
      # Touch first and last slice to ensure file is not truncated
      _ = v[0, 0, 0]
      _ = v[-1, -1, -1]
    return True
  except Exception:
    return False


def download_file_resumable(
  url: str,
  dest_path: Path,
  expected_size: int | None = None,
  retries: int = 5,
  timeout: int = 120,
) -> Path:
  """Downloads a file with HTTP Range resume and retry backoff."""
  dest_path.parent.mkdir(parents=True, exist_ok=True)
  part_path = dest_path.with_suffix(dest_path.suffix + ".part")

  if dest_path.exists() and (expected_size is None or dest_path.stat().st_size == expected_size):
    return dest_path

  for attempt in range(1, retries + 1):
    try:
      existing = part_path.stat().st_size if part_path.exists() else 0
      headers = {}
      mode = "wb"
      if existing > 0:
        headers["Range"] = f"bytes={existing}-"
        mode = "ab"

      with requests.get(url, headers=headers, stream=True, timeout=timeout, allow_redirects=True) as resp:
        if resp.status_code == 416:
          # Range not satisfiable; restart
          if part_path.exists():
            part_path.unlink()
          existing = 0
          mode = "wb"
          resp.close()
          continue
        if resp.status_code == 200 and existing > 0:
          # Server ignored Range; overwrite from 0
          mode = "wb"
          existing = 0
        elif resp.status_code not in (200, 206):
          resp.raise_for_status()

        total_len = resp.headers.get("Content-Length")
        with open(part_path, mode) as f:
          written = 0
          for chunk in resp.iter_content(chunk_size=256 * 1024):
            if chunk:
              f.write(chunk)
              written += len(chunk)
              if written >= 4 * 1024 * 1024:
                f.flush()
                written = 0

      if expected_size is not None and part_path.stat().st_size != expected_size:
        raise IOError(f"Size mismatch for {dest_path.name}: got {part_path.stat().st_size}, expected {expected_size}")

      if dest_path.exists():
        dest_path.unlink()
      part_path.replace(dest_path)
      return dest_path
    except Exception as exc:
      if attempt == retries:
        raise RuntimeError(f"Failed downloading {url} after {retries} attempts: {exc}") from exc
      time.sleep(min(30, 2 ** attempt))

  return dest_path


def fetch_trace21k_manifest(cfg: dict) -> list[dict]:
  """Fetches TraCE-21ka file manifest from GDEX API and builds download tasks for TREFHT, PRECC, PRECL, PHIS."""
  src = cfg["sources"]
  api_url = src["gdex_filelist_api"]
  base_url = src["trace21k_base_url"].rstrip("/")
  min_seg = int(src["trace21k_min_segment"])
  max_seg = int(src["trace21k_max_segment"])
  vars_needed = list(src["trace21k_vars"])

  resp = requests.get(api_url, timeout=60)
  resp.raise_for_status()
  payload = resp.json()
  rows = payload["data"]["groups"][0]["rows"]

  phis_by_seg = {}
  for row in rows:
    finfo = row[0]
    fname = finfo["value"]
    size = int(row[1]["value"])
    m = re.match(r"trace\.(\d+)\.(.+?)\.cam2\.h0\.([A-Z0-9_]+)\.(\d+-\d+)\.nc$", fname)
    if not m:
      continue
    seg = int(m.group(1))
    period = m.group(2)
    var = m.group(3)
    dates = m.group(4)
    if min_seg <= seg <= max_seg and var == "PHIS":
      phis_by_seg[seg] = {
        "seg": seg,
        "period": period,
        "dates": dates,
        "phis_size": size,
        "data_path": finfo["data_path"],
        "osdf_url": finfo["url"],
      }

  tasks = []
  for seg in range(min_seg, max_seg + 1):
    if seg not in phis_by_seg:
      raise KeyError(f"Missing segment {seg} in GDEX filelist")
    meta = phis_by_seg[seg]
    for var in vars_needed:
      fname = f"trace.{seg:02d}.{meta['period']}.cam2.h0.{var}.{meta['dates']}.nc"
      rel_path = f"TraCE/TraCE-Main/atm/proc/tseries/monthly/{var}/{fname}"
      primary_url = f"{base_url}/{rel_path}"
      fallback_url = f"https://osdf-director.osg-htc.org/ncar/gdex/d651050/{rel_path}"
      tasks.append({
        "seg": seg,
        "var": var,
        "fname": fname,
        "primary_url": primary_url,
        "fallback_url": fallback_url,
      })
  return tasks


def download_trace_task(task: dict, out_dir: Path) -> str:
  dest = out_dir / task["var"] / task["fname"]
  if verify_netcdf_file(dest, task["var"]):
    return f"SKIP {task['fname']}"

  if dest.exists():
    dest.unlink()

  try:
    download_file_resumable(task["primary_url"], dest)
  except Exception:
    download_file_resumable(task["fallback_url"], dest)

  if not verify_netcdf_file(dest, task["var"]):
    dest.unlink(missing_ok=True)
    raise RuntimeError(f"NetCDF integrity check failed for {dest}")
  return f"OK   {task['fname']}"


def acquire_worldclim(cfg: dict, wc_dir: Path) -> None:
  """Downloads and extracts WorldClim 2.1 5-arcminute tavg, prec, and elev GeoTIFFs."""
  wc_dir.mkdir(parents=True, exist_ok=True)
  src = cfg["sources"]
  items = [
    ("tavg", src["worldclim_tavg_url"], [f"wc2.1_5m_tavg_{m:02d}.tif" for m in range(1, 13)]),
    ("prec", src["worldclim_prec_url"], [f"wc2.1_5m_prec_{m:02d}.tif" for m in range(1, 13)]),
    ("elev", src["worldclim_elev_url"], ["wc2.1_5m_elev.tif"]),
  ]
  for key, url, expected_files in items:
    all_ok = all((wc_dir / f).exists() and (wc_dir / f).stat().st_size > 1000 for f in expected_files)
    if all_ok:
      print(f"[WorldClim] {key}: already extracted ({len(expected_files)} files)")
      continue
    zip_path = wc_dir / Path(url).name
    print(f"[WorldClim] Downloading {zip_path.name}...")
    download_file_resumable(url, zip_path)
    print(f"[WorldClim] Extracting {zip_path.name}...")
    with zipfile.ZipFile(zip_path, "r") as zf:
      zf.extractall(wc_dir)
    for f in expected_files:
      fpath = wc_dir / f
      if not fpath.exists():
        raise FileNotFoundError(f"Expected {fpath} after extracting {zip_path.name}")
      with tifffile.TiffFile(fpath) as tf:
        _ = tf.pages[0].shape
    print(f"[WorldClim] {key}: verified {len(expected_files)} GeoTIFFs")


def acquire_cru_one(var: str, url: str, cru_dir: Path) -> str:
  """Downloads and gunzips a single CRU TS 4.10 variable NetCDF file."""
  cru_dir.mkdir(parents=True, exist_ok=True)
  gz_name = Path(url).name
  nc_name = gz_name[:-3] if gz_name.endswith(".gz") else gz_name
  gz_path = cru_dir / gz_name
  nc_path = cru_dir / nc_name
  if verify_netcdf_file(nc_path, var):
    return f"SKIP {nc_name}"
  if nc_path.exists():
    nc_path.unlink()
  print(f"[CRU] Downloading {gz_name}...", flush=True)
  download_file_resumable(url, gz_path)
  print(f"[CRU] Gunzipping {gz_name} -> {nc_name}...", flush=True)
  part_nc = nc_path.with_suffix(".nc.part")
  with gzip.open(gz_path, "rb") as f_in, open(part_nc, "wb") as f_out:
    shutil.copyfileobj(f_in, f_out, length=4 * 1024 * 1024)
  part_nc.replace(nc_path)
  if not verify_netcdf_file(nc_path, var):
    nc_path.unlink(missing_ok=True)
    raise RuntimeError(f"CRU NetCDF verification failed for {nc_path}")
  return f"OK   {nc_name}"


def main() -> None:
  paths = hc.resolve_paths("climate")
  cfg = hc.load_json5("climate_config.json5")
  mcfg = hc.load_json5("model_config.json5")
  max_workers = min(8, int(mcfg["compute"]["max_workers"]))

  raw_dir = paths["raw_dir"]
  trace_dir = raw_dir / "trace21k"
  wc_dir = raw_dir / "worldclim"
  cru_dir = raw_dir / "cru"

  print("=== Step A: Acquiring WorldClim 2.1 ===", flush=True)
  acquire_worldclim(cfg, wc_dir)

  print("=== Step A: Fetching TraCE-21ka manifest ===", flush=True)
  tasks = fetch_trace21k_manifest(cfg)
  src = cfg["sources"]
  cru_items = [("tmp", src["cru_tmp_url"]), ("pre", src["cru_pre_url"])]
  total_jobs = len(tasks) + len(cru_items)
  print(f"[Parallel] {total_jobs} jobs queued (CRU=2, TraCE-21ka={len(tasks)}), workers={max_workers}", flush=True)

  done = 0
  with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as pool:
    futs = [pool.submit(acquire_cru_one, var, url, cru_dir) for var, url in cru_items]
    futs.extend(pool.submit(download_trace_task, t, trace_dir) for t in tasks)
    for fut in concurrent.futures.as_completed(futs):
      msg = fut.result()
      done += 1
      if done % 5 == 0 or done == total_jobs or msg.startswith("OK"):
        print(f"  [{done}/{total_jobs}] {msg}", flush=True)

  print("=== Step A complete ===", flush=True)


if __name__ == "__main__":
  main()
