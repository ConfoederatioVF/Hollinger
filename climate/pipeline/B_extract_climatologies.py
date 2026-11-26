"""Step B: Extract CCSM3 TraCE-21ka, CRU TS 4.10, and gap-free WorldClim 2.1 5' climatologies."""

from __future__ import annotations

import re
import sys
from pathlib import Path

import netCDF4 as nc
import numpy as np
from PIL import Image
from scipy.ndimage import distance_transform_edt

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "core"))
import hollinger_core as hc

Image.MAX_IMAGE_PIXELS = None

DAYS_IN_MONTH = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31], dtype=np.float32)


def parse_trace_segment_years(fname: str) -> tuple[int, int, int]:
  """Parses segment number and signed calendar start/end years (CE, 0 = 1 AD) from TraCE filename."""
  m = re.match(r"trace\.(\d+)\.(\d+)-(\d+)BP\.cam2\.h0\.", fname)
  if m:
    seg = int(m.group(1))
    bp_start = int(m.group(2))
    bp_end = int(m.group(3))
    return seg, 1950 - bp_start, 1950 - bp_end

  m_ce = re.match(r"trace\.(\d+)\.(\d+)BP-(\d+)CE\.cam2\.h0\.", fname)
  if m_ce:
    seg = int(m_ce.group(1))
    bp_start = int(m_ce.group(2))
    ce_end = int(m_ce.group(3))
    return seg, 1950 - bp_start, ce_end - 1

  raise ValueError(f"Unrecognised TraCE-21ka filename format: {fname}")


def inpaint_periodic_2d(arr: np.ndarray, valid_mask: np.ndarray, n_smooth: int = 8) -> np.ndarray:
  """Inpaints invalid cells using lon-periodic nearest-neighbour + Laplace relaxation."""
  h, w = arr.shape
  if np.all(valid_mask):
    return arr.astype(np.float32, copy=True)

  # Tile 3x in longitude for periodic distance transform
  pad_w = w // 4
  ext_arr = np.concatenate([arr[:, -pad_w:], arr, arr[:, :pad_w]], axis=1)
  ext_mask = np.concatenate([valid_mask[:, -pad_w:], valid_mask, valid_mask[:, :pad_w]], axis=1)

  _, indices = distance_transform_edt(~ext_mask, return_indices=True)
  filled_ext = ext_arr[indices[0], indices[1]]
  out = filled_ext[:, pad_w : pad_w + w].astype(np.float32, copy=True)

  miss = ~valid_mask
  for _ in range(n_smooth):
    left = np.roll(out, 1, axis=1)
    right = np.roll(out, -1, axis=1)
    up = np.empty_like(out)
    up[0, :] = out[0, :]
    up[1:, :] = out[:-1, :]
    down = np.empty_like(out)
    down[-1, :] = out[-1, :]
    down[:-1, :] = out[1:, :]
    avg = 0.25 * (left + right + up + down)
    out[miss] = avg[miss]

  return out


def bilinear_upsample_periodic(
  src: np.ndarray,
  src_lat: np.ndarray,
  src_lon: np.ndarray,
  dst_lat: np.ndarray,
  dst_lon: np.ndarray,
) -> np.ndarray:
  """Lon-periodic bilinear interpolation from (..., nlat, nlon) N->S grid to (..., H, W) N->S grid."""
  # src_lat is descending (N->S); flip for np.interp / searchsorted ascending
  if src_lat[0] > src_lat[-1]:
    lat_asc = src_lat[::-1]
    src_use = src[..., ::-1, :]
  else:
    lat_asc = src_lat
    src_use = src

  nlat = len(lat_asc)
  nlon = len(src_lon)

  # Latitude indices & weights (clamped at poles)
  y_pos = np.interp(dst_lat, lat_asc, np.arange(nlat, dtype=np.float64))
  y0 = np.clip(np.floor(y_pos).astype(np.int32), 0, nlat - 1)
  y1 = np.clip(y0 + 1, 0, nlat - 1)
  wy = (y_pos - y0).astype(np.float32)

  # Longitude indices & weights (periodic wrap around -180/+180)
  lon_ext = np.concatenate([src_lon[-1:] - 360.0, src_lon, src_lon[:1] + 360.0])
  x_pos = np.interp(dst_lon, lon_ext, np.arange(nlon + 2, dtype=np.float64)) - 1.0
  x0 = np.floor(x_pos).astype(np.int32)
  wx = (x_pos - x0).astype(np.float32)
  x0_mod = np.mod(x0, nlon)
  x1_mod = np.mod(x0 + 1, nlon)

  # Separable 1D lon-then-lat bilinear (avoids 2D fancy indexing; ~10x faster)
  wy_b = wy[:, None]
  wx_b = wx[None, :]
  if src_use.ndim == 3:
    lon_interp = src_use[:, :, x0_mod] * (1.0 - wx_b[None, :, :]) + src_use[:, :, x1_mod] * wx_b[None, :, :]
    return (lon_interp[:, y0, :] * (1.0 - wy_b[None, :, :]) + lon_interp[:, y1, :] * wy_b[None, :, :]).astype(
      np.float32
    )

  lon_interp = src_use[..., x0_mod] * (1.0 - wx_b) + src_use[..., x1_mod] * wx_b
  return (lon_interp[..., y0, :] * (1.0 - wy_b) + lon_interp[..., y1, :] * wy_b).astype(np.float32)


def extract_ccsm3_climatologies(cfg: dict, years: list[int], raw_trace_dir: Path, out_npz: Path) -> dict:
  """Extracts centred monthly T, P, Z climatologies and 1960-1990 normal from TraCE-21ka segments 16..36."""
  ds_cfg = cfg["downscaling"]
  min_seg = int(cfg["sources"]["trace21k_min_segment"])
  max_seg = int(cfg["sources"]["trace21k_max_segment"])
  ancient_cutoff = int(ds_cfg["ancient_cutoff_year"])
  ancient_half = int(ds_cfg["ancient_window_half_yrs"])
  modern_half = int(ds_cfg["modern_window_half_yrs"])
  ref_start = int(ds_cfg["ccsm3_ref_start"])
  ref_end = int(ds_cfg["ccsm3_ref_end"])

  ccsm_years = [y for y in years if y <= 1990]
  n_targets = len(ccsm_years)
  windows = []
  for y in ccsm_years:
    hw = ancient_half if y <= ancient_cutoff else modern_half
    windows.append((y - hw, y + hw))

  # Collect segment file paths
  seg_files = {}
  for seg in range(min_seg, max_seg + 1):
    seg_files[seg] = {}
    for var in ("TREFHT", "PRECC", "PRECL", "PHIS"):
      matches = sorted((raw_trace_dir / var).glob(f"trace.{seg:02d}.*.nc"))
      if len(matches) != 1:
        raise FileNotFoundError(f"Expected 1 file for seg {seg} {var}, found {len(matches)}")
      seg_files[seg][var] = matches[0]

  # Inspect grid from first file
  with nc.Dataset(str(seg_files[min_seg]["TREFHT"]), "r") as ds0:
    raw_lat = np.asarray(ds0.variables["lat"][:], dtype=np.float64)
    raw_lon = np.asarray(ds0.variables["lon"][:], dtype=np.float64)

  flip_lat = bool(raw_lat[0] < raw_lat[-1])
  lat = raw_lat[::-1].copy() if flip_lat else raw_lat.copy()
  lon_wrapped = np.where(raw_lon >= 180.0, raw_lon - 360.0, raw_lon)
  lon_order = np.argsort(lon_wrapped)
  lon = lon_wrapped[lon_order].copy()
  nlat, nlon = len(lat), len(lon)

  t_sum = np.zeros((n_targets, 12, nlat, nlon), dtype=np.float64)
  p_sum = np.zeros((n_targets, 12, nlat, nlon), dtype=np.float64)
  z_sum = np.zeros((n_targets, nlat, nlon), dtype=np.float64)
  counts = np.zeros(n_targets, dtype=np.int32)

  t_ref_sum = np.zeros((12, nlat, nlon), dtype=np.float64)
  p_ref_sum = np.zeros((12, nlat, nlon), dtype=np.float64)
  z_ref_sum = np.zeros((nlat, nlon), dtype=np.float64)
  ref_count = 0

  prev_end_ce = None
  p_scale = (DAYS_IN_MONTH * 86400.0 * 1000.0)[None, :, None, None]

  for seg in range(min_seg, max_seg + 1):
    fname = seg_files[seg]["TREFHT"].name
    _, start_ce, end_ce = parse_trace_segment_years(fname)
    if prev_end_ce is not None and start_ce != prev_end_ce + 1:
      raise ValueError(f"Non-contiguous TraCE segment {seg}: {start_ce} != {prev_end_ce + 1}")
    prev_end_ce = end_ce
    n_seg_years = end_ce - start_ce + 1

    with (
      nc.Dataset(str(seg_files[seg]["TREFHT"]), "r") as ds_t,
      nc.Dataset(str(seg_files[seg]["PRECC"]), "r") as ds_pc,
      nc.Dataset(str(seg_files[seg]["PRECL"]), "r") as ds_pl,
      nc.Dataset(str(seg_files[seg]["PHIS"]), "r") as ds_z,
    ):
      n_months = ds_t.variables["TREFHT"].shape[0]
      if n_months != n_seg_years * 12:
        raise ValueError(f"Segment {seg} month count {n_months} != {n_seg_years * 12}")

      t_arr = np.asarray(ds_t.variables["TREFHT"][:], dtype=np.float32) - 273.15
      p_rate = np.asarray(ds_pc.variables["PRECC"][:], dtype=np.float32) + np.asarray(
        ds_pl.variables["PRECL"][:], dtype=np.float32
      )
      z_arr = np.asarray(ds_z.variables["PHIS"][:], dtype=np.float32) / 9.80665

    if flip_lat:
      t_arr = t_arr[:, ::-1, :]
      p_rate = p_rate[:, ::-1, :]
      z_arr = z_arr[:, ::-1, :]
    t_arr = t_arr[:, :, lon_order].reshape(n_seg_years, 12, nlat, nlon)
    p_arr = (p_rate[:, :, lon_order].reshape(n_seg_years, 12, nlat, nlon)) * p_scale
    z_ann = z_arr[:, :, lon_order].reshape(n_seg_years, 12, nlat, nlon).mean(axis=1)

    seg_years = np.arange(start_ce, end_ce + 1, dtype=np.int32)
    for ti, (w0, w1) in enumerate(windows):
      mask = (seg_years >= w0) & (seg_years <= w1)
      if np.any(mask):
        t_sum[ti] += t_arr[mask].sum(axis=0, dtype=np.float64)
        p_sum[ti] += p_arr[mask].sum(axis=0, dtype=np.float64)
        z_sum[ti] += z_ann[mask].sum(axis=0, dtype=np.float64)
        counts[ti] += int(mask.sum())

    ref_mask = (seg_years >= ref_start) & (seg_years <= ref_end)
    if np.any(ref_mask):
      t_ref_sum += t_arr[ref_mask].sum(axis=0, dtype=np.float64)
      p_ref_sum += p_arr[ref_mask].sum(axis=0, dtype=np.float64)
      z_ref_sum += z_ann[ref_mask].sum(axis=0, dtype=np.float64)
      ref_count += int(ref_mask.sum())

    print(f"  [CCSM3] seg {seg:02d}: {start_ce}..{end_ce} CE ({n_seg_years} yrs)")

  if np.any(counts == 0):
    missing = [ccsm_years[i] for i in range(n_targets) if counts[i] == 0]
    raise RuntimeError(f"Zero-count CCSM3 windows for years: {missing}")

  t_clim = (t_sum / counts[:, None, None, None]).astype(np.float32)
  p_clim = (p_sum / counts[:, None, None, None]).astype(np.float32)
  z_clim = (z_sum / counts[:, None, None]).astype(np.float32)
  t_ref = (t_ref_sum / ref_count).astype(np.float32)
  p_ref = (p_ref_sum / ref_count).astype(np.float32)
  z_ref = (z_ref_sum / ref_count).astype(np.float32)

  out_npz.parent.mkdir(parents=True, exist_ok=True)
  np.savez_compressed(
    out_npz,
    years=np.asarray(ccsm_years, dtype=np.int32),
    lat=lat.astype(np.float32),
    lon=lon.astype(np.float32),
    T_clim=t_clim,
    P_clim=p_clim,
    Z_clim=z_clim,
    T_ref=t_ref,
    P_ref=p_ref,
    Z_ref=z_ref,
  )
  print(f"[CCSM3] Saved {out_npz} ({n_targets} timesteps, {nlat}x{nlon})")
  return {
    "years": np.asarray(ccsm_years, dtype=np.int32),
    "lat": lat.astype(np.float32),
    "lon": lon.astype(np.float32),
    "T_clim": t_clim,
    "P_clim": p_clim,
    "Z_clim": z_clim,
    "T_ref": t_ref,
    "P_ref": p_ref,
    "Z_ref": z_ref,
  }


def extract_cru_climatologies(cfg: dict, years: list[int], raw_cru_dir: Path, out_npz: Path) -> dict:
  """Extracts 30-yr centred monthly CRU TS 4.10 climatologies (1900..2025) + 1970-2000 & 1960-1990 normals."""
  ds_cfg = cfg["downscaling"]
  hw = int(ds_cfg["modern_window_half_yrs"])
  wc_s, wc_e = int(ds_cfg["worldclim_ref_start"]), int(ds_cfg["worldclim_ref_end"])
  cc_s, cc_e = int(ds_cfg["ccsm3_ref_start"]), int(ds_cfg["ccsm3_ref_end"])

  cru_years = [y for y in years if y >= int(ds_cfg["cru_splice_start"])]
  n_targets = len(cru_years)

  tmp_nc = sorted(raw_cru_dir.glob("*tmp*.nc"))[0]
  pre_nc = sorted(raw_cru_dir.glob("*pre*.nc"))[0]

  with nc.Dataset(str(tmp_nc), "r") as ds_t, nc.Dataset(str(pre_nc), "r") as ds_p:
    raw_lat = np.asarray(ds_t.variables["lat"][:], dtype=np.float64)
    lon = np.asarray(ds_t.variables["lon"][:], dtype=np.float64)
    flip_lat = bool(raw_lat[0] < raw_lat[-1])
    lat = raw_lat[::-1].copy() if flip_lat else raw_lat.copy()
    nlat, nlon = len(lat), len(lon)

    n_months = ds_t.variables["tmp"].shape[0]
    n_cru_yrs = n_months // 12
    cru_all_yrs = np.arange(1901, 1901 + n_cru_yrs, dtype=np.int32)

    t_raw = np.ma.filled(ds_t.variables["tmp"][:], np.nan).astype(np.float32)
    p_raw = np.ma.filled(ds_p.variables["pre"][:], np.nan).astype(np.float32)

  if flip_lat:
    t_raw = t_raw[:, ::-1, :]
    p_raw = p_raw[:, ::-1, :]

  t_raw = t_raw.reshape(n_cru_yrs, 12, nlat, nlon)
  p_raw = p_raw.reshape(n_cru_yrs, 12, nlat, nlon)
  valid_land = np.isfinite(t_raw[0, 0]) & np.isfinite(p_raw[0, 0])

  def window_mean_and_fill(arr: np.ndarray, mask_yrs: np.ndarray) -> np.ndarray:
    with np.errstate(all="ignore"):
      sub = np.nanmean(arr[mask_yrs], axis=0).astype(np.float32)
    out = np.empty_like(sub)
    for m in range(12):
      v = valid_land & np.isfinite(sub[m])
      out[m] = inpaint_periodic_2d(np.where(v, sub[m], 0.0), v, n_smooth=6)
    return out

  t_clim = np.empty((n_targets, 12, nlat, nlon), dtype=np.float32)
  p_clim = np.empty((n_targets, 12, nlat, nlon), dtype=np.float32)
  for ti, y in enumerate(cru_years):
    ymask = (cru_all_yrs >= y - hw) & (cru_all_yrs <= y + hw)
    t_clim[ti] = window_mean_and_fill(t_raw, ymask)
    p_clim[ti] = np.maximum(0.0, window_mean_and_fill(p_raw, ymask))

  m_70_00 = (cru_all_yrs >= wc_s) & (cru_all_yrs <= wc_e)
  m_60_90 = (cru_all_yrs >= cc_s) & (cru_all_yrs <= cc_e)
  t_ref_70_00 = window_mean_and_fill(t_raw, m_70_00)
  p_ref_70_00 = np.maximum(0.0, window_mean_and_fill(p_raw, m_70_00))
  t_ref_60_90 = window_mean_and_fill(t_raw, m_60_90)
  p_ref_60_90 = np.maximum(0.0, window_mean_and_fill(p_raw, m_60_90))

  out_npz.parent.mkdir(parents=True, exist_ok=True)
  np.savez_compressed(
    out_npz,
    years=np.asarray(cru_years, dtype=np.int32),
    lat=lat.astype(np.float32),
    lon=lon.astype(np.float32),
    T_clim=t_clim,
    P_clim=p_clim,
    T_ref_70_00=t_ref_70_00,
    P_ref_70_00=p_ref_70_00,
    T_ref_60_90=t_ref_60_90,
    P_ref_60_90=p_ref_60_90,
  )
  print(f"[CRU] Saved {out_npz} ({n_targets} timesteps, {nlat}x{nlon})")
  return {
    "years": np.asarray(cru_years, dtype=np.int32),
    "lat": lat.astype(np.float32),
    "lon": lon.astype(np.float32),
    "T_clim": t_clim,
    "P_clim": p_clim,
    "T_ref_70_00": t_ref_70_00,
    "P_ref_70_00": p_ref_70_00,
    "T_ref_60_90": t_ref_60_90,
    "P_ref_60_90": p_ref_60_90,
  }


def build_worldclim_reference(
  cfg: dict,
  wc_dir: Path,
  ccsm_data: dict,
  out_npz: Path,
) -> None:
  """Builds gap-free 4320x2160x12 WorldClim 2.1 modern reference over land and emergent paleoshelves."""
  ds_cfg = cfg["downscaling"]
  lapse = float(ds_cfg["lapse_rate_c_per_m"])
  p_floor = float(ds_cfg["precip_floor_mm_month"])
  r_min = float(ds_cfg["precip_ratio_min"])
  r_max = float(ds_cfg["precip_ratio_max"])

  lat_5m, lon_5m = hc.get_grid_coords(4320, 2160)
  h, w = 2160, 4320

  elev_raw = np.array(Image.open(wc_dir / "wc2.1_5m_elev.tif"), dtype=np.float32)
  valid_elev = elev_raw > -30000.0
  z_base = np.where(valid_elev, elev_raw, 0.0).astype(np.float32)

  # Coarse CCSM3 1960-1990 precipitation upsampled to 5' for precipitation ratio inpainting
  p_ccsm_5m = bilinear_upsample_periodic(
    ccsm_data["P_ref"],
    ccsm_data["lat"],
    ccsm_data["lon"],
    lat_5m,
    lon_5m,
  )

  t_base = np.empty((12, h, w), dtype=np.float32)
  p_base = np.empty((12, h, w), dtype=np.float32)

  for m in range(1, 13):
    t_raw = np.array(Image.open(wc_dir / f"wc2.1_5m_tavg_{m:02d}.tif"), dtype=np.float32)
    p_raw = np.array(Image.open(wc_dir / f"wc2.1_5m_prec_{m:02d}.tif"), dtype=np.float32)

    vt = (t_raw > -1e20) & np.isfinite(t_raw) & valid_elev
    t_sl = np.where(vt, t_raw + lapse * z_base, 0.0).astype(np.float32)
    t_sl_filled = inpaint_periodic_2d(t_sl, vt, n_smooth=8)
    t_base[m - 1] = t_sl_filled - lapse * z_base

    vp = (p_raw > -30000.0) & np.isfinite(p_raw)
    denom = np.maximum(p_ccsm_5m[m - 1], p_floor)
    ratio = np.where(vp, np.clip(p_raw / denom, r_min, r_max), 1.0).astype(np.float32)
    ratio_filled = inpaint_periodic_2d(ratio, vp, n_smooth=8)
    p_base[m - 1] = np.maximum(0.0, np.where(vp, p_raw, ratio_filled * p_ccsm_5m[m - 1]))
    print(f"  [WorldClim] month {m:02d}: T=[{t_base[m-1].min():.1f}, {t_base[m-1].max():.1f}] C, P_max={p_base[m-1].max():.0f} mm")

  if np.isnan(t_base).any() or np.isnan(p_base).any() or np.isnan(z_base).any():
    raise RuntimeError("NaN detected in gap-free WorldClim 2.1 reference!")

  out_npz.parent.mkdir(parents=True, exist_ok=True)
  np.savez(
    out_npz,
    T_base_5m=t_base,
    P_base_5m=p_base,
    Z_base_5m=z_base,
  )
  print(f"[WorldClim] Saved gap-free 5' reference {out_npz} (ZERO NaNs)")


def main() -> None:
  paths = hc.resolve_paths("climate")
  cfg = hc.load_json5("climate_config.json5")
  years = hc.get_timestep_years()

  raw_dir = paths["raw_dir"]
  inter_dir = paths["intermediate_dir"]

  print("=== Step B1: Extracting CCSM3 TraCE-21ka climatologies ===")
  ccsm_npz = inter_dir / "ccsm3_climatologies.npz"
  ccsm_data = extract_ccsm3_climatologies(cfg, years, raw_dir / "trace21k", ccsm_npz)

  print("=== Step B2: Extracting CRU TS 4.10 climatologies ===")
  cru_npz = inter_dir / "cru_climatologies.npz"
  extract_cru_climatologies(cfg, years, raw_dir / "cru", cru_npz)

  print("=== Step B3: Building gap-free WorldClim 2.1 5' modern reference ===")
  wc_npz = inter_dir / "worldclim_base_5m.npz"
  build_worldclim_reference(cfg, raw_dir / "worldclim", ccsm_data, wc_npz)

  print("=== Step B complete ===")


if __name__ == "__main__":
  main()
