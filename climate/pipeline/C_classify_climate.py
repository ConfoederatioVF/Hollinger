"""Step C: Downscale monthly T/P with paleo lapse adjustment and classify Köppen-Geiger, Trewartha, and Holdridge Life Zones."""

from __future__ import annotations

import concurrent.futures
import os
import shutil
import sys
import time
from pathlib import Path

import netCDF4 as nc
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "core"))
import hollinger_core as hc
from B_extract_climatologies import bilinear_upsample_periodic

_WORKER_CTX: dict = {}


def load_glac1d_elevation(glac_nc: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
  """Loads GLAC-1D HDC elevation (masl) oriented N->S and wrapped to [-180, 180)."""
  with nc.Dataset(str(glac_nc), "r") as ds:
    t_ka = np.asarray(ds.variables["T122KP11"][:], dtype=np.float64)
    raw_lat = np.asarray(ds.variables["YLATGLOBP5"][:], dtype=np.float64)
    raw_lon = np.asarray(ds.variables["XLONGLOB1"][:], dtype=np.float64)
    hdc = np.asarray(ds.variables["HDC"][:], dtype=np.float32)

  if raw_lat[0] < raw_lat[-1]:
    lat = raw_lat[::-1].copy()
    hdc = hdc[:, ::-1, :]
  else:
    lat = raw_lat.copy()

  lon_w = np.where(raw_lon >= 180.0, raw_lon - 360.0, raw_lon)
  lon_order = np.argsort(lon_w)
  lon = lon_w[lon_order].copy()
  hdc = hdc[:, :, lon_order]
  return t_ka, lat.astype(np.float32), lon.astype(np.float32), hdc


def interp_glac1d_dz(
  year: int,
  t_ka: np.ndarray,
  hdc: np.ndarray,
  hdc_modern: np.ndarray,
) -> np.ndarray:
  """Returns GLAC-1D elevation anomaly (m) relative to 0 ka BP (1950 CE) for a signed year."""
  if year >= 1950:
    return np.zeros_like(hdc_modern, dtype=np.float32)
  ka = max(float(t_ka[0]), min(0.0, (float(year) - 1950.0) / 1000.0))
  idx = int(np.searchsorted(t_ka, ka))
  if idx <= 0:
    h_y = hdc[0]
  elif idx >= len(t_ka):
    h_y = hdc[-1]
  else:
    t0, t1 = t_ka[idx - 1], t_ka[idx]
    w = (ka - t0) / (t1 - t0) if t1 > t0 else 0.0
    h_y = (1.0 - w) * hdc[idx - 1] + w * hdc[idx]
  return (np.maximum(0.0, h_y) - np.maximum(0.0, hdc_modern)).astype(np.float32)


def compute_seasonal_metrics(
  t_m: np.ndarray,
  p_m: np.ndarray,
  nh_mask: np.ndarray,
) -> dict[str, np.ndarray]:
  """Computes common monthly summary metrics for Köppen-Geiger and Trewartha."""
  mat = t_m.mean(axis=0)
  map_ann = p_m.sum(axis=0)
  t_hot = t_m.max(axis=0)
  t_cold = t_m.min(axis=0)
  t_mon10 = (t_m >= 10.0).sum(axis=0)
  p_dry = p_m.min(axis=0)

  p_amjjas = p_m[3:9]
  p_ondjfm = np.concatenate([p_m[9:12], p_m[0:3]], axis=0)

  p_sdry = np.where(nh_mask, p_amjjas.min(axis=0), p_ondjfm.min(axis=0))
  p_swet = np.where(nh_mask, p_amjjas.max(axis=0), p_ondjfm.max(axis=0))
  p_ssum = np.where(nh_mask, p_amjjas.sum(axis=0), p_ondjfm.sum(axis=0))

  p_wdry = np.where(nh_mask, p_ondjfm.min(axis=0), p_amjjas.min(axis=0))
  p_wwet = np.where(nh_mask, p_ondjfm.max(axis=0), p_amjjas.max(axis=0))
  p_wsum = np.where(nh_mask, p_ondjfm.sum(axis=0), p_amjjas.sum(axis=0))

  return {
    "mat": mat,
    "map": map_ann,
    "t_hot": t_hot,
    "t_cold": t_cold,
    "t_mon10": t_mon10,
    "p_dry": p_dry,
    "p_sdry": p_sdry,
    "p_swet": p_swet,
    "p_ssum": p_ssum,
    "p_wdry": p_wdry,
    "p_wwet": p_wwet,
    "p_wsum": p_wsum,
  }


def classify_koeppen_geiger(
  t_m: np.ndarray,
  p_m: np.ndarray,
  m: dict[str, np.ndarray],
  water_mask: np.ndarray,
) -> np.ndarray:
  """Classifies Köppen-Geiger (Peel et al. 2007; Beck et al. 2018, 2023) into uint8 codes 0..30."""
  h, w = water_mask.shape
  out = np.zeros((h, w), dtype=np.uint8)

  mat = m["mat"]
  map_ann = m["map"]
  t_hot = m["t_hot"]
  t_cold = m["t_cold"]
  t_mon10 = (t_m > 10.0).sum(axis=0)
  p_dry = m["p_dry"]
  p_sdry, p_swet, p_ssum = m["p_sdry"], m["p_swet"], m["p_ssum"]
  p_wdry, p_wwet, p_wsum = m["p_wdry"], m["p_wwet"], m["p_wsum"]

  p_thresh = np.where(
    p_wsum >= 0.7 * map_ann,
    2.0 * mat,
    np.where(p_ssum >= 0.7 * map_ann, 2.0 * mat + 28.0, 2.0 * mat + 14.0),
  )

  is_s = (p_sdry < 40.0) & (p_sdry < (p_wwet / 3.0))
  is_w = (~is_s) & (p_wdry < (p_swet / 10.0))

  is_e = t_hot < 10.0
  out = np.where(is_e & (t_hot >= 0.0), 29, out)
  out = np.where(is_e & (t_hot < 0.0), 30, out)

  is_b = (~is_e) & (map_ann < 10.0 * p_thresh)
  is_bw = map_ann < 5.0 * p_thresh
  is_bh = mat >= 18.0
  out = np.where(is_b & is_bw & is_bh, 4, out)
  out = np.where(is_b & is_bw & (~is_bh), 5, out)
  out = np.where(is_b & (~is_bw) & is_bh, 6, out)
  out = np.where(is_b & (~is_bw) & (~is_bh), 7, out)

  rem = (~is_e) & (~is_b)
  is_a = rem & (t_cold >= 18.0)
  is_af = p_dry >= 60.0
  is_am = (~is_af) & (p_dry >= (100.0 - map_ann / 25.0))
  out = np.where(is_a & is_af, 1, out)
  out = np.where(is_a & is_am, 2, out)
  out = np.where(is_a & (~is_af) & (~is_am), 3, out)

  is_c = rem & (t_cold > 0.0) & (t_cold < 18.0)
  hot_a = t_hot >= 22.0
  warm_b = (~hot_a) & (t_mon10 >= 4)
  cold_c = (~hot_a) & (~warm_b)

  out = np.where(is_c & is_s & hot_a, 8, out)
  out = np.where(is_c & is_s & warm_b, 9, out)
  out = np.where(is_c & is_s & cold_c, 10, out)
  out = np.where(is_c & is_w & hot_a, 11, out)
  out = np.where(is_c & is_w & warm_b, 12, out)
  out = np.where(is_c & is_w & cold_c, 13, out)
  is_f = (~is_s) & (~is_w)
  out = np.where(is_c & is_f & hot_a, 14, out)
  out = np.where(is_c & is_f & warm_b, 15, out)
  out = np.where(is_c & is_f & cold_c, 16, out)

  is_d = rem & (t_cold <= 0.0)
  very_cold_d = cold_c & (t_cold < -38.0)
  cold_dc = cold_c & (~very_cold_d)

  out = np.where(is_d & is_s & hot_a, 17, out)
  out = np.where(is_d & is_s & warm_b, 18, out)
  out = np.where(is_d & is_s & cold_dc, 19, out)
  out = np.where(is_d & is_s & very_cold_d, 20, out)
  out = np.where(is_d & is_w & hot_a, 21, out)
  out = np.where(is_d & is_w & warm_b, 22, out)
  out = np.where(is_d & is_w & cold_dc, 23, out)
  out = np.where(is_d & is_w & very_cold_d, 24, out)
  out = np.where(is_d & is_f & hot_a, 25, out)
  out = np.where(is_d & is_f & warm_b, 26, out)
  out = np.where(is_d & is_f & cold_dc, 27, out)
  out = np.where(is_d & is_f & very_cold_d, 28, out)

  out[water_mask] = 0
  return out


def classify_trewartha(
  t_m: np.ndarray,
  p_m: np.ndarray,
  m: dict[str, np.ndarray],
  water_mask: np.ndarray,
) -> np.ndarray:
  """Classifies Trewartha (Trewartha & Horn 1980; Patton 1962; Belda et al. 2014) into uint8 codes 0..17."""
  h, w = water_mask.shape
  out = np.zeros((h, w), dtype=np.uint8)

  mat = m["mat"]
  map_ann = m["map"]
  t_hot = m["t_hot"]
  t_cold = m["t_cold"]
  t_mon10 = m["t_mon10"]
  p_dry = m["p_dry"]
  p_sdry, p_ssum = m["p_sdry"], m["p_ssum"]
  p_wdry, p_wwet, p_swet, p_wsum = m["p_wdry"], m["p_wwet"], m["p_swet"], m["p_wsum"]

  pw_pct = 100.0 * p_wsum / np.maximum(map_ann, 1e-3)
  r_cm = 2.3 * mat - 0.64 * pw_pct + 41.0
  p_thresh_mm = 10.0 * r_cm

  is_f = t_hot < 10.0
  out = np.where(is_f & (t_hot >= 0.0), 16, out)
  out = np.where(is_f & (t_hot < 0.0), 17, out)

  is_b = (~is_f) & (map_ann < p_thresh_mm)
  is_bw = map_ann < (0.5 * p_thresh_mm)
  is_bh = t_mon10 >= 8
  out = np.where(is_b & is_bw & is_bh, 5, out)
  out = np.where(is_b & is_bw & (~is_bh), 6, out)
  out = np.where(is_b & (~is_bw) & is_bh, 7, out)
  out = np.where(is_b & (~is_bw) & (~is_bh), 8, out)

  rem = (~is_f) & (~is_b)

  is_a = rem & (t_cold >= 18.0)
  n_dry_months = (p_m < 60.0).sum(axis=0)
  is_ar = n_dry_months <= 2
  is_am = (~is_ar) & (p_dry >= (100.0 - map_ann / 25.0))
  is_as = (~is_ar) & (~is_am) & (p_ssum < p_wsum)
  out = np.where(is_a & is_ar, 1, out)
  out = np.where(is_a & is_am, 2, out)
  out = np.where(is_a & is_as, 4, out)
  out = np.where(is_a & (~is_ar) & (~is_am) & (~is_as), 3, out)

  is_c = rem & (t_cold < 18.0) & (t_mon10 >= 8)
  is_cs = (p_sdry < 30.0) & (p_sdry < (p_wwet / 3.0)) & (p_ssum < p_wsum) & (map_ann < 890.0)
  is_cw = (~is_cs) & (p_wdry < (p_swet / 10.0)) & (p_wsum < p_ssum)
  out = np.where(is_c & is_cs, 9, out)
  out = np.where(is_c & is_cw, 10, out)
  out = np.where(is_c & (~is_cs) & (~is_cw), 11, out)

  is_d = rem & (t_cold < 18.0) & (t_mon10 >= 4) & (t_mon10 <= 7)
  out = np.where(is_d & (t_cold > 0.0), 12, out)
  out = np.where(is_d & (t_cold <= 0.0), 13, out)

  is_e = rem & (t_cold < 18.0) & (t_mon10 >= 1) & (t_mon10 <= 3)
  out = np.where(is_e & (t_cold > 0.0), 14, out)
  out = np.where(is_e & (t_cold <= 0.0), 15, out)

  out[water_mask] = 0
  return out


def build_holdridge_centroids(cfg: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
  """Builds canonical Holdridge centroids (excluding Polar_Desert and Subtrop duplicate rows)."""
  hcfg = cfg["holdridge_lifezones"]
  ids = []
  log_t = []
  log_p = []
  for cls in hcfg["classes"]:
    cid = int(cls["id"])
    if cid <= 1 or (24 <= cid <= 30):
      continue
    ids.append(cid)
    log_t.append(np.log2(float(cls["t_c"])))
    log_p.append(np.log2(float(cls["p_c"])))
  return (
    np.asarray(ids, dtype=np.uint8),
    np.asarray(log_t, dtype=np.float32),
    np.asarray(log_p, dtype=np.float32),
  )


def classify_holdridge(
  t_m: np.ndarray,
  p_m: np.ndarray,
  t_cold: np.ndarray,
  water_mask: np.ndarray,
  c_ids: np.ndarray,
  c_log_t: np.ndarray,
  c_log_p: np.ndarray,
  polar_tbio_c: float,
  frost_tcold_c: float,
) -> np.ndarray:
  """Classifies Holdridge Life Zones (0..38) via log2 hexagonal centroid nearest-neighbour + frost line."""
  h, w = water_mask.shape
  out = np.zeros((h, w), dtype=np.uint8)

  t_bio = np.clip(t_m, 0.0, 30.0).mean(axis=0)
  p_ann = np.maximum(p_m.sum(axis=0), 1.0)

  polar_mask = (~water_mask) & (t_bio < polar_tbio_c)
  out[polar_mask] = 1

  active_mask = (~water_mask) & (~polar_mask)
  lt = np.log2(np.maximum(t_bio[active_mask], 1e-3))
  lp = np.log2(p_ann[active_mask])

  n_act = lt.size
  best_ids = np.empty(n_act, dtype=np.uint8)
  chunk = 500_000
  for i0 in range(0, n_act, chunk):
    i1 = min(n_act, i0 + chunk)
    dt = lt[i0:i1, None] - c_log_t[None, :]
    dp = lp[i0:i1, None] - c_log_p[None, :]
    d2 = dt * dt + dp * dp
    best_ids[i0:i1] = c_ids[np.argmin(d2, axis=1)]

  tc_act = t_cold[active_mask]
  is_warm_subtrop = (best_ids >= 17) & (best_ids <= 23)
  best_ids = np.where(is_warm_subtrop & (tc_act >= frost_tcold_c), best_ids + 7, best_ids)

  out[active_mask] = best_ids
  out[water_mask] = 0
  return out


def _init_worker(shards_dir_str: str) -> None:
  """Initialises per-process read-only climate arrays."""
  paths = hc.resolve_paths("climate")
  cfg = hc.load_json5("climate_config.json5")
  ds_cfg = cfg["downscaling"]
  inter_dir = paths["intermediate_dir"]

  lat_5m, lon_5m = hc.get_grid_coords(4320, 2160)
  nh_mask = (lat_5m >= 0.0)[:, None]
  ant_mask = (lat_5m < float(ds_cfg["antarctica_lat_cutoff"]))[:, None]

  wc = np.load(inter_dir / "worldclim_base_5m.npz")
  ccsm = np.load(inter_dir / "ccsm3_climatologies.npz")
  cru = np.load(inter_dir / "cru_climatologies.npz")

  p_floor = float(ds_cfg["precip_floor_mm_month"])
  r_min = float(ds_cfg["precip_ratio_min"])
  r_max = float(ds_cfg["precip_ratio_max"])

  dt_splice_5m = bilinear_upsample_periodic(
    cru["T_ref_60_90"] - cru["T_ref_70_00"], cru["lat"], cru["lon"], lat_5m, lon_5m
  )
  rp_splice_5m = bilinear_upsample_periodic(
    np.clip(cru["P_ref_60_90"] / np.maximum(cru["P_ref_70_00"], p_floor), r_min, r_max),
    cru["lat"],
    cru["lon"],
    lat_5m,
    lon_5m,
  )
  dt_splice_5m = np.where(ant_mask[None, :, :], 0.0, dt_splice_5m).astype(np.float32)
  rp_splice_5m = np.where(ant_mask[None, :, :], 1.0, rp_splice_5m).astype(np.float32)

  glac_nc = paths["base_dir"] / "landarea" / "1.raw_data" / "TOPicemsk.GLACD26kN9894GE90227A6005GGrBgic.nc"
  t_ka, glac_lat, glac_lon, glac_hdc = load_glac1d_elevation(glac_nc)

  c_ids, c_log_t, c_log_p = build_holdridge_centroids(cfg)
  hcfg = cfg["holdridge_lifezones"]

  _WORKER_CTX.update({
    "cfg": cfg,
    "paths": paths,
    "shards_dir": Path(shards_dir_str),
    "lat_5m": lat_5m,
    "lon_5m": lon_5m,
    "nh_mask": nh_mask,
    "ant_mask": ant_mask,
    "lapse": float(ds_cfg["lapse_rate_c_per_m"]),
    "p_floor": p_floor,
    "r_min": r_min,
    "r_max": r_max,
    "splice_start": int(ds_cfg["cru_splice_start"]),
    "splice_end": int(ds_cfg["cru_splice_end"]),
    "t_base_5m": wc["T_base_5m"],
    "p_base_5m": wc["P_base_5m"],
    "ccsm_idx": {int(y): i for i, y in enumerate(ccsm["years"].tolist())},
    "ccsm_lat": ccsm["lat"],
    "ccsm_lon": ccsm["lon"],
    "ccsm_t": ccsm["T_clim"],
    "ccsm_p": ccsm["P_clim"],
    "ccsm_z": ccsm["Z_clim"],
    "ccsm_t_ref": ccsm["T_ref"],
    "ccsm_p_ref": ccsm["P_ref"],
    "ccsm_z_ref": ccsm["Z_ref"],
    "cru_idx": {int(y): i for i, y in enumerate(cru["years"].tolist())},
    "cru_lat": cru["lat"],
    "cru_lon": cru["lon"],
    "cru_t": cru["T_clim"],
    "cru_p": cru["P_clim"],
    "cru_t_70_00": cru["T_ref_70_00"],
    "cru_p_70_00": cru["P_ref_70_00"],
    "dt_splice_5m": dt_splice_5m,
    "rp_splice_5m": rp_splice_5m,
    "t_ka": t_ka,
    "glac_lat": glac_lat,
    "glac_lon": glac_lon,
    "glac_hdc": glac_hdc,
    "glac_hdc_modern": glac_hdc[-1],
    "c_ids": c_ids,
    "c_log_t": c_log_t,
    "c_log_p": c_log_p,
    "polar_tbio_c": float(hcfg["polar_tbio_threshold_c"]),
    "frost_tcold_c": float(hcfg["subtropical_frost_tcold_c"]),
    "land_nc": nc.Dataset(str(paths["base_dir"] / "landarea" / "landarea.nc"), "r"),
  })


def _process_timestep(args: tuple[int, int]) -> tuple[int, int, float]:
  """Processes a single timestep: downscales, classifies, writes GeoPNGs + shard."""
  ti, y = args
  t0 = time.time()
  ctx = _WORKER_CTX
  lat_5m, lon_5m = ctx["lat_5m"], ctx["lon_5m"]
  ant_mask = ctx["ant_mask"]
  p_floor, r_min, r_max = ctx["p_floor"], ctx["r_min"], ctx["r_max"]

  water_mask = np.asarray(ctx["land_nc"].variables["land_area"][ti, :, :], dtype=np.float32) == 0.0

  if y < ctx["splice_start"]:
    ci = ctx["ccsm_idx"][y]
    dt_ccsm_5m = bilinear_upsample_periodic(
      ctx["ccsm_t"][ci] - ctx["ccsm_t_ref"], ctx["ccsm_lat"], ctx["ccsm_lon"], lat_5m, lon_5m
    )
    rp_ccsm_5m = bilinear_upsample_periodic(
      np.clip(ctx["ccsm_p"][ci] / np.maximum(ctx["ccsm_p_ref"], p_floor), r_min, r_max),
      ctx["ccsm_lat"],
      ctx["ccsm_lon"],
      lat_5m,
      lon_5m,
    )
    dt_5m = dt_ccsm_5m + ctx["dt_splice_5m"]
    rp_5m = np.clip(rp_ccsm_5m * ctx["rp_splice_5m"], r_min, r_max)
  elif y <= ctx["splice_end"]:
    w = float(y - ctx["splice_start"]) / float(ctx["splice_end"] - ctx["splice_start"])
    ci = ctx["ccsm_idx"][y]
    ri = ctx["cru_idx"][y]

    dt_ccsm_5m = bilinear_upsample_periodic(
      ctx["ccsm_t"][ci] - ctx["ccsm_t_ref"], ctx["ccsm_lat"], ctx["ccsm_lon"], lat_5m, lon_5m
    )
    rp_ccsm_5m = bilinear_upsample_periodic(
      np.clip(ctx["ccsm_p"][ci] / np.maximum(ctx["ccsm_p_ref"], p_floor), r_min, r_max),
      ctx["ccsm_lat"],
      ctx["ccsm_lon"],
      lat_5m,
      lon_5m,
    )
    dt_ccsm_anc = dt_ccsm_5m + ctx["dt_splice_5m"]
    rp_ccsm_anc = np.clip(rp_ccsm_5m * ctx["rp_splice_5m"], r_min, r_max)

    dt_cru_5m = bilinear_upsample_periodic(
      ctx["cru_t"][ri] - ctx["cru_t_70_00"], ctx["cru_lat"], ctx["cru_lon"], lat_5m, lon_5m
    )
    rp_cru_5m = bilinear_upsample_periodic(
      np.clip(ctx["cru_p"][ri] / np.maximum(ctx["cru_p_70_00"], p_floor), r_min, r_max),
      ctx["cru_lat"],
      ctx["cru_lon"],
      lat_5m,
      lon_5m,
    )
    dt_5m = np.where(ant_mask[None, :, :], dt_ccsm_anc, (1.0 - w) * dt_ccsm_anc + w * dt_cru_5m)
    rp_5m = np.where(ant_mask[None, :, :], rp_ccsm_anc, (1.0 - w) * rp_ccsm_anc + w * rp_cru_5m)
  else:
    ri = ctx["cru_idx"][y]
    ci = ctx["ccsm_idx"][min(y, 1990)]
    dt_cru_5m = bilinear_upsample_periodic(
      ctx["cru_t"][ri] - ctx["cru_t_70_00"], ctx["cru_lat"], ctx["cru_lon"], lat_5m, lon_5m
    )
    rp_cru_5m = bilinear_upsample_periodic(
      np.clip(ctx["cru_p"][ri] / np.maximum(ctx["cru_p_70_00"], p_floor), r_min, r_max),
      ctx["cru_lat"],
      ctx["cru_lon"],
      lat_5m,
      lon_5m,
    )
    dt_ccsm_5m = bilinear_upsample_periodic(
      ctx["ccsm_t"][ci] - ctx["ccsm_t_ref"], ctx["ccsm_lat"], ctx["ccsm_lon"], lat_5m, lon_5m
    )
    rp_ccsm_5m = bilinear_upsample_periodic(
      np.clip(ctx["ccsm_p"][ci] / np.maximum(ctx["ccsm_p_ref"], p_floor), r_min, r_max),
      ctx["ccsm_lat"],
      ctx["ccsm_lon"],
      lat_5m,
      lon_5m,
    )
    dt_5m = np.where(ant_mask[None, :, :], dt_ccsm_5m, dt_cru_5m)
    rp_5m = np.where(ant_mask[None, :, :], rp_ccsm_5m, rp_cru_5m)

  if y < 1950:
    ci_z = ctx["ccsm_idx"][y]
    dz_ccsm_5m = bilinear_upsample_periodic(
      ctx["ccsm_z"][ci_z] - ctx["ccsm_z_ref"], ctx["ccsm_lat"], ctx["ccsm_lon"], lat_5m, lon_5m
    )
    dz_glac_5m = bilinear_upsample_periodic(
      interp_glac1d_dz(y, ctx["t_ka"], ctx["glac_hdc"], ctx["glac_hdc_modern"]),
      ctx["glac_lat"],
      ctx["glac_lon"],
      lat_5m,
      lon_5m,
    )
    dz_res = dz_glac_5m - dz_ccsm_5m
    t_m = ctx["t_base_5m"] + dt_5m - ctx["lapse"] * dz_res[None, :, :]
  else:
    t_m = ctx["t_base_5m"] + dt_5m

  p_m = np.maximum(0.0, ctx["p_base_5m"] * rp_5m)

  land_idx = ~water_mask
  n_land = int(land_idx.sum())
  t_land = t_m[:, land_idx][:, None, :]
  p_land = p_m[:, land_idx][:, None, :]
  nh_land = np.broadcast_to(ctx["nh_mask"], water_mask.shape)[land_idx][None, :]
  w_false = np.zeros((1, n_land), dtype=bool)

  metrics = compute_seasonal_metrics(t_land, p_land, nh_land)
  kg_1d = classify_koeppen_geiger(t_land, p_land, metrics, w_false)[0]
  tr_1d = classify_trewartha(t_land, p_land, metrics, w_false)[0]
  hl_1d = classify_holdridge(
    t_land,
    p_land,
    metrics["t_cold"],
    w_false,
    ctx["c_ids"],
    ctx["c_log_t"],
    ctx["c_log_p"],
    ctx["polar_tbio_c"],
    ctx["frost_tcold_c"],
  )[0]

  kg = np.zeros(water_mask.shape, dtype=np.uint8)
  tr = np.zeros(water_mask.shape, dtype=np.uint8)
  hl = np.zeros(water_mask.shape, dtype=np.uint8)
  kg[land_idx] = kg_1d
  tr[land_idx] = tr_1d
  hl[land_idx] = hl_1d

  rasters_dir = ctx["paths"]["rasters_dir"]
  cfg = ctx["cfg"]
  for key, arr in (("koeppen_geiger", kg), ("trewartha", tr), ("holdridge_lifezones", hl)):
    png_path = rasters_dir / f"{cfg[key]['geopng_prefix']}{y}.png"
    hc.encode_uint32_geopng(png_path, arr.astype(np.uint32))

  shard_path = ctx["shards_dir"] / f"t_{ti:03d}.npy"
  np.save(shard_path, np.stack([kg, tr, hl], axis=0))
  return ti, y, time.time() - t0


def open_climate_outputs(
  cfg: dict,
  paths: dict,
  years: list[int],
  lat: np.ndarray,
  lon: np.ndarray,
  cell_area: np.ndarray,
) -> dict:
  """Opens all 3 CF-1.8 NetCDF-4 datasets for slice-by-slice writing."""
  enc = cfg["encoding"]
  chunk_shape = tuple(int(x) for x in enc["netcdf_chunk_shape"])
  complevel = int(enc["netcdf_compression_level"])

  out_handles = {}
  specs = [
    (
      "koeppen_geiger",
      "Köppen-Geiger climate classification (Peel et al. 2007; Beck et al. 2018, 2023)",
      "CCSM3 TraCE-21ka (d651050) + WorldClim 2.1 (5') + CRU TS 4.10 + GLAC-1D + landarea.nc",
    ),
    (
      "trewartha",
      "Trewartha climate classification (Trewartha & Horn 1980; Patton 1962; Belda et al. 2014)",
      "CCSM3 TraCE-21ka (d651050) + WorldClim 2.1 (5') + CRU TS 4.10 + GLAC-1D + landarea.nc",
    ),
    (
      "holdridge_lifezones",
      "Holdridge Life Zones classification (Holdridge 1947, 1967; Leemans 1990; Lugo et al. 1999)",
      "CCSM3 TraCE-21ka (d651050) + WorldClim 2.1 (5') + CRU TS 4.10 + GLAC-1D + landarea.nc",
    ),
  ]

  for key, long_name, source_str in specs:
    sub_cfg = cfg[key]
    classes = sub_cfg["classes"]
    flag_values = [int(c["id"]) for c in classes]
    flag_meanings = " ".join(c["code"] for c in classes)
    nc_path = paths["indicator_dir"] / sub_cfg["nc_filename"]
    ds, var = hc.create_categorical_netcdf4(
      nc_path=nc_path,
      var_name=sub_cfg["var_name"],
      years=years,
      lat=lat,
      lon=lon,
      cell_area_km2=cell_area,
      flag_values=flag_values,
      flag_meanings=flag_meanings,
      global_attrs={
        "title": long_name,
        "institution": "SVEA / landuse_Hollinger",
        "source": source_str,
        "history": "Generated by climate/pipeline/C_classify_climate.py",
      },
      var_long_name=long_name,
      chunk_shape=chunk_shape,
      complevel=complevel,
    )
    out_handles[key] = {"ds": ds, "var": var}
  return out_handles


def main() -> None:
  t_start = time.time()
  paths = hc.resolve_paths("climate")
  cfg = hc.load_json5("climate_config.json5")
  mcfg = hc.load_json5("model_config.json5")
  n_workers = min(8, int(mcfg["compute"]["max_workers"]), os.cpu_count() or 4)
  years = hc.get_timestep_years()

  lat_5m, lon_5m = hc.get_grid_coords(4320, 2160)
  cell_area = hc.compute_cell_area_by_row(4320, 2160)

  rasters_dir = paths["rasters_dir"]
  rasters_dir.mkdir(parents=True, exist_ok=True)
  shards_dir = paths["intermediate_dir"] / "shards"
  shards_dir.mkdir(parents=True, exist_ok=True)

  with nc.Dataset(str(paths["base_dir"] / "landarea" / "landarea.nc"), "r") as ds_land:
    land_years = list(np.asarray(ds_land.variables["time"][:], dtype=np.int32).tolist())
    if land_years != years:
      raise ValueError("landarea.nc time axis does not match canonical 218 timesteps")

  print(f"=== Step C: Classifying {len(years)} timesteps with {n_workers} worker processes ===", flush=True)
  done = 0
  with concurrent.futures.ProcessPoolExecutor(
    max_workers=n_workers,
    initializer=_init_worker,
    initargs=(str(shards_dir),),
  ) as pool:
    for ti, y, dt in pool.map(_process_timestep, list(enumerate(years))):
      done += 1
      if done % 20 == 0 or done == 1 or done == len(years):
        print(f"  [{done:3d}/{len(years)}] year={y:6d} ({dt:.2f}s/step)", flush=True)

  print("=== Step C: Assembling 3 compressed CF-1.8 NetCDF-4 datasets ===", flush=True)
  handles = open_climate_outputs(cfg, paths, years, lat_5m, lon_5m, cell_area)
  try:
    for ti in range(len(years)):
      stack = np.load(shards_dir / f"t_{ti:03d}.npy")
      handles["koeppen_geiger"]["var"][ti, :, :] = stack[0]
      handles["trewartha"]["var"][ti, :, :] = stack[1]
      handles["holdridge_lifezones"]["var"][ti, :, :] = stack[2]
  finally:
    for h_obj in handles.values():
      h_obj["ds"].close()

  shutil.rmtree(shards_dir, ignore_errors=True)
  print(f"=== Step C complete in {time.time() - t_start:.1f}s ===", flush=True)


if __name__ == "__main__":
  main()
