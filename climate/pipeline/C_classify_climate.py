"""Step C: Downscale monthly T/P with paleo lapse adjustment and classify Köppen-Geiger, Trewartha, and Holdridge Life Zones."""

from __future__ import annotations

import argparse
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

  raw_s = (p_sdry < 40.0) & (p_sdry < (p_wwet / 3.0))
  raw_w = p_wdry < (p_swet / 10.0)
  # Peel et al. (2007, p. 1636): when both dry-summer (s) and dry-winter (w) criteria
  # hold (e.g. monsoonal northern India, Mexico, Sub-Saharan Africa), assign w if total
  # summer half-year precipitation exceeds winter, and s only when winter >= summer.
  is_s = raw_s & ((~raw_w) | (p_wsum >= p_ssum))
  is_w = raw_w & (~is_s)

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


def correct_trace_hudson_bay_artifact(
  years: np.ndarray,
  lat: np.ndarray,
  lon: np.ndarray,
  t_clim: np.ndarray,
  p_clim: np.ndarray,
  z_clim: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
  """Smooths the unphysical TraCE-21ka 8.0-7.6 ka BP (-6000..-5700) dry-soil Hudson Bay heatwave."""
  yr_list = [int(y) for y in years.tolist()]
  if -6100 not in yr_list or -5600 not in yr_list:
    return t_clim, p_clim, z_clim

  i_pre = yr_list.index(-6100)
  i_post = yr_list.index(-5600)
  t_out = np.array(t_clim, dtype=np.float32, copy=True)
  p_out = np.array(p_clim, dtype=np.float32, copy=True)
  z_out = np.array(z_clim, dtype=np.float32, copy=True)

  lat_arr = np.asarray(lat, dtype=np.float32)[:, None]
  lon_arr = np.asarray(lon, dtype=np.float32)[None, :]

  w_lat = np.ones_like(lat_arr)
  w_lat = np.where(lat_arr < 42.0, 0.0, w_lat)
  w_lat = np.where((lat_arr >= 42.0) & (lat_arr < 50.0), 0.5 * (1.0 - np.cos(np.pi * (lat_arr - 42.0) / 8.0)), w_lat)
  w_lat = np.where((lat_arr > 74.0) & (lat_arr <= 82.0), 0.5 * (1.0 + np.cos(np.pi * (lat_arr - 74.0) / 8.0)), w_lat)
  w_lat = np.where(lat_arr > 82.0, 0.0, w_lat)

  w_lon = np.ones_like(lon_arr)
  w_lon = np.where(lon_arr < -115.0, 0.0, w_lon)
  w_lon = np.where((lon_arr >= -115.0) & (lon_arr < -105.0), 0.5 * (1.0 - np.cos(np.pi * (lon_arr + 115.0) / 10.0)), w_lon)
  w_lon = np.where((lon_arr > -62.0) & (lon_arr <= -52.0), 0.5 * (1.0 + np.cos(np.pi * (lon_arr + 62.0) / 10.0)), w_lon)
  w_lon = np.where(lon_arr > -52.0, 0.0, w_lon)

  w_spat = (w_lat * w_lon).astype(np.float32)

  for y in (-6000, -5900, -5800, -5700):
    if y not in yr_list:
      continue
    iy = yr_list.index(y)
    alpha = float(y - (-6100)) / float((-5600) - (-6100))
    t_interp = (1.0 - alpha) * t_clim[i_pre] + alpha * t_clim[i_post]
    p_interp = (1.0 - alpha) * p_clim[i_pre] + alpha * p_clim[i_post]
    z_interp = (1.0 - alpha) * z_clim[i_pre] + alpha * z_clim[i_post]
    t_out[iy] = (1.0 - w_spat[None, :, :]) * t_clim[iy] + w_spat[None, :, :] * t_interp
    p_out[iy] = (1.0 - w_spat[None, :, :]) * p_clim[iy] + w_spat[None, :, :] * p_interp
    z_out[iy] = (1.0 - w_spat) * z_clim[iy] + w_spat * z_interp

  return t_out, p_out, z_out


def _smoothstep_cos(x: np.ndarray | float, e0: np.ndarray | float, e1: np.ndarray | float) -> np.ndarray:
  """Smooth cosine ramp clipped cleanly to [0.0, 1.0]."""
  t = np.clip((x - e0) / (e1 - e0), 0.0, 1.0)
  return (0.5 * (1.0 - np.cos(np.pi * t))).astype(np.float32)


def apply_beyer2020_africa_integration(
  y: int,
  t_m: np.ndarray,
  p_m: np.ndarray,
  ctx: dict,
) -> tuple[np.ndarray, np.ndarray]:
  """Blends WorldClim-downscaled Beyer et al. (2020) Holocene climatologies over Africa-Arabia for y < end_year_cutoff."""
  bcfg = ctx["cfg"].get("beyer2020_africa_integration", {})
  bey = ctx.get("beyer")
  if not bey or not bcfg.get("enabled", False) or y >= int(bcfg["end_year_cutoff"]):
    return t_m, p_m

  fade_start = float(bcfg["fade_start_year"])
  end_cut = float(bcfg["end_year_cutoff"])
  w_time = 1.0 if float(y) <= fade_start else 1.0 - float(_smoothstep_cos(float(y), fade_start, end_cut))
  if w_time <= 0.0:
    return t_m, p_m

  by_yrs = bey["years_cal"]
  if y <= int(by_yrs[0]):
    t_b = bey["T_clim"][0]
    p_b = bey["P_clim"][0]
  elif y >= int(by_yrs[-1]):
    t_b = bey["T_clim"][-1]
    p_b = bey["P_clim"][-1]
  else:
    idx = int(np.searchsorted(by_yrs, y, side="right")) - 1
    alpha = float(y - int(by_yrs[idx])) / float(int(by_yrs[idx + 1]) - int(by_yrs[idx]))
    t_b = (1.0 - alpha) * bey["T_clim"][idx] + alpha * bey["T_clim"][idx + 1]
    p_b = (1.0 - alpha) * bey["P_clim"][idx] + alpha * bey["P_clim"][idx + 1]

  sub_lat = bey["sub_lat"]
  sub_lon = bey["sub_lon"]
  b_lat = bey["lat"]
  b_lon = bey["lon"]
  r0, r1, c0, c1 = bey["r0"], bey["r1"], bey["c0"], bey["c1"]

  dt_beyer_5m = bilinear_upsample_periodic(t_b - bey["T_ref"], b_lat, b_lon, sub_lat, sub_lon)
  p_beyer_raw = bilinear_upsample_periodic(p_b, b_lat, b_lon, sub_lat, sub_lon)

  t_base_sub = ctx["t_base_5m"][:, r0:r1, c0:c1]
  p_base_sub = ctx["p_base_5m"][:, r0:r1, c0:c1]
  p_ref_5m = bey["p_ref_5m"]

  t_beyer_down = t_base_sub + dt_beyer_5m
  p_dry = p_base_sub * (p_beyer_raw / np.maximum(p_ref_5m, 0.1))
  p_add = p_base_sub + (p_beyer_raw - p_ref_5m)
  p_rat = p_base_sub * np.minimum(p_beyer_raw / np.maximum(p_ref_5m, 1.0), float(bcfg["precip_ratio_max"]))
  w_wet = bey["w_wet_5m"]
  p_wet = (1.0 - w_wet) * p_add + w_wet * p_rat
  p_beyer_down = np.maximum(0.0, np.where(p_beyer_raw <= p_ref_5m, p_dry, p_wet))

  w_eff = (np.float32(w_time) * bey["w_spat"])[None, :, :]
  t_m[:, r0:r1, c0:c1] = (1.0 - w_eff) * t_m[:, r0:r1, c0:c1] + w_eff * t_beyer_down
  p_m[:, r0:r1, c0:c1] = (1.0 - w_eff) * p_m[:, r0:r1, c0:c1] + w_eff * p_beyer_down
  return t_m, p_m


def _init_worker(shards_dir_str: str, active_keys: tuple[str, ...]) -> None:
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
  ccsm_t_corr, ccsm_p_corr, ccsm_z_corr = correct_trace_hudson_bay_artifact(
    ccsm["years"], ccsm["lat"], ccsm["lon"], ccsm["T_clim"], ccsm["P_clim"], ccsm["Z_clim"]
  )

  beyer_data = None
  bcfg = cfg.get("beyer2020_africa_integration", {})
  beyer_path = inter_dir / "beyer2020_holocene_climatologies.npz"
  if bcfg.get("enabled", False) and beyer_path.exists():
    lat_min = float(bcfg["domain_lat_min"])
    lat_max = float(bcfg["domain_lat_max"])
    lon_min = float(bcfg["domain_lon_min"])
    lon_max = float(bcfg["domain_lon_max"])
    apron_lat = float(bcfg["apron_lat_deg"])
    apron_lon = float(bcfg["apron_lon_deg"])

    r0 = int(np.where(lat_5m <= lat_max)[0][0])
    r1 = int(np.where(lat_5m >= lat_min)[0][-1]) + 1
    c0 = int(np.where(lon_5m >= lon_min)[0][0])
    c1 = int(np.where(lon_5m <= lon_max)[0][-1]) + 1
    sub_lat = lat_5m[r0:r1].astype(np.float32)
    sub_lon = lon_5m[c0:c1].astype(np.float32)

    lat_g = sub_lat[:, None]
    lon_g = sub_lon[None, :]
    w_lat = _smoothstep_cos(lat_g, lat_min, lat_min + apron_lat) * (
      1.0 - _smoothstep_cos(lat_g, lat_max - apron_lat, lat_max)
    )
    w_lon = _smoothstep_cos(lon_g, lon_min, lon_min + apron_lon) * (
      1.0 - _smoothstep_cos(lon_g, lon_max - apron_lon, lon_max)
    )
    w_spat = (w_lat * w_lon).astype(np.float32)

    with np.load(beyer_path) as bf:
      b_lat = np.asarray(bf["lat"], dtype=np.float32)
      b_lon = np.asarray(bf["lon"], dtype=np.float32)
      p_ref_arr = np.asarray(bf["P_ref"], dtype=np.float32)
      p_ref_5m = bilinear_upsample_periodic(p_ref_arr, b_lat, b_lon, sub_lat, sub_lon)
      w_wet_5m = _smoothstep_cos(
        p_ref_5m,
        float(bcfg["arid_blend_lo_mm_month"]),
        float(bcfg["arid_blend_hi_mm_month"]),
      )
      beyer_data = {
        "r0": r0,
        "r1": r1,
        "c0": c0,
        "c1": c1,
        "sub_lat": sub_lat,
        "sub_lon": sub_lon,
        "w_spat": w_spat,
        "years_cal": np.asarray(bf["years_cal"], dtype=np.int32),
        "lat": b_lat,
        "lon": b_lon,
        "T_clim": np.asarray(bf["T_clim"], dtype=np.float32),
        "P_clim": np.asarray(bf["P_clim"], dtype=np.float32),
        "T_ref": np.asarray(bf["T_ref"], dtype=np.float32),
        "p_ref_5m": p_ref_5m,
        "w_wet_5m": w_wet_5m,
      }

  _WORKER_CTX.update({
    "beyer": beyer_data,
    "cfg": cfg,
    "paths": paths,
    "active_keys": tuple(active_keys),
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
    "ccsm_t": ccsm_t_corr,
    "ccsm_p": ccsm_p_corr,
    "ccsm_z": ccsm_z_corr,
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

  land_arr = np.asarray(ctx["land_nc"].variables["land_area"][ti, :, :], dtype=np.float32)
  water_mask = (land_arr == 0.0)
  ice_mask = (land_arr < 0.0)

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
  if y < int(ctx["cfg"].get("beyer2020_africa_integration", {}).get("end_year_cutoff", -100000)):
    t_m, p_m = apply_beyer2020_africa_integration(y, t_m, p_m, ctx)

  land_idx = ~water_mask
  n_land = int(land_idx.sum())
  t_land = t_m[:, land_idx][:, None, :]
  p_land = p_m[:, land_idx][:, None, :]
  nh_land = np.broadcast_to(ctx["nh_mask"], water_mask.shape)[land_idx][None, :]
  w_false = np.zeros((1, n_land), dtype=bool)

  metrics = compute_seasonal_metrics(t_land, p_land, nh_land)
  active_keys = ctx["active_keys"]
  rasters_dir = ctx["paths"]["rasters_dir"]
  cfg = ctx["cfg"]
  shard_arrays = []

  if "koeppen_geiger" in active_keys:
    kg = np.zeros(water_mask.shape, dtype=np.uint8)
    kg[land_idx] = classify_koeppen_geiger(t_land, p_land, metrics, w_false)[0]
    kg[ice_mask] = 30
    hc.encode_uint32_geopng(rasters_dir / f"{cfg['koeppen_geiger']['geopng_prefix']}{y}.png", kg.astype(np.uint32))
    shard_arrays.append(kg)

  if "trewartha" in active_keys:
    tr = np.zeros(water_mask.shape, dtype=np.uint8)
    tr[land_idx] = classify_trewartha(t_land, p_land, metrics, w_false)[0]
    tr[ice_mask] = 17
    hc.encode_uint32_geopng(rasters_dir / f"{cfg['trewartha']['geopng_prefix']}{y}.png", tr.astype(np.uint32))
    shard_arrays.append(tr)

  if "holdridge_lifezones" in active_keys:
    hl = np.zeros(water_mask.shape, dtype=np.uint8)
    hl[land_idx] = classify_holdridge(
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
    hl[ice_mask] = 1
    hc.encode_uint32_geopng(rasters_dir / f"{cfg['holdridge_lifezones']['geopng_prefix']}{y}.png", hl.astype(np.uint32))
    shard_arrays.append(hl)

  shard_path = ctx["shards_dir"] / f"t_{ti:03d}.npy"
  np.save(shard_path, np.stack(shard_arrays, axis=0))
  return ti, y, time.time() - t0


def open_climate_outputs(
  cfg: dict,
  paths: dict,
  years: list[int],
  lat: np.ndarray,
  lon: np.ndarray,
  cell_area: np.ndarray,
  active_keys: tuple[str, ...],
) -> dict:
  """Opens requested CF-1.8 NetCDF-4 datasets for slice-by-slice writing."""
  enc = cfg["encoding"]
  chunk_shape = tuple(int(x) for x in enc["netcdf_chunk_shape"])
  complevel = int(enc["netcdf_compression_level"])

  out_handles = {}
  specs = [
    (
      "koeppen_geiger",
      "Köppen-Geiger climate classification (Peel et al. 2007; Beck et al. 2018, 2023)",
      "CCSM3 TraCE-21ka (d651050) + Beyer et al. (2020) Africa Holocene + WorldClim 2.1 (5') + CRU TS 4.10 + GLAC-1D + landarea.nc",
    ),
    (
      "trewartha",
      "Trewartha climate classification (Trewartha & Horn 1980; Patton 1962; Belda et al. 2014)",
      "CCSM3 TraCE-21ka (d651050) + Beyer et al. (2020) Africa Holocene + WorldClim 2.1 (5') + CRU TS 4.10 + GLAC-1D + landarea.nc",
    ),
    (
      "holdridge_lifezones",
      "Holdridge Life Zones classification (Holdridge 1947, 1967; Leemans 1990; Lugo et al. 1999)",
      "CCSM3 TraCE-21ka (d651050) + Beyer et al. (2020) Africa Holocene + WorldClim 2.1 (5') + CRU TS 4.10 + GLAC-1D + landarea.nc",
    ),
  ]

  for key, long_name, source_str in specs:
    if key not in active_keys:
      continue
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
  all_keys = ("koeppen_geiger", "trewartha", "holdridge_lifezones")
  parser = argparse.ArgumentParser(description="Classify Hollinger climate indicators.")
  parser.add_argument("--only", type=str, default="", help="Comma-separated target keys (default: all 3)")
  args = parser.parse_args()
  if args.only.strip():
    req = [k.strip() for k in args.only.split(",") if k.strip()]
    unknown = [k for k in req if k not in all_keys]
    if unknown:
      raise ValueError(f"Unknown --only target(s) {unknown}; valid: {all_keys}")
    active_keys = tuple(k for k in all_keys if k in req)
  else:
    active_keys = all_keys

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

  print(f"=== Step C: Classifying {active_keys} ({len(years)} timesteps) with {n_workers} worker processes ===", flush=True)
  done = 0
  with concurrent.futures.ProcessPoolExecutor(
    max_workers=n_workers,
    initializer=_init_worker,
    initargs=(str(shards_dir), active_keys),
  ) as pool:
    for ti, y, dt in pool.map(_process_timestep, list(enumerate(years))):
      done += 1
      if done % 20 == 0 or done == 1 or done == len(years):
        print(f"  [{done:3d}/{len(years)}] year={y:6d} ({dt:.2f}s/step)", flush=True)

  print(f"=== Step C: Assembling {len(active_keys)} compressed CF-1.8 NetCDF-4 dataset(s) ===", flush=True)
  handles = open_climate_outputs(cfg, paths, years, lat_5m, lon_5m, cell_area, active_keys)
  try:
    for ti in range(len(years)):
      stack = np.load(shards_dir / f"t_{ti:03d}.npy")
      for ki, key in enumerate(active_keys):
        handles[key]["var"][ti, :, :] = stack[ki]
  finally:
    for h_obj in handles.values():
      h_obj["ds"].close()

  shutil.rmtree(shards_dir, ignore_errors=True)
  print(f"=== Step C complete in {time.time() - t_start:.1f}s ===", flush=True)


if __name__ == "__main__":
  main()
