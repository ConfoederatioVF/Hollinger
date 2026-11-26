"""Step B: Build dynamic biome masks, extended pot_vill paleoshelf mask, HYDE land_area, and HYDE time mapping."""

from __future__ import annotations

import csv
import json
import sys
import time
from pathlib import Path

import netCDF4 as nc
import numpy as np
from scipy.ndimage import distance_transform_edt
import tifffile

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "core"))
import hollinger_core as hc

if hasattr(sys.stdout, "reconfigure"):
  sys.stdout.reconfigure(line_buffering=True)


def get_hyde_years() -> list[int]:
  """Returns the canonical 128 HYDE 3.5 signed calendar years (-10000..2025, 0 = 1 AD)."""
  return (
    [-10000 + i * 1000 for i in range(10)]
    + [i * 100 for i in range(18)]
    + [1710 + i * 10 for i in range(25)]
    + [1951 + i for i in range(75)]
  )


def load_esri_ascii_grid(asc_path: Path) -> np.ndarray:
  """Loads an ESRI ASCII (.asc) grid, returning float32 with NODATA replaced by 0.0."""
  header = {}
  with open(asc_path, "r", encoding="utf-8") as f:
    for _ in range(6):
      parts = f.readline().strip().split()
      header[parts[0].lower()] = float(parts[1])
  ncols = int(header["ncols"])
  nrows = int(header["nrows"])
  nodata = header.get("nodata_value", -9999.0)
  arr = np.loadtxt(str(asc_path), skiprows=6, dtype=np.float32)
  if arr.shape != (nrows, ncols):
    raise ValueError(f"Expected {(nrows, ncols)} in {asc_path}, got {arr.shape}")
  arr = np.where((arr == nodata) | (arr < 0.0) | ~np.isfinite(arr), 0.0, arr).astype(np.float32)
  return arr


def load_geotiff_2d(tif_path: Path) -> np.ndarray:
  """Loads a 2D GeoTIFF and validates its shape is (2160, 4320)."""
  arr = tifffile.imread(str(tif_path))
  if arr.ndim != 2 or arr.shape != (2160, 4320):
    raise ValueError(f"Expected (2160, 4320) in {tif_path}, got {arr.shape}")
  return arr


def build_hyde_time_mapping(hollinger_years: list[int], hyde_years: list[int]) -> list[dict]:
  """Maps each of the 218 Hollinger years to bounding HYDE 3.5 indices and linear weight."""
  hyde_idx = {y: i for i, y in enumerate(hyde_years)}
  mapping = []
  for y in hollinger_years:
    if y in hyde_idx:
      idx = hyde_idx[y]
      mapping.append({"year": int(y), "idx0": int(idx), "idx1": int(idx), "weight": 0.0})
    else:
      pos = int(np.searchsorted(hyde_years, y))
      if pos <= 0 or pos >= len(hyde_years):
        raise ValueError(f"Hollinger year {y} is outside HYDE range [{hyde_years[0]}, {hyde_years[-1]}]")
      i0 = pos - 1
      i1 = pos
      y0 = hyde_years[i0]
      y1 = hyde_years[i1]
      w = float(y - y0) / float(y1 - y0)
      mapping.append({"year": int(y), "idx0": int(i0), "idx1": int(i1), "weight": w})
  return mapping


def main() -> None:
  t0 = time.time()
  cfg = hc.load_json5("anthromes_config.json5")
  paths = hc.resolve_paths("anthromes")
  raw_dir = paths["raw_dir"]
  grids_dir = raw_dir / "supporting_5m_grids"
  inter_dir = paths["intermediate_dir"]
  inter_dir.mkdir(parents=True, exist_ok=True)
  base_dir = paths["base_dir"]
  cls_cfg = cfg["classification"]
  src_cfg = cfg["sources"]

  # 1. Extended pot_vill mask via Euclidean nearest-neighbour on valid potvill20.tif pixels
  print("Building extended pot_vill mask from potvill20.tif...")
  potvill_raw = load_geotiff_2d(grids_dir / "potvill20.tif")
  potvill_valid = (potvill_raw == 0) | (potvill_raw == 1)
  potvill_bool = potvill_raw == 1

  # Ensure Greenland (excluding Iceland) and Antarctica (+ their shelves) are False
  lat_1d, lon_1d = hc.get_grid_coords(4320, 2160)
  lat_2d = lat_1d[:, None]
  lon_2d = lon_1d[None, :]
  iceland_box = (lat_2d >= 63.0) & (lat_2d <= 67.0) & (lon_2d >= -25.0) & (lon_2d <= -13.0)
  greenland_box = (lat_2d >= 59.0) & (lon_2d >= -74.0) & (lon_2d <= -11.0) & (~iceland_box)
  antarctica_rows = lat_1d <= -60.0
  potvill_bool[greenland_box] = False
  potvill_valid[antarctica_rows, :] = True
  potvill_bool[antarctica_rows, :] = False

  _, nearest_idx = distance_transform_edt(~potvill_valid, return_indices=True)
  pot_vill = potvill_bool[nearest_idx[0], nearest_idx[1]].astype(bool)
  pot_vill[greenland_box] = False
  pot_vill[antarctica_rows, :] = False
  print(f"  potvill20 valid pixels: {int(potvill_valid.sum()):,}; extended pot_vill True pixels: {int(pot_vill.sum()):,}")

  # 2. potveg15 modern woodland/forest (1..8) and ice/barren (15)
  print("Loading potveg15.tif...")
  potveg_raw = load_geotiff_2d(grids_dir / "potveg15.tif")
  potveg15_valid = np.isfinite(potveg_raw) & (potveg_raw >= 1) & (potveg_raw <= 15)
  potveg15_u8 = np.where(potveg15_valid, potveg_raw, 0).astype(np.uint8)
  tree_ids_pv = np.asarray(cls_cfg["potveg15_tree_ids"], dtype=np.uint8)
  trees_base_modern = potveg15_valid & np.isin(potveg15_u8, tree_ids_pv)
  ice_barren_modern = potveg15_valid & (potveg15_u8 == int(cls_cfg["potveg15_ice_id"]))
  print(
    f"  potveg15 valid: {int(potveg15_valid.sum()):,}, "
    f"trees_base_modern: {int(trees_base_modern.sum()):,}, "
    f"ice_barren_modern: {int(ice_barren_modern.sum()):,}"
  )

  # 3. Holdridge Life Zones modern anchor (1985)
  anchor_year = int(cls_cfg["biome_anchor_year"])
  hlz_path = base_dir / "climate" / "holdridge_lifezones.nc"
  print(f"Extracting Holdridge Life Zones anchor at {anchor_year}...")
  with nc.Dataset(str(hlz_path), "r") as ds_hlz:
    hlz_years = list(np.asarray(ds_hlz.variables["time"][:], dtype=np.int32).tolist())
    if anchor_year not in hlz_years:
      raise ValueError(f"Anchor year {anchor_year} not in holdridge_lifezones.nc")
    t_anchor = hlz_years.index(anchor_year)
    hlz_modern = np.asarray(ds_hlz.variables["holdridge_lifezones"][t_anchor, :, :], dtype=np.uint8)

  hlz_tree_ids = np.asarray(cls_cfg["holdridge_tree_ids"], dtype=np.uint8)
  trees_hlz_modern = np.isin(hlz_modern, hlz_tree_ids)
  print(f"  trees_hlz_modern ({anchor_year}): {int(trees_hlz_modern.sum()):,}")

  # 4. HYDE land_area_per_cell.asc
  hyde_meta_dir = (base_dir / src_cfg["hyde_meta_rel_dir"]).resolve()
  asc_candidates = [
    hyde_meta_dir / "land_area_per_cell.asc",
    hyde_meta_dir.parent / "land_area_per_cell.asc",
  ]
  asc_path = next((p for p in asc_candidates if p.exists()), None)
  if asc_path is None:
    raise FileNotFoundError(f"land_area_per_cell.asc not found in {[str(p) for p in asc_candidates]}")
  print(f"Loading HYDE land_area_per_cell.asc from {asc_path}...")
  la_hyde = load_esri_ascii_grid(asc_path)
  print(f"  la_hyde > 0 cells: {int((la_hyde > 0).sum()):,}, sum: {float(la_hyde.sum()):,.1f} km2")

  # Save static intermediates + nearest HYDE land distance/indices
  print("Computing Euclidean nearest-neighbour distance to HYDE land_area_per_cell...")
  hyde_dist_px, hyde_nn = distance_transform_edt(la_hyde <= 0.0, return_indices=True)
  hyde_dist_px = hyde_dist_px.astype(np.float32)
  hyde_nn_row = hyde_nn[0].astype(np.int16)
  hyde_nn_col = hyde_nn[1].astype(np.int16)

  npz_path = inter_dir / "static_masks.npz"
  np.savez_compressed(
    npz_path,
    pot_vill=pot_vill,
    potveg15_raw=potveg15_u8,
    potveg15_valid=potveg15_valid,
    trees_base_modern=trees_base_modern,
    trees_hlz_modern=trees_hlz_modern,
    ice_barren_modern=ice_barren_modern,
    la_hyde=la_hyde,
    hyde_dist_px=hyde_dist_px,
    hyde_nn_row=hyde_nn_row,
    hyde_nn_col=hyde_nn_col,
  )
  print(f"  [OK] Saved {npz_path}")

  # 5. Precompute HYDE 3.5 timestep mapping for all 218 Hollinger years
  hollinger_years = hc.get_timestep_years()
  hyde_years = get_hyde_years()
  hyde_dir = (base_dir / src_cfg["hyde35_rel_dir"]).resolve()
  with nc.Dataset(str(hyde_dir / "cropland.nc"), "r") as ds_c:
    if len(ds_c.variables["time"]) != len(hyde_years):
      raise ValueError(f"HYDE cropland.nc has {len(ds_c.variables['time'])} steps, expected {len(hyde_years)}")

  mapping = build_hyde_time_mapping(hollinger_years, hyde_years)
  exact_count = sum(1 for m in mapping if m["idx0"] == m["idx1"])
  interp_count = len(mapping) - exact_count
  map_path = inter_dir / "hyde_time_mapping.json"
  map_path.write_text(json.dumps(mapping, indent=2), encoding="utf-8")
  print(f"  [OK] Saved {map_path} ({exact_count} exact + {interp_count} interpolated = {len(mapping)} timesteps)")

  # 6. Precompute GHSL + Gapminder coastal & island population supplement (for la_hyde == 0 cells)
  fusion_cfg = cfg["coastal_and_island_fusion"]
  ghsl_epochs = np.asarray(fusion_cfg["ghsl_epochs"], dtype=np.int32)
  ghsl_scale = float(fusion_cfg["ghsl_pop_scale_factor"])
  ghsl_pop_dir = (base_dir / src_cfg["ghsl_pop_geotiffs_rel_dir"]).resolve()
  ghsl_urb_dir = (base_dir / src_cfg["ghsl_urban_rasters_rel_dir"]).resolve()

  print(f"Loading {len(ghsl_epochs)} GHSL population & urban epochs for island/coastal fusion...")
  n_ep = len(ghsl_epochs)
  ghsl_pop_stack = np.zeros((n_ep, 2160, 4320), dtype=np.float32)
  ghsl_urb_stack = np.zeros((n_ep, 2160, 4320), dtype=bool)
  for ei, ep_yr in enumerate(ghsl_epochs):
    tif_path = ghsl_pop_dir / f"GHS_POP_E{int(ep_yr)}.tif"
    raw_g = tifffile.imread(str(tif_path))
    if raw_g.shape == (2138, 4320):
      # E1975.tif starts at 89.099583°N (row offset 11 from 90°N)
      aligned = np.zeros((2160, 4320), dtype=np.float32)
      aligned[11:11 + 2138, :] = np.where(np.isfinite(raw_g) & (raw_g > 0.0), raw_g, 0.0).astype(np.float32)
    elif raw_g.shape == (2160, 4320):
      aligned = np.where(np.isfinite(raw_g) & (raw_g > 0.0), raw_g, 0.0).astype(np.float32)
    else:
      raise ValueError(f"Unexpected GHSL shape {raw_g.shape} in {tif_path}")
    ghsl_pop_stack[ei] = aligned * ghsl_scale

    urb_path = ghsl_urb_dir / f"GHS_urban_{int(ep_yr)}.png"
    urb_u32 = hc.decode_uint32_geopng(urb_path)
    ghsl_urb_stack[ei] = urb_u32 > 0

  supp_mask = (la_hyde <= 0.0) & ((ghsl_pop_stack.max(axis=0) > 0.0) | ghsl_urb_stack.any(axis=0))
  supp_rows, supp_cols = np.where(supp_mask)
  supp_rows = supp_rows.astype(np.int16)
  supp_cols = supp_cols.astype(np.int16)
  supp_pop = ghsl_pop_stack[:, supp_rows, supp_cols].astype(np.float32)
  supp_urban = ghsl_urb_stack[:, supp_rows, supp_cols].astype(bool)
  del ghsl_pop_stack, ghsl_urb_stack
  print(f"  GHSL supplement pixels (la_hyde==0 & GHSL>0): {len(supp_rows):,}")

  # Load Gapminder country series & HYDE regional population totals for pre-1975 backcasting
  gm_path = (base_dir / src_cfg["gapminder_pop_csv_rel_path"]).resolve()
  gm_series: dict[str, tuple[np.ndarray, np.ndarray]] = {}
  raw_gm: dict[str, list[tuple[int, float]]] = {}
  with open(gm_path, "r", encoding="utf-8-sig") as f_gm:
    for row in csv.DictReader(f_gm):
      iso3 = row["geo"].strip().lower()
      yr_gm = int(row["year"])
      pop_gm = float(row["Total population"])
      raw_gm.setdefault(iso3, []).append((yr_gm, pop_gm))
  for iso3, pairs in raw_gm.items():
    pairs.sort(key=lambda p: p[0])
    gm_series[iso3] = (
      np.asarray([p[0] for p in pairs], dtype=np.float64),
      np.asarray([p[1] for p in pairs], dtype=np.float64),
    )

  # Load HYDE 3.5 population grid at 1975 and all Hollinger years for regional backcast scaling
  with nc.Dataset(str(hyde_dir / "population.nc"), "r") as ds_hpop:
    hyde_1975_idx = hyde_years.index(1975)
    hpop_all = np.zeros((len(hyde_years), 2160, 4320), dtype=np.float32)
    for hi in range(len(hyde_years)):
      v_hp = ds_hpop.variables["population"][hi, :, :]
      a_hp = np.where(v_hp.mask, 0.0, v_hp.data) if hasattr(v_hp, "mask") else np.asarray(v_hp, dtype=np.float32)
      hpop_all[hi] = np.where((a_hp < 0.0) | ~np.isfinite(a_hp), 0.0, a_hp).astype(np.float32)

  # Precompute coarse 5-degree macro-region HYDE population ratio for fallback backcasting
  # Bin size: 60 rows x 60 cols (5 deg x 5 deg)
  n_supp = len(supp_rows)
  supp_backcast_ratio = np.ones((len(hollinger_years), n_supp), dtype=np.float32)
  supp_lat = lat_1d[supp_rows]
  supp_lon = lon_1d[supp_cols]

  # Precompute 5-degree window HYDE pop time series per unique bin (with global fallback for open-ocean atolls)
  global_hyde_ts = hpop_all.sum(axis=(1, 2)).astype(np.float64)
  global_hyde_ratio = np.clip(global_hyde_ts / max(float(global_hyde_ts[hyde_1975_idx]), 1.0), 0.0, 2.0)

  bin_r = (supp_rows // 60).astype(np.int32)
  bin_c = (supp_cols // 60).astype(np.int32)
  unique_bins = set(zip(bin_r.tolist(), bin_c.tolist()))
  bin_ratios: dict[tuple[int, int], np.ndarray] = {}
  for br, bc in unique_bins:
    r0, r1 = max(0, br * 60 - 60), min(2160, (br + 1) * 60 + 60)
    c0, c1 = max(0, bc * 60 - 60), min(4320, (bc + 1) * 60 + 60)
    win_ts = hpop_all[:, r0:r1, c0:c1].sum(axis=(1, 2)).astype(np.float64)
    if float(win_ts[hyde_1975_idx]) < 100.0:
      ratio_hyde = global_hyde_ratio
    else:
      ref_75 = max(float(win_ts[hyde_1975_idx]), 1.0)
      ratio_hyde = np.clip(win_ts / ref_75, 0.0, 2.0)
    h_ratios = np.ones(len(hollinger_years), dtype=np.float32)
    for ti, m_t in enumerate(mapping):
      if m_t["year"] < 1975:
        i0, i1, w = m_t["idx0"], m_t["idx1"], m_t["weight"]
        h_ratios[ti] = float((1.0 - w) * ratio_hyde[i0] + w * ratio_hyde[i1])
    bin_ratios[(br, bc)] = h_ratios

  for si in range(n_supp):
    supp_backcast_ratio[:, si] = bin_ratios[(int(bin_r[si]), int(bin_c[si]))]

  # Override with Gapminder + historical settlement onset for explicit island regions
  hollinger_yrs_arr = np.asarray(hollinger_years, dtype=np.float64)
  for reg in fusion_cfg["island_regions"]:
    in_reg = np.zeros(n_supp, dtype=bool)
    for la0, la1, lo0, lo1 in reg["boxes"]:
      in_reg |= (supp_lat >= la0) & (supp_lat <= la1) & (supp_lon >= lo0) & (supp_lon <= lo1)
    if not np.any(in_reg):
      continue
    settle_yr = float(reg["settlement_year"])
    iso3 = reg["iso3"].lower()
    if iso3 in gm_series:
      gm_y, gm_p = gm_series[iso3]
      p_1975 = float(np.interp(1975.0, gm_y, gm_p))
      gm_min_y = float(gm_y[0])
      p_min = float(gm_p[0])
      r_min = min(1.0, max(0.0, p_min / max(p_1975, 1.0)))
    else:
      gm_y, gm_p, p_1975, gm_min_y, r_min = None, None, 1.0, 1975.0, 1.0

    # Representative regional HYDE curve from first matching pixel's macro-bin
    first_idx = int(np.where(in_reg)[0][0])
    rep_hyde_curve = bin_ratios[(int(bin_r[first_idx]), int(bin_c[first_idx]))].astype(np.float64)
    if gm_y is not None:
      hyde_at_gm_min = max(float(np.interp(gm_min_y, hollinger_yrs_arr, rep_hyde_curve)), 1e-4)

    reg_curve = np.ones(len(hollinger_years), dtype=np.float32)
    for ti, hyr in enumerate(hollinger_years):
      if hyr >= 1975:
        reg_curve[ti] = 1.0
      elif hyr < settle_yr:
        reg_curve[ti] = 0.0
      elif gm_y is not None and hyr >= gm_min_y:
        reg_curve[ti] = float(np.clip(np.interp(float(hyr), gm_y, gm_p) / max(p_1975, 1.0), 0.0, 2.0))
      else:
        # Between settle_yr and gm_min_y: smooth ramp from 0 at settle_yr scaled by regional HYDE growth
        span = max(1.0, gm_min_y - settle_yr)
        ramp = min(1.0, (float(hyr) - settle_yr) / min(200.0, span))
        if gm_y is not None:
          val = r_min * min(1.0, float(rep_hyde_curve[ti]) / hyde_at_gm_min) * ramp
        else:
          val = float(rep_hyde_curve[ti]) * ramp
        reg_curve[ti] = float(np.clip(val, 0.0, 1.0))

    supp_backcast_ratio[:, in_reg] = reg_curve[:, None]

  del hpop_all
  supp_npz_path = inter_dir / "ghsl_supplement.npz"
  np.savez_compressed(
    supp_npz_path,
    supp_rows=supp_rows,
    supp_cols=supp_cols,
    ghsl_epochs=ghsl_epochs,
    supp_pop=supp_pop,
    supp_urban=supp_urban,
    supp_backcast_ratio=supp_backcast_ratio,
  )
  print(f"  [OK] Saved {supp_npz_path}")
  print(f"=== Step B complete in {time.time() - t0:.1f}s ===")


if __name__ == "__main__":
  main()
