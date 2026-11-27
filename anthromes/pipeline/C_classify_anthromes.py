"""Step C: Classify Anthromes 12K across 218 timesteps with dynamic landarea and Holdridge biomes."""

from __future__ import annotations

import concurrent.futures
import json
import os
import shutil
import sys
import time
from pathlib import Path

import netCDF4 as nc
import numpy as np
from scipy.ndimage import label as nd_label

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "core"))
import hollinger_core as hc

if hasattr(sys.stdout, "reconfigure"):
  sys.stdout.reconfigure(line_buffering=True)

_WORKER_CTX: dict = {}


def _read_hyde_slice(ds: nc.Dataset, var_name: str, idx0: int, idx1: int, weight: float) -> np.ndarray:
  """Reads and optionally linearly interpolates a HYDE 3.5 2D slice, clamping invalid/negative to 0."""
  v = ds.variables[var_name]
  raw0 = v[idx0, :, :]
  arr0 = np.where(raw0.mask, 0.0, raw0.data) if hasattr(raw0, "mask") else np.asarray(raw0, dtype=np.float32)
  arr0 = np.where((arr0 < 0.0) | ~np.isfinite(arr0), 0.0, arr0).astype(np.float32)
  if idx0 == idx1 or weight <= 0.0:
    return arr0
  raw1 = v[idx1, :, :]
  arr1 = np.where(raw1.mask, 0.0, raw1.data) if hasattr(raw1, "mask") else np.asarray(raw1, dtype=np.float32)
  arr1 = np.where((arr1 < 0.0) | ~np.isfinite(arr1), 0.0, arr1).astype(np.float32)
  return ((1.0 - weight) * arr0 + weight * arr1).astype(np.float32)


def _init_worker(shards_dir_str: str) -> None:
  """Initialises per-process read-only datasets and static masks."""
  cfg = hc.load_json5("anthromes_config.json5")
  paths = hc.resolve_paths("anthromes")
  base_dir = paths["base_dir"]
  inter_dir = paths["intermediate_dir"]
  hyde_dir = (base_dir / cfg["sources"]["hyde35_rel_dir"]).resolve()

  npz = np.load(inter_dir / "static_masks.npz")
  supp = np.load(inter_dir / "ghsl_supplement.npz")
  mapping = json.loads((inter_dir / "hyde_time_mapping.json").read_text(encoding="utf-8"))

  ds_la = nc.Dataset(str(base_dir / "landarea" / "landarea.nc"), "r")
  ds_hlz = nc.Dataset(str(base_dir / "climate" / "holdridge_lifezones.nc"), "r")
  hyde_ds = {
    "cropland": nc.Dataset(str(hyde_dir / "cropland.nc"), "r"),
    "grazing_land": nc.Dataset(str(hyde_dir / "grazing_land.nc"), "r"),
    "urban_area": nc.Dataset(str(hyde_dir / "urban_area.nc"), "r"),
    "total_rice": nc.Dataset(str(hyde_dir / "total_rice.nc"), "r"),
    "total_irrigated": nc.Dataset(str(hyde_dir / "total_irrigated.nc"), "r"),
    "population": nc.Dataset(str(hyde_dir / "population.nc"), "r"),
  }

  _WORKER_CTX["cfg"] = cfg
  _WORKER_CTX["paths"] = paths
  _WORKER_CTX["shards_dir"] = Path(shards_dir_str)
  _WORKER_CTX["mapping"] = mapping
  _WORKER_CTX["pot_vill"] = npz["pot_vill"]
  _WORKER_CTX["potveg15_valid"] = npz["potveg15_valid"]
  _WORKER_CTX["trees_base_modern"] = npz["trees_base_modern"]
  _WORKER_CTX["trees_hlz_modern"] = npz["trees_hlz_modern"]
  _WORKER_CTX["ice_barren_modern"] = npz["ice_barren_modern"]
  _WORKER_CTX["la_hyde"] = npz["la_hyde"]
  _WORKER_CTX["hyde_dist_px"] = npz["hyde_dist_px"]
  _WORKER_CTX["hyde_nn_row"] = npz["hyde_nn_row"]
  _WORKER_CTX["hyde_nn_col"] = npz["hyde_nn_col"]
  _WORKER_CTX["supp_rows"] = supp["supp_rows"]
  _WORKER_CTX["supp_cols"] = supp["supp_cols"]
  _WORKER_CTX["ghsl_epochs"] = supp["ghsl_epochs"]
  _WORKER_CTX["supp_pop"] = supp["supp_pop"]
  _WORKER_CTX["supp_urban"] = supp["supp_urban"]
  _WORKER_CTX["supp_backcast_ratio"] = supp["supp_backcast_ratio"]
  _WORKER_CTX["hlz_tree_ids"] = np.asarray(cfg["classification"]["holdridge_tree_ids"], dtype=np.uint8)
  _WORKER_CTX["hlz_ice_id"] = int(cfg["classification"]["holdridge_ice_id"])
  _WORKER_CTX["ds_la"] = ds_la
  _WORKER_CTX["ds_hlz"] = ds_hlz
  _WORKER_CTX["hyde_ds"] = hyde_ds


def _classify_timestep(ti: int) -> tuple[int, int, float]:
  """Classifies a single timestep, writes its GeoPNG, and saves a shard for NetCDF assembly."""
  t0 = time.time()
  ctx = _WORKER_CTX
  cfg = ctx["cfg"]
  cls_cfg = cfg["classification"]
  m = ctx["mapping"][ti]
  y = int(m["year"])
  idx0 = int(m["idx0"])
  idx1 = int(m["idx1"])
  weight = float(m["weight"])

  la = np.asarray(ctx["ds_la"].variables["land_area"][ti, :, :], dtype=np.float32)
  hlz = np.asarray(ctx["ds_hlz"].variables["holdridge_lifezones"][ti, :, :], dtype=np.uint8)

  out = np.zeros(la.shape, dtype=np.uint8)
  out[la < 0.0] = 63

  land_mask = la > 0.0
  if np.any(land_mask):
    la_hyde_2d = ctx["la_hyde"]
    hyde_ds = ctx["hyde_ds"]
    crops_2d_km2 = _read_hyde_slice(hyde_ds["cropland"], "cropland", idx0, idx1, weight)
    grazing_2d_km2 = _read_hyde_slice(hyde_ds["grazing_land"], "grazing_land", idx0, idx1, weight)
    urban_2d_km2 = _read_hyde_slice(hyde_ds["urban_area"], "urban_area", idx0, idx1, weight)
    rice_2d_km2 = _read_hyde_slice(hyde_ds["total_rice"], "total_rice", idx0, idx1, weight)
    irri_2d_km2 = _read_hyde_slice(hyde_ds["total_irrigated"], "total_irrigated", idx0, idx1, weight)
    pop_2d_cap = _read_hyde_slice(hyde_ds["population"], "population", idx0, idx1, weight)

    la_1d = la[land_mask]
    la_hyde_1d = la_hyde_2d[land_mask]
    denom = np.maximum(la_1d, la_hyde_1d)
    valid_denom = denom > 0.0
    safe_denom = np.where(valid_denom, denom, 1.0)

    crops = np.where(valid_denom, np.clip(crops_2d_km2[land_mask] / safe_denom, 0.0, 1.0), 0.0)
    grazing = np.where(valid_denom, np.clip(grazing_2d_km2[land_mask] / safe_denom, 0.0, 1.0), 0.0)
    urban = np.where(valid_denom, np.clip(urban_2d_km2[land_mask] / safe_denom, 0.0, 1.0), 0.0)
    rice = np.where(valid_denom, np.clip(rice_2d_km2[land_mask] / safe_denom, 0.0, 1.0), 0.0)
    irrigation = np.where(valid_denom, np.clip(irri_2d_km2[land_mask] / safe_denom, 0.0, 1.0), 0.0)
    popd = np.where(valid_denom, pop_2d_cap[land_mask] / safe_denom, 0.0)

    # Coastal & reclaimed intensive-density extrapolation + GHSL/Gapminder island fusion on la_hyde == 0
    fusion_cfg = cfg["coastal_and_island_fusion"]
    nc_mask_2d = land_mask & (la_hyde_2d <= 0.0)

    ghsl_pop_2d = np.zeros(la.shape, dtype=np.float32)
    ghsl_urb_2d = np.zeros(la.shape, dtype=bool)
    supp_r = ctx["supp_rows"]
    supp_c = ctx["supp_cols"]
    ghsl_ep = ctx["ghsl_epochs"]
    if y >= int(ghsl_ep[0]):
      pos = int(np.searchsorted(ghsl_ep, y))
      if pos <= 0:
        sp_vals = ctx["supp_pop"][0]
        su_vals = ctx["supp_urban"][0]
      elif pos >= len(ghsl_ep):
        sp_vals = ctx["supp_pop"][-1]
        su_vals = ctx["supp_urban"][-1]
      elif int(ghsl_ep[pos]) == y:
        sp_vals = ctx["supp_pop"][pos]
        su_vals = ctx["supp_urban"][pos]
      else:
        e0, e1 = int(ghsl_ep[pos - 1]), int(ghsl_ep[pos])
        wt = float(y - e0) / float(e1 - e0)
        sp_vals = (1.0 - wt) * ctx["supp_pop"][pos - 1] + wt * ctx["supp_pop"][pos]
        su_vals = ctx["supp_urban"][pos - 1] | ctx["supp_urban"][pos]
    else:
      sp_vals = ctx["supp_pop"][0] * ctx["supp_backcast_ratio"][ti]
      su_vals = np.zeros(len(supp_r), dtype=bool)

    ghsl_pop_2d[supp_r, supp_c] = sp_vals
    ghsl_urb_2d[supp_r, supp_c] = su_vals

    if np.any(nc_mask_2d):
      inv_la_hyde = np.where(la_hyde_2d > 0.0, 1.0 / np.maximum(la_hyde_2d, 1e-6), 0.0)
      crops_d2d = np.clip(crops_2d_km2 * inv_la_hyde, 0.0, 1.0)
      grazing_d2d = np.clip(grazing_2d_km2 * inv_la_hyde, 0.0, 1.0)
      urban_d2d = np.clip(urban_2d_km2 * inv_la_hyde, 0.0, 1.0)
      rice_d2d = np.clip(rice_2d_km2 * inv_la_hyde, 0.0, 1.0)
      irri_d2d = np.clip(irri_2d_km2 * inv_la_hyde, 0.0, 1.0)
      popd_d2d = pop_2d_cap * inv_la_hyde

      nn_r = ctx["hyde_nn_row"][nc_mask_2d]
      nn_c = ctx["hyde_nn_col"][nc_mask_2d]
      d_px = ctx["hyde_dist_px"][nc_mask_2d]

      comp_labels, _ = nd_label(land_mask, structure=np.ones((3, 3), dtype=np.int8))
      hyde_comps = np.unique(comp_labels[la_hyde_2d > 0.0])
      same_landmass = np.isin(comp_labels[nc_mask_2d], hyde_comps)

      imm_px = float(fusion_cfg["coastal_immediate_px"])
      fringe_px = float(fusion_cfg["coastal_fringe_px"])
      paleo_max = float(fusion_cfg["paleoshelf_max_px"])
      paleo_decay = float(fusion_cfg["paleoshelf_decay_px"])
      nc_ghsl_pop = ghsl_pop_2d[nc_mask_2d]

      decay = np.zeros(d_px.shape, dtype=np.float32)
      decay[d_px <= imm_px] = 1.0
      m_fringe = (d_px > imm_px) & (d_px <= fringe_px)
      if y >= 1975:
        decay[m_fringe & (nc_ghsl_pop > 0.0)] = 0.65
      else:
        decay[m_fringe] = 0.65

      if y < 1975:
        m_paleo = (d_px > fringe_px) & (d_px <= paleo_max) & same_landmass
        decay[m_paleo] = np.exp(-(d_px[m_paleo] - imm_px) / paleo_decay).astype(np.float32)

      nc_in_land1d = la_hyde_1d <= 0.0
      crops[nc_in_land1d] = crops_d2d[nn_r, nn_c] * decay
      grazing[nc_in_land1d] = grazing_d2d[nn_r, nn_c] * decay
      urban[nc_in_land1d] = urban_d2d[nn_r, nn_c] * decay
      rice[nc_in_land1d] = rice_d2d[nn_r, nn_c] * decay
      irrigation[nc_in_land1d] = irri_d2d[nn_r, nn_c] * decay
      popd[nc_in_land1d] = popd_d2d[nn_r, nn_c] * decay

      min_island_la = float(fusion_cfg["island_min_land_area_km2"])
      has_ghsl = nc_ghsl_pop > 0.0
      if np.any(has_ghsl):
        popd_ghsl = nc_ghsl_pop / np.maximum(la_1d[nc_in_land1d], min_island_la)
        cur_popd = popd[nc_in_land1d]
        cur_popd[has_ghsl] = np.maximum(cur_popd[has_ghsl], popd_ghsl[has_ghsl])
        popd[nc_in_land1d] = cur_popd

      nc_ghsl_urb = ghsl_urb_2d[nc_mask_2d]
      if np.any(nc_ghsl_urb):
        cur_urb = urban[nc_in_land1d]
        cur_urb[nc_ghsl_urb] = np.maximum(cur_urb[nc_ghsl_urb], 0.25)
        urban[nc_in_land1d] = cur_urb

    # Proportionally rescale crops + grazing + urban if sum exceeds 1.0
    total_use = crops + grazing + urban
    over = total_use > 1.0
    if np.any(over):
      scale = 1.0 / total_use[over]
      crops[over] *= scale
      grazing[over] *= scale
      urban[over] *= scale

    rice = np.minimum(rice, crops)
    irrigation = np.minimum(irrigation, crops)
    used = np.minimum(1.0, crops + grazing + urban)

    # Canonical Ellis et al. (2020, 2021) potveg15 woodland mask on valid non-ice land;
    # cleaned Holdridge tree mask strictly on exposed paleoshelves (~pv_valid) & deglaciated modern-ice margins
    hlz_1d = hlz[land_mask]
    hlz_ice_id = ctx["hlz_ice_id"]
    trees_hlz_1d = np.isin(hlz_1d, ctx["hlz_tree_ids"])
    pv_valid_1d = ctx["potveg15_valid"][land_mask]
    trees_base_1d = ctx["trees_base_modern"][land_mask]
    ice_mod_1d = ctx["ice_barren_modern"][land_mask]

    trees = np.where(pv_valid_1d & (~ice_mod_1d), trees_base_1d, trees_hlz_1d)
    ice_barren = (
      (ice_mod_1d & (hlz_1d == hlz_ice_id))
      | ((~pv_valid_1d) & (hlz_1d == hlz_ice_id))
    )
    pv_1d = ctx["pot_vill"][land_mask]

    # Anthromes 12K decision tree (Ellis et al. 2020, 2021) in strict top-to-bottom order
    urban_min = float(cls_cfg["urban_frac_min"])
    urban_popd_min = float(cls_cfg["urban_popd_min"])
    dense_popd_min = float(cls_cfg["dense_popd_min"])
    res_popd_min = float(cls_cfg["residential_popd_min"])
    pop_popd_min = float(cls_cfg["populated_popd_min"])
    used_min = float(cls_cfg["used_frac_min"])
    rice_min = float(cls_cfg["rice_frac_min"])
    irri_min = float(cls_cfg["irrigation_frac_min"])
    crops_min = float(cls_cfg["crops_frac_min"])
    grazing_min = float(cls_cfg["grazing_frac_min"])

    out_1d = np.zeros(la_1d.shape, dtype=np.uint8)

    # 1. Urban (11)
    out_1d[(urban >= urban_min) | (popd >= urban_popd_min)] = 11

    # 2. Villages & Mixed Settlements (12, 21..24)
    m_dense = (out_1d == 0) & (popd >= dense_popd_min)
    m_vill = m_dense & pv_1d
    out_1d[m_dense & (~pv_1d)] = 12
    out_1d[m_vill & (rice >= rice_min)] = 21
    out_1d[(out_1d == 0) & m_vill & (irrigation >= irri_min)] = 22
    out_1d[(out_1d == 0) & m_vill & (crops >= crops_min)] = 23
    out_1d[(out_1d == 0) & m_vill & (grazing >= grazing_min)] = 24
    out_1d[(out_1d == 0) & m_vill] = 12

    # 3. Croplands (31..34)
    m_crop = (out_1d == 0) & (crops >= crops_min)
    out_1d[m_crop & (popd >= res_popd_min) & (irrigation >= irri_min)] = 31
    out_1d[m_crop & (popd >= res_popd_min) & (irrigation < irri_min)] = 32
    out_1d[m_crop & (popd >= pop_popd_min) & (popd < res_popd_min)] = 33
    out_1d[m_crop & (popd < pop_popd_min)] = 34

    # 4. Rangelands (41..43)
    m_range = (out_1d == 0) & (grazing >= grazing_min)
    out_1d[m_range & (popd >= res_popd_min)] = 41
    out_1d[m_range & (popd >= pop_popd_min) & (popd < res_popd_min)] = 42
    out_1d[m_range & (popd < pop_popd_min)] = 43

    # 5. Mixed used >= 20% (routes to Remote croplands 34 or Remote rangelands 43)
    m_used = (out_1d == 0) & (used >= used_min)
    out_1d[m_used & (crops > grazing)] = 34
    out_1d[m_used & (crops <= grazing)] = 43

    # 6. Cultured (51..54)
    m_cult = (out_1d == 0) & (popd > 0.0)
    out_1d[m_cult & trees & (popd >= res_popd_min)] = 51
    out_1d[m_cult & trees & (popd >= pop_popd_min) & (popd < res_popd_min)] = 52
    out_1d[m_cult & trees & (popd < pop_popd_min)] = 53
    out_1d[m_cult & (~trees)] = 54

    # 7. Wildlands (61..63)
    m_wild = out_1d == 0
    out_1d[m_wild & ice_barren] = 63
    out_1d[m_wild & (~ice_barren) & trees] = 61
    out_1d[m_wild & (~ice_barren) & (~trees)] = 62

    out[land_mask] = out_1d

  rasters_dir = ctx["paths"]["rasters_dir"]
  prefix = cfg["encoding"]["geopng_prefix"]
  png_path = rasters_dir / f"{prefix}{y}.png"
  hc.encode_uint32_geopng(png_path, out.astype(np.uint32))

  shard_path = ctx["shards_dir"] / f"t_{ti:03d}.npy"
  np.save(shard_path, out)
  return ti, y, time.time() - t0


def main() -> None:
  t_start = time.time()
  paths = hc.resolve_paths("anthromes")
  cfg = hc.load_json5("anthromes_config.json5")
  mcfg = hc.load_json5("model_config.json5")
  n_workers = min(8, int(mcfg["compute"]["max_workers"]), os.cpu_count() or 4)
  years = hc.get_timestep_years()

  lat_5m, lon_5m = hc.get_grid_coords(4320, 2160)
  cell_area = hc.compute_cell_area_by_row(4320, 2160)

  rasters_dir = paths["rasters_dir"]
  rasters_dir.mkdir(parents=True, exist_ok=True)
  shards_dir = paths["intermediate_dir"] / "shards"
  shards_dir.mkdir(parents=True, exist_ok=True)

  print(f"=== Step C: Classifying {len(years)} timesteps with {n_workers} worker processes ===")
  done = 0
  with concurrent.futures.ProcessPoolExecutor(
    max_workers=n_workers,
    initializer=_init_worker,
    initargs=(str(shards_dir),),
  ) as pool:
    for ti, y, dt in pool.map(_classify_timestep, range(len(years))):
      done += 1
      if done % 20 == 0 or done == 1 or done == len(years):
        print(f"  [{done:3d}/{len(years)}] year={y:6d} ({dt:.2f}s/step)")

  print("=== Step C: Assembling compressed CF-1.8 NetCDF-4 dataset ===")
  enc = cfg["encoding"]
  classes = cfg["classes"]
  flag_values = [int(c["id"]) for c in classes]
  flag_meanings = " ".join(c["code"] for c in classes)
  nc_path = paths["indicator_dir"] / enc["nc_filename"]

  ds, data_var = hc.create_categorical_netcdf4(
    nc_path=nc_path,
    var_name=enc["var_name"],
    years=years,
    lat=lat_5m,
    lon=lon_5m,
    cell_area_km2=cell_area,
    flag_values=flag_values,
    flag_meanings=flag_meanings,
    global_attrs={
      "title": "Anthromes 12K classification (Ellis et al. 2020, 2021) coupled with dynamic landarea and Holdridge biomes",
      "institution": "SVEA / landuse_Hollinger",
      "source": "Anthromes 12K DGG v1 (Ellis et al. 2020, 2021) + HYDE 3.5 + landarea.nc + holdridge_lifezones.nc",
      "history": "Generated by anthromes/pipeline/C_classify_anthromes.py",
    },
    var_long_name="Anthromes 12K anthropogenic biome classification",
    chunk_shape=tuple(int(x) for x in enc["netcdf_chunk_shape"]),
    complevel=int(enc["netcdf_compression_level"]),
  )
  try:
    for ti in range(len(years)):
      data_var[ti, :, :] = np.load(shards_dir / f"t_{ti:03d}.npy")
  finally:
    ds.close()

  shutil.rmtree(shards_dir, ignore_errors=True)
  print(f"=== Step C complete in {time.time() - t_start:.1f}s ===")


if __name__ == "__main__":
  main()
