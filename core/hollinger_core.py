"""Shared WGS84 grid, GeoPNG I/O, CF-1.8 NetCDF-4 writer, and spot-check renderer for landuse_Hollinger."""

import math
import os
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import netCDF4 as nc
import numpy as np
from PIL import Image


def get_base_dir() -> Path:
  return Path(__file__).resolve().parent.parent


def resolve_paths(indicator: str = "landarea") -> dict:
  base_dir = get_base_dir()
  common_dir = base_dir / "common"
  core_dir = base_dir / "core"
  tests_dir = base_dir / "tests"
  ind_dir = base_dir / indicator
  return {
    "base_dir": base_dir,
    "common_dir": common_dir,
    "core_dir": core_dir,
    "tests_dir": tests_dir,
    "indicator_dir": ind_dir,
    "raw_dir": ind_dir / "1.raw_data",
    "intermediate_dir": ind_dir / "2.intermediates",
    "rasters_dir": ind_dir / "3.rasters",
    "pipeline_dir": ind_dir / "pipeline",
    "nc_path": ind_dir / f"{indicator}.nc",
    "diagnostics_dir": tests_dir / "diagnostics" / indicator,
    "exports_dir": base_dir / "exports",
    "renders_dir": base_dir / "exports" / indicator,
    "stop_hooks_dir": tests_dir / ".stop_hooks",
  }


def load_json5(filename: str) -> dict:
  """Minimal JSON5 parser (strips // and /* */ comments, trailing commas, unquoted keys)."""
  import json
  p = get_base_dir() / "common" / filename
  txt = p.read_text(encoding="utf-8")
  # Remove /* ... */ block comments
  txt = re.sub(r"/\*[\s\S]*?\*/", "", txt)
  # Remove // line comments outside strings
  out_lines = []
  for line in txt.splitlines():
    in_str = False
    quote_ch = None
    esc = False
    cut = len(line)
    for idx, ch in enumerate(line):
      if esc:
        esc = False
        continue
      if ch == "\\":
        esc = True
        continue
      if in_str:
        if ch == quote_ch:
          in_str = False
      else:
        if ch in ('"', "'"):
          in_str = True
          quote_ch = ch
        elif ch == "/" and idx + 1 < len(line) and line[idx + 1] == "/":
          cut = idx
          break
    out_lines.append(line[:cut])
  txt = "\n".join(out_lines)
  # Quote unquoted object keys
  txt = re.sub(r'([{\[,]\s*)([A-Za-z_][A-Za-z0-9_]*)\s*:', r'\1"\2":', txt)
  # Remove trailing commas
  txt = re.sub(r",\s*([}\]])", r"\1", txt)
  return json.loads(txt)


def get_timestep_years() -> list[int]:
  cfg = load_json5("model_config.json5")
  t = cfg["timesteps"]
  years = []
  years.extend(range(t["bce_start"], t["bce_end"] + 1, t["bce_step"]))
  years.extend(range(t["ce_century_start"], t["ce_century_end"] + 1, t["ce_century_step"]))
  years.extend(range(t["ce_decade_start"], t["ce_decade_end"] + 1, t["ce_decade_step"]))
  years.extend(range(t["ce_annual_start"], t["ce_annual_end"] + 1, t["ce_annual_step"]))
  return years


def get_grid_coords(width: int = 4320, height: int = 2160) -> tuple[np.ndarray, np.ndarray]:
  lat_step = 180.0 / height
  lon_step = 360.0 / width
  lat = np.array([90.0 - (y + 0.5) * lat_step for y in range(height)], dtype=np.float32)
  lon = np.array([-180.0 + (x + 0.5) * lon_step for x in range(width)], dtype=np.float32)
  return lat, lon


def compute_cell_area_by_row(width: int = 4320, height: int = 2160) -> np.ndarray:
  a = 6378137.0
  f = 1.0 / 298.257223563
  b = a * (1.0 - f)
  dlat_rad = math.pi / height
  dlon_rad = (2.0 * math.pi) / width
  e2 = 1.0 - (b * b) / (a * a)
  e = math.sqrt(e2)
  factor = (b * b * dlon_rad) / (2.0 * 1e6)

  def compute_q(phi: float) -> float:
    s = math.sin(phi)
    es = e * s
    return (s / (1.0 - es * es)) + (math.log((1.0 + es) / (1.0 - es)) / (2.0 * e))

  areas = np.zeros(height, dtype=np.float64)
  for y in range(height):
    phi_north = (math.pi / 2.0) - y * dlat_rad
    phi_south = (math.pi / 2.0) - (y + 1) * dlat_rad
    areas[y] = abs(factor * (compute_q(phi_north) - compute_q(phi_south)))
  return areas


def encode_float32_geopng(file_path: str | Path, arr: np.ndarray) -> None:
  p = Path(file_path)
  p.parent.mkdir(parents=True, exist_ok=True)
  arr_f32 = np.asarray(arr, dtype=np.float32)
  rgba = arr_f32.astype(">f4", copy=False).view(np.uint8).reshape((*arr_f32.shape, 4))
  Image.fromarray(rgba, mode="RGBA").save(p, compress_level=6)


def decode_float32_geopng(file_path: str | Path) -> np.ndarray:
  im = Image.open(file_path)
  return np.frombuffer(im.tobytes(), dtype=">f4").reshape((im.size[1], im.size[0])).astype(np.float32)


def encode_uint32_geopng(file_path: str | Path, arr: np.ndarray) -> None:
  """Encodes integer array as big-endian uint32 RGBA GeoPNG (matching SVEA GeoPNG int32 format)."""
  p = Path(file_path)
  p.parent.mkdir(parents=True, exist_ok=True)
  arr_u32 = np.asarray(arr, dtype=np.uint32)
  rgba = arr_u32.astype(">u4", copy=False).view(np.uint8).reshape((*arr_u32.shape, 4))
  Image.fromarray(rgba, mode="RGBA").save(p, compress_level=6)


def decode_uint32_geopng(file_path: str | Path) -> np.ndarray:
  im = Image.open(file_path)
  return np.frombuffer(im.tobytes(), dtype=">u4").reshape((im.size[1], im.size[0])).astype(np.uint32)


def create_categorical_netcdf4(
  nc_path: str | Path,
  var_name: str,
  years: list[int],
  lat: np.ndarray,
  lon: np.ndarray,
  cell_area_km2: np.ndarray,
  flag_values: list[int],
  flag_meanings: str,
  global_attrs: dict,
  var_long_name: str,
  chunk_shape: tuple[int, int, int] = (1, 216, 432),
  complevel: int = 4,
):
  """Creates and returns an open CF-1.8 compressed NetCDF-4 Dataset ready for slice writing."""
  p = Path(nc_path)
  p.parent.mkdir(parents=True, exist_ok=True)
  if p.exists():
    p.unlink()
  ds = nc.Dataset(str(p), "w", format="NETCDF4")
  ds.Conventions = "CF-1.8"
  for k, v in global_attrs.items():
    setattr(ds, k, v)

  ds.createDimension("time", len(years))
  ds.createDimension("lat", len(lat))
  ds.createDimension("lon", len(lon))

  t_var = ds.createVariable("time", "i4", ("time",))
  t_var[:] = np.asarray(years, dtype=np.int32)
  t_var.standard_name = "time"
  t_var.long_name = "Signed calendar year (negative = BC, 0 = 1 AD matching HYDE, positive = AD)"
  t_var.units = "years"
  t_var.axis = "T"

  lat_var = ds.createVariable("lat", "f4", ("lat",))
  lat_var[:] = np.asarray(lat, dtype=np.float32)
  lat_var.standard_name = "latitude"
  lat_var.long_name = "Latitude of gridcell centre"
  lat_var.units = "degrees_north"
  lat_var.axis = "Y"

  lon_var = ds.createVariable("lon", "f4", ("lon",))
  lon_var[:] = np.asarray(lon, dtype=np.float32)
  lon_var.standard_name = "longitude"
  lon_var.long_name = "Longitude of gridcell centre"
  lon_var.units = "degrees_east"
  lon_var.axis = "X"

  area_var = ds.createVariable("cell_area_km2", "f4", ("lat",))
  area_var[:] = np.asarray(cell_area_km2, dtype=np.float32)
  area_var.long_name = "Exact WGS84 ellipsoidal area of each 5-arcminute gridcell by latitude row"
  area_var.units = "km2"

  data_var = ds.createVariable(
    var_name,
    "u1",
    ("time", "lat", "lon"),
    zlib=True,
    complevel=complevel,
    chunksizes=chunk_shape,
  )
  data_var.long_name = var_long_name
  data_var.coordinates = "time lat lon"
  data_var.flag_values = np.asarray(flag_values, dtype=np.uint8)
  data_var.flag_meanings = flag_meanings
  return ds, data_var


def render_categorical_spotcheck_sheet(
  out_png: str | Path,
  panels: list[dict],
  palette_rgb: dict[int, tuple[int, int, int]],
  legend_items: list[tuple[int, str]],
  sheet_title: str,
  ncols: int = 3,
) -> None:
  """Renders a crisp multi-panel categorical spot-check sheet with graticule and legend."""
  p = Path(out_png)
  p.parent.mkdir(parents=True, exist_ok=True)
  n = len(panels)
  nrows = math.ceil(n / ncols)
  fig = plt.figure(figsize=(ncols * 5.6, nrows * 3.8 + 1.4), dpi=140, facecolor="#14181D")
  fig.suptitle(sheet_title, color="#F0F4F8", fontsize=14, fontweight="bold", y=0.98)

  lut = np.zeros((256, 3), dtype=np.uint8)
  lut[:] = (162, 210, 238)
  for code, rgb in palette_rgb.items():
    if 0 <= code < 256:
      lut[code] = rgb

  for idx, pan in enumerate(panels):
    ax = fig.add_subplot(nrows, ncols, idx + 1)
    lat_min, lat_max = pan["lat_min"], pan["lat_max"]
    lon_min, lon_max = pan["lon_min"], pan["lon_max"]
    r0 = max(0, int(math.floor(((90.0 - lat_max) / 180.0) * 2160)))
    r1 = min(2160, int(math.ceil(((90.0 - lat_min) / 180.0) * 2160)))
    c0 = max(0, int(math.floor(((lon_min + 180.0) / 360.0) * 4320)))
    c1 = min(4320, int(math.ceil(((lon_max + 180.0) / 360.0) * 4320)))
    sub = pan["data"][r0:r1, c0:c1].astype(np.uint8)
    rgb = lut[sub]
    ax.imshow(
      rgb,
      extent=[lon_min, lon_max, lat_min, lat_max],
      interpolation="nearest",
      aspect="auto",
    )
    ax.set_title(pan["title"], color="#E2E8F0", fontsize=9.5, pad=4)
    ax.tick_params(colors="#94A3B8", labelsize=7)
    ax.grid(color="#000000", alpha=0.22, linewidth=0.5)
    for spine in ax.spines.values():
      spine.set_edgecolor("#475569")

  if legend_items:
    handles = [
      mpatches.Patch(
        facecolor=tuple(c / 255.0 for c in palette_rgb.get(code, (128, 128, 128))),
        edgecolor="#334155",
        label=f"{code}: {label}",
      )
      for code, label in legend_items
    ]
    fig.legend(
      handles=handles,
      loc="lower center",
      ncol=min(6, max(3, len(legend_items) // 3 + 1)),
      fontsize=7.2,
      frameon=False,
      labelcolor="#E2E8F0",
      bbox_to_anchor=(0.5, 0.01),
    )

  plt.tight_layout(rect=[0.01, 0.08 if legend_items else 0.02, 0.99, 0.95])
  fig.savefig(str(p), facecolor=fig.get_facecolor())
  plt.close(fig)
