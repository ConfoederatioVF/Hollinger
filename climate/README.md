# Climate Indicator (`climate/`)

High-resolution (5-arcminute, `4320 × 2160`, WGS84 `EPSG:4326`) paleoclimate and historical climate classification datasets spanning **10000 BC – 2025 AD** across the canonical **218** Hollinger timesteps (`-10000` to `2025`, where `0 = 1 AD`).

---

## 1. Output Artefacts

- `climate/koeppen_geiger.nc` — CF-1.8 compressed NetCDF-4 (`uint8`, codes `0..30`, chunks `[1, 216, 432]`, `zlib=True`, `complevel=4`).
- `climate/trewartha.nc` — CF-1.8 compressed NetCDF-4 (`uint8`, codes `0..17`, chunks `[1, 216, 432]`, `zlib=True`, `complevel=4`).
- `climate/holdridge_lifezones.nc` — CF-1.8 compressed NetCDF-4 (`uint8`, codes `0..38`, chunks `[1, 216, 432]`, `zlib=True`, `complevel=4`).
- `climate/3.rasters/` — `654` (`218 × 3`) SVEA `uint32` (`>u4`) big-endian RGBA GeoPNGs (`koeppen_geiger_<year>.png`, `trewartha_<year>.png`, `holdridge_lifezones_<year>.png`).

All outputs mask water cells strictly where `land_area == 0.0` in `landarea/landarea.nc` (both ice-free land `> 0` and glaciated land / permanent ice shelves `< 0` are classified).

---

## 2. Data Sources & Downscaling Methodology

### 2.1 Primary Sources
1. **CCSM3 TraCE-21ka (`d651050`)**: Monthly `TREFHT` (2 m air temperature, K), `PRECC` + `PRECL` (convective + large-scale precipitation rate, m s⁻¹), and `PHIS` (surface geopotential, m² s⁻²) across segments `trace.16` (`12000–11701 BP`, `-10050 CE`) through `trace.36` (`400 BP–1990 CE`).
2. **WorldClim 2.1 (5-arcminute)**: Modern (`1970–2000`) monthly mean temperature (`wc2.1_5m_tavg`), monthly precipitation (`wc2.1_5m_prec`), and elevation (`wc2.1_5m_elev`).
3. **CRU TS 4.10 (`1901–2025`)**: Monthly `0.5° × 0.5°` land temperature (`tmp`, °C) and precipitation (`pre`, mm month⁻¹).
4. **GLAC-1D (`TOPicemsk.GLACD26k...`) & `landarea/landarea.nc`**: Paleo ice-sheet / GIA topography (`HDC`, m a.s.l.) and dynamic land–sea–ice masks.

### 2.2 Climatology Extraction & Anomaly Downscaling
1. **CCSM3 Centred Windows**: Monthly climatologies are computed using a 50-year centred window (`±25` years) for `year ≤ 1700` and a 30-year centred window (`±15` years) for `1710 ≤ year ≤ 1990`, plus the `1960–1990` reference normal. Longitudes are wrapped from `0..360` to `-180..180` and latitudes flipped to north-to-south (`+89..−89`).
2. **CRU TS 4.10 & Seamless Splice**: Centred 30-year monthly climatologies are computed for `1900..2025` alongside `1970–2000` and `1960–1990` normals. Ocean and island gaps are filled smoothly via longitude-periodic distance-transform nearest-neighbour and Laplace relaxation. Over Antarctica (`lat < −60°`), anomalies fall back to CCSM3. The `1960–1990` vs `1970–2000` CRU offset anchors CCSM3 to the WorldClim baseline, and a linear blend across `1900–1950` ensures zero step-jump at the CCSM3–CRU transition.
3. **Gap-Free 5-Arcminute WorldClim Reference**: Emergent paleoshelves (Doggerland, Sundaland, Sahul, Beringia, Persian Gulf, Yellow Sea) and coastal/Antarctic gaps are inpainted without NaNs using sea-level temperature reduction (`Γ = 0.0065 °C m⁻¹`) + Laplace relaxation and precipitation ratio inpainting against coarse CCSM3 precipitation.
4. **Paleo Elevation Lapse Adjustment**: For `year < 1950`, where GLAC-1D elevation anomalies (`HDC(t) − HDC(0 ka)`) differ from coarse CCSM3 `PHIS / 9.80665` anomalies, a residual lapse correction (`−0.0065 × ΔZ_res`) is applied to 5-arcminute monthly temperatures.

---

## 3. Classification Schemes & Integer Codes

### 3.1 Köppen–Geiger (`koeppen_geiger.nc`, codes `0..30`)
Following Peel et al. (2007) and Beck et al. (2018, 2023). High-sun (summer) is April–September for `lat ≥ 0°` and October–March for `lat < 0°`. Aridity threshold `P_thresh` (cm) is `2×MAT` if `≥70%` of annual precipitation `MAP` falls in winter, `2×MAT + 28` if `≥70%` falls in summer, and `2×MAT + 14` otherwise. Precedence: **E → B → A → C/D**.

| ID | Code | Description | Rule Summary |
|---:|:---|:---|:---|
| 0 | `Water` | Water / Ocean / Lake | `land_area == 0.0` |
| 1 | `Af` | Tropical, rainforest | `T_cold ≥ 18`, `P_dry ≥ 60` |
| 2 | `Am` | Tropical, monsoon | `T_cold ≥ 18`, `P_dry < 60`, `P_dry ≥ 100 − MAP/25` |
| 3 | `Aw` | Tropical, savannah | `T_cold ≥ 18`, `P_dry < 100 − MAP/25` |
| 4 | `BWh` | Arid, desert, hot | `MAP < 5×P_thresh×10`, `MAT ≥ 18` |
| 5 | `BWk` | Arid, desert, cold | `MAP < 5×P_thresh×10`, `MAT < 18` |
| 6 | `BSh` | Arid, steppe, hot | `5×P_thresh×10 ≤ MAP < 10×P_thresh×10`, `MAT ≥ 18` |
| 7 | `BSk` | Arid, steppe, cold | `5×P_thresh×10 ≤ MAP < 10×P_thresh×10`, `MAT < 18` |
| 8–10 | `Csa`, `Csb`, `Csc` | Temperate, dry summer | `0 < T_cold < 18`, `P_sdry < 40 & < P_wwet/3`; `a: T_hot≥22`, `b: T_mon10≥4`, `c: else` |
| 11–13 | `Cwa`, `Cwb`, `Cwc` | Temperate, dry winter | `0 < T_cold < 18`, not `s`, `P_wdry < P_swet/10`; `a/b/c` as above |
| 14–16 | `Cfa`, `Cfb`, `Cfc` | Temperate, no dry season | `0 < T_cold < 18`, neither `s` nor `w`; `a/b/c` as above |
| 17–20 | `Dsa`, `Dsb`, `Dsc`, `Dsd` | Cold, dry summer | `T_cold ≤ 0`, `T_hot ≥ 10`, dry summer; `d: T_cold < −38` |
| 21–24 | `Dwa`, `Dwb`, `Dwc`, `Dwd` | Cold, dry winter | `T_cold ≤ 0`, `T_hot ≥ 10`, dry winter; `d: T_cold < −38` |
| 25–28 | `Dfa`, `Dfb`, `Dfc`, `Dfd` | Cold, no dry season | `T_cold ≤ 0`, `T_hot ≥ 10`, no dry season; `d: T_cold < −38` |
| 29 | `ET` | Polar, tundra | `0 ≤ T_hot < 10` |
| 30 | `EF` | Polar, frost | `T_hot < 0` |

### 3.2 Trewartha (`trewartha.nc`, codes `0..17`)
Following Trewartha & Horn (1980), Patton (1962), and Belda et al. (2014). Uses `T_mon10` (months with `T_m ≥ 10 °C`) and Patton's aridity index `R = 2.3×MAT − 0.64×P_w% + 41` (cm), where `P_w%` is the percentage of annual precipitation falling in the winter six months. Precedence: **F → B → A → C/D/E**.

| ID | Code | Description | Rule Summary |
|---:|:---|:---|:---|
| 0 | `Water` | Water / Ocean / Lake | `land_area == 0.0` |
| 1 | `Ar` | Tropical wet | `T_cold ≥ 18`, `≤ 2` months with `P_m < 60 mm` |
| 2 | `Am` | Tropical monsoon | `T_cold ≥ 18`, `> 2` dry months, `P_dry ≥ 100 − MAP/25` |
| 3 | `Aw` | Tropical wet-dry (dry winter) | `T_cold ≥ 18`, winter drier than summer |
| 4 | `As` | Tropical wet-dry (dry summer) | `T_cold ≥ 18`, summer drier than winter |
| 5 | `BWh` | Subtropical/tropical desert | `MAP < 5×R`, `T_mon10 ≥ 8` |
| 6 | `BWk` | Temperate/boreal desert | `MAP < 5×R`, `T_mon10 < 8` |
| 7 | `BSh` | Subtropical/tropical steppe | `5×R ≤ MAP < 10×R`, `T_mon10 ≥ 8` |
| 8 | `BSk` | Temperate/boreal steppe | `5×R ≤ MAP < 10×R`, `T_mon10 < 8` |
| 9 | `Cs` | Subtropical dry summer | `T_mon10 ≥ 8`, `P_sdry < 30 & < P_wwet/3`, `P_ssum < P_wsum`, `MAP < 890` |
| 10 | `Cw` | Subtropical dry winter | `T_mon10 ≥ 8`, not `Cs`, `P_wdry < P_swet/10`, `P_wsum < P_ssum` |
| 11 | `Cr` | Subtropical humid | `T_mon10 ≥ 8`, neither `Cs` nor `Cw` |
| 12 | `Do` | Temperate oceanic | `4 ≤ T_mon10 ≤ 7`, `T_cold > 0` |
| 13 | `Dc` | Temperate continental | `4 ≤ T_mon10 ≤ 7`, `T_cold ≤ 0` |
| 14 | `Eo` | Boreal oceanic | `1 ≤ T_mon10 ≤ 3`, `T_cold > 0` |
| 15 | `Ec` | Boreal continental | `1 ≤ T_mon10 ≤ 3`, `T_cold ≤ 0` |
| 16 | `Ft` | Polar tundra | `0 ≤ T_hot < 10` |
| 17 | `Fi` | Polar ice cap | `T_hot < 0` |

### 3.3 Holdridge Life Zones (`holdridge_lifezones.nc`, codes `0..38`)
Following Holdridge (1947, 1967), Leemans (1990), and Lugo et al. (1999). Biotemperature is `T_bio = mean(clip(T_m, 0, 30))`. Cells with `T_bio < 1.5 °C` are assigned `1` (`Polar_Desert`). Remaining land cells are assigned via Euclidean nearest-neighbour in `(log2 T_bio, log2 P_ann)` space to the canonical hexagonal centroids defined in `common/climate_config.json5`, with the frost line `T_cold ≥ 6.0 °C` separating Subtropical (`24..30`) from Warm Temperate (`17..23`) zones.

- `0`: `Water`
- `1`: `Polar_Desert`
- `2–5`: `Subpolar_Dry_Tundra`, `Subpolar_Moist_Tundra`, `Subpolar_Wet_Tundra`, `Subpolar_Rain_Tundra`
- `6–10`: `Boreal_Desert`, `Boreal_Dry_Scrub`, `Boreal_Moist_Forest`, `Boreal_Wet_Forest`, `Boreal_Rain_Forest`
- `11–16`: `CoolTemp_Desert`, `CoolTemp_Desert_Scrub`, `CoolTemp_Steppe`, `CoolTemp_Moist_Forest`, `CoolTemp_Wet_Forest`, `CoolTemp_Rain_Forest`
- `17–23`: `WarmTemp_Desert`, `WarmTemp_Desert_Scrub`, `WarmTemp_Thorn_Steppe`, `WarmTemp_Dry_Forest`, `WarmTemp_Moist_Forest`, `WarmTemp_Wet_Forest`, `WarmTemp_Rain_Forest`
- `24–30`: `Subtrop_Desert`, `Subtrop_Desert_Scrub`, `Subtrop_Thorn_Woodland`, `Subtrop_Dry_Forest`, `Subtrop_Moist_Forest`, `Subtrop_Wet_Forest`, `Subtrop_Rain_Forest`
- `31–38`: `Trop_Desert`, `Trop_Desert_Scrub`, `Trop_Thorn_Woodland`, `Trop_Very_Dry_Forest`, `Trop_Dry_Forest`, `Trop_Moist_Forest`, `Trop_Wet_Forest`, `Trop_Rain_Forest`

---

## 4. Pipeline Execution

```bash
python climate/pipeline/run_pipeline.py            # Full run (Steps A, B, C)
python climate/pipeline/run_pipeline.py --skip=A   # Re-run climatology & classification
python climate/pipeline/run_pipeline.py --only=C   # Re-classify only
```

---

## 5. Citations

- Beck, H. E., Zimmermann, N. E., McVicar, T. R., Vergopolan, N., Berg, A., & Wood, E. F. (2018). Present and future Köppen-Geiger climate classification maps at 1-km resolution. *Scientific Data*, 5, 180214.
- Beck, H. E., et al. (2023). High-resolution (1 km) Köppen-Geiger maps for 1901–2099 based on constrained CMIP6 projections. *Scientific Data*, 10, 724.
- Belda, M., Holtanová, E., Halenka, T., & Kalvová, J. (2014). Climate classification revisited: from Köppen to Trewartha. *Climate Research*, 59(1), 1–13.
- Fick, S. E., & Hijmans, R. J. (2017). WorldClim 2: new 1-km spatial resolution climate surfaces for global land areas. *International Journal of Climatology*, 37(12), 4302–4315.
- Harris, I., Osborn, T. J., Jones, P., & Lister, D. (2020). Version 4 of the CRU TS monthly high-resolution gridded multivariate climate dataset. *Scientific Data*, 7, 109.
- Holdridge, L. R. (1947). Determination of world plant formations from simple climatic data. *Science*, 105(2727), 367–368.
- Holdridge, L. R. (1967). *Life Zone Ecology*. Tropical Science Center, San José.
- Leemans, R. (1990). *Possible changes in natural vegetation patterns due to a global warming*. IIASA Working Paper WP-90-08.
- Liu, Z., et al. (2009). Transient simulation of last deglaciation with a new mechanism for Bølling-Allerød warming. *Science*, 325(5938), 310–314. (CCSM3 TraCE-21ka).
- Lugo, A. E., Brown, S. L., Dodson, R., Smith, T. S., & Shugart, H. H. (1999). The Holdridge life zones of the conterminous United States in relation to ecosystem mapping. *Journal of Biogeography*, 26(5), 1025–1038.
- Patton, C. P. (1962). A note on the classification of dry climates in the Köppen system. *California Geographer*, 3, 105–112.
- Peel, M. C., Finlayson, B. L., & McMahon, T. A. (2007). Updated world map of the Köppen-Geiger climate classification. *Hydrology and Earth System Sciences*, 11(5), 1633–1644.
- Trewartha, G. T., & Horn, L. H. (1980). *An Introduction to Climate* (5th ed.). McGraw-Hill.
