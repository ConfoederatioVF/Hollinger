# Anthromes (`anthromes/`)

Global 5-arcminute (`4320 x 2160`, WGS84) Anthromes 12K anthropogenic biome classification across 218 timesteps (`10000 BC` to `2025 AD`, signed years `-10000..2025`, `0 = 1 AD`), modelled on the Anthromes 12K DGG v1 scheme (Ellis et al., 2020, 2021) and coupled with Hollinger's dynamic `landarea/landarea.nc` and Holocene biomes from `climate/holdridge_lifezones.nc`.

## Output Files

- **`anthromes/anthromes.nc`**: CF-1.8 compressed NetCDF-4 (`zlib=True`, `complevel=4`, chunks `[1, 216, 432]`), variable `anthromes` (`uint8`) with dimensions `(time, lat, lon)` (`218 x 2160 x 4320`), `cell_area_km2`, `flag_values`, and `flag_meanings`.
- **`anthromes/3.rasters/anthromes_<year>.png`**: 218 SVEA `uint32` (`>u4`) big-endian RGBA GeoPNGs.

## Classification Scheme

All thresholds, class codes, hierarchy levels, and RGB palettes are defined in `common/anthromes_config.json5`.

| Code | Identifier | Name | Level | Type | RGB |
| ---: | :--- | :--- | :--- | :--- | :--- |
| `0` | `Water` | Water / Ocean / Lake | Water | Water | `(162, 210, 238)` |
| `11` | `Urban` | Urban | Dense settlements | Intensive | `(168, 0, 0)` |
| `12` | `Mixed_settlements` | Mixed settlements | Dense settlements | Intensive | `(255, 0, 0)` |
| `21` | `Rice_villages` | Rice villages | Villages | Intensive | `(0, 112, 255)` |
| `22` | `Irrigated_villages` | Irrigated villages | Villages | Intensive | `(0, 169, 230)` |
| `23` | `Rainfed_villages` | Rainfed villages | Villages | Intensive | `(169, 0, 230)` |
| `24` | `Pastoral_villages` | Pastoral villages | Villages | Intensive | `(255, 115, 223)` |
| `31` | `Residential_irrigated_croplands` | Residential irrigated croplands | Croplands | Intensive | `(0, 255, 197)` |
| `32` | `Residential_rainfed_croplands` | Residential rainfed croplands | Croplands | Intensive | `(230, 230, 0)` |
| `33` | `Populated_croplands` | Populated croplands | Croplands | Intensive | `(255, 255, 115)` |
| `34` | `Remote_croplands` | Remote croplands | Croplands | Intensive | `(255, 255, 190)` |
| `41` | `Residential_rangelands` | Residential rangelands | Rangelands | Intensive | `(230, 152, 0)` |
| `42` | `Populated_rangelands` | Populated rangelands | Rangelands | Intensive | `(255, 211, 127)` |
| `43` | `Remote_rangelands` | Remote rangelands | Rangelands | Intensive | `(255, 235, 175)` |
| `51` | `Residential_woodlands` | Residential woodlands | Cultured | Cultured | `(56, 168, 0)` |
| `52` | `Populated_woodlands` | Populated woodlands | Cultured | Cultured | `(165, 245, 122)` |
| `53` | `Remote_woodlands` | Remote woodlands | Cultured | Cultured | `(211, 255, 178)` |
| `54` | `Inhabited_drylands` | Inhabited drylands | Cultured | Cultured | `(217, 189, 117)` |
| `61` | `Wild_woodlands` | Wild woodlands | Wildlands | Wild | `(218, 242, 234)` |
| `62` | `Wild_drylands` | Wild drylands | Wildlands | Wild | `(225, 225, 225)` |
| `63` | `Ice` | Ice | Wildlands | Wild | `(250, 255, 255)` |

## Decision Tree (Strict Top-to-Bottom Precedence)

1. **Water / Glaciated mask (`landarea.nc`)**:
   - `la == 0.0` -> `0` (`Water`)
   - `la < 0.0` -> `63` (`Ice`)
2. **Dense settlements** (`la > 0.0`):
   - `urban >= 0.2` or `popd >= 2500.0` -> `11` (`Urban`)
   - Else if `popd >= 100.0`:
     - If `not pot_vill` -> `12` (`Mixed_settlements`)
     - Else (`pot_vill` is `True`):
       - `rice >= 0.2` -> `21` (`Rice_villages`)
       - Else if `irrigation >= 0.2` -> `22` (`Irrigated_villages`)
       - Else if `crops >= 0.2` -> `23` (`Rainfed_villages`)
       - Else if `grazing >= 0.2` -> `24` (`Pastoral_villages`)
       - Else -> `12` (`Mixed_settlements`)
3. **Croplands**:
   - `crops >= 0.2`:
     - `popd >= 10.0` and `irrigation >= 0.2` -> `31` (`Residential_irrigated_croplands`)
     - `popd >= 10.0` and `irrigation < 0.2` -> `32` (`Residential_rainfed_croplands`)
     - `1.0 <= popd < 10.0` -> `33` (`Populated_croplands`)
     - `popd < 1.0` -> `34` (`Remote_croplands`)
4. **Rangelands**:
   - `grazing >= 0.2`:
     - `popd >= 10.0` -> `41` (`Residential_rangelands`)
     - `1.0 <= popd < 10.0` -> `42` (`Populated_rangelands`)
     - `popd < 1.0` -> `43` (`Remote_rangelands`)
5. **Mixed intensive use (`used >= 0.2`)**:
   - `crops > grazing` -> `34` (`Remote_croplands`)
   - Else -> `43` (`Remote_rangelands`)
6. **Cultured (`popd > 0.0`)**:
   - If `trees`: `51` (`popd >= 10.0`), `52` (`1.0 <= popd < 10.0`), or `53` (`popd < 1.0`)
   - Else -> `54` (`Inhabited_drylands`)
7. **Wildlands (`popd == 0.0`, `used < 0.2`)**:
   - If `ice_barren` -> `63` (`Ice`)
   - Else if `trees` -> `61` (`Wild_woodlands`)
   - Else -> `62` (`Wild_drylands`)

## Dynamic Coupling & Normalisation

- **Coastal-sliver density normalisation & marginal/island fusion**: On `la > 0.0` where HYDE land exists (`la_hyde > 0.0`), land-use fractions and population density (`popd`, persons km⁻²) are normalised by `denom = max(la, la_hyde)` (where `la` is dynamic `landarea.nc[t]` and `la_hyde` is HYDE `land_area_per_cell.asc`). On marginal, reclaimed, or newly recovered coastal cells (`la > 0.0` and `la_hyde == 0.0`), intensive HYDE land-use/population densities are extrapolated from the nearest valid HYDE land cell (`1.0` within 1.5 px, `0.65` within 3 px where GHSL-populated or pre-1975, and smooth exponential decay onto contiguous pre-1975 Holocene paleoshelves). Small island nations and populated coasts omitted from HYDE 3.5 are fused with GHSL (`GHS_POP_E1975..E2025` + `GHS_urban_1975..2025`) and backcast prior to 1975 using Gapminder country series and historical settlement onset dates. Fractions are clamped to `[0, 1]`; if `crops + grazing + urban > 1.0`, they are rescaled proportionally, with `rice <= crops`, `irrigation <= crops`, and `used = min(1.0, crops + grazing + urban)`.
- **Dynamic Holocene biomes**: `trees_base_modern` comes from `potveg15.tif` (`1..8`), anchored against `holdridge_lifezones.nc` at `1985` (`trees_hlz_modern`, Holdridge classes `8-10, 14-16, 20-23, 27-30, 34-38`). On valid `potveg15` cells, `trees(t)` retains `trees_base_modern` unless `trees_hlz(t) != trees_hlz_modern`, in which case it shifts to `trees_hlz(t)`. On emergent paleoshelves where `potveg15` is NODATA, `trees(t) = trees_hlz(t)`. `ice_barren` combines glaciated `la < 0.0` with polar/ice Holdridge cells (`hlz == 1`) on `potveg15 == 15` or emergent shelves.
- **Extended `pot_vill` mask**: `potvill20.tif` is extended onto emergent paleoshelves via Euclidean nearest-neighbour (`scipy.ndimage.distance_transform_edt`), keeping North America, Australia, New Zealand, Greenland, and Antarctica shelves `False` whilst Sundaland, Doggerland, the Persian Gulf, the Yellow Sea, and New Guinea remain `True`.
- **HYDE 3.5 temporal alignment**: Exact values are used for all `128` native HYDE 3.5 timesteps (`-10000..2025`), and linear interpolation across bounding HYDE millennia is applied for the `90` intermediate BCE centuries.

## Pipeline Usage

```bash
python anthromes/pipeline/run_pipeline.py            # Run full pipeline (A -> B -> C)
python anthromes/pipeline/run_pipeline.py --only=B,C # Run specific steps
python anthromes/pipeline/run_pipeline.py --skip=A   # Skip source download
```

## Citations

- Ellis, E.C., Beusen, A.H.W., & Klein Goldewijk, K. (2020). Anthropogenic biomes: 10,000 BCE to 2015 CE. *Land*, 9(5), 129.
- Ellis, E.C., Gauthier, N., Klein Goldewijk, K., et al. (2021). People have shaped most of terrestrial nature for at least 12,000 years. *PNAS*, 118(17), e2023483118.
- Klein Goldewijk, K., et al. HYDE 3.5 historical population and land-use database.
