# Changelog

All notable changes to the `ndvi2gif` package will be documented in this file.

---

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### ⚠️ Values change for some indices on Sentinel-2 and MODIS

**Sentinel-2 and MODIS bands are now surface reflectance (0-1), as Landsat's always were.** Both collections store reflectance as integers multiplied by 10000 and were used as such, so every index that is not a pure band ratio came out wrong on them. Ratios (NDVI, NDWI, MNDWI, NDMI, NBR, GNDVI, NDRE...) cancel the scale and **do not change**. Indices with a constant, a sum or an inverse **do**, and the old values were wrong — over the same Doñana summer composite:

| Index | Landsat | Sentinel-2 before → now | MODIS before → now |
|---|---|---|---|
| SAVI | 0.14 | 0.27 → 0.13 | 0.34 → 0.16 |
| EVI | 0.14 | 0.40 → 0.13 | 0.45 → 0.15 |
| LAI | 0.38 | 1.38 → 0.35 | 1.49 → 0.42 |
| AVI | 0.25 | 550 → 0.25 | 347 → 0.26 |
| WI2015 | 1.5 | −2255 → 1.5 | −2023 → 1.5 |
| CRI1 | 3.2 | 0.0 → 2.5 | 0.0 → 4.4 |

Also affected: `cri2`, `fai`, `mcari` and `ireci`, and the magnitude (not the sign) of `awei`/`aweinsh`. Ratios of differences such as `wdrvi`, `psri`, `mtci`, `reip` and `s2rep` are scale-free and do not change. Results computed with these indices on Sentinel-2 or MODIS with earlier versions should be recomputed. The raw bands (`index='red'`...) were already rescaled and do not change.

### ⚠️ Values change for Sentinel-1

**Sentinel-1 preprocessing and polarimetric indices now work on linear power.** `COPERNICUS/S1_GRD` is served in dB, and it was used as if it were linear:

- **Terrain correction and speckle filtering were applied to dB values.** The terrain correction multiplies backscatter by a factor between 0.5 and 2 — a ±3 dB adjustment on linear power, but on dB it could turn −13 dB into −26 dB on steep slopes. `S1ARDProcessor` gains `input_format` (`'DB'` by default, for `COPERNICUS/S1_GRD`; `'LINEAR'` for `S1_GRD_FLOAT`) and converts to linear before processing. `vv` and `vh` are still returned in dB; on flat terrain they barely change, on slopes they change a lot and are now correct.
- **Indices were computed on dB**, where ratios and sums have no physical meaning. They now convert to linear power first and follow their published definitions (over a Doñana scene, before → now):

  | Index | Formula (linear power) | Before → now |
  |---|---|---|
  | `rvi` | `4·VH / (VV + VH)` | 2.42 → 0.70 |
  | `vv_vh_ratio` | `VV / VH` | 0.65 → 4.9 |
  | `rfdi` | `(VV − VH) / (VV + VH)` (dual-pol adaptation of Mitchard et al. 2012) | −0.53 → 0.65 |
  | `dpsvi` | `(VV² + VV·VH) / √2`, the per-pixel DPSVIm of dos Santos et al. (2021) | was `(VV − VH)/(VV + VH)` on dB |

  Periasamy's original DPSVI needs the maximum VV of the whole scene, so it is not a per-pixel index; the modified form is used and documented as such.
- `vsdi` is unchanged but now documented as **experimental**: no publication defining it has been found (it was attributed to Periasamy 2018, who does not define it).
- `S1ARDProcessor.process_image()` returned an `ee.Element` instead of an `ee.Image` when called directly on an image.

### Added

- **🛰️ Multi-sensor classification in `LandCoverClassifier`**: pass a list of `NdviSeasonality` instances, one per sensor, and give the indices per sensor as a dict:

  ```python
  clf = LandCoverClassifier([s2, s1])
  clf.create_feature_stack(indices={'S2': ['ndvi', 'ndmi'], 'S1': ['vv', 'vh']})
  ```

  Band names carry the sensor as prefix (`S2_ndvi_2024_winter`, `S1_vh_2024_march`) only when there is more than one sensor, so single-sensor code and band names are unchanged. The typical use is combining optical indices with Sentinel-1 backscatter where clouds leave the optical series full of gaps.
- **Choice of resampling when sensors differ in resolution** (`resample=`): `'coarser'` aggregates the finer sensors by their mean onto the grid of the coarsest one, `'finer'` interpolates the coarser ones bilinearly onto the finest grid, and a number sets the pixel size in meters. There is deliberately no default when resolutions differ — the two options mean different things (losing detail vs. smoothing without adding information) — so the classifier raises a `ValueError` explaining them. The common grid is the UTM zone of the ROI unless `crs=` says otherwise. Verified pixel by pixel: a 30 m Sentinel-2 value on the Landsat grid is exactly the mean of the nine 10 m pixels inside it.

### Changed

- `LandCoverClassifier` samples training data, normalizes and exports at the pixel size of the feature stack (`clf.scale`) instead of fixed values: sampling was hardcoded at 10 m and normalization at 30 m whatever the sensor, and `export_results()` defaulted to 30 m.

### Fixed

- **`index='cig'` (Chlorophyll Index Green) was unreachable.** Its method was in the dispatch dictionary but the index was never registered for any sensor, so the constructor rejected it with `ValueError`. It only needs the green and NIR bands, so it is now available on every optical sensor (Sentinel-2, Landsat, MODIS, Sentinel-3), like `gndvi`. A new test fails if any index in the dispatch dictionary is left without a sensor again.
- **`aweinsh` had the sign of its SWIR2 term flipped**: it computed `4·(Green − SWIR1) − 0.25·NIR + 2.75·SWIR2` instead of Feyisa et al.'s (2014) `4·(Green − SWIR1) − (0.25·NIR + 2.75·SWIR2)`. It is one of the water indices of `HydroperiodAnalyzer`, whose water masks with `index='aweinsh'` change accordingly.
- **`wi2015` did not separate water from land**: its coefficients had been divided by 100 (on the belief that the published ones were for unscaled digital numbers) while the constant 1.7204 was kept, so the index sat around 1.5 everywhere and was always "water" at `HydroperiodAnalyzer`'s threshold of 0. It now uses Fisher et al.'s (2016) `1.7204 + 171·G + 3·R − 70·NIR − 45·SWIR1 − 71·SWIR2` on reflectance in 0-1 (at the Guadalquivir mouth: water +8.5, land −24). A new test checks that every water index of the analyzer is positive on water and negative on land.
- **`nmi` (NMDI) added the two SWIR bands instead of subtracting them**: Wang & Qu (2007) define `(NIR − (SWIR1 − SWIR2)) / (NIR + (SWIR1 − SWIR2))`.
- **Sentinel-3 offered indices it cannot compute.** Twelve indices that need SWIR bands (`mndwi`, `ndmi`, `awei`, `aweinsh`, `nbr`, `nbri`, `ndbi`, `ndsi`, `ndti`, `msi`, `nmi`, `wi2015`) were registered for Sentinel-3, whose OLCI sensor has no SWIR, and failed with an Earth Engine error when computed. They are no longer accepted for `sat='S3'`, which now raises a clear `ValueError` instead.
- **`floating_algae` and `tsi` failed on Sentinel-3**, the only sensor that offers them: they looked up a band called `'NIR'` while it is named `'Nir'`.
- **`lst` and `utfvi` were offered on sensors without thermal bands**: on Sentinel-2 and Sentinel-3 they silently returned fully masked images. `lst` is now available on Landsat and MODIS, and `utfvi` on Landsat only. A new test computes every index registered for each sensor on a real image of it.
- `cloud_filter` is now documented as what it is: `False` disables both the scene-level filter and the pixel-level cloud mask. To keep the pixel mask without discarding scenes, use `cloud_filter=True, max_cloud_cover=100`.
- **Default pixel size for MODIS, ERA5 and CHIRPS.** `_default_scale_for_sat()`, used by `export_to_drive()` / `export_to_asset()` when no `scale` is given, returned 250 m for MODIS although its reflectance comes from MOD09A1 at 500 m, and 30 m for ERA5-Land and CHIRPS, whose grids are about 11 km and 5.5 km: exports were four times too large for MODIS and tens of thousands of times for the climate datasets, with no added information. They now default to 500 m, 11132 m and 5566 m.
- **`LandCoverClassifier.get_accuracy_report()` always raised `KeyError`**: it read `producer_accuracy`/`user_accuracy` while the metrics are stored as `producers_accuracy`/`consumers_accuracy`, and treated them as dicts although Earth Engine returns them as a column and a row indexed by class value. It now returns one row per class present in the validation set plus the overall accuracy and kappa.
- **`LandCoverClassifier.get_feature_importance()` always raised `ValueError`**: it looked for `'RandomForest'` in the Python type of the classifier, which is `ee.Classifier` for every algorithm. It now works after `classify_supervised()` with `'random_forest'`, `'cart'` or `'gradient_tree'`, and returns a plain dict sorted by importance instead of a server-side object.
- **`LandCoverClassifier` ignored `class_property`** outside `add_training_data()`: training and accuracy assessment always used `'class'`, so labels stored under any other name made training fail.
- `LandCoverClassifier` temporal statistics selected an index's bands with the pattern `<index>_.*`, so `vv` also caught the `vv_vh_ratio` bands and mixed them into the `vv` mean, std, max and min. The pattern now requires the year after the index name.

### Added

- Tests for `HydroperiodAnalyzer`, which had none: hydrological-year bounds, index validation, and an Earth Engine integration test over the Doñana marshes checking that same-day tiles are mosaicked, that the midpoint weights tile the cycle from day 0 to 365, that `0 <= hydroperiod <= valid_days <= 365`, that `first_flood_doy <= last_flood_doy`, and that the downloadable mask stack is `uint8` with only the codes 0, 1, 2 and 255.

## [1.5.1] - 2026-09-18

### Fixed

- **`get_year_composite()` shifted the band names of every period after an empty one.** A period without a single image (scenes dropped by the cloud filter, a sensor not yet in orbit, acquisition gaps) was dropped, and the remaining bands were then named after the *first* N periods. With monthly Sentinel-2 composites over Doñana in 2017, where February and March have no scenes, the year came out with 10 bands: April's data was labelled `february`, May's `march`, and so on. Reducing a multi-year collection band by band then mixed different months without any error. Every year now has exactly one band per period, always under its own name: an empty period is a fully masked band (exported as nodata), or zeros for `key='count'`, since zero valid observations is a real value there. `get_period_composite()` returns the same placeholder, so the time-series and classification tools that call it get it too.
- **`LandCoverClassifier.create_feature_stack()` dropped the last year and mislabelled years after a skipped one.** It iterated `range(end_year - start_year)` while `end_year` is inclusive, so a 2018-2025 stack stopped at 2024. It also paired each composite with its year by position in the collection, but `get_year_composite()` skips the years without any image: with Sentinel-2 over 2014-2016, 2015's composite was labelled 2014 and 2016's 2015. Years are now matched by value, through the new `NdviSeasonality.period_scene_counts` (`{year: [scenes per period]}`, filled by `get_year_composite()`). Periods without images are left out of the stack with a message instead of entering as fully masked bands, since a single masked band masks every pixel on sampling and classification.
- **`SpatialTrendAnalyzer.calculate_pixel_trends()` only worked when the composite band was called `nd`.** It selected `'nd'` by name, which fails on Sentinel-1 (`'VH'`, `'RVI'`...) and with `key='percentile'` (`'nd_p90'`). The band is now taken by position.

### Changed

- `get_period_composite()` now always returns a float band (32-bit integer for `key='count'`), so that bands with and without data can be exported together.
- `get_year_composite()` checks data availability with one server call per year instead of one per period (10 calls instead of 120 for 10 years of monthly composites). Periods without images are reported by name.

## [1.5.0] - 2026-08-24

### Added

- **🎨 Raw reflectance bands as selectable indices**: the standardized bands can now be requested directly through `index=`, on Sentinel-2, Landsat and MODIS:
  - `index='blue'`, `'green'`, `'red'`, `'nir'`, `'swir1'`, `'swir2'`, plus `'red_edge1'`, `'red_edge2'` and `'red_edge3'` on Sentinel-2.
  - Values come back as **surface reflectance in 0-1 on every sensor**. Sentinel-2 and MODIS store reflectance as integers multiplied by 10000 and are rescaled on the fly (`self.reflectance_scale`); Landsat is already rescaled by `scale_OLI` / `scale_ETM`. The same numeric threshold therefore means the same thing whichever sensor produced it.
  - Not available on Sentinel-3, whose bands are top-of-atmosphere radiances with a different band set, nor on Sentinel-1, which already exposes `vv` and `vh`.

  Spectral indices are ratios and cancel out multiplicative changes in brightness, so a pixel can keep exactly the same NDVI while its reflectance drifts. Combined with the dispersion reducers added in 1.4.0, the bands map how radiometrically invariant each pixel is over a time series, which is how pseudo-invariant features (PIFs) are selected for relative radiometric normalization.

- **🔢 `key='count'`**: number of valid (non-masked) observations per pixel in each period. Not a value statistic but a quality layer: a standard deviation computed from three observations means very little, so `count` tells which parts of a dispersion map can be trusted. With `periods=1` it gives the yearly observation count per pixel, which also exposes the extra coverage in Sentinel-2 tile overlaps.

### Fixed

- **`SpatialTrendAnalyzer.calculate_pixel_trends()` dropped the last year of every trend map.** It iterated `range(start_year, end_year)` while `get_year_composite()` treats `end_year` as inclusive, so a 2018-2025 request was fitted on 2018-2024. The `magnitude` band was wrong on top of that: the slope was multiplied by `end_year - start_year` while the series only spanned one year less, overestimating the accumulated change by a full year.
- **`method='mann_kendall'` was documented but never implemented** — it fell through to the `else` branch and raised `ValueError`. It now returns Sen's slope, intercept and magnitude plus `tau`, Kendall's rank correlation with time, which is scale-free and so comparable across bands and sensors. `ee.Reducer.kendallsCorrelation()` also offers a `p-value` band, deliberately not returned: it comes back fully masked at every series length tested (n = 8 to 60, including trends `scipy.stats.kendalltau` scores below 1e-20), and a fully masked band silently masks whatever it is combined with. The `ValueError` for an unknown method now lists the valid ones.

### Changed

- `book/reference/indices.md`: the overview table now matches the code (110 variables, not 88) and documents the new bands. The per-section counts and a few listed-but-unimplemented indices further down that page are still out of sync and pending a separate revision.
- Dropped `pycrs` from the dependencies: it was declared but never imported anywhere in the package. Removed the unused `import fiona` from `ndvi2gif.py` as well, though `fiona` stays in the dependencies because geopandas uses it as a file engine.

## [1.4.0] - 2026-08-22

### Added

- **📊 Dispersion reducers**: four new `key` options in `NdviSeasonality` that describe how much an index **varies** inside each period, instead of its typical level:
  - `key='std'` — standard deviation of the observations in the period.
  - `key='variance'` — variance of the observations in the period.
  - `key='range'` — maximum minus minimum (within-period amplitude).
  - `key='cv'` — coefficient of variation (std / mean).

  They work with every sensor and index, and the resulting composites keep the usual period band names (`winter, spring, summer, autumn`, `january…december`, …), so they can be exported, animated and analysed exactly like the existing `max`/`median`/`mean` composites. Useful to map phenological change, disturbances and unstable surfaces such as flooded areas.

  `cv` divides by the mean, so it is only meaningful for indices that stay positive; a warning is printed when it is used with Sentinel-1 backscatter in dB.

- **💧 Water mask download in `HydroperiodAnalyzer`**: the per-date binary water masks can now be exported alongside the hydroperiod bands.
  - `get_water_masks_stack()` flattens the mask collection into a single `uint8` `ee.Image` with one band per acquisition date (`water_YYYY_MM_DD`), encoded as `0` = dry, `1` = water, `2` = observed but discarded as cloud/shadow, `255` = no scene covered the pixel or outside the ROI. Both no-value codes are configurable (`masked_value`, `nodata`).
  - `export_water_masks_to_drive()` and `export_water_masks_to_asset()` send that stack to Google Drive or to an Earth Engine asset.
  - `export_to_drive(..., include_masks=True)` and `export_to_asset(..., include_masks=True)` launch both exports at once and return a `(hydroperiod_task, masks_task)` tuple. The index, threshold and cycle are read from the `'index'`/`'threshold'`/`'hyd_year_start'` properties of the image being exported, so the masks always match it even if `compute_hydroperiod()` has since been called with other settings; the cached values are used only for images that carry no such metadata. The asset version accepts `masks_asset_id` and otherwise derives it from `asset_id` with a `_water_masks` suffix.
  - `get_water_masks()` gained an `add_footprint` flag adding a `'footprint'` band (1 inside the acquisition footprint of that date). This is what separates "cloudy" from "never observed" — the footprint survives the cloud mask, so the slanted edges of Landsat scenes and the gaps between orbits are identified as real no-data instead of being lumped in with clouds. Off by default, so `get_water_masks()` keeps returning a single-band collection.

  The masks go to a **separate** file on purpose: a GeoTIFF holds a single data type, so bundling them with the `int16` hydroperiod bands would promote them to `int16` and cancel out the saving. Earth Engine does not write a nodata tag into the GeoTIFF header, so `255` must be declared as nodata when reading the file.

### Fixed

- `tests/test_basic.py` no longer asserts the old silent fallback for invalid `sat`/`key` values (both raise `ValueError` since v1.2), and passes each sensor an index it actually supports.
- `tests/conftest.py` now initializes Earth Engine before collecting the tests. Building an `NdviSeasonality` already talks to the EE client, so the tests that make no EE calls of their own were failing anyway; a bare `ee.Initialize()` only works when the credentials carry a Cloud project, so `EARTHENGINE_PROJECT`/`GOOGLE_CLOUD_PROJECT` is used as a fallback. Without a session those tests are skipped with an explanatory message instead of erroring out.

## [1.3.1] - 2026-07-07

### Fixed

- Sensor status messages during collection setup are now printed only for the **selected** satellite, instead of for every supported sensor, reducing console verbosity (JOSS review, davemlz #4).

### No Breaking Changes

Full backward compatibility with v1.3.0.

## [1.3.0] - 2026-06-14

### Added

- **🌱 SpatialPhenologyAnalyzer**: New GEE-native module that produces **per-pixel phenology rasters** (Start/Peak/End of Season and derived metrics) for the whole ROI, entirely server-side. Complements the point-based phenology of `TimeSeriesAnalyzer`.
  - Output bands: `sos`, `pos`, `eos`, `los`, `amplitude`, `peak_value`, `baseline`, `growth_rate`, `senescence_rate` (SOS/POS/EOS/LOS in day-of-year).
  - Three server-side methods: `threshold` (amplitude crossing), `derivative` (steepest rate of change) and `harmonic` (per-pixel Fourier regression → smooth curve → threshold extraction). The `harmonic` method is the Earth Engine-native replacement for the client-side double-logistic fit, which relies on `scipy.optimize.curve_fit` and cannot run server-side.
  - Two outputs: `extract_phenology_rasters()` returns one image per year (`ee.ImageCollection`); `phenology_summary()` returns a single multi-year aggregate (`ee.Image`).
  - Export to local GeoTIFF (band names embedded via `rasterio`) or Google Drive (batch task, band names preserved by Earth Engine) through `export_target='local'|'drive'`.
- New documentation section "Per-pixel spatial phenology" in `book/advanced/time_series.md` and example notebook `examples_notebooks/Spatial_Phenology.ipynb`.

### Changed

- Added `rasterio` to the runtime dependencies (used to embed band names in locally exported GeoTIFFs).

### No Breaking Changes

Full backward compatibility with v1.2.x.

---

## [1.2.0] - 2026-04-12

### Added

- **💧 HydroperiodAnalyzer**: New GEE-native module for wetland and floodplain hydroperiod analysis.
  - Computes flood duration per pixel (days/year) entirely server-side using the midpoint temporal weighting method, based on the methodology of [phydroperiod](https://github.com/hectocore/phydroperiod).
  - Output bands: `hydroperiod`, `valid_days`, `normalized`, `first_flood_doy`, `last_flood_doy`.
  - Multi-year support: `compute_all_cycles()` and `compute_anomalies()`.
  - IRT (Irreplaceable Resource Taxonomy) metrics: global (`compute_irt_global()`) and per-pixel (`compute_irt_image()`).
  - Export to Google Drive and Earth Engine Assets.
- **SCL cloud masking for Sentinel-2**: New `scl_mask=True` parameter in `NdviSeasonality` and `mask_s2_scl()` method. Uses the Scene Classification Layer for more accurate cloud, shadow and cirrus detection. Set `scl_mask=False` to restore legacy QA60 behaviour.

### Fixed

- Removed `numpy<2.0` pin — numpy 2.x is now fully supported (`numpy>=1.24`).

### No Breaking Changes

Full backward compatibility with v1.1.0. `scl_mask=True` is the new default for Sentinel-2 cloud masking; set `scl_mask=False` to restore previous behaviour.

---

## [1.0.0] - 2025-12-28

### 🎉 **FIRST STABLE RELEASE - JOSS PUBLICATION**

This milestone release marks ndvi2gif as production-ready with comprehensive climate reanalysis data support (ERA5-Land and CHIRPS), expanding the library beyond vegetation monitoring into integrated climate-vegetation analysis. The library now supports 88 variables across 7 different satellite/reanalysis platforms, with intelligent handling of climate vs. vegetation data in time series analysis.

**NEW in 1.0.0**: Complete Jupyter Book documentation for JOSS (Journal of Open Source Software) submission, including comprehensive API reference, dataset guides, and usage tutorials.

---

## ✨ **New Features**

### 🌡️ **ERA5-Land Climate Reanalysis Support** (NEW)

- **NEW DATASET**: Complete integration with ECMWF ERA5-Land Daily Aggregated dataset
- **NEW**: 47 climate variables spanning 1950-present at ~11km resolution:
  - **Temperature (24 variables)**:
    - Basic: `temperature_2m`, `dewpoint_temperature_2m`, `skin_temperature`, `soil_temperature_level_1`
    - Daily min/max variants (8 variables with `_min` and `_max` suffixes)
    - Celsius conversions (12 variables with `_celsius` suffix)
  - **Precipitation & Water Balance (11 variables)**:
    - Meters: `total_precipitation_sum`, `total_evaporation_sum`, `potential_evaporation_sum`, `runoff_sum`, `surface_runoff_sum`
    - L/m² conversions (6 variables with `_lm2` suffix for intuitive units)
  - **Soil Moisture (4 variables)**: `volumetric_soil_water_layer_1` through `layer_4`
  - **Radiation (3 variables)**: Solar radiation and heat flux measurements
  - **Wind & Pressure (3 variables)**: Wind components and surface pressure
  - **Snow (2 variables)**: Snow depth and snowfall
- **NEW**: Unit conversion functions for user-friendly values:
  - Temperature: Kelvin to Celsius (K - 273.15)
  - Precipitation: Meters to L/m² (m × 1000)
- **NEW**: Support for daily aggregated statistics (mean, min, max, sum, median, percentile)

### 🌧️ **CHIRPS Precipitation Dataset** (NEW)

- **NEW DATASET**: Integration with UCSB Climate Hazards Group CHIRPS Daily precipitation
- **NEW**: High-resolution precipitation monitoring (1981-present, ~5.5km resolution)
- **NEW**: Global coverage from 50°S to 50°N latitude
- **NEW**: Combines satellite imagery with in-situ station data for improved accuracy
- **NEW**: Ideal for drought monitoring, precipitation climatology, and trend analysis

### 📊 **Enhanced Statistical Methods**

- **NEW**: `sum` statistic for temporal aggregation (essential for precipitation totals)
- **NEW**: `min` statistic for minimum value extraction (temperature minimums, etc.)
- **IMPROVED**: All statistical reducers now work with climate variables

### 🔧 **Time Series Analysis Improvements**

- **NEW**: Intelligent climate data detection in `TimeSeriesAnalyzer`
- **NEW**: Climate-specific summary statistics panel (replaces vegetation phenology for ERA5/CHIRPS)
- **NEW**: Seasonal climate statistics (winter, spring, summer, autumn averages)
- **NEW**: Annual mean and range calculations for climate variables
- **IMPROVED**: Dashboard automatically adapts display based on data type (vegetation vs. climate)

---

## 🐛 **Bug Fixes**

### 📅 **Critical: Inclusive Year Range**

- **FIXED**: `end_year` parameter is now **inclusive** instead of exclusive
  - **Before**: `start_year=2023, end_year=2023` would return NO data
  - **After**: `start_year=2023, end_year=2023` correctly includes all of 2023
- **FIXED**: Updated all year range calculations throughout the codebase
- **FIXED**: Corrected documentation to reflect inclusive behavior
- **IMPACT**: This is a **breaking change** for users who worked around the old exclusive behavior

### 🗺️ **Geometry Operations**

- **FIXED**: ROI centroid calculation now includes `maxError=1` parameter to prevent geometry operation errors
- **FIXED**: Improved error handling for centroid-based extractions in time series analysis

### 🛰️ **Sentinel-3 Band Naming**

- **FIXED**: Corrected uppercase/lowercase inconsistency in Sentinel-3 band names
- **FIXED**: Changed `'NIR'` to `'Nir'` for consistency with index calculations
- **IMPACT**: Sentinel-3 indices now work correctly without band selection errors

---

## 📚 **Documentation Updates**

### 📖 **README Enhancements**

- **ADDED**: Complete ERA5-Land variable documentation with units and descriptions
- **ADDED**: CHIRPS dataset documentation with coverage and use cases
- **ADDED**: Unit conversion examples (Celsius, L/m²)
- **UPDATED**: Supported datasets section with climate reanalysis platforms
- **UPDATED**: Project statistics reflecting 7 sensors and 88 total variables

### 🧪 **Testing**

- **ADDED**: Unit tests for ERA5 variable availability
- **ADDED**: Unit tests for CHIRPS precipitation integration
- **ADDED**: Tests for CHIRPS in satellite validation suite
- **UPDATED**: Satellite options tests to include ERA5 and CHIRPS

---

## 📊 **Project Statistics (v1.0.0)**

- **Supported Sensors**: 7 (S1, S2, S3, Landsat, MODIS, ERA5-Land, CHIRPS)
- **Total Variables**: 88
  - 40+ vegetation and environmental indices
  - 47 ERA5-Land climate variables
  - 1 CHIRPS precipitation variable
- **Temporal Coverage**: 1950-present (ERA5) and 1981-present (CHIRPS)
- **Spatial Resolutions**: 10m (S2) to ~11km (ERA5) to ~5.5km (CHIRPS)
- **ML Algorithms**: 8 (5 supervised, 3 unsupervised)

---

## 🔄 **API Changes**

### Breaking Changes

⚠️ **BREAKING**: `end_year` parameter behavior changed from exclusive to inclusive
```python
# Before v1.0.0 (exclusive)
NdviSeasonality(start_year=2020, end_year=2023)  # Processed 2020, 2021, 2022

# After v1.0.0 (inclusive)
NdviSeasonality(start_year=2020, end_year=2023)  # Processes 2020, 2021, 2022, 2023
```

### New Parameters

```python
# ERA5 climate data
processor = NdviSeasonality(
    sat='ERA5',
    index='temperature_2m_celsius',  # or any of 47 ERA5 variables
    key='mean',  # or 'min', 'max', 'sum', 'median', 'percentile'
    start_year=2020,
    end_year=2023  # Now inclusive!
)

# CHIRPS precipitation
chirps = NdviSeasonality(
    sat='CHIRPS',
    index='precipitation',
    key='sum',  # Monthly/seasonal totals
    start_year=2020,
    end_year=2023
)
```

---

## 🚀 **Example Use Cases**

### Climate Analysis with ERA5

```python
import ee
from ndvi2gif import NdviSeasonality
from ndvi2gif.timeseries import TimeSeriesAnalyzer

ee.Initialize()

# Temperature analysis
temp = NdviSeasonality(
    roi=ee.Geometry.Point([-6.48, 37.13]).buffer(5000),
    sat='ERA5',
    index='temperature_2m_celsius',
    periods=12,
    start_year=2020,
    end_year=2023,
    key='mean'
)

# Extract time series and analyze trends
analyzer = TimeSeriesAnalyzer(temp)
df = analyzer.extract_time_series()
trends = analyzer.analyze_trend(df=df)
fig = analyzer.plot_comprehensive_analysis()  # Shows climate stats, not phenology
```

### Precipitation Monitoring with CHIRPS

```python
# Monthly precipitation totals
precip = NdviSeasonality(
    roi=roi,
    sat='CHIRPS',
    index='precipitation',
    periods=12,
    start_year=2020,
    end_year=2023,
    key='sum'  # Sum for monthly totals
)

# Analyze drought patterns
analyzer = TimeSeriesAnalyzer(precip)
df = analyzer.extract_time_series()
trends = analyzer.analyze_trend(df=df, method='mann_kendall')
```

---

## 🙏 **Acknowledgments**

- **ERA5-Land**: ECMWF Climate Reanalysis ([dataset](https://developers.google.com/earth-engine/datasets/catalog/ECMWF_ERA5_LAND_DAILY_AGGR))
- **CHIRPS**: UCSB Climate Hazards Center ([Funk et al., 2015](https://doi.org/10.1038/sdata.2015.66))

---

## [0.6.0] - 2025-09-15

### 🧠 **MACHINE LEARNING & CLASSIFICATION RELEASE**

This release introduces comprehensive land cover classification capabilities and enhanced export functionality, positioning ndvi2gif as a complete remote sensing analysis suite. The library continues to mature toward v1.0.0 with advanced documentation and expanded analytical capabilities.

---

## ✨ **New Features**

### 🧠 **Land Cover Classification Module** (NEW)

- **NEW MODULE**: Complete `LandCoverClassifier` class for supervised and unsupervised classification
- **NEW**: Cloud masking options for Sentinel-2 and Landsat collections
- **NEW**: Multi-temporal feature stack generation with automatic normalization
- **NEW**: Support for multiple classification algorithms:
  - Random Forest (with feature importance)
  - Support Vector Machine (SVM)
  - Classification and Regression Trees (CART)
  - Naive Bayes
  - Gradient Tree Boost
- **NEW**: Unsupervised clustering algorithms:
  - K-means
  - Cascade K-means
  - Latent Dirichlet Allocation (LDA)
- **NEW**: Comprehensive accuracy assessment with confusion matrices
- **NEW**: Training data support from shapefiles, GeoJSON, and point/polygon sampling
- **NEW**: Feature importance analysis for Random Forest models
- **NEW**: Visualization tools for confusion matrices and accuracy reports

### 🚀 **Enhanced Export Capabilities**

- **NEW**: `export_to_drive()` - Batch export to Google Drive with full parameter control
- **NEW**: `export_to_asset()` - Export to Earth Engine Assets with pyramiding policies
- **NEW**: `_default_scale_for_sat()` - Automatic scale selection based on sensor
- **NEW**: Advanced export options including:
  - Custom pyramiding policies for classification data
  - Overwrite protection for assets
  - Format-specific options (compression, file per band)
  - Maximum pixel limits and CRS control

### 📊 **Time Series Analysis Enhancements** (Updated)

- **IMPROVED**: Enhanced documentation with complete examples
- **IMPROVED**: Better error handling and user feedback
- **IMPROVED**: More robust phenology extraction methods
- **IMPROVED**: Advanced visualization capabilities with publication-ready plots

---

## 🔧 **Major Improvements**

### 📚 **Documentation Overhaul**

- **IMPROVED**: Complete Sphinx-style docstrings for all classes and methods
- **IMPROVED**: Comprehensive parameter documentation with types and examples
- **IMPROVED**: Scientific references added to all spectral indices
- **IMPROVED**: Detailed usage examples in docstrings
- **IMPROVED**: Better error descriptions with suggested solutions
- **IMPROVED**: Cross-references between related methods

### 🛰️ **SAR Processing Enhancements**

- **IMPROVED**: Enhanced error handling in `S1ARDProcessor`
- **IMPROVED**: Better documentation for terrain correction parameters
- **IMPROVED**: More detailed method descriptions with scientific references
- **IMPROVED**: Improved parameter validation and user feedback

### 🌍 **API Consistency**

- **IMPROVED**: Consistent parameter naming across all modules
- **IMPROVED**: Standardized return types and error handling
- **IMPROVED**: Better integration between `NdviSeasonality` and new modules
- **IMPROVED**: More informative console output and progress tracking

---

## 🔄 **API Changes & Enhancements**

### 📦 **Module Structure**

```python
# NEW imports available in v0.6.0
from ndvi2gif import (
    NdviSeasonality,        # Core functionality (enhanced)
    S1ARDProcessor,         # SAR preprocessing (improved docs)
    TimeSeriesAnalyzer,     # Time series analysis (enhanced)
    SpatialTrendAnalyzer,   # Spatial analysis (enhanced)
    LandCoverClassifier,    # NEW: Classification workflows
)
```

### 🆕 **New Method Signatures**

```python
# NEW: Enhanced export methods
processor.export_to_drive(
    image=classified_map,
    description="landcover_2023",
    folder="ndvi2gif_results",
    scale=30,
    crs="EPSG:4326"
)

processor.export_to_asset(
    image=classification,
    asset_id="users/yourname/landcover_2023",
    pyramiding_policy={"class": "mode"},
    overwrite=True
)

# NEW: Classification workflow
classifier = LandCoverClassifier(processor)
features = classifier.create_feature_stack(
    indices=['ndvi', 'evi', 'ndwi'],
    include_statistics=True,
    normalize=True
)
classifier.add_training_data('training_points.shp')
result = classifier.classify_supervised('random_forest')
```

---

## 🛠️ **Under the Hood**

### 🔧 **Code Quality**

- **IMPROVED**: Consistent error handling with informative messages
- **IMPROVED**: Better type hints throughout the codebase
- **IMPROVED**: More robust parameter validation
- **IMPROVED**: Enhanced memory efficiency in large-area processing
- **IMPROVED**: Better handling of edge cases and invalid inputs

### 📈 **Performance**

- **OPTIMIZED**: Feature stack generation for classification
- **OPTIMIZED**: Memory usage in multi-temporal processing
- **OPTIMIZED**: Export operations with better chunking strategies

---

## 🐛 **Bug Fixes**

- **FIXED**: Improved error handling when no satellite data is available
- **FIXED**: Better validation of ROI inputs and coordinate systems
- **FIXED**: Enhanced handling of edge cases in temporal compositing
- **FIXED**: More robust processing of incomplete time series
- **FIXED**: Better handling of mixed sensor collections

---

## 📖 **Examples & Use Cases**

### 🌾 **Agricultural Monitoring**

```python
# Multi-temporal crop classification
processor = NdviSeasonality(
    roi='farm_boundaries.shp',
    sat='S2', periods=12,
    start_year=2022, end_year=2024
)

classifier = LandCoverClassifier(processor)
features = classifier.create_feature_stack(['ndvi', 'evi', 'ndre'])
classifier.add_training_data('crop_samples.shp')
crop_map = classifier.classify_supervised('random_forest')
```

### 🌊 **Water Quality Assessment**

```python
# Sentinel-3 water quality with export to Drive
processor = NdviSeasonality(
    roi='lake_boundary.shp',
    sat='S3', index='turbidity',
    periods=24  # Bi-monthly
)

composites = processor.get_year_composite()
processor.export_to_drive(
    image=composites.first(),
    description="lake_turbidity_2024",
    folder="water_quality"
)
```

### 🏔️ **SAR Forest Monitoring**

```python
# Advanced SAR processing with classification
processor = NdviSeasonality(
    sat='S1', index='rvi',
    use_sar_ard=True,
    sar_speckle_filter='REFINED_LEE',
    sar_terrain_correction=True
)

classifier = LandCoverClassifier(processor)
forest_map = classifier.classify_unsupervised('kmeans', n_clusters=5)
```

---

## ⚠️ **No Breaking Changes**

Full backward compatibility maintained with v0.5.x. All existing code continues to work unchanged.

---

## 🔄 **Migration Guide**

### From v0.5.x to v0.6.0

No breaking changes! All existing code will continue to work. New features are additive:

```python
# v0.5.x code continues to work unchanged
processor = NdviSeasonality(sat='S2', index='ndvi')
processor.get_gif('animation.gif')

# v0.6.0 adds new capabilities
classifier = LandCoverClassifier(processor)  # NEW
processor.export_to_drive(image, "export")   # NEW
```

---

## 🎯 **Future Roadmap**

### v1.0.0 (Planned) - Complete Climate Analysis Platform
- **📚 Jupyter Book**: Interactive documentation with comprehensive tutorials and examples
- **🌡️ Climate Datasets**: Integration with ERA5, CHIRPS, TerraClimate, and other climate model datasets
- **🌍 Climate Analysis**: Advanced tools for climate change impact assessment and adaptation planning

---

## 📊 **Statistics**

- **New classes**: 1 (`LandCoverClassifier`)
- **New methods**: 15+ (classification, enhanced exports, utilities)
- **Enhanced methods**: 20+ (improved documentation and error handling)
- **Lines of code**: ~3,500 → ~4,800 (+37%)
- **Documentation coverage**: 95%+ (comprehensive docstrings)

---

## 🙏 **Acknowledgments**

Special thanks to the Google Earth Engine team and the open-source remote sensing community for their continued support and feedback that made this release possible.

---

## 📚 **Documentation**

Complete documentation with tutorials available at: [GitHub Repository](https://github.com/Digdgeo/Ndvi2Gif)

---

**Full Changelog**: https://github.com/Digdgeo/Ndvi2Gif/compare/v0.5.0...v0.6.0

## [0.5.0] - 2025-08-28

### Added

- **🛰️ Sentinel-1 ARD Processor**: New `S1ARDProcessor` module for advanced SAR preprocessing:
  - Radiometric terrain correction (angular method, Vollrath et al. 2020).
  - Configurable speckle filters: Boxcar, Lee, Refined Lee, Gamma-MAP, Lee Sigma.
  - Flexible DEM options (Copernicus 30/90, SRTM 30/90).
- **📈 TimeSeriesAnalyzer**: New module for time series and phenological analysis:
  - Robust extraction of temporal profiles from points or polygons.
  - Trend analysis (Mann-Kendall, Linear regression, Sen’s slope).
  - Comprehensive dashboards (trend, seasonality, autocorrelation, quality).
  - Phenological metrics (SOS, EOS, POS, LOS, amplitude, growth/senescence rates).
- **🌱 NdviSeasonality improvements**:
  - Extended ROI handling: DEIMS sites, Sentinel-2 MGRS tiles, Landsat WRS path/row, shapefiles, GeoJSON.
  - Flexible temporal periods (4, 12, 24, or custom definitions).
  - Optional SAR normalization and enhanced orbit handling.
  - More robust sensor-index validation.

### Changed

- **Visualization**: Unified plotting style with Seaborn/Matplotlib, clearer layouts.
- **Documentation**: Updated examples covering SAR and time series analysis.

### Fixed
- More robust handling of null/NaN values in temporal extraction.
- Minor bug fixes in period generation and export routines.

## [0.4.1] - 2025-07-21

### Added

- Just fixing some bugs in Readme.md

## [0.4.0] - 2025-07-21

### Added

- **🛰️ Sentinel-3 OLCI Support**: Revolutionary addition with 21 spectral bands and daily global coverage
- **🌊 Advanced Water Quality Indices**: 10 specialized aquatic indices including OCI, TSI, CDOM, turbidity, SPM, KD490, floating algae detection
- **🔬 Enhanced Sentinel-2**: Complete Red Edge implementation with Surface Reflectance for superior data quality  
- **💧 Cyanobacteria Detection**: New NDCI index for harmful algal bloom monitoring and water quality assessment
- **⚙️ SAR Orbit Control**: Precise control over Sentinel-1 ascending/descending orbits for geometric consistency
- **🎯 40+ Specialized Indices**: Comprehensive coverage with intelligent sensor-index validation
- **📊 Professional Architecture**: Clean, extensible design with enhanced error handling and documentation

### New Sentinel-3 Indices

- **OCI** - OLCI Chlorophyll Index (optimized for S3)
- **TSI** - Trophic State Index (water quality assessment) 
- **CDOM** - Colored Dissolved Organic Matter Index
- **Turbidity** - Water Turbidity Index (sediment monitoring)
- **SPM** - Suspended Particulate Matter Index
- **KD490** - Diffuse Attenuation Coefficient at 490nm
- **Floating Algae** - Floating Algae Index (bloom detection)
- **Red Edge Position** - OLCI-optimized red edge position
- **Fluorescence Height** - Chlorophyll fluorescence detection
- **Water Leaving Reflectance** - Aquatic reflectance analysis

### New SAR Indices

- **RFDI** - Radar Forest Degradation Index (deforestation monitoring)
- **VSDI** - Vegetation Scattering Diversity Index (structural diversity)

### Enhanced Features

- **Intelligent Validation**: Smart index-sensor compatibility checking prevents invalid combinations
- **Orbit Parameter**: Fine control over Sentinel-1 orbit selection (BOTH/ASCENDING/DESCENDING)
- **Advanced Use Cases**: Support for pseudo-invariant area radiometric normalization workflows

### Changed

- **Sentinel-2 to Surface Reflectance**: Upgraded from TOA to Surface Reflectance for better scientific quality
- **Simplified Architecture**: Removed unnecessary complexity while maintaining full functionality
- **Enhanced Documentation**: Professional-grade docstrings and examples

### Technical Improvements

- Modular sensor setup with clean separation of concerns
- Comprehensive sensor-index mapping and validation
- Enhanced error messages for better user experience
- Support for advanced radiometric normalization workflows

---

## [0.3.0] - 2025-07-17

### Added

- **New SAR Indices**: RVI (Radar Vegetation Index), VV/VH ratio, VH, VV, DPSVI for Sentinel-1
- **Flexible Percentiles**: Support for any percentile value (1-99) instead of fixed 90/95
- **Enhanced Sentinel-1**: VV+VH dual polarization with speckle filtering
- **Robust ROI Handling**: Support for drawn features, lists of features, and improved geometry conversion
- **Incomplete Year Support**: Automatic detection and processing of available periods for current/incomplete years
- **Enhanced Dependencies**: Added pycrs and deims as core dependencies (now available in conda)
- **Example Notebooks**: Comprehensive examples in `examples_notebooks/` folder

### Fixed

- ROI conversion for drawn geometries and feature lists from geemap
- Speckle filter now preserves temporal properties (system:time_start)
- Band naming consistency for SAR indices
- Error handling for missing data periods
- Dependency issues with pycrs and deims

### Improved

- More robust error handling throughout the library
- Better documentation and examples
- Enhanced support for agricultural monitoring workflows
- Simplified installation process

## [0.2.0] - 2025-01-27

### Added

- **Dynamic period generation**: Support for any number of temporal periods (4, 6, 8, 12, 24, 52, or any custom number).
- **Flexible temporal analysis**: Easy configuration from traditional 4 seasons to 52 weekly periods or any custom division.
- **Enhanced extensibility**: Adding new satellites and datasets is now trivial with the unified architecture.

### Changed

- **Major code refactoring**: Eliminated over 90% of code duplication by replacing 40+ individual period functions with a single dynamic system.
- **Improved maintainability**: Reduced codebase from ~3,000 lines to ~400 lines while maintaining all functionality.
- **Enhanced performance**: Streamlined period generation and composite creation.

### Technical Details

- Replaced hardcoded period definitions with dynamic `_generate_periods()` method.
- Consolidated all `get_winter()`, `get_january()`, `get_p1()` through `get_p24()` functions into a single `get_period_composite()` method.
- Maintained full backward compatibility - all existing code works without changes.
- Added comprehensive leap year handling to prevent date-related errors.

### Breaking Changes

- None - this release maintains 100% backward compatibility.

---

## [0.1.5] - 2025-05-26

### Fixed

- Fixing bug with MNDWI index.

---

## [0.1.4] - 2025-05-25

### Fixed

- Nothing really changes, just a f* problem with release version management.

---

## [0.1.3] - 2025-05-25

### Fixed

- Nothing really changes, just a f* problem with release version management.

---

## [0.1.2] - 2025-05-25

### Added

- Complete rework and translation of the README into Markdown format.
- Included new seasonal/statistical methods and updated docstrings in English.
- Added support for region input via Sentinel-2 tiles and Landsat path/row.
- Added rich ROI input documentation with tabular summary.
- Added `deims` dependency as optional to avoid conda forge problems

### Changed

- Clarified the purpose of the library as a broader seasonal analysis tool, not just for GIF generation.
- Cleaned and validated `setup.cfg` and `pyproject.toml`.
- Added extra requirements group for `deims`.

---

## [0.1.1] - 2025-05-21

### Fixed

- Fixed rendering issue in `README.rst` that caused PyPI upload failure.
- Rebuilt and republished the package with correct long description format.

---

## [0.1.0] - 2025-05-21

### Added

- Compatibility with Conda packaging and `conda-forge` ecosystem.
- Included `MANIFEST.in` to ensure `LICENSE` and `README.rst` are bundled in source distribution.
- Improved `README.rst` formatting to comply with PyPI rendering rules.

### Changed

- Switched versioning to semantic 0.x.y style for future compatibility.
- Cleaned and validated metadata to allow upload to both PyPI and Conda Forge.

### Note

- This is a technical release — no changes to the core functionality.

---

## [0.0.9] - 2025-05-20

### Changed

- Version bump to align `setup.cfg`, PyPI and GitHub release.
- No functional changes from version 0.0.7.

---

## [0.0.7] - 2025-05-20

### Added

- New method `get_ndmi()` to compute the Normalized Difference Moisture Index (NDMI).
- New (old) method `get_gif()` to download a gif for the selected index/bands.
- Package structure modernized:
  - Added `setup.cfg` and `pyproject.toml` (PEP 517/518 compliant).
  - Optional removal of legacy `setup.py`.
- Updated dependencies:
  - `geemap` pinned to version `0.29.5`.
  - `numpy` constrained to `<2.0` for compatibility.

### Fixed

- Compatibility issues with recent versions of `geemap`, `xarray`, and `numpy`.
- Resolved import error caused by the removal of `np.unicode_` in NumPy 2.0.

---

## [0.0.6] - 2023-03-10

### Added

- Initial public release of the `ndvi2gif` package.
- Generate seasonal composites and extract statistical summaries from several remote sensing index using Google Earth Engine and geemap.
- Export to animated GIF and GeoTIFF format.
