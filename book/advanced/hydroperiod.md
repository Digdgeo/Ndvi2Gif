# Hydroperiod Analysis with `HydroperiodAnalyzer`

`HydroperiodAnalyzer` computes **hydroperiod** — the number of days a pixel remains flooded during a hydrological cycle — entirely in Google Earth Engine. No data downloads are required: all computation runs server-side as a lazy GEE graph that is only evaluated when you add a layer to a map or trigger an export.

---

## What is hydroperiod?

Hydroperiod is a fundamental ecological indicator for wetland ecosystems. It quantifies how long each pixel is covered by water during a defined hydrological cycle (by default September 1st – August 31st). It drives:

- Plant community composition and distribution
- Habitat suitability for waterbirds, amphibians, and invertebrates
- Carbon dynamics and greenhouse gas fluxes
- Long-term wetland monitoring and conservation planning

---

## Methodology: midpoint temporal weighting

Raw satellite archives contain an uneven number of cloud-free observations over time. Simply counting flood detections ignores this irregularity — a pixel observed only in winter would appear drier than one observed year-round, even if both were equally flooded.

The **midpoint weighting method** solves this by distributing the 365 days of the cycle proportionally among available scenes:

1. Sort all cloud-free scenes by acquisition date within the cycle.
2. For each consecutive pair of scenes, compute the temporal midpoint.
3. Each scene is assigned a **territorial weight** spanning from the midpoint with its predecessor to the midpoint with its successor. The first and last scenes use day 0 and day 365 as outer boundaries.
4. A pixel contributes `weight` flood days if it is water in that scene, and `weight` valid days regardless (as long as it is not cloud-masked).

```
Scenes (days from Sep 1):   0    14    45    120   230   310
Midpoints:                      7    29    82    175   270
Territories:             [0,7] [7,29] [29,82] [82,175] [175,270] [270,365]
Weights (days):           7     22     53      93       95        95
```

The sum of all weights is always exactly 365. The final hydroperiod for a pixel is the sum of weights across scenes where that pixel was classified as water.

### Same-day tile mosaicking

Sentinel-2 and other sensors can produce several overlapping scenes from different tiles on the same calendar day. Without correction, same-day scenes get a midpoint distance of zero → weight of 0 days, silently discarding valid observations.

`HydroperiodAnalyzer` solves this by mosaicking all scenes acquired on the same day with `.max()` before computing weights: if *any* tile shows water for a given pixel on that day, the pixel is classified as water.

---

## Installation and setup

```python
import ee
import geemap
from ndvi2gif import NdviSeasonality, HydroperiodAnalyzer

# Authenticate only the first time
# ee.Authenticate()
ee.Initialize(project='your-project-id')
```

`HydroperiodAnalyzer` requires an existing `NdviSeasonality` instance. All satellite configuration — ROI, sensor, date range, cloud filtering, band standardisation — is reused from it.

---

## Quick start

```python
roi = ee.Geometry.Rectangle([-6.55, 36.85, -6.10, 37.20])  # Doñana, Spain

ns = NdviSeasonality(
    roi=roi,
    sat='S2',
    start_year=2022,
    end_year=2023,
    index='mndwi',
    cloud_filter=True,
    max_cloud_cover=20,
    scl_mask=True,       # pixel-level SCL cloud/shadow mask (recommended for S2)
)

ha = HydroperiodAnalyzer(ns)

result = ha.compute_hydroperiod(index='mndwi', threshold=0.0, hyd_year=2022)
# result is an ee.Image with bands: hydroperiod, valid_days, normalized,
#                                   first_flood_doy, last_flood_doy
```

---

## Step-by-step guide

### 1. Configure `NdviSeasonality`

`NdviSeasonality` is the entry point for all satellite data in ndvi2gif. For hydroperiod analysis, the key parameters are:

| Parameter | Recommendation |
|---|---|
| `sat` | `'S2'` (10 m), `'Landsat'` (30 m), `'MODIS'` (500 m) |
| `start_year` / `end_year` | Define the range of hydrological cycles to analyse |
| `cloud_filter=True` | Filters entire scenes above a cloud cover threshold |
| `max_cloud_cover` | Scene-level threshold (%). Lower = stricter. 10–20% is typical |
| `scl_mask=True` | **Recommended for S2.** Pixel-level SCL cloud/shadow mask; more accurate than QA60 |

```python
ns = NdviSeasonality(
    roi=roi,
    sat='S2',
    start_year=2019,
    end_year=2023,
    index='mndwi',
    cloud_filter=True,
    max_cloud_cover=15,   # stricter than default 20
    scl_mask=True,
)
```

> **Note on `scl_mask`:** The SCL (Scene Classification Layer) is a Sentinel-2 specific band that classifies each pixel at acquisition time. Setting `scl_mask=True` (default) masks classes 1 (saturated), 3 (cloud shadows), 8 (medium clouds), 9 (high clouds), and 10 (thin cirrus). This is more accurate than the QA60 band used by the legacy `scl_mask=False` mode, especially for cloud shadows over water. It benefits all ndvi2gif modules, not just hydroperiod.

---

### 2. Create `HydroperiodAnalyzer`

```python
ha = HydroperiodAnalyzer(ns)
# HydroperiodAnalyzer(sat='S2', hyd_year=2019/2020, hyd_start=09-01)
```

By default the hydrological year starts on **September 1st**. For tropical or monsoon regions with a different rainy season, pass a custom start date:

```python
# April 1st start — suitable for tropical/monsoon regions
ha = HydroperiodAnalyzer(ns, hydrological_year_start=(4, 1))
```

---

### 3. Inspect water masks (optional)

`get_water_masks()` is called internally by `compute_hydroperiod()`, but you can call it manually to inspect the binary water masks and their temporal weights before running the full computation.

```python
masks = ha.get_water_masks(index='mndwi', threshold=0.0, hyd_year=2022)

n_scenes = masks.size().getInfo()
total_weight = sum(masks.aggregate_array('weight').getInfo())

print(f'Scenes in cycle: {n_scenes}')
print(f'Sum of weights: {total_weight:.1f} days')  # should be ~365

# Inspect scene weights
weights = masks.aggregate_array('weight').getInfo()
starts  = masks.aggregate_array('start_doy').getInfo()
ends    = masks.aggregate_array('end_doy').getInfo()

for w, s, e in zip(weights[:5], starts[:5], ends[:5]):
    print(f'  weight={w:.1f} d  [{s:.0f} → {e:.0f}]')
```

Each image in the returned collection has a `water` band (1 = water, 0 = dry, masked = cloud/nodata) and three properties: `weight`, `start_doy`, and `end_doy` (all in days from the hydrological year start).

---

### 4. Compute a single hydrological cycle

```python
result = ha.compute_hydroperiod(
    index='mndwi',
    threshold=0.0,
    hyd_year=2022,           # Sep 2022 – Aug 2023
    normalize=True,
    compute_first_last=True,
    min_flood_days=3,        # mask pixels with < 3 flood days (noise filter)
    permanent_threshold=0.95 # pixels always flooded → first=0, last=365
)
```

#### Output bands

| Band | Type | Description |
|---|---|---|
| `hydroperiod` | int16 | Accumulated weighted flood days. Reflects actual temporal coverage of the satellite archive |
| `valid_days` | int16 | Accumulated weighted days with valid (cloud-free) observations |
| `normalized` | int16 | `(hydroperiod / valid_days) × 365` — flood days rescaled to a full year, correcting for uneven coverage |
| `first_flood_doy` | int16 | Day of the hydrological year (0-indexed from Sep 1) when flooding was first detected |
| `last_flood_doy` | int16 | Day of the hydrological year when flooding was last detected |

> **`hydroperiod` vs `normalized`:** Use `normalized` for ecological interpretation. It corrects for the fact that some pixels may have fewer valid observations than others. A pixel observed only in the wet season would inflate its raw `hydroperiod`; `normalized` accounts for this by dividing by `valid_days` before scaling to 365.

#### `first_flood_doy` and `last_flood_doy` details

These bands are computed from the midpoint weighting structure:

- `first_flood_doy` = `start_doy` of the earliest scene where the pixel was detected as water. This represents the left boundary of the temporal territory of the first flood scene — i.e. "the pixel has been flooded *at least* since this day".
- `last_flood_doy` = `end_doy` of the latest scene where water was detected — "the pixel was still flooded *at least* until this day".
- For a pixel flooded in exactly one scene, `first_flood_doy` and `last_flood_doy` are the `start_doy` and `end_doy` of that scene — exactly the midpoints with the adjacent scenes.

Two filters are applied:

- **`min_flood_days`** (default 3): pixels with fewer accumulated flood days are masked. This removes isolated false detections from residual cloud shadows or atmospheric effects that passed through the cloud filter.
- **`permanent_threshold`** (default 0.95): pixels where the flood fraction (`hydroperiod / valid_days`) ≥ 0.95 are considered permanently inundated. They receive `first_flood_doy = 0` and `last_flood_doy = 365` because onset/recession dates are not ecologically meaningful for permanent water bodies (rivers, deep lagoons).

---

### 5. Visualise on a map

```python
Map = geemap.Map()
Map.centerObject(roi, zoom=10)

hyd_vis = {
    'bands': ['normalized'],
    'min': 0, 'max': 365,
    'palette': ['ffffff', 'cce5f5', '91bfdb', '4393c3', '2166ac', '053061'],
}
Map.addLayer(result, hyd_vis, 'Hydroperiod normalised')
Map.add_colorbar(hyd_vis, label='Flood days', orientation='horizontal')
Map
```

Palette interpretation: **white** = never flooded, **dark blue** = flooded nearly all year.

---

### 6. Water indices

Five water indices are supported. All use the `index > threshold` rule (default `threshold=0.0`):

| Index | Full name | Best for |
|---|---|---|
| `mndwi` | Modified NDWI (Xu 2006) | Open water, marshes — **recommended default** |
| `ndwi` | NDWI (McFeeters 1996) | Clear open water |
| `awei` | AWEI shade (Feyisa 2014) | Areas with mixed land cover and shadow |
| `aweinsh` | AWEI no-shadow (Feyisa 2014) | Areas without significant shadow |
| `wi2015` | Water Index 2015 (Fisher 2016) | General-purpose |

```python
# Compare indices for the same cycle
r_mndwi = ha.compute_hydroperiod(index='mndwi', threshold=0.0, hyd_year=2022)
r_awei  = ha.compute_hydroperiod(index='awei',  threshold=0.0, hyd_year=2022)

# Stricter threshold — maps only deep/permanent water
r_strict = ha.compute_hydroperiod(index='mndwi', threshold=0.3, hyd_year=2022)
```

> **Always pass `image=` explicitly when exporting after comparing indices.** Calling `compute_hydroperiod()` multiple times updates `ha._results` to the last call. Using `export_to_drive()` without an explicit `image` argument will export the *last* result, which may not be the one you intended. A `UserWarning` is raised in this situation.

---

### 7. Multiple hydrological cycles

`compute_all_cycles()` iterates over every year in the `NdviSeasonality` date range and returns a `{year: ee.Image}` dictionary. GEE builds the computation graphs for all cycles instantly — the actual computation happens on render or export.

```python
ns_multi = NdviSeasonality(roi=roi, sat='S2', start_year=2019, end_year=2023,
                           index='mndwi', cloud_filter=True)
ha_multi = HydroperiodAnalyzer(ns_multi)

cycles = ha_multi.compute_all_cycles(index='mndwi', threshold=0.0)
# cycles = {2019: ee.Image, 2020: ee.Image, 2021: ee.Image, 2022: ee.Image, 2023: ee.Image}
```

Visualise all cycles by toggling layers:

```python
Map = geemap.Map()
for yr, img in cycles.items():
    Map.addLayer(img.select('normalized'), hyd_vis, f'Hydroperiod {yr}/{yr+1}')
Map
```

Inter-annual difference between two specific cycles:

```python
diff = cycles[2022].select('normalized').subtract(cycles[2019].select('normalized'))

diff_vis = {
    'min': -150, 'max': 150,
    'palette': ['d73027', 'fc8d59', 'fee090', 'ffffff', 'e0f3f8', '91bfdb', '4575b4'],
}
Map.addLayer(diff, diff_vis, 'Δ 2022-23 vs 2019-20 (days)')
```

---

### 8. Mean hydroperiod and anomalies

`compute_anomalies()` computes a **mean hydroperiod** image and the **anomaly** of each cycle relative to that mean:

```
anomaly = cycle_normalized − mean_normalized
```

Positive values indicate a wetter-than-average year; negative values indicate a drier year.

#### Reference modes

The `reference` parameter controls how the mean is defined:

**`reference='period'` (default)**  
The mean is computed from the cycles in the current analysis window (`ns.start_year` → `ns.end_year`). Use this when you want to characterise inter-annual variability within a specific period.

```python
result_period = ha_multi.compute_anomalies(cycles=cycles, reference='period')

mean_hyd   = result_period['mean']       # ee.Image — band: mean_hydroperiod
anomalies  = result_period['anomalies']  # {year: ee.Image} — band: anomaly
```

**`reference='historical'`**  
The mean is computed from the full satellite archive, from the sensor's first reliable year to the last complete hydrological cycle. This gives a **climatological baseline** independent of the analysis window — useful when the analysis period is short or potentially biased (e.g. a run of wet or dry years).

```python
result_hist = ha_multi.compute_anomalies(cycles=cycles, reference='historical')

# Check the baseline metadata
mean_hist = result_hist['mean']
print(mean_hist.get('reference').getInfo())      # 'historical'
print(mean_hist.get('ref_start_year').getInfo()) # e.g. 2017 (S2)
print(mean_hist.get('ref_end_year').getInfo())   # last complete cycle
print(mean_hist.get('ref_n_cycles').getInfo())   # total cycles in baseline
```

Historical archive starts by sensor:

| Sensor | Start year | Notes |
|---|---|---|
| `S2` | 2017 | S2_SR_HARMONIZED quality data |
| `Landsat` | 1984 | Landsat 5 TM |
| `MODIS` | 2000 | MOD09A1 Terra |
| `S1` | 2014 | Sentinel-1A launch |

> **Performance note:** `reference='historical'` builds computation graphs for all archive cycles lazily. The Python call completes in milliseconds, but rendering on a map or exporting will be slower because GEE must evaluate more imagery. Already-computed cycles (from `compute_all_cycles`) are reused automatically to avoid redundant computation nodes.

---

### 9. Temporal Representativity Index (IRT)

The IRT measures how evenly the cloud-free observations are distributed across the hydrological year:

- **IRT = 1**: observations are perfectly spread across all months — hydroperiod estimates are most reliable.
- **IRT = 0**: all observations cluster in a single period — estimates may be biased.

A low IRT warns that hydroperiod values for that cycle should be interpreted with caution.

#### Global IRT (one value per cycle)

Computed client-side from scene acquisition dates using the Gini-based formula from the *phydroperiod* library:

```python
for yr in range(2019, 2024):
    irt = ha_multi.compute_irt_global(hyd_year=yr)
    print(f'{yr}/{yr+1}: IRT = {irt:.3f}')
```

#### Per-pixel IRT (`ee.Image`)

Estimates temporal coverage quality at each pixel using Simpson's diversity index — fully computable server-side in GEE:

```
IRT_pixel = N² / (n_periods × Σ nᵢ²)
```

where *N* is the total number of valid observations for the pixel and *nᵢ* is the count in monthly period *i*.

```python
irt_img = ha_multi.compute_irt_image(hyd_year=2022)

irt_vis = {
    'bands': ['irt'],
    'min': 0, 'max': 1,
    'palette': ['d7191c', 'fdae61', 'ffffbf', 'abd9e9', '2c7bb6'],
}
Map.addLayer(irt_img, irt_vis, 'Per-pixel IRT 2022-23')
```

> **On tile boundaries in the IRT image:** With Sentinel-2, the IRT image typically shows visible tile footprints because different tiles have different temporal coverage. Pixels at the edge of a tile may have fewer valid observations than pixels covered by two overlapping tiles. This is expected behaviour — it reflects the real data availability pattern.

---

### 10. Long time series with Landsat

Landsat provides data since 1984, enabling multi-decadal hydroperiod analysis. The interface is identical to Sentinel-2; just change `sat='Landsat'`:

```python
ns_ls = NdviSeasonality(
    roi=roi,
    sat='Landsat',
    start_year=1990,
    end_year=2020,
    index='mndwi',
    cloud_filter=True,
)
ha_ls = HydroperiodAnalyzer(ns_ls)

# Note: scl_mask has no effect for Landsat (SCL is Sentinel-2 specific)
ls_cycles = ha_ls.compute_all_cycles(index='mndwi', threshold=0.0)
```

Use `reference='historical'` with Landsat to leverage the full 1984–present archive as a climatological baseline.

---

### 11. Custom hydrological year start

The default September 1st start suits Mediterranean and temperate climates. For other regions, pass `hydrological_year_start` as a `(month, day)` tuple:

```python
# Tropical/monsoon regions: April 1st start
ha_tropical = HydroperiodAnalyzer(ns, hydrological_year_start=(4, 1))
# Cycle: Apr 2022 – Mar 2023

# Northern hemisphere boreal wetlands: May 1st start
ha_boreal = HydroperiodAnalyzer(ns, hydrological_year_start=(5, 1))
```

---

### 12. Exporting results

Both methods accept any `ee.Image` via the `image` parameter. Always pass the image explicitly if you have called `compute_hydroperiod()` more than once on the same `HydroperiodAnalyzer` instance.

#### To Google Drive

```python
# Single cycle
task = ha.export_to_drive(
    image=result,                        # pass explicitly
    folder='my_wetland',
    description='hydroperiod_2022_2023',
    scale=10,                            # 10 m — S2 native resolution
    crs='EPSG:25829',                    # local UTM projection
)

# All cycles
for yr, img in cycles.items():
    ha_multi.export_to_drive(
        image=img,
        folder='my_wetland',
        description=f'hydroperiod_{yr}_{yr+1}',
        scale=10,
        crs='EPSG:25829',
    )

# Mean and anomalies
ha_multi.export_to_drive(image=mean_hyd, folder='my_wetland',
                         description='hydroperiod_mean_2019_2023', scale=10)
for yr, anom in anomalies.items():
    ha_multi.export_to_drive(image=anom, folder='my_wetland',
                             description=f'hydroperiod_anomaly_{yr}_{yr+1}', scale=10)
```

Monitor tasks at [code.earthengine.google.com/tasks](https://code.earthengine.google.com/tasks).

#### To an Earth Engine Asset

```python
task = ha.export_to_asset(
    image=result,
    asset_id='projects/your-project/assets/hydroperiod_2022_2023',
    scale=10,
    crs='EPSG:25829',
)
```

---

## Tips and caveats

### Cloud filtering strategy

With Sentinel-2 and a spatially compact ROI (e.g. a single wetland), there is usually enough data to be strict. A `max_cloud_cover` of 10–15% reduces cloud contamination artefacts in `first_flood_doy` and `last_flood_doy` without sacrificing temporal coverage. The `min_flood_days` filter (default 3) provides an additional safeguard by masking sporadic false detections.

For large ROIs or sensors with lower revisit frequency (Landsat, MODIS), a higher `max_cloud_cover` threshold is often necessary to maintain temporal coverage.

### Choosing an index and threshold

`mndwi` with `threshold=0.0` is the best default for most inland wetlands. For turbid or shallow water (estuaries, reservoirs with high sediment load), MNDWI can yield false negatives (water not detected) because sediment increases NIR reflectance. In those situations, try `awei` or lower the threshold slightly (e.g. `threshold=-0.1`).

### ROI design

For coastal wetlands, exclude the open sea and estuarine zones from the ROI if they are not part of the study area. Turbid coastal waters are genuinely difficult for optical water indices at threshold=0.0, and including them adds noise to statistics.

### GEE lazy evaluation

All `compute_*` methods return `ee.Image` objects immediately — no computation happens until you call `.getInfo()`, add the layer to a map, or start an export. This means:

- Calling `compute_all_cycles()` on a 30-year Landsat range is essentially instant in Python.
- The computation cost is paid at render/export time by GEE's servers.
- `reference='historical'` in `compute_anomalies()` adds more computation nodes to the GEE graph but does not block Python execution.

---

## API reference

| Method | Returns | Description |
|---|---|---|
| `get_water_masks(index, threshold, hyd_year)` | `ee.ImageCollection` | Binary water masks with temporal weight properties |
| `compute_hydroperiod(...)` | `ee.Image` | Full hydroperiod computation for one cycle |
| `compute_all_cycles(...)` | `dict[int, ee.Image]` | Hydroperiod for every cycle in the date range |
| `compute_anomalies(cycles, reference)` | `dict` | Mean hydroperiod + per-cycle anomalies |
| `compute_irt_global(hyd_year)` | `float` | Global IRT (one value, client-side) |
| `compute_irt_image(hyd_year, n_periods)` | `ee.Image` | Per-pixel IRT image |
| `export_to_drive(image, folder, ...)` | `ee.batch.Task` | Export to Google Drive |
| `export_to_asset(asset_id, image, ...)` | `ee.batch.Task` | Export to Earth Engine Asset |

**Supported water indices:** `ndwi`, `mndwi`, `awei`, `aweinsh`, `wi2015`

**Supported sensors:** `S2` (10 m), `Landsat` (30 m), `MODIS` (500 m), `S1`*

*SAR-based water mapping with Sentinel-1 uses backscatter thresholding rather than optical indices.

---

## References

### Related Python library

- **`phydroperiod`** (García Díaz & Bustamante Díaz, 2026) — standalone Python library for hydroperiod computation from water masks. Available on PyPI: `pip install phydroperiod`. Repository: [https://pypi.org/project/phydroperiod/](https://pypi.org/project/phydroperiod/). Use `phydroperiod` when you already have water masks as local raster files or NumPy arrays; use the `HydroperiodAnalyzer` in ndvi2gif when you want to run the full pipeline (satellite collection → water masks → hydroperiod) server-side on Google Earth Engine.

### Methodology and context

Bustamante, J., Aragonés, D., Afán, I., Luque, C.J., Pérez-Vázquez, A., Castellanos, E.M., Díaz-Delgado, R. (2016). Predictive models of floristic composition based on flooding frequency in Mediterranean wetlands. *Remote Sensing*, 8(9), 776. https://doi.org/10.3390/rs8090776

Xu, H. (2006). Modification of normalised difference water index (NDWI) to enhance open water features in remotely sensed imagery. *International Journal of Remote Sensing*, 27(14), 3025–3033.

McFeeters, S.K. (1996). The use of the Normalized Difference Water Index (NDWI) in the delineation of open water features. *International Journal of Remote Sensing*, 17(7), 1425–1432.

Feyisa, G.L., Meilby, H., Fensholt, R., Proud, S.R. (2014). Automated Water Extraction Index: A new technique for surface water mapping using Landsat imagery. *Remote Sensing of Environment*, 140, 23–35.

Fisher, A., Flood, N., Danaher, T. (2016). Comparing Landsat water index methods for automated water classification in eastern Australia. *Remote Sensing of Environment*, 175, 167–182.
