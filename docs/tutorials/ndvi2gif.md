# Temporal Compositing with `NdviSeasonality`

`NdviSeasonality` is the core engine of ndvi2gif. It configures the satellite collection, computes spectral and radar indices, and generates temporal composite images — one per period per year. All other classes (`HydroperiodAnalyzer`, `TimeSeriesAnalyzer`, `LandCoverClassifier`) take an existing `NdviSeasonality` instance and build on top of it.

---

## What are temporal composites?

Rather than working with individual satellite scenes (which are cloudy, noisy, and irregularly spaced), temporal compositing aggregates all cloud-free observations within a time window into a single representative image using a statistical reducer (max, median, mean, percentile).

A **seasonal composite** for NDVI answers: "what was the maximum NDVI in spring 2022 across every pixel?"

ndvi2gif supports:

- **4 periods/year** — meteorological seasons (winter, spring, summer, autumn)
- **12 periods/year** — calendar months (january through december)
- **24 periods/year** — bi-monthly (~15 days each)
- **Any N** — equal divisions of the year (p1, p2, …, pN)

---

## Quick start

```python
import ee
import geemap
from ndvi2gif import NdviSeasonality

ee.Initialize(project='your-project-id')

roi = ee.Geometry.Rectangle([-6.55, 36.85, -6.10, 37.20])

ns = NdviSeasonality(
    roi=roi,
    sat='S2',
    start_year=2019,
    end_year=2023,
    periods=4,
    index='ndvi',
    key='max',
)

# Get composite ImageCollection (one image per year, each with 4 period bands)
collection = ns.get_year_composite()
print(collection.size().getInfo())  # → 5 (one per year)
```

---

## Region of interest (ROI)

The `roi` parameter accepts multiple formats:

```python
# Direct EE geometry
roi = ee.Geometry.Rectangle([-6.55, 36.85, -6.10, 37.20])
roi = ee.Geometry.Point([-5.97, 37.27]).buffer(5000)  # 5 km buffer

# Shapefile or GeoJSON
ns = NdviSeasonality(roi='study_area.shp')
ns = NdviSeasonality(roi='study_area.geojson')

# eLTER DEIMS site ID
ns = NdviSeasonality(roi='deimsid/11696159-444f-4e06-b537-d4c5c0a4e97d')

# Sentinel-2 MGRS tile
ns = NdviSeasonality(roi='s2:30TXN')

# Landsat WRS-2 path/row
ns = NdviSeasonality(roi='wrs:200,32')

# Drawn features from a geemap Map
ns = NdviSeasonality(roi=Map.draw_features)

# None → default region (Andalusia, Spain — useful for quick tests)
ns = NdviSeasonality()
```

---

## Satellite sensors

| `sat=` | Sensor | Resolution | Archive start | Best for |
|---|---|---|---|---|
| `'S2'` | Sentinel-2 MSI | 10–20 m | 2015 | Detailed land cover, vegetation |
| `'Landsat'` | Landsat 4–5–7–8–9 | 30 m | 1982 | Long time series |
| `'MODIS'` | Terra/Aqua | 500 m | 2000 | Large-scale, daily coverage |
| `'S1'` | Sentinel-1 SAR | 10 m | 2014 | All-weather, structural changes |
| `'S3'` | Sentinel-3 OLCI | 300 m | 2016 | Ocean/coastal, water quality |

```python
# Landsat for long-term analysis (1990–2023)
ns_ls = NdviSeasonality(roi=roi, sat='Landsat', start_year=1990, end_year=2023)

# MODIS for continental-scale monitoring
ns_mod = NdviSeasonality(roi=big_roi, sat='MODIS', start_year=2000, end_year=2023)

# Sentinel-1 SAR (index must be SAR-compatible)
ns_s1 = NdviSeasonality(roi=roi, sat='S1', index='vh', periods=12)
```

---

## Statistical reducers

The `key` parameter controls how observations within each period are combined:

| `key=` | Description | Use case |
|---|---|---|
| `'max'` | Maximum value | Vegetation peak, cloud-free guarantee for NDVI |
| `'median'` | Median value | Robust to outliers |
| `'mean'` | Mean value | Smooth temporal profiles |
| `'percentile'` | Custom percentile (set with `percentile=`) | Flexible thresholding |
| `'sum'` | Sum | Precipitation, flood accumulation |

```python
# 90th percentile — similar to max but more robust
ns = NdviSeasonality(roi=roi, key='percentile', percentile=90)

# Median — good for SAR or areas with many outliers
ns = NdviSeasonality(roi=roi, key='median')
```

---

## Spectral and radar indices

Over 40 indices are available. Use `get_available_indices()` to check what is supported for each sensor:

```python
ns = NdviSeasonality(sat='S2')
print(ns.get_available_indices())           # S2 indices
print(ns.get_available_indices('Landsat'))  # Landsat indices

# Overview of all sensors at once
all_indices = ns.get_all_available_indices()
for sensor, indices in all_indices.items():
    print(f"{sensor}: {len(indices)} indices")
```

### Index categories

**Vegetation (optical)**

| Index | Description |
|---|---|
| `ndvi` | Normalized Difference Vegetation Index |
| `evi` | Enhanced Vegetation Index |
| `savi` | Soil-Adjusted Vegetation Index |
| `gndvi` | Green NDVI (chlorophyll-sensitive) |
| `wdrvi` | Wide Dynamic Range Vegetation Index |
| `lai` | Leaf Area Index |

**Water**

| Index | Description |
|---|---|
| `ndwi` | McFeeters water index |
| `mndwi` | Modified NDWI (Xu 2006) |
| `awei` | AWEI with shadow correction |
| `aweinsh` | AWEI without shadow correction |
| `wi2015` | Water Index 2015 |

**Burn / Fire**

| Index | Description |
|---|---|
| `nbr` | Normalized Burn Ratio |
| `nbri` | NBR Index |

**Urban / Built-up**

| Index | Description |
|---|---|
| `ndbi` | Normalized Difference Built-up Index |

**Red Edge (S2 only)**

| Index | Description |
|---|---|
| `ndre` | Normalized Difference Red Edge |
| `ireci` | Inverted Red Edge Chlorophyll Index |
| `mcari` | Modified Chlorophyll Absorption Ratio |
| `reip` | Red Edge Inflection Point |
| `s2rep` | S2 Red Edge Position |
| `cire` | Chlorophyll Index Red Edge |
| `mtci` | MERIS Terrestrial Chlorophyll Index |

**Ocean / Water quality (S3)**

| Index | Description |
|---|---|
| `fai` | Floating Algae Index |
| `ndci` | Normalized Difference Chlorophyll Index |
| `tsi` | Trophic State Index |
| `cdom` | Coloured Dissolved Organic Matter |
| `turbidity` | Turbidity index |
| `spm` | Suspended Particulate Matter |
| `kd490` | Diffuse attenuation coefficient |

**SAR (S1)**

| Index | Description |
|---|---|
| `vv` | VV backscatter |
| `vh` | VH backscatter |
| `rvi` | Radar Vegetation Index |
| `vv_vh_ratio` | VV/VH ratio |
| `rfdi` | Radar Forest Degradation Index |

---

## Cloud masking

Cloud masking is applied during collection setup, before any composite is computed.

```python
ns = NdviSeasonality(
    roi=roi,
    sat='S2',
    cloud_filter=True,     # scene-level: remove scenes with >max_cloud_cover% clouds
    max_cloud_cover=20,    # percentage threshold (default 20)
    scl_mask=True,         # pixel-level: SCL cloud/shadow mask (S2 only, recommended)
)
```

`scl_mask=True` (default for S2) applies a pixel-level mask using the Sentinel-2 Scene Classification Layer, removing clouds (classes 8, 9, 10), cloud shadows (class 3), and saturated pixels (class 1). This is more accurate than the legacy QA60 band used by `scl_mask=False`.

For Landsat, pixel-level cloud masking uses the QA_PIXEL band (Fmask). For MODIS, the StateQA band. `scl_mask` has no effect on Landsat, MODIS, or S1.

---

## Generating composites

### `get_year_composite()` — all years

Returns an `ee.ImageCollection` with one multi-band image per year. Each image has as many bands as `periods`, named after the period (`winter`, `spring`, … or `january`, `february`, …).

```python
collection = ns.get_year_composite()

# First image (first year), first band (first period)
first = ee.Image(collection.first())
print(first.bandNames().getInfo())  # ['winter', 'spring', 'summer', 'autumn']
```

Optionally get the image count per period (useful for quality control):

```python
collection, counts = ns.get_year_composite(
    return_counts=True,
    count_mode='unique_dates',   # 'granules' or 'unique_dates'
)
for c in counts[:4]:
    print(c['year'], c['period_name'], c['images_count'])
```

### `get_period_composite()` — single period

Returns a single `ee.Image` for a specific year and period index:

```python
# Winter 2021 (period_idx=0 for 4-period config)
winter_2021 = ns.get_period_composite(2021, 0)

# July 2020 (period_idx=6 for monthly config)
ns_monthly = NdviSeasonality(periods=12, start_year=2020, end_year=2022)
july_2020 = ns_monthly.get_period_composite(2020, 6)
```

---

## Visualising on a map

```python
Map = geemap.Map()
Map.centerObject(roi, zoom=10)

# Add one year's composite
year_img = ee.Image(collection.first())
vis = {'bands': ['spring'], 'min': 0.2, 'max': 0.8, 'palette': ['white', 'green']}
Map.addLayer(year_img, vis, 'Spring NDVI 2019')
Map
```

---

## Generating GIFs

```python
# Animated GIF of seasonal composites over all years
ns.get_gif('donana_ndvi_seasonal.gif')
```

---

## Exporting to GeoTIFF

```python
# Export to Google Drive
ns.get_export(
    folder='ndvi_exports',
    description='donana_ndvi_2019_2023',
    scale=10,
    crs='EPSG:25829',
)
```

---

## Sentinel-1 SAR configuration

When `sat='S1'`, additional parameters control the SAR preprocessing pipeline:

```python
ns_s1 = NdviSeasonality(
    roi=roi,
    sat='S1',
    index='rvi',              # Radar Vegetation Index
    periods=12,
    orbit='DESCENDING',       # 'BOTH', 'ASCENDING', or 'DESCENDING'
    normalize_sar=False,      # Z-score normalisation
    use_sar_ard=True,         # Enable advanced ARD preprocessing
    sar_speckle_filter='REFINED_LEE',
    sar_terrain_correction=True,
    sar_terrain_model='VOLUME',  # 'VOLUME' (vegetation) or 'SURFACE' (bare soil)
)
```

For detailed SAR preprocessing options, see the [S1ARDProcessor tutorial](s1_ard.md).

---

## API reference

| Method | Returns | Description |
|---|---|---|
| `get_year_composite(...)` | `ee.ImageCollection` | Multi-band composites for all years |
| `get_period_composite(year, period_idx)` | `ee.Image` | Composite for one period of one year |
| `get_available_indices(satellite)` | `list[str]` | Indices available for a sensor |
| `get_all_available_indices()` | `dict` | Index catalogue for all sensors |
| `get_gif(filename)` | — | Export animated GIF |
| `get_export(folder, ...)` | — | Export GeoTIFF to Drive |

---

## References

Gorelick, N., Hancher, M., Dixon, M., Ilyushchenko, S., Thau, D., Moore, R. (2017). Google Earth Engine: Planetary-scale geospatial analysis for everyone. *Remote Sensing of Environment*, 202, 18–27.

Tucker, C.J. (1979). Red and photographic infrared linear combinations for monitoring vegetation. *Remote Sensing of Environment*, 8(2), 127–150.
