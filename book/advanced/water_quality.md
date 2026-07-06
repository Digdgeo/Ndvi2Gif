# Water Quality Monitoring in Reservoirs

This tutorial shows how to use **ndvi2gif** to monitor inland water bodies with Sentinel-2: detect the flooded extent each month, compute a monthly chlorophyll-a proxy, and study its inter-annual evolution. The workflow relies on the standard `NdviSeasonality` composites and on `TimeSeriesAnalyzer` for trend detection.

A complete, runnable companion notebook with a reservoirs shapefile (CHG — Confederación Hidrográfica del Guadalquivir) is available in [`examples_notebooks/Water_Quality_Embalses.ipynb`](https://github.com/Digdgeo/Ndvi2Gif/blob/master/examples_notebooks/Water_Quality_Embalses.ipynb).

---

## Overview

Reservoirs are heterogeneous systems: shorelines move with the water level, turbidity varies seasonally, and eutrophication episodes (cyanobacterial blooms) are concentrated in the warmer months. Any water-quality analysis therefore needs two layers:

1. **Water detection** — a binary mask of "where is water *right now*".
2. **Quality indicator** — a biogeochemical proxy evaluated **only on water pixels** of that same period.

`ndvi2gif` makes this straightforward: water indices (NDWI, MNDWI, AWEI, WI2015) and the Red-Edge chlorophyll index NDCI are all available out of the box, all using the same standardised band names.

| Index | Formula (S2) | Role |
|---|---|---|
| `ndwi` | `(Green − NIR) / (Green + NIR)` | Classic water detection |
| `mndwi` | `(Green − SWIR1) / (Green + SWIR1)` | **Default water mask** — more robust |
| `awei` / `aweinsh` | Multi-band (Blue, Green, NIR, SWIR1, SWIR2) | Automatic water extraction |
| `wi2015` | NIR + SWIR1 + SWIR2 linear combination | Complementary water index |
| `ndci` | `(RedEdge1 − Red) / (RedEdge1 + Red)` | **Chlorophyll-a proxy — Sentinel-2 only** |

NDCI exploits the Red-Edge band B5 (705 nm), exclusive to Sentinel-2 among the sensors ndvi2gif supports, and is sensitive to chlorophyll-a concentration and potential cyanobacterial blooms.

---

## 1. Setup and ROI

The ROI is a reservoir polygon. Any `ee.Geometry`, shapefile, or GeoJSON works. To keep this page fully reproducible, we load a **public Earth Engine asset** with the CHG reservoir polygons (Confederación Hidrográfica del Guadalquivir), pick the largest one — the **Embalse de Iznájar** (~24 km², the largest in Andalusia) — and buffer it by 500 m:

```python
import ee, geemap
from ndvi2gif import NdviSeasonality, TimeSeriesAnalyzer

ee.Initialize(project='your-project-id')

# Public asset: CHG reservoir polygons (Guadalquivir basin).
reservoirs = ee.FeatureCollection(
    'projects/ee-digdgeografo/assets/CHG_CaptacionesEmbalses'
)

# Pick the largest reservoir (area computed server-side) and buffer it by 500 m.
reservoirs = reservoirs.map(lambda f: f.set('area_m2', f.geometry().area(1)))
largest = ee.Feature(reservoirs.sort('area_m2', False).first())
roi = largest.geometry().buffer(500)
```

> **Why the buffer?** In water-quality work, forgetting the buffer silently clips the shoreline when the reservoir is full, and creates a ring of dry pixels at low water that contaminate the mask. Buffering by 200–1000 m fixes both issues.
>
> **Using your own vector file.** Any shapefile or GeoJSON path works in place of the asset — just pass it as `roi=` to `NdviSeasonality`, or load it with `geemap.shp_to_ee(...)` / `geopandas` if you need to select and buffer a specific feature first.

---

## 2. Monthly water composites

`periods=12` produces one band per calendar month per year. The `median` reducer is the right choice for water indices: it is robust to residual cloud shadows and sun-glint outliers that would pull a mean toward false positives.

```python
s2_mndwi = NdviSeasonality(
    roi=roi,
    periods=12,
    start_year=2020,
    end_year=2025,
    sat='S2',
    key='median',
    index='mndwi',
)
composite_mndwi = s2_mndwi.get_year_composite()
# -> ee.ImageCollection of 6 images (2020..2025), each with 12 bands
#    named january..december
```

Repeat with `index='ndwi'`, `'aweinsh'`, or `'wi2015'` to compare detectors. Quick visual check for July, inter-annual median:

```python
water_viz = {
    'bands': ['july'], 'min': -0.3, 'max': 0.5,
    'palette': ['#d73027', '#fc8d59', '#fee090',
                '#e0f3f8', '#74add1', '#313695'],
}

Map = geemap.Map()
Map.centerObject(composite_mndwi, zoom=11)
Map.addLayer(composite_mndwi.select('july').median(), water_viz,
             'MNDWI July (interannual median)')
Map
```

> **Choosing the detector.** For most Iberian reservoirs `mndwi > 0.0` is a very good default. In highly turbid water it may under-detect — try `aweinsh` or lower the threshold (e.g. `-0.1`). AWEInsh is usually preferable over AWEI in reservoirs in flat terrain; AWEI's shadow term is designed for mountainous settings and can over-suppress legitimate water.

---

## 3. Binary water masks

Thresholding the monthly composite gives a binary mask per month and per year. Keep the threshold explicit — it is the single most sensitive parameter in water-quality work.

```python
WATER_THRESHOLD = 0.0  # typical range: 0.0 to 0.1

def build_water_masks(composite, threshold=WATER_THRESHOLD):
    def to_binary(img):
        return img.gt(threshold).rename(img.bandNames())
    return composite.map(to_binary)

water_masks = build_water_masks(composite_mndwi)
```

Each image in the returned collection still has 12 bands (`january`..`december`), but each pixel is now 1 (water) or 0 (not water).

### Flooded surface per month

Summing the mask × pixel area gives the monthly flooded surface in km²:

```python
MONTHS = ['january','february','march','april','may','june',
          'july','august','september','october','november','december']

imgs  = water_masks.toList(water_masks.size())
years = list(range(s2_mndwi.start_year, s2_mndwi.end_year + 1))
records = []

for i, year in enumerate(years):
    sums = ee.Image(imgs.get(i)).reduceRegion(
        reducer=ee.Reducer.sum(),
        geometry=s2_mndwi.roi,
        scale=10, maxPixels=1e9,
    ).getInfo()
    for m in MONTHS:
        if sums.get(m) is not None:
            records.append({'year': year, 'month': m,
                            'area_km2': sums[m] * 10 / 1e6})
# Note: `scale=10` for S2; factor 10/1e6 converts pixel-count at 10 m
# (100 m² per pixel) to km² (divide by 10 000).
```

Plotting `area_km2` per month with one line per year makes the hydrological regime of the reservoir immediately visible: filling in winter/spring, draw-down in summer, and the inter-annual variability of each month.

---

## 4. Chlorophyll-a with NDCI

The Normalised Difference Chlorophyll Index (Mishra & Mishra, 2012) uses Red-Edge 1 (B5) and Red (B4) and is the standard chlorophyll-a proxy for inland waters from Sentinel-2.

```python
s2_ndci = NdviSeasonality(
    roi=roi,
    periods=12,
    start_year=2020,
    end_year=2025,
    sat='S2',
    key='median',
    index='ndci',      # Sentinel-2 only
)
composite_ndci = s2_ndci.get_year_composite()
```

### Applying the water mask to NDCI

NDCI is only meaningful *over water*. Any interpretation of NDCI on dry, vegetated, or mixed pixels is wrong by construction. Apply the monthly water mask before reducing or visualising:

```python
def mask_water_quality(index_coll, mask_coll, threshold=0.0):
    idx_list  = index_coll.toList(index_coll.size())
    mask_list = mask_coll.toList(mask_coll.size())
    n = index_coll.size().min(mask_coll.size())

    def mask_i(i):
        idx  = ee.Image(idx_list.get(i))
        mask = ee.Image(mask_list.get(i)).gt(threshold)
        return idx.updateMask(mask)

    return ee.ImageCollection(
        ee.List.sequence(0, n.subtract(1)).map(mask_i)
    )

ndci_water = mask_water_quality(composite_ndci, composite_mndwi)
```

This pairs each year's NDCI with *that same year's* MNDWI mask — not with a single static water polygon. It correctly handles reservoirs whose shape changes season to season.

### NDCI interpretation thresholds

| NDCI | Interpretation |
|---|---|
| < 0 | Clear water, low productivity |
| 0.0 – 0.1 | Moderate chlorophyll |
| 0.1 – 0.2 | Elevated — monitor |
| > 0.2 | Potential algal bloom — alert |

These are practical reference values for Mediterranean reservoirs; local calibration against *in situ* samples is always recommended if you need absolute concentrations.

### Visualising

```python
chl_palette = ['#053061', '#2166ac', '#4393c3', '#92c5de', '#d1e5f0',
               '#f7f7f7', '#fddbc7', '#f4a582', '#d6604d', '#b2182b', '#67001f']

Map.addLayer(
    ndci_water.select('july').median(),
    {'bands': ['july'], 'min': -0.2, 'max': 0.3, 'palette': chl_palette},
    'NDCI July — chlorophyll-a (water only)'
)
```

July is a good default month to inspect in Iberia: water temperature is high and solar radiation peaks, so any tendency to algal blooms shows up clearly.

---

## 5. Time series of water quality

Extract mean/median/std NDCI per month per year directly from the masked composite:

```python
MONTH_NUM = {m: i + 1 for i, m in enumerate(MONTHS)}
imgs = ndci_water.toList(ndci_water.size())

stats_reducer = (
    ee.Reducer.mean()
    .combine(ee.Reducer.median(), sharedInputs=True)
    .combine(ee.Reducer.stdDev(), sharedInputs=True)
)

rows = []
for i, year in enumerate(years):
    stats = ee.Image(imgs.get(i)).reduceRegion(
        reducer=stats_reducer,
        geometry=s2_ndci.roi,
        scale=20, maxPixels=1e9,
    ).getInfo()
    for m in MONTHS:
        rows.append({
            'year': year, 'month': m, 'month_num': MONTH_NUM[m],
            'mean':   stats.get(f'{m}_mean'),
            'median': stats.get(f'{m}_median'),
            'std':    stats.get(f'{m}_stdDev'),
        })
```

> **Scale.** NDCI is reduced at `scale=20` because the Red-Edge B5 band is natively 20 m. Reducing at 10 m forces a resample and does not improve the signal.

A typical plot stacks two panels: one curve per year (phenology-like view of each year's chlorophyll cycle) and an `year × month` heatmap (highlights anomalous years at a glance). The companion notebook contains the full Matplotlib code.

---

## 6. Trend detection with `TimeSeriesAnalyzer`

`TimeSeriesAnalyzer` wraps an `NdviSeasonality` instance and provides trend tests on its composites:

```python
ts = TimeSeriesAnalyzer(s2_ndci)

series = ts.extract_time_series(reducer='mean', scale=20)
trend  = ts.analyze_trend()

mk = trend['mann_kendall']
print(f"tau = {mk['tau']:.3f}, p = {mk['p_value']:.4f}")
print(trend['interpretation'])

ts.plot_comprehensive_analysis()
```

The Mann-Kendall test is non-parametric and robust to the seasonal structure of NDCI. A positive, significant `tau` on chlorophyll-a means the reservoir is becoming more eutrophic over time — this is the kind of signal reservoir managers need.

---

## 7. Reference table of water indices in ndvi2gif

| Key | Name | Bands | Sensors | Purpose |
|---|---|---|---|---|
| `ndwi` | NDWI (Gao 1996) | Green, NIR | All optical | Water detection (classic) |
| `mndwi` | Modified NDWI (Xu 2006) | Green, SWIR1 | All optical | **Recommended water mask** |
| `awei` | AWEI shade (Feyisa 2014) | Blue, Green, NIR, SWIR1, SWIR2 | All optical | Automatic water extraction, shadowed terrain |
| `aweinsh` | AWEI no-shadow (Feyisa 2014) | Green, NIR, SWIR1, SWIR2 | All optical | Automatic water extraction, flat terrain |
| `wi2015` | Water Index 2015 (Fisher 2016) | NIR, SWIR1, SWIR2 | All optical | Complementary water index |
| `ndci` | NDCI (Mishra & Mishra 2012) | RedEdge1, Red | **S2 only** | Chlorophyll-a / cyanobacteria |
| `ndmi` | NDMI (Gao 1996) | NIR, SWIR1 | All optical | Moisture (riparian/aquatic vegetation) |

Sentinel-3 OLCI (300 m) adds dedicated water-quality indices — `oci`, `tsi`, `cdom`, `turbidity`, `spm`, `kd490`, `floating_algae` — useful for large lakes and reservoirs where the coarser resolution is acceptable.

---

## 8. Exporting results

The binary water masks and the masked NDCI are ordinary `ee.Image` / `ee.ImageCollection` objects, so the standard Earth Engine export API applies:

```python
task = ee.batch.Export.image.toDrive(
    image=water_masks.median(),          # interannual median of monthly masks
    description='water_mask_monthly_mndwi',
    folder='ndvi2gif_outputs',
    fileNamePrefix='water_mask_mndwi_2020_2025',
    region=s2_mndwi.roi,
    scale=10,
    crs='EPSG:4326',
    maxPixels=1e10,
)
task.start()
```

Export the NDCI composite the same way (at `scale=20`). Monitor running exports at [code.earthengine.google.com/tasks](https://code.earthengine.google.com/tasks).

---

## Tips and caveats

- **Always mask NDCI with the *monthly* water mask, not a static polygon.** Reservoirs change shape; a static polygon contaminates the dry season with vegetation signal.
- **Threshold calibration.** MNDWI > 0.0 works for most Iberian reservoirs. Turbid or shallow systems may need 0.05–0.1. Document the threshold you used — results are not comparable across studies otherwise.
- **Median over mean.** For both water indices and NDCI, `key='median'` suppresses cloud-shadow and sun-glint outliers far better than `mean`.
- **Scale factors.** Sentinel-2 scaling (0.0001) is applied automatically by `NdviSeasonality`. You do not need to scale bands manually before computing indices.
- **Buffer your ROI.** A 200–1000 m buffer captures the full shoreline range and the transition zone, and makes the statistics representative of the reservoir as a hydrological unit, not just its polygon.

---

## References

- Gao, B. (1996). NDWI — A normalized difference water index for remote sensing of vegetation liquid water from space. *Remote Sensing of Environment*, 58(3), 257–266.
- Xu, H. (2006). Modification of normalised difference water index (NDWI) to enhance open water features. *International Journal of Remote Sensing*, 27(14), 3025–3033.
- Feyisa, G.L., Meilby, H., Fensholt, R., Proud, S.R. (2014). Automated Water Extraction Index: A new technique for surface water mapping using Landsat imagery. *Remote Sensing of Environment*, 140, 23–35.
- Fisher, A., Flood, N., Danaher, T. (2016). Comparing Landsat water index methods for automated water classification in eastern Australia. *Remote Sensing of Environment*, 175, 167–182.
- Mishra, S., Mishra, D.R. (2012). Normalized difference chlorophyll index: a novel model for remote estimation of chlorophyll-a concentration in turbid productive waters. *Remote Sensing of Environment*, 117, 394–406.
