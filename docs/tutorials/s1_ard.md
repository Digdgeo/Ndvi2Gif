# Sentinel-1 SAR Preprocessing with `S1ARDProcessor`

`S1ARDProcessor` implements an Analysis Ready Data (ARD) preprocessing pipeline for Sentinel-1 GRD imagery directly in Google Earth Engine. It applies radiometric terrain correction and speckle filtering to produce analysis-ready backscatter values.

In ndvi2gif, SAR preprocessing is automatically applied when `NdviSeasonality` is configured with `sat='S1'` and `use_sar_ard=True` (the default). `S1ARDProcessor` can also be used standalone to preprocess custom Sentinel-1 collections.

---

## Background: why preprocess SAR data?

Raw Sentinel-1 GRD backscatter values contain two main artefacts that must be corrected before analysis:

### Speckle noise
SAR images suffer from a multiplicative noise pattern called **speckle** — caused by coherent interference of signals backscattered from multiple scatterers within a resolution cell. Speckle makes individual pixels unreliable. Spatial filtering reduces speckle while attempting to preserve edges and structural features.

### Topographic distortion
In mountainous terrain, the sensor's incidence angle varies with slope and aspect, causing backscatter variations unrelated to surface properties. **Radiometric terrain correction** normalises backscatter values to the local incidence angle, making values comparable across terrain and between ascending/descending passes.

---

## Quick start

### Via `NdviSeasonality` (recommended)

For most workflows, SAR preprocessing is configured through `NdviSeasonality`:

```python
import ee
from ndvi2gif import NdviSeasonality

ee.Initialize(project='your-project-id')

roi = ee.Geometry.Rectangle([-6.55, 36.85, -6.10, 37.20])

ns = NdviSeasonality(
    roi=roi,
    sat='S1',
    index='rvi',                      # Radar Vegetation Index
    periods=12,
    start_year=2019,
    end_year=2023,
    orbit='DESCENDING',               # 'BOTH', 'ASCENDING', or 'DESCENDING'
    use_sar_ard=True,                 # enable ARD pipeline (default)
    sar_speckle_filter='REFINED_LEE', # speckle filter algorithm
    sar_terrain_correction=True,      # radiometric terrain correction
    sar_terrain_model='VOLUME',       # scattering model
)

collection = ns.get_year_composite()
```

### Standalone use

`S1ARDProcessor` can also preprocess a custom `ee.ImageCollection`:

```python
from ndvi2gif import S1ARDProcessor

processor = S1ARDProcessor(
    speckle_filter='REFINED_LEE',
    terrain_correction=True,
    terrain_flattening_model='VOLUME',
    dem='COPERNICUS_30',
    format='LINEAR',
)

s1_collection = (
    ee.ImageCollection('COPERNICUS/S1_GRD')
    .filterBounds(roi)
    .filterDate('2022-09-01', '2023-08-31')
    .filter(ee.Filter.eq('instrumentMode', 'IW'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
)

processed = s1_collection.map(processor.process_image)
```

---

## Speckle filters

Five speckle filter algorithms are available. All operate on the VV and VH polarisation bands independently.

### `REFINED_LEE` (recommended)

An adaptive edge-preserving filter. Detects edges using directional gradient kernels and applies filtering only in homogeneous areas, preserving structural boundaries.

```python
processor = S1ARDProcessor(speckle_filter='REFINED_LEE', speckle_filter_kernel_size=7)
```

Best for: vegetation monitoring, land cover analysis, general-purpose SAR analysis.

### `LEE`

Standard Lee adaptive filter. Uses local mean and variance to compute a weight that balances between the original pixel value and the local mean. Preserves edges better than simple averaging but less so than Refined Lee.

```python
processor = S1ARDProcessor(speckle_filter='LEE', speckle_filter_kernel_size=7)
```

Best for: moderate terrain, when Refined Lee is computationally too expensive.

### `GAMMA_MAP`

Gamma Maximum A Posteriori filter. Assumes a Gamma distribution for SAR intensity data (which is physically motivated). Particularly effective at preserving texture information.

```python
processor = S1ARDProcessor(speckle_filter='GAMMA_MAP')
```

Best for: texture analysis, forest structure assessment.

### `LEE_SIGMA`

Lee Sigma filter. Identifies "point targets" (pixels significantly brighter than their neighbourhood) and preserves them, while averaging the rest within a sigma range.

```python
processor = S1ARDProcessor(speckle_filter='LEE_SIGMA')
```

Best for: mixed scenes with point targets (buildings, infrastructure).

### `BOXCAR`

Simple spatial mean filter. Fast but does not adapt to edges — blurs structural features. Useful as a quick baseline or for very noisy data.

```python
processor = S1ARDProcessor(speckle_filter='BOXCAR', speckle_filter_kernel_size=5)
```

### Kernel size

The `speckle_filter_kernel_size` parameter controls the neighbourhood window size (in pixels, must be odd). Larger kernels produce smoother images but reduce spatial resolution:

- **3×3**: minimal smoothing, preserves fine detail
- **7×7** (default): good balance
- **11×11**: heavy smoothing, suitable for coarse land cover mapping

### Disabling speckle filtering

```python
processor = S1ARDProcessor(speckle_filter=None)  # no filtering
```

---

## Radiometric terrain correction

Terrain correction compensates for backscatter variations caused by topographic slope and aspect.

```python
processor = S1ARDProcessor(
    terrain_correction=True,
    terrain_flattening_model='VOLUME',
    dem='COPERNICUS_30',
)
```

### Scattering models

| `terrain_flattening_model=` | Physical model | Recommended for |
|---|---|---|
| `'VOLUME'` (default) | Volume scattering — energy distributed within canopy | Forests, crops, dense vegetation |
| `'SURFACE'` | Surface scattering — energy reflected from top layer | Bare soil, rock, water, urban |

The correction is clamped between 0.5 and 2.0 to prevent overcorrection in extreme terrain.

### DEMs

| `dem=` | Resolution | Notes |
|---|---|---|
| `'COPERNICUS_30'` (default) | 30 m | Recommended — better quality, global coverage |
| `'COPERNICUS_90'` | 90 m | Faster, suitable for coarse analysis |
| `'SRTM_30'` | 30 m | USGS SRTM, good coverage |
| `'SRTM_90'` | 90 m | Legacy choice |

Terrain correction is most important for areas with >10° slopes. For flat terrain (wetlands, plains), it can be disabled to save computation:

```python
processor = S1ARDProcessor(terrain_correction=False)
```

---

## Output format

```python
# Linear scale (default) — required for most indices (RVI, ratios)
processor = S1ARDProcessor(format='LINEAR')

# Decibel scale — better for visualisation, compresses dynamic range
processor = S1ARDProcessor(format='DB')

# Manual conversion after processing
processed_linear = s1_collection.map(processor.process_image)
processed_db = processed_linear.map(processor.to_db)
```

---

## Processing pipeline order

`process_image()` applies the steps in this fixed order:

1. **Terrain correction** (if `terrain_correction=True`) — must come before speckle filtering to avoid introducing terrain-dependent artefacts into the filter statistics.
2. **Speckle filtering** (if `speckle_filter` is set)
3. **Format conversion** (if `format='DB'`)

```python
# The pipeline is applied via .map()
processed = s1_collection.map(processor.process_image)
```

---

## Orbit selection

In `NdviSeasonality`, the `orbit` parameter controls which Sentinel-1 passes are included:

```python
# Both orbits — maximum temporal coverage (default)
ns = NdviSeasonality(sat='S1', orbit='BOTH')

# Ascending only — geometrically consistent for temporal analysis
ns = NdviSeasonality(sat='S1', orbit='ASCENDING')

# Descending only — different look angle, better for some terrain configurations
ns = NdviSeasonality(sat='S1', orbit='DESCENDING')
```

For temporal consistency in time series analysis, using a single orbit direction is recommended — mixing ascending and descending passes can introduce geometric biases in backscatter values.

---

## SAR indices in `NdviSeasonality`

When `sat='S1'`, the following indices can be used as the `index` parameter:

| Index | Formula | Description |
|---|---|---|
| `vv` | σ°_VV | VV polarisation backscatter |
| `vh` | σ°_VH | VH polarisation backscatter |
| `vv_vh_ratio` | σ°_VV / σ°_VH | Polarisation ratio — sensitive to surface roughness |
| `rvi` | 4σ°_VH / (σ°_VV + σ°_VH) | Radar Vegetation Index — correlates with vegetation volume |
| `rfdi` | (σ°_VV − σ°_VH) / (σ°_VV + σ°_VH) | Radar Forest Degradation Index |

```python
ns_rvi = NdviSeasonality(roi=roi, sat='S1', index='rvi', periods=4)
ns_rfdi = NdviSeasonality(roi=roi, sat='S1', index='rfdi', periods=4)
```

---

## Combining SAR and optical data

A common workflow combines Sentinel-2 optical composites with Sentinel-1 SAR for land cover classification:

```python
from ndvi2gif import NdviSeasonality, LandCoverClassifier

# Optical features
ns_s2 = NdviSeasonality(roi=roi, sat='S2', index='ndvi', periods=4,
                         start_year=2021, end_year=2022)

# SAR features
ns_s1 = NdviSeasonality(roi=roi, sat='S1', index='rvi', periods=4,
                         start_year=2021, end_year=2022, orbit='DESCENDING')

# Build a classifier on optical data, then extend the feature stack manually
# or use the optical NdviSeasonality as the base and add SAR bands as extra features
lcc = LandCoverClassifier(ns_s2)
feature_stack = lcc.create_feature_stack(indices=['ndvi', 'mndwi', 'nbr'])
```

---

## Tips

- **Start with Refined Lee** — it is the most versatile filter for vegetation monitoring.
- **Use single orbit for time series** — mixing orbits introduces look-angle-dependent backscatter changes that can be mistaken for real change.
- **Enable terrain correction for any terrain with >5° slopes** — even moderate topography introduces significant backscatter bias in C-band SAR.
- **Use `format='LINEAR'` for indices** — RVI and ratios require linear-scale values. Convert to dB only for visualisation.
- **Copernicus 30m DEM** is generally superior to SRTM for terrain correction — it has fewer voids and better accuracy at high latitudes.

---

## API reference

| Method | Returns | Description |
|---|---|---|
| `process_image(image)` | `ee.Image` | Apply full ARD pipeline to one image |
| `apply_terrain_correction(image)` | `ee.Image` | Radiometric terrain correction only |
| `apply_speckle_filter(image)` | `ee.Image` | Speckle filter only |
| `to_db(image)` | `ee.Image` | Convert linear scale to dB |

---

## References

Vollrath, A., Mullissa, A., Reiche, J. (2020). Angular-based radiometric slope correction for Sentinel-1 on Google Earth Engine. *Remote Sensing*, 12(11), 1867. https://doi.org/10.3390/rs12111867

Lee, J.S. (1980). Digital image enhancement and noise filtering by use of local statistics. *IEEE Transactions on Pattern Analysis and Machine Intelligence*, 2(2), 165–168.

Lee, J.S. (1983). Digital image smoothing and the sigma filter. *Computer Vision, Graphics, and Image Processing*, 24(2), 255–269.

Lopes, A., Nezry, E., Touzi, R., Laur, H. (1993). Structure detection and statistical adaptive speckle filtering in SAR images. *International Journal of Remote Sensing*, 14(9), 1735–1758.
