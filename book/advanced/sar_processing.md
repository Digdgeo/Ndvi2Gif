# Sentinel-1 SAR Processing with `S1ARDProcessor`

Sentinel-1 C-band SAR sees through clouds and night, and therefore fills the gap left by optical sensors in cloud-prone areas (boreal, tropical, monsoon) and in rapidly changing systems (crops, wetlands, floods). But raw GRD scenes are not directly comparable in time: they are affected by speckle noise and by topographic distortions that depend on the scene's acquisition geometry.

`S1ARDProcessor` turns Sentinel-1 GRD data into **Analysis Ready Data (ARD)** inside Google Earth Engine. It applies radiometric terrain correction and one of five speckle filters, producing time-consistent VV and VH backscatter that can be used for indices, composites, classification, or time-series analysis.

The processor follows the angular-based radiometric slope correction method of Vollrath et al. (2020).

---

## What ARD gives you

| Problem | Consequence | What ARD does |
|---|---|---|
| Speckle noise | High pixel-to-pixel variance, unreliable thresholds | Adaptive spatial filtering (Lee variants, Gamma MAP, Lee Sigma, Boxcar) |
| Terrain distortion | Same land cover with different backscatter on pole-facing vs sun-facing slopes | Radiometric correction from local incidence angle (slope, aspect, orbit) |
| Geometric inconsistency between orbits | Asc vs desc scenes not comparable | Filter by `orbit` keeping only one pass direction |
| Raw linear power is hard to interpret | Non-Gaussian dynamic range | Optional conversion to decibel scale |

In short: without ARD, pixel values between two dates are not directly comparable. With ARD, they are.

---

## Two ways to use it

`S1ARDProcessor` can be used **directly** on any `ee.ImageCollection` of Sentinel-1 GRD scenes, or indirectly **through `NdviSeasonality`** by passing `use_sar_ard=True`. The latter is the usual path: it reuses the same ROI, date range, and orbit logic that you already configure for optical sensors.

```python
import ee
from ndvi2gif import NdviSeasonality, S1ARDProcessor

ee.Initialize(project='your-project-id')
```

---

## Quick start — through `NdviSeasonality`

```python
roi = ee.Geometry.Rectangle([-6.55, 36.85, -6.10, 37.20])

ns = NdviSeasonality(
    roi=roi,
    sat='S1',
    index='vh',
    start_year=2022,
    end_year=2023,
    periods=12,
    key='median',
    orbit='DESCENDING',       # pick one pass direction for temporal consistency
    use_sar_ard=True,          # enable ARD pipeline
    sar_speckle_filter='REFINED_LEE',
    sar_terrain_correction=True,
    sar_terrain_model='VOLUME',
)

composite = ns.get_year_composite()
# ee.ImageCollection: 2 images (2022, 2023), each with 12 bands (january..december)
```

That's the full workflow: orbit-filtered, terrain-corrected, speckle-filtered, reduced to monthly medians of VH backscatter. Everything happens lazily in GEE — no download.

---

## Direct use — any GRD collection

If you need ARD outside of `NdviSeasonality` (e.g. to build a custom pipeline, to feed a classifier directly, or to preprocess a single scene), instantiate the processor and map it:

```python
processor = S1ARDProcessor(
    speckle_filter='REFINED_LEE',
    speckle_filter_kernel_size=7,
    terrain_correction=True,
    terrain_flattening_model='VOLUME',
    dem='COPERNICUS_30',
    format='LINEAR',
)

s1 = (ee.ImageCollection('COPERNICUS/S1_GRD')
        .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
        .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
        .filter(ee.Filter.eq('instrumentMode', 'IW'))
        .filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))
        .filterBounds(roi)
        .filterDate('2022-01-01', '2023-01-01')
        .select(['VV', 'VH', 'angle']))

s1_ard = s1.map(processor.process_image)
```

`process_image` runs the full chain in order: **terrain correction → speckle filter → format conversion**. Terrain correction expects the `angle` band, so make sure to include it in the `select()` before mapping.

> **Angle band is required for terrain correction.** If `terrain_correction=True`, keep `'angle'` in the band selection. You can drop it afterwards with `.select(['VV', 'VH'])`.

---

## Speckle filters

Five algorithms are available, covering the usual speed-vs-quality trade-off:

| Filter | Strengths | When to use |
|---|---|---|
| `REFINED_LEE` | Detects edges, preserves structure, low texture loss | **Recommended default** — vegetation, crops, forests |
| `LEE` | Adaptive, good homogeneity/edge balance | Light, general-purpose filtering |
| `GAMMA_MAP` | Statistical, preserves texture | Areas with meaningful texture (forest canopy, urban) |
| `LEE_SIGMA` | Preserves point targets and bright features | Scenes with infrastructure, isolated bright scatterers |
| `BOXCAR` | Simple mean, very fast | Quick exploration, prototyping |
| `None` | No filter | If you intend to do your own temporal filtering |

All filters assume an Equivalent Number of Looks (ENL) of 4.4 for Sentinel-1 IW GRD. The kernel size (`speckle_filter_kernel_size`) must be odd: 3, 5, 7 (default), 9. Larger kernels smooth more but also blur structure — start at 7 and only increase if the remaining speckle is visible in your application.

> **Caveat — directional kernels for Refined Lee.** The current implementation uses fixed 3×3 directional kernels for edge detection regardless of the `speckle_filter_kernel_size`. In practice this still outperforms the plain Lee filter for S1 IW data, but if you need very large smoothing you may want to run a second pass with a larger kernel.

---

## Terrain correction

Radiometric terrain correction (RTC) compensates for the backscatter anisotropy introduced by local slope and aspect relative to the sensor's look geometry. It is essential in mountainous terrain; without it, a forest on a radar-facing slope appears brighter than the same forest on a radar-averted slope.

### Scattering model

| Model | Assumption | Use for |
|---|---|---|
| `VOLUME` | Scattering inside a 3D medium | Vegetation, crops, forests — **default** |
| `SURFACE` | Scattering from a 2D surface | Bare soil, rock, water |

The VOLUME model uses the full local-incidence-angle formulation (Vollrath et al. 2020); SURFACE uses a simpler cosine-of-slope correction. For mixed landscapes, VOLUME is usually the safer choice.

### DEM choice

| DEM | Resolution | Notes |
|---|---|---|
| `COPERNICUS_30` | 30 m | **Recommended default** — best global quality, covers high latitudes |
| `COPERNICUS_90` | 90 m | Faster, use for very large ROIs |
| `SRTM_30` | 30 m | Classic, no coverage above 60°N / below 56°S |
| `SRTM_90` | 90 m | Classic, same latitude limits |

The correction factor is clamped to the `[0.5, 2.0]` range to prevent overcorrection in very steep terrain where the angular geometry becomes unstable.

---

## Orbit filtering: why it matters

Sentinel-1 sees the same area from ascending (SSW→NNE look direction) and descending (NNE→SSW look direction) orbits. Backscatter depends on look angle, so **mixing orbits in a time series adds artificial variability** that is not related to the ground.

Always pass `orbit='ASCENDING'` **or** `orbit='DESCENDING'` when building a time series. `'BOTH'` is only appropriate when revisit frequency is more important than radiometric consistency (e.g. rapid flood mapping).

```python
# Temporal consistency (time series, trend analysis, classification)
NdviSeasonality(..., sat='S1', orbit='DESCENDING', use_sar_ard=True)

# Maximum coverage (flood mapping, single-date snapshots)
NdviSeasonality(..., sat='S1', orbit='BOTH', use_sar_ard=True)
```

Choose the direction with better coverage over your ROI — the GEE code editor shows per-orbit footprints, and for most Mediterranean sites both directions are well covered.

---

## SAR indices available

Once the collection is ARD-processed, the usual `NdviSeasonality` index interface applies:

| Key | Name | Formula | Use |
|---|---|---|---|
| `vv` | VV polarization | `VV` | Surface roughness, bare soil, water |
| `vh` | VH polarization | `VH` | Vegetation volume scattering — **default for vegetation** |
| `rvi` | Radar Vegetation Index | `4·VH / (VV+VH)` | Vegetation water content, robust to speckle |
| `vv_vh_ratio` | VV/VH ratio | `VV / VH` | Structural change, mowing, ploughing |
| `dpsvi` | Dual-pol SAR Vegetation Index | `(VV−VH) / (VV+VH)` | Dense canopies, crop growth |
| `rfdi` | Radar Forest Degradation Index | `(VV−VH) / VV` | Forest disturbance |
| `vsdi` | Vegetation Scattering Diversity Index | `sqrt((VV−VH)² + (VV+VH)²)` | Scattering diversity, biomass |

All of them accept `normalize_sar=True` at the `NdviSeasonality` level to rescale output to `[0, 1]` — useful for visualisation and for stacking with optical indices in a classifier.

---

## Full example — crop phenology with VH + RVI

```python
from ndvi2gif import NdviSeasonality

# Monthly VH and RVI over an agricultural area
ns_vh = NdviSeasonality(
    roi=roi,
    sat='S1',
    index='vh',
    start_year=2020, end_year=2023,
    periods=12,
    key='median',
    orbit='DESCENDING',
    use_sar_ard=True,
    sar_speckle_filter='REFINED_LEE',
    sar_terrain_correction=True,
    sar_terrain_model='VOLUME',
)
composite_vh = ns_vh.get_year_composite()

ns_rvi = NdviSeasonality(
    roi=roi,
    sat='S1',
    index='rvi',
    start_year=2020, end_year=2023,
    periods=12,
    key='median',
    orbit='DESCENDING',
    use_sar_ard=True,
    normalize_sar=True,
)
composite_rvi = ns_rvi.get_year_composite()
```

Add to a map to compare summer (dry, low VH) and spring (active vegetation, high VH):

```python
import geemap

Map = geemap.Map()
Map.centerObject(roi, zoom=11)

vh_viz  = {'bands': ['july'], 'min': 0.0005, 'max': 0.1,
           'palette': ['#000000', '#7f7f7f', '#ffffff']}
rvi_viz = {'bands': ['july'], 'min': 0.0, 'max': 1.0,
           'palette': ['#c7e9b4', '#41b6c4', '#253494']}

Map.addLayer(composite_vh.select('july').median(),  vh_viz,  'VH July (median)')
Map.addLayer(composite_rvi.select('july').median(), rvi_viz, 'RVI July (median)')
Map
```

---

## Format: linear vs decibel

Keep the internal format in `'LINEAR'` during pipeline processing — SAR indices (`RVI`, `DPSVI`, `VV_VH_ratio`, etc.) must be computed on linear power, not on logarithms. `NdviSeasonality` enforces this automatically.

Convert to decibel at the end if you want to visualise VV or VH directly:

```python
processor = S1ARDProcessor(format='DB')  # applies 10·log10 at the end
```

`format='DB'` compresses dynamic range for display but is **not interchangeable** with linear for downstream arithmetic.

---

## Tips and caveats

### Don't skip terrain correction in mountainous ROIs

Even in gently sloping terrain, RTC improves temporal stability of backscatter. Skip it only for flat coastal/agricultural areas where the computation cost matters and the DEM adds no information.

### Match the kernel to your target size

A `speckle_filter_kernel_size=7` (default) is a 70 m window at S1 native resolution. For very fine-grained work (e.g. small reservoirs, individual field boundaries), drop to 5 or 3 to preserve detail. For coarse monitoring (biomass over forest stands), 9 is reasonable.

### ENL assumption

All Lee-family and Gamma MAP filters assume ENL = 4.4, which is the standard value for S1 IW GRD. If you apply ARD to other SAR modes (EW, SM) you should recompute ENL — but `S1ARDProcessor` is targeted at IW, which is what 99% of ndvi2gif users will use.

### Cloud filter has no effect

`cloud_filter` and `max_cloud_cover` are ignored for `sat='S1'` — SAR is not affected by clouds. They are kept in the API only for interface consistency with optical sensors.

### Computational cost

`use_sar_ard=True` adds meaningful server-side work — terrain correction involves a DEM join plus per-pixel trigonometry, and Refined Lee runs four directional convolutions. For large ROIs or long time spans, expect `getInfo()` and map rendering to be noticeably slower. Export to Drive/Asset instead of pulling via `getInfo()`.

### Debugging

If output backscatter looks off, disable one stage at a time:

```python
# Stage 1: raw + basic filter
ns = NdviSeasonality(..., sat='S1', use_sar_ard=False)

# Stage 2: ARD without terrain correction
ns = NdviSeasonality(..., sat='S1', use_sar_ard=True,
                     sar_terrain_correction=False)

# Stage 3: full ARD
ns = NdviSeasonality(..., sat='S1', use_sar_ard=True,
                     sar_terrain_correction=True)
```

This isolates whether a problem comes from the DEM/angle geometry or from the speckle filter.

---

## API reference

### `S1ARDProcessor(...)`

| Parameter | Type | Default | Description |
|---|---|---|---|
| `speckle_filter` | str or None | `'REFINED_LEE'` | `'BOXCAR'`, `'LEE'`, `'REFINED_LEE'`, `'GAMMA_MAP'`, `'LEE_SIGMA'`, or `None` |
| `speckle_filter_kernel_size` | int (odd) | `7` | Filter kernel size in pixels |
| `terrain_correction` | bool | `True` | Apply radiometric terrain correction |
| `terrain_flattening_model` | str | `'VOLUME'` | `'VOLUME'` (vegetation) or `'SURFACE'` (bare soil/water) |
| `dem` | str | `'COPERNICUS_30'` | `'COPERNICUS_30'`, `'COPERNICUS_90'`, `'SRTM_30'`, `'SRTM_90'` |
| `format` | str | `'LINEAR'` | `'LINEAR'` (for indices) or `'DB'` (for visualisation) |

### Methods

| Method | Returns | Description |
|---|---|---|
| `process_image(image)` | `ee.Image` | Full pipeline: terrain correction → speckle filter → format conversion |
| `apply_terrain_correction(image)` | `ee.Image` | Only the RTC step |
| `apply_speckle_filter(image)` | `ee.Image` | Only the speckle step (routes to the selected filter) |
| `to_db(image)` | `ee.Image` | Convert VV and VH from linear to decibel |

### `NdviSeasonality` SAR-specific parameters

| Parameter | Default | Notes |
|---|---|---|
| `sat='S1'` | — | Select Sentinel-1 as source |
| `orbit` | `'BOTH'` | `'ASCENDING'`, `'DESCENDING'`, or `'BOTH'` |
| `use_sar_ard` | `True` | Enable the `S1ARDProcessor` pipeline |
| `sar_speckle_filter` | `'REFINED_LEE'` | Passed to the processor |
| `sar_terrain_correction` | `True` | Passed to the processor |
| `sar_terrain_model` | `'VOLUME'` | Passed to the processor |
| `normalize_sar` | `False` | Rescale SAR indices to `[0, 1]` |

---

## References

Vollrath, A., Mullissa, A., Reiche, J. (2020). Angular-based radiometric slope correction for Sentinel-1 on Google Earth Engine. *Remote Sensing*, 12(11), 1867.

Lee, J.S. (1980). Digital image enhancement and noise filtering by use of local statistics. *IEEE Transactions on Pattern Analysis and Machine Intelligence*, (2), 165–168.

Lee, J.S. (1983). Digital image smoothing and the sigma filter. *Computer Vision, Graphics, and Image Processing*, 24(2), 255–269.

Lopes, A., Nezry, E., Touzi, R., Laur, H. (1993). Structure detection and statistical adaptive speckle filtering in SAR images. *International Journal of Remote Sensing*, 14(9), 1735–1758.

Kim, Y., Jackson, T., Bindlish, R., Lee, H., Hong, S. (2012). Radar vegetation index for estimating the vegetation water content of rice and soybean. *IEEE Geoscience and Remote Sensing Letters*, 9(4), 564–568.

Mandal, D., Kumar, V., Ratha, D., Dey, S., Bhattacharya, A., Lopez-Sanchez, J.M., McNairn, H., Rao, Y.S. (2020). Dual polarimetric radar vegetation index for crop growth monitoring using Sentinel-1 SAR data. *Remote Sensing of Environment*, 247, 111954.

Periasamy, S. (2018). Significance of dual polarimetric synthetic aperture radar in biomass retrieval: An attempt on Sentinel-1. *Remote Sensing of Environment*, 217, 537–549.
