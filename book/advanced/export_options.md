# Export & Batch Processing

`NdviSeasonality` offers several export paths, each aimed at a different scale of analysis. The guidance below is opinionated: pick the smallest tool that still fits your job.

| You want… | Use this |
|---|---|
| A single GeoTIFF on your laptop, ROI < ~2 500 km² | `get_export_single()` |
| One GeoTIFF per year on your laptop | `get_export()` |
| A large file you can leave running | `export_to_drive()` |
| A reusable layer other GEE scripts can read | `export_to_asset()` |
| An ROI too big for any single export | `export_with_fishnet()` |
| An animated preview | `get_gif()` |

All of these sit on the same Earth Engine export primitives (`Export.image.toDrive`, `Export.image.toAsset`). The wrappers exist to fill in ROI, scale, and sensible defaults from the `NdviSeasonality` configuration.

---

## 1. Local exports through `geemap`

Local exports stream the image from GEE to your disk as a single GeoTIFF. They are convenient but have hard limits: roughly 2 GB per file and a few minutes of server time. If you routinely hit those limits, switch to Drive or Asset exports.

### `get_export_single()` — one image

```python
image = ns.get_year_composite().first()
ns.get_export_single(image, name='ndvi_2022.tif', scale=10, crs='EPSG:4326')
```

| Argument | Default | Description |
|---|---|---|
| `image` | — | Any `ee.Image`. Usually a composite from `get_year_composite()`. |
| `name` | `'mycomposition.tif'` | File name written to the current working directory. |
| `crs` | `'EPSG:4326'` | Output projection. Use a local UTM zone if you care about metric measurements on the result. |
| `scale` | `10` | Pixel size in metres. **Match the sensor**, not a number you like: 10 for S2/S1, 30 for Landsat, 250–500 for MODIS, 300 for S3. |

### `get_export()` — one GeoTIFF per year

```python
ns.get_export(scale=10, crs='EPSG:4326')
# ndvi_max_2020.tif, ndvi_max_2021.tif, ndvi_max_2022.tif, ...
```

Filenames follow the pattern `{sat}_{index}_{stat}_{year}.tif`, where `{stat}` is `max`, `mean`, `median`, or `p{N}` when `key='percentile'`. One multi-band file per year; bands are the periods (`winter, spring, summer, autumn` for `periods=4`; `january…december` for `periods=12`; `p1…pN` otherwise).

> **Don't oversample.** Setting `scale=10` on Landsat or MODIS data does not produce 10 m imagery — it produces a resampled fake that wastes quota. Use the native resolution of the sensor.

---

## 2. Drive exports — `export_to_drive()`

For anything that won't fit comfortably into a local download, go through Drive. This starts an Earth Engine batch task, returns immediately, and writes the result to your Google Drive when the server is ready.

```python
task = ns.export_to_drive(
    image=classified,
    description='landcover_2022_RF',
    folder='ndvi2gif_exports',
    scale=10,
    crs='EPSG:4326',
    file_format='GeoTIFF',
    format_options={'cloudOptimized': True, 'compression': 'LZW'},
    max_pixels=int(1e13),
    file_dimensions=10000,
    clip_region=True,
)
```

| Argument | Purpose |
|---|---|
| `description` | Task name shown in the Tasks panel. Also used as the base filename. |
| `region` | Defaults to `self.roi`. Pass an `ee.Geometry` for a sub-region. |
| `scale` | Defaults to the sensor's native resolution via `_default_scale_for_sat()`. |
| `folder` | Drive subfolder. If `None`, the file lands in the Drive root. |
| `format_options` | `{'cloudOptimized': True}` for COGs, `{'compression': 'LZW'}` for smaller files. |
| `file_dimensions` | Tile output into pieces of this pixel-width. Use `8192` or `10000` for massive rasters — otherwise GEE will split them arbitrarily or fail. |
| `clip_region` | Clips the image to `region` before exporting. Recommended — faster and smaller. |

The returned `ee.batch.Task` lets you poll progress:

```python
task.status()
# {'state': 'RUNNING', 'description': '...', ...}
```

Or just watch the GEE Tasks panel at https://code.earthengine.google.com/tasks.

> **Cloud-Optimized GeoTIFFs are worth the flag.** Add `format_options={'cloudOptimized': True}` once and downstream tools (QGIS, rasterio, Cloud Storage) will read your outputs faster without loading the whole file. There's no downside for consumers that don't support COG.

---

## 3. Asset exports — `export_to_asset()`

Use this when the output will be consumed by *other GEE scripts* — e.g. a monthly classification feeding a trend analysis, or a mask you want to re-apply across notebooks.

```python
task = ns.export_to_asset(
    image=classified,
    asset_id='users/yourname/ndvi2gif/landcover_2022',
    pyramiding_policy={'classification': 'mode'},
    scale=10,
    overwrite=True,
)
```

| Argument | Notes |
|---|---|
| `asset_id` | Full asset path. Parent folders must exist. |
| `pyramiding_policy` | **Critical for class maps.** Default is `mean`, which averages class IDs across zoom levels — meaningless. Pass `{'classification': 'mode'}` for categorical data, or `{'ndvi': 'mean'}` for continuous. |
| `crs_transform` | Use this *instead of* `scale` when you need to align to an existing grid (e.g. reproducing Landsat pixel edges). Mutually exclusive with `scale`. |
| `overwrite` | Deletes any existing asset at `asset_id` before re-exporting. |

Assets live on Google's servers and are near-instant to load in other scripts via `ee.Image('users/yourname/...')`. There's no download step.

> **Asset vs Drive, quick rule.** If the next consumer is Python/QGIS, export to Drive. If the next consumer is another GEE script, export to Asset. Doing both is sometimes right; doing neither and relying on `get_export_single()` beyond ~1 GB is almost never.

---

## 4. Tiled exports — `export_with_fishnet()`

Some ROIs are simply too large for any single export, even on Drive — especially at S2's 10 m resolution. `export_with_fishnet()` divides the ROI into a grid and exports each tile separately:

```python
ns.export_with_fishnet(
    image=composite,
    name_prefix='ndvi_max_2022',
    scale=10,
    crs='EPSG:4326',
)
# ndvi_max_2022_tile_0.tif, ndvi_max_2022_tile_1.tif, ...
```

Tile size adapts to the resolution: 50 km × 50 km at ≥ 30 m, 25 km × 25 km at finer scales. Only tiles that intersect the ROI are written; the rest are skipped.

This runs through `geemap.ee_export_image` per tile, so it streams to your local disk. For very large jobs (e.g. national-scale Landsat mosaics) prefer `export_to_drive()` with `file_dimensions` — it lets GEE orchestrate the tiling server-side and keeps the task in the batch queue.

---

## 5. `get_gif()` — previews, not publications

```python
ns.get_gif(name='seasonal_evolution.gif', bands=['winter', 'spring', 'summer'])
```

Renders one frame per year, using the three selected periods as R/G/B. Two files are written:

- `seasonal_evolution.gif` — raw animation.
- `seasonal_evolution_texted.gif` — the same with year labels burned in.

Visualisation defaults to sensor-appropriate stretches (`-25, 0` dB for S1; `0.15, 0.85` for optical indices). Fine for README demos and quick presentations, not for quantitative figures.

> **Pick RGB periods that actually differ.** Using `january, february, march` on a monthly composite gives you a grey animation — the three months look nearly identical. `winter, summer, autumn` (or months spread across the year) encode the seasonal signal visibly.

---

## 6. Deciding the right route

A rough decision tree:

```text
Output size < 1 GB and you want it on disk now
    └─► get_export_single()  or  get_export()  for multi-year

Output is large OR you don't want to babysit it
    └─► export_to_drive()  (add cloudOptimized=True)

Next consumer is another GEE script
    └─► export_to_asset()  (watch pyramidingPolicy for class maps)

ROI is so large it breaks even Drive exports
    └─► export_with_fishnet()  or  export_to_drive(file_dimensions=10000)

You need a visual preview, not data
    └─► get_gif()
```

---

## Tips and caveats

### Always clip

For every method that accepts `clip_region`, leave it on (the default). Without clipping, GEE processes pixels outside your ROI that are then discarded — pure wasted compute and quota.

### `maxPixels` is a ceiling, not a target

The default `1e13` is generous. If you're hitting it, the correct response is usually "reduce scale" or "tile the export", not "raise the ceiling".

### CRS choice

`EPSG:4326` is the universal default and fine for visual products. For measurements (area, distance, slope), export in the appropriate UTM zone or another equal-area projection. QGIS reprojection is non-lossless for rasters; do it once at export time if you can.

### Parallel tasks

Drive and Asset tasks run on the GEE backend, so you can queue many at once. A typical pattern is a for-loop over years calling `export_to_drive()` with unique descriptions — GEE will schedule them in parallel up to your quota.

### Don't commit exported GeoTIFFs

They're large and derivable from the code. Either put them on cloud storage or in `.gitignore`.

---

## API reference

| Method | Returns | Purpose |
|---|---|---|
| `get_export_single(image, name, crs, scale)` | `None` | Local GeoTIFF via `geemap.ee_export_image` |
| `get_export(crs, scale)` | `None` | One local GeoTIFF per year, auto-named |
| `export_to_drive(image, description, region, scale, crs, folder, file_format, format_options, max_pixels, file_dimensions, clip_region)` | `ee.batch.Task` | Batch export to Google Drive |
| `export_to_asset(image, asset_id, description, region, scale, crs, crs_transform, pyramiding_policy, max_pixels, overwrite, clip_region)` | `ee.batch.Task` | Batch export to an Earth Engine asset |
| `export_with_fishnet(image, name_prefix, scale, crs)` | `None` | Tiled local export for very large ROIs |
| `get_gif(name, bands)` | `None` | Animated GIF preview of the time series |

---

## References

Gorelick, N., Hancher, M., Dixon, M., Ilyushchenko, S., Thau, D., Moore, R. (2017). Google Earth Engine: Planetary-scale geospatial analysis for everyone. *Remote Sensing of Environment*, 202, 18–27.

Earth Engine Export Guide. Google Developers. https://developers.google.com/earth-engine/guides/exporting

Cloud Optimized GeoTIFF. https://www.cogeo.org/
