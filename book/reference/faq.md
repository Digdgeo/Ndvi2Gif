# Frequently Asked Questions

## Installation & Setup

### Q: Which Python version should I use?

**A:** Python 3.8 to 3.12 are supported. We recommend Python 3.11 for best performance.

### Q: Can I use Ndvi2Gif in Google Colab?

**A:** Yes! Install with:
```ipython3
!pip install ndvi2gif
import ee
ee.Authenticate()
ee.Initialize(project='your-project-id')
```

### Q: Do I need a Google Earth Engine account?

**A:** Yes. Sign up at https://signup.earthengine.google.com/ - it's free for research and education.

## Usage Questions

### Q: What's the difference between `periods=4` and `periods=12`?

**A:**
- `periods=4`: Seasonal composites (Winter, Spring, Summer, Fall)
- `periods=12`: Monthly composites (Jan, Feb, ..., Dec)
- `periods=24`: Bi-monthly (every ~15 days)

### Q: Should I use 'median' or 'max' for my analysis?

**A:**
- **median**: Best for cloud-prone areas, more robust
- **max**: Good for vegetation indices (captures peak greenness)
- **percentile**: Flexible (try 75th, 85th, 95th)

For optical sensors (S2, Landsat), start with `median`. For SAR (S1), use `mean`.

### Q: How do I map variability instead of the typical value?

**A:** Use one of the dispersion reducers: `key='std'`, `'variance'`, `'range'`
(max − min) or `'cv'` (std / mean). They aggregate the same observations, but
answer "how much did this pixel move inside the period" rather than "what was
its usual value" — handy for phenological change, disturbances and flooded
areas. `'cv'` is only meaningful for indices that stay positive, so avoid it for
Sentinel-1 backscatter in dB.

### Q: How do I choose the right satellite?

**A:**

| Satellite | Resolution | Revisit | Best For |
|-----------|------------|---------|----------|
| Sentinel-2 | 10m | 5 days | Detailed agriculture, high-res |
| Sentinel-1 | 10m | 6-12 days | All-weather (SAR), structure |
| Landsat 8/9 | 30m | 16 days | Historical analysis, large areas |
| MODIS | 500m | Daily | Global monitoring, daily data |

### Q: Can I analyze my entire country/continent?

**A:** Yes, but use care:
- Large areas require Earth Engine task exports (use `export_to_drive()`)
- Consider processing by tiles or regions
- Use coarser resolution (MODIS) for continental scales

### Q: How much time do I need in `start_year` to `end_year`?

**A:** Minimum:
- **Seasonal analysis**: 1 year minimum
- **Monthly analysis**: 1 year minimum
- **Cloud-prone areas**: 2+ years recommended
- **Trend analysis**: 3-5+ years

## Technical Issues

### Q: My analysis is very slow or timing out

**A:** Try:
1. Reduce ROI size: `roi = roi.buffer(-100)` to shrink slightly
2. Increase scale: `scale=30` instead of `scale=10`
3. Shorten time period: Process one year at a time
4. Use Earth Engine tasks instead of direct computation
5. Simplify geometry: Complex polygons slow processing

### Q: I get "User memory limit exceeded" errors

**A:** Your ROI or time span is too large. Solutions:
```python
# Use export instead of direct computation
ndvi.export_to_drive(
    image=composite,
    description='my_analysis',
    maxPixels=1e13,
    scale=30  # Increase scale
)

# Or split into tiles
# Process smaller regions and mosaic later
```

### Q: Clouds are ruining my results

**A:**
1. Cloud filtering is on by default; lower the scene threshold, e.g. `max_cloud_cover=10`
2. Use `median` or percentiles instead of `max`
3. Extend time period to get more images
4. Try higher percentiles: `percentile=85` or `percentile=90`

### Q: I get all zeros or all NaN

**A:** Common causes:
1. **Wrong ROI projection**: Convert to EPSG:4326
2. **ROI outside sensor coverage**: Check if area has data
3. **Wrong date range**: S3 starts 2016, S2 starts 2015, etc.
4. **Cloud masking too aggressive**: Try `cloud_filter=False` or `scl_mask=False`

### Q: How do I handle incomplete years (e.g., monitoring 2024)?

**A:** Ndvi2Gif automatically handles incomplete years:
```python
# Will use all available data in 2024
ndvi = NdviSeasonality(
    roi=roi,
    sat='S2',
    periods=12,
    start_year=2024,
    end_year=2024  # Incomplete year OK!
)
```

## Data & Indices

### Q: Which index should I use for my application?

**A:**
- **Vegetation health**: NDVI, EVI
- **Crop monitoring**: NDVI, EVI, NDRE (Sentinel-2 only)
- **Water detection**: NDWI, MNDWI
- **Drought**: NDMI, MSI
- **Water quality**: Turbidity, NDCI (Sentinel-3)
- **Burn scars**: NBRI
- **Harvest detection**: RVI, VV/VH (Sentinel-1 SAR)

See [Indices Reference](indices.md) for complete list.

### Q: Can I use multiple indices together?

**A:** Yes! Create separate processors:
```python
# NDVI processor
ndvi = NdviSeasonality(roi=roi, sat='S2', index='ndvi', ...)
ndvi_comp = ndvi.get_year_composite()

# NDWI processor
ndwi = NdviSeasonality(roi=roi, sat='S2', index='ndwi', ...)
ndwi_comp = ndwi.get_year_composite()

# Combine in classification or analysis
```

### Q: Why are Sentinel-2 and Landsat NDVI values different?

**A:** Normal! Different sensors have:
- Different band wavelengths
- Different spatial resolution
- Different atmospheric correction

Values may differ by 0.05-0.1. Focus on relative changes, not absolute values.

### Q: Can I access raw bands instead of indices?

**A:** Yes, use Google Earth Engine directly:
```python
import ee

collection = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED') \
    .filterDate('2023-01-01', '2023-12-31') \
    .filterBounds(roi)

# Select specific bands
bands = collection.select(['B4', 'B8'])  # Red, NIR
```

## Export & Output

### Q: How do I export to GeoTIFF instead of GIF?

**A:**
```python
ndvi.get_export(
    year=2023,
    filename='output.tif'
)
```

### Q: Exports are huge! How do I reduce file size?

**A:**
1. Increase scale: `scale=30` instead of `scale=10`
2. Reduce bands: Export single year instead of multi-year
3. Crop extent: Use smaller ROI
4. Compress: Use DEFLATE compression in GeoTIFF

### Q: Can I export to shapefile/vector?

**A:** For zonal statistics, yes:
```python
# Sample composite at points
points = ee.FeatureCollection('path/to/points')
composite = ndvi.get_year_composite().first()

samples = composite.sampleRegions(
    collection=points,
    scale=10
)

# Export to Drive
task = ee.batch.Export.table.toDrive(
    collection=samples,
    description='ndvi_samples',
    fileFormat='SHP'
)
task.start()
```

### Q: How do I open exported GeoTIFFs?

**A:**
```python
import rasterio
import numpy as np

# Open file
with rasterio.open('ndvi_2023.tif') as src:
    data = src.read()  # All bands
    meta = src.meta

# Or specific band
with rasterio.open('ndvi_2023.tif') as src:
    jan = src.read(1)  # First band (January)
```

## Advanced Topics

### Q: Can I use Ndvi2Gif with my own ImageCollection?

**A:** Yes, but requires modification:
```python
# Create custom collection
custom_col = ee.ImageCollection('YOUR/COLLECTION') \
    .filterDate('2023-01-01', '2023-12-31') \
    .filterBounds(roi)

# You'll need to extend NdviSeasonality class
# or process directly with Earth Engine
```

### Q: How do I add a new index?

**A:** See [Contributing Guide](contributing.md). Basic steps:
1. Add to appropriate index set
2. Implement computation function
3. Register in dispatch dictionary
4. Map to compatible sensors

### Q: Can I parallelize processing for multiple ROIs?

**A:** Yes, using Earth Engine tasks:
```python
# Submit multiple tasks
for roi in roi_list:
    ndvi = NdviSeasonality(roi=roi, sat='S2', ...)
    composite = ndvi.get_year_composite().first()

    # Export (runs in background)
    ndvi.export_to_drive(
        image=composite,
        description=f'ndvi_{roi.id}',
        folder='batch_exports'
    )
```

### Q: How accurate is the phenology analysis?

**A:** Depends on:
- Temporal resolution (more images = better)
- Cloud cover (less clouds = better)
- Index choice (NDVI, EVI work well)
- Smoothing parameters

Typical accuracy: ±7-14 days for SOS/EOS with good data.

## Getting Help

### Q: Where can I get more help?

**A:**
1. Check this FAQ
2. Review [tutorials](../tutorials/basic_ndvi.md)
3. Search [GitHub Issues](https://github.com/Digdgeo/Ndvi2Gif/issues)
4. Ask in [GitHub Discussions](https://github.com/Digdgeo/Ndvi2Gif/discussions)
5. Check [Earth Engine Forum](https://groups.google.com/g/google-earth-engine-developers)

### Q: How do I report a bug?

**A:** Create an issue at https://github.com/Digdgeo/Ndvi2Gif/issues with:
- Minimal code example
- Error message
- Python version, ndvi2gif version
- Operating system

### Q: Can I contribute?

**A:** Yes! See [Contributing Guide](contributing.md). We welcome:
- New indices
- Bug fixes
- Documentation improvements
- Example notebooks
- Dataset support
