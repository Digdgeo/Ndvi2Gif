# Quick Start Guide

Let's create your first multi-seasonal NDVI analysis with Ndvi2Gif!

## Your First NDVI Analysis

This example creates a seasonal NDVI composite for a region using Sentinel-2 data.

### Step 1: Import and Initialize

```python
import ee
import geemap
from ndvi2gif import NdviSeasonality

# Initialize Earth Engine
ee.Initialize()
```

### Step 2: Define Your Region of Interest (ROI)

There are multiple ways to define your ROI:

#### Option A: Draw on a Map

```python
# Create an interactive map
Map = geemap.Map()
Map.centerObject(ee.Geometry.Point([-3.7, 40.4]), 10)  # Madrid, Spain
Map

# Draw a polygon on the map, then:
roi = Map.user_roi
```

#### Option B: Use Coordinates

```python
# Define a bounding box
roi = ee.Geometry.Rectangle([-3.8, 40.3, -3.6, 40.5])

# Or a point with buffer
roi = ee.Geometry.Point([-3.7, 40.4]).buffer(5000)  # 5km radius
```

#### Option C: Use a Shapefile

```python
# Load from a shapefile
roi = 'path/to/your/area.shp'
```

### Step 3: Create an NDVI Composite

```python
# Initialize the processor
ndvi_processor = NdviSeasonality(
    roi=roi,
    sat='S2',              # Sentinel-2
    periods=12,            # Monthly composites
    start_year=2023,
    end_year=2024,
    index='ndvi',          # Normalized Difference Vegetation Index
    key='median'           # Use median for cloud-free composites
)
```

### Step 4: Get the Composite

```python
# get_year_composite() returns an ee.ImageCollection with one image per year.
# Each image has one band per period (january, february, ..., december).
all_composites = ndvi_processor.get_year_composite()

# Select a single year's composite (here, the first year in the range)
composite = all_composites.first()
```

### Step 5: Visualize

```python
# Visualize on a map
Map = geemap.Map()
vis_params = {
    'min': 0,
    'max': 1,
    'palette': ['red', 'yellow', 'green']
}
# A palette can only be applied to a single band, so select one period (e.g. june)
Map.addLayer(composite.select('june'), vis_params, 'NDVI June')
Map.centerObject(roi, 10)
Map
```

### Step 6: Create an Animated GIF

```python
# Generate an RGB GIF: one frame per year, using three periods as R, G, B.
# The frame rate (10 fps) and dimensions are set internally.
ndvi_processor.get_gif(
    name='ndvi.gif',
    bands=['march', 'june', 'september']   # [R, G, B] periods
)
```

## Complete Example

Here's a complete working example:

```python
import ee
import geemap
from ndvi2gif import NdviSeasonality

# Initialize
ee.Initialize()

# Define ROI (Madrid, Spain)
roi = ee.Geometry.Rectangle([-3.8, 40.3, -3.6, 40.5])

# Create processor
ndvi = NdviSeasonality(
    roi=roi,
    sat='S2',
    periods=12,
    start_year=2023,
    end_year=2024,
    index='ndvi',
    key='median'
)

# Get the yearly composites (ee.ImageCollection, one image per year)
composites = ndvi.get_year_composite()

# Create an RGB GIF (one frame per year)
ndvi.get_gif(name='madrid_ndvi.gif', bands=['march', 'june', 'september'])

print("✓ Analysis complete! Check your GIF file.")
```

## Understanding the Parameters

### `sat` - Satellite Sensor

- `'S2'` - Sentinel-2 (10m, 5-day revisit)
- `'S1'` - Sentinel-1 SAR (10m, 6-12 day revisit)
- `'Landsat'` - Landsat 4-9 merged (30m, 16-day revisit)
- `'MODIS'` - MODIS (500m, daily)
- `'S3'` - Sentinel-3 OLCI (300m)
- `'ERA5'` - ERA5-Land climate reanalysis (~11km, daily)
- `'CHIRPS'` - CHIRPS precipitation (~5.5km, daily)

### `periods` - Temporal Resolution

- `4` - Seasonal (Winter, Spring, Summer, Fall)
- `12` - Monthly
- `24` - Bi-monthly
- Custom numbers for specific applications

### `key` - Statistical Method

- `'median'` - Robust to outliers (recommended for optical)
- `'max'` - Maximum value (good for vegetation indices)
- `'mean'` - Average value
- `'percentile'` - Custom percentile (specify with `percentile` parameter)

### `index` - Spectral Index

- `'ndvi'` - Vegetation health
- `'evi'` - Enhanced vegetation index
- `'ndwi'` - Water content
- `'ndmi'` - Moisture stress
- Many more! See [Indices Reference](../reference/indices.md)

## Export Options

### Export as GeoTIFF

```python
# Export every year's composite as a GeoTIFF to the current directory.
# Filenames are generated automatically as {sat}_{index}_{key}_{year}.tif
ndvi.get_export(scale=10)
```

### Export to Google Drive

```python
# Export directly to Google Drive
ndvi.export_to_drive(
    image=composite,
    description='ndvi_madrid_2023',
    folder='earthengine_exports',
    scale=10
)
```

## Common Patterns

### Multi-Year Analysis

```python
# Analyze multiple years
ndvi = NdviSeasonality(
    roi=roi,
    sat='S2',
    periods=12,
    start_year=2020,
    end_year=2024,
    index='ndvi'
)

# Get the yearly composites (one image per year)
composites = ndvi.get_year_composite()

# Thanks to the band-based design, cross-year seasonal queries are one line.
# For example, the maximum June NDVI across all years:
june_max = composites.select('june').max()
```

### High Temporal Resolution

```python
# Bi-monthly for detailed monitoring
ndvi = NdviSeasonality(
    roi=roi,
    sat='S2',
    periods=24,  # Every 15 days
    start_year=2024,
    end_year=2024,
    index='ndvi'
)
```

### Using Percentiles

```python
# Use 85th percentile for cloud contamination
ndvi = NdviSeasonality(
    roi=roi,
    sat='S2',
    periods=12,
    start_year=2023,
    end_year=2024,
    index='ndvi',
    key='percentile',
    percentile=85
)
```

## Next Steps

Now that you understand the basics, explore:

- [ROI Options](roi_options.md) - Learn all the ways to define your study area
- [Indices Reference](../reference/indices.md) - Explore 40+ available indices

## Tips & Tricks

1. **Start small**: Test with small ROIs and short time periods first
2. **Use median for optical**: It's more robust to clouds
3. **Percentiles are flexible**: Try 75th, 85th, or 95th percentiles
4. **SAR uses mean**: For Sentinel-1, use `key='mean'`
5. **Check your projection**: ROIs work best in EPSG:4326 (WGS84)
