# Basic NDVI Analysis

A comprehensive guide to vegetation monitoring with NDVI (Normalized Difference Vegetation Index).

## What is NDVI?

NDVI measures vegetation health by comparing red and near-infrared (NIR) light reflected by plants:

$$
NDVI = \frac{NIR - Red}{NIR + Red}
$$

- **Values**: -1 to +1
- **Interpretation**:
  - < 0: Water, clouds, snow
  - 0 - 0.2: Bare soil, sand
  - 0.2 - 0.5: Sparse vegetation, grassland
  - 0.5 - 0.8: Dense vegetation, healthy crops
  - > 0.8: Very dense vegetation, forests

## Basic NDVI Workflow

### Step 1: Setup

```python
import ee
import geemap
from ndvi2gif import NdviSeasonality
import matplotlib.pyplot as plt

# Initialize
ee.Initialize()
```

### Step 2: Define Study Area

```python
# Example: Agricultural area in Spain
roi = ee.Geometry.Rectangle([-6.2, 37.8, -6.0, 38.0])

# Visualize
Map = geemap.Map()
Map.addLayer(roi, {}, 'Study Area')
Map.centerObject(roi, 11)
Map
```

### Step 3: Create NDVI Processor

```python
ndvi = NdviSeasonality(
    roi=roi,
    sat='S2',              # Sentinel-2 (10m resolution)
    periods=12,            # Monthly composites
    start_year=2023,
    end_year=2024,
    index='ndvi',
    key='median'           # Robust to clouds (cloud filtering is on by default)
)
```

### Step 4: Generate Composite

```python
# get_year_composite() returns an ImageCollection (one image per year)
composites = ndvi.get_year_composite()
composite_2023 = composites.first()

# Check band names
print("Bands:", composite_2023.bandNames().getInfo())
# Output: ['january', 'february', 'march', 'april', 'may', 'june',
#          'july', 'august', 'september', 'october', 'november', 'december']
```

### Step 5: Visualize Results

```python
# Visualization parameters
vis_params = {
    'bands': ['july', 'april', 'january'],  # RGB = Summer, Spring, Winter
    'min': 0,
    'max': 0.8
}

# Add to map (a 3-band RGB composite does not take a palette)
Map = geemap.Map()
Map.addLayer(composite_2023, vis_params, 'NDVI 2023')
Map.addLayer(roi, {}, 'ROI')
Map.centerObject(roi, 11)
Map
```

### Step 6: Create Animated GIF

```python
# Generate temporal animation (one frame per year, RGB from three periods)
ndvi.get_gif(
    name='ndvi_monthly.gif',
    bands=['july', 'april', 'january']   # [R, G, B] periods
)
```

## Advanced NDVI Analysis

### Multi-Year Comparison

Compare vegetation patterns across years:

```python
# Process multiple years
ndvi_multi = NdviSeasonality(
    roi=roi,
    sat='S2',
    periods=12,
    start_year=2020,
    end_year=2023,
    index='ndvi',
    key='median'
)

# One call returns all years (2020–2023) as an ImageCollection
composites = ndvi_multi.get_year_composite()
print("Years processed:", composites.size().getInfo())
```

### Seasonal Analysis

Focus on specific seasons:

```python
# Quarterly analysis
ndvi_seasonal = NdviSeasonality(
    roi=roi,
    sat='S2',
    periods=4,  # Winter, Spring, Summer, Fall
    start_year=2023,
    end_year=2024,
    index='ndvi'
)

composite = ndvi_seasonal.get_year_composite().first()

# Visualize seasons (period names for periods=4: winter, spring, summer, autumn)
Map = geemap.Map()
Map.addLayer(
    composite.select('winter'),
    {'min': 0, 'max': 0.8, 'palette': ['red', 'yellow', 'green']},
    'Winter'
)
Map.addLayer(
    composite.select('summer'),
    {'min': 0, 'max': 0.8, 'palette': ['red', 'yellow', 'green']},
    'Summer'
)
Map.centerObject(roi, 11)
Map
```

### High Temporal Resolution

Bi-monthly for detailed crop monitoring:

```python
ndvi_detailed = NdviSeasonality(
    roi=roi,
    sat='S2',
    periods=24,  # Every ~15 days
    start_year=2024,
    end_year=2024,
    index='ndvi',
    key='percentile',
    percentile=75
)
```

## Statistical Methods for NDVI

### Median (Recommended)

```python
# Best for cloud-prone areas
ndvi_median = NdviSeasonality(
    roi=roi,
    sat='S2',
    periods=12,
    start_year=2023,
    end_year=2024,
    index='ndvi',
    key='median'  # Robust to outliers
)
```

### Maximum

```python
# Maximum greenness
ndvi_max = NdviSeasonality(
    roi=roi,
    sat='S2',
    periods=12,
    start_year=2023,
    end_year=2024,
    index='ndvi',
    key='max'  # Peak values
)
```

### Percentile (Flexible)

```python
# 85th percentile - good balance
ndvi_p85 = NdviSeasonality(
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

## Extracting Statistics

### Zonal Statistics

Extract mean NDVI per month:

```python
# Get time series for ROI
composite = ndvi.get_year_composite().first()

# Extract values for each month
stats = []
for month in range(1, 13):
    band_name = composite.bandNames().get(month-1).getInfo()

    value = composite.select(band_name).reduceRegion(
        reducer=ee.Reducer.mean(),
        geometry=roi,
        scale=10,
        maxPixels=1e9
    ).getInfo()

    stats.append({
        'month': month,
        'ndvi': value[band_name]
    })

# Plot
import pandas as pd
df = pd.DataFrame(stats)
df.plot(x='month', y='ndvi', kind='line', marker='o')
plt.title('Mean NDVI by Month')
plt.xlabel('Month')
plt.ylabel('NDVI')
plt.grid(True)
plt.show()
```

### Point Sampling

Extract NDVI at specific locations:

```python
# Define sample points
points = ee.FeatureCollection([
    ee.Feature(ee.Geometry.Point([-6.1, 37.9]), {'name': 'Site A'}),
    ee.Feature(ee.Geometry.Point([-6.15, 37.95]), {'name': 'Site B'})
])

# Sample composite
composite = ndvi.get_year_composite().first()
samples = composite.sampleRegions(
    collection=points,
    scale=10
).getInfo()

print(samples['features'])
```

## Export Options

### Export as GeoTIFF

```python
# Exports each year's composite; filenames auto-generated as {sat}_{index}_{key}_{year}.tif
ndvi.get_export(scale=10)
```

### Export to Google Drive

```python
composite = ndvi.get_year_composite().first()

ndvi.export_to_drive(
    image=composite,
    description='ndvi_2023',
    folder='earthengine_exports',
    scale=10,
    crs='EPSG:4326'
)
```

## Interpretation Guide

### Agricultural Applications

- **Early Season (Mar-Apr)**: NDVI 0.2-0.4 indicates crop emergence
- **Peak Season (Jun-Jul)**: NDVI > 0.6 shows healthy, mature crops
- **Harvest (Sep-Oct)**: NDVI drops < 0.3 after harvest
- **Bare Soil (Nov-Feb)**: NDVI ~ 0.1-0.2

### Natural Vegetation

- **Evergreen Forests**: Consistent NDVI 0.7-0.9 year-round
- **Deciduous Forests**: NDVI varies 0.3 (winter) to 0.8 (summer)
- **Grasslands**: NDVI 0.3-0.6, responsive to rainfall
- **Wetlands**: Variable, influenced by water levels

## Common Issues & Solutions

### Issue: Low NDVI values everywhere

**Cause**: Cloud contamination or wrong bands

**Solution**:
```python
# Enable cloud masking
ndvi = NdviSeasonality(
    roi=roi,
    sat='S2',
    periods=12,
    start_year=2023,
    end_year=2024,
    index='ndvi',
    cloud_filter=True,     # Cloud filtering (on by default)
    max_cloud_cover=10     # ← Stricter scene cloud threshold
)
```

### Issue: Patchy results

**Cause**: Insufficient images

**Solution**: Use longer time periods or lower percentiles:
```python
# Extend time period
start_year=2022, end_year=2024  # More images

# Or use median instead of max
key='median'
```

### Issue: Negative NDVI on land

**Cause**: Water bodies or very bare soil

**Solution**: This is normal! Mask water if needed:
```python
# Mask water bodies
import ee

# Add water mask (example)
water = ee.Image('JRC/GSW1_3/GlobalSurfaceWater').select('occurrence')
composite_masked = composite.updateMask(water.lt(50))  # Keep only low water occurrence
```

## Next Steps

- [Indices Reference](../reference/indices.md) - Explore other vegetation indices
- [Time Series Analysis](../advanced/time_series.md) - Detect trends and phenology
- [Classification](../advanced/classification.md) - Use NDVI for land cover mapping
