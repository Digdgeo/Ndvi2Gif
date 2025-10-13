# Multi-Sensor Comparison

Compare vegetation indices across different satellite sensors.

## Overview

Different sensors have different characteristics:
- **Spatial resolution**: Detail level
- **Temporal resolution**: Revisit frequency  
- **Spectral bands**: Available wavelengths
- **Data availability**: Historical coverage

## Quick Comparison

| Sensor | Resolution | Revisit | Start Date | Free |
|--------|------------|---------|------------|------|
| Sentinel-2 | 10m | 5 days | 2015 | Yes |
| Sentinel-1 | 10m | 6-12 days | 2014 | Yes |
| Landsat 8/9 | 30m | 16 days | 2013 | Yes |
| MODIS | 500m | Daily | 2000 | Yes |

## Example: Comparing S2 and Landsat

```python
import ee
from ndvi2gif import NdviSeasonality

ee.Initialize()

roi = ee.Geometry.Point([-3.7, 40.4]).buffer(5000)

# Sentinel-2
s2_ndvi = NdviSeasonality(
    roi=roi,
    sat='S2',
    periods=12,
    start_year=2023,
    end_year=2024,
    index='ndvi',
    key='median'
)

# Landsat 8
l8_ndvi = NdviSeasonality(
    roi=roi,
    sat='L8',
    periods=12,
    start_year=2023,
    end_year=2024,
    index='ndvi',
    key='median'
)

# Get composites
s2_composite = s2_ndvi.get_year_composite(year=2023)
l8_composite = l8_ndvi.get_year_composite(year=2023)
```

*Complete tutorial coming soon...*
