# Region of Interest (ROI) Options

Ndvi2Gif provides multiple flexible ways to define your region of interest. This guide covers all available options.

## 1. Drawn Geometry (Interactive)

Draw directly on an interactive map using geemap:

```python
import geemap
import ee

# Create interactive map
Map = geemap.Map()
Map.centerObject(ee.Geometry.Point([-3.7, 40.4]), 10)
Map

# Draw a polygon using the drawing tools, then:
roi = Map.user_roi

# Use in analysis
from ndvi2gif import NdviSeasonality
ndvi = NdviSeasonality(roi=roi, sat='S2', ...)
```

```{tip}
Drawing is the easiest way to explore new areas interactively!
```

## 2. Shapefile

Load from a shapefile on disk:

```python
roi = '/path/to/your/study_area.shp'

ndvi = NdviSeasonality(
    roi=roi,
    sat='S2',
    periods=12,
    start_year=2023,
    end_year=2024
)
```

```{note}
Shapefiles should be in EPSG:4326 (WGS84) for best compatibility.
```

## 3. GeoJSON

Use GeoJSON files or strings:

```python
# From file
roi = 'path/to/area.geojson'

# Or directly
roi = {
    "type": "Polygon",
    "coordinates": [[
        [-3.8, 40.3],
        [-3.6, 40.3],
        [-3.6, 40.5],
        [-3.8, 40.5],
        [-3.8, 40.3]
    ]]
}
```

## 4. Earth Engine Geometry

Use Earth Engine geometry objects directly:

### Rectangle (Bounding Box)

```python
import ee

# Define by coordinates [west, south, east, north]
roi = ee.Geometry.Rectangle([-3.8, 40.3, -3.6, 40.5])
```

### Point with Buffer

```python
# Create circular ROI around a point
roi = ee.Geometry.Point([-3.7, 40.4]).buffer(5000)  # 5km radius
```

### Polygon

```python
# Custom polygon
roi = ee.Geometry.Polygon([[
    [-3.8, 40.3],
    [-3.6, 40.3],
    [-3.6, 40.5],
    [-3.8, 40.5],
    [-3.8, 40.3]
]])
```

## 5. eLTER Site (DEIMS ID)

Use ecological research site boundaries via DEIMS IDs:

```python
# Use eLTER site ID
roi = 'deimsid:ab8278e6-0b71-4b36-a6d2-e8f34aa3df30'

ndvi = NdviSeasonality(
    roi=roi,
    sat='S2',
    periods=12,
    start_year=2023,
    end_year=2024
)
```

Find DEIMS IDs at [DEIMS-SDR](https://deims.org/).

```{tip}
This is perfect for long-term ecological research sites!
```

## 6. Sentinel-2 Tile

Use Sentinel-2 MGRS tile codes:

```python
# Use tile code directly
roi = 'T30TVK'  # Madrid area

ndvi = NdviSeasonality(
    roi=roi,
    sat='S2',  # Works best with S2
    periods=12,
    start_year=2023,
    end_year=2024
)
```

Find tile codes at [Sentinel-2 Tiling Grid](https://maps.eatlas.org.au/index.html?intro=false&z=7&ll=146.90137,-19.07287&l0=ea_ref%3AWorld_ESA_Sentinel-2-tiling-grid_Poly,ea_ea-be%3AWorld_Bright-Earth-e-Atlas-basemap,google_HYBRID,google_TERRAIN,google_SATELLITE,google_ROADMAP&v0=,f,f,f,f,f).

## 7. Landsat Path/Row

Use Landsat WRS-2 path and row:

```python
# Use path/row notation
roi = '198/034'  # Path 198, Row 034

ndvi = NdviSeasonality(
    roi=roi,
    sat='L8',  # Works best with Landsat
    periods=12,
    start_year=2023,
    end_year=2024
)
```

Find path/row at [Landsat WRS-2 Grid](https://landsat.usgs.gov/landsat_acq).

## ROI Comparison Table

| Method | Best For | Requires | Projection |
|--------|----------|----------|------------|
| Drawn Geometry | Quick exploration | geemap, Jupyter | Any |
| Shapefile | Existing boundaries | File on disk | EPSG:4326 recommended |
| GeoJSON | Web integration | File or dict | EPSG:4326 recommended |
| EE Geometry | Programmatic | Coordinates | Any |
| DEIMS ID | Ecological sites | DEIMS ID | Automatic |
| S2 Tile | Sentinel-2 analysis | Tile code | Automatic |
| Landsat Path/Row | Landsat analysis | Path/row | Automatic |

## Working with Large Areas

For large regions, consider using multiple tiles:

```python
# Process multiple Sentinel-2 tiles
tiles = ['T30TVK', 'T30TUK', 'T30TVL']

results = []
for tile in tiles:
    ndvi = NdviSeasonality(
        roi=tile,
        sat='S2',
        periods=12,
        start_year=2023,
        end_year=2024
    )
    results.append(ndvi.get_year_composite())
```

## Converting Between Formats

### GeoDataFrame to EE Geometry

```python
import geopandas as gpd
import geemap

# Load shapefile
gdf = gpd.read_file('area.shp')

# Convert to EE Geometry
roi = geemap.geopandas_to_ee(gdf)
```

### EE Geometry to GeoJSON

```python
# Export EE geometry as GeoJSON
geometry = ee.Geometry.Rectangle([-3.8, 40.3, -3.6, 40.5])
geojson = geometry.getInfo()

# Save to file
import json
with open('roi.geojson', 'w') as f:
    json.dump(geojson, f)
```

## Best Practices

1. **Use EPSG:4326**: Shapefiles and GeoJSON work best in WGS84
2. **Simplify geometries**: Complex shapes can slow processing
3. **Test with buffers**: Start with `ee.Geometry.Point().buffer()` for quick tests
4. **Match sensor and tile**: Use S2 tiles with Sentinel-2, Landsat path/row with Landsat
5. **Check boundaries**: Visualize ROI on a map before processing

## Example: Multiple ROI Approaches

```python
import ee
from ndvi2gif import NdviSeasonality

ee.Initialize()

# Same analysis, different ROI definitions
configs = [
    {'name': 'rectangle', 'roi': ee.Geometry.Rectangle([-3.8, 40.3, -3.6, 40.5])},
    {'name': 'buffer', 'roi': ee.Geometry.Point([-3.7, 40.4]).buffer(5000)},
    {'name': 'tile', 'roi': 'T30TVK'},
    {'name': 'shapefile', 'roi': 'madrid_area.shp'}
]

for config in configs:
    print(f"Processing with {config['name']} ROI...")
    ndvi = NdviSeasonality(
        roi=config['roi'],
        sat='S2',
        periods=12,
        start_year=2023,
        end_year=2024,
        index='ndvi'
    )
    composite = ndvi.get_year_composite(year=2023)
    print(f"✓ {config['name']} complete")
```

## Next Steps

- [Quick Start](quick_start.md) - Create your first analysis
- [Basic NDVI Tutorial](../tutorials/basic_ndvi.md) - Detailed NDVI workflow
- [Multi-Sensor Comparison](../tutorials/multi_sensor.md) - Compare sensors
