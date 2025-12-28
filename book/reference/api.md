# API Reference

Complete API documentation for ndvi2gif v1.0.0.

## Core Classes

### NdviSeasonality

Main class for temporal compositing and analysis.

```python
from ndvi2gif import NdviSeasonality

processor = NdviSeasonality(
    roi=None,
    sat='S2',
    periods=12,
    start_year=2020,
    end_year=2023,
    index='ndvi',
    key='median',
    percentile=75,
    mask_clouds=True,
    verbose=False
)
```

**Parameters:**

- **roi** : str, ee.Geometry, ee.FeatureCollection, or None
  - Region of interest. Can be:
    - Path to shapefile (.shp)
    - Path to GeoJSON (.geojson, .json)
    - Earth Engine Geometry object
    - Earth Engine FeatureCollection
    - DEIMS site ID (e.g., 'https://deims.org/...')
    - Sentinel-2 tile code (e.g., '29SQB')
    - Landsat path/row (e.g., 'L204033')
    - None (uses hand-drawn geometry)

- **sat** : str, default='S2'
  - Satellite/dataset platform. Options:
    - 'S2': Sentinel-2 (2015-present, 10m)
    - 'Landsat': Landsat 4-9 collection (1982-present, 30m)
    - 'MODIS': MODIS Terra+Aqua (2000-present, 500m)
    - 'S1': Sentinel-1 SAR (2014-present, 10m)
    - 'S3': Sentinel-3 OLCI (2016-present, 300m)
    - 'ERA5': ERA5-Land climate reanalysis (1950-present, ~11km)
    - 'CHIRPS': CHIRPS precipitation (1981-present, ~5.5km)

- **periods** : int, default=12
  - Number of temporal periods per year. Options:
    - 4: Seasonal (quarterly)
    - 12: Monthly
    - 24: Bi-monthly
    - Custom: Any integer

- **start_year** : int, default=2020
  - First year of analysis

- **end_year** : int, default=2023
  - Last year of analysis (inclusive)

- **index** : str, default='ndvi'
  - Variable/index to compute. See [Indices Reference](indices.md) for complete list.
  - Optical: 'ndvi', 'evi', 'ndwi', 'mndwi', 'savi', 'gndvi', 'ndmi', etc. (40+ options)
  - SAR: 'rvi', 'dpsvi', 'rfdi', 'vsdi', etc.
  - Climate (ERA5): 'temperature_2m', 'total_precipitation_sum', 'volumetric_soil_water_layer_1', etc. (47 variables)
  - CHIRPS: 'precipitation'

- **key** : str, default='median'
  - Statistical reducer for temporal compositing. Options:
    - 'median': Robust to outliers (recommended for optical)
    - 'max': Maximum value (peak greenness)
    - 'mean': Arithmetic mean
    - 'percentile': Use specified percentile
    - 'min': Minimum value
    - 'sum': Sum (for precipitation/flux variables)

- **percentile** : int, default=75
  - Percentile value when key='percentile' (0-100)

- **mask_clouds** : bool, default=True
  - Apply cloud masking for optical sensors (S2, Landsat, MODIS)

- **verbose** : bool, default=False
  - Print processing information

**Attributes:**

- **roi** : ee.Geometry - Processed region of interest
- **ndvi_col** : ee.ImageCollection - Satellite/climate collection
- **imagelist** : list - List of composites (one per year)
- **period_names** : list - Period labels (e.g., ['Jan', 'Feb', ...])
- **period_dates** : list - Date ranges for each period

**Methods:**

#### get_year_composite()

Generate temporal composite for a single year.

```python
composite = processor.get_year_composite(year=2023)
```

**Parameters:**
- **year** : int, optional - Year to process (default: start_year)

**Returns:**
- **ee.Image** - Multi-band image with one band per period

---

#### get_period_composite()

Generate composite for a specific period across all years.

```python
period_composite = processor.get_period_composite(period=0)
```

**Parameters:**
- **period** : int - Period index (0 to periods-1)

**Returns:**
- **ee.Image** - Multi-band image with one band per year

---

#### get_gif()

Generate animated GIF from temporal composites.

```python
processor.get_gif(
    name='animation.gif',
    fps=2,
    figsize=(12, 10),
    palette='RdYlGn',
    vmin=0,
    vmax=0.8,
    title='NDVI Evolution'
)
```

**Parameters:**
- **name** : str - Output filename
- **fps** : int, default=2 - Frames per second
- **figsize** : tuple, default=(12, 10) - Figure size (width, height)
- **palette** : str, default='RdYlGn' - Matplotlib colormap
- **vmin** : float - Minimum value for color scale
- **vmax** : float - Maximum value for color scale
- **title** : str - Plot title

---

#### get_export()

Export composite as GeoTIFF to local file.

```python
processor.get_export(
    year=2023,
    filename='ndvi_2023.tif',
    scale=10
)
```

**Parameters:**
- **year** : int - Year to export
- **filename** : str - Output filename
- **scale** : int, optional - Export resolution in meters (default: native resolution)

---

#### get_stats()

Calculate zonal statistics for ROI.

```python
stats = processor.get_stats(year=2023)
```

**Parameters:**
- **year** : int - Year to analyze

**Returns:**
- **dict** - Statistics dictionary with mean/min/max/std per period

---

### TimeSeriesAnalyzer

Time series analysis and visualization tools.

```python
from ndvi2gif import TimeSeriesAnalyzer

analyzer = TimeSeriesAnalyzer(processor)
```

**Parameters:**
- **processor** : NdviSeasonality - Configured processor object

**Methods:**

#### extract_time_series()

Extract time series data from composites.

```python
df = analyzer.extract_time_series(
    point=(-6.48, 37.13),
    scale=10
)
```

**Parameters:**
- **point** : tuple, optional - (lon, lat) coordinates for point extraction
- **scale** : int, default=10 - Spatial resolution in meters

**Returns:**
- **pandas.DataFrame** - Time series with columns: date, value, year, period, doy, season, month

---

#### analyze_trend()

Detect statistical trends in time series.

```python
results = analyzer.analyze_trend(
    df=df,
    method='all',
    alpha=0.05
)
```

**Parameters:**
- **df** : pandas.DataFrame - Time series from extract_time_series()
- **method** : str, default='all' - Trend method ('mann_kendall', 'sen_slope', 'linear', or 'all')
- **alpha** : float, default=0.05 - Significance level

**Returns:**
- **dict** - Trend statistics and interpretation

---

#### extract_phenology()

Extract phenology metrics from vegetation time series.

```python
phenology = analyzer.extract_phenology(df=df)
```

**Parameters:**
- **df** : pandas.DataFrame - Time series data

**Returns:**
- **dict** - Phenology metrics: SOS, EOS, POS, length of season, growth/senescence rates

---

#### plot_comprehensive_analysis()

Generate multi-panel analysis dashboard.

```python
fig = analyzer.plot_comprehensive_analysis(
    df=df,
    show_trend=True,
    show_phenology=True
)
```

**Parameters:**
- **df** : pandas.DataFrame, optional - Pre-extracted time series (if None, extracts automatically)
- **show_trend** : bool, default=True - Include trend analysis
- **show_phenology** : bool, default=True - Include phenology metrics (vegetation only)

**Returns:**
- **matplotlib.Figure** - Multi-panel figure

---

### S1ARDProcessor

Sentinel-1 SAR Analysis Ready Data preprocessing.

```python
from ndvi2gif import S1ARDProcessor

s1_processor = S1ARDProcessor(
    speckle_filter='REFINED_LEE',
    terrain_correction=True,
    terrain_flattening_model='VOLUME',
    dem='COPERNICUS_30'
)
```

**Parameters:**

- **speckle_filter** : str, default='REFINED_LEE'
  - Speckle filter algorithm. Options:
    - 'REFINED_LEE': Refined Lee filter (recommended)
    - 'LEE': Lee filter
    - 'LEE_SIGMA': Lee Sigma filter
    - 'GAMMA_MAP': Gamma MAP filter
    - 'BOXCAR': Simple boxcar filter

- **terrain_correction** : bool, default=True
  - Apply radiometric terrain correction

- **terrain_flattening_model** : str, default='VOLUME'
  - Scattering model for terrain flattening. Options:
    - 'VOLUME': Volume scattering (recommended for vegetation)
    - 'SURFACE': Surface scattering

- **dem** : str, default='COPERNICUS_30'
  - Digital elevation model. Options:
    - 'COPERNICUS_30': Copernicus DEM 30m (recommended)
    - 'SRTM': SRTM 30m

**Methods:**

#### preprocess()

Apply full preprocessing pipeline to Sentinel-1 collection.

```python
processed = s1_processor.preprocess(
    collection=s1_collection,
    roi=roi,
    orbit='DESCENDING'
)
```

**Parameters:**
- **collection** : ee.ImageCollection - Raw Sentinel-1 collection
- **roi** : ee.Geometry - Region of interest
- **orbit** : str, optional - Filter by orbit direction ('ASCENDING', 'DESCENDING')

**Returns:**
- **ee.ImageCollection** - Preprocessed SAR collection

---

### LandCoverClassifier

Machine learning classification for land cover mapping.

```python
from ndvi2gif.clasification import LandCoverClassifier

classifier = LandCoverClassifier(
    roi=roi,
    training_data=training_samples,
    classifier_type='randomForest',
    n_trees=100
)
```

**Parameters:**

- **roi** : ee.Geometry - Region of interest
- **training_data** : ee.FeatureCollection - Training samples with 'class' property
- **classifier_type** : str, default='randomForest'
  - Classification algorithm. Options:
    - Supervised: 'randomForest', 'svm', 'cart', 'naiveBayes', 'gradientTreeBoost'
    - Unsupervised: 'kmeans', 'cascadeKmeans', 'lda'
- **n_trees** : int, default=100 - Number of trees (for tree-based algorithms)
- **max_nodes** : int, optional - Maximum nodes per tree

**Methods:**

#### classify()

Train classifier and classify image.

```python
classified = classifier.classify(image=composite)
```

**Parameters:**
- **image** : ee.Image - Multi-band image to classify

**Returns:**
- **ee.Image** - Classified image with 'classification' band

---

#### accuracy_assessment()

Calculate accuracy metrics with confusion matrix.

```python
accuracy = classifier.accuracy_assessment(validation_data=validation_samples)
```

**Parameters:**
- **validation_data** : ee.FeatureCollection - Independent validation samples

**Returns:**
- **dict** - Accuracy statistics (overall accuracy, kappa, confusion matrix)

---

#### feature_importance()

Get feature importance for tree-based classifiers.

```python
importance = classifier.feature_importance()
```

**Returns:**
- **dict** - Feature importance scores

---

## Utility Functions

### scale_OLI()

Apply scale factors to Landsat 8-9 (OLI) imagery.

```python
from ndvi2gif import scale_OLI

scaled_image = scale_OLI(image)
```

---

### scale_ETM()

Apply scale factors to Landsat 4-7 (TM/ETM+) imagery.

```python
from ndvi2gif import scale_ETM

scaled_image = scale_ETM(image)
```

---

## Examples

### Basic Vegetation Monitoring

```python
import ee
from ndvi2gif import NdviSeasonality, TimeSeriesAnalyzer

ee.Initialize()

# Create processor
ndvi = NdviSeasonality(
    roi='path/to/study_area.shp',
    sat='S2',
    periods=12,
    start_year=2020,
    end_year=2023,
    index='ndvi'
)

# Generate GIF
ndvi.get_gif('ndvi_animation.gif', fps=2)

# Analyze trends
analyzer = TimeSeriesAnalyzer(ndvi)
df = analyzer.extract_time_series()
trends = analyzer.analyze_trend(df)
print(trends['interpretation'])
```

### Climate Analysis

```python
# Temperature analysis
temp = NdviSeasonality(
    roi=ee.Geometry.Point([-6.48, 37.13]).buffer(10000),
    sat='ERA5',
    index='temperature_2m_celsius',
    periods=12,
    start_year=2020,
    end_year=2023,
    key='mean'
)

analyzer = TimeSeriesAnalyzer(temp)
df_temp = analyzer.extract_time_series(scale=11000)
fig = analyzer.plot_comprehensive_analysis(df_temp)
```

### SAR Processing

```python
from ndvi2gif import NdviSeasonality, S1ARDProcessor

# Configure SAR processor
s1_proc = S1ARDProcessor(
    speckle_filter='REFINED_LEE',
    terrain_correction=True
)

# Use with NdviSeasonality
sar = NdviSeasonality(
    roi=roi,
    sat='S1',
    periods=12,
    start_year=2023,
    end_year=2023,
    index='rvi'  # SAR vegetation index
)
```

### Land Cover Classification

```python
from ndvi2gif import NdviSeasonality
from ndvi2gif.clasification import LandCoverClassifier

# Create multi-temporal feature stack
processor = NdviSeasonality(
    roi=roi,
    sat='S2',
    periods=12,
    start_year=2023,
    end_year=2023,
    index='ndvi'
)
features = processor.get_year_composite()

# Train and classify
classifier = LandCoverClassifier(
    roi=roi,
    training_data=training_samples,
    classifier_type='randomForest',
    n_trees=100
)
classified = classifier.classify(features)

# Assess accuracy
accuracy = classifier.accuracy_assessment(validation_samples)
print(f"Overall accuracy: {accuracy['overall_accuracy']:.2%}")
```

---

## See Also

- [Indices Reference](indices.md) - Complete list of 88 variables
- [Datasets Reference](datasets.md) - Details on 7 satellite/climate platforms
- [Tutorials](../tutorials/basic_ndvi.md) - Step-by-step guides
- [GitHub Repository](https://github.com/Digdgeo/Ndvi2Gif) - Source code and examples
