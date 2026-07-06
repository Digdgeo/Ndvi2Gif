# Land Cover Classification with `LandCoverClassifier`

`LandCoverClassifier` provides supervised and unsupervised land cover classification workflows built on top of multi-temporal `NdviSeasonality` composites. All computation runs in Google Earth Engine — no data downloads are required.

---

## Concept

The classification approach uses **multi-temporal feature stacks**: instead of classifying a single-date image, the classifier ingests a stack of temporal composites (e.g. one per season per year). This greatly improves separability between classes that look similar in any single observation but differ in their seasonal trajectory — e.g., crops vs. natural grasslands, or deciduous vs. evergreen forests.

---

## Quick start

```python
import ee
import geemap
from ndvi2gif import NdviSeasonality, LandCoverClassifier

ee.Initialize(project='your-project-id')

roi = ee.Geometry.Rectangle([-6.55, 36.85, -6.10, 37.20])

ns = NdviSeasonality(
    roi=roi,
    sat='S2',
    index='ndvi',
    periods=4,             # seasonal composites
    start_year=2021,
    end_year=2022,
)

lcc = LandCoverClassifier(ns)
# LandCoverClassifier initialized for S2
# Period: 2021-2021, 4 periods/year
```

---

## Step-by-step workflow

### 1. Build the feature stack

`create_feature_stack()` generates a multi-band `ee.Image` stacking all index values across periods and years.

```python
feature_stack = lcc.create_feature_stack(
    indices=['ndvi', 'mndwi', 'nbr'],  # indices to stack
    include_statistics=True,           # add per-index mean, std, min, max
    normalize=True,                    # normalise bands to [0, 1]
)
# Feature stack ready: N bands
# (3 indices × 4 periods × 1 year) + (3 indices × 4 statistics) = 24 bands
```

The feature stack is stored in `lcc.feature_stack` and will be used for all subsequent steps.

#### Parameter reference

| Parameter | Default | Description |
|---|---|---|
| `indices` | `[ns.index]` | List of index names to include. Must be available for the configured satellite |
| `include_statistics` | `True` | Add per-index temporal statistics: mean, std, max, min |
| `normalize` | `True` | Min-max normalise all bands to [0, 1] using `reduceRegion` over the ROI |

---

### 2. Add training data

Training data can be provided as pre-labelled point features or as labelled polygons from which points are sampled.

```python
# From a shapefile with a 'class' attribute column
lcc.add_training_data(
    training_points='training_points.shp',
    class_property='class',
)

# From polygons (samples N points per class)
lcc.add_training_data(
    training_polygons='land_cover_polygons.geojson',
    class_property='class',
    points_per_class=200,
)

# From an ee.FeatureCollection
training_fc = ee.FeatureCollection('users/myproject/training_points')
lcc.add_training_data(training_points=training_fc)
```

The training data is **automatically split 70/30** into training and validation subsets using a random column.

---

### 3. Supervised classification

```python
classified = lcc.classify_supervised(algorithm='random_forest')
```

Available algorithms:

| `algorithm=` | GEE classifier | Notes |
|---|---|---|
| `'random_forest'` (default) | `ee.Classifier.smileRandomForest` | Best overall performance; provides feature importance |
| `'svm'` | `ee.Classifier.libsvm` | Good for small training sets; slower |
| `'cart'` | `ee.Classifier.smileCart` | Fast, interpretable decision tree |
| `'naive_bayes'` | `ee.Classifier.smileNaiveBayes` | Very fast baseline |
| `'gradient_tree'` | `ee.Classifier.smileGradientTreeBoost` | High accuracy, slower than RF |

#### Customising algorithm parameters

```python
# Random Forest with custom parameters
classified = lcc.classify_supervised(
    algorithm='random_forest',
    params={
        'numberOfTrees': 200,
        'variablesPerSplit': None,   # default: sqrt(n_features)
        'minLeafPopulation': 1,
        'bagFraction': 0.5,
    }
)

# SVM with RBF kernel
classified = lcc.classify_supervised(
    algorithm='svm',
    params={
        'kernelType': 'RBF',
        'gamma': 0.5,
        'cost': 10,
    }
)
```

After classification, accuracy metrics are computed automatically from the validation split.

---

### 4. Unsupervised classification (clustering)

No training data required — the algorithm groups pixels by spectral-temporal similarity.

```python
# K-means clustering
clustered = lcc.classify_unsupervised(
    algorithm='kmeans',
    n_clusters=8,
    max_iterations=20,
)
```

Available algorithms:

| `algorithm=` | Description |
|---|---|
| `'kmeans'` (default) | Weka K-means — fast and widely used |
| `'cascade_kmeans'` | Cascade K-means — automatically finds number of clusters |
| `'lda'` | Learning Vector Quantization |

> **Note:** Unsupervised results produce cluster IDs (0, 1, 2, …) that must be manually labelled by visual inspection against reference imagery.

---

### 5. Accuracy assessment

Accuracy metrics are computed automatically after supervised classification. Retrieve them programmatically:

```python
# Summary printout (shown automatically after classification)
# Overall Accuracy: 0.876
# Kappa: 0.851

# Detailed report as a DataFrame
report = lcc.get_accuracy_report()
print(report)
#           Class  ProducerAccuracy  UserAccuracy
# 0        Forest             0.921         0.903
# 1     Grassland             0.854         0.882
# 2        Wetland             0.911         0.895
# ...
```

#### Confusion matrix plot

```python
class_labels = ['Forest', 'Grassland', 'Wetland', 'Urban', 'Water', 'Bare soil']
ax = lcc.plot_confusion_matrix(labels=class_labels)
```

---

### 6. Feature importance (Random Forest)

```python
importance = lcc.get_feature_importance()
# Returns a dict {band_name: importance_score}

# Sort and display top 10 features
import pandas as pd
imp_df = pd.Series(importance).sort_values(ascending=False)
print(imp_df.head(10))
```

Feature importance reveals which temporal periods and indices are most discriminative for your classification. This is useful for:
- Reducing the feature stack to the most informative bands
- Understanding seasonal separability of land cover classes
- Selecting the optimal time windows for field campaigns

---

### 7. Visualise on a map

```python
Map = geemap.Map()
Map.centerObject(roi, zoom=10)

n_classes = 6
palette = ['1a9641', 'a6d96a', '74add1', 'd7191c', '0571b0', 'ffffbf']

Map.addLayer(
    lcc.classified_image,
    {'min': 0, 'max': n_classes - 1, 'palette': palette},
    'Land cover classification',
)
Map
```

---

### 8. Export

```python
task = lcc.export_results(
    description='donana_landcover_2021',
    scale=10,
    region=roi,
)
# Monitor at: https://code.earthengine.google.com/tasks
```

---

## Multi-index feature stacks

Including multiple indices greatly improves classification by capturing different land cover properties:

```python
# Combine vegetation, water, and burn indices
feature_stack = lcc.create_feature_stack(
    indices=['ndvi', 'mndwi', 'nbr', 'ndbi'],
    include_statistics=True,
    normalize=True,
)
```

Useful index combinations by target class:

| Target | Recommended indices |
|---|---|
| Wetlands | `mndwi`, `ndvi`, `awei` |
| Forests | `ndvi`, `evi`, `nbr` |
| Crops | `ndvi`, `evi`, `savi` |
| Urban | `ndbi`, `ndvi` |
| Burned areas | `nbr`, `nbri`, `ndvi` |

For S2, Red Edge indices add discriminative power for vegetation types:

```python
feature_stack = lcc.create_feature_stack(
    indices=['ndvi', 'ireci', 'ndre', 'mndwi'],
    include_statistics=True,
)
```

---

## Combining optical and SAR

SAR backscatter adds structural information (canopy volume, surface roughness) that is complementary to optical indices. While `LandCoverClassifier` takes a single `NdviSeasonality` instance, you can manually add SAR bands to the feature stack before training:

```python
# Build optical feature stack
ns_s2 = NdviSeasonality(roi=roi, sat='S2', periods=4, start_year=2021, end_year=2022)
lcc = LandCoverClassifier(ns_s2)
optical_stack = lcc.create_feature_stack(indices=['ndvi', 'mndwi'])

# Build SAR feature stack separately
ns_s1 = NdviSeasonality(roi=roi, sat='S1', index='rvi', periods=4,
                         start_year=2021, end_year=2022)
sar_collection = ns_s1.get_year_composite()
sar_image = ee.Image(sar_collection.first())

# Merge into a combined stack
combined_stack = optical_stack.addBands(sar_image)
lcc.feature_stack = combined_stack  # override with combined stack
```

---

## Tips

### Training data quality

Classification accuracy depends more on training data quality than on algorithm choice. Aim for:
- At least 50–100 samples per class
- Samples that cover the full intra-class variability (different dates, management regimes, vegetation density)
- Balanced class representation — highly imbalanced training sets bias the classifier towards dominant classes

### Number of trees (Random Forest)

The default of 100 trees is a good starting point. Increasing to 200–500 trees rarely hurts accuracy but adds computation time. Use `variablesPerSplit=None` (auto = √n_features) unless you have a specific reason to change it.

### Normalisation

`normalize=True` is recommended when combining indices with very different value ranges (e.g. NDVI ≈ 0–1, AWEI ≈ −2 to +2). Without normalisation, features with larger numeric ranges dominate distance-based metrics in SVM and clustering.

### Unsupervised as a pre-labelling tool

Run K-means clustering first with `n_clusters=15–20`, then visually inspect each cluster against high-resolution imagery to group clusters into meaningful land cover classes. Use those grouped clusters as training polygons for a subsequent supervised run.

---

## API reference

| Method | Returns | Description |
|---|---|---|
| `create_feature_stack(indices, include_statistics, normalize)` | `ee.Image` | Multi-temporal feature stack |
| `add_training_data(training_points, training_polygons, class_property, points_per_class)` | — | Load and split training data |
| `classify_supervised(algorithm, params)` | `ee.Image` | Supervised classification |
| `classify_unsupervised(algorithm, n_clusters, max_iterations)` | `ee.Image` | Unsupervised clustering |
| `get_accuracy_report()` | `pd.DataFrame` | Accuracy metrics table |
| `plot_confusion_matrix(labels)` | `Axes` | Confusion matrix heatmap |
| `get_feature_importance()` | `dict` | RF feature importance scores |
| `export_results(description, scale, region)` | `ee.batch.Task` | Export classified image to Drive |

---

## References

Breiman, L. (2001). Random Forests. *Machine Learning*, 45(1), 5–32.

Vapnik, V. (1998). *Statistical Learning Theory*. Wiley.

Friedman, J.H. (2001). Greedy function approximation: a gradient boosting machine. *Annals of Statistics*, 29(5), 1189–1232.

Foody, G.M. (2002). Status of land cover classification accuracy assessment. *Remote Sensing of Environment*, 80(1), 185–201.
