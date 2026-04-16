# Land Cover Classification with `LandCoverClassifier`

`LandCoverClassifier` is the machine-learning layer on top of `NdviSeasonality`. It takes the multi-temporal composite you already configured — any sensor, any index, any period scheme — turns it into a feature stack, and feeds it into one of the supervised or unsupervised classifiers available in Google Earth Engine. Training, prediction, and accuracy assessment all run server-side.

The design intent: you should not have to rewrite your feature engineering every time you want to try a new classifier. The `NdviSeasonality` configuration *is* the feature engineering.

---

## Quick start

```python
import ee
from ndvi2gif import NdviSeasonality, LandCoverClassifier

ee.Initialize()

ns = NdviSeasonality(
    roi=roi,
    sat='S2',
    index='ndvi',
    start_year=2021, end_year=2023,
    periods=4,
    key='median',
)

clf = LandCoverClassifier(ns)

clf.create_feature_stack(indices=['ndvi', 'ndwi', 'bsi'], include_statistics=True)
clf.add_training_data(training_polygons='training.shp', class_property='class')

classified = clf.classify_supervised(algorithm='random_forest')

print(clf.accuracy_results['overall_accuracy'], clf.accuracy_results['kappa'])
```

---

## 1. Building a feature stack

`create_feature_stack()` stacks every `period × year × index` combination into a single multi-band `ee.Image`. A three-year Sentinel-2 seasonal composite with three indices produces 3 × 4 × 3 = 36 bands.

```python
stack = clf.create_feature_stack(
    indices=['ndvi', 'ndwi', 'bsi'],
    include_statistics=True,   # adds mean, std, max, min per index
    normalize=True,            # rescales all bands to [0, 1]
)
```

| Argument | Purpose |
|---|---|
| `indices` | Indices to stack. If `None`, uses the current `NdviSeasonality.index`. Must be compatible with the sensor. |
| `include_statistics` | Adds four temporal summary bands per index (`mean`, `std`, `max`, `min`) — often the most informative features for vegetation classes. |
| `normalize` | Linear min-max rescaling to `[0, 1]`. Recommended for SVM and Naive Bayes, optional for tree-based methods. |

Band names follow the pattern `{index}_{year}_{period}` — e.g. `ndvi_2021_summer`, `ndwi_2022_winter`, `bsi_mean`.

> **Why multi-temporal helps.** A cropland pixel and a shrubland pixel can have identical summer NDVI. What separates them is the seasonal *trajectory* — the winter trough, the spring rise, the post-harvest drop. Stacking the periods makes that trajectory available to the classifier as features, instead of collapsing it into a single statistic.

### Sensor compatibility

The stack inherits the sensor from `NdviSeasonality`. Invalid indices (not registered for that sensor in `self.sensor_indices`) raise a `ValueError` before any GEE work runs.

---

## 2. Training data

Two input formats are supported: **points** (already labelled) or **polygons** (a class label per polygon, from which points are sampled).

```python
# Option A — labelled points
clf.add_training_data(
    training_points='training_points.shp',
    class_property='class',
)

# Option B — labelled polygons, sample 100 points each
clf.add_training_data(
    training_polygons='training_polygons.geojson',
    class_property='class',
    points_per_class=100,
)

# Option C — pass an ee.FeatureCollection you already built
clf.add_training_data(training_points=my_ee_fc, class_property='class')
```

A few things happen under the hood:

1. If you pass polygons, `sampleRegions()` draws random samples inside each polygon.
2. The feature stack is sampled at every training location — you get one feature vector per sample.
3. A random 70/30 split creates `self.training_data` and `self.validation_data`.
4. Sample counts are printed — inspect these. If a class has fewer than ~30 validation samples, the per-class accuracy metrics downstream will be noisy.

> **Polygon quality matters more than point count.** One hundred carefully-drawn polygons over pure, representative patches of each class will outperform a thousand cheaply-sampled points every time. Avoid mixed polygons — the model cannot un-mix them for you.

---

## 3. Supervised classification

```python
classified = clf.classify_supervised(algorithm='random_forest')
```

Five algorithms are available. All are Earth Engine classifiers (`ee.Classifier.*`) — nothing leaves GEE.

| Algorithm | Key parameters | Strengths | Weaknesses |
|---|---|---|---|
| `random_forest` | `numberOfTrees=100`, `variablesPerSplit`, `bagFraction=0.5` | **Default.** Handles high-dimensional stacks, robust to noise, gives feature importance. | Can overfit with tiny training sets. |
| `svm` | `kernelType='RBF'`, `gamma=0.5`, `cost=10` | Effective on small, clean training sets. | Sensitive to scaling — always `normalize=True`. Slow on big ROIs. |
| `cart` | `maxNodes`, `minLeafPopulation=1` | Fast, interpretable decision tree. | High variance; prefer Random Forest unless you need a single tree. |
| `naive_bayes` | (none) | Extremely fast baseline. | Assumes feature independence — usually violated for temporal stacks. |
| `gradient_tree` | `numberOfTrees=50`, `shrinkage=0.05`, `samplingRate=0.7` | Often the highest accuracy when well-tuned. | More hyperparameters to tune. |

Pass algorithm-specific parameters via the `params` dict:

```python
classified = clf.classify_supervised(
    algorithm='random_forest',
    params={'numberOfTrees': 300, 'minLeafPopulation': 5},
)
```

After training, the classifier is applied to the full feature stack and stored in `self.classified_image`. Accuracy against the 30 % validation split is computed automatically.

### Accuracy assessment

```python
clf.accuracy_results['overall_accuracy']   # 0.89
clf.accuracy_results['kappa']              # 0.86
clf.accuracy_results['producers_accuracy'] # per-class, omission errors
clf.accuracy_results['consumers_accuracy'] # per-class, commission errors
clf.accuracy_results['confusion_matrix']   # raw ndarray

clf.plot_confusion_matrix(labels=['water', 'crop', 'forest', 'urban', 'bare'])
```

Confusion matrix interpretation:

- Rows = reference class, columns = predicted class.
- Diagonal = correctly classified samples.
- Off-diagonal cells tell you *which* classes the model confuses. "Crop → bare" confusion usually means your harvest timing moved out of the composite window.

> **Kappa below 0.6 means reconsider, not retune.** If overall accuracy is high but Kappa lags, the classifier is exploiting class imbalance. Before adding more trees or tweaking `gamma`, check whether the sample distribution matches the real landscape. More often than not the problem is the training data, not the model.

### Feature importance (Random Forest only)

```python
importance = clf.get_feature_importance()
# dict mapping band name → importance score
```

Very useful for pruning: if `ndvi_2022_winter` dominates and `bsi_2021_autumn` contributes nothing, drop it from the next run. Smaller feature stacks train faster and generalise better.

---

## 4. Unsupervised classification

When you don't have training data — or want to explore the structure of the feature space before building it:

```python
clustered = clf.classify_unsupervised(
    algorithm='kmeans',
    n_clusters=8,
    max_iterations=20,
)
```

| Algorithm | Description |
|---|---|
| `kmeans` | Weka k-means. Requires you to choose `n_clusters`. |
| `cascade_kmeans` | Cascade k-means; chooses an optimal number of clusters between 2 and `n_clusters`. |
| `lda` | Weka LVQ (Learning Vector Quantization). |

The clusterer is trained on a random sample of 5 000 pixels drawn from the feature stack, then applied to the whole image. Clusters are unlabelled — interpret them by overlaying with known land cover or high-resolution imagery.

> **Cascade k-means is the right first move on unknown terrain.** You rarely know the "true" number of classes before looking at the data. Cascade will suggest one, then you can iterate with plain k-means once you have a hypothesis.

---

## 5. Exporting results

### Quick local export

```python
from ndvi2gif import NdviSeasonality
ns.get_export_single(classified, 'landcover_2022.tif', scale=10)
```

This downloads through `geemap.ee_export_image` — fine for small ROIs. For larger areas, use Drive or Asset exports:

```python
ns.export_to_drive(
    image=classified,
    description='landcover_2022_RF',
    folder='ndvi2gif_exports',
    scale=10,
)

ns.export_to_asset(
    image=classified,
    asset_id='users/you/ndvi2gif/landcover_2022',
    pyramiding_policy={'classification': 'mode'},
    scale=10,
)
```

**`pyramiding_policy='mode'` is important for class maps.** The default (`mean`) averages class IDs across pyramid levels and produces meaningless values at zoom-out.

`LandCoverClassifier` also exposes a thin wrapper:

```python
task = clf.export_results(description='landcover_2022', scale=10)
```

See the [Export options tutorial](export_options.md) for the full set of export strategies.

---

## 6. A realistic workflow

```python
ns = NdviSeasonality(
    roi=my_roi,
    sat='S2',
    index='ndvi',
    start_year=2021, end_year=2023,
    periods=4,
    key='median',
)

clf = LandCoverClassifier(ns)

# 1. Rich feature stack: three indices, temporal statistics, normalized
clf.create_feature_stack(
    indices=['ndvi', 'ndwi', 'bsi'],
    include_statistics=True,
    normalize=True,
)

# 2. Polygon-based training
clf.add_training_data(
    training_polygons='training.shp',
    class_property='class',
    points_per_class=150,
)

# 3. First pass with Random Forest (default)
classified = clf.classify_supervised(algorithm='random_forest')
print('RF accuracy:', clf.accuracy_results['overall_accuracy'])

# 4. Inspect importance, prune low-value features, retry
importance = clf.get_feature_importance()

# 5. Export final map as a GEE asset (mode pyramiding for class maps)
ns.export_to_asset(
    image=classified,
    asset_id='users/you/ndvi2gif/landcover_2022',
    pyramiding_policy={'classification': 'mode'},
    scale=10,
)
```

---

## Tips and caveats

### Class imbalance

If one class represents 80 % of your training samples, the overall accuracy can look great while minority classes are effectively ignored. Prefer Kappa and per-class producer's/consumer's accuracy over overall accuracy, and consider stratified polygon sampling rather than random point sampling.

### Temporal leakage

Including bands from the *same year* as your reference labels inflates accuracy in a way that doesn't generalise. If your labels are for 2022, train on 2021 features and test on 2022 — or vice versa. The current pipeline does not enforce this split; it's your responsibility.

### Scale vs memory

`create_feature_stack()` normalizes using `reduceRegion` at `scale=30`, and training sampling runs at `scale=10`. For very large ROIs the min-max reduction can exceed the default `maxPixels`. If you hit that wall, coarsen the ROI, tile it, or disable normalization (`normalize=False`) and use a tree-based classifier that is scale-insensitive.

### Random Forest vs Gradient Boosting

Start with Random Forest. Try Gradient Boosting once you have a clean baseline and want to squeeze out the last few accuracy points. Don't spend effort on hyperparameter tuning before your training data is solid.

### Interpretability

Random Forest plus `get_feature_importance()` is the only route to "why did the model decide this?". If that matters for your use case (e.g. stakeholder reports), commit to RF and skip SVM.

---

## API reference

### `LandCoverClassifier(ndvi_seasonality_instance)`

| Method | Returns | Description |
|---|---|---|
| `create_feature_stack(indices, include_statistics, normalize)` | `ee.Image` | Build multi-band feature stack from temporal composites |
| `add_training_data(training_points, training_polygons, class_property, points_per_class)` | `None` | Load labelled samples, auto-split 70/30 train/validation |
| `classify_supervised(algorithm, train_fraction, params)` | `ee.Image` | Train and apply a supervised classifier |
| `classify_unsupervised(algorithm, n_clusters, max_iterations, params)` | `ee.Image` | Apply a clustering algorithm to the feature stack |
| `export_results(description, scale, region)` | `ee.batch.Task` | Export classified image to Google Drive |
| `plot_confusion_matrix(labels)` | `matplotlib.axes.Axes` | Confusion matrix heatmap |
| `get_accuracy_report()` | `pandas.DataFrame` | Tabular accuracy metrics |
| `get_feature_importance()` | `dict` | Feature importance (Random Forest only) |

**Supervised algorithms:** `random_forest`, `svm`, `cart`, `naive_bayes`, `gradient_tree`

**Unsupervised algorithms:** `kmeans`, `cascade_kmeans`, `lda`

---

## References

Breiman, L. (2001). Random Forests. *Machine Learning*, 45(1), 5–32.

Friedman, J.H. (2001). Greedy function approximation: a gradient boosting machine. *Annals of Statistics*, 29(5), 1189–1232.

Cohen, J. (1960). A coefficient of agreement for nominal scales. *Educational and Psychological Measurement*, 20(1), 37–46.

Congalton, R.G., Green, K. (2019). *Assessing the Accuracy of Remotely Sensed Data: Principles and Practices* (3rd ed.). CRC Press.

Gorelick, N., Hancher, M., Dixon, M., Ilyushchenko, S., Thau, D., Moore, R. (2017). Google Earth Engine: Planetary-scale geospatial analysis for everyone. *Remote Sensing of Environment*, 202, 18–27.
