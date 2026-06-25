# Quick start

This page shows a minimal end-to-end example: initialize Earth Engine, create seasonal composites, build a feature stack, classify, and export.

## Installation

```bash
pip install ndvi2gif
# Earth Engine (first time only)
pip install geemap earthengine-api
```

Authenticate once if needed:
```python
import ee
ee.Authenticate()  # only once on a new machine/session
ee.Initialize()
```

## Minimal example (Landsat, unsupervised)

```python
import ee
from ndvi2gif import NdviSeasonality, LandCoverClassifier

ee.Initialize()

# 1) Area of interest (change to your AOI)
roi = ee.Geometry.Rectangle([-4.80, 37.80, -4.60, 37.95])

# 2) Processor configuration (Landsat 8/9, seasonal = 4 periods/year)
proc = NdviSeasonality(
    roi=roi,
    start_year=2020,
    end_year=2023,      # exclusive -> 2020–2022
    periods=4,
    sat="Landsat",
    index="ndvi",
    stats=["max"]
)

# 3) Feature stack (NDVI + MNDWI + LST), normalized to [0,1]
clf = LandCoverClassifier(proc)
feature_img = clf.create_feature_stack(
    indices=["ndvi", "mndwi", "lst"],
    include_statistics=True,
    normalize=True
)

# 4) Unsupervised classification (K-means)
clusters = clf.classify_unsupervised(algorithm="kmeans", n_clusters=8)

# 5) Export (Drive or local/Asset). Example: to Google Drive
clusters = clusters.rename("class").toInt16()
task = proc.export_to_drive(
    image=clusters,
    description="landcover_kmeans_2020_2022",
    folder="ndvi2gif",
    scale=30,           # Landsat native resolution
    crs="EPSG:4326"
)
print("Drive task:", task.id)
```

## Supervised classification (skeleton)

```python
# Prepare training data (choose ONE of the following)

# A) Polygons with a 'class' field (points sampled per class)
# clf.add_training_data(
#     training_polygons="/path/to/polygons.shp",  # or .geojson or ee.FeatureCollection
#     class_property="class",
#     points_per_class=200,
#     scale=30   # Landsat; use 10 for Sentinel-2
# )

# B) Labeled points
# clf.add_training_data(
#     training_points="/path/to/points.geojson",
#     class_property="class",
#     scale=30
# )

# Train & classify (Random Forest by default)
# classified = clf.classify_supervised(
#     algorithm="random_forest",
#     params={"numberOfTrees": 200, "bagFraction": 0.6}
# )

# Optional: accuracy report & confusion matrix
# print(clf.accuracy_results)
# labels = ["Water", "Urban", "Crop", "Forest", "Bare"]  # adjust to your classes
# _ = clf.plot_confusion_matrix(labels)

# Export (Drive/Asset/local)
# classified = classified.rename("class").toInt16()
# task = proc.export_to_asset(
#     image=classified,
#     asset_id="users/your_user/ndvi2gif/landcover_rf_2020_2022",
#     pyramiding_policy={"class": "mode"},
#     scale=30,
#     overwrite=True
# )
```

## Notes

- **Scale**: use the sensor’s native resolution (S2/S1 → 10 m; Landsat → 30 m; MODIS → 250–500 m; S3 → ~300 m).
- **Types**: export classifications as integers (`toInt16()`); set pyramiding policy `{"class": "mode"}` for Assets.
- The full tutorial is available at: `tutorials/landcover_classification.ipynb`.
