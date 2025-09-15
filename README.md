# Ndvi2Gif: Multi-Seasonal Remote Sensing Index Composites

[![PyPI version](https://img.shields.io/pypi/v/ndvi2gif.svg)](https://pypi.org/project/ndvi2gif/)
[![PyPI downloads](https://img.shields.io/pypi/dm/ndvi2gif.svg)](https://pypi.org/project/ndvi2gif/)
[![Conda version](https://img.shields.io/conda/vn/conda-forge/ndvi2gif.svg)](https://anaconda.org/conda-forge/ndvi2gif)
[![Conda downloads](https://img.shields.io/conda/dn/conda-forge/ndvi2gif.svg)](https://anaconda.org/conda-forge/ndvi2gif)
[![Build status](https://github.com/Digdgeo/Ndvi2Gif/actions/workflows/python-publish.yml/badge.svg)](https://github.com/Digdgeo/Ndvi2Gif/actions/workflows/python-publish.yml)

![NDVI2GIF Köln](https://i.imgur.com/Y5dOWIk.jpeg)
*Richter's stained glass in Cologne Cathedral. Inspiration for this library.*

**Ndvi2Gif** is a Python library designed to simplify access to global satellite data through the Google Earth Engine platform. While its name highlights the ability to create seasonal GIF animations, the true power of this tool lies in its capability to compute and export pixel-wise statistics for any region on Earth, across any time span covered by supported remote sensing datasets.

Built on top of [Google Earth Engine](https://github.com/google/earthengine-api) and [Geemap](https://github.com/giswqs/geemap), it allows you to:

- Generate annual or multi-annual composited rasters (e.g., median NDVI per season between 2001 and 2020),
- Apply multiple statistics (mean, max, flexible percentiles) across space and time,
- Export results as GeoTIFFs for further analysis,
- Retrieve zonal statistics over user-defined geometries,
- Monitor vegetation structure with advanced SAR indices,
- Handle incomplete years automatically for real-time monitoring,
- **NEW in v0.6.0:** Perform supervised and unsupervised land cover classification with integrated machine learning,
- **NEW in v0.6.0:** Export directly to Google Drive and Earth Engine Assets,
- And yes — also create colorful GIFs for easy visualization.

Whether you're monitoring crop phenology, detecting harvest events, assessing drought trends, classifying land cover, or preparing input layers for further ecological modeling, `ndvi2gif` makes it easier to extract reliable, multi-temporal remote sensing information at scale.

Ndvi2Gif was updated and extended as part of its integration into the eLTER and SUMHAL projects, which also enabled the use of eLTER site boundaries (via `deimsPy`) as one of its input sources.

![Interface Screenshot](https://i.imgur.com/Sv9LfYj.png)

## ✨ What's New in v0.6.0 - Machine Learning & Classification Release 🧠

The **0.6.0 release** transforms *Ndvi2Gif* into a **complete remote sensing analysis platform** with integrated machine learning capabilities, positioning it as a comprehensive solution for remote sensing workflows.

### 🚀 New Classification Capabilities

- 🧠 **LandCoverClassifier** – Complete supervised and unsupervised classification workflows
- 🎯 **Multiple Algorithms** – Random Forest, SVM, CART, Naive Bayes, Gradient Tree Boost, K-means, LDA
- 📊 **Accuracy Assessment** – Confusion matrices and comprehensive validation reports
- 🔧 **Feature Engineering** – Multi-temporal stacks with automatic normalization
- 📤 **Enhanced Exports** – Direct export to Google Drive and Earth Engine Assets

### 📚 Documentation Overhaul

- **95%+ Documentation Coverage** – Comprehensive Sphinx-style docstrings
- **Scientific References** – All indices now include citations
- **Complete Examples** – Every method includes usage examples
- **Better Error Handling** – Informative messages with suggested solutions

### 🛰️ Enhanced Capabilities

- **export_to_drive()** – Batch export with full parameter control
- **export_to_asset()** – Direct Earth Engine Asset creation with pyramiding policies
- **Automatic scale detection** – Sensor-specific resolution handling
- **Improved SAR processing** – Enhanced error handling and documentation
- **Feature importance analysis** – Understand which indices contribute most to classification

### 🔥 All Previous Features, Better Than Ever

- 🛰️ **Sentinel-1 ARD Processor** – Professional SAR preprocessing with terrain correction
- 📈 **TimeSeriesAnalyzer** – Extract robust time series, test for trends, and visualize dynamics
- 🌱 **Extended NdviSeasonality** – Dynamic temporal periods (4, 12, 24, custom)
- 🎨 **Polished Visualizations** – Publication-ready layouts

## Why use Ndvi2Gif?

Unlike many visualization-oriented tools, Ndvi2Gif is designed as a **remote sensing analytics suite** that abstracts much of the complexity of working directly with Google Earth Engine, while giving you the flexibility to go far beyond GIF creation.

You can:

- **Access pixel-wise statistics** over any Earth location, at any scale and time span.  
  - Example: *Obtain the monthly median of the 85th NDVI percentile per pixel from 1984 to 2024 using Landsat data.*  
  - Example: *Calculate the maximum of the seasonal NDWI maximums between 2017 and 2023 using Sentinel-2.*  
  - Example: *Monitor crop harvest timing with bi-monthly VV/VH ratio analysis using Sentinel-1.*  
  - Example: *Track daily algal blooms with Sentinel-3 OLCI turbidity indices.*  

- **Perform advanced machine learning classification** (NEW in v0.6.0):
  - Multi-temporal land cover mapping with Random Forest
  - Crop type classification with SVM
  - Unsupervised clustering with K-means
  - Feature importance analysis for ecological insights

- **Perform nested aggregations**:  
  First compute temporal summaries (e.g., per-season percentiles or means), then apply a second statistical reduction across years (e.g., median, min, max).

- **Run advanced time series analysis** with the `TimeSeriesAnalyzer`:  
  - Trend detection (Mann-Kendall, Sen's slope, linear regression)  
  - Multi-panel dashboards (seasonal patterns, autocorrelation, data quality)  
  - Phenology metrics such as Start/End of Season, Peak, Length, amplitude, and rates of change  

- **Preprocess Sentinel-1 SAR like a pro** with the `S1ARDProcessor`:  
  - Radiometric terrain correction for mountainous regions  
  - Multiple speckle filtering options (Boxcar, Lee, Refined Lee, Gamma-MAP, Lee Sigma)  
  - Flexible DEM support (Copernicus and SRTM)  

- **Target any ecological or phenological metric** by choosing the appropriate index and analysis pipeline.

- **Work globally**, without needing to download or preprocess raw satellite data — all computations are handled via Earth Engine's cloud infrastructure.

- **Handle real-time monitoring** with automatic detection of available data periods for incomplete years.

In other words: if you can describe a temporal range, a spatial region, an index, and a chain of statistics — `ndvi2gif` can not only generate it, but now also help you **classify, analyze and interpret the changes over time**.

Yes, it makes nice GIFs — but it's much more than that.
![GIF Example](https://i.imgur.com/xvrPYMH.gif)
![RGB Example](https://i.imgur.com/tq4aMBv.jpg)
*Crop pattern dance around Los Palacios y Villafranca (SW Spain) and the palette color combinations shown*

### Supported Input Formats for ROI

| Input Type           | Description                                                 | Example / Notes                                      |
|----------------------|-------------------------------------------------------------|------------------------------------------------------|
| Drawn Geometry       | Use geemap to draw a polygon directly on a map             | Works in Jupyter Notebooks                           |
| Shapefile / GeoJSON  | Provide a file path to a vector dataset                    | EPSG:4326 recommended                                |
| eLTER site ID        | Use `deimsPy` to fetch site boundaries by DEIMS ID         | e.g., `deimsid:ab8278e6-0b71-4b36-a6d2-e8f34aa3df30` |
| Sentinel-2 Tile      | Specify MGRS tile code (e.g., `T30TYN`)                    | Automatically fetches tile geometry                  |
| Landsat Path/Row     | Provide WRS-2 path and row codes (e.g., `198/034`)         | Covers full Landsat archive                          |

## Included Statistics

- **Maximum** - Peak values for cloud-free compositing
- **Mean** - Average values across time period
- **Median** - Robust central tendency, excellent for noisy data
- **Flexible Percentiles** - Any percentile from 1 to 99
  - Custom percentiles like 75th, 85th, or 99th for specific applications
  - Perfect for handling varying cloud contamination levels

## Available Indices

### 🌱 Basic Optical Indices (S2, Landsat, MODIS, S3)
- **NDVI** - Normalized Difference Vegetation Index
- **EVI** - Enhanced Vegetation Index  
- **GNDVI** - Green Normalized Difference Vegetation Index
- **SAVI** - Soil Adjusted Vegetation Index
- **NDWI** - Normalized Difference Water Index
- **MNDWI** - Modified Normalized Difference Water Index
- **AWEI** - Automated Water Extraction Index
- **AEWINSH** - AWEI No Shadow
- **NDSI** - Normalized Difference Snow Index
- **NBRI** - Normalized Burn Ratio Index
- **NDMI** - Normalized Difference Moisture Index

### 🌾 Advanced Optical Indices (S2, Landsat, MODIS, S3)
- **MSI** - Moisture Stress Index (drought monitoring)
- **NMI** - Normalized Multi-band Drought Index
- **NDTI** - Normalized Difference Tillage Index
- **CRI1/CRI2** - Carotenoid Reflectance Indices
- **LAI** - Leaf Area Index approximation
- **PRI** - Photochemical Reflectance Index
- **WDRVI** - Wide Dynamic Range Vegetation Index

### 🔬 Sentinel-2 Exclusive (Red Edge B5-B7)
- **IRECI** - Inverted Red-Edge Chlorophyll Index (high sensitivity chlorophyll)
- **MCARI** - Modified Chlorophyll Absorption Ratio Index
- **NDRE** - Normalized Difference Red Edge (chlorophyll content)
- **REIP** - Red Edge Inflection Point (vegetation stress)
- **PSRI** - Plant Senescence Reflectance Index (crop maturity)
- **CIRE** - Chlorophyll Index Red Edge
- **MTCI** - MERIS Terrestrial Chlorophyll Index
- **S2REP** - Sentinel-2 Red Edge Position
- **NDCI** - Normalized Difference Chlorophyll Index (cyanobacteria/water quality) 🆕
- **CIG** - Chlorophyll Index Green

### 🌊 Sentinel-3 Exclusive (OLCI 21-band)
- **OCI** - OLCI Chlorophyll Index (optimized for S3)
- **TSI** - Trophic State Index (water quality assessment)
- **CDOM** - Colored Dissolved Organic Matter Index
- **Turbidity** - Water Turbidity Index (sediment monitoring)
- **SPM** - Suspended Particulate Matter Index
- **KD490** - Diffuse Attenuation Coefficient at 490nm
- **Floating Algae** - Floating Algae Index (bloom detection)
- **Red Edge Position** - OLCI-optimized red edge position
- **Fluorescence Height** - Chlorophyll fluorescence detection
- **Water Leaving Reflectance** - Aquatic reflectance analysis

### 🛰️ SAR Indices (Sentinel-1)
- **RVI** - Radar Vegetation Index (recommended for vegetation monitoring)
- **VV/VH Ratio** - Polarization ratio (excellent for structural change detection)
- **VH** - Cross-polarization (sensitive to volume scattering from vegetation)
- **VV** - Co-polarization (sensitive to surface roughness)
- **DPSVI** - Dual-pol SAR Vegetation Index (optimized for dense vegetation)
- **RFDI** - Radar Forest Degradation Index (deforestation monitoring) 🆕
- **VSDI** - Vegetation Scattering Diversity Index (structural diversity) 🆕

## 🧠 Machine Learning Classification (NEW in v0.6.0)

### Supervised Classification Algorithms
- **Random Forest** - With feature importance analysis
- **Support Vector Machine (SVM)** - For complex decision boundaries
- **CART** - Classification and Regression Trees
- **Naive Bayes** - Probabilistic classification
- **Gradient Tree Boost** - Advanced ensemble method

### Unsupervised Clustering
- **K-means** - Classic clustering algorithm
- **Cascade K-means** - Hierarchical clustering approach
- **LDA** - Latent Dirichlet Allocation for pattern discovery

### Classification Workflow
```python
from ndvi2gif import NdviSeasonality, LandCoverClassifier

# Create multi-temporal features
processor = NdviSeasonality(
    roi='study_area.shp',
    sat='S2',
    periods=12,
    start_year=2022,
    end_year=2024
)

# Initialize classifier
classifier = LandCoverClassifier(processor)

# Create feature stack with multiple indices
features = classifier.create_feature_stack(
    indices=['ndvi', 'evi', 'ndwi', 'ndre'],
    include_statistics=True,
    normalize=True
)

# Add training data
classifier.add_training_data('training_samples.shp')

# Classify with Random Forest
classification = classifier.classify_supervised('random_forest')

# Assess accuracy
classifier.assess_accuracy()

# Export results
processor.export_to_drive(
    image=classification,
    description="landcover_2024",
    folder="classifications"
)
```

## Supported Datasets

**Sentinel:**

- **[Sentinel-1 (SAR)](https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S1_GRD)** - Enhanced with dual polarization (VV+VH), speckle filtering, and orbit control
- **[Sentinel-2 (Surface Reflectance)](https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S2_SR_HARMONIZED)** - High resolution optical imagery with Red Edge bands
- **[Sentinel-3 OLCI (Level-1B TOA)](https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S3_OLCI)** - 21-band ocean and land color instrument with daily global coverage 🆕

**Landsat (Surface Reflectance):**

- [Landsat 4 TM](https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LT04_C02_T1_L2)
- [Landsat 5 TM](https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LT05_C02_T1_L2)
- [Landsat 7 ETM+](https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LE07_C02_T1_L2)
- [Landsat 8 OLI](https://developers.google.com/earth-engine/datasets/catalog/landsat-8)
- [Landsat 9 OLI](https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LC09_C02_T1_L2)

**MODIS (Surface Reflectance):**

- [MOD09A1 (SR)](https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MOD09A1)

You can combine any of the supported indices, datasets, and statistical methods. By default, the tool uses NDVI with the **maximum** statistic to avoid cloud contamination. However, **median** and **custom percentiles** are often visually better for Landsat datasets and specific applications.

Note: **Sentinel-2** uses Surface Reflectance, **Sentinel-3** uses Level-1B TOA radiance (optimized for aquatic applications), while **Landsat and MODIS** use Surface Reflectance (SR) for superior atmospheric correction and scientific quality.

The tool generates rasters with 4 (seasons), 12 (months), or 24 (custom periods) bands per year.

Beyond creating a nice-looking animated GIF, this multi-seasonal compositing method provides insights into vegetation dynamics, phenology, land cover, and more. High values in all seasons (white tones) typically mean perennial vegetation, while low values (dark tones) might represent water, soil, or impervious surfaces.

## GeoTIFF Export

You can also export seasonal composites as GeoTIFF files for further analysis. Multi-year composites are supported as well. For example, you can export median NDVI per season for all of Africa between 2001–2020, bi-monthly VV/VH ratios for crop monitoring, or daily Sentinel-3 turbidity indices for water quality assessment.

## 📤 Enhanced Export Capabilities (NEW in v0.6.0)

### Google Drive Export
Export any image or classification directly to your Google Drive:

```python
processor.export_to_drive(
    image=classified_map,
    description="landcover_2023",
    folder="ndvi2gif_results",
    scale=30,
    crs="EPSG:4326",
    maxPixels=1e13
)
```

### Earth Engine Asset Export
Save results as Earth Engine Assets for further processing:

```python
processor.export_to_asset(
    image=classification,
    asset_id="users/yourname/landcover_2023",
    pyramiding_policy={"class": "mode"},
    overwrite=True
)
```

### Automatic Scale Detection
The library now automatically selects the appropriate scale based on the sensor:
- Sentinel-2: 10m
- Sentinel-3: 300m
- Landsat: 30m
- MODIS: 500m
- Sentinel-1: 10m

---

## Installation

You can install `ndvi2gif` using either **pip** or **conda**:

### Using pip:

```bash
pip install ndvi2gif
```

### Using conda:

```bash
conda install -c conda-forge ndvi2gif
```

## 📚 Examples & Tutorials

Check out our comprehensive examples:

- **[Comprehensive Example](https://github.com/Digdgeo/Ndvi2Gif/blob/master/examples_notebooks/ndvi2gif%20extended%20version.ipynb)** - Complete guide to all ndvi2gif features
- **[Input Types Guide](https://github.com/Digdgeo/Ndvi2Gif/blob/master/examples_notebooks/NDVI2Gif_InputsTypes.ipynb)** - Different ways to specify your region of interest

*More examples are regularly added to showcase new capabilities and use cases.*

## Quick Usage Example

```python
import ee
import geemap
from ndvi2gif import NdviSeasonality, TimeSeriesAnalyzer, LandCoverClassifier

# Authenticate Earth Engine
ee.Authenticate()
ee.Initialize()

# Basic NDVI analysis
ndvi_analysis = NdviSeasonality(
    roi=your_roi,           # Your region of interest
    periods=12,             # Monthly analysis
    start_year=2023,
    end_year=2024,
    sat='S2',               # Sentinel-2
    key='percentile',       # Use percentile statistic
    percentile=85,          # 85th percentile (flexible!)
    index='ndvi'
)

# Generate composite
composite = ndvi_analysis.get_year_composite()

# Create animated GIF
ndvi_analysis.get_gif(name='ndvi_evolution.gif')

# NEW in v0.6.0: Land Cover Classification
classifier = LandCoverClassifier(ndvi_analysis)

# Create multi-index feature stack
features = classifier.create_feature_stack(
    indices=['ndvi', 'evi', 'ndwi'],
    normalize=True
)

# Add training data and classify
classifier.add_training_data('training_points.shp')
landcover = classifier.classify_supervised('random_forest')

# Export to Drive
ndvi_analysis.export_to_drive(
    image=landcover,
    description="classification_2024",
    folder="results"
)

# Get feature importance
importance = classifier.get_feature_importance()
print(f"Most important features: {importance[:5]}")

# NEW: Sentinel-3 water quality monitoring
water_quality = NdviSeasonality(
    roi=your_lake,
    periods=24,             # Bi-monthly for detailed monitoring  
    start_year=2023,
    end_year=2024,
    sat='S3',               # Sentinel-3 OLCI
    key='median',
    index='turbidity'       # Water turbidity assessment
)

# NEW: Daily algal bloom detection
algae_monitor = NdviSeasonality(
    roi=your_water_body,
    periods=12,             # Monthly analysis
    sat='S3',               # Daily coverage with S3
    index='floating_algae', # Specialized for bloom detection
    key='mean',
    start_year=2024,
    end_year=2024
)

# Advanced: Sentinel-2 Red Edge analysis for precision agriculture
chlorophyll_analysis = NdviSeasonality(
    roi=your_agricultural_field,
    periods=24,             # Bi-monthly for detailed monitoring
    sat='S2',               # Only S2 has Red Edge bands
    index='ireci',          # Highly sensitive to chlorophyll
    key='median',
    start_year=2023,
    end_year=2024
)

# SAR-based crop monitoring with orbit control
sar_analysis = NdviSeasonality(
    roi=your_roi,
    periods=24,             # Bi-monthly for detailed monitoring
    start_year=2023,
    end_year=2024,
    sat='S1',               # Sentinel-1 SAR
    key='mean',
    index='vv_vh_ratio',    # Excellent for harvest detection
    orbit='DESCENDING'      # Use only descending orbits for consistency
)

# Cyanobacteria detection with NDCI
cyano_detection = NdviSeasonality(
    roi=your_lake,
    periods=12,             # Monthly monitoring
    sat='S2',               # NDCI requires Red Edge
    index='ndci',           # Cyanobacteria detection
    key='percentile',
    percentile=75,
    start_year=2023,
    end_year=2024
)

#### TimeSeriesAnalyzer – trend and phenology ####
# Seasonal NDVI composites
ndvi = NdviSeasonality(
    roi=your_roi,
    sat='S2',
    periods=12,   # monthly
    start_year=2018,
    end_year=2024,
    index='ndvi'
)

# Analyze temporal trends and phenology
ts = TimeSeriesAnalyzer(ndvi)
df = ts.extract_time_series()
trend = ts.analyze_trend(df)
ts.plot_comprehensive_analysis()

#### SAR Analysis ####

from ndvi2gif import S1ARDProcessor
import ee

ee.Initialize()

# Configure ARD processor with terrain correction + Refined Lee filter
s1 = S1ARDProcessor(
    speckle_filter='REFINED_LEE',
    terrain_correction=True,
    terrain_flattening_model='VOLUME',
    dem='COPERNICUS_30'
)

# Apply corrections to a Sentinel-1 image
image = ee.Image("COPERNICUS/S1_GRD/...")  # replace with your image ID
processed = s1.apply_speckle_filter(s1.apply_terrain_correction(image))


For complete examples, see the [example notebooks](examples_notebooks/) folder.
```
---

## Use Cases

**🤖 Land Cover Classification (NEW in v0.6.0)**
- Multi-temporal crop type mapping
- Urban expansion monitoring
- Forest change detection
- Wetland classification
- Feature importance analysis for ecological studies

**🌾 Agricultural Monitoring**
- Crop phenology tracking with optical indices
- Crop type classification with Random Forest
- Harvest timing detection with SAR VV/VH ratios
- Irrigation monitoring with NDWI
- Yield prediction with multi-temporal NDVI
- Precision agriculture with Red Edge indices (S2 exclusive)

**🌊 Water Quality & Environmental Monitoring**
- Daily algal bloom detection with Sentinel-3
- Cyanobacteria monitoring with NDCI (S2 Red Edge)
- Lake and coastal water quality assessment
- Turbidity and sediment tracking
- Harmful algal bloom early warning systems

**🌍 Environmental Research**
- Drought assessment with flexible percentile analysis
- Vegetation change detection combining optical and SAR
- Snow cover analysis with NDSI
- Multi-sensor ecosystem monitoring
- Land cover change with machine learning classification

**📊 Operational Applications**
- Real-time monitoring with incomplete year support
- Multi-year trend analysis for climate studies
- Automated reporting with GeoTIFF exports
- Quality assessment with robust statistics
- Geometric consistency with SAR orbit control
- Direct export to Google Drive and Earth Engine Assets
- Generate reference rasters for pseudo-invariant feature normalization ([ProtocoloV2](https://github.com/Digdgeo/ProtocoloV2))

## Roadmap 🗺️ 

**v0.6.0 ✅ Machine Learning & Classification Suite**  
Status: **Released August 2025!**

✅ **LandCoverClassifier** – Complete classification workflows  
✅ **Multiple ML Algorithms** – RF, SVM, CART, K-means, and more  
✅ **Enhanced Exports** – Google Drive and EE Assets  
✅ **95%+ Documentation** – Comprehensive docstrings with examples  
✅ **Feature Engineering** – Multi-temporal stacks with normalization  

**v1.0.0 🎯 Complete Climate Analysis Platform**  
Status: **Planned**

📚 **Jupyter Book** – Interactive documentation and tutorials  
🌡️ **Climate Datasets** – ERA5, CHIRPS, TerraClimate integration  
🌍 **Climate Analysis** – Advanced climate change assessment tools  
🔧 **API Stability** – Long-term support commitment  

## Contributing

We welcome contributions from the community! Whether you're a developer, researcher, or just curious about remote sensing, your input can help improve Ndvi2Gif.

🐛 **Bug reports**: [GitHub Issues](https://github.com/Digdgeo/Ndvi2Gif/issues)

💡 **Feature requests**: [GitHub Discussions](https://github.com/Digdgeo/Ndvi2Gif/discussions)

🤝 **Pull requests**: Always welcome!

📚 **Example contributions**: Share your use cases in the `examples_notebooks/` folder

---

## 📖 Citation

JOSS Manuscript in preparation. For now, please cite this software as:

```bibtex
@software{garcia_diaz_ndvi2gif_2024,
  author = {García Díaz, Diego},
  title = {ndvi2gif: Multi-Seasonal Remote Sensing Analysis Suite},
  url = {https://github.com/Digdgeo/Ndvi2Gif},
  version = {0.6.0},
  year = {2025}
}

## Project Statistics

- **Current Version:** 0.6.0
- **Supported Sensors:** 5 (S1, S2, S3, Landsat, MODIS)
- **Available Indices:** 40+
- **ML Algorithms:** 8 (5 supervised, 3 unsupervised)
- **Lines of Code:** ~4,800
- **Documentation Coverage:** 95%+
- **Test Coverage:** Growing with each release

## License

This project is licensed under the MIT License - see the [LICENSE.txt](LICENSE.txt) file for details.

## Acknowledgments

- Built on [Google Earth Engine](https://earthengine.google.com/) and [geemap](https://geemap.org/)
- Special thanks to Qiusheng Wu and to the Google Earth Engine team and the open-source remote sensing community
