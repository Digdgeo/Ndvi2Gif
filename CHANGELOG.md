# Changelog

All notable changes to the `ndvi2gif` package will be documented in this file.

---

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.6.0] - 2025-09-15

### 🧠 **MACHINE LEARNING & CLASSIFICATION RELEASE**

This release introduces comprehensive land cover classification capabilities and enhanced export functionality, positioning ndvi2gif as a complete remote sensing analysis suite. The library continues to mature toward v1.0.0 with advanced documentation and expanded analytical capabilities.

---

## ✨ **New Features**

### 🧠 **Land Cover Classification Module** (NEW)

- **NEW MODULE**: Complete `LandCoverClassifier` class for supervised and unsupervised classification
- **NEW**: Cloud masking options for Sentinel-2 and Landsat collections
- **NEW**: Multi-temporal feature stack generation with automatic normalization
- **NEW**: Support for multiple classification algorithms:
  - Random Forest (with feature importance)
  - Support Vector Machine (SVM)
  - Classification and Regression Trees (CART)
  - Naive Bayes
  - Gradient Tree Boost
- **NEW**: Unsupervised clustering algorithms:
  - K-means
  - Cascade K-means
  - Latent Dirichlet Allocation (LDA)
- **NEW**: Comprehensive accuracy assessment with confusion matrices
- **NEW**: Training data support from shapefiles, GeoJSON, and point/polygon sampling
- **NEW**: Feature importance analysis for Random Forest models
- **NEW**: Visualization tools for confusion matrices and accuracy reports

### 🚀 **Enhanced Export Capabilities**

- **NEW**: `export_to_drive()` - Batch export to Google Drive with full parameter control
- **NEW**: `export_to_asset()` - Export to Earth Engine Assets with pyramiding policies
- **NEW**: `_default_scale_for_sat()` - Automatic scale selection based on sensor
- **NEW**: Advanced export options including:
  - Custom pyramiding policies for classification data
  - Overwrite protection for assets
  - Format-specific options (compression, file per band)
  - Maximum pixel limits and CRS control

### 📊 **Time Series Analysis Enhancements** (Updated)

- **IMPROVED**: Enhanced documentation with complete examples
- **IMPROVED**: Better error handling and user feedback
- **IMPROVED**: More robust phenology extraction methods
- **IMPROVED**: Advanced visualization capabilities with publication-ready plots

---

## 🔧 **Major Improvements**

### 📚 **Documentation Overhaul**

- **IMPROVED**: Complete Sphinx-style docstrings for all classes and methods
- **IMPROVED**: Comprehensive parameter documentation with types and examples
- **IMPROVED**: Scientific references added to all spectral indices
- **IMPROVED**: Detailed usage examples in docstrings
- **IMPROVED**: Better error descriptions with suggested solutions
- **IMPROVED**: Cross-references between related methods

### 🛰️ **SAR Processing Enhancements**

- **IMPROVED**: Enhanced error handling in `S1ARDProcessor`
- **IMPROVED**: Better documentation for terrain correction parameters
- **IMPROVED**: More detailed method descriptions with scientific references
- **IMPROVED**: Improved parameter validation and user feedback

### 🌍 **API Consistency**

- **IMPROVED**: Consistent parameter naming across all modules
- **IMPROVED**: Standardized return types and error handling
- **IMPROVED**: Better integration between `NdviSeasonality` and new modules
- **IMPROVED**: More informative console output and progress tracking

---

## 🔄 **API Changes & Enhancements**

### 📦 **Module Structure**

```python
# NEW imports available in v0.6.0
from ndvi2gif import (
    NdviSeasonality,        # Core functionality (enhanced)
    S1ARDProcessor,         # SAR preprocessing (improved docs)
    TimeSeriesAnalyzer,     # Time series analysis (enhanced)
    SpatialTrendAnalyzer,   # Spatial analysis (enhanced)
    LandCoverClassifier,    # NEW: Classification workflows
)
```

### 🆕 **New Method Signatures**

```python
# NEW: Enhanced export methods
processor.export_to_drive(
    image=classified_map,
    description="landcover_2023",
    folder="ndvi2gif_results",
    scale=30,
    crs="EPSG:4326"
)

processor.export_to_asset(
    image=classification,
    asset_id="users/yourname/landcover_2023",
    pyramiding_policy={"class": "mode"},
    overwrite=True
)

# NEW: Classification workflow
classifier = LandCoverClassifier(processor)
features = classifier.create_feature_stack(
    indices=['ndvi', 'evi', 'ndwi'],
    include_statistics=True,
    normalize=True
)
classifier.add_training_data('training_points.shp')
result = classifier.classify_supervised('random_forest')
```

---

## 🛠️ **Under the Hood**

### 🔧 **Code Quality**

- **IMPROVED**: Consistent error handling with informative messages
- **IMPROVED**: Better type hints throughout the codebase
- **IMPROVED**: More robust parameter validation
- **IMPROVED**: Enhanced memory efficiency in large-area processing
- **IMPROVED**: Better handling of edge cases and invalid inputs

### 📈 **Performance**

- **OPTIMIZED**: Feature stack generation for classification
- **OPTIMIZED**: Memory usage in multi-temporal processing
- **OPTIMIZED**: Export operations with better chunking strategies

---

## 🐛 **Bug Fixes**

- **FIXED**: Improved error handling when no satellite data is available
- **FIXED**: Better validation of ROI inputs and coordinate systems
- **FIXED**: Enhanced handling of edge cases in temporal compositing
- **FIXED**: More robust processing of incomplete time series
- **FIXED**: Better handling of mixed sensor collections

---

## 📖 **Examples & Use Cases**

### 🌾 **Agricultural Monitoring**

```python
# Multi-temporal crop classification
processor = NdviSeasonality(
    roi='farm_boundaries.shp',
    sat='S2', periods=12,
    start_year=2022, end_year=2024
)

classifier = LandCoverClassifier(processor)
features = classifier.create_feature_stack(['ndvi', 'evi', 'ndre'])
classifier.add_training_data('crop_samples.shp')
crop_map = classifier.classify_supervised('random_forest')
```

### 🌊 **Water Quality Assessment**

```python
# Sentinel-3 water quality with export to Drive
processor = NdviSeasonality(
    roi='lake_boundary.shp',
    sat='S3', index='turbidity',
    periods=24  # Bi-monthly
)

composites = processor.get_year_composite()
processor.export_to_drive(
    image=composites.first(),
    description="lake_turbidity_2024",
    folder="water_quality"
)
```

### 🏔️ **SAR Forest Monitoring**

```python
# Advanced SAR processing with classification
processor = NdviSeasonality(
    sat='S1', index='rvi',
    use_sar_ard=True,
    sar_speckle_filter='REFINED_LEE',
    sar_terrain_correction=True
)

classifier = LandCoverClassifier(processor)
forest_map = classifier.classify_unsupervised('kmeans', n_clusters=5)
```

---

## ⚠️ **No Breaking Changes**

Full backward compatibility maintained with v0.5.x. All existing code continues to work unchanged.

---

## 🔄 **Migration Guide**

### From v0.5.x to v0.6.0

No breaking changes! All existing code will continue to work. New features are additive:

```python
# v0.5.x code continues to work unchanged
processor = NdviSeasonality(sat='S2', index='ndvi')
processor.get_gif('animation.gif')

# v0.6.0 adds new capabilities
classifier = LandCoverClassifier(processor)  # NEW
processor.export_to_drive(image, "export")   # NEW
```

---

## 🎯 **Future Roadmap**

### v1.0.0 (Planned) - Complete Climate Analysis Platform
- **📚 Jupyter Book**: Interactive documentation with comprehensive tutorials and examples
- **🌡️ Climate Datasets**: Integration with ERA5, CHIRPS, TerraClimate, and other climate model datasets
- **🌍 Climate Analysis**: Advanced tools for climate change impact assessment and adaptation planning

---

## 📊 **Statistics**

- **New classes**: 1 (`LandCoverClassifier`)
- **New methods**: 15+ (classification, enhanced exports, utilities)
- **Enhanced methods**: 20+ (improved documentation and error handling)
- **Lines of code**: ~3,500 → ~4,800 (+37%)
- **Documentation coverage**: 95%+ (comprehensive docstrings)

---

## 🙏 **Acknowledgments**

Special thanks to the Google Earth Engine team and the open-source remote sensing community for their continued support and feedback that made this release possible.

---

## 📚 **Documentation**

Complete documentation with tutorials available at: [GitHub Repository](https://github.com/Digdgeo/Ndvi2Gif)

---

**Full Changelog**: https://github.com/Digdgeo/Ndvi2Gif/compare/v0.5.0...v0.6.0

## [0.5.0] - 2025-08-28

### Added

- **🛰️ Sentinel-1 ARD Processor**: New `S1ARDProcessor` module for advanced SAR preprocessing:
  - Radiometric terrain correction (angular method, Vollrath et al. 2020).
  - Configurable speckle filters: Boxcar, Lee, Refined Lee, Gamma-MAP, Lee Sigma.
  - Flexible DEM options (Copernicus 30/90, SRTM 30/90).
- **📈 TimeSeriesAnalyzer**: New module for time series and phenological analysis:
  - Robust extraction of temporal profiles from points or polygons.
  - Trend analysis (Mann-Kendall, Linear regression, Sen’s slope).
  - Comprehensive dashboards (trend, seasonality, autocorrelation, quality).
  - Phenological metrics (SOS, EOS, POS, LOS, amplitude, growth/senescence rates).
- **🌱 NdviSeasonality improvements**:
  - Extended ROI handling: DEIMS sites, Sentinel-2 MGRS tiles, Landsat WRS path/row, shapefiles, GeoJSON.
  - Flexible temporal periods (4, 12, 24, or custom definitions).
  - Optional SAR normalization and enhanced orbit handling.
  - More robust sensor-index validation.

### Changed

- **Visualization**: Unified plotting style with Seaborn/Matplotlib, clearer layouts.
- **Documentation**: Updated examples covering SAR and time series analysis.

### Fixed
- More robust handling of null/NaN values in temporal extraction.
- Minor bug fixes in period generation and export routines.

## [0.4.1] - 2025-07-21

### Added

- Just fixing some bugs in Readme.md

## [0.4.0] - 2025-07-21

### Added

- **🛰️ Sentinel-3 OLCI Support**: Revolutionary addition with 21 spectral bands and daily global coverage
- **🌊 Advanced Water Quality Indices**: 10 specialized aquatic indices including OCI, TSI, CDOM, turbidity, SPM, KD490, floating algae detection
- **🔬 Enhanced Sentinel-2**: Complete Red Edge implementation with Surface Reflectance for superior data quality  
- **💧 Cyanobacteria Detection**: New NDCI index for harmful algal bloom monitoring and water quality assessment
- **⚙️ SAR Orbit Control**: Precise control over Sentinel-1 ascending/descending orbits for geometric consistency
- **🎯 40+ Specialized Indices**: Comprehensive coverage with intelligent sensor-index validation
- **📊 Professional Architecture**: Clean, extensible design with enhanced error handling and documentation

### New Sentinel-3 Indices

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

### New SAR Indices

- **RFDI** - Radar Forest Degradation Index (deforestation monitoring)
- **VSDI** - Vegetation Scattering Diversity Index (structural diversity)

### Enhanced Features

- **Intelligent Validation**: Smart index-sensor compatibility checking prevents invalid combinations
- **Orbit Parameter**: Fine control over Sentinel-1 orbit selection (BOTH/ASCENDING/DESCENDING)
- **Advanced Use Cases**: Support for pseudo-invariant area radiometric normalization workflows

### Changed

- **Sentinel-2 to Surface Reflectance**: Upgraded from TOA to Surface Reflectance for better scientific quality
- **Simplified Architecture**: Removed unnecessary complexity while maintaining full functionality
- **Enhanced Documentation**: Professional-grade docstrings and examples

### Technical Improvements

- Modular sensor setup with clean separation of concerns
- Comprehensive sensor-index mapping and validation
- Enhanced error messages for better user experience
- Support for advanced radiometric normalization workflows

---

## [0.3.0] - 2025-07-17

### Added

- **New SAR Indices**: RVI (Radar Vegetation Index), VV/VH ratio, VH, VV, DPSVI for Sentinel-1
- **Flexible Percentiles**: Support for any percentile value (1-99) instead of fixed 90/95
- **Enhanced Sentinel-1**: VV+VH dual polarization with speckle filtering
- **Robust ROI Handling**: Support for drawn features, lists of features, and improved geometry conversion
- **Incomplete Year Support**: Automatic detection and processing of available periods for current/incomplete years
- **Enhanced Dependencies**: Added pycrs and deims as core dependencies (now available in conda)
- **Example Notebooks**: Comprehensive examples in `examples_notebooks/` folder

### Fixed

- ROI conversion for drawn geometries and feature lists from geemap
- Speckle filter now preserves temporal properties (system:time_start)
- Band naming consistency for SAR indices
- Error handling for missing data periods
- Dependency issues with pycrs and deims

### Improved

- More robust error handling throughout the library
- Better documentation and examples
- Enhanced support for agricultural monitoring workflows
- Simplified installation process

## [0.2.0] - 2025-01-27

### Added

- **Dynamic period generation**: Support for any number of temporal periods (4, 6, 8, 12, 24, 52, or any custom number).
- **Flexible temporal analysis**: Easy configuration from traditional 4 seasons to 52 weekly periods or any custom division.
- **Enhanced extensibility**: Adding new satellites and datasets is now trivial with the unified architecture.

### Changed

- **Major code refactoring**: Eliminated over 90% of code duplication by replacing 40+ individual period functions with a single dynamic system.
- **Improved maintainability**: Reduced codebase from ~3,000 lines to ~400 lines while maintaining all functionality.
- **Enhanced performance**: Streamlined period generation and composite creation.

### Technical Details

- Replaced hardcoded period definitions with dynamic `_generate_periods()` method.
- Consolidated all `get_winter()`, `get_january()`, `get_p1()` through `get_p24()` functions into a single `get_period_composite()` method.
- Maintained full backward compatibility - all existing code works without changes.
- Added comprehensive leap year handling to prevent date-related errors.

### Breaking Changes

- None - this release maintains 100% backward compatibility.

---

## [0.1.5] - 2025-05-26

### Fixed

- Fixing bug with MNDWI index.

---

## [0.1.4] - 2025-05-25

### Fixed

- Nothing really changes, just a f* problem with release version management.

---

## [0.1.3] - 2025-05-25

### Fixed

- Nothing really changes, just a f* problem with release version management.

---

## [0.1.2] - 2025-05-25

### Added

- Complete rework and translation of the README into Markdown format.
- Included new seasonal/statistical methods and updated docstrings in English.
- Added support for region input via Sentinel-2 tiles and Landsat path/row.
- Added rich ROI input documentation with tabular summary.
- Added `deims` dependency as optional to avoid conda forge problems

### Changed

- Clarified the purpose of the library as a broader seasonal analysis tool, not just for GIF generation.
- Cleaned and validated `setup.cfg` and `pyproject.toml`.
- Added extra requirements group for `deims`.

---

## [0.1.1] - 2025-05-21

### Fixed

- Fixed rendering issue in `README.rst` that caused PyPI upload failure.
- Rebuilt and republished the package with correct long description format.

---

## [0.1.0] - 2025-05-21

### Added

- Compatibility with Conda packaging and `conda-forge` ecosystem.
- Included `MANIFEST.in` to ensure `LICENSE` and `README.rst` are bundled in source distribution.
- Improved `README.rst` formatting to comply with PyPI rendering rules.

### Changed

- Switched versioning to semantic 0.x.y style for future compatibility.
- Cleaned and validated metadata to allow upload to both PyPI and Conda Forge.

### Note

- This is a technical release — no changes to the core functionality.

---

## [0.0.9] - 2025-05-20

### Changed

- Version bump to align `setup.cfg`, PyPI and GitHub release.
- No functional changes from version 0.0.7.

---

## [0.0.7] - 2025-05-20

### Added

- New method `get_ndmi()` to compute the Normalized Difference Moisture Index (NDMI).
- New (old) method `get_gif()` to download a gif for the selected index/bands.
- Package structure modernized:
  - Added `setup.cfg` and `pyproject.toml` (PEP 517/518 compliant).
  - Optional removal of legacy `setup.py`.
- Updated dependencies:
  - `geemap` pinned to version `0.29.5`.
  - `numpy` constrained to `<2.0` for compatibility.

### Fixed

- Compatibility issues with recent versions of `geemap`, `xarray`, and `numpy`.
- Resolved import error caused by the removal of `np.unicode_` in NumPy 2.0.

---

## [0.0.6] - 2023-03-10

### Added

- Initial public release of the `ndvi2gif` package.
- Generate seasonal composites and extract statistical summaries from several remote sensing index using Google Earth Engine and geemap.
- Export to animated GIF and GeoTIFF format.
