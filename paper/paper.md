---
title: 'Ndvi2Gif: A Python Package for Multi-Seasonal Remote Sensing Analysis via Google Earth Engine'
tags:
  - Python
  - remote sensing
  - Google Earth Engine
  - time series analysis
  - vegetation monitoring
  - SAR processing
  - phenology
  - NDVI
authors:
  - name: Diego García Díaz
    orcid: 0000-0002-2757-7112
    affiliation: 1
affiliations:
 - name: Estación Biológica de Doñana (EBD-CSIC), Spain
   index: 1
date: 28 December 2025
bibliography: paper.bib
---

# Summary

Google Earth Engine (GEE) is a cloud-based platform for planetary-scale geospatial analysis [@Gorelick2017]. Ndvi2Gif is a Python package built on top of GEE and geemap [@Wu2020] that simplifies the generation and analysis of multi-temporal composite images from satellite and climate reanalysis data.

The core concept of Ndvi2Gif is the creation of **multi-seasonal statistical composites**. The `get_year_composite()` method returns an `ee.ImageCollection` where each image represents one year, with bands for each temporal period (e.g., `periods=12` generates bands named `january`, `february`, ..., `december`). This band-based temporal organization enables queries such as: *What was the maximum flood extent detected in January across 40 years of Landsat (1984–2024)?* by simply calling `collection.select('january').max()`. Built-in zonal statistics extraction via `get_stats()` allows direct computation of statistics per polygon without external GIS processing.

# Statement of need

Generating multi-temporal composites in GEE requires extensive code for collection filtering, cloud masking, scaling to reflectance, and temporal aggregation. The following example shows the standard approach to create seasonal NDVI composites from Landsat:

```python
import ee
ee.Initialize()

roi = ee.Geometry.Point([-5.9, 37.2]).buffer(10000)

# Scale Landsat 8-9 (OLI) to surface reflectance
def scale_OLI(image):
    optical = image.select(['SR_B2','SR_B3','SR_B4','SR_B5','SR_B6','SR_B7']) \
                   .multiply(0.0000275).add(-0.2) \
                   .rename(['Blue','Green','Red','Nir','Swir1','Swir2'])
    return image.addBands(optical, None, True)

# Scale Landsat 4-7 (TM/ETM+) to surface reflectance  
def scale_ETM(image):
    optical = image.select(['SR_B1','SR_B2','SR_B3','SR_B4','SR_B5','SR_B7']) \
                   .multiply(0.0000275).add(-0.2) \
                   .rename(['Blue','Green','Red','Nir','Swir1','Swir2'])
    return image.addBands(optical, None, True)

# Cloud mask using QA_PIXEL band
def maskClouds(image):
    qa = image.select('QA_PIXEL')
    mask = qa.bitwiseAnd(1 << 3).eq(0).And(qa.bitwiseAnd(1 << 4).eq(0))
    return image.updateMask(mask)

def addNDVI(image):
    return image.addBands(
        image.normalizedDifference(['Nir', 'Red']).rename('NDVI'))

def getSeasonalComposite(startDate, endDate):
    l8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2') \
           .filterBounds(roi).filterDate(startDate, endDate) \
           .map(maskClouds).map(scale_OLI).map(addNDVI)
    l5 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2') \
           .filterBounds(roi).filterDate(startDate, endDate) \
           .map(maskClouds).map(scale_ETM).map(addNDVI)
    return l8.merge(l5).select('NDVI').max()

composite = ee.Image.cat([
    getSeasonalComposite('2013-01-01', '2013-03-31'),
    getSeasonalComposite('2013-04-01', '2013-06-30'),
    getSeasonalComposite('2013-07-01', '2013-09-30'),
    getSeasonalComposite('2013-10-01', '2013-12-31')
])
```

The 40 lines above can be reduced to 5 lines using Ndvi2Gif:

```python
from ndvi2gif import NdviSeasonality

ndvi = NdviSeasonality(roi=roi, start_year=2013, end_year=2013, 
                       sat='Landsat', index='ndvi', periods=4)
collection = ndvi.get_year_composite()  # ImageCollection: 1 image × 4 bands
```

This simplification extends across 7 satellite platforms, 88 spectral and climate variables, multiple temporal aggregation methods, and flexible region of interest specifications—all through a unified API that inherits geemap's interactive mapping capabilities. The architectural decisions and design philosophy that enable this simplification, along with how Ndvi2Gif differentiates from other GEE Python packages, are detailed in the Software Design section below. \autoref{fig:workflow} provides an overview of the package's workflow and capabilities.

![Ndvi2Gif workflow architecture and capabilities. Color coding: **Blue** — input data sources and ROI specification; **Green** — processing operations (reducers and spectral indices); **Purple** — main analysis classes (NdviSeasonality for temporal aggregation, TimeSeriesAnalyzer, LandCoverClassifier, S1ARDProcessor); **Orange** — final outputs.\label{fig:workflow}](ndvi2gif_workflow.svg)

# Key features

## Multi-platform data access

Ndvi2Gif provides unified access to optical, SAR, and climate reanalysis datasets through a consistent interface (\autoref{tab:datasets}).

: Supported satellite and reanalysis datasets. \label{tab:datasets}

| Platform | Collection | Resolution | Coverage | Variables |
|:---------|:-----------|:-----------|:---------|:----------|
| Sentinel-2 | S2 SR HARMONIZED | 10/20m | 2015–present | 40+ optical |
| Sentinel-3 | OLCI | 300m | 2016–present | 10 water quality |
| Landsat | L4-L9 Col. 2 | 30m | 1982–present | 40+ optical |
| MODIS | Terra/Aqua | 500m | 2000–present | 40+ optical |
| Sentinel-1 | GRD ARD | 10m | 2014–present | 7 SAR |
| ERA5-Land | Hourly/Monthly | 11km | 1950–present | 47 climate |
| CHIRPS | Daily | 5.5km | 1981–present | Precipitation |

## Flexible region of interest

The library supports multiple ROI specification formats, enabling seamless integration with existing workflows and the eLTER research network (\autoref{tab:roi}).

: Supported region of interest formats. \label{tab:roi}

| Format | Example | Description |
|:-------|:--------|:------------|
| Shapefile | `'study_area.shp'` | Local vector file |
| GeoJSON | `'boundaries.geojson'` | GeoJSON file |
| DEIMS ID | `'deimsid:8eda49e9...'` | eLTER site [@Wohner2019] |
| S2 tile | `'s2:29SQB'` | MGRS tile code |
| Landsat WRS | `'wrs:200,32'` | Path and Row |
| ee.Geometry | `ee.Geometry.Point()` | GEE geometry |
| geemap | `Map.draw_features` | Interactive selection |

## Advanced SAR processing

The `S1ARDProcessor` module implements a complete Analysis Ready Data (ARD) pipeline for Sentinel-1 imagery, based on Vollrath et al. [-@Vollrath2020]. This addresses a critical gap in GEE's standard Sentinel-1 preprocessing:

**Radiometric Terrain Correction**: Angular-based correction with configurable DEMs (Copernicus, SRTM).

**Speckle Filtering**: Five algorithms including Refined Lee, Gamma MAP, and Lee Sigma.

**Scattering Models**: Volume (vegetation) and Surface (bare soil, water) models.

This preprocessing is essential for accurate SAR-based vegetation monitoring in mountainous or heterogeneous terrain.

## Time series and phenology analysis

The `TimeSeriesAnalyzer` module extracts temporal patterns from composite stacks:

- **Trend detection**: Mann-Kendall test, Sen's slope, linear regression with significance testing
- **Phenology metrics**: Start/End/Peak of Season (SOS, EOS, POS), season length, growth and senescence rates
- **Visualization**: Multi-panel dashboards for comprehensive temporal analysis

## Machine learning classification

The `LandCoverClassifier` module provides supervised (Random Forest, SVM, CART, Gradient Tree Boost) and unsupervised (K-means, Cascade K-means) classification with accuracy assessment and feature importance analysis.

# Software Design

Several Python packages extend GEE for remote sensing workflows: eemont [@Montero2021] provides preprocessing and spectral indices computation, while wxee [@Zuspan2021] enables temporal aggregation and xarray export. These tools can be combined for similar tasks. However, Ndvi2Gif (first released in May 2020) adopts a fundamentally different temporal organization strategy. While eemont maintains the native temporal structure of `ee.ImageCollection` objects and wxee aggregates time series into sequential steps (e.g., 480 monthly composites for 40 years), ndvi2gif implements a hierarchical band-based approach where **each image represents a temporal unit (typically a year) with bands encoding sub-periods** (e.g., 40 images/years with 12, 24, or 52 bands each for monthly, biweekly, or weekly composites).

This architectural decision emerged from addressing specific research questions in ecological monitoring that require direct multi-decadal comparisons across fixed seasonal windows. The band-based organization enables queries such as: *What is the mean cyanobacteria index (NDCI) detected from Sentinel-2 in Andalusian reservoirs during August across the last 10 years?* by simply calling `collection.select('august').mean()`, or *In which week of the year does chlorophyll (OCI) detected from Sentinel-3 peak in European lakes?* by computing `collection.max().bandNames()` to identify the peak band—rather than requiring temporal filtering, grouping, and reduction operations across hundreds of separate images.

Beyond this organizational approach, ndvi2gif provides integrated capabilities unavailable in other GEE Python packages: complete Sentinel-1 Analysis Ready Data (ARD) preprocessing based on Vollrath et al. [-@Vollrath2020], automated phenology metrics extraction (Start/End/Peak of Season with trend analysis), machine learning classification, zonal statistics via `get_stats()`, and geemap-based visualization—all within a unified API built around the band-organized composite structure. The package has been applied in calibration and validation studies of Earth Observation products, including phenology and surface water monitoring at the Doñana protected area [@DiazDelgado2024].

Contributing these capabilities to existing packages would require fundamental restructuring of their core temporal handling mechanisms and would conflict with their design philosophies. Eemont focuses on preprocessing and spectral indices computation while preserving native GEE temporal structures, whereas wxee prioritizes xarray integration and sequential time series analysis. The band-based temporal organization in ndvi2gif represents a cohesive design philosophy that addresses a distinct set of use cases in multi-seasonal ecological analysis that would not naturally fit within the existing packages' architectures.

# Research Impact Statement

Ndvi2Gif has demonstrated realized research impact since its initial release in May 2020. The package has been applied in peer-reviewed scientific research, including calibration and validation studies of Earth Observation products for phenology and surface water monitoring at the Doñana protected area [@DiazDelgado2024]. It serves as a foundational component of geeltermap, a Python mapping application that provides operational environmental monitoring tools for the eLTER network [@DiazDelgado2024b]. Geeltermap integrates three applications: PhenoApp for phenology monitoring [@Garcia2023], FloodApp for flood extent analysis, and LSTApp for land surface temperature assessment—all built upon ndvi2gif's temporal composite framework and geemap's interactive mapping capabilities. This integration enables standardized multi-temporal analysis across European research sites through the DEIMS-SDR registry.

The package currently serves as a core analytical tool in two ongoing research projects: the eLTER network infrastructure for long-term ecosystem research and the SUMHAL biodiversity monitoring initiative in Spanish protected areas. It addresses a specific gap in the GEE Python ecosystem by enabling direct multi-decadal queries across fixed seasonal windows through band-based temporal organization, while integrating complete Sentinel-1 ARD preprocessing, automated phenology metrics extraction, and machine learning classification within a unified API. Its application in operational environmental monitoring workflows demonstrates readiness for research use and capacity to enable reproducible multi-temporal remote sensing analysis.

# AI Usage Disclosure

Ndvi2Gif was initially developed in May 2020 and underwent sustained development through 2024 before current-generation AI coding assistants became widely available. The core scientific methodology, architectural design decisions, and algorithmic implementations represent the author's original intellectual contributions. Nowadays, the author uses AI tools frequently in daily development work for code refactoring, documentation generation, debugging support, and test implementation. This AI assistance functions as a productivity tool—analogous to an advanced code editor or linting system—while all strategic decisions regarding software architecture, scientific approach, API design, and feature prioritization remain entirely under the author's direction and expertise. The use of AI tools has not altered the fundamental nature of the software or its scientific contributions.

# Acknowledgements

The author thanks the Google Earth Engine team and community and acknowledges Qiusheng Wu for geemap, which provides the foundation for Ndvi2Gif's interactive capabilities, and for his extensive educational materials and tutorials on Google Earth Engine that were instrumental in the development of this package. The SAR ARD module is based on the work of Vollrath et al. and the gee_s1_ard implementation.

# References
