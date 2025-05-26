# Ndvi2Gif: Multi-Seasonal Remote Sensing Index Composites

[![PyPI version](https://img.shields.io/pypi/v/ndvi2gif.svg)](https://pypi.org/project/ndvi2gif/)
[![PyPI downloads](https://img.shields.io/pypi/dm/ndvi2gif.svg)](https://pypi.org/project/ndvi2gif/)
[![Conda version](https://img.shields.io/conda/vn/conda-forge/ndvi2gif.svg)](https://anaconda.org/conda-forge/ndvi2gif)
[![Conda downloads](https://img.shields.io/conda/dn/conda-forge/ndvi2gif.svg)](https://anaconda.org/conda-forge/ndvi2gif)
[![Build status](https://github.com/Digdgeo/Ndvi2Gif/actions/workflows/python-publish.yml/badge.svg)](https://github.com/Digdgeo/Ndvi2Gif/actions/workflows/python-publish.yml)

![NDVI2GIF Köln](https://i.imgur.com/Y5dOWIk.jpeg)
*Richter's stained glass in Cologne Cathedral. Inspiration for this library.*

**Ndvi2Gif** is a Python library designed to simplify access to global satellite data through the Google Earth Engine platform. While its name highlights the ability to create seasonal GIF animations, the true power of this tool lies in its capability to compute and export pixel-wise statistics for any regioque aun hace gifs
n on Earth, across any time span covered by supported remote sensing datasets.

Built on top of [Google Earth Engine](https://github.com/google/earthengine-api) and [Geemap](https://github.com/giswqs/geemap), it allows you to:

- Generate annual or multi-annual composited rasters (e.g., median NDVI per season between 2001 and 2020),
- Apply multiple statistics (mean, max, percentiles) across space and time,
- Export results as GeoTIFFs for further analysis,
- Retrieve zonal statistics over user-defined geometries,
- And yes — also create colorful GIFs for easy visualization.

Whether you're monitoring crop phenology, assessing drought trends, or preparing input layers for further ecological modeling, `ndvi2gif` makes it easier to extract reliable, multi-temporal remote sensing information at scale.

Ndvi2Gif was updated and extended as part of its integration into the eLTER and SUMHAL projects, which also enabled the use of eLTER site boundaries (via `deimsPy`) as one of its input sources.


![Interface Screenshot](https://i.imgur.com/Sv9LfYj.png)

## Why use Ndvi2Gif?

Unlike many visualization-oriented tools, Ndvi2Gif is designed as a **remote sensing analytics helper** that abstracts much of the complexity of working directly with Google Earth Engine.

You can:

- **Access pixel-wise statistics** over any Earth location, at any scale and time span.
  - Example: *Obtain the monthly median of the 95th NDVI percentile per pixel from 1984 to 2024 using Landsat data.*
  - Example: *Calculate the maximum of the seasonal NDWI maximums between 2017 and 2023 using Sentinel-2.*
- **Perform nested aggregations**:
  - First compute temporal summaries (e.g., per-season percentiles or means), then apply a second statistical reduction across years (e.g., median, min, max).
- **Target any ecological or phenological metric** by choosing the appropriate index and statistical pipeline.
- **Work globally**, without needing to download or preprocess raw satellite data — all computations are handled via Earth Engine's cloud infrastructure.

In other words: if you can describe a temporal range, a spatial region, an index, and a chain of statistics — `ndvi2gif` can likely generate it.

Yes, it makes nice GIFs — but it’s much more than that.
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

- Maximum
- Mean
- Median
- Percentile 90
- Percentile 95

## Available Indices

- NDVI
- EVI
- GNDVI
- SAVI
- NDWI
- MNDWI
- AEWI
- AEWINSH
- NDSI
- NBRI
- NDMI

## Supported Datasets

**Sentinel:**

- [Sentinel 1 (SAR)](https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S1_GRD)
- [Sentinel 2 (TOA)](https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S2_HARMONIZED)

**Landsat:**

- [Landsat 4 TM](https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LT04_C02_T1_L2)
- [Landsat 5 TM](https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LT05_C02_T1_L2)
- [Landsat 7 ETM+](https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LE07_C02_T1_L2)
- [Landsat 8 OLI](https://developers.google.com/earth-engine/datasets/catalog/landsat-8)
- [Landsat 9 OLI](https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LC09_C02_T1_L2)

**MODIS:**

- [MOD09A1 (SR)](https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MOD09A1)

You can combine any of the supported indices, datasets, and statistical methods. By default, the tool uses NDVI with the **maximum** statistic to avoid cloud contamination. However, **median** and **percentile 90** are often visually better for Landsat datasets.

Note: Sentinel-2 uses **TOA reflectance** (Surface Reflectance is only available since 2017), while Landsat and MODIS collections use **Surface Reflectance (SR)**.

The tool generates rasters with 4 (seasons), 12 (months), or 24 (custom periods) bands per year.

Beyond creating a nice-looking animated GIF, this multi-seasonal compositing method provides insights into vegetation dynamics, phenology, land cover, and more. High values in all seasons (white tones) typically mean perennial vegetation, while low values (dark tones) might represent water, soil, or impervious surfaces.

With SAR (Sentinel-1) support added, the tool now enables structural monitoring of vegetation using radar backscatter intensity.

## GeoTIFF Export

You can also export seasonal NDVI composites as GeoTIFF files for further analysis. Multi-year composites are supported as well. For example, you can export median NDVI per season for all of Africa between 2001–2020.

---

## Installation

You can install `ndvi2gif` using either **pip** or **conda**:

### Using pip:

Just run:

```bash
pip install ndvi2gif
```

### Using conda:

```bash
conda install -c conda-forge ndvi2gif
```


## Usage Example

See the [ndvi2gif_extended_version notebook](https://github.com/Digdgeo/Ndvi2Gif/blob/master/ndvi2gif/ndvi2gif_extended_version.ipynb) for full usage examples and reproducible workflows.

---

## Roadmap 🗺️ 

**v0.2.0 🚀 Dynamic Periods**
Status: Ready for release

✅ Dynamic period generation - Support for any number of temporal periods

✅ Code optimization - Eliminated 90% of code duplication

✅ Flexible temporality - Easy configuration from 4 seasons to 52 weeks or custom periods

✅ Full backward compatibility - All existing code works without changes

✅ Enhanced extensibility - Adding new satellites and datasets is now trivial

Breaking changes: None

**v0.2.1 🎯 Custom Year Cycles**
Status: In planning

📅 Custom year start - start_doy parameter for non-calendar years

📅 Hydrological years - October 1st to September 30th cycles

📅 Academic years - September 1st to August 31st cycles

📅 Fiscal years - July 1st to June 30th cycles

📅 Cross-calendar logic - Intelligent handling of periods spanning multiple calendar years

📅 Smart period naming - Context-aware naming for custom year cycles

Use cases: Water management, agricultural monitoring, institutional reporting

**v0.3.0 🔮 High-Frequency Analysis**
Status: Future

🛰️ Sentinel-3 support - Daily global coverage for rapid change detection

📊 High-frequency periods - Weekly (52 periods) and daily analysis capabilities

🎨 Enhanced visualizations - Improved GIF generation and export options

⚡ Performance optimizations - Faster processing for large temporal datasets

🌐 Multi-sensor fusion - Combine data from multiple satellite platforms

Use cases: Disaster monitoring, urban growth tracking, real-time environmental assessment

## Contributing

We welcome contributions from the community! Whether you're a developer, researcher, or just curious about remote sensing, your input can help improve Ndvi2Gif.Have ideas for new features? Found a bug? We'd love to hear from you!

🐛 Bug reports: GitHub Issues

💡 Feature requests: GitHub Discussions

🤝 Pull requests: Always welcome!