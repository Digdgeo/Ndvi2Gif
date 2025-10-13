# Ndvi2Gif Tutorial

Welcome to the comprehensive tutorial for **Ndvi2Gif**, a Python library for multi-seasonal remote sensing analysis with Google Earth Engine.

```{image} https://i.imgur.com/Y5dOWIk.jpeg
:alt: Ndvi2Gif Banner
:width: 800px
:align: center
```

## What is Ndvi2Gif?

**Ndvi2Gif** is a powerful remote sensing analytics suite that simplifies access to global satellite data through Google Earth Engine. While its name highlights the ability to create seasonal GIF animations, the true power of this tool lies in its capability to compute and export pixel-wise statistics for any region on Earth, across any time span covered by supported remote sensing datasets.

### Key Features

- 🛰️ **Multi-Sensor Support**: Sentinel-1/2/3, Landsat 4-9, MODIS
- 📊 **40+ Indices**: Vegetation, water, SAR, and specialized indices
- 🤖 **Machine Learning**: Supervised and unsupervised land cover classification
- 📈 **Time Series Analysis**: Trend detection, phenology metrics
- 🌍 **Global Coverage**: Process any region on Earth
- 📤 **Flexible Export**: GeoTIFF, Google Drive, Earth Engine Assets
- 🎨 **Visualization**: Automated GIF generation

### What You'll Learn

This tutorial will guide you through:

1. **Getting Started**: Installation, authentication, and your first analysis
2. **Core Tutorials**: Step-by-step guides for common workflows
3. **Advanced Features**: SAR processing, time series, classification
4. **Use Cases**: Real-world applications in agriculture, wetlands, drought assessment
5. **Reference**: Complete API documentation and indices catalog

## Quick Example

```python
import ee
from ndvi2gif import NdviSeasonality

# Authenticate and initialize
ee.Initialize()

# Create seasonal NDVI composite
ndvi = NdviSeasonality(
    roi='your_area.shp',
    sat='S2',
    periods=12,
    start_year=2023,
    end_year=2024,
    index='ndvi'
)

# Generate animated GIF
ndvi.get_gif('ndvi_evolution.gif')
```

## Who Should Use This Tutorial?

- 🌾 **Agricultural researchers** monitoring crop phenology
- 🌊 **Environmental scientists** assessing water quality
- 🌳 **Ecologists** studying vegetation dynamics
- 🛰️ **Remote sensing analysts** working with multi-temporal data
- 🎓 **Students** learning satellite image analysis

## Prerequisites

- Basic Python knowledge
- Familiarity with geospatial concepts (optional but helpful)
- Google Earth Engine account (free for research and education)

## Support & Community

- 📖 [Documentation](https://ndvi2gif.readthedocs.io)
- 💻 [GitHub Repository](https://github.com/Digdgeo/Ndvi2Gif)
- 🐛 [Issue Tracker](https://github.com/Digdgeo/Ndvi2Gif/issues)
- 💬 [Discussions](https://github.com/Digdgeo/Ndvi2Gif/discussions)

## Citation

If you use Ndvi2Gif in your research, please cite:

```bibtex
@software{garcia_diaz_ndvi2gif_2025,
  author = {García Díaz, Diego},
  title = {ndvi2gif: Multi-Seasonal Remote Sensing Analysis Suite},
  url = {https://github.com/Digdgeo/Ndvi2Gif},
  version = {0.6.0},
  year = {2025}
}
```

---

Let's get started! 🚀

```{tableofcontents}
```

