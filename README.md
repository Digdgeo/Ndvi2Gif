# Ndvi2Gif: Multi-Seasonal Remote Sensing Index Composites

[![DOI](https://joss.theoj.org/papers/10.21105/joss.10654/status.svg)](https://doi.org/10.21105/joss.10654)
[![PyPI version](https://img.shields.io/pypi/v/ndvi2gif.svg)](https://pypi.org/project/ndvi2gif/)
[![PyPI downloads](https://img.shields.io/pypi/dm/ndvi2gif.svg)](https://pypi.org/project/ndvi2gif/)
[![Conda version](https://img.shields.io/conda/vn/conda-forge/ndvi2gif.svg)](https://anaconda.org/conda-forge/ndvi2gif)
[![Conda downloads](https://img.shields.io/conda/dn/conda-forge/ndvi2gif.svg)](https://anaconda.org/conda-forge/ndvi2gif)
[![Build status](https://github.com/Digdgeo/Ndvi2Gif/actions/workflows/python-publish.yml/badge.svg)](https://github.com/Digdgeo/Ndvi2Gif/actions/workflows/python-publish.yml)

![NDVI2GIF Köln](https://i.imgur.com/Y5dOWIk.jpeg)
*Richter's stained glass in Cologne Cathedral. Inspiration for this library.*

**Ndvi2Gif** is a Python library for multi-temporal remote sensing analysis with Google Earth Engine. It provides seasonal compositing, 40+ vegetation and environmental indices, SAR preprocessing, time series analytics, land cover classification, and hydroperiod analysis — all server-side, without downloading raw data.

Built on top of [Google Earth Engine](https://github.com/google/earthengine-api) and [geemap](https://github.com/giswqs/geemap). Adapted and extended through its use in the [eLTER](https://elter-ri.eu/) and [SUMHAL](https://lifewatcheric-sumhal.csic.es/) projects.

---

## 📚 Documentation

**https://digdgeo.github.io/Ndvi2Gif/**

- [Getting Started](https://digdgeo.github.io/Ndvi2Gif/getting_started/installation.html)
- [API Reference](https://digdgeo.github.io/Ndvi2Gif/reference/api.html)
- [Indices Catalog](https://digdgeo.github.io/Ndvi2Gif/reference/indices.html)
- [Example Notebooks](https://digdgeo.github.io/Ndvi2Gif/tutorials/)

---

## ✨ What's New in v1.3.0

**New module: `SpatialPhenologyAnalyzer`** — GEE-native, **per-pixel phenology rasters** (Start/Peak/End of Season and derived metrics) for the whole ROI, entirely server-side. Three methods (`threshold`, `derivative`, and `harmonic` Fourier regression — the Earth Engine-native replacement for double-logistic fitting), per-year or multi-year aggregate outputs, and export to local GeoTIFF or Google Drive with named bands. See [CHANGELOG](CHANGELOG.md) for details.

Previously, in v1.2.0: `HydroperiodAnalyzer` (GEE-native flood duration analysis), SCL-based cloud masking for Sentinel-2, and numpy 2.x support.

---

## Modules

| Module | Description |
|--------|-------------|
| `NdviSeasonality` | Core engine: seasonal compositing, 40+ indices, 7 sensors, flexible ROI input |
| `HydroperiodAnalyzer` | Wetland flood duration analysis (days/year) with multi-year anomaly detection |
| `TimeSeriesAnalyzer` | Trend detection (Mann-Kendall, Sen's slope), phenology metrics, dashboards |
| `SpatialPhenologyAnalyzer` | Per-pixel phenology rasters (SOS/POS/EOS) via threshold, derivative or harmonic methods |
| `S1ARDProcessor` | Sentinel-1 SAR preprocessing: terrain correction, speckle filtering |
| `LandCoverClassifier` | Supervised (RF, SVM, CART) and unsupervised (K-means, LDA) classification |

---

## Installation

```bash
pip install ndvi2gif
```

```bash
conda install -c conda-forge ndvi2gif
```

---

## Quick Start

```python
import ee
from ndvi2gif import NdviSeasonality

ee.Authenticate()
ee.Initialize(project='your-project-id')

# Monthly NDVI composites from Sentinel-2 (2018–2024)
ndvi = NdviSeasonality(
    roi=your_roi, sat='S2', periods=12,
    start_year=2018, end_year=2024,
    key='percentile', percentile=85, index='ndvi'
)

ndvi.get_gif(name='ndvi_evolution.gif')
```

Yes, it makes nice GIFs — but it's much more than that.

![GIF Example](https://i.imgur.com/xvrPYMH.gif)
*Crop pattern dance around Los Palacios y Villafranca (SW Spain)*

---

## What you can do with it

- **Compute pixel-wise statistics** over any region and time span — seasonal medians, percentiles, multi-year aggregations
- **Monitor 40+ indices** across Sentinel-1/2/3, Landsat (4–9), MODIS, ERA5-Land, and CHIRPS
- **Analyse wetland hydroperiod** and multi-year flood anomalies with `HydroperiodAnalyzer`
- **Detect trends and phenology** (SOS, EOS, POS, Length of Season) with `TimeSeriesAnalyzer`
- **Classify land cover** with multi-temporal feature stacks and Random Forest, SVM, or K-means
- **Preprocess Sentinel-1 SAR** with terrain correction and speckle filtering
- **Export** to GeoTIFF, Google Drive, or Earth Engine Assets
- **Use any ROI**: shapefile, GeoJSON, drawn geometry, eLTER DEIMS ID, Sentinel-2 tile, or Landsat path/row

---

## Supported Sensors

Sentinel-1 (SAR) · Sentinel-2 SR · Sentinel-3 OLCI · Landsat 4–9 SR · MODIS MOD09A1 · ERA5-Land · CHIRPS

---

## Contributing

Bug reports and feature requests: [GitHub Issues](https://github.com/Digdgeo/Ndvi2Gif/issues)

Pull requests are welcome. See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines — it includes step-by-step instructions for adding new indices and datasets.

---

## Citation

If you use Ndvi2Gif in your research, please cite the JOSS paper:

```bibtex
@article{GarciaDiaz2026,
  author    = {García Díaz, Diego},
  title     = {Ndvi2Gif: A Python Package for Multi-Seasonal Remote Sensing Analysis with Google Earth Engine},
  journal   = {Journal of Open Source Software},
  year      = {2026},
  volume    = {11},
  number    = {124},
  pages     = {10654},
  publisher = {The Open Journal},
  doi       = {10.21105/joss.10654},
  url       = {https://doi.org/10.21105/joss.10654}
}
```

## Acknowledgments

Special thanks to [Qiusheng Wu](https://github.com/giswqs) for his invaluable work in developing and promoting open-source geospatial software, to the Google Earth Engine team, and to the broader open-source geospatial community.

## License

MIT — see [LICENSE.txt](LICENSE.txt)
