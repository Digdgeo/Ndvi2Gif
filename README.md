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

## ✨ What's New in v1.5.0

**Raw reflectance bands** — the standardized bands are now selectable through `index=` on Sentinel-2, Landsat and MODIS (`'blue'`, `'green'`, `'red'`, `'nir'`, `'swir1'`, `'swir2'`, plus `'red_edge1-3'` on S2), returned as surface reflectance in 0–1 on every sensor. Spectral indices are ratios and cancel out changes in brightness, so a pixel can hold exactly the same NDVI while its reflectance drifts; radiometric work needs the bands themselves. Combined with the dispersion reducers they map how invariant each pixel is across a series — the basis for picking pseudo-invariant features.

**`key='count'`** — valid observations per pixel: the quality layer that says which parts of a dispersion map can be trusted.

Previously, in v1.4.0: **dispersion reducers** — four `key` options (`'std'`, `'variance'`, `'range'`, `'cv'`) that map how much an index **varies** inside each period instead of its typical level, useful for phenological change, disturbances and unstable surfaces such as flooded areas. They work with every sensor and index, and the composites keep the usual period band names, so they export, animate and analyse like any other.

**Downloadable water masks** — `HydroperiodAnalyzer.get_water_masks_stack()` flattens the per-date binary masks into a single `uint8` image (one band per acquisition date), and `include_masks=True` sends them to Drive or to an Earth Engine asset alongside the hydroperiod. Cloudy pixels and pixels no satellite ever saw get their own codes, so a time series can tell them apart.

And in v1.3.0: `SpatialPhenologyAnalyzer`, GEE-native per-pixel phenology rasters (SOS/POS/EOS) with threshold, derivative and harmonic methods. See [CHANGELOG](CHANGELOG.md) for details.

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

- **Compute pixel-wise statistics** over any region and time span — seasonal medians, percentiles, multi-year aggregations, or dispersion (std, variance, range, CV) to map variability instead of level
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
