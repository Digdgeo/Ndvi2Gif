# API Reference

This page documents the public API of **ndvi2gif**. It is generated directly
from the source-code docstrings, so it always matches the installed version of
the package.

The library exposes seven main classes plus a couple of helper functions:

| Object | Module | Purpose |
| :----- | :----- | :------ |
| [`NdviSeasonality`](#api-ndviseasonality) | `ndvi2gif.ndvi2gif` | Temporal compositing, indices, ROI handling, export and GIFs |
| [`S1ARDProcessor`](#api-s1ardprocessor) | `ndvi2gif.s1_ard` | Sentinel-1 Analysis Ready Data preprocessing |
| [`TimeSeriesAnalyzer`](#api-timeseriesanalyzer) | `ndvi2gif.timeseries` | Point/region time series, trends and phenology |
| [`SpatialTrendAnalyzer`](#api-spatialtrendanalyzer) | `ndvi2gif.timeseries` | Per-pixel trend rasters (server-side) |
| [`SpatialPhenologyAnalyzer`](#api-spatialphenologyanalyzer) | `ndvi2gif.timeseries` | Per-pixel phenology rasters (server-side) |
| [`LandCoverClassifier`](#api-landcoverclassifier) | `ndvi2gif.clasification` | Supervised/unsupervised land cover classification |
| [`HydroperiodAnalyzer`](#api-hydroperiodanalyzer) | `ndvi2gif.hydroperiod` | Surface water / flooding metrics |

```{note}
All classes and helper functions can be imported directly from the top-level
package, e.g. `from ndvi2gif import NdviSeasonality, HydroperiodAnalyzer`.
```

(api-ndviseasonality)=
## NdviSeasonality

```{eval-rst}
.. autoclass:: ndvi2gif.NdviSeasonality
   :members:
   :show-inheritance:
```

(api-s1ardprocessor)=
## S1ARDProcessor

```{eval-rst}
.. autoclass:: ndvi2gif.S1ARDProcessor
   :members:
   :show-inheritance:
```

(api-timeseriesanalyzer)=
## TimeSeriesAnalyzer

```{eval-rst}
.. autoclass:: ndvi2gif.TimeSeriesAnalyzer
   :members:
   :show-inheritance:
```

(api-spatialtrendanalyzer)=
## SpatialTrendAnalyzer

```{eval-rst}
.. autoclass:: ndvi2gif.SpatialTrendAnalyzer
   :members:
   :show-inheritance:
```

(api-spatialphenologyanalyzer)=
## SpatialPhenologyAnalyzer

```{eval-rst}
.. autoclass:: ndvi2gif.SpatialPhenologyAnalyzer
   :members:
   :show-inheritance:
```

(api-landcoverclassifier)=
## LandCoverClassifier

```{eval-rst}
.. autoclass:: ndvi2gif.LandCoverClassifier
   :members:
   :show-inheritance:
```

(api-hydroperiodanalyzer)=
## HydroperiodAnalyzer

```{eval-rst}
.. autoclass:: ndvi2gif.HydroperiodAnalyzer
   :members:
   :show-inheritance:
```

## Utility functions

```{eval-rst}
.. autofunction:: ndvi2gif.scale_OLI

.. autofunction:: ndvi2gif.scale_ETM
```

## See Also

- [Indices Reference](indices.md) — complete list of spectral, SAR and climate variables
- [Datasets Reference](datasets.md) — details on the supported satellite/climate platforms
- [Tutorials](../tutorials/basic_ndvi.md) — step-by-step, worked examples
- [GitHub Repository](https://github.com/Digdgeo/Ndvi2Gif) — source code
