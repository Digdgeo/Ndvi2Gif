"""Top-level package for ndvi2gif.

Provides the main entry points for seasonal compositing, SAR preprocessing,
time series / spatial trend analysis, and land-cover classification
with Google Earth Engine.
"""

__author__ = "Diego García Díaz"
__email__ = "diegogarcia@ebd.csic.es"
__version__ = "0.6.0"

from .ndvi2gif import NdviSeasonality, scale_OLI, scale_ETM
from .s1_ard import S1ARDProcessor
from .timeseries import TimeSeriesAnalyzer, SpatialTrendAnalyzer
from .clasification import LandCoverClassifier

__all__ = [
    "NdviSeasonality",
    "S1ARDProcessor",
    "TimeSeriesAnalyzer",
    "SpatialTrendAnalyzer",
    "LandCoverClassifier",
    "scale_OLI",
    "scale_ETM",
]