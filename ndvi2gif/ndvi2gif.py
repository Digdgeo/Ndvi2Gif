"""
ndvi2gif - Remote Sensing Time Series Analysis and Visualization
================================================================

A Python library for generating temporal composite animations and analyses from 
multiple satellite data sources, supporting various spectral and radar indices.

Key Features
------------
* Multi-sensor support (Sentinel-1/2/3, Landsat, MODIS)
* 40+ vegetation and environmental indices
* Flexible temporal compositing (seasonal, monthly, custom)
* Advanced SAR preprocessing with ARD pipeline
* Automated GIF generation for time series visualization
* Large-area processing with fishnet tiling

Quick Start
-----------
.. code-block:: python

    from ndvi2gif import NdviSeasonality
    
    # Create seasonal NDVI composites
    processor = NdviSeasonality(
        roi='path/to/roi.shp',
        sat='S2',
        index='ndvi',
        periods=4,
        start_year=2020,
        end_year=2023
    )
    
    # Generate animated visualization
    processor.get_gif('ndvi_trend.gif')

Installation
------------
.. code-block:: bash

    pip install ndvi2gif

Requirements
------------
* earthengine-api >= 0.1.320
* geemap >= 0.20.0
* geopandas >= 0.12.0
* numpy >= 1.21.0

Module Contents
---------------

Author: Diego García Díaz
Date: 2024
License: MIT
"""

import os
import ee
import math
import geemap
import requests
import zipfile
import geopandas as gpd
from geemap import zonal_statistics
from io import BytesIO
import calendar
from datetime import datetime, timedelta
from typing import Optional
from .s1_ard import S1ARDProcessor
#from .timeseries import TimeSeriesAnalyzer, SpatialTrendAnalyzer
#from .clasification import LandCoverClassifier

def scale_OLI(image):
    """
    Scale Landsat 8-9 OLI/TIRS sensor data to surface reflectance.
    
    Applies Collection 2 Level-2 scaling factors to convert digital numbers (DN)
    to surface reflectance values in the range [0, 1]. The scaling formula
    (DN * 0.0000275 - 0.2) is specific to Landsat Collection 2 products.
    
    Parameters
    ----------
    image : ee.Image
        Raw Landsat 8 or 9 image from Collection 2 Level-2 with SR bands.
        Required bands: SR_B2, SR_B3, SR_B4, SR_B5, SR_B6, SR_B7
        
    Returns
    -------
    ee.Image
        Image with scaled and renamed optical bands:
        
        * ``Blue`` : SR_B2 scaled (0.45-0.51 μm)
        * ``Green`` : SR_B3 scaled (0.53-0.59 μm)
        * ``Red`` : SR_B4 scaled (0.64-0.67 μm)
        * ``Nir`` : SR_B5 scaled (0.85-0.88 μm)
        * ``Swir1`` : SR_B6 scaled (1.57-1.65 μm)
        * ``Swir2`` : SR_B7 scaled (2.11-2.29 μm)
        
    Examples
    --------
    Apply to a single Landsat 8 image::
    
        l8_image = ee.Image('LANDSAT/LC08/C02/T1_L2/LC08_044034_20140318')
        scaled_image = scale_OLI(l8_image)
        print(scaled_image.bandNames().getInfo())
        # Output: ['Blue', 'Green', 'Red', 'Nir', 'Swir1', 'Swir2', ...]
    
    Apply to an image collection::
    
        l8_collection = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
        scaled_collection = l8_collection.map(scale_OLI)

    Raises
    ------
    KeyError
        If expected OLI surface reflectance bands are missing in `image`
        (e.g., ``SR_B2``–``SR_B7``).
    TypeError
        If `image` is not an ``ee.Image``.
    
    Notes
    -----
    The scaling coefficients are defined by USGS for Collection 2:
    
    * Scale factor: 0.0000275
    * Offset: -0.2
    * Valid range after scaling: typically -0.2 to 1.5
    * Negative values may occur in water or shadow areas
    
    See Also
    --------
    scale_ETM : Scaling function for Landsat 4-5-7
    
    References
    ----------
    USGS (2021). Landsat Collection 2 Level-2 Science Products.
    https://www.usgs.gov/landsat-missions/landsat-collection-2-level-2-science-products
    """
    opticalBands = image.select(['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7']).multiply(0.0000275).add(-0.2).rename(['Blue', 'Green', 'Red', 'Nir', 'Swir1', 'Swir2'])
    return image.addBands(opticalBands, None, True)
    

def _to_reflectance(image):
    """
    Rescale the standardized optical bands of Sentinel-2 or MODIS to 0-1.

    Both collections store surface reflectance as integers multiplied by
    10000. Only the standardized band names are touched, so auxiliary bands
    (MODIS LST, added later) keep their own units; ``addBands`` with
    ``overwrite=True`` keeps the image properties and footprint.
    """
    optical = ['Blue', 'Green', 'Red', 'Red_Edge1', 'Red_Edge2', 'Red_Edge3',
               'Nir', 'Swir1', 'Swir2']
    present = image.bandNames().filter(ee.Filter.inList('item', optical))
    return image.addBands(image.select(present).multiply(0.0001), None, True)


# Band-averaged extraterrestrial solar irradiance of the OLCI bands, in
# W m-2 um-1, keyed by the names ndvi2gif gives the bands. Computed by
# convolving the Thuillier et al. (2003) solar spectrum with ESA's mean OLCI
# spectral response functions (S3A_OL_SRF_20160713_mean_rsr.nc4 and its S3B
# equivalent), which is how a band-averaged F0 is defined.
#
# These numbers matter because Earth Engine serves COPERNICUS/S3/OLCI as
# top-of-atmosphere *radiance*. F0 spans a factor of 2.8 across the sixteen
# bands loaded here (1937.6 at Oa04 against 699.7 at Oa21), and unlike the
# illumination terms it does NOT cancel in a band ratio, so every index
# computed straight from the radiance was wrong -- the same class of bug as
# the Sentinel-2 and MODIS scale factors fixed in 1.6.0.
#
# S3A and S3B differ by 0.10% on average and 1.39% at worst (Oa01 only).
_OLCI_SOLAR_IRRADIANCE = {
    'S3A': {
        'Violet': 1515.868,   # Oa01
        'Blue': 1708.009,   # Oa02
        'Blue2': 1890.814,   # Oa03
        'Blue_Green': 1937.628,   # Oa04
        'Green': 1918.784,   # Oa05
        'Green2': 1796.865,   # Oa06
        'Red': 1649.274,   # Oa07
        'Red2': 1530.055,   # Oa08
        'Red3': 1494.781,   # Oa09
        'Red_Edge1': 1468.962,   # Oa10
        'Red_Edge2': 1402.691,   # Oa11
        'Nir': 1266.557,   # Oa12
        'Nir2': 1173.373,   # Oa16
        'Nir3': 959.221,   # Oa17
        'Nir4': 930.995,   # Oa18
        'Nir5': 699.731,   # Oa21
    },
    'S3B': {
        'Violet': 1536.937,   # Oa01
        'Blue': 1709.147,   # Oa02
        'Blue2': 1891.177,   # Oa03
        'Blue_Green': 1936.118,   # Oa04
        'Green': 1920.022,   # Oa05
        'Green2': 1797.488,   # Oa06
        'Red': 1649.397,   # Oa07
        'Red2': 1530.440,   # Oa08
        'Red3': 1495.450,   # Oa09
        'Red_Edge1': 1469.715,   # Oa10
        'Red_Edge2': 1403.510,   # Oa11
        'Nir': 1266.454,   # Oa12
        'Nir2': 1174.141,   # Oa16
        'Nir3': 959.101,   # Oa17
        'Nir4': 931.341,   # Oa18
        'Nir5': 699.833,   # Oa21
    },
}

_OLCI_BANDS = list(_OLCI_SOLAR_IRRADIANCE['S3A'])


def _solar_geometry(image):
    """
    Cosine of the solar zenith angle and the Earth-Sun distance factor.

    Earth Engine's OLCI collection carries no solar angles: the SAFE product
    holds them in tie-point grids that are not ingested, so they have to be
    computed. This is the NOAA solar position algorithm, good to about a tenth
    of a degree, which is far below the other error terms here.

    Returns a two-band image, ``cos_sza`` per pixel (the angle changes across
    a 1270 km swath, so a scene-wide constant would not do) and ``e0``, the
    inverse squared Earth-Sun distance in AU, constant for the scene.
    """
    date = ee.Date(image.get('system:time_start'))
    doy = ee.Number(date.getRelative('day', 'year')).add(1)
    hour = ee.Number(date.getFraction('day')).multiply(24)

    # fractional year, radians
    gamma = doy.subtract(1).add(hour.subtract(12).divide(24)) \
               .multiply(2 * math.pi / 365)

    def harmonic(*terms):
        """Sum of a0 + sum(an*cos(n*gamma) + bn*sin(n*gamma)), as an ee.Number."""
        total = ee.Number(terms[0])
        for n, (a, b) in enumerate(zip(terms[1::2], terms[2::2]), start=1):
            angle = gamma.multiply(n)
            total = total.add(angle.cos().multiply(a)).add(angle.sin().multiply(b))
        return total

    # equation of time, minutes
    eqtime = harmonic(0.000075, 0.001868, -0.032077,
                      -0.014615, -0.040849).multiply(229.18)
    # solar declination, radians
    decl = harmonic(0.006918, -0.399912, 0.070257,
                    -0.006758, 0.000907, -0.002697, 0.00148)
    # eccentricity correction, 1/d^2 with d in AU
    e0 = harmonic(1.00011, 0.034221, 0.00128, 0.000719, 0.000077)

    coords = ee.Image.pixelLonLat()
    lat = coords.select('latitude').multiply(math.pi / 180)
    # true solar time, minutes, then the hour angle in radians
    true_solar = coords.select('longitude').multiply(4) \
                       .add(eqtime).add(hour.multiply(60))
    hour_angle = true_solar.divide(4).subtract(180).multiply(math.pi / 180)

    cos_sza = lat.sin().multiply(ee.Number(decl.sin())) \
                 .add(lat.cos().multiply(ee.Number(decl.cos())).multiply(hour_angle.cos()))

    return cos_sza.rename('cos_sza').addBands(ee.Image.constant(e0).rename('e0'))


def _s3_to_reflectance(image):
    """
    Convert Sentinel-3 OLCI top-of-atmosphere radiance to TOA reflectance.

    ``rho = pi * L * d^2 / (F0 * cos(theta_s))``, with F0 taken per band and
    per spacecraft from :data:`_OLCI_SOLAR_IRRADIANCE` and the geometry from
    :func:`_solar_geometry`.

    Note what each term does to an index. ``cos(theta_s)`` and ``d^2`` are the
    same for every band, so they cancel in any band ratio and only matter for
    absolute values and for indices that are not ratios. ``F0`` does not
    cancel, which is why the radiance this replaces gave wrong answers even
    for NDVI.

    This is TOA reflectance, not surface reflectance: no atmospheric
    correction is applied, and over water in particular the signal is
    dominated by Rayleigh scattering.
    """
    spacecraft = ee.String(image.get('spacecraft'))
    irradiance = ee.Dictionary(ee.Algorithms.If(
        spacecraft.compareTo('S3B').eq(0),
        _OLCI_SOLAR_IRRADIANCE['S3B'],
        _OLCI_SOLAR_IRRADIANCE['S3A']))

    present = image.bandNames().filter(ee.Filter.inList('item', _OLCI_BANDS))
    f0 = ee.Image.constant(irradiance.values(present)).rename(present)

    geometry = _solar_geometry(image)
    illumination = f0.multiply(geometry.select('cos_sza')) \
                     .multiply(geometry.select('e0'))

    reflectance = image.select(present).multiply(math.pi).divide(illumination)
    return image.addBands(reflectance, None, True)


def scale_ETM(image):
    """
    Scale Landsat 4-5-7 ETM+/TM sensor data to surface reflectance.
    
    Applies Collection 2 Level-2 scaling factors for Enhanced Thematic Mapper Plus
    (ETM+) on Landsat 7 and Thematic Mapper (TM) on Landsat 4-5. Uses the same
    scaling formula as OLI but with different band numbering scheme.
    
    Parameters
    ----------
    image : ee.Image
        Raw Landsat 4, 5, or 7 image from Collection 2 Level-2.
        Required bands: SR_B1, SR_B2, SR_B3, SR_B4, SR_B5, SR_B7
        Note: SR_B6 (thermal) is excluded as it requires different scaling
        
    Returns
    -------
    ee.Image
        Image with scaled and renamed optical bands:
        
        * ``Blue`` : SR_B1 scaled (0.45-0.52 μm)
        * ``Green`` : SR_B2 scaled (0.52-0.60 μm)
        * ``Red`` : SR_B3 scaled (0.63-0.69 μm)
        * ``Nir`` : SR_B4 scaled (0.77-0.90 μm)
        * ``Swir1`` : SR_B5 scaled (1.55-1.75 μm)
        * ``Swir2`` : SR_B7 scaled (2.09-2.35 μm)
        
    Examples
    --------
    Apply to a Landsat 5 image::
    
        l5_image = ee.Image('LANDSAT/LT05/C02/T1_L2/LT05_044034_20110716')
        scaled_image = scale_ETM(l5_image)
    
    Apply to mixed Landsat collection::
    
        landsat_457 = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2')
        scaled_collection = landsat_457.map(scale_ETM)

    Raises
    ------
    KeyError
        If expected ETM+ surface reflectance bands are missing in `image`
        (e.g., ``SR_B1``–``SR_B5`` and ``SR_B7``).
    TypeError
        If `image` is not an ``ee.Image``.
    
    Notes
    -----
    Band numbering differs between ETM/TM and OLI sensors:
    
    * ETM/TM Band 1 (Blue) → OLI Band 2
    * ETM/TM Band 2 (Green) → OLI Band 3
    * ETM/TM Band 3 (Red) → OLI Band 4
    * No Band 6 processing (thermal requires different scaling)
    
    Warnings
    --------
    Landsat 7 ETM+ has scan line corrector failure (SLC-off) after May 2003,
    resulting in data gaps. Consider using gap-filling techniques or
    focusing on Landsat 4-5 for historical analysis.
    
    See Also
    --------
    scale_OLI : Scaling function for Landsat 8-9
    
    References
    ----------
    USGS (2021). Landsat Collection 2 Level-2 Science Products.
    https://www.usgs.gov/landsat-missions/landsat-collection-2-level-2-science-products
    """
    opticalBands = image.select(['SR_B1','SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7']).multiply(0.0000275).add(-0.2).rename(['Blue', 'Green', 'Red', 'Nir', 'Swir1', 'Swir2'])
    return image.addBands(opticalBands, None, True)

class NdviSeasonality:
    """
    Generate remote sensing index seasonal composition GIFs and images.
    
    Comprehensive framework for creating temporal composites from multiple
    satellite data sources. Supports 40+ spectral and radar indices with
    configurable temporal periods and statistical reducers.
    
    Parameters
    ----------
    roi : ee.Geometry, str, list, or None, optional
        Region of interest specification:
        
        * ``None`` : Default region (Andalusia, Spain)
        * ``ee.Geometry`` : Direct Earth Engine geometry
        * ``str`` : Path or special format
            - ``'*.shp'`` : Shapefile path
            - ``'*.geojson'`` : GeoJSON file path
            - ``'deimsid/XXXXX'`` : DEIMS.org site ID
            - ``'wrs:path,row'`` : Landsat WRS-2 tile (e.g., 'wrs:200,32')
            - ``'s2:XXXXX'`` : Sentinel-2 MGRS tile (e.g., 's2:30TXN')
        * ``list`` : Feature collection from Map.draw_features
        
    periods : int, optional
        Number of temporal periods per year:
        
        * ``4`` : Seasonal (winter, spring, summer, autumn)
        * ``12`` : Monthly (january through december)
        * ``24`` : Bi-monthly (~15 days each, p1 through p24)
        * Other : Custom equal division (p1 through pN)
        
        Default is 4.
        
    start_year : int, optional
        Starting year for analysis (inclusive). Must be within satellite
        data availability. Default is 2016.
        
    end_year : int, optional
        Ending year for analysis (inclusive). Analysis includes years
        from start_year to end_year. Default is 2020.
        
    sat : {'S2', 'S1', 'Landsat', 'MODIS', 'S3'}, optional
        Satellite sensor selection:
        
        * ``'S2'`` : Sentinel-2 MSI (optical, 10-20m, 2015-present)
        * ``'S1'`` : Sentinel-1 SAR (radar, 10m, 2014-present)
        * ``'Landsat'`` : Merged L4-5-7-8-9 (optical, 30m, 1982-present)
        * ``'MODIS'`` : Terra/Aqua (optical, 500m, 2000-present)
        * ``'S3'`` : Sentinel-3 OLCI (ocean/land, 300m, 2016-present)
        
        Default is 'S2'.
        
    key : {'max', 'min', 'median', 'mean', 'sum', 'percentile', 'std', 'variance', 'range', 'cv', 'count'}, optional
        Statistical reducer for temporal aggregation:

        * ``'max'`` : Maximum value (vegetation peak detection)
        * ``'min'`` : Minimum value
        * ``'median'`` : Median value (robust to outliers)
        * ``'mean'`` : Mean value (smooth temporal profiles)
        * ``'sum'`` : Total sum (ideal for precipitation, accumulation)
        * ``'percentile'`` : Custom percentile (set with percentile param)
        * ``'std'`` : Standard deviation (within-period variability)
        * ``'variance'`` : Variance (squared variability)
        * ``'range'`` : Maximum minus minimum (amplitude)
        * ``'cv'`` : Coefficient of variation, std divided by mean.
          Only meaningful for indices that stay positive
        * ``'count'`` : Number of valid observations per pixel (data
          availability, not a value statistic)

        Default is 'max'.
        
    index : str, optional
        Spectral or radar index to compute. Available indices depend on
        satellite. See :meth:`get_available_indices` for full list.
        Common indices:
        
        * Vegetation: ``'ndvi'``, ``'evi'``, ``'savi'``, ``'gndvi'``
        * Water: ``'ndwi'``, ``'mndwi'``, ``'awei'``
        * Burn/Fire: ``'nbr'``, ``'nbri'``
        * SAR: ``'vh'``, ``'vv'``, ``'rvi'``, ``'vv_vh_ratio'``
        
        Default is 'ndvi'.
        
    percentile : int, optional
        Percentile value when key='percentile'. Range: 0-100.
        Common values:
        
        * 10-25: Lower percentiles (minimum-like)
        * 50: Median equivalent
        * 75-95: Upper percentiles (maximum-like)
        
        Default is 90.
        
    orbit : {'BOTH', 'ASCENDING', 'DESCENDING'}, optional
        Sentinel-1 orbit direction (only for sat='S1'):
        
        * ``'BOTH'`` : All orbits (maximum temporal coverage)
        * ``'ASCENDING'`` : Single orbit (geometric consistency)
        * ``'DESCENDING'`` : Single orbit (geometric consistency)
        
        Default is 'BOTH'.
        
    normalize_sar : bool, optional
        If True, applies Z-score normalization to SAR indices for
        better comparability with optical indices. Only applies when
        sat='S1'. Default is False.
        
    use_sar_ard : bool, optional
        If True, applies advanced SAR preprocessing (ARD pipeline) including
        terrain correction and sophisticated speckle filtering. Recommended
        for mountainous areas. Default is True.
        
    sar_speckle_filter : str or None, optional
        Speckle filter algorithm for SAR (when use_sar_ard=True):
        
        * ``'REFINED_LEE'`` : Refined Lee with edge preservation
        * ``'LEE'`` : Standard Lee filter
        * ``'GAMMA_MAP'`` : Gamma Maximum A Posteriori
        * ``'LEE_SIGMA'`` : Lee Sigma filter
        * ``'BOXCAR'`` : Simple mean filter
        * ``None`` : No speckle filtering
        
        Default is 'REFINED_LEE'.
        
    sar_terrain_correction : bool, optional
        Enable radiometric terrain correction for SAR data.
        Essential for mountainous regions. Default is True.
        
    sar_terrain_model : {'VOLUME', 'SURFACE'}, optional
        Scattering model for terrain correction:
        
        * ``'VOLUME'`` : Volume scattering (vegetation, crops)
        * ``'SURFACE'`` : Surface scattering (bare soil, water)
        
        Default is 'VOLUME'.
    
    Attributes
    ----------
    roi : ee.Geometry
        Processed region of interest geometry
    periods : int
        Number of temporal periods per year
    start_year : int
        First year of analysis
    end_year : int
        Last year of analysis (inclusive)
    sat : str
        Selected satellite sensor
    index : str
        Selected spectral/radar index
    key : str
        Statistical reducer method
    percentile : int
        Percentile value for percentile reducer
    period_dates : list of list
        Date ranges for each period [[start, end], ...]
    period_names : list of str
        Names for each temporal period
    ndvi_col : ee.ImageCollection
        Configured satellite image collection
    imagelist : list of ee.Image
        Processed composite images
    
    Examples
    --------
    Basic seasonal NDVI analysis::
    
        >>> processor = NdviSeasonality(
        ...     roi='study_area.shp',
        ...     periods=4,
        ...     start_year=2020,
        ...     end_year=2023,
        ...     sat='S2',
        ...     index='ndvi'
        ... )
        >>> processor.get_gif('seasonal_ndvi.gif')
    
    Monthly SAR analysis with percentile::
    
        >>> sar_processor = NdviSeasonality(
        ...     periods=12,
        ...     sat='S1',
        ...     index='vh',
        ...     key='percentile',
        ...     percentile=90,
        ...     orbit='DESCENDING'
        ... )
        >>> collection = sar_processor.get_year_composite()
    
    Using DEIMS site with custom periods::
    
        >>> deims_processor = NdviSeasonality(
        ...     roi='deimsid/11696159-444f-4e06-b537-d4c5c0a4e97d',
        ...     periods=8,  # 8 periods per year
        ...     sat='MODIS',
        ...     index='evi'
        ... )
    
    Raises
    ------
    ValueError
        If satellite not supported, index not available for satellite,
        orbit parameter invalid, or ROI cannot be processed.
    ImportError
        If required dependencies missing (e.g., deims package).
    
    See Also
    --------
    get_available_indices : Query available indices for satellites
    get_year_composite : Generate temporal composite images
    get_gif : Create animated visualization
    get_export : Export composites to GeoTIFF files
    
    Notes
    -----
    The class automatically validates index-satellite compatibility during
    initialization. Default ROI covers part of Andalusia, Spain for testing.
    Temporal periods are generated dynamically to support flexible analysis.
    
    References
    ----------
    .. [1] Gorelick et al. (2017). Google Earth Engine: Planetary-scale 
           geospatial analysis for everyone. Remote Sensing of Environment.
    """
    
    def __init__(self, roi=None, periods=4, start_year=2016, end_year=2020,
                 sat='S2', key='max', index='ndvi', percentile=90, orbit='BOTH', normalize_sar=False,
                 use_sar_ard=True, sar_speckle_filter='REFINED_LEE', sar_terrain_correction=True,
                sar_terrain_model='VOLUME', cloud_filter=True, max_cloud_cover=20, scl_mask=True):
        """
        Initialize NdviSeasonality object for temporal remote sensing analysis.
        
        Sets up the complete analysis framework including ROI processing,
        satellite collection configuration, index validation, and temporal
        period generation. Provides extensive input validation and flexible
        ROI specification options.
        
        Parameters
        ----------
        roi : ee.Geometry, str, or None, optional
            Region of interest specification. Multiple formats supported:
            
            - None: Uses default Andalusia region, Spain
            - ee.Geometry: Direct Earth Engine geometry object
            - str ending with '.shp': Path to shapefile
            - str ending with '.geojson': Path to GeoJSON file
            - str starting with 'deimsid': DEIMS site ID (format: 'deimsid/XXXXX')
            - str starting with 'wrs:': Landsat WRS path/row (format: 'wrs:XXX,YY')
            - str starting with 's2:': Sentinel-2 MGRS tile (format: 's2:XXXXX')
            - list: Feature collection (e.g., from Map.draw_features)
            
            Default is None.
            
        periods : int, optional
            Number of temporal periods to divide each year. Common options:
            
            - 4: Traditional seasons (winter, spring, summer, autumn)
            - 12: Monthly periods (january, february, ...)
            - 24: Bi-monthly periods (~15 days each, named p1, p2, ...)
            - Other: Custom equal division of year (named p1, p2, ...)
            
            Default is 4.
            
        start_year : int, optional
            Starting year for temporal analysis (inclusive).
            Must be within satellite data availability range.
            Default is 2016.
            
        end_year : int, optional
            Ending year for temporal analysis (inclusive).
            Analysis includes years from start_year to end_year.
            Default is 2020.
            
        sat : str, optional
            Satellite sensor selection. Available options:
            
            - 'S2': Sentinel-2 MSI (optical, 10-20m resolution, 2015-present)
            - 'Landsat': Landsat 4-5-7-8-9 merged (optical, 30m, 1982-present)
            - 'MODIS': MODIS Terra/Aqua (optical, 500m, 2000-present)
            - 'S1': Sentinel-1 SAR (radar, 10m resolution, 2014-present)
            - 'S3': Sentinel-3 OLCI (ocean/land color, 300m, 2016-present)
            
            Default is 'S2'.
            
        key : str, optional
            Statistical reducer for temporal aggregation within each period:

            - 'max': Maximum value (default, good for vegetation peak detection)
            - 'min': Minimum value
            - 'median': Median value (robust to outliers and clouds)
            - 'mean': Mean value (smooth temporal profiles)
            - 'sum': Total sum (accumulation variables such as precipitation)
            - 'percentile': Custom percentile (specify with percentile parameter)
            - 'std': Standard deviation of the observations in the period
            - 'variance': Variance of the observations in the period
            - 'range': Maximum minus minimum (amplitude within the period)
            - 'cv': Coefficient of variation (std / mean). Only meaningful for
              indices that stay positive; it is unstable when the mean is close
              to zero and negative for Sentinel-1 backscatter in dB
            - 'count': Number of valid (non-masked) observations per pixel in
              the period. Not a value statistic but a quality layer: dispersion
              computed from three observations means very little, so 'count'
              tells which parts of the result can be trusted

            The dispersion statistics ('std', 'variance', 'range' and 'cv') do
            not describe the typical level of the index but how much it varies
            inside each period, which is useful to detect phenological change,
            disturbances or unstable surfaces such as flooded areas.

            Default is 'max'.
            
        index : str, optional
            Spectral or radar index to compute. Available indices depend on satellite:
            
            **All optical sensors (S2, Landsat, MODIS):**
            ndvi, evi, ndwi, mndwi, savi, gndvi, avi, nbri, ndsi, aweinsh, awei,
            ndmi, msi, nmi, ndti, cri1, cri2, lai, pri, wdrvi
            
            **Sentinel-2 exclusive (Red Edge bands):**
            ireci, mcari, ndre, reip, psri, cire, mtci, s2rep, ndci
            
            **Sentinel-1 SAR indices:**
            vh, vv, rvi, vv_vh_ratio, dpsvi, rfdi, vsdi
            
            **Sentinel-3 OLCI exclusive:**
            oci, tsi, cdom, turbidity, spm, kd490, floating_algae,
            red_edge_position, fluorescence_height, water_leaving_reflectance
            
            Default is 'ndvi'.
            
        percentile : int, optional
            Percentile value when key='percentile'. Must be between 0-100.
            Common values:
            
            - 10, 25: Lower percentiles (minimum-like behavior)
            - 50: Median equivalent
            - 75, 90, 95: Upper percentiles (maximum-like behavior)
            
            Default is 90.
            
        orbit : str, optional
            Sentinel-1 orbit direction (only used when sat='S1'):
            
            - 'BOTH': Use all available orbits (maximum temporal coverage)
            - 'ASCENDING': Use only ascending orbits (better geometric consistency)
            - 'DESCENDING': Use only descending orbits (better geometric consistency)
            
            For optical satellites, this parameter is ignored with a warning.
            Default is 'BOTH'.

        normalize_sar : bool, optional
            If True, normalizes all SAR indices to [0,1] range for better comparability
            with optical indices. Only applies when sat='S1'. Default is False.

        use_sar_ard : bool, optional
        If True, applies advanced SAR preprocessing using ARD pipeline.
        Includes terrain correction and sophisticated speckle filtering.
        Recommended for mountainous areas or when high quality is needed.
        Default is True.
        
        sar_speckle_filter : str, optional
            Speckle filter algorithm for SAR preprocessing. Options:
            
            - 'REFINED_LEE': Refined Lee with edge preservation (recommended)
            - 'LEE': Standard Lee filter
            - 'GAMMA_MAP': Gamma Maximum A Posteriori
            - 'LEE_SIGMA': Lee Sigma filter
            - 'BOXCAR': Simple mean filter
            - None: No speckle filtering
            
            Only used when sat='S1' and use_sar_ard=True.
            Default is 'REFINED_LEE'.
            
        sar_terrain_correction : bool, optional
            Enable radiometric terrain correction for SAR data.
            Essential for mountainous regions to reduce topographic effects.
            Only used when sat='S1' and use_sar_ard=True.
            Default is True.
            
        sar_terrain_model : str, optional
            Scattering model for terrain correction:
            
            - 'VOLUME': Volume scattering (vegetation, crops, grasslands)
            - 'SURFACE': Surface scattering (bare soil, water, rock)
            
            Only used when sat='S1', use_sar_ard=True, and sar_terrain_correction=True.
            Default is 'VOLUME'.

        cloud_filter : bool, optional
            If True, applies cloud filtering to optical sensors (S2, Landsat).
            On Sentinel-2 this is two steps: whole scenes above
            ``max_cloud_cover`` are dropped, then cloudy pixels are masked
            (SCL or QA60, see ``scl_mask``). ``False`` disables **both**, so
            clouds enter the composites. To keep the pixel mask without
            discarding scenes — useful in persistently cloudy regions, where
            few scenes pass a scene-level threshold — leave it True and set
            ``max_cloud_cover=100``. No effect on SAR or other sensors.
            Default is True.
        
        max_cloud_cover : int, optional
            Maximum cloud cover percentage (0-100) for initial filtering.
            Only used when cloud_filter=True. Default is 20.
            Lower values = stricter filtering but less data.
            
        Raises
        ------
        ValueError
            - If satellite is not supported
            - If index is not available for the selected satellite
            - If orbit parameter is invalid
            - If ROI cannot be processed (e.g., invalid WRS coordinates)
            
        ImportError
            If required dependencies are missing (e.g., deims package for DEIMS IDs)
            
        Examples
        --------
        >>> # Basic seasonal NDVI analysis with Sentinel-2
        >>> processor = NdviSeasonality(
        ...     roi=my_geometry,
        ...     periods=4,
        ...     start_year=2020,
        ...     end_year=2023,
        ...     sat='S2',
        ...     index='ndvi'
        ... )
        
        >>> # Monthly SAR analysis with 90th percentile
        >>> sar_processor = NdviSeasonality(
        ...     periods=12,
        ...     sat='S1',
        ...     index='vh',
        ...     key='percentile',
        ...     percentile=90,
        ...     orbit='DESCENDING'
        ... )
        
        >>> # Custom ROI from shapefile with Landsat EVI
        >>> landsat_processor = NdviSeasonality(
        ...     roi='/path/to/study_area.shp',
        ...     sat='Landsat',
        ...     index='evi',
        ...     key='median'
        ... )
        
        >>> # Using DEIMS site ID
        >>> deims_processor = NdviSeasonality(
        ...     roi='deimsid/https://deims.org/11696159-444f-4e06-b537-d4c5c0a4e97d',
        ...     sat='MODIS',
        ...     index='ndvi'
        ... )
        
        >>> # Using Landsat WRS coordinates
        >>> wrs_processor = NdviSeasonality(
        ...     roi='wrs:200,32',
        ...     sat='Landsat',
        ...     index='nbri'
        ... )
        
        Notes
        -----
        - The class automatically validates index-satellite compatibility
        - ROI processing handles multiple input formats with extensive error handling
        - Temporal periods are generated dynamically to support flexible time divisions
        - Satellite collections are configured with appropriate scaling and filtering
        - Default ROI covers part of Andalusia, Spain for testing purposes
        
        See Also
        --------
        get_year_composite : Generate temporal composite images
        get_export : Export all composite images
        get_gif : Create animated GIF from composites
        get_available_indices : Query available indices for a satellite
        """
        print('There we go again...')
        
        # Initialize ROI with comprehensive format support
        self.roi = roi
        if self.roi is None:
            # Default ROI: Andalusia region, Spain
            self.roi = ee.Geometry.Polygon(
                [[[-6.766047, 36.776586], 
                  [-6.766047, 37.202186], 
                  [-5.867729, 37.202186], 
                  [-5.867729, 36.776586], 
                  [-6.766047, 36.776586]]], None, False)
        elif isinstance(self.roi, str):
            # Handle string-based ROI specifications
            if self.roi.endswith('.shp'):
                self.roi = geemap.shp_to_ee(self.roi).geometry()
            elif self.roi.endswith('.geojson'):
                self.roi = geemap.geojson_to_ee(self.roi).geometry()
            elif self.roi.startswith('deimsid'):
                print('Con Deims hemos topado, amigo Sancho...')
                try:
                    import deims
                except ImportError:
                    raise ImportError("To use a DEIMS ID, you must install the `deims` package via pip:\n\n    pip install deims\n")
                id_ = self.roi.split('/')[-1]
                gdf = deims.getSiteBoundaries(id_)
                self.roi = geemap.geopandas_to_ee(gdf).geometry()
            elif self.roi.startswith('wrs:'):
                print('Loading Landsat WRS-2 geometry from GitHub...')
                path, row = map(int, self.roi.replace('wrs:', '').split(','))
                url = 'https://raw.githubusercontent.com/Digdgeo/Ndvi2Gif/master/data/l2tiles.geojson'
                wrs = gpd.read_file(url)
                subset = wrs[(wrs['PATH'] == path) & (wrs['ROW'] == row)]
                if subset.empty:
                    raise ValueError(f"No geometry found for Path {path}, Row {row}")
                print(f'Found Landsat tile for Path {path}, Row {row}')
                self.roi = geemap.geopandas_to_ee(subset).geometry()
            elif self.roi.startswith('s2:'):
                print('Loading Sentinel-2 MGRS tile from GitHub...')
                tile_id = self.roi.replace('s2:', '').strip().upper()
                url = 'https://raw.githubusercontent.com/Digdgeo/Ndvi2Gif/master/data/s2tiles_2d.geojson'
                s2 = gpd.read_file(url)
                subset = s2[s2['Name'] == tile_id]
                if subset.empty:
                    raise ValueError(f"No geometry found for Sentinel-2 tile {tile_id}")
                print(f'Found Sentinel-2 tile for {tile_id}')
                self.roi = geemap.geopandas_to_ee(subset).geometry()
            else:
                print('Invalid ROI path format')
        else:
            # Handle geometry objects and feature collections
            if isinstance(self.roi, list) and len(self.roi) > 0:
                # Handle lists of Features (like Map.draw_features or draw_last_feature)
                first_feature = self.roi[0]
                if hasattr(first_feature, 'geometry'):
                    self.roi = first_feature.geometry()
                else:
                    self.roi = ee.Geometry(first_feature)
            elif hasattr(self.roi, 'geometry'):
                self.roi = self.roi.geometry()
            elif not isinstance(self.roi, ee.Geometry):
                try:
                    self.roi = ee.Geometry(self.roi)
                except Exception as e:
                    print('Could not convert the provided roi to ee.Geometry')
                    print(f'ROI type: {type(self.roi)}')
                    print(f'Error: {e}')
                    # Use default ROI instead of failing
                    self.roi = ee.Geometry.Polygon(
                        [[[-6.766047, 36.776586], 
                        [-6.766047, 37.202186], 
                        [-5.867729, 37.202186], 
                        [-5.867729, 36.776586], 
                        [-6.766047, 36.776586]]], None, False)
                    print("Using default ROI")
        
        # Set temporal and processing parameters
        self.periods = periods
        self.start_year = start_year
        self.end_year = end_year
        valid_keys = ['max', 'min', 'median', 'percentile', 'mean', 'sum',
                      'std', 'variance', 'range', 'cv', 'count']
        if key not in valid_keys:
            raise ValueError(
                f"Statistic '{key}' is not supported. Available statistics (key) are: {valid_keys}"
            )
        if key == 'cv' and sat == 'S1' and not normalize_sar:
            print("Warning: key='cv' divides by the mean, and Sentinel-1 backscatter "
                  "is negative in dB, so the coefficient of variation will be negative "
                  "and hard to interpret. Use normalize_sar=True, or key='std'/'range'.")
        self.key = key
        self.percentile = percentile
        self.imagelist = []
        self.index = index
        # Store SAR preprocessing parameters
        self.use_sar_ard = use_sar_ard
        self.sar_speckle_filter = sar_speckle_filter
        self.sar_terrain_correction = sar_terrain_correction
        self.sar_terrain_model = sar_terrain_model

        # Store cloud filtering parameters
        self.cloud_filter = cloud_filter
        self.max_cloud_cover = max_cloud_cover
        self.scl_mask = scl_mask
        
        # Validate orbit parameter
        valid_orbits = ['BOTH', 'ASCENDING', 'DESCENDING']
        if orbit not in valid_orbits:
            raise ValueError(f"Orbit '{orbit}' is not valid. Available options are: {valid_orbits}")
        
        self.orbit = orbit
        self.normalize_sar = normalize_sar
    
        # Show warning if orbit is used with non-SAR sensors
        if sat != 'S1' and orbit != 'BOTH':
            print(f"Warning: orbit parameter '{orbit}' is only used with Sentinel-1. Ignoring for {sat}.")
            
        # Define basic optical indices (available on all optical sensors).
        # Only visible and NIR bands: Sentinel-3 OLCI has no SWIR or thermal
        self.optical_indices = {
            'ndvi', 'ndwi', 'evi', 'savi', 'gndvi', 'avi',
            'cri1', 'cri2', 'lai', 'pri', 'wdrvi', 'vci', 'cig'
        }

        # Sentinel-2 exclusive indices (Red Edge bands)
        self.s2_exclusive_indices = {
            'ireci', 'mcari', 'ndre', 'reip', 'psri', 'cire', 'mtci', 's2rep', 'ndci'
        }

        # Sentinel-3 OLCI exclusive indices
        self.s3_exclusive_indices = {
            'oci', 'tsi', 'cdom', 'turbidity', 'spm', 'kd490', 'floating_algae',
            'red_edge_position', 'fluorescence_height', 'water_leaving_reflectance'
        }

        # Indices requiring the SWIR bands (Swir1 ~1600 nm, Swir2 ~2200 nm),
        # which Sentinel-3 OLCI does not have
        self.swir_optical_indices = {
            'mndwi', 'nbri', 'ndsi', 'aweinsh', 'awei', 'ndmi', 'msi', 'nmi',
            'ndti', 'nbr', 'wi2015', 'ndbi', 'fai'
        }

        # Thermal indices. They read the surface temperature bands of Landsat
        # (ST_B10 / ST_B6) and MODIS (LST_Day_1km), and would return a fully
        # masked image on any other sensor
        self.thermal_indices = {'lst', 'utfvi'}

        # Raw reflectance bands exposed as selectable "indices". Spectral
        # indices are ratios and cancel out multiplicative brightness changes,
        # so NDVI can stay constant while reflectance drifts. Radiometric work
        # (selecting pseudo-invariant features, for instance) needs the bands
        # themselves. Sentinel-3 is left out on purpose: its bands are TOA
        # radiances with a different band set, so they are not reflectance
        self.raw_band_indices = {
            'blue', 'green', 'red', 'nir', 'swir1', 'swir2'
        }

        # Red edge bands, only present in Sentinel-2
        self.s2_rededge_bands = {'red_edge1', 'red_edge2', 'red_edge3'}

        # SAR indices (Sentinel-1)
        self.s1_indices = {
            'rvi', 'vv', 'vh', 'vv_vh_ratio', 'dpsvi', 'rfdi', 'vsdi'
        }

        # ERA5-Land climate reanalysis variables
        self.era5_variables = {
            # Temperature (Kelvin)
            'temperature_2m', 'dewpoint_temperature_2m', 'skin_temperature',
            'soil_temperature_level_1',
            # Temperature min/max (Kelvin)
            'temperature_2m_min', 'temperature_2m_max',
            'dewpoint_temperature_2m_min', 'dewpoint_temperature_2m_max',
            'skin_temperature_min', 'skin_temperature_max',
            'soil_temperature_level_1_min', 'soil_temperature_level_1_max',
            # Temperature (Celsius) - converted versions
            'temperature_2m_celsius', 'dewpoint_temperature_2m_celsius',
            'skin_temperature_celsius', 'soil_temperature_level_1_celsius',
            # Temperature min/max (Celsius) - converted versions
            'temperature_2m_min_celsius', 'temperature_2m_max_celsius',
            'dewpoint_temperature_2m_min_celsius', 'dewpoint_temperature_2m_max_celsius',
            'skin_temperature_min_celsius', 'skin_temperature_max_celsius',
            'soil_temperature_level_1_min_celsius', 'soil_temperature_level_1_max_celsius',
            # Precipitation & water balance (meters)
            'total_precipitation_sum', 'total_evaporation_sum', 'potential_evaporation_sum',
            'runoff_sum', 'surface_runoff_sum',
            # Precipitation & water balance (L/m²) - converted versions
            'total_precipitation_sum_lm2', 'total_evaporation_sum_lm2',
            'potential_evaporation_sum_lm2', 'runoff_sum_lm2', 'surface_runoff_sum_lm2',
            'snowfall_sum_lm2',
            # Soil moisture
            'volumetric_soil_water_layer_1', 'volumetric_soil_water_layer_2',
            'volumetric_soil_water_layer_3', 'volumetric_soil_water_layer_4',
            # Radiation
            'surface_solar_radiation_downwards_sum', 'surface_net_solar_radiation_sum',
            'surface_latent_heat_flux_sum',
            # Wind & pressure
            'u_component_of_wind_10m', 'v_component_of_wind_10m', 'surface_pressure',
            # Snow
            'snow_depth_water_equivalent', 'snowfall_sum'
        }

        # CHIRPS precipitation dataset variables
        self.chirps_variables = {
            'precipitation'  # Daily precipitation in mm
        }

        # VIIRS nighttime lights variables (monthly composites)
        self.viirs_variables = {
            'avg_rad',  # Average radiance, nW/cm2/sr
            'cf_cvg'    # Cloud-free observations used
        }

        # VIIRS Black Marble daily variables (VNP46A2). Short names, because
        # the band names of the product ('DNB_BRDF_Corrected_NTL') are not
        # usable as an index= argument
        self.viirs_daily_variables = {
            'ntl',               # DNB_BRDF_Corrected_NTL
            'ntl_gap_filled',    # Gap_Filled_DNB_BRDF_Corrected_NTL
            'lunar_irradiance',  # DNB_Lunar_Irradiance
            'quality_flag',      # Mandatory_Quality_Flag
            'cloud_mask',        # QF_Cloud_Mask
            'snow_flag'          # Snow_Flag
        }

        # DMSP-OLS nighttime lights variables (annual composites)
        self.dmsp_variables = {
            'avg_vis',           # Average visible band, 6-bit DN
            'stable_lights',     # Lights cleaned of ephemeral events
            'avg_lights_x_pct',  # avg_vis weighted by the share of lit nights
            'cf_cvg'             # Cloud-free observations used
        }

        # Final sensor-to-indices mapping
        self.sensor_indices = {
            'S2': (self.optical_indices | self.s2_exclusive_indices
                   | self.swir_optical_indices | self.raw_band_indices
                   | self.s2_rededge_bands),
            'Landsat': (self.optical_indices | self.swir_optical_indices
                        | self.raw_band_indices | self.thermal_indices),
            'MODIS': (self.optical_indices | self.swir_optical_indices
                      | self.raw_band_indices | {'lst'}),
            'S1': self.s1_indices,
            'S3': self.optical_indices | self.s3_exclusive_indices,
            'ERA5': self.era5_variables,
            'CHIRPS': self.chirps_variables,
            'VIIRS': self.viirs_variables,
            'VIIRS_DAILY': self.viirs_daily_variables,
            'DMSP': self.dmsp_variables
        }
        
        # Validate satellite
        if sat not in self.sensor_indices:
            available_sats = list(self.sensor_indices.keys())
            raise ValueError(f"Satellite '{sat}' is not supported. Available satellites are: {available_sats}")
        
        self.sat = sat

        # Factor that turns the stored band values into surface reflectance
        # (0-1), used by get_raw_band. Every optical collection is already in
        # reflectance: Landsat through scale_OLI / scale_ETM, Sentinel-2 and
        # MODIS through _to_reflectance (before 1.6.0 these two stayed as
        # integers multiplied by 10000 and this factor was 1e-4 for them)
        self.reflectance_scale = 1.0

        # Validate index for selected satellite
        available_indices = self.sensor_indices[self.sat]
        if index not in available_indices:
            raise ValueError(
                f"Index '{index}' is not available for {sat}. "
                f"Available indices for {sat} are: {sorted(list(available_indices))}"
            )
        
        self.index = index
        
        # Complete dictionary of index calculation methods
        self.d = {
            # Basic optical indices 
            'ndvi': self.get_ndvi, 'ndwi': self.get_ndwi, 'mndwi': self.get_mndwi, 
            'evi': self.get_evi, 'savi': self.get_savi, 'gndvi': self.get_gndvi, 
            'avi': self.get_avi, 'nbri': self.get_nbri, 'ndsi': self.get_ndsi, 
            'aweinsh': self.get_aweinsh, 'awei': self.get_awei, 'ndmi': self.get_ndmi,
            'lst': self.get_lst, 'vci': self.get_vci, 'utfvi': self.get_utfvi, 
            'nbr': self.get_nbr, 'wi2015': self.get_wi2015, 'ndbi': self.get_ndbi,
            
            # Additional optical indices 
            'msi': self.get_msi, 'nmi': self.get_nmi, 'ndti': self.get_ndti,
            'cri1': self.get_cri1, 'cri2': self.get_cri2, 'lai': self.get_lai, 
            'pri': self.get_pri, 'wdrvi': self.get_wdrvi,
            
            # Sentinel-2 specific indices 
            'ireci': self.get_ireci, 'mcari': self.get_mcari, 'reip': self.get_reip,
            'psri': self.get_psri, 'ndre': self.get_ndre, 'cig': self.get_cig,
            'cire': self.get_cire, 'mtci': self.get_mtci, 's2rep': self.get_s2rep,
            'ndci': self.get_ndci,

            # FAI - Floating Algae Index (S2, Landsat, MODIS)
            'fai': self.get_fai,

            # Raw reflectance bands (S2, Landsat, MODIS)
            'blue': lambda image: self.get_raw_band(image, 'Blue'),
            'green': lambda image: self.get_raw_band(image, 'Green'),
            'red': lambda image: self.get_raw_band(image, 'Red'),
            'nir': lambda image: self.get_raw_band(image, 'Nir'),
            'swir1': lambda image: self.get_raw_band(image, 'Swir1'),
            'swir2': lambda image: self.get_raw_band(image, 'Swir2'),
            'red_edge1': lambda image: self.get_raw_band(image, 'Red_Edge1'),
            'red_edge2': lambda image: self.get_raw_band(image, 'Red_Edge2'),
            'red_edge3': lambda image: self.get_raw_band(image, 'Red_Edge3'),

            # Sentinel-3 OLCI specific indices
            'oci': self.get_oci, 'tsi': self.get_tsi, 'cdom': self.get_cdom,
            'turbidity': self.get_turbidity, 'spm': self.get_spm, 'kd490': self.get_kd490,
            'floating_algae': self.get_floating_algae, 'red_edge_position': self.get_red_edge_position,
            'fluorescence_height': self.get_fluorescence_height, 
            'water_leaving_reflectance': self.get_water_leaving_reflectance,

            # SAR indices - usando lambdas para pasar el parámetro normalize
            'rvi': lambda image: self.get_rvi(image, normalize=self.normalize_sar),
            'vv': lambda image: self.get_vv(image, normalize=self.normalize_sar),
            'vh': lambda image: self.get_vh(image, normalize=self.normalize_sar),
            'vv_vh_ratio': lambda image: self.get_vv_vh_ratio(image, normalize=self.normalize_sar),
            'dpsvi': lambda image: self.get_dpsvi(image, normalize=self.normalize_sar),
            'rfdi': lambda image: self.get_rfdi(image, normalize=self.normalize_sar),
            'vsdi': lambda image: self.get_vsdi(image, normalize=self.normalize_sar),

            # ERA5-Land climate reanalysis variables
            # Temperature
            'temperature_2m': self.get_era5_temperature_2m,
            'dewpoint_temperature_2m': self.get_era5_dewpoint_temperature_2m,
            'skin_temperature': self.get_era5_skin_temperature,
            'soil_temperature_level_1': self.get_era5_soil_temperature_level_1,
            # Precipitation & water balance
            'total_precipitation_sum': self.get_era5_total_precipitation_sum,
            'total_evaporation_sum': self.get_era5_total_evaporation_sum,
            'potential_evaporation_sum': self.get_era5_potential_evaporation_sum,
            'runoff_sum': self.get_era5_runoff_sum,
            'surface_runoff_sum': self.get_era5_surface_runoff_sum,
            # Soil moisture
            'volumetric_soil_water_layer_1': self.get_era5_volumetric_soil_water_layer_1,
            'volumetric_soil_water_layer_2': self.get_era5_volumetric_soil_water_layer_2,
            'volumetric_soil_water_layer_3': self.get_era5_volumetric_soil_water_layer_3,
            'volumetric_soil_water_layer_4': self.get_era5_volumetric_soil_water_layer_4,
            # Radiation
            'surface_solar_radiation_downwards_sum': self.get_era5_surface_solar_radiation_downwards_sum,
            'surface_net_solar_radiation_sum': self.get_era5_surface_net_solar_radiation_sum,
            'surface_latent_heat_flux_sum': self.get_era5_surface_latent_heat_flux_sum,
            # Wind & pressure
            'u_component_of_wind_10m': self.get_era5_u_component_of_wind_10m,
            'v_component_of_wind_10m': self.get_era5_v_component_of_wind_10m,
            'surface_pressure': self.get_era5_surface_pressure,
            # Snow
            'snow_depth_water_equivalent': self.get_era5_snow_depth_water_equivalent,
            'snowfall_sum': self.get_era5_snowfall_sum,
            # Temperature in Celsius (converted)
            'temperature_2m_celsius': self.get_era5_temperature_2m_celsius,
            'dewpoint_temperature_2m_celsius': self.get_era5_dewpoint_temperature_2m_celsius,
            'skin_temperature_celsius': self.get_era5_skin_temperature_celsius,
            'soil_temperature_level_1_celsius': self.get_era5_soil_temperature_level_1_celsius,
            # Temperature min/max (Kelvin)
            'temperature_2m_min': self.get_era5_temperature_2m_min,
            'temperature_2m_max': self.get_era5_temperature_2m_max,
            'dewpoint_temperature_2m_min': self.get_era5_dewpoint_temperature_2m_min,
            'dewpoint_temperature_2m_max': self.get_era5_dewpoint_temperature_2m_max,
            'skin_temperature_min': self.get_era5_skin_temperature_min,
            'skin_temperature_max': self.get_era5_skin_temperature_max,
            'soil_temperature_level_1_min': self.get_era5_soil_temperature_level_1_min,
            'soil_temperature_level_1_max': self.get_era5_soil_temperature_level_1_max,
            # Temperature min/max (Celsius)
            'temperature_2m_min_celsius': self.get_era5_temperature_2m_min_celsius,
            'temperature_2m_max_celsius': self.get_era5_temperature_2m_max_celsius,
            'dewpoint_temperature_2m_min_celsius': self.get_era5_dewpoint_temperature_2m_min_celsius,
            'dewpoint_temperature_2m_max_celsius': self.get_era5_dewpoint_temperature_2m_max_celsius,
            'skin_temperature_min_celsius': self.get_era5_skin_temperature_min_celsius,
            'skin_temperature_max_celsius': self.get_era5_skin_temperature_max_celsius,
            'soil_temperature_level_1_min_celsius': self.get_era5_soil_temperature_level_1_min_celsius,
            'soil_temperature_level_1_max_celsius': self.get_era5_soil_temperature_level_1_max_celsius,
            # Precipitation in L/m² (converted)
            'total_precipitation_sum_lm2': self.get_era5_total_precipitation_sum_lm2,
            'total_evaporation_sum_lm2': self.get_era5_total_evaporation_sum_lm2,
            'potential_evaporation_sum_lm2': self.get_era5_potential_evaporation_sum_lm2,
            'runoff_sum_lm2': self.get_era5_runoff_sum_lm2,
            'surface_runoff_sum_lm2': self.get_era5_surface_runoff_sum_lm2,
            'snowfall_sum_lm2': self.get_era5_snowfall_sum_lm2,
            # CHIRPS precipitation
            'precipitation': self.get_chirps_precipitation,
            # Nighttime lights. 'cf_cvg' is a band of both VIIRS and DMSP and
            # means the same thing in each, so one method serves both
            'avg_rad': self.get_viirs_avg_rad,
            'cf_cvg': self.get_cf_cvg,
            'ntl': self.get_viirs_ntl,
            'ntl_gap_filled': self.get_viirs_ntl_gap_filled,
            'lunar_irradiance': self.get_viirs_lunar_irradiance,
            'quality_flag': self.get_viirs_quality_flag,
            'cloud_mask': self.get_viirs_cloud_mask,
            'snow_flag': self.get_viirs_snow_flag,
            'avg_vis': self.get_dmsp_avg_vis,
            'stable_lights': self.get_dmsp_stable_lights,
            'avg_lights_x_pct': self.get_dmsp_avg_lights_x_pct
        }

        # Generate dynamic temporal periods - replaces all hardcoded periods
        self.period_dates, self.period_names = self._generate_periods(periods)
        
        # Initialize satellite collections with appropriate configuration
        self._setup_satellite_collections()
    
    def _generate_periods(self, n_periods):
        """
        Dynamically generate period dates and names based on number of periods.
        
        Creates temporal divisions of the year for composite generation. Handles
        common cases (seasonal, monthly, bi-monthly) with meaningful names and
        provides generic equal division for custom period counts.
        
        Parameters
        ----------
        n_periods : int
            Number of periods to divide the year into:
            
            * ``4`` : Traditional seasons
            * ``12`` : Calendar months
            * ``24`` : Bi-monthly (~15 days)
            * Other : Equal day-of-year division
        
        Returns
        -------
        period_dates : list of list
            Date ranges as [start, end] pairs in '-MM-DD' format.
            Example: ``[['-01-01', '-03-31'], ['-04-01', '-06-30'], ...]``
            
        period_names : list of str
            Descriptive names for each period.
            Example: ``['winter', 'spring', 'summer', 'autumn']``
        
        Notes
        -----
        Date formats are designed to be prepended with year strings.
        February always uses 28 days to avoid leap year complications.
        Custom periods use day-of-year calculation with 365-day year.
        
        Examples
        --------
        Generate seasonal periods::
        
            >>> dates, names = self._generate_periods(4)
            >>> print(dates[0])  # Winter period
            ['-01-01', '-03-31']
            >>> print(names[0])
            'winter'
        
        Generate custom 8-period division::
        
            >>> dates, names = self._generate_periods(8)
            >>> len(dates)  # 8 equal periods
            8
            >>> names[0]  # Generic naming
            'p1'
        
        Use in date filtering::
        
            >>> year = 2020
            >>> start_date = str(year) + dates[0][0]  # '2020-01-01'
            >>> end_date = str(year) + dates[0][1]    # '2020-03-31'
        
        See Also
        --------
        get_period_composite : Uses generated periods for filtering
        get_year_composite : Processes all periods for each year
        """
        if n_periods == 4:
            # Traditional seasons (Northern Hemisphere meteorological seasons)
            period_dates = [
                ['-01-01', '-03-31'],  # Winter (simplified as Jan-Mar)
                ['-04-01', '-06-30'],  # Spring
                ['-07-01', '-09-30'],  # Summer
                ['-10-01', '-12-31']   # Autumn
            ]
            period_names = ['winter', 'spring', 'summer', 'autumn']
            
        elif n_periods == 12:
            # Monthly periods - using fixed days to avoid leap year complications
            period_dates = [
                ['-01-01', '-01-31'],  # January
                ['-02-01', '-02-28'],  # February (always 28 days - works for all years)
                ['-03-01', '-03-31'],  # March
                ['-04-01', '-04-30'],  # April
                ['-05-01', '-05-31'],  # May
                ['-06-01', '-06-30'],  # June
                ['-07-01', '-07-31'],  # July
                ['-08-01', '-08-31'],  # August
                ['-09-01', '-09-30'],  # September
                ['-10-01', '-10-31'],  # October
                ['-11-01', '-11-30'],  # November
                ['-12-01', '-12-31']   # December
            ]
            period_names = ['january', 'february', 'march', 'april', 'may', 'june',
                        'july', 'august', 'september', 'october', 'november', 'december']
        
        elif n_periods == 24:
            # Bi-monthly periods (every ~15 days) - fixed dates to avoid leap year complications
            period_dates = [
                ['-01-01', '-01-15'], ['-01-16', '-01-31'],  # January (2 periods)
                ['-02-01', '-02-15'], ['-02-16', '-02-28'],  # February (always 28 days)
                ['-03-01', '-03-15'], ['-03-16', '-03-31'],  # March (2 periods)
                ['-04-01', '-04-15'], ['-04-16', '-04-30'],  # April (2 periods)
                ['-05-01', '-05-15'], ['-05-16', '-05-31'],  # May (2 periods)
                ['-06-01', '-06-15'], ['-06-16', '-06-30'],  # June (2 periods)
                ['-07-01', '-07-15'], ['-07-16', '-07-31'],  # July (2 periods)
                ['-08-01', '-08-15'], ['-08-16', '-08-31'],  # August (2 periods)
                ['-09-01', '-09-15'], ['-09-16', '-09-30'],  # September (2 periods)
                ['-10-01', '-10-15'], ['-10-16', '-10-31'],  # October (2 periods)
                ['-11-01', '-11-15'], ['-11-16', '-11-30'],  # November (2 periods)
                ['-12-01', '-12-15'], ['-12-16', '-12-31']   # December (2 periods)
            ]
            period_names = [f'p{i+1}' for i in range(24)]
        
        else:
            # Generic periods - divide year equally using day-of-year approach
            period_dates = []
            period_names = []
            days_per_period = 365 // n_periods
            
            for i in range(n_periods):
                # Calculate start and end day of year (1-365)
                start_day = i * days_per_period + 1
                if i == n_periods - 1:  # Last period goes to end of year
                    end_day = 365
                else:
                    end_day = (i + 1) * days_per_period
                
                # Convert day of year to month-day using non-leap year (2021)
                # Using 2021 because it's not a leap year, ensuring consistent 365-day calendar
                start_date = datetime(2021, 1, 1) + timedelta(days=start_day - 1)
                end_date = datetime(2021, 1, 1) + timedelta(days=end_day - 1)
                
                # Format as strings in '-MM-DD' format
                start_str = f'-{start_date.month:02d}-{start_date.day:02d}'
                end_str = f'-{end_date.month:02d}-{end_date.day:02d}'
                
                period_dates.append([start_str, end_str])
                period_names.append(f'p{i+1}')
        
        return period_dates, period_names
    
    # Clouds para que os quiero...
    def mask_s2_clouds(self, image):
        """
        Mask clouds and shadows in Sentinel-2 images using QA60 band.
        
        Parameters
        ----------
        image : ee.Image
            Sentinel-2 SR image with QA60 band
            
        Returns
        -------
        ee.Image
            Cloud-masked image
            
        References
        ----------
        Sentinel-2 Cloud Masking with s2cloudless
        https://developers.google.com/earth-engine/tutorials/community/sentinel-2-s2cloudless
        """
        # QA60 es la banda de calidad de Sentinel-2
        qa = image.select('QA60')
        
        # Bits 10 y 11 son nubes y cirrus, respectivamente
        cloud_bit_mask = 1 << 10
        cirrus_bit_mask = 1 << 11
        
        # Crear máscara: 0 donde hay nubes/cirrus
        mask = qa.bitwiseAnd(cloud_bit_mask).eq(0).And(
            qa.bitwiseAnd(cirrus_bit_mask).eq(0))
        
        # Aplicar máscara y copiar propiedades
        return image.updateMask(mask).copyProperties(
            image, ['system:time_start', 'system:time_end', 'system:index'])

    def mask_s2_scl(self, image):
        """
        Mask clouds, cloud shadows and cirrus in Sentinel-2 images using the
        Scene Classification Layer (SCL band).

        SCL provides per-pixel classification at 20m resolution, offering
        more accurate cloud/shadow detection than the QA60 band (which is
        coarser and misses thin clouds and shadows).

        Masked classes:
        - 1: Saturated / defective pixels
        - 3: Cloud shadows
        - 8: Medium probability cloud
        - 9: High probability cloud
        - 10: Thin cirrus

        Parameters
        ----------
        image : ee.Image
            Sentinel-2 SR image with SCL band (S2_SR_HARMONIZED collection).

        Returns
        -------
        ee.Image
            Cloud/shadow-masked image with SCL band removed.

        References
        ----------
        European Space Agency (2021). Sentinel-2 MSI Level-2A Algorithm
        Theoretical Basis Document. ESA-EOPG-GSCB-TN-0001.
        """
        scl = image.select('SCL')
        mask = (scl.neq(1)    # saturated/defective
                .And(scl.neq(3))   # cloud shadows
                .And(scl.neq(8))   # medium probability cloud
                .And(scl.neq(9))   # high probability cloud
                .And(scl.neq(10))) # thin cirrus
        return image.updateMask(mask).copyProperties(
            image, ['system:time_start', 'system:time_end', 'system:index'])

    def mask_viirs_quality(self, image):
        """
        Keep only the good retrievals of a VIIRS Black Marble daily image.

        ``Mandatory_Quality_Flag`` grades every pixel of VNP46A2:

        - 0: high quality, persistent nighttime lights
        - 1: high quality, ephemeral nighttime lights
        - 2: poor quality, outlier or likely cloud contaminated
        - 255: no retrieval

        Grades 0 and 1 are kept. Flag 2 is where a cloud the mask missed shows
        up as a dark or a bright night, which is exactly the noise that ruins
        a map of the date of the maximum.

        Parameters
        ----------
        image : ee.Image
            VNP46A2 image with its quality bands.

        Returns
        -------
        ee.Image
            The same image with the poor and missing retrievals masked out.

        References
        ----------
        Román, M.O., Wang, Z., Sun, Q. et al. (2018). NASA's Black Marble
        nighttime lights product suite. Remote Sensing of Environment, 210,
        113-143. https://doi.org/10.1016/j.rse.2018.03.017
        """
        quality = image.select('Mandatory_Quality_Flag')
        return image.updateMask(quality.lte(1)).copyProperties(
            image, ['system:time_start', 'system:time_end', 'system:index'])

    def mask_landsat_clouds(self, image):
        """
        Mask clouds and shadows in Landsat Collection 2 images using QA_PIXEL band.
        
        Parameters
        ----------
        image : ee.Image
            Landsat Collection 2 Level-2 image with QA_PIXEL band
            
        Returns
        -------
        ee.Image
            Cloud-masked image
            
        References
        ----------
        Landsat Collection 2 Level-2 Science Products
        https://www.usgs.gov/landsat-missions/landsat-collection-2-level-2-science-products
        """
        # QA_PIXEL contiene información de calidad
        qa = image.select('QA_PIXEL')
        
        # Bits para diferentes tipos de nubes/sombras en Collection 2
        cloud_bit = 3           # Cloud
        cloud_shadow_bit = 4    # Cloud Shadow
        cirrus_bit = 2         # Cirrus (solo Landsat 8/9)
        
        # Crear máscaras para cada condición
        cloud_mask = qa.bitwiseAnd(1 << cloud_bit).eq(0)
        shadow_mask = qa.bitwiseAnd(1 << cloud_shadow_bit).eq(0)
        
        # Para Landsat 8/9, también filtrar cirrus
        # Detectar si es L8/9 por la presencia de ST_B10
        is_oli = image.bandNames().contains('ST_B10')
        
        mask = ee.Algorithms.If(
            is_oli,
            # Landsat 8/9: incluir máscara de cirrus
            cloud_mask.And(shadow_mask).And(qa.bitwiseAnd(1 << cirrus_bit).eq(0)),
            # Landsat 4/5/7: solo nubes y sombras
            cloud_mask.And(shadow_mask)
        )
        
        # Aplicar máscara y copiar propiedades
        return image.updateMask(mask).copyProperties(
        image, ['system:time_start', 'system:time_end', 'system:index'])
    
    def _setup_satellite_collections(self):
        """
        Configure satellite image collections based on selected sensor.
        
        Initializes the appropriate Earth Engine ImageCollection with necessary
        preprocessing, band selection, scaling, and filtering. Each sensor
        requires specific handling due to different band configurations,
        scaling factors, and data characteristics.
        
        Notes
        -----
        Collection configurations:
        
        **Landsat (merged):**
            * Collections: L4, L5, L7, L8, L9 Collection 2 Level-2
            * Resolution: 30m
            * Scaling: Applied via scale_OLI/scale_ETM functions
            * Thermal bands: ST_B10 (L8/9), ST_B6 (L4/5/7)
            * Cloud filtering: Optional using QA_PIXEL band
            
        **Sentinel-2:**
            * Collection: COPERNICUS/S2_SR_HARMONIZED
            * Resolution: 10-20m
            * Bands: B2-B8, B11-B12 (includes Red Edge)
            * Coverage: 2015-present
            * Cloud filtering: Optional using QA60 band
            
        **MODIS:**
            * Products: MOD09A1 (reflectance), MOD11A1/MYD11A1 (LST)
            * Resolution: 500m (reflectance), 1km (LST)
            * Temporal: 8-day composites
            * Smart Terra/Aqua selection based on date range
            
        **Sentinel-1:**
            * Product: GRD (Ground Range Detected)
            * Mode: IW (Interferometric Wide)
            * Polarization: VV+VH dual-pol
            * Optional ARD preprocessing with terrain correction
            
        **Sentinel-3:**
            * Instrument: OLCI (Ocean and Land Color)
            * Bands: 16 spectral bands (400-1020nm)
            * Resolution: 300m
            * Focus: Ocean/coastal applications
        
        Raises
        ------
        ValueError
            If satellite sensor is not recognized.
        
        Warnings
        --------
        * Landsat 7 has SLC-off gaps after May 2003
        * MODIS Aqua only available after July 2002
        * Sentinel-3 thermal bands not available in Earth Engine
        
        See Also
        --------
        scale_OLI : Landsat 8-9 scaling
        scale_ETM : Landsat 4-5-7 scaling
        S1ARDProcessor : Advanced SAR preprocessing
        mask_s2_clouds : Sentinel-2 cloud masking
        mask_landsat_clouds : Landsat cloud masking
        """
        
        # ============= LANDSAT CONFIGURATION =============
        # Initialize Landsat collections
        LC09col = ee.ImageCollection("LANDSAT/LC09/C02/T1_L2").filterBounds(self.roi) 
        LC08col = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2").filterBounds(self.roi) 
        LE07col = ee.ImageCollection("LANDSAT/LE07/C02/T1_L2").filterBounds(self.roi) 
        LT05col = ee.ImageCollection("LANDSAT/LT05/C02/T1_L2").filterBounds(self.roi) 
        LT04col = ee.ImageCollection("LANDSAT/LT04/C02/T1_L2").filterBounds(self.roi) 
        
        # Apply cloud cover metadata filter if enabled for Landsat
        if self.sat == 'Landsat' and self.cloud_filter:
            print(f"Applying cloud filter to Landsat: max {self.max_cloud_cover}% cloud cover")
            LC09col = LC09col.filter(ee.Filter.lte('CLOUD_COVER', self.max_cloud_cover))
            LC08col = LC08col.filter(ee.Filter.lte('CLOUD_COVER', self.max_cloud_cover))
            LE07col = LE07col.filter(ee.Filter.lte('CLOUD_COVER', self.max_cloud_cover))
            LT05col = LT05col.filter(ee.Filter.lte('CLOUD_COVER', self.max_cloud_cover))
            LT04col = LT04col.filter(ee.Filter.lte('CLOUD_COVER', self.max_cloud_cover))
        
        # Merge Landsat collections
        OLI = LC09col.merge(LC08col)  # Landsat 8/9
        ETM = LE07col.merge(LT05col).merge(LT04col)  # Landsat 4/5/7
        
        # Scaling function for OLI (L8/9) with thermal band and optional cloud masking
        def scale_OLI_with_thermal(image):
            optical = image.select(['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7']).multiply(0.0000275).add(-0.2).rename(['Blue', 'Green', 'Red', 'Nir', 'Swir1', 'Swir2'])
            thermal = image.select(['ST_B10'])
            scaled = image.addBands(optical, None, True).addBands(thermal, None, True)
            
            # Apply pixel-level cloud mask if enabled
            if self.cloud_filter and self.sat == 'Landsat':
                scaled = self.mask_landsat_clouds(scaled)
            
            return scaled
        
        # Scaling function for OLI (L8/9) with thermal band
        def scale_OLI_with_thermal(image):
            optical = image.select(['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7']).multiply(0.0000275).add(-0.2).rename(['Blue', 'Green', 'Red', 'Nir', 'Swir1', 'Swir2'])
            thermal = image.select(['ST_B10'])
            return image.addBands(optical, None, True).addBands(thermal, None, True)

        # Scaling function for ETM/TM (L4/5/7) with thermal band
        def scale_ETM_with_thermal(image):
            optical = image.select(['SR_B1','SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7']).multiply(0.0000275).add(-0.2).rename(['Blue', 'Green', 'Red', 'Nir', 'Swir1', 'Swir2'])
            thermal = image.select(['ST_B6'])
            return image.addBands(optical, None, True).addBands(thermal, None, True)

        OLI_ = OLI.map(scale_OLI_with_thermal)
        ETM_ = ETM.map(scale_ETM_with_thermal)

        # Aplicar enmascaramiento de nubes DESPUÉS del escalado si está activado
        if self.sat == 'Landsat' and hasattr(self, 'cloud_filter') and self.cloud_filter:
            OLI_ = OLI_.map(self.mask_landsat_clouds)
            ETM_ = ETM_.map(self.mask_landsat_clouds)

        Landsat = OLI_.merge(ETM_)
        
        # ============= SENTINEL-2 CONFIGURATION =============
        if self.sat == 'S2' and self.cloud_filter:
            print(f"Applying cloud filter to Sentinel-2: max {self.max_cloud_cover}% cloud cover")
            # First filter by scene-level cloud cover metadata
            S2col = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED").filterBounds(self.roi).filter(
                ee.Filter.lte('CLOUDY_PIXEL_PERCENTAGE', self.max_cloud_cover)
            )
            # Then apply pixel-level cloud mask (SCL preferred, QA60 as fallback)
            if self.scl_mask:
                print("Applying pixel-level cloud/shadow mask using SCL band")
                S2col = S2col.map(self.mask_s2_scl).select([
                    'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B11', 'B12'
                ], [
                    'Blue', 'Green', 'Red', 'Red_Edge1', 'Red_Edge2', 'Red_Edge3',
                    'Nir', 'Swir1', 'Swir2'
                ])
            else:
                print("Applying pixel-level cloud mask using QA60 band")
                S2col = S2col.map(self.mask_s2_clouds).select([
                    'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B11', 'B12'
                ], [
                    'Blue', 'Green', 'Red', 'Red_Edge1', 'Red_Edge2', 'Red_Edge3',
                    'Nir', 'Swir1', 'Swir2'
                ])
        else:
            # No cloud filtering
            S2col = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED").select([
                'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B11', 'B12'
            ], [
                'Blue', 'Green', 'Red', 'Red_Edge1', 'Red_Edge2', 'Red_Edge3',
                'Nir', 'Swir1', 'Swir2'
            ]).filterBounds(self.roi)

        # Sentinel-2 L2A stores surface reflectance as integers multiplied by
        # 10000. Ratio indices such as NDVI do not notice, but every index with
        # a constant, a sum or an inverse does (SAVI, EVI, AVI, WI2015, CRI...),
        # so the bands are brought to reflectance in 0-1 like Landsat's
        S2col = S2col.map(_to_reflectance)
        
        # ============= MODIS CONFIGURATION =============
        # MODIS Terra + Aqua - Smart configuration based on time period
        # MOD11A1 (Terra): Feb 2000 - present
        # MYD11A1 (Aqua): Jul 2002 - present
        
        aqua_start_date = ee.Date('2002-07-04')
        period_start = ee.Date(f'{self.start_year}-01-01')
        use_aqua = period_start.millis().gte(aqua_start_date.millis())
        
        # Terra reflectance (always available)
        # Reflectance is stored as integers multiplied by 10000, as in Sentinel-2
        MOD09A1 = ee.ImageCollection("MODIS/061/MOD09A1").select(
            ['sur_refl_b03', 'sur_refl_b04', 'sur_refl_b01', 'sur_refl_b02', 'sur_refl_b06', 'sur_refl_b07'], 
            ['Blue', 'Green', 'Red', 'Nir', 'Swir1', 'Swir2']
        ).filterBounds(self.roi).map(_to_reflectance)
        
        # Terra LST (always available since 2000)
        MOD11A1 = ee.ImageCollection("MODIS/061/MOD11A1").select(['LST_Day_1km']).filterBounds(self.roi)
        
        # Conditional configuration based on time period
        if self.start_year >= 2003:  # Period after Aqua availability
            if self.sat == 'MODIS':
                print("Using MODIS Terra + Aqua LST (maximum coverage)")
            
            # Aqua reflectance
            MYD09A1 = ee.ImageCollection("MODIS/061/MYD09A1").select(
                ['sur_refl_b03', 'sur_refl_b04', 'sur_refl_b01', 'sur_refl_b02', 'sur_refl_b06', 'sur_refl_b07'], 
                ['Blue', 'Green', 'Red', 'Nir', 'Swir1', 'Swir2']
            ).filterBounds(self.roi).map(_to_reflectance)
            
            # Aqua LST
            MYD11A1 = ee.ImageCollection("MODIS/061/MYD11A1").select(['LST_Day_1km']).filterBounds(self.roi)
            
            # Combine both satellites
            MODIS_reflectance = MOD09A1.merge(MYD09A1)
            MODIS_LST = MOD11A1.merge(MYD11A1)
            
        else:  # Period before Aqua (2000-2002)
            if self.sat == 'MODIS':
                print("Using MODIS Terra LST only (Aqua not available before July 2002)")
            
            # Only Terra
            MODIS_reflectance = MOD09A1
            MODIS_LST = MOD11A1
        
        # Function to merge reflectance with LST from same day
        def merge_modis_lst(ref_image):
            date = ref_image.date()
            # Find LST from same day (±1 day tolerance)
            lst_same_day = MODIS_LST.filterDate(
                date.advance(-1, 'day'), 
                date.advance(1, 'day')
            )
            
            # If multiple LST observations on same day, take mean
            lst_composite = lst_same_day.mean()
            
            # Only add LST if valid data exists
            return ee.Algorithms.If(
                lst_same_day.size().gt(0),
                ref_image.addBands(lst_composite),
                ref_image.addBands(ee.Image.constant(-9999).rename('LST_Day_1km').updateMask(ee.Image.constant(0)))
            )
        
        MOD09A1_with_LST = MODIS_reflectance.map(merge_modis_lst)
        
        # ============= SENTINEL-1 SAR CONFIGURATION =============
        s1 = ee.ImageCollection('COPERNICUS/S1_GRD').filter(
            ee.Filter.listContains('transmitterReceiverPolarisation', 'VH')
        ).filter(
            ee.Filter.listContains('transmitterReceiverPolarisation', 'VV')
        ).filter(ee.Filter.eq('instrumentMode', 'IW'))

        # Apply orbit filter
        if self.orbit == 'ASCENDING':
            s1 = s1.filter(ee.Filter.eq('orbitProperties_pass', 'ASCENDING'))
            if self.sat == 'S1':
                print("Using Sentinel-1 ascending orbits only.")
        elif self.orbit == 'DESCENDING':
            s1 = s1.filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))
            if self.sat == 'S1':
                print("Using Sentinel-1 descending orbits only.")
        else:
            if self.sat == 'S1':
                print("Using all Sentinel-1 orbits (ascending + descending).")

        # Filter by ROI
        s1 = s1.filterBounds(self.roi)

        # Define preprocessing function
        if hasattr(self, 'use_sar_ard') and self.use_sar_ard:
            if self.sat == 'S1':
                print(f"Applying S1 ARD preprocessing:")
                print(f"  - Speckle filter: {self.sar_speckle_filter}")
                print(f"  - Terrain correction: {self.sar_terrain_correction}")
                if self.sar_terrain_correction:
                    print(f"  - Terrain model: {self.sar_terrain_model}")
            
            # Create ARD processor
            ard_processor = S1ARDProcessor(
                speckle_filter=self.sar_speckle_filter,
                speckle_filter_kernel_size=7,
                terrain_correction=self.sar_terrain_correction,
                terrain_flattening_model=self.sar_terrain_model,
                dem='COPERNICUS_30',
                # COPERNICUS/S1_GRD comes in dB: the processor converts it to
                # linear power for terrain correction and speckle filtering and
                # returns dB, the usual scale for VV/VH composites. Indices
                # convert back to linear themselves (see _s1_linear)
                input_format='DB',
                format='DB'
            )
            
            def apply_speckle_filter(image):
                return ard_processor.process_image(image)

            selected_bands = ['VV', 'VH', 'angle']
        else:
            if self.sat == 'S1':
                print("Using basic S1 preprocessing (focal_median filter only)")

            def apply_speckle_filter(image):
                filtered = image.focal_median(radius=1, kernelType='square', units='pixels')
                return filtered.copyProperties(image, ['system:time_start', 'system:time_end'])

            selected_bands = ['VV', 'VH']

        # Apply the selected preprocessing
        s1S1 = s1.select(selected_bands).map(apply_speckle_filter)
        
        # ============= SENTINEL-3 OLCI CONFIGURATION =============
        S3col = ee.ImageCollection("COPERNICUS/S3/OLCI").select([
            'Oa01_radiance', 'Oa02_radiance', 'Oa03_radiance', 'Oa04_radiance',
            'Oa05_radiance', 'Oa06_radiance', 'Oa07_radiance', 'Oa08_radiance',
            'Oa09_radiance', 'Oa10_radiance', 'Oa11_radiance', 'Oa12_radiance',
            'Oa16_radiance', 'Oa17_radiance', 'Oa18_radiance', 'Oa21_radiance'
        ], [
            'Violet', 'Blue', 'Blue2', 'Blue_Green', 'Green', 'Green2',
            'Red', 'Red2', 'Red3', 'Red_Edge1', 'Red_Edge2', 'Nir',
            'Nir2', 'Nir3', 'Nir4', 'Nir5'
        ]).filterBounds(self.roi).map(_s3_to_reflectance)

        # ============= ERA5-LAND CLIMATE REANALYSIS CONFIGURATION =============
        # ERA5-Land daily aggregated climate variables (1950-present, ~11km resolution)
        ERA5col = ee.ImageCollection("ECMWF/ERA5_LAND/DAILY_AGGR").filterBounds(self.roi)

        # ============= CHIRPS PRECIPITATION CONFIGURATION =============
        # CHIRPS daily precipitation dataset (1981-present, ~5.5km resolution)
        CHIRPScol = ee.ImageCollection("UCSB-CHG/CHIRPS/DAILY").filterBounds(self.roi)

        # ============= NIGHTTIME LIGHTS CONFIGURATION =============
        # VIIRS DNB monthly composites, stray-light corrected (2014-present,
        # ~500m). One image per month, so periods=12 gives one image per period
        VIIRScol = ee.ImageCollection(
            "NOAA/VIIRS/DNB/MONTHLY_V1/VCMSLCFG").filterBounds(self.roi)

        # VIIRS Black Marble daily, BRDF corrected for the moon and the
        # atmosphere (2012-present, ~500m). NOAA/VIIRS/001/VNP46A2 is the
        # deprecated version of this same product, frozen at the end of 2024
        VIIRSdaily = ee.ImageCollection(
            "NASA/VIIRS/002/VNP46A2").filterBounds(self.roi)
        if self.sat == 'VIIRS_DAILY' and self.cloud_filter:
            VIIRSdaily = VIIRSdaily.map(self.mask_viirs_quality)

        # DMSP-OLS annual composites (1992-2013, ~1km). One image per year:
        # there is no intra-annual variation to composite
        DMSPcol = ee.ImageCollection(
            "NOAA/DMSP-OLS/NIGHTTIME_LIGHTS").filterBounds(self.roi)

        # ============= ASSIGN COLLECTIONS BASED ON SENSOR =============
        if self.sat == 'S2':
            self.ndvi_col = S2col
            if self.cloud_filter:
                print("Sentinel-2 collection configured with cloud filtering")
        elif self.sat == 'Landsat':
            self.ndvi_col = Landsat
            print("Landsat collection includes thermal bands: ST_B10 (L8/9) and ST_B6 (L4/5/7)")
            if self.cloud_filter:
                print("Landsat collection configured with cloud filtering")
        elif self.sat == 'MODIS':
            self.ndvi_col = MOD09A1_with_LST
            print("MODIS collection includes LST_Day_1km from Terra and Aqua satellites")
        elif self.sat == 'S1':
            self.ndvi_col = s1S1
        elif self.sat == 'S3':
            self.ndvi_col = S3col
            print("Sentinel-3 collection uses OLCI data only. LST not available (requires SLSTR).")
        elif self.sat == 'ERA5':
            self.ndvi_col = ERA5col
            print("ERA5-Land daily climate reanalysis (1950-present, ~11km resolution)")
        elif self.sat == 'CHIRPS':
            self.ndvi_col = CHIRPScol
            print("CHIRPS daily precipitation (1981-present, ~5.5km resolution, 50°S-50°N)")
        elif self.sat == 'VIIRS':
            self.ndvi_col = VIIRScol
            print("VIIRS DNB monthly nighttime lights (2014-present, ~500m resolution)")
            if self.periods != 12:
                print(f"  Note: VIIRS is monthly, so periods={self.periods} averages "
                      "several months per period; periods=12 keeps each month")
        elif self.sat == 'VIIRS_DAILY':
            self.ndvi_col = VIIRSdaily
            print("VIIRS Black Marble daily nighttime lights (2012-present, ~500m resolution)")
            if self.cloud_filter:
                print("Applying pixel-level quality mask using Mandatory_Quality_Flag")
            else:
                print("No quality mask: poor retrievals and gap-filled pixels are kept")
        elif self.sat == 'DMSP':
            self.ndvi_col = DMSPcol
            print("DMSP-OLS annual nighttime lights (1992-2013, ~1km resolution)")
            if self.periods != 1:
                print(f"  Note: DMSP is annual, so every period of a year reads the "
                      f"same image; periods={self.periods} cannot resolve timing "
                      "within a year")
        else:
            print('Not a valid satellite')
    
    def _period_date_range(self, year, period_idx):
        """
        Date range of a period, ready for ``filterDate``.

        ``filterDate`` excludes its end date, while ``period_dates`` stores the
        **last day** of each period, so filtering by it dropped that day from
        every composite: January ran to the 30th, and in a leap year February
        never saw the 29th, because the month is stored with 28 days.

        The end of a period is therefore taken as the start of the next one,
        and the end of the last period as the first day of the next year. That
        keeps the periods contiguous, restores the missing day, handles leap
        years on its own, and makes one-day periods (``periods=365``) work
        instead of building an empty range that Earth Engine rejects.

        Parameters
        ----------
        year : int
            Year the period belongs to.
        period_idx : int
            Zero-based index of the period.

        Returns
        -------
        (str, str)
            Start date (inclusive) and end date (exclusive), 'YYYY-MM-DD'.
        """
        init = f"{year}{self.period_dates[period_idx][0]}"
        if period_idx + 1 < len(self.period_dates):
            ends = f"{year}{self.period_dates[period_idx + 1][0]}"
        else:
            ends = f"{year + 1}{self.period_dates[0][0]}"
        return init, ends

    def get_period_composite(self, year, period_idx):
        """
        Generate composite image for a specific temporal period within a year.
        
        Creates a single composite by applying the configured statistical reducer
        to all satellite images within the defined temporal period. Core processing
        function for temporal composite generation.
        
        Parameters
        ----------
        year : int
            Target year for composite generation. Must be within satellite
            data availability range.
            
        period_idx : int
            Zero-based index of temporal period within the year.
            Must be less than self.periods.
            
        Returns
        -------
        ee.Image
            Single-band composite image with the selected index values, cast
            to float (to 32-bit integer for ``key='count'``). Band name depends
            on reducer type (see :meth:`_composite_band_name`). A period with
            no images still returns one band: fully masked, or 0 everywhere
            for ``key='count'``.

        Notes
        -----
        Processing workflow:
        
        1. Extract date range from ``self.period_dates[period_idx]``
        2. Filter satellite collection to date range
        3. Apply index calculation (``self.d[self.index]``)
        4. Apply statistical reducer (max/median/mean/percentile)
        5. Return composite image, or a placeholder band if the period has
           no images

        The method assumes valid inputs as validation occurs during
        initialization.
        
        Examples
        --------
        Generate winter composite for 2020::
        
            >>> processor = NdviSeasonality(periods=4, start_year=2020)
            >>> winter_2020 = processor.get_period_composite(2020, 0)  # 0=winter
        
        Generate July composite for monthly analysis::
        
            >>> monthly = NdviSeasonality(periods=12, sat='S2')
            >>> july_2021 = monthly.get_period_composite(2021, 6)  # 6=July (0-based)
        
        See Also
        --------
        get_year_composite : Calls this method for all periods
        _generate_periods : Defines period date ranges
        """
        # Extract temporal boundaries for the specified period
        init, ends = self._period_date_range(year, period_idx)

        # Pre-compute all statistical composites for efficiency
        # This approach avoids redundant filtering and index calculation
        period_stats = {}
        
        # Filter collection to specified temporal window and apply index calculation
        filtered_collection = self.ndvi_col.filterDate(init, ends).map(self.d[self.index])
        
        # Apply all statistical reducers to the filtered and processed collection
        period_stats['max'] = filtered_collection.max()
        period_stats['min'] = filtered_collection.min()
        period_stats['median'] = filtered_collection.median()
        period_stats['mean'] = filtered_collection.mean()
        period_stats['sum'] = filtered_collection.sum()
        period_stats['percentile'] = filtered_collection.reduce(ee.Reducer.percentile([self.percentile]))

        # Dispersion statistics. ee.ImageCollection.reduce() appends the reducer
        # name to every band ('nd' -> 'nd_stdDev'), so we strip that suffix to
        # keep the same band naming as max/median/mean and reuse the band
        # bookkeeping of get_year_composite
        period_stats['std'] = self._strip_reducer_suffix(
            filtered_collection.reduce(ee.Reducer.stdDev()), 'stdDev')
        period_stats['variance'] = self._strip_reducer_suffix(
            filtered_collection.reduce(ee.Reducer.variance()), 'variance')
        # Band names of an arithmetic operation come from the first operand,
        # which already carries the original names
        period_stats['range'] = period_stats['max'].subtract(period_stats['min'])
        period_stats['cv'] = period_stats['std'].divide(period_stats['mean'])

        # Valid observations per pixel. ee.ImageCollection has no count()
        # shortcut, so reduce() is used and its suffix stripped like above
        period_stats['count'] = self._strip_reducer_suffix(
            filtered_collection.reduce(ee.Reducer.count()), 'count')

        # A period without a single image (clouds filtered out at scene level,
        # sensor not yet in orbit, acquisition gaps...) reduces to an image with
        # no bands. Return a placeholder band instead so every period always
        # yields exactly one band: masked for the value reducers, and 0 for
        # 'count', since zero valid observations is a real value there.
        # Both branches are cast to the same type because Export refuses
        # images whose bands have different data types
        band_name = self._composite_band_name()
        if self.key == 'count':
            composite = period_stats['count'].toInt32()
            empty = ee.Image.constant(0).toInt32().rename(band_name)
        else:
            composite = period_stats[self.key].toFloat()
            empty = (ee.Image.constant(0).toFloat().rename(band_name)
                     .updateMask(ee.Image.constant(0)))

        # Return the composite corresponding to the user-specified statistical method
        return ee.Image(ee.Algorithms.If(
            filtered_collection.size().gt(0), composite, empty))

    def _composite_band_name(self):
        """
        Band name of the single band returned by :meth:`get_period_composite`.

        Index functions name their output ``'nd'`` on the optical and climate
        datasets and after the index on Sentinel-1 (``'VH'``, ``'RVI'``...).
        The percentile reducer appends ``'_p<percentile>'`` to it; the other
        reducers keep the name unchanged.

        Returns
        -------
        str
            Band name, e.g. ``'nd'``, ``'nd_p90'`` or ``'VH_p90'``.
        """
        if self.sat == 'S1':
            sar_names = {
                'vv': 'VV', 'vh': 'VH', 'rvi': 'RVI', 'vv_vh_ratio': 'RATIO',
                'dpsvi': 'DPSVI', 'rfdi': 'RFDI', 'vsdi': 'VSDI',
            }
            base = sar_names.get(self.index, 'VH')
        else:
            base = 'nd'

        if self.key == 'percentile':
            return f'{base}_p{self.percentile}'
        return base

    @staticmethod
    def _strip_reducer_suffix(image, reducer_name):
        """
        Remove the reducer suffix that ee.ImageCollection.reduce() appends.

        Earth Engine renames every band of a reduced collection as
        ``<band>_<reducer>`` (for example ``nd`` becomes ``nd_stdDev``). This
        helper restores the original band names so composites built with the
        dispersion reducers are named exactly like the max/median/mean ones.

        Parameters
        ----------
        image : ee.Image
            Output of ``ee.ImageCollection.reduce()``.
        reducer_name : str
            Suffix added by the reducer, without the leading underscore
            (e.g. ``'stdDev'``, ``'variance'``).

        Returns
        -------
        ee.Image
            Same image with the suffix removed from every band name.
        """
        pattern = '_' + reducer_name + '$'
        renamed = image.bandNames().map(
            lambda band: ee.String(band).replace(pattern, '')
        )
        return image.rename(renamed)
    
    def get_available_indices(self, satellite=None):
        """
        Get list of available spectral/radar indices for a specific satellite.
        
        Parameters
        ----------
        satellite : str or None, optional
            Satellite sensor: 'S2', 'S1', 'Landsat', 'MODIS', 'S3'.
            If None, uses current satellite (self.sat).
            
        Returns
        -------
        list of str
            Sorted list of available index names.
            
        Examples
        --------
        >>> processor = NdviSeasonality(sat='S2')
        >>> indices = processor.get_available_indices()
        >>> print(len(indices))  # S2 has most indices 29
        
        >>> # Check for different satellite
        >>> landsat_indices = processor.get_available_indices('Landsat')
        
        See Also
        --------
        get_all_available_indices : Get indices for all satellites
        """
        # Use current satellite if none specified
        if satellite is None:
            satellite = self.sat
            
        # Validate satellite parameter
        if satellite not in self.sensor_indices:
            available_sats = list(self.sensor_indices.keys())
            raise ValueError(f"Satellite '{satellite}' is not supported. Available satellites are: {available_sats}")
        
        # Return sorted list of available indices for the specified satellite
        return sorted(list(self.sensor_indices[satellite]))


    def get_all_available_indices(self):
        """
        Get all available indices organized by satellite sensor.
        
        Returns a comprehensive dictionary mapping each supported satellite to its
        available spectral or radar indices. This method provides a complete overview
        of the library's capabilities across all supported sensors.
        
        Returns
        -------
        dict
            Dictionary with satellite names as keys and sorted lists of available
            indices as values. Structure::

                {'satellite_name': ['index1', 'index2', ...], ...}
            
        Examples
        --------
        >>> processor = NdviSeasonality()
        >>> all_indices = processor.get_all_available_indices()
        >>> print(all_indices.keys())
        dict_keys(['S2', 'Landsat', 'MODIS', 'S1', 'S3'])
        
        >>> # Check capabilities of each sensor
        >>> for sensor, indices in all_indices.items():
        ...     print(f"{sensor}: {len(indices)} indices")
        S2: 29 indices      # Most comprehensive (basic + Red Edge)
        Landsat: 20 indices # Basic optical only
        MODIS: 20 indices   # Basic optical only  
        S1: 7 indices       # SAR only
        S3: 30 indices      # Basic optical + ocean/coastal
        
        >>> # Find common indices across optical sensors
        >>> optical_sensors = ['S2', 'Landsat', 'MODIS']
        >>> common_indices = set(all_indices[optical_sensors[0]])
        >>> for sensor in optical_sensors[1:]:
        ...     common_indices &= set(all_indices[sensor])
        >>> print(f"Common optical indices: {len(common_indices)}")
        Common optical indices: 20
        
        >>> # Check sensor-specific capabilities
        >>> s2_only = set(all_indices['S2']) - set(all_indices['Landsat'])
        >>> print(f"S2-only indices: {sorted(s2_only)}")
        S2-only indices: ['cire', 'ireci', 'mcari', 'mtci', 'ndci', 'ndre', 'psri', 'reip', 's2rep']
        
        >>> # Validate index availability before processing
        >>> target_index = 'ndre'
        >>> compatible_sensors = [sensor for sensor, indices in all_indices.items() 
        ...                      if target_index in indices]
        >>> print(f"'{target_index}' available on: {compatible_sensors}")
        'ndre' available on: ['S2']
        
        Notes
        -----
        **Sensor Capabilities Summary:**
        
        - **Sentinel-2 (S2)**: Most versatile optical sensor with Red Edge bands enabling advanced vegetation analysis and chlorophyll estimation
        - **Landsat**: Long-term optical observations (1982-present) with consistent band configuration across missions, ideal for time series analysis
        - **MODIS**: Global daily coverage at coarser resolution, excellent for large-scale monitoring and climate studies
        - **Sentinel-1 (S1)**: All-weather SAR observations, unique for detecting structural changes, crop monitoring, and flood mapping
        - **Sentinel-3 (S3)**: Specialized ocean and coastal monitoring with many spectral bands optimized for water quality assessment

        This method is particularly useful for:

        - Sensor capability comparison
        - Multi-sensor analysis planning
        - Index availability validation
        - Documentation and tutorial purposes
            
        See Also
        --------
        get_available_indices : Get indices for a specific satellite
        __init__ : Where sensor-index mappings are defined
        """
        # Create result dictionary with sorted indices for each sensor
        result = {}
        for sensor, indices in self.sensor_indices.items():
            result[sensor] = sorted(list(indices))
        return result
    
    def get_year_composite(self, return_counts: bool = False, count_valid_pixels: bool = False,
        scale_for_valid: int = 10, maxPixels_for_valid: float = 1e9,
        count_mode: str = "granules",  # 'granules' | 'unique_dates'
        return_df: bool = False, df_pivot: bool = False):
        """
        Generate temporal composite images for all years in the time range.
            
        Main processing method that creates multi-band images where each band
        represents a temporal period (season, month, etc.) using the configured
        statistical reducer and spectral/radar index.

        Parameters
        ----------
        return_counts : bool, optional
            If True, also return a list of dictionaries with image counts
            per year and period. Default is False (only the ImageCollection).
        count_valid_pixels : bool, optional
            If True, counts only images that contain at least one valid (non-masked)
            pixel within the ROI for each period. If False, simply counts the number
            of images in the filtered ImageCollection, regardless of whether they
            contribute valid pixels to the ROI. Default is False.
        scale_for_valid : int, optional
            Spatial resolution (meters) used when checking valid pixels with
            ``reduceRegion`` (only applies if ``count_valid_pixels=True``).
            Default is 10.
        maxPixels_for_valid : float, optional
            Maximum number of pixels allowed for the validity check
            (only applies if ``count_valid_pixels=True``). Default is 1e9.
        count_mode : {'granules', 'unique_dates'}, optional
            How to count inputs per period when ``return_counts=True``:
            - 'granules': count all scenes (granules) intersecting the ROI
            after filters (dates, clouds, etc.).
            - 'unique_dates': count unique acquisition dates (YYYY-MM-dd),
            collapsing multiple tiles/orbits from the same day into one.
            Default is 'granules'.

        Returns
        -------
        ee.ImageCollection
            Collection of multi-band composite images, one per year.
            Each image has exactly one band per temporal period, named after
            it, so every year shares the same bands. A period without images
            is kept as a fully masked band (all zeros for ``key='count'``).
        (ee.ImageCollection, list of dict), optional
            If ``return_counts=True``, returns a tuple containing the
            ImageCollection and a list of dictionaries. Each dictionary
            has the following keys:

                * year : int
                * period_idx : int
                * period_name : str
                * images_count : int
                * cloud_filter : bool
                * sat : str
                * index : str
                * key : str
                * percentile : int or None
                * count_mode : str

        Notes
        -----
        Processing workflow:

        1. Clear previous results (``self.imagelist = []``)
        2. For each year in range:

            a. Count the scenes of every period (one server call per year),
               store them in ``self.period_scene_counts`` as
               ``{year: [count per period]}``, and skip the year if none of
               its periods has any
            b. Optionally count images per period (granules or unique dates)
            c. Process all periods using :meth:`get_period_composite`, which
               returns a placeholder band for the periods without images
            d. Rename each period band to its period name
            e. Combine into multi-band image

        3. Return ImageCollection from processed images
        4. Optionally return image counts if ``return_counts=True``

        Each band is renamed to its period name before the periods are
        combined, so a missing period can never shift the names of the
        following ones. Band names use period names:

            * 4 periods: ``['winter', 'spring', 'summer', 'autumn']``
            * 12 periods: ``['january', 'february', ..., 'december']``
            * Custom: ``['p1', 'p2', ..., 'pN']``

        Examples
        --------
        Generate seasonal composites::

            >>> collection = processor.get_year_composite()
            >>> print(collection.size().getInfo())  # Number of years

        Count input images per period (unique dates)::

            >>> collection, counts = processor.get_year_composite(
            ...     return_counts=True,
            ...     count_valid_pixels=False,
            ...     count_mode="unique_dates"
            ... )
            >>> print(counts[0])
            {'year': 2020, 'period_idx': 0, 'period_name': 'january', 'images_count': 6, ...}

        Raises
        ------
        ee.EEException
            If Earth Engine computation fails.
        RuntimeError
            If no valid data found for any year.

        Warnings
        --------
        Years without a single image in any period are skipped with a console
        warning. Periods without images are reported and filled, not dropped.
        Large time ranges may approach computation limits.

        See Also
        --------
        get_period_composite : Generates individual period composites
        get_export : Export composites to files
        get_gif : Create animated visualization
        """

        # --- helpers internos ---
        def _filtered_collection_for_period(year, period_idx):
            init, ends = self._period_date_range(year, period_idx)
            # misma colección que se compone en get_period_composite
            return self.ndvi_col.filterDate(init, ends).map(self.d[self.index])

        def _count_images_for_period(year, period_idx):
            ic = _filtered_collection_for_period(year, period_idx)

            if count_mode == "granules":
                # 1) contar escenas (granulitos)
                if not count_valid_pixels:
                    return ic.size()  # ee.Number
                # versión "solo válidas" (más lenta)
                def flag_has_valid(img):
                    valid = img.reduceRegion(
                        reducer=ee.Reducer.count(),
                        geometry=self.roi,
                        scale=scale_for_valid,
                        maxPixels=maxPixels_for_valid
                    ).values().get(0)
                    return img.set('has_valid', ee.Number(valid).gt(0))
                used = ic.map(flag_has_valid).filter(ee.Filter.eq('has_valid', True))
                return used.size()

            elif count_mode == "unique_dates":
                # Usar la colección original (sin aplicar el índice) para no perder system:time_start
                init, ends = self._period_date_range(year, period_idx)
                ic_raw = self.ndvi_col.filterDate(init, ends)  # sin .map(self.d[self.index])

                ts = ic_raw.aggregate_array('system:time_start')  # ee.List de timestamps (ms)
                unique_days = ee.List(ts).map(
                    lambda t: ee.Date(t).format('YYYY-MM-dd')  # server-side
                ).distinct()
                return ee.Number(ee.List(unique_days).size())

            else:
                # fallback
                return ic.size()

        # limpiar resultados previos
        self.imagelist = []
        self.period_scene_counts = {}
        rows = []

        # recorrer años (end_year INCLUSIVO)
        for year in range(self.start_year, self.end_year + 1):
            # Escenas por periodo en una sola llamada al servidor. Solo se usa
            # para informar de los periodos vacíos y saltar los años sin datos:
            # get_period_composite ya devuelve una banda de relleno cuando un
            # periodo no tiene imágenes
            scene_counts = ee.List([
                self.ndvi_col.filterDate(*self._period_date_range(year, idx)).size()
                for idx in range(self.periods)
            ]).getInfo()
            self.period_scene_counts[year] = scene_counts

            if sum(scene_counts) == 0:
                print(f"Year {year}: No data available, skipping")
                continue

            period_images = []
            for period_idx in range(self.periods):
                period_name = self.period_names[period_idx]

                # contar SIEMPRE (aunque luego el composite no tenga datos)
                if return_counts:
                    n = _count_images_for_period(year, period_idx).getInfo()
                    rows.append({
                        'year': year,
                        'period_idx': period_idx,
                        'period_name': period_name,
                        'images_count': int(n),
                        'cloud_filter': bool(self.cloud_filter),
                        'sat': self.sat,
                        'index': self.index,
                        'key': self.key,
                        'percentile': self.percentile if self.key == 'percentile' else None,
                        'count_mode': count_mode
                    })

                if scene_counts[period_idx] == 0:
                    fill = "zeros" if self.key == 'count' else "a masked band"
                    print(f"No data for period {period_idx + 1} ({period_name}) in year {year}, filled with {fill}")

                # Renombrar cada periodo antes de unirlos, para que el nombre de
                # la banda no dependa de su posición en la imagen
                period_composite = self.get_period_composite(year, period_idx)
                period_images.append(period_composite.select([0]).rename(period_name))

            compositer = ee.Image.cat(period_images).clip(self.roi).set(
                'system:time_start', ee.Date(str(year)).millis()
            )
            self.imagelist.append(compositer)
            n_valid = sum(1 for c in scene_counts if c > 0)
            print(f"Year {year}: {n_valid}/{self.periods} periods with data using {self.index} index")

        collection = ee.ImageCollection.fromImages(self.imagelist)

        if not return_counts:
            return collection

        # return_counts=True a partir de aquí
        if not return_df:
            # Comportamiento anterior: lista de dicts
            return collection, rows

        # Intentar devolver un DataFrame
        try:
            import pandas as pd
        except Exception as e:
            print("Warning: pandas is not available; returning list of dicts instead.")
            return collection, rows

        df = pd.DataFrame(rows)

        if df_pivot:
            # Ancho: filas = year, columnas = period_name, valores = images_count
            if not df.empty:
                df = df.pivot_table(
                    index="year",
                    columns="period_name",
                    values="images_count",
                    aggfunc="first"
                ).sort_index(axis=1)
            else:
                # DataFrame vacío pero con la forma correcta
                df = pd.DataFrame()

        return collection, df

    def get_peak_period(self, per_year=False, across_years='median',
                        return_peak_value=False, mask_below=None,
                        return_ties=False, composite=None):
        """
        Map the period of the year in which each pixel reaches its maximum.

        Answers "*when*", not "*how much*": for every pixel, which of the
        temporal periods holds the highest value of the configured index. The
        result is a timing map — the month of peak greenness, the season of
        deepest flooding, the month a city shines brightest.

        The periods are the ones the instance is already configured with, so
        the output values run from 1 to ``self.periods`` and correspond to
        :attr:`period_names` in order.

        Parameters
        ----------
        per_year : bool, optional
            If False (default), the years are collapsed first (see
            ``across_years``) and a single image is returned: the period in
            which the *typical* year peaks. If True, returns one image per
            year, each tagged with a ``year`` property, so the peak can be
            tracked over time.
        across_years : {'median', 'mean', 'max'}, optional
            How the yearly composites are combined before the peak is looked
            up. Ignored when ``per_year=True``. Default is ``'median'``, which
            is robust to a single anomalous year.
        return_peak_value : bool, optional
            If True, the output keeps a second band, ``peak_value``, with the
            index value at the peak. Default is False.
        mask_below : float, optional
            Mask out the pixels whose peak value does not reach this
            threshold. The timing of a pixel with no signal is noise, so this
            is the usual way to restrict the map to pixels that are actually
            lit, flooded or vegetated. Default is None (keep every pixel).
        return_ties : bool, optional
            If True, the output keeps a band ``ties`` with the number of
            periods that reach that same maximum: 1 where the peak is unique,
            2 or more where the choice of period was arbitrary. Default is
            False.
        composite : ee.Image, optional
            Use this multi-band composite — one band per period, in order —
            instead of building one. For a composite that has already been
            processed: chlorophyll masked to the water pixels, an index
            masked to a land cover, a stack assembled by hand. Its bands are
            renamed to :attr:`period_names` if they differ, so the output
            still means "period number". Ignored when ``per_year=True``.

        Returns
        -------
        ee.Image or ee.ImageCollection
            Band ``peak_period``, an integer from 1 to ``self.periods``, plus
            ``peak_value`` and ``ties`` if requested. An ``ee.Image`` when
            ``per_year=False``, an ``ee.ImageCollection`` of one image per
            year otherwise. Pixels with no valid data in any period are
            masked.

        Notes
        -----
        The peak is resolved with ``qualityMosaic`` over one candidate image
        per period, which is what makes empty periods harmless: a period with
        no data is masked and simply does not compete. Two consequences worth
        knowing:

        * **Ties go to the later period.** If two periods hold exactly the
          same maximum, the second one wins. Rather than picking a winner by
          some other rule — the first, the middle, the average of the tied
          periods, all of which invent information the data does not have —
          ``return_ties=True`` reports how many periods tied, so the pixels
          where the answer is arbitrary can be found and dealt with
          downstream. Exact ties are common with integer values (a count, a
          saturated 6-bit DMSP city) and rare with reflectance.
        * A pixel masked in *every* period comes back masked, rather than
          defaulting to period 1.

        Like every composite of this library, the result carries no fixed
        projection, so anything that asks for it at another scale resamples
        it. That is harmless for a display, but it makes a
        ``frequencyHistogram`` of the periods report fractional counts for
        values no pixel holds: count the periods with ``eq()`` at the native
        scale, or ``reproject()`` the image first.

        This is not the phenological POS of
        :class:`~ndvi2gif.timeseries.SpatialPhenologyAnalyzer`, which fits and
        smooths the series before locating the maximum. Here the maximum is
        taken from the composites as they are: cruder, far cheaper, and
        defined for any index and any sensor.

        Examples
        --------
        Month of peak greenness over six years::

            >>> processor = NdviSeasonality(roi=roi, sat='S2', periods=12,
            ...                             index='ndvi', key='median',
            ...                             start_year=2018, end_year=2023)
            >>> peak = processor.get_peak_period()

        Month of peak nighttime light, restricted to lit pixels::

            >>> lights = NdviSeasonality(roi=city, sat='VIIRS', periods=12,
            ...                          index='avg_rad', key='mean',
            ...                          start_year=2017, end_year=2017)
            >>> peak = lights.get_peak_period(mask_below=2, return_peak_value=True)

        Follow the peak year by year::

            >>> peaks = processor.get_peak_period(per_year=True)
            >>> peaks.aggregate_array('year').getInfo()
        """
        if across_years not in ('median', 'mean', 'max'):
            raise ValueError(
                f"across_years must be 'median', 'mean' or 'max', got {across_years!r}"
            )

        if composite is not None and not per_year:
            composite = ee.Image(composite)
            n_bands = composite.bandNames().size().getInfo()
            if n_bands != self.periods:
                raise ValueError(
                    f"composite has {n_bands} bands but the instance is "
                    f"configured with {self.periods} periods: one band per "
                    "period is required, in order"
                )
            return self._peak_from_composite(
                composite.rename(self.period_names), return_peak_value,
                mask_below, return_ties)

        collection = self.get_year_composite()

        if per_year:
            # The years without a single scene are skipped by
            # get_year_composite, so the years are read back from the scene
            # counts instead of assuming the whole range is present
            years = [year for year, counts in self.period_scene_counts.items()
                     if sum(counts) > 0]
            peaks = [self._peak_from_composite(img, return_peak_value,
                                               mask_below, return_ties)
                     .set('year', year,
                          'system:time_start', img.get('system:time_start'))
                     for img, year in zip(self.imagelist, years)]
            return ee.ImageCollection.fromImages(peaks)

        combined = getattr(collection, across_years)()
        return self._peak_from_composite(combined, return_peak_value,
                                         mask_below, return_ties)

    def _peak_from_composite(self, composite, return_peak_value=False,
                             mask_below=None, return_ties=False):
        """
        Period of the maximum of a multi-band composite, one band per period.

        Builds one candidate image per period — the value and a constant with
        the period number, both carrying the period's own mask — and lets
        ``qualityMosaic`` pick the winner. See :meth:`get_peak_period`.
        """
        candidates = []
        for period_idx, period_name in enumerate(self.period_names):
            value = composite.select([period_name]).rename('peak_value').float()
            if mask_below is not None:
                value = value.updateMask(value.gte(mask_below))
            # The constant carries the period's mask, so an empty period
            # never wins by default
            period = (ee.Image.constant(period_idx + 1).toInt()
                      .rename('peak_period').updateMask(value.mask()))
            candidates.append(value.addBands(period))

        peak = ee.ImageCollection.fromImages(candidates).qualityMosaic('peak_value')

        bands = ['peak_period']
        if return_peak_value:
            bands.append('peak_value')

        if return_ties:
            # How many periods hold that same maximum. Counted over the same
            # candidates the mosaic chose from, not over the source composite:
            # they carry the same float precision, the same mask and the same
            # mask_below, so the winner always counts itself and 'ties' is
            # never 0
            tied = ee.ImageCollection.fromImages(
                [c.select('peak_value').eq(peak.select('peak_value'))
                 for c in candidates]
            ).reduce(ee.Reducer.sum()).toInt().rename('ties')
            # A maximum is reached by at least the period that holds it, so 1
            # is the floor. The mosaic and these comparisons are resampled
            # independently (a composite carries no fixed projection), which
            # otherwise leaves a handful of pixels counting zero
            tied = tied.unmask(0).max(1)
            peak = peak.addBands(
                tied.updateMask(peak.select('peak_value').mask()))
            bands.append('ties')

        return peak.select(bands).clip(self.roi)

    def get_peak_statistics(self, mask_below=None, min_years=2, peaks=None):
        """
        Circular mean and concentration of the peak period across the years.

        :meth:`get_peak_period` with ``per_year=True`` gives one peak period
        per year. Summarising those with an ordinary average is wrong, because
        **the period of the year is a circular variable**: the arithmetic mean
        of December (12) and January (1) is 6.5 — June — when the answer is
        late December. The wrap-around is the obvious failure, but not the only
        one. A pixel that peaks in February in half the years and in June in
        the other half averages to April, a month in which nothing ever
        happened there.

        The statistics that are defined on a circle are computed instead. Each
        period becomes an angle, the unit vectors are averaged, and the mean
        direction is read back as a period number:

        Each period ``p`` becomes an angle ``theta = (p - 1) * 2*pi / n``; the
        unit vectors are averaged into ``C = mean(cos theta)`` and
        ``S = mean(sin theta)``. The **circular mean** is ``atan2(S, C)``, read
        back as a period number, and the **resultant length**
        ``R = sqrt(C**2 + S**2)`` is the concentration: how much the years
        agree.

        ``R`` is the number that makes the mean readable, and it should be
        reported with it. R near 1 means every year peaked at the same time and
        the mean is a real date. R near 0 means the years point in every
        direction, the mean direction is arbitrary, and quoting it alone would
        invent a season. A bimodal pixel — February and June — has a low R,
        which is how it can be told from a genuinely April-peaking pixel that
        the arithmetic mean would have made it look like.

        Parameters
        ----------
        mask_below : float, optional
            Passed to :meth:`get_peak_period`: drop the pixels whose peak value
            does not reach this threshold, in the years where it does not. The
            timing of a pixel with no signal is noise. Ignored when ``peaks``
            is given.
        min_years : int, optional
            Mask the pixels with fewer than this many valid years. Two years
            give a direction but no useful concentration, so raise it when R is
            going to be interpreted. Default 2.
        peaks : ee.ImageCollection, optional
            A collection already produced by ``get_peak_period(per_year=True)``,
            to avoid recomputing it. Only its ``peak_period`` band is read.

        Returns
        -------
        ee.Image
            Band ``circular_mean``, the mean peak period, fractional and
            running from 1 to ``self.periods`` (a value of 12.6 with
            ``periods=12`` means "late December, a fifth of the way into
            January"); band ``concentration``, R between 0 and 1; and band
            ``n_years``, how many years the pixel contributed.

        Notes
        -----
        R has no threshold that makes it "significant"; it is a descriptive
        statistic. As a rough guide over yearly data, above about 0.8 the years
        cluster tightly, below about 0.4 the mean carries little meaning. If a
        test is wanted, R is the statistic the Rayleigh test of uniformity is
        built on.

        Examples
        --------
        Month of peak greenness and how repeatable it is::

            >>> processor = NdviSeasonality(roi=roi, periods=12, sat='S2',
            ...                             start_year=2018, end_year=2024)
            >>> stats = processor.get_peak_statistics()
            >>> reliable = stats.updateMask(stats.select('concentration').gt(0.7))

        Reusing a collection of peaks that is already in hand::

            >>> peaks = processor.get_peak_period(per_year=True, mask_below=0.3)
            >>> stats = processor.get_peak_statistics(peaks=peaks)
        """
        if peaks is None:
            peaks = self.get_peak_period(per_year=True, mask_below=mask_below)

        periods = ee.ImageCollection(peaks).select('peak_period')
        step = 2 * math.pi / self.periods

        # period 1 -> angle 0, so the angles tile the circle exactly
        angles = periods.map(lambda image: image.subtract(1).multiply(step))
        cos_mean = angles.map(lambda image: image.cos()).mean()
        sin_mean = angles.map(lambda image: image.sin()).mean()

        concentration = cos_mean.hypot(sin_mean).rename('concentration')

        # Earth Engine's atan2 takes the receiver as x and the argument as y,
        # the opposite order to math.atan2, so this reads cos.atan2(sin) to
        # mean atan2(S, C). Getting it backwards fails silently: the
        # concentration stays right and only the mean direction is wrong.
        # Then, atan2 returns (-pi, pi]; shifting before the modulo keeps the
        # result positive, so it maps back onto period numbers 1..periods
        mean_angle = cos_mean.atan2(sin_mean).add(2 * math.pi).mod(2 * math.pi)
        circular_mean = mean_angle.divide(step).add(1).rename('circular_mean')

        n_years = periods.count().rename('n_years')

        stats = circular_mean.addBands(concentration).addBands(n_years)
        return stats.updateMask(n_years.gte(min_years)).clip(self.roi)

    def get_water_mask(self, mode='dynamic', water_index='mndwi', threshold=0.0,
                       occurrence=90, return_area=False):
        """
        Water masks for the configured period grid, one band per period.

        Which pixels count as water is a choice, not a fact, and it changes
        the answer: over Iznajar in a single July, NDWI called 11.1 km2 water,
        MNDWI 10.4 and AWEInsh 7.0. Leaving that choice to be improvised in
        each analysis is how two studies of the same reservoir stop being
        comparable, so it is a parameter here.

        Parameters
        ----------
        mode : str, tuple or ee.Geometry, optional
            How the water surface is defined.

            * ``'dynamic'`` (default) — the water detected in *each period of
              each year*. Follows a drawdown, and is the only mode that can.
            * ``'permanent'`` — the pixels that are water in **every period
              that has data**: the minimum extent, the pool that never dries.
              Periods with no scene are skipped rather than counted as dry,
              so one cloudy month does not empty the pool.
            * ``'maximum'`` — the pixels that are water in **any** period: the
              maximum extent.
            * ``('jrc', occurrence)`` or ``'jrc'`` — JRC Global Surface Water
              occurrence at or above ``occurrence`` percent.
            * an ``ee.Geometry``, ``ee.FeatureCollection``, or a path to a
              shapefile or GeoJSON — a zone you supply.
        water_index : str, optional
            Index used to detect water in the first three modes. One of
            ``'ndwi'``, ``'mndwi'``, ``'awei'``, ``'aweinsh'``, ``'wi2015'``.
            Default ``'mndwi'``, which handles built-up surroundings better
            than NDWI.
        threshold : float, optional
            A pixel is water where the index exceeds this. Default 0.0.
        occurrence : int, optional
            Occurrence percentage for the ``'jrc'`` mode. Default 90.
        return_area : bool, optional
            If True, also return the masked area in km2 per year and period.

        Returns
        -------
        ee.ImageCollection or tuple
            One image per year, each with one band per period named after
            :attr:`period_names`, holding 1 on water and masked elsewhere. The
            static modes repeat the same mask in every band and every year, so
            that consumers do not have to special-case them. With
            ``return_area=True``, a ``(collection, areas)`` tuple where
            ``areas`` is ``{year: [km2, ...]}``.

        Notes
        -----
        **Report the area whenever the mask is dynamic.** Each period is then
        averaged over a different population of pixels, so a change in the
        series mixes a change in the quantity with a change in what was
        measured. The area is what separates the two, which is why
        ``return_area`` exists and why :attr:`water_mask_area` is filled in.

        **A dense algal bloom can fall out of a dynamic mask, and that bias
        runs downwards.** Bloom water is spectrally close to vegetation — the
        near infrared rises — so a water index stops calling it water, and the
        pixels dropped are exactly the ones carrying the most chlorophyll. A
        mean taken over what is left is biased low, and nothing in the output
        says so. MNDWI holds up better than NDWI here, because the water under
        a bloom still absorbs the shortwave infrared almost completely, which
        is why it is the default. Check the area series against the index
        series before trusting either: a drop in masked area at the same time
        as a peak in chlorophyll is the signature of this, not of a drawdown.
        The opposite mistake is just as real — a fixed polygon keeps the dry
        tail, whose index in summer is pasture, not algae.

        **``'maximum'`` is the one that misleads.** It includes pixels that
        are dry most of the time, and the chlorophyll index of a dry pixel is
        grassland. It defines a study domain; it does not compute a mean over
        water.

        **The JRC layer ends in December 2021** (GSW 1.4 covers 1984-03 to
        2021-12). It is stable and comparable across sites, which is its
        virtue, but it knows nothing about a drought that came later.

        Examples
        --------
        Chlorophyll over the water actually present in each month::

            >>> ndci = NdviSeasonality(roi=reservoir, sat='S2', index='ndci',
            ...                        periods=12, start_year=2020, end_year=2024)
            >>> masks, area = ndci.get_water_mask(return_area=True)
            >>> composites = ndci.get_year_composite()

        The pool that never dries, as a single image::

            >>> permanent = ee.Image(ndci.get_water_mask('permanent').first())

        A fixed domain from the global layer::

            >>> domain = ndci.get_water_mask(('jrc', 95))
        """
        from .hydroperiod import HydroperiodAnalyzer

        occurrence_pct = occurrence
        if isinstance(mode, (tuple, list)):
            if len(mode) != 2 or mode[0] != 'jrc':
                raise ValueError(
                    "a tuple mode must be ('jrc', occurrence), got "
                    f"{mode!r}"
                )
            mode, occurrence_pct = 'jrc', mode[1]

        years = list(range(self.start_year, self.end_year + 1))

        # --- a zone the user supplies: no detection at all
        if not isinstance(mode, str):
            zone = mode
            if isinstance(zone, str):
                zone = (geemap.shp_to_ee(zone) if zone.endswith('.shp')
                        else geemap.geojson_to_ee(zone))
            if isinstance(zone, ee.FeatureCollection):
                zone = zone.geometry()
            static = ee.Image.constant(1).clip(zone).rename('water')
            return self._broadcast_water_mask(static, years, return_area)

        if mode == 'jrc':
            gsw = ee.Image('JRC/GSW1_4/GlobalSurfaceWater').select('occurrence')
            static = gsw.gte(occurrence_pct).selfMask().rename('water')
            return self._broadcast_water_mask(static, years, return_area)

        if mode not in ('dynamic', 'permanent', 'maximum'):
            raise ValueError(
                f"mode must be 'dynamic', 'permanent', 'maximum', 'jrc', a "
                f"('jrc', occurrence) tuple or a geometry, got {mode!r}"
            )

        if water_index not in HydroperiodAnalyzer.WATER_INDICES:
            raise ValueError(
                f"water_index must be one of "
                f"{sorted(HydroperiodAnalyzer.WATER_INDICES)}, got "
                f"{water_index!r}"
            )

        # Composite the water index without disturbing the instance: the
        # caller is configured for some other index, and get_year_composite
        # resets imagelist and period_scene_counts as it goes
        # these two only exist once get_year_composite has run at least once
        saved = (self.index, getattr(self, 'imagelist', None),
                 getattr(self, 'period_scene_counts', None))
        try:
            self.index = water_index
            water = self.get_year_composite()
        finally:
            self.index = saved[0]
            if saved[1] is None:
                self.imagelist, self.period_scene_counts = [], {}
            else:
                self.imagelist, self.period_scene_counts = saved[1], saved[2]

        detected = water.map(lambda image: image.gt(threshold))

        def counts(image):
            """Periods that are water, and periods that have data at all."""
            return (image.reduce(ee.Reducer.sum()).rename('wet')
                    .addBands(image.reduce(ee.Reducer.count()).rename('valid')))

        if mode == 'permanent':
            # Water in every period *that has data*, not in every period. A
            # period with no scene comes back as a fully masked band, so
            # comparing against self.periods would make the permanent pool
            # empty for any year with a gap -- measured on a Guadalquivir
            # reservoir with no March or April scene in 2022, it returned
            # 0.00 km2 instead of the pool that is there all year
            wet = detected.map(
                lambda image: counts(image).expression(
                    'valid > 0 && wet == valid',
                    {'wet': counts(image).select('wet'),
                     'valid': counts(image).select('valid')}))
        elif mode == 'maximum':
            wet = detected.map(
                lambda image: counts(image).select('wet').gt(0))

        if mode == 'dynamic':
            masks = detected.map(lambda image: image.selfMask())
        else:
            # one band per period so every mode has the same shape downstream
            masks = wet.map(lambda image: ee.Image.cat(
                [image.selfMask()] * self.periods).rename(self.period_names))

        masks = ee.ImageCollection(masks)
        self.water_mask_area = self._water_mask_area(masks, years)
        return (masks, self.water_mask_area) if return_area else masks

    def _broadcast_water_mask(self, static, years, return_area):
        """Repeat one static mask over every period and every year."""
        image = ee.Image.cat([static] * self.periods).rename(self.period_names)
        masks = ee.ImageCollection(
            [image.set('year', year) for year in years])
        self.water_mask_area = self._water_mask_area(masks, years)
        return (masks, self.water_mask_area) if return_area else masks

    def _water_mask_area(self, masks, years):
        """Masked area in km2, per year and period, in a single request."""
        pixel_area = ee.Image.pixelArea().divide(1e6)
        scale = self._default_scale_for_sat()

        def measure(image):
            return ee.Feature(None, image.multiply(pixel_area)
                              .rename(self.period_names)
                              .reduceRegion(ee.Reducer.sum(), self.roi, scale,
                                            maxPixels=1e9, bestEffort=True))

        rows = ee.FeatureCollection(masks.map(measure)).getInfo()['features']
        return {year: [row['properties'].get(name) for name in self.period_names]
                for year, row in zip(years, rows)}

    # DMSP-OLS runs to 2013 and the stray-light corrected VIIRS product the
    # library uses for `sat='VIIRS'` only starts in 2014, so the fit has to
    # reach for the uncorrected one, which begins in April 2012
    _VIIRS_OVERLAP_PRODUCT = 'NOAA/VIIRS/DNB/MONTHLY_V1/VCMCFG'
    _DMSP_SATURATION = 63

    def calibrate_nighttime_lights(self, n_pixels=5000, seed=7,
                                   compare_products=True):
        """
        Fit the conversion from VIIRS radiance to DMSP-like digital numbers.

        The two nighttime lights records do not join up. DMSP-OLS (1992-2013)
        is a 6-bit sensor with no on-board calibration; VIIRS (2012-present)
        reports radiance in nW/cm2/sr. Plotting one after the other produces a
        step at the sensor change that is pure artefact — measured on the
        Andalusian coast, "lit area" fell from 7888 km2 in 2013 to 2909 km2 in
        2014 simply because ``DN >= 20`` and ``>= 5 nW/cm2/sr`` are two
        different definitions of lit.

        The records do overlap, so the conversion can be fitted on this ROI
        rather than borrowed from a paper about somewhere else. **VIIRS is
        converted into DMSP-like values, not the other way round**: DMSP
        saturates at 63 over city centres, and a saturating 6-bit integer
        cannot be inverted into an unbounded radiance.

        Parameters
        ----------
        n_pixels : int, optional
            Pixels sampled for the fit. Default 5000.
        seed : int, optional
            Sampling seed, so a fit can be repeated. Default 7.
        compare_products : bool, optional
            Also measure how far the VIIRS product used for the fit is from
            the one ``sat='VIIRS'`` composites, over the years where both
            exist. Default True. See the notes.

        Returns
        -------
        dict
            ``coefficients`` (``a``, ``b``, the ``ceiling`` and the ``model``
            as a string), ``r2``, ``n`` (pixels that entered the fit),
            ``years`` used, and, with ``compare_products``,
            ``product_difference``: the median relative difference between the
            two VIIRS products over this ROI.

        Notes
        -----
        **The fit rests on two annual DMSP composites, 2012 and 2013**, not on
        twenty-one independent dates. The regression runs across pixels, which
        is how the literature does it, but the ``r2`` describes agreement
        between two images over space, not a time series.

        **The coefficients are fitted on a different VIIRS product from the one
        they will be applied to.** DMSP ends in 2013 and the stray-light
        corrected VCMSLCFG only starts in 2014, so the fit uses VCMCFG.
        The two differ mainly at high latitudes in summer;
        ``product_difference`` measures the gap over *your* ROI so the
        assumption can be judged rather than assumed.

        **The model is saturating, and deliberately so.** Measured over the
        Andalusian coast, the mean digital number rises 10.6, 26.1, 33.0,
        44.9, 51.8, 58.8, 61.2, 62.5 across radiance bins from under 1 to over
        80 nW/cm2/sr: monotone, flattening against the 63 ceiling, with 85% of
        the brightest bin saturated. A polynomial fitted to that turns over —
        a quadratic predicted DN 0 for the centre of Seville — so
        ``63 * (1 - exp(-a * x**b))`` is used instead, which cannot exceed the
        ceiling or fall as radiance rises. It linearises, so the fit is still
        a single linear regression.

        **Saturated pixels are excluded.** A DMSP pixel at 63 records only
        that the sensor ran out of range, and the linearisation has no
        logarithm there.

        Examples
        --------
        >>> lights = NdviSeasonality(roi=city, sat='VIIRS', index='avg_rad',
        ...                          periods=12, start_year=2014, end_year=2024)
        >>> fit = lights.calibrate_nighttime_lights()
        >>> fit['r2'], fit['coefficients']
        >>> dmsp_like = lights.to_dmsp_like(some_viirs_image, fit['coefficients'])
        """
        dmsp_col = ee.ImageCollection('NOAA/DMSP-OLS/NIGHTTIME_LIGHTS')
        viirs_col = ee.ImageCollection(self._VIIRS_OVERLAP_PRODUCT)

        years = [2012, 2013]
        pairs = []
        for year in years:
            dmsp = dmsp_col.filterDate(f'{year}-01-01', f'{year + 1}-01-01') \
                           .select('stable_lights').mean()
            # A collection reduction has no projection of its own, so it would
            # default to WGS84 at one degree and reduceResolution below would
            # silently work on that grid instead of the DMSP one
            grid = ee.Image(dmsp_col.first()).projection()
            dmsp = dmsp.setDefaultProjection(grid)

            viirs = viirs_col.filterDate(f'{year}-01-01', f'{year + 1}-01-01') \
                             .select('avg_rad').mean()
            viirs = viirs.setDefaultProjection(
                ee.Image(viirs_col.first()).projection())
            # VIIRS is finer than DMSP, so it is aggregated onto the DMSP grid
            viirs = viirs.reduceResolution(ee.Reducer.mean(), maxPixels=1024) \
                         .reproject(grid)

            valid = (dmsp.lt(self._DMSP_SATURATION)
                     .And(dmsp.gt(0))
                     .And(viirs.gt(0)))
            pairs.append(viirs.rename('x').addBands(dmsp.rename('y'))
                              .updateMask(valid))

        samples = ee.FeatureCollection([
            pair.sample(region=self.roi, scale=1000, numPixels=n_pixels // 2,
                        seed=seed + i, dropNulls=True)
            for i, pair in enumerate(pairs)]).flatten()

        # The relationship is monotone and saturates towards 63, which a
        # polynomial cannot represent without turning over: a quadratic fitted
        # here curved back down above ~50 nW and predicted DN 0 for the centre
        # of Seville. A saturating model cannot do that,
        #
        #     DN = 63 * (1 - exp(-a * x**b))
        #
        # and it linearises, so the fit is still one linear regression:
        #
        #     log(-log(1 - DN/63)) = log(a) + b * log(x)
        #
        # Saturated pixels have to go, since 1 - 63/63 = 0 has no logarithm --
        # and they carry no information anyway, only that the sensor ran out
        # of range.
        def design(feature):
            x = ee.Number(feature.get('x'))
            dn = ee.Number(feature.get('y'))
            response = ee.Number(1).subtract(dn.divide(self._DMSP_SATURATION)) \
                         .log().multiply(-1).log()
            return feature.set('log_x', x.log(), 'one', 1, 'response', response)

        prepared = samples.map(design)
        fitted = prepared.reduceColumns(
            ee.Reducer.linearRegression(numX=2, numY=1),
            ['one', 'log_x', 'response']).getInfo()

        intercept, b = (row[0] for row in fitted['coefficients'])
        a = math.exp(intercept)
        coefficients = {'a': a, 'b': b, 'ceiling': self._DMSP_SATURATION,
                        'model': '63 * (1 - exp(-a * x**b))'}

        # r2 against the digital numbers themselves, not the linearised form,
        # because that is the scale the answer is read on
        def predict(feature):
            x = ee.Number(feature.get('x'))
            value = ee.Number(1).subtract(
                x.pow(b).multiply(-a).exp()).multiply(self._DMSP_SATURATION)
            return feature.set('fit', value)

        stats = prepared.map(predict).reduceColumns(
            ee.Reducer.pearsonsCorrelation(), ['fit', 'y']).getInfo()
        r2 = stats['correlation'] ** 2

        result = {'coefficients': coefficients, 'r2': r2,
                  'n': samples.size().getInfo(), 'years': years}

        if compare_products:
            result['product_difference'] = self._viirs_product_difference()
        return result

    def _viirs_product_difference(self, year=2015):
        """Median relative gap between the two VIIRS monthly products."""
        uncorrected = ee.ImageCollection(self._VIIRS_OVERLAP_PRODUCT) \
            .filterDate(f'{year}-01-01', f'{year + 1}-01-01').select('avg_rad').mean()
        corrected = ee.ImageCollection('NOAA/VIIRS/DNB/MONTHLY_V1/VCMSLCFG') \
            .filterDate(f'{year}-01-01', f'{year + 1}-01-01').select('avg_rad').mean()

        lit = corrected.gt(1)
        relative = uncorrected.subtract(corrected).abs().divide(corrected) \
                              .updateMask(lit)
        value = relative.reduceRegion(
            ee.Reducer.median(), self.roi, 500,
            maxPixels=1e9, bestEffort=True).values().get(0)
        return ee.Number(value).getInfo()

    def to_dmsp_like(self, image, coefficients, band='avg_rad'):
        """
        Apply a fitted conversion, turning VIIRS radiance into DMSP-like DN.

        See :meth:`calibrate_nighttime_lights`. The output is clamped to the
        DMSP range, 0 to 63, because values outside it do not exist on that
        scale — though a pixel pushed to 63 is at the ceiling for the same
        reason the excluded pixels were, and means only "brighter than DMSP
        could record".

        Parameters
        ----------
        image : ee.Image
            A VIIRS image or composite.
        coefficients : dict
            The ``coefficients`` entry of the fit.
        band : str, optional
            Radiance band to convert. Default ``'avg_rad'``.

        Returns
        -------
        ee.Image
            One band, ``dmsp_like``.
        """
        radiance = ee.Image(image).select(band).max(0)
        a = coefficients['a']
        b = coefficients['b']
        ceiling = coefficients.get('ceiling', self._DMSP_SATURATION)
        converted = ee.Image.constant(1).subtract(
            radiance.pow(b).multiply(-a).exp()).multiply(ceiling)
        return converted.clamp(0, self._DMSP_SATURATION).rename('dmsp_like')

    # Index calculation methods (same as original - keeping all of them)
    def get_raw_band(self, image, band):
        """
        Return a single reflectance band through the regular index machinery.

        Spectral indices are ratios, so they cancel out multiplicative changes
        in brightness: a pixel can keep its NDVI unchanged while its actual
        reflectance drifts. Radiometric work needs the bands themselves. The
        typical case is selecting pseudo-invariant features (PIFs) for relative
        radiometric normalization, where the criterion is per-band reflectance
        stability, not index stability.

        The value is surface reflectance (0-1) on every sensor, so the same
        threshold means the same thing whichever produced it: Landsat is
        rescaled by :func:`scale_OLI` / :func:`scale_ETM`, and Sentinel-2 and
        MODIS, which store reflectance as integers multiplied by 10000, by
        :func:`_to_reflectance` when the collection is built.

        Available as ``index='blue'``, ``'green'``, ``'red'``, ``'nir'``,
        ``'swir1'``, ``'swir2'`` on Sentinel-2, Landsat and MODIS, plus
        ``'red_edge1'``, ``'red_edge2'`` and ``'red_edge3'`` on Sentinel-2.
        Sentinel-3 is not supported: its bands are TOA radiances.

        Parameters
        ----------
        image : ee.Image
            Image with standardized band names (Blue, Green, Red, Nir, Swir1,
            Swir2, Red_Edge1-3).
        band : str
            Standardized band name to extract.

        Returns
        -------
        ee.Image
            Single-band image of surface reflectance, renamed to ``'nd'`` to
            match the band bookkeeping of :meth:`get_year_composite`.

        Examples
        --------
        Map how stable SWIR1 is over a whole Sentinel-2 series, one composite
        per year covering the full year::

            >>> n = NdviSeasonality(roi=roi, start_year=2018, end_year=2025,
            ...                     periods=1, sat='S2', index='swir1',
            ...                     key='std')
            >>> within_year = n.get_year_composite()

        Combine it with the between-year dispersion of the annual means::

            >>> m = NdviSeasonality(roi=roi, start_year=2018, end_year=2025,
            ...                     periods=1, sat='S2', index='swir1',
            ...                     key='mean')
            >>> between_year = m.get_year_composite().reduce(ee.Reducer.stdDev())

        See Also
        --------
        get_available_indices : Lists the indices and bands available per sensor
        """
        return image.select(band).multiply(self.reflectance_scale).rename(['nd'])

    def get_ndvi(self, image):
        """
        Normalized Difference Vegetation Index - Most widely used vegetation index.
        
        References
        ----------
        Rouse, J.W., Haas, R.H., Schell, J.A., Deering, D.W. (1974). 
        Monitoring vegetation systems in the Great Plains with ERTS. 
        Third ERTS Symposium, NASA SP-351 I: 309-317.
        """
        return image.normalizedDifference(['Nir', 'Red'])

    def get_ndwi(self, image):
        """
        Normalized Difference Water Index - Water body detection and monitoring.
        
        References
        ----------
        Gao, B. (1996). NDWI - A normalized difference water index for remote sensing 
        of vegetation liquid water from space. Remote Sensing of Environment, 58(3), 257-266.
        """
        return image.normalizedDifference(['Green', 'Nir'])

    def get_mndwi(self, image):
        """
        Modified Normalized Difference Water Index - Enhanced water detection, reduces built-up noise.
        
        References
        ----------
        Xu, H. (2006). Modification of normalised difference water index (NDWI) to enhance 
        open water features in remotely sensed imagery. International Journal of Remote Sensing, 27(14), 3025-3033.
        """
        return image.normalizedDifference(['Green', 'Swir1'])

    def get_evi(self, image):
        """
        Enhanced Vegetation Index - Improved vegetation monitoring with reduced atmospheric influence.
        
        References
        ----------
        Huete, A., Didan, K., Miura, T., Rodriguez, E.P., Gao, X., Ferreira, L.G. (2002). 
        Overview of the radiometric and biophysical performance of the MODIS vegetation indices. 
        Remote Sensing of Environment, 83(1-2), 195-213.
        """
        return image.expression(
            '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))', {
            'NIR': image.select('Nir'),
            'RED': image.select('Red'),
            'BLUE': image.select('Blue')}).rename(['nd'])

    def get_savi(self, image, L=0.428):
        """
        Soil Adjusted Vegetation Index - Reduces soil background influence on vegetation indices.
        
        References
        ----------
        Huete, A.R. (1988). A soil-adjusted vegetation index (SAVI). 
        Remote Sensing of Environment, 25(3), 295-309.
        """
        return image.expression(
            '((NIR - RED) / (NIR + RED + L) * (1 + L))', {
            'NIR': image.select('Nir'),
            'RED': image.select('Red'),
            'L': L}).rename(['nd'])

    def get_aweinsh(self, image):
        """
        Automated Water Extraction Index (no shadow) - Water detection without shadow pixels.
        
        References
        ----------
        Feyisa, G.L., Meilby, H., Fensholt, R., Proud, S.R. (2014). 
        Automated Water Extraction Index: A new technique for surface water mapping using Landsat imagery. 
        Remote Sensing of Environment, 140, 23-35.
        """
        return image.expression(
            '4.0 * (GREEN - SWIR1) - (0.25 * NIR + 2.75 * SWIR2)', {
            'NIR': image.select('Nir'),
            'GREEN': image.select('Green'),
            'SWIR1':image.select('Swir1'),
            'SWIR2':image.select('Swir2')}).rename(['nd'])

    def get_awei(self, image):
        """
        Automated Water Extraction Index - Water detection with shadow consideration.
        
        References
        ----------
        Feyisa, G.L., Meilby, H., Fensholt, R., Proud, S.R. (2014). 
        Automated Water Extraction Index: A new technique for surface water mapping using Landsat imagery. 
        Remote Sensing of Environment, 140, 23-35.
        """
        return image.expression(
            ('BLUE + 2.5 * GREEN - 1.5 * (NIR + SWIR1) - 0.25 * SWIR2'), {
            'NIR': image.select('Nir'),
            'BLUE': image.select('Blue'),
            'GREEN': image.select('Green'),
            'SWIR1':image.select('Swir1'),
            'SWIR2':image.select('Swir2')}).rename(['nd'])

    def get_gndvi(self, image):
        """
        Green Normalized Difference Vegetation Index - Sensitive to chlorophyll content.
        
        References
        ----------
        Gitelson, A.A., Kaufman, Y.J., Merzlyak, M.N. (1996). 
        Use of a green channel in remote sensing of global vegetation from EOS-MODIS. 
        Remote Sensing of Environment, 58(3), 289-298.
        """
        return image.normalizedDifference(['Nir', 'Green'])
    
    def get_avi(self, image, L=0.428):
        """
        Advanced Vegetation Index - Non-linear vegetation index for dense vegetation.
        
        References
        ----------
        Bannari, A., Asalhi, H., Teillet, P.M. (2002). 
        Transformed difference vegetation index (TDVI) for vegetation cover mapping. 
        IEEE International Geoscience and Remote Sensing Symposium, 5, 3053-3055.
        """
        return image.expression(
            '(NIR * (1.0 - RED) * (NIR - RED)) ** (1/3)', {
            'NIR': image.select('Nir'),
            'RED': image.select('Red')}).rename(['nd'])

    def get_nbri(self, image):
        """
        Normalized Burn Ratio Index - Fire damage and burn severity assessment.
        
        References
        ----------
        Key, C., Benson, N. (2006). Landscape Assessment: Ground measure of severity, 
        the Composite Burn Index; and Remote sensing of severity, the Normalized Burn Ratio. 
        FIREMON: Fire Effects Monitoring and Inventory System, RMRS-GTR-164-CD.
        """
        return image.normalizedDifference(['Nir', 'Swir2'])

    def get_ndsi(self, image):
        """
        Normalized Difference Snow Index - Snow cover detection and monitoring.
        
        References
        ----------
        Dozier, J. (1989). Spectral signature of alpine snow cover from the Landsat 
        Thematic Mapper. Remote Sensing of Environment, 28, 9-22.
        """
        return image.normalizedDifference(['Green', 'Swir1'])

    def get_ndmi(self, image):
        """
        Normalized Difference Moisture Index - Vegetation water content assessment.
        
        References
        ----------
        Hardisky, M.A., Klemas, V., Smart, R.M. (1983). The influence of soil salinity, 
        growth form, and leaf moisture on the spectral radiance of Spartina alterniflora canopies. 
        Photogrammetric Engineering and Remote Sensing, 49(1), 77-83.
        """
        return image.normalizedDifference(['Nir', 'Swir1'])

    def get_msi(self, image):
        """
        Moisture Stress Index - Plant water stress detection.
        
        References
        ----------
        Rock, B.N., Vogelmann, J.E., Williams, D.L., Vogelmann, A.F., Hoshizaki, T. (1986). 
        Remote detection of forest damage. BioScience, 36(7), 439-445.
        """
        return image.expression(
            'SWIR1 / NIR', {
            'NIR': image.select('Nir'),
            'SWIR1': image.select('Swir1')
        }).rename(['nd'])

    def get_nmi(self, image):
        """
        Normalized Multi-band Drought Index - Multi-spectral drought assessment.
        
        References
        ----------
        Wang, L., Qu, J.J. (2007). NMDI: A normalized multi‐band drought index for monitoring 
        soil and vegetation moisture with satellite remote sensing. 
        Geophysical Research Letters, 34(20), L20405.
        """
        return image.expression(
            '(NIR - (SWIR1 - SWIR2)) / (NIR + (SWIR1 - SWIR2))', {
            'NIR': image.select('Nir'),
            'SWIR1': image.select('Swir1'),
            'SWIR2': image.select('Swir2')
        }).rename(['nd'])

    def get_ndti(self, image):
        """
        Normalized Difference Tillage Index - Agricultural tillage and residue detection.
        
        References
        ----------
        Van Deventer, A.P., Ward, A.D., Gowda, P.H., Lyon, J.G. (1997). 
        Using thematic mapper data to identify contrasting soil plains and tillage practices. 
        Photogrammetric Engineering and Remote Sensing, 63(1), 87-93.
        """
        return image.normalizedDifference(['Swir1', 'Swir2'])

    def get_cri1(self, image):
        """
        Carotenoid Reflectance Index 1 - Carotenoid pigment detection.
        
        References
        ----------
        Gitelson, A.A., Zur, Y., Chivkunova, O.B., Merzlyak, M.N. (2002). 
        Assessing carotenoid content in plant leaves with reflectance spectroscopy. 
        Photochemistry and Photobiology, 75(3), 272-281.
        """
        return image.expression(
            '(1 / BLUE) - (1 / GREEN)', {
            'BLUE': image.select('Blue'),
            'GREEN': image.select('Green')
        }).rename(['nd'])

    def get_cri2(self, image):
        """
        Carotenoid Reflectance Index 2 - Alternative carotenoid assessment.
        
        References
        ----------
        Gitelson, A.A., Zur, Y., Chivkunova, O.B., Merzlyak, M.N. (2002). 
        Assessing carotenoid content in plant leaves with reflectance spectroscopy. 
        Photochemistry and Photobiology, 75(3), 272-281.
        """
        return image.expression(
            '(1 / BLUE) - (1 / RED)', {
            'BLUE': image.select('Blue'),
            'RED': image.select('Red')
        }).rename(['nd'])

    def get_lai(self, image):
        """
        Leaf Area Index approximation - Estimate of leaf area per unit ground area.
        
        References
        ----------
        Boegh, E., Soegaard, H., Broge, N., Hasager, C.B., Jensen, N.O., Schelde, K., Thomsen, A. (2002). 
        Airborne multispectral data for quantifying leaf area index, nitrogen concentration, 
        and photosynthetic efficiency in agriculture. Remote Sensing of Environment, 81(2-3), 179-193.
        """
        return image.expression(
            '3.618 * EVI - 0.118', {
            'EVI': self.get_evi(image).select('nd')
        }).rename(['nd'])

    def get_pri(self, image):
        """
        Photochemical Reflectance Index - Plant stress and photosynthetic efficiency.
        
        References
        ----------
        Gamon, J.A., Peñuelas, J., Field, C.B. (1992). 
        A narrow-waveband spectral index that tracks diurnal changes in photosynthetic efficiency. 
        Remote Sensing of Environment, 41(1), 35-44.
        """
        return image.normalizedDifference(['Green', 'Blue'])

    def get_wdrvi(self, image):
        """
        Wide Dynamic Range Vegetation Index - Enhanced vegetation monitoring for dense canopies.
        
        References
        ----------
        Gitelson, A.A. (2004). Wide dynamic range vegetation index for remote quantification 
        of biophysical characteristics of vegetation. Journal of Plant Physiology, 161(2), 165-173.
        """
        return image.expression(
            '(0.1 * NIR - RED) / (0.1 * NIR + RED)', {
            'NIR': image.select('Nir'),
            'RED': image.select('Red')
        }).rename(['nd'])

    def get_cig(self, image):
        """
        Chlorophyll Index Green - Chlorophyll content estimation using green band.
        
        References
        ----------
        Gitelson, A.A., Gritz, Y., Merzlyak, M.N. (2003). 
        Relationships between leaf chlorophyll content and spectral reflectance and algorithms 
        for non-destructive chlorophyll assessment in higher plant leaves. 
        Journal of Plant Physiology, 160(3), 271-282.
        """
        return image.expression(
            '(Nir / Green) - 1', {
            'Green': image.select('Green'),
            'Nir': image.select('Nir')
        }).rename(['nd'])
    
    def get_lst(self, image):
        """
        Land Surface Temperature (LST) - Multi-sensor implementation with robust error handling
        Automatically detects sensor and applies appropriate scaling and conversion
        
        Supported sensors and bands:
        - Landsat 8/9 (OLI/TIRS): ST_B10 (Band 10 - thermal infrared)
        - Landsat 7 (ETM+): ST_B6 (Band 6 - thermal infrared) 
        - Landsat 4/5 (TM): ST_B6 (Band 6 - thermal infrared)
        - MODIS Terra/Aqua: LST_Day_1km (Land Surface Temperature)
        
        Note: Sentinel-3 LST is not available in Google Earth Engine as it requires 
        SLSTR instrument data, but only OLCI is available in GEE.
        
        Returns temperature in Celsius degrees with quality control
        """
        
        # Detectar qué bandas térmicas están disponibles
        band_names = image.bandNames()
        
        # LANDSAT 8/9 (OLI/TIRS) - Usa ST_B10
        st_b10 = ee.List(['ST_B10'])
        has_st_b10 = band_names.containsAll(st_b10)
        
        # LANDSAT 4/5/7 (TM/ETM+) - Usa ST_B6  
        st_b6 = ee.List(['ST_B6'])
        has_st_b6 = band_names.containsAll(st_b6)
        
        # MODIS - Usa LST_Day_1km
        modis_lst = ee.List(['LST_Day_1km'])
        has_modis_lst = band_names.containsAll(modis_lst)
        
        # Aplicar el método correcto según el sensor disponible
        return ee.Algorithms.If(
            has_st_b10,
            # Landsat 8/9: ST_B10 en Kelvin, aplicar factor de escala y convertir a Celsius
            # Factor de escala: 0.00341802, offset: 149.0 (ya en Kelvin)
            image.select('ST_B10').multiply(0.00341802).add(149.0).subtract(273.15).rename(['nd']),
            
            ee.Algorithms.If(
                has_st_b6,
                # Landsat 4/5/7: ST_B6 en Kelvin, aplicar factor de escala y convertir a Celsius  
                # Factor de escala: 0.00341802, offset: 149.0 (ya en Kelvin)
                image.select('ST_B6').multiply(0.00341802).add(149.0).subtract(273.15).rename(['nd']),
                
                ee.Algorithms.If(
                    has_modis_lst,
                    # MODIS: LST_Day_1km - FIXED conversion with quality control
                    self._process_modis_lst(image),
                    
                    # Fallback: devolver imagen con valores nulos si no hay bandas térmicas
                    ee.Image.constant(-9999).rename(['nd']).updateMask(ee.Image.constant(0))
                )
            )
        )

    def _process_modis_lst(self, image):
        """
        Process MODIS LST with robust quality control and error handling
        """
        lst_raw = image.select('LST_Day_1km')
        
        # MODIS LST viene en escala: DN * 0.02 = Kelvin
        # Pero también tiene valores especiales/errores que necesitamos filtrar
        
        # 1. Convertir valores válidos (típicamente 7500-65535 en DN)
        # En Kelvin serían ~150K to 1310K, en Celsius: ~-123°C to +1037°C
        lst_kelvin = lst_raw.multiply(0.02)
        
        # 2. Filtrar valores físicamente imposibles
        # LST válido en Kelvin: ~200K (-73°C) a 373K (100°C) para superficies terrestres normales
        valid_mask = lst_kelvin.gte(200).And(lst_kelvin.lte(373))
        
        # 3. Convertir a Celsius
        lst_celsius = lst_kelvin.subtract(273.15)
        
        # 4. Aplicar máscara de calidad y rango válido
        lst_final = lst_celsius.updateMask(valid_mask)
        
        # 5. Castear explícitamente para evitar problemas de tipo
        return lst_final.toFloat().rename(['nd'])
    
    def get_vci(self, image):
        """
        Vegetation Condition Index (VCI) - Índice de Condición de Vegetación
        
        Compara el NDVI actual con los valores históricos mín/máx para detectar estrés vegetal.
        VCI = (NDVI - NDVImin) / (NDVImax - NDVImin) * 100
        
        Valores:
        - 0-20: Sequía severa
        - 20-40: Sequía moderada  
        - 40-60: Condiciones normales
        - 60-80: Condiciones favorables
        - 80-100: Condiciones muy favorables
        
        Referencias
        -----------
        Kogan, F.N. (1995). Application of vegetation index and brightness temperature 
        for drought detection. Advances in Space Research, 15(11), 91-100.
        
        Note: Esta implementación usa valores aproximados. Para análisis precisos,
        usar estadísticas multi-anuales específicas del área.
        """
        # Calcular NDVI actual
        ndvi = self.get_ndvi(image).select('nd')
        
        # Valores típicos de NDVI mín/máx (aproximados para implementación general)
        # En uso real, estos deberían calcularse de series temporales largas
        ndvi_min = 0.1   # NDVI mínimo histórico típico
        ndvi_max = 0.8   # NDVI máximo histórico típico
        
        # Calcular VCI
        vci = ndvi.subtract(ndvi_min).divide(ndvi_max - ndvi_min).multiply(100)
        
        # Limitar rango 0-100
        vci_clipped = vci.clamp(0, 100)
        
        return vci_clipped.rename(['nd'])

    def get_utfvi(self, image):
        """
        Urban Thermal Field Variance Index (UTFVI) - Versión corregida y simplificada
        
        Fórmula simplificada basada en la relación LST-NDVI:
        UTFVI = (LST - LST_mean) / std_dev
        
        Donde LST_mean se aproxima usando la relación inversa con NDVI.
        
        Valores esperados:
        - > 0.5: Zona urbana fuerte (isla de calor)
        - 0 a 0.5: Zona urbana moderada
        - < 0: Zona vegetada/fresca
        """
        # Verificar banda térmica disponible
        band_names = image.bandNames()
        has_st_b10 = band_names.contains('ST_B10')
        has_st_b6 = band_names.contains('ST_B6')
        
        def process_landsat_89():
            # LST en Celsius
            lst = image.select('ST_B10').multiply(0.00341802).add(149.0).subtract(273.15)
            
            # NDVI
            ndvi = image.normalizedDifference(['Nir', 'Red'])
            
            # UTFVI simplificado: diferencia normalizada entre LST y "temperatura esperada"
            # Temperatura esperada = función inversa del NDVI
            # Alta vegetación (NDVI~0.8) -> baja temperatura esperada
            # Baja vegetación (NDVI~0.2) -> alta temperatura esperada
            
            temp_expected = lst.multiply(0.8).add(ndvi.multiply(-10))  # Relación empírica simple
            
            # UTFVI = (LST_observada - LST_esperada) normalizado
            utfvi = lst.subtract(temp_expected).divide(10)  # Dividir por 10 para normalizar
            
            # Máscara de calidad
            quality_mask = lst.gt(0).And(lst.lt(70)).And(ndvi.gt(-0.2)).And(ndvi.lt(1))
            
            return utfvi.updateMask(quality_mask).rename(['nd'])
        
        def process_landsat_457():
            # LST en Celsius
            lst = image.select('ST_B6').multiply(0.00341802).add(149.0).subtract(273.15)
            
            # NDVI  
            ndvi = image.normalizedDifference(['Nir', 'Red'])
            
            # UTFVI simplificado
            temp_expected = lst.multiply(0.8).add(ndvi.multiply(-10))
            utfvi = lst.subtract(temp_expected).divide(10)
            
            # Máscara de calidad
            quality_mask = lst.gt(0).And(lst.lt(70)).And(ndvi.gt(-0.2)).And(ndvi.lt(1))
            
            return utfvi.updateMask(quality_mask).rename(['nd'])
        
        def fallback_no_thermal():
            return ee.Image.constant(0).rename(['nd']).updateMask(ee.Image.constant(0))
        
        return ee.Algorithms.If(
            has_st_b10,
            process_landsat_89(),
            ee.Algorithms.If(
                has_st_b6,
                process_landsat_457(),
                fallback_no_thermal()
            )
        )

    def get_nbr(self, image):
        """
        Normalized Burn Ratio (NBR) - Ratio Normalizado de Quemadura
        
        Detecta áreas quemadas usando la diferencia entre NIR y SWIR2.
        NBR = (NIR - SWIR2) / (NIR + SWIR2)
        
        Valores típicos:
        - > 0.27: Vegetación densa no quemada
        - 0.1 - 0.27: Vegetación moderada  
        - -0.1 - 0.1: Área quemada reciente
        - < -0.1: Quemadura severa
        
        Para detectar cambios: dNBR = NBR_prefire - NBR_postfire
        
        Referencias
        -----------
        Key, C.H., Benson, N.C. (2006). Landscape Assessment: Ground measure of 
        severity, the Composite Burn Index. FIREMON: Fire effects monitoring and 
        inventory framework. USDA Forest Service.
        """
        return image.normalizedDifference(['Nir', 'Swir2']).rename(['nd'])

    def get_wi2015(self, image):
        """
        Water Index 2015 (WI2015) - Linear water index fitted on Landsat surface reflectance.

        ``1.7204 + 171*Green + 3*Red - 70*NIR - 45*SWIR1 - 71*SWIR2`` on
        surface reflectance in 0-1, as published. Positive values are water.
        Before 1.6.0 the coefficients were divided by 100 on the belief that
        the original ones were meant for unscaled digital numbers; the
        constant was left as is, so the index sat around 1.5 everywhere and
        did not separate water from land.

        References
        ----------
        Fisher, A., Flood, N., Danaher, T. (2016). 
        Comparing Landsat water index methods for automated water classification 
        in eastern Australia. Remote Sensing of Environment, 175, 167-182.
        https://doi.org/10.1016/j.rse.2015.12.055
        """
        return image.expression(
            '1.7204 + 171 * Green + 3 * Red - 70 * NIR - 45 * SWIR1 - 71 * SWIR2', {
            'Green': image.select('Green'),
            'Red': image.select('Red'), 
            'NIR': image.select('Nir'),
            'SWIR1': image.select('Swir1'),
            'SWIR2': image.select('Swir2')
        }).rename(['nd'])

    def get_ndbi(self, image):
        """
        Normalized Difference Built-up Index (NDBI) - Índice Normalizado de Construcción
        
        Detecta áreas urbanas y construidas usando la diferencia entre SWIR1 y NIR.
        NDBI = (SWIR1 - NIR) / (SWIR1 + NIR)
        
        Valores:
        - > 0: Área construida (más alto = más urbano)
        - < 0: Vegetación/agua
        - 0.1 - 0.5: Área urbana típica
        - > 0.5: Área densamente construida
        
        Referencias
        -----------
        Zha, Y., Gao, J., Ni, S. (2003). Use of normalized difference built-up index 
        in automatically mapping urban areas from TM imagery. 
        International Journal of Remote Sensing, 24(3), 583-594.
        """
        return image.normalizedDifference(['Swir1', 'Nir']).rename(['nd'])
    
    # ÍNDICES SENTINEL-2 CON RED EDGE
    def get_ireci(self, image):
        """
        Inverted Red-Edge Chlorophyll Index - Highly sensitive to chlorophyll content.
        
        References
        ----------
        Frampton, W.J., Dash, J., Watmough, G., Milton, E.J. (2013). 
        Evaluating the capabilities of Sentinel-2 for quantitative estimation of biophysical variables in vegetation. 
        ISPRS Journal of Photogrammetry and Remote Sensing, 82, 83-92.
        """
        return image.expression(
            '(Red_Edge3 - Red) / (Red_Edge1 / Red_Edge2)', {
            'Red': image.select('Red'),
            'Red_Edge1': image.select('Red_Edge1'),    # B5
            'Red_Edge2': image.select('Red_Edge2'),    # B6
            'Red_Edge3': image.select('Red_Edge3')     # B7
        }).rename(['nd'])

    def get_mcari(self, image):
        """
        Modified Chlorophyll Absorption Ratio Index - Chlorophyll content with reduced soil influence.
        
        References
        ----------
        Daughtry, C.S.T., Walthall, C.L., Kim, M.S., de Colstoun, E.B., McMurtrey, J.E. (2000). 
        Estimating corn leaf chlorophyll concentration from leaf and canopy reflectance. 
        Remote Sensing of Environment, 74(2), 229-239.
        """
        return image.expression(
            '((Red_Edge1 - Red) - 0.2 * (Red_Edge1 - Green)) * (Red_Edge1 / Red)', {
            'Green': image.select('Green'),
            'Red': image.select('Red'), 
            'Red_Edge1': image.select('Red_Edge1')     # B5
        }).rename(['nd'])

    def get_ndre(self, image):
        """
        Normalized Difference Red Edge - Sensitive to chlorophyll content variations.
        
        References
        ----------
        Gitelson, A., Merzlyak, M.N. (1994). 
        Spectral reflectance changes associated with autumn senescence of Aesculus hippocastanum L. 
        and Acer platanoides L. leaves. Journal of Plant Physiology, 143(3), 286-292.
        """
        return image.normalizedDifference(['Nir', 'Red_Edge1'])  # NIR - Red Edge 1

    def get_reip(self, image):
        """
        Red Edge Inflection Point - Wavelength position of maximum slope in red-edge region.
        
        References
        ----------
        Guyot, G., Baret, F. (1988). 
        Utilisation de la haute resolution spectrale pour suivre l'etat des couverts vegetaux. 
        Proceedings of the 4th International Colloquium on Spectral Signatures of Objects in Remote Sensing, 279-286.
        """
        return image.expression(
            '700 + 40 * ((((Red + Red_Edge3) / 2) - Red_Edge1) / (Red_Edge2 - Red_Edge1))', {
            'Red': image.select('Red'),
            'Red_Edge1': image.select('Red_Edge1'),    # B5
            'Red_Edge2': image.select('Red_Edge2'),    # B6
            'Red_Edge3': image.select('Red_Edge3')     # B7
        }).rename(['nd'])

    def get_psri(self, image):
        """
        Plant Senescence Reflectance Index - Plant senescence and carotenoid/chlorophyll ratio.
        
        References
        ----------
        Merzlyak, M.N., Gitelson, A.A., Chivkunova, O.B., Rakitin, V.Y. (1999). 
        Non‐destructive optical detection of pigment changes during leaf senescence and fruit ripening. 
        Physiologia Plantarum, 106(1), 135-141.
        """
        return image.expression(
            '(Red - Blue) / Red_Edge2', {
            'Blue': image.select('Blue'),
            'Red': image.select('Red'),
            'Red_Edge2': image.select('Red_Edge2')     # B6
        }).rename(['nd'])

    def get_cire(self, image):
        """
        Chlorophyll Index Red Edge - Chlorophyll content estimation using red-edge band.
        
        References
        ----------
        Gitelson, A.A., Gritz, Y., Merzlyak, M.N. (2003). 
        Relationships between leaf chlorophyll content and spectral reflectance and algorithms 
        for non-destructive chlorophyll assessment in higher plant leaves. 
        Journal of Plant Physiology, 160(3), 271-282.
        """
        return image.expression(
            '(Nir / Red_Edge1) - 1', {
            'Red_Edge1': image.select('Red_Edge1'),    # B5
            'Nir': image.select('Nir')
        }).rename(['nd'])

    def get_mtci(self, image):
        """
        MERIS Terrestrial Chlorophyll Index - Chlorophyll content adapted for terrestrial vegetation.
        
        References
        ----------
        Dash, J., Curran, P.J. (2004). 
        The MERIS terrestrial chlorophyll index. International Journal of Remote Sensing, 25(23), 5403-5413.
        """
        return image.expression(
            '(Red_Edge2 - Red_Edge1) / (Red_Edge1 - Red)', {
            'Red': image.select('Red'),
            'Red_Edge1': image.select('Red_Edge1'),    # B5
            'Red_Edge2': image.select('Red_Edge2')     # B6
        }).rename(['nd'])

    def get_s2rep(self, image):
        """
        Sentinel-2 Red Edge Position - Simplified red-edge position estimation for Sentinel-2.
        
        References
        ----------
        Frampton, W.J., Dash, J., Watmough, G., Milton, E.J. (2013). 
        Evaluating the capabilities of Sentinel-2 for quantitative estimation of biophysical variables in vegetation. 
        ISPRS Journal of Photogrammetry and Remote Sensing, 82, 83-92.
        """
        return image.expression(
            '705 + 35 * ((((Red + Red_Edge3) / 2) - Red_Edge1) / (Red_Edge2 - Red_Edge1))', {
            'Red': image.select('Red'),
            'Red_Edge1': image.select('Red_Edge1'),    # B5
            'Red_Edge2': image.select('Red_Edge2'),    # B6
            'Red_Edge3': image.select('Red_Edge3')     # B7
        }).rename(['nd'])

    def get_ndci(self, image):
        """
        Normalized Difference Chlorophyll Index - Optimized for cyanobacteria and chlorophyll-a detection in water.
        
        Formula: (Red_Edge1 - Red) / (Red_Edge1 + Red)
        Uses Red Edge 1 (B5) and Red (B4) - optimized for Sentinel-2.
        
        References
        ----------
        Mishra, S., Mishra, D.R. (2012). 
        Normalized difference chlorophyll index: A novel model for remote estimation of chlorophyll-a 
        concentration in turbid productive waters. Remote Sensing of Environment, 117, 394-406.
        
        Gitelson, A.A., Dall'Olmo, G., Moses, W., Rundquist, D.C., Barrow, T., Fisher, T.R., ... Holz, J. (2008). 
        A simple semi-analytical model for remote estimation of chlorophyll-a in turbid waters: Validation. 
        Remote Sensing of Environment, 112(9), 3582-3593.
        """
        return image.normalizedDifference(['Red_Edge1', 'Red'])
    
    def get_fai(self, image):
        """
        Floating Algae Index (FAI) - Detection of floating algae and cyanobacterial blooms in water bodies.

        Computes a NIR baseline by linear interpolation between the Red and SWIR2 bands
        at the NIR wavelength, then subtracts it from the observed NIR reflectance.
        Positive FAI values indicate floating algae or cyanobacterial surface accumulations;
        negative values correspond to open water.

        Works with surface reflectance (SR) data. Compatible with Sentinel-2, Landsat,
        and MODIS. Not applicable to Sentinel-3 OLCI (no SWIR band available).

        Formula
        -------
        FAI = NIR - NIR_baseline
        NIR_baseline = Red + (SWIR2 - Red) * (lambda_NIR - lambda_red) / (lambda_SWIR2 - lambda_red)

        Sensor-specific center wavelengths (nm) and resulting interpolation factors:
            Sentinel-2 MSI     : lambda_NIR=835.1, lambda_red=664.5, lambda_SWIR2=2202.4 -> factor=0.1109
            Landsat 8/9 OLI    : lambda_NIR=864.7, lambda_red=654.6, lambda_SWIR2=2201.2 -> factor=0.1359
            MODIS Terra/Aqua   : lambda_NIR=858.5, lambda_red=645.0, lambda_SWIR2=2130.0 -> factor=0.1438

        Note: The Landsat collection in ndvi2gif merges OLI (L8/9) and TM/ETM+ (L4-7).
        OLI wavelengths are used as representative values. The interpolation factor for
        TM/ETM+ (~0.113) differs slightly but the impact on FAI is minor.

        References
        ----------
        Hu, C. (2009).
        A novel ocean color index to detect floating algae in the global oceans.
        Remote Sensing of Environment, 113(10), 2118-2129.
        https://doi.org/10.1016/j.rse.2009.05.012
        """
        wavelength_factors = {
            'S2':      (835.1 - 664.5) / (2202.4 - 664.5),  # 0.1109
            'Landsat': (864.7 - 654.6) / (2201.2 - 654.6),  # 0.1359
            'MODIS':   (858.5 - 645.0) / (2130.0 - 645.0),  # 0.1438
        }
        factor = wavelength_factors.get(self.sat, 0.1359)

        return image.expression(
            'NIR - (Red + (Swir2 - Red) * factor)', {
            'NIR':   image.select('Nir'),
            'Red':   image.select('Red'),
            'Swir2': image.select('Swir2'),
            'factor': factor
        }).rename(['nd'])

    # =================================================================
    # ÍNDICES S3 IMPLEMENTADOS (todos basados en radiancia L1B)
    # =================================================================

    def get_oci(self, image):
        """
        OLCI Chlorophyll Index - Custom chlorophyll index using OLCI L1B radiance data.
        
        References
        ----------
        Hu, C., Lee, Z., Franz, B. (2012). 
        Chlorophyll a algorithms for oligotrophic oceans: A novel approach based on three‐band reflectance difference. 
        Journal of Geophysical Research: Oceans, 117(C1), C01011.
        """
        return image.expression(
            '(Red_Edge1 - Red2) / (Red_Edge1 + Red2)', {
            'Red2': image.select('Red2'),              # Oa08 - 665nm
            'Red_Edge1': image.select('Red_Edge1')     # Oa10 - 681.25nm
        }).rename(['nd'])

    def get_tsi(self, image):
        """
        Trophic State Index - Water trophic state classification for eutrophication assessment.
        
        References
        ----------
        Carlson, R.E. (1977). 
        A trophic state index for lakes. Limnology and Oceanography, 22(2), 361-369.
        
        Kratzer, S., Håkansson, B., Sahlin, C. (2003). 
        Assessing Secchi and photic zone depth in the Baltic Sea from satellite data. 
        AMBIO: A Journal of the Human Environment, 32(8), 577-585.
        """
        return image.expression(
            '(Red2 - Blue_Green) / (NIR - Blue_Green)', {
            'Blue_Green': image.select('Blue_Green'),  # Oa04 - 490nm
            'Red2': image.select('Red2'),              # Oa08 - 665nm  
            'NIR': image.select('Nir')                 # Oa12 - 753.75nm
        }).rename(['nd'])

    def get_cdom(self, image):
        """
        Colored Dissolved Organic Matter Index - CDOM absorption assessment in water bodies.
        
        References
        ----------
        Mannino, A., Russ, M.E., Hooker, S.B. (2008). 
        Algorithm development and validation for satellite‐derived distributions of DOC and CDOM in the U.S. Middle Atlantic Bight. 
        Journal of Geophysical Research: Oceans, 113(C7), C07051.
        """
        return image.expression(
            'Blue / Blue_Green', {
            'Blue': image.select('Blue'),              # Oa02 - 412.5nm
            'Blue_Green': image.select('Blue_Green')   # Oa04 - 490nm
        }).rename(['nd'])

    def get_turbidity(self, image):
        """
        Water Turbidity Index - Suspended sediment and water clarity assessment using OLCI bands.
        
        References
        ----------
        Nechad, B., Ruddick, K.G., Park, Y. (2010). 
        Calibration and validation of a generic multisensor algorithm for mapping of total suspended matter in turbid waters. 
        Remote Sensing of Environment, 114(4), 854-866.
        """
        return image.expression(
            'Red2 / Blue_Green', {
            'Red2': image.select('Red2'),              # Oa08 - 665nm
            'Blue_Green': image.select('Blue_Green')   # Oa04 - 490nm
        }).rename(['nd'])

    def get_spm(self, image):
        """
        Suspended Particulate Matter Index - Quantification of suspended particles in water.
        
        References
        ----------
        Binding, C.E., Bowers, D.G., Mitchelson‐Jacob, E.G. (2005). 
        Estimating suspended sediment concentrations from ocean colour measurements in moderately turbid waters; 
        the impact of variable particle scattering properties. Remote Sensing of Environment, 94(3), 373-383.
        """
        return image.expression(
            '(Red3 - Red2) / (Red3 + Red2)', {
            'Red2': image.select('Red2'),    # Oa08 - 665nm
            'Red3': image.select('Red3')     # Oa09 - 673.75nm
        }).rename(['nd'])

    def get_kd490(self, image):
        """
        Diffuse Attenuation Coefficient at 490nm - Water transparency and optical depth assessment.
        
        References
        ----------
        Mueller, J.L. (2000). 
        SeaWiFS algorithm for the diffuse attenuation coefficient K(490) using water-leaving radiances at 490 and 555 nm. 
        SeaWiFS Postlaunch Calibration and Validation Analyses, Part 3, 11, 24-27.
        """
        return image.expression(
            'log(Blue2 / Blue_Green)', {
            'Blue2': image.select('Blue2'),            # Oa03 - 442.5nm
            'Blue_Green': image.select('Blue_Green')   # Oa04 - 490nm
        }).rename(['nd'])

    def get_floating_algae(self, image):
        """
        Floating Algae Index - Detection of floating algae and surface algal blooms.
        
        References
        ----------
        Hu, C. (2009). 
        A novel ocean color index to detect floating algae in the global oceans. 
        Remote Sensing of Environment, 113(10), 2118-2129.
        """
        return image.expression(
            '(NIR - Red_Edge2) / (NIR + Red_Edge2)', {
            'Red_Edge2': image.select('Red_Edge2'),    # Oa11 - 708.75nm
            'NIR': image.select('Nir')                 # Oa12 - 753.75nm
        }).rename(['nd'])

    def get_red_edge_position(self, image):
        """
        Red Edge Position optimized for OLCI - Chlorophyll-sensitive wavelength position indicator.
        
        References
        ----------
        Gower, J., King, S., Borstad, G., Brown, L. (2005). 
        Detection of intense plankton blooms using the 709 nm band of the MERIS imaging spectrometer. 
        International Journal of Remote Sensing, 26(9), 2005-2012.
        """
        return image.expression(
            '681.25 + 27.5 * ((Red2 + Red_Edge2) / 2 - Red_Edge1) / (Red_Edge2 - Red_Edge1)', {
            'Red2': image.select('Red2'),              # Oa08 - 665nm
            'Red_Edge1': image.select('Red_Edge1'),    # Oa10 - 681.25nm
            'Red_Edge2': image.select('Red_Edge2')     # Oa11 - 708.75nm
        }).rename(['nd'])

    def get_fluorescence_height(self, image):
        """
        Chlorophyll Fluorescence Line Height - Natural chlorophyll fluorescence detection.
        
        References
        ----------
        Gower, J., King, S., Borstad, G., Brown, L. (2005). 
        Detection of intense plankton blooms using the 709 nm band of the MERIS imaging spectrometer. 
        International Journal of Remote Sensing, 26(9), 2005-2012.
        """
        return image.expression(
            'Red3 - (Red2 + (Red_Edge1 - Red2) * (673.75 - 665) / (681.25 - 665))', {
            'Red2': image.select('Red2'),              # Oa08 - 665nm
            'Red3': image.select('Red3'),              # Oa09 - 673.75nm
            'Red_Edge1': image.select('Red_Edge1')     # Oa10 - 681.25nm
        }).rename(['nd'])

    def get_water_leaving_reflectance(self, image):
        """
        Water Leaving Reflectance - Simplified approximation of water-leaving radiance contribution.
        
        References
        ----------
        Gordon, H.R., Brown, O.B., Evans, R.H., Brown, J.W., Smith, R.C., Baker, K.S., Clark, D.K. (1988). 
        A semianalytic radiance model of ocean color. Journal of Geophysical Research: Atmospheres, 93(D9), 10909-10924.
        """
        return image.expression(
            'Green / (Blue + Green + Red)', {
            'Blue': image.select('Blue'),      # Oa02 - 412.5nm
            'Green': image.select('Green'),    # Oa05 - 510nm
            'Red': image.select('Red')         # Oa07 - 620nm
        }).rename(['nd'])

    #### ERA5-Land Climate Reanalysis Variables ####

    # Temperature variables
    def get_era5_temperature_2m(self, image):
        """
        Air temperature at 2 meters height.

        Returns temperature in Kelvin. For Celsius: subtract 273.15

        References
        ----------
        ERA5-Land: ECMWF Climate Reanalysis
        https://developers.google.com/earth-engine/datasets/catalog/ECMWF_ERA5_LAND_DAILY_AGGR
        """
        return image.select('temperature_2m').rename('nd')

    def get_era5_dewpoint_temperature_2m(self, image):
        """
        Dewpoint temperature at 2 meters height.

        Temperature at which air becomes saturated with water vapor.
        Units: Kelvin

        References
        ----------
        ERA5-Land: ECMWF Climate Reanalysis
        """
        return image.select('dewpoint_temperature_2m').rename('nd')

    def get_era5_skin_temperature(self, image):
        """
        Skin temperature - Earth surface temperature.

        Temperature of the Earth's surface (land or water).
        Units: Kelvin

        References
        ----------
        ERA5-Land: ECMWF Climate Reanalysis
        """
        return image.select('skin_temperature').rename('nd')

    def get_era5_soil_temperature_level_1(self, image):
        """
        Soil temperature at level 1 (0-7 cm depth).

        Temperature in the topmost soil layer.
        Units: Kelvin

        References
        ----------
        ERA5-Land: ECMWF Climate Reanalysis
        """
        return image.select('soil_temperature_level_1').rename('nd')

    # Precipitation and water balance variables
    def get_era5_total_precipitation_sum(self, image):
        """
        Total daily precipitation (rain + snow combined).

        Units: meters of water equivalent
        Note: Flow band - accumulated daily sum

        References
        ----------
        ERA5-Land: ECMWF Climate Reanalysis
        """
        return image.select('total_precipitation_sum').rename('nd')

    def get_era5_total_evaporation_sum(self, image):
        """
        Total evapotranspiration from land surface.

        Includes evaporation from soil, vegetation, and water bodies.
        Units: meters of water equivalent
        Note: Flow band - accumulated daily sum

        References
        ----------
        ERA5-Land: ECMWF Climate Reanalysis
        """
        return image.select('total_evaporation_sum').rename('nd')

    def get_era5_potential_evaporation_sum(self, image):
        """
        Potential evapotranspiration.

        Maximum evaporation that would occur with unlimited water availability.
        Units: meters
        Note: Flow band - accumulated daily sum

        References
        ----------
        ERA5-Land: ECMWF Climate Reanalysis
        """
        return image.select('potential_evaporation_sum').rename('nd')

    def get_era5_runoff_sum(self, image):
        """
        Total runoff (surface + sub-surface).

        Units: meters
        Note: Flow band - accumulated daily sum

        References
        ----------
        ERA5-Land: ECMWF Climate Reanalysis
        """
        return image.select('runoff_sum').rename('nd')

    def get_era5_surface_runoff_sum(self, image):
        """
        Surface runoff only.

        Water that flows over the land surface.
        Units: meters
        Note: Flow band - accumulated daily sum

        References
        ----------
        ERA5-Land: ECMWF Climate Reanalysis
        """
        return image.select('surface_runoff_sum').rename('nd')

    # Soil moisture variables
    def get_era5_volumetric_soil_water_layer_1(self, image):
        """
        Volumetric soil water content at layer 1 (0-7 cm depth).

        Volume fraction of water in soil (0-1).
        Units: m³/m³

        References
        ----------
        ERA5-Land: ECMWF Climate Reanalysis
        """
        return image.select('volumetric_soil_water_layer_1').rename('nd')

    def get_era5_volumetric_soil_water_layer_2(self, image):
        """
        Volumetric soil water content at layer 2 (7-28 cm depth).

        Volume fraction of water in soil (0-1).
        Units: m³/m³

        References
        ----------
        ERA5-Land: ECMWF Climate Reanalysis
        """
        return image.select('volumetric_soil_water_layer_2').rename('nd')

    def get_era5_volumetric_soil_water_layer_3(self, image):
        """
        Volumetric soil water content at layer 3 (28-100 cm depth).

        Volume fraction of water in soil (0-1).
        Units: m³/m³

        References
        ----------
        ERA5-Land: ECMWF Climate Reanalysis
        """
        return image.select('volumetric_soil_water_layer_3').rename('nd')

    def get_era5_volumetric_soil_water_layer_4(self, image):
        """
        Volumetric soil water content at layer 4 (100-289 cm depth).

        Volume fraction of water in soil (0-1).
        Units: m³/m³

        References
        ----------
        ERA5-Land: ECMWF Climate Reanalysis
        """
        return image.select('volumetric_soil_water_layer_4').rename('nd')

    # Radiation variables
    def get_era5_surface_solar_radiation_downwards_sum(self, image):
        """
        Downward solar radiation at the surface.

        Total incoming shortwave radiation.
        Units: J/m²
        Note: Flow band - accumulated daily sum

        References
        ----------
        ERA5-Land: ECMWF Climate Reanalysis
        """
        return image.select('surface_solar_radiation_downwards_sum').rename('nd')

    def get_era5_surface_net_solar_radiation_sum(self, image):
        """
        Net solar radiation at the surface.

        Incoming minus reflected solar radiation.
        Units: J/m²
        Note: Flow band - accumulated daily sum

        References
        ----------
        ERA5-Land: ECMWF Climate Reanalysis
        """
        return image.select('surface_net_solar_radiation_sum').rename('nd')

    def get_era5_surface_latent_heat_flux_sum(self, image):
        """
        Surface latent heat flux.

        Energy used for evaporation/condensation.
        Units: J/m²
        Note: Flow band - accumulated daily sum

        References
        ----------
        ERA5-Land: ECMWF Climate Reanalysis
        """
        return image.select('surface_latent_heat_flux_sum').rename('nd')

    # Wind and pressure variables
    def get_era5_u_component_of_wind_10m(self, image):
        """
        Eastward wind component at 10 meters height.

        Positive values indicate wind from west to east.
        Units: m/s

        References
        ----------
        ERA5-Land: ECMWF Climate Reanalysis
        """
        return image.select('u_component_of_wind_10m').rename('nd')

    def get_era5_v_component_of_wind_10m(self, image):
        """
        Northward wind component at 10 meters height.

        Positive values indicate wind from south to north.
        Units: m/s

        References
        ----------
        ERA5-Land: ECMWF Climate Reanalysis
        """
        return image.select('v_component_of_wind_10m').rename('nd')

    def get_era5_surface_pressure(self, image):
        """
        Atmospheric pressure at the surface.

        Units: Pascals (Pa)

        References
        ----------
        ERA5-Land: ECMWF Climate Reanalysis
        """
        return image.select('surface_pressure').rename('nd')

    # Snow variables
    def get_era5_snow_depth_water_equivalent(self, image):
        """
        Snow depth in water equivalent.

        Amount of water that would result from melting the snow.
        Units: meters of water equivalent

        References
        ----------
        ERA5-Land: ECMWF Climate Reanalysis
        """
        return image.select('snow_depth_water_equivalent').rename('nd')

    def get_era5_snowfall_sum(self, image):
        """
        Snowfall amount.

        Units: meters of water equivalent
        Note: Flow band - accumulated daily sum

        References
        ----------
        ERA5-Land: ECMWF Climate Reanalysis
        """
        return image.select('snowfall_sum').rename('nd')

    #### ERA5-Land Converted Units (Celsius and L/m²) ####

    # Temperature conversions to Celsius
    def get_era5_temperature_2m_celsius(self, image):
        """
        Air temperature at 2 meters height in Celsius.

        Converted from Kelvin to Celsius (K - 273.15).

        References
        ----------
        ERA5-Land: ECMWF Climate Reanalysis
        https://developers.google.com/earth-engine/datasets/catalog/ECMWF_ERA5_LAND_DAILY_AGGR
        """
        return image.select('temperature_2m').subtract(273.15).rename('nd')

    def get_era5_dewpoint_temperature_2m_celsius(self, image):
        """
        Dewpoint temperature at 2 meters height in Celsius.

        Temperature at which air becomes saturated with water vapor.
        Converted from Kelvin to Celsius (K - 273.15).

        References
        ----------
        ERA5-Land: ECMWF Climate Reanalysis
        """
        return image.select('dewpoint_temperature_2m').subtract(273.15).rename('nd')

    def get_era5_skin_temperature_celsius(self, image):
        """
        Skin temperature (Earth surface) in Celsius.

        Temperature of the Earth's surface (land or water).
        Converted from Kelvin to Celsius (K - 273.15).

        References
        ----------
        ERA5-Land: ECMWF Climate Reanalysis
        """
        return image.select('skin_temperature').subtract(273.15).rename('nd')

    def get_era5_soil_temperature_level_1_celsius(self, image):
        """
        Soil temperature at level 1 (0-7 cm depth) in Celsius.

        Temperature in the topmost soil layer.
        Converted from Kelvin to Celsius (K - 273.15).

        References
        ----------
        ERA5-Land: ECMWF Climate Reanalysis
        """
        return image.select('soil_temperature_level_1').subtract(273.15).rename('nd')

    # Temperature min/max (Kelvin)
    def get_era5_temperature_2m_min(self, image):
        """
        Daily minimum air temperature at 2 meters height.

        Units: Kelvin (K)
        For Celsius use temperature_2m_min_celsius

        References
        ----------
        ERA5-Land: ECMWF Climate Reanalysis
        """
        return image.select('temperature_2m_min').rename('nd')

    def get_era5_temperature_2m_max(self, image):
        """
        Daily maximum air temperature at 2 meters height.

        Units: Kelvin (K)
        For Celsius use temperature_2m_max_celsius

        References
        ----------
        ERA5-Land: ECMWF Climate Reanalysis
        """
        return image.select('temperature_2m_max').rename('nd')

    def get_era5_dewpoint_temperature_2m_min(self, image):
        """
        Daily minimum dewpoint temperature at 2 meters height.

        Units: Kelvin (K)
        For Celsius use dewpoint_temperature_2m_min_celsius

        References
        ----------
        ERA5-Land: ECMWF Climate Reanalysis
        """
        return image.select('dewpoint_temperature_2m_min').rename('nd')

    def get_era5_dewpoint_temperature_2m_max(self, image):
        """
        Daily maximum dewpoint temperature at 2 meters height.

        Units: Kelvin (K)
        For Celsius use dewpoint_temperature_2m_max_celsius

        References
        ----------
        ERA5-Land: ECMWF Climate Reanalysis
        """
        return image.select('dewpoint_temperature_2m_max').rename('nd')

    def get_era5_skin_temperature_min(self, image):
        """
        Daily minimum skin temperature (Earth surface).

        Units: Kelvin (K)
        For Celsius use skin_temperature_min_celsius

        References
        ----------
        ERA5-Land: ECMWF Climate Reanalysis
        """
        return image.select('skin_temperature_min').rename('nd')

    def get_era5_skin_temperature_max(self, image):
        """
        Daily maximum skin temperature (Earth surface).

        Units: Kelvin (K)
        For Celsius use skin_temperature_max_celsius

        References
        ----------
        ERA5-Land: ECMWF Climate Reanalysis
        """
        return image.select('skin_temperature_max').rename('nd')

    def get_era5_soil_temperature_level_1_min(self, image):
        """
        Daily minimum soil temperature at level 1 (0-7 cm depth).

        Units: Kelvin (K)
        For Celsius use soil_temperature_level_1_min_celsius

        References
        ----------
        ERA5-Land: ECMWF Climate Reanalysis
        """
        return image.select('soil_temperature_level_1_min').rename('nd')

    def get_era5_soil_temperature_level_1_max(self, image):
        """
        Daily maximum soil temperature at level 1 (0-7 cm depth).

        Units: Kelvin (K)
        For Celsius use soil_temperature_level_1_max_celsius

        References
        ----------
        ERA5-Land: ECMWF Climate Reanalysis
        """
        return image.select('soil_temperature_level_1_max').rename('nd')

    # Temperature min/max conversions to Celsius
    def get_era5_temperature_2m_min_celsius(self, image):
        """
        Daily minimum air temperature at 2 meters height in Celsius.

        Converted from Kelvin to Celsius (K - 273.15).

        References
        ----------
        ERA5-Land: ECMWF Climate Reanalysis
        """
        return image.select('temperature_2m_min').subtract(273.15).rename('nd')

    def get_era5_temperature_2m_max_celsius(self, image):
        """
        Daily maximum air temperature at 2 meters height in Celsius.

        Converted from Kelvin to Celsius (K - 273.15).

        References
        ----------
        ERA5-Land: ECMWF Climate Reanalysis
        """
        return image.select('temperature_2m_max').subtract(273.15).rename('nd')

    def get_era5_dewpoint_temperature_2m_min_celsius(self, image):
        """
        Daily minimum dewpoint temperature at 2 meters height in Celsius.

        Converted from Kelvin to Celsius (K - 273.15).

        References
        ----------
        ERA5-Land: ECMWF Climate Reanalysis
        """
        return image.select('dewpoint_temperature_2m_min').subtract(273.15).rename('nd')

    def get_era5_dewpoint_temperature_2m_max_celsius(self, image):
        """
        Daily maximum dewpoint temperature at 2 meters height in Celsius.

        Converted from Kelvin to Celsius (K - 273.15).

        References
        ----------
        ERA5-Land: ECMWF Climate Reanalysis
        """
        return image.select('dewpoint_temperature_2m_max').subtract(273.15).rename('nd')

    def get_era5_skin_temperature_min_celsius(self, image):
        """
        Daily minimum skin temperature (Earth surface) in Celsius.

        Converted from Kelvin to Celsius (K - 273.15).

        References
        ----------
        ERA5-Land: ECMWF Climate Reanalysis
        """
        return image.select('skin_temperature_min').subtract(273.15).rename('nd')

    def get_era5_skin_temperature_max_celsius(self, image):
        """
        Daily maximum skin temperature (Earth surface) in Celsius.

        Converted from Kelvin to Celsius (K - 273.15).

        References
        ----------
        ERA5-Land: ECMWF Climate Reanalysis
        """
        return image.select('skin_temperature_max').subtract(273.15).rename('nd')

    def get_era5_soil_temperature_level_1_min_celsius(self, image):
        """
        Daily minimum soil temperature at level 1 (0-7 cm depth) in Celsius.

        Converted from Kelvin to Celsius (K - 273.15).

        References
        ----------
        ERA5-Land: ECMWF Climate Reanalysis
        """
        return image.select('soil_temperature_level_1_min').subtract(273.15).rename('nd')

    def get_era5_soil_temperature_level_1_max_celsius(self, image):
        """
        Daily maximum soil temperature at level 1 (0-7 cm depth) in Celsius.

        Converted from Kelvin to Celsius (K - 273.15).

        References
        ----------
        ERA5-Land: ECMWF Climate Reanalysis
        """
        return image.select('soil_temperature_level_1_max').subtract(273.15).rename('nd')

    # Precipitation and water balance conversions to L/m²
    def get_era5_total_precipitation_sum_lm2(self, image):
        """
        Total daily precipitation in liters per square meter (L/m²).

        Rain + snow combined.
        Converted from meters to L/m² (m × 1000).
        Note: 1 mm = 1 L/m²

        References
        ----------
        ERA5-Land: ECMWF Climate Reanalysis
        """
        return image.select('total_precipitation_sum').multiply(1000).rename('nd')

    def get_era5_total_evaporation_sum_lm2(self, image):
        """
        Total evapotranspiration in liters per square meter (L/m²).

        Includes evaporation from soil, vegetation, and water bodies.
        Converted from meters to L/m² (m × 1000).

        References
        ----------
        ERA5-Land: ECMWF Climate Reanalysis
        """
        return image.select('total_evaporation_sum').multiply(1000).rename('nd')

    def get_era5_potential_evaporation_sum_lm2(self, image):
        """
        Potential evapotranspiration in liters per square meter (L/m²).

        Maximum evaporation with unlimited water availability.
        Converted from meters to L/m² (m × 1000).

        References
        ----------
        ERA5-Land: ECMWF Climate Reanalysis
        """
        return image.select('potential_evaporation_sum').multiply(1000).rename('nd')

    def get_era5_runoff_sum_lm2(self, image):
        """
        Total runoff (surface + sub-surface) in liters per square meter (L/m²).

        Converted from meters to L/m² (m × 1000).

        References
        ----------
        ERA5-Land: ECMWF Climate Reanalysis
        """
        return image.select('runoff_sum').multiply(1000).rename('nd')

    def get_era5_surface_runoff_sum_lm2(self, image):
        """
        Surface runoff in liters per square meter (L/m²).

        Water that flows over the land surface.
        Converted from meters to L/m² (m × 1000).

        References
        ----------
        ERA5-Land: ECMWF Climate Reanalysis
        """
        return image.select('surface_runoff_sum').multiply(1000).rename('nd')

    def get_era5_snowfall_sum_lm2(self, image):
        """
        Snowfall amount in liters per square meter (L/m²).

        Converted from meters of water equivalent to L/m² (m × 1000).

        References
        ----------
        ERA5-Land: ECMWF Climate Reanalysis
        """
        return image.select('snowfall_sum').multiply(1000).rename('nd')

    #### CHIRPS Precipitation Dataset ####

    def get_chirps_precipitation(self, image):
        """
        Daily precipitation from CHIRPS dataset.

        Climate Hazards Group InfraRed Precipitation with Station data (CHIRPS)
        is a quasi-global rainfall dataset combining satellite imagery with
        in-situ station data.

        Units: millimeters per day (mm/d)
        Temporal resolution: Daily (1981-present)
        Spatial resolution: ~5.5 km (0.05°)
        Coverage: 50°S to 50°N latitude

        Use with key='sum' for monthly/seasonal precipitation totals.
        Use with key='mean' for average daily precipitation rates.

        References
        ----------
        Funk, C., Peterson, P., Landsfeld, M. et al. (2015).
        The climate hazards infrared precipitation with stations—a new
        environmental record for monitoring extremes.
        Scientific Data, 2, 150066. https://doi.org/10.1038/sdata.2015.66

        UCSB Climate Hazards Center
        https://developers.google.com/earth-engine/datasets/catalog/UCSB-CHG_CHIRPS_DAILY
        """
        return image.select('precipitation').rename('nd')

    #### Nighttime lights variables ####

    def get_viirs_avg_rad(self, image):
        """
        Average nighttime radiance from the VIIRS Day/Night Band.

        The monthly VCMSLCFG composite: cloud-free average radiance with the
        stray light of the polar summer corrected, which is what makes it
        usable at high latitudes. Moonlight, airglow and fires are *not*
        removed, so a bright pixel is not necessarily a lit one.

        Units: nW/cm²/sr
        Temporal resolution: monthly (2014-present)
        Spatial resolution: ~500 m (15 arc-seconds)

        Use with ``key='mean'`` or ``'median'``; a monthly composite already
        is the average of its cloud-free nights.

        References
        ----------
        Elvidge, C.D., Baugh, K., Zhizhin, M., Hsu, F.C., Ghosh, T. (2017).
        VIIRS night-time lights. International Journal of Remote Sensing,
        38(21), 5860-5879. https://doi.org/10.1080/01431161.2017.1342050

        Earth Observation Group, Colorado School of Mines
        https://developers.google.com/earth-engine/datasets/catalog/NOAA_VIIRS_DNB_MONTHLY_V1_VCMSLCFG
        """
        return image.select('avg_rad').rename('nd')

    def get_viirs_ntl(self, image):
        """
        Daily nighttime radiance, corrected for the moon and the atmosphere.

        The Black Marble product (VNP46A2) removes what the monthly composite
        leaves in: moonlight, atmospheric scattering, terrain and the angle of
        view are all modelled away, so what is left is the light the surface
        actually emitted. That is what makes a *daily* value comparable to the
        next day's, and a map of the date of the maximum meaningful.

        Units: nW/cm²/sr
        Temporal resolution: daily (2012-present)
        Spatial resolution: ~500 m

        Cloudy nights are still cloudy: the quality mask
        (:meth:`mask_viirs_quality`, applied when ``cloud_filter=True``) drops
        the retrievals the algorithm itself grades as poor.

        References
        ----------
        Román, M.O., Wang, Z., Sun, Q. et al. (2018). NASA's Black Marble
        nighttime lights product suite. Remote Sensing of Environment, 210,
        113-143. https://doi.org/10.1016/j.rse.2018.03.017

        https://developers.google.com/earth-engine/datasets/catalog/NASA_VIIRS_002_VNP46A2
        """
        return image.select('DNB_BRDF_Corrected_NTL').rename('nd')

    def get_viirs_ntl_gap_filled(self, image):
        """
        Daily nighttime radiance with the missing nights filled in.

        Same as :meth:`get_viirs_ntl`, but the pixels without a valid
        retrieval are filled from a temporal model instead of being left
        empty. Convenient for a continuous series, dangerous for anything that
        counts observations: a filled value is a prediction, not a
        measurement, and ``quality_flag`` is what tells them apart.

        Units: nW/cm²/sr
        """
        return image.select('Gap_Filled_DNB_BRDF_Corrected_NTL').rename('nd')

    def get_viirs_lunar_irradiance(self, image):
        """
        Lunar irradiance at the time of the VIIRS overpass.

        The moonlight that the BRDF correction of :meth:`get_viirs_ntl`
        removed. Useful to check that a pattern in the lights is not the lunar
        cycle in disguise.

        Units: nW/cm²/sr
        """
        return image.select('DNB_Lunar_Irradiance').rename('nd')

    def get_viirs_quality_flag(self, image):
        """
        Retrieval quality of a VIIRS Black Marble daily pixel.

        Values: 0 high quality (persistent lights), 1 high quality (ephemeral
        lights), 2 poor quality (outlier or cloud contaminated), 255 no
        retrieval. See :meth:`mask_viirs_quality`, which keeps 0 and 1.

        With ``key='count'`` this counts the nights that produced a value in
        each period, which is the coverage behind any daily-lights composite.
        """
        return image.select('Mandatory_Quality_Flag').rename('nd')

    def get_viirs_cloud_mask(self, image):
        """
        Cloud mask flags of a VIIRS Black Marble daily pixel.

        ``QF_Cloud_Mask`` is a bitmask, not a number to average: bits hold day
        or night, cloud confidence, cirrus, snow or ice. Exposed so it can be
        decoded with ``bitwiseAnd`` when the default quality mask is not
        enough.
        """
        return image.select('QF_Cloud_Mask').rename('nd')

    def get_viirs_snow_flag(self, image):
        """
        Snow or ice cover under a VIIRS Black Marble daily pixel.

        Values: 0 no snow or ice, 1 snow or ice. Fresh snow reflects
        artificial light back to the sensor and can double the apparent
        radiance of a city, so a winter peak deserves a look at this band.
        """
        return image.select('Snow_Flag').rename('nd')

    def get_cf_cvg(self, image):
        """
        Number of cloud-free observations behind a nighttime lights composite.

        Present in both VIIRS and DMSP-OLS, with the same meaning in each: how
        many nights actually contributed to the pixel. It is the quality band
        of these datasets — a low count means the radiance of that pixel rests
        on very few nights, which matters in the tropics and in the polar
        winter.

        Units: count of observations
        """
        return image.select('cf_cvg').rename('nd')

    def get_dmsp_avg_vis(self, image):
        """
        Average visible band of the DMSP-OLS annual composite.

        The raw annual average of the cloud-free nights, ephemeral lights
        (fires, aurorae, lightning) included. Use
        :meth:`get_dmsp_stable_lights` for the cleaned version.

        Units: 6-bit digital number (0-63), *not* a physical radiance
        Temporal resolution: annual (1992-2013)
        Spatial resolution: ~1 km (30 arc-seconds)

        The 6-bit sensor saturates over city centres and the DN is not
        calibrated between satellites or years, so a DMSP series needs
        intercalibration before its values can be compared over time.

        References
        ----------
        Elvidge, C.D., Baugh, K.E., Kihn, E.A., Kroehl, H.W., Davis, E.R.
        (1997). Mapping city lights with nighttime data from the DMSP
        Operational Linescan System. Photogrammetric Engineering and Remote
        Sensing, 63(6), 727-734.

        https://developers.google.com/earth-engine/datasets/catalog/NOAA_DMSP-OLS_NIGHTTIME_LIGHTS
        """
        return image.select('avg_vis').rename('nd')

    def get_dmsp_stable_lights(self, image):
        """
        DMSP-OLS lights with the ephemeral ones removed.

        Fires, aurorae and other one-off events are dropped and background
        noise is set to zero, leaving the lights that were there all year:
        cities, industry, gas flares. This is the band most DMSP studies use.

        Units: 6-bit digital number (0-63), background set to 0
        Temporal resolution: annual (1992-2013)
        Spatial resolution: ~1 km (30 arc-seconds)

        See :meth:`get_dmsp_avg_vis` for the saturation and intercalibration
        caveats, which apply here too.
        """
        return image.select('stable_lights').rename('nd')

    def get_dmsp_avg_lights_x_pct(self, image):
        """
        DMSP-OLS average lights weighted by the share of lit nights.

        The product of ``avg_vis`` and the percentage of observations in which
        the pixel was lit, which downweights pixels that only shine
        occasionally. Useful to separate a steadily lit town from a place that
        lights up a few nights a year.

        Units: 6-bit digital number scaled by the lit fraction
        Temporal resolution: annual (1992-2013)
        Spatial resolution: ~1 km (30 arc-seconds)
        """
        return image.select('avg_lights_x_pct').rename('nd')

    #### New methods for SAR indices with normalization option ####

    def _normalize_to_01(self, image):
        """
        Normalize SAR image values using Z-score standardization (mean=0, std=1).
        
        Z-score normalization transforms data to have zero mean and unit variance:
        z = (x - μ) / σ
        
        This method is robust to outliers and preserves the distribution shape
        while making values comparable across different indices and time periods.
        
        Parameters
        ----------
        image : ee.Image
            Single-band image to normalize
            
        Returns
        -------
        ee.Image
            Z-score normalized image (mean≈0, std≈1)
            
        Notes
        -----
        - Values typically range from -3 to +3 (99.7% of data)
        - Negative values indicate below-average
        - Positive values indicate above-average
        - Extreme values (|z| > 3) indicate outliers
        """
        # Compute mean and standard deviation within the ROI
        stats = image.reduceRegion(
            reducer=ee.Reducer.mean().combine(
                reducer2=ee.Reducer.stdDev(),
                sharedInputs=True
            ),
            geometry=self.roi,
            scale=10,
            maxPixels=1e9,
            bestEffort=True
        )
        
        # Extract statistics
        band_name = image.bandNames().get(0)
        mean_val = ee.Number(stats.get(ee.String(band_name).cat('_mean')))
        std_val = ee.Number(stats.get(ee.String(band_name).cat('_stdDev')))
        
        # Convert to constant images for proper image operations
        mean_image = ee.Image.constant(mean_val)
        std_image = ee.Image.constant(std_val)
        
        # Avoid division by zero by adding a small constant where std is 0
        std_image = std_image.where(std_image.eq(0), 0.0001)
        
        # Apply Z-score transformation: (x - mean) / std
        zscore = image.subtract(mean_image).divide(std_image)
        
        return zscore.rename(image.bandNames())
            

    # Funciones SAR actualizadas con parámetro normalize
    @staticmethod
    def _s1_linear(image):
        """
        VV and VH of a Sentinel-1 image as linear power.

        The Sentinel-1 collection is kept in dB, the usual scale for
        backscatter composites, but polarimetric indices are defined on linear
        power: ratios and sums of dB values have no physical meaning.
        """
        return ee.Image(10).pow(image.select(['VV', 'VH']).divide(10))

    def get_rvi(self, image, normalize=False):
        """
        Radar Vegetation Index - More robust vegetation indicator than individual polarizations.
        
        Parameters
        ----------
        image : ee.Image
            Input SAR image with VV and VH bands
        normalize : bool, optional
            If True, normalizes output to [0,1] range. Default is False.
        
        Dual-polarization form ``4 * VH / (VV + VH)``, computed on linear
        power: the collection is in dB, which gives meaningless values here.

        References
        ----------
        Kim, Y., Jackson, T., Bindlish, R., Lee, H., Hong, S. (2012). 
        Radar vegetation index for estimating the vegetation water content of rice and soybean. 
        IEEE Geoscience and Remote Sensing Letters, 9(4), 564-568.

        Nasirzadehdizaji, R., Balik Sanli, F., Abdikan, S., Cakir, Z., Sekertekin,
        A., Ustuner, M. (2019). Sensitivity analysis of multi-temporal Sentinel-1
        SAR parameters to crop height and canopy coverage. Applied Sciences,
        9(4), 655. https://doi.org/10.3390/app9040655
        """
        linear = self._s1_linear(image)
        rvi = linear.expression(
            '4 * VH / (VV + VH)', {
            'VV': linear.select('VV'),
            'VH': linear.select('VH')}).rename(['RVI'])
        
        if normalize:
            rvi = self._normalize_to_01(rvi)
        
        return rvi

    def get_vv(self, image, normalize=False):
        """
        VV Polarization - Vertical transmit, vertical receive. Sensitive to rough surface scattering.
        
        Parameters
        ----------
        image : ee.Image
            Input SAR image with VV band
        normalize : bool, optional
            If True, normalizes output to [0,1] range. Default is False.
        
        References
        ----------
        Ulaby, F.T., Moore, R.K., Fung, A.K. (1986). 
        Microwave Remote Sensing: Active and Passive. Volume 3: From Theory to Applications. 
        Artech House, Norwood, MA.
        """
        vv = image.select('VV').rename(['VV'])
        
        if normalize:
            vv = self._normalize_to_01(vv)
        
        return vv

    def get_vh(self, image, normalize=False):
        """
        VH Polarization - Vertical transmit, horizontal receive. Sensitive to vegetation structure.
        
        Parameters
        ----------
        image : ee.Image
            Input SAR image with VH band
        normalize : bool, optional
            If True, normalizes output to [0,1] range. Default is False.
        
        References
        ----------
        Ulaby, F.T., Moore, R.K., Fung, A.K. (1986). 
        Microwave Remote Sensing: Active and Passive. Volume 3: From Theory to Applications. 
        Artech House, Norwood, MA.
        """
        vh = image.select('VH').rename(['VH'])
        
        if normalize:
            vh = self._normalize_to_01(vh)
        
        return vh

    def get_vv_vh_ratio(self, image, normalize=False):
        """
        VV/VH Ratio - Highly sensitive to structural changes, ideal for crop monitoring and mowing detection.
        
        Parameters
        ----------
        image : ee.Image
            Input SAR image with VV and VH bands
        normalize : bool, optional
            If True, normalizes output to [0,1] range. Default is False.
        
        References
        ----------
        Mascolo, L., Lopez‐Sanchez, J.M., Vicente‐Guijalba, F., Nunziata, F., Migliaccio, M., Mazzarella, G. (2016). 
        A complete procedure for crop phenology estimation with PolSAR data based on the complex Wishart classifier. 
        IEEE Transactions on Geoscience and Remote Sensing, 54(11), 6505-6515.

        Computed on linear power. In dB the same information is the
        difference ``VV - VH``; dividing two dB values, as before 1.6.0,
        gives a number with no physical meaning.
        """
        linear = self._s1_linear(image)
        ratio = linear.expression(
            'VV / VH', {
            'VV': linear.select('VV'),
            'VH': linear.select('VH')}).rename(['RATIO'])
        
        if normalize:
            ratio = self._normalize_to_01(ratio)
        
        return ratio

    def get_dpsvi(self, image, normalize=False):
        """
        Modified Dual-Polarization SAR Vegetation Index (DPSVIm).

        ``(VV**2 + VV * VH) / sqrt(2)`` on linear power. It is the modified
        form of Periasamy's DPSVI proposed by dos Santos et al. (2021): the
        original, ``(VV_max - VV + VH) / sqrt(2) * (VV + VH) / VV * VH``,
        depends on ``VV_max``, the maximum VV of the whole scene, so it is not
        a per-pixel index and changes with the extent; the modified form is
        computed pixel by pixel and compares across scenes and dates. Before
        1.6.0 this method returned ``(VV - VH) / (VV + VH)`` on dB values,
        which is neither.

        Parameters
        ----------
        image : ee.Image
            Input SAR image with VV and VH bands
        normalize : bool, optional
            If True, normalizes output to [0,1] range. Default is False.
        
        References
        ----------
        dos Santos, E.P., da Silva, D.D., do Amaral, C.H. (2021). Vegetation
        cover monitoring in tropical regions using SAR-C dual-polarization
        index: seasonal and spatial influences. International Journal of Remote
        Sensing, 42(19), 7581-7609. https://doi.org/10.1080/01431161.2021.1959955

        Periasamy, S. (2018). Significance of dual polarimetric synthetic
        aperture radar in biomass retrieval: An attempt on Sentinel-1. Remote
        Sensing of Environment, 217, 537-549.
        https://doi.org/10.1016/j.rse.2018.09.003
        """
        linear = self._s1_linear(image)
        dpsvi = linear.expression(
            '(VV ** 2 + VV * VH) / sqrt(2)', {
            'VV': linear.select('VV'),
            'VH': linear.select('VH')}).rename(['DPSVI'])
        
        if normalize:
            dpsvi = self._normalize_to_01(dpsvi)
        
        return dpsvi

    def get_rfdi(self, image, normalize=False):
        """
        Radar Forest Degradation Index - Forest disturbance and degradation monitoring.
        
        Parameters
        ----------
        image : ee.Image
            Input SAR image with VV and VH bands
        normalize : bool, optional
            If True, normalizes output to [0,1] range. Default is False.
        
        References
        ----------
        Ningthoujam, R.K., Balzter, H., Tansey, K., Feldpausch, T.R., Mitchard, E.T., Wani, A.A., Joshi, P.K. (2018). 
        Relationships of S-1 C-band SAR backscatter with forest cover, height and aboveground biomass at multiple spatial scales across four forest types. 
        Remote Sensing, 10(9), 1442.

        Mitchard, E.T.A., Saatchi, S.S., White, L.J.T., Abernethy, K.A., Jeffery,
        K.J., Lewis, S.L., et al. (2012). Mapping tropical forest biomass with
        radar and spaceborne LiDAR in Lopé National Park, Gabon. Biogeosciences,
        9(1), 179-191. https://doi.org/10.5194/bg-9-179-2012

        Defined by Mitchard et al. as ``(HH - HV) / (HH + HV)``; Sentinel-1
        carries VV/VH, so the dual-pol adaptation ``(VV - VH) / (VV + VH)`` is
        used, on linear power (before 1.6.0 it was ``(VV - VH) / VV`` on dB).
        """
        linear = self._s1_linear(image)
        rfdi = linear.expression(
            '(VV - VH) / (VV + VH)', {
            'VV': linear.select('VV'),
            'VH': linear.select('VH')
        }).rename(['RFDI'])
        
        if normalize:
            rfdi = self._normalize_to_01(rfdi)
        
        return rfdi

    def get_vsdi(self, image, normalize=False):
        """
        Vegetation Scattering Diversity Index (experimental).

        ``sqrt((VV - VH)**2 + (VV + VH)**2)``, computed on the dB values of
        the collection. No publication defining this index has been found —
        it was previously attributed to Periasamy (2018), who defines DPDD,
        IDPDD, VDDPI and DPSVI but not this — so treat it as an experimental
        feature without a physical interpretation, and prefer the documented
        indices for anything that has to be reported.

        Parameters
        ----------
        image : ee.Image
            Input SAR image with VV and VH bands (dB)
        normalize : bool, optional
            If True, normalizes output to [0,1] range. Default is False.
        """
        vsdi = image.expression(
            'sqrt((VV - VH) ** 2 + (VV + VH) ** 2)', {
            'VV': image.select('VV'),
            'VH': image.select('VH')
        }).rename(['VSDI'])
        
        if normalize:
            vsdi = self._normalize_to_01(vsdi)
        
        return vsdi

    
    # Export methods (same as original)
    def get_export_single(self, image, name='mycomposition.tif', crs='EPSG:4326', scale=10):
        """
        Export a single Earth Engine image to a GeoTIFF file.
        
        Exports a specific composite image or any Earth Engine image to the local
        file system as a GeoTIFF raster. This method is useful for exporting
        individual images, statistical summaries, or custom analyses derived
        from the temporal composites.
        
        Parameters
        ----------
        image : ee.Image
            Earth Engine image to export. Can be any single or multi-band image,
            including temporal composites, statistical summaries, or processed
            derivatives from get_year_composite().
            
        name : str, optional
            Output filename including extension. Should end with '.tif' for
            GeoTIFF format. The file will be saved in the current working directory.
            Default is 'mycomposition.tif'.
            
        crs : str, optional
            Coordinate Reference System for the output raster in EPSG format.
            Common options:
            - 'EPSG:4326': WGS84 Geographic (lat/lon)
            - 'EPSG:3857': Web Mercator
            - 'EPSG:32633': UTM Zone 33N (example)
            Default is 'EPSG:4326'.
            
        scale : int or float, optional
            Output pixel resolution in meters. Should match or be appropriate
            for the input satellite data:
            - Sentinel-2: 10-20m
            - Landsat: 30m
            - MODIS: 500m
            - Sentinel-1: 10m
            - Sentinel-3: 300m
            Default is 10.
            
        Examples
        --------
        >>> # Export a single year composite
        >>> processor = NdviSeasonality(sat='S2', index='ndvi', periods=4)
        >>> collection = processor.get_year_composite()
        >>> single_image = collection.first()
        >>> processor.get_export_single(single_image, 'ndvi_2020_seasonal.tif')
        
        >>> # Export temporal statistics
        >>> mean_composite = collection.mean()
        >>> processor.get_export_single(mean_composite, 'ndvi_multiyear_mean.tif', scale=20)
        
        >>> # Export specific band or analysis
        >>> summer_only = single_image.select('summer')
        >>> processor.get_export_single(summer_only, 'summer_ndvi_2020.tif')
        
        >>> # Export with custom CRS and higher resolution
        >>> processor.get_export_single(
        ...     image=single_image,
        ...     name='high_res_composite.tif',
        ...     crs='EPSG:32633',  # UTM projection
        ...     scale=10
        ... )
        
        >>> # Export derived analysis
        >>> seasonal_range = single_image.select('summer').subtract(single_image.select('winter'))
        >>> processor.get_export_single(seasonal_range, 'seasonal_amplitude.tif')
        
        Notes
        -----
        **File Output:**
        - Files are saved to the current working directory (os.getcwd())
        - GeoTIFF format with embedded CRS and geotransform information
        - Multi-band images preserve all bands in a single file
        - Pixel values maintain original data type and scaling
        
        **Performance Considerations:**
        - Export time depends on: image size × number of bands × scale resolution
        - Large regions or high resolutions may approach Earth Engine limits
        - Consider using export_with_fishnet() for very large areas
        - Processing occurs on Earth Engine servers, then downloads locally
        
        **Common Use Cases:**
        - Single image export from temporal analysis
        - Statistical summaries (mean, max, std) across time series
        - Specific band extraction for focused analysis
        - Custom mathematical operations on composites
        - Quality control and validation sample export
        
        Raises
        ------
        ee.EEException
            If Earth Engine encounters processing errors, authentication issues,
            or export limitations (memory, timeout, etc.).
            
        OSError
            If there are local file system issues (permissions, disk space, etc.).
            
        ValueError
            If parameters are invalid (unsupported CRS, negative scale, etc.).
            
        See Also
        --------
        get_export : Export entire time series automatically
        get_year_composite : Generate temporal composites for export
        export_with_fishnet : Export large areas using tiled approach
        
        References
        ----------
        Earth Engine Export Documentation:
        https://developers.google.com/earth-engine/guides/exporting
        """
        # Construct full file path in current working directory
        filename = os.path.join(os.getcwd(), name)
        
        # Export image using geemap wrapper for Earth Engine export
        geemap.ee_export_image(
            image, 
            filename=filename, 
            scale=scale, 
            crs=crs, 
            region=self.roi, 
            file_per_band=False
        ) 
        
        # Provide user feedback
        print('Image have been exported')
        
    def get_export(self, crs='EPSG:4326', scale=10):
        """
        Export all temporal composite images to individual GeoTIFF files with descriptive names.
        
        Processes the complete time series analysis and exports each year's composite
        as a separate multi-band GeoTIFF file. Filenames are automatically generated
        using a descriptive pattern that includes the index, statistical method, and year.
        
        This method orchestrates the entire workflow: generates temporal composites,
        creates meaningful filenames, and exports all results in a single operation.
        
        Parameters
        ----------
        crs : str, optional
            Coordinate Reference System for output rasters in EPSG format.
            Applied to all exported files. Common options:
            - 'EPSG:4326': WGS84 Geographic (lat/lon) - universal compatibility
            - 'EPSG:3857': Web Mercator - web mapping applications
            - 'EPSG:32633': UTM Zone 33N - metric measurements (example)
            Default is 'EPSG:4326'.
            
        scale : int or float, optional
            Output pixel resolution in meters for all exports. Should be appropriate
            for the satellite sensor:
            - Sentinel-2: 10-20m (native resolution)
            - Landsat: 30m (native resolution)
            - MODIS: 500m (native resolution)  
            - Sentinel-1: 10m (native resolution)
            - Sentinel-3: 300m (native resolution)
            Default is 10.
            
        Returns
        -------
        None
            Files are saved to disk with progress information printed to console.
            
        Notes
        -----
        **Automatic Filename Generation:**
        Filenames follow the pattern: `{sat}_{index}_{statistic}_{year}.tif`
        
        Examples:
        - 'ndvi_max_2020.tif' (NDVI with maximum reducer)
        - 'evi_median_2021.tif' (EVI with median reducer)  
        - 'vh_p90_2022.tif' (VH SAR with 90th percentile)
        - 'rvi_mean_2019.tif' (RVI SAR with mean reducer)
        
        **Statistical Method Naming:**
        - 'max': Maximum value reducer
        - 'median': Median value reducer
        - 'mean': Mean value reducer
        - 'p{N}': Percentile reducer (e.g., 'p90' for 90th percentile)
        
        **Multi-band Structure:**
        Each exported file contains multiple bands representing temporal periods:
        - 4 periods: ['winter', 'spring', 'summer', 'autumn']
        - 12 periods: ['january', 'february', ..., 'december'] 
        - 24 periods: ['p1', 'p2', ..., 'p24']
        - Custom: ['p1', 'p2', ..., 'pN']
        
        **Processing Workflow:**
        1. Clear previous results and initialize processing
        2. Generate temporal composites using get_year_composite()
        3. For each successful year:
        a. Generate descriptive filename
        b. Export multi-band composite to GeoTIFF
        c. Provide progress feedback
        4. Report completion summary
        
        Examples
        --------
        >>> # Basic export with default settings
        >>> processor = NdviSeasonality(
        ...     sat='S2', 
        ...     index='ndvi', 
        ...     start_year=2020, 
        ...     end_year=2023,
        ...     key='max'
        ... )
        >>> processor.get_export()
        # Exports: ndvi_max_2020.tif, ndvi_max_2021.tif, ndvi_max_2022.tif
        
        >>> # SAR analysis with percentile statistics
        >>> sar_processor = NdviSeasonality(
        ...     sat='S1',
        ...     index='vh', 
        ...     key='percentile',
        ...     percentile=90,
        ...     periods=12
        ... )
        >>> sar_processor.get_export(scale=20)
        # Exports: vh_p90_2020.tif, vh_p90_2021.tif, etc.
        
        >>> # Monthly analysis with UTM projection
        >>> monthly_processor = NdviSeasonality(
        ...     sat='Landsat',
        ...     index='evi',
        ...     periods=12,
        ...     key='median'
        ... )
        >>> monthly_processor.get_export(crs='EPSG:32633', scale=30)
        # Exports: evi_median_2020.tif, evi_median_2021.tif, etc.
        
        >>> # Multi-sensor comparison workflow
        >>> sensors = ['S2', 'Landsat', 'MODIS']
        >>> for sat in sensors:
        ...     proc = NdviSeasonality(sat=sat, index='ndvi', start_year=2020, end_year=2022)
        ...     proc.get_export()
        # Creates separate files for each sensor-year combination
        
        Performance Notes
        -----------------
        **Processing Time Factors:**
        - Time range: (end_year - start_year + 1) × periods
        - Spatial extent: ROI area × scale resolution
        - Sensor complexity: SAR < Optical (due to preprocessing)
        - Statistical method: max ≈ median < mean < percentile
        
        **Optimization Strategies:**
        - Use appropriate scale for sensor (avoid unnecessary upsampling)
        - Consider seasonal periods (4) vs monthly (12) for time vs detail trade-off  
        - For very large areas, use export_with_fishnet() instead
        - Process smaller time ranges for iterative analysis
        
        **Output File Characteristics:**
        - Format: GeoTIFF with embedded georeferencing
        - Compression: Default LZW compression for smaller files
        - Data type: Float32 for index values, preserving fractional precision
        - NoData handling: Masked pixels preserved from input data
        
        Raises
        ------
        ee.EEException
            If Earth Engine processing fails due to memory limits, timeout,
            authentication issues, or invalid parameters.
            
        OSError
            If local file system issues occur (insufficient disk space,
            permission errors, invalid file paths).
            
        RuntimeError  
            If no valid temporal composites can be generated (e.g., no satellite
            data available for specified time range and region).
            
        Warnings
        --------
        - Years with insufficient data are automatically skipped with console warnings
        - Large time ranges or high-resolution analyses may approach computation limits
        - Existing files with same names will be overwritten without warning
        
        See Also
        --------
        get_export_single : Export individual images with custom names
        get_year_composite : Generate the temporal composites that are exported
        export_with_fishnet : Alternative for very large spatial extents
        get_gif : Create animated visualizations instead of static exports
        
        Examples of Post-Export Analysis
        --------------------------------
        >>> # Load exported files for further analysis
        >>> import rasterio
        >>> with rasterio.open('ndvi_max_2020.tif') as src:
        ...     data = src.read()  # Shape: (bands, height, width)
        ...     summer_data = src.read(3)  # Read summer band specifically
        
        >>> # Multi-temporal analysis
        >>> import numpy as np
        >>> years = range(2020, 2023)
        >>> summer_trend = []
        >>> for year in years:
        ...     with rasterio.open(f'ndvi_max_{year}.tif') as src:
        ...         summer_trend.append(src.read(3))  # Summer band
        >>> trend = np.array(summer_trend)
        >>> mean_summer = np.mean(trend, axis=0)
        """
        # Initialize processing by clearing previous results
        self.imagelist = []
        
        # Generate temporal composites for all years
        self.get_year_composite()
        
        # Get count of successfully processed years
        count = len(self.imagelist)
        print(count)
        
        # Export each year's composite with descriptive filename
        for n in range(count):
            year = self.start_year + n
            image = self.imagelist[n]
            
            # Generate descriptive filename: index_statistic_year.tif
            # Handle special case for percentile reducer
            if self.key == 'percentile':
                stat_name = f'p{self.percentile}'
            else:
                stat_name = self.key
                
            # Create filename with descriptive pattern
            name = f'{self.sat}_{self.index}_{stat_name}_{year}.tif'
            filename = os.path.join(os.getcwd(), name)
            
            # Provide progress feedback
            print('Exporting {}'.format(filename), '\n')
            
            # Export multi-band composite to GeoTIFF
            geemap.ee_export_image(
                image, 
                filename=filename, 
                scale=scale, 
                crs=crs, 
                region=self.roi, 
                file_per_band=False
            ) 
        
        # Provide completion summary
        print('All the images in the collection have been exported')

    def _default_scale_for_sat(self) -> int:
        """
        Return a sensible default pixel scale (meters) based on the configured sensor.

        It is the native resolution of the collection each sensor is read
        from: MOD09A1 reflectance is 500 m, ERA5-Land 0.1 degrees and CHIRPS
        0.05 degrees (their size at the equator).
        """
        sat = (self.sat or "").upper()
        if sat.startswith("S2") or sat.startswith("S1"):
            return 10
        if sat.startswith("L"):   # Landsat 5/7/8/9
            return 30
        if sat.startswith("MOD"):
            return 500
        if sat.startswith("S3"):
            return 300
        if sat == "ERA5":
            return 11132
        if sat == "CHIRPS":
            return 5566
        return 30
    
    def export_to_drive(
        self,
        image,
        description: str,
        *,
        region=None,
        scale: Optional[int] = None,
        crs: str = "EPSG:4326",
        folder: Optional[str] = None,
        file_format: str = "GeoTIFF",
        format_options: Optional[dict] = None,
        max_pixels: int = int(1e13),
        file_dimensions: Optional[int] = None,
        clip_region: bool = True):
        """
        Export an ee.Image to Google Drive as a GeoTIFF (batch task).

        Parameters
        ----------
        image : ee.Image
            Earth Engine image to export (e.g., a classified map or a composite).
        description : str
            Export task name shown in the Earth Engine Tasks panel.
        region : ee.Geometry, optional
            Export region. If None, uses ``self.roi``.
        scale : int, optional
            Pixel resolution in meters. If None, inferred from the sensor via
            :meth:`_default_scale_for_sat`.
        crs : str, optional
            Target projection (e.g., ``"EPSG:4326"``). Default is WGS84.
        folder : str, optional
            Google Drive folder name. If None, uses the Drive root.
        file_format : str, optional
            Output format (e.g., ``"GeoTIFF"``). Default ``"GeoTIFF"``.
        format_options : dict, optional
            Additional format options for the exporter (e.g., compression).
            Example: ``{"cloudOptimized": True, "compression": "LZW"}``.
        max_pixels : int, optional
            Maximum number of pixels allowed by Earth Engine. Default ``1e13``.
        file_dimensions : int, optional
            Maximum pixel dimension per side for each output tile. If provided,
            large exports are split into multiple files (e.g., 8192 or 10000).
        clip_region : bool, optional
            If True, ``image`` is clipped to ``region`` before exporting.

        Returns
        -------
        ee.batch.Task
            The started Earth Engine export task.

        Raises
        ------
        ValueError
            If `image` is not provided.
        ee.EEException
            If Earth Engine fails to create the export task.

        Notes
        -----
        - This is a batch export: monitor progress in the Earth Engine Tasks panel.
        - Use ``clip_region=True`` and/or ``file_dimensions`` to keep file sizes manageable.
        - For large regions, consider compression via ``format_options``.
        """
        if image is None:
            raise ValueError("`image` must be an ee.Image.")

        if scale is None:
            scale = self._default_scale_for_sat()

        if region is None:
            region = self.roi

        # Optionally clip to reduce processing/size
        if clip_region and region is not None:
            image = image.clip(region)

        kwargs = {
            "image": image,
            "description": description,
            "region": region,
            "scale": scale,
            "maxPixels": max_pixels,
            "fileFormat": file_format,
            "crs": crs,
        }

        if folder:
            kwargs["folder"] = folder
        if format_options:
            kwargs["formatOptions"] = format_options
        if file_dimensions:
            kwargs["fileDimensions"] = file_dimensions

        task = ee.batch.Export.image.toDrive(**kwargs)
        task.start()
        return task
    
    def export_to_asset(
        self,
        image,
        asset_id: str,
        *,
        description: Optional[str] = None,
        region=None,
        scale: Optional[int] = None,
        crs: Optional[str] = None,
        crs_transform: Optional[list] = None,   # <-- NUEVO (opcional)
        pyramiding_policy: Optional[dict] = None,
        max_pixels: int = int(1e13),
        overwrite: bool = False,
        clip_region: bool = True):
        """
        Export an ee.Image to an Earth Engine Asset (batch task).

        Parameters
        ----------
        image : ee.Image
            Earth Engine image to export.
        asset_id : str
            Full asset path, e.g., "users/yourname/ndvi2gif/landcover_2022".
        description : str, optional
            Task name. If None, a name is derived from `asset_id`.
        region : ee.Geometry, optional
            Export region. If None, uses ``self.roi``.
        scale : int, optional
            Pixel resolution in meters. If None, inferred from the sensor via
            :meth:`_default_scale_for_sat`.
        crs : str, optional
            Target projection. If None, uses the image's default.
        crs_transform : list, optional
            Affine transform (6 or 9 numbers) instead of `scale`. Mutually
            exclusive with `scale`.
        pyramiding_policy : dict, optional
            Dict mapping band name to policy (e.g., {"class": "mode"}).
        max_pixels : int, optional
            Maximum number of pixels allowed by Earth Engine. Default 1e13.
        overwrite : bool, optional
            If True, tries to delete any existing asset with the same ID first.
        clip_region : bool, optional
            If True, the image is clipped to `region` (or self.roi) before export.

        Returns
        -------
        ee.batch.Task
            The started Earth Engine export task.
        """
        if image is None or not asset_id:
            raise ValueError("`image` and `asset_id` are required.")

        if scale is None and crs_transform is None:
            scale = self._default_scale_for_sat()

        if region is None:
            region = self.roi

        if description is None:
            description = asset_id.split("/")[-1]

        # Recorta para reducir tamaño y evitar procesar fuera del ROI
        if clip_region and region is not None:
            image = image.clip(region)

        if overwrite:
            try:
                ee.data.deleteAsset(asset_id)
            except Exception:
                # ok si no existía
                pass

        kwargs = {
            "image": image,
            "description": description,
            "assetId": asset_id,
            "region": region,
            "maxPixels": max_pixels,
        }
        if scale is not None:
            kwargs["scale"] = scale
        if crs is not None:
            kwargs["crs"] = crs
        if crs_transform is not None:
            kwargs["crsTransform"] = crs_transform
        if pyramiding_policy:
            kwargs["pyramidingPolicy"] = pyramiding_policy

        task = ee.batch.Export.image.toAsset(**kwargs)
        task.start()
        return task

    def get_gif(self, name='mygif.gif', bands=None):
        """
        Create an animated GIF showing the temporal evolution of period composites.

        Generates a video animation (downloaded as GIF) where each frame corresponds
        to one year and each RGB image uses selected period composites (e.g.,
        winter-spring-summer). Internally calls :meth:`get_year_composite` and
        exports the video using :func:`geemap.download_ee_video`. Optionally, an
        annotated version is created with year labels and a progress bar using
        :func:`geemap.add_text_to_gif`.

        Parameters
        ----------
        name : str, optional
            Name of the output file (``.gif``). Saved in the current working
            directory. Default: ``'mygif.gif'``.
        bands : list of str or None, optional
            Band names to use as ``[R, G, B]`` in the animation. Must be valid
            periods (e.g., ``'winter'``, ``'spring'`` or ``'p1'``, ``'p2'``).
            If ``None``, the first three period names are used
            (``self.period_names[:3]``).

        Notes
        -----
        * **SAR sensors (S1):** a typical dB scale is applied (``min=-25``, ``max=0``).
        * **Optical sensors:** a typical vegetation index range is applied (``min=0.15``, ``max=0.85``).
        * Output video is set to ``dimensions=768`` pixels and ``framesPerSecond=10`` for a balance between clarity and file size.

        Raises
        ------
        ee.EEException
            If processing or export fails (memory limits, runtime errors,
            authentication issues).
        ValueError
            If `bands` does not contain exactly 3 valid periods present in
            ``self.period_names``.
        OSError
            If a file system error occurs when writing the GIF locally.

        See Also
        --------
        get_year_composite : Generates yearly composites used in the animation.
        get_export : Exports static multi-band images instead of an animation.
        geemap.download_ee_video : Handles Earth Engine video export.
        geemap.add_text_to_gif : Adds annotations (year and progress bar).

        Examples
        --------
        >>> # Seasonal composite with Sentinel-2 (RGB = winter, spring, summer)
        >>> NdviSeasonality(sat='S2', index='ndvi', periods=4, start_year=2020, end_year=2023).get_gif('ndvi_seasons.gif')
        >>> # Monthly SAR composite (using specific months as RGB)
        >>> proc = NdviSeasonality(sat='S1', index='vh', periods=12)
        >>> proc.get_gif('vh_monthly.gif', bands=['march', 'june', 'september'])
        """
        # Set default bands if none specified (first 3 periods)
        if bands is None:
            bands = self.period_names[:3]  # Use first 3 periods by default
            
        # Initialize processing by clearing previous results
        self.imagelist = []
        
        # Generate temporal composites for animation
        self.get_year_composite()
        
        # Create output file path
        out_gif = os.path.join(os.getcwd(), name)
        
        # Configure visualization parameters based on sensor type
        if self.sat == 'S1':
            # SAR sensor - use decibel scaling appropriate for backscatter
            video_args = {
                'dimensions': 768,          # Fixed output resolution
                'region': self.roi,         # Spatial extent
                'framesPerSecond': 10,      # Animation speed
                'bands': bands,             # RGB band selection
                'min': -25,                 # Typical SAR minimum (dB)
                'max': 0,                   # Typical SAR maximum (dB)
                'gamma': [1, 1, 1]         # Linear gamma correction
            }
        else:
            # Optical sensors - use reflectance scaling for vegetation indices
            video_args = {
                'dimensions': 768,          # Fixed output resolution  
                'region': self.roi,         # Spatial extent
                'framesPerSecond': 10,      # Animation speed
                'bands': bands,             # RGB band selection
                'min': 0.15,               # Vegetation index minimum
                'max': 0.85,               # Vegetation index maximum  
                'gamma': [1, 1, 1]         # Linear gamma correction
            }
        
        # Generate basic animated GIF using Earth Engine video export
        geemap.download_ee_video(self.get_year_composite(), video_args, out_gif)
        
        # Create annotated version with year labels and progress bar
        texted_gif = out_gif[:-4] + '_texted.gif'
        geemap.add_text_to_gif(
            out_gif,                    # Input GIF path
            texted_gif,                 # Output annotated GIF path  
            xy=('5%', '90%'),          # Text position (bottom-left)
            text_sequence=self.start_year,  # Starting year for labels
            font_size=30,               # Text size
            font_color='#ffffff',       # White text color
            add_progress_bar=False,     # No progress bar overlay
            duration=300                # Frame duration in milliseconds
        )

    def export_with_fishnet(self, image, name_prefix='composite', scale=10, crs='EPSG:4326'):
        """
        Export large images using a tiled fishnet approach to overcome memory limitations.
        
        Divides the ROI into regular grid tiles and exports each tile separately.
        Useful for very large spatial extents that exceed Earth Engine's 
        single-image export limits.
        
        Parameters
        ----------
        image : ee.Image
            Earth Engine image to export (typically from get_year_composite()).
        name_prefix : str, optional
            Prefix for output filenames. Files named as '{prefix}_tile_{id}.tif'.
            Default is 'composite'.
        scale : int, optional
            Output pixel resolution in meters. Default is 10.
        crs : str, optional
            Coordinate reference system. Default is 'EPSG:4326'.
            
        Notes
        -----
        Tile size: 50km × 50km for scale ≥ 30m, 25km × 25km for finer scales.
        Only exports tiles that intersect with the ROI geometry.
        
        References
        ----------
        Adapted from Earth Engine large-area export strategies.
        """
        import math
        
        # Determine tile size based on scale
        tile_km = 50 if scale >= 30 else 25
        tile_m = tile_km * 1000
        
        # Get ROI bounds
        bounds = self.roi.bounds()
        coords = bounds.coordinates().getInfo()[0]
        xmin, ymin = coords[0]
        xmax, ymax = coords[2]
        
        # Calculate grid dimensions
        x_steps = math.ceil((xmax - xmin) * 111320 / tile_m)
        y_steps = math.ceil((ymax - ymin) * 110540 / tile_m)
        tile_id = 0
        
        # Generate and export tiles
        for i in range(x_steps):
            for j in range(y_steps):
                # Calculate tile bounds
                x0 = xmin + (i * tile_m / 111320)
                y0 = ymin + (j * tile_m / 110540)
                x1 = x0 + tile_m / 111320
                y1 = y0 + tile_m / 110540
                cell = ee.Geometry.Rectangle([x0, y0, x1, y1])
                
                # Export only if tile intersects ROI
                if self.roi.intersects(cell, ee.ErrorMargin(1)).getInfo():
                    region = cell.intersection(self.roi, ee.ErrorMargin(1))
                    filename = f"{name_prefix}_tile_{tile_id}.tif"
                    tile_id += 1
                    print(f'Exporting tile {tile_id} to {filename}')
                    geemap.ee_export_image(
                        image.clip(region),
                        filename=os.path.join(os.getcwd(), filename),
                        scale=scale,
                        region=region,
                        crs=crs,
                        file_per_band=False
                    )
        print('All tiles have been exported.')

    def get_stats(self, image=None, geom=None, name=None, stat='MEDIAN', scale=10, to_file=False):
        """
        Compute zonal statistics for temporal composites within specified geometries.
        
        Automatically processes all years in the time series, computing statistics
        for each year separately. If a single image is provided, processes only that image.
        
        Parameters
        ----------
        image : ee.Image, ee.ImageCollection, or None, optional  
            Input for statistical analysis. Can be:
            - None: Uses get_year_composite() to process entire time series (default)
            - ee.Image: Single image (e.g., one year composite)
            - ee.ImageCollection: Custom collection to process
        geom : ee.Geometry, str, or None, optional
            Geometry defining zones for statistics. Can be:
            - None: Use self.roi (default)
            - str: Path to shapefile (.shp) or GeoJSON (.geojson)
            - ee.Geometry: Earth Engine geometry object
        name : str or None, optional
            Output filename prefix. If None, uses 'zonal_stats'.
            For multi-year: creates '{name}_{year}.shp' files.
        stat : str, optional
            Statistical method: 'MEAN', 'MEDIAN', 'MAX', 'MIN', 'STDDEV'.
            Default is 'MEDIAN'.
        scale : int, optional
            Pixel resolution for analysis in meters. Default is 10.
        to_file : bool, optional
            If True, saves results as shapefile(s). Default is False.
            
        Returns
        -------
        dict or geopandas.GeoDataFrame
            - Single image: Returns GeoDataFrame with statistics
            - Multiple years: Returns dict {year: GeoDataFrame} with results per year
        """
        # Determine input: single image, collection, or generate new collection
        if image is None:
            # Generate complete time series
            collection = self.get_year_composite()
            process_collection = True
        elif hasattr(image, 'size'):  # It's a collection
            collection = image
            process_collection = True
        else:  # It's a single image
            single_image = image
            process_collection = False
        
        # Determine the zones to analyse. reduceRegions wants a
        # FeatureCollection and computes one row per feature, which is the
        # whole point of zonal statistics: calling .geometry() on a collection
        # of polygons dissolves them into a single zone and drops every
        # attribute with them (one row for 49 reservoirs, and no name on it)
        if geom is None:
            roi = self.roi
        elif isinstance(geom, str):
            if geom.endswith('.shp'):
                roi = geemap.shp_to_ee(geom)
            elif geom.endswith('.geojson'):
                roi = geemap.geojson_to_ee(geom)
            else:
                raise ValueError("Path must be to a .shp or .geojson file.")
        else:
            roi = geom

        if isinstance(roi, ee.Geometry):
            roi = ee.FeatureCollection([ee.Feature(roi)])
        elif isinstance(roi, ee.Feature):
            roi = ee.FeatureCollection([roi])

        # Set default output name
        if name is None:
            name = 'zonal_stats'

        # Process single image
        if not process_collection:
            out_shp = os.path.join(os.getcwd(), name + '.shp')
            
            # Use correct Earth Engine syntax: image.reduceRegions(collection=polygons)
            stats = single_image.reduceRegions(
                collection=roi,
                reducer=getattr(ee.Reducer, stat.lower())(),
                scale=scale
            )
            gdf = geemap.ee_to_gdf(stats)
            
            if to_file:
                gdf.to_file(out_shp)
                print(f'Saved as {out_shp}')
            
            return gdf
        
        # Process collection (multiple years)
        results = {}
        collection_list = collection.toList(collection.size())
        collection_size = collection.size().getInfo()
        
        for i in range(collection_size):
            year = self.start_year + i
            yearly_image = ee.Image(collection_list.get(i))
            
            print(f"Processing year {year}...")
            
            # Compute statistics for this year using correct syntax
            try:
                stats = yearly_image.reduceRegions(
                    collection=roi,
                    reducer=getattr(ee.Reducer, stat.lower())(),
                    scale=scale
                )
                gdf = geemap.ee_to_gdf(stats)
                
                # Store results
                results[year] = gdf
                print(f"Year {year} processed successfully: {gdf.shape[0]} features")
                
            except Exception as e:
                print(f"Error processing year {year}: {e}")
                continue
            
            # Save to file if requested
            if to_file:
                year_filename = f"{name}_{year}.shp"
                out_shp = os.path.join(os.getcwd(), year_filename)
                gdf.to_file(out_shp)
                print(f'Saved {year} statistics as {out_shp}')
        
        print(f"Completed processing. Available years: {list(results.keys())}")
        return results