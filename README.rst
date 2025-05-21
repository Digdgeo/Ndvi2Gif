Ndvi2Gif. Multi Seasonal Remote Sensing Indexes Composites
==========================================================

.. image:: https://img.shields.io/pypi/v/ndvi2gif.svg
   :target: https://pypi.org/project/ndvi2gif/

.. image:: https://github.com/Digdgeo/Ndvi2Gif/actions/workflows/python-publish.yml/badge.svg
   :target: https://github.com/Digdgeo/Ndvi2Gif/actions/workflows/python-publish.yml

.. image:: https://i.imgur.com/Ikg9jwP.gif

Ndvi2Gif is a python library to create Seasonal Composites based on several statistics applied to some Remote Sensing datasets.  
This tool uses `Google Earth Engine API <https://github.com/google/earthengine-api>`_ and the amazing  
`Geemap package <https://github.com/giswqs/geemap>`_, to create yearly  
compositions based on different statistics. We also have added `deimsPy <https://pypi.org/project/deims/>`_ to get the boundaries of all eLTER sites.  
So now, you can choose between a shapefile, a map draw, or just use an eLTER DeimsID to get the boundaries for your seasonal composite index. 

This tool has been updated in the framework of the `eLTER H2020 <https://elter-project.eu>`_ and  
`SUMHAL <https://lifewatcheric-sumhal.csic.es/descripcion-del-proyecto/>`_ projects, as the main input to  
`PhenoPy <https://github.com/JavierLopatin/PhenoPY/tree/master>`_,  
which is the library that we use to get the phenometrics derived from the seasonal vegetation composites.

.. image:: https://i.imgur.com/Sv9LfYj.png

The statistics included at this point are:

* Maximum
* Mean
* Median 
* Percentile 90
* Percentile 95 

The indices available at present are:

* NDVI
* EVI
* GNDVI 
* SAVI 
* NDWI 
* AEWI
* AEWINSH
* NDSI
* NBRI
* NDMI

Available datasets:

* **Sentinel**
  * Sentinel 1: https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S1_GRD
  * Sentinel 2: https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S2_HARMONIZED

* **Landsat**
  * Landsat 4 TM: https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LT04_C02_T1_L2   
  * Landsat 5 TM: https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LT05_C02_T1_L2    
  * Landsat 7 ETM+: https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LE07_C02_T1_L2   
  * Landsat 8 OLI: https://developers.google.com/earth-engine/datasets/catalog/landsat-8
  * Landsat 9 OLI: https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LC09_C02_T1_L2

* **MODIS**           
  * MOD09A1: https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MOD09A1            

It is possible to create a combination of any of these statistics, indices, and datasets.  
By default, maximum `NDVI <https://en.wikipedia.org/wiki/Normalized_difference_vegetation_index>`__ is used as a seasonal reducer in order to avoid clouds and cloud shadows.  
However, we have added other statistics to choose from when instantiating the class.  
Max remains the default, but sometimes median gives a better visual result, especially with Landsat 4 and 5 that sometimes have band errors.  
Percentile 90 is a good compromise between max and median. 

Landsat collections and MODIS datasets are Surface Reflectance (SR) data, while Sentinel 2 is Top of Atmosphere Reflectance (TOA) data.  
This is because Surface Reflectance for Sentinel 2 is only available since 2017, while TOA is available since 2015.

This process generates a raster with 4 (Autumn, Winter, Spring, Summer), 12 (January–December) or 24 (p1–p24) bands for every year in the chosen time period.  

Beyond a nice gif, a lot of information can be obtained with this kind of multi-seasonal vegetation index approach.  
Knowing the pair Seasonal Index–Raster band that you chose for your gif, and having color interpretation in mind,  
you can infer phenology, crop type, or detect changes through the years.  
White colors mean high NDVI values in the three seasons (perennial vegetation), black colors mean low NDVI (water, sand, impervious surfaces, etc).

Since we have added SAR data, it may no longer be accurate to call this an NDVI-only tool.  
However, SAR behaves similarly: high return values for large plants and low values for bare soil.  
So it remains a valid approach for multi-temporal crop growth analysis.

.. image:: https://i.imgur.com/tq4aMBv.jpg

You also have the option to download the yearly seasonal index composites as GeoTIFF files for further analysis.  
Google Earth Engine reducers also allow you to create beautiful multi-year composites, even for large areas like Africa.  
For example, you could generate a median seasonal NDVI for all of Africa between 2001 and 2020.  
So, besides exporting each year individually, you can export a multi-year composite in a single file.

Installation
============

This Python package depends on `geemap`, which will be installed automatically.  
It is recommended to first install `geemap` in a dedicated environment (see instructions here: `geemap install <https://github.com/giswqs/geemap#installation>`_)  
and then install ndvi2gif via pip:

.. code:: python

    pip install ndvi2gif

Usage
=====

This is intended to be used in a Jupyter Notebook alongside a geemap Map object,  
so you can navigate, draw your region of interest, and preview your data.  
You can also run it from the command line using a shapefile, GeoJSON, or eLTER DeimsID.

See the `example notebook <https://github.com/Digdgeo/Ndvi2Gif/blob/master/ndvi2gif/ndvi2gif_notebook_example.ipynb>`_.

.. code:: python

    import geemap
    from ndvi2gif import NdviSeasonality
    import ee

    ee.Initialize()

    Map = geemap.Map()
    Map.add_basemap('Google Satellite')
    Map

    roi = Map.draw_last_feature

    myclass = NdviSeasonality(roi)
    myclass2 = NdviSeasonality(roi, 2014, 2020, 'Landsat')
    myclass3 = NdviSeasonality(roi, 2010, 2015, 'MODIS', key='median')

    vizParams = {'bands': ['summer', 'autumn', 'winter'], 'min': 0, 'max': 0.7, 'gamma': [0.95, 1.1, 1]}
    Map.addLayer(myclass.get_year_composite(), vizParams, 'mycropsfirstviz')

    wintermax = myclass.get_year_composite().select('winter').max()
    median = myclass.get_year_composite().median()
    Map.addLayer(wintermax, {'min': 0, 'max': 0.75}, 'winterMax')
    Map.addLayer(median, {'min': 0.1, 'max': 0.8}, 'median')

    myclass.get_gif()
    myclass.get_export()

Contributions
=============

Yes, please!  
Feel free to contribute to this project in any way you like — feedback, issues, and pull requests are all welcome!

