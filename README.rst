====================================================
Multi Seasonal Remote Sensing Indexes Composites
====================================================

.. image:: https://i.imgur.com/Ikg9jwP.gif


Ndvi2Gif is a python library to create Seasonal Composites based on several statistics applied to some Remote Sensings datastes.
This tool uses `Google Earth Engine API <https://github.com/google/earthengine-api>`_ and the amazing
`Geemap package <https://github.com/giswqs/geemap>`_, to create yearly
compositions based on different statistics. We also have added `deimsPy <https://pypi.org/project/deims/>` to get the boundaries of all eLTER sites. So now, you can choose between a shapefile, a map draw or
just use an eLTER DeimsID to get the boundaries for your seasonal composite index. 

This tool have been updated in the framework of `eLTER H2020 <https://github.com/google/earthengine-api>`_ and 
`SUMHAL <https://lifewatcheric-sumhal.csic.es/descripcion-del-proyecto/>`_ projects, as the main input to 
`PhenoPy <https://github.com/JavierLopatin/PhenoPY/tree/master>`_ python package, 
which is the library that we use to get the phenometrics derived from the seasonal vegetation composites.

.. image:: https://camo.githubusercontent.com/5c734dbb4d997c26304b31db1426732e9497e4f9a49acbd0c8bbf0f9a99c462c/68747470733a2f2f692e696d6775722e636f6d2f5376394c66596a2e706e67


The stats includes at this point are:

* Maximun
* Mean
* Median 
* Percentile 90
* Percentile 95 

The indexes available at present are:

* NDVI
* EVI
* GNDVI 
* SAVI 
* NDWI 
* AEWI
* AEWINSH
* NDSI
* NBRI


And last, the available datasets are the following: 

* **Sentinel**

  Sentinel 1: https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S1_GRD

  Sentinel 2 https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S2_HARMONIZED

* **Landsat**

  Landsat 4 TM: https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LT04_C02_T1_L2   
                      
  Landsat 5 TM: https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LT05_C02_T1_L2    
                      
  Landsat 7 ETM+: https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LE07_C02_T1_L2   
                       
  Landsat 8 OLI: https://developers.google.com/earth-engine/datasets/catalog/landsat-8

  Landsat 9 OLI: https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LC09_C02_T1_L2
                      
* **MODIS**           
                      
  MOD09A1: https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MOD09A1            

It is possible to create a combination of any of these statistics, indices and datasets. By default, Maximum `NDVI <https://en.wikipedia.org/wiki/Normalized_difference_vegetation_index>`__ is used 
as seasonal reducer in order to avoid clouds and cloud shadows. However, we have added others statistic to choice when instantiating the class. 
Max remains the default, but sometimes median gives a
better visual result, specially with Landsat 4 and 5 that sometimes have band errors 
that can affect NDVI results. Percentile 90 is a good compromise between max and median. 

Landsat collections and MODIS datasets are Surface Reflectance (SR) data, while
Sentinel 2 is Top of Atmosphere Reflectance (TOA) data. This is
because Surface Reflectance for Sentinel 2, is only available since
2017 but since 2015 for TOA. 

So, this process generates a raster with 4 (Autumn, Winter, Spring, Summer), 12 (january, febreuary, march, ... ) or 24 (p1, p2, ..., p24) 
bands for every year in the chosen time period. 

Beyond a nice gif, a lot of information can be obtained with this kind of multi seasonal Vegetation Indexes approach. 
Knowing the pair Seasonal Index-Raster band that you chose for your gif, and having colour formation in mind (graphic below), 
you could tell which is the phenology, and therfore the crop or every parcel, and even how it changes through the years.  
White colours means high NDVI values for the three periods chosen for the vizParams (perennial vegetation), black colour means low NDVI values, 
such as permanent water bodies, sand, impervious surfaces, etc...

Since we have added SAR data, maybe is no longer correct saying this is an NDVI tool, but with SAR the meaning s very similar to the NDVI approach, in this case we get higher return values when plants are bigger, and very low values for baresoil. So, at the end is another way to have a multi-temporal look at crop growth. 

.. image:: https://i.imgur.com/tq4aMBv.jpg

Last, you have the choice to download the yearly seasonal index composites as tiff files into your computer, 
in case you want the data for further analysis. Also, it have been noticed that Google Earth Engine reducers are 
really nice to create gorgeous multi-year composties, even for very large areas with MODIS, e.g. median seasonal NDVI 
for whole Africa between 2001 and 2020. So, besides the automatic export for each year, you also have the chance to export 
your favourite multi-year compostion in a single file. 



Installation
============


This tiny and humble python class depends on geemap, so geemap will be installed for you. Also it could be a good idea install first geemap in a python environment (you can see the details here: `geemap install) <https://github.com/giswqs/geemap#installation>`_ and later install ndvi2gif in that environment via pip:

.. code:: python

  pip install ndvi2gif
 


Usage
=====


This is intend to be executed in a notebook and in tandem with a geemap Map object, so you could navigate around the map 
and pick up your region of interest just by drawing a shape, and visualizing different dates and band combinations directly on 
the map. However, you could just run it in a command line and pass it a DeimsID, a shapefile or a geojson as roi, and ask for the gif or 
for the geotiff rasters.


Please, see the `example notebook <https://github.com/Digdgeo/Ndvi2Gif/blob/master/ndvi2gif/ndvi2gif_notebook_example.ipynb>`_ 

.. code:: python

    import geemap
    from ndvi2gif import NdviSeasonality
    
    #You could need a first login to sart with python Earth Enginelogin 
    ee.Initialize()
    
    #Create the Map Object to choose he rois
    Map = geemap.Map()
    Map.add_basemap('Google Satellite')
    Map
    
    #Set the roi to last drawn feature
    roi = Map.draw_last_feature
    
    #Instance ndvi2gif
    #Three different examples here to instantiate the class
    myclass = NdviSeasonality(roi)
    myclass2 = NdviSeasonality(roi, 2014, 2020, 'Landsat')
    myclass3 = NdviSeasonality(roi, 2010, 2015, 'MODIS', key='median')
    
    #Maybe you feel like playing with the Map and see different colour/season combination efore generate the gif
    vizParams = {'bands': ['summer', 'autumn', 'winter'], 'min': 0, 'max': 0.7, 'gamma': [0.95, 1.1, 1]}
    Map.addLayer(show, vizParams, 'mycropsfirstviz')
    
    #Notice that you also can use the Earh Engine amazing analysis capabilities
    wintermax = myclass.get_year_composite().select('winter').max()
    median = myclass.get_year_composite().median()
    Map.addLayer(wintermax, {'min': 0, 'max': 0.75}, 'winterMax')
    Map.addLayer(median, {'min': 0.1, 'max': 0.8}, 'median')
    
    #To get the gif, ust use the method. 
    myclass.get_gif()
    
    #Last, you can export your yearly seasonal NDVI composites to your computer
    myclass.get_export() 



ToDo list
=========


* Add masking capablities based on NDVI values to show real color composite in the background. Is it that possible?
* Add seasons dates as parameters that can be easily modified
* Add a method to easily export multi-yearly composites



Contributions
=============


Yes, please! git pulls will be welcome, even those related to my english grammar...

