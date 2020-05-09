==========================
Multi Seasonal NDVI to GIF
==========================

ndvi2gif is a python script for creating seasonal NDVI compositions
gifs. This is just a small python class that uses `Google Earth Engine
API <https://github.com/google/earthengine-api>`_ and the amazing
`Geemap package <https://github.com/giswqs/geemap>`_, to create yearly
compositions based on the maximum
`NDVI <https://en.wikipedia.org/wiki/Normalized_difference_vegetation_index>`__
value reached in every seasons. Maximum NDVI is used as season reducer, in order to avoid clouds and
cloud shadows. However, we have added 'perc_90' (percentile 90) and 'median' as others available statistic choice
when instantiating the class. Max remains the default, but sometimes median gives a
better visual result, specially with Landsat 4 and 5 that sometimes has band errors
that can affect NDVI results. Percentile 90 is a good compromise between max and median. 

So, this process generates a raster with 4 bands (Winter, Spring, Summer and
Autumn) for every year in the chosen time period.  

Satellite collections available are:

* Sentinel 1-C-band Synthetic Aperture Radar 10 m pixel VH cross-polarization ('sar') | 2014-10-03 - Present
* Sentinel 2-MSI 10 m pixel NDVI ('Sentinel') | 2015-06-23 - Present
* Landsat 4-TM, 5-TM, 7-ETM+ and 8-OLI 30 m pixel NDVI ('Landsat') | 1982 - Present 
* MODIS MOD09Q1 v006 250 m pixel NDVI ('MODIS') | 2000-03-05 - Present

Landsat collections datasets and MODIS, are Surface Reflectance (SR) data, while
Sentinel 2 is Top of Atmosphere Reflectance (TOA) dataset. This is
because Surface Reflectance for Sentinel 2 data, is only available since
2017. 

If everything runs well, you should get a GIF similar to those ones that
you can find in the pics folder of this repo. Actually, you will get 2
gifs, one of them named "mygif_texted.gif", which add year as text to
the gif. Here you can see an example close to Seville, where you can
tell the blue colours (blue band in this example is summer) showing paddy
fields over a marsh area (summer crops). Outside the marshes, the colours
green and yellow predominate,showing winter crops such as cereals. You
can also realize the intermediate colours for different crops.

.. image:: https://i.imgur.com/UZDptan.gif


Beyond the nice gif, a lot of information can be obtained with this kind of multi seasonal NDVI approach. Knowing the pair NDVI season-Raster band that you chose for your gif, and having colour formation in mind (graphic below), you could tell which is the phenology, and therfore the crop or every parcel, and even how it changes through the years.  White colours means high NDVI values for the three chosen seasons (perennial vegetation), black colour means low NDVI values, such as permanent water bodies, sand, impervious surfaces, etc...

Since we have added SAR data, maybe is no longer correct saying this is an NDVI tool, but with SAR the meaning s very similar to the NDVI approach, in this case we get higher return values when plants are bigger, and very low values for baresoil. So, at the end is another way to have a multi-temporal look at crop growth. 

.. image:: https://i.imgur.com/tq4aMBv.jpg

Last, you have the choice to download the yearly ndvi composites as tiff files into your computer, in case you want the data for further analysis. Also, it have been noticed that Google Earth Engine reducers are really nice to create gorgeous multi-year composties, even for very large areas with MODIS, e.g. median seasonal NDVI for whole Africa between 2001 and 2020. So, besides the automatic export for each year, you also have the chance to export your favourite multi-year compostion in a single file. 



Installation
============


This tiny and humble python class depends on geemap, so geemap will be installed for you. Also it could be a good idea install first geemap in a python environment (you can see the details here: `geemap install) <https://github.com/giswqs/geemap#installation>`_ and later install ndvi2gif in that environment via pip:

.. code:: python

  pip install ndvi2gif
 


Usage
=====


This is intend to be executed in a notebook and in tandem with a geemap Map object, so you could travel around the map and pick up your region of interest just by drawing a shape, and visualizing different dates and band combinations directly on the map. However, you could just run it in a command line and pass it a shapefile or a geojson as roi, and ask for the gif or for the geotiff rasters.


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

