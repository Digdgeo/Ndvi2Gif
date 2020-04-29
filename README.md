## Multi Seasonal NDVI to GIF

ndvi2gif is a python script for creating nice seasonal NDVI compositions Gifs. This is just a small scrip that uses [Google Earh Engine API](https://github.com/google/earthengine-api) and the amazing [Geemap library](https://github.com/giswqs/geemap), to create yearly compositions based on the maximun NDVI value reached in every seasons. We use the maximun to avoid clouds and cloud shadows.  
So, we will have a raster with 4 bands (Winter, Spring, Summer and Autumn) for every year in the time period that you choose. Right now, you can choose between Sentinel 2-MSI and Landsat satellites (4&5-TM, 7-ETM+ and 8-OLI), so in first case you have data to play with from 2015 until the present, and in case you choose Landsat you could go from 1984 until present. In case of Landsat we use Surface Reflectance (SR) Datasets, in case of Sentinel 2 we use Top of Atmosphere Reflectance (TOA) data, because Surface Reflectance is only available since 2017. Anyway, as long as we are using an index, atmosferic correction is not so important. Besides, the main goal of this script is just to get a fancy gif.

If everythigns runs well, you should get a GIF similar to those ones that you can find in the pics folder of this repo. Actually, you will get 2 gifs, one of them named  "mygif_texted.gif", which add year as text to the gif. 

![Image](https://github.com/Digdgeo/GEE_Playground/blob/master/pics/LosPalacios_Spain.gif "Los Palacios, Seville")


More over the nice gif, you can get a lot of information with this kind of multi seasonal NDVI approach. Knowing the pair NDVI season-Raster band that you chose for your gif, and having color formation in mind (graphic below), you could tell which is the phenology, and therfore the crop or every parcel, and how it changes thru the years.  White colors means high NDVI values for the three seasons that you chose in your gif, black color means low NDVI values, such as permanent water bodies, sand, impervous surfaces, etc...

<p align="center"> 
<img src="https://i.stack.imgur.com/tKETN.png">
</p>

Last, you have the choice to download the yearly ndvi composites to your computer, in case you want that data for further analysis.  

### Installation

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install ndvi2gif.

```bash
pip install ndvi2gif
```
### Usage

See the [example notebook](https://github.com/Digdgeo/GEE_Playground/blob/master/ndvi2gif_notebook_example.ipynb) in this repo.

```python
#Imports libraries
import ee
import geemap
from ndvi2gif import ndvi_seasonality

#You could need a first login to sart with python Earth Enginelogin 
ee.Initialize()

#Create the Map Object to choose he rois
Map = geemap.Map()
Map.add_basemap('Google Satellite')
Map

#Set the roi to last drawn feature
roi = Map.draw_last_feature

#Instance ndvi2gif
myclassname = ndvi_seasonality(roi)

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
```

### ToDo list

* Add masking capablities based on NDVI values to show real color composite in the background. Is it that possible?
* Add MODIS dataset
* Add more tools to this package
* Add seasons dates as parameters that can be easily modified

### Contributing

Pull requests are welcome!

## License

[MIT](https://choosealicense.com/licenses/mit/)
