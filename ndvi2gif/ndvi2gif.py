import os
import ee
import geemap
import deims

def scale_OLI(image):
        opticalBands = image.select(['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7']).multiply(0.0000275).add(-0.2).rename(['Blue', 'Green', 'Red', 'Nir', 'Swir1', 'Swir2'])
        #thermalBands = image.select('ST_B.*').multiply(0.00341802).add(149.0)
        return image.addBands(opticalBands, None, True)#.addBands(thermalBands, None, True)
    
def scale_ETM(image):
    opticalBands = image.select(['SR_B1','SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7']).multiply(0.0000275).add(-0.2).rename(['Blue', 'Green', 'Red', 'Nir', 'Swir1', 'Swir2'])
    #thermalBands = image.select('ST_B.*').multiply(0.00341802).add(149.0)
    return image.addBands(opticalBands, None, True)#.addBands(thermalBands, None, True)

class NdviSeasonality:
    
   
    '''Class to generate NDVI seasonal compoositions gifs. Seasons are defined like fixed parameters and we use ee.
    Reducer Max to get the maximun NDVI reached in every seasons. Then we combine then in a raster with 4 bands,
    one band per season. Color combination will show phenology over the seasons and over the years.

    Args:
        roi (ee.Geometry): Geometry to cross with satellite collections. Can be taken from roi on map, geojson or shapefile.
        start_year (int): First year to look for.
        end_year (int): End year to look for.
        sat (ee.ImageCollection): Available Sentinel 2 (default. 2015-present) and Landsat 4 TM, 5 TM, 7 ETM+, 8 OLI (1984-present).
    Returns:
        gif: Gif file with ndvi seasonal yearly composites.
        images: GeoTIFF with ndvi seasonal yearly composites.
    '''
                                                            
    def __init__(self, roi=None, periods=4, start_year=2016, end_year=2020, sat='S2', key='max', index='ndvi'):


        print('There we go again...')     
        # Here we get the roi. Valid inputs are draws on the map, shapefiles or geojson
        self.roi = roi
        if self.roi is None:
            # When no geometry is passed we use a default area over Donana Natural Space
            self.roi = ee.Geometry.Polygon(
                [[[-6.766047, 36.776586], 
                  [-6.766047, 37.202186], 
                  [-5.867729, 37.202186], 
                  [-5.867729, 36.776586], 
                  [-6.766047, 36.776586]]], None, False)
            
        elif isinstance(self.roi, str):
            if self.roi.endswith('.shp'):
                self.roi = geemap.shp_to_ee(self.roi).geometry()
            elif self.roi.endswith('.geojson'):
                self.roi = geemap.geojson_to_ee(self.roi).geometry()
            # Let's do the DEIMS part
            elif self.roi.startswith('deimsid'):
                # Spanish joke with Don Quijote de a mancha 
                print('Con DEIMS hemos topado amigo Sancho...')
                id_ = self.roi.split('/')[-1]
                gdf = deims.getSiteBoundaries(id_)
                self.roi = geemap.geopandas_to_ee(gdf)
            else:
                print('It seems that your path is broken. Remember that await for shapefiles or geojson')

        else:

            if not isinstance(roi, ee.Geometry):

                try:
                    self.roi = self.roi.geometry()
                except Exception as e:
                    print('Could not convert the provided roi to ee.Geometry')
                    print(e)
                    return
        
        self.periods = periods
        self.start_year = start_year
        self.end_year = end_year
        if sat not in ['S2', 'Landsat', 'MODIS', 'S1']:
            print('You should choose one from Sentinel 2 (S2), Landsat, MODIS or Sentinel 1 (S1)')
        else:
            self.sat = sat

        self.sat = sat
        if key not in ['max', 'median', 'perc_90', 'perc_95', 'mean']:
            print('Please choose between max, median, mean, perc_90 or perc_95 as available stats')
        else:
            self.key=key
        self.imagelist = []
        self.index = index
        self.d = {'ndvi': self.get_ndvi, 'ndwi': self.get_ndwi, 'evi': self.get_evi, 'savi': self.get_savi,
                  'gndvi': self.get_gndvi, 'avi': self.get_avi, 'nbri': self.get_nbri, 'ndsi': self.get_ndsi,
                  'aweinsh': self.get_aweinsh, 'awei': self.get_awei, 'ndmi': self.get_ndmi}

        #Periods!
        # Here we define the periods splitting the year in 4, 12 or 24 equal parts. Feel free to change the dates in case your are looking for different seasons
        
        #Seasons
        self.winter = ['-01-01', '-03-31']
        self.spring = ['-04-01', '-06-30']
        self.summer = ['-07-01', '-09-30']
        self.autumn = ['-10-01', '-12-31']

        #Yearly
        self.january = ['-01-01', '-01-31']
        self.february = ['-02-01', '-02-28']
        self.march = ['-03-01', '-03-31']
        self.april = ['-04-01', '-04-30']
        self.may = ['-05-01', '-05-31']
        self.june = ['-06-01', '-06-30']
        self.july = ['-07-01', '-07-31']
        self.august = ['-08-01', '-08-31']
        self.september = ['-09-01', '-09-30']
        self.october = ['-10-01', '-10-31']
        self.november = ['-11-01', '-11-30']
        self.december = ['-12-01', '-12-31']

        #15 days
        self.p1 = ['-01-01', '-01-15']
        self.p2 = ['-01-16', '-01-31']
        self.p3 = ['-02-01', '-02-15']
        self.p4 = ['-02-16', '-02-28']
        self.p5 = ['-03-01', '-03-15']
        self.p6 = ['-03-16', '-03-31']
        self.p7 = ['-04-01', '-04-15']
        self.p8 = ['-04-16', '-04-30']
        self.p9 = ['-05-01', '-05-15']
        self.p10 = ['-05-16', '-05-31']
        self.p11 = ['-06-01', '-06-15']
        self.p12 = ['-06-16', '-06-30']
        self.p13 = ['-07-01', '-07-15']
        self.p14 = ['-07-16', '-07-31']
        self.p15 = ['-08-01', '-08-15']
        self.p16 = ['-08-16', '-08-31']
        self.p17 = ['-09-01', '-09-15']
        self.p18 = ['-09-16', '-09-30']
        self.p19 = ['-10-01', '-10-15']
        self.p20 = ['-10-16', '-10-31']
        self.p21 = ['-11-01', '-11-15']
        self.p22 = ['-11-16', '-11-30']
        self.p23 = ['-12-01', '-12-15']
        self.p24 = ['-12-16', '-12-31']

        # Here we define one dict for each season, in order to have the choice to choose the selected statistic
        
        self.dwinter = {}
        self.dspring = {}
        self.dsummer = {}
        self.dautumn = {}

        self.djanuary = {}
        self.dfebruary ={}
        self.dmarch = {}
        self.dapril = {}
        self.dmay = {}
        self.djune = {}
        self.djuly = {}
        self.daugust = {}
        self.dseptember = {}
        self.doctober = {}
        self.dnovember = {}
        self.ddecember = {}

        self.dp1 = {}
        self.dp2 = {}
        self.dp3 = {}
        self.dp4 = {}
        self.dp5 = {}
        self.dp6 = {}
        self.dp7 = {}
        self.dp8 = {}
        self.dp9 = {}
        self.dp10 = {}
        self.dp11 = {}
        self.dp12 = {}
        self.dp13 = {}
        self.dp14 = {}
        self.dp15 = {}
        self.dp16 = {}
        self.dp17 = {}
        self.dp18 = {}
        self.dp19 = {}
        self.dp20 = {}
        self.dp21 = {}
        self.dp22 = {}
        self.dp23 = {}
        self.dp24 = {}

        # Here we defined the collections to generate the ndvi
        # Just need Red and NIR bandas, so avoid and rename them to apply the ndvi formula... 
        # ...to all sensors in the same way
        # Get Landsat surface reflectance collections for OLI, ETM+ and TM sensors.
        
        LC09col = ee.ImageCollection("LANDSAT/LC09/C02/T1_L2").filterBounds(self.roi) 
        LC08col = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2").filterBounds(self.roi) 
        LE07col = ee.ImageCollection("LANDSAT/LE07/C02/T1_L2").filterBounds(self.roi) 
        LT05col = ee.ImageCollection("LANDSAT/LT05/C02/T1_L2").filterBounds(self.roi) 
        LT04col = ee.ImageCollection("LANDSAT/LT04/C02/T1_L2").filterBounds(self.roi) 
        
        OLI = LC09col.merge(LC08col)
        ETM = LE07col.merge(LT05col).merge(LT04col)
        OLI_ = OLI.map(scale_OLI) 
        ETM_ = ETM.map(scale_ETM)
        Landsat = OLI_.merge(ETM_)
        
        # Notice that S2 is not in Surface Reflectance but in TOA, this because otherwise only
        # had data starting in 2017. Using BOA we have 2 extra years, starting at 2015
        # We use the wide 8 band instead of the 8A narrower badn, that is also more similar to landsat NIR
        # But this way we have NDVI at 10 m instead of 20. And we are using TOA instead of SR so, who cares?
        
        #"COPERNICUS/S2_SR_HARMONIZED"
        S2col = ee.ImageCollection("COPERNICUS/S2_HARMONIZED").select(['B2', 'B3', 'B4', 'B8', 'B11', 'B12'], 
                                                                         ['Blue', 'Green', 'Red', 'Nir', 'Swir1', 'Swir2']).filterBounds(self.roi)
        
        # Get MODIS MOD09Q1.006 Terra Surface Reflectance 8-Day Global 250m
        # This is going to be the coarse resolution we will add (250 m) but it cold be a good choice
        # for large areas. And we good pretty good data from 2000 until present
        MOD09Q1 = ee.ImageCollection("MODIS/061/MOD09A1").select(['sur_refl_b03', 'sur_refl_b04', 'sur_refl_b01', 'sur_refl_b02', 'sur_refl_b06', 'sur_refl_b07'], 
                                                                 ['Blue', 'Green', 'Red', 'Nir', 'Swir1', 'Swir2']).filterBounds(self.roi)
        #MODIS/006/MOD09Q1 old version
        # Let's try to add Sentinel 1 to have some S1 data analysis capabilities
        s1 = ee.ImageCollection('COPERNICUS/S1_GRD').filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH')).filter(ee.Filter.eq('instrumentMode', 'IW'))
        s1Ascending = s1.filter(ee.Filter.eq('orbitProperties_pass', 'ASCENDING'))
        s1Descending = s1.filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))
        s1S1 = s1Ascending.select('VH').merge(s1Descending.select('VH')).filterBounds(self.roi)

        
        
        # Set the collection that will be used
        if self.sat == 'S2':
            self.ndvi_col = S2col
        elif self.sat == 'Landsat':
            self.ndvi_col = Landsat
        elif self.sat == 'MODIS':
            self.ndvi_col = MOD09Q1
        elif self.sat == 'S1':
            self.ndvi_col = s1S1
        else: 
            print('Not a valid satellite')
            pass
    

    #Indexes
    #Someday it would be good just use Awesome spectral indexes library
    def get_ndvi(self, image):

        '''Here we apply the NDVI calculation'''   

        return image.normalizedDifference(['Nir', 'Red'])
    
    def get_ndwi(self, image):

        '''Here we apply the NDWI calculation'''   

        return image.normalizedDifference(['Green', 'Nir'])
    
    def get_evi(self, image):
    
        '''Here we apply the EVI calculation'''   

        # Compute the EVI using an expression.
        return image.expression(
            '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))', {
            'NIR': image.select('Nir'),
            'RED': image.select('Red'),
            'BLUE': image.select('Blue')}).rename(['nd'])
    
    def get_savi(self, image, L=0.428):
        
        '''Here we apply the SAVI calculation'''   

        # Compute the SAVI using an expression.
        return image.expression(
            '((NIR - RED) / (NIR + RED + L) * (1 + L))', {
            'NIR': image.select('Nir'),
            'RED': image.select('Red'),
            'L': L}).rename(['nd'])
    
    def get_aweinsh(self, image):
        
        '''Here we apply the SAVI calculation'''   

        # Compute the SAVI using an expression.
        return image.expression(
            '4.0 * (GREEN - SWIR1) - 0.25 * NIR + 2.75 * SWIR2', {
            'NIR': image.select('Nir'),
            'GREEN': image.select('Green'),
            'SWIR1':image.select('Swir1'),
            'SWIR2':image.select('Swir2')}).rename(['nd'])
    
    def get_awei(self, image):
        
        '''Here we apply the SAVI calculation'''   

        # Compute the SAVI using an expression.
        return image.expression(
            ('BLUE + 2.5 * GREEN - 1.5 * (NIR + SWIR1) - 0.25 * SWIR2'), {
            'NIR': image.select('Nir'),
            'BLUE': image.select('Blue'),
            'GREEN': image.select('Green'),
            'SWIR1':image.select('Swir1'),
            'SWIR2':image.select('Swir2')}).rename(['nd'])

    def get_gndvi(self, image):
        
        '''Here we apply the SAVI calculation'''   

        # Compute the GNDVI using an expression.
        return image.normalizedDifference(['Nir', 'Green'])
    
    
    def get_avi(self, image, L=0.428):
        
        '''Here we apply the SAVI calculation'''   

        # Compute the SAVI using an expression.
        return image.expression(
            '(NIR * (1.0 - RED) * (NIR - RED)) ** (1/3)', {
            'NIR': image.select('Nir'),
            'RED': image.select('Red')}).rename(['nd'])

    def get_nbri(self, image):
        
        '''Here we apply the SAVI calculation'''   

        # Compute the EVI using an expression.
        return image.normalizedDifference(['Nir', 'Swir2'])
    
    def get_ndsi(self, image):
        
        '''Here we apply the SAVI calculation'''   

        # Compute the EVI using an expression.
        return image.normalizedDifference(['Green', 'Swir1'])

    def get_ndmi(self, image):
        
        '''Here we apply the NDMI calculation'''   

        # Compute the EVI using an expression.
        return image.normalizedDifference(['Nir', 'Swir1'])

    
    # Here we define the functions to get the seasonal periods of time

    def get_winter(self, y):
        
        '''Here comes the funny thing. We generate the winter image for each year'''

        init, ends = str(y) + self.winter[0], str(y) + self.winter[1]

        if self.sat != 'S1':

            self.dwinter['max'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).max()
            self.dwinter['median'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).median()
            self.dwinter['mean'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).mean()
            self.dwinter['perc_90'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([90]))
            self.dwinter['perc_95'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([95]))

        else:

            self.dwinter['max'] = self.ndvi_col.filterDate(init, ends).max()
            self.dwinter['median'] = self.ndvi_col.filterDate(init, ends).median()
            self.dwinter['mean'] = self.ndvi_col.filterDate(init, ends).mean()
            self.dwinter['perc_90'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([90]))
            self.dwinter['perc_95'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([95]))

        return self.dwinter[self.key]
    
    def get_spring(self, y):
        
        '''Here comes the funny thing. We generate the spring image for each year'''

        init, ends = str(y) + self.spring[0], str(y) + self.spring[1]
        
        if self.sat != 'S1':

            self.dspring['max'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).max()
            self.dspring['median'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).median()
            self.dspring['mean'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).mean()
            self.dspring['perc_90'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([90]))
            self.dspring['perc_95'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([95]))

        else:

            self.dspring['max'] = self.ndvi_col.filterDate(init, ends).max()
            self.dspring['median'] = self.ndvi_col.filterDate(init, ends).median()
            self.dspring['mean'] = self.ndvi_col.filterDate(init, ends).mean()
            self.dspring['perc_90'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([90]))
            self.dspring['perc_95'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([95]))

        return self.dspring[self.key]


    def get_summer(self, y):
        
        '''Here comes the funny thing. We generate the summer image for each year'''

        init, ends = str(y) + self.summer[0], str(y) + self.summer[1]
        
        if self.sat != 'S1':

            self.dsummer['max'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).max()
            self.dsummer['median'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).median()
            self.dsummer['mean'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).mean()
            self.dsummer['perc_90'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([90]))
            self.dsummer['perc_95'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([95]))

        else:

            self.dsummer['max'] = self.ndvi_col.filterDate(init, ends).max()
            self.dsummer['median'] = self.ndvi_col.filterDate(init, ends).median()
            self.dsummer['mean'] = self.ndvi_col.filterDate(init, ends).mean()
            self.dsummer['perc_90'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([90]))
            self.dsummer['perc_95'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([95]))

        return self.dsummer[self.key]
        
    def get_autumn(self, y):
                
        '''Here comes the funny thing. We generate the autumn image for each year'''

        init, ends = str(y) + self.autumn[0], str(y) + self.autumn[1]

        if self.sat != 'S1':

            self.dautumn['max'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).max()
            self.dautumn['median'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).median()
            self.dautumn['mean'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).mean()
            self.dautumn['perc_90'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([90]))
            self.dautumn['perc_95'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([95]))

        else:

            self.dautumn['max'] = self.ndvi_col.filterDate(init, ends).max()
            self.dautumn['median'] = self.ndvi_col.filterDate(init, ends).median()
            self.dautumn['mean'] = self.ndvi_col.filterDate(init, ends).mean()
            self.dautumn['perc_90'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([90]))
            self.dautumn['perc_95'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([95]))

        return self.dautumn[self.key]
    
    def get_january(self, y):
                
        '''Here comes the funny thing. We generate the winter image for each year'''

        init, ends = str(y) + self.january[0], str(y) + self.january[1]

        if self.sat != 'S1':

            self.djanuary['max'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).max()
            self.djanuary['median'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).median()
            self.djanuary['mean'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).mean()
            self.djanuary['perc_90'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([90]))
            self.djanuary['perc_95'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([95]))

        else:

            self.djanuary['max'] = self.ndvi_col.filterDate(init, ends).max()
            self.djanuary['median'] = self.ndvi_col.filterDate(init, ends).median()
            self.djanuary['mean'] = self.ndvi_col.filterDate(init, ends).mean()
            self.djanuary['perc_90'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([90]))
            self.djanuary['perc_95'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([95]))

        return self.djanuary[self.key]
                
    def get_february(self, y):
        
        '''Here comes the funny thing. We generate the spring image for each year'''

        init, ends = str(y) + self.february[0], str(y) + self.february[1]
        d = {'ndvi': self.get_ndvi, 'ndwi': self.get_ndwi}

        if self.sat != 'S1':

            self.dfebruary['max'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).max()
            self.dfebruary['median'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).median()
            self.dfebruary['mean'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).mean()
            self.dfebruary['perc_90'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([90]))
            self.dfebruary['perc_95'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([95]))

        else:

            self.dfebruary['max'] = self.ndvi_col.filterDate(init, ends).max()
            self.dfebruary['median'] = self.ndvi_col.filterDate(init, ends).median()
            self.dfebruary['mean'] = self.ndvi_col.filterDate(init, ends).mean()
            self.dfebruary['perc_90'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([90]))
            self.dfebruary['perc_95'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([95]))

        return self.dfebruary[self.key]


    def get_march(self, y):
        
        '''Here comes the funny thing. We generate the summer image for each year'''

        init, ends = str(y) + self.march[0], str(y) + self.march[1]

        if self.sat != 'S1':

            self.dmarch['max'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).max()
            self.dmarch['median'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).median()
            self.dmarch['mean'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).mean()
            self.dmarch['perc_90'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([90]))
            self.dmarch['perc_95'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([95]))

        else:

            self.dmarch['max'] = self.ndvi_col.filterDate(init, ends).max()
            self.dmarch['median'] = self.ndvi_col.filterDate(init, ends).median()
            self.dmarch['mean'] = self.ndvi_col.filterDate(init, ends).mean()
            self.dmarch['perc_90'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([90]))
            self.dmarch['perc_95'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([95]))

        return self.dmarch[self.key]

        
    def get_april(self, y):
        
        '''Here comes the funny thing. We generate the autumn image for each year'''

        init, ends = str(y) + self.april[0], str(y) + self.april[1]

        if self.sat != 'S1':

            self.dapril['max'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).max()
            self.dapril['median'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).median()
            self.dapril['mean'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).mean()
            self.dapril['perc_90'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([90]))
            self.dapril['perc_95'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([95]))

        else:

            self.dapril['max'] = self.ndvi_col.filterDate(init, ends).max()
            self.dapril['median'] = self.ndvi_col.filterDate(init, ends).median()
            self.dapril['mean'] = self.ndvi_col.filterDate(init, ends).mean()
            self.dapril['perc_90'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([90]))
            self.dapril['perc_95'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([95]))

        return self.dapril[self.key]


    def get_may(self, y):
        
        '''Here comes the funny thing. We generate the winter image for each year'''

        init, ends = str(y) + self.may[0], str(y) + self.may[1]

        if self.sat != 'S1':

            self.dmay['max'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).max()
            self.dmay['median'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).median()
            self.dmay['mean'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).mean()
            self.dmay['perc_90'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([90]))
            self.dmay['perc_95'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([95]))

        else:

            self.dmay['max'] = self.ndvi_col.filterDate(init, ends).max()
            self.dmay['median'] = self.ndvi_col.filterDate(init, ends).median()
            self.dmay['mean'] = self.ndvi_col.filterDate(init, ends).mean()
            self.dmay['perc_90'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([90]))
            self.dmay['perc_95'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([95]))

        return self.dmay[self.key]
        
    def get_june(self, y):
        
        '''Here comes the funny thing. We generate the spring image for each year'''

        init, ends = str(y) + self.june[0], str(y) + self.june[1]

        if self.sat != 'S1':

            self.djune['max'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).max()
            self.djune['median'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).median()
            self.djune['mean'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).mean()
            self.djune['perc_90'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([90]))
            self.djune['perc_95'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([95]))

        else:

            self.djune['max'] = self.ndvi_col.filterDate(init, ends).max()
            self.djune['median'] = self.ndvi_col.filterDate(init, ends).median()
            self.djune['mean'] = self.ndvi_col.filterDate(init, ends).mean()
            self.djune['perc_90'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([90]))
            self.djune['perc_95'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([95]))

        return self.djune[self.key]


    def get_july(self, y):
        
        '''Here comes the funny thing. We generate the summer image for each year'''

        init, ends = str(y) + self.july[0], str(y) + self.july[1]

        if self.sat != 'S1':

            self.djuly['max'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).max()
            self.djuly['median'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).median()
            self.djuly['mean'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).mean()
            self.djuly['perc_90'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([90]))
            self.djuly['perc_95'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([95]))

        else:

            self.djuly['max'] = self.ndvi_col.filterDate(init, ends).max()
            self.djuly['median'] = self.ndvi_col.filterDate(init, ends).median()
            self.djuly['mean'] = self.ndvi_col.filterDate(init, ends).mean()
            self.djuly['perc_90'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([90]))
            self.djuly['perc_95'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([95]))

        return self.djuly[self.key]

        
    def get_august(self, y):
        
        '''Here comes the funny thing. We generate the autumn image for each year'''

        init, ends = str(y) + self.august[0], str(y) + self.august[1]

        if self.sat != 'S1':

            self.daugust['max'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).max()
            self.daugust['median'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).median()
            self.daugust['mean'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).mean()
            self.daugust['perc_90'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([90]))
            self.daugust['perc_95'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([95]))

        else:

            self.daugust['max'] = self.ndvi_col.filterDate(init, ends).max()
            self.daugust['median'] = self.ndvi_col.filterDate(init, ends).median()
            self.daugust['mean'] = self.ndvi_col.filterDate(init, ends).mean()
            self.daugust['perc_90'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([90]))
            self.daugust['perc_95'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([95]))

        return self.daugust[self.key]


    def get_september(self, y):
            
            '''Here comes the funny thing. We generate the winter image for each year'''

            init, ends = str(y) + self.september[0], str(y) + self.september[1]

            if self.sat != 'S1':

                self.dseptember['max'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).max()
                self.dseptember['median'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).median()
                self.dseptember['mean'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).mean()
                self.dseptember['perc_90'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([90]))
                self.dseptember['perc_95'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([95]))

            else:

                self.dseptember['max'] = self.ndvi_col.filterDate(init, ends).max()
                self.dseptember['median'] = self.ndvi_col.filterDate(init, ends).median()
                self.dseptember['mean'] = self.ndvi_col.filterDate(init, ends).mean()
                self.dseptember['perc_90'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([90]))
                self.dseptember['perc_95'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([95]))

            return self.dseptember[self.key]


    def get_october(self, y):
        
        '''Here comes the funny thing. We generate the spring image for each year'''

        init, ends = str(y) + self.october[0], str(y) + self.october[1]

        if self.sat != 'S1':

            self.doctober['max'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).max()
            self.doctober['median'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).median()
            self.doctober['mean'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).mean()
            self.doctober['perc_90'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([90]))
            self.doctober['perc_95'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([95]))

        else:

            self.doctober['max'] = self.ndvi_col.filterDate(init, ends).max()
            self.doctober['median'] = self.ndvi_col.filterDate(init, ends).median()
            self.doctober['mean'] = self.ndvi_col.filterDate(init, ends).mean()
            self.doctober['perc_90'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([90]))
            self.doctober['perc_95'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([95]))

        return self.doctober[self.key]


    def get_november(self, y):
        
        '''Here comes the funny thing. We generate the summer image for each year'''

        init, ends = str(y) + self.november[0], str(y) + self.november[1]

        if self.sat != 'S1':

            self.dnovember['max'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).max()
            self.dnovember['median'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).median()
            self.dnovember['mean'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).mean()
            self.dnovember['perc_90'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([90]))
            self.dnovember['perc_95'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([95]))

        else:

            self.dnovember['max'] = self.ndvi_col.filterDate(init, ends).max()
            self.dnovember['median'] = self.ndvi_col.filterDate(init, ends).median()
            self.dnovember['mean'] = self.ndvi_col.filterDate(init, ends).mean()
            self.dnovember['perc_90'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([90]))
            self.dnovember['perc_95'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([95]))

        return self.dnovember[self.key]

        
    def get_december(self, y):
        
        '''Here comes the funny thing. We generate the autumn image for each year'''

        init, ends = str(y) + self.december[0], str(y) + self.december[1]

        if self.sat != 'S1':

            self.ddecember['max'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).max()
            self.ddecember['median'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).median()
            self.ddecember['mean'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).mean()
            self.ddecember['perc_90'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([90]))
            self.ddecember['perc_95'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([95]))

        else:

            self.ddecember['max'] = self.ndvi_col.filterDate(init, ends).max()
            self.ddecember['median'] = self.ndvi_col.filterDate(init, ends).median()
            self.ddecember['mean'] = self.ndvi_col.filterDate(init, ends).mean()
            self.ddecember['perc_90'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([90]))
            self.ddecember['perc_95'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([95]))

        return self.ddecember[self.key]
    

    def get_p1(self, y):
                
        '''Here comes the funny thing. We generate the winter image for each year'''

        init, ends = str(y) + self.p1[0], str(y) + self.p1[1]

        if self.sat != 'S1':

            self.dp1['max'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).max()
            self.dp1['median'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).median()
            self.dp1['mean'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).mean()
            self.dp1['perc_90'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([90]))
            self.dp1['perc_95'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([95]))

        else:

            self.dp1['max'] = self.ndvi_col.filterDate(init, ends).max()
            self.dp1['median'] = self.ndvi_col.filterDate(init, ends).median()
            self.dp1['mean'] = self.ndvi_col.filterDate(init, ends).mean()
            self.dp1['perc_90'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([90]))
            self.dp1['perc_95'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([95]))

        return self.dp1[self.key]
        
    def get_p2(self, y):
        
        '''Here comes the funny thing. We generate the spring image for each year'''

        init, ends = str(y) + self.p2[0], str(y) + self.p2[1]

        if self.sat != 'S1':

            self.dp2['max'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).max()
            self.dp2['median'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).median()
            self.dp2['mean'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).mean()
            self.dp2['perc_90'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([90]))
            self.dp2['perc_95'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([95]))

        else:

            self.dp2['max'] = self.ndvi_col.filterDate(init, ends).max()
            self.dp2['median'] = self.ndvi_col.filterDate(init, ends).median()
            self.dp2['mean'] = self.ndvi_col.filterDate(init, ends).mean()
            self.dp2['perc_90'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([90]))
            self.dp2['perc_95'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([95]))

        return self.dp2[self.key]


    def get_p3(self, y):
        
        '''Here comes the funny thing. We generate the summer image for each year'''

        init, ends = str(y) + self.p3[0], str(y) + self.p3[1]

        if self.sat != 'S1':

            self.dp3['max'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).max()
            self.dp3['median'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).median()
            self.dp3['mean'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).mean()
            self.dp3['perc_90'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([90]))
            self.dp3['perc_95'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([95]))

        else:

            self.dp3['max'] = self.ndvi_col.filterDate(init, ends).max()
            self.dp3['median'] = self.ndvi_col.filterDate(init, ends).median()
            self.dp3['mean'] = self.ndvi_col.filterDate(init, ends).mean()
            self.dp3['perc_90'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([90]))
            self.dp3['perc_95'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([95]))

        return self.dp3[self.key]

        
    def get_p4(self, y):
        
        '''Here comes the funny thing. We generate the autumn image for each year'''

        init, ends = str(y) + self.p4[0], str(y) + self.p4[1]

        if self.sat != 'S1':

            self.dp4['max'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).max()
            self.dp4['median'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).median()
            self.dp4['mean'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).mean()
            self.dp4['perc_90'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([90]))
            self.dp4['perc_95'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([95]))

        else:

            self.dp4['max'] = self.ndvi_col.filterDate(init, ends).max()
            self.dp4['median'] = self.ndvi_col.filterDate(init, ends).median()
            self.dp4['mean'] = self.ndvi_col.filterDate(init, ends).mean()
            self.dp4['perc_90'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([90]))
            self.dp4['perc_95'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([95]))

        return self.dp4[self.key]
        
        
    def get_p5(self, y):
        
        '''Here comes the funny thing. We generate the winter image for each year'''

        init, ends = str(y) + self.p5[0], str(y) + self.p5[1]

        if self.sat != 'S1':

            self.dp5['max'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).max()
            self.dp5['median'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).median()
            self.dp5['mean'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).mean()
            self.dp5['perc_90'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([90]))
            self.dp5['perc_95'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([95]))

        else:

            self.dp5['max'] = self.ndvi_col.filterDate(init, ends).max()
            self.dp5['median'] = self.ndvi_col.filterDate(init, ends).median()
            self.dp5['mean'] = self.ndvi_col.filterDate(init, ends).mean()
            self.dp5['perc_90'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([90]))
            self.dp5['perc_95'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([95]))

        return self.dp5[self.key]
        
    def get_p6(self, y):
        
        '''Here comes the funny thing. We generate the spring image for each year'''

        init, ends = str(y) + self.p6[0], str(y) + self.p6[1]

        if self.sat != 'S1':

            self.dp6['max'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).max()
            self.dp6['median'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).median()
            self.dp6['mean'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).mean()
            self.dp6['perc_90'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([90]))
            self.dp6['perc_95'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([95]))

        else:

            self.dp6['max'] = self.ndvi_col.filterDate(init, ends).max()
            self.dp6['median'] = self.ndvi_col.filterDate(init, ends).median()
            self.dp6['mean'] = self.ndvi_col.filterDate(init, ends).mean()
            self.dp6['perc_90'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([90]))
            self.dp6['perc_95'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([95]))

        return self.dp6[self.key]


    def get_p7(self, y):
        
        '''Here comes the funny thing. We generate the summer image for each year'''

        init, ends = str(y) + self.p7[0], str(y) + self.p7[1]

        if self.sat != 'S1':

            self.dp7['max'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).max()
            self.dp7['median'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).median()
            self.dp7['mean'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).mean()
            self.dp7['perc_90'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([90]))
            self.dp7['perc_95'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([95]))

        else:

            self.dp7['max'] = self.ndvi_col.filterDate(init, ends).max()
            self.dp7['median'] = self.ndvi_col.filterDate(init, ends).median()
            self.dp7['mean'] = self.ndvi_col.filterDate(init, ends).mean()
            self.dp7['perc_90'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([90]))
            self.dp7['perc_95'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([95]))

        return self.dp7[self.key]

        
    def get_p8(self, y):
        
        '''Here comes the funny thing. We generate the autumn image for each year'''

        init, ends = str(y) + self.p8[0], str(y) + self.p8[1]

        if self.sat != 'S1':

            self.dp8['max'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).max()
            self.dp8['median'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).median()
            self.dp8['mean'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).mean()
            self.dp8['perc_90'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([90]))
            self.dp8['perc_95'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([95]))

        else:

            self.dp8['max'] = self.ndvi_col.filterDate(init, ends).max()
            self.dp8['median'] = self.ndvi_col.filterDate(init, ends).median()
            self.dp8['mean'] = self.ndvi_col.filterDate(init, ends).mean()
            self.dp8['perc_90'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([90]))
            self.dp8['perc_95'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([95]))

        return self.dp8[self.key]
        
        
    def get_p9(self, y):
            
        '''Here comes the funny thing. We generate the winter image for each year'''

        init, ends = str(y) + self.p9[0], str(y) + self.p9[1]

        if self.sat != 'S1':

            self.dp9['max'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).max()
            self.dp9['median'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).median()
            self.dp9['mean'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).mean()
            self.dp9['perc_90'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([90]))
            self.dp9['perc_95'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([95]))

        else:

            self.dp9['max'] = self.ndvi_col.filterDate(init, ends).max()
            self.dp9['median'] = self.ndvi_col.filterDate(init, ends).median()
            self.dp9['mean'] = self.ndvi_col.filterDate(init, ends).mean()
            self.dp9['perc_90'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([90]))
            self.dp9['perc_95'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([95]))

        return self.dp9[self.key]

        
    def get_p10(self, y):
        
        '''Here comes the funny thing. We generate the spring image for each year'''

        init, ends = str(y) + self.p10[0], str(y) + self.p10[1]

        if self.sat != 'S1':

            self.dp10['max'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).max()
            self.dp10['median'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).median()
            self.dp10['mean'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).mean()
            self.dp10['perc_90'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([90]))
            self.dp10['perc_95'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([95]))

        else:

            self.dp10['max'] = self.ndvi_col.filterDate(init, ends).max()
            self.dp10['median'] = self.ndvi_col.filterDate(init, ends).median()
            self.dp10['mean'] = self.ndvi_col.filterDate(init, ends).mean()
            self.dp10['perc_90'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([90]))
            self.dp10['perc_95'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([95]))

        return self.dp10[self.key]


    def get_p11(self, y):
        
        '''Here comes the funny thing. We generate the summer image for each year'''

        init, ends = str(y) + self.p11[0], str(y) + self.p11[1]

        if self.sat != 'S1':

            self.dp11['max'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).max()
            self.dp11['median'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).median()
            self.dp11['mean'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).mean()
            self.dp11['perc_90'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([90]))
            self.dp11['perc_95'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([95]))

        else:

            self.dp11['max'] = self.ndvi_col.filterDate(init, ends).max()
            self.dp11['median'] = self.ndvi_col.filterDate(init, ends).median()
            self.dp11['mean'] = self.ndvi_col.filterDate(init, ends).mean()
            self.dp11['perc_90'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([90]))
            self.dp11['perc_95'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([95]))

        return self.dp11[self.key]

        
    def get_p12(self, y):
        
        '''Here comes the funny thing. We generate the autumn image for each year'''

        init, ends = str(y) + self.p12[0], str(y) + self.p12[1]

        if self.sat != 'S1':

            self.dp12['max'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).max()
            self.dp12['median'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).median()
            self.dp12['mean'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).mean()
            self.dp12['perc_90'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([90]))
            self.dp12['perc_95'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([95]))

        else:

            self.dp12['max'] = self.ndvi_col.filterDate(init, ends).max()
            self.dp12['median'] = self.ndvi_col.filterDate(init, ends).median()
            self.dp12['mean'] = self.ndvi_col.filterDate(init, ends).mean()
            self.dp12['perc_90'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([90]))
            self.dp12['perc_95'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([95]))

        return self.dp12[self.key]
        

    def get_p13(self, y):
            
        '''Here comes the funny thing. We generate the winter image for each year'''

        init, ends = str(y) + self.p13[0], str(y) + self.p13[1]

        if self.sat != 'S1':

            self.dp13['max'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).max()
            self.dp13['median'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).median()
            self.dp13['mean'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).mean()
            self.dp13['perc_90'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([90]))
            self.dp13['perc_95'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([95]))

        else:

            self.dp13['max'] = self.ndvi_col.filterDate(init, ends).max()
            self.dp13['median'] = self.ndvi_col.filterDate(init, ends).median()
            self.dp13['mean'] = self.ndvi_col.filterDate(init, ends).mean()
            self.dp13['perc_90'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([90]))
            self.dp13['perc_95'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([95]))

        return self.dp13[self.key]
        
    def get_p14(self, y):
        
        '''Here comes the funny thing. We generate the spring image for each year'''

        init, ends = str(y) + self.p14[0], str(y) + self.p14[1]

        if self.sat != 'S1':

            self.dp14['max'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).max()
            self.dp14['median'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).median()
            self.dp14['mean'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).mean()
            self.dp14['perc_90'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([90]))
            self.dp14['perc_95'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([95]))

        else:

            self.dp14['max'] = self.ndvi_col.filterDate(init, ends).max()
            self.dp14['median'] = self.ndvi_col.filterDate(init, ends).median()
            self.dp14['mean'] = self.ndvi_col.filterDate(init, ends).mean()
            self.dp14['perc_90'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([90]))
            self.dp14['perc_95'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([95]))

        return self.dp14[self.key]


    def get_p15(self, y):
        
        '''Here comes the funny thing. We generate the summer image for each year'''

        init, ends = str(y) + self.p15[0], str(y) + self.p15[1]

        if self.sat != 'S1':

            self.dp3['max'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).max()
            self.dp3['median'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).median()
            self.dp3['mean'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).mean()
            self.dp3['perc_90'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([90]))
            self.dp3['perc_95'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([95]))

        else:

            self.dp3['max'] = self.ndvi_col.filterDate(init, ends).max()
            self.dp3['median'] = self.ndvi_col.filterDate(init, ends).median()
            self.dp3['mean'] = self.ndvi_col.filterDate(init, ends).mean()
            self.dp3['perc_90'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([90]))
            self.dp3['perc_95'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([95]))

        return self.dp3[self.key]
            
                    
    def get_p16(self, y):
        
        '''Here comes the funny thing. We generate the autumn image for each year'''

        init, ends = str(y) + self.p16[0], str(y) + self.p16[1]

        if self.sat != 'S1':

            self.dp16['max'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).max()
            self.dp16['median'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).median()
            self.dp16['mean'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).mean()
            self.dp16['perc_90'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([90]))
            self.dp16['perc_95'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([95]))

        else:

            self.dp16['max'] = self.ndvi_col.filterDate(init, ends).max()
            self.dp16['median'] = self.ndvi_col.filterDate(init, ends).median()
            self.dp16['mean'] = self.ndvi_col.filterDate(init, ends).mean()
            self.dp16['perc_90'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([90]))
            self.dp16['perc_95'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([95]))

        return self.dp16[self.key]


    def get_p17(self, y):
        
        '''Here comes the funny thing. We generate the winter image for each year'''

        init, ends = str(y) + self.p17[0], str(y) + self.p17[1]

        if self.sat != 'S1':

            self.dp17['max'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).max()
            self.dp17['median'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).median()
            self.dp17['mean'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).mean()
            self.dp17['perc_90'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([90]))
            self.dp17['perc_95'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([95]))

        else:

            self.dp17['max'] = self.ndvi_col.filterDate(init, ends).max()
            self.dp17['median'] = self.ndvi_col.filterDate(init, ends).median()
            self.dp17['mean'] = self.ndvi_col.filterDate(init, ends).mean()
            self.dp17['perc_90'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([90]))
            self.dp17['perc_95'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([95]))

        return self.dp17[self.key]
                    
    def get_p18(self, y):
        
        '''Here comes the funny thing. We generate the spring image for each year'''

        init, ends = str(y) + self.p18[0], str(y) + self.p18[1]

        if self.sat != 'S1':

            self.dp18['max'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).max()
            self.dp18['median'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).median()
            self.dp18['mean'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).mean()
            self.dp18['perc_90'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([90]))
            self.dp18['perc_95'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([95]))

        else:

            self.dp18['max'] = self.ndvi_col.filterDate(init, ends).max()
            self.dp18['median'] = self.ndvi_col.filterDate(init, ends).median()
            self.dp18['mean'] = self.ndvi_col.filterDate(init, ends).mean()
            self.dp18['perc_90'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([90]))
            self.dp18['perc_95'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([95]))

        return self.dp18[self.key]


    def get_p19(self, y):
        
        '''Here comes the funny thing. We generate the summer image for each year'''

        init, ends = str(y) + self.p19[0], str(y) + self.p19[1]

        if self.sat != 'S1':

            self.dp19['max'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).max()
            self.dp19['median'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).median()
            self.dp19['mean'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).mean()
            self.dp19['perc_90'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([90]))
            self.dp19['perc_95'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([95]))

        else:

            self.dp19['max'] = self.ndvi_col.filterDate(init, ends).max()
            self.dp19['median'] = self.ndvi_col.filterDate(init, ends).median()
            self.dp19['mean'] = self.ndvi_col.filterDate(init, ends).mean()
            self.dp19['perc_90'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([90]))
            self.dp19['perc_95'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([95]))

        return self.dp19[self.key]

        
    def get_p20(self, y):
        
        '''Here comes the funny thing. We generate the autumn image for each year'''

        init, ends = str(y) + self.p20[0], str(y) + self.p20[1]

        if self.sat != 'S1':

            self.dp20['max'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).max()
            self.dp20['median'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).median()
            self.dp20['mean'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).mean()
            self.dp20['perc_90'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([90]))
            self.dp20['perc_95'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([95]))

        else:

            self.dp20['max'] = self.ndvi_col.filterDate(init, ends).max()
            self.dp20['median'] = self.ndvi_col.filterDate(init, ends).median()
            self.dp20['mean'] = self.ndvi_col.filterDate(init, ends).mean()
            self.dp20['perc_90'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([90]))
            self.dp20['perc_95'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([95]))

        return self.dp20[self.key]
            
            
    def get_p21(self, y):
            
        '''Here comes the funny thing. We generate the winter image for each year'''

        init, ends = str(y) + self.p21[0], str(y) + self.p21[1]

        if self.sat != 'S1':

            self.dp21['max'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).max()
            self.dp21['median'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).median()
            self.dp21['mean'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).mean()
            self.dp21['perc_90'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([90]))
            self.dp21['perc_95'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([95]))

        else:

            self.dp21['max'] = self.ndvi_col.filterDate(init, ends).max()
            self.dp21['median'] = self.ndvi_col.filterDate(init, ends).median()
            self.dp21['mean'] = self.ndvi_col.filterDate(init, ends).mean()
            self.dp21['perc_90'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([90]))
            self.dp21['perc_95'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([95]))

        return self.dp21[self.key]

        
    def get_p22(self, y):
        
        '''Here comes the funny thing. We generate the spring image for each year'''

        init, ends = str(y) + self.p22[0], str(y) + self.p22[1]

        if self.sat != 'S1':

            self.dp22['max'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).max()
            self.dp22['median'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).median()
            self.dp22['mean'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).mean()
            self.dp22['perc_90'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([90]))
            self.dp22['perc_95'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([95]))

        else:

            self.dp22['max'] = self.ndvi_col.filterDate(init, ends).max()
            self.dp22['median'] = self.ndvi_col.filterDate(init, ends).median()
            self.dp22['mean'] = self.ndvi_col.filterDate(init, ends).mean()
            self.dp22['perc_90'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([90]))
            self.dp22['perc_95'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([95]))

        return self.dp22[self.key]


    def get_p23(self, y):
        
        '''Here comes the funny thing. We generate the summer image for each year'''

        init, ends = str(y) + self.p23[0], str(y) + self.p23[1]

        if self.sat != 'S1':

            self.dp23['max'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).max()
            self.dp23['median'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).median()
            self.dp23['mean'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).mean()
            self.dp23['perc_90'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([90]))
            self.dp23['perc_95'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([95]))

        else:

            self.dp23['max'] = self.ndvi_col.filterDate(init, ends).max()
            self.dp23['median'] = self.ndvi_col.filterDate(init, ends).median()
            self.dp23['mean'] = self.ndvi_col.filterDate(init, ends).mean()
            self.dp23['perc_90'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([90]))
            self.dp23['perc_95'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([95]))

        return self.dp23[self.key]
            
                    
    def get_p24(self, y):
        
        '''Here comes the funny thing. We generate the autumn image for each year'''

        init, ends = str(y) + self.p24[0], str(y) + self.p24[1]

        if self.sat != 'S1':

            self.dp24['max'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).max()
            self.dp24['median'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).median()
            self.dp24['mean'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).mean()
            self.dp24['perc_90'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([90]))
            self.dp24['perc_95'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([95]))

        else:

            self.dp24['max'] = self.ndvi_col.filterDate(init, ends).max()
            self.dp24['median'] = self.ndvi_col.filterDate(init, ends).median()
            self.dp24['mean'] = self.ndvi_col.filterDate(init, ends).mean()
            self.dp24['perc_90'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([90]))
            self.dp24['perc_95'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([95]))

        return self.dp24[self.key]

    


    
    #This should be done just once. But we will get into this later
    def get_year_composite(self):
        
        '''Return the composite ndvi for each year'''
        
        # Maybe this should do with .map instead of a loop
        if self.periods == 4:

            bands = ['nd', 'nd_1', 'nd_2', 'nd_3']
            sarbands = ['VH', 'VH_1', 'VH_2', 'VH_3']
            sarbandsp90 = ['VH_p90', 'VH_p90_1', 'VH_p90_2', 'VH_p90_3']
            sarbandsp95 = ['VH_p95', 'VH_p95_1', 'VH_p95_2', 'VH_p95_3']
            bandsp90 = ['nd_p90', 'nd_p90_1', 'nd_p90_2', 'nd_p90_3']
            bandsp95 = ['nd_p95', 'nd_p95_1', 'nd_p95_2', 'nd_p95_3']
            nbands = ['winter', 'spring', 'summer', 'autumn']

            for y in range(self.start_year, self.end_year):
                
                composite = ee.Image.cat(self.get_winter(y), self.get_spring(y), self.get_summer(y), self.get_autumn(y)).clip(self.roi)

                
                if self.sat != 'S1':
                    
                    if self.key not in ['perc_90', 'perc_95']:
                        
                        compositer = composite.select(bands, nbands)
                
                    elif self.key == 'perc_90':

                        compositer = composite.select(bandsp90, nbands)

                    elif self.key == 'perc_95':
                        
                        compositer = composite.select(bandsp95, nbands)

                    else:
                        
                        print('Uhm... it seems that you choose a bad statistic name')

                
                else: 
                    
                    if self.key not in ['perc_90', 'perc_95']:
                        
                        compositer = composite.select(sarbands, nbands)
                
                    elif self.key == 'perc_90':
        
                        compositer = composite.select(sarbandsp90, nbands)

                    elif self.key == 'perc_95':
                        
                        compositer = composite.select(sarbandsp95, nbands)

                    else:
                        
                        print('Uhm... it seems that you choose a bad statistic name')
                
                
                
                self.imagelist.append(compositer)
            
            
            ndvi_comp_coll = ee.ImageCollection.fromImages(self.imagelist)
        
            return ndvi_comp_coll
    
        elif self.periods == 12:
    
            bands = ['nd', 'nd_1', 'nd_2', 'nd_3', 'nd_4', 'nd_5', 'nd_6', 'nd_7', 'nd_8', 'nd_9', 'nd_10', 'nd_11']
            sarbands = ['VH', 'VH_1', 'VH_2', 'VH_3', 'VH_4', 'VH_5', 'VH_6', 'VH_7', 'VH_8', 'VH_9', 'VH_10', 'VH_11']
            sarbandsp90 = ['VH_p90', 'VH_p90_1', 'VH_p90_2', 'VH_p90_3', 'VH_p90_4', 'VH_p90_5', 'VH_p90_6', 'VH_p90_7', 'VH_p90_8', 'VH_p90_9', 'VH_p90_10', 'VH_p90_11']
            sarbandsp95 = ['VH_p95', 'VH_p95_1', 'VH_p95_2', 'VH_p95_3', 'VH_p95_4', 'VH_p95_5', 'VH_p95_6', 'VH_p95_7', 'VH_p95_8', 'VH_p95_9', 'VH_p95_10', 'VH_p95_11']
            bandsp90 = ['nd_p90', 'nd_p90_1', 'nd_p90_2', 'nd_p90_3', 'nd_p90_4', 'nd_p90_5', 'nd_p90_6', 'nd_p90_7', 'nd_p90_8', 'nd_p90_9', 'nd_p90_10', 'nd_p90_11']
            bandsp95 = ['nd_p95', 'nd_p95_1', 'nd_p95_2', 'nd_p95_3', 'nd_p95_4', 'nd_p95_5', 'nd_p95_6', 'nd_p95_7', 'nd_p95_8', 'nd_p95_9', 'nd_p95_10', 'nd_p95_11']
            nbands = ['january', 'february', 'march', 'april', 'may', 'june', 'july', 'august', 'september', 'october', 'november', 'december']

            for y in range(self.start_year, self.end_year):
                
                composite = ee.Image.cat(self.get_january(y), self.get_february(y), self.get_march(y), self.get_april(y), self.get_may(y), 
                                self.get_june(y), self.get_july(y), self.get_august(y), self.get_september(y), self.get_october(y),
                                self.get_november(y), self.get_december(y)).clip(self.roi)

                
                if self.sat != 'S1':
                    
                    if self.key not in ['perc_90', 'perc_95']:
                        
                        compositer = composite.select(bands, nbands)
                
                    elif self.key == 'perc_90':

                        compositer = composite.select(bandsp90, nbands)

                    elif self.key == 'perc_95':
                        
                        compositer = composite.select(bandsp95, nbands)

                    else:
                        
                        print('Uhm... it seems that you choose a bad statistic name')

                
                else: 
                    
                    if self.key not in ['perc_90', 'perc_95']:
                        
                        compositer = composite.select(sarbands, nbands)
                
                    elif self.key == 'perc_90':
        
                        compositer = composite.select(sarbandsp90, nbands)

                    elif self.key == 'perc_95':
                        
                        compositer = composite.select(sarbandsp95, nbands)

                    else:
                        
                        print('Uhm... it seems that you choose a bad statistic name')
                
                
                
                self.imagelist.append(compositer)
            
            
            ndvi_comp_coll = ee.ImageCollection.fromImages(self.imagelist)
        
            return ndvi_comp_coll
        
        elif self.periods == 24:

            bands = ['nd', 'nd_1', 'nd_2', 'nd_3', 'nd_4', 'nd_5', 'nd_6', 'nd_7', 'nd_8', 'nd_9', 'nd_10', 'nd_11', 'nd_12', 'nd_13', 'nd_14', 
                     'nd_15', 'nd_16', 'nd_17', 'nd_18', 'nd_19', 'nd_20', 'nd_21', 'nd_22', 'nd_23']
            sarbands = ['VH', 'VH_1', 'VH_2', 'VH_3', 'VH_4', 'VH_5', 'VH_6', 'VH_7', 'VH_8', 'VH_9', 'VH_10', 'VH_11', 'VH_12', 'VH_13', 'VH_14', 'VH_15', 'VH_16', 
                        'VH_17', 'VH_18', 'VH_19', 'VH_20', 'VH_21', 'VH_22', 'VH_23']
            sarbandsp90 = ['VH_p90', 'VH_p90_1', 'VH_p90_2', 'VH_p90_3', 'VH_p90_4', 'VH_p90_5', 'VH_p90_6', 'VH_p90_7', 'VH_p90_8', 'VH_p90_9', 'VH_p90_10', 'VH_p90_11', 
                        'VH_p90_12', 'VH_p90_13', 'VH_p90_14', 'VH_p90_15', 'VH_p90_16', 'VH_p90_17', 'VH_p90_18', 'VH_p90_19', 'VH_p90_20', 'VH_p90_21', 'VH_p90_22', 'VH_p90_23']
            sarbandsp95 = ['VH_p95', 'VH_p95_1', 'VH_p95_2', 'VH_p95_3', 'VH_p95_4', 'VH_p95_5', 'VH_p95_6', 'VH_p95_7', 'VH_p95_8', 'VH_p95_9', 'VH_p95_10', 'VH_p95_11', 
                        'VH_p95_12', 'VH_p95_13', 'VH_p95_14', 'VH_p95_15', 'VH_p95_16', 'VH_p95_17', 'VH_p95_18', 'VH_p95_19', 'VH_p95_20', 'VH_p95_21', 'VH_p95_22', 'VH_p95_23']
            bandsp90 = ['nd_p90', 'nd_p90_1', 'nd_p90_2', 'nd_p90_3', 'nd_p90_4', 'nd_p90_5', 'nd_p90_6', 'nd_p90_7', 'nd_p90_8', 'nd_p90_9', 'nd_p90_10', 'nd_p90_11',
                       'nd_p90_12', 'nd_p90_13', 'nd_p90_14', 'nd_p90_15', 'nd_p90_16', 'nd_p90_17', 'nd_p90_18', 'nd_p90_19', 'nd_p90_20', 'nd_p90_21', 'nd_p90_22', 'nd_p90_23' ]
            bandsp95 = ['nd_p95', 'nd_p95_1', 'nd_p95_2', 'nd_p95_3', 'nd_p95_4', 'nd_p95_5', 'nd_p95_6', 'nd_p95_7', 'nd_p95_8', 'nd_p95_9', 'nd_p95_10', 'nd_p95_11',
                       'nd_p95_12', 'nd_p95_13', 'nd_p95_14', 'nd_p95_15', 'nd_p95_16', 'nd_p95_17', 'nd_p95_18', 'nd_p95_19', 'nd_p95_20', 'nd_p95_21', 'nd_p95_22', 'nd_p95_23' ]
            nbands = ['p1', 'p2', 'p3', 'p4', 'p5', 'p6', 'p7', 'p8', 'p9', 'p10', 'p11', 'p12', 'p13', 'p14', 'p15', 'p16', 'p17', 'p18', 'p19', 'p20', 'p21', 'p22', 'p23', 'p24']


            for y in range(self.start_year, self.end_year):
                
                composite = ee.Image.cat(self.get_p1(y), self.get_p2(y), self.get_p3(y), self.get_p4(y), self.get_p5(y), 
                                self.get_p6(y), self.get_p7(y), self.get_p8(y), self.get_p9(y), self.get_p10(y),
                                self.get_p11(y), self.get_p12(y), self.get_p13(y), self.get_p14(y), self.get_p15(y),
                                self.get_p16(y), self.get_p17(y), self.get_p18(y), self.get_p19(y), self.get_p20(y),
                                self.get_p21(y), self.get_p22(y), self.get_p23(y), self.get_p24(y)).clip(self.roi)

                
                if self.sat != 'S1':
                    
                    if self.key not in ['perc_90', 'perc_95']:
                        
                        compositer = composite.select(bands, nbands)
                
                    elif self.key == 'perc_90':

                        compositer = composite.select(bandsp90, nbands)

                    elif self.key == 'perc_95':
                        
                        compositer = composite.select(bandsp95, nbands)

                    else:
                        
                        print('Uhm... it seems that you choose a bad statistic name')

                
                else: 
                    
                    if self.key not in ['perc_90', 'perc_95']:
                        
                        compositer = composite.select(sarbands, nbands)
                
                    elif self.key == 'perc_90':
        
                        compositer = composite.select(sarbandsp90, nbands)

                    elif self.key == 'perc_95':
                        
                        compositer = composite.select(sarbandsp95, nbands)

                    else:
                        
                        print('Uhm... it seems that you choose a bad statistic name')
                
                
                
                self.imagelist.append(compositer)
            
            
            ndvi_comp_coll = ee.ImageCollection.fromImages(self.imagelist)
        
            return ndvi_comp_coll

    #Get year composite for each periods or at the end just once?!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    # Let's add the 12 and 24 split periods

    #elif self.periods == 12:


    def get_export_single(self, image, name='mycomposition.tif', crs='EPSG:4326', scale=10):
            
        '''Export single composition as .tif to your local folder. So long as we can do really nice
        multiseasonal composites, e.g. median seasonal NDVI for whole Africa between 2001 and 2020. 
        I thought that would be interesting to have the chance to export this composites.
        This method calls geemap.ee_export_image.
        
        Args:
            crs (string): Coordinate Reference System of the output tifs. Use EPSG style e.g. "EPSG:4326"; "EPSG:32629"
            scale (int): Value for pixel size of your output tifs, default is 10 because default sat is Sentinel 2 NDVI,
                with bands 8 and 4 (10 m pix/resolution). In case you choose Landsat yous should change scale to 30. 
                You can also considerer use this parameter like resample in case you get a size limitation error
                when try to download a big area. 
        Returns:
          object(tif): 4 bands mycompositon.tif, with max, median or perc_90 for the period you chosen, 
          downloaded at your current working directory 
          '''
        
        filename = os.path.join(os.getcwd(), name)
        geemap.ee_export_image(image, filename=filename, scale=scale, crs=crs, region=self.roi, file_per_band=False) 

        print('Image have been exported')
        
        
    def get_export(self, crs='EPSG:4326', scale=10):
        
        '''Export NDVI year compositions as .tif to your local folder.
        This method calls geemap.ee_export_image.
        
        Args:
            crs (string): Coordinate Reference System of the output tifs. Use EPSG style e.g. "EPSG:4326"; "EPSG:32629"
            scale (int): Value for pixel size of your output tifs, default is 10 because default sat is Sentinel 2 NDVI,
                with bands 8 and 4 (10 m pix/resolution). In case you choose Landsat yous should change scale to 30. 
                You can also considerer use this parameter like resample in case you get a size limitation error
                when try to download a big area. 
        Returns:
          object(tif): 4 bands ndvi_year.tif per year in your ndvi collection, downloaded at your current working directory 
          '''

        self.imagelist = []
        self.get_year_composite()
            
        count = (len(self.imagelist))
        print(count)
            
        for n in range(count):
            year = self.start_year + n
            image = self.imagelist[n]
            name = 'ndvi_' + str(year) + '.tif'
            filename = os.path.join(os.getcwd(), name)
            print('Exporting {}'.format(filename), '\n')
            geemap.ee_export_image(image, filename=filename, scale=scale, crs=crs, region=self.roi, file_per_band=False) 

        print('All the images in the collection have been exported')


    def get_gif(self, name='mygif.gif', bands=['winter', 'spring', 'summer']):

        '''Export NDVI year compositions as .gif to your local folder. 
        This method calls geemap.download_ee_video & geemap.add_text_to_gif.

        Args:
            name (string): Name of the output gif. It will be saved at your current working directory.
            bands (list): List where you define the band combination for your gif.
        Returns:
            object(gif): myname.gif and myname_texted.gif downloaded at your current working directory
        '''
        
        self.imagelist = [0]
        self.get_year_composite()
        
        # We can't define an minimun and maximun for sar data, since it depends a lot on the region
        # Compute those values can take a lot of time, so we decided to keep a viz param standar for optical data
        # and just compute them for sar data, so we split viz params in two depending on sar or opical data
        # Of course, sar data it takes more time to be generated. This apply to get the gif not to export the images
        out_gif = os.path.join(os.getcwd(), name)
        self.imagelist = []
        
        if self.sat == 'sar':
        
            d = {'winter':self.dwinter, 'spring': self.dspring, 'summer': self.dsummer, 'autumn': self.dautumn}

            minimos = [self.get_perc(d[i][self.key])[0] for i in bands]
            maximos = [self.get_perc(d[i][self.key])[1] for i in bands]

            #min_ = self.get_perc(self.dwinter[self.key])[0]
            #max_ = self.get_perc(self.dwinter[self.key])[1]

            print(minimos, maximos)

            video_args = {
                'dimensions': 768,
                'region': self.roi, 
                'framesPerSecond': 10,
                'bands': bands, 
                'min': minimos,
                'max': maximos,
                'gamma': [1, 1, 1]}
            
        else:
            
            video_args = {
                'dimensions': 768,
                'region': self.roi, 
                'framesPerSecond': 10,
                'bands': bands, 
                'min': 0.15,
                'max': 0.85,
                'gamma': [1, 1, 1]}
        
        geemap.download_ee_video(self.get_year_composite(), video_args, out_gif)
        texted_gif = out_gif[:-4] + '_texted.gif'
        geemap.add_text_to_gif(out_gif, texted_gif, xy=('5%', '90%'), 
                                text_sequence=self.start_year, font_size=30, font_color='#ffffff', 
                                add_progress_bar=False, duration=300)


def get_stats(image, geom, name, stat='MEDIAN', scale=10):
    
        '''Compute zonal statistics for a local shapefile or Map user Geometries. It could be used like a Point Sampling Tool also.
        This method calls geemap.zonal_statistics.
        
        Args:
            image (ee.Image or ee.ImageCollection): Coordinate Reference System of the output tifs. Use EPSG style e.g. "EPSG:4326"; "EPSG:32629"
            geom (Shapefile or ee.FeatureCollection): feature collection with the polygons or points to get the statistics
            name (string): Name of the output shapefile in your current folder
            stat (string): Name of the desired statistic. Available names are 'MEAN', 'MEDIAN', 'MIN' and 'MAX'
            scale(int): Resample pixel size in case it needed
        Returns:
          object(Shapefile): Shapefile with the statistics for each rasters bands in the Ndvi2Gif composite, 
          downloaded at your current working directory 
          '''
        
        out_shp = os.path.join(os.getcwd(), name + '.shp')
        if geom is None:
            # When no geometry is passed we use a default area over Donana Natural Space
            print('Please select an area to compute the statistics')
            
        elif isinstance(geom, str):

            if geom.endswith('.shp'):
                roi = geemap.shp_to_ee(geom)
            elif roi.endswith('.geojson'):
                roi = geemap.geojson_to_ee(geom)
            else:
                print('It seems that your path is broken. Remember that await for shapefiles or geojson')

            
            geemap.zonal_statistics(image, roi, out_shp, stat, scale)


        else:

            if not isinstance(geom, ee.Geometry):

                try:
                    roi = geom.geometry()
                except Exception as e:
                    print('Could not convert the provided roi to ee.Geometry')
                    print(e)
                    return
                
            else:

                try:
                    roi = geom
                except Exception as e:
                    print('Could not convert the provided roi to ee.Geometry')
                    print(e)
                    return


            out_shp = os.path.join(os.getcwd(), 'stats.shp')
            geemap.zonal_statistics(image, roi=roi, out_shp=out_shp, statistics_type=stat, scale=scale)

