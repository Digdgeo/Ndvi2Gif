import os
import ee
import geemap
import requests
import zipfile
import geopandas as gpd
import fiona
from geemap import zonal_statistics
from io import BytesIO
        

def scale_OLI(image):

    """
    Scale Landsat 8 OLI surface reflectance bands to reflectance values.
    
    Applies scaling factors and offsets to convert Landsat 8 OLI digital numbers
    to surface reflectance values and renames bands to standardized names.
    
    Parameters
    ----------
    image : ee.Image
        Landsat 8 OLI image with surface reflectance bands (SR_B2 through SR_B7).
        
    Returns
    -------
    ee.Image
        Image with scaled optical bands added, renamed to standard band names
        (Blue, Green, Red, Nir, Swir1, Swir2).
        
    Notes
    -----
    The scaling formula applied is: reflectance = DN * 0.0000275 - 0.2
    where DN is the digital number from the original band.
    
    Examples
    --------
    >>> # Scale a Landsat 8 OLI image
    >>> oli_image = ee.Image('LANDSAT/LC08/C02/T1_L2/LC08_044034_20140318')
    >>> scaled_image = scale_OLI(oli_image)
    """

    opticalBands = image.select(['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7']).multiply(0.0000275).add(-0.2).rename(['Blue', 'Green', 'Red', 'Nir', 'Swir1', 'Swir2'])
    #thermalBands = image.select('ST_B.*').multiply(0.00341802).add(149.0)
    return image.addBands(opticalBands, None, True)#.addBands(thermalBands, None, True)
    
def scale_ETM(image):

    """
    Scale Landsat 7 ETM+ surface reflectance bands to reflectance values.
    
    Applies scaling factors and offsets to convert Landsat 7 ETM+ digital numbers
    to surface reflectance values and renames bands to standardized names.
    
    Parameters
    ----------
    image : ee.Image
        Landsat 7 ETM+ image with surface reflectance bands (SR_B1 through SR_B7,
        excluding thermal band SR_B6).
        
    Returns
    -------
    ee.Image
        Image with scaled optical bands added, renamed to standard band names
        (Blue, Green, Red, Nir, Swir1, Swir2).
        
    Notes
    -----
    The scaling formula applied is: reflectance = DN * 0.0000275 - 0.2
    where DN is the digital number from the original band.
    
    Examples
    --------
    >>> # Scale a Landsat 7 ETM+ image
    >>> etm_image = ee.Image('LANDSAT/LE07/C02/T1_L2/LE07_044034_20030316')
    >>> scaled_image = scale_ETM(etm_image)
    """

    opticalBands = image.select(['SR_B1','SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7']).multiply(0.0000275).add(-0.2).rename(['Blue', 'Green', 'Red', 'Nir', 'Swir1', 'Swir2'])
    #thermalBands = image.select('ST_B.*').multiply(0.00341802).add(149.0)
    return image.addBands(opticalBands, None, True)#.addBands(thermalBands, None, True)

class NdviSeasonality:
    
    """
    Generate remote sensing index seasonal composition GIFs and images.
    
    This class creates seasonal composites of vegetation indices using Google Earth Engine.
    Seasons are defined as fixed parameters and statistical reducers (max, median, etc.)
    are used to generate the seasonal value for each pixel. The resulting composites
    show phenological patterns across seasons and years.
    
    Parameters
    ----------
    roi : ee.Geometry, str, or None, optional
        Region of interest for analysis. Can be:
        - ee.Geometry object
        - Path to shapefile (.shp)
        - Path to GeoJSON file (.geojson)
        - DEIMS site ID (format: 'deimsid:site_id')
        - Landsat WRS-2 tile (format: 'wrs:path,row')
        - Sentinel-2 MGRS tile (format: 's2:tile_id')
        If None, uses default area over Doñana Natural Space.
    periods : int, optional
        Number of seasonal periods to divide the year into. Default is 4.
    start_year : int, optional
        First year of the analysis period. Default is 2016.
    end_year : int, optional
        Last year of the analysis period. Default is 2020.
    sat : {'S2', 'Landsat', 'MODIS', 'S1'}, optional
        Satellite collection to use:
        - 'S2': Sentinel-2 (2015-present)
        - 'Landsat': Landsat 4 TM, 5 TM, 7 ETM+, 8 OLI (1984-present)
        - 'MODIS': MODIS Terra/Aqua
        - 'S1': Sentinel-1 SAR
        Default is 'S2'.
    key : {'max', 'median', 'perc_90', 'perc_95', 'mean'}, optional
        Statistical reducer to apply for seasonal composites. Default is 'max'.
    index : str, optional
        Vegetation index to calculate. Default is 'ndvi'. Available indices
        include: ndvi, ndwi, mndwi, evi, savi, gndvi, avi, nbri, ndsi, aweinsh,
        awei, ndmi.
        
    Attributes
    ----------
    roi : ee.Geometry
        Processed region of interest geometry.
    periods : int
        Number of seasonal divisions.
    start_year : int
        Start year of analysis.
    end_year : int
        End year of analysis.
    sat : str
        Selected satellite collection.
    key : str
        Selected statistical reducer.
    imagelist : list
        List to store processed images.
    index : str
        Selected vegetation index.
    d : dict
        Dictionary mapping index names to calculation methods.
        
    Examples
    --------
    >>> # Create NDVI seasonality analysis for default area
    >>> ndvi_season = NdviSeasonality()
    
    >>> # Use custom ROI with shapefile
    >>> ndvi_season = NdviSeasonality(roi='my_area.shp', start_year=2018, end_year=2021)

    >>> # Use DEIMS ID with shapefile
    >>> ndvi_season = NdviSeasonality(roi='deimsid:https://deims.org/bcbc866c-3f4f-47a8-bbbc-0a93df6de7b2', start_year=2018, end_year=2021)
    
    >>> # Use Landsat WRS-2 tile
    >>> ndvi_season = NdviSeasonality(roi='wrs:202,34', sat='Landsat', key='median')
    
    >>> # Use Sentinel-2 MGRS tile
    >>> ndvi_season = NdviSeasonality(roi='s2:29SQA', index='evi', periods=6)
    
    Notes
    -----
    The class combines seasonal composites into multi-band rasters where each band
    represents one season. Color combinations in the resulting images reveal
    phenological patterns over seasons and years.
    
    Seasonal periods are calculated as equal divisions of the year. For example,
    with periods=4: Spring (Mar-May), Summer (Jun-Aug), Autumn (Sep-Nov), 
    Winter (Dec-Feb).
    """
                                                            
    def __init__(self, roi=None, periods=4, start_year=2016, end_year=2020, sat='S2', key='max', index='ndvi'):

        """
        Initialize NdviSeasonality instance.
        
        Parameters
        ----------
        roi : ee.Geometry, str, or None, optional
            Region of interest for analysis. Can be:
            - ee.Geometry object
            - Path to shapefile (.shp)
            - Path to GeoJSON file (.geojson)
            - DEIMS site ID (format: 'deimsid:site_id')
            - Landsat WRS-2 tile (format: 'wrs:path,row')
            - Sentinel-2 MGRS tile (format: 's2:tile_id')
            If None, uses default area over Doñana Natural Space.
        periods : int, optional
            Number of seasonal periods to divide the year into. Default is 4.
        start_year : int, optional
            First year of the analysis period. Default is 2016.
        end_year : int, optional
            Last year of the analysis period. Default is 2020.
        sat : {'S2', 'Landsat', 'MODIS', 'S1'}, optional
            Satellite collection to use. Default is 'S2'.
        key : {'max', 'median', 'perc_90', 'perc_95', 'mean'}, optional
            Statistical reducer to apply for seasonal composites. Default is 'max'.
        index : str, optional
            Vegetation index to calculate. Default is 'ndvi'.
        """

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
                print('Con DEIMS hemos topado amigo Sancho...')
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
                print('It seems that your ROI path is invalid.\n\nAccepted formats:\n'
                '- Shapefile path (e.g. "myfile.shp")\n'
                '- GeoJSON path (e.g. "area.geojson")\n'
                '- DEIMS site ID (e.g. "deimsid:bcbc866c-3f4f-47a8-bbbc-0a93df6de7b2")\n'
                '- Landsat WRS2 tile (e.g. "wrs:202,34")\n'
                '- Drawn geometry on map (if using geemap interactive tools)')

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
        self.d = {'ndvi': self.get_ndvi, 'ndwi': self.get_ndwi, 'mndwi': self.get_mndwi, 'evi': self.get_evi, 
                'savi': self.get_savi, 'gndvi': self.get_gndvi, 'avi': self.get_avi, 'nbri': self.get_nbri, 
                'ndsi': self.get_ndsi, 'aweinsh': self.get_aweinsh, 'awei': self.get_awei, 'ndmi': self.get_ndmi}

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

        """
        Calculate Normalized Difference Vegetation Index (NDVI).
        
        NDVI is calculated as (NIR - Red) / (NIR + Red) and is used to assess
        vegetation health and density. Values range from -1 to +1, where higher
        values indicate healthier vegetation.
        
        Parameters
        ----------
        image : ee.Image
            Input image with 'Nir' and 'Red' bands.
            
        Returns
        -------
        ee.Image
            Single-band image with NDVI values.
            
        Notes
        -----
        NDVI formula: (NIR - Red) / (NIR + Red)
        
        Typical value ranges:
        - Water: -1 to 0
        - Bare soil: 0 to 0.2
        - Sparse vegetation: 0.2 to 0.4
        - Dense vegetation: 0.4 to 1.0
        
        Examples
        --------
        >>> ndvi_image = self.get_ndvi(scaled_image)
        """   

        return image.normalizedDifference(['Nir', 'Red'])
    
    def get_ndwi(self, image):

        """
        Calculate Normalized Difference Water Index (NDWI).
        
        NDWI is used to delineate open water features and enhance their presence
        in remotely sensed imagery. It's calculated using Green and NIR bands.
        
        Parameters
        ----------
        image : ee.Image
            Input image with 'Green' and 'Nir' bands.
            
        Returns
        -------
        ee.Image
            Single-band image with NDWI values.
            
        Notes
        -----
        NDWI formula: (Green - NIR) / (Green + NIR)
        
        Positive values typically indicate water bodies, while negative values
        indicate vegetation and dry surfaces.
        
        Examples
        --------
        >>> ndwi_image = self.get_ndwi(scaled_image)
        """   

        return image.normalizedDifference(['Green', 'Nir'])

    def get_mndwi(self, image):

        """
        Calculate Modified Normalized Difference Water Index (MNDWI).
        
        MNDWI is a modification of NDWI that uses SWIR1 instead of NIR,
        providing better discrimination between water and built-up areas.
        It's particularly effective for urban water body detection.
        
        Parameters
        ----------
        image : ee.Image
            Input image with 'Green' and 'Swir1' bands.
            
        Returns
        -------
        ee.Image
            Single-band image with MNDWI values.
            
        Notes
        -----
        MNDWI formula: (Green - SWIR1) / (Green + SWIR1)
        
        MNDWI typically performs better than NDWI for:
        - Urban water body detection
        - Distinguishing water from built-up areas
        - Areas with mixed water and urban features
        
        Positive values indicate water bodies, while negative values indicate
        non-water surfaces. MNDWI is less affected by built-up noise compared
        to traditional NDWI.
        
        Examples
        --------
        >>> mndwi_image = self.get_mndwi(scaled_image)
        """

        return image.normalizedDifference(['Green', 'Swir1'])
    
    def get_evi(self, image):
    
        """
        Calculate Enhanced Vegetation Index (EVI).
        
        EVI is an optimized vegetation index that minimizes canopy background
        variations and maintains sensitivity over dense vegetation conditions.
        
        Parameters
        ----------
        image : ee.Image
            Input image with 'Nir', 'Red', and 'Blue' bands.
            
        Returns
        -------
        ee.Image
            Single-band image with EVI values.
            
        Notes
        -----
        EVI formula: 2.5 * ((NIR - Red) / (NIR + 6 * Red - 7.5 * Blue + 1))
        
        EVI provides better sensitivity in high biomass regions and improved
        vegetation monitoring through a de-coupling of the canopy background
        signal and a reduction in atmosphere influences.
        
        Examples
        --------
        >>> evi_image = self.get_evi(scaled_image)
        """   

        # Compute the EVI using an expression.
        return image.expression(
            '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))', {
            'NIR': image.select('Nir'),
            'RED': image.select('Red'),
            'BLUE': image.select('Blue')}).rename(['nd'])
    
    def get_savi(self, image, L=0.428):
        
        """
        Calculate Soil Adjusted Vegetation Index (SAVI).
        
        SAVI is a vegetation index that attempts to minimize soil brightness
        influences using a soil brightness correction factor (L).
        
        Parameters
        ----------
        image : ee.Image
            Input image with 'Nir' and 'Red' bands.
            
        Returns
        -------
        ee.Image
            Single-band image with SAVI values.
            
        Notes
        -----
        SAVI formula: ((NIR - Red) / (NIR + Red + L)) * (1 + L)
        where L = 0.5 (soil brightness correction factor)
        
        The L factor varies based on vegetation density:
        - L = 1 for low vegetation cover
        - L = 0.5 for intermediate vegetation cover
        - L = 0.25 for high vegetation cover
        
        Examples
        --------
        >>> savi_image = self.get_savi(scaled_image)
        """   

        # Compute the SAVI using an expression.
        return image.expression(
            '((NIR - RED) / (NIR + RED + L) * (1 + L))', {
            'NIR': image.select('Nir'),
            'RED': image.select('Red'),
            'L': L}).rename(['nd'])
    
    def get_aweinsh(self, image):
        
        """
        Calculate Automated Water Extraction Index - No Shadow (AWEInsh).
        
        AWEInsh is designed to extract water bodies while minimizing shadow
        effects that can cause confusion in water detection algorithms.
        
        Parameters
        ----------
        image : ee.Image
            Input image with 'Blue', 'Green', 'Nir', 'Swir1', and 'Swir2' bands.
            
        Returns
        -------
        ee.Image
            Single-band image with AWEInsh values.
            
        Notes
        -----
        AWEInsh formula: 4 * (Green - SWIR1) - (0.25 * NIR + 2.75 * SWIR2)
        
        Positive values typically indicate water bodies, while negative values
        indicate non-water surfaces. The index is particularly effective in
        urban and mountainous areas where shadows are prevalent.
        
        Examples
        --------
        >>> aweinsh_image = self.get_aweinsh(scaled_image)
        """   

        # Compute the SAVI using an expression.
        return image.expression(
            '4.0 * (GREEN - SWIR1) - 0.25 * NIR + 2.75 * SWIR2', {
            'NIR': image.select('Nir'),
            'GREEN': image.select('Green'),
            'SWIR1':image.select('Swir1'),
            'SWIR2':image.select('Swir2')}).rename(['nd'])
    
    def get_awei(self, image):
        
        """
        Calculate Automated Water Extraction Index (AWEI).
        
        AWEI is designed for automatic water extraction from satellite imagery,
        particularly effective in urban environments and areas with built-up surfaces.
        
        Parameters
        ----------
        image : ee.Image
            Input image with 'Blue', 'Green', 'Nir', 'Swir1', and 'Swir2' bands.
            
        Returns
        -------
        ee.Image
            Single-band image with AWEI values.
            
        Notes
        -----
        AWEI formula: Blue + 2.5 * Green - 1.5 * (NIR + SWIR1) - 0.25 * SWIR2
        
        Positive values typically indicate water bodies. AWEI is particularly
        useful for water extraction in complex environments with mixed land cover.
        
        Examples
        --------
        >>> awei_image = self.get_awei(scaled_image)
        """   

        # Compute the SAVI using an expression.
        return image.expression(
            ('BLUE + 2.5 * GREEN - 1.5 * (NIR + SWIR1) - 0.25 * SWIR2'), {
            'NIR': image.select('Nir'),
            'BLUE': image.select('Blue'),
            'GREEN': image.select('Green'),
            'SWIR1':image.select('Swir1'),
            'SWIR2':image.select('Swir2')}).rename(['nd'])

    def get_gndvi(self, image):
        
        """
        Calculate Green Normalized Difference Vegetation Index (GNDVI).
        
        GNDVI uses the green band instead of the red band and is more sensitive
        to chlorophyll concentration than NDVI.
        
        Parameters
        ----------
        image : ee.Image
            Input image with 'Nir' and 'Green' bands.
            
        Returns
        -------
        ee.Image
            Single-band image with GNDVI values.
            
        Notes
        -----
        GNDVI formula: (NIR - Green) / (NIR + Green)
        
        GNDVI is more sensitive to chlorophyll-a than NDVI and can be useful
        for assessing photosynthetic activity and nitrogen content in vegetation.
        
        Examples
        --------
        >>> gndvi_image = self.get_gndvi(scaled_image)
        """   

        # Compute the GNDVI using an expression.
        return image.normalizedDifference(['Nir', 'Green'])
    
    
    def get_avi(self, image, L=0.428):
        
        """
        Calculate Advanced Vegetation Index (AVI).
        
        AVI is designed to be more sensitive to vegetation than traditional
        vegetation indices, particularly in areas with sparse vegetation cover.
        
        Parameters
        ----------
        image : ee.Image
            Input image with 'Nir' and 'Red' bands.
            
        Returns
        -------
        ee.Image
            Single-band image with AVI values.
            
        Notes
        -----
        AVI formula: (NIR * (1 - Red) * (NIR - Red))^(1/3)
        
        AVI is particularly useful for monitoring vegetation in arid and
        semi-arid environments where vegetation cover is sparse.
        
        Examples
        --------
        >>> avi_image = self.get_avi(scaled_image)
        """   

        # Compute the SAVI using an expression.
        return image.expression(
            '(NIR * (1.0 - RED) * (NIR - RED)) ** (1/3)', {
            'NIR': image.select('Nir'),
            'RED': image.select('Red')}).rename(['nd'])

    def get_nbri(self, image):
        
        """
        Calculate Normalized Burn Ratio Index (NBRI).
        
        NBRI is used to identify burned areas and assess burn severity by
        highlighting the difference between healthy vegetation and burned areas.
        
        Parameters
        ----------
        image : ee.Image
            Input image with 'Nir' and 'Swir2' bands.
            
        Returns
        -------
        ee.Image
            Single-band image with NBRI values.
            
        Notes
        -----
        NBRI formula: (NIR - SWIR2) / (NIR + SWIR2)
        
        High values indicate healthy vegetation, while low values indicate
        burned areas. The index is particularly useful for post-fire assessment
        and monitoring vegetation recovery.
        
        Examples
        --------
        >>> nbri_image = self.get_nbri(scaled_image)
        """   

        # Compute the EVI using an expression.
        return image.normalizedDifference(['Nir', 'Swir2'])
    
    def get_ndsi(self, image):
        
        """
        Calculate Normalized Difference Snow Index (NDSI).
        
        NDSI is used to identify snow cover by exploiting the difference in
        reflectance between visible and shortwave infrared bands.
        
        Parameters
        ----------
        image : ee.Image
            Input image with 'Green' and 'Swir1' bands.
            
        Returns
        -------
        ee.Image
            Single-band image with NDSI values.
            
        Notes
        -----
        NDSI formula: (Green - SWIR1) / (Green + SWIR1)
        
        Values above 0.4 typically indicate snow cover, while values below
        0.4 indicate snow-free surfaces. The index is widely used for snow
        mapping and monitoring seasonal snow cover changes.
        
        Examples
        --------
        >>> ndsi_image = self.get_ndsi(scaled_image)
        """   

        # Compute the EVI using an expression.
        return image.normalizedDifference(['Green', 'Swir1'])

    def get_ndmi(self, image):
        
        """
        Calculate Normalized Difference Moisture Index (NDMI).
        
        NDMI is sensitive to moisture levels in vegetation and can be used to
        monitor drought conditions and vegetation water stress.
        
        Parameters
        ----------
        image : ee.Image
            Input image with 'Nir' and 'Swir1' bands.
            
        Returns
        -------
        ee.Image
            Single-band image with NDMI values.
            
        Notes
        -----
        NDMI formula: (NIR - SWIR1) / (NIR + SWIR1)
        
        Values range from -1 to +1, where:
        - High values indicate high moisture content
        - Low values indicate low moisture content or stress
        - Negative values may indicate bare soil or very dry conditions
        
        Examples
        --------
        >>> ndmi_image = self.get_ndmi(scaled_image)
        """   

        # Compute the EVI using an expression.
        return image.normalizedDifference(['Nir', 'Swir1'])

    
    # Here we define the functions to get the seasonal periods of time

    def get_winter(self, y):
        
        """
        Generates the seasonal image for the winter period of the specified year.

        This function filters the input image collection using the date range defined
        in `self.winter` for the given year. It computes statistical summaries
        (maximum, median, mean, 90th percentile, and 95th percentile) over the images
        in that period. For optical sensors, the selected index function is applied
        before aggregation; for radar data (e.g., Sentinel-1), raw values are used.

        Parameters
        ----------
        y : int or str
            Year for which the winter period image should be computed.

        Returns
        -------
        ee.Image
            The image corresponding to the selected statistical metric,
            defined in `self.key`.
        """

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
        
        """
        Generates the seasonal image for the winter period of the specified year.

        This function filters the input image collection using the date range defined
        in `self.spring` for the given year. It computes statistical summaries
        (maximum, median, mean, 90th percentile, and 95th percentile) over the images
        in that period. For optical sensors, the selected index function is applied
        before aggregation; for radar data (e.g., Sentinel-1), raw values are used.

        Parameters
        ----------
        y : int or str
            Year for which the winter period image should be computed.

        Returns
        -------
        ee.Image
            The image corresponding to the selected statistical metric,
            defined in `self.key`.
        """

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
        
        """
        Generates the seasonal image for the winter period of the specified year.

        This function filters the input image collection using the date range defined
        in `self.summer` for the given year. It computes statistical summaries
        (maximum, median, mean, 90th percentile, and 95th percentile) over the images
        in that period. For optical sensors, the selected index function is applied
        before aggregation; for radar data (e.g., Sentinel-1), raw values are used.

        Parameters
        ----------
        y : int or str
            Year for which the winter period image should be computed.

        Returns
        -------
        ee.Image
            The image corresponding to the selected statistical metric,
            defined in `self.key`.
        """

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
                
        """
        Generates the seasonal image for the winter period of the specified year.

        This function filters the input image collection using the date range defined
        in `self.autumn` for the given year. It computes statistical summaries
        (maximum, median, mean, 90th percentile, and 95th percentile) over the images
        in that period. For optical sensors, the selected index function is applied
        before aggregation; for radar data (e.g., Sentinel-1), raw values are used.

        Parameters
        ----------
        y : int or str
            Year for which the winter period image should be computed.

        Returns
        -------
        ee.Image
            The image corresponding to the selected statistical metric,
            defined in `self.key`.
        """

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
                
        """
        Generates the seasonal image for the custom period p1 of the specified year.

        This function filters the input image collection using the date range defined
        in january for the given year. It computes statistical summaries
        (maximum, median, mean, 90th percentile, and 95th percentile) over the images
        in that period. For optical sensors, the selected index function is applied
        before aggregation; for radar data (e.g., Sentinel-1), raw values are used.

        Parameters
        ----------
        y : int or str
            Year for which the custom period p1 image should be computed.

        Returns
        -------
        ee.Image
            The image corresponding to the selected statistical metric,
            defined in `self.key`.
        """

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
        
        """
        Generates the seasonal image for the custom period p1 of the specified year.

        This function filters the input image collection using the date range defined
        in february for the given year. It computes statistical summaries
        (maximum, median, mean, 90th percentile, and 95th percentile) over the images
        in that period. For optical sensors, the selected index function is applied
        before aggregation; for radar data (e.g., Sentinel-1), raw values are used.

        Parameters
        ----------
        y : int or str
            Year for which the custom period p1 image should be computed.

        Returns
        -------
        ee.Image
            The image corresponding to the selected statistical metric,
            defined in `self.key`.
        """

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
        
        """
        Generates the seasonal image for the custom period p1 of the specified year.

        This function filters the input image collection using the date range defined
        in march for the given year. It computes statistical summaries
        (maximum, median, mean, 90th percentile, and 95th percentile) over the images
        in that period. For optical sensors, the selected index function is applied
        before aggregation; for radar data (e.g., Sentinel-1), raw values are used.

        Parameters
        ----------
        y : int or str
            Year for which the custom period p1 image should be computed.

        Returns
        -------
        ee.Image
            The image corresponding to the selected statistical metric,
            defined in `self.key`.
        """

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
        
        """
        Generates the seasonal image for the custom period p1 of the specified year.

        This function filters the input image collection using the date range defined
        in april for the given year. It computes statistical summaries
        (maximum, median, mean, 90th percentile, and 95th percentile) over the images
        in that period. For optical sensors, the selected index function is applied
        before aggregation; for radar data (e.g., Sentinel-1), raw values are used.

        Parameters
        ----------
        y : int or str
            Year for which the custom period p1 image should be computed.

        Returns
        -------
        ee.Image
            The image corresponding to the selected statistical metric,
            defined in `self.key`.
        """

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
        
        """
        Generates the seasonal image for the custom period p1 of the specified year.

        This function filters the input image collection using the date range defined
        in may for the given year. It computes statistical summaries
        (maximum, median, mean, 90th percentile, and 95th percentile) over the images
        in that period. For optical sensors, the selected index function is applied
        before aggregation; for radar data (e.g., Sentinel-1), raw values are used.

        Parameters
        ----------
        y : int or str
            Year for which the custom period p1 image should be computed.

        Returns
        -------
        ee.Image
            The image corresponding to the selected statistical metric,
            defined in `self.key`.
        """

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
        
        """
        Generates the seasonal image for the custom period p1 of the specified year.

        This function filters the input image collection using the date range defined
        in june for the given year. It computes statistical summaries
        (maximum, median, mean, 90th percentile, and 95th percentile) over the images
        in that period. For optical sensors, the selected index function is applied
        before aggregation; for radar data (e.g., Sentinel-1), raw values are used.

        Parameters
        ----------
        y : int or str
            Year for which the custom period p1 image should be computed.

        Returns
        -------
        ee.Image
            The image corresponding to the selected statistical metric,
            defined in `self.key`.
        """

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
        
        """
        Generates the seasonal image for the custom period p1 of the specified year.

        This function filters the input image collection using the date range defined
        in july for the given year. It computes statistical summaries
        (maximum, median, mean, 90th percentile, and 95th percentile) over the images
        in that period. For optical sensors, the selected index function is applied
        before aggregation; for radar data (e.g., Sentinel-1), raw values are used.

        Parameters
        ----------
        y : int or str
            Year for which the custom period p1 image should be computed.

        Returns
        -------
        ee.Image
            The image corresponding to the selected statistical metric,
            defined in `self.key`.
        """
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
        
        """
        Generates the seasonal image for the custom period p1 of the specified year.

        This function filters the input image collection using the date range defined
        in august for the given year. It computes statistical summaries
        (maximum, median, mean, 90th percentile, and 95th percentile) over the images
        in that period. For optical sensors, the selected index function is applied
        before aggregation; for radar data (e.g., Sentinel-1), raw values are used.

        Parameters
        ----------
        y : int or str
            Year for which the custom period p1 image should be computed.

        Returns
        -------
        ee.Image
            The image corresponding to the selected statistical metric,
            defined in `self.key`.
        """

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
            
            """
        Generates the seasonal image for the custom period p1 of the specified year.

        This function filters the input image collection using the date range defined
        in september for the given year. It computes statistical summaries
        (maximum, median, mean, 90th percentile, and 95th percentile) over the images
        in that period. For optical sensors, the selected index function is applied
        before aggregation; for radar data (e.g., Sentinel-1), raw values are used.

        Parameters
        ----------
        y : int or str
            Year for which the custom period p1 image should be computed.

        Returns
        -------
        ee.Image
            The image corresponding to the selected statistical metric,
            defined in `self.key`.
        """

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
        
        """
        Generates the seasonal image for the custom period p1 of the specified year.

        This function filters the input image collection using the date range defined
        in october for the given year. It computes statistical summaries
        (maximum, median, mean, 90th percentile, and 95th percentile) over the images
        in that period. For optical sensors, the selected index function is applied
        before aggregation; for radar data (e.g., Sentinel-1), raw values are used.

        Parameters
        ----------
        y : int or str
            Year for which the custom period p1 image should be computed.

        Returns
        -------
        ee.Image
            The image corresponding to the selected statistical metric,
            defined in `self.key`.
        """

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
        
        """
        Generates the seasonal image for the custom period p1 of the specified year.

        This function filters the input image collection using the date range defined
        in november for the given year. It computes statistical summaries
        (maximum, median, mean, 90th percentile, and 95th percentile) over the images
        in that period. For optical sensors, the selected index function is applied
        before aggregation; for radar data (e.g., Sentinel-1), raw values are used.

        Parameters
        ----------
        y : int or str
            Year for which the custom period p1 image should be computed.

        Returns
        -------
        ee.Image
            The image corresponding to the selected statistical metric,
            defined in `self.key`.
        """

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
        
        """
        Generates the seasonal image for the custom period p1 of the specified year.

        This function filters the input image collection using the date range defined
        in december for the given year. It computes statistical summaries
        (maximum, median, mean, 90th percentile, and 95th percentile) over the images
        in that period. For optical sensors, the selected index function is applied
        before aggregation; for radar data (e.g., Sentinel-1), raw values are used.

        Parameters
        ----------
        y : int or str
            Year for which the custom period p1 image should be computed.

        Returns
        -------
        ee.Image
            The image corresponding to the selected statistical metric,
            defined in `self.key`.
        """

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
                
        """
        Generates the seasonal image for the custom period p1 of the specified year.

        This function filters the input image collection using the date range defined
        in `self.p1` for the given year. It computes statistical summaries
        (maximum, median, mean, 90th percentile, and 95th percentile) over the images
        in that period. For optical sensors, the selected index function is applied
        before aggregation; for radar data (e.g., Sentinel-1), raw values are used.

        Parameters
        ----------
        y : int or str
            Year for which the custom period p1 image should be computed.

        Returns
        -------
        ee.Image
            The image corresponding to the selected statistical metric,
            defined in `self.key`.
        """

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
        
        """
        Generates the seasonal image for the custom period p1 of the specified year.

        This function filters the input image collection using the date range defined
        in `self.p2` for the given year. It computes statistical summaries
        (maximum, median, mean, 90th percentile, and 95th percentile) over the images
        in that period. For optical sensors, the selected index function is applied
        before aggregation; for radar data (e.g., Sentinel-1), raw values are used.

        Parameters
        ----------
        y : int or str
            Year for which the custom period p1 image should be computed.

        Returns
        -------
        ee.Image
            The image corresponding to the selected statistical metric,
            defined in `self.key`.
        """

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
        
        """
        Generates the seasonal image for the custom period p1 of the specified year.

        This function filters the input image collection using the date range defined
        in `self.p3` for the given year. It computes statistical summaries
        (maximum, median, mean, 90th percentile, and 95th percentile) over the images
        in that period. For optical sensors, the selected index function is applied
        before aggregation; for radar data (e.g., Sentinel-1), raw values are used.

        Parameters
        ----------
        y : int or str
            Year for which the custom period p1 image should be computed.

        Returns
        -------
        ee.Image
            The image corresponding to the selected statistical metric,
            defined in `self.key`.
        """

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
        
        """
        Generates the seasonal image for the custom period p1 of the specified year.

        This function filters the input image collection using the date range defined
        in `self.p4` for the given year. It computes statistical summaries
        (maximum, median, mean, 90th percentile, and 95th percentile) over the images
        in that period. For optical sensors, the selected index function is applied
        before aggregation; for radar data (e.g., Sentinel-1), raw values are used.

        Parameters
        ----------
        y : int or str
            Year for which the custom period p1 image should be computed.

        Returns
        -------
        ee.Image
            The image corresponding to the selected statistical metric,
            defined in `self.key`.
        """

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
        
        """
        Generates the seasonal image for the custom period p1 of the specified year.

        This function filters the input image collection using the date range defined
        in `self.p5` for the given year. It computes statistical summaries
        (maximum, median, mean, 90th percentile, and 95th percentile) over the images
        in that period. For optical sensors, the selected index function is applied
        before aggregation; for radar data (e.g., Sentinel-1), raw values are used.

        Parameters
        ----------
        y : int or str
            Year for which the custom period p1 image should be computed.

        Returns
        -------
        ee.Image
            The image corresponding to the selected statistical metric,
            defined in `self.key`.
        """

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
        
        """
        Generates the seasonal image for the custom period p1 of the specified year.

        This function filters the input image collection using the date range defined
        in `self.p6` for the given year. It computes statistical summaries
        (maximum, median, mean, 90th percentile, and 95th percentile) over the images
        in that period. For optical sensors, the selected index function is applied
        before aggregation; for radar data (e.g., Sentinel-1), raw values are used.

        Parameters
        ----------
        y : int or str
            Year for which the custom period p1 image should be computed.

        Returns
        -------
        ee.Image
            The image corresponding to the selected statistical metric,
            defined in `self.key`.
        """

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
        
        """
        Generates the seasonal image for the custom period p1 of the specified year.

        This function filters the input image collection using the date range defined
        in `self.p7` for the given year. It computes statistical summaries
        (maximum, median, mean, 90th percentile, and 95th percentile) over the images
        in that period. For optical sensors, the selected index function is applied
        before aggregation; for radar data (e.g., Sentinel-1), raw values are used.

        Parameters
        ----------
        y : int or str
            Year for which the custom period p1 image should be computed.

        Returns
        -------
        ee.Image
            The image corresponding to the selected statistical metric,
            defined in `self.key`.
        """

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
        
        """
        Generates the seasonal image for the custom period p1 of the specified year.

        This function filters the input image collection using the date range defined
        in `self.p8` for the given year. It computes statistical summaries
        (maximum, median, mean, 90th percentile, and 95th percentile) over the images
        in that period. For optical sensors, the selected index function is applied
        before aggregation; for radar data (e.g., Sentinel-1), raw values are used.

        Parameters
        ----------
        y : int or str
            Year for which the custom period p1 image should be computed.

        Returns
        -------
        ee.Image
            The image corresponding to the selected statistical metric,
            defined in `self.key`.
        """

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
            
        """
        Generates the seasonal image for the custom period p1 of the specified year.

        This function filters the input image collection using the date range defined
        in `self.p9` for the given year. It computes statistical summaries
        (maximum, median, mean, 90th percentile, and 95th percentile) over the images
        in that period. For optical sensors, the selected index function is applied
        before aggregation; for radar data (e.g., Sentinel-1), raw values are used.

        Parameters
        ----------
        y : int or str
            Year for which the custom period p1 image should be computed.

        Returns
        -------
        ee.Image
            The image corresponding to the selected statistical metric,
            defined in `self.key`.
        """

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
        
        """
        Generates the seasonal image for the custom period p1 of the specified year.

        This function filters the input image collection using the date range defined
        in `self.p10` for the given year. It computes statistical summaries
        (maximum, median, mean, 90th percentile, and 95th percentile) over the images
        in that period. For optical sensors, the selected index function is applied
        before aggregation; for radar data (e.g., Sentinel-1), raw values are used.

        Parameters
        ----------
        y : int or str
            Year for which the custom period p1 image should be computed.

        Returns
        -------
        ee.Image
            The image corresponding to the selected statistical metric,
            defined in `self.key`.
        """

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
        
        """
        Generates the seasonal image for the custom period p1 of the specified year.

        This function filters the input image collection using the date range defined
        in `self.p11` for the given year. It computes statistical summaries
        (maximum, median, mean, 90th percentile, and 95th percentile) over the images
        in that period. For optical sensors, the selected index function is applied
        before aggregation; for radar data (e.g., Sentinel-1), raw values are used.

        Parameters
        ----------
        y : int or str
            Year for which the custom period p1 image should be computed.

        Returns
        -------
        ee.Image
            The image corresponding to the selected statistical metric,
            defined in `self.key`.
        """

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
        
        """
        Generates the seasonal image for the custom period p1 of the specified year.

        This function filters the input image collection using the date range defined
        in `self.p12` for the given year. It computes statistical summaries
        (maximum, median, mean, 90th percentile, and 95th percentile) over the images
        in that period. For optical sensors, the selected index function is applied
        before aggregation; for radar data (e.g., Sentinel-1), raw values are used.

        Parameters
        ----------
        y : int or str
            Year for which the custom period p1 image should be computed.

        Returns
        -------
        ee.Image
            The image corresponding to the selected statistical metric,
            defined in `self.key`.
        """

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
            
        """
        Generates the seasonal image for the custom period p1 of the specified year.

        This function filters the input image collection using the date range defined
        in `self.p13` for the given year. It computes statistical summaries
        (maximum, median, mean, 90th percentile, and 95th percentile) over the images
        in that period. For optical sensors, the selected index function is applied
        before aggregation; for radar data (e.g., Sentinel-1), raw values are used.

        Parameters
        ----------
        y : int or str
            Year for which the custom period p1 image should be computed.

        Returns
        -------
        ee.Image
            The image corresponding to the selected statistical metric,
            defined in `self.key`.
        """

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
        
        """
        Generates the seasonal image for the custom period p1 of the specified year.

        This function filters the input image collection using the date range defined
        in `self.p14` for the given year. It computes statistical summaries
        (maximum, median, mean, 90th percentile, and 95th percentile) over the images
        in that period. For optical sensors, the selected index function is applied
        before aggregation; for radar data (e.g., Sentinel-1), raw values are used.

        Parameters
        ----------
        y : int or str
            Year for which the custom period p1 image should be computed.

        Returns
        -------
        ee.Image
            The image corresponding to the selected statistical metric,
            defined in `self.key`.
        """

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
        
        """
        Generates the seasonal image for the custom period p1 of the specified year.

        This function filters the input image collection using the date range defined
        in `self.p15` for the given year. It computes statistical summaries
        (maximum, median, mean, 90th percentile, and 95th percentile) over the images
        in that period. For optical sensors, the selected index function is applied
        before aggregation; for radar data (e.g., Sentinel-1), raw values are used.

        Parameters
        ----------
        y : int or str
            Year for which the custom period p1 image should be computed.

        Returns
        -------
        ee.Image
            The image corresponding to the selected statistical metric,
            defined in `self.key`.
        """

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
        
        """
        Generates the seasonal image for the custom period p1 of the specified year.

        This function filters the input image collection using the date range defined
        in `self.p16` for the given year. It computes statistical summaries
        (maximum, median, mean, 90th percentile, and 95th percentile) over the images
        in that period. For optical sensors, the selected index function is applied
        before aggregation; for radar data (e.g., Sentinel-1), raw values are used.

        Parameters
        ----------
        y : int or str
            Year for which the custom period p1 image should be computed.

        Returns
        -------
        ee.Image
            The image corresponding to the selected statistical metric,
            defined in `self.key`.
        """

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
        
        """
        Generates the seasonal image for the custom period p1 of the specified year.

        This function filters the input image collection using the date range defined
        in `self.p17` for the given year. It computes statistical summaries
        (maximum, median, mean, 90th percentile, and 95th percentile) over the images
        in that period. For optical sensors, the selected index function is applied
        before aggregation; for radar data (e.g., Sentinel-1), raw values are used.

        Parameters
        ----------
        y : int or str
            Year for which the custom period p1 image should be computed.

        Returns
        -------
        ee.Image
            The image corresponding to the selected statistical metric,
            defined in `self.key`.
        """

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
        
        """
        Generates the seasonal image for the custom period p1 of the specified year.

        This function filters the input image collection using the date range defined
        in `self.p18` for the given year. It computes statistical summaries
        (maximum, median, mean, 90th percentile, and 95th percentile) over the images
        in that period. For optical sensors, the selected index function is applied
        before aggregation; for radar data (e.g., Sentinel-1), raw values are used.

        Parameters
        ----------
        y : int or str
            Year for which the custom period p1 image should be computed.

        Returns
        -------
        ee.Image
            The image corresponding to the selected statistical metric,
            defined in `self.key`.
        """

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
        
        """
        Generates the seasonal image for the custom period p1 of the specified year.

        This function filters the input image collection using the date range defined
        in `self.p19` for the given year. It computes statistical summaries
        (maximum, median, mean, 90th percentile, and 95th percentile) over the images
        in that period. For optical sensors, the selected index function is applied
        before aggregation; for radar data (e.g., Sentinel-1), raw values are used.

        Parameters
        ----------
        y : int or str
            Year for which the custom period p1 image should be computed.

        Returns
        -------
        ee.Image
            The image corresponding to the selected statistical metric,
            defined in `self.key`.
        """

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
        
        """
        Generates the seasonal image for the custom period p1 of the specified year.

        This function filters the input image collection using the date range defined
        in `self.p20` for the given year. It computes statistical summaries
        (maximum, median, mean, 90th percentile, and 95th percentile) over the images
        in that period. For optical sensors, the selected index function is applied
        before aggregation; for radar data (e.g., Sentinel-1), raw values are used.

        Parameters
        ----------
        y : int or str
            Year for which the custom period p1 image should be computed.

        Returns
        -------
        ee.Image
            The image corresponding to the selected statistical metric,
            defined in `self.key`.
        """

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
            
        """
        Generates the seasonal image for the custom period p1 of the specified year.

        This function filters the input image collection using the date range defined
        in `self.p21` for the given year. It computes statistical summaries
        (maximum, median, mean, 90th percentile, and 95th percentile) over the images
        in that period. For optical sensors, the selected index function is applied
        before aggregation; for radar data (e.g., Sentinel-1), raw values are used.

        Parameters
        ----------
        y : int or str
            Year for which the custom period p1 image should be computed.

        Returns
        -------
        ee.Image
            The image corresponding to the selected statistical metric,
            defined in `self.key`.
        """

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
        
        """
        Generates the seasonal image for the custom period p1 of the specified year.

        This function filters the input image collection using the date range defined
        in `self.p22` for the given year. It computes statistical summaries
        (maximum, median, mean, 90th percentile, and 95th percentile) over the images
        in that period. For optical sensors, the selected index function is applied
        before aggregation; for radar data (e.g., Sentinel-1), raw values are used.

        Parameters
        ----------
        y : int or str
            Year for which the custom period p1 image should be computed.

        Returns
        -------
        ee.Image
            The image corresponding to the selected statistical metric,
            defined in `self.key`.
        """

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
        
        """
        Generates the seasonal image for the custom period p1 of the specified year.

        This function filters the input image collection using the date range defined
        in `self.p23` for the given year. It computes statistical summaries
        (maximum, median, mean, 90th percentile, and 95th percentile) over the images
        in that period. For optical sensors, the selected index function is applied
        before aggregation; for radar data (e.g., Sentinel-1), raw values are used.

        Parameters
        ----------
        y : int or str
            Year for which the custom period p1 image should be computed.

        Returns
        -------
        ee.Image
            The image corresponding to the selected statistical metric,
            defined in `self.key`.
        """

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
        
        """
        Generates the seasonal image for the custom period p1 of the specified year.

        This function filters the input image collection using the date range defined
        in `self.p24` for the given year. It computes statistical summaries
        (maximum, median, mean, 90th percentile, and 95th percentile) over the images
        in that period. For optical sensors, the selected index function is applied
        before aggregation; for radar data (e.g., Sentinel-1), raw values are used.

        Parameters
        ----------
        y : int or str
            Year for which the custom period p1 image should be computed.

        Returns
        -------
        ee.Image
            The image corresponding to the selected statistical metric,
            defined in `self.key`.
        """

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
        
        """
        Generate composite NDVI or SAR images for each year across specified time periods.
        
        This method creates temporal composites by combining imagery from different periods
        within each year of the specified date range. It supports seasonal (4 periods),
        monthly (12 periods), and bi-monthly (24 periods) compositing strategies.
        
        The method handles both optical satellite data (NDVI-based) and SAR data (Sentinel-1),
        with support for different statistical aggregations including mean, 90th percentile,
        and 95th percentile values.
        
        Returns
        -------
        ee.ImageCollection
            An Earth Engine ImageCollection containing composite images for each year.
            Each image in the collection represents one year of data with bands corresponding
            to the temporal periods (seasons, months, or bi-monthly periods).
        
        Notes
        -----
        The method behavior depends on several instance attributes:
        
        - `self.periods` : {4, 12, 24}
            Number of temporal periods to composite:
            * 4: Seasonal composites (winter, spring, summer, autumn)
            * 12: Monthly composites (january through december)
            * 24: Bi-monthly composites (p1 through p24)
        
        - `self.sat` : str
            Satellite type. If 'S1', uses SAR band naming convention (VH),
            otherwise uses optical band naming (nd for NDVI)
        
        - `self.key` : str
            Statistical measure to extract:
            * Default: mean values
            * 'perc_90': 90th percentile values take
            * 'perc_95': 95th percentile values
        
        - `self.start_year`, `self.end_year` : int
            Year range for composite generation (inclusive start, exclusive end)
        
        - `self.roi` : ee.Geometry
            Region of interest for clipping the composites
        
        - `self.imagelist` : list
            Instance list that gets populated with composite images
        
        Band Naming Convention
        ----------------------
        Optical data (non-S1):
            - Mean: 'nd', 'nd_1', 'nd_2', ...
            - 90th percentile: 'nd_p90', 'nd_p90_1', 'nd_p90_2', ...
            - 95th percentile: 'nd_p95', 'nd_p95_1', 'nd_p95_2', ...
        
        SAR data (S1):
            - Mean: 'VH', 'VH_1', 'VH_2', ...
            - 90th percentile: 'VH_p90', 'VH_p90_1', 'VH_p90_2', ...
            - 95th percentile: 'VH_p95', 'VH_p95_1', 'VH_p95_2', ...
        
        Output Band Names
        -----------------
        4 periods: ['winter', 'spring', 'summer', 'autumn']
        12 periods: ['january', 'february', ..., 'december']
        24 periods: ['p1', 'p2', ..., 'p24']
        
        Examples
        --------
        >>> # Assuming an instance with periods=4, sat='L8', key='mean'
        >>> composite_collection = instance.get_year_composite()
        >>> # Returns ImageCollection with seasonal composites for each year
        
        >>> # For SAR data with 90th percentile
        >>> # instance.sat = 'S1', instance.key = 'perc_90'
        >>> sar_composites = instance.get_year_composite()
        >>> # Returns ImageCollection with VH 90th percentile seasonal composites
        
        Raises
        ------
        ValueError
            If `self.periods` is not 4, 12, or 24
        
        Warning
        -------
        Prints warning message if `self.key` is not recognized (not in ['perc_90', 'perc_95'])
        
        Dependencies
        ------------
        Requires the following instance methods to be implemented:
        - For 4 periods: get_winter(), get_spring(), get_summer(), get_autumn()
        - For 12 periods: get_january() through get_december()
        - For 24 periods: get_p1() through get_p24()
        
        See Also
        --------
        The method assumes Earth Engine (ee) is initialized and available.
        """
        
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
            
        """
        Exports a single Earth Engine image as a GeoTIFF file to the local filesystem.

        This method is useful to export seasonal or custom composites (e.g., median NDVI for Africa
        from 2001–2020) created within the pipeline. It leverages `geemap.ee_export_image` to save
        the image locally in GeoTIFF format.

        Parameters
        ----------
        image : ee.Image
            The Earth Engine image to export.
        name : str, optional
            Filename of the output image, including the extension (default is 'mycomposition.tif').
        crs : str, optional
            Coordinate Reference System in EPSG format (e.g., 'EPSG:4326', 'EPSG:32629').
        scale : int, optional
            Pixel resolution of the output image, in meters. Defaults to 10, appropriate for
            Sentinel-2 data. For Landsat data, use 30. Can also be reduced to avoid
            export size limitations.

        Returns
        -------
        None
            The image is saved to the current working directory as a GeoTIFF file with four bands
            (e.g., max, median, or 90th percentile).

        Notes
        -----
        The export region is defined by `self.roi`.
        """
        
        filename = os.path.join(os.getcwd(), name)
        geemap.ee_export_image(image, filename=filename, scale=scale, crs=crs, region=self.roi, file_per_band=False) 

        print('Image have been exported')
        
        
    def get_export(self, crs='EPSG:4326', scale=10):
        
        """
        Exports yearly NDVI composites from the current image collection as GeoTIFF files.

        This method generates a yearly composite for each year in the range defined
        by `self.start_year` and `self.end_year`, and exports each image to the local
        working directory using `geemap.ee_export_image`.

        Parameters
        ----------
        crs : str, optional
            Coordinate Reference System in EPSG format (e.g., 'EPSG:4326', 'EPSG:32629').
        scale : int, optional
            Pixel resolution of the output images, in meters. Defaults to 10 (suitable for
            Sentinel-2). For Landsat, use 30. Can also be reduced to work around export
            size limitations when downloading large areas.

        Returns
        -------
        None
            The method saves one GeoTIFF per year in the current working directory,
            named as 'ndvi_<year>.tif'. Each file includes 4 bands (e.g., max, median,
            mean, 90th percentile).

        Notes
        -----
        The export region is defined by `self.roi`.
        """

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

        """
        Exports a seasonal NDVI composite animation as a GIF file.

        This method generates a GIF animation from yearly NDVI composites using the selected
        seasonal bands (e.g., winter, spring, summer). The animation is created with
        `geemap.download_ee_video` and annotated using `geemap.add_text_to_gif`.

        Parameters
        ----------
        name : str, optional
            Name of the output GIF file (default is 'mygif.gif'). The file will be saved
            to the current working directory.
        bands : list of str, optional
            List of seasonal or statistical band names to include in the RGB channels
            of the GIF (e.g., ['winter', 'spring', 'summer']).

        Returns
        -------
        None
            Two files are saved:
            - `<name>.gif`: The unannotated composite animation.
            - `<name>_texted.gif`: The same animation with year text overlay.

        Notes
        -----
        - Visualization parameters are automatically adapted based on the satellite type.
        For SAR data, dynamic percentiles are computed; for optical data, standard values
        are used (`min=0.15`, `max=0.85`).
        - The region used for export is defined by `self.roi`.
        """
        
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


    def export_with_fishnet(self, image, name_prefix='composite', scale=10, crs='EPSG:4326'):
    
        """
        Exports an image by dividing the region of interest (ROI) into adaptive tiles using a fishnet grid.

        This method splits the ROI into smaller tiles based on pixel resolution and exports each tile
        individually as a GeoTIFF file. Useful for handling large areas or avoiding export size limits
        in Earth Engine.

        Parameters
        ----------
        image : ee.Image
            The Earth Engine image to export.
        name_prefix : str, optional
            Prefix used for naming each exported tile (default is 'composite').
        scale : int, optional
            Pixel resolution of the output tiles in meters. Default is 10 (suitable for Sentinel-2).
        crs : str, optional
            Coordinate Reference System in EPSG format (e.g., 'EPSG:4326').

        Returns
        -------
        None
            Each tile is exported as a separate GeoTIFF file, saved to the current working directory.

        Notes
        -----
        - The method dynamically adjusts tile size: ~25 km² for high-resolution (<=10 m) and ~50 km² for lower resolution.
        - Only tiles that intersect the ROI are exported.
        - Uses `geemap.ee_export_image` for exporting.
        """

        import math

        # Determinar tamaño del tile según la resolución
        tile_km = 50 if scale >= 30 else 25
        tile_m = tile_km * 1000

        bounds = self.roi.bounds()
        coords = bounds.coordinates().getInfo()[0]
        xmin, ymin = coords[0]
        xmax, ymax = coords[2]

        # Calcular pasos en longitud/latitud aproximando metros
        x_steps = math.ceil((xmax - xmin) * 111320 / tile_m)
        y_steps = math.ceil((ymax - ymin) * 110540 / tile_m)

        tile_id = 0
        for i in range(x_steps):
            for j in range(y_steps):
                x0 = xmin + (i * tile_m / 111320)
                y0 = ymin + (j * tile_m / 110540)
                x1 = x0 + tile_m / 111320
                y1 = y0 + tile_m / 110540
                cell = ee.Geometry.Rectangle([x0, y0, x1, y1])

                # Solo exportar si hay intersección
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


    def get_stats(self, image, geom=None, name=None, stat='MEDIAN', scale=10, to_file=False):

        """
        Computes zonal statistics for a given image and geometry, and returns the result as a GeoDataFrame.

        This method supports both local file paths (e.g., shapefiles, GeoJSON) and Earth Engine geometry
        objects. It calculates summary statistics (e.g., mean, median) over the specified regions using
        `geemap.zonal_statistics`.

        Parameters
        ----------
        image : ee.Image
            The Earth Engine image from which to compute zonal statistics.
        geom : str or ee.FeatureCollection or ee.Geometry, optional
            The geometry to use for zonal statistics. Can be:
            - A file path to a `.shp` or `.geojson` file
            - An Earth Engine `FeatureCollection` or `Geometry`
            - `None` to use `self.roi`
        name : str, optional
            Output filename prefix for saving the results (default is 'zonal_stats').
            Used only if `to_file=True`.
        stat : str, optional
            Type of statistic to compute (e.g., 'MEAN', 'MEDIAN', 'MIN', 'MAX').
        scale : int, optional
            Pixel resolution in meters used for sampling (default is 10).
        to_file : bool, optional
            If True, the result will be saved as a shapefile in the current working directory.

        Returns
        -------
        geopandas.GeoDataFrame
            A GeoDataFrame containing the computed zonal statistics for each feature.

        Raises
        ------
        ValueError
            If the input geometry path is not a `.shp` or `.geojson` file.

        Notes
        -----
        - This method uses `geemap.zonal_statistics` under the hood.
        - Exported shapefiles include the computed statistics and geometry.
        """

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
            roi = geom.geometry() if hasattr(geom, 'geometry') else geom

        if name is None:
            name = 'zonal_stats'

        # Ruta temporal si se desea exportar
        out_shp = os.path.join(os.getcwd(), name + '.shp')

        # Ejecutar análisis y capturar como GeoDataFrame
        gdf = geemap.zonal_statistics(
            image=image,
            roi=roi,
            statistics_type=stat,
            scale=scale,
            return_fc=True  # <- Devuelve FeatureCollection
        )

        gdf = gpd.GeoDataFrame.from_features(gdf.getInfo()['features'])

        if to_file:
            gdf.to_file(out_shp)
            print(f'Saved as {out_shp}')

        return gdf