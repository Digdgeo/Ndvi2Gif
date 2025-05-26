import os
import ee
import geemap
import requests
import zipfile
import geopandas as gpd
import fiona
from geemap import zonal_statistics
from io import BytesIO
import calendar
from datetime import datetime, timedelta

def scale_OLI(image):
    opticalBands = image.select(['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7']).multiply(0.0000275).add(-0.2).rename(['Blue', 'Green', 'Red', 'Nir', 'Swir1', 'Swir2'])
    return image.addBands(opticalBands, None, True)
    
def scale_ETM(image):
    opticalBands = image.select(['SR_B1','SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7']).multiply(0.0000275).add(-0.2).rename(['Blue', 'Green', 'Red', 'Nir', 'Swir1', 'Swir2'])
    return image.addBands(opticalBands, None, True)

class NdviSeasonality:
    """
    Generate remote sensing index seasonal composition GIFs and images with dynamic period generation.
    
    This refactored version eliminates code duplication by dynamically generating periods
    and using a single generic function for all temporal composites.
    """
    
    def __init__(self, roi=None, periods=4, start_year=2016, end_year=2020, sat='S2', key='max', index='ndvi'):
        print('There we go again...')
        
        # Initialize ROI (same as original)
        self.roi = roi
        if self.roi is None:
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
                print('Invalid ROI path format')
        else:
            if not isinstance(roi, ee.Geometry):
                try:
                    self.roi = self.roi.geometry()
                except Exception as e:
                    print('Could not convert the provided roi to ee.Geometry')
                    print(e)
                    return
        
        # Set parameters
        self.periods = periods
        self.start_year = start_year
        self.end_year = end_year
        self.sat = sat if sat in ['S2', 'Landsat', 'MODIS', 'S1'] else 'S2'
        self.key = key if key in ['max', 'median', 'perc_90', 'perc_95', 'mean'] else 'max'
        self.imagelist = []
        self.index = index
        
        # Index calculation methods
        self.d = {
            'ndvi': self.get_ndvi, 'ndwi': self.get_ndwi, 'mndwi': self.get_mndwi, 
            'evi': self.get_evi, 'savi': self.get_savi, 'gndvi': self.get_gndvi, 
            'avi': self.get_avi, 'nbri': self.get_nbri, 'ndsi': self.get_ndsi, 
            'aweinsh': self.get_aweinsh, 'awei': self.get_awei, 'ndmi': self.get_ndmi
        }
        
        # **DYNAMIC PERIOD GENERATION** - This replaces all the hardcoded periods!
        self.period_dates, self.period_names = self._generate_periods(periods)
        
        # Initialize satellite collections (same as original)
        self._setup_satellite_collections()
    
    def _generate_periods(self, n_periods):
        """
        Dynamically generate period dates and names based on the number of periods.
        
        Parameters
        ----------
        n_periods : int
            Number of periods to divide the year into (4, 12, 24, or any number)
        
        Returns
        -------
        tuple
            (period_dates, period_names) where:
            - period_dates: list of [start_date, end_date] pairs
            - period_names: list of period names
        """
        if n_periods == 4:
            # Traditional seasons
            period_dates = [
                ['-01-01', '-03-31'],  # Winter
                ['-04-01', '-06-30'],  # Spring  
                ['-07-01', '-09-30'],  # Summer
                ['-10-01', '-12-31']   # Autumn
            ]
            period_names = ['winter', 'spring', 'summer', 'autumn']
            
        elif n_periods == 12:
            # Monthly periods - using fixed days to avoid leap year issues
            period_dates = [
                ['-01-01', '-01-31'],  # January
                ['-02-01', '-02-28'],  # February (always 28 - works for all years)
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
            # Bi-monthly periods (every ~15 days) - fixed dates to avoid leap year issues
            period_dates = [
                ['-01-01', '-01-15'], ['-01-16', '-01-31'],  # January
                ['-02-01', '-02-15'], ['-02-16', '-02-28'],  # February (always 28)
                ['-03-01', '-03-15'], ['-03-16', '-03-31'],  # March
                ['-04-01', '-04-15'], ['-04-16', '-04-30'],  # April
                ['-05-01', '-05-15'], ['-05-16', '-05-31'],  # May
                ['-06-01', '-06-15'], ['-06-16', '-06-30'],  # June
                ['-07-01', '-07-15'], ['-07-16', '-07-31'],  # July
                ['-08-01', '-08-15'], ['-08-16', '-08-31'],  # August
                ['-09-01', '-09-15'], ['-09-16', '-09-30'],  # September
                ['-10-01', '-10-15'], ['-10-16', '-10-31'],  # October
                ['-11-01', '-11-15'], ['-11-16', '-11-30'],  # November
                ['-12-01', '-12-15'], ['-12-16', '-12-31']   # December
            ]
            period_names = [f'p{i+1}' for i in range(24)]
        
        else:
            # Generic periods - divide year equally using day-of-year approach
            period_dates = []
            period_names = []
            days_per_period = 365 // n_periods
            
            for i in range(n_periods):
                # Calculate start and end day of year
                start_day = i * days_per_period + 1
                if i == n_periods - 1:  # Last period goes to end of year
                    end_day = 365
                else:
                    end_day = (i + 1) * days_per_period
                
                # Convert day of year to month-day (using non-leap year)
                start_date = datetime(2021, 1, 1) + timedelta(days=start_day - 1)  # 2021 is not leap
                end_date = datetime(2021, 1, 1) + timedelta(days=end_day - 1)
                
                start_str = f'-{start_date.month:02d}-{start_date.day:02d}'
                end_str = f'-{end_date.month:02d}-{end_date.day:02d}'
                
                period_dates.append([start_str, end_str])
                period_names.append(f'p{i+1}')
        
        return period_dates, period_names
    
    def _setup_satellite_collections(self):
        """Setup satellite collections (same as original)"""
        # Landsat collections
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
        
        # Sentinel-2
        S2col = ee.ImageCollection("COPERNICUS/S2_HARMONIZED").select(
            ['B2', 'B3', 'B4', 'B8', 'B11', 'B12'], 
            ['Blue', 'Green', 'Red', 'Nir', 'Swir1', 'Swir2']
        ).filterBounds(self.roi)
        
        # MODIS
        MOD09Q1 = ee.ImageCollection("MODIS/061/MOD09A1").select(
            ['sur_refl_b03', 'sur_refl_b04', 'sur_refl_b01', 'sur_refl_b02', 'sur_refl_b06', 'sur_refl_b07'], 
            ['Blue', 'Green', 'Red', 'Nir', 'Swir1', 'Swir2']
        ).filterBounds(self.roi)
        
        # Sentinel-1
        s1 = ee.ImageCollection('COPERNICUS/S1_GRD').filter(
            ee.Filter.listContains('transmitterReceiverPolarisation', 'VH')
        ).filter(ee.Filter.eq('instrumentMode', 'IW'))
        s1Ascending = s1.filter(ee.Filter.eq('orbitProperties_pass', 'ASCENDING'))
        s1Descending = s1.filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))
        s1S1 = s1Ascending.select('VH').merge(s1Descending.select('VH')).filterBounds(self.roi)
        
        # Set the collection
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
    
    def get_period_composite(self, year, period_idx):
        """
        SINGLE GENERIC FUNCTION that replaces ALL the get_winter, get_january, get_p1, etc. functions!
        
        Parameters
        ----------
        year : int
            Year for the composite
        period_idx : int
            Index of the period (0 to n_periods-1)
        
        Returns
        -------
        ee.Image
            Composite image for the specified period and year
        """
        start_date, end_date = self.period_dates[period_idx]
        init = str(year) + start_date
        ends = str(year) + end_date
        
        # Dictionary to store results for all statistics
        period_stats = {}
        
        if self.sat != 'S1':
            # Optical satellites - apply index calculation first
            period_stats['max'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).max()
            period_stats['median'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).median()
            period_stats['mean'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).mean()
            period_stats['perc_90'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([90]))
            period_stats['perc_95'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([95]))
        else:
            # SAR satellite - use raw values
            period_stats['max'] = self.ndvi_col.filterDate(init, ends).max()
            period_stats['median'] = self.ndvi_col.filterDate(init, ends).median()
            period_stats['mean'] = self.ndvi_col.filterDate(init, ends).mean()
            period_stats['perc_90'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([90]))
            period_stats['perc_95'] = self.ndvi_col.filterDate(init, ends).reduce(ee.Reducer.percentile([95]))
        
        return period_stats[self.key]
    
    def get_year_composite(self):
        """
        DRAMATICALLY SIMPLIFIED version that works with any number of periods!
        
        Returns
        -------
        ee.ImageCollection
            Collection of yearly composites
        """
        # Generate band names dynamically
        if self.sat != 'S1':
            if self.key not in ['perc_90', 'perc_95']:
                base_bands = ['nd'] + [f'nd_{i}' for i in range(1, self.periods)]
            elif self.key == 'perc_90':
                base_bands = ['nd_p90'] + [f'nd_p90_{i}' for i in range(1, self.periods)]
            else:  # perc_95
                base_bands = ['nd_p95'] + [f'nd_p95_{i}' for i in range(1, self.periods)]
        else:
            if self.key not in ['perc_90', 'perc_95']:
                base_bands = ['VH'] + [f'VH_{i}' for i in range(1, self.periods)]
            elif self.key == 'perc_90':
                base_bands = ['VH_p90'] + [f'VH_p90_{i}' for i in range(1, self.periods)]
            else:  # perc_95
                base_bands = ['VH_p95'] + [f'VH_p95_{i}' for i in range(1, self.periods)]
        
        # Clear previous results
        self.imagelist = []
        
        # Generate composites for each year
        for year in range(self.start_year, self.end_year):
            # Get composites for all periods in this year
            period_images = []
            for period_idx in range(self.periods):
                period_composite = self.get_period_composite(year, period_idx)
                period_images.append(period_composite)
            
            # Combine all periods into a single multi-band image
            composite = ee.Image.cat(period_images).clip(self.roi)
            
            # Rename bands to meaningful names
            compositer = composite.select(base_bands, self.period_names)
            
            self.imagelist.append(compositer)
        
        return ee.ImageCollection.fromImages(self.imagelist)
    
    # Index calculation methods (same as original - keeping all of them)
    def get_ndvi(self, image):
        return image.normalizedDifference(['Nir', 'Red'])
    
    def get_ndwi(self, image):
        return image.normalizedDifference(['Green', 'Nir'])

    def get_mndwi(self, image):
        return image.normalizedDifference(['Green', 'Swir1'])
    
    def get_evi(self, image):
        return image.expression(
            '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))', {
            'NIR': image.select('Nir'),
            'RED': image.select('Red'),
            'BLUE': image.select('Blue')}).rename(['nd'])
    
    def get_savi(self, image, L=0.428):
        return image.expression(
            '((NIR - RED) / (NIR + RED + L) * (1 + L))', {
            'NIR': image.select('Nir'),
            'RED': image.select('Red'),
            'L': L}).rename(['nd'])
    
    def get_aweinsh(self, image):
        return image.expression(
            '4.0 * (GREEN - SWIR1) - 0.25 * NIR + 2.75 * SWIR2', {
            'NIR': image.select('Nir'),
            'GREEN': image.select('Green'),
            'SWIR1':image.select('Swir1'),
            'SWIR2':image.select('Swir2')}).rename(['nd'])
    
    def get_awei(self, image):
        return image.expression(
            ('BLUE + 2.5 * GREEN - 1.5 * (NIR + SWIR1) - 0.25 * SWIR2'), {
            'NIR': image.select('Nir'),
            'BLUE': image.select('Blue'),
            'GREEN': image.select('Green'),
            'SWIR1':image.select('Swir1'),
            'SWIR2':image.select('Swir2')}).rename(['nd'])

    def get_gndvi(self, image):
        return image.normalizedDifference(['Nir', 'Green'])
    
    def get_avi(self, image, L=0.428):
        return image.expression(
            '(NIR * (1.0 - RED) * (NIR - RED)) ** (1/3)', {
            'NIR': image.select('Nir'),
            'RED': image.select('Red')}).rename(['nd'])

    def get_nbri(self, image):
        return image.normalizedDifference(['Nir', 'Swir2'])
    
    def get_ndsi(self, image):
        return image.normalizedDifference(['Green', 'Swir1'])

    def get_ndmi(self, image):
        return image.normalizedDifference(['Nir', 'Swir1'])
    
    # Export methods (same as original)
    def get_export_single(self, image, name='mycomposition.tif', crs='EPSG:4326', scale=10):
        filename = os.path.join(os.getcwd(), name)
        geemap.ee_export_image(image, filename=filename, scale=scale, crs=crs, region=self.roi, file_per_band=False) 
        print('Image have been exported')
        
    def get_export(self, crs='EPSG:4326', scale=10):
        self.imagelist = []
        self.get_year_composite()
        count = len(self.imagelist)
        print(count)
        
        for n in range(count):
            year = self.start_year + n
            image = self.imagelist[n]
            name = 'ndvi_' + str(year) + '.tif'
            filename = os.path.join(os.getcwd(), name)
            print('Exporting {}'.format(filename), '\n')
            geemap.ee_export_image(image, filename=filename, scale=scale, crs=crs, region=self.roi, file_per_band=False) 
        print('All the images in the collection have been exported')

    def get_gif(self, name='mygif.gif', bands=None):
        if bands is None:
            bands = self.period_names[:3]  # Use first 3 periods by default
            
        self.imagelist = []
        self.get_year_composite()
        
        out_gif = os.path.join(os.getcwd(), name)
        
        if self.sat == 'S1':
            # For SAR, compute dynamic ranges (simplified for this example)
            video_args = {
                'dimensions': 768,
                'region': self.roi, 
                'framesPerSecond': 10,
                'bands': bands, 
                'min': -25,  # Typical SAR values
                'max': 0,
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
        import math
        tile_km = 50 if scale >= 30 else 25
        tile_m = tile_km * 1000
        bounds = self.roi.bounds()
        coords = bounds.coordinates().getInfo()[0]
        xmin, ymin = coords[0]
        xmax, ymax = coords[2]
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

        out_shp = os.path.join(os.getcwd(), name + '.shp')
        gdf = geemap.zonal_statistics(
            image=image,
            roi=roi,
            statistics_type=stat,
            scale=scale,
            return_fc=True
        )
        gdf = gpd.GeoDataFrame.from_features(gdf.getInfo()['features'])

        if to_file:
            gdf.to_file(out_shp)
            print(f'Saved as {out_shp}')
        return gdf