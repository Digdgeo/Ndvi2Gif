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
    
    def __init__(self, roi=None, periods=4, start_year=2016, end_year=2020, 
             sat='S2', key='max', index='ndvi', percentile=90, orbit='BOTH'):
        """
        Initialize NdviSeasonality object for temporal remote sensing analysis.
        
        Parameters
        ----------
        roi : ee.Geometry, str, or None, optional
            Region of interest. Can be ee.Geometry, path to shapefile/geojson,
            DEIMS ID, WRS path/row, or S2 tile ID. Default is None (uses Andalusia).
        periods : int, optional
            Number of temporal periods. Default is 4 (seasons).
        start_year : int, optional
            Start year for analysis. Default is 2016.
        end_year : int, optional
            End year for analysis. Default is 2020.
        sat : str, optional
            Satellite sensor. Options: 'S2', 'Landsat', 'MODIS', 'S1'. Default is 'S2'.
        key : str, optional
            Statistical reducer. Options: 'max', 'median', 'mean', 'percentile'. Default is 'max'.
        index : str, optional
            Spectral index to compute. Available indices depend on satellite. Default is 'ndvi'.
        percentile : int, optional
            Percentile value when key='percentile'. Default is 90.
        orbit : str, optional
            Sentinel-1 orbit direction. Options: 'BOTH', 'ASCENDING', 'DESCENDING'. 
            Default is 'BOTH'. Only used when sat='S1'.
            - 'BOTH': Use all available orbits (maximum data coverage)
            - 'DESCENDING': Use only descending orbits (better geometric consistency)
            - 'ASCENDING': Use only ascending orbits (better geometric consistency)
        
        Examples
        --------
        >>> # Default: use both orbits for maximum coverage
        >>> sar_gif = NdviSeasonality(sat='S1', index='vh')
        
        >>> # Use only descending orbits for consistency
        >>> sar_gif = NdviSeasonality(sat='S1', index='vh', orbit='DESCENDING')
        
        >>> # Optical satellites ignore orbit parameter
        >>> ndvi_gif = NdviSeasonality(sat='S2', index='ndvi', orbit='DESCENDING')  # orbit ignored
        """
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
            # Para ROIs dibujados o Features, convertir a Geometry
            if isinstance(self.roi, list) and len(self.roi) > 0:
                # Manejar listas de Features (como Map.draw_features o draw_last_feature)
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
                    # Usar ROI por defecto en lugar de return
                    self.roi = ee.Geometry.Polygon(
                        [[[-6.766047, 36.776586], 
                        [-6.766047, 37.202186], 
                        [-5.867729, 37.202186], 
                        [-5.867729, 36.776586], 
                        [-6.766047, 36.776586]]], None, False)
                    print("Using default ROI")
        
        # Set parameters
        self.periods = periods
        self.start_year = start_year
        self.end_year = end_year
        self.key = key if key in ['max', 'median', 'percentile', 'mean'] else 'max'
        self.percentile = percentile
        self.imagelist = []
        self.index = index
        
        # Validar parámetro orbit
        valid_orbits = ['BOTH', 'ASCENDING', 'DESCENDING']
        if orbit not in valid_orbits:
            raise ValueError(f"Orbit '{orbit}' is not valid. Available options are: {valid_orbits}")
        
        self.orbit = orbit
    
        # Advertencia si se usa orbit con sensores no-SAR
        if sat != 'S1' and orbit != 'BOTH':
            print(f"Warning: orbit parameter '{orbit}' is only used with Sentinel-1. Ignoring for {sat}.")
            
        # Índices básicos (todos los sensores ópticos)
        self.optical_indices = {
            'ndvi', 'ndwi', 'mndwi', 'evi', 'savi', 'gndvi', 'avi', 
            'nbri', 'ndsi', 'aweinsh', 'awei', 'ndmi', 'msi', 'nmi', 
            'ndti', 'cri1', 'cri2', 'lai', 'pri', 'wdrvi'
        }

        # Exclusivos S2 (Red Edge)
        self.s2_exclusive_indices = {
            'ireci', 'mcari', 'ndre', 'reip', 'psri', 'cire', 'mtci', 's2rep', 'ndci'
        }

        # Exclusivos S3 (OLCI)
        self.s3_exclusive_indices = {
            'oci', 'tsi', 'cdom', 'turbidity', 'spm', 'kd490', 'floating_algae', 
            'red_edge_position', 'fluorescence_height', 'water_leaving_reflectance'
        }

        # SAR (S1)
        self.s1_indices = {
            'rvi', 'vv', 'vh', 'vv_vh_ratio', 'dpsvi', 'rfdi', 'vsdi'
        }

        # Mapeo final SIMPLE
        self.sensor_indices = {
            'S2': self.optical_indices | self.s2_exclusive_indices,
            'Landsat': self.optical_indices,
            'MODIS': self.optical_indices,
            'S1': self.s1_indices,
            'S3': self.optical_indices | self.s3_exclusive_indices
        }
        
        # Validar satélite
        if sat not in self.sensor_indices:
            available_sats = list(self.sensor_indices.keys())
            raise ValueError(f"Satellite '{sat}' is not supported. Available satellites are: {available_sats}")
        
        self.sat = sat
        
        # Validar índice para el satélite seleccionado
        available_indices = self.sensor_indices[self.sat]
        if index not in available_indices:
            raise ValueError(
                f"Index '{index}' is not available for {sat}. "
                f"Available indices for {sat} are: {sorted(list(available_indices))}"
            )
        
        self.index = index
        
        # Diccionario COMPLETO de métodos de cálculo (todos los índices)
        self.d = {
        # Índices básicos ópticos
        'ndvi': self.get_ndvi, 'ndwi': self.get_ndwi, 'mndwi': self.get_mndwi, 
        'evi': self.get_evi, 'savi': self.get_savi, 'gndvi': self.get_gndvi, 
        'avi': self.get_avi, 'nbri': self.get_nbri, 'ndsi': self.get_ndsi, 
        'aweinsh': self.get_aweinsh, 'awei': self.get_awei, 'ndmi': self.get_ndmi,
        
        # Índices adicionales ópticos (compartidos)
        'msi': self.get_msi, 'nmi': self.get_nmi, 'ndti': self.get_ndti,
        'cri1': self.get_cri1, 'cri2': self.get_cri2, 'lai': self.get_lai, 
        'pri': self.get_pri, 'wdrvi': self.get_wdrvi,
        
        # Índices específicos S2 (Red Edge)
        'ireci': self.get_ireci, 'mcari': self.get_mcari, 'reip': self.get_reip,
        'psri': self.get_psri, 'ndre': self.get_ndre, 'cig': self.get_cig,
        'cire': self.get_cire, 'mtci': self.get_mtci, 's2rep': self.get_s2rep,
        'ndci': self.get_ndci,

        # Índices específicos S3 (OLCI 21 bandas)
        'oci': self.get_oci, 'tsi': self.get_tsi, 'cdom': self.get_cdom,
        'turbidity': self.get_turbidity, 'spm': self.get_spm, 'kd490': self.get_kd490,
        'floating_algae': self.get_floating_algae, 'red_edge_position': self.get_red_edge_position,
        'fluorescence_height': self.get_fluorescence_height, 
        'water_leaving_reflectance': self.get_water_leaving_reflectance,

        # Índices SAR
        'rvi': self.get_rvi, 'vv': self.get_vv, 'vh': self.get_vh, 
        'vv_vh_ratio': self.get_vv_vh_ratio, 'dpsvi': self.get_dpsvi,
        'rfdi': self.get_rfdi, 'vsdi': self.get_vsdi
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
        """Setup satellite collections - SIMPLE y DIRECTO"""
        
        # Landsat collections (igual que antes)
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
        
        # Sentinel-2 - Surface Reflectance con todas las bandas
        S2col = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED").select([
            'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B11', 'B12'
        ], [
            'Blue', 'Green', 'Red', 'Red_Edge1', 'Red_Edge2', 'Red_Edge3', 
            'Nir', 'Swir1', 'Swir2'
        ]).filterBounds(self.roi)
        
        # MODIS
        MOD09Q1 = ee.ImageCollection("MODIS/061/MOD09A1").select(
            ['sur_refl_b03', 'sur_refl_b04', 'sur_refl_b01', 'sur_refl_b02', 'sur_refl_b06', 'sur_refl_b07'], 
            ['Blue', 'Green', 'Red', 'Nir', 'Swir1', 'Swir2']
        ).filterBounds(self.roi)
        
        # Sentinel-1 - SAR con control de órbitas
        s1 = ee.ImageCollection('COPERNICUS/S1_GRD').filter(
            ee.Filter.listContains('transmitterReceiverPolarisation', 'VH')
        ).filter(
            ee.Filter.listContains('transmitterReceiverPolarisation', 'VV')
        ).filter(ee.Filter.eq('instrumentMode', 'IW'))
        
        # Aplicar filtro de órbita según parámetro del usuario
        if self.orbit == 'ASCENDING':
            s1 = s1.filter(ee.Filter.eq('orbitProperties_pass', 'ASCENDING'))
            print("Using Sentinel-1 ascending orbits only.")
        elif self.orbit == 'DESCENDING':
            s1 = s1.filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))
            print("Using Sentinel-1 descending orbits only.")
        else:  # orbit == 'BOTH'
            print("Using all Sentinel-1 orbits (ascending + descending).")
        
        def apply_speckle_filter(image):
            filtered = image.focal_median(radius=1, kernelType='square', units='pixels')
            return filtered.copyProperties(image, ['system:time_start', 'system:time_end'])
        
        s1S1 = s1.select(['VV', 'VH']).map(apply_speckle_filter).filterBounds(self.roi)
        
        # Sentinel-3 - OLCI Level-1B (21 bandas espectrales)
        S3col = ee.ImageCollection("COPERNICUS/S3/OLCI").select([
            'Oa01_radiance',  # 400nm - Violet
            'Oa02_radiance',  # 412.5nm - Blue
            'Oa03_radiance',  # 442.5nm - Blue2
            'Oa04_radiance',  # 490nm - Blue_Green
            'Oa05_radiance',  # 510nm - Green
            'Oa06_radiance',  # 560nm - Green2
            'Oa07_radiance',  # 620nm - Red
            'Oa08_radiance',  # 665nm - Red2
            'Oa09_radiance',  # 673.75nm - Red3
            'Oa10_radiance',  # 681.25nm - Red_Edge1
            'Oa11_radiance',  # 708.75nm - Red_Edge2
            'Oa12_radiance',  # 753.75nm - NIR
            'Oa16_radiance',  # 778.75nm - NIR2
            'Oa17_radiance',  # 865nm - NIR3
            'Oa18_radiance',  # 885nm - NIR4
            'Oa21_radiance'   # 1020nm - NIR5
        ], [
            'Violet', 'Blue', 'Blue2', 'Blue_Green', 'Green', 'Green2',
            'Red', 'Red2', 'Red3', 'Red_Edge1', 'Red_Edge2', 'NIR',
            'NIR2', 'NIR3', 'NIR4', 'NIR5'
        ]).filterBounds(self.roi)
        
        # Set the collection - SIMPLE if/elif chain
        if self.sat == 'S2':
            self.ndvi_col = S2col
        elif self.sat == 'Landsat':
            self.ndvi_col = Landsat
        elif self.sat == 'MODIS':
            self.ndvi_col = MOD09Q1
        elif self.sat == 'S1':
            self.ndvi_col = s1S1
        elif self.sat == 'S3':
            self.ndvi_col = S3col
        else: 
            print('Not a valid satellite')
    
    def get_period_composite(self, year, period_idx):
        """
        Función simplificada - ya no necesita validación porque se hace en __init__
        """
        start_date, end_date = self.period_dates[period_idx]
        init = str(year) + start_date
        ends = str(year) + end_date
        
        # Dictionary to store results for all statistics
        period_stats = {}
        
        # Aplicar el método del índice directamente (ya está validado)
        period_stats['max'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).max()
        period_stats['median'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).median()
        period_stats['mean'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).mean()
        period_stats['percentile'] = self.ndvi_col.filterDate(init, ends).map(self.d[self.index]).reduce(ee.Reducer.percentile([self.percentile]))
        
        return period_stats[self.key]
    
    def get_available_indices(self, satellite=None):
        """
        Método público para consultar índices disponibles
        """
        if satellite is None:
            satellite = self.sat
            
        if satellite not in self.sensor_indices:
            available_sats = list(self.sensor_indices.keys())
            raise ValueError(f"Satellite '{satellite}' is not supported. Available satellites are: {available_sats}")
        
        return sorted(list(self.sensor_indices[satellite]))
    
    def get_all_available_indices(self):
        """
        Devuelve todos los índices disponibles organizados por sensor
        """
        result = {}
        for sensor, indices in self.sensor_indices.items():
            result[sensor] = sorted(list(indices))
        return result
    
    def get_year_composite(self):
        """
        Generate temporal composite images for each year in the specified time range.
        
        Creates multi-band images where each band represents a temporal period 
        (seasons, months, etc.) using the specified statistical reducer (max, median, 
        mean, or percentile) and the selected spectral index.
        
        The function automatically handles:
        - Dynamic band naming based on satellite type and index
        - Flexible temporal periods (4, 12, 24, or custom)
        - Multiple statistical reducers with configurable percentiles
        - Proper band naming for both optical and SAR indices
        - Data availability validation for each period
        
        Returns
        -------
        ee.ImageCollection
            Collection of multi-band composite images, one per year. Each image
            contains bands named after the temporal periods (e.g., 'winter', 'spring'
            for seasonal composites, or 'january', 'february' for monthly composites).
            
        Examples
        --------
        >>> # Seasonal NDVI composites
        >>> ndvi_gif = NdviSeasonality(sat='S2', index='ndvi', periods=4)
        >>> collection = ndvi_gif.get_year_composite()
        
        >>> # Monthly SAR composites with 90th percentile
        >>> sar_gif = NdviSeasonality(sat='S1', index='vh', periods=12, 
        ...                          key='percentile', percentile=90)
        >>> collection = sar_gif.get_year_composite()
        """
        # Generate band names dynamically
        if self.sat != 'S1':
            # Optical satellites
            if self.key == 'percentile':
                base_bands = [f'nd_p{self.percentile}'] + [f'nd_p{self.percentile}_{i}' for i in range(1, self.periods)]  # DINÁMICO
            else:
                base_bands = ['nd'] + [f'nd_{i}' for i in range(1, self.periods)]
        else:
            # SAR satellite - CORREGIDO para todos los índices SAR
            if self.index == 'vv':
                band_prefix = 'VV'
            elif self.index == 'vh':
                band_prefix = 'VH'
            elif self.index == 'rvi':
                band_prefix = 'RVI'
            elif self.index == 'vv_vh_ratio':
                band_prefix = 'RATIO'
            elif self.index == 'dpsvi':
                band_prefix = 'DPSVI'
            elif self.index == 'rfdi':
                band_prefix = 'RFDI'
            elif self.index == 'vsdi':
                band_prefix = 'VSDI'
            else:
                # Fallback para compatibilidad hacia atrás
                band_prefix = 'VH'
                print(f"Warning: Unknown SAR index '{self.index}', using VH as fallback")
            
            # Generar nombres de bandas con el prefijo correcto
            if self.key == 'percentile':
                base_bands = [f'{band_prefix}_p{self.percentile}'] + [f'{band_prefix}_p{self.percentile}_{i}' for i in range(1, self.periods)]
            else:
                base_bands = [band_prefix] + [f'{band_prefix}_{i}' for i in range(1, self.periods)]
        
        # Clear previous results
        self.imagelist = []
        
        # Generate composites for each year
        for year in range(self.start_year, self.end_year):
            # Get composites for all periods in this year
            period_images = []
            successful_periods = 0  # Track how many periods actually have data
            
            for period_idx in range(self.periods):
                try:
                    period_composite = self.get_period_composite(year, period_idx)
                    # Check if the period has any data by trying to get band names
                    band_count = period_composite.bandNames().size()
                    
                    # Only add if there's actually data
                    if band_count.getInfo() > 0:
                        period_images.append(period_composite)
                        successful_periods += 1
                    else:
                        print(f"No data for period {period_idx + 1} in year {year}")
                        break
                except Exception as e:
                    print(f"Error processing period {period_idx + 1} in year {year}: {str(e)}")
                    break
            
            # Only proceed if we have at least one period with data
            if successful_periods > 0:
                # Combine all available periods into a single multi-band image
                composite = ee.Image.cat(period_images).clip(self.roi)
                
                # Adjust band names and period names to match actual data
                actual_base_bands = base_bands[:successful_periods]
                actual_period_names = self.period_names[:successful_periods]
                
                # Rename bands to meaningful names
                compositer = composite.select(actual_base_bands, actual_period_names)
                
                self.imagelist.append(compositer)
                
                print(f"Year {year}: Successfully processed {successful_periods} periods using {self.index} index")
            else:
                print(f"Year {year}: No data available, skipping")
        
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
    
    def get_msi(self, image):
        """Moisture Stress Index"""
        return image.expression(
            'SWIR1 / NIR', {
            'NIR': image.select('Nir'),
            'SWIR1': image.select('Swir1')
        }).rename(['nd'])

    def get_nmi(self, image):
        """Normalized Multi-band Drought Index"""
        return image.expression(
            '(NIR - (SWIR1 + SWIR2)) / (NIR + (SWIR1 + SWIR2))', {
            'NIR': image.select('Nir'),
            'SWIR1': image.select('Swir1'),
            'SWIR2': image.select('Swir2')
        }).rename(['nd'])

    def get_ndti(self, image):
        """Normalized Difference Tillage Index"""
        return image.normalizedDifference(['Swir1', 'Swir2'])

    def get_cri1(self, image):
        """Carotenoid Reflectance Index 1"""
        return image.expression(
            '(1 / BLUE) - (1 / GREEN)', {
            'BLUE': image.select('Blue'),
            'GREEN': image.select('Green')
        }).rename(['nd'])

    def get_cri2(self, image):
        """Carotenoid Reflectance Index 2"""
        return image.expression(
            '(1 / BLUE) - (1 / RED)', {
            'BLUE': image.select('Blue'),
            'RED': image.select('Red')
        }).rename(['nd'])

    def get_lai(self, image):
        """Leaf Area Index approximation"""
        return image.expression(
            '3.618 * EVI - 0.118', {
            'EVI': self.get_evi(image).select('nd')
        }).rename(['nd'])

    def get_pri(self, image):
        """Photochemical Reflectance Index"""
        return image.normalizedDifference(['Green', 'Blue'])

    def get_wdrvi(self, image):
        """Wide Dynamic Range Vegetation Index"""
        return image.expression(
            '(0.1 * NIR - RED) / (0.1 * NIR + RED)', {
            'NIR': image.select('Nir'),
            'RED': image.select('Red')
        }).rename(['nd'])

    def get_cig(self, image):
        """Chlorophyll Index Green"""
        return image.expression(
            '(Nir / Green) - 1', {
            'Green': image.select('Green'),
            'Nir': image.select('Nir')
        }).rename(['nd'])

    # ÍNDICES SENTINEL-2 CON RED EDGE
    def get_ireci(self, image):
        """Inverted Red-Edge Chlorophyll Index - Muy sensible a clorofila"""
        return image.expression(
            '(Red_Edge3 - Red) / (Red_Edge1 / Red_Edge2)', {
            'Red': image.select('Red'),
            'Red_Edge1': image.select('Red_Edge1'),    # B5
            'Red_Edge2': image.select('Red_Edge2'),    # B6
            'Red_Edge3': image.select('Red_Edge3')     # B7
        }).rename(['nd'])

    def get_mcari(self, image):
        """Modified Chlorophyll Absorption Ratio Index"""
        return image.expression(
            '((Red_Edge1 - Red) - 0.2 * (Red_Edge1 - Green)) * (Red_Edge1 / Red)', {
            'Green': image.select('Green'),
            'Red': image.select('Red'), 
            'Red_Edge1': image.select('Red_Edge1')     # B5
        }).rename(['nd'])

    def get_ndre(self, image):
        """Normalized Difference Red Edge - Sensible a contenido de clorofila"""
        return image.normalizedDifference(['Nir', 'Red_Edge1'])  # NIR - Red Edge 1

    def get_reip(self, image):
        """Red Edge Inflection Point - Aproximación"""
        return image.expression(
            '700 + 40 * ((((Red + Red_Edge3) / 2) - Red_Edge1) / (Red_Edge2 - Red_Edge1))', {
            'Red': image.select('Red'),
            'Red_Edge1': image.select('Red_Edge1'),    # B5
            'Red_Edge2': image.select('Red_Edge2'),    # B6
            'Red_Edge3': image.select('Red_Edge3')     # B7
        }).rename(['nd'])

    def get_psri(self, image):
        """Plant Senescence Reflectance Index"""
        return image.expression(
            '(Red - Blue) / Red_Edge2', {
            'Blue': image.select('Blue'),
            'Red': image.select('Red'),
            'Red_Edge2': image.select('Red_Edge2')     # B6
        }).rename(['nd'])

    def get_cire(self, image):
        """Chlorophyll Index Red Edge"""
        return image.expression(
            '(Nir / Red_Edge1) - 1', {
            'Red_Edge1': image.select('Red_Edge1'),    # B5
            'Nir': image.select('Nir')
        }).rename(['nd'])

    def get_mtci(self, image):
        """MERIS Terrestrial Chlorophyll Index"""
        return image.expression(
            '(Red_Edge2 - Red_Edge1) / (Red_Edge1 - Red)', {
            'Red': image.select('Red'),
            'Red_Edge1': image.select('Red_Edge1'),    # B5
            'Red_Edge2': image.select('Red_Edge2')     # B6
        }).rename(['nd'])

    def get_s2rep(self, image):
        """Sentinel-2 Red Edge Position - Simplificado"""
        return image.expression(
            '705 + 35 * ((((Red + Red_Edge3) / 2) - Red_Edge1) / (Red_Edge2 - Red_Edge1))', {
            'Red': image.select('Red'),
            'Red_Edge1': image.select('Red_Edge1'),    # B5
            'Red_Edge2': image.select('Red_Edge2'),    # B6
            'Red_Edge3': image.select('Red_Edge3')     # B7
        }).rename(['nd'])
    
    def get_ndci(self, image):
        """
        Normalized Difference Chlorophyll Index (NDCI)
        Especialmente útil para detectar cianobacterias y clorofila-a en agua
        Usa Red Edge 1 (B5) y Red (B4) - optimizado para Sentinel-2
        
        Fórmula: (Red_Edge1 - Red) / (Red_Edge1 + Red)
        
        Referencias:
        - Mishra & Mishra (2012) - Cyanobacteria detection
        - Gitelson et al. (2007) - Chlorophyll-a estimation
        """
        return image.normalizedDifference(['Red_Edge1', 'Red'])
    
    # =================================================================
    # ÍNDICES S3 IMPLEMENTADOS (todos basados en radiancia L1B)
    # =================================================================

    def get_oci(self, image):
        """OLCI Chlorophyll Index - Custom usando radiancia L1B"""
        return image.expression(
            '(Red_Edge1 - Red2) / (Red_Edge1 + Red2)', {
            'Red2': image.select('Red2'),              # Oa08 - 665nm
            'Red_Edge1': image.select('Red_Edge1')     # Oa10 - 681.25nm
        }).rename(['nd'])

    def get_tsi(self, image):
        """Trophic State Index - Estado trófico de agua"""
        return image.expression(
            '(Red2 - Blue_Green) / (NIR - Blue_Green)', {
            'Blue_Green': image.select('Blue_Green'),  # Oa04 - 490nm
            'Red2': image.select('Red2'),              # Oa08 - 665nm  
            'NIR': image.select('NIR')                 # Oa12 - 753.75nm
        }).rename(['nd'])

    def get_cdom(self, image):
        """Colored Dissolved Organic Matter Index"""
        return image.expression(
            'Blue / Blue_Green', {
            'Blue': image.select('Blue'),              # Oa02 - 412.5nm
            'Blue_Green': image.select('Blue_Green')   # Oa04 - 490nm
        }).rename(['nd'])

    def get_turbidity(self, image):
        """Water Turbidity Index usando bandas OLCI"""
        return image.expression(
            'Red2 / Blue_Green', {
            'Red2': image.select('Red2'),              # Oa08 - 665nm
            'Blue_Green': image.select('Blue_Green')   # Oa04 - 490nm
        }).rename(['nd'])

    def get_spm(self, image):
        """Suspended Particulate Matter - Materia particulada suspendida"""
        return image.expression(
            '(Red3 - Red2) / (Red3 + Red2)', {
            'Red2': image.select('Red2'),    # Oa08 - 665nm
            'Red3': image.select('Red3')     # Oa09 - 673.75nm
        }).rename(['nd'])

    def get_kd490(self, image):
        """Diffuse Attenuation Coefficient at 490nm"""
        return image.expression(
            'log(Blue2 / Blue_Green)', {
            'Blue2': image.select('Blue2'),            # Oa03 - 442.5nm
            'Blue_Green': image.select('Blue_Green')   # Oa04 - 490nm
        }).rename(['nd'])

    def get_floating_algae(self, image):
        """Floating Algae Index - Detección de algas flotantes"""
        return image.expression(
            '(NIR - Red_Edge2) / (NIR + Red_Edge2)', {
            'Red_Edge2': image.select('Red_Edge2'),    # Oa11 - 708.75nm
            'NIR': image.select('NIR')                 # Oa12 - 753.75nm
        }).rename(['nd'])

    def get_red_edge_position(self, image):
        """Red Edge Position optimizado para OLCI"""
        return image.expression(
            '681.25 + 27.5 * ((Red2 + Red_Edge2) / 2 - Red_Edge1) / (Red_Edge2 - Red_Edge1)', {
            'Red2': image.select('Red2'),              # Oa08 - 665nm
            'Red_Edge1': image.select('Red_Edge1'),    # Oa10 - 681.25nm
            'Red_Edge2': image.select('Red_Edge2')     # Oa11 - 708.75nm
        }).rename(['nd'])

    def get_fluorescence_height(self, image):
        """Chlorophyll Fluorescence Height"""
        return image.expression(
            'Red3 - (Red2 + (Red_Edge1 - Red2) * (673.75 - 665) / (681.25 - 665))', {
            'Red2': image.select('Red2'),              # Oa08 - 665nm
            'Red3': image.select('Red3'),              # Oa09 - 673.75nm
            'Red_Edge1': image.select('Red_Edge1')     # Oa10 - 681.25nm
        }).rename(['nd'])

    def get_water_leaving_reflectance(self, image):
        """Water Leaving Reflectance (simplified approximation)"""
        return image.expression(
            'Green / (Blue + Green + Red)', {
            'Blue': image.select('Blue'),      # Oa02 - 412.5nm
            'Green': image.select('Green'),    # Oa05 - 510nm
            'Red': image.select('Red')         # Oa07 - 620nm
        }).rename(['nd'])
        
    # Nuevos métodos para índices SAR
    def get_rvi(self, image):
        """Radar Vegetation Index - Más robusto que VH solo"""
        return image.expression(
            '4 * VH / (VV + VH)', {
            'VV': image.select('VV'),
            'VH': image.select('VH')}).rename(['RVI'])

    def get_vv(self, image):
        """VV polarization - Sensible a superficie rugosa"""
        return image.select('VV').rename(['VV'])

    def get_vh(self, image):
        """VH polarization - Sensible a estructura vegetal"""
        return image.select('VH').rename(['VH'])

    def get_vv_vh_ratio(self, image):
        """VV/VH ratio - Muy sensible a cambios estructurales (ideal para siega)"""
        return image.expression(
            'VV / VH', {
            'VV': image.select('VV'),
            'VH': image.select('VH')}).rename(['RATIO'])

    def get_dpsvi(self, image):
        """Dual-pol SAR Vegetation Index - Optimizado para vegetación densa"""
        return image.expression(
            '(VV - VH) / (VV + VH)', {
            'VV': image.select('VV'),
            'VH': image.select('VH')}).rename(['DPSVI'])
    
    def get_rfdi(self, image):
        """Radar Forest Degradation Index"""
        return image.expression(
            '(VV - VH) / VV', {
            'VV': image.select('VV'),
            'VH': image.select('VH')
        }).rename(['RFDI'])

    def get_vsdi(self, image):
        """Vegetation Scattering Diversity Index"""
        return image.expression(
            'sqrt((VV - VH) ** 2 + (VV + VH) ** 2)', {
            'VV': image.select('VV'),
            'VH': image.select('VH')
        }).rename(['VSDI'])

    
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