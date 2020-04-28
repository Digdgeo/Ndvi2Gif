
class ndvi_seasonality(object):
                               
    def __init__(self, roi=None, start_year=2016, end_year=2020, sat='Sentinel'):
        
        '''Start the instance with the parameters for start and end year. Also we set the fixed values for
        the seasons dates and generate the collection as ndvi 32 days merge of all landsats'''
        
        self.roi = roi
        if self.roi is None:
            self.roi = ee.Geometry.Polygon(
                [[[-6.766047, 36.776586], 
                  [-6.766047, 37.202186], 
                  [-5.867729, 37.202186], 
                  [-5.867729, 36.776586], 
                  [-6.766047, 36.776586]]], None, False)
            
        else:

            if not isinstance(roi, ee.Geometry):

                try:
                    self.roi = self.roi.geometry()
                except Exception as e:
                    print('Could not convert the provided roi to ee.Geometry')
                    print(e)
                    return

                
        self.start_year = start_year
        self.end_year = end_year
        sat_list = ['Sentinel', 'Landsat']
        self.sat = sat
        #self.name = name
        #self.crs = crs
        self.imagelist = []
        #self.out_gif = os.path.join(os.getcwd(), self.name)
        # Here we define the periods, feel free to change in case your are looking for seasons
        self.winter = ['-01-01', '-03-31']
        self.spring = ['-04-01', '-06-30']
        self.summer = ['-07-01', '-09-30']
        self.autumn = ['-10-01', '-12-31']

        # Here we defined the collections to generate the ndvi
        # Just need Red and NIR bandas, so avoid and rename them to apply the ndvi formula... 
        # ...to all sensors in the same way
        # Get Landsat surface reflectance collections for OLI, ETM+ and TM sensors.
        
        LC08col = ee.ImageCollection("LANDSAT/LC08/C01/T1_SR").select(['B4', 'B5'], ['Red', 'Nir'])
        LE07col = ee.ImageCollection("LANDSAT/LE07/C01/T1_SR").select(['B3', 'B4'], ['Red', 'Nir'])
        LT05col = ee.ImageCollection("LANDSAT/LT05/C01/T1_SR").select(['B3', 'B4'], ['Red', 'Nir'])
        LT04col = ee.ImageCollection("LANDSAT/LT04/C01/T1_SR").select(['B3', 'B4'], ['Red', 'Nir'])
        
        # Notice that S2 is not in Surface Reflectance but in TOA, this because otherwise only... 
        # ...had data starting in 2017. Using BOA we have 2 extra years, starting at 2015
        # We use the wide 8 band instead of the 8A narrower badn, that is also more similar to landsat NIR
        # But this way we have NDVI at 10 m instead of 20. And we are using TOA instead of SR, so who cares?
        
        S2col = ee.ImageCollection("COPERNICUS/S2").select(['B4', 'B8'], ['Red', 'Nir'])
        
        # Set the collection that will be used
        if self.sat == 'Sentinel':
            self.ndvi_col = S2col
        elif self.sat == 'Landsat':
            self.ndvi_col = LC08col.merge(LE07col).merge(LT05col).merge(LT04col)
        else: 
            print('You should choose between Sentinel and Landsat data')
            pass
    
    def get_ndvi(self, image):
        '''Here we apply the NDVI calculation'''   
        return image.normalizedDifference(['Nir', 'Red'])
        
    
    def get_winter(self, y):
        
        '''Here comes the funny thing. We have to generate the winter image for each year'''
        init, ends = str(y) + self.winter[0], str(y) + self.winter[1]
        return self.ndvi_col.filterDate(init, ends).map(self.get_ndvi).max()
        
    def get_spring(self, y):
        
        '''Here comes the funny thing. We have to generate the spring image for each year'''
        init, ends = str(y) + self.spring[0], str(y) + self.spring[1]
        return self.ndvi_col.filterDate(init, ends).map(self.get_ndvi).max()
        
    def get_summer(self, y):
        
        '''Here comes the funny thing. We have to generate the summer image for each year'''
        init, ends = str(y) + self.summer[0], str(y) + self.summer[1]
        return self.ndvi_col.filterDate(init, ends).map(self.get_ndvi).max()
        
    def get_autumn(self, y):
        
        '''Here comes the funny thing. We have to generate the autumn image for each year'''
        init, ends = str(y) + self.autumn[0], str(y) + self.autumn[1]
        return self.ndvi_col.filterDate(init, ends).map(self.get_ndvi).max()
    
    def get_year_composite(self):
        
        '''return the composite ndvi for each year'''
        for y in range(self.start_year, self.end_year):
            
            composite = ee.Image.cat(self.get_winter(y), self.get_spring(y), self.get_summer(y), self.get_autumn(y)).clip(self.roi)

            compositer = composite.select(['nd', 'nd_1', 'nd_2', 'nd_3'], ['winter', 'spring', 'summer', 'autumn'])
            #self.col = ee.ImageCollection(compositer)
            self.imagelist.append(compositer)
        
        
        ndvi_comp_coll = ee.ImageCollection.fromImages(self.imagelist)
        #ndvi_comp_coll_masked = ndvi_comp_coll.gt(0.1) We could apply masks in case we need it
        
        return ndvi_comp_coll
    
    def get_gif(self, name='mygif.gif', bands=['winter', 'spring', 'summer']):
        
        out_gif = os.path.join(os.getcwd(), name)
        self.imagelist = []
        video_args = {
          'dimensions': 768,
          'region': self.roi, 
          'framesPerSecond': 10,
          'bands': bands, 
          'min': 0.1,
          'max': 0.8,
          'gamma': [1, 1, 1]
        }
        
        geemap.download_ee_video(self.get_year_composite(), video_args, out_gif)
        texted_gif = out_gif[:-4] + '_texted.gif'
        geemap.add_text_to_gif(out_gif, texted_gif, xy=('5%', '90%'), text_sequence=self.start_year, font_size=30, font_color='#ffffff', add_progress_bar=False)
        
        
    def get_export(self, crs='EPSG:4326'):
        
        
        if len(self.imagelist) == 0:
            self.get_year_composite()
            
        count = (len(self.imagelist))
            
        #return geemap.ee_export_image_collection(self.get_year_composite(), os.getcwd())
        for n in range(count):
            year = self.start_year + n
            #print(self.imagelist[n])
            image = self.imagelist[n]
            name = 'ndvi_' + str(year) + '.tif'
            filename = os.path.join(os.getcwd(), name)
            print('Exporting {}'.format(filename), '\n')
            geemap.ee_export_image(image, filename=filename, scale=10,
                            crs=crs, region=self.roi, file_per_band=False)
        print('All the images in the ndvi collection have been exported')
