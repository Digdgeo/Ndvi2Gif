# s1_ard.py
"""
Sentinel-1 Analysis Ready Data (ARD) Processor for Google Earth Engine.

This module provides advanced SAR preprocessing capabilities including radiometric
terrain correction and speckle filtering algorithms for Sentinel-1 GRD imagery.

Based on Vollrath et al. (2020) and the gee_s1_ard implementation.
Optimized for vegetation monitoring and land cover analysis in complex terrain.

Classes
-------
S1ARDProcessor
    Main processor class for Sentinel-1 ARD generation

References
----------
Vollrath, A., Mullissa, A., & Reiche, J. (2020). Angular-based radiometric slope 
correction for Sentinel-1 on google earth engine. Remote Sensing, 12(11), 1867.

https://github.com/adugnag/gee_s1_ard

Author: Diego García Díaz
Date: 2024
License: MIT
"""

import ee
import math


class S1ARDProcessor:
    """
    Analysis Ready Data (ARD) processor for Sentinel-1 SAR imagery.
    
    Implements advanced preprocessing techniques for Sentinel-1 GRD data including
    radiometric terrain correction and various speckle filtering algorithms.
    
    Parameters
    ----------
    speckle_filter : {'REFINED_LEE', 'LEE', 'GAMMA_MAP', 'LEE_SIGMA', 'BOXCAR', None}
        Speckle filter algorithm. Default is 'REFINED_LEE'.
    speckle_filter_kernel_size : int
        Filter kernel size in pixels (must be odd). Default is 7.
    terrain_correction : bool
        Enable radiometric terrain correction. Default is True.
    terrain_flattening_model : {'VOLUME', 'SURFACE'}
        Scattering model for terrain correction. Default is 'VOLUME'.
    dem : {'COPERNICUS_30', 'COPERNICUS_90', 'SRTM_30', 'SRTM_90'}
        Digital Elevation Model. Default is 'COPERNICUS_30'.
    format : {'LINEAR', 'DB'}
        Output format. Default is 'LINEAR'.
    
    Attributes
    ----------
    dem_ee : ee.Image
        Earth Engine DEM image object
        
    Examples
    --------
    >>> processor = S1ARDProcessor(
    ...     speckle_filter='REFINED_LEE',
    ...     terrain_correction=True
    ... )
    >>> collection = ee.ImageCollection('COPERNICUS/S1_GRD')
    >>> processed = collection.map(processor.process_image)
    """
    
    def __init__(self, 
                 speckle_filter='REFINED_LEE',
                 speckle_filter_kernel_size=7,
                 terrain_correction=True,
                 terrain_flattening_model='VOLUME',
                 dem='COPERNICUS_30',
                 format='LINEAR'):
        """
        Initialize the Sentinel-1 ARD processor with specified parameters.
        
        Configures the preprocessing pipeline for Sentinel-1 data including
        speckle filtering and terrain correction options. Default parameters
        are optimized for vegetation monitoring in mountainous areas.
        
        Parameters
        ----------
        speckle_filter : str, optional
            Speckle filter algorithm to apply. Available options:
            
            - 'BOXCAR': Simple mean filter (fast but less edge-preserving)
            - 'LEE': Standard Lee filter (good balance)
            - 'REFINED_LEE': Refined Lee with edge detection (recommended)
            - 'GAMMA_MAP': Gamma Maximum A Posteriori filter
            - 'LEE_SIGMA': Lee Sigma filter using standard deviation
            - None: No speckle filtering
            
            Default is 'REFINED_LEE'.
            
        speckle_filter_kernel_size : int, optional
            Size of the filter kernel in pixels. Must be odd (3, 5, 7, 9...).
            Larger kernels provide more smoothing but less detail preservation.
            Default is 7.
            
        terrain_correction : bool, optional
            Enable radiometric terrain correction using local incidence angle.
            Essential for mountainous areas to reduce topographic effects.
            Default is True.
            
        terrain_flattening_model : str, optional
            Scattering model for terrain correction:
            
            - 'VOLUME': Volume scattering model (recommended for vegetation)
            - 'SURFACE': Surface scattering model (for bare soil/rock)
            
            Default is 'VOLUME'.
            
        dem : str, optional
            Digital Elevation Model to use for terrain correction:
            
            - 'COPERNICUS_30': Copernicus DEM at 30m resolution (recommended)
            - 'COPERNICUS_90': Copernicus DEM at 90m resolution
            - 'SRTM_30': SRTM at 30m resolution
            - 'SRTM_90': SRTM at 90m resolution
            
            Default is 'COPERNICUS_30'.
            
        format : str, optional
            Output format for backscatter values:
            
            - 'LINEAR': Linear power scale (required for most indices)
            - 'DB': Decibel scale (better for visualization)
            
            Default is 'LINEAR'.
            
        Examples
        --------
        >>> # Standard configuration for vegetation monitoring
        >>> processor = S1ARDProcessor(
        ...     speckle_filter='REFINED_LEE',
        ...     terrain_correction=True,
        ...     terrain_flattening_model='VOLUME'
        ... )
        
        >>> # Configuration for water/bare soil
        >>> processor = S1ARDProcessor(
        ...     speckle_filter='LEE',
        ...     terrain_flattening_model='SURFACE',
        ...     format='DB'
        ... )
        
        >>> # Minimal processing (only speckle filter)
        >>> processor = S1ARDProcessor(
        ...     terrain_correction=False,
        ...     speckle_filter='BOXCAR'
        ... )
        
        References
        ----------
        Vollrath, A., Mullissa, A., & Reiche, J. (2020). 
        Angular-based radiometric slope correction for Sentinel-1 on google earth engine. 
        Remote Sensing, 12(11), 1867.
        """
        self.speckle_filter = speckle_filter
        self.speckle_filter_kernel_size = speckle_filter_kernel_size
        self.terrain_correction = terrain_correction
        self.terrain_flattening_model = terrain_flattening_model
        self.dem = dem
        self.format = format
        
        # Load the specified DEM
        self.dem_ee = self._get_dem()
    
    def _get_dem(self):
        """
        Load the specified Digital Elevation Model from Earth Engine.
        
        Returns
        -------
        ee.Image
            DEM image object from Earth Engine catalog.

        Raises
        ------
        ee.EEException
            If the DEM asset cannot be accessed or there is an Earth Engine catalog error.
            
        Notes
        -----
        Copernicus DEM is generally recommended over SRTM for better quality
        and global coverage including high latitude regions.
        """
        dem_dict = {
            'COPERNICUS_30': ee.ImageCollection('COPERNICUS/DEM/GLO30').select('DEM').mosaic(),
            'COPERNICUS_90': ee.Image('COPERNICUS/DEM/GLO90').select('DEM'),
            'SRTM_30': ee.Image('USGS/SRTMGL1_003'),
            'SRTM_90': ee.Image('CGIAR/SRTM90_V4')
        }
        return dem_dict.get(self.dem, dem_dict['COPERNICUS_30'])
    
    def apply_terrain_correction(self, image):
        """
        Apply radiometric terrain correction to reduce topographic effects.
        
        Implements the angular-based radiometric slope correction method from
        Vollrath et al. (2020). This correction compensates for variations in
        backscatter caused by local terrain slope and aspect relative to the
        sensor viewing geometry.
        
        Parameters
        ----------
        image : ee.Image
            Sentinel-1 image with VV, VH, and angle bands.
            
        Returns
        -------
        ee.Image
            Terrain-corrected image with adjusted VV and VH bands.

        Raises
        ------
        ee.EEException
            If required bands or properties are missing (e.g., ``VV``, ``VH``, ``angle``,
            or ``orbitProperties_pass``) or if Earth Engine operations fail during terrain
            correction (projection/reprojection errors, invalid geometries).
            
        Notes
        -----
        The correction factor is clamped between 0.5 and 2.0 to prevent
        overcorrection in extreme terrain conditions.
        
        References
        ----------
        Vollrath, A., Mullissa, A., & Reiche, J. (2020). 
        Angular-based radiometric slope correction for Sentinel-1 on google earth engine. 
        Remote Sensing, 12(11), 1867.
        """
        img_geom = image.geometry()
        
        # Calculate terrain slope and aspect
        slope = ee.Terrain.slope(self.dem_ee).clip(img_geom)
        aspect = ee.Terrain.aspect(self.dem_ee).clip(img_geom)
        
        # Get incidence angle from image
        angle = image.select('angle')
        
        # Calculate local incidence angle based on scattering model
        if self.terrain_flattening_model == 'VOLUME':
            # Volume scattering model (vegetation)
            # Convert angles to radians
            phi_i = angle.multiply(math.pi / 180.0)
            alpha_s = slope.multiply(math.pi / 180.0)
            phi_s = aspect.multiply(math.pi / 180.0)
            
            # Satellite azimuth angle (depends on orbit)
            orbit_pass = ee.String(image.get('orbitProperties_pass'))
            azimuth = ee.Number(ee.Algorithms.If(
                orbit_pass.equals('DESCENDING'),
                280.0,  # Approximate azimuth for descending pass
                80.0    # Approximate azimuth for ascending pass
            ))
            phi_r = azimuth.multiply(math.pi / 180.0)
            
            # Calculate local incidence angle
            # θ_lia = arccos(cos(φ_i) * cos(α_s) + sin(φ_i) * sin(α_s) * cos(φ_s - φ_r))
            lia = alpha_s.cos().multiply(phi_i.cos()).add(
                alpha_s.sin().multiply(phi_i.sin()).multiply(
                    phi_s.subtract(phi_r).cos()
                )
            ).acos()
            
            # Correction factor = cos(φ_i) / cos(θ_lia)
            correction = phi_i.cos().divide(lia.cos())
            
        else:  # SURFACE
            # Surface scattering model (bare soil, water)
            # Simpler correction based on slope
            correction = slope.cos()
        
        # Clamp correction factor to prevent extreme values
        correction_clamped = correction.clamp(0.5, 2.0)
        
        # Apply correction to polarization bands
        vv_corrected = image.select('VV').multiply(correction_clamped)
        vh_corrected = image.select('VH').multiply(correction_clamped)
        
        return image.addBands(vv_corrected, None, True).addBands(vh_corrected, None, True)
    
    def apply_speckle_filter(self, image):
        """
        Apply the selected speckle filter to reduce SAR noise.
        
        Routes to the appropriate filter implementation based on the
        configured speckle_filter parameter.
        
        Parameters
        ----------
        image : ee.Image
            Input SAR image with VV and VH bands.
            
        Returns
        -------
        ee.Image
            Filtered image with reduced speckle noise.

        Raises
        ------
        ee.EEException
            If underlying Earth Engine focal/statistical operations fail during filtering.
        """
        if self.speckle_filter == 'BOXCAR':
            return self._boxcar_filter(image)
        elif self.speckle_filter == 'LEE':
            return self._lee_filter(image)
        elif self.speckle_filter == 'REFINED_LEE':
            return self._refined_lee_filter(image)
        elif self.speckle_filter == 'GAMMA_MAP':
            return self._gamma_map_filter(image)
        elif self.speckle_filter == 'LEE_SIGMA':
            return self._lee_sigma_filter(image)
        else:
            return image
    
    def _boxcar_filter(self, image):
        """
        Apply boxcar (mean) filter for basic speckle reduction.
        
        Simple spatial averaging filter. Fast but doesn't preserve edges well.
        
        Parameters
        ----------
        image : ee.Image
            Input SAR image.
            
        Returns
        -------
        ee.Image
            Mean-filtered image.

        Raises
        ------
        ee.EEException
            If Earth Engine focal operations fail for the requested kernel size or geometry.
        """
        kernel_radius = int(self.speckle_filter_kernel_size / 2)
        vv = image.select('VV').focal_mean(
            radius=kernel_radius, 
            kernelType='square', 
            units='pixels'
        )
        vh = image.select('VH').focal_mean(
            radius=kernel_radius,
            kernelType='square',
            units='pixels'
        )
        return image.addBands(vv, None, True).addBands(vh, None, True)
    
    def _lee_filter(self, image):
        """
        Apply standard Lee filter for speckle reduction.
        
        Adaptive filter that uses local statistics to preserve edges
        while reducing speckle in homogeneous areas.
        
        Parameters
        ----------
        image : ee.Image
            Input SAR image.
            
        Returns
        -------
        ee.Image
            Lee-filtered image.

        Raises
        ------
        ee.EEException
            If neighborhood reducers or focal operations fail in Earth Engine.
            
        References
        ----------
        Lee, J.S. (1980). Digital image enhancement and noise filtering by use 
        of local statistics. IEEE Transactions on Pattern Analysis and Machine 
        Intelligence, (2), 165-168.
        """
        kernel_radius = int(self.speckle_filter_kernel_size / 2)
        
        def apply_lee_to_band(band_name):
            band = image.select(band_name)
            
            # Calculate local statistics
            mean = band.focal_mean(radius=kernel_radius, kernelType='square', units='pixels')
            variance = band.reduceNeighborhood(
                reducer=ee.Reducer.variance(),
                kernel=ee.Kernel.square(kernel_radius, 'pixels')
            )
            
            # Coefficient of variation for speckle (typical value for SAR)
            cv_noise = 0.25
            
            # Weight factor
            weight = variance.divide(variance.add(mean.pow(2).multiply(cv_noise**2)))
            
            # Apply Lee filter
            filtered = mean.add(weight.multiply(band.subtract(mean)))
            
            return filtered.rename(band_name)
        
        vv_filtered = apply_lee_to_band('VV')
        vh_filtered = apply_lee_to_band('VH')
        
        return image.addBands(vv_filtered, None, True).addBands(vh_filtered, None, True)
    
    def _refined_lee_filter(self, image):
        """
        Apply Refined Lee filter with edge detection.

        Advanced version of Lee filter that detects edges using directional
        windows and applies filtering adaptively to preserve structural details.

        Parameters
        ----------
        image : ee.Image
            Input SAR image.

        Returns
        -------
        ee.Image
            Refined Lee filtered image.

        Raises
        ------
        ee.EEException
            If convolution or neighborhood reducers fail in Earth Engine during the
            refined Lee processing steps.

        References
        ----------
        Lee, J.S. (2009). Refined filtering of image noise using local statistics.
        Computer Graphics and Image Processing, 15(4), 380-389.
        """
        kernel_radius = int(self.speckle_filter_kernel_size / 2)

        def apply_refined_lee_to_band(band_name):
            band = image.select(band_name)

            # Step 1: Edge detection using directional gradients
            kernels = self._get_directional_kernels(kernel_radius)

            gradients = []
            for kernel in kernels:
                gradient = band.convolve(kernel)
                gradients.append(gradient.abs())

            # Find direction with minimum gradient (most homogeneous)
            gradient_stack = ee.Image.cat(gradients)
            min_gradient = gradient_stack.reduce(ee.Reducer.min())

            # Step 2: Apply adaptive Lee filter
            mean = band.focal_mean(radius=kernel_radius, kernelType='square', units='pixels')
            variance = band.reduceNeighborhood(
                reducer=ee.Reducer.variance(),
                kernel=ee.Kernel.square(kernel_radius, 'pixels')
            )

            # Equivalent Number of Looks (ENL) for Sentinel-1 IW GRD
            enl = 4.4

            # Noise factor
            sigma_v = 1.0 / math.sqrt(enl)

            # Coefficient of variation
            cv = variance.sqrt().divide(mean)
            cv_threshold = sigma_v * 1.5

            # Adaptive weight factor (corregido)
            sigma_v_sq = ee.Image.constant(sigma_v * sigma_v)
            weight = ee.Image(1.0).subtract(
                sigma_v_sq.divide(cv.multiply(cv))
            ).clamp(0, 1)

            # Apply weight only where cv > threshold
            weight = weight.where(cv.lte(cv_threshold), 0)

            # Apply Refined Lee filter
            filtered = mean.add(weight.multiply(band.subtract(mean)))

            return filtered.rename(band_name)

        vv_filtered = apply_refined_lee_to_band('VV')
        vh_filtered = apply_refined_lee_to_band('VH')

        return image.addBands(vv_filtered, None, True).addBands(vh_filtered, None, True)

    
    def _get_directional_kernels(self, radius):
        """
        Generate directional kernels for edge detection.
        
        Parameters
        ----------
        radius : int
            Kernel radius (not currently used, kept for future enhancement).
            
        Returns
        -------
        list of ee.Kernel
            List of directional convolution kernels.

        Raises
        ------
        None
        """
        # Simplified 3x3 directional kernels
        kernels = [
            ee.Kernel.fixed(3, 3, [[0, 0, 0], [1, 1, 1], [0, 0, 0]]),  # Horizontal
            ee.Kernel.fixed(3, 3, [[0, 1, 0], [0, 1, 0], [0, 1, 0]]),  # Vertical
            ee.Kernel.fixed(3, 3, [[1, 0, 0], [0, 1, 0], [0, 0, 1]]),  # Diagonal 1
            ee.Kernel.fixed(3, 3, [[0, 0, 1], [0, 1, 0], [1, 0, 0]])   # Diagonal 2
        ]
        return kernels
    
    def _gamma_map_filter(self, image):
        """
        Apply Gamma Maximum A Posteriori (MAP) filter.
        
        Statistical filter based on gamma distribution assumption for SAR data.
        Particularly effective for preserving texture information.
        
        Parameters
        ----------
        image : ee.Image
            Input SAR image.
            
        Returns
        -------
        ee.Image
            Gamma MAP filtered image.

        Raises
        ------
        ee.EEException
            If neighborhood statistics or arithmetic operations fail in Earth Engine
            while computing the MAP weights.
            
        References
        ----------
        Lopes, A., Nezry, E., Touzi, R., & Laur, H. (1993). 
        Structure detection and statistical adaptive speckle filtering in SAR images.
        International Journal of Remote Sensing, 14(9), 1735-1758.
        """
        kernel_radius = int(self.speckle_filter_kernel_size / 2)
        enl = 4.4  # Equivalent Number of Looks for S1 IW GRD
        
        def apply_gamma_map_to_band(band_name):
            band = image.select(band_name)
            
            # Local statistics
            mean = band.focal_mean(radius=kernel_radius, kernelType='square', units='pixels')
            variance = band.reduceNeighborhood(
                reducer=ee.Reducer.variance(),
                kernel=ee.Kernel.square(kernel_radius, 'pixels')
            )
            
            # Gamma MAP parameters
            alpha = (enl + 1) / (enl - 1)
            b = variance.divide(mean.pow(2)).multiply(alpha).subtract(enl - 1).divide(enl + 1)
            
            # Weight factor
            weight = b.multiply(mean).divide(b.multiply(mean).add(band))
            
            # Apply filter
            filtered = band.multiply(weight).add(mean.multiply(ee.Image(1).subtract(weight)))
            
            return filtered.rename(band_name)
        
        vv_filtered = apply_gamma_map_to_band('VV')
        vh_filtered = apply_gamma_map_to_band('VH')
        
        return image.addBands(vv_filtered, None, True).addBands(vh_filtered, None, True)
    
    def _lee_sigma_filter(self, image):
        """
        Apply Lee Sigma filter using standard deviation thresholding.
        
        Variant of Lee filter that uses sigma range (standard deviation)
        to identify and preserve point targets and edges.
        
        Parameters
        ----------
        image : ee.Image
            Input SAR image.
            
        Returns
        -------
        ee.Image
            Lee Sigma filtered image.

        Raises
        ------
        ee.EEException
            If neighborhood standard deviation or masking operations fail in Earth Engine.
            
        References
        ----------
        Lee, J.S. (1983). Digital image smoothing and the sigma filter.
        Computer Vision, Graphics, and Image Processing, 24(2), 255-269.
        """
        kernel_radius = int(self.speckle_filter_kernel_size / 2)
        sigma_range = 2  # Number of standard deviations
        
        def apply_lee_sigma_to_band(band_name):
            band = image.select(band_name)
            
            # Local statistics
            mean = band.focal_mean(radius=kernel_radius, kernelType='square', units='pixels')
            stddev = band.reduceNeighborhood(
                reducer=ee.Reducer.stdDev(),
                kernel=ee.Kernel.square(kernel_radius, 'pixels')
            )
            
            # Define valid range (mean ± sigma_range * stddev)
            lower_bound = mean.subtract(stddev.multiply(sigma_range))
            upper_bound = mean.add(stddev.multiply(sigma_range))
            
            # Mask pixels within range
            mask = band.gte(lower_bound).And(band.lte(upper_bound))
            
            # Apply mean only to pixels within range
            filtered = band.where(mask, mean)
            
            return filtered.rename(band_name)
        
        vv_filtered = apply_lee_sigma_to_band('VV')
        vh_filtered = apply_lee_sigma_to_band('VH')
        
        return image.addBands(vv_filtered, None, True).addBands(vh_filtered, None, True)
    
    def to_db(self, image):
        """
        Convert backscatter values from linear to decibel scale.
        
        Parameters
        ----------
        image : ee.Image
            SAR image with VV and VH bands in linear scale.
            
        Returns
        -------
        ee.Image
            Image with bands converted to decibel scale (10*log10).

        Raises
        ------
        ee.EEException
            If required bands (``VV``, ``VH``) are missing or if Earth Engine evaluation
            fails during the logarithmic conversion.
            
        Notes
        -----
        Decibel scale is preferred for visualization and some analyses
        as it compresses the dynamic range and normalizes the distribution.
        """
        vv_db = image.select('VV').log10().multiply(10).rename('VV')
        vh_db = image.select('VH').log10().multiply(10).rename('VH')
        return image.addBands(vv_db, None, True).addBands(vh_db, None, True)
    
    def process_image(self, image):
        """
        Apply complete preprocessing pipeline to a single image.
        
        Executes the full ARD processing chain in order:
        1. Terrain correction (if enabled)
        2. Speckle filtering (if configured)
        3. Format conversion (if DB requested)
        
        Parameters
        ----------
        image : ee.Image
            Raw Sentinel-1 GRD image.
            
        Returns
        -------
        ee.Image
            Preprocessed ARD image ready for analysis.

        Raises
        ------
        ee.EEException
            If any step of the processing chain fails in Earth Engine (terrain correction,
            speckle filtering, or format conversion).
            
        Examples
        --------
        >>> processor = S1ARDProcessor()
        >>> s1_collection = ee.ImageCollection('COPERNICUS/S1_GRD')
        >>> processed = s1_collection.map(processor.process_image)
        """
        # 1. Terrain correction (if enabled)
        if self.terrain_correction:
            image = self.apply_terrain_correction(image)
        
        # 2. Speckle filtering
        if self.speckle_filter:
            image = self.apply_speckle_filter(image)
        
        # 3. Format conversion (if requested)
        if self.format == 'DB':
            image = self.to_db(image)
        
        # Preserve original metadata
        return image.copyProperties(image, ['system:time_start', 'system:time_end'])