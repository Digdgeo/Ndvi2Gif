"""
classification.py - Land Cover Classification Module for ndvi2gif v1.0.0

This module provides supervised and unsupervised classification capabilities
using multi-temporal composite images from NdviSeasonality.

Author: Diego García Díaz
Date: 2024
License: MIT
"""

import ee
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Optional, List, Dict, Union, Tuple, Any
import geemap
import geopandas as gpd

class LandCoverClassifier:
    """
    Land cover classification workflow based on temporal NDVI composites.

    This class integrates seasonal NDVI metrics (from an
    ``NdviSeasonality`` instance) with supervised and unsupervised
    classification methods in Google Earth Engine. It supports feature
    stack generation, training data ingestion, model fitting, and accuracy
    assessment.

    Attributes
    ----------
    processor : NdviSeasonality
        First (or only) processor; ROI and dates are taken from it.
    processors : list of NdviSeasonality
        Every processor, one per sensor.
    multi_sensor : bool
        True when more than one sensor is combined. Band names then carry
        the sensor as prefix (``S2_ndvi_2023_january``).
    scale : float
        Pixel size in meters of the feature stack, used for sampling,
        normalization and export.
    crs : str or None
        CRS of the common grid when the sensors are resampled, else None.
    feature_stack : ee.Image or None
        Image containing stacked features (NDVI indices, temporal metrics).
    training_data : ee.FeatureCollection or None
        FeatureCollection with labeled training samples.
    validation_data : ee.FeatureCollection or None
        FeatureCollection with labeled validation samples.
    classifier : ee.Classifier or None
        Trained Earth Engine classifier.
    classified_image : ee.Image or None
        Output land cover classification map.
    accuracy_results : dict or None
        Accuracy metrics computed from validation data.
    roi : ee.Geometry
        Region of interest inherited from ``processor``.
    periods : int
        Number of temporal periods per year.
    start_year : int
        First year of the analysis.
    end_year : int
        Last year (inclusive) of the analysis.
    sat : str
        Satellite name used (e.g., 'S2', 'L8').
    """
    
    # Accepted values of the ``resample`` argument besides a scale in meters
    RESAMPLE_MODES = ('coarser', 'finer')

    def __init__(self, ndvi_seasonality_instance, resample=None, crs=None):
        """
        Initialize the classifier with one or more ``NdviSeasonality`` instances.

        Parameters
        ----------
        ndvi_seasonality_instance : NdviSeasonality or list of NdviSeasonality
            Processor providing the temporal composites, ROI and analysis
            configuration. Pass a list to combine several sensors in one
            feature stack (e.g. Sentinel-2 optical indices with Sentinel-1
            backscatter); each processor must use a different sensor. The
            ROI is taken from the first one.
        resample : {'coarser', 'finer'} or float, optional
            Common pixel size of the feature stack when the sensors have
            different resolutions. Required in that case, since the choice
            changes what the classification means:

            - ``'coarser'``: finer sensors are aggregated by their **mean**
              onto the grid of the coarsest one. Nothing is invented, but
              detail is lost (S2 + Landsat -> 30 m).
            - ``'finer'``: coarser sensors are **bilinearly interpolated**
              onto the grid of the finest one. The map keeps the finer
              detail, but the interpolated bands are only smoothed, they do
              not gain information (S2 + Landsat -> 10 m).
            - a number: pixel size in meters. Each sensor is aggregated or
              interpolated to it, whichever applies.

            With a single sensor the default ``None`` keeps its native
            resolution, as in previous versions.
        crs : str, optional
            CRS of the common grid when resampling (e.g. ``'EPSG:32650'``).
            Defaults to the UTM zone of the ROI centroid, so pixel sizes are
            true meters.

        Raises
        ------
        ValueError
            If two processors use the same sensor, if ``resample`` is not a
            valid option, or if the sensors have different resolutions and
            no ``resample`` was given.

        Notes
        -----
        The constructor inherits spatial and temporal parameters from the
        first processor:
        - ROI (region of interest)
        - Number of periods per year
        - Start and end years
        - Satellite identifier

        Native resolutions come from
        :meth:`NdviSeasonality._default_scale_for_sat` (10 m for S1/S2, 30 m
        for Landsat, 250 m for MODIS, 300 m for S3...).

        Examples
        --------
        Sentinel-2 indices plus Sentinel-1 backscatter, both at 10 m::

            >>> s2 = NdviSeasonality(roi=roi, sat='S2', periods=12,
            ...                      start_year=2023, end_year=2023)
            >>> s1 = NdviSeasonality(roi=roi, sat='S1', periods=12,
            ...                      start_year=2023, end_year=2023)
            >>> clf = LandCoverClassifier([s2, s1])
            >>> stack = clf.create_feature_stack(
            ...     indices={'S2': ['ndvi', 'mndwi'], 'S1': ['vh', 'rvi']})

        Sentinel-2 and Landsat on the 30 m Landsat grid::

            >>> clf = LandCoverClassifier([s2, landsat], resample='coarser')
        """
        if isinstance(ndvi_seasonality_instance, (list, tuple)):
            processors = list(ndvi_seasonality_instance)
        else:
            processors = [ndvi_seasonality_instance]
        if not processors:
            raise ValueError("Pass at least one NdviSeasonality instance")

        sats = [p.sat for p in processors]
        repeated = sorted({s for s in sats if sats.count(s) > 1})
        if repeated:
            raise ValueError(
                f"Each processor must use a different sensor; repeated: {repeated}. "
                "Put all the indices of a sensor in the same processor through "
                "create_feature_stack(indices={...})."
            )

        self.processors = processors
        self.processor = processors[0]
        self.multi_sensor = len(processors) > 1
        self.feature_stack = None
        self.training_data = None
        self.validation_data = None
        self.classifier = None
        self.classified_image = None
        self.accuracy_results = None
        self.class_property = 'class'
        self.algorithm = None

        # Inherit parameters
        self.roi = self.processor.roi
        self.periods = self.processor.periods
        self.start_year = self.processor.start_year
        self.end_year = self.processor.end_year
        self.sat = self.processor.sat

        # Pixel size of the stack and whether bands have to be put on a common grid
        self.native_scales = {p.sat: p._default_scale_for_sat() for p in processors}
        self.resample = resample
        self.scale = self._resolve_scale(resample)
        self.crs = crs
        needs_grid = resample is not None
        if needs_grid and self.crs is None:
            self.crs = self._utm_crs_for_roi()

        print(f"LandCoverClassifier initialized for {', '.join(sats)}")
        print(f"Period: {self.start_year}-{self.end_year}, {self.periods} periods/year")
        if self.multi_sensor or needs_grid:
            natives = ', '.join(f"{s} {m} m" for s, m in self.native_scales.items())
            grid = f" on {self.crs}" if needs_grid else ""
            print(f"Feature stack at {self.scale} m{grid} (native: {natives})")

    def _resolve_scale(self, resample):
        """Pixel size in meters of the feature stack for a ``resample`` option."""
        natives = sorted(set(self.native_scales.values()))

        if resample is None:
            if len(natives) > 1:
                listed = ', '.join(f"{s} {m} m" for s, m in self.native_scales.items())
                raise ValueError(
                    f"The sensors have different resolutions ({listed}). Choose "
                    "how to combine them with resample='coarser' (mean "
                    "aggregation to the coarsest grid), resample='finer' "
                    "(bilinear interpolation to the finest grid) or "
                    "resample=<meters>."
                )
            return natives[0]

        if resample == 'coarser':
            return natives[-1]
        if resample == 'finer':
            return natives[0]
        if isinstance(resample, (int, float)) and not isinstance(resample, bool) and resample > 0:
            return resample
        raise ValueError(
            f"resample must be one of {self.RESAMPLE_MODES}, a scale in meters, "
            f"or None; got {resample!r}"
        )

    def _utm_crs_for_roi(self):
        """EPSG code of the UTM zone containing the ROI centroid."""
        lon, lat = self.roi.centroid(1).coordinates().getInfo()
        zone = int((lon + 180) // 6) + 1
        return f"EPSG:{(32600 if lat >= 0 else 32700) + zone}"

    def _to_common_grid(self, image, native_scale):
        """
        Put a single-sensor stack on the common grid of the feature stack.

        Composites reduced from an ``ee.ImageCollection`` carry no fixed
        projection, so the sensor's native resolution is declared first.
        From there the image is aggregated by its mean when the target pixel
        is coarser, bilinearly interpolated when it is finer, and simply
        aligned otherwise.
        """
        if self.resample is None:
            return image

        target = ee.Projection(self.crs).atScale(self.scale)
        image = image.setDefaultProjection(crs=self.crs, scale=native_scale)

        if self.scale > native_scale:
            # Input pixels per output pixel, with margin for grid misalignment
            ratio = int(np.ceil(self.scale / native_scale)) + 1
            image = image.reduceResolution(
                reducer=ee.Reducer.mean(), maxPixels=min(ratio * ratio, 65535)
            )
        elif self.scale < native_scale:
            image = image.resample('bilinear')

        return image.reproject(target)

    def create_feature_stack(self,
                           indices: Union[List[str], Dict[str, List[str]]] = None,
                           include_statistics: bool = True,
                           normalize: bool = True) -> ee.Image:
        """
        Create multi-temporal feature stack for classification.
        
        Parameters
        ----------
        indices : list of str or dict, optional
            Indices to stack. With a single sensor, a list such as
            ``['ndvi', 'mndwi']``. With several sensors, a dict mapping each
            sensor to its indices, such as
            ``{'S2': ['ndvi', 'mndwi'], 'S1': ['vh', 'rvi']}``. If None, uses
            the current index of each processor.
        include_statistics : bool
            Add temporal statistics (mean, std, max, min)
        normalize : bool
            Normalize to [0,1] range
            
        Returns
        -------
        ee.Image
            Multi-band feature stack. Bands are named
            ``<index>_<year>_<period>`` with a single sensor and
            ``<sensor>_<index>_<year>_<period>`` with several.

        Raises
        ------
        ValueError
            If `indices` contains unsupported names, names sensors that were
            not passed to the classifier, or is a list with several sensors.
        ee.EEException
            If Earth Engine image processing fails when computing the stack.
        """
        print("Creating feature stack...")

        indices_by_sat = self._indices_by_sat(indices)

        sensor_stacks = []
        stat_groups = []

        for processor in self.processors:
            sat = processor.sat
            prefix = f"{sat}_" if self.multi_sensor else ""
            feature_bands = []

            # Process each index
            for idx_name in indices_by_sat[sat]:
                print(f"  Processing {prefix}{idx_name}...")

                # Temporarily set processor index
                original_index = processor.index
                processor.index = idx_name
                try:
                    feature_bands.extend(self._index_bands(processor, idx_name, prefix))
                finally:
                    # Restore original index
                    processor.index = original_index
                stat_groups.append(f"{prefix}{idx_name}")

            sensor_stack = ee.Image.cat(feature_bands)
            sensor_stacks.append(
                self._to_common_grid(sensor_stack, self.native_scales[sat])
            )

        # Create feature stack
        self.feature_stack = ee.Image.cat(sensor_stacks)
        
        # Add statistics if requested
        if include_statistics:
            print("  Adding temporal statistics...")
            stats = self._compute_statistics(stat_groups)
            self.feature_stack = self.feature_stack.addBands(stats)
        
        # Normalize if requested
        if normalize:
            print("  Normalizing features...")
            self.feature_stack = self._normalize_image(self.feature_stack)
        
        # Clip to ROI
        self.feature_stack = self.feature_stack.clip(self.roi)
        
        band_count = self.feature_stack.bandNames().size().getInfo()
        print(f"Feature stack ready: {band_count} bands")
        
        return self.feature_stack

    def _indices_by_sat(self, indices):
        """Normalize the ``indices`` argument to ``{sat: [index, ...]}`` and validate it."""
        sats = [p.sat for p in self.processors]

        if indices is None:
            indices_by_sat = {p.sat: [p.index] for p in self.processors}
        elif isinstance(indices, dict):
            unknown = sorted(set(indices) - set(sats))
            if unknown:
                raise ValueError(
                    f"No processor for sensor(s) {unknown}; the classifier has {sats}"
                )
            missing = [s for s in sats if s not in indices]
            if missing:
                raise ValueError(
                    f"Give the indices of every sensor; missing: {missing}"
                )
            indices_by_sat = {s: list(indices[s]) for s in sats}
        else:
            if self.multi_sensor:
                raise ValueError(
                    "With several sensors pass indices as a dict, e.g. "
                    "{'S2': ['ndvi'], 'S1': ['vh']}"
                )
            indices_by_sat = {self.sat: list(indices)}

        for processor in self.processors:
            sat = processor.sat
            if not indices_by_sat[sat]:
                raise ValueError(f"No indices given for {sat}")
            available = processor.sensor_indices[sat]
            invalid = [idx for idx in indices_by_sat[sat] if idx not in available]
            if invalid:
                raise ValueError(f"Invalid indices for {sat}: {invalid}")

        return indices_by_sat

    def _index_bands(self, processor, idx_name, prefix):
        """One band per year and period of an index, named after both."""
        feature_bands = []

        # Generate composites
        processor.get_year_composite()

        # get_year_composite skips the years without any image, so the
        # position of an image in the collection is not its year offset.
        # Pair each image with its year from the per-year scene counts,
        # which list every year in the same order
        years_with_data = [
            year for year, counts in processor.period_scene_counts.items()
            if sum(counts) > 0
        ]

        # Stack all years (end_year is inclusive, like get_year_composite)
        for year, year_image in zip(years_with_data, processor.imagelist):
            counts = processor.period_scene_counts[year]

            # Rename bands with descriptive names
            for period_idx, period_name in enumerate(processor.period_names):
                # An empty period comes back as a fully masked band, and a
                # single masked band masks every pixel on sampling and
                # classification, so it is left out of the stack
                if counts[period_idx] == 0:
                    print(f"    Skipping {prefix}{idx_name} {year} {period_name}: no images")
                    continue
                band = year_image.select(period_name)
                band_name = f"{prefix}{idx_name}_{year}_{period_name}"
                feature_bands.append(band.rename(band_name))

        return feature_bands

    def _compute_statistics(self, groups: List[str]) -> ee.Image:
        """
        Compute per-pixel temporal statistics of each index.

        Reduces all the year/period bands of an index to their mean,
        standard deviation, maximum and minimum.

        Parameters
        ----------
        groups : list of str
            Band name prefixes, one per index: ``'ndvi'`` with a single
            sensor, ``'S2_ndvi'`` with several.

        Returns
        -------
        ee.Image
            Four bands per group: ``<group>_mean``, ``<group>_std``,
            ``<group>_max`` and ``<group>_min``.
        """
        stats_bands = []
        
        for group in groups:
            # Select all bands for this index. The year digits keep 'vv' from
            # also matching the 'vv_vh_ratio' bands
            idx_pattern = f"{group}_[0-9]{{4}}_.*"
            idx_bands = self.feature_stack.select(idx_pattern)
            
            # Calculate statistics
            mean = idx_bands.reduce(ee.Reducer.mean()).rename(f"{group}_mean")
            std = idx_bands.reduce(ee.Reducer.stdDev()).rename(f"{group}_std")
            max_val = idx_bands.reduce(ee.Reducer.max()).rename(f"{group}_max")
            min_val = idx_bands.reduce(ee.Reducer.min()).rename(f"{group}_min")
            
            stats_bands.extend([mean, std, max_val, min_val])
        
        return ee.Image.cat(stats_bands)
    
    def _normalize_image(self, image: ee.Image) -> ee.Image:
        """
        Normalize image bands to the [0, 1] range.

        Uses provided statistics (min, max) to scale each band,
        applying (value - min) / (max - min).

        Parameters
        ----------
        image : ee.Image
            Earth Engine image to be normalized.
        stats : dict
            Dictionary of per-band statistics (min, max).

        Returns
        -------
        ee.Image
            Normalized image with values in [0, 1].
        """
        # Get min/max per band
        minMax = image.reduceRegion(
            reducer=ee.Reducer.minMax(),
            geometry=self.roi,
            scale=self.scale,
            maxPixels=1e9,
            bestEffort=True
        )
        
        # Function to normalize a band
        def normalize_band(band_name):
            """
            Normalize a single band to the [0, 1] range.

            Applies (value - min) / (max - min) using band-specific
            statistics from ``Reducer.minMax()``.

            Parameters
            ----------
            band_name : str or ee.String
                Name of the band to normalize.

            Returns
            -------
            ee.Image
                Single-band image with values scaled to [0, 1].

            Notes
            -----
            If min == max for the band, a unit range is used to avoid
            division by zero, resulting in a band of zeros.
            """
            band_name = ee.String(band_name)
            min_key = band_name.cat('_min')
            max_key = band_name.cat('_max')
            
            min_val = ee.Number(minMax.get(min_key))
            max_val = ee.Number(minMax.get(max_key))
            range_val = max_val.subtract(min_val)
            
            # Avoid division by zero
            range_val = ee.Number(ee.Algorithms.If(
                range_val.eq(0), 1, range_val
            ))
            
            normalized = image.select([band_name]).subtract(min_val).divide(range_val)
            return normalized
        
        # Apply to all bands
        band_names = image.bandNames()
        normalized_bands = band_names.map(lambda b: normalize_band(b))
        
        return ee.ImageCollection(normalized_bands).toBands().rename(band_names)
    
    def add_training_data(self,
                         training_points: Union[str, ee.FeatureCollection] = None,
                         training_polygons: Union[str, ee.FeatureCollection] = None,
                         class_property: str = 'class',
                         points_per_class: int = 100,
                         train_fraction: float = 0.7,
                         seed: int = 0) -> None:
        """
        Add training data for supervised classification.
        
        Parameters
        ----------
        training_points : str or ee.FeatureCollection
            Point features with class labels (shapefile path or ee.FeatureCollection)
        training_polygons : str or ee.FeatureCollection
            Polygon features to sample points from
        class_property : str
            Property containing class values
        points_per_class : int
            If using polygons, number of points to sample per class
        train_fraction : float
            Share of the samples used for training; the rest is kept for
            validation. Default 0.7.
        seed : int
            Seed of the random train/validation split. The split is drawn on
            the samples before reading the feature stack, so the same points
            and seed give the same split whatever the stack — which is what
            makes classifiers built on different stacks comparable.

        Raises
        ------
        ValueError
            If no feature stack has been created or if neither `points`
            nor `polygons` are provided.
        ee.EEException
            If Earth Engine sampling fails when extracting training data.
        """
        print("Loading training data...")
        self.class_property = class_property
        
        if self.feature_stack is None:
            raise ValueError("Create feature stack first using create_feature_stack()")
        
        # Load points
        if training_points is not None:
            if isinstance(training_points, str):
                if training_points.endswith('.shp'):
                    gdf = gpd.read_file(training_points)
                    training_fc = geemap.geopandas_to_ee(gdf)
                elif training_points.endswith('.geojson'):
                    training_fc = geemap.geojson_to_ee(training_points)
            else:
                training_fc = training_points
                
        # Load polygons and sample
        elif training_polygons is not None:
            if isinstance(training_polygons, str):
                if training_polygons.endswith('.shp'):
                    gdf = gpd.read_file(training_polygons)
                    polygons_fc = geemap.geopandas_to_ee(gdf)
                elif training_polygons.endswith('.geojson'):
                    polygons_fc = geemap.geojson_to_ee(training_polygons)
            else:
                polygons_fc = training_polygons
            
            # Sample points from polygons
            training_fc = self.feature_stack.sampleRegions(
                collection=polygons_fc,
                properties=[class_property],
                scale=self.scale,
                numPixels=points_per_class,
                geometries=True
            )
        else:
            raise ValueError("Provide either training_points or training_polygons")
        
        # Draw the train/validation split on the samples themselves, before
        # reading the stack: drawn afterwards (as before 1.6.0) it depended on
        # the stack's bands, and two stacks sampled at the same points got
        # different splits
        training_fc = ee.FeatureCollection(training_fc).randomColumn('random', seed)

        # Sample feature values at training locations
        self.training_data = self.feature_stack.sampleRegions(
            collection=training_fc,
            properties=[class_property, 'random'],
            scale=self.scale
        )
        
        # Get sample count
        sample_count = self.training_data.size().getInfo()
        print(f"Training data loaded: {sample_count} samples")
        
        # Split train/validation
        training_split = self.training_data.filter(ee.Filter.lt('random', train_fraction))
        validation_split = self.training_data.filter(ee.Filter.gte('random', train_fraction))
        
        self.training_data = training_split
        self.validation_data = validation_split
        
        train_size = training_split.size().getInfo()
        val_size = validation_split.size().getInfo()
        print(f"Split: {train_size} training, {val_size} validation")

    def classify_supervised(self,
                          algorithm: str = 'random_forest',
                          train_fraction: float = 0.7,
                          params: Dict = None) -> ee.Image:
        """
        Perform supervised classification.
        
        Parameters
        ----------
        algorithm : str
            Classification algorithm:
            - 'random_forest': Random Forest (default)
            - 'svm': Support Vector Machine
            - 'cart': Classification and Regression Trees
            - 'naive_bayes': Naive Bayes
            - 'gradient_tree': Gradient Tree Boost
        train_fraction : float
            Has no effect, and never had: the split is made when the samples
            are loaded, with ``add_training_data(train_fraction=...)``. Kept
            so existing calls do not break; a value other than 0.7 prints a
            warning.
        params : dict
            Algorithm-specific parameters
            
        Returns
        -------
        ee.Image
            Classified image

        Raises
        ------
        ValueError
            If training data has not been added or if the classifier
            `algorithm` is not supported.
        ee.EEException
            If supervised classification fails in Earth Engine.
        """
        if self.training_data is None:
            raise ValueError("Add training data first using add_training_data()")

        if train_fraction != 0.7:
            print("Warning: classify_supervised(train_fraction=...) has no effect; "
                  "pass it to add_training_data() instead.")
        
        print(f"Training {algorithm} classifier...")
        
        # Get input bands
        bands = self.feature_stack.bandNames()
        
        # Select classifier
        if algorithm == 'random_forest':
            default_params = {
                'numberOfTrees': 100,
                'variablesPerSplit': None,  # sqrt(n_features)
                'minLeafPopulation': 1,
                'bagFraction': 0.5,
                'maxNodes': None
            }
            if params:
                default_params.update(params)
            
            self.classifier = ee.Classifier.smileRandomForest(**default_params)
            
        elif algorithm == 'svm':
            default_params = {
                'kernelType': 'RBF',
                'gamma': 0.5,
                'cost': 10
            }
            if params:
                default_params.update(params)
            
            self.classifier = ee.Classifier.libsvm(**default_params)
            
        elif algorithm == 'cart':
            default_params = {
                'maxNodes': None,
                'minLeafPopulation': 1
            }
            if params:
                default_params.update(params)
            
            self.classifier = ee.Classifier.smileCart(**default_params)
            
        elif algorithm == 'naive_bayes':
            self.classifier = ee.Classifier.smileNaiveBayes()
            
        elif algorithm == 'gradient_tree':
            default_params = {
                'numberOfTrees': 50,
                'shrinkage': 0.05,
                'samplingRate': 0.7,
                'maxNodes': None
            }
            if params:
                default_params.update(params)
            
            self.classifier = ee.Classifier.smileGradientTreeBoost(**default_params)
            
        else:
            raise ValueError(f"Unknown algorithm: {algorithm}")
        
        self.algorithm = algorithm

        # Train classifier
        self.classifier = self.classifier.train(
            features=self.training_data,
            classProperty=self.class_property,
            inputProperties=bands
        )
        
        # Apply classifier
        self.classified_image = self.feature_stack.classify(self.classifier)
        
        print(f"{algorithm} classification complete")
        
        # Calculate accuracy if validation data exists
        if self.validation_data is not None:
            self._calculate_accuracy()
        
        return self.classified_image
    
    def classify_unsupervised(self,
                            algorithm: str = 'kmeans',
                            n_clusters: int = 10,
                            max_iterations: int = 20,
                            params: Dict = None) -> ee.Image:
        """
        Perform unsupervised classification (clustering).
        
        Parameters
        ----------
        algorithm : str
            Clustering algorithm:
            - 'kmeans': K-means clustering (default)
            - 'cascade_kmeans': Cascade K-means
            - 'lda': Latent Dirichlet Allocation
        n_clusters : int
            Number of clusters
        max_iterations : int
            Maximum iterations
        params : dict
            Algorithm-specific parameters
            
        Returns
        -------
        ee.Image
            Clustered image

        Raises
        ------
        ValueError
            If no feature stack has been created or if `algorithm`
            is not one of {'kmeans', 'gmm'}.
        ee.EEException
            If unsupervised classification fails in Earth Engine.
        """
        if self.feature_stack is None:
            raise ValueError("Create feature stack first")
        
        print(f"Performing {algorithm} clustering with {n_clusters} clusters...")
        
        # Sample input data for clustering
        training_data = self.feature_stack.sample(
            region=self.roi,
            scale=self.scale,
            numPixels=5000,
            geometries=True
        )
        
        # Select clusterer
        if algorithm == 'kmeans':
            clusterer = ee.Clusterer.wekaKMeans(
                nClusters=n_clusters,
                maxIterations=max_iterations
            )
            
        elif algorithm == 'cascade_kmeans':
            clusterer = ee.Clusterer.wekaCascadeKMeans(
                minClusters=2,
                maxClusters=n_clusters
            )
            
        elif algorithm == 'lda':
            clusterer = ee.Clusterer.wekaLVQ(
                numClusters=n_clusters
            )
            
        else:
            raise ValueError(f"Unknown algorithm: {algorithm}")
        
        # Train clusterer
        clusterer = clusterer.train(training_data)
        
        # Apply to image
        self.classified_image = self.feature_stack.cluster(clusterer)
        
        print(f"Clustering complete: {n_clusters} clusters")
        
        return self.classified_image
    
    def _calculate_accuracy(self):
        """
        Calculate accuracy metrics from validation samples.

        Computes confusion matrix, overall accuracy, kappa coefficient,
        producer's and user's accuracy per class.

        Parameters
        ----------
        validation : ee.FeatureCollection
            FeatureCollection containing validation samples with
            reference and predicted labels.

        Returns
        -------
        dict
            Accuracy metrics, including:
            - ``'confusion_matrix'`` : numpy.ndarray
            - ``'overall_accuracy'`` : float
            - ``'kappa'`` : float
            - ``'producer_accuracy'`` : dict
            - ``'user_accuracy'`` : dict
        """
        if self.validation_data is None or self.classifier is None:
            return
        
        # Classify validation data
        validated = self.validation_data.classify(self.classifier)
        
        # Create confusion matrix
        confusion_matrix = validated.errorMatrix(self.class_property, 'classification')
        
        # Calculate metrics
        self.accuracy_results = {
            'overall_accuracy': confusion_matrix.accuracy().getInfo(),
            'kappa': confusion_matrix.kappa().getInfo(),
            'producers_accuracy': confusion_matrix.producersAccuracy().getInfo(),
            'consumers_accuracy': confusion_matrix.consumersAccuracy().getInfo(),
            'confusion_matrix': confusion_matrix.array().getInfo()
        }
        
        print(f"Overall Accuracy: {self.accuracy_results['overall_accuracy']:.3f}")
        print(f"Kappa: {self.accuracy_results['kappa']:.3f}")

    def export_results(self, description: str, scale: Optional[float] = None, region: Optional[ee.Geometry] = None):
        """
        Export the classified image to Google Drive or Earth Engine Asset.

        Parameters
        ----------
        description : str
            Name of the export task.
        scale : float, optional
            Spatial resolution in meters. Defaults to the pixel size of the
            feature stack (:attr:`scale`), so the map is exported at the
            resolution it was classified at. Before 1.6.0 it was always 30.
        region : ee.Geometry, optional
            Geometry defining the export area. If None, uses the full image extent.

        Returns
        -------
        ee.batch.Task
            The Earth Engine export task object.

        Raises
        ------
        ValueError
            If no classified image is available.
        ee.EEException
            If the export task could not be created.
        """
        if self.classified_image is None:
            raise ValueError("No classified image to export. Run classification first.")

        if scale is None:
            scale = self.scale

        export_args = dict(
            image=self.classified_image,
            description=description,
            scale=scale,
            region=region
        )
        # Keep the grid the stack was built on, so exported pixels match
        # the classified ones
        if self.crs is not None:
            export_args['crs'] = self.crs

        task = ee.batch.Export.image.toDrive(**export_args)
        task.start()
        return task
    
    def plot_confusion_matrix(self, labels: List[str]):
        """
        Plot the confusion matrix of the classification results.

        Parameters
        ----------
        labels : list of str
            List of class names in the same order as the matrix.

        Returns
        -------
        matplotlib.axes.Axes
            Axis object containing the confusion matrix plot.

        Raises
        ------
        ValueError
            If no confusion matrix is available (classification or accuracy not run).
        """
        cm = self.accuracy_results.get('confusion_matrix')
        if cm is None:
            raise ValueError("Confusion matrix not available. Run classification and accuracy first.")

        fig, ax = plt.subplots(figsize=(6, 5))
        sns.heatmap(cm, annot=True, fmt="d", cmap="Blues",
                    xticklabels=labels, yticklabels=labels, ax=ax)
        ax.set_xlabel("Predicted")
        ax.set_ylabel("Reference")
        ax.set_title("Confusion Matrix")
        return ax
    
    def get_accuracy_report(self) -> pd.DataFrame:
        """
        Return accuracy metrics as a pandas DataFrame.

        Returns
        -------
        pandas.DataFrame
            Table with overall accuracy, kappa, producer's and user's
            accuracy for each class.

        Raises
        ------
        ValueError
            If no accuracy metrics are available.
        """
        if not self.accuracy_results:
            raise ValueError("No accuracy metrics available. Run classification first.")

        # Earth Engine returns producer's accuracy as a column (N x 1) and
        # user's (consumer's) accuracy as a row (1 x N), both indexed by class
        # value; classes missing from the validation set are left out
        producers = [row[0] for row in self.accuracy_results['producers_accuracy']]
        consumers = self.accuracy_results['consumers_accuracy'][0]
        matrix = np.array(self.accuracy_results['confusion_matrix'])

        rows = []
        for cls, (pa, ua) in enumerate(zip(producers, consumers)):
            if matrix[cls, :].sum() == 0 and matrix[:, cls].sum() == 0:
                continue
            rows.append({
                'Class': cls,
                'ProducerAccuracy': pa,
                'UserAccuracy': ua
            })
        df = pd.DataFrame(rows)
        df.loc[len(df)] = {
            'Class': 'Overall',
            'ProducerAccuracy': self.accuracy_results.get('overall_accuracy'),
            'UserAccuracy': self.accuracy_results.get('kappa')
        }
        return df
    
    def get_feature_importance(self) -> Dict[str, float]:
        """
        Get feature importance scores of the trained classifier.

        Available for the tree-based algorithms of :meth:`classify_supervised`
        (``'random_forest'``, ``'cart'`` and ``'gradient_tree'``). Scores are
        relative within a model: compare features, not models.

        Returns
        -------
        dict
            Mapping of feature (band) names to importance scores, sorted from
            most to least important.

        Raises
        ------
        ValueError
            If no tree-based classifier has been trained.
        """
        # The Python type of every trained classifier is ee.Classifier, so the
        # algorithm has to come from what classify_supervised() was asked for
        if self.classifier is None or self.algorithm not in ('random_forest', 'cart', 'gradient_tree'):
            raise ValueError(
                "Feature importance is only available after classify_supervised() "
                "with algorithm='random_forest', 'cart' or 'gradient_tree'."
            )

        importance = self.classifier.explain().get('importance').getInfo()
        return dict(sorted(importance.items(), key=lambda kv: kv[1], reverse=True))
