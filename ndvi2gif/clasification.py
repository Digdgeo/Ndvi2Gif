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
        Instance providing temporal NDVI composites and configuration.
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
        Last year (exclusive) of the analysis.
    sat : str
        Satellite name used (e.g., 'S2', 'L8').
    """
    
    def __init__(self, ndvi_seasonality_instance):
        """
        Initialize the classifier with an ``NdviSeasonality`` instance.

        Parameters
        ----------
        ndvi_seasonality_instance : NdviSeasonality
            An initialized instance of the NdviSeasonality processor, providing
            temporal composites, ROI, and analysis configuration.

        Notes
        -----
        The constructor inherits spatial and temporal parameters directly from
        the provided processor:
        - ROI (region of interest)
        - Number of periods per year
        - Start and end years
        - Satellite identifier
        """
        self.processor = ndvi_seasonality_instance
        self.feature_stack = None
        self.training_data = None
        self.validation_data = None
        self.classifier = None
        self.classified_image = None
        self.accuracy_results = None
        
        # Inherit parameters
        self.roi = self.processor.roi
        self.periods = self.processor.periods
        self.start_year = self.processor.start_year
        self.end_year = self.processor.end_year
        self.sat = self.processor.sat
        
        print(f"LandCoverClassifier initialized for {self.sat}")
        print(f"Period: {self.start_year}-{self.end_year-1}, {self.periods} periods/year")
    
    def create_feature_stack(self,
                           indices: List[str] = None,
                           include_statistics: bool = True,
                           normalize: bool = True) -> ee.Image:
        """
        Create multi-temporal feature stack for classification.
        
        Parameters
        ----------
        indices : list of str, optional
            Indices to stack. If None, uses current index.
        include_statistics : bool
            Add temporal statistics (mean, std, max, min)
        normalize : bool
            Normalize to [0,1] range
            
        Returns
        -------
        ee.Image
            Multi-band feature stack

        Raises
        ------
        ValueError
            If `indices` contains unsupported names or is empty.
        ee.EEException
            If Earth Engine image processing fails when computing the stack.
        """
        print("Creating feature stack...")
        
        # Use current index if none specified
        if indices is None:
            indices = [self.processor.index]
        
        # Validate indices
        available = self.processor.sensor_indices[self.sat]
        invalid = [idx for idx in indices if idx not in available]
        if invalid:
            raise ValueError(f"Invalid indices for {self.sat}: {invalid}")
        
        feature_bands = []
        
        # Process each index
        for idx_name in indices:
            print(f"  Processing {idx_name}...")
            
            # Temporarily set processor index
            original_index = self.processor.index
            self.processor.index = idx_name
            
            # Generate composites
            collection = self.processor.get_year_composite()
            
            # Convert to list for iteration
            img_list = collection.toList(100)
            
            # Stack all years
            for year_idx in range(self.end_year - self.start_year):
                year = self.start_year + year_idx
                year_image = ee.Image(img_list.get(year_idx))
                
                # Rename bands with descriptive names
                for period_name in self.processor.period_names:
                    band = year_image.select(period_name)
                    band_name = f"{idx_name}_{year}_{period_name}"
                    feature_bands.append(band.rename(band_name))
            
            # Restore original index
            self.processor.index = original_index
        
        # Create feature stack
        self.feature_stack = ee.Image.cat(feature_bands)
        
        # Add statistics if requested
        if include_statistics:
            print("  Adding temporal statistics...")
            stats = self._compute_statistics(indices)
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
    
    def _compute_statistics(self, indices: List[str]) -> ee.Image:
        """
        Compute basic statistics for a given image.

        Calculates per-band mean, standard deviation, minimum and maximum
        values to support normalization and feature scaling.

        Parameters
        ----------
        image : ee.Image
            Earth Engine image for which statistics will be computed.

        Returns
        -------
        dict
            Dictionary containing per-band statistics:
            ``{'mean': ..., 'std': ..., 'min': ..., 'max': ...}``.
        """
        stats_bands = []
        
        for idx_name in indices:
            # Select all bands for this index
            idx_pattern = f"{idx_name}_.*"
            idx_bands = self.feature_stack.select(idx_pattern)
            
            # Calculate statistics
            mean = idx_bands.reduce(ee.Reducer.mean()).rename(f"{idx_name}_mean")
            std = idx_bands.reduce(ee.Reducer.stdDev()).rename(f"{idx_name}_std")
            max_val = idx_bands.reduce(ee.Reducer.max()).rename(f"{idx_name}_max")
            min_val = idx_bands.reduce(ee.Reducer.min()).rename(f"{idx_name}_min")
            
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
            scale=30,
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
                         points_per_class: int = 100) -> None:
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

        Raises
        ------
        ValueError
            If no feature stack has been created or if neither `points`
            nor `polygons` are provided.
        ee.EEException
            If Earth Engine sampling fails when extracting training data.
        """
        print("Loading training data...")
        
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
                scale=10,
                numPixels=points_per_class,
                geometries=True
            )
        else:
            raise ValueError("Provide either training_points or training_polygons")
        
        # Sample feature values at training locations
        self.training_data = self.feature_stack.sampleRegions(
            collection=training_fc,
            properties=[class_property],
            scale=10
        )
        
        # Get sample count
        sample_count = self.training_data.size().getInfo()
        print(f"Training data loaded: {sample_count} samples")
        
        # Split train/validation (70/30)
        self.training_data = self.training_data.randomColumn('random')
        training_split = self.training_data.filter(ee.Filter.lt('random', 0.7))
        validation_split = self.training_data.filter(ee.Filter.gte('random', 0.7))
        
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
            Fraction of data for training (rest for validation)
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
        
        # Train classifier
        self.classifier = self.classifier.train(
            features=self.training_data,
            classProperty='class',
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
            scale=30,
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
        confusion_matrix = validated.errorMatrix('class', 'classification')
        
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

    def export_results(self, description: str, scale: int = 30, region: Optional[ee.Geometry] = None):
        """
        Export the classified image to Google Drive or Earth Engine Asset.

        Parameters
        ----------
        description : str
            Name of the export task.
        scale : int, optional
            Spatial resolution in meters. Default is 30.
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

        task = ee.batch.Export.image.toDrive(
            image=self.classified_image,
            description=description,
            scale=scale,
            region=region
        )
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

        rows = []
        for cls, pa in self.accuracy_results['producer_accuracy'].items():
            ua = self.accuracy_results['user_accuracy'].get(cls, None)
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
        Get feature importance scores from a Random Forest classifier.

        Returns
        -------
        dict
            Mapping of feature names to importance scores.

        Raises
        ------
        ValueError
            If classifier is not a Random Forest or not trained.
        """
        if not hasattr(self, 'classifier') or 'RandomForest' not in str(type(self.classifier)):
            raise ValueError("Feature importance is only available for Random Forest classifiers.")

        importance = self.classifier.explain()['importance']
        return importance
