"""
Time Series Analysis Module for ndvi2gif v0.5.0 - CLEANED VERSION

This module provides comprehensive temporal analysis capabilities for remote sensing data.
It's designed to work seamlessly with the NdviSeasonality class.

Author: Diego García Díaz
Date: 2024
License: MIT
"""

import os
import ee
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import seaborn as sns
from scipy import stats, signal, interpolate
from scipy.optimize import curve_fit
from datetime import datetime, timedelta
from typing import Optional, List, Dict, Union, Tuple, Any
from matplotlib.ticker import MaxNLocator
import warnings

# Set style for better plots
plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette("husl")

# ========== CONFIGURACIÓN GLOBAL UNIFICADA ==========
plt.rcParams.update({
    'font.size': 11,              # Tamaño base mejorado
    'axes.titlesize': 13,         # Títulos más grandes
    'axes.labelsize': 11,         # Etiquetas de ejes
    'xtick.labelsize': 10,        # Ticks X
    'ytick.labelsize': 10,        # Ticks Y
    'legend.fontsize': 10,        # Leyendas
    'legend.title_fontsize': 11,
    'figure.titlesize': 15,       # Título principal
    'lines.linewidth': 2.5,       # Líneas más gruesas
    'lines.markersize': 7,        # Marcadores más grandes
    'axes.grid': True,
    'grid.alpha': 0.3,
    'axes.spines.top': False,     # Quitar borde superior
    'axes.spines.right': False,   # Quitar borde derecho
})

class TimeSeriesAnalyzer:
    """
    Advanced time series analysis for Earth Engine remote sensing data.
    
    This class provides methods for extracting, analyzing, and visualizing
    temporal patterns in satellite data processed by NdviSeasonality.
    
    Version: 0.5.0 - Improved visualization without overlapping elements
    """
    
    def __init__(self, ndvi_seasonality_instance):
        """Initialize TimeSeriesAnalyzer with an NdviSeasonality instance."""
        self.processor = ndvi_seasonality_instance
        self.time_series_cache = {}
        
        # Inherit key attributes
        self.roi = self.processor.roi
        self.periods = self.processor.periods
        self.start_year = self.processor.start_year
        self.end_year = self.processor.end_year
        self.sat = self.processor.sat
        self.index = self.processor.index
        self.key = self.processor.key
        self.period_names = self.processor.period_names
        self.period_dates = self.processor.period_dates
    
    # ========== CORE METHODS ==========
    
    def extract_time_series(self, 
                       point: Optional[Union[ee.Geometry.Point, Tuple[float, float]]] = None,
                       reducer: str = 'mean',
                       scale: int = 30,
                       use_cache: bool = True) -> pd.DataFrame:
        """
        Extract complete time series for a point or region.
        
        This method combines all temporal periods across all years into a single
        continuous time series suitable for analysis.
        
        Parameters
        ----------
        point : ee.Geometry.Point, tuple, or None
            Location for extraction:
            - None: uses ROI centroid
            - tuple: (longitude, latitude)
            - ee.Geometry.Point: direct point geometry
            - ee.Geometry.Polygon: for spatial averaging
        reducer : str
            Spatial reducer if using polygon: 'mean', 'median', 'max', 'min', 'stdDev'
        scale : int
            Scale in meters for spatial reduction
        use_cache : bool
            Whether to use cached results if available
            
        Returns
        -------
        pd.DataFrame
            DataFrame with columns:
            - date: datetime index
            - value: index values  
            - year: year
            - period: period name
            - doy: day of year
            - season: meteorological season
            - month: month number
        """
        # Create cache key
        cache_key = f"{point}_{reducer}_{scale}"
        
        if use_cache and cache_key in self.time_series_cache:
            print("Using cached time series data")
            return self.time_series_cache[cache_key]
        
        # Handle different point inputs - SOLUCIÓN INTEGRADA
        extraction_geometry = None
        geometry_type = 'point'  # Default assumption
        
        if point is None:
            # Use ROI centroid
            extraction_geometry = self.roi.centroid()
            coords_info = extraction_geometry.coordinates().getInfo()
            print(f"Using ROI centroid: {coords_info}")
            geometry_type = 'point'
            
        elif isinstance(point, tuple) and len(point) == 2:
            # Handle tuple coordinates (lon, lat)
            try:
                lon, lat = float(point[0]), float(point[1])
                extraction_geometry = ee.Geometry.Point([lon, lat])
                print(f"Using point coordinates: ({lon}, {lat})")
                geometry_type = 'point'
            except (ValueError, TypeError) as e:
                raise ValueError(f"Invalid coordinate tuple: {point}. Must be (lon, lat) as numbers.")
                
        else:
            # Assume it's an Earth Engine geometry object
            # Avoid isinstance checks that can fail with EE objects
            extraction_geometry = point
            try:
                # Try to get coordinates to determine if it's a point
                coords = extraction_geometry.coordinates().getInfo()
                if isinstance(coords, list) and len(coords) == 2:
                    geometry_type = 'point'
                    print(f"Using EE Point: {coords}")
                else:
                    geometry_type = 'polygon'
                    print(f"Using EE Polygon with {reducer} reducer")
            except:
                # If we can't determine, assume polygon
                geometry_type = 'polygon'
                print(f"Using EE geometry with {reducer} reducer")
        
        # Initialize containers
        time_series_data = []
        
        # Progress tracking
        total_periods = (self.end_year - self.start_year) * self.periods
        current = 0
        
        print(f"Extracting {total_periods} temporal periods...")
        
        # Get reducer function
        try:
            reducer_func = getattr(ee.Reducer, reducer)()
        except AttributeError:
            print(f"Warning: Unknown reducer '{reducer}', using 'mean'")
            reducer_func = ee.Reducer.mean()
        
        # Process each year and period
        for year in range(self.start_year, self.end_year):
            for period_idx in range(self.periods):
                current += 1
                
                # Progress indicator
                if current % 10 == 0 or current == total_periods:
                    print(f"Progress: {current}/{total_periods} ({current*100/total_periods:.1f}%)", end='\r')
                
                try:
                    # Get period composite
                    period_composite = self.processor.get_period_composite(year, period_idx)
                    
                    # Extract value based on geometry type
                    if geometry_type == 'point':
                        # Point extraction with small buffer for stability
                        buffer_geom = extraction_geometry.buffer(scale/2)
                        extracted = period_composite.reduceRegion(
                            reducer=reducer_func,
                            geometry=buffer_geom,
                            scale=scale,
                            maxPixels=1e9,
                            bestEffort=True
                        )
                    else:
                        # Polygon extraction
                        extracted = period_composite.reduceRegion(
                            reducer=reducer_func,
                            geometry=extraction_geometry,
                            scale=scale,
                            maxPixels=1e9,
                            bestEffort=True
                        )
                    
                    # Get value - IMPROVED ERROR HANDLING
                    value_dict = extracted.getInfo()
                    
                    if value_dict and len(value_dict) > 0:
                        # Extract first value (handles different band names)
                        value = list(value_dict.values())[0]
                        
                        # Enhanced validation
                        if value is not None and not (isinstance(value, float) and np.isnan(value)):
                            # Convert to float for consistency
                            try:
                                value = float(value)
                            except (ValueError, TypeError):
                                print(f"\nWarning: Could not convert value to float: {value}")
                                continue
                                
                            # Calculate temporal attributes
                            period_start, period_end = self.period_dates[period_idx]
                            start_date = datetime.strptime(f"{year}{period_start}", "%Y-%m-%d")
                            end_date = datetime.strptime(f"{year}{period_end}", "%Y-%m-%d")
                            mid_date = start_date + (end_date - start_date) / 2
                            
                            # Determine season
                            month = mid_date.month
                            if month in [12, 1, 2]:
                                season = 'winter'
                            elif month in [3, 4, 5]:
                                season = 'spring'
                            elif month in [6, 7, 8]:
                                season = 'summer'
                            else:
                                season = 'autumn'
                            
                            # Append to results
                            time_series_data.append({
                                'date': mid_date,
                                'value': value,
                                'year': year,
                                'period': self.period_names[period_idx],
                                'doy': mid_date.timetuple().tm_yday,
                                'season': season,
                                'month': month
                            })
                        else:
                            # Value is None or NaN
                            continue
                    else:
                        # Empty result
                        continue
                        
                except Exception as e:
                    print(f"\nWarning: Failed to process {year} {self.period_names[period_idx]}: {str(e)}")
                    continue
        
        print(f"\nSuccessfully extracted {len(time_series_data)} data points")
        
        # Create DataFrame
        df = pd.DataFrame(time_series_data)
        
        if len(df) > 0:
            df = df.sort_values('date').reset_index(drop=True)
            df['date'] = pd.to_datetime(df['date'])
            
            # Cache results
            if use_cache:
                self.time_series_cache[cache_key] = df
        else:
            print("Warning: No valid data points extracted")
        
        return df
    
    def analyze_trend(self, 
                     df: Optional[pd.DataFrame] = None,
                     method: str = 'mann_kendall',
                     alpha: float = 0.05) -> Dict[str, Any]:
        """
        Perform comprehensive trend analysis on time series.
        
        Parameters
        ----------
        df : pd.DataFrame or None
            Time series dataframe. If None, extracts from ROI centroid
        method : str
            Trend test method: 'mann_kendall', 'linear', 'sen_slope', 'all'
        alpha : float
            Significance level for statistical tests
            
        Returns
        -------
        dict
            Dictionary containing trend statistics
        """
        if df is None:
            df = self.extract_time_series()
        
        if len(df) < 3:
            return {'error': 'Insufficient data for trend analysis'}
        
        results = {}
        values = df['value'].values
        n = len(values)
        
        # Create time index
        time_index = np.arange(n)
        
        # Mann-Kendall Test
        if method in ['mann_kendall', 'all']:
            mk_result = self._mann_kendall_test(values)
            results['mann_kendall'] = mk_result
        
        # Linear Regression
        if method in ['linear', 'all']:
            slope, intercept, r_value, p_value, std_err = stats.linregress(time_index, values)
            
            # Calculate confidence interval
            t_stat = stats.t.ppf(1 - alpha/2, n - 2)
            ci = t_stat * std_err
            
            results['linear'] = {
                'slope': slope,
                'intercept': intercept,
                'r_squared': r_value**2,
                'p_value': p_value,
                'std_error': std_err,
                'confidence_interval': (slope - ci, slope + ci),
                'yearly_change': slope * self.periods  # Change per year
            }
        
        # Sen's Slope (Theil-Sen)
        if method in ['sen_slope', 'all']:
            sen_result = self._sen_slope(values)
            results['sen_slope'] = sen_result
        
        # Overall interpretation
        if 'mann_kendall' in results:
            mk = results['mann_kendall']
            if mk['p_value'] < alpha:
                if mk['trend'] > 0:
                    interpretation = f"Significant increasing trend (p={mk['p_value']:.3f})"
                else:
                    interpretation = f"Significant decreasing trend (p={mk['p_value']:.3f})"
            else:
                interpretation = f"No significant trend (p={mk['p_value']:.3f})"
            results['interpretation'] = interpretation
        
        return results
    
    # ========== MAIN VISUALIZATION METHODS ==========
    
    def plot_comprehensive_analysis(self,
                                   point=None,
                                   figsize=(22, 14),
                                   save_path=None):
        """
        Create comprehensive time series analysis dashboard with improved layout.
        """
        df = self.extract_time_series(point)
        
        if len(df) == 0:
            print("No data available for visualization")
            return None
        
        plt.close('all')
        
        # Create figure with optimized layout
        fig = plt.figure(figsize=figsize, constrained_layout=False)
        
        gs = fig.add_gridspec(3, 4,
                            height_ratios=[2.0, 1.4, 1.4],
                            width_ratios=[1.2, 1, 1, 1],
                            hspace=0.40,  # More vertical space
                            wspace=0.35,  # More horizontal space
                            top=0.90,
                            bottom=0.08,
                            left=0.06,
                            right=0.96)
        
        # Create all subplots
        ax1 = fig.add_subplot(gs[0, :])
        self._plot_time_series_with_trend(df, ax1)
        
        ax2 = fig.add_subplot(gs[1, 0])
        self._plot_seasonal_pattern(df, ax2)
        
        ax3 = fig.add_subplot(gs[1, 1])
        self._plot_annual_comparison(df, ax3)
        
        ax4 = fig.add_subplot(gs[1, 2])
        self._plot_trend_summary(df, ax4)
        
        ax5 = fig.add_subplot(gs[1, 3])
        self._plot_autocorrelation(df, ax5)
        
        ax6 = fig.add_subplot(gs[2, 0])
        self._plot_distribution(df, ax6)
        
        ax7 = fig.add_subplot(gs[2, 1])
        self._plot_phenology_summary(df, ax7)
        
        ax8 = fig.add_subplot(gs[2, 2])
        self._plot_data_quality(df, ax8)
        
        ax9 = fig.add_subplot(gs[2, 3])
        self._plot_seasonal_statistics(df, ax9)
        
        # Main title
        title = f"{self.sat} {self.index.upper()} - Time Series Analysis"
        subtitle = f"Period: {self.start_year}-{self.end_year-1} | " \
                  f"Resolution: {self.periods} periods/year | " \
                  f"Total: {len(df)} observations"
        
        fig.suptitle(title, fontsize=16, fontweight='bold', y=0.96)
        fig.text(0.5, 0.92, subtitle, ha='center', fontsize=12, style='italic')
        
        plt.tight_layout(rect=[0, 0.08, 1, 0.90])
        
        if save_path:
            fig.savefig(save_path, dpi=300, bbox_inches='tight',
                       facecolor='white', edgecolor='none')
            print(f"Figure saved to {save_path}")
        
        plt.show()
        return fig
    
    def plot_phenology_analysis(self, 
                               point=None,
                               method='threshold',
                               threshold_percentile=50,
                               figsize=(24, 16),
                               save_path=None):
        """
        Create comprehensive phenological analysis dashboard.
        """
        df = self.extract_time_series(point)
        
        if len(df) == 0:
            print("No data available for phenological analysis")
            return None
        
        # Extract phenology metrics
        if method == 'derivative':
            phenology_results = self.extract_phenology_metrics(df, method=method)
        else:
            phenology_results = self.extract_phenology_metrics(
                df, method=method, threshold_percentile=threshold_percentile
            )
        
        if not phenology_results:
            print("No phenological metrics could be extracted")
            return None
        
        plt.close('all')
        
        fig = plt.figure(figsize=figsize, constrained_layout=False)
        
        gs = fig.add_gridspec(3, 6,
                            height_ratios=[1.8, 1.3, 1.3],
                            width_ratios=[1.5, 1, 1, 1, 1, 1.3],
                            hspace=0.55,
                            wspace=0.25,
                            top=0.85,
                            bottom=0.15,
                            left=0.05,
                            right=0.97)
        
        # Row 1: Main time series and info
        ax1 = fig.add_subplot(gs[0, :4])
        self._plot_time_series_with_phenology(df, phenology_results, ax1, method)
        
        ax_info = fig.add_subplot(gs[0, 4:])
        self._plot_phenology_info(ax_info, method, threshold_percentile, len(phenology_results))
        
        # Row 2: 6 metric panels
        ax2 = fig.add_subplot(gs[1, 0])
        self._plot_phenology_timing(phenology_results, ax2)
        
        ax3 = fig.add_subplot(gs[1, 1])
        self._plot_phenology_amplitude(phenology_results, ax3)
        
        ax4 = fig.add_subplot(gs[1, 2])
        self._plot_phenology_rates(phenology_results, ax4)
        
        ax5 = fig.add_subplot(gs[1, 3])
        self._plot_phenology_duration(phenology_results, ax5)
        
        ax6 = fig.add_subplot(gs[1, 4])
        self._plot_phenology_statistics(phenology_results, ax6)
        
        ax7 = fig.add_subplot(gs[1, 5])
        self._plot_phenology_data_availability(phenology_results, ax7)
        
        # Row 3: Comparison and summaries
        ax8 = fig.add_subplot(gs[2, :3])
        self._plot_annual_phenology_comparison(df, phenology_results, ax8)
        
        ax9 = fig.add_subplot(gs[2, 3:5])
        self._plot_phenology_metrics_summary(phenology_results, ax9)
        
        ax10 = fig.add_subplot(gs[2, 5])
        self._plot_phenology_quality_summary(phenology_results, ax10)
        
        # Titles
        title = f"Phenological Analysis - {self.sat} {self.index.upper()}"
        if method == 'threshold':
            subtitle = f"Method: Threshold | Threshold: {threshold_percentile}% | " \
                      f"{self.start_year}-{self.end_year-1} | {len(phenology_results)} years"
        else:
            subtitle = f"Method: {method.capitalize()} | " \
                      f"{self.start_year}-{self.end_year-1} | {len(phenology_results)} years"
        
        fig.suptitle(title, fontsize=15, fontweight='bold', y=0.92)
        fig.text(0.5, 0.88, subtitle, ha='center', fontsize=11, style='italic')
        
        if save_path:
            fig.savefig(save_path, dpi=300, bbox_inches='tight',
                       facecolor='white', edgecolor='none')
            print(f"Phenological analysis saved to {save_path}")
        
        plt.show()
        return fig
    
    # ========== PHENOLOGY EXTRACTION METHODS ==========
    
    def extract_phenology_metrics(self, 
                             df: Optional[pd.DataFrame] = None,
                             method: str = 'threshold',
                             threshold_percentile: float = 50,
                             smoothing: bool = True,
                             min_season_length: int = 60) -> Dict[str, Any]:
        """
        Extract comprehensive phenological metrics from time series data.
        
        Calculates key phenological parameters including Start of Season (SOS),
        End of Season (EOS), Peak of Season (POS), Length of Season (LOS),
        and various amplitude and timing metrics for each year.
        
        Parameters
        ----------
        df : pd.DataFrame or None
            Time series dataframe. If None, extracts from ROI centroid
        method : str
            Phenology extraction method:
            - 'threshold': Uses percentile-based thresholds
            - 'derivative': Uses first/second derivatives
            - 'logistic': Fits logistic curves to data
        threshold_percentile : float
            Percentile for threshold method (0-100). Default 50 (median).
        smoothing : bool
            Apply Savitzky-Golay smoothing before analysis. Default True.
        min_season_length : int
            Minimum season length in days. Default 60.
            
        Returns
        -------
        dict
            Dictionary with phenological metrics for each year:
            {
                year: {
                    'sos': Start of Season (day of year),
                    'eos': End of Season (day of year), 
                    'pos': Peak of Season (day of year),
                    'los': Length of Season (days),
                    'amplitude': Seasonal amplitude (max - baseline),
                    'baseline': Seasonal baseline (minimum value),
                    'peak_value': Peak vegetation value,
                    'auc': Area Under Curve (seasonal integral),
                    'sos_value': Vegetation value at SOS,
                    'eos_value': Vegetation value at EOS,
                    'growth_rate': Rate of vegetation increase,
                    'senescence_rate': Rate of vegetation decrease
                }
            }
        """
        if df is None:
            df = self.extract_time_series()
        
        if len(df) == 0:
            return {'error': 'No data available for phenology analysis'}
        
        phenology_results = {}
        
        # Process each year separately
        for year in df['year'].unique():
            year_data = df[df['year'] == year].copy()
            
            if len(year_data) < 4:  # Need minimum data points
                print(f"Insufficient data for year {year} ({len(year_data)} points)")
                continue
                
            # Sort by day of year for proper temporal order
            year_data = year_data.sort_values('doy').reset_index(drop=True)
            
            try:
                if method == 'threshold':
                    metrics = self._phenology_threshold_method(
                        year_data, threshold_percentile, smoothing, min_season_length
                    )
                elif method == 'derivative':
                    metrics = self._phenology_derivative_method(
                        year_data, smoothing, min_season_length
                    )
                elif method == 'logistic':
                    metrics = self._phenology_logistic_method(
                        year_data, smoothing, min_season_length
                    )
                else:
                    raise ValueError(f"Unknown method: {method}")
                    
                if metrics is not None:
                    phenology_results[year] = metrics
                    
            except Exception as e:
                print(f"Error processing phenology for year {year}: {e}")
                continue
        
        return phenology_results
    
    # ========== HELPER METHODS ==========
    
    def _mann_kendall_test(self, x: np.ndarray) -> Dict[str, Any]:
        """Mann-Kendall trend test implementation."""
        n = len(x)
        s = 0
        
        # Calculate S statistic
        for i in range(n-1):
            for j in range(i+1, n):
                s += np.sign(x[j] - x[i])
        
        # Calculate variance
        var_s = n * (n - 1) * (2 * n + 5) / 18
        
        # Calculate Z-score
        if s > 0:
            z = (s - 1) / np.sqrt(var_s)
        elif s < 0:
            z = (s + 1) / np.sqrt(var_s)
        else:
            z = 0
        
        # Calculate p-value (two-tailed)
        p_value = 2 * (1 - stats.norm.cdf(abs(z)))
        
        # Calculate Kendall's tau
        tau = s / (n * (n - 1) / 2)
        
        return {
            'trend': np.sign(s),
            'tau': tau,
            'z_score': z,
            'p_value': p_value,
            's_statistic': s,
            'var_s': var_s
        }
    
    def _sen_slope(self, x: np.ndarray) -> Dict[str, Any]:
        """Calculate Sen's slope estimator."""
        n = len(x)
        slopes = []
        
        # Calculate all pairwise slopes
        for i in range(n-1):
            for j in range(i+1, n):
                slopes.append((x[j] - x[i]) / (j - i))
        
        # Sen's slope is the median of all slopes
        sen_slope = np.median(slopes)
        
        # Calculate confidence interval using Kendall's method
        # Simplified version
        ci_lower = np.percentile(slopes, 2.5)
        ci_upper = np.percentile(slopes, 97.5)
        
        return {
            'slope': sen_slope,
            'confidence_interval': (ci_lower, ci_upper),
            'n_slopes': len(slopes)
        }
    
    # ========== PLOT METHODS - COMPREHENSIVE ANALYSIS ==========
    
    def _plot_time_series_with_trend(self, df, ax):
        """Serie temporal con trend mejorada"""
        # Datos observados
        ax.plot(df['date'], df['value'],
                'o-', alpha=0.6, label='Observed',
                markersize=4, color='steelblue',
                linewidth=1.5, markerfacecolor='lightblue',
                markeredgecolor='steelblue', markeredgewidth=0.5)
        
        # Línea de tendencia
        x_numeric = np.arange(len(df))
        z = np.polyfit(x_numeric, df['value'], 1)
        p = np.poly1d(z)
        
        ax.plot(df['date'], p(x_numeric),
                '--', color='red', linewidth=2.5, alpha=0.8,
                label=f'Trend ({z[0]:.2e}/period)')
        
        # Intervalo de confianza
        residuals = df['value'] - p(x_numeric)
        std_resid = np.std(residuals)
        ax.fill_between(df['date'],
                        p(x_numeric) - 1.96*std_resid,
                        p(x_numeric) + 1.96*std_resid,
                        alpha=0.15, color='red', label='95% CI')
        
        # Box de estadísticas MEJORADO (más compacto)
        mk = self._mann_kendall_test(df['value'].values)
        
        if mk['p_value'] < 0.05:
            trend_symbol = "↗" if mk['trend'] > 0 else "↘"
            trend_text = "Significant"
        else:
            trend_symbol = "→"
            trend_text = "No trend"
        
        stats_text = f"Mann-Kendall\nτ={mk['tau']:.3f}\np={mk['p_value']:.3f}\n{trend_symbol} {trend_text}"
        
        # Posicionar mejor el box
        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes,
                verticalalignment='top', fontsize=10,
                bbox=dict(boxstyle='round,pad=0.5', facecolor='wheat',
                         alpha=0.8, edgecolor='orange', linewidth=1))
        
        # Configuración
        ax.set_xlabel('Date', fontsize=12, fontweight='bold')
        ax.set_ylabel(f'{self.index.upper()} Value', fontsize=12, fontweight='bold')
        ax.set_title('Time Series with Trend Analysis',
                    fontsize=14, fontweight='bold', pad=15)
        
        # Leyenda optimizada
        ax.legend(loc='upper right', frameon=True, fancybox=True,
                 shadow=True, fontsize=10, framealpha=0.9)
        
        ax.grid(True, alpha=0.3)
        ax.set_axisbelow(True)
        
        # Formato de fechas
        import matplotlib.dates as mdates
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
        ax.xaxis.set_major_locator(mdates.MonthLocator(interval=4))
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')
    
    def _plot_seasonal_pattern(self, df, ax):
        """Patrón estacional mejorado"""
        if self.periods <= 12:
            period_data = []
            period_labels = []
            
            for p in self.period_names:
                if p in df['period'].unique():
                    data = df[df['period'] == p]['value'].values
                    if len(data) > 0:
                        period_data.append(data)
                        # Abreviar nombres largos
                        label = p[:3] if len(p) > 5 else p
                        period_labels.append(label)
            
            if period_data:
                bp = ax.boxplot(period_data, labels=period_labels, 
                               patch_artist=True, widths=0.7)
                
                # Colores mejorados
                colors = plt.cm.Set3(np.linspace(0, 1, len(bp['boxes'])))
                for box, color in zip(bp['boxes'], colors):
                    box.set_facecolor(color)
                    box.set_alpha(0.7)
                    box.set_linewidth(1.5)
                
                # Estilo de whiskers y outliers
                for whisker in bp['whiskers']:
                    whisker.set_color('black')
                    whisker.set_linewidth(1.2)
                    whisker.set_linestyle('-')
                
                for cap in bp['caps']:
                    cap.set_color('black')
                    cap.set_linewidth(1.2)
                
                for median in bp['medians']:
                    median.set_color('red')
                    median.set_linewidth(2)
                
                for outlier in bp['fliers']:
                    outlier.set_marker('o')
                    outlier.set_markersize(4)
                    outlier.set_alpha(0.5)
        
        ax.set_xlabel('Period', fontsize=11, fontweight='bold')
        ax.set_ylabel('Value', fontsize=11, fontweight='bold')
        ax.set_title('Seasonal Pattern', fontsize=12, fontweight='bold', pad=10)
        
        if self.periods > 6:
            plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')
        
        ax.grid(True, alpha=0.3, axis='y')
        ax.set_axisbelow(True)
    
    def _plot_annual_comparison(self, df, ax):
        """Comparación anual mejorada"""
        annual_stats = df.groupby('year')['value'].agg(['mean', 'std', 'count'])
        years = annual_stats.index
        
        # Plot con barras de error mejoradas
        ax.errorbar(years, annual_stats['mean'],
                   yerr=annual_stats['std'],
                   fmt='o-', capsize=6, capthick=2,
                   alpha=0.8, markersize=8, linewidth=2.5,
                   color='darkgreen', markerfacecolor='lightgreen',
                   markeredgecolor='darkgreen', markeredgewidth=1.5,
                   ecolor='gray', elinewidth=1.5)
        
        # Valores sobre los puntos (sin solaparse)
        for i, (year, mean_val) in enumerate(zip(years, annual_stats['mean'])):
            offset = annual_stats['std'].iloc[i] + 0.01
            ax.annotate(f'{mean_val:.3f}',
                       (year, mean_val + offset),
                       ha='center', fontsize=9,
                       bbox=dict(boxstyle='round,pad=0.3',
                                facecolor='white', alpha=0.7))
        
        ax.set_xlabel('Year', fontsize=11, fontweight='bold')
        ax.set_ylabel('Mean ± Std Dev', fontsize=11, fontweight='bold')
        ax.set_title('Annual Variation', fontsize=12, fontweight='bold', pad=10)
        ax.grid(True, alpha=0.3)
        ax.set_axisbelow(True)
        
        # Ajustar ticks
        from matplotlib.ticker import MaxNLocator # type: ignore
        ax.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=min(len(years), 6)))
    
    def _plot_trend_summary(self, df, ax):
        """Resumen de tendencias mejorado"""
        trend_results = self.analyze_trend(df, method='all')
        
        # Formato más compacto y claro
        summary_lines = []
        summary_lines.append("TREND ANALYSIS")
        summary_lines.append("═" * 15)
        summary_lines.append("")
        
        if 'linear' in trend_results:
            lr = trend_results['linear']
            summary_lines.append("LINEAR")
            summary_lines.append(f"Slope: {lr['slope']:.2e}")
            summary_lines.append(f"R²: {lr['r_squared']:.3f}")
            summary_lines.append(f"p: {lr['p_value']:.3f}")
            
            if lr['p_value'] < 0.05:
                summary_lines.append("✓ Significant")
            else:
                summary_lines.append("✗ Not signif.")
            summary_lines.append("")
        
        if 'mann_kendall' in trend_results:
            mk = trend_results['mann_kendall']
            summary_lines.append("MANN-KENDALL")
            summary_lines.append(f"τ: {mk['tau']:.3f}")
            summary_lines.append(f"p: {mk['p_value']:.3f}")
            
            if mk['p_value'] < 0.05:
                if mk['trend'] > 0:
                    summary_lines.append("↗ Increasing")
                else:
                    summary_lines.append("↘ Decreasing")
            else:
                summary_lines.append("→ No trend")
        
        ax.text(0.1, 0.95, '\n'.join(summary_lines),
               transform=ax.transAxes, verticalalignment='top',
               fontsize=10, family='monospace',
               bbox=dict(boxstyle='round,pad=0.5', facecolor='lightblue',
                        alpha=0.1, edgecolor='blue', linewidth=1))
        
        ax.set_title('Trend Summary', fontsize=12, fontweight='bold', pad=10)
        ax.axis('off')
    
    def _plot_autocorrelation(self, df, ax):
        """Autocorrelación mejorada"""
        try:
            values = df['value'].values
            n_lags = min(20, len(values) // 3)
            
            correlations = []
            for lag in range(n_lags):
                if lag == 0:
                    correlations.append(1.0)
                else:
                    correlation = np.corrcoef(values[:-lag], values[lag:])[0, 1]
                    correlations.append(correlation)
            
            # Barras con colores según significancia
            colors = ['red' if abs(c) > 1.96/np.sqrt(len(values)) else 'blue' 
                     for c in correlations]
            ax.bar(range(n_lags), correlations, color=colors, alpha=0.7,
                  edgecolor='black', linewidth=1)
            
            # Líneas de significancia
            significance = 1.96/np.sqrt(len(values))
            ax.axhline(y=0, color='k', linestyle='-', linewidth=1)
            ax.axhline(y=significance, color='r', linestyle='--', 
                      alpha=0.5, label=f'95% CI')
            ax.axhline(y=-significance, color='r', linestyle='--', alpha=0.5)
            
            ax.set_xlabel('Lag', fontsize=11, fontweight='bold')
            ax.set_ylabel('Correlation', fontsize=11, fontweight='bold')
            ax.set_title('Autocorrelation', fontsize=12, fontweight='bold', pad=10)
            ax.legend(fontsize=9, loc='upper right')
            ax.grid(True, alpha=0.3)
            ax.set_axisbelow(True)
            
        except Exception as e:
            ax.text(0.5, 0.5, 'ACF not available',
                   transform=ax.transAxes, ha='center', va='center')
            ax.set_title('Autocorrelation', fontsize=12, fontweight='bold')
            ax.axis('off')
    
    def _plot_distribution(self, df, ax):
        """Distribución mejorada"""
        values = df['value'].values
        
        # Histograma mejorado
        n, bins, patches = ax.hist(values, bins=20, density=True,
                                   alpha=0.7, edgecolor='black',
                                   linewidth=1.2, color='skyblue')
        
        # Ajuste normal
        mu, std = stats.norm.fit(values)
        xmin, xmax = ax.get_xlim()
        x = np.linspace(xmin, xmax, 100)
        p = stats.norm.pdf(x, mu, std)
        ax.plot(x, p, 'r-', linewidth=2.5, label='Normal fit')
        
        # Box de estadísticas más compacto
        cv = std/mu if mu != 0 else 0
        stats_text = f'μ={mu:.3f}\nσ={std:.3f}\nCV={cv:.3f}'
        ax.text(0.72, 0.95, stats_text, transform=ax.transAxes,
               verticalalignment='top', fontsize=10,
               bbox=dict(boxstyle='round,pad=0.4', facecolor='wheat',
                        alpha=0.7, edgecolor='orange'))
        
        ax.set_xlabel('Value', fontsize=11, fontweight='bold')
        ax.set_ylabel('Density', fontsize=11, fontweight='bold')
        ax.set_title('Value Distribution', fontsize=12, fontweight='bold', pad=10)
        ax.legend(fontsize=10, loc='upper left')
        ax.grid(True, alpha=0.3, axis='y')
        ax.set_axisbelow(True)
    
    def _plot_data_quality(self, df, ax):
        """Calidad de datos mejorada"""
        total_points = len(df)
        valid_points = df['value'].notna().sum()
        completeness = valid_points/total_points*100
        
        # Texto compacto y claro
        quality_text = []
        quality_text.append("DATA QUALITY")
        quality_text.append("─" * 12)
        quality_text.append(f"Total: {total_points}")
        quality_text.append(f"Valid: {valid_points}")
        quality_text.append(f"Complete: {completeness:.0f}%")
        quality_text.append("")
        quality_text.append("COVERAGE")
        quality_text.append("─" * 12)
        quality_text.append(f"Years: {df['year'].nunique()}")
        quality_text.append(f"Days: {(df['date'].max() - df['date'].min()).days}")
        
        ax.text(0.1, 0.95, '\n'.join(quality_text),
               transform=ax.transAxes, fontsize=10,
               family='monospace', verticalalignment='top',
               bbox=dict(boxstyle='round,pad=0.4',
                        facecolor='lightyellow', alpha=0.2))
        
        ax.set_title('Data Quality', fontsize=12, fontweight='bold', pad=10)
        ax.axis('off')
    
    def _plot_seasonal_statistics(self, df, ax):
        """Estadísticas estacionales mejoradas"""
        try:
            if len(df) < 12:
                ax.text(0.5, 0.5, 'Insufficient\ndata',
                       ha='center', va='center', transform=ax.transAxes)
                ax.set_title('Seasonal Stats', fontsize=12, fontweight='bold')
                ax.axis('off')
                return
            
            seasonal_means = df.groupby('season')['value'].mean()
            
            if len(seasonal_means) > 0:
                # Colores temáticos para estaciones
                colors = {
                    'spring': '#90EE90',  # Verde claro
                    'summer': '#FFD700',  # Dorado
                    'autumn': '#FF8C00',  # Naranja oscuro
                    'winter': '#87CEEB'   # Azul cielo
                }
                
                season_colors = [colors.get(s, 'gray') for s in seasonal_means.index]
                
                # Barras mejoradas
                bars = ax.bar(range(len(seasonal_means)), seasonal_means.values,
                             color=season_colors, alpha=0.8,
                             edgecolor='black', linewidth=1.5)
                
                # Valores sobre las barras
                for i, (season, value) in enumerate(seasonal_means.items()):
                    ax.text(i, value + 0.005, f'{value:.3f}',
                           ha='center', fontsize=9, fontweight='bold')
                
                # Labels mejorados
                season_labels = {
                    'spring': 'Spring',
                    'summer': 'Summer',
                    'autumn': 'Autumn',
                    'winter': 'Winter'
                }
                
                ax.set_xticks(range(len(seasonal_means)))
                ax.set_xticklabels([season_labels.get(s, s.capitalize()) 
                                   for s in seasonal_means.index],
                                  fontsize=10)
                
                ax.set_ylabel('Mean Value', fontsize=11, fontweight='bold')
                ax.set_title('Seasonal Statistics', fontsize=12, fontweight='bold', pad=10)
                ax.grid(True, alpha=0.3, axis='y')
                ax.set_axisbelow(True)
                
        except Exception as e:
            ax.text(0.5, 0.5, 'Error in\nanalysis',
                   ha='center', va='center', transform=ax.transAxes)
            ax.set_title('Seasonal Statistics', fontsize=12, fontweight='bold')
            ax.axis('off')
    
    def _plot_phenology_summary(self, df, ax):
        """Resumen fenológico mejorado"""
        try:
            phenology_results = self.extract_phenology_metrics(df, method='threshold')
            
            if phenology_results:
                pheno_df = pd.DataFrame(phenology_results).T
                
                # Texto más compacto y organizado
                summary_text = []
                summary_text.append("PHENOLOGY")
                summary_text.append("─" * 12)
                
                metrics = {
                    'sos': ('SOS', 'day'),
                    'pos': ('POS', 'day'),
                    'eos': ('EOS', 'day'),
                    'los': ('LOS', 'days'),
                    'amplitude': ('Ampl', 'val')
                }
                
                for metric, (label, unit) in metrics.items():
                    if metric in pheno_df.columns:
                        values = pheno_df[metric].dropna()
                        if len(values) > 0:
                            mean_val = values.mean()
                            if unit == 'day' or unit == 'days':
                                summary_text.append(f"{label}: {mean_val:.0f} {unit}")
                            else:
                                summary_text.append(f"{label}: {mean_val:.3f}")
                
                ax.text(0.1, 0.95, '\n'.join(summary_text),
                       transform=ax.transAxes, verticalalignment='top',
                       fontsize=10, family='monospace',
                       bbox=dict(boxstyle='round,pad=0.4',
                                facecolor='lightgreen', alpha=0.1))
            else:
                ax.text(0.5, 0.5, 'No phenology\ndata',
                       transform=ax.transAxes, ha='center', va='center')
        except:
            ax.text(0.5, 0.5, 'Phenology\nnot available',
                   transform=ax.transAxes, ha='center', va='center')
        
        ax.set_title('Phenology Summary', fontsize=12, fontweight='bold', pad=10)
        ax.axis('off')
    
    # ========== PLOT METHODS - PHENOLOGY ANALYSIS ==========
    
    def _plot_time_series_with_phenology(self, df: pd.DataFrame, phenology_results: Dict,
                                     ax: plt.Axes, method: str):
        """Plot time series with phenological markers and smoothed curves."""
        
        # 1. OBSERVED DATA
        ax.plot(df['date'], df['value'], 'o-', alpha=0.6, label='Observed', 
            markersize=3, color='lightblue', linewidth=1.2)
        
        # 2. SMOOTHED CURVES for each year
        colors = plt.cm.Set3(np.linspace(0, 1, len(phenology_results)))
        
        # Add one example smoothed curve to legend
        legend_added = False
        
        for i, (year, metrics) in enumerate(phenology_results.items()):
            year_data = df[df['year'] == year].sort_values('doy')
            
            if len(year_data) > 3:
                values = year_data['value'].values
                dates = year_data['date'].values
                
                # Smoothing
                try:
                    from scipy import signal
                    window_length = min(len(values), 7)
                    if window_length % 2 == 0:
                        window_length -= 1
                    if window_length >= 3:
                        values_smooth = signal.savgol_filter(values, window_length, 3)
                        
                        # Add label only to first smoothed curve
                        label = 'Smoothed curves' if not legend_added else ""
                        if not legend_added:
                            legend_added = True
                        
                        ax.plot(dates, values_smooth, '-', 
                            color=colors[i], linewidth=2.5, alpha=0.8, 
                            label=label, zorder=3)
                    else:
                        values_smooth = values
                except:
                    values_smooth = values
                
                # PHENOLOGICAL MARKERS
                if not np.isnan(metrics.get('sos', np.nan)):
                    sos_date = datetime(year, 1, 1) + timedelta(days=int(metrics['sos'])-1)
                    sos_idx = np.argmin(np.abs(year_data['doy'].values - metrics['sos']))
                    if sos_idx < len(year_data):
                        # Add SOS to legend only once
                        sos_label = 'SOS (Start)' if i == 0 else ""
                        ax.plot(sos_date, year_data.iloc[sos_idx]['value'], 
                            marker='o', markersize=6, color=colors[i], 
                            markeredgecolor='white', markeredgewidth=1, 
                            zorder=4, label=sos_label)
                
                if not np.isnan(metrics.get('pos', np.nan)):
                    pos_date = datetime(year, 1, 1) + timedelta(days=int(metrics['pos'])-1)
                    pos_idx = np.argmin(np.abs(year_data['doy'].values - metrics['pos']))
                    if pos_idx < len(year_data):
                        # Add POS to legend only once
                        pos_label = 'POS (Peak)' if i == 0 else ""
                        ax.plot(pos_date, year_data.iloc[pos_idx]['value'],
                            marker='s', markersize=6, color=colors[i],
                            markeredgecolor='white', markeredgewidth=1, 
                            zorder=4, label=pos_label)
                
                if not np.isnan(metrics.get('eos', np.nan)):
                    eos_date = datetime(year, 1, 1) + timedelta(days=int(metrics['eos'])-1)
                    eos_idx = np.argmin(np.abs(year_data['doy'].values - metrics['eos']))
                    if eos_idx < len(year_data):
                        # Add EOS to legend only once
                        eos_label = 'EOS (End)' if i == 0 else ""
                        ax.plot(eos_date, year_data.iloc[eos_idx]['value'],
                            marker='^', markersize=6, color=colors[i],
                            markeredgecolor='white', markeredgewidth=1, 
                            zorder=4, label=eos_label)
        
        # CONFIGURATION (English)
        ax.set_xlabel('Date', fontsize=12)
        ax.set_ylabel(f'{self.index.upper()} Value', fontsize=12)
        ax.set_title(f'Time Series with Phenological Analysis ({method.capitalize()})', 
                    fontsize=14, fontweight='bold')
        
        # LEGEND INSIDE - UPPER LEFT (won't interfere with data)
        ax.legend(loc='upper left', fontsize=9, frameon=True, framealpha=0.9)
        ax.grid(True, alpha=0.3)
        
        # INFO BOX in bottom right
        info_text = f"Method: {method.capitalize()}\nSmoothing: Savitzky-Golay\n○ SOS  □ POS  △ EOS"
        ax.text(0.98, 0.02, info_text, transform=ax.transAxes,
            verticalalignment='bottom', horizontalalignment='right', 
            fontsize=9, bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.8))
    
    def _plot_phenology_info(self, ax, method, threshold, n_years):
        """Plot phenology analysis information panel."""
        info_lines = []
        info_lines.append("PHENOLOGICAL ANALYSIS")
        info_lines.append("")
        info_lines.append(f"Method: {method.upper()}")
        
        if method == 'threshold':
            info_lines.append(f"Threshold: {threshold}%")
        
        info_lines.append(f"Years: {n_years}")
        info_lines.append(f"Periods: {self.periods}/year")
        info_lines.append("")
        info_lines.append("SYMBOLS:")
        info_lines.append("○ SOS (Start)")
        info_lines.append("□ POS (Peak)")
        info_lines.append("△ EOS (End)")
        
        info_text = '\n'.join(info_lines)
        
        ax.text(0.1, 0.95, info_text, transform=ax.transAxes,
                fontsize=8, verticalalignment='top',
                bbox=dict(boxstyle='round,pad=0.4', facecolor='lightgray', alpha=0.2))
        
        ax.set_title('Analysis Info', fontsize=11, fontweight='bold', pad=10)
        ax.axis('off')
    
    def _plot_phenology_timing(self, phenology_results: Dict, ax: plt.Axes):
        """Plot phenological timing metrics."""
        if not phenology_results:
            ax.text(0.5, 0.5, 'No timing\ndata', ha='center', va='center', transform=ax.transAxes)
            ax.set_title('Phenological Timing')
            ax.axis('off')
            return
        
        years = list(phenology_results.keys())
        sos_values = [phenology_results[year].get('sos', np.nan) for year in years]
        pos_values = [phenology_results[year].get('pos', np.nan) for year in years]
        eos_values = [phenology_results[year].get('eos', np.nan) for year in years]
        
        ax.plot(years, sos_values, 'o-', label='SOS', markersize=6, linewidth=2)
        ax.plot(years, pos_values, 's-', label='POS', markersize=6, linewidth=2)  
        ax.plot(years, eos_values, '^-', label='EOS', markersize=6, linewidth=2)
        
        ax.set_xlabel('Year', fontsize=10, fontweight='bold')
        ax.set_ylabel('Day of Year', fontsize=10, fontweight='bold')
        ax.set_title('Phenological Timing', fontsize=11, fontweight='bold')
        ax.legend(fontsize=8, loc='best')
        ax.grid(True, alpha=0.3)
        ax.set_xticks(years)
        ax.tick_params(axis='both', labelsize=9)
        
        if len(years) > 4:
            plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')
    
    def _plot_phenology_amplitude(self, phenology_results: Dict, ax: plt.Axes):
        """Plot amplitude and peak value metrics."""
        years = list(phenology_results.keys())
        amplitude = [phenology_results[year].get('amplitude', np.nan) for year in years]
        peak_values = [phenology_results[year].get('peak_value', np.nan) for year in years]
        
        ax2 = ax.twinx()
        
        line1 = ax.plot(years, amplitude, 'o-', color='green', label='Amplitude', markersize=6)
        line2 = ax2.plot(years, peak_values, 's-', color='darkgreen', label='Peak Value', markersize=6)
        
        ax.set_xlabel('Year')
        ax.set_ylabel('Amplitude', color='green')
        ax2.set_ylabel('Peak Value', color='darkgreen')
        ax.set_title('Seasonal Amplitude & Peak')
        
        # Combine legends
        lines = line1 + line2
        labels = [l.get_label() for l in lines]
        ax.legend(lines, labels, loc='upper left')
        
        ax.grid(True, alpha=0.3)
    
    def _plot_phenology_rates(self, phenology_results: Dict, ax: plt.Axes):
        """Plot growth and senescence rates."""
        years = list(phenology_results.keys())
        growth_rates = [phenology_results[year].get('growth_rate', np.nan) for year in years]
        senescence_rates = [abs(phenology_results[year].get('senescence_rate', np.nan)) for year in years]
        
        ax.plot(years, growth_rates, 'o-', label='Growth Rate', markersize=6, color='lightgreen')
        ax.plot(years, senescence_rates, 's-', label='Senescence Rate', markersize=6, color='orange')
        
        ax.set_xlabel('Year')
        ax.set_ylabel('Rate (value/day)')
        ax.set_title('Growth & Senescence Rates')
        ax.legend()
        ax.grid(True, alpha=0.3)
    
    def _plot_phenology_duration(self, phenology_results: Dict, ax: plt.Axes):
        """
        NUEVO: Duración de temporada (Length of Season - LOS)
        """
        if not phenology_results:
            ax.text(0.5, 0.5, 'No duration\ndata', 
                ha='center', va='center', transform=ax.transAxes)
            ax.set_title('Season Duration')
            ax.axis('off')
            return
        
        years = list(phenology_results.keys())
        los_values = [phenology_results[year].get('los', np.nan) for year in years]
        
        valid_data = [(year, los) for year, los in zip(years, los_values) if not np.isnan(los)]
        
        if not valid_data:
            ax.text(0.5, 0.5, 'No valid\nduration data', 
                ha='center', va='center', transform=ax.transAxes)
            ax.set_title('Season Duration')
            ax.axis('off')
            return
        
        valid_years, valid_los = zip(*valid_data)
        
        bars = ax.bar(range(len(valid_years)), valid_los, alpha=0.7, 
                    color='lightcoral', edgecolor='darkred')
        
        ax.set_xticks(range(len(valid_years)))
        ax.set_xticklabels([str(y) for y in valid_years], fontsize=8, rotation=45)
        
        # Valores encima de las barras
        for i, los in enumerate(valid_los):
            ax.annotate(f'{los:.0f}d', (i, los), 
                    textcoords="offset points", 
                    xytext=(0,3), ha='center', fontsize=7)
        
        # Línea promedio
        mean_los = np.mean(valid_los)
        ax.axhline(y=mean_los, color='red', linestyle='--', alpha=0.7, linewidth=1.5)
        ax.text(0.02, 0.98, f'Average: {mean_los:.0f} days', 
            transform=ax.transAxes, fontsize=7,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        ax.set_ylabel('Days', fontsize=9, fontweight='bold')
        ax.set_title('Season Duration (LOS)', fontsize=10, fontweight='bold')
        ax.grid(True, alpha=0.3, axis='y')
        ax.tick_params(axis='both', labelsize=8)
    
    def _plot_phenology_statistics(self, phenology_results: Dict, ax: plt.Axes):
        """Plot phenological statistics summary."""
        # Create summary statistics
        pheno_df = pd.DataFrame(phenology_results).T
        
        # Calculate means for key metrics
        metrics_to_show = ['sos', 'pos', 'eos', 'los', 'amplitude']
        means = []
        stds = []
        metric_labels = []
        
        for metric in metrics_to_show:
            if metric in pheno_df.columns:
                values = pheno_df[metric].dropna()
                if len(values) > 0:
                    means.append(values.mean())
                    stds.append(values.std())
                    metric_labels.append(metric.upper())
        
        if means:
            y_pos = np.arange(len(metric_labels))
            ax.barh(y_pos, means, xerr=stds, alpha=0.7, capsize=5)
            ax.set_yticks(y_pos)
            ax.set_yticklabels(metric_labels)
            ax.set_xlabel('Value')
            ax.set_title('Mean ± Std')
            ax.grid(True, alpha=0.3, axis='x')
        else:
            ax.text(0.5, 0.5, 'No statistics\navailable', 
                ha='center', va='center', transform=ax.transAxes)
            ax.set_title('Phenological Statistics')
    
    def _plot_phenology_data_availability(self, phenology_results: Dict, ax: plt.Axes):
        """
        RENOMBRADO: Disponibilidad de datos (antes "completeness").
        Muestra en cuántos años se pudo calcular cada métrica.
        """
        if not phenology_results:
            ax.text(0.5, 0.5, 'No data\navailable', 
                ha='center', va='center', transform=ax.transAxes)
            ax.set_title('Data Availability')
            ax.axis('off')
            return
        
        pheno_df = pd.DataFrame(phenology_results).T
        total_years = len(phenology_results)
        
        # Calcular disponibilidad para métricas clave
        metrics = ['sos', 'pos', 'eos', 'los']
        metric_labels = ['SOS', 'POS', 'EOS', 'LOS']
        availability = []
        
        for metric in metrics:
            if metric in pheno_df.columns:
                valid_count = pheno_df[metric].notna().sum()
                availability.append(valid_count / total_years * 100)
            else:
                availability.append(0)
        
        # Gráfico de barras horizontal
        y_pos = np.arange(len(metric_labels))
        colors = ['red' if a < 50 else 'orange' if a < 80 else 'green' for a in availability]
        
        bars = ax.barh(y_pos, availability, color=colors, alpha=0.7, edgecolor='black')
        
        ax.set_yticks(y_pos)
        ax.set_yticklabels(metric_labels, fontsize=8)
        ax.set_xlabel('% Valid Years', fontsize=8)  # MÁS CLARO
        ax.set_xlim(0, 100)
        ax.set_title('Data Availability', fontsize=10, fontweight='bold')  # MÁS CLARO
        ax.grid(True, alpha=0.3, axis='x')
        
        # Porcentajes en las barras
        for i, avail in enumerate(availability):
            if avail > 10:
                ax.annotate(f'{avail:.0f}%', (avail/2, i), 
                        ha='center', va='center', fontsize=7, 
                        color='white', weight='bold')
        
        ax.tick_params(axis='both', labelsize=8)
    
    def _plot_annual_phenology_comparison(self, df: pd.DataFrame, phenology_results: Dict, ax: plt.Axes):
        """Plot annual phenological curves for comparison."""
        if not phenology_results or len(df) == 0:
            ax.text(0.5, 0.5, 'No data for\ncomparison', 
                ha='center', va='center', transform=ax.transAxes)
            ax.set_title('Annual Comparison')
            ax.axis('off')
            return
    
        # Colors for each year
        years = sorted(phenology_results.keys())
        colors = plt.cm.viridis(np.linspace(0, 1, len(years)))
        
        # Plot curve for each year
        for i, year in enumerate(years):
            year_data = df[df['year'] == year].sort_values('doy')
            
            if len(year_data) > 0:
                # Main line
                ax.plot(year_data['doy'], year_data['value'], 
                    'o-', color=colors[i], label=str(year), 
                    alpha=0.8, markersize=4, linewidth=2)
                
                # Phenological markers
                metrics = phenology_results[year]
                
                # SOS marker
                if not np.isnan(metrics.get('sos', np.nan)):
                    sos_idx = np.argmin(np.abs(year_data['doy'].values - metrics['sos']))
                    if sos_idx < len(year_data):
                        ax.plot(metrics['sos'], year_data.iloc[sos_idx]['value'], 
                            'o', color=colors[i], markersize=8, 
                            markeredgecolor='white', markeredgewidth=2)
                
                # POS marker
                if not np.isnan(metrics.get('pos', np.nan)):
                    pos_idx = np.argmin(np.abs(year_data['doy'].values - metrics['pos']))
                    if pos_idx < len(year_data):
                        ax.plot(metrics['pos'], year_data.iloc[pos_idx]['value'], 
                            's', color=colors[i], markersize=8, 
                            markeredgecolor='white', markeredgewidth=2)
                
                # EOS marker
                if not np.isnan(metrics.get('eos', np.nan)):
                    eos_idx = np.argmin(np.abs(year_data['doy'].values - metrics['eos']))
                    if eos_idx < len(year_data):
                        ax.plot(metrics['eos'], year_data.iloc[eos_idx]['value'], 
                            '^', color=colors[i], markersize=8, 
                            markeredgecolor='white', markeredgewidth=2)
        
        ax.set_xlabel('Day of Year', fontsize=11, fontweight='bold')
        ax.set_ylabel(f'{self.index.upper()} Value', fontsize=11, fontweight='bold')
        ax.set_title('Annual Phenological Curves Comparison', fontsize=12, fontweight='bold')
        
        # LEGEND INSIDE - UPPER RIGHT
        ax.legend(loc='upper right', fontsize=9, frameon=True, framealpha=0.9)
        ax.grid(True, alpha=0.3)
        
        # Info about symbols in bottom left
        ax.text(0.02, 0.02, 'Symbols: ○ SOS, □ POS, △ EOS', 
            transform=ax.transAxes, fontsize=8,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    def _plot_phenology_metrics_summary(self, phenology_results: Dict, ax: plt.Axes):
        """Panel compacto con métricas clave."""
        if not phenology_results:
            ax.text(0.5, 0.5, 'No phenology\ndata available', 
                ha='center', va='center', transform=ax.transAxes)
            ax.set_title('Metrics Summary')
            ax.axis('off')
            return
        
        # Convertir a DataFrame
        pheno_df = pd.DataFrame(phenology_results).T
        
        # Crear resumen compacto
        summary_lines = []
        summary_lines.append("KEY METRICS")
        summary_lines.append("=" * 15)
        
        # Métricas principales
        metrics = ['sos', 'pos', 'eos', 'amplitude']
        labels = ['SOS (day)', 'POS (day)', 'EOS (day)', 'Amplitude']
        
        for metric, label in zip(metrics, labels):
            if metric in pheno_df.columns:
                values = pheno_df[metric].dropna()
                if len(values) > 0:
                    mean_val = values.mean()
                    if 'day' in label:
                        summary_lines.append(f"{label}: {mean_val:.0f}")
                    else:
                        summary_lines.append(f"{label}: {mean_val:.3f}")
        
        # Añadir años analizados
        summary_lines.append("")
        summary_lines.append("YEARS")
        summary_lines.append("=" * 15)
        years = sorted(phenology_results.keys())
        for year in years:
            summary_lines.append(f"✓ {year}")
        
        ax.text(0.05, 0.95, '\n'.join(summary_lines), 
            transform=ax.transAxes, fontsize=9, family='monospace',
            verticalalignment='top')
        ax.set_title('Summary', fontsize=11, fontweight='bold')
        ax.axis('off')

    def _plot_phenology_quality_summary(self, phenology_results: Dict, ax: plt.Axes):
        """
        SIMPLIFICADO: Resumen de calidad fenológica en inglés.
        """
        if not phenology_results:
            ax.text(0.5, 0.5, 'No quality\ndata', ha='center', va='center', transform=ax.transAxes)
            ax.set_title('Quality Summary')
            ax.axis('off')
            return
        
        pheno_df = pd.DataFrame(phenology_results).T
        total_years = len(phenology_results)
        
        # Texto en inglés
        quality_text = []
        quality_text.append("ANALYSIS QUALITY")
        quality_text.append("=" * 16)
        quality_text.append(f"Years analyzed: {total_years}")
        quality_text.append("")
        
        # Años con datos completos (SOS, POS, EOS disponibles)
        complete_years = 0
        for year in phenology_results:
            year_metrics = phenology_results[year]
            if (not pd.isna(year_metrics.get('sos', np.nan)) and 
                not pd.isna(year_metrics.get('pos', np.nan)) and
                not pd.isna(year_metrics.get('eos', np.nan))):
                complete_years += 1
        
        quality_text.append("COMPLETE CYCLES")
        quality_text.append("=" * 16)
        quality_text.append(f"Full cycles: {complete_years}/{total_years}")
        quality_text.append(f"Success rate: {complete_years/total_years*100:.0f}%")
        quality_text.append("")
        
        # Variabilidad de duración de temporada
        if 'los' in pheno_df.columns:
            los_values = pheno_df['los'].dropna()
            if len(los_values) > 1:
                cv_los = los_values.std() / los_values.mean() * 100
                quality_text.append("VARIABILITY")
                quality_text.append("=" * 16)
                quality_text.append(f"LOS CV: {cv_los:.1f}%")
        
        # Mostrar texto
        ax.text(0.05, 0.95, '\n'.join(quality_text), 
            transform=ax.transAxes, fontsize=8, family='monospace',
            verticalalignment='top')
        
        ax.set_title('Quality Summary', fontsize=10, fontweight='bold')
        ax.axis('off')
    
    # ========== PHENOLOGY EXTRACTION HELPER METHODS ==========
    
    def _phenology_threshold_method(self, year_data: pd.DataFrame, 
                                threshold_percentile: float,
                                smoothing: bool,
                                min_season_length: int) -> Dict[str, float]:
        """
        Threshold-based phenology extraction using percentile values.
        
        This method identifies phenological phases by finding when vegetation
        values cross specific thresholds relative to the annual range.
        """
        values = year_data['value'].values
        doys = year_data['doy'].values
        
        # Apply smoothing if requested
        if smoothing and len(values) > 3:
            try:
                # Savitzky-Golay filter for smoothing
                window_length = min(len(values), 7)
                if window_length % 2 == 0:
                    window_length -= 1
                values_smooth = signal.savgol_filter(values, window_length, 3)
            except:
                values_smooth = values
        else:
            values_smooth = values
        
        # Calculate basic statistics
        min_val = np.min(values_smooth)
        max_val = np.max(values_smooth)
        amplitude = max_val - min_val
        
        if amplitude < 0.01:  # Very low amplitude - skip analysis
            return None
        
        # Calculate threshold based on percentile
        threshold = min_val + (amplitude * threshold_percentile / 100)
        
        # Find crossings
        above_threshold = values_smooth >= threshold
        
        # Find Start of Season (first crossing above threshold)
        sos_idx = None
        for i in range(len(above_threshold) - 1):
            if not above_threshold[i] and above_threshold[i + 1]:
                sos_idx = i + 1
                break
        
        # Find End of Season (last crossing below threshold)
        eos_idx = None
        for i in range(len(above_threshold) - 1, 0, -1):
            if above_threshold[i - 1] and not above_threshold[i]:
                eos_idx = i - 1
                break
        
        # Find Peak of Season
        pos_idx = np.argmax(values_smooth)
        
        # Calculate metrics
        metrics = {
            'baseline': float(min_val),
            'peak_value': float(max_val),
            'amplitude': float(amplitude),
            'pos': float(doys[pos_idx]),
            'auc': float(np.trapz(values_smooth, doys))
        }
        
        if sos_idx is not None:
            metrics['sos'] = float(doys[sos_idx])
            metrics['sos_value'] = float(values_smooth[sos_idx])
        else:
            metrics['sos'] = np.nan
            metrics['sos_value'] = np.nan
        
        if eos_idx is not None:
            metrics['eos'] = float(doys[eos_idx])  
            metrics['eos_value'] = float(values_smooth[eos_idx])
        else:
            metrics['eos'] = np.nan
            metrics['eos_value'] = np.nan
        
        # Calculate Length of Season
        if not np.isnan(metrics['sos']) and not np.isnan(metrics['eos']):
            metrics['los'] = float(metrics['eos'] - metrics['sos'])
            
            # Calculate growth and senescence rates
            if sos_idx is not None and pos_idx is not None:
                growth_days = doys[pos_idx] - doys[sos_idx]
                if growth_days > 0:
                    metrics['growth_rate'] = float((values_smooth[pos_idx] - values_smooth[sos_idx]) / growth_days)
                else:
                    metrics['growth_rate'] = np.nan
            else:
                metrics['growth_rate'] = np.nan
                
            if pos_idx is not None and eos_idx is not None:
                senescence_days = doys[eos_idx] - doys[pos_idx]
                if senescence_days > 0:
                    metrics['senescence_rate'] = float((values_smooth[eos_idx] - values_smooth[pos_idx]) / senescence_days)
                else:
                    metrics['senescence_rate'] = np.nan
            else:
                metrics['senescence_rate'] = np.nan
        else:
            metrics['los'] = np.nan
            metrics['growth_rate'] = np.nan
            metrics['senescence_rate'] = np.nan
        
        return metrics
    
    def _phenology_derivative_method(self, year_data: pd.DataFrame,
                               smoothing: bool,
                               min_season_length: int) -> Dict[str, float]:
        """
        Método derivative sin usar threshold_percentile.
        """
        values = year_data['value'].values
        doys = year_data['doy'].values
        
        if len(values) < 5:
            print(f"Insufficient data for derivative method: {len(values)} points")
            return None
        
        try:
            # Smoothing
            if smoothing:
                try:
                    window_length = min(len(values), 7)
                    if window_length % 2 == 0:
                        window_length -= 1
                    if window_length >= 3:
                        values_smooth = signal.savgol_filter(values, window_length, 3)
                    else:
                        values_smooth = values
                except Exception as e:
                    print(f"Smoothing failed: {e}")
                    values_smooth = values
            else:
                values_smooth = values
            
            # Calculate derivatives
            try:
                first_deriv = np.gradient(values_smooth, doys)
            except Exception as e:
                print(f"Derivative calculation failed: {e}")
                return None
            
            # Validate derivatives
            if np.any(np.isnan(first_deriv)) or np.any(np.isinf(first_deriv)):
                print("Invalid derivatives calculated")
                return None
            
            # Umbral basado en desviación estándar (NO porcentaje)
            deriv_std = np.std(first_deriv)
            pos_threshold = deriv_std * 0.25
            neg_threshold = -pos_threshold
            
            # SOS: First significant positive derivative
            sos_idx = None
            for i in range(len(first_deriv)):
                if first_deriv[i] > pos_threshold:
                    sos_idx = i
                    break
            
            # POS: Maximum
            pos_idx = np.argmax(values_smooth)
            
            # EOS: First significant negative derivative after peak
            eos_idx = None
            for i in range(pos_idx, len(first_deriv)):
                if first_deriv[i] < neg_threshold:
                    eos_idx = i
                    break
            
            # Calculate metrics
            min_val = float(np.min(values_smooth))
            max_val = float(np.max(values_smooth))
            
            metrics = {
                'baseline': min_val,
                'peak_value': max_val,
                'amplitude': float(max_val - min_val),
                'pos': float(doys[pos_idx]),
                'auc': float(np.trapz(values_smooth, doys))
            }
            
            # Add SOS and EOS
            if sos_idx is not None:
                metrics['sos'] = float(doys[sos_idx])
                metrics['sos_value'] = float(values_smooth[sos_idx])
            else:
                metrics['sos'] = np.nan
                metrics['sos_value'] = np.nan
            
            if eos_idx is not None:
                metrics['eos'] = float(doys[eos_idx])
                metrics['eos_value'] = float(values_smooth[eos_idx])
            else:
                metrics['eos'] = np.nan
                metrics['eos_value'] = np.nan
            
            # Derived metrics
            if not np.isnan(metrics['sos']) and not np.isnan(metrics['eos']):
                metrics['los'] = float(metrics['eos'] - metrics['sos'])
                
                if sos_idx is not None and pos_idx is not None:
                    growth_days = doys[pos_idx] - doys[sos_idx]
                    if growth_days > 0:
                        metrics['growth_rate'] = float((values_smooth[pos_idx] - values_smooth[sos_idx]) / growth_days)
                    else:
                        metrics['growth_rate'] = np.nan
                else:
                    metrics['growth_rate'] = np.nan
                    
                if pos_idx is not None and eos_idx is not None:
                    senescence_days = doys[eos_idx] - doys[pos_idx]
                    if senescence_days > 0:
                        metrics['senescence_rate'] = float((values_smooth[eos_idx] - values_smooth[pos_idx]) / senescence_days)
                    else:
                        metrics['senescence_rate'] = np.nan
                else:
                    metrics['senescence_rate'] = np.nan
            else:
                metrics['los'] = np.nan
                metrics['growth_rate'] = np.nan
                metrics['senescence_rate'] = np.nan
            
            return metrics
            
        except Exception as e:
            print(f"Error in derivative method: {e}")
            return None
    
    def _phenology_logistic_method(self, year_data: pd.DataFrame,
                                smoothing: bool,
                                min_season_length: int) -> Dict[str, float]:
        """
        Logistic curve fitting for phenology extraction.
        
        Fits double logistic curves to represent vegetation growth and senescence,
        providing smooth phenological parameter estimates.
        """
        values = year_data['value'].values
        doys = year_data['doy'].values
        
        if len(values) < 6:
            return None
        
        # Apply smoothing
        if smoothing:
            try:
                window_length = min(len(values), 7)
                if window_length % 2 == 0:
                    window_length -= 1
                values_smooth = signal.savgol_filter(values, window_length, 3)
            except:
                values_smooth = values
        else:
            values_smooth = values
        
        try:
            # Fit double logistic function
            # This is a simplified approach - in practice you'd use more sophisticated fitting
            from scipy.optimize import curve_fit
            
            def double_logistic(x, baseline, amplitude, sos, growth_rate, eos, senescence_rate):
                """Double logistic function for vegetation phenology"""
                growth = amplitude / (1 + np.exp(-growth_rate * (x - sos)))
                senescence = amplitude / (1 + np.exp(senescence_rate * (x - eos)))
                return baseline + growth - senescence
            
            # Initial parameter guess
            baseline_guess = np.min(values_smooth)
            amplitude_guess = np.max(values_smooth) - baseline_guess
            sos_guess = doys[len(doys)//4]  # First quarter
            eos_guess = doys[3*len(doys)//4]  # Last quarter
            
            initial_guess = [baseline_guess, amplitude_guess, sos_guess, 0.1, eos_guess, -0.1]
            
            # Fit the curve
            popt, _ = curve_fit(double_logistic, doys, values_smooth, 
                            p0=initial_guess, maxfev=1000)
            
            baseline, amplitude, sos, growth_rate, eos, senescence_rate = popt
            
            # Calculate peak (maximum of fitted curve)
            fitted_curve = double_logistic(doys, *popt)
            pos_idx = np.argmax(fitted_curve)
            pos = doys[pos_idx]
            peak_value = fitted_curve[pos_idx]
            
            metrics = {
                'baseline': float(baseline),
                'peak_value': float(peak_value),
                'amplitude': float(amplitude),
                'pos': float(pos),
                'sos': float(sos),
                'eos': float(eos),
                'los': float(eos - sos),
                'growth_rate': float(growth_rate),
                'senescence_rate': float(senescence_rate),
                'sos_value': float(double_logistic(sos, *popt)),
                'eos_value': float(double_logistic(eos, *popt)),
                'auc': float(np.trapz(fitted_curve, doys))
            }
            
            return metrics
            
        except Exception as e:
            print(f"Logistic fitting failed: {e}")
            # Fallback to threshold method
            return self._phenology_threshold_method(year_data, 50, smoothing, min_season_length)
    
    # ========== COMPARISON AND ANALYSIS METHODS ==========
    
    def compare_phenology_years(self, 
                            point: Optional[Union[tuple, ee.Geometry.Point]] = None,
                            reference_year: Optional[int] = None) -> Dict[str, Any]:
        """
        Compare phenological metrics across years to identify anomalies and trends.
        
        Parameters
        ----------
        point : location for extraction
        reference_year : year to use as reference (if None, uses mean of all years)
        
        Returns
        -------
        dict
            Comparison results with anomalies and trends
        """
        # Extract time series and phenology
        df = self.extract_time_series(point)
        phenology_results = self.extract_phenology_metrics(df)
        
        if len(phenology_results) < 2:
            return {'error': 'Need at least 2 years for comparison'}
        
        # Convert to DataFrame for easier analysis
        pheno_df = pd.DataFrame(phenology_results).T
        pheno_df.index = pheno_df.index.astype(int)
        
        # Calculate statistics
        stats = {}
        comparison = {}
        
        for metric in pheno_df.columns:
            if pheno_df[metric].notna().sum() > 1:  # Need at least 2 valid values
                values = pheno_df[metric].dropna()
                
                stats[metric] = {
                    'mean': float(values.mean()),
                    'std': float(values.std()),
                    'min': float(values.min()),
                    'max': float(values.max()),
                    'cv': float(values.std() / values.mean() * 100) if values.mean() != 0 else np.nan
                }
                
                # Calculate anomalies (z-scores)
                z_scores = (values - values.mean()) / values.std() if values.std() > 0 else values * 0
                
                comparison[metric] = {
                    'values': values.to_dict(),
                    'anomalies': z_scores.to_dict(),
                    'trend': self._calculate_trend(values.index, values.values)
                }
        
        # Identify extreme years
        extreme_years = {}
        for metric in comparison:
            anomalies = comparison[metric]['anomalies']
            extreme_years[metric] = {
                'earliest': min(anomalies, key=anomalies.get) if 'sos' in metric else None,
                'latest': max(anomalies, key=anomalies.get) if 'sos' in metric or 'eos' in metric else None,
                'highest': max(anomalies, key=anomalies.get) if 'amplitude' in metric or 'peak' in metric else None,
                'lowest': min(anomalies, key=anomalies.get) if 'amplitude' in metric or 'peak' in metric else None
            }
        
        return {
            'statistics': stats,
            'comparisons': comparison,
            'extreme_years': extreme_years,
            'summary': self._generate_phenology_summary(stats, comparison)
        }
    
    def _calculate_trend(self, years: np.ndarray, values: np.ndarray) -> Dict[str, float]:
        """Calculate trend statistics for phenological metrics."""
        if len(values) < 3:
            return {'slope': np.nan, 'p_value': np.nan, 'r_squared': np.nan}
        
        try:
            slope, intercept, r_value, p_value, std_err = stats.linregress(years, values)
            return {
                'slope': float(slope),
                'p_value': float(p_value),
                'r_squared': float(r_value**2)
            }
        except:
            return {'slope': np.nan, 'p_value': np.nan, 'r_squared': np.nan}
    
    def _generate_phenology_summary(self, stats: Dict, comparisons: Dict) -> str:
        """Generate text summary of phenological analysis."""
        summary_lines = []
        
        # Start of Season analysis
        if 'sos' in stats:
            sos_stats = stats['sos']
            sos_trend = comparisons['sos']['trend']
            
            summary_lines.append("START OF SEASON (SOS):")
            summary_lines.append(f"  Average: Day {sos_stats['mean']:.1f} (±{sos_stats['std']:.1f})")
            summary_lines.append(f"  Range: Day {sos_stats['min']:.1f} to {sos_stats['max']:.1f}")
            summary_lines.append(f"  Variability: {sos_stats['cv']:.1f}% CV")
            
            if not np.isnan(sos_trend['slope']) and sos_trend['p_value'] < 0.05:
                trend_direction = "earlier" if sos_trend['slope'] < 0 else "later"
                summary_lines.append(f"  Trend: {abs(sos_trend['slope']):.2f} days/year {trend_direction} (p={sos_trend['p_value']:.3f})")
            else:
                summary_lines.append("  Trend: No significant trend")
            summary_lines.append("")
        
        # Peak of Season analysis
        if 'pos' in stats:
            pos_stats = stats['pos']
            pos_trend = comparisons['pos']['trend']
            
            summary_lines.append("PEAK OF SEASON (POS):")
            summary_lines.append(f"  Average: Day {pos_stats['mean']:.1f} (±{pos_stats['std']:.1f})")
            summary_lines.append(f"  Range: Day {pos_stats['min']:.1f} to {pos_stats['max']:.1f}")
            
            if not np.isnan(pos_trend['slope']) and pos_trend['p_value'] < 0.05:
                trend_direction = "earlier" if pos_trend['slope'] < 0 else "later"
                summary_lines.append(f"  Trend: {abs(pos_trend['slope']):.2f} days/year {trend_direction}")
            summary_lines.append("")
        
        # End of Season analysis
        if 'eos' in stats:
            eos_stats = stats['eos']
            eos_trend = comparisons['eos']['trend']
            
            summary_lines.append("END OF SEASON (EOS):")
            summary_lines.append(f"  Average: Day {eos_stats['mean']:.1f} (±{eos_stats['std']:.1f})")
            summary_lines.append(f"  Range: Day {eos_stats['min']:.1f} to {eos_stats['max']:.1f}")
            
            if not np.isnan(eos_trend['slope']) and eos_trend['p_value'] < 0.05:
                trend_direction = "earlier" if eos_trend['slope'] < 0 else "later"
                summary_lines.append(f"  Trend: {abs(eos_trend['slope']):.2f} days/year {trend_direction}")
            summary_lines.append("")
        
        # Length of Season analysis
        if 'los' in stats:
            los_stats = stats['los']
            los_trend = comparisons['los']['trend']
            
            summary_lines.append("LENGTH OF SEASON (LOS):")
            summary_lines.append(f"  Average: {los_stats['mean']:.1f} days (±{los_stats['std']:.1f})")
            summary_lines.append(f"  Range: {los_stats['min']:.1f} to {los_stats['max']:.1f} days")
            
            if not np.isnan(los_trend['slope']) and los_trend['p_value'] < 0.05:
                trend_direction = "lengthening" if los_trend['slope'] > 0 else "shortening"
                summary_lines.append(f"  Trend: {abs(los_trend['slope']):.2f} days/year {trend_direction}")
            summary_lines.append("")
        
        # Amplitude analysis
        if 'amplitude' in stats:
            amp_stats = stats['amplitude']
            amp_trend = comparisons['amplitude']['trend']
            
            summary_lines.append("SEASONAL AMPLITUDE:")
            summary_lines.append(f"  Average: {amp_stats['mean']:.3f} (±{amp_stats['std']:.3f})")
            summary_lines.append(f"  Range: {amp_stats['min']:.3f} to {amp_stats['max']:.3f}")
            
            if not np.isnan(amp_trend['slope']) and amp_trend['p_value'] < 0.05:
                trend_direction = "increasing" if amp_trend['slope'] > 0 else "decreasing"
                summary_lines.append(f"  Trend: {trend_direction} ({abs(amp_trend['slope']):.4f}/year)")
        
        return "\n".join(summary_lines)
    

class SpatialTrendAnalyzer:
    """
    Spatial trend analysis for generating pixel-wise trend maps.
    
    This class complements TimeSeriesAnalyzer by providing spatial analysis
    capabilities using Earth Engine's distributed computing.
    """
    
    def __init__(self, ndvi_seasonality_instance):
        """
        Initialize SpatialTrendAnalyzer.
        
        Parameters
        ----------
        ndvi_seasonality_instance : NdviSeasonality
            Configured instance with ROI, periods, years, etc.
        """
        self.processor = ndvi_seasonality_instance
        
    def calculate_pixel_trends(self,
                            method: str = 'linear',
                            min_observations: int = 5,
                            export: bool = False,
                            scale: int = 30) -> ee.Image:
        """
        Calculate per-pixel temporal trends across the ROI.
        
        Parameters
        ----------
        method : str
            Trend calculation method: 'linear', 'sen', 'mann_kendall'
        min_observations : int
            Minimum valid observations required per pixel
        export : bool
            Whether to export the result
        scale : int
            Output scale in meters
            
        Returns
        -------
        ee.Image
            Multi-band image with trend statistics
        """
        print(f"Calculating {method} trend map...")
        
        # Create collection with time bands
        images_list = []
        
        for year_idx, year in enumerate(range(self.processor.start_year, self.processor.end_year)):
            for period_idx in range(self.processor.periods):
                # Get composite
                composite = self.processor.get_period_composite(year, period_idx)
                
                # Add time band (fractional years)
                time_value = year_idx + (period_idx / self.processor.periods)
                time_band = ee.Image.constant(time_value).float().rename('time')
                
                # Add constant band for regression
                constant = ee.Image.constant(1).rename('constant')
                
                # Combine bands
                composite_with_vars = composite.select(['nd']).addBands([constant, time_band])
                
                images_list.append(composite_with_vars)
        
        # Create collection
        collection = ee.ImageCollection.fromImages(images_list)
        
        if method == 'linear':
            # Linear regression
            linear_fit = collection.select(['constant', 'time', 'nd']).reduce(
                ee.Reducer.linearRegression(
                    numX=2,  # constant and time
                    numY=1   # nd value
                )
            )
            
            # Extract coefficients
            coefficients = linear_fit.select('coefficients').arrayProject([0]).arrayFlatten([['constant', 'time']])
            
            # Rename for clarity
            intercept = coefficients.select('constant').rename('intercept')
            slope = coefficients.select('time').rename('slope')
            
            # Calculate magnitude of change
            total_time = self.processor.end_year - self.processor.start_year
            magnitude = slope.multiply(total_time).rename('magnitude')
            
            # Combine all bands
            trend_image = slope.addBands([intercept, magnitude])
            
        elif method == 'sen':
            # Sen's slope
            trend_image = collection.select(['time', 'nd']).reduce(
                ee.Reducer.sensSlope()
            )
            
            # Rename bands
            trend_image = trend_image.select(
                ['offset', 'slope'],
                ['intercept', 'slope']
            )
            
            # Add magnitude
            total_time = self.processor.end_year - self.processor.start_year
            magnitude = trend_image.select('slope').multiply(total_time).rename('magnitude')
            trend_image = trend_image.addBands(magnitude)
            
        else:
            raise ValueError(f"Method '{method}' not supported")
        
        # Mask pixels with insufficient observations
        observation_count = collection.select('nd').count()
        trend_image = trend_image.updateMask(observation_count.gte(min_observations))
        
        # Clip to ROI
        trend_image = trend_image.clip(self.processor.roi)
        
        # Export if requested
        if export:
            filename = f'{self.processor.sat}_{self.processor.index}_trend_{method}_{self.processor.start_year}_{self.processor.end_year-1}.tif'
            print(f"Exporting trend map as {filename}")
            self.processor.get_export_single(trend_image, name=filename, scale=scale)
        
        return trend_image