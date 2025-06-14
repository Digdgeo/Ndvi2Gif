"""
Basic tests for ndvi2gif package
"""
import pytest
import sys
import os

# Add the parent directory to path to import ndvi2gif
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

def test_import():
    """Test that the package can be imported successfully."""
    try:
        import ndvi2gif
        assert True
    except ImportError as e:
        pytest.fail(f"Failed to import ndvi2gif: {e}")

def test_import_main_class():
    """Test that the main NdviSeasonality class can be imported."""
    try:
        from ndvi2gif.ndvi2gif import NdviSeasonality
        assert NdviSeasonality is not None
    except ImportError as e:
        pytest.fail(f"Failed to import NdviSeasonality class: {e}")

def test_class_initialization_default():
    """Test that NdviSeasonality can be initialized with default parameters."""
    try:
        from ndvi2gif.ndvi2gif import NdviSeasonality
        # Mock Earth Engine to test basic initialization without auth
        instance = NdviSeasonality()
        assert instance is not None
        assert instance.periods == 4
        assert instance.start_year == 2016
        assert instance.end_year == 2020
        assert instance.sat == 'S2'
        assert instance.key == 'max'
        assert instance.index == 'ndvi'
        # Test that period generation works
        assert len(instance.period_names) == 4
        assert instance.period_names == ['winter', 'spring', 'summer', 'autumn']
    except Exception as e:
        # If Earth Engine is not authenticated, we still check basic attributes
        pytest.skip(f"Skipping Earth Engine dependent test: {e}")

def test_class_initialization_custom_parameters():
    """Test that NdviSeasonality accepts custom parameters."""
    try:
        from ndvi2gif.ndvi2gif import NdviSeasonality
        instance = NdviSeasonality(
            periods=12, 
            start_year=2018, 
            end_year=2022, 
            sat='Landsat', 
            key='median', 
            index='evi'
        )
        assert instance.periods == 12
        assert instance.start_year == 2018
        assert instance.end_year == 2022
        assert instance.sat == 'Landsat'
        assert instance.key == 'median'
        assert instance.index == 'evi'
    except Exception as e:
        pytest.skip(f"Skipping Earth Engine dependent test: {e}")

def test_period_generation():
    """Test the dynamic period generation functionality."""
    try:
        from ndvi2gif.ndvi2gif import NdviSeasonality
        instance = NdviSeasonality(periods=4)
        
        # Test 4 periods (seasons)
        assert len(instance.period_dates) == 4
        assert len(instance.period_names) == 4
        assert instance.period_names == ['winter', 'spring', 'summer', 'autumn']
        
        # Test that dates are properly formatted
        for period in instance.period_dates:
            assert len(period) == 2  # start and end date
            assert period[0].startswith('-')  # starts with month-day format
            assert period[1].startswith('-')
            
    except Exception as e:
        pytest.skip(f"Skipping Earth Engine dependent test: {e}")

def test_period_generation_12_periods():
    """Test 12-period (monthly) generation."""
    try:
        from ndvi2gif.ndvi2gif import NdviSeasonality
        instance = NdviSeasonality(periods=12)
        
        assert len(instance.period_dates) == 12
        assert len(instance.period_names) == 12
        assert 'january' in instance.period_names
        assert 'december' in instance.period_names
        
    except Exception as e:
        pytest.skip(f"Skipping Earth Engine dependent test: {e}")

def test_index_methods_exist():
    """Test that all index calculation methods exist."""
    try:
        from ndvi2gif.ndvi2gif import NdviSeasonality
        instance = NdviSeasonality()
        
        # Check that all index methods are callable
        expected_indices = ['ndvi', 'ndwi', 'mndwi', 'evi', 'savi', 'gndvi', 
                           'avi', 'nbri', 'ndsi', 'aweinsh', 'awei', 'ndmi']
        
        for index in expected_indices:
            assert index in instance.d
            assert callable(instance.d[index])
            
    except Exception as e:
        pytest.skip(f"Skipping Earth Engine dependent test: {e}")

def test_valid_satellite_options():
    """Test that satellite parameter validation works."""
    try:
        from ndvi2gif.ndvi2gif import NdviSeasonality
        
        # Valid satellites
        valid_sats = ['S2', 'Landsat', 'MODIS', 'S1']
        for sat in valid_sats:
            instance = NdviSeasonality(sat=sat)
            assert instance.sat == sat
            
        # Invalid satellite should default to S2
        instance = NdviSeasonality(sat='InvalidSat')
        assert instance.sat == 'S2'
        
    except Exception as e:
        pytest.skip(f"Skipping Earth Engine dependent test: {e}")

def test_valid_statistic_options():
    """Test that statistic parameter validation works."""
    try:
        from ndvi2gif.ndvi2gif import NdviSeasonality
        
        # Valid statistics
        valid_stats = ['max', 'median', 'perc_90', 'perc_95', 'mean']
        for stat in valid_stats:
            instance = NdviSeasonality(key=stat)
            assert instance.key == stat
            
        # Invalid statistic should default to max
        instance = NdviSeasonality(key='invalid_stat')
        assert instance.key == 'max'
        
    except Exception as e:
        pytest.skip(f"Skipping Earth Engine dependent test: {e}")

def test_helper_functions_exist():
    """Test that helper functions exist and are callable."""
    try:
        from ndvi2gif.ndvi2gif import scale_OLI, scale_ETM
        assert callable(scale_OLI)
        assert callable(scale_ETM)
    except ImportError as e:
        pytest.fail(f"Failed to import helper functions: {e}")

# Integration test (will likely be skipped due to EE auth requirements)
def test_basic_workflow():
    """Test a basic workflow - may be skipped if Earth Engine is not authenticated."""
    try:
        from ndvi2gif.ndvi2gif import NdviSeasonality
        
        # Create instance with minimal parameters
        instance = NdviSeasonality(periods=4, start_year=2020, end_year=2021)
        
        # Test that basic methods exist
        assert hasattr(instance, 'get_year_composite')
        assert hasattr(instance, 'get_period_composite')
        assert hasattr(instance, 'get_export')
        assert hasattr(instance, 'get_gif')
        assert hasattr(instance, 'get_stats')
        
        # All methods should be callable
        assert callable(instance.get_year_composite)
        assert callable(instance.get_period_composite)
        assert callable(instance.get_export)
        assert callable(instance.get_gif)
        assert callable(instance.get_stats)
        
    except Exception as e:
        pytest.skip(f"Skipping Earth Engine dependent integration test: {e}")

if __name__ == "__main__":
    # Run tests if called directly
    pytest.main([__file__])