"""
Basic tests for ndvi2gif package (v0.5.0)

- Verifies public API imports
- Smoke tests for NdviSeasonality defaults
- Non-EE tests for TimeSeriesAnalyzer.analyze_trend (using synthetic DataFrame)
- Smoke test for S1ARDProcessor constructor (no EE ops)
- Integration test marked with @pytest.mark.ee (runs only when enabled)
"""

import os
import pytest
import pandas as pd
import numpy as np


# ---------------------------------------------------------------------
# Public API imports
# ---------------------------------------------------------------------

def test_public_api_imports():
    """The top-level package should expose the main classes."""
    import ndvi2gif as n
    assert hasattr(n, "NdviSeasonality")
    assert hasattr(n, "S1ARDProcessor")
    assert hasattr(n, "TimeSeriesAnalyzer")


def test_direct_class_imports():
    """Classes should be importable from their modules too."""
    from ndvi2gif.ndvi2gif import NdviSeasonality, scale_ETM, scale_OLI
    from ndvi2gif.s1_ard import S1ARDProcessor
    from ndvi2gif.timeseries import TimeSeriesAnalyzer

    assert NdviSeasonality is not None
    assert S1ARDProcessor is not None
    assert TimeSeriesAnalyzer is not None
    assert callable(scale_ETM) and callable(scale_OLI)


# ---------------------------------------------------------------------
# NdviSeasonality (pure construction; no EE calls)
# ---------------------------------------------------------------------

def test_ndvi_seasonality_defaults():
    """Defaults should be set at construction without needing EE auth."""
    from ndvi2gif.ndvi2gif import NdviSeasonality
    inst = NdviSeasonality()

    assert inst.periods >= 4
    assert inst.start_year <= inst.end_year
    assert inst.sat in {"S2", "Landsat", "MODIS", "S1", "S3"}
    assert inst.key in {"max", "median", "mean", "percentile"}
    assert isinstance(inst.index, str)

    # Period definitions are created
    assert isinstance(inst.period_names, list) and len(inst.period_names) == inst.periods
    assert isinstance(inst.period_dates, list) and len(inst.period_dates) == inst.periods


def test_ndvi_seasonality_custom_params():
    """Custom params should be accepted and reflected."""
    from ndvi2gif.ndvi2gif import NdviSeasonality
    inst = NdviSeasonality(
        periods=12,
        start_year=2018,
        end_year=2022,
        sat="Landsat",
        key="percentile",
        percentile=85,
        index="evi",
    )
    assert inst.periods == 12
    assert inst.start_year == 2018
    assert inst.end_year == 2022
    assert inst.sat == "Landsat"
    assert inst.key == "percentile"
    assert getattr(inst, "percentile", None) in (85, "85", 85.0)
    assert inst.index == "evi"


def test_valid_satellite_options_updated():
    """Accept supported sats including S3; invalid falls back to default (S2)."""
    from ndvi2gif.ndvi2gif import NdviSeasonality

    for sat in ["S2", "Landsat", "MODIS", "S1", "S3"]:
        assert NdviSeasonality(sat=sat).sat == sat

    assert NdviSeasonality(sat="InvalidSat").sat == "S2"


def test_statistic_key_validation_updated():
    """Support modern keys; invalid falls back to max."""
    from ndvi2gif.ndvi2gif import NdviSeasonality

    for key in ["max", "median", "mean", "percentile"]:
        inst = NdviSeasonality(key=key, percentile=90)
        assert inst.key == key

    assert NdviSeasonality(key="invalid_stat").key == "max"


def test_core_indices_available():
    """Check presence of a core subset of indices (avoid over-constraining)."""
    from ndvi2gif.ndvi2gif import NdviSeasonality
    inst = NdviSeasonality()
    expected = {"ndvi", "ndwi", "mndwi", "evi", "savi", "gndvi", "ndmi"}
    available = set(inst.d.keys())
    missing = expected - available
    assert not missing, f"Missing index methods: {missing}"


# ---------------------------------------------------------------------
# S1ARDProcessor (constructor only; no EE calls)
# ---------------------------------------------------------------------

def test_s1_ard_smoke_constructor():
    """S1ARDProcessor should be constructible without EE initialization."""
    from ndvi2gif.s1_ard import S1ARDProcessor

    proc = S1ARDProcessor(
        speckle_filter="REFINED_LEE",
        terrain_correction=True,
        terrain_flattening_model="VOLUME",
        dem="COPERNICUS_30",
    )
    assert proc is not None
    assert hasattr(proc, "speckle_filter")


# ---------------------------------------------------------------------
# TimeSeriesAnalyzer (no EE): analyze_trend on synthetic DataFrame
# ---------------------------------------------------------------------

def test_timeseries_analyzer_analyze_trend_on_dataframe():
    """
    analyze_trend should work with a pre-made DataFrame (no EE required).
    """
    from ndvi2gif.timeseries import TimeSeriesAnalyzer
    from collections import namedtuple

    # Minimal stub to satisfy TSA constructor without EE usage
    DummyProc = namedtuple(
        "DummyProc",
        [
            "roi", "periods", "start_year", "end_year", "sat",
            "index", "key", "period_names", "period_dates",
        ],
    )
    dummy = DummyProc(
        roi=None,
        periods=12,
        start_year=2018,
        end_year=2020,
        sat="S2",
        index="ndvi",
        key="median",
        period_names=[f"p{i+1}" for i in range(12)],
        period_dates=[("-01-01", "-01-31")] * 12,
    )

    tsa = TimeSeriesAnalyzer(dummy)

    # Synthetic, reproducible time series with positive trend
    np.random.seed(0)
    n = 60
    dates = pd.date_range("2020-01-01", periods=n, freq="7D")
    values = np.linspace(0.2, 0.8, n) + np.random.normal(0, 0.03, n)
    df = pd.DataFrame({"date": dates, "value": values})

    res = tsa.analyze_trend(df=df, method="all", alpha=0.05)
    assert "mann_kendall" in res
    assert "linear" in res
    assert "sen_slope" in res
    assert "interpretation" in res
    assert isinstance(res["linear"]["slope"], float)
    assert isinstance(res["mann_kendall"]["p_value"], float)


# ---------------------------------------------------------------------
# Integration test with Earth Engine
#   - Marked with @pytest.mark.ee
#   - conftest.py/pytest.ini can skip it unless enabled
# ---------------------------------------------------------------------

@pytest.mark.ee
def test_integration_basic_workflow():
    """Runs only when EE tests are enabled (see pytest.ini/conftest.py)."""
    try:
        import ee
        from ndvi2gif.ndvi2gif import NdviSeasonality
        ee.Initialize()
    except Exception as e:
        pytest.skip(f"Earth Engine not initialized: {e}")

    inst = NdviSeasonality(periods=4, start_year=2020, end_year=2021)
    assert hasattr(inst, "get_year_composite")
    assert hasattr(inst, "get_period_composite")
    assert hasattr(inst, "get_export")
    assert hasattr(inst, "get_gif")
    assert hasattr(inst, "get_stats")


if __name__ == "__main__":
    pytest.main([__file__])
