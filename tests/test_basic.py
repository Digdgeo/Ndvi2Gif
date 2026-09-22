"""
Basic tests for ndvi2gif package (v0.7.0)

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


# Statistical reducers accepted by NdviSeasonality(key=...)
VALID_KEYS = {
    "max", "min", "median", "mean", "sum", "percentile",
    "std", "variance", "range", "cv", "count",
}

# Raw reflectance bands selectable through index=...
RAW_BANDS = {"blue", "green", "red", "nir", "swir1", "swir2"}
S2_REDEDGE_BANDS = {"red_edge1", "red_edge2", "red_edge3"}


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
    assert inst.sat in {"S2", "Landsat", "MODIS", "S1", "S3", "ERA5"}
    assert inst.key in VALID_KEYS
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
    """Accept supported sats including S3, ERA5, and CHIRPS; invalid falls back to default (S2)."""
    from ndvi2gif.ndvi2gif import NdviSeasonality

    # Each sensor needs an index it actually supports (S1 has no 'ndvi')
    sat_index = {
        "S2": "ndvi", "Landsat": "ndvi", "MODIS": "ndvi", "S1": "vv",
        "S3": "ndvi", "ERA5": "temperature_2m", "CHIRPS": "precipitation",
        "VIIRS": "avg_rad", "DMSP": "stable_lights",
    }
    for sat, index in sat_index.items():
        assert NdviSeasonality(sat=sat, index=index).sat == sat

    with pytest.raises(ValueError):
        NdviSeasonality(sat="InvalidSat")


def test_statistic_key_validation_updated():
    """Every documented key is accepted; an unknown one raises ValueError."""
    from ndvi2gif.ndvi2gif import NdviSeasonality

    for key in sorted(VALID_KEYS):
        inst = NdviSeasonality(key=key, percentile=90)
        assert inst.key == key

    with pytest.raises(ValueError):
        NdviSeasonality(key="invalid_stat")


def test_core_indices_available():
    """Check presence of a core subset of indices (avoid over-constraining)."""
    from ndvi2gif.ndvi2gif import NdviSeasonality
    inst = NdviSeasonality()
    expected = {"ndvi", "ndwi", "mndwi", "evi", "savi", "gndvi", "ndmi"}
    available = set(inst.d.keys())
    missing = expected - available
    assert not missing, f"Missing index methods: {missing}"


def test_every_index_is_reachable():
    """Every index in the dispatch dict is registered for at least one sensor.

    An index missing from sensor_indices is rejected by the constructor, so its
    method can never run ('cig' was in that state until 1.6.0).
    """
    from ndvi2gif.ndvi2gif import NdviSeasonality
    inst = NdviSeasonality()

    reachable = set().union(*inst.sensor_indices.values())
    assert set(inst.d) - reachable == set()
    assert reachable - set(inst.d) == set()

    for sat in ("S2", "Landsat", "MODIS", "S3"):
        assert NdviSeasonality(sat=sat, index="cig").index == "cig"


def test_raw_bands_available_per_sensor():
    """Raw reflectance bands are selectable on S2/Landsat/MODIS, not on S1/S3."""
    from ndvi2gif.ndvi2gif import NdviSeasonality
    inst = NdviSeasonality(index="swir1")

    assert RAW_BANDS | S2_REDEDGE_BANDS <= set(inst.d.keys())
    assert RAW_BANDS | S2_REDEDGE_BANDS <= inst.sensor_indices["S2"]

    for sat in ("Landsat", "MODIS"):
        assert RAW_BANDS <= inst.sensor_indices[sat]
        # Red edge only exists in Sentinel-2
        assert not (S2_REDEDGE_BANDS & inst.sensor_indices[sat])

    # Sentinel-3 carries TOA radiances and Sentinel-1 backscatter, so neither
    # should expose bands that claim to be surface reflectance
    for sat in ("S1", "S3"):
        assert not (RAW_BANDS & inst.sensor_indices[sat])

    for sat, index in [("Landsat", "red_edge1"), ("S3", "swir1"), ("S1", "red")]:
        with pytest.raises(ValueError):
            NdviSeasonality(sat=sat, index=index)


def test_raw_band_reflectance_scale():
    """Every optical collection is already reflectance, so raw bands need no factor.

    Since 1.6.0 Sentinel-2 and MODIS are rescaled when the collection is built
    (they store integers x 10000), like Landsat in scale_OLI / scale_ETM.
    """
    from ndvi2gif.ndvi2gif import NdviSeasonality

    for sat in ("S2", "MODIS", "Landsat"):
        assert NdviSeasonality(sat=sat, index="red").reflectance_scale == 1.0


def test_era5_variables_available():
    """ERA5-Land climate variables should be available in dispatch dictionary."""
    from ndvi2gif.ndvi2gif import NdviSeasonality
    inst = NdviSeasonality(sat="ERA5", index="temperature_2m")

    # Check that ERA5 is in the sensor mapping
    assert "ERA5" in inst.sensor_indices

    # Check core ERA5 variables are in dispatch dict
    era5_expected = {
        "temperature_2m", "total_precipitation_sum", "total_evaporation_sum",
        "volumetric_soil_water_layer_1", "surface_pressure"
    }
    available = set(inst.d.keys())
    missing = era5_expected - available
    assert not missing, f"Missing ERA5 variable methods: {missing}"

    # Check that ERA5 variables are mapped to the satellite
    era5_vars = inst.sensor_indices["ERA5"]
    assert "temperature_2m" in era5_vars
    assert "total_precipitation_sum" in era5_vars


def test_chirps_precipitation_available():
    """CHIRPS precipitation variable should be available."""
    from ndvi2gif.ndvi2gif import NdviSeasonality
    inst = NdviSeasonality(sat="CHIRPS", index="precipitation")

    # Check that CHIRPS is in the sensor mapping
    assert "CHIRPS" in inst.sensor_indices

    # Check precipitation variable is in dispatch dict
    assert "precipitation" in inst.d

    # Check that CHIRPS variables are mapped to the satellite
    chirps_vars = inst.sensor_indices["CHIRPS"]
    assert "precipitation" in chirps_vars


def test_nighttime_lights_variables_available():
    """VIIRS and DMSP-OLS variables are registered and reach their methods."""
    from ndvi2gif.ndvi2gif import NdviSeasonality

    viirs = NdviSeasonality(sat="VIIRS", index="avg_rad")
    assert viirs.sensor_indices["VIIRS"] == {"avg_rad", "cf_cvg"}

    dmsp = NdviSeasonality(sat="DMSP", index="stable_lights", periods=1)
    assert dmsp.sensor_indices["DMSP"] == {
        "avg_vis", "stable_lights", "avg_lights_x_pct", "cf_cvg"
    }

    # 'cf_cvg' is a band of both datasets and shares a single method
    assert viirs.d["cf_cvg"].__name__ == dmsp.d["cf_cvg"].__name__ == "get_cf_cvg"

    # An optical index is not available on a lights dataset, and the other way round
    with pytest.raises(ValueError):
        NdviSeasonality(sat="VIIRS", index="ndvi")
    with pytest.raises(ValueError):
        NdviSeasonality(sat="S2", index="avg_rad")


# ---------------------------------------------------------------------
# HydroperiodAnalyzer (construction and validation; no EE computation)
# ---------------------------------------------------------------------

def test_hydroperiod_year_bounds():
    """Cycles run from the configured start day to the same day a year later."""
    from ndvi2gif import NdviSeasonality, HydroperiodAnalyzer
    ns = NdviSeasonality(sat="S2", start_year=2022, end_year=2022)

    # Default hydrological year starts on 1 September; the end is exclusive
    assert HydroperiodAnalyzer(ns)._hyd_year_bounds(2022) == (
        "2022-09-01", "2023-09-01"
    )
    custom = HydroperiodAnalyzer(ns, hydrological_year_start=(10, 1))
    assert custom._hyd_year_bounds(2022) == ("2022-10-01", "2023-10-01")


def test_hydroperiod_rejects_invalid_indices():
    """Only water indices the sensor actually provides are accepted."""
    from ndvi2gif import NdviSeasonality, HydroperiodAnalyzer
    ns = NdviSeasonality(sat="S2", start_year=2022, end_year=2022)
    analyzer = HydroperiodAnalyzer(ns)

    for index in HydroperiodAnalyzer.WATER_INDICES:
        analyzer._validate_index(index)

    # A vegetation index is not a water index
    with pytest.raises(ValueError, match="not a supported water index"):
        analyzer.get_water_masks(index="ndvi")

    # Sentinel-1 has none of the optical water indices
    s1 = HydroperiodAnalyzer(NdviSeasonality(sat="S1", index="vv"))
    with pytest.raises(ValueError, match="not available for sensor"):
        s1.get_water_masks(index="mndwi")


# ---------------------------------------------------------------------
# LandCoverClassifier (multi-sensor argument validation; no EE computation)
# ---------------------------------------------------------------------

def test_classifier_multisensor_validation():
    """Resampling must be chosen explicitly when resolutions differ."""
    from ndvi2gif import NdviSeasonality, LandCoverClassifier

    s2 = NdviSeasonality(sat="S2", index="ndvi")
    s1 = NdviSeasonality(sat="S1", index="vh")
    landsat = NdviSeasonality(sat="Landsat", index="ndvi")

    # Same resolution: nothing to choose
    clf = LandCoverClassifier([s2, s1])
    assert clf.multi_sensor and clf.scale == 10 and clf.crs is None

    # A single sensor keeps working as before
    assert LandCoverClassifier(landsat).scale == 30

    with pytest.raises(ValueError, match="different resolutions"):
        LandCoverClassifier([s2, landsat])
    with pytest.raises(ValueError, match="different sensor"):
        LandCoverClassifier([s2, s2])
    with pytest.raises(ValueError, match="resample must be"):
        LandCoverClassifier([s2, s1], resample="max")

    with pytest.raises(ValueError, match="as a dict"):
        clf.create_feature_stack(indices=["ndvi"])
    with pytest.raises(ValueError, match="missing"):
        clf.create_feature_stack(indices={"S2": ["ndvi"]})
    with pytest.raises(ValueError, match="No processor"):
        clf.create_feature_stack(indices={"S2": ["ndvi"], "S1": ["vh"], "MODIS": ["ndvi"]})
    with pytest.raises(ValueError, match="Invalid indices for S1"):
        clf.create_feature_stack(indices={"S2": ["ndvi"], "S1": ["ndvi"]})


def test_classifier_accuracy_report_and_importance_guard():
    """The accuracy report reads what _calculate_accuracy() stores."""
    from ndvi2gif import NdviSeasonality, LandCoverClassifier
    clf = LandCoverClassifier(NdviSeasonality(sat="S2", index="ndvi"))

    # Shapes as returned by ee.ConfusionMatrix: producers N x 1, consumers
    # 1 x N, indexed by class value. Class 0 is not used here
    clf.accuracy_results = {
        "overall_accuracy": 0.8,
        "kappa": 0.7,
        "producers_accuracy": [[0], [0.9], [0.6]],
        "consumers_accuracy": [[0, 0.75, 0.8]],
        "confusion_matrix": [[0, 0, 0], [0, 9, 1], [0, 3, 4]],
    }
    report = clf.get_accuracy_report()
    assert list(report["Class"]) == [1, 2, "Overall"]
    assert list(report["ProducerAccuracy"]) == [0.9, 0.6, 0.8]
    assert list(report["UserAccuracy"]) == [0.75, 0.8, 0.7]

    # No importance without a trained tree-based classifier
    with pytest.raises(ValueError, match="only available"):
        clf.get_feature_importance()
    clf.classifier, clf.algorithm = object(), "svm"
    with pytest.raises(ValueError, match="only available"):
        clf.get_feature_importance()


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

def _require_ee():
    """Return a working ee module, or skip the test.

    A bare ``ee.Initialize()`` only works when the credentials carry a Cloud
    project, so an already initialized session is reused instead of being
    overwritten.
    """
    try:
        import ee
        try:
            ee.Number(1).getInfo()
        except Exception:
            ee.Initialize()
            ee.Number(1).getInfo()
        return ee
    except Exception as e:
        pytest.skip(f"Earth Engine not initialized: {e}")


@pytest.mark.ee
def test_integration_basic_workflow():
    """Runs only when EE tests are enabled (see pytest.ini/conftest.py)."""
    _require_ee()
    from ndvi2gif.ndvi2gif import NdviSeasonality

    inst = NdviSeasonality(periods=4, start_year=2020, end_year=2021)
    assert hasattr(inst, "get_year_composite")
    assert hasattr(inst, "get_period_composite")
    assert hasattr(inst, "get_export")
    assert hasattr(inst, "get_gif")
    assert hasattr(inst, "get_stats")


@pytest.mark.ee
def test_integration_dispersion_reducers():
    """std/variance/range/cv build composites with the usual period band names."""
    ee = _require_ee()
    from ndvi2gif.ndvi2gif import NdviSeasonality

    roi = ee.Geometry.Rectangle([-6.30, 36.95, -6.25, 37.00])
    expected_bands = ["winter", "spring", "summer", "autumn"]
    means = {}

    for key in ["max", "min", "std", "variance", "range", "cv"]:
        inst = NdviSeasonality(
            roi=roi, periods=4, start_year=2021, end_year=2021,
            sat="S2", key=key, index="ndvi",
        )
        composite = inst.get_year_composite().first()
        # Dispersion reducers must not leak the '_stdDev'/'_variance' suffix
        # that ee.ImageCollection.reduce() appends to every band
        assert composite.bandNames().getInfo() == expected_bands
        means[key] = composite.reduceRegion(
            ee.Reducer.mean(), roi, 100, maxPixels=1e9
        ).getInfo()

    for band in expected_bands:
        assert means["range"][band] == pytest.approx(
            means["max"][band] - means["min"][band], abs=1e-6
        )
        assert means["std"][band] > 0
        assert means["variance"][band] > 0


@pytest.mark.ee
def test_integration_raw_bands_and_count():
    """Raw bands come back as reflectance and key='count' as observation counts."""
    ee = _require_ee()
    from ndvi2gif.ndvi2gif import NdviSeasonality

    roi = ee.Geometry.Rectangle([-6.30, 36.95, -6.25, 37.00])

    def median_of(sat, index, key):
        inst = NdviSeasonality(
            roi=roi, periods=1, start_year=2021, end_year=2021,
            sat=sat, index=index, key=key,
        )
        composite = inst.get_year_composite().first()
        assert composite.bandNames().getInfo() == ["p1"]
        scale = 60 if sat == "S2" else 90
        return composite.reduceRegion(
            ee.Reducer.median(), roi, scale, maxPixels=1e9
        ).getInfo()["p1"]

    # Reflectance, not the raw 0-10000 integers Sentinel-2 stores
    s2_red = median_of("S2", "red", "mean")
    assert 0.0 < s2_red < 1.0

    # Landsat is scaled elsewhere, so both sensors must land on the same range
    landsat_red = median_of("Landsat", "red", "mean")
    assert 0.0 < landsat_red < 1.0

    # Dispersion of a band is small compared with its level on a one-year window
    s2_red_std = median_of("S2", "red", "std")
    assert 0.0 < s2_red_std < s2_red

    # count is a number of images, so it must be a positive integer
    n_obs = median_of("S2", "red", "count")
    assert n_obs >= 1
    assert n_obs == pytest.approx(round(n_obs), abs=1e-6)


@pytest.mark.ee
def test_integration_empty_periods_keep_band_order():
    """A period without images is kept as a band and does not shift the rest.

    Sentinel-2 has no scene over this Doñana ROI in February and March 2017,
    while January and April onwards do have data.
    """
    ee = _require_ee()
    from ndvi2gif.ndvi2gif import NdviSeasonality

    roi = ee.Geometry.Rectangle([-6.30, 36.95, -6.25, 37.00])

    def pixel_count(image):
        return image.reduceRegion(
            ee.Reducer.count(), roi, 200, maxPixels=1e9
        ).getInfo()

    inst = NdviSeasonality(
        roi=roi, periods=12, start_year=2017, end_year=2017,
        sat="S2", index="ndvi", key="median",
    )
    composite = inst.get_year_composite().first()
    assert composite.bandNames().getInfo() == inst.period_names

    counts = pixel_count(composite)
    assert counts["january"] > 0
    assert counts["february"] == 0
    assert counts["march"] == 0
    assert counts["april"] > 0

    # The 'april' band holds April, not the next period with data moved up
    april = inst.get_period_composite(2017, 3)
    diff = composite.select("april").subtract(april).abs()
    max_diff = diff.reduceRegion(
        ee.Reducer.max(), roi, 200, maxPixels=1e9
    ).getInfo()
    assert max_diff["april"] == pytest.approx(0, abs=1e-6)

    # With key='count' an empty period is zero observations, not nodata
    inst = NdviSeasonality(
        roi=roi, periods=12, start_year=2017, end_year=2017,
        sat="S2", index="ndvi", key="count",
    )
    composite = inst.get_year_composite().first()
    assert composite.bandNames().getInfo() == inst.period_names
    feb = composite.select("february").reduceRegion(
        ee.Reducer.minMax(), roi, 200, maxPixels=1e9
    ).getInfo()
    assert feb == {"february_min": 0, "february_max": 0}


@pytest.mark.ee
def test_integration_classifier_stack_years_and_empty_periods():
    """The feature stack covers end_year, tracks skipped years, drops empty periods.

    Sentinel-2 has no scene over this ROI in 2014 and none in winter or spring
    2015, so 2014 is skipped by get_year_composite and the 2015 image sits at
    position 0 of the collection.
    """
    ee = _require_ee()
    from ndvi2gif.ndvi2gif import NdviSeasonality
    from ndvi2gif.clasification import LandCoverClassifier

    roi = ee.Geometry.Rectangle([-6.30, 36.95, -6.25, 37.00])
    inst = NdviSeasonality(
        roi=roi, periods=4, start_year=2014, end_year=2016,
        sat="S2", index="ndvi", key="median",
    )
    clf = LandCoverClassifier(inst)
    stack = clf.create_feature_stack(
        indices=["ndvi"], include_statistics=False, normalize=False
    )

    assert stack.bandNames().getInfo() == [
        "ndvi_2015_summer", "ndvi_2015_autumn",
        "ndvi_2016_winter", "ndvi_2016_spring",
        "ndvi_2016_summer", "ndvi_2016_autumn",
    ]

    # Each band holds the year in its name, not the next year with data
    expected = inst.get_period_composite(2016, 2)
    diff = stack.select("ndvi_2016_summer").subtract(expected).abs()
    max_diff = diff.reduceRegion(
        ee.Reducer.max(), roi, 200, maxPixels=1e9
    ).getInfo()
    assert max_diff["ndvi_2016_summer"] == pytest.approx(0, abs=1e-6)


@pytest.mark.ee
def test_integration_pixel_trends_percentile_band_name():
    """Trend maps work when the composite band is not called 'nd'."""
    ee = _require_ee()
    from ndvi2gif.ndvi2gif import NdviSeasonality
    from ndvi2gif.timeseries import SpatialTrendAnalyzer

    roi = ee.Geometry.Rectangle([-6.30, 36.95, -6.25, 37.00])
    inst = NdviSeasonality(
        roi=roi, periods=4, start_year=2019, end_year=2021,
        sat="S2", index="ndvi", key="percentile", percentile=90,
    )
    assert inst.get_period_composite(2019, 0).bandNames().getInfo() == ["nd_p90"]

    trend = SpatialTrendAnalyzer(inst).calculate_pixel_trends(
        method="linear", min_observations=5
    )
    assert trend.bandNames().getInfo() == ["slope", "intercept", "magnitude"]
    stats = trend.reduceRegion(
        ee.Reducer.count(), roi, 200, maxPixels=1e9
    ).getInfo()
    assert stats["slope"] > 0


@pytest.mark.ee
def test_integration_hydroperiod_invariants():
    """Midpoint weights tile the cycle and the hydroperiod bands stay consistent.

    Doñana marshes, Sentinel-2, hydrological year Sep 2022 - Aug 2023.
    """
    ee = _require_ee()
    from ndvi2gif import NdviSeasonality, HydroperiodAnalyzer

    roi = ee.Geometry.Rectangle([-6.40, 36.93, -6.33, 36.98])
    ns = NdviSeasonality(roi=roi, sat="S2", start_year=2022, end_year=2022)
    analyzer = HydroperiodAnalyzer(ns)

    # --- Water masks: one per day, weights covering the whole cycle ---------
    masks = analyzer.get_water_masks(hyd_year=2022)
    props = ee.Dictionary({
        "time": masks.aggregate_array("system:time_start"),
        "weight": masks.aggregate_array("weight"),
        "start": masks.aggregate_array("start_doy"),
        "end": masks.aggregate_array("end_doy"),
    }).getInfo()

    n_dates = len(props["time"])
    assert n_dates > 10
    # Same-day tiles are mosaicked, so no two masks share a day
    days = {t // 86_400_000 for t in props["time"]}
    assert len(days) == n_dates

    # Each scene owns the days up to the midpoints with its neighbours:
    # contiguous spans from 0 to 365 whose weights add up to the year
    assert props["start"][0] == 0
    assert props["end"][-1] == 365
    assert props["start"][1:] == props["end"][:-1]
    assert sum(props["weight"]) == pytest.approx(365)

    # --- Hydroperiod bands --------------------------------------------------
    result = analyzer.compute_hydroperiod(hyd_year=2022)
    assert result.bandNames().getInfo() == [
        "hydroperiod", "valid_days", "normalized",
        "first_flood_doy", "last_flood_doy",
    ]
    assert result.getInfo()["properties"]["index"] == "mndwi"

    checks = ee.Image.cat([
        result.select("hydroperiod").rename("flood"),
        result.select("valid_days").rename("valid"),
        result.select("normalized").rename("norm"),
        result.select("valid_days").subtract(result.select("hydroperiod"))
        .rename("valid_minus_flood"),
        result.select("last_flood_doy").subtract(result.select("first_flood_doy"))
        .rename("last_minus_first"),
    ])
    stats = checks.reduceRegion(
        ee.Reducer.minMax(), roi, 100, maxPixels=1e9
    ).getInfo()

    assert stats["flood_min"] >= 0
    assert stats["valid_max"] <= 365
    assert stats["valid_minus_flood_min"] >= 0
    assert 0 <= stats["norm_min"] and stats["norm_max"] <= 365
    assert stats["last_minus_first_min"] >= 0
    # The marsh floods part of the year and dries out elsewhere
    assert stats["flood_max"] > 30
    assert stats["flood_min"] == 0

    # --- Downloadable mask stack --------------------------------------------
    stack = analyzer.get_water_masks_stack(hyd_year=2022)
    assert stack.bandNames().size().getInfo() == n_dates
    band_type = stack.bandTypes().values().get(0).getInfo()
    assert (band_type["min"], band_type["max"]) == (0, 255)

    histograms = stack.reduceRegion(
        ee.Reducer.frequencyHistogram(), roi, 200, maxPixels=1e9
    ).values().getInfo()
    seen = {int(value) for hist in histograms for value in hist}
    assert seen <= {0, 1, 2, 255}


@pytest.mark.ee
def test_integration_classifier_multisensor_resampling():
    """Sensors are put on one grid: mean aggregation or bilinear interpolation."""
    ee = _require_ee()
    import math
    from ndvi2gif import NdviSeasonality, LandCoverClassifier

    roi = ee.Geometry.Rectangle([-6.30, 36.95, -6.25, 37.00])
    kw = dict(roi=roi, periods=4, start_year=2021, end_year=2021)
    s2 = NdviSeasonality(sat="S2", index="ndvi", **kw)
    landsat = NdviSeasonality(sat="Landsat", index="ndvi", **kw)
    indices = {"S2": ["ndvi"], "Landsat": ["ndvi"]}

    # --- Coarser: S2 averaged onto the 30 m Landsat grid ---------------------
    clf = LandCoverClassifier([s2, landsat], resample="coarser")
    assert clf.scale == 30
    assert clf.crs == "EPSG:32629"  # UTM zone of Doñana

    stack = clf.create_feature_stack(
        indices=indices, include_statistics=False, normalize=False
    )
    assert stack.bandNames().getInfo() == [
        f"{sat}_ndvi_2021_{season}"
        for sat in ("S2", "Landsat")
        for season in ("winter", "spring", "summer", "autumn")
    ]
    band = stack.select("S2_ndvi_2021_summer")
    assert band.projection().crs().getInfo() == clf.crs
    assert band.projection().nominalScale().getInfo() == pytest.approx(30)

    # The 30 m value is the mean of the nine 10 m pixels inside the cell
    point = ee.Geometry.Point([-6.27, 36.975])
    x, y = point.transform(clf.crs, 1).coordinates().getInfo()
    x0, y0 = math.floor(x / 30) * 30, math.floor(y / 30) * 30
    cell = ee.Geometry.Rectangle([x0, y0, x0 + 30, y0 + 30], clf.crs, False)
    native = s2.get_period_composite(2021, 2)
    inside = native.reduceRegion(
        ee.Reducer.mean().combine(ee.Reducer.count(), None, True),
        cell, 10, crs=clf.crs,
    ).getInfo()
    aggregated = band.reduceRegion(ee.Reducer.first(), point, 30).getInfo()
    assert inside["nd_count"] == 9
    assert aggregated["S2_ndvi_2021_summer"] == pytest.approx(inside["nd_mean"], abs=1e-6)

    # --- Finer: Landsat interpolated onto the 10 m grid ----------------------
    clf = LandCoverClassifier([s2, landsat], resample="finer")
    stack = clf.create_feature_stack(
        indices=indices, include_statistics=False, normalize=False
    )
    band = stack.select("Landsat_ndvi_2021_summer")
    assert band.projection().nominalScale().getInfo() == pytest.approx(10)

    # Interpolation smooths values but keeps the regional mean
    native_mean = landsat.get_period_composite(2021, 2).reduceRegion(
        ee.Reducer.mean(), roi, 30, maxPixels=1e9).getInfo()["nd"]
    fine_mean = band.reduceRegion(
        ee.Reducer.mean(), roi, 10, maxPixels=1e9
    ).getInfo()["Landsat_ndvi_2021_summer"]
    assert fine_mean == pytest.approx(native_mean, abs=0.01)


@pytest.mark.ee
def test_integration_every_registered_index_computes():
    """Each index a sensor accepts can be computed on a real image of it.

    Until 1.6.0 Sentinel-3 accepted twelve SWIR indices (OLCI has no SWIR),
    two of its own water-quality indices looked up a band called 'NIR'
    instead of 'Nir', and 'lst' was offered on sensors without thermal
    bands, where it returned an empty image.
    """
    ee = _require_ee()
    import contextlib, io
    from ndvi2gif import NdviSeasonality

    roi = ee.Geometry.Rectangle([-6.30, 36.95, -6.25, 37.00])
    failures = {}
    # DMSP-OLS ended in 2013, so it needs a date range of its own
    for sat, default, window in [("S2", "ndvi", ("2021-06-01", "2021-09-01")),
                                 ("Landsat", "ndvi", ("2021-06-01", "2021-09-01")),
                                 ("MODIS", "ndvi", ("2021-06-01", "2021-09-01")),
                                 ("S3", "ndvi", ("2021-06-01", "2021-09-01")),
                                 ("S1", "vv", ("2021-06-01", "2021-09-01")),
                                 ("VIIRS", "avg_rad", ("2021-06-01", "2021-09-01")),
                                 ("DMSP", "stable_lights", ("2012-01-01", "2013-01-01"))]:
        with contextlib.redirect_stdout(io.StringIO()):
            inst = NdviSeasonality(roi=roi, sat=sat, index=default,
                                   start_year=2021, end_year=2021)
        image = inst.ndvi_col.filterDate(*window).first()
        for index in sorted(inst.sensor_indices[sat]):
            try:
                ee.Image(inst.d[index](image)).bandNames().getInfo()
            except Exception as e:  # noqa: BLE001 - report every failure at once
                failures[f"{sat}:{index}"] = str(e)[:80]

    assert failures == {}

    inst = NdviSeasonality(sat="S3", index="ndvi")
    assert not {"mndwi", "ndmi", "awei", "lst"} & inst.sensor_indices["S3"]
    assert "lst" not in inst.sensor_indices["S2"]
    assert {"lst", "utfvi"} <= inst.sensor_indices["Landsat"]


@pytest.mark.ee
def test_integration_index_formulas_on_known_reflectance():
    """Index formulas match their published definitions on fixed reflectances."""
    ee = _require_ee()
    from ndvi2gif import NdviSeasonality

    refl = {"Blue": 0.05, "Green": 0.08, "Red": 0.06, "Nir": 0.30,
            "Swir1": 0.15, "Swir2": 0.08}
    image = ee.Image.constant(list(refl.values())).rename(list(refl.keys()))
    b, g, r, n, s1, s2 = refl.values()

    expected = {
        "ndvi": (n - r) / (n + r),
        # ndvi2gif uses L = 0.428 by default (Huete suggested 0.5)
        "savi": 1.428 * (n - r) / (n + r + 0.428),
        "evi": 2.5 * (n - r) / (n + 6 * r - 7.5 * b + 1),
        # Feyisa et al. (2014): the SWIR2 term is subtracted, not added
        "aweinsh": 4 * (g - s1) - (0.25 * n + 2.75 * s2),
        "awei": b + 2.5 * g - 1.5 * (n + s1) - 0.25 * s2,
        "mndwi": (g - s1) / (g + s1),
        # Fisher et al. (2016) coefficients, on reflectance in 0-1
        "wi2015": 1.7204 + 171 * g + 3 * r - 70 * n - 45 * s1 - 71 * s2,
        # Wang & Qu (2007) use the SWIR difference, not the sum
        "nmi": (n - (s1 - s2)) / (n + (s1 - s2)),
    }
    inst = NdviSeasonality(sat="S2", index="ndvi")
    point = ee.Geometry.Point([0, 0])
    for index, value in expected.items():
        got = ee.Image(inst.d[index](image)).reduceRegion(
            ee.Reducer.first(), point, 10).values().get(0).getInfo()
        assert got == pytest.approx(value, abs=1e-5), index


@pytest.mark.ee
def test_integration_optical_sensors_share_reflectance_scale():
    """S2, Landsat and MODIS bands are all reflectance, so indices agree.

    Before 1.6.0 Sentinel-2 and MODIS kept their integer x 10000 values and
    every index with a constant (SAVI, EVI, LAI...) came out several times
    too large on them.
    """
    ee = _require_ee()
    import contextlib, io
    from ndvi2gif import NdviSeasonality

    roi = ee.Geometry.Rectangle([-6.30, 36.95, -6.25, 37.00])
    bands = ["Blue", "Green", "Red", "Nir", "Swir1", "Swir2"]
    values = {}
    for sat in ("S2", "Landsat", "MODIS"):
        with contextlib.redirect_stdout(io.StringIO()):
            inst = NdviSeasonality(roi=roi, sat=sat, start_year=2021, end_year=2021)
        summer = inst.ndvi_col.filterDate("2021-06-01", "2021-09-01").select(bands).median()
        stack = ee.Image.cat([
            summer.select("Red").rename("red"),
            ee.Image(inst.d["savi"](summer)).rename("savi"),
            ee.Image(inst.d["evi"](summer)).rename("evi"),
        ])
        values[sat] = stack.reduceRegion(
            ee.Reducer.median(), roi, 60, maxPixels=1e9).getInfo()

    for sat, v in values.items():
        assert 0 < v["red"] < 1, sat
    for index in ("savi", "evi"):
        ref = values["Landsat"][index]
        for sat in ("S2", "MODIS"):
            assert values[sat][index] == pytest.approx(ref, abs=0.05), (sat, index)


@pytest.mark.ee
def test_integration_sar_indices_on_linear_power():
    """SAR indices use linear power although the collection is in dB."""
    ee = _require_ee()
    import math
    from ndvi2gif import NdviSeasonality, S1ARDProcessor

    vv_db, vh_db = -12.0, -19.0
    vv, vh = 10 ** (vv_db / 10), 10 ** (vh_db / 10)
    image = ee.Image.constant([vv_db, vh_db, 35.0]).rename(["VV", "VH", "angle"])
    point = ee.Geometry.Point([0, 0])

    def first(img):
        return ee.Image(img).reduceRegion(
            ee.Reducer.first(), point, 10).values().get(0).getInfo()

    inst = NdviSeasonality(sat="S1", index="vv")
    expected = {
        "rvi": 4 * vh / (vv + vh),
        "vv_vh_ratio": vv / vh,
        "rfdi": (vv - vh) / (vv + vh),
        "dpsvi": (vv ** 2 + vv * vh) / math.sqrt(2),
        "vv": vv_db,  # single polarizations stay in dB
        "vh": vh_db,
    }
    for index, value in expected.items():
        assert first(inst.d[index](image)) == pytest.approx(value, rel=1e-5), index

    # The ARD processor takes dB, works in linear power and can return dB:
    # with nothing else to do, the round trip gives back the input
    ard = S1ARDProcessor(speckle_filter=None, terrain_correction=False,
                         input_format="DB", format="DB")
    out = ard.process_image(image)
    assert first(out.select("VV")) == pytest.approx(vv_db, abs=1e-4)
    assert first(ard.from_db(image).select("VH")) == pytest.approx(vh, rel=1e-6)

    with pytest.raises(ValueError, match="input_format"):
        S1ARDProcessor(input_format="dB")


@pytest.mark.ee
def test_integration_water_indices_split_water_and_land():
    """Every HydroperiodAnalyzer water index is positive on water, negative on land.

    The analyzer thresholds them at 0. Before 1.6.0 wi2015 sat around 1.5
    everywhere and aweinsh had a flipped term.
    """
    ee = _require_ee()
    import contextlib, io
    from ndvi2gif import NdviSeasonality, HydroperiodAnalyzer

    # Guadalquivir mouth: sea, river and land
    roi = ee.Geometry.Rectangle([-6.45, 36.75, -6.25, 36.95])
    worldcover = ee.ImageCollection("ESA/WorldCover/v200").first()
    with contextlib.redirect_stdout(io.StringIO()):
        inst = NdviSeasonality(roi=roi, sat="S2", start_year=2021, end_year=2021)
    summer = inst.ndvi_col.filterDate("2021-06-01", "2021-09-01").median()

    for index in HydroperiodAnalyzer.WATER_INDICES:
        value = ee.Image(inst.d[index](summer)).rename("v")
        medians = ee.Dictionary({
            name: value.updateMask(mask).reduceRegion(
                ee.Reducer.median(), roi, 60, maxPixels=1e9).get("v")
            for name, mask in (("water", worldcover.eq(80)),
                               ("land", worldcover.neq(80)))
        }).getInfo()
        assert medians["water"] > 0 > medians["land"], (index, medians)


@pytest.mark.ee
def test_integration_s1_ard_dem_has_real_slopes():
    """The terrain-correction DEM keeps its 30 m grid, so slopes are real.

    Mosaicking the Copernicus tiles without their projection, as before
    1.6.0, gave a slope of about 0.07 degrees everywhere and masked part of
    each scene, so terrain correction did nothing.
    """
    ee = _require_ee()
    from ndvi2gif import S1ARDProcessor

    sierra_nevada = ee.Geometry.Rectangle([-3.35, 37.03, -3.28, 37.08])
    slope = ee.Terrain.slope(S1ARDProcessor().dem_ee)
    stats = slope.reduceRegion(
        ee.Reducer.mean().combine(ee.Reducer.max(), None, True),
        sierra_nevada, 90, maxPixels=1e9,
    ).getInfo()
    assert stats["slope_max"] > 30
    assert stats["slope_mean"] > 10

    with pytest.raises(ValueError, match="Unknown DEM"):
        S1ARDProcessor(dem="COPERNICUS_90")


@pytest.mark.ee
def test_integration_training_split_independent_of_stack():
    """The same points and seed give the same train/validation split on any stack."""
    ee = _require_ee()
    import contextlib, io
    from ndvi2gif import NdviSeasonality, LandCoverClassifier

    roi = ee.Geometry.Rectangle([-6.45, 36.95, -6.25, 37.10])
    worldcover = ee.ImageCollection("ESA/WorldCover/v200").first().rename("landcover")
    points = worldcover.stratifiedSample(
        numPoints=30, classBand="landcover", region=roi, scale=30, seed=1,
        geometries=True)

    splits = []
    for indices in (["ndvi"], ["ndvi", "ndwi", "ndmi"]):
        with contextlib.redirect_stdout(io.StringIO()):
            inst = NdviSeasonality(roi=roi, sat="S2", periods=4, start_year=2021,
                                   end_year=2021, key="median")
            clf = LandCoverClassifier(inst)
            clf.create_feature_stack(indices=indices, include_statistics=False,
                                     normalize=False)
            clf.add_training_data(training_points=points, class_property="landcover")
        splits.append(sorted(clf.validation_data.aggregate_array("random").getInfo()))

    assert splits[0] and splits[0] == splits[1]


@pytest.mark.ee
def test_integration_export_model(tmp_path):
    """export_model writes the configuration, the model and the samples."""
    ee = _require_ee()
    import contextlib, io, json
    from ndvi2gif import NdviSeasonality, LandCoverClassifier

    roi = ee.Geometry.Rectangle([-6.45, 36.95, -6.25, 37.10])
    worldcover = ee.ImageCollection("ESA/WorldCover/v200").first().rename("landcover")
    points = worldcover.stratifiedSample(
        numPoints=30, classBand="landcover", region=roi, scale=30, seed=1,
        geometries=True)

    with contextlib.redirect_stdout(io.StringIO()):
        inst = NdviSeasonality(roi=roi, sat="S2", periods=4, start_year=2021,
                               end_year=2021, key="median")
        clf = LandCoverClassifier(inst)
        clf.create_feature_stack(indices=["ndvi", "ndwi"], include_statistics=False,
                                 normalize=False)
        clf.add_training_data(training_points=points, class_property="landcover",
                              seed=3)
        clf.classify_supervised(algorithm="random_forest",
                                params={"numberOfTrees": 10})
        paths = clf.export_model(str(tmp_path / "out" / "model"))

    meta = json.load(open(paths["model"]))
    assert meta["processors"][0]["sat"] == "S2"
    assert meta["features"] == [f"{i}_2021_{s}" for i in ("ndvi", "ndwi")
                                for s in ("winter", "spring", "summer", "autumn")]
    assert meta["training"] == {"class_property": "landcover",
                                "train_fraction": 0.7, "seed": 3}
    assert meta["classifier"]["parameters"]["numberOfTrees"] == 10
    assert len(meta["classifier"]["explain"]["trees"]) == 10
    assert meta["accuracy"]["overall_accuracy"] > 0

    samples = pd.read_csv(paths["samples"])
    assert list(samples.columns[:4]) == ["lon", "lat", "landcover", "split"]
    assert set(samples["split"]) == {"train", "validation"}
    assert set(meta["features"]) <= set(samples.columns)
    assert samples["lon"].between(-6.45, -6.25).all()

    # The CSV is enough to fit an equivalent model outside Earth Engine
    from sklearn.ensemble import RandomForestClassifier
    train = samples[samples["split"] == "train"]
    rf = RandomForestClassifier(n_estimators=10, random_state=0)
    rf.fit(train[meta["features"]], train["landcover"])


def test_period_date_ranges_are_contiguous():
    """A period ends where the next begins, so no day falls between two.

    filterDate excludes its end date, and period_dates stores the last day of
    each period, so filtering by it dropped that day from every composite:
    January ran to the 30th, and February never saw the 29th of a leap year.
    """
    from ndvi2gif.ndvi2gif import NdviSeasonality

    for periods in (4, 12, 24, 365):
        inst = NdviSeasonality(periods=periods)
        ranges = [inst._period_date_range(2016, i) for i in range(periods)]

        assert ranges[0][0] == "2016-01-01"
        assert ranges[-1][1] == "2017-01-01"  # the year is covered to the end
        for (_, end), (start, _) in zip(ranges, ranges[1:]):
            assert end == start
        for start, end in ranges:
            assert start < end  # never an empty range, not even at periods=365

    monthly = NdviSeasonality(periods=12)
    assert monthly._period_date_range(2016, 1) == ("2016-02-01", "2016-03-01")


@pytest.mark.ee
def test_integration_period_covers_every_day():
    """Every day of the year lands in exactly one period.

    CHIRPS is daily, so counting valid observations per period counts days.
    2016 is a leap year: February must have 29.
    """
    ee = _require_ee()
    from ndvi2gif.ndvi2gif import NdviSeasonality

    roi = ee.Geometry.Rectangle([-6.10, 37.25, -5.85, 37.45])
    inst = NdviSeasonality(roi=roi, sat="CHIRPS", index="precipitation",
                           periods=12, key="count",
                           start_year=2016, end_year=2016)
    composite = ee.Image(inst.get_year_composite().first())
    days = composite.reduceRegion(ee.Reducer.max(), roi, 5500,
                                  maxPixels=1e9).getInfo()

    assert days["january"] == 31
    assert days["february"] == 29
    assert days["december"] == 31
    assert sum(days.values()) == 366


@pytest.mark.ee
def test_integration_get_stats_one_row_per_zone():
    """Zonal statistics return one row per feature, attributes kept.

    A path to a shapefile always did. An ee.FeatureCollection passed directly
    fell through to `geom.geometry()`, which dissolves the polygons into one
    zone and drops their attributes: three reservoirs came back as a single
    unnamed row.
    """
    ee = _require_ee()
    from ndvi2gif.ndvi2gif import NdviSeasonality

    zones = ee.FeatureCollection([
        ee.Feature(ee.Geometry.Rectangle([-6.30, 36.95, -6.27, 37.00]), {"zone": "a"}),
        ee.Feature(ee.Geometry.Rectangle([-6.26, 36.95, -6.23, 37.00]), {"zone": "b"}),
        ee.Feature(ee.Geometry.Rectangle([-6.22, 36.95, -6.19, 37.00]), {"zone": "c"}),
    ])

    inst = NdviSeasonality(roi=zones.geometry().bounds(), sat="S2", index="ndvi",
                           key="median", periods=4,
                           start_year=2021, end_year=2021)
    image = ee.Image(inst.get_year_composite().first())

    gdf = inst.get_stats(image=image, geom=zones, stat="MEAN", scale=100)
    assert len(gdf) == 3
    assert sorted(gdf["zone"]) == ["a", "b", "c"]
    assert "summer" in gdf.columns

    # A bare geometry is a single zone
    single = inst.get_stats(image=image, geom=zones.geometry(), stat="MEAN", scale=200)
    assert len(single) == 1


@pytest.mark.ee
def test_integration_peak_period():
    """get_peak_period() maps the period of the maximum, skipping empty ones.

    Same Doñana ROI as the empty-periods test: Sentinel-2 has no scene there
    in February and March 2017, so neither month can be the peak of any pixel.
    """
    ee = _require_ee()
    from ndvi2gif.ndvi2gif import NdviSeasonality

    roi = ee.Geometry.Rectangle([-6.30, 36.95, -6.25, 37.00])
    inst = NdviSeasonality(
        roi=roi, periods=12, start_year=2017, end_year=2017,
        sat="S2", index="ndvi", key="median",
    )

    peak = inst.get_peak_period(return_peak_value=True)
    assert peak.bandNames().getInfo() == ["peak_period", "peak_value"]

    stats = peak.select("peak_period").reduceRegion(
        ee.Reducer.minMax(), roi, 100, maxPixels=1e9
    ).getInfo()
    assert 1 <= stats["peak_period_min"] <= stats["peak_period_max"] <= 12

    # February (2) and March (3) have no data at all, so they can never win.
    # Compared with eq(), not with a frequency histogram: the composite has no
    # fixed projection, and a histogram at a scale other than the native one
    # weights resampled pixels and invents fractional counts for values that
    # no pixel actually holds
    empty_periods = peak.select("peak_period").eq(2).Or(
        peak.select("peak_period").eq(3))
    assert empty_periods.selfMask().reduceRegion(
        ee.Reducer.count(), roi, 10, maxPixels=1e9
    ).getInfo()["peak_period"] == 0

    # The peak value is the maximum across the periods. Both sides are built
    # from the same composite: comparing against a second evaluation of
    # get_year_composite() would differ by the resampling of the borders
    composite = inst.get_year_composite().first()
    direct = inst._peak_from_composite(composite, return_peak_value=True)
    diff = direct.select("peak_value").subtract(
        composite.reduce(ee.Reducer.max())).abs()
    assert diff.reduceRegion(
        ee.Reducer.max(), roi, 10, maxPixels=1e9
    ).getInfo()["peak_value"] == pytest.approx(0, abs=1e-6)

    # mask_below drops the pixels that never reach the threshold; with a
    # threshold above every NDVI value nothing survives
    empty = inst.get_peak_period(mask_below=2)
    assert empty.reduceRegion(
        ee.Reducer.count(), roi, 100, maxPixels=1e9
    ).getInfo()["peak_period"] == 0

    # One image per year, tagged with its year
    yearly = NdviSeasonality(
        roi=roi, periods=4, start_year=2019, end_year=2021,
        sat="S2", index="ndvi", key="median",
    ).get_peak_period(per_year=True)
    assert yearly.aggregate_array("year").getInfo() == [2019, 2020, 2021]

    with pytest.raises(ValueError):
        inst.get_peak_period(across_years="first")


@pytest.mark.ee
def test_integration_peak_period_ties_and_composite():
    """The tie count never drops below 1, and composite= is honoured.

    A maximum is reached by at least the period holding it, so 'ties' starts
    at 1. Integer values are where real ties happen: with key='count' several
    periods routinely share the same number of observations.
    """
    ee = _require_ee()
    from ndvi2gif.ndvi2gif import NdviSeasonality

    roi = ee.Geometry.Rectangle([-6.02, 37.35, -5.93, 37.42])
    inst = NdviSeasonality(roi=roi, sat="VIIRS_DAILY", index="quality_flag",
                           key="count", periods=12,
                           start_year=2017, end_year=2017)

    peak = inst.get_peak_period(return_ties=True)
    bounds = peak.select("ties").reduceRegion(
        ee.Reducer.minMax(), roi, 500, maxPixels=1e9).getInfo()
    assert bounds["ties_min"] >= 1
    assert bounds["ties_max"] > 1      # integer counts do tie

    # composite=: the peak of an image built elsewhere, masked or not
    composite = inst.get_year_composite().first()
    from_composite = inst.get_peak_period(composite=composite)
    diff = from_composite.subtract(inst.get_peak_period(across_years="median")).abs()
    assert diff.reduceRegion(
        ee.Reducer.max(), roi, 500, maxPixels=1e9
    ).getInfo()["peak_period"] == 0    # one year: median of it is itself

    with pytest.raises(ValueError):
        inst.get_peak_period(composite=composite.select(["january", "february"]))


@pytest.mark.ee
def test_integration_peak_statistics_is_circular():
    """The peak month is a circular variable, and the statistics must know it.

    Feeding known peak periods straight in, bypassing the imagery, lets the
    answers be checked against ones worked out by hand. The two that matter
    are the cases an arithmetic mean gets wrong: December and January average
    to June instead of late December, and a bimodal pixel averages to a month
    in which nothing ever happened there.
    """
    ee = _require_ee()
    from ndvi2gif.ndvi2gif import NdviSeasonality

    roi = ee.Geometry.Rectangle([-6.0, 37.0, -5.9, 37.1])
    inst = NdviSeasonality(roi=roi, periods=12, sat="S2", index="ndvi",
                           start_year=2020, end_year=2022)

    def stats(periods):
        peaks = ee.ImageCollection([
            ee.Image.constant(p).rename("peak_period").toInt().clip(roi)
            for p in periods])
        return inst.get_peak_statistics(peaks=peaks, min_years=1).reduceRegion(
            ee.Reducer.first(), roi.centroid(1), 1000).getInfo()

    # every year at the same period: the mean is that period, R is 1
    for period in (1, 6, 12):
        out = stats([period] * 3)
        assert out["circular_mean"] == pytest.approx(period, abs=1e-6)
        assert out["concentration"] == pytest.approx(1.0, abs=1e-6)
        assert out["n_years"] == 3

    # the wrap-around: an arithmetic mean would answer 6.5, June
    out = stats([12, 1])
    assert out["circular_mean"] == pytest.approx(12.5, abs=1e-6)
    assert out["concentration"] > 0.9

    out = stats([11, 12, 1, 2])
    assert out["circular_mean"] == pytest.approx(12.5, abs=1e-6)

    # bimodal: the mean lands in April, and R is what says not to trust it
    out = stats([2, 6])
    assert out["circular_mean"] == pytest.approx(4.0, abs=1e-6)
    assert out["concentration"] == pytest.approx(0.5, abs=1e-6)

    # four periods evenly around the circle cancel exactly
    assert stats([1, 4, 7, 10])["concentration"] == pytest.approx(0.0, abs=1e-6)

    # min_years masks the pixels with too few valid years
    peaks = ee.ImageCollection([
        ee.Image.constant(6).rename("peak_period").toInt().clip(roi)])
    masked = inst.get_peak_statistics(peaks=peaks, min_years=2).reduceRegion(
        ee.Reducer.first(), roi.centroid(1), 1000).getInfo()
    assert masked["circular_mean"] is None


@pytest.mark.ee
def test_integration_water_mask_modes_are_ordered():
    """permanent <= dynamic <= maximum, whatever the reservoir does.

    The three detected modes are nested by construction: a pixel that is water
    in every period is water in some period. If that ordering ever breaks, the
    masks are not measuring what they claim. The regression guarded here is
    'permanent', which used to come back empty for any year with a period that
    had no scene, because it compared against the period count instead of the
    periods that actually carried data.
    """
    ee = _require_ee()
    from ndvi2gif.ndvi2gif import NdviSeasonality

    # Doñana marshes: seasonal flooding, so the three modes differ a lot
    roi = ee.Geometry.Rectangle([-6.30, 36.95, -6.25, 37.00])
    inst = NdviSeasonality(roi=roi, sat="S2", index="ndci", key="median",
                           periods=4, start_year=2020, end_year=2020)

    areas = {}
    for mode in ("permanent", "dynamic", "maximum"):
        _, area = inst.get_water_mask(mode=mode, return_area=True)
        areas[mode] = [v for v in area[2020] if v is not None]
        assert areas[mode], f"{mode} returned no area at all"

    assert max(areas["permanent"]) <= max(areas["dynamic"]) + 1e-6
    assert max(areas["dynamic"]) <= max(areas["maximum"]) + 1e-6

    # the static modes repeat one value, the dynamic one does not have to
    assert len(set(round(v, 6) for v in areas["maximum"])) == 1

    # the mask carries one band per period, named like the periods
    masks = inst.get_water_mask(mode="dynamic")
    assert ee.Image(masks.first()).bandNames().getInfo() == inst.period_names

    # and the instance keeps the areas, as period_scene_counts does
    assert set(inst.water_mask_area) == {2020}

    # compositing the configured index still works afterwards: the temporary
    # switch to the water index must leave the instance as it found it
    assert inst.index == "ndci"

    with pytest.raises(ValueError):
        inst.get_water_mask(mode="whatever")
    with pytest.raises(ValueError):
        inst.get_water_mask(water_index="ndvi")


@pytest.mark.ee
def test_integration_dmsp_conversion_is_monotone_and_bounded():
    """The VIIRS to DMSP conversion may never fall, nor leave 0..63.

    This is the regression that nearly shipped. A quadratic fitted to the
    overlap curved back down above about 50 nW/cm2/sr, so the centre of
    Seville at 107 nW came out as digital number 0 — a saturated city read as
    unlit. The saturating model cannot do that, and this pins it down.
    """
    ee = _require_ee()
    from ndvi2gif.ndvi2gif import NdviSeasonality

    roi = ee.Geometry.Rectangle([-6.1, 37.3, -5.9, 37.5])
    inst = NdviSeasonality(roi=roi, sat="VIIRS", index="avg_rad", key="mean",
                           periods=12, start_year=2015, end_year=2015)

    # coefficients of the shape the fit returns, from the Andalusian coast
    coefficients = {"a": 0.3489, "b": 0.6810, "ceiling": 63}

    radiances = [0.0, 0.1, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 1000]
    converted = [
        inst.to_dmsp_like(
            ee.Image.constant(value).rename("avg_rad"), coefficients
        ).reduceRegion(ee.Reducer.first(), roi.centroid(1), 1000
                       ).getInfo()["dmsp_like"]
        for value in radiances
    ]

    assert all(0 <= v <= 63 for v in converted), converted
    assert converted == sorted(converted), converted
    assert converted[0] == pytest.approx(0.0, abs=1e-9)   # no light, no lights
    assert converted[-1] == pytest.approx(63.0, abs=0.1)  # a city saturates

    # the bright end must not collapse: this is what the quadratic got wrong
    assert converted[radiances.index(100)] > converted[radiances.index(10)]


if __name__ == "__main__":
    pytest.main([__file__])
