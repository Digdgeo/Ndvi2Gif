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
    for sat, default in [("S2", "ndvi"), ("Landsat", "ndvi"), ("MODIS", "ndvi"),
                         ("S3", "ndvi"), ("S1", "vv")]:
        with contextlib.redirect_stdout(io.StringIO()):
            inst = NdviSeasonality(roi=roi, sat=sat, index=default,
                                   start_year=2021, end_year=2021)
        image = inst.ndvi_col.filterDate("2021-06-01", "2021-09-01").first()
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


if __name__ == "__main__":
    pytest.main([__file__])
