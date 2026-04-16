# Time Series Analysis with `TimeSeriesAnalyzer`

`TimeSeriesAnalyzer` turns the period-by-period composites produced by `NdviSeasonality` into a full time-series analysis pipeline: point or polygon extraction, trend detection, phenology metrics, and diagnostic dashboards. Everything runs as a thin Python layer on top of GEE — extraction uses `reduceRegion` lazily, and only the resulting table is brought back to the client for statistics and plotting.

For **per-pixel** trend maps over an ROI (rather than a single point or polygon), ndvi2gif also ships `SpatialTrendAnalyzer`, covered at the end of this page.

---

## Quick start

```python
import ee
from ndvi2gif import NdviSeasonality, TimeSeriesAnalyzer

ee.Initialize()

ns = NdviSeasonality(
    roi=roi,
    sat='S2',
    index='ndvi',
    start_year=2018, end_year=2024,
    periods=12,
    key='median',
)

ts = TimeSeriesAnalyzer(ns)

df     = ts.extract_time_series()                # ROI centroid, mean reducer
trend  = ts.analyze_trend(method='mann_kendall')
pheno  = ts.extract_phenology_metrics(method='threshold')

fig = ts.plot_comprehensive_analysis()
```

`TimeSeriesAnalyzer` reuses the `NdviSeasonality` configuration — ROI, sensor, date range, periods, index. You do not need to repeat any of it.

---

## 1. Extracting a time series

`extract_time_series()` walks every period × year combination in the `NdviSeasonality` date range and calls `reduceRegion` once per period. The result is a `pandas.DataFrame` with one row per period:

| Column | Description |
|---|---|
| `date` | Mid-date of the period (`datetime`) |
| `value` | Reduced index value (float) |
| `year` | Year |
| `period` | Period name (`january`, `winter`, `p3`, ...) |
| `doy` | Day of year of the mid-date |
| `season` | Meteorological season |
| `month` | Month number |

### Three ways to specify location

```python
# 1) ROI centroid (default)
df = ts.extract_time_series()

# 2) A specific lon/lat point
df = ts.extract_time_series(point=(-5.5, 37.1))

# 3) An ee.Geometry (Point or Polygon)
df = ts.extract_time_series(point=ee.Geometry.Polygon([...]),
                            reducer='median', scale=20)
```

For points, a small buffer of `scale/2` is added automatically so the extraction is stable against sub-pixel registration issues. For polygons, the `reducer` argument selects the spatial statistic (`'mean'`, `'median'`, `'max'`, `'min'`, `'stdDev'`).

> **Scale matters.** Pass `scale=10` for S2, `30` for Landsat, `500` for MODIS, `20` for S2 Red-Edge bands. Over-resampling (small scale on a coarse sensor) wastes GEE quota; under-resampling (large scale on a fine sensor) mixes the signal with neighbouring land cover.

### Caching

Results are cached per `(point, reducer, scale)` key. Subsequent calls with the same arguments are instant. Pass `use_cache=False` to force re-extraction after changing the underlying `NdviSeasonality` configuration.

---

## 2. Trend analysis

`analyze_trend()` runs one or more statistical trend tests on the time series:

```python
trend = ts.analyze_trend(method='mann_kendall', alpha=0.05)
# or all three at once
trend_all = ts.analyze_trend(method='all')
```

| Method | What it returns | When to use |
|---|---|---|
| `mann_kendall` | `tau`, `p_value`, `trend` direction | **Default** — non-parametric, robust to seasonality and outliers |
| `linear` | `slope`, `intercept`, `r_squared`, CI, `yearly_change` | When you need an absolute rate of change per year |
| `sen_slope` | Median slope, CI | Robust slope estimator, pair with Mann-Kendall |
| `all` | All three + `interpretation` text | When you want a direct comparison |

The returned `interpretation` field is a one-line human-readable summary (e.g. `"Significant increasing trend (p=0.003)"`) intended for reports or quick-look printouts.

> **Mann-Kendall is the right default for vegetation indices.** NDVI, NDWI and friends have strong seasonality and non-Gaussian residuals. The rank-based Mann-Kendall statistic ignores the seasonal cycle's shape and only tests monotonic change, which is what "is this area greening or browning?" actually asks.

---

## 3. Phenology metrics

`extract_phenology_metrics()` computes start-of-season (SOS), end-of-season (EOS), peak-of-season (POS), and derived metrics for each year in the series. Three extraction methods are available:

| Method | Definition of SOS/EOS | Pros | Cons |
|---|---|---|---|
| `threshold` | Crossing of a per-year percentile | Simple, robust | Threshold is arbitrary — calibrate |
| `derivative` | Inflection points of the smoothed curve | Driven by the shape of the curve | Sensitive to noise; needs smoothing |
| `logistic` | Double-logistic fit, derived SOS/EOS | Best fit for unimodal vegetation cycles | Fails on multi-modal cycles |

```python
pheno = ts.extract_phenology_metrics(
    method='threshold',
    threshold_percentile=50,     # mid-amplitude
    smoothing=True,              # Savitzky-Golay
    smoothing_window=7,
    smoothing_order=3,
    min_season_length=60,        # days; filters false short seasons
)
```

Output is a `dict` keyed by year. Each value contains:

- `sos`, `eos`, `pos` — day of year
- `length_of_season` — days
- `amplitude` — value range
- `greenup_rate`, `senescence_rate` — value/day
- `smoothed` — whether smoothing was applied

### When smoothing matters

Savitzky-Golay smoothing is **on by default**. For S2 monthly composites it typically reduces noise enough that `method='derivative'` becomes usable. For raw Landsat or MODIS series with irregular gaps, consider calling `compare_smoothing_impact()` to see how much the metrics move:

```python
impact = ts.compare_smoothing_impact(method='threshold', threshold_percentile=50)
# returns dict with raw vs smoothed SOS/EOS/POS and quantified differences
```

### Comparing years

```python
comparison = ts.compare_phenology_years(reference_year=2020)
# anomalies of every year relative to 2020 (or relative to multi-year mean if None)
```

---

## 4. Diagnostic dashboards

Two ready-made multi-panel figures are available:

### `plot_comprehensive_analysis()`

Nine-panel layout: time series + trend line, seasonal pattern, year-on-year comparison, trend summary, autocorrelation, distribution, data quality, seasonal statistics, and annual boxplots.

```python
fig = ts.plot_comprehensive_analysis(
    point=None,              # ROI centroid if None
    figsize=(22, 14),
    save_path='ndvi_dashboard.png',
)
```

### `plot_phenology_analysis()`

Eight-panel phenology dashboard: smoothed series with SOS/EOS/POS markers, annual SOS/EOS timing, amplitude, greenup/senescence rates, duration, per-year comparison, summary statistics, and a data-quality panel.

```python
fig = ts.plot_phenology_analysis(
    method='threshold',
    threshold_percentile=50,
    save_path='phenology_dashboard.png',
)
```

These are designed to fit a paper figure or a report page without further editing. When you need publication-level control, use `extract_time_series()` + your own Matplotlib code.

---

## 5. Per-pixel spatial trends with `SpatialTrendAnalyzer`

`TimeSeriesAnalyzer` reduces the ROI down to a single series. When you want **a trend map** — per-pixel slope or Mann-Kendall statistic across the ROI — use `SpatialTrendAnalyzer`:

```python
from ndvi2gif import SpatialTrendAnalyzer

st = SpatialTrendAnalyzer(ns)

trend_map = st.calculate_pixel_trends(
    method='linear',          # or 'sen', 'mann_kendall'
    min_observations=5,
    scale=30,
)
# trend_map is an ee.Image with per-pixel slope / tau / p-value bands
```

Visualise it like any other `ee.Image`:

```python
import geemap
Map = geemap.Map()
Map.centerObject(ns.roi, zoom=10)
Map.addLayer(
    trend_map.select('slope'),
    {'min': -0.01, 'max': 0.01,
     'palette': ['#67001f','#d6604d','#f7f7f7','#4393c3','#053061']},
    'NDVI trend slope (yr⁻¹)',
)
Map
```

This runs entirely on GEE servers — no client-side loop over pixels, and it scales to large ROIs that would be unmanageable in pandas.

---

## Tips and caveats

### `extract_time_series()` is lazy-per-period but eager-per-call

Each period triggers one `reduceRegion().getInfo()` round-trip. For long time series (e.g. Landsat 1984–2024 with `periods=12` = 492 calls) this can be slow. Strategies:

- Limit extraction to the range you need.
- Prefer larger `scale` if sub-pixel accuracy is not required.
- Use cached results (`use_cache=True`, default).

For pixel-level maps over long ranges, use `SpatialTrendAnalyzer` — it does the whole thing server-side.

### Minimum data for trends and phenology

- Trends: at least 3 points are required; in practice ≥20 is needed for Mann-Kendall to be meaningful.
- Phenology: ≥4 points per year; `threshold_percentile=50` is the most forgiving default.

The analyzer prints quality warnings for years that fail these thresholds rather than silently dropping them.

### Phenology methods are not interchangeable

SOS/EOS from `threshold`, `derivative`, and `logistic` will not match. Pick one and stick with it across the whole analysis; compare across sites with the same method, not across methods on the same site.

### Smoothing changes the answer

Savitzky-Golay smoothing (window 7, order 3) is a reasonable default for monthly S2 series, but it widens the transition and can shift SOS/EOS by up to a week. Use `compare_smoothing_impact()` to quantify this for your data before publishing numbers.

### Point vs polygon extraction

Points are more sensitive to sub-pixel registration; buffered polygons (≥ 3×3 pixels) give cleaner series for trend detection. For phenology of a single field, polygon + `'mean'` is usually the safer choice.

---

## API reference

### `TimeSeriesAnalyzer(ndvi_seasonality_instance)`

| Method | Returns | Description |
|---|---|---|
| `extract_time_series(point, reducer, scale, use_cache)` | `pd.DataFrame` | Per-period time series at a point or polygon |
| `analyze_trend(df, method, alpha)` | `dict` | Mann-Kendall, linear, and/or Sen's-slope trend tests |
| `extract_phenology_metrics(df, method, threshold_percentile, smoothing, ...)` | `dict[year → metrics]` | SOS/EOS/POS + derived metrics per year |
| `compare_phenology_years(point, reference_year)` | `dict` | Anomalies vs a reference year or multi-year mean |
| `compare_smoothing_impact(point, method, threshold_percentile)` | `dict` | Quantifies impact of smoothing on phenology metrics |
| `plot_comprehensive_analysis(point, figsize, save_path)` | `matplotlib.Figure` | 9-panel time-series dashboard |
| `plot_phenology_analysis(method, threshold_percentile, save_path)` | `matplotlib.Figure` | 8-panel phenology dashboard |

### `SpatialTrendAnalyzer(ndvi_seasonality_instance)`

| Method | Returns | Description |
|---|---|---|
| `calculate_pixel_trends(method, min_observations, scale, export)` | `ee.Image` | Per-pixel trend map (slope / tau / p-value) |

**Supported trend methods:** `linear`, `sen`, `mann_kendall` (and `all` for `TimeSeriesAnalyzer.analyze_trend`)

**Supported phenology methods:** `threshold`, `derivative`, `logistic`

---

## References

Mann, H.B. (1945). Nonparametric tests against trend. *Econometrica*, 13, 245–259.

Kendall, M.G. (1975). *Rank Correlation Methods*. Griffin, London.

Sen, P.K. (1968). Estimates of the regression coefficient based on Kendall's tau. *Journal of the American Statistical Association*, 63(324), 1379–1389.

Savitzky, A., Golay, M.J.E. (1964). Smoothing and differentiation of data by simplified least squares procedures. *Analytical Chemistry*, 36(8), 1627–1639.

Jönsson, P., Eklundh, L. (2002). Seasonality extraction by function fitting to time-series of satellite sensor data. *IEEE Transactions on Geoscience and Remote Sensing*, 40(8), 1824–1832.

Zhang, X., Friedl, M.A., Schaaf, C.B., Strahler, A.H., Hodges, J.C.F., Gao, F., Reed, B.C., Huete, A. (2003). Monitoring vegetation phenology using MODIS. *Remote Sensing of Environment*, 84(3), 471–475.
