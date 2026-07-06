# Time Series and Trend Analysis

ndvi2gif provides two complementary classes for temporal analysis:

- **`TimeSeriesAnalyzer`** — extracts time series for a point or region, runs statistical trend tests, and generates multi-panel analysis dashboards with phenological metrics.
- **`SpatialTrendAnalyzer`** — computes per-pixel trend maps across the entire ROI using GEE's server-side linear regression.

Both classes take a configured `NdviSeasonality` instance and inherit all its settings (ROI, satellite, date range, index).

---

## Quick start

```python
import ee
import geemap
from ndvi2gif import NdviSeasonality, TimeSeriesAnalyzer, SpatialTrendAnalyzer

ee.Initialize(project='your-project-id')

roi = ee.Geometry.Rectangle([-6.55, 36.85, -6.10, 37.20])

ns = NdviSeasonality(
    roi=roi,
    sat='S2',
    index='ndvi',
    periods=12,          # monthly composites
    start_year=2017,
    end_year=2023,
)

analyzer = TimeSeriesAnalyzer(ns)
df = analyzer.extract_time_series()
print(df.head())
```

---

## `TimeSeriesAnalyzer`

### Extracting a time series

`extract_time_series()` queries each temporal period across all years and returns a `pandas.DataFrame`.

```python
# Extract from ROI centroid (default)
df = analyzer.extract_time_series()

# Extract from a specific point (lon, lat)
df = analyzer.extract_time_series(point=(-5.97, 37.27))

# Extract mean value over the full ROI (polygon reduction)
df = analyzer.extract_time_series(
    point=roi,
    reducer='mean',
    scale=20,
)
```

The returned DataFrame has columns: `date`, `value`, `year`, `period`, `doy` (day of year), `season`, `month`.

> **Performance note:** `extract_time_series()` calls `reduceRegion` for each period of each year — it sends one GEE request per period. For a 6-year monthly analysis that is 72 requests. The results are cached automatically; subsequent calls with the same parameters return the cached DataFrame instantly.

---

### Trend analysis

`analyze_trend()` runs one or more statistical tests on the extracted time series and returns a summary dictionary.

```python
# Mann-Kendall non-parametric test (recommended for environmental data)
trends = analyzer.analyze_trend(method='mann_kendall')
print(trends['interpretation'])
# → "Significant increasing trend (p=0.012)"

# All methods at once
trends = analyzer.analyze_trend(method='all', alpha=0.05)

# Mann-Kendall result
mk = trends['mann_kendall']
print(f"tau={mk['tau']:.3f}, p={mk['p_value']:.4f}, trend={'↑' if mk['trend'] > 0 else '↓'}")

# Linear regression
lr = trends['linear']
print(f"slope={lr['slope']:.4f}/period, R²={lr['r_squared']:.3f}")
print(f"Annual change: {lr['yearly_change']:.4f}")

# Sen's slope (robust estimate)
sen = trends['sen_slope']
print(f"Median slope: {sen['median_slope']:.4f}")
```

Available methods:

| `method=` | Test | Notes |
|---|---|---|
| `'mann_kendall'` | Mann-Kendall non-parametric test | Robust to non-normality, recommended for NDVI |
| `'linear'` | Ordinary least squares regression | Assumes linearity and normality |
| `'sen_slope'` | Theil-Sen estimator | Robust to outliers, median-based |
| `'all'` | All three methods | |

---

### Comprehensive analysis dashboard

`plot_comprehensive_analysis()` generates a 9-panel figure covering the time series, trend, seasonal patterns, distribution, autocorrelation, and data quality:

```python
fig = analyzer.plot_comprehensive_analysis()

# Save to file
fig = analyzer.plot_comprehensive_analysis(save_path='ndvi_analysis.png')

# For a specific location
fig = analyzer.plot_comprehensive_analysis(point=(-5.97, 37.27))
```

Dashboard panels:
1. Full time series with trend line
2. Seasonal pattern boxplot (by season or month)
3. Annual comparison (year-over-year)
4. Trend summary statistics
5. Autocorrelation function
6. Value distribution histogram
7. Phenology summary (SOS, EOS, peak)
8. Data quality (% missing per period)
9. Seasonal statistics table

---

### Phenology metrics

`extract_phenology_metrics()` extracts Start of Season (SOS), End of Season (EOS), Peak of Season (POS), Length of Season (LOS), and rates of green-up and senescence for each year:

```python
# Threshold method (most common)
phenology = analyzer.extract_phenology_metrics(
    method='threshold',
    threshold_percentile=50,    # SOS/EOS when value crosses 50th percentile
    smoothing=True,             # Savitzky-Golay smoothing before extraction
    smoothing_window=7,
    smoothing_order=3,
    min_season_length=60,       # ignore seasons shorter than 60 days
)

# Results per year
for year, metrics in phenology.items():
    print(f"{year}: SOS={metrics['sos_doy']}, EOS={metrics['eos_doy']}, "
          f"POS={metrics['pos_doy']}, LOS={metrics['los_days']}")
```

Available phenology methods:

| `method=` | Description | Best for |
|---|---|---|
| `'threshold'` | SOS/EOS when signal crosses a percentile threshold | Croplands, deciduous forests |
| `'derivative'` | SOS/EOS at maximum rate of change (inflection point) | Multi-modal or gradual signals |
| `'logistic'` | Fits logistic curve to growing/declining season | Clean unimodal signals |

#### Phenology dashboard

```python
fig = analyzer.plot_phenology_analysis(
    method='threshold',
    threshold_percentile=50,
    save_path='phenology_dashboard.png',
)
```

The phenology dashboard includes: time series with SOS/EOS/POS markers, timing trends across years, amplitude analysis, green-up and senescence rates, season length variability, and data quality indicators.

---

### Comparing phenology years

```python
# Compare phenological metrics across all years
comparison = analyzer.compare_phenology_years(phenology)
print(comparison['summary'])  # text interpretation
```

---

## `SpatialTrendAnalyzer`

`SpatialTrendAnalyzer` computes pixel-wise trend maps across the entire ROI — each pixel gets its own slope, intercept, and total magnitude of change. Computation runs server-side in GEE.

```python
spatial = SpatialTrendAnalyzer(ns)

# Linear trend map
trend_map = spatial.calculate_pixel_trends(
    method='linear',
    min_observations=5,   # mask pixels with fewer valid observations
    scale=20,
)

# trend_map is an ee.Image with bands: slope, intercept, magnitude
```

Output bands:

| Band | Description |
|---|---|
| `slope` | Trend slope (index units per time step) |
| `intercept` | Y-intercept of the fitted line |
| `magnitude` | Total change over the full analysis period (`slope × n_years`) |

```python
# Visualise slope map
Map = geemap.Map()
Map.centerObject(roi, zoom=10)

slope_vis = {
    'bands': ['slope'],
    'min': -0.01, 'max': 0.01,
    'palette': ['d73027', 'fee090', 'ffffff', 'e0f3f8', '4575b4'],
}
Map.addLayer(trend_map, slope_vis, 'NDVI trend (slope)')
Map.add_colorbar(slope_vis, label='NDVI change per period  (blue = greening)')
Map
```

```python
# Export trend map to Google Drive
trend_map_export = spatial.calculate_pixel_trends(
    method='linear',
    export=True,
    scale=10,
)
```

Available methods:

| `method=` | Description |
|---|---|
| `'linear'` | Per-pixel linear regression via `ee.Reducer.linearRegression` |
| `'sen'` | Theil-Sen estimator (robust to outliers) |

---

## Combining time series and spatial analysis

A common workflow: first explore the signal at a representative point, then map it spatially.

```python
ns = NdviSeasonality(roi=roi, sat='S2', index='ndvi', periods=4,
                     start_year=2017, end_year=2023)

# Step 1: explore point time series and check for trends
analyzer = TimeSeriesAnalyzer(ns)
df = analyzer.extract_time_series(point=(-5.97, 37.27))
trends = analyzer.analyze_trend(df, method='all')
fig = analyzer.plot_comprehensive_analysis()

# Step 2: map the trend spatially
spatial = SpatialTrendAnalyzer(ns)
trend_map = spatial.calculate_pixel_trends(method='linear', scale=10)

# Step 3: export
ee.batch.Export.image.toDrive(
    image=trend_map,
    description='ndvi_trend_2017_2023',
    folder='ndvi_analysis',
    region=roi,
    scale=10,
    crs='EPSG:25829',
).start()
```

---

## Tips

### Choosing `periods`

More periods give finer phenological resolution but increase the number of GEE requests in `extract_time_series()`. For trend analysis, 4 periods (seasonal) is usually sufficient. For phenology extraction, 12 or 24 periods give better SOS/EOS estimates.

### Choosing the right trend test

- Use **Mann-Kendall** for environmental indices (NDVI, water quality) — it makes no distributional assumptions and handles autocorrelation better than OLS.
- Use **Sen's slope** for the magnitude estimate — it is more robust to outliers than linear regression.
- Use **linear regression** only when normality and linearity are reasonable assumptions.

### Phenology quality

`extract_phenology_metrics()` prints quality warnings for years with insufficient data or smoothing failures. Low-quality years can be excluded by filtering `phenology.keys()` based on the `quality_warnings` output.

---

## API reference

### `TimeSeriesAnalyzer`

| Method | Returns | Description |
|---|---|---|
| `extract_time_series(point, reducer, scale)` | `pd.DataFrame` | Extract time series for a point or region |
| `analyze_trend(df, method, alpha)` | `dict` | Statistical trend tests |
| `extract_phenology_metrics(df, method, ...)` | `dict` | SOS, EOS, POS, LOS per year |
| `plot_comprehensive_analysis(point, ...)` | `Figure` | 9-panel analysis dashboard |
| `plot_phenology_analysis(point, method, ...)` | `Figure` | Phenology dashboard |
| `compare_phenology_years(phenology)` | `dict` | Inter-annual phenology comparison |

### `SpatialTrendAnalyzer`

| Method | Returns | Description |
|---|---|---|
| `calculate_pixel_trends(method, scale, export)` | `ee.Image` | Per-pixel trend map |

---

## References

Mann, H.B. (1945). Nonparametric tests against trend. *Econometrica*, 13(3), 245–259.

Sen, P.K. (1968). Estimates of the regression coefficient based on Kendall's tau. *Journal of the American Statistical Association*, 63(324), 1379–1389.

Theil, H. (1950). A rank-invariant method of linear and polynomial regression analysis. *Proceedings of the Royal Netherlands Academy of Sciences*, 53, 386–392.

White, M.A., Thornton, P.E., Running, S.W. (1997). A continental phenology model for monitoring vegetation responses to interannual climatic variability. *Global Biogeochemical Cycles*, 11(2), 217–234.
