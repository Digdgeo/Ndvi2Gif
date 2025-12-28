# Climate Data Analysis with ERA5 & CHIRPS

```{admonition} New in v1.0.0 🌡️
:class: tip
This tutorial showcases the new climate reanalysis capabilities introduced in v1.0.0, including ERA5-Land (47 variables, 1950-present) and CHIRPS precipitation (1981-present).
```

## Overview

Ndvi2Gif v1.0.0 extends beyond vegetation monitoring to comprehensive **climate analysis**. You can now analyze:

- 🌡️ **Temperature** (ERA5): Air, surface, soil temperatures with daily min/max
- 🌧️ **Precipitation** (ERA5 & CHIRPS): Multiple high-quality rainfall datasets
- 💧 **Soil Moisture** (ERA5): 4 depth layers
- ☀️ **Radiation** (ERA5): Solar radiation and heat flux
- 💨 **Wind & Pressure** (ERA5): Wind components and atmospheric pressure
- ❄️ **Snow** (ERA5): Snow depth and snowfall

## Dataset Comparison

| Dataset | Variables | Resolution | Coverage | Best For |
|---------|-----------|------------|----------|----------|
| **ERA5-Land** | 47 climate vars | ~11 km | 1950-present, global | Multi-variable climate analysis |
| **CHIRPS** | Precipitation | ~5.5 km | 1981-present, 50°S-50°N | High-res rainfall, drought monitoring |

## ERA5-Land: Temperature Analysis

### Basic Temperature Analysis

```python
import ee
from ndvi2gif import NdviSeasonality
from ndvi2gif.timeseries import TimeSeriesAnalyzer

ee.Initialize()

# Define ROI (e.g., Doñana National Park, Spain)
roi = ee.Geometry.Point([-6.48, 37.13]).buffer(10000)  # 10km buffer

# Monthly mean temperature (in Celsius)
temp = NdviSeasonality(
    roi=roi,
    sat='ERA5',
    index='temperature_2m_celsius',  # Automatic Kelvin→Celsius conversion
    periods=12,                      # Monthly composites
    start_year=2020,
    end_year=2023,
    key='mean'                       # Average daily temperature
)

# Generate composite
composite = temp.get_year_composite()
print(f"Generated {len(temp.imagelist)} annual composites")
```

### Temperature Extremes (Min/Max)

```python
# Daily maximum temperatures (summer heat analysis)
temp_max = NdviSeasonality(
    roi=roi,
    sat='ERA5',
    index='temperature_2m_max_celsius',  # Daily maximum
    periods=12,
    start_year=2020,
    end_year=2023,
    key='max'  # Maximum of daily maxima (hottest day of month)
)

# Daily minimum temperatures (frost analysis)
temp_min = NdviSeasonality(
    roi=roi,
    sat='ERA5',
    index='temperature_2m_min_celsius',  # Daily minimum
    periods=12,
    start_year=2020,
    end_year=2023,
    key='min'  # Minimum of daily minima (coldest night of month)
)
```

### Time Series Analysis

```python
# Extract temperature time series
analyzer = TimeSeriesAnalyzer(temp)
df = analyzer.extract_time_series(point=(-6.48, 37.13), scale=11000)

print(df.head())
#         date      value  year    period  doy  season  month
# 0 2020-01-16  11.2       2020   january   16  winter      1
# 1 2020-02-14  12.5       2020  february   45  winter      2
# 2 2020-03-16  15.8       2020     march   75  spring      3

# Analyze temperature trends
trends = analyzer.analyze_trend(df=df, method='mann_kendall')
print(f"Temperature trend: {trends['interpretation']}")

# Comprehensive dashboard (shows climate stats, not phenology)
fig = analyzer.plot_comprehensive_analysis()
fig.savefig('temperature_analysis.png', dpi=300)
```

## Precipitation Analysis

### ERA5-Land Precipitation

```python
# Monthly precipitation totals (in L/m² = mm)
precip_era5 = NdviSeasonality(
    roi=roi,
    sat='ERA5',
    index='total_precipitation_sum_lm2',  # Automatic m→L/m² conversion
    periods=12,
    start_year=2020,
    end_year=2023,
    key='sum'  # Sum for monthly totals
)

# Extract and analyze
analyzer = TimeSeriesAnalyzer(precip_era5)
df_precip = analyzer.extract_time_series(scale=11000)

print(f"Annual precipitation 2020-2023:")
for year in range(2020, 2024):
    yearly_precip = df_precip[df_precip['year'] == year]['value'].sum()
    print(f"  {year}: {yearly_precip:.1f} mm")
```

### CHIRPS Precipitation (Higher Resolution)

```python
# CHIRPS: Higher spatial resolution, station-calibrated
precip_chirps = NdviSeasonality(
    roi=roi,
    sat='CHIRPS',
    index='precipitation',  # mm/day
    periods=12,
    start_year=2020,
    end_year=2023,
    key='sum'  # Monthly totals
)

# Compare ERA5 vs CHIRPS
analyzer_chirps = TimeSeriesAnalyzer(precip_chirps)
df_chirps = analyzer_chirps.extract_time_series(scale=5566)

# Correlation analysis
import numpy as np
from scipy.stats import pearsonr

# Align time series
common_dates = set(df_precip['date']).intersection(df_chirps['date'])
era5_vals = df_precip[df_precip['date'].isin(common_dates)]['value'].values
chirps_vals = df_chirps[df_chirps['date'].isin(common_dates)]['value'].values

r, p = pearsonr(era5_vals, chirps_vals)
print(f"ERA5 vs CHIRPS correlation: r={r:.3f}, p={p:.4f}")
```

## Drought Monitoring Example

```python
# Combine precipitation and evapotranspiration for drought index
# Step 1: Get precipitation
precip = NdviSeasonality(
    roi=roi,
    sat='CHIRPS',
    index='precipitation',
    periods=12,
    start_year=2015,
    end_year=2023,
    key='sum'
)

# Step 2: Get potential evapotranspiration
pet = NdviSeasonality(
    roi=roi,
    sat='ERA5',
    index='potential_evaporation_sum_lm2',
    periods=12,
    start_year=2015,
    end_year=2023,
    key='sum'
)

# Step 3: Calculate water balance
analyzer_precip = TimeSeriesAnalyzer(precip)
analyzer_pet = TimeSeriesAnalyzer(pet)

df_precip = analyzer_precip.extract_time_series(scale=5566)
df_pet = analyzer_pet.extract_time_series(scale=11000)

# Simple drought indicator: P - PET
df_precip['water_balance'] = df_precip['value'] - df_pet['value'].abs()
df_precip['drought'] = df_precip['water_balance'] < -50  # Deficit > 50mm

print(f"Drought months (2015-2023): {df_precip['drought'].sum()}")
```

## Soil Moisture Analysis

```python
# Root zone soil moisture (0-100cm)
soil_moisture = NdviSeasonality(
    roi=roi,
    sat='ERA5',
    index='volumetric_soil_water_layer_3',  # 28-100cm depth
    periods=12,
    start_year=2020,
    end_year=2023,
    key='mean'
)

analyzer = TimeSeriesAnalyzer(soil_moisture)
df_soil = analyzer.extract_time_series(scale=11000)

# Identify dry/wet periods
dry_threshold = df_soil['value'].quantile(0.25)
wet_threshold = df_soil['value'].quantile(0.75)

df_soil['condition'] = 'normal'
df_soil.loc[df_soil['value'] < dry_threshold, 'condition'] = 'dry'
df_soil.loc[df_soil['value'] > wet_threshold, 'condition'] = 'wet'

print(df_soil['condition'].value_counts())
```

## Combined Climate-Vegetation Analysis

```python
# Correlate vegetation (NDVI) with climate variables

# 1. Get NDVI time series
ndvi = NdviSeasonality(
    roi=roi,
    sat='S2',
    index='ndvi',
    periods=12,
    start_year=2020,
    end_year=2023,
    key='median'
)

# 2. Get temperature
temp = NdviSeasonality(
    roi=roi,
    sat='ERA5',
    index='temperature_2m_celsius',
    periods=12,
    start_year=2020,
    end_year=2023,
    key='mean'
)

# 3. Get precipitation
precip = NdviSeasonality(
    roi=roi,
    sat='CHIRPS',
    index='precipitation',
    periods=12,
    start_year=2020,
    end_year=2023,
    key='sum'
)

# 4. Extract and merge
from ndvi2gif.timeseries import TimeSeriesAnalyzer
import pandas as pd

df_ndvi = TimeSeriesAnalyzer(ndvi).extract_time_series(scale=20)
df_temp = TimeSeriesAnalyzer(temp).extract_time_series(scale=11000)
df_prec = TimeSeriesAnalyzer(precip).extract_time_series(scale=5566)

# Merge on date
df_combined = df_ndvi[['date', 'value']].rename(columns={'value': 'ndvi'})
df_combined = df_combined.merge(
    df_temp[['date', 'value']].rename(columns={'value': 'temperature'}),
    on='date'
)
df_combined = df_combined.merge(
    df_prec[['date', 'value']].rename(columns={'value': 'precipitation'}),
    on='date'
)

# Correlation matrix
print(df_combined[['ndvi', 'temperature', 'precipitation']].corr())

# Lagged correlations (vegetation response to climate)
df_combined['precip_lag1'] = df_combined['precipitation'].shift(1)
df_combined['temp_lag1'] = df_combined['temperature'].shift(1)

print("\nLagged correlations (1-month lag):")
print(df_combined[['ndvi', 'precip_lag1', 'temp_lag1']].corr())
```

## Export Climate Data

```python
# Export temperature rasters to Google Drive
temp = NdviSeasonality(
    roi=roi,
    sat='ERA5',
    index='temperature_2m_celsius',
    periods=12,
    start_year=2023,
    end_year=2023,
    key='mean'
)

# Generate composites
temp.get_year_composite()

# Export all 12 monthly images
for i, image in enumerate(temp.imagelist[0]):
    task = ee.batch.Export.image.toDrive(
        image=image,
        description=f'temp_2023_month_{i+1:02d}',
        folder='ERA5_Temperature',
        region=roi,
        scale=11000,
        crs='EPSG:4326'
    )
    task.start()
    print(f"Started export: month {i+1}")
```

## Available ERA5 Variables

### Temperature (24 variables)
- `temperature_2m`, `dewpoint_temperature_2m`, `skin_temperature`, `soil_temperature_level_1`
- Add `_min` or `_max` for daily extremes (8 variables)
- Add `_celsius` for Celsius conversion (12 variables)

### Precipitation & Water Balance (11 variables)
- `total_precipitation_sum`, `total_evaporation_sum`, `potential_evaporation_sum`
- `runoff_sum`, `surface_runoff_sum`, `snowfall_sum`
- Add `_lm2` for L/m² units (6 variables)

### Soil Moisture (4 variables)
- `volumetric_soil_water_layer_1` (0-7 cm)
- `volumetric_soil_water_layer_2` (7-28 cm)
- `volumetric_soil_water_layer_3` (28-100 cm)
- `volumetric_soil_water_layer_4` (100-289 cm)

### Radiation (3 variables)
- `surface_solar_radiation_downwards_sum`
- `surface_net_solar_radiation_sum`
- `surface_latent_heat_flux_sum`

### Wind & Pressure (3 variables)
- `u_component_of_wind_10m`, `v_component_of_wind_10m`
- `surface_pressure`

### Snow (2 variables)
- `snow_depth_water_equivalent`, `snowfall_sum`

## Key Differences: Climate vs Vegetation Data

| Aspect | Vegetation (S2, Landsat) | Climate (ERA5, CHIRPS) |
|--------|-------------------------|----------------------|
| **Statistics** | max, median, percentile (cloud-free) | mean, sum, min, max (all valid) |
| **Time Series** | Phenology metrics (SOS, EOS, POS) | Seasonal climate statistics |
| **Temporal Res** | 5-16 days | Daily (aggregated to periods) |
| **Spatial Res** | 10-30m | 5.5-11 km |
| **Coverage** | 2013-present (S2) | 1950-present (ERA5) |

## Tips & Best Practices

```{admonition} Best Practices
:class: tip
1. **Use `sum` for precipitation/flux variables** (total accumulation)
2. **Use `mean` for temperature/pressure** (average conditions)
3. **Use `min`/`max` for extremes** (frost risk, heat stress)
4. **CHIRPS for precipitation** (better station calibration)
5. **ERA5 for multi-variable analysis** (temperature + precip + soil moisture)
6. **Match scales**: ERA5 ~11km, CHIRPS ~5.5km
7. **Climate dashboards** automatically replace phenology panels
```

## References

- **ERA5-Land**: Muñoz Sabater, J. (2019). [DOI:10.24381/cds.e2161bac](https://doi.org/10.24381/cds.e2161bac)
- **CHIRPS**: Funk et al. (2015). Scientific Data. [DOI:10.1038/sdata.2015.66](https://doi.org/10.1038/sdata.2015.66)

---

```{seealso}
- {doc}`time_series` - Advanced time series analysis
- {doc}`statistical_methods` - Statistical reducers explained
- {doc}`../reference/datasets` - Complete dataset documentation
```
