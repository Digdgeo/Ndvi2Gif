# Datasets Reference

Comprehensive guide to the 7 satellite and climate reanalysis platforms supported by ndvi2gif v1.0.0.

## Overview

| Platform | Type | Coverage | Resolution | Variables | Best For |
|----------|------|----------|------------|-----------|----------|
| **Sentinel-2** | Optical | 2015-present | 10-30m | 40+ indices | High-res vegetation |
| **Sentinel-3** | Optical | 2016-present | 300m | 40+ indices | Large-scale ocean/land |
| **Landsat** | Optical | 1982-present | 30m | 40+ indices | Long-term change |
| **MODIS** | Optical | 2000-present | 500m | 40+ indices | Daily global coverage |
| **Sentinel-1** | SAR | 2014-present | 10m | 7 indices | All-weather monitoring |
| **ERA5-Land** | Climate | 1950-present | ~11km | 47 variables | Climate analysis |
| **CHIRPS** | Climate | 1981-present | ~5.5km | 1 variable | Precipitation monitoring |

---

## Optical Sensors

### Sentinel-2 (S2)

**European Space Agency multispectral imaging mission**

```python
processor = NdviSeasonality(
    roi=roi,
    sat='S2',
    index='ndvi',
    start_year=2020,
    end_year=2023
)
```

**Technical Specifications:**
- **Satellites**: Sentinel-2A (launched 2015), Sentinel-2B (launched 2017)
- **Revisit Time**: 5 days (combined constellation)
- **Swath Width**: 290 km
- **Spatial Resolution**:
  - 10m: Blue, Green, Red, NIR
  - 20m: Red Edge (3 bands), SWIR (2 bands)
  - 60m: Coastal Aerosol, Water Vapor, Cirrus
- **Temporal Coverage**: June 2015 - present
- **Earth Engine Collection**: `COPERNICUS/S2_SR_HARMONIZED`

**Spectral Bands:**
| Band | Wavelength (nm) | Resolution | Description |
|------|----------------|------------|-------------|
| B1 | 443 | 60m | Coastal Aerosol |
| B2 | 490 | 10m | Blue |
| B3 | 560 | 10m | Green |
| B4 | 665 | 10m | Red |
| B5 | 705 | 20m | Red Edge 1 |
| B6 | 740 | 20m | Red Edge 2 |
| B7 | 783 | 20m | Red Edge 3 |
| B8 | 842 | 10m | NIR |
| B8A | 865 | 20m | NIR Narrow |
| B9 | 940 | 60m | Water Vapor |
| B11 | 1610 | 20m | SWIR 1 |
| B12 | 2190 | 20m | SWIR 2 |

**Best For:**
- High-resolution vegetation monitoring (10m)
- Agricultural field-level analysis
- Red Edge indices for detailed vegetation health
- Water quality assessment (chlorophyll, turbidity)
- Urban green space mapping

**Advantages:**
- Excellent spatial resolution (10m)
- Free and open data
- Frequent revisit (5 days)
- Red Edge bands for advanced vegetation analysis
- Active mission with ongoing data collection

**Limitations:**
- Shorter temporal archive (2015+)
- Cloud contamination in tropical/temperate regions
- Data volume requires cloud processing

**Recommended Statistical Methods:**
- `key='median'`: Robust to clouds (recommended)
- `key='percentile', percentile=75`: Good cloud-free balance
- `key='max'`: Maximum greenness/NDVI

---

### Sentinel-3 OLCI (S3)

**Ocean and Land Color Instrument for large-scale monitoring**

```python
processor = NdviSeasonality(
    roi=roi,
    sat='S3',
    index='ndvi',
    start_year=2020,
    end_year=2023
)
```

**Technical Specifications:**
- **Satellites**: Sentinel-3A (2016), Sentinel-3B (2018)
- **Revisit Time**: <2 days (combined)
- **Swath Width**: 1270 km
- **Spatial Resolution**: 300m
- **Temporal Coverage**: 2016 - present
- **Earth Engine Collection**: `COPERNICUS/S3/OLCI`

**Spectral Bands:**
- 21 bands from 400-1020 nm
- Optimized for ocean color and land applications
- Includes atmospheric correction bands

**Best For:**
- Regional to global-scale monitoring
- Coastal and inland water quality
- Rapid coverage applications
- Complementing Sentinel-2 for large areas

**Advantages:**
- Very high temporal frequency (<2 days)
- Wide swath (1270 km)
- Excellent for large-scale monitoring
- Free and open data

**Limitations:**
- Lower spatial resolution (300m)
- Not suitable for field-level agriculture
- Fewer band options than Sentinel-2

---

### Landsat (4-9)

**USGS/NASA long-term Earth observation program**

```python
processor = NdviSeasonality(
    roi=roi,
    sat='Landsat',
    index='ndvi',
    start_year=1990,  # Long-term analysis
    end_year=2023
)
```

**Technical Specifications:**
- **Satellites**:
  - Landsat 4-5 (TM): 1982-2013
  - Landsat 7 (ETM+): 1999-present (SLC-off after 2003)
  - Landsat 8-9 (OLI): 2013-present
- **Revisit Time**: 16 days (8 days with L8+L9)
- **Swath Width**: 185 km
- **Spatial Resolution**: 30m (15m panchromatic)
- **Temporal Coverage**: 1982 - present (40+ years!)
- **Earth Engine Collection**:
  - `LANDSAT/LC08/C02/T1_L2` (Landsat 8)
  - `LANDSAT/LC09/C02/T1_L2` (Landsat 9)
  - `LANDSAT/LE07/C02/T1_L2` (Landsat 7)
  - `LANDSAT/LT05/C02/T1_L2` (Landsat 5)
  - `LANDSAT/LT04/C02/T1_L2` (Landsat 4)

**Spectral Bands (Landsat 8/9 OLI):**
| Band | Wavelength (nm) | Resolution | Description |
|------|----------------|------------|-------------|
| B1 | 435-451 | 30m | Coastal Aerosol |
| B2 | 452-512 | 30m | Blue |
| B3 | 533-590 | 30m | Green |
| B4 | 636-673 | 30m | Red |
| B5 | 851-879 | 30m | NIR |
| B6 | 1566-1651 | 30m | SWIR 1 |
| B7 | 2107-2294 | 30m | SWIR 2 |

**Best For:**
- **Long-term change detection** (1982-present)
- Historical vegetation trends
- Decadal climate impact studies
- Consistent 30m resolution time series
- Comparison with modern sensors

**Advantages:**
- **Longest continuous record** (40+ years)
- Free and open data
- Well-documented and validated
- Consistent 30m resolution across sensors
- Global coverage

**Limitations:**
- 16-day revisit (less frequent than Sentinel-2)
- Cloud contamination
- Landsat 7 SLC-off gaps (2003+)
- Lower spatial resolution than Sentinel-2

**Recommended Use Cases:**
- Deforestation and land use change (1980s-present)
- Agricultural expansion tracking
- Urban growth analysis
- Glacier retreat monitoring
- Climate change impact assessment

---

### MODIS

**NASA Terra + Aqua daily global coverage**

```python
processor = NdviSeasonality(
    roi=roi,
    sat='MODIS',
    index='ndvi',
    start_year=2010,
    end_year=2023,
    key='mean'  # Daily data - mean works well
)
```

**Technical Specifications:**
- **Satellites**: Terra (2000), Aqua (2002)
- **Revisit Time**: 1-2 days (daily global coverage)
- **Swath Width**: 2330 km
- **Spatial Resolution**: 500m (surface reflectance)
- **Temporal Coverage**: 2000 - present
- **Earth Engine Collection**: `MODIS/061/MCD43A4` (NBAR/BRDF-adjusted)

**Spectral Bands:**
- 7 bands optimized for land applications (500m)
- 36 total bands (various resolutions)

**Best For:**
- Global-scale monitoring
- Daily vegetation dynamics
- Rapid response applications
- Coarse-resolution regional analysis
- Phenology studies (daily resolution)

**Advantages:**
- Daily global coverage
- Long archive (2000+)
- Excellent temporal resolution
- Free and open data
- BRDF-corrected products

**Limitations:**
- Coarse spatial resolution (500m)
- Not suitable for field-level analysis
- Mixed pixels in heterogeneous landscapes

**Recommended Statistical Methods:**
- `key='mean'`: Good for daily data
- `key='max'`: MVC (Maximum Value Composite) for vegetation

---

## Synthetic Aperture Radar (SAR)

### Sentinel-1 (S1)

**ESA C-band SAR for all-weather monitoring**

```python
from ndvi2gif import NdviSeasonality, S1ARDProcessor

# Configure SAR preprocessing
s1_proc = S1ARDProcessor(
    speckle_filter='REFINED_LEE',
    terrain_correction=True,
    dem='COPERNICUS_30'
)

# Use with NdviSeasonality
sar = NdviSeasonality(
    roi=roi,
    sat='S1',
    index='rvi',  # SAR vegetation index
    periods=12,
    start_year=2020,
    end_year=2023
)
```

**Technical Specifications:**
- **Satellites**: Sentinel-1A (2014), Sentinel-1B (2016-2021, decommissioned)
- **Revisit Time**: 6-12 days (S1A only after 2021)
- **Swath Width**: 250 km (IW mode)
- **Spatial Resolution**: 10m (IW mode)
- **Frequency**: C-band (5.405 GHz)
- **Polarizations**: VV, VH (dual-pol)
- **Temporal Coverage**: 2014 - present
- **Earth Engine Collection**: `COPERNICUS/S1_GRD`

**SAR Indices Available:**
| Index | Formula | Description |
|-------|---------|-------------|
| RVI | 4×VH / (VV+VH) | Radar Vegetation Index |
| DPSVI | VV+VH | Dual-Polarization SAR Vegetation Index |
| RFDI | (VV-VH) / (VV+VH) | Radar Forest Degradation Index |
| VSDI | VV-VH | Vegetation Scattering Difference Index |

**Best For:**
- **All-weather monitoring** (clouds, darkness, smoke)
- Flood mapping and water extent
- Wetland monitoring
- Rice paddy detection
- Soil moisture proxy
- Forest structure analysis

**Advantages:**
- **Cloud-independent** (microwave penetrates clouds)
- Day/night acquisition
- Sensitive to vegetation structure and moisture
- Free and open data
- 10m spatial resolution

**Limitations:**
- More complex preprocessing (speckle, terrain effects)
- Less intuitive than optical imagery
- Geometric distortions in mountains
- Requires understanding of SAR backscatter physics

**Preprocessing Features:**
- **Speckle Filtering**: REFINED_LEE (recommended), LEE, GAMMA_MAP
- **Terrain Correction**: Radiometric and geometric correction
- **Orbit Filtering**: ASCENDING/DESCENDING for temporal consistency

**Recommended Statistical Methods:**
- `key='median'`: Reduces speckle (recommended)
- `key='mean'`: Alternative for smoother results

---

## Climate Reanalysis

### ERA5-Land

**ECMWF high-resolution land surface reanalysis**

```python
processor = NdviSeasonality(
    roi=roi,
    sat='ERA5',
    index='temperature_2m_celsius',  # Or any of 47 variables
    periods=12,
    start_year=1980,  # Can go back to 1950!
    end_year=2023,
    key='mean'  # For temperature
)
```

**Technical Specifications:**
- **Provider**: European Centre for Medium-Range Weather Forecasts (ECMWF)
- **Temporal Resolution**: Daily aggregates (from hourly data)
- **Spatial Resolution**: ~11 km (0.1° grid)
- **Temporal Coverage**: **1950 - present** (75+ years!)
- **Global Coverage**: Complete global land surface
- **Earth Engine Collection**: `ECMWF/ERA5_LAND/DAILY_AGGR`
- **Update Frequency**: ~5 days behind real-time

**Variable Categories (47 total):**

**Temperature (8 core + 16 variants = 24 total):**
- `temperature_2m` - Air temperature at 2m (Kelvin)
- `temperature_2m_celsius` - Auto-converted to °C
- `temperature_2m_min` / `temperature_2m_max` - Daily extremes
- `temperature_2m_min_celsius` / `temperature_2m_max_celsius` - Extremes in °C
- `dewpoint_temperature_2m` (+ celsius/min/max variants)
- `skin_temperature` (+ celsius/min/max variants)
- `soil_temperature_level_1` (+ celsius/min/max variants)

**Precipitation & Water Balance (5 core + 6 unit variants = 11 total):**
- `total_precipitation_sum` - Total precipitation (meters)
- `total_precipitation_sum_lm2` - In liters/m² (mm)
- `total_evaporation_sum` (+ _lm2 variant)
- `potential_evaporation_sum` (+ _lm2 variant)
- `runoff_sum` (+ _lm2 variant)
- `surface_runoff_sum` (+ _lm2 variant)
- `snowfall_sum` (+ _lm2 variant)

**Soil Moisture (4 layers):**
- `volumetric_soil_water_layer_1` - 0-7 cm depth
- `volumetric_soil_water_layer_2` - 7-28 cm depth
- `volumetric_soil_water_layer_3` - 28-100 cm depth
- `volumetric_soil_water_layer_4` - 100-289 cm depth

**Radiation (3 variables):**
- `surface_solar_radiation_downwards_sum`
- `surface_net_solar_radiation_sum`
- `surface_latent_heat_flux_sum`

**Wind & Pressure (3 variables):**
- `u_component_of_wind_10m` - East-west wind (m/s)
- `v_component_of_wind_10m` - North-south wind (m/s)
- `surface_pressure` - Atmospheric pressure (Pa)

**Snow (2 variables):**
- `snow_depth_water_equivalent` - Snow water content (m)
- `snowfall_sum` - Daily snowfall (m)

**Best For:**
- **Long-term climate analysis** (1950-2023)
- Temperature trend detection
- Precipitation patterns and drought
- Soil moisture dynamics
- Climate-vegetation interactions
- Heat/cold wave identification
- Water balance studies

**Advantages:**
- **Exceptional temporal depth** (1950+)
- **Gap-free**: No missing data (modeled reanalysis)
- 47 climate variables in one dataset
- Daily resolution
- Celsius conversion built-in
- Global coverage
- Free access

**Limitations:**
- Coarse resolution (~11 km)
- Not observed data (reanalysis/model)
- Some variables have known biases
- Evaporation bands have data issues (swapped values)

**Recommended Statistical Methods:**
- **Temperature**: `key='mean'` (daily average)
- **Temperature Extremes**: `key='max'` or `key='min'`
- **Precipitation**: `key='sum'` (monthly/seasonal totals)
- **Soil Moisture**: `key='mean'`
- **Radiation**: `key='sum'` (energy accumulation)

**Usage Tips:**
```python
# Temperature in Celsius (automatic conversion)
temp = NdviSeasonality(
    roi=roi,
    sat='ERA5',
    index='temperature_2m_celsius',  # ← Kelvin → Celsius
    key='mean'
)

# Precipitation in mm (L/m²)
precip = NdviSeasonality(
    roi=roi,
    sat='ERA5',
    index='total_precipitation_sum_lm2',  # ← meters → L/m²
    key='sum'  # Sum for totals
)

# Daily temperature extremes
temp_max = NdviSeasonality(
    roi=roi,
    sat='ERA5',
    index='temperature_2m_max_celsius',  # Daily max
    key='max'  # Maximum of daily maxima
)
```

---

### CHIRPS

**Climate Hazards Group InfraRed Precipitation with Station data**

```python
processor = NdviSeasonality(
    roi=roi,
    sat='CHIRPS',
    index='precipitation',
    periods=12,
    start_year=2010,
    end_year=2023,
    key='sum'  # Monthly totals
)
```

**Technical Specifications:**
- **Provider**: UC Santa Barbara Climate Hazards Center
- **Temporal Resolution**: Daily
- **Spatial Resolution**: ~5.5 km (0.05°)
- **Temporal Coverage**: **1981 - present** (40+ years)
- **Coverage**: 50°S to 50°N (quasi-global tropics/subtropics)
- **Earth Engine Collection**: `UCSB-CHG/CHIRPS/DAILY`
- **Update Frequency**: ~3 weeks behind real-time

**Variable:**
- `precipitation` - Daily precipitation (mm/day)

**Data Sources:**
- Satellite infrared cold cloud duration (IR-CCD)
- In-situ station observations (blended)
- TRMM precipitation estimates (calibration)

**Best For:**
- **High-resolution precipitation monitoring**
- Drought early warning systems
- Agricultural rainfall assessment
- Hydrological modeling inputs
- Tropical/subtropical precipitation analysis
- Comparison with ERA5 precipitation

**Advantages:**
- **Higher resolution than ERA5** (~5.5 km vs ~11 km)
- **Station-calibrated** (better accuracy than satellite-only)
- Specifically designed for drought monitoring
- Long archive (1981+)
- Quasi-global tropical coverage
- Free access

**Limitations:**
- **Single variable** (precipitation only)
- Limited to 50°S - 50°N (no polar regions)
- ~3-week latency
- No other climate variables

**Recommended Statistical Methods:**
- `key='sum'`: Monthly/seasonal precipitation totals (recommended)
- `key='mean'`: Average daily rainfall

**Comparison: CHIRPS vs ERA5 Precipitation**

| Aspect | CHIRPS | ERA5 |
|--------|--------|------|
| Resolution | ~5.5 km | ~11 km |
| Method | Satellite + stations | Reanalysis model |
| Coverage | 50°S-50°N | Global |
| Accuracy | High (station blend) | Moderate (model) |
| Variables | Precipitation only | 47 variables |
| Latency | ~3 weeks | ~5 days |

**Usage Example:**
```python
# CHIRPS for detailed precipitation
precip_chirps = NdviSeasonality(
    roi=roi,
    sat='CHIRPS',
    index='precipitation',
    periods=12,
    start_year=2015,
    end_year=2023,
    key='sum'  # Monthly totals
)

# Compare with ERA5
precip_era5 = NdviSeasonality(
    roi=roi,
    sat='ERA5',
    index='total_precipitation_sum_lm2',
    periods=12,
    start_year=2015,
    end_year=2023,
    key='sum'
)
```

---

## Dataset Selection Guide

### By Application

**Agricultural Monitoring:**
- **Field-level (< 50ha)**: Sentinel-2 (10m)
- **Regional**: Landsat (30m) or MODIS (500m)
- **All-weather**: Sentinel-1 SAR
- **Climate context**: ERA5 (temperature, precipitation, soil moisture)

**Forest/Vegetation Change Detection:**
- **High-detail**: Sentinel-2 (10m, 2015+)
- **Long-term trends**: Landsat (30m, 1982+)
- **Rapid monitoring**: MODIS (daily)
- **Forest structure**: Sentinel-1 SAR

**Water Quality & Wetlands:**
- **High-resolution**: Sentinel-2 (chlorophyll, turbidity)
- **Large water bodies**: Sentinel-3 OLCI
- **Water extent**: Sentinel-1 SAR (all-weather)
- **Hydrological drivers**: ERA5 (precipitation), CHIRPS

**Climate & Drought Analysis:**
- **Long-term trends**: ERA5-Land (1950+)
- **High-res precipitation**: CHIRPS (1981+)
- **Vegetation response**: Sentinel-2 or Landsat NDVI
- **Soil moisture**: ERA5 volumetric soil water

**Phenology & Seasonality:**
- **Daily dynamics**: MODIS (500m)
- **Detailed phenology**: Sentinel-2 (10m, 5-day)
- **Multi-decadal**: Landsat (30m, 40+ years)
- **Climate drivers**: ERA5 temperature + precipitation

### By Temporal Coverage

| Need | Dataset(s) |
|------|-----------|
| **1950-1980** | ERA5-Land only |
| **1981-2000** | ERA5-Land + CHIRPS + Landsat 4-5 |
| **2000-2015** | All except S2/S3 |
| **2015-present** | All 7 platforms |

### By Resolution

| Resolution | Dataset(s) | Best Use |
|------------|-----------|----------|
| **10m** | Sentinel-1, Sentinel-2 | Field-level precision |
| **30m** | Landsat | Regional monitoring |
| **300m** | Sentinel-3 | Ocean/large-scale land |
| **500m** | MODIS | Global/continental |
| **5.5 km** | CHIRPS | Regional precipitation |
| **11 km** | ERA5-Land | Climate analysis |

---

## Earth Engine Collection IDs

For direct Earth Engine access:

```python
import ee

# Optical
s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
s3 = ee.ImageCollection('COPERNICUS/S3/OLCI')
l8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
l9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2')
l7 = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2')
l5 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
modis = ee.ImageCollection('MODIS/061/MCD43A4')

# SAR
s1 = ee.ImageCollection('COPERNICUS/S1_GRD')

# Climate
era5 = ee.ImageCollection('ECMWF/ERA5_LAND/DAILY_AGGR')
chirps = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')
```

---

## References

- **Sentinel-2**: [ESA Sentinel-2](https://sentinels.copernicus.eu/web/sentinel/missions/sentinel-2)
- **Sentinel-3**: [ESA Sentinel-3](https://sentinels.copernicus.eu/web/sentinel/missions/sentinel-3)
- **Landsat**: [USGS Landsat Missions](https://www.usgs.gov/landsat-missions)
- **MODIS**: [NASA MODIS](https://modis.gsfc.nasa.gov/)
- **Sentinel-1**: [ESA Sentinel-1](https://sentinels.copernicus.eu/web/sentinel/missions/sentinel-1)
- **ERA5-Land**: Muñoz-Sabater, J. (2019). [DOI:10.24381/cds.e2161bac](https://doi.org/10.24381/cds.e2161bac)
- **CHIRPS**: Funk et al. (2015). Scientific Data. [DOI:10.1038/sdata.2015.66](https://doi.org/10.1038/sdata.2015.66)

---

## See Also

- [Indices Reference](indices.md) - Complete list of 88 variables
- [API Reference](api.md) - Detailed class documentation
- [Climate Tutorial](../tutorials/climate_analysis.md) - ERA5 & CHIRPS examples
