# Indices & Variables Reference

Complete catalog of the variables supported by ndvi2gif, organized by category and satellite platform.

```{note}
Many of the reference DOIs listed for the spectral indices below were sourced from the
[**Awesome Spectral Indices**](https://github.com/awesome-spectral-indices/awesome-spectral-indices)
catalogue (Montero *et al.*, *A standardized catalogue of spectral indices to advance the use of
remote sensing in Earth system research*, **Scientific Data**, 2023). We gratefully acknowledge this
open community resource.
```

## Overview by Category

| Category | Count | Platforms | Description |
|----------|-------|-----------|-------------|
| **Optical Vegetation** | 40 | S2, S3, Landsat, MODIS | NDVI, EVI, SAVI, and derivatives |
| **SAR Vegetation** | 7 | Sentinel-1 | Radar-based vegetation indices |
| **ERA5 Climate** | 40 | ERA5-Land | Temperature, precipitation, soil, radiation |
| **CHIRPS** | 1 | CHIRPS | High-resolution precipitation |
| **Total** | **88** | 7 platforms | Complete variable library |

---

## Optical Vegetation Indices (40)

Compatible with: **Sentinel-2**, **Sentinel-3**, **Landsat 4-9**, **MODIS**

### Core Vegetation Indices (7)

#### NDVI - Normalized Difference Vegetation Index
```python
index='ndvi'
```
**Formula**: `(NIR - Red) / (NIR + Red)`

**Range**: -1 to +1

**Description**: Most widely used vegetation index. Measures greenness and photosynthetic activity.

**Interpretation**:
- < 0: Water, clouds, snow
- 0 - 0.2: Bare soil, rock
- 0.2 - 0.5: Sparse vegetation, grassland
- 0.5 - 0.8: Dense vegetation, crops
- > 0.8: Very dense vegetation, forests

**Reference**: [Rouse et al. (1974)](https://ntrs.nasa.gov/citations/19740022614)

---

#### EVI - Enhanced Vegetation Index
```python
index='evi'
```
**Formula**: `2.5 * (NIR - Red) / (NIR + 6*Red - 7.5*Blue + 1)`

**Range**: -1 to +1

**Description**: Improved vegetation index with atmospheric and soil corrections. Less sensitive to atmospheric conditions than NDVI.

**Advantages**: Better in high biomass areas, reduced soil brightness influence

**Reference**: [Huete et al. (2002)](https://doi.org/10.1016/S0034-4257(96)00112-5)

---

#### SAVI - Soil Adjusted Vegetation Index
```python
index='savi'
```
**Formula**: `((NIR - Red) / (NIR + Red + L)) * (1 + L)` where L=0.5

**Range**: -1 to +1

**Description**: Minimizes soil brightness influence in areas with sparse vegetation.

**Use Cases**: Arid/semi-arid regions, early crop growth

**Reference**: [Huete (1988)](https://doi.org/10.1016/0034-4257(88)90106-X)

---

#### GNDVI - Green Normalized Difference Vegetation Index
```python
index='gndvi'
```
**Formula**: `(NIR - Green) / (NIR + Green)`

**Range**: -1 to +1

**Description**: More sensitive to chlorophyll content than NDVI. Better for later growth stages.

**Use Cases**: Nitrogen assessment, crop health monitoring

**Reference**: [Original reference](https://doi.org/10.1016/S0034-4257(96)00072-7)

---

#### NDMI - Normalized Difference Moisture Index
```python
index='ndmi'
```
**Formula**: `(NIR - SWIR1) / (NIR + SWIR1)`

**Range**: -1 to +1

**Description**: Sensitive to vegetation water content. Useful for drought monitoring.

**Use Cases**: Irrigation management, drought stress detection, fire risk

**Reference**: [Original reference](https://doi.org/10.1016/S0034-4257(01)00318-2)

---

#### EVI2 - Two-Band Enhanced Vegetation Index
```python
index='evi2'
```
**Formula**: `2.5 * (NIR - Red) / (NIR + 2.4*Red + 1)`

**Description**: EVI without blue band requirement. Compatible with more sensors.

**Advantages**: Works with Landsat TM/ETM+ without blue band

**Reference**: [Original reference](https://doi.org/10.1016/j.rse.2008.06.006)

---

#### ARVI - Atmospherically Resistant Vegetation Index
```python
index='arvi'
```
**Formula**: `(NIR - (2*Red - Blue)) / (NIR + (2*Red - Blue))`

**Description**: Reduces atmospheric effects (aerosols, haze) using blue band.

**Use Cases**: Areas with high atmospheric contamination

**Reference**: [Kaufman & Tanré (1992)](https://doi.org/10.1109/36.134076)

---

### Water Indices (5)

#### NDWI - Normalized Difference Water Index
```python
index='ndwi'
```
**Formula**: `(Green - NIR) / (Green + NIR)`

**Range**: -1 to +1

**Description**: Delineates open water features. Opposite of NDVI pattern.

**Interpretation**:
- > 0.3: Water bodies
- 0 - 0.3: Mixed pixels
- < 0: Land

**Reference**: [McFeeters (1996)](https://doi.org/10.1080/01431169608948714)

---

#### MNDWI - Modified Normalized Difference Water Index
```python
index='mndwi'
```
**Formula**: `(Green - SWIR1) / (Green + SWIR1)`

**Description**: Improved water index that suppresses built-up areas better than NDWI.

**Advantages**: Better separation of water from urban areas

**Reference**: [Xu (2006)](https://doi.org/10.1080/01431160600589179)

---

#### NDMI (Water Content)
```python
index='ndmi'
```
*See Moisture section above - dual use for vegetation water content and water body detection*

---

#### LSWI - Land Surface Water Index
```python
index='lswi'
```
**Formula**: `(NIR - SWIR1) / (NIR + SWIR1)`

**Description**: Detects surface water and soil moisture. Similar to NDMI.

**Reference**: [Original reference](https://doi.org/10.1016/j.rse.2003.11.008)

---

#### AWEInsh - Automated Water Extraction Index (no shadow)
```python
index='aweinsh'
```
**Formula**: `4 * (Green - SWIR1) - (0.25*NIR + 2.75*SWIR2)`

**Description**: Eliminates shadow pixels that contaminate water extraction.

**Advantages**: Robust in urban and mountainous areas

**Reference**: [Feyisa et al. (2014)](https://doi.org/10.1016/j.rse.2013.08.029)

---

### Advanced Vegetation Indices (10)

#### OSAVI - Optimized Soil Adjusted Vegetation Index
```python
index='osavi'
```
**Formula**: `(NIR - Red) / (NIR + Red + 0.16)`

**Description**: Optimized L parameter for SAVI. Works across wider vegetation range.

**Reference**: [Original reference](https://doi.org/10.1016/0034-4257(95)00186-7)

---

#### MSAVI - Modified Soil Adjusted Vegetation Index
```python
index='msavi'
```
**Formula**: `(2*NIR + 1 - sqrt((2*NIR+1)^2 - 8*(NIR-Red))) / 2`

**Description**: Self-adjusting L parameter based on NIR-Red relationship.

**Reference**: [Original reference](https://doi.org/10.1016/0034-4257(94)90134-1)

---

#### VARI - Visible Atmospherically Resistant Index
```python
index='vari'
```
**Formula**: `(Green - Red) / (Green + Red - Blue)`

**Description**: Uses only visible bands. Good for RGB sensors without NIR.

**Reference**: [Original reference](https://doi.org/10.1016/S0034-4257(01)00289-9)

---

#### GCI - Green Chlorophyll Index
```python
index='gci'
```
**Formula**: `(NIR / Green) - 1`

**Description**: Sensitive to leaf chlorophyll content variations.

**Use Cases**: Nitrogen status, crop health assessment

**Reference**: [Original reference](https://doi.org/10.1078/0176-1617-00887)

---

#### SIPI - Structure Insensitive Pigment Index
```python
index='sipi'
```
**Formula**: `(NIR - Blue) / (NIR - Red)`

**Description**: Ratio of carotenoids to chlorophyll. Stress indicator.

**Reference**: [Original reference](https://eurekamag.com/research/009/395/009395053.php)

---

#### NBR - Normalized Burn Ratio
```python
index='nbr'
```
**Formula**: `(NIR - SWIR2) / (NIR + SWIR2)`

**Description**: Highlights burned areas and fire severity.

**Use Cases**: Post-fire assessment, burn scar mapping

**Reference**: [Original reference](https://doi.org/10.3133/ofr0211)

---

#### BAIS2 - Burned Area Index for Sentinel-2
```python
index='bais2'
```
**Formula**: Complex formula using Red, NIR, SWIR1, SWIR2

**Description**: Optimized for Sentinel-2 burn severity mapping.

**Platform**: Sentinel-2 only (requires specific bands)

**Reference**: [Original reference](https://doi.org/10.3390/ecrs-2-05177)

---

#### BSI - Bare Soil Index
```python
index='bsi'
```
**Formula**: `((SWIR1 + Red) - (NIR + Blue)) / ((SWIR1 + Red) + (NIR + Blue))`

**Description**: Identifies bare soil and unvegetated areas.

**Reference**: [Original reference](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.465.8749&rep=rep1&type=pdf)

---

#### NDSI - Normalized Difference Snow Index
```python
index='ndsi'
```
**Formula**: `(Green - SWIR1) / (Green + SWIR1)`

**Description**: Identifies snow cover vs clouds.

**Interpretation**: > 0.4 typically indicates snow

**Reference**: [Original reference](https://doi.org/10.1109/IGARSS.1994.399618)

---

#### CIG - Chlorophyll Index Green
```python
index='cig'
```
**Formula**: `(NIR / Green) - 1`

**Description**: Chlorophyll content using green band. Similar to GCI.

**Reference**: [Original reference](https://doi.org/10.1078/0176-1617-00887)

---

### Sentinel-2 Exclusive Indices (Red Edge) (12)

These indices require Sentinel-2's red edge bands (B5, B6, B7).

#### NDRE - Normalized Difference Red Edge
```python
index='ndre'  # Sentinel-2 only
```
**Formula**: `(NIR - RedEdge1) / (NIR + RedEdge1)`

**Description**: Uses red edge for detailed vegetation assessment.

**Advantages**: More sensitive to chlorophyll variations than NDVI

---

#### MCARI - Modified Chlorophyll Absorption Ratio Index
```python
index='mcari'  # Sentinel-2 only
```
**Formula**: `((RedEdge1 - Red) - 0.2*(RedEdge1 - Green)) * (RedEdge1 / Red)`

**Description**: Chlorophyll content with reduced soil/atmospheric effects.

**Reference**: [Original reference](http://dx.doi.org/10.1016/S0034-4257(00)00113-9)

---

#### CIre - Chlorophyll Index Red Edge
```python
index='cire'  # Sentinel-2 only
```
**Formula**: `(NIR / RedEdge1) - 1`

**Description**: Chlorophyll-sensitive index using red edge.

**Reference**: [Original reference](https://doi.org/10.1078/0176-1617-00887)

---

#### IRECI - Inverted Red Edge Chlorophyll Index
```python
index='ireci'  # Sentinel-2 only
```
**Formula**: `(RedEdge3 - Red) / (RedEdge2 / RedEdge1)`

**Description**: ESA's chlorophyll index for Sentinel-2.

**Reference**: [Original reference](https://doi.org/10.1016/j.isprsjprs.2013.04.007)

---

#### RECI - Red Edge Chlorophyll Index
```python
index='reci'  # Sentinel-2 only
```
**Formula**: `(NIR / RedEdge1) - 1`

**Description**: Alternative formulation for chlorophyll content.

**Reference**: [Original reference](https://doi.org/10.1078/0176-1617-00887)

---

#### S2REP - Sentinel-2 Red Edge Position
```python
index='s2rep'  # Sentinel-2 only
```
**Formula**: `705 + 35 * ((((NIR + Red) / 2) - RedEdge1) / (RedEdge2 - RedEdge1))`

**Description**: Estimates position of red edge inflection point (nm).

**Use Cases**: Nitrogen status, crop health

**Reference**: [Original reference](https://doi.org/10.1016/j.isprsjprs.2013.04.007)

---

#### NDVI_re1, NDVI_re2, NDVI_re3
```python
index='ndvi_re1'  # (NIR - RedEdge1) / (NIR + RedEdge1)
index='ndvi_re2'  # (NIR - RedEdge2) / (NIR + RedEdge2)
index='ndvi_re3'  # (NIR - RedEdge3) / (NIR + RedEdge3)
```
**Description**: NDVI variants using each of the 3 red edge bands.

---

#### MTCI - MERIS Terrestrial Chlorophyll Index
```python
index='mtci'  # Sentinel-2 compatible
```
**Formula**: `(RedEdge2 - RedEdge1) / (RedEdge1 - Red)`

**Description**: Chlorophyll concentration in dense canopies.

**Reference**: [Original reference](https://doi.org/10.1080/0143116042000274015)

---

#### MSAVI_re - Modified SAVI with Red Edge
```python
index='msavi_re'  # Sentinel-2 only
```
**Description**: MSAVI using red edge band instead of red.

---

#### PSRI - Plant Senescence Reflectance Index
```python
index='psri'  # Sentinel-2 only
```
**Formula**: `(Red - Green) / RedEdge2`

**Description**: Detects plant senescence and maturity stages.

**Reference**: [Original reference](https://doi.org/10.1034/j.1399-3054.1999.106119.x)

---

### Water Quality Indices (6)

#### NDTI - Normalized Difference Turbidity Index
```python
index='ndti'
```
**Formula**: `(Red - Green) / (Red + Green)`

**Description**: Estimates water turbidity and suspended sediment.

**Reference**: [Original reference](https://doi.org/10.1016/j.rse.2006.07.012)

---

#### Chlorophyll-a (Simple Ratio)
```python
index='chlorophyll'
```
**Formula**: `(Blue / Green)`

**Description**: Estimates chlorophyll-a concentration in water.

**Use Cases**: Lake/reservoir monitoring, eutrophication

---

#### KIVU - Water Clarity Index
```python
index='kivu'
```
**Description**: Specific index for lake water clarity assessment.

---

#### NDCI - Normalized Difference Chlorophyll Index
```python
index='ndci'
```
**Formula**: `(RedEdge1 - Red) / (RedEdge1 + Red)`

**Description**: Chlorophyll concentration in turbid waters.

**Platform**: Sentinel-2 (requires red edge)

**Reference**: [Original reference](https://doi.org/10.1016/j.rse.2011.10.016)

---

#### Turbidity
```python
index='turbidity'
```
**Formula**: Empirical formula using Red band

**Description**: Direct turbidity estimation from red reflectance.

---

#### Cyanobacteria
```python
index='cyanobacteria'
```
**Formula**: Uses Red and NIR bands

**Description**: Detects harmful algal blooms (cyanobacteria).

---

## SAR Indices (7)

Compatible with: **Sentinel-1** (VV + VH polarizations)

### RVI - Radar Vegetation Index
```python
sat='S1', index='rvi'
```
**Formula**: `(4 * VH) / (VV + VH)`

**Range**: 0 to ~4

**Description**: Most common SAR vegetation index. Sensitive to vegetation structure.

**Interpretation**:
- < 0.5: Bare soil, water
- 0.5 - 1.5: Sparse vegetation
- > 1.5: Dense vegetation

**Reference**: [Original reference](https://doi.org/10.2134/agronj1968.00021962006000060016x)

---

### DPSVI - Dual Polarization SAR Vegetation Index
```python
sat='S1', index='dpsvi'
```
**Formula**: `VV + VH`

**Description**: Simple combination sensitive to total vegetation biomass.

---

### RFDI - Radar Forest Degradation Index
```python
sat='S1', index='rfdi'
```
**Formula**: `(VV - VH) / (VV + VH)`

**Description**: Detects forest degradation and deforestation.

**Reference**: [Original reference](https://doi.org/10.5194/bg-9-179-2012)

---

### VSDI - Vegetation Scattering Difference Index
```python
sat='S1', index='vsdi'
```
**Formula**: `VV - VH`

**Description**: Difference between polarizations. Vegetation structure indicator.

---

### VV - Vertical-Vertical Polarization
```python
sat='S1', index='vv'
```
**Description**: Single VV polarization backscatter. Sensitive to surface roughness.

---

### VH - Vertical-Horizontal Polarization
```python
sat='S1', index='vh'
```
**Description**: Cross-polarization. Sensitive to volume scattering (vegetation canopy).

---

### VV/VH Ratio
```python
sat='S1', index='vv_vh_ratio'
```
**Formula**: `VV / VH`

**Description**: Ratio between polarizations. Discriminates land cover types.

---

## ERA5-Land Climate Variables (40)

Compatible with: **ERA5-Land** (ECMWF reanalysis, 1950-present)

### Temperature Variables (18)

#### Core Temperature (Kelvin)
```python
sat='ERA5', index='temperature_2m'
```
**Units**: Kelvin

**Description**: Air temperature at 2 meters height (daily mean).

**Coverage**: 1950-present, ~11km resolution

---

#### Core Temperature (Celsius) - Auto-converted
```python
sat='ERA5', index='temperature_2m_celsius'
```
**Units**: °C (automatically converted from Kelvin)

**Description**: Same as above but in Celsius. Recommended for user-friendly output.

---

#### Daily Temperature Minimum
```python
sat='ERA5', index='temperature_2m_min'  # Kelvin
sat='ERA5', index='temperature_2m_min_celsius'  # °C
```
**Description**: Daily minimum temperature. Use `key='min'` for coldest day of period.

---

#### Daily Temperature Maximum
```python
sat='ERA5', index='temperature_2m_max'  # Kelvin
sat='ERA5', index='temperature_2m_max_celsius'  # °C
```
**Description**: Daily maximum temperature. Use `key='max'` for hottest day of period.

---

#### Dewpoint Temperature
```python
sat='ERA5', index='dewpoint_temperature_2m'  # Kelvin
sat='ERA5', index='dewpoint_temperature_2m_celsius'  # °C
sat='ERA5', index='dewpoint_temperature_2m_min'  # Daily min
sat='ERA5', index='dewpoint_temperature_2m_max'  # Daily max
```
**Description**: Dewpoint temperature at 2m. Humidity indicator.

---

#### Skin Temperature
```python
sat='ERA5', index='skin_temperature'  # Kelvin
sat='ERA5', index='skin_temperature_celsius'  # °C
sat='ERA5', index='skin_temperature_min'  # Daily min
sat='ERA5', index='skin_temperature_max'  # Daily max
```
**Description**: Land surface temperature.

---

#### Soil Temperature (Layer 1: 0-7cm)
```python
sat='ERA5', index='soil_temperature_level_1'  # Kelvin
sat='ERA5', index='soil_temperature_level_1_celsius'  # °C
sat='ERA5', index='soil_temperature_level_1_min'  # Daily min
sat='ERA5', index='soil_temperature_level_1_max'  # Daily max
```
**Description**: Soil temperature at 0-7cm depth.

---

### Precipitation & Water Balance (10)

#### Total Precipitation
```python
sat='ERA5', index='total_precipitation_sum'  # meters
sat='ERA5', index='total_precipitation_sum_lm2'  # L/m² (mm)
```
**Units**: meters or L/m² (1 L/m² = 1 mm)

**Description**: Daily accumulated precipitation (rain + snow).

**Recommended**: Use `key='sum'` for monthly/seasonal totals.

---

#### Total Evaporation
```python
sat='ERA5', index='total_evaporation_sum'  # meters
sat='ERA5', index='total_evaporation_sum_lm2'  # L/m²
```
**Description**: Daily evaporation from surface (negative values = upward flux).

**Note**: Negative in original data; absolute value recommended.

---

#### Potential Evaporation
```python
sat='ERA5', index='potential_evaporation_sum'  # meters
sat='ERA5', index='potential_evaporation_sum_lm2'  # L/m²
```
**Description**: Potential evapotranspiration (atmospheric demand).

---

#### Total Runoff
```python
sat='ERA5', index='runoff_sum'  # meters
sat='ERA5', index='runoff_sum_lm2'  # L/m²
```
**Description**: Total runoff (surface + subsurface).

---

#### Surface Runoff
```python
sat='ERA5', index='surface_runoff_sum'  # meters
sat='ERA5', index='surface_runoff_sum_lm2'  # L/m²
```
**Description**: Surface runoff only (excludes subsurface drainage).

---

#### Snowfall
```python
sat='ERA5', index='snowfall_sum'  # meters water equivalent
sat='ERA5', index='snowfall_sum_lm2'  # L/m²
```
**Description**: Daily snowfall in water equivalent.

---

### Soil Moisture (4 layers)

```python
sat='ERA5', index='volumetric_soil_water_layer_1'  # 0-7 cm
sat='ERA5', index='volumetric_soil_water_layer_2'  # 7-28 cm
sat='ERA5', index='volumetric_soil_water_layer_3'  # 28-100 cm
sat='ERA5', index='volumetric_soil_water_layer_4'  # 100-289 cm
```
**Units**: m³/m³ (volumetric water content)

**Range**: 0-1

**Description**: Soil moisture at different depths.

**Use Cases**:
- Layer 1: Surface moisture, infiltration
- Layer 3: Root zone (crops)
- Layer 4: Deep groundwater interaction

---

### Radiation (3)

#### Surface Solar Radiation Downwards
```python
sat='ERA5', index='surface_solar_radiation_downwards_sum'
```
**Units**: J/m²

**Description**: Solar energy reaching surface. Use `key='sum'`.

---

#### Surface Net Solar Radiation
```python
sat='ERA5', index='surface_net_solar_radiation_sum'
```
**Units**: J/m²

**Description**: Net solar radiation (incoming - outgoing).

---

#### Surface Latent Heat Flux
```python
sat='ERA5', index='surface_latent_heat_flux_sum'
```
**Units**: J/m²

**Description**: Energy used for evaporation.

---

### Wind (2)

```python
sat='ERA5', index='u_component_of_wind_10m'  # East-West
sat='ERA5', index='v_component_of_wind_10m'  # North-South
```
**Units**: m/s

**Description**: Wind velocity components at 10m height.

**Calculate wind speed**: `sqrt(u² + v²)`

---

### Pressure (1)

```python
sat='ERA5', index='surface_pressure'
```
**Units**: Pa

**Description**: Atmospheric pressure at surface.

---

### Snow (2)

#### Snow Depth (Water Equivalent)
```python
sat='ERA5', index='snow_depth_water_equivalent'
```
**Units**: meters of water

**Description**: Snow pack water content.

---

#### Snowfall (Daily)
```python
sat='ERA5', index='snowfall_sum'
```
**Units**: meters of water

**Description**: Daily snowfall accumulation.

---

## CHIRPS Precipitation (1)

Compatible with: **CHIRPS** (1981-present, ~5.5km, 50°S-50°N)

### Precipitation
```python
sat='CHIRPS', index='precipitation'
```
**Units**: mm/day

**Description**: Satellite + station blended precipitation.

**Recommended**: Use `key='sum'` for monthly/seasonal totals.

**Advantages**: Higher resolution (~5.5km) and station-calibrated vs ERA5 precipitation.

---

## Usage Examples

### Optical Vegetation
```python
from ndvi2gif import NdviSeasonality

# Basic NDVI
ndvi = NdviSeasonality(roi=roi, sat='S2', index='ndvi')

# Advanced Sentinel-2 chlorophyll
chlorophyll = NdviSeasonality(roi=roi, sat='S2', index='ndre')

# Water bodies
water = NdviSeasonality(roi=roi, sat='S2', index='mndwi')
```

### SAR Vegetation
```python
# All-weather vegetation monitoring
sar_veg = NdviSeasonality(
    roi=roi,
    sat='S1',
    index='rvi',
    key='median'  # Reduces speckle
)
```

### Climate Variables
```python
# Temperature analysis (Celsius)
temp = NdviSeasonality(
    roi=roi,
    sat='ERA5',
    index='temperature_2m_celsius',
    start_year=1980,
    end_year=2023,
    key='mean'
)

# Precipitation totals
precip = NdviSeasonality(
    roi=roi,
    sat='ERA5',
    index='total_precipitation_sum_lm2',
    key='sum'  # Monthly/seasonal totals
)

# High-resolution precipitation
precip_chirps = NdviSeasonality(
    roi=roi,
    sat='CHIRPS',
    index='precipitation',
    key='sum'
)

# Soil moisture (root zone)
soil = NdviSeasonality(
    roi=roi,
    sat='ERA5',
    index='volumetric_soil_water_layer_3',  # 28-100cm
    key='mean'
)
```

---

## Platform Compatibility Matrix

| Index Type | S2 | S3 | Landsat | MODIS | S1 | ERA5 | CHIRPS |
|------------|----|----|---------|-------|----|-|--------|
| Optical (basic) | ✓ | ✓ | ✓ | ✓ | ✗ | ✗ | ✗ |
| Red Edge | ✓ | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ |
| SAR | ✗ | ✗ | ✗ | ✗ | ✓ | ✗ | ✗ |
| Climate | ✗ | ✗ | ✗ | ✗ | ✗ | ✓ | ✓ |

---

## References

**Vegetation Indices:**
- Rouse et al. (1974) - NDVI
- Huete et al. (2002) - EVI
- McFeeters (1996) - NDWI
- Xu (2006) - MNDWI
- Feyisa et al. (2014) - AWEInsh

**Climate Datasets:**
- Muñoz-Sabater, J. (2019) - ERA5-Land [DOI:10.24381/cds.e2161bac]
- Funk et al. (2015) - CHIRPS [DOI:10.1038/sdata.2015.66]

---

## See Also

- [Datasets Reference](datasets.md) - Detailed platform documentation
- [API Reference](api.md) - Complete class documentation
- [Tutorials](../tutorials/basic_ndvi.md) - Step-by-step index usage
