# Indices & Variables Reference

Complete catalog of the variables supported by ndvi2gif, organized by the bands they
need and therefore by the sensors that can compute them.

```{note}
Many of the reference DOIs listed for the spectral indices below were sourced from the
[**Awesome Spectral Indices**](https://github.com/awesome-spectral-indices/awesome-spectral-indices)
catalogue (Montero *et al.*, *A standardized catalogue of spectral indices to advance the use of
remote sensing in Earth system research*, **Scientific Data**, 2023). We gratefully acknowledge this
open community resource.
```

```{tip}
This page is kept in sync with the code, but the authoritative list for a sensor is always
`NdviSeasonality(sat=...).get_available_indices()`. Asking for an index a sensor cannot
compute raises a `ValueError` listing the ones it can.
```

## Overview

| Group | Count | S2 | Landsat | MODIS | S3 | S1 | ERA5 | CHIRPS |
|-------|------:|:--:|:-------:|:-----:|:--:|:--:|:----:|:------:|
| [Visible and NIR indices](#visible-and-nir-indices) | 13 | ✓ | ✓ | ✓ | ✓ | | | |
| [SWIR indices](#swir-indices) | 13 | ✓ | ✓ | ✓ | | | | |
| [Thermal indices](#thermal-indices) | 2 | | ✓ | `lst` only | | | | |
| [Sentinel-2 red edge](#sentinel-2-red-edge-indices) | 9 | ✓ | | | | | | |
| [Sentinel-3 OLCI water quality](#sentinel-3-olci-water-quality-indices) | 10 | | | | ✓ | | | |
| [Raw reflectance bands](#raw-reflectance-bands) | 9 | ✓ | 6 | 6 | | | | |
| [SAR](#sar-indices) | 7 | | | | | ✓ | | |
| [ERA5-Land climate](#era5-land-climate-variables) | 47 | | | | | | ✓ | |
| [CHIRPS](#chirps-precipitation) | 1 | | | | | | | ✓ |
| **Total** | **111** | **44** | **34** | **33** | **23** | **7** | **47** | **1** |

Usage is the same for all of them:

```python
from ndvi2gif import NdviSeasonality

ns = NdviSeasonality(roi=roi, sat='S2', index='evi', key='median',
                     periods=12, start_year=2023, end_year=2023)
```

### Units: surface reflectance on every optical sensor

Formulas below use the standardized band names (`Blue`, `Green`, `Red`, `Nir`, `Swir1`,
`Swir2`, `Red_Edge1-3`) and assume **surface reflectance in 0–1**, which is what every
optical collection holds: Landsat is rescaled by `scale_OLI` / `scale_ETM`, and Sentinel-2
and MODIS — which store reflectance as integers multiplied by 10000 — are rescaled when the
collection is built.

```{warning}
Before **v1.6.0** Sentinel-2 and MODIS were left as integers × 10000. Pure band ratios
(NDVI, NDWI, MNDWI, NDMI...) were not affected, but every index with a constant, a sum or an
inverse was: `savi`, `evi`, `lai`, `avi`, `wi2015`, `cri1`, `cri2`, `fai`, `mcari`, `ireci`,
and the magnitude of `awei`/`aweinsh`. Results computed with those indices on Sentinel-2 or
MODIS with earlier versions should be recomputed. See the [changelog](changelog.md).
```

Sentinel-3 is the exception: the OLCI collection available in Earth Engine holds
**top-of-atmosphere radiances**, not reflectance, so its indices are relative indicators
rather than calibrated quantities.

---

## Visible and NIR indices

**13 indices** · Sentinel-2, Landsat, MODIS, Sentinel-3

They only need visible and near-infrared bands, so they are available on every optical sensor.

| `index=` | Name | Formula | Reference |
|---|---|---|---|
| `ndvi` | Normalized Difference Vegetation Index | `(Nir − Red) / (Nir + Red)` | [Rouse et al. (1974)](https://ntrs.nasa.gov/citations/19740022614) |
| `evi` | Enhanced Vegetation Index | `2.5 · (Nir − Red) / (Nir + 6·Red − 7.5·Blue + 1)` | [Huete et al. (2002)](https://doi.org/10.1016/S0034-4257(02)00096-2) |
| `savi` | Soil Adjusted Vegetation Index | `(1 + L) · (Nir − Red) / (Nir + Red + L)`, L = 0.428 | [Huete (1988)](https://doi.org/10.1016/0034-4257(88)90106-X) |
| `gndvi` | Green NDVI | `(Nir − Green) / (Nir + Green)` | [Gitelson et al. (1996)](https://doi.org/10.1016/S0034-4257(96)00072-7) |
| `wdrvi` | Wide Dynamic Range Vegetation Index | `(0.1·Nir − Red) / (0.1·Nir + Red)` | [Gitelson (2004)](https://doi.org/10.1078/0176-1617-01176) |
| `avi` | Advanced Vegetation Index | `(Nir · (1 − Red) · (Nir − Red))^(1/3)` | Bannari, Asalhi & Teillet (2002), IGARSS |
| `lai` | Leaf Area Index (empirical) | `3.618 · EVI − 0.118` | [Boegh et al. (2002)](https://doi.org/10.1016/S0034-4257(01)00342-X) |
| `cig` | Chlorophyll Index Green | `Nir / Green − 1` | [Gitelson et al. (2003)](https://doi.org/10.1078/0176-1617-00887) |
| `cri1` | Carotenoid Reflectance Index 1 | `1/Blue − 1/Green` | [Gitelson et al. (2002)](https://doi.org/10.1562/0031-8655(2002)075%3C0272:ACCIPL%3E2.0.CO;2) |
| `cri2` | Carotenoid Reflectance Index 2 | `1/Blue − 1/Red` | [Gitelson et al. (2002)](https://doi.org/10.1562/0031-8655(2002)075%3C0272:ACCIPL%3E2.0.CO;2) |
| `pri` | Photochemical Reflectance Index (broadband) | `(Green − Blue) / (Green + Blue)` | [Gamon et al. (1992)](https://doi.org/10.1016/0034-4257(92)90059-S) |
| `ndwi` | Normalized Difference Water Index | `(Green − Nir) / (Green + Nir)` | [McFeeters (1996)](https://doi.org/10.1080/01431169608948714) |
| `vci` | Vegetation Condition Index (simplified) | `100 · (NDVI − 0.1) / (0.8 − 0.1)`, clamped to 0–100 | [Kogan (1995)](https://doi.org/10.1016/0273-1177(95)00079-T) |

Notes on how some of them are implemented:

- **`savi`** uses L = 0.428 rather than the 0.5 Huete suggested as a general value.
- **`lai`** is an empirical regression on EVI calibrated for agricultural crops; treat it as
  a relative indicator elsewhere.
- **`cri2`** uses the red band where the original index uses 700 nm, so that it can be
  computed on sensors without a red edge band.
- **`pri`** is a broadband proxy: the original index contrasts two narrow bands at 531 and
  570 nm that no sensor here has.
- **`vci`** normally scales NDVI between its long-term minimum and maximum *for each pixel*;
  here fixed values (0.1 and 0.8) are used, so it is a rescaled NDVI. For the real VCI, compute
  the per-pixel extremes from a long `NdviSeasonality` series.

---

## SWIR indices

**13 indices** · Sentinel-2, Landsat, MODIS (not Sentinel-3, whose OLCI sensor has no SWIR)

| `index=` | Name | Formula | Reference |
|---|---|---|---|
| `mndwi` | Modified NDWI | `(Green − Swir1) / (Green + Swir1)` | [Xu (2006)](https://doi.org/10.1080/01431160600589179) |
| `awei` | Automated Water Extraction Index (shadow) | `Blue + 2.5·Green − 1.5·(Nir + Swir1) − 0.25·Swir2` | [Feyisa et al. (2014)](https://doi.org/10.1016/j.rse.2013.08.029) |
| `aweinsh` | Automated Water Extraction Index (no shadow) | `4·(Green − Swir1) − (0.25·Nir + 2.75·Swir2)` | [Feyisa et al. (2014)](https://doi.org/10.1016/j.rse.2013.08.029) |
| `wi2015` | Water Index 2015 | `1.7204 + 171·Green + 3·Red − 70·Nir − 45·Swir1 − 71·Swir2` | [Fisher et al. (2016)](https://doi.org/10.1016/j.rse.2015.12.055) |
| `ndmi` | Normalized Difference Moisture Index | `(Nir − Swir1) / (Nir + Swir1)` | Hardisky et al. (1983); [Wilson & Sader (2002)](https://doi.org/10.1016/S0034-4257(01)00318-2) |
| `msi` | Moisture Stress Index | `Swir1 / Nir` | [Hunt & Rock (1989)](https://doi.org/10.1016/0034-4257(89)90046-1) |
| `nmi` | Normalized Multi-band Drought Index | `(Nir − (Swir1 − Swir2)) / (Nir + (Swir1 − Swir2))` | [Wang & Qu (2007)](https://doi.org/10.1029/2007GL031021) |
| `nbr` | Normalized Burn Ratio | `(Nir − Swir2) / (Nir + Swir2)` | Key & Benson (2006), FIREMON, USDA RMRS-GTR-164-CD |
| `nbri` | Normalized Burn Ratio (alias of `nbr`) | `(Nir − Swir2) / (Nir + Swir2)` | as `nbr` |
| `ndbi` | Normalized Difference Built-up Index | `(Swir1 − Nir) / (Swir1 + Nir)` | [Zha et al. (2003)](https://doi.org/10.1080/01431160304987) |
| `ndsi` | Normalized Difference Snow Index | `(Green − Swir1) / (Green + Swir1)` | [Dozier (1989)](https://doi.org/10.1016/0034-4257(89)90101-6) |
| `ndti` | Normalized Difference **Tillage** Index | `(Swir1 − Swir2) / (Swir1 + Swir2)` | Van Deventer et al. (1997), *PE&RS* 63(1) |
| `fai` | Floating Algae Index | `Nir − (Red + (Swir2 − Red) · f)` | [Hu (2009)](https://doi.org/10.1016/j.rse.2009.05.012) |

- **Water indices** (`mndwi`, `ndwi`, `awei`, `aweinsh`, `wi2015`) are positive over water and
  negative over land. `HydroperiodAnalyzer` uses them with a default threshold of 0.
- **`ndti`** is the *tillage* index (crop residue on soil), not the *turbidity* index that
  shares its acronym. For water turbidity see `turbidity` and `spm` on Sentinel-3.
- **`mndwi`** and **`ndsi`** have the same formula; the name only says what it is used for.
- **`fai`** uses SWIR2 as the baseline band, with a wavelength factor `f` per sensor (0.1109
  for Sentinel-2, 0.1359 for Landsat 8/9, 0.1438 for MODIS). Hu (2009) used the 1240 nm MODIS
  band, which is not among the standardized bands.

---

## Thermal indices

**2 indices** · Landsat (both), MODIS (`lst` only)

| `index=` | Name | Implementation | Sensors |
|---|---|---|---|
| `lst` | Land Surface Temperature (°C) | Landsat Collection 2 `ST_B10` / `ST_B6` × 0.00341802 + 149.0 − 273.15; MODIS `LST_Day_1km` × 0.02 − 273.15, masked outside 200–373 K | Landsat, MODIS |
| `utfvi` | Urban Thermal Field Variance Index (simplified) | `(LST − (0.8·LST − 10·NDVI)) / 10`, masked outside 0–70 °C | Landsat |

`utfvi` is a simplified empirical version, not the published UTFVI (which normalizes LST by
its mean over the scene). Both were offered on Sentinel-2 and Sentinel-3 before v1.6.0, where
they returned empty images; they are now only accepted on sensors with thermal bands.

---

## Sentinel-2 red edge indices

**9 indices** · Sentinel-2 only (bands B5, B6, B7 → `Red_Edge1`, `Red_Edge2`, `Red_Edge3`)

| `index=` | Name | Formula | Reference |
|---|---|---|---|
| `ndre` | Normalized Difference Red Edge | `(Nir − Red_Edge1) / (Nir + Red_Edge1)` | [Gitelson & Merzlyak (1994)](https://doi.org/10.1016/1011-1344(93)06963-4) |
| `cire` | Chlorophyll Index Red Edge | `Nir / Red_Edge1 − 1` | [Gitelson et al. (2003)](https://doi.org/10.1078/0176-1617-00887) |
| `mcari` | Modified Chlorophyll Absorption Ratio Index | `((Red_Edge1 − Red) − 0.2·(Red_Edge1 − Green)) · (Red_Edge1 / Red)` | [Daughtry et al. (2000)](https://doi.org/10.1016/S0034-4257(00)00113-9) |
| `ireci` | Inverted Red Edge Chlorophyll Index | `(Red_Edge3 − Red) / (Red_Edge1 / Red_Edge2)` | [Frampton et al. (2013)](https://doi.org/10.1016/j.isprsjprs.2013.04.007) |
| `mtci` | MERIS Terrestrial Chlorophyll Index | `(Red_Edge2 − Red_Edge1) / (Red_Edge1 − Red)` | [Dash & Curran (2004)](https://doi.org/10.1080/0143116042000274015) |
| `psri` | Plant Senescence Reflectance Index | `(Red − Blue) / Red_Edge2` | [Merzlyak et al. (1999)](https://doi.org/10.1034/j.1399-3054.1999.106119.x) |
| `reip` | Red Edge Inflection Point (nm) | `700 + 40 · ((Red + Red_Edge3)/2 − Red_Edge1) / (Red_Edge2 − Red_Edge1)` | Guyot & Baret (1988) |
| `s2rep` | Sentinel-2 Red Edge Position (nm) | `705 + 35 · ((Red + Red_Edge3)/2 − Red_Edge1) / (Red_Edge2 − Red_Edge1)` | [Frampton et al. (2013)](https://doi.org/10.1016/j.isprsjprs.2013.04.007) |
| `ndci` | Normalized Difference Chlorophyll Index (water) | `(Red_Edge1 − Red) / (Red_Edge1 + Red)` | [Mishra & Mishra (2012)](https://doi.org/10.1016/j.rse.2011.10.016) |

`ndci` is meant for chlorophyll-a in turbid and productive waters; the rest are vegetation
chlorophyll and senescence indices.

---

## Sentinel-3 OLCI water quality indices

**10 indices** · Sentinel-3 only

OLCI band names are standardized as `Blue` (Oa02, 412.5 nm), `Blue2` (Oa03, 442.5 nm),
`Blue_Green` (Oa04, 490 nm), `Green` (Oa05, 510 nm), `Red` (Oa07, 620 nm), `Red2` (Oa08,
665 nm), `Red3` (Oa09, 673.75 nm), `Red_Edge1` (Oa10, 681.25 nm), `Red_Edge2` (Oa11,
708.75 nm) and `Nir` (Oa12, 753.75 nm). The collection holds top-of-atmosphere radiances,
so these are relative indicators for comparing places and dates, not calibrated
concentrations.

| `index=` | Name | Formula | Reference |
|---|---|---|---|
| `oci` | OLCI chlorophyll index | `(Red_Edge1 − Red2) / (Red_Edge1 + Red2)` | [Hu et al. (2012)](https://doi.org/10.1029/2011JC007395) |
| `fluorescence_height` | Fluorescence Line Height | `Red3 − (Red2 + (Red_Edge1 − Red2) · (673.75 − 665) / (681.25 − 665))` | [Gower et al. (2005)](https://doi.org/10.1080/01431160500075857) |
| `red_edge_position` | Red edge position (nm) | `681.25 + 27.5 · ((Red2 + Red_Edge2)/2 − Red_Edge1) / (Red_Edge2 − Red_Edge1)` | [Gower et al. (2005)](https://doi.org/10.1080/01431160500075857) |
| `floating_algae` | Floating algae (OLCI) | `(Nir − Red_Edge2) / (Nir + Red_Edge2)` | [Hu (2009)](https://doi.org/10.1016/j.rse.2009.05.012) |
| `turbidity` | Turbidity | `Red2 / Blue_Green` | [Nechad et al. (2010)](https://doi.org/10.1016/j.rse.2009.11.022) |
| `spm` | Suspended Particulate Matter | `(Red3 − Red2) / (Red3 + Red2)` | [Binding et al. (2005)](https://doi.org/10.1016/j.rse.2004.11.002) |
| `tsi` | Trophic State (spectral proxy) | `(Red2 − Blue_Green) / (Nir − Blue_Green)` | [Carlson (1977)](https://doi.org/10.4319/lo.1977.22.2.0361) |
| `cdom` | Coloured Dissolved Organic Matter | `Blue / Blue_Green` | [Mannino et al. (2008)](https://doi.org/10.1029/2007JC004493) |
| `kd490` | Diffuse attenuation at 490 nm (proxy) | `log(Blue2 / Blue_Green)` | Mueller (2000), SeaWiFS Postlaunch Tech. Report |
| `water_leaving_reflectance` | Water-leaving reflectance (proxy) | `Green / (Blue + Green + Red)` | [Gordon et al. (1988)](https://doi.org/10.1029/JD093iD09p10909) |

`floating_algae` and `tsi` failed on every Sentinel-3 image before v1.6.0 (they looked up a
band that does not exist); they work from v1.6.0 on.

---

## Raw Reflectance Bands

**9 bands** · Sentinel-2 (9), Landsat (6), MODIS (6)

Spectral indices are ratios, so they cancel out multiplicative changes in
brightness: a pixel can keep exactly the same NDVI while its reflectance
drifts. When the question is about the radiometry itself rather than about
vegetation, the bands are available as regular index names.

```python
index='blue'   # or 'green', 'red', 'nir', 'swir1', 'swir2'
index='red_edge1'   # or 'red_edge2', 'red_edge3' — Sentinel-2 only
```

**Units**: surface reflectance, 0 to 1, on every sensor, so the same numeric threshold means
the same thing whichever sensor produced it.

**Not available on Sentinel-3**, whose bands are top-of-atmosphere radiances
with a different band set, and not on Sentinel-1, which has `vv` and `vh`.

**Typical use — radiometric stability**: combined with the dispersion reducers,
the bands map how invariant each pixel is over a time series, which is how
pseudo-invariant features (PIFs) are chosen for relative radiometric
normalization:

```python
# Within-year dispersion of SWIR1, one composite per year
stability = NdviSeasonality(
    roi=roi, start_year=2018, end_year=2025,
    periods=1, sat='S2', index='swir1', key='std'
).get_year_composite()

# How many observations each pixel is based on
n_obs = NdviSeasonality(
    roi=roi, start_year=2018, end_year=2025,
    periods=1, sat='S2', index='swir1', key='count'
).get_year_composite()
```

| Name | Band | S2 | Landsat | MODIS |
|------|------|----|---------|-------|
| `blue` | Blue (~490 nm) | ✓ | ✓ | ✓ |
| `green` | Green (~560 nm) | ✓ | ✓ | ✓ |
| `red` | Red (~665 nm) | ✓ | ✓ | ✓ |
| `nir` | NIR (~840 nm) | ✓ | ✓ | ✓ |
| `swir1` | SWIR 1 (~1610 nm) | ✓ | ✓ | ✓ |
| `swir2` | SWIR 2 (~2190 nm) | ✓ | ✓ | ✓ |
| `red_edge1` | Red Edge 1 (~705 nm) | ✓ | ✗ | ✗ |
| `red_edge2` | Red Edge 2 (~740 nm) | ✓ | ✗ | ✗ |
| `red_edge3` | Red Edge 3 (~783 nm) | ✓ | ✗ | ✗ |

---

## SAR Indices

**7 variables** · Sentinel-1 (VV + VH, IW mode)

`vv` and `vh` are returned in **dB**, the usual scale for backscatter. The polarimetric
indices are defined on **linear power**, so each of them converts VV and VH to linear
(`10^(dB/10)`) before applying its formula. The default preprocessing (`use_sar_ard=True`)
applies terrain flattening and speckle filtering on linear power too.

| `index=` | Name | Formula | Reference |
|---|---|---|---|
| `vv` | VV backscatter (dB) | `VV` | — |
| `vh` | VH backscatter (dB) | `VH` | — |
| `rvi` | Dual-pol Radar Vegetation Index | `4·VH / (VV + VH)` (linear) | [Kim et al. (2012)](https://doi.org/10.1109/LGRS.2011.2174772); [Nasirzadehdizaji et al. (2019)](https://doi.org/10.3390/app9040655) |
| `vv_vh_ratio` | VV/VH ratio | `VV / VH` (linear) | [Mascolo et al. (2016)](https://doi.org/10.1109/TGRS.2016.2585744) |
| `rfdi` | Radar Forest Degradation Index (dual-pol) | `(VV − VH) / (VV + VH)` (linear) | [Mitchard et al. (2012)](https://doi.org/10.5194/bg-9-179-2012) |
| `dpsvi` | Modified Dual-Pol SAR Vegetation Index (DPSVIm) | `(VV² + VV·VH) / √2` (linear) | [dos Santos et al. (2021)](https://doi.org/10.1080/01431161.2021.1959955) |
| `vsdi` | Vegetation Scattering Diversity Index — **experimental** | `√((VV − VH)² + (VV + VH)²)` (dB) | none found |

- **`rvi`** ranges roughly from 0 (bare soil, water) to 1 (dense vegetation).
- **`vv_vh_ratio`** in linear power is equivalent to the difference `VV − VH` in dB.
- **`rfdi`** was defined by Mitchard et al. with HH/HV; Sentinel-1 carries VV/VH, so the
  common dual-pol adaptation is used.
- **`dpsvi`** is the modified form by dos Santos et al. (2021). Periasamy's original DPSVI
  ([2018](https://doi.org/10.1016/j.rse.2018.09.003)) multiplies by a term that depends on the
  maximum VV of the whole scene, so it is not a per-pixel index and changes with the extent.
- **`vsdi`** has no known published definition. It is kept for compatibility, but prefer the
  documented indices for anything that has to be reported.

```{warning}
Before **v1.6.0** the Sentinel-1 collection (served in dB by Earth Engine) was processed and
indexed as if it were linear power: terrain correction and speckle filtering ran on dB values,
and `rvi`, `vv_vh_ratio`, `rfdi` and `dpsvi` were computed on dB. Their values change from
v1.6.0 (e.g. `rvi` 2.42 → 0.70 over the same scene); `vv` and `vh` change on sloping terrain,
where the old terrain correction distorted them.
```

The `normalize=True` option of the SAR index methods returns a z-score (mean 0, standard
deviation 1) instead of the raw value.

---

## ERA5-Land Climate Variables

**47 variables** · ERA5-Land daily aggregates (ECMWF reanalysis, 1950–present, ~11 km)

Use the reducer that matches the variable: `key='mean'` for states (temperature, soil
moisture, pressure), `key='sum'` for fluxes accumulated over the period (precipitation,
evaporation, runoff, radiation), `key='min'`/`'max'` for extremes.

### Temperature (24)

Four variables, each as daily mean, daily minimum and daily maximum, in Kelvin or converted
to °C (`_celsius` suffix):

| Variable | Kelvin | °C |
|---|---|---|
| Air temperature at 2 m | `temperature_2m`, `temperature_2m_min`, `temperature_2m_max` | `temperature_2m_celsius`, `temperature_2m_min_celsius`, `temperature_2m_max_celsius` |
| Dewpoint temperature at 2 m | `dewpoint_temperature_2m`, `_min`, `_max` | `dewpoint_temperature_2m_celsius`, `_min_celsius`, `_max_celsius` |
| Skin (surface) temperature | `skin_temperature`, `_min`, `_max` | `skin_temperature_celsius`, `_min_celsius`, `_max_celsius` |
| Soil temperature, layer 1 (0–7 cm) | `soil_temperature_level_1`, `_min`, `_max` | `soil_temperature_level_1_celsius`, `_min_celsius`, `_max_celsius` |

### Precipitation & water balance (10)

Each in meters of water, or in L/m² (= mm) with the `_lm2` suffix:

| Variable | meters | L/m² |
|---|---|---|
| Total precipitation | `total_precipitation_sum` | `total_precipitation_sum_lm2` |
| Total evaporation (negative = upward flux) | `total_evaporation_sum` | `total_evaporation_sum_lm2` |
| Potential evaporation | `potential_evaporation_sum` | `potential_evaporation_sum_lm2` |
| Total runoff (surface + subsurface) | `runoff_sum` | `runoff_sum_lm2` |
| Surface runoff | `surface_runoff_sum` | `surface_runoff_sum_lm2` |

### Soil moisture (4)

Volumetric water content (m³/m³) by depth: `volumetric_soil_water_layer_1` (0–7 cm),
`volumetric_soil_water_layer_2` (7–28 cm), `volumetric_soil_water_layer_3` (28–100 cm) and
`volumetric_soil_water_layer_4` (100–289 cm).

### Radiation (3)

In J/m²: `surface_solar_radiation_downwards_sum`, `surface_net_solar_radiation_sum` and
`surface_latent_heat_flux_sum`.

### Wind & pressure (3)

`u_component_of_wind_10m` (east–west) and `v_component_of_wind_10m` (north–south), in m/s —
wind speed is `sqrt(u² + v²)` — and `surface_pressure`, in Pa.

### Snow (3)

`snow_depth_water_equivalent` (m of water), and daily snowfall as `snowfall_sum` (m of water)
or `snowfall_sum_lm2` (L/m²).

---

## CHIRPS Precipitation

**1 variable** · CHIRPS daily (1981–present, ~5.5 km, 50°S–50°N)

```python
sat='CHIRPS', index='precipitation'
```

Satellite and station blended daily precipitation in mm/day. Use `key='sum'` for
monthly or seasonal totals. Finer grid than ERA5-Land and calibrated against stations, but
limited to 50°S–50°N.

---

## Usage Examples

### Optical

```python
from ndvi2gif import NdviSeasonality

ndvi = NdviSeasonality(roi=roi, sat='S2', index='ndvi')
evi = NdviSeasonality(roi=roi, sat='MODIS', index='evi', key='mean')
chlorophyll = NdviSeasonality(roi=roi, sat='S2', index='ndre')   # red edge, S2 only
water = NdviSeasonality(roi=roi, sat='Landsat', index='mndwi')
surface_temp = NdviSeasonality(roi=roi, sat='Landsat', index='lst', key='max')
```

### SAR

```python
# All-weather vegetation monitoring; the median also dampens residual speckle
rvi = NdviSeasonality(roi=roi, sat='S1', index='rvi', key='median', periods=12)
```

### Climate

```python
temperature = NdviSeasonality(roi=roi, sat='ERA5', index='temperature_2m_celsius',
                              start_year=1980, end_year=2023, key='mean')
rain_totals = NdviSeasonality(roi=roi, sat='CHIRPS', index='precipitation', key='sum')
root_zone = NdviSeasonality(roi=roi, sat='ERA5',
                            index='volumetric_soil_water_layer_3', key='mean')
```

### Several sensors in one classification

Indices from different sensors can be combined in a single feature stack with
`LandCoverClassifier` (v1.6.0):

```python
from ndvi2gif import LandCoverClassifier

s2 = NdviSeasonality(roi=roi, sat='S2', periods=4, start_year=2024, end_year=2024)
s1 = NdviSeasonality(roi=roi, sat='S1', periods=12, start_year=2024, end_year=2024)

clf = LandCoverClassifier([s2, s1])
clf.create_feature_stack(indices={'S2': ['ndvi', 'ndmi', 'ndre'], 'S1': ['vv', 'vh', 'rvi']})
```

---

## References

**Climate datasets:**
- Muñoz-Sabater, J. (2019). ERA5-Land hourly data from 1981 to present. Copernicus Climate
  Change Service (C3S) Climate Data Store. [DOI:10.24381/cds.e2161bac](https://doi.org/10.24381/cds.e2161bac)
- Funk, C. et al. (2015). The climate hazards infrared precipitation with stations — a new
  environmental record for monitoring extremes. *Scientific Data*, 2, 150066.
  [DOI:10.1038/sdata.2015.66](https://doi.org/10.1038/sdata.2015.66)

The references of each index are linked in its table.

---

## See Also

- [Datasets Reference](datasets.md) - Detailed platform documentation
- [API Reference](api.md) - Complete class documentation
- [Tutorials](../tutorials/basic_ndvi.md) - Step-by-step index usage
