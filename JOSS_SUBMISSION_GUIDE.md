# JOSS Submission Guide for ndvi2gif

## 📋 JOSS Requirements Checklist

### Essential Requirements ✅

- [x] **Open Source License** → MIT License ✓
- [x] **Public Repository** → GitHub (https://github.com/Digdgeo/Ndvi2Gif) ✓
- [x] **Version Control** → Git ✓
- [x] **Documentation** → Jupyter Book (comprehensive) ✓
- [x] **Installation Instructions** → README + book/getting_started/installation.md ✓
- [x] **Example Usage** → Multiple tutorials and notebooks ✓
- [x] **Community Guidelines** → Contributing, Code of Conduct ✓
- [x] **Substantial Scholarly Effort** → 88 variables, 7 platforms, ML, time series ✓
- [x] **Research Applications** → Remote sensing, climate, ecology ✓

### Recommended (Strengthen Submission) 📝

- [x] **Automated Tests** → pytest suite ✓
- [x] **CI/CD** → GitHub Actions ✓
- [x] **PyPI Distribution** → Published v1.0.0 ✓
- [x] **Conda Distribution** → conda-forge ✓
- [ ] **Statement of Need** → Draft paper (TO DO)
- [ ] **Paper Draft** → paper.md (TO DO)
- [ ] **paper.bib** → References file (TO DO)
- [ ] **Author ORCID** → Add to paper.md (RECOMMENDED)
- [ ] **Reviewers Suggestions** → Optional but helpful

---

## 📄 Paper Structure (paper.md)

JOSS papers are typically 1-2 pages. Create `paper/paper.md`:

```markdown
---
title: 'Ndvi2Gif: A Python Library for Multi-Seasonal Remote Sensing and Climate Analysis with Google Earth Engine'
tags:
  - Python
  - remote sensing
  - Google Earth Engine
  - time series analysis
  - climate reanalysis
  - vegetation monitoring
  - machine learning
  - NDVI
authors:
  - name: Diego García Díaz
    orcid: 0000-0000-0000-0000  # ADD YOUR ORCID
    affiliation: "1, 2"
affiliations:
 - name: Estación Biológica de Doñana (EBD-CSIC), Spain
   index: 1
 - name: eLTER Network / SUMHAL Project
   index: 2
date: DD Month 2025
bibliography: paper.bib
---

# Summary

Ndvi2Gif is a comprehensive Python library that democratizes access to global satellite and climate reanalysis data through Google Earth Engine (GEE). While originally designed for vegetation monitoring via seasonal NDVI composites, the library has evolved into a complete remote sensing and climate analysis platform, supporting 88 variables across 7 satellite and reanalysis platforms spanning from 1950 to present.

The library addresses a critical gap in the remote sensing community: the need for user-friendly, reproducible tools that bridge the gap between Earth Engine's powerful data catalog and practical scientific workflows. By providing high-level abstractions for common remote sensing tasks—temporal compositing, index calculation, time series analysis, and machine learning classification—Ndvi2Gif enables researchers to focus on scientific questions rather than technical implementation.

# Statement of Need

Remote sensing and climate reanalysis datasets offer unprecedented opportunities for Earth system monitoring, yet accessing and processing these data remains technically challenging. Google Earth Engine provides a powerful cloud-based platform, but requires significant programming expertise and domain knowledge. Ndvi2Gif fills this gap by providing:

1. **Unified Multi-Platform Access**: Seamless integration of optical (Sentinel-2/3, Landsat, MODIS), SAR (Sentinel-1), and climate reanalysis data (ERA5-Land, CHIRPS) through a consistent API.

2. **Temporal Analysis Workflows**: Automated generation of seasonal/monthly composites with multiple statistical reducers (max, median, percentile, sum, min), essential for cloud-contaminated optical imagery and precipitation accumulation.

3. **Climate-Vegetation Integration**: First tool to seamlessly combine vegetation indices with climate variables (temperature, precipitation, soil moisture) from ERA5-Land (47 variables, 1950-present) and high-resolution CHIRPS precipitation (1981-present, ~5.5km).

4. **Intelligent Data Handling**: Automatic detection of data types (vegetation vs. climate) with appropriate statistical methods and visualizations—phenology metrics for vegetation, seasonal statistics for climate.

5. **Production-Ready ML Pipeline**: Complete land cover classification workflows with 8 algorithms (Random Forest, SVM, CART, Naive Bayes, Gradient Tree Boost, K-means, Cascade K-means, LDA) including accuracy assessment and feature importance.

6. **Reproducible Science**: Comprehensive Jupyter Book documentation with executable examples, ensuring reproducibility and lowering the barrier to entry for students and researchers.

The library has been validated through integration into the eLTER (European Long-Term Ecosystem Research) network and the SUMHAL project, demonstrating its utility for long-term ecological monitoring across diverse ecosystems.

# Key Features

## Multi-Sensor Data Access

Ndvi2Gif provides standardized access to 7 satellite and reanalysis platforms:

- **Optical**: Sentinel-2 (2015-present, 10m), Sentinel-3 OLCI (2016-present, 300m), Landsat 4-9 (1982-present, 30m), MODIS (2000-present, 500m)
- **SAR**: Sentinel-1 (2014-present, 10m) with advanced preprocessing (speckle filtering, terrain correction, orbit control)
- **Climate**: ERA5-Land (1950-present, ~11km, 47 variables), CHIRPS precipitation (1981-present, ~5.5km)

## Extensive Variable Library

The library supports 88 variables spanning vegetation, water quality, SAR, and climate domains:

- **40+ Optical Indices**: NDVI, EVI, NDWI, MNDWI, SAVI, GNDVI, NDMI, chlorophyll indices, water quality parameters
- **7 SAR Indices**: RVI, DPSVI, RFDI, VSDI with dual polarization (VV+VH)
- **47 ERA5 Climate Variables**: Temperature (with min/max/Celsius variants), precipitation (meters/L/m²), soil moisture (4 layers), radiation, wind, pressure, snow
- **1 CHIRPS Variable**: High-quality satellite + station precipitation

## Advanced Analytics

### Time Series Analysis
- Trend detection (Mann-Kendall, Sen's slope, linear regression)
- Phenology metrics (SOS, EOS, POS, length of season, growth/senescence rates)
- Climate statistics (seasonal means, extremes, water balance)
- Interactive multi-panel dashboards

### Machine Learning Classification
- Supervised: Random Forest, SVM, CART, Naive Bayes, Gradient Tree Boost
- Unsupervised: K-means, Cascade K-means, LDA
- Multi-temporal feature stacks
- Accuracy assessment with confusion matrices
- Feature importance analysis

### Flexible ROI Options
- Shapefiles and GeoJSON
- eLTER DEIMS site IDs (via deimsPy integration)
- Sentinel-2 tile codes
- Landsat path/row
- Hand-drawn geometries
- Point coordinates with buffers

## Validation & Testing

The library includes comprehensive test coverage:
- Unit tests for all major components
- Integration tests for Earth Engine operations
- Continuous integration via GitHub Actions
- Validated through eLTER and SUMHAL projects

# Research Applications

Ndvi2Gif has been designed for and used in:

- **Agricultural Monitoring**: Crop phenology, yield estimation, irrigation assessment
- **Wetland Dynamics**: Water extent changes, vegetation health, hydrological cycles
- **Climate Analysis**: Temperature trends, precipitation patterns, drought monitoring, soil moisture dynamics
- **Ecological Studies**: Vegetation phenology, habitat monitoring, biodiversity indicators
- **Land Cover Mapping**: Multi-temporal classification, change detection
- **Water Quality**: Chlorophyll-a estimation, turbidity monitoring, eutrophication assessment

The integration with ERA5-Land and CHIRPS (v1.0.0) enables novel research on climate-vegetation interactions, drought impacts, and phenological responses to climate variables at regional to global scales.

# Performance & Scalability

Ndvi2Gif leverages Google Earth Engine's cloud-based infrastructure:

- **Global Coverage**: Process any region from local (10m) to continental scales
- **Temporal Depth**: Access data from 1950 (ERA5) to near-real-time
- **Parallel Processing**: Earth Engine's distributed computing
- **Export Options**: GeoTIFF, Google Drive, Earth Engine Assets
- **Memory Efficient**: Server-side operations avoid local memory constraints

# Community & Documentation

The project maintains high standards for documentation and community engagement:

- **Jupyter Book**: Comprehensive interactive tutorials with executable examples
- **API Documentation**: Complete NumPy-style docstrings with scientific references
- **Example Notebooks**: Real-world use cases (agriculture, wetlands, drought, phenology)
- **Active Development**: Regular releases with semantic versioning
- **Open Community**: GitHub discussions, issue tracking, contribution guidelines

# Comparison with Existing Tools

While several Python libraries exist for Earth Engine (geemap, eemont, geetools), Ndvi2Gif is unique in:

1. **Domain Focus**: Specifically designed for time series and climate analysis workflows common in ecology and remote sensing research
2. **Climate Integration**: First to seamlessly combine vegetation indices with ERA5 and CHIRPS climate data
3. **ML Pipeline**: Complete classification workflow from training data to accuracy assessment
4. **Phenology Analysis**: Built-in phenology extraction with multiple methods
5. **Publication-Ready**: Automated generation of analysis-ready visualizations

# Future Development

Planned enhancements include:

- Additional climate datasets (TerraClimate, MODIS LST daily)
- Spatial trend analysis expansion
- Harmonic regression for phenology
- Deep learning classification options
- Interactive web-based visualization dashboard

# Acknowledgements

This work was supported by the eLTER (European Long-Term Ecosystem Research) network and the SUMHAL project. The author thanks the Google Earth Engine team for providing the platform, and the ECMWF and UCSB Climate Hazards Center for making ERA5-Land and CHIRPS datasets publicly available.

# References
```

---

## 📚 Bibliography (paper.bib)

Create `paper/paper.bib`:

```bibtex
@article{Gorelick2017,
  title={Google Earth Engine: Planetary-scale geospatial analysis for everyone},
  author={Gorelick, Noel and Hancher, Matt and Dixon, Mike and Ilyushchenko, Simon and Thau, David and Moore, Rebecca},
  journal={Remote Sensing of Environment},
  volume={202},
  pages={18--27},
  year={2017},
  doi={10.1016/j.rse.2017.06.031}
}

@article{Rouse1974,
  title={Monitoring vegetation systems in the Great Plains with ERTS},
  author={Rouse Jr, John W and Haas, Robert H and Schell, John A and Deering, Donald W},
  journal={NASA Special Publication},
  volume={351},
  pages={309},
  year={1974}
}

@article{Huete2002,
  title={Overview of the radiometric and biophysical performance of the MODIS vegetation indices},
  author={Huete, Alfredo and Didan, Kamel and Miura, Tomoaki and Rodriguez, E Patricia and Gao, Xiang and Ferreira, Laerte G},
  journal={Remote Sensing of Environment},
  volume={83},
  number={1-2},
  pages={195--213},
  year={2002},
  doi={10.1016/S0034-4257(02)00096-2}
}

@article{McFeeters1996,
  title={The use of the Normalized Difference Water Index (NDWI) in the delineation of open water features},
  author={McFeeters, Stuart K},
  journal={International Journal of Remote Sensing},
  volume={17},
  number={7},
  pages={1425--1432},
  year={1996},
  doi={10.1080/01431169608948714}
}

@article{Funk2015,
  title={The climate hazards infrared precipitation with stations—a new environmental record for monitoring extremes},
  author={Funk, Chris and Peterson, Pete and Landsfeld, Martin and Pedreros, Diego and Verdin, James and Shukla, Shraddhanand and Husak, Gregory and Rowland, James and Harrison, Laura and Hoell, Andrew and others},
  journal={Scientific Data},
  volume={2},
  number={1},
  pages={1--21},
  year={2015},
  doi={10.1038/sdata.2015.66}
}

@article{Munoz2021,
  title={ERA5-Land: A state-of-the-art global reanalysis dataset for land applications},
  author={Mu{\~n}oz-Sabater, Joaqu{\'\i}n and Dutra, Emanuel and Agust{\'\i}-Panareda, Anna and Albergel, Cl{\'e}ment and Arduini, Gabriele and Balsamo, Gianpaolo and Boussetta, Souhail and Choulga, Margarita and Harrigan, Shaun and Hersbach, Hans and others},
  journal={Earth System Science Data},
  volume={13},
  number={9},
  pages={4349--4383},
  year={2021},
  doi={10.5194/essd-13-4349-2021}
}

@software{geemap2020,
  title={geemap: A Python package for interactive mapping with Google Earth Engine},
  author={Wu, Qiusheng},
  year={2020},
  publisher={GitHub},
  url={https://github.com/gee-community/geemap},
  doi={10.21105/joss.02305}
}

@article{deims2024,
  title={deimsPy: A Python package for accessing the DEIMS-SDR API},
  author={Wohner, Christoph and Peterseil, Johannes},
  journal={SoftwareX},
  volume={25},
  pages={101602},
  year={2024},
  doi={10.1016/j.softx.2023.101602}
}
```

---

## 🚀 Submission Checklist

### Before Submitting

1. **Create paper folder**:
   ```bash
   mkdir paper
   ```

2. **Write paper.md** (use template above, customize)

3. **Create paper.bib** (add your specific references)

4. **Add ORCID** to paper.md

5. **Build and test Jupyter Book**:
   ```bash
   cd book
   jupyter-book build .
   ```

6. **Verify all links work** in Jupyter Book

7. **Create codemeta.json** (optional but recommended):
   ```bash
   pip install codemetapy
   codemetapy setup.cfg > codemeta.json
   ```

8. **Update README** to mention JOSS submission

### Submission Process

1. **Pre-submission inquiry** (optional):
   - Open issue at https://github.com/openjournals/joss-reviews/issues
   - Tag with `pre-review`
   - Get editor feedback before formal submission

2. **Formal submission**:
   - Go to https://joss.theoj.org/papers/new
   - Fill in repository URL
   - System will check requirements
   - Submit

3. **Review process**:
   - Editor assigns reviewers (usually 2)
   - Reviewers check functionality, documentation, paper
   - Address reviewer comments via GitHub
   - Iterate until accepted

### Expected Timeline

- Pre-submission discussion: 1-2 weeks (optional)
- Review assignment: 1-2 weeks
- Review process: 2-8 weeks
- Revisions: Variable
- **Total**: 1-3 months typically

---

## 📊 Strengthening the Submission

### Highlight Unique Contributions

1. **First climate-vegetation integration tool** for Earth Engine
2. **Largest variable collection** (88 vars) with unified API
3. **Intelligent data type handling** (climate vs vegetation)
4. **eLTER network integration** (DEIMS-SDR)
5. **Publication-ready ML pipeline**

### Potential Reviewer Suggestions

Consider suggesting reviewers who:
- Work with Earth Engine / geemap
- Do remote sensing + climate analysis
- Have JOSS review experience
- Are NOT collaborators

Examples (don't know them personally):
- Dr. Qiusheng Wu (geemap developer)
- Remote sensing researchers from ecology/climate communities

### Additional Materials

- **YouTube demo** (optional but impactful)
- **Binder/Colab links** for interactive examples
- **Published studies** using ndvi2gif (if available)
- **User testimonials** from eLTER/SUMHAL

---

## 📝 Tips for Success

1. **Clear statement of need**: Emphasize gap filled
2. **Target audience**: Researchers + students + practitioners
3. **Research applications**: Cite actual use cases
4. **Documentation**: Highlight Jupyter Book quality
5. **Community**: Show active development, CI/CD, testing
6. **Comparison**: Differentiate from geemap/eemont
7. **Future-proof**: Mention planned features

## 🎯 Expected Outcome

**Strong submission** given:
- ✅ Comprehensive functionality (88 variables, 7 platforms)
- ✅ Excellent documentation (Jupyter Book)
- ✅ Active development (v1.0.0 stable release)
- ✅ Real-world validation (eLTER, SUMHAL)
- ✅ Testing + CI/CD
- ✅ PyPI + Conda distribution
- ✅ Clear research applications

JOSS acceptance rate is ~70%, and ndvi2gif meets or exceeds all criteria.

---

## 📧 Questions?

- JOSS FAQ: https://joss.readthedocs.io/en/latest/submitting.html
- JOSS Review Criteria: https://joss.readthedocs.io/en/latest/review_criteria.html
- Example Papers: Browse accepted JOSS papers in remote sensing/GIS

Good luck! 🍀
