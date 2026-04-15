# Installation

This guide will help you install Ndvi2Gif and set up your environment.

## Installation Options

### Option 1: Install with pip (Recommended)

The simplest way to install Ndvi2Gif is using pip:

```bash
pip install ndvi2gif
```

### Option 2: Install with conda

If you prefer conda/mamba:

```bash
conda install -c conda-forge ndvi2gif
```

### Option 3: Install from source (Development)

For the latest development version or to contribute:

```bash
# Clone the repository
git clone https://github.com/Digdgeo/Ndvi2Gif.git
cd Ndvi2Gif

# Create a conda environment
conda create -n ndvi2gif python=3.11 -y
conda activate ndvi2gif

# Install in editable mode
pip install -e ".[dev]"
```

## Setting Up Your Environment

### Create a Dedicated Environment

We recommend creating a dedicated environment for remote sensing work:

```bash
# Using conda
conda create -n geo_analysis python=3.11 -y
conda activate geo_analysis
pip install ndvi2gif

# Or using venv
python -m venv ndvi2gif_env
source ndvi2gif_env/bin/activate  # On Windows: ndvi2gif_env\Scripts\activate
pip install ndvi2gif
```

### Installing with JupyterLab

If you want to use Ndvi2Gif in Jupyter notebooks:

```bash
conda create -n ndvi2gif python=3.11 -y
conda activate ndvi2gif
conda install -c conda-forge jupyterlab
pip install ndvi2gif
```

## Dependencies

Ndvi2Gif will automatically install the following key dependencies:

- **earthengine-api** (≥0.1.347): Google Earth Engine Python API
- **geemap** (≥0.29.5): Interactive mapping with Earth Engine
- **geopandas**: Geospatial data processing
- **numpy** (>=1.24): Numerical operations
- **pandas**: Data manipulation
- **matplotlib**: Visualization
- **scipy**: Scientific computing
- **scikit-learn** (≥1.0): Machine learning

## Verify Installation

Test that everything is installed correctly:

```python
import ndvi2gif
print(f"Ndvi2Gif version: {ndvi2gif.__version__}")

# Check all modules are available
from ndvi2gif import (
    NdviSeasonality,
    S1ARDProcessor,
    TimeSeriesAnalyzer,
    LandCoverClassifier,
    HydroperiodAnalyzer,
)

print("✓ All modules loaded successfully!")
```

Expected output:
```
Ndvi2Gif version: 1.2.2
✓ All modules loaded successfully!
```

## Troubleshooting

### Common Issues

#### Issue: GDAL/Fiona installation problems

On some systems, installing geopandas can fail due to GDAL dependencies. Try:

```bash
# On Ubuntu/Debian
sudo apt-get install libgdal-dev

# Then install geopandas separately
pip install geopandas

# Finally install ndvi2gif
pip install ndvi2gif
```

Or use conda which handles binary dependencies better:

```bash
conda install -c conda-forge geopandas ndvi2gif
```

#### Issue: Permission errors on Linux/Mac

Use the `--user` flag:

```bash
pip install --user ndvi2gif
```

### Getting Help

If you encounter issues:

1. Check the [FAQ](../reference/faq.md)
2. Search [GitHub Issues](https://github.com/Digdgeo/Ndvi2Gif/issues)
3. Ask in [GitHub Discussions](https://github.com/Digdgeo/Ndvi2Gif/discussions)

## Next Steps

Once installed, proceed to [Authentication](authentication.md) to set up your Google Earth Engine credentials.
