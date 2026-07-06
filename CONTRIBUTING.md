# Contributing to Ndvi2Gif

First off, thanks for taking the time to contribute! 🎉  
This project grew from a small NDVI-to-GIF helper into a **remote sensing analytics suite**. We warmly welcome contributions of all sizes—from typo fixes to new indices, datasets, and examples.

---

## How Can I Contribute?

- **Fix bugs** or improve error messages.
- **Add a new index** (optical or SAR) — the easiest and most impactful contribution.
- **Add a new dataset** (an Earth Engine ImageCollection properly integrated into the pipeline).
- **Improve documentation** (README, tutorials, notebooks).
- **Add tests** or simplify existing code.
- **Share use-case notebooks** in `examples_notebooks/`.

---

## Contribution Workflow

We follow the standard GitHub **fork → branch → pull request** model. You do
**not** need write access to the main repository.

1. **Fork** the repository to your own GitHub account (click *Fork* at the top
   right of [the repo page](https://github.com/Digdgeo/Ndvi2Gif)).

2. **Clone your fork** and add the original repo as an `upstream` remote so you
   can keep your fork in sync:

   ```bash
   git clone https://github.com/<your-username>/Ndvi2Gif.git
   cd Ndvi2Gif
   git remote add upstream https://github.com/Digdgeo/Ndvi2Gif.git
   ```

3. **Create a feature branch** off `master` (never commit directly to `master`):

   ```bash
   git checkout -b my-new-index      # e.g. add-savi-index, fix-export-crs
   ```

4. **Make your changes**, then commit them with a clear message:

   ```bash
   git add .
   git commit -m "Add SAVI index for S2 and Landsat"
   ```

5. **Keep your branch up to date** with upstream before opening the PR:

   ```bash
   git fetch upstream
   git rebase upstream/master
   ```

6. **Push** the branch to your fork and **open a Pull Request** against
   `Digdgeo/Ndvi2Gif:master`:

   ```bash
   git push origin my-new-index
   ```

   Then open the PR from GitHub and describe what you changed and why.

---

## Development Setup

Set up a development environment inside your cloned fork:

```bash
conda create -n ndvi2gif-dev python=3.11 -y
conda activate ndvi2gif-dev
pip install -e ".[dev]"
```

Authenticate Earth Engine (a Google Cloud **project ID is required**):

```python
import ee
ee.Authenticate()
ee.Initialize(project='your-project-id')
```

---

## Project Structure

- `ndvi2gif.py` → **NdviSeasonality**: core seasonal/statistical engine.  
- `s1_ard.py` → **S1ARDProcessor**: Sentinel-1 preprocessing (ARD, speckle filters, terrain correction).  
- `timeseries.py` → **TimeSeriesAnalyzer**: time series extraction, trend analysis, phenology.  
- `examples_notebooks/` → contributed notebooks and tutorials.  

---

## Adding a New Index

1. Add the index name to the correct set (`optical_indices`, `s1_indices`, etc.).  
2. Implement the function in `NdviSeasonality` or `S1ARDProcessor`.  
3. Register it in the `self.d` dispatch dictionary.  
4. Map it to the appropriate sensor in `self.sensor_indices`.  
5. Add a minimal test and a short example (README or notebook).

**Minimal example**:

```python
def get_myindex(self, image):
    num = image.select("Nir").subtract(image.select("Red"))
    den = image.select("Nir").add(image.select("Red"))
    return num.divide(den).rename("MYINDEX")

# Register it
self.optical_indices.add("myindex")
self.d["myindex"] = self.get_myindex
self.sensor_indices["S2"].add("myindex")
self.sensor_indices["Landsat"].add("myindex")
```

---

## Adding a New Dataset

Contributing new **Earth Engine datasets** is one of the most useful ways to expand Ndvi2Gif.

1. **Find the dataset** in the [Earth Engine catalog](https://developers.google.com/earth-engine/datasets/).  
   Example: `COPERNICUS/S5P/NRTI/L3_NO2` or `NASA/ORNL/DAYMET_V4`.  

2. **Define its ImageCollection** inside `_setup_satellite_collections()` in `NdviSeasonality`.  
   - Give it a short, clear `sat` code (e.g., `"S3"`, `"L8"`, `"MOD09"`).  
   - Apply standard filters (`date`, `bounds`).  

3. **Standardize bands** so that indices can use them.  
   - Rename bands consistently: `Red`, `Nir`, `Green`, `Swir1`, `Blue`, etc.  
   - Apply scale factors if the dataset requires (see Sentinel-2 or MODIS examples).  

4. **Register the dataset** in `sensor_indices` to tell the system which indices it supports.  

5. **Document it**:  
   - Add to README under *Supported Datasets*.  
   - Include the Earth Engine catalog link.  
   - Optionally add a short example in a notebook.  

6. **Test it**:  
   - Run a quick seasonal composite with a small ROI and check that it executes.  
   - Add a smoke test in `tests/` if possible.  

**Tip**: use existing datasets (Sentinel-2 SR, MODIS SR, Landsat C2 L2) as templates — they show the correct pattern for scaling, band mapping, and validation.

---

### Example: Adding a New Dataset (MODIS MOD13Q1)

Suppose you want to add the [MOD13Q1](https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MOD13Q1) dataset (16-day NDVI/EVI, 250m).

1. **Locate the dataset** in the Earth Engine catalog:  
   ID: `MODIS/061/MOD13Q1`  
   Bands: `NDVI`, `EVI`, `sur_refl_b01`, `sur_refl_b02`, `sur_refl_b03`, etc.  

2. **Create a scaling helper** (if needed):  
   ```python
   def _scale_mod13(self, img):
       # MODIS scale factor = 0.0001
       return img.multiply(0.0001).copyProperties(img, ["system:time_start"])
   ```

3. **Update `_setup_satellite_collections()`** in `NdviSeasonality`:  
   ```python
   elif self.sat == "MOD13":
       col = ee.ImageCollection("MODIS/061/MOD13Q1") \
                 .filterDate(self.start_date, self.end_date) \
                 .filterBounds(self.roi)
       col = col.map(self._scale_mod13)
       # Standardize band names
       col = col.select(
           ["sur_refl_b01", "sur_refl_b02", "sur_refl_b03", "sur_refl_b07"],
           ["Red", "Nir", "Blue", "Swir2"]
       )
       self.collection = col
   ```

4. **Register it in `sensor_indices`** (at class init):  
   ```python
   self.sensor_indices["MOD13"] = {"ndvi", "evi", "ndwi", "msi"}
   ```

5. **Test it quickly** in a notebook:  
   ```python
   proc = NdviSeasonality(
       roi=ee.Geometry.Point([-3.7, 40.4]).buffer(5000),  # Madrid area
       sat="MOD13",
       periods=12,
       start_year=2020,
       end_year=2021,
       index="ndvi"
   )
   img = proc.get_year_composite(2020)
   print(img.bandNames().getInfo())
   ```

6. **Document it**:  
   - Add MOD13Q1 to README’s *Supported Datasets* table with the EE link.  
   - Add a short notebook in `examples_notebooks/` showing usage.

---

## Testing

We use **pytest**. Add lightweight smoke tests:
```python
def test_ndvi_smoke():
    proc = NdviSeasonality(roi=ee.Geometry.Point([0,0]).buffer(1000),
                           periods=12, start_year=2023, end_year=2024,
                           sat='S2', index='ndvi')
    img = proc.get_period_composite(2023, 0)
    assert isinstance(img, ee.Image)
```

Run tests with:
```bash
pytest -q
```

---

## Pull Request Checklist

Before opening the PR (see [Contribution Workflow](#contribution-workflow) above),
make sure you have:

- [ ] Formatted with `black` + `isort`.
- [ ] Linted with `ruff` or `flake8`.
- [ ] Added/updated docstrings and a `CHANGELOG` entry.
- [ ] Added or updated a notebook if your change introduces a new feature.
- [ ] Run the test suite (`pytest -q`) and confirmed it passes.

---

## Issues

Please include:
- Code snippet (or notebook cell).  
- Dataset/index/parameters used.  
- Expected vs. actual behavior.  
- Environment (Python version, OS, Earth Engine version).  

Labels used: `bug`, `index`, `dataset`, `docs`, `enhancement`.

---

✨ Thanks for contributing to Ndvi2Gif! ✨
