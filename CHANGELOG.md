# Changelog

All notable changes to the `ndvi2gif` package will be documented in this file.

---

## [0.1.5] - 2025-05-26

- Fixing bug with MNDWI index.


## [0.1.4] - 2025-05-25

- Nothing really changes, just a f* problem with release version management.

## [0.1.3] - 2025-05-25

- Nothing really changes, just a f* problem with release version management.

## [0.1.2] - 2025-05-25

### Added

- Complete rework and translation of the README into Markdown format.
- Included new seasonal/statistical methods and updated docstrings in English.
- Added support for region input via Sentinel-2 tiles and Landsat path/row.
- Added rich ROI input documentation with tabular summary.
- Added `deims` dependency as optional to avoid conda forge problems

### Changed

- Clarified the purpose of the library as a broader seasonal analysis tool, not just for GIF generation.
- Cleaned and validated `setup.cfg` and `pyproject.toml`.
- Added extra requirements group for `deims`.

---

## [0.1.1] - 2025-05-21

### Fixed

- Fixed rendering issue in `README.rst` that caused PyPI upload failure.
- Rebuilt and republished the package with correct long description format.

---

## [0.1.0] - 2025-05-21

### Added

- Compatibility with Conda packaging and `conda-forge` ecosystem.
- Included `MANIFEST.in` to ensure `LICENSE` and `README.rst` are bundled in source distribution.
- Improved `README.rst` formatting to comply with PyPI rendering rules.

### Changed

- Switched versioning to semantic 0.x.y style for future compatibility.
- Cleaned and validated metadata to allow upload to both PyPI and Conda Forge.

### Note

- This is a technical release — no changes to the core functionality.

---

## [0.0.9] - 2025-05-20

### Changed

- Version bump to align `setup.cfg`, PyPI and GitHub release.
- No functional changes from version 0.0.7.

---

## [0.0.7] - 2025-05-20

### Added

- New method `get_ndmi()` to compute the Normalized Difference Moisture Index (NDMI).
- New (old) method `get_gif()` to download a gif for the selected index/bands.
- Package structure modernized:
  - Added `setup.cfg` and `pyproject.toml` (PEP 517/518 compliant).
  - Optional removal of legacy `setup.py`.
- Updated dependencies:
  - `geemap` pinned to version `0.29.5`.
  - `numpy` constrained to `<2.0` for compatibility.

### Fixed

- Compatibility issues with recent versions of `geemap`, `xarray`, and `numpy`.
- Resolved import error caused by the removal of `np.unicode_` in NumPy 2.0.

---

## [0.0.6] - 2023-03-10

### Added

- Initial public release of the `ndvi2gif` package.
- Generate seasonal composites and extract statistical summaries from several remote sensing index using Google Earth Engine and geemap.
- Export to animated GIF and GeoTIFF format.

