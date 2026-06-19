# Changelog

All notable changes to LoopDetect will be documented in this file.

## [0.2.0] - 2026-06-20

### Changed

* Modernized the Python package layout by moving the source package into `src/loopdetect/`.
* Standardized the import package name to lowercase `loopdetect`.
* Updated package configuration from the older `setup.py`-based structure to `pyproject.toml`.
* Updated Sphinx documentation configuration so the documentation can import the package from the new `src/` layout.
* Updated documentation examples to use the lowercase import style:

  ```python
  import loopdetect as ld
  ```

  or:

  ```python
  import loopdetect.core as ld
  ```

### Added

* Added `CHANGELOG.md` to track release history.
* Added optional test structure under `tests/` for package import checks, core-function checks, and package-data checks.
* Added development/testing environment guidance for creating a clean Conda environment.

### Preserved

* Kept the core LoopDetect functions available through:

  ```python
  import loopdetect.core as ld
  ```
* Preserved package data under the installable package directory:

  ```text
  src/loopdetect/data/
  ```

### Notes

* `LoopDetect` remains the project display name.
* The Python import name should be lowercase:

  ```python
  import loopdetect
  ```
