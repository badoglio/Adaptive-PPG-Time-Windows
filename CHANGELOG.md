# Changelog

All notable changes to this project are listed here. Versions follow [Semantic Versioning](https://semver.org/).

## [0.2.0] – 2026-09-30

A rewrite of the analysis pipeline as an installable package. Results are not comparable with 0.1.0:
features are now averaged over the adaptive windows (ensemble templates), not computed beat by beat.

### Upgrading from 0.1.0

* The flat scripts (`app.py`, `adaptive_engine.py`, `feature_extraction.py`, `io_loader.py`,
  `preprocessing.py`, `visualization.py`) are replaced by the `adaptive_ppg` package in `src/` and by
  `app/streamlit_app.py`.
* Install with `pip install -e ".[app]"` instead of installing the dependencies by hand, then start the
  dashboard with `adaptive-ppg dashboard` (or `streamlit run app/streamlit_app.py`) instead of
  `streamlit run app.py`. The step-by-step guides are in [`docs/`](docs/INSTALL_en.md).
* An existing `venv/` folder from 0.1.0 can be deleted; the guides create `.venv/`.

### Added

* `adaptive_ppg` package with `pyproject.toml`, the `adaptive-ppg` command (`run`, `benchmark`, `sweep`,
  `init-config`, `dashboard`), batch runs, JSON/YAML configuration and presets.
* Ensemble averaging per feature category: each category (macro, time_volume, derivatives) has its own
  adaptive window; features are extracted from the ensemble template, with CI95 and template SQI.
* Resampling to a fixed processing rate (`target_fs = 125` Hz by default, `--target-fs`), sub-sample
  fiducial points and a fixed 500 Hz template grid, so that results do not depend on the input rate.
* Automatic detection of inverted signals (`polarity = "auto"`, `--polarity`).
* CSV dialect detection (`;` separator, decimal comma) and a least-squares sampling-rate estimate for
  rounded time columns.
* Synthetic signals with known ground truth and a benchmark against fixed, EMA and per-beat estimators.
* Pulse arrival time (PAT) from an ECG channel.
* Sample recordings `Sample1.CSV` and `Sample2.CSV` (Pulse Transit Time PPG Dataset, ODbL 1.0; IR light,
  distal phalanx, seated) with provenance in `sample_data/README.md`.
* Test suite (pytest), ruff configuration and GitHub Actions CI on Linux and Windows.
* Type hints (`py.typed`), `LICENSE` (GPL-3.0-or-later), `.gitattributes`, this changelog, and
  installation guides in English, Italian and Japanese.

### Changed

* README (English and Japanese) rewritten to match the code: method, validation, configuration, package
  layout, sample data and license.

### Removed

* `sample_data/S01.csv`: its provenance and license are unknown, so it is no longer distributed.

## [0.1.0] – initial release

Streamlit dashboard and flat analysis scripts (commit `548ca91`).

[0.2.0]: https://github.com/badoglio/Adaptive-PPG-Time-Windows/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/badoglio/Adaptive-PPG-Time-Windows/releases/tag/v0.1.0
