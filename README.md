# Adaptive PPG Windows

> **Morphological analysis of the photoplethysmogram (PPG) with variable-width ensemble averaging, driven by the non-stationarity of the pulse shape**

[日本語版 README](README_JA.md)

This repository contains a Python package (`adaptive_ppg`) and a Streamlit dashboard. Together they extract 33 morphological features from PPG signals. Features are not averaged over fixed intervals. Each feature is averaged over an **ensemble window whose length follows how fast the pulse shape is changing**:

* long windows (up to `T_max`) while the morphology is stable, for noise suppression;
* short windows (down to `T_min` beats) during transients, for temporal resolution.

The window length, the overlap and the causal anchoring are computed per feature category. The analysis does not depend on the input sampling rate: every recording is resampled to a common processing rate. The result is validated against synthetic recordings with a known ground truth, and on real recordings from two public datasets.

---

## Rationale

The trade-off between **stationarity** and **temporal resolution** is usually settled once, by picking a fixed window (e.g. 30 s, or the 5 min of the 1996 HRV Task Force). Autonomic responses are neither stationary nor synchronous:

* vagal effects act within seconds;
* sympathetic vasoconstriction develops over tens of seconds;
* recovery can take minutes.

A single fixed window blurs fast transients, or it stays noisy during long stable phases. This project inverts the choice: **the measured non-stationarity of the pulse morphology sets the window length `T_w` and overlap `O_w`**.

---

## Method

```
raw PPG ─► resampling to 125 Hz ─► polarity check
        ─► dual-band filtering ─► beat detection (Elgendi) + sub-sample tangent onsets ─► per-beat SQI
        ─► normalized beats (peak-aligned resampling to 128 points)
        ─► morphological variability (drift_z) ─► H ∈ [0, 1] ─► T_w^c, O_w per category
        ─► ensemble templates per window ─► features (+ CI95) ─► 1 Hz feature grid
```

### 1. Preprocessing and signal quality (`preprocessing.py`)

* **Common processing rate.** The signal and its auxiliary channels (e.g. ECG) are resampled to `target_fs` (125 Hz by default) with a polyphase anti-aliasing filter (`scipy.signal.resample_poly`, rational ratio, so the output rate is exact).
  * Rates within 1% of the target are left untouched; `target_fs=None` keeps the native rate.
  * Clipping is detected on the native signal, and the mask is carried over to the new time base.
  * When the input is upsampled, the filter cut-offs stay limited by the *input* Nyquist frequency, so upsampling does not widen the band.
  * Why 125 Hz: it covers the 30 Hz morphology band with margin, it is the rate of BIDMC/MIMIC, and recordings from different devices are then analysed on the same time base. The signals held in memory and exported are smaller by the ratio of the rates. The run time gains little, because it is dominated by the template stage, which works on a fixed grid (below): `Sample1.CSV` takes 0.33 s at the native 500 Hz and 0.32 s at 125 Hz.
  * How much the rate matters: on `Sample1.CSV` and `Sample2.CSV`, analysis at 125 Hz differs from the native 500 Hz analysis by a median of 0.10 per-beat SD per feature (0.03–0.05 SD on 10 s windows). At 250 Hz the difference is 0.02 SD, at 100 Hz 0.14 SD and at 64 Hz 0.3 SD. Rates below 100 Hz are not recommended for derivative features.
* **Polarity.** Some sensors store raw photodiode counts (more light = higher counts), which gives an inverted pulse. With `polarity="auto"` (default), the pipeline uses the sign of the skewness of the first derivative: an upright PPG rises fast and decays slowly, so its skewness is positive. If it is negative, the signal is flipped and a warning is recorded. `"normal"` and `"inverted"` force the polarity. `raw_signal` keeps the original polarity, so the DC level and the perfusion index are unaffected.
* **Two zero-phase Butterworth branches (SOS):**
  * detection: 0.5–8 Hz, tuned for the Elgendi (2013) systolic-peak detector;
  * morphology: 0.5–30 Hz, which preserves the upstroke and the second-derivative waves.
  * Cut-offs above Nyquist are clamped, and a warning is recorded.
* **Onsets:** by the **tangent method**, the intersection of the maximum-slope tangent with the horizontal through the preceding minimum. This is more reproducible than the minimum itself, which drifts with diastolic run-off. `onset_method="minimum"` is available.
* **Sub-sample fiducials.** The tangent onset is kept as a fractional sample position, and peaks and minima are refined by parabolic interpolation (`onset_positions`, `peak_positions`, `end_positions`; `onsets`/`peaks`/`ends` remain the integer indices). Beat times, IBIs and timing features are therefore not quantized to `1 / fs` (8 ms at 125 Hz).
* **Per-beat causal SQI.** A beat is rejected if:
  * its rate lies outside `[min_bpm, max_bpm]`;
  * its IBI jumps against the trailing median by more than `max_ibi_jump`;
  * its correlation with a running reference beat is below `min_template_corr`;
  * its amplitude ratio falls outside `amplitude_ratio_limits`;
  * it is clipped.
* **Beat normalization:** each beat is scaled to `[0, 1]` and resampled piecewise (PCHIP). The systolic peak always sits at 30% of the 128 points, so the variability reflects changes of shape, not of timing.
* **Detector:** the default is a native Elgendi detector. `detector="biosppy"` is optional (`pip install -e ".[biosppy]"`); BioSPPy imports `peakutils` without declaring it, and the extra installs it.

### 2. Morphological variability and H (`engine.py`)

The default metric, **`drift_z`**, measures **non-stationarity**, not beat-to-beat scatter:

```
drift_k   = RMS( mean(beats k-N+1..k) − mean(beats k-2N+1..k-N) )          N = N_past (10)
drift_z_k = drift_k / ( sqrt(2/N) · pooled within-block RMS deviation )
```

When the signal is stationary, the two block means differ only by noise, and `drift_z ≈ 1` *at any noise level*. A sustained change of shape raises it in proportion to the change measured in noise units.

`drift_z` is mapped to `H ∈ [0, 1]` with a **fixed scale**: `H = 0` at `drift_z = 1` and `H = 1` at `drift_z = 3` (`fixed_scale=(1, 3)`). `H` is held at 0 until `2N` valid beats are available (warm-up).

Other options:
* metrics `drift`, `successive` (RMS difference of consecutive beats) and `dispersion`;
* data-driven normalizations: `calibration` (p5–p95 of the first 60 s), `rolling` and `global`;
* regressing out the IBI-driven part of the variability (`ibi_correction`).

`successive` with `calibration` is the method this README originally described; it is kept as the `legacy` preset. It has two weaknesses:
* it mostly measures beat-to-beat noise (on the synthetic benchmark, 0.030 during a transient against 0.026 at baseline);
* the percentile normalization places a resting recording at `H ≈ 0.5` by construction.

### 3. Window length and overlap

For each category `c`:

```
s(H)   = 1 / (1 + exp(k · (H − θ)))                   (sigmoid; "linear": 1 − H)      k = 10, θ = 0.3
T_w^c  = T_min^c + (T_max^c − T_min^c) · s(H)          T_max^c [beats] = T_max_sec^c / local IBI
T_w^c  = max(round(T_w^c), T_crit^c)                   counted in valid beats
O_w    = O_max − (O_max − O_min) · s(H)                shorter windows overlap more (0.25 → 0.85)
```

* **Rate limit:** window growth is limited to `max_expand_per_beat` (1 beat per beat), while contraction is immediate. The window therefore reacts to a transient at once without oscillating when `H` is noisy.
* **Placement:** consecutive anchors are `round(T_w · (1 − O_w))` beats apart.
* **Anchoring:** windows are **trailing** by default. Every value at time `t` depends only on beats up to `t`, so the output is causal. `StreamingAdaptiveEngine` reproduces the batch result exactly, beat by beat. `anchor="centered"` gives offline, zero-lag windows.
* **Coverage:** the last beat is always an anchor, so the end of the recording is covered.

| Category | Features | `T_min` | `T_max` | `T_crit` |
| :--- | :--- | :---: | :---: | :---: |
| **macro** (17) | peak amplitude, pulse widths at 10/25/50/75 %, total/systolic/diastolic area, area ratio, inflection-point area ratio, inflection amplitude, reflection index, max systolic / decay slope, slope ratio, DC level, perfusion index | 3 beats | 30 s | 3 |
| **time_volume** (10) | pulse duration, heart rate, crest time, decay time, duty cycle, notch time, ΔT_DVP, stiffness index (needs subject height), PAT to foot / to peak (needs ECG) | 3 beats | 30 s | 3 |
| **derivatives** (6) | SDPTG b/a, c/a, d/a, e/a, aging index (b−c−d−e)/a, (b−e)/a | 5 beats | 60 s | 10 |

### 4. Ensemble templates and features (`features.py`)

For each window, the valid beats are:
1. aligned on their maximum-slope point, located with sub-sample precision;
2. sampled by cubic B-spline interpolation on a **fixed template grid** (`template_fs`, 500 Hz by default), whatever the processing rate;
3. corrected by removing a linear foot-to-foot baseline;
4. scaled to unit peak and averaged;
5. rescaled by the mean amplitude.

The **template SQI** is the mean correlation of the beats with their template. Peak, notch and diastolic points on the template are refined by parabolic interpolation.

**Derivatives** are Savitzky–Golay (70 ms, order 3), applied on a margin of one filter window around the template. Because the template grid is fixed, the filter is the same at every input rate. Without it, the window rounds to an odd number of samples (e.g. 76 ms at 250 Hz), and SDPTG b/a shifted from −1.12 to −1.59 on the same synthetic signal.

> **SDPTG ratios depend strongly on the derivative smoothing.** On the same signal, b/a is −0.65 with a 50 ms window, −1.12 with 70 ms and −1.64 with 80 ms. Compare b/a, c/a, d/a, e/a and the aging index only between analyses that use the same `sg_window_sec`, and report it.

**SDPTG points** are detected under constraints: the order a < b < c < d < e, and c and d must reach a prominence of at least 5% of a.

**The dicrotic point** is the notch (a minimum followed by a maximum) or, if there is no notch, the inflection point.

**Uncertainty:** each window value comes with `CI95 = 1.96 · SD / sqrt(n)` of the per-beat values.

**Beat-sequence features:**
* the DC level, the perfusion index and PAT are averaged over the window's beats;
* PAT uses R peaks from an optional ECG channel;
* per-beat values are available for every feature.

**Output grid:** the windows are resampled on a 1 Hz grid. Trailing windows use a zero-order hold, which keeps the grid causal; centered windows use linear interpolation.

Fixed-window baselines (e.g. 10 s and 30 s, 50% overlap) are computed with the same code, for comparison.

---

## Validation

`synthetic.py` generates PPG with a known ground truth:
* **Beat model:** three Gaussian waves per beat (systolic, tidal and diastolic reflection).
* **Scenario:** rest 120 s → acute stress (τ = 4 s) → slow recovery (τ = 25 s). During stress:
  * HR rises from 70 to 100 bpm and PAT falls from 250 to 190 ms;
  * the amplitude falls by 40%;
  * the tidal and reflected waves shrink.
* **Corruption:**
  * white noise (30 dB) and baseline wander;
  * respiratory AM and RSA, and IBI jitter;
  * motion-artefact bursts (about 1 per minute);
  * a synthetic ECG.

`benchmark.py` runs the pipeline and scores each estimator against the ground truth. The ground truth is the per-beat features of an ideal, noise-free regeneration of the same scenario, interpolated on the 1 Hz grid. The estimators are:
* the adaptive windows;
* fixed windows of 10, 30 and 60 s;
* a 10 s EMA of the per-beat values;
* the raw per-beat values.

The metrics:
* **NRMSE:** RMSE divided by the ground-truth SD in the stress scenario; lower is better.
* **Stationary NSD:** noise in the settled segments.
* **Latency:** delay of the 50% crossing after the stress onset.

Median over features, **stress scenario, 3 seeds** (`adaptive-ppg benchmark --seeds 1 2 3`):

| Estimator | NRMSE | Stationary NSD | Latency (s) |
| :--- | :---: | :---: | :---: |
| fixed 10 s | **0.374** | 0.158 | **8.2** |
| **adaptive (default)** | 0.414 | **0.116** | 9.3 |
| EMA 10 s | 0.459 | 0.261 | 11.2 |
| per beat | 0.637 | 0.506 | – |
| fixed 30 s | 0.647 | 0.383 | 24.7 |
| fixed 60 s | 0.971 | 0.641 | 40.5 |

**Rest scenario** (stationary, `--scenario rest`, seed 1), NRMSE: fixed 60 s 0.079, EMA 10 s 0.081, fixed 30 s 0.091, **adaptive 0.093**, fixed 10 s 0.116.

**Interpretation.** No fixed window is good in both regimes:
* 10 s windows track transients but are the noisiest at rest;
* 60 s windows are best at rest and worst under stress.

The adaptive windows are close to the best fixed window in each regime, and they need no prior choice of window. They have the lowest stationary noise during the stress protocol. They are **not** the winner of either regime, and across the two scenarios their mean NRMSE (0.254) is close to, but not below, that of fixed 10 s windows (0.245). Their advantage is in the stationary segments: lower noise than 10 s windows, at a latency cost of about 1 s.

**Fragile features.** Short (3-beat) windows make fragile fiducial points noisier. Sub-sample fiducials removed the largest case: ΔT_DVP now has an adaptive NRMSE of 0.46 against 0.59 for fixed 10 s (it was 0.86 against 0.39 with integer fiducials). The inflection-point area ratio has NRMSE > 1 with every estimator and should be interpreted with care. For fragile features, raise `T_crit` of the category.

### Sampling-rate independence

`tests/test_sampling.py` checks, among other things, that the same 500 Hz recording gives the same fixed-window features whether it is analysed natively or at 125 Hz. On the synthetic rest scenario, the stationary NSD at 125 Hz fell from 0.104 with integer fiducials to 0.036 with sub-sample fiducials and the fixed template grid (0.041 and 0.039 when analysed natively at 500 and 1000 Hz). The per-beat b/a median is now between −1.12 and −1.17 at every processing rate from 1000 Hz down to 100 Hz.

### Real recordings

Two excerpts of the Pulse Transit Time PPG Dataset (`Sample1.CSV`, `Sample2.CSV`, 500 Hz, infrared light at the distal phalanx, seated; see [Sample data](#sample-data)) and BIDMC record 01 were analysed with the default configuration. `tests/test_real_data.py` checks that the results stay within plausible ranges.

| | `Sample1.CSV` | `Sample2.CSV` | `bidmc_01` |
| :--- | :---: | :---: | :---: |
| rate (input → processing) | 500 → 125 Hz | 500 → 125 Hz | 125 Hz |
| polarity (derivative skewness) | inverted (−2.45) | inverted (−2.39) | upright |
| beats (valid) | 335 (100 %) | 277 (99.6 %) | 718 (98.2 %) |
| heart rate (median) | 71 bpm | 74 bpm | 91 bpm |
| median drift_z / beats with H > θ | 1.18 / 26 % | 1.06 / 21 % | 0.83 / 7 % |
| macro window, median [p10–p90] | 17 [4–33] beats | 18 [5–34] beats | 35 beats |
| derivatives window, median | 24 beats | 26 beats | 59 beats |
| template SQI, macro / derivatives | 0.985 / 0.958 | 0.989 / 0.982 | 0.993 |
| SDPTG b/a | −0.72 | −0.70 | −0.87 |

* Both PTT excerpts were stored inverted; the automatic polarity check detected it.
* drift_z stays close to 1 on all three recordings, as expected for stable morphology. The PTT excerpts have more beats with H > θ than BIDMC, so their windows are shorter. Both were recorded seated at rest, so activity does not explain the difference; the sensor (raw MAX30101 counts at the fingertip, against the processed pleth channel of a clinical monitor) and contact pressure are more likely causes. This has not been tested.
* On the 1 Hz grid, the SDPTG features are about half as rough with adaptive windows as with fixed 10 s windows (Sample1 b/a: 0.09 against 0.21 per-beat SD per second), with similar or smaller CI95. The notch time and the reflection index are rougher on Sample1 with adaptive windows, because they follow the shorter macro/time_volume windows.

Use `adaptive-ppg sweep` for sensitivity analyses. The defaults (`drift_z`, `fixed_scale=(1, 3)`, `θ = 0.3`, `N_past = 10`, `T_min = 3`) come from such a sweep, over 3 stress seeds and 1 rest seed.

---

## Installation

Python ≥ 3.10 is recommended; the core package supports 3.9. Step-by-step guides for users new to Python: [English](docs/INSTALL_en.md) · [Italiano](docs/INSTALL_it.md) · [日本語](docs/INSTALL_ja.md).

```bash
git clone https://github.com/badoglio/Adaptive-PPG-Time-Windows.git
cd Adaptive-PPG-Time-Windows
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -e ".[app,edf,yaml]"   # add "dev" for pytest/ruff, "biosppy" for the alternative detector
```

The core package needs only numpy, scipy, pandas and plotly. The extras add:
* `edf`: pyedflib, for EDF/BDF files;
* `app`: streamlit and openpyxl;
* `yaml`: YAML configuration files;
* `biosppy`: the BioSPPy detector (installs peakutils as well);
* `dev`: pytest and ruff.

## Usage

### Dashboard

```bash
adaptive-ppg dashboard             # or: streamlit run app/streamlit_app.py
```

`adaptive-ppg dashboard` works from a source checkout (editable install), where it finds `app/streamlit_app.py`; arguments after it are passed to Streamlit.

**Data sources:**
* the bundled `sample_data`;
* your own CSV/TSV/TXT/MAT/EDF/BDF file. The sampling rate is inferred from a time column or header, or entered manually. The CSV separator (`,` `;` tab), the decimal comma and a UTF-8 BOM are detected automatically. Leave the column empty to use the first non-time column;
* a synthetic demo with ground truth.

**Sidebar:**
* presets;
* variability metric, normalization, θ, k and N_past;
* window anchor;
* per-category `T_min`/`T_max`/`T_crit`;
* filter band, onset method and detector;
* resampling to a common processing rate, and the signal polarity (auto/normal/inverted);
* subject height;
* fixed-window baselines.

**Tabs:**
* signal and beats, with rejected beats and artefacts shaded;
* variability, H, `T_w` and `O_w`;
* features: adaptive with CI95 against fixed and per-beat values, plus the ground truth in the demo;
* ensemble templates;
* benchmark;
* export as ZIP of CSVs plus `metadata.json`, as an Excel workbook, or as the configuration JSON.

### Command line

```bash
adaptive-ppg run sample_data/bidmc_01_Signals.csv --signal-column PLETH --ecg-column II --out results/
adaptive-ppg run "data/*.edf" --channel Pleth --preset acute_stress --fixed 10 30 --templates --out results/
adaptive-ppg run my_recording.csv --signal-column PPG --fs 500 --out results/mine   # no time column: give fs
adaptive-ppg run sample_data/Sample1.CSV --target-fs 250 --out results/sample1   # --target-fs 0: native rate
adaptive-ppg init-config my_config.json --preset rest        # edit, then pass with --config
adaptive-ppg benchmark --seeds 1 2 3 --out bench/
adaptive-ppg sweep --param engine.theta 0.2 0.3 0.5 --param engine.N_past 6 10 --with-rest --out sweep.csv
```

Batch runs write one folder per recording and a `batch_summary.csv`. Every export embeds, in `metadata.json`:
* the full configuration;
* the package version and git revision;
* the input sampling rate and its source, and the processing rate;
* the detected polarity and the derivative skewness it was based on;
* the warnings.

### Python

```python
from adaptive_ppg.io import BioSignalLoader
from adaptive_ppg.config import preset
from adaptive_ppg.pipeline import run

sd = BioSignalLoader.auto_load("sample_data/bidmc_01_Signals.csv", signal_column="PLETH", ecg_column="II")
res = run(sd, preset("default"), fixed_window_secs=[30])
res.features.grid            # 1 Hz feature grid (adaptive)
res.window_table()           # every adaptive window with features and CI95
res.save("results/bidmc01", include_templates=True)
```

For real-time use, `adaptive_ppg.engine.StreamingAdaptiveEngine.push(...)` takes one beat at a time and returns the windows as they close.

### Configuration

All parameters are in one serializable `PipelineConfig`, with the sections `preprocessing`, `engine` (including `categories`) and `features`. Partial JSON/YAML files are merged over the defaults, unknown keys raise an error, and a category set to `null` is disabled.

Parameters related to the sampling rate:

| Parameter | Default | Meaning |
| :--- | :---: | :--- |
| `preprocessing.target_fs` | 125.0 | processing rate in Hz (≥ 20); `null` keeps the native rate |
| `preprocessing.polarity` | `"auto"` | `"auto"`, `"normal"` or `"inverted"` |
| `features.template_fs` | 500.0 | template grid in Hz (≥ 50); `null` uses the processing rate |
| `features.sg_window_sec` | 0.07 | Savitzky–Golay window of the derivatives (report it with SDPTG results) |

The presets:
* `default`;
* `rest`: longer windows, θ = 0.5;
* `acute_stress`: shorter windows, θ = 0.2;
* `wearable_64hz`: narrower morphology band (15 Hz), longer derivative smoothing (110 ms) and a larger `T_crit` for the derivatives;
* `legacy`: the original successive/calibration method.

---

## Package layout

```
src/adaptive_ppg/
  io.py             CSV/TSV/TXT (dialect sniffing), MAT (v5–v7.2), EDF/BDF loading; robust fs estimation
  preprocessing.py  resampling, polarity, filtering, Elgendi detection, sub-sample fiducials, SQI,
                    beat normalization
  engine.py         variability metrics, H, per-category schedules; batch and streaming engines
  features.py       ensemble templates, 33 features, PAT, fixed-window baselines, output grid
  synthetic.py      synthetic PPG/ECG with ground truth
  benchmark.py      ground truth, metrics, benchmark and parameter sweeps
  pipeline.py       run() and PipelineResult (tables, ZIP/Excel/CSV export, provenance)
  visualization.py  Plotly figures used by the dashboard
  config.py         PipelineConfig, validation, presets
  cli.py            command-line interface
  _typing.py        array type aliases (FloatArray, IntArray, BoolArray); the package ships py.typed
app/streamlit_app.py  Streamlit dashboard
sample_data/        example recordings and their licenses (sample_data/README.md)
tests/              pytest suite (run: pytest -q)
```

---

## Sample data

The recordings are third-party data under their own licenses; see [`sample_data/README.md`](sample_data/README.md).

* **`sample_data/Sample1.CSV`, `Sample2.CSV`**: excerpts (282 s and 226 s at 500 Hz) of the **Pulse Transit Time PPG Dataset** (PhysioNet, v1.1.0), released under the Open Database License (ODbL 1.0). Columns: `Time (s)`, `Pleth`, with `;` as separator and decimal commas.
  * Acquisition: infrared (IR) light at the distal phalanx, subject seated at rest.
  * `Pleth` holds raw MAX30101 counts, so the pulse is inverted; the automatic polarity check flips it.
  * About 23 % of the samples repeat the previous value exactly. The repeats are spread evenly and are mostly single samples, which suggests resynchronization to the 500 Hz output clock. They have no visible effect after resampling.
  * Citation: Mehrgardt, P., Khushi, M., Poon, S., & Withana, A. (2022). *Pulse Transit Time PPG Dataset* (version 1.1.0). PhysioNet. https://doi.org/10.13026/jpan-6n92 — and Goldberger, A. L. et al., "PhysioBank, PhysioToolkit, and PhysioNet", *Circulation* 101(23):e215–e220, 2000.
* **`sample_data/bidmc_01_Signals.csv`**: record 01 of the **BIDMC PPG and Respiration Dataset** (PhysioNet), 8 min at 125 Hz. Columns: `Time [s]`, `RESP`, `PLETH` (PPG), `V`, `AVR`, `II` (ECG).
  * The time column is rounded to 2–3 decimals. The loader therefore estimates fs by a least-squares fit of time against sample index; the median difference would give 100 Hz.
  * PAT on this record is about 0.53 s. That is longer than physiological PAT, and consistent with the known signal-processing delay of the MIMIC monitor's pleth channel. Treat absolute PAT from BIDMC/MIMIC with caution.
  * Citation: Pimentel, M. A. F., Charlton, P. H., & Clifton, D. A. (2016). *BIDMC PPG and Respiration Dataset* (v1.0.0). PhysioNet. https://doi.org/10.13026/C2GG69 — and Pimentel et al., "Toward a robust estimation of respiratory rate from pulse oximeters", *IEEE Trans. Biomed. Eng.* 64(8):1914–1923, 2017.

## Development

```bash
pip install -e ".[dev,edf,yaml,app]"
ruff check src tests app
pytest -q
```

GitHub Actions runs ruff and the test suite on Linux and Windows (`.github/workflows/ci.yml`). Changes between versions are listed in [`CHANGELOG.md`](CHANGELOG.md).

## License

The code is released under the **GNU General Public License v3.0 or later** (`GPL-3.0-or-later`, see [`LICENSE`](LICENSE)). The sample recordings keep their own licenses (ODC-By 1.0 for BIDMC, ODbL 1.0 for the Pulse Transit Time PPG Dataset), listed in [`sample_data/README.md`](sample_data/README.md).
