# Operative guide to the dashboard

[English](Operative_Guide_en.md) · [Italiano](Operative_Guide_it.md) · [日本語](Operative_Guide_ja.md)

This guide explains every control, plot and table of the web dashboard of **adaptive_ppg 0.2.0**. For each
item it says what it means, how to choose its value for a given study, and how to read the results in
physiological terms. To install and start the dashboard, see the [installation guide](INSTALL_en.md). The
method is summarised in the [README](../README.md).

* The dashboard is a research tool. It is not a medical device, and its outputs are not validated for
  diagnosis.
* The guide describes what the code of version 0.2.0 does. Where the README differs from the code, the guide
  follows the code and says so ([section 20](#20-known-issues)).
* GUI labels are written in **bold**, exactly as they appear on screen. Configuration keys (for the JSON
  file and the command line) are written as `code`.
* A **beat** is one pulse, from its foot (onset) to the foot of the next pulse. A **valid beat** is a beat
  that passed all the quality checks. **a.u.** means arbitrary units, that is, the units of the input
  signal.

---

## Contents

1. [How the analysis works](#1-how-the-analysis-works)
2. [Recommended workflow](#2-recommended-workflow)
3. [Sidebar: 1. Data](#3-sidebar-1-data)
4. [Sidebar: 2. Method](#4-sidebar-2-method)
5. [Sidebar: 3. Windows per category](#5-sidebar-3-windows-per-category)
6. [Sidebar: 4. Signal processing](#6-sidebar-4-signal-processing)
7. [Sidebar: 5. Fixed-window baselines](#7-sidebar-5-fixed-window-baselines)
8. [Page header, warnings and KPIs](#8-page-header-warnings-and-kpis)
9. [Tab: Signal & beats](#9-tab-signal--beats)
10. [Tab: Variability & windows](#10-tab-variability--windows)
11. [Tab: Features](#11-tab-features)
12. [Tab: Templates](#12-tab-templates)
13. [Tab: Benchmark](#13-tab-benchmark)
14. [Tab: Export](#14-tab-export)
15. [Feature reference](#15-feature-reference)
16. [Settings by use case](#16-settings-by-use-case)
17. [Quality-control checklist](#17-quality-control-checklist)
18. [What to report](#18-what-to-report)
19. [Reference values from the sample data](#19-reference-values-from-the-sample-data)
20. [Known issues](#20-known-issues)
21. [References](#21-references)

---

## 1. How the analysis works

The dashboard runs the same pipeline as the `adaptive-ppg run` command.

1. **Loading.** The PPG (and, optionally, an ECG) is read from the file. Missing values (NaN, Inf) are filled
   by linear interpolation. Recordings shorter than 10 s are refused.
2. **Saturation (clipping).** At the native rate, runs of at least 40 ms (and at least 2 samples) that stay
   within 1 % of the signal range from the global minimum or maximum are marked as clipped.
3. **Resampling** to a common processing rate (125 Hz by default), with a polyphase filter.
4. **Filtering** in two branches, both zero-phase (Butterworth of order 4 applied forward and backward to the
   signal minus its median):
   * a *detection* branch, 0.5–8 Hz, used only to find the systolic peaks;
   * a *morphology* branch, 0.5–30 Hz by default, used for everything else.
5. **Polarity.** If the derivative of the detection signal has negative skewness, the signal is flipped (see
   **Signal polarity** in [section 6](#6-sidebar-4-signal-processing)).
6. **Peak detection** (Elgendi et al. 2013 by default) and **beat segmentation**. Each beat runs from its
   onset to the onset of the next beat, so the last detected peak, which has no following onset, is dropped.
7. **Beat quality.** Each beat is compared with a reference built from the preceding 30 beats and is marked
   valid or invalid.
8. **Morphological variability.** Each beat is normalised: the straight line from its onset to its end is
   subtracted, the peak is scaled to 1, and the beat is resampled to 128 points with the systolic peak at
   30 % of the beat. A variability metric is computed on the last valid normalised beats and mapped to an
   index **H** between 0 (stable morphology) and 1 (changing morphology).
9. **Adaptive windows.** For each feature category, H sets the window length (in valid beats) and the
   overlap between consecutive windows: long windows when the morphology is stable, short windows when it
   changes.
10. **Ensemble templates and features.** The beats of each window are aligned and averaged into a template,
    and the features are measured on the template. The same features are also measured on every single beat
    and on fixed-duration windows, for comparison.

### Equations

Default variability metric, `drift_z`, computed on the normalised valid beats (`N` = **N_past**):

```
drift_k   = RMS over the 128 samples of [ mean(beats k−N+1 … k) − mean(beats k−2N+1 … k−N) ]
drift_z_k = drift_k / sqrt( (2/N) · σ²_k )
σ²_k      = within-block variance of the 2N beats, pooled over the two blocks (2N − 2 degrees of freedom)
            and averaged over the 128 samples
```

If the morphology is stationary and the beat-to-beat fluctuations are independent, the expected value of
`drift²` is `2σ²/N`, so `drift_z ≈ 1` whatever the noise level. A step change of RMS size Δ placed between the
two blocks gives approximately `drift_z ≈ sqrt(1 + N·Δ²/(2σ²))`. The metric therefore measures the change in
units of the beat-to-beat scatter, and its sensitivity grows with √N.

Normalisation, window length and overlap (the same H drives all categories; each category `c` has its own
`T_min`, `T_max` and `T_crit`):

```
H      = clip( (v − lo) / (hi − lo), 0, 1 )           v = variability; lo, hi set by H normalization
s(H)   = 1 / (1 + exp( k · (H − θ) ))                  k = Sigmoid slope, θ = threshold
O      = O_max − (O_max − O_min) · s(H)                O_min = 0.25, O_max = 0.85
T_max  = max( T_max_sec / IBI_local , T_min )          in beats; IBI_local = median of the last 10 valid beat durations
T*     = T_min + (T_max − T_min) · s(H)
L_k    = min( T*_k , L_(k−1) + Δ )                     growth limited to Δ beats per beat (Δ = 1); shrinking is immediate
T_k    = max( round(L_k), T_crit )                     window length in valid beats
step_k = max( 1, round( T_k · (1 − O_k) ) )            beats between two consecutive window anchors
```

Values of the shape factor with the defaults (k = 10, θ = 0.3):

| H | 0 | 0.1 | 0.2 | 0.3 | 0.4 | 0.5 | 1 |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| s(H) | 0.953 | 0.881 | 0.731 | 0.500 | 0.269 | 0.119 | 0.001 |

Worked example (`Sample1.CSV`, 71 bpm, IBI ≈ 0.84 s, macro category). `T_max` = 30 s / 0.84 s ≈ 35.5 beats.
With H = 0, T ≈ 3 + 32.5 × 0.953 ≈ 34 beats and O ≈ 0.28, so a new window starts about every 25 beats. With
H = 1 the window shrinks at once to 3 beats with O = 0.85, and a new window starts at every beat. After the
change, the window grows back by at most one beat per beat (about 31 beats to go from 3 back to 34).

### What H can and cannot see

Because each beat is normalised before the variability is computed:

* H does **not** respond to a pure change of pulse amplitude (all the samples of the beat scaled by the same
  factor).
* H does **not** respond to a uniform stretching of the beat in time, and the piecewise resampling (peak always
  at 30 %) also removes a change of the rise-to-decay proportion. Changes of crest time and duty cycle alone
  therefore do not shrink the windows.
* H **does** respond to changes of shape within the upstroke and within the decay: for example, a change of
  the height or position of the late-systolic/diastolic wave and of the dicrotic notch relative to the decay,
  or of the curvature of the upstroke.
* A change of heart rate can still change H indirectly, because the diastolic part of the pulse does not
  scale with the cycle length.

Practical consequence: when the response of interest is mainly a change of amplitude or of timing (for
example the fall in pulse amplitude with sympathetic vasoconstriction), H may stay low and the windows stay
long. In such protocols, set a shorter **T_max (s)** for the categories of interest
([section 16](#16-settings-by-use-case)).

---

## 2. Recommended workflow

1. **Load** the recording ([section 3](#3-sidebar-1-data)). Read the info line under the title: channel,
   sampling rate, duration, ECG, polarity.
2. Read the **warnings** (yellow boxes). If a warning reports a sampling-rate mismatch or a wrong column, fix
   it before going on.
3. In **Signal & beats**, zoom on a few beats. The filtered pulse must point upwards with a steep upstroke.
   The green dots must sit on the systolic peaks and the orange triangles at the feet. Red spans are rejected
   beats.
4. Check the four **KPIs**: valid beats, share of beats with H > θ, median window and template SQI.
5. In **Variability & windows**, check that H rises where you expect a change (for example at the start of a
   task) and stays low at rest.
6. In **Features**, look at the features of interest with their CI95 and compare them with the fixed-window
   baselines.
7. In **Templates**, check that the ensemble pulses have a plausible shape.
8. In **Export**, download the ZIP and the configuration JSON. Write down the fixed windows and the preset:
   they are not saved in the JSON.

Every change of a control re-runs the analysis. Results are cached, so going back to an earlier setting is
fast.

---

## 3. Sidebar: 1. Data

**Source** selects where the signal comes from:

* **Synthetic demo (known ground truth)**: a simulated recording whose true features are known. It is needed
  by the **Benchmark** tab.
* **Sample data**: the recordings in the `sample_data` folder. This is the default when the folder exists.
* **Upload a file**: your own recording.

### Synthetic demo

* **Scenario**:
  * `stress`: 120 s of rest, then 120 s of acute stress reached quickly (time constant 4 s), then 120 s of
    slow recovery (time constant 25 s);
  * `rest`: 300 s, stationary.
* **Sampling rate (Hz)**: 25–1000 Hz, default 125. Use it to see how the acquisition rate affects the
  features, in particular the derivative (SDPTG) features.
* **Noise seed**: 0–10000, default 1. It changes the random noise; the underlying physiology stays the same.

Each simulated beat is the sum of three Gaussian waves: systolic, late-systolic and diastolic. The generator
changes these parameters between rest and stress (amplitudes relative to the systolic wave, delays after the
systolic wave):

| Parameter | Rest | Stress |
| :--- | :---: | :---: |
| Heart rate (bpm) | 70 | 100 |
| Systolic wave: time / SD (s) | 0.16 / 0.055 | 0.13 / 0.045 |
| Late-systolic wave: amplitude / delay (s) / SD (s) | 0.35 / 0.09 / 0.045 | 0.12 / 0.08 / 0.045 |
| Diastolic wave: amplitude / delay (s) / SD (s) | 0.45 / 0.24 / 0.09 | 0.22 / 0.19 / 0.08 |
| Pulse arrival time (s) | 0.25 | 0.19 |
| Pulse amplitude | 1 | 0.6 |

Noise and confounders: white noise at 30 dB SNR; baseline wander (0.3 × pulse amplitude, 0.05–0.15 Hz);
respiration at 0.25 Hz with 8 % amplitude modulation and 3 bpm of respiratory sinus arrhythmia; 1 % random
beat-to-beat interval jitter; motion artefacts about once per minute (1–4 s, 2 × pulse amplitude). A
synthetic ECG is included, so PAT is available. The DC level is 2.0, so the synthetic perfusion index (about
50 %) is not physiological.

### Sample data

**File** lists the files in `sample_data`, sorted by name. `Sample1.CSV` is the default.

| File | Source | Content |
| :--- | :--- | :--- |
| `Sample1.CSV`, `Sample2.CSV` | Pulse Transit Time PPG Dataset 1.1.0 (Mehrgardt et al. 2022), ODbL licence | 500 Hz. Infrared PPG at the distal phalanx, subject seated at rest. Raw MAX30101 counts, which are inverted (more blood, less light). Separator `;`, decimal comma. About 23 % of consecutive samples are repeated values. |
| `bidmc_01_Signals.csv` | BIDMC dataset (Pimentel et al. 2017), derived from MIMIC-II (Goldberger et al. 2000), ODC-By licence | 125 Hz, 8 min, critically ill patient. Columns `Time [s]`, `RESP`, `PLETH`, `V`, `AVR`, `II`. The dashboard uses `PLETH` as the PPG and lead `II` as the ECG. |

### Upload a file

**PPG recording** accepts CSV, TXT, TSV, MAT, EDF and BDF files up to 200 MB.

* **CSV, TXT, TSV.** The separator (tab, `;`, `,`, `|`, otherwise spaces) and a decimal comma are detected
  automatically. If the first field of the first row is a number, the file is read without a header and the
  columns are named `col_0`, `col_1`, … A time column is recognised by its name: `time`, `t`, `sec`, `secs`,
  `seconds`, `tempo` or `timestamp` (any case), followed by a word boundary or `_`, as in `Time [s]` or
  `time_ms`. Milliseconds are recognised from `[ms]`, `(ms)`, `_ms`, `msec` or `millisec` in the name.
* **MAT.** A name selects the variable that holds the PPG. A number selects the channel of the first numeric
  variable (for 2-D arrays). A variable called `time` or `t` (or another time name), otherwise one called
  `fs`, `sampling_rate`, `srate` or `sample_rate`, gives the sampling rate. The ECG field is ignored for MAT
  files.
* **EDF, BDF.** The sampling rate is read from the header. The ECG channel must have the same rate as the PPG.
  Reading needs the `pyedflib` package (install option `edf`) or, failing that, `mne`.

**PPG column / channel** — the column name or a 0-based index. If blank, the first non-time column (or
channel 0) is used. Names are matched exactly first, then without regard to case; a number is taken as an
index over all the columns, the time column included. For `bidmc_01_Signals.csv` the default is `PLETH`.

**ECG column / channel (optional, enables PAT)** — if given, R peaks are detected on it and the pulse arrival
time is computed ([section 15](#15-feature-reference)). The R-peak detector is simple: 5–20 Hz band-pass,
squared slope smoothed over 80 ms, one threshold at 30 % of its 99th percentile for the whole recording,
peaks at least 0.3 s apart. It is adequate for clean monitoring leads such as BIDMC lead II. For ambulatory
or noisy ECG, check the PAT values carefully. For `bidmc_01_Signals.csv` the default is `II`.

**Infer the sampling rate from the file** — when ticked, the rate comes from the time column or, for EDF/BDF,
from the header (always used for EDF/BDF). The rate is fitted by least squares on the whole time column, so
time stamps rounded to 2–3 decimals still give the right rate. When unticked, type the rate in
**Sampling rate (Hz)** (1–10000, default 125).

> **Tip.** A wrong sampling rate scales every time feature and the heart rate. If the rate you type differs
> by more than 1 % from the one implied by the time column, a warning gives the scaling factor. Trust the
> time column unless you know it is wrong.

---

## 4. Sidebar: 2. Method

### Preset

**Preset** loads a set of starting values. The controls below show them and can override them.

| Preset | Changes from `default` | Intended use |
| :--- | :--- | :--- |
| `default` | — | General use |
| `rest` | θ = 0.5; growth limit Δ = 2 beats per beat; **T_max (s)** 60 / 60 / 90 (macro / time_volume / derivatives) | Resting recordings: long windows that react only to large changes |
| `acute_stress` | θ = 0.2; Δ = 0.5; **T_max (s)** 20 / 20 / 40 | Protocols with fast autonomic transients |
| `wearable_64hz` | **Morphology low-pass (Hz)** = 15; derivative smoothing 110 ms; derivatives **T_crit** = 15 | Wrist wearables sampled at about 64 Hz |
| `legacy` | **Variability metric** = successive; **H normalization** = calibration; fixed scale 0.02–0.15; θ = 0.5 | The original method, for comparison |

> **Important.** The growth limit Δ and the derivative smoothing window have no control in the GUI, so the
> preset value always applies. The controls in **3. Windows per category** do **not** follow a change of
> preset after the first run: type the preset values of **T_max (s)** and **T_crit** by hand
> ([section 20](#20-known-issues)). The other controls do follow the preset.

### Variability metric

| Metric | Definition | Behaviour |
| :--- | :--- | :--- |
| `drift_z` (default) | `drift` divided by its expected value under stationarity | ≈ 1 on a stationary signal whatever the noise level; designed for the fixed scale 1–3 |
| `drift` | RMS difference between the mean of the last N beats and the mean of the N beats before | Measures sustained change, in units of the normalised pulse (fractions of the peak); depends on the noise level |
| `successive` | Trailing mean over N beats of the RMS difference between consecutive beats | Mostly measures beat-to-beat noise; used by the `legacy` preset |
| `dispersion` | Mean RMS deviation of the last N beats from their mean | Scatter within the block, not change over time |

`drift` and `drift_z` are undefined until 2N valid beats are available (H is held at 0 meanwhile).
`successive` and `dispersion` are defined from the second valid beat.

> **Tip.** `drift`, `successive` and `dispersion` are expressed as fractions of the peak (typically a few
> hundredths). With **H normalization** = fixed and the default scale 1–3, H stays at 0 and the windows never
> shrink. Use them with calibration, rolling or global normalisation, or with a fixed scale of the same order
> (the `legacy` preset uses 0.02 and 0.15).

### H normalization

* **fixed** (default): H = 0 at **H = 0 at** and H = 1 at **H = 1 at**. The scale is the same for every
  recording, so H can be compared between recordings and subjects. There is no warm-up other than the one of
  the metric itself (the first 2N − 1 valid beats for `drift_z`).
* **calibration**: lo and hi are the 5th and 95th percentiles of the metric during the first 60 s. This
  assumes that the first 60 s are a representative baseline. H is held at 0 during these 60 s (warm-up). If
  fewer than 3 values are available, the fixed scale is used.
* **rolling**: lo and hi are the 5th and 95th percentiles of the last 120 valid beats. The scale follows
  slow changes. H is held at 0 during the first 60 s.
* **global**: percentiles of the whole recording. This uses future data, so it is offline only.

> **Tip.** The relative normalisations (calibration, rolling, global) stretch whatever variability the
> recording has onto 0–1. By construction, about 5 % of the reference beats reach H = 1 and a sizeable share
> exceeds θ, even when nothing changes. Prefer **fixed** with `drift_z` unless you have a reason to rescale
> each recording.

**H = 0 at** and **H = 1 at** (shown only with fixed normalisation; three decimals; **H = 1 at** is kept
above **H = 0 at**). With `drift_z`, the default 1 means "no more change than expected from the beat-to-beat
scatter", and 3 means a change about three times larger. On the resting sample recordings the median
`drift_z` is 0.83–1.18 ([section 19](#19-reference-values-from-the-sample-data)). If a resting recording shows
many contractions, look at its `variability_used` values in the **Beat quality table** and consider raising
**H = 0 at** slightly above their median.

### θ and slope

**θ (H at which windows start shrinking)** — 0.05–0.95, default 0.3. At H = θ the window length is halfway
between `T_max` and `T_min`. The corresponding value of the metric is `lo + θ·(hi − lo)`: with `drift_z` and
the scale 1–3, θ = 0.3 corresponds to `drift_z` = 1.6. A lower θ makes the windows shrink for smaller changes
(more reactive, but noisier features). A higher θ makes them shrink only for large changes.

**Sigmoid slope k** — 2–30, default 10. It sets how sharp the transition between long and short windows is.
With k = 10 and θ = 0.3, s(0) = 0.953, so the windows are almost at `T_max` when nothing changes. With a low
slope the transition is gradual, but the windows never reach `T_max` (k = 2 gives s(0) = 0.646). With
k ≥ 20 the behaviour is nearly binary (long or short).

### N_past

**N_past (beats per drift block)** — 4–30, default 10. The number of beats in each of the two blocks compared
by `drift` and `drift_z`. It is also the averaging length of `successive`, of `dispersion` and of the IBI
regressor. The choice is a trade-off:

* **Larger N**: more beats in each mean, so smaller sustained changes are detected (the response of
  `drift_z` grows as √N) and periodic modulations average out better.
* **Smaller N**: faster response and shorter warm-up. The metric needs 2N valid beats: with N = 10 the
  warm-up is 19 valid beats, about 16 s at 70 bpm. For a step change, the trailing `drift` peaks about N
  beats after the step and stays raised for about 2N beats.

Physiological guidance: each block should cover at least about **two respiratory cycles**, so that the
respiratory modulation of the pulse averages out in both blocks. In beats, `N ≥ 2 × heart rate / breathing
rate`.

* At 15 breaths/min and 70 bpm a breathing cycle lasts about 4.7 beats, so N = 10 covers about two cycles.
* With slow or paced breathing (6 breaths/min at 60 bpm, 10 beats per cycle), use N ≥ 20.
* Mayer waves (about 0.1 Hz; Julien 2006) last about 12 beats at 70 bpm and are only partly averaged out with
  the default N.

### Window anchor

**Window anchor** — `trailing` (default) or `centered`.

* **trailing**: each window ends at its anchor beat. Every value at time t depends only on beats up to t
  (causal, as in real time). A trailing mean lags a trend by about half its length, so feature values lag
  the physiology by about T/2 (about 15 s for a 30 s window at rest).
* **centered**: each window is centred on its anchor beat. The feature value is not delayed relative to the
  anchor, which is better for offline event-related analysis. H and the local IBI are read N_past//2 beats
  later, but H still reacts about N/2 beats after a change (instead of about N beats). Centred windows use
  future data, so they are offline only.

### Regress out IBI-driven variability

**Regress out IBI-driven variability** — when ticked, the part of the metric that follows the fluctuations of
the beat intervals is removed. The regressor is the trailing mean over N beats of `|ΔIBI| / IBI`. A causal
expanding least-squares fit starts after 30 pairs, and the metric becomes
`v − slope · (x − mean(x))`. The slope is saved in the metadata (`ibi_regression_slope`).

Use it when beat-interval fluctuations, for example a strong respiratory sinus arrhythmia or irregular
rhythm, make H follow the intervals rather than the shape. Leave it off when a real change of shape is
expected together with a change of the beat intervals, because part of the real change would be removed.

---

## 5. Sidebar: 3. Windows per category

The features are grouped into three categories. Each category has its own window limits:

| Category | Features | **T_min (beats)** | **T_max (s)** | **T_crit** |
| :--- | :--- | :---: | :---: | :---: |
| **macro** (17) | pulse amplitude, pulse widths, areas, area ratios, inflection amplitude, reflection index, slopes, DC level, perfusion index | 3 | 30 | 3 |
| **time_volume** (10) | pulse duration, heart rate, crest time, decay time, duty cycle, notch time, ΔT_DVP, stiffness index, PAT to foot and to peak | 3 | 30 | 3 |
| **derivatives** (6) | SDPTG b/a, c/a, d/a, e/a, aging index, (b − e)/a | 5 | 60 | 10 |

* **T_min (beats)** (1–60): window length, in valid beats, when the morphology changes (H near 1). It sets
  the time resolution during transients: 3 beats are about 2.5 s at 70 bpm.
* **T_max (s)** (5–300, step 5): window length, in seconds, when the morphology is stable. It is converted to
  beats with the local IBI (median of the last 10 valid beats), so at 100 bpm a 30 s window holds about 50
  beats. It sets the smoothing at rest and the lag of trailing windows (about T_max/2).
* **T_crit** (1–60): the minimum number of valid beats in a window. Windows with fewer beats are not produced,
  and the window length never falls below T_crit.

Guidance:

* The derivative features have longer windows because the second derivative amplifies noise.
* For fragile fiducial points (dicrotic notch, inflection-point area ratio, SDPTG c and d waves), raise
  **T_crit** of the category.
* Choose **T_max (s)** shorter than the shortest change you want to follow without H noticing it. This
  matters for changes of amplitude or timing alone, which H does not see
  ([section 1](#what-h-can-and-cannot-see)).
* At the start of the recording a trailing window holds all the valid beats available so far, so `n_valid`
  can be smaller than `T_beats`.
* The growth limit Δ (`max_expand_per_beat`) has no GUI control: 1 by default, 2 with `rest`, 0.5 with
  `acute_stress`. To set another value, use a configuration file with the command line.

---

## 6. Sidebar: 4. Signal processing

**Resample to a common rate** (default on) and **Processing rate (Hz)** (50–1000, step 25, default 125) — the
whole analysis runs at this rate. Templates are always built on a 500 Hz grid (cubic B-spline interpolation),
so the derivative filter is the same at every rate. Upsampling does not add bandwidth: the filter cut-offs
are limited to 45 % of the lower of the acquisition and processing rates. When unticked, the native rate is
used (slower for 500–1000 Hz recordings).

> **Tip.** The README reports, on `Sample1.CSV` and `Sample2.CSV`, a median difference from the native
> 500 Hz analysis of 0.02 per-beat SD at 250 Hz, 0.10 SD at 125 Hz, 0.14 SD at 100 Hz and 0.3 SD at 64 Hz.
> 125 Hz is adequate for most features. For derivative features, use 250 Hz or more when the acquisition
> rate allows it. Rates below 100 Hz are not recommended for derivative features.

**Signal polarity** — `auto` (default), `normal` or `inverted`. In both transmission and reflectance PPG,
more blood absorbs more light, so the detected light falls during systole. Raw photodetector counts are
therefore upside down with respect to blood volume (as in `Sample1.CSV` and `Sample2.CSV`). Clinical monitors
usually show the pleth already inverted, volume up. `auto` flips the signal when the derivative of the
detection signal has negative skewness: in a correctly oriented pulse the upstroke is shorter and steeper than
the decay, so the derivative is positively skewed. Choose `normal` or `inverted` when you know the
convention, or when `auto` fails (damped pulses with nearly symmetric rise and fall, very noisy recordings).

**Morphology high-pass (Hz)** — 0.05–5, default 0.5. It removes the DC, the baseline wander and the
respiratory baseline. Because the filter is applied forward and backward, the attenuation is −6 dB at the
cut-off.

| Frequency (Hz) | 0.3 | 0.5 | 0.667 (40 bpm) | 0.75 | 1 |
| :--- | :---: | :---: | :---: | :---: | :---: |
| Attenuation with high-pass 0.5 Hz | −36 dB | −6.0 dB | −0.76 dB | −0.29 dB | −0.02 dB |

* At 0.5 Hz, a pulse rate of 40 bpm loses 0.76 dB at its fundamental and 30 bpm loses 6 dB. For bradycardia
  (below about 45 bpm), lower the high-pass to 0.3 Hz (−0.13 dB at 0.5 Hz).
* A high-pass close to the pulse rate changes the pulse shape. Avoid values above about 0.7 Hz.
* Lower values keep more baseline wander. The linear onset-to-end baseline removed from each beat compensates
  for part of it.

**Morphology low-pass (Hz)** — 5–60, default 30 (−0.33 dB at 22.5 Hz, −6 dB at 30 Hz). The second
derivative amplifies high frequencies, so the SDPTG waves depend on this value together with the derivative
smoothing (70 ms) and the processing rate. Use 15 Hz for low-rate wearables (preset `wearable_64hz`; −0.56 dB
at 11.25 Hz). Lower values smooth the upstroke and change the SDPTG ratios. Keep the same value within a
study and report it. If the cut-off is above 45 % of the rate, it is lowered automatically and a warning is
shown.

**Bypass the morphology filter** — the morphology branch becomes the signal minus its median. The detection
branch is still filtered. Use it only for inputs that are already band-limited and free of baseline wander,
or to see the effect of the filter. Without the filter, baseline wander enters the templates (partly removed
by the per-beat linear baseline) and the derivative features are noisier.

**Onset method** — `tangent` (default) or `minimum`.

* **tangent**: the intersection of the tangent at the point of maximum upstroke slope with the horizontal
  line through the minimum before the upstroke (sub-sample precision). The intersecting-tangent foot is less
  sensitive to the shape of the foot and to noise, and is a common choice for pulse transit and arrival
  times (Mukkamala et al. 2015).
* **minimum**: the minimum between consecutive peaks, refined by a parabolic fit.

The tangent foot usually lies after the minimum. With `tangent`, crest time is shorter and `pat_foot` is
longer than with `minimum`. Do not mix the two methods within a study.

**Min. correlation with the reference beat** — 0.5–0.99, default 0.86. A beat passes `corr_ok` if the
Pearson correlation between its normalised shape and the reference (median of the previous 30 beats) is at
least this value. 0.86 follows the PPG template-matching threshold of Orphanidou et al. (2015), who applied
it to the mean correlation over 10 s segments; here it is applied to each beat. Li & Clifford (2012) discuss
other template-matching quality indices for pulsatile signals.

* A higher threshold rejects more atypical beats (artefacts, ectopic beats). It can also reject beats during
  a real, fast change of shape, which biases the analysis towards stationarity.
* A lower threshold keeps more beats but lets artefacts into the templates.

The other quality checks are fixed: pulse rate 30–220 bpm; beat duration within ±30 % of the reference
duration; amplitude between 1/3 and 3 times the reference amplitude; no clipping; finite values. For the
first 5 beats the reference is the median of the first 30 beats.

**Peak detector** — `elgendi` (default) or `biosppy`. `elgendi` is the two-moving-average detector of
Elgendi et al. (2013), with windows of 111 ms and 667 ms and an offset of 0.02, applied to the 0.5–8 Hz
detection branch. Its minimum peak distance of 0.3 s limits it to 200 bpm. `biosppy` runs the BioSPPy
implementation of the same algorithm on the same branch; it needs the optional `biosppy` package (and
`peakutils`). With both detectors the peaks are refined on the morphology signal within ±50 ms, with
parabolic interpolation. The **Benchmark** tab always uses `elgendi`.

**Subject height (m, 0 = unknown)** — enables the stiffness index, `SI = height / ΔT_DVP`
(Millasseau et al. 2002). Height is used as a proxy for the arterial path length.

---

## 7. Sidebar: 5. Fixed-window baselines

**Fixed windows (s, comma-separated)** — default `10, 30`. Each value creates a schedule of fixed-duration
windows, computed for all 33 features. Commas and semicolons are both accepted. An empty field means no fixed
windows. Text that cannot be read shows "Use numbers separated by commas." and 30 s is used.

Fixed windows use the same **Window anchor** as the adaptive windows:

* anchors every `W · (1 − overlap)` seconds, starting at the first peak plus W (trailing) or plus W/2
  (centred); the last peak is always an anchor;
* a window holds the valid beats whose peak lies in `(t_anchor − W, t_anchor]` (trailing) or within
  `t_anchor ± W/2` (centred);
* windows with fewer than 3 valid beats are skipped.

**Overlap of fixed windows** — 0–0.9, default 0.5. A higher overlap gives more frequent values but the same
smoothing.

Use the fixed windows as references. Short fixed windows (10 s) follow fast changes but are noisy. Long ones
(30–60 s) are smooth but lag and blur transients. The adaptive windows aim to combine the two.

---

## 8. Page header, warnings and KPIs

### Header and info line

The title "Adaptive-window PPG morphology" is followed by the package version. The info line gives the file
name, the PPG channel, the acquisition rate and, in brackets, the processing rate, the duration, and the notes
", with ECG (PAT enabled)" and ", inverted" when they apply.

### Warnings and errors

| Message | Meaning | What to do |
| :--- | :--- | :--- |
| "Sampling rate X Hz is low: derivative features (SDPTG) will be unreliable." | Acquisition below 50 Hz | Do not interpret SDPTG features |
| "Detection / Morphology high cutoff X Hz is too close to Nyquist …; using Y Hz." | A low-pass cut-off was above 45 % of the rate and was lowered | Report the cut-off actually used (`effective_highcut_hz` in the metadata) |
| "Signal inverted automatically (derivative skewness S < 0)." | `auto` polarity flipped the signal | Check in **Signal & beats** that the pulse now points upwards |
| "Only X% of beats passed quality checks." | Less than 50 % of valid beats | Check the column, the rate, the polarity and the signal quality |
| "Sampling rate set to X Hz but the time column implies Y Hz. …" | Typed rate and time column disagree by more than 1 % | Use the rate from the file, unless the time column is wrong |
| "No beats detected." | No peak found | Wrong column, flat signal or wrong rate |
| "Fewer than 20 valid beats: check the column, the sampling rate and the signal polarity." | Too few beats to analyse (error) | As above |
| "Could not analyse the recording: …" | The pipeline stopped, for example "Signal too short (x s); at least 10 s are required." | Read the message |
| "Invalid configuration: …" | Inconsistent settings, for example high-pass ≥ low-pass | Fix the setting named in the message |

### KPIs

* **Valid beats** — valid beats / detected beats, with the percentage. The sample recordings give 98–100 %.
  Below about 90 %, look at the rejected spans before interpreting the features.
* **Beats with H > θ** — the share of beats, after the warm-up, in which H exceeds θ, that is, in which the
  windows are shorter than halfway between `T_max` and `T_min`. Invalid beats are counted with the H of the
  last valid beat. With the default settings, the resting sample recordings give 7–26 % and the synthetic
  stress recording (seed 1) gives 14.8 %. A high value at rest suggests real non-stationarity (posture,
  breathing pattern, slow vasomotor waves, sensor contact) or a θ or **H = 0 at** that is too low.
* **Median macro window** — the median macro window length over all the beats (warm-up included), in beats,
  and in seconds (median length × median local IBI). Compare it with **T_max (s)**.
* **Median template SQI** — the median, over the macro windows, of the mean Pearson correlation between each
  beat and the window template. Clean recordings give values close to 1 (0.985–0.993 on the sample data).
  Each beat contributes to the template it is compared with, so the SQI is biased upwards in short windows.

---

## 9. Tab: Signal & beats

**Show the raw signal** (default on) adds a second panel with the raw signal.

The plot shows:

* **Filtered (morphology branch)**: the signal used for all the features, after the polarity correction. The
  title gives the morphology band, the detection band and the processing rate ("PPG — unfiltered" when the
  filter is bypassed).
* **Raw**: the resampled signal in its original polarity. It is **not** flipped, so with raw counts it points
  downwards.
* Green dots: systolic peaks of valid beats. Red dots: peaks of rejected beats. Orange triangles: onsets.
* Red shaded spans: rejected beats. Purple spans ("artefact"): simulated artefacts (synthetic data only).
* Lines with more than 20 000 points are decimated for display (min–max). The markers are exact.

### Beat quality table

One row per detected beat.

| Column | Meaning |
| :--- | :--- |
| `onset_idx`, `peak_idx`, `end_idx` | Onset, systolic peak and end (= next onset), as sample indices at the processing rate |
| `t_onset`, `t_peak` | Onset and peak times (s) |
| `duration_s` | Beat duration, onset to next onset (s): the pulse interval |
| `bpm` | 60 / `duration_s` |
| `amplitude` | Morphology signal at the peak minus at the onset (a.u.) |
| `template_corr` | Pearson correlation of the normalised beat with the reference beat |
| `amplitude_ratio` | Amplitude divided by the reference amplitude |
| `finite_ok` | No missing values |
| `rate_ok` | 30–220 bpm |
| `ibi_ok` | Duration within ±30 % of the reference duration |
| `corr_ok` | `template_corr` ≥ **Min. correlation with the reference beat** |
| `amp_ok` | Amplitude positive and `amplitude_ratio` between 1/3 and 3 |
| `clip_ok` | No clipped samples in the beat |
| `valid` | All the checks passed |
| `variability_raw` | Variability metric before the IBI regression |
| `variability_used` | Metric used for H (equal to `variability_raw` when the regression is off) |
| `H` | Normalised variability (in centred mode, the value read N_past//2 beats later) |
| `overlap` | Overlap O |
| `warmup` | True while H is held at 0 |
| `local_ibi` | Median of the last 10 valid beat durations (s) |
| `T_beats_macro`, `T_beats_time_volume`, `T_beats_derivatives` | Window length of each category at this beat (valid beats) |

Invalid beats carry the values of the last valid beat in the variability columns.

Reading the table:

* `ibi_ok` false: usually an ectopic beat, the beat after it, or a missed or extra peak.
* `corr_ok` false: an artefact or a beat of atypical shape.
* `amp_ok` false: motion, a sudden change of perfusion or of sensor contact.
* `clip_ok` false: sensor or amplifier saturation.
* `rate_ok` false: almost always a detection error.
* `duration_s` of the valid beats is a pulse-interval series. Pulse rate variability is not the same as heart
  rate variability, especially during stress, movement and low perfusion (Schäfer & Vagedes 2013).

---

## 10. Tab: Variability & windows

Four rows share the time axis (beat peak times):

1. **Morphological variability (metric)**: `variability_used`, with a dotted line at the H = 0 level and a
   dashed line at the H = 1 level.
2. **Normalized variability H**: H, with a dashed line at θ. Warm-up spans are shaded in rows 1 and 2.
3. **Window length T_w (valid beats)**: the window length of each category (step lines), with markers at the
   beats where a window was produced.
4. **Overlap O_w**: the overlap.

How to read it:

* At rest, `drift_z` fluctuates around 1 and H stays near 0. The windows sit at `T_max` and the markers are
  sparse.
* A sustained change of shape raises `drift_z` for about 2N beats, with a peak about N beats after the change
  (trailing). The window length drops at once and then grows back by one beat per beat (with the default
  growth limit).
* Events that can change the pulse shape include changes of posture, the start and end of mental or physical
  tasks, vasoactive stimuli (cold, pain), changes of sensor contact pressure, and motion artefacts that passed
  the quality checks. Compare the H peaks with the protocol markers and with the rejected spans.
* A periodic H at the breathing rate means that N_past is too short for the breathing pattern
  ([section 4](#n_past)).

---

## 11. Tab: Features

### Controls

* **Features**: the features to plot, shown as "name [unit]". Defaults: heart_rate, peak_amplitude,
  crest_time, sdptg_b_a. Only features with at least one finite value are listed: the stiffness index needs
  **Subject height**, PAT needs an ECG.
* **Show per-beat values**: grey dots, one per valid beat.
* **Fixed-window baselines**: the fixed-window schedules to overlay (all by default).

### Plot

Each feature is plotted in its own panel, titled "feature [unit] — category: cat".

* **adaptive (category)**: the template feature of each adaptive window, at `t_anchor` (the peak time of the
  anchor beat). With trailing windows the line is drawn as steps (the value holds until the next window);
  with centred windows the points are joined by straight lines. Hovering shows `n_valid` and `T_w`.
* **adaptive CI95**: a band of ± `1.96 · SD / √n` of the per-beat values of the window.
* **fixed_Ws**: the template feature of each fixed window.
* **ground truth** (synthetic data only): the noise-free feature (dotted black line).

About the CI95:

* It is the 95 % interval of the **mean of the per-beat values**, while the plotted value is the feature of
  the **template**. For non-linear features and fragile fiducials the two can differ.
* It uses 1.96 at every n. For small windows the Student t quantile is larger (4.30 at n = 3, 2.26 at
  n = 10), so the band is too narrow in short windows.
* Successive beats are not independent (respiratory and vasomotor modulation), so the effective n is smaller
  than the number of beats.

Treat the band as a lower bound of the uncertainty.

### Adaptive windows table

One row per adaptive window, all categories together. Columns:

* `category`, `window_id`;
* `t_anchor`: peak time of the anchor beat;
* `t_start`, `t_end`: onset of the first beat and end of the last beat;
* `n_valid`: number of beats averaged;
* `template_sqi`: mean correlation of the beats with the template;
* `anchor_beat`, `T_beats`, `overlap`, `H`;
* each feature and its `_ci95`. The columns of the other categories are empty (NaN).

> **Tip.** Short windows are produced much more often than long ones (one per beat during a transient, one
> every ~25 beats at rest). A plain mean or median of this table is therefore weighted towards transients.
> For summaries over time, use `features_grid.csv` from the **Export** tab, which holds the adaptive values
> on a regular 1 Hz grid.

### Per-beat features table

One row per detected beat (invalid beats included; filter on `valid`): `beat_index`, `t_peak`, `valid`, and
the 33 features measured on that single beat with the same procedure as the templates.

> **Tip.** Template and per-beat values are different estimators. On the sample data, per-beat c/a, d/a, e/a
> and aging index differ systematically from the template values, and per-beat notch-related features of
> BIDMC 01 are unstable ([section 19](#19-reference-values-from-the-sample-data)). Compare template values
> with template values and per-beat values with per-beat values.

---

## 12. Tab: Templates

**Category** selects the category. **Templates shown** (2–30, default 12) selects how many templates are drawn,
evenly spaced in time and coloured from blue (early) to red (late). The legend gives the anchor time and the
number of beats ("t=XXs (n=NN)").

How a template is built:

1. The beats of the window are sampled on a 500 Hz grid by cubic B-spline interpolation of the morphology
   signal.
2. They are aligned on the point of maximum systolic slope (sub-sample precision).
3. The straight foot-to-foot baseline of each beat is removed, and each beat is divided by its amplitude.
4. The beats are averaged, and the mean is multiplied by the mean amplitude.

The template keeps the average shape and the average amplitude, in real time (the beats are not stretched to
a common duration). Its duration is the median time from the foot to the point of maximum slope plus the
median time from that point to the end of the beat. The y axis is the amplitude in a.u. above the
foot-to-foot baseline.

> **Known issue.** The time axis ("Time from foot (s)") is stretched by 500 / processing rate (4 times at
> 125 Hz). The feature values are not affected. Read times from the features, not from this axis
> ([section 20](#20-known-issues)).

What to look for:

* A systolic peak, then a dicrotic notch or a late-systolic or diastolic inflection.
* With age and arterial stiffening, the diastolic wave tends to merge with the systolic part and the notch
  becomes less visible (Millasseau et al. 2002, 2006). In this case the program uses the inflection point.
* Templates that change gradually from blue to red show a slow drift of shape. A sudden jump shows a
  transient or an artefact.

---

## 13. Tab: Benchmark

This tab works only with **Synthetic demo**. With other sources it shows a note that suggests the
`adaptive-ppg benchmark` and `adaptive-ppg sweep` commands for larger studies.

The ground truth is the same scenario without noise, wander, respiration, jitter or artefacts, analysed at the
native rate. Its per-beat features are interpolated linearly on the 1 Hz grid of the analysis. Each
estimator is compared with it after `max(60 s, the longest fixed window)`. Estimators:

* `adaptive`: the adaptive windows;
* `fixed_Ws`: the fixed windows;
* `per_beat`: the last valid beat (zero-order hold);
* `ema_10s`: a causal exponential moving average of the per-beat values with a 10 s time constant.

The features scored are those whose ground truth varies enough in the stress scenario (SD at least 2 % of the
mean absolute value) and which are available at least 80 % of the time. Scales are the ground-truth SDs of
the stress scenario (seed 0).

### Summary table

| Column | Meaning |
| :--- | :--- |
| `estimator` | Estimator name |
| `median_nrmse` | Median over features of RMSE / ground-truth SD (lower is better) |
| `median_stationary_nsd` | Median over features of the SD of the estimate within stationary segments, normalised by the same scale (noise at rest; lower is better) |
| `median_latency_s` | Median over features of the delay (s) of the 50 % crossing at the transitions |
| `mean_coverage` | Mean share of the grid where the estimator has a value |

The table is sorted by `median_nrmse`. The stationary segments of the stress scenario are 60–120 s and
140–240 s. The recovery is not included because it does not settle within 120 s (time constant 25 s). The
latency of a transition is counted only for features that change by at least 0.5 SD at that transition. Its
pre and post levels are the means of the 30 s before the event and the 30 s before the next boundary.

**Metric** (nrmse, stationary_nsd, latency_s, bias) draws one box per estimator over the features, sorted by
the median. Bias is in the feature units, so it cannot be compared across features.

**Per-feature metrics** lists `coverage`, `rmse`, `bias`, `nrmse`, `stationary_sd`, `stationary_nsd`,
`latency_s` and `n_transitions` for every estimator and feature.

Example (stress, seed 1, 125 Hz, defaults):

| Estimator | median_nrmse | median_stationary_nsd | median_latency_s |
| :--- | :---: | :---: | :---: |
| fixed_10s | 0.340 | 0.127 | 8.0 |
| ema_10s | 0.361 | 0.201 | 10.0 |
| adaptive | 0.422 | 0.090 | 8.5 |
| per_beat | 0.494 | 0.474 | −2.0 |
| fixed_30s | 0.619 | 0.342 | 23.0 |

Coverage is 1.0 for every estimator. The adaptive windows have the lowest noise at rest and a latency similar
to 10 s fixed windows. A 10 s window has a lower overall error in this scenario. The per-beat latency of
−2 s is not a real lead: it is produced by noise at the crossing.

> **Caveats.** The synthetic model is simple. The benchmark compares estimators under known conditions; it
> does not validate the features on real physiology. The README reports the results over three seeds and
> for the rest scenario.

---

## 14. Tab: Export

* **Include ensemble templates**: adds the templates of each category to the exports.
* **Download ZIP (CSV + metadata.json)**:

| File | Content |
| :--- | :--- |
| `features_grid.csv` | Adaptive features on a 1 Hz grid, from the first whole second after the first beat to the last beat. Trailing: the value of the latest window (empty before the first one). Centred: linear interpolation between anchors. |
| `windows.csv` | The adaptive windows table ([section 11](#adaptive-windows-table)) |
| `beats.csv` | The beat quality table ([section 9](#beat-quality-table)) |
| `features_per_beat.csv` | The per-beat features |
| `feature_dictionary.csv` | Name, category and unit of each feature |
| `windows_fixed_<W>s.csv` | One table for each fixed window |
| `templates_<category>.csv` | The templates (only with **Include ensemble templates**) |
| `metadata.json` | Run metadata (see below) |

* **Download Excel workbook**: the same tables as sheets (names shortened to 31 characters) plus a
  `metadata` sheet. It needs the `openpyxl` package; otherwise "Excel export needs openpyxl." is shown.
* **Download configuration (JSON)**: `config.json`, with all the parameter values. It can be reused with
  `adaptive-ppg run --config config.json`. It does **not** contain the preset name or the list of fixed
  windows (it holds `fixed_overlap_frac` and a `fixed_window_sec` of 30 that the GUI does not use).
* **Run metadata**: shows `metadata.json`:
  * package, version, git revision (with `-dirty` when tracked files have uncommitted changes), creation
    time (UTC), peak detector;
  * input: channel, format, sampling rate, duration, and the loader details (resolved column, separator,
    decimal mark, time column, estimated rate, maximum time-stamp deviation, source of the rate, rate
    warning);
  * summary: number of beats and of valid beats, processing rate, polarity and skewness, effective low-pass
    cut-off, number of windows per category, IBI regression slope;
  * warnings and the full configuration.

---

## 15. Feature reference

All features are measured on a template (or on a single beat) sampled at 500 Hz. The template is smoothed
with a Savitzky–Golay filter (70 ms, order 3; Savitzky & Golay 1964), which also gives the first (`dy`) and
second (`ddy`) derivatives. Notation: onset `o`, systolic peak `p` (maximum between onset and end, with
parabolic refinement), end `e`, baseline `base` = smoothed value at the onset, `amp` = value at the peak −
`base`. Times use sub-sample positions.

The physiological readings below summarise the cited literature. They describe associations reported in
specific populations and set-ups, not calibrated measurements.

### Macro (pulse amplitude and shape)

| Feature | Definition as implemented | Unit | Physiological reading | Caveats |
| :--- | :--- | :---: | :--- | :--- |
| `peak_amplitude` | `amp` | a.u. | Pulsatile (AC) blood-volume change in the tissue under the sensor. It falls with local vasoconstriction (for example sympathetic activation, cold, pain) and rises with vasodilation (Allen 2007). | Depends on site, wavelength, contact pressure and device gain. Not comparable between subjects or devices. Measured on the filtered signal. |
| `pulse_width_10`, `_25`, `_50`, `_75` | Time between the rising and the falling crossing of 10, 25, 50, 75 % of `amp` | s | The width at half height has been reported to correlate with systemic vascular resistance (Awad et al. 2007). | Shortens with heart rate. The 10 % width is affected by the diastolic wave and the notch. |
| `total_area` | Area of the pulse above `base`, onset to end | a.u.·s | Scales with the pulsatile volume and the beat duration | Same dependence on amplitude as `peak_amplitude` |
| `systolic_area`, `diastolic_area` | Area from onset to systolic peak, and from systolic peak to end | a.u.·s | — | Split at the systolic **peak**, not at the notch |
| `area_ratio` | `diastolic_area / systolic_area` | – | Balance between the late and the early part of the pulse | Depends strongly on the duty cycle, so on heart rate (diastole shortens more than systole as heart rate rises) |
| `inflection_point_area_ratio` (IPA) | Area from the notch to the end divided by the area from the onset to the notch | – | Proposed as an index related to total peripheral resistance (Wang et al. 2009) | Depends on the notch position. Noisy: in the README benchmark its NRMSE exceeds 1 with every estimator. |
| `inflection_amplitude` | Smoothed value at the diastolic point − `base` | a.u. | Height of the diastolic (reflected) wave | The diastolic point is a peak when a notch exists, otherwise an inflection |
| `reflection_index` (RI) | `inflection_amplitude / amp` | – | Mainly reflects the tone of small arteries. It falls with vasodilators such as glyceryl trinitrate and with the β2-agonist salbutamol (Chowienczyk et al. 1999; Millasseau et al. 2006). | Ratio, not percent (multiply by 100 to compare with studies in %). Unstable per beat when there is no clear notch. |
| `max_systolic_slope` | Maximum of `dy` between onset and peak | a.u./s | Steepness of the upstroke | Scales with amplitude |
| `max_decay_slope` | Minimum of `dy` between peak and end (negative) | a.u./s | Steepness of the decay | Scales with amplitude |
| `slope_ratio` | abs(`max_systolic_slope / max_decay_slope`) | – | Asymmetry between rise and fall | — |
| `dc_level` | Mean of the raw signal over the beat (resampled, not filtered, not flipped) | a.u. | With raw light intensity: the non-pulsatile light level (tissue, venous and non-pulsatile arterial blood, and source intensity) | Meaningless when the device removes or rescales the DC |
| `perfusion_index` (PI) | `100 · amp / dc_level`; in windows, template amplitude / mean `dc_level` | % | Ratio of pulsatile to non-pulsatile signal. Low values indicate poor peripheral perfusion (Lima et al. 2002). | Meaningful only for DC-coupled raw intensity from the same channel. Here the AC is measured on the filtered signal, so values are not identical to pulse-oximeter PI. Not physiological for the synthetic data and for BIDMC. |

### Time_volume (timing)

| Feature | Definition as implemented | Unit | Physiological reading | Caveats |
| :--- | :--- | :---: | :--- | :--- |
| `pulse_duration` | Onset to end (next onset) | s | Pulse interval | For a template, the median foot-to-maximum-slope span plus the median maximum-slope-to-end span of its beats, close to the median beat duration |
| `heart_rate` | 60 / `pulse_duration` | bpm | Pulse rate | Pulse rate, not ECG heart rate (Schäfer & Vagedes 2013) |
| `crest_time` | Onset to systolic peak | s | Duration of the upstroke. Used with the stiffness index to classify arterial stiffness (Alty et al. 2007). | Depends on the **Onset method** (shorter with `tangent`) |
| `decay_time` | Systolic peak to end | s | Includes diastole, which shortens most as heart rate rises | — |
| `duty_cycle` | `crest_time / pulse_duration` | – | Fraction of the cycle spent rising | Rises with heart rate |
| `notch_time` | Onset to dicrotic notch; without a notch, onset to the maximum of `ddy` between the peak and the inflection point | s | Timing landmark of late systole at the measurement site | Not a measure of left-ventricular ejection time without validation |
| `delta_t_dvp` (ΔT_DVP) | Systolic peak to diastolic peak (or inflection point) | s | Time between the direct and the reflected wave. Shortens with large-artery stiffness and age (Millasseau et al. 2002). | Needs a resolvable diastolic point |
| `stiffness_index` (SI) | Subject height / `delta_t_dvp` | m/s | Index of large-artery stiffness, correlated with carotid–femoral pulse wave velocity (Millasseau et al. 2002) | Needs **Subject height**. With no diastolic peak the inflection point is used, as in Millasseau et al. (2002). |
| `pat_foot`, `pat_peak` | Preceding ECG R peak to the PPG onset, and to the systolic peak | s | Pulse arrival time = pre-ejection period + pulse transit time (Mukkamala et al. 2015). It shortens with higher blood pressure (shorter transit time) and with higher contractility (shorter pre-ejection period). | Needs an ECG. Measured on each beat at the processing rate; window values are means of the per-beat values. Values outside 0.05–0.6 s (foot) are discarded. Device delays between the ECG and pleth channels add an offset. |

### Derivatives (second derivative, SDPTG)

The SDPTG waves (Takazawa et al. 1998; Elgendi 2012) are found on `ddy` as follows:

* **a**: maximum of `ddy` from 50 ms before the onset to the point of maximum slope;
* **b**: minimum of `ddy` after a, up to the systolic peak;
* **e**: maximum of `ddy` after the peak, within the first 60 % of the beat;
* **c**: the highest local maximum between b and e with a prominence of at least 5 % of a;
* **d**: the lowest local minimum between c and e with the same prominence. If d is not found, c is discarded
  too.

The ratios are computed only when a > 0 and b < 0. When c or d cannot be resolved (common in older subjects),
c/a, d/a and the aging index are empty (NaN).

| Feature | Definition | Physiological reading |
| :--- | :--- | :--- |
| `sdptg_b_a` | b/a | Negative. It increases (towards 0) with age and arterial stiffness (Takazawa et al. 1998). |
| `sdptg_c_a` | c/a | Decreases with age (Takazawa et al. 1998) |
| `sdptg_d_a` | d/a | Decreases with age and changes with vasoactive agents (Takazawa et al. 1998) |
| `sdptg_e_a` | e/a | Decreases with age (Takazawa et al. 1998) |
| `aging_index` | (b − c − d − e)/a | SDPTG aging index; increases with age (Takazawa et al. 1998) |
| `aging_index_be` | (b − e)/a | Alternative index when c and d cannot be resolved (reviewed in Elgendi 2012) |

> **Important.** The SDPTG ratios depend strongly on the derivative smoothing, the low-pass cut-off and the
> sampling rate. On the same synthetic signal the README reports b/a = −0.65 with a 50 ms smoothing window,
> −1.12 with 70 ms and −1.64 with 80 ms. Compare SDPTG values only between analyses with identical settings,
> and compare them with published norms only if the processing is the same. For a review of these indices
> and their limits, see Charlton et al. (2022).

---

## 16. Settings by use case

These are starting points. Check the results with the tools of [section 17](#17-quality-control-checklist).

### Resting recordings and baseline characterisation

* **Preset** `rest`. Type **T_max (s)** 60 / 60 / 90 in section 3 by hand ([section 20](#20-known-issues)).
* **H normalization** fixed with `drift_z`.
* **Window anchor** `centered` for offline analysis.
* **Fixed windows** `30, 60` for comparison.
* Summarise each recording with the median of `features_grid.csv`.
* With the default preset, the resting sample recordings have H > θ in 21–26 % of the beats (`Sample1`,
  `Sample2`). The `rest` preset (θ = 0.5) makes the windows shrink only for larger changes.

### Acute autonomic and stress protocols

Examples: mental arithmetic, cold pressor, Valsalva manoeuvre, handgrip, active standing or tilt.

* **Preset** `acute_stress`. Type **T_max (s)** 20 / 20 / 40 by hand.
* H does not see changes of amplitude or timing alone. To follow pulse amplitude, perfusion index or pulse
  rate within seconds, lower **T_max (s)** of macro and time_volume further (for example 10–15 s).
* Avoid **calibration** normalisation unless the first 60 s are a clean baseline. Prefer **fixed** with
  `drift_z`.
* A smaller **N_past** (for example 6–8) reacts faster but is less sensitive and more affected by breathing.
* For event-locked averages, use `centered` offline, remembering that H still lags about N/2 beats. Use
  `trailing` when causality matters, and remember the lag of about T/2.
* Add an ECG for PAT. Its pre-ejection component shortens with sympathetic activation, so PAT changes
  combine cardiac and vascular effects.
* Add **Fixed windows** `10` for comparison.

### Causal or real-time-like analysis

* **Window anchor** `trailing`. **H normalization** fixed, calibration or rolling (not global).
* The IBI regression is causal. For the first 5 beats only, the quality reference uses the first 30 beats.
* The dashboard analyses whole files. For real-time use, the package provides `StreamingAdaptiveEngine`,
  which takes one beat at a time and returns the adaptive windows as they close (see the README).

### Wrist and low-rate wearables (about 64 Hz)

* **Preset** `wearable_64hz`. Type the derivatives **T_crit** 15 by hand.
* Keep **Processing rate (Hz)** at 125. Upsampling adds no bandwidth: the cut-offs are limited to 45 % of
  64 Hz (28.8 Hz) anyway.
* Below 100 Hz the derivative features are not recommended (README: 0.3 SD difference at 64 Hz). Below 50 Hz
  a warning is shown.
* The pulse shape depends on the measurement site. Notch-based features and reference values from finger
  studies do not transfer directly to the wrist. Motion artefacts are more frequent: check the share of valid
  beats.

### Slow or paced breathing, respiratory sinus arrhythmia, biofeedback

* **N_past** ≥ 2 × heart rate / breathing rate: for 6 breaths/min at 60 bpm, N ≥ 20.
* Consider **Regress out IBI-driven variability** if H follows the breathing.
* Short windows show the respiratory modulation of the features. Long windows average it out.

### Ectopic beats, arrhythmias and atrial fibrillation

* The quality checks (`ibi_ok`, `corr_ok`) usually reject ectopic beats and the beats after them. Do not
  lower **Min. correlation with the reference beat** to keep them.
* In atrial fibrillation the intervals are irregular and the pulse amplitude varies with the preceding
  interval. Many beats are rejected and ensemble averaging is not meaningful. The method is not designed for
  this rhythm: analyse sinus-rhythm segments only.

### Clinical monitor exports (MIMIC, BIDMC)

* The pleth is already filtered, scaled and shown volume up, usually at 125 Hz. The defaults work.
* `perfusion_index` and `dc_level` are meaningless (BIDMC 01 gives a PI of about 37 %).
* PAT may include processing delays between the ECG and the pleth channels of the monitor. The BIDMC 01
  value (`pat_foot` ≈ 0.53 s) is longer than the values usually reported for the finger, roughly 0.2–0.4 s.
  Interpret changes within a recording rather than absolute values.

### Bradycardia and tachycardia

* Below about 45 bpm, lower **Morphology high-pass (Hz)** to 0.3 Hz. Beats below 30 bpm are always rejected.
* The Elgendi detector cannot find peaks closer than 0.3 s, so it is limited to 200 bpm. Beats above 220 bpm
  are always rejected. The defaults have not been validated for neonates or small children.
* **T_max (s)** is in seconds: at high heart rates the same window holds more beats.

### Comparing subjects, groups or sessions

* Use the same configuration (save the JSON) and the same fixed windows.
* Keep the processing rate, the low-pass, the onset method and the derivative smoothing the same.
* Compare template values with template values, and per-beat values with per-beat values.
* Summarise with time-weighted values (`features_grid.csv`), not with the windows table.
* Amplitude, areas, slopes, DC and PI are in arbitrary units and depend on the device and on contact. Prefer
  within-subject changes, or dimensionless and timing features.
* Timing features and the area ratio depend on heart rate. Consider heart rate as a covariate.
* Report the share of valid beats, and set an exclusion rule in advance.

---

## 17. Quality-control checklist

1. The info line shows the right channel, sampling rate and duration. No rate-mismatch warning.
2. The filtered pulse points upwards. If not, set **Signal polarity** by hand.
3. Peaks and onsets are on the right points when you zoom on some beats.
4. **Valid beats** is above about 90 %. The rejected spans match visible artefacts.
5. **Median template SQI** is close to 1 (the sample data give 0.985–0.993 for macro).
6. The warm-up is short compared with the recording (19 valid beats with the defaults).
7. **Beats with H > θ** is plausible: low at rest, higher during tasks.
8. The templates have a plausible shape. Check whether a notch is visible.
9. For SDPTG features: check how many beats have empty c/a and d/a in **Per-beat features**.
10. For PAT: the values are not piled up at the limits (0.05 or 0.6 s), and the R-peak detection is
    plausible.

---

## 18. What to report

* Package version and git revision (`metadata.json`).
* The preset and every changed parameter (`config.json`), plus the fixed windows, which are not in the JSON.
* Processing rate, filter bands (and the effective low-pass cut-off), onset method, peak detector.
* Derivative smoothing (Savitzky–Golay 70 ms, order 3, unless changed) and the 500 Hz template grid.
* Quality thresholds and the share of valid beats.
* Variability metric, H normalisation, θ, slope, N_past, anchor and window limits per category.
* Whether the reported values are template or per-beat values, and how they were summarised over time.

---

## 19. Reference values from the sample data

Default settings, processing at 125 Hz, fixed windows 10 and 30 s. "Template" is the median over the adaptive
windows of the category. "Per beat" is the median over the valid beats. These values describe the three
sample recordings. They are not normative values.

| | `Sample1.CSV` | `Sample2.CSV` | `bidmc_01_Signals.csv` |
| :--- | :---: | :---: | :---: |
| Valid beats | 335 / 335 | 276 / 277 | 705 / 718 |
| Pulse rate (bpm) | 71.1 | 73.7 | 90.9 |
| Median `drift_z` | 1.18 | 1.06 | 0.83 |
| Beats with H > θ | 25.9 % | 21.0 % | 7.2 % |
| Macro window, median over all beats (p10–p90 after warm-up) | 18 (4–33) | 17 (5–34) | 36 (13–44) |
| Macro window, median over windows | 8 | 10.5 | 29 |
| Number of macro windows | 47 | 32 | 40 |
| Derivatives window, median over all beats / over windows | 25 / 16.5 | 25 / 24 | 60 / 46 |
| Template SQI, macro / derivatives | 0.985 / 0.958 | 0.989 / 0.982 | 0.993 / 0.993 |
| Warm-up | 19 beats (15.8 s) | 20 beats (15.0 s) | 19 beats (11.9 s) |

The medians over all beats include the warm-up, as the **Median macro window** KPI does. The README gives the
medians after the warm-up (macro: 17, 18 and 35 beats; derivatives: 24, 26 and 59 beats).

Features (template / per beat):

| Feature | `Sample1.CSV` | `Sample2.CSV` | `bidmc_01_Signals.csv` |
| :--- | :---: | :---: | :---: |
| `crest_time` (s) | 0.095 / 0.096 | 0.093 / 0.095 | 0.111 / 0.111 |
| `duty_cycle` | 0.112 / 0.115 | 0.116 / 0.116 | 0.167 / 0.167 |
| `notch_time` (s) | 0.285 / 0.302 | 0.261 / 0.294 | 0.224 / 0.224 |
| `delta_t_dvp` (s) | 0.242 / 0.247 | 0.249 / 0.247 | 0.155 / 0.339 |
| `reflection_index` | 0.657 / 0.665 | 0.628 / 0.665 | 0.289 / −0.033 |
| `pulse_width_50` (s) | 0.429 / 0.424 | 0.404 / 0.416 | 0.166 / 0.166 |
| `area_ratio` | 5.63 / 5.29 | 5.36 / 5.31 | 1.60 / 1.64 |
| `inflection_point_area_ratio` | 0.81 / 0.70 | 0.88 / 0.71 | 0.147 / 0.161 |
| `slope_ratio` | 6.65 / 4.99 | 6.38 / 4.78 | 1.82 / 1.68 |
| `sdptg_b_a` | −0.724 / −0.732 | −0.702 / −0.709 | −0.874 / −0.898 |
| `sdptg_c_a` | 0.058 / 0.116 | 0.080 / 0.119 | −0.413 / −0.307 |
| `sdptg_d_a` | −0.099 / −0.156 | −0.105 / −0.165 | −0.484 / −0.527 |
| `sdptg_e_a` | 0.100 / 0.190 | 0.093 / 0.187 | 0.451 / 0.499 |
| `aging_index` | −0.787 / −0.876 | −0.777 / −0.860 | −0.443 / −0.580 |
| `aging_index_be` | −0.831 / −0.922 | −0.789 / −0.904 | −1.326 / −1.397 |
| `perfusion_index` (%) | 0.337 / 0.381 | 0.353 / 0.357 | 36.7 / 37.6 (not meaningful) |

BIDMC 01 with ECG lead II: `pat_foot` ≈ 0.525 s, `pat_peak` ≈ 0.637 s (medians).

Observations:

* Per-beat c/a, d/a and the aging index are empty in 15 % (`Sample1`), 16 % (`Sample2`) and 22 % (BIDMC 01)
  of the beats.
* Per-beat c/a, d/a, e/a and aging index differ systematically from the template values: averaging smooths
  the small c and d waves.
* In BIDMC 01 the per-beat ΔT_DVP and RI are unstable (the per-beat median RI is negative). Use the template
  values.
* Per-check pass rates in BIDMC 01: `rate_ok` 99.6 %, `ibi_ok` 98.5 %, `corr_ok` 99.3 %, `amp_ok` 99.0 %.
* The very different `area_ratio` and `pulse_width_50` of BIDMC 01 reflect a different pulse contour (higher
  heart rate, critically ill patient, monitor filtering). They should not be read as a difference between
  sensors.

---

## 20. Known issues

These issues were found while checking this guide against the code of version 0.2.0. They are not fixed in
that version.

1. **Presets do not update section 3.** The **T_min (beats)**, **T_max (s)** and **T_crit** controls keep
   their earlier values when you change **Preset** after the first run. Type the `rest`, `acute_stress` and
   `wearable_64hz` values of section 3 by hand ([section 4](#preset)). The other controls follow the preset,
   and the values without a control (growth limit, derivative smoothing) always apply.
2. **Template time axis.** In the **Templates** tab, and in the `t_from_foot_s` column of the exported
   templates, the time axis is divided by the processing rate instead of the 500 Hz template rate. It is
   stretched by 500 / processing rate (4 times at 125 Hz; correct when processing at 500 Hz). Feature values
   are not affected.
3. **README: beat normalisation.** The README says that each beat is scaled to [0, 1]. The code subtracts the
   straight line from the onset to the end and scales the peak to 1, so values below 0 are possible.
4. **README: centred windows.** The README says that `centered` gives zero-lag windows. The feature windows
   are symmetric, but H still reacts about N_past/2 beats after a change.
5. **Settings not saved.** `config.json` does not contain the preset name or the list of fixed windows. The
   **Benchmark** tab always uses the Elgendi detector, whatever **Peak detector** says.

---

## 21. References

* Allen J. Photoplethysmography and its application in clinical physiological measurement. *Physiol Meas*
  2007;28:R1–R39.
* Alty SR, et al. Predicting arterial stiffness from the digital volume pulse waveform. *IEEE Trans Biomed
  Eng* 2007;54:2268–2275.
* Awad AA, et al. The relationship between the photoplethysmographic waveform and systemic vascular
  resistance. *J Clin Monit Comput* 2007;21:365–372.
* Charlton PH, et al. Assessing hemodynamics from the photoplethysmogram to gain insights into vascular age:
  a review from VascAgeNet. *Am J Physiol Heart Circ Physiol* 2022;322:H493–H522.
* Chowienczyk PJ, et al. Photoplethysmographic assessment of pulse wave reflection: blunted response to
  endothelium-dependent beta2-adrenergic vasodilation in type II diabetes mellitus. *J Am Coll Cardiol*
  1999;34:2007–2014.
* Elgendi M. On the analysis of fingertip photoplethysmogram signals. *Curr Cardiol Rev* 2012;8:14–25.
* Elgendi M, et al. Systolic peak detection in acceleration photoplethysmograms measured from emergency
  responders in tropical conditions. *PLoS ONE* 2013;8:e76585.
* Goldberger AL, et al. PhysioBank, PhysioToolkit, and PhysioNet: components of a new research resource for
  complex physiologic signals. *Circulation* 2000;101:e215–e220.
* Julien C. The enigma of Mayer waves: facts and models. *Cardiovasc Res* 2006;70:12–21.
* Li Q, Clifford GD. Dynamic time warping and machine learning for signal quality assessment of pulsatile
  signals. *Physiol Meas* 2012;33:1491–1501.
* Lima AP, Beelen P, Bakker J. Use of a peripheral perfusion index derived from the pulse oximetry signal as a
  noninvasive indicator of perfusion. *Crit Care Med* 2002;30:1210–1213.
* Mehrgardt P, Khushi M, Poon S, Withana A. Pulse Transit Time PPG Dataset (version 1.1.0). PhysioNet 2022.
  doi:10.13026/jpan-6n92.
* Millasseau SC, et al. Determination of age-related increases in large artery stiffness by digital pulse
  contour analysis. *Clin Sci* 2002;103:371–377.
* Millasseau SC, et al. Contour analysis of the photoplethysmographic pulse measured at the finger.
  *J Hypertens* 2006;24:1449–1456.
* Mukkamala R, et al. Toward ubiquitous blood pressure monitoring via pulse transit time: theory and
  practice. *IEEE Trans Biomed Eng* 2015;62:1879–1901.
* Orphanidou C, et al. Signal-quality indices for the electrocardiogram and photoplethysmogram: derivation and
  applications to wireless monitoring. *IEEE J Biomed Health Inform* 2015;19:832–838.
* Pimentel MAF, et al. Toward a robust estimation of respiratory rate from pulse oximeters. *IEEE Trans Biomed
  Eng* 2017;64:1914–1923.
* Savitzky A, Golay MJE. Smoothing and differentiation of data by simplified least squares procedures.
  *Anal Chem* 1964;36:1627–1639.
* Schäfer A, Vagedes J. How accurate is pulse rate variability as an estimate of heart rate variability? A
  review on studies comparing photoplethysmographic technology with an electrocardiogram. *Int J Cardiol*
  2013;166:15–29.
* Takazawa K, et al. Assessment of vasoactive agents and vascular aging by the second derivative of
  photoplethysmogram waveform. *Hypertension* 1998;32:365–370.
* Wang L, et al. Noninvasive cardiac output estimation using a novel photoplethysmogram index. *Proc IEEE Eng
  Med Biol Soc (EMBC)* 2009:1746–1749.
