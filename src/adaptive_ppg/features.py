"""Morphological features from ensemble-averaged PPG templates.

For every window produced by the adaptive engine, the valid beats are
**aligned on their maximum systolic slope**, baseline-corrected (linear
foot-to-foot baseline), scaled to unit amplitude and **averaged into a
template**; the features are then extracted *from the template*. Averaging
waveforms (ensemble averaging) is not the same as averaging per-beat
features, f(mean x) != mean f(x): the template has a higher SNR, which is what
makes second-derivative (SDPTG) fiducials usable.

Beats are aligned with sub-sample precision (fractional fiducials and cubic
B-spline interpolation of the signal) and sampled on a fixed template grid
(``FeatureConfig.template_fs``, 500 Hz by default) whatever the processing
rate; timing fiducials on the template are refined with parabolic
interpolation. The Savitzky-Golay derivatives therefore use the same filter
at every input rate, and timing features are not quantized to ``1 / fs``.

Each feature belongs to one category (:data:`FEATURE_CATEGORIES`) and is
computed with that category's window schedule. Per-beat features (single-beat
"templates") are kept as a secondary output and give the per-window 95 %
confidence interval ``1.96 * SD / sqrt(n)``.
"""

from __future__ import annotations

import functools
import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.integrate import trapezoid
from scipy.ndimage import map_coordinates, spline_filter1d
from scipy.signal import butter, find_peaks, savgol_coeffs, savgol_filter, sosfiltfilt

from ._typing import FloatArray, IntArray
from .config import FeatureConfig
from .engine import AdaptiveResult, Window
from .preprocessing import ProcessedPPG, parabolic_offset

logger = logging.getLogger(__name__)

FEATURE_CATEGORIES: Dict[str, Tuple[str, ...]] = {
    "macro": (
        "peak_amplitude",
        "pulse_width_10",
        "pulse_width_25",
        "pulse_width_50",
        "pulse_width_75",
        "total_area",
        "systolic_area",
        "diastolic_area",
        "area_ratio",
        "inflection_point_area_ratio",
        "inflection_amplitude",
        "reflection_index",
        "max_systolic_slope",
        "max_decay_slope",
        "slope_ratio",
        "dc_level",
        "perfusion_index",
    ),
    "time_volume": (
        "pulse_duration",
        "heart_rate",
        "crest_time",
        "decay_time",
        "duty_cycle",
        "notch_time",
        "delta_t_dvp",
        "stiffness_index",
        "pat_foot",
        "pat_peak",
    ),
    "derivatives": (
        "sdptg_b_a",
        "sdptg_c_a",
        "sdptg_d_a",
        "sdptg_e_a",
        "aging_index",
        "aging_index_be",
    ),
}

FEATURE_TO_CATEGORY: Dict[str, str] = {f: c for c, fs in FEATURE_CATEGORIES.items() for f in fs}
ALL_FEATURES: Tuple[str, ...] = tuple(f for fs in FEATURE_CATEGORIES.values() for f in fs)

FEATURE_UNITS: Dict[str, str] = {
    "peak_amplitude": "a.u.", "pulse_width_10": "s", "pulse_width_25": "s", "pulse_width_50": "s",
    "pulse_width_75": "s", "total_area": "a.u.·s", "systolic_area": "a.u.·s", "diastolic_area": "a.u.·s",
    "area_ratio": "-", "inflection_point_area_ratio": "-", "inflection_amplitude": "a.u.",
    "reflection_index": "-", "max_systolic_slope": "a.u./s", "max_decay_slope": "a.u./s", "slope_ratio": "-",
    "dc_level": "a.u.", "perfusion_index": "%", "pulse_duration": "s", "heart_rate": "bpm", "crest_time": "s",
    "decay_time": "s", "duty_cycle": "-", "notch_time": "s", "delta_t_dvp": "s", "stiffness_index": "m/s",
    "pat_foot": "s", "pat_peak": "s", "sdptg_b_a": "-", "sdptg_c_a": "-", "sdptg_d_a": "-", "sdptg_e_a": "-",
    "aging_index": "-", "aging_index_be": "-",
}

# Features that are properties of a beat sequence rather than of the waveform shape:
# they are averaged over the window's beats instead of being read from the template.
_SEQUENCE_FEATURES = ("dc_level", "perfusion_index", "pat_foot", "pat_peak")


# ---------------------------------------------------------------------- #
# Templates
# ---------------------------------------------------------------------- #
@dataclass
class Template:
    """Ensemble-averaged beat.

    ``y`` is in signal units with the foot-to-foot baseline removed; it spans
    ``margin`` extra samples on both sides of ``[onset, end]`` so that the
    Savitzky-Golay derivatives have no edge effects inside the beat.
    ``onset`` / ``end`` are sample indices; ``onset_pos`` / ``end_pos`` are the
    fractional positions used for the timing features.
    """

    y: FloatArray
    fs: float
    onset: int
    end: int
    n_beats: int
    sqi: float
    amplitude: float
    beat_indices: IntArray = field(default_factory=lambda: np.array([], dtype=int))
    onset_pos: float = float("nan")
    end_pos: float = float("nan")

    def __post_init__(self) -> None:
        if not np.isfinite(self.onset_pos):
            self.onset_pos = float(self.onset)
        if not np.isfinite(self.end_pos):
            self.end_pos = float(self.end)

    @property
    def time(self) -> FloatArray:
        return (np.arange(len(self.y)) - self.onset_pos) / self.fs

    @property
    def core(self) -> FloatArray:
        return self.y[self.onset : self.end + 1]


def _sg_window(fs: float, cfg: FeatureConfig) -> int:
    w = int(round(cfg.sg_window_sec * fs))
    w = max(w, cfg.sg_polyorder + 2)
    return w if w % 2 == 1 else w + 1


@functools.lru_cache(maxsize=16)
def _sg_coeffs(window: int, polyorder: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Savitzky-Golay convolution kernels for the signal and its first two derivatives.

    Templates carry a margin of one filter window on each side, so plain
    convolution (no edge fitting) is exact inside the beat.
    """
    return tuple(savgol_coeffs(window, polyorder, deriv=d, use="conv") for d in range(3))


class TemplateBuilder:
    """Builds aligned ensemble templates from a :class:`ProcessedPPG`.

    Templates are sampled at ``cfg.template_fs`` (the processing rate when
    ``None``); ``fs``, ``margin`` and every template index refer to that grid.
    """

    def __init__(self, ppg: ProcessedPPG, cfg: FeatureConfig):
        self.ppg = ppg
        self.cfg = cfg
        fs = ppg.sampling_rate
        self.fs = float(cfg.template_fs or fs)
        self._ratio = self.fs / fs  # template samples per signal sample
        self.win = _sg_window(self.fs, cfg)
        self.margin = self.win
        y = ppg.filtered_signal
        self._pad = int(np.ceil(self.margin / self._ratio)) + int(fs)
        # Cubic B-spline coefficients: beats are sampled at fractional positions.
        self._coef = spline_filter1d(np.pad(y, self._pad, mode="edge"), order=3, mode="mirror")
        dy = savgol_filter(y, _sg_window(fs, cfg), cfg.sg_polyorder, deriv=1)
        # Alignment point: maximum systolic slope between foot and peak, refined to sub-sample precision.
        ms = [o + int(np.argmax(dy[o : p + 1])) for o, p in zip(ppg.onsets, ppg.peaks)]
        self.max_slope = np.array([k + parabolic_offset(dy, k) for k in ms], dtype=float)
        self._cache: Dict[Tuple[int, ...], Template] = {}

    def _sample(self, positions: np.ndarray) -> np.ndarray:
        """Filtered signal at fractional sample ``positions`` (cubic B-spline)."""
        flat = np.asarray(positions, dtype=float).ravel() + self._pad
        return map_coordinates(self._coef, flat[None], order=3, prefilter=False, mode="mirror").reshape(
            np.shape(positions))

    def build(self, beat_indices: Sequence[int]) -> Template:
        key = tuple(int(i) for i in beat_indices)
        if key in self._cache:
            return self._cache[key]
        idx = np.asarray(key, dtype=int)
        ppg, m = self.ppg, self.margin
        on, pk, en = ppg.onset_positions[idx], ppg.peak_positions[idx], ppg.end_positions[idx]
        ms = self.max_slope[idx]
        # Median foot-to-max-slope and max-slope-to-end spans, in template samples.
        pre_f = float(np.median(ms - on)) * self._ratio
        post_f = float(np.median(en - ms)) * self._ratio
        pre, post = int(np.ceil(pre_f)), int(np.ceil(post_f))
        length = pre + post + 2 * m + 1
        # Row r samples beat r (in signal-sample units) on a grid that puts its max-slope point at index m + pre.
        absolute = ms[:, None] + ((np.arange(length) - pre - m) / self._ratio)[None, :]
        segs = self._sample(absolute)
        y_o, y_e, y_p = self._sample(on), self._sample(en), self._sample(pk)
        span = np.maximum(en - on, 1.0)
        segs -= y_o[:, None] + (y_e - y_o)[:, None] * (absolute - on[:, None]) / span[:, None]
        amps = y_p - (y_o + (y_e - y_o) * (pk - on) / span)
        # The peak of an atypical beat may not rise above its baseline.
        bad = ~np.isfinite(amps) | (amps <= 0)
        amps[bad] = np.maximum(np.ptp(segs[bad], axis=1), 1e-12)
        segs /= amps[:, None]
        mean_norm = segs.mean(axis=0)
        core = slice(m, m + pre + post + 1)
        if len(idx) > 1:
            ref = mean_norm[core] - mean_norm[core].mean()
            c = segs[:, core] - segs[:, core].mean(axis=1, keepdims=True)
            denom = np.linalg.norm(c, axis=1) * np.linalg.norm(ref)
            sqi = float(np.mean(np.where(denom > 0, c @ ref / np.where(denom > 0, denom, 1), 0.0)))
        else:
            sqi = 1.0
        amplitude = float(np.mean(amps))
        onset_pos, end_pos = m + pre - pre_f, m + pre + post_f
        tpl = Template(mean_norm * amplitude, self.fs, int(round(onset_pos)), int(round(end_pos)), len(idx), sqi,
                       amplitude, idx, onset_pos, end_pos)
        self._cache[key] = tpl
        return tpl


# ---------------------------------------------------------------------- #
# Feature extraction from one template
# ---------------------------------------------------------------------- #
def _crossings(y: np.ndarray, o: int, p: int, e: int, level: float) -> Tuple[float, float]:
    """Fractional sample positions where the rising / falling edges cross ``level``."""
    rise = np.nan
    for i in range(p, o, -1):
        if y[i - 1] < level <= y[i]:
            rise = i - 1 + (level - y[i - 1]) / (y[i] - y[i - 1])
            break
    fall = np.nan
    for i in range(p, e):
        if y[i] >= level > y[i + 1]:
            fall = i + (y[i] - level) / (y[i] - y[i + 1])
            break
    return rise, fall


def find_dicrotic(ys: FloatArray, dy: FloatArray, ddy: FloatArray, o: int, p: int, e: int, fs: float,
                  amplitude: float) -> Tuple[Optional[int], Optional[int], str]:
    """Locates the dicrotic notch and the diastolic point.

    Search window: from 50 ms after the systolic peak to 75 % of the beat.
    If the smoothed template has a notch (local minimum followed by a local
    maximum) they are returned directly (``"notch"``). Otherwise the diastolic
    point is the inflection point where the decay slows most (local maximum of
    the first derivative) and the notch is the SDPTG e-wave before it
    (``"inflection"``).
    """
    lo = p + max(int(round(0.05 * fs)), 1)
    hi = min(o + int(round(0.75 * (e - o))), e - 1)
    if hi - lo < 3:
        return None, None, "none"
    seg = ys[lo : hi + 1]
    prom = 0.005 * amplitude
    minima, _ = find_peaks(-seg, prominence=prom)
    for mi in minima:
        maxima, _ = find_peaks(seg[mi:], prominence=prom)
        if len(maxima):
            return lo + int(mi), lo + int(mi + maxima[0]), "notch"
    dseg = dy[lo : hi + 1]
    peaks, _ = find_peaks(dseg)
    if len(peaks):
        infl = lo + int(peaks[np.argmax(dseg[peaks])])
        notch = p + int(np.argmax(ddy[p : infl + 1]))
        return notch, infl, "inflection"
    return None, None, "none"


def sdptg_points(ddy: FloatArray, dy: FloatArray, o: int, p: int, e: int, fs: float) -> Dict[str, Optional[int]]:
    """Constrained a-e wave detection on the second derivative.

    * a: maximum of ddy in the upstroke, from 50 ms before the foot to the max-slope point;
    * b: minimum of ddy after a, up to the systolic peak;
    * e: maximum of ddy after the peak, within 60 % of the beat;
    * c, d: the highest local maximum between b and e and the lowest minimum
      between c and e, each with a prominence of at least 5 % of a (``None``
      when they are not resolvable, as is common in older subjects).
    """
    ms = o + int(np.argmax(dy[o : p + 1]))
    w0 = max(o - int(round(0.05 * fs)), 0)
    a = w0 + int(np.argmax(ddy[w0 : ms + 1]))
    if p <= a:
        return {"a": a, "b": None, "c": None, "d": None, "e": None}
    b = a + 1 + int(np.argmin(ddy[a + 1 : p + 1]))
    e_hi = min(o + int(round(0.6 * (e - o))), len(ddy) - 1)
    if e_hi <= p + 1:
        return {"a": a, "b": b, "c": None, "d": None, "e": None}
    ew = p + 1 + int(np.argmax(ddy[p + 1 : e_hi + 1]))
    c = d = None
    prom = 0.05 * max(float(ddy[a]), 0.0)
    between = ddy[b : ew + 1]
    maxima, _ = find_peaks(between, prominence=prom)
    if len(maxima):
        c = b + int(maxima[np.argmax(between[maxima])])
        after_c = ddy[c : ew + 1]
        minima, _ = find_peaks(-after_c, prominence=prom)
        if len(minima):
            d = c + int(minima[np.argmin(after_c[minima])])
        else:
            c = None
    return {"a": a, "b": b, "c": c, "d": d, "e": ew}


def template_features(tpl: Template, cfg: FeatureConfig) -> Dict[str, float]:
    """Extracts the shape features of one template (single beat or ensemble)."""
    fs = tpl.fs
    win = _sg_window(fs, cfg)
    y = tpl.y
    nan = float("nan")
    out = {f: nan for f in ALL_FEATURES if f not in _SEQUENCE_FEATURES}
    if len(y) <= win:
        return out
    ys, dy, ddy = (np.convolve(y, c, mode="same") * fs ** d for d, c in enumerate(_sg_coeffs(win, cfg.sg_polyorder)))
    o, e = tpl.onset, tpl.end
    o_f, e_f = tpl.onset_pos, tpl.end_pos
    p = o + int(np.argmax(ys[o : e + 1]))
    p_f = p + parabolic_offset(ys, p)
    base = float(np.interp(o_f, np.arange(len(ys)), ys))
    amp = ys[p] - base
    T = (e_f - o_f) / fs
    if amp <= 0 or p <= o or p >= e:
        return out

    out["pulse_duration"] = T
    out["heart_rate"] = 60.0 / T
    out["crest_time"] = (p_f - o_f) / fs
    out["decay_time"] = (e_f - p_f) / fs
    out["duty_cycle"] = (p_f - o_f) / (e_f - o_f)
    out["peak_amplitude"] = amp
    for pct in (10, 25, 50, 75):
        rise, fall = _crossings(ys - base, o, p, e, amp * pct / 100.0)
        out[f"pulse_width_{pct}"] = (fall - rise) / fs if np.isfinite(rise) and np.isfinite(fall) else nan

    yc = ys[o : e + 1] - base
    out["total_area"] = float(trapezoid(yc, dx=1.0 / fs))
    sys_area = float(trapezoid(yc[: p - o + 1], dx=1.0 / fs))
    out["systolic_area"] = sys_area
    out["diastolic_area"] = out["total_area"] - sys_area
    out["area_ratio"] = out["diastolic_area"] / sys_area if sys_area > 0 else nan

    out["max_systolic_slope"] = float(np.max(dy[o : p + 1]))
    out["max_decay_slope"] = float(np.min(dy[p : e + 1]))
    out["slope_ratio"] = abs(out["max_systolic_slope"] / out["max_decay_slope"]) if out["max_decay_slope"] < 0 else nan

    notch, dia, kind = find_dicrotic(ys, dy, ddy, o, p, e, fs, amp)
    if notch is not None and dia is not None:
        # Refine on the curve whose extremum defines each point.
        notch_f = notch + parabolic_offset(ys if kind == "notch" else ddy, notch)
        dia_f = dia + parabolic_offset(ys if kind == "notch" else dy, dia)
        out["notch_time"] = (notch_f - o_f) / fs
        out["inflection_amplitude"] = ys[dia] - base
        out["reflection_index"] = out["inflection_amplitude"] / amp
        dt_dvp = (dia_f - p_f) / fs
        out["delta_t_dvp"] = dt_dvp
        if cfg.subject_height_m and dt_dvp > 0:
            out["stiffness_index"] = cfg.subject_height_m / dt_dvp
        a1 = float(trapezoid(yc[: notch - o + 1], dx=1.0 / fs))
        a2 = out["total_area"] - a1
        out["inflection_point_area_ratio"] = a2 / a1 if a1 > 0 else nan

    pts = sdptg_points(ddy, dy, o, p, e, fs)
    a_val = ddy[pts["a"]]
    b_val = ddy[pts["b"]] if pts["b"] is not None else nan
    if a_val > 0 and np.isfinite(b_val) and b_val < 0:
        out["sdptg_b_a"] = b_val / a_val
        c_val = ddy[pts["c"]] if pts["c"] is not None else nan
        d_val = ddy[pts["d"]] if pts["d"] is not None else nan
        e_val = ddy[pts["e"]] if pts["e"] is not None else nan
        out["sdptg_c_a"] = c_val / a_val
        out["sdptg_d_a"] = d_val / a_val
        out["sdptg_e_a"] = e_val / a_val
        out["aging_index"] = (b_val - c_val - d_val - e_val) / a_val
        out["aging_index_be"] = (b_val - e_val) / a_val
    return out


# ---------------------------------------------------------------------- #
# ECG R peaks and pulse arrival time
# ---------------------------------------------------------------------- #
def detect_r_peaks(ecg: FloatArray, fs: float) -> IntArray:
    """Simple R-peak detector (band-pass 5-20 Hz, squared slope, adaptive height).

    Adequate for clean monitoring leads such as BIDMC lead II; use a dedicated
    detector for ambulatory ECG.
    """
    high = min(20.0, 0.45 * fs)
    sos = butter(3, [5.0 / (fs / 2), high / (fs / 2)], btype="band", output="sos")
    f = sosfiltfilt(sos, ecg - np.median(ecg))
    energy = np.gradient(f) ** 2
    energy = np.convolve(energy, np.ones(max(int(0.08 * fs), 1)) / max(int(0.08 * fs), 1), mode="same")
    peaks, _ = find_peaks(energy, height=0.3 * np.percentile(energy, 99), distance=int(0.3 * fs))
    polarity = 1.0 if np.median(f[peaks]) >= 0 else -1.0 if len(peaks) else 1.0
    half = max(int(0.05 * fs), 1)
    refined = [max(p - half, 0) + int(np.argmax(polarity * f[max(p - half, 0) : p + half + 1])) for p in peaks]
    return np.unique(np.asarray(refined, dtype=int))


def pulse_arrival_times(ppg: ProcessedPPG, r_peaks: IntArray) -> Tuple[FloatArray, FloatArray]:
    """PAT from the preceding R peak to the PPG foot and to the systolic peak (s).

    Values outside 50-600 ms (foot) are set to NaN (missed or spurious R peaks).
    """
    fs = ppg.sampling_rate
    pat_foot = np.full(ppg.n_beats, np.nan)
    pat_peak = np.full(ppg.n_beats, np.nan)
    if len(r_peaks) == 0:
        return pat_foot, pat_peak
    pos = np.searchsorted(r_peaks, ppg.onset_positions, side="right") - 1
    ok = pos >= 0
    r = np.where(ok, r_peaks[np.maximum(pos, 0)], 0)
    pf = (ppg.onset_positions - r) / fs
    pp = (ppg.peak_positions - r) / fs
    good = ok & (pf >= 0.05) & (pf <= 0.6)
    pat_foot[good] = pf[good]
    pat_peak[good] = pp[good]
    return pat_foot, pat_peak


# ---------------------------------------------------------------------- #
# Results
# ---------------------------------------------------------------------- #
@dataclass
class WindowSet:
    """Features of one window schedule (one category, or a fixed-window baseline)."""

    name: str
    features: Tuple[str, ...]
    table: pd.DataFrame
    templates: List[Template]
    anchor: str


@dataclass
class FeatureResult:
    per_beat: pd.DataFrame
    adaptive: Dict[str, WindowSet]
    grid: pd.DataFrame
    fixed: Dict[str, WindowSet]
    fixed_grid: pd.DataFrame
    r_peaks: Optional[IntArray] = None


def series_on_grid(times: FloatArray, values: FloatArray, grid: FloatArray, anchor: str) -> FloatArray:
    """Resamples a window series on ``grid``.

    Trailing windows use a zero-order hold (the latest window whose anchor is
    ``<= t``; NaN before the first) so that the grid stays causal. Centered
    windows are linearly interpolated between anchors.
    """
    out = np.full(len(grid), np.nan)
    ok = np.isfinite(values)
    if ok.sum() == 0:
        return out
    t, v = times[ok], values[ok]
    if anchor == "centered":
        inside = (grid >= t[0]) & (grid <= t[-1])
        out[inside] = np.interp(grid[inside], t, v)
    else:
        pos = np.searchsorted(t, grid, side="right") - 1
        has = pos >= 0
        out[has] = v[pos[has]]
    return out


class PPGFeatureExtractor:
    """Computes template features on adaptive and fixed window schedules."""

    def __init__(self, config: Optional[FeatureConfig] = None, t_crit_fixed: int = 3):
        self.config = config or FeatureConfig()
        self.t_crit_fixed = t_crit_fixed

    # -- per beat ------------------------------------------------------ #
    def per_beat(self, ppg: ProcessedPPG, builder: TemplateBuilder, pat: Tuple[np.ndarray, np.ndarray],
                 dc: np.ndarray) -> pd.DataFrame:
        rows = []
        for k in range(ppg.n_beats):
            feats = template_features(builder.build([k]), self.config)
            feats["dc_level"] = dc[k]
            feats["perfusion_index"] = 100.0 * feats["peak_amplitude"] / dc[k] if dc[k] > 0 else np.nan
            feats["pat_foot"], feats["pat_peak"] = pat[0][k], pat[1][k]
            rows.append(feats)
        df = pd.DataFrame(rows, columns=list(ALL_FEATURES))
        df.insert(0, "valid", ppg.valid_beats_mask)
        df.insert(0, "t_peak", ppg.beat_times)
        df.insert(0, "beat_index", np.arange(ppg.n_beats))
        return df

    # -- windows ------------------------------------------------------- #
    def window_table(self, name: str, windows: Sequence[Tuple[float, float, float, np.ndarray, Dict]],
                     features: Tuple[str, ...], builder: TemplateBuilder, per_beat: pd.DataFrame,
                     anchor: str) -> WindowSet:
        """``windows`` items: (t_anchor, t_start, t_end, beat_indices, extra columns)."""
        rows, templates = [], []
        pb = per_beat[list(features)].to_numpy()
        seq_cols = [i for i, f in enumerate(features) if f in _SEQUENCE_FEATURES]
        for wid, (t_anchor, t_start, t_end, idx, extra) in enumerate(windows):
            tpl = builder.build(idx)
            feats = template_features(tpl, self.config)
            beat_vals = pb[idx]
            row = {"window_id": wid, "t_anchor": t_anchor, "t_start": t_start, "t_end": t_end,
                   "n_valid": len(idx), "template_sqi": tpl.sqi, **extra}
            with np.errstate(invalid="ignore"), _quiet():
                sd = np.nanstd(beat_vals, axis=0, ddof=1) if len(idx) > 1 else np.full(len(features), np.nan)
                n_finite = np.isfinite(beat_vals).sum(axis=0)
                means = np.nanmean(beat_vals, axis=0)
            for i, f in enumerate(features):
                row[f] = means[i] if i in seq_cols else feats.get(f, np.nan)
                row[f"{f}_ci95"] = 1.96 * sd[i] / np.sqrt(n_finite[i]) if n_finite[i] > 1 else np.nan
            rows.append(row)
            templates.append(tpl)
        table = pd.DataFrame(rows)
        if "perfusion_index" in features and len(table):
            dc = table["dc_level"].to_numpy()
            table["perfusion_index"] = np.where(dc > 0, 100.0 * table["peak_amplitude"] / np.where(dc > 0, dc, 1),
                                                np.nan)
        return WindowSet(name, features, table, templates, anchor)

    def fixed_windows(self, ppg: ProcessedPPG, window_sec: float, overlap_frac: float,
                      anchor: str) -> List[Tuple[float, float, float, np.ndarray, Dict]]:
        """Fixed-duration windows with the same anchor semantics as the adaptive engine.

        Anchors are spaced by ``window_sec * (1 - overlap_frac)``; the last
        systolic peak is always an anchor, so the tail is covered.
        """
        fs = ppg.sampling_rate
        t_peak = ppg.beat_times
        valid = np.flatnonzero(ppg.valid_beats_mask)
        if len(valid) == 0:
            return []
        step = max(window_sec * (1.0 - overlap_frac), 1.0 / fs)
        t_first, t_last = t_peak[0], t_peak[-1]
        start = t_first + (window_sec if anchor == "trailing" else window_sec / 2)
        anchors = list(np.arange(start, t_last, step))
        if not anchors or anchors[-1] < t_last:
            anchors.append(t_last)
        out = []
        for ta in anchors:
            lo, hi = (ta - window_sec, ta) if anchor == "trailing" else (ta - window_sec / 2, ta + window_sec / 2)
            idx = valid[(t_peak[valid] > lo) & (t_peak[valid] <= hi)]
            if len(idx) < self.t_crit_fixed:
                continue
            out.append((float(ta), float(ppg.onset_positions[idx[0]] / fs), float(ppg.end_positions[idx[-1]] / fs),
                        idx, {}))
        return out

    # -- pipeline ------------------------------------------------------ #
    def process(self, ppg: ProcessedPPG, adaptive: AdaptiveResult,
                fixed_window_secs: Optional[Sequence[float]] = None) -> FeatureResult:
        cfg = self.config
        fs = ppg.sampling_rate
        builder = TemplateBuilder(ppg, cfg)

        r_peaks = None
        pat = (np.full(ppg.n_beats, np.nan), np.full(ppg.n_beats, np.nan))
        if "ecg" in ppg.aux_signals:
            r_peaks = detect_r_peaks(ppg.aux_signals["ecg"], fs)
            pat = pulse_arrival_times(ppg, r_peaks)
        dc = np.array([np.mean(ppg.raw_signal[o : e + 1]) for o, e in zip(ppg.onsets, ppg.ends)])

        per_beat = self.per_beat(ppg, builder, pat, dc)

        adaptive_sets: Dict[str, WindowSet] = {}
        for cat, feats in FEATURE_CATEGORIES.items():
            if cat not in adaptive.schedules:
                continue
            wins = [(w.t_anchor, w.t_start, w.t_end, w.beat_indices,
                     {"anchor_beat": w.anchor_beat, "T_beats": w.T_beats, "overlap": w.overlap, "H": w.H})
                    for w in adaptive.schedules[cat].windows]
            adaptive_sets[cat] = self.window_table(cat, wins, feats, builder, per_beat, adaptive.anchor)

        grid_t = self._grid(ppg)
        grid = self._to_grid(adaptive_sets, grid_t)

        fixed_sets: Dict[str, WindowSet] = {}
        secs = list(fixed_window_secs) if fixed_window_secs is not None else [cfg.fixed_window_sec]
        for sec in secs:
            wins = self.fixed_windows(ppg, sec, cfg.fixed_overlap_frac, adaptive.anchor)
            fixed_sets[f"fixed_{sec:g}s"] = self.window_table(f"fixed_{sec:g}s", wins, ALL_FEATURES, builder,
                                                              per_beat, adaptive.anchor)
        fixed_grid = self._to_grid(
            {k: v for k, v in fixed_sets.items() if k == f"fixed_{secs[0]:g}s"}, grid_t
        ) if secs else pd.DataFrame({"time": grid_t})
        return FeatureResult(per_beat, adaptive_sets, grid, fixed_sets, fixed_grid, r_peaks)

    def _grid(self, ppg: ProcessedPPG) -> np.ndarray:
        if ppg.n_beats == 0:
            return np.array([])
        step = 1.0 / self.config.grid_hz
        return np.arange(np.ceil(ppg.beat_times[0] / step) * step, ppg.beat_times[-1] + 1e-9, step)

    @staticmethod
    def _to_grid(sets: Dict[str, WindowSet], grid_t: np.ndarray) -> pd.DataFrame:
        cols: Dict[str, np.ndarray] = {"time": grid_t}
        for ws in sets.values():
            if ws.table.empty:
                for f in ws.features:
                    cols[f] = np.full(len(grid_t), np.nan)
                continue
            t = ws.table["t_anchor"].to_numpy()
            for f in ws.features:
                cols[f] = series_on_grid(t, ws.table[f].to_numpy(dtype=float), grid_t, ws.anchor)
        return pd.DataFrame(cols)


class _quiet:
    """Silences numpy 'mean of empty slice' warnings inside a block."""

    def __enter__(self):
        import warnings

        self._ctx = warnings.catch_warnings()
        self._ctx.__enter__()
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return self

    def __exit__(self, *exc):
        return self._ctx.__exit__(*exc)


__all__ = [
    "ALL_FEATURES",
    "FEATURE_CATEGORIES",
    "FEATURE_TO_CATEGORY",
    "FEATURE_UNITS",
    "FeatureResult",
    "PPGFeatureExtractor",
    "Template",
    "TemplateBuilder",
    "WindowSet",
    "Window",
    "detect_r_peaks",
    "find_dicrotic",
    "pulse_arrival_times",
    "sdptg_points",
    "series_on_grid",
    "template_features",
]
