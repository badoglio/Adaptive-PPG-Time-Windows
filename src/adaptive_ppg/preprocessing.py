"""PPG preprocessing: filtering, beat detection, segmentation and signal quality.

Two zero-phase filter branches are used:

* a **detection** branch (default 0.5-8 Hz) that feeds the Elgendi et al. (2013)
  systolic peak detector, where a narrow band maximises detection robustness;
* a **morphology** branch (default 0.5-30 Hz) that preserves the dicrotic notch
  and the high-frequency content needed by the first/second derivatives.

A beat is the segment between two consecutive onsets (feet). Each beat is
assigned a signal-quality flag combining physiological rate limits, IBI jumps,
correlation with a trailing median template, amplitude outliers and clipping.
The per-beat checks only look at *past* beats, so the pipeline stays causal.

Reference: M. Elgendi et al., "Systolic Peak Detection in Acceleration
Photoplethysmograms Measured from Emergency Responders in Tropical Conditions",
PLoS ONE 8(10): e76585, 2013.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from fractions import Fraction
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator
from scipy.ndimage import uniform_filter1d
from scipy.signal import butter, resample_poly, sosfiltfilt

from ._typing import BoolArray, FloatArray, IntArray
from .config import PreprocessingConfig
from .io import SignalData

logger = logging.getLogger(__name__)


@dataclass
class ProcessedPPG:
    """Output of :class:`PPGPreprocessor`.

    Beat ``i`` spans samples ``onsets[i]`` (inclusive) to ``ends[i]`` (the next
    onset, inclusive) and has its systolic peak at ``peaks[i]``.

    Attributes:
        raw_signal: Input signal (NaN-free), resampled to ``sampling_rate`` but with
            its original polarity and DC level (used for DC level and perfusion index).
        filtered_signal: Morphology-branch signal (used for features), upright.
        detection_signal: Detection-branch signal (used for peak detection only), upright.
        sampling_rate: Sampling frequency of all arrays, in Hz (after resampling).
        onsets, peaks, ends: Per-beat sample indices (int arrays of equal length).
        beat_times: Time of each systolic peak in seconds (sub-sample precision).
        beats_normalized: ``(n_beats, L)`` beats resampled with PCHIP, baseline-removed
            and scaled to a unit peak. Used only for the variability metric and SQI.
        valid_beats_mask: Boolean quality flag per beat.
        beat_quality: Per-beat table with every SQI component.
        effective_highcut: Morphology low-pass cutoff actually applied (Hz).
        warnings: Human-readable warnings raised during processing.
        aux_signals: Auxiliary channels passed through from the loader (resampled).
        original_sampling_rate: Sampling frequency of the input, in Hz.
        polarity: +1 if the input was upright, -1 if it was inverted before filtering.
        polarity_skewness: Skewness of the derivative of the detection-band input,
            the statistic used by ``polarity="auto"``.
        onset_positions, peak_positions, end_positions: Sub-sample fiducial
            positions (float sample indices): fractional tangent onset and
            parabolic peak. They keep timing features independent of the
            sampling rate; default to the integer indices.
    """

    raw_signal: FloatArray
    filtered_signal: FloatArray
    detection_signal: FloatArray
    sampling_rate: float
    onsets: IntArray
    peaks: IntArray
    ends: IntArray
    beat_times: FloatArray
    beats_normalized: FloatArray
    valid_beats_mask: BoolArray
    beat_quality: pd.DataFrame
    effective_highcut: float
    warnings: List[str] = field(default_factory=list)
    aux_signals: Dict[str, FloatArray] = field(default_factory=dict)
    original_sampling_rate: float = float("nan")
    polarity: int = 1
    polarity_skewness: float = float("nan")
    # None -> filled from the integer indices in __post_init__
    onset_positions: FloatArray = None  # type: ignore[assignment]
    peak_positions: FloatArray = None  # type: ignore[assignment]
    end_positions: FloatArray = None  # type: ignore[assignment]

    def __post_init__(self) -> None:
        if not np.isfinite(self.original_sampling_rate):
            self.original_sampling_rate = float(self.sampling_rate)
        if self.onset_positions is None:
            self.onset_positions = np.asarray(self.onsets, dtype=float)
        if self.peak_positions is None:
            self.peak_positions = np.asarray(self.peaks, dtype=float)
        if self.end_positions is None:
            self.end_positions = np.asarray(self.ends, dtype=float)

    @property
    def n_beats(self) -> int:
        return int(len(self.peaks))

    @property
    def time(self) -> FloatArray:
        return np.arange(len(self.raw_signal)) / self.sampling_rate

    @property
    def beat_durations(self) -> FloatArray:
        return (self.end_positions - self.onset_positions) / self.sampling_rate


# ---------------------------------------------------------------------- #
# Sampling-rate normalization and polarity
# ---------------------------------------------------------------------- #
def resample_signal(
    x: FloatArray, fs: float, target_fs: Optional[float], tol: float = 0.01, max_denominator: int = 1000
) -> Tuple[FloatArray, float]:
    """Resamples ``x`` from ``fs`` to (approximately) ``target_fs`` with a polyphase FIR.

    The ratio is approximated by a fraction ``up/down`` with ``down <=
    max_denominator``; the returned rate is the exact ``fs * up / down``, so
    the time base stays exact even when it differs slightly from
    ``target_fs``. ``resample_poly`` applies its own anti-aliasing low-pass;
    the ends are padded with a line fit to avoid edge transients on signals
    with a large DC offset (raw optical counts). Rates within ``tol`` of the
    target, and ``target_fs=None``, leave the signal untouched.
    """
    if target_fs is None or abs(fs / target_fs - 1.0) <= tol:
        return np.asarray(x, dtype=np.float64), float(fs)
    frac = Fraction(target_fs / fs).limit_denominator(max_denominator)
    up, down = frac.numerator, frac.denominator
    y = resample_poly(np.asarray(x, dtype=np.float64), up, down, padtype="line")
    return y, float(fs) * up / down


def _resample_mask(mask: BoolArray, n_out: int) -> BoolArray:
    """Maps a boolean sample mask onto a grid of ``n_out`` samples spanning the same time."""
    if len(mask) == 0 or n_out == len(mask) or not mask.any():
        return np.resize(mask, n_out) if len(mask) != n_out else mask
    src = np.linspace(0.0, len(mask) - 1.0, n_out)
    return np.interp(src, np.arange(len(mask)), mask.astype(np.float64)) > 0


def derivative_skewness(x: FloatArray) -> float:
    """Skewness of the first difference of ``x``.

    A pulse wave rises quickly and decays slowly, so the derivative of an
    upright PPG has large, brief positive values and a positive skewness; raw
    optical counts from reflectance sensors, which fall as blood volume rises,
    give a negative one.
    """
    d = np.diff(np.asarray(x, dtype=np.float64))
    sd = d.std()
    if len(d) < 3 or sd == 0:
        return 0.0
    return float(np.mean(((d - d.mean()) / sd) ** 3))


# ---------------------------------------------------------------------- #
# Filtering
# ---------------------------------------------------------------------- #
def bandpass(signal: FloatArray, fs: float, lowcut: float, highcut: float, order: int) -> FloatArray:
    """Zero-phase Butterworth band-pass in second-order sections."""
    nyq = 0.5 * fs
    sos = butter(order, [lowcut / nyq, highcut / nyq], btype="band", output="sos")
    return sosfiltfilt(sos, signal - np.median(signal))


def _safe_highcut(highcut: float, fs: float, what: str, warnings: List[str]) -> float:
    limit = 0.45 * fs
    if highcut >= limit:
        msg = f"{what} high cutoff {highcut:g} Hz is too close to Nyquist ({fs / 2:g} Hz); using {limit:.2f} Hz."
        logger.warning(msg)
        warnings.append(msg)
        return limit
    return highcut


# ---------------------------------------------------------------------- #
# Beat detection
# ---------------------------------------------------------------------- #
def elgendi_peaks(
    x: np.ndarray,
    fs: float,
    peakwindow: float = 0.111,
    beatwindow: float = 0.667,
    beatoffset: float = 0.02,
    mindelay: float = 0.3,
) -> np.ndarray:
    """Systolic peak detector of Elgendi et al. (2013), two moving-average version.

    Matches the parameterisation of ``biosppy.signals.ppg.find_onsets_elgendi2013``
    (which, despite its name, returns systolic peaks).
    """
    y = np.clip(x, 0.0, None) ** 2
    w_peak = max(int(round(peakwindow * fs)), 1)
    w_beat = max(int(round(beatwindow * fs)), 1)
    ma_peak = uniform_filter1d(y, size=w_peak, mode="nearest")
    ma_beat = uniform_filter1d(y, size=w_beat, mode="nearest")
    threshold = ma_beat + beatoffset * np.mean(y)

    blocks = (ma_peak > threshold).astype(np.int8)
    edges = np.diff(np.concatenate(([0], blocks, [0])))
    starts, stops = np.flatnonzero(edges == 1), np.flatnonzero(edges == -1)

    peaks: List[int] = []
    min_dist = int(round(mindelay * fs))
    for a, b in zip(starts, stops):
        if b - a < w_peak:
            continue
        p = a + int(np.argmax(x[a:b]))
        if peaks and p - peaks[-1] < min_dist:
            if x[p] > x[peaks[-1]]:
                peaks[-1] = p
            continue
        peaks.append(p)
    return np.asarray(peaks, dtype=int)


def _biosppy_peaks(x: np.ndarray, fs: float) -> np.ndarray:
    from biosppy.signals.ppg import find_onsets_elgendi2013  # type: ignore

    out = find_onsets_elgendi2013(signal=x, sampling_rate=fs)
    return np.asarray(out[0], dtype=int)  # biosppy calls them "onsets" but they are systolic peaks


def _round_index(pos: FloatArray) -> IntArray:
    return np.rint(pos).astype(int)


def _refine_to_local_max(y: np.ndarray, idx: np.ndarray, half_width: int) -> np.ndarray:
    out = np.empty_like(idx)
    n = len(y)
    for k, p in enumerate(idx):
        a, b = max(p - half_width, 0), min(p + half_width + 1, n)
        out[k] = a + int(np.argmax(y[a:b]))
    return out


def parabolic_offset(y: FloatArray, i: int) -> float:
    """Sub-sample offset (in [-0.5, 0.5]) of the extremum at ``y[i]`` from a parabola through 3 samples."""
    if i <= 0 or i >= len(y) - 1:
        return 0.0
    a, b, c = float(y[i - 1]), float(y[i]), float(y[i + 1])
    den = a - 2.0 * b + c
    if den == 0 or not np.isfinite(den):
        return 0.0
    return float(np.clip(0.5 * (a - c) / den, -0.5, 0.5))


def _tangent_onset(y: np.ndarray, foot: int, peak: int) -> float:
    """Intersection of the max-slope tangent with the horizontal through the foot.

    Returns a fractional sample position: the intersection is a continuous
    quantity, so it is not rounded (rounding would quantize every timing
    feature to ``1 / fs``).
    """
    if peak - foot < 3:
        return float(foot)
    dy = np.gradient(y[foot : peak + 1])
    ms = int(np.argmax(dy))
    if dy[ms] <= 0:
        return float(foot)
    t = ms - (y[foot + ms] - y[foot]) / dy[ms]
    return foot + float(np.clip(t, 0, ms))


# ---------------------------------------------------------------------- #
# Resampling
# ---------------------------------------------------------------------- #
def resample_beat(
    segment: np.ndarray, length: int, peak_pos: Optional[int] = None, peak_fraction: float = 0.3
) -> np.ndarray:
    """Resamples one beat to ``length`` points with shape-preserving PCHIP.

    If ``peak_pos`` is given (piecewise mode), the rising edge is mapped onto
    the first ``peak_fraction`` of the output and the decay onto the rest, so
    systolic peaks line up across beats of different duration.
    """
    n = len(segment)
    x = np.arange(n, dtype=float)
    interp = PchipInterpolator(x, segment)
    if peak_pos is None or peak_pos <= 0 or peak_pos >= n - 1:
        return interp(np.linspace(0, n - 1, length))
    n_rise = max(int(round(peak_fraction * (length - 1))), 1)
    rise = np.linspace(0, peak_pos, n_rise + 1)[:-1]
    fall = np.linspace(peak_pos, n - 1, length - n_rise)
    return interp(np.concatenate((rise, fall)))


def normalize_beat(segment: FloatArray, peak_pos: int) -> FloatArray:
    """Removes the onset-to-end linear baseline and scales the peak to 1."""
    n = len(segment)
    baseline = np.linspace(segment[0], segment[-1], n)
    detr = segment - baseline
    scale = detr[peak_pos]
    if not np.isfinite(scale) or abs(scale) < 1e-12:
        scale = np.ptp(detr) or 1.0
    return detr / scale


# ---------------------------------------------------------------------- #
# Preprocessor
# ---------------------------------------------------------------------- #
class PPGPreprocessor:
    """Filters the PPG, detects beats and computes per-beat quality."""

    def __init__(self, config: Optional[PreprocessingConfig] = None, detector: str = "elgendi"):
        self.config = config or PreprocessingConfig()
        if detector not in ("elgendi", "biosppy"):
            raise ValueError("detector must be 'elgendi' or 'biosppy'")
        self.detector = detector

    # -- filtering ------------------------------------------------------ #
    def filter_signals(
        self, signal: FloatArray, fs: float, warnings: List[str], band_limit_fs: Optional[float] = None
    ) -> Tuple[FloatArray, FloatArray, float]:
        """Returns ``(detection, morphology, effective_highcut)``.

        ``band_limit_fs`` is the rate that limits the usable bandwidth (the
        original rate of an upsampled signal); cutoffs are clamped to 45 % of it.
        """
        cfg = self.config
        limit_fs = min(fs, band_limit_fs) if band_limit_fs else fs
        det_high = _safe_highcut(cfg.detection_highcut, limit_fs, "Detection", warnings)
        detection = bandpass(signal, fs, cfg.detection_lowcut, det_high, cfg.filter_order)
        if cfg.bypass_filter:
            morph = signal - np.median(signal)
            highcut = float("nan")
        else:
            highcut = _safe_highcut(cfg.highcut, limit_fs, "Morphology", warnings)
            morph = bandpass(signal, fs, cfg.lowcut, highcut, cfg.filter_order)
        return detection, morph, highcut

    # -- detection ------------------------------------------------------ #
    def detect_beats(
        self, detection: FloatArray, morph: FloatArray, fs: float
    ) -> Tuple[IntArray, IntArray, IntArray]:
        """Returns ``(onsets, peaks, ends)`` of complete beats (integer sample indices)."""
        onsets, peaks, ends = self.detect_beat_positions(detection, morph, fs)
        return _round_index(onsets), _round_index(peaks), _round_index(ends)

    def detect_beat_positions(
        self, detection: FloatArray, morph: FloatArray, fs: float
    ) -> Tuple[FloatArray, FloatArray, FloatArray]:
        """Returns ``(onsets, peaks, ends)`` of complete beats as fractional sample positions."""
        raw_peaks = _biosppy_peaks(detection, fs) if self.detector == "biosppy" else elgendi_peaks(detection, fs)
        if len(raw_peaks) < 2:
            empty = np.array([], dtype=float)
            return empty, empty, empty
        peaks = _refine_to_local_max(morph, raw_peaks, max(int(round(0.05 * fs)), 1))
        peaks = np.unique(peaks)

        # Foot of beat k: minimum of the morphology signal between peak k-1 and peak k.
        ibi = np.diff(peaks)
        lookback = int(np.median(ibi)) if len(ibi) else int(fs)
        prev = np.concatenate(([max(peaks[0] - lookback, 0)], peaks[:-1]))
        feet = np.array([a + int(np.argmin(morph[a:b])) if b > a else a for a, b in zip(prev, peaks)], dtype=int)
        if self.config.onset_method == "tangent":
            feet_pos = np.array([_tangent_onset(morph, f, p) for f, p in zip(feet, peaks)])
        else:
            feet_pos = feet + np.array([parabolic_offset(morph, f) for f in feet])
        peaks_pos = peaks + np.array([parabolic_offset(morph, k) for k in peaks])

        # A beat needs a following onset to close it; drop beats with degenerate rising edges.
        onsets, beat_peaks, ends = feet_pos[:-1], peaks_pos[:-1], feet_pos[1:]
        keep = (beat_peaks - onsets >= 2) & (ends - beat_peaks >= 2)
        return onsets[keep], beat_peaks[keep], ends[keep]

    # -- quality -------------------------------------------------------- #
    def _flat_clipped_mask(self, raw: FloatArray, fs: float) -> BoolArray:
        """Samples belonging to flat runs at the signal extremes (ADC saturation)."""
        n_min = max(int(round(self.config.clip_min_run_sec * fs)), 2)
        rng = np.ptp(raw)
        mask = np.zeros(len(raw), dtype=bool)
        if rng == 0:
            mask[:] = True
            return mask
        flat = np.concatenate(([False], np.diff(raw) == 0))
        edges = np.diff(np.concatenate(([0], flat.astype(np.int8), [0])))
        lo, hi = raw.min() + 0.01 * rng, raw.max() - 0.01 * rng
        for a, b in zip(np.flatnonzero(edges == 1), np.flatnonzero(edges == -1)):
            start = a - 1
            if b - start >= n_min and (raw[start] <= lo or raw[start] >= hi):
                mask[start:b] = True
        return mask

    def assess_quality(
        self,
        clipped: BoolArray,
        morph: FloatArray,
        fs: float,
        onsets: np.ndarray,
        peaks: np.ndarray,
        ends: np.ndarray,
        beats_norm: np.ndarray,
    ) -> pd.DataFrame:
        cfg = self.config
        n = len(peaks)
        durations = (ends - onsets) / fs
        bpm = 60.0 / durations
        amplitude = morph[peaks] - morph[onsets]
        m_ref = cfg.sqi_reference_beats

        finite = np.array([np.all(np.isfinite(morph[a : b + 1])) for a, b in zip(onsets, ends)], dtype=bool)
        rate_ok = (bpm >= cfg.min_bpm) & (bpm <= cfg.max_bpm)

        # Trailing references (past beats only), bootstrapped on the first beats.
        boot = slice(0, min(n, max(5, m_ref)))
        ref_dur, ref_amp = np.empty(n), np.empty(n)
        corr = np.empty(n)
        boot_template = np.median(beats_norm[boot], axis=0) if n else None
        for k in range(n):
            lo = max(0, k - m_ref)
            if k - lo >= 5:
                ref_dur[k] = np.median(durations[lo:k])
                ref_amp[k] = np.median(amplitude[lo:k])
                template = np.median(beats_norm[lo:k], axis=0)
            else:
                ref_dur[k] = np.median(durations[boot])
                ref_amp[k] = np.median(amplitude[boot])
                template = boot_template
            b = beats_norm[k]
            if np.std(b) > 0 and np.std(template) > 0:
                corr[k] = float(np.corrcoef(b, template)[0, 1])
            else:
                corr[k] = 0.0

        ibi_ok = np.abs(durations - ref_dur) <= cfg.max_ibi_jump * ref_dur
        corr_ok = corr >= cfg.min_template_corr
        lo_r, hi_r = cfg.amplitude_ratio_limits
        with np.errstate(divide="ignore", invalid="ignore"):
            amp_ratio = amplitude / ref_amp
        amp_ok = (amp_ratio >= lo_r) & (amp_ratio <= hi_r) & (amplitude > 0)

        clip_ok = np.array([not clipped[a : b + 1].any() for a, b in zip(onsets, ends)], dtype=bool)

        valid = finite & rate_ok & ibi_ok & corr_ok & amp_ok & clip_ok
        return pd.DataFrame(
            {
                "onset_idx": onsets,
                "peak_idx": peaks,
                "end_idx": ends,
                "t_onset": onsets / fs,
                "t_peak": peaks / fs,
                "duration_s": durations,
                "bpm": bpm,
                "amplitude": amplitude,
                "template_corr": corr,
                "amplitude_ratio": amp_ratio,
                "finite_ok": finite,
                "rate_ok": rate_ok,
                "ibi_ok": ibi_ok,
                "corr_ok": corr_ok,
                "amp_ok": amp_ok,
                "clip_ok": clip_ok,
                "valid": valid,
            }
        )

    # -- pipeline ------------------------------------------------------- #
    def process(self, signal_data: SignalData) -> ProcessedPPG:
        cfg = self.config
        fs_in = float(signal_data.sampling_rate)
        raw_in = np.asarray(signal_data.signal, dtype=np.float64)
        warnings: List[str] = []
        if fs_in < 50:
            msg = f"Sampling rate {fs_in:g} Hz is low: derivative features (SDPTG) will be unreliable."
            logger.warning(msg)
            warnings.append(msg)
        if len(raw_in) < 10 * fs_in:
            raise ValueError(f"Signal too short ({len(raw_in) / fs_in:.1f} s); at least 10 s are required.")

        # Saturation is detected at the native rate: resampling smooths flat runs.
        clipped_in = self._flat_clipped_mask(raw_in, fs_in)
        raw, fs = resample_signal(raw_in, fs_in, cfg.target_fs)
        aux = dict(signal_data.aux_signals)
        if fs != fs_in:
            aux = {k: resample_signal(v, fs_in, cfg.target_fs)[0][: len(raw)] for k, v in aux.items()}
            clipped = _resample_mask(clipped_in, len(raw))
            logger.info("Resampled from %g Hz to %g Hz.", fs_in, fs)
        else:
            clipped = clipped_in

        detection, morph, highcut = self.filter_signals(raw, fs, warnings, band_limit_fs=fs_in)
        skewness = derivative_skewness(detection)
        polarity = {"normal": 1, "inverted": -1}.get(cfg.polarity, 1 if skewness >= 0 else -1)
        if polarity < 0:
            detection, morph = -detection, -morph
            if cfg.polarity == "auto":
                msg = f"Signal inverted automatically (derivative skewness {skewness:.2f} < 0)."
                logger.info(msg)
                warnings.append(msg)
        onset_pos, peak_pos, end_pos = self.detect_beat_positions(detection, morph, fs)
        onsets, peaks, ends = _round_index(onset_pos), _round_index(peak_pos), _round_index(end_pos)
        n = len(peaks)
        L = cfg.target_beat_length
        if n == 0:
            warnings.append("No beats detected.")
            empty = np.array([], dtype=int)
            return ProcessedPPG(raw, morph, detection, fs, empty, empty, empty, np.array([]),
                                np.empty((0, L)), np.array([], dtype=bool), pd.DataFrame(), highcut, warnings,
                                aux, fs_in, polarity, skewness)

        beats_norm = np.empty((n, L))
        for k, (a, p, b) in enumerate(zip(onsets, peaks, ends)):
            seg = normalize_beat(morph[a : b + 1], p - a)
            piece = p - a if cfg.resample_mode == "piecewise" else None
            beats_norm[k] = resample_beat(seg, L, piece, cfg.piecewise_peak_fraction)

        quality = self.assess_quality(clipped, morph, fs, onsets, peaks, ends, beats_norm)
        valid = quality["valid"].to_numpy()
        if valid.mean() < 0.5:
            msg = f"Only {100 * valid.mean():.0f}% of beats passed quality checks."
            logger.warning(msg)
            warnings.append(msg)

        return ProcessedPPG(
            raw_signal=raw,
            filtered_signal=morph,
            detection_signal=detection,
            sampling_rate=fs,
            onsets=onsets,
            peaks=peaks,
            ends=ends,
            beat_times=peak_pos / fs,
            beats_normalized=beats_norm,
            valid_beats_mask=valid,
            beat_quality=quality,
            effective_highcut=highcut,
            warnings=warnings,
            aux_signals=aux,
            original_sampling_rate=fs_in,
            polarity=polarity,
            polarity_skewness=skewness,
            onset_positions=onset_pos,
            peak_positions=peak_pos,
            end_positions=end_pos,
        )
