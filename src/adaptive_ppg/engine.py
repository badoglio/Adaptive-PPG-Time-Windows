"""Adaptive window engine.

For every beat the engine computes a **morphological variability** index
``H_morph`` (it is a variability, not an entropy: the mean RMS difference
between consecutive normalized beats, or their dispersion around the local
mean), normalizes it causally to ``[0, 1]`` and maps it to a target window
length and overlap for each feature category::

    s(H)   = 1 / (1 + exp(k * (H - theta)))        (sigmoid)   or   1 - H  (linear)
    T_w^c  = T_min^c + (T_max^c - T_min^c) * s(H)             T_max^c in beats = T_max_sec^c / local IBI
    T_w^c  = max(T_w^c, T_crit^c)                              counted in *valid* beats
    O_w    = O_max - (O_max - O_min) * s(H)                    shorter windows overlap more

Window growth is rate-limited (``max_expand_per_beat``) while contraction is
immediate, so the window reacts at once to a transient but does not oscillate
beat by beat when ``H`` is noisy.

All statistics are computed on the sequence of *valid* beats; invalid beats
inherit the last value. With ``anchor="trailing"`` and a non-global
normalization, every quantity at beat ``k`` depends only on beats ``<= k``;
:class:`StreamingAdaptiveEngine` reproduces the batch result incrementally.
"""

from __future__ import annotations

import collections
from dataclasses import dataclass, field
from typing import Deque, Dict, List, Optional, Tuple

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view

from ._typing import BoolArray, FloatArray, IntArray
from .config import CategoryConfig, EngineConfig
from .preprocessing import ProcessedPPG

IBI_HISTORY_BEATS = 10


@dataclass
class Window:
    """One ensemble-averaging window of a feature category."""

    category: str
    anchor_beat: int
    t_anchor: float
    beat_indices: IntArray
    t_start: float
    t_end: float
    T_target: float
    T_beats: int
    overlap: float
    H: float

    @property
    def n_valid(self) -> int:
        return int(len(self.beat_indices))

    def to_dict(self) -> Dict[str, object]:
        return {
            "category": self.category,
            "anchor_beat": int(self.anchor_beat),
            "t_anchor": float(self.t_anchor),
            "t_start": float(self.t_start),
            "t_end": float(self.t_end),
            "n_valid": self.n_valid,
            "T_target": float(self.T_target),
            "T_beats": int(self.T_beats),
            "overlap": float(self.overlap),
            "H": float(self.H),
            "beat_indices": [int(i) for i in self.beat_indices],
        }


@dataclass
class CategorySchedule:
    """Per-beat window curves and the resulting windows of one category."""

    name: str
    config: CategoryConfig
    T_max_beats: FloatArray
    T_target: FloatArray
    T_beats: IntArray
    windows: List[Window] = field(default_factory=list)


@dataclass
class AdaptiveResult:
    """Output of :class:`AdaptiveEngine`. All per-beat arrays have length ``n_beats``."""

    beat_times: FloatArray
    valid_mask: BoolArray
    variability_raw: FloatArray
    ibi_variability: FloatArray
    variability_used: FloatArray
    H: FloatArray
    overlap: FloatArray
    norm_lo: FloatArray
    norm_hi: FloatArray
    warmup_mask: BoolArray
    local_ibi: FloatArray
    ibi_slope: float
    anchor: str
    schedules: Dict[str, CategorySchedule]

    @property
    def windows(self) -> Dict[str, List[Window]]:
        return {name: s.windows for name, s in self.schedules.items()}

    @property
    def categories(self) -> List[str]:
        return list(self.schedules)


# ---------------------------------------------------------------------- #
# Building blocks (shared by the batch and the streaming engine)
# ---------------------------------------------------------------------- #
def _rms_rows(x: np.ndarray) -> np.ndarray:
    return np.sqrt(np.mean(x * x, axis=-1))


def trailing_nanmean(x: FloatArray, n: int) -> FloatArray:
    """Mean of the last ``n`` finite values up to and including each index."""
    padded = np.concatenate((np.full(n - 1, np.nan), x))
    view = sliding_window_view(padded, n)
    counts = np.isfinite(view).sum(axis=1)
    sums = np.nansum(view, axis=1)
    out = np.full(len(x), np.nan)
    ok = counts > 0
    out[ok] = sums[ok] / counts[ok]
    return out


def successive_differences(beats: FloatArray) -> FloatArray:
    """RMS difference between each beat and the previous one (NaN for the first)."""
    d = np.full(len(beats), np.nan)
    if len(beats) > 1:
        d[1:] = _rms_rows(beats[1:] - beats[:-1])
    return d


def dispersion(beats: FloatArray, n: int) -> FloatArray:
    """Mean RMS deviation of the last ``n`` beats from their mean template."""
    out = np.full(len(beats), np.nan)
    for j in range(1, min(n - 1, len(beats))):
        blk = beats[: j + 1]
        out[j] = float(np.mean(_rms_rows(blk - blk.mean(axis=0))))
    if len(beats) >= n:
        view = sliding_window_view(beats, n, axis=0)  # (m, L, n)
        dev = view - view.mean(axis=2, keepdims=True)
        out[n - 1 :] = np.mean(np.sqrt(np.mean(dev * dev, axis=1)), axis=1)
    return out


def drift(beats: FloatArray, n: int, standardized: bool = False) -> FloatArray:
    """RMS difference between the mean template of the last ``n`` beats and of the ``n`` before them.

    Beat-to-beat noise and respiratory modulation largely cancel in the two
    means (their contribution shrinks as ``sqrt(2 / n)``), while a sustained
    change of morphology appears with its full size: the metric measures
    non-stationarity rather than beat-to-beat scatter. NaN until ``2 n``
    beats are available.

    With ``standardized=True`` the value is divided by its expectation under
    stationarity, ``sqrt(2 / n)`` times the pooled within-block RMS
    deviation: it is ~1 for a stationary signal *whatever its noise level*
    and grows with the size of a morphology change relative to the noise
    (``drift_z`` metric).
    """
    out = np.full(len(beats), np.nan)
    for j in range(2 * n - 1, len(beats)):
        out[j] = _drift_value(beats[j - 2 * n + 1 : j + 1], n, standardized)
    return out


def _drift_value(block: np.ndarray, n: int, standardized: bool = False) -> float:
    recent, previous = block[n:], block[:n]
    m_recent, m_previous = np.mean(recent, axis=0), np.mean(previous, axis=0)
    value = float(_rms_rows(m_recent - m_previous))
    if not standardized:
        return value
    pooled = (np.sum((recent - m_recent) ** 2, axis=0) + np.sum((previous - m_previous) ** 2, axis=0)) / (2 * n - 2)
    noise = float(np.sqrt(np.mean(pooled) * 2.0 / n))
    return value / noise if noise > 1e-12 else np.nan


def shape_factor(H: FloatArray, cfg: EngineConfig) -> FloatArray:
    """Maps ``H`` in [0, 1] to a window-length factor in [0, 1] (1 = longest window)."""
    if cfg.mapping == "linear":
        return 1.0 - H
    z = np.clip(cfg.sigmoid_slope * (H - cfg.theta), -700, 700)
    return 1.0 / (1.0 + np.exp(z))


class RunningIBIRegression:
    """Expanding (causal) OLS of the variability on the IBI variability.

    ``correct(v, x)`` first adds the pair to the running sums and then returns
    ``v - slope * (x - mean_x)`` with the slope fitted on all pairs seen so far,
    i.e. the part of the morphological variability not explained by heart-rate
    changes. Until ``min_samples`` pairs are available the value is returned
    unchanged. Batch and streaming engines share this class, so both produce
    bit-identical results.
    """

    def __init__(self, min_samples: int = 30):
        self.min_samples = min_samples
        self.n = 0
        self.sx = self.sy = self.sxx = self.sxy = 0.0

    @property
    def slope(self) -> float:
        if self.n < self.min_samples:
            return 0.0
        mx, my = self.sx / self.n, self.sy / self.n
        var = self.sxx / self.n - mx * mx
        return 0.0 if var <= 1e-15 else (self.sxy / self.n - mx * my) / var

    def correct(self, v: float, x: float) -> float:
        if not (np.isfinite(v) and np.isfinite(x)):
            return v
        self.n += 1
        self.sx += x
        self.sy += v
        self.sxx += x * x
        self.sxy += x * v
        if self.n < self.min_samples:
            return v
        return v - self.slope * (x - self.sx / self.n)


def _normalize(value: float, lo: float, hi: float) -> float:
    if not (np.isfinite(value) and np.isfinite(lo) and np.isfinite(hi)) or hi - lo <= 1e-12:
        return 0.0
    return float(min(max((value - lo) / (hi - lo), 0.0), 1.0))


def _step(T_beats: int, overlap: float) -> int:
    return max(1, int(round(T_beats * (1.0 - overlap))))


def _rate_limit(target: np.ndarray, max_expand: Optional[float]) -> np.ndarray:
    if max_expand is None or len(target) == 0:
        return target.copy()
    out = np.empty_like(target)
    out[0] = target[0]
    for k in range(1, len(target)):
        out[k] = min(target[k], out[k - 1] + max_expand)
    return out


def _category_curves(
    s: np.ndarray, local_ibi: np.ndarray, c: CategoryConfig, max_expand: Optional[float]
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    with np.errstate(divide="ignore", invalid="ignore"):
        t_max = np.where(local_ibi > 0, c.T_max_sec / local_ibi, c.T_min_beats)
    t_max = np.maximum(t_max, c.T_min_beats)
    target = c.T_min_beats + (t_max - c.T_min_beats) * s
    limited = _rate_limit(target, max_expand)
    t_beats = np.maximum(np.round(limited), c.T_crit_beats).astype(int)
    return t_max, target, t_beats


# ---------------------------------------------------------------------- #
# Batch engine
# ---------------------------------------------------------------------- #
class AdaptiveEngine:
    """Computes the morphological variability and the adaptive window schedules."""

    def __init__(self, config: Optional[EngineConfig] = None):
        self.config = config or EngineConfig()

    # -- variability ----------------------------------------------------- #
    def variability(self, beats_valid: FloatArray) -> FloatArray:
        cfg = self.config
        if cfg.variability_metric == "dispersion":
            return dispersion(beats_valid, cfg.N_past)
        if cfg.variability_metric in ("drift", "drift_z"):
            return drift(beats_valid, cfg.N_past, standardized=cfg.variability_metric == "drift_z")
        return trailing_nanmean(successive_differences(beats_valid), cfg.N_past)

    def ibi_variability(self, durations_valid: FloatArray) -> FloatArray:
        """Trailing mean of |ΔIBI| / IBI between consecutive valid beats."""
        d = np.full(len(durations_valid), np.nan)
        if len(durations_valid) > 1:
            d[1:] = np.abs(np.diff(durations_valid)) / durations_valid[1:]
        return trailing_nanmean(d, self.config.N_past)

    # -- normalization ------------------------------------------------- #
    def _normalize_valid(
        self, raw: np.ndarray, dibi: np.ndarray, t: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, float]:
        """Returns (used, H, lo, hi, warmup, slope) on the valid-beat sequence."""
        cfg = self.config
        mode = cfg.normalization
        n = len(raw)
        p_lo, p_hi = cfg.percentiles
        cal = t < cfg.calibration_sec

        used = raw.copy()
        slope = 0.0
        if cfg.ibi_correction:
            reg = RunningIBIRegression()
            used = np.array([reg.correct(float(v), float(x)) for v, x in zip(raw, dibi)])
            slope = reg.slope

        lo = np.full(n, np.nan)
        hi = np.full(n, np.nan)
        warmup = np.zeros(n, bool)
        if mode == "global":
            lo[:], hi[:] = np.nanpercentile(used, [p_lo, p_hi]) if np.isfinite(used).any() else (np.nan, np.nan)
        elif mode == "fixed":
            lo[:], hi[:] = cfg.fixed_scale
        elif mode == "calibration":
            warmup = cal.copy()
            ref = used[cal]
            if np.isfinite(ref).sum() >= 3:
                lo[~cal], hi[~cal] = np.nanpercentile(ref, [p_lo, p_hi])
            else:
                lo[~cal], hi[~cal] = cfg.fixed_scale
        else:  # rolling
            warmup = cal.copy()
            for j in np.flatnonzero(~cal):
                buf = used[max(0, j - cfg.rolling_beats + 1) : j + 1]
                if np.isfinite(buf).sum() >= 3:
                    lo[j], hi[j] = np.nanpercentile(buf, [p_lo, p_hi])

        warmup |= ~np.isfinite(used)  # not enough history for the metric yet
        H = np.array([0.0 if warmup[j] else _normalize(used[j], lo[j], hi[j]) for j in range(n)])
        return used, H, lo, hi, warmup, slope

    # -- windows ------------------------------------------------------- #
    @staticmethod
    def _trailing_window(valid_idx: np.ndarray, cum_valid: np.ndarray, k: int, T: int) -> np.ndarray:
        count = int(cum_valid[k])
        take = min(T, count)
        return valid_idx[count - take : count]

    @staticmethod
    def _centered_window(valid_idx: np.ndarray, cum_valid: np.ndarray, k: int, T: int) -> np.ndarray:
        count = int(cum_valid[k])  # valid beats with index <= k
        before, after = count, len(valid_idx) - count
        n_before = min(before, (T + 1) // 2)
        n_after = min(after, T - n_before)
        n_before = min(before, T - n_after)
        return valid_idx[count - n_before : count + n_after]

    def build_windows(
        self,
        name: str,
        c: CategoryConfig,
        T_target: np.ndarray,
        T_beats: np.ndarray,
        overlap: np.ndarray,
        H: np.ndarray,
        valid_mask: np.ndarray,
        onsets_t: np.ndarray,
        ends_t: np.ndarray,
        peak_t: np.ndarray,
    ) -> List[Window]:
        """Places windows at anchor beats; every window has at least ``T_crit`` valid beats.

        The step between anchors is ``round(T_w * (1 - O_w))`` beats. The last
        beat is always used as a final anchor, so the tail of the recording is
        covered and no truncated window is produced.
        """
        n = len(valid_mask)
        valid_idx = np.flatnonzero(valid_mask)
        cum_valid = np.cumsum(valid_mask)
        pick = self._centered_window if self.config.anchor == "centered" else self._trailing_window
        windows: List[Window] = []
        k = 0
        while k < n:
            idx = pick(valid_idx, cum_valid, k, int(T_beats[k]))
            if len(idx) < c.T_crit_beats:
                k += 1
                continue
            windows.append(
                Window(
                    category=name,
                    anchor_beat=int(k),
                    t_anchor=float(peak_t[k]),
                    beat_indices=idx.astype(int),
                    t_start=float(onsets_t[idx[0]]),
                    t_end=float(ends_t[idx[-1]]),
                    T_target=float(T_target[k]),
                    T_beats=int(T_beats[k]),
                    overlap=float(overlap[k]),
                    H=float(H[k]),
                )
            )
            if k == n - 1:
                break
            k = min(k + _step(int(T_beats[k]), float(overlap[k])), n - 1)
        return windows

    # -- pipeline ------------------------------------------------------ #
    def process(self, ppg: ProcessedPPG) -> AdaptiveResult:
        cfg = self.config
        n = ppg.n_beats
        valid = np.asarray(ppg.valid_beats_mask, dtype=bool)
        valid_idx = np.flatnonzero(valid)
        fs = ppg.sampling_rate
        durations = ppg.beat_durations
        t_peak = ppg.beat_times

        raw_v = self.variability(ppg.beats_normalized[valid_idx])
        dibi_v = self.ibi_variability(durations[valid_idx])
        used_v, H_v, lo_v, hi_v, warm_v, slope = self._normalize_valid(raw_v, dibi_v, t_peak[valid_idx])
        ibi_v = np.array([np.median(durations[valid_idx][max(0, j - IBI_HISTORY_BEATS + 1) : j + 1])
                          for j in range(len(valid_idx))])

        # Map valid-beat quantities onto all beats (invalid beats inherit the last value).
        pos = np.cumsum(valid) - 1
        def ffill(values: np.ndarray, default: float) -> np.ndarray:
            out = np.full(n, default, dtype=float)
            has = pos >= 0
            out[has] = values[pos[has]]
            return out

        raw = ffill(raw_v, np.nan)
        dibi = ffill(dibi_v, np.nan)
        used = ffill(used_v, np.nan)
        H = ffill(H_v, 0.0)
        lo = ffill(lo_v, np.nan)
        hi = ffill(hi_v, np.nan)
        warmup = ffill(warm_v.astype(float), 1.0).astype(bool)
        local_ibi = ffill(ibi_v, np.nan)  # NaN before the first valid beat -> T_max = T_min

        if cfg.anchor == "centered":
            # Offline mode: evaluate H and IBI in the middle of the N_past history.
            shift = cfg.N_past // 2
            src = np.minimum(np.arange(n) + shift, n - 1) if n else np.arange(0)
            H, local_ibi = H[src], local_ibi[src]

        s = shape_factor(H, cfg)
        overlap = cfg.O_max - (cfg.O_max - cfg.O_min) * s

        schedules: Dict[str, CategorySchedule] = {}
        for name, c in cfg.categories.items():
            t_max, target, t_beats = _category_curves(s, local_ibi, c, cfg.max_expand_per_beat)
            windows = self.build_windows(name, c, target, t_beats, overlap, H, valid,
                                         ppg.onset_positions / fs, ppg.end_positions / fs, t_peak)
            schedules[name] = CategorySchedule(name, c, t_max, target, t_beats, windows)

        return AdaptiveResult(
            beat_times=t_peak,
            valid_mask=valid,
            variability_raw=raw,
            ibi_variability=dibi,
            variability_used=used,
            H=H,
            overlap=overlap,
            norm_lo=lo,
            norm_hi=hi,
            warmup_mask=warmup,
            local_ibi=local_ibi,
            ibi_slope=slope,
            anchor=cfg.anchor,
            schedules=schedules,
        )


# ---------------------------------------------------------------------- #
# Streaming engine
# ---------------------------------------------------------------------- #
class StreamingAdaptiveEngine:
    """Incremental (online) version of :class:`AdaptiveEngine`.

    Feed beats one at a time with :meth:`push`; each call returns the windows
    that close at that beat. Call :meth:`finalize` at the end of the stream to
    emit the last window. Only causal settings are supported
    (``anchor="trailing"`` and a non-global normalization); with those, the
    windows are identical to the batch engine's.
    """

    def __init__(self, config: Optional[EngineConfig] = None):
        cfg = config or EngineConfig()
        if cfg.anchor != "trailing":
            raise ValueError("StreamingAdaptiveEngine requires anchor='trailing'.")
        if cfg.normalization == "global":
            raise ValueError("StreamingAdaptiveEngine does not support global normalization.")
        self.config = cfg
        n = cfg.N_past
        self._k = -1
        is_drift = cfg.variability_metric in ("drift", "drift_z")
        self._beats: Deque[np.ndarray] = collections.deque(maxlen=2 * n if is_drift else n)
        self._diffs: Deque[float] = collections.deque([np.nan] * (n - 1), maxlen=n)
        self._dibis: Deque[float] = collections.deque([np.nan] * (n - 1), maxlen=n)
        self._durations: Deque[float] = collections.deque(maxlen=IBI_HISTORY_BEATS)
        self._history: List[Tuple[float, float]] = []  # (used variability, t_peak) of valid beats
        self._valid_beats: List[Tuple[int, float, float]] = []  # (beat index, t_onset, t_end)
        self._regression = RunningIBIRegression() if cfg.ibi_correction else None
        self._cal_bounds: Optional[Tuple[float, float]] = None
        self._H = 0.0
        self._ibi: Optional[float] = None
        self._limited: Dict[str, float] = {}
        self._next_anchor = {name: 0 for name in cfg.categories}
        self._last_emitted = {name: -1 for name in cfg.categories}
        self._last_state: Dict[str, Tuple[float, int, float]] = {}
        self._last_t_peak = 0.0

    # -- internal ------------------------------------------------------ #
    def _close_calibration(self) -> None:
        cfg = self.config
        ref = np.array([u for u, t in self._history if t < cfg.calibration_sec])
        if np.isfinite(ref).sum() >= 3:
            lo, hi = np.nanpercentile(ref, cfg.percentiles)
            self._cal_bounds = (float(lo), float(hi))
        else:
            self._cal_bounds = (float(cfg.fixed_scale[0]), float(cfg.fixed_scale[1]))

    def _update_valid(self, beat: np.ndarray, duration: float, t_peak: float) -> None:
        cfg = self.config
        if cfg.variability_metric in ("drift", "drift_z"):
            self._beats.append(beat)
            full = len(self._beats) == self._beats.maxlen
            standardized = cfg.variability_metric == "drift_z"
            raw = _drift_value(np.array(self._beats), cfg.N_past, standardized) if full else np.nan
        elif cfg.variability_metric == "dispersion":
            if len(self._beats) == self._beats.maxlen:
                self._beats.popleft()
            self._beats.append(beat)
            blk = np.array(self._beats)
            raw = float(np.mean(_rms_rows(blk - blk.mean(axis=0)))) if len(blk) > 1 else np.nan
        else:
            d = float(_rms_rows(beat - self._beats[-1])) if self._beats else np.nan
            self._beats.append(beat)
            self._diffs.append(d)
            arr = np.array(self._diffs)
            raw = float(np.nansum(arr) / np.isfinite(arr).sum()) if np.isfinite(arr).any() else np.nan
        prev = self._durations[-1] if self._durations else None
        self._dibis.append(abs(duration - prev) / duration if prev is not None else np.nan)
        arr = np.array(self._dibis)
        dibi = float(np.nansum(arr) / np.isfinite(arr).sum()) if np.isfinite(arr).any() else np.nan
        self._durations.append(duration)
        self._ibi = float(np.median(np.array(self._durations)))
        used = self._regression.correct(raw, dibi) if self._regression is not None else raw
        self._history.append((used, t_peak))

        in_cal = t_peak < cfg.calibration_sec
        if cfg.normalization == "fixed":
            self._H = _normalize(used, *cfg.fixed_scale)
        elif in_cal:
            self._H = 0.0
        elif cfg.normalization == "calibration":
            if self._cal_bounds is None:
                self._close_calibration()
            self._H = _normalize(used, *self._cal_bounds)
        else:  # rolling
            buf = np.array([u for u, _ in self._history[-cfg.rolling_beats:]])
            if np.isfinite(buf).sum() >= 3:
                lo, hi = np.nanpercentile(buf, cfg.percentiles)
                self._H = _normalize(used, float(lo), float(hi))
            else:
                self._H = 0.0

    def _window(self, name: str, k: int, t_peak: float) -> Optional["Window"]:
        target, t_beats, overlap = self._last_state[name]
        take = min(t_beats, len(self._valid_beats))
        if take < self.config.categories[name].T_crit_beats:
            return None
        sel = self._valid_beats[len(self._valid_beats) - take :]
        return Window(name, k, t_peak, np.array([b[0] for b in sel], dtype=int), sel[0][1], sel[-1][2],
                      target, t_beats, overlap, self._H)

    # -- public API ---------------------------------------------------- #
    def push(self, beat_normalized: FloatArray, t_onset: float, t_peak: float, t_end: float,
             valid: bool) -> List[Window]:
        """Adds one beat and returns the windows anchored at it."""
        cfg = self.config
        self._k += 1
        k = self._k
        self._last_t_peak = t_peak
        if valid:
            self._valid_beats.append((k, t_onset, t_end))
            self._update_valid(np.asarray(beat_normalized, dtype=float), t_end - t_onset, t_peak)

        s = float(shape_factor(np.array([self._H]), cfg)[0])
        overlap = cfg.O_max - (cfg.O_max - cfg.O_min) * s
        ibi = self._ibi
        out: List[Window] = []
        for name, c in cfg.categories.items():
            t_max = c.T_max_sec / ibi if ibi and ibi > 0 else c.T_min_beats
            t_max = max(t_max, c.T_min_beats)
            target = c.T_min_beats + (t_max - c.T_min_beats) * s
            prev = self._limited.get(name)
            limited = target if prev is None or cfg.max_expand_per_beat is None else \
                min(target, prev + cfg.max_expand_per_beat)
            self._limited[name] = limited
            t_beats = max(int(np.round(limited)), c.T_crit_beats)
            self._last_state[name] = (target, t_beats, overlap)
            if k == self._next_anchor[name]:
                w = self._window(name, k, t_peak)
                if w is None:
                    self._next_anchor[name] = k + 1
                else:
                    out.append(w)
                    self._last_emitted[name] = k
                    self._next_anchor[name] = k + _step(t_beats, overlap)
        return out

    def finalize(self) -> List[Window]:
        """Emits the closing window at the last beat, as the batch engine does."""
        out: List[Window] = []
        k = self._k
        for name in self.config.categories:
            if k >= 0 and self._next_anchor[name] > k and self._last_emitted[name] < k:
                w = self._window(name, k, self._last_t_peak)
                if w is not None:
                    out.append(w)
                    self._last_emitted[name] = k
        return out

    @property
    def H(self) -> float:
        return self._H
