"""Validation against synthetic ground truth.

The ground truth of a synthetic recording is obtained by generating the same
scenario **without** noise, artefacts, respiratory modulation or heart-rate
jitter and extracting the per-beat features from it: it is the trajectory of
the underlying morphology, which a windowed estimator should follow.

Estimators compared on a common time grid (all causal when
``anchor="trailing"``):

* ``adaptive``: the variable-width ensemble-averaging pipeline;
* ``fixed_<N>s``: fixed-duration ensemble windows (same anchor rules);
* ``ema_<N>s``: exponential moving average of the valid per-beat features;
* ``per_beat``: the raw valid per-beat features (zero-order hold).

Metrics per feature: RMSE and bias on the evaluation grid, SD of the error
in stationary segments (precision), RMSE normalized by the ground-truth SD
(to aggregate across features) and latency after each transition: time for
the estimate to cover 50 % of the ground-truth change, minus the same time
for the ground truth.
"""

from __future__ import annotations

import copy
import itertools
from dataclasses import dataclass
from typing import Any, Dict, Iterable, List, Optional, Sequence

import numpy as np
import pandas as pd

from ._typing import FloatArray
from .config import PipelineConfig
from .engine import AdaptiveEngine
from .features import ALL_FEATURES, PPGFeatureExtractor, series_on_grid
from .pipeline import PipelineResult, run
from .preprocessing import PPGPreprocessor
from .synthetic import NoiseConfig, SyntheticRecording, generate

DEFAULT_FIXED = (10.0, 30.0, 60.0)
DEFAULT_EMA = (10.0,)


# ---------------------------------------------------------------------- #
# Ground truth
# ---------------------------------------------------------------------- #
def ground_truth(rec: SyntheticRecording, config: Optional[PipelineConfig] = None) -> pd.DataFrame:
    """Per-beat ground-truth features (``t_peak`` + one column per feature).

    The noise-free signal is analysed at the generator's native rate
    (``target_fs`` is ignored), so the truth does not inherit the timing
    resolution of the configuration under test.
    """
    cfg = (config or PipelineConfig()).copy()
    cfg.preprocessing.target_fs = None
    quiet = NoiseConfig(snr_db=np.inf, wander_amp=0.0, resp_am=0.0, resp_hr=0.0, hr_jitter=0.0,
                        artifact_rate_per_min=0.0)
    ideal = generate(rec.segments, fs=rec.fs, noise=quiet, with_ecg="ecg" in rec.signal_data.aux_signals, seed=0)
    ppg = PPGPreprocessor(cfg.preprocessing).process(ideal.signal_data)
    feats = PPGFeatureExtractor(cfg.features).process(ppg, AdaptiveEngine(cfg.engine).process(ppg),
                                                       fixed_window_secs=[])
    return feats.per_beat.drop(columns=["beat_index", "valid"])


def _gt_on_grid(gt: pd.DataFrame, grid: np.ndarray, features: Sequence[str]) -> pd.DataFrame:
    out = {"time": grid}
    t = gt["t_peak"].to_numpy()
    for f in features:
        v = gt[f].to_numpy(dtype=float)
        ok = np.isfinite(v)
        col = np.full(len(grid), np.nan)
        if ok.sum() >= 2:
            inside = (grid >= t[ok][0]) & (grid <= t[ok][-1])
            col[inside] = np.interp(grid[inside], t[ok], v[ok])
        out[f] = col
    return pd.DataFrame(out)


# ---------------------------------------------------------------------- #
# Estimators
# ---------------------------------------------------------------------- #
def ema_series(t: FloatArray, v: FloatArray, tau: float) -> FloatArray:
    """Causal EMA for irregular samples (NaNs are skipped)."""
    out = np.full(len(v), np.nan)
    state, t_prev = np.nan, None
    for i, (ti, vi) in enumerate(zip(t, v)):
        if np.isfinite(vi):
            if not np.isfinite(state):
                state = vi
            else:
                alpha = 1.0 - np.exp(-(ti - t_prev) / tau)
                state = state + alpha * (vi - state)
            t_prev = ti
        out[i] = state
    return out


def estimator_grids(result: PipelineResult, grid: FloatArray, ema_secs: Sequence[float] = DEFAULT_EMA
                    ) -> Dict[str, pd.DataFrame]:
    """Every estimator of ``result`` resampled on ``grid``."""
    fr = result.features
    out = {"adaptive": PPGFeatureExtractor._to_grid(fr.adaptive, grid)}
    for name, ws in fr.fixed.items():
        out[name] = PPGFeatureExtractor._to_grid({name: ws}, grid)
    pb = fr.per_beat[fr.per_beat["valid"]]
    t = pb["t_peak"].to_numpy()
    raw = {"time": grid}
    for f in ALL_FEATURES:
        raw[f] = series_on_grid(t, pb[f].to_numpy(dtype=float), grid, "trailing")
    out["per_beat"] = pd.DataFrame(raw)
    for tau in ema_secs:
        cols = {"time": grid}
        for f in ALL_FEATURES:
            cols[f] = series_on_grid(t, ema_series(t, pb[f].to_numpy(dtype=float), tau), grid, "trailing")
        out[f"ema_{tau:g}s"] = pd.DataFrame(cols)
    return out


# ---------------------------------------------------------------------- #
# Metrics
# ---------------------------------------------------------------------- #
def _crossing_time(t: np.ndarray, x: np.ndarray, start: float, stop: float, level: float, rising: bool) -> float:
    sel = (t >= start) & (t < stop) & np.isfinite(x)
    ts, xs = t[sel], x[sel]
    hit = xs >= level if rising else xs <= level
    return float(ts[np.argmax(hit)]) if hit.any() else np.nan


def feature_metrics(
    est: pd.DataFrame,
    gt: pd.DataFrame,
    rec: SyntheticRecording,
    features: Sequence[str],
    eval_start: float,
    scales: Optional[Dict[str, float]] = None,
) -> pd.DataFrame:
    """Per-feature metrics of one estimator.

    Normalized metrics divide by ``scales[feature]`` (default: the SD of the
    ground truth on the evaluation grid). Pass explicit scales for scenarios
    in which the ground truth is constant (e.g. rest).
    """
    t = gt["time"].to_numpy()
    ev = t >= eval_start
    stat = np.zeros(len(t), bool)
    for a, b in rec.stationary:
        stat |= (t >= max(a, eval_start)) & (t <= b)
    rows = []
    boundaries = list(rec.events) + [t[-1] + 1.0]
    for f in features:
        g = gt[f].to_numpy(dtype=float)
        e = est[f].to_numpy(dtype=float)
        ok = ev & np.isfinite(g) & np.isfinite(e)
        err = e[ok] - g[ok]
        gsd = scales[f] if scales and f in scales else np.nanstd(g[ev])
        row: Dict[str, Any] = {
            "feature": f,
            "coverage": float(np.mean(np.isfinite(e[ev & np.isfinite(g)]))) if (ev & np.isfinite(g)).any() else np.nan,
            "rmse": float(np.sqrt(np.mean(err ** 2))) if len(err) else np.nan,
            "bias": float(np.mean(err)) if len(err) else np.nan,
            "nrmse": float(np.sqrt(np.mean(err ** 2)) / gsd) if len(err) and gsd > 0 else np.nan,
        }
        sok = stat & np.isfinite(g) & np.isfinite(e)
        row["stationary_sd"] = float(np.std(e[sok] - g[sok])) if sok.sum() > 2 else np.nan
        row["stationary_nsd"] = row["stationary_sd"] / gsd if gsd > 0 else np.nan
        lats = []
        for i, t_ev in enumerate(rec.events):
            pre = (t >= t_ev - 30.0) & (t < t_ev)
            nxt = boundaries[i + 1]
            post = (t >= nxt - 30.0) & (t < nxt)
            g_pre, g_post = np.nanmean(g[pre]) if pre.any() else np.nan, np.nanmean(g[post]) if post.any() else np.nan
            if not (np.isfinite(g_pre) and np.isfinite(g_post)) or abs(g_post - g_pre) < 0.5 * np.nanstd(g[ev]):
                continue
            level = 0.5 * (g_pre + g_post)
            rising = g_post > g_pre
            tg = _crossing_time(t, g, t_ev, nxt, level, rising)
            te = _crossing_time(t, e, t_ev, nxt, level, rising)
            lats.append(te - tg if np.isfinite(tg) and np.isfinite(te) else np.nan)
        row["latency_s"] = float(np.nanmean(lats)) if lats and np.isfinite(lats).any() else np.nan
        row["n_transitions"] = len(lats)
        rows.append(row)
    return pd.DataFrame(rows)


@dataclass
class BenchmarkResult:
    metrics: pd.DataFrame           # one row per (estimator, feature)
    summary: pd.DataFrame           # one row per estimator
    grids: Dict[str, pd.DataFrame]  # estimator -> feature grid
    ground_truth: pd.DataFrame      # ground truth on the grid
    result: PipelineResult


def evaluable_features(gt_grid: pd.DataFrame, min_coverage: float = 0.8, min_cv: float = 0.02) -> List[str]:
    """Features whose ground truth is defined on most of the grid and varies by at least ``min_cv``.

    A feature that the scenario barely changes (SD below ``min_cv`` times its
    mean magnitude) cannot measure tracking, and its normalized RMSE would be
    dominated by the estimation bias.
    """
    out = []
    for f in ALL_FEATURES:
        g = gt_grid[f].to_numpy(dtype=float)
        if np.mean(np.isfinite(g)) < min_coverage:
            continue
        sd, scale = np.nanstd(g), np.nanmean(np.abs(g))
        if sd > 1e-9 and sd >= min_cv * scale:
            out.append(f)
    return out


def feature_scales(rec: SyntheticRecording, config: Optional[PipelineConfig] = None,
                   gt: Optional[pd.DataFrame] = None) -> Dict[str, float]:
    """SD of each evaluable ground-truth feature of ``rec`` (typically a stress scenario).

    Used as fixed normalization scales so that scores are comparable across
    scenarios.
    """
    gt = ground_truth(rec, config) if gt is None else gt
    t = gt["t_peak"].to_numpy()
    grid = np.arange(np.ceil(t[0]), t[-1], 1.0)
    gt_grid = _gt_on_grid(gt, grid, ALL_FEATURES)
    return {f: float(np.nanstd(gt_grid[f])) for f in evaluable_features(gt_grid)}


def summarize(metrics: pd.DataFrame) -> pd.DataFrame:
    agg = metrics.groupby("estimator").agg(
        median_nrmse=("nrmse", "median"),
        median_stationary_nsd=("stationary_nsd", "median"),
        median_latency_s=("latency_s", "median"),
        mean_coverage=("coverage", "mean"),
    )
    return agg.sort_values("median_nrmse").reset_index()


def benchmark(
    rec: SyntheticRecording,
    config: Optional[PipelineConfig] = None,
    fixed_secs: Sequence[float] = DEFAULT_FIXED,
    ema_secs: Sequence[float] = DEFAULT_EMA,
    features: Optional[Sequence[str]] = None,
    eval_start: Optional[float] = None,
    gt: Optional[pd.DataFrame] = None,
    scales: Optional[Dict[str, float]] = None,
) -> BenchmarkResult:
    """Runs the pipeline on ``rec`` and scores every estimator against the ground truth."""
    cfg = (config or PipelineConfig()).copy()
    result = run(rec.signal_data, cfg, fixed_window_secs=fixed_secs)
    grid = result.features.grid["time"].to_numpy()
    gt = ground_truth(rec, cfg) if gt is None else gt
    gt_grid = _gt_on_grid(gt, grid, ALL_FEATURES)
    feats = list(features) if features is not None else (list(scales) if scales else evaluable_features(gt_grid))
    if eval_start is None:
        eval_start = max(cfg.engine.calibration_sec, max(fixed_secs, default=0.0))
    grids = estimator_grids(result, grid, ema_secs)
    frames = []
    for name, est in grids.items():
        m = feature_metrics(est, gt_grid, rec, feats, eval_start, scales)
        m.insert(0, "estimator", name)
        frames.append(m)
    metrics = pd.concat(frames, ignore_index=True)
    return BenchmarkResult(metrics, summarize(metrics), grids, gt_grid, result)


# ---------------------------------------------------------------------- #
# Sensitivity analysis
# ---------------------------------------------------------------------- #
def _set_path(cfg: PipelineConfig, path: str, value: Any) -> None:
    obj: Any = cfg
    parts = path.split(".")
    for p in parts[:-1]:
        obj = obj[p] if isinstance(obj, dict) else getattr(obj, p)
    last = parts[-1]
    if isinstance(obj, dict):
        obj[last] = value
    else:
        current = getattr(obj, last)
        if isinstance(current, tuple) and isinstance(value, list):
            value = tuple(value)
        setattr(obj, last, value)


def sweep(
    recordings: Sequence[SyntheticRecording],
    grid: Dict[str, Iterable[Any]],
    base: Optional[PipelineConfig] = None,
    fixed_secs: Sequence[float] = (),
    ema_secs: Sequence[float] = (),
    scales: Optional[Dict[str, float]] = None,
) -> pd.DataFrame:
    """Benchmarks the adaptive estimator for every combination in ``grid``.

    ``grid`` maps dotted config paths (e.g. ``"engine.theta"``,
    ``"engine.categories.macro.T_max_sec"``) to candidate values. Returns one
    row per combination and recording, with the adaptive summary metrics
    (normalized by ``scales`` when given, see :func:`feature_scales`).
    """
    base = (base or PipelineConfig()).copy()
    keys = list(grid)
    gts = [ground_truth(r, base) for r in recordings]
    rows = []
    for values in itertools.product(*(list(grid[k]) for k in keys)):
        cfg = copy.deepcopy(base)
        for k, v in zip(keys, values):
            _set_path(cfg, k, v)
        cfg.validate()
        for i, (rec, gt) in enumerate(zip(recordings, gts)):
            res = benchmark(rec, cfg, fixed_secs=fixed_secs, ema_secs=ema_secs, gt=gt, scales=scales,
                            eval_start=max(cfg.engine.calibration_sec, 60.0))
            s = res.summary.set_index("estimator").loc["adaptive"].to_dict()
            rows.append({**dict(zip(keys, values)), "recording": i, **s})
    return pd.DataFrame(rows)


__all__ = [
    "BenchmarkResult",
    "benchmark",
    "ema_series",
    "estimator_grids",
    "evaluable_features",
    "feature_metrics",
    "feature_scales",
    "ground_truth",
    "summarize",
    "sweep",
]
