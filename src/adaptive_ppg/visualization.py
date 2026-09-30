"""Plotly figures for pipeline and benchmark results (no Streamlit dependency)."""

from __future__ import annotations

from typing import List, Optional, Sequence, Tuple

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from .features import FEATURE_TO_CATEGORY, FEATURE_UNITS
from .pipeline import PipelineResult

CATEGORY_COLORS = {"macro": "#1f77b4", "time_volume": "#2ca02c", "derivatives": "#d62728"}
FIXED_COLORS = ("#7f7f7f", "#bcbd22", "#8c564b", "#17becf")
MAX_POINTS = 20_000


def _decimate(t: np.ndarray, y: np.ndarray, max_points: int = MAX_POINTS) -> Tuple[np.ndarray, np.ndarray]:
    """Min-max decimation: keeps the envelope of the signal with at most ~``max_points`` points."""
    n = len(y)
    if n <= max_points:
        return t, y
    bins = max_points // 2
    edges = np.linspace(0, n, bins + 1).astype(int)
    idx = []
    for a, b in zip(edges[:-1], edges[1:]):
        if b <= a:
            continue
        seg = y[a:b]
        i1, i2 = a + int(np.argmin(seg)), a + int(np.argmax(seg))
        idx.extend(sorted((i1, i2)))
    idx = np.asarray(idx)
    return t[idx], y[idx]


def _scatter(n_points: int):
    return go.Scattergl if n_points > 5_000 else go.Scatter


def invalid_spans(result: PipelineResult) -> List[Tuple[float, float]]:
    """Merged time spans covered by invalid beats."""
    ppg = result.ppg
    fs = ppg.sampling_rate
    spans: List[Tuple[float, float]] = []
    for k in np.flatnonzero(~ppg.valid_beats_mask):
        a, b = ppg.onsets[k] / fs, ppg.ends[k] / fs
        if spans and a <= spans[-1][1] + 1e-9:
            spans[-1] = (spans[-1][0], max(b, spans[-1][1]))
        else:
            spans.append((a, b))
    return spans


def _add_spans(fig: go.Figure, spans: Sequence[Tuple[float, float]], color: str, row: Optional[int] = None,
               label: str = "") -> None:
    kwargs = {"row": row, "col": 1} if row is not None else {}
    for i, (a, b) in enumerate(spans):
        fig.add_vrect(x0=a, x1=b, fillcolor=color, opacity=0.18, line_width=0, layer="below",
                      annotation_text=label if i == 0 and label else None, annotation_position="top left",
                      **kwargs)


def plot_signal(result: PipelineResult, show_raw: bool = True,
                artifacts: Optional[Sequence[Tuple[float, float]]] = None,
                max_points: int = MAX_POINTS) -> go.Figure:
    """Filtered (and optionally raw) signal with onsets, peaks and invalid-beat spans."""
    ppg = result.ppg
    fs = ppg.sampling_rate
    t = ppg.time
    pp = result.config.preprocessing
    title = (f"PPG — morphology band {pp.lowcut:g}–{ppg.effective_highcut:g} Hz "
             f"(detection {pp.detection_lowcut:g}–{pp.detection_highcut:g} Hz), fs = {fs:g} Hz")
    if pp.bypass_filter:
        title = f"PPG — unfiltered, fs = {fs:g} Hz"
    rows = 2 if show_raw else 1
    fig = make_subplots(rows=rows, cols=1, shared_xaxes=True, vertical_spacing=0.06,
                        row_heights=[0.65, 0.35][:rows] if show_raw else None,
                        subplot_titles=("Filtered (morphology branch)", "Raw")[:rows])
    Sc = _scatter(min(len(t), max_points))
    td, yd = _decimate(t, ppg.filtered_signal, max_points)
    fig.add_trace(Sc(x=td, y=yd, mode="lines", name="filtered", line=dict(width=1, color="#1f77b4")), row=1, col=1)
    valid = ppg.valid_beats_mask
    for mask, name, color in ((valid, "valid peaks", "#2ca02c"), (~valid, "rejected peaks", "#d62728")):
        fig.add_trace(go.Scatter(x=ppg.peaks[mask] / fs, y=ppg.filtered_signal[ppg.peaks[mask]], mode="markers",
                                 name=name, marker=dict(size=5, color=color)), row=1, col=1)
    fig.add_trace(go.Scatter(x=ppg.onsets / fs, y=ppg.filtered_signal[ppg.onsets], mode="markers", name="onsets",
                             marker=dict(size=4, color="#ff7f0e", symbol="triangle-up")), row=1, col=1)
    if show_raw:
        tr, yr = _decimate(t, ppg.raw_signal, max_points)
        fig.add_trace(Sc(x=tr, y=yr, mode="lines", name="raw", line=dict(width=1, color="#7f7f7f")), row=2, col=1)
    spans = invalid_spans(result)
    for r in range(1, rows + 1):
        _add_spans(fig, spans, "#d62728", row=r, label="rejected" if r == 1 else "")
        if artifacts:
            _add_spans(fig, artifacts, "#9467bd", row=r, label="artefact" if r == 1 else "")
    fig.update_layout(title=title, height=520 if show_raw else 380, hovermode="x unified",
                      margin=dict(l=40, r=20, t=70, b=40), legend=dict(orientation="h", y=-0.12))
    fig.update_xaxes(title_text="Time (s)", row=rows, col=1)
    return fig


def plot_variability(result: PipelineResult) -> go.Figure:
    """Variability, H, window length per category and overlap versus time."""
    a = result.adaptive
    e = result.config.engine
    t = a.beat_times
    fig = make_subplots(rows=4, cols=1, shared_xaxes=True, vertical_spacing=0.04,
                        subplot_titles=(f"Morphological variability ({e.variability_metric})",
                                        "Normalized variability H", "Window length T_w (valid beats)",
                                        "Overlap O_w"))
    fig.add_trace(go.Scatter(x=t, y=a.variability_used, mode="lines", name="variability",
                             line=dict(color="#9467bd")), row=1, col=1)
    if np.isfinite(a.norm_lo).any():
        fig.add_trace(go.Scatter(x=t, y=a.norm_lo, mode="lines", name="H = 0 level",
                                 line=dict(color="#7f7f7f", dash="dot")), row=1, col=1)
        fig.add_trace(go.Scatter(x=t, y=a.norm_hi, mode="lines", name="H = 1 level",
                                 line=dict(color="#7f7f7f", dash="dash")), row=1, col=1)
    fig.add_trace(go.Scatter(x=t, y=a.H, mode="lines", name="H", line=dict(color="#e377c2")), row=2, col=1)
    if e.mapping == "sigmoid":
        fig.add_hline(y=e.theta, line=dict(color="#333", dash="dash"), annotation_text=f"θ = {e.theta:g}",
                      row=2, col=1)
    for name, sched in a.schedules.items():
        color = CATEGORY_COLORS.get(name, None)
        fig.add_trace(go.Scatter(x=t, y=sched.T_beats, mode="lines", name=f"T_w {name}",
                                 line=dict(color=color, shape="hv")), row=3, col=1)
        wt = [w.t_anchor for w in sched.windows]
        wn = [w.T_beats for w in sched.windows]
        fig.add_trace(go.Scatter(x=wt, y=wn, mode="markers", name=f"anchors {name}", showlegend=False,
                                 marker=dict(color=color, size=4)), row=3, col=1)
    fig.add_trace(go.Scatter(x=t, y=a.overlap, mode="lines", name="O_w", line=dict(color="#8c564b")), row=4, col=1)
    warm = a.warmup_mask
    if warm.any():
        spans, start = [], None
        for i, w in enumerate(warm):
            if w and start is None:
                start = t[i]
            if (not w or i == len(warm) - 1) and start is not None:
                spans.append((start, t[i]))
                start = None
        for r in (1, 2):
            _add_spans(fig, spans, "#bcbd22", row=r, label="warm-up" if r == 1 else "")
    fig.update_layout(height=760, hovermode="x unified", margin=dict(l=40, r=20, t=50, b=40),
                      legend=dict(orientation="h", y=-0.08))
    fig.update_yaxes(range=[-0.05, 1.05], row=2, col=1)
    fig.update_xaxes(title_text="Time (s)", row=4, col=1)
    return fig


def plot_feature(result: PipelineResult, feature: str, fixed: Optional[Sequence[str]] = None,
                 show_beats: bool = True, ground_truth: Optional[Tuple[np.ndarray, np.ndarray]] = None) -> go.Figure:
    """One feature: adaptive windows (with CI95), fixed windows, per-beat values and optional ground truth."""
    fr = result.features
    cat = FEATURE_TO_CATEGORY[feature]
    unit = FEATURE_UNITS.get(feature, "")
    fig = go.Figure()
    if show_beats:
        pb = fr.per_beat[fr.per_beat["valid"]]
        Sc = _scatter(len(pb))
        fig.add_trace(Sc(x=pb["t_peak"], y=pb[feature], mode="markers", name="per beat",
                         marker=dict(size=3, color="#aaaaaa"), opacity=0.6))
    if ground_truth is not None:
        fig.add_trace(go.Scatter(x=ground_truth[0], y=ground_truth[1], mode="lines", name="ground truth",
                                 line=dict(color="black", width=2, dash="dot")))
    for i, name in enumerate(fixed if fixed is not None else list(fr.fixed)):
        tb = fr.fixed[name].table
        if tb.empty:
            continue
        fig.add_trace(go.Scatter(x=tb["t_anchor"], y=tb[feature], mode="lines+markers", name=name,
                                 line=dict(color=FIXED_COLORS[i % len(FIXED_COLORS)], width=1.5),
                                 marker=dict(size=4)))
    tb = fr.adaptive[cat].table if cat in fr.adaptive else None
    if tb is not None and not tb.empty:
        color = CATEGORY_COLORS.get(cat, "#1f77b4")
        y, ci = tb[feature].to_numpy(float), tb[f"{feature}_ci95"].to_numpy(float)
        ok = np.isfinite(y) & np.isfinite(ci)
        if ok.any():
            x = tb["t_anchor"].to_numpy()[ok]
            fig.add_trace(go.Scatter(x=np.concatenate([x, x[::-1]]),
                                     y=np.concatenate([(y + ci)[ok], (y - ci)[ok][::-1]]),
                                     fill="toself", fillcolor=color, opacity=0.15, line=dict(width=0),
                                     name="adaptive CI95", hoverinfo="skip"))
        fig.add_trace(go.Scatter(x=tb["t_anchor"], y=y, mode="lines+markers", name=f"adaptive ({cat})",
                                 line=dict(color=color, width=2.5, shape="hv" if result.adaptive.anchor == "trailing"
                                           else "linear"),
                                 marker=dict(size=5), customdata=np.stack([tb["n_valid"], tb["T_beats"]], axis=1),
                                 hovertemplate="%{y:.4g}<br>n_valid=%{customdata[0]}<br>T_w=%{customdata[1]}"))
    fig.update_layout(title=f"{feature} [{unit}] — category: {cat}", height=420, hovermode="x unified",
                      xaxis_title="Time (s)", yaxis_title=f"{feature} ({unit})",
                      margin=dict(l=40, r=20, t=60, b=40), legend=dict(orientation="h", y=-0.2))
    return fig


def plot_templates(result: PipelineResult, category: str, max_templates: int = 12) -> go.Figure:
    """Ensemble templates of a category, coloured by time."""
    ws = result.features.adaptive[category]
    fs = result.ppg.sampling_rate
    n = len(ws.templates)
    fig = go.Figure()
    if n == 0:
        fig.update_layout(title=f"No {category} windows")
        return fig
    pick = np.unique(np.linspace(0, n - 1, min(max_templates, n)).astype(int))
    for j, i in enumerate(pick):
        tp = ws.templates[i]
        row = ws.table.iloc[i]
        shade = j / max(len(pick) - 1, 1)
        color = f"rgba({int(30 + 200 * shade)}, {int(120 - 60 * shade)}, {int(220 - 180 * shade)}, 0.9)"
        fig.add_trace(go.Scatter(x=np.arange(len(tp.core)) / fs, y=tp.core, mode="lines",
                                 name=f"t={row['t_anchor']:.0f}s (n={int(row['n_valid'])})", line=dict(color=color)))
    fig.update_layout(title=f"Ensemble templates — {category}", xaxis_title="Time from foot (s)",
                      yaxis_title="Amplitude (a.u.)", height=420, margin=dict(l=40, r=20, t=60, b=40))
    return fig


def plot_benchmark(metrics, value: str = "nrmse") -> go.Figure:
    """Box plot of a benchmark metric across features, per estimator."""
    fig = go.Figure()
    order = metrics.groupby("estimator")[value].median().sort_values().index
    for name in order:
        m = metrics[metrics["estimator"] == name]
        fig.add_trace(go.Box(y=m[value], name=name, boxpoints="all", jitter=0.4, pointpos=0,
                             text=m["feature"], hovertemplate="%{text}: %{y:.3f}"))
    fig.update_layout(title=f"Benchmark — {value} across features", yaxis_title=value, height=420,
                      showlegend=False, margin=dict(l=40, r=20, t=60, b=40))
    return fig


__all__ = [
    "invalid_spans",
    "plot_benchmark",
    "plot_feature",
    "plot_signal",
    "plot_templates",
    "plot_variability",
]
