"""Streamlit dashboard for the adaptive-window PPG pipeline.

Run from the repository root with ``streamlit run app/streamlit_app.py``
(or ``adaptive-ppg dashboard``; both require the ``app`` extra:
``pip install -e ".[app]"``). All computation is done by the
:mod:`adaptive_ppg` package; this file only builds the user interface.
"""

from __future__ import annotations

import io
import json
import pathlib
import sys

import numpy as np
import streamlit as st

try:
    import adaptive_ppg  # noqa: F401
except ImportError:  # running from a source checkout without installing the package
    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from adaptive_ppg import __version__
from adaptive_ppg.benchmark import benchmark, feature_scales
from adaptive_ppg.config import FEATURE_CATEGORIES_ORDER, PRESET_NAMES, PipelineConfig, preset
from adaptive_ppg.features import FEATURE_CATEGORIES, FEATURE_UNITS
from adaptive_ppg.io import BioSignalLoader, BioSignalLoaderError
from adaptive_ppg.pipeline import run
from adaptive_ppg.synthetic import generate, rest_scenario, stress_scenario
from adaptive_ppg.visualization import plot_benchmark, plot_feature, plot_signal, plot_templates, plot_variability

ROOT = pathlib.Path(__file__).resolve().parents[1]
SAMPLE_DIR = ROOT / "sample_data"
# Known layout of the bundled recordings: (PPG column, ECG column, sampling rate or None to infer).
SAMPLE_DEFAULTS = {
    "bidmc_01_Signals.csv": ("PLETH", "II", None),
    "Sample1.CSV": ("", None, None),  # PhysioNet PTT PPG dataset: 500 Hz time column, inverted MAX30101 counts
    "Sample2.CSV": ("", None, None),
}
SAMPLE_NOTES = {
    "Sample1.CSV": "Pulse Transit Time PPG Dataset (Mehrgardt et al., PhysioNet 2022, ODbL 1.0). "
                   "IR channel, distal phalanx, seated at rest. "
                   "Raw MAX30101 counts are inverted: the automatic polarity check flips them.",
}
SAMPLE_NOTES["Sample2.CSV"] = SAMPLE_NOTES["Sample1.CSV"]
POLARITIES = ("auto", "normal", "inverted")
METRICS = ("drift_z", "drift", "successive", "dispersion")
NORMALIZATIONS = ("fixed", "calibration", "rolling", "global")


# ---------------------------------------------------------------------- #
# Cached computation
# ---------------------------------------------------------------------- #
@st.cache_resource(max_entries=4, show_spinner="Loading signal...")
def load_file(data: bytes, name: str, signal_column: str, ecg_column: str, fs: float | None):
    col = signal_column.strip() or None  # blank: first non-time column (CSV) or channel 0
    if col is not None and col.isdigit():
        col = int(col)
    ecg = None
    if ecg_column.strip():
        ecg = int(ecg_column) if ecg_column.strip().isdigit() else ecg_column.strip()
    buf = io.BytesIO(data)
    buf.name = name
    return BioSignalLoader.auto_load(buf, file_name=name, signal_column=col, sampling_rate=fs, ecg_column=ecg)


@st.cache_resource(max_entries=4, show_spinner="Generating synthetic recording...")
def synthetic_recording(scenario: str, fs: float, seed: int):
    segments = stress_scenario() if scenario == "stress" else rest_scenario()
    return generate(segments, fs=fs, seed=seed)


@st.cache_resource(max_entries=6, show_spinner="Running the pipeline...")
def analyse(source_key: str, _signal_data, config_json: str, fixed: tuple, detector: str):
    return run(_signal_data, PipelineConfig.from_dict(json.loads(config_json)), fixed_window_secs=list(fixed),
               detector=detector)


@st.cache_resource(max_entries=4, show_spinner="Scoring against the ground truth...")
def analyse_synthetic(scenario: str, fs: float, seed: int, config_json: str, fixed: tuple):
    cfg = PipelineConfig.from_dict(json.loads(config_json))
    rec = synthetic_recording(scenario, fs, seed)
    scales = feature_scales(synthetic_recording("stress", fs, 0), cfg)
    return rec, benchmark(rec, cfg, fixed_secs=fixed, ema_secs=(10.0,), scales=scales)


# ---------------------------------------------------------------------- #
# Sidebar
# ---------------------------------------------------------------------- #
def sidebar_source():
    """Returns ``(kind, payload)``: ('file', (bytes, name, col, ecg, fs)) or ('synthetic', (scenario, fs, seed))."""
    with st.sidebar.expander("1. Data", expanded=True):
        options = ["Synthetic demo (known ground truth)", "Sample data", "Upload a file"]
        choice = st.radio("Source", options, index=1 if SAMPLE_DIR.exists() else 0)
        if choice.startswith("Synthetic"):
            scenario = st.selectbox("Scenario", ("stress", "rest"),
                                    help="stress: rest → fast acute stress → slow recovery; rest: stationary")
            fs = st.number_input("Sampling rate (Hz)", 25.0, 1000.0, 125.0, 25.0)
            seed = int(st.number_input("Noise seed", 0, 10_000, 1))
            return "synthetic", (scenario, float(fs), seed)

        if choice == "Sample data":
            files = sorted(p.name for p in SAMPLE_DIR.glob("*") if p.suffix.lower() in (".csv", ".txt", ".tsv",
                                                                                         ".mat", ".edf", ".bdf"))
            if not files:
                st.warning("No files in sample_data/.")
                return None, None
            name = st.selectbox("File", files)
            data = (SAMPLE_DIR / name).read_bytes()
        else:
            up = st.file_uploader("PPG recording", type=["csv", "txt", "tsv", "mat", "edf", "bdf"])
            if up is None:
                return None, None
            name, data = up.name, up.getvalue()

        d_col, d_ecg, d_fs = SAMPLE_DEFAULTS.get(name, ("", None, None))
        if choice == "Sample data" and name in SAMPLE_NOTES:
            st.caption(SAMPLE_NOTES[name])
        signal_column = st.text_input("PPG column / channel", d_col,
                                      help="Name or 0-based index; blank = first non-time column (or channel 0)")
        ecg_column = st.text_input("ECG column / channel (optional, enables PAT)", d_ecg or "")
        infer = st.checkbox("Infer the sampling rate from the file", value=d_fs is None,
                            help="Uses the time column (CSV) or the header (EDF/BDF, always used).")
        fs = None
        if not infer:
            fs = float(st.number_input("Sampling rate (Hz)", 1.0, 10_000.0, d_fs or 125.0, 1.0))
        return "file", (data, name, signal_column, ecg_column, fs)


def sidebar_config() -> tuple[PipelineConfig, tuple, str]:
    with st.sidebar.expander("2. Method", expanded=True):
        name = st.selectbox("Preset", PRESET_NAMES, help="Starting values; the controls below override them.")
        cfg = preset(name)
        e, p, f = cfg.engine, cfg.preprocessing, cfg.features
        e.variability_metric = st.selectbox("Variability metric", METRICS, index=METRICS.index(e.variability_metric),
                                            help="drift_z ≈ 1 on a stationary signal, whatever the noise level")
        e.normalization = st.selectbox("H normalization", NORMALIZATIONS, index=NORMALIZATIONS.index(e.normalization))
        if e.normalization == "fixed":
            lo_default, hi_default = e.fixed_scale
            c1, c2 = st.columns(2)
            lo = c1.number_input("H = 0 at", value=float(lo_default), format="%.3f")
            hi = c2.number_input("H = 1 at", value=float(hi_default), format="%.3f")
            e.fixed_scale = (lo, max(hi, lo + 1e-6))
        e.theta = st.slider("θ (H at which windows start shrinking)", 0.05, 0.95, float(e.theta), 0.05)
        e.sigmoid_slope = st.slider("Sigmoid slope k", 2.0, 30.0, float(e.sigmoid_slope), 1.0)
        e.N_past = st.slider("N_past (beats per drift block)", 4, 30, int(e.N_past))
        e.anchor = st.radio("Window anchor", ("trailing", "centered"), horizontal=True,
                            index=0 if e.anchor == "trailing" else 1, help="trailing = causal (real-time)")
        e.ibi_correction = st.checkbox("Regress out IBI-driven variability", value=e.ibi_correction)

    with st.sidebar.expander("3. Windows per category", expanded=False):
        for cat in FEATURE_CATEGORIES_ORDER:
            cc = e.categories[cat]
            st.markdown(f"**{cat}**")
            c1, c2, c3 = st.columns(3)
            cc.T_min_beats = int(c1.number_input("T_min (beats)", 1, 60, cc.T_min_beats, key=f"tmin_{cat}"))
            cc.T_max_sec = float(c2.number_input("T_max (s)", 5.0, 300.0, cc.T_max_sec, 5.0, key=f"tmax_{cat}"))
            cc.T_crit_beats = int(c3.number_input("T_crit", 1, 60, cc.T_crit_beats, key=f"tcrit_{cat}"))

    with st.sidebar.expander("4. Signal processing", expanded=False):
        resample = st.checkbox("Resample to a common rate", value=p.target_fs is not None,
                               help="Makes the analysis independent of the acquisition rate and faster")
        if resample:
            p.target_fs = float(st.number_input("Processing rate (Hz)", 50.0, 1000.0,
                                                float(p.target_fs or 125.0), 25.0))
        else:
            p.target_fs = None
        p.polarity = st.radio("Signal polarity", POLARITIES, horizontal=True, index=POLARITIES.index(p.polarity),
                              help="auto: flips signals whose derivative is negatively skewed "
                                   "(e.g. raw absorbance counts)")
        c1, c2 = st.columns(2)
        p.lowcut = c1.number_input("Morphology high-pass (Hz)", 0.05, 5.0, float(p.lowcut), 0.05)
        p.highcut = c2.number_input("Morphology low-pass (Hz)", 5.0, 60.0, float(p.highcut), 1.0)
        p.bypass_filter = st.checkbox("Bypass the morphology filter", value=p.bypass_filter)
        p.onset_method = st.radio("Onset method", ("tangent", "minimum"), horizontal=True,
                                  index=0 if p.onset_method == "tangent" else 1)
        p.min_template_corr = st.slider("Min. correlation with the reference beat", 0.5, 0.99,
                                        float(p.min_template_corr), 0.01)
        detector = st.radio("Peak detector", ("elgendi", "biosppy"), horizontal=True,
                            help="biosppy needs the 'biosppy' extra (and peakutils)")
        height = st.number_input("Subject height (m, 0 = unknown)", 0.0, 2.5, 0.0, 0.01,
                                 help="Enables the stiffness index")
        f.subject_height_m = height or None

    with st.sidebar.expander("5. Fixed-window baselines", expanded=False):
        text = st.text_input("Fixed windows (s, comma-separated)", "10, 30")
        try:
            fixed = tuple(sorted({float(x) for x in text.replace(";", ",").split(",") if x.strip()}))
        except ValueError:
            st.error("Use numbers separated by commas.")
            fixed = (30.0,)
        f.fixed_overlap_frac = st.slider("Overlap of fixed windows", 0.0, 0.9, float(f.fixed_overlap_frac), 0.1)

    try:
        cfg.validate()
    except ValueError as exc:
        st.sidebar.error(f"Invalid configuration: {exc}")
        st.stop()
    return cfg, fixed, detector


# ---------------------------------------------------------------------- #
# Main panels
# ---------------------------------------------------------------------- #
def kpis(res) -> None:
    ppg, a, fr = res.ppg, res.adaptive, res.features
    n, n_valid = ppg.n_beats, int(ppg.valid_beats_mask.sum())
    ok = np.isfinite(a.H) & ~a.warmup_mask
    transient = float(np.mean(a.H[ok] > res.config.engine.theta)) if ok.any() else np.nan
    macro = a.schedules.get("macro")
    med_T = float(np.nanmedian(macro.T_beats)) if macro is not None else np.nan
    med_ibi = float(np.nanmedian(a.local_ibi)) if np.isfinite(a.local_ibi).any() else np.nan
    sqi = float(fr.adaptive["macro"].table["template_sqi"].median()) if "macro" in fr.adaptive else np.nan
    c = st.columns(4)
    c[0].metric("Valid beats", f"{n_valid} / {n}", f"{100 * n_valid / max(n, 1):.0f} %", delta_color="off")
    c[1].metric("Beats with H > θ", f"{100 * transient:.1f} %" if np.isfinite(transient) else "–",
                help="Share of beats (after warm-up) where the windows contract")
    c[2].metric("Median macro window", f"{med_T:.0f} beats",
                f"≈ {med_T * med_ibi:.0f} s" if np.isfinite(med_T * med_ibi) else None, delta_color="off")
    c[3].metric("Median template SQI", f"{sqi:.3f}" if np.isfinite(sqi) else "–",
                help="Mean correlation of the beats with their ensemble template (macro windows)")


def tab_features(res, gt_grid=None) -> None:
    all_feats = [f for cat in FEATURE_CATEGORIES_ORDER for f in FEATURE_CATEGORIES[cat]]
    available = [f for f in all_feats if f in res.features.per_beat and res.features.per_beat[f].notna().any()]
    default = [f for f in ("heart_rate", "peak_amplitude", "crest_time", "sdptg_b_a") if f in available]
    chosen = st.multiselect("Features", available, default=default,
                            format_func=lambda f: f"{f} [{FEATURE_UNITS.get(f, '')}]")
    c1, c2 = st.columns(2)
    show_beats = c1.checkbox("Show per-beat values", value=True)
    fixed_names = list(res.features.fixed)
    fixed_sel = c2.multiselect("Fixed-window baselines", fixed_names, default=fixed_names)
    for feat in chosen:
        gt = None
        if gt_grid is not None and feat in gt_grid:
            gt = (gt_grid["time"].to_numpy(), gt_grid[feat].to_numpy())
        st.plotly_chart(plot_feature(res, feat, fixed=fixed_sel, show_beats=show_beats, ground_truth=gt),
                        width="stretch")
    with st.expander("Adaptive windows table"):
        st.dataframe(res.window_table(), width="stretch")
    with st.expander("Per-beat features"):
        st.dataframe(res.features.per_beat, width="stretch")


def tab_templates(res) -> None:
    cats = [c for c in FEATURE_CATEGORIES_ORDER if c in res.features.adaptive]
    c1, c2 = st.columns([2, 1])
    cat = c1.selectbox("Category", cats)
    n = c2.slider("Templates shown", 2, 30, 12)
    st.plotly_chart(plot_templates(res, cat, max_templates=n), width="stretch")


def tab_benchmark(bench) -> None:
    if bench is None:
        st.info("The benchmark needs a known ground truth: choose **Synthetic demo** as the data source. "
                "For larger studies use `adaptive-ppg benchmark` and `adaptive-ppg sweep`.")
        return
    st.markdown("Scores of each estimator against the noise-free ground truth, normalized by the "
                "ground-truth SD in the stress scenario (lower is better).")
    st.dataframe(bench.summary.round(3), width="stretch", hide_index=True)
    metric = st.radio("Metric", ("nrmse", "stationary_nsd", "latency_s", "bias"), horizontal=True)
    st.plotly_chart(plot_benchmark(bench.metrics, metric), width="stretch")
    with st.expander("Per-feature metrics"):
        st.dataframe(bench.metrics, width="stretch")


def export_ui(res, source_name: str) -> None:
    stem = pathlib.Path(source_name).stem or "recording"
    c1, c2, c3 = st.columns(3)
    templates = c1.checkbox("Include ensemble templates", value=False)
    c1.download_button("Download ZIP (CSV + metadata.json)", res.to_zip(include_templates=templates),
                       file_name=f"{stem}_adaptive_ppg.zip", mime="application/zip")
    try:
        xlsx = res.to_excel()
        c2.download_button("Download Excel workbook", xlsx, file_name=f"{stem}_adaptive_ppg.xlsx",
                           mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
    except ImportError:
        c2.caption("Excel export needs openpyxl.")
    c3.download_button("Download configuration (JSON)", res.config.to_json(), file_name="config.json",
                       mime="application/json")
    with st.expander("Run metadata"):
        meta = res.metadata()
        meta.pop("config", None)
        st.json(meta)


def main() -> None:
    st.set_page_config(page_title="Adaptive PPG", page_icon="🫀", layout="wide")
    st.title("Adaptive-window PPG morphology")
    st.caption(f"adaptive_ppg {__version__} — variable-width ensemble averaging driven by morphological variability")

    kind, payload = sidebar_source()
    cfg, fixed, detector = sidebar_config()
    if kind is None:
        st.info("Choose a data source in the sidebar.")
        return
    config_json = cfg.to_json()

    rec = bench = None
    try:
        if kind == "synthetic":
            scenario, fs, seed = payload
            rec, bench = analyse_synthetic(scenario, fs, seed, config_json, fixed)
            res, source_name = bench.result, f"synthetic_{scenario}_seed{seed}"
        else:
            data, name, col, ecg, fs = payload
            sd = load_file(data, name, col, ecg, fs)
            key = f"{name}|{len(data)}|{hash(data)}|{col}|{ecg}|{fs}"
            res, source_name = analyse(key, sd, config_json, fixed, detector), name
    except (BioSignalLoaderError, ValueError, KeyError, ImportError) as exc:
        st.error(f"Could not analyse the recording: {exc}")
        return

    sd = res.signal_data
    ppg = res.ppg
    rate = f"{sd.sampling_rate:g} Hz"
    if ppg.sampling_rate != sd.sampling_rate:
        rate += f" (processed at {ppg.sampling_rate:g} Hz)"
    st.markdown(f"**{source_name}** — channel `{sd.channel_name}`, {rate}, {sd.duration:.0f} s"
                + (", with ECG (PAT enabled)" if "ecg" in sd.aux_signals else "")
                + (", inverted" if ppg.polarity < 0 else ""))
    for w in res.ppg.warnings:
        st.warning(w)
    if sd.metadata.get("sampling_rate_warning"):
        st.warning(sd.metadata["sampling_rate_warning"])
    if res.ppg.valid_beats_mask.sum() < 20:
        st.error("Fewer than 20 valid beats: check the column, the sampling rate and the signal polarity.")

    kpis(res)
    tabs = st.tabs(["Signal & beats", "Variability & windows", "Features", "Templates", "Benchmark", "Export"])
    with tabs[0]:
        show_raw = st.checkbox("Show the raw signal", value=True)
        st.plotly_chart(plot_signal(res, show_raw=show_raw, artifacts=rec.artifacts if rec else None),
                        width="stretch")
        with st.expander("Beat quality table"):
            st.dataframe(res.beat_table(), width="stretch")
    with tabs[1]:
        st.plotly_chart(plot_variability(res), width="stretch")
    with tabs[2]:
        tab_features(res, bench.ground_truth if bench else None)
    with tabs[3]:
        tab_templates(res)
    with tabs[4]:
        tab_benchmark(bench)
    with tabs[5]:
        export_ui(res, source_name)


if __name__ == "__main__":
    main()
