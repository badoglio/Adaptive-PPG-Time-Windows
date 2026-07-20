"""Streamlit orchestrator dashboard for Adaptive PPG Analysis.

This file integrates data loading, preprocessing, adaptive window calculation,
and morphological feature extraction in a responsive web app.
"""

import io
import os
import pathlib
import warnings
from typing import Dict, List, Tuple, Union, Any, Optional

import numpy as np
import pandas as pd
import streamlit as st

# Safe import of local pipeline modules with Mock fallbacks
try:
    from io_loader import BioSignalLoader, SignalData, BioSignalLoaderError
except ImportError:
    st.warning("Failed to import io_loader. Using mock implementation.")
    from dataclasses import dataclass, field

    @dataclass
    class SignalData:
        signal: np.ndarray
        sampling_rate: float
        time: np.ndarray = field(default=None)  # type: ignore
        channel_name: str = "Channel_Mock"
        file_format: str = "MOCK"
        metadata: Optional[dict] = field(default_factory=dict)
        def __post_init__(self):
            if self.time is None:
                self.time = np.arange(len(self.signal)) / self.sampling_rate

    class BioSignalLoader:
        @staticmethod
        def auto_load(*args, **kwargs):
            return None

    class BioSignalLoaderError(Exception):
        pass

try:
    from preprocessing import PPGPreprocessor, ProcessedPPG, PPGPreprocessingError
except ImportError:
    st.warning("Failed to import preprocessing. Using mock implementation.")
    @dataclass
    class ProcessedPPG:
        raw_signal: np.ndarray
        filtered_signal: np.ndarray
        sampling_rate: float
        peaks: np.ndarray
        onsets: np.ndarray
        beats_raw: List[np.ndarray]
        beats_normalized: np.ndarray
        beat_times: np.ndarray
        valid_beats_mask: np.ndarray

    class PPGPreprocessor:
        def __init__(self, **kwargs): pass
        def process(self, *args, **kwargs): return None

    class PPGPreprocessingError(Exception):
        pass

try:
    from adaptive_engine import AdaptiveEngine, AdaptiveConfig, AdaptiveWindowResult, AdaptiveEngineError
except ImportError:
    st.warning("Failed to import adaptive_engine. Using mock implementation.")
    @dataclass
    class AdaptiveWindowResult:
        beat_times: np.ndarray
        entropy_raw: np.ndarray
        entropy_normalized: np.ndarray
        target_window_beats: np.ndarray
        overlap_ratio: np.ndarray
        feature_windows_beats: Dict[str, np.ndarray]
        window_slices: Dict[str, List[Tuple[int, int]]]

    class AdaptiveEngine:
        def __init__(self, *args): pass
        def process(self, *args): return None
        def apply_critical_sensitivity(self, *args): return {}
        def generate_window_slices(self, *args, **kwargs): return {}

    class AdaptiveConfig:
        def __init__(self, **kwargs): pass

    class AdaptiveEngineError(Exception):
        pass

try:
    from feature_extraction import PPGFeatureExtractor, FeatureExtractionResult, FeatureExtractionError
except ImportError:
    st.warning("Failed to import feature_extraction. Using mock implementation.")
    @dataclass
    class FeatureExtractionResult:
        single_beat_df: pd.DataFrame
        adaptive_features_df: pd.DataFrame
        fixed_features_df: pd.DataFrame
        comparison_metrics: dict

    class PPGFeatureExtractor:
        def process(self, *args, **kwargs): return None

    class FeatureExtractionError(Exception):
        pass

try:
    from visualization import PPGVisualizer, render_export_ui
except ImportError:
    st.warning("Failed to import visualization. Using mock implementation.")
    class PPGVisualizer:
        @staticmethod
        def plot_signal_and_beats(*args, **kwargs): return None
        @staticmethod
        def plot_entropy_and_windows(*args, **kwargs): return None
        @staticmethod
        def plot_feature_comparison(*args, **kwargs): return None
        @staticmethod
        def plot_multi_feature_dashboard(*args, **kwargs): return None
    def render_export_ui(*args, **kwargs): pass


# ==========================================
# CASCADE CACHED FUNCTIONS
# ==========================================

@st.cache_data
def cached_load_data(
    uploaded_file_bytes: Optional[bytes],
    file_name: str,
    file_format: str,
    channel_index: Union[int, str],
    sampling_rate_override: float
) -> SignalData:
    """Loads and formats input biosegnal data. Re-evaluates only on file/format changes.

    Args:
        uploaded_file_bytes (Optional[bytes]): File byte stream, or None for synthetic.
        file_name (str): Original file name.
        file_format (str): CSV, MAT, EDF, or SYNTHETIC.
        channel_index (Union[int, str]): Channel offset index or label.
        sampling_rate_override (float): Frequency in Hz used if metadata doesn't provide one.

    Returns:
        SignalData: Standardized biosegnal struct.
    """
    if uploaded_file_bytes is None:
        # Generate clean synthetic demonstration PPG (20 seconds, 100 Hz)
        fs = float(sampling_rate_override)
        t = np.arange(int(20.0 * fs)) / fs
        hr = 75.0
        period = 60.0 / hr
        phase = (t % period) / period
        systolic_wave = np.exp(-((phase - 0.15) / 0.07)**2)
        diastolic_wave = 0.35 * np.exp(-((phase - 0.45) / 0.10)**2)
        clean_ppg = systolic_wave + diastolic_wave
        baseline_wander = 0.25 * np.sin(2 * np.pi * 0.15 * t)
        noise = 0.04 * np.random.randn(len(t))
        sig = clean_ppg + baseline_wander + noise
        return SignalData(
            signal=sig,
            sampling_rate=fs,
            time=t,
            channel_name="PPG_Synthetic",
            file_format="CSV"
        )

    # Convert bytes into file-like object
    buf = io.BytesIO(uploaded_file_bytes)
    ext = pathlib.Path(file_name).suffix.lower()

    # Route parser depending on resolved format
    if ext in (".csv", ".txt", ".tsv"):
        # If channel_index is digits, cast to int
        resolved_col: Union[int, str] = channel_index
        try:
            resolved_col = int(channel_index)
        except ValueError:
            pass
        return BioSignalLoader.load_csv(
            buf,
            signal_column=resolved_col,
            sampling_rate=sampling_rate_override
        )
    elif ext == ".mat":
        # Resolve index
        resolved_idx = 0
        try:
            resolved_idx = int(channel_index)
        except ValueError:
            pass
        return BioSignalLoader.load_mat(
            buf,
            key=None,
            sampling_rate=sampling_rate_override,
            channel_index=resolved_idx
        )
    elif ext in (".edf", ".bdf"):
        resolved_idx = 0
        try:
            resolved_idx = int(channel_index)
        except ValueError:
            pass
        return BioSignalLoader.load_edf(buf, channel_index=resolved_idx)
    else:
        raise ValueError(f"Format extension '{ext}' not supported.")


@st.cache_data
def cached_preprocessing(
    signal_data: SignalData,
    lowcut: float,
    highcut: float,
    target_length: int
) -> ProcessedPPG:
    """Filters, segments, and normalizes signals. Re-evaluates only if cutoffs change.

    Args:
        signal_data (SignalData): Raw signals struct.
        lowcut (float): High-pass frequency cutoff in Hz.
        highcut (float): Low-pass frequency cutoff in Hz.
        target_length (int): Fixed beat sample width (typically 128).

    Returns:
        ProcessedPPG: Normalized segmented dataset.
    """
    preprocessor = PPGPreprocessor(
        lowcut=lowcut,
        highcut=highcut,
        target_beat_length=target_length
    )
    return preprocessor.process(signal_data)


@st.cache_data
def cached_adaptive_engine(
    processed_ppg: ProcessedPPG,
    N_past: int,
    T_min: int,
    T_max: int,
    O_min: float,
    O_max: float,
    theta: float
) -> AdaptiveWindowResult:
    """Computes morphological complexity curves. Re-evaluates on window/entropy configs changes.

    Args:
        processed_ppg (ProcessedPPG): Segmented beat arrays.
        N_past (int): Rolling beat count for entropy.
        T_min (int): Min beats window size.
        T_max (int): Max beats window size.
        O_min (float): Min overlap ratio.
        O_max (float): Max overlap ratio.
        theta (float): Slicing sensitivity threshold.

    Returns:
        AdaptiveWindowResult: Entropy indices and scheduled target sizing curves.
    """
    config = AdaptiveConfig(
        N_past=N_past,
        T_min_beats=T_min,
        T_max_beats=T_max,
        O_min=O_min,
        O_max=O_max,
        theta_sensitivity=theta
    )
    engine = AdaptiveEngine(config)
    return engine.process(processed_ppg)


@st.cache_data
def cached_feature_extraction(
    processed_ppg: ProcessedPPG,
    adaptive_result: AdaptiveWindowResult,
    T_crit_macro: int,
    T_crit_derivatives: int,
    fixed_window_sec: float
) -> FeatureExtractionResult:
    """Extracts features and aggregates them. Re-evaluates on T_crit bounds changes.

    Does NOT recompute morphological entropy.

    Args:
        processed_ppg (ProcessedPPG): Segmented beat arrays.
        adaptive_result (AdaptiveWindowResult): Sizing curves from engine.
        T_crit_macro (int): Minimum beats threshold for macro features.
        T_crit_derivatives (int): Minimum beats threshold for derivative features.
        fixed_window_sec (float): Benchmark window width in seconds (default 30s).

    Returns:
        FeatureExtractionResult: Combined beat-by-beat and comparative dataframes.
    """
    # Recalculate sensitivity bounds and schedules based on new constraints
    config = AdaptiveConfig(
        T_crit_dict={'macro': T_crit_macro, 'time_volume': T_crit_macro, 'derivatives': T_crit_derivatives}
    )
    engine = AdaptiveEngine(config)
    
    # Recalculate critical window limits and schedules
    res_windows = engine.apply_critical_sensitivity(adaptive_result.target_window_beats)
    res_slices = engine.generate_window_slices(res_windows, adaptive_result.overlap_ratio)

    updated_result = AdaptiveWindowResult(
        beat_times=adaptive_result.beat_times,
        entropy_raw=adaptive_result.entropy_raw,
        entropy_normalized=adaptive_result.entropy_normalized,
        target_window_beats=adaptive_result.target_window_beats,
        overlap_ratio=adaptive_result.overlap_ratio,
        feature_windows_beats=res_windows,
        window_slices=res_slices
    )

    extractor = PPGFeatureExtractor()
    return extractor.process(processed_ppg, updated_result, fixed_window_sec=fixed_window_sec)


# ==========================================
# STREAMLIT UI LAYOUT AND ORCHESTRATION
# ==========================================

def main():
    st.set_page_config(
        page_title="Adaptive PPG Analyzer",
        layout="wide",
        initial_sidebar_state="expanded"
    )

    # CSS Injection for Premium Visual Theme
    st.markdown(
        """
        <style>
        .reportview-container {
            background: #0f1116;
        }
        .main-title {
            color: #4f8bf9;
            font-size: 2.5rem;
            font-weight: 700;
            margin-bottom: 0.1rem;
        }
        .subtitle {
            color: #8a99ad;
            font-size: 1.1rem;
            margin-bottom: 1.5rem;
        }
        .badge-status {
            background-color: #1e293b;
            color: #38bdf8;
            padding: 0.4rem 0.8rem;
            border-radius: 0.5rem;
            font-weight: 600;
            display: inline-block;
            margin-bottom: 1rem;
            border: 1px solid #334155;
        }
        div[data-testid="metric-container"] {
            background-color: #1e293b;
            border: 1px solid #334155;
            padding: 0.8rem 1.2rem;
            border-radius: 0.8rem;
            box-shadow: 0 4px 6px -1px rgb(0 0 0 / 0.1);
        }
        </style>
        """,
        unsafe_allow_html=True
    )

    st.markdown("<div class='main-title'>Adaptive PPG Windowing Dashboard</div>", unsafe_allow_html=True)
    st.markdown("<div class='subtitle'>Variable-width ensemble averaging based on morphological state complexity</div>", unsafe_allow_html=True)

    # ----------------- SIDEBAR CONTROLS -----------------
    st.sidebar.markdown("### Parameters & Controls")

    # Sidebar 1: Data Loading
    with st.sidebar.expander("📂 1. Data Loading", expanded=True):
        use_synth = st.checkbox("Use Demo Synthetic Data", value=True)
        uploaded_file = None
        file_bytes = None
        file_name = "synthetic"
        
        if not use_synth:
            uploaded_file = st.file_uploader(
                "Upload a biosegnal file",
                type=["csv", "txt", "mat", "edf"]
            )
            if uploaded_file is not None:
                file_bytes = uploaded_file.getvalue()
                file_name = uploaded_file.name
        
        st.markdown("---")
        channel_index = st.text_input("Channel or Column (Name/Index)", value="0")
        sampling_rate_override = st.number_input("Sampling Frequency (Hz)", min_value=1.0, value=100.0, step=1.0)

    # Sidebar 2: Bandpass Filter
    with st.sidebar.expander("𝄠 2. Bandpass Filter", expanded=False):
        lowcut = st.slider("High-Pass Cutoff (Hz)", min_value=0.1, max_value=5.0, value=0.5, step=0.1)
        highcut = st.slider("Low-Pass Cutoff (Hz)", min_value=5.0, max_value=40.0, value=30.0, step=1.0)
        target_length = st.number_input("Beat Resampling Points", min_value=16, max_value=512, value=128, step=16)

    # Sidebar 3: Entropy & Windowing Parameters
    with st.sidebar.expander("🌀 3. Entropy & Windowing", expanded=False):
        N_past = st.slider("N_past (Historical Beats)", min_value=5, max_value=30, value=10, step=1)
        t_limits = st.slider("Tw Window Limits (Beats)", min_value=3, max_value=60, value=(3, 30), step=1)
        T_min, T_max = t_limits
        o_limits = st.slider("Ow Overlap Ratio Limits", min_value=0.0, max_value=1.0, value=(0.25, 0.85), step=0.05)
        O_min, O_max = o_limits
        theta = st.slider("Entropy Sensitivity θ", min_value=0.05, max_value=1.0, value=0.50, step=0.05)

    # Sidebar 4: Critical Limits & Benchmark
    with st.sidebar.expander("🎯 4. Critical Limits & Benchmark", expanded=False):
        T_crit_macro = st.number_input("T_crit Macro-Features (Beats)", min_value=1, max_value=15, value=3)
        T_crit_derivatives = st.number_input("T_crit Derivative SDPTG (Beats)", min_value=1, max_value=30, value=10)
        fixed_window_sec = st.number_input("Fixed Benchmark Window (Seconds)", min_value=5.0, max_value=120.0, value=30.0, step=5.0)

    # Warn if no signal is active
    if not use_synth and uploaded_file is None:
        st.info("💡 Upload a biosegnal file (.csv, .mat, .edf) from the sidebar or enable synthetic data to start.")
        return

    # ----------------- CASCADE PIPELINE EXECUTION -----------------
    try:
        # Step 1: Load Data
        with st.spinner("Loading biosegnal data..."):
            signal_data = cached_load_data(
                uploaded_file_bytes=file_bytes,
                file_name=file_name,
                file_format=pathlib.Path(file_name).suffix.upper().replace(".", ""),
                channel_index=channel_index,
                sampling_rate_override=sampling_rate_override
            )

        # Step 2: Preprocessing
        with st.spinner("Preprocessing and segmenting beats..."):
            processed_ppg = cached_preprocessing(
                signal_data=signal_data,
                lowcut=lowcut,
                highcut=highcut,
                target_length=target_length
            )

        # Step 3: Adaptive Windowing
        with st.spinner("Computing entropy and adaptive window slices..."):
            adaptive_result = cached_adaptive_engine(
                processed_ppg=processed_ppg,
                N_past=N_past,
                T_min=T_min,
                T_max=T_max,
                O_min=O_min,
                O_max=O_max,
                theta=theta
            )

        # Step 4: Feature Extraction and Aggregations
        with st.spinner("Extracting morphological features and benchmark averages..."):
            feature_result = cached_feature_extraction(
                processed_ppg=processed_ppg,
                adaptive_result=adaptive_result,
                T_crit_macro=T_crit_macro,
                T_crit_derivatives=T_crit_derivatives,
                fixed_window_sec=fixed_window_sec
            )

    except BioSignalLoaderError as e:
        st.error(f"❌ Error loading file: {e}")
        return
    except PPGPreprocessingError as e:
        st.error(f"❌ Preprocessing and segmentation error: {e}")
        return
    except AdaptiveEngineError as e:
        st.error(f"❌ Adaptive engine computation error: {e}")
        return
    except FeatureExtractionError as e:
        st.error(f"❌ Feature extraction error: {e}")
        return
    except Exception as e:
        st.error(f"❌ An unexpected error occurred: {e}")
        return

    # Check for empty beat outputs
    n_beats = len(processed_ppg.beats_raw)
    if n_beats == 0:
        st.warning("⚠️ No heartbeats were detected in the signal. Check the channel index or sampling rate.")
        return

    # ----------------- MAIN PANEL RENDERING -----------------
    # Status Badge
    sqi_val = "High" if np.sum(processed_ppg.valid_beats_mask) / n_beats > 0.8 else "Medium/Low"
    st.markdown(
        f"<div class='badge-status'>Frequency: {signal_data.sampling_rate} Hz | "
        f"Beats Detected: {n_beats} | SQI Quality: {sqi_val}</div>",
        unsafe_allow_html=True
    )

    # KPI Summary Row
    kpi1, kpi2, kpi3, kpi4 = st.columns(4)
    kpi1.metric(
        label="Detected Beats",
        value=f"{n_beats} (Valid: {np.sum(processed_ppg.valid_beats_mask)})"
    )
    
    mean_entropy = float(np.nanmean(adaptive_result.entropy_normalized)) if len(adaptive_result.entropy_normalized) > 0 else 0.0
    kpi2.metric(
        label="Mean H_morph (Complexity)",
        value=f"{mean_entropy:.3f}"
    )
    
    mean_tw = float(np.nanmean(adaptive_result.target_window_beats)) if len(adaptive_result.target_window_beats) > 0 else 0.0
    kpi3.metric(
        label="Mean Tw Window (Beats)",
        value=f"{mean_tw:.1f}",
        delta=f"Fixed benchmark: {fixed_window_sec} s"
    )
    
    # Calculate transient resolution gain proxy
    beats_per_30s = 30.0 / (processed_ppg.beat_times[1] - processed_ppg.beat_times[0]) if len(processed_ppg.beat_times) > 1 else 37.5
    gain_pct = max(0.0, (1.0 - (mean_tw / beats_per_30s)) * 100.0)
    kpi4.metric(
        label="Transient Resolution (Gain)",
        value=f"+{gain_pct:.1f}%",
        delta="vs Fixed Window"
    )

    # Tabs
    tab1, tab2, tab3 = st.tabs([
        "📈 Preprocessed Signal & Beats",
        "🔄 Entropy Dynamics & Tw Windows",
        "📊 Extracted Features & Benchmark (Adaptive vs Fixed)"
    ])

    # TAB 1: PREPROCESSED SIGNAL & BEATS
    with tab1:
        fig_signal = PPGVisualizer.plot_signal_and_beats(processed_ppg)
        if fig_signal is not None:
            st.plotly_chart(fig_signal, use_container_width=True)
        else:
            st.write("Signal visualization is not available.")

        st.markdown("### Normalized Aligned Beats Matrix")
        if len(processed_ppg.beats_normalized) > 0:
            avg_beat = np.mean(processed_ppg.beats_normalized, axis=0)
            df_beats = pd.DataFrame({
                "Sample Point": np.arange(len(avg_beat)),
                "Average BeatTemplate": avg_beat
            })
            for k in range(min(5, len(processed_ppg.beats_normalized))):
                df_beats[f"Beat {k}"] = processed_ppg.beats_normalized[k]
            st.line_chart(df_beats.set_index("Sample Point"))
        else:
            st.write("Segmentation data is not available.")

    # TAB 2: ENTROPY DYNAMICS & ADAPTIVE WINDOWS
    with tab2:
        fig_entropy = PPGVisualizer.plot_entropy_and_windows(adaptive_result, theta_sensitivity=theta)
        if fig_entropy is not None:
            st.plotly_chart(fig_entropy, use_container_width=True)
        else:
            st.write("Entropy dynamics visualization is not available.")

    # TAB 3: EXTRACTED FEATURES & BENCHMARK
    with tab3:
        st.subheader("Beat-by-Beat Morphological Features")
        st.dataframe(feature_result.single_beat_df)

        # Plotly comparative feature chart
        st.markdown("### Transient Features Dynamics (Adaptive vs Fixed)")
        fig_compare = PPGVisualizer.plot_feature_comparison(feature_result, "max_systolic_slope")
        if fig_compare is not None:
            st.plotly_chart(fig_compare, use_container_width=True)

        col_left, col_right = st.columns(2)
        with col_left:
            st.markdown("#### Aggregated over Adaptive Windows (Tw)")
            st.dataframe(feature_result.adaptive_features_df)
        with col_right:
            st.markdown("#### Aggregated over Fixed Windows (30s)")
            st.dataframe(feature_result.fixed_features_df)

        st.markdown("### Comparative Time Series Summary Metrics")
        st.write(feature_result.comparison_metrics)

        # Render export options using the consolidated component
        config_log = {
            "lowcut": lowcut,
            "highcut": highcut,
            "target_beat_length": target_length,
            "N_past": N_past,
            "T_min": T_min,
            "T_max": T_max,
            "O_min": O_min,
            "O_max": O_max,
            "theta": theta,
            "T_crit_macro": T_crit_macro,
            "T_crit_derivatives": T_crit_derivatives,
            "fixed_window_sec": fixed_window_sec
        }
        render_export_ui(feature_result, adaptive_result, config_log)


if __name__ == '__main__':
    main()
