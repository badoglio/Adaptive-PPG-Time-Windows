"""PPG analysis visualization module.

This module provides the `PPGVisualizer` class to generate interactive
Plotly figures for raw/filtered signal traces, beat detection markers,
entropy windowing curves, and adaptive vs. fixed window feature comparisons.
It also includes utility exports for Streamlit dashboard reports.
"""

import io
import json
import warnings
from typing import Dict, List, Tuple, Union, Any, Optional

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px
import streamlit as st

# Safe import of data structures from preceding modules
try:
    from preprocessing import ProcessedPPG
except ImportError:
    from dataclasses import dataclass, field
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

try:
    from adaptive_engine import AdaptiveWindowResult
except ImportError:
    @dataclass
    class AdaptiveWindowResult:
        beat_times: np.ndarray
        entropy_raw: np.ndarray
        entropy_normalized: np.ndarray
        target_window_beats: np.ndarray
        overlap_ratio: np.ndarray
        feature_windows_beats: Dict[str, np.ndarray]
        window_slices: Dict[str, List[Tuple[int, int]]]

try:
    from feature_extraction import FeatureExtractionResult
except ImportError:
    @dataclass
    class FeatureExtractionResult:
        single_beat_df: pd.DataFrame
        adaptive_features_df: pd.DataFrame
        fixed_features_df: pd.DataFrame
        comparison_metrics: dict


class VisualizationError(Exception):
    """Custom exception raised for errors within the PPGVisualizer."""
    pass


class PPGVisualizer:
    """Plotly rendering engine for biosegnal analysis dashboards."""

    @staticmethod
    def _create_empty_figure(message: str) -> go.Figure:
        """Helper to generate an empty figure showing a central alert annotation.

        Args:
            message (str): Text to show in the center of the chart.

        Returns:
            go.Figure: Interactive blank canvas.
        """
        fig = go.Figure()
        fig.add_annotation(
            text=message,
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False,
            font=dict(size=16, color="gray")
        )
        fig.update_layout(
            template="plotly_white",
            xaxis=dict(showticklabels=False, showgrid=False, zeroline=False),
            yaxis=dict(showticklabels=False, showgrid=False, zeroline=False)
        )
        return fig

    @classmethod
    def plot_signal_and_beats(cls, processed_ppg: ProcessedPPG) -> go.Figure:
        """Generates a 2-subplot chart showing raw vs filtered signal and peak-valley markers.

        High-performance synchronization allows matched zooming across both axes.

        Args:
            processed_ppg (ProcessedPPG): Preprocessed PPG signal.

        Returns:
            go.Figure: Synchronized Plotly trace.
        """
        if len(processed_ppg.raw_signal) == 0:
            return cls._create_empty_figure("No signal data available.")

        fs = processed_ppg.sampling_rate
        t = np.arange(len(processed_ppg.raw_signal)) / fs

        fig = make_subplots(
            rows=2, cols=1,
            shared_xaxes=True,
            vertical_spacing=0.08,
            subplot_titles=(
                "Raw vs Filtered PPG (0.5 - 30 Hz)",
                "Detected Systolic Peaks & Beat Onsets"
            )
        )

        # Subplot 1: Raw vs Filtered
        fig.add_trace(
            go.Scatter(x=t, y=processed_ppg.raw_signal, name="Raw", line=dict(color="#94a3b8", width=1)),
            row=1, col=1
        )
        fig.add_trace(
            go.Scatter(x=t, y=processed_ppg.filtered_signal, name="Filtered", line=dict(color="#3b82f6", width=1.5)),
            row=1, col=1
        )

        # Subplot 2: Filtered with Peak/Onset markers
        fig.add_trace(
            go.Scatter(x=t, y=processed_ppg.filtered_signal, name="PPG Wave", line=dict(color="#1d4ed8", width=1.5)),
            row=2, col=1
        )

        # Add systolic peaks
        peaks = processed_ppg.peaks
        if len(peaks) > 0:
            fig.add_trace(
                go.Scatter(
                    x=t[peaks], y=processed_ppg.filtered_signal[peaks],
                    mode="markers",
                    name="Peak (Systolic)",
                    marker=dict(color="#22c55e", size=8, symbol="circle")
                ),
                row=2, col=1
            )

        # Add onsets/feet
        onsets = processed_ppg.onsets
        if len(onsets) > 0:
            fig.add_trace(
                go.Scatter(
                    x=t[onsets], y=processed_ppg.filtered_signal[onsets],
                    mode="markers",
                    name="Onset (Foot)",
                    marker=dict(color="#ef4444", size=8, symbol="triangle-up")
                ),
                row=2, col=1
            )

        # Highlight invalid beats with gray background rectangles
        valid_mask = processed_ppg.valid_beats_mask
        for i in range(len(valid_mask)):
            if not valid_mask[i] and i+1 < len(onsets):
                t_start = onsets[i] / fs
                t_end = onsets[i+1] / fs
                fig.add_vrect(
                    x0=t_start, x1=t_end,
                    fillcolor="#6b7280", opacity=0.15,
                    layer="below", line_width=0,
                    row=2, col=1
                )

        fig.update_layout(
            template="plotly_white",
            height=600,
            hovermode="x unified",
            xaxis2_title="Time (seconds)",
            yaxis_title="Amplitude",
            yaxis2_title="Amplitude",
            legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1)
        )

        return fig

    @classmethod
    def plot_entropy_and_windows(
        cls,
        adaptive_result: AdaptiveWindowResult,
        theta_sensitivity: float = 0.5
    ) -> go.Figure:
        """Generates a 3-subplot synchronized chart showing H_morph, Tw, and Ow curves.

        Args:
            adaptive_result (AdaptiveWindowResult): Sizing curves from engine.
            theta_sensitivity (float): Inflection threshold to highlight. Default 0.5.

        Returns:
            go.Figure: Synchronized Plotly trace.
        """
        if len(adaptive_result.beat_times) == 0:
            return cls._create_empty_figure("No beats analyzed for complexity.")

        beat_times = adaptive_result.beat_times
        
        fig = make_subplots(
            rows=3, cols=1,
            shared_xaxes=True,
            vertical_spacing=0.08,
            subplot_titles=(
                "Normalized Morphological Complexity H_morph[n]",
                "Adaptive Window Length Tw[n] per Feature Category",
                "Adaptive Overlap Ratio Ow[n]"
            )
        )

        # Subplot 1: H_morph
        fig.add_trace(
            go.Scatter(x=beat_times, y=adaptive_result.entropy_normalized, name="H_morph", line=dict(color="#8b5cf6", width=2)),
            row=1, col=1
        )
        # Add threshold line
        fig.add_shape(
            type="line",
            x0=beat_times[0], x1=beat_times[-1],
            y0=theta_sensitivity, y1=theta_sensitivity,
            line=dict(color="#ec4899", width=1.5, dash="dash"),
            row=1, col=1
        )

        # Subplot 2: Windows Tw curves per category
        colors_dict = {'macro': '#10b981', 'time_volume': '#3b82f6', 'derivatives': '#f43f5e'}
        for category, curve in adaptive_result.feature_windows_beats.items():
            color = colors_dict.get(category, '#6b7280')
            fig.add_trace(
                go.Scatter(x=beat_times, y=curve, name=f"Window: {category.replace('_', ' ').title()}", line=dict(color=color, width=1.8)),
                row=2, col=1
            )

        # Subplot 3: Overlap Ow percentage
        fig.add_trace(
            go.Scatter(x=beat_times, y=adaptive_result.overlap_ratio * 100.0, name="Overlap Ow (%)", line=dict(color="#f59e0b", width=2)),
            row=3, col=1
        )

        fig.update_layout(
            template="plotly_white",
            height=700,
            hovermode="x unified",
            xaxis3_title="Time (seconds)",
            yaxis_title="H_morph (0-1)",
            yaxis2_title="Beats (N)",
            yaxis3_title="Overlap (%)",
            legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1)
        )

        return fig

    @classmethod
    def plot_feature_comparison(
        cls,
        extraction_result: FeatureExtractionResult,
        feature_name: str = 'max_systolic_slope'
    ) -> go.Figure:
        """Compares adaptive variables trends to standard 30s fixed window profiles.

        Includes filled semitransparent bands showing standard deviations.

        Args:
            extraction_result (FeatureExtractionResult): Extraction outcome container.
            feature_name (str): Label of the target metric. Default 'max_systolic_slope'.

        Returns:
            go.Figure: Benchmarking Plotly trace.
        """
        adaptive_df = extraction_result.adaptive_features_df
        fixed_df = extraction_result.fixed_features_df

        if adaptive_df.empty and fixed_df.empty:
            return cls._create_empty_figure("No feature data available.")

        fig = go.Figure()

        # Helper to plot line and shaded error band
        def add_series_with_band(df, color, label):
            mean_col = f"{feature_name}_mean"
            std_col = f"{feature_name}_std"
            if mean_col not in df.columns or std_col not in df.columns:
                return

            sub = df.dropna(subset=[mean_col, std_col])
            if sub.empty:
                return

            t_val = sub["timestamp"].to_numpy()
            mean_val = sub[mean_col].to_numpy()
            std_val = sub[std_col].to_numpy()

            # Error band shade trace
            fig.add_trace(go.Scatter(
                x=np.concatenate([t_val, t_val[::-1]]),
                y=np.concatenate([mean_val + std_val, (mean_val - std_val)[::-1]]),
                fill='toself',
                fillcolor=color.replace("1)", "0.15)"),
                line=dict(color='rgba(255,255,255,0)'),
                hoverinfo="skip",
                showlegend=False,
                name=f"{label} Std-Dev"
            ))

            # Main average line
            fig.add_trace(go.Scatter(
                x=t_val, y=mean_val,
                name=label,
                line=dict(color=color, width=2)
            ))

        # Add Adaptive (Light Blue/Sapphire)
        add_series_with_band(adaptive_df, "rgba(59, 130, 246, 1)", "Adaptive Window")
        # Add Fixed (Red/Coral)
        add_series_with_band(fixed_df, "rgba(239, 68, 68, 1)", "Fixed 30s Window")

        fig.update_layout(
            template="plotly_white",
            title=f"Aggregation Comparison: {feature_name.replace('_', ' ').title()}",
            xaxis_title="Time (seconds)",
            yaxis_title="Aggregated Value (Mean ± Std-Dev)",
            height=500,
            hovermode="x unified"
        )

        return fig

    @classmethod
    def plot_multi_feature_dashboard(
        cls,
        extraction_result: FeatureExtractionResult,
        selected_features: List[str]
    ) -> go.Figure:
        """Plots a comparative grid for multiple selected features.

        Args:
            extraction_result (FeatureExtractionResult): Extraction outcomes.
            selected_features (List[str]): List of metric labels.

        Returns:
            go.Figure: Synchronized comparative grid.
        """
        n_feats = len(selected_features)
        if n_feats == 0:
            return cls._create_empty_figure("No features selected.")

        fig = make_subplots(
            rows=n_feats, cols=1,
            shared_xaxes=True,
            vertical_spacing=0.06,
            subplot_titles=[f.replace('_', ' ').title() for f in selected_features]
        )

        for i, feat in enumerate(selected_features):
            row = i + 1
            # Add Adaptive line (no band to avoid cluttering in grid)
            ad_df = extraction_result.adaptive_features_df.dropna(subset=[f"{feat}_mean"])
            if not ad_df.empty:
                fig.add_trace(
                    go.Scatter(
                        x=ad_df["timestamp"], y=ad_df[f"{feat}_mean"],
                        name=f"{feat} (Adaptive)",
                        line=dict(color="#3b82f6", width=2),
                        showlegend=False
                    ),
                    row=row, col=1
                )

            # Add Fixed line
            fx_df = extraction_result.fixed_features_df.dropna(subset=[f"{feat}_mean"])
            if not fx_df.empty:
                fig.add_trace(
                    go.Scatter(
                        x=fx_df["timestamp"], y=fx_df[f"{feat}_mean"],
                        name=f"{feat} (Fixed 30s)",
                        line=dict(color="#ef4444", width=2, dash="dash"),
                        showlegend=False
                    ),
                    row=row, col=1
                )

        fig.update_layout(
            template="plotly_white",
            height=250 * n_feats,
            hovermode="x unified",
            xaxis_title="Time (seconds)"
        )

        return fig


# ==========================================
# STREAMLIT DATA EXPORTER COMPONENT
# ==========================================

def render_export_ui(
    extraction_result: FeatureExtractionResult,
    adaptive_result: AdaptiveWindowResult,
    config_dict: Optional[dict] = None
):
    """Renders a Streamlit widgets layout panel with CSV and JSON download triggers.

    Args:
        extraction_result (FeatureExtractionResult): Features outcome.
        adaptive_result (AdaptiveWindowResult): Slices and entropy curves.
        config_dict (Optional[dict]): Settings log.
    """
    with st.expander("📥 Export Data and Analysis Report", expanded=False):
        st.markdown("Select a format below to download the processed outcomes:")

        col1, col2, col3 = st.columns(3)

        # 1. Beat-by-Beat features
        df_single = extraction_result.single_beat_df
        csv_single = df_single.to_csv(index=False).encode('utf-8')
        col1.download_button(
            label="Download Beat-by-Beat CSV",
            data=csv_single,
            file_name="ppg_single_beat_features.csv",
            mime="text/csv",
            help="Download the complete matrix of 21 features computed beat-by-beat."
        )

        # 2. Consolidated window summary
        # Merge adaptive and fixed aggregates into a single printable comparison frame
        ad_summary = extraction_result.adaptive_features_df.copy()
        fx_summary = extraction_result.fixed_features_df.copy()

        # Add labels to distinguish
        ad_summary["aggregation_type"] = "Adaptive"
        fx_summary["aggregation_type"] = "Fixed_30s"
        df_combined = pd.concat([ad_summary, fx_summary], ignore_index=True)
        csv_combined = df_combined.to_csv(index=False).encode('utf-8')
        
        col2.download_button(
            label="Download Consolidated CSV",
            data=csv_combined,
            file_name="ppg_consolidated_aggregates.csv",
            mime="text/csv",
            help="Download the consolidated temporal averages comparing adaptive to fixed-interval trends."
        )

        # 3. JSON Configurations log
        report = {
            "execution_timestamp": pd.Timestamp.now().isoformat(),
            "summary_metrics": {
                "total_beats": len(df_single),
                "valid_beats": int(df_single["pulse_duration"].dropna().count()),
                "mean_morphological_entropy": float(np.nanmean(adaptive_result.entropy_normalized)) if len(adaptive_result.entropy_normalized) > 0 else 0.0,
                "mean_window_beats": float(np.nanmean(adaptive_result.target_window_beats)) if len(adaptive_result.target_window_beats) > 0 else 0.0
            },
            "configuration": config_dict or {}
        }
        json_bytes = json.dumps(report, indent=4).encode('utf-8')
        col3.download_button(
            label="Download JSON Configurations Log",
            data=json_bytes,
            file_name="ppg_analysis_log.json",
            mime="application/json",
            help="Download the configuration parameters log and overall descriptive stats."
        )


if __name__ == '__main__':
    print("--- Running autonomous tests for PPGVisualizer ---")

    # 1. Generate coherent synthetic structures
    fs = 100.0
    L = 2000
    t = np.arange(L) / fs
    sig_raw = np.sin(2 * np.pi * 1.0 * t) + 0.1 * np.random.randn(L)
    sig_filt = np.sin(2 * np.pi * 1.0 * t)

    # 10 peaks and 10 onsets
    peaks = np.arange(100, 1900, 200)
    onsets = np.arange(50, 1850, 200)
    valid_mask = np.ones(len(peaks), dtype=bool)
    valid_mask[3] = False  # Mark one beat as invalid for validation testing

    mock_ppg = ProcessedPPG(
        raw_signal=sig_raw,
        filtered_signal=sig_filt,
        sampling_rate=fs,
        peaks=peaks,
        onsets=onsets,
        beats_raw=[np.ones(10) for _ in range(len(peaks))],
        beats_normalized=np.ones((len(peaks), 128)),
        beat_times=t[onsets[:-1]],
        valid_beats_mask=valid_mask
    )

    mock_adaptive = AdaptiveWindowResult(
        beat_times=t[onsets[:-1]],
        entropy_raw=np.linspace(0.1, 0.9, len(peaks) - 1),
        entropy_normalized=np.linspace(0.1, 0.9, len(peaks) - 1),
        target_window_beats=np.linspace(30, 5, len(peaks) - 1),
        overlap_ratio=np.linspace(0.25, 0.85, len(peaks) - 1),
        feature_windows_beats={
            'macro': np.linspace(30, 5, len(peaks) - 1),
            'time_volume': np.linspace(30, 5, len(peaks) - 1),
            'derivatives': np.linspace(30, 10, len(peaks) - 1)
        },
        window_slices={'derivatives': [(0, 3), (3, 6)]}
    )

    # DataFrame mock features
    single_data = {
        "beat_index": np.arange(len(peaks) - 1),
        "timestamp": t[onsets[:-1]],
        "pulse_duration": np.ones(len(peaks) - 1) * 0.8,
        "crest_time": np.ones(len(peaks) - 1) * 0.15,
        "total_area": np.ones(len(peaks) - 1) * 5.0,
        "max_systolic_slope": np.linspace(10.0, 15.0, len(peaks) - 1),
        "sdptg_b_a": np.ones(len(peaks) - 1) * -0.8,
        "aging_index": np.ones(len(peaks) - 1) * -1.2
    }
    df_single = pd.DataFrame(single_data)

    ad_data = {
        "window_index": [0, 1],
        "timestamp": [2.0, 6.0],
        "max_systolic_slope_mean": [11.0, 14.0],
        "max_systolic_slope_std": [0.5, 0.4],
        "aging_index_mean": [-0.8, -0.8],
        "aging_index_std": [0.05, 0.05]
    }
    df_ad = pd.DataFrame(ad_data)

    fx_data = {
        "window_index": [0, 1],
        "timestamp": [3.0, 7.0],
        "max_systolic_slope_mean": [11.5, 13.8],
        "max_systolic_slope_std": [0.6, 0.5],
        "aging_index_mean": [-0.8, -0.8],
        "aging_index_std": [0.06, 0.06]
    }
    df_fx = pd.DataFrame(fx_data)

    mock_extraction = FeatureExtractionResult(
        single_beat_df=df_single,
        adaptive_features_df=df_ad,
        fixed_features_df=df_fx,
        comparison_metrics={"test": 1.0}
    )

    # 2. Run plot tests
    print("Testing plot generation...")
    fig1 = PPGVisualizer.plot_signal_and_beats(mock_ppg)
    fig2 = PPGVisualizer.plot_entropy_and_windows(mock_adaptive)
    fig3 = PPGVisualizer.plot_feature_comparison(mock_extraction, "max_systolic_slope")
    fig4 = PPGVisualizer.plot_multi_feature_dashboard(mock_extraction, ["max_systolic_slope", "aging_index"])

    assert isinstance(fig1, go.Figure), "fig1 type mismatch"
    assert isinstance(fig2, go.Figure), "fig2 type mismatch"
    assert isinstance(fig3, go.Figure), "fig3 type mismatch"
    assert isinstance(fig4, go.Figure), "fig4 type mismatch"
    print("-> Plot structures instantiation: SUCCESS")

    # 3. Write check to HTML
    fig1.write_html("test_plot.html")
    print("-> Successfully saved test plot trace 'test_plot.html'")

    print("\n--- All visualization tests completed successfully! ---")
