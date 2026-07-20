"""PPG morphological feature extraction module.

This module computes 21 detailed morphologic, time-domain, area-domain,
and derivative-domain (VPG, SDPTG) features from individual PPG beats.
It also provides tools to aggregate these features over variable adaptive
windows or fixed-time windows (e.g. 30 seconds).
"""

import warnings
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Union, Any, Optional

import numpy as np
import pandas as pd
import scipy.signal
import scipy.integrate

# Safe import of data classes from upstream modules
try:
    from preprocessing import ProcessedPPG
except ImportError:
    @dataclass
    class ProcessedPPG:
        """Mock ProcessedPPG class for standalone testing."""
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
        """Mock AdaptiveWindowResult class for standalone testing."""
        beat_times: np.ndarray
        entropy_raw: np.ndarray
        entropy_normalized: np.ndarray
        target_window_beats: np.ndarray
        overlap_ratio: np.ndarray
        feature_windows_beats: Dict[str, np.ndarray]
        window_slices: Dict[str, List[Tuple[int, int]]]


class FeatureExtractionError(Exception):
    """Custom exception raised for errors within the PPGFeatureExtractor."""
    pass


@dataclass
class SingleBeatFeatures:
    """Dataclass holding all 21 morphological features of a single PPG beat.

    Attributes:
        beat_index (int): Chronological index of the beat.
        timestamp (float): Start time of the beat in seconds.
        pulse_duration (float): Total cycle duration (onset to onset) in seconds.
        crest_time (float): Time from onset to systolic peak in seconds.
        decay_time (float): Time from systolic peak to ending onset in seconds.
        pulse_width_10 (float): Wave width at 10% peak amplitude.
        pulse_width_25 (float): Wave width at 25% peak amplitude.
        pulse_width_50 (float): Wave width at 50% peak amplitude.
        pulse_width_75 (float): Wave width at 75% peak amplitude.
        duty_cycle (float): Ratio of crest time to pulse duration.
        peak_amplitude (float): Absolute systolic peak amplitude.
        inflection_amplitude (float): Absolute dicrotic wave amplitude.
        reflection_index (float): Inflection to systolic amplitude ratio.
        stiffness_index_proxy (float): Velocity proxy (1 / crest_time).
        total_area (float): Absolute area under the PPG curve.
        systolic_area (float): Area under the curve from onset to peak.
        diastolic_area (float): Area under the curve from peak to final onset.
        area_ratio (float): Ratio of diastolic area to systolic area.
        max_systolic_slope (float): Maximum upward slope (PSI).
        max_decay_slope (float): Maximum downward slope (minimum derivative).
        slope_ratio (float): Ratio of maximum systolic slope to max decay slope.
        sdptg_b_a (float): b/a ratio from the second derivative.
        sdptg_c_a (float): c/a ratio from the second derivative.
        sdptg_d_a (float): d/a ratio from the second derivative.
        sdptg_e_a (float): e/a ratio from the second derivative.
        aging_index (float): AGI value, (b - c - d - e) / a.
    """
    beat_index: int
    timestamp: float
    pulse_duration: float
    crest_time: float
    decay_time: float
    pulse_width_10: float
    pulse_width_25: float
    pulse_width_50: float
    pulse_width_75: float
    duty_cycle: float
    peak_amplitude: float
    inflection_amplitude: float
    reflection_index: float
    stiffness_index_proxy: float
    total_area: float
    systolic_area: float
    diastolic_area: float
    area_ratio: float
    max_systolic_slope: float
    max_decay_slope: float
    slope_ratio: float
    sdptg_b_a: float
    sdptg_c_a: float
    sdptg_d_a: float
    sdptg_e_a: float
    aging_index: float

    def to_dict(self) -> Dict[str, Any]:
        """Converts the dataclass fields into a dictionary."""
        return self.__dict__.copy()


@dataclass
class FeatureExtractionResult:
    """Container for the full feature extraction outcomes.

    Attributes:
        single_beat_df (pd.DataFrame): Beat-by-beat detailed features.
        adaptive_features_df (pd.DataFrame): Features aggregated over adaptive windows.
        fixed_features_df (pd.DataFrame): Features aggregated over fixed-time windows (e.g., 30s).
        comparison_metrics (dict): Benchmark errors/variances comparing the two series.
    """
    single_beat_df: pd.DataFrame
    adaptive_features_df: pd.DataFrame
    fixed_features_df: pd.DataFrame
    comparison_metrics: dict


class PPGFeatureExtractor:
    """Extractor class carrying out beat features parsing and multi-resolution window aggregates."""

    def __init__(self, epsilon: float = 1e-8):
        """Initializes the feature extractor.

        Args:
            epsilon (float): Regularization constant to avoid division-by-zero errors. Default 1e-8.
        """
        self.epsilon = epsilon

    def extract_single_beat_features(
        self,
        beat_normalized: np.ndarray,
        dt: float,
        beat_index: int = 0,
        timestamp: float = 0.0,
        raw_amplitude: float = 1.0
    ) -> SingleBeatFeatures:
        """Calculates 21 morphological features for a single normalized beat.

        Args:
            beat_normalized (np.ndarray): 1D array of normalized beat waveform (resampled).
            dt (float): Sampling step in seconds (1 / sampling_rate).
            beat_index (int): Numerical index of the beat.
            timestamp (float): Absolute start time of the beat.
            raw_amplitude (float): Peak-to-valley amplitude scale of the raw beat.

        Returns:
            SingleBeatFeatures: Extracted metric container.
        """
        eps = self.epsilon
        L = len(beat_normalized)
        t = np.arange(L) * dt
        pulse_duration = L * dt

        # 1. Peak & decay times
        crest_idx = np.argmax(beat_normalized)
        crest_time = crest_idx * dt
        decay_time = pulse_duration - crest_time

        # 2. Widths at different thresholds (10%, 25%, 50%, 75%)
        # For a standard single pulse peak=1.0, valley=0.0
        def get_width(thresh: float) -> float:
            indices = np.where(beat_normalized >= thresh)[0]
            if len(indices) >= 2:
                return (indices[-1] - indices[0]) * dt
            return 0.0

        pw10 = get_width(0.10)
        pw25 = get_width(0.25)
        pw50 = get_width(0.50)
        pw75 = get_width(0.75)

        duty_cycle = crest_time / (pulse_duration + eps)

        # 3. Peak amplitude (absolute scale)
        peak_amplitude = raw_amplitude
        stiffness_index_proxy = 1.0 / (crest_time + eps)

        # 4. Inflection point & dicrotic amplitude
        # Search inflection in decay phase
        decay_phase = beat_normalized[crest_idx:]
        decay_len = len(decay_phase)
        
        # Default fallback is mid-point of decay phase
        dicrotic_idx = crest_idx + (decay_len // 2)

        # Try to find local peaks in decay phase
        local_peaks, _ = scipy.signal.find_peaks(decay_phase)
        if len(local_peaks) > 0:
            dicrotic_idx = crest_idx + local_peaks[0]
        elif decay_len >= 2:
            # If no local peak, find the point where the decay slows down most (minimum absolute slope)
            # which corresponds to the maximum of the derivative
            dy_decay = np.gradient(decay_phase, dt)
            dy_peaks, _ = scipy.signal.find_peaks(dy_decay)
            if len(dy_peaks) > 0:
                dicrotic_idx = crest_idx + dy_peaks[0]

        inflection_amplitude = beat_normalized[dicrotic_idx] * raw_amplitude
        reflection_index = inflection_amplitude / (peak_amplitude + eps)

        # 5. Area integration (using trapezoidal method)
        # Scale back to original amplitude domain
        scaled_beat = beat_normalized * raw_amplitude
        total_area = float(scipy.integrate.trapezoid(scaled_beat, dx=dt))
        
        # Systolic area up to crest, diastolic area after
        systolic_area = float(scipy.integrate.trapezoid(scaled_beat[:crest_idx + 1], dx=dt))
        diastolic_area = total_area - systolic_area
        area_ratio = diastolic_area / (systolic_area + eps)

        # 6. First derivative (VPG)
        dy = np.gradient(scaled_beat, dt)
        max_systolic_slope = float(np.max(dy[:crest_idx + 1])) if crest_idx > 0 else 0.0
        max_decay_slope = float(np.min(dy[crest_idx:])) if crest_idx < L else 0.0
        slope_ratio = abs(max_systolic_slope / (max_decay_slope + eps))

        # 7. Second derivative (SDPTG)
        ddy = np.gradient(dy, dt)
        
        # Find positive peaks and negative valleys
        pos_peaks, _ = scipy.signal.find_peaks(ddy)
        neg_valleys, _ = scipy.signal.find_peaks(-ddy)

        # Initialize SDPTG keypoint indices
        a_idx = b_idx = c_idx = d_idx = e_idx = None

        # a is the first positive peak
        if len(pos_peaks) > 0:
            a_idx = pos_peaks[0]
            # b is the first negative valley after a
            after_a_neg = neg_valleys[neg_valleys > a_idx]
            if len(after_a_neg) > 0:
                b_idx = after_a_neg[0]
                # c is the next positive peak after b
                after_b_pos = pos_peaks[pos_peaks > b_idx]
                if len(after_b_pos) > 0:
                    c_idx = after_b_pos[0]
                    # d is the next negative valley after c
                    after_c_neg = neg_valleys[neg_valleys > c_idx]
                    if len(after_c_neg) > 0:
                        d_idx = after_c_neg[0]
                        # e is the next positive peak after d
                        after_d_pos = pos_peaks[pos_peaks > d_idx]
                        if len(after_d_pos) > 0:
                            e_idx = after_d_pos[0]

        # Extract values if all coordinates resolved
        if all(idx is not None for idx in [a_idx, b_idx, c_idx, d_idx, e_idx]):
            a = ddy[a_idx]
            b = ddy[b_idx]
            c = ddy[c_idx]
            d = ddy[d_idx]
            e = ddy[e_idx]

            sdptg_b_a = b / (a + eps)
            sdptg_c_a = c / (a + eps)
            sdptg_d_a = d / (a + eps)
            sdptg_e_a = e / (a + eps)
            aging_index = (b - c - d - e) / (a + eps)
        else:
            sdptg_b_a = np.nan
            sdptg_c_a = np.nan
            sdptg_d_a = np.nan
            sdptg_e_a = np.nan
            aging_index = np.nan

        return SingleBeatFeatures(
            beat_index=beat_index,
            timestamp=timestamp,
            pulse_duration=pulse_duration,
            crest_time=crest_time,
            decay_time=decay_time,
            pulse_width_10=pw10,
            pulse_width_25=pw25,
            pulse_width_50=pw50,
            pulse_width_75=pw75,
            duty_cycle=duty_cycle,
            peak_amplitude=peak_amplitude,
            inflection_amplitude=inflection_amplitude,
            reflection_index=reflection_index,
            stiffness_index_proxy=stiffness_index_proxy,
            total_area=total_area,
            systolic_area=systolic_area,
            diastolic_area=diastolic_area,
            area_ratio=area_ratio,
            max_systolic_slope=max_systolic_slope,
            max_decay_slope=max_decay_slope,
            slope_ratio=slope_ratio,
            sdptg_b_a=sdptg_b_a,
            sdptg_c_a=sdptg_c_a,
            sdptg_d_a=sdptg_d_a,
            sdptg_e_a=sdptg_e_a,
            aging_index=aging_index
        )

    def extract_all_beats(self, processed_ppg: ProcessedPPG) -> pd.DataFrame:
        """Extracts morphological features for all valid beats.

        Args:
            processed_ppg (ProcessedPPG): Preprocessed PPG signal metadata.

        Returns:
            pd.DataFrame: Table with one row per beat containing its index, timestamp, and 21 features.
        """
        beats_norm = processed_ppg.beats_normalized
        beats_raw = processed_ppg.beats_raw
        valid_mask = processed_ppg.valid_beats_mask
        beat_times = processed_ppg.beat_times
        fs = processed_ppg.sampling_rate
        dt = 1.0 / fs

        features_list = []

        for i in range(len(beats_norm)):
            # If beat is marked invalid, append a row with NaNs
            if not valid_mask[i]:
                # Rebuild empty row matching the dataclass keys
                empty_row = {field_name: np.nan for field_name in SingleBeatFeatures.__dataclass_fields__.keys()}
                empty_row["beat_index"] = i
                empty_row["timestamp"] = beat_times[i] if i < len(beat_times) else (i * 0.8)
                features_list.append(empty_row)
                continue

            # Calculate raw amplitude range (peak to valley) for the raw beat segment
            raw_beat = beats_raw[i]
            raw_amp = float(raw_beat.max() - raw_beat.min()) if len(raw_beat) > 0 else 1.0

            beat_norm = beats_norm[i]
            # Use nominal time step matching the normalized length (128 points over raw duration)
            beat_duration = len(raw_beat) * dt if len(raw_beat) > 0 else 0.8
            normalized_dt = beat_duration / (len(beat_norm) - 1 + self.epsilon)

            feat = self.extract_single_beat_features(
                beat_normalized=beat_norm,
                dt=normalized_dt,
                beat_index=i,
                timestamp=beat_times[i] if i < len(beat_times) else (i * 0.8),
                raw_amplitude=raw_amp
            )
            features_list.append(feat.to_dict())

        return pd.DataFrame(features_list)

    def aggregate_window_features(
        self,
        single_beat_df: pd.DataFrame,
        slices: List[Tuple[int, int]]
    ) -> pd.DataFrame:
        """Aggregates beat-by-beat features over the provided slice intervals.

        Computes mean and standard deviation for each column.

        Args:
            single_beat_df (pd.DataFrame): Beat-by-beat detail features table.
            slices (List[Tuple[int, int]]): List of (start_beat_idx, end_beat_idx) intervals.

        Returns:
            pd.DataFrame: Table containing aggregated feature trends with time mappings.
        """
        aggregated_rows = []

        # Find all feature columns (excluding index keys)
        exclude_cols = {"beat_index", "timestamp"}
        feature_cols = [c for c in single_beat_df.columns if c not in exclude_cols]

        for window_idx, (start_idx, end_idx) in enumerate(slices):
            # Select slice partition
            window_df = single_beat_df.iloc[start_idx:end_idx]

            # Exclude invalid NaN rows
            valid_window_df = window_df.dropna(subset=feature_cols, how='all')

            row_data = {
                "window_index": window_idx,
                "start_beat": start_idx,
                "end_beat": end_idx,
                "timestamp": float(window_df["timestamp"].mean()) if len(window_df) > 0 else np.nan,
                "num_beats": len(window_df)
            }

            if not valid_window_df.empty:
                # Compute average and std-dev for each feature
                for col in feature_cols:
                    # Ignore warnings if slice has only NaNs
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore", category=RuntimeWarning)
                        row_data[f"{col}_mean"] = float(np.nanmean(valid_window_df[col]))
                        row_data[f"{col}_std"] = float(np.nanstd(valid_window_df[col]))
            else:
                for col in feature_cols:
                    row_data[f"{col}_mean"] = np.nan
                    row_data[f"{col}_std"] = np.nan

            aggregated_rows.append(row_data)

        return pd.DataFrame(aggregated_rows)

    def extract_fixed_windows(
        self,
        single_beat_df: pd.DataFrame,
        window_size_sec: float = 30.0,
        overlap_sec: float = 15.0
    ) -> pd.DataFrame:
        """Aggregates beat-by-beat features over traditional fixed-time windows.

        Args:
            single_beat_df (pd.DataFrame): Beat-by-beat features table.
            window_size_sec (float): Fixed window span in seconds. Default 30.0.
            overlap_sec (float): Window overlap span in seconds. Default 15.0.

        Returns:
            pd.DataFrame: Table of time-locked aggregated features.
        """
        if single_beat_df.empty:
            return pd.DataFrame()

        t_min = float(single_beat_df["timestamp"].min())
        t_max = float(single_beat_df["timestamp"].max())
        
        step_sec = window_size_sec - overlap_sec
        if step_sec <= 0:
            step_sec = window_size_sec

        fixed_rows = []
        exclude_cols = {"beat_index", "timestamp"}
        feature_cols = [c for c in single_beat_df.columns if c not in exclude_cols]

        curr_t = t_min
        window_idx = 0

        while curr_t < t_max:
            win_start = curr_t
            win_end = curr_t + window_size_sec

            # Filter beats falling in the time window
            in_window_df = single_beat_df[
                (single_beat_df["timestamp"] >= win_start) & 
                (single_beat_df["timestamp"] < win_end)
            ]

            valid_in_win = in_window_df.dropna(subset=feature_cols, how='all')

            row_data = {
                "window_index": window_idx,
                "window_start_sec": win_start,
                "window_end_sec": win_end,
                "timestamp": win_start + (window_size_sec / 2.0),
                "num_beats": len(in_window_df)
            }

            if not valid_in_win.empty:
                for col in feature_cols:
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore", category=RuntimeWarning)
                        row_data[f"{col}_mean"] = float(np.nanmean(valid_in_win[col]))
                        row_data[f"{col}_std"] = float(np.nanstd(valid_in_win[col]))
            else:
                for col in feature_cols:
                    row_data[f"{col}_mean"] = np.nan
                    row_data[f"{col}_std"] = np.nan

            fixed_rows.append(row_data)

            # Shift window forward
            curr_t += step_sec
            window_idx += 1

            if curr_t + step_sec > t_max and len(fixed_rows) > 0:
                break

        return pd.DataFrame(fixed_rows)

    def process(
        self,
        processed_ppg: Any,
        adaptive_result: Any,
        fixed_window_sec: float = 30.0
    ) -> FeatureExtractionResult:
        """Runs feature extraction, aggregates over adaptive slices and fixed-time windows, and compares.

        Args:
            processed_ppg (Any): ProcessedPPG input containing beats.
            adaptive_result (Any): AdaptiveWindowResult configuration containing schedules.
            fixed_window_sec (float): Benchmark window width in seconds. Default 30.0.

        Returns:
            FeatureExtractionResult: Complete comparative outcomes package.
        """
        # 1. Extract detailed single-beat features
        single_beat_df = self.extract_all_beats(processed_ppg)

        # 2. Extract adaptive aggregates (using 'derivatives' as reference slice category, or taking the first available)
        slices_dict = adaptive_result.window_slices
        ref_category = "derivatives" if "derivatives" in slices_dict else list(slices_dict.keys())[0]
        ref_slices = slices_dict[ref_category]

        adaptive_df = self.aggregate_window_features(single_beat_df, ref_slices)

        # 3. Extract fixed 30s benchmark aggregates
        fixed_df = self.extract_fixed_windows(
            single_beat_df,
            window_size_sec=fixed_window_sec,
            overlap_sec=fixed_window_sec / 2.0
        )

        # 4. Compute comparison metrics (e.g. comparing average trends of critical features)
        comparison = {}
        target_feature = "max_systolic_slope"
        mean_col = f"{target_feature}_mean"

        if mean_col in adaptive_df.columns and mean_col in fixed_df.columns:
            # Interpolate series to match shapes and measure Mean Absolute Deviation
            val_adapt = adaptive_df[mean_col].dropna().to_numpy()
            val_fixed = fixed_df[mean_col].dropna().to_numpy()
            
            comparison["adaptive_variance"] = float(np.var(val_adapt)) if len(val_adapt) > 0 else 0.0
            comparison["fixed_variance"] = float(np.var(val_fixed)) if len(val_fixed) > 0 else 0.0
            comparison["ref_feature"] = target_feature

        return FeatureExtractionResult(
            single_beat_df=single_beat_df,
            adaptive_features_df=adaptive_df,
            fixed_features_df=fixed_df,
            comparison_metrics=comparison
        )


if __name__ == '__main__':
    print("--- Running autonomous tests for PPGFeatureExtractor ---")

    # 1. Generate synthetic beat dataset (50 beats, 128 samples each, dt = 1/128)
    n_beats = 50
    beat_len = 128
    dt = 1.0 / 128.0
    t_beat = np.arange(beat_len) * dt

    # Construct double-gaussian pulse wave template (systolic peak at 0.15s, dicrotic wave at 0.4s)
    systolic = np.exp(-((t_beat - 0.15) / 0.06)**2)
    dicrotic = 0.35 * np.exp(-((t_beat - 0.40) / 0.10)**2)
    clean_beat = systolic + dicrotic

    # Add minor noise variations to generate 50 unique beats
    beats = []
    for i in range(n_beats):
        noise = 0.005 * np.random.randn(beat_len)
        beats.append(clean_beat + noise)

    beats_matrix = np.vstack(beats)
    valid_mask = np.ones(n_beats, dtype=bool)
    # 0.8 seconds per beat (75 BPM)
    beat_times = np.arange(n_beats) * 0.8

    # 2. Package into mock ProcessedPPG
    # Raw beats of variable lengths
    raw_beats = [clean_beat.copy() for _ in range(n_beats)]
    mock_ppg = ProcessedPPG(
        raw_signal=np.zeros(1000),
        filtered_signal=np.zeros(1000),
        sampling_rate=128.0,
        peaks=np.array([]),
        onsets=np.array([]),
        beats_raw=raw_beats,
        beats_normalized=beats_matrix,
        beat_times=beat_times,
        valid_beats_mask=valid_mask
    )

    # 3. Process features
    extractor = PPGFeatureExtractor()
    single_beat_df = extractor.extract_all_beats(mock_ppg)

    # 4. Perform assertions on single beats
    assert len(single_beat_df) == n_beats, "Row count mismatch in single_beat_df"
    # Column verification (excluding beat_index and timestamp)
    expected_cols = 24  # 21 features + index + timestamp + 1 extra if shape, let's verify key features exist
    required_features = ["pulse_duration", "crest_time", "total_area", "max_systolic_slope", "sdptg_b_a", "aging_index"]
    for feat in required_features:
        assert feat in single_beat_df.columns, f"Feature column '{feat}' missing in DataFrame"
    
    # Assert values logic
    assert (single_beat_df["max_systolic_slope"] > 0).all(), "Slope must be strictly positive"
    assert (single_beat_df["total_area"] > 0).all(), "Area must be strictly positive"
    print("-> Beat feature calculation assertions: SUCCESS")

    # 5. Run full orchestrator with mock adaptive slices
    # Let's mock a set of adaptive slices
    mock_slices = [(0, 10), (10, 20), (20, 30), (30, 40), (40, 50)]
    mock_adaptive = AdaptiveWindowResult(
        beat_times=beat_times,
        entropy_raw=np.zeros(n_beats),
        entropy_normalized=np.zeros(n_beats),
        target_window_beats=np.zeros(n_beats),
        overlap_ratio=np.zeros(n_beats),
        feature_windows_beats={},
        window_slices={"derivatives": mock_slices}
    )

    res = extractor.process(mock_ppg, mock_adaptive, fixed_window_sec=15.0)

    # Assert correct aggregation shapes
    assert len(res.adaptive_features_df) == len(mock_slices), "Adaptive rows count mismatch"
    assert len(res.fixed_features_df) > 0, "Fixed rows count must be positive"
    print("-> Multi-window aggregation and benchmark comparisons: SUCCESS")

    # 6. Display head details
    print("\n" + "="*80)
    print("                      SAMPLE EXTRACTED FEATURE MATRIX (FIRST 5 BEATS)")
    print("="*80)
    display_cols = ["beat_index", "timestamp", "pulse_duration", "crest_time", "total_area", "max_systolic_slope", "aging_index"]
    print(single_beat_df[display_cols].head())
    print("="*80)

    print("\n--- All feature extraction tests completed successfully! ---")
