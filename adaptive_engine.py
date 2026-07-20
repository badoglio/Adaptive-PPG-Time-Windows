"""Morphological entropy calculation and adaptive windowing engine.

This module implements the calculation of Morphological Entropy (H_morph)
over segmented beat matrices, maps it to adaptive temporal windows (T_w)
and overlaps (O_w), applies feature-specific critical thresholds (T_crit),
and schedules variable-step index slices for feature extraction.
"""

import warnings
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Union, Any, Optional

import numpy as np
import scipy.stats
import scipy.signal

# Safe import of ProcessedPPG from preprocessing.py
try:
    from preprocessing import ProcessedPPG
except ImportError:
    # Mock ProcessedPPG for standalone test execution
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


class AdaptiveEngineError(Exception):
    """Custom exception raised for errors within the AdaptiveEngine."""
    pass


@dataclass
class AdaptiveConfig:
    """Configuration parameters for the adaptive windowing engine.

    Attributes:
        N_past (int): Number of past beats used to calculate H_morph. Default 10.
        T_min_beats (int): Minimum window length in number of beats. Default 3.
        T_max_beats (int): Maximum window length in number of beats. Default 30.
        O_min (float): Minimum overlap ratio (e.g. 0.25 for 25%). Default 0.25.
        O_max (float): Maximum overlap ratio (e.g. 0.85 for 85%). Default 0.85.
        theta_sensitivity (float): Inflection sensitivity threshold for entropy. Default 0.5.
        T_crit_dict (Dict[str, int]): Minimum valid beats threshold for each feature category.
            Defaults to {'macro': 3, 'time_volume': 3, 'derivatives': 10}.
    """
    N_past: int = 10
    T_min_beats: int = 3
    T_max_beats: int = 30
    O_min: float = 0.25
    O_max: float = 0.85
    theta_sensitivity: float = 0.5
    T_crit_dict: Dict[str, int] = field(
        default_factory=lambda: {'macro': 3, 'time_volume': 3, 'derivatives': 10}
    )


@dataclass
class AdaptiveWindowResult:
    """Output results of the adaptive windowing engine.

    Attributes:
        beat_times (np.ndarray): 1D array of timestamps in seconds for each beat.
        entropy_raw (np.ndarray): 1D array of raw morphological entropy values per beat.
        entropy_normalized (np.ndarray): 1D array of normalized H_morph values in [0, 1].
        target_window_beats (np.ndarray): 1D array of target window sizes in beats.
        overlap_ratio (np.ndarray): 1D array of overlap percentages per beat.
        feature_windows_beats (Dict[str, np.ndarray]): Dict of capped window curves per category.
        window_slices (Dict[str, List[Tuple[int, int]]]): Dict of start/end indices lists for each category.
    """
    beat_times: np.ndarray
    entropy_raw: np.ndarray
    entropy_normalized: np.ndarray
    target_window_beats: np.ndarray
    overlap_ratio: np.ndarray
    feature_windows_beats: Dict[str, np.ndarray]
    window_slices: Dict[str, List[Tuple[int, int]]]


class AdaptiveEngine:
    """Class containing methods to calculate entropy, map windows, and produce slices."""

    def __init__(self, config: Optional[AdaptiveConfig] = None):
        """Initializes the adaptive engine.

        Args:
            config (Optional[AdaptiveConfig]): Engine configurations. If None, uses default settings.
        """
        self.config = config if config is not None else AdaptiveConfig()

    def compute_morphological_entropy(
        self,
        beats_matrix: np.ndarray,
        valid_mask: np.ndarray
    ) -> np.ndarray:
        """Computes the Morphological Entropy (H_morph) over a rolling window of past valid beats.

        Uses the mean Root Mean Squared (RMS) consecutive difference of normalized beat waveforms
        as a robust divergence index.

        Args:
            beats_matrix (np.ndarray): 2D array of shape [n_beats, target_length] containing normalized beats.
            valid_mask (np.ndarray): 1D boolean array indicating valid beats.

        Returns:
            np.ndarray: 1D array of length n_beats containing normalized H_morph values in [0, 1].
        """
        n_beats = len(beats_matrix)
        entropy_raw = np.zeros(n_beats, dtype=np.float64)

        if n_beats < self.config.N_past:
            # Insufficient beats to establish baseline, return all zeros
            return entropy_raw

        # Calculate raw divergence for each beat position
        for i in range(n_beats):
            # Extract previous N_past beats up to index i
            start_window = max(0, i - self.config.N_past + 1)
            end_window = i + 1

            # Get indices within window that are valid
            win_mask = valid_mask[start_window:end_window]
            win_beats = beats_matrix[start_window:end_window][win_mask]

            if len(win_beats) < 2:
                # If we don't have at least 2 valid beats to compare, divergence is 0
                entropy_raw[i] = 0.0
                continue

            # Calculate consecutive Root Mean Square differences between successive beats
            diffs = []
            for j in range(1, len(win_beats)):
                rms_diff = np.sqrt(np.mean((win_beats[j] - win_beats[j-1])**2))
                diffs.append(rms_diff)

            entropy_raw[i] = np.mean(diffs)

        # Handle early innesco (priming) phase beats by assigning the first valid calculation's value
        first_valid_idx = self.config.N_past - 1
        if first_valid_idx < n_beats:
            initial_val = entropy_raw[first_valid_idx]
            entropy_raw[:first_valid_idx] = initial_val

        # Normalize raw entropy values into [0, 1] range
        raw_min = entropy_raw.min()
        raw_max = entropy_raw.max()

        if raw_max > raw_min:
            entropy_normalized = (entropy_raw - raw_min) / (raw_max - raw_min)
        else:
            entropy_normalized = np.zeros_like(entropy_raw)

        return entropy_normalized

    def map_entropy_to_window(self, entropy_norm: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Maps normalized entropy values to target temporal windows and overlap ratios.

        T_w decreases linearly as entropy increases (shorter windows for high entropy/transients).
        O_w increases linearly as entropy increases (higher overlap to trace fast transitions).

        Args:
            entropy_norm (np.ndarray): 1D array of normalized H_morph values in [0, 1].

        Returns:
            Tuple[np.ndarray, np.ndarray]: (target_window_beats, overlap_ratio)
        """
        # Linear decresing mapping for T_w:
        # T_w = T_max - (T_max - T_min) * H_morph
        t_span = self.config.T_max_beats - self.config.T_min_beats
        target_window_beats = self.config.T_max_beats - (t_span * entropy_norm)

        # Linear increasing mapping for O_w:
        # O_w = O_min + (O_max - O_min) * H_morph
        o_span = self.config.O_max - self.config.O_min
        overlap_ratio = self.config.O_min + (o_span * entropy_norm)

        return target_window_beats, overlap_ratio

    def apply_critical_sensitivity(self, target_window_beats: np.ndarray) -> Dict[str, np.ndarray]:
        """Applies feature-specific critical constraints (T_crit) and rounds to integers.

        T_w^k = max(T_target, T_crit^k)

        Args:
            target_window_beats (np.ndarray): 1D array of float target window sizes.

        Returns:
            Dict[str, np.ndarray]: Dict of integer window curves for 'macro', 'time_volume', and 'derivatives'.
        """
        curves = {}
        for category, t_crit in self.config.T_crit_dict.items():
            # Apply lower bound constraint
            capped = np.maximum(target_window_beats, t_crit)
            # Round to nearest integer and cast
            curves[category] = np.round(capped).astype(int)
        return curves

    def generate_window_slices(
        self,
        feature_windows_beats: Dict[str, np.ndarray],
        overlap_ratio: np.ndarray
    ) -> Dict[str, List[Tuple[int, int]]]:
        """Converts window size curves and overlap trends into concrete slicing intervals of beat indices.

        Each interval is defined as (start_beat_idx, end_beat_idx) and steps forward by:
        Step = max(1, int(T_w * (1 - O_w)))

        Args:
            feature_windows_beats (Dict[str, np.ndarray]): Dict of integer window curves per category.
            overlap_ratio (np.ndarray): 1D array of overlap percentages.

        Returns:
            Dict[str, List[Tuple[int, int]]]: Slicing schedules grouped by category name.
        """
        schedules = {}
        n_beats = len(overlap_ratio)

        for category, tw_curve in feature_windows_beats.items():
            slices = []
            curr_idx = 0

            while curr_idx < n_beats:
                tw = tw_curve[curr_idx]
                ow = overlap_ratio[curr_idx]

                # Segment starts at current index and spans tw beats
                end_idx = curr_idx + tw
                # Cap the end index to the total number of beats
                if end_idx > n_beats:
                    end_idx = n_beats

                slices.append((curr_idx, end_idx))

                # Step index forward based on adaptive overlap
                step = int(np.round(tw * (1.0 - ow)))
                step = max(1, step)

                # Break if we have reached the end of the signal
                if curr_idx + step >= n_beats or end_idx >= n_beats:
                    break

                curr_idx += step

            schedules[category] = slices

        return schedules

    def process(self, processed_ppg: Any) -> AdaptiveWindowResult:
        """Runs the entire adaptive windowing pipeline on a ProcessedPPG input.

        Args:
            processed_ppg (Any): ProcessedPPG object containing beats and validation mask.

        Returns:
            AdaptiveWindowResult: Computed curves, raw/norm entropy, and slice schedules.

        Raises:
            AdaptiveEngineError: If calculation parameters are invalid.
        """
        # 1. Extract signals and metadata from ProcessedPPG
        if not hasattr(processed_ppg, 'beats_normalized') or not hasattr(processed_ppg, 'valid_beats_mask'):
            raise AdaptiveEngineError(
                "Input must contain 'beats_normalized' matrix and 'valid_beats_mask' array."
            )

        beats_matrix = processed_ppg.beats_normalized
        valid_mask = processed_ppg.valid_beats_mask
        beat_times = processed_ppg.beat_times

        n_beats = len(beats_matrix)

        # 2. Check for empty or very short beat lists
        if n_beats == 0:
            # Return empty arrays safely
            empty_arr = np.array([], dtype=np.float64)
            empty_int = np.array([], dtype=int)
            empty_slices: Dict[str, List[Tuple[int, int]]] = {k: [] for k in self.config.T_crit_dict.keys()}
            empty_windows: Dict[str, np.ndarray] = {k: empty_int for k in self.config.T_crit_dict.keys()}
            return AdaptiveWindowResult(
                beat_times=empty_arr,
                entropy_raw=empty_arr,
                entropy_normalized=empty_arr,
                target_window_beats=empty_arr,
                overlap_ratio=empty_arr,
                feature_windows_beats=empty_windows,
                window_slices=empty_slices
            )

        # Handle edge case where number of beats is less than N_past
        if n_beats < self.config.N_past:
            # Force target window to maximum size, overlap to minimum, and empty entropy
            entropy_raw = np.zeros(n_beats, dtype=np.float64)
            entropy_normalized = np.zeros(n_beats, dtype=np.float64)
            target_window_beats = np.full(n_beats, float(self.config.T_max_beats))
            overlap_ratio = np.full(n_beats, self.config.O_min)
        else:
            # Compute entropy and map values
            entropy_normalized = self.compute_morphological_entropy(beats_matrix, valid_mask)
            # Since we normalized, let's keep the raw array for completeness
            # For the raw, we rebuild the values before normalization
            entropy_raw = np.zeros(n_beats, dtype=np.float64)
            for i in range(n_beats):
                start_window = max(0, i - self.config.N_past + 1)
                win_mask = valid_mask[start_window:i + 1]
                win_beats = beats_matrix[start_window:i + 1][win_mask]
                if len(win_beats) >= 2:
                    diffs = [np.sqrt(np.mean((win_beats[j] - win_beats[j-1])**2)) for j in range(1, len(win_beats))]
                    entropy_raw[i] = np.mean(diffs)
            
            # Match priming phase for raw values too
            first_valid_idx = self.config.N_past - 1
            if first_valid_idx < n_beats:
                entropy_raw[:first_valid_idx] = entropy_raw[first_valid_idx]

            target_window_beats, overlap_ratio = self.map_entropy_to_window(entropy_normalized)

        # 3. Apply constraints
        feature_windows_beats = self.apply_critical_sensitivity(target_window_beats)

        # 4. Generate slices
        window_slices = self.generate_window_slices(feature_windows_beats, overlap_ratio)

        return AdaptiveWindowResult(
            beat_times=beat_times,
            entropy_raw=entropy_raw,
            entropy_normalized=entropy_normalized,
            target_window_beats=target_window_beats,
            overlap_ratio=overlap_ratio,
            feature_windows_beats=feature_windows_beats,
            window_slices=window_slices
        )


if __name__ == '__main__':
    print("--- Running autonomous tests for AdaptiveEngine ---")

    # 1. Generate synthetic beat dataset (100 beats, 128 samples each)
    n_beats = 100
    beat_len = 128
    
    # Base template (sine-like cycle)
    x = np.linspace(0, 2 * np.pi, beat_len)
    base_beat = np.sin(x)

    beats = []
    # Phase 1: Beats 0 to 49 -> morphologically identical
    for i in range(50):
        beats.append(base_beat + 0.01 * np.random.randn(beat_len))

    # Phase 2: Beats 50 to 99 -> morphologically variable (modulated phase and amplitude)
    for i in range(50):
        phase_shift = 0.5 * np.sin(2 * np.pi * i / 10.0)
        amp_mod = 1.0 + 0.3 * np.cos(2 * np.pi * i / 5.0)
        modulated_beat = amp_mod * np.sin(x + phase_shift)
        beats.append(modulated_beat + 0.01 * np.random.randn(beat_len))

    beats_matrix = np.vstack(beats)
    valid_mask = np.ones(n_beats, dtype=bool)
    beat_times = np.arange(n_beats) * 0.8  # Assume 0.8 seconds per beat (75 BPM)

    # 2. Package into mock ProcessedPPG
    mock_ppg = ProcessedPPG(
        raw_signal=np.zeros(1000),
        filtered_signal=np.zeros(1000),
        sampling_rate=100.0,
        peaks=np.array([]),
        onsets=np.array([]),
        beats_raw=[],
        beats_normalized=beats_matrix,
        beat_times=beat_times,
        valid_beats_mask=valid_mask
    )

    # 3. Process with AdaptiveEngine
    config = AdaptiveConfig(
        N_past=10,
        T_min_beats=3,
        T_max_beats=30,
        O_min=0.25,
        O_max=0.85
    )
    engine = AdaptiveEngine(config)
    res = engine.process(mock_ppg)

    # 4. Perform assertions
    # Verify outputs have correct shapes
    assert len(res.entropy_normalized) == n_beats, "Entropy shape mismatch"
    assert len(res.target_window_beats) == n_beats, "Target windows shape mismatch"
    assert len(res.overlap_ratio) == n_beats, "Overlap ratio shape mismatch"
    print("-> Data dimensions: OK")

    # Verify transition in entropy (higher entropy in second half)
    mean_ent_first_half = np.mean(res.entropy_normalized[10:50])
    mean_ent_second_half = np.mean(res.entropy_normalized[50:])
    print(f"-> Mean H_morph: First Half={mean_ent_first_half:.4f}, Second Half={mean_ent_second_half:.4f}")
    assert mean_ent_second_half > mean_ent_first_half, "Entropy should be higher in Phase 2"
    print("-> Morphological Entropy transition validation: SUCCESS")

    # Verify transition in window size (smaller windows in second half due to higher entropy)
    mean_w_first_half = np.mean(res.target_window_beats[10:50])
    mean_w_second_half = np.mean(res.target_window_beats[50:])
    print(f"-> Mean Window Size: First Half={mean_w_first_half:.2f} beats, Second Half={mean_w_second_half:.2f} beats")
    assert mean_w_first_half > mean_w_second_half, "Window size should be smaller in Phase 2"
    print("-> Adaptive Window Size transition validation: SUCCESS")

    # Verify critical sensitivity constraints (T_crit)
    # derivatives category constraint should be >= 10
    assert np.all(res.feature_windows_beats['derivatives'] >= 10), "Derivatives window violates T_crit >= 10 constraint"
    # macro category constraint should be >= 3
    assert np.all(res.feature_windows_beats['macro'] >= 3), "Macro window violates T_crit >= 3 constraint"
    print("-> T_crit sensitivity constraints: VALID")

    # Verify slices boundary safety
    for category, slices in res.window_slices.items():
        assert len(slices) > 0, f"No slices generated for {category}"
        for start, end in slices:
            assert start >= 0 and end <= n_beats, f"Index out of bounds in slice ({start}, {end})"
            assert start < end, f"Invalid slice length ({start}, {end})"
    print("-> Slicing interval index bounds: SAFE")

    # 5. Output synthesis table
    print("\n" + "="*50)
    print("                ADAPTIVE ENGINE SUMMARY")
    print("="*50)
    print(f"{'Metric':<25} | {'Phase 1 (Stab)':<10} | {'Phase 2 (Mod)':<10}")
    print("-"*50)
    print(f"{'Mean H_morph':<25} | {mean_ent_first_half:<10.4f} | {mean_ent_second_half:<10.4f}")
    print(f"{'Mean Target Window (beats)':<25} | {mean_w_first_half:<10.2f} | {mean_w_second_half:<10.2f}")
    print(f"{'Mean Overlap Ratio':<25} | {np.mean(res.overlap_ratio[10:50]):<10.2%} | {np.mean(res.overlap_ratio[50:]):<10.2%}")
    print(f"{'Derivatives Win (beats)':<25} | {np.mean(res.feature_windows_beats['derivatives'][10:50]):<10.1f} | {np.mean(res.feature_windows_beats['derivatives'][50:]):<10.1f}")
    print(f"{'Total Slices (derivatives)':<25} | {len(res.window_slices['derivatives']):<10} | (Total signal)")
    print("="*50)

    print("\n--- All adaptive engine tests passed successfully! ---")
