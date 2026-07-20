"""PPG signal preprocessing and beat segmentation module.

This module provides the `PPGPreprocessor` class and `ProcessedPPG` dataclass
for filtering raw PPG signals, detecting systolic peaks and onsets, segmenting
individual heartbeats, performing amplitude and temporal normalization,
and validating beats based on physiological criteria.
"""

import warnings
from dataclasses import dataclass, field
from typing import List, Tuple, Union, Any, Optional

import numpy as np
import scipy.signal
import scipy.interpolate
import biosppy.signals.ppg as ppg

# Safe import of SignalData from io_loader.py to allow standalone testing
try:
    from io_loader import SignalData
except ImportError:
    # Standalone mock implementation of SignalData
    @dataclass
    class SignalData:
        """Mock SignalData class for standalone testing."""
        signal: np.ndarray
        sampling_rate: float
        time: np.ndarray = field(default=None)  # type: ignore
        channel_name: str = "Channel_0"
        file_format: str = ""
        metadata: Optional[dict] = field(default_factory=dict)

        def __post_init__(self):
            if not isinstance(self.signal, np.ndarray):
                self.signal = np.array(self.signal, dtype=np.float64)
            else:
                self.signal = self.signal.astype(np.float64)
            if self.signal.ndim != 1:
                raise ValueError("Signal must be 1D.")
            if self.time is None:
                self.time = np.arange(len(self.signal)) / self.sampling_rate
            if self.metadata is None:
                self.metadata = {}


class PPGPreprocessingError(Exception):
    """Custom exception raised for errors during PPG preprocessing."""
    pass


@dataclass
class ProcessedPPG:
    """Dataclass holding processed PPG signal and segmented beat information.

    Attributes:
        raw_signal (np.ndarray): 1D array of the original raw PPG signal.
        filtered_signal (np.ndarray): 1D array of the bandpass filtered PPG signal.
        sampling_rate (float): Sampling frequency in Hz.
        peaks (np.ndarray): 1D array of integers containing systolic peak indices.
        onsets (np.ndarray): 1D array of integers containing beat onset (valley) indices.
        beats_raw (List[np.ndarray]): List of 1D arrays representing raw variable-length beats (onset to onset).
        beats_normalized (np.ndarray): 2D array of shape [n_beats, target_length] containing normalized beats.
        beat_times (np.ndarray): 1D array of floats containing start times (seconds) of each beat.
        valid_beats_mask (np.ndarray): 1D boolean array indicating valid beats.
    """
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
    """Preprocessor class that orchestrates the PPG filtering, detection, and segmentation pipeline."""

    def __init__(
        self,
        lowcut: float = 0.5,
        highcut: float = 30.0,
        filter_order: int = 4,
        target_beat_length: int = 128,
        bypass_filter: bool = False
    ):
        """Initializes the preprocessor with filter parameters.

        Args:
            lowcut (float): High-pass cutoff frequency in Hz. Default 0.5.
            highcut (float): Low-pass cutoff frequency in Hz. Default 30.0.
            filter_order (int): Order of the Butterworth filter. Default 4.
            target_beat_length (int): Fixed number of samples for temporal normalization. Default 128.
            bypass_filter (bool): If True, skips the filtering phase during processing. Default False.
        """
        self.lowcut = lowcut
        self.highcut = highcut
        self.filter_order = filter_order
        self.target_beat_length = target_beat_length
        self.bypass_filter = bypass_filter

    def filter_signal(self, signal: np.ndarray, sampling_rate: float) -> np.ndarray:
        """Applies a zero-phase Butterworth bandpass filter after detrending.

        Args:
            signal (np.ndarray): 1D raw PPG signal.
            sampling_rate (float): Sampling rate of the signal in Hz.

        Returns:
            np.ndarray: Filtered signal of the same shape as input.

        Raises:
            PPGPreprocessingError: If the signal is too short to apply the filter.
        """
        if len(signal) < 3 * self.filter_order:
            raise PPGPreprocessingError(
                f"Signal length ({len(signal)}) is too short for filter order {self.filter_order}. "
                f"Minimum length is {3 * self.filter_order}."
            )

        # 1. Detrend signal to avoid edge transients
        detrended = scipy.signal.detrend(signal)

        # 2. Design Butterworth bandpass filter
        nyq = 0.5 * sampling_rate
        low = self.lowcut / nyq
        high = self.highcut / nyq

        # Validate normalized frequency ranges
        if low <= 0 or low >= 1:
            raise PPGPreprocessingError(
                f"Normalized lowcut frequency {low:.4f} is invalid. Check lowcut ({self.lowcut} Hz) "
                f"relative to Nyquist frequency ({nyq} Hz)."
            )
        if high <= 0 or high >= 1:
            raise PPGPreprocessingError(
                f"Normalized highcut frequency {high:.4f} is invalid. Check highcut ({self.highcut} Hz) "
                f"relative to Nyquist frequency ({nyq} Hz)."
            )
        if low >= high:
            raise PPGPreprocessingError(
                f"Lowcut frequency ({self.lowcut} Hz) must be smaller than highcut frequency ({self.highcut} Hz)."
            )

        b, a = scipy.signal.butter(self.filter_order, [low, high], btype='band')

        # 3. Apply zero-phase filter
        try:
            filtered = scipy.signal.filtfilt(b, a, detrended)
        except Exception as e:
            raise PPGPreprocessingError(f"Filtering operation failed: {e}") from e

        return filtered

    def detect_beats(self, filtered_signal: np.ndarray, sampling_rate: float) -> Tuple[np.ndarray, np.ndarray]:
        """Detects pulse peaks and onsets using BioSPPy detectors and local minima checks.

        Args:
            filtered_signal (np.ndarray): 1D filtered PPG signal.
            sampling_rate (float): Sampling frequency in Hz.

        Returns:
            Tuple[np.ndarray, np.ndarray]: (peaks, onsets) as 1D arrays of integer indices.
        """
        if len(filtered_signal) == 0:
            return np.array([], dtype=int), np.array([], dtype=int)

        # Check for flat signal to avoid numeric noise causing false detections
        sig_range = np.max(filtered_signal) - np.min(filtered_signal)
        if sig_range < 1e-5:
            return np.array([], dtype=int), np.array([], dtype=int)

        try:
            # 1. Find initial peaks using Elgendi's algorithm (stored in the 'onsets' field of the return tuple)
            res = ppg.find_onsets_elgendi2013(signal=filtered_signal, sampling_rate=sampling_rate)
            peaks = res['onsets']

            if len(peaks) == 0:
                return np.array([], dtype=int), np.array([], dtype=int)

            # 2. Local minima extraction for robust onset matching (avoiding crashing on edge cases)
            minima = (np.diff(np.sign(np.diff(filtered_signal))) > 0).nonzero()[0]

            onsets_list = []
            peaks_list = []
            for peak_idx in peaks:
                prev_minima = minima[minima < peak_idx]
                if len(prev_minima) > 0:
                    onsets_list.append(prev_minima[-1])
                    peaks_list.append(peak_idx)

            if len(onsets_list) == 0:
                return np.array([], dtype=int), np.array([], dtype=int)

            return np.array(peaks_list, dtype=int), np.array(onsets_list, dtype=int)

        except Exception as e:
            warnings.warn(f"BioSPPy peak/onset detection failed or raised an exception: {e}")
            return np.array([], dtype=int), np.array([], dtype=int)

    def segment_and_normalize_beats(
        self,
        filtered_signal: np.ndarray,
        onsets: np.ndarray,
        peaks: np.ndarray
    ) -> Tuple[List[np.ndarray], np.ndarray, np.ndarray]:
        """Segments the signal from onset to onset, and normalizes beat amplitude and duration.

        Args:
            filtered_signal (np.ndarray): 1D array of filtered PPG.
            onsets (np.ndarray): 1D array of onset indices.
            peaks (np.ndarray): 1D array of peak indices.

        Returns:
            Tuple[List[np.ndarray], np.ndarray, np.ndarray]:
                - List of 1D arrays representing raw variable-length beats.
                - 2D array of shape [n_beats, target_beat_length] representing normalized beats.
                - 1D array of beat start times in seconds.
        """
        beats_raw: List[np.ndarray] = []
        beats_norm_list: List[np.ndarray] = []
        beat_times_list: List[float] = []

        n_beats = len(onsets) - 1
        sampling_rate = 100.0  # Safe default if not computed, overridden in process

        for i in range(n_beats):
            start_idx = int(onsets[i])
            end_idx = int(onsets[i+1])
            beat = filtered_signal[start_idx:end_idx]

            beats_raw.append(beat)

            # Amplitude Normalization (Min-Max Scaling to [0, 1])
            b_min = beat.min()
            b_max = beat.max()
            if b_max > b_min:
                beat_scaled = (beat - b_min) / (b_max - b_min)
            else:
                beat_scaled = np.zeros_like(beat)

            # Temporal Normalization (Resampling to target_beat_length points)
            if len(beat_scaled) >= 2:
                x_old = np.linspace(0, 1, len(beat_scaled))
                x_new = np.linspace(0, 1, self.target_beat_length)
                f_interp = scipy.interpolate.interp1d(x_old, beat_scaled, kind='linear')
                beat_norm = f_interp(x_new)
            else:
                beat_norm = np.zeros(self.target_beat_length)

            beats_norm_list.append(beat_norm)

        # Build final structures
        if n_beats > 0:
            beats_normalized = np.vstack(beats_norm_list)
        else:
            beats_normalized = np.empty((0, self.target_beat_length))

        return beats_raw, beats_normalized, onsets[:-1]

    def validate_beats(self, beats_raw: List[np.ndarray], sampling_rate: float) -> np.ndarray:
        """Validates beats based on physiological BPM ranges (30-220 BPM) and signal quality (NaN/Inf check).

        Args:
            beats_raw (List[np.ndarray]): List of 1D arrays representing raw beats.
            sampling_rate (float): Sampling frequency in Hz.

        Returns:
            np.ndarray: 1D boolean array indicating valid (True) or invalid (False) beats.
        """
        valid_mask = []
        for beat in beats_raw:
            if len(beat) == 0:
                valid_mask.append(False)
                continue

            # 1. Check for non-finite values (NaN / Inf)
            if not np.isfinite(beat).all():
                valid_mask.append(False)
                continue

            # 2. Check physiological duration constraints
            duration_sec = len(beat) / sampling_rate
            bpm = 60.0 / duration_sec
            if bpm < 30.0 or bpm > 220.0:
                valid_mask.append(False)
            else:
                valid_mask.append(True)

        return np.array(valid_mask, dtype=bool)

    def process(self, signal_data: Any, sampling_rate: Optional[float] = None, bypass_filter: Optional[bool] = None) -> ProcessedPPG:
        """Executes the full preprocessing pipeline on the input signal.

        Args:
            signal_data (Any): SignalData object or 1D array of raw signal values.
            sampling_rate (Optional[float]): Sampling rate in Hz. Required if signal_data is an array.
            bypass_filter (Optional[bool]): If provided, overrides self.bypass_filter for this run.

        Returns:
            ProcessedPPG: A container for the processed signals and extracted/normalized beats.

        Raises:
            PPGPreprocessingError: If input signals are invalid or processing fails.
        """
        # 1. Resolve raw signal and sampling frequency
        if hasattr(signal_data, 'signal') and hasattr(signal_data, 'sampling_rate'):
            raw_signal = signal_data.signal
            fs = float(signal_data.sampling_rate)
        else:
            raw_signal = np.asarray(signal_data)
            if sampling_rate is None:
                raise PPGPreprocessingError(
                    "A sampling_rate must be supplied when passing a raw signal array."
                )
            fs = float(sampling_rate)

        if raw_signal.ndim > 1:
            warnings.warn("Raw signal has more than 1 dimension. Extracting the first channel.")
            raw_signal = raw_signal.squeeze()
            if raw_signal.ndim > 1:
                # Down-mix multi-dimensional inputs by taking the first column/row
                if raw_signal.shape[0] <= raw_signal.shape[1]:
                    raw_signal = raw_signal[0, :]
                else:
                    raw_signal = raw_signal[:, 0]

        # Convert to float64 for computational stability
        raw_signal = raw_signal.astype(np.float64)

        # Determine filter bypass setting
        should_bypass = self.bypass_filter if bypass_filter is None else bypass_filter

        # 2. Apply filtering
        if should_bypass:
            filtered_signal = raw_signal.copy()
        else:
            filtered_signal = self.filter_signal(raw_signal, fs)

        # 3. Detect peaks and onsets
        peaks, onsets = self.detect_beats(filtered_signal, fs)

        # 4. Handle case with no detected beats
        if len(onsets) < 2:
            warnings.warn("No beats (or fewer than 2 onsets) were detected in the signal.")
            empty_beats: List[np.ndarray] = []
            empty_norm = np.empty((0, self.target_beat_length))
            empty_times = np.array([], dtype=float)
            empty_mask = np.array([], dtype=bool)
            return ProcessedPPG(
                raw_signal=raw_signal,
                filtered_signal=filtered_signal,
                sampling_rate=fs,
                peaks=peaks,
                onsets=onsets,
                beats_raw=empty_beats,
                beats_normalized=empty_norm,
                beat_times=empty_times,
                valid_beats_mask=empty_mask
            )

        # 5. Segment and normalize
        beats_raw, beats_normalized, onset_indices = self.segment_and_normalize_beats(
            filtered_signal, onsets, peaks
        )
        beat_times = onset_indices / fs

        # 6. Validate segments
        valid_beats_mask = self.validate_beats(beats_raw, fs)

        return ProcessedPPG(
            raw_signal=raw_signal,
            filtered_signal=filtered_signal,
            sampling_rate=fs,
            peaks=peaks,
            onsets=onsets,
            beats_raw=beats_raw,
            beats_normalized=beats_normalized,
            beat_times=beat_times,
            valid_beats_mask=valid_beats_mask
        )


if __name__ == '__main__':
    print("--- Running autonomous tests for PPGPreprocessor ---")

    # 1. Generate a realistic synthetic PPG signal
    # Duration: 20 seconds, sampling rate: 100 Hz, HR: 75 BPM (period 0.8s)
    fs = 100.0
    duration = 20.0
    t = np.arange(int(duration * fs)) / fs
    hr = 75.0
    period = 60.0 / hr  # 0.8 seconds

    # Generate pulse template based on Gaussian curves
    phase = (t % period) / period
    systolic_wave = np.exp(-((phase - 0.15) / 0.07)**2)
    diastolic_wave = 0.35 * np.exp(-((phase - 0.45) / 0.10)**2)
    clean_ppg = systolic_wave + diastolic_wave

    # Add realistic noise and baseline wander
    baseline_wander = 0.25 * np.sin(2 * np.pi * 0.15 * t)  # slow respiratory drift
    noise = 0.05 * np.random.randn(len(t))
    synthetic_raw_signal = clean_ppg + baseline_wander + noise

    # 2. Build mock/real SignalData
    signal_data = SignalData(
        signal=synthetic_raw_signal,
        sampling_rate=fs,
        channel_name="PPG_Synth",
        file_format="SYNTHETIC"
    )

    # 3. Instantiate Preprocessor and process
    preprocessor = PPGPreprocessor(lowcut=0.5, highcut=30.0, target_beat_length=128)
    processed = preprocessor.process(signal_data)

    # 4. Perform assertions
    # Verify filtered signal length
    assert len(processed.filtered_signal) == len(synthetic_raw_signal), "Length mismatch in filtered signal"
    print("-> Filtered signal length: MATCHED")

    # Verify beat count and shapes
    n_beats = len(processed.beats_raw)
    print(f"-> Detected beats: {n_beats}")
    assert n_beats > 0, "No beats were detected from the synthetic signal"
    assert processed.beats_normalized.shape == (n_beats, 128), f"Normalized shape mismatch: {processed.beats_normalized.shape}"
    print(f"-> Beats normalized matrix shape: {processed.beats_normalized.shape} (MATCHED)")

    # Verify normalization range
    assert np.nanmin(processed.beats_normalized) >= 0.0, "Normalized values contain items below 0.0"
    assert np.nanmax(processed.beats_normalized) <= 1.0, "Normalized values contain items above 1.0"
    print("-> Amplitude normalization limits: VALID [0.0, 1.0]")

    # Verify time vector and times array size
    assert len(processed.beat_times) == n_beats, "Beat times list size mismatch with beats count"
    assert len(processed.valid_beats_mask) == n_beats, "Validation mask size mismatch with beats count"
    print("-> Mask and beat times alignment: OK")

    # 5. Test filter bypass ("take-over / bypass")
    bypassed_processed = preprocessor.process(signal_data, bypass_filter=True)
    assert np.allclose(bypassed_processed.filtered_signal, synthetic_raw_signal), "Bypass filter failed to retain raw signal"
    print("-> Filter bypass option: SUCCESS")

    # 6. Test short signal error handling
    short_signal = np.sin(np.linspace(0, 1, 10))
    try:
        preprocessor.process(short_signal, sampling_rate=fs)
        raise AssertionError("Should have raised PPGPreprocessingError for too short signal")
    except PPGPreprocessingError:
        print("-> Short signal exception handling: SUCCESS")

    # 7. Test empty detection handling (pass flat signal)
    flat_signal = np.ones(1000)
    flat_processed = preprocessor.process(flat_signal, sampling_rate=fs)
    assert len(flat_processed.beats_raw) == 0, "Flat signal should result in zero beats"
    assert flat_processed.beats_normalized.shape == (0, 128), "Normalized shape must have 0 rows for zero beats"
    print("-> Zero beats warning & empty return handling: SUCCESS")

    print("\n--- All preprocessing tests passed successfully! ---")
    print(f"Extracted beats: {n_beats}")
    print(f"Normalized Matrix: {processed.beats_normalized.shape}")
