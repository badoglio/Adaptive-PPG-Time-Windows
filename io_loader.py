"""Biosegnal loading and parsing module.

This module provides classes and structures to load biosegnal data (e.g. PPG, ECG)
from various file formats (CSV, MAT, EDF) into a standardized, strongly-typed
structure `SignalData`.
"""

import io
import os
import pathlib
import tempfile
import warnings
from dataclasses import dataclass, field
from typing import Dict, Union, Optional, Any, List

import numpy as np
import pandas as pd
import scipy.io

# Try to import pyedflib or mne for EDF support
try:
    import pyedflib
except ImportError:
    pyedflib = None

try:
    import mne
except ImportError:
    mne = None


class BioSignalLoaderError(Exception):
    """Custom exception raised for errors occurring during biosegnal loading or parsing."""
    pass


@dataclass
class SignalData:
    """Dataclass holding standard biosegnal data.

    Attributes:
        signal (np.ndarray): 1D array of float64 containing the biosegnal.
        sampling_rate (float): Sampling frequency in Hz.
        time (np.ndarray): 1D array of float64 containing the time vector in seconds.
        channel_name (str): Name of the extracted channel (e.g. "PPG", "PLETH", "ECG", or "Channel_0").
        file_format (str): Identified file format/extension (e.g. "CSV", "MAT", "EDF").
        metadata (Optional[dict]): Additional metadata dictionary from the file header.
    """
    signal: np.ndarray
    sampling_rate: float
    time: np.ndarray = field(default=None)  # type: ignore
    channel_name: str = "Channel_0"
    file_format: str = ""
    metadata: Optional[dict] = field(default_factory=dict)

    def __post_init__(self):
        # 1. Ensure signal is a numpy array
        if not isinstance(self.signal, np.ndarray):
            self.signal = np.array(self.signal, dtype=np.float64)
        else:
            self.signal = self.signal.astype(np.float64)

        # 2. Check signal dimensions
        if self.signal.ndim != 1:
            raise ValueError(f"Signal must be a 1D array, got shape {self.signal.shape}")

        # 3. Handle NaN and Inf values
        non_finite_mask = ~np.isfinite(self.signal)
        if np.any(non_finite_mask):
            warnings.warn("Signal contains NaN or Inf values. Imputing via linear interpolation.")
            x = np.arange(len(self.signal))
            finite_mask = ~non_finite_mask
            if not np.any(finite_mask):
                # All values are non-finite, fall back to zeros
                self.signal = np.zeros_like(self.signal)
            else:
                self.signal[non_finite_mask] = np.interp(
                    x[non_finite_mask],
                    x[finite_mask],
                    self.signal[finite_mask]
                )

        # 4. Handle time vector
        if self.time is None:
            self.time = np.arange(len(self.signal)) / self.sampling_rate
        else:
            if not isinstance(self.time, np.ndarray):
                self.time = np.array(self.time, dtype=np.float64)
            else:
                self.time = self.time.astype(np.float64)

            if self.time.ndim != 1:
                raise ValueError("Time vector must be a 1D array.")
            if len(self.time) != len(self.signal):
                raise ValueError(
                    f"Time vector length ({len(self.time)}) does not match signal length ({len(self.signal)})."
                )

        if self.metadata is None:
            self.metadata = {}


class BioSignalLoader:
    """Loader class containing static methods to parse biosegnals from multiple formats."""

    @staticmethod
    def _is_stream(file_input: Any) -> bool:
        """Helper to check if the input is a stream/buffer."""
        return hasattr(file_input, "read")

    @staticmethod
    def _to_file_path_for_edf(file_input: Any) -> tuple[str, bool]:
        """Ensures a physical file path for EDF reading.

        If file_input is a buffer/stream, writes its contents to a temporary file.

        Returns:
            tuple[str, bool]: (file_path, is_temporary)
        """
        if isinstance(file_input, (str, pathlib.Path)):
            return str(file_input), False

        if BioSignalLoader._is_stream(file_input):
            # Read from buffer and write to temp file
            try:
                if hasattr(file_input, "seek"):
                    try:
                        file_input.seek(0)
                    except Exception:
                        pass
                content = file_input.read()
                if isinstance(content, str):
                    content = content.encode("utf-8")
                
                # Create a temporary file
                temp_file = tempfile.NamedTemporaryFile(delete=False, suffix=".edf")
                temp_file.write(content)
                temp_file.close()
                return temp_file.name, True
            except Exception as e:
                raise BioSignalLoaderError(f"Failed to copy stream content to a temporary EDF file: {e}") from e

        raise BioSignalLoaderError(f"Unsupported input type: {type(file_input)}")

    @classmethod
    def load_csv(
        cls,
        file_input: Union[str, pathlib.Path, io.IOBase, Any],
        signal_column: Union[int, str] = 0,
        sampling_rate: float = 100.0,
        delimiter: str = ','
    ) -> SignalData:
        """Loads a biosegnal from a CSV/TXT file.

        Args:
            file_input (Union[str, pathlib.Path, io.IOBase, Any]): Path to file or file-like buffer.
            signal_column (Union[int, str]): Name or index of the signal column (default 0).
            sampling_rate (float): Sampling frequency in Hz (default 100.0).
            delimiter (str): Column delimiter (default ',').

        Returns:
            SignalData: Loaded and validated signal data.

        Raises:
            BioSignalLoaderError: If parsing or validation fails.
        """
        try:
            if cls._is_stream(file_input):
                if hasattr(file_input, "seek"):
                    try:
                        file_input.seek(0)
                    except Exception:
                        pass
            
            # Read CSV with pandas, let it infer headers initially
            df = pd.read_csv(file_input, sep=delimiter, header='infer')
            
            if df.empty:
                raise BioSignalLoaderError("The CSV file is empty.")

            # Identify if header was correctly inferred or if first row was actually numeric data
            is_numeric_header = False
            if len(df.columns) > 0:
                try:
                    float(df.columns[0])
                    is_numeric_header = True
                except ValueError:
                    pass

            if is_numeric_header:
                # Re-read without header
                if cls._is_stream(file_input):
                    if hasattr(file_input, "seek"):
                        try:
                            file_input.seek(0)
                        except Exception:
                            pass
                df = pd.read_csv(file_input, sep=delimiter, header=None)

            col_list = list(df.columns)

            # Extract signal column
            if isinstance(signal_column, str):
                if signal_column not in df.columns:
                    raise BioSignalLoaderError(
                        f"Column name '{signal_column}' not found. Available: {col_list}"
                    )
                col_name = signal_column
                col_data = df[col_name]
            elif isinstance(signal_column, int):
                if signal_column < 0 or signal_column >= df.shape[1]:
                    raise BioSignalLoaderError(
                        f"Column index {signal_column} out of bounds. CSV has {df.shape[1]} columns."
                    )
                col_name = df.columns[signal_column]
                col_data = df[col_name]
            else:
                raise BioSignalLoaderError(f"signal_column must be int or str, got {type(signal_column)}")

            # Handle duplicated column names if pd.Series is returned as pd.DataFrame
            if isinstance(col_data, pd.DataFrame):
                warnings.warn(f"Duplicate column name '{col_name}' found. Selecting the first one.")
                col_data = col_data.iloc[:, 0]

            signal = col_data.to_numpy()

            # Attempt to extract a time column if present
            time_vector = None
            time_candidates = ["time", "t", "sec", "seconds", "tempo"]
            time_cols = [c for c in df.columns if str(c).strip().lower() in time_candidates]
            
            if time_cols:
                # Make sure the time column is different from the signal column
                time_col = time_cols[0]
                if time_col != col_name:
                    time_data = df[time_col]
                    if isinstance(time_data, pd.DataFrame):
                        time_data = time_data.iloc[:, 0]
                    time_vector = time_data.to_numpy()

            # Build metadata
            metadata = {
                "columns": [str(c) for c in df.columns],
                "shape": df.shape,
                "signal_column_resolved": str(col_name)
            }

            return SignalData(
                signal=signal,
                sampling_rate=sampling_rate,
                time=time_vector,
                channel_name=str(col_name),
                file_format="CSV",
                metadata=metadata
            )
        except BioSignalLoaderError:
            raise
        except Exception as e:
            raise BioSignalLoaderError(f"Error loading CSV file: {e}") from e

    @classmethod
    def load_mat(
        cls,
        file_input: Union[str, pathlib.Path, io.BytesIO, Any],
        key: Optional[str] = None,
        sampling_rate: float = 100.0,
        channel_index: int = 0
    ) -> SignalData:
        """Loads a biosegnal from a MATLAB (.mat) file.

        Args:
            file_input (Union[str, pathlib.Path, io.BytesIO, Any]): Path to file or in-memory BytesIO buffer.
            key (Optional[str]): The dictionary key of the signal in the MAT file.
                                 If None, automatically picks the first numeric array key.
            sampling_rate (float): Sampling frequency in Hz (default 100.0).
            channel_index (int): Index of the channel to extract if the key points to a 2D array.

        Returns:
            SignalData: Standardized biosegnal data.

        Raises:
            BioSignalLoaderError: If loading or parsing fails.
        """
        try:
            if cls._is_stream(file_input):
                if hasattr(file_input, "seek"):
                    try:
                        file_input.seek(0)
                    except Exception:
                        pass

            try:
                mat_dict = scipy.io.loadmat(file_input)
            except Exception as e:
                raise BioSignalLoaderError(f"scipy.io.loadmat failed to read MAT file: {e}") from e

            # Filter out MATLAB system headers/keys
            user_keys = {k: v for k, v in mat_dict.items() if not k.startswith("__")}

            if not user_keys:
                raise BioSignalLoaderError("No non-system variables found in MAT file.")

            selected_key = key
            if selected_key is None:
                # Auto-select the first key containing a numeric array
                for k, v in user_keys.items():
                    if isinstance(v, np.ndarray) and np.issubdtype(v.dtype, np.number):
                        selected_key = k
                        break
                if selected_key is None:
                    raise BioSignalLoaderError("Could not find any numeric numpy arrays in MAT file.")
            else:
                if selected_key not in user_keys:
                    raise BioSignalLoaderError(
                        f"Specified key '{selected_key}' not found. Available keys: {list(user_keys.keys())}"
                    )

            data_array = user_keys[selected_key]
            if not isinstance(data_array, np.ndarray):
                raise BioSignalLoaderError(f"Data for key '{selected_key}' is not a NumPy array.")

            # If array has more than 1 dimension, down-mix or isolate the channel
            if data_array.ndim > 1:
                shape = data_array.shape
                if len(shape) == 2:
                    # Determine whether channels are rows or columns
                    if shape[0] <= shape[1]:
                        # Row-wise channels (e.g. shape 2x1000)
                        if channel_index < 0 or channel_index >= shape[0]:
                            raise BioSignalLoaderError(
                                f"Channel index {channel_index} out of bounds for row-wise MAT array of shape {shape}."
                            )
                        signal = data_array[channel_index, :]
                    else:
                        # Column-wise channels (e.g. shape 1000x2)
                        if channel_index < 0 or channel_index >= shape[1]:
                            raise BioSignalLoaderError(
                                f"Channel index {channel_index} out of bounds for column-wise MAT array of shape {shape}."
                            )
                        signal = data_array[:, channel_index]
                else:
                    # 3D or higher, flatten extra dimensions
                    warnings.warn(f"MAT array for key '{selected_key}' has shape {shape}. Flattening extra dimensions.")
                    flat_data = data_array.squeeze()
                    if flat_data.ndim > 1:
                        # Select first slice along extra dimensions
                        idx = (0,) * (flat_data.ndim - 1) + (slice(None),)
                        signal = flat_data[idx]
                    else:
                        signal = flat_data
            else:
                signal = data_array.squeeze()

            # Attempt to extract time vector
            time_vector = None
            time_candidates = ["time", "t", "sec", "seconds", "tempo"]
            time_keys = [k for k in user_keys.keys() if k.lower() in time_candidates]
            
            if time_keys:
                t_arr = user_keys[time_keys[0]]
                if isinstance(t_arr, np.ndarray):
                    time_vector = t_arr.squeeze()

            metadata = {
                "all_keys": list(user_keys.keys()),
                "selected_key": selected_key,
                "original_shape": data_array.shape,
                "channel_index": channel_index
            }

            return SignalData(
                signal=signal,
                sampling_rate=sampling_rate,
                time=time_vector,
                channel_name=f"{selected_key}_ch{channel_index}" if data_array.ndim > 1 else selected_key,
                file_format="MAT",
                metadata=metadata
            )
        except BioSignalLoaderError:
            raise
        except Exception as e:
            raise BioSignalLoaderError(f"Error loading MAT file: {e}") from e

    @classmethod
    def load_edf(cls, file_input: Union[str, pathlib.Path, io.BytesIO, Any], channel_index: int = 0) -> SignalData:
        """Loads a biosegnal from an EDF/BDF file.

        Args:
            file_input (Union[str, pathlib.Path, io.BytesIO, Any]): Path to file or in-memory BytesIO buffer.
            channel_index (int): Index of the channel to extract (default 0).

        Returns:
            SignalData: Standardized biosegnal data.

        Raises:
            BioSignalLoaderError: If loading or parsing fails.
        """
        file_path = ""
        is_temp = False
        try:
            file_path, is_temp = cls._to_file_path_for_edf(file_input)

            if pyedflib is not None:
                try:
                    reader = pyedflib.EdfReader(file_path)
                except Exception as e:
                    raise BioSignalLoaderError(f"pyedflib failed to open EDF file: {e}") from e

                try:
                    num_channels = reader.signals_in_file
                    if channel_index < 0 or channel_index >= num_channels:
                        raise BioSignalLoaderError(
                            f"Channel index {channel_index} out of bounds. File has {num_channels} channels."
                        )

                    channel_name = reader.getSignalLabels()[channel_index]
                    sampling_rate = reader.getSampleFrequency(channel_index)
                    signal = reader.readSignal(channel_index)

                    metadata = {
                        "patient_name": reader.getPatientName(),
                        "patient_code": reader.getPatientCode(),
                        "gender": reader.getSex() if hasattr(reader, "getSex") else reader.getGender(),
                        "birthdate": reader.getBirthdate(),
                        "admincode": reader.getAdmincode(),
                        "technician": reader.getTechnician(),
                        "equipment": reader.getEquipment(),
                        "recording_additional": reader.getRecordingAdditional(),
                        "startdate": reader.getStartdatetime().isoformat() if reader.getStartdatetime() else None,
                        "signal_headers": reader.getSignalHeaders()[channel_index]
                    }
                finally:
                    reader.close()
            elif mne is not None:
                # Fallback implementation using mne
                try:
                    raw = mne.io.read_raw_edf(file_path, preload=True, verbose=False)
                    ch_names = raw.info['ch_names']
                    if channel_index < 0 or channel_index >= len(ch_names):
                        raise BioSignalLoaderError(
                            f"Channel index {channel_index} out of bounds. File has {len(ch_names)} channels."
                        )
                    channel_name = ch_names[channel_index]
                    sampling_rate = raw.info['sfreq']
                    signal = raw.get_data(picks=[channel_name])[0]
                    metadata = {
                        "ch_names": ch_names,
                        "sfreq": sampling_rate,
                        "meas_date": str(raw.info['meas_date']) if raw.info['meas_date'] else None
                    }
                except Exception as e:
                    raise BioSignalLoaderError(f"MNE fallback failed to read EDF: {e}") from e
            else:
                raise BioSignalLoaderError(
                    "EDF parsing libraries not available. Install pyedflib or mne."
                )

            return SignalData(
                signal=signal,
                sampling_rate=float(sampling_rate),
                channel_name=channel_name,
                file_format="EDF",
                metadata=metadata
            )
        except BioSignalLoaderError:
            raise
        except Exception as e:
            raise BioSignalLoaderError(f"Error loading EDF file: {e}") from e
        finally:
            if is_temp and file_path and os.path.exists(file_path):
                try:
                    os.unlink(file_path)
                except Exception as e:
                    warnings.warn(f"Failed to delete temporary file '{file_path}': {e}")

    @classmethod
    def auto_load(
        cls,
        file_input: Union[str, pathlib.Path, io.IOBase, Any],
        file_name: Optional[str] = None,
        **kwargs
    ) -> SignalData:
        """Dispatcher method that identifies file format and loads data.

        Args:
            file_input (Union[str, pathlib.Path, io.IOBase, Any]): File path or in-memory buffer.
            file_name (Optional[str]): Explicit file name (useful if file_input is a buffer).
            **kwargs: Arguments passed to specific load methods.

        Returns:
            SignalData: Standardized biosegnal data.

        Raises:
            BioSignalLoaderError: If format cannot be identified or loading fails.
        """
        resolved_name = ""
        if file_name is not None:
            resolved_name = file_name
        elif isinstance(file_input, (str, pathlib.Path)):
            resolved_name = str(file_input)
        elif hasattr(file_input, 'name') and isinstance(file_input.name, str):
            resolved_name = file_input.name

        if not resolved_name:
            raise BioSignalLoaderError(
                "Cannot identify file format because file name/path was not provided "
                "and the buffer has no 'name' attribute. Please provide 'file_name'."
            )

        ext = pathlib.Path(resolved_name).suffix.lower()

        if ext in (".csv", ".txt", ".tsv"):
            return cls.load_csv(file_input, **kwargs)
        elif ext == ".mat":
            return cls.load_mat(file_input, **kwargs)
        elif ext in (".edf", ".bdf"):
            return cls.load_edf(file_input, **kwargs)
        else:
            raise BioSignalLoaderError(f"Unsupported file extension '{ext}' in file '{resolved_name}'")


if __name__ == '__main__':
    # 1. Generate synthetic data (100 Hz, 10 seconds sine wave)
    fs = 100.0
    duration = 10.0
    t = np.arange(int(fs * duration)) / fs
    freq = 1.0  # 1 Hz
    signal = np.sin(2 * np.pi * freq * t)
    
    print("--- Running autonomous tests for BioSignalLoader ---")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = pathlib.Path(tmpdir)
        csv_path = tmp_path / "test_signal.csv"
        mat_path = tmp_path / "test_signal.mat"
        
        # Save CSV with time and PPG columns
        df_synthetic = pd.DataFrame({"time": t, "PPG": signal})
        df_synthetic.to_csv(csv_path, index=False)
        print(f"Created temporary CSV file at: {csv_path}")
        
        # Save MAT file with PPG and time
        scipy.io.savemat(str(mat_path), {"PPG": signal, "time": t})
        print(f"Created temporary MAT file at: {mat_path}")
        
        # 2. Test CSV Loading
        print("\nTesting CSV Loading...")
        # Load from path
        data_csv_path = BioSignalLoader.load_csv(csv_path, signal_column="PPG", sampling_rate=fs)
        assert np.allclose(data_csv_path.signal, signal), "CSV signal mismatch (path)"
        assert np.allclose(data_csv_path.time, t), "CSV time mismatch (path)"
        assert data_csv_path.sampling_rate == fs, "CSV sampling rate mismatch (path)"
        assert data_csv_path.channel_name == "PPG", "CSV channel name mismatch (path)"
        print("-> Path loading: SUCCESS")
        
        # Load from StringIO (buffer)
        with open(csv_path, "r") as f:
            csv_buffer = io.StringIO(f.read())
        data_csv_buf = BioSignalLoader.load_csv(csv_buffer, signal_column=1, sampling_rate=fs)
        assert np.allclose(data_csv_buf.signal, signal), "CSV signal mismatch (buffer)"
        assert np.allclose(data_csv_buf.time, t), "CSV time mismatch (buffer)"
        print("-> StringIO buffer loading: SUCCESS")
        
        # 3. Test MAT Loading
        print("\nTesting MAT Loading...")
        # Load from path (auto-selecting key)
        data_mat_path = BioSignalLoader.load_mat(mat_path, key=None, sampling_rate=fs)
        # Note: auto-select will pick 'PPG' or 'time'. Let's specify key to ensure PPG is picked
        data_mat_path_explicit = BioSignalLoader.load_mat(mat_path, key="PPG", sampling_rate=fs)
        assert np.allclose(data_mat_path_explicit.signal, signal), "MAT signal mismatch (path)"
        assert np.allclose(data_mat_path_explicit.time, t), "MAT time mismatch (path)"
        print("-> Path loading: SUCCESS")
        
        # Load from BytesIO (buffer)
        with open(mat_path, "rb") as f:
            mat_buffer = io.BytesIO(f.read())
        data_mat_buf = BioSignalLoader.load_mat(mat_buffer, key="PPG", sampling_rate=fs)
        assert np.allclose(data_mat_buf.signal, signal), "MAT signal mismatch (buffer)"
        assert np.allclose(data_mat_buf.time, t), "MAT time mismatch (buffer)"
        print("-> BytesIO buffer loading: SUCCESS")
        
        # 4. Test EDF Loading (if pyedflib is available)
        if pyedflib is not None:
            print("\nTesting EDF Loading (pyedflib is available)...")
            edf_path = tmp_path / "test_signal.edf"
            # Write synthetic EDF
            writer = None
            try:
                writer = pyedflib.EdfWriter(str(edf_path), 1, file_type=pyedflib.FILETYPE_EDFPLUS)
                channel_info = {
                    'label': 'PPG',
                    'dimension': 'uV',
                    'sample_frequency': int(fs),
                    'physical_max': 10.0,
                    'physical_min': -10.0,
                    'digital_max': 32767,
                    'digital_min': -32768,
                    'prefilter': '',
                    'transducer': ''
                }
                writer.setSignalHeader(0, channel_info)
                writer.writeSamples([signal])
            finally:
                if writer is not None:
                    writer.close()
            print(f"Created temporary EDF file at: {edf_path}")
            
            # Load from path
            data_edf_path = BioSignalLoader.load_edf(edf_path, channel_index=0)
            # Digital quantization causes minor discrepancies, use atol=1e-3
            assert np.allclose(data_edf_path.signal, signal, atol=1e-3), "EDF signal mismatch (path)"
            assert np.allclose(data_edf_path.time, t), "EDF time mismatch (path)"
            assert data_edf_path.sampling_rate == fs, "EDF sampling rate mismatch"
            assert data_edf_path.channel_name == "PPG", "EDF channel label mismatch"
            print("-> Path loading: SUCCESS")
            
            # Load from BytesIO (buffer)
            with open(edf_path, "rb") as f:
                edf_buffer = io.BytesIO(f.read())
            data_edf_buf = BioSignalLoader.load_edf(edf_buffer, channel_index=0)
            assert np.allclose(data_edf_buf.signal, signal, atol=1e-3), "EDF signal mismatch (buffer)"
            assert np.allclose(data_edf_buf.time, t), "EDF time mismatch (buffer)"
            print("-> BytesIO buffer loading: SUCCESS")
        else:
            print("\nSkipping EDF tests (pyedflib not available).")
            
        # 5. Test Auto-load dispatcher
        print("\nTesting Auto-load dispatcher...")
        data_auto_csv = BioSignalLoader.auto_load(csv_path, signal_column="PPG", sampling_rate=fs)
        assert np.allclose(data_auto_csv.signal, signal), "Auto-load CSV mismatch"
        
        data_auto_mat = BioSignalLoader.auto_load(mat_path, key="PPG", sampling_rate=fs)
        assert np.allclose(data_auto_mat.signal, signal), "Auto-load MAT mismatch"
        print("-> Auto-load routing: SUCCESS")
        
        # 6. Test Edge Cases: NaN / Inf handling
        print("\nTesting NaN/Inf imputation...")
        signal_with_nan = signal.copy()
        signal_with_nan[50] = np.nan
        signal_with_nan[150] = np.inf
        # Instantiating SignalData with bad values should trigger imputation and a warning
        import warnings
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            data_nan = SignalData(signal=signal_with_nan, sampling_rate=fs)
            assert len(w) >= 1, "NaN/Inf warning was not raised"
            assert "NaN or Inf" in str(w[-1].message), "Incorrect warning message"
            
        assert np.isfinite(data_nan.signal).all(), "Imputed signal still contains non-finite values"
        # Verify that imputed value is close to the original (linear interpolation between index 49 and 51)
        expected_50 = (signal[49] + signal[51]) / 2.0
        assert np.allclose(data_nan.signal[50], expected_50, atol=1e-3), f"Imputed value mismatch at 50: {data_nan.signal[50]} vs {expected_50}"
        print("-> NaN/Inf Imputation: SUCCESS")
        
        # 7. Test Edge Cases: Multi-channel down-mixing for MAT
        print("\nTesting multi-channel down-mixing for MAT...")
        multi_signal = np.vstack([signal, -signal])  # 2 x 1000 array
        multi_mat_path = tmp_path / "multi_channel.mat"
        scipy.io.savemat(str(multi_mat_path), {"data": multi_signal})
        
        # Load first channel
        data_ch0 = BioSignalLoader.load_mat(multi_mat_path, key="data", channel_index=0)
        assert np.allclose(data_ch0.signal, signal), "Multi-channel MAT load ch0 mismatch"
        # Load second channel
        data_ch1 = BioSignalLoader.load_mat(multi_mat_path, key="data", channel_index=1)
        assert np.allclose(data_ch1.signal, -signal), "Multi-channel MAT load ch1 mismatch"
        print("-> Multi-channel down-mixing: SUCCESS")
        
        # 8. Test Error raising
        print("\nTesting error cases...")
        try:
            BioSignalLoader.load_csv(csv_path, signal_column="NON_EXISTENT")
            raise AssertionError("Should have raised BioSignalLoaderError for invalid column")
        except BioSignalLoaderError:
            pass
            
        try:
            BioSignalLoader.auto_load(csv_path, file_name="test.unsupported")
            raise AssertionError("Should have raised BioSignalLoaderError for unsupported format")
        except BioSignalLoaderError:
            pass
        print("-> Error handling: SUCCESS")
        
    print("\n--- All tests completed successfully! ---")
    
    # Print metrics summary as requested
    print(f"Signal Length: {len(signal)} samples")
    print(f"Sampling Frequency: {fs} Hz")
    print(f"Duration: {duration} seconds")
    print(f"Verification matched synth data to within 1e-3 tolerance.")
