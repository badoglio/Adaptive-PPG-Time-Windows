"""Biosignal loading and parsing.

Loads PPG (and optionally ECG) channels from CSV/TXT/TSV, MATLAB and EDF/BDF
files into a :class:`SignalData` container.

Sampling frequency handling: when a time column is present the sampling rate
is estimated from the median sample spacing. If the caller also passes an
explicit ``sampling_rate`` that disagrees by more than 1 %, the explicit value
is kept but a warning is logged and both values are stored in ``metadata``.
Passing ``sampling_rate=None`` means "infer it from the file".
"""

from __future__ import annotations

import io
import logging
import os
import pathlib
import re
import tempfile
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
import scipy.io

from ._typing import FloatArray

logger = logging.getLogger(__name__)

FileInput = Union[str, pathlib.Path, io.IOBase, Any]

_TIME_COLUMN_RE = re.compile(r"^\s*(time|t|sec|secs|seconds|tempo|timestamp)(\b|_)", re.IGNORECASE)
_MS_UNIT_RE = re.compile(r"(\[\s*ms\s*\]|\(\s*ms\s*\)|_ms\b|\bmsec\b|millisec)", re.IGNORECASE)
_FS_MISMATCH_TOL = 0.01
_SEPARATOR_CANDIDATES = ("\t", ";", ",", "|")
_DECIMAL_COMMA_RE = re.compile(r"\d,\d")


class BioSignalLoaderError(Exception):
    """Raised for errors occurring during biosignal loading or parsing."""


@dataclass
class SignalData:
    """Standardized biosignal container.

    Attributes:
        signal: 1D float64 array with the PPG signal.
        sampling_rate: Sampling frequency in Hz.
        time: 1D float64 time vector in seconds (derived from ``sampling_rate`` if omitted).
        channel_name: Name of the extracted channel.
        file_format: Identified file format (e.g. ``"CSV"``, ``"MAT"``, ``"EDF"``).
        metadata: Non-identifying information about the source file.
        aux_signals: Optional auxiliary channels sampled on the same clock
            (e.g. ``{"ecg": array}``), used for pulse arrival time.
    """

    signal: FloatArray
    sampling_rate: float
    time: FloatArray = None  # type: ignore[assignment]  # None -> uniform time base
    channel_name: str = "Channel_0"
    file_format: str = ""
    metadata: Dict[str, Any] = field(default_factory=dict)
    aux_signals: Dict[str, FloatArray] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.signal = _to_finite_1d(self.signal, "Signal")
        if self.sampling_rate is None or not np.isfinite(self.sampling_rate) or self.sampling_rate <= 0:
            raise ValueError(f"sampling_rate must be a positive number, got {self.sampling_rate!r}")
        self.sampling_rate = float(self.sampling_rate)

        if self.time is None:
            self.time = np.arange(len(self.signal)) / self.sampling_rate
        else:
            self.time = np.asarray(self.time, dtype=np.float64)
            if self.time.ndim != 1:
                raise ValueError("Time vector must be a 1D array.")
            if len(self.time) != len(self.signal):
                raise ValueError(
                    f"Time vector length ({len(self.time)}) does not match signal length ({len(self.signal)})."
                )

        if self.metadata is None:
            self.metadata = {}
        clean_aux: Dict[str, np.ndarray] = {}
        for name, values in (self.aux_signals or {}).items():
            arr = _to_finite_1d(values, f"Auxiliary signal '{name}'")
            if len(arr) != len(self.signal):
                raise ValueError(f"Auxiliary signal '{name}' length ({len(arr)}) differs from the PPG signal.")
            clean_aux[name] = arr
        self.aux_signals = clean_aux

    @property
    def duration(self) -> float:
        return len(self.signal) / self.sampling_rate


def _to_finite_1d(values: Any, what: str) -> np.ndarray:
    arr = np.array(values, dtype=np.float64)
    if arr.ndim != 1:
        raise ValueError(f"{what} must be a 1D array, got shape {arr.shape}")
    bad = ~np.isfinite(arr)
    if np.any(bad):
        logger.warning("%s contains %d NaN/Inf values; imputing via linear interpolation.", what, int(bad.sum()))
        if np.all(bad):
            return np.zeros_like(arr)
        x = np.arange(len(arr))
        arr[bad] = np.interp(x[bad], x[~bad], arr[~bad])
    return arr


def estimate_sampling_rate(time: FloatArray) -> Tuple[float, float]:
    """Estimates fs from a time vector.

    The sample period is the least-squares slope of time against sample
    index. Unlike the median of the differences, this is unaffected by
    rounding of the time stamps (e.g. BIDMC stores 0.008 s steps with only
    two decimals after ~76 s, so the median difference would give 100 Hz
    instead of 125 Hz).

    Returns:
        (fs, max_gap): the sampling rate and the largest deviation of a time
        stamp from the fitted uniform grid, in sample periods. Values above
        ~1 indicate gaps or non-uniform sampling. Rates within 1 ppm of an
        integer are rounded to it.
    """
    t = np.asarray(time, dtype=np.float64)
    ok = np.isfinite(t)
    if ok.sum() < 2:
        raise BioSignalLoaderError("Time column has fewer than two finite samples.")
    idx = np.flatnonzero(ok).astype(np.float64)
    slope, intercept = np.polyfit(idx, t[ok], 1)
    if slope <= 0:
        raise BioSignalLoaderError("Time column is not increasing.")
    max_gap = float(np.max(np.abs(t[ok] - (slope * idx + intercept))) / slope)
    fs = 1.0 / float(slope)
    if abs(fs - round(fs)) < 1e-6 * fs:  # drop floating-point residue of the fit (500.0000000000016 Hz)
        fs = float(round(fs))
    return fs, max_gap


def _resolve_sampling_rate(
    requested: Optional[float], time_vector: Optional[np.ndarray], metadata: Dict[str, Any]
) -> float:
    """Reconciles an explicit sampling rate with one estimated from a time column."""
    estimated = None
    if time_vector is not None and len(time_vector) >= 2:
        estimated, max_gap = estimate_sampling_rate(time_vector)
        metadata["sampling_rate_estimated"] = estimated
        metadata["time_max_deviation_samples"] = max_gap
        if max_gap > 1.0:
            logger.warning("Time stamps deviate from a uniform grid by up to %.1f samples (gaps?).", max_gap)

    if requested is None:
        if estimated is None:
            raise BioSignalLoaderError(
                "No time column found: the sampling rate cannot be inferred. Pass sampling_rate explicitly."
            )
        metadata["sampling_rate_source"] = "time_column"
        return estimated

    requested = float(requested)
    metadata["sampling_rate_source"] = "user"
    if estimated is not None and abs(estimated - requested) / estimated > _FS_MISMATCH_TOL:
        message = (
            f"Sampling rate set to {requested:g} Hz but the time column implies {estimated:.3f} Hz. "
            f"All time-based features will be scaled by {estimated / requested:.3f}."
        )
        metadata["sampling_rate_warning"] = message
        logger.warning(message)
    return requested


def _find_column(columns: List[str], wanted: Union[int, str], role: str) -> str:
    """Finds a column by index, exact (stripped) name or case-insensitive name."""
    if isinstance(wanted, (int, np.integer)):
        if wanted < 0 or wanted >= len(columns):
            raise BioSignalLoaderError(f"{role} column index {wanted} out of bounds. File has {len(columns)} columns.")
        return columns[int(wanted)]
    if isinstance(wanted, str):
        key = wanted.strip()
        if key in columns:
            return key
        lowered = {c.lower(): c for c in columns}
        if key.lower() in lowered:
            return lowered[key.lower()]
        if key.isdigit():
            return _find_column(columns, int(key), role)
        raise BioSignalLoaderError(f"{role} column '{wanted}' not found. Available: {columns}")
    raise BioSignalLoaderError(f"{role} column must be int or str, got {type(wanted)}")


def _read_head(file_input: Any, n_bytes: int = 65536) -> str:
    """Returns the first ``n_bytes`` of a path or buffer as text, rewinding buffers."""
    if isinstance(file_input, (str, pathlib.Path)):
        with open(file_input, "rb") as fh:
            raw = fh.read(n_bytes)
    else:
        _rewind(file_input)
        raw = file_input.read(n_bytes)
        _rewind(file_input)
    text = raw.decode("utf-8", errors="replace") if isinstance(raw, bytes) else str(raw)
    return text.lstrip("\ufeff")


def sniff_csv_dialect(head: str, delimiter: Optional[str] = None) -> Tuple[str, str]:
    r"""Guesses the column separator and the decimal mark of a delimited text file.

    The separator is the first of tab, ``;``, ``,`` and ``|`` that occurs the
    same (non-zero) number of times on every sampled line; files with none of
    them are split on whitespace. A comma between digits marks a decimal comma
    (e.g. ``0,79;73689`` from spreadsheets with an Italian or German locale),
    which is only possible when the separator is not itself a comma.

    Args:
        head: The first lines of the file.
        delimiter: A known separator; only the decimal mark is then guessed.

    Returns:
        ``(separator, decimal)``; the separator may be the regex ``r"\s+"``.
    """
    lines = [ln for ln in head.splitlines() if ln.strip()]
    if lines and not head.endswith("\n") and len(lines) > 1:
        lines = lines[:-1]  # the last line may be cut in the middle
    lines = lines[:50]
    sep = delimiter
    if sep is None:
        sep = r"\s+"
        for cand in _SEPARATOR_CANDIDATES:
            counts = {ln.count(cand) for ln in lines}
            if len(counts) == 1 and counts.pop() > 0:
                sep = cand
                break
        else:
            if len(lines) > 0 and all(len(ln.split()) == 1 for ln in lines):
                sep = ","  # single column
    decimal = "."
    if sep != ",":
        data_lines = lines[1:] if len(lines) > 1 else lines
        if any(_DECIMAL_COMMA_RE.search(ln) for ln in data_lines):
            decimal = ","
    return sep, decimal


def _rewind(file_input: Any) -> None:
    if hasattr(file_input, "seek"):
        try:
            file_input.seek(0)
        except Exception:  # noqa: BLE001 - some streams are not seekable
            pass


class BioSignalLoader:
    """Static loaders for the supported file formats."""

    @staticmethod
    def _is_stream(file_input: Any) -> bool:
        return hasattr(file_input, "read")

    # ------------------------------------------------------------------ #
    # CSV / TXT / TSV
    # ------------------------------------------------------------------ #
    @classmethod
    def load_csv(
        cls,
        file_input: FileInput,
        signal_column: Optional[Union[int, str]] = None,
        sampling_rate: Optional[float] = None,
        delimiter: Optional[str] = None,
        ecg_column: Optional[Union[int, str]] = None,
        time_column: Optional[Union[int, str]] = None,
    ) -> SignalData:
        """Loads a biosignal from a delimited text file.

        Args:
            file_input: Path or file-like buffer.
            signal_column: Name or index of the PPG column. Names are matched after
                stripping whitespace, then case-insensitively. ``None`` selects the
                first column that is not a time column.
            sampling_rate: Sampling frequency in Hz, or ``None`` to infer it from the time column.
            delimiter: Column separator; ``None`` (default) sniffs it together with the
                decimal mark (see :func:`sniff_csv_dialect`).
            ecg_column: Optional name or index of an ECG column (for pulse arrival time).
            time_column: Optional explicit time column; otherwise detected by name.
        """
        try:
            sep, decimal = sniff_csv_dialect(_read_head(file_input), delimiter)
            _rewind(file_input)
            read_kwargs: Dict[str, Any] = {"sep": sep, "decimal": decimal}
            if sep == r"\s+":
                read_kwargs["engine"] = "python"
            df = pd.read_csv(file_input, header="infer", **read_kwargs)
            if df.empty:
                raise BioSignalLoaderError("The file contains no data rows.")

            # A numeric first "header" means the file has no header row.
            try:
                float(str(df.columns[0]))
                has_header = False
            except ValueError:
                has_header = True
            if not has_header:
                _rewind(file_input)
                df = pd.read_csv(file_input, header=None, **read_kwargs)
                df.columns = [f"col_{i}" for i in range(df.shape[1])]

            df.columns = [str(c).strip() for c in df.columns]
            if df.columns.duplicated().any():
                logger.warning("Duplicate column names found; keeping the first occurrence of each.")
                df = df.loc[:, ~df.columns.duplicated()]
            columns = list(df.columns)

            if signal_column is None:
                candidates = [c for c in columns if not _TIME_COLUMN_RE.match(c)] or columns
                col_name = candidates[0]
            else:
                col_name = _find_column(columns, signal_column, "Signal")
            signal = pd.to_numeric(df[col_name], errors="coerce").to_numpy()

            if time_column is not None:
                time_col: Optional[str] = _find_column(columns, time_column, "Time")
            else:
                time_col = next((c for c in columns if _TIME_COLUMN_RE.match(c) and c != col_name), None)
            time_vector = None
            metadata: Dict[str, Any] = {"columns": columns, "shape": df.shape, "signal_column_resolved": col_name,
                                        "separator": sep, "decimal": decimal}
            if time_col is not None and time_col != col_name:
                time_vector = pd.to_numeric(df[time_col], errors="coerce").to_numpy(dtype=np.float64)
                if _MS_UNIT_RE.search(time_col):
                    time_vector = time_vector / 1000.0
                if not np.all(np.isfinite(time_vector)):
                    logger.warning("Time column '%s' has non-numeric values; ignoring it.", time_col)
                    time_vector = None
                else:
                    metadata["time_column"] = time_col

            fs = _resolve_sampling_rate(sampling_rate, time_vector, metadata)
            if time_vector is not None:
                time_vector = time_vector - time_vector[0]

            aux: Dict[str, np.ndarray] = {}
            if ecg_column is not None:
                ecg_name = _find_column(columns, ecg_column, "ECG")
                aux["ecg"] = pd.to_numeric(df[ecg_name], errors="coerce").to_numpy()
                metadata["ecg_column_resolved"] = ecg_name

            return SignalData(
                signal=signal,
                sampling_rate=fs,
                time=time_vector,
                channel_name=col_name,
                file_format="CSV",
                metadata=metadata,
                aux_signals=aux,
            )
        except BioSignalLoaderError:
            raise
        except Exception as e:
            raise BioSignalLoaderError(f"Error loading delimited text file: {e}") from e

    # ------------------------------------------------------------------ #
    # MAT
    # ------------------------------------------------------------------ #
    @classmethod
    def load_mat(
        cls,
        file_input: FileInput,
        key: Optional[str] = None,
        sampling_rate: Optional[float] = None,
        channel_index: int = 0,
    ) -> SignalData:
        """Loads a biosignal from a MATLAB (.mat, v5-v7.2) file.

        Args:
            file_input: Path or in-memory buffer.
            key: Variable holding the signal. If ``None``, the first numeric array
                that is not a time vector is used.
            sampling_rate: Sampling frequency in Hz, or ``None`` to infer it from a
                ``time``/``t`` variable, or from an ``fs``/``Fs``/``sampling_rate`` scalar.
            channel_index: Channel to extract when the variable is 2D.
        """
        try:
            _rewind(file_input)
            try:
                mat_dict = scipy.io.loadmat(file_input)
            except Exception as e:
                raise BioSignalLoaderError(f"scipy.io.loadmat failed to read MAT file: {e}") from e

            user_keys = {k: v for k, v in mat_dict.items() if not k.startswith("__")}
            if not user_keys:
                raise BioSignalLoaderError("No non-system variables found in MAT file.")

            time_keys = [k for k in user_keys if _TIME_COLUMN_RE.fullmatch(k)]
            fs_keys = [k for k in user_keys if k.lower() in ("fs", "sampling_rate", "srate", "sample_rate")]

            selected_key = key
            if selected_key is None:
                for k, v in user_keys.items():
                    if k in time_keys or k in fs_keys:
                        continue
                    if isinstance(v, np.ndarray) and np.issubdtype(v.dtype, np.number) and v.size > 1:
                        selected_key = k
                        break
                if selected_key is None:
                    raise BioSignalLoaderError("Could not find any numeric signal array in MAT file.")
            elif selected_key not in user_keys:
                raise BioSignalLoaderError(
                    f"Specified key '{selected_key}' not found. Available keys: {list(user_keys.keys())}"
                )

            data_array = user_keys[selected_key]
            if not isinstance(data_array, np.ndarray):
                raise BioSignalLoaderError(f"Data for key '{selected_key}' is not a NumPy array.")

            squeezed = data_array.squeeze()
            if squeezed.ndim == 1:
                signal = squeezed
            elif squeezed.ndim == 2:
                rows, cols = squeezed.shape
                along_rows = rows <= cols  # channels are the shorter dimension
                n_channels = rows if along_rows else cols
                if channel_index < 0 or channel_index >= n_channels:
                    raise BioSignalLoaderError(
                        f"Channel index {channel_index} out of bounds for MAT array of shape {squeezed.shape}."
                    )
                signal = squeezed[channel_index, :] if along_rows else squeezed[:, channel_index]
            else:
                raise BioSignalLoaderError(
                    f"MAT variable '{selected_key}' has shape {data_array.shape}; only 1D/2D arrays are supported."
                )

            metadata: Dict[str, Any] = {
                "all_keys": list(user_keys.keys()),
                "selected_key": selected_key,
                "original_shape": data_array.shape,
                "channel_index": channel_index,
            }
            time_vector = None
            if time_keys:
                t_arr = np.asarray(user_keys[time_keys[0]], dtype=np.float64).squeeze()
                if t_arr.ndim == 1 and len(t_arr) == len(signal):
                    time_vector = t_arr
                    metadata["time_key"] = time_keys[0]
            if sampling_rate is None and time_vector is None and fs_keys:
                sampling_rate = float(np.asarray(user_keys[fs_keys[0]]).squeeze())
                metadata["sampling_rate_key"] = fs_keys[0]

            fs = _resolve_sampling_rate(sampling_rate, time_vector, metadata)
            if time_vector is not None:
                time_vector = time_vector - time_vector[0]

            return SignalData(
                signal=signal,
                sampling_rate=fs,
                time=time_vector,
                channel_name=f"{selected_key}_ch{channel_index}" if squeezed.ndim > 1 else selected_key,
                file_format="MAT",
                metadata=metadata,
            )
        except BioSignalLoaderError:
            raise
        except Exception as e:
            raise BioSignalLoaderError(f"Error loading MAT file: {e}") from e

    # ------------------------------------------------------------------ #
    # EDF / BDF
    # ------------------------------------------------------------------ #
    @classmethod
    def load_edf(
        cls,
        file_input: FileInput,
        channel_index: Union[int, str] = 0,
        ecg_channel: Optional[Union[int, str]] = None,
        suffix: str = ".edf",
    ) -> SignalData:
        """Loads a channel from an EDF/EDF+/BDF file.

        Patient-identifying header fields (name, code, birth date, admin code,
        technician, recording notes and start date) are never read, so they
        cannot leak into exports.

        Args:
            file_input: Path or in-memory buffer.
            channel_index: Index or label of the PPG channel.
            ecg_channel: Optional index or label of an ECG channel with the same sampling rate.
            suffix: File suffix used for the temporary copy of in-memory buffers.
        """
        file_path, is_temp = "", False
        try:
            file_path, is_temp = cls._to_physical_path(file_input, suffix)
            try:
                import pyedflib  # type: ignore
            except ImportError:
                pyedflib = None  # type: ignore[assignment]

            if pyedflib is not None:
                try:
                    reader = pyedflib.EdfReader(file_path)
                except Exception as e:
                    raise BioSignalLoaderError(f"pyedflib failed to open file: {e}") from e
                try:
                    labels = [str(lbl).strip() for lbl in reader.getSignalLabels()]
                    idx = cls._edf_channel(labels, channel_index)
                    signal = reader.readSignal(idx)
                    fs = float(reader.getSampleFrequency(idx))
                    header = reader.getSignalHeader(idx)
                    aux: Dict[str, np.ndarray] = {}
                    if ecg_channel is not None:
                        ecg_idx = cls._edf_channel(labels, ecg_channel)
                        if float(reader.getSampleFrequency(ecg_idx)) != fs:
                            raise BioSignalLoaderError("ECG and PPG channels have different sampling rates.")
                        aux["ecg"] = reader.readSignal(ecg_idx)
                    metadata = {
                        "n_channels": len(labels),
                        "channel_labels": labels,
                        "signal_header": {k: header.get(k) for k in ("label", "dimension", "sample_frequency",
                                                                     "sample_rate", "prefilter", "transducer")
                                          if k in header},
                        "duration_sec": float(reader.getFileDuration()),
                    }
                finally:
                    reader.close()
            else:
                try:
                    import mne  # type: ignore
                except ImportError as e:
                    raise BioSignalLoaderError("EDF parsing requires 'pyedflib' or 'mne'.") from e
                try:
                    raw = mne.io.read_raw_bdf(file_path, preload=True, verbose=False) if suffix.lower() == ".bdf" \
                        else mne.io.read_raw_edf(file_path, preload=True, verbose=False)
                except Exception as e:
                    raise BioSignalLoaderError(f"MNE failed to read file: {e}") from e
                labels = list(raw.info["ch_names"])
                idx = cls._edf_channel(labels, channel_index)
                fs = float(raw.info["sfreq"])
                signal = raw.get_data(picks=[idx])[0]
                aux = {}
                if ecg_channel is not None:
                    aux["ecg"] = raw.get_data(picks=[cls._edf_channel(labels, ecg_channel)])[0]
                metadata = {"n_channels": len(labels), "channel_labels": labels}

            return SignalData(
                signal=signal,
                sampling_rate=fs,
                channel_name=labels[idx],
                file_format="BDF" if suffix.lower() == ".bdf" else "EDF",
                metadata=metadata,
                aux_signals=aux,
            )
        except BioSignalLoaderError:
            raise
        except Exception as e:
            raise BioSignalLoaderError(f"Error loading EDF/BDF file: {e}") from e
        finally:
            if is_temp and file_path and os.path.exists(file_path):
                try:
                    os.unlink(file_path)
                except OSError as e:
                    logger.warning("Failed to delete temporary file '%s': %s", file_path, e)

    @staticmethod
    def _edf_channel(labels: List[str], wanted: Union[int, str]) -> int:
        if isinstance(wanted, str) and not wanted.strip().isdigit():
            name = _find_column(labels, wanted, "Channel")
            return labels.index(name)
        idx = int(wanted)
        if idx < 0 or idx >= len(labels):
            raise BioSignalLoaderError(f"Channel index {idx} out of bounds. File has {len(labels)} channels.")
        return idx

    @classmethod
    def _to_physical_path(cls, file_input: FileInput, suffix: str) -> Tuple[str, bool]:
        """Returns a file path, writing in-memory buffers to a temporary file."""
        if isinstance(file_input, (str, pathlib.Path)):
            return str(file_input), False
        if cls._is_stream(file_input):
            _rewind(file_input)
            content = file_input.read()
            if isinstance(content, str):
                content = content.encode("utf-8")
            with tempfile.NamedTemporaryFile(delete=False, suffix=suffix) as tmp:
                tmp.write(content)
            return tmp.name, True
        raise BioSignalLoaderError(f"Unsupported input type: {type(file_input)}")

    # ------------------------------------------------------------------ #
    # Dispatcher
    # ------------------------------------------------------------------ #
    @classmethod
    def auto_load(
        cls,
        file_input: FileInput,
        file_name: Optional[str] = None,
        signal_column: Optional[Union[int, str]] = None,
        sampling_rate: Optional[float] = None,
        ecg_column: Optional[Union[int, str]] = None,
    ) -> SignalData:
        """Identifies the format from the file extension and loads the signal.

        Args:
            file_input: Path or buffer.
            file_name: Name used to detect the format when ``file_input`` is a buffer.
            signal_column: Column (CSV), channel index (MAT) or channel index/label (EDF/BDF).
                ``None`` selects the first non-time column (CSV) or channel 0.
            sampling_rate: Sampling rate in Hz, or ``None`` to infer it. Ignored for EDF/BDF,
                whose header always carries the rate.
            ecg_column: Optional ECG column/channel.
        """
        resolved_name = file_name or (str(file_input) if isinstance(file_input, (str, pathlib.Path)) else
                                      getattr(file_input, "name", ""))
        if not resolved_name:
            raise BioSignalLoaderError(
                "Cannot identify the file format: no file name was provided and the buffer has no 'name'."
            )
        ext = pathlib.Path(str(resolved_name)).suffix.lower()

        if ext in (".csv", ".tsv", ".txt"):
            delimiter = "\t" if ext == ".tsv" else None
            return cls.load_csv(file_input, signal_column=signal_column, sampling_rate=sampling_rate,
                                delimiter=delimiter, ecg_column=ecg_column)
        if signal_column is None and ext not in (".csv", ".tsv", ".txt"):
            signal_column = 0
        if ext == ".mat":
            channel = int(signal_column) if str(signal_column).strip().isdigit() else 0
            key = None if str(signal_column).strip().isdigit() else str(signal_column)
            return cls.load_mat(file_input, key=key, sampling_rate=sampling_rate, channel_index=channel)
        if ext in (".edf", ".bdf"):
            return cls.load_edf(file_input, channel_index=signal_column, ecg_channel=ecg_column, suffix=ext)
        raise BioSignalLoaderError(f"Unsupported file extension '{ext}' in file '{resolved_name}'")


SUPPORTED_EXTENSIONS: Tuple[str, ...] = ("csv", "tsv", "txt", "mat", "edf", "bdf")
