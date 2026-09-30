import io

import numpy as np
import pandas as pd
import pytest

from adaptive_ppg.io import BioSignalLoader, BioSignalLoaderError, SignalData, estimate_sampling_rate, sniff_csv_dialect
from conftest import DATA


def test_bidmc_sampling_rate_despite_rounded_time_column(bidmc):
    # The time column is rounded to 2-3 decimals: the median difference would give 100 Hz.
    assert bidmc.sampling_rate == pytest.approx(125.0, rel=1e-4)
    assert bidmc.channel_name == "PLETH"
    assert "ecg" in bidmc.aux_signals and len(bidmc.aux_signals["ecg"]) == len(bidmc.signal)
    assert bidmc.duration == pytest.approx(480.0, abs=1.0)


def test_estimate_sampling_rate_detects_gaps():
    t = np.arange(1000) / 100.0
    fs, gap = estimate_sampling_rate(t)
    assert fs == pytest.approx(100.0) and gap < 1e-6
    t[500:] += 0.5  # 50-sample gap
    assert estimate_sampling_rate(t)[1] > 1.0


def test_csv_without_time_column_needs_fs(tmp_path):
    path = tmp_path / "no_time.csv"
    pd.DataFrame({"PPG": np.sin(np.arange(6000) / 50.0)}).to_csv(path, index=False)
    with pytest.raises(BioSignalLoaderError, match="sampling rate"):
        BioSignalLoader.auto_load(path, signal_column="PPG")
    sd = BioSignalLoader.auto_load(path, signal_column="PPG", sampling_rate=500.0)
    assert sd.sampling_rate == 500.0 and len(sd.signal) == 6000


def test_csv_buffer_ms_time_and_headerless():
    t = np.arange(500) * 4.0  # ms
    x = np.sin(2 * np.pi * 1.2 * t / 1000.0)
    text = "time_ms,ppg\n" + "\n".join(f"{a},{b}" for a, b in zip(t, x))
    buf = io.BytesIO(text.encode())
    buf.name = "rec.csv"
    sd = BioSignalLoader.auto_load(buf, signal_column="ppg")
    assert sd.sampling_rate == pytest.approx(250.0)
    headerless = io.StringIO("\n".join(f"{b}" for b in x))
    sd2 = BioSignalLoader.load_csv(headerless, signal_column=0, sampling_rate=250.0)
    assert len(sd2.signal) == 500


def test_user_fs_mismatch_is_reported():
    text = "time,ppg\n" + "\n".join(f"{i / 100},{np.sin(i / 10)}" for i in range(300))
    sd = BioSignalLoader.load_csv(io.StringIO(text), signal_column="ppg", sampling_rate=125.0)
    assert "sampling_rate_warning" in sd.metadata


def test_mat_round_trip(tmp_path):
    scipy_io = pytest.importorskip("scipy.io")
    x = np.random.default_rng(0).standard_normal(1000)
    path = tmp_path / "rec.mat"
    scipy_io.savemat(path, {"ppg": x, "fs": 64.0})
    sd = BioSignalLoader.auto_load(path)
    np.testing.assert_allclose(sd.signal, x)
    assert sd.sampling_rate == pytest.approx(64.0)


def test_edf_round_trip(tmp_path):
    highlevel = pytest.importorskip("pyedflib.highlevel")
    x = np.sin(np.arange(2560) / 20.0) * 100
    path = str(tmp_path / "rec.edf")
    headers = highlevel.make_signal_headers(["Pleth"], sample_frequency=128, physical_min=-200, physical_max=200)
    highlevel.write_edf(path, [x], headers)
    sd = BioSignalLoader.auto_load(path, signal_column="Pleth")
    assert sd.sampling_rate == pytest.approx(128.0)
    np.testing.assert_allclose(sd.signal, x, atol=0.05)


def test_unsupported_extension_and_bad_signal():
    with pytest.raises(BioSignalLoaderError):
        BioSignalLoader.auto_load(io.BytesIO(b""), file_name="x.wav")
    with pytest.raises(ValueError):
        SignalData(signal=np.array([1.0, 2.0]), sampling_rate=0.0)


# ---------------------------------------------------------------- CSV dialects
@pytest.mark.parametrize("head, expected", [
    ("Time (s);Pleth\n0,79;73689\n0,792;73689\n", (";", ",")),
    ("t,ppg\n0.0,1.5\n0.01,1.6\n", (",", ".")),
    ("t\tppg\n0,0\t1,5\n0,01\t1,6\n", ("\t", ",")),
    ("t ppg\n0.0 1.5\n0.01 1.6\n", (r"\s+", ".")),
    ("1.5\n1.6\n1.7\n", (",", ".")),
])
def test_sniff_csv_dialect(head, expected):
    assert sniff_csv_dialect(head) == expected


def test_semicolon_decimal_comma_csv_and_default_column():
    lines = ["Time (s);Pleth"] + [f"{i / 200:.3f};{1000 + 10 * np.sin(i / 20):.2f}".replace(".", ",")
                                  for i in range(4000)]
    buf = io.BytesIO(("\ufeff" + "\n".join(lines)).encode("utf-8"))
    buf.name = "REC.CSV"  # upper-case extension
    sd = BioSignalLoader.auto_load(buf)
    assert sd.channel_name == "Pleth"  # the time column is skipped by default
    assert sd.sampling_rate == pytest.approx(200.0)
    assert sd.metadata["separator"] == ";" and sd.metadata["decimal"] == ","
    assert sd.signal.max() == pytest.approx(1010.0, abs=0.01)


@pytest.mark.parametrize("name, n", [("Sample1.CSV", 140932), ("Sample2.CSV", 113094)])
def test_ptt_dataset_samples_load(name, n):
    sd = BioSignalLoader.auto_load(DATA / name)
    assert sd.sampling_rate == pytest.approx(500.0, rel=1e-3)
    assert sd.channel_name == "Pleth" and len(sd.signal) == n
