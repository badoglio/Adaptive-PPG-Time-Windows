import json

import pandas as pd
import pytest

from adaptive_ppg.cli import main
from conftest import DATA


def test_init_config_and_run(tmp_path):
    cfg = tmp_path / "cfg.json"
    assert main(["init-config", str(cfg), "--preset", "rest"]) == 0
    assert json.loads(cfg.read_text())["engine"]["theta"] == 0.5
    out = tmp_path / "out"
    rc = main(["run", str(DATA / "bidmc_01_Signals.csv"), "--signal-column", "PLETH", "--ecg-column", "II",
               "--config", str(cfg), "--fixed", "30", "--templates", "--out", str(out)])
    assert rc == 0
    assert (out / "features_grid.csv").exists() and (out / "templates_macro.csv").exists()
    meta = json.loads((out / "metadata.json").read_text())
    assert meta["config"]["engine"]["theta"] == 0.5


def test_batch_run_reports_failures(tmp_path):
    other = tmp_path / "other.csv"
    pd.DataFrame({"Time [s]": [0.0, 0.01, 0.02], "PPG": [1.0, 2.0, 3.0]}).to_csv(other, index=False)
    out = tmp_path / "batch"
    rc = main(["run", str(DATA / "bidmc_01_Signals.csv"), str(other), "--signal-column", "PLETH",
               "--out", str(out)])
    assert rc == 1  # other.csv has no PLETH column
    summary = pd.read_csv(out / "batch_summary.csv")
    assert list(summary["status"].str.startswith("ok")) == [True, False]


def test_dashboard_passes_options_to_streamlit(monkeypatch):
    import adaptive_ppg.cli as cli

    calls = []
    monkeypatch.setattr(cli.subprocess, "call", lambda cmd: calls.append(cmd) or 0)
    assert main(["dashboard", "--server.port", "8502"]) == 0
    assert calls[0][-2:] == ["--server.port", "8502"]


def test_unknown_option_is_rejected_outside_dashboard():
    with pytest.raises(SystemExit):
        main(["init-config", "x.json", "--server.port", "8502"])


def test_run_overrides_processing_rate_and_polarity(tmp_path):
    out = tmp_path / "native"
    rc = main(["run", str(DATA / "Sample1.CSV"), "--target-fs", "0", "--polarity", "inverted", "--fixed", "10",
               "--out", str(out)])
    assert rc == 0
    meta = json.loads((out / "metadata.json").read_text())
    assert meta["config"]["preprocessing"]["target_fs"] is None
    assert meta["config"]["preprocessing"]["polarity"] == "inverted"
    assert meta["summary"]["processing_sampling_rate_hz"] == 500.0
    assert meta["summary"]["polarity"] == "inverted"
