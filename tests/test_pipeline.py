import io
import json
import zipfile

import pytest

from adaptive_ppg.pipeline import run


@pytest.fixture(scope="module")
def result(short_rec):
    return run(short_rec.signal_data, fixed_window_secs=[10.0, 30.0])


def test_zip_export(result):
    zf = zipfile.ZipFile(io.BytesIO(result.to_zip(include_templates=True)))
    names = set(zf.namelist())
    for n in ("features_grid.csv", "windows.csv", "beats.csv", "features_per_beat.csv", "feature_dictionary.csv",
              "windows_fixed_10s.csv", "templates_macro.csv", "metadata.json"):
        assert n in names
    meta = json.loads(zf.read("metadata.json"))
    assert meta["version"] and meta["config"]["engine"]["variability_metric"] == "drift_z"
    assert meta["summary"]["n_beats"] == result.ppg.n_beats


def test_excel_export(result):
    pytest.importorskip("openpyxl")
    assert result.to_excel()[:2] == b"PK"


def test_grid_is_causal_and_complete(result):
    grid = result.features.grid
    assert grid["time"].is_monotonic_increasing
    first = result.features.adaptive["macro"].table["t_anchor"].min()
    for f in ("heart_rate", "peak_amplitude"):
        col = grid[f]
        assert col[grid["time"] < first].isna().all()
        assert col[grid["time"] > first + 1].notna().mean() > 0.95
