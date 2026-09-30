"""Regression on the bundled real recordings (PhysioNet PTT PPG dataset, BIDMC)."""

import numpy as np
import pytest

from adaptive_ppg.io import BioSignalLoader
from adaptive_ppg.pipeline import run
from conftest import DATA


@pytest.fixture(scope="module", params=["Sample1.CSV", "Sample2.CSV"])
def ptt_result(request):
    return run(BioSignalLoader.auto_load(DATA / request.param), fixed_window_secs=[30])


def test_ptt_samples_are_inverted_and_processed_at_125_hz(ptt_result):
    ppg = ptt_result.ppg
    assert ppg.polarity == -1 and ppg.polarity_skewness < -1
    assert ppg.original_sampling_rate == pytest.approx(500.0, rel=1e-3)
    assert ppg.sampling_rate == pytest.approx(125.0)


def test_ptt_samples_quality_and_plausible_features(ptt_result):
    ppg, fr = ptt_result.ppg, ptt_result.features
    assert ppg.valid_beats_mask.mean() >= 0.99
    macro = fr.adaptive["macro"].table
    assert macro["template_sqi"].median() > 0.95
    g = fr.grid
    assert 55 < g["heart_rate"].median() < 90
    assert 0.08 < g["crest_time"].median() < 0.3
    assert 0.1 < g["duty_cycle"].median() < 0.45
    assert g["sdptg_b_a"].median() < 0 < g["sdptg_e_a"].median()
    assert g["dc_level"].median() > 1e4  # raw counts: the DC level keeps the original sign and offset
    assert 0 < g["perfusion_index"].median() < 20


def test_bidmc_is_upright_and_not_resampled(bidmc):
    res = run(bidmc, fixed_window_secs=[])
    assert res.ppg.polarity == 1 and res.ppg.sampling_rate == pytest.approx(125.0)
    assert np.isfinite(res.features.grid["pat_foot"]).mean() > 0.9
