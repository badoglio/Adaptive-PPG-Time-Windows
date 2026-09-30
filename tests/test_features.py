import numpy as np
import pytest

from adaptive_ppg.config import PipelineConfig
from adaptive_ppg.features import ALL_FEATURES, FEATURE_TO_CATEGORY, FEATURE_UNITS, detect_r_peaks
from adaptive_ppg.pipeline import run
from adaptive_ppg.synthetic import NoiseConfig, generate, stress_scenario


@pytest.fixture(scope="module")
def bidmc_result(bidmc):
    return run(bidmc, fixed_window_secs=[30.0])


def test_feature_dictionary_is_consistent():
    assert set(ALL_FEATURES) == set(FEATURE_TO_CATEGORY) == set(FEATURE_UNITS)


def test_bidmc_morphology_is_plausible(bidmc_result):
    pb = bidmc_result.features.per_beat
    v = pb[pb["valid"]]
    assert 60 < v["heart_rate"].median() < 110
    assert v["sdptg_b_a"].median() < 0          # the b wave is negative
    assert 0.08 < v["crest_time"].median() < 0.35
    assert v["stiffness_index"].isna().all()    # no subject height
    macro = bidmc_result.features.adaptive["macro"].table
    assert macro["template_sqi"].median() > 0.95
    assert macro["t_anchor"].iloc[-1] > bidmc_result.signal_data.duration - 10  # windows cover the tail


def test_pat_from_ecg(bidmc_result):
    fr = bidmc_result.features
    assert len(fr.r_peaks) > 0.9 * bidmc_result.ppg.n_beats
    pat = fr.per_beat["pat_foot"].dropna()
    assert len(pat) > 0.8 * len(fr.per_beat) and 0.05 < pat.median() < 0.6


def test_r_peak_detector_on_synthetic_ecg():
    fs = 250.0
    rec = generate(stress_scenario(), fs=fs, seed=2)
    r = detect_r_peaks(rec.signal_data.aux_signals["ecg"], fs) / fs
    true = (rec.beats["t_onset"] - rec.beats["pat"]).to_numpy()
    assert abs(len(r) - len(true)) <= 2
    assert np.median(np.min(np.abs(r[None, :] - true[:, None]), axis=1)) < 0.01


def test_features_follow_the_scenario():
    quiet = NoiseConfig(snr_db=40.0, wander_amp=0.0, artifact_rate_per_min=0.0)
    res = run(generate(stress_scenario(), fs=125.0, noise=quiet, seed=0).signal_data, fixed_window_secs=[])
    pb = res.features.per_beat
    rest = pb[(pb["t_peak"] > 60) & (pb["t_peak"] < 115)]
    stress = pb[(pb["t_peak"] > 170) & (pb["t_peak"] < 235)]
    assert rest["heart_rate"].median() == pytest.approx(70, abs=3)
    assert stress["heart_rate"].median() == pytest.approx(100, abs=4)
    # PAT falls by 60 ms and the pulse amplitude by 40 %
    assert rest["pat_foot"].median() - stress["pat_foot"].median() == pytest.approx(0.06, abs=0.015)
    assert stress["peak_amplitude"].median() / rest["peak_amplitude"].median() == pytest.approx(0.6, abs=0.08)
    assert stress["sdptg_b_a"].median() < rest["sdptg_b_a"].median()


def test_stiffness_index_with_height(bidmc):
    cfg = PipelineConfig()
    cfg.features.subject_height_m = 1.75
    res = run(bidmc, cfg, fixed_window_secs=[])
    si = res.features.adaptive["time_volume"].table["stiffness_index"]
    assert si.notna().mean() > 0.5 and (si.dropna() > 0).all()
