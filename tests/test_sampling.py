"""Sampling-rate normalization, polarity detection and sub-sample fiducials."""

import numpy as np
import pytest

from adaptive_ppg.config import PipelineConfig
from adaptive_ppg.io import SignalData
from adaptive_ppg.pipeline import run
from adaptive_ppg.preprocessing import (
    PPGPreprocessor,
    _resample_mask,
    derivative_skewness,
    parabolic_offset,
    resample_signal,
)
from adaptive_ppg.synthetic import generate, rest_scenario


def _cfg(target_fs=125.0, polarity="auto"):
    cfg = PipelineConfig()
    cfg.preprocessing.target_fs = target_fs
    cfg.preprocessing.polarity = polarity
    return cfg


@pytest.fixture(scope="module")
def rest_500():
    return generate(rest_scenario(120.0), fs=500.0, seed=4)


# ---------------------------------------------------------------- resampling
@pytest.mark.parametrize("fs, target", [(500.0, 125.0), (1000.0, 125.0), (256.0, 125.0), (100.0, 125.0)])
def test_resample_rate_is_exact(fs, target):
    x = np.sin(2 * np.pi * 1.3 * np.arange(int(20 * fs)) / fs)
    y, fs_out = resample_signal(x, fs, target)
    assert fs_out == pytest.approx(target, rel=2e-3)
    assert len(y) == pytest.approx(len(x) * fs_out / fs, abs=1)
    t = np.arange(len(y)) / fs_out
    inner = (t > 1) & (t < 19)  # the time base is exact: the sine is preserved sample by sample
    assert np.max(np.abs(y[inner] - np.sin(2 * np.pi * 1.3 * t[inner]))) < 1e-2


def test_resample_is_identity_within_tolerance_or_disabled():
    x = np.random.default_rng(0).standard_normal(1000)
    for target in (None, 125.0, 126.0):
        y, fs_out = resample_signal(x, 125.0, target)
        assert fs_out == 125.0 and np.array_equal(y, x)


def test_resample_mask_keeps_clipped_segments():
    mask = np.zeros(5000, dtype=bool)
    mask[1000:1100] = True
    out = _resample_mask(mask, 1250)
    assert out.dtype == bool and len(out) == 1250
    assert out[255:270].all() and not out[:240].any() and not out[290:].any()


def test_processing_rate_and_aux_channels_follow_target(rest_500):
    ppg = PPGPreprocessor(_cfg().preprocessing).process(rest_500.signal_data)
    assert ppg.sampling_rate == pytest.approx(125.0)
    assert ppg.original_sampling_rate == 500.0
    n = len(ppg.filtered_signal)
    assert n == pytest.approx(len(rest_500.signal_data.signal) / 4, abs=1)
    assert len(ppg.raw_signal) == n and len(ppg.aux_signals["ecg"]) == n


def test_upsampling_does_not_widen_the_filters():
    rec = generate(rest_scenario(60.0), fs=50.0, seed=2)
    cfg = _cfg(target_fs=125.0)
    cfg.preprocessing.highcut = 30.0  # above the 25 Hz Nyquist of the input
    ppg = PPGPreprocessor(cfg.preprocessing).process(rec.signal_data)
    assert ppg.sampling_rate == pytest.approx(125.0)
    assert ppg.effective_highcut <= 0.45 * 50.0 + 1e-9


def test_features_do_not_depend_on_the_processing_rate(rest_500):
    """A 500 Hz recording analysed natively and resampled to 125 Hz gives the same fixed-window features."""
    tables = {t: run(rest_500.signal_data, _cfg(t), fixed_window_secs=[20]).features.fixed["fixed_20s"].table
              for t in (None, 125.0)}
    assert len(tables[None]) == len(tables[125.0])
    for feat, tol in [("heart_rate", 0.1), ("crest_time", 0.001), ("pulse_width_50", 0.001),
                      ("duty_cycle", 0.005), ("sdptg_b_a", 0.05), ("reflection_index", 0.01),
                      ("notch_time", 0.002)]:
        a, b = tables[None][feat].to_numpy(), tables[125.0][feat].to_numpy()
        assert np.nanmedian(np.abs(a - b)) < tol, feat


def test_timing_is_not_quantized_to_the_sample_period(rest_500):
    """Crest time on a 125 Hz signal takes values between multiples of 8 ms."""
    res = run(rest_500.signal_data, _cfg(), fixed_window_secs=[])
    ct = res.features.per_beat["crest_time"].dropna().to_numpy() * 125.0
    assert np.mean(np.abs(ct - np.round(ct)) > 0.05) > 0.5
    ppg = res.ppg
    assert np.all(np.abs(ppg.onset_positions - ppg.onsets) <= 0.5)
    assert np.all(np.abs(ppg.peak_positions - ppg.peaks) <= 0.5)
    assert np.allclose(ppg.beat_times, ppg.peak_positions / ppg.sampling_rate)


def test_parabolic_offset_recovers_the_vertex():
    x = np.arange(10.0)
    y = -((x - 4.3) ** 2)
    assert 4 + parabolic_offset(y, 4) == pytest.approx(4.3)
    assert parabolic_offset(y, 0) == 0.0 and parabolic_offset(np.ones(5), 2) == 0.0


# ---------------------------------------------------------------- polarity
def test_derivative_skewness_sign(rest_500):
    x = rest_500.signal_data.signal
    assert derivative_skewness(x) > 0 > derivative_skewness(-x)


def test_inverted_signal_is_flipped_automatically(rest_500):
    sd = rest_500.signal_data
    inverted = SignalData(-sd.signal + 1e5, sd.sampling_rate, channel_name="counts", aux_signals=sd.aux_signals)
    upright = run(sd, _cfg(), fixed_window_secs=[])
    flipped = run(inverted, _cfg(), fixed_window_secs=[])
    assert upright.ppg.polarity == 1 and flipped.ppg.polarity == -1
    assert any("inverted automatically" in w for w in flipped.ppg.warnings)
    assert flipped.ppg.valid_beats_mask.mean() > 0.95
    assert np.array_equal(upright.ppg.onsets, flipped.ppg.onsets)
    assert np.allclose(upright.ppg.filtered_signal, flipped.ppg.filtered_signal, atol=1e-6)
    # raw_signal keeps the original polarity and DC level (used for DC level and perfusion index)
    assert np.median(flipped.ppg.raw_signal) > 5e4
    assert flipped.metadata()["summary"]["polarity"] == "inverted"


def test_forced_polarity_is_honoured(rest_500):
    forced = PPGPreprocessor(_cfg(polarity="inverted").preprocessing).process(rest_500.signal_data)
    assert forced.polarity == -1 and not any("automatically" in w for w in forced.warnings)
    normal = PPGPreprocessor(_cfg(polarity="normal").preprocessing).process(rest_500.signal_data)
    assert normal.polarity == 1
    assert np.allclose(forced.filtered_signal, -normal.filtered_signal)


def test_config_validation_of_new_fields():
    cfg = PipelineConfig()
    cfg.preprocessing.target_fs = 10.0
    with pytest.raises(ValueError, match="target_fs"):
        cfg.validate()
    cfg = PipelineConfig()
    cfg.preprocessing.polarity = "upside-down"
    with pytest.raises(ValueError, match="polarity"):
        cfg.validate()
    cfg = PipelineConfig()
    cfg.features.template_fs = 20.0
    with pytest.raises(ValueError, match="template_fs"):
        cfg.validate()
