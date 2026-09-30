import numpy as np
import pytest

from adaptive_ppg.config import PipelineConfig
from adaptive_ppg.engine import AdaptiveEngine, StreamingAdaptiveEngine, drift
from adaptive_ppg.preprocessing import PPGPreprocessor


@pytest.fixture(scope="module")
def short_ppg(short_rec):
    return PPGPreprocessor(PipelineConfig().preprocessing).process(short_rec.signal_data)


def _stream(cfg, ppg):
    se = StreamingAdaptiveEngine(cfg)
    out = {c: [] for c in cfg.categories}
    fs = ppg.sampling_rate
    for k in range(ppg.n_beats):
        for w in se.push(ppg.beats_normalized[k], ppg.onset_positions[k] / fs, ppg.peak_positions[k] / fs,
                         ppg.end_positions[k] / fs,
                         bool(ppg.valid_beats_mask[k])):
            out[w.category].append(w)
    for w in se.finalize():
        out[w.category].append(w)
    return out


@pytest.mark.parametrize("metric", ["drift_z", "drift", "successive", "dispersion"])
@pytest.mark.parametrize("normalization", ["fixed", "calibration", "rolling"])
def test_streaming_reproduces_batch(short_ppg, metric, normalization):
    cfg = PipelineConfig().engine
    cfg.variability_metric = metric
    cfg.normalization = normalization
    cfg.calibration_sec = 30.0
    if metric != "drift_z":
        cfg.fixed_scale = (0.02, 0.15)
    batch = AdaptiveEngine(cfg).process(short_ppg)
    stream = _stream(cfg, short_ppg)
    for cat in cfg.categories:
        a = [(w.anchor_beat, tuple(w.beat_indices)) for w in batch.windows[cat]]
        b = [(w.anchor_beat, tuple(w.beat_indices)) for w in stream[cat]]
        assert a == b, cat


def test_drift_z_is_about_one_when_stationary():
    rng = np.random.default_rng(0)
    for noise in (0.01, 0.1):
        beats = 0.5 + noise * rng.standard_normal((400, 64))
        z = drift(beats, 10, standardized=True)
        assert np.nanmedian(z) == pytest.approx(1.0, abs=0.15)
    # a step change in shape is detected
    beats = 0.01 * rng.standard_normal((200, 64))
    beats[100:] += 0.05
    z = drift(beats, 10, standardized=True)
    assert np.nanmax(z[100:120]) > 5 * np.nanmedian(z[:90])


def test_windows_contract_during_transient_and_respect_limits(short_ppg):
    cfg = PipelineConfig().engine
    res = AdaptiveEngine(cfg).process(short_ppg)
    assert np.all((res.H >= 0) & (res.H <= 1))
    assert np.all(res.H[res.warmup_mask] == 0)
    t = res.beat_times
    sched = res.schedules["macro"]
    rest = np.median(sched.T_beats[(t > 30) & (t < 60)])
    transient = np.min(sched.T_beats[(t > 60) & (t < 75)])
    assert transient < rest
    lims = cfg.categories["macro"]
    assert sched.T_beats.min() >= lims.T_min_beats
    for w in sched.windows:
        assert w.n_valid >= lims.T_crit_beats
        assert short_ppg.valid_beats_mask[w.beat_indices].all()
        assert w.t_start <= w.t_anchor + 1e-9  # trailing windows are causal
