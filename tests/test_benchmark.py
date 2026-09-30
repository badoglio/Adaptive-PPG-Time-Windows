import pytest

from adaptive_ppg.benchmark import benchmark, feature_scales


@pytest.fixture(scope="module")
def bench(stress_rec):
    return benchmark(stress_rec, fixed_secs=(10.0, 60.0), ema_secs=(10.0,), scales=feature_scales(stress_rec))


def test_adaptive_tracks_transients_better_than_long_windows(bench):
    s = bench.summary.set_index("estimator")
    assert s.loc["adaptive", "median_nrmse"] < s.loc["fixed_60s", "median_nrmse"]
    assert s.loc["adaptive", "median_latency_s"] < s.loc["fixed_60s", "median_latency_s"]
    assert s.loc["adaptive", "median_stationary_nsd"] < s.loc["per_beat", "median_stationary_nsd"]


def test_metrics_table(bench):
    m = bench.metrics
    assert {"estimator", "feature", "rmse", "bias", "nrmse", "stationary_nsd", "latency_s", "coverage"} <= set(m)
    assert set(m["estimator"]) == {"adaptive", "fixed_10s", "fixed_60s", "ema_10s", "per_beat"}
    assert (m["coverage"] > 0.5).all()
