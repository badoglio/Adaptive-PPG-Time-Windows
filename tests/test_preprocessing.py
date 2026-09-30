import numpy as np

from adaptive_ppg.config import PipelineConfig
from adaptive_ppg.preprocessing import PPGPreprocessor


def test_bidmc_beats_are_mostly_valid(bidmc_ppg):
    assert bidmc_ppg.n_beats > 650
    assert bidmc_ppg.valid_beats_mask.mean() >= 0.97
    assert np.all(bidmc_ppg.onsets < bidmc_ppg.peaks)
    assert np.all(bidmc_ppg.peaks < bidmc_ppg.ends)


def test_synthetic_detection_matches_ground_truth(stress_rec):
    ppg = PPGPreprocessor(PipelineConfig().preprocessing).process(stress_rec.signal_data)
    gt = stress_rec.beats
    assert abs(ppg.n_beats - len(gt)) <= 0.03 * len(gt)
    # every clean ground-truth systolic peak has a detected peak within 40 ms
    clean = gt[~gt["in_artifact"]]["t_peak"].to_numpy()
    detected = ppg.peaks / ppg.sampling_rate
    nearest = np.min(np.abs(detected[None, :] - clean[:, None]), axis=1)
    assert np.mean(nearest < 0.04) > 0.97


def test_artefacts_are_rejected(stress_rec):
    ppg = PPGPreprocessor(PipelineConfig().preprocessing).process(stress_rec.signal_data)
    t_peak = ppg.peaks / ppg.sampling_rate
    in_art = np.zeros(ppg.n_beats, dtype=bool)
    for a, b in stress_rec.artifacts:
        in_art |= (t_peak >= a) & (t_peak <= b)
    assert in_art.any()
    assert ppg.valid_beats_mask[in_art].mean() < ppg.valid_beats_mask[~in_art].mean()


def test_highcut_above_nyquist_is_clamped(short_rec):
    cfg = PipelineConfig().preprocessing
    cfg.highcut = 100.0  # above Nyquist at 125 Hz
    ppg = PPGPreprocessor(cfg).process(short_rec.signal_data)
    assert ppg.effective_highcut < 62.5
    assert ppg.warnings
