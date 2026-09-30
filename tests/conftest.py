import pathlib

import pytest

from adaptive_ppg.config import PipelineConfig
from adaptive_ppg.io import BioSignalLoader
from adaptive_ppg.preprocessing import PPGPreprocessor
from adaptive_ppg.synthetic import BeatParameters, Segment, generate, stress_scenario

DATA = pathlib.Path(__file__).resolve().parents[1] / "sample_data"


@pytest.fixture(scope="session")
def bidmc():
    return BioSignalLoader.auto_load(DATA / "bidmc_01_Signals.csv", signal_column="PLETH", ecg_column="II")


@pytest.fixture(scope="session")
def bidmc_ppg(bidmc):
    return PPGPreprocessor(PipelineConfig().preprocessing).process(bidmc)


@pytest.fixture(scope="session")
def stress_rec():
    return generate(stress_scenario(), fs=125.0, seed=1)


@pytest.fixture(scope="session")
def short_rec():
    """Two minutes with one transient, for the faster tests."""
    stress = BeatParameters(hr=95.0, sys_time=0.13, tidal_amp=0.15, refl_amp=0.25, pat=0.2, amplitude=0.7)
    return generate([Segment(60.0, BeatParameters(), 0.0, "rest"), Segment(60.0, stress, 4.0, "stress")],
                    fs=125.0, seed=3)
