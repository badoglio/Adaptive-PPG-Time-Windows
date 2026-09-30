"""Synthetic PPG with known ground truth, for validation and benchmarking.

Each beat is the sum of three Gaussian waves (systolic, late-systolic "tidal"
and diastolic reflection, as in pulse-decomposition models), defined in
seconds from the pulse foot, with a linear foot-to-foot baseline removed so
that consecutive beats join continuously. The beat parameters
follow a scenario of stationary segments and autonomic transients (e.g.
rest -> acute stress -> recovery), so that every morphological feature has a
known trajectory. The clean signal is then corrupted with white noise,
baseline wander, respiratory amplitude modulation and motion-artefact bursts.

An optional synthetic ECG (one narrow QRS per beat, placed a known pulse
arrival time before the PPG foot) allows the PAT features to be validated.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from ._typing import BoolArray, FloatArray
from .io import SignalData


@dataclass
class BeatParameters:
    """Morphology of one synthetic beat (times in seconds from the foot)."""

    hr: float = 70.0            # heart rate (bpm)
    sys_time: float = 0.16      # centre of the systolic wave
    sys_width: float = 0.055    # SD of the systolic wave
    tidal_amp: float = 0.35     # late-systolic (tidal) wave / systolic amplitude
    tidal_delay: float = 0.09   # delay of the tidal wave after the systolic one
    tidal_width: float = 0.045  # SD of the tidal wave
    refl_amp: float = 0.45      # diastolic / systolic amplitude
    refl_delay: float = 0.24    # delay of the diastolic wave after the systolic one
    refl_width: float = 0.09    # SD of the diastolic wave
    pat: float = 0.25           # pulse arrival time, ECG R peak -> PPG foot (s)
    amplitude: float = 1.0      # pulse amplitude (AC)

    def as_array(self) -> np.ndarray:
        return np.array([getattr(self, f) for f in _PARAMS])


_PARAMS = ("hr", "sys_time", "sys_width", "tidal_amp", "tidal_delay", "tidal_width",
           "refl_amp", "refl_delay", "refl_width", "pat", "amplitude")


@dataclass
class Segment:
    """One scenario segment: parameters move from the previous state to ``target``.

    The transition starts at the segment start and follows an exponential
    approach with time constant ``tau`` (s); ``tau=0`` is a step.
    """

    duration: float
    target: BeatParameters
    tau: float = 0.0
    label: str = ""


def stress_scenario() -> List[Segment]:
    """Rest (120 s) -> acute stress (fast onset, 120 s) -> slow recovery (120 s)."""
    rest = BeatParameters()
    stress = BeatParameters(hr=100.0, sys_time=0.13, sys_width=0.045, tidal_amp=0.12, tidal_delay=0.08,
                            refl_amp=0.22, refl_delay=0.19,
                            refl_width=0.08, pat=0.19, amplitude=0.6)
    return [
        Segment(120.0, rest, 0.0, "rest"),
        Segment(120.0, stress, 4.0, "stress"),
        Segment(120.0, rest, 25.0, "recovery"),
    ]


def rest_scenario(duration: float = 300.0) -> List[Segment]:
    """Stationary recording (only respiratory modulation)."""
    return [Segment(duration, BeatParameters(), 0.0, "rest")]


@dataclass
class NoiseConfig:
    snr_db: float = 30.0                 # white noise, relative to the pulse RMS
    wander_amp: float = 0.3              # baseline wander amplitude (x pulse amplitude)
    wander_hz: Tuple[float, float] = (0.05, 0.15)
    resp_hz: float = 0.25                # respiration frequency
    resp_am: float = 0.08                # respiratory amplitude modulation (fraction)
    resp_hr: float = 3.0                 # respiratory sinus arrhythmia (bpm)
    hr_jitter: float = 0.01              # random beat-to-beat IBI jitter (fraction)
    artifact_rate_per_min: float = 1.0
    artifact_duration: Tuple[float, float] = (1.0, 4.0)
    artifact_amp: float = 2.0            # x pulse amplitude
    dc_level: float = 2.0


@dataclass
class SyntheticRecording:
    """A synthetic recording and its ground truth.

    Attributes:
        signal_data: Noisy signal, ready for the pipeline (``aux_signals['ecg']`` if requested).
        clean: Noise-free PPG (same DC level, no wander or artefacts).
        beats: Per-beat ground-truth table (foot/peak times and parameters).
        artifacts: List of ``(t_start, t_end)`` motion-artefact intervals.
        events: Times of the scenario transitions (s).
        stationary: ``(t_start, t_end)`` intervals where the parameters are settled.
    """

    signal_data: SignalData
    clean: FloatArray
    beats: pd.DataFrame
    artifacts: List[Tuple[float, float]]
    events: List[float]
    stationary: List[Tuple[float, float]]
    segments: List[Segment] = field(default_factory=list)

    @property
    def fs(self) -> float:
        return self.signal_data.sampling_rate

    def artifact_mask(self) -> BoolArray:
        t = self.signal_data.time
        mask = np.zeros(len(t), dtype=bool)
        for a, b in self.artifacts:
            mask |= (t >= a) & (t <= b)
        return mask


def _beat_waveform(tt: np.ndarray, p: np.ndarray) -> np.ndarray:
    hr, st, sw, ta, td, tw, ra, rd, rw, _pat, amp = p
    w = (np.exp(-0.5 * ((tt - st) / sw) ** 2)
         + ta * np.exp(-0.5 * ((tt - st - td) / tw) ** 2)
         + ra * np.exp(-0.5 * ((tt - st - rd) / rw) ** 2))
    return amp * w


def _parameter_trajectory(segments: Sequence[Segment], t: np.ndarray) -> np.ndarray:
    """Parameter vector at each time in ``t`` (rows) following the scenario."""
    out = np.empty((len(t), len(_PARAMS)))
    state = segments[0].target.as_array()
    start = 0.0
    for seg in segments:
        target = seg.target.as_array()
        sel = (t >= start) & (t < start + seg.duration)
        dt = t[sel] - start
        if seg.tau > 0:
            out[sel] = target + (state - target) * np.exp(-dt / seg.tau)[:, None]
        else:
            out[sel] = target
        # state at the end of the segment
        state = target + (state - target) * (np.exp(-seg.duration / seg.tau) if seg.tau > 0 else 0.0)
        start += seg.duration
    out[t >= start] = state
    return out


def generate(
    segments: Optional[Sequence[Segment]] = None,
    fs: float = 125.0,
    noise: Optional[NoiseConfig] = None,
    with_ecg: bool = True,
    seed: Optional[int] = 0,
) -> SyntheticRecording:
    """Generates a synthetic recording following ``segments`` (default: :func:`stress_scenario`)."""
    segments = list(segments or stress_scenario())
    noise = noise or NoiseConfig()
    rng = np.random.default_rng(seed)
    duration = float(sum(s.duration for s in segments))
    n = int(round(duration * fs))
    t = np.arange(n) / fs

    # -- beat train --------------------------------------------------- #
    feet = [0.3]
    params = []
    while True:
        tf = feet[-1]
        p = _parameter_trajectory(segments, np.array([tf]))[0]
        p[0] += noise.resp_hr * np.sin(2 * np.pi * noise.resp_hz * tf)
        p[_PARAMS.index("amplitude")] *= 1.0 + noise.resp_am * np.sin(2 * np.pi * noise.resp_hz * tf + np.pi / 3)
        ibi = 60.0 / p[0] * (1.0 + noise.hr_jitter * rng.standard_normal())
        if tf + ibi >= duration - 0.3:
            break
        params.append(p)
        feet.append(tf + ibi)
    feet_arr = np.array(feet)
    params_arr = np.array(params)
    n_beats = len(params_arr)

    clean_ac = np.zeros(n)
    sys_peak_t = np.empty(n_beats)
    for k in range(n_beats):
        a, b = feet_arr[k], feet_arr[k + 1]
        i0, i1 = int(np.ceil(a * fs)), int(np.floor(b * fs))
        tt = t[i0 : i1 + 1] - a
        w = _beat_waveform(tt, params_arr[k])
        w0 = _beat_waveform(np.array([0.0]), params_arr[k])[0]
        w1 = _beat_waveform(np.array([b - a]), params_arr[k])[0]
        clean_ac[i0 : i1 + 1] = w - (w0 + (w1 - w0) * tt / (b - a))
        fine = np.linspace(0, b - a, 2000)
        wf = _beat_waveform(fine, params_arr[k]) - (w0 + (w1 - w0) * fine / (b - a))
        sys_peak_t[k] = a + fine[np.argmax(wf)]

    clean = noise.dc_level + clean_ac

    # -- noise ---------------------------------------------------------- #
    pulse_rms = np.sqrt(np.mean(clean_ac ** 2))
    noisy = clean.copy()
    if np.isfinite(noise.snr_db):
        noisy += rng.standard_normal(n) * pulse_rms / (10 ** (noise.snr_db / 20))
    if noise.wander_amp > 0:
        f1, f2 = noise.wander_hz
        for f in (f1, 0.5 * (f1 + f2), f2):
            noisy += noise.wander_amp / 3 * np.sin(2 * np.pi * f * t + rng.uniform(0, 2 * np.pi))
    artifacts: List[Tuple[float, float]] = []
    if noise.artifact_rate_per_min > 0:
        n_art = rng.poisson(noise.artifact_rate_per_min * duration / 60.0)
        for _ in range(n_art):
            d = rng.uniform(*noise.artifact_duration)
            s0 = rng.uniform(5.0, max(duration - d - 5.0, 5.0))
            sel = (t >= s0) & (t <= s0 + d)
            m = int(sel.sum())
            if m < 4:
                continue
            burst = np.cumsum(rng.standard_normal(m))
            burst = burst - np.linspace(burst[0], burst[-1], m)
            burst = burst / (np.max(np.abs(burst)) + 1e-12)
            burst += 0.5 * np.sin(2 * np.pi * rng.uniform(1.0, 3.0) * np.arange(m) / fs)
            noisy[sel] += noise.artifact_amp * burst * np.hanning(m)
            artifacts.append((float(s0), float(s0 + d)))
        artifacts.sort()

    aux: Dict[str, np.ndarray] = {}
    if with_ecg:
        ecg = 0.02 * rng.standard_normal(n)
        r_times = feet_arr[:-1] - params_arr[:, _PARAMS.index("pat")]
        for rt, ibi in zip(r_times, np.diff(feet_arr)):
            ecg += np.exp(-0.5 * ((t - rt) / 0.008) ** 2)
            ecg -= 0.2 * np.exp(-0.5 * ((t - rt - 0.025) / 0.01) ** 2)
            ecg += 0.25 * np.exp(-0.5 * ((t - rt - 0.3 * ibi) / 0.04) ** 2)
        aux["ecg"] = ecg

    beats = pd.DataFrame(params_arr, columns=list(_PARAMS))
    beats.insert(0, "t_end", feet_arr[1:])
    beats.insert(0, "t_peak", sys_peak_t)
    beats.insert(0, "t_onset", feet_arr[:-1])
    beats["in_artifact"] = [any(a - 0.2 <= tp <= b + 0.2 for a, b in artifacts) for tp in sys_peak_t]

    events, stationary, start = [], [], 0.0
    for i, seg in enumerate(segments):
        if i > 0:
            events.append(start)
        settle = 5.0 * seg.tau if i > 0 else 0.0
        if settle < seg.duration:
            stationary.append((start + settle, start + seg.duration))
        start += seg.duration

    sd = SignalData(signal=noisy, sampling_rate=fs, channel_name="synthetic_ppg", file_format="SYNTHETIC",
                    metadata={"generator": "adaptive_ppg.synthetic", "seed": seed,
                              "segments": [s.label for s in segments]},
                    aux_signals=aux)
    return SyntheticRecording(sd, clean, beats, artifacts, events, stationary, segments)


__all__ = [
    "BeatParameters",
    "NoiseConfig",
    "Segment",
    "SyntheticRecording",
    "generate",
    "rest_scenario",
    "stress_scenario",
]
