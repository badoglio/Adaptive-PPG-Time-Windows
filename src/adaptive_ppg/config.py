"""Single serializable configuration for the whole pipeline.

Every tunable parameter of preprocessing, the adaptive engine and feature
extraction lives in :class:`PipelineConfig`. The configuration can be written
to / read from JSON (always available) or YAML (if ``pyyaml`` is installed) and
is embedded, together with the package version and git hash, in every export.
"""

from __future__ import annotations

import copy
import json
import pathlib
from dataclasses import asdict, dataclass, field, fields, is_dataclass
from typing import Any, Dict, Optional, Tuple, Union

FEATURE_CATEGORIES_ORDER: Tuple[str, ...] = ("macro", "time_volume", "derivatives")


@dataclass
class PreprocessingConfig:
    """Filtering, beat detection, segmentation and signal-quality settings.

    Two filter branches are used: a narrow *detection* band tuned for the
    Elgendi (2013) systolic-peak detector, and a wide *morphology* band that
    preserves the systolic upstroke and the second-derivative waves.
    """

    target_fs: Optional[float] = 125.0      # every input is resampled to this rate (Hz); None keeps it
    polarity: str = "auto"                 # "auto" (sign of the derivative skewness) | "normal" | "inverted"
    lowcut: float = 0.5
    highcut: float = 30.0
    filter_order: int = 4
    detection_lowcut: float = 0.5
    detection_highcut: float = 8.0
    bypass_filter: bool = False
    onset_method: str = "tangent"          # "tangent" (max-slope tangent) | "minimum" (min between peaks)
    target_beat_length: int = 128
    resample_mode: str = "piecewise"       # "piecewise" (peak-aligned) | "uniform"
    piecewise_peak_fraction: float = 0.3
    # Signal-quality index (SQI) thresholds
    min_bpm: float = 30.0
    max_bpm: float = 220.0
    max_ibi_jump: float = 0.3              # relative deviation from the trailing median IBI
    min_template_corr: float = 0.86
    amplitude_ratio_limits: Tuple[float, float] = (1.0 / 3.0, 3.0)
    sqi_reference_beats: int = 30
    clip_min_run_sec: float = 0.04


@dataclass
class CategoryConfig:
    """Window limits of one feature category.

    ``T_min_beats`` and ``T_crit_beats`` are counted in *valid* beats;
    ``T_max_sec`` is converted to beats with the local heart rate.
    """

    T_min_beats: int = 3
    T_max_sec: float = 30.0
    T_crit_beats: int = 3


def _default_categories() -> Dict[str, CategoryConfig]:
    return {
        "macro": CategoryConfig(T_min_beats=3, T_max_sec=30.0, T_crit_beats=3),
        "time_volume": CategoryConfig(T_min_beats=3, T_max_sec=30.0, T_crit_beats=3),
        "derivatives": CategoryConfig(T_min_beats=5, T_max_sec=60.0, T_crit_beats=10),
    }


@dataclass
class EngineConfig:
    """Adaptive windowing engine settings."""

    N_past: int = 10
    # "drift_z" (template change / its stationary expectation, ~1 when stationary) |
    # "drift" | "successive" | "dispersion" (raw RMS values, need data-driven normalization)
    variability_metric: str = "drift_z"
    ibi_correction: bool = False             # regress out the IBI-driven part of the variability
    normalization: str = "fixed"             # "fixed" | "calibration" | "rolling" | "global"
    calibration_sec: float = 60.0
    rolling_beats: int = 120
    percentiles: Tuple[float, float] = (5.0, 95.0)
    fixed_scale: Tuple[float, float] = (1.0, 3.0)   # variability mapped to H = 0 and H = 1
    mapping: str = "sigmoid"                 # "sigmoid" | "linear"
    theta: float = 0.3
    sigmoid_slope: float = 10.0
    O_min: float = 0.25
    O_max: float = 0.85
    anchor: str = "trailing"                 # "trailing" (causal) | "centered"
    max_expand_per_beat: Optional[float] = 1.0
    categories: Dict[str, CategoryConfig] = field(default_factory=_default_categories)


@dataclass
class FeatureConfig:
    """Feature extraction and output-grid settings."""

    sg_window_sec: float = 0.07
    sg_polyorder: int = 3
    template_fs: Optional[float] = 500.0   # template grid (Hz); None = processing rate
    subject_height_m: Optional[float] = None
    fixed_window_sec: float = 30.0
    fixed_overlap_frac: float = 0.5
    grid_hz: float = 1.0


@dataclass
class PipelineConfig:
    """Top-level configuration container."""

    preprocessing: PreprocessingConfig = field(default_factory=PreprocessingConfig)
    engine: EngineConfig = field(default_factory=EngineConfig)
    features: FeatureConfig = field(default_factory=FeatureConfig)

    # ------------------------------------------------------------------ #
    # Validation
    # ------------------------------------------------------------------ #
    def validate(self) -> "PipelineConfig":
        """Checks parameter consistency and raises ``ValueError`` on problems."""
        p, e, f = self.preprocessing, self.engine, self.features
        _check(p.target_fs is None or p.target_fs >= 20, "preprocessing.target_fs must be None or >= 20 Hz")
        _check(p.polarity in ("auto", "normal", "inverted"),
               "preprocessing.polarity must be 'auto', 'normal' or 'inverted'")
        _check(0 < p.lowcut < p.highcut, "preprocessing: need 0 < lowcut < highcut")
        _check(0 < p.detection_lowcut < p.detection_highcut,
               "preprocessing: need 0 < detection_lowcut < detection_highcut")
        _check(p.onset_method in ("minimum", "tangent"), "preprocessing.onset_method must be 'minimum' or 'tangent'")
        _check(p.resample_mode in ("uniform", "piecewise"),
               "preprocessing.resample_mode must be 'uniform' or 'piecewise'")
        _check(0 < p.piecewise_peak_fraction < 1, "preprocessing.piecewise_peak_fraction must be in (0, 1)")
        _check(p.target_beat_length >= 16, "preprocessing.target_beat_length must be >= 16")
        _check(0 < p.min_bpm < p.max_bpm, "preprocessing: need 0 < min_bpm < max_bpm")

        _check(e.N_past >= 2, "engine.N_past must be >= 2")
        _check(e.variability_metric in ("successive", "dispersion", "drift", "drift_z"),
               "engine.variability_metric must be 'successive', 'dispersion', 'drift' or 'drift_z'")
        _check(isinstance(e.ibi_correction, bool), "engine.ibi_correction must be true or false")
        _check(e.normalization in ("calibration", "rolling", "fixed", "global"),
               "engine.normalization must be one of calibration/rolling/fixed/global")
        _check(0 <= e.percentiles[0] < e.percentiles[1] <= 100, "engine.percentiles must satisfy 0 <= lo < hi <= 100")
        _check(e.fixed_scale[0] < e.fixed_scale[1], "engine.fixed_scale must satisfy lo < hi")
        _check(e.mapping in ("sigmoid", "linear"), "engine.mapping must be 'sigmoid' or 'linear'")
        _check(e.sigmoid_slope > 0, "engine.sigmoid_slope must be > 0")
        _check(0.0 <= e.O_min <= e.O_max < 1.0, "engine: need 0 <= O_min <= O_max < 1")
        _check(e.anchor in ("trailing", "centered"), "engine.anchor must be 'trailing' or 'centered'")
        _check(e.max_expand_per_beat is None or e.max_expand_per_beat > 0,
               "engine.max_expand_per_beat must be None or > 0")
        _check(len(e.categories) > 0, "engine.categories must not be empty")
        for name, c in e.categories.items():
            _check(c.T_min_beats >= 1, f"category '{name}': T_min_beats must be >= 1")
            _check(c.T_crit_beats >= 1, f"category '{name}': T_crit_beats must be >= 1")
            _check(c.T_max_sec > 0, f"category '{name}': T_max_sec must be > 0")

        _check(f.sg_polyorder >= 2, "features.sg_polyorder must be >= 2 (second derivative)")
        _check(f.sg_window_sec > 0, "features.sg_window_sec must be > 0")
        _check(f.template_fs is None or f.template_fs >= 50, "features.template_fs must be None or >= 50 Hz")
        _check(f.fixed_window_sec > 0, "features.fixed_window_sec must be > 0")
        _check(0 <= f.fixed_overlap_frac < 1, "features.fixed_overlap_frac must be in [0, 1)")
        _check(f.grid_hz > 0, "features.grid_hz must be > 0")
        return self

    # ------------------------------------------------------------------ #
    # Serialization
    # ------------------------------------------------------------------ #
    def to_dict(self) -> Dict[str, Any]:
        return _to_jsonable(asdict(self))

    @classmethod
    def from_dict(cls, data: Optional[Dict[str, Any]]) -> "PipelineConfig":
        """Builds a config from a (possibly partial) nested dictionary.

        Missing keys keep their defaults; unknown keys raise ``ValueError`` so
        that typos in configuration files do not go unnoticed.
        """
        data = copy.deepcopy(data or {})
        _check_keys(data, {"preprocessing", "engine", "features"}, "config")
        cfg = cls()
        cfg.preprocessing = _update_dataclass(PreprocessingConfig(), data.get("preprocessing", {}), "preprocessing")
        engine_data = dict(data.get("engine", {}))
        categories_data = engine_data.pop("categories", None)
        cfg.engine = _update_dataclass(EngineConfig(), engine_data, "engine")
        if categories_data is not None:
            # Categories are merged over the defaults; ``null`` disables one.
            _check_keys(categories_data, set(FEATURE_CATEGORIES_ORDER), "engine.categories")
            for name, cat in categories_data.items():
                if cat is None:
                    cfg.engine.categories.pop(name, None)
                else:
                    base = cfg.engine.categories.get(name, CategoryConfig())
                    cfg.engine.categories[name] = _update_dataclass(base, cat, f"engine.categories.{name}")
        cfg.features = _update_dataclass(FeatureConfig(), data.get("features", {}), "features")
        return cfg.validate()

    def to_json(self, indent: int = 2) -> str:
        return json.dumps(self.to_dict(), indent=indent)

    def save(self, path: Union[str, pathlib.Path]) -> None:
        """Writes the config as JSON or YAML depending on the file extension."""
        path = pathlib.Path(path)
        if path.suffix.lower() in (".yaml", ".yml"):
            yaml = _require_yaml()
            path.write_text(yaml.safe_dump(self.to_dict(), sort_keys=False), encoding="utf-8")
        else:
            path.write_text(self.to_json(), encoding="utf-8")

    @classmethod
    def load(cls, path: Union[str, pathlib.Path]) -> "PipelineConfig":
        """Reads a JSON or YAML configuration file."""
        path = pathlib.Path(path)
        text = path.read_text(encoding="utf-8")
        if path.suffix.lower() in (".yaml", ".yml"):
            data = _require_yaml().safe_load(text)
        else:
            data = json.loads(text)
        return cls.from_dict(data)

    def copy(self) -> "PipelineConfig":
        return copy.deepcopy(self)


# ---------------------------------------------------------------------- #
# Presets
# ---------------------------------------------------------------------- #
def preset(name: str) -> PipelineConfig:
    """Returns a named parameter preset.

    * ``default``: the values of :class:`PipelineConfig`.
    * ``rest``: long windows, slow reaction; for resting-state recordings.
    * ``acute_stress``: shorter windows and a lower θ so that the window
      contracts earlier; for protocols with fast autonomic transients.
    * ``legacy``: the original method (successive differences normalized by
      the p5-p95 range of a 60 s calibration period, θ = 0.5).
    * ``wearable_64hz``: narrower morphology band and longer derivative
      smoothing for wrist wearables sampled at ~64 Hz.
    """
    cfg = PipelineConfig()
    if name == "default":
        pass
    elif name == "rest":
        cfg.engine.theta = 0.5
        cfg.engine.max_expand_per_beat = 2.0
        cfg.engine.categories["macro"].T_max_sec = 60.0
        cfg.engine.categories["time_volume"].T_max_sec = 60.0
        cfg.engine.categories["derivatives"].T_max_sec = 90.0
    elif name == "acute_stress":
        cfg.engine.theta = 0.2
        cfg.engine.max_expand_per_beat = 0.5
        cfg.engine.categories["macro"].T_max_sec = 20.0
        cfg.engine.categories["time_volume"].T_max_sec = 20.0
        cfg.engine.categories["derivatives"].T_max_sec = 40.0
    elif name == "legacy":
        cfg.engine.variability_metric = "successive"
        cfg.engine.normalization = "calibration"
        cfg.engine.fixed_scale = (0.02, 0.15)
        cfg.engine.theta = 0.5
    elif name == "wearable_64hz":
        cfg.preprocessing.highcut = 15.0
        cfg.features.sg_window_sec = 0.11
        cfg.engine.categories["derivatives"].T_crit_beats = 15
    else:
        raise ValueError(f"Unknown preset '{name}'. Available: {', '.join(PRESET_NAMES)}")
    return cfg.validate()


PRESET_NAMES: Tuple[str, ...] = ("default", "rest", "acute_stress", "wearable_64hz", "legacy")


# ---------------------------------------------------------------------- #
# Helpers
# ---------------------------------------------------------------------- #
def _check(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _check_keys(data: Dict[str, Any], allowed: set, where: str) -> None:
    unknown = set(data) - set(allowed)
    if unknown:
        raise ValueError(f"Unknown key(s) in {where}: {sorted(unknown)}. Allowed: {sorted(allowed)}")


def _update_dataclass(obj: Any, values: Dict[str, Any], where: str) -> Any:
    names = {f.name: f for f in fields(obj)}
    _check_keys(values, set(names), where)
    for key, value in values.items():
        current = getattr(obj, key)
        if isinstance(current, tuple) and isinstance(value, (list, tuple)):
            value = tuple(value)
        setattr(obj, key, value)
    return obj


def _to_jsonable(obj: Any) -> Any:
    if is_dataclass(obj):
        return _to_jsonable(asdict(obj))
    if isinstance(obj, dict):
        return {k: _to_jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_to_jsonable(v) for v in obj]
    return obj


def _require_yaml():
    try:
        import yaml  # type: ignore
    except ImportError as exc:  # pragma: no cover - depends on the environment
        raise ImportError("YAML configuration files require 'pyyaml' (pip install pyyaml).") from exc
    return yaml
