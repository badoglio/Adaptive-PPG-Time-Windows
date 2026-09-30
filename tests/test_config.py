import pytest

from adaptive_ppg.config import PRESET_NAMES, PipelineConfig, preset


@pytest.mark.parametrize("name", PRESET_NAMES)
def test_presets_validate_and_round_trip(name, tmp_path):
    cfg = preset(name)
    path = tmp_path / "cfg.json"
    cfg.save(path)
    assert PipelineConfig.load(path).to_dict() == cfg.to_dict()


def test_yaml_round_trip(tmp_path):
    pytest.importorskip("yaml")
    cfg = preset("acute_stress")
    path = tmp_path / "cfg.yaml"
    cfg.save(path)
    assert PipelineConfig.load(path).to_dict() == cfg.to_dict()


def test_partial_dict_keeps_defaults():
    cfg = PipelineConfig.from_dict({"engine": {"theta": 0.4, "categories": {"macro": {"T_max_sec": 12}}}})
    default = PipelineConfig()
    assert cfg.engine.theta == 0.4
    assert cfg.engine.categories["macro"].T_max_sec == 12
    assert cfg.engine.categories["derivatives"].T_crit_beats == default.engine.categories["derivatives"].T_crit_beats


@pytest.mark.parametrize("bad", [
    {"preprocessing": {"lowcut": 10.0, "highcut": 5.0}},
    {"engine": {"variability_metric": "entropy"}},
    {"engine": {"fixed_scale": [3.0, 1.0]}},
    {"engine": {"O_max": 1.0}},
    {"engine": {"categories": {"macro": {"T_min_beats": 0}}}},
    {"features": {"sg_polyorder": 1}},
])
def test_invalid_values_raise(bad):
    with pytest.raises(ValueError):
        PipelineConfig.from_dict(bad).validate()


def test_unknown_key_raises():
    with pytest.raises(ValueError):
        PipelineConfig.from_dict({"engine": {"tehta": 0.3}})
    with pytest.raises(ValueError):
        PipelineConfig.from_dict({"engine": {"categories": {"micro": {"T_min_beats": 3}}}})


def test_category_can_be_disabled():
    cfg = PipelineConfig.from_dict({"engine": {"categories": {"derivatives": None}}})
    assert set(cfg.engine.categories) == {"macro", "time_volume"}
