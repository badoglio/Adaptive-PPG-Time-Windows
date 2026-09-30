"""End-to-end pipeline: loading -> preprocessing -> adaptive windows -> features.

:func:`run` is the single entry point used by the CLI, the Streamlit app, the
benchmark and the tests. Its :class:`PipelineResult` carries everything
needed to reproduce an analysis (configuration, package version, git
revision) and can export all tables.
"""

from __future__ import annotations

import datetime as _dt
import io
import json
import pathlib
import subprocess
import zipfile
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Union

import numpy as np
import pandas as pd

from . import __version__
from .config import PipelineConfig
from .engine import AdaptiveEngine, AdaptiveResult
from .features import FEATURE_TO_CATEGORY, FEATURE_UNITS, FeatureResult, PPGFeatureExtractor
from .io import SignalData
from .preprocessing import PPGPreprocessor, ProcessedPPG


def git_revision() -> Optional[str]:
    """Short git hash of the package checkout (None outside a git repository)."""
    try:
        here = pathlib.Path(__file__).resolve().parent
        out = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=here, capture_output=True, text=True,
                             timeout=5)
        if out.returncode != 0:
            return None
        rev = out.stdout.strip()
        dirty = subprocess.run(["git", "status", "--porcelain", "--untracked-files=no"], cwd=here,
                               capture_output=True, text=True, timeout=5)
        return rev + ("-dirty" if dirty.stdout.strip() else "")
    except (OSError, subprocess.SubprocessError):
        return None


@dataclass
class PipelineResult:
    signal_data: SignalData
    config: PipelineConfig
    ppg: ProcessedPPG
    adaptive: AdaptiveResult
    features: FeatureResult
    provenance: Dict[str, object] = field(default_factory=dict)

    # -- tables -------------------------------------------------------- #
    def beat_table(self) -> pd.DataFrame:
        """Per-beat SQI, variability, H and window lengths."""
        q = self.ppg.beat_quality.copy()
        a = self.adaptive
        q["variability_raw"] = a.variability_raw
        q["variability_used"] = a.variability_used
        q["H"] = a.H
        q["overlap"] = a.overlap
        q["warmup"] = a.warmup_mask
        q["local_ibi"] = a.local_ibi
        for name, sched in a.schedules.items():
            q[f"T_beats_{name}"] = sched.T_beats
        return q

    def window_table(self) -> pd.DataFrame:
        """All adaptive windows of all categories, with their features and CI95."""
        frames = []
        for name, ws in self.features.adaptive.items():
            t = ws.table.copy()
            t.insert(0, "category", name)
            frames.append(t)
        return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()

    def feature_dictionary(self) -> pd.DataFrame:
        return pd.DataFrame(
            [{"feature": f, "category": c, "unit": FEATURE_UNITS.get(f, "")} for f, c in FEATURE_TO_CATEGORY.items()]
        )

    def metadata(self) -> Dict[str, object]:
        sd = self.signal_data
        return {
            **self.provenance,
            "input": {
                "channel": sd.channel_name,
                "format": sd.file_format,
                "sampling_rate_hz": sd.sampling_rate,
                "duration_s": sd.duration,
                **{k: v for k, v in sd.metadata.items() if isinstance(v, (str, int, float, bool, type(None)))},
            },
            "summary": {
                "n_beats": int(self.ppg.n_beats),
                "n_valid_beats": int(self.ppg.valid_beats_mask.sum()),
                "processing_sampling_rate_hz": float(self.ppg.sampling_rate),
                "polarity": "inverted" if self.ppg.polarity < 0 else "normal",
                "polarity_skewness": float(self.ppg.polarity_skewness),
                "effective_highcut_hz": self.ppg.effective_highcut,
                "windows": {k: len(v.windows) for k, v in self.adaptive.schedules.items()},
                "ibi_regression_slope": self.adaptive.ibi_slope,
            },
            "warnings": list(self.ppg.warnings),
            "config": self.config.to_dict(),
        }

    def tables(self) -> Dict[str, pd.DataFrame]:
        out = {
            "features_grid": self.features.grid,
            "windows": self.window_table(),
            "beats": self.beat_table(),
            "features_per_beat": self.features.per_beat,
            "feature_dictionary": self.feature_dictionary(),
        }
        for name, ws in self.features.fixed.items():
            out[f"windows_{name}"] = ws.table
        return out

    def templates_table(self, category: str) -> pd.DataFrame:
        """Ensemble templates of one category, one column per window (time from the foot, s)."""
        ws = self.features.adaptive[category]
        cols: Dict[str, np.ndarray] = {}
        n = max((len(tp.core) for tp in ws.templates), default=0)
        fs = self.ppg.sampling_rate
        cols["t_from_foot_s"] = np.arange(n) / fs
        for wid, tp in zip(ws.table["window_id"], ws.templates):
            v = np.full(n, np.nan)
            v[: len(tp.core)] = tp.core
            cols[f"w{int(wid)}"] = v
        return pd.DataFrame(cols)

    # -- export -------------------------------------------------------- #
    def to_zip(self, include_templates: bool = False) -> bytes:
        """All tables as CSV plus ``metadata.json`` in a ZIP archive."""
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
            for name, df in self.tables().items():
                zf.writestr(f"{name}.csv", df.to_csv(index=False))
            if include_templates:
                for cat in self.features.adaptive:
                    zf.writestr(f"templates_{cat}.csv", self.templates_table(cat).to_csv(index=False))
            zf.writestr("metadata.json", json.dumps(self.metadata(), indent=2, default=_json_default))
        return buf.getvalue()

    def to_excel(self) -> bytes:
        """All tables in one workbook (requires ``openpyxl``)."""
        buf = io.BytesIO()
        with pd.ExcelWriter(buf, engine="openpyxl") as xw:
            for name, df in self.tables().items():
                df.to_excel(xw, sheet_name=name[:31], index=False)
            meta = pd.DataFrame({"metadata_json": [json.dumps(self.metadata(), indent=2, default=_json_default)]})
            meta.to_excel(xw, sheet_name="metadata", index=False)
        return buf.getvalue()

    def save(self, out_dir: Union[str, pathlib.Path], include_templates: bool = False) -> List[pathlib.Path]:
        """Writes every table as CSV and ``metadata.json`` into ``out_dir``."""
        out = pathlib.Path(out_dir)
        out.mkdir(parents=True, exist_ok=True)
        paths = []
        for name, df in self.tables().items():
            p = out / f"{name}.csv"
            df.to_csv(p, index=False)
            paths.append(p)
        if include_templates:
            for cat in self.features.adaptive:
                p = out / f"templates_{cat}.csv"
                self.templates_table(cat).to_csv(p, index=False)
                paths.append(p)
        p = out / "metadata.json"
        p.write_text(json.dumps(self.metadata(), indent=2, default=_json_default), encoding="utf-8")
        paths.append(p)
        return paths


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


def run(
    signal_data: SignalData,
    config: Optional[PipelineConfig] = None,
    fixed_window_secs: Optional[Sequence[float]] = None,
    detector: str = "elgendi",
) -> PipelineResult:
    """Runs the full pipeline on one recording."""
    cfg = (config or PipelineConfig()).copy().validate()
    ppg = PPGPreprocessor(cfg.preprocessing, detector=detector).process(signal_data)
    adaptive = AdaptiveEngine(cfg.engine).process(ppg)
    features = PPGFeatureExtractor(cfg.features).process(ppg, adaptive, fixed_window_secs=fixed_window_secs)
    provenance = {
        "package": "adaptive_ppg",
        "version": __version__,
        "git_revision": git_revision(),
        "created_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
        "detector": detector,
    }
    return PipelineResult(signal_data, cfg, ppg, adaptive, features, provenance)


__all__ = ["PipelineResult", "git_revision", "run"]
