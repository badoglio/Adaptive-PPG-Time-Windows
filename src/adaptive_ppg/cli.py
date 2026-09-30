"""Command-line interface.

Examples::

    adaptive-ppg run data/*.csv --signal-column PLETH --ecg-column II --out results/
    adaptive-ppg run recording.edf --channel Pleth --config my_config.json --preset rest
    adaptive-ppg run sample_data/Sample1.CSV --target-fs 250 --out results/sample1
    adaptive-ppg benchmark --seeds 1 2 3 --out bench/
    adaptive-ppg sweep --param engine.theta 0.2 0.3 0.5 --param engine.N_past 6 10 --out sweep.csv
    adaptive-ppg init-config config.json --preset acute_stress
    adaptive-ppg dashboard            # source checkout with the 'app' extra
"""

from __future__ import annotations

import argparse
import glob
import json
import logging
import pathlib
import subprocess
import sys
from typing import List, Optional, Sequence

import pandas as pd

from . import __version__
from .config import PRESET_NAMES, PipelineConfig, preset


def _column(value: Optional[str]):
    if value is None:
        return None
    return int(value) if value.strip().lstrip("-").isdigit() else value


def _load_config(args) -> PipelineConfig:
    cfg = preset(args.preset) if getattr(args, "preset", None) else PipelineConfig()
    if getattr(args, "config", None):
        data = PipelineConfig.load(args.config).to_dict()
        base = cfg.to_dict()
        _deep_update(base, data)
        cfg = PipelineConfig.from_dict(base)
    return cfg


def _deep_update(base: dict, new: dict) -> None:
    for k, v in new.items():
        if isinstance(v, dict) and isinstance(base.get(k), dict):
            _deep_update(base[k], v)
        else:
            base[k] = v


def _expand(paths: Sequence[str]) -> List[pathlib.Path]:
    out: List[pathlib.Path] = []
    for p in paths:
        matches = sorted(glob.glob(p)) or [p]
        out.extend(pathlib.Path(m) for m in matches)
    return out


# ---------------------------------------------------------------------- #
def cmd_run(args) -> int:
    from .io import BioSignalLoader
    from .pipeline import run

    cfg = _load_config(args)
    if args.height is not None:
        cfg.features.subject_height_m = args.height
    if args.target_fs is not None:
        cfg.preprocessing.target_fs = args.target_fs if args.target_fs > 0 else None
    if args.polarity is not None:
        cfg.preprocessing.polarity = args.polarity
    cfg.validate()
    files = _expand(args.files)
    out_root = pathlib.Path(args.out)
    rows = []
    for path in files:
        logging.info("Processing %s", path)
        try:
            with open(path, "rb") as fh:
                sd = BioSignalLoader.auto_load(fh, file_name=path.name, signal_column=_column(args.signal_column),
                                               sampling_rate=args.fs, ecg_column=_column(args.ecg_column))
            res = run(sd, cfg, fixed_window_secs=args.fixed, detector=args.detector)
            out_dir = out_root / path.stem if len(files) > 1 else out_root
            res.save(out_dir, include_templates=args.templates)
            meta = res.metadata()["summary"]
            rows.append({"file": str(path), "status": "ok", "fs": sd.sampling_rate, "duration_s": sd.duration,
                         **{k: v for k, v in meta.items() if not isinstance(v, dict)},
                         **{f"windows_{k}": v for k, v in meta["windows"].items()},
                         "warnings": " | ".join(res.ppg.warnings)})
            print(f"{path}: {meta['n_valid_beats']}/{meta['n_beats']} valid beats -> {out_dir}")
        except Exception as exc:  # keep going in batch mode
            logging.exception("Failed on %s", path)
            rows.append({"file": str(path), "status": f"error: {exc}"})
            print(f"{path}: ERROR {exc}", file=sys.stderr)
    if len(files) > 1:
        out_root.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(rows).to_csv(out_root / "batch_summary.csv", index=False)
        print(f"Batch summary: {out_root / 'batch_summary.csv'}")
    return 0 if all(r["status"] == "ok" for r in rows) else 1


def cmd_benchmark(args) -> int:
    from .benchmark import benchmark, feature_scales
    from .synthetic import generate, rest_scenario, stress_scenario

    cfg = _load_config(args)
    scenario = rest_scenario() if args.scenario == "rest" else stress_scenario()
    scales = feature_scales(generate(stress_scenario(), fs=args.fs, seed=0), cfg)
    out = pathlib.Path(args.out) if args.out else None
    all_metrics = []
    for seed in args.seeds:
        rec = generate(scenario, fs=args.fs, seed=seed)
        b = benchmark(rec, cfg, fixed_secs=args.fixed, ema_secs=args.ema, scales=scales)
        b.metrics.insert(0, "seed", seed)
        all_metrics.append(b.metrics)
        print(f"\n== seed {seed} ({args.scenario})\n{b.summary.round(3).to_string(index=False)}")
    metrics = pd.concat(all_metrics, ignore_index=True)
    if len(args.seeds) > 1:
        summary = metrics.groupby(["seed", "estimator"]).agg(nrmse=("nrmse", "median"),
                                                              stationary_nsd=("stationary_nsd", "median"),
                                                              latency_s=("latency_s", "median"))
        print("\n== mean over seeds\n" + summary.groupby("estimator").mean().sort_values("nrmse").round(3).to_string())
    if out:
        out.mkdir(parents=True, exist_ok=True)
        metrics.to_csv(out / "benchmark_metrics.csv", index=False)
        (out / "config.json").write_text(cfg.to_json(), encoding="utf-8")
        print(f"\nMetrics written to {out / 'benchmark_metrics.csv'}")
    return 0


def cmd_sweep(args) -> int:
    from .benchmark import feature_scales, sweep
    from .synthetic import generate, rest_scenario, stress_scenario

    cfg = _load_config(args)
    grid = {}
    for spec in args.param:
        key, *values = spec
        grid[key] = [json.loads(v) if _is_json(v) else v for v in values]
    recs = [generate(stress_scenario(), fs=args.fs, seed=s) for s in args.seeds]
    if args.with_rest:
        recs.append(generate(rest_scenario(), fs=args.fs, seed=max(args.seeds) + 1))
    scales = feature_scales(generate(stress_scenario(), fs=args.fs, seed=0), cfg)
    df = sweep(recs, grid, cfg, scales=scales)
    keys = list(grid)
    df[keys] = df[keys].astype(str)
    summary = df.groupby(keys)[["median_nrmse", "median_stationary_nsd", "median_latency_s"]].mean()
    print(summary.sort_values("median_nrmse").round(3).to_string())
    if args.out:
        df.to_csv(args.out, index=False)
        print(f"\nResults written to {args.out}")
    return 0


def _is_json(v: str) -> bool:
    try:
        json.loads(v)
        return True
    except ValueError:
        return False


APP_PATH = pathlib.Path(__file__).resolve().parents[2] / "app" / "streamlit_app.py"


def cmd_dashboard(args) -> int:
    if not APP_PATH.exists():
        print(f"Dashboard not found at {APP_PATH}: it ships with the source checkout "
              "(install it with 'pip install -e \".[app]\"').", file=sys.stderr)
        return 1
    return subprocess.call([sys.executable, "-m", "streamlit", "run", str(APP_PATH), *args.streamlit_args])


def cmd_init_config(args) -> int:
    cfg = preset(args.preset)
    cfg.save(args.path)
    print(f"Wrote {args.path} (preset '{args.preset}')")
    return 0


# ---------------------------------------------------------------------- #
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="adaptive-ppg", description="Adaptive-window PPG morphology analysis")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    p.add_argument("-v", "--verbose", action="store_true")
    sub = p.add_subparsers(dest="command", required=True)

    def config_args(sp):
        sp.add_argument("--config", help="JSON/YAML configuration file (partial files are merged over the preset)")
        sp.add_argument("--preset", choices=PRESET_NAMES, default=None)

    r = sub.add_parser("run", help="analyse one or more recordings")
    r.add_argument("files", nargs="+", help="input files or glob patterns (.csv .tsv .txt .mat .edf .bdf)")
    r.add_argument("--signal-column", "--channel", dest="signal_column", default=None,
                   help="PPG column / channel (name or index; MAT: channel index). "
                        "Default: first non-time column, or channel 0")
    r.add_argument("--ecg-column", default=None, help="optional ECG column / channel, enables PAT")
    r.add_argument("--fs", type=float, default=None, help="sampling rate (Hz); inferred when omitted")
    r.add_argument("--fixed", type=float, nargs="*", default=[30.0], help="fixed-window baselines (s)")
    r.add_argument("--detector", choices=("elgendi", "biosppy"), default="elgendi")
    r.add_argument("--target-fs", type=float, default=None,
                   help="processing rate (Hz); 0 keeps the native rate. Default: from the configuration (125 Hz)")
    r.add_argument("--polarity", choices=("auto", "normal", "inverted"), default=None,
                   help="signal polarity. Default: from the configuration (auto)")
    r.add_argument("--height", type=float, default=None, help="subject height (m), enables the stiffness index")
    r.add_argument("--templates", action="store_true", help="also export the ensemble templates")
    r.add_argument("--out", default="results")
    config_args(r)
    r.set_defaults(func=cmd_run)

    b = sub.add_parser("benchmark", help="adaptive vs fixed windows on synthetic ground truth")
    b.add_argument("--scenario", choices=("stress", "rest"), default="stress")
    b.add_argument("--seeds", type=int, nargs="+", default=[1, 2, 3])
    b.add_argument("--fs", type=float, default=125.0)
    b.add_argument("--fixed", type=float, nargs="*", default=[10.0, 30.0, 60.0])
    b.add_argument("--ema", type=float, nargs="*", default=[10.0])
    b.add_argument("--out", default=None)
    config_args(b)
    b.set_defaults(func=cmd_benchmark)

    s = sub.add_parser("sweep", help="sensitivity analysis of configuration parameters")
    s.add_argument("--param", nargs="+", action="append", required=True, metavar=("PATH", "VALUE"),
                   help="dotted config path followed by candidate values (JSON literals), repeatable")
    s.add_argument("--seeds", type=int, nargs="+", default=[1, 2, 3])
    s.add_argument("--with-rest", action="store_true", help="add a stationary (rest) recording")
    s.add_argument("--fs", type=float, default=125.0)
    s.add_argument("--out", default=None)
    config_args(s)
    s.set_defaults(func=cmd_sweep)

    c = sub.add_parser("init-config", help="write a configuration file to edit")
    c.add_argument("path")
    c.add_argument("--preset", choices=PRESET_NAMES, default="default")
    c.set_defaults(func=cmd_init_config)

    d = sub.add_parser("dashboard", help="start the Streamlit dashboard (source checkout, 'app' extra)")
    d.add_argument("streamlit_args", nargs=argparse.REMAINDER, help="extra arguments passed to 'streamlit run'")
    d.set_defaults(func=cmd_dashboard)
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    # argparse.REMAINDER does not capture options such as "--server.port" when they come first,
    # so unknown arguments are passed to Streamlit for 'dashboard' and rejected otherwise.
    args, unknown = parser.parse_known_args(argv)
    if unknown:
        if args.command != "dashboard":
            parser.error(f"unrecognized arguments: {' '.join(unknown)}")
        args.streamlit_args = [*unknown, *args.streamlit_args]
    logging.basicConfig(level=logging.INFO if args.verbose else logging.WARNING, format="%(levelname)s %(message)s")
    return args.func(args)


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
