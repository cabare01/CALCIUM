#!/usr/bin/env python3
"""Run the three wavelet-analysis stages for every cell in each input CSV (v2).

Required environment variables:
  CALCIUM_RAW_TRACE_ROOT   directory of per-embryo "* - Not_self_Merged.csv" files
                           (searched recursively when CALCIUM_GROUP is set)
  CALCIUM_OUTPUT_ROOT      new (non-existing) output directory
Optional:
  CALCIUM_GROUP            stage_13 | stage_13_devitellinated | stage_14 |
                           stage_14_devitellinated. Embryos are selected from
                           CALCIUM_REF_TABLE (default data10/Embryo_Ref_Table.csv),
                           the definitive grouping, regardless of folder or name.
  CALCIUM_RAW_TRACE_PATTERN  glob for input files (default *.csv)
  CALCIUM_FMIN / CALCIUM_FMAX / CALCIUM_NFREQ   axis (default 0.005 / 0.08 / 80)
  CALCIUM_MIN_VALID_CYCLES   COI-valid cycles required at every frequency (default 1)
  CALCIUM_DETREND            linear | none, per-segment dF/F detrend (default linear)
  CALCIUM_PYWT_WAVELET       (default cmor1.5-1.0)
  CALCIUM_MAX_MISSING_FRAMES (default 3)    CALCIUM_TRIM_AT_SPLIT (default 1)
  CALCIUM_BASELINE_FRACTION  (default 0.10) CALCIUM_GAP_SPLIT_SECONDS (fallback, 10)
  CALCIUM_MIN_CELL_FRAMES    (default 100)
  CALCIUM_MAKE_TRACE_PLOTS   true/false (default true)

Per embryo the pipeline writes:
  tables/singlecell_results_final.csv
  tables/skipped_cells.csv
  SpectralLibrary.npz
and, if enabled, one dF/F PNG per cell.

A run_manifest.json at the output root records every parameter, the package
versions and the git commit.
"""
from __future__ import annotations

import importlib.util
import json
import os
import platform
import re
import subprocess
from collections import Counter
from datetime import datetime
from pathlib import Path

import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import embryo_ref  # noqa: E402


def _dist_version(dist: str) -> str:
    """Installed distribution version. PyWavelets 1.9.0's wheel reports
    pywt.__version__ == '1.8.0' (upstream labelling bug), so don't use that."""
    import importlib.metadata
    return importlib.metadata.version(dist)


def _load_stage(filename: str, module_name: str):
    path = Path(__file__).with_name(filename)
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load pipeline stage {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_stage("01_preprocess_data.py", "wavelets_preprocess")
spectral = _load_stage("02_spectral_features.py", "wavelets_spectral")
results = _load_stage("03_build_results.py", "wavelets_results")


def _safe_name(value: object) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]", "_", str(value)).strip("._") or "unnamed"


def _git_commit() -> str:
    try:
        root = Path(__file__).resolve().parents[2]
        sha = subprocess.run(["git", "-C", str(root), "rev-parse", "HEAD"],
                             capture_output=True, text=True, check=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", str(root), "status", "--porcelain", "scripts/wavelets"],
                               capture_output=True, text=True, check=True).stdout.strip()
        return sha + (" (scripts/wavelets has uncommitted changes)" if dirty else "")
    except Exception as exc:  # git missing or not a repository
        return f"unavailable: {exc}"


def _short_reason(message: str) -> str:
    return re.sub(r"\d+(\.\d+)?", "N", message)[:80]


def main() -> None:
    source_raw = os.environ.get("CALCIUM_RAW_TRACE_ROOT")
    output_raw = os.environ.get("CALCIUM_OUTPUT_ROOT")
    if not source_raw or not output_raw:
        raise ValueError("Set CALCIUM_RAW_TRACE_ROOT and CALCIUM_OUTPUT_ROOT")
    source, output = Path(source_raw).expanduser(), Path(output_raw).expanduser()
    if not source.is_dir():
        raise NotADirectoryError(source)
    if output.exists():
        raise FileExistsError(f"Refusing to reuse existing output directory: {output}")
    pattern = os.environ.get("CALCIUM_RAW_TRACE_PATTERN", "*.csv")
    group = os.environ.get("CALCIUM_GROUP")
    missing_ref, ignored = [], []
    if group:
        wanted = embryo_ref.group_keys(group)
        by_key: dict[str, list[Path]] = {}
        for path in sorted(source.rglob(pattern)):
            by_key.setdefault(embryo_ref.embryo_key(path.name), []).append(path)
        duplicates = {k: v for k, v in by_key.items() if k in wanted and len(v) > 1}
        if duplicates:
            raise ValueError(f"Several raw files for one reference embryo: {duplicates}")
        files = sorted(by_key[k][0] for k in wanted if k in by_key)
        missing_ref = sorted(wanted - set(by_key))
        ignored = sorted(str(p.relative_to(source)) for k, v in by_key.items()
                         if k not in wanted for p in v
                         if p.parent.name == group)  # in this group's folder but not in ref
    else:
        files = sorted(source.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No input files found in {source}")
    make_plots = os.environ.get("CALCIUM_MAKE_TRACE_PLOTS", "true").lower() in {"1", "true", "yes"}
    output.mkdir(parents=True)

    import pywt  # noqa: F401  fail early if PyWavelets is missing
    axis = spectral.frequency_axis()
    coi = spectral.coi_factor()
    run_id = datetime.now().strftime("%Y-%m-%d_%H%M%S")
    manifest = {
        "run_id": run_id, "pipeline_version": "v3",
        "raw_trace_root": str(source), "pattern": pattern, "n_input_files": len(files),
        "group": group, "reference_table": str(embryo_ref.REF_TABLE) if group else None,
        "reference_embryos_without_raw_file": missing_ref,
        "files_in_group_folder_not_in_reference": ignored,
        "input_files": [str(p) for p in files],
        "wavelet": spectral.WAVELET, "frequency_min_hz": float(axis[0]),
        "frequency_max_hz": float(axis[-1]), "n_frequencies": int(axis.size),
        "power": "|W|^2 / scale (Liu et al. 2007), cone-of-influence masked, per-segment CWT",
        "coi_seconds_times_hz": coi,
        "min_segment_seconds": preprocess.min_segment_seconds(float(axis[0]), coi),
        "min_valid_cycles_per_frequency": spectral.MIN_VALID_CYCLES,
        "detrend": preprocess.DETREND,
        "gap_rule": f"bridge <= {preprocess.MAX_MISSING_FRAMES} missing frames (ND.T); split at "
                    "refocus pauses; fallback without ND.T: split > "
                    f"{preprocess.GAP_SPLIT_SECONDS:g} s",
        "trim_at_split": preprocess.TRIM_AT_SPLIT,
        "min_cell_frames": preprocess.MIN_CELL_FRAMES,
        "baseline": f"per-segment mean of lowest {preprocess.BASELINE_FRACTION:.0%} of intensities "
                    "(paper baseline (a)); dF/F=(F-F0)/F0",
        "time_source": "embryo-wide ND.T->time linear map within pause-free runs "
                       "(falls back to parsed Time column if ND.T missing)",
        "make_trace_plots": make_plots,
        "git_commit": _git_commit(),
        "versions": {"python": platform.python_version(), "numpy": np.__version__,
                     "pandas": pd.__version__, "pywt": _dist_version("PyWavelets")},
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
    }
    (output / "run_manifest.json").write_text(json.dumps(manifest, indent=2))
    if missing_ref:
        print(f"WARNING: reference embryos with no raw file (not processed): {missing_ref}", flush=True)
    if ignored:
        print(f"NOTE: files in the {group} folder that the reference table does not list "
              f"(not processed): {ignored}", flush=True)
    print(json.dumps(manifest, indent=2), flush=True)

    if make_plots:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

    totals = Counter()
    for file_number, path in enumerate(files, 1):
        embryo_name = _safe_name(path.stem.removesuffix(" - Not_self_Merged"))
        embryo_dir = output / embryo_name
        tables_dir, figures_dir = embryo_dir / "tables", embryo_dir / "figures"
        print(f"Processing file {file_number}/{len(files)}: {path.name}", flush=True)
        try:
            table = pd.read_csv(path, low_memory=False)
        except Exception as exc:
            print(f"SKIPPED FILE {path.name}: read failed: {exc}", flush=True)
            continue
        if table.empty or "ID" not in table.columns:
            print(f"SKIPPED FILE {path.name}: empty table or missing ID column", flush=True)
            continue
        tables_dir.mkdir(parents=True)
        if make_plots:
            figures_dir.mkdir()
        try:
            frame_info = preprocess.embryo_frame_times(table)
        except Exception as exc:
            print(f"  frame-time map failed ({exc}); using time-based gap rule", flush=True)
            frame_info = None
        if frame_info is not None:
            table["_time_s"] = frame_info["time"]
            print(f"  frame interval {frame_info['frame_dt']:.3f} s; "
                  f"{len(frame_info['pause_frames'])} refocus pause(s)", flush=True)
        rows, skipped, cell_ids, powers, counts = [], [], [], [], []
        for cell_id, cell_rows in table.groupby("ID", sort=False, dropna=True):
            try:
                cell_info = (None if frame_info is None else
                             {**frame_info, "time": cell_rows["_time_s"].to_numpy()})
                data = preprocess.preprocess_data(cell_rows, float(axis[0]), coi, cell_info)
                features, spectrum = spectral.spectral_features(data, axis)
            except Exception as exc:
                skipped.append({"CellID": str(cell_id), "Reason": str(exc)})
                continue
            rows.append(results.build_result_row(data, features, cell_id, path.name, run_id))
            cell_ids.append(str(cell_id))
            powers.append(spectrum["MeanPower"])
            counts.append(spectrum["ValidSamples"])
            if make_plots:
                fig, ax = plt.subplots(figsize=(9, 3))
                for segment in data["segments"]:
                    ax.plot(segment["time"], segment["dFF"], linewidth=1.0)
                ax.set(xlabel="Time (s)", ylabel="ΔF/F",
                       title=f"{path.name} | Cell {cell_id} | {len(data['segments'])} segment(s)")
                ax.grid(True)
                fig.tight_layout()
                fig.savefig(figures_dir / f"{embryo_name}_Cell_{_safe_name(cell_id)}.png", dpi=120)
                plt.close(fig)
        pd.DataFrame.from_records(rows, columns=results.RESULT_COLUMNS).to_csv(
            tables_dir / "singlecell_results_final.csv", index=False)
        pd.DataFrame.from_records(skipped, columns=["CellID", "Reason"]).to_csv(
            tables_dir / "skipped_cells.csv", index=False)
        np.savez_compressed(
            embryo_dir / "SpectralLibrary.npz", CellID=np.asarray(cell_ids), Frequency=axis,
            Power=np.asarray(powers, dtype=float).reshape((-1, axis.size)),
            ValidSamples=np.asarray(counts, dtype=int).reshape((-1, axis.size)))
        reasons = Counter(_short_reason(s["Reason"]) for s in skipped)
        totals.update({"kept": len(rows), "skipped": len(skipped)})
        totals.update(reasons)
        print(f"Finished {embryo_name}: {len(rows)} usable cells, {len(skipped)} skipped "
              f"{dict(reasons)}", flush=True)
    print(f"All files complete. Totals: {dict(totals)}", flush=True)


if __name__ == "__main__":
    main()
