#!/usr/bin/env python3
"""Run the three wavelet-analysis stages for every cell in each input CSV.

Set CALCIUM_RAW_TRACE_ROOT to the directory containing per-embryo CSV files
and CALCIUM_OUTPUT_ROOT to a new output directory before submitting the
companion SLURM script. CALCIUM_RAW_TRACE_PATTERN defaults to *.csv.
"""
from __future__ import annotations

import importlib.util
import os
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd


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
    files = sorted(source.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No files matching {pattern} in {source}")
    output.mkdir(parents=True)
    frequency_axis = np.linspace(0.02, 0.20, 80)
    import pywt
    print(f"PyWavelets {pywt.__version__}; wavelet={spectral.WAVELET}; files={len(files)}", flush=True)

    for file_number, path in enumerate(files, 1):
        embryo_name = _safe_name(path.stem.removesuffix(" - Not_self_Merged"))
        embryo_dir = output / embryo_name
        tables_dir = embryo_dir / "tables"
        figures_dir = embryo_dir / "figures"
        tables_dir.mkdir(parents=True)
        figures_dir.mkdir()
        print(f"Processing file {file_number}/{len(files)}: {path.name}", flush=True)
        try:
            table = pd.read_csv(path)
        except Exception as exc:
            warnings.warn(f"Skipping file {path.name} (read failed): {exc}")
            tables_dir.rmdir()
            figures_dir.rmdir()
            embryo_dir.rmdir()
            continue
        if table.empty or "ID" not in table.columns:
            warnings.warn(f"Skipping {path.name}: empty table or missing ID")
            continue
        rows, cell_ids, powers = [], [], []
        for cell_id, group in table.groupby("ID", sort=False, dropna=True):
            try:
                data = preprocess.preprocess_data(group)
                spectral_feat, spectrum = spectral.spectral_features(data, frequency_axis)
                rows.append(results.build_result_row(data, spectral_feat, cell_id, path.name))
                cell_ids.append(str(cell_id))
                powers.append(spectrum["MeanPower"])
                import matplotlib.pyplot as plt
                fig, ax = plt.subplots(figsize=(9, 3))
                ax.plot(data["time"], data["dFF"], linewidth=1.2)
                ax.set(xlabel="Time (s)", ylabel="ΔF/F",
                       title=f"{path.name} | Cell {cell_id} | Preprocessed ΔF/F")
                ax.grid(True)
                fig.tight_layout()
                fig.savefig(figures_dir / f"{embryo_name}_Cell_{_safe_name(cell_id)}.png", dpi=200)
                plt.close(fig)
            except Exception as exc:
                warnings.warn(f"Skipping embryo {embryo_name} cell {cell_id}: {exc}")
        columns = ["CellID", "RunID", "SourceFile", "DominantFreqHz", "SpectralCentroid",
                   "SpectralSpread", "SpectralEntropy", "TotalBandPower", "PeakPowerFraction",
                   "Duration_min", "MeanX", "MeanY"]
        pd.DataFrame.from_records(rows, columns=columns).to_csv(
            tables_dir / "singlecell_results_final.csv", index=False)
        np.savez_compressed(embryo_dir / "SpectralLibrary.npz",
                            CellID=np.asarray(cell_ids), Frequency=frequency_axis,
                            Power=np.asarray(powers, dtype=float).reshape((-1, frequency_axis.size)))
        print(f"Finished {embryo_name}: {len(rows)} usable cells", flush=True)
    print("All files complete.", flush=True)


if __name__ == "__main__":
    main()
