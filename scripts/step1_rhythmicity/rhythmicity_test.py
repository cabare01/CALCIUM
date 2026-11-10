#!/usr/bin/env python3
"""Step 1: per-cell red-noise (AR(1)) significance test of the wavelet spectrum.

Question: does a cell's time-averaged wavelet spectrum contain more power than
a red-noise (AR(1)) process with the same lag-1 autocorrelation and variance
would give? If not, DominantFreqHz only locates the largest noise peak.
(Torrence & Compo 1998, Bull. Am. Meteorol. Soc. 79:61, section 4; Monte Carlo
form, because power here is |W|^2/scale with COI masking per segment, so the
closed-form chi-square levels do not apply directly.)

For every kept cell of every embryo in --group:
  1. Rebuild the exact dF/F segments with scripts/wavelets/01_preprocess_data.py
     (same raw file as the v3 run, read from run_manifest.json).
  2. Estimate AR(1) alpha from the pooled, detrended segments:
     alpha = (r1 + sqrt(r2)) / 2 (Torrence & Compo bias-reduced form), clipped
     to [0, 0.99].
  3. Simulate --n-surrogates AR(1) series per segment (same length, dt,
     linear detrend and standard deviation as the observed segment) and pass
     them through the identical CWT, scale rectification and COI masking.
  4. Per cell:
       P_global   = Monte Carlo p of the max over frequencies of
                    observed / mean-surrogate power (controls for testing 80
                    frequencies). Minimum p = 1/(n+1).
       DomPeakAbove95 = observed power at DominantFreqHz exceeds the
                    pointwise 95th surrogate percentile.
       MaxRatioFreqHz = frequency where observed/red-noise ratio is largest.
       P_maxT     = Westfall-Young max-T: log power z-scored per frequency
                    against the surrogates, max over frequencies, calibrated
                    on the surrogates' max. Weights frequencies ~equally
                    (P_global's max-ratio is dominated by the noisy
                    low-frequency edge). MaxT / MaxTFreqHz = the statistic
                    and where it peaks.
  Self-check: the observed spectrum recomputed here must equal the saved
  SpectralLibrary.npz power for the first cells; otherwise the run stops.

Outputs (per embryo, resumable; finished embryos are skipped on resubmit):
  <out>/<group>/<embryo>_rhythmicity.csv
  <out>/<group>/<embryo>_rednoise.npz   (CellID, Frequency, SurrogateMean, Surrogate95)
FDR correction and the Rhythmic flag are added by build_rhythmic_tables.py.

Run a pilot first (see PROJECT_PLAN): --embryos <one> --max-cells 50 to time it.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
import zlib
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "analysis_common"))
import common  # noqa: E402

WAVELET_DIR = common.SCRIPTS / "wavelets"
preprocess = common.load_module(WAVELET_DIR / "01_preprocess_data.py", "wl_preprocess")
spectral = common.load_module(WAVELET_DIR / "02_spectral_features.py", "wl_spectral")
N_SELF_CHECK = 3
CHUNK = 100  # surrogates per CWT call (memory ~ 80 scales x CHUNK x samples x 16 B)


def check_manifest(manifest: dict, axis: np.ndarray) -> None:
    """The modules read env vars at import; they must match the v3 run."""
    expected = {
        "wavelet": spectral.WAVELET, "n_frequencies": int(axis.size),
        "detrend": preprocess.DETREND, "trim_at_split": preprocess.TRIM_AT_SPLIT,
        "min_cell_frames": preprocess.MIN_CELL_FRAMES,
    }
    bad = {k: (manifest.get(k), v) for k, v in expected.items() if manifest.get(k) != v}
    if not np.isclose(manifest["frequency_min_hz"], axis[0]) or \
            not np.isclose(manifest["frequency_max_hz"], axis[-1]):
        bad["frequency_axis"] = ((manifest["frequency_min_hz"], manifest["frequency_max_hz"]),
                                 (axis[0], axis[-1]))
    if bad:
        raise RuntimeError(f"Pipeline settings differ from the v3 manifest (manifest, now): {bad}. "
                           "Unset the CALCIUM_* overrides or match them to the manifest.")


def batch_mean_power(batch: list[np.ndarray], dts: list[float], axis: np.ndarray) -> np.ndarray:
    """mean_power_spectrum() from 02_spectral_features.py for B series at once.

    batch[k] has shape (B, n_k) = B realisations of segment k. Returns (B, nfreq).
    Mirrors the masking and normalisation of the original line by line.
    """
    import pywt
    wavelet = pywt.ContinuousWavelet(spectral.WAVELET)
    center = pywt.central_frequency(wavelet)
    coi_seconds = spectral.coi_factor() / axis
    n_batch = batch[0].shape[0]
    power_sum = np.zeros((n_batch, axis.size))
    counts = np.zeros(axis.size, dtype=int)
    for signal, dt in zip(batch, dts):
        scales = center / (axis * dt)
        coef, _ = pywt.cwt(signal, scales, wavelet, sampling_period=dt, method="fft", axis=-1)
        power = np.abs(coef) ** 2 / scales[:, None, None]          # (scales, B, n)
        elapsed = np.arange(signal.shape[-1]) * dt
        edge = np.minimum(elapsed, elapsed[-1] - elapsed)
        valid = edge[None, :] >= coi_seconds[:, None]               # (scales, n)
        power_sum += np.einsum("sbn,sn->bs", power, valid)
        counts += valid.sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = np.where(counts > 0, power_sum / np.maximum(counts, 1), np.nan)
    return np.maximum(mean, 0.0)  # spectral_features() clips at 0 before features


def ar1_alpha(segments: list[dict]) -> tuple[float, float, float]:
    num1 = num2 = den = 0.0
    for seg in segments:
        x = np.asarray(seg["dFF"], float)
        x = x - x.mean()
        den += x @ x
        num1 += x[1:] @ x[:-1]
        num2 += x[2:] @ x[:-2]
    r1, r2 = (num1 / den, num2 / den) if den > 0 else (0.0, 0.0)
    alpha = (r1 + np.sqrt(r2)) / 2 if r2 > 0 else r1
    return float(np.clip(alpha, 0.0, 0.99)), float(r1), float(r2)


def surrogates(seg: dict, alpha: float, n: int, rng: np.random.Generator) -> np.ndarray:
    from scipy.signal import detrend, lfilter
    length = len(seg["dFF"])
    e = rng.standard_normal((n, length))
    scale = np.sqrt(1 - alpha ** 2)
    zi = (alpha * rng.standard_normal(n))[:, None]   # stationary start
    x, _ = lfilter([scale], [1, -alpha], e, axis=1, zi=zi)
    if preprocess.DETREND == "linear":
        x = detrend(x, axis=1, type="linear")
    target = np.std(seg["dFF"])
    sd = x.std(axis=1, keepdims=True)
    return x / np.where(sd > 0, sd, 1) * target


def process_embryo(job: dict) -> dict:
    group, embryo, raw_path = job["group"], job["embryo"], Path(job["raw"])
    out_csv = job["out"] / f"{embryo}_rhythmicity.csv"
    if out_csv.exists() and not job["overwrite"]:
        return {"embryo": embryo, "status": "exists"}
    t0 = time.time()
    lib = np.load(common.SPECTRAL_ROOT / common.GROUP_FOLDER[group] / embryo / "SpectralLibrary.npz")
    ids, axis, observed = lib["CellID"].astype(str), lib["Frequency"], lib["Power"]
    if job["max_cells"]:
        ids, observed = ids[:job["max_cells"]], observed[:job["max_cells"]]
    coi = spectral.coi_factor()
    table = pd.read_csv(raw_path, low_memory=False)
    try:
        frame_info = preprocess.embryo_frame_times(table)
    except Exception:
        frame_info = None
    if frame_info is not None:
        table["_time_s"] = frame_info["time"]
    by_id = {str(k): v for k, v in table.groupby("ID", sort=False, dropna=True)}
    rng = np.random.default_rng([job["seed"], zlib.crc32(embryo.encode())])
    rows, s_mean, s_95 = [], [], []
    for i, cell_id in enumerate(ids):
        cell_rows = by_id[cell_id]
        info = None if frame_info is None else {**frame_info, "time": cell_rows["_time_s"].to_numpy()}
        data = preprocess.preprocess_data(cell_rows, float(axis[0]), coi, info)
        segs = data["segments"]
        p_obs = observed[i]
        if i < N_SELF_CHECK:
            again = batch_mean_power([np.asarray(s["dFF"], float)[None, :] for s in segs],
                                     [s["dt"] for s in segs], axis)[0]
            if not np.allclose(again, np.maximum(p_obs, 0), rtol=1e-6, atol=1e-12):
                raise RuntimeError(f"{embryo} cell {cell_id}: recomputed spectrum differs from "
                                   "SpectralLibrary.npz; preprocessing is not reproducible.")
        alpha, r1, r2 = ar1_alpha(segs)
        parts = []
        for start in range(0, job["n_surr"], CHUNK):
            n = min(CHUNK, job["n_surr"] - start)
            parts.append(batch_mean_power([surrogates(s, alpha, n, rng) for s in segs],
                                          [s["dt"] for s in segs], axis))
        sur = np.vstack(parts)                                   # (n_surr, nfreq)
        mean_s = sur.mean(axis=0)
        p95 = np.quantile(sur, 0.95, axis=0)
        ratio = np.maximum(p_obs, 0) / mean_s
        r_obs = float(np.nanmax(ratio))
        r_sur = np.nanmax(sur / mean_s, axis=1)
        # Max-T (Westfall-Young): log power standardised per frequency by the
        # surrogates' mean/SD, max over frequencies, calibrated on the
        # surrogates' max. Unlike the max-ratio statistic (P_global), every
        # frequency gets ~equal weight: the ratio's null spread is ~4-5x wider
        # at 0.005 Hz than at 0.08 Hz. (Rank-based min-P was rejected: its
        # discreteness puts the p floor near 0.02-0.05 even with 1000 surrogates.)
        tiny = np.finfo(float).tiny
        log_sur = np.log(np.maximum(sur, tiny))
        mu, sd = log_sur.mean(axis=0), log_sur.std(axis=0, ddof=1)
        sd = np.where(sd > 0, sd, np.inf)
        t_obs = (np.log(np.maximum(p_obs, tiny)) - mu) / sd
        t_sur_max = ((log_sur - mu) / sd).max(axis=1)
        t_max = float(np.max(t_obs))
        dom = int(np.nanargmax(p_obs))
        rows.append({
            "Group": group, "EmbryoID": embryo, "CellID": cell_id,
            "Alpha": alpha, "R1": r1, "R2": r2, "NSegments": len(segs),
            "NSurrogates": job["n_surr"],
            "DominantFreqHz": float(axis[dom]),
            "DomPeakRatio": float(ratio[dom]),
            "DomPeakAbove95": bool(p_obs[dom] > p95[dom]),
            "MaxRatio": r_obs, "MaxRatioFreqHz": float(axis[int(np.nanargmax(ratio))]),
            "FracFreqAbove95": float(np.mean(p_obs > p95)),
            "P_global": float((1 + np.sum(r_sur >= r_obs)) / (job["n_surr"] + 1)),
            "MaxT": t_max, "MaxTFreqHz": float(axis[int(np.argmax(t_obs))]),
            "P_maxT": float((1 + np.sum(t_sur_max >= t_max)) / (job["n_surr"] + 1)),
        })
        s_mean.append(mean_s)
        s_95.append(p95)
    tmp = out_csv.with_suffix(".tmp")
    pd.DataFrame(rows).to_csv(tmp, index=False)
    np.savez_compressed(job["out"] / f"{embryo}_rednoise.npz", CellID=np.asarray(ids),
                        Frequency=axis, SurrogateMean=np.asarray(s_mean),
                        Surrogate95=np.asarray(s_95))
    tmp.rename(out_csv)  # written last so a partial embryo is never mistaken for done
    return {"embryo": embryo, "status": "done", "cells": len(rows),
            "seconds": round(time.time() - t0, 1)}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--group", choices=common.GROUPS, required=True)
    ap.add_argument("--n-surrogates", type=int, default=1000,
                    help="1000 gives a p floor of 0.001 (needed for BH-FDR to work when "
                         "only a few %% of cells are rhythmic).")
    ap.add_argument("--processes", type=int,
                    default=int(os.environ.get("SLURM_CPUS_PER_TASK", "1")))
    ap.add_argument("--embryos", default="", help="Comma-separated embryo folder names (pilot).")
    ap.add_argument("--max-cells", type=int, default=0, help="Cells per embryo (pilot); 0 = all.")
    ap.add_argument("--seed", type=int, default=20261102)
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--output-dir", type=Path,
                    default=common.ANALYSIS_ROOT / "step1_rhythmicity")
    args = ap.parse_args()

    group_dir = common.SPECTRAL_ROOT / common.GROUP_FOLDER[args.group]
    manifest = json.loads((group_dir / "run_manifest.json").read_text())
    check_manifest(manifest, spectral.frequency_axis())
    raw_by_embryo = {common.safe_name(Path(p).stem.removesuffix(" - Not_self_Merged")): p
                     for p in manifest["input_files"]}
    embryos = sorted(p.name for p in group_dir.iterdir() if (p / "SpectralLibrary.npz").exists())
    if args.embryos:
        wanted = {e.strip() for e in args.embryos.split(",") if e.strip()}
        unknown = wanted - set(embryos)
        if unknown:
            raise ValueError(f"Unknown embryos: {sorted(unknown)}")
        embryos = [e for e in embryos if e in wanted]
    out = args.output_dir / args.group
    if args.max_cells:
        out = args.output_dir / "pilot" / args.group
    out.mkdir(parents=True, exist_ok=True)
    jobs = [{"group": args.group, "embryo": e, "raw": raw_by_embryo[e], "out": out,
             "n_surr": args.n_surrogates, "max_cells": args.max_cells, "seed": args.seed,
             "overwrite": args.overwrite} for e in embryos]
    print(f"{args.group}: {len(jobs)} embryos, {args.n_surrogates} surrogates, "
          f"{args.processes} processes -> {out}", flush=True)
    with Pool(args.processes) as pool:
        for result in pool.imap_unordered(process_embryo, jobs):
            print(result, flush=True)
    (out / "settings.json").write_text(json.dumps(
        {**vars(args), "output_dir": str(args.output_dir), "manifest_run_id": manifest["run_id"]},
        indent=2, default=str))


if __name__ == "__main__":
    main()
