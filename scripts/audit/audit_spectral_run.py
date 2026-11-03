#!/usr/bin/env python3
"""Audit one spectral pipeline output root (e.g. data10/HPC_spectral_v3).

Usage: python audit_spectral_run.py <spectral_root> [<compare_root>]

Checks, per group and embryo:
  * reference-table membership: every reference embryo has exactly one output
    folder, no extras, no embryo in two groups
  * provenance: every table's RunID matches the group's run_manifest.json and
    every SourceFile is one of the manifest's input files
  * NPZ <-> CSV consistency, non-finite/negative power, NaN features,
    DominantFreqHz == argmax(Power)
  * frequency axis vs. the wavelet's spread and each embryo's Nyquist limit
  * band-edge pile-up (dominant frequency in the lowest/highest 3 bins)
  * template-mappable cells (in *-Boundary_SuperEmbryo.csv) kept vs. skipped,
    honouring scripts/mapping/mapping_exclusions.csv
  * optional comparison of kept/mappable counts with an earlier run
Writes a text report to stdout and per-embryo CSV next to the report path
given by CALCIUM_AUDIT_CSV (optional). Never modifies its inputs.
"""
import json
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import embryo_ref  # noqa: E402

BOUNDARY_ROOT = Path("/sciclone/home/cabare01/data10/boundary_info/Boundaries")
EXCLUSIONS = Path(__file__).resolve().parents[1] / "mapping" / "mapping_exclusions.csv"
GROUPS = ("stage_13", "stage_13_devitellinated", "stage_14", "stage_14_devitellinated")
WAVELET_SPREAD = 1 + 2 / (np.pi * np.sqrt(2 * 1.5))  # cmor1.5: f + 2 sigma ~= 1.37 f


def boundary_ids(key):
    path = BOUNDARY_ROOT / f"{key}-Boundary_SuperEmbryo.csv"
    if not path.exists():
        return None
    return set(pd.read_csv(path, usecols=["Calcium_ID"]).Calcium_ID.astype(str))


def audit_root(root, excluded, issues, label):
    rows, seen = [], {}
    for group in GROUPS:
        gdir = root / group
        if not gdir.is_dir():
            issues.append(f"[{label}] {group}: missing group folder")
            continue
        manifest = json.loads((gdir / "run_manifest.json").read_text())
        run_id = manifest["run_id"]
        inputs = {Path(p).name for p in manifest["input_files"]}
        wanted = embryo_ref.group_keys(group)
        dirs = sorted(d for d in gdir.iterdir() if d.is_dir())
        keys = [embryo_ref.embryo_key(d.name) for d in dirs]
        for k in sorted(wanted):
            if keys.count(k) != 1:
                issues.append(f"[{label}] {group} {k}: {keys.count(k)} output folders")
        for d, k in zip(dirs, keys):
            seen.setdefault(k, []).append(group)
            if k not in wanted:
                issues.append(f"[{label}] {group} {d.name}: not in reference column")
            sc, sk, npz = d / "tables/singlecell_results_final.csv", d / "tables/skipped_cells.csv", d / "SpectralLibrary.npz"
            if not (sc.exists() and sk.exists() and npz.exists()):
                issues.append(f"[{label}] {group} {d.name}: missing output file(s)")
                continue
            df, skipped = pd.read_csv(sc), pd.read_csv(sk)
            z = np.load(npz)
            F, P = z["Frequency"], z["Power"]
            r = {"run": label, "group": group, "embryo": d.name, "key": k,
                 "kept": len(df), "skipped": len(skipped)}
            if set(df.RunID.astype(str)) - {run_id}:
                issues.append(f"[{label}] {group} {k}: RunID differs from manifest (stale output?)")
            if set(df.SourceFile) - inputs:
                issues.append(f"[{label}] {group} {k}: SourceFile not in manifest inputs")
            if list(z["CellID"].astype(str)) != list(df.CellID.astype(str)):
                issues.append(f"[{label}] {group} {k}: NPZ/CSV CellID order mismatch")
            if len(df) and (~np.isfinite(P)).any():
                issues.append(f"[{label}] {group} {k}: non-finite power")
            if len(df) and (P < 0).any():
                issues.append(f"[{label}] {group} {k}: negative power")
            nan_rows = int(df[["DominantFreqHz", "SpectralCentroid", "TotalBandPower"]].isna().any(axis=1).sum())
            if nan_rows:
                issues.append(f"[{label}] {group} {k}: {nan_rows} result rows with NaN features")
            if len(df):
                mism = int((~np.isclose(F[np.nanargmax(P, 1)], df.DominantFreqHz.values)).sum())
                if mism:
                    issues.append(f"[{label}] {group} {k}: {mism} DominantFreqHz != argmax(Power)")
                dt = float(df.FrameInterval_s.median())
                nyq = 1 / (2 * dt)
                r.update(frame_dt=dt, nyquist=nyq,
                         top_reach_over_nyq=F[-1] * WAVELET_SPREAD / nyq,
                         edge_low=float((df.DominantFreqHz <= F[2] + 1e-12).mean()),
                         edge_high=float((df.DominantFreqHz >= F[-3] - 1e-12).mean()),
                         domf_median=float(df.DominantFreqHz.median()),
                         used_min_median=float(df.UsedDuration_min.median()))
                if r["top_reach_over_nyq"] > 1:
                    issues.append(f"[{label}] {group} {k}: f_max x wavelet spread exceeds Nyquist "
                                  f"({F[-1]:.3f} Hz, dt {dt:.2f} s)")
            r["fmin"], r["fmax"], r["nfreq"] = float(F[0]), float(F[-1]), int(F.size)
            b = boundary_ids(k)
            r["excluded_from_mapping"] = k in excluded
            if b is None:
                issues.append(f"[{label}] {group} {k}: no boundary file")
            else:
                kept_ids = set(df.CellID.astype(str))
                skip_ids = set(skipped.CellID.astype(str))
                r["mappable"] = len(b)
                r["mappable_kept"] = len(b & kept_ids)
                r["mappable_skipped"] = len(b & skip_ids)
                r["mappable_absent"] = len(b - kept_ids - skip_ids)
            reasons = skipped.Reason.astype(str).str.replace(r"\d+(\.\d+)?", "N", regex=True)
            for reason, n in reasons.value_counts().items():
                r[f"skip: {reason[:60]}"] = int(n)
            rows.append(r)
    for k, groups in seen.items():
        if len(groups) > 1:
            issues.append(f"[{label}] {k}: present in several groups {groups}")
    return pd.DataFrame(rows)


def main():
    root = Path(sys.argv[1])
    compare = Path(sys.argv[2]) if len(sys.argv) > 2 else None
    excluded = set()
    if EXCLUSIONS.exists():
        excluded = set(pd.read_csv(EXCLUSIONS).EmbryoKey.astype(str).str.strip())
    issues = []
    cur = audit_root(root, excluded, issues, root.name)
    pd.set_option("display.width", 200)
    pd.set_option("display.max_columns", 30)
    print(f"# Audit of {root}\n")
    print(f"Mapping exclusions ({EXCLUSIONS.name}): {sorted(excluded)}\n")
    ax = cur[["fmin", "fmax", "nfreq"]].drop_duplicates()
    print("Frequency axis:", ax.to_dict("records"))
    print(f"Frame interval range: {cur.frame_dt.min():.2f}-{cur.frame_dt.max():.2f} s; "
          f"max (f_max x wavelet spread)/Nyquist = {cur.top_reach_over_nyq.max():.2f} (must be < 1)\n")

    print(f"## Issues ({len(issues)})")
    for line in issues or ["none"]:
        print(" ", line)

    m = cur[~cur.excluded_from_mapping]
    summary = m.groupby("group").agg(
        embryos=("embryo", "size"), kept=("kept", "sum"), skipped=("skipped", "sum"),
        mappable=("mappable", "sum"), mappable_kept=("mappable_kept", "sum"),
        mappable_skipped=("mappable_skipped", "sum"), mappable_absent=("mappable_absent", "sum"),
        edge_low_median=("edge_low", "median"), edge_high_median=("edge_high", "median"),
        domf_median=("domf_median", "median"))
    summary["mappable_kept_pct"] = 100 * summary.mappable_kept / summary.mappable
    print("\n## Per group (mapping-excluded embryos left out)")
    print(summary.round(3).to_string())
    tot = summary[["mappable", "mappable_kept"]].sum()
    print(f"\nTotal template-mappable cells kept: {int(tot.mappable_kept):,} of {int(tot.mappable):,} "
          f"({100 * tot.mappable_kept / tot.mappable:.1f}%); manuscript n = 201,409 (179 embryos)")

    print("\n## Embryos with >50% of dominant frequencies in the lowest 3 bins")
    hi = m[m.edge_low > 0.5][["group", "key", "edge_low", "mappable_kept", "used_min_median"]]
    print(hi.round(3).to_string(index=False) if len(hi) else "  none")

    print("\n## Lowest template-mappable retention (bottom 10 embryos)")
    m2 = m.assign(pct=100 * m.mappable_kept / m.mappable)
    print(m2.nsmallest(10, "pct")[["group", "key", "mappable", "mappable_kept", "pct"]].round(1).to_string(index=False))

    skip_cols = [c for c in cur.columns if c.startswith("skip: ")]
    print("\n## Skip reasons, all cells (per group)")
    print(cur.groupby("group")[skip_cols].sum().T.to_string())

    if compare is not None and compare.is_dir():
        prev_issues = []
        prev = audit_root(compare, excluded, prev_issues, compare.name)
        both = m.merge(prev[~prev.excluded_from_mapping][["group", "key", "kept", "mappable_kept"]],
                       on=["group", "key"], suffixes=("", "_prev"))
        cmp = both.groupby("group")[["kept", "kept_prev", "mappable_kept", "mappable_kept_prev"]].sum()
        cmp["mappable_change_pct"] = 100 * (cmp.mappable_kept / cmp.mappable_kept_prev - 1)
        print(f"\n## Compared with {compare}")
        print(cmp.round(1).to_string())

    out_csv = os.environ.get("CALCIUM_AUDIT_CSV")
    if out_csv:
        cur.to_csv(out_csv, index=False)
        print(f"\nPer-embryo table: {out_csv}")


if __name__ == "__main__":
    main()
