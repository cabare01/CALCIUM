#!/usr/bin/env python3
"""Step 1b: merge per-embryo red-noise results, add FDR, write rhythmic-only tables.

Run after all four run_rhythmicity.sbatch jobs (login node is fine, ~minutes).

Writes to ANALYSIS_ROOT/step1_rhythmicity/:
  rhythmicity_all_groups.csv    every tested cell + Q_global (BH within embryo)
                                + Rhythmic (Q_global < --fdr)
  rhythmicity_by_embryo.csv     per embryo: n cells, % rhythmic, % DomPeakAbove95,
                                median alpha; plus Region split when mapped
  rhythmicity_by_group.csv      the same, summarised over embryos (median, IQR)
  tables_rhythmic/<group>_CW.csv  the data/group_tables (v3) rows of rhythmic cells
                                only, same columns -> input for map_rhythmic.sbatch
  missing_cells.csv             cells in the group tables with no test result

Why Rhythmic uses P_global (decided 2026-11-04, after the full-run audit in
CALCIUM/logs/calcium_iii_rednoise_audit_2026-11-04.md):
  * P_global asks "is this cell's spectrum, relative to its AR(1) red-noise
    expectation, higher anywhere than red noise itself ever gets?" It uses the
    ratio observed / mean-AR(1), i.e. power in units of the null, so the
    statistic is a physical excess (2x, 3x) that can be compared across
    frequencies and cells.
  * P_maxT divides by the null's spread at each frequency instead. Where the
    null is very tight (the 0.08 Hz edge: many cycles averaged) a tiny,
    biologically meaningless mismatch, e.g. a white measurement-noise floor
    that AR(1) does not model, becomes a huge z. In the full run 51% of S14
    P_maxT hits sat at 0.08 Hz, nearly all from the two embryos imaged at a
    4.4-4.5 s frame interval; P_global had 0.1-2.4% edge hits per group.
  * P_global's known weakness is low power at 0.03-0.08 Hz. The real excess
    seen in the data is at 0.01-0.015 Hz, where P_global has power, and once
    P_maxT's edge hits are removed the two tests agree (S13 1.4 vs 1.8%,
    S13dev 10.6 vs 11.4%, S14 1.4 vs 1.7%, S14dev 7.2 vs 8.4%).
  * DomPeakAbove95 is not a test: it compares the already-selected highest
    peak to a pointwise 95% band, so it is inflated (25-34% of cells).
  Caveat for any report: the SD-matched AR(1) null is conservative, and more
  so for low-alpha (low-amplitude) cells, so Rhythmic=False means "not
  detectable", not "aperiodic".

Embryos in scripts/mapping/mapping_exclusions.csv are kept here on purpose
(flagged MappingExcluded) so their test results are ready once their mapping
files are traced; spatial_spectral_mapping.py drops them when mapping.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "analysis_common"))
import common  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fdr", type=float, default=0.05)
    args = ap.parse_args()
    root = common.ANALYSIS_ROOT / "step1_rhythmicity"

    parts = []
    for g in common.GROUPS:
        files = sorted((root / g).glob("*_rhythmicity.csv"))
        if not files:
            raise FileNotFoundError(f"No results for {g} in {root / g}")
        parts.extend(pd.read_csv(f) for f in files)
    res = pd.concat(parts, ignore_index=True)
    res["CellID"] = res["CellID"].astype(str)
    res["Q_global"] = res.groupby(["Group", "EmbryoID"])["P_global"].transform(lambda p: common.bh(p))
    res["Rhythmic"] = res["Q_global"] < args.fdr
    excluded_keys = tuple(pd.read_csv(common.MAPPING_EXCLUSIONS)["EmbryoKey"].str.strip())
    res["MappingExcluded"] = res["EmbryoID"].str.startswith(excluded_keys)
    print(f"Kept {res.loc[res['MappingExcluded'], 'EmbryoID'].nunique()} mapping-excluded embryos "
          f"({int(res['MappingExcluded'].sum())} cells) in all tables (flag MappingExcluded).")
    floor = 1 / (res["NSurrogates"] + 1)
    print(f"{(res['P_global'] <= floor).mean():.1%} of cells are at the Monte Carlo p floor; "
          "if this is large, increase --n-surrogates.")

    master = common.load_master()
    merged = master.merge(res[["Group", "EmbryoID", "CellID", "Rhythmic", "Q_global"]],
                          on=["Group", "EmbryoID", "CellID"], how="left", validate="1:1")
    missing = merged[merged["Rhythmic"].isna()][["Group", "EmbryoID", "CellID"]]
    missing.to_csv(root / "missing_cells.csv", index=False)
    if len(missing):
        print(f"WARNING: {len(missing)} table cells have no test result "
              f"(embryos: {missing['EmbryoID'].nunique()}). See missing_cells.csv.")

    mapped = common.load_mapped_cells()
    mapped["Region"] = common.assign_region(mapped)
    res = res.merge(mapped[["Group", "EmbryoID", "CellID", "Region"]],
                    on=["Group", "EmbryoID", "CellID"], how="left")
    res["Region"] = res["Region"].fillna("unmapped")
    res.to_csv(root / common.RHYTHM_FILE, index=False)

    def summary(d):
        return pd.Series({"NCells": len(d), "PctRhythmic": 100 * d["Rhythmic"].mean(),
                          "PctDomPeakAbove95": 100 * d["DomPeakAbove95"].mean(),
                          "MedianAlpha": d["Alpha"].median()})
    emb = res.groupby(["Group", "EmbryoID"]).apply(summary, include_groups=False)
    reg = (res[res["Region"] != "unmapped"].groupby(["Group", "EmbryoID", "Region"])
           .apply(summary, include_groups=False).unstack("Region"))
    reg.columns = [f"{a}_{b}" for a, b in reg.columns]
    emb = emb.join(reg).reset_index()
    emb.insert(2, "MappingExcluded", emb["EmbryoID"].str.startswith(excluded_keys))
    emb.to_csv(root / "rhythmicity_by_embryo.csv", index=False)
    num = [c for c in emb.columns if c not in ("Group", "EmbryoID", "MappingExcluded")]
    grp = emb.groupby("Group")[num].describe(percentiles=[.25, .5, .75])
    grp.to_csv(root / "rhythmicity_by_group.csv")
    print(emb.groupby("Group")[["PctRhythmic", "PctDomPeakAbove95"]].median().round(1).to_string())

    out = root / "tables_rhythmic"
    out.mkdir(exist_ok=True)
    for g in common.GROUPS:
        keep = merged[(merged["Group"] == g) & (merged["Rhythmic"] == True)]  # noqa: E712
        keep = keep.drop(columns=["Rhythmic", "Q_global"])
        keep.to_csv(out / f"{common.GROUP_FOLDER[g]}_CW.csv", index=False)
        print(f"{g}: {len(keep)} rhythmic cells from {keep['EmbryoID'].nunique()} embryos")


if __name__ == "__main__":
    main()
