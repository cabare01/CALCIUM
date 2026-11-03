#!/usr/bin/env python3
"""Step 2a: how many cells survive each stage, per embryo and per group.

Stages (one row per embryo):
  Tracked          unique track IDs in the raw merged CSV (= kept + skipped)
  PassedSpectralQC kept by the wavelet pipeline (singlecell_results_final.csv)
  Skipped_<reason> wavelet-pipeline skip reasons (skipped_cells.csv)
  ExcludedEmbryo   1 if the embryo is in mapping_exclusions.csv
  Mapped           unique cells joined to the SuperEmbryo boundary file
  InGrid           mapped cells inside the legacy CNPU grid (used in maps)
  Rhythmic         mapped cells flagged rhythmic in step 1 (if available)

Outputs in ANALYSIS_ROOT/step2_cell_flow/:
  cell_flow_by_embryo.csv, cell_flow_by_group.csv, skip_reasons_by_group.csv
No HPC needed; reads only small per-embryo tables.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "analysis_common"))
import common  # noqa: E402


def reason_category(text: str) -> str:
    text = str(text)
    if text.startswith("No segment longer"):
        return "SegmentTooShort"
    if "cone-of-influence-valid" in text:
        return "TooLittleCOIData"
    if text.startswith("Fewer than two"):
        return "TooFewObservations"
    if "observed frames" in text:
        return "Under100Frames"
    return "Other:" + re.sub(r"\d+(\.\d+)?", "N", text)[:40]


def main():
    out = common.out_dir("step2_cell_flow")
    excluded = set(pd.read_csv(common.MAPPING_EXCLUSIONS)["EmbryoKey"].str.strip())
    mapped = common.load_mapped_cells()
    try:
        rhythm = common.load_rhythmicity()
    except FileNotFoundError:
        rhythm = None
        print("Step 1 output not found; Rhythmic column left empty.")

    rows, reasons = [], []
    for group in common.GROUPS:
        gdir = common.SPECTRAL_ROOT / common.GROUP_FOLDER[group]
        for edir in sorted(p for p in gdir.iterdir() if (p / "tables").is_dir()):
            embryo = edir.name
            kept = pd.read_csv(edir / "tables/singlecell_results_final.csv", usecols=["CellID"])
            skipped = pd.read_csv(edir / "tables/skipped_cells.csv")
            cats = skipped["Reason"].map(reason_category).value_counts()
            m = mapped[(mapped.Group == group) & (mapped.EmbryoID == embryo)]
            row = {"Group": group, "EmbryoID": embryo,
                   "EmbryoKey": common.embryo_ref.embryo_key(embryo),
                   "Tracked": len(kept) + len(skipped), "PassedSpectralQC": len(kept),
                   "ExcludedEmbryo": int(common.embryo_ref.embryo_key(embryo) in excluded),
                   "Mapped": m["CellID"].nunique(),
                   "InGrid": m.loc[m.InGrid, "CellID"].nunique()}
            if rhythm is not None:
                r = rhythm[(rhythm.Group == group) & (rhythm.EmbryoID == embryo) & rhythm.Rhythmic]
                row["Rhythmic"] = int(r["CellID"].isin(m["CellID"]).sum())
            row.update({f"Skipped_{k}": int(v) for k, v in cats.items()})
            rows.append(row)
            reasons.extend({"Group": group, "Reason": k, "N": int(v)} for k, v in cats.items())

    emb = pd.DataFrame(rows).fillna(0)
    for col in ("PassedSpectralQC", "Mapped", "InGrid"):
        emb[f"Pct{col}_of_Tracked"] = 100 * emb[col] / emb["Tracked"]
    emb["PctMapped_of_Passed"] = 100 * emb["Mapped"] / emb["PassedSpectralQC"]
    emb.to_csv(out / "cell_flow_by_embryo.csv", index=False)

    counts = [c for c in emb.columns if c in ("Tracked", "PassedSpectralQC", "Mapped", "InGrid", "Rhythmic")
              or c.startswith("Skipped_")]
    grp = emb.groupby("Group")[counts].sum()
    grp["NEmbryos"] = emb.groupby("Group").size()
    grp["NEmbryosMapped"] = emb[emb.Mapped > 0].groupby("Group").size()
    grp.loc["All"] = grp.sum()
    for col in ("PassedSpectralQC", "Mapped", "InGrid", "Rhythmic"):
        if col in grp:
            grp[f"Pct{col}_of_Tracked"] = (100 * grp[col] / grp["Tracked"]).round(1)
    grp.to_csv(out / "cell_flow_by_group.csv")
    (pd.DataFrame(reasons).groupby(["Group", "Reason"])["N"].sum().unstack(fill_value=0)
     .to_csv(out / "skip_reasons_by_group.csv"))
    print(grp.T.to_string())


if __name__ == "__main__":
    main()
