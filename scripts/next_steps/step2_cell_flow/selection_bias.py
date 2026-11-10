#!/usr/bin/env python3
"""Step 2b: are the cells that reach the maps spatially representative?

Each SuperEmbryo boundary file gives template coordinates for every cell it
lists, including cells the wavelet pipeline skipped. So within one embryo we
can compare, on the same template:
  Kept     cells with spectral results (in the group table)
  Skipped  cells the wavelet pipeline rejected (skipped_cells.csv)
Per embryo: % neural (Sox2 polygon), median template Y (anterior-posterior),
median |X - template midline| (medio-lateral). Kept minus Skipped is tested
across embryos with a Wilcoxon signed-rank test per group (embryo = unit).
A consistent difference means the maps over- or under-represent a region.

Outputs in ANALYSIS_ROOT/step2_cell_flow/:
  selection_bias_by_embryo.csv, selection_bias_tests.csv
"""
from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "analysis_common"))
import common  # noqa: E402

MIN_CELLS = 30


def boundary_index() -> dict[str, Path]:
    found = defaultdict(list)
    for p in common.BOUNDARY_ROOT.rglob("*-Boundary_SuperEmbryo.csv"):
        found[common.embryo_ref.embryo_key(p.name)].append(p)
    dup = {k: v for k, v in found.items() if len(v) > 1}
    if dup:
        print(f"WARNING: several boundary files per key (embryos skipped): {sorted(dup)}")
    return {k: v[0] for k, v in found.items() if len(v) == 1}


def main():
    out = common.out_dir("step2_cell_flow")
    excluded = set(pd.read_csv(common.MAPPING_EXCLUSIONS)["EmbryoKey"].str.strip())
    index = boundary_index()
    midline = {t: float(np.mean(common.load_template(p)[:, 0])) for t, p in common.SOX2.items()}
    rows = []
    for group in common.GROUPS:
        gdir = common.SPECTRAL_ROOT / common.GROUP_FOLDER[group]
        for edir in sorted(p for p in gdir.iterdir() if (p / "tables").is_dir()):
            key = common.embryo_ref.embryo_key(edir.name)
            if key in excluded or key not in index:
                continue
            b = pd.read_csv(index[key], usecols=lambda c: c in ("Calcium_ID", "CellID",
                                                                "Centroid_1", "Centroid_2"))
            idc = "Calcium_ID" if "Calcium_ID" in b else "CellID"
            b = b.rename(columns={idc: "CellID"}).dropna(subset=["Centroid_1", "Centroid_2"])
            b["CellID"] = b["CellID"].astype(str).str.strip()
            b = b.drop_duplicates("CellID")
            kept = set(pd.read_csv(edir / "tables/singlecell_results_final.csv",
                                   usecols=["CellID"])["CellID"].astype(str))
            skipped = set(pd.read_csv(edir / "tables/skipped_cells.csv")["CellID"].astype(str))
            b["Status"] = np.where(b.CellID.isin(kept), "Kept",
                                   np.where(b.CellID.isin(skipped), "Skipped", "NotTracked"))
            b["Group"] = group
            b["Region"] = common.assign_region(b)
            b["ML"] = (b["Centroid_1"] - midline[common.TEMPLATE_FOR_GROUP[group]]).abs()
            row = {"Group": group, "EmbryoID": edir.name}
            for status in ("Kept", "Skipped", "NotTracked"):
                s = b[b.Status == status]
                row[f"N_{status}"] = len(s)
                ok = len(s) >= MIN_CELLS
                row[f"PctNeural_{status}"] = 100 * (s.Region == "neural").mean() if ok else np.nan
                row[f"MedianY_{status}"] = s["Centroid_2"].median() if ok else np.nan
                row[f"MedianML_{status}"] = s["ML"].median() if ok else np.nan
            rows.append(row)
    emb = pd.DataFrame(rows)
    emb.to_csv(out / "selection_bias_by_embryo.csv", index=False)

    tests = []
    for group in (*common.GROUPS, "All"):
        e = emb if group == "All" else emb[emb.Group == group]
        for metric in ("PctNeural", "MedianY", "MedianML"):
            d = (e[f"{metric}_Kept"] - e[f"{metric}_Skipped"]).dropna()
            w = stats.wilcoxon(d) if len(d) >= 6 and d.abs().sum() > 0 else None
            tests.append({"Group": group, "Metric": metric, "NEmbryos": len(d),
                          "MedianKeptMinusSkipped": d.median() if len(d) else np.nan,
                          "Wilcoxon_p": w.pvalue if w else np.nan})
    t = pd.DataFrame(tests)
    t["p_holm"] = common.holm(t["Wilcoxon_p"])
    t.to_csv(out / "selection_bias_tests.csv", index=False)
    print(t.round(4).to_string())


if __name__ == "__main__":
    main()
