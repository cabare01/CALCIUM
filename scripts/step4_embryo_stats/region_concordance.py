#!/usr/bin/env python3
"""Step 4c: reconcile the cell-level (paper rule) and embryo-level region tests.

Reads <mapping-root>/combined/Region_Comparison_CellLevel.csv,
Region_Comparison_EmbryoLevel.csv and Embryo_Region_Summary.csv, adds a paired
embryo-level effect size (Cohen's dz of neural - non-neural embryo medians),
and labels every group x feature:
  Supported          paper rule met AND embryo Wilcoxon p_bonf < 0.05 AND same sign
  Cell-level only    paper rule met, embryo test not significant (pseudoreplication risk)
  Embryo-level only  embryo test significant, |d| <= 0.2 or K-S not significant
  Not supported      neither
Only "Supported" rows should be stated as findings on the poster/manuscript.

Run once on mapping_v3 (all cells) and, after step 1, on mapping_v3_rhythmic.
Output: ANALYSIS_ROOT/step4_embryo_stats/region_concordance_<label>.csv
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
    ap.add_argument("--mapping-root", type=Path, default=common.MAPPING_ROOT)
    ap.add_argument("--label", default="all_cells")
    ap.add_argument("--min-cells", type=int, default=10, help="Same as combine_group_maps.py")
    args = ap.parse_args()
    comb = args.mapping_root / "combined"
    cell = pd.read_csv(comb / "Region_Comparison_CellLevel.csv")
    emb = pd.read_csv(comb / "Region_Comparison_EmbryoLevel.csv")
    summ = pd.read_csv(comb / "Embryo_Region_Summary.csv")

    dz = []
    for (g, f), d in summ.groupby(["Group", "Feature"]):
        wide = d.pivot(index="EmbryoID", columns="Region", values="Median")
        n = d.pivot(index="EmbryoID", columns="Region", values="NCells").fillna(0)
        ok = n.reindex(columns=["neural", "non-neural"]).fillna(0).ge(args.min_cells).all(axis=1)
        diff = (wide.reindex(columns=["neural", "non-neural"])[ok].dropna()
                .pipe(lambda w: w["neural"] - w["non-neural"]))
        dz.append({"Group": g, "Feature": f,
                   "PairedDz": diff.mean() / diff.std(ddof=1) if len(diff) > 2 and diff.std(ddof=1) > 0 else np.nan,
                   "PctEmbryosNeuralHigher": 100 * (diff > 0).mean() if len(diff) else np.nan})
    t = cell.merge(emb, on=["Group", "Feature"]).merge(pd.DataFrame(dz), on=["Group", "Feature"])
    cell_sig = t["Significant_paper_rule"].astype(bool)
    emb_sig = t["Wilcoxon_p_bonferroni"] < 0.05
    same_sign = np.sign(t["CohensD"]) == np.sign(t["MedianOfDifferences"])
    t["Concordance"] = np.select(
        [cell_sig & emb_sig & same_sign, cell_sig & ~emb_sig, ~cell_sig & emb_sig],
        ["Supported", "Cell-level only", "Embryo-level only"], default="Not supported")
    out = common.out_dir("step4_embryo_stats")
    cols = ["Group", "Feature", "Concordance", "CohensD", "KS_p_bonferroni", "NEmbryosPaired",
            "MedianOfDifferences", "PairedDz", "PctEmbryosNeuralHigher", "Wilcoxon_p_bonferroni"]
    t[cols].to_csv(out / f"region_concordance_{args.label}.csv", index=False)
    print(t[cols].round(4).to_string(index=False))


if __name__ == "__main__":
    main()
