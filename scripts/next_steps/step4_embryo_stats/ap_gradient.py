#!/usr/bin/env python3
"""Step 4a: per-embryo spatial gradients inside the neural plate.

The shared-scale DominantFreqHz maps show high values anteriorly and low
values posteriorly in all four groups. This turns that into a test with the
embryo as the unit:
  1. For each embryo, take its neural cells (Sox2 polygon) and compute
     Spearman rho between each feature and position along --axis:
       ap  anterior-posterior; rho > 0 means higher toward ANTERIOR
       ml  distance from the template midline; rho > 0 means higher LATERALLY
     plus a Theil-Sen slope (feature units per 100 px).
     Embryos need >= --min-cells neural cells spanning >= --min-span px
     (5th-95th percentile) along the axis; a small imaging patch cannot
     estimate a gradient.
  2. Per group: one-sample Wilcoxon of rho against 0 (is there a gradient?).
  3. Across groups: stage x membrane permutation model on rho (does the
     gradient change with stage or devitellination?).
  Holm correction across the six features within each test family.

VERIFY ORIENTATION FIRST: --anterior top assumes anterior is at small
template Y (top of the maps). Confirm against the Sox2 template / PDT paper.

Outputs in ANALYSIS_ROOT/step4_embryo_stats/<cells>/:
  gradient_<axis>_by_embryo.csv, gradient_<axis>_group_tests.csv,
  gradient_<axis>_stage_membrane.csv
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "analysis_common"))
import common  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--axis", choices=("ap", "ml"), default="ap")
    ap.add_argument("--anterior", choices=("top", "bottom"), default="top")
    ap.add_argument("--cells", choices=("all", "rhythmic"), default="all")
    ap.add_argument("--min-cells", type=int, default=100)
    ap.add_argument("--min-span", type=float, default=200.0)
    ap.add_argument("--n-perm", type=int, default=10000)
    args = ap.parse_args()
    out = common.out_dir("step4_embryo_stats", args.cells)

    cells = common.filter_cells(common.load_mapped_cells(), args.cells)
    cells["Region"] = common.assign_region(cells)
    cells = cells[cells.Region == "neural"].copy()
    if args.axis == "ap":
        sign = -1 if args.anterior == "top" else 1
        cells["Pos"] = sign * cells["Centroid_2"]
    else:
        mid = {t: float(np.mean(common.load_template(p)[:, 0])) for t, p in common.SOX2.items()}
        cells["Pos"] = (cells["Centroid_1"] - cells["Group"].map(common.TEMPLATE_FOR_GROUP).map(mid)).abs()

    rows = []
    for (g, e), d in cells.groupby(["Group", "EmbryoID"]):
        span = np.subtract(*np.percentile(d["Pos"], [95, 5]))
        base = {"Group": g, "EmbryoID": e, "NNeuralCells": len(d), "Span_px": span,
                "Stage": common.GROUP_STAGE[g], "Membrane": common.GROUP_MEMBRANE[g]}
        usable = len(d) >= args.min_cells and span >= args.min_span
        for f in common.FEATURES:
            v = d[[f, "Pos"]].dropna()
            rho = slope = np.nan
            if usable and v[f].nunique() > 1:
                rho = stats.spearmanr(v["Pos"], v[f]).statistic
                slope = stats.theilslopes(v[f], v["Pos"] / 100.0).slope
            rows.append({**base, "Feature": f, "Usable": usable, "Rho": rho, "TheilSenPer100px": slope})
    emb = pd.DataFrame(rows)
    emb.to_csv(out / f"gradient_{args.axis}_by_embryo.csv", index=False)

    tests = []
    for g in common.GROUPS:
        for f in common.FEATURES:
            r = emb[(emb.Group == g) & (emb.Feature == f)]["Rho"].dropna()
            w = stats.wilcoxon(r) if len(r) >= 6 else None
            tests.append({"Group": g, "Feature": f, "NEmbryos": len(r), "MedianRho": r.median(),
                          "PctEmbryosPositive": 100 * (r > 0).mean() if len(r) else np.nan,
                          "Wilcoxon_p": w.pvalue if w else np.nan})
    t = pd.DataFrame(tests)
    t["p_holm"] = t.groupby("Group")["Wilcoxon_p"].transform(lambda p: common.holm(p))
    t.to_csv(out / f"gradient_{args.axis}_group_tests.csv", index=False)
    print(t.round(4).to_string(index=False))

    mods = []
    for f in common.FEATURES:
        e = emb[(emb.Feature == f) & emb.Rho.notna()]
        for r in common.two_by_two(e["Rho"], e["Stage"], e["Membrane"], n_perm=args.n_perm):
            mods.append({"Feature": f, "Outcome": f"rho_{args.axis}", **r})
    m = pd.DataFrame(mods)
    m["p_holm"] = m.groupby("Term")["p_perm"].transform(lambda p: common.holm(p))
    m.to_csv(out / f"gradient_{args.axis}_stage_membrane.csv", index=False)
    print(m.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
