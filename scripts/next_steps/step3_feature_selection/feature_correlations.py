#!/usr/bin/env python3
"""Step 3: which of the six spectral features carry independent information?

Two views, because they answer different questions:
  within-embryo   Spearman rho between features across cells, computed in each
                  embryo separately, then the median over embryos (cell-to-cell
                  redundancy, free of embryo-level confounds)
  between-embryo  Spearman rho between embryo medians (what the group
                  comparisons in step 4 actually test)
A pair with |rho| >= --redundant in both views is redundant: report one of
the two. Run with --cells rhythmic after step 1 to repeat on rhythmic cells.

Outputs in ANALYSIS_ROOT/step3_feature_selection/<cells>/:
  rho_within_embryo_median.csv, rho_between_embryo.csv, redundant_pairs.csv,
  feature_correlations.png/.pdf (two heatmaps, diverging -1..1)
"""
from __future__ import annotations

import argparse
import itertools
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "analysis_common"))
import common  # noqa: E402

MIN_CELLS = 50


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cells", choices=("all", "rhythmic"), default="all")
    ap.add_argument("--redundant", type=float, default=0.7)
    args = ap.parse_args()
    out = common.out_dir("step3_feature_selection", args.cells)
    feats = list(common.FEATURES)
    cells = common.filter_cells(common.load_mapped_cells(), args.cells)

    per_embryo = []
    for (g, e), d in cells.groupby(["Group", "EmbryoID"]):
        if len(d) >= MIN_CELLS:
            per_embryo.append(d[feats].corr(method="spearman").to_numpy())
    within = pd.DataFrame(np.nanmedian(np.stack(per_embryo), axis=0), index=feats, columns=feats)
    q1 = pd.DataFrame(np.nanpercentile(np.stack(per_embryo), 25, axis=0), index=feats, columns=feats)
    q3 = pd.DataFrame(np.nanpercentile(np.stack(per_embryo), 75, axis=0), index=feats, columns=feats)
    medians = cells.groupby(["Group", "EmbryoID"])[feats].median()
    between = medians.corr(method="spearman")
    within.to_csv(out / "rho_within_embryo_median.csv")
    between.to_csv(out / "rho_between_embryo.csv")

    pairs = []
    for a, b in itertools.combinations(feats, 2):
        pairs.append({"FeatureA": a, "FeatureB": b, "RhoWithinMedian": within.loc[a, b],
                      "RhoWithinIQR": f"{q1.loc[a, b]:.2f} to {q3.loc[a, b]:.2f}",
                      "RhoBetween": between.loc[a, b],
                      "Redundant": abs(within.loc[a, b]) >= args.redundant
                      and abs(between.loc[a, b]) >= args.redundant})
    pairs = pd.DataFrame(pairs).sort_values("RhoWithinMedian", key=np.abs, ascending=False)
    pairs.to_csv(out / "redundant_pairs.csv", index=False)
    print(f"{len(per_embryo)} embryos with >= {MIN_CELLS} cells")
    print(pairs.round(2).to_string(index=False))

    mpl = common.plot_style()
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), constrained_layout=True)
    short = [f.replace("Spectral", "").replace("Fraction", "Frac") for f in feats]
    for ax, mat, title in ((axes[0], within, "Within embryo (median over embryos)"),
                           (axes[1], between, "Between embryos (embryo medians)")):
        im = ax.imshow(mat.to_numpy(), cmap="RdBu_r", vmin=-1, vmax=1)
        ax.set_xticks(range(len(feats)), short, rotation=45, ha="right")
        ax.set_yticks(range(len(feats)), short)
        ax.grid(False)
        for i, j in itertools.product(range(len(feats)), repeat=2):
            v = mat.iloc[i, j]
            ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=9,
                    color="white" if abs(v) > 0.6 else common.TEXT)
        ax.set_title(title)
    fig.colorbar(im, ax=axes, label="Spearman rho", shrink=0.8)
    fig.suptitle(f"Spectral feature redundancy ({args.cells} mapped cells)")
    for ext in ("png", "pdf"):
        fig.savefig(out / f"feature_correlations.{ext}")
    plt.close(fig)


if __name__ == "__main__":
    main()
