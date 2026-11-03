#!/usr/bin/env python3
"""Cross-group spatial outputs on shared scales, plus neural vs non-neural summaries.

Run after the four per-group runs of spatial_spectral_mapping.py. Reads
<mapping-root>/<group>/All_MappedCells_with_CNPU.csv and
Group_Cnpu_Spectral_Summary.csv and writes to --output-dir:

  Registration_QC.csv
      Per embryo: R^2 of an affine fit from tracked image position
      (MeanX/MeanY) to template centroid (Centroid_1/2). The fit should be
      near-perfect for a correct CellID <-> Calcium_ID join; low values flag
      boundary files that do not match the recording.
  SharedScale_<feature>_<metric>.png
      2x2 group maps on one colour scale, 9-px CNPUs, Sox2 boundary overlay.
      z-score maps use the paper's fixed -1..1 range (saturating beyond).
  Cell_Region_Assignments.csv
      Every mapped cell with Region = neural | non-neural, from the
      stage-matched Sox2 template polygon (Stage 13 groups -> Stage 14
      template; Stage 14 groups -> Stage 15 template, as in the mapping).
  Embryo_Region_Summary.csv
      Per embryo x region x feature: n cells, median, mean.
  Region_Comparison_CellLevel.csv
      Paper-style (Methods 2.3/2.5): per group x feature, two-sided
      two-sample K-S test and Cohen's d, neural vs non-neural cells.
      Bonferroni across all group x feature tests. Significant if
      p_adj < 0.05 and |d| > 0.2.
  Region_Comparison_EmbryoLevel.csv
      Robustness check that treats the embryo as the unit: Wilcoxon
      signed-rank test on paired embryo medians (embryos with >= --min-cells
      in both regions), Bonferroni-adjusted.
  Region_EmbryoMedians_<feature>.png
      Paired embryo medians, neural vs non-neural, per group.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

FEATURES = ("DominantFreqHz", "SpectralCentroid", "SpectralSpread",
            "SpectralEntropy", "TotalBandPower", "PeakPowerFraction")
GROUPS = ("S13", "S13_devitellinated", "S14", "S14_devitellinated")
TEMPLATE_FOR_GROUP = {"S13": "stage14", "S13_devitellinated": "stage14",
                      "S14": "stage15", "S14_devitellinated": "stage15"}
X_RANGE, Y_RANGE = (120, 1122), (1, 1190)
METRICS = {  # metric column -> (label, fixed limits or None)
    "GroupMedian": ("Embryo-balanced group median", None),
    "PooledMean": ("Pooled cell mean (paper CNP method)", None),
    "SmoothedZ": ("Smoothed group-median z-score", (-1.0, 1.0)),
    "PooledSmoothedZ": ("Smoothed pooled-mean z-score (paper CNP method)", (-1.0, 1.0)),
    "CoveragePct": ("Embryo coverage (%)", (0.0, 100.0)),
}


def load_template(path):
    table = pd.read_csv(path)
    lower = {c.lower(): c for c in table.columns}
    xcol = next((lower[c] for c in ("x", "v1", "col", "boundary_1", "var1") if c in lower), None)
    ycol = next((lower[c] for c in ("y", "v2", "row", "boundary_2", "var2") if c in lower), None)
    if xcol is None or ycol is None:
        if table.shape[1] != 2:
            raise ValueError(f"Cannot identify X/Y columns in {path}")
        xcol, ycol = table.columns
    pts = table[[xcol, ycol]].apply(pd.to_numeric, errors="coerce").dropna().to_numpy()
    if len(pts) < 3:
        raise ValueError(f"Fewer than three boundary points in {path}")
    if not np.allclose(pts[0], pts[-1]):
        pts = np.vstack([pts, pts[0]])
    return pts


def affine_r2(df):
    if len(df) < 4:
        return np.nan
    a = np.c_[df["MeanX"], df["MeanY"], np.ones(len(df))]
    out = []
    for target in ("Centroid_1", "Centroid_2"):
        y = df[target].to_numpy(float)
        coef = np.linalg.lstsq(a, y, rcond=None)[0]
        ss = ((y - y.mean()) ** 2).sum()
        out.append(1 - ((y - a @ coef) ** 2).sum() / ss if ss > 0 else np.nan)
    return float(np.nanmin(out))


def cohens_d(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 2 or len(b) < 2:
        return np.nan
    pooled = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1))
                     / (len(a) + len(b) - 2))
    return (a.mean() - b.mean()) / pooled if pooled > 0 else np.nan


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--mapping-root", type=Path, required=True)
    parser.add_argument("--stage14-sox2", type=Path, required=True)
    parser.add_argument("--stage15-sox2", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--exclude-embryos", default="",
                        help="Comma-separated EmbryoIDs (or a file, one per line) to drop.")
    parser.add_argument("--min-embryos", type=int, default=1,
                        help="Hide CNPUs supported by fewer embryos in shared-scale maps.")
    parser.add_argument("--min-cells", type=int, default=10,
                        help="Min cells per region for an embryo in the embryo-level test.")
    args = parser.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.path import Path as MplPath
    from scipy import stats

    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    exclude = set()
    if args.exclude_embryos:
        source = Path(args.exclude_embryos)
        text = source.read_text() if source.is_file() else args.exclude_embryos.replace(",", "\n")
        exclude = {line.strip() for line in text.splitlines() if line.strip()}

    templates = {"stage14": load_template(args.stage14_sox2),
                 "stage15": load_template(args.stage15_sox2)}
    cells, grids = [], []
    for group in GROUPS:
        gdir = args.mapping_root / group
        cells.append(pd.read_csv(gdir / "All_MappedCells_with_CNPU.csv", low_memory=False))
        grids.append(pd.read_csv(gdir / "Group_Cnpu_Spectral_Summary.csv"))
    cells, grids = pd.concat(cells, ignore_index=True), pd.concat(grids, ignore_index=True)

    # 1. Registration QC (all embryos, before exclusion).
    qc = (cells.groupby(["Group", "EmbryoID"])
          .apply(lambda d: pd.Series({"NMappedCells": len(d), "AffineR2": affine_r2(d)}),
                 include_groups=False).reset_index())
    qc["Flag_R2_below_0.9"] = qc["AffineR2"] < 0.9
    qc["Excluded"] = qc["EmbryoID"].isin(exclude)
    qc.sort_values("AffineR2").to_csv(out / "Registration_QC.csv", index=False)
    if exclude:
        unknown = exclude - set(cells["EmbryoID"])
        if unknown:
            print(f"Warning: exclusion IDs not found in mapped cells: {sorted(unknown)}")
        cells = cells[~cells["EmbryoID"].isin(exclude)]
        print(f"Excluded {len(exclude)} embryo(s) from region summaries. CNPU maps use the "
              "per-group summaries as computed; rerun per-group mapping without these "
              "embryos to drop them from the maps too.")

    # 2. Shared-scale maps.
    for feature in FEATURES:
        sub = grids[grids["Feature"] == feature]
        for metric, (label, fixed) in METRICS.items():
            if metric not in sub.columns:
                print(f"Skipping {metric}: not in group summaries (rerun v2 mapping).")
                continue
            shown = sub[(sub["NEmbryos"] >= args.min_embryos) & sub[metric].notna()]
            if shown.empty:
                continue
            limits = fixed or tuple(np.nanpercentile(shown[metric], [1, 99]))
            fig, axes = plt.subplots(2, 2, figsize=(12, 10), constrained_layout=True)
            last = None
            for ax, group in zip(axes.flat, GROUPS):
                g = shown[shown["Group"] == group]
                last = ax.scatter(g["CNPU_CenterX"], g["CNPU_CenterY"], c=g[metric], s=6,
                                  marker="s", cmap="jet" if fixed and fixed[0] < 0 else "viridis",
                                  vmin=limits[0], vmax=limits[1])
                poly = templates[TEMPLATE_FOR_GROUP[group]]
                ax.plot(poly[:, 0], poly[:, 1], "k:", linewidth=1.5)
                ax.set(xlim=X_RANGE, ylim=Y_RANGE, title=group, xlabel="SuperEmbryo X",
                       ylabel="SuperEmbryo Y")
                ax.invert_yaxis()
                ax.set_aspect("equal")
            fig.colorbar(last, ax=axes, label=label, shrink=0.8,
                         extend="both" if fixed and fixed[0] < 0 else "neither")
            fig.suptitle(f"{feature} — {label} (shared scale, CNPUs with ≥{args.min_embryos} embryo(s))")
            fig.savefig(out / f"SharedScale_{feature}_{metric}.png", dpi=200)
            plt.close(fig)

    # 3. Neural vs non-neural assignment.
    cells = cells.copy()
    cells["Region"] = "non-neural"
    for group in GROUPS:
        mask = cells["Group"] == group
        path = MplPath(templates[TEMPLATE_FOR_GROUP[group]])
        inside = path.contains_points(cells.loc[mask, ["Centroid_1", "Centroid_2"]].to_numpy(float))
        cells.loc[mask, "Region"] = np.where(inside, "neural", "non-neural")
    keep_cols = ["Group", "EmbryoID", "CellID", "Centroid_1", "Centroid_2", "CNPU_ID",
                 "Region", *FEATURES]
    cells[keep_cols].to_csv(out / "Cell_Region_Assignments.csv", index=False)

    long = cells.melt(id_vars=["Group", "EmbryoID", "Region"], value_vars=list(FEATURES),
                      var_name="Feature", value_name="Value").dropna(subset=["Value"])
    embryo = (long.groupby(["Group", "EmbryoID", "Region", "Feature"])["Value"]
              .agg(NCells="size", Median="median", Mean="mean").reset_index())
    embryo.to_csv(out / "Embryo_Region_Summary.csv", index=False)

    cell_rows, embryo_rows = [], []
    for group in GROUPS:
        for feature in FEATURES:
            sel = long[(long["Group"] == group) & (long["Feature"] == feature)]
            a = sel.loc[sel["Region"] == "neural", "Value"].to_numpy()
            b = sel.loc[sel["Region"] == "non-neural", "Value"].to_numpy()
            ks = stats.ks_2samp(a, b) if len(a) and len(b) else None
            cell_rows.append({
                "Group": group, "Feature": feature, "NNeuralCells": len(a), "NNonNeuralCells": len(b),
                "MeanNeural": np.mean(a) if len(a) else np.nan,
                "MeanNonNeural": np.mean(b) if len(b) else np.nan,
                "KS_D": ks.statistic if ks else np.nan, "KS_p": ks.pvalue if ks else np.nan,
                "CohensD": cohens_d(a, b)})
            e = embryo[(embryo["Group"] == group) & (embryo["Feature"] == feature)]
            wide = e.pivot(index="EmbryoID", columns="Region", values="Median")
            counts = e.pivot(index="EmbryoID", columns="Region", values="NCells")
            ok = counts.reindex(columns=["neural", "non-neural"]).fillna(0).ge(args.min_cells).all(axis=1)
            paired = wide.reindex(columns=["neural", "non-neural"])[ok].dropna()
            if len(paired) >= 6:
                w = stats.wilcoxon(paired["neural"], paired["non-neural"])
                stat, p = w.statistic, w.pvalue
            else:
                stat = p = np.nan
            embryo_rows.append({
                "Group": group, "Feature": feature, "NEmbryosPaired": len(paired),
                "MedianOfDifferences": float((paired["neural"] - paired["non-neural"]).median())
                if len(paired) else np.nan,
                "Wilcoxon_W": stat, "Wilcoxon_p": p})
    for rows, pcol in ((cell_rows, "KS_p"), (embryo_rows, "Wilcoxon_p")):
        n_tests = sum(np.isfinite(r[pcol]) for r in rows)
        for r in rows:
            r[pcol + "_bonferroni"] = min(1.0, r[pcol] * n_tests) if np.isfinite(r[pcol]) else np.nan
    cell_df = pd.DataFrame(cell_rows)
    cell_df["Significant_paper_rule"] = (cell_df["KS_p_bonferroni"] < 0.05) & (cell_df["CohensD"].abs() > 0.2)
    cell_df.to_csv(out / "Region_Comparison_CellLevel.csv", index=False)
    pd.DataFrame(embryo_rows).to_csv(out / "Region_Comparison_EmbryoLevel.csv", index=False)

    for feature in FEATURES:
        e = embryo[embryo["Feature"] == feature]
        fig, axes = plt.subplots(1, 4, figsize=(14, 4), sharey=True, constrained_layout=True)
        for ax, group in zip(axes, GROUPS):
            wide = e[e["Group"] == group].pivot(index="EmbryoID", columns="Region", values="Median")
            wide = wide.reindex(columns=["neural", "non-neural"]).dropna()
            for _, r in wide.iterrows():
                ax.plot([0, 1], [r["neural"], r["non-neural"]], color="#999999", linewidth=0.7)
            ax.boxplot([wide["neural"], wide["non-neural"]], positions=[0, 1], widths=0.4)
            ax.set_xticks([0, 1], ["neural", "non-neural"])
            ax.set_title(f"{group} (n={len(wide)})")
        axes[0].set_ylabel(f"Embryo median {feature}")
        fig.savefig(out / f"Region_EmbryoMedians_{feature}.png", dpi=200)
        plt.close(fig)

    print(cells.groupby(["Group", "Region"]).size().unstack().to_string())
    print(f"Outputs: {out}")


if __name__ == "__main__":
    main()
