#!/usr/bin/env python3
"""Step 4b: stage x membrane effects on spectral features, embryo as the unit.

Primary analysis (inference):
  Embryo medians per feature in three scopes (all mapped cells, neural,
  non-neural; >= --min-cells per embryo and scope), then the 2x2 model
  value ~ stage * membrane with Freedman-Lane permutation p-values
  (common.two_by_two). Holm across features within each scope x term.
  If step 1 exists, PctRhythmic (per embryo) is added as a seventh outcome.

Planned pairwise contrasts on embryo medians (Mann-Whitney U + Hedges g, Holm):
  S13 vs S14 (stage, intact)           S13_dev vs S14_dev (stage, devitellinated)
  S13 vs S13_dev (membrane, stage 13)  S14 vs S14_dev (membrane, stage 14)

Paper-parity table (descriptive, cells treated as independent as in the PDT
paper Methods 2.5): the same four contrasts on cells with two-sample K-S,
Cohen's d, Bonferroni, paper rule p_adj < 0.05 and |d| > 0.2.

Optional --lmm: cell-level linear mixed model
  value ~ stage * membrane, random intercept per embryo (statsmodels;
  log10 for TotalBandPower). Requires `conda install -n wavelet-analysis statsmodels`.

Outputs in ANALYSIS_ROOT/step4_embryo_stats/<cells>/:
  embryo_medians.csv, stage_membrane_effects.csv, pairwise_embryo_level.csv,
  pairwise_paper_style_cells.csv, [lmm_effects.csv]
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

CONTRASTS = (("S13", "S14", "Stage (intact)"),
             ("S13_devitellinated", "S14_devitellinated", "Stage (devitellinated)"),
             ("S13", "S13_devitellinated", "Membrane (stage 13)"),
             ("S14", "S14_devitellinated", "Membrane (stage 14)"))


def embryo_medians(cells: pd.DataFrame, min_cells: int) -> pd.DataFrame:
    scopes = {"all": cells, "neural": cells[cells.Region == "neural"],
              "non-neural": cells[cells.Region == "non-neural"]}
    rows = []
    for scope, d in scopes.items():
        grouped = d.groupby(["Group", "EmbryoID"])
        med = grouped[list(common.FEATURES)].median()
        n = grouped.size().rename("NCells")
        med = med[n >= min_cells]
        long = med.reset_index().melt(id_vars=["Group", "EmbryoID"], var_name="Outcome",
                                      value_name="Value")
        long = long.merge(n.reset_index(), on=["Group", "EmbryoID"], how="left")
        long["Scope"] = scope
        rows.append(long)
    return common.add_design(pd.concat(rows, ignore_index=True))


def rhythm_outcome(cells: pd.DataFrame, min_cells: int) -> pd.DataFrame | None:
    try:
        r = common.load_rhythmicity()[["Group", "EmbryoID", "CellID", "Rhythmic"]]
    except FileNotFoundError:
        return None
    d = cells.merge(r, on=["Group", "EmbryoID", "CellID"], how="inner")
    rows = []
    for scope, s in (("all", d), ("neural", d[d.Region == "neural"]),
                     ("non-neural", d[d.Region == "non-neural"])):
        g = s.groupby(["Group", "EmbryoID"])["Rhythmic"].agg(["mean", "size"])
        g = g[g["size"] >= min_cells].reset_index()
        rows.append(pd.DataFrame({"Group": g.Group, "EmbryoID": g.EmbryoID, "Outcome": "PctRhythmic",
                                  "Value": 100 * g["mean"], "Scope": scope, "NCells": g["size"]}))
    return common.add_design(pd.concat(rows, ignore_index=True))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cells", choices=("all", "rhythmic"), default="all")
    ap.add_argument("--min-cells", type=int, default=20)
    ap.add_argument("--n-perm", type=int, default=10000)
    ap.add_argument("--lmm", action="store_true")
    args = ap.parse_args()
    out = common.out_dir("step4_embryo_stats", args.cells)

    cells = common.filter_cells(common.load_mapped_cells(), args.cells)
    cells["Region"] = common.assign_region(cells)
    med = embryo_medians(cells, args.min_cells)
    if args.cells == "all":
        extra = rhythm_outcome(cells, args.min_cells)
        if extra is not None:
            med = pd.concat([med, extra], ignore_index=True)
    med.to_csv(out / "embryo_medians.csv", index=False)

    effects = []
    for (scope, outcome), d in med.groupby(["Scope", "Outcome"]):
        for r in common.two_by_two(d["Value"], d["Stage"], d["Membrane"], n_perm=args.n_perm):
            cell_means = d.groupby("Group")["Value"].median()
            effects.append({"Scope": scope, "Outcome": outcome, **r,
                            **{f"Median_{g}": cell_means.get(g, np.nan) for g in common.GROUPS}})
    eff = pd.DataFrame(effects)
    eff["p_holm"] = eff.groupby(["Scope", "Term"])["p_perm"].transform(lambda p: common.holm(p))
    eff.to_csv(out / "stage_membrane_effects.csv", index=False)
    print(eff[eff.Scope == "all"].round(4).to_string(index=False))

    pair = []
    for (scope, outcome), d in med.groupby(["Scope", "Outcome"]):
        for a, b, label in CONTRASTS:
            x, y = d.loc[d.Group == a, "Value"], d.loc[d.Group == b, "Value"]
            u = stats.mannwhitneyu(x, y) if len(x) >= 3 and len(y) >= 3 else None
            pair.append({"Scope": scope, "Outcome": outcome, "Contrast": label, "A": a, "B": b,
                         "nA": len(x), "nB": len(y), "MedianA": x.median(), "MedianB": y.median(),
                         "HedgesG_BminusA": common.hedges_g(y, x), "MWU_p": u.pvalue if u else np.nan})
    pw = pd.DataFrame(pair)
    pw["p_holm"] = pw.groupby(["Scope", "Contrast"])["MWU_p"].transform(lambda p: common.holm(p))
    pw.to_csv(out / "pairwise_embryo_level.csv", index=False)

    paper = []
    for f in common.FEATURES:
        for a, b, label in CONTRASTS:
            x = cells.loc[cells.Group == a, f].dropna()
            y = cells.loc[cells.Group == b, f].dropna()
            ks = stats.ks_2samp(x, y)
            paper.append({"Feature": f, "Contrast": label, "nCellsA": len(x), "nCellsB": len(y),
                          "KS_D": ks.statistic, "KS_p": ks.pvalue,
                          "CohensD_BminusA": common.cohens_d(y, x)})
    pp = pd.DataFrame(paper)
    pp["KS_p_bonferroni"] = np.minimum(1, pp["KS_p"] * len(pp))
    pp["Significant_paper_rule"] = (pp.KS_p_bonferroni < 0.05) & (pp.CohensD_BminusA.abs() > 0.2)
    pp.to_csv(out / "pairwise_paper_style_cells.csv", index=False)

    if args.lmm:
        try:
            import statsmodels.formula.api as smf
        except ImportError:
            sys.exit("statsmodels not installed: conda install -n wavelet-analysis statsmodels")
        rows = []
        for f in common.FEATURES:
            d = cells[["EmbryoID", "Stage", "Membrane", f]].dropna().copy()
            d["y"] = np.log10(d[f].clip(lower=1e-12)) if f in common.LOG_FEATURES else d[f]
            d["s"] = np.where(d.Stage == 14, 0.5, -0.5)
            d["m"] = np.where(d.Membrane == "devitellinated", 0.5, -0.5)
            fit = smf.mixedlm("y ~ s * m", d, groups=d["EmbryoID"]).fit(reml=True)
            for term, name in (("s", "Stage"), ("m", "Membrane"), ("s:m", "Interaction")):
                rows.append({"Feature": f, "Scale": "log10" if f in common.LOG_FEATURES else "raw",
                             "Term": name, "Estimate": fit.params[term], "SE": fit.bse[term],
                             "p_wald": fit.pvalues[term], "NCells": len(d),
                             "NEmbryos": d.EmbryoID.nunique(), "Converged": fit.converged})
        lm = pd.DataFrame(rows)
        lm["p_holm"] = lm.groupby("Term")["p_wald"].transform(lambda p: common.holm(p))
        lm.to_csv(out / "lmm_effects.csv", index=False)
        print(lm.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
