#!/usr/bin/env python3
"""Step 5c: PDT event metrics on the same CNPU grid as the spectral maps.

1. Maps: for each --event-columns entry, the paper CNP method (pooled mean per
   9-px CNPU, 49-CNPU smoothing over occupied CNPUs, z-score) per group, drawn
   as 2x2 shared-scale maps on the -1..1 jet scale with the Sox2 outline and
   CNPUs supported by >= --min-embryos embryos (same rules as the spectral
   maps in mapping_v3/combined, so the two can be shown side by side).
2. Spatial concordance, embryo = unit: inside each embryo, pooled means per
   CNPU (CNPUs with >= 2 cells) for the event metric and each spectral
   feature, Spearman rho across CNPUs (>= --min-cnpus). Per group: median rho
   and Wilcoxon vs 0. (Correlating the two group maps directly is reported
   too, labelled descriptive: neighbouring CNPUs are not independent.)

Needs events_spectral_cells.csv from event_spectral_join.py.
Outputs in ANALYSIS_ROOT/step5_event_integration/maps/:
  <event>_cnpu.csv, SharedScale_<event>_PooledSmoothedZ.png/.pdf,
  spatial_concordance_by_embryo.csv, spatial_concordance_tests.csv,
  group_map_correlation_descriptive.csv
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


def draw_map(grid: pd.DataFrame, column: str, title: str, path: Path, min_embryos: int):
    common.plot_style()
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(11, 10), constrained_layout=True)
    last = None
    for ax, g in zip(axes.flat, common.GROUPS):
        d = grid[(grid.Group == g) & (grid.NEmbryos >= min_embryos)]
        last = ax.scatter(d.CNPU_CenterX, d.CNPU_CenterY, c=d[column], s=6, marker="s",
                          cmap="jet", vmin=-1, vmax=1)
        poly = common.load_template(common.SOX2[common.TEMPLATE_FOR_GROUP[g]])
        ax.plot(poly[:, 0], poly[:, 1], ":", color=common.TEXT, linewidth=1.5)
        ax.set(xlim=common.X_RANGE, ylim=common.Y_RANGE, title=g,
               xlabel="SuperEmbryo X", ylabel="SuperEmbryo Y")
        ax.invert_yaxis()
        ax.set_aspect("equal")
        ax.grid(False)
    fig.colorbar(last, ax=axes, label="Smoothed CNPU mean z-score", shrink=0.8, extend="both")
    fig.suptitle(title)
    for ext in ("png", "pdf"):
        fig.savefig(path.with_suffix(f".{ext}"))
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--event-columns", nargs="+", required=True,
                    help="e.g. Th_SD_1 MeanSpikeMag (confirm definitions first)")
    ap.add_argument("--min-embryos", type=int, default=2)
    ap.add_argument("--min-cnpus", type=int, default=20)
    args = ap.parse_args()
    base = common.ANALYSIS_ROOT / "step5_event_integration"
    out = common.out_dir("step5_event_integration", "maps")
    cells = pd.read_csv(base / "events_spectral_cells.csv", low_memory=False)
    cells = cells[cells["CNPU_X"].notna()].copy()
    cells["InGrid"] = cells["InGrid"].astype(bool)

    spectral_grid = pd.concat([pd.read_csv(common.MAPPING_ROOT / g / "Group_Cnpu_Spectral_Summary.csv")
                               for g in common.GROUPS])
    map_corr = []
    for ec in args.event_columns:
        grid = common.cnpu_pooled_map(cells, ec)
        grid.to_csv(out / f"{ec}_cnpu.csv", index=False)
        draw_map(grid, "PooledSmoothedZ",
                 f"{ec} — smoothed CNPU mean z-score (CNPUs with ≥{args.min_embryos} embryos; "
                 "embryos with event files only)", out / f"SharedScale_{ec}_PooledSmoothedZ",
                 args.min_embryos)
        for f in common.FEATURES:
            s = spectral_grid[spectral_grid.Feature == f][["Group", "CNPU_X", "CNPU_Y", "PooledSmoothedZ"]]
            m = grid.merge(s, on=["Group", "CNPU_X", "CNPU_Y"], suffixes=("_event", "_spectral"))
            m = m[m.NEmbryos >= args.min_embryos].dropna(subset=["PooledSmoothedZ_event",
                                                                  "PooledSmoothedZ_spectral"])
            for g, d in m.groupby("Group"):
                map_corr.append({"Group": g, "EventColumn": ec, "Feature": f, "NCNPUs": len(d),
                                 "SpearmanRho_descriptive": stats.spearmanr(
                                     d.PooledSmoothedZ_event, d.PooledSmoothedZ_spectral).statistic})
    pd.DataFrame(map_corr).to_csv(out / "group_map_correlation_descriptive.csv", index=False)

    rows = []
    in_grid = cells[cells.InGrid]
    for (g, e), d in in_grid.groupby(["Group", "EmbryoID"]):
        per = d.groupby(["CNPU_X", "CNPU_Y"])
        n = per.size()
        means = per[[*args.event_columns, *common.FEATURES]].mean()[n >= 2]
        if len(means) < args.min_cnpus:
            continue
        for ec in args.event_columns:
            for f in common.FEATURES:
                rows.append({"Group": g, "EmbryoID": e, "EventColumn": ec, "Feature": f,
                             "NCNPUs": len(means),
                             "Rho": stats.spearmanr(means[ec], means[f], nan_policy="omit").statistic})
    emb = pd.DataFrame(rows)
    emb.to_csv(out / "spatial_concordance_by_embryo.csv", index=False)
    tests = []
    for (g, ec, f), d in pd.concat([emb, emb.assign(Group="All")]).groupby(["Group", "EventColumn", "Feature"]):
        r = d["Rho"].dropna()
        w = stats.wilcoxon(r) if len(r) >= 6 else None
        tests.append({"Group": g, "EventColumn": ec, "Feature": f, "NEmbryos": len(r),
                      "MedianRho": r.median(), "Wilcoxon_p": w.pvalue if w else np.nan})
    t = pd.DataFrame(tests)
    t["q_bh"] = common.bh(t["Wilcoxon_p"])
    t.to_csv(out / "spatial_concordance_tests.csv", index=False)
    print(t[t.Group == "All"].round(4).to_string(index=False))


if __name__ == "__main__":
    main()
