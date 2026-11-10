#!/usr/bin/env python3
"""Step 6: poster / manuscript figure panels from the outputs of steps 1-5.

Each panel is written as PDF (vector, for the poster layout) and PNG to
ANALYSIS_ROOT/step6_figures/. A panel whose inputs are missing is skipped
with a message, so this can be rerun as steps finish.

Panels (--panels, default all):
  flow        cell flow per group: lost at spectral QC / not mapped /
              outside grid / mapped (non-rhythmic) / rhythmic
  spectra     example cells: observed spectrum vs AR(1) mean and 95% band
              (one rhythmic, one not; --example-embryo)
  maps        2x2 shared-scale smoothed z-maps for --features (paper -1..1 jet
              scale, CNPUs with >= --min-embryos embryos, Sox2 outline)
  rhythm_map  % rhythmic cells per CNPU (smoothed), sequential blue scale
  gradient    per-embryo anterior-posterior rho for --features (step 4a)
  effects     per-embryo medians for --features + PctRhythmic (step 4b)
  events      per-embryo event-spectral rho for --event-column (step 5b)

Encoding: x position = stage, colour = membrane (blue intact, orange
devitellinated), legend always present; one dot per embryo, black bar =
group median. Z-maps keep the paper's jet scale for consistency with the
PDT figures (project rule).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "analysis_common"))
import common  # noqa: E402

STEP = common.ANALYSIS_ROOT
FLOW_ORDER = ["Lost at spectral QC", "Not in boundary file", "Outside CNPU grid",
              "Mapped, not rhythmic", "Mapped, rhythmic"]
FLOW_COLORS = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#0d366b"]  # ordinal blue ramp
XPOS = {("intact", 13): 0.0, ("devitellinated", 13): 0.35,
        ("intact", 14): 1.2, ("devitellinated", 14): 1.55}


def save(fig, out: Path, name: str):
    for ext in ("pdf", "png"):
        fig.savefig(out / f"{name}.{ext}")
    import matplotlib.pyplot as plt
    plt.close(fig)
    print(f"wrote {name}")


def need(path: Path) -> bool:
    if not path.exists():
        print(f"  skipped: {path} not found")
        return False
    return True


def dot_axes(ax, df, value, ylabel, zero_line=False):
    """One dot per embryo at (stage, membrane) positions, median bar, legend."""
    rng = np.random.default_rng(0)
    for membrane in ("intact", "devitellinated"):
        for stage in (13, 14):
            v = df[(df.Membrane == membrane) & (df.Stage == stage)][value].dropna()
            x0 = XPOS[(membrane, stage)]
            ax.scatter(x0 + rng.uniform(-0.09, 0.09, len(v)), v, s=22,
                       color=common.MEMBRANE_COLOR[membrane], edgecolor="white", linewidth=0.6,
                       label=membrane if stage == 13 else None, zorder=3)
            if len(v):
                ax.plot([x0 - 0.13, x0 + 0.13], [v.median()] * 2, color=common.TEXT, lw=2, zorder=4)
    if zero_line:
        ax.axhline(0, color=common.TEXT_MUTED, lw=0.8, ls="--", zorder=1)
    ax.set_xticks([0.175, 1.375], ["Stage 13", "Stage 14"])
    ax.set_xlim(-0.35, 1.9)
    ax.set_ylabel(ylabel)
    ax.grid(axis="x", visible=False)


def panel_flow(out, args):
    path = STEP / "step2_cell_flow/cell_flow_by_group.csv"
    if not need(path):
        return
    import matplotlib.pyplot as plt
    g = pd.read_csv(path, index_col=0).drop(index="All", errors="ignore").reindex(common.GROUPS)
    rhythmic = g["Rhythmic"] if "Rhythmic" in g else pd.Series(0, index=g.index)
    parts = pd.DataFrame({
        FLOW_ORDER[0]: g.Tracked - g.PassedSpectralQC,
        FLOW_ORDER[1]: g.PassedSpectralQC - g.Mapped,
        FLOW_ORDER[2]: g.Mapped - g.InGrid,
        FLOW_ORDER[3]: g.InGrid - rhythmic,
        FLOW_ORDER[4]: rhythmic})
    if "Rhythmic" not in g:
        parts = parts.drop(columns=FLOW_ORDER[4]).rename(columns={FLOW_ORDER[3]: "Mapped (in grid)"})
    pct = 100 * parts.div(g.Tracked, axis=0)
    fig, ax = plt.subplots(figsize=(9, 3.6))
    left = np.zeros(len(pct))
    for col, color in zip(pct.columns, FLOW_COLORS[:len(pct.columns)]):
        ax.barh(range(len(pct)), pct[col], left=left, color=color, edgecolor="white",
                linewidth=2, label=col, height=0.6)
        left += pct[col].to_numpy()
    ax.set_yticks(range(len(pct)), [f"{i}\n(n = {int(n):,} tracked)" for i, n in zip(pct.index, g.Tracked)])
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("% of tracked cells")
    ax.grid(axis="y", visible=False)
    ax.legend(ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.22))
    ax.set_title("Cells retained at each analysis step")
    save(fig, out, "fig_cell_flow")


def panel_spectra(out, args):
    root = STEP / "step1_rhythmicity"
    table = root / common.RHYTHM_FILE
    if not need(table):
        return
    import matplotlib.pyplot as plt
    r = pd.read_csv(table, low_memory=False)
    r["CellID"] = r["CellID"].astype(str)
    emb = args.example_embryo or r.groupby("EmbryoID")["Rhythmic"].mean().idxmax()
    d = r[r.EmbryoID == emb]
    group = d.Group.iloc[0]
    lib = np.load(common.SPECTRAL_ROOT / common.GROUP_FOLDER[group] / emb / "SpectralLibrary.npz")
    red = np.load(root / group / f"{emb}_rednoise.npz")
    lib_idx = {c: i for i, c in enumerate(lib["CellID"].astype(str))}
    red_idx = {c: i for i, c in enumerate(red["CellID"].astype(str))}
    picks = [("Rhythmic cell", d[d.Rhythmic].sort_values("MaxRatio").iloc[-1].CellID),
             ("Non-rhythmic cell", d[~d.Rhythmic].sort_values("Q_global").iloc[-1].CellID)]
    f = lib["Frequency"] * 1000
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.6), sharey=False, constrained_layout=True)
    for ax, (title, cid) in zip(axes, picks):
        ax.fill_between(f, 0, red["Surrogate95"][red_idx[cid]], color=common.GRID,
                        label="AR(1) 95th percentile")
        ax.plot(f, red["SurrogateMean"][red_idx[cid]], color=common.TEXT_MUTED, lw=1.5,
                ls="--", label="AR(1) mean")
        ax.plot(f, lib["Power"][lib_idx[cid]], color=common.MEMBRANE_COLOR["intact"], lw=2,
                label="Observed")
        ax.set(xlabel="Frequency (mHz)", ylabel="Mean wavelet power", title=f"{title} ({cid})")
    axes[0].legend()
    fig.suptitle(f"Red-noise test, {emb}")
    save(fig, out, "fig_example_spectra")


def map_axes(fig, axes, grid, value, cmap, vmin, vmax, label, extend):
    last = None
    for ax, g in zip(axes.flat, common.GROUPS):
        d = grid[grid.Group == g]
        last = ax.scatter(d.CNPU_CenterX, d.CNPU_CenterY, c=d[value], s=5, marker="s",
                          cmap=cmap, vmin=vmin, vmax=vmax)
        poly = common.load_template(common.SOX2[common.TEMPLATE_FOR_GROUP[g]])
        ax.plot(poly[:, 0], poly[:, 1], ":", color=common.TEXT, lw=1.5)
        ax.set(xlim=common.X_RANGE, ylim=common.Y_RANGE, title=g.replace("_", " "))
        ax.invert_yaxis()
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.grid(False)
        for s in ax.spines.values():
            s.set_visible(False)
    fig.colorbar(last, ax=axes, label=label, shrink=0.75, extend=extend)


def panel_maps(out, args):
    import matplotlib.pyplot as plt
    files = [args.mapping_root / g / "Group_Cnpu_Spectral_Summary.csv" for g in common.GROUPS]
    if not all(need(p) for p in files):
        return
    grid = pd.concat([pd.read_csv(p) for p in files])
    for feat in args.features:
        d = grid[(grid.Feature == feat) & (grid.NEmbryos >= args.min_embryos)]
        fig, axes = plt.subplots(2, 2, figsize=(8, 8), constrained_layout=True)
        map_axes(fig, axes, d, "PooledSmoothedZ", "jet", -1, 1, f"{feat}, smoothed z-score", "both")
        fig.suptitle(f"{feat} (CNPUs with ≥{args.min_embryos} embryos)")
        save(fig, out, f"fig_map_{feat}_{args.mapping_root.name}")


def panel_rhythm_map(out, args):
    if not need(STEP / "step1_rhythmicity" / common.RHYTHM_FILE):
        return
    import matplotlib.pyplot as plt
    cells = common.load_mapped_cells()
    r = common.load_rhythmicity()[["Group", "EmbryoID", "CellID", "Rhythmic"]]
    cells = cells.merge(r, on=["Group", "EmbryoID", "CellID"], how="inner")
    cells["PctRhythmic"] = 100 * cells["Rhythmic"].astype(float)
    grid = common.cnpu_pooled_map(cells, "PctRhythmic")
    grid = grid[grid.NEmbryos >= args.min_embryos]
    vmax = float(np.nanpercentile(grid.PooledSmoothed, 99))
    fig, axes = plt.subplots(2, 2, figsize=(8, 8), constrained_layout=True)
    map_axes(fig, axes, grid, "PooledSmoothed", common.sequential_cmap(), 0, vmax,
             "% rhythmic cells (49-CNPU smoothed)", "max")
    fig.suptitle(f"Rhythmic cells (red-noise test; CNPUs with ≥{args.min_embryos} embryos)")
    save(fig, out, "fig_map_pct_rhythmic")


def panel_gradient(out, args):
    path = STEP / f"step4_embryo_stats/{args.cells}/gradient_ap_by_embryo.csv"
    if not need(path):
        return
    import matplotlib.pyplot as plt
    e = pd.read_csv(path)
    fig, axes = plt.subplots(1, len(args.features), figsize=(4.2 * len(args.features), 3.8),
                             constrained_layout=True, squeeze=False)
    for ax, feat in zip(axes[0], args.features):
        dot_axes(ax, e[e.Feature == feat], "Rho", f"{feat}\nrho with anterior position", zero_line=True)
        ax.set_title(feat)
    axes[0][0].legend(loc="lower left")
    fig.suptitle("Anterior-posterior gradient within the neural plate (one dot per embryo)")
    save(fig, out, f"fig_ap_gradient_{args.cells}")


def panel_effects(out, args):
    path = STEP / f"step4_embryo_stats/{args.cells}/embryo_medians.csv"
    if not need(path):
        return
    import matplotlib.pyplot as plt
    m = pd.read_csv(path)
    m = m[m.Scope == "all"]
    outcomes = [o for o in [*args.features, "PctRhythmic"] if o in set(m.Outcome)]
    fig, axes = plt.subplots(1, len(outcomes), figsize=(4.2 * len(outcomes), 3.8),
                             constrained_layout=True, squeeze=False)
    for ax, o in zip(axes[0], outcomes):
        d = m[m.Outcome == o]
        dot_axes(ax, d, "Value", "% rhythmic cells" if o == "PctRhythmic" else f"Embryo median {o}")
        if o in common.LOG_FEATURES:
            ax.set_yscale("log")
        ax.set_title(o)
    axes[0][0].legend(loc="best")
    fig.suptitle("Stage and devitellination effects (one dot per embryo)")
    save(fig, out, f"fig_group_effects_{args.cells}")


def panel_events(out, args):
    path = STEP / "step5_event_integration/event_spectral_rho_by_embryo.csv"
    if not need(path):
        return
    import matplotlib.pyplot as plt
    rho = common.add_design(pd.read_csv(path))
    rho = rho[rho.EventColumn == args.event_column]
    feats = [f for f in common.FEATURES if f in set(rho.Feature)]
    fig, axes = plt.subplots(1, len(feats), figsize=(3.4 * len(feats), 3.6),
                             constrained_layout=True, squeeze=False, sharey=True)
    for ax, f in zip(axes[0], feats):
        dot_axes(ax, rho[rho.Feature == f], "Rho", "Spearman rho (cells within embryo)", zero_line=True)
        ax.set_title(f.replace("Spectral", ""))
    axes[0][0].legend(loc="lower left")
    fig.suptitle(f"PDT {args.event_column} vs spectral features")
    save(fig, out, f"fig_events_{args.event_column}")


PANELS = {"flow": panel_flow, "spectra": panel_spectra, "maps": panel_maps,
          "rhythm_map": panel_rhythm_map, "gradient": panel_gradient,
          "effects": panel_effects, "events": panel_events}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panels", nargs="+", default=list(PANELS), choices=list(PANELS))
    ap.add_argument("--features", nargs="+", default=["DominantFreqHz", "TotalBandPower"],
                    help="Choose from step 3 (two non-redundant features).")
    ap.add_argument("--cells", choices=("all", "rhythmic"), default="all")
    ap.add_argument("--mapping-root", type=Path, default=common.MAPPING_ROOT,
                    help="mapping_v3 or mapping_v3_rhythmic")
    ap.add_argument("--min-embryos", type=int, default=2)
    ap.add_argument("--example-embryo", default="")
    ap.add_argument("--event-column", default="Th_SD_1")
    args = ap.parse_args()
    common.plot_style()
    out = common.out_dir("step6_figures")
    for name in args.panels:
        print(f"[{name}]")
        PANELS[name](out, args)


if __name__ == "__main__":
    main()
