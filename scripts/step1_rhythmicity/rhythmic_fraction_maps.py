#!/usr/bin/env python3
"""Step 1d: where are rhythmic cells concentrated? Rhythmic-fraction CNPU maps.

For every CNPU the map value is the fraction of ALL mapped cells that are
rhythmic (Rhythmic = Q_global < 0.05 from build_rhythmic_tables.py). So every
CNPU carries the full all-cell embryo support of the v3 maps, unlike the
rhythmic-only feature maps (mapping_v3_rhythmic), where the median CNPU holds
1 cell from 1 embryo.

Method = the v3 / PDT-paper CNP method, applied to the 0/1 Rhythmic indicator:
  * cells and CNPU assignments: the audited all-cell v3 mapping
    (mapping_v3/<group>/All_MappedCells_with_CNPU.csv; mapping exclusions
    already applied there), repeated boundary IDs deduplicated, cells
    outside the CNPU grid dropped (as in the v3 maps);
  * PooledFraction = mean of the indicator over all cells of all embryos in
    the CNPU (the paper's pooled CNPU mean);
  * 7x7 (49-CNPU) mean over occupied CNPUs only, z-scored within group over
    occupied CNPUs (sample SD) -> PooledSmoothedZ, shown on the fixed -1..1
    jet scale; combined maps show CNPUs supported by >= --min-embryos.
  A self-check recomputes the production TotalBandPower PooledSmoothedZ with
  the same code and stops unless it matches mapping_v3.

Additions (not in the paper; secondary):
  * PooledSmoothedPct on one shared sequential scale for all groups: z-maps
    are relative within a group, so only this map compares groups.
  * EmbryoBalancedFraction = mean over embryos of each embryo's fraction in
    the CNPU (the v3 "embryo-balanced median" is degenerate for a 0/1
    indicator, so the mean is used), smoothed and z-scored the same way.
  * Support QC per CNPU over its 7x7 window: rhythmic cells, embryos that
    contribute rhythmic cells, and the share of those cells from the single
    top embryo (high share = a hot spot driven by one embryo).

Writes to --output-dir (must not exist yet, unless --overwrite):
  Group_Cnpu_RhythmicFraction.csv, Embryo_Cnpu_RhythmicFraction.csv,
  RhythmicFraction_*.png, run_settings.json
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "analysis_common"))
import common  # noqa: E402

NX, NY = common.MAX_CNPU_X + 1, common.MAX_CNPU_Y + 1
HALF = 3  # 7x7 window


def window_sum(a: np.ndarray) -> np.ndarray:
    """Sum over the 7x7 window (zero outside the grid), same shape as a."""
    from scipy.signal import convolve2d
    return convolve2d(a, np.ones((2 * HALF + 1, 2 * HALF + 1)), mode="same", boundary="fill")


def smooth_z(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """v3 smoother: mean of occupied CNPUs in the 7x7 window, only at occupied
    CNPUs; z-score over occupied CNPUs with the sample SD (statistics.stdev)."""
    occ = np.isfinite(values)
    total = window_sum(np.where(occ, values, 0.0))
    count = window_sum(occ.astype(float))
    smoothed = np.where(occ, total / np.maximum(count, 1), np.nan)
    obs = smoothed[occ]
    sd = obs.std(ddof=1) if obs.size > 1 else 0.0
    z = (smoothed - obs.mean()) / sd if sd else np.full_like(smoothed, np.nan)
    return smoothed, np.where(occ, z, np.nan)


def to_grid(df: pd.DataFrame, col: str) -> np.ndarray:
    a = np.full((NX, NY), np.nan)
    a[df["CNPU_X"].to_numpy(int), df["CNPU_Y"].to_numpy(int)] = df[col].to_numpy(float)
    return a


def self_check() -> dict:
    """Recompute the production TotalBandPower PooledSmoothedZ; must match mapping_v3.
    Uses the raw mapped rows (production does not dedupe repeated boundary IDs)."""
    cells = common.load_mapped_cells(dedupe=False)
    cells = cells[cells["InGrid"]]
    out = {}
    for g in common.GROUPS:
        c = cells[cells["Group"] == g]
        pooled = c.groupby(["CNPU_X", "CNPU_Y"])["TotalBandPower"].mean().reset_index()
        _, z = smooth_z(to_grid(pooled, "TotalBandPower"))
        prod = pd.read_csv(common.MAPPING_ROOT / g / "Group_Cnpu_Spectral_Summary.csv")
        prod = prod[prod["Feature"] == "TotalBandPower"]
        ref = to_grid(prod, "PooledSmoothedZ")
        same_mask = bool(np.array_equal(np.isfinite(ref), np.isfinite(z)))
        ok = same_mask and bool(np.allclose(ref[np.isfinite(ref)], z[np.isfinite(z)],
                                            rtol=1e-9, atol=1e-9))
        out[g] = {"occupied_cnpus": int(np.isfinite(z).sum()), "matches_mapping_v3": ok}
        if not ok:
            raise RuntimeError(f"{g}: recomputed TotalBandPower PooledSmoothedZ differs from "
                               "mapping_v3; smoothing is not identical to production.")
    return out


def group_table(c: pd.DataFrame, n_group_embryos: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    emb = (c.groupby(["EmbryoID", "CNPU_X", "CNPU_Y"])["Rhythmic"]
           .agg(NCells="size", NRhythmic="sum").reset_index())
    emb["Fraction"] = emb["NRhythmic"] / emb["NCells"]
    cnpu = (emb.groupby(["CNPU_X", "CNPU_Y"])
            .agg(NCells=("NCells", "sum"), NRhythmic=("NRhythmic", "sum"),
                 NEmbryos=("EmbryoID", "size"),
                 NEmbryosWithRhythmic=("NRhythmic", lambda s: int((s > 0).sum())),
                 EmbryoBalancedFraction=("Fraction", "mean"))
            .reset_index())
    cnpu["PooledFraction"] = cnpu["NRhythmic"] / cnpu["NCells"]
    cnpu["NGroupEmbryos"] = n_group_embryos
    cnpu["CoveragePct"] = 100 * cnpu["NEmbryos"] / n_group_embryos

    grids = {}
    for src, pre in (("PooledFraction", "Pooled"), ("EmbryoBalancedFraction", "EmbryoBalanced")):
        s, z = smooth_z(to_grid(cnpu, src))
        grids[f"{pre}Smoothed"], grids[f"{pre}SmoothedZ"] = s, z
    # Support QC over the same 7x7 window.
    per_emb = []
    for _, d in emb[emb["NRhythmic"] > 0].groupby("EmbryoID"):
        a = np.zeros((NX, NY))
        a[d["CNPU_X"].to_numpy(int), d["CNPU_Y"].to_numpy(int)] = d["NRhythmic"].to_numpy(float)
        per_emb.append(window_sum(a))
    stack = np.stack(per_emb) if per_emb else np.zeros((1, NX, NY))
    tot = stack.sum(axis=0)
    grids["WindowRhythmicCells"] = tot
    grids["WindowNEmbryosWithRhythmic"] = (stack > 0.5).sum(axis=0).astype(float)
    with np.errstate(invalid="ignore", divide="ignore"):
        grids["WindowTopEmbryoShare"] = np.where(tot > 0, stack.max(axis=0) / tot, np.nan)
    ix, iy = cnpu["CNPU_X"].to_numpy(int), cnpu["CNPU_Y"].to_numpy(int)
    for k, a in grids.items():
        cnpu[k] = a[ix, iy]
    cnpu["PooledSmoothedPct"] = 100 * cnpu["PooledSmoothed"]
    cnpu["CNPU_ID"] = cnpu["CNPU_X"].astype(str) + "_" + cnpu["CNPU_Y"].astype(str)
    cnpu["CNPU_CenterX"] = common.ORIGIN_X + cnpu["CNPU_X"] * common.CNPU_SIZE + common.CNPU_SIZE / 2
    cnpu["CNPU_CenterY"] = common.ORIGIN_Y + cnpu["CNPU_Y"] * common.CNPU_SIZE + common.CNPU_SIZE / 2
    return cnpu, emb


def plot_2x2(table, metric, label, path, *, limits, cmap, extend, min_embryos, title, counts):
    mpl = common.plot_style()
    import matplotlib.pyplot as plt
    mpl.rcParams.update({"axes.grid": False})
    shown = table[(table["NEmbryos"] >= min_embryos) & table[metric].notna()]
    fig, axes = plt.subplots(2, 2, figsize=(12, 10.5), constrained_layout=True)
    last = None
    for ax, g in zip(axes.flat, common.GROUPS):
        d = shown[shown["Group"] == g]
        last = ax.scatter(d["CNPU_CenterX"], d["CNPU_CenterY"], c=d[metric], s=6, marker="s",
                          cmap=cmap, vmin=limits[0], vmax=limits[1], linewidths=0)
        poly = common.load_template(common.SOX2[common.TEMPLATE_FOR_GROUP[g]])
        ax.plot(poly[:, 0], poly[:, 1], color=common.TEXT, linestyle=":", linewidth=1.5)
        n_emb, n_cells, n_rh = counts[g]
        ax.set(xlim=common.X_RANGE, ylim=common.Y_RANGE, xlabel="SuperEmbryo X",
               ylabel="SuperEmbryo Y")
        ax.set_title(f"{g.replace('_', ' ')}\n{n_emb} embryos · {n_cells:,} mapped cells · "
                     f"{n_rh:,} rhythmic", fontsize=11)
        ax.invert_yaxis()
        ax.set_aspect("equal")
        ax.text(0.02, 0.02, "anterior ↑", transform=ax.transAxes, fontsize=9,
                color=common.TEXT_MUTED)
    fig.colorbar(last, ax=axes, label=label, shrink=0.8, extend=extend)
    fig.suptitle(f"{title}\nCNPUs supported by ≥{min_embryos} embryos; dotted line = Sox2 "
                 "neural boundary", fontsize=12)
    fig.savefig(path, dpi=200)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--output-dir", type=Path,
                    default=common.ANALYSIS_ROOT / "step1_rhythmicity" / "rhythmic_fraction_maps")
    ap.add_argument("--min-embryos", type=int, default=2,
                    help="Display mask for the combined maps (v3 combined maps use 2).")
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()
    out = args.output_dir
    if out.exists() and any(out.iterdir()):
        if not args.overwrite:
            raise FileExistsError(f"{out} exists; pass --overwrite or choose a new --output-dir.")
        shutil.rmtree(out)
    out.mkdir(parents=True, exist_ok=True)

    cells = common.load_mapped_cells()                     # dedupes repeated boundary IDs
    n_all = len(cells)
    cells = cells[cells["InGrid"]].copy()
    rh = common.load_rhythmicity()[["Group", "EmbryoID", "CellID", "Rhythmic"]]
    cells = cells.merge(rh, on=["Group", "EmbryoID", "CellID"], how="left", validate="1:1")
    if cells["Rhythmic"].isna().any():
        raise RuntimeError(f"{int(cells['Rhythmic'].isna().sum())} mapped cells lack a Rhythmic flag")
    cells["Rhythmic"] = cells["Rhythmic"].astype(bool).astype(int)
    excluded = tuple(pd.read_csv(common.MAPPING_EXCLUSIONS)["EmbryoKey"].str.strip())
    if cells["EmbryoID"].str.startswith(excluded).any():
        raise RuntimeError("mapping-excluded embryos found in mapping_v3 cells")

    check = self_check()
    print("Self-check (TotalBandPower PooledSmoothedZ vs mapping_v3):", check)

    tables, emb_tables, counts = [], [], {}
    for g in common.GROUPS:
        c = cells[cells["Group"] == g]
        t, e = group_table(c, c["EmbryoID"].nunique())
        t.insert(0, "Group", g)
        e.insert(0, "Group", g)
        tables.append(t)
        emb_tables.append(e)
        counts[g] = (c["EmbryoID"].nunique(), len(c), int(c["Rhythmic"].sum()))
        assert t["NCells"].sum() == len(c) and t["NRhythmic"].sum() == c["Rhythmic"].sum()
    table = pd.concat(tables, ignore_index=True)
    lead = ["Group", "CNPU_ID", "CNPU_X", "CNPU_Y", "CNPU_CenterX", "CNPU_CenterY"]
    table = table[lead + [c for c in table.columns if c not in lead]]
    table.to_csv(out / "Group_Cnpu_RhythmicFraction.csv", index=False)
    pd.concat(emb_tables, ignore_index=True).to_csv(out / "Embryo_Cnpu_RhythmicFraction.csv", index=False)

    shown = table[table["NEmbryos"] >= args.min_embryos]
    pct_max = float(np.ceil(np.nanpercentile(shown["PooledSmoothedPct"], 99)))
    jet = "jet"
    seq = common.sequential_cmap()
    plot_2x2(table, "PooledSmoothedZ", "Smoothed rhythmic-fraction z-score (paper CNP method)",
             out / "RhythmicFraction_PooledSmoothedZ.png", limits=(-1, 1), cmap=jet, extend="both",
             min_embryos=args.min_embryos, counts=counts,
             title="Where are rhythmic cells concentrated? Fraction of mapped cells that are "
                   "rhythmic, smoothed z-score within each group")
    plot_2x2(table, "PooledSmoothedPct", "Smoothed % of cells rhythmic (shared scale)",
             out / "RhythmicFraction_PooledSmoothedPct.png", limits=(0, pct_max), cmap=seq,
             extend="max", min_embryos=args.min_embryos, counts=counts,
             title="Rhythmic cells as % of mapped cells (7×7-CNPU smoothed), one scale for all groups")
    plot_2x2(table, "EmbryoBalancedSmoothedZ", "Smoothed embryo-balanced rhythmic-fraction z-score",
             out / "RhythmicFraction_EmbryoBalancedSmoothedZ.png", limits=(-1, 1), cmap=jet,
             extend="both", min_embryos=args.min_embryos, counts=counts,
             title="Embryo-balanced version: mean of per-embryo rhythmic fractions, smoothed z-score")
    plot_2x2(table, "WindowNEmbryosWithRhythmic", "Embryos contributing rhythmic cells (7×7 window)",
             out / "RhythmicFraction_Support_NEmbryos.png", limits=(0, 10), cmap=seq, extend="max",
             min_embryos=args.min_embryos, counts=counts,
             title="Support QC: how many embryos contribute rhythmic cells to each CNPU's smoothing window")
    plot_2x2(table, "WindowTopEmbryoShare", "Share of window's rhythmic cells from the top embryo",
             out / "RhythmicFraction_Support_TopEmbryoShare.png", limits=(0, 1), cmap=seq,
             extend="neither", min_embryos=args.min_embryos, counts=counts,
             title="Support QC: 1.0 = all rhythmic cells near this CNPU come from one embryo")

    settings = {
        "created": datetime.now().isoformat(timespec="seconds"),
        "rhythmic_definition": "Q_global < 0.05 (P_global, BH within embryo)",
        "inputs": {"mapping": str(common.MAPPING_ROOT),
                   "rhythmicity": str(common.ANALYSIS_ROOT / "step1_rhythmicity" / common.RHYTHM_FILE)},
        "mapped_cells_dedup": n_all, "cells_in_grid": len(cells),
        "per_group_embryos_cells_rhythmic": counts, "min_embryos_display": args.min_embryos,
        "shared_pct_scale_max": pct_max, "self_check": check,
    }
    (out / "run_settings.json").write_text(json.dumps(settings, indent=2, default=str))
    print(json.dumps({k: settings[k] for k in ("cells_in_grid", "per_group_embryos_cells_rhythmic",
                                               "shared_pct_scale_max")}, default=str))
    print(f"Outputs: {out}")


if __name__ == "__main__":
    main()
