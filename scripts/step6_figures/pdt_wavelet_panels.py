#!/usr/bin/env python3
"""
SfN poster panels 2 and 3: peak detection setups (PDTs) vs wavelet features.

Inputs (all per cell, joined by embryo + cell ID):
  * wavelet features (v3 cell-level group tables, CALCIUM/data/group_tables/stage_*_CW.csv)
  * the red-noise result (step 1, rhythmicity_all_groups.csv; Rhythmic = Q_global < 0.05)
  * Sudip's per-hour event counts for every PDT setup:
      Stage_<13|14>[_devitelined]_<36|54>_Conditions_Wave_<CW|EW|P>[_fpf]_Unique_LastFrameIDs_Perhour_Counts.csv
    one row per cell, columns Total_<code> (events per cell for setup <code>),
    plus CY5_Intensity (Sox2) for each cell
  * AllStage_Change_in_Mag_Data.csv: magnitude of every local peak (read in chunks)

Group membership comes only from data10/Embryo_Ref_Table.csv (via embryo_ref.py).
Embryo names are normalised with embryo_ref.embryo_key plus the six explicit PDT
aliases from the 2026-11-03 PDT_data audit.

Panel 2  Correlation of each PDT setup's per-cell event count with each wavelet
         feature and with Sox2, WITHIN each embryo, averaged across embryos.
           panel2_family_summary   compact, annotated (poster)
           panel2_all_setups       every setup, one row each (supplement / zoom)
           supp_pdt_x_pdt          how much the setups agree with each other
Panel 3  A  Cells above red noise vs the rest: event counts (each setup) and magnitude
         B  Membrane removal: devitellinated / intact fold for every setup's event
            count, against band power and event magnitude

Usage on the HPC (conda env wavelet-analysis; see pdt_wavelet_panels.sbatch):
    python pdt_wavelet_panels.py --inspect          # check files, IDs, matching; no figures
    python pdt_wavelet_panels.py                    # default --out in analysis_v3/step6_figures
    python pdt_wavelet_panels.py --skip-mag         # skip the 2.9 GB magnitude file

Original draft generated in a claude.ai chat (tested on mock data only); data
loading adapted to the cluster layout on 2026-11-08.
"""
import argparse
import glob
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy import stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from embryo_ref import GROUP_COLUMNS, embryo_key, load_reference  # noqa: E402

# ============================================================================
# CONFIG
# ============================================================================
HOME = "/sciclone/home/cabare01"
PDT_DIR = os.environ.get("PDT_DIR", "/sciclone/data10/cabare01/PDT_data")
ANALYSIS = os.path.join(HOME, "scr10/CALCIUM/data/analysis_v3")
CONFIG = {
    # Group membership (definitive)
    "ref_table": os.path.join(HOME, "data10/Embryo_Ref_Table.csv"),

    # Per-cell wavelet features (v3 group tables; *_CW = cell-level rows, *_EW = embryo summaries)
    "spectral_glob": os.path.join(HOME, "CALCIUM/data/group_tables/stage_*_CW.csv"),
    "spec_embryo_col": "EmbryoID",
    "spec_cell_col": "CellID",

    # Per-cell red-noise result (step 1). None -> panel 3A skipped
    "rednoise_glob": os.path.join(ANALYSIS, "step1_rhythmicity/rhythmicity_all_groups.csv"),
    "rn_embryo_col": "EmbryoID",
    "rn_cell_col": "CellID",
    "rn_pass_col": "Rhythmic",           # True = exceeds red noise (Q_global < 0.05)

    # PDT per-hour count files
    "pdt_count_glob": os.path.join(PDT_DIR, "Stage_*_Conditions_Wave_*_Unique_LastFrameIDs_Perhour_Counts.csv"),
    # Which PDT column holds the cell ID that matches the wavelet cell ID.
    # "auto" tries ID and Calcium_ID and keeps the one that matches more cells.
    "pdt_cell_col": "auto",
    # Per-setup measure: Total_<code> (events) or CellWiseActiveIncidence_<code>
    "pdt_measure_prefix": "Total_",

    # Peak magnitudes (every local peak). None -> skipped
    "mag_file": os.path.join(PDT_DIR, "AllStage_Change_in_Mag_Data.csv"),
    "mag_cell_col": "ID",
    "mag_value_col": "Peak_Mag_absDiff",

    "min_cells": 30,                     # joined cells for an embryo to enter correlations
    "min_cells_per_class": 5,            # cells in EACH red-noise class for panel 3A
}
DEFAULT_OUT = os.path.join(ANALYSIS, "step6_figures/pdt_wavelet_panels")

# PDT names that embryo_key cannot resolve (PDT_data audit 2026-11-03).
# No date + E-number fallback: 190918SP_E1 would collide with 190918RN_E1.
PDT_ALIASES = {
    "180509P_E7_": "180509SP_E7",
    "180509P_E6_": "180509SP_E6",
    "191030SPRN_": "191030SPRN_E3",
    "191020RNSP_": "191020RNSP_E6",
    "200313SP-E1": "200313SP_E1",
    "200313SP-E2": "200313SP_E2",
}

WAVELET_FEATURES = [
    ("DominantFreqHz",    "Dominant frequency",  False),
    ("SpectralCentroid",  "Spectral centroid",   False),
    ("SpectralSpread",    "Spectral spread",     False),
    ("SpectralEntropy",   "Spectral entropy",    False),
    ("TotalBandPower",    "Band power",          True),
    ("PeakPowerFraction", "Peak power fraction", True),
]
SOX2 = ("Sox2", "Sox2 level", False)
MAG = ("EventMag", "Event magnitude", True)

GROUP_ORDER = ["S13", "S13_devitellinated", "S14", "S14_devitellinated"]
FAMILY_ORDER = ["36-CW", "36-EW", "36-P", "54-CW", "54-EW", "54-P"]
INTACT_COLOR = "#2a6fdb"
DEVIT_COLOR = "#e8833a"
PDT_DOT = "#7d8796"
ACCENT = "#1f1f1f"

METRIC_RE = re.compile(r"(ActivePeaksPerHr|Total|CellWiseActiveIncidence|OnlySpikes|OnlyRosettes|OnlyWaves)_\d+$")
FILE_RE = re.compile(r"Stage_(13|14)(_devitelined)?_(36|54)_Conditions_Wave_(CW|EW|P)(_fpf)?_Unique")

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 11, "axes.titlesize": 12,
    "axes.labelsize": 11, "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": "#888888", "axes.linewidth": 0.8, "xtick.color": "#444444",
    "ytick.color": "#444444", "svg.fonttype": "none", "pdf.fonttype": 42,
})


# ============================================================================
# Loading
# ============================================================================
def norm_embryo(x):
    """Any embryo / file name -> reference key, e.g. 180118SPE7_ and 180118SP_E7_st_13_Q -> 180118SP_E7."""
    s = str(x).strip()
    return PDT_ALIASES.get(s) or embryo_key(s)


def as_id(series):
    """Cell IDs as clean integer strings ('18', not '18.0'); unparseable IDs -> NaN (callers drop them)."""
    num = pd.to_numeric(series, errors="coerce")
    num = num.where(num % 1 == 0)
    return num.astype("Int64").astype(str).where(num.notna())


def drop_bad_ids(d, col, label):
    bad = d[col].isna()
    if bad.any():
        print(f"  [{label}] dropped {bad.sum()} rows with an unparseable cell ID")
    return d[~bad]


def coerce_counts(d, cols, label):
    """Numeric setup columns. A row with any unparseable value (corrupt bytes in the 2022
    files shift or merge fields) has ALL its setup values in this file set to NaN."""
    num = d[cols].apply(pd.to_numeric, errors="coerce")
    bad = (num.isna() & d[cols].notna()).any(axis=1)
    if bad.any():
        who = d.loc[bad, ["Name", "ID"]].astype(str).agg(" cell ".join, axis=1).tolist()
        print(f"  [{label}] {bad.sum()} row(s) with corrupt values -> counts set to NaN: {who[:10]}")
        num.loc[bad] = np.nan
    out = d.drop(columns=cols)
    return pd.concat([out, num], axis=1)


def report_unmatched(raw_names, label):
    raw = pd.Series(pd.unique(pd.Series(raw_names).astype(str)))
    bad = raw[raw.map(norm_embryo).isna()]
    if len(bad):
        print(f"  [{label}] {len(bad)} names could not be normalised (dropped): {sorted(bad)[:10]}")


def need(pattern, label):
    files = sorted(glob.glob(pattern, recursive=True))
    if not files:
        sys.exit(f"[{label}] no files match: {pattern}\n  -> edit CONFIG")
    return files


def load_ref():
    assignment = load_reference(CONFIG["ref_table"])   # {key: ref column}
    ref = pd.DataFrame({"embryo": list(assignment), "column": list(assignment.values())})
    ref["group"] = ref["column"].map(lambda c: GROUP_COLUMNS[c][1])
    ref["stage"] = ref["column"].map(lambda c: GROUP_COLUMNS[c][2])
    ref["membrane"] = ref["column"].map(lambda c: GROUP_COLUMNS[c][3])
    ref["embryo_raw"] = ref["embryo"]
    return ref[["embryo", "embryo_raw", "group", "stage", "membrane"]]


def load_spectral():
    d = pd.concat([pd.read_csv(f) for f in need(CONFIG["spectral_glob"], "spectral")], ignore_index=True)
    d = d.rename(columns={CONFIG["spec_embryo_col"]: "embryo_raw", CONFIG["spec_cell_col"]: "cell"})
    report_unmatched(d["embryo_raw"], "spectral")
    d["embryo"] = d["embryo_raw"].map(norm_embryo)
    multi = d.groupby("embryo")["embryo_raw"].nunique()
    if (multi > 1).any():
        print(f"  [spectral] WARNING: several runs share a key: {multi[multi > 1].index.tolist()}")
    d = d.dropna(subset=["embryo"])
    d["cell"] = as_id(d["cell"])
    d = drop_bad_ids(d, "cell", "spectral")
    keep = ["embryo", "cell"] + [c for c, _, _ in WAVELET_FEATURES if c in d.columns]
    return d[keep].drop_duplicates(["embryo", "cell"])


def load_rednoise():
    pat = CONFIG["rednoise_glob"]
    files = sorted(glob.glob(pat, recursive=True)) if pat else []
    if not files:
        print(f"[red noise] no files ({pat}); panel 3A skipped")
        return None
    cols = [CONFIG["rn_embryo_col"], CONFIG["rn_cell_col"], CONFIG["rn_pass_col"]]
    d = pd.concat([pd.read_csv(f, usecols=cols) for f in files], ignore_index=True)
    d = d.rename(columns={CONFIG["rn_embryo_col"]: "embryo", CONFIG["rn_cell_col"]: "cell",
                          CONFIG["rn_pass_col"]: "exceeds_rn"})
    d["embryo"] = d["embryo"].map(norm_embryo)
    d = d.dropna(subset=["embryo"])
    d["cell"] = as_id(d["cell"])
    d = drop_bad_ids(d, "cell", "red noise")
    d["exceeds_rn"] = d["exceeds_rn"].astype(str).str.lower().isin(["1", "true", "t", "yes"])
    return d[["embryo", "cell", "exceeds_rn"]].drop_duplicates(["embryo", "cell"])


def load_pdt_counts():
    """Wide per-cell table: embryo, ID, Calcium_ID, Sox2, and one column per setup '<family>.<code>'."""
    prefix = CONFIG["pdt_measure_prefix"]
    wide, setups, base = None, [], None
    for f in need(CONFIG["pdt_count_glob"], "PDT counts"):
        m = FILE_RE.search(os.path.basename(f))
        if not m:
            print(f"  skipping unrecognised file name: {os.path.basename(f)}")
            continue
        fam = f"{m.group(3)}-{m.group(4)}"
        fgroup = f"S{m.group(1)}" + ("_devitellinated" if m.group(2) else "")
        # Read every per-setup metric column, not only <prefix>: a corrupt byte that eats a comma
        # anywhere in a row shifts all later fields, so the whole row must be checked.
        d = pd.read_csv(f, low_memory=False,
                        usecols=lambda c: c in ("Name", "ID", "Calcium_ID", "CY5_Intensity")
                        or bool(METRIC_RE.match(c)))
        report_unmatched(d["Name"], os.path.basename(f))
        d["embryo"] = d["Name"].map(norm_embryo)
        d = d.dropna(subset=["embryo"])
        metric = [c for c in d.columns if METRIC_RE.match(c)]
        d = coerce_counts(d, metric, os.path.basename(f))
        d = d.drop(columns=[c for c in metric if not c.startswith(prefix)])
        d["CY5_Intensity"] = pd.to_numeric(d["CY5_Intensity"], errors="coerce")
        d["ID"] = as_id(d["ID"])
        d = drop_bad_ids(d, "ID", os.path.basename(f))
        d["Calcium_ID"] = as_id(d["Calcium_ID"]) if "Calcium_ID" in d else d["ID"]
        d["file_group"] = fgroup
        codes = [c[len(prefix):] for c in d.columns if c.startswith(prefix)]
        if fam == "54-EW":
            # code = d1 d3 d4 d5; baseline 1/3 with %-of-baseline (d4 = 2) has no EW form
            # (Flowchart 4.1). Those 12 columns hold values but are invalid -> 42 EW setups.
            bad = [c for c in codes if len(c) == 4 and c[1] in "13" and c[2] == "2"]
            d = d.drop(columns=[f"{prefix}{c}" for c in bad])
            codes = [c for c in codes if c not in bad]
        ren = {f"{prefix}{c}": f"{fam}.{c}" for c in codes}
        d = d.rename(columns=ren).drop_duplicates(["embryo", "ID"])
        setups += list(ren.values())
        b = d[["embryo", "ID", "Calcium_ID", "CY5_Intensity", "file_group"]]
        base = b if base is None else pd.concat([base, b]).drop_duplicates(["embryo", "ID"])
        part = d[["embryo", "ID"] + list(ren.values())]
        if wide is None:
            wide = part
        else:
            # same family appears in several group files -> stack rows; new family -> new columns
            shared = [c for c in ren.values() if c in wide.columns]
            if shared:
                wide = wide.set_index(["embryo", "ID"])
                part = part.set_index(["embryo", "ID"])
                wide = wide.combine_first(part).reset_index().copy()
            else:
                wide = wide.merge(part, on=["embryo", "ID"], how="outer")
        print(f"  {os.path.basename(f)}: {fam}, {len(codes)} setups, {d.embryo.nunique()} embryos, {len(d):,} cells")
    out = base.merge(wide, on=["embryo", "ID"], how="left").rename(columns={"CY5_Intensity": "Sox2"})
    setups = sorted(set(setups), key=lambda s: (FAMILY_ORDER.index(s.split(".")[0])
                                                if s.split(".")[0] in FAMILY_ORDER else 99, s))
    return out, setups


def load_magnitude():
    f = CONFIG["mag_file"]
    if not f or not os.path.exists(f):
        print(f"[magnitude] {f} not found; skipped")
        return None
    cc, vc = CONFIG["mag_cell_col"], CONFIG["mag_value_col"]
    sums = []
    for i, ch in enumerate(pd.read_csv(f, usecols=["Name", cc, vc], dtype={cc: str}, chunksize=5_000_000)):
        ch[vc] = pd.to_numeric(ch[vc], errors="coerce")
        sums.append(ch.groupby(["Name", cc])[vc].agg(["sum", "count"]))
        print(f"  magnitude chunk {i + 1} ({(i + 1) * 5}M rows)", flush=True)
    s = pd.concat(sums).groupby(level=[0, 1]).sum().reset_index()
    report_unmatched(s["Name"], "magnitude")
    s["embryo"] = s["Name"].map(norm_embryo)
    s = s.dropna(subset=["embryo"])
    s["mag_cell"] = as_id(s[cc])
    s = drop_bad_ids(s, "mag_cell", "magnitude")
    s["EventMag"] = s["sum"] / s["count"]
    # one embryo name form per key; guard against duplicate keys after aliasing
    s = s.groupby(["embryo", "mag_cell"], as_index=False)["EventMag"].mean()
    return s[["embryo", "mag_cell", "EventMag"]]


# ============================================================================
# Helpers
# ============================================================================
def log_transform(df, cols):
    out = df.copy()
    for c in cols:
        v = out[c].astype(float)
        pos = v[v > 0]
        eps = pos.quantile(0.01) / 2 if len(pos) else 1e-6
        out[c] = np.log10(v.clip(lower=0) + eps)
    return out


def bh(p):
    p = np.asarray(p, float)
    q = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    if ok.sum():
        pv = p[ok]
        o = np.argsort(pv)
        r = pv[o] * len(pv) / (np.arange(len(pv)) + 1)
        r = np.minimum.accumulate(r[::-1])[::-1]
        t = np.empty_like(pv)
        t[o] = np.clip(r, 0, 1)
        q[ok] = t
    return q


def save(fig, out, name):
    for ext in ("png", "pdf", "svg"):
        fig.savefig(os.path.join(out, f"{name}.{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def fam_of(setup):
    return setup.split(".")[0]


def within_embryo_r(data, xcols, ycols, min_cells):
    """Per-embryo Pearson r between every x and y column -> Fisher-z array (embryos, x, y)."""
    zs, embs = [], []
    allc = list(dict.fromkeys(xcols + ycols))
    for emb, g in data.groupby("embryo"):
        g = g[allc]
        if g.dropna(how="all").shape[0] < min_cells:
            continue
        c = g.corr(min_periods=min_cells)
        r = c.loc[xcols, ycols].to_numpy(float)
        zs.append(np.arctanh(np.clip(r, -0.999999, 0.999999)))
        embs.append(emb)
    return (np.stack(zs) if zs else None), embs


def summarise_z(Z):
    r = np.tanh(np.nanmean(Z, axis=0))
    p = np.full(r.shape, np.nan)
    same = np.full(r.shape, np.nan)
    for idx in np.ndindex(r.shape):
        z = Z[(slice(None),) + idx]
        z = z[~np.isnan(z)]
        if len(z) >= 6 and np.any(z != 0):
            p[idx] = stats.wilcoxon(z).pvalue
            same[idx] = np.mean(np.sign(z) == np.sign(np.median(z)))
    q = bh(p.ravel()).reshape(p.shape)
    n = np.sum(~np.isnan(Z), axis=0)
    return r, q, same, n


def legend_membrane(fig, inside_ax=None):
    from matplotlib.lines import Line2D
    h = [Line2D([], [], marker="o", ls="", ms=7, mfc=INTACT_COLOR, mec="white", label="Membrane intact"),
         Line2D([], [], marker="o", ls="", ms=7, mfc=DEVIT_COLOR, mec="white", label="Devitellinated"),
         Line2D([], [], color=ACCENT, lw=2, label="Median")]
    if inside_ax is not None:
        inside_ax.legend(handles=h, frameon=False, fontsize=10, loc="lower left")
    else:
        fig.legend(handles=h, frameon=False, fontsize=10, loc="upper center", bbox_to_anchor=(0.5, 0.0), ncol=3)


def fold_axis(ax, lim):
    ticks = [t for t in (-3, -2, -1, 0, 1, 2, 3) if abs(t) <= lim]
    ax.set_yticks(ticks)
    ax.set_yticklabels([f"{2.0 ** t:g}×" for t in ticks])
    ax.set_ylim(-lim, lim)
    ax.axhline(0, color="#999999", lw=0.8, zorder=1)


# ============================================================================
# Panel 2
# ============================================================================
def panel2(tdata, setups, out, lines):
    wcols = [c for c, _, _ in WAVELET_FEATURES if c in tdata]
    ycols = wcols + (["Sox2"] if "Sox2" in tdata else [])
    ylab = {c: l for c, l, _ in WAVELET_FEATURES + [SOX2]}
    Z, embs = within_embryo_r(tdata, setups, ycols, CONFIG["min_cells"])
    if Z is None:
        lines.append("Panel 2: no embryo had enough joined cells.")
        return None
    r, q, same, n = summarise_z(Z)
    long = []
    for i, s in enumerate(setups):
        for j, y in enumerate(ycols):
            long.append(dict(setup=s, family=fam_of(s), code=s.split(".", 1)[1], feature=ylab[y],
                             mean_r=r[i, j], q_bh=q[i, j], frac_same_sign=same[i, j], n_embryos=n[i, j]))
    long = pd.DataFrame(long)
    long.to_csv(os.path.join(out, "panel2_setup_x_feature_r.csv"), index=False)

    # --- 2a  family summary (poster) ---------------------------------------
    fams = [f for f in FAMILY_ORDER if any(fam_of(s) == f for s in setups)]
    S = np.array([[long[(long.family == f) & (long.feature == ylab[y])].mean_r.median() for y in ycols] for f in fams])
    LO = np.array([[long[(long.family == f) & (long.feature == ylab[y])].mean_r.quantile(.1) for y in ycols] for f in fams])
    HI = np.array([[long[(long.family == f) & (long.feature == ylab[y])].mean_r.quantile(.9) for y in ycols] for f in fams])
    def agree(v):   # share of setups with BH q < .05 AND the sign of the family median
        return ((np.sign(v.mean_r) == np.sign(v.mean_r.median())) & (v.q_bh < 0.05)).mean()
    SIG = np.array([[agree(long[(long.family == f) & (long.feature == ylab[y])]) for y in ycols] for f in fams])
    nset = [sum(fam_of(s) == f for s in setups) for f in fams]
    fig, ax = plt.subplots(figsize=(1.05 * len(ycols) + 2.6, 0.62 * len(fams) + 2.0))
    im = ax.imshow(S, cmap="RdBu_r", norm=TwoSlopeNorm(0, -1, 1), aspect="auto")
    for i in range(len(fams)):
        for j in range(len(ycols)):
            txt = f"{S[i, j]:.2f}".replace("0.", ".").replace("-.", "−.")
            ax.text(j, i - 0.12, txt, ha="center", va="center", fontsize=10.5,
                    color="white" if abs(S[i, j]) > 0.6 else "#222222",
                    fontweight="bold" if SIG[i, j] >= 0.9 else "normal")
            ax.text(j, i + 0.22, f"{LO[i, j]:.2f} to {HI[i, j]:.2f}".replace("0.", ".").replace("-.", "−."),
                    ha="center", va="center", fontsize=7, color="white" if abs(S[i, j]) > 0.6 else "#555555")
    ax.set_xticks(range(len(ycols)))
    ax.set_xticklabels([ylab[y] for y in ycols], rotation=35, ha="right")
    ax.set_yticks(range(len(fams)))
    ax.set_yticklabels([f"{f}  ({k} setups)" for f, k in zip(fams, nset)])
    ax.tick_params(length=0)
    for s_ in ax.spines.values():
        s_.set_visible(False)
    if "Sox2" in ycols:
        ax.axvline(len(wcols) - 0.5, color="#333333", lw=1.2)
    cb = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
    cb.set_label("Within-embryo Pearson r")
    cb.outline.set_visible(False)
    ax.set_title("Event counts from every peak-detection setup vs wavelet features", loc="left")
    ax.text(0.0, -0.40, f"Cell-level r within each embryo, averaged across embryos (n = {len(embs)}); "
             "value = median over setups in the family, small text = 10th–90th percentile.\n"
             "Bold = >= 90% of the family's setups agree in sign with BH q < 0.05 (Wilcoxon on Fisher z). "
             "Counts, band power and peak power fraction log-transformed.", fontsize=8.5, color="#555555",
            transform=ax.transAxes, va="top")
    save(fig, out, "panel2_family_summary")

    # --- 2b  every setup -----------------------------------------------------
    h = max(4.0, 0.045 * len(setups) + 1.5)
    fig, ax = plt.subplots(figsize=(0.8 * len(ycols) + 2.8, h))
    im = ax.imshow(r, cmap="RdBu_r", norm=TwoSlopeNorm(0, -1, 1), aspect="auto", interpolation="nearest")
    ax.set_xticks(range(len(ycols)))
    ax.set_xticklabels([ylab[y] for y in ycols], rotation=35, ha="right")
    starts = [i for i, s in enumerate(setups) if i == 0 or fam_of(s) != fam_of(setups[i - 1])]
    ends = starts[1:] + [len(setups)]
    ax.set_yticks([(a + b - 1) / 2 for a, b in zip(starts, ends)])
    ax.set_yticklabels([fam_of(setups[a]) for a in starts])
    for a in starts[1:]:
        ax.axhline(a - 0.5, color="white", lw=2)
    ax.tick_params(length=0)
    for s_ in ax.spines.values():
        s_.set_visible(False)
    cb = fig.colorbar(im, ax=ax, shrink=0.5, pad=0.02)
    cb.set_label("Within-embryo Pearson r")
    cb.outline.set_visible(False)
    ax.set_title(f"All {len(setups)} setups", loc="left")
    save(fig, out, "panel2_all_setups")

    # numbers
    lines.append(f"Panel 2: {len(embs)} embryos, {len(setups)} PDT setups.")
    for y in ycols:
        v = long[long.feature == ylab[y]]
        lines.append(f"   {ylab[y]:<20} median r over setups {v.mean_r.median():+.2f} "
                     f"(range {v.mean_r.min():+.2f} to {v.mean_r.max():+.2f}); "
                     f"{(v.q_bh < 0.05).mean():.0%} of setups consistent across embryos")
    for f in fams:
        v = long[(long.family == f) & (long.feature == "Band power")]
        if len(v):
            lines.append(f"   Band power vs {f}: median r {v.mean_r.median():+.2f}")
    return long


def supp_pdt_x_pdt(tdata, setups, out, lines):
    Z, embs = within_embryo_r(tdata, setups, setups, CONFIG["min_cells"])
    if Z is None:
        return
    r = np.tanh(np.nanmean(Z, axis=0))
    pd.DataFrame(r, index=setups, columns=setups).to_csv(os.path.join(out, "supp_pdt_x_pdt_r.csv"))
    fig, ax = plt.subplots(figsize=(7.5, 6.6))
    im = ax.imshow(r, cmap="RdBu_r", norm=TwoSlopeNorm(0, -1, 1), interpolation="nearest")
    starts = [i for i, s in enumerate(setups) if i == 0 or fam_of(s) != fam_of(setups[i - 1])]
    ends = starts[1:] + [len(setups)]
    mids = [(a + b - 1) / 2 for a, b in zip(starts, ends)]
    for a in starts[1:]:
        ax.axhline(a - 0.5, color="white", lw=1.5)
        ax.axvline(a - 0.5, color="white", lw=1.5)
    ax.set_xticks(mids)
    ax.set_xticklabels([fam_of(setups[a]) for a in starts], rotation=45, ha="right")
    ax.set_yticks(mids)
    ax.set_yticklabels([fam_of(setups[a]) for a in starts])
    ax.tick_params(length=0)
    for s_ in ax.spines.values():
        s_.set_visible(False)
    cb = fig.colorbar(im, ax=ax, shrink=0.7, pad=0.02)
    cb.set_label("Within-embryo Pearson r")
    cb.outline.set_visible(False)
    ax.set_title("How much do peak-detection setups agree with each other?", loc="left")
    save(fig, out, "supp_pdt_x_pdt")
    off = r[~np.eye(len(r), dtype=bool)]
    lines.append(f"Setup-vs-setup agreement: median r {np.nanmedian(off):.2f} "
                 f"(10th–90th percentile {np.nanpercentile(off, 10):.2f} to {np.nanpercentile(off, 90):.2f}).")


# ============================================================================
# Panel 3A: cells above red noise vs the rest
# ============================================================================
def panel3a(data, setups, out, lines):
    mincls = CONFIG["min_cells_per_class"]
    rows = []
    for emb, g in data.groupby("embryo"):
        a, b = g[g.exceeds_rn], g[~g.exceeds_rn]
        if len(a) < mincls or len(b) < mincls:
            continue
        row = dict(embryo=emb, membrane=g.membrane.iloc[0], group=g.group.iloc[0], n_exceed=len(a), n_not=len(b))
        for s in setups:                                   # counts: mean events per cell
            ma, mb = a[s].mean(), b[s].mean()
            row[s] = np.log2(ma / mb) if (ma > 0 and mb > 0) else np.nan
        if "EventMag" in g:
            ma, mb = a.EventMag.median(), b.EventMag.median()
            row["EventMag"] = np.log2(ma / mb) if (ma > 0 and mb > 0) else np.nan
        rows.append(row)
    t = pd.DataFrame(rows)
    if t.empty:
        lines.append("Panel 3A: no embryo had enough cells in both red-noise classes.")
        return
    t.to_csv(os.path.join(out, "panel3A_per_embryo_log2_fold.csv"), index=False)

    per_setup = pd.DataFrame([dict(setup=s, family=fam_of(s), median_log2=t[s].median(),
                                   p=stats.wilcoxon(t[s].dropna()).pvalue if t[s].notna().sum() >= 6 else np.nan)
                              for s in setups])
    per_setup["q_bh"] = bh(per_setup.p)
    per_setup.to_csv(os.path.join(out, "panel3A_per_setup.csv"), index=False)

    has_mag = "EventMag" in t and t.EventMag.notna().any()
    fams = [f for f in FAMILY_ORDER if f in set(per_setup.family)]
    fig, axes = plt.subplots(1, 2 if has_mag else 1, figsize=(9.0 if has_mag else 7.0, 3.8),
                             gridspec_kw={"width_ratios": [len(fams), 1.6]} if has_mag else None, sharey=True)
    axes = np.atleast_1d(axes)
    vals = np.concatenate([per_setup.median_log2.dropna().to_numpy()] +
                          ([t.EventMag.dropna().to_numpy()] if has_mag else []))
    lim = min(3, max(1, np.ceil(np.nanmax(np.abs(vals)) * 1.15 * 2) / 2)) if len(vals) else 1
    rng = np.random.default_rng(0)
    ax = axes[0]
    for i, f in enumerate(fams):
        v = per_setup.loc[per_setup.family == f]
        ax.scatter(i + rng.uniform(-0.2, 0.2, len(v)), v.median_log2, s=14, c=PDT_DOT,
                   alpha=0.85, edgecolor="none", zorder=3)
        ax.plot([i - 0.3, i + 0.3], [v.median_log2.median()] * 2, color=ACCENT, lw=2, zorder=4)
    ax.set_xticks(range(len(fams)))
    ax.set_xticklabels(fams)
    ax.set_xlabel("Peak-detection family (one dot per setup)")
    fold_axis(ax, lim)
    ax.set_ylabel("Cells above red noise vs the rest\n(fold, log scale)")
    ax.set_title("A  Cells above red noise vs other cells: event counts", loc="left")
    if has_mag:
        ax2 = axes[1]
        v = t.EventMag.dropna()
        cols = [DEVIT_COLOR if m == "devitellinated" else INTACT_COLOR for m in t.loc[v.index, "membrane"]]
        ax2.scatter(rng.uniform(-0.2, 0.2, len(v)), v, s=20, c=cols, edgecolor="white", linewidth=0.6, zorder=3)
        ax2.plot([-0.3, 0.3], [v.median()] * 2, color=ACCENT, lw=2, zorder=4)
        ax2.set_xticks([0])
        ax2.set_xticklabels(["Event\nmagnitude"])
        ax2.set_xlim(-0.6, 0.6)
        ax2.set_xlabel("One dot per embryo")
        fold_axis(ax2, lim)
        ax2.set_title("Event size", loc="left")
        legend_membrane(fig)
    fig.tight_layout()
    save(fig, out, "panel3A_rednoise_vs_pdt")

    lines.append(f"Panel 3A: {len(t)} embryos with >= {mincls} cells in each class "
                 f"({(t.membrane == 'devitellinated').sum()} devitellinated, {(t.membrane == 'intact').sum()} intact).")
    lines.append(f"   Event counts: median fold over setups {2 ** per_setup.median_log2.median():.2f}x; "
                 f"{(per_setup.median_log2 > 0).mean():.0%} of setups > 1x; "
                 f"{(per_setup.q_bh < 0.05).mean():.0%} consistent across embryos (BH q < 0.05)")
    for f in fams:
        v = per_setup[per_setup.family == f]
        lines.append(f"     {f}: median {2 ** v.median_log2.median():.2f}x")
    if has_mag:
        v = t.EventMag.dropna()
        p = stats.wilcoxon(v).pvalue if len(v) >= 6 else np.nan
        lines.append(f"   Event magnitude: median {2 ** v.median():.2f}x (n={len(v)} embryos, Wilcoxon p={p:.2g})")


# ============================================================================
# Panel 3B: membrane effect, every setup vs band power and magnitude
# ============================================================================
def panel3b(data, setups, out, lines):
    cols = setups + [c for c in ("TotalBandPower", "EventMag") if c in data]
    agg = {s: "mean" for s in setups}
    agg.update({c: "median" for c in ("TotalBandPower", "EventMag") if c in data})
    emb = data.groupby(["embryo", "group", "stage", "membrane"]).agg(agg).reset_index()
    emb.to_csv(os.path.join(out, "panel3B_embryo_values.csv"), index=False)

    rows = []
    for stage in (13, 14):
        a = emb[(emb.stage == stage) & (emb.membrane == "intact")]
        b = emb[(emb.stage == stage) & (emb.membrane == "devitellinated")]
        for c in cols:
            x, y = a[c].dropna(), b[c].dropna()
            if len(x) >= 3 and len(y) >= 3 and x.median() > 0 and y.median() > 0:
                rows.append(dict(stage=stage, measure=c, family=fam_of(c) if c in setups else c,
                                 n_intact=len(x), n_devit=len(y),
                                 log2_fold=np.log2(y.median() / x.median()),
                                 p=stats.mannwhitneyu(x, y).pvalue))
    st = pd.DataFrame(rows)
    if st.empty:
        lines.append("Panel 3B: not enough embryos per group.")
        return
    st["q_bh"] = np.nan
    for stage in (13, 14):
        m = st.stage == stage
        st.loc[m, "q_bh"] = bh(st.loc[m, "p"])
    st.to_csv(os.path.join(out, "panel3B_stats.csv"), index=False)

    fams = [f for f in FAMILY_ORDER if f in set(st.family)]
    refs = [(c, l) for c, l in (("EventMag", "Event\nmagnitude"), ("TotalBandPower", "Band\npower")) if c in set(st.measure)]
    lim = min(3, max(1, np.ceil(np.nanmax(np.abs(st.log2_fold)) * 1.15 * 2) / 2))
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.9), sharey=True)
    rng = np.random.default_rng(1)
    for ax, stage in zip(axes, (13, 14)):
        s = st[st.stage == stage]
        for i, f in enumerate(fams):
            v = s[s.family == f]
            ax.scatter(i + rng.uniform(-0.2, 0.2, len(v)), v.log2_fold, s=14, c=PDT_DOT,
                       alpha=0.85, edgecolor="none", zorder=3)
            ax.plot([i - 0.3, i + 0.3], [v.log2_fold.median()] * 2, color=ACCENT, lw=2, zorder=4)
        for k, (c, l) in enumerate(refs):
            v = s[s.measure == c]
            x = len(fams) + 0.4 + k
            if len(v):
                ax.scatter([x], v.log2_fold, s=70, marker="D", c=DEVIT_COLOR, edgecolor=ACCENT, linewidth=0.8, zorder=4)
                ax.text(x, v.log2_fold.iloc[0] + lim * 0.09, f"{2 ** v.log2_fold.iloc[0]:.1f}×",
                        ha="center", fontsize=9)
        ax.axvline(len(fams) - 0.3, color="#cccccc", lw=0.8)
        ax.set_xticks(list(range(len(fams))) + [len(fams) + 0.4 + k for k in range(len(refs))])
        ax.set_xticklabels(fams + [l for _, l in refs], rotation=35, ha="right")
        fold_axis(ax, lim)
        ax.set_title(f"Stage {stage}", fontsize=11)
    axes[0].set_ylabel("Devitellinated vs intact\n(fold in embryo medians, log scale)")
    fig.suptitle("B  Membrane removal: event counts from every setup vs event size and wavelet power",
                 x=0.01, ha="left", fontsize=12)
    from matplotlib.lines import Line2D
    fig.legend(handles=[Line2D([], [], marker="o", ls="", ms=6, mfc=PDT_DOT, mec="none", label="One peak-detection setup (event count)"),
                        Line2D([], [], color=ACCENT, lw=2, label="Median of setups"),
                        Line2D([], [], marker="D", ls="", ms=8, mfc=DEVIT_COLOR, mec=ACCENT, label="Wavelet power / event magnitude")],
               frameon=False, fontsize=10, loc="upper center", bbox_to_anchor=(0.5, 0.0), ncol=3)
    fig.tight_layout()
    save(fig, out, "panel3B_membrane_all_setups")

    for stage in (13, 14):
        s = st[st.stage == stage]
        c = s[s.measure.isin(setups)]
        lines.append(f"Panel 3B, stage {stage}: event counts devitellinated/intact median {2 ** c.log2_fold.median():.2f}x "
                     f"over setups; {(c.log2_fold > 0).mean():.0%} of setups higher in devitellinated, "
                     f"{((c.log2_fold > 0) & (c.q_bh < 0.05)).mean():.0%} significantly higher, "
                     f"{((c.log2_fold < 0) & (c.q_bh < 0.05)).mean():.0%} significantly lower (BH q < 0.05)")
        for m, l in (("TotalBandPower", "Band power"), ("EventMag", "Event magnitude")):
            v = s[s.measure == m]
            if len(v):
                v = v.iloc[0]
                lines.append(f"   {l}: {2 ** v.log2_fold:.2f}x (n={v.n_intact}/{v.n_devit} embryos, "
                             f"Mann-Whitney p={v.p:.2g})")


# ============================================================================
# Build + main
# ============================================================================
def build(inspect=False, skip_mag=False):
    ref = load_ref()
    spec = load_spectral()
    print("Reading PDT count files ...", flush=True)
    pdt, setups = load_pdt_counts()
    rn = load_rednoise()
    mag = None if skip_mag else load_magnitude()

    # which PDT ID column matches the wavelet cell IDs?
    common = sorted(set(spec.embryo) & set(pdt.embryo))
    keys = set(zip(spec.embryo, spec.cell))
    rate = {}
    for col in ("ID", "Calcium_ID"):
        sub = pdt[pdt.embryo.isin(common)]
        rate[col] = np.mean([k in keys for k in zip(sub.embryo, sub[col])]) if len(sub) else 0.0
    col = CONFIG["pdt_cell_col"] if CONFIG["pdt_cell_col"] != "auto" else max(rate, key=rate.get)
    pdt["cell"] = pdt[col]

    # diagnostics (always printed so the job log records the join)
    print("\n== Setups ==")
    print(f"{len(setups)} setups:", pd.Series([fam_of(s) for s in setups]).value_counts().reindex(FAMILY_ORDER).dropna().astype(int).to_dict())
    print("example setup names:", setups[:4])
    print("\n== Embryos ==")
    print(f"reference {ref.embryo.nunique()} | spectral {spec.embryo.nunique()} | PDT {pdt.embryo.nunique()} | both spectral+PDT {len(common)}")
    print("PDT embryos not in reference:", sorted(set(pdt.embryo) - set(ref.embryo))[:20])
    print("reference embryos without PDT:", sorted(set(ref.embryo) - set(pdt.embryo))[:20])
    print("\n== Cell-ID matching (share of PDT cells found among wavelet cells) ==")
    for k, v in rate.items():
        print(f"  {k}: {v:.0%}")
    print(f"  -> using {col}")
    if rn is not None:
        print(f"red-noise cells matching wavelet cells: {np.mean([k in keys for k in zip(rn.embryo, rn.cell)]):.0%}")
    if mag is not None:
        mk = set(zip(pdt.embryo, pdt.cell))
        print(f"magnitude cells matching PDT cells: {np.mean([k in mk for k in zip(mag.embryo, mag.mag_cell)]):.0%}")

    data = spec.merge(pdt.drop(columns=["ID", "Calcium_ID"]), on=["embryo", "cell"], how="inner")
    if mag is not None:
        data = data.merge(mag.rename(columns={"mag_cell": "cell"}), on=["embryo", "cell"], how="left")
    if rn is not None:
        data = data.merge(rn, on=["embryo", "cell"], how="left")
        data["exceeds_rn"] = data["exceeds_rn"].astype("boolean").fillna(False).astype(bool)
    data = data.merge(ref, on="embryo", how="inner")
    mism = (data.drop_duplicates("embryo").pipe(lambda d: d[d.file_group != d.group]))
    if len(mism):
        print(f"note: {len(mism)} embryos sit in a PDT file for a different group than the reference table; "
              "the reference table is used.")

    print("\n== Joined ==")
    print(f"embryos {data.embryo.nunique()}, cells {len(data):,}")
    print(data.groupby("group").embryo.nunique().reindex(GROUP_ORDER).to_string())
    per_emb = (data.groupby("embryo").size() / spec.groupby("embryo").size()).dropna()
    if len(per_emb):
        print(f"per-embryo share of wavelet cells joined: median {per_emb.median():.0%}, "
              f"min {per_emb.min():.0%} ({per_emb.idxmin()})")
    print("(wavelet tables hold all QC-passed cells; PDT tables only last-frame cells, so <100% is expected)")
    if "EventMag" in data:
        print(f"joined cells with an event magnitude: {data.EventMag.notna().mean():.0%}")
    if data.embryo.nunique() == 0 or max(rate.values()) < 0.5:
        print("WARNING: fewer than half of PDT cells matched a wavelet cell -> check cell-ID columns in CONFIG.")
    return data, setups


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inspect", action="store_true")
    ap.add_argument("--skip-mag", action="store_true", help="skip the large magnitude file")
    ap.add_argument("--skip-supp", action="store_true", help="skip the setup-vs-setup matrix")
    ap.add_argument("--out", default=DEFAULT_OUT)
    args = ap.parse_args()

    data, setups = build(inspect=args.inspect, skip_mag=args.skip_mag)
    if args.inspect:
        return
    os.makedirs(args.out, exist_ok=True)
    lines = [f"Joined: {data.embryo.nunique()} embryos, {len(data):,} cells, {len(setups)} PDT setups."]

    logc = setups + [c for c, _, lg in WAVELET_FEATURES + [MAG] if lg and c in data]
    tdata = log_transform(data, logc)

    panel2(tdata, setups, args.out, lines)
    if not args.skip_supp:
        supp_pdt_x_pdt(tdata, setups, args.out, lines)
    if "exceeds_rn" in data:
        panel3a(data, setups, args.out, lines)
    panel3b(data, setups, args.out, lines)

    with open(os.path.join(args.out, "numbers_for_poster.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n" + "\n".join(lines))
    print(f"\nFigures and tables in {args.out}/")


if __name__ == "__main__":
    main()
