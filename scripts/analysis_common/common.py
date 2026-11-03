#!/usr/bin/env python3
"""Shared paths, constants and helpers for the post-mapping analysis steps.

Imported by scripts/step1_* ... scripts/step6_*. Every step writes under
ANALYSIS_ROOT/<step name>/ and never modifies spectral, table or mapping outputs.

Environment overrides (all optional):
  CALCIUM_SPECTRAL_ROOT  default data10/HPC_spectral_v3
  CALCIUM_TABLE_DIR      default CALCIUM/data/master_tables/group_tables_v3
  CALCIUM_MAPPING_ROOT   default scr10/CALCIUM/data/mapping_v3
  CALCIUM_ANALYSIS_ROOT  default scr10/CALCIUM/data/analysis_v3
"""
from __future__ import annotations

import importlib.util
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HOME = Path("/sciclone/home/cabare01")
REPO = HOME / "CALCIUM"
SCRIPTS = REPO / "scripts"
sys.path.insert(0, str(SCRIPTS))
import embryo_ref  # noqa: E402

SPECTRAL_ROOT = Path(os.environ.get("CALCIUM_SPECTRAL_ROOT", HOME / "data10/HPC_spectral_v3"))
TABLE_DIR = Path(os.environ.get("CALCIUM_TABLE_DIR", REPO / "data/master_tables/group_tables_v3"))
MAPPING_ROOT = Path(os.environ.get("CALCIUM_MAPPING_ROOT", HOME / "scr10/CALCIUM/data/mapping_v3"))
ANALYSIS_ROOT = Path(os.environ.get("CALCIUM_ANALYSIS_ROOT", HOME / "scr10/CALCIUM/data/analysis_v3"))
BOUNDARY_ROOT = HOME / "data10/boundary_info/Boundaries"
SOX2_DIR = Path("/sciclone/data10/cabare01/boundary_info/SuperEmbryoBoundaries")
SOX2 = {"stage14": SOX2_DIR / "1_new_Stage_14_Sox2_Boundary.csv",
        "stage15": SOX2_DIR / "1_new_Stage_15_Sox2_Boundary.csv"}
MAPPING_EXCLUSIONS = SCRIPTS / "mapping" / "mapping_exclusions.csv"

GROUPS = ("S13", "S13_devitellinated", "S14", "S14_devitellinated")
GROUP_FOLDER = {"S13": "stage_13", "S13_devitellinated": "stage_13_devitellinated",
                "S14": "stage_14", "S14_devitellinated": "stage_14_devitellinated"}
GROUP_STAGE = {"S13": 13, "S13_devitellinated": 13, "S14": 14, "S14_devitellinated": 14}
GROUP_MEMBRANE = {"S13": "intact", "S13_devitellinated": "devitellinated",
                  "S14": "intact", "S14_devitellinated": "devitellinated"}
TEMPLATE_FOR_GROUP = {"S13": "stage14", "S13_devitellinated": "stage14",
                      "S14": "stage15", "S14_devitellinated": "stage15"}
FEATURES = ("DominantFreqHz", "SpectralCentroid", "SpectralSpread",
            "SpectralEntropy", "TotalBandPower", "PeakPowerFraction")
# Heavy-tailed features are analysed on a log10 scale in parametric models.
LOG_FEATURES = {"TotalBandPower"}

# Legacy CNPU grid (identical to scripts/mapping/spatial_spectral_mapping.py).
CNPU_SIZE = 9
ORIGIN_X, ORIGIN_Y = 120, 1
X_RANGE, Y_RANGE = (120, 1122), (1, 1190)
MAX_CNPU_X = int(np.ceil((X_RANGE[1] - ORIGIN_X) / CNPU_SIZE)) - 1
MAX_CNPU_Y = int(np.ceil((Y_RANGE[1] - ORIGIN_Y) / CNPU_SIZE)) - 1


# --------------------------------------------------------------------------- I/O
def out_dir(step: str, *parts: str) -> Path:
    path = ANALYSIS_ROOT.joinpath(step, *parts)
    path.mkdir(parents=True, exist_ok=True)
    return path


def load_module(path: Path, name: str):
    """Import a pipeline stage whose file name starts with a digit."""
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def safe_name(value: object) -> str:
    """Embryo folder name used by scripts/wavelets/pipeline.py."""
    return re.sub(r"[^A-Za-z0-9_.-]", "_", str(value)).strip("._") or "unnamed"


def add_design(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["Stage"] = df["Group"].map(GROUP_STAGE)
    df["Membrane"] = df["Group"].map(GROUP_MEMBRANE)
    return df


def load_master(groups=GROUPS) -> pd.DataFrame:
    """Cell-level spectral tables (all kept cells, mapped or not)."""
    files = {g: TABLE_DIR / f"{GROUP_FOLDER[g]}_CW.csv" for g in groups}
    df = pd.concat([pd.read_csv(p, low_memory=False) for p in files.values()],
                   ignore_index=True)
    df["CellID"] = df["CellID"].astype(str)
    df["EmbryoID"] = df["EmbryoID"].astype(str).str.strip()
    return df


def load_mapped_cells(groups=GROUPS, dedupe: bool = True) -> pd.DataFrame:
    """Mapped cells with template coordinates and CNPU indices.

    The mapping mirrors MATLAB's innerjoin, so a Calcium_ID repeated in a
    boundary file yields more than one row (9 rows in 200217SP_E3). For
    statistics each cell must count once: dedupe keeps the first row.
    """
    frames = [pd.read_csv(MAPPING_ROOT / g / "All_MappedCells_with_CNPU.csv", low_memory=False)
              for g in groups]
    df = pd.concat(frames, ignore_index=True)
    df["CellID"] = df["CellID"].astype(str)
    df["EmbryoID"] = df["EmbryoID"].astype(str).str.strip()
    if dedupe:
        dup = df.duplicated(["Group", "EmbryoID", "CellID"])
        if dup.any():
            print(f"load_mapped_cells: dropped {int(dup.sum())} repeated boundary-ID rows")
        df = df[~dup].copy()
    df["InGrid"] = df["CNPU_X"].between(0, MAX_CNPU_X) & df["CNPU_Y"].between(0, MAX_CNPU_Y)
    return add_design(df)


def load_template(path: Path) -> np.ndarray:
    table = pd.read_csv(path)
    lower = {c.lower(): c for c in table.columns}
    xcol = lower.get("x", table.columns[0])
    ycol = lower.get("y", table.columns[1])
    pts = table[[xcol, ycol]].apply(pd.to_numeric, errors="coerce").dropna().to_numpy(float)
    if not np.allclose(pts[0], pts[-1]):
        pts = np.vstack([pts, pts[0]])
    return pts


def assign_region(df: pd.DataFrame, x="Centroid_1", y="Centroid_2") -> pd.Series:
    """neural / non-neural from the stage-matched Sox2 polygon (as in combine_group_maps.py)."""
    from matplotlib.path import Path as MplPath
    region = pd.Series("non-neural", index=df.index)
    for group, template in TEMPLATE_FOR_GROUP.items():
        mask = df["Group"] == group
        if mask.any():
            inside = MplPath(load_template(SOX2[template])).contains_points(
                df.loc[mask, [x, y]].to_numpy(float))
            region[mask] = np.where(inside, "neural", "non-neural")
    return region


RHYTHM_FILE = "rhythmicity_all_groups.csv"


def load_rhythmicity() -> pd.DataFrame:
    path = ANALYSIS_ROOT / "step1_rhythmicity" / RHYTHM_FILE
    if not path.exists():
        raise FileNotFoundError(f"{path} not found: run step1 build_rhythmic_tables.py first")
    df = pd.read_csv(path, low_memory=False)
    df["CellID"] = df["CellID"].astype(str)
    return df


def filter_cells(df: pd.DataFrame, which: str) -> pd.DataFrame:
    """which = all | rhythmic | nonrhythmic (needs step 1 output for the latter two)."""
    if which == "all":
        return df
    r = load_rhythmicity()[["Group", "EmbryoID", "CellID", "Rhythmic"]]
    merged = df.merge(r, on=["Group", "EmbryoID", "CellID"], how="left", validate="m:1")
    missing = merged["Rhythmic"].isna().sum()
    if missing:
        print(f"filter_cells: {missing} cells have no rhythmicity result and are dropped")
    keep = merged["Rhythmic"].eq(which == "rhythmic") & merged["Rhythmic"].notna()
    return merged[keep].drop(columns="Rhythmic")


# -------------------------------------------------------------------- statistics
def bh(p) -> np.ndarray:
    """Benjamini-Hochberg q-values; NaN stays NaN."""
    p = np.asarray(p, float)
    q = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if ok.sum():
        pv = p[ok]
        order = np.argsort(pv)
        ranked = pv[order] * ok.sum() / np.arange(1, ok.sum() + 1)
        ranked = np.minimum.accumulate(ranked[::-1])[::-1]
        out = np.empty_like(pv)
        out[order] = np.minimum(ranked, 1.0)
        q[ok] = out
    return q


def holm(p) -> np.ndarray:
    p = np.asarray(p, float)
    q = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    m = ok.sum()
    if m:
        pv = p[ok]
        order = np.argsort(pv)
        adj = np.maximum.accumulate((m - np.arange(m)) * pv[order])
        out = np.empty_like(pv)
        out[order] = np.minimum(adj, 1.0)
        q[ok] = out
    return q


def hedges_g(a, b) -> float:
    a, b = np.asarray(a, float), np.asarray(b, float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    na, nb = len(a), len(b)
    if na < 2 or nb < 2:
        return np.nan
    sp = np.sqrt(((na - 1) * a.var(ddof=1) + (nb - 1) * b.var(ddof=1)) / (na + nb - 2))
    if sp == 0:
        return np.nan
    return (a.mean() - b.mean()) / sp * (1 - 3 / (4 * (na + nb) - 9))


def cohens_d(a, b) -> float:
    """Pooled-SD Cohen's d, identical to combine_group_maps.py (paper rule)."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 2 or len(b) < 2:
        return np.nan
    pooled = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1))
                     / (len(a) + len(b) - 2))
    return (a.mean() - b.mean()) / pooled if pooled > 0 else np.nan


def _ols(X, y):
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    dof = len(y) - X.shape[1]
    sigma2 = resid @ resid / dof if dof > 0 else np.nan
    se = np.sqrt(np.diag(sigma2 * np.linalg.pinv(X.T @ X)))
    return beta, se, X @ beta, resid


def two_by_two(values, stage, membrane, n_perm=10000, seed=1) -> list[dict]:
    """Stage x membrane effects on embryo-level values (embryo = unit).

    Effect coding (stage 13 = -0.5, 14 = +0.5; intact = -0.5, devitellinated
    = +0.5) so that each coefficient is in the outcome's units:
      Stage       = mean(S14 groups) - mean(S13 groups), averaged over membrane
      Membrane    = mean(devitellinated) - mean(intact), averaged over stage
      Interaction = (S14dev - S14) - (S13dev - S13)
    p-values: Freedman-Lane residual permutation of the t statistic, valid for
    unbalanced groups and non-normal residuals.
    """
    y = np.asarray(values, float)
    s = np.where(np.asarray(stage) == 14, 0.5, -0.5)
    m = np.where(np.asarray(membrane) == "devitellinated", 0.5, -0.5)
    ok = np.isfinite(y)
    y, s, m = y[ok], s[ok], m[ok]
    # With +-0.5 coding the raw product s*m (+-0.25) gives a coefficient equal
    # to the difference of differences, and b1/b2 equal the averaged main effects.
    X = np.column_stack([np.ones_like(y), s, m, s * m])
    rng = np.random.default_rng(seed)
    beta, se, _, _ = _ols(X, y)
    out = []
    for j, name in ((1, "Stage"), (2, "Membrane"), (3, "Interaction")):
        t_obs = beta[j] / se[j]
        reduced = np.delete(X, j, axis=1)
        _, _, fitted_r, resid_r = _ols(reduced, y)
        exceed = 0
        for _ in range(n_perm):
            y_star = fitted_r + rng.permutation(resid_r)
            b, e, _, _ = _ols(X, y_star)
            exceed += abs(b[j] / e[j]) >= abs(t_obs)
        out.append({"Term": name, "Estimate": float(beta[j]), "t": float(t_obs),
                    "p_perm": (exceed + 1) / (n_perm + 1), "NEmbryos": int(len(y))})
    return out


# ------------------------------------------------------------------- CNPU maps
def cnpu_pooled_map(cells: pd.DataFrame, value: str, by=("Group",)) -> pd.DataFrame:
    """Paper CNP method (Methods 2.3) for any cell-level column.

    Pooled mean of all cells per CNPU, 7x7 (49-CNPU) mean over CNPUs that
    contain cells, z-scored over those CNPUs within each `by` level. Empty
    CNPUs are never filled (matches the 2026-11-02 smoother change in
    spatial_spectral_mapping.py). Out-of-grid cells are dropped.
    """
    d = cells[cells["InGrid"] & cells[value].notna()]
    agg = (d.groupby([*by, "CNPU_X", "CNPU_Y"])
           .agg(PooledMean=(value, "mean"), NCells=(value, "size"),
                NEmbryos=("EmbryoID", "nunique")).reset_index())
    pieces = []
    for _, g in agg.groupby(list(by)):
        lookup = {(x, y): v for x, y, v in zip(g.CNPU_X, g.CNPU_Y, g.PooledMean)}
        sm = [np.mean([lookup[(x + dx, y + dy)] for dx in range(-3, 4) for dy in range(-3, 4)
                       if (x + dx, y + dy) in lookup]) for x, y in zip(g.CNPU_X, g.CNPU_Y)]
        g = g.copy()
        g["PooledSmoothed"] = sm
        sd = g["PooledSmoothed"].std(ddof=1)
        g["PooledSmoothedZ"] = (g["PooledSmoothed"] - g["PooledSmoothed"].mean()) / sd if sd else np.nan
        pieces.append(g)
    out = pd.concat(pieces, ignore_index=True)
    out["CNPU_CenterX"] = ORIGIN_X + out.CNPU_X * CNPU_SIZE + CNPU_SIZE / 2
    out["CNPU_CenterY"] = ORIGIN_Y + out.CNPU_Y * CNPU_SIZE + CNPU_SIZE / 2
    return out


# ------------------------------------------------------------------------ style
TEXT = "#0b0b0b"
TEXT_MUTED = "#52514e"
GRID = "#e4e3df"
MEMBRANE_COLOR = {"intact": "#2a78d6", "devitellinated": "#eb6834"}
BLUE_RAMP = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]


def plot_style():
    import matplotlib as mpl
    mpl.use("Agg")
    mpl.rcParams.update({
        "font.size": 11, "axes.titlesize": 12, "axes.labelsize": 11,
        "axes.edgecolor": TEXT_MUTED, "axes.labelcolor": TEXT, "text.color": TEXT,
        "xtick.color": TEXT_MUTED, "ytick.color": TEXT_MUTED,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
        "axes.axisbelow": True, "legend.frameon": False,
        "savefig.dpi": 300, "savefig.bbox": "tight", "pdf.fonttype": 42,
    })
    return mpl


def sequential_cmap():
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list("blue_seq", BLUE_RAMP)
