#!/usr/bin/env python3
"""SfN poster figures 1, 2A, 2B, 2C, 3A and 3B (poster draft of 2026-11-04).

  Fig 1   one cell: dF/F trace and wavelet power (cone of influence shaded),
          stacked on the left; time-averaged spectrum and its six features on
          the right. The cell is the Fig 2A "exceeds red noise" example.
  Fig 2A  two cells from one devitellinated embryo: observed spectrum vs AR(1)
          red-noise mean and 95th percentile ("Exceeds" / "Within red noise").
  Fig 2B  % of cells exceeding red noise per embryo (all tested cells).
  Fig 2C  smoothed % of mapped cells exceeding red noise per CNPU, one shared
          0-44% scale (from step1_rhythmicity/rhythmic_fraction_maps).
  Fig 3A  TotalBandPower smoothed pooled-mean z-maps (-1..1 jet, z-scored
          within each group), all mapped cells, from mapping_v3.
  Fig 3B  embryo medians (mapped cells) of total band power (log), spectral
          entropy (max 6.32) and peak power fraction.

"Exceeds red noise" = the step 1 Rhythmic flag (Q_global < 0.05). Poster text
never says rhythmic / oscillator / quasi-periodic or gives periods in seconds:
the 2026-11-04 transient check showed these cells carry irregular slow
transients, not a rhythm.

Poster style (author's spec, 2026-11-04): Arial, else Liberation Sans (bundled
in ./fonts, SIL OFL); text >= 24 pt at the placed size, per-panel n lines
>= 20 pt; intact #004E38, devitellinated #6AA84F, axes/secondary text #595959;
Hz, not mHz. Each figure is saved at exactly its placed size, 300 dpi,
transparent background, as PDF and PNG. A QA pass checks fonts, sizes, words,
clipping and overlaps before saving, and the saved files' size and alpha after.
figure_facts.json holds every number a caption might quote.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "analysis_common"))
import common  # noqa: E402

RT = common.load_module(common.SCRIPTS / "step1_rhythmicity" / "rhythmicity_test.py", "rt")
RHYTHM_ROOT = common.ANALYSIS_ROOT / "step1_rhythmicity"
FONT_DIR = Path(__file__).resolve().with_name("fonts")
DPI = 300
SIZES = {"Fig1": (14.6, 9.0), "Fig2A": (14.6, 5.0), "Fig2B": (14.6, 5.0), "Fig2C": (14.6, 6.0),
         "Fig3A": (14.6, 9.5), "Fig3B": (14.6, 6.8)}
INTACT, DEVIT, GREY, BLACK = "#004E38", "#6AA84F", "#595959", "#000000"
LIGHT = "#E3E3E3"  # red-noise band / spread band fill (not text)
MIN_PT, MIN_N_PT, N_PT, BASE_PT = 24, 20, 20, 24
FORBIDDEN = re.compile(r"rhythm|quasi|oscillat|mHz|\b\d+(\.\d+)?\s?s\b", re.I)
GROUP_ORDER = [("intact", 13), ("devitellinated", 13), ("intact", 14), ("devitellinated", 14)]
XPOS = {("intact", 13): 0.0, ("devitellinated", 13): 0.6, ("intact", 14): 1.55, ("devitellinated", 14): 2.15}
MEMB_COLOR = {"intact": INTACT, "devitellinated": DEVIT}
PANEL_TITLE = {"S13": "Stage 13 · intact", "S13_devitellinated": "Stage 13 · devitellinated",
               "S14": "Stage 14 · intact", "S14_devitellinated": "Stage 14 · devitellinated"}
FREQ_TICKS = [0.02, 0.04, 0.06, 0.08]
FONT = None


# ------------------------------------------------------------------ style + QA
def style() -> str:
    import matplotlib as mpl
    mpl.use("Agg")
    from matplotlib import font_manager as fm
    names = {f.name for f in fm.fontManager.ttflist}
    if "Arial" in names:
        family = "Arial"
    else:
        for f in sorted(FONT_DIR.glob("LiberationSans-*.ttf")):
            fm.fontManager.addfont(str(f))
        family = "Liberation Sans"
        fm.findfont(family, fallback_to_default=False)  # raises if missing
    mpl.rcParams.update({
        "font.family": family, "font.size": BASE_PT, "axes.titlesize": BASE_PT,
        "axes.labelsize": BASE_PT, "xtick.labelsize": BASE_PT, "ytick.labelsize": BASE_PT,
        "legend.fontsize": BASE_PT, "legend.title_fontsize": BASE_PT, "figure.titlesize": 26,
        "mathtext.fontset": "custom", "mathtext.rm": family, "mathtext.it": f"{family}:italic",
        "mathtext.bf": f"{family}:bold", "mathtext.sf": family, "mathtext.cal": family,
        "mathtext.default": "regular",
        "text.color": BLACK, "axes.titlecolor": BLACK, "axes.labelcolor": GREY,
        "axes.edgecolor": GREY, "xtick.color": GREY, "ytick.color": GREY,
        "axes.linewidth": 1.5, "xtick.major.width": 1.5, "ytick.major.width": 1.5,
        "xtick.major.size": 7, "ytick.major.size": 7, "xtick.minor.width": 1, "ytick.minor.width": 1,
        "axes.spines.top": False, "axes.spines.right": False, "axes.grid": False,
        "legend.frameon": False, "lines.linewidth": 3,
        "savefig.dpi": DPI, "savefig.bbox": "standard", "savefig.transparent": True,
        "pdf.fonttype": 42, "ps.fonttype": 42, "axes.unicode_minus": True,
    })
    return family


def new_figure(key: str):
    import matplotlib.pyplot as plt
    return plt.figure(figsize=SIZES[key], layout="constrained")


def n_text(t):
    """Mark a Text as a per-panel n line (allowed down to 20 pt)."""
    t.set_gid("nline")
    return t


def qa(fig, key: str) -> None:
    from matplotlib.text import Text
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    W, H = fig.bbox.width, fig.bbox.height
    # Tick labels that are never drawn: ticks outside the view interval, hidden
    # labels, and every tick of an axes that is switched off.
    undrawn = set()
    for ax in fig.axes:
        for axis in (ax.xaxis, ax.yaxis):
            lo, hi = sorted(axis.get_view_interval())
            tol = 1e-9 * max(abs(lo), abs(hi), 1)
            for tick in axis.get_major_ticks() + axis.get_minor_ticks():
                for lab in (tick.label1, tick.label2):
                    if (not ax.axison or not axis.get_visible() or not lab.get_visible()
                            or not (lo - tol <= tick.get_loc() <= hi + tol)):
                        undrawn.add(id(lab))
    problems, boxes = [], []
    for t in fig.findobj(Text):
        s = t.get_text()
        if not s.strip() or not t.get_visible() or id(t) in undrawn:
            continue
        try:
            bb = t.get_window_extent(renderer)
        except Exception:
            continue
        if bb.width == 0 and bb.height == 0:
            continue
        floor = MIN_N_PT if t.get_gid() == "nline" else MIN_PT
        if t.get_fontsize() < floor - 1e-6:
            problems.append(f"{s!r}: {t.get_fontsize():.1f} pt < {floor} pt")
        if t.get_fontname() != FONT:
            problems.append(f"{s!r}: font {t.get_fontname()}")
        plain = s.replace("$\\mathdefault{", "").replace("}$", "")
        if FORBIDDEN.search(plain):
            problems.append(f"{s!r}: forbidden wording")
        if bb.x0 < -1 or bb.y0 < -1 or bb.x1 > W + 1 or bb.y1 > H + 1:
            problems.append(f"{s!r}: outside the figure (clipped)")
        boxes.append((s, bb))
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            a, b = boxes[i][1], boxes[j][1]
            if a.x0 < b.x1 - 2 and b.x0 < a.x1 - 2 and a.y0 < b.y1 - 2 and b.y0 < a.y1 - 2:
                problems.append(f"overlap: {boxes[i][0]!r} / {boxes[j][0]!r}")
    if problems:
        raise RuntimeError(f"{key} QA failed:\n  " + "\n  ".join(problems))


def save(fig, out: Path, key: str, name: str) -> dict:
    import matplotlib.pyplot as plt
    from PIL import Image
    qa(fig, key)
    w, h = SIZES[key]
    for ext in ("pdf", "png"):
        fig.savefig(out / f"{name}.{ext}", dpi=DPI, transparent=True)
    plt.close(fig)
    img = Image.open(out / f"{name}.png")
    px = (round(w * DPI), round(h * DPI))
    alpha = np.asarray(img.convert("RGBA"))[..., 3]
    box = re.search(rb"/MediaBox\s*\[\s*0 0 ([\d.]+) ([\d.]+)", (out / f"{name}.pdf").read_bytes())
    pdf_in = (float(box.group(1)) / 72, float(box.group(2)) / 72) if box else None
    check = {"png_px": list(img.size), "expected_px": list(px), "png_dpi": [round(d) for d in img.info.get("dpi", (0, 0))],
             "transparent_px_pct": round(100 * float((alpha == 0).mean()), 1),
             "pdf_in": [round(v, 3) for v in pdf_in] if pdf_in else None}
    if tuple(img.size) != px or img.mode != "RGBA" or alpha[0, 0] != 0 or \
            pdf_in is None or not np.allclose(pdf_in, (w, h), atol=1e-3):
        raise RuntimeError(f"{name}: saved file check failed {check}")
    print(f"wrote {name}  {w}x{h} in  {img.size[0]}x{img.size[1]} px  OK")
    return check


def hz_axis(ax):
    ax.set_xticks(FREQ_TICKS, [f"{v:g}" for v in FREQ_TICKS])
    ax.set_xlabel("Frequency (Hz)")


def white_green():
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list("white_green", ["#FFFFFF", INTACT])


# ------------------------------------------------------------------ example cells
def pick_examples(r: pd.DataFrame, group: str) -> dict:
    """Embryo: the one in `group` whose % of cells exceeding red noise is
    closest to the group's embryo median, among embryos with >= 20 such
    single-segment cells. "Exceeds" cell: single segment, >= 55 min used, peak
    excess at 0.01-0.015 Hz, largest power ratio closest to the group median
    ratio of exceeding cells. "Within" cell: same embryo, single segment,
    >= 55 min, P_global closest to the group median P_global of cells within
    red noise."""
    g = r[r.Group == group]
    ok = (g.NSegments == 1) & (g.UsedDuration_min >= 55)
    pct = g.groupby("EmbryoID")["Rhythmic"].mean() * 100
    n_ok = g[ok & g.Rhythmic].groupby("EmbryoID").size()
    cand = pct[n_ok.reindex(pct.index).fillna(0) >= 20]
    emb = (cand - pct.median()).abs().idxmin()
    d = g[(g.EmbryoID == emb) & ok]
    rh = d[d.Rhythmic & d.MaxRatioFreqHz.between(0.01, 0.015)]
    target_ratio = g.loc[g.Rhythmic, "MaxRatio"].median()
    rcell = rh.iloc[(rh.MaxRatio - target_ratio).abs().argmin()]
    nr = d[~d.Rhythmic]
    target_p = g.loc[~g.Rhythmic, "P_global"].median()
    ncell = nr.iloc[(nr.P_global - target_p).abs().argmin()]
    return {"group": group, "embryo": emb, "embryo_pct_rhythmic": float(pct[emb]),
            "group_median_pct_rhythmic": float(pct.median()),
            "rhythmic": rcell, "nonrhythmic": ncell}


def cell_segments(group: str, embryo: str, cell_id: str):
    gdir = common.SPECTRAL_ROOT / common.GROUP_FOLDER[group]
    manifest = json.loads((gdir / "run_manifest.json").read_text())
    raw = {common.safe_name(Path(p).stem.removesuffix(" - Not_self_Merged")): p
           for p in manifest["input_files"]}[embryo]
    table = pd.read_csv(raw, low_memory=False)
    try:
        info = RT.preprocess.embryo_frame_times(table)
    except Exception:
        info = None
    if info is not None:
        table["_time_s"] = info["time"]
    rows = table[table["ID"].astype(str) == cell_id]
    cinfo = None if info is None else {**info, "time": rows["_time_s"].to_numpy()}
    axis = RT.spectral.frequency_axis()
    data = RT.preprocess.preprocess_data(rows, float(axis[0]), RT.spectral.coi_factor(), cinfo)
    return data["segments"], axis


# ------------------------------------------------------------------ Figure 1
def fig1(out, ex, master, facts, checks):
    import pywt
    from matplotlib.ticker import MaxNLocator
    group, emb, cell = ex["group"], ex["embryo"], str(ex["rhythmic"].CellID)
    segs, axis = cell_segments(group, emb, cell)
    lib = np.load(common.SPECTRAL_ROOT / common.GROUP_FOLDER[group] / emb / "SpectralLibrary.npz")
    i = list(lib["CellID"].astype(str)).index(cell)
    saved = np.maximum(lib["Power"][i], 0)
    again = RT.batch_mean_power([np.asarray(s["dFF"], float)[None, :] for s in segs],
                                [s["dt"] for s in segs], axis)[0]
    if not np.allclose(again, saved, rtol=1e-6, atol=1e-12):
        raise RuntimeError("Fig 1: recomputed spectrum differs from SpectralLibrary.npz")
    feats = master.set_index(["Group", "EmbryoID", "CellID"]).loc[(group, emb, cell)]

    fig = new_figure("Fig1")
    outer = fig.add_gridspec(1, 2, width_ratios=[1.12, 1])
    left = outer[0].subgridspec(2, 1, height_ratios=[1, 1.45])
    right = outer[1].subgridspec(2, 1, height_ratios=[1.5, 1])
    ax_a = fig.add_subplot(left[0])
    ax_b = fig.add_subplot(left[1], sharex=ax_a)
    ax_c = fig.add_subplot(right[0])
    ax_t = fig.add_subplot(right[1])

    wavelet = pywt.ContinuousWavelet(RT.spectral.WAVELET)
    center = pywt.central_frequency(wavelet)
    coi_s = RT.spectral.coi_factor() / axis
    t0 = segs[0]["time"][0]
    pw, vmax = [], 0.0
    for s in segs:
        scales = center / (axis * s["dt"])
        coef, _ = pywt.cwt(np.asarray(s["dFF"], float), scales, wavelet, sampling_period=s["dt"],
                           method="fft")
        p = np.abs(coef) ** 2 / scales[:, None]
        pw.append(p)
        vmax = max(vmax, float(np.percentile(p, 99.5)))
    mesh = None
    for s, p in zip(segs, pw):
        tm = (s["time"] - t0) / 60
        ax_a.plot(tm, s["dFF"], color=BLACK, lw=1.6)
        mesh = ax_b.pcolormesh(tm, axis, p, cmap=white_green(), vmin=0, vmax=vmax, shading="auto",
                               rasterized=True)
        el = s["time"] - s["time"][0]
        invalid = (np.minimum(el, el[-1] - el)[None, :] < coi_s[:, None]).astype(float)
        ax_b.contourf(tm, axis, invalid, levels=[0.5, 1.5], colors=GREY, alpha=0.35)
        ax_b.contour(tm, axis, invalid, levels=[0.5], colors=GREY, linewidths=2, linestyles="--")
    ax_a.set_ylabel("ΔF/F")
    ax_a.set_title("ΔF/F of one cell")
    ax_a.set_title("A", loc="left", fontweight="bold", fontsize=28)
    ax_a.yaxis.set_major_locator(MaxNLocator(3))
    ax_a.tick_params(labelbottom=False)
    ax_b.set(xlabel="Time (min)", ylabel="Frequency (Hz)", ylim=(axis[0], axis[-1]))
    ax_b.set_yticks(FREQ_TICKS, [f"{v:g}" for v in FREQ_TICKS])
    ax_b.set_xticks([0, 20, 40, 60])
    ax_b.set_title("Wavelet power")
    ax_b.set_title("B", loc="left", fontweight="bold", fontsize=28)
    cb = fig.colorbar(mesh, ax=ax_b, extend="max", pad=0.02, aspect=12)
    cb.set_label("Power", color=GREY)
    cb.locator = MaxNLocator(3)
    cb.update_ticks()
    cb.outline.set_edgecolor(GREY)

    dom, cen, spr = float(feats.DominantFreqHz), float(feats.SpectralCentroid), float(feats.SpectralSpread)
    ax_c.axvspan(cen - spr, cen + spr, color=LIGHT, lw=0, zorder=0)
    ax_c.axvline(cen, color=GREY, ls="--", lw=2.5)
    ax_c.plot(axis, saved, color=BLACK, lw=3, zorder=3)
    ax_c.plot([dom], [saved.max()], "o", ms=15, color=BLACK, mec="white", mew=2, zorder=4)
    ax_c.set(ylabel="Mean power", xlim=(axis[0], axis[-1]), ylim=(0, saved.max() * 1.08))
    hz_axis(ax_c)
    ax_c.yaxis.set_major_locator(MaxNLocator(3))
    ax_c.set_title("Time-averaged spectrum")
    ax_c.set_title("C", loc="left", fontweight="bold", fontsize=28)
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    ax_c.legend([Line2D([], [], ls="", marker="o", ms=15, color=BLACK, mec="white", mew=2),
                 Line2D([], [], color=GREY, ls="--", lw=2.5), Patch(fc=LIGHT, ec=GREY, lw=1.2)],
                ["Dominant frequency", "Centroid", "Centroid ± spread"], loc="center right",
                bbox_to_anchor=(1.0, 0.55), handlelength=1.4, borderaxespad=0.2, labelspacing=0.3)

    rows = [("Dominant frequency", f"{dom:.4f} Hz"), ("Spectral centroid", f"{cen:.4f} Hz"),
            ("Spectral spread", f"{spr:.4f} Hz"),
            ("Entropy (max 6.32 bits)", f"{feats.SpectralEntropy:.2f} bits"),
            ("Total band power", f"{feats.TotalBandPower:.3f}"),
            ("Peak power fraction", f"{feats.PeakPowerFraction:.3f}")]
    ax_t.set_axis_off()
    for k, (lab, val) in enumerate(rows):
        y = 1 - (k + 0.5) / len(rows)
        ax_t.text(0.0, y, lab, transform=ax_t.transAxes, ha="left", va="center", color=GREY)
        ax_t.text(1.0, y, val, transform=ax_t.transAxes, ha="right", va="center", color=BLACK)
    checks["Fig1"] = save(fig, out, "Fig1", "Fig1_single_cell_wavelet")
    facts["fig1"] = {"group": group, "embryo": emb, "cell": cell, "n_segments": len(segs),
                     "duration_min": float(feats.UsedDuration_min),
                     "frame_interval_s": float(feats.FrameInterval_s),
                     **{k: float(feats[k]) for k in common.FEATURES}}


# ------------------------------------------------------------------ Figure 2A
def fig2a(out, ex, facts, checks):
    from matplotlib.ticker import MaxNLocator
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    group, emb = ex["group"], ex["embryo"]
    lib = np.load(common.SPECTRAL_ROOT / common.GROUP_FOLDER[group] / emb / "SpectralLibrary.npz")
    red = np.load(RHYTHM_ROOT / group / f"{emb}_rednoise.npz")
    li = {c: k for k, c in enumerate(lib["CellID"].astype(str))}
    ri = {c: k for k, c in enumerate(red["CellID"].astype(str))}
    f = lib["Frequency"]
    fig = new_figure("Fig2A")
    axes = fig.subplots(1, 2)
    facts["fig2a"] = {"group": group, "embryo": emb}
    for ax, key, title in zip(axes, ("rhythmic", "nonrhythmic"), ("Exceeds red noise", "Within red noise")):
        c = ex[key]
        cid = str(c.CellID)
        top = max(red["Surrogate95"][ri[cid]].max(), lib["Power"][li[cid]].max())
        ax.fill_between(f, 0, red["Surrogate95"][ri[cid]], color=LIGHT, lw=0)
        ax.plot(f, red["SurrogateMean"][ri[cid]], color=GREY, ls="--", lw=2.5)
        ax.plot(f, np.maximum(lib["Power"][li[cid]], 0), color=BLACK, lw=3)
        ax.set(ylabel="Mean power", xlim=(f[0], f[-1]), ylim=(0, top * 1.05))
        hz_axis(ax)
        ax.yaxis.set_major_locator(MaxNLocator(3))
        ax.set_title(title)
        q = f"q = {c.Q_global:.3f}" if c.Q_global >= 0.001 else "q < 0.001"
        ax.text(0.98, 0.96, f"max {c.MaxRatio:.1f}× red noise\nat {c.MaxRatioFreqHz:.3f} Hz; {q}",
                transform=ax.transAxes, ha="right", va="top", color=BLACK, linespacing=1.15)
        facts["fig2a"][key] = {"cell": cid, "P_global": float(c.P_global), "Q_global": float(c.Q_global),
                               "MaxRatio": float(c.MaxRatio), "MaxRatioFreqHz": float(c.MaxRatioFreqHz),
                               "Alpha": float(c.Alpha)}
    axes[0].set_title("A", loc="left", fontweight="bold", fontsize=28)
    handles = [Line2D([], [], color=BLACK, lw=3), Line2D([], [], color=GREY, lw=2.5, ls="--"),
               Patch(color=LIGHT)]
    fig.legend(handles, ["Observed", "Red-noise mean", "Red-noise 95th percentile"],
               loc="outside lower center", ncol=3, handlelength=1.6, columnspacing=1.5)
    checks["Fig2A"] = save(fig, out, "Fig2A", "Fig2A_example_spectra")


# ------------------------------------------------------------------ dot plots (2B, 3B)
def dot_panel(ax, e: pd.DataFrame, value: str, log: bool = False, rng=None, legend: bool = False):
    rng = rng or np.random.default_rng(1)
    meds = {}
    for memb, st in GROUP_ORDER:
        d = e[(e.Membrane == memb) & (e.Stage == st)][value].dropna()
        x = XPOS[memb, st]
        ax.scatter(x + rng.uniform(-0.14, 0.14, len(d)), d, s=120, color=MEMB_COLOR[memb], alpha=0.9,
                   edgecolor="white", linewidth=1, zorder=3,
                   label=("Intact" if memb == "intact" else "Devitellinated") if (legend and st == 13) else None)
        med = float(d.median())
        ax.plot([x - 0.2, x + 0.2], [med, med], color=BLACK, lw=4.5, zorder=4, solid_capstyle="butt")
        meds[f"{memb}_{st}"] = {"n_embryos": int(len(d)), "median": med,
                                "q25": float(d.quantile(0.25)), "q75": float(d.quantile(0.75))}
    if log:
        ax.set_yscale("log")
    lo, hi = ax.get_ylim()
    if log:
        span = np.log10(hi) - np.log10(lo)
        ax.set_ylim(10 ** (np.log10(lo) - 0.22 * span), hi)
    else:
        ax.set_ylim(lo - 0.22 * (hi - lo), hi)
    for memb, st in GROUP_ORDER:
        n = meds[f"{memb}_{st}"]["n_embryos"]
        n_text(ax.text(XPOS[memb, st], 0.02, f"n={n}", transform=ax.get_xaxis_transform(), ha="center",
                       va="bottom", fontsize=N_PT, color=GREY))
    ax.set_xticks([(XPOS["intact", 13] + XPOS["devitellinated", 13]) / 2,
                   (XPOS["intact", 14] + XPOS["devitellinated", 14]) / 2], ["Stage 13", "Stage 14"])
    ax.set_xlim(-0.4, 2.55)
    ax.tick_params(axis="x", length=0)
    return meds


def fig2b(out, r, facts, checks):
    from scipy.stats import mannwhitneyu
    e = (r.groupby(["Group", "EmbryoID"])["Rhythmic"].mean() * 100).rename("Pct").reset_index()
    e = common.add_design(e)
    fig = new_figure("Fig2B")
    ax = fig.subplots()
    meds = dot_panel(ax, e, "Pct", legend=True)
    ax.set_ylabel("% of cells exceeding\nred noise")
    ax.set_title("Cells exceeding red noise, per embryo (bar = median)")
    ax.set_title("B", loc="left", fontweight="bold", fontsize=28)
    ax.legend(title="Vitelline membrane", loc="upper left", bbox_to_anchor=(1.01, 1.0),
              handletextpad=0.3, alignment="left")
    facts["fig2b"] = {k: {"n_embryos": v["n_embryos"], "median_pct": v["median"],
                          "n_with_any": int((e[(e.Membrane == k.split("_")[0]) &
                                               (e.Stage == int(k.split("_")[1]))].Pct > 0).sum())}
                      for k, v in meds.items()}
    for st in (13, 14):
        a = e[(e.Stage == st) & (e.Membrane == "intact")].Pct
        b = e[(e.Stage == st) & (e.Membrane == "devitellinated")].Pct
        facts["fig2b"][f"MWU_p_stage{st}"] = float(mannwhitneyu(a, b).pvalue)
    checks["Fig2B"] = save(fig, out, "Fig2B", "Fig2B_pct_rhythmic_per_embryo")


def fig3b(out, facts, checks):
    cells = common.load_mapped_cells()
    feats = ["TotalBandPower", "SpectralEntropy", "PeakPowerFraction"]
    e = cells.groupby(["Group", "EmbryoID"])[feats].median().reset_index()
    e = common.add_design(e)
    fig = new_figure("Fig3B")
    axes = fig.subplots(1, 3)
    titles = {"TotalBandPower": "Total band power", "SpectralEntropy": "Spectral entropy",
              "PeakPowerFraction": "Peak power fraction"}
    facts["fig3b"] = {"cells": "mapped cells (all-cell v3 mapping, repeated boundary IDs deduplicated)",
                      "n_cells": int(len(cells)), "n_embryos": int(e.EmbryoID.nunique())}
    for k, (ax, feat) in enumerate(zip(axes, feats)):
        meds = dot_panel(ax, e, feat, log=feat == "TotalBandPower", legend=k == 0)
        ax.set_title(titles[feat])
        if k == 0:
            ax.set_ylabel("Embryo median (log scale)")
        if feat == "SpectralEntropy":
            lo, hi = ax.get_ylim()
            ax.set_ylim(lo, 6.32 + 0.12 * (6.32 - lo))
            ax.axhline(6.32, color=GREY, ls="--", lw=2.5, zorder=1)
            ax.text(2.53, 6.32, "max 6.32", ha="right", va="bottom", color=GREY)
        facts["fig3b"][feat] = meds
        for st in (13, 14):
            facts["fig3b"][feat][f"devit_over_intact_median_stage{st}"] = \
                meds[f"devitellinated_{st}"]["median"] / meds[f"intact_{st}"]["median"]
    axes[0].set_title("B", loc="left", fontweight="bold", fontsize=28)
    fig.legend(*axes[0].get_legend_handles_labels(), title="Vitelline membrane", loc="outside upper center",
               ncol=2, handletextpad=0.3, columnspacing=1.5)
    checks["Fig3B"] = save(fig, out, "Fig3B", "Fig3B_embryo_medians")


# ------------------------------------------------------------------ maps (2C, 3A)
def grid_image(d: pd.DataFrame, col: str) -> np.ndarray:
    img = np.full((common.MAX_CNPU_Y + 1, common.MAX_CNPU_X + 1), np.nan)
    img[d.CNPU_Y.to_numpy(int), d.CNPU_X.to_numpy(int)] = d[col].to_numpy(float)
    return img


EXTENT = (common.ORIGIN_X, common.ORIGIN_X + (common.MAX_CNPU_X + 1) * common.CNPU_SIZE,
          common.ORIGIN_Y + (common.MAX_CNPU_Y + 1) * common.CNPU_SIZE, common.ORIGIN_Y)


def draw_map(ax, img, cmap, vmin, vmax, group, title, nline, footprint=False):
    im = ax.imshow(img, extent=EXTENT, cmap=cmap, vmin=vmin, vmax=vmax, interpolation="nearest")
    if footprint:  # outline of the area with data, so 0% (white) is not mistaken for no data
        from scipy import ndimage
        mask = ndimage.binary_fill_holes(ndimage.binary_closing(np.isfinite(img), iterations=2))
        yc = common.ORIGIN_Y + (np.arange(mask.shape[0]) + 0.5) * common.CNPU_SIZE
        xc = common.ORIGIN_X + (np.arange(mask.shape[1]) + 0.5) * common.CNPU_SIZE
        ax.contour(xc, yc, mask.astype(float), levels=[0.5], colors=GREY, linewidths=1.2)
    poly = common.load_template(common.SOX2[common.TEMPLATE_FOR_GROUP[group]])
    ax.plot(poly[:, 0], poly[:, 1], ":", color=BLACK, lw=2.2)
    ax.set(xlim=(150, 1100), ylim=(1190, 110))
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title(title)
    n_text(ax.set_xlabel(nline, fontsize=N_PT, color=GREY, linespacing=1.15))
    return im


def fig2c(out, min_emb, facts, checks):
    path = RHYTHM_ROOT / "rhythmic_fraction_maps" / "Group_Cnpu_RhythmicFraction.csv"
    t = pd.read_csv(path)
    t = t[t.NEmbryos >= min_emb]
    settings = json.loads((path.parent / "run_settings.json").read_text())
    vmax = settings["shared_pct_scale_max"]
    cmap = white_green().copy()
    cmap.set_bad((0, 0, 0, 0))
    fig = new_figure("Fig2C")
    axes = fig.subplots(1, 4)
    im = None
    for ax, g in zip(axes, common.GROUPS):
        n, c, k = settings["per_group_embryos_cells_rhythmic"][g]
        im = draw_map(ax, grid_image(t[t.Group == g], "PooledSmoothedPct"), cmap, 0, vmax, g,
                      PANEL_TITLE[g].replace(" · ", "\n"),
                      f"{n} embryos\n{c:,} mapped cells\n{k:,} exceed red noise", footprint=True)
    cb = fig.colorbar(im, ax=axes, extend="max", shrink=0.85, aspect=14, pad=0.015)
    cb.set_label("% of cells (smoothed)", color=GREY)
    cb.set_ticks([0, 10, 20, 30, 40])
    cb.outline.set_edgecolor(GREY)
    fig.suptitle("% of mapped cells exceeding red noise")
    fig.text(0.005, 0.995, "C", ha="left", va="top", fontsize=28, fontweight="bold")
    checks["Fig2C"] = save(fig, out, "Fig2C", "Fig2C_pct_rhythmic_maps")
    facts["fig2c"] = {"scale_max_pct": vmax, "min_embryos": min_emb,
                      "median_smoothed_pct": t.groupby("Group").PooledSmoothedPct.median().round(2).to_dict(),
                      "n_cnpus_shown": t.groupby("Group").size().to_dict()}


def fig3a(out, min_emb, facts, checks):
    import matplotlib.pyplot as plt
    parts = [pd.read_csv(common.MAPPING_ROOT / g / "Group_Cnpu_Spectral_Summary.csv") for g in common.GROUPS]
    t = pd.concat(parts)
    t = t[(t.Feature == "TotalBandPower") & (t.NEmbryos >= min_emb)]
    cells = common.load_mapped_cells()
    cells = cells[cells.InGrid]
    cmap = plt.get_cmap("jet").copy()
    cmap.set_bad((0, 0, 0, 0))
    fig = new_figure("Fig3A")
    axes = fig.subplots(2, 2)
    layout = [["S13", "S13_devitellinated"], ["S14", "S14_devitellinated"]]
    im = None
    for r_i, row in enumerate(layout):
        for c_i, g in enumerate(row):
            d = cells[cells.Group == g]
            im = draw_map(axes[r_i, c_i], grid_image(t[t.Group == g], "PooledSmoothedZ"), cmap, -1, 1, g,
                          PANEL_TITLE[g], f"{d.EmbryoID.nunique()} embryos · {len(d):,} cells")
    cb = fig.colorbar(im, ax=axes, extend="both", shrink=0.7, aspect=18, pad=0.02)
    cb.set_label("Smoothed z-score", color=GREY)
    cb.set_ticks([-1, -0.5, 0, 0.5, 1])
    cb.outline.set_edgecolor(GREY)
    fig.suptitle("Total band power (0.005–0.08 Hz), z-scored within each group")
    fig.text(0.005, 0.995, "A", ha="left", va="top", fontsize=28, fontweight="bold")
    checks["Fig3A"] = save(fig, out, "Fig3A", "Fig3A_band_power_zmaps")
    facts["fig3a"] = {"min_embryos": min_emb, "n_cnpus_shown": t.groupby("Group").size().to_dict(),
                      "z_within_group": True}


def main():
    global FONT
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--output-dir", type=Path,
                    default=common.ANALYSIS_ROOT / "step6_figures" / "sfn_poster")
    ap.add_argument("--example-group", default="S13_devitellinated")
    ap.add_argument("--min-embryos", type=int, default=2)
    args = ap.parse_args()
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    FONT = style()

    r = common.load_rhythmicity()
    master = common.load_master()[["Group", "EmbryoID", "CellID", "UsedDuration_min", "FrameInterval_s",
                                   *common.FEATURES]]
    r = r.merge(master[["Group", "EmbryoID", "CellID", "UsedDuration_min"]],
                on=["Group", "EmbryoID", "CellID"], validate="1:1")
    facts = {"rhythmic_definition": "Q_global < 0.05 (P_global, BH within embryo); "
                                    "poster wording: 'exceeds red noise'"}
    ex = pick_examples(r, args.example_group)
    facts["examples"] = {"group": ex["group"], "embryo": ex["embryo"],
                         "embryo_pct_rhythmic": ex["embryo_pct_rhythmic"],
                         "group_median_pct_rhythmic": ex["group_median_pct_rhythmic"],
                         "rule": " ".join(pick_examples.__doc__.split())}
    checks = {}
    fig1(out, ex, master, facts, checks)
    fig2a(out, ex, facts, checks)
    fig2b(out, r, facts, checks)
    fig2c(out, args.min_embryos, facts, checks)
    fig3a(out, args.min_embryos, facts, checks)
    fig3b(out, facts, checks)
    facts["style"] = {"font": FONT, "dpi": DPI, "transparent": True,
                      "placed_size_in": {k: list(v) for k, v in SIZES.items()},
                      "min_text_pt": MIN_PT, "min_n_line_pt": MIN_N_PT,
                      "colors": {"intact": INTACT, "devitellinated": DEVIT, "axes_secondary_text": GREY},
                      "file_checks": checks}
    (out / "figure_facts.json").write_text(json.dumps(facts, indent=2, default=str))
    print("facts written")


if __name__ == "__main__":
    main()
