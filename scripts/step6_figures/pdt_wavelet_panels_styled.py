#!/usr/bin/env python3
"""Poster-style (v2) versions of the PDT x wavelet panels 2, 3A and 3B.

Reads the tables written by pdt_wavelet_panels.py (no raw data) and redraws them
with the SfN poster style of rednoise_bandpower_figures.py: Liberation Sans / Arial,
text >= 24 pt (n lines >= 20 pt), intact #004E38, devitellinated #6AA84F,
grey #595959, exact placed size, 300 dpi, transparent PDF + PNG, and the same
QA pass (fonts, sizes, forbidden words, clipping, overlaps).

  PDT_Fig2   within-embryo r: PDT family (rows) x wavelet feature + Sox2 (columns),
             median over the family's setups. Bold = >= 90% of the family's setups
             are significant (BH q < 0.05, Wilcoxon on per-embryo Fisher z) AND share
             the sign of the family median. (The old rule - half the setups with
             q < 0.05 - marked every cell, because with 179 embryos almost every
             setup passes.)
  PDT_Fig3A  cells exceeding vs within red noise: per-setup median fold in event
             counts, plus per-embryo fold in event magnitude.
  PDT_Fig3B  devitellinated / intact fold of embryo medians: every setup's event
             count vs event magnitude and band power, per stage.

Wording: "exceeds / within red noise" (never rhythmic). PDT families are named
by input signal and threshold: "Magnitude" = change-in-magnitude PDTs (36-set),
"Intensity" = literature PDTs on intensity (54-set).

Usage:
    python pdt_wavelet_panels_styled.py                       # reads the default pdt_wavelet_panels folder
    python pdt_wavelet_panels_styled.py --in DIR --out DIR2
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import rednoise_bandpower_figures as rbf  # noqa: E402  (style, QA and save helpers)

DEFAULT_IN = Path("/sciclone/home/cabare01/scr10/CALCIUM/data/analysis_v3/step6_figures/pdt_wavelet_panels")
# Placed sizes (inches). Width matches the other poster figures; change here if the layout changes.
SIZES = {"PDT_Fig2": (14.6, 8.8), "PDT_Fig3A": (14.6, 8.4), "PDT_Fig3B": (14.6, 9.0)}
FAMILY_ORDER = ["36-CW", "36-EW", "36-P", "54-CW", "54-EW"]
FAMILY_LABEL = {"36-CW": ("Magnitude", "cell SD"), "36-EW": ("Magnitude", "embryo SD"),
                "36-P": ("Magnitude", "% baseline"), "54-CW": ("Intensity", "cell-wise"),
                "54-EW": ("Intensity", "embryo-wise")}
FEATURES = ["Dominant frequency", "Spectral centroid", "Spectral spread", "Spectral entropy",
            "Band power", "Peak power fraction", "Sox2 level"]
AGREE = 0.90
SETUP_DOT = "#8A8A8A"


def fam_label(f, n=None):
    a, b = FAMILY_LABEL.get(f, (f, ""))
    return f"{a}\n{b}" + (f"\n{n} setups" if n is not None else "")


def fold_ticks(ax, lo, hi, axis="x", step=1):
    """Log2 fold axis: labelled ticks every `step` doublings (anchored at 1x), minor ticks at each doubling."""
    ticks = [t for t in range(-6, 6) if lo <= t <= hi and t % step == 0]
    minor = [t for t in range(-6, 6) if lo <= t <= hi and t % step != 0]
    (ax.set_xticks if axis == "x" else ax.set_yticks)(minor, minor=True)
    labels = [f"{2.0 ** t:g}×" for t in ticks]
    if axis == "x":
        ax.set_xticks(ticks, labels)
        ax.set_xlim(lo, hi)
        ax.axvline(0, color=rbf.GREY, lw=1.5, zorder=1)
    else:
        ax.set_yticks(ticks, labels)
        ax.set_ylim(lo, hi)
        ax.axhline(0, color=rbf.GREY, lw=1.5, zorder=1)


def span(values, pad=0.3):
    v = np.asarray(values, float)
    v = v[np.isfinite(v)]
    return np.floor(v.min() * 2) / 2 - pad, np.ceil(v.max() * 2) / 2 + pad


def strip(ax, y, x, rng, color=SETUP_DOT, size=90, jitter=0.2, median=True, edge="none"):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    ax.scatter(x, y + rng.uniform(-jitter, jitter, len(x)), s=size, c=color, alpha=0.8 if edge == "none" else 1,
               edgecolor=edge, linewidth=1.2, zorder=3)
    if median and len(x):
        ax.plot([np.median(x)] * 2, [y - 0.36, y + 0.36], color=rbf.BLACK, lw=4, zorder=4)


# ----------------------------------------------------------------------------- Fig 2
def fig2(src, out, facts):
    from matplotlib.colors import TwoSlopeNorm
    long = pd.read_csv(src / "panel2_setup_x_feature_r.csv")
    fams = [f for f in FAMILY_ORDER if f in set(long.family)]
    feats = [f for f in FEATURES if f in set(long.feature)]
    med = long.groupby(["family", "feature"]).mean_r.median().unstack().loc[fams, feats]

    def agree(g):
        s = np.sign(g.mean_r.median())
        return float(((np.sign(g.mean_r) == s) & (g.q_bh < 0.05)).mean())

    agr = long.groupby(["family", "feature"]).apply(agree).unstack().loc[fams, feats]
    nset = long.groupby("family").setup.nunique()
    n_emb = int(long.n_embryos.max())
    lim = 0.8
    if np.nanmax(np.abs(med.to_numpy())) > lim:
        raise RuntimeError("Fig 2: |median r| exceeds the colour range")

    fig = rbf.new_figure("PDT_Fig2")
    ax = fig.add_subplot()
    # rows = features (horizontal labels), columns = PDT families
    im = ax.imshow(med.T.to_numpy(), cmap="RdBu_r", norm=TwoSlopeNorm(0, -lim, lim), aspect="auto")
    for j, f in enumerate(fams):
        for i, y in enumerate(feats):
            v = med.loc[f, y]
            ax.text(j, i, f"{v:.2f}".replace("-", "−"), ha="center", va="center",
                    color="white" if abs(v) > 0.45 else rbf.BLACK,
                    fontweight="bold" if agr.loc[f, y] >= AGREE else "normal")
    ax.set_yticks(range(len(feats)), feats)
    ax.set_xticks(range(len(fams)), [fam_label(f, int(nset[f])) for f in fams])
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    if "Sox2 level" in feats:
        ax.axhline(feats.index("Sox2 level") - 0.5, color=rbf.BLACK, lw=2.5)
    cb = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.02, ticks=[-0.8, -0.4, 0, 0.4, 0.8])
    cb.ax.set_yticklabels(["−0.8", "−0.4", "0", "0.4", "0.8"])
    cb.set_label("Within-embryo r")
    cb.outline.set_visible(False)
    fig.suptitle("Peak-detection event counts vs wavelet features", x=0.01, ha="left")
    rbf.n_text(ax.set_title(f"n = {n_emb} embryos · {int(nset.sum())} peak-detection setups; "
                            "value = median over the family's setups", loc="left", fontsize=rbf.N_PT))
    ax.set_xlabel(f"Bold: ≥ {AGREE:.0%} of the family's setups agree in sign (BH q < 0.05)", labelpad=12)
    facts["PDT_Fig2"] = {
        "file": "PDT_Fig2_family_summary", "n_embryos": n_emb, "n_setups": int(nset.sum()),
        "setups_per_family": {f: int(nset[f]) for f in fams},
        "median_r": {f: {y: round(float(med.loc[f, y]), 3) for y in feats} for f in fams},
        "p10_r": {f: {y: round(float(long[(long.family == f) & (long.feature == y)].mean_r.quantile(.1)), 3)
                      for y in feats} for f in fams},
        "p90_r": {f: {y: round(float(long[(long.family == f) & (long.feature == y)].mean_r.quantile(.9)), 3)
                      for y in feats} for f in fams},
        "frac_setups_agree": {f: {y: round(float(agr.loc[f, y]), 3) for y in feats} for f in fams},
        "bold_rule": f">= {AGREE:.0%} of the family's setups have BH q < 0.05 (Wilcoxon signed-rank on "
                     "per-embryo Fisher z, BH over all setup x feature pairs) and the sign of the family median",
        "method": "Pearson r between cells within each embryo (embryos with >= 30 joined cells), Fisher z "
                  "averaged across embryos and back-transformed. Event counts (Total_<code>), band power and "
                  "peak power fraction log10-transformed (offset = half the 1st percentile of positive values); "
                  "Sox2 = CY5_Intensity (Sox2 FISH), untransformed.",
    }
    return rbf.save(fig, out, "PDT_Fig2", "PDT_Fig2_family_summary")


# ----------------------------------------------------------------------------- Fig 3A
def fig3a(src, out, facts):
    from matplotlib.lines import Line2D
    per = pd.read_csv(src / "panel3A_per_setup.csv")
    emb = pd.read_csv(src / "panel3A_per_embryo_log2_fold.csv")
    fams = [f for f in FAMILY_ORDER if f in set(per.family)]
    mag = emb.dropna(subset=["EventMag"]) if "EventMag" in emb else emb.iloc[:0]
    lo, hi = span(np.concatenate([per.median_log2.to_numpy(), mag.EventMag.to_numpy() if len(mag) else []]))
    rng = np.random.default_rng(0)

    fig = rbf.new_figure("PDT_Fig3A")
    gs = fig.add_gridspec(2, 1, height_ratios=[len(fams), 1.25 if len(mag) else 0.01])
    ax = fig.add_subplot(gs[0])
    for k, f in enumerate(fams):
        strip(ax, -k, per.loc[per.family == f, "median_log2"], rng)
    ax.set_yticks([-k for k in range(len(fams))], [fam_label(f) for f in fams])
    ax.set_ylim(-len(fams) + 0.45, 0.55)
    fold_ticks(ax, lo, hi)
    ax.tick_params(axis="x", labelbottom=False)
    ax.set_title("Event counts: one dot per peak-detection setup (median over embryos)", loc="left")
    n_i, n_d = int((emb.membrane == "intact").sum()), int((emb.membrane == "devitellinated").sum())
    fig.suptitle("Cells exceeding red noise vs cells within red noise", x=0.01, ha="left")
    facts["PDT_Fig3A"] = {
        "file": "PDT_Fig3A_red_noise_vs_pdt", "n_embryos": len(emb), "n_intact": n_i, "n_devitellinated": n_d,
        "rule": "embryos with >= 5 cells in each class; per embryo log2(mean count exceeding / mean count within) "
                "per setup; dot = median over embryos; exceeds red noise = Q_global < 0.05",
        "median_fold_all_setups": round(float(2 ** per.median_log2.median()), 3),
        "frac_setups_fold_gt_1": round(float((per.median_log2 > 0).mean()), 3),
        "frac_setups_q_lt_0.05": round(float((per.q_bh < 0.05).mean()), 3),
        "median_fold_by_family": {f: round(float(2 ** per.loc[per.family == f, "median_log2"].median()), 3)
                                  for f in fams},
    }
    if len(mag):
        from scipy import stats
        ax2 = fig.add_subplot(gs[1], sharex=ax)
        cols = [rbf.DEVIT if m == "devitellinated" else rbf.INTACT for m in mag.membrane]
        strip(ax2, 0, mag.EventMag, rng, color=cols, size=110, jitter=0.25, edge="white")
        ax2.set_yticks([0], ["Event\nmagnitude"])
        ax2.set_ylim(-0.55, 0.55)
        fold_ticks(ax2, lo, hi)
        m = float(2 ** mag.EventMag.median())
        p = float(stats.wilcoxon(mag.EventMag).pvalue)
        ax2.set_title(f"Event size: one dot per embryo · median {m:.2f}×", loc="left")
        ax2.set_xlabel("Fold difference, exceeding ÷ within red noise (log scale)")
        facts["PDT_Fig3A"].update({"event_magnitude_median_fold": round(m, 3), "event_magnitude_n": len(mag),
                                   "event_magnitude_wilcoxon_p": p,
                                   "event_magnitude_rule": "per embryo log2(median of per-cell mean "
                                                           "Peak_Mag_absDiff, exceeding / within)"})
    else:
        ax.tick_params(axis="x", labelbottom=True)
        ax.set_xlabel("Fold difference, exceeding ÷ within red noise (log scale)")
    handles = [Line2D([], [], marker="o", ls="", ms=11, mfc=SETUP_DOT, mec="none", label="Peak-detection setup"),
               Line2D([], [], color=rbf.BLACK, lw=4, label="Median"),
               Line2D([], [], marker="o", ls="", ms=12, mfc=rbf.INTACT, mec="white", label="Intact embryo"),
               Line2D([], [], marker="o", ls="", ms=12, mfc=rbf.DEVIT, mec="white", label="Devitellinated embryo")]
    leg = fig.legend(handles=handles, loc="outside lower center", ncol=2, columnspacing=2.5, handletextpad=0.5)
    leg.set_title(f"n = {len(emb)} embryos with ≥ 5 cells in each class ({n_i} intact, {n_d} devitellinated)",
                  prop={"size": rbf.N_PT})
    rbf.n_text(leg.get_title())
    return rbf.save(fig, out, "PDT_Fig3A", "PDT_Fig3A_red_noise_vs_pdt")


# ----------------------------------------------------------------------------- Fig 3B
def fig3b(src, out, facts):
    from matplotlib.lines import Line2D
    st = pd.read_csv(src / "panel3B_stats.csv")
    fams = [f for f in FAMILY_ORDER if f in set(st.family)]
    refs = [(c, l) for c, l in (("EventMag", "Event\nmagnitude"), ("TotalBandPower", "Band\npower"))
            if c in set(st.measure)]
    rows = [("fam", f) for f in fams] + [("ref", c) for c, _ in refs]
    ypos = {r: -(k + (0.6 if r[0] == "ref" else 0)) for k, r in enumerate(rows)}
    lo, hi = span(st.log2_fold)
    hi += 0.9                                               # room for the diamond labels
    rng = np.random.default_rng(1)

    fig = rbf.new_figure("PDT_Fig3B")
    axes = fig.subplots(1, 2, sharey=True)
    facts["PDT_Fig3B"] = {"file": "PDT_Fig3B_membrane_all_setups",
                          "rule": "fold = median over devitellinated embryos / median over intact embryos of the "
                                  "embryo value (mean count per cell; median band power / event magnitude); "
                                  "Mann-Whitney per measure, BH within stage", "stages": {}}
    for ax, stage in zip(axes, (13, 14)):
        s = st[st.stage == stage]
        for f in fams:
            strip(ax, ypos[("fam", f)], s.loc[s.family == f, "log2_fold"], rng)
        for c, _ in refs:
            v = s[s.measure == c]
            if len(v):
                x = float(v.log2_fold.iloc[0])
                ax.scatter([x], [ypos[("ref", c)]], s=260, marker="D", c=rbf.BLACK, zorder=4)
                ax.text(x + 0.22, ypos[("ref", c)], f"{2 ** x:.1f}×", va="center", ha="left")
        ax.axhline((ypos[("fam", fams[-1])] + ypos[("ref", refs[0][0])]) / 2 if refs else -99,
                   color=rbf.LIGHT, lw=2)
        fold_ticks(ax, lo, hi, step=2)
        cnt = s[s.measure.isin(st.measure[st.family.isin(fams)])]
        ni, nd = int(cnt.n_intact.max()), int(cnt.n_devit.max())
        ax.set_title(f"Stage {stage}", loc="left")
        rbf.n_text(ax.set_title(f"n = {ni} / {nd} embryos", loc="right", fontsize=rbf.N_PT))
        facts["PDT_Fig3B"]["stages"][str(stage)] = {
            "n_intact": ni, "n_devitellinated": nd,
            "median_fold_setups": round(float(2 ** cnt.log2_fold.median()), 3),
            "frac_setups_higher_devit": round(float((cnt.log2_fold > 0).mean()), 3),
            "frac_setups_sig_higher": round(float(((cnt.log2_fold > 0) & (cnt.q_bh < 0.05)).mean()), 3),
            "frac_setups_sig_lower": round(float(((cnt.log2_fold < 0) & (cnt.q_bh < 0.05)).mean()), 3),
            "median_fold_by_family": {f: round(float(2 ** s.loc[s.family == f, "log2_fold"].median()), 3)
                                      for f in fams},
            **{c: {"fold": round(float(2 ** v.log2_fold.iloc[0]), 3), "p": float(v.p.iloc[0]),
                   "n_intact": int(v.n_intact.iloc[0]), "n_devit": int(v.n_devit.iloc[0])}
               for c, _ in refs for v in [s[s.measure == c]] if len(v)},
        }
    axes[0].set_yticks([ypos[r] for r in rows],
                       [fam_label(r[1]) if r[0] == "fam" else dict(refs)[r[1]] for r in rows])
    axes[0].set_ylim(min(ypos.values()) - 0.55, 0.55)
    axes[1].tick_params(axis="y", length=0)
    fig.suptitle("Membrane removal: event counts vs event size and wavelet power", x=0.01, ha="left")
    for ax in axes:
        ax.set_xlabel("Devitellinated ÷ intact")
    handles = [Line2D([], [], marker="o", ls="", ms=11, mfc=SETUP_DOT, mec="none",
                      label="One peak-detection setup"),
               Line2D([], [], color=rbf.BLACK, lw=4, label="Median of setups"),
               Line2D([], [], marker="D", ls="", ms=13, mfc=rbf.BLACK, mec="none", label="Single measure")]
    leg = fig.legend(handles=handles, loc="outside lower center", ncol=3, columnspacing=2.0, handletextpad=0.5)
    leg.set_title("n = intact / devitellinated embryos", prop={"size": rbf.N_PT})
    rbf.n_text(leg.get_title())
    return rbf.save(fig, out, "PDT_Fig3B", "PDT_Fig3B_membrane_all_setups")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in", dest="src", type=Path, default=DEFAULT_IN)
    ap.add_argument("--out", type=Path, default=None, help="default: <in>/poster_v2")
    args = ap.parse_args()
    out = args.out or args.src / "poster_v2"
    out.mkdir(parents=True, exist_ok=True)
    rbf.FONT = rbf.style()
    rbf.SIZES.update(SIZES)
    facts, checks = {"source_tables": str(args.src),
                     "wording": "exceeds / within red noise = step 1 Rhythmic flag (Q_global < 0.05)"}, {}
    checks["PDT_Fig2"] = fig2(args.src, out, facts)
    checks["PDT_Fig3A"] = fig3a(args.src, out, facts)
    checks["PDT_Fig3B"] = fig3b(args.src, out, facts)
    facts["style"] = {"font": rbf.FONT, "dpi": rbf.DPI, "transparent": True,
                      "placed_size_in": {k: list(v) for k, v in SIZES.items()},
                      "min_text_pt": rbf.MIN_PT, "min_n_line_pt": rbf.MIN_N_PT,
                      "colors": {"intact": rbf.INTACT, "devitellinated": rbf.DEVIT, "setup": SETUP_DOT,
                                 "axes_secondary_text": rbf.GREY}, "file_checks": checks}
    (out / "pdt_figure_facts.json").write_text(json.dumps(facts, indent=2, default=str))
    print(f"facts written -> {out / 'pdt_figure_facts.json'}")


if __name__ == "__main__":
    main()
