#!/usr/bin/env python3
"""Step 5b: join PDT event metrics to spectral features per cell, and test agreement.

Join: event file `ID` == spectral `CellID` within the same embryo (both come
from the same tracking). The match rate is reported per embryo; a low rate
(< --min-match %) means the event file was built from a different
recording/tracking (as found for two boundary files) and the embryo is
excluded from the tests.

CONFIRM COLUMN MEANINGS with the PDT paper before interpreting: the event
files hold ID, Cy5Intensity, Cy5_Norm1/2, Counts_0, Counts_M, Th_SD_1..5,
MeanSpikeMag, VarSpikeMag, MeanLocalPeakMag, VarLocalPeakMag,
MeanSpikeDuration, VarSpikeDuration. Pass the ones to test with
--event-columns (default: every numeric column except ID).

Tests (embryo = unit):
  1. Per embryo, Spearman rho between each event column and each spectral
     feature over matched cells (>= --min-cells). Per group: median rho,
     % embryos with rho > 0, Wilcoxon of rho vs 0. BH-FDR over all
     group x event x feature tests.
     Expected if the two methods agree: event rate/magnitude correlates
     positively with TotalBandPower and PeakPowerFraction.
  2. If step 1 exists: per embryo, median event value in rhythmic minus
     non-rhythmic cells; Wilcoxon across embryos per group.

Outputs in ANALYSIS_ROOT/step5_event_integration/:
  events_spectral_cells.csv   matched cells (event + spectral + template columns)
  event_match_by_embryo.csv
  event_spectral_rho_by_embryo.csv, event_spectral_rho_tests.csv
  event_rhythmic_contrast.csv (when step 1 output exists)
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


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--event-columns", nargs="*", default=None)
    ap.add_argument("--min-match", type=float, default=50.0)
    ap.add_argument("--min-cells", type=int, default=50)
    args = ap.parse_args()
    out = common.out_dir("step5_event_integration")
    index = pd.read_csv(out / "event_file_index.csv", keep_default_na=False)
    index = index[index.EventFile != ""]

    master = common.load_master()
    master["EmbryoKey"] = master["EmbryoID"].map(common.embryo_ref.embryo_key)
    mapped = common.load_mapped_cells()[["Group", "EmbryoID", "CellID", "Centroid_1", "Centroid_2",
                                         "CNPU_X", "CNPU_Y", "InGrid"]]

    joined, match = [], []
    for r in index.itertuples():
        ev = pd.read_csv(r.EventFile)
        ev["CellID"] = ev["ID"].astype(str).str.strip()
        ev = ev.drop(columns="ID").drop_duplicates("CellID")
        sp = master[master.EmbryoKey == r.EmbryoKey]
        both = sp.merge(ev, on="CellID", how="inner")
        pct = 100 * len(both) / len(sp) if len(sp) else np.nan
        match.append({"Group": r.Group, "EmbryoKey": r.EmbryoKey,
                      "EmbryoID": sp.EmbryoID.iloc[0] if len(sp) else "",
                      "NEventCells": len(ev), "NSpectralCells": len(sp), "NMatched": len(both),
                      "PctSpectralMatched": pct, "Usable": bool(pct >= args.min_match)})
        if pct >= args.min_match:
            joined.append(both)
    match = pd.DataFrame(match)
    match.to_csv(out / "event_match_by_embryo.csv", index=False)
    print(match.groupby("Group")[["Usable"]].agg(["sum", "size"]))
    cells = pd.concat(joined, ignore_index=True).merge(
        mapped, on=["Group", "EmbryoID", "CellID"], how="left")
    ev_cols = args.event_columns or [c for c in ev.columns if c != "CellID"
                                     and pd.api.types.is_numeric_dtype(cells[c])]
    cells.to_csv(out / "events_spectral_cells.csv", index=False)

    rows = []
    for (g, e), d in cells.groupby(["Group", "EmbryoID"]):
        if len(d) < args.min_cells:
            continue
        for ec in ev_cols:
            for f in common.FEATURES:
                v = d[[ec, f]].dropna()
                rho = stats.spearmanr(v[ec], v[f]).statistic if v[ec].nunique() > 1 else np.nan
                rows.append({"Group": g, "EmbryoID": e, "EventColumn": ec, "Feature": f,
                             "NCells": len(v), "Rho": rho})
    rho = pd.DataFrame(rows)
    rho.to_csv(out / "event_spectral_rho_by_embryo.csv", index=False)
    tests = []
    for (g, ec, f), d in pd.concat([rho, rho.assign(Group="All")]).groupby(["Group", "EventColumn", "Feature"]):
        r = d["Rho"].dropna()
        w = stats.wilcoxon(r) if len(r) >= 6 else None
        tests.append({"Group": g, "EventColumn": ec, "Feature": f, "NEmbryos": len(r),
                      "MedianRho": r.median(), "PctPositive": 100 * (r > 0).mean() if len(r) else np.nan,
                      "Wilcoxon_p": w.pvalue if w else np.nan})
    t = pd.DataFrame(tests)
    t["q_bh"] = common.bh(t["Wilcoxon_p"])
    t.to_csv(out / "event_spectral_rho_tests.csv", index=False)
    print(t[t.Group == "All"].sort_values("Wilcoxon_p").head(20).round(4).to_string(index=False))

    try:
        rh = common.load_rhythmicity()[["Group", "EmbryoID", "CellID", "Rhythmic"]]
    except FileNotFoundError:
        return
    d = cells.merge(rh, on=["Group", "EmbryoID", "CellID"], how="inner")
    contrast = []
    for ec in ev_cols:
        per = d.groupby(["Group", "EmbryoID", "Rhythmic"])[ec].median().unstack("Rhythmic")
        per = per.reindex(columns=[True, False]).dropna()
        groups = list(per.index.get_level_values("Group").unique())
        for g in [*groups, "All"]:
            p = per if g == "All" else per.xs(g, level="Group")
            diff = p[True] - p[False]
            w = stats.wilcoxon(diff) if len(diff) >= 6 and diff.abs().sum() > 0 else None
            contrast.append({"Group": g, "EventColumn": ec, "NEmbryos": len(diff),
                             "MedianRhythmicMinusNon": diff.median(),
                             "Wilcoxon_p": w.pvalue if w else np.nan})
    c = pd.DataFrame(contrast)
    c["q_bh"] = common.bh(c["Wilcoxon_p"])
    c.to_csv(out / "event_rhythmic_contrast.csv", index=False)


if __name__ == "__main__":
    main()
