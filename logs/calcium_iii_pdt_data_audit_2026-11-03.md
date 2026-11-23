# Audit of `data10/PDT_data` (Sudip Paudel's 2022 PDT outputs) for step 5 — 2026-11-03

Read-only audit. Nothing in `PDT_data`, `boundary_info` or the v3 outputs was modified. Audit scripts and intermediate CSVs are in the session scratchpad (`audit_percell.py`, `out/percell_*.csv`, `out/quick_rho_by_embryo.csv`); nothing was written to `scr10/.../analysis_v3`.

## 1. Coverage and joins: Established

**Per-cell tables** (`*_Unique_LastFrameIDs_Perhour_Counts.csv`, 20 files = 4 groups × {36-CW, 36-EW, 36-P, 54-CW, 54-EW}):
- All 179 reference-table embryos are present in the right group (46/46/46/41). No embryo is in the wrong group file.
- Row counts are 53,016 / 48,247 / 55,563 / 44,583 (S13 / S13dev / S14 / S14dev) = **201,409**, exactly the manuscript's CNP n.
- `ID == Calcium_ID` in every row. There are no duplicate `(Name, Calcium_ID)` pairs.

**Name variants:** `embryo_ref.embryo_key` alone does not resolve these PDT names:

| PDT `Name` | Reference key | Problem |
|---|---|---|
| `180509P_E7_`, `180509P_E6_` | `180509SP_E7`, `180509SP_E6` | Initials truncated |
| `191030SPRN_` | `191030SPRN_E3` | **No embryo number** (unique by date + initials) |
| `191020RNSP_` | `191020RNSP_E6` | **No embryo number** (unique by date + initials) |
| `200313SP-E1`, `200313SP-E2` | `200313SP_E1`, `200313SP_E2` | Hyphen instead of underscore |

Step 5 code needs this explicit alias table. A date + E-number fallback is unsafe: `190918SP_E1` (a boundary file not in the reference table) collides with `190918RN_E1`.

**PDT per-cell vs boundary files:**
- The PDT per-cell rows are the boundary-file cells (median 97.6% of boundary cells), with identical Cy5 and template centroids.
- Exception: **190213SP_E2** has the same IDs and Cy5, but its centroids are offset from the current boundary file (median 5 px, max 46 px). The PDT outputs used a different registration of this embryo. v3 matches the current boundary file exactly.

**PDT per-cell vs v3 mapped cells:**
- 100% of the 191,527 v3 mapped cells have a PDT row.
- Centroids are exactly equal in 175 of 177 embryos. The two exceptions are 190213SP_E2 (above) and 200217SP_E3 (98.7%; the known repeated `Calcium_ID`s).
- **Step 5 can use all 177 mapped embryos, not 79.**

**Excluded embryos:** the PDT per-cell tables for 180405SP_E4 and 200114SP_E1 contain only 396 and 518 cells (36% and 38% of their boundary cells). So the paper's CNP was built from the same mismatched boundary files. This supports keeping the exclusion.

## 2. Event-level files (`*_ActivityWiseDetails_*`): Established

**Inventory:**
- 36-set: CW, EW and P for all 4 groups.
- 54-set: EW for all 4 groups; **CW for the devitellinated groups only.** There are no 54-CW event files for intact S13 or S14.

**Structure:**
- **One row per event.** `ActivityID` is unique within (embryo, condition), and `MarkStart` = 1 in every row.
- Each row names **one** cell (`ID`) and gives `No_of_PPcells`; rosette and wave participants are not listed. For example, in S13dev 36-EW, wave events have median 8 and max 1,611 participating cells per row.
- So **per-cell event times exist only for spikes** (single-cell events) and for whichever cell is listed as a multi-cell event's `ID`.

**ISI / CV_ISI feasibility** (S13dev; share of last-frame cells with ≥3 spike rows):
- 36-EW: median PDT 6%, range 0.4–56%.
- 54-EW: median PDT 6%, range 0.3–96%.

CV_ISI from PDT files works only for the most permissive PDTs, and only for spikes.

## 3. Condition codes: Established from dissertation Flowchart 4.1 (corrected same day)

> **Correction.** The first version of this section inferred the codes from the data and got **digit 1 wrong** (it said 1 = none, 2 = min-subtract). It also misread the 3rd digit of the 36-set as the threshold approach and left `P` unexplained. Sudip Paudel's dissertation (`data10/Paudel_wmgrad_0261E_10696.pdf`, Flowchart 4.1, printed p. 92) settles all three.

**Five-digit code:**
1. Preprocessing: **1 background subtraction** (frame-wise minimum) · **2 none** · 3 Savitzky–Golay
2. Input signal: 1 Local · 2 Modified; "–" for literature PDTs
3. Baseline: 1 lower 10% · 2 mean* · 3 moving mean
4. Threshold technique: 1 SD* · 2 % of baseline · 3 findpeaks*
5. Parameter: 1/2 = 1SD/2SD, 110%/200%, or mean+1SD/mean+2SD

\* = also computed embryo-wise.

**Literature PDTs (54-set):**
- File code `Split_W_<d1 d3 d4 d5>`; paper name `d1-d3d4d5.CW|EW`.
- The 12 codes with baseline 1 or 3 and %-of-baseline thresholding have no EW form, giving 54 CW + 42 EW.
- **Those 12 `_EW` columns still hold values in the per-cell 54-EW files and must not be used.**

**Change-in-magnitude PDTs (36-set):**
- File code `Split_W_<d1 d2 d3 d5>` = preprocessing, Local/Modified, **baseline**, Th1/Th2.
- **The threshold technique is carried by `ThType`:** CW = SD cell-wise, EW = SD embryo-wise, **P = percentage of baseline** (cell-wise only, hence labelled `.CW` in the paper's names). findpeaks is not used.
- Local gives about 6× more events than Modified. The event files start one frame later (ND.T 2 vs 1).

**`fpf`:** not mentioned in the dissertation. The Protocol 4 run script checks for 18 findpeaks columns and reruns them with `Get_FindPeaks_Revised`, so "fpf" plausibly means "findpeaks fixed" `[CONFIRM]`.

**Per-cell columns** (dissertation `Get_CellWise_Activity_Counts`):
- Each count column = raw count × `TotalFrames`/`CellLengths`, i.e. scaled to the full recording. Verified: value × `CellLengths`/`TotalFrames` is an integer for every value.
- `Total` = the number of collective-activity IDs the cell joined, including rosettes and waves. `Only*` = the same, split by event type.
- `ActivePeaksPerHr` = the number of above-threshold runs.
- `CellWiseActiveIncidence` = % of frames active.
- These are the paper's NNDist `Category` columns.

## 4. Paper CNPU maps (`AllGroups_Data_for_Plotting.csv`): Established

- Groups present: 4 stages × {Spikes, Rosettes, Waves}, with 13–37 PDT groups each. Which PDTs belong to each `GroupID` is **not** in the file, only `GroupSize`.
- **Grid phase differs from v3 and differs per group.** CenX mod 9 = 3.80 / 0.72 / 8.64 / 6.14, versus v3's fixed origin (120, 1). The nearest v3 CNPU centre is a median 3.6–4.5 px away. CNPU-by-CNPU comparison needs re-gridding or a nearest-centre match; `bID` ≠ `CNPU_ID`.
- **The paper keeps only CNPUs with ≥3 cells:**
  - That is 4,909 / 4,738 / 5,501 / 4,651 CNPUs, against 7.4k–9.6k occupied.
  - These CNPUs hold 70–93% of cells.
  - A further 283–533 CNPUs with ≥3 cells are also dropped.
  - v3 combined maps use ≥2 embryos instead.
  - **Update (same day, from Supp. Text S6):** the text gives the rule as CNPUs with **≥3 embryos** inside a **rubber-band boundary of n = 5 embryos**.
  - **Resolved (same day, from Sudip's R code, `Activity_Per_hour_Counts_Parallel_220305SP.R`):**
    1. Take the convex hull of template pixels covered by ≥5 embryos.
    2. Keep cells inside the hull.
    3. Build the grid with its origin at those cells' minimum centroid (verified: each group's grid origin equals a real cell centroid).
    4. Keep CNPUs with **≥3 cells**.
    - So the "≥3 cells" reading from the data was right, and the text's "≥3 embryos" is not what the code did.

## 6. Addendum (same day): Sudip's R scripts (`PDT_data/R_ScriptsSP/`)

- **`fpf` = "findpeaks fixed":** the 18 findpeaks PDTs were rerun on the 54-set files (`Findpeaks_iPeaks_Only_{CW,EW}.R`).
- **P = percentage-of-baseline** (`1_Functions_36_cond_P.R`).
- **The cell-minimum subtraction is EW-only**, not applied after every preprocessing step.
- **Event files:** the listed `ID` is a leader cell at the start frame, in image coordinates; the files include all tracked cells.
- **Grouping:** complete-linkage `hclust` on 1 − r, cut at 0.5.
- **Smoothing:** 49 nearest CNPUs within 39 px, padded by repetition at edges.
- **Sudip's spike CV** = sd(gap)/(mean(gap) − min(gap)), in frames. His output is not on the cluster.
- **Not on the cluster** (need W: drive):
  - `MethodGroups_RealData.csv`
  - `*_Spatial_Moving_Average.csv`
  - `*MatrixForm_Only_IDs.csv`
  - the CV outputs
  - the frame-level `*_Conditions_Wave_*` files
  - `All_NeuralPlate_Regional_Boundaries.csv`
  - Supp. Table S3's 27 retained groups match `AllGroups` `GroupID` numbering; all have GroupSize ≥ 11.

## 5. Quick descriptive check: PDT rate vs spectral features: Preliminary, not a step 5 result

Method: per embryo, the Spearman ρ over mapped cells between `ActivePeaksPerHr_<code>` and each spectral feature. The table gives the median of embryo-level ρ (177 embryos), with the median over PDTs in each set.

| Feature | 36-CW | 36-P | 54-CW | 54-EW |
|---|---|---|---|---|
| TotalBandPower | −0.22 | −0.10 | −0.08 | **+0.39** |
| PeakPowerFraction | −0.13 | −0.08 | −0.09 | +0.15 |
| SpectralCentroid | **+0.38** | +0.20 | +0.28 | **−0.31** |
| DominantFreqHz | +0.29 | +0.16 | +0.16 | −0.23 |

Individual PDTs span about ±0.8 for every feature. Examples:

| PDT | Feature | S13dev ρ | S14dev ρ | Intact ρ |
|---|---|---|---|---|
| 54-EW 1112 = `1-112.EW` (background-subtracted, lower-10% baseline, SD, 2SD) | TotalBandPower | +0.70 | +0.63 | +0.37 / +0.41 |
| 54-EW 2112 = `2-112.EW` (no preprocessing; otherwise the same) | TotalBandPower | +0.71 | +0.66 | +0.39 / +0.45 |
| 54-CW 1111 | TotalBandPower | −0.59 | −0.52 | about 0 |

**The sign of PDT–spectral agreement depends on the PDT**, mainly on CW vs EW:
- **EW PDTs** (one absolute threshold per embryo) count large transients, which carry band power.
- **CW PDTs** normalise each cell by its own SD. They count small crossings in low-variance cells, which have flat, high-centroid spectra.

This matches the paper's CW/EW-dependent wave bias, and it supports the caveat that frequency and centroid maps partly reflect SNR. The step 5 docstring's expectation ("event rate correlates positively with TotalBandPower") holds only for EW PDTs.

## Consequences for the plan

- The step 5 inputs exist for all 177 mapped embryos. `find_event_files.py` (which globs ConvergeSpikeCy5) and the "79/179" fallback are obsolete.
- `event_spectral_join.py` must read the PDT per-cell tables and use the alias table above.
- Pre-specify the PDT panel; do not pick a single PDT post hoc.
- Report agreement across all 204 PDTs (the distribution of embryo-level ρ), split by CW/EW/P and Th1/Th2.
- CV_ISI from PDT data is limited to spikes in permissive PDTs.

## 7. Addendum (same day): CNPU inclusion rebuilt from boundary files

**What's in `boundary_info/SuperEmbryoBoundaries/SuperEmbryo_MatrixForm/`:**
- Only older (Jan 2021) single-group matrices, not the 2022 per-group `*_Only_IDs.csv` files.
- `1_SuperEmbryo_MatrixFormID.csv` has 41 embryo columns: 40 S13dev + 190320SP_E6, now Stage 14.
- Its ≥5-embryo hull gives a different grid origin from the paper's, and 4,351 of 4,738 S13dev CNPUs. Not used.

**Rebuild:**
1. Rasterise cell outlines from the boundary files (`Boundary_*` = Y coordinates then X).
2. Count embryos per pixel.
3. Take the convex hull of pixels with >4 embryos.
4. Keep centroids inside the hull.
5. Grid at the paper origin and keep CNPUs with ≥3 cells.

**Overlap with `AllGroups`:**

| Group | Overlap | Extras |
|---|---|---|
| S13 | 4,900/4,909 | 12 |
| S13dev | 4,613/4,738 | 2 |
| S14 | 5,487/5,501 | 2 |
| S14dev | 4,641/4,651 | 0 |

- Every paper CNPU has ≥3 cells.
- Using the origin from the rebuilt hull, S14 matches exactly; the other groups are offset by up to ~16 px.
- The prototype is the scratch script `rebuild_cnpu.py`, not in the repo.
