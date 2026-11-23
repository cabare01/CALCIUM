# Calcium III spectral analysis: project summary for poster and manuscript work

*Prepared 2026-11-02 from the analysis repository (`CALCIUM`) on the W&M HPC cluster. **Updated 2026-11-03** after an audit of the full PDT outputs and a reading of three sources: the manuscript's supplementary text, Sudip Paudel's PhD dissertation, and Sudip's original R scripts. **Updated 2026-11-04** with the results of the red-noise rhythmicity test (step 1) and the rhythmic-only maps (step 1c); see Section 6.8. This document is self-contained: the reader cannot open the cluster files named here, which are listed only so the author can find them.*

---

## 0. How to use this document (read first)

You are helping write (1) a **Society for Neuroscience (SfN) poster** and (2) the **spectral-analysis sections of a calcium-imaging manuscript**. The manuscript already contains event-detection methods: **PDT = peak detection technique**, from the earlier analysis in the same lab by Sudip Paudel. The draft is called *CB Calcium III Working*. The spectral results go in its Results §3.6, "Spectral analysis reveals distinct temporal characteristics of calcium activity beyond peak frequency", which is currently a heading only, plus a matching Methods subsection.

**Authors** (SfN abstract): Catherine Bare (submitting and presenting author), Sudip Paudel, Margaret Saha (senior author). The PDT analysis files originate from the `W:\SudipP\1Calcium III\…` project directory.

Ground rules:

1. **Do not invent results.** Use only numbers stated here. Where a result is pending, write a clearly marked placeholder such as `[RESULT: % rhythmic cells, step 1]`.
2. **Status labels used below:**
   - **Established:** computed, audited and safe to state.
   - **Preliminary:** computed, but with a caveat that must travel with it.
   - **Pending:** not yet computed; a script exists.
3. **Methods must parallel the PDT paper** wherever a choice exists. Section 4 lists the rules. New embryo-level statistics are *added* alongside the paper's methods and never replace them.
4. **Only claims supported with the embryo as the unit of analysis are findings.** Cell-level significance alone is not enough (Section 6.4).
5. **Items marked `[CONFIRM]`** are details not recorded in the analysis repository (for example, full citations). Ask the author for them rather than guessing.
6. **The SfN abstract was written before the current (v3) analysis.** Section 12 lists where it no longer matches the data. The poster must be accurate to the data while staying within the abstract's topic.
7. **Section 13 summarises the PDT study** (methods, codes, findings) from the manuscript, its supplement, the dissertation and the original code. Use it for background, citations and parallel methods text.

**What changed on 2026-11-04** (details in Section 6.8):
- **The red-noise (AR(1)) rhythmicity test has been run** on all 300,023 analysed cells from 179 embryos. Its integrity is audited.
- **"Rhythmic" is defined by the max-ratio global test (P_global), BH-FDR within embryo, q < 0.05.** This is the author's decision. The rationale is in 6.8.
- **Result (Preliminary):** rhythmicity is rare in intact embryos and common in devitellinated ones. The median embryo has 0% rhythmic cells (intact, both stages) vs 7.9% (S13 devitellinated) and 2.3% (S14 devitellinated).
- **Rhythmic cells share one spectral signature:** about 2–2.7× more power than red noise at **0.01–0.015 Hz**, and very high band power.
- **Correction (11-04, exploratory transient check):** these cells do **not** oscillate every 80–100 s. They carry a few large, slow calcium transients: median 4 per hour, about 26 s wide at half height, about 6 min apart, at irregular intervals (CV about 0.7). The 0.01–0.015 Hz excess reflects the *duration* of these transients, not a repeat period. **Do not write "quasi-periodic" or "80–100 s rhythm".** Write "slow, large calcium transients that exceed red noise" (6.8, caveat 3).
- **The 0.08 Hz dominant frequency of most intact cells is noise.** It is where a red-noise-like spectrum happens to be highest, which resolves the question in 6.2.
- **The rhythmic-only maps are too sparse to support spatial claims.** The median CNPU (composite neural plate unit) holds 1 cell from 1 embryo. **Rhythmic-fraction maps were built instead** (6.8). Intact hot spots are single-embryo. S13 devitellinated shows a multi-embryo rise in rhythmic % toward lateral non-neural tissue, but the embryo-level test does not support it yet.

**What changed on 2026-11-03** (details in Sections 6.6, 6.7 and 13):
- **The full PDT outputs are now on the cluster for all 179 embryos.** They join 100% of the 191,527 v3 mapped cells, so the PDT ↔ spectral comparison (step 5) can use all 177 mapped embryos, not 79.
- **Whether PDT event rate agrees with spectral power depends on the PDT** (Preliminary):
  - PDTs with embryo-wide thresholds ("EW") give a **positive** correlation with band power.
  - PDTs with per-cell thresholds ("CW") give a **negative** one in devitellinated embryos.
  - So step 5 must report across all 204 PDTs, plus a pre-specified headline PDT.
- **The PDT naming code is decoded** from the dissertation's Flowchart 4.1.
- **Template orientation is confirmed:** anterior = small template Y = top of the maps.
- **The paper's real CNPU inclusion rule** (from Sudip's code) differs from the text. It has been reproduced from the boundary files at 97–99.8% agreement.
- **The devitellination effect in the PDT study is PDT-dependent:** "fewer, longer-lasting events", not simply "more activity".
- **Sudip computed a spike CV with a non-standard formula.** Its output is not on the cluster.

---

## 1. Project in one paragraph

*Xenopus laevis* embryos expressing the calcium indicator GCaMP6m and membrane-targeted RFP (memRFP) were imaged at early neurula stages (stages 13 and 14), to record calcium activity in every tracked neural plate cell. The PDT (peak detection technique) analysis detected discrete calcium events per cell with 204 detection variants, classified them as spikes, rosettes or waves, and mapped them onto a composite neural plate ("SuperEmbryo" template). This project adds a **time–frequency (continuous wavelet) analysis** of each cell's ΔF/F trace. It asks **where in the neural plate activity is rhythmic, at what frequency and power, and whether this differs by developmental stage (13 vs 14) or by removal of the vitelline membrane (intact vs devitellinated)**. Spectral features are mapped onto the same template, with the same composite-neural-plate method as the PDT paper, so the two analyses can be integrated in one manuscript.

---

## 2. Experimental design and data

| Group label | Stage | Vitelline membrane | Embryos in reference table | Embryos mapped |
|---|---|---|---|---|
| S13 | 13 | intact | 46 | 44 (2 excluded, see 2.1) |
| S13_devitellinated | 13 | removed | 46 | 46 |
| S14 | 14 | intact | 46 | 46 |
| S14_devitellinated | 14 | removed | 41 | 41 |
| **Total** | | | **179** | **177** |

- **Design:** 2 × 2, stage × membrane. The embryo is the biological unit.
- **Group membership** comes from a curated reference table (`Embryo_Ref_Table.csv`), which is definitive. Some file names suggest a different stage; the author has confirmed these embryos are correctly grouped. Do not flag them.
- **Species and imaging protocol** (the author's text, usable in Methods):

  > Embryos that expressed memRFP and GCaMP6m were sorted using a dissecting scope (Nikon, SMZ800N). Each embryo was oriented with its dorsal side facing down in a 35 mm Glass Bottom Microwell Dish (MatTek Corporation, P35G-1.5-10-C) and then immobilized using a cover glass (Goldseal Cover Glass, Thermo Scientific, Reorder No 3306). Two regions of each neural plate that expressed both GCaMP6m (using FITC channel as described in Ablondi et al., 2020) and memRFP (TRITC channel, emission wavelength 595 nm, laser wavelength 561.5 nm) were imaged at 20X magnification using the multi-point acquisition module of Nikon Eclipse Ti-E CLSM with NIS Elements software for 60 min at the acquisition rate of 0.33Hz. In addition, 10X magnification snapshot images of each embryo were acquired at the beginning and end of each imaging session.

  - **Species:** *Xenopus laevis*.
  - **Staging:** the PDT supplement gives the timing:
    - embryos imaged at stage 13 and fixed at stage 14 were 14 h 45 min – 16 h 15 min post-fertilisation;
    - embryos imaged at stage 14 and fixed at stage 15 were 16 h 15 min – 17 h 30 min.

    This is consistent with Nieuwkoop–Faber timing, but the staging system is not named `[CONFIRM]`.
  - **Injections** (PDT Methods): 1 ng GCaMP6m mRNA + 0.5 ng memRFP mRNA.
  - **Nominal vs measured frame interval:** acquisition was nominally 0.33 Hz (3.0 s). The measured interval in the tracked data has a median of 3.57 s (range 2.5–4.5 s; 84% of cells > 3.4 s), most likely because two regions were acquired per cycle with multi-point acquisition. The spectral pipeline uses the measured per-embryo interval, so Methods should state both.
  - **Refocus pauses:** when the image blurred, acquisition was paused for 28.9 ± 14.6 s and the region refocused (PDT Methods).
  - The median usable tracking duration per cell is about 52–59 min of the 60 min recording.
- **Neural vs non-neural territory** is defined by a Sox2-based boundary on the template:
  - stage-13 groups are overlaid with the **Stage 14 Sox2 template**;
  - stage-14 groups are overlaid with the **Stage 15 Sox2 template**;
  - this matches the PDT paper's convention.

  The per-cell Sox2 signal is a **whole-mount Sox2 FISH** (DIG probe, Cy5 tyramide) measured after imaging, not an immunostain.
- **Template orientation (Established 11-03): anterior = small template Y = top of the maps;** the blastopore is at the bottom. Three checks agree:
  - the Sox2 template is widest at small Y and tapers toward large Y;
  - the template's normalised AP coordinate is 0 at the smallest Y;
  - the PDT supplement draws the forebrain at the top.
- **Recording year is unbalanced across groups** (embryos per group):

  | Group | 2018 | 2019 | 2020 |
  |---|---|---|---|
  | S13 | 31 | 7 | 6 |
  | S13_devitellinated | 10 | 15 | 21 |
  | S14 | 17 | 11 | 18 |
  | S14_devitellinated | 5 | 14 | 22 |

  This is a potential batch confound; see 6.5.

### 2.1 Cell counts through the pipeline (v3)

| Group | Cells tracked | Passed spectral QC | Mapped to template | Mapped cells inside CNPU grid |
|---|---|---|---|---|
| S13 | 148,341 | 78,647 | 49,575 (44 embryos) | 49,553 |
| S13_devitellinated | 127,977 | 68,413 | 46,409 | 45,986 |
| S14 | 159,424 | 83,052 | 53,243 | 52,405 |
| S14_devitellinated | 141,201 | 69,911 | 42,300 | 41,768 |
| **Total** | **576,943** | **300,023 (52%)** | **191,527 (33% of tracked)** | **189,712** |

- The PDT study's composite neural plate (CNP) holds **201,409 cells from 179 embryos**: every cell present on the last frame. The v3 mapped set (191,527) is the subset that also passed spectral QC.
- **Why cells fail spectral QC** (all groups): no gap-free segment longer than about 490 s (about 71% of skips); too little cone-of-influence-valid data at 0.005 Hz (about 14%); fewer than two complete observations (about 15%); under 100 frames (about 1%).
- **Why passed cells don't map:** they don't appear in the template boundary file (about 35% of passed cells).
  - **Established** (PDT Methods and dissertation): cell outlines come from the **last frame** of each recording, so only cells present at the end can be mapped. The maps therefore over-represent cells tracked until the end.
  - The PDT paper states the same limitation. A bias check is pending (step 2).
- **Excluded from mapping:** 180405SP_E4 and 200114SP_E1 (S13). Their boundary files were built from a different recording than the tracked calcium data, so cell IDs do not correspond (affine registration R² 0.26 and 0.09).
  - The PDT outputs for these two embryos contain only 36–38% of their boundary cells. The PDT study therefore used the same mismatched files, which supports the exclusion.
- **Registration QC for the rest:** an affine fit of image position to template position gives a median R² of 0.989. Three embryos fall below 0.9 (0.78, 0.84, 0.89); they are kept but noted.

---

## 3. Methods as implemented (manuscript-ready facts)

### 3.1 Preprocessing (shared with the PDT methods)
- Cells tracked in the merged tracking tables. A trace needs **≥ 100 observed frames** (the paper's tracking criterion).
- **Gaps of up to 3 missing frames** are bridged by linear interpolation. This mirrors the tracking software's 3-frame maximum gap.
- **Refocus pauses** split a trace into segments, which are processed separately. The point immediately before and after each pause is dropped. A pause is a frame step of 1 with a time jump above 1.5× the frame interval; observed pauses were 23–32 s.
  - The PDT code detected pauses as a >4 s gap in the microscope's timestamps and worked in frame numbers.
- **Time base:** an embryo-wide frame → time map, regressed within pause-free runs. This removes timestamp rounding jitter.
- **Baseline F0** = mean of the lowest 10% of intensities, per segment (the paper's baseline definition (a)). ΔF/F = (F − F0)/F0.
- Each segment is resampled to a uniform grid and **linearly detrended** to remove photobleaching drift.
  - The PDT pipeline ran on raw intensities. It computed 2nd/3rd-order polynomial detrends but did not use them in the final analysis.
- A segment must be longer than 2 × the cone of influence at the lowest frequency (about 490 s) to be used.

### 3.2 Continuous wavelet transform
- **Wavelet:** complex Morlet `cmor1.5-1.0` (PyWavelets 1.9.0, which self-reports as 1.8.0; Python 3.12).
- **Frequency axis:** **0.005–0.08 Hz**, 80 linear bins (periods 12.5–200 s).
  - The top stays below the lowest Nyquist frequency in the dataset (0.111 Hz) with margin for the wavelet's bandwidth.
  - Measured calcium events last about 8 frames (about 30 s) and carry little power above 0.06 Hz.
  - The bottom was lowered to 0.005 Hz because event rates are about 0.005–0.02 Hz.
- **Power** = |W|² / scale (scale-rectified; Liu, Liang & Weisberg 2007). Unrectified power grows with scale and had pushed 53–71% of cells to the lowest bin in v1.
- Each segment is transformed separately. **Cone-of-influence (COI) masking**: only samples at least √B·C/f seconds from a segment edge are used (COI ≈ 1.22/f s). The time-averaged spectrum pools COI-valid samples across segments.
- A cell is kept only if, at every frequency, it has **at least 1 full cycle of COI-valid data** (200 s at 0.005 Hz).
- **Validation:**
  - Synthetic sinusoids sampled at 3.3 s with rounding jitter, 3% dropped frames and a 30 s pause were recovered within one frequency bin.
  - Equal-amplitude sinusoids give roughly flat peak power across frequency.
  - White noise does not pile up at the lowest frequency.
  - These checks were run on the v2 axis (0.02–0.15 Hz); repeating them on the v3 axis is recommended.

### 3.3 Six spectral features per cell (from the time-averaged spectrum P(f))
| Feature | Definition |
|---|---|
| DominantFreqHz | Frequency of maximum P(f) |
| SpectralCentroid | Σ f·P / Σ P |
| SpectralSpread | √(Σ (f − centroid)²·P / Σ P) |
| SpectralEntropy | Shannon entropy (log₂) of P normalised to sum 1. Maximum for a flat spectrum = log₂80 = 6.32 |
| TotalBandPower | Σ P over 0.005–0.08 Hz |
| PeakPowerFraction | max P / Σ P |

### 3.4 Mapping to the composite neural plate (as in PDT paper Methods 2.3)
- Cells are joined to the template by cell ID and assigned to **9 × 9-pixel CNPUs** (composite neural plate units) on a fixed grid: origin (120, 1), X 120–1122, Y 1–1190. The same grid is used for every group.
- **Primary map (paper method):**
  1. Take the pooled mean of all cells from all embryos in each CNPU.
  2. Smooth with a 7 × 7 (49-CNPU) mean **over occupied CNPUs only**; empty CNPUs are not filled.
  3. Z-score within each group and feature.
  4. Display on a fixed **−1…1 jet colour scale** that saturates beyond ±1.
- Combined 2 × 2 maps share one colour scale across groups and **show only CNPUs supported by ≥ 2 embryos**.
- Secondary maps: an embryo-balanced median (median of per-embryo CNPU medians), SD, and coverage (% of a group's embryos contributing).
- Cells outside the grid (1,815, mostly from four laterally or posteriorly imaged embryos) are excluded from maps but kept in region statistics.
- **How the PDT paper's maps were actually built (Established 11-03 from Sudip's code).** Use this for parity and in Methods text that claims to "follow the PDT paper":
  - **Inclusion.** The text says "CNPUs with at least three embryos within a rubber-band boundary of n = 5 embryos". The code did this:
    1. Take template pixels covered by ≥5 embryos and draw their **convex hull** (the "rubber band").
    2. Keep cells whose centroid falls inside the hull.
    3. Keep CNPUs with **≥3 cells** (not ≥3 embryos).

    This leaves 4,909 / 4,738 / 5,501 / 4,651 CNPUs (S13 / S13dev / S14 / S14dev) and includes surrounding non-neural ectoderm.
  - **Grid.** The paper's grid starts at the lowest cell coordinate inside the hull, so **its origin differs per group and from the v3 grid**: a sub-CNPU offset, with a median 3.6–4.5 px to the nearest v3 CNPU centre. Cell-by-CNPU comparisons with the paper's maps need re-gridding.
  - **Smoothing.** Mean of the 49 nearest CNPU centres within 39 px. At map edges, where fewer exist, the found neighbours are repeated to make 49. Then z = (local mean − overall mean) / overall SD.
  - **Reproduction.** Rebuilding the hull from the cell outlines in the boundary files and gridding at the paper's origins reproduces 4,900/4,909, 4,613/4,738, 5,487/5,501 and 4,641/4,651 of the paper's CNPUs. A "paper-CNPU-set" version of the spectral maps is therefore possible; it is not yet in the pipeline.
  - **What to write.** The v3 maps follow the paper's method (9-px CNPUs, pooled mean, 49-CNPU smoothing, z-score, −1…1 jet). They differ in two documented ways: a fixed grid shared by all groups, and a ≥2-embryo display mask instead of the hull + ≥3-cells rule. Say so in Methods, or add the parity version as a supplementary figure.

### 3.5 Statistics
- **Paper-style (descriptive, matches PDT Methods 2.5):** two-sample K-S test plus Cohen's d on cells, Bonferroni across tests. "Significant" = p_adj < 0.05 **and** |d| > 0.2.
- **Embryo-level (inferential, added):**
  - Paired Wilcoxon signed-rank tests on embryo medians (neural vs non-neural; embryos with ≥ 10 cells per region), Bonferroni.
  - **Pending:** a stage × membrane permutation model on embryo medians (Freedman–Lane, Holm), plus a cell-level linear mixed model with a random intercept per embryo.

---

## 4. Rules the author has set (keep these in all text and figures)
- Follow the PDT paper wherever a method choice exists: 9-px CNPUs, pooled CNPU mean, 49-CNPU smoothing, z-score, −1…1 jet scale, K-S + Cohen's d + Bonferroni, ≥ 100-frame tracks, 3-missing-frame gap rule, refocus splits with edge points dropped, lowest-10% baseline. Present alternatives only as secondary or supplementary.
- Group membership comes only from the reference table.
- State results at the embryo level. Present cell-level paper-style statistics as descriptive and for parity with the PDT paper.
- **Settled analysis choices** (made by the author after the audits; do not re-propose):
  - the 0.005–0.08 Hz band with 80 linear bins;
  - linear detrend;
  - ≥1 COI-valid cycle per frequency;
  - the two excluded embryos.

---

## 5. Completed work (chronological)
1. **MATLAB → Python/Slurm port** of the original wavelet workflow (Oct 26–27). Audits of the spectral tables and of an extreme-power/gap issue.
2. **Full pipeline audit** (Nov 1). It found:
   - frequencies above Nyquist;
   - a scale bias in power that put most cells at the lowest frequency;
   - CWT run across trace splices;
   - dropped frames splitting traces;
   - mis-grouped embryos;
   - two boundary files that don't match their recordings;
   - a recording-year confound.
3. **Pipeline v2 → v3** (Nov 2): every fix in Section 3. v3 is the current production run, and every reference embryo has raw data.
4. **Group tables v3**, then **template mapping v3** for the four groups plus a combined step. The mapping was audited:
   - aggregation recomputed exactly;
   - group membership matches the reference table;
   - boundary files all matched to the correct embryo.
5. **Map fixes after the audit:** the smoother no longer fills empty CNPUs, and the combined maps hide CNPUs with < 2 embryos. Rerun and re-audited the same day.
6. **Analysis plan and scripts written for steps 1–6** (Section 7). Step 1 has now been run (item 8); steps 2–6 have not.
7. **PDT data audit and source review** (Nov 3; read-only). This covered the full 2022 PDT outputs (47 CSVs, about 52 GB), the manuscript supplement, the Paudel dissertation and Sudip's R scripts. Results are in Sections 6.6, 6.7 and 13. No spectral or map outputs were changed.
8. **Step 1, red-noise rhythmicity test** (Nov 3–4).
   - A 50-cell pilot, then the full run: 4 jobs, 1,000 surrogates per cell, about 2–6 h per embryo.
   - Integrity audit: every analysed cell tested exactly once; the pilot was reproduced exactly; the recomputed spectra equal the production spectra.
   - Then the merged tables, the choice of test (P_global), and the **rhythmic-only maps (step 1c)**. The maps were audited as an exact subset of the v3 maps. Results are in 6.8.

---

## 6. Current results and their status

### 6.1 Group-level spectral descriptors (Preliminary; read 6.2 first)
Median across embryos of embryo-median values (mapped cells), with IQR across embryos:

| Feature | S13 (n=44) | S13 devitellinated (n=46) | S14 (n=46) | S14 devitellinated (n=41) |
|---|---|---|---|---|
| DominantFreqHz | 0.080 [0.053–0.080] | 0.012 [0.008–0.069] | 0.080 [0.078–0.080] | 0.013 [0.009–0.074] |
| SpectralCentroid (Hz) | 0.050 [0.047–0.051] | 0.036 [0.026–0.049] | 0.051 [0.049–0.051] | 0.041 [0.028–0.050] |
| SpectralEntropy (max 6.32) | 6.20 [6.17–6.22] | 6.11 [5.71–6.17] | 6.18 [6.16–6.20] | 6.16 [5.87–6.18] |
| TotalBandPower | 0.0055 [0.0037–0.0091] | 0.0167 [0.0124–0.0364] | 0.0069 [0.0046–0.0095] | 0.0159 [0.0105–0.0350] |
| PeakPowerFraction | 0.020 [0.019–0.020] | 0.031 [0.021–0.054] | 0.020 [0.020–0.020] | 0.027 [0.020–0.048] |

**Apparent pattern:** devitellinated embryos have about 2–3× higher band power, lower dominant frequency and centroid, lower entropy and higher peak fraction at both stages. Their spectra are more structured, with more slow (about 0.01 Hz, about 100 s period) activity. Stage differences look small. **This has not yet been tested at the embryo level** (step 4b).

**How this relates to the PDT devitellination result** (Section 13.6): the PDT study found devitellination gives larger events (significant at stage 13) and "fewer long-lasting events", with the details depending on the PDT. Fewer, longer, larger events would shift power toward low frequencies and raise band power. So the spectral pattern is *consistent* with the PDT result, but it is not the same claim as "more events".

### 6.2 Major caveat: many spectra may be noise-like (central open question)
- **In intact embryos, 42–44% of cells have their dominant frequency at the upper band edge (0.08 Hz),** and the median embryo's dominant frequency *is* that edge. Another 6–14% sit at the lower edge (0.005 Hz).
- Median entropy (6.15–6.20) is close to the flat-spectrum maximum (6.32).
- So **for many cells, especially in intact embryos, DominantFreqHz probably does not mark a real oscillation.** It marks where a nearly flat or noise-dominated spectrum happens to be highest.
- **Consequences:**
  - "Intact embryos oscillate at 0.08 Hz" is **not** a valid claim.
  - The intact vs devitellinated frequency difference may partly reflect activity level or signal-to-noise rather than rhythm frequency.
  - The maps' anterior-high / posterior-low DominantFreqHz pattern may likewise reflect where activity is weak.
- **New supporting evidence (6.7):** the PDTs that count many small threshold crossings in low-variance cells (CW thresholds) correlate *positively* with SpectralCentroid and DominantFreqHz. That fits the view that high-frequency, flat spectra mark quiet, noise-dominated cells.
- **Resolution (step 1, done 11-04; Section 6.8):** the red-noise test confirms the concern.
  - Non-rhythmic cells still have their dominant frequency at the 0.08 Hz edge in 43–45% of intact and 28% of devitellinated cells. Rhythmic cells almost never do (0.1–0.4%).
  - Non-rhythmic intact cells have *less* power than red noise across 0.005–0.05 Hz and reach the red-noise level only at about 0.08 Hz.
  - **Present frequency features only for rhythmic cells.** For all cells, present power and structure features (TotalBandPower, PeakPowerFraction, entropy).
- **Why this matters for the manuscript:** the PDT Discussion and the dissertation argue that wavelet and Fourier methods suit periodic cells, whereas pre-synaptic neural plate traces "lack periodicity and predictability" (citing Marken et al. 2016). The red-noise test is the direct, quantitative answer to that claim, whichever way it comes out.

### 6.3 Spatial maps (Established as maps; interpretation Preliminary)
- Shared-scale smoothed z-maps exist for all six features and four groups, plus coverage maps.
- **DominantFreqHz:** in all four groups, the **anterior** neural plate (top of the maps; orientation now confirmed) is high (z ≥ +1) and the posterior and lateral regions are low.
  - The pattern repeats across four independent groups, which makes it the most robust *descriptive* observation.
  - But see 6.2 and 6.7: it may reflect activity or signal-to-noise rather than frequency.
  - **Update 11-04 (6.8):** in intact embryos, 98–99% of cells are not rhythmic, and their dominant frequency marks the noise-dominated top of the band. So an "anterior-high DominantFreqHz" map most likely shows *where cells are quiet*, not where they oscillate faster. Do not label it a frequency gradient.
  - It has not yet been tested per embryo (step 4a).
- The map edges, where only 2–4 embryos contribute, are still often saturated; interpret the cores of the maps.
- **Context from the PDT study:**
  - PDT maps found **no clusters of high- or low-spiking CNPUs inside the neural plate**; cells there were "equally capable of spiking".
  - Wave patterns had an anterior–posterior bias that flipped with the PDT: posterior with EW or Th2 PDTs, anterior with CW or Th1 PDTs.
  - A spectral gradient would be new, but it should be presented alongside this PDT-dependence (6.7).

### 6.4 Neural vs non-neural (Preliminary, none supported yet)
- About 66–72% of mapped cells fall inside the Sox2 neural boundary.
- **Cell-level paper rule** (K-S Bonferroni p < 0.05 and |d| > 0.2) is met in five tests:
  - S13_devitellinated: SpectralCentroid (d = +0.25), SpectralEntropy (+0.23), TotalBandPower (−0.22), PeakPowerFraction (−0.20);
  - S14_devitellinated: DominantFreqHz (+0.22).
- **None of these five survives the embryo-level paired Wilcoxon test.** The three embryo-level results that do survive Bonferroni fail the paper rule (|d| < 0.2): S13 SpectralCentroid, S14 SpectralCentroid, S14 SpectralSpread.
- **No neural vs non-neural difference is currently supported at both levels.** With about 50,000 cells per group, every K-S p-value is effectively 0, so cell-level tests overstate significance. A concordance table is pending (step 4c).

### 6.5 Recording-year confound (Preliminary)
- Embryo-median TotalBandPower rises with recording year inside the intact groups. For example, S13 by year: 0.0044 (2018), 0.0070 (2019), 0.0090 (2020).
- Devitellinated groups are mostly 2019–2020 and S13 mostly 2018.
- **However, the devitellinated > intact power difference and the lower dominant frequency hold within each year.** For example, in 2018, S13 is 0.0044 and S13 devitellinated is 0.0121. So the membrane effect is not explained by year alone.
- The final models should include recording date or year (planned as an addition to step 4b; see 7). Of 67 recording dates, 33 span more than one group, so a date random effect is estimable.

### 6.6 PDT outputs now available (audited 2026-11-03; Established unless marked)

This replaces the 11-02 inventory. That inventory knew only the partial `ConvergeSpikeCy5` files (79 embryos), which are now superseded.

**A. Per-cell PDT tables** (one file per group × PDT set × threshold type; 20 files).
- **Coverage:** all **179 embryos** (46/46/46/41) and **201,409 cells**, exactly the manuscript's composite-neural-plate count. One row per last-frame cell. Each row carries the cell's template centroid and Sox2 (Cy5) intensity, identical to the boundary files.
- **Join to the spectral data:**
  - **100% of the 191,527 v3 mapped cells have a PDT row.** Template centroids are identical in 175 of 177 embryos.
  - **190213SP_E2:** the PDT study used a slightly different registration (median offset 5 px).
  - **200217SP_E3:** a 98.7% match, from duplicate cell IDs.
  - **Step 5 can use all 177 mapped embryos.**
- **Columns per PDT:**
  - **Scaling:** every count is scaled to the full recording, i.e. raw count × (total frames / frames the cell was tracked). The recording is ≈60 min, so "per hour" is approximately events per hour, computed in frames.
  - `Total`: the number of distinct events (spikes, rosettes or waves) the cell took part in.
  - `OnlySpikes` / `OnlyRosettes` / `OnlyWaves`: the same, split by event type.
  - `ActivePeaksPerHr`: the number of above-threshold runs; about 1.36× `Total` (why is `[CONFIRM]`).
  - `CellWiseActiveIncidence`: % of the cell's frames that were active.
  - These are the columns the paper's spatial tests used.
- **Name pitfalls:** six PDT embryo names need explicit aliases. Two lack the embryo number (`191030SPRN_`, `191020RNSP_`) and two use a hyphen (`200313SP-E1`, `-E2`). The analysis code handles this.

**B. Event-level PDT tables** (one row per detected event; 0.8–4 GB each).
- **Columns:** start and end frame, duration, number of participating cells, number of leader cells, event type (1 spike / 2 rosette / 3 wave), and the position of the event's leader cells.
- **Each row names only one cell**, a leader cell at the event's first frame. Rosette and wave participants are not listed.
- **Per-cell event times are therefore available only for spikes.** At a typical PDT only about 6% of cells have ≥3 spike events (range 0.3–96% across PDTs). So **CV_ISI from these files is feasible only for spikes under the most lenient PDTs.**
- **Inventory gap:** the literature PDTs have no cell-wise (CW) event files for the intact S13 and S14 groups.

**C. Paper map and summary tables.**
- Group-averaged PDT maps for all PDT groups (the 27 retained groups are identifiable; see 13.5).
- Per-embryo event counts and incidence for all 204 PDTs.
- The paper's spatial-pattern (NNDist) test tables.

**D. Sudip's original R scripts.** These define every column above; Section 13 summarises them.

**E. Still not on the cluster** (on Sudip's lab drive, `W:\SudipP\1Calcium III\…`; worth requesting from Sudip or Margaret):
- the list of which PDTs belong to each group (`MethodGroups_RealData.csv`);
- the per-PDT smoothed maps;
- Sudip's spike-CV results;
- the per-embryo, per-frame event files, which would allow per-cell event times for every event type;
- the a priori neural plate region polygons (`All_NeuralPlate_Regional_Boundaries.csv`).

**F. Older material (superseded for integration):**
- `ConvergeSpikeCy5` per-cell summaries (79 embryos) from an early, pre-204-PDT pipeline. `Counts_0` = all "sandwich" spikes; `Counts_M` = spikes larger than the cell's smallest; `Th_SD_k` = spikes above the cell minimum + k × the embryo's SD of spike size. All are scaled to the recording length.
- Per-embryo Sox2-vs-activity correlations (`1_R and C.csv`).
- A Jan 2021 pixel-level matrix for an early 41-embryo S13dev build.

### 6.7 PDT ↔ spectral agreement depends on the PDT (Preliminary, 2026-11-03)

**Method:** a quick descriptive check, not the step 5 analysis.
- For each embryo, the Spearman ρ across mapped cells between each PDT's `ActivePeaksPerHr` and each spectral feature.
- Summarised as the median across embryos (177 embryos).

| Feature | Change-in-magnitude PDTs, cell-wise SD (36 PDTs) | Change-in-magnitude PDTs, %-of-baseline "P" (36) | Literature PDTs, cell-wise (54) | Literature PDTs, embryo-wise (42) |
|---|---|---|---|---|
| TotalBandPower | −0.22 | −0.10 | −0.08 | **+0.39** |
| PeakPowerFraction | −0.13 | −0.08 | −0.09 | +0.15 |
| SpectralCentroid | **+0.38** | +0.20 | +0.28 | **−0.31** |
| DominantFreqHz | +0.29 | +0.16 | +0.16 | −0.23 |

*Values are the median over the PDTs in each set of the embryo-median ρ. Individual PDTs range from about −0.8 to +0.8 for every feature.*

**Methods-matched example PDTs:** literature PDTs with no preprocessing, the lowest-10% baseline (the same baseline as the spectral pipeline) and an SD threshold. Paper codes `2-111` (1 SD) and `2-112` (2 SD). Embryo-median ρ with TotalBandPower:

| PDT | S13 | S13 devitellinated | S14 | S14 devitellinated |
|---|---|---|---|---|
| `2-112.EW` (embryo-wide SD, 2 SD) | +0.39 | +0.71 | +0.45 | +0.66 |
| `2-111.EW` (embryo-wide SD, 1 SD) | +0.47 | +0.51 | +0.48 | +0.55 |
| `2-111.CW` (cell-wise SD, 1 SD) | +0.03 | −0.64 | −0.01 | −0.52 |
| `2-112.CW` (cell-wise SD, 2 SD) | −0.08 | −0.64 | −0.10 | −0.56 |

Background subtraction vs none (codes `1-…` vs `2-…`) barely changes these values; **CW vs EW sets the sign.**

**Interpretation (Preliminary):**
- **EW PDTs** apply one absolute threshold per embryo, so they count large transients. Those carry band power, hence the positive ρ.
- **CW PDTs** scale the threshold to each cell's own variability, so quiet cells register many small crossings. Those cells have flat, high-centroid spectra.
- This mirrors the PDT study's own findings:
  - the wave AP bias flips between CW/Th1 and EW/Th2 PDTs;
  - in simulations, CW/Th1 PDTs favour low-count ("REV") patterns and EW/Th2 favour high-count ("FWD") patterns.

**Consequences for writing:**
- "Spectral power tracks PDT-detected events" is true only for **embryo-wide-threshold PDTs**. Write it that way.
- Step 5 must report the **distribution of agreement across all 204 PDTs**, split by threshold type and Th1/Th2. This follows the PDT supplement's own "% of PDTs" reporting style.
- It should also name a **pre-specified headline PDT set**. Suggested: `2-111`/`2-112` in all four CW/EW × Th1/Th2 combinations; this choice is the author's `[CONFIRM]`.
- **Never choose a single PDT after seeing results:** the sign can be made either way.

### 6.8 Red-noise rhythmicity test: results (run 2026-11-03/04)

**What the test does (Methods-ready).**
- **Question:** for each cell, does its time-averaged wavelet spectrum contain more power than a red-noise process would give? Red noise here is a first-order autoregressive, AR(1), process: random fluctuations with short memory and no rhythm.
- **Fitting the null:** the AR(1) coefficient is estimated from the cell's own detrended ΔF/F with the bias-reduced form of Torrence & Compo (1998): α = (r₁ + √r₂)/2, clipped to [0, 0.99].
- **Surrogates:** 1,000 AR(1) series are simulated per cell, matched segment by segment in length, frame interval, linear detrend and standard deviation.
- **Same processing as the data:** every surrogate goes through the identical CWT, scale rectification and cone-of-influence masking.
- **Rhythmic** = the global test P_global, BH-FDR within each embryo, q < 0.05.
  - P_global is the Monte Carlo p-value of the largest observed / mean-surrogate power ratio over the 80 frequencies.
  - Taking the maximum over frequencies corrects for testing 80 frequencies.
  - The minimum possible p is 1/1001.
- This is the Monte Carlo form of the Torrence & Compo red-noise significance test.

**Integrity (Established).**
- All 300,023 cells passing spectral QC, from all 179 embryos, were tested exactly once; there are no missing values.
- The spectra recomputed inside the test equal the production spectra.
- A 50-cell pilot was reproduced exactly.
- The two embryos excluded from mapping are kept in all rhythmicity tables, flagged, so their results are ready if their mapping files are recovered.

**Why P_global rather than the alternatives** (the author's decision, 11-04):
- **P_global measures effect size.** It works in units of the red-noise level: "this cell has 3× the power red noise would give somewhere in the band". That means the same thing at every frequency and in every cell.
- **The alternative max-T test (P_maxT) measures surprise.** It works in units of the null's *spread* at each frequency. Toward 0.08 Hz the spread is very small, because many cycles are averaged: the null's 95th percentile is 1.27× its mean at 0.08 Hz vs 1.72× at 0.005 Hz.
  - So a small, consistent mismatch between data and model at 0.08 Hz becomes "highly significant", e.g. a measurement-noise floor that AR(1) lacks.
  - In the data, 51% of P_maxT hits in S14 sat exactly at 0.08 Hz. Almost all came from two embryos imaged at a 4.4–4.5 s frame interval instead of about 3.5 s, an artifact.
  - For P_global, 0.0–0.1% of hits peak above 0.05 Hz.
- **P_global's known weakness doesn't matter here.** It has low power for rhythms at 0.03–0.08 Hz. But the excess in the data is at 0.01–0.015 Hz, and once P_maxT's edge artifacts are removed the two tests agree to within about 1 percentage point per group.
- **A third column is not a valid test.** "Dominant peak above the pointwise 95% band" compares the highest of 80 peaks with a cutoff meant for a pre-chosen frequency, so 22–35% of cells "pass".

**Results (cell counts Established; interpretation Preliminary).**

| | S13 | S13 devitellinated | S14 | S14 devitellinated |
|---|---|---|---|---|
| Cells tested | 78,647 | 68,413 | 83,052 | 69,911 |
| Rhythmic, pooled cells | 1.4% (1,132) | 10.6% (7,231) | 1.4% (1,176) | 7.2% (5,041) |
| Rhythmic, **median embryo** | **0%** | **7.9%** | **0%** | **2.3%** |
| Embryos with ≥ 1 rhythmic cell | 5 / 46 | 33 / 46 | 7 / 46 | 22 / 41 |
| Rhythmic among mapped neural / mapped non-neural / unmapped cells | 1.9 / 2.3 / 0.4% | 11.9 / 17.3 / 4.4% | 2.2 / 1.2 / 0.6% | 9.5 / 11.8 / 2.6% |

- **Membrane effect.** Devitellinated > intact at both stages.
  - The embryo-level Mann–Whitney comparison of per-embryo % rhythmic is p = 2 × 10⁻⁸ at stage 13 (n = 44 vs 46) and p = 1 × 10⁻⁴ at stage 14 (46 vs 41).
  - This is a quick exploratory test; the planned stage × membrane permutation model with recording year (step 4b) is the confirmatory one.
  - The effect also holds among cells with similar α. For example, at α 0.90–0.95, rhythmic cells are 7.4% vs 23.1% (stage 13) and 5.4% vs 17.6% (stage 14).
  - It fits 6.1 (more slow, high-power activity in devitellinated embryos) and the PDT finding of "fewer, longer-lasting events" after devitellination (13.6).
- **Intact embryos.** Rhythmicity is essentially absent at the embryo level: the median is 0%, and 5–7 embryos per group carry all rhythmic cells. A few outlier embryos (e.g. one S13 embryo with about 37% rhythmic cells) dominate the pooled intact percentages.
- **Rhythmic cells' spectral signature** (all four groups alike):
  - **Shape:** a band-limited excess over red noise of about 1.9–2.7× (group means) at 0.01–0.015 Hz, falling *below* red noise above about 0.03 Hz.
  - **Where the excess peaks:** median 0.010–0.013 Hz. This is a spectral position that reflects transient width, not a repeat period (caveat 3). 96–97% of rhythmic cells peak at 0.005–0.02 Hz; none are at the band edges.
  - **Features vs non-rhythmic cells:** median TotalBandPower 0.17–0.27 vs 0.006–0.013 (about 20–40× higher). Entropy 5.1–5.2 vs 6.2. PeakPowerFraction 0.06–0.08 vs 0.02. DominantFreqHz 0.008–0.009 Hz.
  - **Largest power ratio:** a median of about 3× red noise.
- **Where rhythmic cells are** (exploratory, embryo level):
  - **Neural vs non-neural:** no difference in the rhythmic fraction (S13 devitellinated: 15/32 embryos higher in neural tissue, p = 0.45; S14 devitellinated: 11/20, p = 0.39).
  - **Front to back:** in the S13 devitellinated neural plate, the posterior third has a higher rhythmic fraction than the anterior third in 12 of 14 embryos (median 24.6% vs 5.6%; Wilcoxon p = 0.004). S14 devitellinated points the same way (4/6, p = 0.44).
  - This was found after looking and uses only embryos with cells in all three thirds. Treat it as a hypothesis to pre-specify for step 4a, not a finding.
- **Selection.** Unmapped cells (those not present on the last frame) are 3–5× less often rhythmic than mapped cells, so mapped percentages are higher than whole-population percentages.

**Caveats that must accompany any rhythmicity statement.**
1. **The null is conservative.** It is variance-matched to the cell, but in the 50-cell pilot, real cells had a median 53% of their variance below 0.005 Hz, against 19% in AR(1) surrogates. So the null predicts *more* in-band power than typical cells have: 36–49% of cells have p > 0.9. **Rhythmic percentages are lower bounds.**
2. **Detection depends on amplitude and α.** 62% of low-α cells (α ≤ 0.5) have p > 0.9, against 14% of high-α cells (α > 0.95). "Not rhythmic" means **"not detectable"**, not "aperiodic", especially for low-amplitude cells.
3. **Rejecting red noise is not proof of oscillation, and here it is not oscillation** (exploratory check, 11-04).
   - **Sample:** about 1,400 sampled cells from 6 devitellinated embryos (3 per group, chosen as closest to their group's median % rhythmic) and the one high-rhythmic intact embryo. Transients were counted with a simple prominence rule: ≥ 6× the robust noise SD and ≥ 0.1 ΔF/F.
   - **Rhythmic devitellinated cells (n = 515):** median 4 transients per hour (IQR 3–6), median width 26 s at half height, median interval 6.2 min. The intervals are irregular: median CV_ISI 0.72, and only 5% of cells with ≥ 6 events have a CV below 0.5.
   - **Non-rhythmic cells (n = 900):** median 2 transients per hour, width 20 s, interval 6.8 min, CV 0.67.
   - **Meaning:** "rhythmic" here means *more, larger and broader slow transients than red noise predicts*, not periodicity. The 0.01–0.015 Hz power excess corresponds to transient duration.
   - **Wording:** "slow, large calcium transients that exceed red noise" or "structured activity beyond red noise". Not "quasi-periodic", "80–100 s rhythm" or "oscillators".
   - **Files:** script `CALCIUM/scripts/step1_rhythmicity/transient_check.py`; output `analysis_v3/step1_rhythmicity/transient_check/transients.csv`. This is a quick rule, not a validated event detector.
4. **Weak power for faster rhythms.** P_global has little power at 0.03–0.08 Hz, so faster rhythms would be under-detected. The absence of rhythmic cells there is uninformative.
5. **Cell-level p-values are descriptive.** Group claims rest on per-embryo rates.

**Rhythmic-only maps (step 1c, audited).** The v3 mapping was rerun on rhythmic cells only.
- **Integrity:** the output is an exact subset of the v3 maps: identical cells, coordinates, CNPUs and feature values. Mapped: 1,008 / 6,259 / 993 / 4,337 cells from 5 / 33 / 7 / 22 embryos.
- **Coverage is too sparse for spatial conclusions:**

  | | S13 | S13 dev | S14 | S14 dev |
  |---|---|---|---|---|
  | CNPUs with ≥ 2 embryos (shown in combined maps) | 56 | 1,371 | 45 | 708 |
  | All-cell v3 maps, same rule | 5,724 | 5,867 | 6,638 | 6,093 |
  | Max embryos in any CNPU | 2 | 6 | 3 | 5 |

  - The median occupied CNPU holds 1 cell from 1 embryo.
  - **Intact maps** reflect 1–3 embryos; one S13 embryo supplies 61% of its rhythmic cells.
  - **Devitellinated maps:** the hot and cold patches are driven by a few embryos. The top embryo supplies about half the cells in each patch's smoothing window, and 4 embryos lead 73–88% of patches.
  - Also, these maps z-score features *among rhythmic cells*, which doesn't show where rhythmic cells are concentrated.
- **For the poster and manuscript:**
  - Don't use the rhythmic-only feature maps as main figures; at most, show them as supplementary.
  - The informative spatial figure is the **rhythmic-fraction map** below (built 11-04).
  - Within rhythmic cells, one neural vs non-neural feature difference passes the embryo-level test (S13 devitellinated SpectralSpread), but its size is negligible (0.0005 Hz).

**Rhythmic-fraction maps (step 1d, built 2026-11-04): where rhythmic cells are concentrated.**

*Method (Methods-ready).*
- For each CNPU: the fraction of **all** mapped cells (from all embryos) that are rhythmic. That is the paper's pooled CNPU mean applied to a 0/1 rhythmic indicator.
- Then the paper's 7 × 7-CNPU smoothing over occupied CNPUs and a z-score within each group, shown on the −1…1 jet scale. Combined maps show CNPUs supported by ≥2 embryos.
- Every CNPU therefore keeps the full embryo support of the all-cell maps.
- The code reproduces the production v3 maps exactly when applied to an existing feature (built-in check).
- 189,703 mapped cells inside the grid.

*Extra panels:*
- **A smoothed "% of cells rhythmic" map on one shared scale (0–44%).** z-scores are relative within a group, so this is the only version that compares groups.
- Support QC maps.
- An embryo-balanced version, which is redundant: r ≥ 0.993 with the pooled map, because an embryo contributes a median of 1 cell per CNPU.

| | S13 | S13 devitellinated | S14 | S14 devitellinated |
|---|---|---|---|---|
| CNPUs shown (≥2 embryos) | 5,724 | 5,867 | 6,638 | 6,093 |
| Median smoothed % rhythmic (90th percentile) | 0.4% (5.9%) | 12.2% (30.1%) | 0.3% (6.3%) | 8.4% (19.0%) |
| CNPUs whose rhythmic cells come from ≥3 embryos (7 × 7 window) | 9% | 77% | 7% | 51% |
| "Hot" CNPUs (z > 1): median embryos contributing rhythmic cells; top embryo's share | 1; 100% | 3; 56% | 1; 100% | 1; 100% |

*Interpretation (Preliminary; descriptive maps, cell-pooled).*
- **Intact embryos.** About 0% rhythmic almost everywhere. Every z > 1 hot spot comes from a single embryo. **Do not read or show the intact z-maps as spatial patterns.** On the shared % scale they are uniformly pale, which is the honest picture.
- **S13 devitellinated** (best-supported map). Rhythmic cells are widespread and supported by multiple embryos. The rhythmic fraction rises from the centre of the neural plate outward:

  | Location | Median smoothed % rhythmic |
  |---|---|
  | Deep inside the Sox2 boundary (> 50 px) | 10.6% |
  | Just inside | 13.3% |
  | Just outside | 15.6% |
  | Lateral non-neural tissue (> 50 px outside) | 28.7% |

  - 88% of multi-embryo hot spots lie outside the Sox2 boundary, mainly mid-AP and lateral.
  - **But** the embryo-level paired test of neural vs non-neural rhythmic fraction is **not significant** (p = 0.45). Lateral tissue is covered by fewer embryos.
  - So "rhythmic cells concentrate lateral to the neural plate" is a map-level observation and a hypothesis for step 4a, not a finding.
  - Within the neural plate, the AP thirds are similar (10.9 / 11.5 / 13.2%, anterior → posterior).
- **S14 devitellinated.** Weaker structure. The anterior neural plate is lowest (6.0% vs 15.4% mid and 12.6% posterior). Only 3% of CNPUs have z > 1, nearly all single-embryo speckles at the map edge.
- **Leads to pre-specify for step 4a:**
  - lateral/non-neural enrichment (S13 devitellinated);
  - posterior > anterior (S13 devitellinated, 12/14 embryos by AP thirds per embryo, found post hoc);
  - an anterior-low neural plate (S14 devitellinated).

---

## 7. Planned analyses (step 1 run; steps 2–6 scripts written, not run)

| Step | Question | Output used for | Status |
|---|---|---|---|
| 0 | Housekeeping: git commit; install statsmodels. ~~Confirm template orientation~~ (done: anterior = top). ~~Event-column definitions~~ (done; Section 13) | — | Partly done |
| 1 | **Which cells are rhythmic?** Per-cell Monte Carlo AR(1) red-noise test (1,000 surrogates through the identical CWT/COI pipeline; global p from the max ratio over frequencies; BH-FDR within embryo, q < 0.05). Then maps of rhythmic cells only | % rhythmic per group and region; rhythmic-only maps; example spectra vs red-noise band | **Done 11-04** (6.8). Rhythmic-only maps too sparse; **rhythmic-fraction maps built** (step 1d) |
| 2 | Cell flow, and whether skipped cells sit in different template positions from kept ones (selection bias) | Cell-flow figure; Methods counts; limitations | Pending |
| 3 | Which features are redundant (within-embryo and between-embryo Spearman ρ) → choose 2 for the main figures | Feature choice | Pending |
| 4a | Anterior–posterior and medio-lateral gradient per embryo (Spearman ρ with position in the neural plate), tested across embryos; stage × membrane effects on the gradient. Also consider the paper's a priori AP regions (forebrain, midbrain, hindbrain, spinal cord) and DV regions for parity | Main quantitative claim about the maps | Pending |
| 4b | Stage × membrane effects on embryo medians (permutation model, Holm), planned pairwise contrasts, paper-style cell table, mixed model. **Add recording year/date as a covariate or random effect** (not yet in the script) | Main group-comparison claims | Pending |
| 4c | Region concordance: label each neural vs non-neural result Supported / Cell-level only / Embryo-level only / Not supported | Which region claims can be stated | Pending |
| 5 | Link to PDT events. Per-embryo correlation of event rate with spectral features **across all 204 PDTs** plus the pre-specified headline set (6.7); events in rhythmic vs non-rhythmic cells; event maps on the same grid and per-embryo spatial agreement | Integration section; validation panel | Pending. **Inputs now complete for all 177 mapped embryos**; scripts need adapting to the new per-cell tables |
| 6 | Figure panels: cell flow; example spectra; 2 × 2 maps; % rhythmic map; per-embryo dot plots (x = stage, colour = membrane); events vs spectral | Poster and manuscript figures | Pending |

**Decision rule after step 1:**

| % of cells rhythmic | Headline |
|---|---|
| ≥ about 30% | Frequency features are meaningful (rhythmic cells in the main figures) |
| about 5–30% | Rhythmicity is the headline (% rhythmic and where); frequency given for rhythmic cells only |
| < about 5% | Headline is power and structure; the absence of rhythm is stated explicitly. This would agree with the PDT study's "lack periodicity" view, quantitatively |

**Outcome (11-04): the result splits by membrane.**
- **Intact embryos fall in the < 5% row:** median embryo 0%, pooled 1.4%. This quantitatively supports the PDT study's "lack periodicity" view for intact embryos, with the conservative-null caveat.
- **Devitellinated embryos fall in the 5–30% row:** median embryo 2.3–7.9%, pooled 7–11%.
- **Suggested headline (the author's call `[CONFIRM]`):** the membrane contrast in rhythmicity, with frequency features reported only for rhythmic cells, and power and structure features for all cells.

---

## 8. Poster plan (SfN)

**Fixed details:**
- **Title (from the abstract):** *Decoding Calcium Activity in Early Neural Development using Wavelet Analysis*.
- **Category:** A. Development.
- **Keywords:** Calcium, Wavelets, Xenopus.
- **Authors:** Catherine Bare, Sudip Paudel, Margaret Saha.
- **Board:** 6 ft wide × 4 ft high (1.8 m × 1.2 m), landscape. A 72 × 48 in layout works well with 4 columns, a title band about 6–8 in tall, body text about 28–32 pt, panel headings about 44 pt, title about 90 pt.

**Single message (draft, to be finalised after steps 1 and 4):** *"Calcium activity across the neural plate has spatially organised spectral structure; removing the vitelline membrane increases slow, high-power activity at both stages; spectral power tracks events detected with embryo-wide thresholds, and the agreement depends on how events are detected."* Each clause becomes a claim only once its embryo-level test (steps 1, 4, 5) supports it; otherwise rewrite it.

**Suggested revision after step 1** (the author's call `[CONFIRM]`):
> *"Most neural plate calcium activity is indistinguishable from red noise, but removing the vitelline membrane increases the share of cells whose slow, large calcium transients exceed red noise, at both stages."*

*(Corrected 11-04: the earlier draft's "quasi-periodic, about 80–100 s component" is not supported. Rhythmic cells have irregular transients about 6 min apart; see 6.8 caveat 3. The poster's Figure 2 title and conclusions need the same change.)*

- The step 1 clause is now supported by the exploratory embryo-level test. The stage × membrane model (step 4b) is still needed to state it as confirmed.
- The "spatially organised" clause is **not yet supported**: the rhythmic-only maps are too sparse (6.8), and the all-cell frequency maps mostly reflect noise (6.3).
- The PDT clause is unchanged (6.7; step 5).

**Layout** (about 6 panels, read left to right):
1. **Question and background:**
   - neural plate calcium activity in *Xenopus*; GCaMP6m imaging;
   - what peak detection (PDT) captures and misses (timing structure, rhythmicity, multiscale activity);
   - PDT findings to cite (Section 13.6): results depend strongly on the PDT (up to about 10-fold); no clusters within the neural plate; wave AP bias that depends on the PDT; devitellination effects.
   - A useful hook: the PDT literature survey (254 papers) found no PDT combination used by more than 2.4% of studies.
2. **Approach:** schematic of tracking → ΔF/F → wavelet → red-noise test → template mapping. Include the cell-flow bar (577k tracked → 300k analysed → 191k mapped; the composite neural plate is 201k last-frame cells).
3. **Validation and rhythmicity** (data now available, 6.8):
   - an example rhythmic and an example non-rhythmic spectrum, each plotted against its red-noise mean and 95% band;
   - the average observed / red-noise ratio curve: rhythmic cells show about 2–2.7× at 0.01–0.015 Hz; non-rhythmic cells sit below red noise;
   - **% rhythmic per embryo** as one dot per embryo (x = stage, colour = membrane), with n embryos;
   - the conservative-null caveat in the caption.
4. **Main maps:** 2 × 2 (stage × membrane) smoothed z-maps. Use the paper's −1…1 jet scale and the Sox2 outline, anterior at the top.
   - **Updated 11-04:** don't use rhythmic-only feature maps as main panels; they're too sparse (6.8).
   - Use all-cell maps of a power or structure feature, e.g. TotalBandPower.
   - **Rhythmicity panel (built, 6.8):** the 2 × 2 **shared-scale "% of cells rhythmic" map**, with n embryos, mapped cells and rhythmic cells in each panel title. This is the honest cross-group comparison: intact panels are pale, devitellinated panels structured.
     - Optionally add the paper-style z-map for S13 devitellinated only.
     - Don't show intact z-map hot spots; each is a single embryo.
     - Caption: "descriptive, cell-pooled; region differences not yet tested at the embryo level" (S13 devitellinated lateral enrichment: embryo-level p = 0.45).
   - Avoid an all-cell DominantFreqHz map, or label it as reflecting quiet, noise-dominated cells (6.3).
5. **Quantification:** one dot per embryo, median bar (x = stage, colour = membrane: blue intact, orange devitellinated). Show the gradient strength and the group effects, with permutation p-values.
6. **Events vs spectral:**
   - the PDT-dependence of agreement: distribution of embryo-level ρ across 204 PDTs, CW vs EW, plus the headline PDT set;
   - then conclusions, then limitations:
     - selection toward end-of-recording cells;
     - 2 embryos excluded for registration;
     - recording-year imbalance;
     - many spectra noise-like;
     - the red-noise null is conservative and less sensitive for low-amplitude cells ("not rhythmic" = "not detectable");
     - rejecting red noise ≠ strict periodicity;
     - CV_ISI limited to spikes (if shown).

**Figure conventions:** one dot per embryo; state n embryos (and cells) on every panel; avoid dual axes; keep the jet scale only for the paper-style z-maps.

---

## 9. Manuscript plan (spectral sections inside the PDT draft)
- **Methods:** a new subsection after event detection, "Time–frequency analysis", using Section 3 nearly verbatim:
  - preprocessing shared with PDT;
  - CWT parameters and rationale for the frequency band;
  - scale rectification and COI;
  - the six features;
  - the red-noise test (Torrence & Compo 1998, Monte Carlo form), described as in 6.8 "What the test does": the AR(1) fit, 1,000 variance-matched surrogates through the identical pipeline, P_global, and BH-FDR within embryo;
  - mapping by the paper's composite-neural-plate method, stating the grid and display-mask differences (3.4);
  - statistics: paper K-S/d plus embryo-level permutation model, mixed model, Holm/BH;
  - the PDT ↔ spectral comparison across all 204 PDTs.
- **Results** (§3.6): spectral results follow the event results and tie back to them:
  1. spectral power agrees with event rate for embryo-wide-threshold PDTs, and the agreement depends on the PDT (validation; 6.7);
  2. how many cells are rhythmic, and where. Numbers are in 6.8: rare in intact embryos, common after devitellination, with a shared signature of slow, large, irregular transients (about 26 s wide, about 4 per hour; 6.8 caveat 3). "Where": the rhythmic-fraction maps (6.8) are descriptive, with S13 devitellinated rising toward lateral non-neural tissue; an embryo-level test is pending (step 4a);
  3. spatial organisation (gradient);
  4. stage and devitellination effects, related to the PDT devitellination result (13.6);
  5. neural vs non-neural (only Supported rows).
- **Discussion:** engage directly with the PDT Discussion's claim that neural plate traces "lack periodicity" (6.2).
  - The red-noise result *supports* that claim for intact embryos.
  - After devitellination, more cells show slow transients that exceed red noise. But these are irregular (CV_ISI about 0.7), so they also *lack periodicity*. The spectral layer refines the claim rather than contradicting it.
  - State the conservative-null caveat.
- **Supplement:**
  - cell flow and selection bias;
  - registration R² table and exclusions;
  - coverage maps;
  - feature correlations;
  - all-cell versions of the rhythmic-only maps;
  - cell-level vs embryo-level concordance tables;
  - recording-year analysis;
  - optionally, spectral maps on the paper's exact CNPU set (3.4).
- **Key references:**
  - Torrence & Compo 1998, *Bull. Am. Meteorol. Soc.* 79:61 (wavelet analysis, red-noise significance, COI);
  - Liu, Liang & Weisberg 2007, *J. Atmos. Oceanic Technol.* 24:2093 (scale rectification of wavelet power);
  - PyWavelets (Lee et al. 2019, *J. Open Source Softw.* 4:1237);
  - Ablondi et al. 2020, *JoVE* 156:e60726 (fluorescent calcium imaging and subsequent in situ hybridization in *Xenopus*; cited in the imaging protocol);
  - Paudel et al. 2019, *Int. J. Mol. Sci.* 20:1880 ("Calcium activity dynamics correlate with neuronal phenotype at a single cell level and in a threshold-dependent manner"; used Hurst exponent and Markovian entropy, a precedent for SpectralEntropy);
  - Sudip Paudel's PhD dissertation, William & Mary, 2022: Chapter 4 is the Calcium III PDT study `[CONFIRM exact title and year]`;
  - Marken et al. 2016 (cited by the PDT study for the "lack of periodicity" in developing neural cells) `[CONFIRM full citation]`;
  - Neugornet et al. 2021 (wavelet + ridge-walking peak detection; cited in the dissertation's limitations) `[CONFIRM full citation]`;
  - the PDT manuscript itself (*CB Calcium III Working*, Paudel et al.) `[CONFIRM citation and status]`.

---

## 10. Open questions for the author
1. Staging system (Nieuwkoop–Faber? The hpf windows are known, 2), imaging temperature, and how embryos were devitellinated (for Methods).
2. **Which PDTs to pre-specify** as the headline event reference for step 5. Suggested: `2-111`/`2-112` in all four CW/EW × Th1/Th2 combinations.
3. **How to handle CV_ISI** (Section 12). Options:
   - spikes only, from the PDT event files, for lenient PDTs;
   - request Sudip's per-frame files or his CV outputs;
   - compute it from the v3 traces with a stated detection rule;
   - omit it.

   If it is meant to match Sudip's, note his non-standard formula (13.4).
4. **Files to request from `W:\SudipP\1Calcium III\…`:**
   - PDT group membership (`MethodGroups_RealData.csv`);
   - per-PDT smoothed maps;
   - spike-CV outputs;
   - per-frame event files;
   - region polygons.
5. Whether to add a "paper-CNPU-set" parity version of the spectral maps (3.4).
6. How to handle the abstract's temporal-class claims, which the current pipeline doesn't yet produce (Section 12).
7. Citation details marked `[CONFIRM]` in Section 9.
8. **After step 1 (11-04):**
   - ~~Whether to build the rhythmic-fraction map~~ (built; 6.8). Which version to show: the shared-scale % map, and/or the S13 devitellinated z-map.
   - Whether to pre-specify the **lateral/non-neural enrichment** test (S13 devitellinated) and the **anterior-low** test (S14 devitellinated) for step 4a.
   - Whether to **pre-specify the posterior > anterior rhythmic-fraction test** (seen post hoc in S13 devitellinated) for step 4a.
   - Which poster headline to adopt (Section 7 outcome; Section 8 suggested revision).

*Resolved on 11-03 (no longer open):*
- template orientation (anterior = top);
- PDT column definitions;
- `Counts_0`/`Counts_M`;
- the PDT naming code;
- what "P" and "fpf" mean;
- which embryos the PDT outputs cover (all 179).

---

## 11. File glossary (on the cluster, for the author's reference)
| What | Where |
|---|---|
| Reference table (groups) | `data10/Embryo_Ref_Table.csv` |
| Spectral run v3 | `data10/HPC_spectral_v3/<group>/` (`run_manifest.json`, per-embryo tables, `SpectralLibrary.npz`) |
| Group tables v3 | `CALCIUM/data/group_tables/` |
| Maps v3 | `scr10/CALCIUM/data/mapping_v3/<group>/` and `.../combined/` |
| Analysis plan (steps 0–7, exact commands) | `CALCIUM/PROJECT_PLAN_spectral_next_steps.md` |
| Project notes for Claude Code (most detailed) | `CALCIUM/CLAUDE.md` |
| Analysis scripts | `CALCIUM/scripts/step1_rhythmicity` … `step6_figures`, `analysis_common` |
| Analysis outputs (when run) | `scr10/CALCIUM/data/analysis_v3/<step>/` |
| Red-noise test results | `scr10/CALCIUM/data/analysis_v3/step1_rhythmicity/` (`rhythmicity_all_groups.csv` per cell; `rhythmicity_by_embryo.csv`; `rhythmicity_by_group.csv`; `tables_rhythmic/`; per-embryo `<group>/*_rhythmicity.csv` and `*_rednoise.npz` with surrogate mean and 95% band per cell) |
| Rhythmic-only maps | `scr10/CALCIUM/data/mapping_v3_rhythmic/<group>/` and `.../combined/` |
| Rhythmic-fraction maps (step 1d) | `scr10/CALCIUM/data/analysis_v3/step1_rhythmicity/rhythmic_fraction_maps/` (`RhythmicFraction_PooledSmoothedPct.png` shared scale; `RhythmicFraction_PooledSmoothedZ.png` paper z-map; support QC maps; `Group_Cnpu_RhythmicFraction.csv`; `Embryo_Cnpu_RhythmicFraction.csv`). Script: `CALCIUM/scripts/step1_rhythmicity/rhythmic_fraction_maps.py` |
| **SfN poster figures** (1, 2A, 2B, 2C, 3A, 3B) | `scr10/CALCIUM/data/analysis_v3/step6_figures/sfn_poster/` (PDF + PNG at exact placed size, 300 dpi, transparent, Liberation Sans ≥24 pt; `figure_facts.json` holds every number for the captions). Script: `CALCIUM/scripts/step6_figures/rednoise_bandpower_figures.py` |
| Red-noise audit and P_global rationale | `CALCIUM/logs/calcium_iii_rednoise_audit_2026-11-04.md` |
| Pipeline notes and audits | `CALCIUM/logs/calcium_iii_*.md` (PDT data audit: `calcium_iii_pdt_data_audit_2026-11-03.md`) |
| Methods write-up (author's) | `CALCIUM/Calcium III Wavelet Methods (1).pdf` (not reviewed for this summary) |
| PDT manuscript draft | `data10/CalciumIII_Working.pdf` |
| PDT supplementary text | `data10/CalciumIII_Supplementary_text.pdf` |
| Paudel dissertation | `data10/Paudel_wmgrad_0261E_10696.pdf` (Flowchart 4.1 = PDF p. 106; R code in Appendix 3) |
| **Full PDT outputs (2022)** | `data10/PDT_data/` (per-cell `*_Unique_LastFrameIDs_Perhour_Counts.csv`, event-level `*_ActivityWiseDetails_*.csv`, `AllGroups_Data_for_Plotting.csv`, summary and NNDist tables) |
| Sudip's R scripts | `data10/PDT_data/R_ScriptsSP/` |
| Template boundaries per embryo | `data10/boundary_info/Boundaries/<embryo>-Boundary_SuperEmbryo.csv` |
| Old PDT per-cell summaries (superseded) | `data10/boundary_info/Boundaries/Stage_*/ConvergeSpikeCy5/*.csv` |
| PDT Sox2-vs-activity correlations | `data10/boundary_info/Boundaries/Stage_*/1_R and C.csv` |
| Old pixel matrices (Jan 2021, early S13dev build) | `data10/boundary_info/SuperEmbryoBoundaries/SuperEmbryo_MatrixForm/` |
| Not on the cluster | `W:\SudipP\1Calcium III\…`: group membership, per-PDT maps, CV outputs, per-frame event files, region polygons |

---

## 12. SfN abstract (as submitted) and how it compares with the current data

**Title:** Decoding Calcium Activity in Early Neural Development using Wavelet Analysis. **Category:** A. Development. **Keywords:** Calcium, Wavelets, Xenopus.

> Calcium signaling is a key regulator of neural induction, patterning, and early excitability, yet the fine‑scale temporal structure of calcium activity during neurulation remains poorly understood in contrast to adult neuron activity. To determine how individual neural plate cells encode developmental information at this stage, we developed a continuous wavelet transform (CWT)–based pipeline to complement traditional peak detection methods that extracts multiscale temporal features from high‑speed calcium imaging of Xenopus laevis embryos at early neurula stages, when neuroepithelial identity is being established and cells begin diverging toward distinct fates. The framework quantifies dominant frequency, spectral bandwidth, band‑limited power, and inter‑event interval variability for thousands of single cells.
> Analysis of >2,000 neural plate cells in early neurula embryos using the genetically encoded calcium indicator GCaMP revealed multiple excitability phenotypes rather than a uniform field of activity. Peak event frequencies spanned more than an order of magnitude (~0.01–0.13 Hz), and band‑limited power varied across six orders of magnitude. Inter‑event interval variability (CVISI) ranged from near‑zero to >1.4, distinguishing highly regular oscillators from irregular spikers and burst‑like activity patterns. These features delineated distinct temporal classes: regular narrowband oscillators, broadband irregular spikers, high‑power burst‑prone cells, and near‑silent ROIs. Although visually similar in raw traces, these classes differ in timing statistics that are known to selectively activate downstream pathways such as CaMKII, CREB, calcineurin, and NFAT. The temporal phenotypes likely reflect underlying developmental states—ranging from progenitors with variable responsiveness to inductive cues, to mechanically active cells experiencing morphogenetic strain, to early differentiating neurons exhibiting stabilized oscillatory calcium entry. Because calcium dynamics influence gene expression, cytoskeletal remodeling, and responsiveness to embryonic development signals, these classes provide insight into how neural plate cells may differentially interpret patterning cues and commit to distinct fates. By linking phenotypes to potential developmental roles, wavelet analysis at the single-cell level offers a foundation for connecting excitability states to fate specification, tissue mechanics, and environmental sensitivity across neurulation.

**Where the abstract and the v3 analysis line up, and where they don't.** The abstract preceded the corrected v3 pipeline. Its numbers appear to come from an earlier pilot of about 2,000 cells on the v1 method, which had frequencies above Nyquist and a scale bias.

| Abstract statement | Current (v3) status | Recommendation for the poster |
|---|---|---|
| CWT pipeline complementing peak detection | ✔ Implemented (Section 3); PDT outputs now joinable for all mapped cells (6.6) | Core of the poster; show the PDT comparison (6.7) |
| Dominant frequency, spectral bandwidth, band-limited power | ✔ DominantFreqHz, SpectralSpread (bandwidth), TotalBandPower; plus centroid, entropy, peak fraction | Use as stated |
| ">2,000 neural plate cells" | Now 300,023 analysed cells (191,527 mapped) from 177 embryos | Report the larger n; it strengthens the abstract |
| "Peak event frequencies ~0.01–0.13 Hz" | v3 band is 0.005–0.08 Hz. Frequencies above about 0.08 Hz approach Nyquist for some embryos (measured frame interval up to 4.5 s; Nyquist 0.11 Hz), and v1 included above-Nyquist bins. Many v3 cells peak at the band edges (Section 6.2) | Report the v3 range with the reason for the band limits; avoid the 0.13 Hz figure. **After step 1:** frequency is meaningful only for rhythmic cells. Their excess over red noise peaks at 0.005–0.02 Hz for 96–97% of them (median about 0.010–0.013 Hz). That reflects transient duration (about 26 s wide), not a repeat period (intervals about 6 min, irregular), so report it as "slow transients", not a peak oscillation frequency |
| "Band-limited power varied across six orders of magnitude" | v3: the central 98% of cells spans about 2.8 orders of magnitude (the full min–max span is dominated by near-zero, silent cells) | Report about 3 orders (1st–99th percentile); describe near-silent cells separately |
| "Inter-event interval variability (CV_ISI) from near-zero to >1.4" | **Not computed in v3.** The PDT event files give per-cell event times for spikes only, and ≥3 spikes per cell only for lenient PDTs (6.6 B). Sudip computed a spike CV per PDT, but with CV = sd(gap) / (mean(gap) − min(gap)) in frames, not sd/mean (13.4), and his outputs are not on the cluster | Choose one, stating the definition: (a) spike-only CV_ISI from the PDT event files for a stated lenient PDT; (b) request Sudip's per-frame files or CV outputs; (c) detect events on the v3 ΔF/F traces with a stated rule; (d) present CV_ISI as preliminary or omit it |
| "Temporal classes: regular narrowband oscillators, broadband irregular spikers, high-power burst-prone, near-silent" | **Not yet derived in v3.** Step 1 (done) gives rhythmic vs not (6.8). Rhythmic cells are narrowband, with low entropy and high power, but they are **not** regular oscillators: their transients are irregular (CV_ISI about 0.7; 6.8 caveat 3). The abstract's "regular narrowband oscillators" class is not supported by the current data. Full classes would need clustering on (rhythmicity, PeakPowerFraction/entropy, TotalBandPower, CV_ISI) | A simple, pre-specified rule now looks feasible: rhythmic / non-rhythmic high-power / non-rhythmic low-power ("near-silent"), with class proportions per embryo. Report it as exploratory |
| "High-speed calcium imaging" | Acquisition 0.33 Hz (frame interval about 3–3.6 s) | Use "time-lapse" rather than "high-speed" |
| Links to CaMKII, CREB, calcineurin, NFAT and to developmental states | Interpretive, not tested | Keep in the Discussion box as hypotheses |

Analyses that would close these gaps are not in the step 1–6 scripts yet:
1. a CV_ISI step: from the PDT spike events, Sudip's files if retrieved, or event detection on the v3 segments;
2. a temporal-class step: an exploratory classification on per-cell features, with the class map on the template and class proportions per group with the embryo as the unit.

---

## 13. The PDT study: background for writing

*Sources: the PDT manuscript (Methods §2.1–2.5) and its supplementary text, Sudip Paudel's dissertation (Chapter 4 and Appendix 3), and his R scripts. Where the code and the text disagree, the code is what produced the data, and the difference is noted.*

### 13.1 Imaging, segmentation and tracking
- **Segmentation:** NIS-Elements AR on the memRFP channel, using six binary layers: "Segment Tight Borders" at thresholds 15/30/45/60/75 (inverted and eroded) plus "Spot Detection" (12 µm diameter).
- **Tracking:** objects tracked with a **maximum gap of 3 frames**. New objects (daughter cells, cells drifting in) are kept if present for **≥100 frames**.
- **Merging layers:** objects from different layers are the same cell if their centroids are closer than the larger object's equivalent diameter. The shorter track fills gaps in the longer one.
- **Coverage:** tracking found >85% of imaged cells (versus manual counts).
- **Cell outlines** come from the **last frame**, so only cells present at the end are mapped (a stated limitation).

### 13.2 The 204 peak detection techniques (PDTs)
- **Literature-based PDTs (96).** These come from a survey of 254 in vivo calcium-imaging papers (2019–2021). They combine:
  - **preprocessing:** none, frame-wise minimum (background) subtraction, or Savitzky–Golay filter (order 3, 7 frames);
  - **baseline:** mean of the lowest 10% of intensities, mean intensity, or a 20-frame trailing moving mean;
  - **threshold technique:** SD-based, % of baseline (110% / 200%), or R `findpeaks` (minimum height = mean + 1 SD or + 2 SD);
  - **threshold value:** Th1 (lower) or Th2 (higher).

  Parameters are computed per cell (**CW**, cell-wise) or per embryo (**EW**, embryo-wise) where that makes sense, giving 54 CW + 42 EW.
- **Change-in-magnitude PDTs (108).** These keep negative deviations.
  - **"Local":** ΔF from a local minimum to the next local maximum.
  - **"Modified" ("sandwich"):** a base is a local minimum lower than its neighbouring minima; the peak is the maximum between two bases.
  - The magnitudes are then thresholded by cell-wise SD, embryo-wise SD or % of baseline ("P").
  - 3 preprocessing × 2 (Local/Modified) × 3 baselines × 2 threshold values × 3 threshold types = 108.
- **PDT naming code** (dissertation Flowchart 4.1; five digits):

  | Digit | Step | Values |
  |---|---|---|
  | 1 | Preprocessing | 1 background subtraction · 2 none · 3 Savitzky–Golay filter |
  | 2 | Input signal | 1 Local · 2 Modified; "–" = intensity values (literature PDTs) |
  | 3 | Baseline | 1 lowest 10% · 2 mean · 3 moving mean |
  | 4 | Threshold technique | 1 SD · 2 % of baseline · 3 findpeaks |
  | 5 | Threshold value | 1 = Th1 (1 SD / 110% / mean+1SD) · 2 = Th2 (2 SD / 200% / mean+2SD) |

  - **Note:** digit 1's numbering is not the manuscript's (a)/(b)/(c) listing order (the manuscript lists "none" first).
  - The paper writes literature PDTs as e.g. `2-111.CW` (no preprocessing, lowest-10% baseline, SD, 1 SD, cell-wise) and change-in-magnitude PDTs as e.g. `21121.EW`.
  - In the data files, "P" (the %-of-baseline change-in-magnitude PDTs) is labelled cell-wise.
- **Code details worth stating if methods are compared:**
  - **SD thresholds:** a frame is "active" if F > baseline + Th × SD.
  - **%-of-baseline thresholds:** active if ΔF/F > 10% (Th1) or > 100% (Th2).
  - **Refocus segments:** baselines and SDs are computed separately for each segment between refocus pauses.
  - **Embryo-wide (EW) runs:** each cell's minimum is subtracted first.
  - **What counts as an event:** a run of consecutive active frames. The earlier version's findpeaks was found to be erroneous and rerun; the corrected files carry the suffix "fpf", for "findpeaks fixed".
- **Main finding about methods:** each pipeline step changes peak counts by up to about 10-fold. CW gives more peaks than EW, and Th1 more than Th2.

### 13.3 Event classification
- **Neighbours:** two active cells are neighbours if their centroid distance is less than the sum of half their maximum Feret diameters plus 3 µm (the measured gap between cells is 2.78 ± 0.89 µm).
- **Event types:**
  - **spike** = 1 cell with no active neighbours;
  - **rosette** = 2–4 adjacent cells;
  - **wave** = ≥5 adjacent cells propagating from one or a few initiating "leader cells".
- **Linking and splitting:** events are linked across frames. Where two waves meet, the contact cells go to the wave whose member is nearest.
- **Mean rates per cell per hour:** 3.85 ± 5.13 spikes, 1.43 ± 2.10 rosettes and 0.70 ± 1.24 waves (n = 179 embryos, 201,409 cells).
- **Stage:** average event counts per embryo did not differ between stage 13 and stage 14 (all PDTs pooled).

### 13.4 Inter-event intervals in the PDT code
- Sudip computed, per cell and per PDT, a variability score for **spike** intervals only.
- The gaps are runs of frames between events, excluding gaps that touch the start or end of the track.
- **CV = sd(gap) / (mean(gap) − min(gap)).** This is not the standard sd/mean, and it is in frames.
- The outputs are not on the cluster. Any CV_ISI on the poster must state which definition it uses.

### 13.5 Composite neural plate (CNP) and spatial analysis
- **Registration:**
  1. After imaging, embryos were fixed and processed for whole-mount Sox2 FISH (Cy5), anti-memRFP immunostaining and DAPI.
  2. The calcium last frame was registered to the 20× stained image (≥20 control points, local weighted mean transform).
  3. That image was registered to the 10× image (≥6 points, affine).
  4. The 10× image was registered to the template (≥20 landmarks, local weighted mean).
  5. Each step was checked by ≥3 people.
- **Template:** built from published Sox2 in situ images (Xenbase), with landmark distances normalised to AP length.
- **Coverage:** on average 10–14 embryos contribute to each pixel inside the neural plate (14 ± 4.1, 12 ± 3, 13 ± 3.7, 10 ± 3.1 for S13, S13dev, S14, S14dev).
- **CNPUs and maps:**
  - 9 × 9-px CNPUs, about one cell each (mean cell area 84 ± 47 px).
  - Mean event count per CNPU over all cells from all embryos → 49-neighbour smoothing → z-score.
  - For inclusion, grid and smoothing as actually coded, see 3.4.
- **A priori regions:**
  - AP: forebrain (diencephalon–telencephalon), midbrain, hindbrain, spinal cord;
  - DV: central, mid-lateral, lateral;
  - neural vs non-neural.
  - Ten pairwise K-S + Cohen's d comparisons, reported as **the percentage of PDTs** finding A > B, B > A or neither.
- **PDT groups:**
  - PDTs with similar maps were grouped by complete-linkage hierarchical clustering cut at r = 0.5, so every pair in a group correlates at r ≥ 0.5. Their z-maps were averaged.
  - Only groups with ≥5% of PDTs (≥11) whose |z| > 1 CNPUs covered >1% of the map were kept.
  - **27 groups were retained** (Supplementary Table S3).
- **Spatial pattern test (NNDist):**
  1. Take CNPUs with |z| > 1: FWD = z > 1, REV = z < −1.
  2. Compute each one's mean distance to its 3 nearest neighbours.
  3. Compare with 1,000 shuffles (K-S + d).
  4. The pattern is called clustered, dispersed or random if Bonferroni p < 0.05 and d > 0.2 in ≥95% of shuffles.
- **Simulation (325 synthetic "embryos" on a 3,097-cell grid):**
  - 14 synthetic time series: 12 signal types (amplitude × frequency × regularity), plus noise and a randomly spiking trace.
  - 5 spatial patterns and 5 inside/outside mixes.
  - Results:
    - PDTs were FWD in 43 ± 24, REV in 50 ± 24 and neither in 7 ± 4 conditions.
    - **Th1, CW, and randomly spiking surrounding cells favour REV; Th2 and EW favour FWD.**
    - CW counts don't depend on tissue composition; EW counts do.

### 13.6 Key PDT findings to cite
- **PDT choice drives conclusions:** counts change up to about 10-fold between PDTs, and pairs of PDTs can agree in one simulated condition and reverse in another.
- **No clusters of high- or low-spiking CNPUs inside the neural plate;** cells there are "equally capable of spiking".
- **Clusters outside the plate:** PDT-dependent clusters of spikes and rosettes appear in the pre-placodal ectoderm and near the blastopore / chordoneural hinge.
- **Waves are anterior–posteriorly biased in a PDT-dependent way:** posterior (blastopore, neural ridge) with EW or Th2 PDTs; anterior with CW or Th1 PDTs.
- **Devitellination** (Supplementary Text 5):
  - The change in magnitude per event is higher in devitellinated embryos, significant only at **stage 13** with the Modified approach.
  - With CW PDTs, there is little or no difference in event counts or durations.
  - With EW PDTs, *intact* embryos show more spikes and rosettes, of shorter duration.
  - With EW Th2 PDTs, devitellinated embryos show more waves.
  - The authors' summary: removing the membrane "increases the intracellular calcium levels, which resulted in fewer long-lasting events".
  - **Do not cite this as "devitellination increases event frequency".**
- **Framing for the spectral work:**
  - The literature survey found 31.5% of papers used non-spike-counting methods, such as frequency analysis, alone or with peak counting.
  - The dissertation names Fourier and wavelet methods as quantifying calcium-event frequency "beyond simply detecting peaks". It argues they suit periodic cells, whereas developing neural cells "lack periodicity and predictability" (Marken et al. 2016).
  - Sudip also tried an exploratory, unpublished FFT on the traces. It had no Hz axis, no edge-effect correction and no significance test, so it is not comparable to the CWT here.
