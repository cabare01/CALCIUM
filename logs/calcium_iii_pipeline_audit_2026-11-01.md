# Calcium III full pipeline audit (PyWavelets → group tables → SuperEmbryo mapping)

Audit date: 2026-11-01. This was read-only. All checks ran in memory, and no code, data or mapping outputs were changed. This document is the only file written.

**Scope note:** The four earlier audits (2026-10-26/27) cover the **superseded MATLAB-derived** tables under `data/master_tables/spectral_analysis_tables/`, which are now deleted from the working tree. Their counts (417,650 rows, 11 TotalBandPower values above 1,000, and so on) do not describe the current data. This audit covers the current chain:

| Step | Code | Output |
| --- | --- | --- |
| 1. Drop tracking columns | `scripts/raw_prep/preprocess_stage_*.py` | rewrites `data10/merged_data/<group>/*.csv` in place |
| 2. Spectral features | `scripts/wavelets/pipeline.py` + `01_/02_/03_*.py` (PyWavelets 1.8.0, `cmor1.5-1.0`) | `data10/HPC_spectral/<group>/<embryo>/` |
| 3. Group tables | `scripts/create_hpc_group_master_tables.py` | `data/master_tables/group_tables/*_CW.csv`, `*_EW.csv` |
| 4. Spatial mapping | `scripts/mapping/spatial_spectral_mapping.py` | `scr10/CALCIUM/data/mapping/<group>/` |

## 0. Current population counts

| Group | Raw files | Spectral embryos | Spectral cells (CW) | Cells skipped: no ≥100-pt segment | Cells skipped: <2 obs | Mapped embryos | Mapped rows | Unmatched cells |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| S13 | 46 | 46 | 75,029 | 63,216 | 9,614 | 45 | 49,964 | 24,218 |
| S13_devitellinated | 46 | 46 | 67,480 | 51,181 | 9,122 | 46 | 46,181 (46,172 unique) | 21,308 |
| S14 | **45** | 45 | 78,930 | 66,318 | 11,225 | 45 | 52,431 | 26,499 |
| S14_devitellinated | 41 | 41 | 68,377 | 60,770 | 11,146 | 41 | 41,993 | 26,384 |
| **Total** | 178 | 178 | 289,816 | 241,485 | 41,107 | 177 | 190,569 | 98,409 |

- **Confirmed:** No CW row has missing features. The pipeline drops failing cells instead of writing NaN rows.
- **Confirmed:** All four wavelet jobs started after the last edit to the wavelet scripts on 2026-10-29 at 13:54. `sacct` shows starts at 10-29 13:56, 10-30 14:09, 14:22 and 14:22, so every group used the same code.
- **Confirmed:** No biological embryo (date + embryo number) appears in more than one group, so the groups are independent samples.

## 1. Data integrity and group labels (blocking: resolve before any comparison)

1. **Potential issue: stage or membrane tokens in filenames disagree with the group folder.**

   | Embryo ID | Group it sits in |
   | --- | --- |
   | `190320SP_E6_st_13` | S14 |
   | `190730SP_E3_st_13_devitelined` | S14_devitellinated |
   | `190730SP_E5_st_13_devitelined` | S14_devitellinated |
   | `190327SP_E6_st_13_devitelined001` | S14_devitellinated |
   | `190619SP_E2_st_13_devitelined_albino_Q` | S14_devitellinated |
   | `190411SP_E3_st_14_devitelined` | S13_devitellinated |

   The previous student kept a folder named `boundary_info/Boundaries/wrong stage folder need to cross reference/`. It contains spike files for `190320SP_E6_st_13`, `190730SP_E3/E5_st_13_devitelined` and `190831SP_E4_st_13_devitelined`, among others. The cross-referencing was apparently never finished. **Unresolved:** whether the folder (assigned stage) or the filename (stage at recording) is correct. This is a user/lab decision.
2. **Potential issue: probable swap at 200217 E3.**
   - `200217RN_E3_st_13_devitelined_albino001` is in **intact** S13. It matched no boundary file, so it contributes 0 mapped cells.
   - `200217SP_E3_st_13_albino` (no "devitelined" token) is in **S13_devitellinated**.
   - The boundary `200217SP_E3-Boundary_SuperEmbryo.csv` has 47 repeated `Calcium_ID`s (in the 30xxx and 40xxx series). The first S13_dev run crashed on these (job 583208). The rerun (583210) used an inner-join, which created 9 duplicated mapped rows.

   This pattern fits two recordings that were mislabelled, or merged, across groups.
3. **Confirmed:** S14 has 45 raw files, not 46. See `stage14_preprocess-577775.out` ("Updated 45 CSV files") and the mapping QC warning. The 46th embryo's merged CSV is missing from `data10/merged_data/stage_14`.
4. **Potential issue: wrong boundary association for two S13 embryos.** For each embryo, I fitted an affine transform from the raw-image cell mean position (MeanX/MeanY) to the template centroid (Centroid_1/2) for the joined cells. A correct ID join should fit almost perfectly, since it is the same cells under a registration.
   - The median R² is 0.98–0.99 in every group, and every embryo in the other three groups is at least 0.90.
   - **`200114SP_E1_st_13001`: R² = 0.09.** **`180405SP_E4_st_13`: R² = 0.26.** `180612SP_E2_st_13_Q` is borderline at 0.78.

   These embryos' cells are almost certainly placed at the wrong template positions. The boundary files were matched by stem prefix (`200114SP_E1-…`, `180405SP_E4-…`), so either the boundary belongs to a different recording or its Calcium_IDs come from a different tracking run. Inspect them with `scripts/mapping/check_mapped_boundaries.py`, and exclude them from spatial maps until resolved.
5. **Confirmed: mixed time formats in one embryo.** `190724SP_E2_st_13_devitelined` mixes `hh:mm:ss`, `mm:ss.fff` and plain seconds, sometimes within one cell track: 217 IDs mix clock and numeric values. As a result, 567 cells were skipped: 217 with "Cannot parse time column" and 350 with pandas "Mixed timezones detected". The embryo kept 943 cells. No other embryo had parse errors.
6. **Potential issue: in-place raw rewrite.** `raw_prep/preprocess_stage_*.py` replaces the raw CSVs in place with no backup (docstring line 6, `os.replace` at line 97). **Unresolved:** whether an untouched copy of the original exports exists elsewhere. Confirm before any further in-place step.

## 2. Spectral method (affects every feature value)

Acquisition facts from the raw CSVs: about 1-hour recordings, frames roughly 3.3 s apart, and `Time` stored as `hh:mm:ss` **rounded to whole seconds** (e.g. 3, 6, 9, 13, 16, 19, 23, 26 s). In a 300-cell test on `180209SP_E4_st_14`, the computed `dt_med` was 3 s for 273 cells and 4 s for 5 cells.

1. **Confirmed: the frequency axis goes past Nyquist.** `pipeline.py:54` uses `linspace(0.02, 0.20, 80)`. The true frame rate is about 0.30 Hz, so Nyquist is about 0.14–0.15 Hz. Bins above that contain no real signal, only interpolation and wavelet leakage. They hold about 7% of total power (median, 300-cell test), and they enter TotalBandPower, entropy, spread and centroid.
2. **Confirmed: dominant frequency is pinned to the lowest bin.**

   | Group | Cells with DominantFreqHz = 0.02 Hz |
   | --- | ---: |
   | S13 | 58.8% |
   | S13_devitellinated | 71.1% |
   | S14 | 52.8% |
   | S14_devitellinated | 67.0% |

   The median is 0.02 Hz in every group, and no cell has a dominant frequency at or above 0.15 Hz. `02_spectral_features.py:34` averages `|coef|²` without scale normalization. PyWavelets coefficient power grows with scale, which biases every spectrum toward low frequency. In the 300-cell test, the share of cells at 0.02 Hz was:

   | Variant | Share at 0.02 Hz |
   | --- | ---: |
   | As run | 65.5% |
   | Power rectified by scale | 15.1% |
   | 5-min rolling-median detrend only | 50.0% |

   The pile-up is therefore mainly the normalization artifact, plus some slow drift or bleaching. As computed, DominantFreqHz is largely uninformative, and centroid, spread, entropy and PeakPowerFraction share the same bias.
3. **Confirmed: segments are spliced before the CWT.** As in the MATLAB method, retained segments are concatenated (`01_preprocess_data.py:181`) and transformed in a single CWT, with no cone-of-influence masking. In the test, 64% of cells had at least one internal gap. The baseline is one global 10th percentile with an eps floor (`:182`), with no bleaching correction.
4. **Potential issue: the gap rule fragments traces.** A gap is any interval over 2×dt_med (`:160`), which is more than 6 s. Timestamps are rounded and frames are about 3.3 s apart, so one dropped frame (about 6.6 s, recorded as 6 or 7 s) can split the trace. Segments under 100 points (about 5 min) are then discarded (`:170`). Even among tracks with ≥300 frames, 22 of 300 failed this rule. With an overall skip of about 241k cells, the retained population is biased toward long, stable tracks.
5. **Confirmed: provenance is thin.** Outputs record PyWavelets version and wavelet only in the job `.out` file. `RunID` is a per-cell wall-clock timestamp. No git SHA, frequency axis or thresholds are stored with the outputs. A per-cell PNG at 200 dpi (`pipeline.py:92`) accounts for most of the 32 GB in `HPC_spectral`.
6. **Confirmed: the old extreme values are gone.** The maximum TotalBandPower in the new tables is 25.6 (S14). The 10⁹-scale values from the old MATLAB audit no longer appear, and no group has a value above 1,000. TotalBandPower remains an uncalibrated coefficient-power sum.

## 3. Mapping and group maps

1. **Confirmed: the cell join is sound apart from §1.2 and §1.4.** 31–39% of spectral cells lack a boundary centroid. They are mainly short track fragments: unmatched cells have a median duration of about 20 min versus about 59 min for matched cells. Over 92% of cells with Duration ≥ 55 min match, compared with 26–30% of cells under 30 min. This fits boundaries segmented on the last frame. Matched and unmatched cells have similar feature medians.
2. **Confirmed: CNPU coverage is sparse.**
   - No CNPU reaches 50% embryo coverage; the maximum is 43.5% / 37.0% / 48.9% / 43.9% for S13 / S13_dev / S14 / S14_dev.
   - The median is 3–5 contributing embryos per observed CNPU.
   - 1,538 / 1,995 / 2,282 / 2,644 CNPUs rest on a single embryo.
   - The median is 1 cell per embryo per CNPU.

   Group medians at 9-px resolution are therefore noisy. The 7×7 smoother (`spatial_spectral_mapping.py:291-316`) averages CNPUs with equal weight, regardless of how many embryos support each one.
3. **Confirmed: colour scales are not comparable across groups.** Each sbatch runs a single `--group`, so the "shared" colour limits (`:383-397`) cover that group only. `SmoothedZ` (`:316`) is z-scored within each group, which removes level differences between groups. Side-by-side group PNGs therefore cannot be compared by colour.
4. **Confirmed: some cells fall outside the grid.** 23 / 426 / 830 / 517 mapped cells lie outside the legacy grid. They stay in `All_MappedCells_with_CNPU.csv` but are not part of the grid summaries.

## 4. Design confounds for statistics

1. **Potential issue: a strong batch effect by recording year.** Embryo-median TotalBandPower, taken as the median across embryos:

   | Group | 2018 | 2019 | 2020 |
   | --- | ---: | ---: | ---: |
   | S13 | 0.0147 | 0.0229 | 0.0393 |
   | S13_devitellinated | 0.0255 | 0.0449 | 0.0501 |
   | S14 | 0.0134 | 0.0219 | 0.0346 |
   | S14_devitellinated | 0.0295 | 0.0476 | 0.0449 |

   The groups are unbalanced by year. S13 is 32/7/7 embryos across 2018/2019/2020, while S13_dev is 10/14/22 and S14_dev is 5/14/22. The roughly 2× higher TotalBandPower in the devitellinated groups is therefore partly confounded with year (microscope, dye or settings). Albino embryos number 5/6/5/8 per group, and their effect is inconsistent.
2. **Confirmed: mixed models are feasible.** There are 67 recording dates and 33 span more than one group; 15 include both membrane statuses and 32 include both stages. A random effect for date is estimable.
3. **Confirmed: integration data exist.**
   - There are 87 `ConvergeSpikeCy5` per-cell event + Cy5 files under `Boundaries/Stage_*/ConvergeSpikeCy5/` and the "wrong stage" folder. 82 match mapped embryos, covering 94,345 of 190,569 mapped cells. Five files do not match because of naming variants (e.g. `191030SPRN_…`, `191020RNSP_…`).
   - Every boundary CSV carries a per-cell `CY5_Intensity` (Sox2), fully populated in the file checked. Spectral-versus-Sox2 analysis is therefore possible for all mapped cells, not just the spike subset.
   - `Stage_*/1_R and C.csv` holds the previous student's per-embryo Cy5-versus-spike-feature correlations, which can serve as a replication reference.

## 5. Prioritized next analysis steps

**Step 0: resolve inputs (blocking).**
1. Decide group membership for the §1.1 embryos and the §1.2 200217 E3 pair.
2. Locate the missing 46th S14 embryo, and confirm an original raw backup exists.
3. Inspect the boundary files for `200114SP_E1` and `180405SP_E4` with `check_mapped_boundaries.py`; re-source them or exclude those embryos.
4. Make time parsing per-value tolerant of mixed formats for `190724SP_E2`.

**Step 1: revise and validate the spectral method, then rerun.**
1. Cap the axis below Nyquist, at about 0.12–0.14 Hz.
2. Rectify power by scale (divide by scale), or use a normalized CWT.
3. Run the CWT per segment and mask the cone of influence, instead of concatenating.
4. Tolerate single dropped frames by interpolating short gaps (≤2 missing frames, or ≤10 s).
5. Use a rolling or bleaching-corrected baseline with a sane F0 floor.
6. Make the per-cell PNGs optional.
7. Write a run manifest: git SHA, wavelet, axis, thresholds and package versions.
8. Validate on synthetic sinusoids at known frequencies (e.g. 0.03, 0.06, 0.10 Hz) with injected gaps and drift. Then run one pilot embryo and compare with the current outputs.
9. Rerun the four groups and rebuild the CW/EW tables with `create_hpc_group_master_tables.py`.

**Step 2: re-map with fixes.**
1. Compute colour limits across all groups (one combined job, or a second pass), and report z-scores alongside absolute values.
2. Mask CNPUs below a minimum coverage (e.g. at least 5 embryos), and weight the smoothing by embryo count.
3. Add region-level summaries: inside versus outside the stage-matched Sox2 neural-plate polygon, and signed distance to the boundary. At the current coverage these are much more robust than 9-px CNPUs.

**Step 3: statistics with the embryo as the unit.**
1. Fit a mixed model: feature ~ Stage × Membrane + albino + (1 | date) + (1 | embryo) on cells. Alternatively, fit embryo medians with (1 | date).
2. Report a sensitivity analysis by year, and repeat with and without the §1 problem embryos.
3. Keep explicit denominators: spectral cells, mapped cells, and cells with event data.

**Step 4: integration.**
1. Correlate cell-level spectral features with `CY5_Intensity` (Sox2) for all mapped cells.
2. Run event–spectral correlations for the 82 embryos with spike data, reusing the `correlation` / `corr_summary` workflows in `scripts/calcium_iii_pywavelets_workflows.py`.
3. Compare against the previous student's `1_R and C.csv`.

## 6. Decisions needed from the user
1. The authoritative stage and membrane label for the §1.1 and §1.2 embryos.
2. The upper frequency limit; I recommend 0.12–0.14 Hz.
3. The gap tolerance, and whether to keep the ≥100-point (about 5 min) minimum segment length.
4. Whether to keep the legacy 9-px CNPU maps as the primary spatial readout, or move to region-level (Sox2 inside/outside) summaries.

## Evidence inspected
- Scripts:
  - `wavelets/pipeline.py`
  - `wavelets/01_preprocess_data.py`, `wavelets/02_spectral_features.py`, `wavelets/03_build_results.py`
  - all `run_pipeline_stage_*.sbatch`
  - `raw_prep/*`
  - `create_hpc_group_master_tables.py`
  - `mapping/spatial_spectral_mapping.py`, `mapping/check_mapped_boundaries.py`
  - `mapping/*.sbatch`
- Logs:
  - `st*_wavelets-*.out/.err` (skip reasons tallied from the `.err` files)
  - `stage*_preprocess-*.out`
  - `spatial_*-*.out/.err`, including the failed jobs 583190 and 583208
  - `sacct` job timings
- Data:
  - all 8 group tables
  - all 4 `All_MappedCells_with_CNPU.csv` and `Group_Cnpu_Spectral_Summary.csv`
  - raw `merged_data` headers and times, plus a 300-cell in-memory reprocessing of `180209SP_E4_st_14` using the production preprocessing and spectral functions
  - the `boundary_info` tree, including the `ConvergeSpikeCy5` files and the Sox2 templates
