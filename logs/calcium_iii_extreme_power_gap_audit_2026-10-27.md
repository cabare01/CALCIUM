# Calcium III population, TotalBandPower, and CWT gap audit

Audit date: 2026-10-27. This was read-only. No saved data or code were changed; no intermediate outputs were saved; no group-comparison hypothesis tests or exclusion decisions were made. The only output created is this requested audit document.

## 1. Resolved cell-count explanation

- **Confirmed:** The current cell spectral master contains 417,650 rows across 179 embryos. These are all cell records assembled into the master, whether or not all six feature columns are valid.
- **Confirmed:** 415,186 rows have valid values for all six spectral features. The other 2,464 rows are missing all six features. This is a feature-valid population, not the full row count.
- **Confirmed:** The Wavelet Methods PDF's spectral preprocessor retains a continuous segment only when its resampled grid has at least 10 observations (`preprocessCellTrace`, PDF pp. 6–8). The documented event/spectral integration uses an inner join by `CellID`; the correlation workflow reports event, spectral, and merged row counts (PDF pp. 13 and 27). An inner join necessarily retains only shared keys.
- **Confirmed from user-provided context, not independently recounted:** 199,814 is the matched event-detection/spectral population used for integrated spatial analysis. It describes a different analysis population from both the full spectral master and the all-six-features-valid subset.
- **Potential issue:** The stated previous-student eligibility threshold of at least 100 time points would make the matched population selected by trace length as well as event-data availability. The local PDF documents the spectral side's 10-resampled-observation segment rule and the inner-join mechanism, but the previous student's event-detection code and the saved matched-cell table/count were not present in the inspected CALCIUM or `data10/cabare01` trees. Thus the exact 100-point rule and 199,814 count could not be independently verified here.

The three counts should therefore remain separate in future work: 417,650 full spectral table rows; 415,186 complete six-feature rows; 199,814 event/spectral matched cells as reported for integrated spatial analyses. The first two were confirmed from the master audit. The third and its eligibility rule remain user-reported pending the event-pipeline code and merge record. Spectral-only summaries refer to the spectral population; integrated event or composite-template analyses refer to the matched population.

## 2. TotalBandPower investigation

### Source records

All 11 finite cell-master values above 1,000 were re-matched by `EmbryoID`, `CellID`, and `InputFile` to the original per-cell spectral-result CSV. The source CSV line number below includes its header as line 1. Each source-row `TotalBandPower` string exactly equals the corresponding master value.

| EmbryoID | CellID | GlobalID | RunID | SourceFile carried in row | Source result CSV (line) | TotalBandPower |
| --- | ---: | --- | --- | --- | --- | ---: |
| 180125SP_E2_st_13 | 1833 | 180125SP_E2_st_13::1833 | 2026-07-21_151503 | 180125SP_E2_st_13 - Not_self_Merged.csv | `stage_13/180125SP_E2_st_13_spectral_results.csv`: 1882 | 3,017,943,905.64631 |
| 180125SP_E2_st_13 | 1854 | 180125SP_E2_st_13::1854 | 2026-07-21_151514 | 180125SP_E2_st_13 - Not_self_Merged.csv | `stage_13/180125SP_E2_st_13_spectral_results.csv`: 1894 | 45,399,825.841592 |
| 191020SP_E2_st_13 | 1869 | 191020SP_E2_st_13::1869 | 2026-07-23_172107 | 191020SP_E2_st_13 - Not_self_Merged.csv | `stage_13/191020SP_E2_st_13_spectral_results.csv`: 1731 | 12,955,379,820.3898 |
| 190411SP_E1_st_13_devitelined | 1507 | 190411SP_E1_st_13_devitelined::1507 | 2026-07-25_041345 | 190411SP_E1_st_13_devitelined - Not_self_Merged.csv | `stage_13_devitellinated/190411SP_E1_st_13_devitelined_spectral_results.csv`: 1319 | 6,681,046,118.84971 |
| 201007SP_E6_st_13_devitelined | 674 | 201007SP_E6_st_13_devitelined::674 | 2026-07-27_120002 | 201007SP_E6_st_13_devitelined - Not_self_Merged.csv | `stage_13_devitellinated/201007SP_E6_st_13_devitelined_spectral_results.csv`: 1097 | 249,122,157.129723 |
| 191020SP_E3_st_14 | 1284 | 191020SP_E3_st_14::1284 | 2026-07-27_232250 | 191020SP_E3_st_14 - Not_self_Merged.csv | `stage_14/191020SP_E3_st_14_spectral_results.csv`: 1272 | 361,484,135.832097 |
| 201011SP_E5_st_14 | 1050 | 201011SP_E5_st_14::1050 | 2026-07-28_185416 | 201011SP_E5_st_14 - Not_self_Merged.csv | `stage_14/201011SP_E5_st_14_spectral_results.csv`: 1583 | 1,929,635,860.60135 |
| 201011SP_E2_st_14_devitelined | 2158 | 201011SP_E2_st_14_devitelined::2158 | 2026-08-04_113503 | 201011SP_E2_st_14_devitelined - Not_self_Merged.csv | `stage_14_devitellinated/201011SP_E2_st_14_devitelined_spectral_results.csv`: 1887 | 623,675,270.691646 |
| 201014SP_E1_st_14_devitelined | 1456 | 201014SP_E1_st_14_devitelined::1456 | 2026-08-04_121049 | 201014SP_E1_st_14_devitelined - Not_self_Merged.csv | `stage_14_devitellinated/201014SP_E1_st_14_devitelined_spectral_results.csv`: 1292 | 111,733,775.963161 |
| 201014SP_E1_st_14_devitelined | 2318 | 201014SP_E1_st_14_devitelined::2318 | 2026-08-04_121522 | 201014SP_E1_st_14_devitelined - Not_self_Merged.csv | `stage_14_devitellinated/201014SP_E1_st_14_devitelined_spectral_results.csv`: 1536 | 6,035,022,834.77704 |
| 201024SP_E1_st_14_devitelined_albino | 2448 | 201024SP_E1_st_14_devitelined_albino::2448 | 2026-08-05_024308 | 201024SP_E1_st_14_devitelined_albino - Not_self_Merged.csv | `stage_14_devitellinated/201024SP_E1_st_14_devitelined_albino_spectral_results.csv`: 2552 | 412,340,733.315302 |

**Confirmed:** These values pre-exist the master-table assembly. In `scripts/create_spectral_analysis_tables.py`, `INPUT_COLUMNS` includes `TotalBandPower` (lines 42–46); `build_embryo_records` parses input feature values for summaries (lines 125–164); and `write_cell_master` copies each original input row value into the output (lines 167–186). There is no CWT calculation or feature scaling in this builder.

### Formula and limits of the evidence

**Confirmed from the attached PDF:** `preprocessCellTrace` calculates `F0 = prctile(traceUniform,10)`, floors it with `max(F0,eps)`, and calculates `dFF = (traceUniform-F0)/F0` (PDF pp. 7–8). `extractSpectralFeatures` computes MATLAB analytic-Morlet CWT coefficients, squares their magnitudes, averages power over time, linearly interpolates to 80 frequencies from 0.02–0.20 Hz with zero fill, clamps power at zero, and sets `TotalBandPower = sum(meanPower)` for positive totals (PDF p. 10). The normalized spectrum `P` is used for entropy; it is not used to normalize TotalBandPower.

**Confirmed implication:** This is a sum of CWT coefficient power over a fixed set of bins, not a frequency-bin-width integral. The PDF supplies no physical calibration or measurement unit for that sum. Values depend on the signal normalization and CWT implementation/settings; magnitude alone cannot establish an error.

**Potential issue:** If an input trace has a zero, negative, or very small 10th-percentile baseline, the epsilon floor can yield a very large ΔF/F. Squared CWT coefficients and their sum can then be very large. This is a plausible pathway, not a finding about these specific cells.

**Unresolved for all 11 records:** The targeted search under `/sciclone/home/cabare01` and `/sciclone/data10/cabare01` found the nine embryo spectral-result CSVs, but no matching raw `Not_self_Merged.csv` traces, no Excel trace files, and no saved `SpectralLibrary.mat` files. The source-result CSVs hold scalar feature columns, not the 80-bin mean-power vectors. Therefore the pre-epsilon P10 baseline, epsilon-floor status, ΔF/F minimum/maximum/discontinuities, 80 mean-power values, and their per-bin contributions cannot be recovered from the inspected files. The 11 source rows prove that the values were present in the spectral-result outputs, but do not by themselves validate them against fluorescence traces.

**Unresolved on processing versions:** The `RunID`s are timestamps. PDF `buildResultRow` constructs a timestamped `RunID` (p. 12), but no run manifest, executable MATLAB source version, parameter file, or relevant job log was found. The date strings do not identify a code commit or establish that settings changed. The PDF documents one method version; whether every listed output was produced with precisely that version remains unverified.

## 3. Acquisition-gap investigation

### Documented execution sequence

**Confirmed from PDF `preprocessCellTrace` (pp. 5–8):** After removing rows with NaN in time, fluorescence, X, or Y, it calculates `dt = diff(time)`, `dt_med = median(dt)`, and marks `gapIdx = find(dt > 2*dt_med)`. It later validates finite times and positive `dt_med`, calculates `Fs = 1/dt_med`, and makes segment starts/ends from those original gap indices. Each segment is uniqued in time, assigned `expectedN = round(segment_duration/dt_med)+1` (minimum 2), and resampled with `linspace` and linear interpolation. A segment is retained only if its resampled grid has at least 10 points. Fluorescence and X/Y arrays for retained segments are appended to `traceUniformAll`, `xUniformAll`, and `yUniformAll`. The combined fluorescence vector is normalized to ΔF/F after concatenation.

**Confirmed:** The driver calls `preprocessCellTrace` and passes its returned `data` into `extractSpectralFeatures` (PDF pp. 2–3). The feature function calls `cwt(dFF, Fs, 'amor')` once on the entire `data.dFF` vector (PDF p. 10). It receives no segment boundaries or timestamps.

**Confirmed answer to the adjacency question:** Samples on opposite sides of every retained gap are adjacent in the concatenated `dFF` array. The CWT is applied to that concatenated array with one constant `Fs`; it does not transform segments separately. Thus the CWT treats each join as a one-sample-step transition, despite the real elapsed-time gap.

### Potential effects and methodological concerns

- **Confirmed behavior / Potential impact:** The segment-wise interpolation avoids filling fluorescence values across the long gap. However, concatenation does not preserve the gap for the later transform. A value change across a join can act like a sharp transient and plausibly add broadband or high-frequency wavelet power near the join. The size and direction of changes to the six per-cell features are unresolved without the affected traces and segment boundaries.
- **Potential issue:** `dt_med` is computed across all adjacent retained input times, including long gaps. The median is robust when ordinary intervals are the majority; if gaps/irregular intervals are numerous, the representative interval and `Fs` could be shifted, and the `>2*dt_med` detector could classify intervals differently. Per-cell time vectors are needed to assess whether this occurs.
- **Potential issue:** Segment grid length is rounded from duration divided by `dt_med`, then points are spaced by `linspace`. The resulting segment interval can differ slightly from `dt_med`; nevertheless one global `Fs=1/dt_med` is supplied to CWT. The likely impact depends on segment durations and residual sampling jitter and cannot be quantified here.
- **Confirmed behavior / Potential concern:** `data.time` keeps the original segment time coordinates when arrays are appended, so `Duration_min = (time(end)-time(1))/60` includes elapsed gaps, while the CWT sees only concatenated samples at nominal `Fs`. Duration is descriptive and is not used to normalize the documented TotalBandPower, but the two notions of signal duration differ when gaps exist.
- **Potential issue:** The PDF computes `dt`, its median, and `gapIdx` before its later finite/positive/monotonic checks. It checks finite values and positive median interval, but the shown code does not reject decreasing time values. Unsorted timestamps could therefore make gap segmentation unreliable. The missing `convertTimeToSeconds` implementation also leaves timestamp units/conversion unverified.
- **Confirmed scope of the stated goal:** The code does avoid *interpolation across* flagged gaps. It does not preserve discontinuities for the subsequent CWT, so it does not ensure that the CWT treats recording segments as temporally separate.

The pre-existing translated scripts in this project preserve the same sequence: MATLAB port `calcium_iii_hpc_workflows.m:128–151` concatenates segments and calls CWT once at lines 151–153; Python port `calcium_iii_pywavelets_workflows.py:124–165` concatenates arrays and calls PyWavelets once at lines 163–165. These are later ports, not proof of which exact source-code revision generated the archived spectral-result CSVs. The table builder itself only assembles existing features (`create_spectral_analysis_tables.py:167–186`).

No limited in-memory signal check was run because the relevant raw signals are unavailable. With raw traces and the executed code/settings, a targeted validation could compare the existing concatenated-signal calculation with per-segment CWT calculations for affected cells and representative gapped cells, reporting feature differences in memory only. That comparison is necessary to estimate impact; the code inspection alone confirms the join behavior but not its magnitude in current results.

## 4. Prioritized next steps

1. **Before interpreting group-level spectral comparisons:** Locate the raw trace files and exact MATLAB source/run settings for the 11 high-power cells and cells with acquisition gaps. Recover or regenerate in memory their pre-floor P10 baselines, ΔF/F range, detected gap indices, segment lengths, and 80-bin spectra. Use that evidence to judge whether observed values follow the documented calculation; do not classify by magnitude alone.
2. **Before interpreting group-level spectral comparisons:** Quantify the impact of gap joins on current feature estimates using paired in-memory calculations with the original continuous segments kept separate. Report how many cells have gaps and feature-level differences; save no intermediate outputs during this audit.
3. **Before analyses combining event and spectral features:** Locate the previous student's event filter and the saved merge/mapping records to verify the at-least-100-point rule and the 199,814 matched-cell count. Maintain separate labels and denominators for the full spectral, all-six-valid, and event/spectral matched populations.
4. **Document separately:** Establish the time conversion helper and run-version provenance, define calibration/units for summed CWT power if available, and document whether `TotalBandPower` should be interpreted as a discrete 80-bin coefficient-power sum rather than a physically calibrated integral.

### Evidence inspected

- `/sciclone/home/cabare01/CALCIUM/Calcium III Wavelet Methods (1).pdf`, especially `preprocessCellTrace` (pp. 5–8), `extractSpectralFeatures` (p. 10), `buildResultRow` (p. 12), and the event/spectral inner-join description/code (pp. 13 and 27).
- `/sciclone/home/cabare01/CALCIUM/scripts/create_spectral_analysis_tables.py`, especially `INPUT_COLUMNS`, `build_embryo_records`, and `write_cell_master`.
- `/sciclone/home/cabare01/CALCIUM/scripts/calcium_iii_hpc_workflows.m:103–160` and `/sciclone/home/cabare01/CALCIUM/scripts/calcium_iii_pywavelets_workflows.py:100–181` for the already-created ports' corresponding execution sequence.
- Both existing master CSVs and the nine relevant source spectral-result CSVs under `/sciclone/data10/cabare01/embryo_spectral_results/`.
- Targeted read-only file searches under `/sciclone/home/cabare01` and `/sciclone/data10/cabare01` for the nine affected embryo stems and raw trace/library file types.

