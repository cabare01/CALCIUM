# Calcium III spectral tables read-only audit

Audit date: 2026-10-26; methods supplement added 2026-10-27. Calculations were performed in memory. The two tables are under `data/master_tables/spectral_analysis_tables/`, one directory below the path named in the request.

## 1. Overall audit summary

- **Confirmed:** Cell master: 417,650 rows × 20 columns. Embryo summary: 179 rows × 38 columns. The four reported group cell totals and embryo counts match exactly.
- **Confirmed:** Both tables have the same 179 embryos. No duplicate/missing embryo IDs, duplicate `(EmbryoID, CellID)` pairs, or duplicate `GlobalID` values.
- **Confirmed:** 2,464 rows are missing all six spectral features; none have only partial feature missingness. Each feature has 415,186 finite observations.
- **Confirmed:** All 5,370 embryo-feature statistics (valid n, mean, sample SD, median, IQR) in the embryo summary agree with independent calculations from the cell table within numeric serialization tolerance.
- **Confirmed:** There are 11 finite `TotalBandPower` values above 1,000, with maximum 12,955,379,820.3898. All 11 occur identically in the corresponding original per-cell spectral-result CSV rows. The table-assembly script did not create these high values.
- **Confirmed:** The attached Calcium III Wavelet Methods PDF defines the documented calculation: `TotalBandPower` is the sum of time-averaged squared CWT coefficient magnitudes after interpolation to 80 bins from 0.02–0.20 Hz. The table-building script copies the result and does not calculate it.
- **Potential issue:** The documented preprocessing floors the 10th-percentile fluorescence baseline at machine epsilon before computing ΔF/F. If a source trace has a zero/negative or extremely small baseline, this can make ΔF/F and subsequent CWT power very large.
- **Unresolved:** The corresponding raw traces and upstream execution/version records were not available in the inspected source-results tree, so the baseline mechanism cannot be confirmed for the 11 records, and the numerical values cannot yet be judged against their raw signals. The PDF does not state a physically calibrated unit for the summed CWT power.

## 2. Cell-level master table audit

File: `data/master_tables/spectral_analysis_tables/cell_level_spectral_master.csv`.

**Confirmed:** 417,650 rows and 20 columns. Column meanings based on `scripts/create_spectral_analysis_tables.py`, especially `CELL_COLUMNS`, `INPUT_COLUMNS`, and `write_cell_master`:

| Column | Apparent meaning |
| --- | --- |
| Group | Experimental group label |
| Stage | Developmental stage (13 or 14) |
| MembraneStatus | intact or devitellinated |
| EmbryoID | Identifier derived from input result filename |
| InputFile | Per-embryo spectral-result CSV filename |
| GlobalID | Constructed as `EmbryoID::CellID` |
| CellID | Cell ID from input result CSV |
| RunID | Source processing run identifier |
| SourceFile | Source-record filename carried from input row |
| DominantFreqHz | Dominant frequency (Hz by name) |
| SpectralCentroid | Spectral centroid feature |
| SpectralSpread | Spectral spread feature |
| SpectralEntropy | Spectral entropy feature |
| TotalBandPower | Sum of time-averaged squared CWT coefficients on the common frequency bins; physical calibration/units unspecified |
| PeakPowerFraction | Peak power fraction feature |
| Duration_min | Duration in minutes by name |
| MeanX, MeanY | Mean spatial coordinates; units unspecified |
| EmbryoCellCount | Number of input cell rows for this embryo |
| GroupCellCount | Number of cell rows in the group |

All five requested identifier/provenance fields are present and complete. `EmbryoID` has 179 unique values. `CellID` has 19,748 unique values across rows and is reused across embryos, so it is not globally unique on its own. `(EmbryoID, CellID)` has zero duplicates. `GlobalID` has 417,650 unique values and zero duplicates. `RunID` repeats and has 346,947 distinct values; `SourceFile` has 184 distinct values. These are provenance fields, not row-unique keys.

| Group | Embryos | Cell rows | Intended embryos | Match? |
| --- | ---: | ---: | ---: | --- |
| S13 | 46 | 109,831 | 46 | Yes |
| S13_devitellinated | 46 | 90,917 | 46 | Yes |
| S14 | 46 | 116,937 | 46 | Yes |
| S14_devitellinated | 41 | 99,965 | 41 | Yes |
| **Total** | **179** | **417,650** | **179** | **Yes** |

`EmbryoCellCount` agrees with actual row counts for every embryo. `GroupCellCount` is constant per group and matches the table above. No listed cell-level metadata or identifier field is blank.

| Feature | Missing rows | Finite observations |
| --- | ---: | ---: |
| DominantFreqHz | 2,464 | 415,186 |
| SpectralCentroid | 2,464 | 415,186 |
| SpectralEntropy | 2,464 | 415,186 |
| SpectralSpread | 2,464 | 415,186 |
| TotalBandPower | 2,464 | 415,186 |
| PeakPowerFraction | 2,464 | 415,186 |

**Confirmed:** The same 2,464 rows are missing all six features; zero rows are partially missing. No infinities or negative feature values were found. Textual NaN values are treated as missing.

`TotalBandPower` distribution excludes missing values (n=415,186); percentiles use linear interpolation:

| Statistic | Value |
| --- | ---: |
| Q1 (25th percentile) | 0.008476012740617721 |
| Median | 0.01637671167077325 |
| Q3 (75th percentile) | 0.03295146733448433 |
| 90th percentile | 0.07385914324665745 |
| 95th percentile | 0.1419918795387855 |
| 99th percentile | 0.45708837296201493 |
| Maximum | 12,955,379,820.3898 |
| Values > 1 | 1,032 |
| Values > 1,000 | 11 |

## 3. Embryo-level summary audit

File: `data/master_tables/spectral_analysis_tables/embryo_level_spectral_summary.csv`.

**Confirmed:** 179 rows, 38 columns, 179 unique nonmissing `EmbryoID`s, no duplicate IDs, and no missing IDs. Group counts are 46/46/46/41. Group labels, stages, and membrane-status labels match the design and cell master. Metadata columns (`Group`, `Stage`, `MembraneStatus`, `EmbryoID`, `InputFile`, `EmbryoCellCount`, `GroupEmbryoCount`, `GroupCellCount`) have no blank values or inconsistencies.

Each feature has exactly these statistic columns: `<feature>_n_valid`, `<feature>_mean`, `<feature>_sd`, `<feature>_median`, and `<feature>_iqr`. Feature prefixes are `DominantFreqHz`, `SpectralCentroid`, `SpectralEntropy`, `SpectralSpread`, `TotalBandPower`, and `PeakPowerFraction`.

`EmbryoCellCount` matches the cell master for all embryos. Group embryo counts are 46, 46, 46, and 41; group cell counts are 109,831, 90,917, 116,937, and 99,965. They agree with the cell rows and summed embryo counts. All 179 `InputFile` values resolve in the corresponding group folders under `/sciclone/home/cabare01/data10/embryo_spectral_results/`.

## 4. Cross-table consistency checks

**Confirmed:** The tables contain identical sets of 179 `EmbryoID`s. Group, stage, and membrane status agree for every embryo. For every embryo-feature pair, I treated blank/NA/NaN/null and non-finite values as invalid, then calculated valid n, arithmetic mean, sample SD (n−1 denominator), median, and linearly interpolated Q3−Q1. These rules match `numeric_or_none` and `feature_summary` in the builder for the observed data (no infinities were present).

All 5,370 comparisons (179 × 6 × 5) matched within `max(1e-9, 1e-10 × |computed value|)`. No mismatched embryo or statistic was found; thus no observed difference needs an explanation based on rounding or missing-value handling. Every `*_n_valid` matches valid per-cell observations.

## 5. Cell-count discrepancy investigation

**Confirmed:** The cell master row total is 417,650, exactly the sum of the group counts and the previously reported totals. Embryo-level `EmbryoCellCount` values and group counts agree with the cell master.

The 2,464 all-feature-missing rows remain in cell and group row totals but are excluded from each feature's `*_n_valid` and descriptive statistics. Thus 417,650 cell rows and 415,186 valid observations per feature use different counting definitions; the difference is accounted for by those missing rows.

## 6. TotalBandPower anomaly investigation

**Confirmed:** Each finite value above 1,000 was found unchanged in the matching `CellID` row of its `InputFile` spectral-result CSV. The source row also has matching `RunID` and `SourceFile`. These records do not coincide with the all-six-feature-missing rows.

| EmbryoID | CellID | GlobalID | RunID | SourceFile | TotalBandPower |
| --- | ---: | --- | --- | --- | ---: |
| 180125SP_E2_st_13 | 1833 | 180125SP_E2_st_13::1833 | 2026-07-21_151503 | 180125SP_E2_st_13 - Not_self_Merged.csv | 3,017,943,905.64631 |
| 180125SP_E2_st_13 | 1854 | 180125SP_E2_st_13::1854 | 2026-07-21_151514 | 180125SP_E2_st_13 - Not_self_Merged.csv | 45,399,825.841592 |
| 191020SP_E2_st_13 | 1869 | 191020SP_E2_st_13::1869 | 2026-07-23_172107 | 191020SP_E2_st_13 - Not_self_Merged.csv | 12,955,379,820.3898 |
| 190411SP_E1_st_13_devitelined | 1507 | 190411SP_E1_st_13_devitelined::1507 | 2026-07-25_041345 | 190411SP_E1_st_13_devitelined - Not_self_Merged.csv | 6,681,046,118.84971 |
| 201007SP_E6_st_13_devitelined | 674 | 201007SP_E6_st_13_devitelined::674 | 2026-07-27_120002 | 201007SP_E6_st_13_devitelined - Not_self_Merged.csv | 249,122,157.129723 |
| 191020SP_E3_st_14 | 1284 | 191020SP_E3_st_14::1284 | 2026-07-27_232250 | 191020SP_E3_st_14 - Not_self_Merged.csv | 361,484,135.832097 |
| 201011SP_E5_st_14 | 1050 | 201011SP_E5_st_14::1050 | 2026-07-28_185416 | 201011SP_E5_st_14 - Not_self_Merged.csv | 1,929,635,860.60135 |
| 201011SP_E2_st_14_devitelined | 2158 | 201011SP_E2_st_14_devitelined::2158 | 2026-08-04_113503 | 201011SP_E2_st_14_devitelined - Not_self_Merged.csv | 623,675,270.691646 |
| 201014SP_E1_st_14_devitelined | 1456 | 201014SP_E1_st_14_devitelined::1456 | 2026-08-04_121049 | 201014SP_E1_st_14_devitelined - Not_self_Merged.csv | 111,733,775.963161 |
| 201014SP_E1_st_14_devitelined | 2318 | 201014SP_E1_st_14_devitelined::2318 | 2026-08-04_121522 | 201014SP_E1_st_14_devitelined - Not_self_Merged.csv | 6,035,022,834.77704 |
| 201024SP_E1_st_14_devitelined_albino | 2448 | 201024SP_E1_st_14_devitelined_albino::2448 | 2026-08-05_024308 | 201024SP_E1_st_14_devitelined_albino - Not_self_Merged.csv | 412,340,733.315302 |

In `scripts/create_spectral_analysis_tables.py`, `INPUT_COLUMNS` includes `TotalBandPower`; `read_header` requires it; `build_embryo_records` parses it for summary statistics; and `write_cell_master` copies the source row fields into the master. This script has no feature calculation, scaling, or cell-record merge. `GlobalID` is formed from embryo and cell identifiers. Therefore the table-assembly script did not originate the extremes; they are present in the per-cell source outputs.

**Confirmed from the attached PDF, `extractSpectralFeatures`:** The code calls MATLAB CWT with analytic Morlet (`amor`), `power = abs(cfs).^2`, `meanPower_raw = mean(power,2)`, then linearly interpolates to `FrequencyAxis = linspace(0.02,0.20,80)` with out-of-range values set to zero and clamps the interpolated values at zero. It sets `totalPower = sum(meanPower)` and assigns that value to `TotalBandPower` when positive. No frequency-bin-width integration, amplitude rescaling, or duration normalization appears in this documented feature formula. The normalized distribution `P = meanPower/(totalPower+eps)` is used for entropy and is not the TotalBandPower calculation.

**Confirmed from the attached PDF, `preprocessCellTrace`:** `F0 = prctile(traceUniform,10)`, `F0 = max(F0,eps)`, then `dFF = (traceUniform-F0)/F0`. This establishes a plausible amplification path for low baselines. **Potential issue:** very small baseline values could produce very large ΔF/F and CWT power. **Unresolved:** Raw traces for the 11 rows were not located in the inspected data tree, so it is not possible to determine whether their fluorescence baselines trigger this mechanism. The PDF specifies the algorithm but does not supply the run-specific upstream code/version or a physical calibration for the sum. Git history returned no commits for the table-building script; that absence does not establish that no upstream processing-version changes occurred. The 11 high values are confirmed in source result files; they are not thereby proven erroneous.

## 7. Recommended next steps

1. **Before group comparisons:** Locate the raw traces and run-specific feature-generation code/configuration for the 11 extreme source rows. Check their P10 baseline, ΔF/F amplitude, and 80-bin mean power using the documented method; establish any unit/calibration convention. Do not remove or transform them based only on magnitude.
2. **Before group comparisons:** Define and document the denominator for each feature. Preserve the distinction between 417,650 cell rows and 415,186 valid observations per feature, and disclose the 2,464 all-feature-missing rows.
3. **Before group comparisons:** Confirm whether inference will use cells or embryos and pre-specify how any subsequently verified data-quality concern will be handled. This audit performed no hypothesis tests and made no exclusion decisions.
4. **Document separately:** Preserve the verified provenance chain (`InputFile`, `SourceFile`, `RunID`, `EmbryoID`, `CellID`) and record that the tables are internally consistent. Keep the high values flagged for source-method review while their status is unresolved.

Evidence inspected: `scripts/create_spectral_analysis_tables.py` (constants and functions `numeric_or_none`, `percentile`, `feature_summary`, `build_embryo_records`, `write_cell_master`, `write_embryo_table`); `README.md` (no method detail); both existing master CSVs; existence of all 179 `InputFile` source paths; and direct row-level checks of the 11 extreme-value source records. No files other than this requested audit document were intentionally written.
