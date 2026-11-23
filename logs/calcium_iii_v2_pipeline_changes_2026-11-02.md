# Calcium III v2 pipeline: changes, validation, run order

Date: 2026-11-02 (revised the same day). **Guiding rule: wherever a method choice exists, v2 follows the PDT paper (CB Calcium III Working draft, Methods 2.1–2.5)** so the spectral results can be integrated into it. **Group membership comes from `data10/Embryo_Ref_Table.csv`**, the definitive grouping, not from folder or file names. This follows `calcium_iii_pipeline_audit_2026-11-01.md`. The v1 scripts that produced `data10/HPC_spectral`, `data/master_tables/group_tables` and `scr10/CALCIUM/data/mapping` are archived unchanged in `scripts/archive_v1_pre_2026-11-01/`. v2 writes only to new folders (`HPC_spectral_v2`, `group_tables_v2`, `mapping_v2`), so v1 outputs are never overwritten.

## 1. Spectral pipeline (`scripts/wavelets/`)

| Item | v1 | v2 |
| --- | --- | --- |
| Frequency axis | 0.02–0.20 Hz, 80 bins (above Nyquist) | 0.02–**0.15** Hz, 80 bins (`CALCIUM_FMIN/FMAX/NFREQ`) |
| Power | \|W\|² (grows with scale → 53–71% of cells pinned at 0.02 Hz) | **\|W\|²/scale** (Liu et al. 2007 rectification) |
| Segments | concatenated, one CWT across splices | **each segment transformed separately**; mean power uses only **cone-of-influence-valid** samples (COI = √B·C/f s; 61 s at 0.02 Hz) |
| Gap rule | split at any interval > 2×median (> 6 s), so one dropped frame split a trace | **paper rules:** up to **3 missing frames** bridged by interpolation (TrackBinaries' 3-frame gap; ND.T step ≤ 4); **refocus pauses** split the trace, segments processed separately, and the point just before and after each pause is ignored. A pause = frame step of 1 with a time jump above 1.5× the frame interval (23–32 s observed). (`CALCIUM_MAX_MISSING_FRAMES`, `CALCIUM_TRIM_AT_SPLIT`; fallback without ND.T: split > `CALCIUM_GAP_SPLIT_SECONDS` = 10 s) |
| Inclusion | every segment needed ≥ 100 resampled points | a segment must be > 2×COI at 0.02 Hz (~122 s ≈ 37 frames); a cell needs ≥ **100 observed frames** in usable segments, the paper's tracking criterion (`CALCIUM_MIN_CELL_FRAMES`) |
| Time base | rounded `hh:mm:ss` (±0.5 s jitter); median step = 3 s | **embryo-wide frame→time map**: time regressed on `ND.T` within pause-free runs. Removes rounding jitter; estimated frame interval is 3.3–3.6 s depending on the embryo |
| Baseline | one P10 across concatenated segments | **paper baseline (a): average of the lower 10% of intensities**, per segment (the paper processes refocus segments separately); dF/F = (F − F0)/F0 (`CALCIUM_BASELINE_FRACTION`) |
| Embryo selection | every CSV in the group folder | embryos listed for the group in **`Embryo_Ref_Table.csv`** (`CALCIUM_GROUP`), found anywhere under `merged_data`; unlisted files ignored; missing ones reported in the manifest |
| Out-of-order rows | cell rejected ("not monotonic") | rows sorted by time (merged multi-layer tracks append gap-fill rows) |
| Mixed time formats | 190724SP_E2 lost ~567 cells | parsed value by value, including the `00:00.-01` first-frame artefact |
| Outputs | CSV + NPZ + 200-dpi PNG per cell (~32 GB) | CSV with QC columns (`UsedDuration_min, FrameInterval_s, FramesObserved, FramesUsed, NSegments, NBridgedGaps, NSplitGaps`), `skipped_cells.csv` with reasons, NPZ with `ValidSamples`, **`run_manifest.json`** (all parameters, versions, git commit). PNGs off by default in the sbatch files (`CALCIUM_MAKE_TRACE_PLOTS=true` to enable) |

The six feature definitions (dominant frequency, centroid, spread, Shannon entropy, summed band power, peak power fraction) are unchanged. TotalBandPower values are not comparable to v1, because power is now scale-rectified.

### Validation (in memory, before any production run)

- **Synthetic sinusoids**, sampled at 3.3 s with integer-rounded timestamps, 3% dropped frames and a 30 s pause. True frequencies 0.025–0.12 Hz were recovered within one bin (0.0249, 0.0397, 0.0595, 0.0792, 0.1006, 0.1204). A 0.14 Hz input was misassigned to 0.15 Hz, which is too close to Nyquist; treat dominant frequencies above ~0.13 Hz with caution.
- **Clean equal-amplitude sinusoids**: peak power is roughly flat across frequency (0.0098 at 0.025 Hz, 0.0069 at 0.10 Hz). Before rectification it was proportional to scale.
- **White noise**: 0% of 200 noise traces have their dominant frequency at 0.02 Hz (median 0.097 Hz), so the low-frequency pile-up is gone.
- **Pilot embryos** (real data):

  | Embryo | Cells v1 → v2 | Share at 0.02 Hz, v1 → v2 | Median dominant frequency, v1 → v2 |
  | --- | --- | --- | --- |
  | 180209SP_E4_st_14 | 1,355 → 1,627 | 64% → 2% | 0.020 → 0.099 Hz |
  | 190724SP_E2_st_13_devitelined | 943 → 1,716 | 59% → 20% | 0.020 → 0.081 Hz |

  These were measured with the earlier 10 s gap rule. With the paper's 3-frame rule, cell counts rose to **1,910** and **2,005**. On 180731SP_E1 (2 refocus pauses) both pauses were detected, 1,050 cells were split into 2 segments, and 2,114 cells were kept.

  Run time is about 40 s per embryo without PNGs.

**Interpretation caveat:** in 180209SP_E4 the median spectral entropy (6.28) is close to the flat-spectrum maximum (log2 80 = 6.32), and its median dominant frequency matches what white noise gives. Many cells may have no significant rhythm, in which case DominantFreqHz simply locates the largest noise peak. A red-noise (AR(1)) significance test per cell (Torrence & Compo 1998) is the recommended next addition before dominant frequency is interpreted biologically.

## 2. Group tables (`scripts/create_hpc_group_master_tables.py`)

Reads `HPC_spectral_v2` and writes `data/master_tables/group_tables_v2/` by default, with the v2 QC columns. It refuses to write into a non-empty folder. Override with `CALCIUM_SPECTRAL_ROOT`, `CALCIUM_TABLE_DIR` and `CALCIUM_TABLE_VERSION=v1`.

## 3. Mapping (`scripts/mapping/`)

- The 9-px CNPU grid, origin and extent are unchanged, to match the composite neural plate in the PDT paper.
- **The paper's CNP maps are now the primary maps.** `<feature>_GroupSpatialMaps.png` = CNPU mean over all cells, smoothed over 49 CNPUs, z-scored, on a fixed −1…1 jet scale (paper Figs 2–4). The embryo-balanced version moved to `<feature>_EmbryoBalanced_SmoothedZ_GroupSpatialMaps.png`.
- **Added the paper's CNPU statistic.** The paper's Methods 2.3 averages all cells from all embryos within a CNPU. v1 used embryo medians, then a median across embryos. Both are now written: `PooledMean`, `PooledSmoothed` and `PooledSmoothedZ` (7×7 / 49-CNPU mean, z-scored) sit alongside the existing `GroupMedian`, `SmoothedMedian` and `SmoothedZ`. There is a new per-group map, `<feature>_PooledMean_SmoothedZ_GroupSpatialMaps.png`.
- **New `combine_group_maps.py` + `combine_group_maps.sbatch`**, run after the four group jobs:
  - 2×2 maps on **shared colour scales** across groups. z-maps use the paper's fixed −1…1 scale.
  - Neural vs non-neural assignment from the stage-matched Sox2 polygon.
  - Embryo × region summaries.
  - Paper-style cell-level K-S + Cohen's d with Bonferroni (significant if p_adj < 0.05 and |d| > 0.2).
  - Embryo-level paired Wilcoxon as a robustness check.
  - `Registration_QC.csv`: affine R² of image position → template position per embryo.
  - `--exclude-embryos` for sensitivity runs.
- The four mapping sbatch files now read `group_tables_v2` and write to `scr10/CALCIUM/data/mapping_v2`.

## 4. Data findings and inventory (run `python3 scripts/embryo_inventory.py`)

`scripts/embryo_inventory.py` writes `data/embryo_inventory.csv`. It lists the raw, boundary and spike/Cy5 files for every reference embryo. Rerun it after adding files. Results for 2026-11-02:

- **All 179 reference embryos have a SuperEmbryo boundary file.** The only embryo that previously failed to map, `200217RN_E3_st_13_devitelined_albino001`, is not in the reference table.
- **Two reference embryos have no merged raw CSV** and must be retrieved from the original source:
  - `191020RNSP_E6` (Stage 14);
  - `191030SPRN_E3` (Stage_13_dev). The S13_dev folder only appeared complete at 46 because it held `200217SP_E3`.
- **One raw CSV sits in the wrong folder:** `200217SP_E3_st_13_albino` is a Stage 13 embryo per the reference table but is in `merged_data/stage_13_devitellinated`. v2 selects it correctly by reference anyway; moving it is optional, for tidiness.
- **Not in the reference table:** the raw file `stage_13/200217RN_E3_st_13_devitelined_albino001` (ignored by v2), and boundary files for 190905SP_E2, 191103RN_E7 and 190918SP_E1.
- **Only 79 of 179 embryos have a local `ConvergeSpikeCy5` event file.** The full PDT event outputs will be needed for event–spectral integration.
- **Registration errors** in boundary files made by the previous student (`200114SP_E1` R² 0.09, `180405SP_E4` R² 0.26) cannot be fixed here. Run the combined step with and without them (`EXCLUDE=...`).

## 5. Run order

```bash
# 0. Add the two missing merged CSVs, then: python3 scripts/embryo_inventory.py
# 1. Spectral (one job per group; embryos chosen from Embryo_Ref_Table.csv)
sbatch scripts/wavelets/run_pipeline_stage_13.sbatch
sbatch scripts/wavelets/run_pipeline_stage_13_devitellinated.sbatch
sbatch scripts/wavelets/run_pipeline_stage_14.sbatch
sbatch scripts/wavelets/run_pipeline_stage_14_devitellinated.sbatch
# 2. Tables (after all four finish)
python3 scripts/create_hpc_group_master_tables.py
# 3. Mapping, one job per group, then the combined step
sbatch scripts/mapping/spatial_mapping_stage_13.sbatch   # etc. for the other three
sbatch --dependency=afterok:<ids> scripts/mapping/combine_group_maps.sbatch
```

Check each `run_manifest.json`, and the `Totals:` line at the end of each wavelet `.out` log, which gives kept and skipped cells by reason.
