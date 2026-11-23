# Calcium III spectral analysis: next steps toward the SfN poster and the PDT manuscript

Written 2026-11-02. Inputs at the start: the v3 spectral run (`data10/HPC_spectral_v3`), v3 group tables (`data/master_tables/group_tables_v3`), and the audited v3 maps (`scr10/CALCIUM/data/mapping_v3`, smoother limited to occupied CNPUs, combined maps at ≥ 2 embryos).

**Goal.** Spectral analysis is a second layer on top of the PDT event analysis. It shows where in the neural plate activity is rhythmic, at what frequency, and whether that changes with stage (13 vs 14) or with devitellination. A claim goes on the poster or into the manuscript only when it holds with the embryo as the unit of analysis.

**Rules carried over:**
- Methods follow the PDT paper wherever a choice exists (9-px CNPUs, pooled mean, 49-CNPU smoothing, z-scores on a −1…1 jet scale, K-S plus Cohen's d with Bonferroni). Embryo-level statistics are added alongside the paper's methods; they never replace them.
- `data10/Embryo_Ref_Table.csv` is the only source of group membership.
- No existing output is overwritten. Every new step writes to `scr10/CALCIUM/data/analysis_v3/<step>/`. The rhythmic-only maps go to `mapping_v3_rhythmic`.

## Script layout

All scripts are written but **none has been run**. Every script prints its own documentation with `--help`.

| Folder (`CALCIUM/scripts/…`) | Script | Runs on | Purpose |
|---|---|---|---|
| `analysis_common/` | `common.py` | (imported) | Paths, groups, region assignment, BH/Holm, the 2×2 permutation model, the paper's CNPU method, plot style |
| | `run_step.sbatch` | Slurm | Generic wrapper for any single-script step |
| `step1_rhythmicity/` | `rhythmicity_test.py` + `run_rhythmicity.sbatch` | Slurm, 16 CPUs per group | AR(1) red-noise Monte Carlo test for each cell |
| | `build_rhythmic_tables.py` | Slurm or login node | FDR, the Rhythmic flag, per-embryo % rhythmic, rhythmic-only group tables |
| | `map_rhythmic.sbatch` | Slurm | Unchanged v3 mapping + combine step run on rhythmic cells only |
| `step2_cell_flow/` | `cell_flow.py` | wrapper | Cells retained at each step, per embryo and per group |
| | `selection_bias.py` | wrapper | Whether skipped and kept cells differ in position |
| `step3_feature_selection/` | `feature_correlations.py` | wrapper | Which features are redundant |
| `step4_embryo_stats/` | `ap_gradient.py` | wrapper | Anterior–posterior (and medio-lateral) gradient per embryo |
| | `group_effects.py` | wrapper | Stage × membrane effects on embryo medians, the paper-style cell table, optional mixed model |
| | `region_concordance.py` | wrapper | Neural vs non-neural: cell-level vs embryo-level agreement |
| `step5_event_integration/` | `find_event_files.py` | login node | Index of PDT event files and a list of missing ones |
| | `event_spectral_join.py` | wrapper | Cell-level join of events and spectral features; agreement tests |
| | `event_maps.py` | wrapper | Event metrics on the CNPU grid; spatial agreement with spectral maps |
| `step6_figures/` | `pipeline_summary_panels.py` | wrapper | Poster and manuscript panels (PDF and PNG) |

Shorthand used below, run on Bora:

```bash
module load slurm/bora
S=/sciclone/home/cabare01/CALCIUM/scripts
RUN="$S/analysis_common/run_step.sbatch"
```

## Step 0: Housekeeping (½ day, do first)

1. **Commit the current state** so every result can be traced to a commit. Many v3 scripts are uncommitted, and the wavelet manifest records "scripts/wavelets has uncommitted changes".
   ```bash
   cd /sciclone/home/cabare01/CALCIUM && git checkout -b analysis-v3 && git add scripts PROJECT_PLAN_spectral_next_steps.md && git commit -m "v3 pipeline + analysis steps 1-6"
   ```
2. **Install statsmodels** (needed only for the optional mixed model in step 4b):
   `conda install -n wavelet-analysis statsmodels`.
3. **Confirm the template orientation:** is anterior at small template Y (the top of the maps)? `ap_gradient.py` assumes `--anterior top`. Check against the Sox2 template or the PDT paper's figures.
4. **Confirm the event-file columns** in the PDT paper: which column is the event rate (`Th_SD_1`…`Th_SD_5`? `Counts_M`?) and which is event magnitude. Write the choice here: primary event-rate column = `________`, magnitude column = `________`.
5. **Start retrieving event files now, because it is the slowest dependency.** Run step 5a (seconds on a login node), then send `missing_event_files.csv` to whoever holds the full PDT outputs. Put retrieved files in `data10/pdt_events/` (the folder does not exist yet; create it).
   ```bash
   python $S/step5_event_integration/find_event_files.py
   ```

## Step 1: Which cells are actually rhythmic? (gating step, about 1–2 days of wall time)

**Why.** Median spectral entropy is about 6.15, close to the flat-spectrum maximum of 6.32. For many cells, DominantFreqHz may only locate the largest noise peak. Every frequency claim depends on this step.

**Method.** For each cell, `rhythmicity_test.py`:
1. Rebuilds the exact dF/F segments with the v3 preprocessing. It reads the raw file list from `run_manifest.json`, and it stops if the current pipeline settings differ from the manifest.
2. Estimates the AR(1) coefficient α as (r1 + √r2)/2.
3. Simulates 1,000 AR(1) surrogate series with the same segment lengths, dt, linear detrend and SD.
4. Passes them through the identical CWT, |W|²/scale rectification and cone-of-influence masking.

The p-value tests the maximum observed/red-noise power ratio across all 80 frequencies, so it already corrects for testing many frequencies. `build_rhythmic_tables.py` then applies BH-FDR within each embryo: a cell is **Rhythmic** if q < 0.05. Built-in check: the observed spectrum is recomputed for the first three cells of every embryo and must equal `SpectralLibrary.npz`. If it doesn't, the run stops.

**Run:**
```bash
# 1. Pilot: one embryo, 50 cells, to time it and confirm the self-check passes
sbatch --time=01:00:00 --export=ALL,GROUP=S14,EXTRA="--embryos 180209SP_E4_st_14 --max-cells 50" $S/step1_rhythmicity/run_rhythmicity.sbatch
```
- **Check the pilot:** the `.out` log shows `seconds` for 50 cells. Projected hours per group ≈ (seconds / 50) × about 75,000 cells / 16 processes / 3600. If the projection is over about 40 h, rerun with `N_SURROGATES=500` (p floor 0.002) or raise `--cpus-per-task`.
- Look at `analysis_v3/step1_rhythmicity/pilot/S14/180209SP_E4_st_14_rhythmicity.csv`. α should mostly lie between 0.3 and 0.95, and P_global should vary across cells.

```bash
# 2. Full run (four jobs). Finished embryos are skipped if a job is resubmitted after a timeout
for g in S13 S13_devitellinated S14 S14_devitellinated; do
  sbatch --job-name=rhythm_$g --export=ALL,GROUP=$g $S/step1_rhythmicity/run_rhythmicity.sbatch
done
# 3. After all four finish: FDR, summaries, rhythmic-only tables
sbatch --job-name=step1b --export=ALL,CMD="python $S/step1_rhythmicity/build_rhythmic_tables.py" $RUN
# 4. Rhythmic-only maps (new folder mapping_v3_rhythmic)
sbatch $S/step1_rhythmicity/map_rhythmic.sbatch
```

**Check:**
- `missing_cells.csv` should be empty.
- The printed share of cells "at the Monte Carlo p floor" should be small.
- Look at `rhythmicity_by_group.csv`.

**Decide:**

| Result | What it means | What to do |
|---|---|---|
| ≥ about 30% of cells rhythmic | Frequency features are meaningful | Report frequency features for rhythmic cells (main), all cells (supplement) |
| about 5–30% | Rhythmicity is the headline | Lead with % rhythmic and where it occurs; give frequency for rhythmic cells only |
| < about 5% | Little rhythm above red noise | Headline becomes power and broadband activity (TotalBandPower, PeakPowerFraction); state the negative result about rhythm explicitly |

## Step 2: Cell flow and selection bias (1 hour)

```bash
sbatch --job-name=step2_flow --export=ALL,CMD="python $S/step2_cell_flow/cell_flow.py" $RUN
sbatch --job-name=step2_bias --export=ALL,CMD="python $S/step2_cell_flow/selection_bias.py" $RUN
```

**Outputs:**
- `cell_flow_by_group.csv`: Tracked → passed spectral QC → mapped → in grid → rhythmic. Current totals: about 577k tracked, about 300k passed QC, 191k mapped.
- `skip_reasons_by_group.csv`.
- `selection_bias_tests.csv`: kept minus skipped cells within each embryo for % neural, median AP position and medio-lateral distance (Wilcoxon across embryos, Holm-corrected).

**Decide:** if a bias test is significant with a sizeable median difference (for example > 5 percentage points neural, or > 50 px), state it as a limitation and avoid conclusions about the under-represented region. The flow numbers go into the Methods and the poster's cell-flow panel.

## Step 3: Pick the features to show (30 min)

```bash
sbatch --job-name=step3 --export=ALL,CMD="python $S/step3_feature_selection/feature_correlations.py" $RUN
# after step 1:
sbatch --job-name=step3r --export=ALL,CMD="python $S/step3_feature_selection/feature_correlations.py --cells rhythmic" $RUN
```

**Output:** `redundant_pairs.csv` and `feature_correlations.pdf`. These give the within-embryo median Spearman ρ and the between-embryo ρ.

**Decide:** pick **two** features for the main figures:
- one frequency feature (DominantFreqHz, or SpectralCentroid if the two are redundant and centroid is the more stable of the pair);
- one power feature (TotalBandPower or PeakPowerFraction) with |ρ| < 0.5 against the first.

Record the choice here: `________` and `________`. The other features go to the supplement.

## Step 4: Embryo-level statistics (half a day including review)

Run each analysis on all cells, and again on rhythmic cells once step 1 is done (`--cells rhythmic`).

```bash
sbatch --job-name=step4_ap  --export=ALL,CMD="python $S/step4_embryo_stats/ap_gradient.py --axis ap --anterior top" $RUN
sbatch --job-name=step4_ml  --export=ALL,CMD="python $S/step4_embryo_stats/ap_gradient.py --axis ml" $RUN
sbatch --job-name=step4_eff --export=ALL,CMD="python $S/step4_embryo_stats/group_effects.py --lmm" $RUN
python $S/step4_embryo_stats/region_concordance.py                                   # all cells (seconds)
python $S/step4_embryo_stats/region_concordance.py --mapping-root /sciclone/home/cabare01/scr10/CALCIUM/data/mapping_v3_rhythmic --label rhythmic
```

**4a. Gradient** (`gradient_ap_group_tests.csv`, `gradient_ap_stage_membrane.csv`)
- Tests the anterior-high / posterior-low DominantFreqHz pattern seen in all four maps.
- Embryos need at least 100 neural cells spanning at least 200 px along the axis.
- Claim a gradient in a group if its Holm-corrected p < 0.05 and most embryos (≥ about 70%) agree in sign.
- Report the stage and membrane terms as a change in gradient.

**4b. Group effects**
- `stage_membrane_effects.csv` is the inferential result: embryo medians, a 2×2 model, Freedman–Lane permutation p-values, Holm correction across outcomes.
- `pairwise_embryo_level.csv` gives the four planned contrasts.
- `pairwise_paper_style_cells.csv` is the PDT-paper-style table, labelled descriptive.
- `lmm_effects.csv` is a cell-level mixed model and should agree in direction with the permutation results.
- One pattern to test: devitellinated groups showed lower DominantFreqHz and higher TotalBandPower at both stages in the cell-level means.
- **To add before final results: recording year/date.** Recording year is unbalanced across groups: S13 is 31/7/6 embryos across 2018/19/20, S14_dev is 5/14/22. TotalBandPower rises with year inside the intact groups. The membrane effect holds within each year, but the final model should adjust for it:
  - add year as a covariate in `two_by_two`, or a stratified permutation within year;
  - add `(1|RecordingDate)` to the mixed model.

  `group_effects.py` does not do this yet.

**4c. Region concordance** (`region_concordance_*.csv`)
- Each group × feature is labelled Supported, Cell-level only, Embryo-level only, or Not supported.
- On the current all-cell results, nothing is Supported: the five tests that pass the paper rule are all cell-level only.
- Only Supported rows are stated as findings.

## Step 5: Link spectral results to PDT events (depends on event-file retrieval)

```bash
python $S/step5_event_integration/find_event_files.py          # rerun after adding files to data10/pdt_events
sbatch --job-name=step5_join --export=ALL,CMD="python $S/step5_event_integration/event_spectral_join.py" $RUN
sbatch --job-name=step5_maps --export=ALL,CMD="python $S/step5_event_integration/event_maps.py --event-columns <rate column> <magnitude column>" $RUN
```

1. **Join check.** In `event_match_by_embryo.csv`, embryos below 50% matched are excluded. A low rate means the event file came from a different tracking run, as with 180405SP_E4 and 200114SP_E1; add any such embryos to the plan's exclusion notes.
2. **Validation of the two methods.** In `event_spectral_rho_tests.csv`, event rate and magnitude should correlate positively with TotalBandPower and PeakPowerFraction across embryos. If they do, the two methods support each other, which reviewers will look for. If they don't, find out why before using the spectral results.
3. **Rhythmic vs non-rhythmic.** `event_rhythmic_contrast.csv` tests whether rhythmic cells have more or larger events.
4. **Spatial agreement.** `maps/spatial_concordance_tests.csv` correlates event and spectral values across CNPUs within each embryo. The direct group-map correlation is descriptive only.

**Minimum to proceed:** about 79 of 179 embryos have event files now. Run step 5 on what exists for the poster, labelled with the n per group. Rerun when the remaining files arrive, for the manuscript.

## Step 6: Figures (1 day including layout)

```bash
sbatch --job-name=step6 --export=ALL,CMD="python $S/step6_figures/pipeline_summary_panels.py --features <feat1> <feat2> --event-column <rate column>" $RUN
# rhythmic-only maps:
python $S/step6_figures/pipeline_summary_panels.py --panels maps --mapping-root /sciclone/home/cabare01/scr10/CALCIUM/data/mapping_v3_rhythmic --features <feat1> <feat2>
```

Panels are written to `analysis_v3/step6_figures/` as `fig_*.pdf` (vector) and `.png`:
- cell flow;
- example spectra against red noise;
- 2×2 maps for both features;
- % rhythmic map;
- anterior–posterior gradient dot plots;
- stage × membrane dot plots;
- events vs spectral.

Every dot plot shows one dot per embryo with the group median. Colour encodes membrane state (blue intact, orange devitellinated) and x-position encodes stage. The z-maps keep the paper's jet −1…1 scale.

## Step 7: Synthesis into the poster and the manuscript

**Claims table.** Fill this table first; the poster and manuscript text follow from it. A row is a finding only when the embryo-level column is significant.

| # | Claim | Embryo-level evidence (file, statistic, n embryos) | Paper-style cell evidence | Status |
|---|---|---|---|---|
| 1 | X% of cells show rhythmic activity above red noise | `rhythmicity_by_group.csv` | — | |
| 2 | Anterior–posterior gradient of <feat1> in the neural plate | `gradient_ap_group_tests.csv` | — | |
| 3 | Devitellination changes <feat> | `stage_membrane_effects.csv` (Membrane) | `pairwise_paper_style_cells.csv` | |
| 4 | Stage changes <feat> | `stage_membrane_effects.csv` (Stage) | same | |
| 5 | Spectral power tracks PDT event rate | `event_spectral_rho_tests.csv` | — | |
| 6 | Neural vs non-neural difference | `region_concordance_*.csv` (Supported only) | `Region_Comparison_CellLevel.csv` | |

**Poster structure** (one message; about 6 panels):
1. **Question:** where, and how rhythmically, does the neural plate signal, and does that depend on stage and the vitelline membrane?
2. **Pipeline schematic:** tracking → PDT events → wavelet → red-noise test → template mapping. Include the cell-flow bar.
3. **Validation:** synthetic-sinusoid and white-noise checks from the v2 notes, the example spectra, and % rhythmic.
4. **Main figure:** 2×2 maps of the two chosen features (rhythmic cells), plus the % rhythmic map.
5. **Quantification:** dot plots of the anterior–posterior gradient and the stage × membrane effects.
6. **Events vs spectral panel**, then conclusions and limitations: registration, coverage, two excluded embryos, partial event files, selection bias from step 2.

**Manuscript** (spectral analysis integrated into the PDT draft):
- **Methods:** a new subsection after event detection, "Time–frequency analysis". Cover:
  - preprocessing shared with the PDT methods (gaps, refocus splits, baseline);
  - the CWT (cmor1.5-1.0, 0.005–0.08 Hz, |W|²/scale, cone of influence);
  - the six features;
  - the red-noise test;
  - mapping, using the paper's CNP method;
  - statistics: the paper's K-S/d plus the embryo-level permutation model and mixed model.

  Base the text on `Calcium III Wavelet Methods (1).pdf` and the v2/v3 change notes in `logs/`. That PDF couldn't be read when this plan was written, so check it against this list.
- **Results:** spectral results follow the event results and tie back to them: the step 5 validation first, then rhythmicity, the gradient, and the group effects.
- **Supplement:**
  - registration R² table;
  - coverage maps;
  - `mapping_exclusions.csv`;
  - cell flow and selection bias;
  - feature correlations;
  - all-cell versions of the rhythmic-only maps;
  - cell-level vs embryo-level concordance tables.

## Order and timing

| When | Work |
|---|---|
| Day 1 | Step 0, step 5a (request event files), step 1 pilot, steps 2–3 on all cells |
| Days 1–3 | Step 1 full run (Slurm). In parallel: steps 4a–4c on all cells, draft Methods text |
| Day 3–4 | Step 1b + rhythmic maps; rerun steps 3–4 with `--cells rhythmic`; fill the claims table |
| Day 4–5 | Step 5 on the available event files; step 6 figures; poster layout |
| After the poster | Remaining event files → rerun step 5; manuscript Results and supplement |

Work back from the SfN poster deadline. Steps 1 and 5 are the only ones that can slip by days, so start both on day 1.

## Risks

| Risk | Effect | Mitigation |
|---|---|---|
| Few cells pass the red-noise test | Frequency claims weaken | The step 1 decision table; power and rhythmicity become the headline |
| Step 1 runtime too long | Delays everything after it | Pilot timing; `N_SURROGATES=500`; resubmission skips finished embryos |
| Event files missing or IDs mismatched | Step 5 is partial | Request on day 1; n per group stated on the poster; join-rate check |
| Anterior–posterior orientation assumed wrong | Gradient sign reversed | Step 0.3 check before step 4a |
| Selection bias toward a region | Maps not representative | Step 2b; state as a limitation |
| Cell-level significance from pseudoreplication | False positives | Only embryo-level-supported claims (step 4c) |
