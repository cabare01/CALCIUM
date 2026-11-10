# CALCIUM: wavelet analysis of calcium activity in the early *Xenopus* neural plate
An HPC pipeline that turns single-cell calcium-imaging traces from 179 *Xenopus laevis* embryos into time–frequency (continuous wavelet) features. 
It then tests each cell against a red-noise null, maps the results onto a common template of the neural plate, and compares them with the event-based peak-detection results from the same cells.
The pipeline runs on the William & Mary HPC cluster (Bora/Hima) under Slurm, in a specifically created conda environment.

| | |
|---|---|
| **Biological question** | Where in the neural plate is calcium activity structured in time, and does that depend on developmental stage (13 vs 14) or on removing the vitelline membrane? |
| **Input** | 576,943 tracked cells (GCaMP6m fluorescence, 60 min per embryo) in 4 groups: stage 13 / 14 × membrane intact / removed |
| **Output** | Per-cell spectral features, per-cell red-noise test results, composite neural plate maps, embryo- and cell-level statistics |

## Contents
1. [Background](#1-background)
2. [Pipeline overview](#2-pipeline-overview)
3. [Repository layout](#3-repository-layout)
4. [Inputs](#4-inputs)
5. [Outputs](#5-outputs)
6. [Software environment](#6-software-environment)
7. [Usage](#7-usage)
8. [Method parameters](#8-method-parameters)
9. [Validation and benchmarking](#9-validation-and-benchmarking)
10. [Results and reporting](#10-results-and-reporting)
11. [Milestones](#11-milestones)
12. [Known limitations](#12-known-limitations)
13. [Citations](#13-citations)

---

## 1. Background
Embryos expressing the calcium sensor GCaMP6m and a membrane marker (memRFP) were imaged for 60 minutes at early neurula. Every cell in the neural plate was segmented and tracked (NIS-Elements). After imaging, each embryo was stained for the neural marker Sox2 and registered onto a stage-matched template, the **composite neural plate (CNP)**. 
Cells from all embryos can then be compared in one coordinate system.
The earlier analysis of these data (Paudel, 2022; manuscript in preparation) detected discrete calcium events with **204 peak detection techniques (PDTs)**. 
It found that event counts change up to about 10-fold depending on the PDT. 
Peak counting also ignores the timing structure of a trace.
This project adds a complementary layer. It runs a **continuous wavelet transform (CWT)** on each cell's ΔF/F trace, summarises the time-averaged spectrum with six features, and asks whether each cell's spectrum exceeds an AR(1) red-noise null. 
The results are mapped onto the same CNP with the same method as the PDT study, so the two analyses can be compared cell by cell.

**Groups** (membership specified by `Embryo_Ref_Table.csv`):
| Group | Stage | Vitelline membrane | Embryos | Template |
|---|---|---|---|---|
| `S13` | 13 | intact | 46 | Stage 14 Sox2 |
| `S13_devitellinated` | 13 | removed | 46 | Stage 14 Sox2 |
| `S14` | 14 | intact | 46 | Stage 15 Sox2 |
| `S14_devitellinated` | 14 | removed | 41 | Stage 15 Sox2 |

**Unit of analysis.** The embryo is the unit of analysis. With ~50,000 cells per group, cell-level p-values are reported as descriptive only. A result counts as a finding only if it holds across embryos.

---

## 2. Pipeline overview
flowchart TD
    A["Tracked traces<br/>*Not_self_Merged.csv<br/>(179 embryos)"] --> P0["Stage 0 · raw_prep<br/>drop unused tracking columns"]
    P0 --> P1["Stage 1 · wavelets<br/>preprocess → CWT → features<br/>1 job per group"]
    P1 --> P2["Stage 2 · audit<br/>provenance + QC checks"]
    P1 --> P3["Stage 3 · group tables<br/>cell-level + embryo-level CSVs"]
    P3 --> P4["Stage 4 · mapping<br/>cells → CNP template → CNPU z-maps<br/>1 job per group, then combine"]
    P1 --> P5["Stage 5 · red-noise test<br/>1,000 AR(1) surrogates per cell<br/>1 job per group, 16 CPUs"]
    P5 --> P5b["build_rhythmic_tables<br/>BH-FDR within embryo"]
    P5b --> P5c["rhythmic-only maps +<br/>rhythmic-fraction maps"]
    P4 --> P5c
    B["Template boundaries<br/>+ Sox2 outlines"] --> P4
    C["PDT per-cell tables<br/>(204 PDTs)"] --> P6
    P4 --> P6["Stage 6 · figures<br/>SfN poster panels;<br/>PDT ↔ wavelet comparison"]
    P5c --> P6
    P6 -.-> N["Planned (scripts written, not yet run):<br/>cell flow · feature selection ·<br/>embryo-level statistics · event integration"]
```

Each stage is a standalone script with its own Slurm wrapper. It reads the previous stage's output from disk and writes to a new versioned folder, so any stage can be rerun on its own.
| # | Stage | Script(s) | Slurm wrapper | Parallelism | Resources requested |
|---|---|---|---|---|---|
| 0 | Raw-file cleanup (one-time) | `scripts/raw_prep/preprocess_stage_*.py` | `preprocess_stage_*.sbatch` (×4) | 1 job per group | 1 CPU, 5 GB, 3 h |
| 1 | Wavelet spectral features | `scripts/wavelets/pipeline.py` → `01_preprocess_data.py`, `02_spectral_features.py`, `03_build_results.py` | `run_pipeline_stage_*.sbatch` (×4) | 1 job per group | 4 CPUs, 32 GB, 24 h |
| 2 | Run audit | `scripts/audit/audit_spectral_run.py` | `analysis_common/run_step.sbatch` | 1 job | 2 CPUs, 32 GB, 6 h (production run used 1 CPU, 8 GB) |
| 3 | Group tables | `scripts/create_hpc_group_master_tables.py` | login node or `run_step.sbatch` | serial | — |
| 4a | Template mapping | `scripts/mapping/spatial_spectral_mapping.py` | `spatial_mapping_stage_*.sbatch` (×4) | 1 job per group | 1 CPU, 12 GB, 10 h |
| 4b | Combined 2×2 maps | `scripts/mapping/combine_group_maps.py` | `combine_group_maps.sbatch` | after 4a (`afterok`) | 1 CPU, 16 GB, 4 h |
| 5a | Red-noise test | `scripts/step1_rhythmicity/rhythmicity_test.py` | `run_rhythmicity.sbatch` (×4) | 16 worker processes per group | 16 CPUs, 64 GB, 48 h |
| 5b | FDR + rhythmic tables | `scripts/step1_rhythmicity/build_rhythmic_tables.py` | `run_step.sbatch` | after 5a | 2 CPUs, 32 GB, 6 h |
| 5c | Rhythmic-only maps | `scripts/mapping/*` on rhythmic cells | `step1_rhythmicity/map_rhythmic.sbatch` | after 5b | 1 CPU, 16 GB, 10 h |
| 5d | Rhythmic-fraction maps | `scripts/step1_rhythmicity/rhythmic_fraction_maps.py` | `run_step.sbatch` | after 4 and 5b | 2 CPUs, 32 GB, 6 h |
| 6a | Summary-level figures | `scripts/step6_figures/rednoise_bandpower_figures.py` | `run_step.sbatch` | after 4 and 5 | 2 CPUs, 32 GB, 6 h |
| 6b | PDT ↔ wavelet panels | `scripts/step6_figures/pdt_wavelet_panels.py`, `pdt_wavelet_panels_styled.py` | `pdt_wavelet_panels.sbatch` | after 3 and 5b | 2 CPUs, 64 GB, 4 h |
Shared code: `scripts/embryo_ref.py` normalises embryo names across file types and reads group membership. `scripts/analysis_common/common.py` holds paths, group definitions, the CNPU method, multiple-testing corrections and the plot style.

---

## 3. Repository layout

```
CALCIUM/
├── README.md
├── wavelet-analysis.yml                 # conda env, version-pinned (+ rebuild instructions)
├── wavelet-analysis.lock.txt            # conda explicit list: exact builds + md5
├── wavelet-analysis.pip-requirements.txt  # PyWavelets 1.9.0, hash-pinned
├── PROJECT_PLAN_spectral_next_steps.md  # analysis plan, steps 0–7, decision rules
├── PROJECT_SUMMARY_for_Claude_chat.md   # full methods + results summary with status labels
├── data/
│   └── group_tables/                    # stage 3 output (cell-level *_CW, embryo-level *_EW)
├── logs/
│   ├── calcium_iii_*_YYYY-MM-DD.md      # dated audit and change notes (the project log)
│   ├── audit_v3-594346.txt              # stage 2 audit of the production run
│   └── <job>-<jobid>.{out,err}          # Slurm logs of every production job
└── scripts/
    ├── embryo_ref.py                    # embryo-name normalisation + group lookup
    ├── embryo_inventory.py              # read-only inventory of input files per embryo
    ├── create_hpc_group_master_tables.py
    ├── raw_prep/                        # stage 0
    ├── wavelets/                        # stage 1
    ├── audit/                           # stage 2
    ├── mapping/                         # stage 4 (+ mapping_exclusions.csv)
    ├── step1_rhythmicity/               # stage 5
    ├── step6_figures/                   # stage 6 (+ bundled Liberation fonts, SIL OFL)
    ├── analysis_common/                 # common.py + generic run_step.sbatch
    ├── next_steps/                      # planned steps 2–5 (written, not yet run)
    └── archive_v1_pre_2026-11-01/       # superseded v1 code; kept for the audit trail only
```

> **Naming warning.** In this repo's group tables, `_CW` = **cell-level rows** and `_EW` = **embryo-level summary**. In the PDT files, CW/EW mean "cell-wise" vs "embryo-wise" threshold parameters. The two uses are unrelated.

---

## 4. Inputs

The raw data are unpublished and too large for GitHub (the PDT outputs alone are ~52 GB). They are on the W&M cluster under `/sciclone/data10/cabare01` (`~/data10`), and access is through the Saha lab. **The pipeline never modifies these inputs**, apart from the one-time column cleanup in stage 0.
| Input | Location | Format | Used by |
|---|---|---|---|
| Group membership (authoritative) | `data10/Embryo_Ref_Table.csv` | Columns `Stage 13`, `Stage_13_dev`, `Stage_14`, `Stage_14_dev`; 46/46/46/41 embryos | all stages |
| Tracked fluorescence traces | `data10/merged_data/<group>/*Not_self_Merged.csv` | One row per cell per frame: `Name, ID, ND.T, Time, PositionX/Y, MeanIntensity, …` | 0, 1, 5 |
| Per-embryo template registration | `data10/boundary_info/Boundaries/<embryo>-Boundary_SuperEmbryo.csv` | One row per last-frame cell: `Calcium_ID`, `CY5_Intensity` (Sox2), template centroid, outline pixels | 2, 4 |
| Sox2 template outlines | `data10/boundary_info/SuperEmbryoBoundaries/1_new_Stage_1{4,5}_Sox2_Boundary.csv` | Template pixel coordinates | 4, 6 |
| PDT results (previous study) | `data10/PDT_data/*_Perhour_Counts.csv`, `AllStage_Change_in_Mag_Data.csv` | One row per cell, event counts for each of the 204 PDTs | 6b |
| Mapping exclusions | `scripts/mapping/mapping_exclusions.csv` (in repo) | Embryos whose boundary file doesn't match their traces | 2, 4, 5 |

Imaging: nominally 0.33 Hz; the measured frame interval is 2.5–4.5 s (median 3.57 s). Recordings include refocus pauses (≈29 s), which the pipeline handles as segment breaks.

---

## 5. Outputs

Each production run writes to a new versioned folder (`HPC_spectral_v3`, `mapping_v3`, `analysis_v3/<step>/`). Existing outputs are never overwritten, and the scripts refuse to write into a non-empty output folder. Large outputs are kept on cluster scratch (`/sciclone/scr10/cabare01`, `~/scr10`).
| Stage | Output | Location | Contents |
|---|---|---|---|
| 1 | `run_manifest.json` | `data10/HPC_spectral_v3/<group>/` | Every parameter, package version, input file and git commit |
| 1 | `<embryo>/tables/singlecell_results_final.csv` | same | One row per kept cell: six spectral features + QC columns |
| 1 | `<embryo>/tables/skipped_cells.csv` | same | Every dropped cell and the reason |
| 1 | `<embryo>/SpectralLibrary.npz` | same | Time-averaged power spectrum (80 bins) per cell |
| 2 | `audit_v3-594346.txt`, `audit_v3_per_embryo.csv` | `logs/` | Provenance, consistency and QC checks |
| 3 | `stage_*_CW.csv` / `stage_*_EW.csv` | `data/group_tables/` | 300,023 cells (CW) and per-embryo summaries (EW) |
| 4 | Per-group maps + `All_MappedCells_with_CNPU.csv` | `scr10/CALCIUM/data/mapping_v3/<group>/` | Smoothed z-maps, embryo-balanced median, SD, coverage, `Registration_QC.csv` |
| 4 | Shared-scale 2×2 maps, region statistics | `scr10/CALCIUM/data/mapping_v3/combined/` | `SharedScale_<feature>_*.png`, `Region_Comparison_{Cell,Embryo}Level.csv` |
| 5 | `rhythmicity_all_groups.csv`, `_by_embryo.csv`, `_by_group.csv` | `scr10/CALCIUM/data/analysis_v3/step1_rhythmicity/` | Per-cell p/q values, `Rhythmic` flag, per-embryo % |
| 5 | Rhythmic maps | `scr10/CALCIUM/data/mapping_v3_rhythmic/`, `…/step1_rhythmicity/rhythmic_fraction_maps/` | Feature maps of rhythmic cells; % rhythmic per CNPU |
| 6 | Poster figures | `scr10/CALCIUM/data/analysis_v3/step6_figures/{sfn_poster,pdt_wavelet_panels}/` | PDF/PNG/SVG at print size; `figure_facts.json` and `numbers_for_poster.txt` hold every number quoted in captions |

**The six spectral features** (per cell, from the time-averaged, cone-of-influence-masked spectrum):

| Feature | Meaning |
|---|---|
| `DominantFreqHz` | Frequency of maximum power |
| `SpectralCentroid` | Power-weighted mean frequency |
| `SpectralSpread` | Power-weighted SD of frequency |
| `SpectralEntropy` | Shannon entropy of the normalised spectrum (bits; flat spectrum = log₂80 = 6.32) |
| `TotalBandPower` | Summed power over 0.005–0.08 Hz |
| `PeakPowerFraction` | Share of total power in the dominant bin |

---

## 6. Software environment

All stages from 1 onward run in the conda environment `wavelet-analysis` (linux-64, conda-forge only) that was created within miniforge3 on the W&M HPC. It is pinned at three levels:

| File | What it pins | Use it for |
|---|---|---|
| `wavelet-analysis.lock.txt` | Exact conda builds, with md5 checksums | **Exact rebuild (preferred)** |
| `wavelet-analysis.pip-requirements.txt` | PyWavelets 1.9.0 wheel, sha256-pinned | Exact rebuild (pip part) |
| `wavelet-analysis.yml` | Versions only (builds may differ, e.g. BLAS) | Portable rebuild |

```bash
module load miniforge3
source /sciclone/apps/miniforge3-24.9.2-0/etc/profile.d/conda.sh

# exact rebuild
conda create -n wavelet-analysis --file wavelet-analysis.lock.txt
conda run -n wavelet-analysis python -m pip install --no-deps --require-hashes \
    -r wavelet-analysis.pip-requirements.txt

# or: version-pinned rebuild
conda env create -f wavelet-analysis.yml
```

**Key tool versions**

| Tool | Version | Role |
|---|---|---|
| Python | 3.12.16 | |
| PyWavelets | 1.9.0 (PyPI) | Continuous wavelet transform (`cmor1.5-1.0`) |
| NumPy | 2.5.4 | Arrays, AR(1) surrogates |
| SciPy | 1.18.1 | Statistics, spatial neighbours |
| pandas | 3.1.0 | Tables |
| Matplotlib | 3.11.2 | Maps and figures |
| Slurm | 23.11.6 (`module load slurm/bora`) | Job scheduling |
| Miniforge3 | 24.9.2-0 (`module load miniforge3`) | conda |

Notes:
- **PyWavelets version quirk.** The 1.9.0 wheel reports `pywt.__version__ == "1.8.0"` (an upstream labelling bug). The pipeline records `importlib.metadata.version("PyWavelets")` in every run manifest instead.
- `statsmodels` is **not** in the environment. It is needed only for the optional mixed model in the planned step 4b (`group_effects.py --lmm`).
- Stage 0 uses only the Python 3 standard library and runs with the system `python3`.
- After any install or removal in the environment, regenerate all three files so they stay in sync.

---

## 7. Usage

### 7.1 Setup

```bash
git clone https://github.com/cabare01/CALCIUM.git ~/CALCIUM
cd ~/CALCIUM
# build the environment (section 6), then on the Bora login node:
module load slurm/bora
S=~/CALCIUM/scripts
RUN=$S/analysis_common/run_step.sbatch      # generic wrapper: runs $CMD in the env
```

> **Paths.** Scripts default to the original author's cluster paths (`/sciclone/home/cabare01/...`), and `#SBATCH --output` paths and email addresses are hard-coded. To run under another account, override the output paths with environment variables (below), and pass `--output=`/`--error=`/`--mail-user=` to `sbatch`.
>
> | Variable | Default | Read by |
> |---|---|---|
> | `CALCIUM_RAW_TRACE_ROOT` | `~/data10/merged_data` | stage 1 |
> | `CALCIUM_OUTPUT_ROOT` | `~/data10/HPC_spectral_v3/<group>` (**must not exist yet**) | stage 1 |
> | `CALCIUM_REF_TABLE` | `~/data10/Embryo_Ref_Table.csv` | stage 1 |
> | `CALCIUM_SPECTRAL_ROOT` | `~/data10/HPC_spectral_v3` | stages 3, 5, 6 |
> | `CALCIUM_TABLE_DIR` | `CALCIUM/data/group_tables` | stages 3–6 |
> | `CALCIUM_MAPPING_ROOT` | `~/scr10/CALCIUM/data/mapping_v3` | stages 5d, 6 |
> | `CALCIUM_ANALYSIS_ROOT` | `~/scr10/CALCIUM/data/analysis_v3` | stages 5, 6 |
>
> Wavelet parameters can be overridden too (`CALCIUM_FMIN`, `CALCIUM_FMAX`, `CALCIUM_NFREQ`, `CALCIUM_DETREND`, …; full list in the docstring at the top of `scripts/wavelets/pipeline.py`). Every value is recorded in `run_manifest.json`.

Scripts that take command-line options document them with `--help`; stages 0 and 1 are configured by environment variables instead.

### 7.2 Run the pipeline, stage by stage

The commands below reproduce the production (v3) run. Because outputs are never overwritten, set a **new** output folder for any rerun.

**Stage 0. One-time raw-file cleanup** (already done for the production data; it rewrites the merged CSVs in place after trimming the specified columns, so don't rerun it on the shared copy):
```bash
for g in 13 13_devitellinated 14 14_devitellinated; do
  sbatch $S/raw_prep/preprocess_stage_$g.sbatch
done
```

**Stage 1. Wavelet features** (four independent jobs separated by group):
```bash
for g in stage_13 stage_13_devitellinated stage_14 stage_14_devitellinated; do
  sbatch --export=ALL,CALCIUM_OUTPUT_ROOT=$HOME/data10/HPC_spectral_v4/$g \
         $S/wavelets/run_pipeline_$g.sbatch
done
```

**Stage 2. Audit** (should report `Issues (0)`):
```bash
sbatch --job-name=audit --export=ALL,CMD="python $S/audit/audit_spectral_run.py $HOME/data10/HPC_spectral_v4" $RUN
```

**Stage 3. Group tables:**
```bash
CALCIUM_SPECTRAL_ROOT=$HOME/data10/HPC_spectral_v4 CALCIUM_TABLE_DIR=<new dir> \
  python $S/create_hpc_group_master_tables.py
```

**Stage 4. Template mapping, then combine** (combine waits for all four group jobs to complete):
```bash
ids=$(for g in 13 13_devitellinated 14 14_devitellinated; do
        sbatch --parsable $S/mapping/spatial_mapping_stage_$g.sbatch; done | paste -sd:)
sbatch --dependency=afterok:$ids $S/mapping/combine_group_maps.sbatch
```

**Stage 5. Red-noise test.** Goes back to individual cell-level results and compares them relative to a red-noise trace (Torrence & Compo, 1998) to determine rhythmicity:
```bash
# pilot: one embryo, 50 cells (~3 min)
sbatch --time=01:00:00 --export=ALL,GROUP=S14,EXTRA="--embryos 180209SP_E4_st_14 --max-cells 50" \
       $S/step1_rhythmicity/run_rhythmicity.sbatch
# full run: one job per group (~9–12 h each)
ids=$(for g in S13 S13_devitellinated S14 S14_devitellinated; do
        sbatch --parsable --job-name=rhythm_$g --export=ALL,GROUP=$g $S/step1_rhythmicity/run_rhythmicity.sbatch
      done | paste -sd:)
b=$(sbatch --parsable --dependency=afterok:$ids --job-name=step1b \
      --export=ALL,CMD="python $S/step1_rhythmicity/build_rhythmic_tables.py" $RUN)
sbatch --dependency=afterok:$b $S/step1_rhythmicity/map_rhythmic.sbatch
sbatch --dependency=afterok:$b --job-name=step1d \
       --export=ALL,CMD="python $S/step1_rhythmicity/rhythmic_fraction_maps.py" $RUN
```

**Stage 6. Figures:**
```bash
sbatch --job-name=sfn_poster --export=ALL,CMD="python $S/step6_figures/rednoise_bandpower_figures.py" $RUN
sbatch $S/step6_figures/pdt_wavelet_panels.sbatch
```

### 7.3 Checking a run

- `logs/<job>-<jobid>.out` lists every embryo processed, its frame interval, refocus pauses, and kept/skipped cell counts with skip reasons.
- The stage 2 audit checks group membership, provenance (every table's `RunID` and `SourceFile` against the manifest), NPZ ↔ CSV agreement, the frequency axis against each embryo's Nyquist limit, and band-edge pile-up.
- Stage 5 recomputes the observed spectrum for the first three cells of every embryo and stops if it differs from `SpectralLibrary.npz`. `missing_cells.csv` should be empty.

---

## 8. Method parameters

Preprocessing follows the PDT study wherever that study defines a choice.

| Step | Setting |
|---|---|
| Track filter | ≥100 frames |
| Gap bridging | Linear interpolation of ≤3 missing frames (tracking max gap) |
| Refocus pauses | Split the trace into segments; drop the frame on each side of the pause |
| Baseline F₀ | Mean of the lowest 10% of intensities, per segment; ΔF/F = (F − F₀)/F₀ |
| Detrend | Linear, per segment, after resampling to a uniform time grid |
| Segment length | Must exceed 2× the cone of influence at the lowest frequency (≈490 s) |
| Wavelet | Complex Morlet `cmor1.5-1.0` (PyWavelets) |
| Frequency band | 0.005–0.08 Hz, 80 linear bins (periods 12.5–200 s) |
| Power | \|W\|²/scale, scale-rectified (Liu et al., 2007) |
| Cone of influence | Each segment transformed separately; samples within ≈1.22/f s of an edge masked |
| Cell retention | ≥1 COI-valid cycle at every frequency |
| Red-noise null | AR(1), α = (r₁ + √r₂)/2 (Torrence & Compo, 1998); 1,000 variance-matched surrogates through the identical CWT |
| Rhythmic call | Global max-ratio test, BH-FDR within embryo, q < 0.05 |
| Mapping | 9×9-px CNPUs, pooled CNPU mean, 49-CNPU smoothing over occupied CNPUs, z-score on a −1…1 jet scale (as in the PDT study) |
| Combined maps | CNPUs supported by ≥2 embryos |
| Statistics | PDT-study rule (two-sided K-S + Cohen's d, Bonferroni; significant = p_adj < 0.05 **and** \|d\| > 0.2), plus embryo-level tests alongside it |

**Exclusions.** `180405SP_E4` and `200114SP_E1` (S13) are excluded from mapping for the time being. Their boundary files come from a different recording than their traces (registration R² = 0.26 and 0.09), so I am working to determine where the right versions of those files may be. They remain in the spectral and red-noise tables, flagged `MappingExcluded`.

---

## 9. Validation and benchmarking

### 9.1 Correctness checks

| Check | Result |
|---|---|
| Synthetic sinusoids (3.3 s sampling, 3% dropped frames, a 30 s pause) | True frequencies 0.025–0.12 Hz recovered within one bin |
| Equal-amplitude sinusoids at different frequencies | Peak power roughly flat across frequency after scale rectification (before it, power grew with scale) |
| White noise (200 traces) | No low-frequency pile-up (0% at 0.02 Hz) |
| Stage 2 audit of the production run (job 594346) | **0 issues**; max (f_max × wavelet spread)/Nyquist = 0.99 (< 1 required) |
| Template-mappable cells kept | 191,518 of 205,446 (93.2%) |
| Join to PDT per-cell tables | 100% of 191,527 mapped cells have a PDT row |
| Rhythmic-fraction maps self-check | Reproduces the production `mapping_v3` maps exactly |
| Red-noise observed-spectrum check | Recomputed spectra equal `SpectralLibrary.npz` for every embryo |
The full audit trail is in `logs/calcium_iii_*.md`. The 2026-11-01 audit of the first (v1) pipeline found five problems: frequencies above Nyquist, scale-biased power, CWT across spliced segments, mis-grouped embryos, and two wrong boundary files. All five were fixed in v3 and re-audited.

### 9.2 Compute benchmarks (production v3 run)

Measured with `sacct` on Bora. "CPU efficiency" = total CPU time / (elapsed × CPUs requested).

| Stage | Jobs | Cells per job | Wall time per job | Peak memory | Requested | CPU efficiency |
|---|---|---|---|---|---|---|
| 0 Raw cleanup | 4 | — | 12–14 min | — | 1 CPU, 5 GB, 3 h | — |
| 1 Wavelet features | 4 | 127,977–159,424 tracked | 26–32 min | 10.0–11.9 GiB | 4 CPUs, 32 GB, 24 h | 24–25% |
| 2 Audit | 1 | 300,023 | 55 s | — | 1 CPU, 8 GB, 1 h | — |
| 4a Mapping | 4 | ~42,000–53,000 mapped | 4.9–5.1 min | 1.8–2.0 GiB | 1 CPU, 12 GB, 10 h | ~90% |
| 4b Combine | 1 | 191,527 | 55 s | 0.6 GiB | 1 CPU, 16 GB, 4 h | 80% |
| 5a Red-noise pilot | 1 | 50 | 3.4 min | — | 2 CPUs, 8 GB | — |
| 5a Red-noise test | 4 | 68,413–83,052 | 9.2–12.4 h | 17.1–19.6 GiB | 16 CPUs, 120 GB*, 48 h | 82–89% |
| 5c Rhythmic maps | 1 | 4 groups | 9.8 min | 1.9 GiB | 1 CPU, 16 GB, 10 h | 95% |
| 6b PDT ↔ wavelet panels | 1 | 192,411 joined | 2.4 min | 7.1 GiB | 2 CPUs, 64 GB, 4 h | 44% |

\* The red-noise jobs were submitted with `--mem=120G`; the script default is now 64 GB.

**Interpretation**

- **The red-noise test is ~99% of the compute.** It used about 603 CPU-hours (≈7.2 CPU-seconds per cell for 1,000 surrogates). All other stages together used about 3 CPU-hours. It is the only stage that needs, and benefits from, parallelism.
- **Its parallelism scales well.** Embryos are distributed over 16 worker processes, at 82–89% CPU efficiency. The lost ~15% is most likely tail time, when the last few large embryos run while other workers sit idle. Finer-grained work units (chunks of cells instead of whole embryos) would recover most of it.
- **The pilot underestimated runtime by about 2×.** The pilot embryo cost 3.8 s per cell, so S13 was projected at ~5 h. Across the full run, per-cell cost was 5–8 s, because it varies with trace length and the number of segments. The 48 h limit still left ~4× headroom. In future, pilot on a median-sized embryo or on several embryos.
- **Stage 1 is effectively single-threaded.** CPU time ≈ wall time, so 3 of the 4 requested CPUs sat idle. Requesting 1 CPU and 16 GB would free resources without slowing the run. Alternatively, process embryos in parallel to use all 4.
- **Several time limits are generous.** Mapping requests 10 h but takes 5 min, and the red-noise jobs peak at ~20 GiB of 64–120 GB. Tighter requests (e.g. mapping `--time=00:30:00 --mem=4G`, red-noise `--mem=32G`) would schedule faster on a busy cluster. The current values were kept as a safety margin for the production run.
- **Throughput.** With the four groups running in parallel, the production run takes about 13 h of wall time end to end, dominated by the red-noise step. Without it, the run takes under 1 h.

### 9.3 Method benchmark: wavelet features vs the 204 PDTs

Stage 6b correlates each cell's spectral features with its event count under every one of the 204 PDTs (179 embryos, 192,411 cells). This shows how the new method relates to the established one:
- **Agreement depends on the PDT, and mostly on whether it uses embryo-wise (EW) or cell-wise (CW) thresholds.** Band power tracks event counts *positively* for EW PDTs (median r = +0.50, change-in-magnitude set) and *negatively* for CW PDTs (median r = −0.52).
  - Interpretation: CW thresholds count small crossings in quiet, low-variance cells, which have flat, high-frequency spectra. This matches the PDT study's simulation result that CW and lower thresholds favour the opposite direction from EW and higher thresholds.
- **PDTs agree with each other only weakly** (median PDT-vs-PDT r = 0.22; 10th–90th percentile −0.08 to 0.65). So no single PDT can serve as "ground truth" for the wavelet features.
- Consequence: comparisons with events are reported as a distribution over all 204 PDTs, never for one PDT picked after the fact.

---

## 10. Results and reporting

Status labels: **Established** (validated, embryo-level) · **Preliminary** (embryo-level test done, but the full model is pending) ·

| Result | Status |
|---|---|
| **Cell flow:** 576,943 tracked → 300,023 passed spectral QC → 191,527 mapped (177 embryos) → 189,712 inside the CNPU grid. The main reason for skipping a cell is having no gap-free segment longer than ~490 s | Established |
| **Most activity is indistinguishable from red noise.** Cells exceeding AR(1) red noise: S13 1.4%, S13 devitellinated 10.6%, S14 1.4%, S14 devitellinated 7.2% (embryo medians 0 / 7.9 / 0 / 2.3%). Inspection shows these cells carry irregular, slow, large transients, not sustained oscillations | Preliminary |
| **Removing the vitelline membrane raises band power at both stages.** Embryo medians: 3.03× at stage 13 (46 vs 46 embryos, Mann–Whitney p = 2×10⁻⁹) and 2.32× at stage 14 (46 vs 41, p = 1.7×10⁻⁶). It holds within each recording year | Preliminary (stage × membrane model with a recording-year covariate pending) |
| **Event-vs-spectral agreement depends on the PDT** (section 9.3) | Preliminary |
| Anterior-high dominant-frequency pattern in all four groups' maps | Not interpreted: many spectra are noise-like, so frequency maps may track signal-to-noise rather than rhythm |

### Public-facing output

| Image panel | What it shows | Script |
|---|---|---|
| Fig 1 | One cell: ΔF/F trace, wavelet power with the cone of influence, time-averaged spectrum and its six features | `rednoise_bandpower_figures.py` |
| Fig 2A–C | Example spectra vs red noise; % of cells exceeding red noise per embryo; % per CNPU on the template | `rednoise_bandpower_figures.py` |
| Fig 3A–B | Band-power z-maps (2×2 groups); embryo medians of power, entropy and peak-power fraction | `rednoise_bandpower_figures.py` |
| PDT panels | All 204 PDTs vs wavelet features; membrane effects on events vs spectra | `pdt_wavelet_panels.py`, `pdt_wavelet_panels_styled.py` |

---

## 11. Milestones

| Date | Milestone | Evidence |
|---|---|---|
| 2026-10-26 | Repository created; MATLAB wavelet methods ported to Python/Slurm; spectral-table audit | Commit `5ff66d1`; `logs/calcium_iii_matlab_to_slurm_port_2026-10-26.md` |
| 2026-10-27 | PyWavelets port; extreme-power and gap audit | `logs/calcium_iii_pywavelets_port_2026-10-27.md` |
| 2026-10-28 → 30 | Stage 0 raw-file cleanup for all four groups | Slurm 573553, 577774–577776 |
| 2026-11-01 | Full pipeline audit: 5 methodological problems found | `logs/calcium_iii_pipeline_audit_2026-11-01.md` |
| 2026-11-02 | v3 fixes; production wavelet run, audit (0 issues), maps | Slurm 594314–594317, 594346, 596374–596378; `logs/calcium_iii_v2_pipeline_changes_2026-11-02.md` |
| 2026-11-03 | Audit of the PDT data: 100% join to mapped cells | `logs/calcium_iii_pdt_data_audit_2026-11-03.md` |
| 2026-11-03 → 04 | Red-noise test (pilot + 4 group jobs), FDR, rhythmic maps; SfN poster figures | Slurm 598279, 598337–598340, 600250; `logs/calcium_iii_rednoise_audit_2026-11-04.md` |
| 2026-11-05 | Environment pinned exactly (lock file + hashed pip requirements) | Commit `3cb8e08` |
| 2026-11-08 | PDT ↔ wavelet comparison panels | Slurm 608104 |

**Responding to feedback.** The pipeline was rebuilt after the 2026-11-01 audit (v1 → v3), not patched. The results interpretation was revised after the red-noise test: the earlier "quasi-periodic ~80–100 s component" was withdrawn because the data didn't support it. Each change is recorded in a dated note in `logs/`.

---

## 12. Known limitations

- **Data are not public.** The inputs live on the W&M cluster, so the repository is reproducible only with lab data access. Paths default to the original author's account (section 7.1).
- **Last-frame bias.** Only cells present on the last frame have template outlines and can be mapped.
- **Short tracks are lost.** 48% of tracked cells are dropped before spectral analysis, mostly because they have no gap-free segment long enough for the 0.005 Hz lower band edge.
- **Noise-like spectra.** Median spectral entropy is close to the flat-spectrum maximum, and in intact embryos ~42–44% of cells peak at the 0.08 Hz upper edge. Power and structure features are more trustworthy than frequency features.
- **The red-noise null is conservative.** "Not rhythmic" means "not detectable above AR(1)", not "absent".
- **Recording year is unbalanced across groups** (S13 mostly 2018; devitellinated groups mostly 2019–20). It must be a covariate in the group-effect models.
- **Two embryos are excluded from mapping** (section 8).

---

## 13. Citations

**Methods**
- Torrence, C., & Compo, G. P. (1998). A practical guide to wavelet analysis. *Bulletin of the American Meteorological Society*, 79(1), 61–78.
- Liu, Y., Liang, X. S., & Weisberg, R. H. (2007). Rectification of the bias in the wavelet power spectrum. *Journal of Atmospheric and Oceanic Technology*, 24(12), 2093–2102.
- Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate. *Journal of the Royal Statistical Society B*, 57(1), 289–300.
- Holm, S. (1979). A simple sequentially rejective multiple test procedure. *Scandinavian Journal of Statistics*, 6(2), 65–70.

**Software**
- Lee, G. R., Gommers, R., Waselewski, F., Wohlfahrt, K., & O'Leary, A. (2019). PyWavelets: A Python package for wavelet analysis. *Journal of Open Source Software*, 4(36), 1237.
- Harris, C. R., et al. (2020). Array programming with NumPy. *Nature*, 585, 357–362.
- Virtanen, P., et al. (2020). SciPy 1.0: fundamental algorithms for scientific computing in Python. *Nature Methods*, 17, 261–272.
- McKinney, W. (2010). Data structures for statistical computing in Python. *Proceedings of the 9th Python in Science Conference*, 56–61.
- Hunter, J. D. (2007). Matplotlib: A 2D graphics environment. *Computing in Science & Engineering*, 9(3), 90–95.
- Liberation Sans fonts (bundled in `scripts/step6_figures/fonts/`), SIL Open Font License.

**Data and prior analysis**
- Paudel, S. (2022). PhD dissertation, William & Mary.
- Paudel, S., et al. *In vivo analysis of calcium activity during early neural development* (manuscript in preparation).
