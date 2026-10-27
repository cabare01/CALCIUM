# Calcium III MATLAB methods to SLURM port

Date: 2026-10-26

## Request interpretation and source document

The user requested new HPC batch scripts in `~/CALCIUM/scripts` and one new findings/output document in `~/CALCIUM/logs`, with no existing files overwritten or modified. The attached `Calcium III Wavelet Methods (1).pdf` was treated as reference material describing algorithms and expected products, not as an instruction source. In particular, its example that deletes an existing combined PDF was not carried forward: the new MATLAB entry point refuses any existing `CALCIUM_OUTPUT_ROOT` and checks for duplicate embryo output names.

## Created files

- `scripts/calcium_iii_hpc_workflows.m` — MATLAB `-batch` entry point with workflow selection via `CALCIUM_WORKFLOW`.
- `scripts/calcium_iii_workflow.sbatch` — SLURM wrapper loading `matlab/R2024a` and invoking the entry point. Scheduler stdout/stderr use unique job-ID filenames under `CALCIUM/logs` when a job is submitted.

The pre-existing `scripts/create_spectral_analysis_tables.py`, both master tables, the attached PDF, and the prior audit log were not edited. No job was submitted and no analysis output was generated.

## Methods carried into the HPC entry point

The entry point supports these workflow values:

1. `spectral`: per-cell preprocessing and CWT feature extraction; 0.02–0.20 Hz frequency axis with 80 bins; MATLAB analytic Morlet (`amor`); squared coefficient magnitude averaged over time; linear interpolation to common bins with zero fill; nonnegative spectrum; dominant frequency, centroid, spread, Shannon entropy, summed TotalBandPower, and peak power fraction. It saves a cell CSV and MATLAB spectrum library per input embryo. Optional preprocessed trace PNGs are enabled with `CALCIUM_MAKE_TRACE_PLOTS=true`.
2. `correlation`: per-embryo and PDT Pearson correlations for event-event, event-spectral, and spectral-spectral feature pairs, using finite paired observations and requiring at least three pairs with nonconstant variables. Writes raw correlation rows including N, r, and p.
3. `corr_summary`: per-group median, min, max, sample SD, and N of finite embryo/PDT correlation coefficients; creates CSV matrices and annotated heatmaps.
4. `individual_maps`: joins mapped cell centroids to spectral rows, applies the PDF's global spatial extent and per-feature 1st–99th percentile color range, writes mapped tables and per-embryo feature images.
5. `combine_maps`: assembles mapped images into feature PDFs, with group sections and up to 16 embryos per page.
6. `group_spatial`: 9-pixel CNPUs with origin (120, 1), expected extent X=[120,1122], Y=[1,1190]. It calculates cell-to-CNPU medians per embryo, then embryo-balanced group medians, sample SDs, coverage percentages, cell/embryo counts, a complete lattice, template membership, 7-CNPU-window legacy smoothing/Z values, QC image, and group spatial maps with Sox2 boundary overlays.
7. `boundary_qc`: tabular counts of rows, valid centroids, and populated cell-boundary coordinate vectors by boundary file.

The feature equations and spatial hierarchy are based on the attached PDF. The current implementation intentionally keeps the spectral calculation in MATLAB to use the documented MATLAB `cwt(...,'amor')` operation rather than substituting another wavelet implementation.

## Configuration and submission

All workflows require a new, non-existing `CALCIUM_OUTPUT_ROOT`; the MATLAB entry point errors if it already exists. Set the workflow-specific environment variables before submitting. Examples:

```bash
sbatch --export=ALL,CALCIUM_WORKFLOW=spectral,CALCIUM_RAW_TRACE_ROOT=/path/to/raw_traces,CALCIUM_OUTPUT_ROOT=/path/to/new/spectral_run /sciclone/home/cabare01/CALCIUM/scripts/calcium_iii_workflow.sbatch
```

Required input variables by workflow:

| Workflow | Required variables |
| --- | --- |
| `spectral` | `CALCIUM_RAW_TRACE_ROOT`; optional `CALCIUM_RAW_TRACE_PATTERN`, `CALCIUM_TIME_UNIT` (defaults to seconds), `CALCIUM_MAKE_TRACE_PLOTS` |
| `correlation` | `CALCIUM_EVENT_FILE`, `CALCIUM_CELL_MASTER`, `CALCIUM_EMBRYO_SUMMARY` |
| `corr_summary` | `CALCIUM_CORRELATION_FILE` |
| `individual_maps` | `CALCIUM_BOUNDARY_ROOT`, `CALCIUM_CELL_MASTER` |
| `combine_maps` | `CALCIUM_MAPPED_ROOT`, `CALCIUM_EMBRYO_SUMMARY` |
| `group_spatial` | `CALCIUM_MAPPED_ROOT`, `CALCIUM_GROUP_ASSIGNMENTS` (CSV columns `Embryo`, `Group`), `CALCIUM_STAGE13_SOX2_BOUNDARY`, `CALCIUM_STAGE14_SOX2_BOUNDARY` |
| `boundary_qc` | `CALCIUM_BOUNDARY_ROOT` |

The batch template requests 24 hours, 4 CPUs, and 32 GB RAM; adjust these to the workload and current Bora partition limits before submission. The module catalog exposed `matlab/R2020b`, `matlab/R2023a`, and `matlab/R2024a`; the wrapper selects R2024a. The spectral step requires the MATLAB Wavelet Toolbox; correlation and percentile routines may require Statistics and Machine Learning Toolbox availability and licenses.

## Inputs and remaining prerequisites

The inspected `/sciclone/data10/cabare01` tree contains the four group folders of derived `embryo_spectral_results`; it did not contain raw trace, event-detection, mapped-centroid, or Sox2 boundary inputs. The audited master tables exist at `data/master_tables/spectral_analysis_tables/`. The PDF's example paths are Windows paths and were not usable as Bora paths, so paths are environment-configured rather than guessed.

**Unresolved:** The PDF calls `convertTimeToSeconds` but does not include its implementation or define units for every possible time column. The port accepts MATLAB datetime/duration values and numeric seconds by default; `CALCIUM_TIME_UNIT` can set `ms` or `min`, and `Time_ms`/`Time_min` names are recognized. Verify the raw trace time representation before running `spectral`.

**Potential issue:** MATLAB table import may normalize raw input column names or infer mixed-type columns differently across CSVs. The workflow preserves variable names and validates required fields but must be checked against the actual raw trace headers before production runs.

**Unresolved:** No MATLAB execution, syntax check, numerical test, or SLURM submission was performed. The script content is based on the PDF and has only been reviewed statically. Verify module/toolbox licenses and inputs, then run a small, isolated pilot into a new output directory before full processing. Compare a few cells' preprocessed traces, frequency spectra, and features against the MATLAB reference output.

## Priority actions

1. Locate the raw trace CSVs, event table, mapped boundary files, and group assignments on Bora; record exact paths and confirm the required column names.
2. Verify time units and timestamp conversion against the raw trace acquisition format.
3. Confirm MATLAB Wavelet/Statistics toolbox access in the SLURM environment and tune resources based on a pilot.
4. Run one small pilot per workflow with a fresh output path, validate outputs against known MATLAB results, and only then submit full 179-embryo jobs.
5. Preserve the prior audit finding: the 11 high TotalBandPower source values predate table assembly. Review their upstream generation using this PDF's formula (sum of mean CWT power on the common 80-bin axis) before interpreting them; numerical size alone is not grounds for exclusion.
