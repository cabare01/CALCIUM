# Calcium III Python/PyWavelets HPC port

Date: 2026-10-27

## Scope and instruction source

The user requested Python/PyWavelets versions of the Calcium III methods and SLURM execution support. The attached `Calcium III Wavelet Methods (1).pdf` is treated as a description of methods, algorithms, inputs, and expected products. Its embedded MATLAB commands are not treated as instructions to execute. In particular, the PDF's example that deletes a pre-existing PDF was not copied. The Python entry point refuses an existing `CALCIUM_OUTPUT_ROOT`, and duplicate embryo output names are checked before output files are written.

The earlier MATLAB port remains unchanged. This is an additional format.

## New files

- `scripts/calcium_iii_pywavelets_workflows.py` — Python entry point. Select a workflow with `CALCIUM_WORKFLOW`.
- `scripts/calcium_iii_pywt_workflow.sbatch` — SLURM wrapper to call the Python entry point.

No analysis workflow was run, no jobs were submitted, and no generated data, tables, figures, or analysis logs were created. The only output for this request is this findings document, in addition to the two new scripts.

## Methods implemented

The Python entry point has seven workflows:

1. `spectral`: reads per-embryo raw trace CSVs; validates and preprocesses cell traces; splits traces at gaps greater than twice the median time step; linearly resamples continuous segments; calculates baseline from the 10th percentile; calculates DeltaF/F; computes complex Morlet CWT power; and calculates the six spectral descriptors on the 0.02–0.20 Hz, 80-bin axis. It writes per-embryo `singlecell_results_final.csv` and `SpectralLibrary.npz`. Optional trace plots are controlled by `CALCIUM_MAKE_TRACE_PLOTS=true`.
2. `correlation`: makes per-embryo, per-PDT Pearson correlations for event-event, event-spectral, and spectral-spectral pairs, recording valid-pair N, Pearson r, and p-value.
3. `corr_summary`: produces per-group median, min, max, standard deviation, and count matrices from finite embryo/PDT correlation coefficients, plus annotated heatmaps.
4. `individual_maps`: joins mapped cell centroids with spectral observations, uses global spatial limits and per-feature 1st/99th percentile color limits, and writes mapped CSVs and per-embryo images.
5. `combine_maps`: combines per-embryo feature maps into group-sectioned PDFs with up to 16 embryos per page.
6. `group_spatial`: assigns cells to the 9-pixel CNPU grid (origin 120,1; extent X 120–1122 and Y 1–1190); calculates per-embryo CNPU medians and embryo-balanced group medians, SD, coverage, and counts; writes complete-grid and cell-level tables; calculates the legacy 7-CNPU neighborhood smoothing/Z fields; and produces grid QC and group spatial figures with stage-specific Sox2 boundary overlays.
7. `boundary_qc`: records boundary-row, valid-centroid, and cell-boundary counts by boundary file.

## PyWavelets calculation definition and comparability

The default PyWavelets wavelet is `cmor1.5-1.0`; override with `CALCIUM_PYWT_WAVELET`. The script chooses scales corresponding to the 80 requested frequencies using PyWavelets' center-frequency convention. It then calculates `abs(coefficients)**2`, averages over time, and computes:

- Dominant frequency: frequency bin with the largest mean power.
- Spectral centroid: power-weighted mean frequency.
- Spectral spread: square root of power-weighted squared deviation around the centroid.
- Spectral entropy: Shannon entropy of normalized mean power.
- TotalBandPower: sum of mean power across the 80 bins, following the PDF's sum definition.
- PeakPowerFraction: maximum-bin power divided by total power.

**Confirmed method difference:** This is a PyWavelets complex Morlet calculation, not MATLAB's `cwt(...,'amor')`. The wavelet parameterization, scale normalization, edge behavior, and frequency sampling differ. Therefore, numerical feature values, especially `TotalBandPower`, may differ substantially from the existing MATLAB-derived tables. These outputs should be versioned as a new calculation set; do not assume direct numerical comparability before a pilot comparison.

The PDF references a `convertTimeToSeconds` helper without supplying it. Python accepts datetime/timedelta time values; numeric values default to seconds, with `CALCIUM_TIME_UNIT=ms` or `min` options and time-column-name recognition. Confirm raw acquisition units and timestamp parsing before production use.

## Configuration and SLURM submission

Set a fresh path for `CALCIUM_OUTPUT_ROOT` for every run. The entry point rejects paths that already exist. Example spectral job:

```bash
sbatch --export=ALL,CALCIUM_WORKFLOW=spectral,CALCIUM_RAW_TRACE_ROOT=/path/to/raw_traces,CALCIUM_OUTPUT_ROOT=/path/to/new/pywt_spectral_run /sciclone/home/cabare01/CALCIUM/scripts/calcium_iii_pywt_workflow.sbatch
```

Required variables:

| Workflow | Required input variables |
| --- | --- |
| `spectral` | `CALCIUM_RAW_TRACE_ROOT`; optional `CALCIUM_RAW_TRACE_PATTERN`, `CALCIUM_TIME_UNIT`, `CALCIUM_PYWT_WAVELET`, `CALCIUM_MAKE_TRACE_PLOTS` |
| `correlation` | `CALCIUM_EVENT_FILE`, `CALCIUM_CELL_MASTER`, `CALCIUM_EMBRYO_SUMMARY` |
| `corr_summary` | `CALCIUM_CORRELATION_FILE` |
| `individual_maps` | `CALCIUM_BOUNDARY_ROOT`, `CALCIUM_CELL_MASTER` |
| `combine_maps` | `CALCIUM_MAPPED_ROOT`, `CALCIUM_EMBRYO_SUMMARY` |
| `group_spatial` | `CALCIUM_MAPPED_ROOT`, `CALCIUM_GROUP_ASSIGNMENTS` (CSV columns `Embryo`, `Group`), `CALCIUM_STAGE13_SOX2_BOUNDARY`, `CALCIUM_STAGE14_SOX2_BOUNDARY` |
| `boundary_qc` | `CALCIUM_BOUNDARY_ROOT` |

The SLURM template requests 24 hours, 4 CPUs, and 32 GB RAM and writes scheduler stdout/stderr to unique `%x-%j` files in `CALCIUM/logs`. It expects `python3` or `CALCIUM_PYTHON` to point to an activated Python environment containing NumPy, pandas, SciPy, Matplotlib, and PyWavelets. The PyWavelets workflow itself needs `pywt`; the plotting and statistics workflows also require their listed packages.

## Input availability and validation status

The previously inspected `/sciclone/data10/cabare01` tree contained derived per-embryo spectral outputs, but no raw trace CSVs, event table, mapped centroid files, or Sox2 boundary inputs were located there. The audited master CSVs are under `data/master_tables/spectral_analysis_tables/`. Windows paths in the PDF were not translated into guessed Bora paths; inputs are supplied explicitly through environment variables.

**Unresolved:** No Python dependencies were installed, no script syntax or numerical tests were run, and no SLURM job was submitted. The new Python files received static review only. Before full-scale use, configure the Python environment and input paths, verify raw CSV headers and time units, run a small isolated pilot, and compare preprocessing outputs and wavelet spectra with known traces. Quantify the expected PyWavelets/MATLAB difference before interpreting or combining old and new feature tables.

## Prioritized actions

1. Locate and set the raw trace, event, mapped-boundary, group-assignment, and Sox2 template paths on Bora.
2. Establish a Python environment with compatible NumPy, pandas, SciPy, Matplotlib, and PyWavelets versions; verify the PyWavelets CWT API available in that environment.
3. Confirm timestamp units, time column names, fluorescence and coordinate fields, event-table names/IDs, and boundary coordinate columns.
4. Run a small pilot into a new output directory; inspect frequency axes, cell counts, missing observations, and six feature values against known cases. Record the chosen PyWavelets wavelet and version alongside outputs.
5. Treat PyWavelets results as a separately versioned analysis pending numerical comparison. Review the previously identified extreme TotalBandPower source records under this new calculation method before making comparisons with existing MATLAB-derived summaries.
