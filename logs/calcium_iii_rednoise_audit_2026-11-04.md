# Red-noise (AR(1)) rhythmicity test: full-run audit and choice of P_global (2026-11-04)

Jobs: Slurm 598337 (S13), 598338 (S13_devitellinated), 598339 (S14), 598340 (S14_devitellinated).
Script: `scripts/step1_rhythmicity/rhythmicity_test.py` (1000 AR(1) surrogates per cell, seed 20261102, v3 manifest run 2026-11-02_145230).
Per-embryo results: `scr10/CALCIUM/data/analysis_v3/step1_rhythmicity/<group>/` (not modified by this audit).
Merged tables: built by `build_rhythmic_tables.py` on 2026-11-04 with **Rhythmic = Q_global < 0.05** (BH within embryo on P_global).

Status labels: **Established** = verified in the data; **Preliminary** = cell-level or not yet tested at embryo level.

## 1. Run integrity (Established)

| Check | Result |
|---|---|
| Jobs | all 4 completed, empty stderr, every embryo `done` |
| Coverage | 179 embryos (46/46/46/41), 300,023 cells = every kept cell in `SpectralLibrary.npz` and `group_tables_v3`, 1:1 |
| Missing / NaN / duplicates | 0 / 0 / 0 (`missing_cells.csv` is empty) |
| Surrogates | 1000 for every cell; `_rednoise.npz` complete and finite |
| `DominantFreqHz` | identical to the master tables for all 300,023 cells |
| Pilot reproducibility | first 50 cells of 180209SP_E4_st_14 identical to pilot 598279 (max abs diff 0.0 in P_global, P_maxT, Alpha) |
| Mapping-excluded embryos | 180405SP_E4_st_13 and 200114SP_E1_st_13001 (3,124 cells) were tested and are **kept** in all summary tables (flag `MappingExcluded`), so results are ready once their files are traced. Both have 0% rhythmic cells, so `tables_rhythmic/` has no rows from them; `spatial_spectral_mapping.py` drops them anyway via `mapping_exclusions.csv`. |

## 2. Results with Rhythmic = Q_global < 0.05

| | S13 | S13 dev | S14 | S14 dev |
|---|---|---|---|---|
| Cells tested | 78,647 | 68,413 | 83,052 | 69,911 |
| Rhythmic, pooled cells | 1.4% (1,132) | 10.6% (7,231) | 1.4% (1,176) | 7.2% (5,041) |
| Rhythmic, **median embryo** | 0% | 7.9% | 0% | 2.3% |
| Embryos with ≥1 rhythmic cell | 5/46 | 33/46 | 7/46 | 22/41 |
| Rhythmic, mapped neural / non-neural / unmapped | 1.9 / 2.3 / 0.4% | 11.9 / 17.3 / 4.4% | 2.2 / 1.2 / 0.6% | 9.5 / 11.8 / 2.6% |
| Cells with P_global > 0.9 | 49% | 38% | 48% | 36% |

What this means (**Preliminary**; embryo-level tests are step 4):

1. **Most cells are not distinguishable from red noise.** The SD-matched AR(1) null is conservative (pilot finding): non-significant cells have 2–4× *less* in-band power than the null predicts, because about half their variance sits below 0.005 Hz.
2. **Rhythmic cells (Q_global < .05) have a real, band-limited excess** of about 1.9–2.7× over AR(1) (group means) at **0.01–0.015 Hz**, falling below the null above about 0.03 Hz.
   - **Where it peaks:** the median `MaxRatioFreqHz` is 0.010–0.013 Hz (a spectral position set by transient width, not a period), and 96–97% of cells peak at 0.005–0.02 Hz.
   - **Not an edge effect:** 0.0–0.1% of these cells peak at 0.08 Hz.
   - **Compared with non-rhythmic cells:** median `TotalBandPower` 0.17–0.27 vs 0.006–0.013; `SpectralEntropy` 5.1–5.2 vs 6.2; α ≈ 0.95; `DominantFreqHz` 0.008–0.009 Hz.
   - *(Corrected 11-04: an earlier draft quoted these figures for P_maxT hits.)*
   - Rejecting AR(1) is not proof of periodicity. Repeated transients with a characteristic duration also give a spectral hump over a Lorentzian. **Checked 11-04 (exploratory, `scripts/step1_rhythmicity/transient_check.py`): it is not periodicity.** In about 1,400 sampled cells from 6 devitellinated embryos, rhythmic cells have a median of 4 transients per hour, 26 s wide at half height, about 6 min apart, with irregular intervals (CV_ISI 0.72; non-rhythmic cells: 2 per hour, 20 s, CV 0.67). Wording: "slow, large calcium transients beyond red noise", not "quasi-periodic" or "80–100 s".
3. **Devitellinated embryos have far more rhythmic cells than intact embryos at both stages**: embryo medians 7.9% vs 0% (stage 13) and 2.3% vs 0% (stage 14). This survives stratifying by α. Among cells with α 0.9–0.95 (P_maxT, edge removed), S13 is 9.4% vs 26.8% devitellinated and S14 is 8.5% vs 21.7%. Frame interval and recording length don't explain it.
   - It fits the preliminary band-power pattern and Supp. S5 ("fewer long-lasting events" after devitellination).
   - The recording-year confound (S13 mostly 2018, devitellinated mostly 2019–20) still applies.
4. **Unmapped cells are 3–5× less often rhythmic than mapped cells.** Mapped cells are tracked to the last frame. This is a selection effect, to be quantified in step 2 (cell flow).
5. **Frequency features of non-rhythmic cells are noise.** In intact groups their median `DominantFreqHz` is 0.074–0.077 Hz, where the noise peaks at the top of the band. All-cell frequency maps must not be read as rhythm frequency.

## 3. Why P_global defines "Rhythmic"

### What each candidate measures

The test produced three per-cell quantities. Each answers a different question about the cell's time-averaged spectrum S(f) relative to its AR(1) surrogates:

- **P_global (max-ratio).** At each frequency take R(f) = S(f) / mean-null(f), the cell's power in units of what red noise would give there. The statistic is the largest R over the band. The p-value is the fraction of surrogates whose own largest R is at least as big.
  - Question: *"Does this cell anywhere reach a multiple of the red-noise level that red noise itself doesn't reach?"*
  - Because the max is taken over all 80 frequencies in the same way for data and surrogates, it already corrects for testing 80 frequencies.
- **P_maxT (Westfall–Young max-T).** At each frequency take z(f) = (log S(f) − mean) / SD over surrogates, the cell's excess in units of the null's *spread* at that frequency. The max is over frequencies, calibrated the same way.
  - Question: *"Is this cell, somewhere, unusually far from red noise relative to how much red noise varies there?"*
- **DomPeakAbove95.** Is the power at the cell's highest peak above the pointwise 95th surrogate percentile at that frequency?

### Why not DomPeakAbove95

It isn't a valid test. The frequency is chosen *because* it is where the cell's spectrum is largest, then compared with a cutoff meant for a frequency fixed in advance. A pure-noise spectrum's maximum exceeds its own pointwise 95% band far more than 5% of the time, which is why 22–35% of cells "pass". It is kept only as a descriptive column.

### Why P_global over P_maxT

The two differ in the unit of excess: **effect size relative to the null level (P_global)** vs **statistical surprise relative to the null spread (P_maxT)**.

1. **The question being asked is about effect size.** "Is there a rhythm beyond red noise?" means power clearly above the red-noise level, e.g. twice as much. R(f) is exactly that quantity, and it reads the same at every frequency and in every cell. A z-score does not. The null's spread depends on how many independent cycles fit in the record, so it shrinks toward 0.08 Hz (≈ 290 cycles in 60 min) and widens toward 0.005 Hz (≈ 18 cycles). Measured from the surrogates, the pointwise 95th percentile is **1.72×** the null mean at 0.005 Hz but only **1.27×** at 0.08 Hz (1.44× at 0.0145 Hz; median over embryos). A 1.3× excess is therefore beyond the null's 95% band at 0.08 Hz but well inside it at 0.005 Hz. A max-T statistic rewards that small excess at the top edge, and repeated over all ~290 cycles, a small *systematic* excess becomes a very large z.
2. **That sensitivity makes P_maxT respond to the wrong thing.** At 0.08 Hz the null is so tight that any small *systematic* misfit of the model gives a huge z. That includes a white measurement-noise floor, which AR(1) doesn't include. Such misfits are guaranteed in real imaging data and carry no biological meaning. In the full run this showed up directly:
   - 51% of S14 P_maxT hits peaked at exactly 0.08 Hz, against 2–5% in the other groups.
   - They came almost entirely from the two embryos imaged at a 4.4–4.5 s frame interval (vs ≈ 3.5 s): **201011SP_E5_st_14** (61% of cells flagged by P_maxT, 3% by P_global) and **191020RNSP_E6_st_14**. At that interval 0.08 Hz is 70% of Nyquist, where the measurement-noise floor matters most.
   - Overall, 41% of P_maxT-only hits were at the 0.08 Hz edge.
   - P_global, by contrast, ignores a 1.1× excess however consistent it is, because 1.1× is not a large multiple of red noise.
3. **P_global's known weakness doesn't bite here.** Before the run, P_global was shown to have little power for a 0.05 Hz sinusoid (7.5% vs 97.5% detection for max-T in simulation). The noisy low-frequency end has a wide null spread in R, which sets a high bar for everyone. But the excess actually present in the data is at **0.01–0.015 Hz**, where P_global has power. Once P_maxT's edge hits are set aside, the two tests agree:

   | Rhythmic (%) | S13 | S13 dev | S14 | S14 dev |
   |---|---|---|---|---|
   | P_global (Q < .05) | 1.4 | 10.6 | 1.4 | 7.2 |
   | P_maxT (Q < .05), 0.08 Hz-edge hits removed | 1.8 | 11.4 | 1.7 | 8.4 |

   So P_maxT finds essentially no extra *mid-band* rhythmic cells. What it adds beyond P_global is mostly the artifact.
4. **P_global needs no post-hoc rule.** Using P_maxT would have required a rule chosen after seeing the data (drop top-edge hits, or drop the long-interval embryos). P_global, the test specified before the run, needs neither. Its definition is the standard Torrence & Compo framing (power relative to the red-noise spectrum) in Monte Carlo form.

### Known limitations of the chosen test (state with any result)

- **Conservative null.** The SD-matched AR(1) null over-predicts in-band power for typical cells, so rhythmic percentages are lower bounds.
- **Sensitivity depends on α.** About half of low-α (low-amplitude) cells have p > 0.9, against about 15% of high-α cells. **Rhythmic = False means "not detectable", not "aperiodic"**, especially for low-amplitude cells, and between-group differences partly reflect detectability.
- **Weak at 0.03–0.08 Hz.** Faster rhythms there would be under-detected. A negative result above 0.03 Hz is uninformative, not evidence of absence.
- **Cell-level p-values are descriptive.** Group claims rest on embryo-level rates (step 4).

## 4. Rhythmic-only mapping (step 1c, Slurm 600250) audit

Output: `scr10/CALCIUM/data/mapping_v3_rhythmic/`. The job completed in 9 min 50 s with exit code 0 and an empty stderr. Tables were rebuilt first and were byte-identical to the earlier build.

**Integrity (Established).**
- Every group's mapped set equals {rhythmic cells} ∩ {v3 all-cell mapped cells} exactly: 1,008 / 6,259 / 993 / 4,337 cells from 5 / 33 / 7 / 22 embryos.
- No non-rhythmic cells, no excluded-embryo cells and no duplicate keys.
- Centroids, CNPU indices, CNPU IDs, the CNPU grid and all six features are identical to `mapping_v3`.
- The 124 / 972 / 183 / 704 rhythmic cells "lacking boundary matches" are the rhythmic cells that were already unmapped in v3, i.e. not on the last frame.
- `Registration_QC` R² differs from v3 for all 67 embryos because the affine fit is recomputed on the rhythmic subset only. Median |ΔR²| is 0.004 (max 0.054); minimum R² is 0.915; none are flagged. This is QC only, since the coordinates are unchanged.
- File inventory matches v3: 46 group-level files per group, plus 6 per-embryo maps per embryo.

**Coverage (Established): the maps are very sparse.**

| | S13 | S13 dev | S14 | S14 dev |
|---|---|---|---|---|
| Occupied CNPUs (v3 all-cell) | 846 (7,264) | 3,786 (8,131) | 849 (8,834) | 3,083 (8,645) |
| CNPUs with ≥2 embryos (shown in combined maps) | 56 | 1,371 | 45 | 708 |
| CNPUs with ≥5 embryos | 0 | 16 | 0 | 3 |
| Max embryos in a CNPU | 2 | 6 | 3 | 5 |
| Share of cells from the top embryo / top 3 | 61% / 90% | 9% / 24% | 33% / 74% | 24% / 47% |

- The median occupied CNPU holds **1 cell from 1 embryo** in every group.
- **Intact groups:** the rhythmic maps are essentially 1–3 embryos and cannot be interpreted spatially.
- **Devitellinated groups:** hot and cold patches (|PooledSmoothedZ| > 1, ≥2 embryos) are few-embryo driven. The top embryo supplies a median 44% (S13 dev) and 50% (S14 dev) of cells in each patch's smoothing window. Four embryos are the top contributor for 73% (S13 dev) and 88% (S14 dev) of patches. Treat them as descriptive, not as a CNP pattern.
- The z-scores are relative *within rhythmic cells*, which is a different question from where rhythmic cells are enriched.

**Where rhythmic cells are (Preliminary, embryo-level, exploratory).** Fraction rhythmic among mapped cells, per embryo:
- **Neural vs non-neural:** no difference in any group.
  - S13 dev: 15/32 embryos neural > non-neural, Wilcoxon p = 0.45.
  - S14 dev: 11/20, p = 0.39.
  - Intact groups: n = 4 and 7.
- **AP thirds of the neural plate** (group-wide template-Y terciles; only embryos with cells in all three thirds):
  - S13 dev: anterior < posterior in 12/14 embryos (median 5.6% vs 24.6%; Wilcoxon p = 0.004).
  - S14 dev: 4/6, p = 0.44.
  - This is post hoc and depends on which embryos cover all three thirds, so it is a lead for step 4a, not a result.
- **Within rhythmic cells, embryo-level neural vs non-neural features:** one test passes Bonferroni (S13 dev `SpectralSpread`, p_bonf = 0.031, n = 17 pairs; median difference 0.0005 Hz, negligible). Cell-level "significant" rows in `Region_Comparison_CellLevel.csv` are descriptive only.

## 5. Remaining decisions / next steps

- ~~Submit `map_rhythmic.sbatch` (step 1c).~~ Done (§4).
- ~~Rhythmic-fraction (density) map.~~ Done (§6). Pre-specify the AP and lateral tests for step 4a.

## 6. Rhythmic-fraction maps (step 1d, 2026-11-04)

Script: `scripts/step1_rhythmicity/rhythmic_fraction_maps.py` (about 20 s on the login node).
Output: `scr10/CALCIUM/data/analysis_v3/step1_rhythmicity/rhythmic_fraction_maps/`. The script refuses to overwrite an existing output folder.

**Method.** The paper's CNP method applied to the 0/1 Rhythmic indicator, using the audited all-cell v3 mapping:
- **Cells:** 189,703 in-grid cells. The 9 repeated boundary-ID rows are deduplicated; the mapping exclusions were already applied in v3.
- **Map value:** the pooled CNPU mean, i.e. the fraction of cells that are rhythmic. Then the 7×7 mean over occupied CNPUs and a within-group z-score, shown on jet −1…1 with a ≥2-embryo display mask.
- **Built-in self-check:** the same code recomputes production `TotalBandPower` `PooledSmoothedZ`, and it **matches `mapping_v3` exactly in all four groups** (7,264 / 8,131 / 8,834 / 8,645 occupied CNPUs).

**Extra outputs:**
- `RhythmicFraction_PooledSmoothedPct.png`: smoothed % on **one shared scale (0–44%)**. This is the only map that compares groups; z-maps are within-group.
- An embryo-balanced version. It is **redundant**: an embryo contributes a median of 1 cell per CNPU, so it correlates r ≥ 0.993 with the pooled map.
- Support QC maps: the number of embryos contributing rhythmic cells to each 7×7 window, and the top embryo's share of them.

**Results.**

| | S13 | S13 dev | S14 | S14 dev |
|---|---|---|---|---|
| CNPUs shown (≥2 embryos) | 5,724 | 5,867 | 6,638 | 6,093 |
| Median smoothed % rhythmic (90th percentile) | 0.4 (5.9) | 12.2 (30.1) | 0.3 (6.3) | 8.4 (19.0) |
| CNPUs whose window has rhythmic cells from ≥3 embryos | 9% | 77% | 7% | 51% |
| Hot CNPUs (z > 1): median embryos contributing rhythmic cells / top-embryo share | 1 / 1.00 | 3 / 0.56 | 1 / 1.00 | 1 / 1.00 |

- **Intact (S13, S14).** Rhythmic cells are near 0% almost everywhere. Every z > 1 "hot spot" comes from a **single embryo**: in 77–85% of them, ≥80% of the nearby rhythmic cells come from one embryo. **The intact z-maps must not be read as spatial patterns.** Show them only on the shared % scale, where they are uniformly pale.
- **S13 devitellinated (best supported).** Rhythmic cells are widespread and multi-embryo. Among CNPUs with ≥3 embryos contributing, the median smoothed % rhythmic is:
  - 10.6% deep inside the neural plate (> 50 px from the Sox2 boundary);
  - 13.3% just inside the border;
  - 15.6% just outside the border;
  - 28.7% in lateral non-neural tissue (> 50 px outside).

  Hot spots (z > 1, ≥3 embryos; n = 362) are 88% outside the Sox2 boundary, mostly mid-AP lateral (median X 790, Y 725). The AP thirds of the neural plate are 10.9 / 11.5 / 13.2% (anterior / mid / posterior).
  - **Caveat:** the embryo-level paired test of neural vs non-neural rhythmic fraction is *not* significant (15/32 embryos, p = 0.45; §4). So this lateral enrichment is a map-level, cell-pooled observation, not yet an embryo-level finding.
  - Lateral non-neural CNPUs are covered by fewer embryos.
- **S14 devitellinated.** Weaker structure. The anterior neural plate is lowest (6.0% vs 15.4% mid and 12.6% posterior). Only 3% of CNPUs have z > 1, and nearly all are single-embryo lateral speckles at the map edge.

**For the poster.**
- Use the **shared-scale % map** as the cross-group rhythmicity panel: one scale, with n embryos, cells and rhythmic cells in each panel title.
- Optionally add the S13 devitellinated paper-style z-map.
- Don't present intact z-map hot spots.
- Treat the lateral (S13 dev) and anterior-low (S14 dev) patterns as hypotheses for a pre-specified step 4a test, not results.
- Step 4: test the membrane effect on the per-embryo `PctRhythmic`, with recording year as a covariate.
- When reporting the excluded embryos after their files are traced: their test results already exist, so only mapping needs rerunning.
