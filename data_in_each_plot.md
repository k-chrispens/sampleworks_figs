# Data in Each Plot

This document describes, for every figure produced by the analysis notebooks, exactly what data goes in, how it is filtered and transformed, and what is left out.

---

## Global Filters (applied at data-loading time)

These filters are baked into the data-loading functions in `analysis_helpers.py` and apply to **every** downstream plot that uses that loader. Any data excluded here never reaches any figure.

### 1. Occupancy filter — path-tagged CSVs (`load_and_tag_path_csv` / `tag_path_dataframe`)

Used by: clashscore, bond angle, and bond length data.

- The file path of each row is parsed with a regex to extract `occ_spec` (e.g. `0.5occA_0.5occB`).
- **Only rows where `occ_spec == "0.5occA_0.5occB"` are kept.**
- All other occupancy configurations (e.g. `0.25occA_0.75occB`, `1.0occA`) are discarded.

### 2. Occupancy filter — RSCC data (`load_rscc`)

- The RSCC CSV has a string column `altloc_occupancies`.
- **Only rows where `altloc_occupancies == "{'A': 0.5, 'B': 0.5}"` are kept.**
- All other occupancy splits are discarded.

### 3. LDDT data (`load_lddt`)

- The LDDT CSV also has a string column `altloc_occupancies`; **only rows where `altloc_occupancies == "{'A': 0.5, 'B': 0.5}"` are kept** (same filter as RSCC).
- Columns retained: `protein`, `model_name`, `guidance_weight`, `selection`, `occupancies`, `avg_silhouette`, `avg_silhouette_to_ref`.

### 4. Model name aliasing (all loaders)

After loading, model names are canonicalized via `MODEL_NAME_ALIAS`:

| Raw name in CSV | Canonical name used in plots |
|---|---|
| `boltz2_X-RAY_DIFFRACTION` | `boltz2X-RAY` |
| `boltz2_MD` | `boltz2MD` |
| `boltz2xray` | `boltz2X-RAY` |
| `boltz2md` | `boltz2MD` |
| `rf3` (raw) | `rf3-old` |
| `rf3-fixed` | `rf3` |

For Boltz2 models in the RSCC and LDDT loaders, the model name is constructed by concatenating the `model` and `method` columns (e.g. `boltz2` + `X-RAY_DIFFRACTION`), then aliased.

### 5. Guidance strength bucketing (all loaders)

Every loader adds a `guidance_strength` column mapped from the numeric `guidance_weight`:

| guidance_weight | guidance_strength |
|---|---|
| 0.0 | baseline |
| 0.005 | weak |
| 0.01 | weak |
| 0.02 | medium |
| 0.035 | medium |
| 0.05 | strong |
| 0.1 | weak |
| 0.2 | medium |
| 0.35 | medium |
| 0.5 | strong |

### 6. Selection name mapping

When a `selection` column is present, a `selection_name` column is derived: if the selection string contains `"chain_id"`, it is labeled `"all altloc"`; otherwise the raw selection string is used as-is.

### 7. Model groups

Several figures are rendered in three variants via a UI dropdown. The model groups are:

- **Main models**: `boltz2X-RAY`, `boltz2MD`, `protenix`, `rf3` (i.e. everything NOT in `MODELS_ALT`)
- **MODELS_ALT** (alternative/variant models): `protenix-mini`, `protenix-tiny`, `rf3-old`, `rf3-chiral-off`
- **Protenix group**: `protenix`, `protenix-mini`, `protenix-tiny`
- **RF3 group**: `rf3`, `rf3-old`, `rf3-chiral-off`

---

## Notebook 1: Geometry Quality (`01_geometry.py`)


### Source data

| File | Metric columns | Rows (raw) |
|---|---|---|
| `initial40results/clashscore_metrics.csv` | `clashscore`, `num_clashes` | ~75k |
| `initial40results/bond_angle_violation_fractions.csv` | `outlier_fraction` (renamed `outlier_fraction_bond_angles`) | ~75k |
| `initial40results/bond_length_violation_fractions.csv` | `outlier_fraction` (renamed `outlier_fraction_bond_lengths`) | ~75k |

Clashscores are loaded via `load_and_tag_path_csv` with a hardcoded `selection = "all residues"` added to every row.

Bond angle and bond length CSVs are loaded separately, then merged on `["model", "model_n"]` to produce a single `bond_stats` dataframe with both `outlier_fraction_bond_angles` and `outlier_fraction_bond_lengths`. This merged frame is then tagged via `tag_path_dataframe`.

Both loaders apply the global occupancy filter (`occ_spec == "0.5occA_0.5occB"`).

---

### Figure 01 — Clashscore Distribution

**What is plotted:** Ridge plot (KDE per row) of raw clashscore values. One row per model, hue = guidance strength.

**Data included:**
- All clashscore rows that survive the global occupancy filter.
- Filtered to `guidance_strength != 'strong'` (baseline, weak, and medium only).
- Filtered to `selection_name != 'all altloc'`. (In practice this is a no-op because clashscores are loaded with a hardcoded `selection = "all residues"`, so `selection_name` is always `"all residues"` and never `"all altloc"`.)
- Further filtered by the active model group dropdown (main / protenix / rf3).

**What is excluded:**
- All occupancy configurations other than 50/50.
- **Strong guidance** (`guidance_weight` 0.05 or 0.5) is excluded.
- Models outside the selected group are excluded.

---

### Figure 02 — Delta Clashscore Distribution

**What is plotted:** Ridge plot of the *change* in clashscore relative to the no-guidance baseline, per model.

**Delta computation** (`make_delta_frame`):
1. Baseline rows (`guidance_strength == "baseline"`, i.e. `guidance_weight == 0.0`) are extracted per `(protein, model_name)`.
2. Non-baseline rows are merged onto their corresponding baseline.
3. `delta_clashscore = clashscore_guided - clashscore_baseline` is computed.
4. Baseline rows themselves are excluded from the result (only guided rows remain).

**Data included:**
- All (protein, model, guidance_weight) combinations that have both a baseline and a guided measurement.
- Filtered by model group dropdown.

**What is excluded:**
- All occupancy configurations other than 50/50 (global filter).
- **Rows with `delta_clashscore >= 100`** are excluded as outliers for readability.
- Models outside the selected group.
- Any (protein, model) pair missing a baseline measurement has no delta and is absent.

---

### Figure 03 — Delta Bond Angle Outlier Fraction

**What is plotted:** Ridge plot of the change in bond-angle outlier fraction relative to baseline, per model.

**Delta computation:** Same `make_delta_frame` logic as Figure 02, but on the `outlier_fraction_bond_angles` column.

**Data included:**
- All (protein, model, guidance_weight) combinations from the merged bond stats that have a baseline.
- Filtered by model group dropdown.

**What is excluded:**
- All occupancy configurations other than 50/50 (global filter).
- **Rows with `delta_outlier_fraction_bond_angles >= 0.05`** are excluded as extreme geometry distortions.
- Models outside the selected group.
- Any (protein, model) pair missing a baseline measurement.

---

## Notebook 2: Fit Quality (`02_fit_quality.py`)

### Source data

| File | Key columns |
|---|---|
| `initial40results/rscc_results.csv` | `protein`, `model`, `method`, `altloc_occupancies`, `guidance_weight`, `rscc`, `selection`, `ensemble_size` |
| `initial40results/lddt_results.csv` | `protein`, `model`, `method`, `altloc_occupancies`, `guidance_weight`, `selection`, `occupancies`, `avg_silhouette`, `avg_silhouette_to_ref` |

**RSCC aggregation step:** After loading, RSCC rows are grouped by `(protein, model_name, guidance_weight)` and averaged across all per-altloc-region selections (`.mean()`), yielding one RSCC value per (protein, model, gw) triple.

**LDDT baseline merge:** Baseline LDDT rows (`guidance_strength == "baseline"`) are extracted with their `avg_silhouette_to_ref` value, then merged back onto non-baseline rows by `(protein, model_name, selection_name)` to compute `delta_silhouette_to_ref = guided - baseline`.

---

### Figure 04 — RSCC Baseline Distribution

**What is plotted:** Combined histogram + KDE of RSCC values with all models overlaid on a single axis (one color per model). Baseline only.

**Data included:**
- Aggregated RSCC data (one row per protein/model/gw).
- Filtered to `guidance_strength == 'baseline'` only (gw = 0.0).
- Filtered by model group dropdown.

**What is excluded:**
- All occupancy configurations other than 50/50 (global RSCC filter).
- All guided rows (any `guidance_weight > 0.0`).
- Models outside the selected group.

---

### Figure 05 — RSCC Baseline vs Weak Guidance

**What is plotted:** Ridge plot (one row per model, hue = guidance strength) comparing RSCC distributions at baseline and weak guidance.

**Data included:**
- Aggregated RSCC data.
- Filtered to `guidance_strength in ('baseline', 'weak')`.
- Filtered by model group dropdown.

**What is excluded:**
- All occupancy configurations other than 50/50.
- **Medium and strong guidance** are excluded.
- Models outside the selected group.

---

### Figure 06 — Delta RSCC Distribution

**What is plotted:** Ridge plot of the change in RSCC from baseline, across all non-baseline guidance strengths, per model.

**Delta computation:** Same `make_delta_frame` logic — baseline RSCC is subtracted from each guided RSCC per (protein, model).

**Data included:**
- All non-baseline (protein, model, guidance_weight) triples that have a matching baseline.
- Filtered by model group dropdown.

**What is excluded:**
- All occupancy configurations other than 50/50.
- **Rows where `delta_rscc` falls outside (-0.3, +0.3)** are clipped/excluded for readability.
- Models outside the selected group.

---

### Figure 07 — Delta RSCC (Weak Guidance, All Main Models Combined)

**What is plotted:** Combined histogram + KDE of delta RSCC, with all main models overlaid on one axis. Weak guidance only.

**Data included:**
- Delta RSCC values at `guidance_strength == 'weak'` only.
- Main models only (excludes `MODELS_ALT`).
- Clipped to `delta_rscc` in (-0.3, +0.3).

**What is excluded:**
- All occupancy configurations other than 50/50.
- Baseline, medium, and strong guidance rows.
- Alternative models (`protenix-mini`, `protenix-tiny`, `rf3-old`, `rf3-chiral-off`).
- Delta RSCC values outside +/-0.3.

---

### Figure 08 — Delta RSCC (Weak Guidance, Protenix Group)

**What is plotted:** Same as Figure 07 but restricted to the Protenix group (`protenix`, `protenix-mini`, `protenix-tiny`).

**Data included:**
- Delta RSCC values at `guidance_strength == 'weak'`.
- Protenix group models only.
- Clipped to (-0.3, +0.3).

**What is excluded:**
- Same as Figure 07, except model filter is Protenix group instead of main models.

---

### Figure 09 — LDDT-to-Reference Baseline Distribution

**What is plotted:** Combined histogram + KDE of `avg_silhouette_to_ref` at baseline, all models overlaid.

**Data included:**
- LDDT data at `guidance_strength == 'baseline'`.
- Filtered to `selection_name != 'all altloc'` (only non-trivial residue selections).
- Main models only (excludes `MODELS_ALT`). Variant-specific versions also exist for Protenix and RF3 groups.

**What is excluded:**
- **Rows with `selection_name == 'all altloc'`** — the full-chain altloc selection is excluded.
- All non-baseline guidance rows.
- Alternative models (in the main version).

---

### Figure 10 — LDDT-to-Reference Across All Guidance Strengths

**What is plotted:** Ridge plot (one row per model, hue = guidance strength) of `avg_silhouette_to_ref` across all guidance levels.

**Data included:**
- All LDDT rows (baseline + all guided).
- Filtered to `selection_name != 'all altloc'`.
- Main models only (with variant versions for Protenix/RF3 groups).

**What is excluded:**
- Rows with `selection_name == 'all altloc'`.
- Alternative models (in the main version).

---

### Figure 11 — Delta LDDT-to-Reference Distribution

**What is plotted:** Ridge plot of the change in `avg_silhouette_to_ref` from baseline, per model.

**Delta computation:** The `lddt_merged` dataframe is built by:
1. Extracting baseline rows (one per protein/model/selection_name).
2. Merging onto non-baseline rows by `(protein, model_name, selection_name)`.
3. Computing `delta_silhouette_to_ref = avg_silhouette_to_ref - avg_silhouette_to_ref_baseline`.

**Data included:**
- All non-baseline (protein, model, selection, gw) combinations with a matching baseline.
- Main models only (with variant versions). The `selection_name != 'all altloc'` filter is NOT explicitly applied at the plot level, but it is inherited from the merge since baseline/guided rows for "all altloc" both exist and produce valid deltas — so "all altloc" selections ARE included in Figure 11, unlike Figures 09 and 10.

**What is excluded:**
- Alternative models (in the main version).
- Any (protein, model, selection) triple missing a baseline.

---

### RSCC Scatter Plots — Baseline vs Optimal Guidance (per model)

**What is plotted:** Scatter plot with one point per **(protein, selection)** pair for the selected model. X-axis = baseline RSCC, Y-axis = RSCC at the model's optimal guidance weight. Diagonal `y = x` spans the full axis range. Both axes fixed to `[0, 1]`. Marginal KDE density panel below.

**Highlighting** (capped so the figure stays readable):
- Top-3 most-improved selections: green ring.
- Top-3 most-decreased selections (most-negative `delta_rscc`): red ring.
- Loop-of-interest selections: blue ring. Currently `6B8X 181-183` and `5MHX 162-164` (the altloc-resolved residues inside each protein's most flexible large loop).
- Labels use the compact form `"<PDB> <start>-<end>"`.

**Optimal guidance weights:**

| Model | Optimal gw |
|---|---|
| boltz2X-RAY | 0.1 |
| boltz2MD | 0.1 |
| protenix | 0.1 |
| protenix-mini | 0.1 |
| protenix-tiny | 0.1 |
| rf3 | 0.01 |
| rf3-old | 0.01 |
| rf3-chiral-off | 0.01 |

**Data construction:**
1. Load RSCC filtered to `selection_name != "all altloc"` (drop chain-level rows).
2. Baseline: rows with `guidance_weight == 0.0`, keeping `(protein, selection, model_name, rscc)` renamed to `rscc_baseline`.
3. Optimal: for each model, rows at that model's optimal gw, same columns renamed to `rscc_optimal`.
4. Merge on `(protein, selection, model_name)`.
5. `delta_rscc = rscc_optimal - rscc_baseline`.
6. Add `label = "<PDB> <start>-<end>"` for annotations.

**Data included:**
- One point per (protein, selection) pair for the selected model, requiring both baseline and optimal-gw measurements.

**What is excluded:**
- `selection_name == "all altloc"` (chain-level rows).
- All occupancy configurations other than 50/50.
- All guidance weights other than 0.0 (baseline) and the model's designated optimal weight.
- (protein, selection) pairs missing from either slice.

---

### Aggregated RSCC Scatter — All Models Combined

**What is plotted:** Scatter with one point per **(protein, selection)** pair. Baseline and optimal RSCC are collapsed across main models using the **median** (no error bars; per-model spread is shown by the per-model RSCC scatter figures). Diagonal y = x line. Marginal KDE. Both axes fixed to `[0, 1]` so the figure is directly comparable across runs and model groups.

**Data construction:**
1. Load RSCC filtered to `selection_name != "all altloc"` (drop the chain-level rows).
2. Baseline: rows with `guidance_weight == 0.0` and `model_name not in MODELS_ALT`, keeping `(protein, selection, model_name, rscc)`.
3. Optimal: rows at each main model's optimal gw (from `OPTIMAL_GW`), same columns.
4. Merge on `(protein, selection, model_name)`, then group by `(protein, selection)` and take median / count across models.

**Highlighting** (capped so the figure stays readable):
- Top-3 most-improved (protein, selection) pairs: green ring.
- Top-3 most-decreased pairs (most-negative delta): red ring.
- Loop-of-interest pairs: blue ring. Currently: `6B8X 181-183` (altloc-resolved residues inside the most flexible large loop `chain A and resi 178-193`, backbone lDDT 0.676) and `5MHX 162-164` (altloc-resolved residues inside the only large loop `chain A and resi 158-167`, backbone lDDT 0.582).
- Labels use the compact form `"<PDB> <start>-<end>"` (e.g., `6B8X 181-183`).

**Data included:**
- One point per (protein, selection) pair that has `n >= 3` main models contributing both a baseline and an optimal measurement.

**What is excluded:**
- Alternative models (`MODELS_ALT`).
- `selection_name == "all altloc"` (chain-level rows).
- (protein, selection) pairs with fewer than 3 main models.
- Everything excluded at the upstream RSCC loader level (non-50/50 occupancies, etc.).

---

### RSCC vs LDDT Cross-Metric Scatter Plots (per model)

**What is plotted:** Scatter comparing RSCC (x) vs LDDT silhouette-to-ref (y), both at the optimal guidance weight. One point per protein. Regression line + Pearson r displayed.

**Data construction:**
1. RSCC at optimal gw per model (same as scatter plot above).
2. LDDT at optimal gw: the LDDT data is filtered to `selection_name != 'all altloc'`, then for each model at its optimal gw, `avg_silhouette_to_ref` is **averaged across all remaining selections** per protein.
3. Merged on `(protein, model_name)`.

**Data included:**
- Proteins that have both an RSCC and an LDDT measurement at the optimal gw for the selected model.

**What is excluded:**
- LDDT rows with `selection_name == 'all altloc'`.
- All non-optimal guidance weights.
- All occupancy configurations other than 50/50 (RSCC side).
- Proteins missing from either metric.

---

### Aggregated RSCC vs LDDT Scatter — All Models Combined

**What is plotted:** Same cross-metric scatter, collapsed across models using the **median** (no error bars; per-model spread is shown by the per-model cross-metric scatter figures). Pearson r.

**Data included:**
- The per-model `cross_metric_df`, grouped by protein. Median across models.
- Only proteins with `n >= 3` models.

**What is excluded:**
- Proteins in fewer than 3 models.
- Same upstream exclusions as the per-model cross-metric scatter.

---

### Aggregated LDDT Scatter — All Models Combined

**What is plotted:** Scatter with one point per **(protein, selection)** pair. Baseline and optimally-guided LDDT silhouette-to-ref collapsed across main models using the **median** (no error bars; per-model spread is shown by the per-model LDDT scatter figures). Diagonal y = x. Marginal KDE. Both axes fixed to `[0, 1]` for cross-figure comparability.

**Data construction:**
1. LDDT is filtered to `selection_name != 'all altloc'`.
2. Baseline: rows at `guidance_weight == 0.0` and `model_name not in MODELS_ALT`, keeping `(protein, selection, model_name, avg_silhouette_to_ref)`.
3. Optimal: rows at each main model's optimal gw (from `OPTIMAL_GW`), same columns.
4. Merge on `(protein, selection, model_name)`, then group by `(protein, selection)` and take median / count across models.
5. (protein, selection) pairs with `n >= 3` models kept.

**Highlighting:** Same as the Aggregated RSCC Scatter — top-3 improved (green), top-3 most-decreased (red), `6B8X 181-183` and `5MHX 162-164` (blue). Labels use `"<PDB> <start>-<end>"`.

**What is excluded:**
- Alternative models (`MODELS_ALT`).
- LDDT rows with `selection_name == 'all altloc'`.
- Non-baseline / non-optimal guidance weights.
- (protein, selection) pairs with fewer than 3 main models.

**Note:** The previous version of this figure averaged `avg_silhouette_to_ref` across selections per (protein, model) before aggregating across models. That collapse step has been removed so each residue-level selection becomes its own point.

---

## Notebook 3: Guidance-Weight Spiral Trajectories (`03_tradeoffs.py`)

### Source data

All four metric files are loaded (same as Notebooks 1 and 2):
- `initial40results/clashscore_metrics.csv`
- `initial40results/bond_angle_violation_fractions.csv`
- `initial40results/bond_length_violation_fractions.csv`
- `initial40results/rscc_results.csv`

### Tradeoff dataframe construction

This is the most complex data pipeline. The steps are:

1. **Load and tag** clashscores, bond stats, and RSCC (same loaders/filters as Notebooks 1 and 2).

2. **Aggregate per (protein, model, gw, guidance_strength):**
   - `clashscore`: mean across rows within each group.
   - `outlier_fraction_bond_angles`: mean across rows.
   - `outlier_fraction_bond_lengths`: mean across rows.
   - `rscc`: mean across rows.

3. **Compute deltas** for each metric via `make_delta_frame`:
   - `delta_clashscore = clashscore_guided - clashscore_baseline`
   - `delta_outlier_fraction_bond_angles = guided - baseline`
   - `delta_outlier_fraction_bond_lengths = guided - baseline`
   - `delta_rscc = rscc_guided - rscc_baseline`

4. **Merge** all four delta frames on `(protein, model_name, guidance_weight, guidance_strength)`. The RSCC merge uses `how="left"`, so rows with geometry data but missing RSCC are kept (with `delta_rscc = NaN`).

5. **Convert bond fractions to percentage points:**
   - `delta_outlier_fraction_bond_lengths_pct = 100 * delta_outlier_fraction_bond_lengths`
   - `delta_outlier_fraction_bond_angles_pct = 100 * delta_outlier_fraction_bond_angles`

6. **Compute geometry penalty** (unweighted mean of three terms):
   ```
   geometry_penalty = mean(
       delta_clashscore,
       delta_outlier_fraction_bond_lengths_pct,
       delta_outlier_fraction_bond_angles_pct
   )
   ```

The final `tradeoff_df` has one row per (protein, model, guidance_weight) at every non-baseline guidance level, with columns `delta_rscc` and `geometry_penalty`.

---

### Figures 12–19 — Per-Model Spiral Trajectories

**Figure numbering:** 12 = boltz2X-RAY, 13 = boltz2MD, 14 = protenix, 15 = rf3, 16 = protenix-mini, 17 = protenix-tiny, 18 = rf3-old, 19 = rf3-chiral-off.

**What is plotted:** For a single model, one trajectory line per protein traces through (delta RSCC, geometry penalty) space as guidance weight increases. Each trajectory starts at the origin (0, 0) representing baseline. Line segments and dots are colored by guidance weight using the "flare" colormap (cool-to-warm). Arrowheads indicate direction of increasing guidance. The top-5 most extreme trajectories (by total arc length in tradeoff space) are labeled with the protein name. Y-axis uses symmetric log scale (`symlog`, `linthresh=1.0`).

**Data included:**
- All rows in `tradeoff_df` for the selected model.
- All non-baseline guidance weights present for that model.
- All proteins.

**What is excluded:**
- All occupancy configurations other than 50/50 (global filter).
- Baseline rows (they define the origin, not plotted as explicit data points).
- `inf`/`-inf` values in delta_rscc or geometry_penalty are replaced with NaN and excluded from axis-limit computation.
- **No outlier clipping** is applied — unlike Figures 02/03/06, all values are shown.

**Shared axis limits:** X and Y limits are computed across ALL models simultaneously to enable direct visual comparison between figures 12–19.

---

### Figure 20 — Aggregated Spiral (All Models Overlaid)

**What is plotted:** Multi-model overlay. For each (model, guidance_weight), all proteins are collapsed to a single point at the **median** delta RSCC and **median** geometry penalty. IQR error bars (Q1–Q3) in both X and Y. Trajectory line from origin (0,0) through each guidance weight, with arrowhead. Y-axis uses symlog scale. Fixed `ylim=[-1, 1000]`.

**Models included in Figure 20:** `boltz2X-RAY`, `boltz2MD`, `protenix`, `rf3`, `rf3-chiral-off` (hardcoded model_order in the notebook cell). Note: `protenix-mini`, `protenix-tiny`, and `rf3-old` are **not included** in this specific figure.

**Data included:**
- Aggregated statistics from `tradeoff_df` grouped by `(model_name, guidance_weight)`.
- Only groups with **n >= 3 proteins** (the `min_n` default).

**What is excluded:**
- `protenix-mini`, `protenix-tiny`, `rf3-old` (not in the hardcoded model order for this figure).
- Groups with fewer than 3 proteins.
- All upstream exclusions (occupancy filter, baseline-only subtraction).

---

### Figures 21–28 — Per-Model Aggregated Spirals

**Figure numbering:** 21 = boltz2X-RAY, 22 = boltz2MD, 23 = protenix, 24 = rf3, 25 = protenix-mini, 26 = protenix-tiny, 27 = rf3-old, 28 = rf3-chiral-off.

**What is plotted:** Single-model aggregated spiral. Same aggregation as Figure 20 but for one model at a time. Dots colored by guidance weight (flare colormap). IQR error bars. Colorbar legend for guidance weight.

**Data included:**
- All (model, guidance_weight) groups for the selected model with n >= 3 proteins.

**What is excluded:**
- Groups with fewer than 3 proteins.
- All upstream exclusions.

---

### Figures 29–36 — Training-Distribution Bias

**Common inputs:**
- `lddt_assignment_counts.csv` — pre-cutoff (PDB depositions before 2021-09-30) per-selection counts of cluster members assigned to altloc A vs altloc B by lDDT comparison, for each protein's sequence cluster. Loaded via `load_training_cluster_counts` (lowercases `protein`, renames `selection_part` → `selection`, drops the multi-part `selection` column to avoid name collision, computes `frac_a_cluster = n_A / (n_A + n_B)` and `x_one_minus_fa = 1 − frac_a_cluster`).
- `initial40results/Sampleworks_initial_structures.csv` — per-protein deposition date and `Pre`/`Post` cutoff flag. Loaded via `load_deposition_cutoffs` (lowercases PDB ID, dedupes the duplicate `9B22` row, fills the four missing proteins `1VME, 4OLE, 5SOP, 6B8X` as `Pre` since all are 2004–2017 depositions).
- `initial40results/lddt_results.csv` — predicted occupancies; loaded via `load_lddt` (so the global `OCC_FILTER_ALTLOC` of `{'A': 0.5, 'B': 0.5}` and `MODEL_NAME_ALIAS` are already applied), then narrowed to main models via `filter_by_model_group(df, "main")` (drops `protenix-mini`, `protenix-tiny`, `rf3-old`, `rf3-chiral-off`).

**Common transforms:** `build_training_bias_df` parses the **second** float of each `occupancies` cell (numpy-style `"[0.375 0.625]"` → `0.625`) into a `b_occupancy` column, inner-joins lddt to the cluster frame on `(protein, selection)`, and left-joins the cutoff flag on `protein`. Warnings are printed for: NaN occupancies (~46 rows under the OCC filter), lddt rows with no matching cluster row (~22 keys, mostly `7z0e` chain P), and ambiguous cluster rows where `n_cluster_evaluated > 0` but `A == B == 0` (21 rows, dropped — X is undefined for plot 1, eval > 0 disqualifies plot 2).

**Models plotted:** boltz2X-RAY, boltz2MD, protenix, rf3.

#### Figure 29 — Training-Cluster vs Predicted State-B Fraction (Small Multiples)

**Output:** `figures/29_training_bias_small_multiples.svg`.

**What is plotted:** 4-column × 2-row grid of small scatter panels. Columns: Boltz2-XRD, Boltz2-MD, Protenix, RF3 (one per main model). Rows: top = Baseline (`guidance_weight == 0.0`), bottom = Optimal guidance (each model at its `OPTIMAL_GW`: `0.1` for Boltz2-XRD/Boltz2-MD/Protenix, `0.01` for RF3). Axes are shared across all 8 panels: `xlim = ylim = [-0.02, 1.02]`. Reference line `y = 0.5` (light gray dotted) drawn behind the points. Model identity is encoded by column position only — points are monochrome dark gray.

**Axes:**
- X = `frac_b_cluster = n_assigned_to_B / (n_assigned_to_A + n_assigned_to_B)` — fraction of pre-cutoff cluster assigned to altloc B by lDDT comparison. Computed directly from the loader's count columns (not via `1 − frac_a_cluster`).
- Y = `b_occupancy` — model-predicted fraction in state B, parsed from the second float of each `occupancies` cell.

**Point encoding:**
- Color: monochrome `#333333` (no model color).
- Alpha: `alpha = log(s + 1) / log(s_max_global + 1)` where `s = n_assigned_to_A + n_assigned_to_B` and `s_max_global` is the maximum `s` over the union of the baseline and optimal frames. Alpha is comparable across all 8 panels by construction.
- Jitter: x is perturbed by uniform `±0.01` (seed = 0) so boundary pileups at `frac_b_cluster ∈ {0, 1}` are visible.

**In-panel decoration:** No in-panel text. Column titles (top row only) use `MODEL_DISPLAY`. Row labels "Baseline" / "Optimal guidance" are placed to the left of the leftmost panel of each row. Figure-level X and Y labels: `fraction in state B (training cluster)` and `predicted fraction in state B`.

**Data included:**
- One point per `(protein, selection_part, model_name, condition)` combination with `s > 0`.

**What is excluded:**
- `s = 0` cases (no cluster member within RMSD threshold of either altloc) are dropped from this figure — they are part of the pre-cutoff small-cluster subset in Figures 30–33.
- All upstream exclusions (occupancy filter, model alias, non-main models).

#### Figure 34 — Calibration Overlay (Redesign A, 1x4)

**Output:** `figures/34_training_bias_calibration_overlay.svg`.

**What is plotted:** 1-row × 4-column grid of panels, one per main model (Boltz2-XRD, Boltz2-MD, Protenix, RF3). Each panel pairs a scatter axis (left, width 4) with a marginal y-histogram strip (right, width 1) sharing y, via nested `add_gridspec` / `subgridspec`. Within each scatter, **Baseline** (gray, `#888888`) and **Optimal guidance** (dark, `#333333`) points are overlaid. Both reference lines are drawn behind the points: `y = x` (dashed `#AAAAAA`, lw=1, "training prior") and `y = 0.5` (dotted `#CCCCCC`, lw=0.8, "guidance target"). Per-condition binned conditional-mean curves are overlaid: baseline as a thin gray line through bin-center dots (`color="#888888", lw=1.5, alpha=0.8, marker="o", markersize=4`); optimal guidance as a thicker dark line (`color="#222222", lw=2.5, alpha=0.95, marker="o", markersize=5`). Conditional means use 10 **equal-count (quantile)** bins of `frac_b_cluster`; duplicate quantile edges from ties at 0 / 1 are collapsed (so the effective bin count can be less than 10), and each bin's dot sits at the within-bin **mean of x** rather than the bin midpoint. Bins with fewer than 5 points are skipped (NaN sentinel breaks the line). Slope of `b_occupancy` on `frac_b_cluster` and Pearson r are annotated for the **Optimal-guidance** condition only, in the upper-left corner of each scatter axis with the rounded white box style used in `make_scatter_comparison`. Each marginal axis overlays both conditions' `b_occupancy` distributions as horizontal histograms (`bins=np.linspace(0, 1, 21)`, baseline `alpha=0.45`, guidance `alpha=0.65`). A single horizontal legend above the panel row names the four artists ("Baseline (points & mean)", "Optimal guidance (points & mean)", "y = x (training prior)", "y = 0.5 (guidance target)").

**Axes:**
- X = `frac_b_cluster` (same definition as Figure 29).
- Y = `b_occupancy` (same definition as Figure 29).
- `xlim = ylim = [-0.02, 1.02]` for every panel; ticks at `[0, 0.5, 1]`.

**Point encoding:**
- Color: baseline `#888888`, optimal guidance `_TRAINING_BIAS_POINT_COLOR` = `#333333` (no per-model color).
- Alpha: `alpha = 0.5 * log(s + 1) / log(s_max_global + 1)` then clipped to [0, 1], with `s_max_global` taken across the **union** of the baseline and optimal prepared frames so alpha is comparable across all panels. The half-strength factor lets the conditional-mean curves dominate visually.
- Marker size: `s=12` (smaller than Figure 29's `s=18`).
- Jitter: x is perturbed by uniform `±0.01` (seed = 0), drawn from a single shared RNG so the pattern is reproducible.

**In-panel decoration:** Direct text labels "training prior" (at data coords `(0.78, 0.86)`) and "guidance target" (at `(0.78, 0.55)`) are placed in the **leftmost panel only** (`fontsize=8, color="#666666"`). Column titles use `MODEL_DISPLAY`. The rightmost marginal axis carries an `x` label `"count"` (`fontsize=8`); other marginals hide the x spine and ticks.

**Data included:**
- Same as Figure 29: one point per `(protein, selection_part, model_name, condition)` combination with `s > 0`, for both Baseline and Optimal-guidance conditions.

**What is excluded:**
- Same as Figure 29 (`s = 0` rows; all upstream filters: occupancy, model alias, non-main models).

#### Figure 35 — Calibration Small Multiples v2 (Redesign B, 2x4)

**Output:** `figures/35_training_bias_small_multiples_v2.svg`.

**What is plotted:** 2-row × 4-column grid mirroring Figure 29's geometry (top = Baseline, bottom = Optimal guidance, columns = main models). Each scatter panel is paired with a single-condition marginal y-histogram strip on its right edge (width ratio 4:1 via nested `subgridspec`). Unlike Redesign A, the marginals are **per-panel** — each marginal shows only its own panel's condition. Every panel draws both reference lines (`y = x` dashed `#AAAAAA` and `y = 0.5` dotted `#CCCCCC`) — including the Baseline row, where the contrast against a near-zero slope is part of the rhetorical comparison. Each panel overlays its own binned conditional-mean curve (`color="#222222", lw=2.0, marker="o", markersize=4`, same equal-count-quantile / `min_count=5` logic as Figure 34) on top of a `mean ± std` shadow band (`color="#222222", alpha=0.15`, `ddof=1` sample std per bin, NaN where the bin is below the count threshold). Slope and Pearson r for the panel's own condition are annotated upper-left (same bbox style as Figure 34). The marginal histogram uses `bins=np.linspace(0, 1, 21), color=_TRAINING_BIAS_POINT_COLOR, alpha=0.6`.

**Axes:**
- X = `frac_b_cluster` (same definition as Figure 29).
- Y = `b_occupancy` (same definition as Figure 29).
- All scatter panels share x and y. `xlim = ylim = [-0.02, 1.02]`; ticks at `[0, 0.5, 1]`.

**Point encoding:**
- Color: monochrome `_TRAINING_BIAS_POINT_COLOR` (`#333333`); model identity is encoded by column position (Figure 29 convention).
- Alpha: same formula as Figure 34 (half-strength log-scaled with `s_max` over the union of both frames).
- Marker size: `s=12`. Jitter: uniform `±0.01` (seed = 0).

**In-panel decoration:** No inline reference-line labels and no inset legend in this figure (they will be added separately). Column titles on the top row only, using `MODEL_DISPLAY`. Row labels "Baseline" / "Optimal guidance" are placed to the left of the leftmost panel of each row via the same `axes[row, 0].annotate(label, xy=(-0.32, 0.5), xycoords="axes fraction", rotation=90, fontsize=11)` pattern as `make_training_bias_scatter`. Marginal axes hide their x spine and tick labels.

**Data included:**
- Same as Figure 29 / Figure 34: one point per `(protein, selection_part, model_name, condition)` combination with `s > 0`, plotted in its panel's row × column cell.

**What is excluded:**
- Same as Figure 29 (`s = 0` rows; all upstream filters: occupancy, model alias, non-main models).

#### Figure 36 — Calibration Explainer (Empty Panel)

**Output:** `figures/36_training_bias_calibration_explainer.svg`.

**What is plotted:** A single 4×4 inch panel with the same axis style as the v2 small multiples (`xlim = ylim = [-0.02, 1.02]`, ticks at `[0, 0.5, 1]`, axis labels `fraction in state B (training cluster)` / `predicted fraction in state B`). Three lines, no data:
- `y = x` (training prior) — dashed `#AAAAAA`, lw=1.
- `y = 0.5` (guidance target) — dashed `#AAAAAA`, lw=1 (note: dashed here, not dotted as in Figures 34/35).
- Solid bisector at `y = 0.5*x + 0.25` — `color="#222222", lw=2.0`. Slope 0.5 through the intersection of the two reference lines at `(0.5, 0.5)`; runs from `(-0.02, 0.24)` to `(1.02, 0.76)`.

No annotations, no title. Intended as a minimal explanatory companion to Figure 35 for downstream labeling.

**Data included / excluded:** N/A (no data).

#### Figures 30, 31 — Calibration Overlay on Out-of-Training-Distribution Subsets

**Output filenames:**
- `figures/30_training_bias_calibration_overlay_low_cluster.svg`
- `figures/31_training_bias_calibration_overlay_post_cutoff.svg`

**What is plotted:** Both figures use the exact Figure 34 layout (`make_training_bias_calibration_overlay`): a 1×4 row of panels, one per main model (`boltz2X-RAY`, `boltz2MD`, `protenix`, `rf3`), with Baseline and Optimal guidance overlaid in the same axes, binned conditional-mean curves for each condition, `y = x` (training prior) and `y = 0.5` (guidance target) reference lines, slope and Pearson r annotated for both conditions, and a per-panel marginal y-histogram on the right edge. The only difference from Figure 34 is the subset filter applied to both conditions before plotting (passed via the function's `max_cluster_size` and `cutoff_status` parameters).

- **Figure 30 (low-cluster subset):** `max_cluster_size=5` — keeps only rows with `n_cluster_evaluated ≤ 5` (sparse training prior on both sides of the cutoff). Pre/post-cutoff is not filtered here.
- **Figure 31 (post-cutoff subset):** `cutoff_status="Post"` — keeps only rows deposited after 2021-09-30 (`cutoff_status == "Post"`), regardless of cluster size.

**Data included:**
- Same upstream pipeline as Figure 34 (`build_training_bias_df` → `_training_bias_prepare_scatter_frame` requiring `s_cluster > 0`), plus the subset filter for each figure.
- For each panel, baseline (`guidance_weight == 0.0`) and optimal-guidance (`OPTIMAL_GW[model]`) conditions are both drawn.

**What is excluded:**
- Same upstream exclusions as Figure 34 (rows missing `b_occupancy` or `frac_b_cluster`; rows with `s_cluster ≤ 0`; non-main models).
- Figure 30: rows with `n_cluster_evaluated > 5`. Note: with only ~9 surviving points per panel, the binned mean uses `n_bins=3, min_bin_count=2` (overrides forwarded to `_training_bias_binned_means`) so a coarse trend curve is still drawn.
- Figure 31: rows with `cutoff_status != "Post"` (the small set with missing `cutoff_status` is also dropped).

#### Figure 32 — Empty-Training-Cluster Strip Plot

**Output:** `figures/32_training_bias_no_cluster_strip.svg`.

**What is plotted:** 1×4 row of panels (one per main model). For each panel, a jittered strip plot overlays Baseline (`#888888`) and Optimal-guidance (`MODEL_PALETTE[model]`) predicted B-occupancy values; the X axis carries no quantitative meaning (uniform jitter on `[-0.35, 0.35]`, seed=0) and its ticks are hidden. Each scatter axis has a marginal y-histogram on its right edge using the same 21-edge bins as Figures 34/35, with both conditions overlaid (`alpha=0.45` for Baseline, `0.65` for Optimal). A `y = 0.5` reference line is drawn behind the points (dashed, `#CCCCCC`). Each panel's Optimal-condition sample size is annotated in the top-left corner. Panel title = `MODEL_DISPLAY[model_name]`.

This figure exists because the rows it shows have `frac_b_cluster = N_B / (N_A + N_B)` undefined (both numerator and denominator zero), so they cannot appear on the X-axis of Figures 34, 30, 31, or 35. They are nonetheless the most extreme out-of-distribution case (no pre-cutoff PDB cluster member assigned to either altloc by lDDT comparison) and worth visualizing separately.

**Data included:**
- All `build_training_bias_df` rows where `n_assigned_to_A == 0` AND `n_assigned_to_B == 0`, for the four main models, at Baseline (gw=0.0) and Optimal guidance (`OPTIMAL_GW[model]`).
- At the time of writing this is ~117 selections per model panel for protenix and ~134 for the other three main models, spanning 14-17 distinct proteins.

**What is excluded:**
- Any row with at least one non-zero altloc count (those appear in Figures 29-31, 34, 35 instead).
- Rows missing `b_occupancy` (NaN-occupancy rows are dropped upstream by `build_training_bias_df`).
- Non-main models.

---

## Summary: What Is Universally Excluded

The following data **never appears in any plot**:

1. **Non-50/50 occupancy configurations** — only the equal-occupancy (A=0.5, B=0.5) conformer split is analyzed.
2. **Proteins or models missing from the source CSVs** — the analysis is limited to whatever was computed in `initial40results/`.
3. **Baseline rows in delta plots** — deltas are defined as guided minus baseline, so the baseline itself (gw=0.0) does not appear as a data point in any delta/spiral figure (it is the implicit origin).
4. **(protein, model) pairs missing a baseline measurement** — the `make_delta_frame` merge is an inner join, so if there is no gw=0.0 row for a (protein, model), that pair is absent from all delta-based figures.

## Summary: Per-Figure Exclusion Table

| Figure | Excludes strong guidance | Excludes "all altloc" selection | Outlier/range clipping | Model subset |
|---|---|---|---|---|
| 01 Clashscore Distribution | Yes | Yes (no-op) | None | By dropdown |
| 02 Delta Clashscore | No | No | `delta >= 100` excluded | By dropdown |
| 03 Delta Bond Angles | No | No | `delta >= 0.05` excluded | By dropdown |
| 04 RSCC Baseline | N/A (baseline only) | No | None | By dropdown |
| 05 RSCC Baseline vs Weak | Medium+strong excluded | No | None | By dropdown |
| 06 Delta RSCC | No | No | `abs(delta) >= 0.3` excluded | By dropdown |
| 07 Delta RSCC Weak All | Medium+strong excluded | No | `abs(delta) >= 0.3` excluded | Main only |
| 08 Delta RSCC Weak Protenix | Medium+strong excluded | No | `abs(delta) >= 0.3` excluded | Protenix group |
| 09 LDDT Baseline | N/A (baseline only) | Yes | None | Main (+ variants) |
| 10 LDDT All Guidance | No | Yes | None | Main (+ variants) |
| 11 Delta LDDT | No | No (included) | None | Main (+ variants) |
| 12–19 Per-Model Spirals | No | N/A (geometry+RSCC) | None | One model each |
| 20 Aggregated Spiral All | No | N/A | n < 3 groups dropped | 5 models (see above) |
| 21–28 Aggregated Spirals | No | N/A | n < 3 groups dropped | One model each |
| RSCC Scatter (per model) | N/A (specific gw) | No | None | One model each |
| RSCC Scatter Aggregated | N/A (specific gw) | No | n < 3 proteins dropped | All 8 models |
| RSCC vs LDDT (per model) | N/A (specific gw) | Yes (LDDT side) | None | One model each |
| RSCC vs LDDT Aggregated | N/A (specific gw) | Yes (LDDT side) | n < 3 proteins dropped | All 8 models |
| LDDT Scatter Aggregated | N/A (specific gw) | Yes | n < 3 proteins dropped | All 8 models |
| 34 Calibration Overlay | No | N/A (cluster-frame join) | None | Main only |
| 35 Calibration Small Multiples v2 | No | N/A (cluster-frame join) | None | Main only |
| 36 Calibration Explainer | N/A (no data) | N/A (no data) | None | N/A (no data) |
