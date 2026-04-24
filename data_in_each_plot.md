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

**What is plotted:** Scatter plot with one point per protein. X-axis = baseline RSCC, Y-axis = RSCC at optimal guidance weight. Diagonal line at y = x. Top-3 most-improved proteins highlighted green, any decreased proteins highlighted red. Marginal KDE density panel below.

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
1. Baseline: all rows with `guidance_weight == 0.0`, extracting `(protein, model_name, rscc)` renamed to `rscc_baseline`.
2. Optimal: for each model, rows at that model's optimal gw, extracting `rscc` renamed to `rscc_optimal`.
3. Merged on `(protein, model_name)`.
4. `delta_rscc = rscc_optimal - rscc_baseline`.

**Data included:**
- One point per protein for the selected model.
- Only proteins that have BOTH a baseline and an optimal-gw measurement for that model.

**What is excluded:**
- All occupancy configurations other than 50/50.
- All guidance weights other than 0.0 (baseline) and the model's designated optimal weight.
- Proteins missing from either the baseline or optimal slice.

---

### Aggregated RSCC Scatter — All Models Combined

**What is plotted:** Scatter with one point per protein. Baseline and optimal RSCC are collapsed across all 8 models using **median**. Error bars show IQR (Q1–Q3) across models. Diagonal y = x line. Marginal KDE.

**Data included:**
- The per-model `rscc_scatter_df` (baseline vs optimal for each model), grouped by protein.
- Statistics: median, Q1, Q3 of both `rscc_baseline` and `rscc_optimal` across models.
- **Only proteins present in >= 3 models** (`n >= 3` filter).

**What is excluded:**
- Proteins appearing in fewer than 3 models.
- Everything excluded at the per-model scatter level (non-50/50 occupancies, non-optimal gw, etc.).

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

**What is plotted:** Same cross-metric scatter, collapsed across models. Median +/- IQR error bars. Pearson r.

**Data included:**
- The per-model `cross_metric_df`, grouped by protein. Median and IQR across models.
- Only proteins with `n >= 3` models.

**What is excluded:**
- Proteins in fewer than 3 models.
- Same upstream exclusions as the per-model cross-metric scatter.

---

### Aggregated LDDT Scatter — All Models Combined

**What is plotted:** Scatter with one point per protein. Baseline and optimally-guided LDDT silhouette-to-ref collapsed across models (median +/- IQR). Diagonal y = x. Marginal KDE.

**Data construction:**
1. LDDT is filtered to `selection_name != 'all altloc'`.
2. Baseline: rows at `guidance_weight == 0.0`, `avg_silhouette_to_ref` averaged across selections per (protein, model).
3. Optimal: rows at each model's optimal gw, same averaging.
4. Merged per (protein, model), then aggregated across models per protein (median, Q1, Q3).
5. Proteins with `n >= 3` models kept.

**What is excluded:**
- LDDT rows with `selection_name == 'all altloc'`.
- Non-baseline/non-optimal guidance weights.
- Proteins in fewer than 3 models.

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
