# Data in Each Plot

This document describes, for every figure produced by the analysis notebooks, exactly what data goes in, how it is filtered and transformed, and what is left out.

---

## Global Filters (applied at data-loading time)

These filters are baked into the data-loading functions in `loaders.py` and apply to **every** downstream plot that uses that loader. Any data excluded here never reaches any figure.

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
- **MODELS_ALT** (alternative/variant models): `protenix-mini`, `protenix-tiny`, `rf3-chiral-off`
- **Protenix group**: `protenix`, `protenix-mini`, `protenix-tiny`
- **RF3 group**: `rf3`, `rf3-chiral-off`

---

## Notebook 1: Fit Quality (`02_fit_quality.py`)

### Source data

| File | Key columns |
|---|---|
| `initial40results/rscc_results.csv` | `protein`, `model`, `method`, `altloc_occupancies`, `guidance_weight`, `rscc`, `selection`, `ensemble_size` |
| `initial40results/lddt_results.csv` | `protein`, `model`, `method`, `altloc_occupancies`, `guidance_weight`, `selection`, `occupancies`, `avg_silhouette`, `avg_silhouette_to_ref` |
| `initial40results/rscc_results_domain_shift.csv` | Same schema as `rscc_results.csv`; clean extract of the domain-shift whole-domain rows from the legacy `rscc_results.csv.bak`, with the legacy `rf3` (= `rf3-old`) rows dropped. |
| `small-altlocs-classified-rmsd.csv` | `protein`, `selection`, `chain`, `start_res`, `end_res`, `classification` |

**RSCC aggregation step:** After loading, RSCC rows are grouped by `(protein, model_name, guidance_weight)` and averaged across all per-altloc-region selections (`.mean()`), yielding one RSCC value per (protein, model, gw) triple.

**LDDT baseline merge:** Baseline LDDT rows (`guidance_strength == "baseline"`) are extracted with their `avg_silhouette_to_ref` value, then merged back onto non-baseline rows by `(protein, model_name, selection_name)` to compute `delta_silhouette_to_ref = guided - baseline`.

---

### Figures 04–09 — Combined Baseline Distributions (`04_09_baseline_distributions[_variant].svg`)

**What is plotted:** A two-panel figure: RSCC baseline distribution on top, LDDT-to-reference baseline distribution on the bottom, both with all main models overlaid on a single axis (one color per model). A single shared model legend in the upper-right of the top panel. Supports a `cumulative` switch (KDE or ECDF cumulative views via `cumulative_style`).

**Data included:**
- RSCC panel: aggregated RSCC at `guidance_strength == 'baseline'`, restricted to the current dropdown model group.
- LDDT panel: LDDT at `guidance_strength == 'baseline'` with `selection_name != 'all altloc'` and `model_name not in MODELS_ALT`.

**What is excluded:**
- All occupancy configurations other than 50/50 (global loaders).
- All guided rows.
- For LDDT: chain-level "all altloc" selections.

---

### RSCC Scatter Plots — Baseline vs Optimal Guidance (per model)

**What is plotted:** Scatter plot with one point per **(protein, selection)** pair for the selected model. X-axis = baseline RSCC, Y-axis = improvement (residual = guided − baseline). Reference lines `y = 0` (no change) and `y = 1 − x` (max-improvement ceiling). X axis pinned to `[0, 1]`; Y pinned to `[min(0, residual_min) - 0.02, 1.0]`. Marginal KDE density panel below.

**Highlighting** (capped so the figure stays readable):
- Top-3 most-improved selections: green ring.
- Every decreased selection (negative `delta_rscc`): red ring.
- 5 selections furthest below the `y = 1 − x` ceiling: orange ring.
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

### Aggregated RSCC Scatter — All Main Models Combined

**What is plotted:** Same residual scatter as the per-model figures, but each point is one **(protein, selection)** pair collapsed across main models using the **median** (no error bars; per-model spread is shown by the per-model figures). Marginal KDE.

**Data construction:**
1. Same per-model construction as above.
2. Restrict to `model_name not in MODELS_ALT` (main models only).
3. Group by `(protein, selection)`; take median of `rscc_baseline`, `rscc_optimal`, and count of contributing models.
4. Drop groups with fewer than 3 main models contributing.

**What is excluded:**
- Alternative models (`MODELS_ALT`).
- `selection_name == "all altloc"`.
- (protein, selection) pairs with fewer than 3 main models.

---

### RSCC vs LDDT Cross-Metric Scatter Plots (per model)

**What is plotted:** Scatter comparing RSCC (x) vs LDDT silhouette-to-ref (y), both at the optimal guidance weight. One point per protein. Regression line + Pearson r displayed.

**Data construction:**
1. RSCC at optimal gw per model (same as scatter plot above).
2. LDDT at optimal gw: filtered to `selection_name != 'all altloc'`, then for each model at its optimal gw, `avg_silhouette_to_ref` is **averaged across all remaining selections** per protein.
3. Merged on `(protein, model_name)`.

**Data included:**
- Proteins that have both an RSCC and an LDDT measurement at the optimal gw for the selected model.

**What is excluded:**
- LDDT rows with `selection_name == 'all altloc'`.
- All non-optimal guidance weights.
- All occupancy configurations other than 50/50 (RSCC side).
- Proteins missing from either metric.

---

### Aggregated LDDT Scatter — All Main Models Combined

**What is plotted:** Same residual scatter style as the aggregated RSCC scatter, but for LDDT silhouette-to-ref. One point per **(protein, selection)** pair; baseline and optimally-guided LDDT collapsed across main models using the **median**.

**Data construction:**
1. LDDT filtered to `selection_name != 'all altloc'`.
2. Baseline: `guidance_weight == 0.0` and `model_name not in MODELS_ALT`.
3. Optimal: each main model at its `OPTIMAL_GW`.
4. Merge on `(protein, selection, model_name)`, group by `(protein, selection)`, take median / count.
5. Drop groups with fewer than 3 main models.

**Highlighting:** Same as the Aggregated RSCC scatter — top-3 improved (green), top-5 most-decreased (red), top-3 furthest below ceiling (orange), `6B8X 181-183` and `5MHX 162-164` (blue).

**What is excluded:**
- Alternative models.
- LDDT rows with `selection_name == 'all altloc'`.
- Non-baseline / non-optimal guidance weights.
- (protein, selection) pairs with fewer than 3 main models.

---

### Combined By-Class Grid (`rscc_by_class_combined.svg`)

**Output:** `figures/rscc_scatterplots/rscc_by_class_combined.svg`.

**What is plotted:** A 4 × 4 grid: rows = main models (Boltz2-XRD, Boltz2-MD, Protenix, RF3), columns = conformational-change classification (`side_chain_only`, `small_loop`, `large_loop`, `domain_shift`). All 16 panels share x and y axes. X is pinned to `[0, 1]`; y is pinned globally to `[min(0, global_min_residual) − 0.02, 1.0]` so every panel exposes the same `y = 1 − x` improvement ceiling. Each panel draws (a) a light-gray context shadow consisting of **that row's per-model pooled (baseline, residual) cloud**, (b) the per-class points overlaid in the classification color, (c) `y = 0` (no-change) and `y = 1 − x` (max-improvement ceiling) reference lines. Column titles on the top row report `CLASSIFICATION_DISPLAY[cls]` plus the unique-`(protein, selection)` count for that class. Each model row carries a centered bold banner above the row (row 0's banner sits above the column titles so it doesn't collide with them). One figure-level `supxlabel("Baseline RSCC")` and one `supylabel("Improvement (Guided − Baseline RSCC)")`.

**Data construction (in `02_fit_quality.py`'s by-class cell):**
1. **Small-segment classes** (`side_chain_only`, `small_loop`, `large_loop`): load `initial40results/rscc_results.csv`, drop `selection_name == "all altloc"`, run `merge_classification_by_containment` against `small-altlocs-classified-rmsd.csv`, exclude `classification == "domain_shift"`. Baseline = `guidance_weight == 0.0`; Optimal = each model at its `OPTIMAL_GW`. Merge on `(protein, selection, model_name)`.
2. **Domain-shift class**: load `initial40results/rscc_results_domain_shift.csv` (the clean extract; the live `rscc_results.csv` only contains 1–3 residue altloc-resolved selections inside each domain, which under-samples the conformational change). Inner-join with the `domain_shift` rows of `small-altlocs-classified-rmsd.csv` on `(protein, selection)`. Same baseline / optimal slicing.
3. Concat the two frames; restrict to main models (`model_name not in MODELS_ALT`).
4. `make_by_class_combined_grid(merged_df)` does the rest.

**Data included:**
- One point per (protein, selection, model) triple in the merged frame.

**What is excluded:**
- Alternative models (`MODELS_ALT`).
- `selection_name == "all altloc"` (small-segment branch).
- All occupancy configurations other than 50/50.
- All guidance weights other than 0.0 (baseline) and each model's optimal weight.
- The legacy `rf3` (= `rf3-old`) rows in the domain-shift CSV are dropped at extraction time.

---

## Notebook 2: Guidance-Weight Spiral Trajectories (`03_tradeoffs.py`)

### Source data

All four metric files are loaded:
- `initial40results/clashscore_metrics.csv`
- `initial40results/bond_angle_violation_fractions.csv`
- `initial40results/bond_length_violation_fractions.csv`
- `initial40results/rscc_results.csv`

### Tradeoff dataframe construction

1. **Load and tag** clashscores, bond stats, and RSCC.
2. **Aggregate per (protein, model, gw, guidance_strength):** mean of `clashscore`, `outlier_fraction_bond_angles`, `outlier_fraction_bond_lengths`, `rscc`.
3. **Compute deltas** for each metric via `make_delta_frame` (`delta_X = X_guided − X_baseline`).
4. **Merge** all four delta frames on `(protein, model_name, guidance_weight, guidance_strength)`. The RSCC merge uses `how="left"`.
5. **Convert bond fractions to percentage points** (`*_pct = 100 * *`).
6. **Compute geometry penalty** = unweighted mean of `delta_clashscore`, `delta_outlier_fraction_bond_lengths_pct`, `delta_outlier_fraction_bond_angles_pct`.

The final `tradeoff_df` has one row per (protein, model, guidance_weight) at every non-baseline guidance level, with columns `delta_rscc` and `geometry_penalty`.

---

### Figures 12–17, 19 — Per-Model Spiral Trajectories

**Figure numbering:** 12 = boltz2X-RAY, 13 = boltz2MD, 14 = protenix, 15 = rf3, 16 = protenix-mini, 17 = protenix-tiny, 19 = rf3-chiral-off.

**What is plotted:** Per-model trajectory plot. One line per protein traces through (delta RSCC, geometry penalty) space as guidance weight increases. Each trajectory starts at the origin (baseline). Line segments and dots colored by guidance weight using the "flare" colormap. Arrowheads indicate the direction of increasing guidance. Top-5 most extreme trajectories (by total arc length) are labeled with the protein name. Y-axis uses symmetric log scale (`symlog`, `linthresh=1.0`).

**Data included:**
- All rows in `tradeoff_df` for the selected model.

**What is excluded:**
- All occupancy configurations other than 50/50 (global filter).
- Baseline rows (they define the origin).
- `inf`/`-inf` in delta_rscc / geometry_penalty replaced with NaN and excluded from axis-limit computation.

**Shared axis limits:** X and Y limits are computed across **all models** simultaneously to enable direct visual comparison between figures 12–17 and 19.

---

### Figure 20 — Aggregated Spiral (All Main Models Overlaid)

**What is plotted:** Multi-model overlay. For each (model, guidance_weight), proteins are collapsed to a single point at the **median** delta RSCC and **median** geometry penalty. IQR error bars in both X and Y. Trajectory line from origin through each guidance weight, with arrowhead. Y-axis on symlog. Fixed `ylim=[-1, 1000]`.

**Models included:** `boltz2X-RAY`, `boltz2MD`, `protenix`, `rf3` (hardcoded model_order in the notebook cell).

**Data included:**
- Aggregated statistics from `tradeoff_df` grouped by `(model_name, guidance_weight)`.
- Only groups with `n ≥ 3 proteins`.

**What is excluded:**
- `protenix-mini`, `protenix-tiny`, `rf3-chiral-off` (not in the hardcoded model order for this figure).
- Groups with fewer than 3 proteins.

---

### Figures 21–26, 28 — Per-Model Aggregated Spirals

**Figure numbering:** 21 = boltz2X-RAY, 22 = boltz2MD, 23 = protenix, 24 = rf3, 25 = protenix-mini, 26 = protenix-tiny, 28 = rf3-chiral-off.

**What is plotted:** Single-model aggregated spiral. Same aggregation as Figure 20 but for one model. Dots colored by guidance weight (flare colormap). IQR error bars. Colorbar legend for guidance weight.

**Data included:**
- All (model, guidance_weight) groups for the selected model with `n ≥ 3 proteins`.

**What is excluded:**
- Groups with fewer than 3 proteins.

---

### Figures 29–36 — Training-Distribution Bias

**Common inputs:**
- `lddt_assignment_counts.csv` — pre-cutoff (PDB depositions before 2021-09-30) per-selection counts of cluster members assigned to altloc A vs altloc B by lDDT comparison. Loaded via `load_training_cluster_counts` (lowercases `protein`, renames `selection_part` → `selection`, drops the multi-part `selection` column, computes `frac_a_cluster = n_A / (n_A + n_B)` and `x_one_minus_fa = 1 − frac_a_cluster`).
- `initial40results/Sampleworks_initial_structures.csv` — per-protein deposition date and `Pre`/`Post` cutoff flag. Loaded via `load_deposition_cutoffs` (dedupes `9B22`, fills `1VME, 4OLE, 5SOP, 6B8X` as `Pre`).
- `initial40results/lddt_results.csv` — predicted occupancies; loaded via `load_lddt` then narrowed to main models via `filter_by_model_group(df, "main")`.

**Common transforms:** `build_training_bias_df` parses the **second** float of each `occupancies` cell (numpy-style `"[0.375 0.625]"` → `0.625`) into a `b_occupancy` column, inner-joins lddt to the cluster frame on `(protein, selection)`, and left-joins the cutoff flag on `protein`.

**Models plotted:** boltz2X-RAY, boltz2MD, protenix, rf3.

#### Figure 29 — Training-Cluster vs Predicted State-B Fraction (Small Multiples)

**Output:** `figures/29_training_bias_small_multiples.svg`.

**What is plotted:** 4-column × 2-row scatter grid. Columns = main models, rows = condition (top: Baseline gw=0.0; bottom: Optimal guidance per `OPTIMAL_GW`). Shared axes (`xlim = ylim = [-0.02, 1.02]`). `y = 0.5` reference line drawn behind the points. Points are monochrome dark gray; alpha scales with `log(cluster_size + 1) / log(s_max + 1)` so larger-cluster cases dominate visually.

**Axes:** X = `frac_b_cluster`, Y = `b_occupancy`.

**Data included:** one point per `(protein, selection, model_name, condition)` with `s_cluster > 0`.

**What is excluded:** `s_cluster == 0` rows (these go to Figure 32 instead); all upstream filters.

#### Figure 34 — Calibration Overlay (Redesign A, 1×4)

**Output:** `figures/34_training_bias_calibration_overlay.svg` + legend at `figures/34_training_bias_calibration_overlay_legend.svg`.

**What is plotted:** 1-row × 4-column grid of panels, one per main model. Each panel pairs a scatter axis (left, width 4) with a marginal y-histogram strip (right, width 1) sharing y. Within each scatter, Baseline (gray, `#888888`) and Optimal guidance (model color from `MODEL_PALETTE`) points are overlaid. Reference lines `y = x` and `y = 0.5` are drawn behind the points. Per-condition binned conditional-mean curves use 10 equal-count quantile bins (or `n_bins=3, min_bin_count=2` for sparse subsets); each bin's marker sits at the within-bin mean of x. Slope of `b_occupancy` on `frac_b_cluster` and Pearson r are annotated for both conditions. Marginal axes overlay both conditions' `b_occupancy` distributions.

The main F34 cell calls the builder with `annotation_layout="split"`: the baseline slope/r block sits at top-left `(0.03, 0.97)` and the optimal-guidance slope/r block at bottom-right `(0.97, 0.03)`, so the two annotations never visually collide. The subset variants (Figures 30, 31) intentionally keep the default `"stacked"` layout (both annotations top-left).

Legend (in `make_training_bias_calibration_overlay_legend`): two-row vertical legend — "Baseline" with a single gray marker, "Optimal guidance" with four model-colored markers shown side-by-side via `HandlerTuple`.

#### Figure 35 — Calibration Small Multiples v2 (Redesign B, 2×4)

**Output:** `figures/35_training_bias_small_multiples_v2.svg`.

**What is plotted:** 2-row × 4-column grid mirroring Figure 29's geometry (top = Baseline, bottom = Optimal guidance, columns = main models). Each scatter panel is paired with a single-condition marginal y-histogram (per-panel, not per-condition like Figure 34). Every panel draws both reference lines, its own binned conditional-mean curve on top of a `mean ± std` shadow band, and annotates slope / Pearson r for the panel's own condition. Points monochrome (Figure 29 convention).

#### Figure 36 — Calibration Explainer (Empty Panel)

**Output:** `figures/36_training_bias_calibration_explainer.svg`.

**What is plotted:** A single 4×4 inch square-locked panel (`set_box_aspect(1)`) with the same axis style as the v2 small multiples and three lines: `y = x` (training prior, dashed), `y = 0.5` (guidance target, dashed), and a solid bisector at `y = 0.5*x + 0.25` (slope 0.5 through `(0.5, 0.5)`). Two italic black labels name the references — "training set prior" in the upper-right region rotated 45° to follow the `y = x` diagonal, and "perfect guidance" along the `y = 0.5` line left of center, unrotated. Two black double-headed arrows sit at the two ends of the solid bisector, perpendicular to it (purely visual, no text). No data; intended as a minimal explanatory companion to Figure 35 for downstream labeling.

#### Figures 30, 31 — Calibration Overlay on OOD Subsets

**Output filenames:**
- `figures/30_training_bias_calibration_overlay_low_cluster.svg`
- `figures/31_training_bias_calibration_overlay_post_cutoff.svg`

**What is plotted:** Both use the Figure 34 layout (`make_training_bias_calibration_overlay`). The only difference from Figure 34 is the subset filter applied to both conditions before plotting:

- **Figure 30 (low-cluster subset):** `max_cluster_size=5` — keeps only rows with `n_cluster_evaluated ≤ 5`. With only ~9 surviving points per panel, the binned mean uses `n_bins=3, min_bin_count=2`.
- **Figure 31 (post-cutoff subset):** `cutoff_status="Post"` — keeps only rows deposited after 2021-09-30.

#### Figure 32 — Empty-Training-Cluster Strip Plot

**Output:** `figures/32_training_bias_no_cluster_strip.svg`.

**What is plotted:** 1×4 row of panels (one per main model). Each panel pairs a half-violin (Gaussian KDE, clipped to `[0, 1]`, Scott bandwidth) on the left with a jittered scatter on the right; the X axis carries no quantitative meaning. Baseline (`#888888`) and Optimal-guidance (`MODEL_PALETTE[model]`) are overlaid in both halves (KDEs share the left region; jittered points share the right region). The KDE is skipped for any condition with fewer than 3 points (or zero variance); jittered points are always drawn. A `y = 0.5` reference line is drawn behind the data. Spines are stripped to a single left spine on the leftmost panel only; model names appear below each panel as the x-axis label. Each panel's Optimal-condition sample size is annotated.

This figure exists because its rows have `frac_b_cluster = N_B / (N_A + N_B)` undefined (both numerator and denominator zero), so they cannot appear on the X-axis of Figures 29, 30, 31, 34, or 35.

**Data included:** rows where `n_assigned_to_A == 0` AND `n_assigned_to_B == 0`, for the four main models, at Baseline and Optimal guidance.

---

## Summary: What Is Universally Excluded

The following data **never appears in any plot**:

1. **Non-50/50 occupancy configurations** — only the equal-occupancy (A=0.5, B=0.5) conformer split is analyzed.
2. **Proteins or models missing from the source CSVs** — the analysis is limited to whatever was computed in `initial40results/`.
3. **Baseline rows in delta plots** — deltas are defined as guided minus baseline, so the baseline itself (gw=0.0) does not appear as a data point in any delta/spiral figure (it is the implicit origin).
4. **(protein, model) pairs missing a baseline measurement** — the `make_delta_frame` merge is an inner join, so if there is no gw=0.0 row for a (protein, model), that pair is absent from all delta-based figures.

## Summary: Per-Figure Exclusion Table

| Figure | Excludes "all altloc" selection | Outlier/range clipping | Model subset |
|---|---|---|---|
| 04_09 Combined Baseline | Yes (LDDT side) | None | By dropdown |
| 12–17, 19 Per-Model Spirals | N/A (geometry+RSCC) | None | One model each |
| 20 Aggregated Spiral All | N/A | n < 3 groups dropped | 4 main models |
| 21–26, 28 Aggregated Spirals | N/A | n < 3 groups dropped | One model each |
| RSCC Scatter (per model) | Yes | None | One model each |
| RSCC Scatter Aggregated | Yes | n < 3 main models dropped | Main only |
| RSCC vs LDDT (per model) | Yes (LDDT side) | None | One model each |
| Aggregated LDDT Scatter | Yes | n < 3 main models dropped | Main only |
| By-Class Combined Grid | Yes (small-segment branch) | None | Main only |
| 29 Training Bias Small Multiples | N/A (cluster-frame join) | s_cluster == 0 dropped | Main only |
| 30 Calibration Overlay Low Cluster | N/A | n_cluster_evaluated > 5 dropped | Main only |
| 31 Calibration Overlay Post Cutoff | N/A | cutoff_status != Post dropped | Main only |
| 32 No-Cluster Strip | N/A | A+B > 0 rows dropped | Main only |
| 34 Calibration Overlay | N/A | None | Main only |
| 34L Calibration Overlay Legend | N/A (no data) | N/A | Main only |
| 35 Calibration Small Multiples v2 | N/A | None | Main only |
| 36 Calibration Explainer | N/A (no data) | None | N/A (no data) |
