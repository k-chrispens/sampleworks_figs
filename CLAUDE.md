# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Analysis codebase for a structural biology paper evaluating AI protein structure prediction models (Boltz2, Protenix, RoseTTAFold3 variants) on alternate conformations with experimental density guidance. Compares model outputs against ground truth crystallographic data across metrics like RSCC, lDDT, clash scores, and bond geometry.

## Commands

Run notebooks interactively (the primary workflow):
```
uv run marimo edit 02_fit_quality.py
```

Run a notebook headlessly:
```
uv run marimo run 02_fit_quality.py
```

Run standalone scripts:
```
uv run python merge_new_rscc_lddt.py
```

Install dependencies:
```
uv sync
```

## Architecture

**Marimo notebooks** (`01_geometry.py`, `02_fit_quality.py`, `03_tradeoffs.py`) are the analysis entry points. Each is a marimo `App` with `@app.cell` decorated functions — not standalone scripts. They import everything from `analysis_helpers.py` and should be run via `marimo edit` or `marimo run`, not `python` directly.

**`analysis_helpers.py`** is the shared library. It contains:
- All data loading (`load_and_tag_path_csv`, `load_and_tag_path_csv_fast` (polars), `load_rscc`, `load_lddt`, `load_classification`, `merge_classification_by_containment`, `get_rscc_nan_proteins`, `load_deposition_cutoffs`, `load_training_cluster_counts`)
- All plotting functions:
  - Distributions: `make_ridge_plot`, `make_combined_baseline_distributions`
  - Per-model RSCC scatters: `make_rscc_scatter`, `make_scatter_comparison`
  - Classification-aware scatters: `make_classified_scatter`, `make_classified_scatter_grid`, `make_classified_aggregated_scatter`
  - Aggregated (mean ± SD over proteins): `make_aggregated_scatter`, `make_aggregated_spiral_figure`
  - Spiral / tradeoff plots: `make_spiral_figure`, `make_tradeoff_grid`, `annotate_tradeoff_highlights`, `make_conceptual_tradeoff_figure`
  - Training-distribution bias: `make_training_bias_scatter`, `make_training_bias_calibration_overlay` (supports `max_cluster_size` / `cutoff_status` subset filters and `n_bins` / `min_bin_count` overrides for sparse subsets), `make_training_bias_calibration_overlay_legend`, `make_training_bias_no_cluster_strip` (for the A=B=0 / undefined-X subset), `make_training_bias_small_multiples_v2`
- Save helpers: `save_fig`, `save_model_legend` (writes `figures/model_legend.svg`), `save_per_model_figures`
- Constants: `OCC_FILTER` (`"0.5occA_0.5occB"`), `MODEL_PALETTE` / `MODEL_DISPLAY` / `MODEL_NAME_ALIAS`, `GUIDANCE_PALETTE` / `GUIDANCE_MARKERS`, `GW_TO_STRENGTH`, `OPTIMAL_GW`, `CLASSIFICATION_PALETTE` / `CLASSIFICATION_DISPLAY` / `CLASSIFICATION_ORDER`, `NICE_XLABELS`, `MODELS_ALT` / `MODELS_PROTENIX_GROUP` / `MODELS_RF3_GROUP`, `FIGURE_DIR`, `TRAINING_CUTOFF_DATE`
- Data transformations: `extract_metadata_from_path`, `add_simplified_metadata_columns` (adds `guidance_strength` bucket), `tag_path_dataframe`, `make_delta_frame` (guided − baseline), `build_training_bias_df` (joins lddt + training-cluster + cutoff frames), `robust_z_score`

**Data flow**: CSV files in `initial40results/` (plus the root-level `lddt_assignment_counts.csv` for training-cluster bias) → loaded/filtered by helpers → plotted in notebooks → SVGs saved to `figures/`.

**Key data conventions**:
- All analysis uses the `OCC_FILTER` constant (`0.5occA_0.5occB`) — equal-occupancy alternate conformations only
- Raw model names are aliased at load time (e.g., `rf3` → `rf3-old`, `rf3-fixed` → `rf3`, `boltz2_X-RAY_DIFFRACTION` → `boltz2X-RAY`). See `MODEL_NAME_ALIAS`.
- Guidance weights are bucketed into strengths via `GW_TO_STRENGTH`: 0.0=baseline, {0.005,0.01,0.1}=weak, {0.02,0.035,0.2,0.35}=medium, {0.05,0.5}=strong
- Model groups for split plotting: `MODELS_ALT` (the four "alternate" model variants), `MODELS_PROTENIX_GROUP`, `MODELS_RF3_GROUP`. `filter_by_model_group(df, "main")` drops `MODELS_ALT`.
- Notebooks expose toggles for cumulative distributions (`cumulative` switch + `cumulative_style` radio: KDE vs ECDF) and for excluding proteins with NaN RSCC (`exclude_nan_rscc`, sourced via `get_rscc_nan_proteins`).
- Training-distribution sources: `lddt_assignment_counts.csv` holds per-selection counts of pre-cutoff PDB cluster members assigned to altloc A vs B by lDDT comparison (use `load_training_cluster_counts`); `initial40results/Sampleworks_initial_structures.csv` holds the per-protein deposition date and `Pre`/`Post` flag relative to `TRAINING_CUTOFF_DATE` (use `load_deposition_cutoffs`, which auto-fills the four pre-2018 PDBs missing from the CSV: `1VME, 4OLE, 5SOP, 6B8X`).

**Merge scripts** (`merge_combined_into_initial.py`, `merge_new_rscc_lddt.py`) are one-time data integration scripts that append new model results into `initial40results/`. They are not part of the regular analysis workflow.

**`data_in_each_plot.md`** documents exactly what data, filters, and transformations go into each figure — consult it when modifying plots.

**Key helpers for notebook code** (used to reduce duplication across cells):
- `filter_by_model_group(df, group, extra_query)` — filters by model group ("main", "protenix", "rf3") and returns `(filtered_df, filename_suffix)`
- `model_slug(model_name)` — returns filesystem-safe slug from MODEL_DISPLAY name
- `save_per_model_figures(df, model_col, make_fn, path_template)` — batch-saves one figure per model
- `OPTIMAL_GW` — dict mapping model names to their optimal guidance weights

## Important Conventions

- Do not duplicate data pipelines from the notebooks into standalone scripts — use the notebook's own execution via marimo.
- Figures output as SVG to `figures/` directory. The `save_fig` helper handles path construction.
- The codebase uses both pandas and polars: notebooks use pandas via the helpers; merge scripts use polars.
- Each notebook has a **figure catalog table** in its intro cell listing every figure it produces.

## Adding a New Figure

**Ridge plot with model-group variants:**
1. Add a cell with `filter_by_model_group()` + `make_ridge_plot()` + `save_fig()`
2. Update the figure table in the notebook's intro cell

**Per-model scatter:**
1. Build the dataframe (use `OPTIMAL_GW` for optimal guidance extraction)
2. Use `save_per_model_figures()` for batch save
3. Show selected model with the shared `scatter_model` dropdown

**Aggregated scatter:**
1. Group by protein, aggregate with mean +/- SD
2. Call `make_aggregated_scatter()` with the aggregated df

**Classification-aware scatter (side-chain / small-loop / large-loop / domain-shift):**
1. Merge metric df with `load_classification("small-altlocs-classified-rmsd.csv")` (use `merge_classification_by_containment` when the selection-name match needs to be by residue range containment)
2. Per-model: `make_classified_scatter()`; aggregated (mean ± SD across proteins): `make_classified_aggregated_scatter()`; combined model grid: `make_classified_scatter_grid()`

**Training-distribution bias (predicted occupancy vs. pre-cutoff cluster composition):**
1. Build the tidy frame with `build_training_bias_df(lddt_df, cluster_df, cutoffs_df, model_to_gw)` — pass `{m: 0.0 for m in MAIN_MODELS}` for the baseline figure, `{m: OPTIMAL_GW[m] for m in MAIN_MODELS}` for the optimal-guidance figure
2. Plot 1 (full scatter): `make_training_bias_scatter(df)` — X = `1 − fraction A in pre-cutoff training cluster`, Y = predicted B occupancy, all models color-coded; reference lines `y = x` (training prior) and `y = 0.5` (density-only)
3. Plot 2 (calibration overlay, Figure 34): `make_training_bias_calibration_overlay(baseline_df, optimal_df)` — 1×4 row of panels, one per main model. Baseline (gray) and Optimal guidance (dark) overlaid in the same axes; both reference lines (y = x training prior, y = 0.5 guidance target — both drawn as light dashed lines); per-condition binned conditional-mean curves; guidance-only slope/r annotation; per-panel two-condition marginal y-histogram on the right edge. The legend is rendered separately by `make_training_bias_calibration_overlay_legend()` so it can be arranged independently in the final layout. Pass `max_cluster_size=5` (Figure 30, plus `n_bins=3, min_bin_count=2` for the 9-point subset) or `cutoff_status="Post"` (Figure 31) to restrict both conditions to an out-of-distribution subset before plotting.
4. Plot 3 (empty-cluster strip, Figure 32): `make_training_bias_no_cluster_strip(baseline_df, optimal_df)` — 1×4 strip plot for selections where `n_assigned_to_A == 0` AND `n_assigned_to_B == 0` (X = `frac_b_cluster` is undefined; these rows cannot appear on Figure 34). Per-panel jittered strip overlaying Baseline (gray) and Optimal guidance (model color), `y = 0.5` reference line, marginal y-histogram on the right edge, and an `n = …` annotation.
5. Plot 4 (calibration small multiples v2, Figure 35): `make_training_bias_small_multiples_v2(baseline_df, optimal_df)` — 2×4 mirroring Figure 29 with per-panel conditional mean, per-panel slope/r annotation, per-panel single-condition marginal y-histogram, and an inset legend in the top-left panel. Use this for the "every panel earns its place" rhetorical comparison against Figure 29.

## Validation
- For quick syntax/import checks while editing helpers, prefer `uv run python -c "import analysis_helpers"` over running the notebook headlessly. `uv run marimo run …` is slow.
