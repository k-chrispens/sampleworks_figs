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
  - Training-distribution bias: `make_training_bias_scatter`, `make_training_bias_subset_grid`
- Save helpers: `save_fig`, `save_model_legend` (writes `figures/model_legend.svg`), `save_per_model_figures`
- Constants: `OCC_FILTER` (`"0.5occA_0.5occB"`), `MODEL_PALETTE` / `MODEL_DISPLAY` / `MODEL_NAME_ALIAS`, `GUIDANCE_PALETTE` / `GUIDANCE_MARKERS`, `GW_TO_STRENGTH`, `OPTIMAL_GW`, `CLASSIFICATION_PALETTE` / `CLASSIFICATION_DISPLAY` / `CLASSIFICATION_ORDER`, `NICE_XLABELS`, `MODELS_ALT` / `MODELS_PROTENIX_GROUP` / `MODELS_RF3_GROUP`, `FIGURE_DIR`, `TRAINING_CUTOFF_DATE`
- Data transformations: `extract_metadata_from_path`, `add_simplified_metadata_columns` (adds `guidance_strength` bucket), `tag_path_dataframe`, `make_delta_frame` (guided − baseline), `build_training_bias_df` (joins lddt + training-cluster + cutoff frames), `robust_z_score`

**Data flow**: CSV files in `initial40results/` (plus the root-level `max_rmsd_within_rmsd_2021-09-30_1A.csv` for training-cluster bias) → loaded/filtered by helpers → plotted in notebooks → SVGs saved to `figures/`.

**Key data conventions**:
- All analysis uses the `OCC_FILTER` constant (`0.5occA_0.5occB`) — equal-occupancy alternate conformations only
- Raw model names are aliased at load time (e.g., `rf3` → `rf3-old`, `rf3-fixed` → `rf3`, `boltz2_X-RAY_DIFFRACTION` → `boltz2X-RAY`). See `MODEL_NAME_ALIAS`.
- Guidance weights are bucketed into strengths via `GW_TO_STRENGTH`: 0.0=baseline, {0.005,0.01,0.1}=weak, {0.02,0.035,0.2,0.35}=medium, {0.05,0.5}=strong
- Model groups for split plotting: `MODELS_ALT` (the four "alternate" model variants), `MODELS_PROTENIX_GROUP`, `MODELS_RF3_GROUP`. `filter_by_model_group(df, "main")` drops `MODELS_ALT`.
- Notebooks expose toggles for cumulative distributions (`cumulative` switch + `cumulative_style` radio: KDE vs ECDF) and for excluding proteins with NaN RSCC (`exclude_nan_rscc`, sourced via `get_rscc_nan_proteins`).
- Training-distribution sources: `max_rmsd_within_rmsd_2021-09-30_1A.csv` holds per-selection counts of pre-cutoff PDB depositions within 1Å RMSD of altloc A vs B (use `load_training_cluster_counts`); `initial40results/Sampleworks_initial_structures.csv` holds the per-protein deposition date and `Pre`/`Post` flag relative to `TRAINING_CUTOFF_DATE` (use `load_deposition_cutoffs`, which auto-fills the four pre-2018 PDBs missing from the CSV: `1VME, 4OLE, 5SOP, 6B8X`).

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
3. Plot 2 (out-of-distribution 2x2 grid): `make_training_bias_subset_grid(baseline_df, optimal_df, model_name)` — per-model 2x2 strip plot. Rows: post-cutoff depositions (no `n_cluster_evaluated` filter) and pre-cutoff selections with `n_cluster_evaluated <= 5`. Columns: Baseline and Optimal guidance.

## Validation
- For quick syntax/import checks while editing helpers, prefer `uv run python -c "import analysis_helpers"` over running the notebook headlessly. `uv run marimo run …` is slow.
