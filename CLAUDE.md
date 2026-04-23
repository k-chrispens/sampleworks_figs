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
- All data loading (`load_and_tag_path_csv`, `load_rscc`, `load_lddt`, `load_classification`)
- All plotting functions (`make_ridge_plot`, `make_rscc_scatter`, `make_spiral_figure`, `make_tradeoff_grid`, etc.)
- Constants: model palettes, display name mappings, guidance weight→strength bucketing, figure directory
- Data transformations: path metadata extraction, model name aliasing via `MODEL_NAME_ALIAS`, delta computation

**Data flow**: CSV files in `initial40results/` → loaded/filtered by helpers → plotted in notebooks → SVGs saved to `figures/`.

**Key data conventions**:
- All analysis uses `0.5occA_0.5occB` occupancy filter (equal-occupancy alternate conformations only)
- Raw model names are aliased at load time (e.g., `rf3` → `rf3-old`, `rf3-fixed` → `rf3`). See `MODEL_NAME_ALIAS` dict.
- Guidance weights are bucketed into strengths: 0.0=baseline, {0.005,0.01,0.1}=weak, {0.02,0.035,0.2,0.35}=medium, {0.05,0.5}=strong

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
