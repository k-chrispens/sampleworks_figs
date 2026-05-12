# Figures Cleanup & Helper Refactor

## Refined prompt

> I want to clean up this repository.
>
> * **Delete deprecated figures and their generators**:
>   * `figures/01_*` through `figures/11_*` SVGs (including `_protenix_alts`/`_rf3_alts` variants) **except** `figures/04_09_baseline_distributions.svg` (keep that one).
>   * The notebook cells in `01_geometry.py` and `02_fit_quality.py` that produced figs 01–11.
>   * `01_geometry.py` entirely (every figure it produces is being removed).
>   * `figures/rscc_scatterplots/rscc_vs_lddt_aggregated_all_models.svg` and the cell that generates it (the per-protein `cross_metric_df` aggregated panel).
>
> * **Fully purge `rf3-old`**:
>   * Drop from `MODELS_ALT` and `MODELS_RF3_GROUP`; remove any `rf3-old` entries from `MODEL_DISPLAY`, `MODEL_PALETTE`, `OPTIMAL_GW`, and any spiral figure-number maps.
>   * Delete `figures/18_spiral_rf3_old.svg`, `figures/27_aggregated_spiral_rf3_old.svg`, `figures/rscc_scatterplots/rscc_by_class_rf3_old.svg`, `figures/rscc_scatterplots/rscc_scatter_rf3_old.svg`.
>   * The live data is already clean — `MODEL_NAME_ALIAS` has no `rf3 → rf3-old` mapping, and no current CSV produces an `rf3-old` row after load + alias.
>   * Update the stale alias documentation in `data_in_each_plot.md`.
>
> * **Replace the per-model by-class figures with one combined figure**:
>   * `4 main_models × {side_chain_only, small_loop, large_loop, domain_shift}` grid with shared X and Y axes.
>   * Y range pinned across all rows and columns: lower = `min(0, global_min_residual) - 0.02` taken across the **full pooled dataset** (all 4 models and all 4 classes), upper = `1.0`. X = `[0, 1]`.
>   * Each panel keeps both reference lines: `y = 0` (no change) and `y = 1 − x` (max possible improvement ceiling).
>   * Each panel keeps a light-gray context shadow, but the shadow uses **that row's per-model pooled points** (model-specific context), not the full pooled cloud.
>   * Horizontal model banner above each row with proper spacing so banners do not collide with the next row's class titles.
>   * Saved to `figures/rscc_scatterplots/rscc_by_class_combined.svg`. Replaces all 13 existing variants (6 per-model small-segment, 1 aggregated, 1 grid, 5 domain-shift variants).
>   * Validate by rendering temp PNGs during development and visually checking for text overlap.
>   * **Extract whole-domain RSCC source**: extract the `domain_shift` whole-domain rows from the legacy `initial40results/rscc_results.csv.bak` into a new clean `initial40results/rscc_results_domain_shift.csv`. Update the by-class cell to load from this clean file instead of the `.bak`. The `.bak` file itself is left on disk (scope confirmation 3a).
>
> * **PNG-render audit on kept figures**:
>   * Render PNGs of every kept figure (spirals 12–17/19/20–26/28, training-bias 29–36, 04_09 baseline, scatters in `rscc_scatterplots/` and `cross_metric_scatterplots/`).
>   * Fix any glaring text-overlap regressions with minimal adjustments (margins, spacing, font sizes); no broader styling rework.
>
> * **Split `analysis_helpers.py` (4149 lines) into 8 flat top-level modules**:
>   * `style.py` — matplotlib rcParams, colors/palettes, `MODEL_DISPLAY`/`MODEL_PALETTE`/`MODELS_*` groups, `CLASSIFICATION_*` constants, `OPTIMAL_GW`, `OCC_FILTER*`, `TRAINING_CUTOFF_DATE`, `NICE_XLABELS`, `GW_TO_STRENGTH`, `MODEL_NAME_ALIAS`, plus `save_fig`, `save_model_legend`, `save_per_model_figures`.
>   * `loaders.py` — `load_rscc`, `load_lddt`, `load_classification`, `load_and_tag_path_csv`, `load_and_tag_path_csv_fast`, `load_deposition_cutoffs`, `load_training_cluster_counts`, `extract_metadata_from_path`, `tag_path_dataframe`, `add_simplified_metadata_columns`, `get_rscc_nan_proteins`, `merge_classification_by_containment`, `_parse_b_occupancy`.
>   * `transforms.py` — `make_delta_frame`, `build_training_bias_df`, `filter_by_model_group`, `model_slug`, `format_residue_range`, `robust_z_score`.
>   * `plot_distributions.py` — `make_combined_baseline_distributions` and its `_make_combined_hist_plot` internal helper.
>   * `plot_scatter.py` — `make_scatter_comparison`, `make_rscc_scatter`, `make_aggregated_scatter`.
>   * `plot_classified.py` — the new `make_by_class_combined_grid` function and helpers.
>   * `plot_tradeoff.py` — `make_spiral_figure`, `make_aggregated_spiral_figure`, `make_tradeoff_grid`, `annotate_tradeoff_highlights`, `make_conceptual_tradeoff_figure`, `add_tradeoff_guidance_regions`, `compute_conceptual_tradeoff_limits`.
>   * `plot_training_bias.py` — `make_training_bias_scatter`, `make_training_bias_calibration_overlay`, `make_training_bias_calibration_overlay_legend`, `make_training_bias_no_cluster_strip`, `make_training_bias_small_multiples_v2`, `make_training_bias_calibration_explainer`, and `_training_bias_*` private helpers.
>   * Delete `analysis_helpers.py` after the split; update `02_fit_quality.py` and `03_tradeoffs.py` to import from the new modules directly.
>   * After moving code, aggressively purge any function that ends up with zero callers across the surviving notebooks (`02_fit_quality.py` + `03_tradeoffs.py`); iterate to fixed-point so cascading dead helpers (private functions only called by other deleted helpers) are also removed.
>
> * **Phasing**:
>   * **Phase 1 (split first)**: create the 8 new modules, move code, rewire imports, delete `analysis_helpers.py`. One commit: `refactor: split analysis_helpers into flat modules`.
>   * **Phase 2 (deletions)**: figs 01–11 SVGs + their cells; `01_geometry.py`; `rscc_vs_lddt_aggregated_all_models.svg` + cell; rf3-old purge across constants and SVGs; cascading dead-code purge. One commit: `chore: remove deprecated figures 01-11 and rf3-old`.
>   * **Phase 3 (combined by-class figure)**: replace per-model/aggregated/grid/domain-shift by-class helpers with one combined function in `plot_classified.py`; update `02_fit_quality.py`; delete the 13 existing by-class SVGs; render+inspect temp PNG. One commit: `feat: combined by-class figure with domain shift`.
>   * **Phase 4 (audit)**: render PNGs of kept figures and fix glaring overlaps. One commit: `fix: text-overlap regressions on kept figures`.
>   * **Phase 5 (docs)**: update `data_in_each_plot.md`, `CLAUDE.md`, and the figure-catalog markdown in `02_fit_quality.py`'s intro cell. One commit: `docs: update plot catalog and conventions`.
>
> * **Scope**: only `figures/` subfolder, helper code, the two surviving notebooks (`02_fit_quality.py`, `03_tradeoffs.py`), and the two markdown docs. Leave repo-root paper composites (`fig1.svg`/`.af`, `fig2.svg`, …, `fig5new.af`, `*.jpg`, `*.png`) and merge scripts (`merge_*.py`) alone. Leave intermediate root-level data CSVs alone.

---

## Overview

- The codebase has accumulated deprecated figures: the entire 01–11 series in `figures/` was superseded by the combined `04_09_baseline_distributions.svg`, and the by-class scatter family has 13 overlapping variants where one grid would do.
- `rf3-old` is a phantom model: it appears in `MODELS_ALT`/`MODELS_RF3_GROUP` and produces 4 stale SVGs, but no live data path produces an `rf3-old` model name after alias. Purging it removes dead branches in the per-model batch loops and the variant dropdown.
- `analysis_helpers.py` is 4149 lines in one file. Splitting into 8 flat modules along natural seams (constants/IO, data loading, data transforms, plot families) makes individual plots locate-and-tweak in seconds rather than requiring a scroll.
- The by-class scatter figures have the worst text-overlap issues in the kept figure set — axis labels collide with tick labels, per-axis xlabels run together at column boundaries, and suptitle/row-banner/class-title text overlaps in the grid. Visible in PNG renders of the existing SVGs.
- The cleanup is phased so each phase produces a working state: split first (no behavior change), then delete (smaller code), then build the new combined figure (replaces 13 SVGs), then audit kept figures, then update docs.

## Expected behavior

- `uv run marimo edit 02_fit_quality.py` and `uv run marimo edit 03_tradeoffs.py` continue to launch; all surviving cells render their figures.
- Imports in the surviving notebooks come from new flat modules (`from style import ...`, `from loaders import load_rscc`, etc.); `analysis_helpers.py` no longer exists.
- `figures/` contains exactly `04_09_baseline_distributions.svg` from the 01–11 range; figs 12–17, 19, 20–26, 28, 29–36 are unchanged in content; `figures/model_legend.svg` is regenerated without `rf3-old`.
- `figures/rscc_scatterplots/` contains `rscc_by_class_combined.svg` (new), the per-model and aggregated `rscc_scatter_*.svg` (kept), `silhouette_scatter_aggregated_all_models.svg` (kept), and no by-class or rf3-old artifacts.
- `figures/cross_metric_scatterplots/` keeps the per-model `rscc_vs_lddt_{model}.svg` files; no aggregated cross-metric panel.
- The new combined by-class figure shows 4 model rows × 4 class columns; Y range identical in every panel and pinned to the global residual extent; `y = 0` and `y = 1 − x` reference lines visible in every panel; a model-specific gray context shadow under each row; no text collisions between suptitle, model banners, class titles, axis ticks, and axis labels.
- The domain-shift column in the combined by-class figure reads from a new clean `initial40results/rscc_results_domain_shift.csv` (extracted from the legacy `.bak`), so the production code path no longer references a `.bak` file.
- The `Model group` dropdown in `02_fit_quality.py` is kept — `fig_04_09` still produces protenix/rf3 baseline-distribution variants via `filter_by_model_group`, and the dropdown drives that switch.
- The `exclude_nan_rscc` switch and `get_rscc_nan_proteins` helper are kept; the cross-notebook NaN-protein blacklist remains available.
- `data_in_each_plot.md` reflects the post-cleanup figure inventory; `CLAUDE.md` references the new module layout; the figure-catalog markdown table in `02_fit_quality.py`'s intro cell matches what the notebook still produces.
- The `auto_save` switch in `02_fit_quality.py` saves only the figures that survive the cleanup (combined baseline, per-model RSCC and RSCC-vs-LDDT scatters, per-model and aggregated RSCC scatters, the new combined by-class figure).
- `03_tradeoffs.py` figure outputs are unchanged except for the removal of figures 18 (rf3-old spiral) and 27 (rf3-old aggregated spiral) from its per-model loops.

## Implementation plan

### New modules (replacing `analysis_helpers.py`)

#### `style.py`
- **Constants moved verbatim from `analysis_helpers.py`**: `OCC_FILTER`, `OCC_FILTER_ALTLOC`, `TRAINING_CUTOFF_DATE`, `GW_TO_STRENGTH`, `OPTIMAL_GW`, `MODEL_PALETTE`, `MODEL_DISPLAY`, `MODEL_NAME_ALIAS`, `MODELS_ALT`, `MODELS_PROTENIX_GROUP`, `MODELS_RF3_GROUP`, `GUIDANCE_PALETTE`, `GUIDANCE_MARKERS`, `CLASSIFICATION_PALETTE`, `CLASSIFICATION_DISPLAY`, `CLASSIFICATION_ORDER`, `NICE_XLABELS`, `FIGURE_DIR`.
- **rcParams block** at module top (current "Publication styling" section).
- **Functions moved**: `save_fig`, `save_model_legend`, `save_per_model_figures`.
- **Constant edits**: remove `"rf3-old"` from `MODELS_ALT`; remove `"rf3-old"` from `MODELS_RF3_GROUP`; remove any `MODEL_PALETTE["rf3-old"]`, `MODEL_DISPLAY["rf3-old"]`, `OPTIMAL_GW["rf3-old"]` entries.

#### `loaders.py`
- `extract_metadata_from_path`, `add_simplified_metadata_columns`, `tag_path_dataframe`, `load_and_tag_path_csv`, `load_and_tag_path_csv_fast`, `load_rscc`, `get_rscc_nan_proteins`, `_parse_b_occupancy`, `load_lddt`, `load_classification`, `merge_classification_by_containment`, `load_deposition_cutoffs`, `load_training_cluster_counts`.
- Imports from `style` for `MODEL_NAME_ALIAS`, `OCC_FILTER`, `OCC_FILTER_ALTLOC`, `GW_TO_STRENGTH`, `TRAINING_CUTOFF_DATE`.

#### `transforms.py`
- `make_delta_frame`, `build_training_bias_df`, `filter_by_model_group`, `model_slug`, `format_residue_range`, `robust_z_score`.
- Imports from `style` for `MODELS_ALT`, `MODELS_PROTENIX_GROUP`, `MODELS_RF3_GROUP`, `MODEL_DISPLAY`.

#### `plot_distributions.py`
- `make_combined_baseline_distributions` and its `_make_combined_hist_plot` helper.
- Imports from `style` for palettes and rcParams influence.

#### `plot_scatter.py`
- `make_scatter_comparison`, `make_rscc_scatter`, `make_aggregated_scatter`.

#### `plot_classified.py` (rewritten)
- **Delete**: `make_classified_scatter`, `make_classified_scatter_grid`, `make_classified_aggregated_scatter` (all 3 existing variants).
- **New function**: `make_by_class_combined_grid(merged_df, *, model_order=None) -> plt.Figure`.
  - Input: single tidy frame with columns `protein, selection, model_name, rscc_baseline, rscc_optimal, classification`. The caller concatenates per-altloc-region rows (small-segment classes) and whole-domain rows (`domain_shift`) upstream; the helper is presentation-only.
  - Derives `residual = rscc_optimal − rscc_baseline`.
  - Iterates model rows (default `[boltz2X-RAY, boltz2MD, protenix, rf3]` from `MODEL_DISPLAY`, restricted to main models present in the data) × class columns (`CLASSIFICATION_ORDER`).
  - `sharex=True, sharey=True` so all 16 panels share `xlim=[0,1]`, `ylim=[min(0, global_min_residual)-0.02, 1.0]`.
  - Per-panel light-gray shadow uses `(rscc_baseline, residual)` of every row in **that model's row** (not the full pooled cloud).
  - Per-panel class points overlaid on top, colored by `CLASSIFICATION_PALETTE[cls]`.
  - Reference lines drawn behind points: `y = 0` (no-change) and `y = 1 − x` (max-improvement ceiling).
  - Column titles (top row only) show `CLASSIFICATION_DISPLAY[cls]` plus `n = unique (protein, selection) count for that class across all models`.
  - Model row label: horizontal banner above each row via `fig.text(...)`. Banner Y position placed above the row's top edge with `hspace` enlarged so the banner sits in clear whitespace between rows.
  - Figure-level `supxlabel("Baseline RSCC")` and `supylabel("Improvement (Guided − Baseline RSCC)")` only; no per-axis xlabels (the source of the run-on text).
  - `subplots_adjust(left=0.09, right=0.99, top=0.92, bottom=0.08, hspace=0.55, wspace=0.08)` — wider `hspace` than the current `0.45` so banners fit cleanly; narrower `wspace` keeps columns tight.
- **Helpers** (private): `_compute_global_residual_range(df)`, `_draw_class_panel(ax, x_bg, y_bg, x_pts, y_pts, color)`.

#### `plot_tradeoff.py`
- `make_spiral_figure`, `make_aggregated_spiral_figure`, `make_tradeoff_grid`, `annotate_tradeoff_highlights`, `make_conceptual_tradeoff_figure`, `add_tradeoff_guidance_regions`, `compute_conceptual_tradeoff_limits`.

#### `plot_training_bias.py`
- `make_training_bias_scatter`, `make_training_bias_calibration_overlay`, `make_training_bias_no_cluster_strip`, `make_training_bias_calibration_overlay_legend`, `make_training_bias_small_multiples_v2`, `make_training_bias_calibration_explainer`, plus `_training_bias_*` helpers.

### Notebook updates

#### `02_fit_quality.py`
- Replace the single `from analysis_helpers import (...)` block with explicit imports from each new module.
- **Delete cells** that generate (per current line ranges): `fig_05` (~L240–263), `fig_06` (~L266–291), `fig_07`/`fig_08` (~L294–321), `fig_10` (~L789–821), `fig_11` (~L824–853).
- **Delete the `rscc_vs_lddt_aggregated_all_models` block** inside the current "Aggregated Scatter Plots" cell (~L747–770 — the `_cross_agg_df` block and its `mo.tabs` entry). Keep the LDDT scatter aggregated tab.
- **Replace the by-class cell** (~L506–682) with a much smaller block:
  - Build `_rscc_by_class_df` (small-segment data, exclude `domain_shift`) — keep the existing data-prep logic.
  - Build `_rscc_ds_df` (domain-shift data) by reading the new `initial40results/rscc_results_domain_shift.csv` (no longer `.csv.bak`).
  - `pd.concat([_rscc_by_class_df, _rscc_ds_df], ignore_index=True)` to produce a single merged frame.
  - Call `make_by_class_combined_grid(merged_df)` once.
  - Save (when `auto_save`) to `rscc_scatterplots/rscc_by_class_combined.svg`.
  - Drop the `save_all_scatter`-driven per-model batch and the per-model "show current model" branch for by-class.
- **Update the figure-catalog markdown table** in the intro cell to list only surviving figures.
- **`scatter_model` dropdown**: remove `"rf3-old"` from the `_options` dict.
- Decide whether `model_group` dropdown stays — see Open Questions.

#### `03_tradeoffs.py`
- Replace single `from analysis_helpers import (...)` with explicit imports.
- **Remove `rf3-old` entries** from the per-model figure-number maps: `"rf3-old": 18` in `_spiral_nums` and `"rf3-old": 27` in `_agg_spiral_nums`.
- The per-model loops iterate over `df["model_name"].unique()`; once `rf3-old` is gone from `MODELS_ALT`, this naturally stops producing figs 18/27.

### Figure deletions

Delete from `figures/`:
- `01_clashscore_distribution.svg`
- `02_delta_clashscore_distribution.svg`
- `03_delta_bond_angle_outlier_fraction_distribution.svg`
- `04_rscc_baseline_all_models.svg`
- `05_rscc_baseline_vs_weak_distribution.svg`
- `06_delta_rscc_distribution.svg`
- `07_delta_rscc_weak_all_models.svg`
- `08_delta_rscc_weak_protenix.svg`
- `09_lddt_to_ref_baseline_distribution.svg`
- `10_lddt_to_ref_all_guidance_distribution.svg`, `10_lddt_to_ref_all_guidance_distribution_protenix_alts.svg`, `10_lddt_to_ref_all_guidance_distribution_rf3_alts.svg`
- `11_delta_lddt_to_ref_distribution.svg`, `11_delta_lddt_to_ref_distribution_protenix_alts.svg`, `11_delta_lddt_to_ref_distribution_rf3_alts.svg`
- `18_spiral_rf3_old.svg`
- `27_aggregated_spiral_rf3_old.svg`

Delete from `figures/rscc_scatterplots/`:
- `rscc_by_class_aggregated.svg`
- `rscc_by_class_boltz2_md.svg`, `rscc_by_class_boltz2_xrd.svg`, `rscc_by_class_protenix.svg`, `rscc_by_class_protenix_mini.svg`, `rscc_by_class_protenix_tiny.svg`, `rscc_by_class_rf3.svg`, `rscc_by_class_rf3_chiral_off.svg`, `rscc_by_class_rf3_old.svg`
- `rscc_by_class_domain_shift_aggregated.svg`
- `rscc_by_class_domain_shift_boltz2_md.svg`, `rscc_by_class_domain_shift_boltz2_xrd.svg`, `rscc_by_class_domain_shift_protenix.svg`, `rscc_by_class_domain_shift_protenix_mini.svg`, `rscc_by_class_domain_shift_protenix_tiny.svg`, `rscc_by_class_domain_shift_rf3.svg`, `rscc_by_class_domain_shift_rf3_chiral_off.svg`
- `rscc_by_class_domain_shift_grid_all_models.svg`
- `rscc_scatter_rf3_old.svg`
- `rscc_vs_lddt_aggregated_all_models.svg`

Delete `01_geometry.py`.

### Dead-code purge (fixed-point iteration)

Re-running the purge until no further changes occur. Likely deletes (verify per pass):
- `make_ridge_plot` — only caller was figs 01/02/03/05/06/07/08/10/11.
- Possibly `_make_combined_hist_plot` if kept private to `make_combined_baseline_distributions` only.
- Possibly `get_rscc_nan_proteins` if both surviving notebooks drop the `exclude_nan_rscc` switch (verify; likely kept).
- Possibly `extract_metadata_from_path` if only called by `tag_path_dataframe` and `tag_path_dataframe` is unused — confirm.
- Final list determined empirically per pass.

### Docs

- `data_in_each_plot.md`: remove the "Figure 01" through "Figure 11" sections; remove the stale `MODEL_NAME_ALIAS` row claiming `rf3 → rf3-old`; replace the by-class scatter section with one section describing the combined figure; drop the `rscc_vs_lddt_aggregated_all_models.svg` mention.
- `CLAUDE.md`: update the helper-module listing to reflect the 8 flat modules; remove `make_ridge_plot` and other deleted functions from the "All plotting functions" enumeration; replace the "RSCC by classification" entries with the single combined figure entry; remove `01_geometry.py` from the "marimo notebooks" list.
- `02_fit_quality.py` intro markdown table: list only kept figures.

## Implementation phases

### Phase 1 — Split `analysis_helpers.py` into 8 flat modules

- Create `style.py`, `loaders.py`, `transforms.py`, `plot_distributions.py`, `plot_scatter.py`, `plot_classified.py`, `plot_tradeoff.py`, `plot_training_bias.py`. Move code from the corresponding sections of `analysis_helpers.py`; preserve function signatures and behavior.
- Update internal cross-imports so each module imports what it needs from `style`, `loaders`, etc.
- Update `02_fit_quality.py` and `03_tradeoffs.py` to import from the new modules instead of `analysis_helpers`.
- Delete `analysis_helpers.py`.
- Verify: `uv run python -c "import style, loaders, transforms, plot_distributions, plot_scatter, plot_classified, plot_tradeoff, plot_training_bias"` runs clean.
- Verify: `uv run python -c "from plot_classified import make_classified_scatter_grid"` (or other functions that still exist at this stage) runs clean.
- Commit: `refactor: split analysis_helpers into flat modules`.
- **Working state**: both notebooks still produce the same set of figures as before, just sourced from the new modules.

### Phase 2 — Delete deprecated figures, cells, file, and rf3-old; cascading dead-code purge

- Delete the SVGs listed in **Figure deletions** above (figs 01–11 except 04_09 and the four `rf3_old` variants and `rscc_vs_lddt_aggregated_all_models.svg`).
- Delete the corresponding cells in `02_fit_quality.py` (figs 05/06/07/08/10/11, the `rscc_vs_lddt_aggregated_all_models` block, and the `rf3-old` dropdown option).
- Delete `01_geometry.py`.
- Edit `style.py`: remove `"rf3-old"` from `MODELS_ALT` and `MODELS_RF3_GROUP`; remove any `MODEL_PALETTE["rf3-old"]`, `MODEL_DISPLAY["rf3-old"]`, `OPTIMAL_GW["rf3-old"]` entries.
- Edit `03_tradeoffs.py`: remove `"rf3-old": 18` and `"rf3-old": 27` from the figure-number maps.
- Regenerate `figures/model_legend.svg` by triggering `save_model_legend()` (manually via `auto_save` or a one-line invocation) so the legend no longer references `rf3-old`.
- Run fixed-point dead-code purge:
  1. `grep -rE '^def [a-zA-Z_]+' *.py` to enumerate all module functions.
  2. For each function, `grep -rE 'name\(|from [^ ]+ import [^,)]*\bname\b' *.py 02_fit_quality.py 03_tradeoffs.py`.
  3. Delete any function with zero callers outside its own definition.
  4. Repeat from step 1 until no further deletions.
- Verify: `uv run marimo edit 02_fit_quality.py` and `uv run marimo edit 03_tradeoffs.py` launch and every surviving cell executes.
- Commit: `chore: remove deprecated figures 01-11 and rf3-old`.
- **Working state**: every surviving notebook cell still runs; the by-class cell still produces the old per-model + aggregated + grid SVGs (those get replaced in Phase 3).

### Phase 3 — Build the combined by-class figure

- **Extract domain-shift source CSV**: one-shot script (or `uv run python -c "..."`) that reads `initial40results/rscc_results.csv.bak`, joins with the `domain_shift` rows of `small-altlocs-classified-rmsd.csv` by `(protein, selection)`, and writes the resulting frame to `initial40results/rscc_results_domain_shift.csv` (same column schema as `rscc_results.csv` so `load_rscc` can read it). The `.bak` file is left on disk untouched.
- Implement `make_by_class_combined_grid` in `plot_classified.py` per the spec above; remove `make_classified_scatter`, `make_classified_scatter_grid`, `make_classified_aggregated_scatter`.
- Update `02_fit_quality.py`'s by-class cell to load from the new `initial40results/rscc_results_domain_shift.csv` and call `make_by_class_combined_grid(...)` once.
- Render a temp PNG **directly via matplotlib** so fonts match the production SVG output: `fig.savefig('/tmp/by_class_combined.png', dpi=140)` from a tiny test harness. Avoid `inkscape` for this validation step — Inkscape may substitute fonts and produce overlap diagnostics that don't reflect the final SVG.
- Inspect the PNG for residual text overlap (Y-axis label vs ticks, suptitle vs banner, banner vs class title, x-tick run-on between columns). Tune `subplots_adjust` parameters and banner Y position until clean; commit only the final SVG (no PNG snapshot needs to be checked in).
- Delete the 13 by-class SVG variants enumerated in **Figure deletions**.
- Re-save the final SVG via the notebook (auto_save path) so the committed artifact comes from the production code path.
- Verify: `figures/rscc_scatterplots/rscc_by_class_combined.svg` exists and contains 16 panels (4 model rows × 4 class columns) with both reference lines, row banners, shared axes.
- Commit: `feat: combined by-class figure with domain shift`.

### Phase 4 — PNG-render audit of kept figures

- For each kept figure, re-execute its production code path and call `fig.savefig('/tmp/audit_<name>.png', dpi=120)` directly via matplotlib so the audit fonts match the SVG output exactly. Reserve `inkscape --export-type=png ...` only for SVGs whose generating code is not easily isolatable (worst case).
- Visually inspect each PNG for:
  - Tick labels overlapping figure-level supxy labels.
  - Suptitle overlapping the top row of subplots.
  - Per-axis labels colliding at column boundaries (the same family of issues seen in the by-class plots).
  - Legend collisions with data.
- For each overlap found, apply the minimal fix (adjust `subplots_adjust`, `tight_layout(pad=...)`, `bbox_to_anchor`, or font size).
- No new figure designs; no semantic changes.
- Re-render until clean; re-save final SVGs via the notebooks.
- Commit: `fix: text-overlap regressions on kept figures`.

### Phase 5 — Documentation update

- Edit `data_in_each_plot.md`: drop figs 01–11 sections, replace the by-class section, drop `rf3 → rf3-old` alias row, drop `rscc_vs_lddt_aggregated_all_models` mention, refresh the per-figure exclusion table.
- Edit `CLAUDE.md`: update the helpers section to list the 8 new modules; remove deleted-function names from the plotting-functions enumeration; replace by-class entries with the single combined figure; remove `01_geometry.py`.
- Edit the figure-catalog markdown table in `02_fit_quality.py`'s intro cell to match the post-cleanup state.
- Verify: every figure listed in the docs exists on disk; every helper named in `CLAUDE.md` exists in a module.
- Commit: `docs: update plot catalog and conventions`.

## Testing strategy

- **Phase 1 import smoke test**: `uv run python -c "import style, loaders, transforms, plot_distributions, plot_scatter, plot_classified, plot_tradeoff, plot_training_bias"`. No errors means the split is internally consistent.
- **Phase 1 notebook smoke test**: launch `uv run marimo edit 02_fit_quality.py` and `uv run marimo edit 03_tradeoffs.py`; confirm every cell runs without ImportError.
- **Phase 2 cell-by-cell test**: re-run each surviving cell in both notebooks; confirm the expected SVGs are produced and no cell raises.
- **Phase 2 dead-code purge correctness**: after each deletion pass, re-run the notebook smoke tests. If anything breaks, the function was not actually unreferenced — restore and add a unit-callsite check.
- **Phase 3 combined figure visual test**:
  - Render PNG of the new combined SVG; eyeball for (a) all 4 rows present with banners, (b) all 4 class columns present, (c) shared axes (every panel shows the same Y range), (d) no text collisions, (e) the domain-shift column shows just the n=1 point.
  - Spot-check data correctness: for one (model, class) cell, compare a few points against the old per-model SVG to confirm the same residuals.
- **Phase 4 audit**: render-then-eyeball loop for every kept SVG. Acceptance: no labels overlap data, no titles overlap subplots, no legend overlaps data.
- **Phase 5 docs sanity**: `ls figures/` and `ls figures/*/` to confirm every figure named in the docs exists; `grep -l '<function_name>' *.py` to confirm every helper named in `CLAUDE.md` exists.
- **Cross-phase regression**: after all phases, run `uv run marimo run 02_fit_quality.py` and `uv run marimo run 03_tradeoffs.py` once headlessly (acknowledging the user's memory says this is slow — use only as a final acceptance check, not for iterative validation).

## Open questions

- **Domain-shift CSV schema** — `rscc_results.csv.bak` has the same column set as `rscc_results.csv` so `load_rscc` should accept the extracted file without changes. Verify during Phase 3 that no column is silently dropped by the existing load path.

---

The plan is ready. Path:

```
blueprint/figures-cleanup-refactor/plan-figures-cleanup-refactor.md
```

You can chat to refine the plan freely, or ask "what are the open questions?" to get a list of the issues to resolve.
