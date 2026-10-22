# CLAUDE.md

Guidance for Claude Code (claude.ai/code) working in this repository.

## Project Overview

Analysis codebase for a structural biology paper evaluating AI protein structure prediction models (Boltz2, Protenix, RoseTTAFold3 variants) on alternate conformations with experimental density guidance. Compares model outputs against ground-truth crystallographic data across RSCC, lDDT, clashscores, bond geometry, and minimum RMSD to each altloc.

The full analysis is three marimo notebooks (`02_fit_quality.py`, `03_tradeoffs.py`, `04_recycles_steps_ensembles.py`) on top of flat pure-Python modules. Figures save as SVG into `figures/`. Focused reconstructions are `01_baseline.py` (Figure 2C, supported by `baseline.py`), `02_guidance.py` (Figure 3B, supported by `guidance.py`), `03_rscc_improvement.py` / `04_clustering_improvement.py` (Figures 3C/D, supported by `region_improvement.py`), and `05_min_rmsd.py` (Figure 3E, supported by `min_rmsd.py`).

## Commands

```sh
uv sync                                # install
uv run marimo edit 01_baseline.py      # first figure reconstruction + exploration
uv run marimo edit 02_guidance.py      # guidance/geometry reconstruction + exploration
uv run marimo edit 03_rscc_improvement.py  # region-level RSCC improvement + model pairs
uv run marimo edit 04_clustering_improvement.py  # clustering improvement + score components
uv run marimo edit 05_min_rmsd.py       # ensemble minima + across-model mean/median
uv run marimo edit 02_fit_quality.py   # primary workflow (interactive)
uv run marimo run  02_fit_quality.py   # headless render (slow — usually overkill)
```

Notebooks **must** be invoked via marimo, not `python file.py` — they are `marimo.App`s of `@app.cell` functions, not scripts.

## Module Map

Each module imports only from those above it; no cycles.

| Module | Role | Most-used surface |
|---|---|---|
| `style.py` | rcParams, palettes, constants, save helpers | `MODEL_PALETTE`, `MODEL_DISPLAY`, `MODEL_NAME_ALIAS`, `MODELS_ALT`, `OCC_FILTER[_ALTLOC]`, `OPTIMAL_GW`, `GW_TO_STRENGTH`, `save_fig`, `save_per_model_figures`, `save_model_legend` |
| `loaders.py` | CSV → tidy DataFrame | `load_rscc`, `load_rscc_max_rmsd`, `load_rscc_recyc_diff`, `load_lddt`, `load_min_altloc_rmsd`, `load_classification`, `merge_classification_by_containment`, `load_training_cluster_counts`, `load_deposition_cutoffs`, `get_rscc_nan_proteins`, `add_simplified_metadata_columns` |
| `baseline.py` | Figure 2C snapshot preparation, ECDFs, rendering, provenance | `load_baseline_snapshot`, `prepare_baseline`, `baseline_ecdf`, `make_baseline_static`, `make_baseline_interactive`, `baseline_provenance` |
| `transforms.py` | DataFrame → DataFrame (pure) | `make_delta_frame`, `build_training_bias_df`, `filter_by_model_group`, `format_residue_range`, `model_slug` |
| `guidance.py` | Figure 3B snapshot, protein deltas, median/IQR, renderers | `load_guidance_snapshot`, `prepare_guidance`, `make_guidance_static`, `make_guidance_interactive` |
| `region_improvement.py` | Figures 3C/D paired model results, region medians, highlights/KDEs | `load_improvement_snapshot`, `prepare_improvement`, `prepare_density`, `highlighted_regions`, `make_improvement_static`, `make_improvement_interactive` |
| `min_rmsd.py` | Figure 3E ensemble minima, mean/median aggregation, two-panel rendering | `load_min_rmsd_snapshot`, `prepare_min_rmsd`, `highlighted_min_rmsd`, `make_min_rmsd_static`, `make_min_rmsd_interactive` |
| `plot_interaction.py` | Persistent Plotly widget, stable plot/table selection, display-only outlines | `make_linked_widget`, `highlight_table_rows`, `highlight_observations`, `select_observations` |
| `plot_distributions.py` | Baseline KDE/ECDF | `make_combined_baseline_distributions` |
| `plot_scatter.py` | Scatter family | `make_scatter_comparison`, `make_rscc_scatter`, `make_aggregated_scatter`, `make_min_rmsd_scatter` |
| `plot_classified.py` | 4×4 model-by-class grid | `make_by_class_combined_grid` |
| `plot_tradeoff.py` | Spiral trajectories | `make_spiral_figure`, `make_aggregated_spiral_figure` |
| `plot_training_bias.py` | Figures 29–36 (training-cluster vs predicted occupancy) | see module docstrings + `data_in_each_plot.md` |

Read module docstrings for full APIs and design notes. **`data_in_each_plot.md` is the authoritative per-figure spec** — consult it when modifying any plot.

## Data Sources

All figure inputs live in **`results/`** (the old top-level + `initial40results/` CSVs are gone). The four formerly-monolithic sources (LDDT, base RSCC, clashscore, bond-violations) are now **five per-model-group files each, concatenated in-memory inline in the notebook cells** with hardcoded `results/…` paths — there is no central `RESULTS_DIR` constant and no `load_*_all()` helper.

**Focused reconstruction exception:** the five focused notebooks accept an explicit
results folder (default `results`, or `SAMPLEWORKS_RESULTS`) and read it without
writing. They fingerprint their source files and reuse existing loader declarations.
Protein-mean RSCC remains notebook 01's default; per-region RSCC is a labeled
alternative pending paper adjudication. Figures 2C/3B use the native legend for
model visibility; Figure 3C aggregates a fixed main-model population and its legend
controls display categories only, as do Figures 3D/E. Figure 3E defaults to the
existing across-model median and exposes the caption-described mean as an explicit
unverified alternative. Circled regions have no fixed segment text; their static
reconstruction ZIPs bundle the SVG, companion highlight table, and recipe.
Toolbar SVGs reflect the visible view; separate
reconstruction SVGs and CSVs have explicit fixed scopes. Both renderers consume the
same prepared coordinates (and, for Figure 3C, density curves).
Main observation/region tables are bidirectionally linked through a persistent
Plotly `FigureWidget` and `mo.ui.anywidget`. Table callbacks update marker styles
in place; do not recreate the chart on row selection (that loses zoom/legend state).
Black inspection outlines never enter fixed reconstruction exports.
Tests: `uv run python -m unittest discover -s tests -v`. See
`open_problems.md`; merged-artifact storage remains unapproved.

| Path (all under `results/`) | Contents |
|---|---|
| `rscc_results_max_rmsd_{boltz_protenix,rf3_fixed,protenix_mini,protenix_tiny}.csv` + `rscc_results_rf3_chiral_off.csv` | **Primary RSCC source for all non-classification figures** (max-RMSD altloc region). Concatenated in-memory by `load_rscc_max_rmsd` (paths in `RSCC_MAX_RMSD_SOURCES`); covers boltz2 MD/XRAY, protenix, rf3, protenix-mini/-tiny, rf3-chiral-off (last from its non-max CSV — no max-RMSD variant exists). Sparser gw sweep than base (boltz2/protenix 6 gws; rf3 7 gws). Also backs `get_rscc_nan_proteins` |
| `rscc_results_{boltz_protenix,rf3_fixed,protenix_mini,protenix_tiny,rf3_chiral_off}.csv` | Base per-altloc-region RSCC (all 10 gws). **Used only for the by-class / conformational-class grid** — concatenated inline in the nb-02 by-class cell via `load_rscc` (mini/-tiny/chiral-off via `model_name_override`), then region-containment classification. The lone domain-shift region (`2e7z`) is folded in from this same concat — no separate domain-shift CSV |
| `lddt_results_{boltz_protenix,rf3_fixed,protenix_mini,protenix_tiny,rf3_chiral_off}.csv` | lDDT silhouette + predicted occupancies. **Concatenated inline** in nb 02 + 03 via `load_lddt(..., model_name_override=…)` (mini/-tiny/chiral-off disambiguated); loader derives `clustering_score = min(predicted_occ) * avg_silhouette_to_ref / min(map_occ)` |
| `Sampleworks_initial_structures.csv` | Per-protein deposition date + Pre/Post cutoff flag (`load_deposition_cutoffs` default) |
| `{clashscore_metrics,bond_angle_violation_fractions,bond_length_violation_fractions}_{boltz_protenix,rf3_fixed,protenix_mini,protenix_tiny,rf3_chiral_off}.csv` | Geometry metrics for the spiral figures. **Concatenated inline** in nb 03; the path-tagged variant files (mini/-tiny/chiral-off) get a post-tag `model_name` override because their raw run-dir collapses to `protenix_X-RAY_DIFFRACTION` / `rf3_X-RAY_DIFFRACTION` (see Conventions) |
| `min_altloc_rmsd_results_{boltz_protenix,protenix_mini,protenix_tiny,rf3_fixed}.csv` | Per-(protein, model, gw, selection) minimum RMSD across the 8-member ensemble to each altloc — loaded one-CSV-per-call by `load_min_altloc_rmsd` and concatenated in nb 02 (no rf3-chiral-off CSV — that model is absent from min-RMSD figures) |
| `lddt_assignment_counts.csv` | Pre-cutoff PDB cluster members assigned to altloc A vs B by lDDT (`load_training_cluster_counts` default) |
| `rscc_results_5mhx_recyc_diff.csv` | **5MHX-only** recycles × diffusion-steps grid-search RSCC (recycles `{1,3,10}` × steps `{50,100,200,400}`, four main models, optimal gw). Loaded by `load_rscc_recyc_diff` (keeps `n_recycles`/`n_diffusion_steps`); backs nb 04 only |
| `rscc_results_6b8x_ensemble_size.csv` | **6B8X-only** ensemble-size sweep RSCC (sizes `{1,2,4,8,16}`, four main models, optimal gw, 32 selections). Loaded via plain `load_rscc` (it already keeps `ensemble_size`); single-source (every model already at its optimal gw, no default-point merge); backs nb 04 only |
| `classify_altloc_regions.csv` | Conformational-change classification per (protein, selection) by backbone lDDT (`side_chain_only` / `small_loop` / `large_loop` / `domain_shift`). Loaded via `load_classification(CLASSIFICATION_CSV)`; the single swap-point constant `CLASSIFICATION_CSV` lives in `style.py`. Replaced the retired `small-altlocs-classified-rmsd.csv` |

There are no on-disk merge scripts (the legacy `merge_combined_into_initial.py` / `merge_new_rscc_lddt.py` were removed). **Present but intentionally unplotted:** `results/bond_{length,angle}_outliers_*.csv` and `results/max_rmsd_within_rmsd_*.csv` have no loader/figure yet (`results/rscc_results_6b8x_ensemble_size.csv` is now wired into nb 04 — see its row above). `results/boltz2_msa_depth.tsv~` is kept at the user's request — never delete it during cleanup.

## Conventions That Bite

- **50/50 occupancy only.** Every loader filters to `OCC_FILTER == "0.5occA_0.5occB"` (or `OCC_FILTER_ALTLOC == "{'A': 0.5, 'B': 0.5}"`). Other splits never reach any figure.
- **Model name aliasing.** Raw CSV names go through `MODEL_NAME_ALIAS` (e.g. `boltz2_X-RAY_DIFFRACTION → boltz2X-RAY`, `rf3-fixed → rf3`). Always work with the canonical form downstream and look it up via `MODEL_DISPLAY[name]` for figure titles.
- **RSCC source split.** `load_rscc_max_rmsd` (the max-RMSD-region CSVs) is the primary RSCC frame for every figure *except* the by-class grid, which concatenates the five **base** (non-max) `results/rscc_results_*.csv` inline (it needs the base's full per-region selection coverage). All four conformational classes — including the lone `domain_shift` region (`2e7z`, ~44 region rows across the 7 models) — now flow through **one** `merge_classification_by_containment` path; the old `rscc_results_domain_shift.csv` whole-domain source is retired. The max-RMSD CSVs have a sparser gw sweep than base (boltz2/protenix lack `0.005–0.035`; rf3 lacks `0.2–0.5`), so the rf3 spiral shows only its low-gw arm.
- **Disambiguating shared `model` columns.** Several CSVs store the same raw `model` value — both protenix-mini and protenix-tiny CSVs store `model == "protenix"`, and the rf3 / rf3-chiral-off CSVs store `model == "rf3"`. Disambiguate **at load time** with the `model_name_override` arg: `load_rscc(...)`, `load_lddt(...)`, `load_min_altloc_rmsd(...)` all take it. `load_rscc_max_rmsd` applies the per-source overrides from `RSCC_MAX_RMSD_SOURCES` (in `loaders.py`). **The path-tagged geometry files are a separate trap:** `load_and_tag_path_csv` / `tag_path_dataframe` extract `model_name` from the raw run-dir, which is `protenix_X-RAY_DIFFRACTION` for **both** mini and tiny and `rf3_X-RAY_DIFFRACTION` for chiral-off — so the clashscore/bond concats in nb 03 set `model_name` **after** tagging (mini/-tiny/chiral-off), or mini+tiny silently merge into one model.
- **Model groups + the dropdown.** `MODELS_ALT = {"protenix-mini", "protenix-tiny", "rf3-chiral-off"}` are the alternative variants. `MODELS_PROTENIX_GROUP` / `MODELS_RF3_GROUP` define the two variant families. The `model_group` dropdown (Main / Protenix variants / RF3 variants) lives in both nb 02 and nb 03 and drives **only the per-model-colored figures**: baseline distributions + by-class grid (02), the Figure-20 overlay spiral + the training-bias suite (03). `filter_by_model_group(df, group)` returns `(filtered_df, suffix)`; non-main saves get the `_protenix_alts` / `_rf3_alts` suffix so they never overwrite Main outputs. The **median-aggregated scatters, printed summary tables, per-model spirals, and nb 04** stay main-models-only / per-`selected_model` and ignore the group dropdown.
- **Guidance buckets.** `GW_TO_STRENGTH` maps numeric gw → `{baseline, weak, medium, strong}`. `OPTIMAL_GW[model]` is the per-model "best" gw used by every "optimal vs baseline" figure (`0.1` for the three main models, `0.01` for RF3).
- **Pinned cross-figure highlights.** `6B8X 181-183` and `5MHX 162-164` are circled in blue (`#1565C0`) on every aggregated scatter (RSCC, LDDT, min-RMSD). Pin new always-circled selections in the same color so they read as the same family across figures.
- **`exclude_nan_rscc` toggle.** When on, the notebook blacklists proteins with any NaN RSCC across the four main models (via `get_rscc_nan_proteins`) and threads the blacklist through every loader output. Don't bypass it.
- **Pre-cutoff training cluster source.** `load_training_cluster_counts()` (default `results/lddt_assignment_counts.csv`) adds `frac_a_cluster` / `x_one_minus_fa`. `load_deposition_cutoffs()` (default `results/Sampleworks_initial_structures.csv`) auto-fills `1VME, 4OLE, 5SOP, 6B8X` (pre-2018, absent from the CSV) as Pre.
- **Training-bias grids size to the group.** The `make_training_bias_*` helpers take `model_order` and build their subplot grid with `n_cols = len(model_order)` (4 main, 3 protenix, 2 rf3) — width scales so the 4-model case is byte-identical. `make_training_bias_calibration_explainer` (Fig 36) has no data/models and is group-independent.
- **All pandas.** Analysis code is pandas throughout (the removed merge scripts were the only polars users).

## Marimo Gotchas

- **`pd.query("@VAR")` references are invisible to marimo's static analysis.** A name only used inside a query string (e.g. `df.query("model_name not in @MODELS_ALT")`) is *not* added to the cell's argument list. Don't add it manually — marimo rewrites the cell function on save and any manually-added arg without a detected reference is stripped, producing a runtime `TypeError`. The `@VAR` lookup still resolves because marimo injects all defined names into the cell's globals.
- **Use leading-underscore locals** (`_df`, `_fig`, …) inside cells so marimo doesn't promote them to inter-cell variables. Only return names you want downstream cells to consume.
- **Cell signatures auto-rewrite on save.** Stick to the marimo-generated arg list; if you need a name that's missed by static analysis (the `@VAR` case above) rely on the runtime globals injection.
- **Shared UI controls** (`auto_save`, `exclude_nan_rscc`, `scatter_model`, `save_all_scatter`, etc.) are produced by upstream cells and consumed by every figure cell that needs them. New figure cells should take them as parameters rather than introducing parallel toggles.

## Adding a New Figure

The dominant pattern is **per-model figure + aggregated companion**, used by RSCC, LDDT, and min-RMSD scatters. Recipe:

1. **Loader** (`loaders.py`): one function per CSV family. Apply the 50/50 filter, alias model names through `MODEL_NAME_ALIAS`, drop NaN rows in the metric of interest with a one-line `print` warning per CSV. Mirror `load_rscc` / `load_min_altloc_rmsd` for the docstring style and column ordering.
2. **Plot helper** (`plot_*.py`): one signature serving both per-model and aggregated cells. Caller passes `title` and (optionally) `extra_highlight_groups: dict[color, [labels]]` for ad-hoc highlight overlays on top of the pinned blue ones. Wrap `adjustText.adjust_text` in a defensive `try/except` so labels degrade gracefully on overlap. `make_min_rmsd_scatter` and `make_aggregated_scatter` are the reference shapes.
3. **Notebook cells** (`02_fit_quality.py`):
   - **Per-model frame cell:** inner-join with `OPTIMAL_GW`, restrict to `guidance_weight == optimal_gw`, drop `selection_name == "all altloc"`, attach `label = "<PDB> <residue range>"` via `format_residue_range`. Print an audit when (protein, selection) pairs are silently dropped.
   - **Per-model display cell:** drive from the shared `scatter_model` dropdown; render `mo.callout(..., kind="warn")` in place of the figure when the model has no data; iterate via `save_per_model_figures` when `save_all_scatter` is on.
   - **Aggregated cell:** restrict to `model_name not in MODELS_ALT`, group by `(protein, selection, label)`, take median, drop `n < 3`, save under a metric-specific subdirectory of `figures/`.
4. **Document it:** add a row to the figure-catalog table in the notebook's intro cell, and a new section in `data_in_each_plot.md` (data sources, filters, what's included / excluded).

For non-scatter families look at the existing implementations: `make_by_class_combined_grid` (models × 4-class grid from a single per-region concat, sized by `model_order`), `make_spiral_figure` / `make_aggregated_spiral_figure` (trajectories through guidance space; the aggregated one is group-aware), and the `make_training_bias_*` family (figures 29–36, group-aware grids, fully documented in `data_in_each_plot.md`).

## Summary Tables & Shared Frames (`02_fit_quality.py`)

Beyond the figures, the notebook ends with a block of **printed summary tables** (rendered as marimo tables, no saved output) comparing baseline vs each model's optimal gw at the per-(protein, selection, model) level. They sit on shared infrastructure introduced to kill duplication — reuse it instead of re-deriving:

- **Helpers cell** (top of the notebook): `optimal_gw_df` (the `OPTIMAL_GW` frame to merge on), `baseline_vs_optimal(region, value_cols)` (pairs each metric's `gw == 0.0` value with its optimal-gw value → `<col>_baseline` / `<col>_optimal`), `add_label(df)` (the `"<PDB> <range>"` label column), `per_model_counts(df)` (the `"Model: n; …"` summary string).
- **Two delta frames**, each per (protein, selection, model): `rscc_scatter_df` (RSCC; also feeds the RSCC scatters) and `cs_scatter_df` (`clustering_score` + `avg_silhouette_to_ref` + `min_predicted_occupancy`, each baseline/optimal, plus `delta_clustering_score`). `rscc_cs_delta_df` is their inner join (RSCC ∩ clustering coverage); the trade-off / collapse / reverse tables filter it.
- **Column names follow `baseline_vs_optimal`'s scheme** (`<metric>_baseline` / `<metric>_optimal`, e.g. `clustering_score_optimal`, `min_rmsd_to_A_optimal`) — longer but self-describing.
- The tables themselves (threshold fractions, single-residue enrichment, RSCC↔clustering trade-offs + collapse, RSCC vs min-RMSD) are documented in `data_in_each_plot.md`.

## Validation

- **Fast import smoke test** — use before launching marimo:
  ```sh
  uv run python -c "import style, loaders, transforms, plot_distributions, plot_scatter, plot_classified, plot_tradeoff, plot_training_bias"
  ```
- **End-to-end cell execution without the UI** — much faster than `marimo run` and catches cell-signature `TypeError`s and runtime errors immediately:
  ```sh
  uv run marimo export script 02_fit_quality.py > _smoke.py && uv run python _smoke.py; rm _smoke.py
  ```
  Write the temp script into the project directory (not `/tmp`) so relative imports resolve.
- **Single-helper smoke test with real data** — prefer an inline `uv run python -c "..."` that loads the CSV, builds the frame, calls the new helper, and writes a PNG to `/tmp/`. Faster than the notebook when iterating on a plot function.

## Reference

- **`data_in_each_plot.md`** — exhaustive per-figure spec. Update it whenever you add or modify a figure.
- **Notebook intro cells** — figure catalog table; every figure produced by the notebook must be listed there with its output path.
