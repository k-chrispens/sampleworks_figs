# CLAUDE.md

Guidance for Claude Code (claude.ai/code) working in this repository.

## Project Overview

Analysis codebase for a structural biology paper evaluating AI protein structure prediction models (Boltz2, Protenix, RoseTTAFold3 variants) on alternate conformations with experimental density guidance. Compares model outputs against ground-truth crystallographic data across RSCC, lDDT, clashscores, bond geometry, and minimum RMSD to each altloc.

The analysis is two marimo notebooks (`02_fit_quality.py`, `03_tradeoffs.py`) on top of eight flat pure-Python modules. Figures save as SVG into `figures/`.

## Commands

```sh
uv sync                                # install
uv run marimo edit 02_fit_quality.py   # primary workflow (interactive)
uv run marimo run  02_fit_quality.py   # headless render (slow — usually overkill)
uv run python merge_new_rscc_lddt.py   # one-off data integration scripts only
```

Notebooks **must** be invoked via marimo, not `python file.py` — they are `marimo.App`s of `@app.cell` functions, not scripts.

## Module Map

Each module imports only from those above it; no cycles.

| Module | Role | Most-used surface |
|---|---|---|
| `style.py` | rcParams, palettes, constants, save helpers | `MODEL_PALETTE`, `MODEL_DISPLAY`, `MODEL_NAME_ALIAS`, `MODELS_ALT`, `OCC_FILTER[_ALTLOC]`, `OPTIMAL_GW`, `GW_TO_STRENGTH`, `save_fig`, `save_per_model_figures`, `save_model_legend` |
| `loaders.py` | CSV → tidy DataFrame | `load_rscc`, `load_lddt`, `load_min_altloc_rmsd`, `load_classification`, `merge_classification_by_containment`, `load_training_cluster_counts`, `load_deposition_cutoffs`, `get_rscc_nan_proteins`, `add_simplified_metadata_columns` |
| `transforms.py` | DataFrame → DataFrame (pure) | `make_delta_frame`, `build_training_bias_df`, `filter_by_model_group`, `format_residue_range`, `model_slug` |
| `plot_distributions.py` | Baseline KDE/ECDF | `make_combined_baseline_distributions` |
| `plot_scatter.py` | Scatter family | `make_scatter_comparison`, `make_rscc_scatter`, `make_aggregated_scatter`, `make_min_rmsd_scatter` |
| `plot_classified.py` | 4×4 model-by-class grid | `make_by_class_combined_grid` |
| `plot_tradeoff.py` | Spiral trajectories | `make_spiral_figure`, `make_aggregated_spiral_figure` |
| `plot_training_bias.py` | Figures 29–36 (training-cluster vs predicted occupancy) | see module docstrings + `data_in_each_plot.md` |

Read module docstrings for full APIs and design notes. **`data_in_each_plot.md` is the authoritative per-figure spec** — consult it when modifying any plot.

## Data Sources

| Path | Contents |
|---|---|
| `initial40results/rscc_results.csv` | RSCC for main models, all 50/50 occupancies, all gws (per altloc-region) |
| `initial40results/rscc_results_domain_shift.csv` | Whole-domain RSCC for the domain-shift class only (clean extract from legacy `.bak`; `rf3-old` excluded) |
| `initial40results/lddt_results.csv` | lDDT silhouette + predicted occupancies (loader derives `clustering_score = min(predicted_occ) * avg_silhouette_to_ref / min(map_occ)`) |
| `initial40results/Sampleworks_initial_structures.csv` | Per-protein deposition date + Pre/Post cutoff flag |
| `initial40results/{clashscore_metrics,bond_*_violation_fractions}.csv` | Geometry metrics for the spiral figures |
| `min_altloc_rmsd_results_{boltz_protenix,protenix_mini,protenix_tiny,rf3_fixed}.csv` | Per-(protein, model, gw, selection) minimum RMSD across the 8-member ensemble to each altloc — loaded one-CSV-per-call by `load_min_altloc_rmsd` and concatenated in the notebook |
| `lddt_assignment_counts.csv` | Pre-cutoff PDB cluster members assigned to altloc A vs B by lDDT |
| `small-altlocs-classified-rmsd.csv` | Conformational-change classification per (protein, selection) |

`merge_combined_into_initial.py` and `merge_new_rscc_lddt.py` are one-shot integration scripts, not part of the analysis loop.

## Conventions That Bite

- **50/50 occupancy only.** Every loader filters to `OCC_FILTER == "0.5occA_0.5occB"` (or `OCC_FILTER_ALTLOC == "{'A': 0.5, 'B': 0.5}"`). Other splits never reach any figure.
- **Model name aliasing.** Raw CSV names go through `MODEL_NAME_ALIAS` (e.g. `boltz2_X-RAY_DIFFRACTION → boltz2X-RAY`, `rf3-fixed → rf3`). Always work with the canonical form downstream and look it up via `MODEL_DISPLAY[name]` for figure titles.
- **Disambiguating shared `model` columns.** When several CSVs store the same raw `model` value, two patterns coexist:
  - **At merge time** — `merge_new_rscc_lddt.py` rewrites `model` to a canonical name (e.g. `"rf3-chiral-off"`). Use this for new merge scripts.
  - **At load time** — `load_min_altloc_rmsd(csv_path, model_name_override="protenix-mini")` because both protenix-mini and protenix-tiny CSVs store `model == "protenix"`. Use this for ad-hoc per-source loads where the merge step is skipped.
- **Model groups.** `MODELS_ALT = {"protenix-mini", "protenix-tiny", "rf3-chiral-off"}` are the alternative variants. `filter_by_model_group(df, "main")` drops them; aggregated figures are main-models-only by convention.
- **Guidance buckets.** `GW_TO_STRENGTH` maps numeric gw → `{baseline, weak, medium, strong}`. `OPTIMAL_GW[model]` is the per-model "best" gw used by every "optimal vs baseline" figure (`0.1` for the three main models, `0.01` for RF3).
- **Pinned cross-figure highlights.** `6B8X 181-183` and `5MHX 162-164` are circled in blue (`#1565C0`) on every aggregated scatter (RSCC, LDDT, min-RMSD). Pin new always-circled selections in the same color so they read as the same family across figures.
- **`exclude_nan_rscc` toggle.** When on, the notebook blacklists proteins with any NaN RSCC across the four main models (via `get_rscc_nan_proteins`) and threads the blacklist through every loader output. Don't bypass it.
- **Pre-cutoff training cluster source.** `load_training_cluster_counts("lddt_assignment_counts.csv")` adds `frac_a_cluster` / `x_one_minus_fa`. `load_deposition_cutoffs(...)` auto-fills `1VME, 4OLE, 5SOP, 6B8X` (pre-2018, absent from the CSV) as Pre.
- **Pandas + polars split.** Analysis code is pandas; only the one-shot merge scripts use polars.

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

For non-scatter families look at the existing implementations: `make_by_class_combined_grid` (4×4 grid with a two-source merge for the domain-shift column), `make_spiral_figure` / `make_aggregated_spiral_figure` (trajectories through guidance space), and the `make_training_bias_*` family (figures 29–36, fully documented in `data_in_each_plot.md`).

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
