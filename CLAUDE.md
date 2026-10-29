# Project guide

Analysis of AI structure predictions against crystallographic density and
alternate conformations. Python/pandas modules prepare data; marimo notebooks
present figures. `README.md` explains local/pod launch. `data_in_each_plot.md`
is the authoritative source/population/figure contract; consult it before
changing any plot. `open_problems.md` contains scientific approval gates.

## Entry points

| Notebook | Role |
|---|---|
| `01_baseline.py` | Figure 2C baseline ECDFs; `baseline.py` |
| `02_guidance.py` | Figure 3B aggregate guidance/geometry; `guidance.py` |
| `03_rscc_improvement.py` | Figure 3C aggregate RSCC improvement |
| `04_clustering_improvement.py` | Figure 3D aggregate clustering improvement |
| `05_min_rmsd.py` | Figure 3E aggregate minima; `min_rmsd.py` |
| `06_fit_supplements.py` | Per-model, cross-metric, class-grid and summary analyses |
| `07_training_tradeoffs.py` | Individual-protein spirals and training-bias suite |
| `08_sampling_budgets.py` | Supplementary 5MHX/6B8X budget sweeps |

Notebooks 03/04 share `region_improvement.py`. Legacy notebook names and duplicate
aggregate renderers were retired; their last full snapshot is commit `4667be8`.
Do not recreate alternate implementations for the same panel.

Notebooks are `marimo.App`s, not Python scripts:

```sh
uv sync --locked
uv run marimo edit 01_baseline.py
uv run marimo run 01_baseline.py --host 127.0.0.1 --token
```

Use the same commands with the other notebook filenames. Source stays local;
when using actl, run code/tests on the existing pod through normal source sync.

## Modules and boundaries

- `style.py`: rcParams, model palettes/display/aliases, guidance constants, save helpers.
- `loaders.py`: CSV parsing, source declarations, historical 50/50 filters,
  canonical names, occupancy metrics, classification and training metadata.
- `transforms.py`: DataFrame transformations and labels.
- `baseline.py`, `guidance.py`, `region_improvement.py`, `min_rmsd.py`: focused
  preparation, audit/provenance, static/interactive renderers. Both renderers
  consume the same final coordinates; do not independently aggregate in a renderer.
- `plot_interaction.py`: persistent Plotly widget, stable observation lookup,
  display-only table selection outlines.
- `plot_scatter.py`, `plot_classified.py`, `plot_tradeoff.py`,
  `plot_training_bias.py`, `plot_budgets.py`: supplementary/per-model Plotly figures;
  prepared values and summaries live in figure metadata, shared by viewing/export.
- `supplementary_plot.py`: observation traces, persistent linked panels, inspection
  tables and explicit fixed SVG/data bundles for 06–08. Callback-bearing chart,
  button and table objects must remain strongly referenced. Selection generations
  prevent stale filters from leaking into a rebuilt chart.

Keep dependencies acyclic. Reuse current helpers and source declarations rather
than introducing parallel loaders or a general notebook framework.

## Scientific rules

- Runtime inputs are under **`results/`**, ignored by Git. Focused notebooks 01–05
  accept a read-only explicit folder (`results` or `SAMPLEWORKS_RESULTS` by default),
  fingerprint sources, and retain original data-record ordinals. 06–08 currently
  use hardcoded `results/...` paths. Never write to an input folder.
- All current loaders filter to 50/50 A/B map occupancy. Canonicalize model names
  through `MODEL_NAME_ALIAS`, display via `MODEL_DISPLAY`. The three alternatives
  are `MODELS_ALT`; main-model aggregation must not accidentally include them.
- Use max-RMSD-region RSCC everywhere **except the class grid**, which needs all
  base selections. Main and alternate source declarations are in
  `RSCC_MAX_RMSD_SOURCES`. Chiral-off lacks a max-RMSD source and uses its base CSV.
- Mini/tiny share raw `model == "protenix"`; chiral-off shares raw `model == "rf3"`.
  Set `model_name_override` at load time. For path-tagged geometry, apply variant
  overrides **after tagging** or variants silently collapse into main models.
- Guided comparisons use `OPTIMAL_GW`, not a newly selected optimum.
  Guidance-strength buckets come from `GW_TO_STRENGTH`, not weight ordering.
- Clustering uses the minimum of **all** predicted occupancy entries times
  silhouette/reference-map occupancy. Training-bias B occupancy uses the **second**
  entry instead. Do not conflate those reductions.
- 3B preserves marginal protein medians/IQRs and differing finite cohorts.
  3C/D subtract separately computed medians and count paired rows, not finite pairs.
  3E aggregates complete A/B model pairs; median is the reconstruction default,
  mean is an explicitly unverified caption-described alternative.
- Respect the NaN-RSCC blacklist when enabled in 06/07: it covers all four main
  models and must be threaded through every loaded metric frame in that notebook.
- Classification uses the narrowest containing protein/chain residue interval;
  unmatched/unparseable regions are reported and dropped. The current classifier
  is backbone-lDDT; publication all-atom correspondence is unresolved.
- Per-model RMSD (including mini/tiny and 6NI5/6NI6 rings) is distinct from the
  fixed main-model aggregate. Do not delete it as a duplicate of Figure 3E.
- Model-group controls affect class/training-bias grids, with distinct variant
  output suffixes. Per-model controls are separate; never overwrite main outputs
  with a variant-family export.
- Circled-region static reconstructions have no fixed segment text. ZIPs include
  the SVG, companion highlighted-region table and recipe. Toolbar SVGs follow the
  visible view; fixed reconstruction/CSV scopes are explicit and independent.
- Keep authored molecular panel assets and the user-retained
  `results/boltz2_msa_depth.tsv~`. Merged-artifact storage still needs user approval.

## Marimo and interaction rules

- Use leading-underscore locals within cells; expose only real inter-cell values.
- Cell arguments auto-rewrite from static references. `pd.query("@VAR")` does not
  count as a static reference: rely on marimo's injected globals rather than
  manually adding an argument that marimo will remove.
- Keep callback-bearing UI elements in returned/global variables; anonymous
  elements can lose their Python callback when collected.
- Main tables use a persistent `FigureWidget` through `mo.ui.anywidget`.
  Update markers in place on table selection; recreating the chart loses zoom
  and legend state. Stable IDs, never sorted table row positions, identify points.
- Empty anchors and density vertices are not observations. Black inspection
  outlines never enter fixed reconstruction exports. A new analysis resets inspection.

## Validation

Run sequentially on the pod when resources are small:

```sh
uv run python -m unittest discover -s tests -v
uv run marimo check --strict [0-9][0-9]_*.py
uv run python -c "import style, loaders, transforms, plot_scatter, plot_classified, plot_tradeoff, plot_training_bias, plot_budgets"
```

Execute each notebook through marimo's exported-script path, not `python notebook.py`:

```sh
uv run marimo export script 06_fit_supplements.py > _smoke.py
uv run python _smoke.py
rm _smoke.py
```

Keep the temporary script in the project root for relative imports. Run this for
all changed notebooks, including retained supplements. Browser behavior is separate:
`uv run --with playwright python tests/browser_linked_selection.py "$FIGURE_URL"`.
Test a single-panel figure and 3E's linked panels after interaction changes.
