# Sampleworks figures

## Focused reconstructions

Focused notebooks reconstruct the existing results before adding exploration:

| Notebook | Paper panel |
|---|---|
| `01_baseline.py` | Figure 2C: baseline RSCC and clustering-score distributions |
| `02_guidance.py` | Figure 3B: guidance-strength versus geometry tradeoff |
| `03_rscc_improvement.py` | Figure 3C: aggregated RSCC improvement and model-pair inspection |
| `04_clustering_improvement.py` | Figure 3D: clustering-score improvement and occupancy/silhouette inspection |
| `05_min_rmsd.py` | Figure 3E: minimum RMSD to A/B, with explicit across-model mean/median |

The larger analysis notebooks remain available.

### Linked plot and table selection

All five focused notebooks link the main observation/region table to the plot:

- Box/lasso selection in the plot narrows the table; use **Clear plot selection**
  below the plot to restore all observations.
- Selecting table rows adds **black inspection outlines** to their points.
  Clearing the rows removes those outlines. Sorting/paging uses stable observation
  IDs, not row positions. For Figure 3E, the region is outlined in both panels.
- Table clicks preserve chart zoom and legend visibility. A legend-hidden model
  stays hidden; enable it in the legend to see its selected observations.
- Changing the snapshot or analysis controls starts a fresh inspection.
- The camera export includes the visible inspection outlines. Fixed reconstruction
  SVG/ZIP/CSV downloads do not change with these transient selections.

The plots use Plotly `FigureWidget` through marimo's `anywidget` support so marker
updates happen in place. Run `uv sync --locked` after updating this checkout.
The `cpu-small` pod example is intended for one active notebook at a time. Use a
larger-memory profile if you want several loaded figure notebooks concurrently.

### Inputs for either launch method

Run from the checkout/worktree containing `01_baseline.py`. Put the existing
results snapshot in `results/`, or choose another folder in the notebook.
The folder must be accessible **where marimo runs**: a laptop path is not
automatically available inside a pod.

`results/` is ignored by Git but is **not excluded by `.actlignore`**, so a snapshot
placed there is included in normal `actl` sync. The initial reconstruction used
exact copies from the primary checkout's
`/Users/kyroyo/Downloads/sampleworks_paper/for_team/results`. That machine-specific
path is not a default in the code. Use the flat results snapshot, not an individual
raw `occ_sweep_results*` folder: the general raw-export importer is not built yet.

### Run on the pod, view in your local browser

This follows the pattern in
`~/Projects/latent-testing/scripts/notebook_on_pod.sh`: the server binds to the
pod's loopback address, authentication uses a fresh token, and a port forward
makes it available in the laptop's browser. No GPU or local `actl mount` is needed.

**1. In a local terminal, choose one way to open a forwarded pod shell.**

If `sampleworks-figs` and its code-sync session are already running:

```sh
actl pod exec sampleworks-figs --namespace diffuse --port 2719
```

This opens another remote shell and forwards local port 2719 to pod port 2719.
It can run alongside the existing sync session. It does **not** start sync itself;
make sure the active sync is using this checkout.

If the workspace or sync is not running, use this **instead**, from this repo:

```sh
actl pod up sampleworks-figs --namespace diffuse \
  --profile cpu-small --image python --pvc-size 10Gi --port 2719
```

This creates/resumes the CPU workspace, syncs the current checkout, forwards the
port, and opens a remote shell. Do not start a second `pod up`/`pod sync` while
another sync session for this workspace is active.

**2. In that remote shell, run the following pod-side commands.**

```sh
cd /home/dev/workspace
export PATH="$PWD/.tools:$HOME/.local/bin:$PATH"

# One-time bootstrap if the minimal image does not already have uv.
if ! command -v uv >/dev/null 2>&1; then
  curl -LsSf https://astral.sh/uv/install.sh |
    UV_UNMANAGED_INSTALL="$PWD/.tools" sh
fi

# Keep dependency installation off the slower shared filesystem.
export UV_CACHE_DIR=/tmp/figure-uv-cache
export UV_PROJECT_ENVIRONMENT=/tmp/figure-venv
uv sync --locked
uv run marimo run 01_baseline.py \
  --host 127.0.0.1 --port 2719 --headless --token
```

`uv sync --locked` restores the environment after a pod restart; `/tmp` is
disposable. The commands do not run structure prediction or change input CSVs.
If your snapshot is elsewhere on the pod, set `SAMPLEWORKS_RESULTS` to that
**pod-side path** before launching marimo.

**3. Open the token-bearing URL printed by marimo in your local browser.**

It has the form `http://127.0.0.1:2719/?access_token=...`. Click **Load snapshot**,
then click model names in the plot legend to show/hide them, or double-click a visible model to
isolate one. `marimo run` prevents editing the
notebook's Python source in the browser; it does not disable interactive controls.
Keep the terminal open while using the notebook.

Port 2719 avoids the earlier development preview on 2718. If it is occupied, use
another free port in both `actl --port` and `marimo --port` (or follow any remapped
local port printed by `actl`). Do not omit `--token`, bind to `0.0.0.0`, or reuse
an old session's token URL.

**Stopping and editing:** Ctrl-C stops marimo and returns to the pod shell;
`exit` closes that shell and its forward. Edit source files locally, let sync
propagate them, then restart marimo to pick up helper-module changes. When finished
with the entire workspace, `actl pod down sampleworks-figs --namespace diffuse`
releases its CPU resources and keeps its home PVC; it also stops other sessions
using that workspace. No broad `pkill` or process-group cleanup is needed.

### Run entirely locally

No `actl` session or port forward is needed. From this checkout:

```sh
uv sync --locked
uv run marimo run 01_baseline.py --host 127.0.0.1 --port 2720 --token
```

Marimo opens the local browser and prints its authenticated URL. Port 2720 lets
this run alongside the remote example. Click **Load snapshot** to load `results`.
For inputs elsewhere on your laptop:

```sh
SAMPLEWORKS_RESULTS=/absolute/path/to/results \
  uv run marimo run 01_baseline.py --host 127.0.0.1 --port 2720 --token
```

If you want to edit notebook code locally as well as use its controls, substitute
`marimo edit` for `marimo run`, retaining the loopback address and token:

```sh
uv run marimo edit 01_baseline.py --host 127.0.0.1 --port 2720 --token
```

Stop a local server with Ctrl-C. Both launch methods use the same data preparation
and downloads; neither writes into the selected results folder.

For the other plots, substitute any focused notebook from the table above in the same launch commands. If
running both notebooks simultaneously, use a separate free port for the second
server and, remotely, a matching `actl --port` forward.

### What the first figure does

- Starts with the four main models; variants are available as hidden legend entries.
- Click a legend entry to toggle that model in both panels; double-click a visible
  model to isolate it.
- Uses exact empirical CDFs, as in the saved `04_09_baseline_distributions.svg`.
- Defaults to **protein-mean RSCC**, reproducing notebook 02's calculation. A
  separate **per-region RSCC** option explores the population described by the
  paper caption. Neither choice is silently presented as a verified paper recipe.
- Retains the existing per-region clustering-score definition and 50/50 occupancy
  filter. Three/four-entry occupancy arrays remain a documented scientific question.
- Shows source-value exclusions and observation counts.
- Supports hover and chart selection. Select observation-table rows to see their
  contributing CSV records; the table is also a keyboard-accessible alternative.
- The plot toolbar's camera button downloads the **currently visible view as SVG**.
- Separate downloads provide the **fixed main-model reconstruction SVG**, all-model
  observations/curve vertices, and a recipe with input/code checksums. Legend choices
  are browser view state, not a filter on those all-model CSV downloads.

There is no separate model-selection block. Legend interaction does not reload
the snapshot or recompute the statistics. Without a point selection, the inspection
table includes all loaded models, even hidden ones. Source reload is explicit;
inputs are never overwritten and downloads are built in memory.

The snapshot must contain the five RSCC sources listed in
`loaders.RSCC_MAX_RMSD_SOURCES` and the corresponding five `lddt_results_*.csv`
files. Missing files and ambiguous duplicate baseline records stop loading rather
than silently omit a model. This is a reconstruction adapter for today's exports,
not the future general input contract or merged-artifact storage implementation.

### Figure 3B: guidance and geometry

`02_guidance.py` reads the five existing RSCC sources plus each group's clashscore
and bond-angle/bond-length summary CSVs. It reproduces notebook 03's sequence:
member/region means per protein, guided-minus-baseline changes, then marginal
medians and IQRs across proteins. The symmetric-log vertical scale preserves zero
and negative penalties; hover and tables report untransformed values.

Use the same native legend controls. Marker size increases with guidance within
each model; black outlines mark the existing chosen guidance weights. Select
points and table rows to inspect the contributing proteins. Coverage tables show
missing baseline pairs, finite contributor counts, and undrawn groups. Different
RSCC/geometry populations are disclosed, not silently paired or filled.

The plot toolbar exports the visible SVG. Separate downloads contain the fixed
main-model SVG, all-model summary points, protein contributions, and provenance.
Initial reconstruction files are under `figures/fig3b/`. Like Figure 2C, this is
a current-code reconstruction, not an independently verified historical paper recipe.

### Figure 3C: RSCC improvement

`03_rscc_improvement.py` reads only the two main-model max-RMSD RSCC files and
pairs each model's baseline with its existing chosen guidance weight. Each region
requires at least three paired model rows. The plotted improvement is **median
guided RSCC minus median baseline RSCC**, not the median of individual changes.

Hover and select regions to inspect the contributing models, scores, guidance
weights, and original CSV records. The legend controls region/highlight/density
visibility; it does **not** change the fixed four-model aggregation. The original
green/red/orange/blue highlight rules and marginal Gaussian KDE settings are
retained, without fixed segment text labels in either renderer. Region identities
remain available through hover and selection. Finite-pair counts, unpaired regions,
and degenerate density estimates are reported explicitly.

The camera button exports the visible view. Separate downloads provide the full
reconstruction, region/model-pair tables, density coordinates, and provenance.
The reconstruction ZIP includes the unlabeled SVG, a companion circled-region
table with categories/colors and values, and the recipe. That table is also
displayed in the notebook and downloadable separately. Initial reconstruction
files are under `figures/fig3c/`.

### Figure 3D: clustering-score improvement

`04_clustering_improvement.py` reads the two main-model lDDT exports. It uses the
same pairing, marginal-median subtraction, minimum-three-model-row rule, and
density settings as the existing notebook. Green/red/blue circles follow the
same categories as 3C; orange uses three far-from-ceiling candidates, as in the
original 3D calculation. Hover/select to inspect the raw occupancy arrays,
silhouette components, scores, and source records.

Both 3C and 3D use `region_improvement.py` for their common preparation/rendering.
Neither prints fixed segment labels. Reconstruction ZIPs contain the SVG,
circled-region table, and recipe; the toolbar saves the visible view. The
multi-entry-occupancy interpretation remains an open scientific question.
Reconstruction files are under `figures/fig3d/`.

### Figure 3E: minimum RMSD to each altloc

`05_min_rmsd.py` reads the two main-model minimum-RMSD exports. Each input value
is already a minimum over an eight-member ensemble. The notebook combines those
values **across models**, requiring three complete A/B pairs per region at the
existing chosen guidance weights.

The default **median** reproduces the original notebook; the **mean** option
matches the caption's stated across-model reduction. This does not establish
which historical computation produced the published panel. Both views retain
the full-range and 0–3 Å panels, original curated circles (red in zoom only),
hover/selection inspection, coverage reports, and unlabeled static reconstructions
with companion tables. The selected aggregation is recorded in filenames and
recipes. Reconstruction files are under `figures/fig3e/`.

### Checks

```sh
uv run python -m unittest discover -s tests -v
uv run marimo check 01_baseline.py 02_guidance.py 03_rscc_improvement.py 04_clustering_improvement.py 05_min_rmsd.py
uv run marimo export script 01_baseline.py > _smoke.py
uv run python _smoke.py
rm _smoke.py
```

Repeat the export-script check with each other focused notebook.

The tests cover population weighting, missing values, exact ECDF ties,
interactive/static/reference parity, model aliases, source-record tracing, and
source-file immutability. Real-data reconstruction and browser inspection are
separate checks, not implied by the unit tests.

For a running focused notebook, the browser regression checks both directions,
clearing, preserved view state, visible SVG export, and unchanged
fixed reconstruction downloads:

```sh
uv run --with playwright playwright install chromium  # once in the test environment
uv run --with playwright python tests/browser_linked_selection.py "$FIGURE_URL"
```

Set `FIGURE_URL` to its URL including the access token. Run this where that URL is
reachable (on the pod, use its loopback address and port).

The exported script loads `SAMPLEWORKS_RESULTS` (or `results`) with the initial
controls; browser sessions still require an explicit snapshot-load submission.

See `data_in_each_plot.md` for figure inputs and `open_problems.md` for unresolved
scientific correspondence to the paper. The merged-artifact format is still
deferred and needs explicit approval before implementation.
