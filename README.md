# Sampleworks figures

## Focused reconstructions

These five figures have **initial interactive drafts**, with hover, linked
plot/table selection, and SVG/data downloads. They reconstruct the existing
results; they are not yet verified historical publication recipes.

| Notebook | Paper panel |
|---|---|
| `01_baseline.py` | Figure 2C: baseline RSCC and clustering-score distributions |
| `02_guidance.py` | Figure 3B: guidance-strength versus geometry tradeoff |
| `03_rscc_improvement.py` | Figure 3C: aggregated RSCC improvement and model-pair inspection |
| `04_clustering_improvement.py` | Figure 3D: clustering-score improvement and occupancy/silhouette inspection |
| `05_min_rmsd.py` | Figure 3E: minimum RMSD to A/B, with explicit across-model mean/median |

### Supplementary interactive figures

| Notebook | Remaining scope |
|---|---|
| `06_fit_supplements.py` | Per-model residual/RMSD panels, cross-metric comparisons, class grid, and summary tables |
| `07_training_tradeoffs.py` | Individual-protein trajectories and training-bias figures |
| `08_sampling_budgets.py` | Recycles, diffusion steps, and ensemble-size sweeps |

These now use real Plotly figures with hover, zoom, legend controls, and linked
inspection tables. They load `results/` immediately, without the focused
notebooks' source form. Numerical summary tables remain tables; legends and the
data-free calibration explainer do not invent selectable observations.
Mol* ensemble viewing and the general raw-output importer are still pending.

The old `02_fit_quality.py`, `03_tradeoffs.py`, and
`04_recycles_steps_ensembles.py` entry points are retired, not compatibility
wrappers. Their last full snapshot is commit `4667be8`. Duplicate baseline and
aggregate panels/renderers were deleted, while unique analyses moved to 06–08.
Historical figure assets remain intact.

### Linked plot and table selection

All eight notebooks link observation/region tables to the plots:

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

In 06–08, expand a figure's named **inspection table near the bottom**. Each plot
has its own clear button; a lasso in one plot does not filter another plot.
Budget rows represent mean/SD points; their second table lists the contributing
selection measurements. Densities, regression lines, uncertainty bands, and
class-grid background clouds are not additional observations.

The toolbar camera downloads the visible view as SVG. For a fixed reconstruction,
choose a figure near the bottom and click **Build SVG + data bundle**, then
download its ZIP. The ZIP includes the uninspected SVG, prepared points/summaries,
a circled-region table when applicable, budget contributors when applicable, and
a recipe containing the figure specification. These supplementary recipes record
prepared values, **not a verified historical input snapshot**; upstream
provenance and publication-correspondence gates remain open.

### Fixed Plotly SVG exports: Chrome prerequisite

Notebooks 06–08 use Kaleido for Python-side SVG export (fixed bundles and optional
Auto-save). Viewing figures and browser-toolbar SVG export do not need Chrome
on the Python host. For fixed exports, install Chrome/Chromium **where Python
runs**, not merely on the laptop viewing the notebook:

```sh
uv sync --locked
uv run plotly_get_chrome -y
```

Alternatively, point Kaleido at an existing installation:

```sh
export BROWSER_PATH=/absolute/path/to/chrome-or-chromium
```

Set this before starting marimo, in the local shell for local execution or the
pod shell for remote execution. A missing Chrome installation blocks fixed
exports, not interactive viewing. Fixed exports are explicit because launching
Chrome for every panel on notebook startup is wasteful on a small pod.

## Demo setup

The recommended demo is **one read-only marimo server for this directory**.
Its chooser lists all eight notebooks; no separate server per figure is needed.
Use the current checkout containing notebooks 01–08, and run commands from its
root. The demo uses port 2730 below. If it is already running, reuse its live
token-bearing URL and forward rather than starting a second server on that port.

### Inputs and prerequisites

Provide the existing flat evaluation-results snapshot in `results/`. These inputs
are **not** installed by `uv sync` or supplied by a Git clone. Required filenames
and populations are listed in [`data_in_each_plot.md`](data_in_each_plot.md).
To demo all eight notebooks, provide inputs for the retained analyses too:
base RSCC, variant RMSD, classification, training metadata, and budget sweeps.

`results/` is ignored by Git but is **not excluded by `.actlignore`**, so a snapshot
placed there is included in normal `actl` sync. Use the flat results snapshot,
not a raw `occ_sweep_results*` folder: the raw-export importer is not built yet.

The folder must be accessible **where marimo runs**. Notebooks 01–05 can select
another folder, or default to `SAMPLEWORKS_RESULTS`; 06–08 currently require
`results/` relative to the project root. A laptop path is not automatically a pod
path. Inputs are read-only; keep **Auto-save figures** off in 06–08 for a
demo without figure-file exports (their optional saves go to `figures/`, never
the input folder).

For pod hosting, install/configure `actl` and authenticate to the `diffuse`
namespace. A CPU workspace is sufficient; no GPU or laptop PVC mount is needed.
For local hosting, install `uv`; the lockfile supplies the Python dependencies.

### Run on the pod, view in your local browser

The server binds to the pod's loopback address. A local forward makes it
reachable in your browser, with a fresh marimo access token.

**1. In a local terminal, choose one way to open a forwarded pod shell.**

If `sampleworks-figs` and its code-sync session are already running:

```sh
actl pod exec sampleworks-figs --namespace diffuse --port 2730
```

This opens another remote shell and forwards local port 2730 to pod port 2730.
It can run alongside the existing sync session. It does **not** start sync itself;
make sure the active sync is using this checkout.

Open the forwarded shell first, then start marimo **inside** it. Command-mode
`actl pod exec ... -- <command>` does not establish the same port-forwarding
session; do not use it as a replacement for this interactive-shell step.

If the workspace or sync is not running, use this **instead**, from this repo:

```sh
actl pod up sampleworks-figs --namespace diffuse \
  --profile cpu-small --image python --pvc-size 10Gi --port 2730
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
uv run marimo run . \
  --host 127.0.0.1 --port 2730 --headless --token --session-ttl 30
```

`uv sync --locked` restores the environment after a pod restart; `/tmp` is
disposable. The commands do not run structure prediction or change input CSVs.
The `.` is intentional: it opens the notebook chooser. The 30-second session
TTL releases a notebook session shortly after its tab disconnects; it does not
terminate the server or an actively connected notebook.

**3. Open the token-bearing URL printed by marimo in your local browser.**

It has the form `http://127.0.0.1:2730/?access_token=...`. Open that URL to
authenticate and see the chooser. Click a notebook; it opens in a separate tab.
In 01–05, click **Load snapshot** with `results` selected. Notebooks 06–08 load
their inputs immediately.

`marimo run` prevents editing Python source in the browser; figure controls still
work. Keep the terminal open while using the notebook. Close finished notebook
tabs before loading another on the 4 GiB `cpu-small` pod; the chooser is not a
memory-free way to keep eight loaded sessions active. Use a larger-memory profile
for concurrent loaded notebooks.

If port 2730 is occupied by another service, use
another free port in both `actl --port` and `marimo --port` (or follow any remapped
local port printed by `actl`). Do not omit `--token`, bind to `0.0.0.0`, or reuse
an old session's token URL. A `localhost` link works only on the machine with the
forward; use screen sharing for a remote demo rather than treating it as a public URL.

**Stopping and editing:** Ctrl-C stops marimo and returns to the pod shell;
`exit` closes that shell and its forward. Edit source files locally, let sync
propagate them, then restart marimo to pick up helper-module changes. If the page
becomes unreachable after closing the terminal, reopen the forward and check
whether marimo is still running before starting another server. When finished
with the entire workspace, `actl pod down sampleworks-figs --namespace diffuse`
releases its CPU resources and keeps its home PVC; it also stops other sessions
using that workspace. No broad `pkill` or process-group cleanup is needed.

### Run entirely locally

No `actl` session or port forward is needed. From this checkout:

```sh
uv sync --locked
uv run marimo run . --host 127.0.0.1 --port 2731 --token --session-ttl 30
```

Marimo opens the local chooser and prints its authenticated URL. Port 2731 lets
this run alongside the remote example. Choose a notebook and use the same
snapshot-loading steps as above. `SAMPLEWORKS_RESULTS` or the source form can point
01–05 at another local folder; 06–08 still need the project-root `results/` inputs.

If you want to edit notebook code locally as well as use its controls, substitute
`marimo edit` for `marimo run`, retaining the loopback address and token. This is
for development, not the read-only demo:

```sh
uv run marimo edit 01_baseline.py --host 127.0.0.1 --port 2731 --token
```

Stop a local server with Ctrl-C. Both launch methods use the same data preparation
and downloads; neither writes into the selected results folder.

For a single-notebook demo, replace `.` with its filename in either launch command.
You do not need another server/port just to switch figures in the chooser.

### Suggested demo flow

1. Open **03 RSCC Improvement**, then **Load snapshot**. The current snapshot
   produces 791 plotted regions; counts are data-dependent, not a schema guarantee.
2. Hover for region identity and values. Box/lasso-select points to narrow the table.
3. Select a table row to add a black outline; change rows after zooming or toggling
   a legend entry to show that view state is preserved.
4. Use **Clear plot selection** to restore the full table; uncheck rows to remove
   inspection outlines.
5. Use the camera for the visible SVG, or the reconstruction ZIP for the fixed
   SVG plus circled-region table and recipe.

In **01/02**, legend entries select model visibility; in **03–05**, they control
display categories, **not** the fixed contributing-model population. Figure 3E
also exposes the existing median reconstruction and the explicitly unverified
caption-described mean alternative. These scientific choices are not final paper
adjudications; see `open_problems.md`.

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
uv run marimo check --strict [0-9][0-9]_*.py
uv run marimo export script 01_baseline.py > _smoke.py
uv run python _smoke.py
rm _smoke.py
```

Repeat the export-script check with every notebook, including 06–08 after changing
retained analyses. Run checks sequentially on small pods.

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

For 06–08, use the directory chooser's authenticated URL and a notebook filename:

```sh
uv run --with playwright python tests/browser_supplements.py \
  "$DEMO_URL" 06_fit_supplements.py
```

Repeat for 07 and 08. This checks every observation-bearing panel's row outlines,
native box selection, preserved view state, visible camera SVG, and fixed bundles.
Restart `marimo run` after changing shared Python plotting modules: a running
server can retain imported modules even when source sync has delivered new files.

Set `FIGURE_URL` to its URL including the access token. Run this where that URL is
reachable (on the pod, use its loopback address and port).

The exported script loads `SAMPLEWORKS_RESULTS` (or `results`) with the initial
controls; browser sessions still require an explicit snapshot-load submission.

See `data_in_each_plot.md` for figure inputs and `open_problems.md` for unresolved
scientific correspondence to the paper. The merged-artifact format is still
deferred and needs explicit approval before implementation.
