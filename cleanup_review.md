# Cleanup review and retirement decisions

Reviewed from `4667be8` on `refactor/prune-legacy-figures`, using three parallel
scout audits, an independent safety/test review, and direct caller/data checks.
The intended workload is a researcher exploring one dataset, not a multi-user
service. This pass removes duplicate execution paths rather than building a
general notebook framework.

At the retirement checkpoint (before the Plotly conversion below), measured
production Python: **8,129 → 7,124 lines**, a net removal of **1,005
lines (12.4%)**, including the new 06–08 entry points and excluding tests/cache.
The active figure specification shrank from 843 to 242 lines. Regression coverage
grew rather than being counted as a compression opportunity.

## Implemented

### One implementation per reconstructed main-paper panel

| Retired path | Surviving implementation |
|---|---|
| Old fit-quality baseline-distribution cell and `plot_distributions.py` | `01_baseline.py` / `baseline.py`, exact ECDFs |
| Old aggregate guidance overlay and per-model aggregate spirals | `02_guidance.py` / `guidance.py`; all seven models, native legend subsets/isolation |
| Old aggregate RSCC residual cell and aggregate scatter renderer | `03_rscc_improvement.py` / `region_improvement.py` |
| Old aggregate clustering residual cell and aggregate scatter renderer | `04_clustering_improvement.py` / `region_improvement.py` |
| Old aggregate minimum-RMSD cell | `05_min_rmsd.py` / `min_rmsd.py` |

Before retiring aggregate guidance rendering, all seven old per-model median
trajectories were compared with focused preparation on the current snapshot:
5 points each for Boltz2-MD/XRD, Protenix/mini/tiny, 6 each for RF3/chiral-off;
coordinates matched to 1e-12. Existing reconstruction work had already checked
the aggregate RSCC, clustering and median-RMSD coordinates against their old paths.

Old-only convenience modes—baseline KDE controls and batch generation of multiple
aggregate spiral files—are removed, not reproduced behind compatibility wrappers.
The focused notebooks own the explicit current reconstruction/export scopes.

### Preserve unique analyses, retire misleading entry-point names

| Old entry point | Remaining notebook |
|---|---|
| `02_fit_quality.py` | `06_fit_supplements.py`: per-model residual/RMSD, cross-metric, class-grid, summary tables |
| `03_tradeoffs.py` | `07_training_tradeoffs.py`: individual-protein spirals and training bias |
| `04_recycles_steps_ensembles.py` | `08_sampling_budgets.py`: both supplementary budget sweeps |

This is not just renaming: duplicate output cells, obsolete controls/imports,
the 346-line aggregate scatter function, the 213-line aggregate spiral function,
and the 107-line distribution module are deleted.

Important corrections to initial audit suggestions:

- **Per-model RMSD is not superseded by Figure 3E.** Keep mini/tiny, the original
  per-panel highlights, and their shared inputs. Summary tables also still need
  baseline/guided RMSD data.
- **Aggregate clustering is superseded by Figure 3D.** More detailed inspection
  UI does not make the old scientific panel a distinct analysis.
- Aggregate guidance curves are not the individual-protein trajectory family.
  Only the former moved to focused Figure 3B.

### Remove dead infrastructure without weakening scientific checks

- Replace the baseline test's dependency on the retired distribution renderer
  with explicit independent ECDF values, retaining ties and static/interactive parity.
- Remove unused `polars` (and its resolved runtime package). Keep `pyarrow` for
  marimo's table/data workflows, and SciPy/seaborn for live analyses. AdjustText
  was retained at this checkpoint and removed in the subsequent Plotly conversion.
  A missing direct import alone is not sufficient to delete a runtime/tool dependency.
- Remove the tracked `__marimo__/session/initial_40_analysis.py.json` cache.
  Its notebook no longer exists; the directory is already ignored.
- Replace the contradictory historical source inventory with one current figure
  contract. Previous data reorganizations remain available in Git, not as
  conflicting “canonical” instructions.

### Three concrete retained-analysis bugs fixed

1. Classification containment crashed on a parseable unmatched region because
   a NumPy boolean array was accessed through `.values`. It now reports/drops that
   row as documented. Tests cover exact/nested matches, narrowest containment,
   domain shift, unmatched and unparseable selections.
2. The default NaN-RSCC blacklist omitted RF3 despite the documented four-main-model
   rule. RF3 is now included. A regression fixture distinguishes RF3-only failure
   from a failure confined to chiral-off.
3. Training-bias inputs independently reloaded lDDT without applying the enabled
   NaN-RSCC blacklist. The lDDT, cluster and cutoff frames now follow the same
   exclusion as the geometry/RSCC branch. Disabled-toggle results are unchanged.

Additional fixtures protect training B occupancy versus minimum occupancy, empty/
one-sided/mixed training clusters, cutoff labels, unaggregated per-model RMSD,
and the three-versus-five orange-candidate distinction in Figures 3D/3C.

## Next compression targets, in order

1. **Completed in the follow-up Plotly conversion:** 06–08 now use real Plotly
   figures, persistent widgets, linked inspection and explicit fixed SVG/data
   bundles. Static label placement, coincident-label grouping/expansion and the
   final AdjustText dependency are removed. Per-model populations and colored-ring
   membership remain in the prepared data. Shared helper code handles presentation
   and observation IDs, not scientific population selection.
2. **Adjudicate training-bias redesigns against the paper.** IDs 29, 34 and 35
   include alternative presentations; delete rejected designs once their exact
   publication role is agreed. Do not treat every historical redesign as mandatory.
3. **Consolidate mechanical snapshot bookkeeping.** Four modules repeat
   hash-before/read/hash-after, source ordinals and manifest assembly. One small
   read helper can remove about 35–55 net lines. Keep family selection/validation
   and differing scientific missingness rules explicit at callers.
4. **Share reconstruction bundle assembly, not a notebook framework.** 3C/3D
   repeat SVG/highlight CSV/recipe/ZIP assembly; about 30–40 net lines can be removed.
   Keep source forms, scope labels and marimo dependencies visible in notebook cells.
5. **Retire machine-specific copy scripts before rebuilding molecular assets.**
   Existing scripts hardcode destinations outside the current checkout and use
   forced copies/truncated manifests. Their frozen MANIFEST/CIF/map inputs should
   survive; automatic copying should not become part of the read-only figure path.
6. **Make transfer failures explicit before using the old CSV pull script again.**
   `pull_remote_results.sh` suppresses SSH errors and can report success with no
   matching files when the source is unreachable. Separate listing failure from
   a successful listing with zero regex matches, or retire the script when the
   mounted-folder importer replaces it.

## Deliberately not compressed

- Different scientific reductions and unresolved publication choices in
  `open_problems.md`. Replacing them with one uniform aggregate would change results.
  In particular, Figure 32 still includes evaluated-but-unassigned A = B = 0
  rows, and budget fallback still operates after NaN removal. Their actual behavior
  is documented rather than silently “corrected” by this preservation pass.
- Small explicit marimo state/table callbacks: hiding them behind a framework
  saves few lines and obscures reactivity and export scope.
- Authored molecular graphics, CIFs, maps and paper layouts. They dominate the
  roughly 982 MiB tracked checkout, but are not dead code. Similarly named 5SOP
  layout/image pairs have different hashes; names do not establish duplication.
  Several existing `examples_panels/` ensembles/maps may support future Mol*
  testing after their provenance is verified.
- Git history, publication artifacts, or the deferred merged-artifact storage
  decision. No history rewrite, data deletion, storage choice, or upstream edit
  is part of this cleanup.

## Retirement-checkpoint verification gates

- **Completed on the existing pod:** 42 unit tests, locked dependency sync,
  strict marimo checks on all eight notebooks, lint on changed code, and export/
  execute of all eight notebooks. Named-cell tests exercise the training producer
  for all three model groups with exclusions off/on and the budget overlap rule.
- **Not repeated in this pass:** browser interaction/pixel checks and historical
  publication-archive equivalence. The focused interaction implementation is unchanged.
  Supplement execution still reports existing layout/label-placement warnings and
  missing training-cluster joins; these were not hidden or filled.

- Run the complete unit suite and strict marimo checks on **all eight** notebooks.
- Export and execute every notebook from the project root, sequentially on small pods.
  The supplement inputs must be present; focused-only smoke tests are insufficient.
- Keep numerical population/pairing/highlight tests, not SVG byte equality as
  a scientific oracle. Explicitly record the RF3 blacklist correction as a bug fix.
- Preserve the linked-widget browser tests and the static reconstruction scope.
  Browser sorting/paging, visible SVG contents, all CSV scopes and analysis-control
  resets need stronger automation; the current browser test does not establish them.

## Follow-up Plotly conversion verification

- The current suite has **53 tests**, including independent residual/KDE,
  quantile-bin/sample-SD, regression, subset, budget-contributor, legend-group,
  observation-identity, clear-event and stale-generation checks.
- All eight notebooks pass strict marimo validation and exported-script execution.
- `tests/browser_supplements.py` exercises the 18 supplementary plots, selecting
  rows in every observation-bearing panel; legends/explainer contain no observations.
  It checks native box selection, in-place outlines, preserved view, visible camera
  SVG, fixed bundles without inspection, and model/group control changes.
- Class-grid and calibration screenshots were inspected; full-width notebooks and
  shared grid headings replace the overlapping headings seen at the default width.
- Historical publication equivalence and comprehensive browser sorting/paging
  remain unverified. Recipes for 06–08 describe prepared data, not fingerprinted
  atomic source snapshots. Paper provenance and scientific approval gates remain open.
- Browser/Kaleido testing shares the small pod's memory with the server. Use one
  active notebook, and restart the server after shared-module changes; otherwise
  cached imports can make a browser check test stale code.
