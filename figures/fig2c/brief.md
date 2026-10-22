# Figure 2C — first reconstruction

claim: Baseline RSCC and conformational-recovery distributions differ among predictors.
venue: Paper panel reconstruction and interactive marimo notebook.
concept: Two stacked empirical CDFs, sharing model identity and colors.
data: Ten CSVs copied exactly from the user's existing local results snapshot;
the snapshot hash and individual source/code checksums are in the versioned
`v01-provenance.json` / `v02-provenance.json` files.

## Controls and inspection

- Results folder + Load snapshot: explicitly choose which read-only inputs are used.
- Plot legend: click to toggle a model in both panels; double-click a visible model to isolate it.
  All variants remain available, initially hidden. There is no separate model picker.
- RSCC observation: reconstruct protein means or explore per-region weighting.
- Hover/box/lasso and observation-table row selection: inspect values and contributors.
- The plot toolbar exports the visible SVG. Separate downloads provide a fixed
  main-model reconstruction SVG, all-model CSVs, and a recipe with explicit scopes.

## States

- **default:** four main predictors, protein-mean RSCC, exact ECDFs. There are 40
  RSCC means per model; finite clustering counts are 791/785/791/782 for
  Boltz2-MD/Boltz2-XRD/Protenix/RF3.
- **region:** same models with per-region RSCC. This is an explicit alternative,
  not an adjudicated correction to the publication.
- **cross-family:** arbitrary subsets work, e.g. Boltz2-MD plus RF3-Chiral-Off.
- **all hidden:** legend clicks can hide all traces without changing the loaded data.

## Evidence and limits

- Both panels/all main models match the existing helper's observed scores and
  cumulative proportions exactly. A finite display anchor replaces its initial
  negative-infinity sentinel; it is not a plotted observation.
- Unit tests cover missingness, weighting, ties, source tracing/immutability,
  renderer parity, and native-selection payloads after table sorting.
- Browser checks cover grouped legend toggles/isolation, current-view SVG export,
  both populations, and chart/table selection. Marimo 0.23.2 remains sufficient.
- See `../../open_problems.md`: exact historical paper identity, RSCC population,
  and multi-entry occupancy semantics are not resolved by this reconstruction.

log:
- v01 — reconstruct current results and add exploration without changing metric
  definitions; move the static legend outside the plot.
- v02 — remove the model-selection form; use native legend interaction and
  distinguish visible-view SVG export from fixed reconstruction/all-model exports.
