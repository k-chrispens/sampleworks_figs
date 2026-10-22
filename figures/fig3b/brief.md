# Figure 3B — guidance versus geometry

claim: Guidance improves median density fit but increases geometric penalties.
venue: Paper panel reconstruction and interactive marimo notebook.
concept: Model-colored trajectories through median delta RSCC / geometry space,
with marginal IQRs, increasing marker sizes, and selected-weight outlines.
data: Twenty files from the existing local results snapshot. SHA-256 fingerprints
and analysis settings are recorded in `v01-provenance.json`.

## Interaction

- Load snapshot: choose a read-only results folder.
- Click model names in the native legend to toggle them; double-click a visible model to isolate.
  Main models start visible and variants remain available but hidden.
- Hover reports guidance, raw penalty, and protein/finite contributor counts.
- Select points/table rows to inspect the protein-level metric means and deltas.
- The camera button exports the visible SVG. Separate exports preserve the fixed
  main-model reconstruction and all-model data/provenance.

## Scientific reconstruction

- Match notebook 03's member/region means, baseline pairing, marginal medians/IQRs,
  and NaN-skipping geometry penalty. No new optimal guidance weights are chosen.
- Use the same prepared coordinates for static and interactive renderers, including
  the symmetric-log display transformation. Raw values remain in tables and hover.
- Main-model contributions and transformed trajectories match the old implementation
  numerically (tolerance 1e-12).
- The main models have 40 finite RSCC/geometry protein contributors at each point.
  Mini/tiny have 35/34 finite RSCC contributors versus 40 geometry contributors;
  this legacy mismatch is exposed, not silently changed.

log: v01 — reconstruct the next generated paper panel with native legend exploration,
point-to-protein inspection, coverage diagnostics, and explicit export scopes.
