# Figure 3C — aggregated RSCC improvement

claim: Density guidance changes region-level fit unevenly; baseline fit limits the available improvement.
venue: Paper panel reconstruction and interactive marimo notebook.
concept: Region residual scatter above paired baseline/guided marginal densities.
data: The two main-model max-RMSD RSCC files from the existing local snapshot.
Checksums and settings are in `v01-provenance.json`.

## Interpretation and interaction

- The model population is fixed to the four main models at their existing chosen
  guidance weights; require three paired model rows as in notebook 02.
- Plot separate marginal medians and subtract them. Do not substitute a median
  of per-model deltas. Non-finite contributor counts remain visible.
- Use the original top-improved/decreased/far-from-ceiling/curated ring rules.
- Hover/select regions and inspect model pairs plus original source record numbers.
- Native legend controls display only, not model membership. Density samples are
  not individual regions.
- Toolbar SVG exports the visible view; reconstruction/data downloads retain the
  full fixed analysis.

## Checks and limits

- All 791 region coordinates and residuals match notebook 02 exactly.
- Marginal KDE coordinates match the existing renderer to 1e-12.
- All four models contribute 791 paired rows; one Boltz2-MD pair is non-finite.
  The legacy row-count rule and NaN-skipping medians are preserved and reported.
- Historical paper identity/aggregation interpretation remain open; see
  `../../open_problems.md`.

log: v01 — reconstruct Figure 3C with shared final coordinates/densities and
region-to-model-pair inspection, without changing the scientific calculation.

log: v02 — remove fixed segment labels while retaining all ring categories;
bundle the reconstruction SVG with a highlighted-region table and recipe.
