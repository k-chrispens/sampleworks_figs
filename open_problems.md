# Open scientific and pipeline questions

These are review gates, not permission to silently change published calculations.

## Figure 2C: what counts as an RSCC observation?

- **Evidence:** notebook 02 averages region RSCC per protein before drawing the
  baseline distribution; the publication caption describes altloc selections.
- **Current behavior:** `01_baseline.py` reconstructs protein means by default and
  labels the per-region alternative explicitly. Both retain the existing inputs
  and score formula. The tracked saved baseline SVG uses ECDFs.
- **Resolution needed:** approve the intended paper population after comparing
  the reconstruction with the published panel. Do not silently change weighting.

## Clustering-score state correspondence

- Some 50/50 A/B map rows have three/four predicted occupancy entries.
- The current loader takes the minimum over every entry and divides by map
  occupancy 0.5. This behavior is preserved in the first reconstruction.
- Establish the reference-state labels and intended treatment of extra entries,
  and reconcile the publication's formula with the upstream silhouette definition.

## Later panels

- **Figure 3C:** the existing code uses a difference of marginal medians, while
  “median improvement” could be read as a median of paired model differences.
  The reconstruction preserves the former and exposes model pairs. Its minimum
  three-model rule counts rows, not finite pairs; confirm the intended paper rule
  separately rather than silently change the calculation.
- **Figure 3B:** reconstruction preserves marginal medians/IQRs. Mini/tiny have
  35/34 finite RSCC contributors versus 40 geometry contributors; decide separately
  whether future comparisons should use a common complete-case population.
  Guidance-weight outlines retain the existing constants, not a new optimization.
- **Figure 3E:** publication says mean across predictors; current aggregated
  minimum-RMSD notebook uses median. `05_min_rmsd.py` now defaults to that median
  reconstruction and offers a labeled across-model mean alternative. The
  within-ensemble minimum is unchanged. Establish the historical computation
  before declaring either view the pinned paper recipe.
- **Classification:** publication describes all-atom lDDT; the local source is
  documented as backbone-lDDT. Confirm definition/thresholds before pinning it.
- **Sampling budgets:** existing 5MHX/6B8X CSVs provide budget axes. Document their
  generation history and whether Supplement S6 also requires clustering data.
- **Historical reference:** source checksums identify the current reconstruction;
  equivalence to the published archive/code revision is not yet verified.

## Implementation decisions that remain deferred

- Ask the user to approve the merged-artifact storage format before implementing
  persistence. The first figure does not make that choice.
- General import/conflict resolution follows the first working figure.
- Mol* requires supplied ensemble/reference/map examples; integrated sync is later.
