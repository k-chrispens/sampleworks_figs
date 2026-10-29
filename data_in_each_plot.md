# Figure data and analysis contract

This is the current figure inventory, not a history of previous CSV layouts.
Historical sources and implementations remain in Git (`4667be8` is the last
pre-retirement snapshot). `open_problems.md` records unresolved scientific choices.
Figures below reconstruct the available results; historical publication identity
has not been established.

## Sources

All runtime CSVs live in `results/`, which is ignored by Git. Focused notebooks
01–05 accept another folder via the source form or `SAMPLEWORKS_RESULTS`.
Supplement notebooks 06–08 currently use explicit `results/...` paths.

The five source groups are `boltz_protenix`, `rf3_fixed`, `protenix_mini`,
`protenix_tiny`, and `rf3_chiral_off`. Apply the mini/tiny/chiral-off model override
when loading: their raw model/path identities overlap the main models.

| CSV family, relative to `results/` | Used for |
|---|---|
| `rscc_results_max_rmsd_{boltz_protenix,rf3_fixed,protenix_mini,protenix_tiny}.csv` plus `rscc_results_rf3_chiral_off.csv` | Primary max-RMSD-region RSCC; all non-classification analyses |
| `rscc_results_{boltz_protenix,rf3_fixed,protenix_mini,protenix_tiny,rf3_chiral_off}.csv` | Base per-region RSCC for the conformational-class grid only |
| `lddt_results_<group>.csv` | Silhouette, predicted cluster occupancies, derived clustering score |
| `clashscore_metrics_<group>.csv` | Ensemble-member clashscores for geometry tradeoffs |
| `bond_{angle,length}_violation_fractions_<group>.csv` | Ensemble-member geometry fractions, joined on original CIF path and `model_n` |
| `min_altloc_rmsd_results_{boltz_protenix,rf3_fixed,protenix_mini,protenix_tiny}.csv` | Ensemble minima to A/B; no chiral-off source |
| `classify_altloc_regions.csv` | Backbone-lDDT conformational classifications; single path declaration is `style.CLASSIFICATION_CSV` |
| `lddt_assignment_counts.csv` | Training-cluster assignments to reference states A/B |
| `Sampleworks_initial_structures.csv` | Deposition/training cutoff status; missing pre-2018 `1VME`, `4OLE`, `5SOP`, `6B8X` are filled as Pre |
| `rscc_results_5mhx_recyc_diff.csv` | 5MHX recycles × diffusion-step grid |
| `rscc_results_6b8x_ensemble_size.csv` | 6B8X ensemble-size sweep |

Bond outlier-detail and within-RMSD files have no active figure. They are not
substitutes for the summary inputs above. Do not delete the user-retained
`boltz2_msa_depth.tsv~`.

## Shared rules

- Loaders retain only the declared **50/50 A/B map occupancy**, currently matched
  by its exact serialized representation. Arbitrary occupancy support is future work.
- Canonical main models: `boltz2X-RAY`, `boltz2MD`, `protenix`, `rf3`.
  Alternatives: `protenix-mini`, `protenix-tiny`, `rf3-chiral-off`.
  Display names/colors come from `style.py`, not filenames.
- Baseline is guidance weight 0. The chosen guided condition is `OPTIMAL_GW`:
  0.1 for Boltz2/Protenix and 0.01 for RF3, including their corresponding variants.
  These are existing constants, not a new optimization.
- Per-region analyses exclude `selection_name == "all altloc"`. Existing mapping
  marks selections containing `chain_id` as "all altloc"; do not silently broaden it.
- Clustering score is
  `min(all predicted occupancy entries) * avg_silhouette_to_ref / min(map occupancies)`.
  The denominator is 0.5 here. Three/four-entry predicted arrays are retained;
  their state correspondence is unresolved.
- Training-bias B occupancy is the **second array entry**, not the minimum.
- The optional NaN-RSCC blacklist in 06/07 uses **all four main models** and is
  threaded through each notebook's metric frames. Focused recipes explicitly
  default to no protein blacklist; failures and coverage remain visible.
- Do not equate missing exports, absent condition rows, and recorded non-finite
  values. Do not infer a complete experiment grid without a manifest.

## Focused notebooks

### 01 — baseline distributions (Figure 2C)

`01_baseline.py`, prepared by `baseline.py`.

- Read/fingerprint the five RSCC and corresponding five lDDT sources without writes.
- Default RSCC observation: the NaN-skipping region mean per `(protein, model)`
  at baseline, including the existing whole-selection rows. The explicitly labeled
  region alternative excludes "all altloc" and changes population weighting.
- Clustering observations are finite per-region baseline scores.
- Exact ECDF vertices retain ties and stable observation IDs. A repeated first x
  anchors proportion 0; that anchor is not an observation. No KDE alternative.
- Main models start visible; variants remain available in the shared native legend.
- Fixed reconstruction SVG is main-model-only. Observation/curve/contributor CSVs
  retain all loaded models and expose original CSV record ordinals.
- Versioned artifacts: `figures/fig2c/`.

### 02 — guidance versus geometry (Figure 3B)

`02_guidance.py`, prepared by `guidance.py`.

- Read/fingerprint five RSCC plus five each clashscore/angle/length sources.
  Disambiguate variants after path tagging. Bond member joins are one-to-one;
  unmatched members are reported, not fabricated or paired by renumbered IDs.
- Per protein/model/weight: mean each metric across its member/region rows, pair
  with baseline, subtract. Inner-join geometry deltas; left-join RSCC.
- Geometry penalty is the row-wise NaN-skipping mean of delta clashscore,
  100 × delta bond-angle fraction, and 100 × delta bond-length fraction.
- Summaries are separate marginal medians and 25th/75th percentiles across
  protein rows. Require at least three protein rows; report finite counts per
  metric. Do not substitute a common complete-case cohort.
- Both renderers consume prepared symmetric-log coordinates (threshold/scale 1,
  base 10); hover/tables retain raw penalties. Marker area scales within a model.
  Black chosen-weight outlines are fixed constants.
- All seven models are loaded. Native legend subsets/isolation replace the old
  aggregate overlay and per-model aggregate spiral renderers.
- Fixed SVG is main-model-only; summary/protein CSVs retain all loaded models.
  Artifacts: `figures/fig3b/`.

### 03/04 — RSCC and clustering improvement (Figures 3C/D)

`03_rscc_improvement.py` and `04_clustering_improvement.py` share
`region_improvement.py`.

- Read only the two main-model max-RMSD RSCC files (3C), or the two main-model
  lDDT files (3D). Require eight-member ensembles.
- Pair baseline/chosen-guidance rows on `(protein, selection, model_name)`.
  Reject missing identities and duplicate pairing keys.
- Group each region across main models: x = baseline median; y = guided median
  minus baseline median, **not median of paired improvements**.
- Require three paired model **rows**, including non-finite rows. Report finite
  contributors separately. Draw only finite resulting coordinates.
- Green: top three positive improvements. Red: all decreases. Orange: take
  five (3C) or three (3D) largest `1 - guided_median` candidates, then remove
  green/red candidates without backfilling. Blue overlays: `6B8X 181-183`,
  `5MHX 162-164`, `5SOP 131-133`, including overlaps.
- Lower marginal Gaussian KDEs use independently finite marginal values among
  count-eligible regions: Scott bandwidth × 0.5, cut 3, 200 support values.
  Degenerate marginals are explicitly unavailable.
- Reference guides: no change `y=0` and guided-score ceiling `y=1-x`; x range [0,1].
- Contributions retain original CSV record ordinals; 3D also retains occupancy/
  silhouette components. No scientific recomputation in either renderer.
- Legend toggles display categories, not the contributing-model population.
  Artifacts: `figures/fig3c/` and `figures/fig3d/`.

### 05 — minimum RMSD to each altloc (Figure 3E)

`05_min_rmsd.py`, prepared by `min_rmsd.py`.

- Read two main-model RMSD exports with failure retention enabled.
- Input A/B values are already separate minima over each model's eight-member
  ensemble. The best member can differ between A and B.
- Select the chosen guided weight per model; exclude "all altloc"; reject
  duplicate selected identities and negative RMSD.
- A model contributes only if **both** A/B values are finite. Require three
  contributing main models per region.
- Default across-model **median** reconstructs the original analysis. Explicit
  **mean** is the caption-described alternative, not a verified historical recipe.
- Full-range and 0–3 Å equal-aspect panels use identical prepared coordinates.
  Blue: `6B8X 181-183`, `5MHX 162-164`, `5SOP 131-133`.
  Orange: `2QCV 2-4`, `4OLE 424-426`, `6RP1 5`.
  Red, zoom only: the curated eight `7Z0E` regions in `min_rmsd.RED_REGIONS`.
  No 6NI5/6NI6 rings in this aggregate; those remain in per-model supplements.
- Coverage distinguishes missing chosen conditions and non-finite pairs.
  Filenames and recipes identify the chosen reduction. Artifacts: `figures/fig3e/`.

### Linked inspection and export scopes

Persistent Plotly `FigureWidget` instances are displayed through `mo.ui.anywidget`.
Plot selection filters the primary table by stable rendered observation IDs.
Table callbacks add black inspection outlines in place without resetting zoom/
legend state. Density vertices and empty anchors are not observations.
Figure 3E outlines both rendered copies. Hidden models stay hidden.
Clear plot selection restores the full table; clearing rows restores marker styles.
Changing analysis inputs creates a fresh widget.

Toolbar SVGs reflect the visible view, including inspection. Fixed reconstruction
downloads do not depend on transient selection/legend state. Circled-region
reconstructions contain no fixed segment text; ZIPs bundle SVG, a companion CSV
of every circled region/categories/colors/values, and provenance. Separate data
downloads expose prepared coordinates, contributors, coverage, and source hashes.

## Remaining paper and supplementary analyses

### 06 — fit-quality supplements and summaries

`06_fit_supplements.py` retains analyses not replaced by 01–05:

| Family | Population and calculation | Output under `figures/` |
|---|---|---|
| Per-model RSCC/clustering residuals | Baseline/chosen-guidance pairs per region. Clustering rings mirror each model's RSCC-selected regions; hover identifies individual observations, including coincident points | `rscc_scatterplots/{rscc,clustering_score}_scatter_<model>.svg` |
| Per-model RMSD | Six available models, including mini/tiny. Original A/B ensemble minima at chosen weights, no across-model reduction; full-range + 0–3 Å panels. Preserve per-panel extra rings including 6NI5/6NI6 | `rmsd_scatterplots/min_rmsd_scatter_<model>.svg` |
| RSCC versus clustering | Chosen guidance, per-protein marginal means or per-selection pairs; inner-join the two metric populations | `cross_metric_scatterplots/rscc_vs_lddt[_per_selection]_<model>.svg` |
| Conformational-class grid | **Base**, not max-RMSD, per-region RSCC; baseline/chosen conditions; narrowest-containing protein/chain interval class. Unparseable/unmatched selections are reported and dropped; domain shift comes through the same join | `rscc_scatterplots/rscc_by_class_combined[_variant].svg` |

The class grid follows Main/Protenix/RF3 model groups and sizes to their members.
Per-model panels follow the model dropdown, not that group control. Alternatives
use non-overwriting output suffixes. Fixed segment text and the merged-label
expansion table are replaced by hover and linked observation tables.

Summary tables are printed, not saved; use per-region/model baseline/chosen pairs:

- Inclusive threshold fractions: RSCC ≥ .8/.9, clustering ≥ .7/.5.
- Single-residue enrichment among baseline RSCC ≥ .8 versus all baseline regions.
- Strict RSCC-up/clustering-down tradeoff and the strict reverse, on the inner
  metric intersection. Collapse means chosen minimum predicted occupancy equals 0;
  report silhouette components for collapsed cases.
- RSCC-up with chosen RMSD to **both** A/B strictly < 1 Å; RSCC-up with either
  RMSD strictly increasing, on the RSCC/RMSD intersection.
- Model and missingness coverage follows actual joins; do not hardcode historical
  claims that mini/tiny lack lDDT data. Chiral-off has no RMSD source.

### 07 — individual-protein tradeoffs and training bias

`07_training_tradeoffs.py` keeps:

- Individual-protein guidance trajectories (old IDs 12–19), not aggregate medians.
  Geometry/delta preparation retains the existing per-protein means and baseline
  joins. Lines run from the baseline origin through increasing guidance weight.
  Hover and linked rows identify proteins/conditions; no fixed protein labels.
  Geometry uses the same base-10 symlog display transform as Figure 3B, with raw
  values retained for inspection. Line colors indicate guidance weight.
- Training-bias IDs 29–36, using `build_training_bias_df`, the five lDDT sources,
  training-cluster counts, and deposition cutoffs. Baseline and chosen guided
  conditions use the same keys; variant groups are disambiguated at loading.
- Training x = fraction of reference-cluster assignments to B (`1 - frac_a_cluster`);
  y = second predicted occupancy entry. Missing predictions/cluster keys are
  reported and excluded by the required join/finite filters, not filled.
- 29: baseline/guided small-multiple scatters. 34: overlaid calibration with
  quantile-bin conditional means, slope/r and marginal histograms; 34L is its legend.
  35: separate baseline/guided calibration rows with mean ± SD summaries.
  36: data-free calibration explainer.
- 30: calibration restricted to evaluated training clusters ≤ 5 (including 5).
  31: Post-cutoff calibration. 32: cluster occupancy strip for **A = B = 0**.
  33: one-sided nonempty A-only/B-only cluster strips.
  Undefined training x is not zero. Figure 32 currently also includes rows with
  positive evaluated counts but zero assignments; despite its historical “empty”
  label, its filter does not require zero evaluated members. This remains an
  explicit scientific review gate, not a cleanup correction.

Outputs retain `12`–`19` spiral and `29`–`36` training-bias filenames; group suffixes
are `_protenix_alts` and `_rf3_alts`. Per-model spirals use their own dropdown.
Training-bias redesign alternatives remain until publication-panel adjudication.

### 08 — sampling budgets

`08_sampling_budgets.py` retains the separate supplementary sweeps:

- **5MHX** recycles `{1,3,10}` × diffusion steps `{50,100,200,400}`.
  Recycles panel holds steps at 200; steps panel holds recycles at production
  defaults (Boltz2 3, Protenix/RF3 10).
  A default max-RMSD-run cell is included only when the **loaded finite grid**
  has no row at that `(model, recycles, steps)` key; finite grid cells win overlaps.
  All-NaN exported grid cells can therefore fall back to the default run. This
  historical rule is not an assessment of which evaluations were attempted.
  No row-wise backfill.
  Show mean ± sample SD across selections, not across configurations.
  Output: `recycles_steps/5mhx_rscc_vs_recycles_steps.svg`.
- **6B8X** ensemble sizes `{1,2,4,8,16}` at chosen guidance, one source and no
  default-point merge. Mean ± sample SD across the 32 selections; evenly spaced
  log2 x positions with dodged capped error bars.
  Output: `ensemble_size/6b8x_rscc_vs_ensemble_size.svg`.

These historical budget columns are not yet a new upstream export contract.
No inference or structural evaluation runs in any of these figure notebooks.

### Plotly workflow in 06–08

All retained figure families render as Plotly vectors, not embedded static images.
Final coordinates, jittered display coordinates, raw observations, and computed
summary curves are prepared once; exports consume that same figure specification.
Calibration bins/regressions use raw fractions, never display jitter.
One-sided-strip jitter is now deterministic (seed 0); it affects presentation only.

Each chart is a persistent widget with its own stable-ID inspection table at the
bottom of the notebook. A control change creates a fresh chart. Class-grid shadows,
density vertices, regression lines, histograms, and uncertainty bands are not
observations. Budget summary points are selectable; contributor rows are associated
by the summary point ID. Table outlines never enter fixed reconstructions.

Explicit fixed-export ZIPs contain SVG, point/summaries CSVs, highlighted-region
CSV where applicable, budget contributor CSV where applicable, and the complete
figure specification. They require Kaleido and host-side Chrome. The recipes state
their prepared-value scope; unlike 01–05, they do not claim fingerprinted atomic
input snapshots or verified historical paper provenance.
