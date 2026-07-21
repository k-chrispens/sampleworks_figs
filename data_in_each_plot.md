# Data in Each Plot

This document describes, for every figure produced by the analysis notebooks, exactly what data goes in, how it is filtered and transformed, and what is left out.

---

## Data Source Inventory

> **Maintenance note (canonical inventory — keep in sync).** This section is the reference source for *every* CSV in the project: what it holds, whether it is used, and where. The tables below were last audited **2026-05-29**. **If any CSV listed here is deleted, renamed, or added, update this section** — otherwise the "Status" column silently goes stale. The per-figure "Source data" subsections elsewhere in this document describe individual figures; this section is the project-wide index.
>
> **2026-05-29 cleanup.** Deleted the Unused raw-detail files (`bond_angle_outliers.csv`, `bond_length_outliers.csv`, `initial_40_tortoize_residues.csv`, `initial_40_tortoize_protein_stats.csv` — ~11.2 GB), the Superseded small files (`clashscore_w_baseline.csv`, `clashscore_metrics_unpatched_cifs.csv`, `boltz2_msa_depth.csv`), all 7 `*.csv.bak` snapshots, and the `initial40results_backup/` directory. The `rscc_results.csv.bak` was confirmed reproducible first: all 66 of its per-region selections are present in `small-altlocs-classified-rmsd.csv` (40 proteins, 0 missing), and its 40 per-protein "all altloc" aggregates are just the union of each protein's regions. Kept: all Core files, the 6 Rebuild-source CSVs, the 4 Orphaned protenix-mini/tiny files, and `combined_occ_sweep_results/` (still the only rebuild source for the geometry Core files).
>
> **2026-06-01 update (RSCC rewiring — supersedes parts of the inventory below).** The RSCC figures now load the more-complete **max-RMSD-region** CSVs instead of the merged `rscc_results.csv`. New **Core** RSCC source (base directory, git-tracked), concatenated in-memory by `load_rscc_max_rmsd`: `rscc_results_max_rmsd_boltz_protenix.csv` (boltz2 MD/XRD + protenix), `rscc_results_max_rmsd_rf3_fixed.csv` (rf3), `rscc_results_max_rmsd_protenix_mini.csv` + `_tiny.csv` (protenix-mini/-tiny via `model_name_override`), and `rscc_results_rf3_chiral_off.csv` (rf3-chiral-off — no max-RMSD variant exists, so its non-max CSV is the source). This **closes the protenix-mini/-tiny RSCC coverage gap** described below. Consequences: (1) `initial40results/rscc_results.csv` is now **Core only for the nb-02 by-class grid** (region-containment classification) — every other RSCC figure (nb 02 + 03) uses the max-RMSD frame; the domain-shift column still uses `rscc_results_domain_shift.csv`. (2) `get_rscc_nan_proteins` now derives from `load_rscc_max_rmsd`. (3) Both merge scripts (`merge_combined_into_initial.py`, `merge_new_rscc_lddt.py`) were **deleted** — the RSCC "Rebuild source" rows below no longer have a consumer. (4) The superseded non-max per-model-group RSCC CSVs (`rscc_results_boltz_protenix.csv`, `rscc_results_rf3_fixed.csv`, `rscc_results_protenix_mini.csv`, `rscc_results_protenix_tiny.csv`) were renamed to `*.csv.bak`. Per-model coverage caveat (pre-existing upstream gaps, **identical in base and max_rmsd — no regression from this switch**): 33 of 40 proteins appear for all 7 models. The protenix family (protenix / -mini / -tiny) has no predictions for `2I6H`, `4OLE`, `6RP1`, `7Z0E`, `8Z76` (35/40); boltz2-XRD/MD has none for `6NI5`/`6NI6` (38/40); rf3 / rf3-chiral-off cover all 40.
>
> **2026-06-04 reorg (results/ rewiring — supersedes the path columns below).** All figure inputs now live in **`results/`**; `initial40results/` and the old top-level CSVs are gone. The four formerly-monolithic sources are now **five per-model-group files each, concatenated in-memory inline in the notebook cells** (hardcoded `results/…` paths, no `RESULTS_DIR` constant, no `load_*_all()` helpers): (1) **LDDT** — `load_lddt("results/lddt_results_*.csv", model_name_override=…)` in nb 02 + 03; (2) **base RSCC** (by-class grid only) — five `results/rscc_results_*.csv` via `load_rscc`; (3) **clashscore** — `load_and_tag_path_csv("results/clashscore_metrics_*.csv")`; (4) **bond violations** — `results/bond_{angle,length}_violation_fractions_*.csv`. All four take the same per-model `model_name_override` for protenix-mini/-tiny/rf3-chiral-off; the geometry (clashscore/bond) overrides are applied **after** path-tagging because the raw run-dir collapses mini+tiny to `protenix_X-RAY_DIFFRACTION` and chiral-off to `rf3_X-RAY_DIFFRACTION`. **Classification** now reads `results/classify_altloc_regions.csv` (backbone-lDDT scheme) via the single swap-point constant `style.CLASSIFICATION_CSV`, replacing `small-altlocs-classified-rmsd.csv`. **Domain-shift** is folded into the by-class per-region concat (the lone `2e7z` region survives containment), retiring `rscc_results_domain_shift.csv`. The **`model_group` dropdown** now lives in nb 03 too and drives the per-model-colored figures in both notebooks (baseline distributions + by-class grid in 02; the Figure-20 overlay spiral + the training-bias suite in 03); non-main saves get the `_protenix_alts` / `_rf3_alts` suffix. The genuinely-new-but-unplotted families (`bond_{length,angle}_outliers_*`, `max_rmsd_within_rmsd_*`) get no loader/figure in this pass. (`rscc_results_6b8x_ensemble_size` was since wired into nb 04 — see the Ensemble-Size figure below — and is loaded via `load_rscc`, which already keeps `ensemble_size`.)

### Structural facts

- **`.gitignore` excludes `initial40results*/` and `combined_occ_sweep_results/`.** So the base-directory CSVs are the **git-tracked source inputs**, while everything in `initial40results/` (the derived/merged products + multi-GB raw detail) and the two backup/source dirs are **local-only** and not version-controlled.
- **The two merge scripts were removed 2026-06-01** (`merge_combined_into_initial.py`, `merge_new_rscc_lddt.py`). RSCC is now concatenated in-memory by `load_rscc_max_rmsd` (no on-disk merge step). Provenance of the still-present derived files: the geometry Core files (`clashscore_metrics.csv` + the bond-violation-fraction files) were last written by the old `merge_combined_into_initial.py` (so they still carry the `protenix-mini`/`-tiny`/`rf3-fixed`/`rf3-chiral-off` variants) and are now frozen artifacts; `rscc_results.csv` / `lddt_results.csv` were last rebuilt by `merge_new_rscc_lddt.py` and are likewise kept as-is.

**Status legend:** **Core** = loaded directly by a notebook, feeds figures · **Rebuild source** = not loaded by notebooks but consumed by a merge script to build a Core file · **Orphaned** = present but neither loaded nor consumed by the current pipeline (contains usable data) · **Superseded** = replaced by another file the pipeline now uses · **Unused** = referenced nowhere (raw detail or exploratory).

### `initial40results/` (local-only, gitignored)

| File | Size | What it actually holds | Status — used where |
|---|---|---|---|
| `rscc_results.csv` | 68 MB | Per-altloc-region RSCC, all occupancies/gws → boltz2 MD/XRAY, protenix, rf3, rf3-chiral-off (no protenix-mini/tiny) | **Core — by-class grid only (as of 2026-06-01)** — `load_rscc`; nb 02 by-class grid. Every other RSCC figure now uses `load_rscc_max_rmsd` (see 2026-06-01 note) |
| `lddt_results.csv` | 53 MB | lDDT silhouette + predicted occupancies; same model set as above | **Core** — `load_lddt`; nb 02 + 03 |
| `rscc_results_domain_shift.csv` | 0.1 MB | Whole-domain RSCC, domain-shift class only (clean extract from legacy `.bak`, `rf3-old` dropped) | **Core** — nb 02 by-class grid (domain-shift column) |
| `clashscore_metrics.csv` | 11.5 MB | Per-model clashscore / num_clashes (path-encoded). From the *old* merge → still carries protenix-mini/tiny/rf3-fixed/rf3-chiral-off | **Core** — nb 03 spiral |
| `bond_angle_violation_fractions.csv` | 11.6 MB | Aggregated bond-angle outlier fraction per model | **Core** — nb 03 spiral |
| `bond_length_violation_fractions.csv` | 11.9 MB | Aggregated bond-length outlier fraction per model | **Core** — nb 03 spiral |
| `Sampleworks_initial_structures.csv` | 2 KB | Per-protein deposit date, Pre/Post-cutoff flag, resolution, Rfree, cluster size | **Core** — `load_deposition_cutoffs` |

*Deleted 2026-05-29 (see cleanup note above):* `bond_angle_outliers.csv`, `bond_length_outliers.csv`, `initial_40_tortoize_residues.csv`, `initial_40_tortoize_protein_stats.csv` (Unused raw per-atom/per-residue detail behind the aggregated fraction files; no figure touched tortoize); `clashscore_metrics_unpatched_cifs.csv`, `clashscore_w_baseline.csv` (Superseded by `clashscore_metrics.csv`); `boltz2_msa_depth.csv` (Unused exploratory); and all 7 `*.csv.bak` snapshots (Superseded pre-2026-04-24 `rscc`/`lddt` backups). The geometry detail is regenerable from `combined_occ_sweep_results/` if ever needed.

### Base directory (git-tracked source inputs)

| File | What it is | Status — used where |
|---|---|---|
| `small-altlocs-classified-rmsd.csv` | Region classification by **RMSD** (large_loop 450 / small_loop 274 / side_chain 66 / domain_shift 1) | **Core** — `load_classification`; nb 02 |
| `min_altloc_rmsd_results_boltz_protenix.csv` | Min RMSD to each altloc — boltz2 XRD/MD, protenix | **Core** — `load_min_altloc_rmsd`; nb 02 |
| `min_altloc_rmsd_results_protenix_mini.csv` / `_tiny.csv` | Same, protenix-mini / -tiny (`model_name_override`) | **Core** — nb 02 (only place mini/tiny get a figure) |
| `min_altloc_rmsd_results_rf3_fixed.csv` | Same, rf3 | **Core** — nb 02 |
| `lddt_assignment_counts.csv` | Pre-cutoff cluster members assigned to altloc A/B by **lDDT** | **Core** — `load_training_cluster_counts`; nb 03 training bias |
| `rscc_results_5mhx_recyc_diff.csv` | **5MHX** recycles × diffusion-steps grid-search RSCC (recycles `{1,3,10}` × steps `{50,100,200,400}`, per model, at optimal gw) | **Core** — `load_rscc_recyc_diff`; nb 04 |
| `rscc_results_6b8x_ensemble_size.csv` | **6B8X** ensemble-size sweep RSCC (sizes `{1,2,4,8,16}`, per model, at optimal gw, 32 selections) | **Core** — `load_rscc`; nb 04 |
| `rscc_results_other_models.csv` | boltz2 + protenix RSCC (raw run) | **Rebuild source** → `rscc_results.csv` |
| `rscc_results_rf3_fixed.csv` | rf3 RSCC (raw run) | **Rebuild source** → `rscc_results.csv` |
| `rscc_results_rf3_fixed_chiral_off.csv` | rf3-chiral-off RSCC | **Rebuild source** → `rscc_results.csv` |
| `lddt_results_other_models.csv` | boltz2 + protenix lDDT | **Rebuild source** → `lddt_results.csv` |
| `lddt_results_rf3_fixed.csv` | rf3 lDDT | **Rebuild source** → `lddt_results.csv` |
| `lddt_results_rf3_fixed_chiral_off.csv` | rf3-chiral-off lDDT | **Rebuild source** → `lddt_results.csv` |
| `rscc_results_protenix_mini.csv` / `_tiny.csv` | protenix-mini/-tiny **RSCC** (3,081 / 4,290 non-null at 50/50, 35 proteins) | **Orphaned** — not in merge sources, never loaded |
| `lddt_results_protenix_mini.csv` / `_tiny.csv` | protenix-mini/-tiny **lDDT** (≈4,277 non-null at 50/50, 35 proteins) | **Orphaned** — docstring calls them an "older run," not merged |
| `classify_altloc_regions.csv` | Same 791 regions classified by **backbone lDDT** instead of RMSD | **Superseded** by the RMSD version above |
| `max_rmsd_within_rmsd_2021-09-30.csv` / `_1A.csv` | Cluster assignment by **RMSD threshold** (2 Å / 1 Å) | **Superseded** — RMSD twin of `lddt_assignment_counts.csv` |
| `classify_altloc_selections.csv` | Upstream pipeline metadata: selection strings, map patterns | **Unused** — provenance only |
| `classify_input_proteins.csv` | Upstream pipeline metadata: input structure/density paths | **Unused** — provenance only |

### Two cross-cutting patterns

1. **Every classification/assignment file has an lDDT twin and an RMSD twin, and the notebooks mix which they use.** Region classification uses **RMSD** (`small-altlocs-classified-rmsd.csv`) and ignores the lDDT twin (`classify_altloc_regions.csv`); training-cluster assignment uses **lDDT** (`lddt_assignment_counts.csv`) and ignores the RMSD twins (`max_rmsd_within_rmsd_*.csv`). The unused twins are the road-not-taken metric for each — archival, not junk.
2. **Model-coverage gap (not just unused files).** The current `rscc_results.csv` / `lddt_results.csv` deliberately omit protenix-mini/tiny, but the geometry CSVs (clashscore + bond, from the old merge) *and* the min-RMSD CSVs include them. Net: **protenix-mini/-tiny appear in spiral and min-RMSD figures but cannot appear in any RSCC or LDDT figure** (the per-model dropdown is empty for them; aggregated figures hide this by dropping `MODELS_ALT`). The four **Orphaned** base-dir files hold the data that would close this gap if folded into `merge_new_rscc_lddt.py`'s source lists. **(RSCC half closed 2026-06-01:** the rewiring to `load_rscc_max_rmsd` — see the update note above — pulls protenix-mini/-tiny RSCC from the max-RMSD CSVs, so they now appear in the per-model RSCC scatter. LDDT still omits them.)**

---

## Global Filters (applied at data-loading time)

These filters are baked into the data-loading functions in `loaders.py` and apply to **every** downstream plot that uses that loader. Any data excluded here never reaches any figure.

### 1. Occupancy filter — path-tagged CSVs (`load_and_tag_path_csv` / `tag_path_dataframe`)

Used by: clashscore, bond angle, and bond length data.

- The file path of each row is parsed with a regex to extract `occ_spec` (e.g. `0.5occA_0.5occB`).
- **Only rows where `occ_spec == "0.5occA_0.5occB"` are kept.**
- All other occupancy configurations (e.g. `0.25occA_0.75occB`, `1.0occA`) are discarded.

### 2. Occupancy filter — RSCC data (`load_rscc`)

- The RSCC CSV has a string column `altloc_occupancies`.
- **Only rows where `altloc_occupancies == "{'A': 0.5, 'B': 0.5}"` are kept.**
- All other occupancy splits are discarded.

### 3. LDDT data (`load_lddt`)

- The LDDT CSV also has a string column `altloc_occupancies`; **only rows where `altloc_occupancies == "{'A': 0.5, 'B': 0.5}"` are kept** (same filter as RSCC).
- Columns retained: `protein`, `model_name`, `guidance_weight`, `selection`, `occupancies`, `avg_silhouette`, `avg_silhouette_to_ref`, `min_predicted_occupancy`, `min_map_occupancy`, `clustering_score`.
- **Clustering score** is computed at load time as
  ``clustering_score = min(predicted_occupancies) * avg_silhouette_to_ref / min(map_occupancies)``,
  where ``predicted_occupancies`` is parsed from the ``occupancies`` column
  (numpy array repr, may have 2+ entries) and ``map_occupancies`` is parsed
  from the ``altloc_occupancies`` column (Python dict repr). Because every
  loader output is filtered to ``altloc_occupancies == "{'A': 0.5, 'B': 0.5}"``,
  ``min_map_occupancy`` is always ``0.5`` in the current analysis; the
  divisor is kept explicit so non-50/50 splits can be supported without
  code changes downstream.

### 4. Model name aliasing (all loaders)

After loading, model names are canonicalized via `MODEL_NAME_ALIAS`:

| Raw name in CSV | Canonical name used in plots |
|---|---|
| `boltz2_X-RAY_DIFFRACTION` | `boltz2X-RAY` |
| `boltz2_MD` | `boltz2MD` |
| `boltz2xray` | `boltz2X-RAY` |
| `boltz2md` | `boltz2MD` |
| `rf3-fixed` | `rf3` |

For Boltz2 models in the RSCC and LDDT loaders, the model name is constructed by concatenating the `model` and `method` columns (e.g. `boltz2` + `X-RAY_DIFFRACTION`), then aliased.

### 5. Guidance strength bucketing (all loaders)

Every loader adds a `guidance_strength` column mapped from the numeric `guidance_weight`:

| guidance_weight | guidance_strength |
|---|---|
| 0.0 | baseline |
| 0.005 | weak |
| 0.01 | weak |
| 0.02 | medium |
| 0.035 | medium |
| 0.05 | strong |
| 0.1 | weak |
| 0.2 | medium |
| 0.35 | medium |
| 0.5 | strong |

### 6. Selection name mapping

When a `selection` column is present, a `selection_name` column is derived: if the selection string contains `"chain_id"`, it is labeled `"all altloc"`; otherwise the raw selection string is used as-is.

### 7. Model groups

Several figures are rendered in three variants via a UI dropdown. The model groups are:

- **Main models**: `boltz2X-RAY`, `boltz2MD`, `protenix`, `rf3` (i.e. everything NOT in `MODELS_ALT`)
- **MODELS_ALT** (alternative/variant models): `protenix-mini`, `protenix-tiny`, `rf3-chiral-off`
- **Protenix group**: `protenix`, `protenix-mini`, `protenix-tiny`
- **RF3 group**: `rf3`, `rf3-chiral-off`

---

## Notebook 1: Fit Quality (`02_fit_quality.py`)

### Source data

| File (all under `results/`) | Key columns |
|---|---|
| `rscc_results_max_rmsd_*.csv` + `rscc_results_rf3_chiral_off.csv` (concatenated by `load_rscc_max_rmsd`) | `protein`, `model`, `method`, `altloc_occupancies`, `guidance_weight`, `rscc`, `selection`, `ensemble_size` — **primary RSCC source** for every nb-02 figure except the by-class grid |
| `rscc_results_{boltz_protenix,rf3_fixed,protenix_mini,protenix_tiny,rf3_chiral_off}.csv` | Same schema — **by-class grid only**, concatenated inline via `load_rscc` (mini/-tiny/chiral-off via `model_name_override`); region-containment classification |
| `lddt_results_{boltz_protenix,rf3_fixed,protenix_mini,protenix_tiny,rf3_chiral_off}.csv` | `protein`, `model`, `method`, `altloc_occupancies`, `guidance_weight`, `selection`, `occupancies`, `avg_silhouette`, `avg_silhouette_to_ref` — concatenated inline via `load_lddt(..., model_name_override=…)` |
| `classify_altloc_regions.csv` | `protein`, `selection`, `chain`, `start_res`, `end_res`, `classification` (backbone-lDDT scheme; via `load_classification(CLASSIFICATION_CSV)`). The lone `domain_shift` region (`2e7z`) is folded into the by-class per-region concat — no separate domain-shift CSV |

**RSCC aggregation step:** After loading, RSCC rows are grouped by `(protein, model_name, guidance_weight)` and averaged across all per-altloc-region selections (`.mean()`), yielding one RSCC value per (protein, model, gw) triple.

**Shared baseline-vs-optimal frames:** A helper `baseline_vs_optimal(region, value_cols)` pairs each metric's `guidance_weight == 0.0` value with its value at the model's `OPTIMAL_GW`, per `(protein, selection, model_name)`, naming columns `<col>_baseline` / `<col>_optimal`. Two frames built this way drive the per-(protein, selection, model) scatters and summary tables: `rscc_scatter_df` (RSCC, `selection_name != "all altloc"`) and `cs_scatter_df` (`clustering_score` + `avg_silhouette_to_ref` + `min_predicted_occupancy`, baseline/optimal, plus `delta_clustering_score`). (2026-06-03)

---

### Figures 04–09 — Combined Baseline Distributions (`04_09_baseline_distributions[_variant].svg`)

**What is plotted:** A two-panel figure: RSCC baseline distribution on top, LDDT-to-reference baseline distribution on the bottom, both with all main models overlaid on a single axis (one color per model). A single shared model legend in the upper-right of the top panel. Supports a `cumulative` switch (KDE or ECDF cumulative views via `cumulative_style`).

**Data included:**
- RSCC panel: aggregated RSCC at `guidance_strength == 'baseline'`, restricted to the current dropdown model group.
- LDDT panel: LDDT at `guidance_strength == 'baseline'` with `selection_name != 'all altloc'`, **also restricted to the current dropdown model group** (both panels now honor the dropdown via `filter_by_model_group`). Non-main groups save with the `_protenix_alts` / `_rf3_alts` suffix.

**What is excluded:**
- All occupancy configurations other than 50/50 (global loaders).
- All guided rows.
- For LDDT: chain-level "all altloc" selections.

---

### RSCC Scatter Plots — Baseline vs Optimal Guidance (per model)

**What is plotted:** Scatter plot with one point per **(protein, selection)** pair for the selected model. X-axis = baseline RSCC, Y-axis = improvement (residual = guided − baseline). Reference lines `y = 0` (no change) and `y = 1 − x` (max-improvement ceiling). X axis pinned to `[0, 1]`; Y pinned to `[min(0, residual_min) - 0.02, 1.0]`. Marginal KDE density panel below.

**Highlighting** (capped so the figure stays readable):
- Top-3 most-improved selections: green ring.
- Every decreased selection (negative `delta_rscc`): red ring.
- 5 selections furthest below the `y = 1 − x` ceiling: orange ring.
- Loop-of-interest selections: blue ring. Currently `6B8X 181-183` and `5MHX 162-164` (the altloc-resolved residues inside each protein's most flexible large loop).
- Labels use the compact form `"<PDB> <start>-<end>"`.
- **Ring style:** every highlight is an open (hollow) circle sized larger than the data marker, so the ring sits with a visible gap around its point instead of hugging it — the same convention as the min-RMSD scatters. Shared across the per-model RSCC and clustering panels, the cross-metric per-selection scatter, and the aggregated RSCC/clustering scatters.

**Optimal guidance weights:**

| Model | Optimal gw |
|---|---|
| boltz2X-RAY | 0.1 |
| boltz2MD | 0.1 |
| protenix | 0.1 |
| protenix-mini | 0.1 |
| protenix-tiny | 0.1 |
| rf3 | 0.01 |
| rf3-chiral-off | 0.01 |

**Data construction:**
1. Load RSCC filtered to `selection_name != "all altloc"` (drop chain-level rows).
2. Baseline: rows with `guidance_weight == 0.0`, keeping `(protein, selection, model_name, rscc)` renamed to `rscc_baseline`.
3. Optimal: for each model, rows at that model's optimal gw, same columns renamed to `rscc_optimal`.
4. Merge on `(protein, selection, model_name)`.
5. `delta_rscc = rscc_optimal - rscc_baseline`.
6. Add `label = "<PDB> <start>-<end>"` for annotations.

**Data included:**
- One point per (protein, selection) pair for the selected model, requiring both baseline and optimal-gw measurements.

**What is excluded:**
- `selection_name == "all altloc"` (chain-level rows).
- All occupancy configurations other than 50/50.
- All guidance weights other than 0.0 (baseline) and the model's designated optimal weight.
- (protein, selection) pairs missing from either slice.

---

### Clustering Score Scatter — Baseline vs Optimal Guidance (per model)

**What is plotted:** The clustering-score analog of the per-model RSCC scatter above. One point per **(protein, selection)** pair for the selected model. X-axis = baseline clustering score, Y-axis = improvement (residual = guided − baseline). Reference lines `y = 0` (no change) and `y = 1 − x` (max-improvement ceiling). X axis pinned to `[0, 1]`; Y pinned to `[min(0, residual_min) - 0.02, 1.0]`. Marginal KDE density panel below.

**Highlighting:** The labels are **pinned to the exact same selections, in the exact same colors, as that model's RSCC residual panel** (fig 4 panel C) — *not* recomputed on the clustering metric. A selection that is green (top-3 improved) in the RSCC panel is green here; red (decreased) stays red; orange (far-from-ceiling) stays orange; blue (curated) stays blue — even where the clustering values would have sorted them differently. This is what makes the two panels directly comparable label-for-label.

Mechanically: the green / red / orange sets come from `rscc_residual_label_override(rscc_scatter_df[model])`, which returns the `top` / `decreased` / `far` label-string lists that the RSCC panel's residual mode selects (top-3 by RSCC improvement; every negative-RSCC-residual selection; 5 furthest below the RSCC `y = 1 − x` ceiling). These are passed to `make_clustering_scatter(..., residual_label_override=...)`, which forwards them to `make_scatter_comparison`; that function then rings/labels those exact selections (with the same green→red→blue→orange draw order and seen-dedup, so overlap colors — e.g. a curated point that is also far-from-ceiling shows orange — resolve identically to panel C). The blue curated set is the model-agnostic `RSCC_HIGHLIGHT_PROTEINS` (`6B8X 181-183`, `5MHX 162-164`, `5SOP 131-133`, `6DUR 389-391`, `2QCV 2-4`, `6RP1 5`, and the `6NI5` / `6NI6` altloc segments), shared with `make_rscc_scatter`. Labels use the compact form `"<PDB> <start>-<end>"`.

**Coincident labels merged automatically.** On the clustering metric almost all of these selections have baseline clustering ≈ 0, so they stack on the y-axis (x≈0) and several land on the *exact* same coordinate (e.g. for RF3, `5SOP 131-133`, `6DUR 389-391`, `3T94 229`, `5MHX 162-164`, `3L4P 18` are all at `(0, 0)`), which no placement algorithm can separate. `make_clustering_scatter` passes `group_coincident_labels=True`, so `make_scatter_comparison` merges labels that round to the same `coincident_tol` (default `1e-3`, data coords) cell into one combined label — the first `coincident_max_members` (default 2) joined by `, ` plus a `(+N)` count for the rest, e.g. `6DUR 389-391, 3T94 229 (+3)`. Only the **text** is merged; each selection keeps its own colored ring. Selections with a NaN clustering score (e.g. RF3 `5RBR 164`) are dropped from text since they cannot be plotted. The old hand-maintained `suppress_label_proteins` / `CLUSTERING_NO_LABEL_PROTEINS` set is retired (the param remains as an optional manual escape hatch, default `None`). The full membership behind each `(+N)` is recoverable from the printed **Merged Clustering Labels — (+N) Expansion** table that follows the clustering panel in `02_fit_quality.py` (built from the same `coincident_groups_out` grouping `make_clustering_scatter` populates during the draw, so it cannot drift from the panels).

**Deterministic placement.** After merging, the surviving labels are placed in-plot by the same `adjust_text` the RSCC panels use (so the two metrics' panels stay visually consistent). It is made reproducible by seeding NumPy with `adjust_seed` (default `0`) and using a fixed `iter_lim` (400) with `time_lim=None` instead of the old wall-clock `time_lim` — so the figure is byte-identical run-to-run. The earlier `static_labels` fixed-offset hack is removed; the merge step (not suppression or hand-annotation) is what resolves the former exact-`(0,0)` smear. If a model's cluster ever stays too dense in-plot, the intended fallback is an opt-in right-margin leader-line gutter (not currently enabled).

Because clustering scores concentrate at low values where the RSCC selections are spread out, the merged labels still sit near the left edge, fanned off the y-axis by `adjust_text`; this is accepted as the cost of an exact label-for-label match with panel C. A selection labeled in the RSCC panel that has no clustering-score row for that model simply does not appear (it cannot be drawn).

**Optimal guidance weights:** same table as the per-model RSCC scatter above.

**Data construction:** reuses the shared `cs_scatter_df` frame (the clustering analog of `rscc_scatter_df`, built by `baseline_vs_optimal` from the in-memory LDDT concat, filtered to `selection_name != "all altloc"`), with the `label` column added in the figure cell:
1. Baseline: `clustering_score` at `guidance_weight == 0.0` → `clustering_score_baseline`.
2. Optimal: `clustering_score` at each model's optimal gw → `clustering_score_optimal`.
3. `delta_clustering_score = clustering_score_optimal - clustering_score_baseline`.
4. Add `label = "<PDB> <start>-<end>"` for annotations.

**Data included:**
- One point per (protein, selection) pair for the selected model, requiring both baseline and optimal-gw clustering measurements.

**What is excluded:**
- `selection_name == "all altloc"`, non-50/50 occupancies, and all guidance weights other than 0.0 and the model's optimal weight (inherited from the LDDT loader / `baseline_vs_optimal`).
- Models with no clustering-score coverage render a `mo.callout(kind="warn")` in place of the figure rather than an empty plot.

---

### Aggregated RSCC Scatter — All Main Models Combined

**What is plotted:** Same residual scatter as the per-model figures, but each point is one **(protein, selection)** pair collapsed across main models using the **median** (no error bars; per-model spread is shown by the per-model figures). Marginal KDE.

**Data construction:**
1. Same per-model construction as above.
2. Restrict to `model_name not in MODELS_ALT` (main models only).
3. Group by `(protein, selection)`; take median of `rscc_baseline`, `rscc_optimal`, and count of contributing models.
4. Drop groups with fewer than 3 main models contributing.

**What is excluded:**
- Alternative models (`MODELS_ALT`).
- `selection_name == "all altloc"`.
- (protein, selection) pairs with fewer than 3 main models.

---

### Aggregated Min-RMSD Scatter — Min RMSD to Altloc A vs B

**Output:** `figures/rmsd_scatterplots/min_rmsd_scatter_aggregated_all_models.svg`.

**What is plotted:** Two side-by-side panels in one SVG. Each point is one **(protein, selection)** pair; X = median across main models of `min_rmsd_to_A` (the per-protein/selection minimum RMSD across the 8-member ensemble to altloc A), Y = same for altloc B. Both panels share a dashed gray `y = x` diagonal across the visible range, italic gray corner cues (`"closer to A"` upper-left, `"closer to B"` lower-right), equal aspect, and the same point cloud / highlights:

- **Left panel (Auto-fit):** axes auto-scaled to the data range (with a 5% padding) so global outliers stay visible.
- **Right panel (Zoomed):** axes pinned to `[0, 3]` Å × `[0, 3]` Å so the dense low-RMSD cluster is readable (range still controlled by `zoom_xylim`, just no longer printed in the title). Off-scale points clip naturally — no special marker.

**Highlights (split per panel):**
- Pinned (blue, `#1565C0`): `6B8X 181-183` and `5MHX 162-164` on **both** panels (the same selections circled in the aggregated RSCC scatter, for direct cross-reference). Pinned highlights are not per-panel.
- RSCC-matched extra highlights (defined as `MIN_RMSD_EXTRA_HIGHLIGHTS` in the notebook, passed via `extra_highlight_groups` to both the aggregated and per-model figures). Colors mirror the per-model RSCC scatter (`make_rscc_scatter`):
  - Orange `#F57C00` (RSCC "far-from-ceiling" family): `2QCV 2-4`, `4OLE 424-426`, `6RP1 5` — **both panels**.
  - Blue `#1565C0` (RSCC manually-curated biological highlights): `5SOP 131-133` — **both panels**. The `6NI5` / `6NI6` altloc segments that the per-model RMSD figures circle blue are **intentionally omitted from this aggregated scatter**: the aggregated cell strips every `6NI5`/`6NI6` label out of each panel/color group of `MIN_RMSD_EXTRA_HIGHLIGHTS` before plotting, leaving the shared dict — and therefore the per-model figure — untouched.
  - Red `#D32F2F` (RSCC "decreased" family): a curated subset of `7Z0E` altloc segments — `11-13`, `18-20`, `54-56`, `82-84`, `122-124`, `181-183`, `206-208`, `214-216` — **zoom panel only**, intentionally kept off the auto-fit panel where the dense cluster crowds the global view. Intentionally a subset (not all 37 7Z0E selections) to keep the zoom panel readable; segments listed here are the ones flagged in the RSCC scatter.
- Optional `extra_highlight_groups` knob accepts either the flat form `dict[str, list[str]]` (same extras on both panels, backward-compatible) or the nested per-panel form `{"autofit": {color: [labels]}, "zoom": {color: [labels]}}`. Nested form is detected by an explicit key-set check (outer keys must be a subset of `{"autofit", "zoom"}`); a missing panel key means no extras on that panel (pinned highlights still render on both). Highlights are restricted to rows visible inside each panel before label placement so `adjust_text` doesn't fight off-axes points.

**Source data:**
- `min_altloc_rmsd_results_boltz_protenix.csv` — Boltz2-XRD, Boltz2-MD, Protenix.
- `min_altloc_rmsd_results_protenix_mini.csv` — Protenix-Mini (`model_name_override="protenix-mini"` because the CSV stores `model == "protenix"`).
- `min_altloc_rmsd_results_protenix_tiny.csv` — Protenix-Tiny (same override).
- `min_altloc_rmsd_results_rf3_fixed.csv` — RF3.

Each loader applies the standard 50/50 occupancy filter (`altloc_occupancies == "{'A': 0.5, 'B': 0.5}"`), aliases the model name through `MODEL_NAME_ALIAS`, and **drops rows where either `min_rmsd_to_A` or `min_rmsd_to_B` is NaN**, printing one warning per CSV with the dropped count.

**Data construction:**
1. Concat the four loader outputs into `min_rmsd_all`; apply the `nan_rscc_proteins` exclusion (when the toggle is on).
2. `min_rmsd_per_model_df`: drop `selection_name == "all altloc"`, inner-join on `model_name` with the `OPTIMAL_GW` table, restrict to `guidance_weight == optimal_gw`, attach `label = "<PDB> <residue range>"` via `format_residue_range`. Print a one-line audit per model when any (protein, selection) pair existed at some gw but lacked a row at the model's optimal gw.
3. Aggregated frame: restrict to main models (`model_name not in MODELS_ALT`), group by `(protein, selection, label)`, take median of both RMSDs, count contributing models, drop groups with `n < 3`.

**What is excluded:**
- Alternative models (`MODELS_ALT`) for the aggregated frame only.
- `selection_name == "all altloc"` (chain-level rows).
- All occupancy configurations other than 50/50 (loader filter).
- All guidance weights other than each model's optimal weight.
- (protein, selection) pairs with fewer than 3 main-model contributions.
- Rows with NaN `min_rmsd_to_A` or `min_rmsd_to_B` (loader filter).

---

### Per-Model Min-RMSD Scatter Plots

**Output:** `figures/rmsd_scatterplots/min_rmsd_scatter_{model_slug}.svg` per model that has a min-RMSD CSV.

**What is plotted:** Same two-panel layout as the aggregated figure but for one model at a time, drawn from `min_rmsd_per_model_df` filtered to the dropdown-selected model. Title is `"Minimum RMSD to Each Altloc: {ModelDisplay}"`. No baseline overlay; X / Y are the model's raw `min_rmsd_to_A` / `min_rmsd_to_B` at its optimal guidance weight. Passes the **full** `MIN_RMSD_EXTRA_HIGHLIGHTS` overlay (the aggregated figure strips the `6NI5`/`6NI6` segments first, so this per-model figure keeps them while the aggregated one does not): pinned blues plus the orange `2QCV 2-4` / `4OLE 424-426` / `6RP1 5` and blue `5SOP 131-133` + `6NI5`/`6NI6` segment groups appear on both panels, while the red curated `7Z0E` subset appears on the zoom panel only — selections absent from a particular model's frame, or sitting outside a panel's limits (e.g. 6NI5/6NI6 `213-215` on the zoom panel), simply don't get circled in that panel.

**Models with data:** `boltz2X-RAY`, `boltz2MD`, `protenix`, `protenix-mini`, `protenix-tiny`, `rf3` (the alternative models `protenix-mini` and `protenix-tiny` *do* appear here; only the aggregated figure drops them).

**Empty case:** When the dropdown selects a model with no min-RMSD CSV (currently only `rf3-chiral-off`), the cell renders `mo.callout("No min-RMSD data available for {ModelDisplay}.", kind="warn")` in place of the figure — no traceback, no empty axes. The `save_all_scatter` switch iterates only over models with data via `save_per_model_figures`.

**What is excluded:** same as the aggregated figure, except the `MODELS_ALT` and `n < 3` constraints don't apply (each model's own (protein, selection) pairs are plotted as-is).

---

### RSCC vs LDDT Cross-Metric Scatter Plots (per model)

**What is plotted:** Scatter comparing RSCC (x) vs clustering score (y), both at the optimal guidance weight. One point per protein (RSCC and clustering both averaged across selections). Regression line + Pearson r (upper-left) displayed.

**Data construction:**
1. RSCC at optimal gw per model (same as scatter plot above) — `rscc_all` is the per-(protein, model, gw) RSCC averaged across selections.
2. LDDT at optimal gw: filtered to `selection_name != 'all altloc'`, then for each model at its optimal gw, `clustering_score` is **averaged across all remaining selections** per protein.
3. Merged on `(protein, model_name)`.

**Data included:**
- Proteins that have both an RSCC and an LDDT measurement at the optimal gw for the selected model.

**What is excluded:**
- LDDT rows with `selection_name == 'all altloc'`.
- All non-optimal guidance weights.
- All occupancy configurations other than 50/50 (RSCC side).
- Proteins missing from either metric.

---

### RSCC vs LDDT Cross-Metric Scatter Plots — Per Selection (per model)

**Output:** `figures/cross_metric_scatterplots/rscc_vs_lddt_per_selection_{model}.svg`

**What is plotted:** The same RSCC (x) vs clustering score (y) comparison at the optimal guidance weight, but with **one point per (protein, selection)** altloc region — neither metric averaged across selections. This is the finer-grained companion to the per-protein plot above; with ~700–800 regions per model it gives a far more populated cloud and a more stable Pearson r (annotated upper-left). Regression line drawn; `6B8X 181-183` and `5MHX 162-164` circled blue + labeled (other points unlabeled to keep the dense cloud readable).

**Data construction:**
1. RSCC per `(protein, selection, model)` at each model's optimal gw, from the max-RMSD source (`selection_name != 'all altloc'`).
2. Clustering score per `(protein, selection, model)` at each model's optimal gw.
3. Merged on `(protein, selection, model_name)` (inner — keeps regions present in both metrics).

**Data included:**
- (protein, selection) regions that have both an RSCC and a clustering-score measurement at the optimal gw for the selected model.

**What is excluded:**
- LDDT/RSCC rows with `selection_name == 'all altloc'`.
- All non-optimal guidance weights.
- All occupancy configurations other than 50/50 (RSCC side).
- Regions missing from either metric (e.g. mini/tiny/chiral-off lack the pinned blue regions, so those points are simply absent there).

---

### Aggregated Clustering Score Scatter — All Main Models Combined

**What is plotted:** Same residual scatter style as the aggregated RSCC scatter, but for the LDDT-derived clustering score. One point per **(protein, selection)** pair; baseline and optimally-guided clustering score collapsed across main models using the **median**.

**Data construction:**
1. LDDT filtered to `selection_name != 'all altloc'`.
2. Baseline: `guidance_weight == 0.0` and `model_name not in MODELS_ALT`.
3. Optimal: each main model at its `OPTIMAL_GW`.
4. Merge on `(protein, selection, model_name)`, group by `(protein, selection)`, take median / count.
5. Drop groups with fewer than 3 main models.

**Highlighting:** Same as the Aggregated RSCC scatter — top-3 improved (green), top-5 most-decreased (red), top-3 furthest below ceiling (orange), `6B8X 181-183` and `5MHX 162-164` (blue).

**What is excluded:**
- Alternative models.
- LDDT rows with `selection_name == 'all altloc'`.
- Non-baseline / non-optimal guidance weights.
- (protein, selection) pairs with fewer than 3 main models.

---

### Combined By-Class Grid (`rscc_by_class_combined.svg`)

**Output:** `figures/rscc_scatterplots/rscc_by_class_combined[_variant].svg` (the `_protenix_alts` / `_rf3_alts` suffix on non-main groups).

**What is plotted:** An *N* × 4 grid: rows = the selected model group's members (Main → the 4 main models; Protenix variants → protenix / -mini / -tiny; RF3 variants → rf3 / -chiral-off), columns = conformational-change classification (`side_chain_only`, `small_loop`, `large_loop`, `domain_shift`). All panels share x and y axes. X is pinned to `[0, 1]`; y is pinned globally to `[min(0, global_min_residual) − 0.02, 1.0]` so every panel exposes the same `y = 1 − x` improvement ceiling. Each panel draws (a) a light-gray context shadow consisting of **that row's per-model pooled (baseline, residual) cloud**, (b) the per-class points overlaid in the classification color, (c) `y = 0` (no-change) and `y = 1 − x` (max-improvement ceiling) reference lines. Column titles on the top row report `CLASSIFICATION_DISPLAY[cls]` plus the unique-`(protein, selection)` count for that class **within the current group's frame** (so the counts differ between groups — e.g. Large Loop is ≈19 for Main but ≈36 for the protenix family). Each model row carries a centered bold banner above the row. One figure-level `supxlabel("Baseline RSCC")` and one `supylabel("Improvement (Guided − Baseline RSCC)")`.

**Data construction (in `02_fit_quality.py`'s by-class cell):**
1. **Single per-region concat:** concatenate the five base (non-max) `results/rscc_results_*.csv` via `load_rscc` (mini/-tiny/chiral-off via `model_name_override`), drop `selection_name == "all altloc"`, run `merge_classification_by_containment` against `results/classify_altloc_regions.csv`. **All four classes — including `domain_shift` — flow through this one path** (the lone `domain_shift` classification range, `2e7z` `chain A and resi 204-324`, contains ≈2 per-region RSCC selections that survive the containment merge). Baseline = `guidance_weight == 0.0`; Optimal = each model at its `OPTIMAL_GW`. Merge on `(protein, selection, model_name)`.
2. `filter_by_model_group(by_class_df, model_group.value)` restricts to the selected group and returns the save suffix; the group's members (canonical order) are passed as `model_order`.
3. `make_by_class_combined_grid(merged_df, model_order=…)` does the rest.

**Data included:**
- One point per (protein, selection, model) triple in the selected group's merged frame.

**What is excluded:**
- Models outside the selected group.
- `selection_name == "all altloc"`.
- All occupancy configurations other than 50/50.
- All guidance weights other than 0.0 (baseline) and each model's optimal weight.

**Note on domain-shift sourcing:** this now reads the per-altloc-region RSCC for the residues inside the domain (the old whole-domain `rscc_results_domain_shift.csv` is retired), so the Domain Shift column is sparse (n≈1 for Main, n≈2 for the protenix family). Accepted trade-off of the single-source rewiring.

---

### Summary Tables — Baseline vs Optimal (printed, not saved) (2026-06-03)

The notebook closes with a block of tabular analyses (rendered as marimo tables, **no SVG output**). All operate at the **per-(protein, selection, model)** level on the two shared delta frames `rscc_scatter_df` and `cs_scatter_df` (and their inner join `rscc_cs_delta_df`), comparing each metric's baseline (`gw == 0.0`) with its value at the model's `OPTIMAL_GW`. RSCC comes from the max-RMSD CSVs, clustering score from `lddt_results.csv`, min RMSD from `min_altloc_rmsd_results_*.csv`. `selection_name == "all altloc"` is excluded throughout.

| Table | Scope | What it reports |
|---|---|---|
| Threshold-clearing fractions | main models | fraction of predictions clearing RSCC ≥ 0.8, RSCC ≥ 0.9, clustering score ≥ 0.7, and clustering score ≥ 0.5, baseline vs optimal, over the identical point set |
| …by model | all 7 models | the same four fractions per model — columns `rscc_ge08_*` / `rscc_ge09_*` / `cs_ge07_*` / `cs_ge05_*` (variants without clustering data show NaN in the `cs_*` columns) |
| Single-residue enrichment | main models | fraction of the baseline RSCC ≥ 0.8 set that are single-residue (`resi N-N`) selections vs the full baseline universe, with the overrepresentation factor |
| RSCC↑ / clustering↓ trade-off | RSCC ∩ clustering | selections where optimal gw raises RSCC but lowers clustering score (strict), per model |
| Clustering collapse vs partial | trade-off set | within that set, how many collapse clustering to 0 (`min_predicted_occupancy_optimal == 0`) vs partially decrease, plus the `avg_silhouette_to_ref` of the collapsed cases |
| RSCC↓ / clustering↑ (reverse) | RSCC ∩ clustering | the mirror direction (optimal gw lowers RSCC, raises clustering) |
| RSCC↑ & min RMSD to both altlocs < 1.0 Å | RSCC ∩ min-RMSD | RSCC-up selections whose guided ensemble stays within 1.0 Å of both A and B |
| RSCC↑ & min RMSD to A or B increases | RSCC ∩ min-RMSD | RSCC-up selections that drift away from at least one altloc vs baseline |

**Coverage / threshold notes:**
- Clustering-score tables use the RSCC ∩ clustering inner join, which drops `protenix-mini` / `protenix-tiny` (no clustering data) but keeps `rf3-chiral-off`.
- Min-RMSD tables use the RSCC ∩ min-RMSD inner join, covering the 4 main models + `protenix-mini` / `protenix-tiny`; `rf3-chiral-off` has no min-RMSD CSV yet.
- The threshold-fraction tables report two cutoffs per metric — RSCC ≥ 0.8 / ≥ 0.9 and clustering ≥ 0.7 / ≥ 0.5, all inclusive; the trade-off increase/decrease and the min-RMSD drift / < 1.0 Å tests are strict (`>` / `<`).

---

## Notebook 2: Guidance-Weight Spiral Trajectories (`03_tradeoffs.py`)

### Source data

All four metric families are loaded (each concatenated inline from the five per-model-group `results/` files):
- `results/clashscore_metrics_{boltz_protenix,rf3_fixed,protenix_mini,protenix_tiny,rf3_chiral_off}.csv` (via `load_and_tag_path_csv`; variant `model_name` overridden after tagging — see the 2026-06-04 note).
- `results/bond_angle_violation_fractions_*.csv` (same five sources, merged with the bond-length files on `[model, model_n]` per source, then variant override).
- `results/bond_length_violation_fractions_*.csv`.
- `results/rscc_results_max_rmsd_*.csv` + `results/rscc_results_rf3_chiral_off.csv` (RSCC, concatenated by `load_rscc_max_rmsd`; sparser gw sweep than base, so rf3 spirals show only their low-gw arm).

### Tradeoff dataframe construction

1. **Load and tag** clashscores, bond stats, and RSCC.
2. **Aggregate per (protein, model, gw, guidance_strength):** mean of `clashscore`, `outlier_fraction_bond_angles`, `outlier_fraction_bond_lengths`, `rscc`.
3. **Compute deltas** for each metric via `make_delta_frame` (`delta_X = X_guided − X_baseline`).
4. **Merge** all four delta frames on `(protein, model_name, guidance_weight, guidance_strength)`. The RSCC merge uses `how="left"`.
5. **Convert bond fractions to percentage points** (`*_pct = 100 * *`).
6. **Compute geometry penalty** = unweighted mean of `delta_clashscore`, `delta_outlier_fraction_bond_lengths_pct`, `delta_outlier_fraction_bond_angles_pct`.

The final `tradeoff_df` has one row per (protein, model, guidance_weight) at every non-baseline guidance level, with columns `delta_rscc` and `geometry_penalty`.

---

### Figures 12–17, 19 — Per-Model Spiral Trajectories

**Figure numbering:** 12 = boltz2X-RAY, 13 = boltz2MD, 14 = protenix, 15 = rf3, 16 = protenix-mini, 17 = protenix-tiny, 19 = rf3-chiral-off.

**What is plotted:** Per-model trajectory plot. One line per protein traces through (delta RSCC, geometry penalty) space as guidance weight increases. Each trajectory starts at the origin (baseline). Line segments and dots colored by guidance weight using the "flare" colormap. Arrowheads indicate the direction of increasing guidance. Top-5 most extreme trajectories (by total arc length) are labeled with the protein name. Y-axis uses symmetric log scale (`symlog`, `linthresh=1.0`).

**Data included:**
- All rows in `tradeoff_df` for the selected model.

**What is excluded:**
- All occupancy configurations other than 50/50 (global filter).
- Baseline rows (they define the origin).
- `inf`/`-inf` in delta_rscc / geometry_penalty replaced with NaN and excluded from axis-limit computation.

**Shared axis limits:** X and Y limits are computed across **all models** simultaneously to enable direct visual comparison between figures 12–17 and 19.

---

### Figure 20 — Aggregated Spiral (selected model group overlaid)

**Output:** `figures/20_aggregated_spiral_all[_variant].svg` (the `_protenix_alts` / `_rf3_alts` suffix on non-main groups).

**What is plotted:** Multi-model overlay. For each (model, guidance_weight), proteins are collapsed to a single point at the **median** delta RSCC and **median** geometry penalty. IQR error bars in both X and Y. Trajectory line from origin through each guidance weight, with arrowhead. Y-axis on symlog. Fixed `ylim=[-1, 1000]`, `xlim=[-0.05, 0.25]`.

**Models included:** the selected model group — Main (`boltz2X-RAY`, `boltz2MD`, `protenix`, `rf3`), Protenix variants (`protenix`, `protenix-mini`, `protenix-tiny`), or RF3 variants (`rf3`, `rf3-chiral-off`). The model list comes from `filter_by_model_group(tradeoff_df, model_group.value)`.

**Geometry-coverage guard:** when a variant group is selected, the cell prints each model's per-gw point count (gw cells with `n ≥ 3`) and warns if any variant's gw sweep is **sparser** than its family's main model — a degenerate trajectory would otherwise be drawn silently. Current data: the variants match their family's sweep (protenix-mini/-tiny = 5 gws like protenix; rf3-chiral-off = 6 gws like rf3), so no warning fires.

**Data included:**
- Aggregated statistics from `tradeoff_df` grouped by `(model_name, guidance_weight)`, restricted to the selected group.
- Only groups with `n ≥ 3 proteins`.

**What is excluded:**
- Models outside the selected group.
- Groups with fewer than 3 proteins.

---

### Figures 21–26, 28 — Per-Model Aggregated Spirals

**Figure numbering:** 21 = boltz2X-RAY, 22 = boltz2MD, 23 = protenix, 24 = rf3, 25 = protenix-mini, 26 = protenix-tiny, 28 = rf3-chiral-off.

**What is plotted:** Single-model aggregated spiral. Same aggregation as Figure 20 but for one model. Dots colored by guidance weight (flare colormap). IQR error bars. Colorbar legend for guidance weight.

**Data included:**
- All (model, guidance_weight) groups for the selected model with `n ≥ 3 proteins`.

**What is excluded:**
- Groups with fewer than 3 proteins.

---

### Figures 29–36 — Training-Distribution Bias

**Common inputs:**
- `results/lddt_assignment_counts.csv` — pre-cutoff (PDB depositions before 2021-09-30) per-selection counts of cluster members assigned to altloc A vs altloc B by lDDT comparison. Loaded via `load_training_cluster_counts` (lowercases `protein`, renames `selection_part` → `selection`, drops the multi-part `selection` column, computes `frac_a_cluster = n_A / (n_A + n_B)` and `x_one_minus_fa = 1 − frac_a_cluster`).
- `results/Sampleworks_initial_structures.csv` — per-protein deposition date and `Pre`/`Post` cutoff flag. Loaded via `load_deposition_cutoffs` (dedupes `9B22`, fills `1VME, 4OLE, 5SOP, 6B8X` as `Pre`).
- `results/lddt_results_*.csv` — predicted occupancies; concatenated inline via `load_lddt(..., model_name_override=…)` (same five-file concat as nb 02), then narrowed to the **selected model group** via `filter_by_model_group(df, model_group.value)`.

**Common transforms:** `build_training_bias_df` parses the **second** float of each `occupancies` cell (numpy-style `"[0.375 0.625]"` → `0.625`) into a `b_occupancy` column, inner-joins lddt to the cluster frame on `(protein, selection)`, and left-joins the cutoff flag on `protein`. It builds for whichever models are in `bias_model_order` (the selected group); `OPTIMAL_GW` covers the variants.

**Models plotted:** the selected model group (Main → boltz2X-RAY, boltz2MD, protenix, rf3; Protenix variants → protenix, -mini, -tiny; RF3 variants → rf3, -chiral-off). The `make_training_bias_*` grids size to `len(model_order)` (4 / 3 / 2 columns); the Main 4-column layout is byte-identical to before. Non-main groups save with the `_protenix_alts` / `_rf3_alts` suffix. The Figure-36 explainer has no data and is group-independent.

#### Figure 29 — Training-Cluster vs Predicted State-B Fraction (Small Multiples)

**Output:** `figures/29_training_bias_small_multiples.svg`.

**What is plotted:** 4-column × 2-row scatter grid. Columns = main models, rows = condition (top: Baseline gw=0.0; bottom: Optimal guidance per `OPTIMAL_GW`). Shared axes (`xlim = ylim = [-0.02, 1.02]`). `y = 0.5` reference line drawn behind the points. Points are monochrome dark gray; alpha scales with `log(cluster_size + 1) / log(s_max + 1)` so larger-cluster cases dominate visually.

**Axes:** X = `frac_b_cluster`, Y = `b_occupancy`.

**Data included:** one point per `(protein, selection, model_name, condition)` with `s_cluster > 0`.

**What is excluded:** `s_cluster == 0` rows (these go to Figure 32 instead); all upstream filters.

#### Figure 34 — Calibration Overlay (Redesign A, 1×4)

**Output:** `figures/34_training_bias_calibration_overlay.svg` + legend at `figures/34_training_bias_calibration_overlay_legend.svg`.

**What is plotted:** 1-row × 4-column grid of panels, one per main model. Each panel pairs a scatter axis (left, width 4) with a marginal y-histogram strip (right, width 1) sharing y. Within each scatter, Baseline (gray, `#888888`) and Optimal guidance (model color from `MODEL_PALETTE`) points are overlaid. Reference lines `y = x` and `y = 0.5` are drawn behind the points. Per-condition binned conditional-mean curves use 10 equal-count quantile bins (or `n_bins=3, min_bin_count=2` for sparse subsets); each bin's marker sits at the within-bin mean of x. Slope of `b_occupancy` on `frac_b_cluster` and Pearson r are annotated for both conditions. Marginal axes overlay both conditions' `b_occupancy` distributions.

The main F34 cell calls the builder with `annotation_layout="split"`: the baseline slope/r block sits at top-left `(0.03, 0.97)` and the optimal-guidance slope/r block at bottom-right `(0.97, 0.03)`, so the two annotations never visually collide. The subset variants (Figures 30, 31) intentionally keep the default `"stacked"` layout (both annotations top-left).

Legend (in `make_training_bias_calibration_overlay_legend`): two-row vertical legend — "Baseline" with a single gray marker, "Optimal guidance" with four model-colored markers shown side-by-side via `HandlerTuple`.

#### Figure 35 — Calibration Small Multiples v2 (Redesign B, 2×4)

**Output:** `figures/35_training_bias_small_multiples_v2.svg`.

**What is plotted:** 2-row × 4-column grid mirroring Figure 29's geometry (top = Baseline, bottom = Optimal guidance, columns = main models). Each scatter panel is paired with a single-condition marginal y-histogram (per-panel, not per-condition like Figure 34). Every panel draws both reference lines, its own binned conditional-mean curve on top of a `mean ± std` shadow band, and annotates slope / Pearson r for the panel's own condition. Points monochrome (Figure 29 convention).

#### Figure 36 — Calibration Explainer (Empty Panel)

**Output:** `figures/36_training_bias_calibration_explainer.svg`.

**What is plotted:** A single 4×4 inch square-locked panel (`set_box_aspect(1)`) with the same axis style as the v2 small multiples and three lines: `y = x` (training prior, dashed), `y = 0.5` (guidance target, dashed), and a solid bisector at `y = 0.5*x + 0.25` (slope 0.5 through `(0.5, 0.5)`). Two italic black labels name the references — "training set prior" in the upper-right region rotated 45° to follow the `y = x` diagonal, and "perfect guidance" along the `y = 0.5` line left of center, unrotated. Two black double-headed arrows sit at the two ends of the solid bisector, perpendicular to it (purely visual, no text). No data; intended as a minimal explanatory companion to Figure 35 for downstream labeling.

#### Figures 30, 31 — Calibration Overlay on OOD Subsets

**Output filenames:**
- `figures/30_training_bias_calibration_overlay_low_cluster.svg`
- `figures/31_training_bias_calibration_overlay_post_cutoff.svg`

**What is plotted:** Both use the Figure 34 layout (`make_training_bias_calibration_overlay`). The only difference from Figure 34 is the subset filter applied to both conditions before plotting:

- **Figure 30 (low-cluster subset):** `max_cluster_size=5` — keeps only rows with `n_cluster_evaluated ≤ 5`. With only ~9 surviving points per panel, the binned mean uses `n_bins=3, min_bin_count=2`.
- **Figure 31 (post-cutoff subset):** `cutoff_status="Post"` — keeps only rows deposited after 2021-09-30.

#### Figure 32 — Empty-Training-Cluster Strip Plot

**Output:** `figures/32_training_bias_no_cluster_strip.svg`.

**What is plotted:** 1×4 row of panels (one per main model). Each panel pairs a half-violin (Gaussian KDE, clipped to `[0, 1]`, Scott bandwidth) on the left with a jittered scatter on the right; the X axis carries no quantitative meaning. Baseline (`#888888`) and Optimal-guidance (`MODEL_PALETTE[model]`) are overlaid in both halves (KDEs share the left region; jittered points share the right region). The KDE is skipped for any condition with fewer than 3 points (or zero variance); jittered points are always drawn. A `y = 0.5` reference line is drawn behind the data. Spines are stripped to a single left spine on the leftmost panel only; model names appear below each panel as the x-axis label. Each panel's Optimal-condition sample size is annotated.

This figure exists because its rows have `frac_b_cluster = N_B / (N_A + N_B)` undefined (both numerator and denominator zero), so they cannot appear on the X-axis of Figures 29, 30, 31, 34, or 35.

**Data included:** rows where `n_assigned_to_A == 0` AND `n_assigned_to_B == 0`, for the four main models, at Baseline and Optimal guidance.

#### Figure 33 — One-sided Training-Cluster Strip Plot

**Output:** `figures/33_training_bias_one_sided_strip.svg`.

**What is plotted:** 1×4 row of panels (one per main model, model name as panel title). Each panel holds two side-by-side jittered strips: left = "only-A" (`n_assigned_to_A > 0 AND n_assigned_to_B == 0`, `frac_b_cluster == 0`), right = "only-B" (`n_assigned_to_A == 0 AND n_assigned_to_B > 0`, `frac_b_cluster == 1`). X-tick labels are `A` and `B`; the X axis carries no quantitative meaning beyond the altloc label. Point placement uses `sns.stripplot` with `jitter=0.25` and two overlaid calls per panel — Baseline (`#888888`, alpha 0.7) drawn first, Optimal-guidance (`MODEL_PALETTE[model]`, alpha 0.7) drawn on top — so points alpha-blend across conditions. Each strip carries two horizontal mean ticks (gray for Baseline, model color for Optimal) at the strip's x-center, half-width ±0.25. A `y = 0.5` reference line is drawn behind the data. Per-strip `n = …` annotated just above the strip (Baseline and Optimal share row counts via the join, so a single `n` per strip is correct); empty strips render `n = 0` and skip the mean tick.

This figure is a sibling of Figure 32: Figure 32 covers `A == B == 0` (`frac_b_cluster` undefined), and Figure 33 covers the adjacent subset where the training cluster exists but is fully one-sided (`frac_b_cluster ∈ {0, 1}`). On Figure 34 these rows pile up at the `x ∈ {0, 1}` endpoints; the strip plot exposes their per-model spread instead.

**Data included:** rows where `(n_assigned_to_A > 0 AND n_assigned_to_B == 0)` OR `(n_assigned_to_A == 0 AND n_assigned_to_B > 0)`, for the four main models, at Baseline and Optimal guidance.

---

## Notebook 3: Recycles, Diffusion Steps & Ensemble Size (`04_recycles_steps_ensembles.py`)

### Source data

- `results/rscc_results_5mhx_recyc_diff.csv` — the **5MHX** recycles × diffusion-steps grid search, loaded by `load_rscc_recyc_diff` (same 50/50 occupancy filter + model-name aliasing as `load_rscc`, but keeps the `n_recycles` / `n_diffusion_steps` columns; NaN-`rscc` rows dropped with a print). Covers recycles `{1, 3, 10}` × steps `{50, 100, 200, 400}` for each of `boltz2MD`, `boltz2X-RAY`, `protenix`, `rf3`, 20 selections each, at each model's optimal guidance weight (`0.1` main / `0.01` rf3), ensemble size 8.
- `results/rscc_results_max_rmsd_boltz_protenix.csv` + `results/rscc_results_max_rmsd_rf3_fixed.csv` — the production **default** config, loaded via `load_rscc`, filtered to 5MHX at `OPTIMAL_GW` per model. Carries no recycles/steps columns, so it is pinned to `(DEFAULT_RECYCLES[model], 200)` where `DEFAULT_RECYCLES = {boltz2 → 3, protenix/rf3 → 10}`. (nb 04 stays on the four main models — no group dropdown.)
- `results/rscc_results_6b8x_ensemble_size.csv` — the **6B8X** ensemble-size sweep, loaded via `load_rscc` (it already keeps `ensemble_size`; same 50/50 filter + model-name aliasing). Covers ensemble sizes `{1, 2, 4, 8, 16}` for each of `boltz2MD`, `boltz2X-RAY`, `protenix`, `rf3`, 32 selections each, at each model's optimal guidance weight (`0.1` main / `0.01` rf3). Single-source: every model is already at its optimal gw, so there is **no default-point merge** (unlike the recycles/steps figure). NaN-free and error-free as delivered.

### Combine logic

The grid and default frames are concatenated, but a default `(model, recycles, steps)` cell is added **only when the grid does not already cover it** ("missing points come from the default"). The combine cell prints how many default cells it added. With the current full grid it adds **0** (every cell is covered), so the figure is **grid-only and self-consistent**; the default-merge is a fallback that activated only while the grid was incomplete. (Cross-check: at the default cell, grid RSCC matches the separate default CSV within ≈0.00–0.01 for protenix/rf3 and ≈0.05–0.07 for boltz2.)

### Figure — RSCC vs recycles ‖ RSCC vs diffusion steps (`recycles_steps/5mhx_rscc_vs_recycles_steps.svg`)

1×2 panel, mean RSCC ± SD **across the 20 selections** at each single config (so the SD band is over selections, not over configs). One line per model, colored by `MODEL_PALETTE`, labelled by `MODEL_DISPLAY`; markers + a translucent ±SD `fill_between`. X positions are evenly spaced categorical, labelled with the real values.

- **Left — RSCC vs recycles:** diffusion steps held at **200**; x = `{1, 3, 10}`.
- **Right — RSCC vs diffusion steps:** recycles held at each model's **default** (`boltz2 = 3`, `protenix`/`rf3 = 10`); x = `{50, 100, 200, 400}`.

**Data included:** 5MHX only; 50/50 occupancy; `selection_name == "all altloc"` excluded; main four models only (the grid has no protenix-mini/-tiny or rf3-chiral-off). Saved only when the `auto_save` switch is on.

### Figure — RSCC vs ensemble size (`ensemble_size/6b8x_rscc_vs_ensemble_size.svg`)

Single panel, mean RSCC ± SD **across the 32 selections** at each ensemble size (so the SD is over selections, not over configs). One line per model, colored by `MODEL_PALETTE`, labelled by `MODEL_DISPLAY`. X positions are `log2(size)` — for powers of two this is identical to even categorical spacing — labelled with the real sizes `{1, 2, 4, 8, 16}`. SD is drawn as **dodged capped error bars** (each model nudged ±off the shared tick by `0.06 * (i − (n−1)/2)`) rather than a `fill_between` band: here the per-selection spread is ~4× the signal range, so overlapping fills muddy the trend, whereas dodged bars keep the mean lines crisp and the models separable.

- **x = ensemble size**; held config = each model's optimal gw (no recycles/steps axis here).
- Mean RSCC rises then plateaus around size 4–8; the SD tightens as size grows (more members → more consistent across selections). Model ordering: Boltz2-MD ≈ Protenix > Boltz2-XRD > RF3.

**Data included:** 6B8X only; 50/50 occupancy; `selection_name == "all altloc"` excluded; main four models only (the CSV has no protenix-mini/-tiny or rf3-chiral-off). Saved only when the `auto_save` switch is on.

---

## Summary: What Is Universally Excluded

The following data **never appears in any plot**:

1. **Non-50/50 occupancy configurations** — only the equal-occupancy (A=0.5, B=0.5) conformer split is analyzed.
2. **Proteins or models missing from the source CSVs** — the analysis is limited to whatever was computed in the `results/` source CSVs.
3. **Baseline rows in delta plots** — deltas are defined as guided minus baseline, so the baseline itself (gw=0.0) does not appear as a data point in any delta/spiral figure (it is the implicit origin).
4. **(protein, model) pairs missing a baseline measurement** — the `make_delta_frame` merge is an inner join, so if there is no gw=0.0 row for a (protein, model), that pair is absent from all delta-based figures.

## Summary: Per-Figure Exclusion Table

_Saved figures only. The nb-02 printed summary tables are documented in their own section above ("Summary Tables — Baseline vs Optimal")._

| Figure | Excludes "all altloc" selection | Outlier/range clipping | Model subset |
|---|---|---|---|
| 04_09 Combined Baseline | Yes (LDDT side) | None | By dropdown |
| 12–17, 19 Per-Model Spirals | N/A (geometry+RSCC) | None | One model each |
| 20 Aggregated Spiral All | N/A | n < 3 groups dropped | By dropdown |
| 21–26, 28 Aggregated Spirals | N/A | n < 3 groups dropped | One model each |
| RSCC Scatter (per model) | Yes | None | One model each |
| Clustering Score Scatter (per model) | Yes (LDDT side) | None; warn callout if model has no clustering data | One model each |
| RSCC Scatter Aggregated | Yes | n < 3 main models dropped | Main only |
| Min-RMSD Scatter Aggregated | Yes | n < 3 main models dropped; NaN RMSD pre-dropped | Main only |
| Min-RMSD Scatter (per model) | Yes | NaN RMSD pre-dropped | One model each (data exists for 6 of 7) |
| RSCC vs LDDT (per model) | Yes (LDDT side) | None | One model each |
| Aggregated Clustering Score Scatter | Yes | n < 3 main models dropped | Main only |
| By-Class Combined Grid | Yes | None | By dropdown |
| 29 Training Bias Small Multiples | N/A (cluster-frame join) | s_cluster == 0 dropped | By dropdown |
| 30 Calibration Overlay Low Cluster | N/A | n_cluster_evaluated > 5 dropped | By dropdown |
| 31 Calibration Overlay Post Cutoff | N/A | cutoff_status != Post dropped | By dropdown |
| 32 No-Cluster Strip | N/A | A+B > 0 rows dropped | By dropdown |
| 33 One-sided Cluster Strip | N/A | A and B both > 0, or A+B == 0, rows dropped | By dropdown |
| 34 Calibration Overlay | N/A | None | By dropdown |
| 34L Calibration Overlay Legend | N/A (no data) | N/A | By dropdown |
| 35 Calibration Small Multiples v2 | N/A | None | By dropdown |
| 36 Calibration Explainer | N/A (no data) | None | N/A (no data) |
| 04 Recycles/Steps (5MHX) | Yes | None | 4 main models (5MHX only) |
| 04 Ensemble Size (6B8X) | Yes | None | 4 main models (6B8X only) |
