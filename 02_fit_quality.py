import marimo

__generated_with = "0.23.2"
app = marimo.App()


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell
def _():
    import matplotlib.pyplot as plt
    import pandas as pd
    from loaders import (
        add_simplified_metadata_columns,
        get_rscc_nan_proteins,
        load_classification,
        load_lddt,
        load_min_altloc_rmsd,
        load_rscc,
        load_rscc_max_rmsd,
        merge_classification_by_containment,
    )
    from plot_classified import make_by_class_combined_grid
    from plot_distributions import make_combined_baseline_distributions
    from plot_scatter import (
        make_aggregated_scatter,
        make_clustering_scatter,
        make_min_rmsd_scatter,
        make_rscc_scatter,
        make_scatter_comparison,
        rscc_residual_label_override,
    )
    from style import (
        CLASSIFICATION_CSV,
        MODEL_DISPLAY,
        MODELS_ALT,
        OPTIMAL_GW,
        save_fig,
        save_model_legend,
        save_per_model_figures,
    )
    from transforms import (
        filter_by_model_group,
        format_residue_range,
        model_slug,
    )

    return (
        CLASSIFICATION_CSV,
        MODELS_ALT,
        MODEL_DISPLAY,
        OPTIMAL_GW,
        add_simplified_metadata_columns,
        filter_by_model_group,
        format_residue_range,
        get_rscc_nan_proteins,
        load_classification,
        load_lddt,
        load_min_altloc_rmsd,
        load_rscc,
        load_rscc_max_rmsd,
        make_aggregated_scatter,
        make_by_class_combined_grid,
        make_clustering_scatter,
        make_combined_baseline_distributions,
        make_min_rmsd_scatter,
        make_rscc_scatter,
        make_scatter_comparison,
        merge_classification_by_containment,
        model_slug,
        pd,
        rscc_residual_label_override,
        save_fig,
        save_model_legend,
        save_per_model_figures,
    )


@app.cell
def _(OPTIMAL_GW, format_residue_range, pd):
    # ---- Shared helpers for the baseline-vs-optimal summary tables below ----
    # Each model's optimal guidance weight, as a frame to merge on.
    optimal_gw_df = pd.DataFrame(
        list(OPTIMAL_GW.items()), columns=["model_name", "optimal_gw"]
    )

    def baseline_vs_optimal(region, value_cols):
        """One row per (protein, selection, model) pairing each metric's baseline
        (guidance_weight == 0) value with its value at that model's optimal weight,
        named ``<col>_baseline`` / ``<col>_optimal`` for every col in ``value_cols``."""
        keys = ["protein", "selection", "model_name"]
        baseline = (
            region.query("guidance_weight == 0.0")[keys + value_cols]
            .rename(columns={c: f"{c}_baseline" for c in value_cols})
        )
        optimal = (
            region.merge(optimal_gw_df, on="model_name")
            .query("guidance_weight == optimal_gw")[keys + value_cols]
            .rename(columns={c: f"{c}_optimal" for c in value_cols})
        )
        return baseline.merge(optimal, on=keys)

    def add_label(df):
        """Return ``df`` with a ``label`` column like ``6B8X 181-183``."""
        df = df.copy()
        df["label"] = [
            f"{str(p).upper()} {format_residue_range(s)}"
            for p, s in zip(df["protein"], df["selection"])
        ]
        return df

    def per_model_counts(df, col="model"):
        """``"Model: n; ..."`` row-count summary, or ``"none"`` if empty."""
        return "; ".join(f"{m}: {n}" for m, n in df.groupby(col).size().items()) or "none"

    return add_label, baseline_vs_optimal, optimal_gw_df, per_model_counts


@app.cell
def _(mo):
    mo.md(r"""
    # Fit Quality Analysis

    How well do guided structure predictions fit experimental data? Two complementary
    metrics, each compared at **baseline** (no guidance) vs each model's **optimal
    guidance weight** (`OPTIMAL_GW`):

    1. **RSCC** — agreement between predicted and observed electron density.
    2. **Clustering score** — lDDT-based measure of how well the predicted ensemble
       clusters near the reference, weighted by recovered occupancy:
       `min(predicted_occ) * avg_silhouette_to_ref / min(map_occ)`.

    ## Figures (saved under `figures/`)

    | Figure | Output path |
    |--------|-------------|
    | Model color legend | `model_legend.svg` |
    | Baseline distributions (RSCC + clustering) | `04_09_baseline_distributions[_variant].svg` |
    | RSCC scatter, per model | `rscc_scatterplots/rscc_scatter_{model}.svg` |
    | Clustering score scatter, per model | `rscc_scatterplots/clustering_score_scatter_{model}.svg` |
    | RSCC vs clustering scatter, per model | `cross_metric_scatterplots/rscc_vs_lddt_{model}.svg` |
    | RSCC vs clustering scatter, per selection | `cross_metric_scatterplots/rscc_vs_lddt_per_selection_{model}.svg` |
    | Aggregated RSCC scatter | `rscc_scatterplots/rscc_scatter_aggregated_all_models.svg` |
    | Aggregated / per-model min-RMSD scatter | `rmsd_scatterplots/min_rmsd_scatter_*.svg` |
    | RSCC by conformational class (models × 4 grid, by dropdown) | `rscc_scatterplots/rscc_by_class_combined[_variant].svg` |
    | Aggregated clustering-score scatter | `rscc_scatterplots/clustering_score_scatter_aggregated_all_models.svg` |

    ## Summary tables (printed, not saved)

    Baseline-vs-optimal tables at the foot of the notebook, all at the
    per-(protein, selection, model) level: threshold-clearing fractions (overall and
    per model), single-residue enrichment, RSCC↔clustering trade-offs in both
    directions with a clustering-collapse breakdown, and RSCC vs min-RMSD-to-altloc.
    A **Merged Clustering Labels -- (+N) Expansion** table (just below the clustering
    scatter) lists every selection collapsed into each merged `(+N)` label.
    """)

    model_group = mo.ui.dropdown(
        options={"Main models": "main", "Protenix variants": "protenix", "RF3 variants": "rf3"},
        value="Main models",
        label="Model group",
    )
    smoothing = mo.ui.slider(
        start=0.05, stop=1.0, step=0.05,
        value=0.25,
        label="KDE smoothing",
    )
    auto_save = mo.ui.switch(label="Auto-save figures", value=False)
    exclude_nan_rscc = mo.ui.switch(label="Exclude NaN-RSCC proteins", value=False)
    cumulative = mo.ui.switch(label="Cumulative (fig 04_09)", value=False)
    cumulative_style = mo.ui.radio(
        options={"Smoothed (KDE)": "kde", "ECDF (step)": "ecdf"},
        value="Smoothed (KDE)",
        label="Cumulative style",
        inline=True,
    )
    mo.hstack(
        [model_group, smoothing, auto_save, exclude_nan_rscc, cumulative, cumulative_style],
        justify="start", gap=1,
    )
    return (
        auto_save,
        cumulative,
        cumulative_style,
        exclude_nan_rscc,
        model_group,
        smoothing,
    )


@app.cell
def _(exclude_nan_rscc, get_rscc_nan_proteins):
    nan_rscc_proteins = (
        get_rscc_nan_proteins()
        if exclude_nan_rscc.value else set()
    )
    return (nan_rscc_proteins,)


@app.cell
def _(add_simplified_metadata_columns, load_rscc_max_rmsd, nan_rscc_proteins):
    rscc = load_rscc_max_rmsd()
    if nan_rscc_proteins:
        rscc = rscc[~rscc["protein"].isin(nan_rscc_proteins)]
    # Collapse the per-altloc-region rows to a single RSCC per (protein, model, gw)
    # by averaging across all selections.
    rscc_all = (
        rscc.groupby(["protein", "model_name", "guidance_weight"], as_index=False)["rscc"]
        .mean()
    )
    rscc_all = add_simplified_metadata_columns(rscc_all)
    return rscc, rscc_all


@app.cell
def _(add_simplified_metadata_columns, load_lddt, nan_rscc_proteins, pd):
    # LDDT concatenated in-memory from the five per-model-group CSVs in results/.
    # protenix-mini/-tiny store model == "protenix" and rf3-chiral-off stores
    # model == "rf3", so those three are disambiguated via model_name_override at
    # load time; boltz_protenix and rf3_fixed carry distinguishable model columns.
    lddt_results = pd.concat(
        [
            load_lddt("results/lddt_results_boltz_protenix.csv"),
            load_lddt("results/lddt_results_rf3_fixed.csv"),
            load_lddt("results/lddt_results_protenix_mini.csv", model_name_override="protenix-mini"),
            load_lddt("results/lddt_results_protenix_tiny.csv", model_name_override="protenix-tiny"),
            load_lddt("results/lddt_results_rf3_chiral_off.csv", model_name_override="rf3-chiral-off"),
        ],
        ignore_index=True,
    )
    lddt_results_1 = add_simplified_metadata_columns(lddt_results)
    if nan_rscc_proteins:
        lddt_results_1 = lddt_results_1[~lddt_results_1["protein"].isin(nan_rscc_proteins)]
    return (lddt_results_1,)


@app.cell
def _(CLASSIFICATION_CSV, load_classification, nan_rscc_proteins):
    classification_df = load_classification(CLASSIFICATION_CSV)
    if nan_rscc_proteins:
        classification_df = classification_df[~classification_df["protein"].isin(nan_rscc_proteins)]
    return (classification_df,)


@app.cell
def _(load_min_altloc_rmsd, nan_rscc_proteins, pd):
    # Concatenate the four per-source min-RMSD CSVs into one tidy frame.
    # protenix-mini / protenix-tiny CSVs both store ``model == "protenix"``
    # so we override model_name at load time to disambiguate.
    _pieces = [
        load_min_altloc_rmsd("results/min_altloc_rmsd_results_boltz_protenix.csv"),
        load_min_altloc_rmsd(
            "results/min_altloc_rmsd_results_protenix_mini.csv",
            model_name_override="protenix-mini",
        ),
        load_min_altloc_rmsd(
            "results/min_altloc_rmsd_results_protenix_tiny.csv",
            model_name_override="protenix-tiny",
        ),
        load_min_altloc_rmsd("results/min_altloc_rmsd_results_rf3_fixed.csv"),
    ]
    min_rmsd_all = pd.concat(_pieces, ignore_index=True)
    if nan_rscc_proteins:
        min_rmsd_all = min_rmsd_all[~min_rmsd_all["protein"].isin(nan_rscc_proteins)]
    return (min_rmsd_all,)


@app.cell
def _(OPTIMAL_GW, add_label, min_rmsd_all, optimal_gw_df):
    # Per-(protein, selection, model) min RMSD at each model's optimal gw.
    _region = min_rmsd_all[min_rmsd_all["selection_name"] != "all altloc"]

    # Audit: which (protein, selection) pairs exist for each model in any gw
    # but lack a row at that model's optimal gw?
    _all_keys_per_model = {
        m: set(zip(g["protein"], g["selection"]))
        for m, g in _region.groupby("model_name")
    }

    min_rmsd_per_model_df = add_label(
        _region.merge(optimal_gw_df, on="model_name")
        .query("guidance_weight == optimal_gw")
        [["protein", "selection", "model_name", "guidance_weight",
          "min_rmsd_to_A", "min_rmsd_to_B"]]
        .reset_index(drop=True)
    )

    _kept_keys_per_model = {
        m: set(zip(g["protein"], g["selection"]))
        for m, g in min_rmsd_per_model_df.groupby("model_name")
    }
    for _m in sorted(_all_keys_per_model):
        _all = _all_keys_per_model[_m]
        _kept = _kept_keys_per_model.get(_m, set())
        _dropped = len(_all) - len(_all & _kept)
        if _dropped > 0:
            print(
                f"\u2139\ufe0f  min-RMSD per-model join: {_m} dropped {_dropped} of "
                f"{len(_all)} (protein, selection) pair(s) lacking a row at "
                f"optimal gw = {OPTIMAL_GW.get(_m)}."
            )

    # Cross-figure highlight groups for the min-RMSD scatter family,
    # nested per panel. Colors match the RSCC scatter convention in
    # plot_scatter.py:
    #   #F57C00 orange = far-from-ceiling family (2QCV 2-4, 4OLE 424-426,
    #                    6RP1 5),
    #   #1565C0 blue   = manually-curated biological highlights: 5SOP 131-133
    #                    plus the 6NI5 / 6NI6 altloc segments circled in blue
    #                    in the RSCC per-model figures (make_rscc_scatter's
    #                    highlight_proteins). Same color as the pinned blue
    #                    rings, so they read as one curated family across
    #                    figures.
    #   #D32F2F red    = decreased-RSCC family (curated 7Z0E altloc segments
    #                    flagged in the RSCC scatter; intentionally not the
    #                    full 7Z0E set to keep the panel readable).
    # The orange + blue groups are duplicated verbatim under both "autofit"
    # and "zoom" (no magic "both" key — explicit duplication keeps the
    # mapping readable). The 6NI5/6NI6 segments are listed on both panels;
    # the two 213-215 segments sit beyond the 3 A zoom range, so the panel
    # visibility filter drops them from the zoom panel automatically and they
    # appear on auto-fit only. The red 7Z0E set is intentionally zoom-only to
    # keep the global auto-fit panel readable.
    _ni56_segments = [
        "6NI5 151-153", "6NI5 201-202", "6NI5 204-206",
        "6NI5 209-210", "6NI5 213-215",
        "6NI6 151-153", "6NI6 201-202", "6NI6 204-206",
        "6NI6 209-210", "6NI6 213-215",
    ]
    MIN_RMSD_EXTRA_HIGHLIGHTS = {
        "autofit": {
            "#F57C00": ["2QCV 2-4", "4OLE 424-426", "6RP1 5"],
            "#1565C0": ["5SOP 131-133", *_ni56_segments],
        },
        "zoom": {
            "#F57C00": ["2QCV 2-4", "4OLE 424-426", "6RP1 5"],
            "#1565C0": ["5SOP 131-133", *_ni56_segments],
            "#D32F2F": [
                "7Z0E 11-13",
                "7Z0E 18-20",
                "7Z0E 54-56",
                "7Z0E 82-84",
                "7Z0E 122-124",
                "7Z0E 181-183",
                "7Z0E 206-208",
                "7Z0E 214-216",
            ],
        },
    }
    return MIN_RMSD_EXTRA_HIGHLIGHTS, min_rmsd_per_model_df


@app.cell
def _(auto_save, save_model_legend):
    if auto_save.value:
        save_model_legend()
    return


@app.cell
def _(
    auto_save,
    cumulative,
    cumulative_style,
    filter_by_model_group,
    lddt_results_1,
    make_combined_baseline_distributions,
    mo,
    model_group,
    rscc_all,
    save_fig,
    smoothing,
):
    _rscc_data, _suffix = filter_by_model_group(
        rscc_all, model_group.value, extra_query="guidance_strength == 'baseline'",
    )
    # Group-filter the LDDT panel too (was main-only) so both panels honor the
    # dropdown; the suffix is shared with the RSCC frame above.
    _lddt_data, _ = filter_by_model_group(
        lddt_results_1, model_group.value,
        extra_query="guidance_strength == 'baseline' and selection_name != 'all altloc'",
    )
    fig_04_09 = make_combined_baseline_distributions(
        _rscc_data, _lddt_data,
        smoothing=smoothing.value,
        cumulative=cumulative.value,
        cumulative_style=cumulative_style.value,
    )
    if auto_save.value:
        save_fig(fig_04_09, f"04_09_baseline_distributions{_suffix}.png")
    mo.vstack([
        mo.md("""## Figures 04 & 09 -- Baseline Distributions (RSCC + LDDT)
        All models at no-guidance baseline.
        `figures/04_09_baseline_distributions[_variant].svg`"""),
        fig_04_09,
    ])
    return


@app.cell
def _(MODEL_DISPLAY, mo):
    # Build dropdown from MODEL_DISPLAY (consistent labels across all scatter sections)
    _options = {
        v: k for k, v in MODEL_DISPLAY.items()
        if k in ("boltz2X-RAY", "boltz2MD", "protenix", "rf3",
                  "protenix-mini", "protenix-tiny", "rf3-chiral-off")
    }
    scatter_model = mo.ui.dropdown(options=_options, value="Boltz2-XRD", label="Model")
    save_all_scatter = mo.ui.switch(label="Save all models at once", value=False)
    mo.hstack([scatter_model, save_all_scatter], justify="start", gap=1)
    return save_all_scatter, scatter_model


@app.cell
def _(add_label, baseline_vs_optimal, rscc):
    # Per-(protein, selection, model) RSCC: baseline vs optimal guidance weight.
    # Feeds the RSCC scatters and most summary tables.
    _region = rscc[rscc["selection_name"] != "all altloc"]
    rscc_scatter_df = baseline_vs_optimal(_region, ["rscc"])
    rscc_scatter_df["delta_rscc"] = (
        rscc_scatter_df["rscc_optimal"] - rscc_scatter_df["rscc_baseline"]
    )
    rscc_scatter_df = add_label(rscc_scatter_df)
    return (rscc_scatter_df,)


@app.cell
def _(
    OPTIMAL_GW,
    auto_save,
    make_rscc_scatter,
    mo,
    model_slug,
    rscc_scatter_df,
    save_all_scatter,
    save_fig,
    save_per_model_figures,
    scatter_model,
):
    if save_all_scatter.value:
        save_per_model_figures(
            rscc_scatter_df, "model_name",
            lambda df, m: make_rscc_scatter(df, m, OPTIMAL_GW[m]),
            "rscc_scatterplots/rscc_scatter_{slug}.svg",
        )

    _cur = scatter_model.value
    _cur_df = rscc_scatter_df[rscc_scatter_df["model_name"] == _cur]
    _fig = make_rscc_scatter(_cur_df, _cur, OPTIMAL_GW[_cur])
    if auto_save.value:
        save_fig(_fig, f"rscc_scatterplots/rscc_scatter_{model_slug(_cur)}.svg")
    mo.vstack([
        mo.md("""## RSCC Scatter -- Baseline vs Optimal Guidance
        Per-model scatter comparing baseline RSCC with RSCC at optimal guidance weight.
        `figures/rscc_scatterplots/rscc_scatter_{model}.svg`"""),
        _fig,
    ])
    return


@app.cell
def _(
    MODEL_DISPLAY,
    OPTIMAL_GW,
    add_label,
    auto_save,
    cs_scatter_df,
    make_clustering_scatter,
    mo,
    model_slug,
    rscc_residual_label_override,
    rscc_scatter_df,
    save_all_scatter,
    save_fig,
    save_per_model_figures,
    scatter_model,
):
    # Clustering-score analog of the RSCC per-model scatter above. cs_scatter_df
    # is the clustering baseline-vs-optimal frame defined with the summary-table
    # frames below (marimo's DAG resolves the forward reference); it omits the
    # label column to avoid a collision in the rscc_cs_delta_df merge, so we add
    # the label locally here for annotations.
    _cs_labeled = add_label(cs_scatter_df)

    def _make_cluster_fig(_df, _m):
        # Pin the green/red/orange label sets to this model's RSCC residual panel
        # (fig 4 panel C) so the clustering panel labels the SAME selections in the
        # SAME colors. Blue curated highlights are shared inside the helper.
        _override = rscc_residual_label_override(
            rscc_scatter_df[rscc_scatter_df["model_name"] == _m]
        )
        return make_clustering_scatter(
            _df, _m, OPTIMAL_GW[_m], residual_label_override=_override
        )

    if save_all_scatter.value:
        save_per_model_figures(
            _cs_labeled, "model_name",
            _make_cluster_fig,
            "rscc_scatterplots/clustering_score_scatter_{slug}.svg",
        )

    _cur = scatter_model.value
    _cur_df = _cs_labeled[_cs_labeled["model_name"] == _cur]
    if _cur_df.empty:
        _content = mo.callout(
            f"No clustering-score data available for {MODEL_DISPLAY.get(_cur, _cur)}.",
            kind="warn",
        )
    else:
        _fig = _make_cluster_fig(_cur_df, _cur)
        if auto_save.value:
            save_fig(_fig, f"rscc_scatterplots/clustering_score_scatter_{model_slug(_cur)}.svg")
        _content = _fig
    mo.vstack([
        mo.md("""## Clustering Score Scatter -- Baseline vs Optimal Guidance
        Per-model scatter comparing baseline clustering score with the clustering score
        at optimal guidance weight. Green / red / orange labels are pinned to the same
        selections (in the same colors) as that model's RSCC residual panel (fig 4
        panel C); blue rings mark the shared curated highlights. Models without
        clustering data render a warning callout in place of the figure.
        `figures/rscc_scatterplots/clustering_score_scatter_{model}.svg`"""),
        _content,
    ])
    return


@app.cell
def _(
    MODEL_DISPLAY,
    OPTIMAL_GW,
    add_label,
    cs_scatter_df,
    make_clustering_scatter,
    mo,
    pd,
    rscc_residual_label_override,
    rscc_scatter_df,
):
    # Expand the (+N) merged clustering labels: which selections collapsed into each
    # combined "A, B (+N)" marker, per model. Built from the SAME grouping the
    # figures use -- make_clustering_scatter populates coincident_groups_out during
    # the draw -- so this table can never drift from the panels. A throwaway figure
    # is drawn and closed per model purely to harvest the groups.
    import matplotlib.pyplot as _plt

    _cs_labeled = add_label(cs_scatter_df)
    _rows = []
    for _m in [
        _name
        for _name in (
            "boltz2X-RAY", "boltz2MD", "protenix", "rf3",
            "protenix-mini", "protenix-tiny", "rf3-chiral-off",
        )
        if _name in set(_cs_labeled["model_name"])
    ]:
        _d = _cs_labeled[_cs_labeled["model_name"] == _m]
        if _d.empty:
            continue
        _override = rscc_residual_label_override(
            rscc_scatter_df[rscc_scatter_df["model_name"] == _m]
        )
        _groups: list = []
        _fig = make_clustering_scatter(
            _d, _m, OPTIMAL_GW[_m],
            residual_label_override=_override,
            coincident_groups_out=_groups,
        )
        _plt.close(_fig)
        for _gx, _gy, _members in sorted(_groups, key=lambda _g: -len(_g[2])):
            _rows.append({
                "model": MODEL_DISPLAY.get(_m, _m),
                "baseline_cs": round(float(_gx), 3),
                "improvement": round(float(_gy), 3),
                "n": len(_members),
                "selections": ", ".join(_members),
            })

    merged_label_groups_df = pd.DataFrame(
        _rows, columns=["model", "baseline_cs", "improvement", "n", "selections"]
    )
    mo.vstack([
        mo.md("""## Merged Clustering Labels -- (+N) Expansion
        Each row is a coincident group on a clustering panel whose member labels were
        merged into one combined `"A, B (+N)"` marker (selections landing on the same
        `coincident_tol` cell). Lists every selection in the group so the `(+N)` count
        is recoverable. Built from the same `coincident_groups_out` grouping the figures
        use, so it cannot drift from the panels. Not saved."""),
        merged_label_groups_df,
    ])
    return (merged_label_groups_df,)


@app.cell
def _(lddt_results_1, optimal_gw_df, rscc_all):
    # RSCC and clustering score per (protein, model) at each model's optimal gw,
    # clustering averaged across selections (excl. "all altloc").
    _rscc_at_opt = (
        rscc_all.merge(optimal_gw_df, on="model_name")
        .query("guidance_weight == optimal_gw")
        [["protein", "model_name", "rscc"]]
    )
    _lddt_at_opt = (
        lddt_results_1[lddt_results_1["selection_name"] != "all altloc"]
        .merge(optimal_gw_df, on="model_name")
        .query("guidance_weight == optimal_gw")
        .groupby(["protein", "model_name"], as_index=False)
        .agg(clustering_score=("clustering_score", "mean"))
    )
    cross_metric_df = _rscc_at_opt.merge(_lddt_at_opt, on=["protein", "model_name"])
    return (cross_metric_df,)


@app.cell
def _(
    auto_save,
    cross_metric_df,
    make_scatter_comparison,
    mo,
    model_slug,
    save_all_scatter,
    save_fig,
    save_per_model_figures,
    scatter_model,
):
    _scatter_kwargs = dict(
        x_col="rscc", y_col="clustering_score",
        x_label="RSCC",
        y_label="Clustering Score",
        draw_regression=True, show_correlation=True, label_all=True, equal_aspect=False,
    )

    if save_all_scatter.value:
        save_per_model_figures(
            cross_metric_df, "model_name",
            make_scatter_comparison,
            "cross_metric_scatterplots/rscc_vs_lddt_{slug}.svg",
            make_kwargs=_scatter_kwargs,
        )

    _cur = scatter_model.value
    _cur_df = cross_metric_df[cross_metric_df["model_name"] == _cur]
    _fig = make_scatter_comparison(_cur_df, _cur, **_scatter_kwargs)
    if auto_save.value:
        save_fig(_fig, f"cross_metric_scatterplots/rscc_vs_lddt_{model_slug(_cur)}.svg")
    mo.vstack([
        mo.md("""## RSCC vs LDDT Scatter -- Cross-Metric Comparison
        Per-model scatter comparing RSCC with LDDT clustering score at optimal guidance,
        clustering averaged across selections (one point per protein). Pearson r annotated.
        `figures/cross_metric_scatterplots/rscc_vs_lddt_{model}.svg`"""),
        _fig,
    ])
    return


@app.cell
def _(add_label, lddt_results_1, optimal_gw_df, rscc):
    # RSCC and clustering score per (protein, selection, model) at each model's optimal
    # gw -- the same cross-metric comparison but NOT averaged across selections, so each
    # point is one altloc region rather than one protein.
    _rscc_at_opt = (
        rscc[rscc["selection_name"] != "all altloc"]
        .merge(optimal_gw_df, on="model_name")
        .query("guidance_weight == optimal_gw")
        [["protein", "selection", "model_name", "rscc"]]
    )
    _lddt_at_opt = (
        lddt_results_1[lddt_results_1["selection_name"] != "all altloc"]
        .merge(optimal_gw_df, on="model_name")
        .query("guidance_weight == optimal_gw")
        [["protein", "selection", "model_name", "clustering_score"]]
    )
    cross_metric_per_sel_df = add_label(
        _rscc_at_opt.merge(_lddt_at_opt, on=["protein", "selection", "model_name"])
    )
    return (cross_metric_per_sel_df,)


@app.cell
def _(
    auto_save,
    cross_metric_per_sel_df,
    make_scatter_comparison,
    mo,
    model_slug,
    save_all_scatter,
    save_fig,
    save_per_model_figures,
    scatter_model,
):
    _scatter_kwargs = dict(
        x_col="rscc", y_col="clustering_score",
        x_label="RSCC",
        y_label="Clustering Score",
        draw_regression=True, show_correlation=True, equal_aspect=False,
        label_col="label", highlight_proteins=["6B8X 181-183", "5MHX 162-164"],
    )

    if save_all_scatter.value:
        save_per_model_figures(
            cross_metric_per_sel_df, "model_name",
            make_scatter_comparison,
            "cross_metric_scatterplots/rscc_vs_lddt_per_selection_{slug}.svg",
            make_kwargs=_scatter_kwargs,
        )

    _cur = scatter_model.value
    _cur_df = cross_metric_per_sel_df[cross_metric_per_sel_df["model_name"] == _cur]
    _fig = make_scatter_comparison(_cur_df, _cur, **_scatter_kwargs)
    if auto_save.value:
        save_fig(_fig, f"cross_metric_scatterplots/rscc_vs_lddt_per_selection_{model_slug(_cur)}.svg")
    mo.vstack([
        mo.md("""## RSCC vs LDDT Scatter -- Per Selection (Cross-Metric)
        Same RSCC vs clustering comparison as above but with **one point per (protein,
        selection)** altloc region instead of one per protein -- clustering is not averaged
        across selections. Pearson r annotated; `6B8X 181-183` / `5MHX 162-164` circled blue.
        `figures/cross_metric_scatterplots/rscc_vs_lddt_per_selection_{model}.svg`"""),
        _fig,
    ])
    return


@app.cell
def _(auto_save, make_aggregated_scatter, mo, rscc_scatter_df, save_fig):
    _rscc_agg_df = (
        rscc_scatter_df.query("model_name not in @MODELS_ALT")
        .groupby(["protein", "selection", "label"], as_index=False)
        .agg(
            x_med=("rscc_baseline", "median"),
            y_med=("rscc_optimal", "median"),
            n=("model_name", "count"),
        )
    )
    _rscc_agg_df = _rscc_agg_df[_rscc_agg_df["n"] >= 3]

    _fig = make_aggregated_scatter(
        _rscc_agg_df,
        x_col="x_med", y_col="y_med",
        x_label="Baseline RSCC (median across models)",
        y_label="Median Improvement (Guided - Baseline)",
        title="RSCC Improvement Over Baseline",
        label_col="label",
        residual=True,
        show_marginal_kde=True,
        show_errors=False,
        point_color="#377EB8",
        label_top_n=3, label_far_from_ceiling_n=5,
        highlight_proteins=["6B8X 181-183", "5MHX 162-164", "5SOP 131-133"],
        xlim=(0.0, 1.0),
    )
    if auto_save.value:
        save_fig(_fig, "rscc_scatterplots/rscc_scatter_aggregated_all_models.png")
    mo.vstack([
        mo.md("""## Aggregated RSCC Scatter -- Collapsed Across Main Models
        Each point is one (protein, selection); baseline and guided RSCC collapsed across main models using the median (per-model spread shown in the per-model figures above).
        `figures/rscc_scatterplots/rscc_scatter_aggregated_all_models.svg`"""),
        _fig,
    ])
    return


@app.cell
def _(
    MIN_RMSD_EXTRA_HIGHLIGHTS,
    auto_save,
    make_min_rmsd_scatter,
    min_rmsd_per_model_df,
    mo,
    save_fig,
):
    _min_rmsd_agg_df = (
        min_rmsd_per_model_df.query("model_name not in @MODELS_ALT")
        .groupby(["protein", "selection", "label"], as_index=False)
        .agg(
            min_rmsd_to_A=("min_rmsd_to_A", "median"),
            min_rmsd_to_B=("min_rmsd_to_B", "median"),
            n=("model_name", "count"),
        )
    )
    _min_rmsd_agg_df = _min_rmsd_agg_df[_min_rmsd_agg_df["n"] >= 3]

    # The 6NI5/6NI6 altloc segments stay on the per-model RMSD figures but are
    # dropped from this aggregated scatter (per user request). Filter them out
    # of every panel/color group of the shared highlight dict so the per-model
    # cell, which still passes MIN_RMSD_EXTRA_HIGHLIGHTS unchanged, is untouched.
    _agg_highlights = {
        _panel: {
            _color: [
                _lbl for _lbl in _labels
                if not _lbl.startswith(("6NI5", "6NI6"))
            ]
            for _color, _labels in _groups.items()
        }
        for _panel, _groups in MIN_RMSD_EXTRA_HIGHLIGHTS.items()
    }

    _fig_min_rmsd_agg = make_min_rmsd_scatter(
        _min_rmsd_agg_df,
        title="Minimum RMSD to Each Altloc",
        extra_highlight_groups=_agg_highlights,
    )
    if auto_save.value:
        save_fig(
            _fig_min_rmsd_agg,
            "rmsd_scatterplots/min_rmsd_scatter_aggregated_all_models.svg",
        )
    mo.vstack([
        mo.md("""## Aggregated Min-RMSD Scatter -- Min RMSD to Altloc A vs B
        Each point is one (protein, selection); X / Y are the median across main models of the minimum RMSD across the 8-member ensemble at each model's optimal guidance weight. Pinned blue rings on `6B8X 181-183` and `5MHX 162-164` appear on both panels (matching the aggregated RSCC scatter). Extra RSCC-matched highlight rings split per panel: orange (`#F57C00`) `2QCV 2-4`, `4OLE 424-426` and `6RP1 5` and blue (`#1565C0`) `5SOP 131-133` appear on **both** the auto-fit and zoom panels; the red (`#D32F2F`) curated `7Z0E` subset (`11-13`, `18-20`, `54-56`, `82-84`, `122-124`, `181-183`, `206-208`, `214-216`) appears on the **zoom panel only**, keeping the global auto-fit panel readable. The 6NI5 / 6NI6 altloc segments circled blue on the per-model RMSD figures are intentionally omitted from this aggregated scatter.
        `figures/rmsd_scatterplots/min_rmsd_scatter_aggregated_all_models.svg`"""),
        _fig_min_rmsd_agg,
    ])
    return


@app.cell
def _(
    MIN_RMSD_EXTRA_HIGHLIGHTS,
    MODEL_DISPLAY,
    auto_save,
    make_min_rmsd_scatter,
    min_rmsd_per_model_df,
    mo,
    model_slug,
    save_all_scatter,
    save_fig,
    save_per_model_figures,
    scatter_model,
):
    def _make_for_model(_df, _m):
        return make_min_rmsd_scatter(
            _df,
            title=f"Minimum RMSD to Each Altloc: {MODEL_DISPLAY.get(_m, _m)}",
            extra_highlight_groups=MIN_RMSD_EXTRA_HIGHLIGHTS,
        )

    if save_all_scatter.value:
        save_per_model_figures(
            min_rmsd_per_model_df, "model_name",
            _make_for_model,
            "rmsd_scatterplots/min_rmsd_scatter_{slug}.svg",
        )

    _cur = scatter_model.value
    _cur_df = min_rmsd_per_model_df[min_rmsd_per_model_df["model_name"] == _cur]
    if _cur_df.empty:
        _content = mo.callout(
            f"No min-RMSD data available for {MODEL_DISPLAY.get(_cur, _cur)}.",
            kind="warn",
        )
    else:
        _fig_min_rmsd_per = _make_for_model(_cur_df, _cur)
        if auto_save.value:
            save_fig(
                _fig_min_rmsd_per,
                f"rmsd_scatterplots/min_rmsd_scatter_{model_slug(_cur)}.svg",
            )
        _content = _fig_min_rmsd_per
    mo.vstack([
        mo.md("""## Per-Model Min-RMSD Scatter -- Min RMSD to Altloc A vs B
        Per-model scatter at each model's optimal guidance weight (no baseline overlay). Highlight split: pinned blues plus the orange + blue extra groups (blue here includes the 6NI5 / 6NI6 altloc segments matching the RSCC per-model figures -- unlike the aggregated figure, which omits them) appear on both panels, while the red curated `7Z0E` subset appears on the zoom panel only. Per-panel visibility is filtered to in-range points, so segments outside a panel's limits (e.g. 6NI5/6NI6 `213-215` on the zoom panel) simply don't render there. Models without a min-RMSD CSV (e.g. `rf3-chiral-off`) render a warning callout in place of the figure.
        `figures/rmsd_scatterplots/min_rmsd_scatter_{model}.svg`"""),
        _content,
    ])
    return


@app.cell
def _(
    auto_save,
    classification_df,
    filter_by_model_group,
    load_rscc,
    make_by_class_combined_grid,
    merge_classification_by_containment,
    mo,
    model_group,
    nan_rscc_proteins,
    optimal_gw_df,
    pd,
    save_fig,
):
    # Per-altloc-region RSCC from the BASE (non-max) per-model CSVs in results/,
    # concatenated in-memory. The by-class grid keeps the base frame (not the
    # max-RMSD one) because classification needs the base's full per-region
    # selection coverage. All four conformational classes — including the lone
    # domain_shift region (2e7z) — now flow through a single containment merge;
    # the retired rscc_results_domain_shift.csv is no longer needed. Variant
    # files are disambiguated via model_name_override, mirroring the max-RMSD
    # sources (RSCC_MAX_RMSD_SOURCES).
    _base_rscc = pd.concat(
        [
            load_rscc("results/rscc_results_boltz_protenix.csv"),
            load_rscc("results/rscc_results_rf3_fixed.csv"),
            load_rscc("results/rscc_results_protenix_mini.csv", model_name_override="protenix-mini"),
            load_rscc("results/rscc_results_protenix_tiny.csv", model_name_override="protenix-tiny"),
            load_rscc("results/rscc_results_rf3_chiral_off.csv", model_name_override="rf3-chiral-off"),
        ],
        ignore_index=True,
    )
    if nan_rscc_proteins:
        _base_rscc = _base_rscc[~_base_rscc["protein"].isin(nan_rscc_proteins)]
    _region = _base_rscc[_base_rscc["selection_name"] != "all altloc"].copy()
    _region_cls = merge_classification_by_containment(_region, classification_df)

    _baseline = (
        _region_cls.query("guidance_weight == 0.0")
        [["protein", "selection", "model_name", "rscc", "classification"]]
        .rename(columns={"rscc": "rscc_baseline"})
    )
    _optimal = (
        _region_cls.merge(optimal_gw_df, on="model_name")
        .query("guidance_weight == optimal_gw")
        [["protein", "selection", "model_name", "rscc"]]
        .rename(columns={"rscc": "rscc_optimal"})
    )
    _by_class_df = _baseline.merge(_optimal, on=["protein", "selection", "model_name"])

    # The model-group dropdown picks which family fills the grid rows: Main →
    # the four main models; Protenix / RF3 variants → that family's members and
    # colors. Non-main saves get the _protenix_alts / _rf3_alts suffix.
    _merged, _suffix = filter_by_model_group(_by_class_df, model_group.value)
    _canonical = [
        "boltz2X-RAY", "boltz2MD", "protenix", "rf3",
        "protenix-mini", "protenix-tiny", "rf3-chiral-off",
    ]
    _present = set(_merged["model_name"].unique())
    _model_order = [m for m in _canonical if m in _present]

    _fig_combined = make_by_class_combined_grid(_merged, model_order=_model_order)
    if auto_save.value:
        save_fig(_fig_combined, f"rscc_scatterplots/rscc_by_class_combined{_suffix}.svg")

    mo.vstack([
        mo.md("""## RSCC by Conformational Type — Combined Grid (models × 4 classes)
        Single combined figure: rows = the selected model group's members, columns =
        conformational-change classes. Per-panel gray shadow shows that row's model
        pooled cloud. All four classes — including the lone `domain_shift` region
        (`2e7z`) — come from one in-memory concat of the per-altloc-region base
        `results/rscc_results_*.csv` files, classified by containment against
        `results/classify_altloc_regions.csv`. Honors the model-group dropdown;
        non-main groups save with the `_protenix_alts` / `_rf3_alts` suffix.
        `figures/rscc_scatterplots/rscc_by_class_combined[_variant].svg`"""),
        _fig_combined,
    ])
    return


@app.cell
def _(
    add_label,
    auto_save,
    lddt_results_1,
    make_aggregated_scatter,
    mo,
    optimal_gw_df,
    save_fig,
):
    # Per-(protein, selection, model) baseline + optimal clustering score (main models),
    # collapsed to one (protein, selection) median.
    _region = lddt_results_1[lddt_results_1["selection_name"] != "all altloc"]
    _baseline = (
        _region.query("guidance_weight == 0.0 and model_name not in @MODELS_ALT")
        [["protein", "selection", "model_name", "clustering_score"]]
        .rename(columns={"clustering_score": "lddt_baseline"})
    )
    _optimal = (
        _region.merge(optimal_gw_df, on="model_name")
        .query("guidance_weight == optimal_gw and model_name not in @MODELS_ALT")
        [["protein", "selection", "model_name", "clustering_score"]]
        .rename(columns={"clustering_score": "lddt_optimal"})
    )
    _agg = (
        _baseline.merge(_optimal, on=["protein", "selection", "model_name"])
        .groupby(["protein", "selection"], as_index=False)
        .agg(x_med=("lddt_baseline", "median"), y_med=("lddt_optimal", "median"),
             n=("model_name", "count"))
    )
    _agg = add_label(_agg[_agg["n"] >= 3])

    _fig_lddt = make_aggregated_scatter(
        _agg,
        x_col="x_med", y_col="y_med",
        x_label="Baseline Clustering Score (median across models)",
        y_label="Median Improvement (Guided - Baseline)",
        title="Clustering Score Improvement Over Baseline",
        label_col="label",
        residual=True,
        show_marginal_kde=True, show_errors=False, point_color="#377EB8",
        label_top_n=3, label_bottom_n=5, label_far_from_ceiling_n=3,
        highlight_proteins=["6B8X 181-183", "5MHX 162-164", "5SOP 131-133"],
        xlim=(0.0, 1.0),
    )
    if auto_save.value:
        save_fig(_fig_lddt, "rscc_scatterplots/clustering_score_scatter_aggregated_all_models.png")

    mo.vstack([
        mo.md("""## Aggregated Clustering Score Scatter
        Each point is one (protein, selection), median across main models of baseline vs
        guided clustering score (per-model spread is in the per-model figures)."""),
        _fig_lddt,
    ])
    return


@app.cell
def _(baseline_vs_optimal, lddt_results_1):
    # Per-(protein, selection, model) clustering score: baseline vs optimal gw, plus the
    # silhouette and min predicted occupancy that define it (used by the collapse
    # analysis). Parallel to rscc_scatter_df; feeds the summary tables below.
    _region = lddt_results_1[lddt_results_1["selection_name"] != "all altloc"]
    cs_scatter_df = baseline_vs_optimal(
        _region, ["clustering_score", "avg_silhouette_to_ref", "min_predicted_occupancy"]
    )
    cs_scatter_df["delta_clustering_score"] = (
        cs_scatter_df["clustering_score_optimal"] - cs_scatter_df["clustering_score_baseline"]
    )
    return (cs_scatter_df,)


@app.cell
def _(MODELS_ALT, cs_scatter_df, mo, pd, rscc_scatter_df):
    # Fraction of predictions at or above each quality threshold, baseline vs optimal.
    # Unit = (protein, selection, model), main models only; baseline and optimal share
    # the identical point set. Thresholds inclusive (>=).
    _rscc = rscc_scatter_df[~rscc_scatter_df["model_name"].isin(MODELS_ALT)].dropna(
        subset=["rscc_baseline", "rscc_optimal"]
    )
    _cs = cs_scatter_df[~cs_scatter_df["model_name"].isin(MODELS_ALT)].dropna(
        subset=["clustering_score_baseline", "clustering_score_optimal"]
    )
    _summary = pd.DataFrame([
        {
            "metric": "RSCC >= 0.8",
            "baseline_fraction": (_rscc["rscc_baseline"] >= 0.8).mean(),
            "optimal_fraction": (_rscc["rscc_optimal"] >= 0.8).mean(),
            "n": len(_rscc),
        },
        {
            "metric": "RSCC >= 0.9",
            "baseline_fraction": (_rscc["rscc_baseline"] >= 0.9).mean(),
            "optimal_fraction": (_rscc["rscc_optimal"] >= 0.9).mean(),
            "n": len(_rscc),
        },
        {
            "metric": "Clustering score >= 0.7",
            "baseline_fraction": (_cs["clustering_score_baseline"] >= 0.7).mean(),
            "optimal_fraction": (_cs["clustering_score_optimal"] >= 0.7).mean(),
            "n": len(_cs),
        },
        {
            "metric": "Clustering score >= 0.5",
            "baseline_fraction": (_cs["clustering_score_baseline"] >= 0.5).mean(),
            "optimal_fraction": (_cs["clustering_score_optimal"] >= 0.5).mean(),
            "n": len(_cs),
        },
    ])
    mo.vstack([
        mo.md("""## Threshold-Clearing Fractions -- Baseline vs Optimal
        Fraction of per-(protein, selection, model) predictions (main models, "all altloc"
        excluded) at or above each threshold, over the identical baseline/optimal point set."""),
        _summary,
    ])
    return


@app.cell
def _(MODEL_DISPLAY, cs_scatter_df, mo, rscc_scatter_df):
    # Per-model version of the threshold-clearing fractions (same metrics/thresholds as
    # the overall table, but per model and including the variant models).
    _rscc = rscc_scatter_df.dropna(subset=["rscc_baseline", "rscc_optimal"])
    _rscc_pm = _rscc.groupby("model_name").agg(
        rscc_ge08_baseline=("rscc_baseline", lambda s: (s >= 0.8).mean()),
        rscc_ge08_optimal=("rscc_optimal", lambda s: (s >= 0.8).mean()),
        rscc_ge09_baseline=("rscc_baseline", lambda s: (s >= 0.9).mean()),
        rscc_ge09_optimal=("rscc_optimal", lambda s: (s >= 0.9).mean()),
        rscc_n=("rscc_baseline", "size"),
    )
    _cs = cs_scatter_df.dropna(subset=["clustering_score_baseline", "clustering_score_optimal"])
    _cs_pm = _cs.groupby("model_name").agg(
        cs_ge07_baseline=("clustering_score_baseline", lambda s: (s >= 0.7).mean()),
        cs_ge07_optimal=("clustering_score_optimal", lambda s: (s >= 0.7).mean()),
        cs_ge05_baseline=("clustering_score_baseline", lambda s: (s >= 0.5).mean()),
        cs_ge05_optimal=("clustering_score_optimal", lambda s: (s >= 0.5).mean()),
        cs_n=("clustering_score_baseline", "size"),
    )
    _per_model = _rscc_pm.join(_cs_pm, how="outer").sort_index().reset_index()
    _per_model.insert(
        0, "model", _per_model.pop("model_name").map(lambda m: MODEL_DISPLAY.get(m, m))
    )

    mo.vstack([
        mo.md("""## Threshold-Clearing Fractions by Model
        Per-model version of the table above. `rscc_ge08_*` / `rscc_ge09_*` = fraction
        with RSCC >= 0.8 / >= 0.9; `cs_ge07_*` / `cs_ge05_*` = fraction with clustering
        score >= 0.7 / >= 0.5; `*_n` = predictions counted (RSCC and clustering have
        different selection coverage). Variant models without clustering data show NaN in
        the `cs_*` columns."""),
        _per_model,
    ])
    return


@app.cell
def _(MODELS_ALT, format_residue_range, mo, pd, rscc_scatter_df):
    # Are single-residue selections (range resi N-N) overrepresented among baseline
    # RSCC >= 0.8 predictions vs the selection universe? Unit = (protein, selection,
    # model) baseline prediction, main models only.
    _base = rscc_scatter_df[~rscc_scatter_df["model_name"].isin(MODELS_ALT)].dropna(
        subset=["rscc_baseline"]
    ).copy()
    _base["is_single"] = _base["selection"].map(lambda s: "-" not in format_residue_range(s))
    _pass = _base[_base["rscc_baseline"] >= 0.8]

    _summary = pd.DataFrame([
        {"set": "Baseline RSCC >= 0.8", "n_predictions": len(_pass),
         "n_single_residue": int(_pass["is_single"].sum()),
         "frac_single_residue": _pass["is_single"].mean()},
        {"set": "All baseline predictions", "n_predictions": len(_base),
         "n_single_residue": int(_base["is_single"].sum()),
         "frac_single_residue": _base["is_single"].mean()},
    ])
    _ratio = _summary.loc[0, "frac_single_residue"] / _summary.loc[1, "frac_single_residue"]

    mo.vstack([
        mo.md(f"""## Single-Residue Selections in the Baseline RSCC >= 0.8 Set
        Single-residue predictions are **{_summary.loc[0, 'frac_single_residue']:.1%}** of the
        baseline RSCC >= 0.8 set vs **{_summary.loc[1, 'frac_single_residue']:.1%}** of all
        baseline predictions -- an overrepresentation factor of **{_ratio:.2f}x**."""),
        _summary,
    ])
    return


@app.cell
def _(MODEL_DISPLAY, cs_scatter_df, mo, per_model_counts, rscc_scatter_df):
    # RSCC and clustering deltas joined per (protein, selection, model). The inner join
    # keeps only selections with BOTH metrics (drops protenix-mini/-tiny, which have no
    # clustering data). rscc_cs_delta_df also feeds the reverse-direction table below.
    rscc_cs_delta_df = rscc_scatter_df.merge(
        cs_scatter_df, on=["protein", "selection", "model_name"]
    )
    rscc_cs_delta_df["model"] = rscc_cs_delta_df["model_name"].map(
        lambda m: MODEL_DISPLAY.get(m, m)
    )

    tradeoff_df = (
        rscc_cs_delta_df[
            (rscc_cs_delta_df["delta_rscc"] > 0)
            & (rscc_cs_delta_df["delta_clustering_score"] < 0)
        ]
        .sort_values(["model", "delta_clustering_score"])
        .reset_index(drop=True)
    )
    _display = tradeoff_df[
        ["model", "label", "rscc_baseline", "rscc_optimal", "delta_rscc",
         "clustering_score_baseline", "clustering_score_optimal", "delta_clustering_score"]
    ]

    mo.vstack([
        mo.md(f"""## RSCC Up / Clustering Down -- Trade-off Selections (per model)
        Optimal guidance **raises RSCC but lowers clustering score** (strict). RSCC from
        the max-RMSD CSVs, clustering from base lDDT, inner-joined. Sorted by model then
        clustering drop. Per model: {per_model_counts(tradeoff_df)}."""),
        _display,
    ])
    return rscc_cs_delta_df, tradeoff_df


@app.cell
def _(mo, pd, tradeoff_df):
    # Within the trade-off set: how many collapse clustering fully to 0 (min predicted
    # occupancy == 0 -> all ensemble weight on one conformation) vs a partial decrease?
    # In the collapsed cases, is avg_silhouette_to_ref still high (a tight single
    # conformation near the reference, abandoning the other altloc)?
    _sil = "avg_silhouette_to_ref_optimal"
    _df = tradeoff_df.assign(collapsed=tradeoff_df["min_predicted_occupancy_optimal"] == 0.0)
    _coll = _df[_df["collapsed"]]

    _counts = _df.groupby("model").agg(
        n_tradeoff=("collapsed", "size"), n_collapsed=("collapsed", "sum"),
    )
    _counts["n_partial"] = _counts["n_tradeoff"] - _counts["n_collapsed"]
    _counts["frac_collapsed"] = _counts["n_collapsed"] / _counts["n_tradeoff"]
    _sil_stats = _coll.groupby("model").agg(
        collapsed_sil_mean=(_sil, "mean"),
        collapsed_sil_median=(_sil, "median"),
        collapsed_frac_sil_ge_05=(_sil, lambda s: (s >= 0.5).mean()),
    )
    _per_model = _counts.join(_sil_stats, how="left").reset_index()
    _overall = pd.DataFrame([{
        "model": "All models",
        "n_tradeoff": len(_df), "n_collapsed": int(_df["collapsed"].sum()),
        "n_partial": int((~_df["collapsed"]).sum()), "frac_collapsed": _df["collapsed"].mean(),
        "collapsed_sil_mean": _coll[_sil].mean(),
        "collapsed_sil_median": _coll[_sil].median(),
        "collapsed_frac_sil_ge_05": (_coll[_sil] >= 0.5).mean(),
    }])
    _summary = pd.concat([_per_model, _overall], ignore_index=True)

    mo.vstack([
        mo.md(f"""## Clustering Collapse vs Partial Decrease (within the trade-off set)
        Of {len(_df)} RSCC-up / clustering-down selections,
        **{int(_df['collapsed'].sum())} ({_df['collapsed'].mean():.0%})** collapse clustering
        to 0 (zero occupancy on one altloc -> a single conformation); the rest partially
        decrease. In the collapsed cases, mean `avg_silhouette_to_ref` =
        **{_coll[_sil].mean():.2f}** (median **{_coll[_sil].median():.2f}**,
        **{(_coll[_sil] >= 0.5).mean():.0%}** at >= 0.5): a high silhouette with a zero
        clustering score means the model collapsed onto a single tight conformation near
        the reference, missing the other altloc (0.5 cutoff illustrative)."""),
        _summary,
    ])
    return


@app.cell
def _(mo, per_model_counts, rscc_cs_delta_df):
    # Opposite trade-off: guidance lowers RSCC but raises clustering score (strict) --
    # the mirror of the table above. Sorted by model then RSCC drop.
    _reverse = (
        rscc_cs_delta_df[
            (rscc_cs_delta_df["delta_rscc"] < 0)
            & (rscc_cs_delta_df["delta_clustering_score"] > 0)
        ]
        .sort_values(["model", "delta_rscc"])
        .reset_index(drop=True)
    )
    _display = _reverse[
        ["model", "label", "rscc_baseline", "rscc_optimal", "delta_rscc",
         "clustering_score_baseline", "clustering_score_optimal", "delta_clustering_score"]
    ]
    mo.vstack([
        mo.md(f"""## RSCC Down / Clustering Up -- Reverse Trade-off Selections (per model)
        Optimal guidance **lowers RSCC but raises clustering score** (strict) -- the
        mirror of the table above. Per model: {per_model_counts(_reverse)}."""),
        _display,
    ])
    return


@app.cell
def _(
    MODEL_DISPLAY,
    baseline_vs_optimal,
    min_rmsd_all,
    mo,
    per_model_counts,
    rscc_scatter_df,
):
    # RSCC delta vs ensemble min-RMSD to each altloc, per (protein, selection, model).
    # min-RMSD coverage is partial (no rf3-chiral-off yet); the inner join keeps only
    # models with both RSCC and min-RMSD data.
    _region = min_rmsd_all[min_rmsd_all["selection_name"] != "all altloc"]
    _rmsd = baseline_vs_optimal(_region, ["min_rmsd_to_A", "min_rmsd_to_B"])
    rscc_minrmsd_df = rscc_scatter_df[
        ["protein", "selection", "model_name", "rscc_baseline", "rscc_optimal", "delta_rscc", "label"]
    ].merge(_rmsd, on=["protein", "selection", "model_name"])
    rscc_minrmsd_df["model"] = rscc_minrmsd_df["model_name"].map(lambda m: MODEL_DISPLAY.get(m, m))
    rscc_minrmsd_df["delta_min_rmsd_to_A"] = (
        rscc_minrmsd_df["min_rmsd_to_A_optimal"] - rscc_minrmsd_df["min_rmsd_to_A_baseline"]
    )
    rscc_minrmsd_df["delta_min_rmsd_to_B"] = (
        rscc_minrmsd_df["min_rmsd_to_B_optimal"] - rscc_minrmsd_df["min_rmsd_to_B_baseline"]
    )

    # Analysis 1: RSCC increases AND guided min RMSD to BOTH altlocs < 1.0 A.
    _both_close = (
        rscc_minrmsd_df[
            (rscc_minrmsd_df["delta_rscc"] > 0)
            & (rscc_minrmsd_df["min_rmsd_to_A_optimal"] < 1.0)
            & (rscc_minrmsd_df["min_rmsd_to_B_optimal"] < 1.0)
        ]
        .sort_values(["model", "delta_rscc"], ascending=[True, False])
        .reset_index(drop=True)
    )
    _display = _both_close[
        ["model", "label", "rscc_baseline", "rscc_optimal", "delta_rscc",
         "min_rmsd_to_A_optimal", "min_rmsd_to_B_optimal"]
    ]
    _n_up = int((rscc_minrmsd_df["delta_rscc"] > 0).sum())
    _cov = ", ".join(sorted(rscc_minrmsd_df["model"].unique()))

    mo.vstack([
        mo.md(f"""## RSCC Up & Min RMSD to Both Altlocs < 1.0 A (per model)
        Optimal guidance **raises RSCC** and the guided ensemble's min RMSD to **both**
        altlocs is under 1.0 A (fits density while staying near both conformations; min
        RMSD = min across the 8-member ensemble). {len(_both_close)} of {_n_up} RSCC-up
        selections qualify. Coverage: {_cov} (no rf3-chiral-off min-RMSD CSV yet).
        Per model: {per_model_counts(_both_close)}."""),
        _display,
    ])
    return (rscc_minrmsd_df,)


@app.cell
def _(mo, per_model_counts, rscc_minrmsd_df):
    # Analysis 2: RSCC increases but guided min RMSD to A or B grows vs baseline
    # (better fit while drifting away from at least one conformation).
    _worse = rscc_minrmsd_df[
        (rscc_minrmsd_df["delta_rscc"] > 0)
        & ((rscc_minrmsd_df["delta_min_rmsd_to_A"] > 0)
           | (rscc_minrmsd_df["delta_min_rmsd_to_B"] > 0))
    ].copy()
    _worse["_drift"] = _worse[["delta_min_rmsd_to_A", "delta_min_rmsd_to_B"]].max(axis=1)
    _worse = _worse.sort_values(["model", "_drift"], ascending=[True, False]).reset_index(drop=True)
    _display = _worse[
        ["model", "label", "delta_rscc", "delta_min_rmsd_to_A", "delta_min_rmsd_to_B",
         "min_rmsd_to_A_optimal", "min_rmsd_to_B_optimal"]
    ]
    mo.vstack([
        mo.md(f"""## RSCC Up & Min RMSD to A or B Increases (per model)
        Optimal guidance **raises RSCC** but the guided ensemble's min RMSD to **at least
        one** altloc grows vs baseline -- a better fit bought by moving away from a
        conformation. Sorted by model then largest RMSD increase.
        {len(_worse)} selections qualify. Per model: {per_model_counts(_worse)}."""),
        _display,
    ])
    return


if __name__ == "__main__":
    app.run()
