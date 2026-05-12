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
        load_rscc,
        merge_classification_by_containment,
    )
    from plot_classified import (
        make_classified_aggregated_scatter,
        make_classified_scatter,
        make_classified_scatter_grid,
    )
    from plot_distributions import make_combined_baseline_distributions
    from plot_scatter import (
        make_aggregated_scatter,
        make_rscc_scatter,
        make_scatter_comparison,
    )
    from style import (
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
        MODELS_ALT,
        MODEL_DISPLAY,
        OPTIMAL_GW,
        add_simplified_metadata_columns,
        filter_by_model_group,
        format_residue_range,
        get_rscc_nan_proteins,
        load_classification,
        load_lddt,
        load_rscc,
        merge_classification_by_containment,
        make_aggregated_scatter,
        make_classified_aggregated_scatter,
        make_classified_scatter,
        make_classified_scatter_grid,
        make_combined_baseline_distributions,
        make_rscc_scatter,
        make_scatter_comparison,
        model_slug,
        pd,
        plt,
        save_fig,
        save_model_legend,
        save_per_model_figures,
    )


@app.cell
def _(mo):
    mo.md(r"""
    # Fit Quality Analysis

    This notebook evaluates how well guided structure predictions fit
    experimental data, using two complementary metrics:

    1. **RSCC** — measures how well predicted density matches observed electron density.
    2. **LDDT-based silhouette score** — measures structural clustering quality
       relative to a reference structure.

    ## Figures produced by this notebook

    | ID | Description | Output |
    |----|-------------|--------|
    | legend | Model color legend | `model_legend.svg` |
    | 04_09 | Combined baseline distributions (RSCC + LDDT) | `04_09_baseline_distributions[_variant].svg` |
    | -- | RSCC scatter per model | `rscc_scatterplots/rscc_scatter_{model}.svg` |
    | -- | RSCC vs LDDT scatter per model | `cross_metric_scatterplots/rscc_vs_lddt_{model}.svg` |
    | -- | Aggregated RSCC scatter | `rscc_scatterplots/rscc_scatter_aggregated_all_models.svg` |
    | -- | RSCC by classification (combined 4×4 grid) | `rscc_scatterplots/rscc_by_class_combined.svg` |
    | -- | Aggregated LDDT scatter | `rscc_scatterplots/silhouette_scatter_aggregated_all_models.svg` |
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
    return auto_save, cumulative, cumulative_style, exclude_nan_rscc, model_group, smoothing


@app.cell
def _(exclude_nan_rscc, get_rscc_nan_proteins):
    nan_rscc_proteins = (
        get_rscc_nan_proteins("initial40results/rscc_results.csv")
        if exclude_nan_rscc.value else set()
    )
    return (nan_rscc_proteins,)


@app.cell
def _(add_simplified_metadata_columns, load_rscc, nan_rscc_proteins):
    rscc = load_rscc("initial40results/rscc_results.csv")
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
def _(add_simplified_metadata_columns, load_lddt, nan_rscc_proteins):
    lddt_results = load_lddt("initial40results/lddt_results.csv")
    lddt_results_1 = add_simplified_metadata_columns(lddt_results)
    if nan_rscc_proteins:
        lddt_results_1 = lddt_results_1[~lddt_results_1["protein"].isin(nan_rscc_proteins)]

    _lddt_baseline = lddt_results_1.query('guidance_strength == "baseline"')[
        ["protein", "model_name", "selection_name", "avg_silhouette_to_ref"]
    ]
    lddt_merged = _lddt_baseline.merge(
        lddt_results_1.query('guidance_strength != "baseline"'),
        on=["protein", "model_name", "selection_name"],
        suffixes=["_baseline", ""],
    )
    lddt_merged["delta_silhouette_to_ref"] = (
        lddt_merged.avg_silhouette_to_ref - lddt_merged.avg_silhouette_to_ref_baseline
    )
    return lddt_merged, lddt_results_1


@app.cell
def _(load_classification, nan_rscc_proteins):
    classification_df = load_classification("small-altlocs-classified-rmsd.csv")
    if nan_rscc_proteins:
        classification_df = classification_df[~classification_df["protein"].isin(nan_rscc_proteins)]
    return (classification_df,)


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
    _lddt_data = lddt_results_1.query(
        "guidance_strength == 'baseline' and selection_name != 'all altloc'"
        " and model_name not in @MODELS_ALT"
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
def _(OPTIMAL_GW, format_residue_range, pd, rscc):
    _gw_df = pd.DataFrame(list(OPTIMAL_GW.items()), columns=["model_name", "optimal_gw"])
    _region = rscc[rscc["selection_name"] != "all altloc"]

    _baseline = (
        _region.query("guidance_weight == 0.0")
        [["protein", "selection", "model_name", "rscc"]]
        .rename(columns={"rscc": "rscc_baseline"})
    )
    _optimal = (
        _region.merge(_gw_df, on="model_name")
        .query("guidance_weight == optimal_gw")
        [["protein", "selection", "model_name", "rscc"]]
        .rename(columns={"rscc": "rscc_optimal"})
    )
    rscc_scatter_df = _baseline.merge(
        _optimal, on=["protein", "selection", "model_name"]
    )
    rscc_scatter_df["delta_rscc"] = (
        rscc_scatter_df["rscc_optimal"] - rscc_scatter_df["rscc_baseline"]
    )

    rscc_scatter_df["label"] = rscc_scatter_df.apply(
        lambda r: f"{str(r['protein']).upper()} {format_residue_range(r['selection'])}", axis=1,
    )
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
def _(OPTIMAL_GW, lddt_results_1, pd, rscc_all):
    _gw_df = pd.DataFrame(list(OPTIMAL_GW.items()), columns=["model_name", "optimal_gw"])

    # RSCC at optimal guidance weight per model
    _rscc_at_opt = (
        rscc_all.merge(_gw_df, on="model_name")
        .query("guidance_weight == optimal_gw")
        [["protein", "model_name", "rscc"]]
    )

    # LDDT at optimal guidance weight, averaged across selections (excl. "all altloc")
    _lddt_filtered = lddt_results_1[lddt_results_1["selection_name"] != "all altloc"]
    _lddt_at_opt = (
        _lddt_filtered.merge(_gw_df, on="model_name")
        .query("guidance_weight == optimal_gw")
        .groupby(["protein", "model_name"], as_index=False)
        .agg(avg_silhouette_to_ref=("avg_silhouette_to_ref", "mean"))
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
        x_col="rscc", y_col="avg_silhouette_to_ref",
        x_label="RSCC",
        y_label="Clustering Score",
        draw_regression=True, label_all=True, equal_aspect=False,
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
        Per-model scatter comparing RSCC with LDDT silhouette-to-ref at optimal guidance.
        `figures/cross_metric_scatterplots/rscc_vs_lddt_{model}.svg`"""),
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
        highlight_proteins=["6B8X 181-183", "5MHX 162-164"],
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
    OPTIMAL_GW,
    auto_save,
    classification_df,
    load_rscc,
    make_classified_aggregated_scatter,
    make_classified_scatter,
    make_classified_scatter_grid,
    merge_classification_by_containment,
    mo,
    model_slug,
    nan_rscc_proteins,
    pd,
    plt,
    rscc,
    save_all_scatter,
    save_fig,
    scatter_model,
):
    # Build rscc_by_class_df for non-domain-shift classes (small-segment data)
    _gw_df = pd.DataFrame(list(OPTIMAL_GW.items()), columns=["model_name", "optimal_gw"])
    _region = rscc[rscc["selection_name"] != "all altloc"].copy()
    _region_cls = merge_classification_by_containment(_region, classification_df)
    _region_cls = _region_cls[_region_cls["classification"] != "domain_shift"]

    _baseline = (
        _region_cls.query("guidance_weight == 0.0")
        [["protein", "selection", "model_name", "rscc", "classification"]]
        .rename(columns={"rscc": "rscc_baseline"})
    )
    _optimal = (
        _region_cls.merge(_gw_df, on="model_name")
        .query("guidance_weight == optimal_gw")
        [["protein", "selection", "model_name", "rscc"]]
        .rename(columns={"rscc": "rscc_optimal"})
    )
    _rscc_by_class_df = _baseline.merge(_optimal, on=["protein", "selection", "model_name"])
    _rscc_by_class_df["delta_rscc"] = _rscc_by_class_df["rscc_optimal"] - _rscc_by_class_df["rscc_baseline"]

    # Domain shift uses whole-domain RSCC from the legacy backup CSV; the
    # small-segment file in initial40results/ only contains 1-3 residue
    # altloc-resolved selections inside the domain, which under-samples the
    # conformational change for this category.
    _rscc_legacy = load_rscc("initial40results/rscc_results.csv.bak")
    if nan_rscc_proteins:
        _rscc_legacy = _rscc_legacy[~_rscc_legacy["protein"].isin(nan_rscc_proteins)]
    _ds_cls = classification_df.query("classification == 'domain_shift'")[
        ["protein", "selection", "classification"]
    ]
    _legacy_region = _rscc_legacy.merge(_ds_cls, on=["protein", "selection"], how="inner")

    _baseline_ds = (
        _legacy_region.query("guidance_weight == 0.0")
        [["protein", "selection", "model_name", "rscc", "classification"]]
        .rename(columns={"rscc": "rscc_baseline"})
    )
    _optimal_ds = (
        _legacy_region.merge(_gw_df, on="model_name")
        .query("guidance_weight == optimal_gw")
        [["protein", "selection", "model_name", "rscc"]]
        .rename(columns={"rscc": "rscc_optimal"})
    )
    _rscc_ds_df = _baseline_ds.merge(_optimal_ds, on=["protein", "selection", "model_name"])
    _rscc_ds_df["delta_rscc"] = _rscc_ds_df["rscc_optimal"] - _rscc_ds_df["rscc_baseline"]

    # Per-model batch save (non-domain-shift classes)
    if save_all_scatter.value:
        for _model in sorted(_rscc_by_class_df["model_name"].unique()):
            _mdf = _rscc_by_class_df[_rscc_by_class_df["model_name"] == _model]
            if _mdf.empty:
                continue
            _fig = make_classified_scatter(_mdf, _model)
            save_fig(_fig, f"rscc_scatterplots/rscc_by_class_{model_slug(_model)}.svg")
            plt.close(_fig)
        for _model in sorted(_rscc_ds_df["model_name"].unique()):
            _mdf_ds = _rscc_ds_df[_rscc_ds_df["model_name"] == _model]
            if _mdf_ds.empty:
                continue
            _fig_ds = make_classified_scatter(_mdf_ds, _model)
            save_fig(_fig_ds, f"rscc_scatterplots/rscc_by_class_domain_shift_{model_slug(_model)}.svg")
            plt.close(_fig_ds)

    # Show current model
    _cur = scatter_model.value
    _cur_df = _rscc_by_class_df[_rscc_by_class_df["model_name"] == _cur]
    _fig_cur = make_classified_scatter(_cur_df, _cur)
    if auto_save.value:
        save_fig(_fig_cur, f"rscc_scatterplots/rscc_by_class_{model_slug(_cur)}.svg")

    _cur_ds_df = _rscc_ds_df[_rscc_ds_df["model_name"] == _cur]
    _fig_cur_ds = (
        make_classified_scatter(_cur_ds_df, _cur) if not _cur_ds_df.empty else None
    )
    if auto_save.value and _fig_cur_ds is not None:
        save_fig(_fig_cur_ds, f"rscc_scatterplots/rscc_by_class_domain_shift_{model_slug(_cur)}.svg")

    # Aggregated scatter (non-domain-shift)
    _rscc_by_class_agg = (
        _rscc_by_class_df
        .groupby(["protein", "selection", "classification"], as_index=False)
        .agg(
            x_mean=("rscc_baseline", "median"),
            y_mean=("rscc_optimal", "median"),
            n=("model_name", "count"),
        )
    )
    _rscc_by_class_agg = _rscc_by_class_agg[_rscc_by_class_agg["n"] >= 3]
    _fig_agg = make_classified_aggregated_scatter(
        _rscc_by_class_agg,
        x_err_lo_col=None, x_err_hi_col=None,
        y_err_lo_col=None, y_err_hi_col=None,
        x_label="Baseline RSCC (median across models)",
        y_label="Median Improvement (Guided - Baseline)",
        show_errors=False,
    )
    if auto_save.value:
        save_fig(_fig_agg, "rscc_scatterplots/rscc_by_class_aggregated.png")

    # Aggregated scatter (domain-shift)
    _rscc_ds_agg = (
        _rscc_ds_df
        .groupby(["protein", "selection", "classification"], as_index=False)
        .agg(
            x_mean=("rscc_baseline", "median"),
            y_mean=("rscc_optimal", "median"),
            n=("model_name", "count"),
        )
    )
    _rscc_ds_agg = _rscc_ds_agg[_rscc_ds_agg["n"] >= 3]
    _fig_agg_ds = (
        make_classified_aggregated_scatter(
            _rscc_ds_agg,
            x_err_lo_col=None, x_err_hi_col=None,
            y_err_lo_col=None, y_err_hi_col=None,
            x_label="Baseline RSCC (median across models)",
            y_label="Median Improvement (Guided - Baseline)",
            show_errors=False,
        ) if not _rscc_ds_agg.empty else None
    )
    if auto_save.value and _fig_agg_ds is not None:
        save_fig(_fig_agg_ds, "rscc_scatterplots/rscc_by_class_domain_shift_aggregated.png")

    # Combined grid: rows = main models only, cols = classes, shared x and y axes.
    _fig_grid = make_classified_scatter_grid(
        _rscc_by_class_df.query("model_name not in @MODELS_ALT"),
    )
    if auto_save.value:
        save_fig(_fig_grid, "rscc_scatterplots/rscc_by_class_grid_all_models.svg")

    _ds_grid_data = _rscc_ds_df.query("model_name not in @MODELS_ALT")
    _fig_grid_ds = (
        make_classified_scatter_grid(_ds_grid_data) if not _ds_grid_data.empty else None
    )
    if auto_save.value and _fig_grid_ds is not None:
        save_fig(_fig_grid_ds, "rscc_scatterplots/rscc_by_class_domain_shift_grid_all_models.svg")

    _stack_items = [
        mo.md("""## RSCC by Conformational Type -- Baseline vs Optimal Guidance
        Per-selection scatter colored by conformational change classification.
        Domain-shift entries use the whole-domain selection RSCC (from the
        legacy backup CSV) and are plotted separately below.
        `figures/rscc_scatterplots/rscc_by_class_{model}.svg`, `rscc_by_class_aggregated.svg`,
        `rscc_by_class_grid_all_models.svg`,
        `rscc_by_class_domain_shift_{model}.svg`, `rscc_by_class_domain_shift_aggregated.svg`,
        and `rscc_by_class_domain_shift_grid_all_models.svg`"""),
        _fig_cur,
        _fig_agg,
        _fig_grid,
    ]
    if _fig_cur_ds is not None:
        _stack_items.append(_fig_cur_ds)
    if _fig_agg_ds is not None:
        _stack_items.append(_fig_agg_ds)
    if _fig_grid_ds is not None:
        _stack_items.append(_fig_grid_ds)
    mo.vstack(_stack_items)
    return


@app.cell
def _(
    OPTIMAL_GW,
    auto_save,
    format_residue_range,
    lddt_results_1,
    make_aggregated_scatter,
    mo,
    pd,
    save_fig,
):
    _gw_df = pd.DataFrame(list(OPTIMAL_GW.items()), columns=["model_name", "optimal_gw"])

    # Build per-(protein, selection, model) baseline + optimal LDDT
    _lddt_nontrivial = lddt_results_1[lddt_results_1["selection_name"] != "all altloc"]
    _baseline_lddt = (
        _lddt_nontrivial.query("guidance_weight == 0.0 and model_name not in @MODELS_ALT")
        [["protein", "selection", "model_name", "avg_silhouette_to_ref"]]
        .rename(columns={"avg_silhouette_to_ref": "lddt_baseline"})
    )
    _optimal_lddt = (
        _lddt_nontrivial.merge(_gw_df, on="model_name")
        .query("guidance_weight == optimal_gw and model_name not in @MODELS_ALT")
        [["protein", "selection", "model_name", "avg_silhouette_to_ref"]]
        .rename(columns={"avg_silhouette_to_ref": "lddt_optimal"})
    )
    _lddt_per_sel = _baseline_lddt.merge(
        _optimal_lddt, on=["protein", "selection", "model_name"]
    )

    # Aggregated LDDT scatter (per protein, selection)
    _lddt_agg_df = (
        _lddt_per_sel.groupby(["protein", "selection"], as_index=False)
        .agg(
            x_med=("lddt_baseline", "median"),
            y_med=("lddt_optimal", "median"),
            n=("model_name", "count"),
        )
    )
    _lddt_agg_df = _lddt_agg_df[_lddt_agg_df["n"] >= 3]

    _lddt_agg_df["label"] = _lddt_agg_df.apply(
        lambda r: f"{str(r['protein']).upper()} {format_residue_range(r['selection'])}", axis=1,
    )

    _fig_lddt = make_aggregated_scatter(
        _lddt_agg_df,
        x_col="x_med", y_col="y_med",
        x_label="Baseline Clustering Score (median across models)",
        y_label="Median Improvement (Guided - Baseline)",
        title="Clustering Score Improvement Over Baseline",
        label_col="label",
        residual=True,
        show_marginal_kde=True, show_errors=False, point_color="#377EB8",
        label_top_n=3, label_bottom_n=5, label_far_from_ceiling_n=3,
        highlight_proteins=["6B8X 181-183", "5MHX 162-164"],
        xlim=(0.0, 1.0),
    )
    if auto_save.value:
        save_fig(_fig_lddt, "rscc_scatterplots/silhouette_scatter_aggregated_all_models.png")

    mo.vstack([
        mo.md("""## Aggregated Clustering Score Scatter
        Collapsed across main models using the per-(protein, selection) median; per-model spread is shown in the per-model figures."""),
        _fig_lddt,
    ])
    return


if __name__ == "__main__":
    app.run()
