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
    from analysis_helpers import (
        add_simplified_metadata_columns,
        filter_by_model_group,
        get_rscc_nan_proteins,
        load_rscc,
        load_lddt,
        load_classification,
        make_ridge_plot,
        make_rscc_scatter,
        make_scatter_comparison,
        make_aggregated_scatter,
        make_classified_scatter,
        make_classified_aggregated_scatter,
        make_combined_baseline_distributions,
        make_delta_frame,
        model_slug,
        save_fig,
        save_model_legend,
        save_per_model_figures,
        MODEL_DISPLAY,
        MODELS_ALT,
        MODELS_PROTENIX_GROUP,
        MODELS_RF3_GROUP,
        OPTIMAL_GW,
    )

    return (
        MODELS_PROTENIX_GROUP,
        MODELS_RF3_GROUP,
        MODEL_DISPLAY,
        OPTIMAL_GW,
        add_simplified_metadata_columns,
        filter_by_model_group,
        get_rscc_nan_proteins,
        load_classification,
        load_lddt,
        load_rscc,
        make_aggregated_scatter,
        make_classified_aggregated_scatter,
        make_classified_scatter,
        make_combined_baseline_distributions,
        make_delta_frame,
        make_ridge_plot,
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
    | 04 | RSCC baseline distribution | `04_rscc_baseline_all_models[_variant].svg` |
    | 05 | RSCC baseline vs weak guidance | `05_rscc_baseline_vs_weak_distribution[_variant].svg` |
    | 06 | Delta RSCC distribution | `06_delta_rscc_distribution[_variant].svg` |
    | 07 | Delta RSCC weak (all models) | `07_delta_rscc_weak_all_models.svg` |
    | 08 | Delta RSCC weak (protenix) | `08_delta_rscc_weak_protenix.svg` |
    | -- | RSCC scatter per model | `rscc_scatterplots/rscc_scatter_{model}.svg` |
    | -- | RSCC vs LDDT scatter per model | `cross_metric_scatterplots/rscc_vs_lddt_{model}.svg` |
    | -- | Aggregated RSCC scatter | `rscc_scatterplots/rscc_scatter_aggregated_all_models.svg` |
    | -- | RSCC by classification per model | `rscc_scatterplots/rscc_by_class_{model}.svg` |
    | -- | RSCC by classification aggregated | `rscc_scatterplots/rscc_by_class_aggregated.svg` |
    | -- | Aggregated LDDT scatter | `rscc_scatterplots/lddt_scatter_aggregated_all_models.svg` |
    | -- | Aggregated RSCC vs LDDT scatter | `rscc_scatterplots/rscc_vs_lddt_aggregated_all_models.svg` |
    | 09 | LDDT baseline distribution | `09_lddt_to_ref_baseline_distribution[_variant].svg` |
    | 10 | LDDT all guidance distribution | `10_lddt_to_ref_all_guidance_distribution[_variant].svg` |
    | 11 | Delta LDDT distribution | `11_delta_lddt_to_ref_distribution[_variant].svg` |
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
    mo.hstack([model_group, smoothing, auto_save, exclude_nan_rscc], justify="start", gap=1)
    return auto_save, exclude_nan_rscc, model_group, smoothing


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
    classification_df = load_classification("classify_altloc_regions.csv")
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
        _rscc_data, _lddt_data, smoothing=smoothing.value,
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
def _(
    auto_save,
    filter_by_model_group,
    make_ridge_plot,
    mo,
    model_group,
    rscc_all,
    save_fig,
):
    _data, _suffix = filter_by_model_group(
        rscc_all, model_group.value,
        extra_query="guidance_strength in ('baseline', 'weak')",
    )
    fig_05 = make_ridge_plot(_data, "rscc")
    if auto_save.value:
        save_fig(fig_05, f"05_rscc_baseline_vs_weak_distribution{_suffix}.png")
    mo.vstack([
        mo.md("""## Figure 05 -- RSCC Baseline vs Weak Guidance
        Comparing baseline and weak guidance distributions.
        `figures/05_rscc_baseline_vs_weak_distribution[_variant].svg`"""),
        fig_05,
    ])
    return


@app.cell
def _(
    auto_save,
    filter_by_model_group,
    make_delta_frame,
    make_ridge_plot,
    mo,
    model_group,
    rscc_all,
    save_fig,
):
    _rscc_delta = make_delta_frame(rscc_all, "rscc")
    _data, _suffix = filter_by_model_group(
        _rscc_delta, model_group.value,
        extra_query="delta_rscc < 0.3 and delta_rscc > -0.3",
    )
    fig_06 = make_ridge_plot(_data, "delta_rscc", is_delta=True)
    if auto_save.value:
        save_fig(fig_06, f"06_delta_rscc_distribution{_suffix}.png")
    mo.vstack([
        mo.md("""## Figure 06 -- Delta RSCC Distribution
        Change in RSCC from baseline (clipped to +/-0.3).
        `figures/06_delta_rscc_distribution[_variant].svg`"""),
        fig_06,
    ])
    return


@app.cell
def _(auto_save, make_delta_frame, make_ridge_plot, mo, rscc_all, save_fig):
    _delta = make_delta_frame(rscc_all, "rscc")
    _weak_filt = "guidance_strength == 'weak' and delta_rscc < 0.3 and delta_rscc > -0.3"
    _xlabel = "Change in RSCC (Guided vs. Baseline)"

    fig_07 = make_ridge_plot(
        _delta.query(f"{_weak_filt} and model_name not in @MODELS_ALT"),
        "delta_rscc", custom_xlabel=_xlabel, is_delta=True, combined=True,
    )
    if auto_save.value:
        save_fig(fig_07, "07_delta_rscc_weak_all_models.png")

    fig_08 = make_ridge_plot(
        _delta.query(f"{_weak_filt} and model_name in @MODELS_PROTENIX_GROUP"),
        "delta_rscc", custom_xlabel=_xlabel, is_delta=True, combined=True,
    )
    if auto_save.value:
        save_fig(fig_08, "08_delta_rscc_weak_protenix.png")

    mo.vstack([
        mo.md("""## Figures 07-08 -- Delta RSCC at Weak Guidance
        Fig 07: all main models combined. Fig 08: Protenix variants.
        `figures/07_delta_rscc_weak_all_models.svg`, `figures/08_delta_rscc_weak_protenix.svg`"""),
        fig_07,
        fig_08,
    ])
    return


@app.cell
def _(MODEL_DISPLAY, mo):
    # Build dropdown from MODEL_DISPLAY (consistent labels across all scatter sections)
    _options = {
        v: k for k, v in MODEL_DISPLAY.items()
        if k in ("boltz2X-RAY", "boltz2MD", "protenix", "rf3",
                  "protenix-mini", "protenix-tiny", "rf3-old", "rf3-chiral-off")
    }
    scatter_model = mo.ui.dropdown(options=_options, value="Boltz2-XRD", label="Model")
    save_all_scatter = mo.ui.switch(label="Save all models at once", value=False)
    mo.hstack([scatter_model, save_all_scatter], justify="start", gap=1)
    return save_all_scatter, scatter_model


@app.cell
def _(OPTIMAL_GW, pd, rscc_all):
    _gw_df = pd.DataFrame(list(OPTIMAL_GW.items()), columns=["model_name", "optimal_gw"])

    _baseline = rscc_all.query("guidance_weight == 0.0")[["protein", "model_name", "rscc"]].rename(
        columns={"rscc": "rscc_baseline"}
    )
    _optimal = (
        rscc_all.merge(_gw_df, on="model_name")
        .query("guidance_weight == optimal_gw")
        [["protein", "model_name", "rscc"]]
        .rename(columns={"rscc": "rscc_optimal"})
    )
    rscc_scatter_df = _baseline.merge(_optimal, on=["protein", "model_name"])
    rscc_scatter_df["delta_rscc"] = rscc_scatter_df["rscc_optimal"] - rscc_scatter_df["rscc_baseline"]
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
        x_label="RSCC at Optimal Guidance",
        y_label="Avg. Silhouette to Ref",
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
        .groupby("protein", as_index=False)
        .agg(
            x_med=("rscc_baseline", "mean"),
            x_q1=("rscc_baseline", lambda s: s.mean() - s.std()),
            x_q3=("rscc_baseline", lambda s: s.mean() + s.std()),
            y_med=("rscc_optimal", "mean"),
            y_q1=("rscc_optimal", lambda s: s.mean() - s.std()),
            y_q3=("rscc_optimal", lambda s: s.mean() + s.std()),
            n=("model_name", "count"),
        )
    )
    _rscc_agg_df = _rscc_agg_df[_rscc_agg_df["n"] >= 3]

    _fig = make_aggregated_scatter(
        _rscc_agg_df,
        x_col="x_med", y_col="y_med",
        x_q1_col="x_q1", x_q3_col="x_q3",
        y_q1_col="y_q1", y_q3_col="y_q3",
        x_label="Baseline RSCC (mean +/- SD)",
        y_label="Guided RSCC at optimal gw (mean +/- SD)",
        title="RSCC Improvement Over Baseline",
        draw_diagonal=True, equal_aspect=True,
        show_marginal_kde=True,
        point_color="#377EB8",
        highlight_top_n=3, highlight_decreased=True,
        highlight_proteins=["6b8x", "5mhx"],
    )
    if auto_save.value:
        save_fig(_fig, "rscc_scatterplots/rscc_scatter_aggregated_all_models.png")
    mo.vstack([
        mo.md("""## Aggregated RSCC Scatter -- Collapsed Across Main Models
        Each point is one protein; baseline and guided RSCC collapsed across main models (mean +/- SD).
        `figures/rscc_scatterplots/rscc_scatter_aggregated_all_models.svg`"""),
        _fig,
    ])
    return


@app.cell
def _(
    OPTIMAL_GW,
    auto_save,
    classification_df,
    make_classified_aggregated_scatter,
    make_classified_scatter,
    mo,
    model_slug,
    pd,
    plt,
    rscc,
    save_all_scatter,
    save_fig,
    scatter_model,
):
    # Build rscc_by_class_df
    _gw_df = pd.DataFrame(list(OPTIMAL_GW.items()), columns=["model_name", "optimal_gw"])
    _region = rscc[rscc["selection_name"] != "all altloc"].copy()
    _region_cls = _region.merge(classification_df, on=["protein", "selection"], how="inner")

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

    # Per-model batch save
    if save_all_scatter.value:
        for _model in sorted(_rscc_by_class_df["model_name"].unique()):
            _mdf = _rscc_by_class_df[_rscc_by_class_df["model_name"] == _model]
            if _mdf.empty:
                continue
            _fig = make_classified_scatter(_mdf, _model, y_label=f"RSCC at gw={OPTIMAL_GW[_model]}")
            save_fig(_fig, f"rscc_scatterplots/rscc_by_class_{model_slug(_model)}.svg")
            plt.close(_fig)

    # Show current model
    _cur = scatter_model.value
    _cur_df = _rscc_by_class_df[_rscc_by_class_df["model_name"] == _cur]
    _fig_cur = make_classified_scatter(_cur_df, _cur, y_label=f"RSCC at gw={OPTIMAL_GW[_cur]}")
    if auto_save.value:
        save_fig(_fig_cur, f"rscc_scatterplots/rscc_by_class_{model_slug(_cur)}.svg")

    # Aggregated scatter
    _rscc_by_class_agg = (
        _rscc_by_class_df
        .groupby(["protein", "selection", "classification"], as_index=False)
        .agg(
            x_mean=("rscc_baseline", "mean"),
            x_err_lo=("rscc_baseline", lambda s: s.mean() - s.std()),
            x_err_hi=("rscc_baseline", lambda s: s.mean() + s.std()),
            y_mean=("rscc_optimal", "mean"),
            y_err_lo=("rscc_optimal", lambda s: s.mean() - s.std()),
            y_err_hi=("rscc_optimal", lambda s: s.mean() + s.std()),
            n=("model_name", "count"),
        )
    )
    _rscc_by_class_agg = _rscc_by_class_agg[_rscc_by_class_agg["n"] >= 3]
    _fig_agg = make_classified_aggregated_scatter(_rscc_by_class_agg)
    if auto_save.value:
        save_fig(_fig_agg, "rscc_scatterplots/rscc_by_class_aggregated.png")

    mo.vstack([
        mo.md("""## RSCC by Conformational Type -- Baseline vs Optimal Guidance
        Per-selection scatter colored by conformational change classification.
        `figures/rscc_scatterplots/rscc_by_class_{model}.svg` and `rscc_by_class_aggregated.svg`"""),
        _fig_cur,
        _fig_agg,
    ])
    return


@app.cell
def _(
    OPTIMAL_GW,
    auto_save,
    cross_metric_df,
    lddt_results_1,
    make_aggregated_scatter,
    mo,
    pd,
    save_fig,
):
    _gw_df = pd.DataFrame(list(OPTIMAL_GW.items()), columns=["model_name", "optimal_gw"])

    # Build lddt_scatter_df
    _lddt_nontrivial = lddt_results_1[lddt_results_1["selection_name"] != "all altloc"]
    _baseline_lddt = (
        _lddt_nontrivial.query("guidance_weight == 0.0")
        .groupby(["protein", "model_name"], as_index=False)
        .agg(lddt_baseline=("avg_silhouette_to_ref", "mean"))
    )
    _optimal_lddt = (
        _lddt_nontrivial.merge(_gw_df, on="model_name")
        .query("guidance_weight == optimal_gw")
        .groupby(["protein", "model_name"], as_index=False)
        .agg(lddt_optimal=("avg_silhouette_to_ref", "mean"))
    )
    _lddt_scatter_df = _baseline_lddt.merge(_optimal_lddt, on=["protein", "model_name"])
    _lddt_scatter_df["delta_lddt"] = _lddt_scatter_df["lddt_optimal"] - _lddt_scatter_df["lddt_baseline"]

    # Aggregated LDDT scatter
    _lddt_agg_df = (
        _lddt_scatter_df.query("model_name not in @MODELS_ALT")
        .groupby("protein", as_index=False)
        .agg(
            x_med=("lddt_baseline", "mean"),
            x_q1=("lddt_baseline", lambda s: s.mean() - s.std()),
            x_q3=("lddt_baseline", lambda s: s.mean() + s.std()),
            y_med=("lddt_optimal", "mean"),
            y_q1=("lddt_optimal", lambda s: s.mean() - s.std()),
            y_q3=("lddt_optimal", lambda s: s.mean() + s.std()),
            n=("model_name", "count"),
        )
    )
    _lddt_agg_df = _lddt_agg_df[_lddt_agg_df["n"] >= 3]
    _fig_lddt = make_aggregated_scatter(
        _lddt_agg_df,
        x_col="x_med", y_col="y_med",
        x_q1_col="x_q1", x_q3_col="x_q3",
        y_q1_col="y_q1", y_q3_col="y_q3",
        x_label="Baseline Avg Silhouette to Ref (mean +/- SD)",
        y_label="Guided Avg Silhouette to Ref at optimal gw (mean +/- SD)",
        title="Avg Silhouette to Ref Improvement Over Baseline",
        draw_diagonal=True, equal_aspect=True,
        show_marginal_kde=True, point_color="#377EB8",
        highlight_top_n=3, highlight_decreased=True,
        highlight_proteins=["6b8x", "5mhx"],
    )
    if auto_save.value:
        save_fig(_fig_lddt, "rscc_scatterplots/silhouette_scatter_aggregated_all_models.png")

    # Aggregated cross-metric scatter
    _cross_agg_df = (
        cross_metric_df.query("model_name not in @MODELS_ALT")
        .groupby("protein", as_index=False)
        .agg(
            x_med=("rscc", "mean"),
            x_q1=("rscc", lambda s: s.mean() - s.std()),
            x_q3=("rscc", lambda s: s.mean() + s.std()),
            y_med=("avg_silhouette_to_ref", "mean"),
            y_q1=("avg_silhouette_to_ref", lambda s: s.mean() - s.std()),
            y_q3=("avg_silhouette_to_ref", lambda s: s.mean() + s.std()),
            n=("model_name", "count"),
        )
    )
    _cross_agg_df = _cross_agg_df[_cross_agg_df["n"] >= 3]
    _fig_cross = make_aggregated_scatter(
        _cross_agg_df,
        x_col="x_med", y_col="y_med",
        x_q1_col="x_q1", x_q3_col="x_q3",
        y_q1_col="y_q1", y_q3_col="y_q3",
        x_label="Guided RSCC at optimal gw (mean +/- SD)",
        y_label="Guided Avg Silhouette to Ref (mean +/- SD)",
        title="RSCC vs Avg Silhouette to Ref at optimal guidance",
        draw_diagonal=False, equal_aspect=False,
        highlight_top_n=3, highlight_decreased=True,
        highlight_proteins=["6b8x", "5mhx"],
    )
    if auto_save.value:
        save_fig(_fig_cross, "rscc_scatterplots/rscc_vs_lddt_aggregated_all_models.png")

    mo.vstack([
        mo.md("""## Aggregated Scatter Plots -- LDDT and Cross-Metric
        Collapsed across main models with mean +/- SD whiskers."""),
        mo.tabs({
            "LDDT Scatter": _fig_lddt,
            "RSCC vs LDDT": _fig_cross,
        }),
    ])
    return


@app.cell
def _():
    # Fig 09 (standalone LDDT baseline) merged into the combined cell above (fig_04_09).
    return


@app.cell
def _(
    MODELS_PROTENIX_GROUP,
    MODELS_RF3_GROUP,
    auto_save,
    lddt_results_1,
    make_ridge_plot,
    mo,
    save_fig,
    smoothing,
):
    _all_lddt = lddt_results_1.query("selection_name != 'all altloc'")

    fig_10 = make_ridge_plot(
        _all_lddt.query("model_name not in @MODELS_ALT"),
        "avg_silhouette_to_ref", smoothing=smoothing.value,
    )
    if auto_save.value:
        save_fig(fig_10, "10_lddt_to_ref_all_guidance_distribution.png")
        for _group, _tag in [(MODELS_PROTENIX_GROUP, "protenix_alts"), (MODELS_RF3_GROUP, "rf3_alts")]:
            _fig = make_ridge_plot(
                _all_lddt[_all_lddt["model_name"].isin(_group)],
                "avg_silhouette_to_ref", smoothing=smoothing.value,
            )
            save_fig(_fig, f"10_lddt_to_ref_all_guidance_distribution_{_tag}.png")

    mo.vstack([
        mo.md("""## Figure 10 -- LDDT-to-Reference Across All Guidance Strengths
        Silhouette scores at baseline, weak, and strong guidance.
        `figures/10_lddt_to_ref_all_guidance_distribution[_variant].svg`"""),
        fig_10,
    ])
    return


@app.cell
def _(
    MODELS_PROTENIX_GROUP,
    MODELS_RF3_GROUP,
    auto_save,
    lddt_merged,
    make_ridge_plot,
    mo,
    save_fig,
):
    fig_11 = make_ridge_plot(
        lddt_merged.query("model_name not in @MODELS_ALT"),
        "delta_silhouette_to_ref", is_delta=True,
    )
    if auto_save.value:
        save_fig(fig_11, "11_delta_lddt_to_ref_distribution.png")
        for _group, _tag in [(MODELS_PROTENIX_GROUP, "protenix_alts"), (MODELS_RF3_GROUP, "rf3_alts")]:
            _fig = make_ridge_plot(
                lddt_merged[lddt_merged["model_name"].isin(_group)],
                "delta_silhouette_to_ref", is_delta=True,
            )
            save_fig(_fig, f"11_delta_lddt_to_ref_distribution_{_tag}.png")

    mo.vstack([
        mo.md("""## Figure 11 -- Delta LDDT-to-Reference Distribution
        Change in silhouette score from baseline.
        `figures/11_delta_lddt_to_ref_distribution[_variant].svg`"""),
        fig_11,
    ])
    return


if __name__ == "__main__":
    app.run()
