import marimo

__generated_with = "0.23.0"
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
        load_rscc,
        load_lddt,
        make_ridge_plot,
        make_rscc_scatter,
        make_scatter_comparison,
        make_aggregated_scatter,
        make_delta_frame,
        save_fig,
        save_model_legend,
        MODEL_PALETTE,
        MODEL_DISPLAY,
        MODELS_ALT,
        MODELS_PROTENIX_GROUP,
        MODELS_RF3_GROUP,
    )

    return (
        MODELS_PROTENIX_GROUP,
        MODELS_RF3_GROUP,
        MODEL_DISPLAY,
        add_simplified_metadata_columns,
        load_lddt,
        load_rscc,
        make_aggregated_scatter,
        make_delta_frame,
        make_ridge_plot,
        make_rscc_scatter,
        make_scatter_comparison,
        pd,
        plt,
        save_fig,
        save_model_legend,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Fit Quality Analysis

    This notebook evaluates how well guided structure predictions fit
    experimental data, using two complementary metrics:

    1. **RSCC (Real Space Correlation Coefficient)** — measures how well the
       predicted model density matches the observed electron density map.
       Values closer to 1.0 indicate better fit.
    2. **LDDT-based silhouette score** — measures structural clustering quality
       relative to a reference structure, capturing whether guidance pushes
       predictions closer to the experimentally determined conformation.

    Each metric is examined at baseline (no guidance), across guidance strengths,
    and as deltas from baseline to quantify the effect of guidance.
    """)
    return


@app.cell
def _(mo):
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
    mo.hstack([model_group, smoothing, auto_save], justify="start", gap=1)
    return auto_save, model_group, smoothing


@app.cell
def _(load_rscc):
    rscc = load_rscc("initial40results/rscc_results.csv")
    return (rscc,)


@app.cell
def _(add_simplified_metadata_columns, rscc):
    rscc_all = rscc.groupby(["protein", "model_name", "guidance_weight"]).last().reset_index()
    rscc_all = add_simplified_metadata_columns(rscc_all)
    return (rscc_all,)


@app.cell
def _(load_lddt):
    lddt_results = load_lddt("initial40results/lddt_results_w_baseline.csv")
    return (lddt_results,)


@app.cell
def _(add_simplified_metadata_columns, lddt_results):
    lddt_results_1 = add_simplified_metadata_columns(lddt_results)
    return (lddt_results_1,)


@app.cell
def _(lddt_results_1):
    lddt_baseline = lddt_results_1.query('guidance_strength == "baseline"')[
        ["protein", "model_name", "selection_name", "avg_silhouette_to_ref"]
    ]
    lddt_merged = lddt_baseline.merge(
        lddt_results_1.query('guidance_strength != "baseline"'),
        on=["protein", "model_name", "selection_name"],
        suffixes=["_baseline", ""],
    )
    lddt_merged["delta_silhouette_to_ref"] = (
        lddt_merged.avg_silhouette_to_ref - lddt_merged.avg_silhouette_to_ref_baseline
    )
    return (lddt_merged,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Model Legend

    A standalone color legend mapping each model family to its plot color,
    used as a reference across all figures.

    **Output:** `figures/model_legend.svg`
    """)
    return


@app.cell
def _(auto_save, save_model_legend):
    if auto_save.value:
        save_model_legend()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Figure 04 — RSCC Baseline Distribution

    Ridge plot of RSCC values at the no-guidance baseline, showing all models
    overlaid on a single combined axis. This establishes the starting quality
    of each model before any guidance is applied.

    **Output:** `figures/04_rscc_baseline_all_models[_variant].svg`
    """)
    return


@app.cell
def _(auto_save, make_ridge_plot, model_group, rscc_all, save_fig):
    _group_val = model_group.value
    if _group_val == "main":
        _data = rscc_all.query("guidance_strength == 'baseline' and model_name not in @MODELS_ALT")
        _suffix = ""
    elif _group_val == "protenix":
        _data = rscc_all.query("guidance_strength == 'baseline' and model_name in @MODELS_PROTENIX_GROUP")
        _suffix = "_protenix_alts"
    else:
        _data = rscc_all.query("guidance_strength == 'baseline' and model_name in @MODELS_RF3_GROUP")
        _suffix = "_rf3_alts"
    fig_04 = make_ridge_plot(_data, "rscc", combined=True)
    if auto_save.value:
        save_fig(fig_04, f"04_rscc_baseline_all_models{_suffix}.png")
    fig_04
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Figure 05 — RSCC Baseline vs Weak Guidance

    Ridge plot comparing RSCC distributions at baseline and weak guidance for
    each model. This reveals how much (if any) density-fit improvement
    comes from low-strength guidance.

    **Output:** `figures/05_rscc_baseline_vs_weak_distribution[_variant].svg`
    """)
    return


@app.cell
def _(auto_save, make_ridge_plot, model_group, rscc_all, save_fig):
    _group_val = model_group.value
    if _group_val == "main":
        _data = rscc_all.query("guidance_strength in ('baseline', 'weak') and model_name not in @MODELS_ALT")
        _suffix = ""
    elif _group_val == "protenix":
        _data = rscc_all.query("guidance_strength in ('baseline', 'weak') and model_name in @MODELS_PROTENIX_GROUP")
        _suffix = "_protenix_alts"
    else:
        _data = rscc_all.query("guidance_strength in ('baseline', 'weak') and model_name in @MODELS_RF3_GROUP")
        _suffix = "_rf3_alts"
    fig_05 = make_ridge_plot(_data, "rscc")
    if auto_save.value:
        save_fig(fig_05, f"05_rscc_baseline_vs_weak_distribution{_suffix}.png")
    fig_05
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Figure 06 — Delta RSCC Distribution

    Ridge plot of the change in RSCC from baseline across all guidance strengths.
    Positive delta means guidance improved density fit. Values are clipped to
    ±0.3 for readability.

    **Output:** `figures/06_delta_rscc_distribution[_variant].svg`
    """)
    return


@app.cell
def _(
    auto_save,
    make_delta_frame,
    make_ridge_plot,
    model_group,
    rscc_all,
    save_fig,
):
    rscc_delta = make_delta_frame(rscc_all, "rscc")
    _group_val = model_group.value
    if _group_val == "main":
        _data = rscc_delta.query("delta_rscc < 0.3 and delta_rscc > -0.3 and model_name not in @MODELS_ALT")
        _suffix = ""
    elif _group_val == "protenix":
        _data = rscc_delta.query("delta_rscc < 0.3 and delta_rscc > -0.3 and model_name in @MODELS_PROTENIX_GROUP")
        _suffix = "_protenix_alts"
    else:
        _data = rscc_delta.query("delta_rscc < 0.3 and delta_rscc > -0.3 and model_name in @MODELS_RF3_GROUP")
        _suffix = "_rf3_alts"
    fig_06 = make_ridge_plot(_data, "delta_rscc", is_delta=True)
    if auto_save.value:
        save_fig(fig_06, f"06_delta_rscc_distribution{_suffix}.png")
    fig_06
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Figure 07 — Delta RSCC (Weak Guidance, All Models Combined)

    Ridge plot of RSCC change at weak guidance only, with all main models shown
    on a single combined axis. This focuses on the most practically relevant
    guidance level.

    **Output:** `figures/07_delta_rscc_weak_all_models.svg`
    """)
    return


@app.cell
def _(auto_save, make_delta_frame, make_ridge_plot, rscc_all, save_fig):
    _delta = make_delta_frame(rscc_all, "rscc")
    fig_07 = make_ridge_plot(
        _delta.query(
            "guidance_strength == 'weak'"
            " and delta_rscc < 0.3 and delta_rscc > -0.3"
            " and model_name not in @MODELS_ALT"
        ),
        "delta_rscc",
        custom_xlabel="∆ Real Space Correlation Coefficient (Guided vs. Baseline)",
        is_delta=True,
        combined=True,
    )
    if auto_save.value:
        save_fig(fig_07, "07_delta_rscc_weak_all_models.png")
    fig_07
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Figure 08 — Delta RSCC (Weak Guidance, Protenix Group)

    Same as Figure 07 but restricted to Protenix and its size variants,
    highlighting how model capacity affects density-fit improvement under
    weak guidance.

    **Output:** `figures/08_delta_rscc_weak_protenix.svg`
    """)
    return


@app.cell
def _(auto_save, make_delta_frame, make_ridge_plot, rscc_all, save_fig):
    _delta = make_delta_frame(rscc_all, "rscc")
    fig_08 = make_ridge_plot(
        _delta.query(
            "guidance_strength == 'weak'"
            " and delta_rscc < 0.3 and delta_rscc > -0.3"
            " and model_name in @MODELS_PROTENIX_GROUP"
        ),
        "delta_rscc",
        custom_xlabel="∆ Real Space Correlation Coefficient (Guided vs. Baseline)",
        is_delta=True,
        combined=True,
    )
    if auto_save.value:
        save_fig(fig_08, "08_delta_rscc_weak_protenix.png")
    fig_08
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## RSCC Scatter Plots — Baseline vs Optimal Guidance

    Per-model scatter plots comparing baseline RSCC (x-axis) with RSCC at the
    optimal guidance weight (y-axis). The optimal weight is **0.1** for
    Boltz2 and Protenix models and **0.01** for RF3 models.

    Points above the diagonal improved with guidance; points below degraded.
    The top-3 most-improved proteins and any that decreased are highlighted.

    **Output:** `figures/rscc_scatterplots/rscc_scatter_{model}.svg`
    """)
    return


@app.cell
def _(mo):
    scatter_model = mo.ui.dropdown(
        options={
            "Boltz2 X-RAY": "boltz2X-RAY",
            "Boltz2 MD": "boltz2MD",
            "Protenix": "protenix",
            "RF3": "rf3",
            "Protenix-Mini": "protenix-mini",
            "Protenix-Tiny": "protenix-tiny",
            "RF3-Old": "rf3-old",
            "RF3-Chiral-Off": "rf3-chiral-off",
        },
        value="Boltz2 X-RAY",
        label="Model",
    )
    save_all_scatter = mo.ui.switch(label="Save all models at once", value=False)
    mo.hstack([scatter_model, save_all_scatter], justify="start", gap=1)
    return save_all_scatter, scatter_model


@app.cell
def _(pd, rscc_all):
    OPTIMAL_GW = {
        "boltz2X-RAY": 0.1,
        "boltz2MD": 0.1,
        "protenix": 0.1,
        "protenix-mini": 0.1,
        "protenix-tiny": 0.1,
        "rf3": 0.01,
        "rf3-old": 0.01,
        "rf3-chiral-off": 0.01,
    }

    _baseline = rscc_all.query("guidance_weight == 0.0")[["protein", "model_name", "rscc"]].rename(
        columns={"rscc": "rscc_baseline"}
    )

    _rows = []
    for _model, _gw in OPTIMAL_GW.items():
        _mask = (rscc_all["model_name"] == _model) & (rscc_all["guidance_weight"] == _gw)
        _opt = rscc_all.loc[_mask, ["protein", "model_name", "rscc"]].rename(
            columns={"rscc": "rscc_optimal"}
        )
        _rows.append(_opt)
    _optimal = pd.concat(_rows, ignore_index=True)

    rscc_scatter_df = _baseline.merge(_optimal, on=["protein", "model_name"])
    rscc_scatter_df["delta_rscc"] = rscc_scatter_df["rscc_optimal"] - rscc_scatter_df["rscc_baseline"]
    return OPTIMAL_GW, rscc_scatter_df


@app.cell
def _(
    MODEL_DISPLAY,
    OPTIMAL_GW,
    auto_save,
    make_rscc_scatter,
    plt,
    rscc_scatter_df,
    save_all_scatter,
    save_fig,
    scatter_model,
):
    _all_models = sorted(rscc_scatter_df["model_name"].unique())

    if save_all_scatter.value:
        for _model in _all_models:
            _model_df = rscc_scatter_df[rscc_scatter_df["model_name"] == _model]
            if _model_df.empty:
                continue
            _slug = MODEL_DISPLAY.get(_model, _model).lower().replace("-", "_").replace(" ", "_")
            _fig = make_rscc_scatter(_model_df, _model, OPTIMAL_GW[_model])
            save_fig(_fig, f"rscc_scatterplots/rscc_scatter_{_slug}.png")
            plt.close(_fig)

    _cur = scatter_model.value
    _cur_df = rscc_scatter_df[rscc_scatter_df["model_name"] == _cur]
    rscc_scatter_fig = make_rscc_scatter(_cur_df, _cur, OPTIMAL_GW[_cur])
    if auto_save.value:
        _slug = MODEL_DISPLAY.get(_cur, _cur).lower().replace("-", "_").replace(" ", "_")
        save_fig(rscc_scatter_fig, f"rscc_scatterplots/rscc_scatter_{_slug}.png")
    rscc_scatter_fig
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## RSCC vs LDDT Scatter Plots — Cross-Metric Comparison

    Per-model scatter plots comparing RSCC (x-axis) with LDDT silhouette-to-ref
    (y-axis), both measured at the optimal guidance weight. This reveals whether
    proteins that improve in density fit also improve in structural similarity.

    LDDT values are averaged across non-trivial residue selections (excluding
    "all altloc") per protein.

    **Output:** `figures/cross_metric_scatterplots/rscc_vs_lddt_{model}.svg`
    """)
    return


@app.cell
def _(mo):
    cross_scatter_model = mo.ui.dropdown(
        options={
            "Boltz2 X-RAY": "boltz2X-RAY",
            "Boltz2 MD": "boltz2MD",
            "Protenix": "protenix",
            "RF3": "rf3",
            "Protenix-Mini": "protenix-mini",
            "Protenix-Tiny": "protenix-tiny",
            "RF3-Old": "rf3-old",
            "RF3-Chiral-Off": "rf3-chiral-off",
        },
        value="Boltz2 X-RAY",
        label="Model",
    )
    save_all_cross_scatter = mo.ui.switch(label="Save all models at once", value=False)
    mo.hstack([cross_scatter_model, save_all_cross_scatter], justify="start", gap=1)
    return cross_scatter_model, save_all_cross_scatter


@app.cell
def _(OPTIMAL_GW, lddt_results_1, pd, rscc_all):
    # RSCC at optimal guidance weight per model
    _rscc_rows = []
    for _model, _gw in OPTIMAL_GW.items():
        _mask = (rscc_all["model_name"] == _model) & (rscc_all["guidance_weight"] == _gw)
        _rscc_rows.append(rscc_all.loc[_mask, ["protein", "model_name", "rscc"]])
    _rscc_at_opt = pd.concat(_rscc_rows, ignore_index=True)

    # LDDT at optimal guidance weight, averaged across selections (excl. "all altloc")
    _lddt_filtered = lddt_results_1[lddt_results_1["selection_name"] != "all altloc"]
    _lddt_rows = []
    for _model, _gw in OPTIMAL_GW.items():
        _mask = (_lddt_filtered["model_name"] == _model) & (_lddt_filtered["guidance_weight"] == _gw)
        _agg = (
            _lddt_filtered.loc[_mask]
            .groupby(["protein", "model_name"], as_index=False)
            .agg(avg_silhouette_to_ref=("avg_silhouette_to_ref", "mean"))
        )
        _lddt_rows.append(_agg)
    _lddt_at_opt = pd.concat(_lddt_rows, ignore_index=True)

    cross_metric_df = _rscc_at_opt.merge(_lddt_at_opt, on=["protein", "model_name"])
    return (cross_metric_df,)


@app.cell
def _(
    MODEL_DISPLAY,
    auto_save,
    cross_metric_df,
    cross_scatter_model,
    make_scatter_comparison,
    plt,
    save_all_cross_scatter,
    save_fig,
):
    _all_models = sorted(cross_metric_df["model_name"].unique())

    if save_all_cross_scatter.value:
        for _model in _all_models:
            _model_df = cross_metric_df[cross_metric_df["model_name"] == _model]
            if _model_df.empty:
                continue
            _slug = MODEL_DISPLAY.get(_model, _model).lower().replace("-", "_").replace(" ", "_")
            _fig = make_scatter_comparison(
                _model_df, _model,
                x_col="rscc", y_col="avg_silhouette_to_ref",
                x_label="RSCC at Optimal Guidance",
                y_label="Avg. Silhouette to Ref",
                draw_regression=True, show_correlation=True,
                label_all=True, equal_aspect=False,
            )
            save_fig(_fig, f"cross_metric_scatterplots/rscc_vs_lddt_{_slug}.png")
            plt.close(_fig)

    _cur = cross_scatter_model.value
    _cur_df = cross_metric_df[cross_metric_df["model_name"] == _cur]
    cross_scatter_fig = make_scatter_comparison(
        _cur_df, _cur,
        x_col="rscc", y_col="avg_silhouette_to_ref",
        x_label="RSCC at Optimal Guidance",
        y_label="Avg. Silhouette to Ref",
        draw_regression=True, show_correlation=True,
        label_all=True, equal_aspect=False,
    )
    if auto_save.value:
        _slug = MODEL_DISPLAY.get(_cur, _cur).lower().replace("-", "_").replace(" ", "_")
        save_fig(cross_scatter_fig, f"cross_metric_scatterplots/rscc_vs_lddt_{_slug}.png")
    cross_scatter_fig
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Aggregated RSCC Scatter — collapsed across all models

    Each point is one protein; baseline and optimally-guided RSCC values are
    collapsed across all 8 models. The marker sits at the median across models
    and the asymmetric whiskers span the interquartile range (Q1–Q3) in both
    axes — identical error-bar semantics to the aggregated spiral figures.
    Proteins above the y=x diagonal improve under guidance on median across
    models.

    **Output:** `figures/rscc_scatterplots/rscc_scatter_aggregated_all_models.svg`
    """)
    return


@app.cell
def _(auto_save, make_aggregated_scatter, rscc_scatter_df, save_fig):
    rscc_agg_df = (
        rscc_scatter_df
        .groupby("protein", as_index=False)
        .agg(
            x_med=("rscc_baseline", "median"),
            x_q1=("rscc_baseline", lambda s: s.quantile(0.25)),
            x_q3=("rscc_baseline", lambda s: s.quantile(0.75)),
            y_med=("rscc_optimal", "median"),
            y_q1=("rscc_optimal", lambda s: s.quantile(0.25)),
            y_q3=("rscc_optimal", lambda s: s.quantile(0.75)),
            n=("model_name", "count"),
        )
    )
    rscc_agg_df = rscc_agg_df[rscc_agg_df["n"] >= 3]

    fig_rscc_agg = make_aggregated_scatter(
        rscc_agg_df,
        x_col="x_med", y_col="y_med",
        x_q1_col="x_q1", x_q3_col="x_q3",
        y_q1_col="y_q1", y_q3_col="y_q3",
        x_label="Baseline RSCC (median across models)",
        y_label="Guided RSCC at optimal gw (median across models)",
        title="RSCC Improvement Over Baseline Aggregated Across Models",
        draw_diagonal=True,
        equal_aspect=True,
        label_points=True,
    )
    if auto_save.value:
        save_fig(fig_rscc_agg, "rscc_scatterplots/rscc_scatter_aggregated_all_models.png")
    fig_rscc_agg
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Aggregated LDDT Scatter — collapsed across all models

    Analogue of the aggregated RSCC scatter for LDDT silhouette-to-ref.
    For each (protein, model) we average the silhouette across non-trivial
    residue selections (excluding ``all altloc``), take the baseline value
    and the value at the model's optimal guidance weight, then aggregate
    across models per protein with median ± IQR whiskers. Proteins above
    the y=x diagonal improve under guidance on median across models.

    Not all values exist in the source data yet; proteins missing from too
    many models (n<3) are dropped automatically.

    **Output:** `figures/rscc_scatterplots/lddt_scatter_aggregated_all_models.svg`
    """)
    return


@app.cell
def _(OPTIMAL_GW, lddt_results_1, pd):
    # One row per (protein, model): baseline and optimally-guided LDDT
    # silhouette-to-ref, averaged across non-trivial residue selections.
    _lddt_nontrivial = lddt_results_1[lddt_results_1["selection_name"] != "all altloc"]

    _baseline_lddt = (
        _lddt_nontrivial.query("guidance_weight == 0.0")
        .groupby(["protein", "model_name"], as_index=False)
        .agg(lddt_baseline=("avg_silhouette_to_ref", "mean"))
    )

    _opt_rows = []
    for _model, _gw in OPTIMAL_GW.items():
        _mask = (
            (_lddt_nontrivial["model_name"] == _model)
            & (_lddt_nontrivial["guidance_weight"] == _gw)
        )
        _opt = (
            _lddt_nontrivial.loc[_mask]
            .groupby(["protein", "model_name"], as_index=False)
            .agg(lddt_optimal=("avg_silhouette_to_ref", "mean"))
        )
        _opt_rows.append(_opt)
    _optimal_lddt = pd.concat(_opt_rows, ignore_index=True)

    lddt_scatter_df = _baseline_lddt.merge(_optimal_lddt, on=["protein", "model_name"])
    lddt_scatter_df["delta_lddt"] = (
        lddt_scatter_df["lddt_optimal"] - lddt_scatter_df["lddt_baseline"]
    )
    return (lddt_scatter_df,)


@app.cell
def _(auto_save, lddt_scatter_df, make_aggregated_scatter, save_fig):
    lddt_agg_df = (
        lddt_scatter_df
        .groupby("protein", as_index=False)
        .agg(
            x_med=("lddt_baseline", "median"),
            x_q1=("lddt_baseline", lambda s: s.quantile(0.25)),
            x_q3=("lddt_baseline", lambda s: s.quantile(0.75)),
            y_med=("lddt_optimal", "median"),
            y_q1=("lddt_optimal", lambda s: s.quantile(0.25)),
            y_q3=("lddt_optimal", lambda s: s.quantile(0.75)),
            n=("model_name", "count"),
        )
    )
    lddt_agg_df = lddt_agg_df[lddt_agg_df["n"] >= 3]

    fig_lddt_agg = make_aggregated_scatter(
        lddt_agg_df,
        x_col="x_med", y_col="y_med",
        x_q1_col="x_q1", x_q3_col="x_q3",
        y_q1_col="y_q1", y_q3_col="y_q3",
        x_label="Baseline LDDT silhouette-to-ref (median across models)",
        y_label="Guided LDDT at optimal gw (median across models)",
        title="LDDT Improvement Over Baseline Aggregated Across Models",
        draw_diagonal=True,
        equal_aspect=True,
        label_points=True,
    )
    if auto_save.value:
        save_fig(fig_lddt_agg, "rscc_scatterplots/lddt_scatter_aggregated_all_models.png")
    fig_lddt_agg
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Aggregated RSCC vs LDDT Scatter — collapsed across all models

    Each point is one protein; guided RSCC and guided LDDT silhouette-to-ref
    are taken at each model's optimal guidance weight and then collapsed across
    all models (median ± IQR across models). Same error-bar styling as the
    aggregated spiral figures.

    **Output:** `figures/rscc_scatterplots/rscc_vs_lddt_aggregated_all_models.svg`
    """)
    return


@app.cell
def _(auto_save, cross_metric_df, make_aggregated_scatter, save_fig):
    cross_agg_df = (
        cross_metric_df
        .groupby("protein", as_index=False)
        .agg(
            x_med=("rscc", "median"),
            x_q1=("rscc", lambda s: s.quantile(0.25)),
            x_q3=("rscc", lambda s: s.quantile(0.75)),
            y_med=("avg_silhouette_to_ref", "median"),
            y_q1=("avg_silhouette_to_ref", lambda s: s.quantile(0.25)),
            y_q3=("avg_silhouette_to_ref", lambda s: s.quantile(0.75)),
            n=("model_name", "count"),
        )
    )
    cross_agg_df = cross_agg_df[cross_agg_df["n"] >= 3]

    fig_cross_agg = make_aggregated_scatter(
        cross_agg_df,
        x_col="x_med", y_col="y_med",
        x_q1_col="x_q1", x_q3_col="x_q3",
        y_q1_col="y_q1", y_q3_col="y_q3",
        x_label="Guided RSCC at optimal gw (median across models)",
        y_label="Guided LDDT silhouette-to-ref (median across models)",
        title="RSCC vs LDDT at optimal guidance — aggregated across all models",
        draw_diagonal=False,
        equal_aspect=False,
        label_points=True,
        show_correlation=True,
    )
    if auto_save.value:
        save_fig(fig_cross_agg, "rscc_scatterplots/rscc_vs_lddt_aggregated_all_models.png")
    fig_cross_agg
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Figure 09 — LDDT-to-Reference Baseline Distribution

    Ridge plot of average silhouette scores (LDDT-based clustering relative to
    the deposited reference structure) at baseline. Higher values indicate the
    model's predicted ensemble is closer to the experimental structure.

    **Output:** `figures/09_lddt_to_ref_baseline_distribution.svg`
    """)
    return


@app.cell
def _(auto_save, lddt_results_1, make_ridge_plot, save_fig, smoothing):
    _base = lddt_results_1.query("guidance_strength == 'baseline' and selection_name != 'all altloc'")
    fig_09 = make_ridge_plot(
        _base.query("model_name not in @MODELS_ALT"),
        plot_column_name="avg_silhouette_to_ref",
        custom_xlabel="Average Silhouette Score to Reference (LDDT Clustering)",
        smoothing=smoothing.value,
        combined=True,
    )
    if auto_save.value:
        save_fig(fig_09, "09_lddt_to_ref_baseline_distribution.png")
    fig_09
    return


@app.cell
def _(
    MODELS_PROTENIX_GROUP,
    MODELS_RF3_GROUP,
    auto_save,
    lddt_results_1,
    make_ridge_plot,
    save_fig,
    smoothing,
):
    _base = lddt_results_1.query("guidance_strength == 'baseline' and selection_name != 'all altloc'")
    for _group, _tag in [(MODELS_PROTENIX_GROUP, "protenix_alts"), (MODELS_RF3_GROUP, "rf3_alts")]:
        _fig = make_ridge_plot(
            _base.query("model_name in @_group"),
            plot_column_name="avg_silhouette_to_ref",
            custom_xlabel="Average Silhouette Score to Reference (LDDT Clustering)",
            smoothing=smoothing.value,
            combined=True,
        )
        if auto_save.value:
            save_fig(_fig, f"09_lddt_to_ref_baseline_distribution_{_tag}.png")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Figure 10 — LDDT-to-Reference Across All Guidance Strengths

    Ridge plot of silhouette-to-reference scores across baseline, weak, and
    strong guidance for every model. Shows whether increasing guidance
    drives predictions closer to or further from the reference.

    **Output:** `figures/10_lddt_to_ref_all_guidance_distribution.svg`
    """)
    return


@app.cell
def _(auto_save, lddt_results_1, make_ridge_plot, save_fig, smoothing):
    fig_10 = make_ridge_plot(
        lddt_results_1.query("selection_name != 'all altloc' and model_name not in @MODELS_ALT"),
        "avg_silhouette_to_ref",
        smoothing=smoothing.value,
    )
    if auto_save.value:
        save_fig(fig_10, "10_lddt_to_ref_all_guidance_distribution.png")
    fig_10
    return


@app.cell
def _(
    MODELS_PROTENIX_GROUP,
    MODELS_RF3_GROUP,
    auto_save,
    lddt_results_1,
    make_ridge_plot,
    save_fig,
    smoothing,
):
    for _group, _tag in [(MODELS_PROTENIX_GROUP, "protenix_alts"), (MODELS_RF3_GROUP, "rf3_alts")]:
        _fig = make_ridge_plot(
            lddt_results_1.query("selection_name != 'all altloc' and model_name in @_group"),
            "avg_silhouette_to_ref",
            smoothing=smoothing.value,
        )
        if auto_save.value:
            save_fig(_fig, f"10_lddt_to_ref_all_guidance_distribution_{_tag}.png")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Figure 11 — Delta LDDT-to-Reference Distribution

    Ridge plot of the change in silhouette-to-reference score from baseline.
    Positive delta means guidance moved predictions closer to the experimental
    reference structure.

    **Output:** `figures/11_delta_lddt_to_ref_distribution.svg`
    """)
    return


@app.cell
def _(auto_save, lddt_merged, make_ridge_plot, save_fig):
    fig_11 = make_ridge_plot(
        lddt_merged.query("model_name not in @MODELS_ALT"),
        "delta_silhouette_to_ref", is_delta=True,
    )
    if auto_save.value:
        save_fig(fig_11, "11_delta_lddt_to_ref_distribution.png")
    fig_11
    return


@app.cell
def _(
    MODELS_PROTENIX_GROUP,
    MODELS_RF3_GROUP,
    auto_save,
    lddt_merged,
    make_ridge_plot,
    save_fig,
):
    for _group, _tag in [(MODELS_PROTENIX_GROUP, "protenix_alts"), (MODELS_RF3_GROUP, "rf3_alts")]:
        _fig = make_ridge_plot(
            lddt_merged.query("model_name in @_group"),
            "delta_silhouette_to_ref", is_delta=True,
        )
        if auto_save.value:
            save_fig(_fig, f"11_delta_lddt_to_ref_distribution_{_tag}.png")
    return


if __name__ == "__main__":
    app.run()
