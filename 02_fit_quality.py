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
        merge_classification_by_containment,
    )
    from plot_classified import make_by_class_combined_grid
    from plot_distributions import make_combined_baseline_distributions
    from plot_scatter import (
        make_aggregated_scatter,
        make_min_rmsd_scatter,
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
        load_min_altloc_rmsd,
        load_rscc,
        merge_classification_by_containment,
        make_aggregated_scatter,
        make_by_class_combined_grid,
        make_combined_baseline_distributions,
        make_min_rmsd_scatter,
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
    | -- | Aggregated min-RMSD scatter (A vs B) | `rmsd_scatterplots/min_rmsd_scatter_aggregated_all_models.svg` |
    | -- | Per-model min-RMSD scatter (A vs B) | `rmsd_scatterplots/min_rmsd_scatter_{model}.svg` |
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
def _(load_min_altloc_rmsd, nan_rscc_proteins, pd):
    # Concatenate the four per-source min-RMSD CSVs into one tidy frame.
    # protenix-mini / protenix-tiny CSVs both store ``model == "protenix"``
    # so we override model_name at load time to disambiguate.
    _pieces = [
        load_min_altloc_rmsd("min_altloc_rmsd_results_boltz_protenix.csv"),
        load_min_altloc_rmsd(
            "min_altloc_rmsd_results_protenix_mini.csv",
            model_name_override="protenix-mini",
        ),
        load_min_altloc_rmsd(
            "min_altloc_rmsd_results_protenix_tiny.csv",
            model_name_override="protenix-tiny",
        ),
        load_min_altloc_rmsd("min_altloc_rmsd_results_rf3_fixed.csv"),
    ]
    min_rmsd_all = pd.concat(_pieces, ignore_index=True)
    if nan_rscc_proteins:
        min_rmsd_all = min_rmsd_all[~min_rmsd_all["protein"].isin(nan_rscc_proteins)]
    return (min_rmsd_all,)


@app.cell
def _(OPTIMAL_GW, format_residue_range, min_rmsd_all, pd):
    # Per-(protein, selection, model) min RMSD at each model's optimal gw.
    _gw_df = pd.DataFrame(list(OPTIMAL_GW.items()), columns=["model_name", "optimal_gw"])
    _region = min_rmsd_all[min_rmsd_all["selection_name"] != "all altloc"]

    # Audit: which (protein, selection) pairs exist for each model in any gw
    # but lack a row at that model's optimal gw?
    _all_keys_per_model = {
        m: set(zip(g["protein"], g["selection"]))
        for m, g in _region.groupby("model_name")
    }

    min_rmsd_per_model_df = (
        _region.merge(_gw_df, on="model_name")
        .query("guidance_weight == optimal_gw")
        [["protein", "selection", "model_name", "guidance_weight",
          "min_rmsd_to_A", "min_rmsd_to_B"]]
        .reset_index(drop=True)
    )
    min_rmsd_per_model_df["label"] = min_rmsd_per_model_df.apply(
        lambda r: f"{str(r['protein']).upper()} {format_residue_range(r['selection'])}",
        axis=1,
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

    # Cross-figure highlight groups for the min-RMSD scatter family.
    # Colors match the RSCC scatter convention in plot_scatter.py:
    #   #F57C00 orange = far-from-ceiling family (2QCV 2-4, 4OLE 424-426),
    #   #1565C0 blue   = manually-curated biological highlight (5SOP 131-133),
    #   #D32F2F red    = decreased-RSCC family (curated 7Z0E altloc segments
    #                    flagged in the RSCC scatter; intentionally not the
    #                    full 7Z0E set to keep the panel readable).
    MIN_RMSD_EXTRA_HIGHLIGHTS = {
        "#F57C00": ["2QCV 2-4", "4OLE 424-426"],
        "#1565C0": ["5SOP 131-133"],
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

    _fig_min_rmsd_agg = make_min_rmsd_scatter(
        _min_rmsd_agg_df,
        title="Minimum RMSD to Each Altloc",
        extra_highlight_groups=MIN_RMSD_EXTRA_HIGHLIGHTS,
    )
    if auto_save.value:
        save_fig(
            _fig_min_rmsd_agg,
            "rmsd_scatterplots/min_rmsd_scatter_aggregated_all_models.svg",
        )
    mo.vstack([
        mo.md("""## Aggregated Min-RMSD Scatter -- Min RMSD to Altloc A vs B
        Each point is one (protein, selection); X / Y are the median across main models of the minimum RMSD across the 8-member ensemble at each model's optimal guidance weight. Pinned highlights ring `6B8X 181-183` and `5MHX 162-164` in blue (matching the aggregated RSCC scatter). Extra RSCC-matched highlight rings (per-model RSCC scatter colors): orange (`#F57C00`) `2QCV 2-4` and `4OLE 424-426`, blue (`#1565C0`) `5SOP 131-133`, red (`#D32F2F`) a curated `7Z0E` subset (`11-13`, `18-20`, `54-56`, `82-84`, `122-124`, `181-183`, `206-208`, `214-216`).
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
        Per-model scatter at each model's optimal guidance weight (no baseline overlay). Models without a min-RMSD CSV (e.g. `rf3-chiral-off`) render a warning callout in place of the figure.
        `figures/rmsd_scatterplots/min_rmsd_scatter_{model}.svg`"""),
        _content,
    ])
    return


@app.cell
def _(
    OPTIMAL_GW,
    auto_save,
    classification_df,
    load_rscc,
    make_by_class_combined_grid,
    merge_classification_by_containment,
    mo,
    nan_rscc_proteins,
    pd,
    rscc,
    save_fig,
):
    _gw_df = pd.DataFrame(list(OPTIMAL_GW.items()), columns=["model_name", "optimal_gw"])

    # Small-segment classes (per-altloc-region RSCC from rscc_results.csv).
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

    # Domain-shift class uses whole-domain RSCC from the extracted clean CSV.
    _rscc_ds = load_rscc("initial40results/rscc_results_domain_shift.csv")
    if nan_rscc_proteins:
        _rscc_ds = _rscc_ds[~_rscc_ds["protein"].isin(nan_rscc_proteins)]
    _ds_cls = classification_df.query("classification == 'domain_shift'")[
        ["protein", "selection", "classification"]
    ]
    _ds_region = _rscc_ds.merge(_ds_cls, on=["protein", "selection"], how="inner")

    _baseline_ds = (
        _ds_region.query("guidance_weight == 0.0")
        [["protein", "selection", "model_name", "rscc", "classification"]]
        .rename(columns={"rscc": "rscc_baseline"})
    )
    _optimal_ds = (
        _ds_region.merge(_gw_df, on="model_name")
        .query("guidance_weight == optimal_gw")
        [["protein", "selection", "model_name", "rscc"]]
        .rename(columns={"rscc": "rscc_optimal"})
    )
    _rscc_ds_df = _baseline_ds.merge(_optimal_ds, on=["protein", "selection", "model_name"])

    _merged = pd.concat(
        [_rscc_by_class_df, _rscc_ds_df],
        ignore_index=True,
    ).query("model_name not in @MODELS_ALT")

    _fig_combined = make_by_class_combined_grid(_merged)
    if auto_save.value:
        save_fig(_fig_combined, "rscc_scatterplots/rscc_by_class_combined.svg")

    mo.vstack([
        mo.md("""## RSCC by Conformational Type — Combined Grid (4 models × 4 classes)
        Single combined figure: rows = main models, columns = conformational-change
        classes. Per-panel gray shadow shows that row's model pooled cloud.
        Domain-shift column reads from
        `initial40results/rscc_results_domain_shift.csv` (whole-domain RSCC
        extracted from the legacy `.bak`); the other three columns use the
        per-altloc-region RSCC from `rscc_results.csv`.
        `figures/rscc_scatterplots/rscc_by_class_combined.svg`"""),
        _fig_combined,
    ])
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
