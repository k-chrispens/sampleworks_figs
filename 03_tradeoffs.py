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
    import numpy as np
    import pandas as pd
    from loaders import (
        add_simplified_metadata_columns,
        get_rscc_nan_proteins,
        load_and_tag_path_csv,
        load_deposition_cutoffs,
        load_lddt,
        load_rscc_max_rmsd,
        load_training_cluster_counts,
        tag_path_dataframe,
    )
    from plot_tradeoff import (
        make_aggregated_spiral_figure,
        make_spiral_figure,
    )
    from plot_training_bias import (
        make_training_bias_calibration_explainer,
        make_training_bias_calibration_overlay,
        make_training_bias_calibration_overlay_legend,
        make_training_bias_no_cluster_strip,
        make_training_bias_one_sided_strip,
        make_training_bias_scatter,
        make_training_bias_small_multiples_v2,
    )
    from style import (
        MODEL_DISPLAY,
        OPTIMAL_GW,
        save_fig,
    )
    from transforms import (
        build_training_bias_df,
        filter_by_model_group,
        make_delta_frame,
        model_slug,
    )

    return (
        MODEL_DISPLAY,
        OPTIMAL_GW,
        add_simplified_metadata_columns,
        build_training_bias_df,
        filter_by_model_group,
        get_rscc_nan_proteins,
        load_and_tag_path_csv,
        load_deposition_cutoffs,
        load_lddt,
        load_rscc_max_rmsd,
        load_training_cluster_counts,
        make_aggregated_spiral_figure,
        make_delta_frame,
        make_spiral_figure,
        make_training_bias_calibration_explainer,
        make_training_bias_calibration_overlay,
        make_training_bias_calibration_overlay_legend,
        make_training_bias_no_cluster_strip,
        make_training_bias_one_sided_strip,
        make_training_bias_scatter,
        make_training_bias_small_multiples_v2,
        model_slug,
        np,
        pd,
        plt,
        save_fig,
        tag_path_dataframe,
    )


@app.cell
def _(MODEL_DISPLAY, mo):
    mo.md(r"""
    # Guidance-Weight Spiral Trajectories

    Each line traces a single protein through **(delta RSCC, geometry penalty)**
    space as guidance weight increases from baseline to the strongest setting.

    ## Figures produced by this notebook

    | ID | Description | Output |
    |----|-------------|--------|
    | 12-19 | Spiral trajectory per model | `{nn}_spiral_{model}.svg` |
    | 20 | Aggregated spiral (selected model group) | `20_aggregated_spiral_all[_variant].svg` |
    | 21-28 | Aggregated spiral per model | `{nn}_aggregated_spiral_{model}.svg` |
    | 29 | Training-cluster vs predicted state-B fraction -- 2×N small multiples | `29_training_bias_small_multiples[_variant].svg` |
    | 30 | Calibration overlay on low-training-cluster subset (n_cluster_evaluated ≤ 5) -- 1×N | `30_training_bias_calibration_overlay_low_cluster[_variant].svg` |
    | 31 | Calibration overlay on post-cutoff subset -- 1×N | `31_training_bias_calibration_overlay_post_cutoff[_variant].svg` |
    | 32 | Empty-training-cluster strip (`A==B==0`): half-violin + jittered strip of predicted B occupancy -- 1×N | `32_training_bias_no_cluster_strip[_variant].svg` |
    | 33 | One-sided training-cluster strip (`A>0 XOR B>0`): paired A/B swarms of predicted B occupancy per model -- 1×N | `33_training_bias_one_sided_strip[_variant].svg` |
    | 34 | Calibration overlay (Redesign A) -- 1×N | `34_training_bias_calibration_overlay[_variant].svg` |
    | 34L | Calibration overlay legend (standalone) | `34_training_bias_calibration_overlay_legend[_variant].svg` |
    | 35 | Calibration small multiples v2 (Redesign B) -- 2×N | `35_training_bias_small_multiples_v2[_variant].svg` |
    | 36 | Calibration explainer -- empty panel with y=x, y=0.5, and their bisector | `36_training_bias_calibration_explainer.svg` |
    """)

    _options = {
        v: k for k, v in MODEL_DISPLAY.items()
        if k in ("boltz2X-RAY", "boltz2MD", "protenix", "rf3",
                  "protenix-mini", "protenix-tiny", "rf3-chiral-off")
    }
    selected_model = mo.ui.dropdown(options=_options, value="Boltz2-XRD", label="Model")
    # Model-group selector, identical to notebook 02's. Drives the per-model-
    # colored multi-model figures (the Figure-20 overlay spiral and the
    # training-bias suite); the per-model spirals stay on `selected_model`.
    model_group = mo.ui.dropdown(
        options={"Main models": "main", "Protenix variants": "protenix", "RF3 variants": "rf3"},
        value="Main models",
        label="Model group",
    )
    auto_save = mo.ui.switch(label="Auto-save figures", value=False)
    save_all = mo.ui.switch(label="Save all models at once", value=False)
    exclude_nan_rscc = mo.ui.switch(label="Exclude NaN-RSCC proteins", value=False)
    mo.hstack([selected_model, model_group, auto_save, save_all, exclude_nan_rscc], justify="start", gap=1)
    return auto_save, exclude_nan_rscc, model_group, save_all, selected_model


@app.cell
def _(exclude_nan_rscc, get_rscc_nan_proteins):
    nan_rscc_proteins = (
        get_rscc_nan_proteins()
        if exclude_nan_rscc.value else set()
    )
    return (nan_rscc_proteins,)


@app.cell
def _(
    add_simplified_metadata_columns,
    load_and_tag_path_csv,
    load_rscc_max_rmsd,
    nan_rscc_proteins,
    pd,
    tag_path_dataframe,
):
    # Geometry metrics come per-model-group from results/. The clashscore files
    # carry the run path in `model_name` and the bond-violation files in `model`;
    # path extraction yields a canonical name for boltz_protenix / rf3_fixed but
    # collapses the variants to their raw run dir (protenix_X-RAY_DIFFRACTION for
    # BOTH mini and tiny, rf3_X-RAY_DIFFRACTION for chiral-off). So the three
    # variant files get an explicit model_name override AFTER tagging, mirroring
    # the RSCC/LDDT concats — without it mini+tiny silently merge into one model.
    _GEOM_SOURCES = [
        ("boltz_protenix", None),
        ("rf3_fixed", None),
        ("protenix_mini", "protenix-mini"),
        ("protenix_tiny", "protenix-tiny"),
        ("rf3_chiral_off", "rf3-chiral-off"),
    ]

    def _load_clash(_slug, _override):
        _c = load_and_tag_path_csv(
            f"results/clashscore_metrics_{_slug}.csv",
            filename_col="model_name",
            extra_cols={"selection": "all residues"},
        )
        if _override is not None:
            _c["model_name"] = _override
        return _c

    def _load_bond(_slug, _override):
        _a = pd.read_csv(f"results/bond_angle_violation_fractions_{_slug}.csv")
        _l = pd.read_csv(f"results/bond_length_violation_fractions_{_slug}.csv")
        _m = _a.merge(_l, on=["model", "model_n"], suffixes=["_bond_angles", "_bond_lengths"])
        _m = tag_path_dataframe(_m, filename_col="model")
        if _override is not None:
            _m["model_name"] = _override
        return _m

    clashes_1 = pd.concat(
        [_load_clash(_slug, _ov) for _slug, _ov in _GEOM_SOURCES], ignore_index=True
    )
    bond_stats = pd.concat(
        [_load_bond(_slug, _ov) for _slug, _ov in _GEOM_SOURCES], ignore_index=True
    )

    _rscc = load_rscc_max_rmsd()
    # One RSCC per (protein, model, gw) by averaging across per-altloc-region rows.
    rscc_all = (
        _rscc.groupby(["protein", "model_name", "guidance_weight"], as_index=False)["rscc"]
        .mean()
    )
    rscc_all = add_simplified_metadata_columns(rscc_all)
    if nan_rscc_proteins:
        clashes_1 = clashes_1[~clashes_1["protein"].isin(nan_rscc_proteins)]
        bond_stats = bond_stats[~bond_stats["protein"].isin(nan_rscc_proteins)]
        rscc_all = rscc_all[~rscc_all["protein"].isin(nan_rscc_proteins)]
    return bond_stats, clashes_1, rscc_all


@app.cell
def _(bond_stats, clashes_1, make_delta_frame, rscc_all):
    _group_cols = ["protein", "model_name", "guidance_weight", "guidance_strength"]

    _clash_stats = clashes_1.groupby(_group_cols, as_index=False).agg(
        clashscore=("clashscore", "mean")
    )
    _bond_stats_agg = bond_stats.groupby(_group_cols, as_index=False).agg(
        outlier_fraction_bond_angles=("outlier_fraction_bond_angles", "mean"),
        outlier_fraction_bond_lengths=("outlier_fraction_bond_lengths", "mean"),
    )
    _rscc_stats = rscc_all.groupby(_group_cols, as_index=False).agg(
        rscc=("rscc", "mean")
    )

    _clash_delta = make_delta_frame(_clash_stats, "clashscore")[_group_cols + ["delta_clashscore"]]
    _bond_angle_delta = make_delta_frame(_bond_stats_agg, "outlier_fraction_bond_angles")[
        _group_cols + ["delta_outlier_fraction_bond_angles"]
    ]
    _bond_length_delta = make_delta_frame(_bond_stats_agg, "outlier_fraction_bond_lengths")[
        _group_cols + ["delta_outlier_fraction_bond_lengths"]
    ]
    _rscc_delta = make_delta_frame(_rscc_stats, "rscc")[_group_cols + ["delta_rscc"]]

    _merge_on = _group_cols
    tradeoff_df = (
        _clash_delta
        .merge(_bond_angle_delta, on=_merge_on)
        .merge(_bond_length_delta, on=_merge_on)
        .merge(_rscc_delta, on=_merge_on, how="left")
    )
    tradeoff_df["delta_outlier_fraction_bond_lengths_pct"] = (
        100.0 * tradeoff_df["delta_outlier_fraction_bond_lengths"]
    )
    tradeoff_df["delta_outlier_fraction_bond_angles_pct"] = (
        100.0 * tradeoff_df["delta_outlier_fraction_bond_angles"]
    )
    tradeoff_df["geometry_penalty"] = tradeoff_df[[
        "delta_clashscore",
        "delta_outlier_fraction_bond_lengths_pct",
        "delta_outlier_fraction_bond_angles_pct",
    ]].mean(axis=1)
    return (tradeoff_df,)


@app.cell
def _(
    auto_save,
    make_spiral_figure,
    mo,
    model_slug,
    np,
    plt,
    save_all,
    save_fig,
    selected_model,
    tradeoff_df,
):
    _spiral_nums = {
        "boltz2X-RAY": 12, "boltz2MD": 13, "protenix": 14, "rf3": 15,
        "protenix-mini": 16, "protenix-tiny": 17, "rf3-chiral-off": 19,
    }
    _all_models = [
        m for m in list(_spiral_nums.keys())
        if m in tradeoff_df["model_name"].unique()
    ]

    # Shared axis limits across all models
    _mx = tradeoff_df["delta_rscc"].replace([np.inf, -np.inf], np.nan).dropna()
    _my = tradeoff_df["geometry_penalty"].replace([np.inf, -np.inf], np.nan).dropna()
    _x_lo, _x_hi = min(0.0, float(_mx.min())), max(0.0, float(_mx.max()))
    _y_lo, _y_hi = min(0.0, float(_my.min())), max(0.0, float(_my.max()))
    _x_pad = max(0.08 * (_x_hi - _x_lo), 0.02)
    _y_pad = max(0.08 * (_y_hi - _y_lo), 0.50)
    _shared_xlim = (_x_lo - _x_pad, _x_hi + _x_pad)
    _shared_ylim = (_y_lo - _y_pad, _y_hi + _y_pad)

    if save_all.value and auto_save.value:
        for _model in _all_models:
            _model_df = tradeoff_df[tradeoff_df["model_name"] == _model].copy()
            if _model_df.empty:
                continue
            _fig = make_spiral_figure(_model_df, _model, xlim=_shared_xlim, ylim=_shared_ylim)
            save_fig(_fig, f"{_spiral_nums[_model]:02d}_spiral_{model_slug(_model)}.png")
            plt.close(_fig)

    _cur_model = selected_model.value
    _spiral_fig = None
    if _cur_model in tradeoff_df["model_name"].unique():
        _model_df = tradeoff_df[tradeoff_df["model_name"] == _cur_model].copy()
        _spiral_fig = make_spiral_figure(_model_df, _cur_model, xlim=_shared_xlim, ylim=_shared_ylim)
        if auto_save.value and not save_all.value:
            save_fig(_spiral_fig, f"{_spiral_nums[_cur_model]:02d}_spiral_{model_slug(_cur_model)}.png")

    mo.vstack([
        mo.md("""## Figures 12-19 -- Spiral Trajectory per Model
        Each protein traces through (delta RSCC, geometry penalty) space.
        `figures/12_spiral_boltz2_xrd.svg` through `figures/19_spiral_rf3_chiral_off.svg`"""),
        _spiral_fig,
    ])
    return


@app.cell
def _(
    auto_save,
    filter_by_model_group,
    make_aggregated_spiral_figure,
    mo,
    model_group,
    pd,
    save_fig,
    tradeoff_df,
):
    # Figure 20 honors the model-group dropdown. Main → the four main models;
    # Protenix / RF3 variants → that family's members and colors. Non-main saves
    # get the _protenix_alts / _rf3_alts suffix.
    _grp_df, _suffix = filter_by_model_group(tradeoff_df, model_group.value)
    _canonical = [
        "boltz2X-RAY", "boltz2MD", "protenix", "rf3",
        "protenix-mini", "protenix-tiny", "rf3-chiral-off",
    ]
    _present = set(_grp_df["model_name"].unique())
    _model_order = [m for m in _canonical if m in _present]

    # Geometry-coverage guard (report, don't fix): the spiral drops any
    # (model, gw) cell with < 3 proteins, so a variant whose geometry sweep is
    # sparser than its family's main model would draw a degenerate trajectory.
    # Surface that to the user rather than silently shipping it; the underlying
    # data gap is resolved upstream, not in the notebook.
    if model_group.value != "main":
        _min_n = 3
        _surv = (
            _grp_df.groupby(["model_name", "guidance_weight"]).size()
            .reset_index(name="n")
        )
        _surv = _surv[_surv["n"] >= _min_n]
        _gws = {m: sorted(g["guidance_weight"].unique()) for m, g in _surv.groupby("model_name")}
        print(f"Figure 20 geometry coverage [{model_group.value}] — gw points with n>={_min_n}:")
        for _m in _model_order:
            print(f"  {_m}: {_gws.get(_m, [])}")
        _family_main = "protenix" if model_group.value == "protenix" else "rf3"
        _main_gws = set(_gws.get(_family_main, []))
        for _m in _model_order:
            if _m == _family_main:
                continue
            _missing = _main_gws - set(_gws.get(_m, []))
            if _missing:
                print(
                    f"⚠️  geometry coverage: variant {_m} is missing gw(s) "
                    f"{sorted(_missing)} that {_family_main} has (n>={_min_n} after "
                    f"aggregation). Figure 20's {_m} arm is sparser than the main "
                    f"model's — a data gap to resolve upstream, not in the notebook."
                )

    _fig = make_aggregated_spiral_figure(
        tradeoff_df,
        model_order=_model_order,
        ylim=[-1., 1000.],
        xlim=[-0.05, 0.25],
    )
    if auto_save.value:
        save_fig(_fig, f"20_aggregated_spiral_all{_suffix}.svg")
    mo.vstack([
        mo.md("""## Figure 20 -- Aggregated Spiral (selected model group)
        Median delta RSCC and geometry penalty per guidance weight, with IQR error bars.
        Honors the model-group dropdown (Main / Protenix variants / RF3 variants);
        non-main groups save with the `_protenix_alts` / `_rf3_alts` suffix.
        `figures/20_aggregated_spiral_all[_variant].svg`"""),
        _fig,
    ])
    return


@app.cell
def _(
    auto_save,
    make_aggregated_spiral_figure,
    mo,
    model_slug,
    plt,
    save_all,
    save_fig,
    selected_model,
    tradeoff_df,
):
    _agg_spiral_nums = {
        "boltz2X-RAY": 21, "boltz2MD": 22, "protenix": 23, "rf3": 24,
        "protenix-mini": 25, "protenix-tiny": 26, "rf3-chiral-off": 28,
    }
    _all_models_agg = [
        m for m in list(_agg_spiral_nums.keys())
        if m in tradeoff_df["model_name"].unique()
    ]

    if save_all.value and auto_save.value:
        for _model in _all_models_agg:
            _fig = make_aggregated_spiral_figure(tradeoff_df, model_key=_model)
            save_fig(_fig, f"{_agg_spiral_nums[_model]:02d}_aggregated_spiral_{model_slug(_model)}.svg")
            plt.close(_fig)

    _cur_model = selected_model.value
    _agg_fig = None
    if _cur_model in tradeoff_df["model_name"].unique():
        _agg_fig = make_aggregated_spiral_figure(tradeoff_df, model_key=_cur_model)
        if auto_save.value and not save_all.value:
            save_fig(_agg_fig, f"{_agg_spiral_nums[_cur_model]:02d}_aggregated_spiral_{model_slug(_cur_model)}.svg")

    mo.vstack([
        mo.md("""## Figures 21-28 -- Aggregated Spiral per Model
        Per-model aggregated trajectories with median + IQR.
        `figures/21_aggregated_spiral_boltz2_xrd.svg` through `figures/28_aggregated_spiral_rf3_chiral_off.svg`"""),
        _agg_fig,
    ])
    return


@app.cell
def _(
    filter_by_model_group,
    load_deposition_cutoffs,
    load_lddt,
    load_training_cluster_counts,
    model_group,
    pd,
):
    cluster_df = load_training_cluster_counts()
    cutoffs_df = load_deposition_cutoffs()
    # LDDT concatenated from the five per-model-group CSVs in results/ (same
    # inline concat + variant overrides as notebook 02's LDDT cell).
    lddt_all = pd.concat(
        [
            load_lddt("results/lddt_results_boltz_protenix.csv"),
            load_lddt("results/lddt_results_rf3_fixed.csv"),
            load_lddt("results/lddt_results_protenix_mini.csv", model_name_override="protenix-mini"),
            load_lddt("results/lddt_results_protenix_tiny.csv", model_name_override="protenix-tiny"),
            load_lddt("results/lddt_results_rf3_chiral_off.csv", model_name_override="rf3-chiral-off"),
        ],
        ignore_index=True,
    )
    # The training-bias suite honors the model-group dropdown. Resolve the
    # group's members (canonical order) and save suffix once here; every
    # training-bias figure cell consumes bias_model_order / bias_suffix.
    lddt_for_bias, bias_suffix = filter_by_model_group(lddt_all, model_group.value)
    _canonical = [
        "boltz2X-RAY", "boltz2MD", "protenix", "rf3",
        "protenix-mini", "protenix-tiny", "rf3-chiral-off",
    ]
    bias_model_order = [
        m for m in _canonical if m in set(lddt_for_bias["model_name"].unique())
    ]
    return bias_model_order, bias_suffix, cluster_df, cutoffs_df, lddt_for_bias


@app.cell
def _(
    OPTIMAL_GW,
    bias_model_order,
    build_training_bias_df,
    cluster_df,
    cutoffs_df,
    lddt_for_bias,
):
    # Models follow the selected group (bias_model_order); OPTIMAL_GW covers the
    # variants, so the baseline (gw=0) and optimal frames build for whichever
    # family is selected.
    bias_baseline_df = build_training_bias_df(
        lddt_for_bias, cluster_df, cutoffs_df,
        {m: 0.0 for m in bias_model_order},
    )
    bias_optimal_df = build_training_bias_df(
        lddt_for_bias, cluster_df, cutoffs_df,
        {m: OPTIMAL_GW[m] for m in bias_model_order},
    )
    return bias_baseline_df, bias_optimal_df


@app.cell
def _(
    auto_save,
    bias_model_order,
    bias_suffix,
    bias_baseline_df,
    bias_optimal_df,
    make_training_bias_scatter,
    mo,
    save_fig,
):
    _fig = make_training_bias_scatter(
        bias_baseline_df, bias_optimal_df, model_order=bias_model_order,
    )
    if auto_save.value:
        save_fig(_fig, f"29_training_bias_small_multiples{bias_suffix}.svg")
    mo.vstack([
        mo.md("""## Figure 29 -- Training-Cluster vs Predicted State-B Fraction
        2-row small multiples: columns = the selected model group's members,
        rows = condition (Baseline, Optimal guidance). One point per case with
        cluster sample size s > 0; point alpha scales with log(s + 1) normalized
        across the full dataset. `y = 0.5` reference is drawn behind the points.
        Honors the model-group dropdown; non-main groups save with the
        `_protenix_alts` / `_rf3_alts` suffix.
        `figures/29_training_bias_small_multiples[_variant].svg`"""),
        _fig,
    ])
    return


@app.cell
def _(
    auto_save,
    bias_model_order,
    bias_suffix,
    bias_baseline_df,
    bias_optimal_df,
    make_training_bias_calibration_overlay,
    mo,
    save_fig,
):
    _fig = make_training_bias_calibration_overlay(
        bias_baseline_df, bias_optimal_df,
        model_order=bias_model_order,
        annotation_layout="split",
    )
    if auto_save.value:
        save_fig(_fig, f"34_training_bias_calibration_overlay{bias_suffix}.svg")
    mo.vstack([
        mo.md("""## Figure 34 -- Redesign A: Calibration Overlay (1×N)
        Single row of panels, one per model in the selected group. Baseline and Optimal guidance
        overlaid in the same axes, with the y = x diagonal (training prior), y = 0.5
        (guidance target), binned conditional means for each condition, and a
        marginal y-histogram on the right edge of each panel. Guidance slope and
        Pearson r annotated per panel. Honors the model-group dropdown; non-main
        groups save with the `_protenix_alts` / `_rf3_alts` suffix.
        `figures/34_training_bias_calibration_overlay[_variant].svg`"""),
        _fig,
    ])
    return


@app.cell
def _(
    auto_save,
    bias_model_order,
    bias_suffix,
    make_training_bias_calibration_overlay_legend,
    mo,
    save_fig,
):
    _fig = make_training_bias_calibration_overlay_legend(model_order=bias_model_order)
    if auto_save.value:
        save_fig(_fig, f"34_training_bias_calibration_overlay_legend{bias_suffix}.svg")
    mo.vstack([
        mo.md("""## Figure 34L -- Calibration Overlay Legend (standalone)
        Standalone legend for Figure 34, separated so it can be arranged
        independently of the plot panels in the final layout. Shows the selected
        group's optimal-guidance markers; non-main groups save with the
        `_protenix_alts` / `_rf3_alts` suffix.
        `figures/34_training_bias_calibration_overlay_legend[_variant].svg`"""),
        _fig,
    ])
    return


@app.cell
def _(
    auto_save,
    bias_model_order,
    bias_suffix,
    bias_baseline_df,
    bias_optimal_df,
    make_training_bias_small_multiples_v2,
    mo,
    save_fig,
):
    _fig = make_training_bias_small_multiples_v2(
        bias_baseline_df, bias_optimal_df, model_order=bias_model_order,
    )
    if auto_save.value:
        save_fig(_fig, f"35_training_bias_small_multiples_v2{bias_suffix}.svg")
    mo.vstack([
        mo.md("""## Figure 35 -- Redesign B: Calibration Small Multiples v2 (2×N)
        2-row grid mirroring Figure 29's geometry (rows = condition, columns = model).
        Each panel adds: both reference lines (y = x training prior, y = 0.5 guidance
        target), a binned conditional-mean curve, an OLS slope and Pearson r
        annotation, and a single-condition marginal y-histogram strip on the right
        edge. Inset legend in the top-left panel names the curve and reference lines.
        Honors the model-group dropdown; non-main groups save with the
        `_protenix_alts` / `_rf3_alts` suffix.
        `figures/35_training_bias_small_multiples_v2[_variant].svg`"""),
        _fig,
    ])
    return


@app.cell
def _(
    auto_save,
    make_training_bias_calibration_explainer,
    mo,
    save_fig,
):
    _fig = make_training_bias_calibration_explainer()
    if auto_save.value:
        save_fig(_fig, "36_training_bias_calibration_explainer.svg")
    mo.vstack([
        mo.md("""## Figure 36 -- Calibration Explainer (Empty Panel)
        Intentionally minimal companion to Figure 35. Same axes as the v2 small
        multiples panels: `y = x` (training prior, dashed) and `y = 0.5` (guidance
        target, dashed) reference lines, plus a solid bisector at `y = 0.5*x + 0.25`
        passing through their intersection at `(0.5, 0.5)`. No data, no annotations
        -- intended for downstream labeling.
        `figures/36_training_bias_calibration_explainer.svg`"""),
        _fig,
    ])
    return


@app.cell
def _(
    auto_save,
    bias_model_order,
    bias_suffix,
    bias_baseline_df,
    bias_optimal_df,
    make_training_bias_calibration_overlay,
    mo,
    save_fig,
):
    # F30/F31 intentionally keep the default stacked annotation layout — the
    # split-corner layout is only for the main F34 cell.
    _fig_low = make_training_bias_calibration_overlay(
        bias_baseline_df, bias_optimal_df,
        model_order=bias_model_order,
        max_cluster_size=5,
        n_bins=3, min_bin_count=2,
    )
    _fig_post = make_training_bias_calibration_overlay(
        bias_baseline_df, bias_optimal_df,
        model_order=bias_model_order, cutoff_status="Post",
    )
    if auto_save.value:
        save_fig(_fig_low, f"30_training_bias_calibration_overlay_low_cluster{bias_suffix}.svg")
        save_fig(_fig_post, f"31_training_bias_calibration_overlay_post_cutoff{bias_suffix}.svg")

    mo.vstack([
        mo.md("""## Figures 30, 31 -- Calibration Overlay on OOD Subsets
        Figure-34 layout (1x4, baseline + optimal overlaid with binned mean curves,
        `y = x` and `y = 0.5` references, slope/r annotations, marginal y-histogram per
        panel) restricted to two out-of-distribution subsets.

        - Figure 30: selections with `n_cluster_evaluated ≤ 5` (sparse training prior).
        - Figure 31: `cutoff_status == "Post"` (deposited after the training cutoff).

        Honors the model-group dropdown; non-main groups save with the
        `_protenix_alts` / `_rf3_alts` suffix.
        `figures/30_training_bias_calibration_overlay_low_cluster[_variant].svg`
        `figures/31_training_bias_calibration_overlay_post_cutoff[_variant].svg`"""),
        _fig_low,
        _fig_post,
    ])
    return


@app.cell
def _(
    auto_save,
    bias_model_order,
    bias_suffix,
    bias_baseline_df,
    bias_optimal_df,
    make_training_bias_no_cluster_strip,
    mo,
    save_fig,
):
    _fig = make_training_bias_no_cluster_strip(
        bias_baseline_df, bias_optimal_df, model_order=bias_model_order,
    )
    if auto_save.value:
        save_fig(_fig, f"32_training_bias_no_cluster_strip{bias_suffix}.svg")
    mo.vstack([
        mo.md("""## Figure 32 -- Empty-Training-Cluster Strip
        1×N vertical strip (one panel per model in the selected group) for selections where the
        pre-cutoff training cluster has zero structures assigned to either
        altloc by lDDT comparison (`n_assigned_to_A == 0` AND
        `n_assigned_to_B == 0`). These rows have `frac_b_cluster`
        undefined and cannot appear on Figure 34, so they are shown here
        instead: each panel renders a half-violin (KDE) on the left and a
        jittered strip on the right, overlaying Baseline and Optimal-guidance.
        `y = 0.5` reference line; per-panel `n` annotated. No marginal histogram.
        Honors the model-group dropdown; non-main groups save with the
        `_protenix_alts` / `_rf3_alts` suffix.
        `figures/32_training_bias_no_cluster_strip[_variant].svg`"""),
        _fig,
    ])
    return


@app.cell
def _(
    auto_save,
    bias_model_order,
    bias_suffix,
    bias_baseline_df,
    bias_optimal_df,
    make_training_bias_one_sided_strip,
    mo,
    save_fig,
):
    _fig = make_training_bias_one_sided_strip(
        bias_baseline_df, bias_optimal_df, model_order=bias_model_order,
    )
    if auto_save.value:
        save_fig(_fig, f"33_training_bias_one_sided_strip{bias_suffix}.svg")
    mo.vstack([
        mo.md("""## Figure 33 -- One-sided Training-Cluster Strip
        1×N paired-swarm strip (one panel per model in the selected group) for selections where
        the pre-cutoff training cluster is strictly one-sided -- either
        `n_assigned_to_A > 0 AND n_assigned_to_B == 0` (only-A, left strip) or
        `n_assigned_to_A == 0 AND n_assigned_to_B > 0` (only-B, right strip).
        Equivalent to `frac_b_cluster in {0, 1}`. These rows pile up at
        `x in {0, 1}` on Figure 34; this companion breaks them out by altloc.
        Baseline (gray, alpha 0.7) drawn first; Optimal guidance (model color,
        alpha 0.7) drawn on top, so points alpha-blend across conditions but
        do not overlap within a condition. Horizontal mean tick per condition.
        `y = 0.5` reference line; per-strip `n` annotated. Honors the model-group
        dropdown; non-main groups save with the `_protenix_alts` / `_rf3_alts` suffix.
        `figures/33_training_bias_one_sided_strip[_variant].svg`"""),
        _fig,
    ])
    return


if __name__ == "__main__":
    app.run()
