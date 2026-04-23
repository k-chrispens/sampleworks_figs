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
    from analysis_helpers import (
        add_simplified_metadata_columns,
        load_and_tag_path_csv,
        tag_path_dataframe,
        load_rscc,
        make_delta_frame,
        make_spiral_figure,
        make_aggregated_spiral_figure,
        model_slug,
        save_fig,
        MODEL_DISPLAY,
    )

    return (
        MODEL_DISPLAY,
        add_simplified_metadata_columns,
        load_and_tag_path_csv,
        load_rscc,
        make_aggregated_spiral_figure,
        make_delta_frame,
        make_spiral_figure,
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
    | 20 | Aggregated spiral (all models) | `20_aggregated_spiral_all.svg` |
    | 21-28 | Aggregated spiral per model | `{nn}_aggregated_spiral_{model}.svg` |
    """)

    _options = {
        v: k for k, v in MODEL_DISPLAY.items()
        if k in ("boltz2X-RAY", "boltz2MD", "protenix", "rf3",
                  "protenix-mini", "protenix-tiny", "rf3-old", "rf3-chiral-off")
    }
    selected_model = mo.ui.dropdown(options=_options, value="Boltz2-XRD", label="Model")
    auto_save = mo.ui.switch(label="Auto-save figures", value=False)
    save_all = mo.ui.switch(label="Save all models at once", value=False)
    mo.hstack([selected_model, auto_save, save_all], justify="start", gap=1)
    return auto_save, save_all, selected_model


@app.cell
def _(
    add_simplified_metadata_columns,
    load_and_tag_path_csv,
    load_rscc,
    pd,
    tag_path_dataframe,
):
    clashes_1 = load_and_tag_path_csv(
        "initial40results/clashscore_metrics.csv",
        filename_col="model_name",
        extra_cols={"selection": "all residues"},
    )
    _ba = pd.read_csv("initial40results/bond_angle_violation_fractions.csv")
    _bl = pd.read_csv("initial40results/bond_length_violation_fractions.csv")
    bond_stats = _ba.merge(_bl, on=["model", "model_n"], suffixes=["_bond_angles", "_bond_lengths"])
    bond_stats = tag_path_dataframe(bond_stats, filename_col="model")

    _rscc = load_rscc("initial40results/rscc_results.csv")
    rscc_all = _rscc.groupby(["protein", "model_name", "guidance_weight"]).last().reset_index()
    rscc_all = add_simplified_metadata_columns(rscc_all)
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
        "protenix-mini": 16, "protenix-tiny": 17, "rf3-old": 18, "rf3-chiral-off": 19,
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
def _(auto_save, make_aggregated_spiral_figure, mo, save_fig, tradeoff_df):
    _fig = make_aggregated_spiral_figure(
        tradeoff_df,
        model_order=["boltz2X-RAY", "boltz2MD", "protenix", "rf3"],
        ylim=[-1., 1000.],
        xlim=[-0.05, 0.2],
    )
    if auto_save.value:
        save_fig(_fig, "20_aggregated_spiral_all.svg")
    mo.vstack([
        mo.md("""## Figure 20 -- Aggregated Spiral (All Models)
        Median delta RSCC and geometry penalty per guidance weight, with IQR error bars.
        `figures/20_aggregated_spiral_all.svg`"""),
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
        "protenix-mini": 25, "protenix-tiny": 26, "rf3-old": 27, "rf3-chiral-off": 28,
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


if __name__ == "__main__":
    app.run()
