import marimo

__generated_with = "0.22.4"
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
        load_lddt,
        make_delta_frame,
        make_spiral_figure,
        make_aggregated_spiral_figure,
        save_fig,
        MODEL_DISPLAY,
        MODELS_ALT,
    )

    return (
        MODEL_DISPLAY,
        add_simplified_metadata_columns,
        load_and_tag_path_csv,
        load_rscc,
        make_aggregated_spiral_figure,
        make_delta_frame,
        make_spiral_figure,
        np,
        pd,
        plt,
        save_fig,
        tag_path_dataframe,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Guidance-Weight Spiral Trajectories

    This notebook visualizes how each protein's fit quality and geometry evolve
    as guidance weight increases, producing one spiral trajectory plot per model.

    Each line in a spiral plot traces a single protein through
    **(delta RSCC, geometry penalty)** space as guidance weight increases from
    baseline (the origin) to the strongest setting tested. Dots are colored
    from cool to warm with increasing guidance weight.

    **Common trajectory shapes:**

    - **Rightward sweep** — RSCC improves with tolerable geometry cost (correctable).
    - **Upward hook / spiral** — geometry degrades faster than RSCC improves at strong guidance.
    - **Stuck at origin** — neither metric moves much regardless of guidance.

    The **geometry penalty** is the unweighted mean of three baseline-relative terms:
    delta clashscore, delta bond-length outlier fraction (pp), and delta bond-angle
    outlier fraction (pp).
    """)
    return


@app.cell
def _(mo):
    selected_model = mo.ui.dropdown(
        options={
            "Boltz2 X-RAY": "boltz2X-RAY",
            "Boltz2 MD": "boltz2MD",
            "Protenix": "protenix",
            "RF3": "rf3",
            "Protenix Mini": "protenix-mini",
            "Protenix Tiny": "protenix-tiny",
            "RF3 Old": "rf3-old",
            "RF3 Chiral Off": "rf3-chiral-off",
        },
        value="Boltz2 X-RAY",
        label="Model",
    )
    auto_save = mo.ui.switch(label="Auto-save figures", value=False)
    save_all = mo.ui.switch(label="Save all models at once", value=False)
    mo.hstack([selected_model, auto_save, save_all], justify="start", gap=1)
    return auto_save, save_all, selected_model


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Data Loading

    We load clash scores, bond violation fractions, RSCC, and LDDT metrics,
    then compute baseline-relative deltas and merge them into a single
    tradeoff dataframe used to build spiral trajectories.
    """)
    return


@app.cell
def _(load_and_tag_path_csv):
    clashes_1 = load_and_tag_path_csv(
        "initial40results/clashscore_metrics.csv",
        filename_col="model_name",
        extra_cols={"selection": "all residues"},
    )
    return (clashes_1,)


@app.cell
def _(pd, tag_path_dataframe):
    _ba = pd.read_csv("initial40results/bond_angle_violation_fractions.csv")
    _bl = pd.read_csv("initial40results/bond_length_violation_fractions.csv")
    bond_stats = _ba.merge(_bl, on=["model", "model_n"], suffixes=["_bond_angles", "_bond_lengths"])
    bond_stats = tag_path_dataframe(bond_stats, filename_col="model")
    return (bond_stats,)


@app.cell
def _(add_simplified_metadata_columns, load_rscc):
    _rscc = load_rscc("initial40results/rscc_results.csv")
    rscc_all = _rscc.groupby(["protein", "model_name", "guidance_weight"]).last().reset_index()
    rscc_all = add_simplified_metadata_columns(rscc_all)
    return (rscc_all,)


@app.cell
def _(bond_stats, clashes_1, make_delta_frame, rscc_all):
    _group_cols = ["protein", "model_name", "guidance_weight", "guidance_strength"]

    clash_stats = clashes_1.groupby(_group_cols, as_index=False).agg(
        clashscore=("clashscore", "mean")
    )
    bond_stats_agg = bond_stats.groupby(_group_cols, as_index=False).agg(
        outlier_fraction_bond_angles=("outlier_fraction_bond_angles", "mean"),
        outlier_fraction_bond_lengths=("outlier_fraction_bond_lengths", "mean"),
    )
    rscc_stats = rscc_all.groupby(_group_cols, as_index=False).agg(
        rscc=("rscc", "mean")
    )

    clash_delta = make_delta_frame(clash_stats, "clashscore")[
        _group_cols + ["delta_clashscore"]
    ]
    bond_angle_delta = make_delta_frame(bond_stats_agg, "outlier_fraction_bond_angles")[
        _group_cols + ["delta_outlier_fraction_bond_angles"]
    ]
    bond_length_delta = make_delta_frame(bond_stats_agg, "outlier_fraction_bond_lengths")[
        _group_cols + ["delta_outlier_fraction_bond_lengths"]
    ]
    rscc_delta_tradeoff = make_delta_frame(rscc_stats, "rscc")[
        _group_cols + ["delta_rscc"]
    ]
    return (
        bond_angle_delta,
        bond_length_delta,
        clash_delta,
        rscc_delta_tradeoff,
    )


@app.cell
def _(bond_angle_delta, bond_length_delta, clash_delta, rscc_delta_tradeoff):
    _merge_on = ["protein", "model_name", "guidance_weight", "guidance_strength"]

    tradeoff_df = (
        clash_delta
        .merge(bond_angle_delta, on=_merge_on)
        .merge(bond_length_delta, on=_merge_on)
        .merge(rscc_delta_tradeoff, on=_merge_on, how="left")
    )

    tradeoff_df["delta_outlier_fraction_bond_lengths_pct"] = (
        100.0 * tradeoff_df["delta_outlier_fraction_bond_lengths"]
    )
    tradeoff_df["delta_outlier_fraction_bond_angles_pct"] = (
        100.0 * tradeoff_df["delta_outlier_fraction_bond_angles"]
    )

    _geom_terms = [
        "delta_clashscore",
        "delta_outlier_fraction_bond_lengths_pct",
        "delta_outlier_fraction_bond_angles_pct",
    ]
    tradeoff_df["geometry_penalty"] = tradeoff_df[_geom_terms].mean(axis=1)
    return (tradeoff_df,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Figures 12–19 — Spiral Trajectory per Model

    Each figure traces every protein through (delta RSCC, geometry penalty) space
    as guidance weight increases. The top-N most extreme trajectories are labeled.

    Select a model from the dropdown above to view its spiral. Toggle
    **Save all models at once** to batch-save all eight spirals in one run.

    **Output:** `figures/12_spiral_boltz2_xrd.svg` through `figures/19_spiral_rf3_chiral_off.svg`
    """)
    return


@app.cell
def _(
    MODEL_DISPLAY,
    auto_save,
    make_spiral_figure,
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
        m for m in [
            "boltz2X-RAY", "boltz2MD", "protenix", "rf3",
            "protenix-mini", "protenix-tiny", "rf3-old", "rf3-chiral-off",
        ]
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
        # Batch save all models
        for _model in _all_models:
            _model_df = tradeoff_df[tradeoff_df["model_name"] == _model].copy()
            if _model_df.empty:
                continue
            _slug = MODEL_DISPLAY.get(_model, _model).lower().replace("-", "_")
            _fig = make_spiral_figure(_model_df, _model, xlim=_shared_xlim, ylim=_shared_ylim)
            save_fig(_fig, f"{_spiral_nums[_model]:02d}_spiral_{_slug}.png")
            plt.close(_fig)

    # Always show the currently selected model
    _cur_model = selected_model.value
    if _cur_model in tradeoff_df["model_name"].unique():
        _model_df = tradeoff_df[tradeoff_df["model_name"] == _cur_model].copy()
        _slug = MODEL_DISPLAY.get(_cur_model, _cur_model).lower().replace("-", "_")
        spiral_fig = make_spiral_figure(_model_df, _cur_model, xlim=_shared_xlim, ylim=_shared_ylim)
        if auto_save.value and not save_all.value:
            save_fig(spiral_fig, f"{_spiral_nums[_cur_model]:02d}_spiral_{_slug}.png")
        spiral_fig
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Figures 20–28 — Aggregated Spiral Trajectories

    These figures aggregate over proteins. Each guidance weight becomes one
    point at the **median** Δ RSCC and **median** geometry penalty across
    proteins; error bars show the **interquartile range (IQR)**. The Y axis
    uses a symmetric log scale so differences at small guidance weights
    remain visible.

    Median + IQR (rather than mean ± SD) keeps the central tendency robust
    to a few extreme proteins in the heavy tail of `geometry_penalty`.

    **Output:**
    - `figures/20_aggregated_spiral_all.svg` — all models overlaid
    - `figures/21_aggregated_spiral_boltz2_xrd.svg` … `28_aggregated_spiral_rf3_chiral_off.svg` — per model
    """)
    return


@app.cell
def _(auto_save, make_aggregated_spiral_figure, save_fig, tradeoff_df):
    agg_spiral_all_fig = make_aggregated_spiral_figure(tradeoff_df)
    if auto_save.value:
        save_fig(agg_spiral_all_fig, "20_aggregated_spiral_all.svg")
    agg_spiral_all_fig
    return


@app.cell
def _(
    MODEL_DISPLAY,
    auto_save,
    make_aggregated_spiral_figure,
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
        m for m in [
            "boltz2X-RAY", "boltz2MD", "protenix", "rf3",
            "protenix-mini", "protenix-tiny", "rf3-old", "rf3-chiral-off",
        ]
        if m in tradeoff_df["model_name"].unique()
    ]

    if save_all.value and auto_save.value:
        for _model in _all_models_agg:
            _slug = MODEL_DISPLAY.get(_model, _model).lower().replace("-", "_")
            _fig = make_aggregated_spiral_figure(tradeoff_df, model_key=_model)
            save_fig(_fig, f"{_agg_spiral_nums[_model]:02d}_aggregated_spiral_{_slug}.svg")
            plt.close(_fig)

    _cur_model = selected_model.value
    if _cur_model in tradeoff_df["model_name"].unique():
        _slug = MODEL_DISPLAY.get(_cur_model, _cur_model).lower().replace("-", "_")
        agg_spiral_fig = make_aggregated_spiral_figure(tradeoff_df, model_key=_cur_model)
        if auto_save.value and not save_all.value:
            save_fig(agg_spiral_fig, f"{_agg_spiral_nums[_cur_model]:02d}_aggregated_spiral_{_slug}.svg")
        agg_spiral_fig
    return


if __name__ == "__main__":
    app.run()
