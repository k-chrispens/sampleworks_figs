import marimo

__generated_with = "0.23.1"
app = marimo.App()


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell
def _():
    import pandas as pd
    from analysis_helpers import (
        filter_by_model_group,
        get_rscc_nan_proteins,
        load_and_tag_path_csv,
        tag_path_dataframe,
        make_ridge_plot,
        make_delta_frame,
        save_fig,
    )

    return (
        filter_by_model_group,
        get_rscc_nan_proteins,
        load_and_tag_path_csv,
        make_delta_frame,
        make_ridge_plot,
        pd,
        save_fig,
        tag_path_dataframe,
    )


@app.cell
def _(mo):
    mo.md(r"""
    # Geometry Quality Analysis

    This notebook analyzes the geometric quality of protein structure predictions
    across different models and guidance strengths.

    ## Figures produced by this notebook

    | ID | Description | Output |
    |----|-------------|--------|
    | 01 | Clashscore Distribution | `01_clashscore_distribution[_variant].svg` |
    | 02 | Delta Clashscore Distribution | `02_delta_clashscore_distribution[_variant].svg` |
    | 03 | Delta Bond Angle Outlier Fraction | `03_delta_bond_angle_outlier_fraction_distribution[_variant].svg` |
    """)

    model_group = mo.ui.dropdown(
        options={"Main models": "main", "Protenix variants": "protenix", "RF3 variants": "rf3"},
        value="Main models",
        label="Model group",
    )
    smoothing = mo.ui.slider(
        start=0.05, stop=1.0, step=0.05,
        value=0.15,
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
def _(load_and_tag_path_csv, nan_rscc_proteins, pd, tag_path_dataframe):
    clashes = load_and_tag_path_csv(
        "initial40results/clashscore_metrics.csv",
        filename_col="model_name",
        extra_cols={"selection": "all residues"},
    )
    _bond_angles = pd.read_csv("initial40results/bond_angle_violation_fractions.csv")
    _bond_lengths = pd.read_csv("initial40results/bond_length_violation_fractions.csv")
    bond_stats = _bond_angles.merge(
        _bond_lengths, on=["model", "model_n"], suffixes=["_bond_angles", "_bond_lengths"]
    )
    bond_stats = tag_path_dataframe(bond_stats, filename_col="model")
    if nan_rscc_proteins:
        clashes = clashes[~clashes["protein"].isin(nan_rscc_proteins)]
        bond_stats = bond_stats[~bond_stats["protein"].isin(nan_rscc_proteins)]
    return bond_stats, clashes


@app.cell
def _(
    auto_save,
    clashes,
    filter_by_model_group,
    make_ridge_plot,
    mo,
    model_group,
    save_fig,
    smoothing,
):
    _filt = clashes.query("guidance_strength != 'strong' and selection_name != 'all altloc'")
    _data, _suffix = filter_by_model_group(_filt, model_group.value)
    fig_01 = make_ridge_plot(_data, "clashscore", smoothing=smoothing.value)
    if auto_save.value:
        save_fig(fig_01, f"01_clashscore_distribution{_suffix}.png")
    mo.vstack([
        mo.md("""## Figure 01 -- Clashscore Distribution
        Ridge plot of Phenix clashscores (strong guidance excluded).
        `figures/01_clashscore_distribution[_variant].svg`"""),
        fig_01,
    ])
    return


@app.cell
def _(
    auto_save,
    clashes,
    filter_by_model_group,
    make_delta_frame,
    make_ridge_plot,
    mo,
    model_group,
    save_fig,
    smoothing,
):
    _diff = make_delta_frame(clashes, "clashscore")
    _data, _suffix = filter_by_model_group(_diff, model_group.value, extra_query="delta_clashscore < 100")
    fig_02 = make_ridge_plot(_data, "delta_clashscore", is_delta=True, smoothing=smoothing.value)
    if auto_save.value:
        save_fig(fig_02, f"02_delta_clashscore_distribution{_suffix}.png")
    mo.vstack([
        mo.md("""## Figure 02 -- Delta Clashscore Distribution
        Change in clashscore from baseline (outliers > 100 excluded).
        `figures/02_delta_clashscore_distribution[_variant].svg`"""),
        fig_02,
    ])
    return


@app.cell
def _(
    auto_save,
    bond_stats,
    filter_by_model_group,
    make_delta_frame,
    make_ridge_plot,
    mo,
    model_group,
    save_fig,
    smoothing,
):
    _delta = make_delta_frame(bond_stats, "outlier_fraction_bond_angles")
    _data, _suffix = filter_by_model_group(
        _delta, model_group.value, extra_query="delta_outlier_fraction_bond_angles < 0.05",
    )
    fig_03 = make_ridge_plot(
        _data, "delta_outlier_fraction_bond_angles",
        is_delta=True, smoothing=smoothing.value,
    )
    if auto_save.value:
        save_fig(fig_03, f"03_delta_bond_angle_outlier_fraction_distribution{_suffix}.png")
    mo.vstack([
        mo.md("""## Figure 03 -- Delta Bond Angle Outlier Fraction
        Change in bond-angle outlier fraction from baseline (outliers >= 0.05 excluded).
        `figures/03_delta_bond_angle_outlier_fraction_distribution[_variant].svg`"""),
        fig_03,
    ])
    return


if __name__ == "__main__":
    app.run()
