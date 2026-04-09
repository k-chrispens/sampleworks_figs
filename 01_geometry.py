import marimo

__generated_with = "0.23.0"
app = marimo.App()


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell
def _():
    import pandas as pd
    from analysis_helpers import (
        load_and_tag_path_csv,
        tag_path_dataframe,
        make_ridge_plot,
        make_delta_frame,
        save_fig,
        MODELS_ALT,
        MODELS_PROTENIX_GROUP,
        MODELS_RF3_GROUP,
    )

    return (
        load_and_tag_path_csv,
        make_delta_frame,
        make_ridge_plot,
        pd,
        save_fig,
        tag_path_dataframe,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Geometry Quality Analysis

    This notebook analyzes the geometric quality of protein structure predictions
    across different models and guidance strengths. It examines two key aspects:

    1. **Clash scores** — steric clashes detected by Phenix, measuring how often
       atoms overlap in sterically forbidden ways.
    2. **Bond angle outlier fractions** — the fraction of bond angles that fall
       outside expected geometric ranges, indicating distortions in molecular geometry.

    Each metric is shown as a distribution across proteins, and as a delta
    (change from the no-guidance baseline) to reveal the effect of guidance.
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
        value=0.15,
        label="KDE smoothing",
    )
    auto_save = mo.ui.switch(label="Auto-save figures", value=False)
    mo.hstack([model_group, smoothing, auto_save], justify="start", gap=1)
    return auto_save, model_group, smoothing


@app.cell
def _(load_and_tag_path_csv):
    clashes = load_and_tag_path_csv(
        "initial40results/clashscore_metrics.csv",
        filename_col="model_name",
        extra_cols={"selection": "all residues"},
    )
    return (clashes,)


@app.cell
def _(pd, tag_path_dataframe):
    _bond_angles = pd.read_csv("initial40results/bond_angle_violation_fractions.csv")
    _bond_lengths = pd.read_csv("initial40results/bond_length_violation_fractions.csv")
    bond_stats = _bond_angles.merge(
        _bond_lengths, on=["model", "model_n"], suffixes=["_bond_angles", "_bond_lengths"]
    )
    bond_stats = tag_path_dataframe(bond_stats, filename_col="model")
    return (bond_stats,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Figure 01 — Clashscore Distribution

    Ridge plot showing the distribution of Phenix clashscores across all proteins
    for each model at baseline and weak guidance strengths (strong guidance excluded).
    Lower clashscores indicate better steric geometry.

    **Output:** `figures/01_clashscore_distribution[_variant].svg`
    """)
    return


@app.cell
def _(auto_save, clashes, make_ridge_plot, model_group, save_fig, smoothing):
    _filt = clashes.query("guidance_strength != 'strong' and selection_name != 'all altloc'")
    _group_val = model_group.value
    if _group_val == "main":
        _data = _filt.query("model_name not in @MODELS_ALT")
        _suffix = ""
    elif _group_val == "protenix":
        _data = _filt.query("model_name in @MODELS_PROTENIX_GROUP")
        _suffix = "_protenix_alts"
    else:
        _data = _filt.query("model_name in @MODELS_RF3_GROUP")
        _suffix = "_rf3_alts"
    fig_01 = make_ridge_plot(_data, "clashscore", smoothing=smoothing.value)
    if auto_save.value:
        save_fig(fig_01, f"01_clashscore_distribution{_suffix}.png")
    fig_01
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Figure 02 — Delta Clashscore Distribution

    Ridge plot of the *change* in clashscore relative to the no-guidance baseline
    for each model. Positive values mean guidance increased clashes; negative means
    guidance reduced them. Outliers with delta > 100 are excluded for readability.

    **Output:** `figures/02_delta_clashscore_distribution[_variant].svg`
    """)
    return


@app.cell
def _(
    auto_save,
    clashes,
    make_delta_frame,
    make_ridge_plot,
    model_group,
    save_fig,
    smoothing,
):
    _diff = make_delta_frame(clashes, "clashscore")
    _group_val = model_group.value
    if _group_val == "main":
        _data = _diff.query("delta_clashscore < 100 and model_name not in @MODELS_ALT")
        _suffix = ""
    elif _group_val == "protenix":
        _data = _diff.query("delta_clashscore < 100 and model_name in @MODELS_PROTENIX_GROUP")
        _suffix = "_protenix_alts"
    else:
        _data = _diff.query("delta_clashscore < 100 and model_name in @MODELS_RF3_GROUP")
        _suffix = "_rf3_alts"
    fig_02 = make_ridge_plot(_data, "delta_clashscore", is_delta=True, smoothing=smoothing.value)
    if auto_save.value:
        save_fig(fig_02, f"02_delta_clashscore_distribution{_suffix}.png")
    fig_02
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Figure 03 — Delta Bond Angle Outlier Fraction

    Ridge plot of the change in bond-angle outlier fraction relative to the
    no-guidance baseline. Positive values indicate guidance caused more bond-angle
    distortions. Extreme outliers (delta >= 0.05) are excluded for clarity.

    **Output:** `figures/03_delta_bond_angle_outlier_fraction_distribution[_variant].svg`
    """)
    return


@app.cell
def _(
    auto_save,
    bond_stats,
    make_delta_frame,
    make_ridge_plot,
    model_group,
    save_fig,
    smoothing,
):
    _delta = make_delta_frame(bond_stats, "outlier_fraction_bond_angles")
    _group_val = model_group.value
    if _group_val == "main":
        _data = _delta.query("delta_outlier_fraction_bond_angles < 0.05 and model_name not in @MODELS_ALT")
        _suffix = ""
    elif _group_val == "protenix":
        _data = _delta.query("delta_outlier_fraction_bond_angles < 0.05 and model_name in @MODELS_PROTENIX_GROUP")
        _suffix = "_protenix_alts"
    else:
        _data = _delta.query("delta_outlier_fraction_bond_angles < 0.05 and model_name in @MODELS_RF3_GROUP")
        _suffix = "_rf3_alts"
    fig_03 = make_ridge_plot(
        _data, "delta_outlier_fraction_bond_angles",
        is_delta=True, smoothing=smoothing.value,
    )
    if auto_save.value:
        save_fig(fig_03, f"03_delta_bond_angle_outlier_fraction_distribution{_suffix}.png")
    fig_03
    return


if __name__ == "__main__":
    app.run()
