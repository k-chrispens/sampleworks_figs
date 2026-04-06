import marimo

__generated_with = "0.20.4"
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
    import seaborn as sns
    from analysis_helpers import (
        add_simplified_metadata_columns,
        load_rscc,
        load_lddt,
        make_ridge_plot,
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
        MODELS_ALT,
        MODELS_PROTENIX_GROUP,
        MODELS_RF3_GROUP,
        add_simplified_metadata_columns,
        load_lddt,
        load_rscc,
        make_delta_frame,
        make_ridge_plot,
        pd,
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
def _(
    MODELS_ALT,
    MODELS_PROTENIX_GROUP,
    MODELS_RF3_GROUP,
    auto_save,
    make_ridge_plot,
    model_group,
    rscc_all,
    save_fig,
):
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
    return (fig_04,)


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
def _(
    MODELS_ALT,
    MODELS_PROTENIX_GROUP,
    MODELS_RF3_GROUP,
    auto_save,
    make_ridge_plot,
    model_group,
    rscc_all,
    save_fig,
):
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
    return (fig_05,)


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
    MODELS_ALT,
    MODELS_PROTENIX_GROUP,
    MODELS_RF3_GROUP,
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
    return (fig_06,)


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
    return (fig_07,)


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
def _(
    MODELS_PROTENIX_GROUP,
    auto_save,
    make_delta_frame,
    make_ridge_plot,
    rscc_all,
    save_fig,
):
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
    return (fig_08,)


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
    fig_09 = make_ridge_plot(
        lddt_results_1.query("guidance_strength == 'baseline' and selection_name != 'all altloc'"),
        plot_column_name="avg_silhouette_to_ref",
        custom_xlabel="Average Silhouette Score to Reference (LDDT Clustering)",
        smoothing=smoothing.value,
        combined=True,
    )
    if auto_save.value:
        save_fig(fig_09, "09_lddt_to_ref_baseline_distribution.png")
    fig_09
    return (fig_09,)


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
        lddt_results_1.query("selection_name != 'all altloc'"),
        "avg_silhouette_to_ref",
        smoothing=smoothing.value,
    )
    if auto_save.value:
        save_fig(fig_10, "10_lddt_to_ref_all_guidance_distribution.png")
    fig_10
    return (fig_10,)


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
    fig_11 = make_ridge_plot(lddt_merged, "delta_silhouette_to_ref", is_delta=True)
    if auto_save.value:
        save_fig(fig_11, "11_delta_lddt_to_ref_distribution.png")
    fig_11
    return (fig_11,)


if __name__ == "__main__":
    app.run()
