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


@app.cell
def _(mo):
    smoothing = mo.ui.slider(
        start=0.05, stop=1.0, step=0.05,
        value=0.15,
        label="KDE smoothing",
    )
    mo.vstack([smoothing])
    return (smoothing,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Clashes

    Let's take a look first at the clashes detected using phenix. We'll look at clashscores first.
    """)
    return


@app.cell
def _(load_and_tag_path_csv):
    clashes = load_and_tag_path_csv(
        "initial40results/clashscore_metrics.csv",
        filename_col="model_name",
        extra_cols={"selection": "all residues"},
    )
    return (clashes,)


@app.cell
def _(clashes):
    clash_summary = clashes.groupby(
        ["model_name", "pdb_id", "guidance_weight"]
    ).agg({"num_clashes": "mean", "clashscore": "mean"})
    return (clash_summary,)


@app.cell
def _(clash_summary):
    clash_summary
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Make sure all proteins are represented well across the models

    We want to be careful that our results aren't biased by any missing results.
    """)
    return


@app.cell
def _(clashes):
    clashes.groupby("pdb_id").agg({"model_name": "count"})
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    First let's see whether all the models are available (all should have 8, they do):
    """)
    return


@app.cell
def _(clashes):
    (clashes.groupby(["model_name", "pdb_id", "guidance_weight"]).model_id.count() == 8).all()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Pick the first one, check whether all models are represented for each protein:

    On my first pass through, the answer was an emphatic "no". There were missing results from Boltz for a large number of proteins, most of them, in fact. Some tinkering indicated that phenix.clashscore (and presumably all of phenix) wouldn't recognize a CIF file that lacked "B_iso_or_equiv", and possibly also "occupancy". I will modify my CIF patching script to set B factors to 20.0 if missing, and set occupancies to 1.0 for all atoms.

    I'm not sure if these matter for clashes or not. I suppose they could go into how "bad" a clash is.
    """)
    return


@app.cell
def _(clashes):
    clashes.query("model_id == 1").pivot_table(
        index="pdb_id",
        columns=["model_name", "guidance_weight"],
        values="clashscore",
        aggfunc="count",
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Clash scores and number of clashes by model and guidance weight

    The first question I have is whether higher guidance weight leads to more clashes. The highest guidance we tried for each model apparently leads to 3-10x more clashes and similarly worse Phenix clashscores. The results validate our choice to use 10x lower guidance weight for RF3, although even lower guidance may lead to better final results.

    Based on these results, it appears protenix results in the smallest clashscores and number of clashes, by a small margin (probably not statistically significant) over Boltz2-XRD although the worst-case models for Boltz2-XRD are vastly better than the worst-case for Protenix. Boltz2-MD does not perform as well as the X-ray version. We need to understand the clashscore a little better too, but it appears that while Boltz2-XRD has better median clashscores while Protenix has lower median number of clashes.
    """)
    return


@app.cell
def _(clash_summary):
    clash_summary.groupby(["model_name", "guidance_weight"]).describe(percentiles=[0.5])
    return


@app.cell
def _(clashes, make_ridge_plot, save_fig, smoothing):
    _filt = clashes.query("guidance_strength != 'strong' and selection_name != 'all altloc'")
    _g = make_ridge_plot(_filt.query("model_name not in @MODELS_ALT"), "clashscore", smoothing=smoothing.value)
    save_fig(_g, "01_clashscore_distribution.png")
    return


@app.cell
def _(clashes, make_ridge_plot, save_fig, smoothing):
    _g = make_ridge_plot(
        clashes.query("guidance_strength != 'strong' and selection_name != 'all altloc' and model_name in @MODELS_PROTENIX_GROUP"),
        "clashscore", smoothing=smoothing.value,
    )
    save_fig(_g, "01_clashscore_distribution_protenix_alts.png")
    return


@app.cell
def _(clashes, make_ridge_plot, save_fig, smoothing):
    _g = make_ridge_plot(
        clashes.query("guidance_strength != 'strong' and selection_name != 'all altloc' and model_name in @MODELS_RF3_GROUP"),
        "clashscore", smoothing=smoothing.value,
    )
    save_fig(_g, "01_clashscore_distribution_rf3_alts.png")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Generally it looks like RF3 is having a harder time that the others with respect to clashes.

    Kind of hard to see trends with guidance, so let's look at the difference between each guidance level and the baseline.
    """)
    return


@app.cell
def _(clashes, make_delta_frame, make_ridge_plot, save_fig, smoothing):
    clashes_diff = make_delta_frame(clashes, "clashscore")
    _g = make_ridge_plot(
        clashes_diff.query("delta_clashscore < 100 and model_name not in @MODELS_ALT"),
        "delta_clashscore", is_delta=True, smoothing=smoothing.value,
    )
    save_fig(_g, "02_delta_clashscore_distribution.png")
    return


@app.cell
def _(clashes, make_delta_frame, make_ridge_plot, save_fig, smoothing):
    _diff = make_delta_frame(clashes, "clashscore")
    _g = make_ridge_plot(
        _diff.query("delta_clashscore < 100 and model_name in @MODELS_PROTENIX_GROUP"),
        "delta_clashscore", is_delta=True, smoothing=smoothing.value,
    )
    save_fig(_g, "02_delta_clashscore_distribution_protenix_alts.png")
    return


@app.cell
def _(clashes, make_delta_frame, make_ridge_plot, save_fig, smoothing):
    _diff = make_delta_frame(clashes, "clashscore")
    _g = make_ridge_plot(
        _diff.query("delta_clashscore < 100 and model_name in @MODELS_RF3_GROUP"),
        "delta_clashscore", is_delta=True, smoothing=smoothing.value,
    )
    save_fig(_g, "02_delta_clashscore_distribution_rf3_alts.png")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Let's look at whether there are particular proteins that are easier/harder.

    There is a lot of variability in clash scores across proteins, and it seems to be a function of protein, not model, in that it is rare for some models but not others to do poorly (or vice versa). However, the spread can be very wide, e.g. 4GMU the lowest clashscore is about 6 and the highest almost 30, or the even more extreme 5I09, where it ranges from about 4 to 117 (in both cases Boltz2-XRD is best, RF3 the worst).
    """)
    return


@app.cell
def _(clash_summary):
    by_pdb_id = (
        clash_summary.reset_index()
        .query(
            "(model_name not in ['rf3','rf3-old','rf3-chiral-off'] and guidance_weight == 0.1) or "
            "(model_name in ['rf3','rf3-old','rf3-chiral-off'] and guidance_weight == 0.01)"
        )
        .pivot(index="pdb_id", columns="model_name", values=["num_clashes", "clashscore"])
    )
    by_pdb_id
    return (by_pdb_id,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    In more than half the cases, Boltz2-XRD produces the best clash score for a given PDB, while Protenix wins about a third of the time. We hypothesize that this is due to Boltz2's physical geometry restraints. RF3 is not always the worst, but it is never the best.
    """)
    return


@app.cell
def _(by_pdb_id):
    by_pdb_id.num_clashes.idxmin(axis=1).reset_index(name="model_name").groupby("model_name").size()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Bond Length and Angle Violations

    Next we look at bond lengths and angles to understand whether the models have trouble with this.
    """)
    return


@app.cell
def _(pd, tag_path_dataframe):
    _bond_angles = pd.read_csv("initial40results/bond_angle_violation_fractions.csv")
    _bond_lengths = pd.read_csv("initial40results/bond_length_violation_fractions.csv")
    bond_stats = _bond_angles.merge(
        _bond_lengths, on=["model", "model_n"], suffixes=["_bond_angles", "_bond_lengths"]
    )
    bond_stats = tag_path_dataframe(bond_stats, filename_col="model")
    return (bond_stats,)


@app.cell
def _(bond_stats):
    bond_stats.groupby("pdb_id").agg({"model_name": "count"})
    return


@app.cell
def _(bond_stats):
    bond_stats_summary = bond_stats.groupby(
        ["model_name", "pdb_id", "guidance_strength"]
    ).agg({
        "outlier_fraction_bond_angles": "mean",
        "outlier_fraction_bond_lengths": "mean",
    })
    return (bond_stats_summary,)


@app.cell
def _(bond_stats_summary):
    bond_stats_summary.groupby(["model_name", "guidance_strength"]).describe(percentiles=[0.5])
    return


@app.cell
def _(bond_stats, make_delta_frame, make_ridge_plot, save_fig, smoothing):
    _delta = make_delta_frame(bond_stats, "outlier_fraction_bond_angles")
    _g = make_ridge_plot(
        _delta.query("delta_outlier_fraction_bond_angles < 0.05 and model_name not in @MODELS_ALT"),
        "delta_outlier_fraction_bond_angles",
        is_delta=True,
        smoothing=smoothing.value,
    )
    save_fig(_g, "03_delta_bond_angle_outlier_fraction_distribution.png")
    print(f"total bond angles in dataframe: {len(_delta)}")
    _delta.query("delta_outlier_fraction_bond_angles >= 0.05").groupby(
        ["protein", "model_name", "guidance_strength"]
    ).agg({
        "outlier_fraction_bond_angles": ["mean", "size"],
        "delta_outlier_fraction_bond_angles": ["mean", "max"],
    })
    return


@app.cell
def _(bond_stats, make_delta_frame, make_ridge_plot, save_fig, smoothing):
    _delta = make_delta_frame(bond_stats, "outlier_fraction_bond_angles")
    _g = make_ridge_plot(
        _delta.query("delta_outlier_fraction_bond_angles < 0.05 and model_name in @MODELS_PROTENIX_GROUP"),
        "delta_outlier_fraction_bond_angles", is_delta=True, smoothing=smoothing.value,
    )
    save_fig(_g, "03_delta_bond_angle_outlier_fraction_distribution_protenix_alts.png")
    return


@app.cell
def _(bond_stats, make_delta_frame, make_ridge_plot, save_fig, smoothing):
    _delta = make_delta_frame(bond_stats, "outlier_fraction_bond_angles")
    _g = make_ridge_plot(
        _delta.query("delta_outlier_fraction_bond_angles < 0.05 and model_name in @MODELS_RF3_GROUP"),
        "delta_outlier_fraction_bond_angles", is_delta=True, smoothing=smoothing.value,
    )
    save_fig(_g, "03_delta_bond_angle_outlier_fraction_distribution_rf3_alts.png")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Again the story seems to be that high guidance leads to more bond angle and length outliers, although the effect is quite small I would say. Protenix and RF3 are unsurprisingly the most susceptible to these effects. Protenix is especially susceptible relative to the no-guidance baseline. In this case I would venture (but haven't yet calculated) that the effect is statistically significant.

    We should discuss whether this is a sign that the models aren't able to generate _distributions_ of idealized conformers, but are "trying" to fit the guidance density by distorting individual conformers.

    The worst cases here are pretty spectacular: for RF3 the worst case scenario has nearly 70% bond length and bond angle outliers, and note this is even with minimal guidance. That worst case scenario is limited to just one of 40 proteins (5I09, which also had clash problems).
    """)
    return


@app.cell
def _(bond_stats_summary):
    (
        bond_stats_summary.reset_index()
        .query("guidance_strength == 'weak'")
        .pivot(
            index="pdb_id",
            columns="model_name",
            values=["outlier_fraction_bond_angles", "outlier_fraction_bond_lengths"],
        )
    )
    return


if __name__ == "__main__":
    app.run()
