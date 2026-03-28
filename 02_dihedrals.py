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
        load_and_tag_path_csv_fast,
        save_fig,
    )
    return (
        load_and_tag_path_csv,
        load_and_tag_path_csv_fast,
        np,
        pd,
        plt,
        save_fig,
        sns,
    )


@app.cell
def _(mo):
    load_residues = mo.ui.switch(label="Load per-residue data (1.8 GB)", value=False)
    mo.vstack([load_residues])
    return (load_residues,)


# ── Backbone and Sidechain Dihedrals ─────────────────────────────────────────


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Backbone and side chain dihedrals
    """)
    return


@app.cell
def _(load_and_tag_path_csv):
    dihedrals = load_and_tag_path_csv(
        "initial40results/initial_40_tortoize_protein_stats.csv",
        filename_col="path",
        extra_renames={"model": "model_n"},
    )
    dihedrals
    return (dihedrals,)


@app.cell
def _(dihedrals):
    dihedrals_summary = dihedrals.groupby(
        ["model_name", "pdb_id", "guidance_weight"]
    ).agg({"ramachandran_z_score": "mean", "torsion_z_score": "mean"})
    dihedrals_summary.groupby(["model_name", "guidance_weight"]).describe(
        percentiles=[0.1, 0.5, 0.9]
    )
    return (dihedrals_summary,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Wow that looks... awful. By eyeball it looks like the Boltz models again do better, but there are still quite a few bad proteins in there. As before it seems things get worse with stronger guidance. Curiously stronger guidance seems to shift the z-scores to be more negative. It isn't clear to me why that should be, we should discuss. Let's see if the distributions shift a lot. One thing we notice is that there are NaN rows--something went seriously wrong but fortunately in just a few cases.

    The torsion scores are bonkers--I need to dig in to this more, but it seems that in aggregate across each protein, the distribution must be biased enough to cause the individual values not to seem too nutty, but taken together they are extremely unlikely.
    """)
    return


@app.cell
def _(dihedrals_summary):
    (
        dihedrals_summary.reset_index()
        .query(
            "(model_name != 'rf3' and guidance_weight == 0.1) or "
            "(model_name == 'rf3' and guidance_weight == 0.01)"
        )
        .pivot(
            index="pdb_id",
            columns="model_name",
            values=["ramachandran_z_score", "torsion_z_score"],
        )
    )
    return


@app.cell
def _(dihedrals_summary, pd, plt, save_fig, sns):
    _dsw = (
        dihedrals_summary.reset_index()
        .query(
            "(model_name != 'rf3' and guidance_weight == 0.1) or "
            "(model_name == 'rf3' and guidance_weight == 0.01)"
        )
        .dropna()
        .assign(guidance="weak")
    )
    _dss = (
        dihedrals_summary.reset_index()
        .query(
            "(model_name != 'rf3' and guidance_weight == 0.5) or "
            "(model_name == 'rf3' and guidance_weight == 0.05)"
        )
        .dropna()
        .assign(guidance="strong")
    )
    _dsb = (
        dihedrals_summary.reset_index()
        .query("guidance_weight == 0.0")
        .dropna()
        .assign(guidance="baseline")
    )
    ds = pd.concat([_dsw, _dss, _dsb])

    _hue_order = ["weak", "strong", "baseline"]
    _palette = {"weak": "#377EB8", "strong": "#E41A1C", "baseline": "#888888"}
    _g = sns.FacetGrid(
        ds, row="model_name", hue="guidance", hue_order=_hue_order,
        aspect=5, height=1.5, palette=_palette,
    )
    _g.map(sns.kdeplot, "ramachandran_z_score", bw_adjust=0.5, clip_on=False, fill=True, alpha=0.75, linewidth=1.5)
    _g.map(sns.kdeplot, "ramachandran_z_score", clip_on=False, color="w", lw=2, bw_adjust=0.5)
    _g.refline(y=0, linewidth=2, linestyle="-", color=None, clip_on=False)

    def _label(*x, color=None, label=None):
        text = x[0].iloc[0]
        ax = plt.gca()
        ax.text(0, 0.2, text, fontweight="bold", color="black", ha="left", va="center", transform=ax.transAxes)
        ax.text(0, 0.4 if label == "strong" else 0.6, label, fontweight="bold", color=color, ha="left", va="center", transform=ax.transAxes)

    _g.map(_label, "model_name")
    _g.set(xlabel="Ramachandran Z-Score")
    _g.set_titles("")
    _g.set(yticks=[0.0, 0.2, 0.4], ylabel="")
    _g.set(xticks=[-9, -6, -3, 0, 3, 6])
    _g.despine(bottom=True, left=True)
    save_fig(_g, "04_ramachandran_zscore_kde_by_model.png")
    return (ds,)


@app.cell
def _(ds, plt, save_fig, sns):
    _hue_order = ["weak", "strong"]
    _palette = {"weak": "#377EB8", "strong": "#E41A1C"}
    _g = sns.FacetGrid(
        ds, row="model_name", hue="guidance", hue_order=_hue_order,
        aspect=5, height=1.5, palette=_palette,
    )
    _g.map(sns.kdeplot, "torsion_z_score", bw_adjust=0.5, clip_on=False, fill=True, alpha=0.75, linewidth=1.5)
    _g.map(sns.kdeplot, "torsion_z_score", clip_on=False, color="w", lw=2, bw_adjust=0.5)
    _g.refline(y=0, linewidth=2, linestyle="-", color=None, clip_on=False)

    def _label(x, color, label):
        text = x.iloc[0]
        ax = plt.gca()
        ax.text(0, 0.2, text, fontweight="bold", color="black", ha="left", va="center", transform=ax.transAxes)
        ax.text(0, 0.4 if label == "strong" else 0.6, label, fontweight="bold", color=color, ha="left", va="center", transform=ax.transAxes)

    _g.map(_label, "model_name")
    _g.set(xlabel="Torsion Z-Score")
    _g.set_titles("")
    _g.set(yticks=[0.0, 0.2, 0.4], ylabel="")
    _g.set(xticks=[-9, -6, -3, 0, 3, 6])
    _g.despine(bottom=True, left=True)
    save_fig(_g, "05_torsion_zscore_kde_by_model.png")
    return


# ── Per-Residue Analysis (optional heavy data) ───────────────────────────────


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Torsion outlier counts

    Super weird. Let's take a look at the actual number of outliers. Turns out they are quite small.

    **Toggle the switch above to load per-residue data (1.8 GB).**
    """)
    return


@app.cell
def _(load_and_tag_path_csv_fast, load_residues, mo):
    mo.stop(not load_residues.value)
    torsion_res = load_and_tag_path_csv_fast(
        "initial40results/initial_40_tortoize_residues.csv",
    )
    torsion_res
    return (torsion_res,)


@app.cell
def _(torsion_res):
    torsion_res.query("pdb_id == '1VME'").groupby(["model_name", "guidance_weight"])[
        ["ramachandran_z_score", "torsion_z_score"]
    ].describe()
    return


@app.cell
def _(torsion_res):
    torsion_res.query("pdb_id == '1VME'").groupby(["model_name", "guidance_weight"]).agg(
        rama_outlier_frac=("ramachandran_z_score", lambda x: (abs(x) > 3).sum() / len(x)),
        torsion_outlier_frac=("torsion_z_score", lambda x: (abs(x) > 3).sum() / len(x)),
    )
    return


@app.cell
def _(torsion_res):
    torsion_outliers = torsion_res.groupby(
        ["model_n", "model_name", "guidance_weight", "pdb_id"]
    ).agg(
        rama_outlier_frac=("ramachandran_z_score", lambda x: (abs(x) > 3).sum() / len(x)),
        torsion_outlier_frac=("torsion_z_score", lambda x: (abs(x) > 3).sum() / len(x)),
    )
    torsion_outliers
    return (torsion_outliers,)


@app.cell
def _(torsion_outliers):
    torsion_summary = torsion_outliers.groupby(
        ["model_name", "pdb_id", "guidance_weight"]
    ).agg({"rama_outlier_frac": "mean", "torsion_outlier_frac": "mean"})
    torsion_summary.groupby(["model_name", "guidance_weight"]).describe(
        percentiles=[0.1, 0.5, 0.9]
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Well, I'm not sure what to make of this. The Ramachandran outliers seem to be fewer (not necessarily better) with more guidance. This appears to be mostly true for torsion (side chain X_1, X_2) angles. But it is clear that the whole-protein scores get _much_ worse, if I'm correctly interpreting the whole protein scores. This implies that while there are fewer outliers with higher guidance, the outliers there are must get much worse. Let's see what the distributions look like:
    """)
    return


@app.cell
def _(torsion_res):
    dihedral_means = torsion_res.groupby(
        ["model_name", "guidance_weight", "pdb_id"]
    ).agg({"ramachandran_z_score": "mean", "torsion_z_score": "mean"})
    return (dihedral_means,)


@app.cell
def _(dihedral_means, pd, plt, save_fig, sns):
    _dsw = (
        dihedral_means.reset_index()
        .query(
            "(model_name != 'rf3' and guidance_weight == 0.1) or "
            "(model_name == 'rf3' and guidance_weight == 0.01)"
        )
        .dropna()
        .assign(guidance="weak")
    )
    _dss = (
        dihedral_means.reset_index()
        .query(
            "(model_name != 'rf3' and guidance_weight == 0.5) or "
            "(model_name == 'rf3' and guidance_weight == 0.05)"
        )
        .dropna()
        .assign(guidance="strong")
    )
    ds_1 = pd.concat([_dsw, _dss])

    _hue_order = ["weak", "strong"]
    _palette = {"weak": "#377EB8", "strong": "#E41A1C"}
    _g = sns.FacetGrid(
        ds_1, row="model_name", hue="guidance", hue_order=_hue_order,
        aspect=5, height=1.5, palette=_palette,
    )
    _g.map(sns.kdeplot, "ramachandran_z_score", bw_adjust=0.5, clip_on=False, fill=True, alpha=0.75, linewidth=1.5)
    _g.map(sns.kdeplot, "ramachandran_z_score", clip_on=False, color="w", lw=2, bw_adjust=0.5)
    _g.refline(y=0, linewidth=2, linestyle="-", color=None, clip_on=False)

    def _label(x, color, label):
        text = x.iloc[0]
        ax = plt.gca()
        ax.text(0, 0.2, text, fontweight="bold", color="black", ha="left", va="center", transform=ax.transAxes)
        ax.text(0, 0.4 if label == "strong" else 0.6, label, fontweight="bold", color=color, ha="left", va="center", transform=ax.transAxes)

    _g.map(_label, "model_name")
    _g.set(xlabel="Ramachandran Z-Score")
    _g.set_titles("")
    _g.set(yticks=[0, 2, 4], ylabel="")
    _g.despine(bottom=True, left=True)
    save_fig(_g, "06_ramachandran_zscore_residue_mean_kde.png")
    return (ds_1,)


@app.cell
def _(ds_1, plt, save_fig, sns):
    _hue_order = ["weak", "strong"]
    _palette = {"weak": "#377EB8", "strong": "#E41A1C"}
    _g = sns.FacetGrid(
        ds_1, row="model_name", hue="guidance", hue_order=_hue_order,
        aspect=5, height=1.5, palette=_palette,
    )
    _g.map(sns.kdeplot, "torsion_z_score", bw_adjust=0.5, clip_on=False, fill=True, alpha=0.75, linewidth=1.5)
    _g.map(sns.kdeplot, "torsion_z_score", clip_on=False, color="w", lw=2, bw_adjust=0.5)
    _g.refline(y=0, linewidth=2, linestyle="-", color=None, clip_on=False)

    def _label(x, color, label):
        text = x.iloc[0]
        ax = plt.gca()
        ax.text(0, 0.2, text, fontweight="bold", color="black", ha="left", va="center", transform=ax.transAxes)
        ax.text(0, 0.4 if label == "strong" else 0.6, label, fontweight="bold", color=color, ha="left", va="center", transform=ax.transAxes)

    _g.map(_label, "model_name")
    _g.set(xlabel="Torsion Z-Score")
    _g.set_titles("")
    _g.set(yticks=[0, 2, 4], ylabel="")
    _g.despine(bottom=True, left=True)
    save_fig(_g, "07_torsion_zscore_residue_mean_kde.png")
    return


if __name__ == "__main__":
    app.run()
