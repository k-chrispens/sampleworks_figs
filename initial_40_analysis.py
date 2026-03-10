import marimo

__generated_with = "0.20.4"
app = marimo.App()


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    TODO:
    check if geometry stats are already bad at 0.5 guidance for boltz/protenix see if we can go higher weight.
    fix metadata storage
    """)
    return


@app.cell
def _():
    import os
    os.listdir()
    return


@app.cell
def _():
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import seaborn as sns
    import re
    from pathlib import Path

    path_pattern = r"([0-9A-Za-z]{4})_native_occ/(.*?)/(.*?)/ens8_gw([0-9\.]+)/"
    def extract_metadata_from_path(p: str) -> pd.Series:
        m = re.search(path_pattern, p)
        return pd.Series(
            {"pdb_id": m.group(1), "model_name": m.group(2), "steering": m.group(3), "guidance_weight": m.group(4)}
        )


    # ── Consistent visual identity ────────────────────────────────────────────────
    # ColorBrewer Set1 — designed for maximum discriminability
    _GUIDANCE_PALETTE = {
        "baseline": "#888888",  # neutral grey
        "weak":     "#377EB8",  # blue
        "medium":   "#FF7F00",  # orange
        "strong":   "#E41A1C",  # red
    }

    _MODEL_DISPLAY = {
        # path-extracted names
        "boltz2_X-RAY_DIFFRACTION": "Boltz2-XRD",
        "boltz2_MD":                "Boltz2-MD",
        "rf3":                      "RF3",
        "protenix":                 "Protenix",
        # combined names (used in rscc/lddt cells)
        "boltz2X-RAY_DIFFRACTION":  "Boltz2-XRD",
        "boltz2X-RAY":              "Boltz2-XRD",
        "boltz2xray":               "Boltz2-XRD",
        "boltz2MD":                 "Boltz2-MD",
        "boltz2md":                 "Boltz2-MD",
    }

    _NICE_XLABELS = {
        "delta_clashscore":                    "Δ Clash Score",
        "delta_outlier_fraction_bond_angles":  "Δ Bond Angle Outlier Fraction",
        "delta_outlier_fraction_bond_lengths": "Δ Bond Length Outlier Fraction",
        "delta_rscc":                          "Δ RSCC",
        "delta_silhouette_to_ref":             "Δ Avg. Silhouette to Ref",
        "clashscore":                          "Clash Score",
        "rscc":                                "RSCC",
        "avg_silhouette_to_ref":               "Avg. Silhouette to Ref",
        "ramachandran_z_score":                "Ramachandran Z-Score",
        "torsion_z_score":                     "Torsion Z-Score",
        "outlier_fraction_bond_angles":        "Bond Angle Outlier Fraction",
        "outlier_fraction_bond_lengths":       "Bond Length Outlier Fraction",
    }
    # ─────────────────────────────────────────────────────────────────────────────


    def make_model_comparison_plot(df, plot_column_name, smoothing=0.15, is_delta=False, xlim=None):
        """Ridge/joy-plot comparing distributions across model types and guidance strengths."""
        from matplotlib.patches import Patch

        hue_order = [k for k in ["baseline", "weak", "medium", "strong"]
                     if k in df["guidance_strength"].unique()]
        palette   = [_GUIDANCE_PALETTE[k] for k in hue_order]
        model_order = sorted(df["model_name"].unique())

        g = sns.FacetGrid(
            df,
            row="model_name", row_order=model_order,
            hue="guidance_strength", hue_order=hue_order,
            aspect=5, height=1.8,
            palette=palette, sharex=True, sharey=False,
        )

        g.map(sns.kdeplot, plot_column_name,
              bw_adjust=smoothing, clip_on=False,
              fill=True, alpha=0.75, linewidth=1.5)
        g.map(sns.kdeplot, plot_column_name,
              clip_on=False, color="w", lw=2, bw_adjust=smoothing)

        g.refline(y=0, linewidth=1.5, linestyle="-", color=None, clip_on=False)

        _xlim = xlim
        if is_delta:
            for ax in g.axes.flat:
                ax.axvline(x=0, color="#333333", lw=1.5, linestyle="--", alpha=0.6, zorder=3)

            if _xlim is None:
                _values = df[plot_column_name].replace([np.inf, -np.inf], np.nan).dropna()
                if not _values.empty:
                    _max_abs = float(np.nanmax(np.abs(_values.to_numpy())))
                    if _max_abs == 0:
                        _max_abs = 0.01
                    _pad = max(0.06 * _max_abs, 0.005)
                    _xlim = (-_max_abs - _pad, _max_abs + _pad)

        if _xlim is not None:
            g.set(xlim=_xlim)

        for ax in g.axes.flat:
            ax.spines["bottom"].set_visible(True)
            ax.tick_params(axis="x", which="both", bottom=True, labelbottom=True)

        g.despine(left=True)

        # Model name labels inside each row at the upper-left, with a white backing
        for ax, model in zip(g.axes.flat, model_order):
            display = _MODEL_DISPLAY.get(model, model)
            ax.text(0.02, 0.90, display, transform=ax.transAxes,
                    fontweight="bold", fontsize=10, va="top", ha="left", color="#222222",
                    bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.8, ec="none"))

        legend_handles = [
            Patch(facecolor=_GUIDANCE_PALETTE[k], alpha=0.75, label=k.capitalize())
            for k in hue_order
        ]
        g.figure.legend(
            handles=legend_handles,
            loc="upper right", frameon=True, framealpha=0.9, edgecolor="none",
            title="Guidance", fontsize=9, title_fontsize=9,
            bbox_to_anchor=(0.99, 0.99),
        )

        xlabel = _NICE_XLABELS.get(plot_column_name,
                                   plot_column_name.replace("_", " ").title())
        g.set(xlabel=xlabel)

        g.figure.subplots_adjust(hspace=0.6, right=0.95)
        g.set_titles("")
        g.set(yticks=[], yticklabels=[], ylabel="")

        return g


    def make_model_comparison_delta_frame(df, plot_column_name, baseline_name="baseline"):
        merge_cols = ["protein", "model_name"]
        if "selection_name" in df.columns:
            merge_cols.append("selection_name")
        baseline = df.query('guidance_strength == @baseline_name')[merge_cols + [plot_column_name, ]]
        merged = baseline.merge(df.query('guidance_strength != @baseline_name'), on=merge_cols, suffixes=["_baseline", ""])
        merged["delta_" + plot_column_name] = merged[plot_column_name] - merged[plot_column_name + "_baseline"]
        return merged

    def add_simplified_metadata_columns(df):
        if "selection" in df.columns:
            df["selection_name"] = df.selection.apply(lambda x: 'all altloc' if 'chain_id' in x else x)
        df["guidance_weight"] = df.guidance_weight.astype(float)
        df["guidance_strength"] = df.guidance_weight.apply(lambda x: {
            0.5: "strong", 0.05: "strong", 0.2: "medium", 0.02: "medium", 0.1: "weak", 0.01: "weak", 0.0: "baseline"}[x])
        return df

    def add_tradeoff_guidance_regions(ax, df, x_col, y_col, guidance_order, guidance_palette):
        """Overlay smooth guidance-specific KDE regions on a tradeoff axis.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axis to draw on.
        df : pandas.DataFrame
            Data containing `guidance_strength`, `x_col`, and `y_col`.
        x_col : str
            Column name for the x-axis variable.
        y_col : str
            Column name for the y-axis variable.
        guidance_order : list[str]
            Guidance strengths to draw in order.
        guidance_palette : dict[str, str]
            Mapping from guidance strength to color.
        """
        for _guidance_strength in guidance_order:
            _subset = df.loc[
                df["guidance_strength"] == _guidance_strength,
                [x_col, y_col],
            ].replace([np.inf, -np.inf], np.nan).dropna()
            if (
                len(_subset) < 8
                or _subset[x_col].nunique() < 3
                or _subset[y_col].nunique() < 3
            ):
                continue
            try:
                sns.kdeplot(
                    data=_subset,
                    x=x_col,
                    y=y_col,
                    ax=ax,
                    fill=True,
                    levels=4,
                    thresh=0.18,
                    bw_adjust=1.05,
                    cut=1.6,
                    color=guidance_palette[_guidance_strength],
                    alpha=0.10,
                    warn_singular=False,
                    zorder=0,
                )
                sns.kdeplot(
                    data=_subset,
                    x=x_col,
                    y=y_col,
                    ax=ax,
                    fill=False,
                    levels=2,
                    thresh=0.24,
                    bw_adjust=1.05,
                    cut=1.6,
                    color=guidance_palette[_guidance_strength],
                    alpha=0.35,
                    linewidths=0.8,
                    warn_singular=False,
                    zorder=0.2,
                )
            except Exception:
                continue

    _figure_dir = Path("figures")
    _figure_dir.mkdir(parents=True, exist_ok=True)

    def save_plot_png(plot_obj, filename, dpi=220):
        figure = plot_obj.figure if hasattr(plot_obj, "figure") else plot_obj
        output_path = _figure_dir / filename
        output_path.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(output_path, dpi=dpi, bbox_inches="tight")
        return output_path

    return (
        add_simplified_metadata_columns,
        add_tradeoff_guidance_regions,
        extract_metadata_from_path,
        make_model_comparison_delta_frame,
        make_model_comparison_plot,
        pd,
        plt,
        save_plot_png,
        sns,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Clashes

    Let's take a look first at the clashes detected using phenix. We'll look at clashscores first.
    """)
    return


@app.cell
def _(pd):
    clashes = pd.read_csv("initial40results/clashscore_w_baseline.csv")
    clashes = clashes.rename(columns={"model_name": "filename"})
    return (clashes,)


@app.cell
def _(clashes):
    clashes.loc[0,'filename']
    return


@app.cell
def _(
    add_simplified_metadata_columns,
    clashes,
    extract_metadata_from_path,
    pd,
):
    clashes_1 = pd.merge(clashes, clashes.filename.apply(extract_metadata_from_path), left_index=True, right_index=True)
    clashes_1['selection'] = 'all residues'
    clashes_1['guidance_weight'] = clashes_1.guidance_weight.astype(float)
    clashes_1 = add_simplified_metadata_columns(clashes_1)
    clashes_1['protein'] = clashes_1.pdb_id.str.lower()
    return (clashes_1,)


@app.cell
def _(clashes_1):
    clash_summary = clashes_1.groupby(['model_name', 'pdb_id', 'guidance_weight']).agg({'num_clashes': 'mean', 'clashscore': 'mean'})
    return (clash_summary,)


@app.cell
def _(clash_summary):
    clash_summary
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Make sure all proteins are represented well across the models (there were certainly errors for several)

    We want to be careful that our results aren't biased by any missing results.
    """)
    return


@app.cell
def _(clashes_1):
    clashes_1.groupby('pdb_id').agg({'model_name': 'count'})
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    First let's see whether all the models are available (all should have 8, they do):
    """)
    return


@app.cell
def _(clashes_1):
    (clashes_1.groupby(['model_name', 'pdb_id', 'guidance_weight']).model_id.count() == 8).all()
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
def _(clashes_1):
    clashes_1.query('model_id == 1').pivot_table(index='pdb_id', columns=['model_name', 'guidance_weight'], values='clashscore', aggfunc='count')
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
def _(clashes_1, make_model_comparison_plot, save_plot_png):
    _g = make_model_comparison_plot(
        clashes_1.query("guidance_strength != 'strong' and selection_name != 'all altloc'"),
        "clashscore",
    )
    save_plot_png(_g, "01_clashscore_distribution.png")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Generally it looks like RF3 is having a harder time that the others with respect to clashes.

    Kind of hard to see trends with guidance, so let's look at the difference between each guidance level and the baseline.
    """)
    return


@app.cell
def _(
    clashes_1,
    make_model_comparison_delta_frame,
    make_model_comparison_plot,
    save_plot_png,
):
    clashes_diff_to_baseline = make_model_comparison_delta_frame(clashes_1, 'clashscore')
    clashes_diff_to_baseline['protein'] = clashes_diff_to_baseline.pdb_id
    g = make_model_comparison_plot(clashes_diff_to_baseline.query('delta_clashscore < 100'), 'delta_clashscore', is_delta=True)
    save_plot_png(g, "02_delta_clashscore_distribution.png")
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
    by_pdb_id = (clash_summary.reset_index()
     .query("guidance_weight == 0.1 or (model_name == 'rf3' and guidance_weight == 0.01)")
     .pivot(index="pdb_id", columns="model_name", values=['num_clashes', 'clashscore'])
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
def _(add_simplified_metadata_columns, extract_metadata_from_path, pd):
    bond_angle_stats = pd.read_csv("initial40results/bond_angle_violation_fractions.csv")

    bond_length_stats = pd.read_csv("initial40results/bond_length_violation_fractions.csv")

    bond_stats = bond_angle_stats.merge(bond_length_stats, on=["model", "model_n"], suffixes=["_bond_angles", "_bond_lengths"])
    bond_stats.rename(columns={"model": "filename"}, inplace=True)
    bond_stats = pd.merge(bond_stats, bond_stats.filename.apply(extract_metadata_from_path), left_index=True, right_index=True)
    bond_stats = add_simplified_metadata_columns(bond_stats)
    return (bond_stats,)


@app.cell
def _(bond_stats):
    bond_stats.groupby("pdb_id").agg({"model_name": "count"})  # check everything is there. Yep!
    return


@app.cell
def _(bond_stats):
    bond_stats_summary = bond_stats.groupby(["model_name", "pdb_id", "guidance_strength"]).agg(
        {"outlier_fraction_bond_angles": 'mean', "outlier_fraction_bond_lengths": "mean"}
    )
    return (bond_stats_summary,)


@app.cell
def _(bond_stats_summary):
    bond_stats_summary.groupby(["model_name", "guidance_strength"]).describe(percentiles=[0.5])
    return


@app.cell
def _(
    bond_stats,
    make_model_comparison_delta_frame,
    make_model_comparison_plot,
    save_plot_png,
):
    bond_stats['protein'] = bond_stats.pdb_id.str.lower()
    bond_angles_delta = make_model_comparison_delta_frame(bond_stats, 'outlier_fraction_bond_angles')
    g_1 = make_model_comparison_plot(bond_angles_delta.query('delta_outlier_fraction_bond_angles < 0.05'), 'delta_outlier_fraction_bond_angles', is_delta=True)
    save_plot_png(g_1, "03_delta_bond_angle_outlier_fraction_distribution.png")
    print(f'total bond angles in dataframe: {len(bond_angles_delta)}')
    # print the outliers
    bond_angles_delta.query('delta_outlier_fraction_bond_angles >= 0.05').groupby(['protein', 'model_name', 'guidance_strength']).agg({'outlier_fraction_bond_angles': ['mean', 'size'], 'delta_outlier_fraction_bond_angles': ['mean', 'max']})
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
    bonds_by_pdb_id = (bond_stats_summary.reset_index()
     .query("guidance_strength == 'weak'")
     .pivot(index="pdb_id", columns="model_name", values=['outlier_fraction_bond_angles', 'outlier_fraction_bond_lengths'])
    )
    bonds_by_pdb_id
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Backbone and side chain dihedrals
    """)
    return


@app.cell
def _(add_simplified_metadata_columns, extract_metadata_from_path, pd):
    dihedrals = pd.read_csv("initial40results/initial_40_tortoize_protein_stats.csv")
    dihedrals.rename(columns={"path": "filename", "model": "model_n"}, inplace=True)
    dihedrals = pd.merge(dihedrals, dihedrals.filename.apply(extract_metadata_from_path), left_index=True, right_index=True)
    dihedrals = add_simplified_metadata_columns(dihedrals)  # converts guidance_weight to float
    dihedrals["protein"] = dihedrals.pdb_id.str.lower()
    dihedrals
    return (dihedrals,)


@app.cell
def _(dihedrals):
    dihedrals_summary = dihedrals.groupby(["model_name", "pdb_id", "guidance_weight"]).agg(
        {"ramachandran_z_score": 'mean', "torsion_z_score": "mean"}
    )
    dihedrals_summary.groupby(["model_name", "guidance_weight"]).describe(percentiles=[0.1, 0.5, 0.9])
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
    dihedrals_by_pdb_id = (dihedrals_summary.reset_index()
     .query("guidance_weight == 0.1 or (model_name == 'rf3' and guidance_weight == 0.01)")
     .pivot(index="pdb_id", columns="model_name", values=['ramachandran_z_score', 'torsion_z_score'])
    )
    dihedrals_by_pdb_id
    return


@app.cell
def _(dihedrals_summary, pd, plt, save_plot_png, sns):
    dsw = dihedrals_summary.reset_index().query("guidance_weight == 0.1 or (model_name == 'rf3' and guidance_weight == 0.01)").dropna()
    dsw['guidance'] = 'weak'
    dss = dihedrals_summary.reset_index().query("guidance_weight == 0.5 or (model_name == 'rf3' and guidance_weight == 0.05)").dropna()
    dss['guidance'] = 'strong'
    dsb = dihedrals_summary.reset_index().query("guidance_weight == 0.0").dropna()
    dsb['guidance'] = 'baseline'
    ds = pd.concat([dsw, dss, dsb])
    _hue_order = ['weak', 'strong', 'baseline']
    _palette = {'weak': '#377EB8', 'strong': '#E41A1C', 'baseline': '#888888'}
    g_2 = sns.FacetGrid(ds, row='model_name', hue='guidance', hue_order=_hue_order, aspect=5, height=1.5, palette=_palette)
    g_2.map(sns.kdeplot, 'ramachandran_z_score', bw_adjust=0.5, clip_on=False, fill=True, alpha=0.75, linewidth=1.5)
    g_2.map(sns.kdeplot, 'ramachandran_z_score', clip_on=False, color='w', lw=2, bw_adjust=0.5)
    g_2.refline(y=0, linewidth=2, linestyle='-', color=None, clip_on=False)

    def label(*x, color=None, label=None):
        text = x[0].iloc[0]
        ax = plt.gca()
        ax.text(0, 0.2, text, fontweight='bold', color='black', ha='left', va='center', transform=ax.transAxes)
        ax.text(0, 0.4 if label == 'strong' else 0.6, label, fontweight='bold', color=color, ha='left', va='center', transform=ax.transAxes)
    g_2.map(label, 'model_name')
    g_2.set(xlabel='Ramachandran Z-Score')
    g_2.set_titles('')
    g_2.set(yticks=[0.0, 0.2, 0.4], ylabel='')
    g_2.set(xticks=[-9, -6, -3, 0, 3, 6])
    g_2.despine(bottom=True, left=True)
    save_plot_png(g_2, "04_ramachandran_zscore_kde_by_model.png")
    return (ds,)


@app.cell
def _(ds, plt, save_plot_png, sns):
    _hue_order = ['weak', 'strong']
    _palette = {'weak': '#377EB8', 'strong': '#E41A1C'}
    g_3 = sns.FacetGrid(ds, row='model_name', hue='guidance', hue_order=_hue_order, aspect=5, height=1.5, palette=_palette)
    g_3.map(sns.kdeplot, 'torsion_z_score', bw_adjust=0.5, clip_on=False, fill=True, alpha=0.75, linewidth=1.5)
    g_3.map(sns.kdeplot, 'torsion_z_score', clip_on=False, color='w', lw=2, bw_adjust=0.5)
    g_3.refline(y=0, linewidth=2, linestyle='-', color=None, clip_on=False)

    def label_1(x, color, label):
        text = x.iloc[0]
        ax = plt.gca()
        ax.text(0, 0.2, text, fontweight='bold', color='black', ha='left', va='center', transform=ax.transAxes)
        ax.text(0, 0.4 if label == 'strong' else 0.6, label, fontweight='bold', color=color, ha='left', va='center', transform=ax.transAxes)
    g_3.map(label_1, 'model_name')
    g_3.set(xlabel='Torsion Z-Score')
    g_3.set_titles('')
    g_3.set(yticks=[0.0, 0.2, 0.4], ylabel='')
    g_3.set(xticks=[-9, -6, -3, 0, 3, 6])
    g_3.despine(bottom=True, left=True)
    save_plot_png(g_3, "05_torsion_zscore_kde_by_model.png")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Torsion outlier counts

    Super weird. Let's take a look at the actual number of outliers. Turns out they are quite small.
    """)
    return


@app.cell
def _(add_simplified_metadata_columns, extract_metadata_from_path, pd):
    # this is a very large file for this kind of thing: almost 200MB
    torsion_res = pd.read_csv("initial40results/initial_40_tortoize_residues.csv")
    torsion_res.rename(columns={"path": "filename", "model": "model_n"}, inplace=True)
    torsion_res = pd.merge(torsion_res, torsion_res.filename.apply(extract_metadata_from_path), left_index=True, right_index=True)
    torsion_res = add_simplified_metadata_columns(torsion_res)
    torsion_res
    return (torsion_res,)


@app.cell
def _(torsion_res):
    torsion_res.query("pdb_id == '1VME'").groupby(["model_name", "guidance_weight"])[
        ['ramachandran_z_score', 'torsion_z_score']
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
    torsion_outliers = torsion_res.groupby(["model_n", "model_name", "guidance_weight", "pdb_id"]).agg(
        rama_outlier_frac=("ramachandran_z_score", lambda x: (abs(x) > 3).sum() / len(x)),
        torsion_outlier_frac=("torsion_z_score", lambda x: (abs(x) > 3).sum() / len(x)),
    )
    torsion_outliers
    return (torsion_outliers,)


@app.cell
def _(torsion_outliers):
    torsion_summary = torsion_outliers.groupby(["model_name", "pdb_id", "guidance_weight"]).agg(
        {"rama_outlier_frac": 'mean', "torsion_outlier_frac": "mean"}
    )
    torsion_summary.groupby(["model_name", "guidance_weight"]).describe(percentiles=[0.1, 0.5, 0.9])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Well, I'm not sure what to make of this. The Ramachandran outliers seem to be fewer (not necessarily better) with more guidance. This appears to be mostly true for torsion (side chain X_1, X_2) angles. But it is clear that the whole-protein scores get _much_ worse, if I'm correctly interpreting the whole protein scores. This implies that while there are fewer outliers with higher guidance, the outliers there are must get much worse. Let's see what the distributions look like:
    """)
    return


@app.cell
def _(torsion_res):
    dihedral_means = torsion_res.groupby(["model_name", "guidance_weight", "pdb_id"]).agg({
        'ramachandran_z_score': 'mean', 'torsion_z_score': 'mean'
    })
    return (dihedral_means,)


@app.cell
def _(dihedral_means, pd, plt, save_plot_png, sns):
    dsw_1 = dihedral_means.reset_index().query("guidance_weight == 0.1 or (model_name == 'rf3' and guidance_weight == 0.01)").dropna()
    dsw_1['guidance'] = 'weak'
    dss_1 = dihedral_means.reset_index().query("guidance_weight == 0.5 or (model_name == 'rf3' and guidance_weight == 0.05)").dropna()
    dss_1['guidance'] = 'strong'
    ds_1 = pd.concat([dsw_1, dss_1])
    _hue_order = ['weak', 'strong']
    _palette = {'weak': '#377EB8', 'strong': '#E41A1C'}
    g_4 = sns.FacetGrid(ds_1, row='model_name', hue='guidance', hue_order=_hue_order, aspect=5, height=1.5, palette=_palette)
    g_4.map(sns.kdeplot, 'ramachandran_z_score', bw_adjust=0.5, clip_on=False, fill=True, alpha=0.75, linewidth=1.5)
    g_4.map(sns.kdeplot, 'ramachandran_z_score', clip_on=False, color='w', lw=2, bw_adjust=0.5)
    g_4.refline(y=0, linewidth=2, linestyle='-', color=None, clip_on=False)

    def label_2(x, color, label):
        text = x.iloc[0]
        ax = plt.gca()
        ax.text(0, 0.2, text, fontweight='bold', color='black', ha='left', va='center', transform=ax.transAxes)
        ax.text(0, 0.4 if label == 'strong' else 0.6, label, fontweight='bold', color=color, ha='left', va='center', transform=ax.transAxes)
    g_4.map(label_2, 'model_name')
    g_4.set(xlabel='Ramachandran Z-Score')
    g_4.set_titles('')
    g_4.set(yticks=[0, 2, 4], ylabel='')
    g_4.despine(bottom=True, left=True)
    save_plot_png(g_4, "06_ramachandran_zscore_residue_mean_kde.png")
    return (ds_1,)


@app.cell
def _(ds_1, plt, save_plot_png, sns):
    _hue_order = ['weak', 'strong']
    _palette = {'weak': '#377EB8', 'strong': '#E41A1C'}
    g_5 = sns.FacetGrid(ds_1, row='model_name', hue='guidance', hue_order=_hue_order, aspect=5, height=1.5, palette=_palette)
    g_5.map(sns.kdeplot, 'torsion_z_score', bw_adjust=0.5, clip_on=False, fill=True, alpha=0.75, linewidth=1.5)
    g_5.map(sns.kdeplot, 'torsion_z_score', clip_on=False, color='w', lw=2, bw_adjust=0.5)
    g_5.refline(y=0, linewidth=2, linestyle='-', color=None, clip_on=False)

    def label_3(x, color, label):
        text = x.iloc[0]
        ax = plt.gca()
        ax.text(0, 0.2, text, fontweight='bold', color='black', ha='left', va='center', transform=ax.transAxes)
        ax.text(0, 0.4 if label == 'strong' else 0.6, label, fontweight='bold', color=color, ha='left', va='center', transform=ax.transAxes)
    g_5.map(label_3, 'model_name')
    g_5.set(xlabel='Torsion Z-Score')
    g_5.set_titles('')
    g_5.set(yticks=[0, 2, 4], ylabel='')
    g_5.despine(bottom=True, left=True)
    save_plot_png(g_5, "07_torsion_zscore_residue_mean_kde.png")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # RSCC results
    """)
    return


@app.cell
def _(add_simplified_metadata_columns, pd):
    rscc = pd.read_csv("initial40results/rscc_results_aligned_260302.csv")

    rscc["model_name"] = rscc.apply(lambda x: x["model"] + x["method"] if x["model"] == "boltz2" else x["model"], axis=1)
    rscc = rscc[["protein", "model_name", "ensemble_size", "guidance_weight", "rscc", "selection"]]
    rscc = add_simplified_metadata_columns(rscc)
    return (rscc,)


@app.cell
def _(add_simplified_metadata_columns, rscc):
    rscc_all = rscc.groupby(["protein", "model_name", "guidance_weight"]).last().reset_index()
    rscc_all = add_simplified_metadata_columns(rscc_all)
    rscc_all
    return (rscc_all,)


@app.cell
def _(make_model_comparison_plot, rscc_all, save_plot_png):
    _g = make_model_comparison_plot(rscc_all.query("guidance_strength in ('baseline', 'weak')"), "rscc")
    save_plot_png(_g, "08_rscc_baseline_vs_weak_distribution.png")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    what's with all the really low RSCC structures? Let's list some:
    """)
    return


@app.cell
def _(pd, rscc_all):
    # rscc_all.query("guidance_strength == 'weak' and rscc < 0.15")
    boltz2_msa_depth = pd.read_csv("initial40results/boltz2_msa_depth.csv")
    boltz2_msa_depth.pdb_id = boltz2_msa_depth.pdb_id.str.lower()
    w_msa_depth = rscc_all.merge(boltz2_msa_depth, left_on="protein", right_on="pdb_id")

    protein_stats = pd.read_csv("initial40results/Sampleworks_initial_structures.csv")
    protein_stats.PDB_ID = protein_stats.PDB_ID.str.lower()
    w_all_stats = w_msa_depth.merge(protein_stats, left_on="pdb_id", right_on="PDB_ID")
    w_all_stats.groupby("model_name")[[
        'rscc', 'msa_depth', "Resolution", "Rfree", 
        "Cluster_Size", "Num_Structure_100_identity"
    ]].corr().reset_index().query("level_1 == 'rscc'")
    return


@app.cell
def _(
    make_model_comparison_delta_frame,
    make_model_comparison_plot,
    rscc_all,
    save_plot_png,
):
    rscc_delta = make_model_comparison_delta_frame(rscc_all, "rscc")
    _g = make_model_comparison_plot(
        rscc_delta.query("delta_rscc < 0.3 and delta_rscc > -0.3"),
        "delta_rscc",
        is_delta=True,
    )
    save_plot_png(_g, "09_delta_rscc_distribution.png")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Do any zero-guidance models achieve any meaningful RSCC?
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Do any zero-guidance models achieve any meaningful RSCC?
    """)
    return


@app.cell
def _(rscc):
    rscc.query("guidance_weight == 0.0").pivot(index=["protein", "selection_name"], columns="model_name", values="rscc").describe()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    As we can see, the scores are mostly 0.8-0.9, which is not bad exactly but we would expect better. Boltz and Protenix RSCC fall apart at strong (0.5) guidance, but RF3 does not. Unfortunately, this is still not very discriminating. However, from the metrics above, we know that RF3 had the worst clash scores. It appears to be favoring RSCC over maintaining proper side-chain geometry. Let's see some examples:
    """)
    return


@app.cell
def _(rscc):
    rscc.query("protein == '5i09'")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # LDDT
    """)
    return


@app.cell
def _(pd):
    lddt_results = pd.read_csv("initial40results/lddt_results_w_baseline.csv")
    lddt_results["model_name"] = lddt_results.apply(lambda x: x["model"] + x["method"] if x["model"] == "boltz2" else x["model"], axis=1)
    lddt_results = lddt_results[["protein", "model_name", "guidance_weight", "selection", "occupancies", "avg_silhouette", "avg_silhouette_to_ref"]]
    lddt_results.head()
    return (lddt_results,)


@app.cell
def _(lddt_results):
    lddt_results.query("protein == '6b8x'")
    return


@app.cell
def _(add_simplified_metadata_columns, lddt_results):
    lddt_results_1 = add_simplified_metadata_columns(lddt_results)
    return (lddt_results_1,)


@app.cell
def _(lddt_results_1, make_model_comparison_plot, save_plot_png):
    _g = make_model_comparison_plot(
        lddt_results_1.query("guidance_strength == 'baseline' and selection_name != 'all altloc'"),
        "avg_silhouette_to_ref",
        smoothing=0.25,
    )
    save_plot_png(_g, "10_lddt_to_ref_baseline_distribution.png")
    return


@app.cell
def _(lddt_results_1, make_model_comparison_plot, save_plot_png):
    _g = make_model_comparison_plot(
        lddt_results_1.query("selection_name != 'all altloc'"),
        "avg_silhouette_to_ref",
        smoothing=0.25,
    )
    save_plot_png(_g, "11_lddt_to_ref_all_guidance_distribution.png")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Yikes
    """)
    return


@app.cell
def _(lddt_results_1):
    lddt_by_pdb_id = lddt_results_1.query("guidance_strength == 'weak' and selection_name != 'all altloc'").pivot(index=['protein', 'selection_name'], columns='model_name', values=['avg_silhouette_to_ref', 'avg_silhouette'])
    lddt_by_pdb_id
    return (lddt_by_pdb_id,)


@app.cell
def _(lddt_by_pdb_id):
    lddt_by_pdb_id.describe()
    return


@app.cell
def _(lddt_by_pdb_id):
    lddt_by_pdb_id.avg_silhouette_to_ref.idxmax(axis=1).reset_index(name="model_name").groupby("model_name").size()
    return


@app.cell
def _(lddt_by_pdb_id):
    pdb_id, sel = lddt_by_pdb_id[("avg_silhouette_to_ref", "boltz2X-RAY")].sort_values().dropna().index[-1]
    lddt_by_pdb_id.loc[(pdb_id, sel)]
    return


@app.cell
def _(lddt_by_pdb_id):
    lddt_by_pdb_id.query('protein in ("5sop", "8z76", "7z0e", "5rbr", "7frd")')
    return


@app.cell
def _(lddt_by_pdb_id):
    (lddt_by_pdb_id.avg_silhouette_to_ref.max(axis=1) - lddt_by_pdb_id.avg_silhouette_to_ref.min(axis=1)).sort_values()
    return


@app.cell
def _(lddt_results_1):
    # we'd like to know what is the change in LDDT statistics relative to no guidance.
    lddt_baseline = lddt_results_1.query('guidance_strength == "baseline"')[['protein', 'model_name', 'selection_name', 'avg_silhouette_to_ref']]
    lddt_merged = lddt_baseline.merge(lddt_results_1.query('guidance_strength != "baseline"'), on=['protein', 'model_name', 'selection_name'], suffixes=['_baseline', ''])
    lddt_merged['delta_silhouette_to_ref'] = lddt_merged.avg_silhouette_to_ref - lddt_merged.avg_silhouette_to_ref_baseline
    lddt_merged
    return (lddt_merged,)


@app.cell
def _(lddt_merged, make_model_comparison_plot, save_plot_png):
    _g = make_model_comparison_plot(lddt_merged, "delta_silhouette_to_ref", is_delta=True)
    save_plot_png(_g, "12_delta_lddt_to_ref_distribution.png")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Geometry vs fit trade-offs (across guidance strengths)

    We combine the baseline-relative geometry metrics we can measure directly
    into a simple **geometry penalty** score, then compare that penalty against
    ΔRSCC and ΔLDDT-to-reference.

    The simple geometry score averages Δ clashscore and the percentage-point
    changes in bond-length and bond-angle outlier fractions. Tortoize terms are
    shown separately in a companion plot.
    """)
    return


@app.cell
def _(
    bond_stats,
    clashes_1,
    dihedrals,
    lddt_merged,
    make_model_comparison_delta_frame,
    pd,
    rscc_all,
):
    model_name_alias = {
        "boltz2_X-RAY_DIFFRACTION": "boltz2X-RAY",
        "boltz2_MD": "boltz2MD",
        "boltz2xray": "boltz2X-RAY",
        "boltz2md": "boltz2MD",
    }

    def robust_z_score(series):
        _median = series.median()
        _mad = (series - _median).abs().median()
        _scale = 1.4826 * _mad + 1e-9
        return (series - _median) / _scale

    clash_stats = (
        clashes_1.groupby(
            ["protein", "model_name", "guidance_weight", "guidance_strength"],
            as_index=False,
        )
        .agg(clashscore=("clashscore", "mean"))
    )
    bond_stats_agg = (
        bond_stats.assign(protein=bond_stats.pdb_id.str.lower())
        .groupby(
            ["protein", "model_name", "guidance_weight", "guidance_strength"],
            as_index=False,
        )
        .agg(
            outlier_fraction_bond_angles=("outlier_fraction_bond_angles", "mean"),
            outlier_fraction_bond_lengths=("outlier_fraction_bond_lengths", "mean"),
        )
    )
    rscc_stats = (
        rscc_all.groupby(
            ["protein", "model_name", "guidance_weight", "guidance_strength"],
            as_index=False,
        )
        .agg(rscc=("rscc", "mean"))
    )
    dihedral_stats = (
        dihedrals.groupby(
            ["protein", "model_name", "guidance_weight", "guidance_strength"],
            as_index=False,
        )
        .agg(
            ramachandran_z_score=("ramachandran_z_score", "mean"),
            torsion_z_score=("torsion_z_score", "mean"),
        )
        .dropna(subset=["ramachandran_z_score", "torsion_z_score"])
    )
    dihedral_stats["abs_ramachandran_z_score"] = dihedral_stats["ramachandran_z_score"].abs()
    dihedral_stats["abs_torsion_z_score"] = dihedral_stats["torsion_z_score"].abs()

    clash_delta = make_model_comparison_delta_frame(clash_stats, "clashscore")[
        ["protein", "model_name", "guidance_weight", "guidance_strength", "delta_clashscore"]
    ]
    bond_angle_delta = make_model_comparison_delta_frame(
        bond_stats_agg, "outlier_fraction_bond_angles"
    )[
        [
            "protein",
            "model_name",
            "guidance_weight",
            "guidance_strength",
            "delta_outlier_fraction_bond_angles",
        ]
    ]
    bond_length_delta = make_model_comparison_delta_frame(
        bond_stats_agg, "outlier_fraction_bond_lengths"
    )[
        [
            "protein",
            "model_name",
            "guidance_weight",
            "guidance_strength",
            "delta_outlier_fraction_bond_lengths",
        ]
    ]
    rama_z_term = dihedral_stats[
        [
            "protein",
            "model_name",
            "guidance_weight",
            "guidance_strength",
            "abs_ramachandran_z_score",
        ]
    ]
    torsion_z_term = dihedral_stats[
        [
            "protein",
            "model_name",
            "guidance_weight",
            "guidance_strength",
            "abs_torsion_z_score",
        ]
    ]
    rscc_delta_tradeoff = make_model_comparison_delta_frame(rscc_stats, "rscc")[
        ["protein", "model_name", "guidance_weight", "guidance_strength", "delta_rscc"]
    ]
    lddt_delta = (
        lddt_merged.query("selection_name != 'all altloc'")
        .groupby(
            ["protein", "model_name", "guidance_weight", "guidance_strength"],
            as_index=False,
        )
        .agg(delta_silhouette_to_ref=("delta_silhouette_to_ref", "mean"))
    )

    for _df in [
        clash_delta,
        bond_angle_delta,
        bond_length_delta,
        rama_z_term,
        torsion_z_term,
        rscc_delta_tradeoff,
        lddt_delta,
    ]:
        _df["model_name"] = _df["model_name"].map(lambda x: model_name_alias.get(x, x))

    tradeoff_df = (
        clash_delta.merge(
            bond_angle_delta,
            on=["protein", "model_name", "guidance_weight", "guidance_strength"],
        )
        .merge(
            bond_length_delta,
            on=["protein", "model_name", "guidance_weight", "guidance_strength"],
        )
        .merge(
            rama_z_term,
            on=["protein", "model_name", "guidance_weight", "guidance_strength"],
        )
        .merge(
            torsion_z_term,
            on=["protein", "model_name", "guidance_weight", "guidance_strength"],
        )
        .merge(
            rscc_delta_tradeoff,
            on=["protein", "model_name", "guidance_weight", "guidance_strength"],
        )
        .merge(
            lddt_delta,
            on=["protein", "model_name", "guidance_weight", "guidance_strength"],
        )
    )

    tradeoff_df["guidance_order"] = tradeoff_df["guidance_strength"].map(
        {"weak": 1.0, "medium": 2.0, "strong": 3.0}
    )

    tradeoff_df["delta_outlier_fraction_bond_lengths_pct"] = (
        100.0 * tradeoff_df["delta_outlier_fraction_bond_lengths"]
    )
    tradeoff_df["delta_outlier_fraction_bond_angles_pct"] = (
        100.0 * tradeoff_df["delta_outlier_fraction_bond_angles"]
    )

    simple_geometry_terms = [
        "delta_clashscore",
        "delta_outlier_fraction_bond_lengths_pct",
        "delta_outlier_fraction_bond_angles_pct",
    ]
    geometry_weights = pd.DataFrame(
        {
            "geometry_term": [
                "delta_clashscore",
                "delta_outlier_fraction_bond_lengths_pct",
                "delta_outlier_fraction_bond_angles_pct",
            ],
            "display_name": [
                "Δ clashscore",
                "Δ bond-length outlier fraction (percentage points)",
                "Δ bond-angle outlier fraction (percentage points)",
            ],
            "weight": [1 / 3, 1 / 3, 1 / 3],
        }
    )

    tradeoff_df["geometry_penalty"] = tradeoff_df[simple_geometry_terms].mean(axis=1)
    tradeoff_df["tortoize_combined_score"] = 0.5 * (
        tradeoff_df["abs_ramachandran_z_score"]
        + tradeoff_df["abs_torsion_z_score"]
    )
    tradeoff_df["geometry_penalty_z_model"] = tradeoff_df.groupby("model_name")[
        "geometry_penalty"
    ].transform(robust_z_score).clip(-5, 5)

    for _col in ["delta_rscc", "delta_silhouette_to_ref"]:
        tradeoff_df[f"{_col}_z_model"] = tradeoff_df.groupby("model_name")[_col].transform(
            robust_z_score
        ).clip(-5, 5)

    tradeoff_df["balance_score"] = (
        tradeoff_df["delta_rscc_z_model"]
        + tradeoff_df["delta_silhouette_to_ref_z_model"]
        - tradeoff_df["geometry_penalty_z_model"]
    )

    weak_tradeoff = tradeoff_df.query("guidance_strength == 'weak'").copy()

    tradeoff_df["peer_model_count"] = tradeoff_df.groupby(
        ["protein", "guidance_strength"]
    )["model_name"].transform("nunique")
    _peer_balance_sum = tradeoff_df.groupby(["protein", "guidance_strength"])[
        "balance_score"
    ].transform("sum")
    tradeoff_df["peer_balance_mean_other"] = (
        (_peer_balance_sum - tradeoff_df["balance_score"])
        / (tradeoff_df["peer_model_count"] - 1)
    ).where(tradeoff_df["peer_model_count"] > 1)
    tradeoff_df["balance_vs_peers"] = (
        tradeoff_df["balance_score"] - tradeoff_df["peer_balance_mean_other"]
    )

    tradeoff_balance_best = (
        tradeoff_df.sort_values(["model_name", "balance_score"], ascending=[True, False])
        .groupby("model_name", group_keys=False)
        .head(1)
        .copy()
    )
    tradeoff_balance_best["highlight_kind"] = "best"
    tradeoff_balance_best["case_type"] = "Best balance"

    tradeoff_balance_worst = (
        tradeoff_df.sort_values(["model_name", "balance_score"], ascending=[True, True])
        .groupby("model_name", group_keys=False)
        .head(1)
        .copy()
    )
    tradeoff_balance_worst["highlight_kind"] = "worst"
    tradeoff_balance_worst["case_type"] = "Worst balance"

    tradeoff_peer_best = (
        tradeoff_df.query("peer_model_count > 1")
        .sort_values(["model_name", "balance_vs_peers"], ascending=[True, False])
        .groupby("model_name", group_keys=False)
        .head(1)
        .copy()
    )
    tradeoff_peer_best["highlight_kind"] = "peer+"
    tradeoff_peer_best["case_type"] = "Differentially better than peers"

    tradeoff_peer_worst = (
        tradeoff_df.query("peer_model_count > 1")
        .sort_values(["model_name", "balance_vs_peers"], ascending=[True, True])
        .groupby("model_name", group_keys=False)
        .head(1)
        .copy()
    )
    tradeoff_peer_worst["highlight_kind"] = "peer-"
    tradeoff_peer_worst["case_type"] = "Differentially worse than peers"

    tradeoff_case_summary = pd.concat(
        [
            tradeoff_balance_best,
            tradeoff_balance_worst,
            tradeoff_peer_best,
            tradeoff_peer_worst,
        ],
        ignore_index=True,
    )[
        [
            "protein",
            "model_name",
            "guidance_strength",
            "guidance_weight",
            "highlight_kind",
            "case_type",
            "delta_rscc",
            "delta_silhouette_to_ref",
            "geometry_penalty",
            "tortoize_combined_score",
            "balance_score",
            "balance_vs_peers",
        ]
    ].copy()
    tradeoff_case_summary["case_note"] = (
        tradeoff_case_summary["protein"].str.upper()
        + " ("
        + tradeoff_case_summary["guidance_strength"].str[0].str.upper()
        + ")"
    )

    tradeoff_case_table = (
        tradeoff_case_summary.pivot(
            index="model_name",
            columns="case_type",
            values="case_note",
        )
        .reindex(["boltz2X-RAY", "boltz2MD", "protenix", "rf3"])
        .reset_index()
    )

    tradeoff_highlight_points = (
        tradeoff_case_summary.groupby(
            [
                "protein",
                "model_name",
                "guidance_strength",
                "guidance_weight",
                "delta_rscc",
                "delta_silhouette_to_ref",
                "geometry_penalty",
                "tortoize_combined_score",
                "balance_score",
                "balance_vs_peers",
            ],
            as_index=False,
            dropna=False,
        )
        .agg(
            highlight_kinds=(
                "highlight_kind",
                lambda s: tuple(dict.fromkeys(s.tolist())),
            ),
            case_types=(
                "case_type",
                lambda s: tuple(dict.fromkeys(s.tolist())),
            ),
        )
    )
    tradeoff_highlight_points["is_balance_best"] = tradeoff_highlight_points[
        "highlight_kinds"
    ].map(lambda values: "best" in values)
    tradeoff_highlight_points["is_balance_worst"] = tradeoff_highlight_points[
        "highlight_kinds"
    ].map(lambda values: "worst" in values)
    tradeoff_highlight_points["is_peer_best"] = tradeoff_highlight_points[
        "highlight_kinds"
    ].map(lambda values: "peer+" in values)
    tradeoff_highlight_points["is_peer_worst"] = tradeoff_highlight_points[
        "highlight_kinds"
    ].map(lambda values: "peer-" in values)
    tradeoff_highlight_points["point_label"] = tradeoff_highlight_points.apply(
        lambda row: f"{row['protein'].upper()} ({row['guidance_strength'][0].upper()})",
        axis=1,
    )

    focus_5i09 = weak_tradeoff.query("protein == '5i09'")[
        [
            "protein",
            "model_name",
            "delta_rscc",
            "delta_silhouette_to_ref",
            "abs_ramachandran_z_score",
            "abs_torsion_z_score",
            "delta_clashscore",
            "delta_outlier_fraction_bond_angles",
            "delta_outlier_fraction_bond_lengths",
            "geometry_penalty",
            "balance_score",
        ]
    ].sort_values("model_name")

    well_balanced_examples = (
        weak_tradeoff.query(
            "delta_rscc > 0 and delta_silhouette_to_ref > 0 and geometry_penalty <= 0"
        )
        .sort_values(["model_name", "balance_score"], ascending=[True, False])
        .groupby("model_name")
        .head(3)[
            [
                "protein",
                "model_name",
                "delta_rscc",
                "delta_silhouette_to_ref",
                "geometry_penalty",
                "balance_score",
            ]
        ]
    )

    rscc_tradeoff_examples = (
        weak_tradeoff.query(
            "delta_rscc > 0 and (delta_silhouette_to_ref < 0 or geometry_penalty > 0)"
        )
        .sort_values(["model_name", "delta_rscc"], ascending=[True, False])
        .groupby("model_name")
        .head(3)[
            [
                "protein",
                "model_name",
                "delta_rscc",
                "delta_silhouette_to_ref",
                "geometry_penalty",
            ]
        ]
    )

    lddt_tradeoff_examples = (
        weak_tradeoff.query(
            "delta_silhouette_to_ref > 0 and (delta_rscc < 0 or geometry_penalty > 0)"
        )
        .sort_values(["model_name", "delta_silhouette_to_ref"], ascending=[True, False])
        .groupby("model_name")
        .head(3)[
            [
                "protein",
                "model_name",
                "delta_rscc",
                "delta_silhouette_to_ref",
                "geometry_penalty",
            ]
        ]
    )
    return (
        focus_5i09,
        geometry_weights,
        lddt_tradeoff_examples,
        rscc_tradeoff_examples,
        tradeoff_case_summary,
        tradeoff_case_table,
        tradeoff_df,
        tradeoff_highlight_points,
        weak_tradeoff,
        well_balanced_examples,
    )


@app.cell
def _(plt, sns):
    _CONCEPTUAL_MODEL_PALETTE = {
        "boltz2X-RAY": "#1B9E77",
        "boltz2MD": "#D95F02",
        "protenix": "#7570B3",
        "rf3": "#E7298A",
    }

    def compute_conceptual_tradeoff_limits(
        tradeoff_df, guidance_strengths=("weak", "strong")
    ):
        """Compute shared RSCC/geometry axis limits for conceptual figures.

        Parameters
        ----------
        tradeoff_df : pandas.DataFrame
            Tradeoff table containing `guidance_strength`, `delta_rscc`, and
            `geometry_penalty`.
        guidance_strengths : tuple[str, ...], default=("weak", "strong")
            Guidance buckets whose combined range should define the limits.

        Returns
        -------
        tuple[tuple[float, float], tuple[float, float]]
            Shared `(xlim, ylim)` pair for conceptual RSCC tradeoff plots.

        Raises
        ------
        ValueError
            If the requested subset has no finite rows.
        """
        _subset = (
            tradeoff_df[
                tradeoff_df["guidance_strength"].isin(guidance_strengths)
            ][["delta_rscc", "geometry_penalty"]]
            .replace([float("inf"), float("-inf")], float("nan"))
            .dropna()
        )
        if _subset.empty:
            raise ValueError("No finite rows available for conceptual limits.")

        _x_min = float(_subset["delta_rscc"].min())
        _x_max = float(_subset["delta_rscc"].max())
        _y_min = float(_subset["geometry_penalty"].min())
        _y_max = float(_subset["geometry_penalty"].max())

        _x_pad = max(0.08 * (_x_max - _x_min), 0.02)
        _y_pad = max(0.10 * (_y_max - _y_min), 0.50)

        return (
            (_x_min - _x_pad, _x_max + _x_pad),
            (_y_min - _y_pad, _y_max + _y_pad),
        )

    def make_conceptual_tradeoff_figure(
        tradeoff_df,
        guidance_strength,
        *,
        xlim,
        ylim,
        title,
        zoom_quantiles=None,
    ):
        """Create a simplified RSCC-versus-geometry tradeoff figure.

        Parameters
        ----------
        tradeoff_df : pandas.DataFrame
            Tradeoff table containing `model_name`, `guidance_strength`,
            `delta_rscc`, and `geometry_penalty`.
        guidance_strength : str
            Guidance bucket to visualize.
        xlim : tuple[float, float]
            X-axis limits to apply when `zoom_quantiles` is not used.
        ylim : tuple[float, float]
            Y-axis limits to apply when `zoom_quantiles` is not used.
        title : str
            Figure title.
        zoom_quantiles : tuple[float, float] | None, default=None
            Optional quantile range used to create a zoomed-in weak-guidance
            view. When provided, rows outside the padded quantile window are
            omitted from the plot and the axes are tightened to that window.

        Returns
        -------
        matplotlib.figure.Figure
            Figure with one shared RSCC-versus-geometry panel and no model
            labels.

        Raises
        ------
        ValueError
            If the requested guidance subset is empty.
        """

        def _add_model_regions(ax, df, x_col, y_col, model_order, model_palette):
            """Overlay model-specific KDE regions on a shared tradeoff axis.

            Parameters
            ----------
            ax : matplotlib.axes.Axes
                Axis to draw on.
            df : pandas.DataFrame
                Data containing `model_name`, `x_col`, and `y_col`.
            x_col : str
                Column name for the x-axis variable.
            y_col : str
                Column name for the y-axis variable.
            model_order : list[str]
                Model names to draw in order.
            model_palette : dict[str, str]
                Mapping from model name to color.
            """
            for _model in model_order:
                _subset = (
                    df.loc[df["model_name"] == _model, [x_col, y_col]]
                    .replace([float("inf"), float("-inf")], float("nan"))
                    .dropna()
                )
                if (
                    len(_subset) < 8
                    or _subset[x_col].nunique() < 3
                    or _subset[y_col].nunique() < 3
                ):
                    continue
                try:
                    sns.kdeplot(
                        data=_subset,
                        x=x_col,
                        y=y_col,
                        ax=ax,
                        fill=True,
                        levels=4,
                        thresh=0.16,
                        bw_adjust=1.0,
                        cut=1.6,
                        color=model_palette[_model],
                        alpha=0.09,
                        warn_singular=False,
                        zorder=0,
                    )
                    sns.kdeplot(
                        data=_subset,
                        x=x_col,
                        y=y_col,
                        ax=ax,
                        fill=False,
                        levels=2,
                        thresh=0.22,
                        bw_adjust=1.0,
                        cut=1.6,
                        color=model_palette[_model],
                        alpha=0.30,
                        linewidths=0.8,
                        warn_singular=False,
                        zorder=0.2,
                    )
                except Exception:
                    continue

        _guidance_df = tradeoff_df[
            tradeoff_df["guidance_strength"] == guidance_strength
        ].copy()
        if _guidance_df.empty:
            raise ValueError(
                f"No rows found for guidance_strength={guidance_strength!r}."
            )

        _plot_df = _guidance_df.copy()
        _plot_xlim = xlim
        _plot_ylim = ylim
        if zoom_quantiles is not None:
            _low_q, _high_q = zoom_quantiles
            _x_low, _x_high = _guidance_df["delta_rscc"].quantile([_low_q, _high_q])
            _y_low, _y_high = _guidance_df["geometry_penalty"].quantile([_low_q, _high_q])
            _x_pad = max(0.10 * (_x_high - _x_low), 0.01)
            _y_pad = max(0.16 * (_y_high - _y_low), 0.25)
            _plot_xlim = (float(_x_low - _x_pad), float(_x_high + _x_pad))
            _plot_ylim = (float(_y_low - _y_pad), float(_y_high + _y_pad))
            _plot_df = _guidance_df[
                _guidance_df["delta_rscc"].between(*_plot_xlim)
                & _guidance_df["geometry_penalty"].between(*_plot_ylim)
            ].copy()

        _model_order = [
            _model
            for _model in ["boltz2X-RAY", "boltz2MD", "protenix", "rf3"]
            if _model in _plot_df["model_name"].unique()
        ]
        if not _model_order:
            _model_order = sorted(_plot_df["model_name"].unique())

        _fig, _ax = plt.subplots(figsize=(5.8, 4.2))
        _add_model_regions(
            _ax,
            _plot_df,
            "delta_rscc",
            "geometry_penalty",
            _model_order,
            _CONCEPTUAL_MODEL_PALETTE,
        )

        for _model in _model_order:
            _model_df = _plot_df[_plot_df["model_name"] == _model]
            if _model_df.empty:
                continue

            _color = _CONCEPTUAL_MODEL_PALETTE.get(_model, "#4D4D4D")
            _ax.scatter(
                _model_df["delta_rscc"],
                _model_df["geometry_penalty"],
                s=20,
                c=_color,
                alpha=0.24,
                edgecolors="white",
                linewidths=0.25,
                zorder=1.2,
            )
            _ax.scatter(
                [_model_df["delta_rscc"].median()],
                [_model_df["geometry_penalty"].median()],
                s=95,
                c=_color,
                edgecolors="black",
                linewidths=0.9,
                zorder=4,
            )

        _ax.axvline(0, color="#555555", linestyle="--", linewidth=1.0, alpha=0.6)
        _ax.axhline(0, color="#555555", linestyle="--", linewidth=1.0, alpha=0.6)
        _ax.set_xlim(_plot_xlim)
        _ax.set_ylim(_plot_ylim)
        _ax.set_xlabel("Δ RSCC vs baseline")
        _ax.set_ylabel("Geometry penalty")
        _ax.set_title(title, fontsize=13)
        _fig.tight_layout()
        return _fig

    return (compute_conceptual_tradeoff_limits, make_conceptual_tradeoff_figure)


@app.cell
def _(add_tradeoff_guidance_regions, plt, save_plot_png, tradeoff_df):
    from matplotlib.lines import Line2D as _Line2D

    _model_display_tradeoff_plot = {
        "boltz2X-RAY": "Boltz2-XRD",
        "boltz2MD": "Boltz2-MD",
        "protenix": "Protenix",
        "rf3": "RF3",
    }
    _guidance_palette = {
        "weak": "#377EB8",
        "medium": "#FF7F00",
        "strong": "#E41A1C",
    }
    _guidance_markers = {
        "weak": "o",
        "medium": "s",
        "strong": "^",
    }

    _y_values = tradeoff_df["geometry_penalty"].replace(
        [float("inf"), float("-inf")], float("nan")
    ).dropna()
    _y_pad = max(0.16 * (_y_values.max() - _y_values.min()), 0.30)
    _y_lim = (float(_y_values.min() - _y_pad), float(_y_values.max() + _y_pad))

    _metric_rows = [
        ("delta_rscc", "Δ RSCC vs baseline"),
        ("delta_silhouette_to_ref", "Δ LDDT silhouette-to-ref vs baseline"),
    ]
    _model_order = ["boltz2X-RAY", "boltz2MD", "protenix", "rf3"]
    _guidance_order = ["weak", "medium", "strong"]

    _row_xlim = {}
    for _metric, _ in _metric_rows:
        _values = tradeoff_df[_metric].replace([float("inf"), float("-inf")], float("nan")).dropna()
        if _values.empty:
            _row_xlim[_metric] = None
            continue
        _max_abs = float(_values.abs().max())
        if _max_abs == 0:
            _max_abs = 0.01
        _pad = max(0.12 * _max_abs, 0.02)
        _row_xlim[_metric] = (-_max_abs - _pad, _max_abs + _pad)

    _fig, _axes = plt.subplots(
        2, len(_model_order), figsize=(17, 7.6), sharex="row", sharey=True
    )

    for _row_idx, (_metric, _xlabel) in enumerate(_metric_rows):
        for _col_idx, _model in enumerate(_model_order):
            _ax = _axes[_row_idx, _col_idx]
            _model_df = tradeoff_df[tradeoff_df["model_name"] == _model]

            if _model_df.empty:
                _ax.axis("off")
                continue

            add_tradeoff_guidance_regions(
                _ax,
                _model_df,
                _metric,
                "geometry_penalty",
                _guidance_order,
                _guidance_palette,
            )

            for _guidance_strength in _guidance_order:
                _guidance_df = _model_df[
                    _model_df["guidance_strength"] == _guidance_strength
                ]
                if _guidance_df.empty:
                    continue

                _ax.scatter(
                    _guidance_df[_metric],
                    _guidance_df["geometry_penalty"],
                    s=16,
                    c=_guidance_palette[_guidance_strength],
                    alpha=0.24,
                    marker=_guidance_markers[_guidance_strength],
                    edgecolors="white",
                    linewidths=0.2,
                    zorder=1.2,
                )

                _centroid_x = _guidance_df[_metric].median()
                _centroid_y = _guidance_df["geometry_penalty"].median()
                _ax.scatter(
                    [_centroid_x],
                    [_centroid_y],
                    s=70,
                    c=_guidance_palette[_guidance_strength],
                    marker=_guidance_markers[_guidance_strength],
                    edgecolors="black",
                    linewidths=0.8,
                    zorder=4,
                )

            _focus = _model_df[_model_df["protein"] == "5i09"].sort_values(
                "guidance_order"
            )
            if not _focus.empty:
                _ax.plot(
                    _focus[_metric],
                    _focus["geometry_penalty"],
                    color="black",
                    linewidth=1.0,
                    alpha=0.75,
                    zorder=5,
                )
                for _focus_row in _focus.itertuples(index=False):
                    _ax.scatter(
                        [getattr(_focus_row, _metric)],
                        [_focus_row.geometry_penalty],
                        s=70,
                        facecolors="white",
                        edgecolors="black",
                        linewidths=1.0,
                        marker=_guidance_markers.get(_focus_row.guidance_strength, "o"),
                        zorder=6,
                    )

            _ax.axvline(0, color="#555555", linestyle="--", linewidth=1.0, alpha=0.6)
            _ax.axhline(0, color="#555555", linestyle="--", linewidth=1.0, alpha=0.6)

            if _row_xlim[_metric] is not None:
                _ax.set_xlim(_row_xlim[_metric])
            _ax.set_ylim(_y_lim)

            if _row_idx == 0:
                _ax.set_title(_model_display_tradeoff_plot.get(_model, _model), fontsize=11)
            if _col_idx == 0:
                _ax.set_ylabel("Simple geometry score")
            _ax.set_xlabel(_xlabel)

    _legend_handles = [
        _Line2D(
            [0],
            [0],
            marker=_guidance_markers[_guidance],
            color="none",
            markerfacecolor=_guidance_palette[_guidance],
            markeredgecolor="black",
            markersize=8,
            label=f"{_guidance.capitalize()} guidance (median)",
        )
        for _guidance in _guidance_order
    ]
    _legend_handles.append(
        _Line2D(
            [0],
            [0],
            color="black",
            linewidth=1.0,
            label="5I09 trajectory",
        )
    )

    _fig.legend(
        handles=_legend_handles,
        loc="upper center",
        ncol=4,
        frameon=False,
        bbox_to_anchor=(0.5, 1.01),
    )
    _fig.suptitle(
        "Guidance trade-offs by model (simple geometry penalty; weak/medium/strong)",
        y=1.05,
        fontsize=13,
    )
    _fig.tight_layout()
    save_plot_png(_fig, "13_tradeoff_scatter_all_guidance.png")
    plt.show()
    _tradeoff_plot = _fig
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        "\n".join(
            [
                "### Conceptual overlaid RSCC tradeoff views",
                "",
                "These companion figures collapse the conceptual view to a single",
                "RSCC-versus-geometry axis, with color indicating model family and",
                "no model names shown in the plot itself.",
            ]
        )
    )
    return


@app.cell
def _(compute_conceptual_tradeoff_limits, tradeoff_df):
    conceptual_tradeoff_xlim, conceptual_tradeoff_ylim = (
        compute_conceptual_tradeoff_limits(tradeoff_df)
    )
    return conceptual_tradeoff_xlim, conceptual_tradeoff_ylim


@app.cell
def _(
    conceptual_tradeoff_xlim,
    conceptual_tradeoff_ylim,
    make_conceptual_tradeoff_figure,
    plt,
    save_plot_png,
    tradeoff_df,
):
    _conceptual_tradeoff_weak = make_conceptual_tradeoff_figure(
        tradeoff_df,
        "weak",
        xlim=conceptual_tradeoff_xlim,
        ylim=conceptual_tradeoff_ylim,
        title="Tradeoff Landscape (Weak Guidance)",
    )
    save_plot_png(_conceptual_tradeoff_weak, "16_tradeoff_conceptual_weak.png")
    plt.show()
    return


@app.cell
def _(
    conceptual_tradeoff_xlim,
    conceptual_tradeoff_ylim,
    make_conceptual_tradeoff_figure,
    plt,
    save_plot_png,
    tradeoff_df,
):
    _conceptual_tradeoff_weak_zoomed = make_conceptual_tradeoff_figure(
        tradeoff_df,
        "weak",
        xlim=conceptual_tradeoff_xlim,
        ylim=conceptual_tradeoff_ylim,
        title="Tradeoff Landscape (Weak Guidance, Zoomed)",
        zoom_quantiles=(0.02, 0.95),
    )
    save_plot_png(
        _conceptual_tradeoff_weak_zoomed,
        "17_tradeoff_conceptual_weak_zoomed.png",
    )
    plt.show()
    return


@app.cell
def _(
    conceptual_tradeoff_xlim,
    conceptual_tradeoff_ylim,
    make_conceptual_tradeoff_figure,
    plt,
    save_plot_png,
    tradeoff_df,
):
    _conceptual_tradeoff_strong = make_conceptual_tradeoff_figure(
        tradeoff_df,
        "strong",
        xlim=conceptual_tradeoff_xlim,
        ylim=conceptual_tradeoff_ylim,
        title="Tradeoff Landscape (Strong Guidance)",
    )
    save_plot_png(_conceptual_tradeoff_strong, "18_tradeoff_conceptual_strong.png")
    plt.show()
    return


@app.cell
def _(
    add_tradeoff_guidance_regions,
    plt,
    save_plot_png,
    tradeoff_df,
    tradeoff_highlight_points,
):
    from matplotlib.lines import Line2D as _Line2D

    _model_display_tradeoff_plot = {
        "boltz2X-RAY": "Boltz2-XRD",
        "boltz2MD": "Boltz2-MD",
        "protenix": "Protenix",
        "rf3": "RF3",
    }
    _guidance_markers = {
        "weak": "o",
        "medium": "s",
        "strong": "^",
    }
    _guidance_palette = {
        "weak": "#377EB8",
        "medium": "#FF7F00",
        "strong": "#E41A1C",
    }
    _highlight_styles = {
        "is_balance_best": {
            "color": "#1B9E77",
            "marker": "o",
            "label": "Best balance",
            "size": 180,
        },
        "is_balance_worst": {
            "color": "#CC3366",
            "marker": "o",
            "label": "Worst balance",
            "size": 180,
        },
        "is_peer_best": {
            "color": "#1F78B4",
            "marker": "D",
            "label": "Differentially better than peers",
            "size": 135,
        },
        "is_peer_worst": {
            "color": "#6A3D9A",
            "marker": "D",
            "label": "Differentially worse than peers",
            "size": 135,
        },
    }
    _label_templates = {
        ("right", "up"): [(8, 10), (8, 22), (8, -2), (8, 34)],
        ("right", "down"): [(8, -8), (8, -20), (8, 4), (8, -32)],
        ("left", "up"): [(-10, 10), (-10, 22), (-10, -2), (-10, 34)],
        ("left", "down"): [(-10, -8), (-10, -20), (-10, 4), (-10, -32)],
    }

    _y_values = tradeoff_df["geometry_penalty"].replace(
        [float("inf"), float("-inf")], float("nan")
    ).dropna()
    _y_pad = max(0.16 * (_y_values.max() - _y_values.min()), 0.30)
    _y_lim = (float(_y_values.min() - _y_pad), float(_y_values.max() + _y_pad))

    _metric_rows = [
        ("delta_rscc", "Δ RSCC vs baseline"),
        ("delta_silhouette_to_ref", "Δ LDDT silhouette-to-ref vs baseline"),
    ]
    _model_order = ["boltz2X-RAY", "boltz2MD", "protenix", "rf3"]
    _guidance_order = ["weak", "medium", "strong"]

    _row_xlim = {}
    for _metric, _ in _metric_rows:
        _values = tradeoff_df[_metric].replace([float("inf"), float("-inf")], float("nan")).dropna()
        if _values.empty:
            _row_xlim[_metric] = None
            continue
        _max_abs = float(_values.abs().max())
        if _max_abs == 0:
            _max_abs = 0.01
        _pad = max(0.12 * _max_abs, 0.02)
        _row_xlim[_metric] = (-_max_abs - _pad, _max_abs + _pad)

    _fig, _axes = plt.subplots(
        2, len(_model_order), figsize=(18, 8.2), sharex="row", sharey=True
    )

    for _row_idx, (_metric, _xlabel) in enumerate(_metric_rows):
        for _col_idx, _model in enumerate(_model_order):
            _ax = _axes[_row_idx, _col_idx]
            _model_df = tradeoff_df[tradeoff_df["model_name"] == _model]

            if _model_df.empty:
                _ax.axis("off")
                continue

            add_tradeoff_guidance_regions(
                _ax,
                _model_df,
                _metric,
                "geometry_penalty",
                _guidance_order,
                _guidance_palette,
            )

            for _guidance_strength in _guidance_order:
                _guidance_df = _model_df[
                    _model_df["guidance_strength"] == _guidance_strength
                ]
                if _guidance_df.empty:
                    continue
                _ax.scatter(
                    _guidance_df[_metric],
                    _guidance_df["geometry_penalty"],
                    s=14,
                    c=_guidance_palette[_guidance_strength],
                    alpha=0.22,
                    marker=_guidance_markers[_guidance_strength],
                    edgecolors="white",
                    linewidths=0.2,
                    zorder=1.2,
                )

            _highlight_df = tradeoff_highlight_points[
                tradeoff_highlight_points["model_name"] == _model
            ].copy()
            if not _highlight_df.empty:
                _highlight_df["primary_kind"] = _highlight_df["highlight_kinds"].map(
                    lambda values: values[0] if len(values) > 0 else "best"
                )
                _highlight_df = _highlight_df.sort_values(
                    ["geometry_penalty", _metric], ascending=[False, True]
                )

                for _flag, _style in _highlight_styles.items():
                    _subset = _highlight_df[_highlight_df[_flag]]
                    if _subset.empty:
                        continue
                    _ax.scatter(
                        _subset[_metric],
                        _subset["geometry_penalty"],
                        s=_style["size"],
                        facecolors="none",
                        edgecolors=_style["color"],
                        linewidths=2.0,
                        marker=_style["marker"],
                        zorder=4,
                    )

                for _guidance_strength in _guidance_order:
                    _highlight_guidance_df = _highlight_df[
                        _highlight_df["guidance_strength"] == _guidance_strength
                    ]
                    if _highlight_guidance_df.empty:
                        continue
                    _ax.scatter(
                        _highlight_guidance_df[_metric],
                        _highlight_guidance_df["geometry_penalty"],
                        s=42,
                        c="black",
                        marker=_guidance_markers[_guidance_strength],
                        zorder=5,
                    )

                if _row_idx == 0 and _row_xlim[_metric] is not None:
                    _x_min, _x_max = _row_xlim[_metric]
                    _y_min, _y_max = _y_lim
                    _highlight_df["x_frac"] = (
                        (_highlight_df[_metric] - _x_min) / (_x_max - _x_min)
                    )
                    _highlight_df["y_frac"] = (
                        (_highlight_df["geometry_penalty"] - _y_min) / (_y_max - _y_min)
                    )
                    _highlight_df["label_side"] = _highlight_df["x_frac"].map(
                        lambda value: (
                            "right"
                            if value < 0.22
                            else "left"
                            if value > 0.78
                            else "right"
                            if value < 0.5
                            else "left"
                        )
                    )
                    _highlight_df["label_vertical"] = _highlight_df["y_frac"].map(
                        lambda value: "down" if value > 0.72 else "up"
                    )
                    for (_label_side, _label_vertical), _side_df in _highlight_df.groupby(
                        ["label_side", "label_vertical"], sort=False
                    ):
                        _side_df = _side_df.sort_values(
                            ["geometry_penalty", _metric], ascending=[False, True]
                        ).reset_index(drop=True)
                        _side_offsets = _label_templates[(_label_side, _label_vertical)]
                        for _label_idx, _highlight_row in enumerate(
                            _side_df.itertuples(index=False)
                        ):
                            _offset = _side_offsets[min(_label_idx, len(_side_offsets) - 1)]
                            _ax.annotate(
                                _highlight_row.point_label,
                                xy=(
                                    getattr(_highlight_row, _metric),
                                    _highlight_row.geometry_penalty,
                                ),
                                xytext=_offset,
                                textcoords="offset points",
                                fontsize=6.4,
                                ha="left" if _label_side == "right" else "right",
                                va="center",
                                bbox={
                                    "boxstyle": "round,pad=0.14",
                                    "fc": "white",
                                    "ec": "#666666",
                                    "alpha": 0.92,
                                    "lw": 0.65,
                                },
                                arrowprops={
                                    "arrowstyle": "-",
                                    "color": "#666666",
                                    "lw": 0.55,
                                    "shrinkA": 0,
                                    "shrinkB": 4,
                                },
                                zorder=6,
                            )

            _ax.axvline(0, color="#555555", linestyle="--", linewidth=1.0, alpha=0.6)
            _ax.axhline(0, color="#555555", linestyle="--", linewidth=1.0, alpha=0.6)

            if _row_xlim[_metric] is not None:
                _ax.set_xlim(_row_xlim[_metric])
            _ax.set_ylim(_y_lim)

            if _row_idx == 0:
                _ax.set_title(_model_display_tradeoff_plot.get(_model, _model), fontsize=11)
            if _col_idx == 0:
                _ax.set_ylabel("Simple geometry score")
            _ax.set_xlabel(_xlabel)

    _legend_handles = [
        _Line2D(
            [0],
            [0],
            marker=_guidance_markers[_guidance],
            color="none",
            markerfacecolor=_guidance_palette[_guidance],
            markeredgecolor="black",
            linestyle="none",
            markersize=7,
            label=f"{_guidance.capitalize()} guidance",
        )
        for _guidance in _guidance_order
    ]
    _legend_handles.extend(
        [
            _Line2D(
                [0],
                [0],
                marker=_highlight_styles[_flag]["marker"],
                color="none",
                markerfacecolor="none",
                markeredgecolor=_highlight_styles[_flag]["color"],
                markeredgewidth=2.0,
                markersize=8,
                label=_highlight_styles[_flag]["label"],
            )
            for _flag in [
                "is_balance_best",
                "is_balance_worst",
                "is_peer_best",
                "is_peer_worst",
            ]
        ]
    )

    _fig.legend(
        handles=_legend_handles,
        loc="upper center",
        ncol=4,
        frameon=False,
        bbox_to_anchor=(0.5, 1.03),
    )
    _fig.suptitle(
        "Tradeoff highlights with simple geometry penalty",
        y=1.08,
        fontsize=13,
    )
    _fig.tight_layout()
    save_plot_png(_fig, "14_tradeoff_highlight_cases.png")
    plt.show()
    _tradeoff_highlight_plot = _fig
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        "\n".join(
            [
                "### Same highlighted cases with combined Tortoize score",
                "",
                "This companion figure uses the same highlighted cases, but the y-axis",
                "is the **combined Tortoize score**:",
                "`0.5 × (|Ramachandran Z| + |Torsion Z|)`.",
            ]
        )
    )
    return


@app.cell
def _(
    add_tradeoff_guidance_regions,
    plt,
    save_plot_png,
    tradeoff_df,
    tradeoff_highlight_points,
):
    from matplotlib.lines import Line2D as _Line2D

    _model_display_tradeoff_plot = {
        "boltz2X-RAY": "Boltz2-XRD",
        "boltz2MD": "Boltz2-MD",
        "protenix": "Protenix",
        "rf3": "RF3",
    }
    _guidance_markers = {
        "weak": "o",
        "medium": "s",
        "strong": "^",
    }
    _guidance_palette = {
        "weak": "#377EB8",
        "medium": "#FF7F00",
        "strong": "#E41A1C",
    }
    _highlight_styles = {
        "is_balance_best": {
            "color": "#1B9E77",
            "marker": "o",
            "label": "Best balance",
            "size": 180,
        },
        "is_balance_worst": {
            "color": "#CC3366",
            "marker": "o",
            "label": "Worst balance",
            "size": 180,
        },
        "is_peer_best": {
            "color": "#1F78B4",
            "marker": "D",
            "label": "Differentially better than peers",
            "size": 135,
        },
        "is_peer_worst": {
            "color": "#6A3D9A",
            "marker": "D",
            "label": "Differentially worse than peers",
            "size": 135,
        },
    }
    _label_templates = {
        ("right", "up"): [(8, 10), (8, 22), (8, -2), (8, 34)],
        ("right", "down"): [(8, -8), (8, -20), (8, 4), (8, -32)],
        ("left", "up"): [(-10, 10), (-10, 22), (-10, -2), (-10, 34)],
        ("left", "down"): [(-10, -8), (-10, -20), (-10, 4), (-10, -32)],
    }

    _metric_rows = [
        ("delta_rscc", "Δ RSCC vs baseline"),
        ("delta_silhouette_to_ref", "Δ LDDT silhouette-to-ref vs baseline"),
    ]
    _model_order = ["boltz2X-RAY", "boltz2MD", "protenix", "rf3"]
    _guidance_order = ["weak", "medium", "strong"]
    _y_col = "tortoize_combined_score"
    _y_values = tradeoff_df[_y_col].replace([float("inf"), float("-inf")], float("nan")).dropna()
    _y_pad = max(0.16 * (_y_values.max() - _y_values.min()), 0.30)
    _y_lim = (float(_y_values.min() - _y_pad), float(_y_values.max() + _y_pad))

    _row_xlim = {}
    for _metric, _ in _metric_rows:
        _values = tradeoff_df[_metric].replace([float("inf"), float("-inf")], float("nan")).dropna()
        if _values.empty:
            _row_xlim[_metric] = None
            continue
        _max_abs = float(_values.abs().max())
        if _max_abs == 0:
            _max_abs = 0.01
        _pad = max(0.12 * _max_abs, 0.02)
        _row_xlim[_metric] = (-_max_abs - _pad, _max_abs + _pad)

    _fig, _axes = plt.subplots(
        2, len(_model_order), figsize=(18, 8.2), sharex="row", sharey=True
    )

    for _row_idx, (_metric, _xlabel) in enumerate(_metric_rows):
        for _col_idx, _model in enumerate(_model_order):
            _ax = _axes[_row_idx, _col_idx]
            _model_df = tradeoff_df[tradeoff_df["model_name"] == _model]

            if _model_df.empty:
                _ax.axis("off")
                continue

            add_tradeoff_guidance_regions(
                _ax,
                _model_df,
                _metric,
                _y_col,
                _guidance_order,
                _guidance_palette,
            )

            for _guidance_strength in _guidance_order:
                _guidance_df = _model_df[
                    _model_df["guidance_strength"] == _guidance_strength
                ]
                if _guidance_df.empty:
                    continue
                _ax.scatter(
                    _guidance_df[_metric],
                    _guidance_df[_y_col],
                    s=14,
                    c=_guidance_palette[_guidance_strength],
                    alpha=0.22,
                    marker=_guidance_markers[_guidance_strength],
                    edgecolors="white",
                    linewidths=0.2,
                    zorder=1.2,
                )

            _highlight_df = tradeoff_highlight_points[
                tradeoff_highlight_points["model_name"] == _model
            ].copy()
            if not _highlight_df.empty:
                _highlight_df = _highlight_df.sort_values(
                    [_y_col, _metric], ascending=[False, True]
                )

                for _flag, _style in _highlight_styles.items():
                    _subset = _highlight_df[_highlight_df[_flag]]
                    if _subset.empty:
                        continue
                    _ax.scatter(
                        _subset[_metric],
                        _subset[_y_col],
                        s=_style["size"],
                        facecolors="none",
                        edgecolors=_style["color"],
                        linewidths=2.0,
                        marker=_style["marker"],
                        zorder=4,
                    )

                for _guidance_strength in _guidance_order:
                    _highlight_guidance_df = _highlight_df[
                        _highlight_df["guidance_strength"] == _guidance_strength
                    ]
                    if _highlight_guidance_df.empty:
                        continue
                    _ax.scatter(
                        _highlight_guidance_df[_metric],
                        _highlight_guidance_df[_y_col],
                        s=42,
                        c="black",
                        marker=_guidance_markers[_guidance_strength],
                        zorder=5,
                    )

                if _row_idx == 0 and _row_xlim[_metric] is not None:
                    _x_min, _x_max = _row_xlim[_metric]
                    _y_min, _y_max = _y_lim
                    _highlight_df["x_frac"] = (
                        (_highlight_df[_metric] - _x_min) / (_x_max - _x_min)
                    )
                    _highlight_df["y_frac"] = (
                        (_highlight_df[_y_col] - _y_min) / (_y_max - _y_min)
                    )
                    _highlight_df["label_side"] = _highlight_df["x_frac"].map(
                        lambda value: (
                            "right"
                            if value < 0.22
                            else "left"
                            if value > 0.78
                            else "right"
                            if value < 0.5
                            else "left"
                        )
                    )
                    _highlight_df["label_vertical"] = _highlight_df["y_frac"].map(
                        lambda value: "down" if value > 0.72 else "up"
                    )
                    for (_label_side, _label_vertical), _side_df in _highlight_df.groupby(
                        ["label_side", "label_vertical"], sort=False
                    ):
                        _side_df = _side_df.sort_values(
                            [_y_col, _metric], ascending=[False, True]
                        ).reset_index(drop=True)
                        _side_offsets = _label_templates[(_label_side, _label_vertical)]
                        for _label_idx, _highlight_row in enumerate(
                            _side_df.itertuples(index=False)
                        ):
                            _offset = _side_offsets[min(_label_idx, len(_side_offsets) - 1)]
                            _ax.annotate(
                                _highlight_row.point_label,
                                xy=(
                                    getattr(_highlight_row, _metric),
                                    getattr(_highlight_row, _y_col),
                                ),
                                xytext=_offset,
                                textcoords="offset points",
                                fontsize=6.4,
                                ha="left" if _label_side == "right" else "right",
                                va="center",
                                bbox={
                                    "boxstyle": "round,pad=0.14",
                                    "fc": "white",
                                    "ec": "#666666",
                                    "alpha": 0.92,
                                    "lw": 0.65,
                                },
                                arrowprops={
                                    "arrowstyle": "-",
                                    "color": "#666666",
                                    "lw": 0.55,
                                    "shrinkA": 0,
                                    "shrinkB": 4,
                                },
                                zorder=6,
                            )

            _ax.axvline(0, color="#555555", linestyle="--", linewidth=1.0, alpha=0.6)
            _ax.axhline(0, color="#555555", linestyle="--", linewidth=1.0, alpha=0.6)

            if _row_xlim[_metric] is not None:
                _ax.set_xlim(_row_xlim[_metric])
            _ax.set_ylim(_y_lim)

            if _row_idx == 0:
                _ax.set_title(_model_display_tradeoff_plot.get(_model, _model), fontsize=11)
            if _col_idx == 0:
                _ax.set_ylabel("Combined Tortoize")
            _ax.set_xlabel(_xlabel)

    _legend_handles = [
        _Line2D(
            [0],
            [0],
            marker=_guidance_markers[_guidance],
            color="none",
            markerfacecolor=_guidance_palette[_guidance],
            markeredgecolor="black",
            linestyle="none",
            markersize=7,
            label=f"{_guidance.capitalize()} guidance",
        )
        for _guidance in _guidance_order
    ]
    _legend_handles.extend(
        [
            _Line2D(
                [0],
                [0],
                marker=_highlight_styles[_flag]["marker"],
                color="none",
                markerfacecolor="none",
                markeredgecolor=_highlight_styles[_flag]["color"],
                markeredgewidth=2.0,
                markersize=8,
                label=_highlight_styles[_flag]["label"],
            )
            for _flag in [
                "is_balance_best",
                "is_balance_worst",
                "is_peer_best",
                "is_peer_worst",
            ]
        ]
    )

    _fig.legend(
        handles=_legend_handles,
        loc="upper center",
        ncol=4,
        frameon=False,
        bbox_to_anchor=(0.5, 1.03),
    )
    _fig.suptitle(
        "Highlighted tradeoff cases with combined Tortoize score",
        y=1.08,
        fontsize=13,
    )
    _fig.tight_layout()
    save_plot_png(_fig, "15_tradeoff_highlight_cases_tortoize.png")
    save_plot_png(_fig, "15_tradeoff_highlight_cases_absolute_geometry.png")
    plt.show()
    _tradeoff_highlight_plot_tortoize = _fig
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        "\n".join(
            [
                "### Highlighted tradeoff cases by model",
                "",
                "The new companion plot focuses on four case types for each model:",
                "",
                "- **Best balance**: highest balance score for that model.",
                "- **Worst balance**: lowest balance score for that model.",
                "- **Differentially better than peers**: largest `balance_score - mean(other models)`",
                "  for the same protein and guidance-strength bucket.",
                "- **Differentially worse than peers**: smallest value of that same peer-relative score.",
            ]
        )
    )
    return


@app.cell
def _(tradeoff_case_table):
    _model_display = {
        "boltz2X-RAY": "Boltz2-XRD",
        "boltz2MD": "Boltz2-MD",
        "protenix": "Protenix",
        "rf3": "RF3",
    }
    tradeoff_case_table.assign(
        model_name=lambda df: df["model_name"].map(lambda x: _model_display.get(x, x))
    ).rename(
        columns={
            "model_name": "Model",
            "Best balance": "Best",
            "Worst balance": "Worst",
            "Differentially better than peers": "Differentially better",
            "Differentially worse than peers": "Differentially worse",
        }
    )
    return


@app.cell
def _(tradeoff_case_summary):
    _model_display = {
        "boltz2X-RAY": "Boltz2-XRD",
        "boltz2MD": "Boltz2-MD",
        "protenix": "Protenix",
        "rf3": "RF3",
    }
    tradeoff_case_summary.assign(
        model_name=lambda df: df["model_name"].map(lambda x: _model_display.get(x, x)),
        guidance_strength=lambda df: df["guidance_strength"].str.capitalize(),
    ).rename(
        columns={
            "protein": "Protein",
            "model_name": "Model",
            "guidance_strength": "Guidance",
            "guidance_weight": "Guidance weight",
            "case_type": "Case",
            "delta_rscc": "ΔRSCC",
            "delta_silhouette_to_ref": "ΔLDDT-to-ref",
            "geometry_penalty": "Simple geometry penalty",
            "tortoize_combined_score": "Combined Tortoize score",
            "balance_score": "Balance score",
            "balance_vs_peers": "Balance vs peers",
        }
    ).sort_values(["Model", "Case"])
    return


@app.cell(hide_code=True)
def _(geometry_weights, mo, weak_tradeoff):
    _model_display_score_docs = {
        "boltz2X-RAY": "Boltz2-XRD",
        "boltz2MD": "Boltz2-MD",
        "protenix": "Protenix",
        "rf3": "RF3",
    }

    _summary = (
        weak_tradeoff.groupby("model_name", as_index=False)
        .agg(
            geometry_penalty_median=("geometry_penalty", "median"),
            geometry_penalty_mean=("geometry_penalty", "mean"),
            tortoize_combined_median=("tortoize_combined_score", "median"),
            tortoize_combined_mean=("tortoize_combined_score", "mean"),
            balance_score_median=("balance_score", "median"),
            balance_score_mean=("balance_score", "mean"),
            rscc_gain_fraction=("delta_rscc", lambda s: (s > 0).mean()),
            lddt_gain_fraction=("delta_silhouette_to_ref", lambda s: (s > 0).mean()),
        )
        .sort_values("model_name")
    )

    _summary_rows = []
    for _row in _summary.itertuples(index=False):
        _summary_rows.append(
            "| "
            + " | ".join(
                [
                    _model_display_score_docs.get(_row.model_name, _row.model_name),
                    f"{_row.geometry_penalty_median:.2f}",
                    f"{_row.geometry_penalty_mean:.2f}",
                    f"{_row.tortoize_combined_median:.2f}",
                    f"{_row.tortoize_combined_mean:.2f}",
                    f"{_row.balance_score_median:.2f}",
                    f"{_row.balance_score_mean:.2f}",
                    f"{100 * _row.rscc_gain_fraction:.0f}%",
                    f"{100 * _row.lddt_gain_fraction:.0f}%",
                ]
            )
            + " |"
        )

    _component_rows = []
    for _row in geometry_weights.itertuples(index=False):
        _component_rows.append(
            "| "
            + " | ".join(
                [
                    _row.display_name,
                    f"{100 * _row.weight:.1f}%",
                ]
            )
            + " |"
        )

    mo.md(
        "\n".join(
            [
                "### Balance and geometry score definitions",
                "",
                "- **Simple geometry penalty** is the unweighted mean of the baseline-relative terms we can compute directly:",
                "  - `Δ clashscore`",
                "  - `100 × Δ bond-length outlier fraction` (percentage-point change)",
                "  - `100 × Δ bond-angle outlier fraction` (percentage-point change)",
                "- So `geometry_penalty = mean(Δ clashscore, Δ bond-length pp, Δ bond-angle pp)`.",
                "- **Combined Tortoize score** is shown separately because we do not currently have baseline Tortoize values:",
                "  - `tortoize_combined_score = 0.5 × (|Ramachandran Z| + |Torsion Z|)`",
                "- **Balance score** is used only to rank highlighted cases:",
                "  - `robust_z_model(ΔRSCC) + robust_z_model(ΔLDDT_to_ref) - robust_z_model(simple geometry penalty)`",
                "",
                "#### Simple geometry-penalty components",
                "",
                "| Component | Weight in simple score |",
                "|---|---:|",
                *_component_rows,
                "",
                "#### Weak-guidance summary by model",
                "",
                "| Model | Simple geometry penalty (median) | Simple geometry penalty (mean) | Combined Tortoize score (median) | Combined Tortoize score (mean) | Balance score (median) | Balance score (mean) | % with ΔRSCC > 0 | % with ΔLDDT-to-ref > 0 |",
                "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
                *_summary_rows,
            ]
        )
    )
    return


@app.cell(hide_code=True)
def _(
    focus_5i09,
    lddt_tradeoff_examples,
    mo,
    rscc_tradeoff_examples,
    well_balanced_examples,
):
    _model_display_callouts = {
        "boltz2X-RAY": "Boltz2-XRD",
        "boltz2MD": "Boltz2-MD",
        "protenix": "Protenix",
        "rf3": "RF3",
    }

    lines = ["### Example callouts", ""]

    if not focus_5i09.empty:
        lines.append("**5I09 at weak guidance:**")
        for _row in focus_5i09.itertuples(index=False):
            lines.append(
                f"- **{_model_display_callouts.get(_row.model_name, _row.model_name)}**: "
                f"ΔRSCC={_row.delta_rscc:.3f}, "
                f"ΔLDDT-to-ref={_row.delta_silhouette_to_ref:.3f}, "
                f"|Rama Z|={_row.abs_ramachandran_z_score:.3f}, "
                f"|Torsion Z|={_row.abs_torsion_z_score:.3f}, "
                f"simple geometry penalty={_row.geometry_penalty:.2f}."
            )
        lines.append("")

    lines.append("**Examples that are relatively well-balanced across metrics:**")
    best_per_model = (
        well_balanced_examples.sort_values(["model_name", "balance_score"], ascending=[True, False])
        .groupby("model_name")
        .head(1)
    )
    for _row in best_per_model.itertuples(index=False):
        lines.append(
            f"- **{_model_display_callouts.get(_row.model_name, _row.model_name)} / {_row.protein.upper()}**: "
            f"ΔRSCC={_row.delta_rscc:.3f}, "
            f"ΔLDDT-to-ref={_row.delta_silhouette_to_ref:.3f}, "
            f"simple geometry penalty={_row.geometry_penalty:.2f}."
        )

    lines.append("")
    lines.append("**Examples with single-metric gains but trade-offs elsewhere:**")

    rscc_focus = (
        rscc_tradeoff_examples.sort_values(["model_name", "delta_rscc"], ascending=[True, False])
        .groupby("model_name")
        .head(1)
    )
    for _row in rscc_focus.itertuples(index=False):
        lines.append(
            f"- RSCC-heavy: **{_model_display_callouts.get(_row.model_name, _row.model_name)} / {_row.protein.upper()}** "
            f"(ΔRSCC={_row.delta_rscc:.3f}, ΔLDDT-to-ref={_row.delta_silhouette_to_ref:.3f}, "
            f"simple geometry penalty={_row.geometry_penalty:.2f})."
        )

    lddt_focus = (
        lddt_tradeoff_examples.sort_values(["model_name", "delta_silhouette_to_ref"], ascending=[True, False])
        .groupby("model_name")
        .head(1)
    )
    for _row in lddt_focus.itertuples(index=False):
        lines.append(
            f"- LDDT-heavy: **{_model_display_callouts.get(_row.model_name, _row.model_name)} / {_row.protein.upper()}** "
            f"(ΔLDDT-to-ref={_row.delta_silhouette_to_ref:.3f}, ΔRSCC={_row.delta_rscc:.3f}, "
            f"simple geometry penalty={_row.geometry_penalty:.2f})."
        )

    mo.md("\n".join(lines))
    return


@app.cell
def _(well_balanced_examples):
    well_balanced_examples
    return


@app.cell
def _(rscc_tradeoff_examples):
    rscc_tradeoff_examples
    return


@app.cell
def _(lddt_tradeoff_examples):
    lddt_tradeoff_examples
    return


if __name__ == "__main__":
    app.run()
