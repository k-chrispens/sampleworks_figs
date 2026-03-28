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
        load_and_tag_path_csv,
        tag_path_dataframe,
        load_rscc,
        load_lddt,
        make_delta_frame,
        make_tradeoff_grid,
        annotate_tradeoff_highlights,
        compute_conceptual_tradeoff_limits,
        make_conceptual_tradeoff_figure,
        make_spiral_figure,
        save_fig,
        robust_z_score,
        MODEL_DISPLAY,
        MODEL_NAME_ALIAS,
    )

    return (
        MODEL_DISPLAY,
        MODEL_NAME_ALIAS,
        add_simplified_metadata_columns,
        annotate_tradeoff_highlights,
        compute_conceptual_tradeoff_limits,
        load_and_tag_path_csv,
        load_lddt,
        load_rscc,
        make_conceptual_tradeoff_figure,
        make_delta_frame,
        make_spiral_figure,
        make_tradeoff_grid,
        pd,
        plt,
        robust_z_score,
        save_fig,
        tag_path_dataframe,
    )


@app.cell
def _(mo):
    mo.md(r"""
    # Geometry vs fit trade-offs (across guidance strengths)

    We combine the baseline-relative geometry metrics we can measure directly
    into a simple **geometry penalty** score, then compare that penalty against
    ΔRSCC and ΔLDDT-to-reference.

    The simple geometry score averages Δ clashscore and the percentage-point
    changes in bond-length and bond-angle outlier fractions. Tortoize terms are
    shown separately in a companion plot.
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
def _(load_and_tag_path_csv):
    dihedrals = load_and_tag_path_csv(
        "initial40results/initial_40_tortoize_protein_stats.csv",
        filename_col="path",
        extra_renames={"model": "model_n"},
    )
    return (dihedrals,)


@app.cell
def _(add_simplified_metadata_columns, load_rscc):
    _rscc = load_rscc("initial40results/rscc_results.csv")
    rscc_all = _rscc.groupby(["protein", "model_name", "guidance_weight"]).last().reset_index()
    rscc_all = add_simplified_metadata_columns(rscc_all)
    return (rscc_all,)


@app.cell
def _(add_simplified_metadata_columns, load_lddt):
    lddt_results_1 = add_simplified_metadata_columns(
        load_lddt("initial40results/lddt_results_w_baseline.csv")
    )
    return (lddt_results_1,)


@app.cell
def _(lddt_results_1):
    _baseline = lddt_results_1.query('guidance_strength == "baseline"')[
        ["protein", "model_name", "selection_name", "avg_silhouette_to_ref"]
    ]
    lddt_merged = _baseline.merge(
        lddt_results_1.query('guidance_strength != "baseline"'),
        on=["protein", "model_name", "selection_name"],
        suffixes=["_baseline", ""],
    )
    lddt_merged["delta_silhouette_to_ref"] = (
        lddt_merged.avg_silhouette_to_ref - lddt_merged.avg_silhouette_to_ref_baseline
    )
    return (lddt_merged,)


@app.cell
def _(
    bond_stats,
    clashes_1,
    dihedrals,
    lddt_merged,
    make_delta_frame,
    rscc_all,
):
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
    dihedral_stats = (
        dihedrals.groupby(_group_cols, as_index=False)
        .agg(ramachandran_z_score=("ramachandran_z_score", "mean"),
             torsion_z_score=("torsion_z_score", "mean"))
        .dropna(subset=["ramachandran_z_score", "torsion_z_score"])
    )
    dihedral_stats["abs_ramachandran_z_score"] = dihedral_stats["ramachandran_z_score"].abs()
    dihedral_stats["abs_torsion_z_score"] = dihedral_stats["torsion_z_score"].abs()

    clash_delta = make_delta_frame(clash_stats, "clashscore")[
        _group_cols + ["delta_clashscore"]
    ]
    bond_angle_delta = make_delta_frame(bond_stats_agg, "outlier_fraction_bond_angles")[
        _group_cols + ["delta_outlier_fraction_bond_angles"]
    ]
    bond_length_delta = make_delta_frame(bond_stats_agg, "outlier_fraction_bond_lengths")[
        _group_cols + ["delta_outlier_fraction_bond_lengths"]
    ]
    rama_z_term = dihedral_stats[_group_cols + ["abs_ramachandran_z_score"]]
    torsion_z_term = dihedral_stats[_group_cols + ["abs_torsion_z_score"]]
    rscc_delta_tradeoff = make_delta_frame(rscc_stats, "rscc")[
        _group_cols + ["delta_rscc"]
    ]
    lddt_delta = (
        lddt_merged.query("selection_name != 'all altloc'")
        .groupby(_group_cols, as_index=False)
        .agg(delta_silhouette_to_ref=("delta_silhouette_to_ref", "mean"))
    )
    return (
        bond_angle_delta,
        bond_length_delta,
        clash_delta,
        lddt_delta,
        rama_z_term,
        rscc_delta_tradeoff,
        torsion_z_term,
    )


@app.cell
def _(
    MODEL_NAME_ALIAS,
    bond_angle_delta,
    bond_length_delta,
    clash_delta,
    lddt_delta,
    pd,
    rama_z_term,
    robust_z_score,
    rscc_delta_tradeoff,
    torsion_z_term,
):
    _merge_on = ["protein", "model_name", "guidance_weight", "guidance_strength"]
    _alias = lambda df: df.assign(
        model_name=df["model_name"].map(lambda x: MODEL_NAME_ALIAS.get(x, x))
    )

    tradeoff_df = (
        _alias(clash_delta)
        .merge(_alias(bond_angle_delta), on=_merge_on)
        .merge(_alias(bond_length_delta), on=_merge_on)
        .merge(_alias(rama_z_term), on=_merge_on)
        .merge(_alias(torsion_z_term), on=_merge_on)
        .merge(_alias(rscc_delta_tradeoff), on=_merge_on)
        .merge(_alias(lddt_delta), on=_merge_on)
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

    _geom_terms = [
        "delta_clashscore",
        "delta_outlier_fraction_bond_lengths_pct",
        "delta_outlier_fraction_bond_angles_pct",
    ]
    geometry_weights = pd.DataFrame({
        "geometry_term": _geom_terms,
        "display_name": [
            "Δ clashscore",
            "Δ bond-length outlier fraction (percentage points)",
            "Δ bond-angle outlier fraction (percentage points)",
        ],
        "weight": [1 / 3, 1 / 3, 1 / 3],
    })

    tradeoff_df["geometry_penalty"] = tradeoff_df[_geom_terms].mean(axis=1)
    tradeoff_df["tortoize_combined_score"] = 0.5 * (
        tradeoff_df["abs_ramachandran_z_score"] + tradeoff_df["abs_torsion_z_score"]
    )

    # Balance scores — must be computed before weak_tradeoff is extracted
    tradeoff_df["geometry_penalty_z_model"] = (
        tradeoff_df.groupby("model_name")["geometry_penalty"]
        .transform(robust_z_score).clip(-5, 5)
    )
    for _col in ["delta_rscc", "delta_silhouette_to_ref"]:
        tradeoff_df[f"{_col}_z_model"] = (
            tradeoff_df.groupby("model_name")[_col]
            .transform(robust_z_score).clip(-5, 5)
        )
    tradeoff_df["balance_score"] = (
        tradeoff_df["delta_rscc_z_model"]
        + tradeoff_df["delta_silhouette_to_ref_z_model"]
        - tradeoff_df["geometry_penalty_z_model"]
    )
    tradeoff_df["peer_model_count"] = tradeoff_df.groupby(
        ["protein", "guidance_strength"]
    )["model_name"].transform("nunique")
    _peer_sum = tradeoff_df.groupby(
        ["protein", "guidance_strength"]
    )["balance_score"].transform("sum")
    tradeoff_df["peer_balance_mean_other"] = (
        (_peer_sum - tradeoff_df["balance_score"])
        / (tradeoff_df["peer_model_count"] - 1)
    ).where(tradeoff_df["peer_model_count"] > 1)
    tradeoff_df["balance_vs_peers"] = (
        tradeoff_df["balance_score"] - tradeoff_df["peer_balance_mean_other"]
    )

    weak_tradeoff = tradeoff_df.query("guidance_strength == 'weak'").copy()
    return geometry_weights, tradeoff_df, weak_tradeoff


@app.cell
def _(pd, tradeoff_df):
    _best = (
        tradeoff_df.sort_values(["model_name", "balance_score"], ascending=[True, False])
        .groupby("model_name", group_keys=False).head(1).copy()
    )
    _best["highlight_kind"] = "best"
    _best["case_type"] = "Best balance"

    _worst = (
        tradeoff_df.sort_values(["model_name", "balance_score"], ascending=[True, True])
        .groupby("model_name", group_keys=False).head(1).copy()
    )
    _worst["highlight_kind"] = "worst"
    _worst["case_type"] = "Worst balance"

    _peer_best = (
        tradeoff_df.query("peer_model_count > 1")
        .sort_values(["model_name", "balance_vs_peers"], ascending=[True, False])
        .groupby("model_name", group_keys=False).head(1).copy()
    )
    _peer_best["highlight_kind"] = "peer+"
    _peer_best["case_type"] = "Differentially better than peers"

    _peer_worst = (
        tradeoff_df.query("peer_model_count > 1")
        .sort_values(["model_name", "balance_vs_peers"], ascending=[True, True])
        .groupby("model_name", group_keys=False).head(1).copy()
    )
    _peer_worst["highlight_kind"] = "peer-"
    _peer_worst["case_type"] = "Differentially worse than peers"

    _summary_cols = [
        "protein", "model_name", "guidance_strength", "guidance_weight",
        "highlight_kind", "case_type", "delta_rscc", "delta_silhouette_to_ref",
        "geometry_penalty", "tortoize_combined_score", "balance_score", "balance_vs_peers",
    ]
    tradeoff_case_summary = pd.concat(
        [_best, _worst, _peer_best, _peer_worst], ignore_index=True
    )[_summary_cols].copy()
    tradeoff_case_summary["case_note"] = (
        tradeoff_case_summary["protein"].str.upper()
        + " (" + tradeoff_case_summary["guidance_strength"].str[0].str.upper() + ")"
    )

    tradeoff_case_table = (
        tradeoff_case_summary.pivot(
            index="model_name", columns="case_type", values="case_note"
        )
        .reindex(["boltz2X-RAY", "boltz2MD", "protenix", "rf3"])
        .reset_index()
    )

    tradeoff_highlight_points = (
        tradeoff_case_summary.groupby(
            [
                "protein", "model_name", "guidance_strength", "guidance_weight",
                "delta_rscc", "delta_silhouette_to_ref", "geometry_penalty",
                "tortoize_combined_score", "balance_score", "balance_vs_peers",
            ],
            as_index=False, dropna=False,
        ).agg(
            highlight_kinds=("highlight_kind", lambda s: tuple(dict.fromkeys(s.tolist()))),
            case_types=("case_type", lambda s: tuple(dict.fromkeys(s.tolist()))),
        )
    )
    for _flag, _kind in [
        ("is_balance_best", "best"), ("is_balance_worst", "worst"),
        ("is_peer_best", "peer+"), ("is_peer_worst", "peer-"),
    ]:
        tradeoff_highlight_points[_flag] = tradeoff_highlight_points["highlight_kinds"].map(
            lambda vals, k=_kind: k in vals
        )
    tradeoff_highlight_points["point_label"] = tradeoff_highlight_points.apply(
        lambda r: f"{r['protein'].upper()} ({r['guidance_strength'][0].upper()})", axis=1
    )
    return (
        tradeoff_case_summary,
        tradeoff_case_table,
        tradeoff_highlight_points,
    )


@app.cell
def _(weak_tradeoff):
    focus_5i09 = weak_tradeoff.query("protein == '5i09'")[
        [
            "protein", "model_name", "delta_rscc", "delta_silhouette_to_ref",
            "abs_ramachandran_z_score", "abs_torsion_z_score",
            "delta_clashscore", "delta_outlier_fraction_bond_angles",
            "delta_outlier_fraction_bond_lengths",
            "geometry_penalty", "balance_score",
        ]
    ].sort_values("model_name")

    well_balanced_examples = (
        weak_tradeoff.query(
            "delta_rscc > 0 and delta_silhouette_to_ref > 0 and geometry_penalty <= 0"
        )
        .sort_values(["model_name", "balance_score"], ascending=[True, False])
        .groupby("model_name").head(3)[
            ["protein", "model_name", "delta_rscc", "delta_silhouette_to_ref",
             "geometry_penalty", "balance_score"]
        ]
    )

    rscc_tradeoff_examples = (
        weak_tradeoff.query(
            "delta_rscc > 0 and (delta_silhouette_to_ref < 0 or geometry_penalty > 0)"
        )
        .sort_values(["model_name", "delta_rscc"], ascending=[True, False])
        .groupby("model_name").head(3)[
            ["protein", "model_name", "delta_rscc", "delta_silhouette_to_ref",
             "geometry_penalty"]
        ]
    )

    lddt_tradeoff_examples = (
        weak_tradeoff.query(
            "delta_silhouette_to_ref > 0 and (delta_rscc < 0 or geometry_penalty > 0)"
        )
        .sort_values(["model_name", "delta_silhouette_to_ref"], ascending=[True, False])
        .groupby("model_name").head(3)[
            ["protein", "model_name", "delta_rscc", "delta_silhouette_to_ref",
             "geometry_penalty"]
        ]
    )
    return (
        focus_5i09,
        lddt_tradeoff_examples,
        rscc_tradeoff_examples,
        well_balanced_examples,
    )


@app.cell
def _(make_tradeoff_grid, plt, save_fig, tradeoff_df):
    _fig, _axes = make_tradeoff_grid(
        tradeoff_df,
        "geometry_penalty",
        "Simple geometry score",
        trajectory_protein="5i09",
        title="Guidance trade-offs by model (simple geometry penalty; weak/medium/strong)",
    )
    save_fig(_fig, "14_tradeoff_scatter_all_guidance.png")
    plt.show()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("\n".join([
        "### Guidance-weight spiral trajectories",
        "",
        "Each line traces a single protein through (Δ RSCC, geometry penalty)",
        "space as guidance weight increases from baseline (origin) to the",
        "strongest setting. Dots are colored from cool to warm with increasing",
        "guidance weight. Common trajectory shapes include:",
        "",
        "- **Rightward sweep** — RSCC improves with tolerable geometry cost ('correctable').",
        "- **Upward hook / spiral** — geometry degrades faster than RSCC improves at strong guidance.",
        "- **Stuck at origin** — neither metric moves much regardless of guidance.",
    ]))
    return


@app.cell
def _(MODEL_DISPLAY, make_spiral_figure, plt, save_fig, tradeoff_df):
    _spiral_nums = {"boltz2X-RAY": 21, "boltz2MD": 22, "protenix": 23, "rf3": 24}
    _model_order = [
        m for m in ["boltz2X-RAY", "boltz2MD", "protenix", "rf3"]
        if m in tradeoff_df["model_name"].unique()
    ]
    for _model in _model_order:
        _model_df = tradeoff_df[tradeoff_df["model_name"] == _model].copy()
        if _model_df.empty:
            continue
        _slug = MODEL_DISPLAY.get(_model, _model).lower().replace("-", "_")
        _fig = make_spiral_figure(_model_df, _model)
        save_fig(_fig, f"{_spiral_nums[_model]:02d}_spiral_{_slug}.png")
        plt.show()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("\n".join([
        "### Conceptual overlaid RSCC tradeoff views",
        "",
        "These companion figures collapse the conceptual view to a single",
        "RSCC-versus-geometry axis, with color indicating model family and",
        "no model names shown in the plot itself.",
    ]))
    return


@app.cell
def _(compute_conceptual_tradeoff_limits, tradeoff_df):
    conceptual_xlim, conceptual_ylim = compute_conceptual_tradeoff_limits(tradeoff_df)
    return conceptual_xlim, conceptual_ylim


@app.cell
def _(
    conceptual_xlim,
    conceptual_ylim,
    make_conceptual_tradeoff_figure,
    plt,
    save_fig,
    tradeoff_df,
):
    _fig = make_conceptual_tradeoff_figure(
        tradeoff_df, "weak",
        xlim=conceptual_xlim, ylim=conceptual_ylim,
        title="Tradeoff Landscape (Weak Guidance)",
    )
    save_fig(_fig, "18_tradeoff_conceptual_weak.png")
    plt.show()
    return


@app.cell
def _(
    conceptual_xlim,
    conceptual_ylim,
    make_conceptual_tradeoff_figure,
    plt,
    save_fig,
    tradeoff_df,
):
    _fig = make_conceptual_tradeoff_figure(
        tradeoff_df, "weak",
        xlim=conceptual_xlim, ylim=conceptual_ylim,
        title="Tradeoff Landscape (Weak Guidance, Zoomed)",
        zoom_quantiles=(0.02, 0.95),
    )
    save_fig(_fig, "19_tradeoff_conceptual_weak_zoomed.png")
    plt.show()
    return


@app.cell
def _(
    conceptual_xlim,
    conceptual_ylim,
    make_conceptual_tradeoff_figure,
    plt,
    save_fig,
    tradeoff_df,
):
    _fig = make_conceptual_tradeoff_figure(
        tradeoff_df, "strong",
        xlim=conceptual_xlim, ylim=conceptual_ylim,
        title="Tradeoff Landscape (Strong Guidance)",
    )
    save_fig(_fig, "20_tradeoff_conceptual_strong.png")
    plt.show()
    return


@app.cell
def _(
    annotate_tradeoff_highlights,
    make_tradeoff_grid,
    plt,
    save_fig,
    tradeoff_df,
    tradeoff_highlight_points,
):
    _fig, _axes = make_tradeoff_grid(
        tradeoff_df, "geometry_penalty", "Simple geometry score",
        figsize=(18, 8.2), scatter_size=14, scatter_alpha=0.22, show_centroids=False,
    )
    annotate_tradeoff_highlights(
        _fig, _axes, tradeoff_highlight_points, "geometry_penalty",
        title="Tradeoff highlights with simple geometry penalty",
    )
    save_fig(_fig, "15_tradeoff_highlight_cases.png")
    plt.show()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("\n".join([
        "### Same highlighted cases with combined Tortoize score",
        "",
        "This companion figure uses the same highlighted cases, but the y-axis",
        "is the **combined Tortoize score**:",
        "`0.5 * (|Ramachandran Z| + |Torsion Z|)`.",
    ]))
    return


@app.cell
def _(
    annotate_tradeoff_highlights,
    make_tradeoff_grid,
    plt,
    save_fig,
    tradeoff_df,
    tradeoff_highlight_points,
):
    _fig, _axes = make_tradeoff_grid(
        tradeoff_df, "tortoize_combined_score", "Combined Tortoize",
        figsize=(18, 8.2), scatter_size=14, scatter_alpha=0.22, show_centroids=False,
    )
    annotate_tradeoff_highlights(
        _fig, _axes, tradeoff_highlight_points, "tortoize_combined_score",
        title="Highlighted tradeoff cases with combined Tortoize score",
    )
    save_fig(_fig, "16_tradeoff_highlight_cases_tortoize.png")
    save_fig(_fig, "17_tradeoff_highlight_cases_absolute_geometry.png")
    plt.show()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("\n".join([
        "### Highlighted tradeoff cases by model",
        "",
        "The new companion plot focuses on four case types for each model:",
        "",
        "- **Best balance**: highest balance score for that model.",
        "- **Worst balance**: lowest balance score for that model.",
        "- **Differentially better than peers**: largest `balance_score - mean(other models)`",
        "  for the same protein and guidance-strength bucket.",
        "- **Differentially worse than peers**: smallest value of that same peer-relative score.",
    ]))
    return


@app.cell
def _(MODEL_DISPLAY, tradeoff_case_table):
    tradeoff_case_table.assign(
        model_name=lambda df: df["model_name"].map(
            lambda x: MODEL_DISPLAY.get(x, x)
        )
    ).rename(columns={
        "model_name": "Model",
        "Best balance": "Best",
        "Worst balance": "Worst",
        "Differentially better than peers": "Differentially better",
        "Differentially worse than peers": "Differentially worse",
    })
    return


@app.cell
def _(MODEL_DISPLAY, tradeoff_case_summary):
    tradeoff_case_summary.assign(
        model_name=lambda df: df["model_name"].map(
            lambda x: MODEL_DISPLAY.get(x, x)
        ),
        guidance_strength=lambda df: df["guidance_strength"].str.capitalize(),
    ).rename(columns={
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
    }).sort_values(["Model", "Case"])
    return


@app.cell(hide_code=True)
def _(MODEL_DISPLAY, geometry_weights, mo, weak_tradeoff):
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
    for _r in _summary.itertuples(index=False):
        _summary_rows.append(
            "| " + " | ".join([
                MODEL_DISPLAY.get(_r.model_name, _r.model_name),
                f"{_r.geometry_penalty_median:.2f}",
                f"{_r.geometry_penalty_mean:.2f}",
                f"{_r.tortoize_combined_median:.2f}",
                f"{_r.tortoize_combined_mean:.2f}",
                f"{_r.balance_score_median:.2f}",
                f"{_r.balance_score_mean:.2f}",
                f"{100 * _r.rscc_gain_fraction:.0f}%",
                f"{100 * _r.lddt_gain_fraction:.0f}%",
            ]) + " |"
        )

    _component_rows = []
    for _r in geometry_weights.itertuples(index=False):
        _component_rows.append(
            "| " + " | ".join([_r.display_name, f"{100 * _r.weight:.1f}%"]) + " |"
        )

    mo.md("\n".join([
        "### Balance and geometry score definitions",
        "",
        "- **Simple geometry penalty** is the unweighted mean of the baseline-relative terms we can compute directly:",
        "  - `Δ clashscore`",
        "  - `100 * Δ bond-length outlier fraction` (percentage-point change)",
        "  - `100 * Δ bond-angle outlier fraction` (percentage-point change)",
        "- So `geometry_penalty = mean(Δ clashscore, Δ bond-length pp, Δ bond-angle pp)`.",
        "- **Combined Tortoize score** is shown separately because we do not currently have baseline Tortoize values:",
        "  - `tortoize_combined_score = 0.5 * (|Ramachandran Z| + |Torsion Z|)`",
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
    ]))
    return


@app.cell(hide_code=True)
def _(
    MODEL_DISPLAY,
    focus_5i09,
    lddt_tradeoff_examples,
    mo,
    rscc_tradeoff_examples,
    well_balanced_examples,
):
    lines = ["### Example callouts", ""]

    if not focus_5i09.empty:
        lines.append("**5I09 at weak guidance:**")
        for _r in focus_5i09.itertuples(index=False):
            lines.append(
                f"- **{MODEL_DISPLAY.get(_r.model_name, _r.model_name)}**: "
                f"ΔRSCC={_r.delta_rscc:.3f}, "
                f"ΔLDDT-to-ref={_r.delta_silhouette_to_ref:.3f}, "
                f"|Rama Z|={_r.abs_ramachandran_z_score:.3f}, "
                f"|Torsion Z|={_r.abs_torsion_z_score:.3f}, "
                f"simple geometry penalty={_r.geometry_penalty:.2f}."
            )
        lines.append("")

    lines.append("**Examples that are relatively well-balanced across metrics:**")
    _best = (
        well_balanced_examples.sort_values(["model_name", "balance_score"], ascending=[True, False])
        .groupby("model_name").head(1)
    )
    for _r in _best.itertuples(index=False):
        lines.append(
            f"- **{MODEL_DISPLAY.get(_r.model_name, _r.model_name)} / {_r.protein.upper()}**: "
            f"ΔRSCC={_r.delta_rscc:.3f}, ΔLDDT-to-ref={_r.delta_silhouette_to_ref:.3f}, "
            f"simple geometry penalty={_r.geometry_penalty:.2f}."
        )

    lines.append("")
    lines.append("**Examples with single-metric gains but trade-offs elsewhere:**")
    _rscc = (
        rscc_tradeoff_examples.sort_values(["model_name", "delta_rscc"], ascending=[True, False])
        .groupby("model_name").head(1)
    )
    for _r in _rscc.itertuples(index=False):
        lines.append(
            f"- RSCC-heavy: **{MODEL_DISPLAY.get(_r.model_name, _r.model_name)} / {_r.protein.upper()}** "
            f"(ΔRSCC={_r.delta_rscc:.3f}, ΔLDDT-to-ref={_r.delta_silhouette_to_ref:.3f}, "
            f"simple geometry penalty={_r.geometry_penalty:.2f})."
        )
    _lddt = (
        lddt_tradeoff_examples.sort_values(["model_name", "delta_silhouette_to_ref"], ascending=[True, False])
        .groupby("model_name").head(1)
    )
    for _r in _lddt.itertuples(index=False):
        lines.append(
            f"- LDDT-heavy: **{MODEL_DISPLAY.get(_r.model_name, _r.model_name)} / {_r.protein.upper()}** "
            f"(ΔLDDT-to-ref={_r.delta_silhouette_to_ref:.3f}, ΔRSCC={_r.delta_rscc:.3f}, "
            f"simple geometry penalty={_r.geometry_penalty:.2f})."
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
