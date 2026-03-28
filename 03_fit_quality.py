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
    )

    return (
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
    # RSCC results
    """)
    return


@app.cell
def _(load_rscc):
    rscc = load_rscc("initial40results/rscc_results.csv")
    return (rscc,)


@app.cell
def _(add_simplified_metadata_columns, rscc):
    rscc_all = rscc.groupby(["protein", "model_name", "guidance_weight"]).last().reset_index()
    rscc_all = add_simplified_metadata_columns(rscc_all)
    rscc_all
    return (rscc_all,)


@app.cell
def _(make_ridge_plot, rscc_all, save_fig):
    _fig = make_ridge_plot(
        rscc_all.query("guidance_strength == 'baseline'"),
        "rscc",
        combined=True,
    )
    save_fig(_fig, "08_rscc_baseline_all_models.png")
    return


@app.cell
def _(save_model_legend):
    save_model_legend()
    return


@app.cell
def _(make_ridge_plot, rscc_all, save_fig):
    _g = make_ridge_plot(
        rscc_all.query("guidance_strength in ('baseline', 'weak')"),
        "rscc",
    )
    save_fig(_g, "09_rscc_baseline_vs_weak_distribution.png")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    what's with all the really low RSCC structures? Let's list some:
    """)
    return


@app.cell
def _(pd, rscc_all):
    _msa = pd.read_csv("initial40results/boltz2_msa_depth.csv")
    _msa.pdb_id = _msa.pdb_id.str.lower()
    _w_msa = rscc_all.merge(_msa, left_on="protein", right_on="pdb_id")

    _pstats = pd.read_csv("initial40results/Sampleworks_initial_structures.csv")
    _pstats.PDB_ID = _pstats.PDB_ID.str.lower()
    _w_all = _w_msa.merge(_pstats, left_on="pdb_id", right_on="PDB_ID")
    _w_all.groupby("model_name")[[
        "rscc", "msa_depth", "Resolution", "Rfree",
        "Cluster_Size", "Num_Structure_100_identity",
    ]].corr().reset_index().query("level_1 == 'rscc'")
    return


@app.cell
def _(make_delta_frame, make_ridge_plot, rscc_all, save_fig):
    rscc_delta = make_delta_frame(rscc_all, "rscc")
    _g = make_ridge_plot(
        rscc_delta.query("delta_rscc < 0.3 and delta_rscc > -0.3"),
        "delta_rscc",
        is_delta=True,
    )
    save_fig(_g, "10_delta_rscc_distribution.png")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Do any zero-guidance models achieve any meaningful RSCC?
    """)
    return


@app.cell
def _(rscc):
    rscc.query("guidance_weight == 0.0").pivot(
        index=["protein", "selection_name"],
        columns="model_name",
        values="rscc",
    ).describe()
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
def _(load_lddt):
    lddt_results = load_lddt("initial40results/lddt_results_w_baseline.csv")
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
def _(lddt_results_1, make_ridge_plot, save_fig):
    _fig = make_ridge_plot(
        lddt_results_1.query("guidance_strength == 'baseline' and selection_name != 'all altloc'"),
        plot_column_name="avg_silhouette_to_ref",
        custom_xlabel="Average Silhouette Score to Reference (LDDT Clustering)",
        smoothing=0.25,
        combined=True,
    )
    save_fig(_fig, "11_lddt_to_ref_baseline_distribution.png")
    return


@app.cell
def _(lddt_results_1, make_ridge_plot, save_fig):
    _g = make_ridge_plot(
        lddt_results_1.query("selection_name != 'all altloc'"),
        "avg_silhouette_to_ref",
        smoothing=0.25,
    )
    save_fig(_g, "12_lddt_to_ref_all_guidance_distribution.png")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Yikes
    """)
    return


@app.cell
def _(lddt_results_1):
    lddt_by_pdb_id = (
        lddt_results_1.query("guidance_strength == 'weak' and selection_name != 'all altloc'")
        .pivot(
            index=["protein", "selection_name"],
            columns="model_name",
            values=["avg_silhouette_to_ref", "avg_silhouette"],
        )
    )
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
    _pdb_id, _sel = lddt_by_pdb_id[("avg_silhouette_to_ref", "boltz2X-RAY")].sort_values().dropna().index[-1]
    lddt_by_pdb_id.loc[(_pdb_id, _sel)]
    return


@app.cell
def _(lddt_by_pdb_id):
    lddt_by_pdb_id.query('protein in ("5sop", "8z76", "7z0e", "5rbr", "7frd")')
    return


@app.cell
def _(lddt_by_pdb_id):
    (
        lddt_by_pdb_id.avg_silhouette_to_ref.max(axis=1)
        - lddt_by_pdb_id.avg_silhouette_to_ref.min(axis=1)
    ).sort_values()
    return


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
    lddt_merged
    return (lddt_merged,)


@app.cell
def _(lddt_merged, make_ridge_plot, save_fig):
    _g = make_ridge_plot(lddt_merged, "delta_silhouette_to_ref", is_delta=True)
    save_fig(_g, "13_delta_lddt_to_ref_distribution.png")
    return


if __name__ == "__main__":
    app.run()
