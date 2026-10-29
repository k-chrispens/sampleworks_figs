import marimo

__generated_with = "0.23.2"
app = marimo.App(width="full")


@app.cell
def _(mo):
    from supplementary_plot import linked_plot, inspection_tables, reconstruction_bundle
    from plot_budgets import make_budget_figure
    get_inspection, set_inspection = mo.state({})
    return get_inspection, set_inspection, linked_plot, inspection_tables, reconstruction_bundle, make_budget_figure


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell
def _():
    import numpy as np
    import pandas as pd
    from loaders import load_rscc, load_rscc_recyc_diff
    from style import (
        OPTIMAL_GW,
        save_fig,
    )

    return (
        OPTIMAL_GW,
        load_rscc,
        load_rscc_recyc_diff,
        np,
        pd,
        save_fig,
    )


@app.cell
def _(mo):
    mo.md(r"""
    # 08 — Recycles, diffusion steps, ensemble size

    Sampling-budget sweeps: how mean RSCC responds to the number of **recycles**,
    **diffusion steps** (both **5MHX**), and **ensemble size** (**6B8X**), each at
    every model's optimal guidance weight (`0.1` for the three main models, `0.01`
    for RF3). The recycles/steps grid sweep (`rscc_results_5mhx_recyc_diff.csv`)
    covers recycles `{1, 3, 10}` × steps `{50, 100, 200, 400}` per model; the
    ensemble-size sweep (`rscc_results_6b8x_ensemble_size.csv`) covers sizes
    `{1, 2, 4, 8, 16}` per model.

    | Figure | Output path |
    |--------|-------------|
    | RSCC vs recycles (steps=200) ‖ RSCC vs diffusion steps (recycles=default) | `figures/recycles_steps/5mhx_rscc_vs_recycles_steps.svg` |
    | RSCC vs ensemble size | `figures/ensemble_size/6b8x_rscc_vs_ensemble_size.svg` |

    - **Recycles panel** holds diffusion steps at **200**; x = `{1, 3, 10}`.
    - **Steps panel** holds recycles at each model's production default
      (`boltz2 = 3`, `protenix`/`rf3 = 10`); x = `{50, 100, 200, 400}`.
    - The production **default** config (from the `rscc_results_max_rmsd_*` CSVs) is
      merged only for `(model, recycles, steps)` cells the grid does **not** cover. With
      the current full grid it contributes **0** cells (the combine cell prints this), so
      the figure is grid-only and self-consistent.
    """)
    return


@app.cell
def _(mo):
    auto_save = mo.ui.switch(label="Auto-save figures", value=False)
    auto_save
    return (auto_save,)


@app.cell
def _():
    # Config the default (max-RMSD / occ-sweep) run represents per model — it carries no
    # recycles/diffusion-steps columns, so we pin it here. Recycles axis then reads {1,3,10}.
    DEFAULT_RECYCLES = {"boltz2MD": 3, "boltz2X-RAY": 3, "protenix": 10, "rf3": 10}
    DEFAULT_STEPS = 200
    PROTEIN = "5MHX"
    MODEL_ORDER = ["boltz2MD", "boltz2X-RAY", "protenix", "rf3"]
    return DEFAULT_RECYCLES, DEFAULT_STEPS, MODEL_ORDER, PROTEIN


@app.cell
def _(PROTEIN, load_rscc_recyc_diff):
    _grid = load_rscc_recyc_diff("results/rscc_results_5mhx_recyc_diff.csv")
    _grid = _grid[_grid.protein.str.upper() == PROTEIN]
    grid_df = _grid[_grid.selection_name != "all altloc"].copy()
    return (grid_df,)


@app.cell
def _(DEFAULT_RECYCLES, DEFAULT_STEPS, OPTIMAL_GW, PROTEIN, load_rscc, np, pd):
    # The default point for each model = its optimal-gw max-RMSD run, placed at
    # (DEFAULT_RECYCLES[model], DEFAULT_STEPS).
    _d = pd.concat(
        [
            load_rscc("results/rscc_results_max_rmsd_boltz_protenix.csv"),
            load_rscc("results/rscc_results_max_rmsd_rf3_fixed.csv"),
        ],
        ignore_index=True,
    )
    _d = _d[_d.protein.str.upper() == PROTEIN]
    _d = _d[_d.selection_name != "all altloc"]
    _d = _d[_d.model_name.isin(DEFAULT_RECYCLES)]
    _d = _d[np.isclose(_d.guidance_weight, _d.model_name.map(OPTIMAL_GW))].copy()
    _d["n_recycles"] = _d.model_name.map(DEFAULT_RECYCLES)
    _d["n_diffusion_steps"] = DEFAULT_STEPS
    default_df = _d
    return (default_df,)


@app.cell
def combine_budget_grid(default_df, grid_df, pd):
    _cols = [
        "protein",
        "model_name",
        "guidance_weight",
        "n_recycles",
        "n_diffusion_steps",
        "rscc",
        "selection",
    ]
    # "Missing points come from the default": only add a default (model, recycles, steps)
    # cell when the grid doesn't already cover it. Today the grid has no rows at the default
    # recycle count, so the default fills it; after the grid is re-run with those recycles,
    # the grid's own cell wins (no double-count / cross-pipeline blend at steps=200).
    _grid_keys = set(
        zip(grid_df.model_name, grid_df.n_recycles, grid_df.n_diffusion_steps)
    )
    _def = default_df[
        ~default_df.apply(
            lambda r: (r.model_name, r.n_recycles, r.n_diffusion_steps) in _grid_keys,
            axis=1,
        )
    ]
    print(
        f"combine: {len(_def) // 20 if len(_def) else 0} default model-cells added "
        f"({len(default_df) - len(_def)} already covered by grid)"
    )
    combined_df = pd.concat([grid_df[_cols], _def[_cols]], ignore_index=True)
    return (combined_df,)


@app.cell
def _(
    DEFAULT_RECYCLES,
    DEFAULT_STEPS,
    MODEL_ORDER,
    auto_save,
    combined_df,
    make_budget_figure,
    linked_plot,
    set_inspection,
    mo,
    save_fig,
):
    fig = make_budget_figure(combined_df, MODEL_ORDER, default_recycles=DEFAULT_RECYCLES, default_steps=DEFAULT_STEPS)
    budget_panel = linked_plot(fig, mo, set_inspection, "recycles-steps")
    if auto_save.value:
        save_fig(fig, "recycles_steps/5mhx_rscc_vs_recycles_steps.svg")
    budget_panel.view
    return fig, budget_panel


@app.cell
def _(mo):
    mo.md(r"""
    ## Ensemble size (6B8X)

    Mean RSCC ± SD across the **32 altloc selections** in **6B8X** as a function of
    **ensemble size** (`{1, 2, 4, 8, 16}`), one line per model at its optimal
    guidance weight. Source: `rscc_results_6b8x_ensemble_size.csv`, loaded via
    `load_rscc` (it already keeps `ensemble_size`). Every model is already at its
    optimal gw, so — unlike the recycles/steps figure — there is no default-point
    merge: the figure is single-source.

    - x-axis: ensemble size on **even (log2) spacing** — powers of two get equal width.
    - Error bars: mean ± 1 SD across selections, not across configurations.
    """)
    return


@app.cell
def _(OPTIMAL_GW, load_rscc, np):
    # 6B8X ensemble-size sweep at each model's optimal gw. load_rscc already keeps
    # ensemble_size; the CSV is single-protein, single-gw-per-model, so the three
    # filters below are defensive guards against future rows, not load-bearing today.
    _e = load_rscc("results/rscc_results_6b8x_ensemble_size.csv")
    _e = _e[_e.protein.str.upper() == "6B8X"]
    _e = _e[_e.selection_name != "all altloc"]
    _e = _e[np.isclose(_e.guidance_weight, _e.model_name.map(OPTIMAL_GW))]
    ensemble_df = _e.copy()
    return (ensemble_df,)


@app.cell
def _(
    MODEL_ORDER,
    auto_save,
    ensemble_df,
    make_budget_figure,
    linked_plot,
    set_inspection,
    mo,
    save_fig,
):
    ens_fig = make_budget_figure(ensemble_df, MODEL_ORDER)
    ensemble_panel = linked_plot(ens_fig, mo, set_inspection, "ensemble-size")
    if auto_save.value:
        save_fig(ens_fig, "ensemble_size/6b8x_rscc_vs_ensemble_size.svg")
    ensemble_panel.view
    return ens_fig, ensemble_panel


@app.cell
def _(budget_panel, ensemble_panel):
    panels = {"recycles-steps": budget_panel, "ensemble-size": ensemble_panel}
    return (panels,)


@app.cell
def _(panels, get_inspection, inspection_tables, mo):
    inspection = inspection_tables(panels, get_inspection(), mo)
    inspection
    return (inspection,)


@app.cell
def _(panels, mo):
    export_panel = mo.ui.dropdown(options=list(panels), value=next(iter(panels)), label="Fixed reconstruction")
    build_export = mo.ui.run_button(label="Build SVG + data bundle")
    mo.hstack([export_panel, build_export])
    return export_panel, build_export


@app.cell
def _(build_export, export_panel, panels, reconstruction_bundle, mo):
    reconstruction_download = mo.download(
        data=reconstruction_bundle(panels[export_panel.value].figure),
        filename=f"08-{export_panel.value}-reconstruction.zip", label="Download reconstruction",
    ) if build_export.value else None
    reconstruction_download
    return (reconstruction_download,)


if __name__ == "__main__":
    app.run()
