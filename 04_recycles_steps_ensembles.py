import marimo

__generated_with = "0.23.2"
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
    from loaders import load_rscc, load_rscc_recyc_diff
    from style import (
        MODEL_DISPLAY,
        MODEL_PALETTE,
        OPTIMAL_GW,
        save_fig,
    )

    return (
        MODEL_DISPLAY,
        MODEL_PALETTE,
        OPTIMAL_GW,
        load_rscc,
        load_rscc_recyc_diff,
        np,
        pd,
        plt,
        save_fig,
    )


@app.cell
def _(mo):
    mo.md(
        r"""
        # 04 — Recycles, diffusion steps, ensemble size

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
        """
    )
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
def _(default_df, grid_df, pd):
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
    MODEL_DISPLAY,
    MODEL_PALETTE,
    MODEL_ORDER,
    auto_save,
    combined_df,
    plt,
    save_fig,
):
    def _agg(_df, _group):
        return (
            _df.groupby(["model_name", _group])
            .rscc.agg(mean="mean", sd="std", n="count")
            .reset_index()
        )

    # Recycles panel: hold diffusion steps at DEFAULT_STEPS.
    _rec = _agg(
        combined_df[combined_df.n_diffusion_steps == DEFAULT_STEPS], "n_recycles"
    )
    # Steps panel: hold recycles at each model's default.
    _def_rec = combined_df.model_name.map(DEFAULT_RECYCLES)
    _stp = _agg(combined_df[combined_df.n_recycles == _def_rec], "n_diffusion_steps")

    _panels = [
        ("n_recycles", _rec, "Number of recycles", f"RSCC vs recycles (steps = {DEFAULT_STEPS})"),
        ("n_diffusion_steps", _stp, "Number of diffusion steps", "RSCC vs diffusion steps (recycles = default)"),
    ]

    fig, _axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    for _ax, (_xcol, _data, _xlabel, _title) in zip(_axes, _panels):
        # Evenly spaced categorical x positions, labelled by the real values.
        _xvals = sorted(_data[_xcol].unique())
        _pos = {v: i for i, v in enumerate(_xvals)}
        for _m in MODEL_ORDER:
            _sub = _data[_data.model_name == _m].sort_values(_xcol)
            if _sub.empty:
                continue
            _x = [_pos[v] for v in _sub[_xcol]]
            _y = _sub["mean"].to_numpy()
            _e = _sub["sd"].to_numpy()
            _color = MODEL_PALETTE[_m]
            _ax.plot(_x, _y, marker="o", color=_color, label=MODEL_DISPLAY[_m], zorder=3)
            if len(_x) > 1:
                _ax.fill_between(_x, _y - _e, _y + _e, color=_color, alpha=0.2, lw=0, zorder=2)
        _ax.set_xticks(list(_pos.values()))
        _ax.set_xticklabels([str(v) for v in _xvals])
        _ax.set_xlabel(_xlabel)
        _ax.set_title(_title)

    _axes[0].set_ylabel("Mean RSCC ± SD")
    _axes[0].legend(title="Model", frameon=False, loc="lower right")
    fig.suptitle("5MHX — RSCC vs recycles and diffusion steps")
    fig.tight_layout()

    if auto_save.value:
        save_fig(fig, "recycles_steps/5mhx_rscc_vs_recycles_steps.svg")
    fig
    return (fig,)


@app.cell
def _(mo):
    mo.md(
        r"""
        ## Ensemble size (6B8X)

        Mean RSCC ± SD across the **32 altloc selections** in **6B8X** as a function of
        **ensemble size** (`{1, 2, 4, 8, 16}`), one line per model at its optimal
        guidance weight. Source: `rscc_results_6b8x_ensemble_size.csv`, loaded via
        `load_rscc` (it already keeps `ensemble_size`). Every model is already at its
        optimal gw, so — unlike the recycles/steps figure — there is no default-point
        merge: the figure is single-source.

        - x-axis: ensemble size on **even (log2) spacing** — powers of two get equal width.
        - Spread band: mean ± 1 SD across the 32 selections (so the band is over
          selections, not over configs).
        """
    )
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
    MODEL_DISPLAY,
    MODEL_PALETTE,
    MODEL_ORDER,
    auto_save,
    ensemble_df,
    np,
    plt,
    save_fig,
):
    _agg = (
        ensemble_df.groupby(["model_name", "ensemble_size"])
        .rscc.agg(mean="mean", sd="std", n="count")
        .reset_index()
    )

    # Powers-of-two sizes -> even spacing == log2 spacing. Plot at log2(size) and
    # label the ticks with the real ensemble sizes.
    _sizes = sorted(_agg.ensemble_size.unique())
    _pos = {v: np.log2(v) for v in _sizes}

    # Dodge each model off the shared tick so the ±SD error bars stay legible
    # instead of overlapping into a muddy band (the spread is ~4x the signal range).
    _n = len(MODEL_ORDER)
    _dodge = {m: (i - (_n - 1) / 2) * 0.06 for i, m in enumerate(MODEL_ORDER)}

    ens_fig, _ax = plt.subplots(figsize=(7, 5))
    for _m in MODEL_ORDER:
        _sub = _agg[_agg.model_name == _m].sort_values("ensemble_size")
        if _sub.empty:
            continue
        _x = np.array([_pos[v] for v in _sub.ensemble_size]) + _dodge[_m]
        _y = _sub["mean"].to_numpy()
        _e = _sub["sd"].to_numpy()
        _ax.errorbar(
            _x, _y, yerr=_e, marker="o", color=MODEL_PALETTE[_m],
            label=MODEL_DISPLAY[_m], capsize=3, lw=1.6, elinewidth=1.2, zorder=3,
        )
    _ax.set_xticks([_pos[v] for v in _sizes])
    _ax.set_xticklabels([str(v) for v in _sizes])
    _ax.set_xlabel("Ensemble size")
    _ax.set_ylabel("Mean RSCC ± SD")
    _ax.set_title("6B8X — RSCC vs ensemble size")
    _ax.legend(title="Model", frameon=False, loc="lower right")
    ens_fig.tight_layout()

    if auto_save.value:
        save_fig(ens_fig, "ensemble_size/6b8x_rscc_vs_ensemble_size.svg")
    ens_fig
    return (ens_fig,)


if __name__ == "__main__":
    app.run()
