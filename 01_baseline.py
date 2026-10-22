import marimo

__generated_with = "0.23.2"
app = marimo.App()


@app.cell
def _():
    import io
    import json
    import os
    from pathlib import Path

    import marimo as mo
    import matplotlib.pyplot as plt

    from baseline import (
        baseline_ecdf,
        baseline_provenance,
        load_baseline_snapshot,
        make_baseline_interactive,
        make_baseline_static,
        prepare_baseline,
    )
    from plot_interaction import highlight_table_rows, make_linked_widget, select_observations
    from style import MODELS_ALT

    return (
        MODELS_ALT, Path, baseline_ecdf, baseline_provenance,
        io, json, load_baseline_snapshot, make_baseline_interactive,
        make_baseline_static, mo, os, plt, prepare_baseline, select_observations,
        highlight_table_rows, make_linked_widget,
    )


@app.cell
def _(mo, os):
    _intro = mo.md("""
    # Figure 2C — baseline fit and conformational recovery

    Reconstruct the existing results first, then inspect the underlying observations.
    A curve farther right indicates higher scores. These are exact empirical CDFs,
    matching the existing saved baseline figure's ECDF mode—not smoothed estimates.

    **Reconstruction, not a verified historical paper snapshot.** The current notebook
    uses protein-mean RSCC; the paper caption describes regions. Both are available
    explicitly below. Clustering retains the current score formula, including its
    treatment of occupancy arrays with more than two entries.

    Inputs are read-only. Downloads are generated in memory; nothing writes to results.

    | Figure | Output |
    |---|---|
    | Baseline RSCC and clustering ECDFs | Plot-toolbar SVG for the visible view; `fig2c-reconstruction.svg` for main models |
    """)
    source_path = mo.ui.text(
        value=os.environ.get("SAMPLEWORKS_RESULTS", "results"),
        label="Results folder",
        full_width=True,
    ).form(submit_button_label="Load snapshot")
    mo.vstack([_intro, source_path])
    return (source_path,)


@app.cell
def _(Path, load_baseline_snapshot, mo, os, source_path):
    _path = source_path.value
    if _path is None and not mo.running_in_notebook():
        _path = os.environ.get("SAMPLEWORKS_RESULTS", "results")
    mo.stop(_path is None, mo.md("Choose the results folder and load the snapshot."))
    try:
        rscc_baseline, lddt_baseline, source_manifest = load_baseline_snapshot(
            Path(_path)
        )
        _error = None
    except (OSError, ValueError, KeyError) as _exc:
        _error = str(_exc)
    mo.stop(_error is not None, mo.callout(f"Snapshot could not be loaded: {_error}", kind="warn"))
    return lddt_baseline, rscc_baseline, source_manifest


@app.cell
def _(lddt_baseline, mo, rscc_baseline, source_manifest):
    model_ids = sorted(set(rscc_baseline.model_name) | set(lddt_baseline.model_name))
    rscc_unit = mo.ui.radio(
        options={
            "Protein mean — reconstruct notebook 02": "protein",
            "Region — explore the paper caption's population": "region",
        },
        value="Protein mean — reconstruct notebook 02",
        label="RSCC observation",
    )
    mo.vstack([
        mo.md(f"Loaded snapshot: `{source_manifest['snapshot_sha256'][:16]}`. "
              "Changing controls does not reload files."),
        rscc_unit,
    ])
    return model_ids, rscc_unit


@app.cell
def _(
    baseline_ecdf, baseline_provenance, lddt_baseline,
    model_ids, prepare_baseline, rscc_baseline, rscc_unit, source_manifest,
):
    observations, contributors, coverage = prepare_baseline(
        rscc_baseline, lddt_baseline, model_ids, rscc_unit.value
    )
    curves = baseline_ecdf(observations)
    provenance = baseline_provenance(source_manifest, model_ids, rscc_unit.value)
    return contributors, coverage, curves, observations, provenance


@app.cell
def _(curves, make_baseline_interactive, make_linked_widget, mo):
    mo.stop(curves.empty, mo.callout("No finite baseline observations.", kind="warn"))
    base_figure = make_baseline_interactive(curves)
    chart_point_ids = [[_row[0] for _row in _trace.customdata] for _trace in base_figure.data]
    get_plot_selection, _set_plot_selection = mo.state([])
    figure_widget, clear_plot_selection = make_linked_widget(base_figure, _set_plot_selection, "fig2c-visible")
    return base_figure, chart_point_ids, clear_plot_selection, figure_widget, get_plot_selection


@app.cell
def _(clear_plot_selection, figure_widget, mo):
    chart = mo.ui.anywidget(figure_widget)
    clear_button = mo.ui.button(label="Clear plot selection", on_click=lambda _: clear_plot_selection())
    mo.vstack([
        chart,
        clear_button,
        mo.md("**Click a model in the legend to show/hide it in both panels; "
              "double-click to isolate it.** Variants start hidden but remain in the legend. "
              "Hover for values or box/lasso-select observations. Select table rows "
              "to add black inspection outlines; clear those rows to remove them. "
              "Use the plot's camera button to save the **visible view as SVG**. "
              "The initial ECDF anchor is not an observation."),
    ])
    return chart, clear_button


@app.cell
def _(base_figure, chart_point_ids, figure_widget, get_plot_selection, highlight_table_rows, mo, observations, select_observations):
    _shown = select_observations(observations, chart_point_ids, get_plot_selection())
    highlight_table_rows(figure_widget, base_figure, observations.iloc[:0])
    observation_table = mo.ui.table(
        _shown, selection="multi", label="Loaded observations (select rows to inspect source records)",
        page_size=10, on_change=lambda rows: highlight_table_rows(figure_widget, base_figure, rows),
    )
    mo.vstack([
        mo.md(f"## Observations\n{len(_shown)} observations shown. "
              "Without a chart selection this table includes all loaded models, "
              "including legend-hidden ones."),
        observation_table,
    ])
    return (observation_table,)


@app.cell
def _(contributors, mo, observation_table):
    _selected = observation_table.value
    _table = mo.ui.table(
        contributors.loc[contributors.point_id.isin(_selected.point_id)]
        if len(_selected) else contributors.iloc[:0],
        selection=None, page_size=10,
    )
    mo.vstack([mo.md("### Source records for selected observations"), _table])
    return


@app.cell
def _(coverage, mo, provenance):
    mo.vstack([
        mo.md("## Coverage and settings\nMissing source values are counted before plotting. "
              "Protein means retain notebook 02's NaN-skipping behavior. "
              "No missing value is filled from another metric file."),
        mo.ui.table(coverage, selection=None),
        mo.accordion({"Recipe and source checksums": mo.json(provenance)}),
    ])
    return


@app.cell
def _(MODELS_ALT, curves, io, json, make_baseline_static, mo, observations, plt, provenance):
    mo.stop(curves.empty, mo.md("No observations to export."))
    _figure = make_baseline_static(curves.loc[~curves.model_name.isin(MODELS_ALT)])
    _svg = io.BytesIO()
    _figure.savefig(_svg, format="svg", bbox_inches="tight")
    plt.close(_figure)
    mo.hstack([
        mo.download(_svg.getvalue(), filename="fig2c-reconstruction.svg",
                    label="Download main-model reconstruction SVG"),
        mo.download(observations.to_csv(index=False).encode(), filename="fig2c-observations.csv",
                    label="Download all-model values"),
        mo.download(curves.to_csv(index=False).encode(), filename="fig2c-curves.csv",
                    label="Download ECDF vertices"),
        mo.download(json.dumps(provenance, indent=2).encode(), filename="fig2c-provenance.json",
                    label="Download recipe and checksums"),
    ], wrap=True)
    return


if __name__ == "__main__":
    app.run()
