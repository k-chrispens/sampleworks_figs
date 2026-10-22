import marimo

__generated_with = "0.23.2"
app = marimo.App(width="full")


@app.cell
def _():
    import io
    import json
    import os
    import zipfile

    import marimo as mo
    import matplotlib.pyplot as plt

    from min_rmsd import (
        highlighted_min_rmsd, load_min_rmsd_snapshot, make_min_rmsd_interactive,
        make_min_rmsd_static, min_rmsd_provenance, prepare_min_rmsd,
    )
    from plot_interaction import highlight_table_rows, make_linked_widget, select_observations

    return (
        highlighted_min_rmsd, io, json, load_min_rmsd_snapshot, make_min_rmsd_interactive,
        make_min_rmsd_static, min_rmsd_provenance, mo, os, plt, prepare_min_rmsd,
        select_observations, zipfile, highlight_table_rows, make_linked_widget,
    )


@app.cell
def _(mo, os):
    source_path = mo.ui.text(
        value=os.environ.get("SAMPLEWORKS_RESULTS", "results"), label="Results folder",
        full_width=True,
    ).form(submit_button_label="Load snapshot")
    aggregation_control = mo.ui.radio(
        options={"Median — existing notebook reconstruction": "median",
                 "Mean — caption-described alternative, historical computation unverified": "mean"},
        value="Median — existing notebook reconstruction", label="Across-model aggregation",
    )
    mo.vstack([
        mo.md("""
        # Figure 3E — minimum RMSD to each altloc

        **First, within each prediction model:** the evaluation has already taken
        the minimum RMSD over its eight ensemble structures, separately for A/B.
        The best-matching member can differ between A and B.

        **Then, across the four main prediction models:** combine those minima into
        one coordinate per region. The existing notebook uses a **median**; the
        paper caption says **mean**. Compare them below without changing the
        within-ensemble minimum. Neither choice verifies the historical paper recipe.

        Require at least three complete model pairs at their existing chosen guidance
        weights. The full-range and 0–3 Å panels show the same prepared coordinates.
        The diagonal means equal distance to A/B, not a recovery threshold.

        | Figure | Output |
        |---|---|
        | Minimum RMSD to altloc A/B | Visible-view toolbar SVG; two-panel reconstruction ZIP with circled-region table |

        Inputs are read-only; no structure evaluation is run here.
        """),
        source_path, aggregation_control,
    ])
    return aggregation_control, source_path


@app.cell
def _(load_min_rmsd_snapshot, mo, os, source_path):
    _path = source_path.value
    if _path is None and not mo.running_in_notebook():
        _path = os.environ.get("SAMPLEWORKS_RESULTS", "results")
    mo.stop(_path is None, mo.md("Choose the results folder and load the snapshot."))
    try:
        measurements, manifest = load_min_rmsd_snapshot(_path)
        _error = None
    except (OSError, ValueError, KeyError) as _exc:
        _error = str(_exc)
    mo.stop(_error is not None, mo.callout(f"Snapshot could not be loaded: {_error}", kind="warn"))
    return manifest, measurements


@app.cell
def _(aggregation_control, manifest, measurements, min_rmsd_provenance, prepare_min_rmsd):
    aggregation = aggregation_control.value
    points, contributions, coverage = prepare_min_rmsd(measurements, aggregation)
    provenance = min_rmsd_provenance(manifest, aggregation)
    return aggregation, contributions, coverage, points, provenance


@app.cell
def _(aggregation, make_min_rmsd_interactive, make_linked_widget, mo, points):
    mo.stop(not points.plotted.any(), mo.callout("No regions have three complete model pairs.", kind="warn"))
    base_figure = make_min_rmsd_interactive(points, aggregation)
    chart_point_ids = [[_row[0] for _row in _trace.customdata] for _trace in base_figure.data]
    get_plot_selection, _set_plot_selection = mo.state([])
    figure_widget, clear_plot_selection = make_linked_widget(
        base_figure, _set_plot_selection, f"fig3e-{aggregation}-visible"
    )
    return base_figure, chart_point_ids, clear_plot_selection, figure_widget, get_plot_selection


@app.cell
def _(aggregation, clear_plot_selection, figure_widget, mo, provenance):
    chart = mo.ui.anywidget(figure_widget)
    clear_button = mo.ui.button(label="Clear plot selection", on_click=lambda _: clear_plot_selection())
    mo.vstack([
        mo.md(f"Loaded snapshot: `{provenance['snapshot_sha256'][:16]}`. "
              f"Across-model aggregation: **{aggregation}**. Legend controls display, not model membership."),
        chart,
        clear_button,
        mo.md("Hover for identities; box/lasso-select either panel and inspect table rows. "
              "The circles follow the existing RSCC-matched categories, with red examples "
              "in the zoom panel only. There are no fixed segment text labels. Table-selected "
              "regions get black outlines in both panels; clearing rows removes them."),
    ])
    return chart, clear_button


@app.cell
def _(base_figure, chart_point_ids, figure_widget, get_plot_selection, highlight_table_rows, mo, points, select_observations):
    _shown = select_observations(points.loc[points.plotted], chart_point_ids, get_plot_selection())
    highlight_table_rows(figure_widget, base_figure, points.iloc[:0])
    region_table = mo.ui.table(
        _shown, selection="multi", page_size=10, label="Region summaries",
        on_change=lambda rows: highlight_table_rows(figure_widget, base_figure, rows),
    )
    mo.vstack([mo.md(f"## Regions\n{len(_shown)} regions shown; clear selection to restore all."),
               region_table])
    return (region_table,)


@app.cell
def _(contributions, mo, region_table):
    _selected = region_table.value
    _rows = contributions.loc[
        contributions.point_id.isin(_selected.point_id)
    ] if len(_selected) else contributions.iloc[:0]
    mo.vstack([
        mo.md("### Contributing model minima\nPrecomputed ensemble minima, guidance weights, "
              "source record numbers, and inclusion status. No minimum is recomputed here."),
        mo.ui.table(_rows, selection=None, page_size=10),
    ])
    return


@app.cell
def _(contributions, coverage, mo, points, provenance):
    mo.vstack([
        mo.md("## Coverage and recipe\nA model must have finite A **and** B minima to "
              "contribute. Recorded failures and absent chosen-weight rows are distinct. "
              "Coverage does not presume an unexported complete experimental grid."),
        mo.ui.table(coverage, selection=None),
        mo.ui.table(contributions.loc[~contributions.included], selection=None, label="Excluded model pairs"),
        mo.ui.table(points.loc[~points.plotted], selection=None, label="Below minimum model coverage"),
        mo.accordion({"Source checksums and reconstruction recipe": mo.json(provenance)}),
    ])
    return


@app.cell
def _(
    aggregation, contributions, highlighted_min_rmsd, io, json, make_min_rmsd_static,
    mo, plt, points, provenance, zipfile,
):
    mo.stop(not points.plotted.any(), mo.md("No reconstruction to export."))
    _figure = make_min_rmsd_static(points, aggregation)
    _svg = io.BytesIO()
    _figure.savefig(_svg, format="svg", bbox_inches="tight")
    plt.close(_figure)
    _highlights = highlighted_min_rmsd(points)
    _csv = _highlights.to_csv(index=False).encode()
    _recipe = json.dumps(provenance, indent=2).encode()
    _prefix = f"fig3e-{aggregation}"
    _bundle = io.BytesIO()
    with zipfile.ZipFile(_bundle, "w", compression=zipfile.ZIP_DEFLATED) as _archive:
        _archive.writestr(f"{_prefix}-reconstruction.svg", _svg.getvalue())
        _archive.writestr(f"{_prefix}-highlighted-regions.csv", _csv)
        _archive.writestr(f"{_prefix}-provenance.json", _recipe)
    mo.vstack([
        mo.md("## Circled regions and downloads\nThe table identifies every visible circled "
              "region and its categories/colors, including the zoom-only red examples. "
              "ZIP/SVG/CSV downloads use the selected aggregation and full region population; "
              "chart selections and legend toggles do not filter them."),
        mo.ui.table(_highlights, selection=None, page_size=10, label="Circled regions"),
        mo.hstack([
            mo.download(_bundle.getvalue(), filename=f"{_prefix}-reconstruction.zip", label="Download reconstruction bundle"),
            mo.download(_svg.getvalue(), filename=f"{_prefix}-reconstruction.svg", label="Download full SVG"),
            mo.download(_csv, filename=f"{_prefix}-highlighted-regions.csv", label="Download circled-region table"),
            mo.download(points.to_csv(index=False).encode(), filename=f"{_prefix}-regions.csv", label="Download regions"),
            mo.download(contributions.to_csv(index=False).encode(), filename=f"{_prefix}-models.csv", label="Download model contributions"),
            mo.download(_recipe, filename=f"{_prefix}-provenance.json", label="Download recipe"),
        ], wrap=True),
    ])
    return


if __name__ == "__main__":
    app.run()
