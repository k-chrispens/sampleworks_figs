import marimo

__generated_with = "0.23.2"
app = marimo.App()


@app.cell
def _():
    import io
    import json
    import os
    import zipfile

    import marimo as mo
    import matplotlib.pyplot as plt

    from plot_interaction import highlight_table_rows, make_linked_widget, select_observations
    from region_improvement import (
        highlighted_regions, improvement_provenance, load_improvement_snapshot,
        make_improvement_interactive, make_improvement_static,
        prepare_density, prepare_improvement,
    )

    return (
        highlighted_regions, improvement_provenance, io, json, load_improvement_snapshot,
        make_improvement_interactive, make_improvement_static, mo, os, plt,
        prepare_density, prepare_improvement, select_observations, zipfile,
        highlight_table_rows, make_linked_widget,
    )


@app.cell
def _(mo, os):
    source_path = mo.ui.text(
        value=os.environ.get("SAMPLEWORKS_RESULTS", "results"),
        label="Results folder", full_width=True,
    ).form(submit_button_label="Load snapshot")
    mo.vstack([
        mo.md("""
        # Figure 3D — clustering-score improvement

        Each point is one region. Pair baseline with the existing chosen guidance
        weight for each main model, then subtract the **marginal medians**.
        At least three paired model rows are required, as in notebook 02.

        **Existing score:** minimum predicted occupancy × silhouette-to-reference
        / minimum map occupancy. Inputs here are the 50/50 A/B map, eight-member
        ensembles. The minimum uses **all entries** in the predicted occupancy
        array, even when there are more than two. This reconstructs the current
        calculation; state correspondence remains an open scientific question.

        | Figure | Output |
        |---|---|
        | Aggregated clustering-score improvement | Visible-view toolbar SVG; full reconstruction ZIP with circled-region table |

        No structural evaluation is run and source files are never written.
        The exact historical paper snapshot remains unverified.
        """),
        source_path,
    ])
    return (source_path,)


@app.cell
def _(
    improvement_provenance, load_improvement_snapshot, mo, os, prepare_density,
    prepare_improvement, source_path,
):
    _path = source_path.value
    if _path is None and not mo.running_in_notebook():
        _path = os.environ.get("SAMPLEWORKS_RESULTS", "results")
    mo.stop(_path is None, mo.md("Choose the results folder and load the snapshot."))
    try:
        _data, _manifest = load_improvement_snapshot(_path, "clustering_score")
        points, contributions, coverage = prepare_improvement(_data, "clustering_score")
        density, density_coverage = prepare_density(points)
        provenance = improvement_provenance(_manifest, "clustering_score")
        _error = None
    except (OSError, ValueError, KeyError) as _exc:
        _error = str(_exc)
    mo.stop(_error is not None, mo.callout(f"Snapshot could not be loaded: {_error}", kind="warn"))
    return contributions, coverage, density, density_coverage, points, provenance


@app.cell
def _(density, make_improvement_interactive, make_linked_widget, mo, points):
    mo.stop(not points.plotted.any(), mo.callout("No regions meet the plotting requirements.", kind="warn"))
    base_figure = make_improvement_interactive(points, density, "clustering_score")
    chart_point_ids = [[_row[0] for _row in _trace.customdata] for _trace in base_figure.data]
    get_plot_selection, _set_plot_selection = mo.state([])
    figure_widget, clear_plot_selection = make_linked_widget(base_figure, _set_plot_selection, "fig3d-visible")
    return base_figure, chart_point_ids, clear_plot_selection, figure_widget, get_plot_selection


@app.cell
def _(clear_plot_selection, figure_widget, mo, provenance):
    chart = mo.ui.anywidget(figure_widget)
    clear_button = mo.ui.button(label="Clear plot selection", on_click=lambda _: clear_plot_selection())
    mo.vstack([
        mo.md(f"Loaded snapshot: `{provenance['snapshot_sha256'][:16]}`. "
              "**The legend controls display categories, not the fixed contributing models.**"),
        chart,
        clear_button,
        mo.md("Hover for region identity. Select scatter points, then table rows, to inspect "
              "the underlying occupancy/silhouette values and original CSV records. "
              "Segment identities are not printed over the circles. Table-selected regions "
              "get black outlines; clearing rows removes them."),
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
    mo.vstack([mo.md(f"## Regions\n{len(_shown)} regions shown. Clear selection to restore all."),
               region_table])
    return (region_table,)


@app.cell
def _(contributions, mo, region_table):
    _selected = region_table.value
    _rows = contributions.loc[
        contributions.point_id.isin(_selected.point_id)
    ] if len(_selected) else contributions.iloc[:0]
    mo.vstack([
        mo.md("### Contributing model pairs\nValues include the score's occupancy/silhouette "
              "components, chosen guidance weights, and original CSV record numbers."),
        mo.ui.table(_rows, selection=None, page_size=10),
    ])
    return


@app.cell
def _(coverage, density_coverage, mo, points, provenance):
    mo.vstack([
        mo.md("## Coverage and recipe\nThe legacy threshold counts paired rows, not finite "
              "pairs. Missing pairs, undrawn regions, and marginal KDE populations are explicit."),
        mo.ui.table(coverage, selection=None),
        mo.ui.table(points.loc[~points.plotted], selection=None, label="Undrawn regions"),
        mo.ui.table(density_coverage, selection=None),
        mo.accordion({"Source checksums and reconstruction recipe": mo.json(provenance)}),
    ])
    return


@app.cell
def _(
    contributions, density, highlighted_regions, io, json, make_improvement_static,
    mo, plt, points, provenance, zipfile,
):
    mo.stop(not points.plotted.any(), mo.md("No reconstruction to export."))
    _figure = make_improvement_static(points, density, "clustering_score")
    _svg = io.BytesIO()
    _figure.savefig(_svg, format="svg", bbox_inches="tight")
    plt.close(_figure)
    _highlights = highlighted_regions(points)
    _csv = _highlights.to_csv(index=False).encode()
    _recipe = json.dumps(provenance, indent=2).encode()
    _bundle = io.BytesIO()
    with zipfile.ZipFile(_bundle, "w", compression=zipfile.ZIP_DEFLATED) as _archive:
        _archive.writestr("fig3d-reconstruction.svg", _svg.getvalue())
        _archive.writestr("fig3d-highlighted-regions.csv", _csv)
        _archive.writestr("fig3d-provenance.json", _recipe)
    mo.vstack([
        mo.md("## Circled regions and downloads\nThe full reconstruction ZIP contains its "
              "unlabeled SVG, companion table with all highlight categories/colors, and recipe. "
              "These downloads retain the full analysis, independent of transient selections."),
        mo.ui.table(_highlights, selection=None, page_size=10, label="Circled regions"),
        mo.hstack([
            mo.download(_bundle.getvalue(), filename="fig3d-reconstruction.zip", label="Download reconstruction bundle"),
            mo.download(_svg.getvalue(), filename="fig3d-reconstruction.svg", label="Download full SVG"),
            mo.download(_csv, filename="fig3d-highlighted-regions.csv", label="Download circled-region table"),
            mo.download(points.to_csv(index=False).encode(), filename="fig3d-regions.csv", label="Download regions"),
            mo.download(contributions.to_csv(index=False).encode(), filename="fig3d-model-pairs.csv", label="Download model pairs"),
            mo.download(density.to_csv(index=False).encode(), filename="fig3d-density.csv", label="Download density"),
            mo.download(_recipe, filename="fig3d-provenance.json", label="Download recipe"),
        ], wrap=True),
    ])
    return


if __name__ == "__main__":
    app.run()
