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
        highlighted_regions,
        load_improvement_snapshot, make_improvement_interactive,
        make_improvement_static, prepare_density, prepare_improvement,
        improvement_provenance,
    )

    return (
        highlight_table_rows, highlighted_regions, io, json, load_improvement_snapshot, make_improvement_interactive,
        make_improvement_static, mo, os, plt, prepare_density,
        prepare_improvement, improvement_provenance, make_linked_widget, select_observations, zipfile,
    )


@app.cell
def _(mo, os):
    source_path = mo.ui.text(
        value=os.environ.get("SAMPLEWORKS_RESULTS", "results"),
        label="Results folder", full_width=True,
    ).form(submit_button_label="Load snapshot")
    mo.vstack([
        mo.md("""
        # Figure 3C — RSCC improvement over baseline

        Each point is one altloc region, aggregated across the **four main models**
        at the guidance weights used in the existing analysis. At least three paired
        model rows are required.

        **The existing calculation is a difference of medians:** median guided RSCC
        minus median baseline RSCC. It is not the median of per-model improvements.
        Select a region below to inspect those individual pairs and source records.

        The dashed horizontal line is no change; the diagonal is the improvement
        ceiling at guided RSCC = 1. The lower panel shows the same baseline/guided
        Gaussian density estimates as the existing figure.

        | Figure | Output |
        |---|---|
        | Aggregated RSCC improvement | Plot-toolbar visible SVG; full reconstruction ZIP with SVG and highlighted-region table |

        This reconstructs notebook 02, not a newly verified historical paper snapshot.
        Inputs are read-only; no inference or evaluation is run.
        """),
        source_path,
    ])
    return (source_path,)


@app.cell
def _(
    load_improvement_snapshot, mo, os, prepare_density,
    prepare_improvement, improvement_provenance, source_path,
):
    _path = source_path.value
    if _path is None and not mo.running_in_notebook():
        _path = os.environ.get("SAMPLEWORKS_RESULTS", "results")
    mo.stop(_path is None, mo.md("Choose the results folder and load the snapshot."))
    try:
        _rscc, _manifest = load_improvement_snapshot(_path)
        points, contributions, coverage = prepare_improvement(_rscc)
        density, density_coverage = prepare_density(points)
        provenance = improvement_provenance(_manifest)
        _error = None
    except (OSError, ValueError, KeyError) as _exc:
        _error = str(_exc)
    mo.stop(_error is not None, mo.callout(f"Snapshot could not be loaded: {_error}", kind="warn"))
    return contributions, coverage, density, density_coverage, points, provenance


@app.cell
def _(density, make_improvement_interactive, make_linked_widget, mo, points):
    mo.stop(not points.plotted.any(), mo.callout(
        "No regions meet the model-count and finite-coordinate requirements.", kind="warn"
    ))
    base_figure = make_improvement_interactive(points, density)
    chart_point_ids = [[_row[0] for _row in _trace.customdata] for _trace in base_figure.data]
    get_plot_selection, _set_plot_selection = mo.state([])
    figure_widget, clear_plot_selection = make_linked_widget(base_figure, _set_plot_selection, "fig3c-visible")
    return base_figure, chart_point_ids, clear_plot_selection, figure_widget, get_plot_selection


@app.cell
def _(clear_plot_selection, figure_widget, mo, provenance):
    chart = mo.ui.anywidget(figure_widget)
    clear_button = mo.ui.button(label="Clear plot selection", on_click=lambda _: clear_plot_selection())
    mo.vstack([
        mo.md(f"Loaded snapshot: `{provenance['snapshot_sha256'][:16]}`. "
              "The model population is fixed for this aggregated paper panel. "
              "**The legend toggles point/highlight/density display, not model membership.**"),
        chart,
        clear_button,
        mo.md("Hover for region identity and contributor counts. Box/lasso-select the "
              "upper scatter, then select table rows to inspect model pairs. Density "
              "curve vertices are not individual regions. The camera button saves "
              "the current visible view as SVG. Table-selected regions get black "
              "outlines; clearing rows removes them."),
    ])
    return chart, clear_button


@app.cell
def _(base_figure, chart_point_ids, figure_widget, get_plot_selection, highlight_table_rows, mo, points, select_observations):
    _shown = select_observations(points.loc[points.plotted], chart_point_ids, get_plot_selection())
    highlight_table_rows(figure_widget, base_figure, points.iloc[:0])
    region_table = mo.ui.table(
        _shown[["point_id", "label", "selection", "baseline_value", "guided_value", "improvement",
                "n_model_rows", "n_finite_pairs", "top", "decreased", "far", "curated"]],
        selection="multi", page_size=10, label="Region summaries",
        on_change=lambda rows: highlight_table_rows(figure_widget, base_figure, rows),
    )
    mo.vstack([
        mo.md(f"## Regions\n{len(_shown)} regions shown. Clearing selection restores all "
              "regions; toggling highlight visibility does not recalculate the statistics."),
        region_table,
    ])
    return (region_table,)


@app.cell
def _(contributions, mo, region_table):
    _selected = region_table.value
    _rows = contributions.loc[
        contributions.point_id.isin(_selected.point_id)
    ] if len(_selected) else contributions.iloc[:0]
    mo.vstack([
        mo.md("### Contributing model pairs\nEach row shows a model's baseline/guided RSCC, "
              "selected guidance weight, individual improvement, and original CSV record numbers."),
        mo.ui.table(_rows, selection=None, page_size=10),
    ])
    return


@app.cell
def _(coverage, density_coverage, mo, points, provenance):
    mo.vstack([
        mo.md("## Coverage and recipe\nThe legacy threshold counts paired model rows, "
              "not finite pairs. Finite baseline/guided counts and undrawn regions "
              "remain explicit. Marginal density populations are counted independently."),
        mo.ui.table(coverage, selection=None),
        mo.ui.table(points.loc[~points.plotted], selection=None, label="Undrawn regions"),
        mo.ui.table(density_coverage, selection=None),
        mo.accordion({"Source checksums and reconstruction recipe": mo.json(provenance)}),
    ])
    return


@app.cell
def _(
    contributions, density, highlighted_regions, io, json, make_improvement_static, mo,
    plt, points, provenance, zipfile,
):
    mo.stop(not points.plotted.any(), mo.md("No region reconstruction to export."))
    _figure = make_improvement_static(points, density)
    _svg = io.BytesIO()
    _figure.savefig(_svg, format="svg", bbox_inches="tight")
    plt.close(_figure)
    _highlights = highlighted_regions(points)
    _highlight_csv = _highlights.to_csv(index=False).encode()
    _bundle = io.BytesIO()
    with zipfile.ZipFile(_bundle, "w", compression=zipfile.ZIP_DEFLATED) as _archive:
        _archive.writestr("fig3c-reconstruction.svg", _svg.getvalue())
        _archive.writestr("fig3c-highlighted-regions.csv", _highlight_csv)
        _archive.writestr("fig3c-provenance.json", json.dumps(provenance, indent=2))
    mo.vstack([
        mo.md("## Circled regions and reconstruction downloads\n"
              "The static figure has no fixed segment text. This companion table identifies "
              "every circled region, including overlapping categories and their colors. "
              "The ZIP contains the full SVG, this table, and the reconstruction recipe. "
              "It is not filtered by chart selection or legend visibility."),
        mo.ui.table(_highlights, selection=None, page_size=10, label="Circled regions"),
        mo.hstack([
            mo.download(_bundle.getvalue(), filename="fig3c-reconstruction.zip",
                        label="Download SVG + circled-region table + recipe"),
            mo.download(_svg.getvalue(), filename="fig3c-reconstruction.svg",
                        label="Download full reconstruction SVG"),
            mo.download(_highlight_csv, filename="fig3c-highlighted-regions.csv",
                        label="Download circled-region table"),
            mo.download(points.to_csv(index=False).encode(), filename="fig3c-regions.csv",
                        label="Download region summaries"),
            mo.download(contributions.to_csv(index=False).encode(), filename="fig3c-model-pairs.csv",
                        label="Download model pairs"),
            mo.download(density.to_csv(index=False).encode(), filename="fig3c-density.csv",
                        label="Download density coordinates"),
            mo.download(json.dumps(provenance, indent=2).encode(), filename="fig3c-provenance.json",
                        label="Download recipe and checksums"),
        ], wrap=True),
    ])
    return


if __name__ == "__main__":
    app.run()
