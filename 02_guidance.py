import marimo

__generated_with = "0.23.2"
app = marimo.App()


@app.cell
def _():
    import io
    import json
    import os

    import marimo as mo
    import matplotlib.pyplot as plt

    from guidance import (
        guidance_provenance, load_guidance_snapshot, make_guidance_interactive,
        make_guidance_static, prepare_guidance,
    )
    from plot_interaction import highlight_table_rows, make_linked_widget, select_observations
    from style import MODELS_ALT

    return (
        MODELS_ALT, guidance_provenance, io, json, load_guidance_snapshot,
        make_guidance_interactive, make_guidance_static, mo, os, plt,
        prepare_guidance, select_observations, highlight_table_rows, make_linked_widget,
    )


@app.cell
def _(mo, os):
    source_path = mo.ui.text(
        value=os.environ.get("SAMPLEWORKS_RESULTS", "results"),
        label="Results folder", full_width=True,
    ).form(submit_button_label="Load snapshot")
    mo.vstack([
        mo.md("""
        # Figure 3B — guidance versus geometry

        Each point summarizes protein-level changes at a guidance weight:
        **median Δ RSCC versus median geometry penalty, with IQR error bars**.
        Larger markers mean stronger guidance within each model; black outlines
        mark the weights selected in the existing analysis.

        This reconstructs notebook 03, not a newly verified historical paper
        snapshot. It preserves protein-mean RSCC, member-mean geometry, and
        NaN-skipping marginal summaries. Counts below expose incomplete coverage.
        The vertical axis uses the existing symmetric-log scale; hover/tables
        report the original, untransformed penalties.

        | Figure | Output |
        |---|---|
        | Guidance/geometry tradeoff | Plot-toolbar visible SVG; `fig3b-reconstruction.svg` for main models |

        Inputs are read-only. No evaluation or model inference is run here.
        """),
        source_path,
    ])
    return (source_path,)


@app.cell
def _(guidance_provenance, load_guidance_snapshot, mo, os, prepare_guidance, source_path):
    _path = source_path.value
    if _path is None and not mo.running_in_notebook():
        _path = os.environ.get("SAMPLEWORKS_RESULTS", "results")
    mo.stop(_path is None, mo.md("Choose the results folder and load the snapshot."))
    try:
        _frames, _manifest = load_guidance_snapshot(_path)
        proteins, summary, coverage = prepare_guidance(_frames)
        provenance = guidance_provenance(_manifest)
        _error = None
    except (OSError, ValueError, KeyError) as _exc:
        _error = str(_exc)
    mo.stop(_error is not None, mo.callout(f"Snapshot could not be loaded: {_error}", kind="warn"))
    return coverage, proteins, provenance, summary


@app.cell
def _(make_guidance_interactive, make_linked_widget, mo, summary):
    mo.stop(not summary.plotted.any(), mo.callout(
        "No guidance points have sufficient finite coverage.", kind="warn"
    ))
    base_figure = make_guidance_interactive(summary)
    chart_point_ids = [[_row[0] for _row in _trace.customdata] for _trace in base_figure.data]
    get_plot_selection, _set_plot_selection = mo.state([])
    figure_widget, clear_plot_selection = make_linked_widget(base_figure, _set_plot_selection, "fig3b-visible")
    return base_figure, chart_point_ids, clear_plot_selection, figure_widget, get_plot_selection


@app.cell
def _(clear_plot_selection, figure_widget, mo, provenance):
    chart = mo.ui.anywidget(figure_widget)
    clear_button = mo.ui.button(label="Clear plot selection", on_click=lambda _: clear_plot_selection())
    mo.vstack([
        mo.md(f"Loaded snapshot: `{provenance['snapshot_sha256'][:16]}`. "
              "**Click legend entries to show/hide models; double-click to isolate one.** "
              "Variants start hidden. Use the camera button to download the visible SVG."),
        chart,
        clear_button,
        mo.md("Box/lasso-select points, then select table rows to inspect contributing proteins. "
              "Table-selected points get black outlines; clearing rows removes them. "
              "The origin is the baseline reference, not an additional observation. "
              "Marker sizes are scaled within each model, not across models."),
    ])
    return chart, clear_button


@app.cell
def _(base_figure, chart_point_ids, figure_widget, get_plot_selection, highlight_table_rows, mo, select_observations, summary):
    _shown = select_observations(summary.loc[summary.plotted], chart_point_ids, get_plot_selection())
    highlight_table_rows(figure_widget, base_figure, summary.iloc[:0])
    point_table = mo.ui.table(
        _shown[["point_id", "model_name", "guidance_weight", "x_med", "x_q1", "x_q3",
                "y_med", "y_q1", "y_q3", "n", "n_rscc", "n_geometry", "incomplete_geometry"]],
        selection="multi", page_size=10, label="Guidance points",
        on_change=lambda rows: highlight_table_rows(figure_widget, base_figure, rows),
    )
    mo.vstack([
        mo.md(f"## Guidance points\n{len(_shown)} points shown. Without a chart selection, "
              "this table includes all loaded models, including legend-hidden ones."),
        point_table,
    ])
    return (point_table,)


@app.cell
def _(mo, point_table, proteins):
    _selected = point_table.value
    _rows = proteins.loc[proteins.point_id.isin(_selected.point_id)] if len(_selected) else proteins.iloc[:0]
    mo.vstack([
        mo.md("### Contributing proteins\nBaseline and guided metric means, their differences, "
              "and the resulting geometry penalty. NaNs are not filled from another source."),
        mo.ui.table(_rows, selection=None, page_size=10),
    ])
    return


@app.cell
def _(coverage, mo, provenance, summary):
    mo.vstack([
        mo.md("## Coverage and provenance\nPairing counts are shown separately for each metric. "
              "The legacy minimum is three protein rows; finite RSCC and geometry counts may "
              "differ. Groups with non-finite plotted coordinates/IQRs are not drawn."),
        mo.ui.table(coverage, selection=None, page_size=10),
        mo.ui.table(summary.loc[~summary.plotted], selection=None, label="Undrawn groups"),
        mo.accordion({"Source checksums, member joins, and recipe": mo.json(provenance)}),
    ])
    return


@app.cell
def _(MODELS_ALT, io, json, make_guidance_static, mo, plt, proteins, provenance, summary):
    _main = summary.loc[~summary.model_name.isin(MODELS_ALT)]
    mo.stop(not _main.plotted.any(), mo.md("No main-model reconstruction to export."))
    _figure = make_guidance_static(_main)
    _svg = io.BytesIO()
    _figure.savefig(_svg, format="svg", bbox_inches="tight")
    plt.close(_figure)
    mo.hstack([
        mo.download(_svg.getvalue(), filename="fig3b-reconstruction.svg",
                    label="Download main-model reconstruction SVG"),
        mo.download(summary.to_csv(index=False).encode(), filename="fig3b-points.csv",
                    label="Download all-model summaries"),
        mo.download(proteins.to_csv(index=False).encode(), filename="fig3b-proteins.csv",
                    label="Download protein contributions"),
        mo.download(json.dumps(provenance, indent=2).encode(), filename="fig3b-provenance.json",
                    label="Download recipe and checksums"),
    ], wrap=True)
    return


if __name__ == "__main__":
    app.run()
