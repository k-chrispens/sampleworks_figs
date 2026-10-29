"""Small Plotly primitives and linked inspection for the supplementary notebooks."""

from dataclasses import dataclass
from io import BytesIO
import json
from uuid import uuid4
from zipfile import ZipFile

import pandas as pd
import plotly.graph_objects as go

from plot_interaction import highlight_table_rows, make_linked_widget, select_observations


def observations(frame):
    """Identify observations by scientific keys, never their displayed row position."""
    out = frame.copy()
    keys = [c for c in (
        "model_name", "protein", "selection", "condition", "guidance_weight",
        "panel", "budget", "budget_value",
    ) if c in out]
    out["point_id"] = [
        json.dumps(row, separators=(",", ":"), default=str)
        for row in out[keys].itertuples(index=False, name=None)
    ]
    if out.point_id.duplicated().any():
        raise ValueError(f"Duplicate plotted observation keys: {keys}")
    return out


def records(frame):
    """JSON records with real nulls rather than non-standard NaN literals."""
    return json.loads(frame.to_json(orient="records"))


def finish(figure, points, summaries=None):
    figure.update_layout(
        template="plotly_white", dragmode="lasso", hovermode="closest",
        margin=dict(l=65, r=30, t=85, b=65),
        meta={"observations": records(points),
              "summaries": records(summaries) if summaries is not None else []},
    )
    return figure


def scatter(points, *, color="#777777", name="", size=6, opacity=0.65, **kwargs):
    """Only this trace represents selectable observations; guides use empty IDs."""
    fields = [c for c in (
        "model_name", "protein", "selection", "condition", "guidance_weight",
        "budget", "budget_value", "n", "sd", "frac_b_cluster", "b_occupancy",
        "classification", "label", "geometry_penalty",
    ) if c in points]
    custom = [[row.point_id, *[str(row[c]) for c in fields]]
              for _, row in points.iterrows()]
    hover = "<br>".join(f"{c}: %{{customdata[{i + 1}]}}" for i, c in enumerate(fields))
    return go.Scatter(
        x=points.x.tolist(), y=points.y.tolist(), customdata=custom,
        mode="markers", name=name, meta={"observations": True},
        marker=dict(color=color, size=size, opacity=opacity, line=dict(width=0)),
        hovertemplate=hover + "<br>x: %{x:.4g}<br>y: %{y:.4g}<extra>%{fullData.name}</extra>",
        **kwargs,
    )


def line(x, y, *, color="#aaaaaa", name="", **kwargs):
    return go.Scatter(
        x=list(x), y=list(y), customdata=[[""]] * len(x),
        mode="lines", line=dict(color=color), name=name,
        hoverinfo="skip", showlegend=False, **kwargs,
    )


def band(figure, x, low, high, color, *, row=1, col=1, legendgroup=None):
    figure.add_trace(line(x, low, color=color, opacity=0, legendgroup=legendgroup), row=row, col=col)
    figure.add_trace(line(x, high, color=color, fill="tonexty", opacity=0.2, legendgroup=legendgroup),
                     row=row, col=col)


@dataclass
class LinkedPlot:
    figure: go.Figure
    widget: go.FigureWidget
    points: pd.DataFrame
    ids: list
    view: object
    chart: object
    clear_button: object
    generation: str
    table: object = None
    has_outlines: bool = False


def outline_rows(panel, rows):
    if panel.has_outlines or not rows.empty:
        highlight_table_rows(panel.widget, panel.figure, rows)
    panel.has_outlines = not rows.empty


def linked_plot(figure, mo, set_selection, key):
    """Keep the callback-bearing chart/button alive independently of table state."""
    generation = uuid4().hex
    figure.layout.meta["inspection_key"] = key
    points = pd.DataFrame(figure.layout.meta["observations"])
    if "point_id" not in points:
        points["point_id"] = pd.Series(dtype=str)
    def publish(event):
        set_selection(lambda old: {
            **old, key: {"generation": generation, "revision": uuid4().hex, "points": event},
        })
    widget, clear = make_linked_widget(figure, publish, key)
    ids = [[row[0] for row in trace.customdata] for trace in figure.data]
    button = mo.ui.button(label="Clear plot selection", on_click=lambda _: clear())
    chart = mo.ui.anywidget(widget)
    view = mo.vstack([chart, button])
    return LinkedPlot(figure, widget, points, ids, view, chart, button, generation)


def inspection_tables(panels, selections, mo):
    """Called from a separate reactive cell, never from the chart construction cell."""
    tables = {}
    for key, panel in panels.items():
        if panel is None:
            continue
        selected = selections.get(key, {})
        current = selected.get("generation") == panel.generation
        event = selected.get("points", []) if current else []
        shown = select_observations(panel.points, panel.ids, event)
        # Marimo retires cell-created UI on rerun: rebuild tables, never the charts.
        outline_rows(panel, shown.iloc[:0])
        panel.table = mo.ui.table(
            shown, selection="multi", page_size=10,
            on_change=lambda rows, p=panel: outline_rows(p, rows),
        )
        table = panel.table
        tables[key] = mo.vstack([
            mo.md("Plot selection filters these rows. Select rows to outline observations; "
                  "clearing the plot restores all rows. Downloads use the uninspected figure."),
            table,
            mo.download(data=panel.points.to_csv(index=False).encode(), filename=f"{key}-points.csv",
                        label="All prepared points (CSV)"),
            mo.download(data=json.dumps(panel.figure.layout.meta).encode(), filename=f"{key}-data.json",
                        label="Prepared points and summaries (JSON)"),
        ])
        if "contributors" in panel.figure.layout.meta:
            tables[key] = mo.vstack([tables[key], mo.ui.table(
                pd.DataFrame(panel.figure.layout.meta["contributors"]).loc[
                    lambda frame: frame.point_id.isin(shown.point_id)],
                selection=None, page_size=10,
            )])
    return mo.accordion(tables)


def reconstruction_bundle(figure):
    """Use the original figure, not the live widget or transient inspection state."""
    output = BytesIO()
    with ZipFile(output, "w") as archive:
        archive.writestr("figure.svg", figure.to_image(format="svg"))
        archive.writestr("points.csv", pd.DataFrame(figure.layout.meta["observations"]).to_csv(index=False))
        points = pd.DataFrame(figure.layout.meta["observations"])
        if "highlight_colors" in points:
            archive.writestr("highlighted-regions.csv",
                             points.loc[points.highlight_colors != ""].to_csv(index=False))
        if "contributors" in figure.layout.meta:
            archive.writestr("contributors.csv",
                             pd.DataFrame(figure.layout.meta["contributors"]).to_csv(index=False))
        archive.writestr("summaries.csv", pd.DataFrame(figure.layout.meta["summaries"]).to_csv(index=False))
        archive.writestr("recipe.json", json.dumps({
            "scope": "Full prepared figure; no transient legend/inspection state",
            "provenance": "Prepared values only; historical upstream provenance remains unresolved",
            "figure": json.loads(figure.to_json()),
        }, indent=2))
    return output.getvalue()
