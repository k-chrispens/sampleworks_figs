"""Link plot/table inspection using stable observation IDs, not displayed row positions."""


def make_linked_widget(figure, on_select, filename):
    """Keep one live Plotly widget so marker updates preserve zoom and legend state."""
    import plotly.graph_objects as go

    widget = go.FigureWidget(figure)
    widget._config = {
        **widget._config,
        "toImageButtonOptions": {"format": "svg", "filename": filename},
    }
    selected_by_trace = {}

    def selected(trace, points, selector):
        index = points.trace_index
        selected_by_trace[index] = [
            {"curveNumber": index, "pointIndex": i}
            for i in points.point_inds
        ]
        # FigureWidget dispatches every trace in order for one selection event.
        # Publish once, so the table never shows an intermediate partial selection.
        if index == len(widget.data) - 1:
            on_select([point for group in selected_by_trace.values() for point in group])

    def clear():
        selected_by_trace.clear()
        with widget.batch_update():
            for trace in widget.data:
                trace.selectedpoints = None
            widget.layout.selections = ()
        on_select([])

    def deselected(trace, points):
        if points.trace_index == len(widget.data) - 1:
            clear()

    for trace in widget.data:
        trace.on_selection(selected)
        trace.on_deselect(deselected)
    return widget, clear


def highlight_table_rows(widget, figure, rows):
    """Change marker styles only; never send the original layout/visibility back."""
    highlighted = highlight_observations(figure, rows.point_id)
    with widget.batch_update():
        for target, source in zip(widget.data, highlighted.data):
            if source.meta and source.meta.get("observations"):
                target.marker = source.marker


def highlight_observations(figure, point_ids):
    """Copy the display figure and outline selected observations without changing data."""
    import plotly.graph_objects as go

    result = go.Figure(figure)
    wanted = set(point_ids)
    if not wanted:
        return result
    for trace in result.data:
        if not trace.meta or not trace.meta.get("observations"):
            continue
        selected = [bool(row[0]) and row[0] in wanted for row in trace.customdata]
        if not any(selected):
            continue

        def values(value, default):
            if value is None:
                return [default] * len(selected)
            if isinstance(value, (str, int, float)):
                return [value] * len(selected)
            return list(value)

        sizes = values(trace.marker.size, 6)
        opacity = values(trace.marker.opacity, 1)
        widths = values(trace.marker.line.width, 0)
        colors = values(trace.marker.line.color, "#000000")
        trace.marker.size = [max(size, 15) if picked else size for size, picked in zip(sizes, selected)]
        trace.marker.opacity = [1 if picked else alpha for alpha, picked in zip(opacity, selected)]
        trace.marker.line.width = [2.5 if picked else width for width, picked in zip(widths, selected)]
        trace.marker.line.color = ["#000000" if picked else color for color, picked in zip(colors, selected)]
    return result


def select_observations(points, chart_point_ids, selection):
    """Use rendered vertex IDs, never positions in a sorted/filtered table."""
    if not selection:
        return points
    ids = set()
    for point in selection:
        if point.get("customdata"):
            ids.add(point["customdata"][0])
        else:
            # marimo's native event payload may omit Plotly customdata.
            index = point.get("pointIndex", point.get("pointNumber"))
            if index is not None and "curveNumber" in point:
                ids.add(chart_point_ids[point["curveNumber"]][index])
    return points.loc[points.point_id.isin(ids)]
