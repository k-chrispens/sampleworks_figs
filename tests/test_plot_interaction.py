import unittest

import pandas as pd
import plotly.graph_objects as go
from plotly.callbacks import BoxSelector, Points

from plot_interaction import (
    highlight_observations, highlight_table_rows, make_linked_widget, select_observations,
)


def figure():
    result = go.Figure()
    for axis in ("x", "x2"):
        result.add_trace(go.Scatter(
            x=[0, 1, 2], y=[0, 2, 1], xaxis=axis,
            customdata=[[""], ["a"], ["b"]], meta={"observations": True},
            marker={"size": [0, 5, 7], "opacity": [0, .3, .4],
                    "line": {"width": [0, 1, 2], "color": ["white", "red", "blue"]}},
        ))
    result.add_trace(go.Scatter(
        x=[1], y=[2], customdata=[["a"]], meta={"observations": False},
        marker={"symbol": "circle-open", "color": "green", "size": 11},
    ))
    result.add_trace(go.Scatter(x=[1], y=[.2], mode="lines", customdata=[[""]]))
    return result


class PlotInteractionTest(unittest.TestCase):
    def test_outlines_use_ids_on_both_panels_without_changing_scientific_rings(self):
        original = figure()
        before = original.to_json()
        highlighted = highlight_observations(original, ["b"])
        for trace in highlighted.data[:2]:
            self.assertEqual(list(trace.marker.size), [0, 5, 15])
            self.assertEqual(list(trace.marker.opacity), [0, .3, 1])
            self.assertEqual(list(trace.marker.line.width), [0, 1, 2.5])
            self.assertEqual(list(trace.marker.line.color), ["white", "red", "#000000"])
        for i in (2, 3):
            self.assertEqual(highlighted.data[i].to_plotly_json(), original.data[i].to_plotly_json())
        self.assertEqual(original.to_json(), before)
        for left, right in zip(original.data, highlighted.data):
            self.assertEqual(list(left.x), list(right.x))
            self.assertEqual(list(left.y), list(right.y))
            self.assertEqual(list(left.customdata), list(right.customdata))

    def test_no_selection_unknown_id_and_anchor_never_create_outlines(self):
        original = figure()
        for ids in ([], ["missing"], [""]):
            self.assertEqual(highlight_observations(original, ids).to_plotly_json(), original.to_plotly_json())

    def test_in_place_update_and_clear_preserve_view_and_restore_original_markers(self):
        original = figure()
        widget, _ = make_linked_widget(original, lambda value: None, "visible")
        widget.layout.xaxis.range = [.5, 1.5]
        widget.data[0].visible = "legendonly"
        # Sorted table order and index labels must not be used as plot vertex positions.
        rows = pd.DataFrame({"point_id": ["b", "a"]}, index=[99, 20]).iloc[:1]
        highlight_table_rows(widget, original, rows)
        self.assertEqual(widget.data[0].marker.line.color[2], "#000000")
        highlight_table_rows(widget, original, rows.iloc[:0])
        self.assertEqual(widget.data[0].marker.to_plotly_json(), original.data[0].marker.to_plotly_json())
        self.assertEqual(widget.layout.xaxis.range, (.5, 1.5))
        self.assertEqual(widget.data[0].visible, "legendonly")
        self.assertEqual(widget._config["toImageButtonOptions"]["format"], "svg")

    def test_plot_callbacks_combine_traces_and_deduplicate_shared_observations(self):
        updates = []
        widget, clear = make_linked_widget(figure(), updates.append, "visible")
        for index, trace in enumerate(widget.data):
            indices = {0: [1], 2: [0]}.get(index, [])
            trace._dispatch_on_selection(Points(point_inds=indices, trace_index=index), BoxSelector())
        self.assertEqual(len(updates), 1)
        points = pd.DataFrame({"point_id": ["a", "b"]})
        ids = [[row[0] for row in trace.customdata] for trace in widget.data]
        self.assertEqual(select_observations(points, ids, updates[-1]).point_id.tolist(), ["a"])
        for index, trace in enumerate(widget.data):
            trace._dispatch_on_selection(Points(point_inds=[], trace_index=index), BoxSelector())
        self.assertEqual(len(updates), 2)
        self.assertEqual(updates[-1], [])
        for index, trace in enumerate(widget.data):
            trace._dispatch_on_deselect(Points(trace_index=index))
        self.assertEqual(len(updates), 3)
        self.assertEqual(updates[-1], [])
        widget.data[0].selectedpoints = [1]
        widget.layout.selections = [{"x0": 0, "x1": 2, "y0": 0, "y1": 2}]
        clear()
        self.assertTrue(all(trace.selectedpoints is None for trace in widget.data))
        self.assertEqual(tuple(widget.layout.selections), ())

    def test_density_only_selection_does_not_invent_region_observations(self):
        updates = []
        widget, _ = make_linked_widget(figure(), updates.append, "visible")
        for index, trace in enumerate(widget.data):
            trace._dispatch_on_selection(
                Points(point_inds=[0] if index == 3 else [], trace_index=index), BoxSelector()
            )
        ids = [[row[0] for row in trace.customdata] for trace in widget.data]
        points = pd.DataFrame({"point_id": ["a", "b"]})
        self.assertTrue(select_observations(points, ids, updates[-1]).empty)


if __name__ == "__main__":
    unittest.main()
