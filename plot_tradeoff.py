"""Individual-protein guidance trajectories; raw deltas and symlog display values."""

import numpy as np
from plotly.colors import sample_colorscale
from plotly.subplots import make_subplots

from guidance import symlog
from style import MODEL_DISPLAY
from supplementary_plot import finish, line, observations, scatter


def make_spiral_figure(model_df, model_key, *, xlim=None, ylim=None):
    points = observations(model_df.sort_values(["protein", "guidance_weight"]))
    points["x"] = points.delta_rscc
    points["y"] = symlog(points.geometry_penalty)
    fig = make_subplots(rows=1, cols=1)
    maximum = points.guidance_weight.max() if len(points) else 1
    segments = {}
    for _, group in points.groupby("protein", sort=False):
        x = [0, *group.x]
        y = [0, *group.y]
        for i, gw in enumerate(group.guidance_weight, 1):
            xs, ys = segments.setdefault(gw, ([], []))
            xs.extend([x[i-1], x[i], None])
            ys.extend([y[i-1], y[i], None])
    for gw, (xs, ys) in sorted(segments.items()):
        color = sample_colorscale("Plasma", [float(gw / maximum) if maximum > 0 else 0])[0]
        fig.add_trace(line(xs, ys, color=color), row=1, col=1)
    trace = scatter(points, name="Protein / guidance weight", opacity=.85, showlegend=False)
    trace.marker.color = points.guidance_weight.tolist()
    trace.marker.colorscale = "Plasma"
    trace.marker.colorbar = dict(title="Guidance weight")
    trace.marker.showscale = True
    fig.add_trace(trace, row=1, col=1)
    ticks = np.array([-1000, -100, -10, -1, 0, 1, 10, 100, 1000, 10000])
    fig.update_xaxes(title_text="Δ RSCC", range=xlim, zeroline=True)
    fig.update_yaxes(title_text="Δ geometry penalty (symlog)", tickvals=symlog(ticks),
                     ticktext=ticks.astype(str), range=symlog(ylim).tolist() if ylim else None)
    fig.update_layout(title=MODEL_DISPLAY.get(model_key, model_key), height=600)
    return finish(fig, points)
