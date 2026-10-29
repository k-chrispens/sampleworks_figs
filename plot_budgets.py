"""Sampling-budget summaries and their contributing selection measurements."""

import numpy as np
import pandas as pd
from plotly.subplots import make_subplots

from style import MODEL_DISPLAY, MODEL_PALETTE
from supplementary_plot import band, finish, observations, records, scatter


def make_budget_figure(frame, model_order, *, default_recycles=None, default_steps=200):
    ensemble = default_recycles is None
    panels = [("ensemble_size", frame)] if ensemble else [
        ("n_recycles", frame.loc[frame.n_diffusion_steps == default_steps]),
        ("n_diffusion_steps", frame.loc[frame.n_recycles == frame.model_name.map(default_recycles)]),
    ]
    summaries, contributors = [], []
    for axis, subset in panels:
        summary = subset.groupby(["model_name", axis]).rscc.agg(mean="mean", sd="std", n="count").reset_index()
        summary = summary.rename(columns={axis: "budget_value"}).assign(budget=axis)
        summary = observations(summary)
        subset = subset.assign(budget=axis, budget_value=subset[axis])
        contributors.append(subset.merge(summary[["model_name", "budget", "budget_value", "point_id"]],
                                         on=["model_name", "budget", "budget_value"]))
        summaries.append(summary)
    points = pd.concat(summaries, ignore_index=True)
    points["x"], points["y"] = np.nan, points["mean"]
    fig = make_subplots(rows=1, cols=len(panels), subplot_titles=[a for a, _ in panels])
    for col, (axis, _) in enumerate(panels, 1):
        values = sorted(points.loc[points.budget == axis, "budget_value"].unique())
        positions = {v: np.log2(v) if ensemble else i for i, v in enumerate(values)}
        for i, model in enumerate(model_order):
            sub = points.loc[(points.model_name == model) & (points.budget == axis)].sort_values("budget_value").copy()
            dodge = (i - (len(model_order) - 1) / 2) * .06 if ensemble else 0
            sub["x"] = sub.budget_value.map(positions) + dodge
            points.loc[sub.index, "x"] = sub.x
            trace = scatter(sub, color=MODEL_PALETTE[model], name=MODEL_DISPLAY[model],
                            opacity=1, legendgroup=model, showlegend=col == 1)
            trace.mode = "lines+markers"
            trace.line.color = MODEL_PALETTE[model]
            if ensemble:
                trace.error_y = dict(type="data", array=sub.sd.tolist(), visible=True)
            elif len(sub) > 1:
                band(fig, sub.x, sub.y - sub.sd, sub.y + sub.sd, MODEL_PALETTE[model],
                     col=col, legendgroup=model)
            fig.add_trace(trace, row=1, col=col)
        fig.update_xaxes(tickvals=list(positions.values()), ticktext=[str(v) for v in values],
                         title_text=axis, row=1, col=col)
    fig.update_yaxes(title_text="Mean RSCC ± sample SD")
    fig.update_layout(title="6B8X · ensemble size" if ensemble else "5MHX · recycles and diffusion steps",
                      height=550)
    finish(fig, points)
    fig.layout.meta["contributors"] = records(pd.concat(contributors, ignore_index=True))
    return fig
