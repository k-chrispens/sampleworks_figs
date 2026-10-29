"""Conformational-class Plotly grid from final baseline/guided pairs."""

import numpy as np
from plotly.subplots import make_subplots

from style import CLASSIFICATION_DISPLAY, CLASSIFICATION_ORDER, CLASSIFICATION_PALETTE, MODEL_DISPLAY
from supplementary_plot import finish, line, observations, scatter


def make_by_class_combined_grid(df, *, model_order=None):
    models = model_order if model_order is not None else [
        m for m in ("boltz2X-RAY", "boltz2MD", "protenix", "rf3") if m in set(df.model_name)]
    classes = [c for c in CLASSIFICATION_ORDER if c in set(df.classification)]
    points = observations(df)
    points["x"] = points.rscc_baseline
    points["y"] = points.rscc_optimal - points.rscc_baseline
    counts = points.drop_duplicates(["protein", "selection"]).classification.value_counts()
    fig = make_subplots(rows=max(1, len(models)), cols=max(1, len(classes)),
                        row_titles=[MODEL_DISPLAY.get(m, m) for m in models],
                        column_titles=[f"{CLASSIFICATION_DISPLAY.get(c, c)}<br>n={counts[c]}"
                                       for c in classes])
    values = points.y[np.isfinite(points.y)]
    low = min(0, values.min()) - .02 if len(values) else -.05
    for row, model in enumerate(models, 1):
        background = points.loc[points.model_name == model]
        for col, cls in enumerate(classes, 1):
            shadow = scatter(background, opacity=.08, name="All classes", showlegend=False)
            shadow.meta = {"observations": False}
            shadow.customdata = [[""]] * len(background)
            shadow.hoverinfo, shadow.hovertemplate = "skip", None
            fig.add_trace(shadow, row=row, col=col)
            subset = background.loc[background.classification == cls]
            fig.add_trace(scatter(subset, color=CLASSIFICATION_PALETTE[cls], showlegend=False),
                          row=row, col=col)
            fig.add_trace(line([0, 1], [0, 0]), row=row, col=col)
            fig.add_trace(line([0, 1], [1, 0]), row=row, col=col)
    fig.update_xaxes(range=[0, 1], title_text="Baseline RSCC")
    fig.update_yaxes(range=[low, 1], title_text="Δ RSCC")
    fig.update_layout(height=max(350, 260 * len(models)), title="RSCC by conformational class")
    fig.update_annotations(font_size=12)
    finish(fig, points)
    fig.update_layout(margin=dict(r=105))
    return fig
