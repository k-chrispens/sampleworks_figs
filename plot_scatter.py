"""Per-model Plotly residual, cross-metric and minimum-RMSD figures."""

import numpy as np
import pandas as pd
from plotly.subplots import make_subplots
from scipy.stats import gaussian_kde, pearsonr

from style import MODEL_DISPLAY
from supplementary_plot import finish, line, observations, scatter

RSCC_HIGHLIGHT_PROTEINS = [
    "6B8X 181-183", "5MHX 162-164",
    "5SOP 131-133", "6DUR 389-391", "2QCV 2-4", "6RP1 5",
    "6NI5 151-153", "6NI5 201-202", "6NI5 204-206",
    "6NI5 209-210", "6NI5 213-215",
    "6NI6 151-153", "6NI6 201-202", "6NI6 204-206",
    "6NI6 209-210", "6NI6 213-215",
]


def rscc_residual_label_override(df, *, x_col="rscc_baseline", y_col="rscc_optimal",
                                 label_col="label", label_top_n=3, label_far_from_ceiling_n=5):
    delta = df[y_col] - df[x_col]
    top = df.assign(_d=delta).nlargest(label_top_n, "_d")
    far = df.assign(_b=1.0 - df[y_col]).nlargest(label_far_from_ceiling_n, "_b")
    return {"top": top.loc[top._d > 0, label_col].tolist(),
            "decreased": df.loc[delta < 0, label_col].tolist(),
            "far": far.loc[far._b > 0, label_col].tolist()}


def _rings(fig, points, groups, *, row=1, col=1):
    for color, labels in groups.items():
        subset = points.loc[points.label.isin(labels)]
        trace = scatter(subset, color=color, name=color, size=13, opacity=1, showlegend=False)
        trace.marker.symbol = "circle-open"
        trace.marker.line.width = 2
        trace.meta = {"observations": False}
        fig.add_trace(trace, row=row, col=col)


def make_scatter_comparison(df, model_name, x_col, y_col, x_label, y_label,
                            title=None, *, draw_regression=False, show_correlation=False,
                            equal_aspect=True, label_col="protein", highlight_proteins=None,
                            residual=False, residual_label_override=None,
                            marginal_kde=False):
    points = observations(df)
    points["x"] = points[x_col]
    points["y"] = points[y_col] - points[x_col] if residual else points[y_col]
    points["label"] = points[label_col].astype(str)
    points["highlight_colors"] = ""
    fig = make_subplots(rows=1, cols=2 if marginal_kde else 1,
                        column_widths=[0.82, 0.18] if marginal_kde else None,
                        horizontal_spacing=0.04)
    fig.add_trace(scatter(points, name=MODEL_DISPLAY.get(model_name, model_name)), row=1, col=1)
    groups = {"#1565C0": highlight_proteins or []}
    if residual:
        selected = residual_label_override or rscc_residual_label_override(
            df, x_col=x_col, y_col=y_col, label_col=label_col)
        top = selected["top"]
        decreased = [v for v in selected["decreased"] if v not in top]
        far = [v for v in selected["far"] if v not in set(top + decreased)]
        groups = {"#2ca02c": top, "#d62728": decreased,
                  "#1565C0": highlight_proteins or [], "#ff7f0e": far}
        fig.add_trace(line([0, 1], [0, 0]), row=1, col=1)
        fig.add_trace(line([0, 1], [1, 0]), row=1, col=1)
        finite = points.y[np.isfinite(points.y)]
        fig.update_xaxes(range=[0, 1], row=1, col=1)
        fig.update_yaxes(range=[min(0, finite.min()) - .02 if len(finite) else -.05, 1],
                         row=1, col=1)
    _rings(fig, points, groups)
    for color, labels in groups.items():
        points.loc[points.label.isin(labels), "highlight_colors"] += color + " "
    summaries = []
    finite = points.loc[np.isfinite(points.x) & np.isfinite(points.y)]
    if draw_regression and len(finite) >= 2 and finite.x.std() > 0:
        slope, intercept = np.polyfit(finite.x, finite.y, 1)
        support = np.linspace(finite.x.min(), finite.x.max(), 100)
        fitted = support * slope + intercept
        fig.add_trace(line(support, fitted, color="#333333"), row=1, col=1)
        summaries.extend(dict(kind="regression", x=x, y=y) for x, y in zip(support, fitted))
    if show_correlation and len(finite) >= 5:
        r = float(pearsonr(finite.x, finite.y).statistic)
        fig.add_annotation(x=.03, y=.97, xref="paper", yref="paper", text=f"Pearson r = {r:.3f}",
                           showarrow=False)
        summaries.append(dict(kind="correlation", r=r, n=len(finite)))
    if marginal_kde:
        for column, color in ((x_col, "#777777"), (y_col, "#1565C0")):
            values = df[column].to_numpy()
            values = values[np.isfinite(values)]
            if len(values) < 2 or np.std(values) == 0:
                continue
            kde = gaussian_kde(values)
            kde.set_bandwidth(kde.factor * .5)
            bw = np.sqrt(kde.covariance.squeeze())
            support = np.linspace(values.min() - 3 * bw, values.max() + 3 * bw, 200)
            density = kde(support)
            fig.add_trace(line(density, support, color=color), row=1, col=2)
            summaries.extend(dict(kind="density", metric=column, x=x, y=y)
                             for x, y in zip(density, support))
        fig.update_yaxes(range=[0, 1], row=1, col=2)
    fig.update_xaxes(title_text=x_label, row=1, col=1)
    fig.update_yaxes(title_text=y_label, row=1, col=1)
    if equal_aspect and not residual:
        fig.update_yaxes(scaleanchor="x", scaleratio=1, row=1, col=1)
    fig.update_layout(title=title or MODEL_DISPLAY.get(model_name, model_name), height=550,
                      showlegend=False)
    return finish(fig, points, pd.DataFrame(summaries))


def make_rscc_scatter(df, model_name, optimal_gw):
    return make_scatter_comparison(
        df, model_name, "rscc_baseline", "rscc_optimal", "Baseline RSCC", "Δ RSCC",
        title=f"{MODEL_DISPLAY.get(model_name, model_name)} · guidance {optimal_gw:g}",
        label_col="label", residual=True, marginal_kde=True,
        highlight_proteins=RSCC_HIGHLIGHT_PROTEINS)


def make_clustering_scatter(df, model_name, optimal_gw, *, residual_label_override=None):
    return make_scatter_comparison(
        df, model_name, "clustering_score_baseline", "clustering_score_optimal",
        "Baseline clustering score", "Δ clustering score",
        title=f"{MODEL_DISPLAY.get(model_name, model_name)} · guidance {optimal_gw:g}",
        label_col="label", residual=True, marginal_kde=True,
        residual_label_override=residual_label_override,
        highlight_proteins=RSCC_HIGHLIGHT_PROTEINS)


def make_min_rmsd_scatter(df, *, title, x_col="min_rmsd_to_A", y_col="min_rmsd_to_B",
                          label_col="label", pinned_highlights=("6B8X 181-183", "5MHX 162-164"),
                          pinned_color="#1565C0", extra_highlight_groups=None, zoom_xylim=(0, 3)):
    points = observations(df)
    points["x"], points["y"] = points[x_col], points[y_col]
    points["label"] = points[label_col]
    points["highlight_colors"] = ""
    fig = make_subplots(rows=1, cols=2, subplot_titles=["Full range", "Zoom"])
    finite = points[["x", "y"]].to_numpy()
    finite = finite[np.isfinite(finite)]
    limits = [(0, max(.1, finite.max() * 1.05) if len(finite) else .1), zoom_xylim]
    extras = extra_highlight_groups or {}
    for col, (panel, bounds) in enumerate(zip(("autofit", "zoom"), limits), 1):
        fig.add_trace(scatter(points, name="Regions", showlegend=False), row=1, col=col)
        fig.add_trace(line(bounds, bounds), row=1, col=col)
        groups = {pinned_color: list(pinned_highlights)}
        panel_extras = extras.get(panel, {}) if set(extras) <= {"autofit", "zoom"} else extras
        for color, labels in panel_extras.items():
            groups.setdefault(color, []).extend(labels)
        visible = points.loc[points.x.between(*bounds) & points.y.between(*bounds)]
        _rings(fig, visible, groups, col=col)
        for color, labels in groups.items():
            points.loc[points.point_id.isin(visible.loc[visible.label.isin(labels), "point_id"]),
                       "highlight_colors"] += f"{panel}:{color} "
        fig.update_xaxes(range=bounds, title_text="Minimum RMSD to A (Å)", row=1, col=col)
        fig.update_yaxes(range=bounds, title_text="Minimum RMSD to B (Å)",
                         scaleanchor="x" if col == 1 else "x2", scaleratio=1, row=1, col=col)
    fig.update_layout(title=title, height=560)
    return finish(fig, points)
