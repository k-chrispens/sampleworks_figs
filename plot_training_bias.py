"""Training-occupancy Plotly figures; summaries are computed before rendering."""

import numpy as np
import pandas as pd
from plotly.subplots import make_subplots
from scipy.stats import gaussian_kde

from style import MODEL_DISPLAY, MODEL_PALETTE
from supplementary_plot import band, finish, line, observations, scatter

_MAIN = ["boltz2X-RAY", "boltz2MD", "protenix", "rf3"]
_CONDITIONS = ["Baseline", "Optimal guidance"]


def _training_bias_prepare_scatter_frame(df, condition_label):
    out = df.copy()
    s = (out["n_assigned_to_A"] + out["n_assigned_to_B"]).astype(float)
    out["s_cluster"] = s
    out["frac_b_cluster"] = np.where(s > 0, out.n_assigned_to_B / s.replace(0, np.nan), np.nan)
    out["condition"] = condition_label
    out = out.dropna(subset=["frac_b_cluster", "b_occupancy"])
    return out[out.s_cluster > 0]


def _training_bias_compute_alpha(s_values, s_max):
    if s_max <= 0:
        return np.zeros_like(s_values, dtype=float)
    return np.clip(.5 * np.log(s_values + 1) / np.log(s_max + 1), 0, 1)


def _training_bias_binned_means(x, y, *, n_bins=10, min_count=5):
    if len(x) == 0:
        return np.array([]), np.array([]), np.array([])
    edges = np.unique(np.quantile(x, np.linspace(0, 1, n_bins + 1)))
    centers, means, stds = [np.full(max(0, len(edges) - 1), np.nan) for _ in range(3)]
    for i in range(len(edges) - 1):
        mask = (x >= edges[i]) & ((x <= edges[i + 1]) if i == len(edges) - 2 else (x < edges[i + 1]))
        if int(mask.sum()) >= min_count:
            centers[i], means[i], stds[i] = x[mask].mean(), y[mask].mean(), np.std(y[mask], ddof=1)
    return centers, means, stds


def _training_bias_slope_r(x, y):
    if len(x) < 2 or float(np.std(x)) == 0:
        return None
    return float(np.polyfit(x, y, 1)[0]), float(np.corrcoef(x, y)[0, 1])


def _color(model, condition):
    return "#777777" if condition == "Baseline" else MODEL_PALETTE[model]


def _calibration(baseline, optimal, *, model_order=None, jitter=.01, seed=0,
                 separate=False, plain=False, max_cluster_size=None, cutoff_status=None,
                 n_bins=10, min_bin_count=5):
    models = _MAIN if model_order is None else model_order
    prepared = pd.concat([_training_bias_prepare_scatter_frame(frame, condition)
                          for frame, condition in zip((baseline, optimal), _CONDITIONS)],
                         ignore_index=True)
    if max_cluster_size is not None:
        prepared = prepared.loc[prepared.n_cluster_evaluated <= max_cluster_size]
    if cutoff_status is not None:
        prepared = prepared.loc[prepared.cutoff_status == cutoff_status]
    points = observations(prepared)
    points["x"], points["y"] = points.frac_b_cluster, points.b_occupancy
    rows = 2 if separate else 1
    cols = max(1, len(models)) * (1 if plain else 2)
    titles = []
    for condition in (_CONDITIONS if separate else [""]):
        for model in models:
            titles.extend([f"{MODEL_DISPLAY.get(model, model)}<br>{condition}"] + ([] if plain else [""]))
    fig = make_subplots(rows=rows, cols=cols, subplot_titles=titles,
                        column_widths=None if plain else [4, 1] * max(1, len(models)),
                        horizontal_spacing=.025, vertical_spacing=.15)
    rng = np.random.default_rng(seed)
    s_max = float(points.s_cluster.max()) if len(points) else 1
    summaries = []
    # Preserve the old display-jitter iteration order in each layout.
    panels = [(r, c, model, [condition]) for r, condition in enumerate(_CONDITIONS, 1)
              for c, model in enumerate(models, 1)] if separate else [
                  (1, c, model, _CONDITIONS) for c, model in enumerate(models, 1)]
    for row, index, model, conditions in panels:
        col = index if plain else index * 2 - 1
        fig.add_trace(line([-.02, 1.02], [.5, .5]), row=row, col=col)
        if not plain:
            fig.add_trace(line([-.02, 1.02], [-.02, 1.02]), row=row, col=col)
        for condition in conditions:
            sub = points.loc[(points.model_name == model) & (points.condition == condition)].copy()
            sub["x"] += rng.uniform(-jitter, jitter, len(sub))
            points.loc[sub.index, "x"] = sub.x
            color = "#333333" if plain else _color(model, condition)
            legendgroup = f"{model}:{condition}"
            alpha = _training_bias_compute_alpha(sub.s_cluster.to_numpy(), s_max)
            if plain:
                alpha *= 2
            fig.add_trace(scatter(sub, color=color, name=f"{MODEL_DISPLAY.get(model, model)} · {condition}",
                                  opacity=alpha.tolist(), legendgroup=legendgroup,
                                  showlegend=True), row=row, col=col)
            if plain:
                continue
            x, y = sub.frac_b_cluster.to_numpy(), sub.b_occupancy.to_numpy()
            centers, means, stds = _training_bias_binned_means(x, y, n_bins=n_bins,
                                                              min_count=min_bin_count)
            band(fig, centers, means - stds, means + stds, color, row=row, col=col, legendgroup=legendgroup)
            fig.add_trace(line(centers, means, color=color, legendgroup=legendgroup), row=row, col=col)
            summaries.extend(dict(kind="bin", model_name=model, condition=condition,
                                  x=a, y=b, sd=s) for a, b, s in zip(centers, means, stds))
            stats = _training_bias_slope_r(x, y)
            if stats is not None:
                slope, r = stats
                summaries.append(dict(kind="regression", model_name=model, condition=condition,
                                      slope=slope, r=r, n=len(x)))
                fig.add_annotation(x=.05, y=.98 if condition == "Baseline" else .87,
                                   xref=f"x{(row-1)*cols+col} domain" if (row-1)*cols+col > 1 else "x domain",
                                   yref=f"y{(row-1)*cols+col} domain" if (row-1)*cols+col > 1 else "y domain",
                                   text=f"{condition}: slope {slope:.2f}, r {r:.2f}",
                                   showarrow=False, font=dict(size=10, color=color))
            counts, edges = np.histogram(y, bins=np.linspace(0, 1, 21))
            centers_y = (edges[:-1] + edges[1:]) / 2
            # Precomputed counts, not a browser-side histogram aggregation.
            fig.add_trace(line(counts, centers_y, color=color, legendgroup=legendgroup), row=row, col=col+1)
            summaries.extend(dict(kind="histogram", model_name=model, condition=condition,
                                  x=int(n), y=v) for n, v in zip(counts, centers_y))
        fig.update_xaxes(range=[-.02, 1.02], title_text="Training state-B fraction", row=row, col=col)
        fig.update_yaxes(range=[-.02, 1.02], title_text="Predicted state-B fraction", row=row, col=col)
        if not plain:
            fig.update_yaxes(range=[-.02, 1.02], row=row, col=col+1)
            fig.update_xaxes(title_text="Count", row=row, col=col+1)
    fig.update_layout(height=760 if separate else 480, legend=dict(orientation="h", y=-.25))
    return finish(fig, points, pd.DataFrame(summaries))


def make_training_bias_scatter(baseline_df, optimal_df, **kwargs):
    return _calibration(baseline_df, optimal_df, plain=True, separate=True, **kwargs)


def make_training_bias_calibration_overlay(baseline_df, optimal_df, *, annotation_layout="stacked", **kwargs):
    return _calibration(baseline_df, optimal_df, **kwargs)


def make_training_bias_small_multiples_v2(baseline_df, optimal_df, **kwargs):
    return _calibration(baseline_df, optimal_df, separate=True, **kwargs)


def _strips(baseline, optimal, *, model_order=None, one_sided=False, seed=0):
    models = _MAIN if model_order is None else model_order
    prepared = pd.concat([frame.assign(condition=condition) for frame, condition in
                          zip((baseline, optimal), _CONDITIONS)], ignore_index=True)
    prepared = prepared.dropna(subset=["b_occupancy"])
    a, b = prepared.n_assigned_to_A, prepared.n_assigned_to_B
    mask = ((a > 0) & (b == 0)) | ((a == 0) & (b > 0)) if one_sided else (a == 0) & (b == 0)
    points = observations(prepared.loc[mask])
    points["x"], points["y"] = .5, points.b_occupancy
    fig = make_subplots(rows=1, cols=max(1, len(models)),
                        subplot_titles=[MODEL_DISPLAY.get(m, m) for m in models])
    rng, summaries = np.random.default_rng(seed), []
    for col, model in enumerate(models, 1):
        for condition in _CONDITIONS:
            sub = points.loc[(points.model_name == model) & (points.condition == condition)].copy()
            color = _color(model, condition)
            legendgroup = f"{model}:{condition}"
            if one_sided:
                sub["x"] = (sub.n_assigned_to_B > 0).astype(float) + rng.uniform(-.25, .25, len(sub))
                for at, group in sub.groupby(sub.n_assigned_to_B > 0):
                    center = int(at)
                    mean = group.b_occupancy.mean()
                    fig.add_trace(line([center-.25, center+.25], [mean, mean], color=color,
                                       legendgroup=legendgroup), row=1, col=col)
                    summaries.append(dict(kind="mean", model_name=model, condition=condition, x=center, y=mean,
                                          n=len(group)))
            else:
                sub["x"] = .5 + rng.uniform(-.18, .18, len(sub))
                values = sub.y.to_numpy()
                if len(values) >= 3 and np.std(values) > 0:
                    support = np.linspace(0, 1, 200)
                    density = gaussian_kde(values)(support)
                    left = -.05 - .35 * density / density.max()
                    fig.add_trace(line(np.r_[left, np.full(200, -.05)],
                                       np.r_[support, support[::-1]], color=color,
                                       fill="toself", opacity=.25, legendgroup=legendgroup), row=1, col=col)
                    summaries.extend(dict(kind="kde", model_name=model, condition=condition, x=x, y=y)
                                     for x, y in zip(left, support))
            points.loc[sub.index, "x"] = sub.x
            fig.add_trace(scatter(sub, color=color, name=f"{model} · {condition}",
                                  opacity=.7 if one_sided else .5, showlegend=True,
                                  legendgroup=legendgroup),
                          row=1, col=col)
        fig.add_trace(line([-.5, 1.5] if one_sided else [-.45, .73], [.5, .5]), row=1, col=col)
        fig.update_xaxes(range=[-.5, 1.5] if one_sided else [-.45, .73],
                         tickvals=[0, 1] if one_sided else [], ticktext=["Only A", "Only B"] if one_sided else [],
                         row=1, col=col)
    fig.update_yaxes(range=[-.02, 1.1] if one_sided else [-.02, 1.02], title_text="Predicted state-B fraction")
    fig.update_layout(height=500, legend=dict(orientation="h", y=-.2))
    return finish(fig, points, pd.DataFrame(summaries))


def make_training_bias_no_cluster_strip(baseline_df, optimal_df, **kwargs):
    return _strips(baseline_df, optimal_df, **kwargs)


def make_training_bias_one_sided_strip(baseline_df, optimal_df, **kwargs):
    return _strips(baseline_df, optimal_df, one_sided=True, **kwargs)


def make_training_bias_calibration_overlay_legend(*, model_order=None):
    fig = make_subplots(rows=1, cols=1)
    entries = [("Baseline", "#777777")] + [
        (MODEL_DISPLAY[model], MODEL_PALETTE[model])
        for model in (_MAIN if model_order is None else model_order)]
    for name, color in entries:
        trace = scatter(pd.DataFrame(columns=["x", "y", "point_id"]), color=color, name=name)
        trace.x, trace.y, trace.customdata = [None], [None], [[""]]
        trace.meta, trace.showlegend = {"observations": False}, True
        fig.add_trace(trace)
    fig.update_xaxes(visible=False)
    fig.update_yaxes(visible=False)
    fig.update_layout(height=250)
    return finish(fig, pd.DataFrame(columns=["point_id"]))


def make_training_bias_calibration_explainer():
    fig = make_subplots(rows=1, cols=1)
    for x, y, name in (([0, 1], [0, 1], "Training prior"),
                        ([0, 1], [.5, .5], "Guidance target"),
                        ([0, 1], [.25, .75], "Bisector")):
        trace = line(x, y, name=name)
        trace.showlegend = True
        fig.add_trace(trace, row=1, col=1)
    fig.add_annotation(x=.78, y=.81, text="training set prior", textangle=-45, showarrow=False)
    fig.add_annotation(x=.18, y=.52, text="perfect guidance", showarrow=False, yanchor="bottom")
    perpendicular = np.array([-.5, 1]) / np.sqrt(1.25) * .09
    for center in (.1, .9):
        x, y = center, .5 * center + .25
        fig.add_annotation(x=x+perpendicular[0], y=y+perpendicular[1],
                           ax=x-perpendicular[0], ay=y-perpendicular[1],
                           xref="x", yref="y", axref="x", ayref="y",
                           text="", showarrow=True, arrowside="end+start", arrowhead=3,
                           startarrowhead=3, arrowwidth=2, arrowcolor="#222222")
    fig.update_xaxes(range=[-.02, 1.02], title_text="Training state-B fraction")
    fig.update_yaxes(range=[-.02, 1.02], title_text="Predicted state-B fraction", scaleanchor="x")
    fig.update_layout(height=500)
    return finish(fig, pd.DataFrame(columns=["point_id"]))
