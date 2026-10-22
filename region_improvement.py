"""Figures 3C/D: paired measurements, region medians, highlights, and marginal KDEs.

Reconstructs notebook 02: subtract the marginal medians, do not take the median
of per-model differences. Model-row counts retain the legacy minimum-three rule.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

from loaders import RSCC_MAX_RMSD_SOURCES, add_simplified_metadata_columns, load_lddt, load_rscc
from style import MODELS_ALT, OPTIMAL_GW
from transforms import format_residue_range

POINT_COLOR = "#377EB8"
METRIC_LABELS = {"rscc": "RSCC", "clustering_score": "Clustering score"}
CURATED = ("6B8X 181-183", "5MHX 162-164", "5SOP 131-133")
HIGHLIGHTS = (
    ("Most improved", "top", "#2E7D32"),
    ("Decreased", "decreased", "#D32F2F"),
    ("Far from ceiling", "far", "#F57C00"),
    ("Biological examples", "curated", "#1565C0"),
)


def load_improvement_snapshot(results_dir, metric="rscc"):
    """Read/fingerprint only the two main-model files for the requested metric."""
    if metric not in METRIC_LABELS:
        raise ValueError(f"Unknown improvement metric: {metric}")
    root = Path(results_dir).expanduser().resolve()
    frames, sources = [], []
    for source, override in RSCC_MAX_RMSD_SOURCES:
        if override in MODELS_ALT:
            continue
        filename = Path(source).name
        if metric == "clustering_score":
            filename = filename.replace("rscc_results_max_rmsd_", "lddt_results_")
        path = root / filename
        with path.open("rb") as handle:
            before = hashlib.file_digest(handle, "sha256").hexdigest()
        loader = load_rscc if metric == "rscc" else load_lddt
        frame = loader(str(path), model_name_override=override)
        if metric == "clustering_score":
            sizes = pd.read_csv(path, usecols=["ensemble_size"]).ensemble_size
            frame["ensemble_size"] = sizes.loc[frame.index]
            frame = add_simplified_metadata_columns(frame)
        if not frame.ensemble_size.eq(8).all():
            raise ValueError(f"{path.name}: region improvement requires ensemble_size=8")
        frame["source_file"], frame["source_record"] = path.name, frame.index + 1
        with path.open("rb") as handle:
            after = hashlib.file_digest(handle, "sha256").hexdigest()
        if before != after:
            raise ValueError(f"{path.name} changed while loading; reload a stable snapshot")
        frames.append(frame)
        sources.append({"file": path.name, "sha256": before})
    return pd.concat(frames, ignore_index=True), {
        "snapshot_sha256": hashlib.sha256(json.dumps(sources, sort_keys=True).encode()).hexdigest(),
        "source_files": sources,
    }


def prepare_improvement(data, metric="rscc", min_models=3):
    """Return all region summaries, paired model/source contributions, and coverage."""
    if metric not in METRIC_LABELS:
        raise ValueError(f"Unknown improvement metric: {metric}")
    if min_models < 1:
        raise ValueError("min_models must be positive")
    data = data.loc[
        ~data.model_name.isin(MODELS_ALT) & (data.selection_name != "all altloc")
    ].copy()
    weights = data.model_name.map(OPTIMAL_GW)
    if weights.isna().any():
        raise ValueError("A model has no declared guided comparison weight")
    keys = ["protein", "selection", "model_name"]
    if data[keys].isna().any().any():
        raise ValueError("Missing protein, selection, or model identity")
    columns = keys + [metric, "source_file", "source_record"]
    if metric == "clustering_score":
        columns += ["occupancies", "avg_silhouette_to_ref",
                    "min_predicted_occupancy", "min_map_occupancy"]
    baseline = data.loc[data.guidance_weight == 0, columns].rename(
        columns={metric: "baseline_value"}
    )
    guided = data.loc[data.guidance_weight == weights, columns + ["guidance_weight"]].rename(
        columns={metric: "guided_value", "guidance_weight": "selected_guidance_weight"}
    )
    joined = baseline.merge(
        guided, on=keys, how="outer", suffixes=["_baseline", "_guided"],
        indicator=True, validate="one_to_one",
    )
    audit = []
    for model, group in joined.groupby("model_name"):
        both = group.loc[group["_merge"] == "both"]
        audit.append({
            "model_name": model,
            "paired_regions": len(both),
            "missing_guided": int(group["_merge"].eq("left_only").sum()),
            "missing_baseline": int(group["_merge"].eq("right_only").sum()),
            "finite_pairs": int(np.isfinite(both[["baseline_value", "guided_value"]]).all(axis=1).sum()),
        })
    pairs = joined.loc[joined["_merge"] == "both"].drop(columns="_merge").copy()
    pairs["delta_value"] = pairs.guided_value - pairs.baseline_value
    pairs["finite_pair"] = np.isfinite(pairs[["baseline_value", "guided_value"]]).all(axis=1)
    pairs["point_id"] = [
        json.dumps([p, s]) for p, s in pairs[["protein", "selection"]].itertuples(index=False, name=None)
    ]
    pairs["label"] = [
        f"{p.upper()} {format_residue_range(s)}"
        for p, s in pairs[["protein", "selection"]].itertuples(index=False, name=None)
    ]
    points = pairs.groupby(["protein", "selection", "label"], as_index=False).agg(
        point_id=("point_id", "first"),
        baseline_value=("baseline_value", "median"),
        guided_value=("guided_value", "median"),
        n_model_rows=("model_name", "count"),
        n_finite_pairs=("finite_pair", "sum"),
        n_baseline_values=("baseline_value", lambda s: np.isfinite(s).sum()),
        n_guided_values=("guided_value", lambda s: np.isfinite(s).sum()),
    )
    points["improvement"] = points.guided_value - points.baseline_value
    points["meets_minimum"] = points.n_model_rows >= min_models
    points["plotted"] = points.meets_minimum & np.isfinite(
        points[["baseline_value", "guided_value", "improvement"]]
    ).all(axis=1)
    for _, flag, _ in HIGHLIGHTS:
        points[flag] = False
    shown = points.loc[points.plotted]
    top = shown.improvement.nlargest(3)
    top_ids = top.index[top > 0].tolist()
    decreased_ids = shown.index[shown.improvement < 0].tolist()
    # The legacy renderer takes five candidates first, then excludes already-highlighted ones.
    distance = (1 - shown.guided_value).nlargest(5 if metric == "rscc" else 3)
    far_ids = [
        i for i in distance.index
        if distance.loc[i] > 0 and i not in set(top_ids + decreased_ids)
    ]
    curated_ids = shown.index[shown.label.isin(CURATED)].tolist()
    for flag, ids in [("top", top_ids), ("decreased", decreased_ids),
                      ("far", far_ids), ("curated", curated_ids)]:
        points.loc[ids, flag] = True
    return points, pairs, pd.DataFrame(audit)


def highlighted_regions(points):
    """Companion key for the unlabeled rings, including overlapping categories."""
    flags = [flag for _, flag, _ in HIGHLIGHTS]
    table = points.loc[points.plotted & points[flags].any(axis=1)].copy()
    table["rings"] = [
        "; ".join(f"{name} ({color})" for name, flag, color in HIGHLIGHTS if row[flag])
        for _, row in table.iterrows()
    ]
    return table


def prepare_density(points):
    """Prepare the legacy Gaussian KDE curves once; report degenerate marginals."""
    curves, audit = [], []
    eligible = points.loc[points.meets_minimum]
    for series, column in [("Baseline", "baseline_value"), ("Guided", "guided_value")]:
        values = eligible.loc[np.isfinite(eligible[column]), column].to_numpy()
        status = "drawn"
        if len(values) < 2 or np.var(values) == 0:
            status = "unavailable: fewer than two distinct finite values"
        else:
            kernel = gaussian_kde(values)
            kernel.set_bandwidth(kernel.factor * 0.5)
            bandwidth = float(np.sqrt(kernel.covariance.squeeze()))
            support = np.linspace(values.min() - 3 * bandwidth, values.max() + 3 * bandwidth, 200)
            curves.append(pd.DataFrame({
                "series": series, "x": support, "density": kernel(support),
            }))
        audit.append({"series": series, "finite_regions": len(values), "status": status})
    return (
        pd.concat(curves, ignore_index=True) if curves else pd.DataFrame(
            columns=["series", "x", "density"]
        ),
        pd.DataFrame(audit),
    )


def _ylim(points):
    shown = points.loc[points.plotted]
    return (min(0, float(shown.improvement.min())) - 0.02, 1) if len(shown) else (-0.05, 1)


def make_improvement_static(points, density, metric="rscc"):
    """Render prepared coordinates, rings, and density; identities stay in hover/tables."""
    import matplotlib.pyplot as plt

    fig, (ax, marginal) = plt.subplots(
        2, 1, figsize=(6, 7.5), sharex=True, layout="constrained",
        gridspec_kw={"height_ratios": [4, 1], "hspace": 0.05}
    )
    shown = points.loc[points.plotted]
    label = METRIC_LABELS[metric]
    ax.scatter(shown.baseline_value, shown.improvement, s=18, color=POINT_COLOR,
               alpha=0.3, edgecolors="none", zorder=4)
    for _, flag, color in HIGHLIGHTS:
        group = shown.loc[shown[flag]]
        ax.scatter(group.baseline_value, group.improvement, s=80, facecolors="none",
                   edgecolors=color, linewidths=1.6, zorder=6 if flag == "curated" else 5)
    ax.axhline(0, color="#AAAAAA", ls="--", lw=1)
    ax.plot([0, 1], [1, 0], color="#AAAAAA", ls="--", lw=0.8)
    ax.set(xlim=(0, 1), ylim=_ylim(points), ylabel="Improvement (guided median − baseline median)",
           title=f"{label} improvement over baseline")
    ax.tick_params(axis="x", labelbottom=False, length=0)
    for series, color in [("Baseline", "#888888"), ("Guided", POINT_COLOR)]:
        curve = density.loc[density.series == series]
        if len(curve):
            marginal.plot(curve.x, curve.density, color=color, lw=1.5, label=series)
            marginal.fill_between(curve.x, curve.density, color=color, alpha=0.3)
    marginal.set(xlabel=f"Baseline {label} (median across main models)", ylabel="Density", yticks=[])
    if len(density):
        marginal.legend(loc="upper right", framealpha=0.8, edgecolor="none")
    return fig


def make_improvement_interactive(points, density, metric="rscc"):
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    fig = make_subplots(rows=2, cols=1, shared_xaxes=True,
                        row_heights=[4, 1], vertical_spacing=0.06)
    label = METRIC_LABELS[metric]
    shown = points.loc[points.plotted]
    for name, flag, color in [("Regions", None, POINT_COLOR), *HIGHLIGHTS]:
        group = shown if flag is None else shown.loc[shown[flag]]
        if group.empty:
            continue
        fig.add_trace(go.Scatter(
            x=group.baseline_value.tolist(), y=group.improvement.tolist(), mode="markers",
            name=name, legendgroup=name,
            meta={"observations": flag is None},
            marker={"size": 5 if flag is None else 11,
                    "symbol": "circle" if flag is None else "circle-open",
                    "color": color, "opacity": 0.3 if flag is None else 1,
                    "line": {"width": 0 if flag is None else 1.6, "color": color}},
            customdata=group[[
                "point_id", "label", "selection", "guided_value", "n_model_rows", "n_finite_pairs"
            ]].values.tolist(),
            hovertemplate=(
                "%{customdata[1]}<br>%{customdata[2]}<br>Baseline median: %{x:.4f}"
                "<br>Guided median: %{customdata[3]:.4f}<br>Improvement: %{y:.4f}"
                "<br>Model rows / finite pairs: %{customdata[4]} / %{customdata[5]}<extra></extra>"
            ),
        ), row=1, col=1)
    for series, color in [("Baseline", "#888888"), ("Guided", POINT_COLOR)]:
        curve = density.loc[density.series == series]
        if len(curve):
            fig.add_trace(go.Scatter(
                x=curve.x.tolist(), y=curve.density.tolist(), mode="lines", name=series + " density",
                line={"color": color}, fill="tozeroy", opacity=0.35,
                customdata=[[""]] * len(curve), hoverinfo="skip",
            ), row=2, col=1)
    fig.add_hline(y=0, line_dash="dash", line_color="#AAAAAA", row=1, col=1)
    fig.add_shape(type="line", x0=0, y0=1, x1=1, y1=0, xref="x", yref="y",
                  line={"dash": "dash", "color": "#AAAAAA"})
    fig.update_xaxes(range=[0, 1])
    fig.update_xaxes(title_text=f"Baseline {label} (median across main models)", row=2, col=1)
    fig.update_yaxes(title_text="Improvement (difference of medians)", range=_ylim(points), row=1, col=1)
    fig.update_yaxes(title_text="Density", showticklabels=False, row=2, col=1)
    fig.update_layout(template="plotly_white", height=760, margin={"t": 100, "b": 65},
                      legend={"orientation": "h", "y": 1.02, "yanchor": "bottom",
                              "groupclick": "togglegroup"},
                      clickmode="event+select", uirevision=f"{metric}-improvement")
    return fig


def improvement_provenance(manifest, metric="rscc"):
    return {
        **manifest,
        "metric": metric,
        "clustering_formula": "min(all predicted occupancy entries) * avg_silhouette_to_ref / 0.5",
        "status": "Notebook 02 reconstruction; historical paper identity unverified",
        "models": [m for m in OPTIMAL_GW if m not in MODELS_ALT],
        "chosen_guidance_weights": {m: gw for m, gw in OPTIMAL_GW.items() if m not in MODELS_ALT},
        "filters": {"map_occupancy": {"A": 0.5, "B": 0.5}, "ensemble_size": 8,
                    "minimum_model_rows": 3, "exclude_nan_rscc_proteins": False},
        "aggregation": "x = median(baseline); y = median(guided) - median(baseline), not median(delta)",
        "density": "Gaussian KDE per marginal, Scott bandwidth * 0.5, cut=3, gridsize=200",
        "highlights": {"top": 3, "far_candidates": 5 if metric == "rscc" else 3,
                       "decreased": "all", "curated": CURATED},
        "highlight_labels": "No fixed segment text; hover/selection and companion highlighted-region CSV",
        "export_scope": "Legend controls display only, not contributing models. Toolbar SVG: visible "
                        "view; reconstruction SVG/CSVs: full fixed main-model analysis.",
        "caveats": [
            "The legacy minimum counts paired model rows, including rows with NaN scores.",
            "Medians retain pandas' NaN-skipping behavior; finite contributor counts are exposed.",
            "Each marginal KDE excludes non-finite values independently.",
            "Only regions with finite scatter coordinates are drawn.",
        ],
        "analysis_sha256": {
            name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
            for name in ("region_improvement.py", "loaders.py", "style.py", "transforms.py",
                         "plot_interaction.py",
                         "03_rscc_improvement.py" if metric == "rscc" else "04_clustering_improvement.py",
                         "uv.lock")
        },
    }
