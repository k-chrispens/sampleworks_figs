"""Figure 3E: precomputed ensemble minima, explicitly aggregated across predictors."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from loaders import load_min_altloc_rmsd
from style import MODELS_ALT, OPTIMAL_GW
from transforms import format_residue_range

RMSD_COLUMNS = ["min_rmsd_to_A", "min_rmsd_to_B"]
RINGS = (
    ("Biological examples", "curated", "#1565C0"),
    ("RSCC far-from-ceiling examples", "far", "#F57C00"),
    ("RSCC-decreased examples (zoom only)", "decreased", "#D32F2F"),
)
RED_REGIONS = (
    "7Z0E 11-13", "7Z0E 18-20", "7Z0E 54-56", "7Z0E 82-84",
    "7Z0E 122-124", "7Z0E 181-183", "7Z0E 206-208", "7Z0E 214-216",
)


def load_min_rmsd_snapshot(results_dir):
    """Read the two main-model exports, retaining failures and original row ordinals."""
    frames, sources = [], []
    root = Path(results_dir).expanduser().resolve()
    for suffix in ("boltz_protenix", "rf3_fixed"):
        path = root / f"min_altloc_rmsd_results_{suffix}.csv"
        with path.open("rb") as handle:
            before = hashlib.file_digest(handle, "sha256").hexdigest()
        frame = load_min_altloc_rmsd(str(path), keep_failed=True)
        if not frame.ensemble_size.eq(8).all():
            raise ValueError(f"{path.name}: Figure 3E requires ensemble_size=8")
        frame["source_file"], frame["source_record"] = path.name, frame.index + 1
        with path.open("rb") as handle:
            after = hashlib.file_digest(handle, "sha256").hexdigest()
        if before != after:
            raise ValueError(f"{path.name} changed while loading; reload a stable snapshot")
        frames.append(frame)
        sources.append({"file": path.name, "sha256": before})
    return pd.concat(frames, ignore_index=True), {
        "source_files": sources,
        "snapshot_sha256": hashlib.sha256(json.dumps(sources, sort_keys=True).encode()).hexdigest(),
    }


def prepare_min_rmsd(data, aggregation="median", min_models=3):
    """Aggregate complete guided model pairs; expose failures and missing conditions."""
    if aggregation not in ("median", "mean"):
        raise ValueError("aggregation must be median or mean")
    if min_models < 1:
        raise ValueError("min_models must be positive")
    data = data.loc[~data.model_name.isin(MODELS_ALT) & data.selection_name.ne("all altloc")].copy()
    keys = ["protein", "selection", "model_name"]
    if data[keys].isna().any().any():
        raise ValueError("Missing protein, selection, or model identity")
    weights = data.model_name.map(OPTIMAL_GW)
    if weights.isna().any():
        raise ValueError("A model has no declared guided comparison weight")
    selected = data.loc[data.guidance_weight == weights].copy()
    if selected.duplicated(keys).any():
        raise ValueError("Duplicate model/region records at the chosen guidance weight")
    selected["included"] = np.isfinite(selected[RMSD_COLUMNS]).all(axis=1)
    if (selected[RMSD_COLUMNS] < 0).any().any():
        raise ValueError("RMSD cannot be negative")
    selected["point_id"] = [
        json.dumps([p, s]) for p, s in selected[["protein", "selection"]].itertuples(index=False, name=None)
    ]
    selected["label"] = [
        f"{p.upper()} {format_residue_range(s)}"
        for p, s in selected[["protein", "selection"]].itertuples(index=False, name=None)
    ]
    audit = []
    for model, group in data.groupby("model_name"):
        chosen = selected.loc[selected.model_name == model]
        audit.append({
            "model_name": model,
            "regions_at_any_weight": len(group[["protein", "selection"]].drop_duplicates()),
            "regions_at_chosen_weight": len(chosen),
            "missing_chosen_weight": len(group[["protein", "selection"]].drop_duplicates()) - len(chosen),
            "finite_pairs": int(chosen.included.sum()),
            "nonfinite_pairs": int((~chosen.included).sum()),
        })
    complete = selected.loc[selected.included]
    points = complete.groupby(["protein", "selection", "label"], as_index=False).agg(
        point_id=("point_id", "first"),
        min_rmsd_to_A=("min_rmsd_to_A", aggregation),
        min_rmsd_to_B=("min_rmsd_to_B", aggregation),
        n_models=("model_name", "count"),
    )
    points["plotted"] = points.n_models >= min_models
    points["zoom_visible"] = points.plotted & points[RMSD_COLUMNS].le(3).all(axis=1)
    points["curated"] = points.label.isin(("6B8X 181-183", "5MHX 162-164", "5SOP 131-133"))
    points["far"] = points.label.isin(("2QCV 2-4", "4OLE 424-426", "6RP1 5"))
    points["decreased"] = points.label.isin(RED_REGIONS)
    return points, selected, pd.DataFrame(audit)


def highlighted_min_rmsd(points):
    """One row per visible circled region; red circles are confined to the zoom panel."""
    table = points.loc[
        points.plotted & (points.curated | points.far | (points.decreased & points.zoom_visible))
    ].copy()
    table["rings"] = [
        "; ".join(f"{name} ({color})" for name, flag, color in RINGS
                  if getattr(row, flag) and (flag != "decreased" or row.zoom_visible))
        for row in table.itertuples()
    ]
    return table


def _limits(points):
    shown = points.loc[points.plotted, RMSD_COLUMNS]
    return (0, max(0.1, float(shown.max().max()) * 1.05)) if len(shown) else (0, 1)


def make_min_rmsd_static(points, aggregation="median"):
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(12, 6), layout="constrained")
    shown = points.loc[points.plotted]
    for ax, title, limits in zip(axes, ("Auto-fit", "Zoomed"), (_limits(points), (0, 3))):
        ax.scatter(shown.min_rmsd_to_A, shown.min_rmsd_to_B, s=18, color="#888888",
                   alpha=0.45, edgecolors="white", linewidths=0.4)
        visible = shown.loc[shown[RMSD_COLUMNS].le(limits[1]).all(axis=1)]
        for _, flag, color in RINGS:
            if flag == "decreased" and title != "Zoomed":
                continue
            group = visible.loc[visible[flag]]
            ax.scatter(group.min_rmsd_to_A, group.min_rmsd_to_B, s=80, facecolors="none",
                       edgecolors=color, linewidths=1.6)
        ax.plot(limits, limits, ls="--", color="#AAAAAA", lw=1)
        ax.set(xlim=limits, ylim=limits, title=title, aspect="equal",
               xlabel="Min RMSD to altloc A (Å)", ylabel="Min RMSD to altloc B (Å)")
        ax.text(.03, .97, "closer to A", transform=ax.transAxes, va="top",
                fontstyle="italic", color="#888888")
        ax.text(.97, .03, "closer to B", transform=ax.transAxes, ha="right",
                fontstyle="italic", color="#888888")
    fig.suptitle(f"Minimum RMSD to each altloc ({aggregation} across main models)")
    return fig


def make_min_rmsd_interactive(points, aggregation="median"):
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    fig = make_subplots(rows=1, cols=2, subplot_titles=("Auto-fit", "Zoomed"),
                        horizontal_spacing=0.12)
    shown = points.loc[points.plotted]
    for col, limits in enumerate((_limits(points), (0, 3)), start=1):
        for name, flag, color in [("Regions", None, "#888888"), *RINGS]:
            if flag == "decreased" and col == 1:
                continue
            group = shown if flag is None else shown.loc[
                shown[flag] & shown[RMSD_COLUMNS].le(limits[1]).all(axis=1)
            ]
            if group.empty:
                continue
            fig.add_trace(go.Scatter(
                x=group.min_rmsd_to_A.tolist(), y=group.min_rmsd_to_B.tolist(), mode="markers",
                name=name, legendgroup=name, showlegend=(col == 1 or flag == "decreased"),
                meta={"observations": flag is None},
                marker={"size": 5 if flag is None else 11, "color": color,
                        "symbol": "circle" if flag is None else "circle-open",
                        "opacity": .45 if flag is None else 1,
                        "line": {"width": 0 if flag is None else 1.6, "color": color}},
                customdata=group[["point_id", "label", "selection", "n_models"]].values.tolist(),
                hovertemplate=("%{customdata[1]}<br>%{customdata[2]}"
                               "<br>RMSD to A: %{x:.4f} Å<br>RMSD to B: %{y:.4f} Å"
                               "<br>Contributing models: %{customdata[3]}<extra></extra>"),
            ), row=1, col=col)
        fig.add_shape(type="line", x0=limits[0], y0=limits[0], x1=limits[1], y1=limits[1],
                      line={"dash": "dash", "color": "#AAAAAA"}, row=1, col=col)
        fig.update_xaxes(range=limits, constrain="domain", title_text="Min RMSD to altloc A (Å)",
                         row=1, col=col)
        fig.update_yaxes(range=limits, title_text="Min RMSD to altloc B (Å)",
                         scaleanchor="x" if col == 1 else "x2", scaleratio=1,
                         constrain="domain", row=1, col=col)
    fig.update_layout(template="plotly_white", height=650, margin={"t": 120},
                      title={"text": f"Minimum RMSD ({aggregation} across main models)",
                             "x": .5, "xanchor": "center", "y": .99, "yanchor": "top"},
                      legend={"orientation": "h", "y": 1.05, "yanchor": "bottom",
                              "groupclick": "togglegroup"},
                      clickmode="event+select", uirevision=f"min-rmsd-{aggregation}")
    return fig


def min_rmsd_provenance(manifest, aggregation="median"):
    return {
        **manifest, "aggregation": aggregation,
        "status": "Current-results reconstruction; historical paper computation unverified",
        "within_model": "Precomputed minimum over the eight ensemble members, separately for A/B",
        "across_models": f"{aggregation} of complete per-model A/B pairs; at least three main models",
        "chosen_guidance_weights": {m: w for m, w in OPTIMAL_GW.items() if m not in MODELS_ALT},
        "filters": {"map_occupancy": {"A": .5, "B": .5}, "ensemble_size": 8,
                    "minimum_models": 3, "exclude_nan_rscc_proteins": False},
        "caveats": [
            "Median reproduces the existing notebook; mean is the caption-described alternative.",
            "Neither option establishes the historical published calculation.",
            "Non-finite A/B pairs are excluded before aggregation and retained in the contribution table.",
            "Expected evaluations not present at any guidance weight cannot be inferred.",
        ],
        "highlights": {"blue": ["6B8X 181-183", "5MHX 162-164", "5SOP 131-133"],
                       "orange": ["2QCV 2-4", "4OLE 424-426", "6RP1 5"],
                       "red_zoom_only": RED_REGIONS, "fixed_segment_labels": False},
        "export_scope": "Toolbar SVG: visible view. Reconstruction ZIP: both panels, "
                        "circled-region table and recipe for the selected aggregation.",
        "analysis_sha256": {
            name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
            for name in ("min_rmsd.py", "loaders.py", "style.py", "transforms.py",
                         "plot_interaction.py", "05_min_rmsd.py", "uv.lock")
        },
    }
