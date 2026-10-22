"""Figure 2C reconstruction: read a snapshot, prepare observations, then ECDFs.

The default reproduces notebook 02's protein-mean RSCC and per-region clustering.
The per-region RSCC option is a separate, explicit interpretation of the caption.
No function in this module writes to the input snapshot.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from loaders import (
    RSCC_MAX_RMSD_SOURCES,
    add_simplified_metadata_columns,
    load_lddt,
    load_rscc,
)

METRICS = ("RSCC", "Clustering score")


def load_baseline_snapshot(results_dir: str | Path):
    """Return baseline RSCC/lDDT frames and checksums of the exact source files.

    CSV record numbers are one-based data-record positions, not physical lines.
    Existing loaders supply the same 50/50 filter and score formula as notebook 02.
    """
    root = Path(results_dir).expanduser().resolve()
    frames = {"RSCC": [], "Clustering score": []}
    sources = []
    for source, override in RSCC_MAX_RMSD_SOURCES:
        rscc_name = Path(source).name
        lddt_name = rscc_name.replace("rscc_results_max_rmsd", "lddt_results").replace(
            "rscc_results", "lddt_results"
        )
        for metric, name, loader in [
            ("RSCC", rscc_name, load_rscc),
            ("Clustering score", lddt_name, load_lddt),
        ]:
            path = root / name
            with path.open("rb") as handle:
                before = hashlib.file_digest(handle, "sha256").hexdigest()
            frame = loader(str(path), model_name_override=override)
            frame = add_simplified_metadata_columns(frame)
            frame["source_file"] = name
            frame["source_record"] = frame.index + 1
            frame = frame.loc[frame.guidance_weight == 0].copy()
            keys = ["protein", "model_name", "selection"]
            if frame.duplicated(keys).any():
                raise ValueError(f"{name}: ambiguous duplicate baseline records for {keys}")
            with path.open("rb") as handle:
                after = hashlib.file_digest(handle, "sha256").hexdigest()
            if before != after:
                raise ValueError(f"{name} changed while loading; reload a stable snapshot")
            frames[metric].append(frame)
            sources.append({"file": name, "sha256": before, "model_override": override})
    identity = hashlib.sha256(json.dumps(sources, sort_keys=True).encode()).hexdigest()
    manifest = {
        "snapshot_sha256": identity,
        "source_files": sources,
        "filters": {"map_occupancy": {"A": 0.5, "B": 0.5}, "guidance_weight": 0},
        "clustering_formula": "min(predicted_occupancies) * avg_silhouette_to_ref / 0.5",
        "status": "Current-results reconstruction; historical paper identity unverified",
    }
    return (
        pd.concat(frames["RSCC"], ignore_index=True),
        pd.concat(frames["Clustering score"], ignore_index=True),
        manifest,
    )


def prepare_baseline(rscc, lddt, models, rscc_unit="protein"):
    """Return plotted observations, their source records, and coverage counts."""
    if rscc_unit not in {"protein", "region"}:
        raise ValueError("rscc_unit must be 'protein' or 'region'")
    points, records, coverage = [], [], []
    for metric, frame, column in [
        ("RSCC", rscc, "rscc"),
        ("Clustering score", lddt, "clustering_score"),
    ]:
        for model in sorted(models):
            original = frame.loc[frame.model_name == model].copy()
            per_protein = metric == "RSCC" and rscc_unit == "protein"
            data = original if per_protein else original.loc[
                original.selection_name != "all altloc"
            ].copy()
            data = data.rename(columns={column: "value"})
            keys = ["protein"] if per_protein else ["protein", "selection"]
            data["point_id"] = [
                json.dumps([metric, model, *key], separators=(",", ":"))
                for key in data[keys].itertuples(index=False, name=None)
            ]
            # pandas mean intentionally retains notebook 02's NaN-skipping behavior.
            plotted = data.groupby("point_id", sort=False, as_index=False).agg(
                protein=("protein", "first"),
                selection=("selection", lambda s: "protein mean" if per_protein else s.iloc[0]),
                value=("value", "mean"),
                source_rows=("value", "size"),
            )
            plotted["metric"], plotted["model_name"] = metric, model
            plotted = plotted.loc[np.isfinite(plotted.value)].copy()
            points.append(plotted)
            records.append(data.loc[data.point_id.isin(plotted.point_id), [
                "point_id", "protein", "selection", "value", "source_file", "source_record"
            ]].assign(metric=metric, model_name=model))
            coverage.append({
                "metric": metric,
                "model_name": model,
                "unit": "protein mean" if per_protein else "region",
                "baseline_rows": len(original),
                "excluded_chain_rows": len(original) - len(data),
                "nonfinite_source_values": int((~np.isfinite(data.value)).sum()),
                "plotted_observations": len(plotted),
            })
    columns = ["point_id", "protein", "selection", "value", "source_rows", "metric", "model_name"]
    observations = pd.concat(points, ignore_index=True) if points else pd.DataFrame(columns=columns)
    record_columns = [
        "point_id", "protein", "selection", "value", "source_file", "source_record",
        "metric", "model_name",
    ]
    contributors = pd.concat(records, ignore_index=True) if records else pd.DataFrame(
        columns=record_columns
    )
    return observations, contributors, pd.DataFrame(coverage)


def baseline_ecdf(points):
    """Prepare final, finite ECDF vertices once for both renderers, retaining ties."""
    curves = []
    for (metric, model), group in points.groupby(["metric", "model_name"], sort=False):
        ordered = group.sort_values(["value", "point_id"], kind="stable")
        values = ordered.value.to_numpy()
        # A repeated first x draws the initial jump without an unserializable -inf.
        curve = pd.DataFrame({
            "x": np.r_[values[0], values],
            "y": np.r_[0, np.arange(1, len(values) + 1) / len(values)],
            "point_id": ["", *ordered.point_id],
            "protein": ["", *ordered.protein],
            "selection": ["", *ordered.selection],
        })
        curves.append(curve.assign(metric=metric, model_name=model, n=len(values)))
    columns = ["x", "y", "point_id", "protein", "selection", "metric", "model_name", "n"]
    return pd.concat(curves, ignore_index=True) if curves else pd.DataFrame(columns=columns)


def make_baseline_static(curves):
    """Render prepared ECDF vertices; no metric calculations live in the renderer."""
    import matplotlib.pyplot as plt

    from style import MODEL_DISPLAY, MODEL_PALETTE

    fig, axes = plt.subplots(2, 1, figsize=(10, 5.6))
    for axis, metric in zip(axes, METRICS):
        for model, curve in curves.loc[curves.metric == metric].groupby("model_name"):
            axis.step(curve.x, curve.y, where="post", linewidth=1.8,
                      color=MODEL_PALETTE.get(model, "#888888"),
                      label=MODEL_DISPLAY.get(model, model))
        axis.set(xlabel=metric, ylabel="Proportion", ylim=(-0.05, 1.05))
        axis.spines[["top", "right"]].set_visible(False)
    if len(curves):
        axes[0].legend(
            title="Model", loc="upper left", bbox_to_anchor=(1.01, 1),
            framealpha=0.9, edgecolor="none",
        )
    fig.tight_layout()
    return fig


def make_baseline_interactive(curves):
    """Plotly ECDF with observation identity in hover and selection payloads."""
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    from style import MODEL_DISPLAY, MODEL_PALETTE, MODELS_ALT

    fig = make_subplots(rows=2, cols=1, vertical_spacing=0.17)
    for row, metric in enumerate(METRICS, 1):
        for model, curve in curves.loc[curves.metric == metric].groupby("model_name"):
            fig.add_trace(go.Scatter(
                x=curve.x.tolist(), y=curve.y.tolist(), mode="lines+markers",
                line={"shape": "hv", "color": MODEL_PALETTE.get(model, "#888888"), "width": 2},
                # Native point selection keeps subplot identity; line-only fallback does not.
                marker={"size": 3, "opacity": [0, *([0.35] * (len(curve) - 1))]},
                name=MODEL_DISPLAY.get(model, model), legendgroup=model, showlegend=row == 1,
                meta={"observations": True},
                visible="legendonly" if model in MODELS_ALT else True,
                customdata=curve[["point_id", "protein", "selection", "n"]].values.tolist(),
                hovertemplate=(
                    "%{fullData.name}<br>%{customdata[1]} %{customdata[2]}"
                    "<br>Score: %{x:.4f}<br>Cumulative proportion: %{y:.3f}"
                    "<br>Observations: %{customdata[3]}<extra></extra>"
                ),
            ), row=row, col=1)
        fig.update_xaxes(title_text=metric, row=row, col=1)
        fig.update_yaxes(title_text="Proportion", range=[-0.05, 1.05], row=row, col=1)
    fig.update_layout(
        height=650, template="plotly_white", margin={"t": 100, "b": 45},
        legend={"orientation": "h", "y": 1.02, "yanchor": "bottom",
                "x": 0, "xanchor": "left", "groupclick": "togglegroup"},
        clickmode="event+select", uirevision="baseline-models",
    )
    return fig


def baseline_provenance(manifest, models, rscc_unit):
    """Describe the reconstruction settings, including unresolved science choices."""
    from style import MODELS_ALT

    return {
        **manifest,
        "models": list(models),
        "initially_visible_models": [m for m in models if m not in MODELS_ALT],
        "static_reconstruction_models": [m for m in models if m not in MODELS_ALT],
        "export_scope": (
            "CSV downloads include all loaded models. The separate reconstruction SVG "
            "uses main models. Use the plot toolbar SVG for the current legend-visible view."
        ),
        "rscc_unit": rscc_unit,
        "clustering_unit": "region",
        "distribution": "ECDF (ties retained; no smoothing)",
        "exclude_nan_rscc_proteins": False,
        "analysis_sha256": {
            name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
            for name in (
                "baseline.py", "loaders.py", "style.py", "plot_interaction.py",
                "01_baseline.py", "uv.lock",
            )
        },
        "caveats": [
            "Notebook 02 averages RSCC per protein; the Figure 2C caption describes regions.",
            "Clustering uses the existing minimum across every predicted occupancy entry, "
            "including three/four-entry arrays; state correspondence is not yet adjudicated.",
        ],
    }
