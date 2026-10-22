"""Figure 3B: reconstruct notebook 03's protein deltas and median/IQR trajectories.

Scientific reductions happen here, not in the two renderers. Missing-value means
and marginal quantiles intentionally retain notebook 03's behavior and are audited.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from loaders import RSCC_MAX_RMSD_SOURCES, load_rscc, tag_path_dataframe
from style import MODEL_DISPLAY, MODEL_PALETTE, MODELS_ALT, OPTIMAL_GW
from transforms import make_delta_frame

KEYS = ["protein", "model_name", "guidance_weight", "guidance_strength"]
GROUPS = ["model_name", "guidance_weight"]
YTICKS = np.array([-1, 0, 1, 10, 100, 1000])


def symlog(values):
    """Matplotlib's base-10 symlog with linthresh=linscale=1, including zero."""
    values = np.asarray(values, dtype=float)
    magnitude = np.abs(values)
    return np.where(
        magnitude <= 1, values / 0.9,
        np.sign(values) * (1 / 0.9 + np.log10(np.maximum(magnitude, 1))),
    )


def _tag_members(frame, path_column, override):
    # Parse each ensemble path once, not once for every member.
    paths = frame[[path_column]].drop_duplicates().reset_index(drop=True)
    try:
        tagged = tag_path_dataframe(paths, filename_col=path_column)
    except (AttributeError, TypeError) as exc:
        raise ValueError("Geometry paths must use the existing ens8_gw export layout") from exc
    result = frame.rename(columns={path_column: "filename"}).merge(
        tagged, on="filename", validate="many_to_one"
    )
    if override is not None:
        result["model_name"] = override
    return result


def load_guidance_snapshot(results_dir):
    """Load 20 existing CSVs read-only; retain source identities and join diagnostics."""
    root = Path(results_dir).expanduser().resolve()
    sources, joins = [], []
    frames = {"rscc": [], "clashes": [], "bonds": []}

    def read(name, loader=pd.read_csv, **kwargs):
        path = root / name
        with path.open("rb") as handle:
            before = hashlib.file_digest(handle, "sha256").hexdigest()
        frame = loader(str(path), **kwargs)
        frame["source_file"], frame["source_record"] = name, frame.index + 1
        with path.open("rb") as handle:
            after = hashlib.file_digest(handle, "sha256").hexdigest()
        if before != after:
            raise ValueError(f"{name} changed while loading; reload a stable snapshot")
        sources.append({"file": name, "sha256": before})
        return frame

    for source, override in RSCC_MAX_RMSD_SOURCES:
        name = Path(source).name
        suffix = name.removeprefix("rscc_results_max_rmsd").removeprefix("rscc_results")
        rscc = read(name, load_rscc, model_name_override=override)
        if not rscc.ensemble_size.eq(8).all():
            raise ValueError(f"{name}: this reconstruction requires ensemble_size=8")
        if rscc.duplicated(["protein", "model_name", "guidance_weight", "selection"]).any():
            raise ValueError(f"{name}: ambiguous duplicate RSCC observations")
        frames["rscc"].append(rscc)
        clash = read("clashscore_metrics" + suffix)
        if clash.duplicated(["model_name", "model_id"]).any():
            raise ValueError(f"{suffix}: duplicate clashscore members")
        frames["clashes"].append(_tag_members(clash, "model_name", override))
        angles = read("bond_angle_violation_fractions" + suffix)
        lengths = read("bond_length_violation_fractions" + suffix)
        both = angles.merge(
            lengths, on=["model", "model_n"], suffixes=["_bond_angles", "_bond_lengths"],
            how="outer", indicator=True, validate="one_to_one",
        )
        joins.append({
            "source_group": suffix.removesuffix(".csv").lstrip("_"),
            "model_override": override,
            "angle_only_members": int(both["_merge"].eq("left_only").sum()),
            "length_only_members": int(both["_merge"].eq("right_only").sum()),
        })
        frames["bonds"].append(_tag_members(
            both.loc[both["_merge"] == "both"].drop(columns="_merge"), "model", override
        ))
    manifest = {
        "snapshot_sha256": hashlib.sha256(json.dumps(sources, sort_keys=True).encode()).hexdigest(),
        "source_files": sources,
        "bond_member_joins_all_occupancies": joins,
    }
    return {key: pd.concat(parts, ignore_index=True) for key, parts in frames.items()}, manifest


def prepare_guidance(frames, min_n=3):
    """Return protein-level contributions, final plotted summaries, and pairing coverage."""
    metrics = [
        ("clashes", "clashscore"),
        ("bonds", "outlier_fraction_bond_angles"),
        ("bonds", "outlier_fraction_bond_lengths"),
        ("rscc", "rscc"),
    ]
    deltas, coverage = {}, []
    for family, metric in metrics:
        frame = frames[family]
        if frame.guidance_strength.isna().any():
            raise ValueError("Unrecognized guidance weight: define its analysis meaning before plotting")
        stats = frame.groupby(KEYS, as_index=False)[metric].mean()
        delta = make_delta_frame(stats, metric)
        deltas[metric] = delta[KEYS + [metric, metric + "_baseline", "delta_" + metric]]
        candidates = stats.loc[stats.guidance_weight != 0].groupby(GROUPS).size()
        paired = delta.groupby(GROUPS).size()
        finite = delta.loc[np.isfinite(delta["delta_" + metric])].groupby(GROUPS).size()
        audit = pd.DataFrame({"candidate_proteins": candidates, "paired_proteins": paired,
                              "finite_deltas": finite}).fillna(0).astype(int).reset_index()
        audit["missing_baseline"] = audit.candidate_proteins - audit.paired_proteins
        coverage.append(audit.assign(metric=metric))
    proteins = deltas["clashscore"]
    for metric in ["outlier_fraction_bond_angles", "outlier_fraction_bond_lengths", "rscc"]:
        proteins = proteins.merge(
            deltas[metric], on=KEYS, how="left" if metric == "rscc" else "inner",
            validate="one_to_one",
        )
    terms = proteins[[
        "delta_clashscore", "delta_outlier_fraction_bond_lengths", "delta_outlier_fraction_bond_angles"
    ]].mul([1, 100, 100])
    proteins["geometry_penalty"] = terms.mean(axis=1)
    proteins["geometry_terms"] = terms.notna().sum(axis=1)
    proteins["point_id"] = [
        json.dumps([m, float(gw)]) for m, gw in proteins[GROUPS].itertuples(index=False, name=None)
    ]
    summary = proteins.groupby(GROUPS, as_index=False).agg(
        point_id=("point_id", "first"),
        x_med=("delta_rscc", "median"),
        x_q1=("delta_rscc", lambda s: s.quantile(0.25)),
        x_q3=("delta_rscc", lambda s: s.quantile(0.75)),
        y_med=("geometry_penalty", "median"),
        y_q1=("geometry_penalty", lambda s: s.quantile(0.25)),
        y_q3=("geometry_penalty", lambda s: s.quantile(0.75)),
        n=("protein", "count"),
        n_rscc=("delta_rscc", lambda s: np.isfinite(s).sum()),
        n_geometry=("geometry_penalty", lambda s: np.isfinite(s).sum()),
        incomplete_geometry=("geometry_terms", lambda s: (s < 3).sum()),
    ).sort_values(GROUPS)
    summary["plotted"] = (summary.n >= min_n) & np.isfinite(
        summary[["x_med", "x_q1", "x_q3", "y_med", "y_q1", "y_q3"]]
    ).all(axis=1)
    for field in ["y_med", "y_q1", "y_q3"]:
        summary[field + "_plot"] = symlog(summary[field])
    summary["marker_area"] = 0.0
    for _, group in summary.loc[summary.plotted].groupby("model_name"):
        lo, hi = group.guidance_weight.min(), group.guidance_weight.max()
        summary.loc[group.index, "marker_area"] = (
            30 + (group.guidance_weight - lo) / (hi - lo) * 190 if hi > lo else 125
        )
    summary["chosen_weight"] = np.isclose(
        summary.guidance_weight, summary.model_name.map(OPTIMAL_GW)
    )
    return proteins, summary, pd.concat(coverage, ignore_index=True)


def make_guidance_static(summary):
    """Draw final median/IQR coordinates; the y coordinates are already symlog-transformed."""
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 6))
    for model, group in summary.loc[summary.plotted].groupby("model_name"):
        color = MODEL_PALETTE.get(model, "#888888")
        ax.plot(np.r_[0, group.x_med], np.r_[0, group.y_med_plot], color=color,
                label=MODEL_DISPLAY.get(model, model), lw=1.6)
        ax.errorbar(group.x_med, group.y_med_plot,
                    xerr=[group.x_med - group.x_q1, group.x_q3 - group.x_med],
                    yerr=[group.y_med_plot - group.y_q1_plot, group.y_q3_plot - group.y_med_plot],
                    fmt="none", ecolor=color, alpha=0.35, capsize=3)
        ax.scatter(group.x_med, group.y_med_plot, s=group.marker_area, c=color,
                   edgecolors=np.where(group.chosen_weight, "#111111", "white"),
                   linewidths=np.where(group.chosen_weight, 2, 0.5), zorder=4)
        for point in group.itertuples():
            ax.annotate(f"{point.guidance_weight:g}", (point.x_med, point.y_med_plot),
                        xytext=(5, 5), textcoords="offset points", fontsize=9, zorder=5,
                        bbox={"boxstyle": "round,pad=0.12", "fc": "white",
                              "ec": "none", "alpha": 0.7})
    ax.set(xlim=(-0.05, 0.25), ylim=symlog([-1, 1000]),
           xlabel="Δ RSCC vs. baseline", ylabel="Δ Geometry penalty vs. baseline")
    ax.set_yticks(symlog(YTICKS), [str(v) for v in YTICKS])
    ax.axvline(0, color="#555555", ls="--", lw=0.8, alpha=0.5)
    ax.axhline(0, color="#555555", ls="--", lw=0.8, alpha=0.5)
    if summary.plotted.any():
        ax.legend(loc="upper left", frameon=False)
    fig.tight_layout()
    return fig


def make_guidance_interactive(summary):
    """Legend-controlled trajectories with raw-valued hover on the symlog display."""
    import plotly.graph_objects as go

    fig = go.Figure()
    for model, group in summary.loc[summary.plotted].groupby("model_name"):
        fig.add_trace(go.Scatter(
            x=[0, *group.x_med], y=[0, *group.y_med_plot], mode="lines+markers",
            name=MODEL_DISPLAY.get(model, model), legendgroup=model,
            meta={"observations": True},
            visible="legendonly" if model in MODELS_ALT else True,
            line={"color": MODEL_PALETTE.get(model, "#888888")},
            marker={
                "size": [0, *(np.sqrt(group.marker_area) * 4 / 3)],
                "line": {"color": ["white", *np.where(group.chosen_weight, "#111111", "white")],
                         "width": [0, *np.where(group.chosen_weight, 2, 0.5)]},
            },
            error_x={"type": "data", "symmetric": False,
                     "array": [0, *(group.x_q3 - group.x_med)],
                     "arrayminus": [0, *(group.x_med - group.x_q1)]},
            error_y={"type": "data", "symmetric": False,
                     "array": [0, *(group.y_q3_plot - group.y_med_plot)],
                     "arrayminus": [0, *(group.y_med_plot - group.y_q1_plot)]},
            customdata=[["", 0, 0, 0, 0, 0], *group[
                ["point_id", "guidance_weight", "y_med", "n", "n_rscc", "n_geometry"]
            ].values.tolist()],
            hovertemplate=(
                "%{fullData.name}<br>Guidance: %{customdata[1]}<br>Δ RSCC: %{x:.4f}"
                "<br>Geometry penalty: %{customdata[2]:.3f}<br>Protein rows: %{customdata[3]}"
                "<br>Finite RSCC / geometry: %{customdata[4]} / %{customdata[5]}<extra></extra>"
            ),
        ))
    fig.update_layout(
        template="plotly_white", height=650, uirevision="guidance-models",
        margin={"t": 100, "b": 60},
        legend={"orientation": "h", "y": 1.02, "yanchor": "bottom",
                "x": 0, "xanchor": "left", "groupclick": "togglegroup"},
        xaxis={"title": "Δ RSCC vs. baseline", "range": [-0.05, 0.25]},
        yaxis={"title": "Δ Geometry penalty vs. baseline (symlog)",
               "range": symlog([-1, 1000]).tolist(), "tickvals": symlog(YTICKS).tolist(),
               "ticktext": [str(v) for v in YTICKS]},
        clickmode="event+select",
    )
    return fig


def guidance_provenance(manifest):
    return {
        **manifest,
        "status": "Notebook 03 reconstruction; historical paper identity unverified",
        "filters": {"map_occupancy": {"A": 0.5, "B": 0.5}, "ensemble_size": 8,
                    "exclude_nan_rscc_proteins": False, "minimum_protein_rows": 3},
        "geometry_penalty": "mean(delta_clashscore, 100*delta_bond_length_fraction, "
                            "100*delta_bond_angle_fraction), skipping NaN terms as in notebook 03",
        "aggregation": "Marginal medians and 25th/75th percentiles across protein deltas",
        "display": "base-10 symlog, linthresh=1, linscale=1; raw values retained in tables/hover",
        "chosen_guidance_weights": OPTIMAL_GW,
        "initially_visible_models": [m for m in OPTIMAL_GW if m not in MODELS_ALT],
        "static_reconstruction_models": [m for m in OPTIMAL_GW if m not in MODELS_ALT],
        "export_scope": "CSV: all loaded models; reconstruction SVG: main models; "
                        "plot-toolbar SVG: current legend-visible view",
        "caveats": [
            "RSCC is averaged over selections per protein; geometry is averaged over members.",
            "Marginal RSCC/geometry summaries may have different finite contributors.",
            "RF3-Chiral-Off retains notebook 03's base-RSCC fallback.",
            "The baseline origin is a reference, not an additional measured point.",
        ],
        "analysis_sha256": {
            name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
            for name in ("guidance.py", "loaders.py", "transforms.py", "style.py",
                         "plot_interaction.py", "02_guidance.py", "uv.lock")
        },
    }
