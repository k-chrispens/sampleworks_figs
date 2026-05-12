"""Tradeoff and spiral figures (geometry penalty vs delta-RSCC trajectories)."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.cm import ScalarMappable
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D

from style import (
    GUIDANCE_MARKERS,
    GUIDANCE_PALETTE,
    GUIDANCE_PALETTE_NO_BASELINE,
    HIGHLIGHT_STYLES,
    LABEL_TEMPLATES,
    MODEL_DISPLAY,
    MODEL_PALETTE,
    OPTIMAL_GW,
)


def add_tradeoff_guidance_regions(
    ax,
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    guidance_order: list[str],
    guidance_palette: dict[str, str],
) -> None:
    """Overlay smooth guidance-specific KDE regions on a tradeoff axis."""
    for gs in guidance_order:
        subset = (
            df.loc[df["guidance_strength"] == gs, [x_col, y_col]]
            .replace([np.inf, -np.inf], np.nan)
            .dropna()
        )
        if (
            len(subset) < 8
            or subset[x_col].nunique() < 3
            or subset[y_col].nunique() < 3
        ):
            continue
        try:
            sns.kdeplot(
                data=subset, x=x_col, y=y_col, ax=ax,
                fill=True, levels=4, thresh=0.18,
                bw_adjust=1.05, cut=1.6,
                color=guidance_palette[gs], alpha=0.10,
                warn_singular=False, zorder=0,
            )
            sns.kdeplot(
                data=subset, x=x_col, y=y_col, ax=ax,
                fill=False, levels=2, thresh=0.24,
                bw_adjust=1.05, cut=1.6,
                color=guidance_palette[gs], alpha=0.35,
                linewidths=0.8, warn_singular=False, zorder=0.2,
            )
        except Exception:
            continue


def make_tradeoff_grid(
    tradeoff_df: pd.DataFrame,
    y_col: str,
    y_label: str,
    metric_rows: list[tuple[str, str]] | None = None,
    model_order: list[str] | None = None,
    guidance_order: list[str] | None = None,
    figsize: tuple[float, float] | None = None,
    scatter_size: int = 16,
    scatter_alpha: float = 0.24,
    show_centroids: bool = True,
    trajectory_protein: str | None = None,
    title: str | None = None,
):
    """Create a 2×N model subplot grid showing guidance tradeoff scatter plots."""
    if metric_rows is None:
        metric_rows = [
            ("delta_rscc", "Δ RSCC vs baseline"),
            ("delta_silhouette_to_ref", "Δ Clustering Score vs baseline"),
        ]
    if model_order is None:
        model_order = [
            "boltz2X-RAY", "boltz2MD", "protenix", "rf3",
            "protenix-mini", "protenix-tiny", "rf3-old", "rf3-chiral-off",
        ]
    _present = set(tradeoff_df["model_name"].unique())
    model_order = [m for m in model_order if m in _present]
    if guidance_order is None:
        guidance_order = ["weak", "medium", "strong"]
    if figsize is None:
        figsize = (4.25 * len(model_order), 7.6)

    y_vals = tradeoff_df[y_col].replace([np.inf, -np.inf], np.nan).dropna()
    y_pad = max(0.16 * (y_vals.max() - y_vals.min()), 0.30)
    y_lim = (float(y_vals.min() - y_pad), float(y_vals.max() + y_pad))

    row_xlim = {}
    for metric, _ in metric_rows:
        vals = tradeoff_df[metric].replace([np.inf, -np.inf], np.nan).dropna()
        if vals.empty:
            row_xlim[metric] = None
            continue
        max_abs = float(vals.abs().max())
        if max_abs == 0:
            max_abs = 0.01
        pad = max(0.12 * max_abs, 0.02)
        row_xlim[metric] = (-max_abs - pad, max_abs + pad)

    fig, axes = plt.subplots(
        len(metric_rows),
        len(model_order),
        figsize=figsize,
        sharex="row",
        sharey=True,
    )

    for row_idx, (metric, xlabel) in enumerate(metric_rows):
        for col_idx, model in enumerate(model_order):
            ax = axes[row_idx, col_idx]
            model_df = tradeoff_df[tradeoff_df["model_name"] == model]

            if model_df.empty:
                ax.axis("off")
                continue

            add_tradeoff_guidance_regions(
                ax, model_df, metric, y_col,
                guidance_order, GUIDANCE_PALETTE_NO_BASELINE,
            )

            for gs in guidance_order:
                gdf = model_df[model_df["guidance_strength"] == gs]
                if gdf.empty:
                    continue
                ax.scatter(
                    gdf[metric], gdf[y_col],
                    s=scatter_size, c=GUIDANCE_PALETTE_NO_BASELINE[gs],
                    alpha=scatter_alpha, marker=GUIDANCE_MARKERS[gs],
                    edgecolors="white", linewidths=0.2, zorder=1.2,
                )
                if show_centroids:
                    ax.scatter(
                        [gdf[metric].median()], [gdf[y_col].median()],
                        s=70, c=GUIDANCE_PALETTE_NO_BASELINE[gs],
                        marker=GUIDANCE_MARKERS[gs],
                        edgecolors="black", linewidths=0.8,
                        zorder=4,
                    )

            if trajectory_protein is not None:
                focus = model_df[model_df["protein"] == trajectory_protein].sort_values(
                    "guidance_order"
                )
                if not focus.empty:
                    ax.plot(
                        focus[metric], focus[y_col],
                        color="black", linewidth=1.0, alpha=0.75, zorder=5,
                    )
                    for row in focus.itertuples(index=False):
                        ax.scatter(
                            [getattr(row, metric)], [getattr(row, y_col)],
                            s=70, facecolors="white", edgecolors="black",
                            linewidths=1.0, marker=GUIDANCE_MARKERS.get(row.guidance_strength, "o"),
                            zorder=6,
                        )

            ax.axvline(0, color="#555555", linestyle="--", linewidth=1.0, alpha=0.6)
            ax.axhline(0, color="#555555", linestyle="--", linewidth=1.0, alpha=0.6)
            if row_xlim[metric] is not None:
                ax.set_xlim(row_xlim[metric])
            ax.set_ylim(y_lim)
            if row_idx == 0:
                ax.set_title(MODEL_DISPLAY.get(model, model))
            if col_idx == 0:
                ax.set_ylabel(y_label)
            ax.set_xlabel(xlabel)

    label_suffix = " (median)" if show_centroids else ""
    legend_handles = [
        Line2D(
            [0], [0],
            marker=GUIDANCE_MARKERS[gs], color="none",
            markerfacecolor=GUIDANCE_PALETTE_NO_BASELINE[gs],
            markeredgecolor="black", linestyle="none",
            markersize=8 if show_centroids else 7,
            label=f"{gs.capitalize()} guidance{label_suffix}",
        )
        for gs in guidance_order
    ]
    if trajectory_protein is not None:
        legend_handles.append(
            Line2D(
                [0], [0],
                color="black", linewidth=1.0,
                label=f"{trajectory_protein.upper()} trajectory",
            )
        )

    fig.legend(
        handles=legend_handles,
        loc="upper center",
        ncol=len(legend_handles),
        frameon=False,
        bbox_to_anchor=(0.5, 1.01),
    )
    if title is not None:
        fig.suptitle(title, y=1.05)
    fig.tight_layout()
    return fig, axes, model_order


def annotate_tradeoff_highlights(
    fig,
    axes,
    highlight_points: pd.DataFrame,
    y_col: str,
    metric_rows: list[tuple[str, str]] | None = None,
    model_order: list[str] | None = None,
    guidance_order: list[str] | None = None,
    title: str | None = None,
):
    """Add highlight ring markers, guidance dots, and smart labels to a tradeoff grid."""
    if metric_rows is None:
        metric_rows = [
            ("delta_rscc", "Δ RSCC vs baseline"),
            ("delta_silhouette_to_ref", "Δ Clustering Score vs baseline"),
        ]
    if model_order is None:
        model_order = [
            "boltz2X-RAY", "boltz2MD", "protenix", "rf3",
            "protenix-mini", "protenix-tiny", "rf3-old", "rf3-chiral-off",
        ]
    if guidance_order is None:
        guidance_order = ["weak", "medium", "strong"]

    for row_idx, (metric, _) in enumerate(metric_rows):
        for col_idx, model in enumerate(model_order):
            if col_idx >= axes.shape[-1]:
                break
            ax = axes[row_idx, col_idx]
            hdf = highlight_points[highlight_points["model_name"] == model].copy()
            if hdf.empty:
                continue
            hdf = hdf.sort_values([y_col, metric], ascending=[False, True])

            for flag, style in HIGHLIGHT_STYLES.items():
                subset = hdf[hdf[flag]]
                if subset.empty:
                    continue
                ax.scatter(
                    subset[metric], subset[y_col],
                    s=style["size"], facecolors="none",
                    edgecolors=style["color"], linewidths=2.0,
                    marker=style["marker"], zorder=4,
                )

            for gs in guidance_order:
                hgdf = hdf[hdf["guidance_strength"] == gs]
                if hgdf.empty:
                    continue
                ax.scatter(
                    hgdf[metric], hgdf[y_col],
                    s=42, c="black",
                    marker=GUIDANCE_MARKERS[gs], zorder=5,
                )

            if row_idx == 0:
                x_min, x_max = ax.get_xlim()
                y_min, y_max = ax.get_ylim()
                if x_max > x_min and y_max > y_min:
                    hdf["x_frac"] = (hdf[metric] - x_min) / (x_max - x_min)
                    hdf["y_frac"] = (hdf[y_col] - y_min) / (y_max - y_min)
                    hdf["label_side"] = hdf["x_frac"].map(
                        lambda v: (
                            "right" if v < 0.22
                            else "left" if v > 0.78
                            else "right" if v < 0.5
                            else "left"
                        )
                    )
                    hdf["label_vertical"] = hdf["y_frac"].map(
                        lambda v: "down" if v > 0.72 else "up"
                    )
                    for (side, vert), side_df in hdf.groupby(
                        ["label_side", "label_vertical"], sort=False
                    ):
                        side_df = side_df.sort_values(
                            [y_col, metric], ascending=[False, True]
                        ).reset_index(drop=True)
                        offsets = LABEL_TEMPLATES[(side, vert)]
                        for idx, row in enumerate(side_df.itertuples(index=False)):
                            offset = offsets[min(idx, len(offsets) - 1)]
                            ax.annotate(
                                row.point_label,
                                xy=(getattr(row, metric), getattr(row, y_col)),
                                xytext=offset, textcoords="offset points",
                                ha="left" if side == "right" else "right",
                                va="center",
                                bbox={
                                    "boxstyle": "round,pad=0.14",
                                    "fc": "white", "ec": "#666666",
                                    "alpha": 0.92, "lw": 0.65,
                                },
                                arrowprops={
                                    "arrowstyle": "-",
                                    "color": "#666666", "lw": 0.55,
                                    "shrinkA": 0, "shrinkB": 4,
                                },
                                zorder=6,
                            )

    for leg in fig.legends:
        leg.remove()
    legend_handles = [
        Line2D(
            [0], [0],
            marker=GUIDANCE_MARKERS[gs], color="none",
            markerfacecolor=GUIDANCE_PALETTE_NO_BASELINE[gs],
            markeredgecolor="black", linestyle="none",
            markersize=7,
            label=f"{gs.capitalize()} guidance",
        )
        for gs in guidance_order
    ]
    legend_handles.extend(
        [
            Line2D(
                [0], [0],
                marker=HIGHLIGHT_STYLES[f]["marker"], color="none",
                markerfacecolor="none",
                markeredgecolor=HIGHLIGHT_STYLES[f]["color"],
                markeredgewidth=2.0, markersize=8,
                label=HIGHLIGHT_STYLES[f]["label"],
            )
            for f in [
                "is_balance_best",
                "is_balance_worst",
                "is_peer_best",
                "is_peer_worst",
            ]
        ]
    )
    fig.legend(
        handles=legend_handles,
        loc="upper center",
        ncol=4,
        frameon=False,
        bbox_to_anchor=(0.5, 1.03),
    )
    if title is not None:
        fig.suptitle(title, y=1.08)
    fig.tight_layout()
    return fig


def compute_conceptual_tradeoff_limits(
    tradeoff_df: pd.DataFrame,
    guidance_strengths: tuple[str, ...] = ("weak", "strong"),
) -> tuple[tuple[float, float], tuple[float, float]]:
    """Compute shared RSCC/geometry axis limits for conceptual figures."""
    subset = (
        tradeoff_df[tradeoff_df["guidance_strength"].isin(guidance_strengths)][
            ["delta_rscc", "geometry_penalty"]
        ]
        .replace([float("inf"), float("-inf")], float("nan"))
        .dropna()
    )
    if subset.empty:
        raise ValueError("No finite rows available for conceptual limits.")

    x_min, x_max = float(subset["delta_rscc"].min()), float(subset["delta_rscc"].max())
    y_min, y_max = (
        float(subset["geometry_penalty"].min()),
        float(subset["geometry_penalty"].max()),
    )
    x_pad = max(0.08 * (x_max - x_min), 0.02)
    y_pad = max(0.10 * (y_max - y_min), 0.50)
    return (x_min - x_pad, x_max + x_pad), (y_min - y_pad, y_max + y_pad)


def make_conceptual_tradeoff_figure(
    tradeoff_df: pd.DataFrame,
    guidance_strength: str,
    *,
    xlim: tuple[float, float],
    ylim: tuple[float, float],
    title: str,
    zoom_quantiles: tuple[float, float] | None = None,
):
    """Create a simplified RSCC-versus-geometry tradeoff figure."""

    def _add_model_regions(ax, df, x_col, y_col, order, palette):
        for m in order:
            sub = (
                df.loc[df["model_name"] == m, [x_col, y_col]]
                .replace([float("inf"), float("-inf")], float("nan"))
                .dropna()
            )
            if len(sub) < 8 or sub[x_col].nunique() < 3 or sub[y_col].nunique() < 3:
                continue
            try:
                sns.kdeplot(
                    data=sub, x=x_col, y=y_col, ax=ax,
                    fill=True, levels=4, thresh=0.16,
                    bw_adjust=1.0, cut=1.6,
                    color=palette[m], alpha=0.09,
                    warn_singular=False, zorder=0,
                )
                sns.kdeplot(
                    data=sub, x=x_col, y=y_col, ax=ax,
                    fill=False, levels=2, thresh=0.22,
                    bw_adjust=1.0, cut=1.6,
                    color=palette[m], alpha=0.30,
                    linewidths=0.8, warn_singular=False, zorder=0.2,
                )
            except Exception:
                continue

    gdf = tradeoff_df[tradeoff_df["guidance_strength"] == guidance_strength].copy()
    if gdf.empty:
        raise ValueError(f"No rows found for guidance_strength={guidance_strength!r}.")

    plot_df = gdf.copy()
    plot_xlim, plot_ylim = xlim, ylim

    if zoom_quantiles is not None:
        lo, hi = zoom_quantiles
        x_lo, x_hi = gdf["delta_rscc"].quantile([lo, hi])
        y_lo, y_hi = gdf["geometry_penalty"].quantile([lo, hi])
        xp = max(0.10 * (x_hi - x_lo), 0.01)
        yp = max(0.16 * (y_hi - y_lo), 0.25)
        plot_xlim = (float(x_lo - xp), float(x_hi + xp))
        plot_ylim = (float(y_lo - yp), float(y_hi + yp))
        plot_df = gdf[
            gdf["delta_rscc"].between(*plot_xlim)
            & gdf["geometry_penalty"].between(*plot_ylim)
        ].copy()

    m_order = [
        m
        for m in [
            "boltz2X-RAY", "boltz2MD", "protenix", "rf3",
            "protenix-mini", "protenix-tiny", "rf3-old", "rf3-chiral-off",
        ]
        if m in plot_df["model_name"].unique()
    ] or sorted(plot_df["model_name"].unique())

    fig, ax = plt.subplots(figsize=(5.8, 4.2))
    _add_model_regions(
        ax, plot_df, "delta_rscc", "geometry_penalty", m_order, MODEL_PALETTE
    )

    for m in m_order:
        mdf = plot_df[plot_df["model_name"] == m]
        if mdf.empty:
            continue
        c = MODEL_PALETTE.get(m, "#4D4D4D")
        ax.scatter(
            mdf["delta_rscc"], mdf["geometry_penalty"],
            s=20, c=c, alpha=0.24,
            edgecolors="white", linewidths=0.25,
            zorder=1.2,
        )
        ax.scatter(
            [mdf["delta_rscc"].median()], [mdf["geometry_penalty"].median()],
            s=95, c=c,
            edgecolors="black", linewidths=0.9,
            zorder=4,
        )

    ax.axvline(0, color="#555555", linestyle="--", linewidth=1.0, alpha=0.6)
    ax.axhline(0, color="#555555", linestyle="--", linewidth=1.0, alpha=0.6)
    ax.set_xlim(plot_xlim)
    ax.set_ylim(plot_ylim)
    ax.set_xlabel("Δ RSCC vs baseline")
    ax.set_ylabel("Geometry penalty")
    ax.set_title(title)
    fig.tight_layout()
    return fig


def make_spiral_figure(
    model_df: pd.DataFrame,
    model_key: str,
    *,
    label_top_n: int = 5,
    xlim: tuple[float, float] | None = None,
    ylim: tuple[float, float] | None = None,
):
    """Create a single-model guidance spiral trajectory figure."""
    cmap = sns.color_palette("flare", as_cmap=True)
    fig, ax = plt.subplots(figsize=(7, 6))

    gw_values = sorted(model_df["guidance_weight"].unique())
    norm = Normalize(vmin=0.0, vmax=max(gw_values))
    endpoints = []

    for protein, pdf in model_df.groupby("protein"):
        pdf = pdf.sort_values("guidance_weight")
        xs = np.concatenate([[0.0], pdf["delta_rscc"].values])
        ys = np.concatenate([[0.0], pdf["geometry_penalty"].values])
        gws = np.concatenate([[0.0], pdf["guidance_weight"].values])

        arc = float(np.sqrt(np.diff(xs) ** 2 + np.diff(ys) ** 2).sum())
        endpoints.append((protein, xs[-1], ys[-1], arc))

        points = np.column_stack([xs, ys]).reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        seg_colors = cmap(norm(gws[1:]))
        seg_colors[:, 3] = 0.35
        ax.add_collection(
            LineCollection(segments, colors=seg_colors, linewidths=1.1, zorder=1)
        )

        dot_colors = cmap(norm(gws[1:]))
        dot_colors[:, 3] = 0.80
        ax.scatter(
            xs[1:], ys[1:],
            c=dot_colors, s=26,
            edgecolors="white", linewidths=0.3,
            zorder=2,
        )

        if len(xs) >= 2:
            ax.annotate(
                "",
                xy=(xs[-1], ys[-1]),
                xytext=(xs[-2], ys[-2]),
                arrowprops=dict(
                    arrowstyle="-|>",
                    color=cmap(norm(gws[-1])),
                    lw=0.0,
                    mutation_scale=8,
                ),
                zorder=3,
            )

    endpoints.sort(key=lambda t: t[3], reverse=True)
    for protein, xe, ye, _ in endpoints[:label_top_n]:
        ax.annotate(
            protein.upper(),
            xy=(xe, ye),
            xytext=(6, 6),
            textcoords="offset points",
            fontweight="bold",
            color="#333333",
            bbox=dict(
                boxstyle="round,pad=0.15",
                fc="white", ec="#aaaaaa",
                alpha=0.85, lw=0.5,
            ),
            zorder=6,
        )

    ax.axvline(0, color="#555555", ls="--", lw=0.8, alpha=0.5)
    ax.axhline(0, color="#555555", ls="--", lw=0.8, alpha=0.5)

    if xlim is not None and ylim is not None:
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
    else:
        mx = model_df["delta_rscc"].replace([np.inf, -np.inf], np.nan).dropna()
        my = model_df["geometry_penalty"].replace([np.inf, -np.inf], np.nan).dropna()
        x_lo, x_hi = min(0.0, float(mx.min())), max(0.0, float(mx.max()))
        y_lo, y_hi = min(0.0, float(my.min())), max(0.0, float(my.max()))
        x_pad = max(0.08 * (x_hi - x_lo), 0.02)
        y_pad = max(0.08 * (y_hi - y_lo), 0.50)
        ax.set_xlim(x_lo - x_pad, x_hi + x_pad)
        ax.set_ylim(y_lo - y_pad, y_hi + y_pad)

    ax.set_yscale("symlog", linthresh=1.0)
    ax.set_xlabel("Δ RSCC vs. baseline")
    ax.set_ylabel("∆ Geometry penalty vs. baseline")

    display = MODEL_DISPLAY.get(model_key, model_key)
    ax.set_title(
        f"{display}: guidance strength tradeoff", fontweight="bold"
    )

    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.75, pad=0.02)
    cbar.set_label("Guidance weight")

    fig.tight_layout()
    return fig


def make_aggregated_spiral_figure(
    tradeoff_df: pd.DataFrame,
    *,
    model_key: str | None = None,
    model_order: list[str] | None = None,
    xlim: tuple[float, float] | None = None,
    ylim: tuple[float, float] | None = None,
    linthresh: float = 1.0,
    min_n: int = 3,
    optimal_gw: dict[str, float] | None = None,
):
    """Aggregated guidance-weight spiral trajectory figure."""
    if optimal_gw is None:
        optimal_gw = OPTIMAL_GW

    agg = (
        tradeoff_df
        .groupby(["model_name", "guidance_weight"], as_index=False)
        .agg(
            x_med=("delta_rscc", "median"),
            x_q1=("delta_rscc", lambda s: s.quantile(0.25)),
            x_q3=("delta_rscc", lambda s: s.quantile(0.75)),
            y_med=("geometry_penalty", "median"),
            y_q1=("geometry_penalty", lambda s: s.quantile(0.25)),
            y_q3=("geometry_penalty", lambda s: s.quantile(0.75)),
            n=("protein", "count"),
        )
    )
    agg = agg[agg["n"] >= min_n].sort_values(["model_name", "guidance_weight"])

    single_mode = model_key is not None
    if model_key is not None:
        models_to_draw: list[str] = [model_key]
    else:
        if model_order is None:
            model_order = [
                "boltz2X-RAY", "boltz2MD", "protenix", "rf3",
                "protenix-mini", "protenix-tiny", "rf3-old", "rf3-chiral-off",
            ]
        models_to_draw = [m for m in model_order if m in agg["model_name"].unique()]

    fig, ax = plt.subplots(figsize=(8, 6) if not single_mode else (7, 6))

    flare_cmap = sns.color_palette("flare", as_cmap=True)
    if single_mode:
        gw_all = sorted(agg.loc[agg["model_name"] == model_key, "guidance_weight"].unique())
        gw_norm = Normalize(vmin=0.0, vmax=max(gw_all) if gw_all else 1.0)

    for m in models_to_draw:
        sub = agg[agg["model_name"] == m]
        if sub.empty:
            continue

        gws = sub["guidance_weight"].to_numpy()
        x_med = sub["x_med"].to_numpy()
        y_med = sub["y_med"].to_numpy()
        x_q1 = sub["x_q1"].to_numpy()
        x_q3 = sub["x_q3"].to_numpy()
        y_q1 = sub["y_q1"].to_numpy()
        y_q3 = sub["y_q3"].to_numpy()

        xerr = np.vstack([x_med - x_q1, x_q3 - x_med])
        yerr = np.vstack([y_med - y_q1, y_q3 - y_med])

        xs = np.concatenate([[0.0], x_med])
        ys = np.concatenate([[0.0], y_med])

        line_color = MODEL_PALETTE.get(m, "#4D4D4D")

        ax.plot(
            xs, ys,
            color=line_color, lw=1.6, alpha=0.85, zorder=3,
            label=MODEL_DISPLAY.get(m, m) if not single_mode else None,
        )

        ax.errorbar(
            x_med, y_med,
            xerr=xerr, yerr=yerr,
            fmt="none",
            ecolor=line_color, alpha=0.35,
            capsize=3, lw=1.0, zorder=2,
        )

        s_min, s_max = 30.0, 220.0
        gw_hi = float(max(gws)) if len(gws) else 1.0
        gw_lo = float(min(gws)) if len(gws) else 0.0
        if gw_hi > gw_lo:
            sizes = s_min + (gws - gw_lo) / (gw_hi - gw_lo) * (s_max - s_min)
        else:
            sizes = np.full_like(gws, (s_min + s_max) / 2.0, dtype=float)

        opt_gw = optimal_gw.get(m)
        if opt_gw is None:
            edge_colors = ["white"] * len(gws)
            edge_widths = [0.6 if single_mode else 0.5] * len(gws)
        else:
            edge_colors = [
                "#111111" if np.isclose(gw, opt_gw) else "white" for gw in gws
            ]
            edge_widths = [
                2.0 if np.isclose(gw, opt_gw) else (0.6 if single_mode else 0.5)
                for gw in gws
            ]

        if single_mode:
            dot_colors = flare_cmap(gw_norm(gws))
            ax.scatter(
                x_med, y_med,
                c=dot_colors, s=sizes,
                edgecolors=edge_colors, linewidths=edge_widths, zorder=4,
            )
        else:
            ax.scatter(
                x_med, y_med,
                c=[line_color] * len(gws), s=sizes,
                edgecolors=edge_colors, linewidths=edge_widths, zorder=4,
            )

        for gw, x, y in zip(gws, x_med, y_med):
            is_opt = opt_gw is not None and np.isclose(gw, opt_gw)
            ax.annotate(
                f"{gw:g}",
                xy=(x, y),
                xytext=(5, 5),
                textcoords="offset points",
                fontsize=7,
                fontweight="bold" if is_opt else "normal",
                color="#222222",
                bbox=dict(
                    boxstyle="round,pad=0.12",
                    fc="white", ec="none",
                    alpha=0.7,
                ),
                zorder=5,
            )

    ax.axvline(0, color="#555555", ls="--", lw=0.8, alpha=0.5)
    ax.axhline(0, color="#555555", ls="--", lw=0.8, alpha=0.5)

    ax.set_yscale("symlog", linthresh=linthresh)

    if xlim is None or ylim is None:
        agg_for_limits = (
            agg[agg["model_name"].isin(models_to_draw)] if models_to_draw else agg
        )
        all_x = np.concatenate([
            agg_for_limits["x_q1"].to_numpy(dtype=float),
            agg_for_limits["x_q3"].to_numpy(dtype=float),
            np.array([0.0]),
        ])
        all_y = np.concatenate([
            agg_for_limits["y_q1"].to_numpy(dtype=float),
            agg_for_limits["y_q3"].to_numpy(dtype=float),
            np.array([0.0]),
        ])
        all_x = all_x[np.isfinite(all_x)]
        all_y = all_y[np.isfinite(all_y)]
        if all_x.size == 0:
            all_x = np.array([0.0])
        if all_y.size == 0:
            all_y = np.array([0.0])
        x_lo, x_hi = float(all_x.min()), float(all_x.max())
        y_lo, y_hi = float(all_y.min()), float(all_y.max())
        x_pad = max(0.08 * (x_hi - x_lo), 0.02)
        y_pad = max(0.08 * (y_hi - y_lo), 0.50)

    ax.set_xlim(xlim if xlim is not None else (x_lo - x_pad, x_hi + x_pad))
    ax.set_ylim(ylim if ylim is not None else (y_lo - y_pad, y_hi + y_pad))

    ax.set_xlabel("Δ RSCC vs. baseline (median)")
    ax.set_ylabel("∆ Geometry penalty vs. baseline (median)")

    optimal_handle = Line2D(
        [0], [0],
        marker="o", color="w",
        markerfacecolor="#CCCCCC", markeredgecolor="#111111",
        markeredgewidth=2.0, markersize=10,
        label="Optimal guidance",
    )

    if model_key is not None:
        display = MODEL_DISPLAY.get(model_key, model_key)
        ax.set_title(
            f"{display}: aggregated guidance tradeoff",
            fontweight="bold",
        )
        sm = ScalarMappable(cmap=flare_cmap, norm=gw_norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, shrink=0.75, pad=0.02)
        cbar.set_label("Guidance weight")
        ax.legend(
            handles=[optimal_handle],
            loc="upper right",
            frameon=True, framealpha=0.9, edgecolor="none",
        )
    else:
        ax.set_title(
            "Aggregated guidance strength tradeoff",
            fontweight="bold",
        )
        model_handles, model_labels = ax.get_legend_handles_labels()
        ax.legend(
            handles=list(model_handles) + [optimal_handle],
            labels=list(model_labels) + ["Optimal guidance"],
            loc="lower right",
            frameon=True,
            framealpha=0.9,
            edgecolor="none",
            title="Model",
        )

    fig.tight_layout()
    return fig
