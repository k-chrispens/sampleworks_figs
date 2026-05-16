"""Spiral and aggregated-spiral tradeoff figures."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.cm import ScalarMappable
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D

from style import MODEL_DISPLAY, MODEL_PALETTE, OPTIMAL_GW


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
                "protenix-mini", "protenix-tiny", "rf3-chiral-off",
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
                fontsize=10,
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
