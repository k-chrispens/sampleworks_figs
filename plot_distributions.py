"""Distribution-style plots (ridge / combined hist / baseline panels)."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from style import (
    GUIDANCE_PALETTE,
    MODEL_DISPLAY,
    MODEL_PALETTE,
    NICE_XLABELS,
)


def make_ridge_plot(
    df: pd.DataFrame,
    plot_column_name: str,
    custom_xlabel: str | None = None,
    smoothing: float = 0.15,
    is_delta: bool = False,
    xlim: tuple[float, float] | None = None,
    combined: bool = False,
    bin_width: float = 0.02,
    title: str | None = None,
):
    """Ridge/joy-plot comparing distributions across model types and guidance strengths.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with ``guidance_strength``, ``model_name``, and the plot column.
    plot_column_name : str
        Column to plot on the x-axis.
    smoothing : float
        KDE bandwidth adjustment.
    is_delta : bool
        If True, add a vertical reference line at x=0 and auto-center xlim.
    xlim : tuple[float, float] | None
        Manual x-axis limits.
    combined : bool
        If True, plot all models on a single axes as histplots with KDE,
        colored by model using MODEL_PALETTE.
    bin_width : float
        Bin width for histplot when ``combined=True``.
    title : str or None
        Optional suptitle for the figure.

    Returns
    -------
    sns.FacetGrid | plt.Figure
        The completed FacetGrid (default) or Figure (combined).
    """
    xlabel = (
        custom_xlabel
        if custom_xlabel
        else NICE_XLABELS.get(
            plot_column_name,
            plot_column_name.replace("_", " ").title(),
        )
    )

    if combined:
        return _make_combined_hist_plot(
            df,
            plot_column_name,
            xlabel,
            smoothing=smoothing,
            is_delta=is_delta,
            xlim=xlim,
            bin_width=bin_width,
            title=title,
        )

    hue_order = [
        k
        for k in ["baseline", "weak", "medium", "strong"]
        if k in df["guidance_strength"].unique()
    ]
    palette = [GUIDANCE_PALETTE[k] for k in hue_order]
    model_order = sorted(df["model_name"].unique())

    g = sns.FacetGrid(
        df,
        row="model_name",
        row_order=model_order,
        hue="guidance_strength",
        hue_order=hue_order,
        aspect=5,
        height=1.8,
        palette=palette,
        sharex=True,
        sharey=False,
    )
    g.map(
        sns.kdeplot,
        plot_column_name,
        bw_adjust=smoothing,
        clip_on=False,
        fill=True,
        alpha=0.75,
        linewidth=1.5,
    )
    g.map(
        sns.kdeplot,
        plot_column_name,
        clip_on=False,
        color="w",
        lw=2,
        bw_adjust=smoothing,
    )
    g.refline(y=0, linewidth=1.5, linestyle="-", color=None, clip_on=False)

    _xlim = xlim
    if is_delta:
        for ax in g.axes.flat:
            ax.axvline(
                x=0, color="#333333", lw=1.5, linestyle="--", alpha=0.6, zorder=3
            )
        if _xlim is None:
            _values = df[plot_column_name].replace([np.inf, -np.inf], np.nan).dropna()
            if not _values.empty:
                _max_abs = float(np.nanmax(np.abs(_values.to_numpy())))
                if _max_abs == 0:
                    _max_abs = 0.01
                _pad = max(0.06 * _max_abs, 0.005)
                _xlim = (-_max_abs - _pad, _max_abs + _pad)
    if _xlim is not None:
        g.set(xlim=_xlim)

    for ax in g.axes.flat:
        ax.spines["bottom"].set_visible(True)
        ax.tick_params(axis="x", which="both", bottom=True, labelbottom=True)
    g.despine(left=True)

    for ax, model in zip(g.axes.flat, model_order):
        display = MODEL_DISPLAY.get(model, model)
        ax.text(
            0.02,
            0.90,
            display,
            transform=ax.transAxes,
            fontweight="bold",
            va="top",
            ha="left",
            color="#222222",
            bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.8, ec="none"),
        )

    legend_handles = [
        Patch(facecolor=GUIDANCE_PALETTE[k], alpha=0.75, label=k.capitalize())
        for k in hue_order
    ]
    g.figure.legend(
        handles=legend_handles,
        loc="upper right",
        frameon=True,
        framealpha=0.9,
        edgecolor="none",
        title="Guidance",
        bbox_to_anchor=(0.99, 0.99),
    )

    g.set(xlabel=xlabel)
    _top = 0.93 if title else 1.0
    g.figure.subplots_adjust(hspace=0.6, right=0.95, top=_top)
    if title:
        g.figure.suptitle(title, fontweight="bold")
    g.set_titles("")
    g.set(yticks=[], yticklabels=[], ylabel="")
    return g


def _make_combined_hist_plot(
    df: pd.DataFrame,
    plot_column_name: str,
    xlabel: str,
    *,
    smoothing: float,
    is_delta: bool,
    xlim: tuple[float, float] | None,
    bin_width: float,
    title: str | None = None,
) -> plt.Figure:
    """Single-axes KDE with all models overlaid, colored by MODEL_PALETTE."""
    model_order = sorted(df["model_name"].unique())

    fig, ax = plt.subplots(figsize=(10, 3.2))

    for model in model_order:
        data = df[df["model_name"] == model]
        if data.empty:
            continue
        color = MODEL_PALETTE.get(model, "#888888")
        sns.kdeplot(
            data=data,
            x=plot_column_name,
            ax=ax,
            fill=True,
            alpha=0.35,
            linewidth=1.5,
            color=color,
            bw_adjust=smoothing,
        )

    if is_delta:
        ax.axvline(x=0, color="#333333", lw=1.5, linestyle="--", alpha=0.6, zorder=3)

    if xlim is not None:
        ax.set_xlim(xlim)

    ax.set_xlabel(xlabel)
    ax.set_ylabel("Density")
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

    legend_handles = [
        Patch(
            facecolor=MODEL_PALETTE.get(m, "#888888"),
            alpha=0.35,
            label=MODEL_DISPLAY.get(m, m),
        )
        for m in model_order
    ]
    ax.legend(
        handles=legend_handles,
        loc="upper right",
        frameon=True,
        framealpha=0.9,
        edgecolor="none",
        title="Model",
    )

    if title:
        fig.suptitle(title, fontweight="bold")
        fig.tight_layout(rect=[0, 0, 1, 0.95])
    else:
        fig.tight_layout()
    return fig


def make_combined_baseline_distributions(
    rscc_df: pd.DataFrame,
    lddt_df: pd.DataFrame,
    rscc_col: str = "rscc",
    lddt_col: str = "avg_silhouette_to_ref",
    rscc_xlabel: str = "RSCC",
    lddt_xlabel: str = "Clustering Score",
    smoothing: float = 0.15,
    bin_width: float = 0.02,
    cumulative: bool = False,
    cumulative_style: str = "kde",
) -> plt.Figure:
    """Two baseline distribution panels (RSCC on top, LDDT below) sharing one legend.

    When ``cumulative`` is True, panels switch to a cumulative view; ``cumulative_style``
    selects ``"kde"`` (smoothed cumulative KDE, respects ``smoothing``) or ``"ecdf"``
    (empirical CDF step plot).
    """
    if cumulative_style not in ("kde", "ecdf"):
        raise ValueError(f"cumulative_style must be 'kde' or 'ecdf', got {cumulative_style!r}")

    model_order = sorted(
        set(rscc_df["model_name"].unique()) | set(lddt_df["model_name"].unique())
    )

    fig, (ax_rscc, ax_lddt) = plt.subplots(
        2, 1, figsize=(10, 5.6), sharex=False,
    )

    for ax, df, col, xlabel in [
        (ax_rscc, rscc_df, rscc_col, rscc_xlabel),
        (ax_lddt, lddt_df, lddt_col, lddt_xlabel),
    ]:
        for model in model_order:
            data = df[df["model_name"] == model]
            if data.empty:
                continue
            color = MODEL_PALETTE.get(model, "#888888")
            if not cumulative:
                sns.kdeplot(
                    data=data, x=col, ax=ax,
                    fill=True, alpha=0.35, linewidth=1.5,
                    color=color, bw_adjust=smoothing,
                )
            elif cumulative_style == "kde":
                sns.kdeplot(
                    data=data, x=col, ax=ax,
                    cumulative=True, fill=False, linewidth=1.8,
                    color=color, bw_adjust=smoothing,
                )
            else:
                sns.ecdfplot(
                    data=data, x=col, ax=ax,
                    color=color, linewidth=1.8,
                )

        ax.set_xlabel(xlabel)
        if not cumulative:
            ax.set_ylabel("Density")
        elif cumulative_style == "kde":
            ax.set_ylabel("Cumulative Density")
        else:
            ax.set_ylabel("Proportion")
        for spine in ["top", "right"]:
            ax.spines[spine].set_visible(False)

    if cumulative:
        legend_handles = [
            Line2D(
                [0], [0],
                color=MODEL_PALETTE.get(m, "#888888"),
                linewidth=1.8,
                label=MODEL_DISPLAY.get(m, m),
            )
            for m in model_order
        ]
    else:
        legend_handles = [
            Patch(
                facecolor=MODEL_PALETTE.get(m, "#888888"),
                alpha=0.35,
                label=MODEL_DISPLAY.get(m, m),
            )
            for m in model_order
        ]
    ax_rscc.legend(
        handles=legend_handles,
        loc="upper right",
        frameon=True, framealpha=0.9, edgecolor="none",
        title="Model",
    )

    fig.tight_layout()
    return fig
