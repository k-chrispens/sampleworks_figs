"""Distribution-style plots (combined baseline RSCC + LDDT panel)."""

from __future__ import annotations

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from style import MODEL_DISPLAY, MODEL_PALETTE


def make_combined_baseline_distributions(
    rscc_df: pd.DataFrame,
    lddt_df: pd.DataFrame,
    rscc_col: str = "rscc",
    lddt_col: str = "clustering_score",
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
