"""Combined 4 models x 4 classes by-class RSCC residual scatter."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from style import (
    CLASSIFICATION_DISPLAY,
    CLASSIFICATION_ORDER,
    CLASSIFICATION_PALETTE,
    MODEL_DISPLAY,
)

_DEFAULT_MAIN_MODELS = ["boltz2X-RAY", "boltz2MD", "protenix", "rf3"]


def _compute_global_residual_range(df: pd.DataFrame) -> tuple[float, float]:
    """Return ``(y_lo, y_hi)`` for residual y across the full pooled frame.

    ``y_lo = min(0, global_min_residual) - 0.02``; ``y_hi = 1.0``.
    """
    resid = (df["rscc_optimal"] - df["rscc_baseline"]).to_numpy()
    finite = np.isfinite(resid)
    if finite.any():
        y_lo = float(min(0.0, np.nanmin(resid[finite]))) - 0.02
    else:
        y_lo = -0.05
    return y_lo, 1.0


def _draw_class_panel(
    ax: plt.Axes,
    x_bg: np.ndarray, y_bg: np.ndarray,
    x_pts: np.ndarray, y_pts: np.ndarray,
    color: str,
    x_ref_line: np.ndarray,
) -> None:
    """Draw the per-panel shadow + class points + reference lines."""
    # Reference geometry drawn first (behind points).
    ax.axhline(0.0, ls="--", color="#AAAAAA", lw=1, zorder=1)
    ax.plot(x_ref_line, 1.0 - x_ref_line, ls="--", color="#AAAAAA", lw=0.8, zorder=1)
    # Light-gray shadow: that row's model pooled points.
    ax.scatter(
        x_bg, y_bg,
        color="lightgray", alpha=0.25,
        s=14, edgecolors="none", zorder=2,
    )
    # Class points on top.
    ax.scatter(
        x_pts, y_pts,
        c=color, alpha=0.7,
        s=18, edgecolors="none", zorder=3,
    )


def make_by_class_combined_grid(
    merged_df: pd.DataFrame,
    *,
    model_order: list[str] | None = None,
) -> plt.Figure:
    """4 main models × 4 classes by-class residual scatter.

    Parameters
    ----------
    merged_df
        Tidy frame with columns ``protein, selection, model_name,
        rscc_baseline, rscc_optimal, classification``. Rows for the
        domain-shift class use whole-domain RSCC; small-segment classes
        use per-altloc-region RSCC. Caller is responsible for the concat.
    model_order
        Optional explicit row ordering of models. Defaults to the four
        main models in display order, restricted to those present.

    Returns
    -------
    matplotlib.figure.Figure
    """
    if model_order is None:
        model_order = list(_DEFAULT_MAIN_MODELS)
    models = [m for m in model_order if m in set(merged_df["model_name"].unique())]
    classes = [c for c in CLASSIFICATION_ORDER if c in set(merged_df["classification"])]
    n_rows, n_cols = max(1, len(models)), max(1, len(classes))

    # Global y range across the full pooled dataset.
    y_lo, y_hi = _compute_global_residual_range(merged_df)
    x_lo, x_hi = 0.0, 1.0
    x_ref = np.linspace(x_lo, x_hi, 100)

    # n per class = unique (protein, selection) count across all rows in this
    # class — invariant across model rows.
    class_counts = (
        merged_df.drop_duplicates(["protein", "selection"])["classification"]
        .value_counts()
        .to_dict()
    )

    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(3.8 * n_cols, 3.4 * n_rows),
        sharex=True, sharey=True,
        squeeze=False,
    )

    for row_idx, model in enumerate(models):
        sub_df = merged_df[merged_df["model_name"] == model]
        x_bg = sub_df["rscc_baseline"].to_numpy()
        y_bg = (sub_df["rscc_optimal"] - sub_df["rscc_baseline"]).to_numpy()

        for col_idx, cls in enumerate(classes):
            ax = axes[row_idx, col_idx]
            mask = sub_df["classification"] == cls
            x_pts = sub_df.loc[mask, "rscc_baseline"].to_numpy()
            y_pts = (
                sub_df.loc[mask, "rscc_optimal"]
                - sub_df.loc[mask, "rscc_baseline"]
            ).to_numpy()
            _draw_class_panel(
                ax, x_bg, y_bg, x_pts, y_pts,
                color=CLASSIFICATION_PALETTE[cls],
                x_ref_line=x_ref,
            )
            if row_idx == 0:
                ax.set_title(
                    f"{CLASSIFICATION_DISPLAY[cls]} "
                    f"(n={int(class_counts.get(cls, 0))})",
                    pad=10,
                )

    axes[0, 0].set_xlim(x_lo, x_hi)
    axes[0, 0].set_ylim(y_lo, y_hi)

    fig.supxlabel("Baseline RSCC", y=0.02)
    fig.supylabel("Improvement (Guided − Baseline RSCC)", x=0.01)

    # Wider hspace than the legacy grid so the row banner sits in clear
    # whitespace between rows without colliding with the next row's class
    # titles. Narrower wspace keeps columns tight.
    fig.subplots_adjust(
        left=0.09, right=0.99,
        top=0.92, bottom=0.08,
        hspace=0.55, wspace=0.08,
    )

    # Horizontal model banner above each row.
    for row_idx, model in enumerate(models):
        left_box = axes[row_idx, 0].get_position()
        right_box = axes[row_idx, -1].get_position()
        x_center = 0.5 * (left_box.x0 + right_box.x1)
        y_text = left_box.y1 + 0.015
        fig.text(
            x_center, y_text,
            MODEL_DISPLAY.get(model, model),
            ha="center", va="bottom",
            fontsize=13, weight="bold",
        )

    return fig
