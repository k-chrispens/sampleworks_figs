"""Classification-aware scatter plots (side-chain / loop / domain-shift classes).

Phase 1 preserves the existing three classified scatter helpers; Phase 3
replaces them with a single combined-grid helper.
"""

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


def make_classified_scatter(
    df: pd.DataFrame,
    model_name: str,
    *,
    x_col: str = "rscc_baseline",
    y_col: str = "rscc_optimal",
    class_col: str = "classification",
    x_label: str = "Baseline RSCC",
    y_label: str | None = None,
    title: str | None = None,
    draw_diagonal: bool = True,  # ignored
    show_correlation: bool = False,
    equal_aspect: bool = True,  # ignored
) -> plt.Figure:
    """Per-class small-multiples residual plot for a single model."""
    display_name = MODEL_DISPLAY.get(model_name, model_name)
    if title is None:
        title = f"{display_name}: RSCC Improvement by Conformation Type"
    if y_label is None:
        y_label = "Improvement (Guided − Baseline RSCC)"

    classes = [c for c in CLASSIFICATION_ORDER if c in set(df[class_col])]
    n_classes = max(1, len(classes))
    n_cols = min(4, n_classes)
    n_rows = int(np.ceil(n_classes / n_cols))
    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(3.4 * n_cols, 3.6 * n_rows),
        sharex=True, sharey=True,
        squeeze=False,
    )
    flat_axes = axes.flatten()

    x_all = df[x_col].to_numpy()
    y_all_resid = (df[y_col] - df[x_col]).to_numpy()

    x_lo, x_hi = 0.0, 1.0
    y_hi = 1.0
    _finite = np.isfinite(y_all_resid)
    if _finite.any():
        y_lo = float(min(0.0, np.nanmin(y_all_resid[_finite]))) - 0.02
    else:
        y_lo = -0.05
    _xs = np.linspace(x_lo, x_hi, 100)

    for ax_i, cls in zip(flat_axes, classes):
        mask = df[class_col] == cls
        ax_i.scatter(
            x_all, y_all_resid,
            color="lightgray", alpha=0.25,
            s=14, edgecolors="none", zorder=1,
        )
        sub_x = df.loc[mask, x_col].to_numpy()
        sub_y = (df.loc[mask, y_col] - df.loc[mask, x_col]).to_numpy()
        color = CLASSIFICATION_PALETTE[cls]
        ax_i.scatter(
            sub_x, sub_y,
            c=color, alpha=0.7,
            s=18, edgecolors="none", zorder=3,
        )
        ax_i.axhline(0.0, ls="--", color="#AAAAAA", lw=1, zorder=2)
        ax_i.plot(_xs, 1.0 - _xs, ls="--", color="#AAAAAA", lw=0.8, zorder=2)
        ax_i.set_title(
            f"{CLASSIFICATION_DISPLAY[cls]} (n={int(mask.sum())})", pad=10,
        )
        ax_i.set_xlabel(x_label)

    for ax_i in flat_axes[len(classes):]:
        ax_i.set_visible(False)

    for row_idx in range(n_rows):
        flat_axes[row_idx * n_cols].set_ylabel(y_label)

    flat_axes[0].set_xlim(x_lo, x_hi)
    flat_axes[0].set_ylim(y_lo, y_hi)

    fig.suptitle(title, y=0.99)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.92))
    return fig


def make_classified_scatter_grid(
    df: pd.DataFrame,
    *,
    model_col: str = "model_name",
    x_col: str = "rscc_baseline",
    y_col: str = "rscc_optimal",
    class_col: str = "classification",
    x_label: str = "Baseline RSCC",
    y_label: str = "Improvement (Guided − Baseline RSCC)",
    title: str = "RSCC Improvement by Conformation Type",
    model_order: list[str] | None = None,
) -> plt.Figure:
    """Combined per-class residual scatter for many models on one figure."""
    classes = [c for c in CLASSIFICATION_ORDER if c in set(df[class_col])]
    n_cols = max(1, len(classes))

    if model_order is not None:
        models = [m for m in model_order if m in set(df[model_col])]
    else:
        models = sorted(df[model_col].unique())
    n_rows = max(1, len(models))

    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(3.8 * n_cols, 3.4 * n_rows),
        sharex=True, sharey=True,
        squeeze=False,
    )

    y_all_resid = (df[y_col] - df[x_col]).to_numpy()
    x_lo, x_hi = 0.0, 1.0
    y_hi = 1.0
    _finite = np.isfinite(y_all_resid)
    if _finite.any():
        y_lo = float(min(0.0, np.nanmin(y_all_resid[_finite]))) - 0.02
    else:
        y_lo = -0.05
    _xs = np.linspace(x_lo, x_hi, 100)

    _unique_pairs = df.drop_duplicates(["protein", "selection"])
    class_counts = (
        _unique_pairs[class_col].value_counts().to_dict()
    )

    for row_idx, model in enumerate(models):
        sub_df = df[df[model_col] == model]
        x_sub = sub_df[x_col].to_numpy()
        y_sub = (sub_df[y_col] - sub_df[x_col]).to_numpy()

        for col_idx, cls in enumerate(classes):
            ax = axes[row_idx, col_idx]
            ax.scatter(
                x_sub, y_sub,
                color="lightgray", alpha=0.25,
                s=14, edgecolors="none", zorder=1,
            )
            mask = sub_df[class_col] == cls
            cx = sub_df.loc[mask, x_col].to_numpy()
            cy = (sub_df.loc[mask, y_col] - sub_df.loc[mask, x_col]).to_numpy()
            color = CLASSIFICATION_PALETTE[cls]
            ax.scatter(
                cx, cy,
                c=color, alpha=0.7,
                s=18, edgecolors="none", zorder=3,
            )
            ax.axhline(0.0, ls="--", color="#AAAAAA", lw=1, zorder=2)
            ax.plot(_xs, 1.0 - _xs, ls="--", color="#AAAAAA", lw=0.8, zorder=2)
            if row_idx == 0:
                ax.set_title(
                    f"{CLASSIFICATION_DISPLAY[cls]} "
                    f"(n={int(class_counts.get(cls, 0))})",
                    pad=10,
                )

    axes[0, 0].set_xlim(x_lo, x_hi)
    axes[0, 0].set_ylim(y_lo, y_hi)

    fig.supxlabel(x_label, y=0.02)
    fig.supylabel(y_label, x=0.01)
    fig.suptitle(title, y=0.985)
    fig.subplots_adjust(
        left=0.07, right=0.99,
        top=0.93, bottom=0.07,
        hspace=0.45, wspace=0.08,
    )

    for row_idx, model in enumerate(models):
        left_box = axes[row_idx, 0].get_position()
        right_box = axes[row_idx, -1].get_position()
        x_center = 0.5 * (left_box.x0 + right_box.x1)
        y_text = left_box.y1 + 0.005
        fig.text(
            x_center, y_text,
            MODEL_DISPLAY.get(model, model),
            ha="center", va="bottom",
            fontsize=13, weight="bold",
        )

    return fig


def make_classified_aggregated_scatter(
    agg_df: pd.DataFrame,
    *,
    x_col: str = "x_mean",
    y_col: str = "y_mean",
    x_err_lo_col: str | None = "x_err_lo",
    x_err_hi_col: str | None = "x_err_hi",
    y_err_lo_col: str | None = "y_err_lo",
    y_err_hi_col: str | None = "y_err_hi",
    class_col: str = "classification",
    x_label: str = "Baseline RSCC (mean ± SD)",
    y_label: str = "Improvement (Guided − Baseline RSCC)",
    title: str = "RSCC Improvement by Conformation Type",
    draw_diagonal: bool = True,  # ignored
    equal_aspect: bool = True,  # ignored
    show_correlation: bool = False,
    show_errors: bool = True,
    label_col: str | None = None,
) -> plt.Figure:
    """Per-class small-multiples residual plot with propagated error bars."""
    classes = [c for c in CLASSIFICATION_ORDER if c in set(agg_df[class_col])]
    n_classes = max(1, len(classes))
    n_cols = min(4, n_classes)
    n_rows = int(np.ceil(n_classes / n_cols))
    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(3.4 * n_cols, 3.6 * n_rows),
        sharex=True, sharey=True,
        squeeze=False,
    )
    flat_axes = axes.flatten()

    x_all_med = agg_df[x_col].to_numpy()
    y_all_resid = (agg_df[y_col] - agg_df[x_col]).to_numpy()

    x_lo, x_hi = 0.0, 1.0
    y_hi = 1.0
    _finite_y = np.isfinite(y_all_resid)
    if _finite_y.any():
        y_lo = float(min(0.0, np.nanmin(y_all_resid[_finite_y]))) - 0.02
    else:
        y_lo = -0.05
    _xs = np.linspace(x_lo, x_hi, 100)

    for ax_i, cls in zip(flat_axes, classes):
        mask = agg_df[class_col] == cls
        ax_i.scatter(
            x_all_med, y_all_resid,
            color="lightgray", alpha=0.25,
            s=14, edgecolors="none", zorder=1,
        )

        sub = agg_df[mask]
        color = CLASSIFICATION_PALETTE[cls]
        if not sub.empty:
            x_med = sub[x_col].to_numpy()
            y_med = (sub[y_col] - sub[x_col]).to_numpy()
            _err_cols_present = all(
                c is not None for c in (
                    x_err_lo_col, x_err_hi_col, y_err_lo_col, y_err_hi_col,
                )
            )
            draw_errors = show_errors and _err_cols_present
            if draw_errors:
                sd_x = (sub[x_err_hi_col].to_numpy() - sub[x_err_lo_col].to_numpy()) / 2.0
                sd_y = (sub[y_err_hi_col].to_numpy() - sub[y_err_lo_col].to_numpy()) / 2.0
                sd_resid = np.sqrt(sd_x**2 + sd_y**2)
                xerr = np.vstack([
                    np.clip(x_med - sub[x_err_lo_col].to_numpy(), 0, None),
                    np.clip(sub[x_err_hi_col].to_numpy() - x_med, 0, None),
                ])
                yerr = np.vstack([sd_resid, sd_resid])
                ax_i.errorbar(
                    x_med, y_med,
                    xerr=xerr, yerr=yerr,
                    fmt="none", ecolor=color, alpha=0.35, capsize=3, lw=1.0, zorder=2,
                )
            ax_i.scatter(
                x_med, y_med,
                c=color, alpha=0.7,
                s=22, edgecolors="none", zorder=4,
            )

        ax_i.axhline(0.0, ls="--", color="#AAAAAA", lw=1, zorder=2)
        ax_i.plot(_xs, 1.0 - _xs, ls="--", color="#AAAAAA", lw=0.8, zorder=2)
        ax_i.set_title(
            f"{CLASSIFICATION_DISPLAY[cls]} (n={int(mask.sum())})", pad=10,
        )
        ax_i.set_xlabel(x_label)

    for ax_i in flat_axes[len(classes):]:
        ax_i.set_visible(False)

    for row_idx in range(n_rows):
        flat_axes[row_idx * n_cols].set_ylabel(y_label)

    flat_axes[0].set_xlim(x_lo, x_hi)
    flat_axes[0].set_ylim(y_lo, y_hi)

    fig.suptitle(title, y=0.99)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.92))
    return fig
