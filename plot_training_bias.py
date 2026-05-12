"""Training-distribution-bias plots (cluster vs predicted occupancy)."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import to_rgb
from matplotlib.legend_handler import HandlerTuple
from matplotlib.lines import Line2D

from style import MODEL_DISPLAY, MODEL_PALETTE

_TRAINING_BIAS_POINT_COLOR = "#333333"
_TRAINING_BIAS_MAIN_MODELS = ["boltz2X-RAY", "boltz2MD", "protenix", "rf3"]


def _training_bias_prepare_scatter_frame(df: pd.DataFrame, condition_label: str) -> pd.DataFrame:
    out = df.copy()
    n_a = out["n_assigned_to_A"]
    n_b = out["n_assigned_to_B"]
    s = (n_a + n_b).astype("float")
    out["s_cluster"] = s
    out["frac_b_cluster"] = np.where(s > 0, n_b / s.replace(0, np.nan), np.nan)
    out["condition"] = condition_label
    out = out.dropna(subset=["frac_b_cluster", "b_occupancy"])
    out = out[out["s_cluster"] > 0]
    return out


def make_training_bias_scatter(
    baseline_df: pd.DataFrame,
    optimal_df: pd.DataFrame,
    *,
    model_order: list[str] | None = None,
    jitter: float = 0.01,
    seed: int = 0,
) -> plt.Figure:
    if model_order is None:
        model_order = list(_TRAINING_BIAS_MAIN_MODELS)

    baseline_prepared = _training_bias_prepare_scatter_frame(baseline_df, "Baseline")
    optimal_prepared = _training_bias_prepare_scatter_frame(optimal_df, "Optimal guidance")
    combined = pd.concat([baseline_prepared, optimal_prepared], ignore_index=True)
    s_max = float(combined["s_cluster"].max()) if not combined.empty else 1.0
    log_s_max = np.log(s_max + 1.0) if s_max > 0 else 1.0

    fig, axes = plt.subplots(
        2, 4,
        sharex=True, sharey=True,
        figsize=(11, 6),
        constrained_layout=True,
    )

    row_label = {"Baseline": 0, "Optimal guidance": 1}
    rng = np.random.default_rng(seed)

    for condition, row in row_label.items():
        for col, model_name in enumerate(model_order):
            ax = axes[row, col]
            ax.axhline(0.5, ls=":", color="#CCCCCC", lw=0.8, zorder=1)

            sub = combined[
                (combined["condition"] == condition)
                & (combined["model_name"] == model_name)
            ]
            n_points = len(sub)
            if n_points > 0:
                x = sub["frac_b_cluster"].to_numpy() + rng.uniform(-jitter, jitter, n_points)
                y = sub["b_occupancy"].to_numpy()
                alpha = np.log(sub["s_cluster"].to_numpy() + 1.0) / log_s_max
                rgba = np.zeros((n_points, 4))
                rgba[:, :3] = to_rgb(_TRAINING_BIAS_POINT_COLOR)
                rgba[:, 3] = np.clip(alpha, 0.0, 1.0)
                ax.scatter(x, y, s=18, c=rgba, linewidths=0, zorder=3)

            ax.set_xlim(-0.02, 1.02)
            ax.set_ylim(-0.02, 1.02)
            ax.set_xticks([0.0, 0.5, 1.0])
            ax.set_yticks([0.0, 0.5, 1.0])
            if row == 0:
                ax.set_title(MODEL_DISPLAY.get(model_name, model_name), fontsize=11)
            ax.tick_params(labelbottom=(row == 1), labelleft=(col == 0))

    for row, label in [(0, "Baseline"), (1, "Optimal guidance")]:
        axes[row, 0].annotate(
            label,
            xy=(-0.32, 0.5), xycoords="axes fraction",
            ha="center", va="center", rotation=90, fontsize=11,
        )

    fig.supxlabel("fraction in state B (training cluster)")
    fig.supylabel("fraction in state B (predicted)")
    return fig


def _training_bias_compute_alpha(s_values: np.ndarray, s_max: float) -> np.ndarray:
    """Half-strength log-scaled per-point alpha used by the calibration figures."""
    if s_max <= 0:
        return np.zeros_like(s_values, dtype=float)
    return np.clip(0.5 * np.log(s_values + 1.0) / np.log(s_max + 1.0), 0.0, 1.0)


def _training_bias_binned_means(
    x: np.ndarray,
    y: np.ndarray,
    *,
    n_bins: int = 10,
    min_count: int = 5,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per-bin mean and sample std of ``y`` in ``n_bins`` equal-count bins of ``x``."""
    if len(x) == 0:
        return np.array([]), np.array([]), np.array([])
    quantiles = np.linspace(0.0, 1.0, n_bins + 1)
    edges = np.unique(np.quantile(x, quantiles))
    n_actual_bins = len(edges) - 1
    if n_actual_bins <= 0:
        return np.array([]), np.array([]), np.array([])
    centers = np.full(n_actual_bins, np.nan)
    means = np.full(n_actual_bins, np.nan)
    stds = np.full(n_actual_bins, np.nan)
    for i in range(n_actual_bins):
        if i == n_actual_bins - 1:
            mask = (x >= edges[i]) & (x <= edges[i + 1])
        else:
            mask = (x >= edges[i]) & (x < edges[i + 1])
        if int(mask.sum()) >= min_count:
            centers[i] = float(x[mask].mean())
            means[i] = float(y[mask].mean())
            stds[i] = float(np.std(y[mask], ddof=1))
    return centers, means, stds


def _training_bias_slope_r(
    x: np.ndarray,
    y: np.ndarray,
) -> tuple[float, float] | None:
    """OLS slope of ``y`` on ``x`` and Pearson r, or ``None`` when undefined."""
    if len(x) < 2 or float(np.std(x)) == 0.0:
        return None
    slope = float(np.polyfit(x, y, 1)[0])
    r = float(np.corrcoef(x, y)[0, 1])
    return slope, r


def _training_bias_draw_reference_lines(ax: plt.Axes) -> None:
    """Draw y=x (training prior) and y=0.5 (guidance target) reference lines."""
    ax.axhline(0.5, ls="--", color="#CCCCCC", lw=1, zorder=0)
    ax.plot(
        [-0.02, 1.02], [-0.02, 1.02],
        ls="--", color="#AAAAAA", lw=1, zorder=1,
    )


def make_training_bias_calibration_overlay(
    baseline_df: pd.DataFrame,
    optimal_df: pd.DataFrame,
    *,
    model_order: list[str] | None = None,
    jitter: float = 0.01,
    seed: int = 0,
    max_cluster_size: int | None = None,
    cutoff_status: str | None = None,
    n_bins: int = 10,
    min_bin_count: int = 5,
) -> plt.Figure:
    """1x4 calibration overlay of Baseline and Optimal guidance per model."""
    if model_order is None:
        model_order = list(_TRAINING_BIAS_MAIN_MODELS)

    baseline_prepared = _training_bias_prepare_scatter_frame(baseline_df, "Baseline")
    optimal_prepared = _training_bias_prepare_scatter_frame(optimal_df, "Optimal guidance")

    if max_cluster_size is not None:
        baseline_prepared = baseline_prepared[baseline_prepared["n_cluster_evaluated"] <= max_cluster_size]
        optimal_prepared = optimal_prepared[optimal_prepared["n_cluster_evaluated"] <= max_cluster_size]
    if cutoff_status is not None:
        baseline_prepared = baseline_prepared[baseline_prepared["cutoff_status"] == cutoff_status]
        optimal_prepared = optimal_prepared[optimal_prepared["cutoff_status"] == cutoff_status]

    combined = pd.concat([baseline_prepared, optimal_prepared], ignore_index=True)
    s_max = float(combined["s_cluster"].max()) if not combined.empty else 1.0

    fig = plt.figure(figsize=(14, 4.2))
    outer = fig.add_gridspec(
        1, 4, wspace=0.18,
        left=0.06, right=0.995, top=0.92, bottom=0.18,
    )
    rng = np.random.default_rng(seed)

    base_color = "#888888"
    hist_bins = np.linspace(0.0, 1.0, 21)

    anchor_sc: plt.Axes | None = None
    for col, model_name in enumerate(model_order):
        guidance_color = MODEL_PALETTE[model_name]
        inner = outer[0, col].subgridspec(1, 2, width_ratios=[4, 1], wspace=0.08)
        if anchor_sc is None:
            ax_sc = fig.add_subplot(inner[0, 0])
            anchor_sc = ax_sc
        else:
            ax_sc = fig.add_subplot(inner[0, 0], sharex=anchor_sc, sharey=anchor_sc)
        ax_mg = fig.add_subplot(inner[0, 1], sharey=ax_sc)

        _training_bias_draw_reference_lines(ax_sc)

        baseline_panel = baseline_prepared[baseline_prepared["model_name"] == model_name]
        optimal_panel = optimal_prepared[optimal_prepared["model_name"] == model_name]

        for panel_df, color in [(baseline_panel, base_color), (optimal_panel, guidance_color)]:
            n_points = len(panel_df)
            if n_points == 0:
                continue
            x = panel_df["frac_b_cluster"].to_numpy() + rng.uniform(-jitter, jitter, n_points)
            y = panel_df["b_occupancy"].to_numpy()
            alpha = _training_bias_compute_alpha(panel_df["s_cluster"].to_numpy(), s_max)
            rgba = np.zeros((n_points, 4))
            rgba[:, :3] = to_rgb(color)
            rgba[:, 3] = alpha
            ax_sc.scatter(x, y, s=12, c=rgba, linewidths=0, zorder=3)

        for panel_df, color, lw, alpha_curve, ms, line_zorder, band_alpha, band_zorder in [
            (baseline_panel, base_color,     1.5, 0.8,  4, 4, 0.08, 3.4),
            (optimal_panel,  guidance_color, 2.5, 0.95, 5, 5, 0.15, 3.6),
        ]:
            if panel_df.empty:
                continue
            centers, means, stds = _training_bias_binned_means(
                panel_df["frac_b_cluster"].to_numpy(),
                panel_df["b_occupancy"].to_numpy(),
                n_bins=n_bins, min_count=min_bin_count,
            )
            ax_sc.fill_between(
                centers, means - stds, means + stds,
                color=color, alpha=band_alpha, linewidth=0, zorder=band_zorder,
            )
            ax_sc.plot(
                centers, means,
                color=color, lw=lw, alpha=alpha_curve,
                marker="o", markersize=ms, zorder=line_zorder,
            )

        y_text = 0.97
        for panel_df, text_color in [
            (baseline_panel, base_color),
            (optimal_panel, guidance_color),
        ]:
            slope_r = _training_bias_slope_r(
                panel_df["frac_b_cluster"].to_numpy(),
                panel_df["b_occupancy"].to_numpy(),
            )
            if slope_r is None:
                continue
            slope, r_val = slope_r
            ax_sc.text(
                0.03, y_text, f"slope = {slope:.2f}\nr = {r_val:.2f}",
                transform=ax_sc.transAxes, ha="left", va="top",
                fontsize=9, color=text_color, zorder=6,
            )
            y_text -= 0.24

        ax_sc.set_xlim(-0.02, 1.02)
        ax_sc.set_ylim(-0.02, 1.02)
        ax_sc.set_xticks([0.0, 0.5, 1.0])
        ax_sc.set_yticks([0.0, 0.5, 1.0])
        ax_sc.set_box_aspect(1)
        ax_sc.set_title(MODEL_DISPLAY.get(model_name, model_name), fontsize=11)
        ax_sc.tick_params(labelleft=(col == 0))

        if not baseline_panel.empty:
            ax_mg.hist(
                baseline_panel["b_occupancy"].to_numpy(),
                bins=hist_bins, orientation="horizontal",
                color=base_color, alpha=0.45,
            )
        if not optimal_panel.empty:
            ax_mg.hist(
                optimal_panel["b_occupancy"].to_numpy(),
                bins=hist_bins, orientation="horizontal",
                color=guidance_color, alpha=0.65,
            )
        ax_mg.set_box_aspect(4)
        plt.setp(ax_mg.get_yticklabels(), visible=False)
        ax_mg.tick_params(
            axis="x", which="both", bottom=False, labelbottom=False,
        )

    fig.supxlabel("fraction in state B (training cluster)", y=0.04)
    fig.supylabel("fraction in state B (predicted)", x=0.005)
    return fig


def make_training_bias_no_cluster_strip(
    baseline_df: pd.DataFrame,
    optimal_df: pd.DataFrame,
    *,
    model_order: list[str] | None = None,
    jitter: float = 0.35,
    seed: int = 0,
) -> plt.Figure:
    """1x4 strip plot for selections with no PDB neighbors in the training cluster."""
    if model_order is None:
        model_order = list(_TRAINING_BIAS_MAIN_MODELS)

    def _filter_no_cluster(df: pd.DataFrame) -> pd.DataFrame:
        sub = df.dropna(subset=["b_occupancy"]).copy()
        return sub[
            (sub["n_assigned_to_A"] == 0)
            & (sub["n_assigned_to_B"] == 0)
        ]

    baseline_nc = _filter_no_cluster(baseline_df)
    optimal_nc = _filter_no_cluster(optimal_df)

    fig = plt.figure(figsize=(14, 4.2))
    outer = fig.add_gridspec(
        1, 4, wspace=0.18,
        left=0.06, right=0.995, top=0.92, bottom=0.18,
    )
    rng = np.random.default_rng(seed)

    base_color = "#888888"
    hist_bins = np.linspace(0.0, 1.0, 21)

    anchor_sc: plt.Axes | None = None
    for col, model_name in enumerate(model_order):
        guidance_color = MODEL_PALETTE[model_name]
        inner = outer[0, col].subgridspec(1, 2, width_ratios=[4, 1], wspace=0.08)
        if anchor_sc is None:
            ax_sc = fig.add_subplot(inner[0, 0])
            anchor_sc = ax_sc
        else:
            ax_sc = fig.add_subplot(inner[0, 0], sharex=anchor_sc, sharey=anchor_sc)
        ax_mg = fig.add_subplot(inner[0, 1], sharey=ax_sc)

        ax_sc.axhline(0.5, ls="--", color="#CCCCCC", lw=1, zorder=1)

        baseline_panel = baseline_nc[baseline_nc["model_name"] == model_name]
        optimal_panel = optimal_nc[optimal_nc["model_name"] == model_name]

        for panel_df, color in [(baseline_panel, base_color), (optimal_panel, guidance_color)]:
            n_points = len(panel_df)
            if n_points == 0:
                continue
            x = rng.uniform(-jitter, jitter, n_points)
            y = panel_df["b_occupancy"].to_numpy()
            ax_sc.scatter(x, y, s=12, c=color, alpha=0.5, linewidths=0, zorder=3)

        ax_sc.text(
            0.03, 0.97,
            f"n = {len(optimal_panel)}",
            transform=ax_sc.transAxes, ha="left", va="top",
            fontsize=9, color="#333333", zorder=6,
        )

        ax_sc.set_xlim(-0.5, 0.5)
        ax_sc.set_ylim(-0.02, 1.02)
        ax_sc.set_xticks([])
        ax_sc.set_yticks([0.0, 0.5, 1.0])
        ax_sc.set_box_aspect(1)
        ax_sc.set_title(MODEL_DISPLAY.get(model_name, model_name), fontsize=11)
        ax_sc.tick_params(labelleft=(col == 0))

        if not baseline_panel.empty:
            ax_mg.hist(
                baseline_panel["b_occupancy"].to_numpy(),
                bins=hist_bins, orientation="horizontal",
                color=base_color, alpha=0.45,
            )
        if not optimal_panel.empty:
            ax_mg.hist(
                optimal_panel["b_occupancy"].to_numpy(),
                bins=hist_bins, orientation="horizontal",
                color=guidance_color, alpha=0.65,
            )
        ax_mg.set_box_aspect(4)
        plt.setp(ax_mg.get_yticklabels(), visible=False)
        ax_mg.tick_params(
            axis="x", which="both", bottom=False, labelbottom=False,
        )

    fig.supylabel("fraction in state B (predicted)", x=0.005)
    return fig


def make_training_bias_calibration_overlay_legend(
    *,
    model_order: list[str] | None = None,
) -> plt.Figure:
    """Standalone legend figure for ``make_training_bias_calibration_overlay``."""
    if model_order is None:
        model_order = list(_TRAINING_BIAS_MAIN_MODELS)

    fig, ax = plt.subplots(figsize=(2.5, 0.9), constrained_layout=True)
    ax.axis("off")

    baseline_handle = Line2D(
        [], [], color="#888888", lw=1.5, alpha=0.8,
        marker="o", markersize=4,
    )
    optimal_handles = tuple(
        Line2D(
            [], [], color=MODEL_PALETTE[m], lw=2.5, alpha=0.95,
            marker="o", markersize=5,
        )
        for m in model_order
    )

    ax.legend(
        handles=[baseline_handle, optimal_handles],
        labels=["Baseline", "Optimal guidance"],
        handler_map={tuple: HandlerTuple(ndivide=None, pad=0.5)},
        loc="center",
        ncol=1,
        frameon=False,
        fontsize=9,
    )
    return fig


def make_training_bias_small_multiples_v2(
    baseline_df: pd.DataFrame,
    optimal_df: pd.DataFrame,
    *,
    model_order: list[str] | None = None,
    jitter: float = 0.01,
    seed: int = 0,
) -> plt.Figure:
    """2x4 calibrated small multiples with per-panel mean curves and marginals."""
    if model_order is None:
        model_order = list(_TRAINING_BIAS_MAIN_MODELS)

    baseline_prepared = _training_bias_prepare_scatter_frame(baseline_df, "Baseline")
    optimal_prepared = _training_bias_prepare_scatter_frame(optimal_df, "Optimal guidance")
    combined = pd.concat([baseline_prepared, optimal_prepared], ignore_index=True)
    s_max = float(combined["s_cluster"].max()) if not combined.empty else 1.0

    fig = plt.figure(figsize=(13, 7), constrained_layout=True)
    outer = fig.add_gridspec(2, 4, wspace=0.25, hspace=0.18)
    rng = np.random.default_rng(seed)
    hist_bins = np.linspace(0.0, 1.0, 21)

    rows = [
        (0, "Baseline", baseline_prepared),
        (1, "Optimal guidance", optimal_prepared),
    ]
    axes_scatter: list[list[plt.Axes | None]] = [[None] * 4, [None] * 4]

    for row, _condition_label, prepared_frame in rows:
        for col, model_name in enumerate(model_order):
            inner = outer[row, col].subgridspec(1, 2, width_ratios=[4, 1], wspace=0.05)
            anchor = axes_scatter[0][0]
            if anchor is None:
                ax_sc = fig.add_subplot(inner[0, 0])
            else:
                ax_sc = fig.add_subplot(inner[0, 0], sharex=anchor, sharey=anchor)
            ax_mg = fig.add_subplot(inner[0, 1], sharey=ax_sc)
            axes_scatter[row][col] = ax_sc

            _training_bias_draw_reference_lines(ax_sc)

            panel_color = "#888888" if row == 0 else MODEL_PALETTE[model_name]

            panel_df = prepared_frame[prepared_frame["model_name"] == model_name]
            n_points = len(panel_df)
            if n_points > 0:
                x = panel_df["frac_b_cluster"].to_numpy() + rng.uniform(-jitter, jitter, n_points)
                y = panel_df["b_occupancy"].to_numpy()
                alpha = _training_bias_compute_alpha(panel_df["s_cluster"].to_numpy(), s_max)
                rgba = np.zeros((n_points, 4))
                rgba[:, :3] = to_rgb(panel_color)
                rgba[:, 3] = alpha
                ax_sc.scatter(x, y, s=12, c=rgba, linewidths=0, zorder=3)

                centers, means, stds = _training_bias_binned_means(
                    panel_df["frac_b_cluster"].to_numpy(),
                    panel_df["b_occupancy"].to_numpy(),
                )
                ax_sc.fill_between(
                    centers, means - stds, means + stds,
                    color=panel_color, alpha=0.15, linewidth=0, zorder=4,
                )
                ax_sc.plot(
                    centers, means,
                    color=panel_color, lw=2.0,
                    marker="o", markersize=4, zorder=5,
                )

                slope_r = _training_bias_slope_r(
                    panel_df["frac_b_cluster"].to_numpy(),
                    panel_df["b_occupancy"].to_numpy(),
                )
                if slope_r is not None:
                    slope, r_val = slope_r
                    ax_sc.text(
                        0.03, 0.97, f"slope = {slope:.2f}\nr = {r_val:.2f}",
                        transform=ax_sc.transAxes, ha="left", va="top",
                        fontsize=9, zorder=6,
                    )

            ax_sc.set_xlim(-0.02, 1.02)
            ax_sc.set_ylim(-0.02, 1.02)
            ax_sc.set_xticks([0.0, 0.5, 1.0])
            ax_sc.set_yticks([0.0, 0.5, 1.0])
            if row == 0:
                ax_sc.set_title(MODEL_DISPLAY.get(model_name, model_name), fontsize=11)
            ax_sc.tick_params(labelbottom=(row == 1), labelleft=(col == 0))

            if not panel_df.empty:
                ax_mg.hist(
                    panel_df["b_occupancy"].to_numpy(),
                    bins=hist_bins, orientation="horizontal",
                    color=panel_color, alpha=0.6,
                )
            plt.setp(ax_mg.get_yticklabels(), visible=False)
            ax_mg.tick_params(
                axis="x", which="both", bottom=False, labelbottom=False,
            )
            ax_mg.spines["bottom"].set_visible(False)

    for row, label in [(0, "Baseline"), (1, "Optimal guidance")]:
        anchor_ax = axes_scatter[row][0]
        if anchor_ax is not None:
            anchor_ax.annotate(
                label,
                xy=(-0.32, 0.5), xycoords="axes fraction",
                ha="center", va="center", rotation=90, fontsize=11,
            )

    fig.supxlabel("fraction in state B (training cluster)")
    fig.supylabel("fraction in state B (predicted)")
    return fig


def make_training_bias_calibration_explainer() -> plt.Figure:
    """Single empty panel with the y=x and y=0.5 references and their bisector."""
    fig, ax = plt.subplots(figsize=(4, 4), constrained_layout=True)

    ax.plot(
        [-0.02, 1.02], [-0.02, 1.02],
        ls="--", color="#AAAAAA", lw=1, zorder=1,
    )
    ax.plot(
        [-0.02, 1.02], [0.5, 0.5],
        ls="--", color="#AAAAAA", lw=1, zorder=1,
    )
    ax.plot(
        [-0.02, 1.02], [0.5 * -0.02 + 0.25, 0.5 * 1.02 + 0.25],
        ls="-", color="#222222", lw=2.0, zorder=2,
    )

    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xticks([0.0, 0.5, 1.0])
    ax.set_yticks([0.0, 0.5, 1.0])
    ax.set_xlabel("fraction in state B (training cluster)")
    ax.set_ylabel("fraction in state B (predicted)")
    return fig
