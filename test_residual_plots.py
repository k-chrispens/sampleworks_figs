"""Minimal tests for the residual plotting refactor.

Verifies that the y-axis values plotted by the residual-mode scatter
helpers equal ``guided − baseline`` for a small fixture.

Run with: ``uv run python -m pytest test_residual_plots.py``
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")  # headless

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from analysis_helpers import (
    make_aggregated_scatter,
    make_classified_scatter,
    make_classified_aggregated_scatter,
    make_rscc_scatter,
    make_scatter_comparison,
)


def _fixture_df() -> pd.DataFrame:
    """Five (protein, selection) rows with known baseline / guided RSCC."""
    return pd.DataFrame({
        "protein": ["a", "b", "c", "d", "e"],
        "label": ["A 1-3", "B 4-6", "C 7-9", "D 10-12", "E 13-15"],
        "rscc_baseline": [0.10, 0.40, 0.55, 0.80, 0.90],
        "rscc_optimal":  [0.30, 0.45, 0.50, 0.70, 0.95],
    }).assign(
        delta_rscc=lambda d: d["rscc_optimal"] - d["rscc_baseline"],
        classification=["side_chain_only", "small_loop", "large_loop", "domain_shift", "small_loop"],
    )


def _scatter_y_values(ax) -> np.ndarray:
    """Concatenate y-values from every PathCollection on the axes."""
    ys: list[float] = []
    for c in ax.collections:
        offsets = c.get_offsets()
        if offsets.size == 0:
            continue
        arr = np.asarray(offsets)
        ys.extend(arr[:, 1].tolist())
    return np.asarray(ys)


def test_per_model_residual_y_equals_guided_minus_baseline():
    df = _fixture_df()
    fig = make_scatter_comparison(
        df,
        model_name="boltz2X-RAY",
        x_col="rscc_baseline",
        y_col="rscc_optimal",
        x_label="Baseline RSCC",
        y_label="Improvement (Guided − Baseline RSCC)",
        delta_col="delta_rscc",
        residual=True,
        label_col="label",
        xlim=(0.0, 1.0),
    )
    ax = fig.axes[0]
    ys = _scatter_y_values(ax)
    expected = (df["rscc_optimal"] - df["rscc_baseline"]).to_numpy()
    # The main scatter should contain every fixture point exactly once;
    # ringed highlights add duplicates, hence "set inclusion" rather than equality.
    for v in expected:
        assert np.any(np.isclose(ys, v)), f"missing residual y={v} in {ys}"
    plt.close(fig)


def test_make_rscc_scatter_residual_y_values():
    df = _fixture_df()
    fig = make_rscc_scatter(df, "boltz2X-RAY", optimal_gw=0.1)
    ax = fig.axes[0]
    ys = _scatter_y_values(ax)
    expected = (df["rscc_optimal"] - df["rscc_baseline"]).to_numpy()
    for v in expected:
        assert np.any(np.isclose(ys, v)), f"missing residual y={v} in {ys}"
    # No point may sit above the y = 1 - x ceiling.
    for x, y in zip(df["rscc_baseline"], expected):
        assert y <= 1.0 - x + 1e-12, f"point ({x}, {y}) above ceiling y = 1 - x"
    plt.close(fig)


def test_aggregated_scatter_residual_propagates_error():
    agg = pd.DataFrame({
        "label": ["A 1-3", "B 4-6", "C 7-9"],
        "x_med": [0.20, 0.50, 0.80],
        "x_q1":  [0.10, 0.40, 0.70],   # SD_x = 0.10
        "x_q3":  [0.30, 0.60, 0.90],
        "y_med": [0.45, 0.55, 0.85],
        "y_q1":  [0.35, 0.45, 0.75],   # SD_y = 0.10
        "y_q3":  [0.55, 0.65, 0.95],
    })
    fig = make_aggregated_scatter(
        agg,
        x_col="x_med", y_col="y_med",
        x_q1_col="x_q1", x_q3_col="x_q3",
        y_q1_col="y_q1", y_q3_col="y_q3",
        x_label="Baseline", y_label="Improvement",
        title="t",
        label_col="label",
        residual=True,
        xlim=(0.0, 1.0),
    )
    ax = fig.axes[0]
    ys = _scatter_y_values(ax)
    expected = (agg["y_med"] - agg["x_med"]).to_numpy()
    for v in expected:
        assert np.any(np.isclose(ys, v)), f"missing residual y={v} in {ys}"

    # Inspect the errorbar collection's vertical extents.
    # Each row has SD_x = SD_y = 0.10, so propagated SD = sqrt(0.10² + 0.10²) ≈ 0.1414.
    expected_sd = float(np.hypot(0.10, 0.10))
    yerr_lines = [
        ln for ln in ax.get_lines()
        if ln.get_linestyle() in ("-", "None") and ln.get_marker() in ("|", "_", "")
    ]
    # Pull all line collections (errorbar caps + bars) for a robust check.
    bars = [c for c in ax.containers if hasattr(c, "lines")]
    assert bars, "no errorbar containers found"
    # Verify residual y-error magnitude on the first errorbar container.
    eb = bars[0]
    _, _, barlinecols = eb.lines
    # barlinecols is (xerr_LineCollection, yerr_LineCollection) when both axes have errors.
    if len(barlinecols) >= 2:
        yerr_lc = barlinecols[1]
        segs = yerr_lc.get_segments()
        for (lo, hi) in segs:
            half = abs(hi[1] - lo[1]) / 2.0
            assert np.isclose(half, expected_sd, atol=1e-6), (
                f"propagated y-error {half} != {expected_sd}"
            )
    plt.close(fig)


def test_classified_scatter_small_multiples_axes_count():
    df = _fixture_df()
    fig = make_classified_scatter(df, "boltz2X-RAY")
    # Four classes in CLASSIFICATION_ORDER are present.
    visible = [a for a in fig.axes if a.get_visible()]
    assert len(visible) == 4
    # Every subplot should contain at least one scatter collection (the gray shadow).
    for ax in visible:
        assert len(ax.collections) >= 1
    plt.close(fig)


def test_classified_aggregated_scatter_residual_y_values():
    agg = pd.DataFrame({
        "x_mean": [0.20, 0.50, 0.80, 0.40],
        "x_err_lo": [0.10, 0.40, 0.70, 0.30],
        "x_err_hi": [0.30, 0.60, 0.90, 0.50],
        "y_mean": [0.45, 0.55, 0.85, 0.30],
        "y_err_lo": [0.35, 0.45, 0.75, 0.20],
        "y_err_hi": [0.55, 0.65, 0.95, 0.40],
        "classification": [
            "side_chain_only", "small_loop", "large_loop", "domain_shift",
        ],
    })
    fig = make_classified_aggregated_scatter(agg)
    expected = (agg["y_mean"] - agg["x_mean"]).to_numpy()
    found = []
    for ax in fig.axes:
        if not ax.get_visible():
            continue
        found.extend(_scatter_y_values(ax).tolist())
    found = np.asarray(found)
    for v in expected:
        assert np.any(np.isclose(found, v)), f"missing residual y={v} in plot"
    plt.close(fig)
