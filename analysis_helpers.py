"""Shared helpers for the initial-40 analysis notebooks.

Data loading, tagging, plotting utilities, and constants used across
01_geometry, 02_dihedrals, 03_fit_quality, and 04_tradeoffs notebooks.
"""

from __future__ import annotations

import os
import re
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import polars as pl
import seaborn as sns
from matplotlib.cm import ScalarMappable
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

# ── Constants ──────────────────────────────────────────────────────────────────

OCC_FILTER = "0.5occA_0.5occB"
OCC_FILTER_ALTLOC = "{'A': 0.5, 'B': 0.5}"

# PyMOL "good" colors (from ~/.pymolrc, cred. Sam Pellock / RFdiffusion)
COLOR_GOOD_GRAY = "#DCDCDC"
COLOR_GOOD_DARKBLUE = "#4B5FAA"
COLOR_GOOD_TEAL = "#4FB9AF"
COLOR_GOOD_NAVAHO = "#FFE0AC"
COLOR_GOOD_PINK = "#FFACB7"
COLOR_GOOD_LIGHTBLUE = "#8EDCE6"

# Model groups for split plotting
MODELS_ALT = {"protenix-mini", "protenix-tiny", "rf3-old", "rf3-chiral-off"}
MODELS_PROTENIX_GROUP = ["protenix", "protenix-mini", "protenix-tiny"]
MODELS_RF3_GROUP = ["rf3", "rf3-old", "rf3-chiral-off"]

GUIDANCE_PALETTE = {
    "baseline": "#888888",
    "weak": "#377EB8",
    "medium": "#FF7F00",
    "strong": "#E41A1C",
}

GUIDANCE_PALETTE_NO_BASELINE = {
    k: v for k, v in GUIDANCE_PALETTE.items() if k != "baseline"
}

MODEL_PALETTE = {
    "boltz2X-RAY": COLOR_GOOD_DARKBLUE,
    "boltz2MD": COLOR_GOOD_TEAL,
    "protenix": COLOR_GOOD_PINK,
    "rf3": "#FFD700",
    "ground_truth": COLOR_GOOD_GRAY,
    "protenix-mini": "#FF6B9D",
    "protenix-tiny": "#CC3377",
    "rf3-old": COLOR_GOOD_NAVAHO,
    "rf3-chiral-off": "#DAA520",
}

GUIDANCE_MARKERS = {
    "weak": "o",
    "medium": "s",
    "strong": "^",
}

HIGHLIGHT_STYLES = {
    "is_balance_best": {
        "color": "#1B9E77",
        "marker": "o",
        "label": "Best balance",
        "size": 180,
    },
    "is_balance_worst": {
        "color": "#CC3366",
        "marker": "o",
        "label": "Worst balance",
        "size": 180,
    },
    "is_peer_best": {
        "color": "#1F78B4",
        "marker": "D",
        "label": "Differentially better than peers",
        "size": 135,
    },
    "is_peer_worst": {
        "color": "#6A3D9A",
        "marker": "D",
        "label": "Differentially worse than peers",
        "size": 135,
    },
}

LABEL_TEMPLATES = {
    ("right", "up"): [(8, 10), (8, 22), (8, -2), (8, 34)],
    ("right", "down"): [(8, -8), (8, -20), (8, 4), (8, -32)],
    ("left", "up"): [(-10, 10), (-10, 22), (-10, -2), (-10, 34)],
    ("left", "down"): [(-10, -8), (-10, -20), (-10, 4), (-10, -32)],
}

MODEL_DISPLAY = {
    "boltz2_X-RAY_DIFFRACTION": "Boltz2-XRD",
    "boltz2_MD": "Boltz2-MD",
    "rf3": "RF3",
    "rf3-old": "RF3-Old",
    "protenix": "Protenix",
    "boltz2X-RAY_DIFFRACTION": "Boltz2-XRD",
    "boltz2X-RAY": "Boltz2-XRD",
    "boltz2xray": "Boltz2-XRD",
    "boltz2MD": "Boltz2-MD",
    "boltz2md": "Boltz2-MD",
    "ground_truth": "Ground truth",
    "protenix-mini": "Protenix-Mini",
    "protenix-tiny": "Protenix-Tiny",
    "rf3-chiral-off": "RF3-Chiral-Off",
}

MODEL_NAME_ALIAS = {
    "boltz2_X-RAY_DIFFRACTION": "boltz2X-RAY",
    "boltz2_MD": "boltz2MD",
    "boltz2xray": "boltz2X-RAY",
    "boltz2md": "boltz2MD",
    "rf3": "rf3-old",
    "rf3-fixed": "rf3",
}

NICE_XLABELS = {
    "delta_clashscore": "Δ Clash Score",
    "delta_outlier_fraction_bond_angles": "Δ Bond Angle Outlier Fraction",
    "delta_outlier_fraction_bond_lengths": "Δ Bond Length Outlier Fraction",
    "delta_rscc": "Δ RSCC",
    "delta_silhouette_to_ref": "Δ Avg. Silhouette to Ref",
    "clashscore": "Clash Score",
    "rscc": "Real Space Correlation Coefficient",
    "avg_silhouette_to_ref": "Avg. Silhouette to Ref",
    "ramachandran_z_score": "Ramachandran Z-Score",
    "torsion_z_score": "Torsion Z-Score",
    "outlier_fraction_bond_angles": "Bond Angle Outlier Fraction",
    "outlier_fraction_bond_lengths": "Bond Length Outlier Fraction",
}

GW_TO_STRENGTH = {
    0.0: "baseline",
    0.005: "weak",
    0.01: "weak",
    0.02: "medium",
    0.035: "medium",
    0.05: "strong",
    0.1: "weak",
    0.2: "medium",
    0.35: "medium",
    0.5: "strong",
}

FIGURE_DIR = Path("figures")

_PATH_PATTERN = re.compile(
    r"([0-9A-Za-z]{4})_([\d.]+occ[A-Z](?:_[\d.]+occ[A-Z])*)/(.*?)/(.*?)/ens8_gw([0-9\.]+)/"
)


# ── Data loading and tagging ──────────────────────────────────────────────────


def extract_metadata_from_path(p: str) -> pd.Series:
    """Extract pdb_id, occ_spec, model_name, steering, guidance_weight from a path string.

    Parameters
    ----------
    p : str
        File path containing the standard occ-sweep directory structure.

    Returns
    -------
    pd.Series
        Series with keys: pdb_id, occ_spec, model_name, steering, guidance_weight.
    """
    m = _PATH_PATTERN.search(p)
    return pd.Series(
        {
            "pdb_id": m.group(1),
            "occ_spec": m.group(2),
            "model_name": m.group(3),
            "steering": m.group(4),
            "guidance_weight": m.group(5),
        }
    )


def add_simplified_metadata_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Add guidance_strength bucket and optional selection_name to a DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with at least a ``guidance_weight`` column.

    Returns
    -------
    pd.DataFrame
        Same DataFrame with ``guidance_strength`` (and optionally ``selection_name``) added.
    """
    if "selection" in df.columns:
        df["selection_name"] = df.selection.apply(
            lambda x: "all altloc" if "chain_id" in str(x) else x
        )
    df["guidance_weight"] = df.guidance_weight.astype(float)
    df["guidance_strength"] = df.guidance_weight.map(GW_TO_STRENGTH)
    return df


def tag_path_dataframe(
    df: pd.DataFrame,
    filename_col: str = "filename",
    extra_renames: dict[str, str] | None = None,
    extra_cols: dict[str, str] | None = None,
) -> pd.DataFrame:
    """Apply path-based metadata tagging to an already-loaded DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        Raw DataFrame to tag.
    filename_col : str
        Column containing file paths; renamed to ``"filename"``.
    extra_renames : dict[str, str] | None
        Additional column renames (e.g. ``{"model": "model_n"}``).
    extra_cols : dict[str, str] | None
        Columns to add before ``add_simplified_metadata_columns``
        (e.g. ``{"selection": "all residues"}``).

    Returns
    -------
    pd.DataFrame
        Tagged DataFrame with ``protein``, ``guidance_strength``, etc.
    """
    renames = {filename_col: "filename"}
    if extra_renames:
        renames.update(extra_renames)
    df = df.rename(columns=renames)
    df = pd.merge(
        df,
        df.filename.apply(extract_metadata_from_path),
        left_index=True,
        right_index=True,
    )
    _occ = OCC_FILTER
    df = df.query("occ_spec == @_occ").copy()
    if extra_cols:
        for col, val in extra_cols.items():
            df[col] = val
    df = add_simplified_metadata_columns(df)
    df["protein"] = df.pdb_id.str.lower()
    df["model_name"] = df["model_name"].map(lambda x: MODEL_NAME_ALIAS.get(x, x))
    return df


def load_and_tag_path_csv(
    csv_path: str,
    filename_col: str = "filename",
    extra_renames: dict[str, str] | None = None,
    extra_cols: dict[str, str] | None = None,
) -> pd.DataFrame:
    """Load a CSV and apply path-based metadata tagging.

    Parameters
    ----------
    csv_path : str
        Path to the CSV file.
    filename_col : str
        Column containing file paths; renamed to ``"filename"``.
    extra_renames : dict[str, str] | None
        Additional column renames.
    extra_cols : dict[str, str] | None
        Columns to add before ``add_simplified_metadata_columns``.

    Returns
    -------
    pd.DataFrame
        Tagged DataFrame.
    """
    df = pd.read_csv(csv_path)
    return tag_path_dataframe(df, filename_col, extra_renames, extra_cols)


def load_and_tag_path_csv_fast(
    csv_path: str,
    filename_col: str = "path",
    model_col: str = "model",
) -> pd.DataFrame:
    """Load a large path-tagged CSV using Polars for ~10-20x faster I/O.

    Drop-in replacement for ``load_and_tag_path_csv`` on million-row files.
    Uses Polars' multithreaded CSV reader and vectorized regex extraction
    instead of per-row ``pandas.Series.apply``.

    Parameters
    ----------
    csv_path : str
        Path to the CSV file (e.g. per-residue tortoize results).
    filename_col : str
        Column containing file paths (renamed to ``"filename"``).
    model_col : str
        Column containing the model/ensemble index (renamed to ``"model_n"``).

    Returns
    -------
    pd.DataFrame
        Tagged DataFrame with the same columns as ``load_and_tag_path_csv``.
    """

    return (
        pl.scan_csv(csv_path)
        .rename({filename_col: "filename", model_col: "model_n"})
        .with_columns(
            pl.col("filename")
            .str.extract_groups(
                r"(?P<pdb_id>[0-9A-Za-z]{4})"
                r"_(?P<occ_spec>[\d.]+occ[A-Z](?:_[\d.]+occ[A-Z])*)"
                r"/(?P<model_name>.*?)/(?P<steering>.*?)"
                r"/ens8_gw(?P<guidance_weight>[0-9\.]+)/"
            )
            .alias("_meta"),
        )
        .unnest("_meta")
        .filter(pl.col("occ_spec") == OCC_FILTER)
        .with_columns(
            pl.col("guidance_weight").cast(pl.Float64),
            pl.col("pdb_id").str.to_lowercase().alias("protein"),
        )
        .with_columns(
            pl.col("guidance_weight")
            .replace_strict(GW_TO_STRENGTH, default=None, return_dtype=pl.String)
            .alias("guidance_strength"),
        )
        .with_columns(
            pl.col("model_name")
            .replace_strict(MODEL_NAME_ALIAS, default=None, return_dtype=pl.String)
            .fill_null(pl.col("model_name"))
            .alias("model_name"),
        )
        .collect()
        .to_pandas()
    )


def load_rscc(csv_path: str) -> pd.DataFrame:
    """Load RSCC results with altloc filtering and model-name construction.

    Parameters
    ----------
    csv_path : str
        Path to the RSCC CSV file.

    Returns
    -------
    pd.DataFrame
        RSCC DataFrame with ``model_name``, ``guidance_strength``, etc.
    """
    df = pd.read_csv(csv_path)
    _occ = OCC_FILTER_ALTLOC
    df = df.query("altloc_occupancies == @_occ").copy()
    df["model_name"] = df.apply(
        lambda x: x["model"] + x["method"] if x["model"] == "boltz2" else x["model"],
        axis=1,
    )
    df["model_name"] = df["model_name"].map(lambda x: MODEL_NAME_ALIAS.get(x, x))
    df = df[
        [
            "protein",
            "model_name",
            "ensemble_size",
            "guidance_weight",
            "rscc",
            "selection",
        ]
    ]
    df = add_simplified_metadata_columns(df)
    return df


def load_lddt(csv_path: str) -> pd.DataFrame:
    """Load LDDT results with graceful missing-file fallback.

    Parameters
    ----------
    csv_path : str
        Path to the LDDT CSV file.

    Returns
    -------
    pd.DataFrame
        LDDT DataFrame (may be empty if file not found).
    """
    _cols = [
        "protein",
        "model_name",
        "guidance_weight",
        "selection",
        "occupancies",
        "avg_silhouette",
        "avg_silhouette_to_ref",
    ]
    if not os.path.exists(csv_path):
        print(f"⚠️  {csv_path} not found — LDDT sections will be empty")
        return pd.DataFrame(columns=_cols)
    df = pd.read_csv(csv_path)
    df["model_name"] = df.apply(
        lambda x: x["model"] + x["method"] if x["model"] == "boltz2" else x["model"],
        axis=1,
    )
    df["model_name"] = df["model_name"].map(lambda x: MODEL_NAME_ALIAS.get(x, x))
    return df[_cols]


# ── Plotting helpers ──────────────────────────────────────────────────────────


def save_fig(plot_obj, filename: str, fmt: str = "svg") -> Path:
    """Save a matplotlib figure or FacetGrid to the figures directory.

    Parameters
    ----------
    plot_obj : matplotlib.figure.Figure or seaborn.FacetGrid
        Object to save.
    filename : str
        Output filename (relative to ``FIGURE_DIR``).  The extension is
        replaced by *fmt* so existing ``.png`` names are handled automatically.
    fmt : str
        Output format (default ``"svg"``).

    Returns
    -------
    Path
        Path to the saved file.
    """
    figure = plot_obj.figure if hasattr(plot_obj, "figure") else plot_obj
    output_path = (FIGURE_DIR / filename).with_suffix(f".{fmt}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, bbox_inches="tight")
    return output_path


def save_model_legend(filename: str = "model_legend", fmt: str = "svg") -> Path:
    """Create and save a standalone model-colour legend as its own figure."""
    handles = [
        Patch(facecolor=color, alpha=0.8, label=MODEL_DISPLAY.get(key, key))
        for key, color in MODEL_PALETTE.items()
    ]
    fig, ax = plt.subplots(figsize=(2.5, 0.4 * len(handles)))
    ax.axis("off")
    ax.legend(
        handles=handles,
        loc="center",
        frameon=False,
        fontsize=10,
        title="Model",
        title_fontsize=11,
    )
    output_path = (FIGURE_DIR / filename).with_suffix(f".{fmt}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)
    return output_path


def make_scatter_comparison(
    df: pd.DataFrame,
    model_name: str,
    x_col: str,
    y_col: str,
    x_label: str,
    y_label: str,
    title: str | None = None,
    delta_col: str | None = None,
    draw_diagonal: bool = False,
    highlight_top_n: int = 3,
    highlight_decreased: bool = True,
    label_all: bool = False,
    draw_regression: bool = False,
    show_correlation: bool = False,
    equal_aspect: bool = True,
) -> plt.Figure:
    """Generic scatter plot comparing two metrics for a single model.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain *x_col*, *y_col*, and ``protein``.  If *delta_col* is
        given, it must also be present for highlighting.
    model_name : str
        Key into ``MODEL_PALETTE`` / ``MODEL_DISPLAY``.
    x_col, y_col : str
        Column names for the x- and y-axes.
    x_label, y_label : str
        Human-readable axis labels.
    title : str or None
        Plot title.  Defaults to ``"{display_name}: {y_label} vs {x_label}"``.
    delta_col : str or None
        Column used to identify improved / decreased points.  If *None*,
        no delta-based highlighting is applied.
    draw_diagonal : bool
        Draw a dashed y = x reference line (useful when both axes share a
        scale, e.g. baseline vs guided RSCC).
    highlight_top_n : int
        Number of most-improved points to highlight (requires *delta_col*).
    highlight_decreased : bool
        Whether to highlight points where *delta_col* < 0.
    label_all : bool
        Annotate every point with its uppercased protein name.
    draw_regression : bool
        Overlay an OLS best-fit line.
    show_correlation : bool
        Display Pearson *r* in a text box (requires >= 5 data points).
    equal_aspect : bool
        If *True*, enforce equal axis scaling.

    Returns
    -------
    matplotlib.figure.Figure
    """
    color = MODEL_PALETTE.get(model_name, "#888888")
    display_name = MODEL_DISPLAY.get(model_name, model_name)
    if title is None:
        title = f"{display_name}: {y_label} vs {x_label}"

    fig, ax = plt.subplots(figsize=(6, 6))

    x, y = df[x_col], df[y_col]

    # Reference diagonal
    if draw_diagonal:
        lo = min(x.min(), y.min()) - 0.02
        hi = max(x.max(), y.max()) + 0.02
        ax.plot([lo, hi], [lo, hi], ls="--", color="#AAAAAA", lw=1, zorder=1)

    # All points
    ax.scatter(x, y, c=color, s=40, alpha=0.7, edgecolors="white", linewidths=0.5, zorder=2)

    # Delta-based highlighting
    if delta_col is not None and delta_col in df.columns:
        if highlight_decreased:
            dec = df[df[delta_col] < 0]
            if not dec.empty:
                ax.scatter(
                    dec[x_col], dec[y_col],
                    facecolors=color, s=60, edgecolors="#D32F2F", linewidths=1.5, zorder=3,
                )
                for _, r in dec.iterrows():
                    ax.annotate(
                        r["protein"].upper(), (r[x_col], r[y_col]),
                        textcoords="offset points", xytext=(6, -8),
                        fontsize=7, color="#D32F2F", fontweight="bold",
                    )

        if highlight_top_n > 0:
            imp = df.nlargest(highlight_top_n, delta_col)
            imp = imp[imp[delta_col] > 0]
            if not imp.empty:
                ax.scatter(
                    imp[x_col], imp[y_col],
                    facecolors=color, s=60, edgecolors="#2E7D32", linewidths=1.5, zorder=3,
                )
                for _, r in imp.iterrows():
                    ax.annotate(
                        r["protein"].upper(), (r[x_col], r[y_col]),
                        textcoords="offset points", xytext=(6, 6),
                        fontsize=7, color="#2E7D32", fontweight="bold",
                    )

    # Label every point (when no delta highlighting)
    if label_all and delta_col is None:
        for _, r in df.iterrows():
            ax.annotate(
                r["protein"].upper(), (r[x_col], r[y_col]),
                textcoords="offset points", xytext=(6, 4),
                fontsize=6, color="#333333",
            )

    # Drop NaNs for regression / correlation (polyfit SVD fails on NaN)
    _finite = np.isfinite(x) & np.isfinite(y)
    _n_dropped = int((~_finite).sum())
    if _n_dropped > 0:
        warnings.warn(
            f"make_scatter_comparison({display_name}): dropping {_n_dropped} of "
            f"{len(df)} rows with NaN in '{x_col}' or '{y_col}' before regression.",
            stacklevel=2,
        )
    _xf, _yf = x[_finite], y[_finite]

    # Regression line
    if draw_regression and len(_xf) >= 2 and _xf.min() != _xf.max():
        m, b = np.polyfit(_xf, _yf, 1)
        x_range = np.linspace(_xf.min(), _xf.max(), 100)
        ax.plot(x_range, m * x_range + b, ls="-", color=color, lw=1.5, alpha=0.6, zorder=1)

    # Pearson correlation
    if show_correlation and len(_xf) >= 5:
        r_val = np.corrcoef(_xf, _yf)[0, 1]
        ax.text(
            0.97, 0.03, f"r = {r_val:.3f}",
            transform=ax.transAxes, fontsize=9, ha="right", va="bottom",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#CCCCCC", alpha=0.8),
        )

    ax.set_xlabel(x_label, fontsize=11)
    ax.set_ylabel(y_label, fontsize=11)
    ax.set_title(title, fontsize=12)
    if equal_aspect:
        ax.set_aspect("equal", adjustable="datalim")
        pad = 0.01
        ax.set_xlim(lo + pad, hi)
        ax.set_ylim(lo + pad, hi)
    fig.tight_layout()
    return fig


def make_rscc_scatter(
    df: pd.DataFrame,
    model_name: str,
    optimal_gw: float,
) -> plt.Figure:
    """Create a scatter plot of baseline vs optimally-guided RSCC for one model."""
    return make_scatter_comparison(
        df,
        model_name,
        x_col="rscc_baseline",
        y_col="rscc_optimal",
        x_label="Baseline RSCC",
        y_label=f"RSCC at gw={optimal_gw}",
        title=f"{MODEL_DISPLAY.get(model_name, model_name)}: RSCC Improvement Over Baseline",
        delta_col="delta_rscc",
        draw_diagonal=True,
        highlight_top_n=3,
        highlight_decreased=True,
    )


def make_aggregated_scatter(
    agg_df: pd.DataFrame,
    *,
    x_col: str,
    y_col: str,
    x_q1_col: str,
    x_q3_col: str,
    y_q1_col: str,
    y_q3_col: str,
    x_label: str,
    y_label: str,
    title: str,
    label_col: str = "protein",
    draw_diagonal: bool = False,
    equal_aspect: bool = False,
    point_color: str = "#4D4D4D",
    label_points: bool = False,
    show_correlation: bool = False,
) -> plt.Figure:
    """Scatter plot with asymmetric IQR error bars for cross-model aggregations.

    Each row of *agg_df* becomes one point; ``x_col``/``y_col`` are the medians
    and ``*_q1``/``*_q3`` provide the IQR whiskers. Error-bar styling mirrors
    :func:`make_aggregated_spiral_figure` so the two figure families stay
    visually consistent.
    """
    fig, ax = plt.subplots(figsize=(6, 6))

    x_med = agg_df[x_col].to_numpy()
    y_med = agg_df[y_col].to_numpy()
    x_q1 = agg_df[x_q1_col].to_numpy()
    x_q3 = agg_df[x_q3_col].to_numpy()
    y_q1 = agg_df[y_q1_col].to_numpy()
    y_q3 = agg_df[y_q3_col].to_numpy()

    xerr = np.vstack([x_med - x_q1, x_q3 - x_med])
    yerr = np.vstack([y_med - y_q1, y_q3 - y_med])

    if draw_diagonal:
        lo = float(min(np.nanmin(x_q1), np.nanmin(y_q1))) - 0.02
        hi = float(max(np.nanmax(x_q3), np.nanmax(y_q3))) + 0.02
        ax.plot([lo, hi], [lo, hi], ls="--", color="#AAAAAA", lw=1, zorder=1)

    ax.errorbar(
        x_med, y_med,
        xerr=xerr, yerr=yerr,
        fmt="none",
        ecolor=point_color,
        alpha=0.35,
        capsize=3,
        lw=1.0,
        zorder=2,
    )
    ax.scatter(
        x_med, y_med,
        c=point_color, s=55,
        edgecolors="white", linewidths=0.5, zorder=4,
    )

    if label_points and label_col in agg_df.columns:
        for _, r in agg_df.iterrows():
            ax.annotate(
                str(r[label_col]).upper(),
                (r[x_col], r[y_col]),
                textcoords="offset points", xytext=(6, 4),
                fontsize=6, color="#333333",
            )

    if show_correlation:
        _finite = np.isfinite(x_med) & np.isfinite(y_med)
        if int(_finite.sum()) >= 5:
            r_val = float(np.corrcoef(x_med[_finite], y_med[_finite])[0, 1])
            ax.text(
                0.97, 0.03, f"r = {r_val:.3f}",
                transform=ax.transAxes, fontsize=9, ha="right", va="bottom",
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#CCCCCC", alpha=0.8),
            )

    ax.set_xlabel(x_label, fontsize=11)
    ax.set_ylabel(y_label, fontsize=11)
    ax.set_title(title, fontsize=12)

    if equal_aspect and draw_diagonal:
        ax.set_aspect("equal", adjustable="datalim")
        pad = 0.01
        ax.set_xlim(lo + pad, hi)
        ax.set_ylim(lo + pad, hi)
    elif equal_aspect:
        ax.set_aspect("equal", adjustable="datalim")

    fig.tight_layout()
    return fig


def make_ridge_plot(
    df: pd.DataFrame,
    plot_column_name: str,
    custom_xlabel: str | None = None,
    smoothing: float = 0.15,
    is_delta: bool = False,
    xlim: tuple[float, float] | None = None,
    combined: bool = False,
    bin_width: float = 0.02,
) -> sns.FacetGrid | plt.Figure:
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
            fontsize=10,
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
        fontsize=9,
        title_fontsize=9,
        bbox_to_anchor=(0.99, 0.99),
    )

    g.set(xlabel=xlabel)
    g.figure.subplots_adjust(hspace=0.6, right=0.95)
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
) -> plt.Figure:
    """Single-axes histplot+KDE with all models overlaid, colored by MODEL_PALETTE."""
    model_order = sorted(df["model_name"].unique())

    values = df[plot_column_name].replace([np.inf, -np.inf], np.nan).dropna()
    lo = float(np.floor(values.min() / bin_width) * bin_width)
    hi = float(np.ceil(values.max() / bin_width) * bin_width) + bin_width / 2
    bins = np.arange(lo, hi, bin_width)

    fig, ax = plt.subplots(figsize=(10, 3.2))

    for model in model_order:
        data = df[df["model_name"] == model]
        color = MODEL_PALETTE.get(model, "#888888")
        sns.histplot(
            data=data,
            x=plot_column_name,
            ax=ax,
            bins=bins,
            kde=True,
            stat="density",
            alpha=0.35,
            linewidth=1.5,
            color=color,
            edgecolor=color,
            kde_kws={"bw_adjust": smoothing},
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
        fontsize=9,
        title_fontsize=9,
    )

    fig.tight_layout()
    return fig


def make_delta_frame(
    df: pd.DataFrame,
    plot_column_name: str,
    baseline_name: str = "baseline",
) -> pd.DataFrame:
    """Compute baseline-relative deltas for a metric column.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with ``protein``, ``model_name``, ``guidance_strength``, and the metric.
    plot_column_name : str
        Column to compute deltas for.
    baseline_name : str
        Value of ``guidance_strength`` that represents the baseline.

    Returns
    -------
    pd.DataFrame
        Merged DataFrame with a ``delta_<plot_column_name>`` column.
    """
    merge_cols = ["protein", "model_name"]
    if "selection_name" in df.columns:
        merge_cols.append("selection_name")
    baseline = df.query("guidance_strength == @baseline_name")[
        merge_cols + [plot_column_name]
    ]
    merged = baseline.merge(
        df.query("guidance_strength != @baseline_name"),
        on=merge_cols,
        suffixes=["_baseline", ""],
    )
    merged["delta_" + plot_column_name] = (
        merged[plot_column_name] - merged[plot_column_name + "_baseline"]
    )
    return merged


def add_tradeoff_guidance_regions(
    ax,
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    guidance_order: list[str],
    guidance_palette: dict[str, str],
) -> None:
    """Overlay smooth guidance-specific KDE regions on a tradeoff axis.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis to draw on.
    df : pd.DataFrame
        Data containing ``guidance_strength``, ``x_col``, and ``y_col``.
    x_col : str
        Column name for the x-axis variable.
    y_col : str
        Column name for the y-axis variable.
    guidance_order : list[str]
        Guidance strengths to draw in order.
    guidance_palette : dict[str, str]
        Mapping from guidance strength to color.
    """
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
                data=subset,
                x=x_col,
                y=y_col,
                ax=ax,
                fill=True,
                levels=4,
                thresh=0.18,
                bw_adjust=1.05,
                cut=1.6,
                color=guidance_palette[gs],
                alpha=0.10,
                warn_singular=False,
                zorder=0,
            )
            sns.kdeplot(
                data=subset,
                x=x_col,
                y=y_col,
                ax=ax,
                fill=False,
                levels=2,
                thresh=0.24,
                bw_adjust=1.05,
                cut=1.6,
                color=guidance_palette[gs],
                alpha=0.35,
                linewidths=0.8,
                warn_singular=False,
                zorder=0.2,
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
    """Create a 2×N model subplot grid showing guidance tradeoff scatter plots.

    Parameters
    ----------
    tradeoff_df : pd.DataFrame
        Tradeoff table with model_name, guidance_strength, metric columns.
    y_col : str
        Column for the y-axis (e.g. ``"geometry_penalty"``).
    y_label : str
        Display label for the y-axis.
    metric_rows : list[tuple[str, str]] | None
        List of ``(column_name, xlabel)`` tuples defining subplot rows.
    model_order : list[str] | None
        Model names defining subplot columns.
    guidance_order : list[str] | None
        Guidance strengths to plot.
    figsize : tuple[float, float] | None
        Figure size.
    scatter_size : int
        Marker size for individual points.
    scatter_alpha : float
        Alpha for individual points.
    show_centroids : bool
        Whether to draw median centroid markers.
    trajectory_protein : str | None
        If set, draw a trajectory path for this protein across guidance levels.
    title : str | None
        Figure suptitle.

    Returns
    -------
    tuple[matplotlib.figure.Figure, np.ndarray]
        ``(fig, axes)`` for optional downstream annotation.
    """
    if metric_rows is None:
        metric_rows = [
            ("delta_rscc", "Δ RSCC vs baseline"),
            ("delta_silhouette_to_ref", "Δ LDDT silhouette-to-ref vs baseline"),
        ]
    if model_order is None:
        model_order = [
            "boltz2X-RAY", "boltz2MD", "protenix", "rf3",
            "protenix-mini", "protenix-tiny", "rf3-old", "rf3-chiral-off",
        ]
    # Keep only models present in the data
    _present = set(tradeoff_df["model_name"].unique())
    model_order = [m for m in model_order if m in _present]
    if guidance_order is None:
        guidance_order = ["weak", "medium", "strong"]
    if figsize is None:
        figsize = (4.25 * len(model_order), 7.6)

    # y-axis limits from data
    y_vals = tradeoff_df[y_col].replace([np.inf, -np.inf], np.nan).dropna()
    y_pad = max(0.16 * (y_vals.max() - y_vals.min()), 0.30)
    y_lim = (float(y_vals.min() - y_pad), float(y_vals.max() + y_pad))

    # x-axis limits per metric row
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
                ax,
                model_df,
                metric,
                y_col,
                guidance_order,
                GUIDANCE_PALETTE_NO_BASELINE,
            )

            for gs in guidance_order:
                gdf = model_df[model_df["guidance_strength"] == gs]
                if gdf.empty:
                    continue
                ax.scatter(
                    gdf[metric],
                    gdf[y_col],
                    s=scatter_size,
                    c=GUIDANCE_PALETTE_NO_BASELINE[gs],
                    alpha=scatter_alpha,
                    marker=GUIDANCE_MARKERS[gs],
                    edgecolors="white",
                    linewidths=0.2,
                    zorder=1.2,
                )
                if show_centroids:
                    ax.scatter(
                        [gdf[metric].median()],
                        [gdf[y_col].median()],
                        s=70,
                        c=GUIDANCE_PALETTE_NO_BASELINE[gs],
                        marker=GUIDANCE_MARKERS[gs],
                        edgecolors="black",
                        linewidths=0.8,
                        zorder=4,
                    )

            if trajectory_protein is not None:
                focus = model_df[model_df["protein"] == trajectory_protein].sort_values(
                    "guidance_order"
                )
                if not focus.empty:
                    ax.plot(
                        focus[metric],
                        focus[y_col],
                        color="black",
                        linewidth=1.0,
                        alpha=0.75,
                        zorder=5,
                    )
                    for row in focus.itertuples(index=False):
                        ax.scatter(
                            [getattr(row, metric)],
                            [getattr(row, y_col)],
                            s=70,
                            facecolors="white",
                            edgecolors="black",
                            linewidths=1.0,
                            marker=GUIDANCE_MARKERS.get(row.guidance_strength, "o"),
                            zorder=6,
                        )

            ax.axvline(0, color="#555555", linestyle="--", linewidth=1.0, alpha=0.6)
            ax.axhline(0, color="#555555", linestyle="--", linewidth=1.0, alpha=0.6)
            if row_xlim[metric] is not None:
                ax.set_xlim(row_xlim[metric])
            ax.set_ylim(y_lim)
            if row_idx == 0:
                ax.set_title(MODEL_DISPLAY.get(model, model), fontsize=11)
            if col_idx == 0:
                ax.set_ylabel(y_label)
            ax.set_xlabel(xlabel)

    # Legend
    label_suffix = " (median)" if show_centroids else ""
    legend_handles = [
        Line2D(
            [0],
            [0],
            marker=GUIDANCE_MARKERS[gs],
            color="none",
            markerfacecolor=GUIDANCE_PALETTE_NO_BASELINE[gs],
            markeredgecolor="black",
            linestyle="none",
            markersize=8 if show_centroids else 7,
            label=f"{gs.capitalize()} guidance{label_suffix}",
        )
        for gs in guidance_order
    ]
    if trajectory_protein is not None:
        legend_handles.append(
            Line2D(
                [0],
                [0],
                color="black",
                linewidth=1.0,
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
        fig.suptitle(title, y=1.05, fontsize=13)
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
    """Add highlight ring markers, guidance dots, and smart labels to a tradeoff grid.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure returned by ``make_tradeoff_grid``.
    axes : np.ndarray
        Axes array returned by ``make_tradeoff_grid``.
    highlight_points : pd.DataFrame
        DataFrame with boolean flag columns and ``point_label``.
    y_col : str
        Column used for y-axis (must match ``make_tradeoff_grid`` call).
    metric_rows : list[tuple[str, str]] | None
        Same as passed to ``make_tradeoff_grid``.
    model_order : list[str] | None
        Same as passed to ``make_tradeoff_grid``.
    guidance_order : list[str] | None
        Same as passed to ``make_tradeoff_grid``.
    title : str | None
        New suptitle (replaces any existing).

    Returns
    -------
    matplotlib.figure.Figure
    """
    if metric_rows is None:
        metric_rows = [
            ("delta_rscc", "Δ RSCC vs baseline"),
            ("delta_silhouette_to_ref", "Δ LDDT silhouette-to-ref vs baseline"),
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
                    subset[metric],
                    subset[y_col],
                    s=style["size"],
                    facecolors="none",
                    edgecolors=style["color"],
                    linewidths=2.0,
                    marker=style["marker"],
                    zorder=4,
                )

            for gs in guidance_order:
                hgdf = hdf[hdf["guidance_strength"] == gs]
                if hgdf.empty:
                    continue
                ax.scatter(
                    hgdf[metric],
                    hgdf[y_col],
                    s=42,
                    c="black",
                    marker=GUIDANCE_MARKERS[gs],
                    zorder=5,
                )

            if row_idx == 0:
                x_min, x_max = ax.get_xlim()
                y_min, y_max = ax.get_ylim()
                if x_max > x_min and y_max > y_min:
                    hdf["x_frac"] = (hdf[metric] - x_min) / (x_max - x_min)
                    hdf["y_frac"] = (hdf[y_col] - y_min) / (y_max - y_min)
                    hdf["label_side"] = hdf["x_frac"].map(
                        lambda v: (
                            "right"
                            if v < 0.22
                            else "left"
                            if v > 0.78
                            else "right"
                            if v < 0.5
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
                                xytext=offset,
                                textcoords="offset points",
                                fontsize=6.4,
                                ha="left" if side == "right" else "right",
                                va="center",
                                bbox={
                                    "boxstyle": "round,pad=0.14",
                                    "fc": "white",
                                    "ec": "#666666",
                                    "alpha": 0.92,
                                    "lw": 0.65,
                                },
                                arrowprops={
                                    "arrowstyle": "-",
                                    "color": "#666666",
                                    "lw": 0.55,
                                    "shrinkA": 0,
                                    "shrinkB": 4,
                                },
                                zorder=6,
                            )

    # Replace legend
    for leg in fig.legends:
        leg.remove()
    legend_handles = [
        Line2D(
            [0],
            [0],
            marker=GUIDANCE_MARKERS[gs],
            color="none",
            markerfacecolor=GUIDANCE_PALETTE_NO_BASELINE[gs],
            markeredgecolor="black",
            linestyle="none",
            markersize=7,
            label=f"{gs.capitalize()} guidance",
        )
        for gs in guidance_order
    ]
    legend_handles.extend(
        [
            Line2D(
                [0],
                [0],
                marker=HIGHLIGHT_STYLES[f]["marker"],
                color="none",
                markerfacecolor="none",
                markeredgecolor=HIGHLIGHT_STYLES[f]["color"],
                markeredgewidth=2.0,
                markersize=8,
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
        fig.suptitle(title, y=1.08, fontsize=13)
    fig.tight_layout()
    return fig


def compute_conceptual_tradeoff_limits(
    tradeoff_df: pd.DataFrame,
    guidance_strengths: tuple[str, ...] = ("weak", "strong"),
) -> tuple[tuple[float, float], tuple[float, float]]:
    """Compute shared RSCC/geometry axis limits for conceptual figures.

    Parameters
    ----------
    tradeoff_df : pd.DataFrame
        Tradeoff table.
    guidance_strengths : tuple[str, ...]
        Guidance buckets whose combined range defines the limits.

    Returns
    -------
    tuple[tuple[float, float], tuple[float, float]]
        ``(xlim, ylim)`` pair.
    """
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
    """Create a simplified RSCC-versus-geometry tradeoff figure.

    Parameters
    ----------
    tradeoff_df : pd.DataFrame
        Tradeoff table.
    guidance_strength : str
        Guidance bucket to visualize.
    xlim : tuple[float, float]
        X-axis limits (used when ``zoom_quantiles`` is None).
    ylim : tuple[float, float]
        Y-axis limits (used when ``zoom_quantiles`` is None).
    title : str
        Figure title.
    zoom_quantiles : tuple[float, float] | None
        Optional quantile range for a zoomed view.

    Returns
    -------
    matplotlib.figure.Figure
    """

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
                    data=sub,
                    x=x_col,
                    y=y_col,
                    ax=ax,
                    fill=True,
                    levels=4,
                    thresh=0.16,
                    bw_adjust=1.0,
                    cut=1.6,
                    color=palette[m],
                    alpha=0.09,
                    warn_singular=False,
                    zorder=0,
                )
                sns.kdeplot(
                    data=sub,
                    x=x_col,
                    y=y_col,
                    ax=ax,
                    fill=False,
                    levels=2,
                    thresh=0.22,
                    bw_adjust=1.0,
                    cut=1.6,
                    color=palette[m],
                    alpha=0.30,
                    linewidths=0.8,
                    warn_singular=False,
                    zorder=0.2,
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
            mdf["delta_rscc"],
            mdf["geometry_penalty"],
            s=20,
            c=c,
            alpha=0.24,
            edgecolors="white",
            linewidths=0.25,
            zorder=1.2,
        )
        ax.scatter(
            [mdf["delta_rscc"].median()],
            [mdf["geometry_penalty"].median()],
            s=95,
            c=c,
            edgecolors="black",
            linewidths=0.9,
            zorder=4,
        )

    ax.axvline(0, color="#555555", linestyle="--", linewidth=1.0, alpha=0.6)
    ax.axhline(0, color="#555555", linestyle="--", linewidth=1.0, alpha=0.6)
    ax.set_xlim(plot_xlim)
    ax.set_ylim(plot_ylim)
    ax.set_xlabel("Δ RSCC vs baseline")
    ax.set_ylabel("Geometry penalty")
    ax.set_title(title, fontsize=13)
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
    """Create a single-model guidance spiral trajectory figure.

    Parameters
    ----------
    model_df : pd.DataFrame
        Rows from tradeoff_df for one model.
    model_key : str
        Internal model name (for display and filename).
    label_top_n : int
        Number of most-extreme trajectories to label.
    xlim, ylim : tuple[float, float] or None
        Fixed axis limits. When *None* (default) the limits are computed
        from *model_df* alone (original behaviour).

    Returns
    -------
    matplotlib.figure.Figure
    """
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
            xs[1:],
            ys[1:],
            c=dot_colors,
            s=26,
            edgecolors="white",
            linewidths=0.3,
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
            fontsize=7,
            fontweight="bold",
            color="#333333",
            bbox=dict(
                boxstyle="round,pad=0.15",
                fc="white",
                ec="#aaaaaa",
                alpha=0.85,
                lw=0.5,
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
    ax.set_xlabel("Δ RSCC vs. baseline", fontsize=11)
    ax.set_ylabel("∆ Geometry penalty vs. baseline", fontsize=11)

    display = MODEL_DISPLAY.get(model_key, model_key)
    ax.set_title(
        f"{display}: guidance strength tradeoff", fontsize=13, fontweight="bold"
    )

    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.75, pad=0.02)
    cbar.set_label("Guidance weight", fontsize=10)

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
):
    """Aggregated guidance-weight spiral trajectory figure.

    For each (model, guidance_weight) group, collapses across proteins to a
    single (median Δ RSCC, median geometry penalty) point with IQR error bars
    in both X and Y. The Y axis uses symmetric log scaling to keep small
    guidance weights legible.

    Parameters
    ----------
    tradeoff_df : pd.DataFrame
        The full tradeoff dataframe (all models). Must include columns
        ``protein``, ``model_name``, ``guidance_weight``, ``delta_rscc``,
        and ``geometry_penalty``.
    model_key : str or None
        If provided, render a single-model aggregated figure. If ``None``,
        render a multi-model overlay.
    model_order : list[str] or None
        Models to include in multi-model mode (default: all present).
    xlim, ylim : tuple[float, float] or None
        Fixed axis limits. When ``None`` they are derived from the
        aggregated whiskers with padding.
    linthresh : float
        Linear-region threshold for the symlog Y axis.
    min_n : int
        Minimum number of proteins required in a group; smaller groups are
        dropped (IQR is unreliable with very few points).

    Returns
    -------
    matplotlib.figure.Figure
    """
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

    # For single-model mode, use a flare colormap on guidance weight to
    # match the visual language of make_spiral_figure.
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

        # Asymmetric error arrays (lower, upper) for errorbar
        xerr = np.vstack([x_med - x_q1, x_q3 - x_med])
        yerr = np.vstack([y_med - y_q1, y_q3 - y_med])

        # Trajectory anchored at origin (baseline)
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
            ecolor=line_color,
            alpha=0.35,
            capsize=3,
            lw=1.0,
            zorder=2,
        )

        if single_mode:
            dot_colors = flare_cmap(gw_norm(gws))
            ax.scatter(
                x_med, y_med,
                c=dot_colors, s=70,
                edgecolors="white", linewidths=0.6, zorder=4,
            )
        else:
            ax.scatter(
                x_med, y_med,
                c=[line_color], s=55,
                edgecolors="white", linewidths=0.5, zorder=4,
            )

        # Arrowhead on the final segment
        if len(xs) >= 2:
            ax.annotate(
                "",
                xy=(xs[-1], ys[-1]),
                xytext=(xs[-2], ys[-2]),
                arrowprops=dict(
                    arrowstyle="-|>",
                    color=line_color,
                    lw=0.0,
                    mutation_scale=10,
                ),
                zorder=5,
            )

    ax.axvline(0, color="#555555", ls="--", lw=0.8, alpha=0.5)
    ax.axhline(0, color="#555555", ls="--", lw=0.8, alpha=0.5)

    ax.set_yscale("symlog", linthresh=linthresh)

    if xlim is not None and ylim is not None:
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
    else:
        all_x = np.concatenate([
            agg["x_q1"].to_numpy(),
            agg["x_q3"].to_numpy(),
            np.array([0.0]),
        ])
        all_y = np.concatenate([
            agg["y_q1"].to_numpy(),
            agg["y_q3"].to_numpy(),
            np.array([0.0]),
        ])
        x_lo, x_hi = float(all_x.min()), float(all_x.max())
        y_lo, y_hi = float(all_y.min()), float(all_y.max())
        x_pad = max(0.08 * (x_hi - x_lo), 0.02)
        y_pad = max(0.08 * (y_hi - y_lo), 0.50)
        ax.set_xlim(x_lo - x_pad, x_hi + x_pad)
        ax.set_ylim(y_lo - y_pad, y_hi + y_pad)

    ax.set_xlabel("Δ RSCC vs. baseline (median)", fontsize=11)
    ax.set_ylabel("∆ Geometry penalty vs. baseline (median)", fontsize=11)

    if model_key is not None:
        display = MODEL_DISPLAY.get(model_key, model_key)
        ax.set_title(
            f"{display}: aggregated guidance tradeoff",
            fontsize=13, fontweight="bold",
        )
        sm = ScalarMappable(cmap=flare_cmap, norm=gw_norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, shrink=0.75, pad=0.02)
        cbar.set_label("Guidance weight", fontsize=10)
    else:
        ax.set_title(
            "Aggregated guidance strength tradeoff",
            fontsize=13, fontweight="bold",
        )
        ax.legend(
            loc="best",
            frameon=True,
            framealpha=0.9,
            edgecolor="none",
            fontsize=9,
            title="Model",
            title_fontsize=10,
        )

    fig.tight_layout()
    return fig


# ── Utility ───────────────────────────────────────────────────────────────────


def robust_z_score(series: pd.Series) -> pd.Series:
    """MAD-based robust z-score.

    Parameters
    ----------
    series : pd.Series
        Input values.

    Returns
    -------
    pd.Series
        Robust z-scores.
    """
    median = series.median()
    mad = (series - median).abs().median()
    scale = 1.4826 * mad + 1e-9
    return (series - median) / scale
