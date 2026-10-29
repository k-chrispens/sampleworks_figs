"""Publication styling, palettes, constants, and figure-save helpers.

Pure presentation-layer module: rcParams, colors, model/guidance/classification
palettes, model groups, filter constants, and the three figure-save helpers
(``save_fig``, ``save_model_legend``, ``save_per_model_figures``).
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import seaborn as sns

# ── Publication styling ───────────────────────────────────────────────────────
sns.set_context("paper")
sns.set_style("ticks")
plt.rcParams.update({
    "axes.titlesize": 17,
    "axes.labelsize": 15,
    "xtick.labelsize": 13,
    "ytick.labelsize": 13,
    "legend.fontsize": 12,
    "legend.title_fontsize": 13,
    "figure.titlesize": 18,
    "font.size": 12,
})

# ── Constants ──────────────────────────────────────────────────────────────────

OCC_FILTER = "0.5occA_0.5occB"
OCC_FILTER_ALTLOC = "{'A': 0.5, 'B': 0.5}"

# Pre-cutoff training-distribution analysis. The training-cluster CSV
# (``lddt_assignment_counts.csv``) holds, per (protein, selection), the
# number of pre-cutoff PDB cluster members assigned to altloc A vs B by
# lDDT comparison. Four proteins from the initial-40 set are absent from
# the cutoffs CSV; all four are pre-2018 depositions, so they are filled
# in as "Pre".
TRAINING_CUTOFF_DATE = "2021-09-30"

# PyMOL "good" colors (from ~/.pymolrc, cred. Sam Pellock / RFdiffusion)
COLOR_GOOD_GRAY = "#DCDCDC"
COLOR_GOOD_DARKBLUE = "#4B5FAA"
COLOR_GOOD_TEAL = "#4FB9AF"
COLOR_GOOD_NAVAHO = "#FFE0AC"
COLOR_GOOD_PINK = "#FFACB7"
COLOR_GOOD_LIGHTBLUE = "#8EDCE6"

# Model groups for split plotting
MODELS_ALT = {"protenix-mini", "protenix-tiny", "rf3-chiral-off"}
MODELS_PROTENIX_GROUP = ["protenix", "protenix-mini", "protenix-tiny"]
MODELS_RF3_GROUP = ["rf3", "rf3-chiral-off"]

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

CLASSIFICATION_PALETTE = {
    "side_chain_only": "#E7298A",  # magenta (ColorBrewer Dark2)
    "small_loop":      "#1B9E77",  # teal-green
    "large_loop":      "#D95F02",  # orange
    "domain_shift":    "#7570B3",  # purple
}
CLASSIFICATION_DISPLAY = {
    "side_chain_only": "Side-Chain Only",
    "small_loop":      "Small Loop",
    "large_loop":      "Large Loop",
    "domain_shift":    "Domain Shift",
}
CLASSIFICATION_ORDER = ["side_chain_only", "small_loop", "large_loop", "domain_shift"]

# Single swap-point for the by-class grid's classification scheme. To switch to a
# different scheme, repoint this and (if it introduces new labels) extend
# CLASSIFICATION_ORDER / CLASSIFICATION_PALETTE / CLASSIFICATION_DISPLAY above.
CLASSIFICATION_CSV = "results/classify_altloc_regions.csv"

MODEL_NAME_ALIAS = {
    "boltz2_X-RAY_DIFFRACTION": "boltz2X-RAY",
    "boltz2_MD": "boltz2MD",
    "boltz2xray": "boltz2X-RAY",
    "boltz2md": "boltz2MD",
    "rf3-fixed": "rf3",
}

NICE_XLABELS = {
    "delta_clashscore": "Change in Clashscore (Guided − Baseline)",
    "delta_outlier_fraction_bond_angles": "Change in Bond Angle Outlier Fraction (Guided − Baseline)",
    "delta_outlier_fraction_bond_lengths": "Change in Bond Length Outlier Fraction (Guided − Baseline)",
    "delta_rscc": "Change in RSCC (Guided − Baseline)",
    "delta_clustering_score": "Change in Clustering Score (Guided − Baseline)",
    "clashscore": "Clashscore (steric clashes per 1000 atoms)",
    "rscc": "Real-Space Correlation Coefficient (RSCC)",
    "clustering_score": "Clustering Score",
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

OPTIMAL_GW = {
    "boltz2X-RAY": 0.1,
    "boltz2MD": 0.1,
    "protenix": 0.1,
    "protenix-mini": 0.1,
    "protenix-tiny": 0.1,
    "rf3": 0.01,
    "rf3-chiral-off": 0.01,
}

FIGURE_DIR = Path("figures")


def save_fig(
    plot_obj,
    filename: str,
    fmt: str = "svg",
    annotate_source: bool = True,
) -> Path:
    """Save a Matplotlib or Plotly figure to the figures directory.

    Parameters
    ----------
    plot_obj : matplotlib.figure.Figure or seaborn.FacetGrid
        Object to save.
    filename : str
        Output filename (relative to ``FIGURE_DIR``).  The extension is
        replaced by *fmt* so existing ``.png`` names are handled automatically.
    fmt : str
        Output format (default ``"svg"``).
    annotate_source : bool
        If True, add a small source annotation to the figure corner.

    Returns
    -------
    Path
        Path to the saved file.
    """
    figure = plot_obj.figure if hasattr(plot_obj, "figure") else plot_obj
    output_path = (FIGURE_DIR / filename).with_suffix(f".{fmt}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if hasattr(figure, "write_image"):
        figure.write_image(output_path, format=fmt)
    else:
        figure.savefig(output_path, bbox_inches="tight")
    return output_path


def save_model_legend(filename: str = "model_legend", fmt: str = "svg") -> Path:
    """Create and save a standalone model-colour legend as its own figure."""
    import plotly.graph_objects as go

    figure = go.Figure([
        go.Scatter(x=[None], y=[None], mode="markers", name=MODEL_DISPLAY.get(key, key),
                   marker=dict(color=color, size=10))
        for key, color in MODEL_PALETTE.items()
    ])
    figure.update_layout(template="plotly_white", height=350, width=350,
                         xaxis=dict(visible=False), yaxis=dict(visible=False))
    return save_fig(figure, filename, fmt=fmt)


def save_per_model_figures(
    df,
    model_col: str,
    make_fn,
    path_template: str,
    *,
    make_kwargs: dict | None = None,
    close: bool = True,
) -> None:
    """Save one figure per unique model in *df*.

    Parameters
    ----------
    make_fn : callable
        ``make_fn(model_df, model_name, **make_kwargs) -> Figure``.
    path_template : str
        Format string with ``{slug}`` placeholder.
    """
    from transforms import model_slug

    make_kwargs = make_kwargs or {}
    for _model in sorted(df[model_col].unique()):
        _model_df = df[df[model_col] == _model]
        if _model_df.empty:
            continue
        fig = make_fn(_model_df, _model, **make_kwargs)
        save_fig(fig, path_template.format(slug=model_slug(_model)))
        if close and not hasattr(fig, "write_image"):
            plt.close(fig)
