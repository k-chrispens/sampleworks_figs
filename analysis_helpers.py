"""Shared helpers for the initial-40 analysis notebooks.

Data loading, tagging, plotting utilities, and constants used across
01_geometry, 02_fit_quality, and 03_tradeoffs notebooks.
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
from adjustText import adjust_text
from matplotlib.cm import ScalarMappable
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize, to_rgb
from matplotlib.legend_handler import HandlerTuple
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

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
_MISSING_CUTOFF_FILL = {"1vme": "Pre", "4ole": "Pre", "5sop": "Pre", "6b8x": "Pre"}

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


def filter_by_model_group(
    df: pd.DataFrame,
    group: str,
    extra_query: str | None = None,
) -> tuple[pd.DataFrame, str]:
    """Filter a dataframe by model group and return (filtered_df, filename_suffix)."""
    if group == "main":
        result, suffix = df.query("model_name not in @MODELS_ALT"), ""
    elif group == "protenix":
        result, suffix = df.query("model_name in @MODELS_PROTENIX_GROUP"), "_protenix_alts"
    else:
        result, suffix = df.query("model_name in @MODELS_RF3_GROUP"), "_rf3_alts"
    if extra_query:
        result = result.query(extra_query)
    return result, suffix


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


def model_slug(model_name: str) -> str:
    """Return a filesystem-safe slug for a model name."""
    display = MODEL_DISPLAY.get(model_name, model_name)
    return display.lower().replace("-", "_").replace(" ", "_")


def format_residue_range(sel: str) -> str:
    """Extract residue range from a selection string, condensing single-residue ranges.

    ``"chain A and resi 181-181"`` -> ``"181"``;
    ``"chain A and resi 178-193"`` -> ``"178-193"``.
    """
    m = re.search(r"resi\s+(\d+)-(\d+)", str(sel))
    if m:
        s, e = m.group(1), m.group(2)
        return s if s == e else f"{s}-{e}"
    m2 = re.search(r"resi\s+([\d-]+)", str(sel))
    return m2.group(1) if m2 else str(sel)


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
    "delta_silhouette_to_ref": "Change in Clustering Score (Guided − Baseline)",
    "clashscore": "Clashscore (steric clashes per 1000 atoms)",
    "rscc": "Real-Space Correlation Coefficient (RSCC)",
    "avg_silhouette_to_ref": "Clustering Score",
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
    "rf3-old": 0.01,
    "rf3-chiral-off": 0.01,
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
    df = pd.read_csv(csv_path, dtype={"error": str}).fillna({"error": ""})
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


def get_rscc_nan_proteins(
    csv_path: str = "initial40results/rscc_results.csv",
    models: tuple[str, ...] = ("boltz2X-RAY", "boltz2MD", "protenix"),
) -> set[str]:
    """Return proteins with at least one NaN ``rscc`` value for the given models.

    Used by the notebooks' "Exclude NaN-RSCC proteins" toggle to produce a
    consistent cross-notebook blacklist. Reuses ``load_rscc`` so the occupancy
    filter and ``MODEL_NAME_ALIAS`` pass stay in sync.
    """
    df = load_rscc(csv_path)
    mask = df["model_name"].isin(models) & df["rscc"].isna()
    return set(df.loc[mask, "protein"].unique())


def _parse_b_occupancy(s) -> float:
    """Extract the B-altloc occupancy (second value) from an ``occupancies`` cell.

    The ``occupancies`` column in the LDDT results is stored as the numpy
    string repr of a 2-element array, e.g. ``"[0.375 0.625]"`` — first value
    is A, second is B.
    """
    if not isinstance(s, str):
        return np.nan
    parts = s.strip("[] ").split()
    if len(parts) < 2:
        return np.nan
    try:
        return float(parts[1])
    except ValueError:
        return np.nan


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
    df = pd.read_csv(csv_path, dtype={"error": str}).fillna({"error": ""})
    _occ = OCC_FILTER_ALTLOC
    df = df.query("altloc_occupancies == @_occ").copy()
    df["model_name"] = df.apply(
        lambda x: x["model"] + x["method"] if x["model"] == "boltz2" else x["model"],
        axis=1,
    )
    df["model_name"] = df["model_name"].map(lambda x: MODEL_NAME_ALIAS.get(x, x))
    return df[_cols]


def load_classification(csv_path: str) -> pd.DataFrame:
    """Load conformational-change classification with lowercased protein names.

    Parameters
    ----------
    csv_path : str
        Path to a classification CSV (e.g. ``classify_altloc_regions.csv`` or
        ``small-altlocs-classified-rmsd.csv``).

    Returns
    -------
    pd.DataFrame
        Columns: ``protein``, ``selection``, ``classification``, ``chain``,
        ``start_res``, ``end_res``. ``protein`` is lowercased.
    """
    df = pd.read_csv(csv_path)
    df["protein"] = df["protein"].str.lower()
    return df[["protein", "selection", "classification", "chain", "start_res", "end_res"]]


_RSCC_SEL_PATTERN = re.compile(r"chain\s+(\S+)\s+and\s+resi\s+(\d+)-(\d+)")


def merge_classification_by_containment(
    rscc_df: pd.DataFrame,
    classification_df: pd.DataFrame,
) -> pd.DataFrame:
    """Inner-join RSCC rows to classifications by residue-range containment.

    For each RSCC selection ``chain X and resi A-B``, finds the classification
    row with the same protein and chain whose ``[start_res, end_res]`` fully
    contains ``[A, B]``. Exact-match selections (most rows) match through this
    same path. If multiple classifications contain the same range, the
    narrowest one wins.

    Selections that don't parse, or have no containing classification, are
    dropped; a one-line summary is printed.
    """
    parsed = rscc_df["selection"].str.extract(_RSCC_SEL_PATTERN)
    parsed.columns = ["chain", "start_res", "end_res"]
    parseable = parsed["chain"].notna()
    if (~parseable).any():
        print(
            f"⚠️  {(~parseable).sum()} RSCC rows had unparseable selections "
            f"(e.g. {rscc_df.loc[~parseable, 'selection'].iloc[0]!r}); dropped."
        )
    left = rscc_df[parseable].copy()
    left["start_res"] = parsed.loc[parseable, "start_res"].astype(int)
    left["end_res"] = parsed.loc[parseable, "end_res"].astype(int)
    left["chain"] = parsed.loc[parseable, "chain"]

    cls = classification_df.rename(columns={
        "start_res": "_cls_start",
        "end_res": "_cls_end",
        "selection": "_cls_selection",
    })

    merged = left.merge(cls, on=["protein", "chain"], how="left")
    contained = (
        (merged["_cls_start"] <= merged["start_res"])
        & (merged["_cls_end"] >= merged["end_res"])
    )
    merged = merged[contained].copy()
    merged["_cls_span"] = merged["_cls_end"] - merged["_cls_start"]
    merged = (
        merged
        .sort_values("_cls_span")
        .drop_duplicates(subset=list(left.columns), keep="first")
    )

    matched_keys = merged.set_index(list(left.columns)).index
    left_keys = left.set_index(list(left.columns)).index
    unmatched = ~left_keys.isin(matched_keys)
    if unmatched.any():
        n = int(unmatched.sum())
        sample = left.loc[unmatched.values, ["protein", "selection"]].drop_duplicates().head(3)
        print(
            f"⚠️  {n} RSCC rows had no containing classification; dropped. "
            f"Examples: {sample.to_dict(orient='records')}"
        )

    return merged.drop(columns=["_cls_start", "_cls_end", "_cls_selection", "_cls_span"])


def load_deposition_cutoffs(
    path: str = "initial40results/Sampleworks_initial_structures.csv",
) -> pd.DataFrame:
    """Load per-protein deposition dates and pre/post-cutoff classification.

    Lowercases the PDB ID, drops duplicate rows, and fills in any proteins
    missing from the CSV that the analysis still uses (see
    ``_MISSING_CUTOFF_FILL``). Cutoff date is :data:`TRAINING_CUTOFF_DATE`.

    Returns
    -------
    pd.DataFrame
        Columns: ``protein`` (lowercase 4-char), ``deposit_date``,
        ``cutoff_status`` (``"Pre"`` or ``"Post"``).
    """
    df = pd.read_csv(path)
    df = df.rename(columns={
        "PDB_ID": "protein",
        "Date_Deposited": "deposit_date",
        "Pre_post_cutoff": "cutoff_status",
    })
    df["protein"] = df["protein"].str.lower()
    df = df.drop_duplicates(subset=["protein"]).reset_index(drop=True)
    missing = [p for p in _MISSING_CUTOFF_FILL if p not in set(df["protein"])]
    if missing:
        fill_rows = pd.DataFrame({
            "protein": missing,
            "deposit_date": [pd.NaT] * len(missing),
            "cutoff_status": [_MISSING_CUTOFF_FILL[p] for p in missing],
        })
        df = pd.concat([df, fill_rows], ignore_index=True)
        print(
            f"ℹ️  Filled {len(missing)} protein(s) missing from {path}: "
            f"{sorted(missing)} → Pre (all are pre-2018 depositions)."
        )
    return df[["protein", "deposit_date", "cutoff_status"]]


def load_training_cluster_counts(
    path: str = "lddt_assignment_counts.csv",
) -> pd.DataFrame:
    """Load per-selection pre-cutoff training-cluster counts.

    Lowercases ``protein`` and renames ``selection_part`` → ``selection`` so
    the frame joins directly onto LDDT/RSCC results. Adds ``frac_a_cluster``
    (fraction of cluster members assigned to altloc A by lDDT comparison)
    and ``x_one_minus_fa = 1 − frac_a_cluster``. Both are NaN where the
    cluster has no members assigned to either altloc (``A == B == 0``).
    """
    df = pd.read_csv(path)
    df["protein"] = df["protein"].str.lower()
    # The CSV's ``selection`` column is the full multi-part selection string;
    # ``selection_part`` is the per-residue-range piece that joins to LDDT.
    # Drop the multi-part selection first so the rename doesn't collide.
    df = df.drop(columns=["selection"]).rename(columns={"selection_part": "selection"})
    n_a = df["n_assigned_to_A"]
    n_b = df["n_assigned_to_B"]
    total = n_a + n_b
    df["frac_a_cluster"] = np.where(total > 0, n_a / total.replace(0, np.nan), np.nan)
    df["x_one_minus_fa"] = 1.0 - df["frac_a_cluster"]
    return df[[
        "protein", "selection",
        "n_assigned_to_A", "n_assigned_to_B",
        "n_cluster_evaluated",
        "frac_a_cluster", "x_one_minus_fa",
    ]]


# ── Plotting helpers ──────────────────────────────────────────────────────────


def save_fig(
    plot_obj,
    filename: str,
    fmt: str = "svg",
    annotate_source: bool = True,
) -> Path:
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
        title="Model",
    )
    output_path = (FIGURE_DIR / filename).with_suffix(f".{fmt}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)
    return output_path


def save_per_model_figures(
    df: pd.DataFrame,
    model_col: str,
    make_fn: callable,
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
    make_kwargs = make_kwargs or {}
    for _model in sorted(df[model_col].unique()):
        _model_df = df[df[model_col] == _model]
        if _model_df.empty:
            continue
        fig = make_fn(_model_df, _model, **make_kwargs)
        save_fig(fig, path_template.format(slug=model_slug(_model)))
        if close:
            plt.close(fig)


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
    highlight_decreased: bool | int = True,
    highlight_proteins: list[str] | None = None,
    label_col: str = "protein",
    label_all: bool = False,
    draw_regression: bool = False,
    show_correlation: bool = False,
    equal_aspect: bool = True,
    show_marginal_kde: bool = False,
    xlim: tuple[float, float] | None = None,
    ylim: tuple[float, float] | None = None,
    residual: bool = False,
    label_top_n: int = 3,
    label_bottom_n: int = 5,
    label_far_from_ceiling_n: int = 5,
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
        scale, e.g. baseline vs guided RSCC). Ignored when ``residual=True``.
    highlight_top_n : int
        Number of most-improved points to highlight (requires *delta_col*).
    highlight_decreased : bool | int
        ``True`` highlights every point with ``delta_col`` < 0. An integer
        ``n`` highlights only the ``n`` most-decreased points (smallest
        ``delta_col`` values below zero). ``False`` disables the highlight.
    label_col : str
        Column to use for the text labels and ``highlight_proteins`` matching.
        Defaults to ``"protein"`` — pass ``"label"`` (or similar) when
        plotting per-(protein, selection) data.
    label_all : bool
        Annotate every point with its uppercased label.
    draw_regression : bool
        Overlay an OLS best-fit line.
    show_correlation : bool
        Display Pearson *r* in a text box (requires >= 5 data points).
    equal_aspect : bool
        If *True*, enforce equal axis scaling. Forced ``False`` when
        ``residual=True`` because residual y has a different scale than x.
    residual : bool
        If *True*, plot ``y_col - x_col`` on the y-axis instead of ``y_col``,
        replace the y = x diagonal with a y = 0 horizontal reference and a
        y = 1 − x ceiling envelope. In residual mode, four label sets are
        drawn: (a) the ``label_top_n`` most improved (delta > 0, green ring),
        (b) every point with delta < 0 (red ring), (c) the
        ``label_far_from_ceiling_n`` points furthest below the y = 1 − x
        ceiling (orange ring), and (d) any point whose ``label_col`` value
        matches an entry in ``highlight_proteins`` (blue ring overlay,
        stacked on top of the other rings). Marker styling switches to
        small, un-edged, alpha=0.3 dots. Single-row data has no error to
        propagate; the propagation rule sqrt(sd_g^2 + sd_b^2) applies in
        :func:`make_aggregated_scatter`.
    label_top_n : int
        Active only when ``residual=True``: number of most-improved points
        to label.
    label_bottom_n : int
        Retained for API compatibility; ignored in residual mode (every
        decreased point is labeled instead).
    label_far_from_ceiling_n : int
        Active only when ``residual=True``: number of points furthest below
        the y = 1 − x ceiling envelope to label (low guided value).

    Returns
    -------
    matplotlib.figure.Figure
    """
    color = MODEL_PALETTE.get(model_name, "#888888")
    display_name = MODEL_DISPLAY.get(model_name, model_name)
    if title is None:
        title = f"{display_name}: {y_label} vs {x_label}"

    _w = 8 if label_all else 6
    if show_marginal_kde:
        fig = plt.figure(figsize=(_w, _w * 1.25))
        gs = fig.add_gridspec(2, 1, height_ratios=[4, 1], hspace=0.05)
        ax = fig.add_subplot(gs[0])
        ax_kde = fig.add_subplot(gs[1], sharex=ax)
    else:
        fig, ax = plt.subplots(figsize=(_w, _w))
        ax_kde = None

    x_raw = df[x_col]
    y_raw = df[y_col]
    if residual:
        # Residual = guided − baseline (improvement over baseline)
        y_plot = y_raw - x_raw
        # Override y-label with residual phrasing unless caller passed
        # a custom one that already speaks of "Improvement".
        if "Improvement" not in y_label:
            y_label_used = "Improvement (Guided − Baseline RSCC)"
        else:
            y_label_used = y_label
        # Residual delta is the y-value itself; precompute for highlight selection
        _delta_for_highlight = y_plot
        # Equal aspect doesn't make sense across two different scales
        equal_aspect = False
        # Ignore the diagonal flag; residual mode owns the reference geometry
        draw_diagonal_used = False
    else:
        y_plot = y_raw
        y_label_used = y_label
        _delta_for_highlight = df[delta_col] if (delta_col is not None and delta_col in df.columns) else None
        draw_diagonal_used = draw_diagonal

    x = x_raw

    # Reference geometry
    if draw_diagonal_used:
        lo = min(x.min(), y_plot.min()) - 0.02
        hi = max(x.max(), y_plot.max()) + 0.02
        if xlim is not None:
            lo = min(lo, xlim[0])
            hi = max(hi, xlim[1])
        if ylim is not None:
            lo = min(lo, ylim[0])
            hi = max(hi, ylim[1])
        ax.plot([lo, hi], [lo, hi], ls="--", color="#AAAAAA", lw=1, zorder=1)

    if residual:
        # y = 0 horizontal reference (replaces y = x diagonal); same dashed style.
        _x_lo = float(np.nanmin(x)) if xlim is None else xlim[0]
        _x_hi = float(np.nanmax(x)) if xlim is None else xlim[1]
        # Pad a bit so the line spans visible range
        _x_lo = min(_x_lo, 0.0) - 0.02
        _x_hi = max(_x_hi, 1.0) + 0.02
        ax.axhline(0.0, ls="--", color="#AAAAAA", lw=1, zorder=1)
        # Upper envelope y = 1 - x — guided cannot exceed 1.0, so improvement is bounded.
        _xs = np.linspace(_x_lo, _x_hi, 100)
        ax.plot(_xs, 1.0 - _xs, ls="--", color="#AAAAAA", lw=0.8, zorder=1)

    # All points — styling depends on mode
    if residual:
        ax.scatter(x, y_plot, c=color, s=14, alpha=0.3, edgecolors="none", zorder=2)
    else:
        ax.scatter(x, y_plot, c=color, s=40, alpha=0.7, edgecolors="white", linewidths=0.5, zorder=2)

    # Collect every label we want to draw; adjust_text will reposition them
    # all together at the end to avoid mutual overlaps.
    texts: list = []

    if residual:
        # Residual mode label sets:
        #   1. top-N most improved (delta > 0)            — green ring
        #   2. every point with delta < 0                 — red ring
        #   3. N points furthest below the y = 1 − x line — orange ring
        #   4. highlight_proteins                          — blue ring overlay
        # The third criterion ranks by 1 − y_col (large = low guided RSCC,
        # low spot below the ceiling); equivalent to perpendicular distance
        # below the line up to a constant factor.
        _has_label_col = label_col in df.columns
        _seen_idx: set = set()

        if label_top_n > 0 and _delta_for_highlight is not None:
            imp = df.assign(_d=_delta_for_highlight).nlargest(label_top_n, "_d")
            imp = imp[imp["_d"] > 0]
            if not imp.empty:
                ax.scatter(
                    imp[x_col], imp[y_col] - imp[x_col],
                    facecolors=color, s=24, edgecolors="#2E7D32", linewidths=1.2, zorder=3,
                )
                if _has_label_col:
                    for _, r in imp.iterrows():
                        _seen_idx.add(r.name)
                        texts.append(ax.text(
                            r[x_col], r[y_col] - r[x_col], str(r[label_col]).upper(),
                            color="#333333", zorder=10,
                        ))

        if _delta_for_highlight is not None:
            dec = df.assign(_d=_delta_for_highlight)
            dec = dec[dec["_d"] < 0]
            dec = dec.loc[~dec.index.isin(_seen_idx)]
            if not dec.empty:
                ax.scatter(
                    dec[x_col], dec[y_col] - dec[x_col],
                    facecolors=color, s=24, edgecolors="#D32F2F", linewidths=1.2, zorder=3,
                )
                if _has_label_col:
                    for _, r in dec.iterrows():
                        _seen_idx.add(r.name)
                        texts.append(ax.text(
                            r[x_col], r[y_col] - r[x_col], str(r[label_col]).upper(),
                            color="#333333", zorder=10,
                        ))

        if label_far_from_ceiling_n > 0:
            below_dist = 1.0 - df[y_col]
            far = df.assign(_b=below_dist).nlargest(label_far_from_ceiling_n, "_b")
            far = far[far["_b"] > 0]
            far = far.loc[~far.index.isin(_seen_idx)]
            if not far.empty:
                ax.scatter(
                    far[x_col], far[y_col] - far[x_col],
                    facecolors=color, s=24, edgecolors="#F57C00", linewidths=1.2, zorder=3,
                )
                if _has_label_col:
                    for _, r in far.iterrows():
                        _seen_idx.add(r.name)
                        texts.append(ax.text(
                            r[x_col], r[y_col] - r[x_col], str(r[label_col]).upper(),
                            color="#333333", zorder=10,
                        ))

        if highlight_proteins and _has_label_col:
            # Blue ring overlay; stacks on top of any green/red/orange ring.
            wanted = {p.lower() for p in highlight_proteins}
            picked = df[df[label_col].str.lower().isin(wanted)]
            if not picked.empty:
                ax.scatter(
                    picked[x_col], picked[y_col] - picked[x_col],
                    facecolors=color, s=24, edgecolors="#1565C0", linewidths=1.2, zorder=4,
                )
                for _, r in picked.iterrows():
                    if r.name in _seen_idx:
                        continue
                    _seen_idx.add(r.name)
                    texts.append(ax.text(
                        r[x_col], r[y_col] - r[x_col], str(r[label_col]).upper(),
                        color="#333333", zorder=10,
                    ))

    else:
        # Legacy (non-residual) delta-based highlighting
        if delta_col is not None and delta_col in df.columns:
            if highlight_decreased:
                dec = df[df[delta_col] < 0]
                if isinstance(highlight_decreased, bool):
                    pass  # keep all decreased
                else:
                    # Integer count: keep only the N most-decreased rows.
                    dec = dec.nsmallest(int(highlight_decreased), delta_col)
                if not dec.empty:
                    ax.scatter(
                        dec[x_col], dec[y_col],
                        facecolors=color, s=60, edgecolors="#D32F2F", linewidths=1.5, zorder=3,
                    )
                    for _, r in dec.iterrows():
                        texts.append(ax.text(
                            r[x_col], r[y_col], str(r[label_col]).upper(),
                            color="#333333", zorder=10,
                        ))

            if highlight_top_n > 0:
                imp = df.nlargest(highlight_top_n, delta_col)
                imp = imp[imp[delta_col] > 0]
                if not imp.empty:
                    ax.scatter(
                        imp[x_col], imp[y_col],
                        facecolors=color, s=60, edgecolors="#2E7D32", linewidths=1.5, zorder=3,
                    )
                    for _, r in imp.iterrows():
                        texts.append(ax.text(
                            r[x_col], r[y_col], str(r[label_col]).upper(),
                            color="#333333", zorder=10,
                        ))

        if highlight_proteins and label_col in df.columns:
            wanted = {p.lower() for p in highlight_proteins}
            picked = df[df[label_col].str.lower().isin(wanted)]
            if not picked.empty:
                ax.scatter(
                    picked[x_col], picked[y_col],
                    facecolors=color, s=60, edgecolors="#1565C0", linewidths=1.5, zorder=4,
                )
                for _, r in picked.iterrows():
                    texts.append(ax.text(
                        r[x_col], r[y_col], str(r[label_col]).upper(),
                        color="#333333", zorder=10,
                    ))

        # Label every point (when no delta highlighting)
        if label_all and delta_col is None:
            for _, r in df.iterrows():
                texts.append(ax.text(
                    r[x_col], r[y_col], str(r[label_col]).upper(),
                    color="#333333", zorder=10,
                ))

    # Drop NaNs for regression / correlation (polyfit SVD fails on NaN)
    _finite = np.isfinite(x) & np.isfinite(y_plot)
    _n_dropped = int((~_finite).sum())
    if _n_dropped > 0:
        warnings.warn(
            f"make_scatter_comparison({display_name}): dropping {_n_dropped} of "
            f"{len(df)} rows with NaN in '{x_col}' or '{y_col}' before regression.",
            stacklevel=2,
        )
    _xf, _yf = x[_finite], y_plot[_finite]

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
            transform=ax.transAxes, ha="right", va="bottom",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#CCCCCC", alpha=0.8),
        )

    ax.set_ylabel(y_label_used)
    ax.set_title(title, pad=14)

    if ax_kde is not None:
        plt.setp(ax.get_xticklabels(), visible=False)
        ax.tick_params(axis="x", length=0)
        # Marginal KDE always operates on the baseline x-axis with raw
        # (un-residualised) baseline and guided values; this stays valid in
        # residual mode because both metrics share the [0, 1] axis.
        sns.kdeplot(x_raw, ax=ax_kde, color="#888888", fill=True, alpha=0.3,
                    label="Baseline", lw=1.5, bw_adjust=0.5)
        sns.kdeplot(y_raw, ax=ax_kde, color=color, fill=True, alpha=0.3,
                    label="Guided", lw=1.5, bw_adjust=0.5)
        ax_kde.set_xlabel(x_label)
        ax_kde.set_ylabel("Density")
        ax_kde.set_yticks([])
        ax_kde.legend(loc="upper right", frameon=True, framealpha=0.8, edgecolor="none")
    else:
        ax.set_xlabel(x_label)

    _fixed_limits = xlim is not None or ylim is not None
    if equal_aspect:
        ax.set_aspect("equal", adjustable="box" if _fixed_limits else "datalim")
    if equal_aspect and draw_diagonal_used and not _fixed_limits:
        pad = 0.01
        ax.set_xlim(lo + pad, hi)
        ax.set_ylim(lo + pad, hi)
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)

    if residual:
        # X is the baseline RSCC range (pinned to [0, 1]); y has its upper
        # bound pinned at 1 to expose the y = 1 - x ceiling, and a tight
        # lower bound so decreased points (y < 0) remain visible.
        if xlim is None:
            ax.set_xlim(0.0, 1.0)
        if ylim is None:
            _yvals = y_plot[np.isfinite(y_plot)]
            if len(_yvals) > 0:
                _y_lo = float(min(0.0, np.nanmin(_yvals))) - 0.02
            else:
                _y_lo = -0.05
            ax.set_ylim(_y_lo, 1.0)

    if texts:
        try:
            _n = len(texts)
            if _n > 15:
                adjust_text(
                    texts, ax=ax,
                    arrowprops=dict(arrowstyle="-", color="#888888", lw=0.5),
                    expand=(1.5, 1.8),
                    force_text=(1.0, 1.5),
                    force_static=(0.8, 1.2),
                    force_explode=(0.5, 0.8),
                    ensure_inside_axes=False,
                    max_move=(40, 40),
                    time_lim=10,
                )
            else:
                adjust_text(
                    texts, ax=ax,
                    arrowprops=dict(arrowstyle="-", color="#888888", lw=0.5),
                    expand=(1.3, 1.5),
                    force_text=(0.5, 0.8),
                    force_static=(0.5, 0.8),
                    force_explode=(0.3, 0.6),
                    ensure_inside_axes=True,
                    max_move=(20, 20),
                )
        except Exception as exc:
            warnings.warn(
                f"make_scatter_comparison({display_name}): adjust_text failed "
                f"({exc}); dropping labels rather than overlapping.",
                stacklevel=2,
            )
            for t in texts:
                t.remove()

    fig.tight_layout()
    return fig


def make_rscc_scatter(
    df: pd.DataFrame,
    model_name: str,
    optimal_gw: float,
) -> plt.Figure:
    """Per-model residual plot: baseline RSCC vs improvement (guided − baseline).

    ``df`` must contain a ``label`` column formatted ``"<PDB> <range>"``
    as produced by the RSCC scatter cell in ``02_fit_quality.py``. Four
    label sets are drawn: top-3 most improved, every decreased point, the
    5 points furthest below the y = 1 − x ceiling, and the manually
    highlighted ``6B8X 181-183`` and ``5MHX 162-164`` selections.
    """
    return make_scatter_comparison(
        df,
        model_name,
        x_col="rscc_baseline",
        y_col="rscc_optimal",
        x_label="Baseline RSCC",
        y_label="Improvement (Guided − Baseline RSCC)",
        title=f"{MODEL_DISPLAY.get(model_name, model_name)}: RSCC Improvement Over Baseline",
        delta_col="delta_rscc",
        label_col="label",
        show_marginal_kde=True,
        xlim=(0.0, 1.0),
        residual=True,
        label_top_n=3,
        label_far_from_ceiling_n=5,
        highlight_proteins=["6B8X 181-183", "5MHX 162-164"],
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
    draw_diagonal: bool = True,  # ignored — residual mode owns the geometry
    show_correlation: bool = False,
    equal_aspect: bool = True,  # ignored
) -> plt.Figure:
    """Per-class small-multiples residual plot.

    One subplot per classification, sharing x and y axes. In each subplot,
    the full pooled dataset is plotted in light gray as a context shadow,
    and only the points belonging to that subplot's class are overlaid in
    the class colour. Each subplot title is ``"<class> (n=<count>)"``.
    The y-axis is the residual ``y_col − x_col`` (improvement over baseline).
    Every subplot shows the y = 0 reference and the y = 1 − x ceiling.
    """
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

    # X is pinned to [0, 1]; y upper bound pinned at 1 (the y = 1 - x
    # ceiling), lower bound tight to the data so decreased points stay
    # visible.
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
        # Light-gray context shadow: full pooled dataset
        ax_i.scatter(
            x_all, y_all_resid,
            color="lightgray", alpha=0.25,
            s=14, edgecolors="none", zorder=1,
        )
        # Class points on top
        sub_x = df.loc[mask, x_col].to_numpy()
        sub_y = (df.loc[mask, y_col] - df.loc[mask, x_col]).to_numpy()
        color = CLASSIFICATION_PALETTE[cls]
        ax_i.scatter(
            sub_x, sub_y,
            c=color, alpha=0.7,
            s=18, edgecolors="none", zorder=3,
        )
        # Reference geometry
        ax_i.axhline(0.0, ls="--", color="#AAAAAA", lw=1, zorder=2)
        ax_i.plot(_xs, 1.0 - _xs, ls="--", color="#AAAAAA", lw=0.8, zorder=2)
        ax_i.set_title(
            f"{CLASSIFICATION_DISPLAY[cls]} (n={int(mask.sum())})", pad=10,
        )
        ax_i.set_xlabel(x_label)

    # Hide any unused axes (when classes don't fill the grid)
    for ax_i in flat_axes[len(classes):]:
        ax_i.set_visible(False)

    # Single y-label on the leftmost axis of each row
    for row_idx in range(n_rows):
        flat_axes[row_idx * n_cols].set_ylabel(y_label)

    # Apply shared limits
    flat_axes[0].set_xlim(x_lo, x_hi)
    flat_axes[0].set_ylim(y_lo, y_hi)

    # Reserve room above the subplots so the suptitle does not sit flush
    # against the per-axes title.
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
    """Combined per-class residual scatter for many models on one figure.

    Rows are models, columns are classifications. All subplots share a
    single x and y axis pair (one figure-level x label, one figure-level
    y label). Class column titles with per-class ``n`` appear only on the
    top row; each row is preceded by a horizontal banner with the model
    display name.

    The per-class ``n`` is the count of unique ``(protein, selection)``
    pairs in that class — independent of which model is in the row, since
    every (protein, selection) maps to exactly one classification.

    The light-gray context shadow in each subplot is that row's model's
    pooled residual cloud (matching :func:`make_classified_scatter`).
    """
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

    # n per class is the count of unique (protein, selection) pairs in
    # that class — invariant across rows.
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
    x_label: str = "Baseline RSCC (mean \u00b1 SD)",
    y_label: str = "Improvement (Guided \u2212 Baseline RSCC)",
    title: str = "RSCC Improvement by Conformation Type",
    draw_diagonal: bool = True,  # ignored \u2014 residual mode owns the geometry
    equal_aspect: bool = True,  # ignored
    show_correlation: bool = False,
    show_errors: bool = True,
    label_col: str | None = None,
) -> plt.Figure:
    """Per-class small-multiples residual plot with propagated error bars.

    Each row of *agg_df* is one (protein, selection) pair aggregated across
    models. Per-axis error columns are mean \u00b1 SD; y-axis residual error is
    propagated as ``sqrt(sd_g^2 + sd_b^2)`` assuming independence between
    the baseline and guided per-protein estimates.

    The figure is a row of subplots, one per classification. Each subplot
    shows the full pooled dataset in light gray as a context shadow, then
    overlays only the points belonging to that subplot's class in the
    class colour. Subplots share x and y axes; titles are
    ``"<class> (n=<count>)"``. Every subplot shows the y = 0 reference
    and the y = 1 \u2212 x ceiling envelope.
    """
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

    # X pinned to [0, 1]; y upper bound pinned at 1, lower bound tight to
    # data so decreased points stay visible.
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
        # Light-gray context shadow (no error bars; just point cloud)
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


def make_aggregated_scatter(
    agg_df: pd.DataFrame,
    *,
    x_col: str,
    y_col: str,
    x_q1_col: str | None = None,
    x_q3_col: str | None = None,
    y_q1_col: str | None = None,
    y_q3_col: str | None = None,
    x_label: str,
    y_label: str,
    title: str,
    label_col: str = "protein",
    draw_diagonal: bool = False,
    equal_aspect: bool = False,
    point_color: str = "#4D4D4D",
    label_points: bool = False,
    highlight_top_n: int = 3,
    highlight_decreased: bool | int = True,
    highlight_proteins: list[str] | None = None,
    show_correlation: bool = False,
    show_marginal_kde: bool = False,
    show_errors: bool = True,
    xlim: tuple[float, float] | None = None,
    ylim: tuple[float, float] | None = None,
    residual: bool = False,
    label_top_n: int = 3,
    label_bottom_n: int | None = None,
    label_far_from_ceiling_n: int = 5,
) -> plt.Figure:
    """Scatter plot with asymmetric IQR error bars for cross-model aggregations.

    Each row of *agg_df* becomes one point; ``x_col``/``y_col`` are the medians
    and ``*_q1``/``*_q3`` provide the IQR whiskers. Error-bar styling mirrors
    :func:`make_aggregated_spiral_figure` so the two figure families stay
    visually consistent.

    Pass ``show_errors=False`` (or omit the ``*_q1_col``/``*_q3_col``
    arguments) to draw only the central points; per-model spread is shown
    elsewhere in this codebase by per-model scatter figures, so collapsed
    aggregates don't always need whiskers.

    ``xlim``/``ylim`` override the auto-computed axis ranges (including the
    range derived from ``equal_aspect + draw_diagonal``) and also extend the
    diagonal so it spans the full fixed range.

    ``highlight_decreased`` accepts either a bool (``True`` = highlight every
    decreased point, the default) or an integer count (``n`` = only label the
    ``n`` most decreased rows). Use a bounded count on per-selection plots
    where the raw decreased count is too high to be readable.

    When ``residual=True``, the y-axis is ``y_col − x_col`` (improvement over
    baseline). The y = x diagonal is replaced by a y = 0 horizontal reference
    and a y = 1 − x ceiling envelope. Y error bars are propagated as
    sqrt(sd_g^2 + sd_b^2), assuming the baseline and guided per-protein
    estimates are independent — this is an approximation: in practice the two
    share most of the same models, so the true residual SD is somewhat smaller.
    Marker styling switches to small, un-edged, alpha=0.3 dots. In residual
    mode four label sets are drawn: (a) the ``label_top_n`` most improved
    (green ring), (b) decreased points (red ring) — every point with
    delta < 0 by default, or only the ``label_bottom_n`` most decreased when
    that argument is given, (c) the ``label_far_from_ceiling_n`` points
    furthest below the y = 1 − x ceiling (orange ring), and (d) any point
    whose ``label_col`` value matches an entry in ``highlight_proteins``
    (blue ring overlay, stacked on top of the other rings).
    """
    if show_marginal_kde:
        fig = plt.figure(figsize=(6, 7.5))
        gs = fig.add_gridspec(2, 1, height_ratios=[4, 1], hspace=0.05)
        ax = fig.add_subplot(gs[0])
        ax_kde = fig.add_subplot(gs[1], sharex=ax)
    else:
        fig, ax = plt.subplots(figsize=(6, 6))
        ax_kde = None

    x_med = agg_df[x_col].to_numpy()
    y_med_raw = agg_df[y_col].to_numpy()

    _err_cols_present = all(
        c is not None for c in (x_q1_col, x_q3_col, y_q1_col, y_q3_col)
    )
    draw_errors = show_errors and _err_cols_present

    if draw_errors:
        x_q1 = agg_df[x_q1_col].to_numpy()
        x_q3 = agg_df[x_q3_col].to_numpy()
        y_q1 = agg_df[y_q1_col].to_numpy()
        y_q3 = agg_df[y_q3_col].to_numpy()

    if residual:
        y_med = y_med_raw - x_med
        equal_aspect = False
        if "Improvement" not in y_label:
            y_label_used = "Improvement (Guided − Baseline RSCC)"
        else:
            y_label_used = y_label
        if draw_errors:
            # propagate independent errors as sqrt(sd_g^2 + sd_b^2). q1/q3 are
            # mean ± SD here (see notebook aggregation), so SD = (q3 − q1) / 2
            # per axis.
            sd_x = (x_q3 - x_q1) / 2.0
            sd_y = (y_q3 - y_q1) / 2.0
            sd_resid = np.sqrt(sd_y**2 + sd_x**2)
            xerr = np.vstack([x_med - x_q1, x_q3 - x_med])
            yerr = np.vstack([sd_resid, sd_resid])
    else:
        y_med = y_med_raw
        y_label_used = y_label
        if draw_errors:
            xerr = np.vstack([x_med - x_q1, x_q3 - x_med])
            yerr = np.vstack([y_med - y_q1, y_q3 - y_med])

    if draw_diagonal and not residual:
        if draw_errors:
            lo = float(min(np.nanmin(x_q1), np.nanmin(y_q1))) - 0.02
            hi = float(max(np.nanmax(x_q3), np.nanmax(y_q3))) + 0.02
        else:
            lo = float(min(np.nanmin(x_med), np.nanmin(y_med))) - 0.02
            hi = float(max(np.nanmax(x_med), np.nanmax(y_med))) + 0.02
        if xlim is not None:
            lo = min(lo, xlim[0])
            hi = max(hi, xlim[1])
        if ylim is not None:
            lo = min(lo, ylim[0])
            hi = max(hi, ylim[1])
        ax.plot([lo, hi], [lo, hi], ls="--", color="#AAAAAA", lw=1, zorder=1)

    if residual:
        # y = 0 horizontal reference and y = 1 - x envelope
        _x_lo = float(np.nanmin(x_med)) if xlim is None else xlim[0]
        _x_hi = float(np.nanmax(x_med)) if xlim is None else xlim[1]
        _x_lo = min(_x_lo, 0.0) - 0.02
        _x_hi = max(_x_hi, 1.0) + 0.02
        ax.axhline(0.0, ls="--", color="#AAAAAA", lw=1, zorder=1)
        _xs = np.linspace(_x_lo, _x_hi, 100)
        ax.plot(_xs, 1.0 - _xs, ls="--", color="#AAAAAA", lw=0.8, zorder=1)

    if draw_errors:
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
    if residual:
        ax.scatter(
            x_med, y_med,
            c=point_color, s=18, alpha=0.3,
            edgecolors="none", zorder=4,
        )
    else:
        ax.scatter(
            x_med, y_med,
            c=point_color, s=55,
            edgecolors="white", linewidths=0.5, zorder=4,
        )

    texts: list = []
    # Delta is always y_raw - x_med, irrespective of mode (this drives ranking)
    _delta = agg_df[y_col] - agg_df[x_col]
    # Helper to map x-y back to plot space when placing rings/labels.
    def _yp(idx):
        if residual:
            return agg_df.loc[idx, y_col] - agg_df.loc[idx, x_col]
        return agg_df.loc[idx, y_col]

    if residual:
        # Residual mode label sets:
        #   1. top-N most improved (delta > 0)            — green ring
        #   2. every point with delta < 0                 — red ring
        #   3. N points furthest below the y = 1 − x line — orange ring
        #   4. highlight_proteins                          — blue ring overlay
        # Distance below ceiling = 1 − y_col (large = low guided value).
        _seen_idx: set = set()

        if label_top_n > 0:
            imp_idx = _delta.nlargest(label_top_n).index
            imp_idx = [i for i in imp_idx if _delta.loc[i] > 0]
            if imp_idx:
                imp = agg_df.loc[imp_idx]
                ax.scatter(
                    imp[x_col], _yp(imp.index),
                    facecolors=point_color, s=32,
                    edgecolors="#2E7D32", linewidths=1.2, zorder=5,
                )
                if label_col in agg_df.columns:
                    for i in imp_idx:
                        _seen_idx.add(i)
                        texts.append(ax.text(
                            agg_df.loc[i, x_col], _yp(i),
                            str(agg_df.loc[i, label_col]).upper(),
                            color="#333333", zorder=10,
                        ))

        dec_mask = _delta < 0
        if label_bottom_n is not None and label_bottom_n > 0:
            dec_ranked = _delta[dec_mask].nsmallest(label_bottom_n).index
            dec_idx = [i for i in dec_ranked if i not in _seen_idx]
        else:
            dec_idx = [i for i in agg_df.index[dec_mask] if i not in _seen_idx]
        if dec_idx:
            dec = agg_df.loc[dec_idx]
            ax.scatter(
                dec[x_col], _yp(dec.index),
                facecolors=point_color, s=32,
                edgecolors="#D32F2F", linewidths=1.2, zorder=5,
            )
            if label_col in agg_df.columns:
                for i in dec_idx:
                    _seen_idx.add(i)
                    texts.append(ax.text(
                        agg_df.loc[i, x_col], _yp(i),
                        str(agg_df.loc[i, label_col]).upper(),
                        color="#333333", zorder=10,
                    ))

        if label_far_from_ceiling_n > 0:
            below_dist = 1.0 - agg_df[y_col]
            far_idx = below_dist.nlargest(label_far_from_ceiling_n).index
            far_idx = [i for i in far_idx if below_dist.loc[i] > 0 and i not in _seen_idx]
            if far_idx:
                far = agg_df.loc[far_idx]
                ax.scatter(
                    far[x_col], _yp(far.index),
                    facecolors=point_color, s=32,
                    edgecolors="#F57C00", linewidths=1.2, zorder=5,
                )
                if label_col in agg_df.columns:
                    for i in far_idx:
                        _seen_idx.add(i)
                        texts.append(ax.text(
                            agg_df.loc[i, x_col], _yp(i),
                            str(agg_df.loc[i, label_col]).upper(),
                            color="#333333", zorder=10,
                        ))

        if highlight_proteins and label_col in agg_df.columns:
            # Blue ring overlay; stacks on top of any green/red/orange ring.
            wanted = {p.lower() for p in highlight_proteins}
            picked_idx = [
                i for i in agg_df.index
                if str(agg_df.loc[i, label_col]).lower() in wanted
            ]
            if picked_idx:
                picked = agg_df.loc[picked_idx]
                ax.scatter(
                    picked[x_col], _yp(picked.index),
                    facecolors=point_color, s=32,
                    edgecolors="#1565C0", linewidths=1.2, zorder=6,
                )
                for i in picked_idx:
                    if i in _seen_idx:
                        continue
                    _seen_idx.add(i)
                    texts.append(ax.text(
                        agg_df.loc[i, x_col], _yp(i),
                        str(agg_df.loc[i, label_col]).upper(),
                        color="#333333", zorder=10,
                    ))

    else:
        # Legacy non-residual highlighting
        if highlight_decreased and label_col in agg_df.columns:
            dec_mask = _delta < 0
            dec = agg_df[dec_mask]
            if not isinstance(highlight_decreased, bool):
                dec_delta = _delta[dec_mask]
                dec = dec.loc[dec_delta.nsmallest(int(highlight_decreased)).index]
            if not dec.empty:
                ax.scatter(
                    dec[x_col], dec[y_col],
                    facecolors=point_color, s=70,
                    edgecolors="#D32F2F", linewidths=1.5, zorder=5,
                )
                for _, r in dec.iterrows():
                    texts.append(ax.text(
                        r[x_col], r[y_col], str(r[label_col]).upper(),
                        color="#333333", zorder=10,
                    ))

        if highlight_top_n > 0 and label_col in agg_df.columns:
            imp = agg_df.loc[_delta.nlargest(highlight_top_n).index]
            imp = imp[_delta.loc[imp.index] > 0]
            if not imp.empty:
                ax.scatter(
                    imp[x_col], imp[y_col],
                    facecolors=point_color, s=70,
                    edgecolors="#2E7D32", linewidths=1.5, zorder=5,
                )
                for _, r in imp.iterrows():
                    texts.append(ax.text(
                        r[x_col], r[y_col], str(r[label_col]).upper(),
                        color="#333333", zorder=10,
                    ))

        if highlight_proteins and label_col in agg_df.columns:
            wanted = {p.lower() for p in highlight_proteins}
            picked = agg_df[agg_df[label_col].str.lower().isin(wanted)]
            already = {t.get_text().lower() for t in texts}
            picked = picked[~picked[label_col].str.upper().str.lower().isin(already)]
            if not picked.empty:
                ax.scatter(
                    picked[x_col], picked[y_col],
                    facecolors=point_color, s=70,
                    edgecolors="#1565C0", linewidths=1.5, zorder=5,
                )
                for _, r in picked.iterrows():
                    texts.append(ax.text(
                        r[x_col], r[y_col], str(r[label_col]).upper(),
                        color="#333333", zorder=10,
                    ))

    if show_correlation:
        _finite = np.isfinite(x_med) & np.isfinite(y_med)
        if int(_finite.sum()) >= 5:
            r_val = float(np.corrcoef(x_med[_finite], y_med[_finite])[0, 1])
            ax.text(
                0.97, 0.03, f"r = {r_val:.3f}",
                transform=ax.transAxes, ha="right", va="bottom",
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#CCCCCC", alpha=0.8),
            )

    ax.set_ylabel(y_label_used)
    ax.set_title(title, pad=14)

    if ax_kde is not None:
        plt.setp(ax.get_xticklabels(), visible=False)
        ax.tick_params(axis="x", length=0)
        # Marginal panel always shows raw baseline/guided distributions on
        # the shared baseline x-axis, even in residual mode.
        sns.kdeplot(x_med, ax=ax_kde, color="#888888", fill=True, alpha=0.3,
                    label="Baseline", lw=1.5, bw_adjust=0.5)
        sns.kdeplot(y_med_raw, ax=ax_kde, color=point_color, fill=True, alpha=0.3,
                    label="Guided", lw=1.5, bw_adjust=0.5)
        ax_kde.set_xlabel(x_label)
        ax_kde.set_ylabel("Density")
        ax_kde.set_yticks([])
        ax_kde.legend(loc="upper right", frameon=True, framealpha=0.8, edgecolor="none")
    else:
        ax.set_xlabel(x_label)

    # When fixed xlim/ylim are supplied, use adjustable="box" so the axes box
    # is reshaped to honor the requested data ranges instead of matplotlib
    # overriding one axis to preserve aspect ratio.
    _fixed_limits = xlim is not None or ylim is not None
    if equal_aspect:
        ax.set_aspect("equal", adjustable="box" if _fixed_limits else "datalim")
    if equal_aspect and draw_diagonal and not residual and not _fixed_limits:
        pad = 0.01
        ax.set_xlim(lo + pad, hi)
        ax.set_ylim(lo + pad, hi)
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)

    if residual:
        # X pinned to [0, 1]; y upper bound pinned at 1 with a tight lower
        # bound. See note in make_scatter_comparison.
        if xlim is None:
            ax.set_xlim(0.0, 1.0)
        if ylim is None:
            _yvals = y_med[np.isfinite(y_med)]
            if len(_yvals) > 0:
                _y_lo = float(min(0.0, np.nanmin(_yvals))) - 0.02
            else:
                _y_lo = -0.05
            ax.set_ylim(_y_lo, 1.0)

    if texts:
        try:
            adjust_text(
                texts, ax=ax,
                arrowprops=dict(arrowstyle="-", color="#888888", lw=0.5),
                expand=(1.3, 1.5),
                force_text=(0.5, 0.8),
                force_static=(0.5, 0.8),
                force_explode=(0.3, 0.6),
                ensure_inside_axes=True,
                max_move=(20, 20),
            )
        except Exception as exc:
            warnings.warn(
                f"make_aggregated_scatter: adjust_text failed ({exc}); "
                "dropping labels rather than overlapping.",
                stacklevel=2,
            )
            for t in texts:
                t.remove()

    fig.tight_layout()
    return fig


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
    """Half-strength log-scaled per-point alpha used by the calibration figures.

    Parameters
    ----------
    s_values
        Cluster sample sizes (``s_cluster`` values) for each point.
    s_max
        Maximum ``s_cluster`` taken across the union of the baseline and
        optimal-guidance prepared frames so alpha is comparable across panels
        and across both calibration figures.

    Returns
    -------
    np.ndarray
        Per-point alpha in ``[0, 1]``. Returns zeros when ``s_max <= 0``
        (degenerate, should not occur post-filter).
    """
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
    """Per-bin mean and sample std of ``y`` in ``n_bins`` equal-count bins of ``x``.

    Bin edges are quantiles of ``x`` so each bin contains roughly the same
    number of points; duplicate edges (from ties at boundary values) are
    collapsed via ``np.unique`` so the effective bin count can be less than
    ``n_bins``. The dot for each bin is plotted at the within-bin mean of
    ``x`` (data-weighted) rather than the bin-edge midpoint, so the curve
    traces the support of the data instead of an arbitrary bin grid. Bins
    with fewer than ``min_count`` points return ``np.nan`` in both ``means``
    and ``stds`` so callers can pass ``(centers, means)`` straight to
    ``ax.plot`` and ``(centers, means - stds, means + stds)`` to
    ``ax.fill_between``; matplotlib breaks the line / skips the band at NaN.
    """
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
    """1x4 calibration overlay of Baseline and Optimal guidance per model.

    One panel per main model. Baseline (gray, ``#888888``) and Optimal
    guidance (dark, ``#333333``) points are overlaid in the same axes, each
    with a binned conditional-mean curve. Reference lines ``y = x`` (training
    prior) and ``y = 0.5`` (guidance target) are drawn behind the points.
    Each scatter axis has a marginal y-histogram on its right edge that
    overlays both conditions' ``b_occupancy`` distributions. Slope of
    ``b_occupancy`` on ``frac_b_cluster`` and Pearson r are annotated for the
    Optimal-guidance condition only.

    Parameters
    ----------
    baseline_df, optimal_df
        Outputs of ``build_training_bias_df`` for the baseline (gw=0.0) and
        optimal-guidance conditions, respectively.
    model_order
        Optional explicit model ordering. Defaults to
        ``_TRAINING_BIAS_MAIN_MODELS``.
    jitter
        Symmetric uniform jitter half-width applied to x. Defaults to 0.01.
    seed
        Seed for the jitter RNG. Defaults to 0.
    max_cluster_size
        If set, restrict both conditions to rows with
        ``n_cluster_evaluated <= max_cluster_size`` (sparse-training-prior subset).
    cutoff_status
        If set, restrict both conditions to rows with the given
        ``cutoff_status`` (``"Pre"`` or ``"Post"``).
    n_bins, min_bin_count
        Forwarded to ``_training_bias_binned_means`` for both conditions'
        binned conditional-mean curves. Lower these for small subsets so a
        curve can still be drawn (e.g. ``n_bins=3, min_bin_count=2`` for
        ~9-point subsets).
    """
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
    """1x4 strip plot for selections with no PDB neighbors in the training cluster.

    Restricted to ``(protein, selection)`` rows where
    ``n_assigned_to_A == 0`` and ``n_assigned_to_B == 0`` — i.e., no
    pre-cutoff PDB cluster member is assigned to either altloc by lDDT
    comparison. The calibration X-axis (``frac_b_cluster``) is undefined
    for these rows, so they cannot appear on Figure 34; this figure shows
    their predicted B-occupancy distribution per model instead.

    Layout mirrors ``make_training_bias_calibration_overlay``: 1x4 row of
    panels, one per main model. Each panel overlays a jittered strip of
    Baseline (gray) and Optimal-guidance (model-colored) points, with a
    marginal y-histogram on the right edge and a ``y = 0.5`` reference line.
    Per-panel sample sizes are annotated.
    """
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
    """Standalone legend figure for ``make_training_bias_calibration_overlay``.

    Two-row vertical legend: ``Baseline`` (single gray marker) and
    ``Optimal guidance`` (four model-colored markers shown side-by-side via
    ``HandlerTuple``). Reference lines (``y = x``, ``y = 0.5``) are not
    included; they remain on the plot panels themselves.
    """
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
    """2x4 calibrated small multiples with per-panel mean curves and marginals.

    Mirrors the geometry of ``make_training_bias_scatter`` (rows = condition,
    columns = model) but adds, per panel: both reference lines (y = x training
    prior, y = 0.5 guidance target), a binned conditional-mean curve, an OLS
    slope and Pearson r annotation, and a single-condition marginal
    y-histogram strip on the right edge.

    Parameters
    ----------
    baseline_df, optimal_df
        Outputs of ``build_training_bias_df`` for the baseline (gw=0.0) and
        optimal-guidance conditions, respectively.
    model_order
        Optional explicit model ordering. Defaults to
        ``_TRAINING_BIAS_MAIN_MODELS``.
    jitter
        Symmetric uniform jitter half-width applied to x. Defaults to 0.01.
    seed
        Seed for the jitter RNG. Defaults to 0.
    """
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
    """Single empty panel with the y=x and y=0.5 references and their bisector.

    Intentionally minimal explanatory companion to the v2 small multiples:
    same axes (``xlim = ylim = (-0.02, 1.02)``, ticks at ``[0, 0.5, 1]``),
    ``y = x`` (training prior) and ``y = 0.5`` (guidance target) drawn as
    light dashed reference lines, and a solid line bisecting them
    (``y = 0.5 * x + 0.25``, slope 0.5 through the intersection at
    ``(0.5, 0.5)``) for direct annotation in downstream editing.
    """
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


def build_training_bias_df(
    lddt_df: pd.DataFrame,
    cluster_df: pd.DataFrame,
    cutoffs_df: pd.DataFrame,
    model_to_gw: dict[str, float],
    *,
    verbose: bool = True,
) -> pd.DataFrame:
    """Compose LDDT, training-cluster, and cutoff frames into one tidy frame.

    For each ``(model_name, guidance_weight)`` pair in *model_to_gw*, slices
    *lddt_df* at exactly that condition, parses the B occupancy from the
    ``occupancies`` string, inner-joins on (``protein``, ``selection``) with
    *cluster_df* to attach pre-cutoff cluster counts, and left-joins on
    ``protein`` with *cutoffs_df* to attach the Pre/Post deposition flag.

    Warnings (``verbose=True``) are printed for: NaN occupancies, lddt rows
    with no matching cluster row, and cluster rows where ``n_cluster_evaluated
    > 0`` but ``A == B == 0`` (X is undefined for plot 1, eval > 0
    disqualifies plot 2).

    Returns
    -------
    pd.DataFrame
        Columns: ``protein``, ``selection``, ``model_name``,
        ``guidance_weight``, ``b_occupancy``, ``n_assigned_to_A``,
        ``n_assigned_to_B``, ``n_cluster_evaluated``,
        ``frac_a_cluster``, ``x_one_minus_fa``, ``cutoff_status``.
    """
    pieces: list[pd.DataFrame] = []
    for model_name, gw in model_to_gw.items():
        sub = lddt_df.query(
            "model_name == @model_name and guidance_weight == @gw"
        ).copy()
        if sub.empty:
            if verbose:
                print(
                    f"⚠️  No lddt rows for model={model_name!r} gw={gw}; "
                    "skipping in build_training_bias_df."
                )
            continue
        sub["b_occupancy"] = sub["occupancies"].map(_parse_b_occupancy)
        nan_mask = sub["b_occupancy"].isna()
        if verbose and nan_mask.any():
            n_nan = int(nan_mask.sum())
            sample = (
                sub.loc[nan_mask, ["protein", "model_name", "guidance_weight", "selection"]]
                .head(min(n_nan, 50))
                .to_dict(orient="records")
            )
            print(
                f"⚠️  {n_nan} NaN B-occupancy row(s) for model={model_name!r} "
                f"gw={gw}; dropped. Listing up to 50:\n  " + "\n  ".join(repr(r) for r in sample)
            )
        sub = sub.loc[~nan_mask, ["protein", "selection", "model_name", "guidance_weight", "b_occupancy"]]
        pieces.append(sub)
    if not pieces:
        return pd.DataFrame(columns=[
            "protein", "selection", "model_name", "guidance_weight", "b_occupancy",
            "n_assigned_to_A", "n_assigned_to_B", "n_cluster_evaluated",
            "frac_a_cluster", "x_one_minus_fa", "cutoff_status",
        ])
    lddt_subset = pd.concat(pieces, ignore_index=True)

    if verbose:
        unmatched = (
            lddt_subset.merge(
                cluster_df[["protein", "selection"]].assign(_in_cluster=1),
                on=["protein", "selection"], how="left",
            )
            .pipe(lambda d: d[d["_in_cluster"].isna()])
            [["protein", "selection"]]
            .drop_duplicates()
        )
        if not unmatched.empty:
            print(
                f"⚠️  {len(unmatched)} (protein, selection) key(s) in lddt have no "
                f"training-cluster row; dropped. Examples: "
                f"{unmatched.head(3).to_dict(orient='records')}"
            )
        ambiguous = cluster_df[
            (cluster_df["n_cluster_evaluated"] > 0)
            & (cluster_df["n_assigned_to_A"] == 0)
            & (cluster_df["n_assigned_to_B"] == 0)
        ]
        if not ambiguous.empty:
            print(
                f"⚠️  {len(ambiguous)} cluster row(s) have eval > 0 but A == B == 0 "
                "(X undefined, eval > 0 disqualifies the empty-cluster plot); dropped."
            )

    merged = lddt_subset.merge(
        cluster_df, on=["protein", "selection"], how="inner",
    )
    merged = merged.merge(cutoffs_df[["protein", "cutoff_status"]], on="protein", how="left")
    return merged[[
        "protein", "selection", "model_name", "guidance_weight", "b_occupancy",
        "n_assigned_to_A", "n_assigned_to_B", "n_cluster_evaluated",
        "frac_a_cluster", "x_one_minus_fa", "cutoff_status",
    ]]


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
            ("delta_silhouette_to_ref", "Δ Clustering Score vs baseline"),
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
                ax.set_title(MODEL_DISPLAY.get(model, model))
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
        fig.suptitle(title, y=1.08)
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
    optimal_gw : dict[str, float] or None
        Mapping of ``model_name`` → optimal guidance weight. If provided, the
        point corresponding to each model's optimal guidance weight is drawn
        with a thicker black outline.

    Returns
    -------
    matplotlib.figure.Figure
    """
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

        # Scale dot size by guidance weight (baseline omitted here, so use
        # gws directly). Map [min_gw, max_gw] → [s_min, s_max].
        s_min, s_max = 30.0, 220.0
        gw_hi = float(max(gws)) if len(gws) else 1.0
        gw_lo = float(min(gws)) if len(gws) else 0.0
        if gw_hi > gw_lo:
            sizes = s_min + (gws - gw_lo) / (gw_hi - gw_lo) * (s_max - s_min)
        else:
            sizes = np.full_like(gws, (s_min + s_max) / 2.0, dtype=float)

        # Per-point outline: thick black ring on the optimal guidance weight,
        # thin white ring everywhere else.
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
                    fc="white",
                    ec="none",
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
        [0], [0], marker="o", color="w",
        markerfacecolor="#CCCCCC", markeredgecolor="#111111", markeredgewidth=2.0,
        markersize=10, label="Optimal guidance",
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
