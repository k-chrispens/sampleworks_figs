"""Data transformations: delta frames, training-bias join, group filters, slugs.

Pure DataFrame -> DataFrame logic. No matplotlib imports.
"""

from __future__ import annotations

import re

import numpy as np
import pandas as pd

from loaders import _parse_b_occupancy
from style import (
    MODEL_DISPLAY,
    MODELS_ALT,
    MODELS_PROTENIX_GROUP,
    MODELS_RF3_GROUP,
)


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


