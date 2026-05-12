"""CSV loaders and path-metadata tagging utilities.

All functions that read raw CSV files (RSCC, LDDT, classifications,
training-cluster counts, deposition cutoffs) plus the path-pattern
extraction helpers and metadata-column post-processing.
"""

from __future__ import annotations

import os
import re

import numpy as np
import pandas as pd

from style import (
    GW_TO_STRENGTH,
    MODEL_NAME_ALIAS,
    OCC_FILTER,
    OCC_FILTER_ALTLOC,
)

# Proteins missing from the deposition-cutoffs CSV that the analysis still
# references; all four are pre-2018 depositions, so they are filled in as "Pre".
_MISSING_CUTOFF_FILL = {"1vme": "Pre", "4ole": "Pre", "5sop": "Pre", "6b8x": "Pre"}

_PATH_PATTERN = re.compile(
    r"([0-9A-Za-z]{4})_([\d.]+occ[A-Z](?:_[\d.]+occ[A-Z])*)/(.*?)/(.*?)/ens8_gw([0-9\.]+)/"
)


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
