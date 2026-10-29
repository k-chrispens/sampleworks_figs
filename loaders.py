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


def load_rscc(
    csv_path: str,
    model_name_override: str | None = None,
) -> pd.DataFrame:
    """Load RSCC results with altloc filtering and model-name construction.

    Parameters
    ----------
    csv_path : str
        Path to the RSCC CSV file.
    model_name_override : str | None
        If given, force-set ``model_name`` to this value. Used for the
        protenix-mini / protenix-tiny / rf3-chiral-off max-RMSD CSVs whose
        ``model`` column reads ``"protenix"`` / ``"rf3"`` and is not otherwise
        distinguishable. When ``None``, model_name is built from the
        ``model``/``method`` columns and aliased.

    Returns
    -------
    pd.DataFrame
        RSCC DataFrame with ``model_name``, ``guidance_strength``, etc.
    """
    df = pd.read_csv(csv_path, dtype={"error": str}).fillna({"error": ""})
    _occ = OCC_FILTER_ALTLOC
    df = df.query("altloc_occupancies == @_occ").copy()
    if model_name_override is not None:
        df["model_name"] = model_name_override
    else:
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


def load_rscc_recyc_diff(csv_path: str) -> pd.DataFrame:
    """Load the recycles/diffusion-steps grid-search RSCC CSV.

    Same 50/50 altloc filter, model-name construction (``model``/``method`` aliased
    through ``MODEL_NAME_ALIAS``), and metadata columns as :func:`load_rscc`, but
    additionally keeps the ``n_recycles`` and ``n_diffusion_steps`` grid axes. Rows
    with a NaN ``rscc`` are dropped with a one-line warning.

    Parameters
    ----------
    csv_path : str
        Path to the recycles/diffusion grid-search RSCC CSV
        (e.g. ``rscc_results_5mhx_recyc_diff.csv``).

    Returns
    -------
    pd.DataFrame
        Columns: ``protein``, ``model_name``, ``ensemble_size``, ``guidance_weight``,
        ``n_recycles``, ``n_diffusion_steps``, ``rscc``, ``selection`` plus the
        ``selection_name`` / ``guidance_strength`` metadata columns.
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
            "n_recycles",
            "n_diffusion_steps",
            "rscc",
            "selection",
        ]
    ]
    _n_before = len(df)
    df = df[df.rscc.notna()].copy()
    if len(df) < _n_before:
        print(f"load_rscc_recyc_diff: dropped {_n_before - len(df)} NaN-rscc rows from {csv_path}")
    df = add_simplified_metadata_columns(df)
    return df


# Per-model-group max-RMSD-region RSCC CSVs that together form the primary RSCC
# frame for every non-classification figure. Each tuple is
# ``(csv_path, model_name_override)``. protenix-mini/tiny store
# ``model == "protenix"`` and rf3-chiral-off stores ``model == "rf3"``, so those
# need an override to disambiguate; boltz_protenix and rf3_fixed carry a
# distinguishable ``model`` column already. rf3-chiral-off has no max-RMSD
# variant, so its non-max CSV is used (the only RSCC source for that model).
RSCC_MAX_RMSD_SOURCES: list[tuple[str, str | None]] = [
    ("results/rscc_results_max_rmsd_boltz_protenix.csv", None),
    ("results/rscc_results_max_rmsd_rf3_fixed.csv", None),
    ("results/rscc_results_max_rmsd_protenix_mini.csv", "protenix-mini"),
    ("results/rscc_results_max_rmsd_protenix_tiny.csv", "protenix-tiny"),
    ("results/rscc_results_rf3_chiral_off.csv", "rf3-chiral-off"),
]


def load_rscc_max_rmsd() -> pd.DataFrame:
    """Load and concatenate the max-RMSD-region RSCC CSVs into one frame.

    This is the primary RSCC source for all non-classification figures
    (per-model + aggregated scatters, baseline distributions, the RSCC-vs-LDDT
    cross-metric scatter, and the spiral trajectories). The by-class grid
    instead concatenates the five base (non-max) ``results/rscc_results_*.csv``
    inline via :func:`load_rscc` (it needs the base's full per-region selection
    coverage; all four conformational classes, incl. domain-shift, flow through
    one containment merge).

    Each source CSV is loaded through :func:`load_rscc` (same 50/50 filter,
    model-name construction, and metadata columns), with a
    ``model_name_override`` where the raw ``model`` column is ambiguous
    (see :data:`RSCC_MAX_RMSD_SOURCES`). Output contract matches
    :func:`load_rscc`.
    """
    frames = [
        load_rscc(path, model_name_override=override)
        for path, override in RSCC_MAX_RMSD_SOURCES
    ]
    return pd.concat(frames, ignore_index=True)


def load_min_altloc_rmsd(
    csv_path: str,
    model_name_override: str | None = None,
    *,
    keep_failed: bool = False,
) -> pd.DataFrame:
    """Load per-(protein, model, gw, selection) minimum-RMSD-to-each-altloc results.

    Mirrors :func:`load_rscc` structurally: applies the
    :data:`OCC_FILTER_ALTLOC` 50/50 filter on ``altloc_occupancies``,
    constructs ``model_name`` from ``model + method`` for ``boltz2`` rows
    (otherwise from ``model``), and runs the result through
    :data:`MODEL_NAME_ALIAS`. Drops rows where either ``min_rmsd_to_A`` or
    ``min_rmsd_to_B`` is NaN, printing a one-line warning per CSV when any
    rows were dropped.

    Parameters
    ----------
    csv_path : str
        Path to a ``min_altloc_rmsd_results_*.csv``.
    model_name_override : str | None
        If given, force-set ``model_name`` to this value (used for the
        protenix-mini / protenix-tiny CSVs whose ``model`` column reads
        ``"protenix"`` and is not otherwise distinguishable from the
        canonical ``protenix`` rows). When ``None``, model_name is built
        from the ``model``/``method`` columns and aliased.
    keep_failed : bool
        Preserve rows with missing RMSDs for explicit coverage reporting.
        The default retains the original complete-pair-only behavior.

    Returns
    -------
    pd.DataFrame
        Columns: ``protein``, ``model_name``, ``ensemble_size``,
        ``guidance_weight``, ``selection``, ``min_rmsd_to_A``,
        ``min_rmsd_to_B``, plus ``selection_name`` and
        ``guidance_strength`` from :func:`add_simplified_metadata_columns`.
    """
    df = pd.read_csv(csv_path, dtype={"error": str}).fillna({"error": ""})
    _occ = OCC_FILTER_ALTLOC
    df = df.query("altloc_occupancies == @_occ").copy()
    if model_name_override is not None:
        df["model_name"] = model_name_override
    else:
        df["model_name"] = df.apply(
            lambda x: x["model"] + x["method"] if x["model"] == "boltz2" else x["model"],
            axis=1,
        )
        df["model_name"] = df["model_name"].map(lambda x: MODEL_NAME_ALIAS.get(x, x))
    n_before = len(df)
    if not keep_failed:
        df = df.dropna(subset=["min_rmsd_to_A", "min_rmsd_to_B"])
    n_dropped = n_before - len(df)
    if n_dropped > 0:
        print(
            f"\u26a0\ufe0f  load_min_altloc_rmsd({csv_path}): dropped "
            f"{n_dropped} rows with NaN min_rmsd_to_A/B"
        )
    df = df[
        [
            "protein",
            "model_name",
            "ensemble_size",
            "guidance_weight",
            "selection",
            "min_rmsd_to_A",
            "min_rmsd_to_B",
        ]
    ]
    df = add_simplified_metadata_columns(df)
    return df


def get_rscc_nan_proteins(
    models: tuple[str, ...] = ("boltz2X-RAY", "boltz2MD", "protenix", "rf3"),
) -> set[str]:
    """Return proteins with at least one NaN ``rscc`` value for the given models.

    Used by the notebooks' "Exclude NaN-RSCC proteins" toggle to produce a
    consistent cross-notebook blacklist. Derives from
    :func:`load_rscc_max_rmsd` (the primary plotting frame) so the blacklist
    matches the data actually shown.
    """
    df = load_rscc_max_rmsd()
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


def _parse_min_predicted_occupancy(s) -> float:
    """Return the minimum entry of an ``occupancies`` cell.

    Handles the numpy array repr used in ``lddt_results.csv``, which may have
    2 or more entries (e.g. ``"[0.375 0.625]"``, ``"[1. 0.]"``,
    ``"[0.5 0.25 0.25]"``). NaN for unparseable / missing cells.
    """
    if not isinstance(s, str):
        return np.nan
    parts = s.strip("[] ").split()
    if not parts:
        return np.nan
    try:
        return float(min(float(p) for p in parts))
    except ValueError:
        return np.nan


_ALTLOC_OCC_PATTERN = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")


def _parse_min_map_occupancy(s) -> float:
    """Return the minimum occupancy from a stringified ``altloc_occupancies`` dict.

    The ``altloc_occupancies`` column is the Python ``repr`` of a dict, e.g.
    ``"{'A': 0.5, 'B': 0.5}"`` or ``"{'A': 0.25, 'B': 0.75}"``. We pull out
    every numeric literal and return the minimum. NaN for unparseable cells.
    """
    if not isinstance(s, str):
        return np.nan
    nums = [float(m) for m in _ALTLOC_OCC_PATTERN.findall(s)]
    if not nums:
        return np.nan
    return float(min(nums))


def load_lddt(
    csv_path: str,
    model_name_override: str | None = None,
) -> pd.DataFrame:
    """Load LDDT results with graceful missing-file fallback.

    Parameters
    ----------
    csv_path : str
        Path to the LDDT CSV file.
    model_name_override : str | None
        If given, force-set ``model_name`` to this value. Used for the
        per-model LDDT concat to disambiguate ``protenix-mini`` /
        ``protenix-tiny`` (raw ``model == "protenix"``) and ``rf3-chiral-off``
        (raw ``model == "rf3"``), mirroring :func:`load_rscc`. When ``None``,
        model_name is built from the ``model``/``method`` columns and aliased.

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
        "min_predicted_occupancy",
        "min_map_occupancy",
        "clustering_score",
    ]
    if not os.path.exists(csv_path):
        print(f"⚠️  {csv_path} not found — LDDT sections will be empty")
        return pd.DataFrame(columns=_cols)
    df = pd.read_csv(csv_path, dtype={"error": str}).fillna({"error": ""})
    _occ = OCC_FILTER_ALTLOC
    df = df.query("altloc_occupancies == @_occ").copy()
    if model_name_override is not None:
        df["model_name"] = model_name_override
    else:
        df["model_name"] = df.apply(
            lambda x: x["model"] + x["method"] if x["model"] == "boltz2" else x["model"],
            axis=1,
        )
        df["model_name"] = df["model_name"].map(lambda x: MODEL_NAME_ALIAS.get(x, x))
    df["min_predicted_occupancy"] = df["occupancies"].map(_parse_min_predicted_occupancy)
    df["min_map_occupancy"] = df["altloc_occupancies"].map(_parse_min_map_occupancy)
    df["clustering_score"] = (
        df["min_predicted_occupancy"]
        * df["avg_silhouette_to_ref"]
        / df["min_map_occupancy"]
    )
    return df[_cols]


def load_classification(csv_path: str) -> pd.DataFrame:
    """Load conformational-change classification with lowercased protein names.

    Parameters
    ----------
    csv_path : str
        Path to a classification CSV (the ``CLASSIFICATION_CSV`` constant,
        i.e. ``results/classify_altloc_regions.csv``).

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
        sample = left.loc[unmatched, ["protein", "selection"]].drop_duplicates().head(3)
        print(
            f"⚠️  {n} RSCC rows had no containing classification; dropped. "
            f"Examples: {sample.to_dict(orient='records')}"
        )

    return merged.drop(columns=["_cls_start", "_cls_end", "_cls_selection", "_cls_span"])


def load_deposition_cutoffs(
    path: str = "results/Sampleworks_initial_structures.csv",
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
    path: str = "results/lddt_assignment_counts.csv",
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
