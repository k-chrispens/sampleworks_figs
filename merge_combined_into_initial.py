"""One-time script to merge combined_occ_sweep_results/ into initial40results/.

Reads each CSV from combined_occ_sweep_results/, rewrites model identifiers
to unique canonical names, and appends rows to the corresponding initial40results/ CSV.

Usage:
    pixi run python merge_combined_into_initial.py
    # or
    uv run python merge_combined_into_initial.py
"""

from __future__ import annotations

import shutil
from pathlib import Path

import polars as pl

INITIAL_DIR = Path("initial40results")
COMBINED_DIR = Path("combined_occ_sweep_results")
BACKUP_DIR = Path("initial40results_backup")

# source_model → canonical model name
SOURCE_TO_CANONICAL = {
    "protenix_mini": "protenix-mini",
    "protenix_tiny": "protenix-tiny",
    "rf3": "rf3-fixed",
    "rf3_chiral_off": "rf3-chiral-off",
}

# source_model → (old path segment after occ_spec, new path segment)
PATH_MODEL_SEGMENT = {
    "protenix_mini": ("protenix_X-RAY_DIFFRACTION", "protenix-mini"),
    "protenix_tiny": ("protenix_X-RAY_DIFFRACTION", "protenix-tiny"),
    "rf3": ("rf3", "rf3-fixed"),
    "rf3_chiral_off": ("rf3_X-RAY_DIFFRACTION", "rf3-chiral-off"),
}

# Anchored regex to replace only the model directory segment (after occ spec),
# not parent directories like occ_sweep_results_rf3_fixed.
# Captures: .../{occ_spec}/{MODEL_NAME}/{steering}/...
import re
OCC_SPEC_PATTERN = r"([\d.]+occ[A-Z](?:_[\d.]+occ[A-Z])*)/{}/"


def _rewrite_paths(df: pl.LazyFrame, path_col: str) -> pl.LazyFrame:
    """Rewrite model directory in paths based on source_model column."""
    frames = []
    for source_model, (old_seg, new_seg) in PATH_MODEL_SEGMENT.items():
        pattern = OCC_SPEC_PATTERN.format(re.escape(old_seg))
        replacement = f"${{1}}/{new_seg}/"
        chunk = (
            df.filter(pl.col("source_model") == source_model)
            .with_columns(
                pl.col(path_col).str.replace(pattern, replacement)
            )
        )
        frames.append(chunk)
    return pl.concat(frames)


def merge_path_based(
    combined_name: str,
    initial_name: str,
    path_col: str,
    initial_col_order: list[str],
    strip_cols: list[str] | None = None,
) -> None:
    """Merge a path-based CSV from combined into initial."""
    combined_path = COMBINED_DIR / combined_name
    initial_path = INITIAL_DIR / initial_name
    print(f"  Merging {combined_path} → {initial_path}")

    combined = pl.scan_csv(combined_path)
    combined = _rewrite_paths(combined, path_col)
    combined = combined.drop("source_model").select(initial_col_order)

    if strip_cols:
        for col in strip_cols:
            combined = combined.with_columns(
                pl.col(col).cast(pl.String).str.strip_chars().alias(col)
            )

    new_data = combined.collect()
    print(f"    New rows: {len(new_data)}")

    initial = pl.read_csv(initial_path)
    print(f"    Existing rows: {len(initial)}")

    # Ensure matching dtypes before concat
    for col in initial_col_order:
        if col in initial.columns and col in new_data.columns:
            if initial[col].dtype != new_data[col].dtype:
                target_dtype = initial[col].dtype
                new_data = new_data.with_columns(pl.col(col).cast(target_dtype))

    result = pl.concat([initial, new_data])
    print(f"    Merged rows: {len(result)}")
    result.write_csv(initial_path)


def merge_rscc() -> None:
    """Merge RSCC results (column-based, not path-based)."""
    combined_path = COMBINED_DIR / "rscc_results.csv"
    initial_path = INITIAL_DIR / "rscc_results.csv"
    print(f"  Merging {combined_path} → {initial_path}")

    initial = pl.read_csv(initial_path)
    initial_cols = initial.columns
    print(f"    Existing rows: {len(initial)}")

    combined = pl.read_csv(combined_path)

    # Rewrite model column using source_model mapping
    combined = combined.with_columns(
        pl.col("source_model").replace_strict(SOURCE_TO_CANONICAL).alias("model")
    )

    # Rename trial_dir → exp_dir to match initial schema
    if "trial_dir" in combined.columns:
        combined = combined.rename({"trial_dir": "exp_dir"})

    combined = combined.drop("source_model")

    # Select only columns that exist in initial, in the same order
    select_cols = [c for c in initial_cols if c in combined.columns]
    combined = combined.select(select_cols)

    # Add any missing columns as null
    for c in initial_cols:
        if c not in combined.columns:
            combined = combined.with_columns(pl.lit(None).alias(c))
    combined = combined.select(initial_cols)

    # Match dtypes
    for col in initial_cols:
        if initial[col].dtype != combined[col].dtype:
            combined = combined.with_columns(pl.col(col).cast(initial[col].dtype))

    print(f"    New rows: {len(combined)}")
    result = pl.concat([initial, combined])
    print(f"    Merged rows: {len(result)}")
    result.write_csv(initial_path)


def main() -> None:
    # Step 0: Backup
    if not BACKUP_DIR.exists():
        print(f"Backing up {INITIAL_DIR} → {BACKUP_DIR}")
        shutil.copytree(INITIAL_DIR, BACKUP_DIR)
    else:
        print(f"Backup already exists at {BACKUP_DIR}, skipping backup")

    # ── Path-based files ─────────────────────────────────────────────────────

    # 1. clashscore_metrics.csv
    print("\n[1/8] clashscore_metrics.csv")
    merge_path_based(
        "clashscore_metrics.csv",
        "clashscore_metrics.csv",
        path_col="model_name",
        initial_col_order=["model_name", "model_id", "clashscore", "num_clashes"],
        strip_cols=["model_id"],
    )

    # 2. bond_angle_violation_fractions.csv
    print("\n[2/8] bond_angle_violation_fractions.csv")
    merge_path_based(
        "bond_angle_violation_fractions.csv",
        "bond_angle_violation_fractions.csv",
        path_col="model",
        initial_col_order=["outlier_fraction", "model", "model_n"],
    )

    # 3. bond_length_violation_fractions.csv
    print("\n[3/8] bond_length_violation_fractions.csv")
    merge_path_based(
        "bond_length_violation_fractions.csv",
        "bond_length_violation_fractions.csv",
        path_col="model",
        initial_col_order=["outlier_fraction", "model", "model_n"],
    )

    # 4. tortoize_protein_stats.csv
    print("\n[4/8] tortoize_protein_stats.csv")
    merge_path_based(
        "tortoize_protein_stats.csv",
        "initial_40_tortoize_protein_stats.csv",
        path_col="path",
        initial_col_order=[
            "model", "ramachandran_z_score", "ramachandran_jackknife_sd",
            "torsion_z_score", "torsion_jackknife_sd", "path",
        ],
    )

    # 5. tortoize_residues.csv (large ~2.5GB)
    print("\n[5/8] tortoize_residues.csv (large file)")
    merge_path_based(
        "tortoize_residues.csv",
        "initial_40_tortoize_residues.csv",
        path_col="path",
        initial_col_order=[
            "model", "asymID", "compID", "seqID", "strandID",
            "insCode", "ss_type", "ramachandran_z_score", "torsion_z_score", "path",
        ],
    )

    # 6. bond_angle_outliers.csv (large ~2.5GB, column reorder needed)
    print("\n[6/8] bond_angle_outliers.csv (large file, column reorder)")
    merge_path_based(
        "bond_angle_outliers.csv",
        "bond_angle_outliers.csv",
        path_col="model",
        initial_col_order=[
            "atom1_residue_id", "atom1_residue_name", "atom1_chain_id", "atom1_atom_id",
            "atom2_residue_id", "atom2_residue_name", "atom2_chain_id", "atom2_atom_id",
            "center_atom_residue_id", "center_atom_residue_name",
            "center_atom_chain_id", "center_atom_atom_id",
            "bond_length", "lower_bound", "upper_bound",
            "model", "model_n",
        ],
    )

    # 7. bond_length_outliers.csv (large ~3.5GB)
    print("\n[7/8] bond_length_outliers.csv (large file)")
    merge_path_based(
        "bond_length_outliers.csv",
        "bond_length_outliers.csv",
        path_col="model",
        initial_col_order=[
            "atom1_residue_id", "atom1_residue_name", "atom1_chain_id", "atom1_atom_id",
            "atom2_residue_id", "atom2_residue_name", "atom2_chain_id", "atom2_atom_id",
            "bond_length", "lower_bound", "upper_bound",
            "model", "model_n",
        ],
    )

    # ── RSCC (column-based) ──────────────────────────────────────────────────

    # 8. rscc_results.csv
    print("\n[8/8] rscc_results.csv")
    merge_rscc()

    print("\nDone! All files merged successfully.")


if __name__ == "__main__":
    main()
