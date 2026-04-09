"""One-time script to integrate new RSCC and LDDT CSV files into initial40results.

New files:
- rscc_results_rf3_fixed_chiral_off.csv → replaces rf3-chiral-off rows in rscc_results.csv
- lddt_results_protenix_mini.csv → appends to lddt_results_w_baseline.csv
- lddt_results_protenix_tiny.csv → appends to lddt_results_w_baseline.csv
- lddt_results_rf3_fixed_chiral_off.csv → appends to lddt_results_w_baseline.csv

Usage:
    uv run python merge_new_rscc_lddt.py
"""

from __future__ import annotations

import polars as pl
from pathlib import Path

INITIAL_DIR = Path("initial40results")

# (source_file, canonical_model_name)
LDDT_FILES = [
    ("lddt_results_protenix_mini.csv", "protenix-mini"),
    ("lddt_results_protenix_tiny.csv", "protenix-tiny"),
    ("lddt_results_rf3_fixed_chiral_off.csv", "rf3-chiral-off"),
]

RSCC_FILE = ("rscc_results_rf3_fixed_chiral_off.csv", "rf3-chiral-off")


def merge_rscc() -> None:
    """Replace existing rf3-chiral-off RSCC rows with the new superset file."""
    src_file, canonical = RSCC_FILE
    target = INITIAL_DIR / "rscc_results.csv"
    print(f"[RSCC] Merging {src_file} → {target}")

    # Read existing, remove old rf3-chiral-off rows
    existing = pl.read_csv(target)
    before = len(existing)
    existing = existing.filter(pl.col("model") != "rf3-chiral-off")
    removed = before - len(existing)
    print(f"  Removed {removed} existing rf3-chiral-off rows")

    # Read new file, rewrite model, rename columns
    new = pl.read_csv(src_file)
    new = new.with_columns(pl.lit(canonical).alias("model"))
    if "trial_dir" in new.columns:
        new = new.rename({"trial_dir": "exp_dir"})

    # Align columns to existing schema
    for col in existing.columns:
        if col not in new.columns:
            new = new.with_columns(pl.lit(None).alias(col))
    new = new.select(existing.columns)

    # Match dtypes
    for col in existing.columns:
        if existing[col].dtype != new[col].dtype:
            new = new.with_columns(pl.col(col).cast(existing[col].dtype))

    result = pl.concat([existing, new])
    print(f"  {before} - {removed} + {len(new)} = {len(result)} rows")
    result.write_csv(target)


def merge_lddt() -> None:
    """Append new LDDT data for protenix-mini, protenix-tiny, rf3-chiral-off."""
    target = INITIAL_DIR / "lddt_results_w_baseline.csv"
    existing = pl.read_csv(target)
    print(f"\n[LDDT] Target: {target} ({len(existing)} existing rows)")

    frames = [existing]
    for src_file, canonical in LDDT_FILES:
        print(f"  Processing {src_file} → {canonical}")
        new = pl.read_csv(src_file)

        # Rewrite model column to canonical name
        new = new.with_columns(pl.lit(canonical).alias("model"))

        # Rename columns to match existing schema
        renames = {}
        if "altloc_occupancies" in new.columns:
            renames["altloc_occupancies"] = "occ_a"
        if "trial_dir" in new.columns:
            renames["trial_dir"] = "exp_dir"
        if renames:
            new = new.rename(renames)

        # Align columns
        for col in existing.columns:
            if col not in new.columns:
                new = new.with_columns(pl.lit(None).alias(col))
        new = new.select(existing.columns)

        print(f"    {len(new)} rows")
        frames.append(new)

    # Unify dtypes across all frames before concat: for each column, if any
    # frame has a different dtype, cast all to the most permissive type.
    target_dtypes = {}
    for col in existing.columns:
        dtypes = {f[col].dtype for f in frames if col in f.columns}
        if len(dtypes) > 1:
            # If any frame has String, cast all to String; otherwise use Float64
            target_dtypes[col] = pl.String if pl.String in dtypes else pl.Float64
    for i, f in enumerate(frames):
        casts = {col: f.schema[col] != dt for col, dt in target_dtypes.items() if col in f.columns}
        if any(casts.values()):
            frames[i] = f.with_columns([
                pl.col(col).cast(dt) for col, dt in target_dtypes.items()
                if col in f.columns and f.schema[col] != dt
            ])

    result = pl.concat(frames)
    print(f"  Total: {len(result)} rows")
    result.write_csv(target)


def main() -> None:
    merge_rscc()
    merge_lddt()
    print("\nDone!")


if __name__ == "__main__":
    main()
