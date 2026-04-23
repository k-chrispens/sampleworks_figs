"""One-time script to rebuild initial40results/rscc_results.csv from new CSVs.

RSCC rebuild (`merge_rscc`): full replacement of
``initial40results/rscc_results.csv`` from three source files in the project
root. The old RSCC data has been backed up as ``*.csv.bak`` and is abandoned.

Sources:
- rscc_results_other_models.csv               → boltz2 (MD+X-RAY), protenix
- rscc_results_rf3_fixed.csv                  → rf3 (fixed run; model kept as "rf3")
- rscc_results_rf3_fixed_chiral_off.csv       → rf3 (chiral-off run; model
                                                 rewritten to "rf3-chiral-off")

LDDT appends (`merge_lddt`): previously used to splice protenix-mini/tiny and
rf3-chiral-off rows into ``lddt_results_w_baseline.csv``. That merge has
already been run, so ``merge_lddt`` is not invoked by default — kept here for
reference / re-run if the target file is ever regenerated.

Usage:
    uv run python merge_new_rscc_lddt.py
"""

from __future__ import annotations

import polars as pl
from pathlib import Path

INITIAL_DIR = Path("initial40results")

# (source_file, canonical_model_name | None)
# None means "leave the model column as-is".
RSCC_SOURCES: list[tuple[str, str | None]] = [
    ("rscc_results_other_models.csv", None),
    ("rscc_results_rf3_fixed.csv", None),
    ("rscc_results_rf3_fixed_chiral_off.csv", "rf3-chiral-off"),
]

LDDT_FILES = [
    ("lddt_results_protenix_mini.csv", "protenix-mini"),
    ("lddt_results_protenix_tiny.csv", "protenix-tiny"),
    ("lddt_results_rf3_fixed_chiral_off.csv", "rf3-chiral-off"),
]


def merge_rscc() -> None:
    """Rebuild initial40results/rscc_results.csv from the three new CSVs."""
    target = INITIAL_DIR / "rscc_results.csv"
    print(f"[RSCC] Rebuilding {target} from {len(RSCC_SOURCES)} source files")

    frames: list[pl.DataFrame] = []
    for src_file, canonical in RSCC_SOURCES:
        df = pl.read_csv(src_file)
        if canonical is not None:
            df = df.with_columns(pl.lit(canonical).alias("model"))
        if "trial_dir" in df.columns:
            df = df.rename({"trial_dir": "exp_dir"})
        print(f"  {src_file}: {len(df)} rows (model rewrite: {canonical or '—'})")
        frames.append(df)

    result = pl.concat(frames, how="vertical_relaxed")
    print(f"  Total: {len(result)} rows → {target}")
    result.write_csv(target)


def merge_lddt() -> None:
    """Append new LDDT data for protenix-mini, protenix-tiny, rf3-chiral-off.

    Not called by default — already run. Retained for reference.
    """
    target = INITIAL_DIR / "lddt_results_w_baseline.csv"
    existing = pl.read_csv(target)
    print(f"\n[LDDT] Target: {target} ({len(existing)} existing rows)")

    frames = [existing]
    for src_file, canonical in LDDT_FILES:
        print(f"  Processing {src_file} → {canonical}")
        new = pl.read_csv(src_file)

        new = new.with_columns(pl.lit(canonical).alias("model"))

        renames = {}
        if "altloc_occupancies" in new.columns:
            renames["altloc_occupancies"] = "occ_a"
        if "trial_dir" in new.columns:
            renames["trial_dir"] = "exp_dir"
        if renames:
            new = new.rename(renames)

        for col in existing.columns:
            if col not in new.columns:
                new = new.with_columns(pl.lit(None).alias(col))
        new = new.select(existing.columns)

        print(f"    {len(new)} rows")
        frames.append(new)

    target_dtypes = {}
    for col in existing.columns:
        dtypes = {f[col].dtype for f in frames if col in f.columns}
        if len(dtypes) > 1:
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
    print("\nDone!")


if __name__ == "__main__":
    main()
