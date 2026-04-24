"""One-time script to (re)build initial40results/{rscc,lddt}_results.csv.

Both rebuilds are full replacements; the previous CSVs have been backed up as
``*.csv.bak`` in ``initial40results/`` and are no longer used.

RSCC rebuild (`merge_rscc`) — 3 source files in the project root:
- rscc_results_other_models.csv               → boltz2 (MD+X-RAY), protenix
- rscc_results_rf3_fixed.csv                  → rf3 (fixed run; model kept as "rf3")
- rscc_results_rf3_fixed_chiral_off.csv       → rf3 (chiral-off run; model
                                                 rewritten to "rf3-chiral-off")

LDDT rebuild (`merge_lddt`) — 3 source files in the project root:
- lddt_results_other_models.csv               → boltz2 (MD+X-RAY), protenix
- lddt_results_rf3_fixed.csv                  → rf3 (fixed run; "rf3")
- lddt_results_rf3_fixed_chiral_off.csv       → "rf3-chiral-off"

Note: protenix-mini / protenix-tiny are absent from both CSVs by design —
their new RSCC and LDDT data have not been produced yet. The
``lddt_results_protenix_{mini,tiny}.csv`` files in the project root are an
older run and are intentionally not merged here.

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

LDDT_SOURCES: list[tuple[str, str | None]] = [
    ("lddt_results_other_models.csv", None),
    ("lddt_results_rf3_fixed.csv", None),
    ("lddt_results_rf3_fixed_chiral_off.csv", "rf3-chiral-off"),
]


def _rebuild(sources: list[tuple[str, str | None]], target: Path, label: str) -> None:
    print(f"[{label}] Rebuilding {target} from {len(sources)} source files")
    frames: list[pl.DataFrame] = []
    for src_file, canonical in sources:
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


def merge_rscc() -> None:
    """Rebuild initial40results/rscc_results.csv from the 3 new RSCC CSVs."""
    _rebuild(RSCC_SOURCES, INITIAL_DIR / "rscc_results.csv", "RSCC")


def merge_lddt() -> None:
    """Rebuild initial40results/lddt_results.csv from the 5 new LDDT CSVs."""
    _rebuild(LDDT_SOURCES, INITIAL_DIR / "lddt_results.csv", "LDDT")


def main() -> None:
    merge_rscc()
    merge_lddt()
    print("\nDone!")


if __name__ == "__main__":
    main()
