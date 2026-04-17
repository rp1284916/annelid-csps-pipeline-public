#!/usr/bin/env python3
"""Overlay curated display labels and clade assignments onto fresh tip metadata.

This helper is intended for the recovery path where a curated tip metadata file
has replaced the original exact ``leaf_id`` values with prettified labels.
Rather than reverse-engineering the exact ids, rerun ``s10`` to regenerate a
fresh canonical tip table, then copy the curated columns back by row order.
"""

from __future__ import annotations

import argparse
import csv
import shutil
from pathlib import Path


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_rows(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def overlay_by_row(
    fresh_tips: Path,
    curated_tips: Path,
    output: Path,
    backup_suffix: str,
) -> None:
    fresh_rows = read_rows(fresh_tips)
    curated_rows = read_rows(curated_tips)

    if len(fresh_rows) != len(curated_rows):
        raise ValueError(
            f"Row count mismatch: fresh={len(fresh_rows)} curated={len(curated_rows)}. "
            "This recovery path assumes the curated file came from the same exported table."
        )

    fieldnames = list(fresh_rows[0].keys())
    if "leaf_display_label" not in fieldnames:
        insert_at = fieldnames.index("leaf_id") + 1 if "leaf_id" in fieldnames else len(fieldnames)
        fieldnames.insert(insert_at, "leaf_display_label")
    if (
        any("clade_assignment" in row for row in curated_rows)
        or any("clade_assignment" in row for row in fresh_rows)
    ) and "clade_assignment" not in fieldnames:
        fieldnames.append("clade_assignment")

    if output == fresh_tips:
        backup_path = fresh_tips.with_name(fresh_tips.name + backup_suffix)
        shutil.copy2(fresh_tips, backup_path)
        print(f"Backed up fresh tip metadata to: {backup_path}")

    merged_rows: list[dict[str, str]] = []
    for index, (fresh, curated) in enumerate(zip(fresh_rows, curated_rows), start=1):
        if fresh.get("species_abbrev") != curated.get("species_abbrev"):
            raise ValueError(
                f"Row {index} species mismatch: fresh={fresh.get('species_abbrev')} "
                f"curated={curated.get('species_abbrev')}"
            )

        merged = dict(fresh)
        curated_display = (
            curated.get("leaf_display_label")
            or curated.get("leaf_id")
            or fresh.get("leaf_display_label")
            or fresh.get("leaf_id")
            or ""
        )
        merged["leaf_display_label"] = curated_display
        if "clade_assignment" in fresh or "clade_assignment" in curated:
            merged["clade_assignment"] = curated.get("clade_assignment", fresh.get("clade_assignment", ""))
        merged_rows.append(merged)

    write_rows(output, fieldnames, merged_rows)
    print(f"Wrote merged tip metadata to: {output}")
    print(f"Rows merged: {len(merged_rows)}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fresh-tips", required=True, help="Fresh exact tip metadata TSV from rerun s10.")
    parser.add_argument("--curated-tips", required=True, help="Curated tip metadata TSV with clade assignments.")
    parser.add_argument(
        "--output",
        help="Output TSV path. Defaults to overwriting --fresh-tips after making a backup.",
    )
    parser.add_argument(
        "--backup-suffix",
        default=".pre_overlay_backup.tsv",
        help="Suffix appended when backing up --fresh-tips before overwrite.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    fresh_tips = Path(args.fresh_tips)
    output = Path(args.output) if args.output else fresh_tips
    overlay_by_row(
        fresh_tips=fresh_tips,
        curated_tips=Path(args.curated_tips),
        output=output,
        backup_suffix=args.backup_suffix,
    )


if __name__ == "__main__":
    main()
