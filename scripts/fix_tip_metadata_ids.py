#!/usr/bin/env python3
"""Repair curated tip metadata so leaf_id matches exact tree tip names.

This helper is for the common curation workflow where the user edits the
human-readable labels in ``leaf_id`` and later needs to restore the original
Newick tip ids expected by the renderer. The script keeps the readable labels
in ``leaf_display_label`` and rewrites ``leaf_id`` to exact tree names.
"""

from __future__ import annotations

import argparse
import csv
import re
import shutil
from collections import defaultdict
from pathlib import Path

from Bio import Phylo


def normalize(text: str) -> str:
    """Normalize a label for fuzzy matching across punctuation variants."""
    text = text.strip().lower()
    text = text.replace("|", " ")
    text = re.sub(r"[./_-]+", " ", text)
    text = re.sub(r"[^a-z0-9 ]+", " ", text)
    text = re.sub(r"\s+", " ", text).strip()
    return text


def build_candidate_strings(row: dict[str, str]) -> list[str]:
    """Return raw strings that may correspond to the original tree tip."""
    candidates: list[str] = []
    for key in ("leaf_id", "leaf_display_label"):
        value = (row.get(key) or "").strip()
        if value and value not in candidates:
            candidates.append(value)

    species = (row.get("species_abbrev") or "").strip()
    cell_state = (row.get("cell_state_label") or "").strip()
    if species and cell_state:
        combined = f"{species}|{cell_state}"
        if combined not in candidates:
            candidates.append(combined)

    return candidates


def choose_unique_normalized_match(
    row: dict[str, str],
    tree_by_norm: dict[str, list[str]],
    used_tips: set[str],
) -> str | None:
    """Return a single unused tip when fuzzy matching is unambiguous."""
    matches: set[str] = set()
    for candidate in build_candidate_strings(row):
        for tip in tree_by_norm.get(normalize(candidate), []):
            if tip not in used_tips:
                matches.add(tip)
    if len(matches) == 1:
        return next(iter(matches))
    return None


def group_key_for_row(row: dict[str, str], tree_by_norm: dict[str, list[str]]) -> str:
    """Return the display-driven group key for ambiguous duplicate labels."""
    preferred = (row.get("leaf_display_label") or row.get("leaf_id") or "").strip()
    if normalize(preferred) in tree_by_norm:
        return normalize(preferred)

    match = re.match(r"^(.*)\.(\d+)$", preferred)
    if match:
        base = match.group(1)
        if normalize(base) in tree_by_norm:
            return normalize(base)

    return normalize(preferred)


def repair_tip_metadata(tree_path: Path, tips_path: Path, backup_suffix: str) -> None:
    """Rewrite ``leaf_id`` to exact tree tips while preserving display labels."""
    tree_tips = [tip.name for tip in Phylo.read(tree_path, "newick").get_terminals()]
    tree_tip_set = set(tree_tips)
    tree_by_norm: dict[str, list[str]] = defaultdict(list)
    for tip in tree_tips:
        tree_by_norm[normalize(tip)].append(tip)
    for key in tree_by_norm:
        tree_by_norm[key] = sorted(tree_by_norm[key])

    with tips_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError("Tip metadata file is empty.")

    fieldnames = list(rows[0].keys())
    if "leaf_display_label" not in fieldnames:
        insert_at = fieldnames.index("leaf_id") + 1 if "leaf_id" in fieldnames else len(fieldnames)
        fieldnames.insert(insert_at, "leaf_display_label")

    backup_path = tips_path.with_name(tips_path.name + backup_suffix)
    shutil.copy2(tips_path, backup_path)

    used_tips: set[str] = set()
    unresolved_rows: list[dict[str, str]] = []

    for row in rows:
        original_label = (row.get("leaf_id") or "").strip()
        if not row.get("leaf_display_label"):
            row["leaf_display_label"] = original_label

        if original_label in tree_tip_set and original_label not in used_tips:
            used_tips.add(original_label)
            continue

        unique_match = choose_unique_normalized_match(row, tree_by_norm, used_tips)
        if unique_match is not None:
            row["leaf_id"] = unique_match
            used_tips.add(unique_match)
            continue

        unresolved_rows.append(row)

    ambiguous_groups: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in unresolved_rows:
        ambiguous_groups[group_key_for_row(row, tree_by_norm)].append(row)

    for group_key, group_rows in ambiguous_groups.items():
        candidate_tips = [tip for tip in tree_by_norm.get(group_key, []) if tip not in used_tips]
        if len(candidate_tips) != len(group_rows):
            raise ValueError(
                f"Group size mismatch for '{group_key}': "
                f"{len(group_rows)} metadata rows vs {len(candidate_tips)} tree tips. "
                f"Metadata rows: {[row['leaf_id'] for row in group_rows]}; "
                f"Tree tips: {candidate_tips}"
            )

        for row, matched_tip in zip(sorted(group_rows, key=lambda item: item["leaf_id"]), candidate_tips):
            row["leaf_id"] = matched_tip
            used_tips.add(matched_tip)

    fixed_tips = {row["leaf_id"] for row in rows}
    missing = sorted(tree_tip_set - fixed_tips)
    extra = sorted(fixed_tips - tree_tip_set)
    if missing or extra:
        raise ValueError(
            f"Repair incomplete. Missing={len(missing)} Extra={len(extra)}. "
            f"Examples missing: {missing[:5]} extra: {extra[:5]}"
        )

    rows.sort(key=lambda row: row["leaf_id"])
    with tips_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Backed up original metadata to: {backup_path}")
    print(f"Repaired {len(rows)} rows.")
    print("Missing from metadata: 0")
    print("Extra in metadata: 0")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tree", required=True, help="Path to the rooted Newick tree.")
    parser.add_argument("--tips", required=True, help="Path to the curated tip metadata TSV.")
    parser.add_argument(
        "--backup-suffix",
        default=".pretty_backup.tsv",
        help="Suffix appended to the original TSV before rewriting.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    repair_tip_metadata(Path(args.tree), Path(args.tips), args.backup_suffix)


if __name__ == "__main__":
    main()
