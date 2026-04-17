#!/usr/bin/env python3
"""Convert curated clade assignments into node-anchored tree annotations.

This helper bridges the two tree-annotation modes used in the annelid pipeline:

1. a clade-level curation table keyed by biological module names
2. a node-level annotation table keyed by internal UPGMA node ids

For each visible clade in the curated template, the script finds the common
ancestor of all leaves assigned to that clade and writes out a new TSV with an
explicit `anchor_node_id`. This makes it straightforward to render a
manuscript-style tree whose right-hand labels correspond to selected internal
nodes rather than generic visual clusters.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

from Bio import Phylo


def normalize_clade_key(value: str) -> str:
    value = re.sub(r"\s*/\s*", "/", value.strip())
    value = re.sub(r"\s+", " ", value)
    return value.casefold()


def parse_bool(value: str | None, default: bool = False) -> bool:
    if value is None or value == "":
        return default
    return value.strip().lower() in {"1", "true", "yes", "y"}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def ensure_internal_ids(tree) -> None:
    existing = {
        clade.name
        for clade in tree.find_clades(order="preorder")
        if not clade.is_terminal() and clade.name not in {None, ""}
    }
    counter = 1
    for clade in tree.find_clades(order="preorder"):
        if clade.is_terminal():
            continue
        if clade.name in {None, ""}:
            if clade is tree.root:
                clade.name = "Root"
            else:
                while f"Node{counter}" in existing:
                    counter += 1
                clade.name = f"Node{counter}"
                existing.add(clade.name)
                counter += 1


def build_node_template(tree_path: Path, tips_path: Path, clades_path: Path, output_path: Path) -> None:
    tree = Phylo.read(str(tree_path), "newick")
    ensure_internal_ids(tree)
    leaves_by_name = {leaf.name: leaf for leaf in tree.get_terminals()}

    tip_rows = read_tsv(tips_path)
    clade_rows = read_tsv(clades_path)

    leaf_ids_by_clade: dict[str, list[str]] = {}
    for row in tip_rows:
        clade_assignment = (row.get("clade_assignment") or "").strip()
        if not clade_assignment:
            continue
        leaf_ids_by_clade.setdefault(normalize_clade_key(clade_assignment), []).append(row["leaf_id"])

    out_rows: list[dict[str, str]] = []
    for row in clade_rows:
        clade_key = (row.get("clade_key") or row.get("clade_name") or "").strip()
        if not clade_key or not parse_bool(row.get("show"), default=True):
            continue

        normalized_key = normalize_clade_key(clade_key)
        assigned_leaf_ids = [leaf_id for leaf_id in leaf_ids_by_clade.get(normalized_key, []) if leaf_id in leaves_by_name]
        if not assigned_leaf_ids:
            continue

        ancestor = tree.common_ancestor([leaves_by_name[leaf_id] for leaf_id in assigned_leaf_ids])
        out_rows.append(
            {
                "clade_key": clade_key,
                "clade_name": row.get("clade_name") or clade_key,
                "label_color": row.get("label_color") or row.get("color") or "#2f2f2f",
                # Leave tf_text blank so render_tf_program_tree.py can derive
                # shared descendant TFs from UPGMA.markers.tsv.
                "tf_text": "",
                "show": "TRUE",
                "anchor_node_id": ancestor.name or "",
            }
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["clade_key", "clade_name", "label_color", "tf_text", "show", "anchor_node_id"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(out_rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tree", type=Path, required=True)
    parser.add_argument("--tips", type=Path, required=True)
    parser.add_argument("--clades", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    build_node_template(args.tree, args.tips, args.clades, args.output)
    print(f"Wrote node-anchored annotation template to {args.output}")


if __name__ == "__main__":
    main()
