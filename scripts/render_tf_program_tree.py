#!/usr/bin/env python3
"""Render a manuscript-style TF-annotated cell-type tree.

The design is intentionally asymmetric:

1. Left: rooted dendrogram with thin branch lines and support labels.
2. Middle: full leaf labels, with only small species cues as secondary metadata.
3. Right: clade/module labels plus the defining TF program text.

The TF program annotation column is the main interpretive output. Tree topology
is only the scaffold that connects those TF programs to their corresponding
selected clades.
"""

from __future__ import annotations

import argparse
import csv
import io
import re
import textwrap
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence

try:
    from Bio import Phylo
except ImportError as exc:  # pragma: no cover - import guard
    raise SystemExit(
        "Biopython is required for this script. Install it with 'pip install biopython'."
    ) from exc

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError as exc:  # pragma: no cover - import guard
    raise SystemExit(
        "matplotlib is required for this script. Install it with 'pip install matplotlib'."
    ) from exc


EXAMPLE_NEWICK = """(
    ((Ctel|Gut.1,Ofus|Ventral.Pharynx.2):0.8,(Ctel|Gut.3,Ofus|Digestive.4):0.9)Node1:0.7,
    (((Ctel|Neurosecretory.12,Ofus|Neurosecretory.7):0.8,Ctel|Neurosecretory.5):0.7,
      (Ofus|Trunk.6,Ctel|Trunk.1):0.6)Node2:0.5,
    (Ctel|Ciliated.Pharynx,Ofus|Ciliated.Band):0.9
)Root;
"""

EXAMPLE_LEAF_METADATA = {
    "Ctel|Gut.1": {
        "species": "Ctel",
        "cell_state": "Gut.1",
        "label": "Ctel|Gut.1",
        "tip_color": "#d95f02",
        "species_color": "#4c78a8",
    },
    "Ofus|Ventral.Pharynx.2": {
        "species": "Ofus",
        "cell_state": "Ventral.Pharynx.2",
        "label": "Ofus|Ventral.Pharynx.2",
        "tip_color": "#d95f02",
        "species_color": "#e45756",
    },
    "Ctel|Gut.3": {
        "species": "Ctel",
        "cell_state": "Gut.3",
        "label": "Ctel|Gut.3",
        "tip_color": "#d95f02",
        "species_color": "#4c78a8",
    },
    "Ofus|Digestive.4": {
        "species": "Ofus",
        "cell_state": "Digestive.4",
        "label": "Ofus|Digestive.4",
        "tip_color": "#d95f02",
        "species_color": "#e45756",
    },
    "Ctel|Neurosecretory.12": {
        "species": "Ctel",
        "cell_state": "Neurosecretory.12",
        "label": "Ctel|Neurosecretory.12",
        "tip_color": "#7570b3",
        "species_color": "#4c78a8",
    },
    "Ofus|Neurosecretory.7": {
        "species": "Ofus",
        "cell_state": "Neurosecretory.7",
        "label": "Ofus|Neurosecretory.7",
        "tip_color": "#7570b3",
        "species_color": "#e45756",
    },
    "Ctel|Neurosecretory.5": {
        "species": "Ctel",
        "cell_state": "Neurosecretory.5",
        "label": "Ctel|Neurosecretory.5",
        "tip_color": "#7570b3",
        "species_color": "#4c78a8",
    },
    "Ofus|Trunk.6": {
        "species": "Ofus",
        "cell_state": "Trunk.6",
        "label": "Ofus|Trunk.6",
        "tip_color": "#1b9e77",
        "species_color": "#e45756",
    },
    "Ctel|Trunk.1": {
        "species": "Ctel",
        "cell_state": "Trunk.1",
        "label": "Ctel|Trunk.1",
        "tip_color": "#1b9e77",
        "species_color": "#4c78a8",
    },
    "Ctel|Ciliated.Pharynx": {
        "species": "Ctel",
        "cell_state": "Ciliated.Pharynx",
        "label": "Ctel|Ciliated.Pharynx",
        "tip_color": "#66a61e",
        "species_color": "#4c78a8",
    },
    "Ofus|Ciliated.Band": {
        "species": "Ofus",
        "cell_state": "Ciliated.Band",
        "label": "Ofus|Ciliated.Band",
        "tip_color": "#66a61e",
        "species_color": "#e45756",
    },
}

EXAMPLE_INTERNAL_NODE_ANNOTATIONS = {
    "Node2": {
        "label": "Example module A",
        "color": "#009E73",
        "tf_list": ["TF_A", "TF_B", "TF_C"],
    },
    "Node1": {
        "label": "Example module B",
        "color": "#E69F00",
        "tf_list": ["marker set B"],
    },
    "Node3": {
        "label": "Example module C",
        "color": "#CC79A7",
        "tf_list": ["marker set C"],
    },
    "Node4": {
        "label": "Example module D",
        "color": "#56B4E9",
        "tf_list": ["TF_D", "marker set D"],
    },
}


@dataclass(frozen=True)
class TipMetadata:
    leaf_id: str
    species_abbrev: str
    cell_state_label: str
    leaf_display_label: str
    tip_color: str
    species_color: str
    clade_assignment: str


@dataclass(frozen=True)
class CladeAnnotation:
    clade_key: str
    display_name: str
    label_color: str
    tf_text: str
    show: bool
    anchor_node_id: str | None
    support: float | None


@dataclass(frozen=True)
class CladePlacement:
    clade_key: str
    anchor_clade: object
    target_y: float
    support: float | None


@dataclass(frozen=True)
class MarkerRecord:
    orthogroup: str
    orthogroup_name: str
    probability: float
    present: bool
    node_id: str


def normalize_clade_key(value: str) -> str:
    """Normalize clade names so minor spacing and slash differences still match."""
    value = re.sub(r"\s*/\s*", "/", value.strip())
    value = re.sub(r"\s+", " ", value)
    return value.casefold()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Render a TF-annotated rooted tree with manuscript-style layout."
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        help=(
            "Pipeline results directory containing csps.<comparison>.dendrogram.<focus>.* files. "
            "When provided, tree/tips/clades/node-support/markers are inferred automatically unless overridden."
        ),
    )
    parser.add_argument(
        "--comparison",
        default="annelids",
        help="Comparison id used in the pipeline filenames when --results-dir is used.",
    )
    parser.add_argument(
        "--focus-id",
        default="cts",
        help="Focus id used in the pipeline filenames when --results-dir is used.",
    )
    parser.add_argument("--tree", type=Path, help="Rooted Newick tree.")
    parser.add_argument("--tips", type=Path, help="Tip metadata TSV.")
    parser.add_argument("--clades", type=Path, help="Clade annotation TSV.")
    parser.add_argument(
        "--node-support",
        type=Path,
        help="Optional node support TSV exported by s10.",
    )
    parser.add_argument(
        "--output-prefix",
        type=Path,
        default=Path("tf_program_tree"),
        help="Output prefix without extension.",
    )
    parser.add_argument(
        "--example",
        action="store_true",
        help="Render the built-in synthetic example.",
    )
    parser.add_argument(
        "--write-example",
        type=Path,
        help="Write the synthetic example tree and TSV files to this directory and exit.",
    )
    parser.add_argument("--png-dpi", type=int, default=600, help="PNG resolution.")
    parser.add_argument("--figure-width", type=float, default=14.5, help="Figure width in inches.")
    parser.add_argument("--leaf-height", type=float, default=0.44, help="Height contribution per leaf.")
    parser.add_argument("--min-height", type=float, default=8.0, help="Minimum figure height.")
    parser.add_argument("--support-min", type=float, default=50.0, help="Minimum support to display.")
    parser.add_argument("--tf-wrap", type=int, default=48, help="Wrap width for TF program text.")
    parser.add_argument("--annotation-gap", type=float, default=0.7, help="Minimum gap between annotation blocks.")
    parser.add_argument(
        "--markers",
        type=Path,
        help="Optional UPGMA.markers.tsv file used to derive shared descendant TFs for selected nodes.",
    )
    parser.add_argument(
        "--node-probability-min",
        type=float,
        default=0.7,
        help="Minimum node posterior probability for a TF to be considered shared in a selected node.",
    )
    parser.add_argument(
        "--outside-probability-max",
        type=float,
        default=0.45,
        help="Maximum tolerated posterior probability in disjoint outside nodes when deriving node TFs.",
    )
    parser.add_argument(
        "--max-node-tfs",
        type=int,
        default=8,
        help="Maximum number of derived TFs shown per selected internal node.",
    )
    return parser.parse_args()


def resolve_pipeline_paths(args: argparse.Namespace) -> argparse.Namespace:
    """Infer the standard pipeline filenames from a results directory."""
    if not args.results_dir:
        return args

    base = f"csps.{args.comparison}.dendrogram.{args.focus_id}"
    resolved = argparse.Namespace(**vars(args))

    if resolved.tree is None:
        resolved.tree = args.results_dir / f"{base}.UPGMA.newick"
    if resolved.tips is None:
        resolved.tips = args.results_dir / f"{base}.tip_metadata.tsv"
    if resolved.clades is None:
        resolved.clades = args.results_dir / f"{base}.tf_program_template.tsv"
    if resolved.node_support is None:
        candidate = args.results_dir / f"{base}.node_support.tsv"
        resolved.node_support = candidate if candidate.exists() else None
    if resolved.markers is None:
        candidate = args.results_dir / f"{base}.UPGMA.markers.tsv"
        resolved.markers = candidate if candidate.exists() else None
    if resolved.output_prefix == Path("tf_program_tree"):
        resolved.output_prefix = args.results_dir / f"{base}.tf_program_tree"

    return resolved


def parse_bool(value: str | None, default: bool = False) -> bool:
    if value is None or value == "":
        return default
    return value.strip().lower() in {"1", "true", "yes", "y"}


def parse_float(value: str | None) -> float | None:
    if value is None or value == "":
        return None
    return float(value)


def read_tsv_rows(path: Path) -> List[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def add_extension(prefix: Path, extension: str) -> Path:
    """Append an output extension without treating dotted prefixes as suffixes."""
    return Path(f"{prefix}.{extension}")


def load_tip_metadata(rows: Sequence[dict[str, str]]) -> Dict[str, TipMetadata]:
    """Load leaf-level metadata and keep the display label separate from species cues."""
    tip_metadata: Dict[str, TipMetadata] = {}
    for row in rows:
        leaf_id = row["leaf_id"]
        species_abbrev = row.get("species_abbrev") or leaf_id.split("|", 1)[0]
        cell_state_label = row.get("cell_state_label") or row.get("cell_type_label") or leaf_id.split("|", 1)[-1]
        leaf_display_label = row.get("leaf_display_label") or row.get("label") or leaf_id
        clade_assignment = (row.get("clade_assignment") or "").strip()
        tip_metadata[leaf_id] = TipMetadata(
            leaf_id=leaf_id,
            species_abbrev=species_abbrev,
            cell_state_label=cell_state_label,
            leaf_display_label=leaf_display_label,
            tip_color=row.get("tip_color") or row.get("color") or "#6c6c6c",
            species_color=row.get("species_color") or "#555555",
            clade_assignment=clade_assignment,
        )
    return tip_metadata


def load_clade_annotations(rows: Sequence[dict[str, str]]) -> Dict[str, CladeAnnotation]:
    """Load the right-hand annotation system keyed by node or biological module."""
    annotations: Dict[str, CladeAnnotation] = {}
    for row in rows:
        clade_key = row.get("clade_key") or row.get("node_id") or row.get("anchor_node_id") or row.get("clade_name")
        if not clade_key:
            continue
        clade_key = clade_key.strip()
        tf_text = row.get("tf_text") or row.get("tf_program") or row.get("tf_list") or row.get("suggested_tfs") or ""
        display_name = row.get("clade_name") or row.get("label") or clade_key
        annotations[clade_key] = CladeAnnotation(
            clade_key=clade_key,
            display_name=display_name,
            label_color=row.get("label_color") or row.get("color") or "#2f2f2f",
            tf_text=tf_text,
            show=parse_bool(row.get("show"), default=bool(display_name or tf_text)),
            anchor_node_id=row.get("anchor_node_id") or row.get("node_id"),
            support=parse_float(row.get("support")),
        )
    return annotations


def load_node_support(rows: Sequence[dict[str, str]]) -> Dict[str, float]:
    """Load support values keyed by internal node ID."""
    support_lookup: Dict[str, float] = {}
    for row in rows:
        node_id = row.get("node_id")
        support = parse_float(row.get("support"))
        if node_id and support is not None:
            support_lookup[node_id] = support
    return support_lookup


def load_marker_records(rows: Sequence[dict[str, str]]) -> List[MarkerRecord]:
    """Load per-node ancestral TF probabilities exported by s10."""
    records: List[MarkerRecord] = []
    for row in rows:
        node_id = (row.get("node") or row.get("node_id") or "").strip()
        probability = parse_float(row.get("probability"))
        if not node_id or probability is None:
            continue
        records.append(
            MarkerRecord(
                orthogroup=(row.get("orthogroup") or "").strip(),
                orthogroup_name=(row.get("orthogroup_name") or row.get("orthogroup") or "").strip(),
                probability=probability,
                present=parse_bool(row.get("present"), default=probability >= 0.7),
                node_id=node_id,
            )
        )
    return records


def example_tip_rows() -> List[dict[str, str]]:
    rows = []
    for leaf_id, meta in EXAMPLE_LEAF_METADATA.items():
        rows.append(
            {
                "leaf_id": leaf_id,
                "species_abbrev": meta["species"],
                "cell_state_label": meta["cell_state"],
                "leaf_display_label": meta["label"],
                "tip_color": meta["tip_color"],
                "species_color": meta["species_color"],
                "clade_assignment": "",
            }
        )
    return rows


def example_clade_rows() -> List[dict[str, str]]:
    rows = []
    for node_id, meta in EXAMPLE_INTERNAL_NODE_ANNOTATIONS.items():
        rows.append(
            {
                "clade_key": node_id,
                "clade_name": meta["label"],
                "label_color": meta["color"],
                "tf_text": ", ".join(meta["tf_list"]),
                "show": "TRUE",
                "anchor_node_id": node_id,
            }
        )
    return rows


def load_example_dataset():
    tree = Phylo.read(io.StringIO(EXAMPLE_NEWICK), "newick")
    ensure_internal_ids(tree)
    tips = load_tip_metadata(example_tip_rows())
    clades = load_clade_annotations(example_clade_rows())
    return tree, tips, clades


def write_example_dataset(outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "example_tree.newick").write_text(EXAMPLE_NEWICK, encoding="utf-8")
    with (outdir / "example_tips.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "leaf_id",
                "species_abbrev",
                "cell_state_label",
                "leaf_display_label",
                "tip_color",
                "species_color",
                "clade_assignment",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(example_tip_rows())
    with (outdir / "example_clades.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["clade_key", "clade_name", "label_color", "tf_text", "show", "anchor_node_id"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(example_clade_rows())


def read_tree(path: Path):
    tree = Phylo.read(str(path), "newick")
    ensure_internal_ids(tree)
    return tree


def ensure_internal_ids(tree) -> None:
    """Give unnamed internal clades stable IDs so support values remain addressable."""
    existing_names = {
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
                while f"Node{counter}" in existing_names:
                    counter += 1
                clade.name = f"Node{counter}"
                existing_names.add(clade.name)
        if clade is not tree.root:
            counter += 1


def build_parent_lookup(tree) -> Dict[int, object]:
    parent_lookup: Dict[int, object] = {}
    for parent in tree.find_clades(order="preorder"):
        for child in parent.clades:
            parent_lookup[id(child)] = parent
    return parent_lookup


def assign_y_positions(tree) -> Dict[int, float]:
    """Place leaves on evenly spaced rows and propagate midpoints upward."""
    y_positions: Dict[int, float] = {}
    leaves = tree.get_terminals()
    for index, leaf in enumerate(leaves):
        y_positions[id(leaf)] = float(index)

    def visit(clade) -> float:
        if clade.is_terminal():
            return y_positions[id(clade)]
        child_y = [visit(child) for child in clade.clades]
        y_positions[id(clade)] = (min(child_y) + max(child_y)) / 2.0
        return y_positions[id(clade)]

    visit(tree.root)
    return y_positions


def assign_x_positions(tree) -> Dict[int, float]:
    """Accumulate branch lengths from the root, defaulting missing lengths to 1."""
    x_positions: Dict[int, float] = {}

    def visit(clade, current_x: float) -> None:
        x_positions[id(clade)] = current_x
        for child in clade.clades:
            branch_length = child.branch_length if child.branch_length is not None else 1.0
            visit(child, current_x + float(branch_length))

    visit(tree.root, 0.0)
    return x_positions


def scale_x_positions(x_positions: Dict[int, float], tree_width: float) -> Dict[int, float]:
    max_x = max(x_positions.values()) if x_positions else 1.0
    if max_x == 0:
        max_x = 1.0
    return {key: (value / max_x) * tree_width for key, value in x_positions.items()}


def wrap_tf_text(text: str, width: int) -> str:
    if not text:
        return ""
    return textwrap.fill(text, width=width, break_long_words=False, break_on_hyphens=False)


def estimate_half_height(display_name: str, tf_text: str) -> float:
    text_block = display_name + ("\n" if display_name and tf_text else "") + tf_text
    line_count = text_block.count("\n") + 1 if text_block else 1
    return max(0.6, 0.38 * line_count)


def build_support_lookup(tree, clade_annotations: Dict[str, CladeAnnotation]) -> Dict[str, float]:
    support_lookup: Dict[str, float] = {}
    for clade in tree.find_clades(order="preorder"):
        confidence = getattr(clade, "confidence", None)
        if confidence is not None and clade.name:
            support_lookup[clade.name] = float(confidence)
    for annotation in clade_annotations.values():
        if annotation.anchor_node_id and annotation.support is not None:
            support_lookup[annotation.anchor_node_id] = annotation.support
    return support_lookup


def node_leaf_sets(tree) -> Dict[str, set[str]]:
    """Return descendant leaf sets for every internal node with a stable id."""
    leaf_sets: Dict[str, set[str]] = {}
    for clade in tree.find_clades(order="preorder"):
        if clade.is_terminal() or not clade.name:
            continue
        leaf_sets[clade.name] = {leaf.name for leaf in clade.get_terminals()}
    return leaf_sets


def derive_node_tf_text(
    clade_annotations: Dict[str, CladeAnnotation],
    marker_records: Sequence[MarkerRecord],
    leaf_sets_by_node: Dict[str, set[str]],
    node_probability_min: float,
    outside_probability_max: float,
    max_node_tfs: int,
) -> Dict[str, str]:
    """Derive shared descendant TF labels for selected internal nodes.

    The score is based on ancestral-state support within the selected node while
    penalizing TFs that also have high support in disjoint outside nodes. This
    approximates a shared descendant / putative ancestral TF program rather than
    a simple clade-enriched heatmap hit list.
    """
    records_by_node: Dict[str, List[MarkerRecord]] = {}
    records_by_gene: Dict[str, List[MarkerRecord]] = {}
    for record in marker_records:
        records_by_node.setdefault(record.node_id, []).append(record)
        records_by_gene.setdefault(record.orthogroup, []).append(record)

    derived: Dict[str, str] = {}
    for clade_key, annotation in clade_annotations.items():
        if annotation.tf_text.strip():
            derived[clade_key] = annotation.tf_text.strip()
            continue

        anchor_node_id = annotation.anchor_node_id or clade_key
        selected_leaf_set = leaf_sets_by_node.get(anchor_node_id)
        if not selected_leaf_set:
            derived[clade_key] = ""
            continue

        candidates = []
        for record in records_by_node.get(anchor_node_id, []):
            if record.probability < node_probability_min:
                continue

            outside_max = 0.0
            for other_record in records_by_gene.get(record.orthogroup, []):
                if other_record.node_id == anchor_node_id:
                    continue
                other_leaf_set = leaf_sets_by_node.get(other_record.node_id)
                if not other_leaf_set or not other_leaf_set.isdisjoint(selected_leaf_set):
                    continue
                outside_max = max(outside_max, other_record.probability)

            if outside_max > outside_probability_max:
                continue

            specificity_score = record.probability - outside_max
            candidates.append((specificity_score, record.probability, record.orthogroup_name))

        unique_labels: List[str] = []
        seen = set()
        for _, _, label in sorted(candidates, key=lambda item: (-item[0], -item[1], item[2])):
            if not label or label in seen:
                continue
            unique_labels.append(label)
            seen.add(label)
            if len(unique_labels) >= max_node_tfs:
                break

        derived[clade_key] = ", ".join(unique_labels)
    return derived


def compute_clade_placements(
    tree,
    tip_metadata: Dict[str, TipMetadata],
    clade_annotations: Dict[str, CladeAnnotation],
    y_positions: Dict[int, float],
    support_lookup: Dict[str, float],
) -> Dict[str, CladePlacement]:
    """Anchor annotation blocks by the leaves assigned to each biological module."""
    leaves_by_name = {leaf.name: leaf for leaf in tree.get_terminals()}
    clades_by_name = {clade.name: clade for clade in tree.find_clades(order="preorder")}
    placements: Dict[str, CladePlacement] = {}
    tips_by_clade: Dict[str, List[object]] = {}

    for clade_key in clade_annotations:
        normalized_key = normalize_clade_key(clade_key)
        assigned = [
            leaves_by_name[leaf_id]
            for leaf_id, metadata in tip_metadata.items()
            if leaf_id in leaves_by_name and normalize_clade_key(metadata.clade_assignment) == normalized_key
        ]
        tips_by_clade[clade_key] = assigned

    for clade_key, annotation in clade_annotations.items():
        assigned_leaves = tips_by_clade.get(clade_key, [])
        anchor_clade = None
        if annotation.anchor_node_id:
            anchor_clade = clades_by_name.get(annotation.anchor_node_id)
        if anchor_clade is None and assigned_leaves:
            anchor_clade = tree.common_ancestor(assigned_leaves)
        if anchor_clade is None:
            continue
        if not assigned_leaves:
            assigned_leaves = anchor_clade.get_terminals()
        assigned_y = sorted(y_positions[id(leaf)] for leaf in assigned_leaves)
        target_y = (assigned_y[0] + assigned_y[-1]) / 2.0
        placements[clade_key] = CladePlacement(
            clade_key=clade_key,
            anchor_clade=anchor_clade,
            target_y=target_y,
            support=support_lookup.get(anchor_clade.name) if anchor_clade.name else None,
        )
    return placements


def validate_annotation_inputs(
    tip_metadata: Dict[str, TipMetadata],
    clade_annotations: Dict[str, CladeAnnotation],
    clade_placements: Dict[str, CladePlacement],
) -> None:
    """Fail with a useful message when no biological clades can be matched."""
    visible_keys = [key for key, annotation in clade_annotations.items() if annotation.show]
    visible_placed = [key for key in visible_keys if key in clade_placements]
    if visible_keys and not visible_placed:
        tip_clades = sorted({meta.clade_assignment for meta in tip_metadata.values() if meta.clade_assignment})
        clade_keys = sorted(visible_keys)
        raise ValueError(
            "No visible TF annotation nodes/clades could be placed. "
            "This usually means either the tip_metadata 'clade_assignment' values do not match the "
            "clade annotation 'clade_key' values, or the selected 'anchor_node_id' values are missing "
            "from the Newick tree. "
            f"Tip clades: {tip_clades[:10]}. Clade keys: {clade_keys[:10]}."
        )


def resolve_annotation_y_positions(
    clade_annotations: Dict[str, CladeAnnotation],
    clade_placements: Dict[str, CladePlacement],
    wrapped_tf: Dict[str, str],
    lower: float,
    upper: float,
    gap: float,
) -> Dict[str, float]:
    """Push the right-hand annotation blocks apart while preserving clade order."""
    visible = [
        clade_annotations[key]
        for key in clade_annotations
        if clade_annotations[key].show and key in clade_placements
    ]
    ordered = sorted(visible, key=lambda annotation: clade_placements[annotation.clade_key].target_y)
    if not ordered:
        return {}

    positions: Dict[str, float] = {}
    half_heights = {
        annotation.clade_key: estimate_half_height(annotation.display_name, wrapped_tf.get(annotation.clade_key, ""))
        for annotation in ordered
    }

    previous_y = None
    previous_half = 0.0
    for annotation in ordered:
        clade_key = annotation.clade_key
        half_height = half_heights[clade_key]
        minimum_y = lower + half_height
        if previous_y is not None:
            minimum_y = max(minimum_y, previous_y + previous_half + half_height + gap)
        proposed = max(clade_placements[clade_key].target_y, minimum_y)
        positions[clade_key] = proposed
        previous_y = proposed
        previous_half = half_height

    next_y = None
    next_half = 0.0
    for annotation in reversed(ordered):
        clade_key = annotation.clade_key
        half_height = half_heights[clade_key]
        maximum_y = upper - half_height
        if next_y is not None:
            maximum_y = min(maximum_y, next_y - next_half - half_height - gap)
        positions[clade_key] = min(positions[clade_key], maximum_y)
        next_y = positions[clade_key]
        next_half = half_height

    return positions


def draw_branches(ax, tree, x_positions: Dict[int, float], y_positions: Dict[int, float], tip_x: float) -> None:
    """Draw the rooted tree with restrained branch styling."""
    branch_color = "#3d3d3d"
    extension_color = "#8f8f8f"
    for clade in tree.find_clades(order="preorder"):
        x_here = x_positions[id(clade)]
        if clade.clades:
            child_ys = [y_positions[id(child)] for child in clade.clades]
            ax.plot(
                [x_here, x_here],
                [min(child_ys), max(child_ys)],
                color=branch_color,
                linewidth=0.85,
                solid_capstyle="round",
                zorder=1,
            )
        for child in clade.clades:
            child_x = x_positions[id(child)]
            child_y = y_positions[id(child)]
            ax.plot(
                [x_here, child_x],
                [child_y, child_y],
                color=branch_color,
                linewidth=0.85,
                solid_capstyle="round",
                zorder=1,
            )
            if child.is_terminal() and child_x < tip_x:
                ax.plot(
                    [child_x, tip_x],
                    [child_y, child_y],
                    color=extension_color,
                    linewidth=0.55,
                    solid_capstyle="round",
                    zorder=1,
                )


def draw_support_values(
    ax,
    tree,
    parent_lookup: Dict[int, object],
    x_positions: Dict[int, float],
    y_positions: Dict[int, float],
    support_lookup: Dict[str, float],
    support_min: float,
) -> None:
    """Place support values on internal branches without dominating the figure."""
    for clade in tree.find_clades(order="preorder"):
        if clade.is_terminal() or clade is tree.root or not clade.name:
            continue
        support = support_lookup.get(clade.name)
        if support is None or support < support_min:
            continue
        parent = parent_lookup.get(id(clade))
        x_parent = x_positions[id(parent)] if parent is not None else x_positions[id(clade)]
        x_here = x_positions[id(clade)]
        x_text = x_parent + (x_here - x_parent) * 0.52
        y_text = y_positions[id(clade)] - 0.34
        label = f"{int(round(support))}" if abs(support - round(support)) < 1e-6 else f"{support:.1f}"
        ax.text(
            x_text,
            y_text,
            label,
            fontsize=7.1,
            color="#595959",
            ha="center",
            va="bottom",
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.12},
            zorder=3,
        )


def add_tip_markers_and_labels(
    ax,
    leaves: Sequence[object],
    tip_metadata: Dict[str, TipMetadata],
    y_positions: Dict[int, float],
    tip_x: float,
    species_x: float,
    label_x: float,
) -> None:
    """Draw full leaf labels as the primary middle column and species cues as secondary."""
    for leaf in leaves:
        metadata = tip_metadata[leaf.name]
        y_here = y_positions[id(leaf)]
        ax.scatter(
            [tip_x],
            [y_here],
            s=24,
            c=[metadata.tip_color],
            edgecolors="white",
            linewidths=0.5,
            zorder=3,
        )
        ax.text(
            species_x,
            y_here,
            metadata.species_abbrev,
            fontsize=7.6,
            color=metadata.species_color,
            ha="left",
            va="center",
            fontweight="semibold",
            zorder=3,
        )
        ax.text(
            label_x,
            y_here,
            metadata.leaf_display_label,
            fontsize=8.5,
            color="#1f1f1f",
            ha="left",
            va="center",
            zorder=3,
        )


def add_clade_annotations(
    ax,
    clade_annotations: Dict[str, CladeAnnotation],
    clade_placements: Dict[str, CladePlacement],
    x_positions: Dict[int, float],
    y_positions: Dict[int, float],
    annotation_y: Dict[str, float],
    connector_x: float,
    clade_x: float,
    tf_x: float,
    wrapped_tf: Dict[str, str],
) -> None:
    """Draw colored node/clade labels and the aligned TF program column."""
    ordered_keys = sorted(annotation_y, key=lambda key: annotation_y[key])
    for clade_key in ordered_keys:
        annotation = clade_annotations[clade_key]
        placement = clade_placements[clade_key]
        y_anchor = annotation_y[clade_key]
        anchor_x = x_positions[id(placement.anchor_clade)]
        anchor_y = y_positions[id(placement.anchor_clade)]

        ax.plot(
            [anchor_x, connector_x],
            [anchor_y, y_anchor],
            color="#c7c7c7",
            linewidth=0.65,
            zorder=1,
        )
        ax.plot(
            [connector_x, clade_x - 0.8],
            [y_anchor, y_anchor],
            color="#d6d6d6",
            linewidth=0.55,
            zorder=1,
        )

        ax.text(
            clade_x,
            y_anchor,
            annotation.display_name,
            fontsize=10.1,
            color=annotation.label_color,
            ha="left",
            va="center",
            fontweight="semibold",
            zorder=3,
        )
        tf_text = wrapped_tf.get(clade_key, "")
        if tf_text:
            ax.text(
                tf_x,
                y_anchor,
                tf_text,
                fontsize=9.2,
                color="#2a2a2a",
                ha="left",
                va="center",
                linespacing=1.16,
                zorder=3,
            )


def render_tree(
    tree,
    tip_metadata: Dict[str, TipMetadata],
    clade_annotations: Dict[str, CladeAnnotation],
    external_node_support: Dict[str, float],
    marker_records: Sequence[MarkerRecord],
    output_prefix: Path,
    figure_width: float,
    leaf_height: float,
    min_height: float,
    support_min: float,
    tf_wrap: int,
    annotation_gap: float,
    png_dpi: int,
    node_probability_min: float,
    outside_probability_max: float,
    max_node_tfs: int,
) -> None:
    """Compose the final figure and save PDF, SVG, and high-resolution PNG."""
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["DejaVu Sans", "Arial", "Liberation Sans"],
            "axes.facecolor": "white",
            "figure.facecolor": "white",
        }
    )

    leaves = tree.get_terminals()
    missing_tips = sorted({leaf.name for leaf in leaves} - set(tip_metadata))
    if missing_tips:
        raise ValueError(f"Missing tip metadata for: {', '.join(missing_tips)}")

    parent_lookup = build_parent_lookup(tree)
    y_positions = assign_y_positions(tree)
    x_positions = scale_x_positions(assign_x_positions(tree), tree_width=28.5)
    support_lookup = build_support_lookup(tree, clade_annotations)
    support_lookup.update(external_node_support)
    clade_placements = compute_clade_placements(
        tree,
        tip_metadata=tip_metadata,
        clade_annotations=clade_annotations,
        y_positions=y_positions,
        support_lookup=support_lookup,
    )
    validate_annotation_inputs(tip_metadata, clade_annotations, clade_placements)
    leaf_sets_by_node = node_leaf_sets(tree)
    tf_text_by_clade = derive_node_tf_text(
        clade_annotations=clade_annotations,
        marker_records=marker_records,
        leaf_sets_by_node=leaf_sets_by_node,
        node_probability_min=node_probability_min,
        outside_probability_max=outside_probability_max,
        max_node_tfs=max_node_tfs,
    )

    wrapped_tf = {
        key: wrap_tf_text(tf_text_by_clade.get(key, annotation.tf_text), width=tf_wrap)
        for key, annotation in clade_annotations.items()
    }
    annotation_y = resolve_annotation_y_positions(
        clade_annotations=clade_annotations,
        clade_placements=clade_placements,
        wrapped_tf=wrapped_tf,
        lower=-0.5,
        upper=len(leaves) - 0.5,
        gap=annotation_gap,
    )

    tip_x = 31.0
    species_x = 32.8
    label_x = 35.8
    connector_x = 58.0
    clade_x = 60.2
    tf_x = 73.5
    x_max = 128.0

    figure_height = max(min_height, leaf_height * len(leaves) + 2.3)
    fig, ax = plt.subplots(figsize=(figure_width, figure_height))

    draw_branches(ax, tree, x_positions=x_positions, y_positions=y_positions, tip_x=tip_x)
    draw_support_values(
        ax,
        tree,
        parent_lookup=parent_lookup,
        x_positions=x_positions,
        y_positions=y_positions,
        support_lookup=support_lookup,
        support_min=support_min,
    )
    add_tip_markers_and_labels(
        ax,
        leaves=leaves,
        tip_metadata=tip_metadata,
        y_positions=y_positions,
        tip_x=tip_x,
        species_x=species_x,
        label_x=label_x,
    )
    add_clade_annotations(
        ax,
        clade_annotations=clade_annotations,
        clade_placements=clade_placements,
        x_positions=x_positions,
        y_positions=y_positions,
        annotation_y=annotation_y,
        connector_x=connector_x,
        clade_x=clade_x,
        tf_x=tf_x,
        wrapped_tf=wrapped_tf,
    )

    ax.set_xlim(-1.0, x_max)
    ax.set_ylim(len(leaves) - 0.5, -0.5)
    ax.axis("off")
    fig.tight_layout()

    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(add_extension(output_prefix, "pdf"), bbox_inches="tight")
    fig.savefig(add_extension(output_prefix, "svg"), bbox_inches="tight")
    fig.savefig(add_extension(output_prefix, "png"), dpi=png_dpi, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = resolve_pipeline_paths(parse_args())

    if args.write_example:
        write_example_dataset(args.write_example)
        return

    if args.example:
        tree, tip_metadata, clade_annotations = load_example_dataset()
        node_support = {}
        marker_records = []
    else:
        if not args.tree or not args.tips or not args.clades:
            raise SystemExit("Provide --tree, --tips, and --clades, use --results-dir, or use --example.")
        missing = [str(path) for path in (args.tree, args.tips, args.clades) if path is not None and not path.exists()]
        if missing:
            raise SystemExit(f"Required input file(s) not found: {', '.join(missing)}")
        tree = read_tree(args.tree)
        tip_metadata = load_tip_metadata(read_tsv_rows(args.tips))
        clade_annotations = load_clade_annotations(read_tsv_rows(args.clades))
        node_support = load_node_support(read_tsv_rows(args.node_support)) if args.node_support else {}
        marker_records = load_marker_records(read_tsv_rows(args.markers)) if args.markers else []

    render_tree(
        tree=tree,
        tip_metadata=tip_metadata,
        clade_annotations=clade_annotations,
        external_node_support=node_support,
        marker_records=marker_records,
        output_prefix=args.output_prefix,
        figure_width=args.figure_width,
        leaf_height=args.leaf_height,
        min_height=args.min_height,
        support_min=args.support_min,
        tf_wrap=args.tf_wrap,
        annotation_gap=args.annotation_gap,
        png_dpi=args.png_dpi,
        node_probability_min=args.node_probability_min,
        outside_probability_max=args.outside_probability_max,
        max_node_tfs=args.max_node_tfs,
    )


if __name__ == "__main__":
    main()
