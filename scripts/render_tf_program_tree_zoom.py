#!/usr/bin/env python3
"""Render a zoomed companion panel for selected node-based TF tree annotations.

This script reuses the main manuscript-style tree renderer logic but crops the
vertical extent to the band containing the selected annotated nodes. The result
is intended to complement the full-tree overview: the overview shows the whole
cross-species dendrogram, while the zoom panel makes the selected annotated
region legible at manuscript scale.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt

from render_tf_program_tree import (
    add_clade_annotations,
    add_extension,
    add_tip_markers_and_labels,
    assign_x_positions,
    assign_y_positions,
    build_parent_lookup,
    build_support_lookup,
    compute_clade_placements,
    derive_node_tf_text,
    draw_branches,
    draw_support_values,
    load_clade_annotations,
    load_marker_records,
    load_node_support,
    load_tip_metadata,
    node_leaf_sets,
    read_tree,
    read_tsv_rows,
    resolve_annotation_y_positions,
    scale_x_positions,
    validate_annotation_inputs,
    wrap_tf_text,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tree", type=Path, required=True, help="Rooted Newick tree.")
    parser.add_argument("--tips", type=Path, required=True, help="Tip metadata TSV.")
    parser.add_argument("--clades", type=Path, required=True, help="Selected node/clade annotation TSV.")
    parser.add_argument("--node-support", type=Path, required=True, help="Optional node support TSV.")
    parser.add_argument("--markers", type=Path, required=True, help="UPGMA ancestral TF marker TSV.")
    parser.add_argument("--output-prefix", type=Path, required=True, help="Output prefix for the zoom panel.")
    parser.add_argument("--padding-leaves", type=float, default=8.0, help="Vertical padding around the selected annotation band.")
    parser.add_argument("--figure-width", type=float, default=12.5, help="Figure width in inches.")
    parser.add_argument("--leaf-height", type=float, default=0.19, help="Height per visible leaf in inches.")
    parser.add_argument("--min-height", type=float, default=7.5, help="Minimum figure height in inches.")
    parser.add_argument("--support-min", type=float, default=20.0, help="Minimum support value to draw.")
    parser.add_argument("--tf-wrap", type=int, default=28, help="Approximate wrap width for TF labels.")
    parser.add_argument("--annotation-gap", type=float, default=0.8, help="Minimum vertical gap between annotation blocks.")
    parser.add_argument("--png-dpi", type=int, default=300, help="PNG export DPI.")
    parser.add_argument("--node-probability-min", type=float, default=0.7, help="Minimum node probability for derived TF labels.")
    parser.add_argument("--outside-probability-max", type=float, default=0.4, help="Maximum probability allowed in disjoint outside nodes.")
    parser.add_argument("--max-node-tfs", type=int, default=4, help="Maximum number of derived TF labels per selected node.")
    return parser.parse_args()


def render_zoom_panel(args: argparse.Namespace) -> None:
    tree = read_tree(args.tree)
    tip_metadata = load_tip_metadata(read_tsv_rows(args.tips))
    clade_annotations = load_clade_annotations(read_tsv_rows(args.clades))
    node_support = load_node_support(read_tsv_rows(args.node_support))
    marker_records = load_marker_records(read_tsv_rows(args.markers))

    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["DejaVu Sans", "Arial", "Liberation Sans"],
            "axes.facecolor": "white",
            "figure.facecolor": "white",
        }
    )

    leaves = tree.get_terminals()
    parent_lookup = build_parent_lookup(tree)
    y_positions = assign_y_positions(tree)
    x_positions = scale_x_positions(assign_x_positions(tree), tree_width=28.5)
    support_lookup = build_support_lookup(tree, clade_annotations)
    support_lookup.update(node_support)
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
        node_probability_min=args.node_probability_min,
        outside_probability_max=args.outside_probability_max,
        max_node_tfs=args.max_node_tfs,
    )
    wrapped_tf = {
        key: wrap_tf_text(tf_text_by_clade.get(key, annotation.tf_text), width=args.tf_wrap)
        for key, annotation in clade_annotations.items()
    }
    annotation_y = resolve_annotation_y_positions(
        clade_annotations=clade_annotations,
        clade_placements=clade_placements,
        wrapped_tf=wrapped_tf,
        lower=-0.5,
        upper=len(leaves) - 0.5,
        gap=args.annotation_gap,
    )

    if not annotation_y:
        raise SystemExit("No visible annotations were placed; cannot build a zoom panel.")

    min_y = min(annotation_y.values())
    max_y = max(annotation_y.values())
    lower = max(-0.5, min_y - args.padding_leaves)
    upper = min(len(leaves) - 0.5, max_y + args.padding_leaves)
    visible_span = upper - lower

    tip_x = 31.0
    species_x = 32.8
    label_x = 35.8
    connector_x = 58.0
    clade_x = 60.2
    tf_x = 73.5
    x_max = 128.0

    figure_height = max(args.min_height, args.leaf_height * visible_span + 2.5)
    fig, ax = plt.subplots(figsize=(args.figure_width, figure_height))

    draw_branches(ax, tree, x_positions=x_positions, y_positions=y_positions, tip_x=tip_x)
    draw_support_values(
        ax,
        tree,
        parent_lookup=parent_lookup,
        x_positions=x_positions,
        y_positions=y_positions,
        support_lookup=support_lookup,
        support_min=args.support_min,
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
    ax.set_ylim(upper, lower)
    ax.axis("off")
    fig.tight_layout()

    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(add_extension(args.output_prefix, "pdf"), bbox_inches="tight")
    fig.savefig(add_extension(args.output_prefix, "svg"), bbox_inches="tight")
    fig.savefig(add_extension(args.output_prefix, "png"), dpi=args.png_dpi, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    render_zoom_panel(parse_args())


if __name__ == "__main__":
    main()
