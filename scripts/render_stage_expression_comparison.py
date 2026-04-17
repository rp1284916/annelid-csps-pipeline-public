#!/usr/bin/env python3
"""Render developmental stage-expression comparisons for a cross-species orthogroup.

This script is designed for the annelid raw inputs used by the CSPS pipeline.
It reads the original h5ad objects, detects a stage/time column in `.obs`,
aggregates expression by stage, and plots species-specific developmental
expression profiles for a target orthogroup or explicit gene pair.

The intended interpretation is onset relative to a user-defined developmental
milestone within each species. The figure therefore uses one developmental axis
per species rather than implying stage-for-stage equivalence between lineages.

Example:
    python scripts/render_stage_expression_comparison.py \
      --config config/annelid_inputs.json \
      --orthogroup OG_PLACEHOLDER \
      --output-prefix work/results_stage_expression/example_signal
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Iterable

try:
    import anndata as ad
except ImportError as exc:  # pragma: no cover - import guard
    raise SystemExit(
        "anndata is required for this script. Install it with 'pip install anndata h5py pandas matplotlib'."
    ) from exc

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError as exc:  # pragma: no cover - import guard
    raise SystemExit(
        "matplotlib is required for this script. Install it with 'pip install matplotlib'."
    ) from exc

import numpy as np
import pandas as pd


DEFAULT_STAGE_COLUMNS = (
    "stage",
    "Stage",
    "timepoint",
    "Timepoint",
    "time_point",
    "age",
    "Age",
    "sample",
    "Sample",
    "batch",
    "Batch",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True, help="Pipeline config JSON with raw_data_dir and species h5ad paths.")
    parser.add_argument("--orthogroup", default="OG_PLACEHOLDER", help="Orthogroup to render.")
    parser.add_argument("--ctel-gene", default="CTELG_PLACEHOLDER", help="Primary Capitella gene to highlight.")
    parser.add_argument("--ofus-gene", default="OFUSG_PLACEHOLDER", help="Primary Owenia gene to highlight.")
    parser.add_argument("--ctel-milestone", default="MILESTONE_A", help="Capitella milestone label to mark on the developmental axis.")
    parser.add_argument("--ofus-milestone", default="MILESTONE_B", help="Owenia milestone label to mark on the developmental axis.")
    parser.add_argument("--stage-column-ctel", default="", help="Optional explicit stage column for Ctel.")
    parser.add_argument("--stage-column-ofus", default="", help="Optional explicit stage column for Ofus.")
    parser.add_argument("--expression-layer", default="", help="Optional AnnData layer name; defaults to X.")
    parser.add_argument("--output-prefix", type=Path, required=True, help="Prefix for PDF/SVG/PNG outputs.")
    return parser.parse_args()


def load_config(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def detect_stage_column(obs: pd.DataFrame, explicit: str = "") -> str:
    """Choose the obs column that best represents developmental stage/time."""
    if explicit:
        if explicit not in obs.columns:
            raise ValueError(f"Requested stage column '{explicit}' not found. Available columns: {list(obs.columns)}")
        return explicit

    lower_lookup = {col.lower(): col for col in obs.columns}
    for candidate in DEFAULT_STAGE_COLUMNS:
        if candidate in obs.columns:
            return candidate
        if candidate.lower() in lower_lookup:
            return lower_lookup[candidate.lower()]

    fuzzy_hits = [
        col for col in obs.columns
        if any(token in col.lower() for token in ("stage", "time", "age", "hpf", "dpf", "sample", "batch"))
    ]
    if len(fuzzy_hits) == 1:
        return fuzzy_hits[0]

    raise ValueError(
        "Could not detect a stage/time column automatically. "
        f"Available obs columns: {list(obs.columns)}"
    )


def load_orthogroup_members(orthogroup_path: Path, orthogroup: str, gene_prefix: str) -> list[str]:
    df = pd.read_csv(orthogroup_path, sep="\t", header=None, names=["orthogroup", "gene"], dtype=str)
    members = df.loc[df["orthogroup"] == orthogroup, "gene"].dropna().tolist()
    return [gene for gene in members if gene.startswith(gene_prefix)]


def var_names_upper(adata: ad.AnnData) -> pd.Index:
    return pd.Index([str(v).upper() for v in adata.var_names])


def extract_gene_vector(adata: ad.AnnData, gene_ids: Iterable[str], layer: str = "") -> pd.Series:
    gene_ids = [str(g).upper() for g in gene_ids]
    var_index = var_names_upper(adata)
    found = []
    for gene_id in gene_ids:
        matches = np.where((var_index == gene_id) | var_index.str.startswith(f"{gene_id}-"))[0]
        if matches.size:
            found.extend(matches.tolist())
    if not found:
        raise ValueError(f"None of the requested genes were found in var_names: {gene_ids}")
    subset = adata[:, found]
    if hasattr(subset, "to_memory"):
        subset = subset.to_memory()
    if layer:
        matrix = subset.layers[layer]
    else:
        matrix = subset.X
    if hasattr(matrix, "toarray"):
        matrix = matrix.toarray()
    matrix = np.asarray(matrix)
    values = matrix.sum(axis=1)
    return pd.Series(values, index=adata.obs_names, dtype=float)


def simplify_stage_label(label: str) -> str:
    label = str(label)
    label = label.strip()
    label = re.sub(r"\s+", " ", label)
    return label


def stage_sort_key(label: str):
    label_s = simplify_stage_label(label)
    match = re.search(r"(\d+(?:\.\d+)?)\s*(hpf|dpf)", label_s.lower())
    if match:
        value = float(match.group(1))
        unit = match.group(2)
        if unit == "dpf":
            value *= 24
        return (0, value, label_s.lower())
    day_match = re.search(r"day\s*(\d+(?:\.\d+)?)", label_s.lower())
    if day_match:
        return (0, float(day_match.group(1)) * 24, label_s.lower())
    stage_match = re.search(r"st\s*(\d+(?:\.\d+)?)", label_s.lower())
    if stage_match:
        return (0, float(stage_match.group(1)), label_s.lower())
    if re.search(r"embry", label_s.lower()):
        return (1, 0, label_s.lower())
    if re.search(r"gastr", label_s.lower()):
        return (1, 1, label_s.lower())
    if re.search(r"larv", label_s.lower()):
        return (1, 2, label_s.lower())
    if re.search(r"juven", label_s.lower()):
        return (1, 3, label_s.lower())
    return (2, label_s.lower())


def summarise_by_stage(adata: ad.AnnData, stage_column: str, gene_ids: list[str], layer: str = "") -> pd.DataFrame:
    """Aggregate orthogroup-level signal within each developmental stage label."""
    expr = extract_gene_vector(adata, gene_ids, layer=layer)
    df = adata.obs.copy()
    df["_expr"] = expr.reindex(df.index).values
    df["_stage"] = df[stage_column].astype(str).map(simplify_stage_label)
    grouped = (
        df.groupby("_stage", dropna=False)
        .agg(mean_expression=("_expr", "mean"), median_expression=("_expr", "median"), n_cells=("_expr", "size"))
        .reset_index()
        .rename(columns={"_stage": "stage"})
    )
    grouped["relative_expression"] = (
        grouped["mean_expression"] / grouped["mean_expression"].max()
        if grouped["mean_expression"].max() > 0
        else 0
    )
    grouped = grouped.sort_values(by="stage", key=lambda s: s.map(stage_sort_key)).reset_index(drop=True)
    return grouped


def build_output_table(ctel_df: pd.DataFrame, ofus_df: pd.DataFrame) -> pd.DataFrame:
    """Write the summarised stage series back out in a flat, species-tagged table."""
    ctel_table = ctel_df.copy()
    ctel_table.insert(0, "species", "Ctel")
    ctel_table.insert(1, "developmental_order", np.arange(1, len(ctel_table) + 1))

    ofus_table = ofus_df.copy()
    ofus_table.insert(0, "species", "Ofus")
    ofus_table.insert(1, "developmental_order", np.arange(1, len(ofus_table) + 1))

    return pd.concat([ctel_table, ofus_table], ignore_index=True)


def find_stage_index(df: pd.DataFrame, stage_label: str) -> int | None:
    """Locate a named stage so a milestone marker can be added to the panel."""
    normalized_target = simplify_stage_label(stage_label).lower()
    for idx, label in enumerate(df["stage"].tolist()):
        if simplify_stage_label(label).lower() == normalized_target:
            return idx
    return None


def save_all_formats(fig: plt.Figure, output_prefix: Path) -> None:
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "svg", "png"):
        fig.savefig(output_prefix.with_suffix(f".{ext}"), dpi=300, bbox_inches="tight", facecolor="white")


def render_figure(
    ctel_df: pd.DataFrame,
    ofus_df: pd.DataFrame,
    orthogroup: str,
    ctel_label: str,
    ofus_label: str,
    ctel_milestone: str,
    ofus_milestone: str,
    output_prefix: Path,
) -> None:
    """Draw one species-specific panel per lineage with milestone markers."""
    ctel_x = np.arange(len(ctel_df))
    ofus_x = np.arange(len(ofus_df))

    fig, axes = plt.subplots(
        2,
        1,
        figsize=(max(10, max(len(ctel_df), len(ofus_df)) * 1.6), 9.6),
        constrained_layout=False,
    )
    fig.subplots_adjust(top=0.91, bottom=0.14, hspace=0.22)
    species_colors = {"Ctel": "#1f78b4", "Ofus": "#d95f02"}
    milestone_line_color = "#666666"

    species_panels = (
        (axes[0], "Ctel", ctel_df, ctel_x, ctel_label, ctel_milestone),
        (axes[1], "Ofus", ofus_df, ofus_x, ofus_label, ofus_milestone),
    )

    for ax, species, df, xvals, gene_label, milestone_label in species_panels:
        # Relative expression emphasises profile shape and timing within each
        # species while avoiding a misleading absolute scale comparison.
        ax.plot(
            xvals,
            df["relative_expression"],
            marker="o",
            lw=2.5,
            color=species_colors[species],
            label=gene_label,
            zorder=3,
        )
        ax.set_ylabel("Relative expression")
        ax.set_ylim(-0.02, 1.08)
        ax.axhline(1, color="#d9d9d9", lw=1, ls="--", zorder=1)
        ax.set_title(f"{species} developmental profile", loc="left", color=species_colors[species], fontsize=12)
        ax.set_xticks(xvals)
        ax.set_xticklabels(df["stage"].tolist())
        ax.set_xlabel(f"{species} developmental stage / time")

        onset_idx = find_stage_index(df, milestone_label)
        if onset_idx is not None:
            # The dashed line anchors the interpretation to the user-defined
            # developmental milestone for that species.
            ax.axvline(onset_idx, color=milestone_line_color, lw=1.2, ls=(0, (4, 3)), zorder=2)
            ax.text(
                onset_idx + 0.05,
                1.03,
                f"Milestone\n({milestone_label})",
                ha="left",
                va="top",
                fontsize=9,
                color=milestone_line_color,
            )

        for xpos, ypos in zip(xvals, df["relative_expression"].tolist()):
            ax.scatter([xpos], [ypos], s=34, color=species_colors[species], edgecolor="white", linewidth=0.6, zorder=4)

    for ax in axes:
        ax.set_facecolor("white")
        ax.grid(axis="y", color="#e5e5e5", lw=0.8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    fig.text(
        0.01,
        0.035,
        "Stages are shown on species-specific developmental scales; interpret expression onset relative to the marked milestone within each species, not as direct stage equivalence across species.",
        fontsize=9,
        color="#404040",
    )
    fig.suptitle(f"{orthogroup} expression onset relative to developmental milestones", fontsize=14, y=0.975)
    save_all_formats(fig, output_prefix)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    cfg = load_config(args.config)
    raw_data_dir = Path(cfg["raw_data_dir"])
    processed_dir = args.config.parent.parent / cfg["processed_dir"]

    ctel_h5ad = raw_data_dir / cfg["species"]["Ctel"]["h5ad"]
    ofus_h5ad = raw_data_dir / cfg["species"]["Ofus"]["h5ad"]
    orthology_ctel = processed_dir / "orthology" / "Ctel.gene_orthogroups.tsv"
    orthology_ofus = processed_dir / "orthology" / "Ofus.gene_orthogroups.tsv"

    ctel_members = load_orthogroup_members(orthology_ctel, args.orthogroup, cfg["species"]["Ctel"]["gene_prefix"])
    ofus_members = load_orthogroup_members(orthology_ofus, args.orthogroup, cfg["species"]["Ofus"]["gene_prefix"])

    # Keep the requested gene first for labelling, but include all orthogroup
    # members in the aggregated stage signal to reflect the whole OG.
    ctel_gene_ids = [args.ctel_gene] + [g for g in ctel_members if g != args.ctel_gene]
    ofus_gene_ids = [args.ofus_gene] + [g for g in ofus_members if g != args.ofus_gene]

    ctel = ad.read_h5ad(ctel_h5ad, backed="r")
    ofus = ad.read_h5ad(ofus_h5ad, backed="r")

    try:
        ctel_stage_col = detect_stage_column(ctel.obs, args.stage_column_ctel)
        ofus_stage_col = detect_stage_column(ofus.obs, args.stage_column_ofus)

        ctel_stage = summarise_by_stage(ctel, ctel_stage_col, ctel_gene_ids, layer=args.expression_layer)
        ofus_stage = summarise_by_stage(ofus, ofus_stage_col, ofus_gene_ids, layer=args.expression_layer)
    finally:
        if getattr(ctel, "isbacked", False):
            ctel.file.close()
        if getattr(ofus, "isbacked", False):
            ofus.file.close()

    render_figure(
        ctel_stage,
        ofus_stage,
        orthogroup=args.orthogroup,
        ctel_label=f"{args.ctel_gene} (+{len(ctel_gene_ids) - 1} OG members)" if len(ctel_gene_ids) > 1 else args.ctel_gene,
        ofus_label=f"{args.ofus_gene} (+{len(ofus_gene_ids) - 1} OG members)" if len(ofus_gene_ids) > 1 else args.ofus_gene,
        ctel_milestone=args.ctel_milestone,
        ofus_milestone=args.ofus_milestone,
        output_prefix=args.output_prefix,
    )

    output_table = args.output_prefix.with_suffix(".stage_expression.tsv")
    merged = build_output_table(ctel_stage, ofus_stage)
    merged.to_csv(output_table, sep="\t", index=False)
    print(f"Wrote stage-expression comparison to {output_table}")
    print(f"Detected stage columns: Ctel={ctel_stage_col}, Ofus={ofus_stage_col}")
    print(f"Milestone markers: Ctel={args.ctel_milestone}, Ofus={args.ofus_milestone}")


if __name__ == "__main__":
    main()
