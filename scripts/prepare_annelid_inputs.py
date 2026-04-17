#!/usr/bin/env python3
"""Prepare annelid-specific inputs for the reduced cross-species pipeline."""

import argparse
import gzip
import json
import re
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import stats


def canonical_gene_id(raw_id: str, gene_prefix: str) -> str:
    """Collapse transcript-style identifiers onto the expected gene identifier prefix."""
    raw_id = str(raw_id).strip()
    gene_match = re.match(rf"^({re.escape(gene_prefix)}\d+)", raw_id)
    if gene_match:
        return gene_match.group(1)
    return raw_id


def canonical_ortholog_id(raw_id: str, gene_prefix: str) -> str:
    """Normalise orthology identifiers before matching them across source files."""
    raw_id = str(raw_id).strip()
    raw_id = re.sub(r"\.\d+$", "", raw_id)
    return canonical_gene_id(raw_id, gene_prefix)


def parse_cluster_map(path: Path) -> dict[str, str]:
    """Parse the hand-curated cluster label mapping embedded in the Excel export."""
    df = pd.read_excel(path, header=None)
    mapping = {}
    for raw in df.iloc[:, 0].dropna().astype(str):
        text = raw.strip()
        match = re.search(r'"?\s*(\d+)\s*"?\s*:\s*"(\d+)\|(.+?)"\s*,?$', text)
        if match:
            mapping[match.group(1)] = match.group(3)
    if not mapping:
        raise ValueError(f"Could not parse cluster labels from {path}")
    return mapping


def benjamini_hochberg(pvalues: np.ndarray) -> np.ndarray:
    """Apply Benjamini-Hochberg correction without introducing an extra dependency."""
    pvalues = np.asarray(pvalues, dtype=float)
    n = pvalues.size
    order = np.argsort(pvalues)
    ranked = pvalues[order]
    adjusted = np.empty(n, dtype=float)
    running = 1.0
    for i in range(n - 1, -1, -1):
        rank = i + 1
        candidate = ranked[i] * n / rank
        running = min(running, candidate)
        adjusted[i] = running
    out = np.empty(n, dtype=float)
    out[order] = np.minimum(adjusted, 1.0)
    return out


def load_species_metadata(raw_dir: Path, species_id: str, cfg: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load metacell labels and colours from the species h5ad file."""
    h5ad_path = raw_dir / cfg["h5ad"]
    label_map = parse_cluster_map(raw_dir / cfg["cluster_labels_xlsx"])
    adata = ad.read_h5ad(h5ad_path, backed="r")

    leiden_series = adata.obs["leiden_6"].astype(str)
    if hasattr(adata.obs["leiden_6"], "cat"):
        categories = [str(v) for v in adata.obs["leiden_6"].cat.categories]
    else:
        categories = sorted(set(leiden_series))
    colors = list(adata.uns["leiden_6_colors"])
    color_map = {cluster: color for cluster, color in zip(categories, colors)}

    cell_metadata = pd.DataFrame(
        {
            "metacell": adata.obs_names.astype(str),
            "leiden_6": leiden_series.values,
        }
    )
    cell_metadata["cell_type"] = cell_metadata["leiden_6"].map(label_map)
    cell_metadata["color"] = cell_metadata["leiden_6"].map(color_map)
    cell_metadata["cell_type_species"] = [f"{species_id}|{ct}" for ct in cell_metadata["cell_type"]]
    cell_metadata.index = cell_metadata["metacell"]

    annot = (
        cell_metadata[["leiden_6", "cell_type", "color"]]
        .drop_duplicates()
        .sort_values("leiden_6", key=lambda s: s.astype(int))
        .reset_index(drop=True)
    )

    adata.file.close()
    return cell_metadata, annot


def load_gene_annotations(path: Path, gene_prefix: str) -> pd.DataFrame:
    """Reduce transcript-level annotation rows to one record per canonical gene."""
    annot = pd.read_csv(path, sep="\t")
    annot["gene"] = annot["transcript_id"].astype(str).map(lambda x: canonical_gene_id(x, gene_prefix))
    annot["Panther_hit"] = annot["Panther_hit"].fillna("").astype(str)
    annot["Panther_hit"] = annot["Panther_hit"].replace(".", "")
    annot = annot.sort_values("Panther_hit", ascending=False)
    annot = annot.drop_duplicates("gene")
    return annot[["gene", "Panther_hit"]]


def load_tf_genes(path: Path, gene_prefix: str) -> set[str]:
    """Load the species transcription-factor gene list in canonical gene-id form."""
    genes = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            gene = line.strip()
            if gene:
                genes.append(canonical_ortholog_id(gene, gene_prefix))
    return set(genes)


def compute_gene_total_umis(raw_dir: Path, cfg: dict, gene_prefix: str) -> pd.DataFrame:
    """Sum total UMIs per canonical gene from the h5ad layer used upstream."""
    h5ad_path = raw_dir / cfg["h5ad"]
    adata = ad.read_h5ad(h5ad_path, backed="r")
    var_names = [canonical_gene_id(v, gene_prefix) for v in adata.var_names.astype(str)]
    adata.file.close()

    with h5py.File(h5ad_path, "r") as h5f:
        ds = h5f["layers"]["total_umis"]
        sums = np.zeros(ds.shape[1], dtype=np.float64)
        chunk = 256
        for start in range(0, ds.shape[0], chunk):
            stop = min(start + chunk, ds.shape[0])
            sums += ds[start:stop, :].sum(axis=0)

    totals = pd.DataFrame({"gene": var_names, "total_umis": sums})
    totals = totals.groupby("gene", as_index=False)["total_umis"].sum()
    return totals


def prepare_orthology(repo_root: Path, raw_dir: Path, config: dict) -> None:
    """Write orthology tables in the format expected by the downstream R stages."""
    orthology_cfg = config["orthology"]
    orthology_path = raw_dir / orthology_cfg["xlsx"]
    orthology = pd.read_excel(orthology_path)
    out_dir = repo_root / config["processed_dir"] / "orthology"
    out_dir.mkdir(parents=True, exist_ok=True)

    sp_cfg = config["species"]
    columns = []
    for species_id in config["pipeline"]["species_order"]:
        col = sp_cfg[species_id]["orthology_column"]
        prefix = sp_cfg[species_id]["gene_prefix"]
        orthology[species_id] = orthology[col].map(lambda x: canonical_ortholog_id(x, prefix) if pd.notna(x) else np.nan)
        columns.append(species_id)

    pair_df = orthology[[orthology_cfg["orthogroup_column"]] + columns].dropna()
    pair_df = pair_df.drop_duplicates()
    pair_df.to_csv(out_dir / "ortholog_pairs.tsv", sep="\t", index=False)

    for species_id in columns:
        mapping = pair_df[[orthology_cfg["orthogroup_column"], species_id]].rename(
            columns={orthology_cfg["orthogroup_column"]: "orthogroup", species_id: "gene"}
        )
        mapping = mapping.drop_duplicates()
        mapping.to_csv(out_dir / f"{species_id}.gene_orthogroups.tsv", sep="\t", index=False)


def write_dataframe(df: pd.DataFrame, path: Path) -> None:
    """Persist a dataframe as TSV, using gzip when the destination suffix requests it."""
    path.parent.mkdir(parents=True, exist_ok=True)
    compression = "gzip" if str(path).endswith(".gz") else None
    df.to_csv(path, sep="\t", index=False, compression=compression)


def process_footprints(
    repo_root: Path,
    raw_dir: Path,
    species_id: str,
    species_cfg: dict,
    cell_metadata: pd.DataFrame,
    tf_genes: set[str],
    annot: pd.DataFrame,
    config: dict,
    chunk_size: int,
) -> None:
    """Build metacell footprints, cell-type footprints, and TF marker tables for one species."""
    out_dir = repo_root / config["processed_dir"] / species_id
    out_dir.mkdir(parents=True, exist_ok=True)

    annot_order = (
        cell_metadata[["leiden_6", "cell_type", "color"]]
        .drop_duplicates()
        .sort_values("leiden_6", key=lambda s: s.astype(int))
        .reset_index(drop=True)
    )
    cell_types = annot_order["cell_type"].tolist()
    metacells = cell_metadata.index.tolist()
    cell_type_masks = {ct: (cell_metadata.loc[metacells, "cell_type"].values == ct) for ct in cell_types}

    raw_fp_path = raw_dir / species_cfg["footprints_csv"]
    gene_prefix = species_cfg["gene_prefix"]
    pseudocount = float(config["pipeline"]["footprint_pseudocount"])
    fp_regularizer = float(config["pipeline"]["footprint_regularizer"])
    logfc_pseudocount = float(config["pipeline"]["marker_logfc_pseudocount"])

    mcs_out_path = out_dir / f"dat.{species_id}.expression.mcs_fp.tsv.gz"
    cts_out_path = out_dir / f"dat.{species_id}.expression.cts_fp.tsv.gz"

    marker_records = {ct: [] for ct in cell_types}
    annot_map = dict(zip(annot["gene"], annot["Panther_hit"]))

    with gzip.open(mcs_out_path, "wt", encoding="utf-8") as mcs_handle, gzip.open(
        cts_out_path, "wt", encoding="utf-8"
    ) as cts_handle:
        mcs_handle.write("\t".join(["gene"] + metacells) + "\n")
        cts_handle.write("\t".join(["gene"] + cell_types) + "\n")

        reader = pd.read_csv(raw_fp_path, index_col=0, chunksize=chunk_size)
        header_checked = False
        for chunk in reader:
            if not header_checked:
                # The raw footprint matrix must stay aligned to the h5ad metacell order,
                # otherwise all downstream cell-type summaries would be mislabelled.
                csv_columns = [str(c) for c in chunk.columns]
                if csv_columns != metacells:
                    raise ValueError(f"Metacell order mismatch in {raw_fp_path}")
                header_checked = True

            # Collapse transcript-style rows onto canonical gene IDs before any averaging.
            chunk.index = [canonical_gene_id(idx, gene_prefix) for idx in chunk.index.astype(str)]
            chunk = chunk.groupby(level=0, sort=False).mean()

            chunk.to_csv(mcs_handle, sep="\t", header=False)

            mat = chunk.to_numpy(dtype=np.float32)
            geomean_cols = []
            for cell_type in cell_types:
                mask = cell_type_masks[cell_type]
                vals = np.clip(mat[:, mask], pseudocount, None)
                geo = np.exp(np.mean(np.log(vals), axis=1))
                geomean_cols.append(geo)
            geomean = np.column_stack(geomean_cols)
            # Regularise and median-normalise each gene so cell-type footprints are
            # interpretable as relative enrichment instead of raw magnitude.
            medians = np.median(fp_regularizer + geomean, axis=1, keepdims=True)
            medians[medians == 0] = 1.0
            cts_fp = (fp_regularizer + geomean) / medians
            cts_chunk = pd.DataFrame(cts_fp, index=chunk.index, columns=cell_types)
            cts_chunk.to_csv(cts_handle, sep="\t", header=False)

            tf_chunk = chunk.loc[chunk.index.intersection(tf_genes)]
            if tf_chunk.empty:
                continue

            # Reuse the same pass through the footprint matrix to assemble TF markers
            # expected by the ancestral-state stage.
            tf_mat = tf_chunk.to_numpy(dtype=np.float32)
            tf_log = np.log2(np.clip(tf_mat, pseudocount, None))
            genes = tf_chunk.index.to_numpy()
            for cell_type in cell_types:
                mask = cell_type_masks[cell_type]
                in_vals = tf_log[:, mask]
                out_vals = tf_log[:, ~mask]
                pvals = stats.ttest_ind(in_vals, out_vals, axis=1, equal_var=False, nan_policy="omit").pvalue
                pvals = np.nan_to_num(pvals, nan=1.0, posinf=1.0, neginf=1.0)
                mean_in = np.mean(tf_mat[:, mask], axis=1)
                mean_out = np.mean(tf_mat[:, ~mask], axis=1)
                avg_log2fc = np.log2((mean_in + logfc_pseudocount) / (mean_out + logfc_pseudocount))
                marker_records[cell_type].append(
                    pd.DataFrame(
                        {
                            "gene": genes,
                            "focus_node": cell_type,
                            "p_val": pvals,
                            "avg_log2FC": avg_log2fc,
                            "is_tf": True,
                            "annotation": [annot_map.get(g, "") for g in genes],
                        }
                    )
                )

    marker_tables = []
    for cell_type, frames in marker_records.items():
        if not frames:
            continue
        df = pd.concat(frames, ignore_index=True)
        df["p_val_adj"] = benjamini_hochberg(df["p_val"].to_numpy())
        df = df.sort_values(["p_val_adj", "avg_log2FC"], ascending=[True, False])
        marker_tables.append(df)

    markers_df = pd.concat(marker_tables, ignore_index=True) if marker_tables else pd.DataFrame(
        columns=["gene", "focus_node", "p_val", "avg_log2FC", "is_tf", "annotation", "p_val_adj"]
    )
    markers_df.to_csv(out_dir / f"seu.{species_id}.markers.cts.tsv.gz", sep="\t", index=False, compression="gzip")


def main() -> None:
    """Parse configuration and materialise all derived annelid inputs."""
    parser = argparse.ArgumentParser(description="Prepare annelid inputs for the cross-species cell-type pipeline.")
    parser.add_argument("--config", default="config/annelid_inputs.json")
    parser.add_argument("--chunk-size", type=int, default=512)
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    config_path = (repo_root / args.config).resolve()
    if not config_path.exists():
        raise FileNotFoundError(
            f"Missing config file: {config_path}. Copy config/annelid_inputs.example.json "
            "to config/annelid_inputs.json or pass --config explicitly."
        )
    with open(config_path, encoding="utf-8") as handle:
        config = json.load(handle)

    raw_dir = Path(config["raw_data_dir"])
    processed_root = repo_root / config["processed_dir"]
    processed_root.mkdir(parents=True, exist_ok=True)

    prepare_orthology(repo_root, raw_dir, config)

    for species_id, species_cfg in config["species"].items():
        cell_metadata, annot = load_species_metadata(raw_dir, species_id, species_cfg)
        species_out = processed_root / species_id
        species_out.mkdir(parents=True, exist_ok=True)

        write_dataframe(annot, species_out / f"annot.{species_id}.leiden.tsv")
        write_dataframe(cell_metadata.reset_index(drop=True), species_out / f"dat.{species_id}.cell_metadata.tsv")

        gene_annot = load_gene_annotations(raw_dir / species_cfg["gene_annotation_tsv"], species_cfg["gene_prefix"])
        write_dataframe(gene_annot, species_out / f"dat.{species_id}.gene_annotations.tsv")

        tf_genes = load_tf_genes(raw_dir / species_cfg["tf_list"], species_cfg["gene_prefix"])
        pd.DataFrame({"gene": sorted(tf_genes)}).to_csv(
            species_out / f"dat.{species_id}.tf_genes.tsv", sep="\t", index=False
        )

        gene_totals = compute_gene_total_umis(raw_dir, species_cfg, species_cfg["gene_prefix"])
        write_dataframe(gene_totals, species_out / f"dat.{species_id}.gene_total_umis.tsv")

        process_footprints(
            repo_root=repo_root,
            raw_dir=raw_dir,
            species_id=species_id,
            species_cfg=species_cfg,
            cell_metadata=cell_metadata,
            tf_genes=tf_genes,
            annot=gene_annot,
            config=config,
            chunk_size=args.chunk_size,
        )


if __name__ == "__main__":
    main()
