# Annelid CSPS Pipeline

This repository is a sanitised, code-only snapshot of an annelid cross-species cell-state analysis workflow. It is shared as a methods and software portfolio piece, not as a data release.

No raw data, processed matrices, result tables, figures, unpublished interpretations, or workstation-specific configs are included here.

## What This Pipeline Does

At a high level, the workflow:

1. prepares species-specific input tables from local raw-data packages
2. scores ortholog conservation with ICC-based filtering
3. computes cross-species cell-state similarity
4. builds reduced-dimension summaries and tree-oriented downstream outputs
5. renders publication-style visualizations from locally generated result files

The repository is designed to show the code structure and analysis logic behind that workflow without exposing lab data or unpublished results.

## Public Repo Scope

Included here:

- pipeline scripts and helper code
- Appocrita batch runners
- a sanitized config template
- documentation for local setup and execution order

Not included here:

- raw input data
- processed data products or intermediate objects
- result tables, figures, dendrograms, or derived outputs from lab data
- local absolute paths or real machine-specific config files
- internal notes that mix code with unpublished interpretation

## Repository Layout

- `scripts/prepare_annelid_inputs.py`
  prepares processed inputs from a local raw-data package
- `R/csps_minimal.R`
  minimal helper layer used by the downstream R stages
- `results_scatlas/s08_csps_icc_calculation_2026-03-25.R`
  computes ICC-based ortholog scoring
- `results_scatlas/s09_csps_similarity_2026-03-25.R`
  computes cross-species similarity summaries
- `results_scatlas/s10_cell_type_trees_2026-03-25.R`
  builds reduced-dimension summaries, tree-oriented outputs, and node-level marker tables
- `scripts/render_tf_program_tree.py`
  renders a tree with editable annotation columns
- `scripts/render_stage_expression_comparison.py`
  plots species-specific developmental trajectories for a chosen orthogroup or gene pair
- `appocrita/`
  Slurm submission helpers for longer runs

## Requirements

The codebase uses both R and Python.

Typical local requirements:

- R plus the packages installed by `scripts/install_r_deps.R`
- Python 3
- Python packages used by individual scripts, including `matplotlib`, `biopython`, `anndata`, `h5py`, and `pandas` where relevant

Some scripts are lightweight utilities, while others expect a local single-cell analysis environment and a local raw-data package that is not distributed in this repository.

## Local Setup

1. Copy `config/annelid_inputs.example.json` to `config/annelid_inputs.json`.
2. Replace the placeholder paths with local paths on your machine.
3. Keep `config/annelid_inputs.json` untracked.

The tracked template uses fake placeholder filenames on purpose. Adjust it to match your local raw-data package layout.

## Minimal Run Order

From the repo root:

```bash
Rscript scripts/install_r_deps.R
py -3 scripts/prepare_annelid_inputs.py
cd results_scatlas
Rscript s08_csps_icc_calculation_2026-03-25.R
Rscript s09_csps_similarity_2026-03-25.R
Rscript s10_cell_type_trees_2026-03-25.R
```

If you keep multiple local configs, pass one explicitly:

```bash
py -3 scripts/prepare_annelid_inputs.py --config path/to/config.json
```

## Outputs

Local runs write under `work/`, which is intentionally git-ignored. Typical output folders are:

- `work/results/`
- `work/results_cell_type_trees/`
- `work/results_stage_expression/`

Those directories are local analysis artifacts and are not part of this public repository.

## Synthetic Example

The tree renderer includes a built-in synthetic example so the plotting code can be demonstrated without lab data:

```bash
py -3 scripts/render_tf_program_tree.py --example --output-prefix work/results_cell_type_trees/example_tf_program_tree
```

## Running on Appocrita

The `appocrita/` directory contains Slurm submission helpers for longer runs.

Typical usage:

```bash
sbatch appocrita/install_r_deps.sbatch
sbatch appocrita/s08_csps_icc.sbatch
sbatch appocrita/s09_csps_similarity.sbatch
sbatch appocrita/s10_cell_type_trees.sbatch
```

## Git Hygiene for This Public Copy

The following are intentionally excluded from version control:

- `config/annelid_inputs.json`
- `work/`
- `annelid-csps-results/`
- common single-cell data formats such as `.h5ad` and `.rds`

Before sharing any future changes, review new files for data, results, screenshots, or project-specific interpretation.
