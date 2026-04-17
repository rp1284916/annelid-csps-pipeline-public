# Annelid CSPS Pipeline

This repository is a code-only, portfolio-style snapshot of an annelid cross-species cell-state analysis workflow. It is shared as a methods and engineering example rather than as a data release.

No raw data, derived analysis outputs, figures, or workstation-specific configs are included here.

## Scope

This repo contains:

- pipeline scripts and helper code
- Appocrita batch runners
- a sanitized config template
- documentation for local setup and run order

This repo does not contain:

- raw input data
- processed matrices or intermediate objects
- result tables, figures, or trees generated from lab data
- absolute local paths or private machine-specific configs
- project notes that mix code with unpublished interpretation

## Layout

- `scripts/prepare_annelid_inputs.py`
  prepares processed inputs from a local raw-data package
- `R/csps_minimal.R`
  minimal helper layer used by the R stages
- `results_scatlas/s08_csps_icc_calculation_2026-03-25.R`
  computes ICC-based ortholog scoring
- `results_scatlas/s09_csps_similarity_2026-03-25.R`
  computes cross-species similarity summaries
- `results_scatlas/s10_cell_type_trees_2026-03-25.R`
  builds reduced-dimension summaries, trees, and node-level marker outputs
- `scripts/render_tf_program_tree.py`
  renders a tree with editable annotation columns
- `scripts/render_stage_expression_comparison.py`
  plots species-specific developmental trajectories for a chosen orthogroup or gene pair
- `appocrita/`
  Slurm submission helpers for longer runs

## Local Setup

1. Copy `config/annelid_inputs.example.json` to `config/annelid_inputs.json`.
2. Replace the placeholder paths with local paths on your machine.
3. Keep `config/annelid_inputs.json` untracked.

The tracked template uses fake placeholder filenames on purpose. Adjust it to match your local raw-data package layout.

## Run Order

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

Those files are local analysis artifacts and are not part of this public repo.

## Synthetic Example

The tree renderer includes a built-in synthetic example so the plotting code can be demonstrated without lab data:

```bash
py -3 scripts/render_tf_program_tree.py --example --output-prefix work/results_cell_type_trees/example_tf_program_tree
```

## Appocrita

The `appocrita/` directory contains Slurm submission helpers for longer runs on Appocrita.

Typical usage:

```bash
sbatch appocrita/install_r_deps.sbatch
sbatch appocrita/s08_csps_icc.sbatch
sbatch appocrita/s09_csps_similarity.sbatch
sbatch appocrita/s10_cell_type_trees.sbatch
```

## Public Repo Hygiene

The following are intentionally excluded from version control in this public copy:

- `config/annelid_inputs.json`
- `work/`
- `annelid-csps-results/`
- common single-cell data formats such as `.h5ad` and `.rds`

Before publishing or sharing further, review any new files for data, results, screenshots, or project-specific interpretation.
