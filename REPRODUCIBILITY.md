# Reproducibility and Availability Crosswalk

## Purpose

This document clarifies the scope of the public release associated with the submitted annelid dissertation. The complete research workspace includes laboratory-supplied inputs, derived matrices, curation tables, result tables, figures, and analysis-specific notes. Those materials are not automatically suitable for public redistribution.

The public repository therefore prioritises a clear disclosure boundary over literal publication of every file named in the dissertation. It contains sanitised code and templates only. Omission from this repository does not mean that a file was absent from the research workspace.

## Dissertation-to-repository crosswalk

| Dissertation statement or named material | Public status | Reason or public substitute |
| --- | --- | --- |
| Adapted ICC stages s08-s10 | Included | The R stage scripts and reduced helper library are public under `results_scatlas/` and `R/`. |
| Python input preparation | Included | A sanitised, configurable preparation script is public as `scripts/prepare_annelid_inputs.py`. |
| Appocrita execution | Included | Generic Slurm wrappers are public under `appocrita/`; logs, account details, and cluster-specific paths are excluded. |
| Configuration files | Partly included | `config/annelid_inputs.example.json` documents the expected schema with placeholder paths. Real configurations are excluded. |
| `broad_celltype_mapping.tsv` | Not public | This is a project curation table derived from laboratory annotations. Broad categories are reported in Supplementary Table S1 of the dissertation, but the machine-readable working table is withheld. |
| `reciprocal_best_hits_<run>.tsv` | Not public | These files are derived result tables. The reported hits and calibration summaries remain in the dissertation; the machine-readable working outputs are withheld. |
| `monoaminergic_candidate_clusters.tsv` | Not public | This is a derived candidate-selection result. The retained candidates are reported in the dissertation, but the working table is withheld. |
| `monoaminergic_tf_validation.tsv` | Not public | This combines derived measurements and curated evidence calls. The reported validation summary remains in the dissertation, but the working table is withheld. |
| Input footprint matrices, atlas objects, and orthogroup tables | Not public | These were supplied by, or derived from material supplied by, the Martin-Duran laboratory and are not redistributed. |
| SAMap scripts and mapping resources | Not public in this snapshot | The research version depends on local atlas objects, proteomes, reciprocal mapping files, and workspace-specific assumptions. It requires a separate sanitisation and laboratory approval review before publication. |
| Focused, stage-aware, validation, calibration, added-species, and complete figure-generation workflows | Not public in this snapshot | These scripts are interleaved with derived inputs, project-specific labels, and unpublished working material. Only independently sanitised utilities are included here. |
| Exact software environment | Not fully recoverable here | The versions explicitly reported in the dissertation are recorded below. A complete lockfile was not deposited with the public snapshot, so the repository does not claim byte-for-byte environment reconstruction. |

## Software versions reported in the dissertation

The methods section reports the following versions:

- Python 3.12
- Scanpy 1.12.1
- NumPy 2.4.6
- Numba 0.65.1
- SAMap (`sc-samap`) 3.0.1

AnnData, SciPy, pandas, and Matplotlib are also named, but their exact versions are not stated in the submitted text. R, `ape`, and the other R dependencies are named without a complete version-locked environment. The dependency installer in `scripts/install_r_deps.R` records the required R package set but deliberately does not pretend to reconstruct historical package versions.

## What can be reproduced publicly

The repository documents the sanitised pipeline structure, data contracts, main ICC calculations, similarity stage, tree/ancestral-marker stage, generic batch execution, and selected rendering utilities. The synthetic tree-rendering example can be run without laboratory data.

Reproducing the dissertation's numerical results requires the non-public inputs and derived intermediate files. The public snapshot alone is therefore a methods/code release, not a complete computational reproducibility package.

## Access and future additions

Questions about access to omitted material should be directed to the data-owning laboratory and handled under its approval and data-sharing procedures. Future additions to this repository should be limited to files that have been independently reviewed for:

1. laboratory-supplied or unpublished data;
2. derived result values and biological interpretations;
3. sample, cell, gene, barcode, or internal identifier leakage;
4. local paths, usernames, email addresses, cluster accounts, and logs;
5. credentials, tokens, secrets, and private service endpoints;
6. licensing or redistribution restrictions on third-party datasets.

This clarification supersedes broader repository-availability wording in the submitted dissertation; it does not alter the scientific results reported there.
