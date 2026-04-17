#!/bin/bash
# Install the R package set once into the user's home-library on Appocrita.
#SBATCH --job-name=annelid_rdeps
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=%x.o%j

set -euo pipefail

repo_dir="${CSPS_REPO_DIR:-${SLURM_SUBMIT_DIR:-$HOME/annelid-csps-pipeline}}"
module load R
# Appocrita's module library is read-only, so install into the user library.
export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/%p-library/%v}"
mkdir -p "$R_LIBS_USER"
cd "$repo_dir"
Rscript scripts/install_r_deps.R
