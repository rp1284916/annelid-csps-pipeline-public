#!/bin/bash
# Submit the dimensional-reduction, tree-building, and ancestral-marker stage.
#SBATCH --job-name=annelid_s10
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=16:00:00
#SBATCH --output=%x.o%j

set -euo pipefail

repo_dir="${CSPS_REPO_DIR:-${SLURM_SUBMIT_DIR:-$HOME/annelid-csps-pipeline}}"
module load R
export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/%p-library/%v}"
# Keep BLAS/OpenMP thread counts aligned with Slurm to avoid oversubscription.
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export CSPS_THREADS="${CSPS_THREADS:-${SLURM_CPUS_PER_TASK:-1}}"
export CSPS_STAGE_TO_SCRATCH="${CSPS_STAGE_TO_SCRATCH:-1}"
bash "$repo_dir/appocrita/run_stage.sh" "s10_cell_type_trees_2026-03-25.R"
