#!/bin/bash
# Submit the weighted cell-type similarity stage.
#SBATCH --job-name=annelid_s09
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=08:00:00
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
bash "$repo_dir/appocrita/run_stage.sh" "s09_csps_similarity_2026-03-25.R"
