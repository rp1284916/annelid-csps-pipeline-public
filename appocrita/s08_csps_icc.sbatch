#!/bin/bash
# Submit the ICC ortholog-scoring stage with enough memory for large
# coexpression matrices.
#SBATCH --job-name=annelid_s08
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=24:00:00
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
bash "$repo_dir/appocrita/run_stage.sh" "s08_csps_icc_calculation_2026-03-25.R"
