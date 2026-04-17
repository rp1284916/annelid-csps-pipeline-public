#!/bin/bash
# Submit s08 first, then fan out the two stages that only depend on ICC output.
set -euo pipefail

repo_dir="$(cd "$(dirname "$0")/.." && pwd)"
cd "$repo_dir"

job_s08="$(sbatch --parsable appocrita/s08_csps_icc.sbatch)"
job_s09="$(sbatch --parsable --dependency=afterok:${job_s08} appocrita/s09_csps_similarity.sbatch)"
job_s10="$(sbatch --parsable --dependency=afterok:${job_s08} appocrita/s10_cell_type_trees.sbatch)"

printf 'Submitted s08: %s\n' "$job_s08"
printf 'Submitted s09: %s (afterok:%s)\n' "$job_s09" "$job_s08"
printf 'Submitted s10: %s (afterok:%s)\n' "$job_s10" "$job_s08"
