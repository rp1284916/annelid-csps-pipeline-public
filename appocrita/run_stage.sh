#!/bin/bash
# Run one R stage either directly from the repo or from node-local scratch,
# then always sync result files back to the source checkout.
set -euo pipefail

stage_script="${1:?usage: run_stage.sh <stage-script>}"
repo_dir="${CSPS_REPO_DIR:-${SLURM_SUBMIT_DIR:-$HOME/annelid-csps-pipeline}}"
stage_to_tmpdir="${CSPS_STAGE_TO_SCRATCH:-${CSPS_STAGE_TO_TMPDIR:-1}}"
stage_root="${CSPS_STAGE_ROOT:-}"
job_id="${SLURM_JOB_ID:-${JOB_ID:-manual}}"
run_repo="$repo_dir"

sync_back() {
	local status=$?
	trap - EXIT
	if [[ "$stage_to_tmpdir" == "1" && "$run_repo" != "$repo_dir" ]]; then
		# Sync only materialised results back; the scratch copy is disposable.
		mkdir -p "$repo_dir/work/results" "$repo_dir/work/results_cell_type_trees"
		if [[ -d "$run_repo/work/results" ]]; then
			rsync -a "$run_repo/work/results/" "$repo_dir/work/results/"
		fi
		if [[ -d "$run_repo/work/results_cell_type_trees" ]]; then
			rsync -a "$run_repo/work/results_cell_type_trees/" "$repo_dir/work/results_cell_type_trees/"
		fi
	fi
	exit "$status"
}

trap sync_back EXIT

if [[ "$stage_to_tmpdir" == "1" ]]; then
	if [[ -z "$stage_root" ]]; then
		if [[ -n "${SLURM_TMPDIR:-}" ]]; then
			stage_root="$SLURM_TMPDIR"
		elif [[ -n "${TMPDIR:-}" ]]; then
			stage_root="$TMPDIR"
		elif [[ -d "/gpfs/scratch" ]]; then
			stage_root="/gpfs/scratch/$USER"
		else
			stage_root="/tmp"
		fi
	fi
	job_root="${stage_root%/}/csps-${job_id}"
	run_repo="$job_root/repo"
	rm -rf "$run_repo"
	mkdir -p "$run_repo"
	# Stage a clean working copy onto fast local storage to avoid hammering the
	# shared filesystem during large R jobs.
	rsync -a --delete \
		--exclude '.git/' \
		--exclude 'appocrita/logs/' \
		"$repo_dir/" "$run_repo/"
fi

cd "$run_repo/results_scatlas"
Rscript "$stage_script"
