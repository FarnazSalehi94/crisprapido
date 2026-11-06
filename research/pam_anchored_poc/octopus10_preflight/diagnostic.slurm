#!/usr/bin/env bash
#SBATCH --job-name=crispr-octopus10-diagnostic
#SBATCH --partition=workers
#SBATCH --nodelist=octopus10
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:10:00
#SBATCH --no-requeue
#SBATCH --output=crispr-diagnostic-octopus10-%j.slurm.log
# Read-only checks only; do not submit without explicit user approval.
set -euo pipefail
: "${SLURM_JOB_ID:?requires Slurm}" "${SLURMD_NODENAME:?requires Slurm node name}" "${REPO:?source pilot.env}"
[[ "$SLURMD_NODENAME" == octopus10 ]] || { echo 'Wrong node; expected octopus10' >&2; exit 1; }
[[ "$REPO" = /* && -r "$REPO/research/pam_anchored_poc/pilot_preparation/worker_diagnostic.slurm" ]] || { echo 'Diagnostic entry point unavailable' >&2; exit 1; }
exec bash "$REPO/research/pam_anchored_poc/pilot_preparation/worker_diagnostic.slurm"
