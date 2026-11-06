#!/usr/bin/env bash
#SBATCH --job-name=crispr-hg00735-octopus10
#SBATCH --partition=workers
#SBATCH --nodelist=octopus10
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --no-requeue
#SBATCH --output=crispr-pilot-octopus10-%j.slurm.log
# Prepared for review only. The delegated pilot retains all binding, cgroup,
# scratch, single-guide, incomplete-output and independent-audit checks.
set -euo pipefail
: "${SLURM_JOB_ID:?requires Slurm}" "${SLURMD_NODENAME:?requires Slurm node name}" "${REPO:?source pilot.env}"
[[ "$SLURMD_NODENAME" == octopus10 ]] || { echo 'Wrong node; expected octopus10' >&2; exit 1; }
[[ "$REPO" = /* && -r "$REPO/research/pam_anchored_poc/pilot_preparation/one_guide.slurm" ]] || { echo 'Pilot entry point unavailable' >&2; exit 1; }
exec bash "$REPO/research/pam_anchored_poc/pilot_preparation/one_guide.slurm"
