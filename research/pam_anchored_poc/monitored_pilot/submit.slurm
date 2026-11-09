#!/usr/bin/env bash
#SBATCH --job-name=crispr-hg00735-monitored
#SBATCH --partition=workers
#SBATCH --nodelist=octopus10
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --no-requeue
#SBATCH --output=crispr-monitored-%j.slurm.log
set -euo pipefail
: "${SLURM_JOB_ID:?}" "${SLURMD_NODENAME:?}" "${REPO:?}" "${RESULT_ROOT:?}"
[[ "$SLURMD_NODENAME" == octopus10 ]]
monitor_dir="$RESULT_ROOT/monitor-$SLURM_JOB_ID"
mkdir "$monitor_dir"
export MONITORED_PILOT_SUPERVISED="$SLURM_JOB_ID"
cd "$REPO"
exec /usr/bin/python3 research/pam_anchored_poc/monitored_pilot/supervise.py --limit-kib 50331648 --interval .1 --log "$monitor_dir/rss.log" --summary "$monitor_dir/summary.json" -- bash research/pam_anchored_poc/monitored_pilot/workload.sh
