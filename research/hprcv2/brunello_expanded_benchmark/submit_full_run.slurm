#!/usr/bin/env bash
# PREPARED TEMPLATE — do not submit before user authorization of the full run.
# One array task per assembly; the worker skips assemblies that are already
# COMPLETE with matching plan checksums, so re-running a partial array is safe.
#SBATCH --job-name=brunello-full
#SBATCH --partition=workers
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=72G
#SBATCH --time=72:00:00
#SBATCH --no-requeue
#SBATCH --array=0-463%32
set -euo pipefail
umask 077
: "${SLURM_JOB_ID:?}" "${SLURMD_NODENAME:?}" "${SLURM_ARRAY_TASK_ID:?}"
plan=${1:?absolute authorized full-run plan required}
manifest=${2:?464-assembly manifest required}
[[ "$plan" = /* && "$manifest" = /* ]]
export USER="$(id -un)" PYTHONDONTWRITEBYTECODE=1
assembly=$(awk -v id="$SLURM_ARRAY_TASK_ID" 'NR == id+1 {print $1}' "$manifest")
[[ -n "$assembly" ]]
worker=$(/usr/bin/python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["worker"])' "$plan")
[[ "$worker" = /* ]]
# One application thread; the extra allocation exists only so this single-CPU
# job is not packed against memory-heavy neighbours on CR_CPU partitions.
exec /usr/bin/python3 "$worker" --plan "$plan" --assembly "$assembly"
