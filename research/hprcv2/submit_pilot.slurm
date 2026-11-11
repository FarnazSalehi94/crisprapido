#!/usr/bin/env bash
#SBATCH --job-name=hprcv2-guarded-pilot
#SBATCH --partition=1tbmem
#SBATCH --nodelist=octopus02
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --no-requeue
set -euo pipefail
umask 077
: "${SLURM_JOB_ID:?}" "${SLURMD_NODENAME:?}"
plan=${1:?absolute immutable pilot-plan.json required}
task=${2:?explicit pilot configuration required}
[[ "$plan" = /* ]]
export USER="$(id -un)" PYTHONDONTWRITEBYTECODE=1
worker=$(/usr/bin/python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["worker"])' "$plan")
[[ "$worker" = /* ]]
exec /usr/bin/python3 "$worker" --plan "$plan" --task "$task"
