#!/usr/bin/env bash
#SBATCH --job-name=brunello-bench16
#SBATCH --partition=1tbmem
#SBATCH --nodelist=octopus02
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --no-requeue
set -euo pipefail
umask 077
: "${SLURM_JOB_ID:?}" "${SLURMD_NODENAME:?}"
plan=${1:?absolute immutable expanded-benchmark-plan.json required}
assembly=${2:?assembly required}
config=${3:?configuration required}
[[ "$plan" = /* ]]
export USER="$(id -un)" PYTHONDONTWRITEBYTECODE=1
worker=$(/usr/bin/python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["worker"])' "$plan")
[[ "$worker" = /* ]]
exec /usr/bin/python3 "$worker" --plan "$plan" --assembly "$assembly" --config "$config"
