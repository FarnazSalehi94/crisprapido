#!/bin/bash
#SBATCH --job-name=sbds-nrr-pilot
#SBATCH --cpus-per-task=1
#SBATCH --mem=56G
#SBATCH --time=04:00:00
#SBATCH --array=0-3%1
set -euo pipefail
: "${PLAN:?export absolute frozen PLAN}"
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export TF_NUM_INTRAOP_THREADS=1 TF_NUM_INTEROP_THREADS=1 CUDA_VISIBLE_DEVICES=-1
mapfile -t paths < <(python3 - "$PLAN" <<'PY'
import json,sys
p=json.load(open(sys.argv[1]))
for k in ['supervisor','pilot_worker','root']: print(p[k])
PY
)
pams=(GG AG GA AA)
python3 "${paths[0]}" --limit-kib 50331648 --interval .1 \
 --log "${paths[2]}/logs/monitor-${SLURM_JOB_ID}.jsonl" \
 --summary "${paths[2]}/logs/monitor-${SLURM_JOB_ID}.json" \
 -- python3 "${paths[1]}" --plan "$PLAN" --pam "${pams[$SLURM_ARRAY_TASK_ID]}"
