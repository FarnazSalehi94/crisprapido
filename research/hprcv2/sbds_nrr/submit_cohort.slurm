#!/bin/bash
#SBATCH --job-name=sbds-nrr-cohort
#SBATCH --cpus-per-task=1
#SBATCH --mem=56G
#SBATCH --time=12:00:00
# Restartable per-assembly array. Completed assemblies are checksum-verified
# and skipped; failed assemblies keep FAILED.json plus verified search
# checkpoints for resume. Re-submit the same command to resume; the worker
# rejects plan drift. Do not pass comma lists through --export (sbatch splits
# them); use ASSEMBLIES_FILE with one assembly name per line.
#
# sbatch --array=0-465%16 --export=ALL,PLAN=...,ASSEMBLIES_FILE=... submit_cohort.slurm
set -euo pipefail
: "${PLAN:?}" "${ASSEMBLIES_FILE:?}"
mapfile -t NAMES < "$ASSEMBLIES_FILE"
NAME=${NAMES[$SLURM_ARRAY_TASK_ID]}
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export TF_NUM_INTRAOP_THREADS=1 TF_NUM_INTEROP_THREADS=1 CUDA_VISIBLE_DEVICES=-1
ROOT=/moosefs/salehi/HPRCv2_CRISPRapido/runs/sbds_ngg_nrr_v1
mkdir -p "$ROOT/logs"
python3 "$ROOT/00_manifest/nrr-tooling-v1/supervise.py" --limit-kib 50331648 \
  --interval .1 --log "$ROOT/logs/cohort-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.jsonl" \
  --summary "$ROOT/logs/cohort-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.json" \
  -- python3 "$ROOT/00_manifest/nrr-tooling-v1/cohort_worker.py" --plan "$PLAN" --assembly "$NAME"
