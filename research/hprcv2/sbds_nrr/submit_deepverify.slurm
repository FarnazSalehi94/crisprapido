#!/bin/bash
#SBATCH --job-name=sbds-deepverify
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
# Stage-1 deep verification batch. One assembly per task; restartable (committed
# assemblies are rechecked and skipped). ASSEMBLIES_FILE: one name per line.
# Never rebuilds or overwrites originals; failures preserve evidence.
set -euo pipefail
: "${ASSEMBLIES_FILE:?}"
mapfile -t NAMES < "$ASSEMBLIES_FILE"
NAME=${NAMES[$SLURM_ARRAY_TASK_ID]}
ROOT=/moosefs/salehi/HPRCv2_CRISPRapido/runs/sbds_ngg_nrr_v1
mkdir -p "$ROOT/logs"
python3 "$ROOT/00_manifest/deep-verify-v1/deep_verify.py" "$NAME"
