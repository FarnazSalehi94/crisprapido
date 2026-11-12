#!/bin/bash
#SBATCH --job-name=sbds-rebuild-pilot
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=04:00:00
set -euo pipefail
: "${ASSEMBLIES_FILE:?}"
mapfile -t NAMES < "$ASSEMBLIES_FILE"
NAME=${NAMES[$SLURM_ARRAY_TASK_ID]}
ROOT=/moosefs/salehi/HPRCv2_CRISPRapido/runs/sbds_ngg_nrr_v1
python3 "$ROOT/00_manifest/deep-verify-v1/rebuild_pilot.py" --assembly "$NAME" --out "/export/local/home/salehi/sbds_rebuilt_indexes_v1/$NAME"
