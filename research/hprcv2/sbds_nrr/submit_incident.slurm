#!/bin/bash
#SBATCH --job-name=sbds-incident
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=06:00:00
set -euo pipefail
: "${ASSEMBLIES_FILE:?}"
mapfile -t NAMES < "$ASSEMBLIES_FILE"
ROOT=/moosefs/salehi/HPRCv2_CRISPRapido/runs/sbds_ngg_nrr_v1
python3 "$ROOT/00_manifest/deep-verify-v1/incident_scan.py" "${NAMES[$SLURM_ARRAY_TASK_ID]}"
