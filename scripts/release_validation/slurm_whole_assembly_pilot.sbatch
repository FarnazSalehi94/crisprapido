#!/bin/bash
#SBATCH --job-name=crispr_pilot_search
#SBATCH --partition=1tbmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --nice=10000
set -euo pipefail
: "${PILOT_ROOT:?}" "${PILOT_ASSEMBLY:?}" "${PILOT_GUIDE:?}" "${PILOT_REFERENCE:?}" "${PILOT_INDEX:?}"
exec nice -n 19 ionice -c 3 python3 "$PILOT_ROOT/tools/whole_assembly_pilot.py" \
 --root "$PILOT_ROOT" --assembly "$PILOT_ASSEMBLY" --guide "$PILOT_GUIDE" \
 --reference "$PILOT_REFERENCE" --index "$PILOT_INDEX" \
 --binary /export/local/home/salehi/crisprapido-release-pilot-0aa9944-20261102/release/crisprapido \
 --profile /export/local/home/salehi/crisprapido-profile-fast-20261031-target/release/crisprapido \
 --columba /moosefs/salehi/HPRCv2_CRISPRapido/results/PAM_CFD_complete_search_v1_9c2d3f5_2026-10-28/release/columba-runtime \
 --matrices /export/local/home/salehi/crisprapido-release-pilot-0aa9944-20261102/release \
 --samtools /export/local/home/salehi/verkko_env/bin/samtools \
 --time-bin /home/salehi/miniconda3/bin/time \
 --oracle "$PILOT_ROOT/inputs/${PILOT_ASSEMBLY}_${PILOT_GUIDE}.jsonl"
