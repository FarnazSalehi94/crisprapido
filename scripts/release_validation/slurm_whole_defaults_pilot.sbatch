#!/bin/bash
#SBATCH --job-name=crispr_defaults_pilot
#SBATCH --partition=1tbmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --nice=10000
set -euo pipefail
: "${PILOT_ROOT:?}" "${PILOT_ASSEMBLY:?}" "${PILOT_GUIDE:?}" "${PILOT_REFERENCE:?}" "${PILOT_INDEX:?}"
# A single output file may not exceed 2 GiB; abort instead of unbounded growth.
ulimit -f 2097152
exec nice -n 19 ionice -c 3 python3 "$PILOT_ROOT/tools/whole_defaults_pilot.py" \
 --root "$PILOT_ROOT" --assembly "$PILOT_ASSEMBLY" --guide "$PILOT_GUIDE" \
 --reference "$PILOT_REFERENCE" --index "$PILOT_INDEX" \
 --old-root /moosefs/salehi/HPRCv2_CRISPRapido/results/whole_assembly_resource_pilot_0aa9944_20261102_v1 \
 --binary "$PILOT_ROOT/tools/crisprapido" \
 --columba /moosefs/salehi/HPRCv2_CRISPRapido/results/PAM_CFD_complete_search_v1_9c2d3f5_2026-10-28/release/columba-runtime \
 --matrices /export/local/home/salehi/crisprapido-release-pilot-0aa9944-20261102/release \
 --samtools /export/local/home/salehi/verkko_env/bin/samtools --time-bin /home/salehi/miniconda3/bin/time \
 --oracle "$PILOT_ROOT/inputs/${PILOT_ASSEMBLY}_${PILOT_GUIDE}.jsonl"
