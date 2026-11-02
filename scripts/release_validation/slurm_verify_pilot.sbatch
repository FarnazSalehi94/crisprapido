#!/bin/bash
#SBATCH --job-name=crispr_pilot_verify
#SBATCH --partition=1tbmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=04:00:00
#SBATCH --nice=10000
set -euo pipefail
: "${PILOT_ROOT:?must provide PILOT_ROOT}"
exec nice -n 19 ionice -c 3 python3 "$PILOT_ROOT/tools/verify_pilot_manifest.py" \
  --assemblies "$PILOT_ROOT/inputs/verified_assemblies.tsv" \
  --checksums "$PILOT_ROOT/inputs/verified_component_checksums.tsv" \
  --binary /export/local/home/salehi/crisprapido-release-pilot-0aa9944-20261102/release/crisprapido \
  --freeze "$PILOT_ROOT/inputs/freeze_manifest.json" \
  --rate-mib 128 --expected-assemblies 12 --output "$PILOT_ROOT/verification"
