#!/bin/bash
#SBATCH -p workers
#SBATCH -w octopus11
#SBATCH -c 48
#SBATCH -J crisprapido_hg00097

set -euo pipefail

REFERENCE="/lizardfs/old/salehi/crispr2/HG00097/HG00097.fa"
GUIDES="/lizardfs/old/salehi/crispr2/crispr-progress/ref-only/geckov2_gRNAs.txt"
OUTPUT="/lizardfs/old/salehi/crispr2/crispr-progress/crisprapido-pangenome/HG00097/HG00097_results.txt"
THREADS=48
THREAD_MULTIPLIER=3.0
CHANNEL_DEPTH=200000

TMPDIR="/tmp/crisprapido_${SLURM_JOB_ID:-$$}"
mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT

LOCAL_REFERENCE="$TMPDIR/reference.fa"
LOCAL_GUIDES="$TMPDIR/guides.txt"
LOCAL_OUTPUT="$TMPDIR/results.txt"

cp "$REFERENCE" "$LOCAL_REFERENCE"
cp "$GUIDES" "$LOCAL_GUIDES"

crisprapido \
  -r "$LOCAL_REFERENCE" \
  --guides-file "$LOCAL_GUIDES" \
  -p NGG \
  -m 4 \
  --threads "$THREADS" \
  --thread-multiplier "$THREAD_MULTIPLIER" \
  --allow-oversubscribe \
  --channel-depth "$CHANNEL_DEPTH" \
  > "$LOCAL_OUTPUT"

cp "$LOCAL_OUTPUT" "$OUTPUT"
