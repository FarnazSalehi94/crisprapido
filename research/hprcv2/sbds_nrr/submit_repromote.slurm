#!/bin/bash
#SBATCH --job-name=sbds-repromote
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
# One promotion+binding task per promoted-set candidate. Runs AFTER cross-node
# reverify. Promotes (quarantine+rollback), then binds at the FINAL
# authoritative path (fresh binding if none exists; the search campaign picks
# the assembly up in the next driver wave).
set -euo pipefail
: "${ASSEMBLIES_FILE:?}"
mapfile -t NAMES < "$ASSEMBLIES_FILE"
NAME=${NAMES[$SLURM_ARRAY_TASK_ID]}
ROOT=/moosefs/salehi/HPRCv2_CRISPRapido/runs/sbds_ngg_nrr_v1
GZ=$(awk -F'\t' -v a="$NAME" '$1==a {print $2}' /moosefs/salehi/HPRCv2_CRISPRapido/analysis_genomes.tsv)
# Cross-node trust gate: ensure this worker differs from the BUILD worker.
BUILD_HOST=$(cut -d= -f2 "/moosefs/salehi/HPRCv2_CRISPRapido/indices_rebuilt_v1/$NAME/build-host.txt")
if [ "$(hostname)" = "$BUILD_HOST" ]; then echo "promotion must run on a different node than build node $BUILD_HOST" >&2; exit 1; fi
python3 "$ROOT/00_manifest/rebuild-tooling-v1/reverify_sbds.py" "$NAME" "$BUILD_HOST"
python3 "$ROOT/00_manifest/promotion-tooling-v1/promote_index.py" "$NAME"
if [ ! -f "$ROOT/00_manifest/bindings/$NAME/binding.json" ]; then
  python3 "$ROOT/00_manifest/prepare_binding-v1.py" --assembly "$NAME" \
    --index-prefix "/moosefs/salehi/HPRCv2_CRISPRapido/indices/$NAME/$NAME" \
    --fasta-gz "$GZ" --fai "$GZ.fai" \
    --out "$ROOT/00_manifest/bindings/$NAME" \
    --scoring-python /usr/bin/python3 \
    --reference-index-helper /export/local/home/salehi/crisprapido-pam-validation/research/crispr_bulge/reference_index.py
else
  python3 "$ROOT/00_manifest/promotion-tooling-v1/repoint_binding.py" "$NAME"
fi
