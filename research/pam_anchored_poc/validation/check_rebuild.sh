#!/usr/bin/env bash
set -euo pipefail
: "${COLUMBA_BUILD_BIN:?set COLUMBA_BUILD_BIN}"
ref=${1:?reference FASTA}; old=${2:?existing index prefix}; new=${3:?scratch index prefix}
"$COLUMBA_BUILD_BIN" -r "$new" -f "$ref" > "$new.build.log" 2>&1
for ext in bwt brt rev.brt sa.4 sa.bv.4 pos sna txt.bin cct meta; do
 if cmp -s "$old.$ext" "$new.$ext"; then printf '%s identical\n' "$ext"; else printf '%s DIFFERENT\n' "$ext"; exit 1; fi
done
