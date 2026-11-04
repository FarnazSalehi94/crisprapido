#!/usr/bin/env bash
set -euo pipefail
: "${COLUMBA_BUILD_BIN:?set COLUMBA_BUILD_BIN to a Columba index-builder wrapper}"
: "${POC_BIN:?set POC_BIN to a prototype executable wrapper}"
root=$(cd "$(dirname "$0")/../.." && pwd)
out=${1:?provide scratch directory}
mkdir -p "$out/alt"
for name in indexed_synthetic indexed_real_chm13 indexed_12_windows; do
 p="$root/docs/coverage_oracles_2026-10-29_v1/$name"
 "$COLUMBA_BUILD_BIN" -r "$out/$name" -f "$p/ref.fa" >"$out/$name.build.log" 2>&1
 "$POC_BIN" "$out/$name" "$p/ref.fa" "$p/guides.fa" 4 1 0.75 >"$out/$name.tsv" 2>"$out/$name.log"
done
python3 - "$out/alt" <<'PY'
import sys
from pathlib import Path
p=Path(sys.argv[1]); g='AACGGCCTCCCAAAGTGCTG'; near='CGCGGCCTCCCAAAGTGCTG'
rc=lambda s:s.translate(str.maketrans('ACGT','TGCA'))[::-1]
refs={'terminal':'CT'+g[2:]+'AGG','near':'CT'+near[:2]+near[4:]+'AGG','reverse':'CCT'+rc('CT'+g[2:])+'AAAA'}
for base in 'ACGT':
 refs['pam_'+base+'_plus']=g+base+'GG'
 refs['pam_'+base+'_minus']='CC'+rc(base)+rc(g)
(p/'ref.fa').write_text(''.join('>'+k+'\n'+v+'\n' for k,v in refs.items()))
(p/'guides.fa').write_text('>g\n'+g+'\n>near\n'+near+'\n')
PY
"$COLUMBA_BUILD_BIN" -r "$out/alt/index" -f "$out/alt/ref.fa" >"$out/alt/build.log" 2>&1
"$POC_BIN" "$out/alt/index" "$out/alt/ref.fa" "$out/alt/guides.fa" 0 2 0.75 >"$out/alt/out.tsv" 2>"$out/alt/run.log"
python3 "$root/research/pam_anchored_poc/check_fixtures.py" "$out"
