#!/usr/bin/env bash
# Bounded research validation only. Build/install tools separately; no deployment.
set -euo pipefail
[[ $# == 1 && "$1" = /* && ! -e "$1" ]] || { echo 'usage: reproduce.sh NEW_ABSOLUTE_EXTERNAL_DIRECTORY' >&2; exit 2; }
: "${CRISPR_CLI:?set absolute development CLI path}"
: "${CRISPR_OLD_CLI:?set preserved CFD baseline CLI path}"
: "${PAM_SEARCH:?set separately built unpruned PAM search path}"
: "${COLUMBA_BUILDER:?set compatible vanilla/32-bit Columba builder path}"
: "${MODEL_PYTHON:?set absolute inference environment interpreter}"
: "${MODEL_ROOT:?set trusted checksum-pinned model root}"
python=${RESEARCH_PYTHON:-python3}
repo=$(cd "$(dirname "$0")/../.." && pwd)
out=$1
mkdir -m 700 "$out"
cd "$repo"
export PYTHONDONTWRITEBYTECODE=1
extra=()
if [[ -n ${RETAINED_PANEL:-} ]]; then extra+=(--retained-panel "$RETAINED_PANEL"); fi
if [[ -n ${OFFICIAL_PYTHON:-} || -n ${OFFICIAL_ROOT:-} ]]; then
  : "${OFFICIAL_PYTHON:?both official options are required}"
  : "${OFFICIAL_ROOT:?both official options are required}"
  extra+=(--official-python "$OFFICIAL_PYTHON" --official-root "$OFFICIAL_ROOT")
fi
"$python" research/reproducibility/validate_workflow.py \
 --cli "$CRISPR_CLI" --old-cli "$CRISPR_OLD_CLI" --native "$PAM_SEARCH" --builder "$COLUMBA_BUILDER" \
 --model-python "$MODEL_PYTHON" --model-root "$MODEL_ROOT" --out "$out/workflow" "${extra[@]}" \
 > "$out/workflow.stdout" 2> "$out/workflow.stderr"
export CRISPR_SCORING_CLI="$CRISPR_CLI" CRISPR_MODEL_PYTHON="$MODEL_PYTHON" CRISPR_MODEL_ROOT="$MODEL_ROOT"
"$python" -m unittest discover -s research/crispr_bulge -p 'test_*.py' -v > "$out/scoring-tests.log" 2>&1
"$python" -m unittest discover -s research/reproducibility -p 'test_*.py' -v > "$out/validator-tests.log" 2>&1
"$python" -m unittest discover -s research/pam_anchored_poc/cfd_feasibility -p 'test_*.py' -v > "$out/cfd-model-tests.log" 2>&1
for test in test_binding test_options test_mixed_gaps; do
  "$python" "research/pam_anchored_poc/compatibility/$test.py" "$COLUMBA_BUILDER" "$PAM_SEARCH" "$out/$test" \
    > "$out/$test.log" 2>&1
done
echo "Complete bounded validation and Python regressions: $out (Rust tests are a separate build step)"
