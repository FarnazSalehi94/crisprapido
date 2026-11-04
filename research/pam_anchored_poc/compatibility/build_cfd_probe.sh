#!/usr/bin/env bash
set -euo pipefail
: "${CRISPR_LIB:?set CRISPR_LIB to the built crisprapido .rlib}"
: "${DEPS_DIR:?set DEPS_DIR to its dependency directory}"
read -r -a extra <<< "${RUSTC_EXTRA_FLAGS:-}"
"${RUSTC:-rustc}" --edition=2021 "$(dirname "$0")/cfd_probe.rs" --extern "crisprapido=$CRISPR_LIB" -L "dependency=$DEPS_DIR" "${extra[@]}" -o "${OUTPUT:-/tmp/crispr-cfd-probe}"
