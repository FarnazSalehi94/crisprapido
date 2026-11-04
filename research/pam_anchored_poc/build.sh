#!/usr/bin/env bash
set -euo pipefail
# Provide a separately built, unmodified Columba edcc7ce checkout and build.
: "${COLUMBA_DIR:?set COLUMBA_DIR to the Columba source checkout}"
: "${COLUMBA_BUILD:?set COLUMBA_BUILD to its vanilla build directory}"
CXX=${CXX:-g++}
PYTHON=$(command -v "${PYTHON:-python3}")
case "$PYTHON" in /*) ;; *) echo "Python interpreter must resolve to an absolute path" >&2; exit 1;; esac
HELPER=$(cd "$(dirname "$0")/compatibility" && pwd)/index_manifest.py
"$CXX" -std=gnu++14 -O3 -mpopcnt -DALPHABET=5 -DTHIRTY_TWO -DHAVE_ZLIB "-DINDEX_BINDING_HELPER=\"$HELPER\"" "-DINDEX_BINDING_PYTHON=\"$PYTHON\"" -include limits -include cstdint \
 -I"$COLUMBA_DIR/src" -I"$COLUMBA_BUILD/_deps/parallel-hashmap-src" -I"$COLUMBA_BUILD/_deps/fmtlib-src/include" \
 "$(dirname "$0")/poc.cpp" "$COLUMBA_BUILD/libAuxiliaryClasses.a" "$COLUMBA_BUILD/_deps/fmtlib-build/libfmt.a" -lz -pthread -o "${OUTPUT:-/tmp/crispr-pam-poc}"
