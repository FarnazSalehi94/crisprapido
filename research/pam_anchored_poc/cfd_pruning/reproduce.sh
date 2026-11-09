#!/usr/bin/env bash
# Site-specific bounded reproduction. Never overwrites an existing output root.
# No Slurm submission, whole-genome query, or production binary modification.
set -Eeuo pipefail
[[ $# == 1 && "$1" = /* && ! -e "$1" ]] || { echo 'usage: reproduce.sh NEW_ABSOLUTE_OUTPUT_ROOT' >&2; exit 2; }
out=$1
research_python=${RESEARCH_PYTHON:-python3}
"$research_python" -c 'import sys; assert sys.version_info >= (3, 9), "research harness requires Python 3.9+"'
repo=$(git rev-parse --show-toplevel)
profile=/gnu/store/3c2mxm30rbzvnw7qsi235mrkk3m38fym-profile
base=/export/local/home/salehi/crisprapido-pilot-rebuild-f7de425-v1
columba=/moosefs/raid5/salehi/crispr2/crispr-progress/columba-crispr-benchmark/columba
builder=/export/local/home/salehi/crisprapido-index-recovery-v1-2026-10-28/rebuilt_HG00128_hap2_attempt02/tools/columba_build
mkdir -m 700 "$out"
cp -a --reflink=auto "$base/target" "$out/target"
# Do not inherit incompatible system/Guix include paths or an old PATH rustc.
env -u GUIX_PROFILE -u CPLUS_INCLUDE_PATH -u C_INCLUDE_PATH -u LIBRARY_PATH -u CPATH -u RUSTFLAGS bash -c '
 set -eo pipefail
 source "$1/etc/profile"
 export PATH="$1/bin:$PATH" RUSTC="$1/bin/rustc" CARGO_BUILD_JOBS=1
 export LIBRARY_PATH="$3/wfa/lib:$LIBRARY_PATH" CARGO_TARGET_DIR="$2/target"
 export RUSTFLAGS="-C link-arg=-Wl,--dynamic-linker=/gnu/store/ln6hxqjvz6m9gdd9s97pivlqck7hzs99-glibc-2.35/lib/ld-linux-x86-64.so.2 -C link-arg=-Wl,-rpath,/gnu/store/ln6hxqjvz6m9gdd9s97pivlqck7hzs99-glibc-2.35/lib:/gnu/store/gmv6n5vy5qcsn71pkapg2hnknyn1p7g3-gcc-12.3.0-lib/lib"
 cd "$4"
 cargo test --release --locked --offline -- --test-threads=1 > "$2/rust-tests.log" 2>&1
 PYTHONDONTWRITEBYTECODE=1 /usr/bin/python3 research/pam_anchored_poc/cfd_pruning/build.py \
   --out "$2/native-build" --columba-dir "$5" --columba-build "$5/build_Vanilla" --cxx "$1/bin/g++"
' reproduce "$profile" "$out" "$base" "$repo" "$columba"
cd "$repo"
export PYTHONDONTWRITEBYTECODE=1
export CRISPR_CFD_CLI="$out/target/release/crisprapido" CRISPR_CFD_NATIVE="$out/native-build/native"
export CRISPR_CFD_UNPRUNED="$base/native/wrapper" COLUMBA_BUILDER="$builder"
# The native/CLI regression tests also run without retained genomic extracts.
"$research_python" -m unittest discover -s research/pam_anchored_poc/cfd_pruning -p 'test_*.py' -v > "$out/native-tests.log" 2>&1
"$research_python" -m unittest discover -s research/pam_anchored_poc/cfd_feasibility -p 'test_*.py' -v > "$out/model-regressions.log" 2>&1
"$research_python" research/pam_anchored_poc/cfd_pruning/validate.py --out "$out/validation" \
 --cli "$CRISPR_CFD_CLI" --old-cli "$base/target/release/crisprapido" \
 --native "$CRISPR_CFD_NATIVE" --old-native "$CRISPR_CFD_UNPRUNED" --builder "$builder" \
 --real-panel-root "${BOUNDED_REAL_ROOT:-/tmp/crispr-v2-3mb}"
