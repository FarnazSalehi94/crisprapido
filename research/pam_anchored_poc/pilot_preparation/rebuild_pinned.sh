#!/usr/bin/env bash
# Rebuild the development CLI in an isolated new directory. No index is built.
set -euo pipefail
[[ $# == 2 ]] || { echo 'usage: rebuild_pinned.sh NEW_OUTPUT_DIR LOCAL_WFA2_GIT_REPOSITORY' >&2; exit 2; }
out=$1; wfa_repo=$2
[[ "$out" = /* && "$wfa_repo" = /* && ! -e "$out" ]]
repo=$(git rev-parse --show-toplevel)
commit=923e8d7a08e3005b16d18af8ff63d7a83bd7324d
wfa=df3dce4b99c37ac4f34fdcef74b774c80b7fadc4
profile=/gnu/store/3c2mxm30rbzvnw7qsi235mrkk3m38fym-profile
[[ -r "$profile/etc/profile" ]]
git -C "$repo" cat-file -e "$commit^{commit}"
git -C "$wfa_repo" cat-file -e "$wfa^{commit}"
[[ $(git -C "$repo" ls-tree "$commit" WFA2-lib | awk '{print $3}') == "$wfa" ]]
mkdir -m 700 "$out"
mkdir "$out/source" "$out/wfa"
git -C "$repo" archive "$commit" | tar -x -C "$out/source"
git -C "$wfa_repo" archive "$wfa" | tar -x -C "$out/wfa"
# The binary embeds CARGO_MANIFEST_DIR; keep $out/source at this path.
# The Cargo cache must contain all locked dependencies, including lib_wfa2.
env -u GUIX_PROFILE -u CPLUS_INCLUDE_PATH -u C_INCLUDE_PATH -u LIBRARY_PATH -u CPATH -u RUSTFLAGS bash -c '
 set -euo pipefail
 set +u
 source "$1/etc/profile"
 set -u
 cd "$2/wfa"
 make -j2 build BUILD_TOOLS=0 BUILD_EXAMPLES=0 CC=gcc CXX=g++
 export LIBRARY_PATH="$2/wfa/lib:$LIBRARY_PATH"
 export CARGO_TARGET_DIR="$2/target"
 export RUSTFLAGS="-C link-arg=-Wl,--dynamic-linker=/gnu/store/ln6hxqjvz6m9gdd9s97pivlqck7hzs99-glibc-2.35/lib/ld-linux-x86-64.so.2 -C link-arg=-Wl,-rpath,/gnu/store/ln6hxqjvz6m9gdd9s97pivlqck7hzs99-glibc-2.35/lib:/gnu/store/gmv6n5vy5qcsn71pkapg2hnknyn1p7g3-gcc-12.3.0-lib/lib"
 cd "$2/source"
 cargo build --release --locked --offline
' rebuild "$profile" "$out"
sha256sum "$out/source/Cargo.lock" "$out/wfa/lib/libwfa.a" "$out/target/release/crisprapido" | tee "$out/sha256.txt"
printf 'source commit=%s\nwfa commit=%s\n' "$commit" "$wfa" | tee "$out/revisions.txt"
printf 'Keep the source snapshot in place: its runner path is embedded in the executable.\n'
