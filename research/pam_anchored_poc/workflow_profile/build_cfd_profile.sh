#!/usr/bin/env bash
# Builds an isolated benchmark; it links the existing source CFD function.
set -euo pipefail
[[ $# == 2 ]] || { echo 'usage: build_cfd_profile.sh OUTPUT_PATH LAZY_STATIC_RLIB' >&2; exit 2; }
out=$1; lib=$2
[[ "$out" = /* && "$lib" = /* && -r "$lib" && ! -e "$out" ]]
root=$(git rev-parse --show-toplevel)
env -u GUIX_PROFILE -u CPLUS_INCLUDE_PATH -u C_INCLUDE_PATH -u LIBRARY_PATH -u CPATH bash -c '
 set -eo pipefail
 source /gnu/store/3c2mxm30rbzvnw7qsi235mrkk3m38fym-profile/etc/profile
 cd "$1"
 rustc --edition=2021 -O research/pam_anchored_poc/workflow_profile/profile_cfd.rs --extern lazy_static="$3" -L dependency="$(dirname "$3")" -C link-arg=-Wl,--dynamic-linker=/gnu/store/ln6hxqjvz6m9gdd9s97pivlqck7hzs99-glibc-2.35/lib/ld-linux-x86-64.so.2 -C link-arg=-Wl,-rpath,/gnu/store/ln6hxqjvz6m9gdd9s97pivlqck7hzs99-glibc-2.35/lib:/gnu/store/gmv6n5vy5qcsn71pkapg2hnknyn1p7g3-gcc-12.3.0-lib/lib -o "$2"
' profile "$root" "$out" "$lib"
