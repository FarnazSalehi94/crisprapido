# Sourced by the monitored workload under Bash. SECONDS is a Bash builtin.
# Records whole seconds for the command only, on both success and failure.
run_timed() {
    local timing_file=$1 start=$SECONDS rc
    shift
    if "$@"; then rc=0; else rc=$?; fi
    printf 'elapsed_seconds=%s\nexit_code=%s\nclock=bash_SECONDS_wall_seconds\n' "$((SECONDS - start))" "$rc" > "$timing_file"
    return "$rc"
}
