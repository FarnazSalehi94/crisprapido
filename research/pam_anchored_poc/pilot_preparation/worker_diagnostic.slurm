#!/usr/bin/env bash
#SBATCH --job-name=crispr-pilot-diagnostic
#SBATCH --partition=workers
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:10:00
#SBATCH --output=crispr-diagnostic-%j.log
# Read-only environment diagnosis. Does not run a guide or submit another job.
set -u
umask 077
: "${SLURM_JOB_ID:?Slurm allocation required}" "${REPO:?source pilot.env before submission}" "${CRISPR_BIN:?}" "${CLI_PAYLOAD:?}" "${PAM_BIN:?}" "${PAM_PAYLOAD:?}" "${RESULT_ROOT:?}" "${SCRATCH_ROOT:?}"
status=0
fail() { printf 'FAIL: %s\n' "$*"; status=1; }
section() { printf '\n--- %s ---\n' "$*"; }
section 'identity and allocation'
date -u; hostname; id; printf 'job=%s\n' "$SLURM_JOB_ID"; scontrol show job "$SLURM_JOB_ID" || fail 'cannot read job allocation'
section 'cgroup mounts and current membership'
cat /proc/self/cgroup; awk '$0 ~ / - cgroup2? / {print}' /proc/self/mountinfo
v1=$(awk -F: '$2 ~ /(^|,)memory(,|$)/ {print $3}' /proc/self/cgroup)
if [[ -n "$v1" && -r "/sys/fs/cgroup/memory${v1}/memory.limit_in_bytes" ]]; then printf 'v1 memory limit bytes: '; cat "/sys/fs/cgroup/memory${v1}/memory.limit_in_bytes"; else printf 'v1 memory limit not available at conventional path\n'; fi
v2=$(awk -F: '$1=="0" {print $3}' /proc/self/cgroup)
for root in /sys/fs/cgroup /sys/fs/cgroup/unified; do
 if [[ -n "$v2" && -r "$root${v2}/memory.max" ]]; then printf 'v2 memory max at %s: ' "$root${v2}"; cat "$root${v2}/memory.max"; fi
done
printf 'A numeric limit here is observation only; verify that it is job-specific and enforced for all descendants.\n'
section 'candidate paths and dependencies'
for f in "$CRISPR_BIN" "$CLI_PAYLOAD" "$PAM_BIN" "$PAM_PAYLOAD" "$REPO/research/pam_anchored_poc/integration/run_snapshot.py" "$REPO/research/pam_anchored_poc/compatibility/index_manifest.py" /usr/bin/python3.7; do
 if [[ -r "$f" ]]; then stat -Lc '%A %s %n' "$f" || fail "$f stat"; else fail "unreadable $f"; fi
done
for f in "$CRISPR_BIN" "$PAM_BIN"; do [[ -x "$f" ]] || fail "not executable $f"; done
sha256sum "$CRISPR_BIN" "$CLI_PAYLOAD" "$PAM_BIN" "$PAM_PAYLOAD" || fail 'candidate hashing'
[[ "$(sha256sum "$CRISPR_BIN" | cut -d' ' -f1)" == "${CRISPR_SHA256:-missing}" ]] || fail 'CLI wrapper hash mismatch'
[[ "$(sha256sum "$CLI_PAYLOAD" | cut -d' ' -f1)" == "${CLI_PAYLOAD_SHA256:-missing}" ]] || fail 'CLI payload hash mismatch'
[[ "$(sha256sum "$PAM_BIN" | cut -d' ' -f1)" == "${PAM_SHA256:-missing}" ]] || fail 'native wrapper hash mismatch'
[[ "$(sha256sum "$PAM_PAYLOAD" | cut -d' ' -f1)" == "${PAM_PAYLOAD_SHA256:-missing}" ]] || fail 'native payload hash mismatch'
command -v python3; python3 --version
if [[ -x "$CRISPR_BIN" ]]; then timeout 20 "$CRISPR_BIN" --help > /dev/null || fail 'CLI help did not execute'; fi
if [[ -x "$PAM_BIN" ]]; then timeout 20 "$PAM_BIN" > /dev/null 2>&1; rc=$?; [[ "$rc" == 2 ]] || fail "native usage returned $rc, expected 2"; fi
if command -v ldd >/dev/null 2>&1; then ldd "$CLI_PAYLOAD" || true; ldd "$PAM_PAYLOAD" || true; fi
section 'bundle metadata only; no index rebuild or guide search'
stage=/export/local/home/salehi/crisprapido-hprc-import-hg00735-20261105-v1
for f in reference.fa index.crispr-binding.json index.bwt index.brt index.rev.brt index.sa.4 index.sa.bv.4 index.pos index.sna index.txt.bin index.cct index.meta index.fsid index.headerSN.bin; do
 if [[ -r "$stage/$f" ]]; then stat -Lc '%A %s %n' "$stage/$f" || fail "cannot stat $f"; else fail "cannot read $stage/$f"; fi
done
[[ "$(sha256sum "$stage/index.crispr-binding.json" 2>/dev/null | cut -d' ' -f1)" == 2fe715d09c5a8b79f20e59ff6cd85486b52c8f3836ab54b7d566b8a50ad9b0d4 ]] || fail 'binding manifest absent or mismatch'
fai=/moosefs/pangenomes/HPRCv2/HG00735_pat_hprc_r2_v1.0.1.fa.gz.fai
[[ "$(sha256sum "$fai" 2>/dev/null | cut -d' ' -f1)" == 7f6187771a76e8af51b983e2df4e00970dd72722dea0b9bd5c425ef817e3c26d ]] || fail 'FAI absent or mismatch'
section 'storage and permissions'
for d in "$SCRATCH_ROOT" "$RESULT_ROOT" "$stage" /tmp; do
 if [[ -d "$d" ]]; then stat -Lc '%A %U:%G %n' "$d"; df -PB1 "$d"; command -v findmnt >/dev/null 2>&1 && findmnt -T "$d" -o TARGET,SOURCE,FSTYPE,OPTIONS; else fail "directory absent $d"; fi
done
[[ -w "$SCRATCH_ROOT" ]] || fail 'scratch not writable'
[[ -w "$RESULT_ROOT" ]] || fail 'result directory not writable'
if command -v quota >/dev/null 2>&1; then quota -s || printf 'quota command failed or unavailable for this filesystem\n'; else printf 'quota command unavailable; quota remains unverified\n'; fi
section 'result'
printf 'diagnostic exit=%s; cgroup and quota output require administrator interpretation\n' "$status"
exit "$status"
