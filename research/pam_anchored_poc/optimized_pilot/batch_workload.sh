#!/usr/bin/env bash
#SBATCH --job-name=crispr-hg00735-pilot
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --output=crispr-pilot-%j.slurm.log
# Preparation only: do not submit until the report's conditions are reviewed.
set -Eeuo pipefail
SECONDS=0
umask 077
: "${SLURM_JOB_ID:?Slurm job required}" "${REPO:?absolute checkout path required}" "${CRISPR_BIN:?fresh experimental CLI path required}" "${CRISPR_SHA256:?pinned CLI SHA256 required}" "${PAM_BIN:?experimental native executable required}" "${PAM_SHA256:?pinned native wrapper SHA256 required}" "${PAM_PAYLOAD:?native payload required}" "${PAM_PAYLOAD_SHA256:?pinned native payload SHA256 required}" "${CLI_PAYLOAD:?CLI payload required}" "${CLI_PAYLOAD_SHA256:?pinned CLI payload SHA256 required}" "${GUIDE:?one 20-base ACGT guide required}" "${RESULT_ROOT:?existing private output parent required}" "${SCRATCH_ROOT:?existing scratch directory required}"
[[ "$GUIDE" == GCCAGGTACAGGTTCCGTCC ]]
[[ "$REPO" = /* && "$CRISPR_BIN" = /* && "$PAM_BIN" = /* && "$PAM_PAYLOAD" = /* && "$CLI_PAYLOAD" = /* && "$RESULT_ROOT" = /* && "$SCRATCH_ROOT" = /* ]]
STAGE=/export/local/home/salehi/crisprapido-hprc-import-hg00735-20261105-v1
FAI=/moosefs/pangenomes/HPRCv2/HG00735_pat_hprc_r2_v1.0.1.fa.gz.fai
[[ $(sha256sum "$STAGE/index.crispr-binding.json" | cut -d' ' -f1) == 2fe715d09c5a8b79f20e59ff6cd85486b52c8f3836ab54b7d566b8a50ad9b0d4 ]]
[[ $(sha256sum "$FAI" | cut -d' ' -f1) == 7f6187771a76e8af51b983e2df4e00970dd72722dea0b9bd5c425ef817e3c26d ]]
[[ $(sha256sum "$CRISPR_BIN" | cut -d' ' -f1) == "$CRISPR_SHA256" ]]
[[ $(sha256sum "$PAM_BIN" | cut -d' ' -f1) == "$PAM_SHA256" ]]
[[ $(sha256sum "$PAM_PAYLOAD" | cut -d' ' -f1) == "$PAM_PAYLOAD_SHA256" ]]
[[ $(sha256sum "$CLI_PAYLOAD" | cut -d' ' -f1) == "$CLI_PAYLOAD_SHA256" ]]
# Memory is monitored by the required parent supervisor, not a hard cgroup limit.
: "${MONITORED_PILOT_SUPERVISED:?supervisor required}"
[[ $MONITORED_PILOT_SUPERVISED == "$SLURM_JOB_ID" ]]
fail_preflight() { printf 'pilot preflight failed: %s\n' "$*" >&2; exit 1; }
free_bytes=$(df -PB1 "$SCRATCH_ROOT" | awk 'NR==2 {print $4}')
(( free_bytes >= 60000000000 )) || fail_preflight 'less than 60 GB available on scratch filesystem'
out="$RESULT_ROOT/hg00735-batch-$SLURM_JOB_ID"
mkdir "$out"
printf 'INCOMPLETE\n' > "$out/status"
scratch=$(mktemp -d "$SCRATCH_ROOT/crispr-pilot-$SLURM_JOB_ID.XXXXXXXX")
cleanup() { rc=$?; trap - EXIT; printf 'elapsed_script_seconds=%s exit=%s\n' "$SECONDS" "$rc" >> "$out/progress.log"; rm -rf -- "$scratch"; if (( rc != 0 )); then printf 'INCOMPLETE exit=%s\n' "$rc" > "$out/status"; fi; exit "$rc"; }
trap cleanup EXIT
trap 'exit 143' TERM
trap 'exit 130' INT
export TMPDIR="$scratch" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
# Maximum individual output file: 10 GiB; scratch requires separate site quotas.
ulimit -f 10485760
cp "$REPO/research/pam_anchored_poc/optimized_pilot/batch_guides.fa" "$out/guide.fa"
{ date -u; hostname; printf 'job=%s memory=sampled-process-tree-48GiB-not-hard-limit scratch=%s\n' "$SLURM_JOB_ID" "$scratch"; sha256sum "$CRISPR_BIN" "$PAM_BIN" "$CLI_PAYLOAD" "$PAM_PAYLOAD" "$STAGE/index.crispr-binding.json" "$FAI"; } > "$out/progress.log"
cd "$REPO"
sha256sum -c research/pam_anchored_poc/pilot_preparation/dependencies.sha256 > "$out/helper-checksums.txt"
sha256sum -c research/pam_anchored_poc/monitored_pilot/timing.sha256 >> "$out/helper-checksums.txt"
sha256sum -c research/pam_anchored_poc/optimized_pilot/runtime.sha256 >> "$out/helper-checksums.txt"
source research/pam_anchored_poc/monitored_pilot/portable_timing.sh
printf '%s preflight full binding check\n' "$(date -u +%FT%TZ)" >> "$out/progress.log"
python3 research/pam_anchored_poc/compatibility/index_manifest.py verify --reference "$STAGE/reference.fa" --index "$STAGE/index" > "$out/preverify.json"
printf '%s start three-guide workflow\n' "$(date -u +%FT%TZ)" >> "$out/progress.log"
run_timed "$out/workflow.time" "$CRISPR_BIN" --search-backend pam-anchored --pam-search-bin "$PAM_BIN" --columba-index "$STAGE/index" -r "$STAGE/reference.fa" -G "$out/guide.fa" -t 1 -m 5 -b 1 -z 2 -f 0.75 -p GG > "$out/output.paf.incomplete" 2> "$out/workflow.stderr"
printf '%s postflight full binding check\n' "$(date -u +%FT%TZ)" >> "$out/progress.log"
python3 research/pam_anchored_poc/compatibility/index_manifest.py verify --reference "$STAGE/reference.fa" --index "$STAGE/index" > "$out/postverify.json"
[[ $(sha256sum "$STAGE/index.crispr-binding.json" | cut -d' ' -f1) == 2fe715d09c5a8b79f20e59ff6cd85486b52c8f3836ab54b7d566b8a50ad9b0d4 ]]
[[ $(sha256sum "$FAI" | cut -d' ' -f1) == 7f6187771a76e8af51b983e2df4e00970dd72722dea0b9bd5c425ef817e3c26d ]]
printf '%s independent output eligibility audit\n' "$(date -u +%FT%TZ)" >> "$out/progress.log"
python3 research/pam_anchored_poc/optimized_pilot/audit_batch.py --reference "$STAGE/reference.fa" --fai "$FAI" --guides "$out/guide.fa" --targets "$REPO/research/pam_anchored_poc/optimized_pilot/batch_targets.json" --paf "$out/output.paf.incomplete" > "$out/audit.json"
[[ $(sha256sum "$STAGE/reference.fa" | cut -d' ' -f1) == bae8df3c5cf8ca96d3316b59eaa63fa37b7d3b58551c8a0700f2c530dc1e4f8b ]]
[[ $(sha256sum "$FAI" | cut -d' ' -f1) == 7f6187771a76e8af51b983e2df4e00970dd72722dea0b9bd5c425ef817e3c26d ]]
sha256sum -c research/pam_anchored_poc/pilot_preparation/dependencies.sha256 >> "$out/helper-checksums.txt"
sha256sum -c research/pam_anchored_poc/monitored_pilot/timing.sha256 >> "$out/helper-checksums.txt"
sha256sum -c research/pam_anchored_poc/optimized_pilot/runtime.sha256 >> "$out/helper-checksums.txt"
python3 research/pam_anchored_poc/optimized_pilot/phase_metrics.py "$out/workflow.stderr" "$out/workflow.time" "$out/audit.json" > "$out/search-metrics.json"
python3 research/pam_anchored_poc/optimized_pilot/extract_batch_metrics.py "$out/workflow.stderr" "$out/audit.json" "$out/search-metrics.json" > "$out/per-guide-metrics.json"
baseline=/export/local/home/salehi/crisprapido-resource-pilot-results/hg00735-optimized-2871283/output.paf
[[ $(sha256sum "$baseline" | cut -d' ' -f1) == a28d0574a915a57d431851ab681b0602b8f5fb1740d1c8efeb6e3b4a378f23be ]]
awk -F '\t' '$1=="pilot"' "$out/output.paf.incomplete" > "$scratch/pilot.paf"
python3 research/pam_anchored_poc/optimized_pilot/compare_paf.py --baseline "$baseline" --candidate "$scratch/pilot.paf" --out "$out/pilot-comparison.json"
python3 - "$out/pilot-comparison.json" <<'PYCOMPARE'
import json,sys
x=json.load(open(sys.argv[1]));assert x['representative_alignment_differences']==0 and x['cfd_differences']==0, 'pilot results differ from 2871283'
PYCOMPARE
mv "$out/output.paf.incomplete" "$out/output.paf"
printf 'COMPLETE eligibility audited; recall not assessed\n' > "$out/status"
printf '%s finished\n' "$(date -u +%FT%TZ)" >> "$out/progress.log"
