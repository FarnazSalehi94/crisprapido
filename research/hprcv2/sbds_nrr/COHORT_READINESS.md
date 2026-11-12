# SBDS 466-genome readiness: ledger, expansion tooling, multi-assembly test

Baseline `195e2cf`; no search algorithm, scoring formula or production binary
changed; historical outputs preserved. Root:
`/moosefs/salehi/HPRCv2_CRISPRapido/runs/sbds_ngg_nrr_v1/`. **The complete
466-genome search was NOT launched** (explicitly out of scope for this pass).

## 1. Binding ledger audit

`00_manifest/binding-ledger-audit-v1.json` (tooling: `ledger.py`, committed):

- **12 assemblies have trusted full-content index evidence** (recovery
  report 2026-10-28/29): CHM13, GRCh38, HG00128_hap2, HG00735_pat, HG01081_pat,
  HG01496_pat, HG02055_mat, HG02630_mat, HG03369_hap2, HG04184_mat, NA18906_mat,
  NA21309_pat — all 12 index prefixes still present on disk.
- **454 assemblies have NO trusted content evidence** and require deep
  verification (full component reads) before any search. Metadata presence is
  never accepted as verification.
- Native binding sidecars ready today: HG00735_pat, HG01496_pat (prior imports)
  plus the two new bindings below. This audit authorizes nothing by itself;
  every job re-verifies its own binding at execution time.

## 2. Repaired expansion tooling (independently tested)

The withdrawn Brunello-draft defects are fixed in new SBDS-specific tools
(`prepare_binding.py`, `cohort_worker.py`, `make_cohort_plan.py`,
`verify_cohort.py`, `submit_cohort.slurm`, `ledger.py`):

- `prepare_binding.py`: correct file handles (no closed-file fsync, no string
  `.open()`), correct digest semantics (ledger `fasta_sha256` is the
  **compressed** .gz hash; decompressed stream hashed separately), trusted
  component drift comparison, exclusive stage, committed binding.json.
- `cohort_worker.py`: stable per-PAM search outputs with checksummed
  SEARCH_COMPLETE checkpoints; **restart reuses verified searches** and re-runs
  only cheap audit/scoring; uncommitted partial PAFs are removed, never
  trusted; exclusive attempts; SUCCESS markers re-verified with checksums;
  per-assembly pins (a job hashes only its own 3 GB reference, fai, binding,
  scoring index); one CPU; whole-job supervisor memory/disk guards.
- Six committed regression tests (`test_cohort_worker.py`) cover: first-PAM
  failure + clean restart; mid-run failure with checkpoint reuse (proven again
  on real data: v3 rerun reused both assemblies' four searches and finished in
  100 s); tampered completed-output detection; stage-timeout process-group kill
  (no leaked processes); bad plan/binding rejection; same-assembly lock
  exclusivity. Failure detection is real: two array attempts were correctly
  stopped by the scratch guard (octopus05) and by the sbatch comma-export
  issue; their FAILED.json/attempt dirs are preserved.

## 3. Controlled multi-assembly test (verified bindings)

Two new assemblies prepared with the repaired binding tool (jobs 2871770/71,
exit 0): **HG02055_mat** (original MooseFS index) and **NA18906_mat** (rebuilt
/export index). Both passed the trusted component drift comparison and
compressed-FASTA/fai hash checks.

Plan `cohort-plan-test-v3.json` (SHA256
`48f3c3271e5373b9a996d61e8a73b944a3743c73c263f16c815194f96c52b859`) pinned the
frozen NRR scorer bundle (`scorer-v2`), byte-identical to the scorer validated
against official inference on 42,425 pilot predictions + 8,272 domain fixtures.
An earlier v2 plan mistakenly pointed at the pre-NRR attested runtime source;
its outputs are preserved under SUPERSEDED markers, not deleted or reused.

Cohort array **2871780**, one CPU, 56G, exit 0, supervisor completed:

| Assembly | NGG | NAG | NGA | NAA | Total hits | CFD defined | CB scored | CB N/A |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| HG02055_mat | 38,037 | 46,380 | 42,870 | 65,091 | 192,378 | 192,378 | 42,465 | 149,913 |
| NA18906_mat | 38,180 | 46,421 | 42,851 | 64,696 | 192,148 | 192,148 | 42,389 | 149,759 |

All eight (assembly, PAM) sets passed the independent audit with explicit
allowed-PAM suffixes (PASS markers under `04_qc/...`), every row CFD-defined,
neural N/A carrying explicit reasons (24-column limit / gapped mismatch domain
/ multibase bulge), actual three-base PAMs retained in annotations. Full search
runs (2871778) took 22.5–29.3 min per assembly on one CPU; peak sampled tree
RSS in the v3 scoring rerun 656,356 KiB; pilot search peaks were ~13.6–14.4 GiB
(one search process at a time; supervisor limit 48 GiB).

## 4. Measured resource model for the complete 466 (scenario, not a launched cost)

Per assembly, one CPU: ~28–30 min four-PAM search (m5/b1/z2, this guide) +
~2 min audit/scoring; binding preparation adds ~125–135 s (decompress+verify
15 s, **full component read ~9 s warm / up to minutes cold**, scoring index
~100 s). Storage: ~10 MB compressed PAF+annotations per assembly
(≈4.75 GB total), reference bindings ~3.1 GB each (≈1.45 TB on /export if all
are kept; they are reusable and could be re-derived, but deletion would
re-verify nothing).

- **Verification campaign for the 454 unverified assemblies:** ~454 × 12.4 GB
  ≈ **5.6 TB of component reads** (~17–20 CPU-hours at measured rates; wall
  time bounded by aggregate read bandwidth, realistically hours to a day at
  moderate concurrency) plus reference decompression/index per assembly
  (~454 × ~120 s ≈ 15 CPU-hours). **Expected EIO failures:** the historical
  cohort showed damaged/unreadable components on multiple index sets; deep
  verification will surface them assembly by assembly — those assemblies stop
  with explicit failure markers and need rebuild/recovery decisions, never
  zero-filling.
- **Search campaign for 466:** ≈ 466 × ~30 min ≈ **~233 one-CPU job-hours**
  (~15 h wall at %16, ~7.5 h at %32), ~4.75 GB compressed outputs, scratch
  headroom ≥60 GB per job (guards refuse nodes like the full octopus05).

## 5. Recommended Slurm submission settings (NOT submitted)

```bash
# Stage 1 — deep verification + binding preparation (restartable, safe to resubmit):
sbatch --array=0-453%12 --cpus-per-task=1 --mem=8G --time=02:00:00 \
  --nodelist=<nodes with /scratch headroom and /export access> \
  --export=ALL,SPEC_FILE=<per-assembly args>,... prepare_binding-wrapper.slurm
# Stage 2 — four-PAM search per assembly (after ALL bindings committed):
sbatch --array=0-465%16 --cpus-per-task=1 --mem=56G --time=06:00:00 \
  --export=ALL,PLAN=...,ASSEMBLIES_FILE=... submit_cohort.slurm
```

`%12`/`%16` concurrency respects the shared cluster (current load ~40 CPUs in
use); the array is restartable — completed assemblies verify checksums and
skip, failed ones resume from checkpoints. Do not pass comma lists through
sbatch `--export`. Node choice must respect the 60 GB /scratch guard.

## 6. Remaining blockers before the complete 466 launch

1. **454 assemblies lack trusted content evidence** — the deep verification
   campaign above must run first (user authorization for its I/O budget and
   for expected EIO-rebuild decisions required).
2. **Noncanonical biological validity remains unestablished** (computational
   agreement only); reports must keep that caveat.
3. Historical EIO assemblies (HG03784_hap1, NA18983_hap1 SBDS; ten Brunello
   artifacts) and the four historical no-perfect-match cases need explicit
   input investigation when reached; never treated as zero-hit.
4. Historical SBDS comparisons (chr7/chr12/chr17, HG00321 duplication) are
   NOT resolved by this test and remain part of the full-run analysis.
