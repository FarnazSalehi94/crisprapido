# Brunello execution-wrapper review — 2026-11-15

Independent read-only review by `scout` / `openai-codex/gpt-6-luna`, integrated
by the main Astra session. **No benchmark or full cohort submitted.**

## Frozen scope and measured basis

- Library: 77,441 unique guides; manuscript cohort: 464 HPRCv2 assemblies.
- Exact search: m4/b1/z1/f0.75/GG, corrected PAM-aware backend, one CPU.
- Both CFD and CRISPR-Bulge; preserve every eligible hit and explicit N/A reasons.
- Existing 10k plan: `/moosefs/salehi/HPRCv2_CRISPRapido/runs/brunello_pam_aware_v1/00_manifest/timing-10000-v1/plan.json`.
- First 10,000 guides are deterministic first-in-file, not a random sample.
- 1k measured stages: traversal 8,315.85 s, audit 536 s, annotation 3,963 s.
- Full-library linear scenario: 178.89 + 11.53 + 85.25 = 275.67 hours per
  assembly; 127,909 workflow job-hours across 464. These are not measured CPU
  hours or validated full-library runtimes.
- 10k scenario: ~35.6 hours plus setup/output overhead; proposed 72 h limit,
  one CPU, 56 GiB requested, 48 GiB sampled tree-RSS stop, ≥60 GB free scratch,
  ≥150 GB persistent output headroom. No resource slot is reserved yet.
- Cohort storage scenario: 12.96 TB authoritative PAF+annotations and 6.44 TB
  audit records; 19.40 TB total before uncompressed spools/retries/staging.
  Current 34.7 TB free is not a capacity reservation.

## Blocking implementation work

The frozen plan is a command/pin inventory, **not an executable restart plan**.
The withdrawn `brunello_expanded_benchmark/full_run_worker.py` must not be used.
Do not pass the plan to the SBDS worker: schema, root, PAMs, parameters and
scoring timeout differ. Its 1,800 s timeout is insufficient for the proposed
10k annotation stage (~39,630 s); the frozen command requests 43,200 s.

Implement a bounded benchmark-specific worker before any full-cohort worker:

1. Strict schema, explicit benchmark authorization, input/tool/model/worker pins,
   output namespace and exact guide count/configuration; refuse cohort use.
2. Exclusive nonblocking per-assembly lock; immutable per-attempt directories.
3. Stage completion checkpoints bind scientific identity, command/pins, output
   paths, bytes and SHA256. Verify all before reuse; never accept a marker alone.
4. One batched search for the bounded 10k timing measurement. Search, audit and
   annotation checkpoints independently reusable; preserve failed attempts.
   Full-cohort guide chunking needs separate coverage/order/merge validation.
5. Audit through a scratch reference alias with pinned copied `.fai`; verify
   `.fai` provenance explicitly. Do not generate an index beside a frozen input.
6. Whole-process-tree supervisor with timeout, memory and scratch/output floors;
   process CPU user/system accounting as well as wall/native-stage measurements.
7. Independently audited row accounting, scoring stdout/bundle equality, explicit
   score availability, and final checksummed SUCCESS covering all authoritative
   outputs. Final SUCCESS requires supervisor success as well as worker success.

Required failure tests: crash before/during search, crash around checkpoint
publication, changed pins/output tampering, audit/scoring failure, row loss,
concurrent lock acquisition, timeout/low-space/memory shutdown with descendants,
restart reuse without duplicate rows, and final completion checksum verification.
For later chunked operation additionally test exact guide coverage, stable order,
partial merge recovery and byte-identical deterministic reassembly.

## Throughput and manuscript preservation

Prefer independent assembly jobs and bounded large guide batches. Guide-level
sharding may reduce per-assembly wall time but repeats integrity setup and index
memory; benchmark it without changing search algorithms, formulas, or eligible
hit policy. SBDS has priority; no performance/timing job may compete for its
reserved node resources or materially increase shared-storage pressure.

Preserve actual PAMs, coordinates, strands, CIGARs, both scores, explicit missing
score reasons, assembly/reference/index provenance and guide identity. Re-extract
IMPG flanks only from frozen reference identities. Unmapped/deleted loci are
assessability states, never zero-hit alleles. Score thresholds are downstream
views, not allele-frequency denominators. Extended CRISPRme Figure 1C remains
post-cohort analysis, not a justification for filtering authoritative results.
