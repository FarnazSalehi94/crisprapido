# Brunello timing reconciliation and prepared 10,000-guide benchmark

This supersedes the incorrect runtime/storage statements in
`BRUNELLO_READINESS.md`; no full Brunello cohort has been launched.

## 178 hours versus 276 hours

These describe different stages, not conflicting search measurements.
Measured HG00735 m4/b1/z1/GG, first 1,000 library records:

| Stage | Measured 1,000 guides | Linear 77,441-guide scenario |
|---|---:|---:|
| Native traversal | 8,315.85 s | 178.89 h (7.45 days) |
| Independent audit | 536 s | 11.53 h |
| Both-score annotation | 3,963 s | 85.25 h |
| Sum | 12,814.85 s | **275.67 h (11.49 days)** |

The full-library **search-only** scenario is ~179 h. The complete audited,
both-score workflow scenario is ~276 h, **not 7.4 days**. At 464 assemblies:
~83,003 traversal job-hours and ~127,909 full-workflow job-hours. These are
wall-stage × job-count scenarios, not measured CPU consumption, guaranteed
throughput, or a validated full-library measurement. Whole-search wall includes
snapshot/index/reference validation and output work; not all of the 8,698 s
measured whole-search wall can be identified as traversal. Fixed startup should
not be multiplied by the number of 1,000-guide blocks for one batched call.
Guide order, output load, caches, CPUs and concurrent I/O may alter scaling.

## Storage correction

At ~5,855 rows/guide, 77,441 guides, 464 assemblies: ~210.39 billion rows.
Using measured gzip bytes/row (PAF 24.1, annotations 37.5, audit 30.6):

- PAF ~5.07 TB;
- annotations ~7.89 TB;
- audit records ~6.44 TB;
- **total ~19.40 TB**, of which authoritative PAF+annotations ~12.96 TB.

The previous claim that audit records alone occupy ~14 TB and the authoritative
archive ~5.4 TB was incorrect. No evidence may be silently discarded based on
that claim. Storage reservation must also cover uncompressed temporary output,
spools, failed attempts, scoring inputs and reference/index staging. Retention
and compression changes require lossless round-trip validation. Threshold views
must never define allele-frequency denominators; preserve all sites and explicit
missing-score reasons, reference identity, CIGARs, coordinates, strands and actual
PAMs for IMPG/Figure 1C.

## Prepared benchmark (not submitted)

`prepare_brunello10k.py` freezes 10,000 unique, unmodified guides from the
checksum-pinned 77,441-record library. Frozen files:
`brunello_pam_aware_v1/00_manifest/timing-10000-v1/{guides.fa,plan.json}`.
Output namespace: `01_pilot/timing-10000-v1/HG00735_pat/` (currently empty).
Guide sampling is first-in-file, not representative random sampling.

Exact settings: m4/b1/z1/f0.75/GG, one CPU, same verified HG00735 bound
reference, pinned native/CLI/model/scorer, independent audit, separate neural
annotation. The plan records commands and pins; **it is not a validated
restartable execution wrapper**. Do not reuse the defective old Brunello
full-run worker. Before submission test the supervisor/locking/stage-checkpoint
wrapper and verify reference .fai provenance, free scratch ≥60 GB and output
headroom ≥150 GB. Use one node/task/CPU, 56 GiB request, 48 GiB sampled RSS stop,
72-hour outer limit; estimated traversal 23.1 h + audit 1.49 h + annotation
11.01 h (~35.6 h total before further overhead).

Cluster availability was checked while SBDS and recovery jobs were active.
SBDS/rebuilding remain priority, so **no Brunello benchmark job was submitted**.
It should be scheduled only after an isolated resource slot and the remaining
wrapper checks are satisfied. Full cohort remains review-gated.
