# SBDS deep verification — controlled batch 1 (16 assemblies)

Baseline `fd813db`. **No rebuild or overwrite of any original index; no
uncontrolled campaign and no searches launched.** Root:
`/moosefs/salehi/HPRCv2_CRISPRapido/runs/sbds_ngg_nrr_v1/`; batch records in
`00_manifest/verification-batch-1/`, updated ledger
`00_manifest/binding-ledger-audit-v2.json`. All previous results unchanged.

## Storage preflight (task 1)

- `/export` (/dev/md0p1): 11 TB free of 14 TB. `/moosefs`: 45 TB free of 211 TB.
- Worker scratch: octopus07 318 G, octopus08 392 G, octopus06 136 G free;
  **octopus05 only 16 G free — excluded** (guards would refuse it anyway).
  octopus09–11 fully allocated by other users during the batch.

## Batch (tasks 2–7)

Deterministic 16-assembly selection: the 7 historically notable assemblies
(four historical no-perfect-match cases HG02258_mat/HG03225_hap1/HG03195_hap1/
HG03834_mat, the duplicated-context HG00321_hap2, and the two historical EIO
result assemblies HG03784_hap1/NA18983_hap1) plus a 9-assembly spread
(every N-th of the remaining unverified list). Stage-1 array **2871790**,
`--array=0-15%4`, one CPU, 8G, nodes octopus06/07/08, one task per assembly,
restartable; committed assemblies recheck and skip.

**Result: 6 BOUND, 10 FAILED — every failure is EIO (task 7).**

| Outcome | Assemblies |
|---|---|
| BOUND (deep-verified, reusable) | HG00097_hap1, HG01258_mat, HG02027_mat, HG02984_pat, HG03710_pat, NA18983_hap1 |
| EIO (originals untouched) | HG00321_hap2, HG005_pat, HG02258_mat, HG02451_pat, HG03195_hap1, HG03225_hap1, HG03784_hap1, HG03834_mat, NA18747_hap2, NA19443_hap1 |

Per-component characterization on two EIO assemblies (separate bounded reads;
first failing offset per file, not full corruption maps): HG03784_hap1 fails
in `rev.brt` (~0.54 GB read) and `sa.4` (~0.67 GB); HG00321_hap2 fails in
`brt` (~1.28 GB) and `txt.bin` (~2.89 GB). Different components and offsets
per assembly — genuine distributed storage degradation, consistent with the
2026-10-28/29 recovery findings, not a single bad node or file.

**Historical-anomaly correlation:** 6 of the 7 historically notable assemblies
have EIO-corrupted indexes (all four historical no-perfect-match cases plus
HG00321_hap2 and HG03784_hap1). NA18983_hap1's *index* is intact — its
historical EIO was on the results PAF. This is strong evidence the historical
SBDS anomalies were at least partly caused by degraded indexes, but it does
not by itself explain each historical result; final attribution awaits
post-search comparison.

## Measurements (task 4)

- Successful verification per assembly: **12.16 GB read** (12 components +
  compressed FASTA + fai), **147–172 s** (median 154 s): component full read
  ~30 s cold, decompress+hash ~28 s, scoring reference index ~114 s.
  Total batch: 73.0 GB read, 939 s of one-CPU work for the 6 successes.
- Disk consumption per committed binding: **~3.06 GB** new bytes
  (decompressed reference + scoring index + sidecars; index components are
  symlinks, not copies). The batch's measured `written_bytes` of 14.31 GB per
  assembly incorrectly followed symlinks — recorded here as a measurement bug,
  not storage cost.
- EIO assemblies fail fast (57–102 s), preserving partial stages and FAILED
  records with `failure_kind=EIO`; nothing was rebuilt or overwritten.
- Reusability confirmed (task 5): re-running a committed assembly prints
  `RECHECKED` after re-hashing the full 3 GB reference against its binding
  record.

## Extrapolation to the remaining 438 (task 8) — scenarios, not measurements

- Verification I/O: ≤ 438 × 12.16 GB ≈ **5.3 TB reads** worst case; EIO
  assemblies stop early so actual is less. One-CPU time ≈ 438 × ~155 s ≈
  **~19 CPU-hours** (~2.6 h wall at %12 concurrency). Storage for all
  successful bindings: ~3.06 GB × (number bound) — up to **~1.4 TB** on
  MooseFS if every assembly binds (45 TB free: fits).
- **Expected EIO prevalence:** the spread sample (non-notable assemblies)
  failed 4/9 ≈ **44%**. If representative, ~190 of the remaining 438 need full
  index rebuilds from frozen FASTAs before they can be searched at all. The
  2026-10-28 recovery rebuilt 6 indexes and text-repaired 3; its measured
  rebuild costs were not recorded per-assembly, so the rebuild campaign's
  cost is **not yet estimated from data** and requires its own bounded pilot
  before any commitment.
- Historical severity recheck: the recovery's original scan found 17 EIO
  components across 12 known-bad sets; this batch found 10 more bad sets in a
  16-assembly sample, so degradation is more widespread than the original
  12-set inventory implied.

## Updated ledger (task 9)

`binding-ledger-audit-v2.json`: 12 trusted-recovery assemblies + **8 committed
bindings** (6 new + HG02055_mat/NA18906_mat from the cohort test) = assemblies
currently search-ready: **20 of 466**; 10 batch EIO failures recorded; 438
not yet attempted. `ledger.py` now merges trusted evidence, committed
bindings, and batch failures; metadata presence alone still authorizes
nothing.

## Safety verdict (is it safe to proceed with the complete verification campaign?)

**Yes for the campaign itself** — it is bounded (~5.3 TB reads, ~19 CPU-hours,
restartable, preserves all originals and evidence, fits storage). **No
implication that a complete 466-genome search will then be possible**: at the
observed prevalence, a large rebuild campaign (~hundreds of Columba index
builds from frozen FASTAs, cost unmeasured) stands between verification and a
complete search, and MooseFS degradation appears ongoing — verified bindings
should be re-checked at search time (the worker already re-hashes every
component per job). Recommended next steps, each requiring explicit
authorization:

1. Complete the verification campaign in batches with the exact tooling used
   here (same %4–%12 limits; no octopus05).
2. A bounded rebuild pilot (a few EIO assemblies) to measure per-assembly
   rebuild time/storage before committing to hundreds.
3. Raise the storage-degradation issue with the cluster team — 6/7 historical
   anomaly assemblies and ~44% of ordinary ones have unreadable chunks in the
   original index store.
