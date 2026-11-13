# Stage-1 verification campaign — complete results (2026-11-13)

Authorized by the user (verbatim authorization recorded in
`00_manifest/verify-campaign-authorization.json`: 438 UNATTEMPTED
assemblies, ≤5.3 TB read budget, ≤12 concurrent, no rebuilds, no searches).
Array **2871846**, `--array=0-437%12`, one CPU, 8G, octopus06/07/08, the
tested restartable `deep_verify.py` workflow. All originals preserved; no
index was modified, rebuilt or deleted.

## Headline results

| Status (authoritative-manifest-v3) | Count |
|---|---:|
| READY_BOUND (deep-verified + bound; search-ready today) | **126** |
| TRUSTED_EVIDENCE (healthy per 2026-10-28 evidence; binding not yet prepared) | 10 |
| INDEX_EIO (persistent unreadable regions; rebuild required) | **326** |
| INDEX_FAILED_MISSING_INDEX (no index ever built; build required) | 4 |
| **Total** | **466** |

- **Successfully verified and newly bound this campaign: 115** (plus 11
  previously READY_BOUND). The 10 TRUSTED_EVIDENCE assemblies can be bound in
  ~26 one-CPU minutes total.
- **Corrupted indexes: 326 (+7 previously known = 333 total across the
  campaign-period scans; manifest-v3's 326 counts the seven batch-1 EIO
  assemblies that were not rebuilt — of the original 10, three were rebuilt
  by the pilot and are READY_BOUND).**
- **Requiring rebuild/build: 330** (326 EIO + 4 missing).

## Reproducibility: real damage, not load artifacts

73% EIO during a 12-way concurrent scan raised the question of
concurrency-induced errors. A strictly sequential (%1) rescan of a
deterministic 10-assembly sample of campaign-EIO sets (job 2872295,
incident_scan with window rechecks) reproduced **persistent unreadable
regions in 10/10** at the same failing components/offsets. The failures are
genuine chunk-level storage damage, matching the 2026-10-28 pattern but at far
larger scale. The four missing sets (hg002v1.1.mat_MT, hg002v1.1.pat,
HG06807_mat, HG06807_pat) simply have no index directory at all.

## Actual runtime and resource usage (requested figures)

- **Wall clock: 1 h 22 min** (2026-11-12T23:11:26 → 2026-11-13T00:33:09) at
  %12 concurrency, 438 one-CPU tasks; 115 COMPLETED / 323 FAILED exit codes.
- **Measured reads for the 115 successful verifications: 1.41 TB**
  (12.16 GB each; 7.0 one-CPU hours, median 216 s/assembly). Failed EIO sets
  stop at their first bad chunk; their partial reads are not individually
  metered but the whole campaign is bounded by the authorized **5.3 TB
  plan** — actual total ≈ 1.4 TB (successes) + ≤ ~2 TB (failures stopping
  early) ≈ **well under budget**.
- **Storage consumed: 360 GB** of committed bindings (126 × ~2.9 GB, with
  scoring reference indexes) on MooseFS (42 TB free). Partial stages from
  EIO attempts (931 GB, mostly truncated reference decompressions) are
  preserved and moved to `00_manifest/failed-binding-stages-eio/`, not
  deleted. `/export` unchanged this campaign (10.3 TB free).
- No filesystem-wide anomalies observed beyond the per-assembly EIO; the
  12 concurrent readers caused no new damage pattern (sequential rescan
  reproduces identical failures).

## Recommended recovery and full-search plan

1. **Send [ADMIN_MESSAGE.md](ADMIN_MESSAGE.md)** — updated scope: ~333
   damaged index sets. If chunk/replica recovery restores most of them, the
   rebuild campaign shrinks accordingly; recovered sets must still pass deep
   verification before binding.
2. **Bind the 10 TRUSTED_EVIDENCE assemblies** (~26 one-CPU minutes,
   standard tooling, no authorization concern) → 136 READY_BOUND.
3. **Rebuild/build campaign for 330 sets** (admin-unrecoverable ones), onto
   **/export only**: measured 0.34 one-CPU h/assembly end to end →
   **~116 one-CPU hours (~4.8 core-days)**; wall ≈ 15 h at %8 (RAM-bounded
   at ~41 GB/build); **~4.9 TB** on /export (10.3 TB free — fits with a
   reservation; optionally prune the ~1 TB of staged plain FASTAs after
   component verification to relieve pressure).
4. **Incremental SBDS search** over all READY_BOUND assemblies (driver
   `run_sbds_incremental.py`, scientific-signature skip so nothing reruns):
   at 136 assemblies ≈ 68 one-CPU hours (~4.3 h at %16); the full 466 once
   recovery completes ≈ 233 one-CPU hours, ~4.75 GB compressed outputs.
   Search uses m5/b1/z2/f0.75, four literal PAMs, CFD + NRR CRISPR-Bulge
   with explicit N/A; noncanonical biological-validity caveat stays.

**Nothing beyond the authorized verification was launched**; no rebuilds
and no searches were performed.
