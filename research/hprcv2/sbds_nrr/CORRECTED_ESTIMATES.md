# Corrected SBDS resource estimates (2026-11-13)

Baseline `ff30795`. All figures derive from measured pilots on this cluster
(rebuild pilot 2871822, verification batch 2871790, cohort test 2871778/2871780,
NRR pilot 2871759). **No campaign launched.**

## 1. Rebuild costs, recomputed (Task 1)

Measured per assembly: decompress+verify FASTA 30–46 s; columba_build
**955–1,032 s (15:54–17:10)**; component verification 8–9 s; binding
preparation ~155 s. End-to-end ≈ **0.34 one-CPU hours**.

| Rebuild count | One-CPU hours | Wall at %4 | Wall at %8 | Wall at %16 | New storage |
|---:|---:|---:|---:|---:|---:|
| 150 | 52 | 13 h | 6.5 h | 3.2 h | ~2.2 TB |
| 190 | 65 | 17 h | 8 h | 4 h | ~2.8 TB |
| 270 (all-EIO worst case) | 93 | 23 h | 12 h | 6 h | ~4.0 TB |

**Where "55–60 rebuild-days" came from:** it was ~54 **serial one-CPU hours**
(190 × ~17 min ≈ 54 h) mislabeled as "days", then wrongly divided again by
concurrency. Correct serial figure: ~2.2 core-days; at any concurrency ≥4 it
is under a day of wall time. The real constraint is **RAM** (≈41 GB peak per
build), which caps simultaneous builds per node (~4–8 on the 386–516 GB
nodes); Slurm's memory accounting handles this automatically at %8–%16.

Component breakdown (measured): build 93 % of end-to-end time; verification
(component hashing) <1 %; staging/decompression ~4 %; binding preparation
~3 %. Verification campaign per assembly adds 12.16 GB reads / ~155 s.

## 2. Verification campaign (438 UNATTEMPTED) — prepared, NOT launched

- Budget requested: **≤ 5.3 TB reads** (438 × 12.16 GB; EIO sets stop early),
  **~19 one-CPU hours** (~1.6 h wall at %12; realistically 2–4 h with
  cold-cache variance), 1 CPU / 8G per task, up to **~1.4 TB bindings** on
  MooseFS (45 TB free).
- Assembly list frozen: `00_manifest/verify-campaign-unattempted-438.txt`;
  submission validated with `sbatch --test-only` (job would start on
  octopus06–08). **Awaiting explicit user authorization of the ~5.3 TB read
  budget.** Restartable; committed assemblies recheck and skip; failures
  labelled EIO/other with originals untouched.

## 3. SBDS four-PAM search — incremental driver prepared, NOT launched

`run_sbds_incremental.py` (committed + staged) reads the authoritative
manifest, selects READY_BOUND-and-not-COMPLETE assemblies, freezes a fresh
plan and prints the sbatch command; `--submit` additionally requires a
verbatim authorization string and refuses concurrency >32. The worker now
records per-PAM **scientific search signatures** in SUCCESS, so completed
assemblies are **never rerun** even under later expanded plans — they are
checksum-verified and skipped; only a changed scientific identity (guide,
parameters, reference, tooling) is refused. Seven regression tests cover
restart, tamper, timeout, exclusivity and the new expansion semantics.

Current dry run: **11 READY_BOUND, 8 pending** (HG00097_hap1, HG005_pat,
HG01258_mat, HG02027_mat, HG02984_pat, HG03225_hap1, HG03710_pat,
NA18983_hap1; HG02055_mat/NA18906_mat/HG03784_hap1 already complete),
prepared plan `sbds-plan-dryrun-v1.json`
(SHA256 `966f4cb9291646e74810ed396a36108b67330f62a5e251922d0128d6202cb5d9`).
Search cost model: ~30 min/assembly one CPU (measured 22–30 min + audit/
scoring), 56G per job, ~233 one-CPU hours for all 466; ~14.6 h wall at %16.
Storage ~4.75 GB compressed outputs total.

## 4. Brunello separation

Brunello (77,441 guides, provisional m4/b1/z1/GG) remains a separate project
with its own reports; nothing Brunello was touched or launched here.

## 5. Next actions (in order)

1. Send [ADMIN_MESSAGE.md](ADMIN_MESSAGE.md); wait for chunk/replica
   diagnosis — replica recovery would avoid most rebuilds.
2. User authorizes the **verification campaign I/O (~5.3 TB)** → launch
   `submit_deepverify.slurm` over the 438 list at %12.
3. Verification outcome sets the true rebuild count; if needed, user
   authorizes the rebuild campaign (≤~93 one-CPU hours, ≤~4 TB on /export,
   onto /export only).
4. User authorizes incremental SBDS search via the driver (`--submit` with
   verbatim authorization); the driver is re-run after each verification/
   rebuild wave to admit newly READY_BOUND assemblies without rerunning
   completed ones.
