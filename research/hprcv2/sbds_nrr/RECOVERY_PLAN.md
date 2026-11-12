# SBDS recovery strategy and campaign preparation (not launched)

Baseline `07d3581`; originals and historical results unchanged; nothing in
this pass was launched as a campaign. Root:
`/moosefs/salehi/HPRCv2_CRISPRapido/runs/sbds_ngg_nrr_v1/`.

## 1. Corrected ledger (Task 1)

The previously reported "20 search-ready + 10 EIO + 438 = 468" double-counted
HG02055_mat and NA18906_mat (both in the trusted-12 AND the committed-bindings
8). The authoritative manifest
(`00_manifest/authoritative-manifest-v2.json`, tool `reconcile.py`) contains
**exactly 466 unique assemblies with mutually exclusive statuses**:

| Status | Count | Meaning |
|---|---:|---|
| READY_BOUND | 11 | committed binding.json (search-ready today) |
| TRUSTED_EVIDENCE | 10 | trusted 2026-10-28 evidence, binding not yet prepared |
| INDEX_EIO | 7 | deep verification failed: persistent unreadable index regions |
| UNATTEMPTED | 438 | no full-content evidence either way |

(11+10+7+438 = 466; three formerly-EIO assemblies became READY_BOUND via the
rebuild pilot below, leaving 7 EIO.)

## 2. Storage-failure investigation (Task 2)

Full per-component scans of all ten failed sets with reproducible re-reads
(`incident_scan.py`, evidence committed): **all ten are persistent unreadable
regions (errno 5), reproducible across attempts and nodes — not transient
errors; not missing components; and NOT corrupted assemblies**: every frozen
source FASTA.gz decompresses with valid CRC and full SHA256. This is
distributed chunk-level damage in the MooseFS index store, growing since
2026-10-28. A concise administrator-facing report with exact paths, failing
components, offsets and a reproducer is [INCIDENT_REPORT.md](INCIDENT_REPORT.md).
Failed originals are preserved untouched; nothing was rebuilt onto MooseFS.

## 3. Controlled rebuild pilot (Task 3) — measured, end to end

Three failed assemblies rebuilt from their verified frozen FASTAs onto
**healthy local /export RAID** (`rebuild_pilot.py`, jobs 2871822, exit 0;
preserved `columba_build` binary + guix loader from the attested recovery,
original indexes untouched):

| Assembly | Build wall | Peak RSS | Index bytes | Decompress | Verify |
|---|---:|---:|---:|---:|---:|
| HG005_pat | 15:54.6 | 40.15 GB | 11.10 GB | 30 s | 8 s |
| HG03225_hap1 | 16:54.7 | 40.23 GB | 11.12 GB | 46 s | 8 s |
| HG03784_hap1 | 17:10.5 | 41.75 GB | 11.55 GB | 46 s | 9 s |

All 12 components of each rebuilt set were fully read and hashed. All three
were bound via the standard deep verification helper (jobs 2871828–30), and
**HG03784_hap1 then passed a complete four-PAM verification search**
(job 2871832, 30 min, 192,353 hits, 100 % CFD-defined, 42,415 CRISPR-Bulge
defined, explicit N/A elsewhere, all four audits PASS) — proving the full
rebuild → binding → search → audit → scoring chain. The rebuilt indexes live
on /export (`sbds_rebuilt_indexes_v1/`), **not** on the degraded store.

Per-assembly rebuild cost model (measured): ~17–19 min wall (incl. decompress
and verification), ~41 GB peak RSS (limits concurrency to ~4–8 per
high-memory node), ~14.7 GB new storage (3.05 GB plain FASTA + 11.4 GB index
+ binding on MooseFS 3.06 GB). Six historical rebuilds (11–12 min, 41–42 GB)
agree.

## 4. Campaign preparation (Task 4) — ready, awaiting explicit authorization

**Verification campaign (438 UNATTEMPTED, existing tested tooling):**
`submit_deepverify.slurm` arrays over `authoritative-manifest-v2.json` rows
with status UNATTEMPTED. Requested I/O budget: **≤ 438 × 12.16 GB ≈ 5.3 TB
reads** (EIO sets stop early), ~438 × ~155 s ≈ **19 one-CPU hours** (~2.6 h
wall at %12), committed bindings up to ~1.4 TB on MooseFS (45 TB free).
Conservative settings: `%4–%12`, 1 CPU, 8G, nodes with scratch headroom
(octopus05 excluded). Every failure is preserved and separately labelled;
metadata presence authorizes nothing; restart rechecks completed work.

**Rebuild campaign (conditional on verification outcome + cluster-team
advice):** if EIO prevalence stays ≈44% of ordinary assemblies, roughly
**~190 of 438** will need rebuilds: **~55–60 rebuild-days of one-CPU wall
(~19–20 h measured per batch of 60 at %4 on 41 GB-RAM nodes)**, ~2.9 TB of
rebuilt indexes + ~0.6 TB staged FASTAs on /export (11 TB free — fits, but
needs a reservation), plus ~0.6 TB of bindings on MooseFS. RAM-bound
concurrency is the limiter. **Do not rebuild onto MooseFS before the
incident is diagnosed**; keep rebuilt sets on /export as the pilot did.

**Recommended order:** (1) administrator diagnosis of the MooseFS damage
(INCIDENT_REPORT.md) — replica recovery could avoid most rebuilds;
(2) authorized verification campaign (surfaces the true EIO count);
(3) authorization decision for rebuild scale with the measured cost model;
(4) only then the four-PAM SBDS search over all READY_BOUND assemblies,
expanding as sets become ready. Historical outputs remain untouched; missing
inputs are never zero-filled.
