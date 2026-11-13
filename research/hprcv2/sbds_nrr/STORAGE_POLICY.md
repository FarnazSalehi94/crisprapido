# MooseFS storage policy compliance and reliability evidence (2026-11-13)

Implements the user's corrected storage directive: **all new genomic data,
rebuilt indexes, intermediates, logs and results live under
`/moosefs/salehi/HPRCv2_CRISPRapido/`**; `/export` hosts no new persistent
genomic datasets. Rebuilt indexes: `indices_rebuilt_v1/`; SBDS results:
`runs/sbds_ngg_nrr_v1/`; existing `indices/` untouched.

## New permanent layout

```text
/moosefs/salehi/HPRCv2_CRISPRapido/
  indices/              original (damaged) indexes — never modified or deleted
  indices_rebuilt_v1/   rebuilt/trusted sets, one dir per assembly:
                        <asm>.{12 components}, build-host.txt, completion.json
  runs/sbds_ngg_nrr_v1/ all campaign tooling, bindings, results, logs
```

Rebuilds build on node-local `/scratch` (fast; no partial MooseFS exposure),
verify every component on scratch, copy to MooseFS with per-file fsync,
re-read every byte from MooseFS (write-back hash comparison), then a
**separate job on a different worker** re-reads all components and compares
against the completion manifest before a set is `TRUSTED_CROSS_NODE`. Only
then may binding/scoring use it. Checksummed completion manifests record
source/builder/component hashes and timings. Confirmed-corrupted originals
are **never deleted by tooling**; deletion requires a reviewed manifest and
is not implemented at all yet (requirement 7). No healthy index, source
FASTA or historical result is ever touched (requirement 8).

## Reliability evidence gathered this session

| Test | Data | Result |
|---|---|---|
| Synthetic write test | 2 GB written, read back on octopus06 | checksum match |
| Pilot-index migration | 3 × ~11.5 GB copied to MooseFS, fsync + write-back verified | all match; cross-node recheck PASS (octopus06) |
| First rebuild batch | 8 × ~11.5 GB built on scratch, written to MooseFS, write-back verified | all match |
| Cross-node trust gate | all 11 sets re-read on a different worker (octopus06/07) | 11/11 TRUSTED_CROSS_NODE |
| End-to-end relocation check | HG03784_hap1 search rerun after index relocation | results **identical** to pre-relocation run (192,353 hits; 42,415 CB scored; 149,938 N/A) |

~137 GB of new MooseFS writes verified from independent workers with zero
mismatches, plus one full search reproduction. **Empirically, new writes
under the store verify reliably so far** — but the corruption that hit
`indices/` appeared over time, so this cannot prove long-term safety
(user requirement 2 acknowledged). The administrator consultation
([ADMIN_MESSAGE.md](ADMIN_MESSAGE.md), updated) now asks explicitly whether
the filesystem is healthy enough for mass rebuilding and what chunk-level
diagnostics show for the damaged files.

## Migrated pilot indexes

HG005_pat, HG03225_hap1, HG03784_hap1 rebuilt during the /export pilot were
migrated to `indices_rebuilt_v1/` with hash-identical verification, their
bindings regenerated against the new locations, and HG03784_hap1's four-PAM
search rerun with identical results (its old SUCCESS marker is preserved as
`SUPERSEDED-index-relocated-to-moosfe.json`). The /export pilot copies are
preserved untouched pending a deletion review; no new genomic data was
created under /export.

## Ledger: authoritative-manifest-v4.json

**144 READY_BOUND** (126 + 10 newly bound trusted + 8 rebuilt-batch-1),
**318 INDEX_EIO**, **4 MISSING_INDEX** — exactly 466, mutually exclusive.
All 12 trusted-recovery assemblies now have committed bindings. The 124-
assembly incremental search plan is prepared (`sbds-plan-relocation-verify-
v1.json`) and **deliberately not launched**.

## Remaining scale and what waits on whom

- 318 EIO + 4 missing sets need rebuilds (~110 one-CPU hours, ~3.8 TB under
  the measured model; `indices_rebuilt_v1` currently uses ~130 GB).
- **Mass rebuilding is deliberately held** pending the administrators'
  answer on MooseFS health (user requirement: confirm reliability before
  large-scale rebuilding). The 8-assembly batch above is the bounded
  reliability pilot; results support proceeding once admins confirm the
  root cause does not affect new writes.
- The /export pilot copy deletion and the corrupted-original deletion both
  require a reviewed deletion manifest — not implemented, awaiting review.
