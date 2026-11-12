# MooseFS HPRCv2 index-store incident report (for cluster administration)

**Date:** 2026-11-13. **Reporter:** CRISPRapido project (salehi).
**Store:** `mfs#octopus04:9521`, path `/moosefs/salehi/HPRCv2_CRISPRapido/indices/`.
**Impact:** 17 of 28 deep-verified index sets (61%) have persistent unreadable
regions; all affected source assemblies are INTACT. No file has been modified
or deleted by us; all originals are preserved exactly as found.

## Summary

Reading index components returns `OSError [Errno 5] Input/output error` at
specific byte ranges. Errors are **persistent and reproducible**: re-reading
the 4 MiB window at each recorded failing offset fails consistently across
repeated attempts and nodes (octopus06/07/08). They are **not transient**, and
they are **not assembly corruption**: for every affected assembly the frozen
source `*.fa.gz` under `/moosefs/pangenomes/HPRCv2/` decompresses completely
with a valid gzip CRC and full SHA256 (incident-evidence-v1.json).

First observed 2026-10-28 (17 EIO components across 12 index sets during the
index-recovery campaign). A controlled 2026-11-13 batch found the degradation
is far more widespread: **10 of 16 freshly scanned index sets fail**, including
6 of 7 assemblies with historical analysis anomalies.

## Affected files (exact paths, failing components, first failing offsets)

Each row: component file, byte offset where the first read failed (reads use
4 MiB blocks; the true first bad byte is within that block), errno (always 5).
Machine-readable records with per-file recheck evidence:
`00_manifest/incident-evidence-v1/` (copied to
`research/hprcv2/sbds_nrr/results/incident-evidence-v1.json`).

| Assembly | Failing component(s) — offset(s) |
|---|---|
| HG00321_hap2 | `.brt` @ 1,275,068,416; `.txt.bin` @ 2,885,681,152 |
| HG005_pat | `.brt` @ 805,306,368; `.rev.brt` @ 1,409,286,144; `.txt.bin` @ 1,140,850,688 |
| HG02258_mat | `.txt.bin` @ 0 |
| HG02451_pat | `.bwt` @ 738,197,504; `.sa.4` @ 805,306,368; `.txt.bin` @ 134,217,728 |
| HG03195_hap1 | `.bwt` @ 1,006,632,960; `.txt.bin` @ 1,811,939,328 |
| HG03225_hap1 | `.rev.brt` @ 872,415,232 |
| HG03784_hap1 | `.rev.brt` @ 536,870,912; `.sa.4` @ 671,088,640 |
| HG03834_mat | `.bwt` @ 805,306,368; `.rev.brt` @ 1,879,048,192 |
| NA18747_hap2 | `.rev.brt` @ 1,811,939,328; `.txt.bin` @ 1,409,286,144 |
| NA19443_hap1 | `.bwt` @ 671,088,640; `.sa.4` @ 1,073,741,824 |

Historical 2026-10-28 failures additionally affected: HG00128_hap2,
HG01081_pat, HG02630_mat, HG03369_hap2, HG04184_mat, NA18906_mat, NA21309_pat,
CHM13, GRCh38 (all since rebuilt/recovered and verified).

## Reproducer

```python
# any node mounting mfs#octopus04:9521
p = '/moosefs/salehi/HPRCv2_CRISPRapido/indices/HG02258_mat/HG02258_mat.txt.bin'
with open(p, 'rb') as f: f.seek(0); f.read(4*1024*1024)
# -> OSError: [Errno 5] Input/output error
```

## Pattern observations for diagnosis

- Failing components are large binary arrays (`brt`, `rev.brt`, `txt.bin`,
  `bwt`, `sa.4`); tiny components (`meta`, `fsid`, `pos`, `sna`, `cct`,
  `headerSN.bin`) always read fine.
- Failure offsets are multiples of 4 MiB blocks but differ per file/assembly —
  not a single failed disk or chunkserver event; consistent with chunk-level
  damage spread across the store.
- The degradation set grew between 2026-10-28 and 2026-11-13 (new failures in
  files that passed earlier, e.g. HG005_pat) or earlier scans simply did not
  cover these sets; a systematic re-scan would tell.
- Data written to other areas of the same store (run outputs, FASTAs under
  `/moosefs/pangenomes/`) reads back fine so far.

## Request

Please investigate chunk health/replication for the affected files (e.g.
`mfscli` chunk diagnostics), advise whether the damaged chunks are
recoverable from replicas, and whether an fsck/migration is advisable before
we commit to rebuilding ~hundreds of indexes. We will not write anything back
into this index store pending your guidance.
