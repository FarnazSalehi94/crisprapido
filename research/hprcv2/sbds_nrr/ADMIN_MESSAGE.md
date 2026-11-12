# Message to cluster administrators — MooseFS index corruption

**To:** cluster-admin. **From:** salehi (CRISPRapido project).
**Subject:** Persistent EIO on MooseFS index store — request chunk/replica investigation.

Hello,

we need your help with persistent read failures on the MooseFS store
(`mfs#octopus04:9521`) mounted at
`/moosefs/salehi/HPRCv2_CRISPRapido/indices/`.

**What we see:** reading certain large index components fails with
`OSError [Errno 5] Input/output error` at specific byte ranges. The failures
are **persistent and reproducible**: re-reading the same 4 MiB window fails
consistently across repeated attempts and across nodes (octopus06/07/08).
Example reproducer on any node:

```python
p='/moosefs/salehi/HPRCv2_CRISPRapido/indices/HG02258_mat/HG02258_mat.txt.bin'
open(p,'rb').read(4*1024*1024)   # -> OSError: [Errno 5] Input/output error
```

**Scale:** we have so far fully scanned 28 index sets (12 in late October,
16 more this week): **17 have unreadable regions**. Affected components are
the large arrays (`brt`, `rev.brt`, `txt.bin`, `bwt`, `sa.4`); small metadata
components always read fine. Failing offsets differ per file — consistent
with damaged chunks scattered across the store rather than one bad disk.
The full evidence (paths, components, first failing offsets, per-file
re-read confirmations) is in
`research/hprcv2/sbds_nrr/results/incident-evidence-v1.json` and
[INCIDENT_REPORT.md](INCIDENT_REPORT.md) in the project repository.

**Importantly, the biological data is fine:** every affected assembly's
source FASTA (under `/moosefs/pangenomes/HPRCv2/`) decompresses completely
with a valid gzip CRC and full SHA256. Only the derived index files on the
CRISPRapido store are damaged.

**Requests:**
1. Please check chunk health/replication for the affected files (e.g. via
   `mfscli`/master metadata) and tell us whether the damaged chunks can be
   recovered from replicas or backups.
2. Let us know if an fsck/migration is advisable — we have stopped writing
   to this store pending your guidance.
3. If recovery is possible, we would much prefer it over rebuilding: a
   mass rebuild would cost us roughly 65–90 one-CPU hours and ~2.8 TB of new
   storage, and rebuilding back onto the same store would risk recurrence
   if the root cause remains.

We have not modified or deleted anything on the damaged store; all original
files are preserved exactly as found. Happy to provide any further scans or
details. Thank you!
