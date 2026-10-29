# Index recovery and independent auditor report — v1

Started 2026-10-28; finalized 2026-10-29. Branch: `fix/pam-cfd-validation`.

**Index recovery and emitted-row auditor verification passed. The complete search
was NOT launched.** All twelve resolved index sets have matching contig metadata,
complete component checksums and full frozen-FASTA-to-indexed-text verification.
Nine original sets had unreadable components: six were rebuilt and three needed
only text regeneration. Three intact sets were verified without rebuilding.

Application source remains `9c2d3f53673fd135f0b4379812a9a7a4a7f20961`.
Candidate binary SHA256 remains
`f59b1bbf58e70b321359c7bd7b8687aad81bb874b22dc2f083319453726bc876`.
The original production executable remains
`ee543e772ac6e956e8fc3eb88e6ea52517e5192eff9da255902c3bad6073faed`.
No production binary, historical result or original index was changed; no merge
into main or CRISPR-Bulge integration was performed.

## Storage and scope

Original indexes are on MooseFS (`mfs#octopus04:9521`). Available capacity was
about 74.5 TB at report completion; capacity did not prevent unreadable chunks.
The underlying storage cause was not repaired or diagnosed beyond confirmed EIO
reads. Damaged originals are preserved for storage investigation.

Recovery artifacts are on the local `/export` filesystem, `/dev/md0p1`, with
about 11.2 TB free. The recovery directory contains about 130.4 GB, including
source staging, preserved failed attempts, indexes, logs and verification data:

```text
/export/local/home/salehi/crisprapido-index-recovery-v1-2026-10-28
```

[Storage capacity metadata](storage_capacity.json) records exact byte counts.
Large indexes/FASTA/PAF files remain outside Git. This directory is local storage
on the validation host; do not assume its absolute paths are available on every
Slurm worker. Check access or stage checksummed copies before a future job.

## Original index inventory

The original inventory covered twelve required index components plus source
FASTA, FAI and GZI for each assembly: **180 file checks**. Original scans used
one file reader at a time, a 4 MiB block size, and bounded read workers after a
source read stalled. Completed earlier checks were preserved and reused when
the interrupted inventory resumed; they were not replaced with guessed hashes.

- Contig names/order/lengths and index dictionaries matched for all twelve.
- 162 original file checks completed with full SHA256 checksums.
- Seventeen index-component reads returned `errno 5` (Input/output error).
- One HG01081_pat source read timed out and initially had no full checksum.
- A subsequent complete source read passed, SHA256
  `7436cea9f07e120b9efa0838e8f71ba3fab090dfeb05635bd7f9daa991e5d3e7`;
  independent full gzip decompression also passed. The initial timeout remains
  recorded separately and is not described as a sequence mismatch.

[Original file checks](original_file_integrity.tsv) preserve checksums only for
complete reads. [Failure offsets](failure_offsets.tsv) give exact zero-based
read-request starts and requested byte counts, not a claim that every byte in
those ranges is corrupt or that every bad region of each file was enumerated.
A failed file was stopped at its first error, then the next file was checked.
Metadata matching and initial-byte readability were never substituted for these
full-content checks. A targeted retry of a previously failing HG00128_hap2 BRT
range still returned EIO; the source timeout recovery did not make damaged
index components safe to reuse.

## Status of all twelve assemblies

| Assembly | Original unreadable index components | Recovery action | Final result |
|---|---:|---|---|
| HG00128_hap2 | 3 | Full pinned-builder rebuild; successful attempt02 | PASS |
| HG00735_pat | 0 | Original retained; full checks | PASS |
| HG01081_pat | 2 | Full rebuild after successful source retry | PASS |
| HG01496_pat | 0 | Original retained; full checks | PASS |
| HG02055_mat | 0 | Original retained; full checks | PASS |
| HG02630_mat | 3 | Full rebuild | PASS |
| HG03369_hap2 | 3 | Full rebuild | PASS |
| HG04184_mat | 1 | Full rebuild of affected set | PASS |
| NA18906_mat | 1 | Regenerate `.txt.bin`; copy verified other components | PASS |
| NA21309_pat | 1 | Regenerate `.txt.bin`; copy verified other components | PASS |
| CHM13 | 1 | Regenerate `.txt.bin`; copy verified other components | PASS |
| GRCh38 | 2 | Full rebuild | PASS |

Exact original/resolved paths are in [index_status.tsv](index_status.tsv).
**Use [verified_assemblies.tsv](verified_assemblies.tsv) for future execution,
not the original index prefixes for the nine damaged sets.** This is a new
resolved manifest; the earlier frozen manifest is preserved as provenance.
Assembly versions, guide identities and biological sequence content are unchanged.
The nine recovered sets use verified local gzip copies of their frozen FASTAs,
with byte-identical FAI/GZI copies checked against the frozen hashes. The three
intact sets retain their existing original FASTA/index paths.

## Copy search and safe recovery

Existing-copy searches covered the known analysis/software/home locations.
Exact-name candidate search completed without finding a separate healthy frozen
copy; aliases of the same files were not counted as backups. Broader metadata
searches encountered stalls in unrelated trees and were stopped. This is not
an exhaustive claim about every filesystem or every possible renamed copy.

Two additional candidates were explicitly rejected:

- An old HG01081_pat temporary plain FASTA was only 236,150,528 bytes, smaller
  than its 3,034,790,068 required sequence bases. It was truncated and not used.
- A standalone CHM13 benchmark index used `chr1`-style names rather than the frozen
  `CHM13#0#chr1` dictionary. Full version/sequence equivalence was not established;
  it was not silently substituted for the masked-Y/rCRS PanSN derivative.

Six sets with damaged search structures were rebuilt with the existing pinned
Columba VANILLA builder. Its executable SHA256 was
`5987ad90e0c4a061cbf4aac7d6b087a284ccf7fe0f9e1835219c567c5e8d1ad2`;
source checkout commit was `edcc7ce2a1e28ded6784a44afbd8b9348a4805ad`.
All rebuilds used new directories, staged source gzip checksum verification,
full decompression, one thread, default seed length 0 and SA sparseness 4.
At most two full builders ran concurrently on local storage. Each completed
with exit status zero; measured wall times ranged from 11:05.73 to 12:08.76,
with peak RSS about 39.4–40.4 GiB per builder.

The three text-only cases avoided unnecessary suffix-array/BWT rebuilding.
A small exclusive-output helper regenerated the documented `.txt.bin` format
from the verified source: 32-bit length prefix, concatenated uppercased sequence
with deterministic non-ACGT replacement, terminal `$`. Every other component
was copied into the new directory and checksum-matched to its fully readable
original. Recovery-copy/text processing was sequential.

The first HG00128 staging attempt failed before a builder was launched because
`/usr/bin/time` is absent on this host. Its source staging/logs remain preserved;
the successful attempt02 reused those independently checksummed source files.
Use `/home/salehi/miniconda3/bin/time`, not the old plan's `/usr/bin/time` path.
The failed attempt's progress status is not a success marker or a running job.

## Full-content verification

Final verification performed or explicitly reused the complete original scan
for all **144 index components**, with full SHA256s recorded in
[verified_component_checksums.tsv](verified_component_checksums.tsv).
Newly recovered components were fully read and hashed again. All **91** comparisons
between recovered components and their readable original counterparts matched.
No checksum is claimed for the seventeen unreadable original components.

For every assembly, the complete indexed text was compared against the full
frozen FASTA under the pinned Columba normalization, including ambiguous-base
replacement. The source gzip copies for recovered sets checksum-match their
frozen source reads. Contig IDs/order/lengths, offsets, sequence-name arrays,
file-ID arrays and build tag were checked separately. The FAI/GZI copies for
staged sources also match frozen hashes. All twelve passed.

Full verification evidence is in `verification_final/`, with a per-assembly JSON,
comparison stdout/stderr, overall `summary.json` and a PASS marker. Executed
verification source snapshots are archived in `execution_sources/`. The final
repository script also has a tested cardinality guard that prevents an empty
or incomplete panel from falsely claiming twelve assemblies passed, returns a
nonzero exit for blocked verification, and propagates failed rebuild jobs. The actual
panel and the final summary were independently checked to contain exactly twelve
unique assemblies and 144 components.

This establishes complete readability, checksums, frozen-source binding and
recovery consistency. It is not a mathematical proof of every FM-index data
structure, nor a demonstration of genome-wide search sensitivity. No queries
were needed to perform these checks.

## Independent auditors and tests

`scripts/release_validation/audit_emitted.py` is a streaming independent validator
for the corrected Columba/WFA2 PAF contract. It fetches target/PAM sequences from
indexed FASTA with samtools rather than using Rust-extracted fields. It checks:

- Complete canonical NGG PAM and eligibility on both strands; no truncated or
  ambiguous PAM acceptance.
- Guide-oriented negative CIGAR replay, with target/PAM reverse-complemented once.
- Exact match/mismatch assertions, RNA/DNA gaps, query/reference consumption,
  full-guide coordinates, contig lengths and boundary conditions.
- Mismatch/bulge limits, minimum match fraction, match/block counts and tags.
- Independent CFD aligned-input reconstruction, existing DNA-bulge projection,
  matrix multiplication and four-decimal output tolerance.
- Missing/unusable scores, duplicate tags, invalid operations and malformed
  coordinates. Failures retain diagnostics and cannot produce an audit PASS.

**Twenty-three unit/integration tests passed**, including deliberate bad PAM, CFD,
orientation, CIGAR and coordinate cases, known numerical fixtures, terminal
insertions/deletions, contig edges, recovery overwrite guards, truncated index
text, and an empty-panel rejection. Small real-sequence fixtures cover both
strands and mismatch/RNA/DNA-bulge classes.

The auditor rechecked the prior strict 88-row real subset with **zero errors**:
49 positive strand, 39 negative strand, 72 bulged rows and 16 ungapped rows.
It also checked these same 88 alignments at their original genomic coordinates
against the staged full HG00128_hap2 indexed FASTA, again with zero errors. The
coordinate lift used the original saved source-window offsets, not a new search.
All 88 CFD reconstructions agreed; 50 scores were >=0.20 and 38 were below.
The prior permissive-output negative control correctly rejected all 32 non-NGG
rows while validating its 88 canonical rows. That expected FAIL is a successful
negative-control test, not an unresolved discrepancy.

These checks reuse previously validated examples; they are not new complete-search
results or a revalidation of all natural loci. A row-soundness auditor cannot
prove that an empty/missing output contains every possible candidate. Per-guide
zero counts in these subset tests are not genome-wide absence claims.

## Reproduction and metadata

Detailed paths, hashes and auditor summaries are in
[recovery_metadata.json](recovery_metadata.json), with final script hashes in
[script_hashes.json](script_hashes.json). The large artifact directory retains:

- `integrity_original/` and `integrity_resumed/`: sequential scan evidence, exact
  read errors and preserved checksums.
- Copy-search records, source retry/decompression results, original interrupted
  read evidence and rejected candidate details.
- Each `rebuilt_*/` or `recovered_text_*/` directory: source staging, command,
  executable hashes, construction logs, time/RSS measurements and status.
- `verification_initial/` and `verification_final/`: separate preserved passes.
- Auditor test logs, record-level TSVs and errors, including expected negative
  controls. No prior result directory was reused or overwritten.

Run small tests without index scans:

```bash
python3 -m unittest discover -s scripts/release_validation -p 'test_*.py' -v
```

Representative independent emitted-row audit (choose a NEW output directory):

```bash
python3 scripts/release_validation/audit_emitted.py \
  --paf /path/to/corrected.paf --fasta /exact/verified/assembly.fa.gz \
  --guides docs/release_candidate_2026-10-28/guides_100.tsv \
  --matrices /path/to/pinned/score-matrices \
  --samtools /moosefs/guarracino/condatools/samtools/1.22.1/bin/samtools \
  --output /new/audit-output
```

To reproduce the expensive full-content pass, use the original frozen manifest
plus the resolved override/evidence files and a NEW output directory:

```bash
ROOT=/export/local/home/salehi/crisprapido-index-recovery-v1-2026-10-28
python3 scripts/release_validation/verify_recovery_set.py \
  --manifest docs/release_candidate_2026-10-28/assemblies_12.tsv \
  --integrity-summary "$ROOT/integrity_with_source_retry.json" \
  --overrides "$ROOT/overrides_final.json" \
  --compare-bin "$ROOT/tools/compare_index_text" \
  --output /new/full-verification-output
```

Individual build/status files record exact source, builder, loader, parameters
and output prefix. Recovery helpers refuse to reuse an output directory. They do
not invoke CRISPRapido or Columba queries.

## Remaining risks and complete-search resource estimate

No index-recovery blocker remains for the resolved twelve-set manifest. The
original MooseFS faults remain unresolved, and the nine original damaged sets
must not be used. Three intact sets still reside on MooseFS; their complete
checks passed, but verify input access/hashes again at execution time. Local
recovery paths must be accessible or staged with checksum verification on the
chosen execution host.

The **100-guide complete-search validation is still unrun**. The existing plan's
independent exact-match recall and bounded exhaustive mismatch/bulge coverage
oracles are separate from this emitted-row soundness auditor; global sensitivity
has not been established. Complete those coverage checks as required by the
reviewed execution plan, and do not describe index recovery as a release or
CRISPR-Bulge integration pass. Main's previously documented merge conflicts also
remain unchanged; no merge was attempted.

For the future strict-NGG run, retain the frozen 100 guides and parameters
`-p GG -m 4 -b 1 -z 1 -f 0.75 -t 4`. Planning requests remain six allocated CPUs,
64 GiB RAM and two hours per assembly, first one pilot, then at most two
concurrent jobs. Budget 15–60 minutes/search (3–12 job-hours for twelve), plus
independent auditing; these remain estimates, not measured search timings.
Use bounded audit batches (default 500), reserve at least 50 GiB local scratch
per search when staging, and measure pilot output/log growth before expansion.

**Stopped here as requested: no complete search, no merge, no production binary
replacement, no historical-index/result overwrite, and no CRISPR-Bulge work.**
