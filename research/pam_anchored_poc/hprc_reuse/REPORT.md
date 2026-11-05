# Existing HPRC Columba indexes and experimental-backend reuse

2026-11-05. **The existing HPRC indexes have the same observed basic Columba format as the prototype indexes, but the experimental backend cannot currently use them directly.** The blocking interface requirements are the missing experimental reference/index binding manifest and, for the references listed in the resolved manifest, compressed FASTA input. No HPRC index was rebuilt. A metadata and provenance audit is not an index search or a proof that every file remains readable today.

The smallest promising reuse path is an **independently verified import of an existing index**, not a rebuild: verify the complete current reference and component bytes against trusted evidence, verify their text/contig relationship including original Ns, then create the prototype's binding sidecar in a separate staged prefix and use the unchanged checked snapshot path. This import path does not exist today and must be implemented and tested before claiming the existing indexes are usable. It must not accept metadata alone as proof or weaken the current verification. Whole-assembly resource use and large-coordinate correctness also remain untested for this backend.

## Location and inventory

The panel source list is `/moosefs/salehi/HPRCv2_CRISPRapido/hprc_haplotypes.tsv` (464 rows, SHA256 `9147e9c96cc39cb42de2c9420523e88aedd1ab17d7514ce5ee90e2210fdc5fb9`). Original index directories are under:

```
/moosefs/salehi/HPRCv2_CRISPRapido/indices/<index-directory>/<index-prefix>.*
```

That root currently has 466 directories, including the two reference assemblies; directory count is not a health check. The read-only [metadata inventory](inspect.py) followed the 464 panel IDs and attempted to read `.meta` and stat all 12 expected components. It observed nonempty components and `21 / 4 / VANILLA` metadata for **459** direct name matches. Four additional IDs use different on-disk directory/prefix names, and [separate alias checks](inspect_aliases.py) found their 12 nonempty components and the same metadata; their FASTA `.fai` dictionaries and index header/offset arrays also matched:

| Panel ID | Directory and prefix name |
|---|---|
| `hg002v1.1.mat_MT` | `hg002v1.1.mat_MT.PanSN` |
| `hg002v1.1.pat` | `hg002v1.1.pat.PanSN` |
| `HG06807_mat` | `HG06807_mat_v1.0.genbank` |
| `HG06807_pat` | `HG06807_pat_v1.0.genbank` |

Thus **463 of 464** panel entries had the expected component names, nonzero recorded sizes and metadata observed via their respective paths. The remaining `HG02723_mat` returned an **I/O error** during metadata inspection; a later attempt to read its `.meta` also stalled and was interrupted. This does not establish which other components are damaged, or that every byte of the other 463 is readable. The full [panel metadata](results/panel_metadata.json), [alias metadata](results/alias_metadata.json), and [initial inventory summary](results/summary.json) record what was checked. The summary's five initial failures consist of these four name mismatches plus the I/O error; the aliases were checked separately. No binding manifest was found among the 463 successfully observed original prefixes. The failed entry's binding status is unknown.

### Previously resolved 12-assembly subset

The earlier [recovery report](../../../docs/index_recovery_2026-10-29_v1/RECOVERY_REPORT.md) records 10 HPRC haplotypes plus CHM13 and GRCh38. Of the 12 original sets, six were rebuilt previously in separate local directories, three had `.txt.bin` regenerated previously with other verified components copied, and three originals were retained. Those actions happened before this task; **no reconstruction or rebuilding was performed now**. [verified_assemblies.tsv](../../../docs/index_recovery_2026-10-29_v1/verified_assemblies.tsv) gives the resolved prefix/reference per assembly; do not assume a damaged original prefix is equivalent to its recovered prefix.

Examples:

* Original retained `HG00735_pat`: `/moosefs/salehi/HPRCv2_CRISPRapido/indices/HG00735_pat/HG00735_pat`.
* Previously recovered `HG00128_hap2`: `/export/local/home/salehi/crisprapido-index-recovery-v1-2026-10-28/rebuilt_HG00128_hap2_attempt02/index/HG00128_hap2`.

For **all 12 resolved prefixes**, the current read-only metadata pass found matching index header names/order/lengths against the listed FASTA `.fai`, matching `.pos` and `.sna`, correct text length prefix/file length/sentinel, `21 / 4 / VANILLA`, and all 12 component paths present. All 12 listed references are `.fa.gz`; none of those prefixes has the experimental binding sidecar. See [resolved_metadata.json](results/resolved_metadata.json). Only the text length prefix and final byte were read from `.txt.bin`, not its contents, and component bodies were not hashed in this pass. This is a structural observation, not a current whole-file identity/readability certification.

The historical recovery evidence reports a complete pinned-normalization FASTA-to-index-text comparison, complete source and component hashes for all 12 resolved sets, and explicit comparisons with readable originals. The later [whole-assembly pilot](../../../docs/whole_assembly_pilot_0aa9944_2026-11-02_v1/REPORT.md) reports a successful full checksum/size scan of these inputs on 2026-11-02. The historical [component manifest](../../../docs/index_recovery_2026-10-29_v1/verified_component_checksums.tsv) has SHA256 `27d1db97d287b40f9ea9d7f690c38404462d411beaa6c7cdc7bffb429b46fffc`; the [resolved assembly manifest](../../../docs/index_recovery_2026-10-29_v1/verified_assemblies.tsv) has SHA256 `948e0c17f144bf4b30624fba62f593816472463e2e1db8f755b3be611a62b591`. Those are valuable provenance, but do not prove the files have not changed or become unreadable since those scans. They cover the 12 resolved assemblies, not the whole 464-entry panel.

## Format and coordinate comparison

The same 12 suffixes occur in the existing indexes and the prototype's tested vanilla indexes:

```
.bwt .brt .rev.brt .sa.4 .sa.bv.4 .pos .sna .txt.bin
.cct .meta .fsid .headerSN.bin
```

Both record build tag **21**, **4-byte** `length_t`, and **VANILLA** in `.meta`. The `.sa.4` components and prototype constructor use SA sampling four. The native build uses alphabet 5 and `-DTHIRTY_TWO`, linked against the recorded local Columba source checkout `edcc7ce2a1e28ded6784a44afbd8b9348a4805ad`. The untracked runtime build is not independently attested to the source. These observations support *format compatibility*, subject to full file integrity and version/ABI checks, but do not demonstrate an actual HPRC prototype search.

As a concrete comparison, `HG00735_pat` contains **3,029,527,752** reference bases and **82** contigs. Its 12 components total **11,455,409,905 bytes**: `.txt.bin` is 3,029,527,757, `.sa.4` is 3,029,527,756, `.bwt` 1,136,072,928, `.brt` and `.rev.brt` 1,893,454,928 each, and `.sa.bv.4` 473,363,736 bytes; the six remaining metadata components account for the rest. The experimental CHM13 3 Mb index has the same component set and metadata, plus its extra `.crispr-binding.json` sidecar; its 12 index components total 11,345,125 bytes and its sidecar 1,577 bytes. `HG00735_pat`'s `.pos` is 332 bytes, exactly 83 32-bit offsets for 82 contigs plus the final end.

In the current prototype, `idx.findSA` returns an offset in the concatenated indexed text. The prototype uses the supplied FASTA record order and cumulative lengths to locate the contig, rejects spans that cross a contig boundary, and emits zero-based half-open contig coordinates; a plus-strand PAM anchor is the target end and a minus-strand anchor is the target start. Current checks on the 12 resolved indexes establish that their names/order/lengths and index offsets agree **at the metadata level**. They do not validate actual large-index `findSA` behavior, PAM anchors, both strands or boundary arithmetic in this backend. PanSN names such as `HG00735#1#CM087877.1` are retained as contig IDs; any result must be audited against the same exact assembly, not mapped to a similarly named reference by assumption.

## Why direct invocation fails

The integrated [snapshot runner](../integration/run_snapshot.py) expects `PREFIX.crispr-binding.json` and copies it alongside all components. The [native prototype](../poc.cpp) independently invokes the mandatory [index_manifest.py](../compatibility/index_manifest.py) verifier before opening the FM index. The helper requires the `crispr-poc-index-binding-v1` sidecar with SHA256 of the **exact plain FASTA bytes** and all 12 index components; its `build` mode refuses an already-existing prefix and creates a binding as part of a fresh build. Existing HPRC component manifests have different provenance/schema and are not automatically accepted as that sidecar. Do not use the absence of a sidecar to conclude the underlying FM data must be rebuilt.

The native and Python reference parsers also require a **plain** FASTA, with A/C/G/T/N reference symbols and supported naming rules. The resolved manifest points to BGZF/`.fa.gz` inputs. Some previous local staging directories already have plain FASTAs, including `rebuilt_HG00128_hap2/input/HG00128_hap2.fa` (3,075,802,891 bytes) and `rebuilt_HG01081_pat/input/HG01081_pat.fa`. Their current exact identity to their frozen compressed sources was not rechecked here; verify it before reuse. An existing plain staged FASTA may avoid decompression, but neither its existence nor matching metadata is a binding proof. Other input symbols and duplicate IDs must still be rejected or explicitly supported without changing N/PAM semantics.

A deliberately limited **preflight-only** direct-native invocation used the existing `HG00735_pat` prefix, its listed gzip FASTA and one guide, with zero edits. It exited 1 with no stdout before opening the FM index:

```
index binding error: index binding manifest missing; rebuild with index_manifest.py build
mandatory reference/index binding verification failed
```

[direct_gate.json](results/direct_gate.json) records the invocation. The current error text suggests rebuilding because the helper only implements build-time binding; **this investigation does not recommend or authorize rebuilding the HPRC indexes**. The test proves the fail-closed gate, not that the FM-index can be loaded or searched. No actual query was executed. Even a single-guide invocation of the current prototype searches the entire supplied indexed assembly; it has no genomic-window search switch. Running it on an HPRC prefix would be a whole-assembly search and is outside this task's constraint.

## Minimal safe reuse design

1. Add a separate importer for an **existing** supported 32-bit vanilla index. It must require an explicitly identified reference and trusted build/recovery evidence, not infer pairing from filenames or `.fai` alone. Confirm metadata/width/flavor/SA sampling and exact contig names, order, lengths and offsets; reject unsupported symbols or formats clearly.
2. Revalidate **current full bytes**, with bounded streaming and timeout/error handling: hash the original compressed source against a trusted pinned digest, establish the plain FASTA bytes exactly correspond to that source, hash all 12 current index components against trusted component hashes, and independently compare the full FASTA transformed with the proven builder normalization against `.txt.bin`. The historical verifier and recovery records are a starting point, not an automatic current PASS. Preserve original Ns for prototype FASTA replay; the index's substituted text alone cannot recover their identity. Resolve provenance individually for the rest of the 464-panel set, for which the 12-assembly evidence does not apply.
3. Atomically create a **new staged prefix** and the exact current `crispr-poc-index-binding-v1` sidecar, with the plain FASTA SHA256 and verified current component SHA256s. Reuse existing components via verified read-only staging or component links where appropriate; this does not construct a new BWT or SA. The existing runner copies linked targets into a private snapshot and verifies hashes before and after search. A stage that can change without detection or bypasses that verification is not acceptable. Historical original/recovery directories should remain untouched.
4. Test the importer and checked runner on small pre-existing indexes with intentionally mismatched, truncated, changed and N-containing reference/index inputs; include alias names and multiple contigs. A region-restricted controlled search facility or another separately authorized bounded testing strategy is needed before an actual HPRC query. Re-run independent anchor, strand, boundary and CFD checks after the identity contract is satisfied. Do not silently replace a failed HPRC component or use older checksums as proof of current bytes.

The narrow implementation touchpoints would be a new existing-index import function alongside `bind_build`/`verify` in `compatibility/index_manifest.py`, a streaming, width-aware identity checker replacing the whole-file assumptions in `validation/index_check.py` where necessary, and a tested staged-prefix handoff to `integration/run_snapshot.py`. The existing `poc.cpp` `require_index_binding` and `fasta` paths must continue to fail closed; `src/pam_backend.rs::search` must not introduce an unchecked alternate path. A gzip-aware parser or a verified pre-existing plain FASTA is a separate input decision.

The original compressed source plus all index files remain important provenance even if a verified plain FASTA is available. Persistent identities must be guarded against time-of-check/time-of-use changes; the existing trusted-input/same-user limitation, unsigned manifest, helper paths and snapshot costs still apply. A read-only name or `chmod` alone does not make an index immutable. No genome partitioning is proposed or implemented here.

## Large-reference and memory considerations

The existing HPRC example has about 3.03 billion bases: larger than signed 32-bit maximum, but within the prototype's nominal **unsigned** 32-bit text limit including its sentinel. This is not an end-to-end validation. The source builder has an explicit >2^31 suffix-array construction path; there is no reason to invoke it for intact existing indexes. Prototype guide/edit counters use signed `int` but global indexed positions use Columba `length_t` and host `size_t`; practical correctness at billion-base offsets, 32-bit intermediate arithmetic and boundary handling remain to be tested independently. [The earlier index audit](../scalability_v2/INDEX_AUDIT.md) documents the precise assumptions and limitations.

Reading an existing `HG00735_pat` index entails approximately **11.46 GB of index components**, plus the original reference and working allocations. For the locally staged `HG00128_hap2` example, 11,439,717,349 bytes of index plus 3,075,802,891 bytes of plain FASTA is approximately **14.52 GB** of input files; the current runner makes a further full temporary snapshot per invocation. The verifier and search also materialize multiple reference/text copies and buffer output. Exact experimental-backend aggregate peak RAM and temporary-disk needs for a whole HPRC input have **not been measured** and cannot be inferred reliably from 3 Mb runs. Historical recovery builders used around 39.4–40.4 GiB RSS, but that is index-**construction** memory, not the memory needed merely to reuse/search an existing index. Set conservative time, aggregate-memory and temporary-disk limits before any future authorized pilot.

## Reproduction and limits

[inspect.py](inspect.py) takes the index root, panel TSV, resolved TSV and a new output directory. It reads metadata and component file stats, and for the 12 resolved sets reads dictionaries, the text length prefix and final byte. It does not hash or stream full component bodies; network filesystem reads can block or fail. [inspect_aliases.py](inspect_aliases.py) separately checks the four known name differences. [test_results.py](test_results.py) checks the retained evidence without accessing HPRC storage. The panel, resolved manifest and prototype executable hashes used in this inspection are captured in [provenance.json](results/provenance.json). No HPRC index, FASTA, production binary, frozen validation branch or defaults were modified, and no whole-genome query was started.
