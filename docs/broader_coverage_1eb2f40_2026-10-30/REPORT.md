# Broader coverage validation of 1eb2f40

**Decision: frozen preflight does not pass.** Two exact-only whole-genome
comparisons passed for ten guides each, but a broader bounded oracle found 16
missing strict-NGG PAM-anchored sites. Eleven occur in repetitive synthetic
contexts and five have independently valid alignments that the selected WFA2
alignment fails to verify. No unexpected reported anchor was found in these
panels. Independent human review remains **outstanding**. The 100-guide ×
12-assembly search was not run.

## Frozen implementation and dependencies

Source commit `1eb2f4073adaea85f76e0e677edab7b0e076fb92`, tree
`7a9f0099057e1d78f4888f0a0d5c20e08c5a12aa`, branch
`fix/pam-cfd-validation`. The separately copied optimized executable used for
these tests has SHA256
`8a207b61a3cc5ad2b5f208115a19dbf3eaaf76b7064519aab0340a154aa9cbd8`.
Its isolated release build used Rust 1.85.1, `--release --offline --locked`,
Cargo.lock SHA256 `93564e52dcafb54bfdacd2b7c6d69ac806aa74fb00ea1a20a6aee75f13c7f0f2`,
lib_wfa2 revision `c608c436a6753d2c21c97d9f5c338efae99d042b` and WFA2 static
library SHA256 `c19e330629971317d92aad90fa7c6c532441e4e5b3fccc72c18fd4016c927dbf`.
Columba executable SHA256 is
`7e44974c81dad4943292fa64a71e54073908b4909ce327480464d769322c40d1`;
the index builder is
`5987ad90e0c4a061cbf4aac7d6b087a284ccf7fe0f9e1835219c567c5e8d1ad2`.
The Guix 2.35 loader and GCC 12.3 library paths are the ones in the earlier
release manifest; the Columba wrapper SHA256 is
`2367ef7b053a44c83f76b059ed66a8f671f83773d3fe69240cda50893408f28a`.
The FASTA auditor used samtools 1.23.1 (binary SHA256
`68315154fcdd270e119d91ccb1329e265095dcab74a56db68ed9305d82ea9827`).
Score matrix, scanner, oracle and guide hashes are in `hashes.txt`. The two
whole-assembly FASTA/index pairs are the CHM13 and GRCh38 entries in
`docs/index_recovery_2026-10-29_v1/verified_assemblies.tsv`, including their
prior full index/source verification; these comparisons did not rehash the
entire inputs. Tests used `-p GG -m 4 -b 1 -z 1 -f 0.75 -t 1`, preserving
all biological bounds and using one thread for these bounded checks.

## Independent exact comparison

The independent streaming scanner examined both strands of the entire frozen
CHM13 and GRCh38 FASTAs for the first ten guides of the frozen guide list.
Actual Columba `-a all -m edit -e 0` searches used the corresponding recovered
whole-assembly indexes. Every mapped SAM row had its full sequence and
oriented three-base PAM checked against the exact FASTA using samtools.

| Assembly | Oracle exact NGG loci | Recovered in raw exact SAM | Precision | Strands |
|---|---:|---:|---:|---|
| CHM13 | 10 | 10 (100%) | 10/10 (100%) | 6 plus, 4 minus |
| GRCh38 | 10 | 10 (100%) | 10/10 (100%) | 6 plus, 4 minus |

No unexpected strict exact locus or missed exact locus appeared. There were
no overlapping exact sites in these two real guide subsets. Synthetic exact
controls include two distinct overlapping guide/PAM sites on both strands;
all 20 synthetic exact oracle loci were in raw exact SAM and all 20 were
reported at identical intervals. Whole-genome exact results are **raw Columba**
comparisons; they do not prove final WFA/representative output coverage.
Scanner elapsed times under concurrent execution were 157.41 s for CHM13 and
155.03 s for GRCh38, at <=14,760 KiB peak RSS. Columba exact elapsed 4.14 s
and 33.44 s respectively, with peak RSS about 11.0 GiB; cache and concurrent
I/O make these single-run times unsuitable for general extrapolation.

## Broader bounded panels

The deterministic synthetic panel uses five 20-base guides (two frozen,
one repetitive and two overlapping-PAM controls) and 134 separately flanked
contigs. It includes mismatches at guide positions 0, 4, 9, 14 and 19;
two four-mismatch patterns; guide and reference bulges at 0, 5, 10, 15 and
19; mismatch-plus-bulge cases; repetitive sequence; adjacent `AGG`/`GGG`
PAMs; and both strands. The real-window panel uses 16 deliberately designed
20-base queries against four **unchanged** 1 kb windows from CHM13, GRCh38,
HG04184 and HG00128: exact and four-mismatch designs plus gap-oriented query
designs. Guides are chosen to create local tests, so these panels are not
unbiased samples of the genome. The natural real-panel observed hits here are
plus-strand only; both strands are covered by the synthetic panel. Reference
and guide FASTAs, design metadata, raw e=5/e=0 SAMs, PAFs and independent
metrics are retained in `synthetic/` and `real/`.

A locus is keyed by guide, contig, strand and exact genomic position of the
complete oriented PAM, as specified in the prior review. Different PAM
positions remain distinct even when protospacers overlap. The exhaustive
oracle evaluates A/C/G/T targets, every 19-, 20- and 21-base interval and all
single-gap placements under the frozen limits. Recall counts oracle PAM keys
with any valid reported alignment at that key; precision counts reported keys
also present in the oracle. Matching a key does not assert identical edit
structure or biological effect. Both categories and intervals can share a key.

| Panel | Oracle PAM loci | Recovered | Recall | Reported correct / reported | Precision |
|---|---:|---:|---:|---:|---:|
| Synthetic | 237 | 222 | 93.67% | 222/222 | 100% |
| Real windows | 20 | 19 | 95.00% | 19/19 | 100% |

By oracle category, synthetic PAM-key recovery is exact 20/20, mismatch-only
61/67, guide insertion 64/64, reference deletion 68/68, and mismatch+bulge
222/237. Real-window recovery is exact 4/4, mismatch-only 5/5, guide insertion
4/4, reference deletion 10/10 and mismatch+bulge 19/20. These category sets
**overlap**: for example, an exact locus may also have a valid mismatch-plus-
bulge representation. This is not evidence that each edit representation is
reported. `locus_categories.json` preserves the keys.

Strict interval-and-category identity results, independently of PAM-key
recovery, are retained separately:

| Panel | Exact | Mismatch only | I | D | Mismatch + bulge |
|---|---:|---:|---:|---:|---:|
| Synthetic | 20/20 | 45/67 | 36/64 | 43/68 | 78/377 |
| Real windows | 4/4 | 5/5 | 0/4 | 6/10 | 4/31 |

Strict-category precision was 100% for every emitted category; real-window I
precision is undefined because no I category was emitted. Different valid
alignments at the same PAM are not presumed to be identical intervals or edit
representations. The independent all-row auditor passed 222 synthetic and 19
real PAF rows, including CIGAR, reference/PAM, thresholds and CFD; no invalid
or duplicate output PAM anchor was found in either panel.

## Every missing site

`missed_sites.json` has all 16 missing keys, sequences, oracle alignments,
overlapping raw candidates from **both** e=5 and e=0 SAMs, isolated imported
candidate results, and direct WFA2 diagnostic results. It distinguishes the
following confirmed observations, rather than merging overlapping sites:

- **Eleven synthetic repetitive keys:** at least one eligible interval is
  accepted by anchored WFA2 when supplied directly, but every such accepted
  interval starts at least 3–6 bases away from every overlapping raw candidate.
  The current `z=1` fallback considers starts at most two bases away. This
  bounded expansion cannot reach those tested accepted intervals from those
  candidates. Other potentially valid alignments or non-overlapping candidates
  are not assumed equivalent. All eleven are retained as genuine missing PAM
  keys under the stated oracle, not recategorized as duplicate output.
- **Four synthetic keys and one real-window key:** every oracle-valid interval
  at the missing key was separately supplied to the same anchored WFA2 logic;
  the selected WFA2 alignments failed frozen gap limits. For example, synthetic
  `frozen2_I_10_plus` `[45,65)` has a valid four-mismatch ungapped oracle
  alignment and GGG PAM at `[65,68)`, but WFA2 chose an I and a D (two gap
  events), so verification rejected it. Real `r3_I` on GRCh38 `[115,136)`
  has a valid four-mismatch one-D alignment and GGG PAM at `[136,139)`;
  WFA2 chose an I plus two-base D and rejected it. Raw CIGARs and mismatch/gap
  accounting are in `wfa_probe_raw.log`. These direct test-only diagnostics
  used a disposable source copy of the frozen implementation; they did not
  alter the frozen release binary. Exact failure stage in the original
  automatic path can include candidate non-discovery or bounded expansion as
  well; the direct tests establish that supplying the intervals is insufficient.

No unexpected strict PAM key was emitted. Lower strict interval/category
recalls also include alternative intervals and alignments sharing a recovered
PAM; each listed discrepancy is retained in the generated external review, and
the stored metrics keep them distinct. Full verbose alignment reviews and
indexes remain outside Git under
`/export/local/home/salehi/crisprapido-broader-coverage-1eb2f40-20261030/`.

## Expansion, duplicates and resource effects

A diagnostic unit test in a disposable copy of the frozen source replayed the
real e=5/e=0 union and exposed the verifier's existing internal counters; it
is not the timed release executable. In synthetic/real panels respectively:

| Measure | Synthetic | Real |
|---|---:|---:|
| Raw e=5 + e=0 mapped rows | 165 | 55 |
| Unique raw candidates after exact candidate deduplication | 145 | 51 |
| Alternative intervals verified | 1,288 | 624 |
| Alternate hits accepted by WFA2 | 351 | 10 |
| Exact duplicate verified hits removed | 0 | 0 |
| Strict-PAM eligible verified hits | 308 | 19 |
| Distinct output PAM anchors | 222 | 19 |

Thus 86 synthetic strict-eligible alignments shared an already represented
anchor and were reduced to one final representative; no duplicate final anchor
remained. Candidate deduplication removed 20 and four duplicate e=5/e=0
records respectively. The synthetic repetitive panel has a substantial
candidate-expansion multiplier; no threshold has been established to label it
acceptable or excessive for production workloads.

Ten alternating runs of original versus frozen corrected optimized binaries
on identical small indexes gave median 50.94 vs 52.75 ms (synthetic, 1.036×)
and 29.32 vs 29.95 ms (real, 1.022×). Separate one-run peak RSS observations
were 40,292 vs 40,624 KiB (synthetic) and 27,428 vs 27,424 KiB (real). Startup,
cache, and fixture size dominate; no genome-scale runtime or memory conclusion
is supported. `benchmark.json` and the stats log preserve the measurements.

## Tests and review status

All 152 standard Rust tests passed, as did the normally ignored standalone
pinned CFD test; all 46 Python tests passed. An initial Rust test invocation
accidentally reused a test executable from the disposable instrumentation
checkout because it shared a debug target directory and failed on a missing
diagnostic environment variable. The worktree source was then forced to
rebuild and the full suite passed; the frozen release binary hash was checked
unchanged. Future diagnostic builds should use a separate target directory.
The new fixture checks validate overlap, miss evidence and output-anchor
uniqueness.

No independent reviewer was available in this session. This report and the
oracle assumptions, PAM-anchor equivalence, gap/terminal-gap policy, and 16
miss classifications therefore require independent review; self-review is not
counted as satisfying the gate. No production binary, old results, index, or
frozen parameters were overwritten. No merge or CRISPR-Bulge integration was
performed. The preflight remains blocked by the observed misses, missing
independent review, and the uncompleted frozen 100-guide / 12-assembly coverage.

## Subsequent development correction

The later [bounded correction report](../bounded_coverage_correction_2026-10-30_v1/REPORT.md)
traces all 16 sites and records a separately built correction and revalidation.
The measurements above remain the historical baseline. The frozen preflight
still does not pass.
