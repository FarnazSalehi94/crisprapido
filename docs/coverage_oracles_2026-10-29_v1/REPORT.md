# Coverage oracle development report

Status: **INCOMPLETE — frozen preflight gate is not satisfied.** No complete search was launched. The independent small-window prototype and tests are committed as work in progress; no claim of CRISPRapido recall or precision is supported yet.

## Algorithm

The prototype fetches a specified 1-based inclusive interval using `samtools faidx`, checks its exact length, and scans its sequence and reverse complement. It enumerates all target starts and lengths guide length minus one, equal, and plus one; ungapped comparisons and every possible single-base gap position are evaluated directly with character comparisons. It requires a complete oriented A/C/G/T NGG PAM, up to four mismatches, at least 75% guide matches, and no more than one single-base gap. Exact mode restricts to ungapped zero mismatches. Coordinates are genomic zero-based half-open, and overlaps are retained. No WFA2, CIGAR reconstruction, candidate discovery, CFD score, or output score threshold is used.

Each emitted JSON object is an alignment representation including aligned strings and gap position; multiple representations may identify the same `(contig, start, end, strand)` locus. A locus comparison can collapse exact duplicates using this tuple, but shifted overlapping intervals must not be presumed equivalent: retain and inspect the specific representative and selection boundary. Exact versus imperfect classification must be compared separately. The oracle excludes ambiguous target and PAM bases; behavior of indexed ambiguous sequences and potentially terminal gaps requires further investigation against the release implementation.

## Checks performed

Hand-specified synthetic examples test exact plus and minus coordinates, oriented PAM, incomplete and ambiguous PAM, overlapping sites, mismatch with each gap direction on both strands, and ambiguous targets. These tests exercise the oracle independently; they are not a candidate-discovery comparison.

## Missing work and limitations

No comparison against CRISPRapido raw exact SAM or final PAF has been performed on synthetic or real indexed windows. Expected/recovered/missed/unexpected counts, recall, and precision are therefore **not measured**. A genome-wide streaming multi-pattern scanner, explicit raw-candidate comparison and overlap-representative matching, fixed windows across all twelve assemblies, further independent truth cases including four mismatches, and resource measurements are still required. The current script materializes a window and reverse complement, does not stream whole genomes or accept multiple guides, and does not deduplicate equivalent gap placements. Its inability to retrieve flanking PAM beyond a specified window boundary means callers must include context and distinguish window edges from true contig edges. CFD and downstream output filtering remain separate checks covered by the emitted-row auditor, not by this prototype.

The frozen preflight gate remains blocked until these missing comparisons and independent review are complete. Do not run the complete search on the basis of this report.

## Follow-up: fraction filter and controlled imported-candidate probe

Status remains **INCOMPLETE; frozen preflight gate does not pass**. The
multi-guide exact scanner now streams FASTA, includes both strands, overlaps,
NGG and line/record boundaries; tests compare it to the independent window
oracle. It has not yet been benchmarked on whole genomes or compared to raw
Columba SAM. The bounded oracle and final PAF must not be equated: reporting
selects representatives among overlapping same-strand hits.

`-f` does not enter the Columba edit bound or candidate generation. In both
anchored and window verification it checks `100 * matches / non_N_guide_bases`
against `100 * f` (f32), and independently requires at least one match, unless
`--no-filter` bypasses filtering. M and = each increment matches; X does not,
including X at guide N (which also does not count as a mismatch). I consumes a
guide base without a match; D consumes a reference base without a match and
does not consume a guide base. Guide N bases are excluded from the denominator,
but M/= at N still increment the numerator. Thus this is not in general a
strictly bounded fraction for guides containing N. Free reference-end D
operations are trimmed in window alignment; anchored alignment retains them.
The alignment selected by WFA is filtered, not every possible alignment.

For 20-base A/C/G/T guides and `-m 4 -b 1 -z 1`, a full alignment that satisfies
those limits has at least 15 matches: four X and one I give 15/20; D instead
leaves at least 16/20. Consequently 0.75 adds no exclusion to such alignments.
This does not establish complete candidate discovery or claim that the verifier
finds every qualifying alignment. For shorter guides, it can exclude candidates.

`fraction_probe.py` tests each oracle locus independently using an imported SAM
record and the release binary SHA256 `f59b1bbf58e70b321359c7bd7b8687aad81bb874b22dc2f083319453726bc876`.
Full output is in `fraction_probe.jsonl`. Coordinates are zero-based half-open
relative to each window. The synthetic 20-base guide
`ACGTCAGTACGATCGTACGA` produced the same two accepted loci in both modes:
`[0,20,+]` (16 matches, four mismatches) and `[1,20,+]` (16 matches,
three mismatches, one I). CHM13 `CHM13#0#chr1:100001-100500`, guide
`GAACATTTACCATGATAAAG`, produced the same three accepted loci in both modes:
`[0,21,+]` (20M plus D), `[1,21,+]` (20M), `[2,21,+]` (19M plus I).
The 10-base guide `ACGTCAGTAC` and sequence `TGGTAGTACAGG` produced one
mode-dependent locus `[0,9,+]`: 0.0 emits `2X2=1I5=` (7/10 matches,
two mismatches, one I); 0.75 rejects it. The probe uses deliberately supplied
candidate intervals and cannot measure discovery recall or final multi-candidate
overlap behavior. No broad claim is made about other windows or guides.

Still required: indexed Columba comparisons including raw exact SAM, independent
truth and final overlap-boundary analysis across the frozen assembly windows,
category-specific recall/precision and discrepancy investigation, scanner
performance/resource measurements, and independent review. Therefore category
recall and precision are **not measured**, and the frozen preflight remains
blocked. No complete search was launched and the default was not changed.

## Indexed-window follow-up, 2026-10-30

This section supersedes the earlier statements that no indexed comparison or
coverage counts had been performed. **The frozen preflight gate still does
not pass.** The interval-level discrepancies below are real differences in
reported representations, and their downstream significance remains unresolved.

Three fresh small reference sets were indexed with the pinned Columba builder,
then searched with Columba `-a all -m edit -e 5 -t 1`, a separate `-e 0` search,
and the unchanged release CRISPRapido binary with `-p GG -m 4 -b 1 -z 1
-f 0.75 -t 1`. The automatic CRISPRapido invocation also ran both candidate
searches. These were window indexes, not searches of the original whole-assembly
indexes. Input FASTA, guides, raw SAMs, PAF, oracle truth, metrics, discrepancy
records and auditor summaries are retained in `indexed_synthetic/`,
`indexed_real_chm13/` and `indexed_12_windows/`. Builder SHA256
`5987ad90e0c4a061cbf4aac7d6b087a284ccf7fe0f9e1835219c567c5e8d1ad2`;
Columba SHA256 `7e44974c81dad4943292fa64a71e54073908b4909ce327480464d769322c40d1`;
release executable SHA256 `f59b1bbf58e70b321359c7bd7b8687aad81bb874b22dc2f083319453726bc876`.

The synthetic panel contains two 20-base guides and 28 separately flanked
sequences: exact, one mismatch, four mismatches, I, D, four mismatches plus I,
and four mismatches plus D, each on both strands for each guide. All 28 designed
intervals were present as intended in the final PAF; each guide had seven hits
per strand. The real CHM13 panel contains 1 kb each from chr1:100001-101000
and chr2:200001-201000; its two guides were selected from those windows to
ensure exact examples, so they are not an unbiased guide sample. A second real
panel contains 1 kb at positions 100001-101000 on each of the first two FAI
contigs from each of the twelve resolved assemblies, indexed together, and
uses the same two CHM13-derived guides. `windows.json` gives the exact source,
contig and coordinate for every record. Guide `real2` had no eligible site in
that 24-window panel. These windows do not certify the original assembly indexes
or genome-wide coverage and do not sample natural negative-strand sites here.

The following metrics use **strict interval and category identity**. The oracle
enumerates all single-gap positions and eligible intervals. A locus may belong
to several categories because distinct valid alignments exist for the same
interval; thus category denominators overlap. An emitted CIGAR has exactly one
category. `I` means guide-consuming insertion (RNA bulge, 19-base reference
span); `D` means reference-consuming deletion (DNA bulge, 21-base span).
Mismatch + bulge combines both directions. Recall is category-matching emitted
intervals divided by all oracle intervals in that category; precision is those
matches divided by emitted intervals in the category. Undefined denominators
are shown as —, not 100%. These counts do not equate overlapping boundaries.

| Panel | Category | Oracle | Emitted | True positive | Recall | Precision |
|---|---|---:|---:|---:|---:|---:|
| Synthetic | Exact | 4 | 4 | 4 | 100% | 100% |
| Synthetic | Mismatch only | 8 | 8 | 8 | 100% | 100% |
| Synthetic | I | 8 | 4 | 4 | 50% | 100% |
| Synthetic | D | 9 | 4 | 4 | 44.4% | 100% |
| Synthetic | Mismatch + bulge | 45 | 8 | 8 | 17.8% | 100% |
| CHM13 | Exact | 2 | 2 | 2 | 100% | 100% |
| CHM13 | Mismatch only | 0 | 0 | 0 | — | — |
| CHM13 | I | 2 | 0 | 0 | 0% | — |
| CHM13 | D | 3 | 0 | 0 | 0% | — |
| CHM13 | Mismatch + bulge | 5 | 0 | 0 | 0% | — |
| 12 assemblies | Exact | 1 | 1 | 1 | 100% | 100% |
| 12 assemblies | Mismatch only | 1 | 1 | 1 | 100% | 100% |
| 12 assemblies | I | 1 | 0 | 0 | 0% | — |
| 12 assemblies | D | 2 | 0 | 0 | 0% | — |
| 12 assemblies | Mismatch + bulge | 6 | 0 | 0 | 0% | — |

Raw strict-NGG `e=0` exact SAM recovered 4/4, 2/2 and 1/1 oracle exact loci,
respectively, with no unexpected strict exact locus. The multi-guide streaming
exact scanner independently agreed with those window-oracle loci. Raw `e=5`
SAM contained the identical interval for 28/57, 2/7 and 2/8 distinct oracle
loci, respectively. Columba's raw edit SAM is a candidate superset without PAM
filtering, so its unfiltered rows are not precision predictions for strict NGG;
no raw-edit precision is claimed.

All category-level final misses are enumerated in each `discrepancies.json`:
46 in the synthetic panel, ten in CHM13, and nine in the twelve-assembly panel.
A single locus can contribute more than once when alternative alignments place
it in multiple categories. Eight synthetic category misses share an identical
raw candidate interval but the final CIGAR represents another category; the
remaining 38 synthetic category misses and all 19 real category misses lack an
identical raw `e=5` interval. Every miss has an overlapping raw edit candidate
and final PAF interval with the same guide, contig, and strand; the files record
each actual pair and category. For example, on an exact site, an oracle-valid
19-base interval with a terminal guide gap overlaps the 20-base exact PAF hit.
Columba does not return that exact 19-base interval in these cases. Other misses
are alternative suboptimal gap placements at the same interval. Therefore the
evidence identifies a stage and a competing representation but does not prove
that every shifted interval is biologically equivalent, or isolate whether a
hypothetical separately supplied interval would survive later overlap selection.
No unexpected final category/interval occurred in these panels. A separate
emitted-row auditor verified PAM, coordinates, CIGAR, limits and CFD for all 28,
2 and 2 PAF rows respectively, with zero errors.

The comparison and fixture generators are in `scripts/release_validation/`;
`test_indexed_coverage.py` recomputes metrics and checks streaming exact-scan
agreement and discrepancy overlap evidence from retained fixtures. Synthetic
reference construction uses seed 20261030. Coordinate values in retained files
are zero-based half-open within their named window; add the recorded source
window offset for genomic positions. SAM/PAF headers retain the original local
run paths. The independent oracle enumerates A/C/G/T only, requires complete
NGG, and examines all single-base gap placements; it is not WFA2 or Columba.

Remaining gate limitations: the 100-guide genome-wide exact oracle has not been
run against the 12 whole assemblies or their original raw exact SAM; its runtime
and resource envelope are not established. Window inputs cover only tiny slices,
two chosen guides, and few real positive examples, and are not a representative
sensitivity estimate. Independent review of the oracle/comparison and overlap
policy remains outstanding. Resolve the interval-level discrepancies and obtain
that review without changing the frozen criteria before declaring the preflight
passed. No 100-guide × 12-assembly search, binary modification, or merge was
performed.

## PAM-anchored discrepancy review and whole-genome exact pilot

**Decision: the frozen preflight gate does not pass.** The unchanged strict
interval metrics in the preceding section remain valid. A conservative
PAM-anchored analysis exposes missing distinct target sites, and the bounded
whole-genome exact comparison is not the frozen 100-guide, 12-assembly scope.

### Independent comparison criteria

The comparison has two separate levels. An **identical alignment** requires the
same guide, contig, strand, reference interval and reconstructed aligned guide
and target (equivalently, the same position-specific edit structure after
replay). Gap placements are retained individually. Repeated sequence may admit
other normalizations; no gap-shift normalization is assumed here. The existing
strict category/interval metrics are preserved separately and do not require
identical gap placement within their category.

For a **PAM-anchored target site**, the key is `(guide ID, contig, strand,
PAM start, PAM end, oriented PAM sequence)`, with zero-based half-open genomic
coordinates. A positive key must have a full A/C/G/T NGG and at least one
oracle-valid full-guide alignment at the frozen bounds. A reported hit counts
as the same target site only when its independently reconstructed target,
PAM, guide identity and alignment are valid and it has that identical key.
The reference sequence and reconstructed edit structure for every comparison
are retained in `alignment_review.json`; matching a site does not assert that
its alternative edit structures have identical effects or CFD scores. Identical
intervals with distinct valid alignments, and distinct intervals that end at
the same plus-strand PAM or start at the same minus-strand PAM, are classified
**B** (alternative alignments at the same PAM-anchored site), rather than
identical alignments. Overlap with a different PAM position is **C** (a
separate PAM-anchored target site); it is not used to recover the oracle site.
This definition distinguishes genomic target sites, not experimentally proven
cleavage outcomes. No transitive overlap clustering is used.

`resolve_bulge_discrepancies.py` reconstructs every oracle alignment, target
sequence, PAM and coordinates, normalized CIGAR, mismatch count and gap
position; it also replays every overlapping raw SAM CIGAR (reversing the SAM
operation order on the minus strand) and final guide-oriented PAF CIGAR. The
three `alignment_review.json` files record every interval/category discrepancy
and every competing alignment. The earlier `discrepancies.json` files and
strict calculations remain unchanged. Class counts below are **discrepancy
records**, not deduplicated sites: A identical 0; B same PAM with a different
alignment; C different PAM-anchored site; D no reported competitor 0; E no
classification 0. Class C establishes a missing site under the stated key;
it does not establish a distinct in-vivo biological outcome.

| Panel | B | C | Oracle PAM-anchored sites | Recovered | Locus recall |
|---|---:|---:|---:|---:|---:|
| Synthetic | 40 | 6 | 33 | 28 | 84.85% |
| CHM13 windows | 8 | 2 | 3 | 2 | 66.67% |
| Twelve-assembly windows | 6 | 3 | 4 | 2 | 50% |

These are small, deliberately enriched fixtures; the real panels reuse one CHM13
window, so their counts must not be pooled as independent evidence. Six, two,
and three C records correspond to five, one, and two missing distinct
PAM-anchored keys respectively. The repeated CHM13 key appears in both real
panels. There are seven distinct missing keys across these fixtures: five
synthetic, one CHM13 and one HG04184. Each omitted site has an oracle-valid
21-base reference span and adjacent GGG PAM. The other final hit at each site
uses a shifted AGG PAM. Detailed missing keys and competing sequences are
machine-readable in each panel's `alignment_review.json`.

For example, `g1` on `g1_exact_plus` has an oracle-valid `[60,81)` target and
PAM `[81,84)` GGG, with `20=1D` or `19=1D1=` and zero mismatches. The final
hit is `[60,80)`, `20=`, with PAM `[80,83)` AGG. Those PAMs overlap but have
different genomic positions. Calling them one recovered site would hide the
GGG-anchored candidate.

### Stage investigation and unresolved causes

None of the eight panel-specific missing intervals appears at the identical
coordinates in the retained `e=5` raw Columba SAM. Every one has a different
overlapping raw candidate. This establishes non-emission of the exact interval
by Columba, but the production verifier can expand some candidates, so absence
from raw SAM alone does not locate the ultimate loss.

`probe_missing_loci.py` separately supplies each missing interval as a synthetic
SAM candidate to the unchanged release binary. Its strict-mode output and
permissive diagnostic controls are recorded in each panel's
`imported_locus_probes.json`; the controls do not change the frozen results.
Both synthetic minus-strand missing GGG intervals are emitted at their exact
coordinates when supplied independently, showing that anchored WFA2
verification, strict NGG filtering and single-candidate reporting can pass.
The six panel-specific plus-strand probes instead emit an overlapping hit at a
different PAM. A test-only Rust probe, appended to a disposable copy of the
unchanged source, calls the actual anchored WFA2 verifier on all five distinct
plus-strand guide/target sequences: **all five pass** at the frozen limits
with a terminal D, and each associated PAM is strict NGG. In the individual
imported-candidate experiments these verified, PAM-eligible plus intervals
were displaced by overlap representative selection. The test source and log
are `scripts/release_validation/missing_pam_wfa_test.rs` and
`genome_exact_chm13_2guides/wfa_stage_probe.log`. It is a debug test in a
disposable source copy, not a modification of the release binary.

In the original automatic search, raw Columba omitted all eight panel-specific
intervals. The CLI does not expose per-interval expansion or verification
counters, so the retained logs do not establish whether each absent interval
was ever expanded and then removed by selection during that particular run.
The tests establish an independent discovery gap at exact-interval granularity,
and demonstrate selection loss when the five plus intervals are supplied;
they do not claim a complete internal trace of the automatic search.

### Genome-wide exact-match pilot

After a 10 Mb pilot (0.15 s wall, 14,500 KiB peak RSS), the optimized streaming
scanner scanned the **entire verified CHM13 FASTA** for the first two frozen
20-base guide records. This is a bounded two-guide, one-assembly exact-only
validation, not the prohibited 100-guide × 12-assembly search. In parallel,
Columba searched the recovered verified CHM13 index with `-a all -m edit -e 0
-t 1`. Each SAM result was checked by a fresh `samtools faidx` fetch of the
entire target and oriented PAM. The independent scanner found two strict NGG
loci, one on each strand and one per guide; Columba recovered **2/2 (100%
recall)**, with 2/2 strict SAM loci correct (100% precision), no misses and no
unexpected strict loci. Details, exact SAM, guide FASTA, oracle and timing files
are in `genome_exact_chm13_2guides/`.

Scanner elapsed time was 61.11 s and peak RSS 14,772 KiB; Columba elapsed 4.19 s
and peak RSS 11,558,660 KiB. Both exited zero. The pilot source is the verified
CHM13 FASTA identified in the recovery manifest, SHA256
`f267b6e3083a0da98d50a3e84d36a8d2de789665eba9766bac57b529e306c33d`;
scanner SHA256 `ac982938a639a16a1c8cd204f5e3b68ad75da92a84fd13384c2eb57972c2a515`.
Large indexes, compressed genome and 10 Mb pilot stay outside Git under
`/export/local/home/salehi/crisprapido-coverage-20261030-panel1/`.
This exact-only result does not assess final PAF overlap selection, imperfect
search sensitivity, the other frozen guides, other assemblies, or index content
at bases excluded by the oracle's A/C/G/T rule. A larger scope, per-stage
instrumentation of the original automatic searches, and independent review
remain required before the frozen gate can pass.

## Subsequent development correction

The later [PAM-anchor correction report](../pam_anchor_correction_2026-10-30_v1/REPORT.md)
traces the seven sites and records a separately built correction and fresh
indexed revalidation. The earlier measurements in this report remain the
historical baseline. The frozen preflight remains not passed.
