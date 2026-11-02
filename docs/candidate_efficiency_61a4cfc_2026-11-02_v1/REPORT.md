# Candidate efficiency investigation — development commit 61a4cfc

**Recommendation:** develop and independently review an interval-level PAM
precheck before WFA2, preserving the existing candidate and overlapping-interval
expansion behavior. A disposable prototype produced identical full output
records on the bounded panels and reduced instrumented WFA2 calls by 96.3% on
the 8 kb oracle panel. It has **not** been merged into the production source,
and this evidence does not prove whole-genome correctness. Columba e=7
candidate discovery and its memory/SAM volume remain separate bottlenecks.

## Why candidates grow

For `-m 5 -b 1 -z 2`, `candidate_edit_distance_bound` asks Columba for all edit
distance ≤7 hits, plus a separate exact search. Columba does not know the
final NGG rule, maximum one gap event, match-fraction constraint, or
representative-PAM policy. Most raw e=7 candidates therefore cannot become
reported sites. The application reads the SAM into candidate vectors,
deduplicates by `(guide, contig, start, strand, CIGAR)`, and may verify a raw
span and many overlapping alternative spans with WFA2. It does not stream
candidates from Columba into WFA2. `NM` below is the Columba raw edit-distance
tag, not the final mismatch count.

A full external-sort analysis of the saved 1.94 GB CHM13–ELMOD2 e=7 SAM from
the stopped pilot, with no new whole-genome search, found:

| Measurement | Count |
|---|---:|
| Mapped e=7 rows | 17,653,389 |
| Duplicate exact candidate keys within e=7 | 0 |
| e=0 keys already present in e=7 | 1/1 |
| Raw intervals without an immediately adjacent strict NGG | 16,976,693 (96.17%) |
| Raw intervals with an immediately adjacent strict NGG | 676,696 (3.83%) |
| Candidates with no possible strict-NGG 18–22 bp interval overlapping the raw span | 6,502,192 (36.83%) |

Of the 17,653,389 e=7 rows, 14,948,701 (84.68%) have `NM=7` and 2,455,923
(13.91%) have `NM=6`; the remainder have `NM≤5`. The no-overlapping-PAM count
is a **conservative necessary-condition screen** for strict NGG, a 20-base
resolved guide and b1/z2: it considers every feasible 18–22 bp target span
that overlaps the candidate, both strands and complete PAMs at contig edges.
It can reject those candidates before WFA2 without relying on the raw CIGAR's
PAM alone. The actual program also permits noncanonical PAMs and `--no-filter`;
a production implementation must disable or correctly generalize the screen
for those modes, other PAMs, guide lengths and gap settings. The geometry
implementation was checked against an independently brute-forced test with
both strands, ambiguous bases and contig edges. Its full raw and unique counts,
`NM` distribution and external-sort method are in `genome_candidate_summary.json`.

**A raw interval with an invalid adjacent PAM is not necessarily an invalid
candidate.** Another overlapping interval can have a valid, distinct PAM. In
the 8 kb oracle panel, discarding raw candidates solely for their own PAM
lost 21/204 valid anchors. The 6,502,192 figure uses the more conservative
all-overlapping-spans condition; it is not the 16,976,693 direct-PAM count.
No exhaustive final whole-genome run was performed, so additional rejects,
final genome-wide sensitivity and the impact on reported anchors are unknown.

The earlier e=5 SAM contained 248,954 unique candidate keys. Of those, 190
were not present at the same start/strand in the e=7 SAM; different Columba
search bounds must not be assumed to produce nested raw records. All 3,027
prior *reported reference spans* had at least one same-guide, same-contig,
same-strand overlapping e=7 raw candidate (`prior_anchor_overlap.json`).
Overlap is only a necessary discovery relation: it is not proof of equivalent
alignment, same PAM, or final recovery. The new whole-genome verification
remains incomplete.

## Can the edit bound be tightened?

Not from 7 to 6 for the selected defaults. With a 20-base resolved guide,
`-f 0.75` permits 15 matches, five mismatches and a two-base reference bulge:
the valid edit distance is seven. In a direct indexed e=6 versus e=7 test
against the independent oracle, e=6 missed two sites, both with only a valid
five-mismatch plus two-base reference-bulge alignment:

- `GRCh38_nonrepeat_candidate_minus_mm5` on CHM13 chr9, minus strand, PAM
  start 31,315,521 (1-based), PAM `GGG`, guide `TTATATTAGACTAAGTTGCC`, target
  `TTTTATTATACTTTAAGTTTTA`.
- `GRCh38_repeat_candidate_plus_D1` on GRCh38 chr1, plus strand, PAM start
  100,001,610 (1-based), PAM `GGG`, guide `AGTTTTCTCTAGGGTTTTTA`, target
  `AGTTTTTTTAATGGAATCTTTA`.

For resolved 20-base guides and `-f 0.75`, a configured `-m 7` cannot admit
more than five substitutions, so seven is also a conditional upper bound
with a two-base reference gap; the current generic formula asks for nine.
This observation requires careful treatment of guide Ns, float thresholds,
`--no-filter`, gap configurations and Columba completeness before any change.
At `-m 7 -f 0.65`, seven mismatches plus a two-base reference bulge are
eligible, so nine remains necessary as a general bound. Splitting overlapping
e=5/e=6/e=7 searches does not reduce the e=7 set unless a candidate generator
can enforce genuinely disjoint edit classes and prove their completeness.

## Independent 8 kb real-reference oracle comparison

The prior independent [real-reference panel](../real_defaults_validation_49f20fc_2026-11-02_v1/REPORT.md)
contains four CHM13/GRCh38 windows, 96 guides and 204 strict PAM anchors,
including both strands and mismatch-plus-two-base-gap alignments. Experimental
SAM filtering used the original development binary with `--columba-sam`;
the filters are research scripts, not production behavior. Every emitted row
in every mode independently passed PAM, CIGAR, alignment limits and CFD audit.
A successful audit measures reported-row validity, not completeness.

| Mode | Candidate keys | Oracle recall | Precision | Verification seconds | Max RSS KiB |
|---|---:|---:|---:|---:|---:|
| e=7 + e=0 baseline | 6,284 | 204/204 | 204/204 | 13.08 | 5,812 |
| Raw candidate's direct PAM only — **unsafe** | 426 | 183/204 | 183/183 | 0.72 | 7,752 |
| Any overlapping-span PAM screen | 4,460 | 204/204 | 204/204 | 9.28 | 5,488 |
| e=6 + e=0 — **unsafe** | 1,545 | 202/204 | 202/202 | 3.82 | 4,640 |

The raw union contained 6,293 mapped rows; all nine e=0 rows duplicated e=7
keys. The overlap screen excluded 1,824/6,284 unique candidates and produced
the exact same full PAF row multiset as baseline. The direct-PAM screen's 21
misses and e=6 screen's two misses are retained in `oracle_benchmark.json`.
An initial e=6 run was discarded because its scratch output path collided
with its input SAM. The script now rejects that collision, and the table uses
the corrected rerun. Times are single runs, one thread, measuring verification
of imported SAMs; Columba discovery is excluded. Python prescreening took approximately 0.05 s,
but is a research implementation, not production overhead.

A separately built **instrumented** copy produced identical PAF records:
baseline used 1,105,628 WFA2 calls, overlap candidate screen 782,115, and
the interval-precheck prototype 41,077. These counters are from disposable
instrumented binaries, not the production executable. The prototype preserves
all candidates and the existing expansion logic; only immediately before
verifying a particular interval does it skip WFA2 when that interval's
complete PAM cannot pass the strict requested rule. In strict b1/z2 mode it
preserved all 204 oracle anchors, 204/204 precision and the exact PAF rows.
It disables the precheck for `--allow-noncanonical-pam` and `--no-filter`.
The synthetic oracle-backed m5, m6 and m7 fixtures, strict `-p GA`, b=0, z=1,
permissive-PAM and `--no-filter` controls all produced identical baseline and
prototype records (`interval_modes.json`). These checks are finite evidence,
not an independent proof of all flag combinations or sequences.

## Larger bounded regional benchmark

One ELMOD2 guide was run on each previously indexed 2 Mb CHM13 and GRCh38
chr1/chr9 region, including the higher-repeat chr9 proxies. Actual source
metadata identifies the chr9 windows as 30,000,001–32,000,000; an older
scaling report text says 34,000,001–36,000,000, so these measurements use the
actual `windows.json` and FASTA headers. These windows do not have exhaustive
oracle truth; equality to the original full PAF records is measured, and all
69 emitted rows independently passed eligibility and CFD audit.

| Window | Raw unique candidates | No overlapping PAM | Baseline verification s | Interval probe s | Baseline / probe RSS KiB | Reported rows |
|---|---:|---:|---:|---:|---:|---:|
| CHM13 chr1 | 13,033 | 5,064 | 34.97 | 1.67 | 10,156 / 10,328 | 10 |
| CHM13 chr9 | 13,272 | 5,604 | 37.65 | 1.67 | 10,200 / 10,328 | 24 |
| GRCh38 chr1 | 12,987 | 5,059 | 31.14 | 2.00 | 10,116 / 10,184 | 12 |
| GRCh38 chr9 | 13,269 | 5,593 | 34.45 | 1.57 | 10,376 / 10,348 | 23 |

Baseline and disposable prototype output records were byte-identical after
sorting in every window. RSS is GNU time's maximum process RSS, not aggregate
concurrency. These are paired, single-run, one-thread timings on a shared node;
they exclude Columba discovery, indexing and the screen preprocessing. The
2 Mb Columba e=7 phases separately took about 0.47–0.65 s and 65.6–68.3 MiB
maximum RSS. These small indexes do not represent whole-genome memory.
Exact candidate-key deduplication, and collapsing identical primary reference
intervals while ignoring CIGAR, each removed **zero** e=7 records in the
four 2 Mb samples; the 8 kb sample also had no repeated primary intervals.
This says nothing about overlapping *expanded* intervals or other guides.

## Next steps and limits

1. Review the interval-level precheck and its interaction with primary-hit
   control flow, expansion, scoring, strands, contig boundaries, Ns and all
   permissive flags. Then implement it in a separate change with exhaustive
   and randomized differential tests against the unmodified binary, including
   diverse mismatch-plus-two-base gaps and `-m 7`.
2. Keep exact-key deduplication correct, but do not expect it alone to resolve
   this sample's cost. Evaluate caching or deduplicating **expanded verification
   intervals** independently; the primary-interval counts above do not measure
   that opportunity.
3. Investigate streaming SAM ingestion and bounded memory with an exact
   deduplication strategy. This can reduce disk and application-side vector
   materialization, but the prior e=7 Columba phase itself consumed substantial
   memory and wrote its results only after mapping; streaming alone is not
   shown to eliminate Columba's internal cost. No streaming implementation was
   benchmarked here.
4. Research PAM-anchored candidate generation or sound, guide-aware search
   schemes that avoid enumerating most invalid-PAM hits. Any alternative needs
   an explicit completeness argument and independent oracle comparisons for
   both strands, 1–2-base bulges, high mismatch settings and repetitive data.

No production algorithm, default, PAM semantics, frozen branch or binary was
changed, and no new whole-genome search was launched. The prototype and
instrumentation are disposable copies; scripts and small summaries are
retained for review. The real-reference oracle is only 8 kb and the regional
samples use one guide; independent human review, whole-genome sensitivity and
resource requirements remain open.

All 160 regular Rust tests, the pinned standalone CFD comparison and 58 Python
tests passed. Analysis and benchmark scripts use separately versioned output
directories; no large SAM or index is committed.
