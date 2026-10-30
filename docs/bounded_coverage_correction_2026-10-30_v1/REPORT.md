# Bounded candidate and eligibility-aware alignment correction

**Result:** the development build recovers all 16 previously missing PAM
anchors in the tested panels, with no unexpected reported anchor. The frozen
preflight **does not pass**: this is bounded fixture evidence, independent
review is outstanding, and complete search sensitivity and production-scale
performance remain unproven. No complete search or production-binary change
was performed.

## Candidate geometry and the eleven expansion misses

The prior [broader report](../broader_coverage_1eb2f40_2026-10-30/REPORT.md)
records the independent raw Columba e=5/e=0 candidates, source sequence,
PAMs, and all oracle alignments. The eleven repetitive cases have a valid
WFA-accepted target with start offsets of -6, -5, -4, -3, +3 or +4 bases from
an overlapping raw candidate. The prior `2 * z` rule covered only -2..+2 at
`z=1`. `pipeline_trace.tsv` records each of the 16 sites, its raw SAM CIGAR
and span, required interval and signed start offset, oracle and selected WFA2
CIGARs, edit/mismatch/gap counts, costs, PAM coordinates and corrected output.
Offsets are relative to the candidate, on the forward reference coordinate
system on both strands.

A finite bound follows from interval overlap rather than an assumed repeat
length. Let guide length be L, total permitted gap bases B = `b*z`, candidate
reference interval `[s,e)` of span C, and maximum eligible target span
T = L+B. Any overlapping target starts between `s-(T-1)` and `e-1`, clipped
to the contig. Eligible reference spans are between `max(1,L-B)` and `L+B`;
each resulting interval must also overlap the raw candidate and fit the
contig. No smaller start bound follows from mismatch count alone: repeats can
create multiple valid matches with substantially shifted starts. Mismatch limit
m contributes to Columba's candidate edit bound `m+B` and to verification,
not to a safe tighter overlap bound. This geometry is conditional on Columba
returning an overlapping candidate; it cannot guarantee a candidate for every
biological locus or correct a non-overlapping/missing Columba result.

At L=20, b=1, z=1 and C=20, this produces at most 40 possible starts and
three spans, or 120 interval combinations before boundary, overlap, PAM and
original-interval exclusions, versus 15 start/span combinations in the former
five-start rule. The added work is bounded by candidate geometry and configured
limits, not a fixed arbitrary shift. For permissive or unusually large CLI
limits the cost can be much larger; those settings were not benchmarked.

## WFA2 alignment selection and five misses

The old anchored verifier selected one unconstrained affine-optimal WFA2
alignment and rejected the whole interval if that CIGAR violated the independent
mismatch, gap count, gap size or fraction limits. A different alignment of the
same interval could satisfy all of them. For example, `frozen2_I_10_plus`
`[45,65)` with GGG PAM `[65,68)` has a four-mismatch, no-gap oracle alignment
with affine cost 12; WFA2 chose zero mismatches plus I and D, also cost 12,
but two gap events exceed `-b 1`. For real `r3_I` on GRCh38 `[115,136)` with
GGG PAM `[136,139)`, the oracle has four mismatches plus one D, cost 18;
WFA2 chose I plus a two-base D, cost 13, violating both gap-event and gap-size
limits. The other three examples are tied at cost 18 but WFA2 selected three
gap events. These are actual competing alignments, not a relaxed eligibility
rule. The trace and prior `missed_sites.json` contain full aligned sequences,
CIGARs, mismatch and edit distances, gap sizes, strands and PAMs.

The correction keeps WFA2's result whenever it passes. When an anchored WFA2
result fails and the requested limits are exactly one bulge of size one, it
also tests the ungapped alignment or all positions of one I or one D as
appropriate to the two lengths. Each candidate is checked by the existing
CIGAR accounting and fraction/threshold validation. The best eligible
alignment under the same affine penalties (mismatch 3, gap open 5 plus one per
base) is selected with a deterministic CIGAR tie break. There are at most
L+1 gap positions, and no candidate with >1 gap event or gap length >1 is
created. Strict PAM filtering, the guide/reference orientation, coordinate
conventions and CFD scoring code are unchanged. The eligibility fallback is
not a general constrained aligner for other b/z configurations; direct
window scanning is also unchanged. Those modes retain their prior limitations.

## Regression and independent revalidation

`tests/bounded_coverage.rs` checks all 16 prior missing PAM keys using actual
retained indexed Columba SAMs, both strands and no duplicate output key. Unit
regressions exercise the overlap-bound edge cases, rejected WFA2 alignment,
all single-gap positions including terminal gaps, and a threshold-failure
control. An older test that required *no* recovery from an incorrect candidate
was revised to retain its meaningful guard against a **disjoint** shift and
against accepting the incorrect original interval: the wider, explicitly
bounded rule can now verify a valid one-base-overlapping target. Existing
PAM/CFD and coordinate tests remain in force.

The original and broader panels were reindexed and searched with the
separately built optimized development binary under unchanged `-m 4 -b 1 -z 1
-f 0.75 -p GG`; both actual e=5 and e=0 candidate searches ran. PAM-key recall
and precision, under the previously documented exact PAM-position criterion:

| Panel | Before recall | After recall | After precision |
|---|---:|---:|---:|
| Original synthetic | 33/33 | 33/33 | 33/33 |
| Original CHM13 windows | 3/3 | 3/3 | 3/3 |
| Original 12-assembly windows | 4/4 | 4/4 | 4/4 |
| Broader synthetic | 222/237 | 237/237 | 237/237 |
| Broader real windows | 19/20 | 20/20 | 20/20 |

No oracle PAM key remains missing and no unexpected key was emitted in these
panels. All 297 emitted PAF rows passed the independent sequence, CIGAR,
limits, three-base PAM and CFD auditor; each panel had exactly one output row
per guide/contig/strand/PAM anchor. Exact controls remained 20/20 on the
broader synthetic panel and 4/4 on the real panel. This does not replace the
prior bounded whole-genome raw exact comparison or establish final-output
whole-genome recall.

Strict interval/category counts remain separate; one PAM representative can
replace other valid intervals. On the broader synthetic panel the after counts
are exact 20/20, mismatch 53/67, I 37/64, D 43/68 and mismatch+bulge 84/377;
all emitted categories had 100% strict precision. Broader real-window counts
are exact 4/4, mismatch 5/5, I 0/4, D 6/10, mismatch+bulge 5/31; strict I
precision is undefined because no such category was emitted, and the other
emitted-category precisions were 100%. Categories can overlap at the same
interval. Full strict metric files and corrected PAF/SAM are retained in the
panel subdirectories; sequences and previous truth remain in the earlier report.
Multiple equally eligible alignments at one PAM remain possible and are not
asserted equivalent in sequence effect or CFD score; reporting still chooses
one deterministic representative.

## Cost, tests and provenance

In ten alternating-order optimized-binary runs, old/new medians were 53.09/
142.01 ms (2.67×) on the broader synthetic panel and 29.40/80.48 ms (2.74×)
on the real-window panel. Single-run peak RSS was 40,744/42,064 KiB and
27,272/25,896 KiB respectively. These tiny startup-bound panels do not
predict genome-scale runtime or memory. A separate test-only stats replay on
an isolated target reported 8,838 alternative intervals verified for 145
unique synthetic raw candidates, versus 1,288 previously; on real windows it
was 4,713 for 51 candidates versus 624 previously. The new run removed 17
identical verified hits and selected 237 anchors from 360 strict-eligible
synthetic alignments. The growth is material; measure and review before any
larger workload. `benchmark.json` and `stats_probe.log` retain the samples.

All 156 standard Rust tests passed, as did the normally ignored pinned CFD
test, and all 46 Python tests passed. The optimized development binary lives
under `/export/local/home/salehi/crisprapido-bounded-20261030-target/release/`
and has SHA256
`b4cd93f57f512a4727a3561b6c13289f03646b58964117f13058c1884c8c8b45`.
After a test-only source addition, this final binary was rebuilt and all five
panel PAFs matched the independently audited initial build (`11a43942…`)
byte-for-byte; timing and peak RSS above use the final binary.
The prior production-candidate binary remained SHA256
`f59b1bbf58e70b321359c7bd7b8687aad81bb874b22dc2f083319453726bc876`.
The source is based on `789b3c2`; the code change is contained in the commit
that adds this report. Generated indexes and verbose oracle reviews remain in
`/export/local/home/salehi/crisprapido-bounded-validation-20261030/`, outside
Git. Small result files and the trace are preserved here. No independent human
review, 100-guide × 12-assembly search, production binary replacement, merge,
or CRISPR-Bulge integration occurred.
