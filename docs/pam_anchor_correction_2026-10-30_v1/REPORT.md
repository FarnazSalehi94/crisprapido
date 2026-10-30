# PAM-anchor correctness correction

**Result:** the separately built development executable recovers every
PAM-anchored oracle site in the three small indexed panels. This is not a
claim of whole-genome completeness. The frozen preflight remains **not
passed**: the frozen multi-guide, multi-assembly exact coverage and independent
review remain incomplete, and broader imperfect-site sensitivity is unproven.
The 100-guide × 12-assembly search was not launched.

## Root causes and pipeline trace

The independent oracle enumerates full-guide alignments with at most four
mismatches, one single-base gap, >=75% matches and a complete strict NGG PAM.
It identifies five missing PAM keys in the synthetic panel, one in the CHM13
panel and two in the twelve-assembly panel. One CHM13 key is in both panels;
there are seven distinct sites, represented by eight panel-specific records.
The records include both strands. `pipeline_trace.tsv` identifies each guide,
interval and PAM, the raw candidate and the corrected final CIGAR. The prior
`alignment_review.json` files retain the full independent sequences and all
alignment representations.

For **each** of these records, actual Columba `e=5` returned one overlapping
20M candidate at a different PAM, but not the oracle's 21-base interval.
The original WFA2 verifier accepted the 20M candidate (perfect or one
mismatch), and its adjacent PAM was canonical AGG or TGG. The original
`should_expand` condition skipped alternative intervals for such an accepted
canonical primary without D. Consequently the missing interval was never
sent to WFA2, never had its PAM checked, and never entered representative
selection in that original fixture search. This conclusion follows from each
retained raw SAM, independently reconstructed primary alignment, the frozen
limits, and the original control flow; no stage counter was present in the
old binary. The retained raw exact-SAM union does not introduce another
candidate at the missing interval.

Separately, supplied-candidate experiments showed a second defect: the old
reporter grouped overlapping same-strand hits without comparing PAM anchors,
and could discard an independently WFA-verified, strict-NGG plus-strand site
in favor of a hit at another PAM. The prior direct WFA2 probe and imported
SAM records confirm this selection failure when those intervals are supplied;
they do not mean the originally unexpanded intervals entered reporting.
All missing sites have oracle-valid terminal D alignments, and their recovered
final CIGARs are in the trace. The corrected indexed runs now discover a
nearby PAM through expansion, verify the candidate, require the strict PAM,
and output it as a distinct representative.

## Change

For a passing canonical primary that previously skipped expansion, inspect
nearby alternative intervals and WFA-verify those having a different eligible
PAM anchor. The existing fallback remains in force for rejected/noncanonical
primaries and accepted reference-bulge primaries. Candidates are still subject
to the existing mismatch, gap count, gap size, match fraction and PAM checks.
This is bounded local expansion around Columba candidates; it is not an
exhaustive search guarantee.

For each guide, contig and strand, reporting selects one representative per
PAM-adjacent boundary: reference end on `+`, reference start on `-`. Distinct
PAM positions survive even when protospacers overlap; the same anchor is
reported once. Among same-anchor alignments it prefers the larger WFA2 score
(scores are non-positive) and uses position/end/CIGAR for deterministic ties.
Strict NGG eligibility precedes selection. CIGAR orientation, coordinate
conventions, and frozen parameters are unchanged. Alternate alignments at one
anchor can still be omitted by design.

## Independent before/after comparison

Every panel was rebuilt into a fresh small Columba index, searched with e=5 and
separate e=0 recovery, and independently enumerated using the same oracle and
criteria as the historical report. Results below are PAM-anchored site recall;
precision is the fraction of reported anchors present in oracle truth.

| Panel | Before recall | Corrected recall | Corrected precision |
|---|---:|---:|---:|
| Synthetic | 28/33 (84.85%) | 33/33 (100%) | 33/33 (100%) |
| CHM13 windows | 2/3 (66.67%) | 3/3 (100%) | 3/3 (100%) |
| Twelve-assembly windows | 2/4 (50%) | 4/4 (100%) | 4/4 (100%) |

Strict interval **and category** recall is intentionally different: one
representative per PAM anchor need not emit every valid alternate interval or
edit placement. The original and corrected synthetic counts are:

| Category | Before | Corrected | Corrected precision |
|---|---:|---:|---:|
| Exact | 4/4 | 4/4 | 4/4 |
| Mismatch only | 8/8 | 8/8 | 8/8 |
| Guide insertion / RNA bulge | 4/8 | 4/8 | 4/4 |
| Reference deletion / DNA bulge | 4/9 | 5/9 | 5/5 |
| Mismatch + bulge | 8/45 | 12/45 | 12/12 |

The strict category denominators overlap where different valid alignments exist
at the same interval. The real-window category counts, every remaining strict
interval discrepancy, candidate SAM and corrected PAF are retained in the
panel subdirectories. No PAM-anchored oracle site remains missing, and no
unexpected PAM anchor was reported in these fixtures. The independent
emitted-row auditor passed every corrected row: 33 synthetic, three CHM13,
and four twelve-assembly-window rows, with no errors.

## Tests, runtime and provenance

`tests/pam_anchor_coverage.rs` uses retained independent fixture references
and real Columba SAMs to require all seven distinct missing keys (eight panel
records), both strands, strict PAMs and no duplicate PAM anchors. Reporting
unit tests check same-anchor deduplication, best-score selection and preservation
of different anchors on both strands. Existing reverse verification tests now
assert that the original exact hit remains present rather than that no other
valid PAM anchor is recovered. All 152 standard Rust tests passed; the one
normally ignored live CFD test also passed after fetching its pinned source and
pickles. All 43 Python tests passed across release-validation and both mismatch
benchmark suites. The auditor and tests do not establish population sensitivity.

A 10-repeat alternating-order comparison of the original and corrected
**optimized** executables on these tiny window indexes, including both Columba
calls, yielded median elapsed times:

| Panel | Original | Corrected | Ratio |
|---|---:|---:|---:|
| Synthetic | 22.37 ms | 22.76 ms | 1.017× |
| CHM13 windows | 16.45 ms | 16.58 ms | 1.008× |
| Twelve-assembly windows | 43.47 ms | 45.27 ms | 1.041× |

These near-startup-bound measurements are not a genome-scale performance
estimate; extra PAM screening and WFA verifications may cost substantially more
on dense genomic candidates. Samples are in `benchmark.json`; the runner is
`scripts/release_validation/benchmark_pam_anchors.py`.

The corrected optimized executable is isolated under
`/export/local/home/salehi/crisprapido-pam-corrected-20261030-target/release/`
with SHA256 `8a207b61a3cc5ad2b5f208115a19dbf3eaaf76b7064519aab0340a154aa9cbd8`.
The original candidate executable still hashes to
`f59b1bbf58e70b321359c7bd7b8687aad81bb874b22dc2f083319453726bc876`.
Rebuilt indexes and test artifacts remain outside Git in
`/export/local/home/salehi/crisprapido-pam-corrected-validation-20261030-v2/`.
Small PAF/SAM, metrics, trace and audit summaries are retained here; historical
outputs were not overwritten. No merge or production-binary replacement occurred.

## Later broader validation

The subsequent [broader coverage report](../broader_coverage_1eb2f40_2026-10-30/REPORT.md)
finds additional missed PAM-anchored sites in repetitive and gap-alignment
fixtures. The historical results above remain unchanged; the preflight gate
continues to be blocked.
