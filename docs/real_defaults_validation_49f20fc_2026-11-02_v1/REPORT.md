# Real human reference validation of m5/b1/z2 defaults

**Result:** all 204 independently enumerated PAM-anchored sites were recovered
(204/204 recall); all 204 reported anchors were expected (204/204 precision).
There were no missed, unexpected or duplicate anchors. All 204 emitted rows
passed the independent sequence, CIGAR, PAM and eligibility audit. This is a
small, deliberately constructed guide panel over real sequence, not a
whole-genome sensitivity estimate.

## Inputs and method

The development implementation was commit `49f20fc`, compiled separately
with Rust 1.85.1 in release mode, SHA256
`cb10e5d6be970bb6ef90aac4f2fe5ab035ae2b4d4ce643f28ec32f3895d7e42d`.
The original frozen executable and branch were not changed. Reference input
was fetched from the previously verified indexed CHM13 and GRCh38 FASTAs with
samtools; each extracted window was checked byte-for-byte against the earlier
2 Mb reference copy. Paths, exact coordinates and sequence SHA256 values are
in `panel.json`; `ref.fa` contains the four small, unchanged extracts.

| Assembly | Region, 1-based inclusive | Repeated 20-mer samples / 100 |
|---|---|---:|
| CHM13 | chr1:100004001–100006000 | 0 |
| CHM13 | chr9:31314001–31316000 | 54 |
| GRCh38 | chr1:100000001–100002000 | 0 |
| GRCh38 | chr9:31750001–31752000 | 29 |

Window selection uses the existing 2 Mb chr1/chr9 panels. Within each parent
region, nonoverlapping 20-mers were counted; the 2 kb block with the smallest
or largest number of repeated samples was chosen. This is a repeat-content
proxy within each parent window, not a formal repeat annotation or a random
sample. The coordinates are not homologous pairings between assemblies.

At one strict NGG PAM per strand in each window, 12 distinct 20-base guide
sequences were derived from the unmodified reference, yielding 96 unique
guides. Designs include exact, each of 1–5 substitutions, one- and two-base
guide insertions and reference insertions, and two-mismatch plus two-base-gap
combinations. The direct oracle confirms all 96 intended edit structures at
their source PAM anchors. Guides are selected around known sites, so these
results do not represent randomly selected guides. Every guide is searched
against every selected window, including the two repeat-rich windows.

The independent oracle scans both strands and every complete strict-NGG PAM
within each 2 kb extract. For target spans of 18–22 bases it directly
compares each ungapped sequence or each possible placement of one contiguous
one- or two-base gap, enforcing at most five mismatches and at least 15 guide
matches. It does not use Columba, WFA2 or CRISPRapido candidate expansion.
Results are keyed by guide ID, contig, strand and precise genomic PAM anchor;
overlap alone never establishes equivalence. Window-boundary-crossing sites
are outside the panel. `oracle.jsonl` records the expected anchors and all
observed eligible edit classes; `actual.paf` is the indexed Columba search
output and `audit.jsonl` records the independent emitted-row checks.

## Results

| Oracle class | Expected anchors | Recovered |
|---|---:|---:|
| Exact | 9 | 9 |
| 1 / 2 / 3 / 4 / 5 mismatches, ungapped | 12 / 13 / 13 / 19 / 31 | 12 / 13 / 13 / 19 / 31 |
| One-base guide insertion, no mismatch | 25 | 25 |
| Two-base guide insertion, no mismatch | 29 | 29 |
| One-base reference insertion, no mismatch | 19 | 19 |
| Two-base reference insertion, no mismatch | 19 | 19 |
| Any eligible mismatch plus bulge | 204 | 204 |

**Classes overlap:** one PAM anchor can have several eligible alignments, and
the program selects a representative. The table counts anchors having at
least one alignment in each class; it does not claim the representative has
each listed CIGAR. Of the 204 actual representative alignments, 22 have five
mismatches, 37 contain `2I`, and 28 contain `2D`. All 96 designed anchors
were recovered. Expected/recovered counts by window and strand were CHM13
chr1 +18/18 and −15/15, CHM13 chr9 +30/30 and −37/37, GRCh38 chr1 +24/24
and −12/12, and GRCh38 chr9 +36/36 and −32/32.

The emitted-row audit reconstructs CIGARs against the guide and reference,
checks coordinates, strand, PAM, mismatches, gap count and size, and the match
fraction. It does not independently recompute CFD. No discrepancies were
observed, so no correctness change was warranted. A regression test verifies
reference-window hashes, the 96 intended edit structures and expected versus
reported PAM keys.

The one-thread development search used `-p GG -m 5 -b 1 -z 2 -f 0.75`,
Columba e=7 plus e=0, and the existing scoring matrices. The e=7 and e=0 SAMs
contained 6,284 and 9 mapped rows respectively. GNU time recorded 13.44 s
wall time and 35,520 KiB maximum process RSS for the complete batched indexed
search; reference extraction, index construction and the independent oracle
are excluded. Oracle enumeration took 11.28 s separately. Times are from a
shared machine and must not be extrapolated to large references or other guide
distributions. The index, raw SAMs, logs and full run are retained at
`/export/local/home/salehi/crisprapido-real-defaults-49f20fc-20261102-v2`.

## Verification and limitations

The full 160-test Rust suite, pinned standalone CFD comparison and all 54
Python tests passed. No production binary, historical result or main branch
was modified; the 100-guide × 12-genome search was not launched.

This panel covers only 8 kb and 96 selected guides derived from eight real
PAM anchors. Repetitive regions are represented but not comprehensively
sampled. The test does not prove Columba always discovers eligible sites,
does not test PAM anchors crossing window boundaries, and does not measure
whole-genome runtime or memory. Higher gap-event counts are outside this
configuration. Independent human review and the existing release preflight
requirements remain outstanding; previous frozen-parameter performance and
sensitivity results cannot be treated as validation of these new defaults.

## Subsequent whole-assembly attempt

The [controlled whole-assembly pilot](../whole_defaults_pilot_60f9d6c_2026-11-02_v1/REPORT.md)
was stopped after CHM13–ELMOD2 Columba e=7 discovery produced 17,653,389
matches, approximately 70.9 times the previous e=5 count. No final output was
available to audit or compare. The other guide–assembly pairs were not started;
this leaves the whole-assembly correctness and resource evaluation incomplete.
