# Development defaults and bounded validation

This separate development branch starts at validated commit `2c3f544` on
`fix/pam-cfd-validation`. It changes the CLI defaults to `-m 5 -b 1 -z 2`;
all remain configurable. PAM handling and `-f 0.75` are unchanged. For a
20-base fully specified guide the fraction requires at least 15 matching
bases, so six or seven mismatches cannot pass at the default fraction.
The `-m 6` and `-m 7` tests below explicitly use `-f 0.70` and `-f 0.65`.

## WFA selection defect and correction

The previous eligible-alignment fallback ran only when `-b 1 -z 1`. A
regression with guide `TGCCACCAATACCCTTTGGA` and target
`TGCCACCAATCCCTTTGGAA` demonstrates an eligible four-mismatch alignment
rejected at `-m 5 -b 1 -z 2 -f 0.75`: WFA selected a CIGAR that violated
the configured limits and the fallback was skipped. The new branch retains
the optimized zero/one-base path and, when WFA's selected alignment fails,
enumerates the bounded possible placements of a single contiguous longer
gap. It applies the existing mismatch, gap, fraction and affine scoring rules.
No eligibility threshold, PAM rule, or CFD calculation was relaxed. The
frozen executable and prior validation artifacts were not modified.

## Independent small indexed synthetic panel

`scripts/release_validation/validate_new_defaults.py` deterministically
constructs 32 small contigs using two guides, both strands, exact sites,
5/6/7 mismatches, two-base insertions and deletions, and two-base gaps plus
mismatches. An independent direct string enumerator checks every possible
single gap of size 0, 1 or 2 against strict NGG PAM anchors. Actual indexed
Columba searches use its `e = m + b*z` bound and WFA verification. Output
CIGAR, reference and guide consumption, mismatch count, fraction, gap event
and length, strand, coordinates and PAM are independently checked for every
emitted record. This check does not independently recompute CFD scores.

| Settings | Columba edit bound | Oracle anchors recovered | Emitted alignment audit | Missed / unexpected / duplicates |
|---|---:|---:|---:|---:|
| `-m 5 -b 1 -z 2 -f 0.75` | 7 | 31/31 | 31/31 | 0 / 0 / 0 |
| `-m 6 -b 1 -z 2 -f 0.70` | 8 | 36/36 | 36/36 | 0 / 0 / 0 |
| `-m 7 -b 1 -z 2 -f 0.65` | 9 | 41/41 | 41/41 | 0 / 0 / 0 |

The command with default values produced the same row multiset as explicit
`-m 5 -b 1 -z 2 -f 0.75`. Independent eligible planted PAM anchors recovered
were 26/26, 29/29 and 32/32 respectively. Other oracle anchors arise from
additional valid overlapping PAMs in the deterministic fixture. At `-m 5`,
some planted six/seven-mismatch sequences are eligible at an alternative
alignment; the eligibility denominator is determined by the oracle rather
than the intended construction label. `summary.json` records the counts.
Separate direct Columba probes also exited successfully at bounds 7, 8 and 9.
The focused tests cover every position of a two-base guide or reference gap,
and the historical ineligible-WFA example above.

## Tests and limits

All 160 regular Rust tests, the separately run pinned CFD test, and all 52
Python tests passed. The existing frozen-parameter regressions remain passing.
The first version of the new fallback regression failed before the correction,
and passes afterward.

These are small, synthetic bounded tests, not proof of genome-wide discovery
or sensitivity. Columba's all-candidate completeness at these higher bounds,
repetitive regions, large-genome runtime and memory, and more than one gap
event require further independent validation. The alternative enumeration is
restricted to `-b 1`; for higher values WFA's one selected alignment can still
miss a different eligible alignment. The original independent-review and
preflight requirements remain open. The previously validated four-mismatch,
one-base-bulge resource results must not be reused as performance estimates
for the new defaults.

## Subsequent real-reference validation

The [bounded CHM13 and GRCh38 panel](../real_defaults_validation_49f20fc_2026-11-02_v1/REPORT.md)
recovered 204/204 independent PAM anchors, including the requested mismatch and
one- and two-base bulge classes on both strands. These results do not establish
whole-genome sensitivity or resource requirements for the new defaults.
