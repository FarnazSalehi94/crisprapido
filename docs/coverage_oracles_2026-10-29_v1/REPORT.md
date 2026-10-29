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
