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
