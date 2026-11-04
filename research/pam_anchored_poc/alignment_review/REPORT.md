# Representative alignment and native safety review

The 54 differing representatives at the 204 expected PAM anchors are all independently eligible on both sides. They are not interchangeable for scoring: **41/54 pairs have different CFD values**. There is no evidence here for a universal biological preference, so the prototype's existing deterministic selection is unchanged and experimental integration remains a separate decision. The native research executable now requires a reference/index binding even when invoked directly; it fails closed if verification cannot run.

## Independent review of all 54 pairs

[compare_representatives.py](compare_representatives.py) joins by guide, contig, strand and exact PAM-adjacent boundary, then independently extracts each target and PAM from the original FASTA and replays each CIGAR. It checks coordinates, strand, full guide/reference consumption, actual `=`/`X` bases, mismatch count, contiguous gap-event count and combined event size, strict PAM and 15/20 matches. All 54 baseline and all 54 prototype alignments pass. It calculates CFD using the separate matrix oracle in `scripts/release_validation/audit_emitted.py`, and compares all **108** scores with the actual Rust `get_cfd_score`; all agree. The 54 baseline reported four-decimal tags agree within their display precision. [representatives.json](results/representatives.json) records each anchor, both spans/CIGARs, counts, oriented PAM, gap positions, match fraction, scores and comparison flags.

| Difference among the 54 pairs | Anchors |
|---|---:|
| Both alignments independently eligible | 54 |
| Mismatch count differs | 5 |
| Number of contiguous bulge events differs | 5 |
| Maximum gap-event size differs | 5 |
| Target start/end coordinates differ | 5 |
| Gap operation/position differs | 54 |
| CFD differs by more than 1e-9 | 41 |

All 108 alignments in this resolved-reference panel have defined CFD values. For example, at `CHM13_nonrepeat_candidate_minus_D1` / `CHM13_nonrepeat_candidate` / minus boundary 903, baseline `1X8=1X1=2D4=1X3=1X` scores approximately 0.0020202 and prototype `1X7=2D1=1X5=1X3=1X` scores 0.0. A CFD difference is not evidence that either alignment is the biological truth; CFD is not a calibrated cleavage probability.

As an **alignment-model comparison**, independently recomputing the current affine cost `3*mismatches + sum(5 + run_length)` over each I/D run favors the baseline in **4** pairs, the prototype in **1**, and ties in **49**. In the prototype-favored case (`GRCh38_repeat_candidate_plus_D1`, `CHM13_repeat_candidate`, minus boundary 428), baseline has 5 mismatches and no gap (cost 15); the prototype has 2 mismatches and 1D (cost 12). This is not proof that either software emitted the global best alignment for every possible path, nor a biological preference. The case shows why assuming WFA's selected record is necessarily best among all eligible anchor alternatives would be unsafe.

A consistent proposed policy, if product compatibility is desired, is: enumerate eligible alignments at an anchor, minimize a **declared** affine cost, then break ties by `(start, end, expanded guide-oriented CIGAR)` in a documented order. The current pipeline uses its obtained WFA score then `(start,end,cigar)`; the prototype currently minimizes `(mismatches + inserted bases + deleted bases, mismatches, compressed CIGAR, start)`. These are deterministic but not equivalent. We have not changed the prototype to the proposed policy: the existing pipeline may not enumerate every eligible alternative, mixed I/D penalty interpretation and PAF scoring need a decision, and changing representative choice can change CFD. A future backend must select and document an objective rather than maximize CFD or claim one pair member is biologically correct. No CFD-guided pruning was added.

## Mixed-gap completeness

Production `validate_cigar_accounting` treats adjacent I and D columns as one contiguous event and applies `-z` to their combined length. The existing two-strand m=0/b=1/z=2 fixture has an eligible adjacent I/D alignment at each anchor; the prototype reports both, whereas the WFA-first pipeline reports neither. The [fixture and its outputs](../compatibility/results/mixed_gap/) document that difference; neither implementation was treated as an oracle.

[mixed_exhaustive.py](mixed_exhaustive.py) additionally enumerates all A/C/N targets of lengths 1–5, their reverse complements, a four-base guide and every full end-to-end alignment path with at most one contiguous event of total size up to two. It compares all PAM anchors with an independently recursive direct-string oracle under six mismatch/bulge/fraction settings. Across **726 contigs**, the configurations had 8, 150, 212, 556, 110 and 918 expected anchors respectively; every expected anchor was recovered, with zero unexpected or duplicate anchors. These deliberately small exhaustive tests and the previous 20-base panels found **no missed eligible prototype anchor**. They do not prove completeness for arbitrary guide lengths, larger gaps or multi-event settings.

## Mandatory reference/index binding

The built native prototype now calls the build-pinned Python interpreter and [index_manifest.py](../compatibility/index_manifest.py) via `fork`/`exec` before opening its FM index or writing output. Direct invocation and the checked wrapper both require the build-time binding: a missing manifest, changed FASTA including N positions, changed index component, failed/missing verifier or a mismatched index fails closed. The interpreter and helper paths are embedded by [build.sh](../build.sh); PATH substitution cannot replace the interpreter. [test_binding.py](../compatibility/test_binding.py) checks direct invocation as well as the wrapper, including a fake `python3` on PATH. The separately trusted builder must create the manifest from the actual original FASTA and index; the earlier unbound indexes must be rebuilt. The helper and interpreter must remain available at the embedded locations.

This prevents the documented **direct-executable bypass**, not adversarial modification of the executable, helper or binding. The manifest is not signed and verification still has a time-of-check/time-of-use window if files change after hashing. Use trusted, immutable inputs; an atomic file-descriptor-based scheme, versioned index format and scalability of full-file hashing remain integration work. Columba's replacement of N means an unbound old index cannot establish its original N positions from index text alone.

## CLI and rerun

The prototype remains research-only. `checked_search.py` exposes -p, -m, -b, -z and -f with defaults GG/5/1/2/0.75. All 16 resolved two-base PAM suffixes, both strands, m=0/2/5/20, b=0/1, z=0/1/2 and f=0/.75 passed **96 independent option/PAM runs**. Negative/malformed limits, invalid PAMs, b>1, z>2, guide N, non-ACGTN reference symbols, duplicate/empty guides and unsupported wrapper flags fail rather than silently narrowing output. The existing CRISPRapido CLI additionally supports more gap events/sizes, guide ambiguity, noncanonical PAM policy, `--no-filter`, window/threads and imported SAM; these are not implemented by this POC and must remain available through the unchanged backend. Broader f32 thresholds and random long/repetitive references need testing.

After native enforcement, the rerun recovered **204/204** and **172/172** independent anchors with no unexpected or duplicates. The difficult PAM panels recovered 33/33, 3/3 and 4/4; three alternative-anchor examples, all eight wildcard/strand cases, reference-N boundary tests, randomized resolved and N panels, four single-N substitutions and all 64 triple-N substitution indexes passed. Independent CFD tests on 1,024 option-panel alignments covered all PAM suffixes, both strands and all operations: 896 scores agreed and 128 agreed as unavailable. The 172-panel check agreed on 97 scores and 75 unavailable N scores. Historical unbound fixture builders were changed to create build-time bindings for their test indexes; the deliberately substituted-N test builders bind the original N FASTA while passing a controlled replacement copy to the trusted test builder.

The native reference check adds whole-file hashing on every invocation; no performance claim is made for this change or for genome-scale indexes. No production binary, default, main branch, backend integration, whole-genome run or CRISPR-Bulge work was changed or performed. Representative/CFD policy, atomic binding and unsupported modes remain blockers.
