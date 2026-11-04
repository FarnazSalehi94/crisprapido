# Reference N correction in PAM-anchored Columba prototype

The research prototype now recovers all 172 anchors in the independent engineered-N panel, with no extras or duplicates, without WFA2. The correction changes only the standalone prototype and its tests. It does not change CRISPRapido defaults, its production binary or PAM semantics.

## Root cause and semantics

Columba VANILLA 2.0.3's `buildindex.cpp` `replaceNonACGT` replaces input `N` by an A/C/G/T base before indexing; its indexed alphabet does not preserve an N symbol or its identity. The prototype's `walk` explores those four indexed letters and calculates provisional mismatches against the replacements, not the original reference. Previously, `emit` rejected any span containing a reference N. Simply accepting the indexed CIGAR would be wrong: a replacement that happens to equal the guide is represented as `=`, though the original reference N is a mismatch for a resolved guide; replacement-dependent counts and match fraction would also be wrong. On minus strand the original target must be reverse complemented, with N complementing itself.

For a resolved A/C/G/T guide, the existing verifier's `eligible_one_gap_alignment` / `validate_cigar_accounting` and the independent direct oracle treat aligned reference N as an `X`: one mismatch and no match. N consumed by a reference-only `D` is part of the deletion event and does not itself consume a mismatch. A guide-only `I` consumes no reference base. N is not assumed to equal an A, C, G or T. The prototype still rejects guide Ns and other ambiguous target characters; it supports one contiguous gap event only.

A pattern PAM `NGG` uses N as a **wildcard over resolved A/C/G/T** in the reference. The existing `extract_adjacent_pam` requires all three reference PAM bases to be resolved. Therefore reference `NGG`, `ANG` and `AGN` are all ineligible: an actual reference N is not a valid wildcard base or required G. Both strands require the complete correctly oriented PAM. This is distinct from a target-sequence reference N.

## Correction and coverage argument

The FM-index search is unchanged. After locating each indexed alignment, `emit` now fetches the original target, reverses/complements it for minus strand, replays the existing gap structure, and recalculates each `=`/`X`, the mismatch count and match fraction from the true reference. It enforces the configured thresholds and only then selects one deterministic representative per guide/contig/strand/PAM boundary. Different occurrences in the same suffix-array interval are replayed independently. Reference N inside D remains D, with no mismatch. Strict PAM validation remains in place before replay.

For any eligible alignment under this restricted model, a resolved guide paired with a reference N has one true mismatch. Columba's substitute letter either matches the guide (indexed cost 0) or differs (indexed cost 1). Thus each indexed partial path has **no more** mismatches and **no fewer** matches than the true path. A mismatch-limit or optimistic maximum-match pruning condition cannot remove an otherwise eligible true path. The same gap transitions are explored, including Ns in a deletion, and the four PAM seeds cover every valid resolved PAM wildcard. Consequently the real path has an indexed representative; replay filters replacement-induced false positives. This argument depends on the index and reference agreeing at resolved bases, on search traversal and location correctness, on resolved guides and this exact counting model. Tests support it but are not a formal proof of all configurations. The separate index checker in the prior report can check names, offsets and resolved text but is **not enforced by the prototype**. Columba replaces all non-ACGT input, while the prototype intentionally supports only A/C/G/T/N target text; other IUPAC symbols are not covered.

## Independent tests

| Check | Result |
|---|---|
| 172-anchor real-derived/engineered-N panel, strict m5/b1/z2/f.75 | **172/172 recovered**, 0 unexpected, 0 duplicate anchors; baseline also 172/172 |
| Prior original 8 kb panel | **204/204**, 0 unexpected or duplicates |
| Difficult PAM-anchor panels | Synthetic 33/33; CHM13 3/3; 12-window 4/4; no extras or duplicates |
| Existing alternative-alignment examples | 3/3; no CIGAR identity requirement |
| Resolved PAM wildcard bases on both strands | 8/8 |
| Deterministic randomized resolved-reference set | 29/29, no extras/duplicates |
| Deterministic randomized reference-N set | 17/17, no extras/duplicates |
| Targeted N/PAM/boundary regression | 18 oriented cases, 12 valid anchors recovered, six invalid PAM cases rejected, 0 extras/duplicates |

[Targeted test](test_reference_n.py) uses independent exhaustive one-gap enumeration and independently replays every emitted CIGAR, count and eligibility. It includes N within the target, immediately adjacent to PAM, near guide insertion and reference deletion, Ns inside the deletion, consecutive Ns, and both strands. Targets begin at a contig boundary or end at one on the reverse strand. It distinguishes an actual N at the wildcard PAM position and at either G; all are rejected on both strands. A separate zero-mismatch run checks that an aligned N fails but an N consumed in a deletion remains eligible. An additional guide is constructed from the actual indexed replacement at an N: even when the replacement equals its guide base, the true alignment is rejected at m=0 and reported with one X at m=1. The targeted index was built and tested four times; Columba can choose different replacement letters. The randomized N test adds 1–3 target Ns per small contig. The independent 172-panel audit replays all emitted rows and compares exact anchors; the fixed output and comparison are in [n_correction_results](n_correction_results/), while the input and oracle remain in [results](results/). The previous 75 misses are all recovered. Passing these finite tests is not a genome-wide or clinical sensitivity claim.

## Bounded performance

The same-host, same-input, single-thread repeated runs use prebuilt indexes and strict m5/b1/z2/f.75. The baseline is the same separately instrumented existing implementation as the [previous validation](REPORT.md), not the production binary. The prototype's state count counts `walk` entries; its search timer includes indexed traversal, locating and FASTA replay, not input loading or final TSV writing. The instrumented baseline includes Columba discovery, SAM I/O, WFA2, output and CFD. Peak RSS is the maximum individual process value, not aggregate parent/child RSS. These are not equivalent workloads; do not interpret their ratio as a search speedup. Index building is excluded. Measurements are retained as JSON under [n_correction_results](n_correction_results/).

| Reference | Guides | Prototype wall s | Search wall s | Search states | Prototype RSS KiB | Baseline wall s | Baseline RSS KiB | Anchors each |
|---|---:|---|---|---:|---|---|---|---:|
| Engineered-N 8 kb panel | 64 | 0.25, 0.25 | 0.252, 0.258 | 2,930,404 | 4,304, 4,284 | 11.52, 11.48 | 22,512, 22,524 | 172 |
| Original resolved 8 kb panel | 96 | 0.70, 0.72 | 0.705, 0.719 | 8,313,485 | 4,308, 4,380 | 13.75, 13.75 | 30,732, 30,720 | 204 |
| CHM13 chr1 2 Mb | 1 | 0.65, 0.67 | 0.635, 0.660 | 8,406,922 | 13,116, 13,080 | 29.71, 29.31 | 66,888, 67,280 | 10 |
| GRCh38 chr9 2 Mb | 1 | 0.60, 0.60 | 0.582, 0.584 | 7,616,819 | 14,252, 13,468 | 30.05, 29.70 | 67,944, 67,156 | 23 |

The prototype and baseline report identical anchor sets for each run; the emitted CIGARs from both 2 Mb windows were independently replayed in both repetitions and are not exhaustively oracle-enumerated. Neither measurement establishes whole-genome scaling. The replay change does not add a genome scan, but N-heavy and highly repetitive reference behavior needs more testing. No CFD-guided pruning was introduced.

## Reproduction and limitations

Build using `research/pam_anchored_poc/build.sh` and the recorded vanilla Columba checkout. Run the existing 172-panel FASTA and guides against its prebuilt or newly built index; write output as `poc.tsv` beside the retained `oracle.jsonl`, `baseline.paf`, and `index-check.json`, then run `python3 research/pam_anchored_poc/validation/compare.py PANEL_DIR`. Run `test_reference_n.py BUILDER_WRAPPER POC_WRAPPER SCRATCH_DIR` and `randomized.py BUILDER_WRAPPER POC_WRAPPER SCRATCH_DIR --with-n`, plus the existing original-panel and fixture scripts. `benchmark.py --help` documents repeated comparisons. The tested POC binary SHA256 is `ac634f3937df3cb234bd06d9c86c2f5e63f8e608575a9af3fd111862c02a65c1` (the later source-only edit names an unused parameter and adds a comment).

The POC still lacks general PAM/noncanonical modes, multiple gaps, ambiguous guides, general IUPAC target semantics, automatic index/reference binding, production scoring and representative compatibility, and broader exhaustive tests. It remains a research implementation and is not ready for production integration. No frozen branch, main branch, production binary or defaults were modified; no whole-genome or pangenome job and no CRISPR-Bulge integration was run.
