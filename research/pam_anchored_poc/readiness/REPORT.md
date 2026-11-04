# Experimental-backend readiness assessment

**Recommendation: retain the prototype as research code; do not integrate it yet.** The N-substitution argument survived the new bounded exhaustive tests, and no missing prototype site was found. The evidence supports a future explicitly selected experimental backend after the interface and compatibility blockers below are addressed. It does not establish completeness for every reference or configuration. The existing Columba + WFA2 backend must remain selectable and continue to be the default.

## Reference-N search completeness

For a resolved guide base aligned to reference `N`, the validation model charges a mismatch. Columba substitutes one A/C/G/T base during index construction. In the indexed shadow of any eligible alignment, every substituted N has a provisional mismatch cost of either zero or one, never greater than the true cost of one; its provisional match contribution is never less than the true contribution of zero. Thus the prototype's mismatch pruning and optimistic maximum-match pruning cannot remove that alignment merely because of the substitution. An N inside a reference-only deletion has no mismatch cost in either representation. Exact A/C/G/T PAM positions, strand conversion, gap structure and suffix-array location must still be correct; FASTA replay recomputes the true CIGAR and scores. This is a conditional argument, not a formal proof: it assumes a correct index/reference pair, resolved guide, current one-gap model, complete underlying FM traversal, and no resource exhaustion.

We tested that argument independently with [substitution_oracle.py](substitution_oracle.py). Its synthetic references include a complex guide and a repetitive `AC` guide; seven designed structures per guide on both strands give **28 designed cases**. The cases include five total mismatches including N; five mismatches plus 1D or 2D; four plus 1I; three plus 2I; three consecutive Ns; and two Ns consumed by 2D. The test constructs independently enumerable PAM-anchored alignments, checks each intended class, and includes adjacent AGG/GGG anchors. It builds four separate indexes by replacing every N with A, C, G or T while using the unchanged original N-containing FASTA for replay. The exhaustive oracle finds **40 anchors**. An independent audit reconstructs every reported alignment, coordinates, mismatch and gap counts, and the 15/20 match threshold; it compares PAM-start keys and checks duplicates.

| Indexed substitute | Independent expected | Prototype recovered | Prototype unexpected / duplicates | Existing pipeline recovered |
|---|---:|---:|---:|---:|
| A | 40 | 40 | 0 / 0 | 37 |
| C | 40 | 40 | 0 / 0 | 38 |
| G | 40 | 40 | 0 / 0 | 37 |
| T | 40 | 40 | 0 / 0 | 37 |

To test combinations rather than just a uniform substitution, [multi_replacements.py](multi_replacements.py) enumerates **all 64 A/C/G/T triples** for three consecutive reference Ns, on both strands in complex and repetitive synthetic contexts. Every index gave **8/8 oracle anchors**, no extras or duplicates. These are short constructed references; they do not prove arbitrary repeat or pangenome behavior. The existing real-derived, repeat-enriched 8 kb N panel was rerun: **172/172** prototype and baseline anchors, no extras/duplicates, with every emitted record audited. The independent original resolved panel was also rerun: **204/204** prototype anchors, no extras/duplicates. Details and exact keys are in [results](results/).

The baseline discrepancy above is real for these constructed inputs, not a prototype false positive: each baseline-missing site has an independently enumerated eligible full-guide N alignment. For each missed guide/contig, a separate inspection of both Columba e=7 and e=0 SAMs found **no mapped candidate anywhere on that contig**. Independently aligning the indexed substituted target with up to one gap gives edit distance 6 or 7 at the intended PAM, within the requested e=7 candidate bound. Thus the immediate observable failure is before WFA2 verification for these runs. The precise Columba search/overlap/reporting mechanism was not isolated; these short multi-contig fixtures cannot establish its frequency in genomes or attribute all misses to N substitution. Raw SAMs, exact misses and indexed distances are retained. No existing-pipeline fix was attempted. The earlier conclusion that 172/172 held on one panel remains valid for that panel; baseline agreement is not a substitute for an independent oracle.

## Bounded performance and phase attribution

Two repetitions on the same host, one search thread, prebuilt verified indexes, same guides, references and m5/b1/z2/f.75. The baseline is a **disposable instrumented copy** of the existing source; its output on the 172 panel was byte-identical to the uninstrumented baseline. The prototype is a separately instrumented copy whose 172-panel output was byte-identical to the current prototype. Instrumentation and retained SAM add overhead; outputs and summaries are in [results](results/). No production executable was changed. Peak RSS is the maximum individual process value, not aggregate process-tree memory.

| Reference / guides | Prototype total wall s | Prototype traversal excluding locate/replay s | Prototype locate/replay s | States | Prototype peak RSS KiB | Baseline total wall s | Baseline peak RSS KiB | Anchors, each |
|---|---|---|---|---:|---|---|---|---:|
| N-containing 8 kb / 64 | 0.25, 0.25 | 0.252, 0.253 | 0.00217, 0.00217 | 2,930,404 | 4,520, 4,284 | 11.61, 11.59 | 29,460, 29,356 | 172 |
| CHM13 chr1 2 Mb / 1 | 0.66, 0.66 | 0.645, 0.648 | 0.000037, 0.000040 | 8,406,922 | 15,000, 13,292 | 29.50, 29.29 | 70,176, 69,920 | 10 |

For the same runs the baseline phase counters (seconds, two repetitions) were:

| Phase | N-containing 8 kb | CHM13 2 Mb |
|---|---|---|
| Columba mapping-worker elapsed | 0.171, 0.171 | 0.224, 0.225 |
| Complete Columba subprocess wait | 0.190, 0.191 | 0.473, 0.476 |
| SAM read and parse | 0.00398, 0.00408 | 0.00979, 0.00988 |
| WFA alignment calls | 8.702, 8.712 | 22.303, 22.163 |
| Eligible-alignment enumeration | 0.688, 0.680 | 1.699, 1.691 |
| CFD calls | 0.00077, 0.00076 | 0.00004, 0.00004 |
| SAM output bytes, both Columba runs | 665,041 each | 1,395,719 each |

The prototype has no SAM stage, WFA stage or CFD; its locate/replay timer covers suffix-array locating, PAM checks, reference replay and representative insertion. `search_wall - locate/replay` includes DFS and timing overhead, not a hardware-level pure FM-index counter. The baseline mapping-worker timer includes Columba search, internal verification and output-record generation, and may overlap the writer. Subprocess wait includes startup, index loading, SAM generation/writing and synchronization. SAM **reading/parsing** is measured, but SAM **writing time is not isolated**; bytes are recorded. Other timers can overlap and are not an additive breakdown of total wall time. CFD timing measures scoring calls, not a guarantee that N-containing sites receive a score. These phase boundaries are deliberately reported rather than pretending an apples-to-apples indexed-search or I/O measurement exists. Different output contracts, instrumentation and scoring mean the elapsed times do not establish a comparable speedup. The 2 Mb anchors matched and both repetitions' output CIGARs were independently audited, but no exhaustive 2 Mb oracle was run.

## What blocks an optional backend

- **Options and semantics:** prototype hardcodes strict NGG, assumes resolved guides, supports one contiguous gap event, and does not support `--allow-noncanonical-pam` or `--no-filter`. CRISPRapido permits a configurable two-base PAM, mismatch, number and size of bulges, fraction and threads. Any experimental selector must keep all CLI parameters configurable and either implement their semantics or reject unsupported combinations clearly; it must not silently change defaults or discard sites. `N` in a PAM is rejected, and other reference ambiguity codes remain unsupported.
- **Reference/index contract:** the standalone checker is not enforced in the prototype. A backend needs reliable index/reference identity, contig order and names, boundary handling, 32-/64-bit index compatibility, and an explicit failure instead of treating an incomplete/failed search as no hits. Existing checks on original Ns cannot prove how they were replaced in the index.
- **Output and scoring:** `reporting.rs` chooses the WFA-scored representative per exact PAM boundary; the prototype uses edit count and a different tie order. Earlier bounded results had 54 of 204 anchor representatives differ. A backend must provide correct target sequence, strand, CIGAR, spans and score to the reporting layer and decide whether equivalent eligible alternatives should retain existing representative behavior. CFD must remain downstream and optional to search. Current `get_cfd_score` can return no score for N-containing targets; preserve that behavior explicitly and validate against the existing matrix scorer, rather than converting N to a base or zero. CFD is not a calibrated cleavage probability.
- **Scope and operational assurance:** the current 32-bit/sparse-SA-4 build, single-threaded proof of concept, FFI/process interface, licensing and maintenance of Columba primitives, resource limits, timeout/error reporting, performance on N-heavy/repetitive regions, random indel tests, high mismatch/alternative paths, parameter ranges, and representative ties require review. Agree on completeness as PAM-anchor recall versus all alignments. Further bounded multi-reference and independent review are needed before genome-scale use.
- **Baseline discrepancy:** independently investigate the observed missing synthetic candidates and fix or document discovery limitations separately. An opt-in backend that recovers more anchors changes reported results, and must not be presented as bit-for-bit equivalent to the current path.

A future design can add an explicit backend selector while leaving the present Columba + WFA2 path available, and feed verified anchored alignments through the existing PAM policy, representative selection and CFD/reporting layer. That design is not implemented here. Record-scoring equivalence and all supported flags must be established before experimental adoption; keep the existing default path until then.

## Reproduction

`substitution_oracle.py --help` describes the four replacement builds and comparison; `multi_replacements.py --help` reproduces all 64 triples. These use separate scratch indexes and original Ns in the verification FASTA; all resolved reference letters and contig metadata are checked against each index. [instrument.py](instrument.py) creates disposable source copies for timing; it adds the SAM-read guard to a copy produced by the existing profiling script and wraps POC locate/replay. Build against the recorded Columba vanilla checkout and run [benchmark.py](../validation/benchmark.py) with `--keep-columba-logs`; the latter retains SAMs/logs in its scratch output to measure bytes and mapping-worker time. The audited 172 panel, oracle, previous performance data and index identity details are in the [N-correction report](../validation/N_CORRECTION_REPORT.md) and [validation report](../validation/REPORT.md).

No production/default/CLI changes, integration, merge, whole-genome experiment, CFD-guided pruning or CRISPR-Bulge work was performed.
