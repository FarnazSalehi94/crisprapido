# Independent review package — review not yet performed

The purpose is to let a reviewer who did not develop the correction examine
the release evidence. **No independent review has been performed or signed.**
This package is not an approval request to run the full experiment.

- Implementation under test: commit `0aa9944e2eae4b7f8286415ee5c3fb3667cf6097`;
  binary and dependency hashes and frozen settings are in `freeze_manifest.json`.
- Search path: Columba e=5 plus e=0, bounded candidate-overlap expansion,
  anchored WFA2, conditional one-single-gap eligibility enumeration, strict
  three-base NGG and one representative per guide/strand/genomic PAM anchor.
  See `src/columba.rs`, `src/verification.rs` and `src/reporting.rs`.
- Independent oracle: direct character comparisons of ungapped and every
  possible one-base I/D placement on both strands, excluding ambiguous bases,
  at the same mismatch, gap and match-fraction limits. PAM key is guide,
  contig, strand and exact PAM interval; same PAM does not mean identical
  alignment. See `scripts/release_validation/coverage_oracle.py` and
  `resolve_bulge_discrepancies.py`.
- Evidence: `heldout_synthetic/`, `heldout_real/`, the prior
  `docs/performance_validation_a817ab5_2026-10-31_v1/`, this report's
  `scaling/per_guide.tsv`, per-run profile counts and independent audit
  summaries. Large index/reference/raw SAM artifacts are outside Git in the
  documented versioned run directory.
- Main results: synthetic 212/212 and real-window 18/18 PAM keys, precision
  100% in both; no output-invalid rows in the 230 held-out or 43 regional
  pilot outputs. Small sample results are not genome-wide sensitivity.
- Critical uncertainty: no proof Columba returns an overlapping candidate for
  every eligible locus. Different equal-PAM alignments are deliberately
  collapsed to a representative, and gap fallback completeness is restricted
  to b=1,z=1; non-ACGT and whole-genome regimes need separate treatment.
- Performance: 16 individual searches over four selected 2 Mb windows and
  four guides each. Time and memory cannot be directly projected to batched
  whole-assembly jobs; repeats and high candidate density may dominate.

Suggested independent checks: rederive the overlap bound; inspect handling of
reverse-strand and terminal gaps, Ns and float thresholds; challenge PAM-key
identity and source-coordinate assumptions; recompute a stratified oracle and
inspect Columba candidate absence; independently verify all 16 per-guide
profiles and audit summaries; assess whether the scale envelope and old full
plan remain appropriate. Record reviewer identity, exact reviewed revision,
method, findings and disposition before marking this requirement complete.

## Controlled whole-assembly pilot

A separately versioned [pilot report](../whole_assembly_pilot_0aa9944_2026-11-02_v1/REPORT.md)
records the frozen-source verification and four completed serialized Slurm
jobs. The independent exact truth was 4/4 and all 13,067 emitted rows passed
the auditor; individual searches took 320.70–450.58 seconds, and observed
maximum process RSS was 11.1–11.2 GiB. The search-only 1,200-process
extrapolation is roughly 121 serial hours, with substantial uncertainty and
without auditing or index construction. This does not constitute independent
review. Reviewers should also
examine the diagnostic replay equivalence, full-row audit, timing and memory
measurement limits, and non-random guide selection.
