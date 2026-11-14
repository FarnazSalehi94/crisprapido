# Brunello full-cohort readiness review (2026-11-13) — prepared, NOT launched

**Correction:** see [timing/storage reconciliation](TIMING_RECONCILIATION.md).
178 h is traversal only; ~276 h (~11.5 days) includes audit and scoring.
The compressed authoritative PAF+annotations scenario is ~13 TB, not 5.4 TB;
audit records ~6.4 TB, not 14 TB. The 10,000-guide inputs/commands are now
frozen, but its restart/monitor execution wrapper is still launch-gated.
Statements below are historical and must not be used instead of the correction.

Separate manuscript project. Library: 77,441 unique guides, frozen FASTA
SHA256 re-verified `a4bac2dc4c8f5ad0c37713751835b416d1669f899de1cd05ca1a0e828be48971`.
Provisional configuration **m4/b1/z1/f0.75/GG** (the validated candidate-A
from the expanded benchmark: 100 % both-score coverage, perfect nesting,
~100 additional CFD≥0.20 sites per guide vs m3). SBDS rebuilding/searches are
the current priority; Brunello launches only after a final capacity review.

## 1. Plan verification (existing evidence, no rerun needed)

- Configuration decision benchmark: 16 timing-blind guides × 2 assemblies +
  1,000-guide library-mix runs; **11.7 M rows audited, zero errors**;
  m4 = 8.3 s/guide traversal, 5,855 rows/guide, 100 % CFD + CRISPR-Bulge
  (NGG-only targets; NRR changes nothing here — Brunello uses strict GG).
- Both-score NRR-enabled scorer bundle is the same officially validated one
  staged for SBDS (`nrr-tooling-v1/scorer-v2`).
- The cohort worker, checkpoint/restart semantics, per-assembly pins and
  scientific-signature skip are all proven at SBDS scale; the Brunello
  full-run needs a *guide-batching* variant of the same worker (single
  invocation per assembly, 77,441 guides, chunked like the benchmark-1000
  workflow) — the Brunello-specific `full_run_worker.py` draft was withdrawn
  for defects; the corrected chunk/restart logic now exists in the SBDS
  worker and must be ported and failure-tested before launch.

## 2. Resource check (verified today)

- MooseFS free: **41.7 TB** (SBDS campaign will consume: ~3.7 TB rebuilt
  indexes + ~1.4 TB bindings + ~5 GB SBDS outputs; quarantined originals
  ~3.7 TB are recovered after verified deletion). Brunello compressed
  outputs ≈ **19.4 TB** (PAF+annotations+audit records, 24/37.5/30.6 B/row ×
  453 M rows/assembly). Projected total free-space use ≈ 28 TB of 41.7 TB —
  **feasible but tight; requires the deletion program to actually run** and
  a storage reservation decision before cohort launch.
- CPU: SBDS wave (~141 × 30 min) + rebuild campaign (~322 × 0.34 h) finish
  in days at current concurrency; Brunello at ~128,000 CPU-hours
  (464 × 276 h) is the dominant cost: at %16 sustained ≈ 333 days wall —
  **%16 is insufficient**; the campaign needs the full idle cluster
  (hundreds of concurrent 1-CPU jobs when free) and/or explicit scope
  decisions (see §4).

## 3. Full-library single-assembly timing experiment (prepared, NOT launched)

Direct measurement beats extrapolation before any cohort commitment:
one verified assembly (HG00735_pat), m4/b1/z1/GG, all 77,441 guides, one
invocation, chunked with per-chunk SEARCH_COMPLETE checkpoints (restartable).
Measured model predicts ~178 search hours + ~12 h audit + ~86 h scoring
≈ 7.4 days one-CPU. A cheaper linearity validation first: a **10,000-guide
run (~23 h search)** to confirm the 8.3 s/guide marginal at 10× the prior
benchmark scale, then the full run if linear. Both plans are ready to freeze
on authorization. This experiment needs its own explicit authorization.

## 4. The ~128,000 CPU-hour / 19.4 TB estimate and safe reductions

Without changing search correctness, formulas, representative policy or
results:
- **Throughput, not per-row cost, is the lever**: one CPU per assembly ×
  hundreds of concurrent jobs when the cluster is otherwise idle.
- **Separate neural scoring from search** (already the design): search-only
  is ~178 h/assembly ≈ 90,000 CPU-hours cohort-wide; CRISPR-Bulge annotation
  (~86 h/assembly) can run as a second pass and is embarrassingly parallel
  with no search-correctness implications.
- **Audit records (~14 TB of the 19.4)** are validation artifacts; storing
  the full both-score PAF + annotations (≈5.4 TB cohort-wide) as the
  authoritative archive with audit records kept per-assembly only until
  review, or gzip -9/dictionary compression, are storage decisions — never
  filtering (IMPG denominators need every eligible hit).
- No pruning, no CFD thresholding, no representative change is proposed.

## 5. IMPG / Figure 1C preservation requirements (binding on any plan)

Complete eligible both-score PAF + annotations per assembly (actual PAMs,
coordinates, strands, CIGARs, both scores, explicit N/A reasons); thresholds
are derived views only and can never define allele denominators; flanking
context re-extracted later from frozen references; deletions/unmapped are
assessability states, never zero-hits. Figure 1C is generated only after the
cohort completes.

## 6. Launch gate (final review required)

1. SBDS campaign complete (indexes + searches) or explicitly deprioritized.
2. Corrected Brunello guide-batching worker ported and failure-tested.
3. Timing experiment (§3) run and analyzed.
4. Storage reservation signed off (§2) including the deletion program.
5. Explicit user launch authorization with concurrency budget.

Until then: no Brunello cohort job is submitted.
