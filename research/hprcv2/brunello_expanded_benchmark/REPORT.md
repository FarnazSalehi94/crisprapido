# Brunello expanded benchmark — configuration decision and full-run plan

**2026-11-12. SBDS is PAUSED** pending confirmation of the prime editor's PAM;
nothing here changes SBDS parameters or assumes NGG applies to it. This
document is Brunello-only. No full 464-assembly search and no full-library
77,441-guide run was launched; everything below is measured on bounded,
guarded, independently audited runs. External evidence root:
`/moosefs/salehi/HPRCv2_CRISPRapido/runs/brunello_pam_aware_v1/`
(`01_pilot/expanded-benchmark-v1/`, `01_pilot/benchmark-1000-v1/`,
`08_reports/expanded-benchmark-analysis.json`).

Subagent note (Task 8): no subagent/agent tool is installed in this Pi
environment, so the configuration/biology, HPC-cost and IMPG-output reviews
were performed serially in this one report; conclusions are cross-checked
against each other here.

## 1. Benchmark design

- **Backends and integrity:** experimental PAM-anchored Columba backend (no
  WFA2), strict `-p GG` (NGG), `-f 0.75`, one application thread, no
  CFD/neural pruning, no threshold in search. Same attested development
  CLI/native binaries and pinned five-model CRISPR-Bulge ensemble as all prior
  milestones. Every emitted row of every run passed the independent
  samtools/CIGAR/PAM/CFD audit; every run ran under the memory/disk supervisor
  with exclusive output directories and checksum-committed completion markers.
- **Match-fraction compatibility check (Task 1):** the backend's rule is
  `matched_bases / guide_length >= f`, where `matched` counts only `=`
  operations; bulge bases neither count as matches nor inflate the denominator.
  With 20-nt guides and `-f 0.75`, eligibility therefore requires
  `20 - mismatches >= 15`, i.e. **mismatches <= 5 — exactly the widest nominal
  setting tested**. All four configurations are fully compatible with `-f 0.75`:
  no configuration's domain is clipped by the match fraction (confirmed on real
  output: candidate-B emits 5-mismatch rows both with and without 2-base
  bulges, e.g. 260,361 such rows on this panel), and the binding constraint is
  the mismatch budget itself. The audit's independent replay enforces the same
  rule, so search and validation agree by construction.
- **Guides (Task 2):** 16 deterministic guides = 15 cells of
  (historical retained-count stratum x GC-content stratum) + one non-targeting
  control; within each cell the median-count guide with an ID tie-break
  (`00_manifest/panel16.{tsv,fa}`). Selection is timing-blind. The historical
  retained-count proxy is thresholded and from legacy-semantics runs — it
  stratifies workload, it does not measure the corrected backend.
- **Assemblies (Task 2):** HG00735_pat (82 contigs, 3.03 Gb) and HG01496_pat
  (105 contigs, 3.05 Gb). HG01496_pat had no binding at session start; it was
  safely prepared via the trusted import path (verified component checksums,
  normalized-text comparison PASS, binding sidecar, one-time 371 s) and a
  pinned scoring reference index. Both bindings are fully verified on every run.
- **Configurations (Task 1):** historical m2/b1/z1, intermediate m3/b1/z1,
  candidate-A m4/b1/z1, candidate-B m5/b1/z2 — one invocation per
  (assembly, configuration) containing all 16 guides.
- **Scaling check (Task 5):** a bounded 1,000-guide benchmark (first 1,000
  frozen-library records, file order) at intermediate and candidate-A on both
  assemblies.

## 2. Results — 16-guide panel, both assemblies

Traversal is the native search-only time (single thread). "s/guide" is
traversal/16; complete-workflow numbers include the fixed ~310 s per-invocation
verification/snapshot overhead, independent audit and both-score annotation.

| Configuration | s/guide (HG00735/HG01496) | Rows per assembly (16 guides) | Neural defined | Workflow wall s | Peak tree RSS |
|---|---:|---:|---:|---:|---:|
| historical m2/b1/z1 | 0.184 / 0.185 | 499 / 496 | 100 % | 371 / 374 | 14.3 GiB |
| intermediate m3/b1/z1 | 1.534 / 1.457 | 7,552 / 7,579 | 100 % | 397 / 396 | 14.3 GiB |
| candidate-A m4/b1/z1 | 9.194 / 8.960 | 88,462 / 88,807 | 100 % | 584 / 583 | 14.3 GiB |
| candidate-B m5/b1/z2 | 93.628 / 99.042 | 560,138 / 561,390 | 17.3 % | 2,045 / 2,141 | 14.4 GiB |

Every run: independent audit PASS on 100 % of rows (zero errors); CRISPR-Bulge
models loaded exactly once per invocation. Candidate-B's unavailable neural
predictions are explicit N/A with reasons (294,273 exceed the official 24-column
input, 148,095 outside the gapped mismatch domain, 20,658 multi-base RNA bulges
on HG00735) — never zero.

### Nested site recovery (exact anchor keys; HG00735, HG01496 in parentheses)

| Step | Added sites | CFD>=0.20 among added | CFD>=0.05 among added | Shared-representative changes |
|---|---:|---:|---:|---:|
| m2 -> m3 | 7,053 (7,083) | 473 (476) | 2,102 (2,142) | 0 (0) |
| m3 -> m4 | 80,910 (81,228) | 1,927 (1,931) | 13,616 (13,717) | 0 (0) |
| m4 -> m5/b1/z2 | 471,676 (472,583) | 1,526 (1,518) | 19,973 (19,975) | 2,667 (2,687) |

Additional-site CFD profiles (HG00735): m3-only median 0.013, max 0.80;
m4-only median 0.005, max 0.73 (41 sites >= 0.50); m5z2-only median 0.0,
max 0.67. **Broader settings recover large numbers of additional high-CFD
candidates up to m4; the m5/b1/z2 step adds mostly low-CFD sites (0.32 % of its
additions >= 0.20 — fewer new high-CFD sites than m4 adds), while changing
2,667 shared-anchor representatives (2,405 stored-CFD changes) because 2-base
bulges create new edit-cheaper representatives.**

## 3. Biological consequences (Task 3)

- Sites appearing only at broader settings are real, independently verified
  NGG-eligible alignments (all replayed against the frozen reference and
  CFD-reconstructed). They are not artifacts — but their CFD distribution is
  strongly skewed low: median CFD drops from 0.013 (m3-only) to 0.005 (m4-only)
  to ~0 (m5z2-only). Low CFD is NOT evidence of biological irrelevance — CFD is
  an uncalibrated prioritization score, and these exact sites are precisely
  what a pangenome allele-frequency analysis needs for denominators.
- High-CFD recovery is concentrated in the m3 and m4 steps: on this panel m3
  adds 473 and m4 a further 1,927 sites with CFD>=0.20 per assembly. m5/b1/z2
  adds fewer new CFD>=0.20 sites than m4 despite ~10x the cost.
- **Representative stability matters for allele frequencies.** m2->m3->m4
  nesting is perfect on this panel (zero representative changes, zero stored-
  CFD changes at shared anchors). m5/b1/z2 changes 2,667 representatives at
  shared anchors; for those loci the *aligned target sequence* — the very thing
  an allele-frequency analysis counts — changes for policy reasons (bulge
  width), not biology. Frequency results are therefore not comparable across
  the bulge-length setting, and the configuration must be fixed before any
  IMPG/frequency analysis.
- CRISPR-Bulge N/A at m5/b1/z2 (83 % of rows) is recorded with explicit
  reasons; those rows still have CFD and remain fully usable for CFD-based
  frequency views. At m3/m4, 100 % of rows have both scores.

## 4. Output-size control (Task 4)

Measured compressed sizes (gzip level 6, real streams) per 16-guide run and
per row (1,000-guide intermediate run: PAF 25.0 B/row, annotations 37.2 B/row,
audit records 30.6 B/row):

| Configuration | full PAF+annotations gzip (16 guides, both assemblies ~) | views >=0.05 / >=0.20 as % of rows |
|---|---:|---:|
| m2 | ~36 KB | 45 % / 21 % |
| m3 | ~0.52 MB | 30.8 % / 7.6 % |
| m4 | ~5.85 MB | 18.0 % / 2.8 % |
| m5/b1/z2 | ~28.8 MB | 6.4 % / 0.7 % |

**Recommended retention strategy: keep the complete eligible both-score PAF
plus the complete annotations TSV per assembly as the authoritative compressed
archive; treat CFD thresholds strictly as derived, regenerable views.**
A CFD threshold can never define an allele-frequency denominator: a hit below
threshold is an allele *present in that haplotype*, and dropping it would
silently inflate "absence". All IMPG inputs (exact anchor coordinates, complete
observed target+PAM sequences, strand, CIGAR, both scores and explicit
unavailable-reasons) are preserved in the archive; flanking context for IMPG
mapping is re-extracted later from the frozen checksummed references, so no
extra per-hit data needs to be stored now. If the m4 archive proves too large
(see §5), the fallback is a validated dictionary/compact format — which needs
a round-trip test and your approval — never silent filtering.

## 5. Full-library performance estimate (Task 5)

From the 1,000-guide library-mix benchmark (intermediate measured; candidate-A
measured; see pending section for candidate-A details) and the 16-guide panel:

| | m2 | m3 | m4 | m5/b1/z2 |
|---|---:|---:|---:|---:|
| Traversal s/guide (library mix) | 0.18 | 1.46 | ~9.1 | ~96 |
| Rows/guide (library mix) | 33 | 523 | (pending 1k-run) | ~23,500 (panel) |
| Search h per assembly (77,441 guides, 1 core) | ~4 | ~31 | ~196 | ~2,070 |
| Audit h per assembly | ~0.1 | ~1.0 | (pending) | ~14 |
| Both-score annotation h per assembly | ~0.2 | ~8.5 | (pending) | ~33 |
| **Total CPU-h per assembly** | ~4.3 | **~41** | **~210–260** | **~2,120** |
| **CPU-h for 464 assemblies** | ~2,000 | **~19,000** | **~100,000–120,000** | **~980,000 (infeasible)** |
| Wall at 100 concurrent 1-CPU jobs | ~1 day | ~8 days | ~42–50 days | ~1.1 years |
| Compressed archive for 464 assemblies | ~90 GB | **~1.8 TB** | **~8–13 TB (pending)** | ~35+ TB |

Notes: audit ~10,800 rows/s and annotation ~0.68–0.74 ms/row are measured
constants from the 1,000-guide runs; neural annotation is a separate stage and
scales with rows, not guides. Per-assembly estimates anchor on the 1,000-guide
library-mix runs (traversal 8.3 s/guide and 5,855 rows/guide at m4; 1.45 s/guide
and 523 rows/guide at m3), validated against the 16-guide panel (9.0 s/guide,
5,529 rows/guide — the panel slightly over-represents expensive guides). These
remain single-measurement extrapolations on two assemblies with warm caches and
a shared node — the earlier performance investigation showed this backend's
linear per-guide extrapolation is structurally sound (its marginal is pure
traversal), while the old pipeline's was not. SBDS-style fixed overhead
(~310 s/invocation) is amortized by guide-batching and is negligible at
full-library scale.

## 6. Recommendation

1. **Recommended scientific configuration: candidate-A (m4/b1/z1)** — the
   widest *feasible* domain: every combination stays inside the CFD matrices'
   designed mismatch domain AND the CRISPR-Bulge supported domain (100 %
   both-score coverage on all 11.7 M audited 1,000-guide rows and 177,269 panel
   rows — zero N/A); it nests perfectly with narrower settings (zero
   representative/CFD changes at every shared anchor, proven at panel AND
   1,000-guide scale); and it recovers ~100 additional CFD>=0.20 sites per
   guide (~7.7 M per assembly) that m3 omits. Measured campaign cost:
   **~128,000 CPU-hours (~53 days at 100-way, ~23 days at 232-way concurrency)
   and ~19.4 TB compressed archive** for 464 assemblies. Choose this if the
   "all eligible off-target hits" deliverable is non-negotiable and the
   cluster/storage budget is approved.
2. **Pragmatic alternative / staged path: intermediate (m3/b1/z1)** —
   **~19,000 CPU-hours (~8 days at 100-way) and ~1.7 TB**, recovering most
   m2-missed high-CFD candidates (~30/guide CFD>=0.20 retained at m3 vs ~130/guide
   at m4). Because m3 -> m4 nesting is *perfect* (zero representative changes),
   an m3-first campaign is fully forward-compatible: a later m4 campaign
   reproduces every m3 anchor byte-for-byte and only adds sites. A staged
   m3-now, m4-upon-budget-approval plan is scientifically safe, but the m3
   omission (~7.7 M CFD>=0.20 sites per assembly) must be stated in any
   allele-frequency denominators produced from m3-only data.
3. **Not recommended: candidate-B (m5/b1/z2)** — ~993,000 CPU-hours (113
   core-years), ~77 TB (exceeds free space), 83 % of rows outside the neural
   model's supported domain, 2,667 representative changes that would destabilize
   allele definitions, and *fewer* additional CFD>=0.20 sites than m4.
4. **Not recommended: historical (m2/b1/z1)** — misses hundreds of CFD>=0.20
   sites even at panel scale.

**This decision (m4 vs m3, concurrency, storage budget) is the user's to make
before any launch; nothing has been launched.**

## 7. Prepared full-run workflow (Task 6) — NOT launched

All scripts are written, unit-checked and committed; execution requires your
explicit authorization of the configuration and scope:

- `prepare_binding.py`: safe in-place symlink binding for any of the 12
  recovery-verified assemblies (full component re-verification against the
  trusted scan; standard binding sidecar; verified decompressed reference;
  scoring reference index). The remaining 452 assemblies need a full-content
  verification campaign first (~5.3 TB of index reads + reference staging;
  array-able, roughly one day of I/O) — metadata presence alone is never
  accepted.
- `make_full_plan.py`: freezes the authorized plan (77,441 guides x 464
  assemblies, chosen configuration, chunk size, all pins, your
  approval/decision reference recorded verbatim; it cannot be fabricated).
- `full_run_worker.py`: one job per assembly, 1 CPU, guides in checksummed
  chunks (default sized so a restart loses at most a few hours), per-chunk
  SEARCH_COMPLETE checkpoints, merged-PAF full independent audit, separate
  both-score annotation stage, explicit N/A reasons, exclusive directories,
  supervisor memory/disk guards, SUCCESS markers with checksums, safe re-runs
  (completed assemblies are skipped after checksum verification).
- `submit_full_run.slurm`: `--array=0-463%<concurrency>` template; the array
  limit is set from your approved concurrency.
- `verify_full_run.py`: completion ledger; missing/changed/failed assemblies
  are reported INCOMPLETE, never as zero-hit results.

## 8. IMPG allele-frequency readiness (Task 7) — design only, nothing run

The proposed archive preserves everything the ALLELE_PLAN requires: exact
anchor keys (contig/strand/PAM-adjacent boundary), complete observed target
and PAM sequences per hit, CIGAR and aligned sequences (representative
ambiguity is reconstructible), both scores with explicit unavailable reasons,
and assembly identity. IMPG flanking-context extraction happens later from the
frozen checksummed references, so the archive stays minimal. Thresholded
views are derived only. Distinct homologous copies, deletion states,
off-reference loci and assessability states remain specified in
[ALLELE_PLAN.md](../ALLELE_PLAN.md); no mapping or frequency number is claimed
here, and the final Figure 1C is not generated.

## 9. Completed 1,000-guide scaling benchmark

Jobs 2871750/2871752 (intermediate) and 2871754/2871755 (candidate-A), both
assemblies, all exit 0; all 11.7 million emitted rows audited with zero errors:

| Measurement (1,000 guides) | m3 | m4 |
|---|---:|---:|
| Native traversal s/guide (HG00735/HG01496) | 1.460 / 1.435 | 8.316 / 8.241 |
| Search wall s | 1,777 / 1,753 | 8,698 / 8,624 |
| Eligible rows | 522,687 / 524,232 | 5,833,431 / 5,876,550 |
| Neural defined | 100 % | 100 % |
| Audit wall s | 47 / 47 | 536 / 558 |
| Both-score annotation wall s | 398 / 404 | 3,963 / 4,021 |
| Complete workflow s | 2,287 / 2,270 | 13,491 / 13,496 |
| Peak tree RSS | 14.3 / 14.4 GiB | 14.95 / 15.0 GiB |
| CFD views >=0.05 / >=0.20 | 128,178 / 30,616 | 896,989 / 129,756 |
| gzip B/row (PAF / annotations) | 25.0 / 37.2 | 24.1 / 37.5 |

**Nested m3 -> m4 at library mix (the decisive biology number):** the
4-mismatch step adds **5,310,744 / 5,352,318 sites per assembly**, of which
**99,140 / 100,494 have CFD >= 0.20** (max 0.93) and 768,811 / 770,749 have
CFD >= 0.05 — while the 522,687/524,232 shared anchors show **zero
representative changes and zero CFD changes**. Scaled to the full library,
candidate-A recovers on the order of **~100 CFD>=0.20 sites per guide
(~7.7 M per assembly) that the intermediate setting omits**. These are genuine,
independently verified 4-mismatch NGG sites in permissive positions — the class
a pangenome off-target/allele analysis cannot silently lose, and the reason m3
must not be adopted without an explicit scope decision.

## 10. Limitations

- 16-panel guide selection is deterministic but not random; the 1,000-guide
  mix (first records in file order) reduces but does not eliminate guide-mix
  bias; per-guide cost varies by off-target load (panel spread at m3: 6 to
  ~1,051 rows).
- Two assemblies, single fresh calls, shared node, warm caches; no repeated
  medians; timing caveats from the performance investigation apply.
- All accuracy statements are emitted-row soundness/containment at these
  settings, not genome-wide recall or biological calibration.
- No algorithm, default, representative policy, scoring formula, production
  binary or historical artifact was modified; SBDS untouched and paused.
