# Brunello A/B safety pilot — review required, no full run launched

**2026-11-11.** Three fixed GC-stratified guides on one already verified
HG00735_pat assembly. This is a completed safety/feasibility pilot, **not a
representative library/cohort cost model or completed 464-assembly analysis**.

## Library and run provenance

Frozen library: 77,441 unique IDs/sequences, resolved 20-nt protospacers;
19,114 targeted gene symbols plus 1,000 non-targeting controls. FASTA SHA256
`a4bac2dc4c8f5ad0c37713751835b416d1669f899de1cd05ca1a0e828be48971`.
Local raw/Unix/clean metadata and exact original FASTA labels agree under the
explicit pre-existing space-to-underscore control-label mapping; no IDs or
sequences were replaced/deduplicated. Published-release equivalence has not
been independently downloaded/certified.

Native PAM-aware Columba, no WFA2 verification, one application thread,
strict NGG, f0.75; both scores and no CFD filtering/pruning. Same development
binaries/reference/helper/model pins in A and B; three guides share one native
index load and one five-component neural ensemble load per configuration.

| Exact original guide ID | GC bases | A hits / neural defined | B hits / neural defined | A CFD>=0.05 / >=0.20 | B CFD>=0.05 / >=0.20 |
|---|---:|---:|---:|---:|---:|
| QKI_AAATAGAGGCAAGCCCAATT | 8 | 5,718 / 5,718 | 30,016 / 6,463 | 666 / 77 | 1,294 / 113 |
| STX12_AGCAGCTGGTCTCATTTGAC | 10 | 7,133 / 7,133 | 35,660 / 7,738 | 1,192 / 166 | 2,555 / 267 |
| TMEM145_CCAGTCGGCCGTAATCCGAG | 13 | 623 / 623 | 4,688 / 723 | 122 / 17 | 376 / 39 |
| **Total** | | **13,474 / 13,474** | **70,364 / 14,924** | **1,980 / 260** | **4,225 / 419** |

A = **m4/b1/z1**, Slurm **2871690**; B = **m5/b1/z2**, Slurm **2871692**.
Both exit 0. CFD was available for every emitted pilot hit. All **83,838** rows
passed independent samtools/CIGAR/PAM/CFD, identity and uniqueness checks.
Unavailable neural predictions remain explicit `N/A` with reasons, never zero:
B has 23,991 DNA2/input-length cases, 24,747 outside the gapped mismatch domain,
and 6,702 unsupported multi-base RNA bulges. Defined coverage: A **100%**,
B **21.2%**. No neural prediction is extrapolated to unsupported inputs.

Sixty-five stratified/tail predictions across A/B matched unmodified official
inference, max absolute difference **8.8974e-11** at tolerance 1e-6, scoring-only
Slurm **2871694**. Official reference-only trusted pickle loading is separate
from the inference worker's HDF5-only path. This validates sampled implementation
agreement, not biological cleavage accuracy, calibration or efficiency.

## Nested search domains are not nested score-prioritized sets

All **13,474 A anchors occur in B**; B adds **56,890**. However, **303** shared
anchors select different valid representatives in the expanded domain, with
**290** different stored CFD values and **303** defined-to-unavailable neural
status changes. Representatives minimize edit/mismatch/CIGAR/start order, not
CFD. No representative policy or formula was changed by this project.

Consequently, A/B downstream CFD views are not nested: **27** A anchors passing
0.05 do not pass B's 0.05 view; **eight** A anchors passing 0.20 do not pass B's
0.20 view. They were **not lost by search**; their B representatives changed.
B adds 2,272 and 167 other passing anchors at those thresholds. Do not call a
higher-threshold result a recall-complete substitute for the full eligible set.

All thresholds use independently reconstructed full precision, not rounded
stored PAF values or a neural cutoff. All primary rows are retained.

## Observed resources and output sizes

One fresh call/configuration, same octopus02 worker, concurrent SBDS neighbour,
uncontrolled caches and no repetitions: not a universal performance ranking.

| Metric | A | B |
|---|---:|---:|
| Search wall, including mandatory binding/snapshot/load/reporting work | 335.498 s | 575.854 s |
| Native traversal/replay, three guides | 22.9669 s | 262.804 s |
| Mean native traversal/replay per guide (not individual timings) | 7.6556 s | 87.6013 s |
| Post-search both-score annotation | 19.029 s | 25.041 s |
| Supervised whole workload | 362.496 s | 614.106 s |
| Sampled summed-tree peak | 14,313,948 KiB | 14,320,468 KiB |
| Full both-score PAF, measured gzip level 6 | 349,899 bytes | 1,437,846 bytes |
| Full annotations TSV, measured gzip level 6 | 547,256 bytes | 2,530,758 bytes |
| All retained diagnostic artifacts (uncompressed; includes duplicates) | 13,022,612 bytes | 69,688,892 bytes |

Analytical PAF-view compressed bytes: A none/0.05/0.20 =
349,899 / 55,780 / 7,953; B = 1,437,846 / 107,577 / 11,953.
These are real compressed stream measurements, not materialized replacements
for the full primary data. No missing CFD arose here; future unscored CFD must
still be preserved separately.

Most A wall time is outside traversal. The residual includes integrity/staging,
reference/index loading, Python/Rust, CFD and formatting; it is **not isolated
index time**. The attested executable has aggregate, not individual per-guide,
native timings. Larger guide batches can amortize work but have not yet been
memory/time validated. Neither current snapshots nor reference hashes were
bypassed to improve timing.

## Cost sensitivity scenarios — NOT estimates with statistical precision

77,441 guides x 464 haplotypes = **35,932,624 guide-haplotype combinations**.
Simply applying the observed three-guide mean everywhere gives:

| Unvalidated linear scenario | A | B |
|---|---:|---:|
| Native-only single-core hours | **76,413** | **874,374** |
| Full both-score PAF gzip, decimal TB | **4.19** | **17.22** |
| Full annotation TSV gzip, decimal TB | **6.55** | **30.31** |
| PAF + full TSV, decimal TB | **10.75** | **47.53** |

These are **sensitivity calculations**, not confidence intervals, bookings or
launch approval. Three non-random, non-repeat-stress guides and one assembly
cannot characterize tails, assembly variation or library-wide costs. CPU figures
exclude batching/setup, scoring, QC, failures, mapping and allele analysis.
Storage figures exclude references/indexes, transient snapshots/spools, legacy
archives, mappings and derived tables. About 53 TB was free on the shared
filesystem at audit time: B's naive PAF+TSV scenario nearly consumes that alone,
without a safe margin. Keeping pilot-style uncompressed duplicates at scale
would be substantially worse and is not an approved retention policy.

## Recommendation and approval gates

**Do not launch full Brunello now.** A is a plausible lower-cost/full-library
candidate only if the narrower m4/b1/z1 biological search scope is acceptable;
B buys many additional anchors at much greater traversal/output cost and mostly
unsupported neural coverage. A is not an equivalent faster implementation of B,
and 100% model coverage is not evidence of better biological recall/accuracy.
Current CLI defaults and the SBDS primary domain remain m5/b1/z2.

Before choosing a full-library configuration:

1. Validate a modest expanded pilot (proposed 16 guides on two verified HPRC
   assemblies, identical A/B inputs): GC strata, non-targeting controls,
   repeat/low-complexity and historical high-retained-count proxies, plus seeded
   held-out selections. Historical CFD>=0.20 counts are biased proxies, not raw
   candidate-count strata. Resolve new reference/index bindings first; budget
   imports, runtime, disk and stop conditions before submission.
2. Test bounded batches (e.g. 16 then 64) under monitored memory/disk limits and
   stage timing. Do not extrapolate the current three-guide residual as a fixed
   per-batch cost, assume 64/256 fit current two-hour limits, or change native
   algorithms/pruning/representatives to make an estimate favourable.
3. Agree a **lossless full-eligible retention policy**. Candidate options are one
   compressed both-score PAF plus full compressed TSV, or one authoritative
   full PAF plus compact sequence/score/status dictionaries and reproducible
   annotation reconstruction from frozen references. The latter needs an exact
   round-trip test and user approval; its compression benefit is not measured
   here. Pilot duplicates have not been deleted. Thresholded PAFs can be derived
   views, never the sole surviving archive.
4. Freeze manifests, matched hardware/versions, node packing, concurrency,
   allocated-core versus application-core budgets, checkpoints and storage
   reservations. The cluster's historic memory scheduling was CR_CPU; sampled
   RSS is not kernel protection. Do not hide extra reserved CPUs in cost claims.
5. Review A/B scientific domain, coverage, cost ranges and retention with the
   user. Only then authorize an explicitly scoped full Brunello launch.

The IMPG/sequence-allele/assessability/off-reference/extended CRISPRme-style
requirements are in [ALLELE_PLAN.md](ALLELE_PLAN.md). No whole-library mapping,
frequency or population specificity result is claimed here.

## Artifacts

Root: `/moosefs/salehi/HPRCv2_CRISPRapido/runs/brunello_pam_aware_v1/`.

- `00_manifest/{assemblies.tsv,pilot-plan.json,brunello-pilot.fa}`.
- `01_pilot/HG00735_pat/{A,B}/SUCCESS.json`: committed checksums, independent
  audits, full raw/both-score PAF and annotations, stage/tree metrics.
- `02_search`, `03_scoring`, `04_impg`, `05_allele_frequency`: no full-run outputs.
- `06_figures/pilot-v2/pilot-diagnostics.{pdf,png}` and provenance.
- `07_qc/official-sample-v1/agreement.json`.
- `08_reports/pilot-complete-v2.json` and this report.

A/B are independent of the SBDS output tree; large inputs/outputs/models stay
outside Git. Small checksummed evidence is committed under `research/hprcv2/results/`.
