# SBDS: completed guarded pilot and historical comparison

**Update:** standard SpCas9 NGG is confirmed; the new requested NGG/NRR pilot
is complete. See [NGG/NRR results and expansion blockers](sbds_nrr/REPORT.md).
This older report remains the historical NGG/edit-domain comparison.

**2026-11-11 — one assembly, not the requested completed 466-input rerun.**
The other 465 primary searches have not been launched. Full-cohort findings and
six updated population figures remain pending input-integrity/binding gates.

## Inputs and execution

Guide `ACTGAAATCTGTAAGCAGGC`; already bound HG00735_pat, 3,029,527,752 bases,
82 contigs. Native PAM-aware Columba backend, no WFA2 verification, no CFD/neural
pruning, no CFD threshold in search. Same attested development CLI/native source
as the previous clean workflow, not a production rebuild. Both modes retain all
eligible hits. Strict NGG and f0.75 in both configurations.

| Configuration | Slurm | Eligible hits | CFD defined | Neural defined / unavailable | Search wall s | Neural annotation wall s | Native traversal/replay s |
|---|---:|---:|---:|---:|---:|---:|---:|
| Primary m5/b1/z2 | 2871689 | 38,060 | 38,060 | 8,104 / 29,956 | 419.615 | 17.028 | 101.030 |
| Historical-domain control m2/b1/z1 | 2871691 | 61 | 61 | 61 / 0 | 317.455 | 10.017 | 0.196491 |

Both Slurm tasks completed with exit 0. Sampled summed-tree peaks were
14,295,252 / 14,260,628 KiB, below the 48-GiB sampled stop; no guard fired.
RSS is not a hard cgroup limit and can double-count shared pages/miss spikes.
One fresh process/configuration, shared worker and uncontrolled caches: these
are workload observations, not replicated benchmarks or universal speed ratios.

All 38,121 emitted rows passed independent samtools/CIGAR/PAM/CFD checks,
identity/availability checks and exact-anchor uniqueness. All 61 control anchors
occur in the primary output, with unchanged spans/CIGARs/CFDs. This is conditional
recovery of finite emitted-site controls, **not genome-wide recall**.

Primary neural unavailability, with every row retained:

- 13,882 DNA2 alignments exceed the official 24-column input;
- 14,376 one-base-gap/five-mismatch alignments exceed conservative gapped domain;
- 1,698 multi-base RNA bulges are outside the supported policy.

No cropping, realignment, retraining, unsupported predictions or numeric-zero
substitution. Numerical PAF `cb` is absent when unavailable; TSV has `N/A` and
reason. Thirty stratified/tail predictions across both configurations matched
unmodified official inference, max absolute difference **4.9651336e-10**
(tolerance 1e-6), scoring-only Slurm **2871693**. Trusted upstream reference
pickles were used only by the independent reference checker, never the deployed
neural worker. This establishes sampled implementation agreement, not cleavage
accuracy or calibrated probability.

Analytical full-precision CFD views, not search filtering:

| Configuration | None | >=0.05 | >=0.20 |
|---|---:|---:|---:|
| Primary | 38,060 | 4,093 | 595 |
| Control | 61 | 18 | 8 |

Every primary row remains in the authoritative outputs regardless of these views.

## Historical comparison: settings, eligibility, geometry, scoring

Original HG00735 output has **813** rows. Its exact original guide sequence was
verified; the explicit ID mapping is `SBDS_therapeutic -> SBDS`, not a locus
normalization. Current full reference hashes and independent samtools fetches
reconfirmed every original target/PAM against the prior PAM-audit cache.

- **756/813** original observed PAMs are non-NGG; **57** old rows have valid
  strict m2/b1/z1 CIGAR witnesses. All **57/57** occur in both corrected outputs.
- **228** original rows have invalid standard PAF query geometry; **two** CIGARs
  fail actual-sequence replay. These diagnostics are separate from eligible
  CIGAR witnesses; old fields are neither rewritten nor silently repaired.
- Control has four additional anchors beyond those 57 witnesses. This does not
  establish exhaustive new genome-wide recall or diagnose each old engine stage.
- At the 57 shared anchors, **20** current representatives differ from old
  spans/CIGARs; **26** stored CFD values differ. **24** score differences occur
  even with unchanged span/CIGAR. Independently reconstructing the original
  eligible CIGARs with current CFD tables/formula disagrees with **25** old
  numeric scores beyond rounding tolerance. Both strands and RNA/DNA bulges
  occur. These are historical scoring/reconstruction discrepancies, not all
  explainable by representative choice. No current CFD formula was changed;
  exact historical internal causality is not newly source-to-binary attested.
- Example unchanged alignment: `1I2=1X9=1X6=`, observed AGG, old CFD **0.1089**,
  current independently reconstructed **0.4200**. Do not silently treat old
  scores as interchangeable with the current validated implementation.

Selected old patterns, independently re-extracted in this assembly:

| Previous-study locus context | Target sequence | Observed PAM | Old CFD | Corrected exact-anchor status |
|---|---|---|---:|---|
| chr7 target pattern | ACTGAAATCTGTAAGCAGGC | GGG | 1.0000 | retained; one perfect hit, CFD 1.0000 |
| chr12 recurrent pattern | ACTGAATCTGCAAGCAGTC | GGC | 0.2721 | ineligible PAM |
| chr17 recurrent pattern | ACAAAATCTATAAGCAGGC | AGT | 0.2151 | ineligible PAM |

Chromosome context comes from previous GRCh38/CHM13 sequence work; this is not a
new cross-assembly homology assignment. Dropping these two invalid PAM anchors
is not loss of valid NGG-site recall, nor evidence that the entire homologous
region contains no other eligible anchor.

## Full SBDS expansion and figures

The initial audit currently rereads 464/466 historical PAFs; HG03784_hap1 and
NA18983_hap1 return EIO. Earlier complete 461/4/1 perfect-match counts must be
labelled historical, not claimed freshly reproduced. Only 12 index sets have
prior full-content/recovery evidence; one already has the ready native binding.
Another 454 sets still need full reads and sequence binding, not metadata-only
approval. Do not reuse the nine known damaged original index prefixes.

Next preparation should prioritize the four no-perfect cases
HG02258_mat/HG03225_hap1/HG03195_hap1/HG03834_mat, HG00321_hap2, and the exact
GRCh38/CHM13 derivatives. Use verified replacements where available; stop and
report EIO without overwriting originals. Confirm duplicated-context/copy status
before treating HG00321's two old records as a biological duplication. Budget
multiple-terabyte full index reads before authorizing all remaining preparation
and search tasks; never replace an unreadable/missing result with zero hits.

Only after a 466-input completion/assessability ledger passes should the six
historical figure topics be updated: CFD distribution, candidate reduction,
per-assembly candidate distribution, recurrent high-CFD patterns, sequence-level
target variation and perfect-match distribution, with neural coverage and
old/new context supplements. The current PDF/PNG is explicitly an **n=1 pilot
diagnostic**, not any of those completed population figures.

## Artifacts

Root: `/moosefs/salehi/HPRCv2_CRISPRapido/runs/sbds_pam_aware_v1/`.

- `00_manifest/pilot-plan.json`: code/binary/input/model pins and exact limits.
- `01_search/HG00735_pat/{primary,historical-domain-control}/SUCCESS.json`.
- `02_scoring/HG00735_pat/`: full both-score PAF/TSV bundles.
- `03_qc/`: every emitted-row audit; `historical-HG00735-v2/{summary.json,records.tsv}`;
  `official-sample-v1/agreement.json`. Failed comparator v1 is preserved: it
  stopped on an invalid old CIGAR; v2 labels failures without repairing them.
- `04_figures/pilot-v2/pilot-diagnostics.{pdf,png}` and provenance.
- `05_reports/pilot-complete-v2.json` and this report.

Small committed evidence is under `research/hprcv2/results/`; large/full PAFs,
sequence-bearing annotations, models and environments remain external.
