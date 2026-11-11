# Proposed figure program

These are **proposals**, not generated new figures or claims that TODO experiments have run. Preserve the existing `benchmark/figures` files. Use source-data tables with exact revision/settings/denominators and explicit unavailable values; do not silently reuse old captions under the new algorithm.

## Figure 1 — Algorithm and architecture

**A. Component graph:** reference/guides → alternative discovery paths (windows, Columba candidates, PAM-seeded Columba; historical Sassy in a separate diagnostic box). WFA2 appears downstream of windows/candidates, **not** after PAM-seeded traversal. Shared eligibility/reference replay/reporting and CFD, with optional neural worker as a post-search branch.

**B. PAM-oriented traversal:** plus/minus PAM boundary, fixed suffix, mismatch/gap-state bookkeeping, true-reference N handling and contig edges. Label supported b/z/index layout and immutable/bound snapshot checks.

**C. Reporting contract:** one declared canonical representative per PAM anchor; preserve unsupported-scoring hits; separate CFD/CB/N/A. Show inference loads once across chunks/guides. No score-guided pruning arrow.

Source: current source/inventory and [integration report](../research/pam_anchored_poc/integration/REPORT.md). Existing figure1 architecture is a historical two-pipeline Columba/WFA2 schematic and requires redesign, not a new engine label pasted onto WFA2.

**TODO:** final pseudocode/invariant review and vector rendering. Do not draw five independent engines.

## Figure 2 — Search-mode concordance

**A. Primary UpSet/membership panel:** three current complete pipelines, exact guide/reference/strand/PAM-boundary keys, separate strata m3/z1 and m3/z2. Current counts: **20/26/26**, shared **20**; **24/28/30**, shared **24**. Truth expected **26/30**; label both missed and unexpected eligible anchors (zero unexpected here).

**B. Diagnostic lane:** historical WFA2 and Sassy raw/export-valid/independent-eligible counts, hatched/grey, **not equal Venn sets**. Sassy's missing mandatory PAF field, fixed gap span, z2 radius and build-attestation limitations block a fair five-way engine figure. The diagnostic 7/7 all-artifact witness intersections must not be advertised as five-way equivalent pipeline concordance.

**C. Discrepancy mechanism:** preferred alignment versus independent eligible witness; raw candidate support/injection results. Mixed-only truth **2**, WFA-based pipelines **0**, native **2**. Link all [70 cases](../research/search_concordance/results/cases.json), not selected omissions only.

Source: [summary](../research/search_concordance/results/summary.json), [shared keys](../research/search_concordance/results/shared_sites.tsv), [all audited records](../research/search_concordance/results/records.tsv).

**TODO:** profile-attested Sassy/historical rebuild and authorized semantics/exporter reconciliation; larger held-out/randomized bounded panel; independent review. No chromosome/genome query needed to draw current finite-control results.

## Figure 3 — Runtime and peak memory

**A. New matched-input 700-base control:** all three current pipelines, three fresh invocations, median plus individual observations, explicit recovered-anchor counts and different integrity/setup work. **0.004152 / 0.016266 / 0.114388 s**, max individual RSS **7,076 / 20,376 / 15,100 KiB**. Native is not fastest on this tiny input. Do not call this a kernel-speed or equal-work ratio.

**B. Existing 8 kb/64-guide comparison:** historical Columba+WFA2 **12.37/11.49 s**, native **0.42/0.43 s**, same 172 expected anchors; label versions and two replicates. It illustrates workload dependence, not a universal performance crossover law.

**C. Optional larger evidence panel:** existing manuscript scaling series only after strict-PAM/raw-data provenance re-audit. Separate hardware, parameters, baseline replicate counts, cache conditions and memory definitions. Never pool individual-child RSS with summed/aggregate RSS on one unqualified axis.

**TODO:** harmonized held-out bounded scaling; model-free engine/source-profile matched Sassy; larger matched workflow profiling only after presenting/approving the plan. Historical whole-genome timing is reused, not rerun. Existing figure2 legends have mixed hardware/replicate scopes that must remain visible.

## Figure 4 — Mismatch, bulges and alignment-dependent scoring

**A. Eligibility strata:** exact, mismatch-only, one-/two-base I and D, adjacent mixed I/D, both strands, reference N. Report independent anchor denominators and all eligible alignment classes where they overlap; a single canonical output class is not all possible classes at a site.

**B. Two explicit mechanistic cartoons:** mixed witness versus rejected WFA optimum; terminal D versus clipped preferred boundary/nonmatching PAM. These are bounded examples, not biological cleavage observations.

**C. Shared-site alternative alignments:** at `D5_2+`, boundary 22: `[1,22) 1X2=1D2X15=` CFD **0.2700** versus `[3,22) 1I1=3X15=` CFD **0.1736**. Both independently eligible; gap direction is not a uniquely inferred biological variant.

**D. Reused broader representative audit:** **54** valid differing pairs, **41** CFD differences, all **108** independent-versus-Rust CFD checks agreed; label that separate 204-site dataset and version, do not pool with the new 26-site stratum.

Source: [new concordance report](../research/search_concordance/REPORT.md), [existing representative review](../research/pam_anchored_poc/alignment_review/REPORT.md). Existing figure4 mismatch-only result **1,038/1,038** is conditional baseline-valid recovery under historical rules, not an independent genome-wide all-site recall denominator.

**TODO:** category-specific held-out oracle recovery and comparable external tools; exact definitions for overlapping class denominators. Do not label an alternative valid CIGAR “wrong” merely because CFD differs.

## Figure 5 — HPRC / pangenome application

**Currently supported panel:** HG00735_pat, **3,029,527,752 bases / 82 contigs**; single-guide **21,447** reported anchors and three-guide **78,582**, all emitted rows eligibility/CFD-audited, exact controls recovered. Show input/index provenance and native versus setup/total time scopes. Explicitly label one haplotype and three guides; counts include targets.

**Proposed application panels (TODO):** matched reference/haplotype target sets after declared coordinate projection; variants altering candidate sequences; guide-to-guide/haplotype heterogeneity with prespecified denominators and mapping uncertainty; independently inspect apparent haplotype-specific sites. Source/accession and natural-variant interpretation must be verified before plotting.

Do not turn assembly-coordinate count differences into variant-specific off-target gains/losses without projection/sequence evidence. No complete 100-guide × 12-assembly campaign or population-level effect estimate is available. No batch-versus-three-single-job speedup was measured for the three-guide run. HPRC application expansion requires a separate plan and approval; **do not run it for this manuscript scaffold**.

## Figure 6 — CFD and CRISPR-Bulge scoring (optional / likely supplement)

**A. Availability waterfall on saved HG00735:** all **21,447** hits/CFD scores retained; **4,526** CB / both scores, **16,921 N/A**. Reasons: **8,176 DNA2 shape**, **7,624 DNA1/five-mismatch policy**, **1,121 RNA2 policy**. Distinguish hard shape failure from encodable-but-unvalidated policy exclusion.

**B. Paired supported-score distributions:** only **4,526** rows, log-aware axes/zero handling without replacing unavailable with zero. Spearman **0.4281**, Pearson **0.3447** from stored four-decimal CFD; class stratification/ties and subset selection must be shown. No calibration, common threshold or biological predictive superiority claim.

**C. Implementation/reproducibility and overhead:** all **4,526** saved-result predictions matched official inference; complete bounded workflow **1,177** panel-row checks also passed within `1e-6`. Distinguish 13–14 s cached-offset saved-file annotation from ~5 s tiny full-workflow inference and one-time reference-index preparation. New 46-wheel offline environment reproduced the 172-row workflow.

Source: [CLI integration](../research/crispr_bulge/CLI_INTEGRATION.md), [assembly evaluation](../research/crispr_bulge/ASSEMBLY_EVALUATION.md), [reproducibility report](../research/reproducibility/REPORT.md).

**Recommendation:** keep this supplementary unless the paper's primary contribution includes a substantive validated scoring analysis. It is pretrained-model integration, not model development. **TODO biological ground truth/cross-guide generalization/representative sensitivity** before any accuracy claim; no retraining or unsupported-score invention is proposed.

## Submission-layout decision

Six proposed figure topics do not require six main figures. A compact software/method article could combine architecture+concordance, combine detection+alignment effects, and move scoring/reproducibility/historical diagnostics to supplements. **TODO choose final article type and check current journal display/length policies before rendering/submission.**
