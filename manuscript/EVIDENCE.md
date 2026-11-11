# Claim-to-evidence ledger and missing work

## Reviewed existing manuscript assets

- [Workflow/correctness figure package](../benchmark/figures/correctness/README.md): source data, drawing code and draft captions; compares historical window and automatic Columba+WFA2, not five engines.
- [Performance figure package](../benchmark/figures/performance_scaling/README.md): mixed guide counts/hardware and unequal replicate counts explicitly documented; do not remove those caveats.
- [Candidate/mismatch figure package](../benchmark/figures/candidate_characterization/README.md): separates old m0 whole-genome composition from m1/m2 mismatch-only chr22 tests. Its 4,249 old final Columba rows include **701 GG and 3,548 non-GG** records. These are historical outputs from before strict-PAM corrections, not the current strict-NGG pipeline result.
- [chr22 100-guide benchmark](../benchmark/performance/chr22_100_guides/README.md): historical locus-equivalence and PAM limitations; observed 53.3× pilot ratio is not a current strict-PAM-anchor performance claim.
- [Project handoff](../docs/PROJECT_HANDOFF.md), [representative review](../research/pam_anchored_poc/alignment_review/REPORT.md), [integration benchmark](../research/pam_anchored_poc/integration/REPORT.md), [reproducibility](../research/reproducibility/REPORT.md), [scoring](../research/crispr_bulge/CLI_INTEGRATION.md).

No complete prose draft was located in the inspected checkout/material paths. The new outline/figure plan is a proposal. Existing figures/data are preserved and are not silently certified under today's stricter contract.

## Usable claims and boundaries

| Claim / current evidence | Status and safe wording | Missing work / prohibited inference |
|---|---|---|
| New strict-anchor bounded comparison: common z1 **20/26 window**, **26/26 Columba+WFA2**, **26/26 native**; z2 **24/30**, **28/30**, **30/30**; mixed **0/2**, **0/2**, **2/2** | Primary finite-control evidence; all current emitted alignments/CFD independently checked, no unexpected eligible/duplicate anchors | Not random/unbiased recall, not genome-wide completeness, not intrinsic Sassy quality |
| Complete [70-case ledger](../research/search_concordance/results/cases.json), [254 raw-row audits](../research/search_concordance/results/records.tsv), candidate union/injection and preferred-alignment probes | Every bounded missing/raw-additional case investigated for eligibility/export validity; current window/verifier limitations localized | Historical/Sassy internal root cause/source attestation remains TODO; do not fake a fair five-way Venn/rank |
| One new shared valid alignment pair gives CFD **0.2700 vs 0.1736**; old separate 204 panel gives **54** differing valid pairs / **41** different CFD values | Site agreement does not imply alignment/score agreement; both representations can be eligible | No biological preference or score calibration claim; do not pool separate datasets |
| Current full input/output/scoring bounded tests recover **566/566**, **570/570**, **80/80**, **172/172** for PAM-aware CLI scoring modes | Implementation/finite oracle recovery with all rows preserved in CFD/neural modes | These are different panels/settings from new m3 comparison; not proof that all other search pipelines recover every anchor |
| Fresh source/native/Rust builds; 164 Rust tests, one old ignored; 20 scorer tests; offline 46-wheel venv reproduces full 172-row output | Research reproducibility evidence with explicit cache/toolchain/source constraints | Not a clean OS/network bootstrap, portable installed binary or security certification; model network fetch failed |
| Historical old figure baseline-valid **127/127** recovery plus **574** additional oracle-valid candidate loci | May be cited only as conditional historical benchmark under its recorded locus/PAM rules | Not all eligible genomic PAM anchors; old non-GG/raw-output and representative-equivalence rules need re-audit before primary reuse |
| Historical mismatch-only **1,038/1,038** baseline-valid loci (955 plus / 83 minus); m1 **31/31** | Conditional recovery in that chosen chr22/10-guide experiment | Not independent all-locus genomic recall or bulge sensitivity at m5/b1/z2 |
| Existing 8 kb/64-guide/172-anchor workflow: Columba+WFA2 **12.37/11.49 s**, native **0.42/0.43 s** | Two observed historical full-workflow replicates on one bounded workload | Not universal speedup; current 700-base controls have opposite ordering; phases/integrity work differ |
| HG00735_pat jobs 2871278/2871283: **21,447** anchors, exact-target recovery and all emitted CFD/eligibility audits; optimized total **645.6085 s** | Completed single-haplotype/single-guide workflow observations; pilot outputs byte-identical | No independent genome-wide recall denominator; no new assembly-wide comparator run authorized |
| Job 2871284: **78,582** anchors for exactly three guides; all rows independently eligible/CFD-validated; total **837.8747 s**, native search **257.493 s** | Completed three-guide one-haplotype batch, not a future experiment | Not a population/pangenome effect study; not three-single-jobs versus batch speedup |
| Saved HG00735 scoring: **4,526 CB**, **16,921 N/A**, all **21,447 CFD/rows** retained; full official parity within 1e-6 | Coverage/implementation agreement of existing pretrained model on supported representatives | Not biological accuracy/cleavage probabilities or support for all default search alignments |
| Paired score rho **0.4281**, Pearson **0.3447** among supported rows | Descriptive association, conditional on model eligibility and rounded CFD | No direct calibration, transferred cutoff, fusion benefit or ranking superiority |

## Mandatory TODOs before strong manuscript claims

1. **Name/version resolution:** choose exactly what the paper calls original; preserve old figure baseline versus historical standalone versus divergent Sassy distinctions.
2. **Fair Sassy comparison:** pinned dependency/release build, source-to-artifact attestation, standard/semantically faithful exporter and equivalent search/PAM/error budgets. Any fix is a separate reviewed change, not hidden in benchmark normalization.
3. **Held-out validation:** randomized/independent bounded references and guides, Ns/edges/repeats/parameter range; report every missing/extra eligible anchor. Predefine selection and denominators; avoid learning success cases then calling them unbiased accuracy.
4. **External tools / novelty:** verify supported domains and primary literature, then fair bounded competitor configurations. Do not fabricate results for Cas-OFFinder/CRISPRitz/FlashFry/SWOffinder or claim novelty solely from wrapper integration.
5. **Historical figure reuse:** re-audit existing raw outputs under strict observed PAM/anchor/strand definitions and corrected CFD/representative accounting; do not rerun completed whole genomes without approval.
6. **Pangenome application:** explicit cohort/guide sampling, source/assembly accessions, coordinate projection/mapping uncertainty, variant attribution and on-target rules. No full 100-guide × 12-assembly result exists in the validated evidence used here.
7. **Scoring science:** predictive ground truth, representative sensitivity and broader guide/domain generalization if claiming more than implementation parity. Model-domain extension/retraining is not part of this work.
8. **Research distribution:** reliable model provisioning (network helper failed), portable coherent build/relocation/container, rights/security review (Columba AGPL, project license arrangement, model/dependency terms), resource/concurrency limits and public archiving.
9. **Editorial checks:** current article type/length/figure/data/code policies, verified citations, authorship/funding/conflicts and archive DOI/accessions. All remain TODO; no acceptance prediction or production-release implication.

## Larger-experiment gate — plan only

First complete the bounded held-out/source-attested comparison and raw-data re-audit. Then present: fixed assembly/guide manifests and hashes; identical requested PAM/mismatch/gap/fraction rules; supported-budget matrix per tool; expected output/storage/cost limits; matched hardware/threads/replicates; full-workflow versus phase metrics; safe stop conditions; independent recall/eligibility plan and every-discordance ledger. Only after user review may a new chromosome/assembly/cohort run be considered. Do not automatically repeat completed whole-genome experiments, merge main, change defaults or enable pruning.
