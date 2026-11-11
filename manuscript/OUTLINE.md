# Manuscript draft structure — proposal, not a submission-ready manuscript

## Working title

**CRISPRapido: PAM-anchored off-target search and alignment-aware scoring in genome assemblies**

Alternative, if the application becomes adequately supported: **PAM-anchored CRISPR off-target enumeration in haplotype-resolved assemblies with CRISPRapido**.

Do not use “exhaustive genome-wide,” “clinically validated,” “calibrated cleavage prediction” or “complete pangenome characterization” in the title/abstract on current evidence. Authors, affiliations, funding, disclosures, final article type and citations are **TODO**.

### Materials reviewed

Existing manuscript materials in this checkout are figure source tables, plotting scripts and draft legends under `benchmark/figures/{correctness,performance_scaling,candidate_characterization}/`. Their READMEs and source provenance were reviewed along with benchmark and current handoff reports. **No complete prose manuscript was located in the inspected checkout/material paths.** This outline does not replace an unseen author draft. Existing figures were not overwritten; their historical parameters/known PAM issues require a reuse decision in [EVIDENCE.md](EVIDENCE.md).

## Abstract structure

**Motivation (draft):** Assembly-resolved off-target analysis requires accurate handling of PAM-adjacent sites, mismatches and bulges. Candidate discovery, representative alignment and scoring are distinct operations: agreement in one does not establish agreement or completeness in the others.

**Methods (draft):** We introduce an optional PAM-anchored Columba FM-index traversal that reports deterministic eligible representatives without WFA2 verification, retaining the existing CFD implementation and optional pretrained CRISPR-Bulge annotation. We compare exact strand-aware PAM-anchor sets, independently replay reported alignments and scores, and distinguish search-domain exclusions from unavailable model predictions. Original search defaults and alternative pipelines remain available.

**Results skeleton — select only after scope is settled:**

- New 700-base adversarial controls: at m3/b1/z1, indexed pipelines each recover **26/26** expected anchors versus **20/26** for direct-window WFA2. At m3/b1/z2, PAM-aware Columba recovers **30/30**, Columba+WFA2 **28/30**, window WFA2 **24/30**. Candidate/witness probes identify constrained-alignment limitations. These are finite designed controls, **not an unbiased genome-wide recall estimate**.
- Same shared anchor, two independently valid representatives: CFD **0.2700 versus 0.1736**. Existing separate 204-panel audit: **54** different eligible representative pairs, **41** different CFD values.
- Existing HG00735_pat demonstration reports **21,447** single-guide anchors and **78,582** across three guides, all passing emitted-site eligibility/CFD audits; this is **not** genome-wide recall validation or a complete HPRC population study.
- Optional scoring, if kept in the abstract: **4,526/21,447** supported CRISPR-Bulge predictions agree with official inference; **16,921** remain explicitly unavailable. This validates implementation/availability, not biological accuracy.
- **TODO:** final abstract performance claim only from a chosen comparable workload/version/hardware/replicate scope. Tiny controls favor different pipelines than existing 8 kb/64-guide results; no universal speedup should be advertised.

**Conclusion (draft):** Separating PAM-anchor recovery, alignment choice and model availability makes assembly-based off-target analyses more interpretable and reproducible. Broader recall/biological benchmarking, model-domain coverage and portable packaging remain necessary.

**Availability:** TODO stable public revision/archive, installation/container, data/accession records and complete licensing/model-redistribution review. Current research branch and manifests are evidence, not a tagged production release.

## 1. Introduction

1. Why sequence variation/haplotype-resolved assemblies can alter candidate off-target sequences; distinguish assembly-native analysis from projection onto one canonical reference. **TODO verified literature citations**, not new empirical claims about prevalence.
2. Computational requirements: PAM constraints, strand, mismatch/bulge limits, many guides/assemblies and repeat-rich sequence.
3. Why discovery, constraint eligibility, chosen alignment and empirical/model score are different. Explain why one best affine alignment followed by constraints need not find every eligible alignment at an anchor.
4. Related work: WFA2, Columba, Sassy and established off-target tools (e.g. Cas-OFFinder/CRISPRitz/FlashFry/SWOffinder where comparable). **TODO verify primary references, supported domains and novelty positioning; no unmeasured superiority claim.**
5. Contributions: optional PAM-aware indexed traversal; explicit deterministic representative/availability policy; systematic concordance audit; reproducible research workflow; limited assembly demonstration. Pretrained CRISPR-Bulge is an integration, **not a newly trained model or biological validation contribution**.

## 2. Methods

### 2.1 Software taxonomy and provenance

Use [implementation inventory](../research/search_concordance/IMPLEMENTATIONS.md). “Original” is version-dependent: historical window/WFA2, the window branch used in old benchmark figures, and Sassy on divergent main are not synonyms. Report current CLI/native/builder hashes, source commits, dependency/index layout, compiler/profile and all commands. WFA2 is an aligner/verifier, not an independent indexed search engine. No fabricated five-engine taxonomy.

### 2.2 Inputs and search contract

Exact FASTA/guides/PAM/source hashes, source-window offsets, resolved guide requirements, supported parameter ranges and fail-closed unsupported settings. Current defaults: m5/b1/z2/f.75/NGG; comparisons explicitly choose equivalent m3 strata where required. Adjacent I/D columns count as one event of combined size. Reference Ns, ambiguous PAM rejection, contig boundaries and negative orientation must be explicit. Assembly-native coordinates are not automatically comparable across haplotypes.

### 2.3 Pipeline architectures

- Overlapping windows → ends-free WFA2 → eligibility/PAM/reporting.
- Columba `e=m+b*z` plus exact recovery → candidate intervals → anchored WFA2/bounded fallback → reporting.
- PAM-seeded indexed traversal → mismatch/gap tracking → true-reference replay → native canonical representative → existing reporting/CFD. No WFA2 query verifier in this branch.
- Sassy historical diagnostic: unit-edit search, `m+b` budget, explicit indels, distinct PAM/export semantics. Not part of equal-profile primary results until rebuilt/reconciled.

**TODO formal pseudocode, state invariants and independent review for the supported scope.** Finite tests are not a general completeness proof. CFD-guided pruning remains disabled and should not be presented as an integrated optimization.

### 2.4 Reference/index integrity

Vanilla/32-bit/alphabet-5/sparse-SA-4 contract; original N-containing reference binding, component hashes, private snapshot and rechecks. Account for verification/copy costs in workflow timing. Unsigned manifests assume trusted preparation; do not claim adversarial race immunity or graph-native pangenome indexing.

### 2.5 Representative alignments and reporting

Define native comparator `(mismatches+inserted+deleted, mismatches, compressed CIGAR, start)` and the distinct WFA-based selection policy. Neither maximizes CFD. State what one representative per guide/contig/strand/PAM anchor means; all eligible alignments need not be emitted. PAF remains byte-compatible by default; opted-in TSV includes sequence/alignment context and N/A reasons. Reporter/schema defects in old artifacts are measured separately, not silently repaired.

### 2.6 Concordance and independent validation

Exact key: guide sequence + exact reference/contig + strand + PAM-adjacent boundary (`end` plus, `start` minus). No broad overlap equivalence. Report expected/shared/unique/missing/unexpected/duplicates, then start/end/CIGAR/count/score agreement separately. Independent direct-string anchor oracle, CIGAR/reference replay, raw candidate controls, candidate injection and isolated preferred-alignment probes. Label historical schema-adapted witnesses separately from valid pipeline records. Retain every discrepancy, including unresolved internal cause.

### 2.7 Scoring

Existing CFD matrices/formulas/projection conventions, defined versus unavailable scores and four-decimal output precision. `--min-cfd` is explicitly opt-in/disabled by default and is not used in concordance runs. Pretrained five-component GRU-Embd `c_2` classification ensemble, arithmetic mean, pinned upstream/hash-checked HDF5, CPU environment and official-reference agreement. Document the conservative preprocessing-based eligibility policy, not an asserted exhaustive inventory of training data: 24 positions including PAM; DNA2 requires 25; other encodable cases can still be scientifically unvalidated. No clipping/realignment/retraining. Both scores remain uncalibrated and noninterchangeable; N/A is never zero.

### 2.8 Performance and reproducibility

Identical reference/guide/settings per comparison, explicit threads, fresh processes, cache state, full-workflow versus kernel/setup scopes, hardware and repetitions. Distinguish individual-descendant high-water RSS from sampled summed tree RSS and unique physical memory. Different recovery/setup/integrity work prevents naive speed ratios. Record clean source/native/Cargo builds, Python hash lock/offline wheel test, models/rights and known network-download blocker. **TODO portable packaged environment and profile-attested Sassy baseline.**

### 2.9 HPRC application design

Describe only actual HG00735_pat single-guide/three-guide runs and verified index reuse. Counts include exact target controls. **TODO** prespecified additional haplotypes/guides, source accessions, liftover/anchor projection and mapping ambiguity, population denominators, on-target exclusion rules, statistical/biological interpretation and resource plan. No new large experiment is authorized by this outline.

## 3. Results

1. Architecture/provenance and supported scope (Figure 1).
2. Strict PAM-anchor concordance with complete case ledger; current window and verifier limitations; Sassy/historical comparability blocker (Figure 2).
3. Matched-input full-workflow runtime/memory, including tiny-control crossover and honest historical scope distinctions (Figure 3).
4. Mismatch/bulge classes, mixed-event controls and alternative valid alignments with different CFD (Figure 4). **TODO held-out/randomized wider panel and fair external competitors.**
5. Reuse of one validated haplotype and three-guide batch; no population-level variant claim (Figure 5). **TODO broader pangenome application/projection.**
6. Optional scorer coverage, preserved N/A hits, official implementation agreement and resource overhead (Figure 6 or supplement). **TODO biological predictive validation only if accuracy is to be claimed.**
7. Reproducibility/availability and unresolved packaging/model-installation limitations.

## 4. Discussion

- Main methodological lesson: search, eligibility, representative alignment and scoring agreement are separate.
- PAM-aware traversal versus best-alignment-then-filter; limits of candidate-only comparisons and all-valid-site claims.
- Multiple valid representations can change CFD without changing sites; neither software alignment is automatically biological truth.
- Assembly-native analysis versus truly population-wide/pangenome conclusions; coordinate projection and sampling/assembly uncertainty.
- Limited model domain, unavailable versus inactive, no cutoff transfer/score fusion/neural pruning.
- Performance/setup/caching and integrity costs; output-volume and concurrent-writer risks.
- Limitations: bounded/adversarial samples, unresolved historical source/export equivalence, no genome-wide recall oracle, incomplete clean network model install, no broad biological validation.
- Next work requiring review: held-out bounded cases and attested competitors first, then a separately approved large application plan. Do not imply a production release or clinical readiness.

## Required back matter / TODO checklist

Verified citations and comparison baselines; author contributions; funding/conflicts; code/data accessions and licenses (including Columba AGPL and model/dependency terms); reproducible scripts/manifests; supplementary discrepancy tables; explicit failed/unsupported configurations; journal-format/availability policy check. No results/authors/grants should be fabricated to fill these sections.
