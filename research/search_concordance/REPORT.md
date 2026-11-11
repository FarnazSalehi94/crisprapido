# Bounded search-mode concordance and alignment/score audit

## Decision and scope

Baseline **`ce52c15`**, branch `feature/pam-anchored-search-poc`. **No search algorithm, default, scorer, production binary or main branch changed. No large experiment was launched.** The new work is a benchmark/auditor, explicit case/record ledgers, implementation inventory and manuscript planning.

**These names do not denote five independent engines.** See [IMPLEMENTATIONS.md](IMPLEMENTATIONS.md) for exact source lineages, executable hashes, commands and limitations. Three current configurations share one attested release executable: direct-window WFA2, automatic Columba+WFA2, and PAM-anchored Columba without WFA2. Historical standalone WFA2 and divergent Sassy artifacts are separately labelled diagnostics. A fair release/profile/source-attested five-configuration comparison is **not complete**; it must not be advertised as five independent search-engine results.

## 1. Inputs and matching rules

New tests are **700 bases / 28 synthetic contigs** (existing resolved fixture subsets plus an explicit two-PAM dense control), one exact guide `GAAACAGTCGATTTTATCAC`, both strands, and a separate **46-base / two-contig mixed I/D control**. These are adversarial bounded controls, not a random biological sensitivity sample. References/guides are identical within each comparison. [Input provenance](results/implementations.json) records reference/guide/index-component hashes; post-run reference/components match their trusted pre-run build binding, and guide bytes and independent oracle sources are checked.

- **Common stratum:** m3/b1/z1/f.75, strict complete observed NGG, one thread. Sassy's `m+b` radius equals `m+b*z=4` here, with explicit `--allow-indels`; historical PAM restriction is external validation, not a repaired producer.
- **Z2 diagnostic stratum:** same reference at m3/b1/z2/f.75. Modern modes have equivalent requested limits; Sassy's radius 4 is insufficient for maximum budget 5, so Sassy cannot support an equivalent-z2 capability claim.
- **Mixed-only stratum:** m0/b1/z2/f.75, two eligible adjacent I/D events (one contiguous event with combined length 2). Sassy radius 1 versus required budget 2 is again non-equivalent.

A site key is **(exact guide sequence, reference identity/contig, strand, zero-based PAM-adjacent target boundary)**: target `end` on plus, target `start` on minus. Minus PAM bases are the reverse complement of `[start-3,start)`; the key is **not** `start-3`. There is no broad-overlap or shifted-boundary equivalence. Window coordinates must be projected with their recorded source offsets before comparisons across separately named genomic windows; no such cross-reference aliasing is performed here.

Each original CIGAR is independently replayed against the supplied reference/guide. Checks cover coordinate bounds, full consumption, actual `=`/`X`, mismatch and contiguous gap-event counts/combined size, at least 15/20 matches, strand and complete NGG. A direct-string forward oracle enumerates all allowed gap words and representatives; separate anchor-local enumeration supplies witnesses for every discordant/raw-additional anchor. Both standard-export validity and mathematical alignment eligibility are recorded. **No endpoint, query coordinate or CIGAR is silently repaired.**

Independent CFD is computed only for eligible, auditable geometry; reported CFD agreement is evaluated separately. `N/A` / `not_auditable`, absent historical scoring, unavailable mixed-gap CFD and actual numeric zero are distinct. Default CFD-only scoring is used; no threshold or pruning is enabled.

## 2. Site concordance

### Current complete pipelines

| Panel | Independent expected | Window + WFA2 | Columba + WFA2 | PAM-aware Columba | Shared by all three |
|---|---:|---:|---:|---:|---:|
| Common m3/z1 | 26 | **20** | **26** | **26** | **20** |
| m3/z2 | 30 | **24** | **28** | **30** | **24** |
| Mixed-only m0/z2 | 2 | **0** | **0** | **2** | **0** |

**Zero unexpected eligible sites and zero duplicate standard PAM anchors** occurred in these three configurations. Every emitted current record was independently eligible and its defined CFD matched the independent calculation at four-decimal output precision (`0.000051` tolerance). PAM-aware representatives also agree with the independent oracle's declared selection policy. Equal site counts are not alignment/score identity.

The exact **20 / 24 shared-site lists** are in [shared_sites.tsv](results/shared_sites.tsv). Each stratum is separate; shared fixture sites are not independent biological observations.

### Historical / Sassy diagnostics — not another equal Venn set

| Panel | Historical raw / eligible witnesses / standards-valid sites | Sassy raw / eligible witnesses / standards-valid sites |
|---|---|---|
| m3/z1 | 21 / **11** / 10 | 24 / **16** / 0 |
| m3/z2 | 25 / **13** / 11 | 28 / **16** / 0 |
| Mixed-only | 0 / 0 / 0 | 0 / 0 / 0 |

Sassy's known eleven-field schema can be read diagnostically but is not valid 12-column PAF; gapped target spans may disagree with CIGAR consumption. Historical WFA2 exports expose query-coordinate/negative-CIGAR/base-consumption problems. Both older producers lack the same observed-PAM output contract. Invalid/ineligible rows remain in the audit, not discarded without explanation.

The diagnostic intersection across all five artifact configurations is **7 / 7 / 0 eligible witnesses**, but **must not be represented as equivalent five-way complete-pipeline concordance**. Sassy has zero standards-valid PAF sites under this exporter, and its z2 budget, ambiguity policy and build profile are not equivalent. Source-attested rebuilding and a separately reviewed exporter/semantics reconciliation are TODO; no fixes were made under this task.

## 3. Every discordance and raw additional record

[records.tsv](results/records.tsv) records the audit outcome for **all 254 first-replicate raw records** across panels/configurations; repeated common-stratum records were audited and checked for deterministic record content. [cases.json](results/cases.json) contains **all 70 anchor cases**, including all independent expected anchors and raw additional anchors, their witnesses, mode membership, linked row IDs and controlled probes. These are complete bounded ledgers, not a selected example list. Raw extra boundaries that fail requested PAM/geometry are not called additional eligible off-target sites.

### Current window misses: controlled preferred-alignment evidence

The following exhausts current window misses:

- **Terminal DNA gap, both strands:** `D20_1` has eligible `20=1D`, CFD **0.6**, at the true NGG boundary. Diagnostic `--no-filter` exposes preferred `20=` ending one base earlier (or one base later on the forward coordinate for minus), whose observed PAM is **not NGG**. Source `scan_window` uses one ends-free optimum and trims reference overhangs; downstream PAM filtering does not enumerate the valid terminal-D alternative.
- **Constrained alternative at z1, both strands:** `D5_2` has valid one-base-gap/three-mismatch alternatives. The window's preferred CIGAR is `5=2D15=`, violating z1; the valid higher-cost alternative is not enumerated by the window scan. At z2, this preferred gap becomes allowed and both anchors are recovered.
- **Two PAM anchors in one full window, both strands:** the dense controls contain two exact sites. Only one optimum per window/strand is emitted. Controlled preferred records show the competing valid exact anchor: plus 47 instead of missed 20, minus 3 instead of missed 30. Native/Columba candidate enumeration recovers both. No broad overlap merge is used to hide the omission.
- **m3/z2 mixed alternative, both strands:** `suffix_GGGG` has the eligible mixed witness `16=1D1I3X`. Window diagnostics prefer `16=4I`, violating the allowed gap size. The eligible constrained alternative is missed.
- **m0/z2 mixed control, both strands:** valid `5=1D1I14=` has zero mismatches. Window diagnostics prefer `5=1X14=`, which violates m0.

Diagnostic `--no-filter` is **only an isolated causal probe**, never part of primary concordance/timing. Its extra unfiltered records are stored separately and not called eligible predictions.

### Columba+WFA2 omissions: candidate controls and injection

Separately running the same `e=m+b*z` and exact `e=0` Columba searches, merging their raw SAMs, and importing that union produces **byte-identical** output to automatic Columba+WFA2 in every stratum. Thus comparisons did not accidentally substitute a different candidate configuration.

For every missing Columba+WFA2 anchor, an independently eligible oracle interval/CIGAR was explicitly supplied as a valid tagged SAM candidate (`NM` and `AS`, reference-oriented negative CIGAR). **None of the four missing panel-anchor cases was recovered by candidate injection.**

- m3/z2 `suffix_GGGG`, two strands: raw candidates exist on the contig/strand but not at the exact witness interval/anchor. Injecting that exact interval still fails. Unfiltered anchored verification at the requested boundary selects `16=4X` (four mismatches), while `16=1D1I3X` is an eligible three-mismatch alternative. Loss cannot be attributed solely to absent raw discovery.
- m0/z2 mixed control, two strands: the raw candidate interval **and PAM anchor are present**. Injection still fails; unfiltered anchored output selects `5=1X14=` rather than the valid zero-mismatch mixed path.

This agrees with the inspected fallback's equal-length limitation: it does not exhaustively enumerate adjacent mixed I/D alternatives after the preferred WFA alignment fails constraints. Raw support, injected recovery booleans and actual preferred anchored records are retained per case. **WFA cost-optimal alignment plus filtering is not the same operation as testing existence of every constraint-eligible alignment.** No algorithm was changed to resolve this limitation.

### Historical omissions/additions

All are investigated for independent witness eligibility, original field/schema/coordinate/base assertions and common requested constraints in the linked ledgers. Observed m3/z1 historical rejections include six base-assertion failures, two guide-consumption failures and two non-NGG records; counts of reasons can overlap in Sassy records (span mismatch, reference bounds, non-NGG). Query-coordinate defects and the missing Sassy block-length field are recorded separately from alignment eligibility.

**Underlying historical/Sassy engine-level causal attribution is not complete**, because reporting defects can hide a valid internal hit and exact old build attestation is unavailable. Do not mislabel every absent valid export as a search-engine miss, or claim these diagnostics establish intrinsic Sassy recall. A faithful release/profile-matched rebuild and separately authorized reporting/semantics review are required before that claim. Every observed case is recorded; unresolved internal cause is explicitly TODO, not invented.

## 4. CIGAR, mismatch/bulge and CFD comparison

The main comparison separates site membership from representative fields: start/end, CIGAR, mismatch count, I/D base counts, event count, maximum event size and independent CFD. All modern shared records are checked independently, not presumed equivalent because anchors agree.

At m3/z1, the two indexed pipelines share **26/26 sites**, but **one representative and CFD differ**:

| Shared anchor | Columba+WFA2 | PAM-aware Columba |
|---|---|---|
| `D5_2+`, plus boundary 22 | `[1,22)`, `1X2=1D2X15=` | `[3,22)`, `1I1=3X15=` |
| Mismatches / gap | 3 / one-base DNA gap | 3 / one-base RNA gap |
| Independent CFD | **0.2699724514822287** | **0.1735537189834711** |
| Reported CFD | **0.2700** | **0.1736** |

Both consume the full guide, satisfy common limits and have the same complete observed PAM. Neither is declared the biological truth. The native representative minimizes `(edit bases, mismatches, compressed CIGAR, start)`; WFA-based reporting uses its obtained affine alignment/selection policy. Neither selects to maximize CFD. Construction label `D5_2` does not force a unique molecular/represented gap interpretation.

Window and Columba+WFA2 agree in all compared representative fields on their 20 common z1 sites. At z2, Columba+WFA2 and native agree on all fields/CFD for their 28 shared sites; the two extra native mixed representatives have **unavailable CFD**, not zero. Historical original has no CFD field/implementation; all 16 eligible Sassy diagnostic records per resolved stratum agree with independent CFD, but invalid-span outputs cannot be scored reliably from their exported geometry. The nominal Sassy source's unavailable-to-zero conversion is a documented reporting hazard, not inferred inactivity.

### Reused, not rerun

The existing [204-panel representative audit](../pam_anchored_poc/alignment_review/REPORT.md) independently validated all **54 differing pairs**, with **41/54 different CFD values** and all **108** Rust-versus-independent scores agreeing. It documents five pairs with mismatch/event/size/span changes and all 54 with gap-placement changes. The [integrated 204/172-panel audit](../pam_anchored_poc/integration/REPORT.md) recovers all expected anchors for both indexed pipelines, with 204 defined CFD values in the resolved panel and 97 defined / 75 unavailable in the 172 panel. These are distinct datasets/versions; counts are not pooled with the new 700-base controls.

## 5. Runtime and memory

Common 700-base input, m3/b1/z1/f.75/NGG, one thread, **three fresh invocations per artifact**, rotated order. Includes CLI startup, search, applicable index loading/copies/checks, WFA where applicable, CFD/reporting. Excludes one-time tiny-index construction, independent validation and causal probes. Warm/uncontrolled storage caches, shared host, tiny-call setup dominates.

| Current complete configuration | Median wall s | Max individual descendant RSS KiB | Max sampled summed tree RSS KiB | Eligible anchors |
|---|---:|---:|---:|---:|
| Window + WFA2 | **0.004152** | **7,076** | 27,112 | 20/26 |
| Columba + WFA2 | **0.016266** | **20,376** | 27,060 | 26/26 |
| PAM-aware Columba | **0.114388** | **15,100** | 46,740 | 26/26 |

The native route is **not fastest on this tiny panel**. It performs mandatory snapshot/native binding checks that the candidate-generation path does not perform in the same way. Different work/output and security setup prevent interpreting these totals as isolated algorithm speed ratios. The JSON `ranking_applicable` flag only distinguishes current release-profile artifacts from historical diagnostics; it does not certify equal work or completeness.

Historical raw-output timings are retained separately: old WFA2 **0.004249 s / 7,024 KiB**, debug-path Sassy **0.016716 s / 7,108 KiB**. These include different PAM/export/scoring contracts and unverified old build profiles. **No fair Sassy/original speed ranking is claimed.**

Sampler requested 20 ms, 2 GiB sampled tree stop and 120 s child timeout; no stop fired. Sub-20-ms calls can finish between samples, summed RSS can double-count shared pages/miss peaks, and child `ru_maxrss` is an individual-descendant high-water mark, not unique physical memory. Three tiny repetitions provide no large-input scalability guarantee.

Reuse the existing [8 kb / 64-guide / 172-hit complete-workflow benchmark](../pam_anchored_poc/integration/REPORT.md#bounded-complete-workflow-performance): Columba+WFA2 **12.37/11.49 s**, native **0.42/0.43 s**, same expected eligible anchors, different valid representative policies. This historical opposite workload trend and the old manuscript scaling data must be labelled by dataset/version/hardware/integrity scope, not combined into a universal speedup. No chromosome/genome rerun occurred.

## 6. Reproduction, tests and artifacts

Use [benchmark.py](benchmark.py) with explicit current CLI/native/Columba/builder and historical diagnostic paths, plus a **NEW external output directory**. The driver accepts no arbitrary genome input. [export_results.py](export_results.py) publishes small complete record/case ledgers and hashes the full external artifacts. Exact commands are in [summary.json](results/summary.json); source/artifact distinctions are in [implementations.json](results/implementations.json).

**Six new audit tests passed**: strand-aware anchor normalization, labelled legacy schema, no endpoint/query repair, mixed-gap witness/unavailable CFD, required injected SAM tags/negative CIGAR orientation, and duplicate preservation. Also reran **4 existing reproduction-auditor, 9 CFD-model and 20 configured scoring Python tests**, all passed. Core Rust source is unchanged; the prior ce52c15 164-pass record is reused, not presented as a fresh Rust rebuild. Existing core remains unchanged; no new production build or genome query is needed for these tests. Source/executable pre/post hashes and independent witness checks passed. Initial development probes omitted required SAM tags and stopped; those harness attempts are preserved, not pipeline results.

External root: **`/export/local/home/salehi/crisprapido-concordance-v1/`**, authoritative run **`final-v4/`**. It retains original PAFs, raw Columba SAMs, independent oracle injections, unfiltered causal probes, full audited rows, metrics, logs and full summary. Earlier development runs are retained; final-v4 was chosen after completing the diagnostic protocol, not for favorable timing. Published [records](results/records.tsv), [cases](results/cases.json), [shared sites](results/shared_sites.tsv), [summary](results/summary.json), [inventory](results/implementations.json) and [tests](results/tests.log) contain no model/index/binary/genomic assembly payloads.

**Next plan, not executed:** source-attested release/profile-matched Sassy/historical baselines; authorized exporter/semantics reconciliation; larger randomized/held-out bounded windows with per-case witnesses and candidate-stage traces; re-audit existing manuscript raw data under exact strict-PAM-anchor rules before borrowing historical plots. Only after reviewing that plan/results should any new chromosome, assembly or multi-haplotype experiment be proposed. Main merge, algorithm/default changes and production replacement remain out of scope.
