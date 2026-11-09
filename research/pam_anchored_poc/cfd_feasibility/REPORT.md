# CFD-guided indexed searching: bounded feasibility study

## Decision summary

**A safe CFD upper bound is feasible for defined scores under a validated multiplicative scoring contract. It is not sufficient by itself for a drop-in pruning change.** The main obstacles are preserving the current representative alignment and retaining unavailable scores, especially mixed I/D paths and reference-N cases. The recommended future design is optional **score-aware candidate discovery followed by unpruned, anchor-local representative reconstruction**, with conservative unavailable-score fallbacks. Start on small resolved references; do not enable this for HPRC by default.

This study changed **no production/search/CLI/scoring code or defaults**, and ran no whole-genome search, Slurm job or CRISPR-Bulge implementation. `study.py` is an isolated Python occurrence-set model, **not a modified FM-index backend**. It tests the proposed mathematics/semantics; native performance remains unmeasured. Existing experimental limits and higher user-selected mismatch limits remain unchanged.

## 1. Existing scorer and search

Source reviewed at development baseline `7c4e77f`:

- [`src/cfd_score.rs`](../../../src/cfd_score.rs): `calculate_cfd`, `get_cfd_score`, reconstruction and matrix parser.
- [`poc.cpp`](../poc.cpp): `Search::run`, `walk`, `emit` and deterministic representative ordering.
- [`src/pam_backend.rs`](../../../src/pam_backend.rs), [`run_snapshot.py`](../integration/run_snapshot.py), and [`reporting.rs`](../../../src/reporting.rs): integration, mandatory snapshot checks and downstream scoring/formatting.
- Independent scorer: [`audit_emitted.py`](../../../scripts/release_validation/audit_emitted.py). This is separate from the Rust implementation, but uses the same score tables by design.

For a supported alignment projected to exactly 20 columns, the existing model is

```
CFD = PAM_factor(observed final two PAM bases) × product(column factors, positions 1..20)
```

Positions are **guide-oriented, PAM-distal 1 to PAM-proximal 20**. Matches contribute 1. For mismatches the key uses RNA guide base (T→U), **complement of the target base**, and the scored column position. Scores depend on substitution identity and position, not just mismatch count. A mismatch can have factor 1 or 0: there is no generally safe fixed penalty per mismatch or CFD-derived mismatch cutoff.

The supplied tables contain **393 mismatch/gap entries and 16 PAM entries**, all finite in `[0,1]`; 51 mismatch/gap and nine PAM factors are zero. There are 14 mismatch/gap factors equal to one. `GG` has factor 1. The scorer treats matches neutrally without consulting a matrix entry. [Matrix summary](results/matrices.json).

The normal CLI permits alternate matrix files. Its current parser does **not** enforce the stricter finite/unit-interval/completeness contract needed here; a future pruning path must validate the actual loaded tables, bind their digest/values across Rust and native search, and fail clearly or fall back to unpruned search if unsupported. Do not change ordinary scoring behavior silently. Factors above 1 would invalidate the simple suffix-product bound; supporting them would require a separate future-factor-maxima bound. Missing keys are unavailable scores, not zero penalties.

### Bulge and ambiguity semantics — not all bulges are unscorable

For the current resolved **20-base guide** case:

| Alignment class | Current CFD behavior |
|---|---|
| No gaps | Standard 20-column product |
| Pure I (guide base aligned to target gap) | Still 20 columns; gap factors are supported |
| Pure D (extra target/DNA base) | Remove the total D count from the **PAM-distal prefix of both aligned strings**, keeping the last 20 columns; retained gaps/mismatches use projected positions |
| Adjacent mixed I/D | No pure-D projection; 20+D columns, hence `None` |
| Retained target N | Missing matrix key, hence `None` |
| N removed by pure-D projection | Can yield a defined score; not universally `None` |

A gap at scored position 1 is explicitly neutral. A D event can remove a distal mismatch or move a retained column into a different matrix position; naïvely multiplying factors at the original guide indices is unsafe. For non-20-base guides, score support cannot be inferred from this table (some mixed-gap lengths can become 20); conservatively disable score pruning. Malformed input must continue to fail normally, rather than masquerading as an unavailable biological score.

`Search::run` seeds the full indexed PAM and traverses PAM-proximal to distal on both strands (`BACKWARD` plus, `FORWARD` minus). The required PAM suffix is already known, so its factor can initialize the bound; the resolved wildcard base has no separate CFD penalty. On minus traversal, complement each indexed target character before looking up the guide-oriented factor. Mismatch/bulge/fraction eligibility remains an independent constraint at every user setting.

## 2. Safe branch bound

Let `r` be the number of alignment columns already traversed **from the PAM**. For a completion that is scoreable under the 20-column conventions above, its j-th column from the PAM has scored position **`21-j`**. This stays fixed even when a future distal D event changes the total alignment length: pure-D projection always retains the **last 20 columns**. Already traversed columns beyond 20 are discarded for pure D. Pure I also has exactly 20 columns. Thus no speculative guide-index shift is needed if factors are attached to reverse **column count**, not guide-base count.

For a branch with known suffix columns, define

```
U = PAM_factor × product(factors for traversed columns j=1..min(r,20))
```

For each **defined-score completion**, all unvisited scored factors lie in `[0,1]`, so `CFD_completion <= U` over the exact real values of the loaded factors. Ignoring mismatch/gap/fraction constraints only makes this upper bound looser; it never introduces a lower mismatch limit. For score threshold `tau`, prune only when **`U < tau`**, not at equality.

This is a bound on **defined alignment scores**. `None` is not an ordered numeric value. The actual safe prune condition is:

```
upper_bound_on_defined_scores < tau
AND no eligible completion could have an unavailable score
```

If a gap event is still open or can still begin, retain branches that can produce mixed I/D. With `z=2` and a remaining guide base, an unused event can still become I+D or D+I. This blocks early pruning even on resolved references. Once a pure event is closed (or the budget precludes mixing), this particular uncertainty disappears. The model's conservative state guard implements this distinction and retains all mixed paths.

### Numerical safety

The mathematical inequality alone does not guarantee a safe implementation with binary64 multiplication in a different order. The model rounds each reverse product upward using `nextafter` and adds a conservative final relative allowance (`64 × machine epsilon`) and subnormal allowance (`64 × smallest subnormal`). This covers the ordinary IEEE-754 forward scorer's at-most-21 multiplications; equality and near-threshold cases are retained. Tests also compare directed-up products against exact rational products of the loaded binary64 factors. A native implementation must independently verify these numerical assumptions, avoid unsafe fast-math/reassociation/flush-to-zero behavior, and use the same matrices and inclusion threshold as the full scorer. Comparing to four-decimal PAF text is not sufficient.

### Reference N is a separate blocking issue

Columba substitutes Ns in its text; an SA interval can mix true resolved targets with true N-containing targets. The search currently lacks a certified N-overlap summary for an interval. A low score computed from substituted letters must **not** delete an occurrence whose original-reference score is unavailable.

The tested baseline therefore **disables all score pruning if the reference contains any N**, while preserving ordinary eligibility and final threshold filtering. This is safe but yields no search-state reduction on the N panel and would yield no reduction on the existing HG00735_pat bundle. Improving it requires a separately validated interval/position-aware N-overlap mechanism or a complete independent N-overlap fallback search. Pre/post binding checks do not provide that interval annotation. An N in a PAM still fails strict PAM eligibility; an N projected away need not be unscorable. Never bypass reference verification or replay to work around this limitation.

## 3. Optional threshold semantics and representative preservation

Proposed future option, **not added to the CLI**: `--min-cfd TAU`, finite `TAU ∈ [0,1]` for the validated standard-table mode. Default **absent/disabled**, not a new nonzero cutoff. Zero retains all eligible sites, including defined zeros and unavailable scores. Use full-precision scores, inclusive `score >= TAU`; retain unavailable representative scores with an explicit status/reason in any future output design. A threshold-selected result is **not exhaustive output under mismatch/bulge settings alone**.

The compatibility target should be:

1. Run the current unpruned search and choose the current deterministic representative at each anchor.
2. Keep that representative if its CFD is unavailable or at least `TAU`.

This is deliberately **not** “pick the highest-CFD alignment” or “keep an anchor if any alternative scores highly.” Those would change scientific/output semantics. If a defined low-scoring canonical representative has an unavailable alternative, the rule above still filters that site: it is not an unavailable *reported representative*. This distinction must be agreed with the user before implementation. All eligible unavailable-score representatives remain retained.

Simply pruning and applying the old comparator to surviving paths violates this target. The study found a concrete example at `I14_1-`, minus anchor 3, `m5/z1/f.75`, `TAU=.1`:

- Current representative: `11=1I8=`, CFD **0.0307692307692**; should be filtered.
- Naïve pruned-search representative: `10=1I1X8=`, CFD **0.2215384615384**; incorrectly emitted.

The higher-CFD alignment is not the better representative under the declared ordering. With `z=2`, unavailable mixed alternatives also replace pruned low-score representatives, magnifying this issue. [Counterexamples](results/summary.json).

### Recommended future architecture

- First establish the optional post-selection threshold contract against unchanged search; no traversal pruning yet.
- Score-aware **discovery** may discard a path only under the safe condition in section 2. Record candidate anchors reached by qualifying or potentially-unavailable completions; carry uncertainty explicitly.
- For every retained candidate anchor, independently enumerate **all eligible alignments** in its bounded local reference span under the same snapshot, mismatch/bulge/f32-fraction rules and exact existing comparator, **without CFD pruning**. Then score/filter the chosen representative.
- A qualifying/unavailable canonical representative necessarily witnesses discovery of its anchor. Refinement restores representatives lost from discovery, and removes anchors witnessed only by high-score/unknown alternatives whose canonical representative is defined and low. This proves equivalence to the unpruned-then-filter target, conditional on complete discovery/refinement.
- Keep reference/index identity checks, atomic result publication and unsupported-mode rejection intact. Pin score tables as part of the run. Validate targeted refinement against current representatives before considering deployment.

The local span is small (current guide length ±2) but refinement can still dominate for many candidate anchors. Existing Python finite-oracle refinement here uses precomputed paths; it is **not** a production anchor verifier or a measured estimate of its cost.

## 4. Bounded independent comparison actually performed

Commands and artifacts are below. One synthetic guide `GAAACAGTCGATTTTATCAC`, not an HPRC experiment. The resolved panel has **598 contigs / 13,768 bases**; the N panel adds both-strand target-N, projected-distal-N and invalid-PAM-N controls, totaling **606 contigs / 13,956 bases**. Resolved targets exhaust all 256 suffix 4-mers on both strands, plus finite pure/mixed-gap, repetitive and 16-PAM controls. This is a deliberately enriched synthetic panel, not representative genomic sequence.

Three independent layers:

1. **Forward direct-string oracle** enumerates every strict PAM on each strand and every contiguous gap word (`I`, `D`, `II`, `DD`, `ID`, `DI`) at every position, then replays and scores all eligible alignments. This is exhaustive for each supplied finite reference/configuration.
2. **Unchanged native Columba search** uses freshly bound tiny indexes, retains all mandatory binding checks, and is compared to oracle anchor sets and exact representative coordinates/CIGARs/mismatch counts. **25 native comparisons passed**: five NGG configurations on each panel plus the other 15 PAM suffixes on the resolved panel. Tested `(m,z,f)` are `(5,0,.75)`, `(5,1,.75)`, `(5,2,.75)`, `(8,2,.5)`, `(20,0,0)`. Higher mismatch settings were not capped to five.
3. **Reverse occurrence-set model** groups explicit occurrences by indexed target letters instead of using FM-index intervals. It uses a deterministic N→A shadow, true-reference replay, and the proposed bound/unknown guard. Every unpruned alignment path matched the forward oracle. At thresholds **0, .01, .1, .5, 1**, all **125 panel/configuration/threshold comparisons** retained every qualifying or unavailable oracle alignment path. Candidate refinement and final filtering had **zero missing/unexpected anchors and zero representative differences** versus oracle unpruned-then-filter output.

The forward oracle and reverse model use different path enumeration but share the replay/eligibility helper and existing independent Python CFD oracle; they are not wholly independent implementations of every check. Native representative comparisons and the Rust score probe provide additional implementation-independent checks. The native executable has documented measurement-only timers, not changed search transitions. The native builder uses its actual N substitutions; equality to the true-reference oracle does not require it to use the model's N→A shadow. The model is contig-local and excludes cross-contig interval artifacts that native emission rejects. Its states are **not native FM-index states**.

All **15,478 distinct** oracle scoring inputs were also checked with the unchanged Rust `get_cfd_score`: **4,861 defined scores matched, 10,617 unavailable scores agreed, zero differences** (probe reports 12 decimals; numerical comparison tolerance `1e-9`). These include alternative paths, not just one representative per anchor. Nine unit/regression tests passed, including exact-score threshold equality, matrix rejection, pure-D suffix bounds, the representative-change counterexample, mixed-only unavailable representatives at `m=0`, N fallback, non-20-base fallback and higher mismatch configurations. [Tests](test_study.py), [test log](results/tests.log).

These are new small feasibility checks, not repetitions of completed whole-assembly experiments. No native pruned search was built or timed. Finite agreement is not a proof of arbitrary whole-genome recall, nor proof of biological utility of a CFD threshold.

## 5. Measured model costs and likely benefit

At `TAU=.1`, on the **resolved synthetic panel**:

| m / z / f | Unpruned model states | Discovery model states | Reduction | Median Python unpruned / discovery seconds | Discovery / final anchors |
|---|---:|---:|---:|---:|---:|
| 5 / 0 / .75 | 4,748 | 1,417 | 70.156% | .035810 / .014314 | 158 / 158 |
| 5 / 1 / .75 | 54,623 | 5,627 | 89.698% | .304505 / .046683 | 172 / 170 |
| 5 / 2 / .75 | 203,786 | 142,556 | 30.046% | 1.442522 / 1.148075 | 558 / 176 |
| 8 / 2 / .5 | 277,129 | 147,684 | 46.709% | 2.013679 / 1.189526 | 558 / 176 |
| 20 / 0 / 0 | 5,004 | 1,418 | 71.663% | .038363 / .014365 | 158 / 158 |

Timings are medians of three alternating model calls per threshold/configuration on the same shared host, not isolated native search or full-workflow timings. Explicit occurrence grouping, Python scoring and validation differ from Columba costs; refinement, index setup, snapshotting and binding are **excluded**. They establish a plausible state-reduction mechanism, **not a measured CRISPRapido speedup**. At `m5/z2/f.75`, `.5` and `1` thresholds reduce states only about 30.7%, because potentially-unavailable mixed paths must survive. At `.1`, naïve deduplication without refinement differs from the specified result at **382 anchors** in that configuration.

The N panel has **zero state reduction at every threshold**, with small observed overhead. Its two unavailable canonical representatives remain retained even at threshold 1. Other unavailable alternatives are preserved during discovery; separate mixed-only unit tests cover unavailable canonical bulge cases. Zero threshold made no state-count change on either panel.

For perspective only, the saved optimized single-guide job spent about **80.6/645.6 s** in native search, and the three-guide batch about **257.5/837.9 s**. Even eliminating all native search could save at most roughly **12.5% / 30.7%** of those historical totals if everything else were fixed; real pruning adds work and cannot remove all search. This is a hypothetical fixed-overhead ceiling, not a forecast. Larger batches may amortize setup, but are not tested here. For the actual N-containing HPRC reference, the tested conservative fallback produces **no traversal benefit** until N-aware handling is developed.

## 6. Reproduction, provenance and limits

Run the small unit tests (no index or HPRC access):

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover \
  -s research/pam_anchored_poc/cfd_feasibility -p 'test_*.py' -v
```

The study used an isolated probe compiled from the existing `compatibility/cfd_probe.rs` against the already built `f7de425` library. Its CFD source SHA256 equals current `src/cfd_score.rs`; archived/current `poc.cpp` also match. Compilation requires the **absolute pinned rustc 1.94.1**; an initial PATH-selected rustc 1.70.0 rejected library metadata, then the absolute-compiler build succeeded. To rebuild the isolated probe, set `RUSTC` from `results/tools.json`, export its `build_environment` values (`CRISPR_LIB`, `DEPS_DIR`, `RUSTC_EXTRA_FLAGS`), choose a new `OUTPUT` path, and run `bash research/pam_anchored_poc/compatibility/build_cfd_probe.sh`. No existing executable was overwritten. [Tool provenance](results/tools.json).

From the checkout root, reproduce into a **new** results directory:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 research/pam_anchored_poc/cfd_feasibility/run_study.py \
  --out /path/to/new-small-results \
  --builder /export/local/home/salehi/crisprapido-index-recovery-v1-2026-10-28/rebuilt_HG00128_hap2_attempt02/tools/columba_build \
  --native /export/local/home/salehi/crisprapido-pilot-rebuild-f7de425-v1/native/wrapper \
  --cfd-probe /export/local/home/salehi/crisprapido-cfd-feasibility-tools.XanPvlMD/cfd_probe
```

The study imposes a 2 GiB **per-process virtual-address** cap, uses subprocess timeouts and temporary synthetic indexes, and leaves only small logs/JSON in Git. It completed in **234.238 s** on the recorded shared host. [Full results and hashes](results/summary.json). This is not aggregate memory profiling or a native performance benchmark. No HPRC input or production binary was modified, and no experiment output from prior jobs was replaced.

**Recommendation before implementation:** obtain approval for canonical-representative threshold semantics and retained-unknown behavior. Then implement/test optional downstream filtering and a separately isolated resolved-reference native prototype, with identical unpruned controls, targeted-representative equivalence and numerical review. Do not activate pruning on N-containing bundles until a completeness-preserving N strategy is independently validated. Keep default mode unpruned, mismatch settings unrestricted within the current backend contract, and CRISPR-Bulge integration separate. CFD bounds cannot be transferred to a different scorer without deriving that model's own bounds and score-support rules.
