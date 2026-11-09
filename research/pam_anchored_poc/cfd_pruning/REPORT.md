# Optional downstream CFD filtering and isolated native pruning

## Status and scope

**Implemented:** optional `--min-cfd` downstream filtering in the development CLI, plus a separately built native FM-index pruning experiment. **Not implemented/enabled:** pruning in either integrated search backend, whole-genome pruning, CRISPR-Bulge, or a new default search mode. No production binary, existing index, frozen branch or main was changed. No Slurm/whole-assembly job was launched.

This follows the approved semantics and the [feasibility study](../cfd_feasibility/REPORT.md). The integrated default remains Columba + WFA2; the integrated experimental PAM-anchored search remains **unpruned**, with the same configurable `-m 5 -b 1 -z 2 -f .75 -p GG` defaults. Higher selected mismatch limits are preserved.

## 1. Downstream filter

`src/main.rs` accepts optional `--min-cfd SCORE`, finite and in `[0,1]`, **absent by default**. All reporting paths pass it to `src/reporting.rs`. After the existing per-anchor representative selection, the formatter computes CFD once and:

- retains a defined score at or above the threshold (inclusive equality, full precision);
- removes a defined score below the threshold;
- preserves unavailable scores, without promoting any different alignment;
- conservatively preserves non-finite custom-matrix results because they cannot establish that a site is below the threshold.

The comparison uses the same unrounded value that is subsequently formatted as `cf:f` to four decimals. Missing CFD continues to have **no `cf:f` tag**, not a zero tag. No new output tags, representative objective, eligibility rules or score formula were introduced. `--no-filter` continues to bypass sequence/PAM eligibility, but does not cancel an explicitly requested `--min-cfd`. Each backend retains its own prior representative policy; this study does not equate WFA2 and experimental representatives.

With standard matrices, `--min-cfd 0` and absence of the flag produced byte-identical output. The absence-of-flag path was compared against the preserved CLI built from `f7de425`; its `src/` is identical to baseline `252a6bb` before this work. **All 18 bounded configurations matched byte-for-byte**, including PAF order and tags. This is regression evidence on these panels, not a proof for every possible invocation/custom matrix.

Example (add to an existing development CLI invocation):

```bash
crisprapido ... --min-cfd 0.1
```

This flag alone **does not reduce integrated native search states**. Threshold-selected output must not be described as exhaustive under mismatch/bulge settings alone. CFD is not a calibrated probability of cleavage.

## 2. Isolated native algorithm

[`native.cpp`](native.cpp) includes the unchanged `../poc.cpp`, reuses `Search::run` for the unpruned control and `Search::emit` for checked discovery, and adds a **separate** `Discovery` traversal. The research executable uses a distinct positional protocol and TSV format; it is **not a drop-in `--pam-search-bin`**. A regression test confirms that attempting to use it through the integrated adapter fails with no PAF output.

### Bound, mismatches, positions and fallbacks

For scoreable 20-column completions, initialize the bound with the requested/observed two-base PAM factor. Walking from PAM-proximal to distal, the j-th processed **alignment column** maps to CFD position `21-j`, independently of original guide position. Multiply the bound by the actual position- and base-pair-dependent factor for that column. Matches have factor 1; no new mismatch-count cutoff is imposed.

Pure DNA/D bulges retain the **last 20 aligned columns** under the current scorer's projection, so processed suffix positions remain fixed even when distal D columns are added. Ignore processed columns beyond 20 when maintaining the bound. Pure I bulges keep 20 columns and use existing gap penalties. The position-1 gap exception is applied **before** matrix lookup, including a target N at a neutral projected position-1 gap.

Prune only if both conditions hold:

1. The conservative numeric upper bound is **strictly less** than the threshold.
2. No completion could have an unavailable score.

A mixed I/D completion for a 20-base guide is unscorable. If an event can still become mixed, keep exploring. Once the event is closed or cannot mix under the remaining size/base budget, this uncertainty can be removed. The guard can retain unnecessary paths; it must not discard unknown-score paths.

Pruning is disabled, with complete unpruned discovery/filtering, for:

- **any reference containing N** (the index does not expose a certified N-overlap summary for SA intervals);
- guides whose length is not 20;
- thresholds of zero;
- tables that fail the complete finite-`[0,1]` mismatch/gap/PAM contract.

Valid custom tables with a factor above one or a missing key were tested to fall back, not silently narrow the search or change mismatch limits. Malformed/unreadable research inputs fail clearly. The research loader/scorer is not a packaged replacement for the Rust scorer; integration would require a shared pinned table contract and broader parser compatibility review.

Every reverse product is rounded upward with `nextafter`; the comparison adds the feasibility study's conservative relative (`64 × epsilon`) and subnormal (`64 × denorm_min`) allowance. Equality is retained. Build flags explicitly disable fast-math and floating-point contraction. This assumes normal IEEE-754 execution; broader compiler/architecture and extreme custom-table numerical review remains necessary.

### Preserving the deterministic representative

Discovery stores candidate PAM anchors, not final representatives. If any paths were pruned for a guide/strand, reconstruct **all eligible alignments at every retained anchor** directly from its bounded reference span, with no CFD pruning:

- enumerate no-gap and every single contiguous `I`, `D`, `II`, `DD`, `ID`, `DI` event at every guide boundary;
- replay true reference bases on the correct strand;
- enforce the same mismatch limit, positive match count and f32 match-fraction threshold;
- minimize the unchanged `(mismatches + I bases + D bases, mismatches, compressed CIGAR, genomic start)` tuple;
- score and filter only this reconstructed representative.

If no paths were pruned, use the already-complete original representatives and skip reconstruction. Mixed/unavailable alternatives may witness discovery, but a final **defined low-scoring canonical representative** is still removed under the approved semantics.

For any qualifying or unavailable canonical representative, its path cannot be safely pruned, so its anchor is discovered. Complete local reconstruction restores the original representative even if some other paths were pruned. Final filtering therefore reproduces unpruned-then-filter output, conditional on discovery/refinement completeness. The independent finite checks below test that contract; they do not prove arbitrary genome-wide recall.

### Integrity and isolation

The research executable requires full reference/index binding **before loading** and again **before publishing the TSV results**. Its bounded tests use exclusive temporary staged references/indexes. It does not remove any production snapshot, hashing, native or post-search verification. It is not a persistent cache or an adversary-proof immutable reference session: same-user transient mutate-and-restore races remain outside the trusted-input model. Table bytes and executable/source hashes are recorded in the study evidence.

## 3. Independent validation

Authoritative evidence: [`results/summary.json`](results/summary.json). Three outputs were compared:

**A.** New integrated PAM-anchored CLI, still unpruned, with downstream `--min-cfd`.

**B.** Isolated native CFD-pruned search plus unpruned anchor-local reconstruction.

**C.** Direct-string exhaustive forward-template oracle from `cfd_feasibility/study.py`, with the independent Python CFD matrix calculation. It enumerates every strict PAM on both strands, every permitted event word and position, all eligible spans and the canonical representative **before** filtering. Native reconstruction is implemented independently in C++; it does not call the oracle. Rust reporting, native scoring and the Python scorer were checked against one another for reported representatives. All use the same empirical matrices by design.

| Panel | Reference size / contigs | Guide |
|---|---|---|
| Resolved synthetic | 13,768 bases / 598 | `GAAACAGTCGATTTTATCAC` |
| Synthetic N controls | 13,956 bases / 606 | `GAAACAGTCGATTTTATCAC` |
| CHM13 chr1:100000001–100050000 | 50,000 bases / 1 | `TAGCCACTGCACTCCAGCCT` |
| GRCh38 chr9:34000001–34050000 | 50,000 bases / 1 | `CACTGCAACCTCTGCCGCCC` |

The genomic guides are the 20 reference bases immediately preceding the first selected strict PAM boundary at local offset 40; they were extracted, not invented. Both 50 kb references are prefixes of retained 3 Mb extracts whose full FASTA hashes matched the earlier [bounded source metadata](../scalability_v2/results/3mb/regions.json). New tiny indexes were necessary for these restricted inputs; no HPRC index was rebuilt or queried. The GRCh38 parent was previously called a repeat candidate; this is not a formal repeat annotation.

**18 A/B/C configurations × thresholds `0, .1, .2, .5` = 72 successful comparisons.** Synthetic panels used `(m,z,f) = (5,0,.75), (5,1,.75), (5,2,.75), (8,2,.5), (20,0,0), (0,2,.75)`; genomic windows used the first three settings. Bulge events were enabled only when `z>0`. Results had **zero final-anchor or representative-coordinate/CIGAR differences**, no duplicates, and matching score availability/CFD within the applicable output precision.

Coverage includes mismatches, one-/two-base pure bulges, mismatch+bulge alignments, both adjacent mixed-event orders, alternative alignments, both strands, strict and ambiguous PAM controls, original-reference Ns (rather than indexed substitutions) in target alignments, and projected Ns. At `m=0,z=2,f=.75`, the resolved panel has **32 unavailable canonical representatives**; every threshold retained all 32. The N panel has **34** in that configuration, all retained. These include synthetic N edits, not a population-wide genomic sample.

Additional native-versus-oracle checks covered the other **15 PAM suffixes × four thresholds = 60 comparisons**, non-20-base fallback, custom-table fallback (above-one factor/missing key), invalid protocol and altered-reference rejection. These 60 extra comparisons do not all include the integrated CLI; the main 72 do.

Regression tests specifically cover the earlier alternative-alignment trap, threshold equality and neighboring binary64 values, a defined CFD **0** versus unavailable CFD, and a score of **0.19999999** that displays as `0.2000` but must be removed at threshold `.2`. They also cover the neutral projected gap against N and rejection of accidental research-binary integration.

Tests run on the final code:

- **164 Rust tests passed**, one pre-existing standalone-Python CFD test remained ignored because it requires external pinned source/pickle files.
- **5 native/CLI regression tests passed**.
- **9 prior feasibility/model regression tests passed**.

No A/B/C result discrepancy was observed. Build-stage issues (inherited incompatible C++ include paths; an updated test formatter call initially missing its new optional argument) were corrected before the final build. Bounded validation was repeated after reusing the formatter's score and extending native RSS measurement through output publication; the final-v4 evidence below is authoritative, not pooled with earlier runs.

## 4. Measured native performance — not model predictions

The native executable has three modes using the **same reference/index, tables, compiler and threshold**:

- `unpruned`: unchanged `Search::run`, then canonical score filtering;
- `overhead`: compute bounds/guards/counters at the **same unpruned state count**, but do not return early and do not reconstruct;
- `pruned`: actual FM-index pruning, followed by reconstruction if paths were lost.

The unpruned native control implements A's search/filter semantics and is validated against its Rust CLI output, but uses C++ scoring/TSV formatting for a like-for-like native comparison. **These are not end-to-end integrated CLI speedups**: B has not been integrated, and full CLI snapshot/Python/Rust work is not timed here.

There are **108 measured runs**: four panels × three bulge sizes × three modes × three repetitions, all at `m5,f.75,NGG,threshold=.1`. Mode order reverses on alternate repetitions. This is one shared host with warm/uncontrolled caches, one native search thread and a separate process-tree monitor. Small runtimes and three repetitions are not a robust general benchmark.

### Current two-base-bulge setting (`b1,z2`)

Times are **median milliseconds across three runs**, independently summarized by phase. Search includes actual FM traversal, locate and discovery/representative replay. Algorithm total additionally includes reconstruction, filtering, contig lookup and TSV construction; index/reference loading, binding checks, final stdout publication and the external monitor are excluded. Phase medians need not sum to the median total.

| Panel | Unpruned → pruned states | Unpruned → pruned search ms | Pruned reconstruction ms | Unpruned → pruned algorithm ms | Reported anchors, both |
|---|---:|---:|---:|---:|---:|
| Resolved synthetic | 425,148 → 291,423 | 85.099936 → 74.976577 | 15.002871 | 86.330783 → 91.308192 | 176 |
| N synthetic | 427,944 → 427,944 | 85.537228 → 87.918747 | 0 | 86.761551 → 89.149655 | 180 |
| CHM13 50 kb | 610,313 → 294,716 | 55.131933 → 41.789465 | 1.181248 | 55.223578 → 43.066048 | 24 |
| GRCh38 50 kb | 713,591 → 348,562 | 65.689357 → 49.562037 | 2.315961 | 65.806571 → 52.005667 | 19 |

The dense resolved synthetic case explored **31.454% fewer indexed states but was 5.766% slower** over the measured algorithm scope: 558 discovered anchors required local reconstruction before filtering to 176 final sites. The two genomic windows explored about **51% fewer states** and had **22.015% / 20.972% lower algorithm times**, respectively, on these particular runs. N fallback pruned no states and was **2.752% slower**. None of these is a genome-wide forecast.

With bound evaluation enabled but early returns disabled, median algorithm times were **102.100697 ms** (resolved synthetic), **89.264072 ms** (N fallback), **91.150215 ms** (CHM13) and **105.988394 ms** (GRCh38). Comparing these with their unpruned totals measures the added traversal/bound bookkeeping at fixed state count, not the refinement cost or a per-operation CPU cost. The N control skips factor updates entirely but retains some structural traversal overhead.

### Zero-/one-base bulges, resolved panels

| Panel / z | Unpruned → pruned states | Unpruned → pruned algorithm ms | Final anchors |
|---|---:|---:|---:|
| Synthetic / 0 | 9,702 → 2,882 | 3.045817 → 1.265813 | 158 |
| Synthetic / 1 | 115,295 → 11,983 | 18.295278 → 5.080014 | 170 |
| CHM13 / 0 | 15,862 → 2,262 | 1.439775 → 0.408229 | 24 |
| CHM13 / 1 | 226,985 → 9,625 | 20.402549 → 1.812736 | 25 |
| GRCh38 / 0 | 18,305 → 3,131 | 1.728646 → 0.480867 | 18 |
| GRCh38 / 1 | 265,529 → 11,928 | 23.965293 → 1.825681 | 19 |

Unlike the previous Python occurrence-set study, **all state counts and search timings in this report are measured native FM-index values**. Reconstruction trials are recorded separately and must not be confused with indexed states. Pure one-base events cannot become mixed I/D, permitting earlier pruning, but that is not a reason to change the default two-base allowance. N panels remain unpruned at all tested sizes.

### Memory

Peak native process RSS comes from Linux `getrusage(RUSAGE_SELF).ru_maxrss`, measured after result publication. It excludes the verifier and monitoring processes and is not aggregate memory. The additional process-tree metric samples summed RSS every **20 ms**, including the supervisor/verifier; it can double-count shared pages and miss short peaks.

At `z2`, maximum values over each mode's three repetitions were:

| Panel | Native peak KiB, unpruned / pruned | Sampled process-tree peak KiB, unpruned / pruned |
|---|---:|---:|
| Synthetic | 7,484 / 7,536 | 29,876 / 29,800 |
| N | 7,460 / 7,496 | 29,900 / 29,880 |
| CHM13 | 7,452 / 7,540 | 32,824 / 29,912 |
| GRCh38 | 7,524 / 7,520 | 30,192 / 32,732 |

No material native memory saving is established on these tiny inputs. A 2 GiB **per-process virtual-address** cap and 1.5 GiB sampled process-tree stop were used; neither fired. These are not a cgroup or guaranteed aggregate memory limit. They do not replace the documented safeguards for any separately authorized future large run.

## 5. Reproduction and artifacts

Run from `/export/local/home/salehi/crisprapido-pam-validation` on the development branch. The site-specific script below uses a **new exclusive external output root**, pinned compiler profile, locked/offline Cargo, separate target/native binaries, existing pinned WFA2 libraries, and only bounded synthetic/50 kb reference indexes:

```bash
bash research/pam_anchored_poc/cfd_pruning/reproduce.sh /path/to/NEW_ABSOLUTE_OUTPUT_ROOT
```

The research harness needs Python 3.9+; set `RESEARCH_PYTHON` if `python3` resolves to an older interpreter. The build-pinned binding verifier can remain `/usr/bin/python3`. `BOUNDED_REAL_ROOT` may identify the retained 3 Mb source extracts instead of `/tmp/crispr-v2-3mb`; their hashes must match the recorded metadata. If those extracts no longer exist, recover/verify only the documented bounded inputs before reproducing; do not substitute another assembly. The script is provided for reproduction and was syntax-checked; the final experiment used its component commands as described in the retained build/test evidence, rather than executing a further duplicate full run.

To run only the isolated native/CLI regression suite, set `CRISPR_CFD_CLI`, `CRISPR_CFD_NATIVE`, `CRISPR_CFD_UNPRUNED` and `COLUMBA_BUILDER` to the separate development/native/control/builder paths, then:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover \
  -s research/pam_anchored_poc/cfd_pruning -p 'test_*.py' -v
```

Retained small evidence:

- [summary](results/summary.json): all A/B/C comparisons, 108 raw native timing/state/RSS records, monitor summaries, input and binary hashes;
- [native build](results/native-build.json): exact command, compiler/library/source hashes; GCC 12.3.0, `-O3 -fno-fast-math -ffp-contract=off`;
- [source provenance](results/source-provenance.json), [Rust tests](results/rust-tests.log), [native tests](results/native-tests.log), [model regressions](results/model-regressions.log).

Isolated runtime root: `/export/local/home/salehi/crisprapido-cfd-pruning-dev-v1/`. Final native: `native-build-v4/native`; final development CLI: `target/release/crisprapido`; authoritative full logs: `validation-final-v4/`. That validation took **75.546 s** overall, including preparation, oracles, CLI/native checks and performance runs; it is not a search-only timing. Original CLI and native payload hashes remained `9dac869a…afd1ce` and `ce5a4f7b…69269d`. No genomic FASTA/index or large PAF was committed. Columba source/build attestation and packaging limitations from earlier reports remain.

## 6. Stop point and remaining work

**Pruning is deliberately not wired into the integrated backend.** Downstream filtering is available in freshly built development CLIs; absence of `--min-cfd` preserves prior behavior. The original unpruned search remains available with or without filtering.

Before considering integration: review native code and numerical assumptions independently; broaden parameter/randomized/boundary oracles and custom-table compatibility; validate larger bounded guide batches and output volumes; measure actual full-workflow costs with matched scopes; and implement a completeness-preserving N-aware strategy if useful. At present, an HPRC assembly containing Ns would receive **no traversal reduction** under the safe fallback. Any future representative reconstruction must retain the established comparator, not maximize CFD. Genome-wide recall and biological sensitivity remain unproven; CRISPR-Bulge requires its own separate scoring/support/bound validation.
