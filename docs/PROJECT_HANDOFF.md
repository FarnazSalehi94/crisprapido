# CRISPRapido project handoff

## MooseFS incident, reconciliation and recovery pilot (2026-11-13)

[RECOVERY_PLAN.md](../research/hprcv2/sbds_nrr/RECOVERY_PLAN.md) and
[INCIDENT_REPORT.md](../research/hprcv2/sbds_nrr/INCIDENT_REPORT.md).
Ledger corrected: the previous "20+10+438=468" double-counted HG02055_mat
and NA18906_mat. Authoritative manifest v2 has **exactly 466 unique,
mutually exclusive** assemblies: 11 READY_BOUND + 10 TRUSTED_EVIDENCE +
7 INDEX_EIO + 438 UNATTEMPTED. Incident scan of all ten failed sets:
**persistent, reproducible errno-5 unreadable regions** (not transient, not
missing files, **not corrupted assemblies** — every source FASTA.gz verifies
with valid CRC/SHA256). Damage is distributed across the MooseFS index
store and grew since 2026-10-28; admin-facing report with paths/offsets/
reproducer is ready. Controlled rebuild pilot (3 assemblies, onto healthy
/export, jobs 2871822): **~16–17 min build, ~40–42 GB peak RSS, 11.1–11.55 GB
index each**, all components hashed; all three bound (2871828–30) and
**HG03784_hap1 passed a complete four-PAM verification search** (2871832,
30 min, 192,353 audited hits, full scoring). Campaign plans prepared but NOT
launched: verification ≤5.3 TB reads / ~19 CPU-h; rebuild scenario ~190 sets
≈ 55–60 rebuild-days at %4, ~2.9 TB on /export; do NOT rebuild onto MooseFS
before cluster-team diagnosis (replica recovery may avoid most rebuilds).

## SBDS deep-verification batch 1 (2026-11-13; campaign not launched)

[VERIFICATION_BATCH1.md](../research/hprcv2/sbds_nrr/VERIFICATION_BATCH1.md):
16-assembly controlled Stage-1 batch (array 2871790, %4 concurrency, one
CPU, octopus06/07/08; octopus05 excluded for 16 G scratch). **6 deep-verified
bindings committed; 10 EIO failures, originals untouched.** Measured per
successful assembly: 12.16 GB read, 147–172 s; committed bindings ~3.06 GB
each; reusability recheck passes. EIO prevalence: 6/7 historically notable
assemblies (all four historical no-perfect-match cases, HG00321_hap2,
HG03784_hap1) plus 4/9 of an ordinary spread ≈ **44% of ordinary assemblies
may need index rebuilds**. Ledger v2: **20 of 466 currently search-ready**
(12 trusted + 8 committed bindings), 10 EIO recorded, 438 unattempted.
Remaining-campaign scenario: ≤5.3 TB reads, ~19 CPU-hours, ~1.4 TB bindings —
**safe to proceed with verification in batches**; a complete search also
needs an authorized, separately piloted **index-rebuild campaign** and a
cluster-team escalation of the MooseFS degradation. Rebuild cost unmeasured.

## SBDS 466-genome readiness (2026-11-13; cohort NOT launched)

Infrastructure blockers resolved and a controlled multi-assembly test passed:
[COHORT_READINESS.md](../research/hprcv2/sbds_nrr/COHORT_READINESS.md).
Ledger audit: **12 assemblies** have trusted full-content index evidence
(incl. GRCh38/CHM13); **454 require deep verification** (~5.6 TB component
reads, ~20-35 CPU-hours at measured rates; EIO assemblies will stop with
explicit failures, never zero-filled). Expansion tooling rewritten and
unit-tested (restart-safe checkpoints, tamper detection, process-group kills,
per-assembly pins; 6 regression tests). Two new bindings (HG02055_mat,
NA18906_mat) committed via the repaired trusted binding helper. Cohort test
array 2871780 (one CPU, 56G): 192,378/192,148 total hits per assembly, all
CFD-defined, 42,465/42,389 CRISPR-Bulge-defined (NRR scorer byte-identical
to the officially validated bundle), 149,913/149,759 explicit N/A, all eight
(assembly,PAM) audits PASS. Searches ~28-30 min/assembly on one CPU; full
466 scenario ≈233 one-CPU job-hours, ~4.75 GB compressed outputs, bindings
~1.45 TB on /export. Remaining blockers: the 454-assembly deep-verification
campaign (needs I/O budget + EIO-rebuild authorization), noncanonical
biological-validity caveat, historical EIO/no-perfect input investigations.
Recommended submission settings are in the readiness report; nothing submitted
for the cohort.

## SBDS resumed: NGG/NRR pilot completed; NRR scoring validated; cohort still gated

Latest SBDS request supersedes the PAM pause below: standard SpCas9 NGG is
confirmed; NRR coverage is requested for alternative PAMs. See
[NGG/NRR report](../research/hprcv2/sbds_nrr/REPORT.md) and the
[saved-pilot NRR reannotation](../research/hprcv2/sbds_nrr/REANNOTATION.md).
Four unchanged literal-PAM searches on verified HG00735 completed under
one-CPU Slurm array 2871759: NGG 38,060; NAG 46,358; NGA 42,918; NAA 64,712.
All 192,048 hits independently audited; all CFD defined. NGG both-score output
is byte-identical to prior primary.

**Saved-pilot reannotation (searches NOT repeated):** the Python scorer now
supports NRR target suffixes within the unchanged mismatch/gap domain; guide
NGG / actual target PAM. Neural defined **42,425** total: NGG 8,104, NAG
10,452, NGA 9,421, NAA 14,448 (**+34,321 noncanonical**). Remaining 149,623
retain N/A with reasons (>24 columns 70,474; gapped mismatch domain 69,176;
multibase RNA 9,973). All **42,425** supported predictions independently match
official inference (max diff 2.2e-9), plus 8,272 synthetic domain fixtures
(max diff 6.1e-13); categorical features hash-equal. Default CFD pass-through
byte-identical on all four inputs. Reannotation job 2871764 (103.7 s, peak
0.944 GiB sampled); domain job 2871765. **Noncanonical biological validity is
NOT established**; scores are experimental. Full 466-genome rerun remains
BLOCKED and explicitly stopped pending further authorization.

Root: `/moosefs/salehi/HPRCv2_CRISPRapido/runs/sbds_ngg_nrr_v1/`.
Serial workflow 28.7 minutes, maximum sampled tree RSS 13.64 GiB; naive 466-fold
scenario 222.9 one-CPU job-hours and 4.75 GB compressed PAF+annotations, excluding
binding preparation and diagnostic copies. No 466-genome rerun or population
figures yet: complete trusted binding ledger and safe tested expansion tooling
are missing. Earlier claims below that draft full-run helpers are ready are
**withdrawn**: binding helper contains closed-file fsync/string-path/digest
errors; full-run worker has chunk stream/status/checkpoint defects. Do not use
them for this expansion until repaired and failure-tested. Pilot used only the
already verified HG00735 binding. Brunello remains separate and unlaunched.

## Current development state

- Repository: https://github.com/FarnazSalehi94/crisprapido
- Checkout: `/export/local/home/salehi/crisprapido-pam-validation`
- Development branch: `feature/pam-anchored-search-poc`.
- Whole-assembly experiment milestone: **`dfd1636`** (three-guide batch completed); handoff documentation baseline: **`7c4e77f`**; CFD feasibility milestone: **`252a6bb`**. The subsequent downstream filtering/isolated native pruning implementation is linked below; pruning is **not enabled in the integrated backend**. The earlier handoff ending at `bfbb8ac` is superseded; **do not submit the three-guide batch as an unperformed experiment**.
- Experimental backend is integrated but is neither the default nor a production release. The original **Columba candidate search → WFA2 verification → PAM filtering → CFD → PAF** backend remains the default and is preserved. Earlier original-versus-Columba/WFA2 benchmarking is separate from the experimental backend results below; do not repeat it without a scientific reason.
- `--search-backend pam-anchored` uses Columba FM-index primitives, seeds the PAM constraint, tracks mismatches/bulges, handles both strands and reference Ns, and reports one deterministic representative per PAM anchor without WFA2. CFD scoring remains downstream.
- Keep CLI-configurable settings **`-m 5 -b 1 -z 2 -f 0.75 -p GG`**; `-p GG` means strict NGG. Do not change defaults or representative selection without discussion. Representatives minimize `(mismatches + inserted bases + deleted bases, mismatches, compressed CIGAR, start)`, not CFD.

## Active HPRCv2 SBDS / Brunello request (2026-11-11)

See [independent project audit and guarded pilot plan](../research/hprcv2/README.md)
and [downstream allele/IMPG requirements](../research/hprcv2/ALLELE_PLAN.md).
Separate run trees now exist under `/moosefs/salehi/HPRCv2_CRISPRapido/runs/`.
Initial audit Slurm 2871688: frozen 464+2 cohort and 77,441-guide library verified;
two historical SBDS PAFs and ten Brunello small-artifact checks currently return
EIO. Earlier successful historical QC must not be described as current full
readability. Prior resolved index evidence covers 12 sets, not all 466.

The requested SBDS primary search is m5/b1/z2/f0.75/GG with both scores and no
CFD filter. Historical runs used m2/b1/z1, so a separate matched-domain control
is needed. Safety pilots on the already bound HG00735_pat assembly are now
complete: SBDS primary/control Slurm **2871689/2871691**, **38,060/61** eligible
hits, **8,104/61** neural defined. Brunello A=m4/b1/z1 and B=m5/b1/z2 share three
exact frozen-library GC-stratified guides: jobs **2871690/2871692**, **13,474/70,364**
hits, **13,474/14,924** neural defined. All **121,959** rows independently audited;
**95** sampled predictions match official inference within 1e-6 (jobs 2871693/4).
See [SBDS pilot](../research/hprcv2/SBDS_PILOT.md) and
[Brunello review report](../research/hprcv2/BRUNELLO_PILOT.md), including actual
resources/compression and explicitly unvalidated 76k/874k native-core-hour linear
full-library scenarios. Forty-eight Python tests passed; core Rust unchanged.

On HG00735, all **57/57** strict-eligible historical CIGAR witnesses recovered;
old raw had **756 non-NGG/813**. Prior chr12/17 patterns have **GGC/AGT** PAMs
and are correctly excluded. Historical gapped CFD discrepancies also remain:
25 original eligible CIGAR scores differ from independent current reconstruction;
24 shared unchanged representatives have changed CFD. Do not attribute every
score change to representatives or claim new biological calibration. Brunello B
contains every A anchor but changes 303 representatives/neural availabilities;
thresholded A/B views are not nested despite nested eligibility.

**No full Brunello launch is authorized before user review.** Full SBDS expansion
awaits remaining integrity/binding gates; no 466-task search array has been
submitted. Two historical SBDS PAFs and ten Brunello small-artifact checks return
EIO; 454 index sets still lack full-content evidence here. No completed
population rerun/updated six SBDS figures, IMPG mapping or allele-frequency
results are claimed yet. External figures are clearly labelled n=1/three-guide
pilot diagnostics. Completed prior HG00735 experiments below are unchanged and
must not be repeated as if new.

## Brunello expanded benchmark and full-run preparation (2026-11-12; SBDS paused)

**SBDS is PAUSED pending confirmation of the prime editor's PAM; no SBDS
search/parameter work is in progress.** The completed bounded benchmark is in
[brunello_expanded_benchmark/REPORT.md](../research/hprcv2/brunello_expanded_benchmark/REPORT.md):
16 timing-blind stratum-selected guides x 2 verified assemblies (HG01496_pat
safely bound this session, 371 s trusted import) x 4 configurations, plus a
1,000-guide library-mix scaling benchmark at m3 and m4; **11.7 M rows audited
with zero errors; nothing full-scale launched**. Key measured results: m3/b1/z1
≈ 1.45 s/guide, 523 rows/guide, 100 % both scores; m4/b1/z1 ≈ 8.3 s/guide,
5,855 rows/guide, 100 % both scores, and it adds ~100 CFD>=0.20 sites per guide
that m3 omits with *perfect nesting* (zero representative/CFD changes at shared
anchors); m5/b1/z2 ≈ 96 s/guide, 83 % neural N/A, 2,667 representative changes
and fewer additional high-CFD sites than m4 — rejected. Full-HPRCv2 campaign
estimates: m4 ≈ 128,000 CPU-hours/~19.4 TB; m3 ≈ 19,000 CPU-hours/~1.7 TB;
m2 rejected (misses high-CFD sites). Recommended scientific config: candidate-A
(m4); pragmatic staged alternative: m3 (forward-compatible because nesting is
perfect). All four configurations are fully compatible with -f 0.75 (the
match-fraction rule is matched/20 and clips nothing at m<=5).

Prepared but NOT launched (Task 6 of the request): prepare_binding.py (in-place
trusted binding for the 12 recovery-verified assemblies; 452 others need the
verification campaign first), make_full_plan.py (requires a recorded user
authorization token), full_run_worker.py (chunked, checkpointed, restart-safe,
1 CPU/assembly, full audits, explicit N/A), submit_full_run.slurm (464-task
array template), verify_full_run.py (ledger; never zero-fills). Storage/filtering
strategy: keep complete eligible both-score PAF + annotations compressed as
authoritative; CFD thresholds are derived views only and can never define
allele-frequency denominators. IMPG/Figure-1C needs remain specified in
research/hprcv2/ALLELE_PLAN.md; nothing mapped or plotted yet.

## Performance investigation (2026-11-11, bounded; no full reruns launched)

See [performance report](../research/hprcv2/performance/REPORT.md). The apparent
7.5-min-for-one-guide discrepancy is ~310 s of fixed per-invocation overhead
(snapshot copy, pre/post full verification, index/reference loads) that batching
eliminates, plus the genuinely expensive requested m5/b1/z2 domain (~101 s/guide
traversal vs 0.196 s at matched m2/b1/z1). Matched-settings benchmark
(jobs 2871709-2871711/2871718-2871720, one verified assembly, identical first
1/100/1000 frozen-library guides, m2/b1/z1/f0.75/GG): new backend marginal
0.179 s/guide wall, ~20,100 guides/h/assembly on one core, all rows
independently audited; old preserved binary 11.1/37.5/383.2 s with 99/85/79%
of one CPU. Old emits ~575 rows/guide with only ~5% strict-eligible; **all
28,750 old strict-eligible witnesses are contained in the new output, zero
missing**. Old-engine full-library wall is measured at 5,775 s (~1.9 cores);
linear extrapolation from N=1000 overestimates it 5.15x, so new-engine
extrapolations (3.94 h/assembly matched) remain unvalidated at scale.
Requested m5/b1/z2 full Brunello would be ~1,900 core-hours per assembly
(~100 core-years over 464 assemblies) — infeasible without authorized
algorithmic/domain change. Batch per assembly, parallelize across assemblies;
full-library/full-cohort launches remain gated on user review.

## Completed validation and whole-assembly experiments

Independent bounded oracle panels recovered **204/204** and **172/172** anchors, with no unexpected or duplicate anchors. Tests include both strands, mismatch/bulge cases, reference Ns, PAMs, alternative alignments, and supported options. Additional bounded large-coordinate probes exercised HG00735_pat positions around `2^31`, contig boundaries and near the assembly end; these are targeted checks, not an exhaustive recall oracle.

All three jobs below used the imported **HG00735_pat** assembly (3,029,527,752 bases, 82 contigs), pilot guide `GCCAGGTACAGGTTCCGTCC`, and the settings above. The batch additionally used exactly the two guides listed in the next table.

| Slurm job | Scope | Reported anchors | Native search including locate/replay | Total runtime |
|---|---|---:|---:|---|
| 2871278 | First successful single-guide pilot | 21,447 | 82.465 s | Supervisor 1,070.97 s |
| 2871283 | Optimized verification/reference handling, same search/scoring policy | 21,447 | 80.5854 s | Supervisor 645.6085467999801 s; Slurm 10 min 47 s |
| 2871284 | Completed three-guide batch, one snapshot/native process/index load | 78,582 | 257.493 s | Supervisor 837.874695419101 s; Slurm 13 min 59 s |

The optimized standalone run reduced observed total runtime by about 40% versus the first successful pilot. This is a historical same-workload observation on a shared system, not a general speedup guarantee. The two standalone PAFs are byte-identical. All 21,447 records passed independent eligibility and CFD audits in each run, with the known exact target recovered.

### Three-guide batch: job 2871284

| Guide ID | Exact guide sequence | Reported anchors / CFD validated | Search states | Native search s |
|---|---|---:|---:|---:|
| pilot | `GCCAGGTACAGGTTCCGTCC` | 21,447 / 21,447 | 856,229,607 | 80.4914 |
| beginning | `CATTATTCTCGGCTGCAAAG` | 28,835 / 28,835 | 926,942,591 | 87.5259 |
| before_2_31 | `TGAACTTCTGATAACTCCCC` | 28,300 / 28,300 | 931,599,680 | 89.4761 |

**All 78,582 reported sites passed independent eligibility and CFD audits. Genome-wide recall remains unproven.** The audit extracted reference targets/PAMs, replayed CIGARs, checked coordinates, strands, mismatch/bulge/match-fraction limits and independently recalculated CFD values (all were supported and matched at output precision). It found zero duplicate anchors and recovered all three independently known exact targets. Counts include those three targets; they are not counts excluding intended on-target sites.

The pilot subset matched job 2871283: 21,447 shared anchors, zero missing/unexpected anchors, zero representative-alignment or CFD differences, and byte-identical PAF records. This demonstrates output stability and within-batch reuse, not an independent genome-wide sensitivity estimate.

- Native index loading: **3.30994 s**, exactly one logged event; snapshot copy: **97.069908 s**, one event.
- CLI workflow: **654 s**, whole-second Bash timer, distinct from native search and total job time.
- Peak sampled summed process-tree RSS: **14,267,208 KiB (13.61 GiB)**.
- Ran on **octopus10**, one CPU, requested 64 GiB, two-hour limit, no requeue/restarts, sampled 48 GiB stop threshold at 100 ms intervals. **Slurm did not establish a hard 64 GiB memory limit.** Sampling can miss spikes and double-count shared pages.
- The other two guides were not run separately; there is **no measured whole-assembly three-single-job versus batch speedup**.

Authoritative batch report and small evidence: [batch_2871284/REPORT.md](../research/pam_anchored_poc/optimized_pilot/batch_2871284/REPORT.md), [audit.json](../research/pam_anchored_poc/optimized_pilot/batch_2871284/audit.json), [per-guide-metrics.json](../research/pam_anchored_poc/optimized_pilot/batch_2871284/per-guide-metrics.json), and [pilot-comparison.json](../research/pam_anchored_poc/optimized_pilot/batch_2871284/pilot-comparison.json).

## Reuse, artifacts and reproducibility

Do not rebuild existing HPRC indexes unnecessarily. Original indexes are under `/moosefs/salehi/HPRCv2_CRISPRapido/indices/`. The verified imported bundle is:

```
/export/local/home/salehi/crisprapido-hprc-import-hg00735-20261105-v1/
```

It contains `reference.fa`, the existing Columba index components and `index.crispr-binding.json`; it was imported **without rebuilding the FM-index**. Preserve mandatory reference/index verification, original N-containing reference bytes and fail-closed behavior.

Large outputs remain outside Git under `/export/local/home/salehi/crisprapido-resource-pilot-results/`:

- `hg00735-pilot-2871278/`
- `hg00735-optimized-2871283/`
- `hg00735-batch-2871284/`
- Corresponding `monitor-2871278/`, `monitor-optimized-2871283/`, `monitor-batch-2871284/` directories retain sampled resource logs/summaries.

Expected PAF SHA256:

- Both standalone PAFs and the batch pilot subset: `a28d0574a915a57d431851ab681b0602b8f5fb1740d1c8efeb6e3b4a378f23be`.
- Full batch PAF: `a038919786424e8312eebe5eedb64f5c9dd8bee2074615823da577d7987d36a4`.

The isolated CLI was built from `f7de425` at `/export/local/home/salehi/crisprapido-pilot-rebuild-f7de425-v1/target/release/crisprapido`. Measurement-only runner/native instrumentation is separately documented and checksum-pinned; production binaries were not modified. Preserve its archived source because the executable embeds the runner path.

Batch reproducibility files live in `research/pam_anchored_poc/optimized_pilot/`: `batch_guides.fa`, `batch_targets.json`, `batch.env`, `submit_batch.slurm`, `batch_workload.sh`, `runtime.sha256`, `audit_batch.py`, and `extract_batch_metrics.py`. **Use `batch.env`, not just `optimized.env`, for the per-guide instrumented runtime.** Submission is not needed now and requires new approval plus fresh worker/path/capacity/integrity checks for any future run. The historical preparation [report](../research/pam_anchored_poc/optimized_pilot/REPORT.md) is retained with a current-status notice.

## Important commits

| Commit | Milestone |
|---|---|
| `8dcc2ad` | Optional integrated PAM-anchored backend; existing backend retained |
| `671d0b6` | Verified import of existing Columba indexes without rebuilding |
| `53a3b88`, `02bacec` | Bounded large-coordinate and repetitive HG00735 probes |
| `1d3ed09` | Portable timing used for successful monitored pilot |
| `f7de425` | Verification/reference-loading overhead reduction, unchanged search/scoring |
| `bfbb8ac` | Optimized comparison and audited batch preparation |
| `dfd1636` | Completed and audited three-guide HG00735 batch |
| `252a6bb` | CFD bound feasibility study, independent bounded oracles and model tests |
| `178de00` | Optional downstream CFD filtering and isolated native pruning, not enabled in the integrated backend |
| `4756527` | Limited-domain isolated CRISPR-Bulge postprocessor and official-model validation |
| `90d958a` | Evaluated saved 21,447-row HG00735 PAF; initial CLI deferral, subsequently superseded below |
| `6fa96d3` | Optional live/saved-file CLI scoring, rich N/A annotations and pinned reference-offset reuse |
| `ce52c15` | Complete bounded workflow validation, fresh/offline installation evidence and research reproduction setup |

## Remaining limitations and proposed next step

- Emitted-site audits establish eligibility and score correctness for reported alignments, **not absence of missed sites**. Agreement between implementations is not an independent recall oracle. Bounded oracle recovery does not prove arbitrary genome-wide completeness.
- The experimental backend supports one thread, `b=0/1`, `z=0/1/2`, resolved A/C/G/T guides and A/C/G/T/N reference FASTA. Unsupported CLI modes fail explicitly; the existing backend remains available. The native format is the tested 32-bit vanilla, sparse-SA-4 layout, not an unrestricted index API.
- Different valid representatives can yield different CFD scores; neither backend's representative is assumed biologically superior. CFD is not a calibrated cleavage probability. Missing/unsupported scores must not silently become zero.
- One batch amortizes setup, but every invocation still copies the complete bundle and performs repeated full checks. There is no immutable verified cross-job cache. Unsigned manifests/read-only file modes do not protect against malicious same-user mutate-and-restore races; trusted inputs/process isolation are required. Output buffering and large/repetitive batches remain resource risks.
- **Integrated CFD-guided pruning remains disabled.** An isolated native pruning prototype exists. CRISPR-Bulge is now an experimental opt-in Rust CLI scoring mode using the limited-domain batched postprocessor; see the current integration section below. Wider guide/assembly performance and extreme-repeat behavior remain to be studied; three guides do not establish general scalability.

### Completed bounded CFD feasibility investigation

See [CFD feasibility report](../research/pam_anchored_poc/cfd_feasibility/REPORT.md) and its scripts, tests and small retained evidence. No production/search/scoring code or CLI changed; no whole-genome job was run. A PAM-proximal column-product upper bound is safe for defined scores under validated finite `[0,1]` factors, but branches that may yield unavailable scores must survive. Pure-D projection retains the last 20 columns; mixed I/D for a 20-base guide is unavailable. Reference Ns require conservative fallback or separately validated N-aware interval handling.

The independent forward path oracle matched unchanged native representatives in **25 tiny-index configurations** (both strands, all 16 PAM suffixes, mismatch limits 5/8/20). **125 model threshold comparisons** had no lost qualifying/unavailable paths or final anchor/representative differences after unpruned refinement. **15,478** distinct oracle scoring inputs agreed with Rust (4,861 defined, 10,617 unavailable). Nine study tests passed. These are finite synthetic checks, not whole-genome recall evidence. The forward oracle and reverse model share replay/scoring helpers; the native/probe comparisons add separate implementation checks.

At threshold `.1`, resolved-panel model states fell **70.156% without bulges**, but only **30.046% at m5/b1/z2/f.75** because mixed-gap possibilities must survive. These are Python occurrence-set model counts, not native FM-index speedups. The conservative N-containing-reference fallback gave **zero state reduction** and would do so on HG00735_pat. Naïve pruning changed representatives/membership; any future implementation must preserve the current comparator, for example with candidate discovery followed by unpruned anchor-local refinement.

### Implemented downstream filtering and isolated native pruning

The user approved canonical-representative threshold semantics. Fresh development CLIs now accept optional **`--min-cfd SCORE`**, finite in `[0,1]`, disabled by default. The threshold is applied after the existing representative is selected, using the formatter's full-precision CFD value once. Defined scores below the threshold are removed; unavailable scores are kept without promoting another alignment. Omitted `cf:f` still denotes unavailable scoring, not zero. Original CLI options/defaults and higher mismatch limits remain unchanged; integrated searching remains unpruned with or without this flag.

See the [implementation/validation report](../research/pam_anchored_poc/cfd_pruning/REPORT.md), code and small retained evidence. The separate native executable performs real FM-index pruning and unpruned anchor-local reconstruction. It reuses unchanged native search/emission for controls, rounds bounds conservatively, preserves mixed-gap possibilities, and falls back on references containing any N, non-20-base guides, zero threshold or unsupported table contracts. Its distinct protocol is rejected by the integrated adapter; **do not set it as a production `--pam-search-bin`**.

All **72 A/B/C comparisons** (18 configurations × thresholds `0,.1,.2,.5`) matched integrated unpruned+filter, isolated pruned, and exhaustive direct-string oracle final anchors/representatives. Default output was byte-identical to the preserved prior CLI in all 18 cases. Additional native/oracle checks covered the other 15 PAM suffixes at four thresholds, custom-table/non-20 fallbacks, and binding rejection. Both strands, mismatch+bulge cases, two-base mixed events and reference Ns were tested. All 32 unavailable canonical mixed-gap representatives in the resolved m0/z2 panel survived every threshold. **164 Rust tests, 5 native/CLI regressions and 9 prior model regressions passed**; one pre-existing external-Python CFD test remained ignored.

Actual native measurements (three repetitions per mode; threshold `.1`, m5/b1/z2/f.75): CHM13/GRCh38 **50 kb** windows explored about 51% fewer states and had **22.015% / 20.972% lower native algorithm times**, including reconstruction/filtering/TSV construction but excluding loading, verification and full CLI setup. Dense synthetic output explored 31.454% fewer states but was **5.766% slower** after reconstruction. N fallback reduced no states. These are small native benchmarks, **not modeled predictions, whole-workflow speedups or genome-wide guarantees**. Native peaks were about 7.3 MiB; process-tree sampling has the usual shared-page/spike limitations. The report gives exact times, states, memory and scopes.

Runtime root: `/export/local/home/salehi/crisprapido-cfd-pruning-dev-v1/`; final native `native-build-v4/native`, development CLI `target/release/crisprapido`, authoritative logs `validation-final-v4/`. Production binaries and existing indexes remain untouched. Reproduction uses a new output directory via `research/pam_anchored_poc/cfd_pruning/reproduce.sh`; no large genomic input/output is in Git.

**Stop point:** do not enable integrated or HPRC pruning without separate approval and further independent review. Next candidates are broader bounded/randomized oracle coverage, custom-table/numerical review, matched full-workflow profiling and completeness-preserving N-aware discovery. Preserve representative policy and higher mismatch options. Threshold-selected output is not exhaustive under mismatch/bulge settings alone; successful finite tests do not establish genome-wide recall.

### Initial limited CRISPR-Bulge postprocessor milestone (historical)

**The statements in this initial milestone about no core CLI option are superseded by the current CLI integration below.** Original results/limits are retained.

See [CRISPR-Bulge report](../research/crispr_bulge/REPORT.md), scripts, manifest and retained evidence. Official source is `OrensteinLab/CRISPR-Bulge`, pinned **`3eddcd5bfcaff00b2bdf29425116ec9756ace870`**; the selected upstream-recommended model is the five-component GRU-Embd (`c_2`) classification ensemble. Root license is MIT; attribution is retained, while weights, pickles, datasets and the dedicated dependency environment remain outside Git. Inference uses checked HDF5 models without unpickling; only original-reference validation uses the trusted pinned upstream pickle loader.

**Full compatibility with `m5/b1/z2` is not achieved.** Checkpoints accept 24 columns including PAM; a two-base DNA bulge requires 25 and the official encoder rejects it. Published GUIDE-seq preprocessing uses six mismatches without gaps, four with gaps, and one gapped nucleotide (not one arbitrary-length event). The adapter scores only resolved 20-base-guide NGG alignments within that conservative subset. Two-base RNA, mixed/multiple gaps, five-mismatch single-bulge cases, ambiguous reference/guide/PAM, non-NGG and other unsupported cases remain present with an unavailable reason. Some are numerically encodable but scientifically unvalidated, not hard model-shape failures. No cropping, CFD-style projection, realignment, retraining or other workaround was implemented; limitations were reported before implementation.

The isolated `research/crispr_bulge/score_paf.py` provides **`--scoring cfd|crispr-bulge|both`**, default **cfd**. CFD mode passes the existing PAF through byte-for-byte; `both` appends `cb:f` only when scored and `bu:Z` status, retaining existing CFD tags/rows/CIGARs/order. CB-only intentionally omits `cf` tags, although CFD was still calculated upstream. It never selects a different alignment or filters sites. Core Rust CLI/search/scorer code, defaults and optional `--min-cfd` remain unchanged; **there is no core Rust `--scoring` option yet**. Apply `--min-cfd` upstream if wanted, not to rounded PAF values or as a neural-score threshold.

Validation: **542 rows retained** (370 synthetic plus the existing 172-panel output), **385 scored / 157 unavailable**; retained-panel coverage was **43/172 scored**, with all 129 other rows retained. Default CFD output was byte-identical on both panels; all original fields/CFD tags in `both` were preserved. Independent reference/CIGAR extraction and exact categorical-feature parity were checked against unmodified official `ensemble_predict` on identical inputs. Maximum prediction difference was **0** at chunk 512, **2.0116567611694336e-7** at chunk 32 (385 records), within `1e-6`; actual annotated PAF scores matched within **4.990997060083657e-10**. Singleton and three repeated 1,024-row comparisons also passed. **10 stdlib-only regression tests passed**. This is implementation agreement, not biological accuracy or genome-wide recall evidence.

Measured median fresh-process CPU wall time for 1,024 supported repeated alignments: Rust CFD probe **0.008286481 s**, isolated CB worker **5.070746292 s**, official reference wrapper **6.322246979 s** (three repetitions; scopes include startup/loading/I/O, not matched kernels or full search). Maximum individual-child RSS was **7,052 / 562,532 / 650,904 KiB** respectively. Model import/loading dominate cold calls; five models load once per invocation and are reused across record chunks/guides. The report distinguishes exact worker peaks from sampled summed tree RSS and records numerical/batch variation.

Runtime paths: upstream/models `/export/local/home/salehi/crisprapido-crispr-bulge-review-v1/`; isolated Python 3.11.15/TF-CPU 2.12.0 environment `/export/local/home/salehi/crisprapido-bulge-env-v1/`; raw results and preserved small legacy inputs `/export/local/home/salehi/crisprapido-bulge-validation-v1/`. The 172-panel search was not rerun; no whole-genome scoring/query or Slurm job occurred.

**Stop point:** review limited coverage/availability semantics before any core integration or extension outside the documented domain. Complete two-base-DNA support needs an appropriate model/data strategy, not truncating alignments. Biological validation, larger bounded batches, licensing/packaging review and full-workflow profiling remain. Do not transfer CFD pruning bounds to this nonlinear scorer, claim either score is a calibrated cleavage probability, or promote pruning automatically. Manuscript updates and full HPRC analysis still require separate decisions.

### Evaluated CRISPR-Bulge on the saved whole-assembly output (historical milestone)

The deferral below records the earlier `90d958a` decision; the user subsequently approved the optional CLI work documented next.

See [saved HG00735_pat evaluation](../research/crispr_bulge/ASSEMBLY_EVALUATION.md), aggregate evidence and `evaluate_assembly.py`. Starting from `4756527`, applied the **unchanged** batched `both` postprocessor to the existing job **2871283** PAF, not a new search. All **21,447 original records and CFD scores** were preserved; input SHA256 remains `a28d0574a915a57d431851ab681b0602b8f5fb1740d1c8efeb6e3b4a378f23be`.

- **4,526 CB scores / both scores (21.10%)**, **16,921 unavailable**. Supported: 1,098 ungapped, 1,080 DNA1, 2,348 RNA1.
- Unavailable: **8,176 two-base DNA bulges** exceed 24 columns; **7,624 DNA1/five-mismatch** rows exceed conservative with-gap mismatch policy; **1,121 RNA2** rows exceed one-gapped-nucleotide policy. Zero malformed CIGAR/mixed/multiple-gap or ambiguity statuses. No unsupported prediction was invented.
- Independent samtools extraction and separate CIGAR reconstruction for **93 stratified/tail-enriched rows**; official feature hashes matched, official prediction maximum absolute difference **2.5598145825256324e-10**. **14 tests passed**. Implementation agreement, not biological accuracy or recall.
- Paired Spearman **0.4281**, Pearson **0.3447**; score distributions differ markedly, especially gapped versus ungapped. No direct calibration, common thresholds or interchangeability claim.
- One measured CPU annotation: **91.5757 s**, including **81.7841 s reference scan**, **6.0402 s neural worker**; peak individual RSS **561,608 KiB**, sampled summed tree **617,108 KiB**. Not search runtime or a repeated median; sampling caveats remain.
- External raw artifacts: `/export/local/home/salehi/crisprapido-bulge-hg00735-2871283-v1/`, annotated output `both/output.paf`, authoritative audit `evaluation-final.json`, independent check `independent-final/`. Models/reference/original results untouched; no large inputs/outputs committed.

**Readiness/stop point:** this is usable as an isolated research postprocessor and technically supports a future experimental opt-in wrapper, **not unrestricted production CLI promotion**. Core integration remains deferred pending agreement on partial availability, dependency/licensing packaging and CLI lifecycle/error tests. Full default-domain coverage remains unsupported. Findings are reported before any integration; no main CLI/search/default/binary changes, no new search or Slurm job, no main merge. Seek explicit user direction before proceeding.

### Current optional CLI scoring integration

See [CLI_INTEGRATION.md](../research/crispr_bulge/CLI_INTEGRATION.md), [small results](../research/crispr_bulge/results/cli-integration/summary.json), and `src/scoring.rs`. User approved partial-domain CLI integration starting from **`90d958a`**. Development CLI now supports **`--scoring cfd|crispr-bulge|both`**, default CFD unchanged. Neural modes spool original reporting once per invocation and call the existing five-model batched worker once, including multiple guides. They preserve all sites and existing CFD tags, append neural tags and write `annotations.tsv` with sequences/alignment/PAM/context plus explicit **N/A** values/reasons. CB-only intentionally retains existing CFD, like `both`; mode selection is recorded in metadata. No score fusion or implicit threshold is applied.

- **`--score-paf PATH`** annotates existing results without any search/CFD recalculation. Neural modes require a NEW `--scoring-output-dir`, absolute `--crispr-bulge-python` and existing `--crispr-bulge-model-root`; source-bundle/script requirements are documented. Original PAF is read only. Default saved-file CFD is byte pass-through with no model/reference scan. `--min-cfd` and search executable/index/candidate options are rejected with `--score-paf`; during search the existing threshold is still explicitly opt-in/disabled by default.
- Original default PAF bytes, both search algorithms, representative selection, CFD formulas/matrices, mismatch/bulge defaults, optional threshold semantics and integrated-unpruned policy are unchanged. Only reporting destination/lifecycle and a test fixture in verification changed. No production binary, existing index or original PAF was modified; only an external development target was built.
- Existing 21,447 HG00735 job 2871283 rows were checked in all three modes: all CFD tags/rows preserved, **4,526 CB / 16,921 N/A** in each neural mode. Same unavailable reasons: **8,176 DNA2 over-length, 7,624 DNA1/five-mismatch domain exclusions, 1,121 RNA2 unsupported**. Both neural PAFs equal the previously evaluated annotated PAF byte-for-byte.
- Independent samtools extraction/separate CIGAR replay plus **unmodified official ensemble inference for ALL 4,526 supported rows**, not just the prior sample. Feature hashes exact; max PAF-score difference **2.5598145825256324e-10** in each mode, within `1e-6`. This is implementation agreement, not biological accuracy/calibration/recall.
- Optional `--scoring-reference-index` reuses explicit trusted offset manifests prepared by `research/crispr_bulge/reference_index.py`. Requires **both trusted index and reference SHA256 pins**. Full reference hashes still occur before/after scoring and index digest is rechecked; no automatic unsigned cache trust, no weakening search/index verification. Fallback remains the full original Python FASTA scan. Plain regularly wrapped FASTA required in neural modes.
- One fresh-process saved-result measurement per mode on octopus01: **12.929960 s CB / 13.781680 s both**, worker ~6.1 s, initial pinned-index/hash check ~2.23 s; max individual RSS **573,596 / 577,588 KiB**, sampled summed tree **629,832 / 634,612 KiB**. One-time index preparation **75.4562 s** (separate, not included in annotation timing). Prior full-scan annotation was 91.576 s; this is amortized reference-validation improvement, not a search speedup or repeated median. CFD pass-through **0.007562 s** is not CFD calculation time.
- **164 Rust tests passed** (one pre-existing ignored), **20 configured Python tests passed**, including tiny live single/multi-guide reporting, prior default compatibility, all modes, unsupported DNA2/RNA2/mixed/N/A/zero, empty inputs, failure/timeout, model reuse and pinned-index identity checks.
- External root `/export/local/home/salehi/crisprapido-scoring-cli-v1/`: development binary `target/release/crisprapido`, trusted offsets `reference-index.json`, complete scoring bundles under `saved-results/{crispr-bulge,both}/scoring/`, full official checks `saved-results/official-all/`. No genomic/model outputs committed.

**Stop point/limitations:** optional experimental CLI scoring is implemented, not a production packaging release. Full default-domain model coverage remains unsupported; no fabricated scores, truncation, retraining or re-alignment. Source-bundle/TF dependency and licensing/packaging review remain; full pre/post reference hashing and disk spooling still cost resources. No hard CLI memory cap/global timeout, GPU/concurrent-writer protection, broad biological validation or genome-wide recall guarantee. Do not transfer CFD pruning bounds, turn on thresholds/pruning automatically, change defaults, launch a new whole-genome search or merge main without separate direction.

### Complete bounded workflow and reproducible research setup

Current continuation starts from **`6fa96d3`**. See [end-to-end/reproducibility report](../research/reproducibility/REPORT.md), [installation guide](RESEARCH_INSTALLATION.md), [workflow results](../research/reproducibility/results/workflow-summary.json) and [installation provenance](../research/reproducibility/results/installation.json). **No core algorithms/scoring/defaults/Cargo files changed.** Two stale compatibility tests now create proper bound-index fixtures instead of raw unbound indexes; native integrity enforcement was not weakened.

- Actual **PAM-anchored FM-index search → canonical reporting → CFD → batched neural scoring → final output**, not saved-PAF pass-through. All three scoring modes recovered **566/566 resolved**, **570/570 N-containing**, **80/80 mixed-only m0 control**, and **172/172 pinned retained-panel** expected anchors; **zero missing/unexpected/duplicate anchors**. Default CFD PAF was byte-identical to the preserved preintegration CLI in every panel. All CFD scores/availability were independently replayed against tables; all rows/CFD tags/coordinates/CIGAR/order survived neural modes.
- CB defined/unavailable by panel: **550/16**, **552/18**, **32/48**, **43/129**. N/A reasons explicitly retained in TSV (DNA2, RNA2, mixed, ambiguity and mismatch-domain exclusions). Default search limits remain m5/b1/z2/f.75/NGG, with explicit one-thread tests; only the mixed-only control explicitly requests m0. No score threshold/pruning flag is passed.
- **1,177 supported panel-row predictions** matched unmodified official ensemble inference within `1e-6`; maximum absolute difference **1.789291992349007e-7**, allowing known batch/float variation. Panels overlap intentionally; this is not a count of distinct biological sites. Finite oracle recovery/implementation agreement does not prove biological accuracy/calibration or arbitrary genome-wide recall.
- Identical **13,956-base / 606-contig N panel**, three fresh complete CLI calls per mode: median **0.214627 s CFD / 5.020807 s CB-requested / 4.970432 s both**. Maximum individual RSS **15,056 / 538,428 / 536,748 KiB**; sampled summed tree **51,008 / 589,128 / 587,724 KiB**. Includes search/snapshot verification/loading/reporting/scoring/publication, not one-time index construction or independent validation. CPU/shared-host/warm-cache/sampling caveats remain; tiny-run setup dominates, no genome-scale speedup claim. CB-requested mode still retains upstream CFD like both.
- Clean `git archive 6fa96d3` source and **empty Cargo target** built with Rust/Cargo **1.94.1**, coherent GCC/G++ **12.3.0**, fresh WFA2 v2.3 library at **df3dce4**, fresh pinned Columba **edcc7ce/2.0.3** libraries/builder/unpruned native. Trusted Cargo/dependency source caches reused, not target objects. Guix/system header/linker/-lrt pitfalls resolved and documented, not a portable container certification.
- New Python **3.11.15** inference venv (TF-CPU **2.12.0**, NumPy **1.23.5**, Keras **2.12.0**, h5py **3.16.0**) passed full execution. A second new venv installed the **46-wheel hash lock offline**, had identical frozen packages, and produced a byte-identical complete retained172 PAF. New lock: `research/reproducibility/requirements-inference-linux-py311.lock` (Linux x86_64/CPython 3.11 only). Standard-library Python is required for PAM backend integrity/snapshot helpers even in CFD mode; no neural environment is needed for CFD.
- **Installation blocker:** actual network model-fetch helper attempts failed (Git 2.41 internal partial-fetch abort; Git 2.20 HTTP 400/RPC failure). TLS was not disabled. Execution used the five existing reviewed weights copied into a new external directory and SHA256-verified, **not a successful fresh network model install**. The official-reference environment was reused. Reliable model distribution, rights/security review and portable packaging remain release gates.
- **164 Rust tests passed**, one pre-existing ignored; **20 scoring + 4 new validator + 9 prior offline CFD model tests passed**. Also **96 PAM/option oracle comparisons**, five invalid-option/two invalid-guide checks, binding/tamper/PATH-spoof checks and both mixed-gap strands passed. The final bounded reproduction wrapper itself completed successfully.
- External root **`/export/local/home/salehi/crisprapido-research-workflow-v1/`**. Authoritative all-pass run **`final/`**, development binary **`target-final/release/crisprapido`**, clean source **`source/`**, fresh native **`pam-search`**, builder **`columba-build-v2/columba_build`**, model root **`models-verified/`**, venvs **`inference-env/`** and **`offline-env/`**, offline check **`offline-smoke/`**. Retain sources because helpers are path-embedded. Earlier exploratory/failed installation/test attempts are preserved but are not the authoritative run. No genomic/model/wheel/environment outputs committed.

**Release stop point:** bounded research workflow passes; no production release, main merge, default change, new whole-genome experiment or CFD-pruning enablement is authorized. Resolve reliable model setup/rights, portable build/relocation/environment packaging, partial-domain scientific policy, broader independent recall/biological benchmarking and resource/concurrent-writer hardening before release consideration. Large live-search/neural memory overlap and repeat-heavy outputs remain unprofiled; hashes assume trusted preparation and cannot defeat malicious mutate-and-restore races. Existing default Columba/WFA2 completeness limitations remain unchanged.

### Systematic search-mode concordance and manuscript scaffold

Current continuation starts from **`ce52c15`**. See [search concordance](../research/search_concordance/REPORT.md), [implementation taxonomy](../research/search_concordance/IMPLEMENTATIONS.md), [complete cases](../research/search_concordance/results/cases.json), [all raw-record audits](../research/search_concordance/results/records.tsv) and [manuscript outline](../manuscript/OUTLINE.md). **No algorithms/defaults/production binaries changed; no large experiment or main merge.**

- The five requested names are **not five independent engines**. “Original” can mean historical standalone window/WFA2 (`46cf652`), the window branch in older `f75f82aa` figure benchmarks, or Sassy-labelled original on divergent main (`cfb1587`). Current window/WFA2, automatic Columba+WFA2 and PAM-aware Columba are three complete configurations of the same attested release executable (built 6fa96d3, core identical to ce52c15). WFA2 is an aligner/verifier shared by two families; native PAM-aware traversal performs no WFA2 query verification. Sassy uses a different branch/interface; no merge was attempted.
- New **700-base / 28-contig adversarial controls**, exact guide/ref/PAM/strand and equivalent modern limits: at m3/b1/z1/f.75/NGG, **window 20/26, Columba+WFA2 26/26, native 26/26**, shared by all current modes **20**. At m3/b1/z2: **24/30, 28/30, 30/30**, shared **24**. Separate 46-base m0/z2 mixed-only control: **0/2, 0/2, 2/2**, shared **0**. Zero unexpected eligible/duplicate standard anchors in current pipelines; every emitted current alignment/defined CFD independently valid. Exact shared keys are committed.
- **All 70 anchor cases / 254 first-replicate raw records** audited and linked; repeats also checked. Current window misses are localized by preferred-alignment probes: terminal-D/PAM boundary, invalid preferred gap versus valid constrained alternative, and one optimum hiding another PAM in a dense window. Four panel-anchor Columba+WFA2 omissions have valid oracle candidates explicitly injected and still fail: unfiltered anchored `4X` versus eligible mixed/three-X, or `1X` versus eligible mixed/m0. Raw candidates already contain the exact interval for both m0 mixed cases. This is verification/fallback constrained-alignment loss, not simply absent raw discovery. No fix was made under the no-algorithm-change scope.
- Shared anchor `D5_2+`, boundary22: two independently valid representatives `[1,22) 1X2=1D2X15=` vs `[3,22) 1I1=3X15=` give CFD **0.2700 vs 0.1736**. Existing separate 204-site audit reused: **54 valid differing pairs, 41 CFD differences, 108 independent/Rust scores agreed**. Site agreement is not score/alignment identity or biological preference.
- Historical standalone WFA2 and Sassy are **diagnostics, not equal Venn/ranking sets**. Historical eligible-witness counts **11/13/0**; Sassy **16/16/0**, under separate schema/validity labels. Sassy exports eleven mandatory-like fields (missing PAF block length) and fixed target spans for gaps; z2 radius is `m+b`, not `m+b*z`. No coordinate/CIGAR repair was hidden in normalization. Exact old source-to-build attestation unavailable; Sassy release rebuild blocked by missing pinned Git cache after using compatible Cargo. Do not call absent valid exports intrinsic engine recall failures. Historical/Sassy causal attribution remains TODO where reporting obscures it.
- Matched common controls, three fresh calls/current mode: medians **0.004152 s window / 0.016266 s Columba+WFA2 / 0.114388 s native**; max individual RSS **7,076 / 20,376 / 15,100 KiB**. Native is not fastest on this tiny input; setup/integrity work differs and window misses sites. Historical/debug-path timings retained but **no fair five-way speed rank**. Earlier 8 kb/64-guide native-vs-Columba timings are reused separately, not rerun/pooled into a universal speedup. No chromosome/genome search occurred.
- **Six new audit tests passed** (anchors/strand, labelled legacy schema, no geometry repair, mixed witness/N/A, valid injected SAM tags/orientation, duplicates); **4 existing auditor + 9 CFD-model + 20 scoring Python tests** also reran and passed. Prior 164-pass Rust evidence is reused because core is unchanged, not claimed as a fresh rebuild. Executable/source hashes checked before/after. External root **`/export/local/home/salehi/crisprapido-concordance-v1/`**, authoritative **`final-v4/`**, earlier development/failed-tag probes preserved. Existing current/historical/Sassy binaries were read only; no production build replacement.
- Manuscript proposals: [outline/abstract structure](../manuscript/OUTLINE.md), [six figure topics](../manuscript/FIGURES.md), [claim/evidence/TODO ledger](../manuscript/EVIDENCE.md), [journal fit](../manuscript/JOURNALS.md). Existing figure/legend materials reviewed and preserved. No full prose draft was located in inspected material paths. Old manuscript whole-genome/chr figures predate strict PAM fixes (including non-GG outputs) and need re-audit before current-primary reuse. No population-level HPRC, broad recall, neural biological accuracy or universal speed claim is invented.

**Next stop point:** choose the intended paper baseline; obtain source/profile-attested historical/Sassy comparators and separately review exporter/semantic reconciliation; expand held-out bounded validation and re-audit old raw figures before any new large experiment proposal. Current native finite-control recovery is not a general completeness proof. Manuscript figures/outline are proposals with TODOs. Journal assessment is conditional (Bioinformatics methods novelty/rigor, NAR GAB genomics application, Bioinformatics Advances research-software fit); exact current submission policies/citations/availability remain to verify. No release, main merge, algorithm/default change, pruning enablement or new chromosome/whole-genome run is authorized.

## Working rules for future sessions

Check branch, commit and status first; preserve existing experiments, frozen branches, deterministic selection, production binaries and configurable defaults. Never merge into main or submit jobs without approval. Add regression tests for code fixes and run relevant tests. Commit/push completed source and small documentation/evidence only; exclude credentials, FASTA, indexes and large outputs. The latest continuation adds bounded normalized search-mode concordance/forensics and a manuscript/evidence/figure/journal scaffold. It changes no algorithm/default/production binary, enables no pruning and runs no large experiment. Use the taxonomy, full case ledger and explicit comparability/TODO gates before claiming a five-engine comparison or planning larger work.
