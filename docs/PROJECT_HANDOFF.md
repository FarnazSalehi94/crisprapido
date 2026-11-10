# CRISPRapido project handoff

## Current development state

- Repository: https://github.com/FarnazSalehi94/crisprapido
- Checkout: `/export/local/home/salehi/crisprapido-pam-validation`
- Development branch: `feature/pam-anchored-search-poc`.
- Whole-assembly experiment milestone: **`dfd1636`** (three-guide batch completed); handoff documentation baseline: **`7c4e77f`**; CFD feasibility milestone: **`252a6bb`**. The subsequent downstream filtering/isolated native pruning implementation is linked below; pruning is **not enabled in the integrated backend**. The earlier handoff ending at `bfbb8ac` is superseded; **do not submit the three-guide batch as an unperformed experiment**.
- Experimental backend is integrated but is neither the default nor a production release. The original **Columba candidate search → WFA2 verification → PAM filtering → CFD → PAF** backend remains the default and is preserved. Earlier original-versus-Columba/WFA2 benchmarking is separate from the experimental backend results below; do not repeat it without a scientific reason.
- `--search-backend pam-anchored` uses Columba FM-index primitives, seeds the PAM constraint, tracks mismatches/bulges, handles both strands and reference Ns, and reports one deterministic representative per PAM anchor without WFA2. CFD scoring remains downstream.
- Keep CLI-configurable settings **`-m 5 -b 1 -z 2 -f 0.75 -p GG`**; `-p GG` means strict NGG. Do not change defaults or representative selection without discussion. Representatives minimize `(mismatches + inserted bases + deleted bases, mismatches, compressed CIGAR, start)`, not CFD.

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

## Remaining limitations and proposed next step

- Emitted-site audits establish eligibility and score correctness for reported alignments, **not absence of missed sites**. Agreement between implementations is not an independent recall oracle. Bounded oracle recovery does not prove arbitrary genome-wide completeness.
- The experimental backend supports one thread, `b=0/1`, `z=0/1/2`, resolved A/C/G/T guides and A/C/G/T/N reference FASTA. Unsupported CLI modes fail explicitly; the existing backend remains available. The native format is the tested 32-bit vanilla, sparse-SA-4 layout, not an unrestricted index API.
- Different valid representatives can yield different CFD scores; neither backend's representative is assumed biologically superior. CFD is not a calibrated cleavage probability. Missing/unsupported scores must not silently become zero.
- One batch amortizes setup, but every invocation still copies the complete bundle and performs repeated full checks. There is no immutable verified cross-job cache. Unsigned manifests/read-only file modes do not protect against malicious same-user mutate-and-restore races; trusted inputs/process isolation are required. Output buffering and large/repetitive batches remain resource risks.
- **Integrated CFD-guided pruning remains disabled.** An isolated native pruning prototype exists. CRISPR-Bulge is now available only as a limited-domain, isolated PAF postprocessor, not a Rust CLI scoring option. Wider guide/assembly performance and extreme-repeat behavior remain to be studied; three guides do not establish general scalability.

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

### Limited optional CRISPR-Bulge scoring integration

See [CRISPR-Bulge report](../research/crispr_bulge/REPORT.md), scripts, manifest and retained evidence. Official source is `OrensteinLab/CRISPR-Bulge`, pinned **`3eddcd5bfcaff00b2bdf29425116ec9756ace870`**; the selected upstream-recommended model is the five-component GRU-Embd (`c_2`) classification ensemble. Root license is MIT; attribution is retained, while weights, pickles, datasets and the dedicated dependency environment remain outside Git. Inference uses checked HDF5 models without unpickling; only original-reference validation uses the trusted pinned upstream pickle loader.

**Full compatibility with `m5/b1/z2` is not achieved.** Checkpoints accept 24 columns including PAM; a two-base DNA bulge requires 25 and the official encoder rejects it. Published GUIDE-seq preprocessing uses six mismatches without gaps, four with gaps, and one gapped nucleotide (not one arbitrary-length event). The adapter scores only resolved 20-base-guide NGG alignments within that conservative subset. Two-base RNA, mixed/multiple gaps, five-mismatch single-bulge cases, ambiguous reference/guide/PAM, non-NGG and other unsupported cases remain present with an unavailable reason. Some are numerically encodable but scientifically unvalidated, not hard model-shape failures. No cropping, CFD-style projection, realignment, retraining or other workaround was implemented; limitations were reported before implementation.

The isolated `research/crispr_bulge/score_paf.py` provides **`--scoring cfd|crispr-bulge|both`**, default **cfd**. CFD mode passes the existing PAF through byte-for-byte; `both` appends `cb:f` only when scored and `bu:Z` status, retaining existing CFD tags/rows/CIGARs/order. CB-only intentionally omits `cf` tags, although CFD was still calculated upstream. It never selects a different alignment or filters sites. Core Rust CLI/search/scorer code, defaults and optional `--min-cfd` remain unchanged; **there is no core Rust `--scoring` option yet**. Apply `--min-cfd` upstream if wanted, not to rounded PAF values or as a neural-score threshold.

Validation: **542 rows retained** (370 synthetic plus the existing 172-panel output), **385 scored / 157 unavailable**; retained-panel coverage was **43/172 scored**, with all 129 other rows retained. Default CFD output was byte-identical on both panels; all original fields/CFD tags in `both` were preserved. Independent reference/CIGAR extraction and exact categorical-feature parity were checked against unmodified official `ensemble_predict` on identical inputs. Maximum prediction difference was **0** at chunk 512, **2.0116567611694336e-7** at chunk 32 (385 records), within `1e-6`; actual annotated PAF scores matched within **4.990997060083657e-10**. Singleton and three repeated 1,024-row comparisons also passed. **10 stdlib-only regression tests passed**. This is implementation agreement, not biological accuracy or genome-wide recall evidence.

Measured median fresh-process CPU wall time for 1,024 supported repeated alignments: Rust CFD probe **0.008286481 s**, isolated CB worker **5.070746292 s**, official reference wrapper **6.322246979 s** (three repetitions; scopes include startup/loading/I/O, not matched kernels or full search). Maximum individual-child RSS was **7,052 / 562,532 / 650,904 KiB** respectively. Model import/loading dominate cold calls; five models load once per invocation and are reused across record chunks/guides. The report distinguishes exact worker peaks from sampled summed tree RSS and records numerical/batch variation.

Runtime paths: upstream/models `/export/local/home/salehi/crisprapido-crispr-bulge-review-v1/`; isolated Python 3.11.15/TF-CPU 2.12.0 environment `/export/local/home/salehi/crisprapido-bulge-env-v1/`; raw results and preserved small legacy inputs `/export/local/home/salehi/crisprapido-bulge-validation-v1/`. The 172-panel search was not rerun; no whole-genome scoring/query or Slurm job occurred.

**Stop point:** review limited coverage/availability semantics before any core integration or extension outside the documented domain. Complete two-base-DNA support needs an appropriate model/data strategy, not truncating alignments. Biological validation, larger bounded batches, licensing/packaging review and full-workflow profiling remain. Do not transfer CFD pruning bounds to this nonlinear scorer, claim either score is a calibrated cleavage probability, or promote pruning automatically. Manuscript updates and full HPRC analysis still require separate decisions.

## Working rules for future sessions

Check branch, commit and status first; preserve existing experiments, frozen branches, deterministic selection, production binaries and configurable defaults. Never merge into main or submit jobs without approval. Add regression tests for code fixes and run relevant tests. Commit/push completed source and small documentation/evidence only; exclude credentials, FASTA, indexes and large outputs. The latest continuation adds an isolated, limited-domain neural scoring postprocessor and independent model validation; it does not change either search backend, enable integrated pruning, rerun completed whole-assembly jobs or modify production binaries.
