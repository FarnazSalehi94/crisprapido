# CRISPRapido project handoff

## Current development state

- Repository: https://github.com/FarnazSalehi94/crisprapido
- Checkout: `/export/local/home/salehi/crisprapido-pam-validation`
- Development branch: `feature/pam-anchored-search-poc`.
- Latest implementation/experiment commit at this documentation review: **`dfd1636`** (three-guide batch completed). The earlier handoff ending at `bfbb8ac` is superseded; **do not submit the three-guide batch as an unperformed experiment**.
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

## Remaining limitations and proposed next step

- Emitted-site audits establish eligibility and score correctness for reported alignments, **not absence of missed sites**. Agreement between implementations is not an independent recall oracle. Bounded oracle recovery does not prove arbitrary genome-wide completeness.
- The experimental backend supports one thread, `b=0/1`, `z=0/1/2`, resolved A/C/G/T guides and A/C/G/T/N reference FASTA. Unsupported CLI modes fail explicitly; the existing backend remains available. The native format is the tested 32-bit vanilla, sparse-SA-4 layout, not an unrestricted index API.
- Different valid representatives can yield different CFD scores; neither backend's representative is assumed biologically superior. CFD is not a calibrated cleavage probability. Missing/unsupported scores must not silently become zero.
- One batch amortizes setup, but every invocation still copies the complete bundle and performs repeated full checks. There is no immutable verified cross-job cache. Unsigned manifests/read-only file modes do not protect against malicious same-user mutate-and-restore races; trusted inputs/process isolation are required. Output buffering and large/repetitive batches remain resource risks.
- **CFD-guided indexed pruning and CRISPR-Bulge scoring are not implemented.** Wider guide/assembly performance and extreme-repeat behavior remain to be studied; three guides do not establish general scalability.

Recommended next scientific development, **subject to approval**, is a bounded scoring/pruning feasibility study before modifying the indexed search. Define whether a CFD threshold applies to the current deterministic representative or to any eligible alignment; derive a conservative upper bound for partial search branches, including gap/PAM effects, and fall back to no pruning where the bound or score is undefined. Validate against independent exhaustive small-reference oracles and the unchanged unpruned backend, including alternative alignments, both strands, Ns and bulges. Threshold-pruned output must be labeled as such, not advertised as exhaustive under mismatch/bulge settings alone. Measure search states and whole-workflow costs separately.

Plan CRISPR-Bulge as a separately optional downstream scorer first: pin the exact method/model/version, confirm licensing, input and bulge support, define unsupported cases and independently validate reference predictions before integration. Do not assume CFD bounds apply to CRISPR-Bulge or use score agreement as evidence of biological accuracy. Larger batches, manuscript updates and full HPRC analysis require separate decisions; do not start them automatically.

## Working rules for future sessions

Check branch, commit and status first; preserve existing experiments, frozen branches, deterministic selection, production binaries and configurable defaults. Never merge into main or submit jobs without approval. Add regression tests for code fixes and run relevant tests. Commit/push completed source and small documentation/evidence only; exclude credentials, FASTA, indexes and large outputs. This handoff update is documentation-only and did not rerun experiments or change algorithms.
