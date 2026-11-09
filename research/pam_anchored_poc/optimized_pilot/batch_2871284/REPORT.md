# Three-guide HG00735_pat batch: job 2871284

2026-11-09. The single authorized three-guide batch **completed successfully** on octopus10: 78,582 PAM-anchored records passed independent alignment, PAM and CFD validation with zero duplicate anchors. The pilot guide's 21,447 records were **byte-identical** to job 2871283. All three independently known exact targets were recovered. No additional job, retry, guide or assembly was run. The independent audit checks reported-site eligibility, not genome-wide recall.

## Per-guide results

Counts in the table are **all reported anchors**, including the specified known exact target. If that specified target is treated as the intended on-target site, subtract one from each total; the resulting off-target counts are shown as a convention, not an independent determination of experimental intent. Each guide had exactly one reported `20=` exact alignment. Every reported record passed the independent PAM/coordinate/strand/CIGAR/mismatch/bulge audit, and every CFD score independently matched. All duplicate counts are zero.

| Guide ID | Guide sequence | All reported anchors | Off-targets if designated exact target excluded | CFD validated | Duplicates | Search states | Native search s |
|---|---|---:|---:|---:|---:|---:|---:|
| pilot | `GCCAGGTACAGGTTCCGTCC` | 21,447 | 21,446 | 21,447 | 0 | 856,229,607 | 80.4914 |
| beginning | `CATTATTCTCGGCTGCAAAG` | 28,835 | 28,834 | 28,835 | 0 | 926,942,591 | 87.5259 |
| before_2_31 | `TGAACTTCTGATAACTCCCC` | 28,300 | 28,299 | 28,300 | 0 | 931,599,680 | 89.4761 |
| **Total** | | **78,582** | **78,579** | **78,582** | **0** | **2,714,771,878** | **257.493** |

The native search values sum separately instrumented plus/minus traversal-and-replay timings; rounded per-strand values need not sum exactly to the logged aggregate. The additional measurement-only native executable logged guide names, strand, state counts and elapsed time around each unchanged search call. It did not change search transitions, representative selection or CFD. [per-guide-metrics.json](per-guide-metrics.json), [audit.json](audit.json) and [workflow.stderr](workflow.stderr) contain the detailed evidence.

The independent audit checked sequence and strict NGG PAM against the reference, strand and genomic bounds, replayed each reported CIGAR and enforced mismatch, bulge and match-fraction limits. It independently recalculated each supported CFD value, and all 78,582 were supported and matched. The three exact target anchors are recorded in [batch_targets.json](../batch_targets.json). Passing this audit does not prove that the search found every eligible site in the assembly.

The pilot subset's PAF was also compared with the retained job-2871283 PAF: **21,447 shared, zero missing or unexpected, zero duplicates, zero representative-alignment differences and zero CFD differences**. A separate `cmp` confirmed the pilot PAF records were byte-for-byte identical. Both pilot PAF SHA256 values are `a28d0574a915a57d431851ab681b0602b8f5fb1740d1c8efeb6e3b4a378f23be`; the full 78,582-record batch PAF SHA256 is `a038919786424e8312eebe5eedb64f5c9dd8bee2074615823da577d7987d36a4`. The PAF itself is not committed.

## Performance and resources

| Metric | Measured batch value | Scope |
|---|---:|---|
| Total Slurm job | 13 min 59 s | Scheduler reported; supervisor 837.875 s |
| CLI workflow | 654 s | Whole-second Bash timer; includes snapshot, internal checks, search, scoring and output |
| Native search | 257.493 s | Sum of three guides and both strands; includes locate/replay |
| Native index loading | 3.310 s | Exactly one logged native index-load event |
| Snapshot creation | 97.070 s | One logged snapshot-copy event |
| Internal binding verification | 77.633 + 78.215 + 77.420 = 233.268 s | Snapshot precheck, native check, snapshot postcheck |
| Outer binding checks | about 80 s each | Wrapper timestamps before/after the workflow, separate from the 654 s |
| Native / Python reference loading | 18.657 / 27.411 s | Separate phases within the workflow |
| Peak aggregate process-tree RSS | 14,267,208 KiB (13.61 GiB) | 100 ms sampled sum, 48 GiB stop threshold |
| Total reported sites | 78,582 | Includes three designated exact target anchors |

The run used one `-G` invocation, one private snapshot, one native process and one native index load. It retained all pre-existing full checks; the single snapshot was verified repeatedly, rather than verified only once. It requested one CPU, 64 GiB and two hours; `Requeue=0` and `Restarts=0`. Slurm reported `ReqTRES=cpu=1,mem=64G` and `AllocTRES=cpu=1,node=1`, so the requested 64 GiB is **not claimed to be a hard enforced limit**. The 48 GiB stop is based on sampled summed process-tree RSS, which can double-count shared pages and miss short spikes; the threshold was not approached.

## Comparison with separate execution

The only measured whole-assembly standalone comparison is job **2871283**, for `pilot` alone: supervisor elapsed **645.609 s**, native search **80.585 s**, 21,447 sites. The batch took **837.875 s** for all three guides and produced the identical pilot result. The other two guides were **not** run separately, so an actual three-single-job runtime or measured savings is unavailable.

As an explicitly hypothetical model, take the earlier single-guide elapsed time minus its 80.585 s search as a fixed per-run overhead, charge that overhead three times, and add this batch's measured 257.493 s total search. That yields about **1,952.6 s** for three separate runs versus **837.9 s** observed here, a nominal difference of **1,114.7 s** (18 min 35 s). This is **not a measured speedup**: it assumes identical non-search overhead, ignores different site counts, caching and contention, and the standalone job included post-run microbenchmarks not performed by the batch. It should not be used as a guaranteed savings estimate. A separate measured 8 kb/64-guide comparison in the [workflow profile](../../workflow_profile/REPORT.md) demonstrated snapshot amortization on a bounded panel, not a whole-assembly speed ratio.

## Reproducibility and artifacts

The job used the existing imported HG00735_pat bundle with `-m 5 -b 1 -z 2 -f 0.75 -p GG`. [guide.fa](guide.fa), [slurm.log](slurm.log), [slurm-job.txt](slurm-job.txt), [progress.log](progress.log), [rss.log](rss.log), [summary.json](summary.json), phase timings, verification manifests and audit/comparison JSON files are small retained artifacts. Large genomic files, indexes and the 12 MB batch PAF are not committed. Runtime helper hashes were checked before and after the run; [batch.env](../batch.env), [batch_workload.sh](../batch_workload.sh), [submit_batch.slurm](../submit_batch.slurm) and [build_batch_native.py](../build_batch_native.py) record the job's instrumentation and setup. The batch-specific audit and per-guide metrics parser were tested on a bounded three-guide panel before submission. No production binary, original backend, default settings, search algorithm or CFD calculation was changed.
