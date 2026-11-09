# Optimized HG00735_pat comparison and batch preparation

2026-11-09. The experimental development CLI was rebuilt in an isolated directory from commit **`f7de4251473ac5259e1902ca0cb36f3a7ba3ba2c`**; its SHA256 is **`9dac869a187856bff3c8309c3918ce48198d7f164b30f702c7c6f6fa33afd1ce`**. The previous executable and production binaries were not changed. The single-guide comparison and a three-guide batch are prepared for review, but **neither whole-assembly job was submitted**. Consequently there are no new whole-assembly timings, memory results, anchor comparisons, representative-alignment differences or CFD differences to report yet. This report does not claim the optimized whole workflow is faster.

## Build and instrumentation

[rebuild.sh](rebuild.sh) archives the exact source commit, builds separately against the repository-pinned WFA2 commit `df3dce4b99c37ac4f34fdcef74b774c80b7fadc4`, and runs locked offline Cargo in an isolated target. The resulting binary is `/export/local/home/salehi/crisprapido-pilot-rebuild-f7de425-v1/target/release/crisprapido`; its Cargo lock SHA256 is `93564e52dcafb54bfdacd2b7c6d69ac806aa74fb00ea1a20a6aee75f13c7f0f2`. It was built with the pinned Guix profile, rustc/cargo 1.94.1 and GCC 12.3.0. The archived source directory is retained because the Rust executable embeds its absolute runner path. The earlier CLI at `/export/local/home/salehi/crisprapido-pilot-rebuild-923e8d7-v4/release/crisprapido` remains intact.

For phase measurement, [prepare_profile_runner.py](prepare_profile_runner.py) first checks the archived Python runner against `f7de425`, saves an uninstrumented copy, and adds timers to the **separate archived runner**. Thus the CLI binary is compiled from the exact commit, while its runtime runner has explicitly documented measurement-only additions. [build_phase_native.py](build_phase_native.py) builds a separate native executable with timers around the existing mandatory binding check, reference load and index load; search and replay code are unchanged. The profiled runner, native executable and runtime helpers are pinned in [runtime.sha256](runtime.sha256) and [optimized.env](optimized.env). Runtime paths and libraries must exist on the chosen worker; ordinary same-user mutate-and-restore limitations remain.

The phase log records snapshot copying, pre/post snapshot verification, native process total, native binding check, native FASTA loading, index constructor, native search including alignment replay, Python snapshot-reference reconstruction, and runner total. [phase_metrics.py](phase_metrics.py) rejects missing or duplicate phase markers and reports the whole CLI time separately. The CLI timer uses Bash whole seconds, so the difference between CLI and runner time is a coarse combined residual that includes Rust parsing, scoring, formatting, output and startup; it is **not** a separately measured CFD or output phase. After the comparison, a separate benchmark calls the unchanged Rust CFD scorer on the emitted sites and an independent PAF profiler measures replay and output-writing proxy. Those post-run measurements do not claim to be the original in-workflow scorer/output time. The independent audit and comparison must succeed before output is published.

## Bounded verification

The rebuilt CLI, profiled runner and timed native executable were tested only on bounded prebuilt panels:

- The independent N panel recovered **172/172**, no extra or duplicate anchors; 97 supported CFD values matched. Its PAF was byte-identical to the retained pre-change experimental PAF.
- The independent resolved panel recovered **204/204**, no extras or duplicates; all 204 CFD values matched, with byte-identical PAF.
- A separate bounded three-guide fixture produced 13 eligible records, no duplicates and all three fixture exact targets; seven CFD values were supported and matched. Logs show **one** snapshot-copy event and **one** native index-load event. This demonstrates reuse in the small fixture, not a whole-assembly batch result.
- The multiguide auditor was tested to reject duplicate anchors, a missing exact target and incorrect mismatch counts. The single-guide comparator self-compared all 21,447 saved records and was tested to detect missing/duplicate anchors and a changed CFD tag. Supervisor limit and portable-timing tests passed.

Small evidence and phase logs are in [results](results/). The original saved job-2871278 PAF is not committed. Its expected SHA256 is `a28d0574a915a57d431851ab681b0602b8f5fb1740d1c8efeb6e3b4a378f23be`; the prepared single-guide workload refuses a mismatch. The comparison uses PAM-anchored `(guide, contig, strand, anchor)` keys, rejects missing/unexpected or duplicate sites and requires the known exact target. It reports representative-alignment and CFD differences rather than requiring identical CIGARs. A separate independent audit checks every emitted alignment and score; agreement with the baseline does not establish genome-wide sensitivity.

## Prepared one-guide comparison

[submit_single.slurm](submit_single.slurm) and [single_workload.sh](single_workload.sh) are prepared for HG00735_pat, guide `GCCAGGTACAGGTTCCGTCC`, `-m 5 -b 1 -z 2 -f 0.75 -p GG`, using the existing bound HPRC index. They request octopus10, one CPU, 64 GiB and two hours, disable requeue, and retain the 48 GiB sampled process-tree RSS stop from the prior authorized pilot. Sampling can miss spikes or double-count shared pages; Slurm's requested 64 GiB is not a verified hard memory limit. The script preserves full binding verification and the private snapshot; it also verifies runtime helper hashes and keeps results incomplete unless independent audit, baseline comparison and profiling succeed.

After separate approval, the prepared command is:

```bash
cd /export/local/home/salehi/crisprapido-pam-validation
source research/pam_anchored_poc/optimized_pilot/optimized.env
sbatch --parsable --export=ALL research/pam_anchored_poc/optimized_pilot/submit_single.slurm
```

This command was **not run**. Recheck worker availability, paths, capacity, integrity and the accepted memory-monitor limitation at submission time. The prior job's observed values were 758 s workflow, 82.465 s native search, 1,071 s total and 21,447 audited sites; its replay/reference loading/copy/internal-check phases were not separately measured, so not every new phase has an equivalent old measurement. The profiled runtime and post-run microbenchmarks introduce measurement overhead. Compare equivalent scopes and do not attribute a speedup to incomparable timing.

## Prepared three-guide batch

[batch_guides.fa](batch_guides.fa) contains exactly three guides: the pilot guide plus `CATTATTCTCGGCTGCAAAG` and `TGAACTTCTGATAACTCCCC`. Their independent known plus-strand exact targets and PAM anchors are recorded in [batch_targets.json](batch_targets.json), based on earlier source-FASTA and targeted-index probes. [submit_batch.slurm](submit_batch.slurm) and [batch_workload.sh](batch_workload.sh) use one `-G` invocation, one private snapshot, one native process/index load and the same verification policy. [audit_batch.py](audit_batch.py) independently checks per-guide PAF coordinates, strand, PAM, eligibility, CFD where supported, duplicates, and the three known exact targets. The supervisor records total runtime and sampled peak process-tree memory. The number of reported sites, runtime and peak memory are not yet known; no complete HPRC batch was launched.

After separate approval, the prepared batch command is:

```bash
cd /export/local/home/salehi/crisprapido-pam-validation
source research/pam_anchored_poc/optimized_pilot/optimized.env
sbatch --parsable --export=ALL research/pam_anchored_poc/optimized_pilot/submit_batch.slurm
```

One snapshot/batch is not an immutable reusable cache across jobs. The existing original backend and defaults remain available. No whole-assembly comparison, batch or full HPRC analysis was launched in this preparation.
