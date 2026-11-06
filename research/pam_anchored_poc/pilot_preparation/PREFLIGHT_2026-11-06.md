# HG00735_pat single-guide pilot: preflight, 2026-11-06

**Not ready for submission.** The controller reports `ConstrainRAMSpace = no`, `TaskPlugin = (null)` and `SelectTypeParameters = CR_CPU`. A `--mem=64G` request alone therefore has not established an enforced 64 GiB process-tree limit. This host uses a hybrid cgroup layout with its memory controller mounted as cgroup v1; the prepared script deliberately requires a verifiable cgroup-v2 `memory.max` at or below 64 GiB. It will refuse to run if that check cannot be met. Direct noninteractive SSH to `octopus07` failed with `Permission denied (publickey)`, so actual worker-node cgroups, path availability and scratch quotas could not be checked. No Slurm job was submitted.

Before approval, an administrator or documented cluster configuration must provide a verified hard job-level 64 GiB limit and appropriate cgroup visibility to the script, and confirm one CPU and the two-hour wall limit. A move to a suitable cgroup-v2 partition or a carefully implemented and tested v1 equivalent are possible, but neither is assumed available. Do not merely remove the fail-closed check. Confirm scheduler and filesystem behavior on the actual worker using an explicitly authorized preflight allocation if needed; no such allocation was started here. Resource limits in the script are not a promise of completion.

## Scheduler and storage findings

At inspection, default partition `workers` was up with seven 48-CPU nodes, unlimited maximum partition time, and idle nodes including octopus07–11. `octopus07` advertised 386,641 MiB configured memory and no allocated CPUs, although the controller's reported `FreeMem` was only 1,011 MiB and must not be treated as guaranteed available memory. Availability changes. Slurm 24.05.3 reports `proctrack/cgroup`, but process tracking does not itself prove a hard memory ceiling; the controller explicitly reports `ConstrainRAMSpace = no` and does not use memory as a select resource (`CR_CPU`). Its `KillWait` is 30 seconds, but signal and cleanup behavior under an actual allocation were not verified.

On **octopus01**, the imported bundle occupies approximately 14.54 GB and the selected scratch parent `/export/local/home/salehi/crisprapido-resource-pilot-scratch` is on local ext4 `/dev/md0p1`, with 11,203,735,044,096 bytes filesystem-wide free at inspection. That is sufficient aggregate filesystem space for the 14.54 GB snapshot and the script's 60 GB initial-free-space check, but is not proof of per-user quota, reserved capacity, worker visibility, or shared mount identity. `quota`, `mfsgetquota`, and `lfs` were unavailable. `/tmp` had 37,664,862,208 bytes free, below the 60 GB preflight threshold. The `ulimit` setting caps each output file at 10 GiB but does not cap total scratch use. Establish and monitor an enforceable scratch quota and verify all paths and store dependencies on the assigned node before submission; `/export` is a local block device on the controller, so the same path on another node cannot be presumed to contain the bundle.

Full binding verification on octopus01 succeeded for the original imported reference and all 12 index components. The reference SHA256 is `bae8df3c5cf8ca96d3316b59eaa63fa37b7d3b58551c8a0700f2c530dc1e4f8b`; pinned manifest and original FAI hashes are checked again by the job. The original importer establishes reference identity under its documented provenance assumptions, and a full binding check is not an independent proof of FM-index construction. The adapter makes its own snapshot and mandatory checks. Files and unsigned manifests are not protected from a malicious same-user mutate-and-restore race; use a trusted isolated environment.

## Selected guide and known exact target

Exactly one guide is fixed by the job script:

```
GCCAGGTACAGGTTCCGTCC
```

Independent `samtools faidx` retrieval from the original compressed source returned `GCCAGGTACAGGTTCCGTCCAGG` at `HG00735#1#JAHBCH020000014.1:43957005-43957027` (1-based inclusive). The target is on the plus strand; zero-based guide coordinates are `[43957004, 43957024)`, PAM `AGG` occupies `[43957024, 43957027)`, and the PAM-adjacent anchor is 43957024. Earlier bounded FM-index tests located the corresponding global start 2,147,483,653 and independently checked the same source window. The guide's complete assembly-wide site count and runtime remain unknown. The independent PAF auditor now requires this exact target with `20=` and zero mismatches/bulges; failure to report it prevents publication as complete. The auditor also checks all reported rows for strict PAM, strand, bounds, CIGAR, mismatches, gap limits, match fraction, duplicates and independently reproduced CFD where supported. This checks the correctness of emitted rows and one known exact target, not recall of all eligible loci.

## Candidate executables and verification

The separately staged candidate CLI and native wrapper/payload paths and SHA256 pins are recorded in [pilot.env](pilot.env). On octopus01 the native program's dynamic loader and libraries require paths under `/gnu/store`; its wrapper invokes the loader explicitly. The CLI embeds the absolute repository runner path and the native payload embeds the binding helper and `/usr/bin/python3.7`; these and the scorer, matrix files and Python imports must exist and work on the worker. Helper and matrix checksums are checked before and after execution using [dependencies.sha256](dependencies.sha256). These are integrity pins, not authentication against a same-user adversary.

The candidate CLI is the previously benchmarked development binary with SHA256 `3c1033355ecc896243197cbdb544ae5b7a5264b1c424ee81d2b78d43af64dfab`; its exact source/build provenance is not independently attested. A fresh offline build attempt using cargo 1.85.1 failed because host build scripts required newer glibc than the system loader provides. This unresolved reproducibility limitation must be accepted or repaired before a pilot. The separate native payload was built from the current prototype with only an index-constructor timer added; its post-build loader/RPATH adjustment requires the explicit wrapper. The original production binaries were not modified.

On octopus01 the candidate CLI and timed native wrapper ran the existing bounded independent 172-anchor panel with the same limits and produced **172/172**, zero extra or duplicate anchors; 97 CFD values were independently reproduced and N-containing cases without supported scores remained unscored. Native logs contained index loading time, 4,229,884 search states and search time. This is a small multi-guide fixture, not the proposed assembly-wide pilot; it tests candidate compatibility, not worker dependencies or whole-assembly sensitivity.

The [Slurm script](one_guide.slurm) records total script elapsed time, `/usr/bin/time -v` for the workflow, sampled process-tree RSS, cgroup current/peak and memory events when supported, disk free space, native index-loading time, native search time/states, and audited reported-site count. Process-tree RSS is a one-second sum and can miss brief peaks or double-count shared pages; cgroup memory includes charged cache and may cover more than the search processes. A missing native timing or known exact target makes the result incomplete. A hard kill can prevent final logs and cleanup, but the initial `INCOMPLETE` marker is retained until all checks pass. The script has one guide, no retry, and leaves partial PAF under `.incomplete`.

## Exact command for later approval

Only after the memory-enforcement, cgroup-version, worker-path/dependency, quota and provenance blockers above are resolved, the prepared command from the checkout is:

```bash
cd /export/local/home/salehi/crisprapido-pam-validation
source research/pam_anchored_poc/pilot_preparation/pilot.env
sbatch --partition=workers --export=ALL research/pam_anchored_poc/pilot_preparation/one_guide.slurm
```

The command is provided for review and was **not run**. If site remediation changes the partition or paths, update the pinned environment and recheck the exact command rather than submitting it as written. A successful resource pilot would not establish genome-wide sensitivity.

Validation completed for these preparation changes: `bash -n` for the Slurm script; positive and negative tests of known-target, duplicate, mismatch, CIGAR and PAM auditing; metric extraction and missing-timer rejection; one-second process-tree sampling; recheck of all 17 recorded bounded exact cases; and the independent 172-anchor candidate-executable panel described above.
