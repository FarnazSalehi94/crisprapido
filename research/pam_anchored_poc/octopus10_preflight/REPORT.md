# Octopus10 preparation for the HG00735_pat pilot

2026-11-06, 19:41–19:43 UTC. **Octopus10 is not yet verified safe for the pilot.** Slurm showed the node idle, but direct SSH failed with `Permission denied (publickey)`, so the imported index, candidate executables, scratch and actual job memory control could not be checked on that node. The controller reports memory enforcement disabled. No diagnostic or genome-search job was submitted, and no search was run.

## Node and resources

`scontrol show node octopus10` reported `IDLE`, 48 logical CPUs, zero Slurm-allocated CPUs, CPU load 0.01, configured `RealMemory=322130` MiB and `FreeMem=163446` MiB. `squeue -w octopus10` showed no jobs at inspection. `FreeMem` is a transient scheduler observation, not an independently measured `MemAvailable` or a reservation; unmanaged workloads and reclaimable cache are unknown. `TmpDisk=0` is the scheduler's configured value, not a measurement of the worker's scratch filesystem. The node is in `workers`, `allnodes`, `386mem` and `lowmem`. Availability can change before an allocation.

Controller configuration remains `SelectTypeParameters=CR_CPU`, `TaskPlugin=(null)`, `ProctrackType=proctrack/cgroup`, `ConstrainRAMSpace=no`, `ConstrainSwapSpace=no` and `JobAcctGatherType=(null)`. A Slurm `--mem=64G` request therefore **does not establish an enforced 64 GiB limit**. The administrator may know of an additional mechanism, but none is verified here. The controller uses a hybrid v1/v2 cgroup layout; its memory cgroup is not writable by the current user. These observations do not prove octopus10's worker hierarchy or user privileges. Do not assume a per-process `ulimit`, free-memory reading or `timeout` supplies a hard limit for all descendants.

The existing pilot refuses to run unless it can read a numeric cgroup-v2 `memory.max` at or below 64 GiB. That check remains mandatory. If an administrator provides another hard job-wide enforcement mechanism, it needs a separately implemented and tested fail-closed check before use; no bypass or v1 implementation was added here. Confirm time limit enforcement, child cleanup, OOM reporting and a scratch quota with the administrator as well.

## Reference, executables and scratch

On the controller, the staged HG00735_pat directory and chosen scratch directory are readable, the pinned candidate CLI and native-wrapper hashes match [pilot.env](../pilot_preparation/pilot.env), and CLI help/native usage execute. The imported bundle's reference and 12 index components passed full binding verification in the earlier [pilot preflight](../pilot_preparation/PREFLIGHT_2026-11-06.md). The rebuilt candidate passed the independent 172-anchor panel as recorded in the [build review](../pilot_preparation/HPC_RESOLUTION_2026-11-06.md). Neither the controller test nor that earlier bounded panel verifies execution or paths on octopus10.

The pre-existing diagnostic log for **octopus06**, inspected in the [tux comparison](../tux_assessment/REPORT.md), showed the candidate executable and bundle accessible there via NFS from `octopus01:/export`; it also reported an effectively unlimited memory cgroup. It is evidence about octopus06 only. It does not verify the octopus10 mount, library loader, Python helper path or quota.

At inspection, the controller's `/export` filesystem had 11,202,328,801,280 bytes free and contains the chosen scratch parent; this is filesystem-wide space, not a guaranteed quota or evidence that octopus10 sees the same mount. `/tmp` on octopus10 and its quota could not be measured. The pilot requires at least 60 GB initially free beyond the persistent approximately 14.54 GB bundle, an enforceable aggregate scratch/output quota, and enough space for a full snapshot. A 10 GiB per-file limit does not enforce an aggregate scratch limit.

## Prepared Slurm scripts

[pilot.slurm](pilot.slurm) is a node-specific wrapper with `--partition=workers`, `--nodelist=octopus10`, one task/CPU, `--mem=64G`, two-hour wall limit and `--no-requeue`. It checks Slurm's node name and delegates to the unchanged [checked pilot](../pilot_preparation/one_guide.slurm), retaining exactly one guide, reference/index binding checks, the cgroup-v2 guard, scratch precheck, incomplete-output handling, independent alignment/CFD audit and expected exact target. `--no-requeue` is supported by the installed `sbatch`. A successful pilot would measure resources and emitted-site eligibility, not genome-wide sensitivity.

[diagnostic.slurm](diagnostic.slurm) pins a read-only, one-CPU, 1 GiB, ten-minute diagnostic to octopus10 and delegates to the previously prepared worker diagnostic. It does not execute a guide search or build an index. The diagnostic itself must **not be submitted without explicit user approval**. Its observed numeric memory values require administrator interpretation; a diagnostic that reports its observations successfully does not certify enforcement. The output directory must be visible on the assigned worker.

Both wrappers passed `bash -n`. Local shell simulations, **not Slurm jobs**, confirmed that the pilot wrapper rejects a wrong node and that the delegated pilot rejects the controller's unavailable cgroup-v2 limit before snapshot or search. Independent audit and metric-extraction tests passed, and all helper checksum checks passed. These checks do not validate octopus10.

After explicit approval for the diagnostic and confirmation that the output location is appropriate, its exact prepared command is:

```bash
cd /export/local/home/salehi/crisprapido-pam-validation
source research/pam_anchored_poc/pilot_preparation/pilot.env
sbatch --export=ALL research/pam_anchored_poc/octopus10_preflight/diagnostic.slurm
```

After diagnostics, administrator confirmation of a hard 64 GiB process-tree limit and usable scratch/quota, and any required fail-closed implementation changes, the prepared **pilot** command for separate authorization is:

```bash
cd /export/local/home/salehi/crisprapido-pam-validation
source research/pam_anchored_poc/pilot_preparation/pilot.env
sbatch --export=ALL research/pam_anchored_poc/octopus10_preflight/pilot.slurm
```

Neither command was run. Under the currently observed configuration, the pilot should be expected to refuse to start, not to enforce memory merely because 64 GiB was requested. If configuration, partition or paths change, review and revalidate the scripts and environment before submission.
