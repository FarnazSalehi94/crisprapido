# Portable timing fix for monitored HG00735_pat pilot

2026-11-09. The monitored workload now calls [portable_timing.sh](portable_timing.sh), which uses Bash's `SECONDS` variable instead of `/usr/bin/time`. It writes `workflow.time` with whole elapsed wall seconds and the command's exit code on normal success or failure. The helper returns the exact command status, so `set -e` still stops the workflow on errors; a signal or forced kill can interrupt timing output and leaves the attempt incomplete. Whole-second resolution means a subsecond workload may report zero. The supervisor's higher-resolution overall elapsed time and sampled process-tree RSS remain separate.

The previous job **2871276** was executed on octopus10 by Bash and its existing cleanup code already used the same `SECONDS` variable: its worker log reported `elapsed_script_seconds=174 exit=127`. That confirms this Bash clock was available and advanced on octopus10 during the original job. The new helper itself has **not** been run on octopus10: noninteractive SSH still returns `Permission denied (publickey)`, and no diagnostic or pilot job was submitted for this fix. An explicitly approved diagnostic allocation or other authorized worker access is needed to test the exact new launch path there. The original log remains unchanged; the old workload version is available in Git at `971859a`.

On the controller, [test_timing.py](test_timing.py) checked a one-second successful command, exit 37 and an absent command (exit 127), including failure propagation under `set -e`. The existing supervisor tests passed normal exit, child failure and a process-tree memory-limit stop. A bounded 172-anchor panel was run through the candidate CLI using the new timing helper; the independent audit reported **172/172**, zero extras and duplicates, and 97 independently scored CFD values. Its PAF matched the previous bounded output. It completed in under the one-second timer resolution and recorded `elapsed_seconds=0`; native search and index-load timings remain in the workflow stderr. These controller tests are not an octopus10 execution or whole-assembly sensitivity evidence.

The separate monitored job retains the **48 GiB sampled process-tree RSS stop**, 100 ms sampling interval, one guide, octopus10, one CPU, requested 64 GiB, two-hour wall limit, no requeue, reference/index binding checks, 60 GB initial scratch check, post-verification, independent PAF/PAM/eligibility/CFD audit and the known exact target requirement. Its `timing.sha256` is checked before sourcing the helper and again before publishing output. The original fail-closed script was not modified. Sampled RSS can miss short spikes or double-count shared memory and is not a hard memory limit; Slurm's requested 64 GiB is not established as enforced. The failed original attempt remains incomplete.

No Slurm job was submitted for this fix. The exact command to submit **one** monitored pilot, after review and any required worker confirmation, is:

```bash
cd /export/local/home/salehi/crisprapido-pam-validation
source research/pam_anchored_poc/pilot_preparation/pilot.env
sbatch --parsable --export=ALL research/pam_anchored_poc/monitored_pilot/submit.slurm
```

The script pins the same guide `GCCAGGTACAGGTTCCGTCC`, HG00735_pat imported index, and `-m 5 -b 1 -z 2 -f 0.75 -p GG`. An eventual incomplete or failed job must not be retried automatically.
