# HG00735_pat monitored single-guide attempt, 2026-11-09

Slurm job **2871276** failed before the genome search started. It ran on octopus10, completed the full reference/index binding precheck, and exited because `/usr/bin/time` was absent on that node (`workload.sh` line 47, exit 127). The monitor reported `child_exit_127`; Slurm marked the job `FAILED` with exit `1:0`. No retry or additional guide was launched. The empty `output.paf.incomplete` is **not a result**, and the job status remains `INCOMPLETE exit=127`.

The requested guide was `GCCAGGTACAGGTTCCGTCC` with `-m 5 -b 1 -z 2 -f 0.75 -p GG`. The job requested one CPU, 64 GiB and two hours on octopus10; Slurm showed `Requeue=0`, `Restarts=0`, and `ReqTRES=cpu=1,mem=64G` but `AllocTRES=cpu=1,node=1`. No hard 64 GiB memory enforcement is claimed. The original fail-closed script was not modified. This separate attempt used a 100 ms sampled aggregate process-tree RSS supervisor with a 48 GiB stop threshold, no automatic retry, and process-tree termination on error or a sampled threshold breach. Summed RSS can double-count shared pages and sampling can miss brief peaks: monitoring is **not equivalent to a hard memory limit**.

| Observation | Result |
|---|---|
| Slurm total runtime | 2 min 57 s |
| Supervisor elapsed | 174.60 s |
| Sampled peak aggregate tree RSS | 30,116 KiB (29.4 MiB); precheck only |
| 48 GiB monitoring threshold | Not reached |
| Reference/index precheck | Passed; reference SHA256 `bae8df3c5cf8ca96d3316b59eaa63fa37b7d3b58551c8a0700f2c530dc1e4f8b`, all 12 component hashes verified |
| Indexed search / alignment / CFD / output audit | Not run |
| Reported off-target sites | Unknown; search did not run |
| Known exact target recovery | Not assessed |
| PAM, mismatch, bulge and CFD verification | Not assessed for this attempt |

The worker logged the exact preflight error: `research/pam_anchored_poc/monitored_pilot/workload.sh: line 47: /usr/bin/time: No such file or directory`. The executable's compatibility on a small panel and other tests from earlier reports do not turn this into a completed whole-assembly test. The scratch directory was removed by the workload's error cleanup. Runtime and memory above do not estimate the full search.

A future separately authorized attempt must first replace the hard-coded `/usr/bin/time` dependency with a verified worker-available timing mechanism, test the complete launch path and monitor on the worker, and preserve equivalent index verification, independent output audit and the explicitly accepted memory-monitor limitation. The submitted workload is preserved in Git at commit `971859a`; the working workload was later updated for the separately documented timing fix. A new job would be a retry and was not submitted here.

[submit.slurm](submit.slurm) and [supervise.py](supervise.py) are the scripts used; the submitted version of `workload.sh` is in commit `971859a` and its current version includes the later timing fix; [test_supervisor.py](test_supervisor.py) verified a normal child, a failing child and a spawned child exceeding a small synthetic tree-RSS threshold before submission. [Results](results/) contain the small supervisor summary, sampled RSS log, precheck result and failure messages; no genomic data or index is committed. Full logs and the incomplete artifact remain in `/export/local/home/salehi/crisprapido-resource-pilot-results/` under job 2871276. The job requested 64 GiB but the scheduler did not report memory as an allocated TRES, consistent with the earlier configuration inspection.
