# Pilot blocker review and reproducible rebuild

2026-11-06. **Do not submit the whole-assembly pilot yet.** Build provenance has been resolved for a newly staged development CLI and a read-only worker diagnostic has been prepared. Hard 64 GiB job-memory enforcement, worker access to the required paths, and scratch quota remain unverified. Neither the diagnostic nor the pilot was submitted, and no assembly-wide guide search was run.

## Memory and partition

Controller configuration still reports `SelectType=select/cons_tres` with `SelectTypeParameters=CR_CPU`, `TaskPlugin=(null)`, `ProctrackType=proctrack/cgroup`, `ConstrainRAMSpace=no`, `ConstrainSwapSpace=no`, and `JobAcctGatherType=(null)`. CPU-only selection does not establish memory scheduling, and process tracking is not proof of a hard memory limit. No other cluster-wide job-memory enforcement mechanism is visible in this configuration. An external or site-specific mechanism cannot be ruled out without administrator confirmation or an authorized worker diagnostic. `--mem=64G` by itself is insufficient evidence.

On octopus01, `/proc/self/cgroup` and mount information show a hybrid hierarchy: cgroup v1 owns the `memory` controller at `/sys/fs/cgroup/memory`, with a separate cgroup-v2 mount at `/sys/fs/cgroup/unified`. The observed root v1 memory limit `9223372036854771712` is effectively unlimited and describes the controller host, not a worker allocation. Actual worker hierarchy and job-specific limits remain unknown. The pilot script requires a readable cgroup-v2 `memory.max` of at most 64 GiB and refuses otherwise; it does not pretend to support v1. A local fail-closed test using a fictitious job ID exited before any snapshot or search with `pilot preflight failed: a readable job cgroup-v2 memory.max is required`. That test was not a Slurm job.

The `workers` partition was up and advertised nodes large enough for 64 GiB and a two-hour job, but it is **not verified suitable for this pilot**: the visible configuration lacks the necessary enforced job-level memory limit, and direct noninteractive access to octopus07 was denied. Availability, mount identity and configured resources must be reconfirmed at submission time.

Questions for the HPC administrator:

1. What mechanism enforces a requested `--mem=64G` for the entire batch job and descendants, including Python, the native process and verifier children? Is it a hard limit, and how are OOM events and status reported?
2. Which partition/node type has enforced limits and supports a two-hour, one-CPU job? Does enforcement use v1, v2 or another mechanism, and what is the actual job-specific controller path? If v1, the pilot needs a separately implemented and tested fail-closed v1 check; do not remove its existing guard.
3. Are the staged `/export/local/home/salehi/...` bundle, repository, executables, `/gnu/store` runtime libraries and `/moosefs` FAI mounted with the **same content** on the chosen worker? May a non-search diagnostic allocation inspect them?
4. What worker scratch filesystem, user/project quota and output location are appropriate? The pilot needs at least 60 GB initially free in addition to the 14.54 GB persistent bundle, plus an enforceable aggregate quota; `ulimit -f` caps individual files only.
5. Are accounting, timeouts, OOM termination and cleanup guaranteed for the job and all child processes, and can a writable output directory remain available for the incomplete marker and logs?

## Lightweight diagnostic, not submitted

[worker_diagnostic.slurm](worker_diagnostic.slurm) requests one CPU, 1 GiB and ten minutes. It reports allocation information, cgroup v1/v2 memberships and observed limits, executable and helper accessibility, candidate executable hashes and usage/help, reference/index metadata and manifest/FAI hashes, filesystem mounts/free space and quotas where available. It does not search a guide, rebuild an index, run the whole-reference binding verifier, or write the staged bundle. Missing resources produce a nonzero diagnostic result. `ldd` and usage checks are restricted to the trusted candidate files. Its output itself needs a shared Slurm submission/output directory; this is among the points that must be confirmed.

**Only after explicit user approval** and path/partition confirmation, the prepared diagnostic command is:

```bash
cd /export/local/home/salehi/crisprapido-pam-validation
source research/pam_anchored_poc/pilot_preparation/pilot.env
sbatch --partition=workers --export=ALL research/pam_anchored_poc/pilot_preparation/worker_diagnostic.slurm
```

A numeric cgroup value reported by the diagnostic is not, by itself, proof of job-level enforcement: interpret it with the controller configuration and administrator. This command was not run.

## Build provenance

The prior development binary `3c1033355ecc896243197cbdb544ae5b7a5264b1c424ee81d2b78d43af64dfab` was previously tested but its precise source/build provenance could not be established retrospectively. A new isolated CLI build replaces it in [pilot.env](pilot.env):

- CRISPRapido source commit `923e8d7a08e3005b16d18af8ff63d7a83bd7324d`, tree `9c776c5ee89ef910ea6b9e5adb0e8d1a9f9af058`; tracked source and lock were unchanged during the build.
- `Cargo.lock` SHA256 `93564e52dcafb54bfdacd2b7c6d69ac806aa74fb00ea1a20a6aee75f13c7f0f2`; Rust `lib_wfa2` git dependency locked to `c608c436a6753d2c21c97d9f5c338efae99d042b` and verified in the offline Cargo cache.
- Repository-pinned WFA2 submodule `df3dce4b99c37ac4f34fdcef74b774c80b7fadc4`; source retrieved by authenticated GitHub SSH after HTTPS certificate verification failed, with no TLS verification bypass. It was built separately using `make -j2 build BUILD_TOOLS=0 BUILD_EXAMPLES=0 CC=gcc CXX=g++`; the WFA archive in that build has SHA256 `08c636ffb3496e1784520bc69a3b498bcbd88032c6bf7abb1412d4373455ac14`.
- Guix profile `/gnu/store/3c2mxm30rbzvnw7qsi235mrkk3m38fym-profile`: rustc 1.94.1 (`e408947bf`, LLVM 21.1.8), cargo 1.94.1, GCC 12.3.0. The linker was configured with the `/gnu/store/ln6hxqjvz6m9gdd9s97pivlqck7hzs99-glibc-2.35` interpreter and that glibc plus `/gnu/store/gmv6n5vy5qcsn71pkapg2hnknyn1p7g3-gcc-12.3.0-lib/lib` in the runtime path.
- Build command from the clean source checkout, after building WFA2: source the profile, set `LIBRARY_PATH` to the separate WFA library plus the profile value, set `CARGO_TARGET_DIR=/export/local/home/salehi/crisprapido-pilot-rebuild-923e8d7-v4`, set `RUSTFLAGS` to the exact linker arguments in [rebuild_pinned.sh](rebuild_pinned.sh), then run `cargo build --release --locked --offline`.
- Result `/export/local/home/salehi/crisprapido-pilot-rebuild-923e8d7-v4/release/crisprapido`, SHA256 `18b66c801be6f5756932145f435b87cd23a0eea6d0465badeb3f190834354fdf`. The timed native search wrapper and payload remain separately pinned in `pilot.env` and unchanged by this rebuild.

The rebuilt CLI plus the timed native executable recovered **172/172** independent 172-panel anchors, with zero extras and duplicates; the independent auditor reproduced 97 supported CFD scores. Its PAF was byte-identical to the prior candidate on that panel. This validates a bounded test, not worker compatibility or genome-wide completeness.

The separate [rebuild_pinned.sh](rebuild_pinned.sh) archives both exact commits into a new output directory, builds WFA2 there, and runs the locked offline Cargo build using the same profile and linker flags. It requires the pinned WFA2 commit to be available in a local Git repository and the locked Cargo dependencies to be cached. It succeeded on octopus01 with a second binary SHA256 `2e27988a3b6aeaee7e0b386e85cd2f6a83123c7ce138e8ed725a8678e1b52015`; that binary also passed the 172-anchor independent panel and produced the same PAF. The two builds used different source/build locations and have different checksums; this is a repeatable source-pinned build procedure, **not** a claim of bit-for-bit reproducibility. The isolated executable embeds its archived source directory as the helper path; preserve that directory and update and verify all associated helper pins before using that alternative in a pilot. The currently selected binary instead embeds the current checkout's runner path, which must remain trustworthy and available.

## Pilot command after blockers are resolved

The original one-guide, one-assembly pilot remains `-m 5 -b 1 -z 2 -f 0.75 -p GG`, one CPU, requested 64 GiB and two hours, with independent eligibility/CFD audit and a required known exact target. The memory guard remains mandatory and now gives a clear error on unsupported or insufficiently constrained cgroups. After administrator confirmation, any required script changes, diagnostic review, and separate user authorization, the current prepared command is:

```bash
cd /export/local/home/salehi/crisprapido-pam-validation
source research/pam_anchored_poc/pilot_preparation/pilot.env
sbatch --partition=workers --export=ALL research/pam_anchored_poc/pilot_preparation/one_guide.slurm
```

Do not execute this as written while the listed blockers persist; the partition, paths or cgroup implementation may need to change. An eventual resource pilot would not demonstrate genome-wide sensitivity.
