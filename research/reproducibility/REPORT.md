# Complete bounded workflow and research reproducibility

## Status and scope

**The complete explicit PAM-anchored workflow passes the bounded checks below. It is usable for reproducible research with the documented trusted dependencies/models, but is not ready for an unrestricted production release.** A clean network model installation remains blocked, and portability, packaging, coverage and broader scientific/resource validation remain open.

Baseline search/scoring code: **`6fa96d3`**, branch `feature/pam-anchored-search-poc`. This continuation changes **no core search/scoring code, defaults, CFD formulas, representative selection, Cargo files or production binaries**. It adds a bounded reproduction driver, a hash-locked inference wheel set, tests and documentation. Two older compatibility tests were repaired to build proper bound fixtures; native integrity checks were not bypassed. No main merge, Slurm submission, new/repeated whole-genome experiment or default pruning activation occurred.

The tested route is:

**verified bound index → PAM-seeded Columba FM-index traversal → canonical alignment/PAF reporting → existing CFD calculation → optional batched CRISPR-Bulge → final PAF/TSV**.

This is **live search**, not `--score-paf` pass-through or an isolated scoring benchmark. The route explicitly selects `--search-backend pam-anchored`; the original Columba/WFA2 backend remains default and unchanged. Its previously documented candidate/verification limitations are not resolved or recharacterized by this run. Genome-wide completeness of either route remains unproven.

## 1. End-to-end results

Existing fixtures from `cfd_feasibility/study.py` were regenerated without new biological input. For resolved/N/mixed panels, a separate direct-string forward alignment oracle enumerated eligible PAM anchors and selected the existing canonical representative. The retained 172 panel used its checksum-pinned, previously independently validated inputs/PAF. All runs use default `m5/b1/z2/f.75/NGG`, except the **explicit `m0` mixed-only control**; `-t 1` is explicit. No score threshold is passed.

| Existing bounded panel | Reference bases / contigs | Expected / recovered | CFD defined / unavailable | CB defined / unavailable |
|---|---|---|---|---|
| Resolved fixtures | 13,768 / 598 | **566 / 566** | **566 / 0** | **550 / 16** |
| Same fixtures with N controls | 13,956 / 606 | **570 / 570** | **568 / 2** | **552 / 18** |
| Mixed-only control, same resolved reference, m0 | 13,768 / 598 | **80 / 80** | **48 / 32** | **32 / 48** |
| Retained independent 172 panel, 64 guides | 8,000 / 8 | **172 / 172** | **97 / 75** | **43 / 129** |

Each panel ran all three live CLI modes: `cfd`, `crispr-bulge`, `both`. **Zero missing, unexpected or duplicate anchors** occurred. Coordinates, canonical CIGAR and mismatch counts matched the independent expected representatives. Independent bounded sequence/CIGAR replay and matrix evaluation checked every CFD availability/value at PAF precision (absolute tolerance `0.000051`). Default/explicit CFD output matched the preserved preintegration CFD CLI byte-for-byte for every panel.

All original fields, CFD tags (including absence), records and ordering survived both neural modes and repeated runs. Within each panel the neural PAFs were byte-identical across modes/repetitions, so the official-reference comparison applies to each output. Each sidecar has one row per hit; guide/target/aligned sequences, observed PAM and CFD values were independently reconstructed/checked. Numeric CB annotations match the sidecar. Unsupported CB is **`N/A` plus a reason**, not zero or a dropped record. Missing CFD remains `N/A` / `upstream_cfd_unavailable` in the sidecar; default PAF retains its historical absent numeric tag.

Unavailable CB reasons:

- Resolved: **8 DNA2 over-length**, **8 multi-base RNA bulges**.
- N panel: the same plus **2 ambiguous-reference** rows.
- Mixed-only: **8 DNA2**, **8 RNA2**, **32 mixed I/D** rows. All 32 CFD-unavailable mixed representatives were retained.
- Retained172: **75 ambiguous-reference**, **20 RNA2**, **18 DNA2 over-length**, **16 mismatch-domain exclusions**.

These are first-failure categories. No unsupported model extension, truncation, projection, alternate alignment or retraining was used. Counts include exact control targets. Different panels intentionally share fixtures; **their row totals are not counts of distinct biological sites**.

### Official neural prediction agreement

All **1,177 supported panel-row predictions** were reconstructed independently from bounded reference coordinates/CIGARs and passed to unmodified official `ensemble_predict`, pinned upstream **`3eddcd5bfcaff00b2bdf29425116ec9756ace870`**, same five pretrained classification components. Both live CLI modes and repetitions agreed with these official-reference values within **`1e-6`**; maximum absolute difference **`1.789291992349007e-7`**. External chunk sizes/official batching differ, so bitwise float agreement is not required. The official check used the existing reviewed reference environment, while actual CLI inference used the new clean inference venv.

This is independent conversion/model-implementation agreement, **not biological cleavage accuracy, calibration, score interchangeability, or genome-wide recall**. Earlier [model-domain limitations](../crispr_bulge/REPORT.md#2-compatibility-with-mismatches-and-bulges) and [score-comparison caveats](../crispr_bulge/ASSEMBLY_EVALUATION.md#cfd-versus-crispr-bulge) still apply.

## 2. Complete-workflow benchmark

**Identical N-containing bounded input**, one guide, 570 hits, default search limits, one search/compute thread. Three fresh CLI processes per mode; order reversed on the middle repetition. Every measurement includes snapshot copying, index/reference verification, native loading/search, reporting, CFD, and (when selected) neural import/model loading/inference/annotations/output publication. One-time index construction, installation, independent oracle/reference validation and supervisor startup are not CLI work and are outside child wall time. Storage caches were warm/uncontrolled on shared `octopus01`; this is not a genome-scale throughput or speedup claim.

| CLI scoring mode | Median child wall seconds | Max individual descendant RSS KiB | Max sampled summed tree RSS KiB |
|---|---:|---:|---:|
| `cfd` | **0.214627** | **15,056** | **51,008** |
| `crispr-bulge` (upstream CFD retained) | **5.020807** | **538,428** | **589,128** |
| `both` | **4.970432** | **536,748** | **587,724** |

Unlike the previous saved-PAF measurements, the CFD row here includes **real search and CFD computation**. CB-requested mode intentionally retains existing CFD computation/tags, like `both`; it is not a pure neural-only kernel. Their small timing difference is run variability, not evidence that combined scoring is intrinsically faster. Fresh-model setup dominates these tiny fixtures. Exactly five model loads are reused across chunks/guides per invocation; retained172 also exercises 64-guide batching.

Measurement uses fresh-process `RUSAGE_CHILDREN` and a requested **20 ms** process-tree sampler, with a **4 GiB sampled stop** and **300 s outer timeout**; no stop fired. Actual intervals include sampling overhead. Linux maximum child RSS is an individual-descendant high-water mark, not concurrent summed physical memory. Summed RSS can double-count shared pages/miss spikes and includes supervisors. The stop is not a kernel memory cap. Full large-assembly live-search-plus-neural peak memory and huge/repetitive outputs were **not** measured.

## 3. Installation and clean-environment testing

See [RESEARCH_INSTALLATION.md](../../docs/RESEARCH_INSTALLATION.md) for versions, commands, output contracts, safe index/model preparation and platform caveats.

Successfully tested:

1. **Fresh `git archive 6fa96d3` source**, an **empty Cargo target**, Rust/Cargo **1.94.1**, `--locked --offline`; all **164 Rust tests passed**, with one pre-existing external-Python CFD test ignored. Existing trusted Cargo source/download caches were reused, not target objects.
2. WFA2-lib **`df3dce4...` / v2.3** freshly cloned from upstream and compiled; optional example/benchmark tools disabled using supported make flags. No installed production library was overwritten.
3. Fresh local clones of pinned Columba **edcc7ce / 2.0.3**, libsais **2.8.5**, fmt **10.2.1**, parallel-hashmap **e5b892b**; all libraries, builder and unpruned PAM executable compiled afresh with **GCC/G++ 12.3.0**, CMake **3.25.1**, zlib **1.2.13**. Dependency source caches were reused, **not prebuilt libraries/objects**. Vanilla/32-bit/alphabet-5/sparse-SA-4 contract preserved.
4. A **new CPython 3.11.15 venv** installed the original inference requirements, passed `pip check` and the complete workflow. No system-site packages were enabled. It is a fresh package environment, **not a freshly installed OS or Python distribution**.
5. Captured the **46-wheel Linux x86_64/CPython 3.11 hash lock**, downloaded a wheelhouse, installed a **second new venv offline with `--require-hashes`**, passed `pip check`, verified identical frozen versions, and ran the complete retained172 workflow. Its final PAF was **byte-identical** to the first clean environment: 172 rows, 97 CFD defined, 43 CB defined, 129 CB N/A. The wheelhouse stays external; the hash lock is committed.

Resolved installation/test pitfalls (raw failed attempts retained externally):

- A Guix profile under shell `nounset`, mixed system G++ 8 / Guix 12 headers, and missing development `-lrt` caused initial build failures. Explicit coherent compilers and matching glibc 2.35 static compatibility libraries fixed the clean Rust build. WFA optional benchmark tools were unnecessary; library-only build options avoid unrelated tool linkage requirements.
- Two old compatibility scripts built **unbound** raw indexes and correctly failed the current native verifier. Their fixtures now use `index_manifest.py build`; integrity enforcement remains unchanged. All 96 option/PAM comparisons and mixed-gap checks then passed.

**Unresolved network model setup:** the real `fetch_models.py` helper was attempted, not mocked. Git **2.41.0** aborted internally during partial fetch (`config.c:129 kvi ...`, SIGABRT); Git **2.20.1** retry failed HTTP 400/RPC negotiation. No TLS bypass was used. Network model installation is **not certified**. For execution, the five already reviewed checkpoints were copied into a new directory and all pinned SHA256s verified. This is the tested trusted-local-model research path, not a successful clean network bootstrap. The larger official-reference Python environment was reused rather than freshly reinstalled.

## 4. Regression tests and reproduction

Authoritative final script run completed successfully:

- **164 Rust tests passed**, one pre-existing ignored.
- **20 configured scoring Python tests passed**, none skipped: all CLI modes, tiny single/multi-guide reporting, prior CFD compatibility, unsupported DNA2/RNA2/mixed, N/A versus zero/missing CFD, empty input, checksum failures, worker failure/timeout, model reuse and pinned-reference-offset checks.
- **4 new validator tests passed**: duplicate/query identity, missing/unexpected anchors, representative/CFD mutations, unavailable versus numeric zero.
- **9 prior offline CFD model tests passed**; these do not enable production pruning.
- **96 native PAM/option oracle comparisons**, five invalid-option and two invalid-guide-batch checks passed; binding/tamper/PATH-spoof checks and **two mixed-gap strands** passed.

Run [reproduce.sh](reproduce.sh) with separately prepared development tools and a **new absolute external output directory**; exact environment variables are in the installation guide. The driver has no arbitrary genome input option. It creates only tiny synthetic bound indexes plus the optional exact-hash retained172 panel. All three modes include actual live search. Optional official-reference arguments are explicit and their presence/absence is recorded. Repetition does not launch a whole-genome query or use the pruning prototype.

Small committed evidence: [workflow summary](results/workflow-summary.json), [installation/provenance](results/installation.json), [offline install execution](results/offline-validation.json), [Rust tests](results/rust-tests.log), [scoring tests](results/scoring-tests.log), [validator tests](results/validator-tests.log), [CFD model tests](results/cfd-model-tests.log), [PAM/option tests](results/pam-option-tests.json), [binding checks](results/binding-tests.json), [mixed gaps](results/mixed-gap-tests.json).

External root: **`/export/local/home/salehi/crisprapido-research-workflow-v1/`**. Retain:

- `source/`: frozen clean source archive (executables embed helper paths);
- `wfa/`, `columba/`, `columba-build-v2/`, `target-final/release/crisprapido`, `pam-search`;
- `inference-env/`, `offline-env/`, `wheelhouse/`, `models-verified/`;
- `final/`: authoritative complete reproduction, full per-mode PAF/annotations/metrics/logs, oracle/official requests and test logs;
- `offline-smoke/`, `offline-smoke.paf`, `offline-validation.json`;
- `clean-build.sh`, `clean-columba.sh`, failed-attempt logs and earlier exploratory `workflow/`, `reproduced/` directories. The earlier reproduction stopped on the stale unbound-test fixture; **only `final/` is the authoritative all-tests-passing run**.

No environment, wheel, checkpoint, source dependency checkout, FASTA/index, sequence-bearing request, PAF or large output is committed. Core source/algorithm/default/matrix hashes were checked against `6fa96d3` and the clean archive.

## 5. Release blockers / stop point

Research validation is not a release decision. Before considering a release:

1. **Reliable model distribution/setup**, checkpoint/dataset/dependency rights and security/lifecycle review. Fresh network helper installation failed here. Columba's checked root license is AGPL-3.0, this repository declares an academic AGPL/commercial arrangement, and WFA2/CRISPR-Bulge root licenses are MIT; the entire linked/distributed stack needs review, not an assumed blanket MIT grant.
2. **Portable coherent build/install packaging** (including helper relocation, Columba internal API/index layout and tested OS/CPU/toolchain matrix); current source-bundle paths and site-specific Guix linking are not a portable release environment. Rust MSRV is not established.
3. **Explicit partial-score contract and scientific validation**: DNA2 exceeds the checkpoint shape; RNA2/mixed/high-mismatch/ambiguity remain unvalidated or outside conservative policy. Do not fabricate missing scores, equate CFD/CB calibration, or silently change defaults/representatives.
4. **Wider bounded/randomized recall review and biological benchmarking**. Finite oracle recovery and official code agreement do not establish arbitrary genome-wide completeness or biological accuracy. Existing default-backend limitations remain documented and unchanged.
5. **Resource/lifecycle/concurrency hardening**: huge outputs/repeats, live large-reference memory overlap, whole-reference pre/post hashing/snapshot costs, hard memory/global timeout behavior, trusted index/model provisioning and malicious mutate-and-restore races remain open.

Do not enable integrated CFD pruning, implicit score filtering, new model-domain extrapolation, main merge, production binary replacement or further whole-genome experiments automatically.
