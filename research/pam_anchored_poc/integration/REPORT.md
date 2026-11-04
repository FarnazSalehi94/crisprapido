# Optional experimental PAM-anchored backend

The development source now accepts `--search-backend pam-anchored` and keeps `columba-wfa2` as the default. This is an experimental, separately built process adapter. The existing candidate-generation/WFA2 path remains intact. The prototype emits one representative per PAM boundary, chosen by `(mismatches + inserted bases + deleted bases, mismatches, compressed CIGAR, start)`. The adapter preserves that representative and the existing PAF formatter and downstream CFD scorer. It never selects on CFD. There were 54 differing representatives versus the existing pipeline on the 204 panel, and 41 of those had different CFD scores; both sides of all 54 pairs were independently eligible (see [alignment review](../alignment_review/REPORT.md)). CFD is not a calibrated probability of cleavage.

## Usage and supported scope

Build the native research executable with [build.sh](../build.sh), against the recorded Columba checkout. Build a **fresh bound** index with `python3 research/pam_anchored_poc/compatibility/index_manifest.py build --builder COLUMBA_BUILD --reference REF.fa --index PREFIX`. Then run, for example:

```
crisprapido --search-backend pam-anchored --pam-search-bin /path/to/poc \
  --columba-index PREFIX -r REF.fa -G GUIDES.fa -t 1
```

Defaults remain `-m 5 -b 1 -z 2 -f 0.75 -p GG`. Resolved two-base PAM suffixes, nonnegative mismatch limits, b=0/1, z=0/1/2, fraction 0..1, both strands, and single or batch guides are supported. The prototype uses one thread. It rejects multiple threads, `--window-size`, imported SAM, baseline Columba executable options, `--no-filter`, noncanonical PAM mode, more than one contiguous bulge event or bulge size above two, and ambiguous guides. Reference FASTA must be uncompressed with unique record IDs and A/C/G/T/N bases. Reference N within a target is a mismatch unless consumed by a deletion; PAM wildcard must be a resolved A/C/G/T base and required PAM suffix bases must match exactly. Unsupported or malformed input fails rather than returning partial results. The existing backend retains its options and default behavior.

The executable requires Python 3 and the trusted runner, binding helper, and native binary's build-pinned helper and interpreter at their build locations. The runner uses `/usr/bin/python3` by default, or an absolute interpreter path set as `CRISPR_PAM_PYTHON` when building the Rust binary; it does not search PATH. The runner's source path is embedded at Rust build time, so relocating or distributing the executable requires preserving that path or packaging work. Native components must be trusted. The prototype has a 32-bit index build and sparse-SA-4 assumption; this is not a general-purpose or production-ready index interface. The runner currently buffers output and reference sequences and copies the complete index, so storage and memory requirements make genome-scale use inappropriate without further work.

## Binding and changing files

For every invocation, the adapter copies the FASTA, complete manifest and index components into a private temporary directory, verifies their binding, runs the native executable (which independently verifies binding), and re-verifies the snapshot before accepting results. The target sequences used for scoring come from that snapshot. A changed or missing component, missing binding, mismatched FASTA or failed verifier aborts without reporting hits. Copying rather than linking isolates the search from ordinary changes to the source files after copying; a change during copying is caught when the snapshot's hash fails to match the binding. Temporary files are made read-only. This is not an authenticated manifest or a defense against a same-user adversary modifying the temporary files, executable, helper, interpreter or manifest, and a check/use race remains if such an actor can change a snapshot and restore it between checks. Trusted inputs and process isolation remain necessary. Whole-file copies and repeated hashing are included in measured total time.

## Validation

The independent PAF auditor [audit.py](audit.py) extracts reference and PAM, replays each CIGAR, checks full consumption, base identities, coordinates, counts, limits, output fields, and independent CFD matrix results. It compares each PAM-start key to the exhaustive oracle. The 204 resolved panel recovered **204/204**, and the 172 N panel **172/172**, with zero unexpected and duplicate anchors. CFD agreed for all 204 scored resolved records and 97 scored N-panel records; 75 N-containing records correctly had no score. The existing pipeline also recovered and passed audit on both panels. The integrated 204 representative coordinates and CIGARs matched the native prototype exactly.

[test_adapter.py](test_adapter.py) compared the integrated records and native prototype representatives under all 16 PAM suffixes and six supported parameter combinations (96 runs), including both strands, against a fixture already independently enumerated in [test_options.py](../compatibility/test_options.py). All matched. It also checked rejection of invalid or unsupported combinations, an altered FASTA and a tampered index with a fake Python executable on PATH. The independently checked mixed adjacent I/D fixture produced both expected strand anchors with `m=0,b=1,z=2`, gap count 1 and gap size 2. The integrated output also matched the independently checked native representatives on the three difficult PAM panels (33, 3 and 4 sites) and the alternative-alignment/wildcard fixture (26 records), via [test_fixtures.py](test_fixtures.py). Earlier exhaustive random, reference-boundary and N substitution tests cover the native search; see [alignment review](../alignment_review/REPORT.md) and [compatibility report](../compatibility/REPORT.md). The random, boundary and N-substitution suites have not all been repeated end-to-end through the adapter. The Rust release test suite also passed (160 tests, one previously ignored). Bounded testing does not prove general search completeness.

## Bounded complete-workflow performance

Two single-thread runs on the same host used the same prebuilt 8 kb reference/index, 64 guides, m5/b1/z2/f0.75/NGG, and the same PAF and CFD reporter. Both reported the same 172 independently audited anchors. [Measurements](results/benchmark_172.json) and [runner](benchmark.py):

| Backend | Total wall seconds | Native indexed search seconds | Maximum individual process RSS KiB | Anchors |
|---|---|---|---|---:|
| Existing Columba + WFA2 | 12.37, 11.49 | not isolated | 22,444, 22,392 | 172 |
| Experimental | 0.42, 0.43 | 0.3591, 0.3589 | 14,600, 14,552 | 172 |

Each experimental search explored 4,229,884 states. Baseline search states are unavailable. GNU time's RSS is the maximum for an individual process, not aggregate process-tree memory. Total time includes copying, binding checks, search, SAM generation/parsing where applicable, verification, scoring and output; the benchmark excludes the separate post-run audit. The search times are not comparable across algorithms or phases. This panel's end-to-end experimental runs were faster, but the two implementations choose different eligible alignments and internal work, and two repetitions on one bounded panel do not establish a general or genome-scale speedup. Baseline indexed search time was not independently measured in this run.

## Remaining work

The prototype remains experimental: broaden end-to-end oracle coverage and randomized parameter tests, assess sensitivity and resource consumption on larger bounded windows, package and pin its native executable/helpers and dependency versions, decide stable representative policy and consequences for CFD, provide scalable binding/immutable index handling and process-tree memory profiling, address single-thread and index format limitations, and obtain independent review of completeness and licensing before promotion. No CFD-guided pruning or CRISPR-Bulge integration was added.
