# PAM-anchored Columba proof of concept

Research branch based on CRISPRapido `32a56e9`; no production implementation or binary was modified. This is a small strict-NGG experiment, not a replacement for CRISPRapido. It directly invokes Columba's FM-index primitives and does not use WFA2, Columba's alignment routines, or CFD during search or validation.

## Runtime provenance

The release manifest records the [biointec/columba](https://github.com/biointec/columba) source checkout at [`edcc7ce2a1e28ded6784a44afbd8b9348a4805ad`](https://github.com/biointec/columba/commit/edcc7ce2a1e28ded6784a44afbd8b9348a4805ad), 2026-03-27. Its `origin` is that repository, tracked tree is clean, and the untracked `build_Vanilla` directory contains the compiled library and executable. The release executable announces **Columba VANILLA 2.0.3**. Its SHA256 is `7e44974c81dad4943292fa64a71e54073908b4909ce327480464d769322c40d1`; the loader wrapper SHA256 is `2367ef7b053a44c83f76b059ed66a8f671f83773d3fe69240cda50893408f28a`, matching the [release manifest](../../docs/release_candidate_2026-10-28/release_manifest.json). That manifest records the checkout revision, but the untracked build and absent reproducible build attestation mean we cannot prove the binary was compiled from exactly that tree or exclude unrecorded build-time differences. The hash identifies the actual runtime.

The local `origin/main` tracking ref equals `edcc7ce`, but is not evidence of the current remote: a fresh GitHub fetch failed due to local certificate configuration, and a web view did not expose a current head SHA. A precise current-upstream diff is therefore unavailable. The prior feasibility report examined upstream main via the web; the recovered 2.0.3 tree has the same major index/search concepts, but no claim of source identity with present upstream main is made.

The prototype links that checkout's vanilla `libAuxiliaryClasses.a` and uses `FMIndex` (`src/fmindex/fmindex.h`, `.cpp`), `IndexInterface::setDirection` and `addChar` (`src/indexinterface.h`, `.cpp`), `SARangePair` and `findSA`. It does not use search schemes or alter Columba source. The build was 32-bit index coordinates (`THIRTY_TWO`, `ALPHABET=5`, sparse SA 4). The measured POC executable SHA256 is `8aa7327b52e5a8f61f80994f96a302dbd68a88f558830f39c2b42be4503bb9da`.

## Algorithm and output

For each guide and strand, search four exact PAM seeds (AGG, CGG, GGG, TGG); backward plus-strand extension and forward minus-strand extension then consume the guide from its PAM-proximal end. PAM N represents A/C/G/T, with no PAM mismatch or gap transitions. A state carries the paired FM ranges, remaining guide position, mismatches, gap direction/open or closed status, gap length, target span, matches and a traceback. A DFS branches over four resolved reference letters and one contiguous insertion or deletion run. Only empty FM intervals, mismatch limit and an upper bound on attainable match count are pruned. The index collapses repeated substrings; completed intervals are located with `findSA` and checked against the supplied FASTA for contig bounds, resolved PAM and reference/CIGAR consistency, important because Columba indexes substitute non-ACGT bases. Only one deterministic representative is retained per guide, contig, strand and PAM-adjacent boundary. The CIGAR is in guide orientation and excludes PAM. Coordinates are zero-based half-open: plus anchor is protospacer end; minus anchor is protospacer start. `I` consumes guide only; `D` consumes reference only. Output includes target span, observed PAM, mismatch count, and inserted/deleted bases.

Command arguments configure maximum mismatches, maximum length of **one** gap event, and minimum match fraction. The tested settings are 5/2/0.75; 0/2 and 4/1 also exercised. This POC does not expose production defaults or change them. It requires resolved guides and strict NGG, rejects ambiguous target spans, assumes the reference FASTA exactly matches the indexed concatenation/order, and uses Columba's sparse-SA-4 index. It has no generic PAM, noncanonical-PAM or `--no-filter` compatibility, no production CLI integration, no affine-score selection or CFD. Representative choice differs from WFA2; 54/204 panel anchors have a different `(start,end,CIGAR)` representative. Independent row audits, not literal CIGAR equality, determine validity. It does not retain all equally valid alignments in output.

## Correctness results

The test uses retained [8 kb independent panel](../../docs/real_defaults_validation_49f20fc_2026-11-02_v1/REPORT.md), four windows, 96 guides, strict NGG, m=5, b=1, z=2, f=0.75. The independent oracle scans each PAM and enumerates full-guide ungapped and one-/two-base single-gap placements. It identifies anchors using PAM start, so the minus-strand prototype boundary is converted by subtracting three for comparison. The saved baseline is Columba + WFA2 at `49f20fc`, identical limits and same reference, already independently audited in the cited report. The recorded output is [panel.tsv](results/panel.tsv); Columba emits two startup information lines before its TSV header.

| Result | Existing baseline | Prototype | Independent oracle |
|---|---:|---:|---:|
| Unique expected anchors recovered | 204/204 | 204/204 | 204 |
| Unexpected anchors | 0 | 0 | — |
| Duplicate anchors | 0 | 0 | — |
| Plus / minus anchors | 108 / 96 | 108 / 96 | 108 / 96 |

The baseline and prototype anchor sets are identical. `check_panel.py` independently reconstructs every prototype alignment from guide, reference and CIGAR, checks full consumption, observed PAM, anchor, mismatch/bulge limits and match fraction, and compares the oracle set. The independent panel oracle and audit are not derived from this DFS. On the previously retained difficult PAM-anchor panels, the prototype matched the entire saved oracle anchor sets: synthetic 33/33, CHM13 3/3 and 12-window 4/4, no extras or duplicates, using the original m=4, b=1, z=1, f=0.75. Three existing fallback examples (terminal and near-terminal two-base guide insertions on plus, terminal insertion on minus) were reproduced under m=0, b=1, z=2: 3/3. The two near-terminal gap normalizations need not match the WFA CIGAR literally. The noncanonical-PAM fallback example is outside this strict-NGG POC. Additional synthetic contigs exercise A, C, G and T at PAM N on both strands (8/8), including contig-edge PAMs.

These are bounded examples and do not establish general completeness. In particular, Columba's index builder randomly replaces ambiguous bases: checking the original FASTA prevents reporting a substituted base as real, but can also prevent discovery of eligible ambiguous-reference matches. Search only considers A/C/G/T text branches. Completeness is conditional on resolved indexed sequences, correct index/reference pairing, available resources, and the tested one-gap/full-guide model; these assumptions have not been proved for every parameter or arbitrary input. The `findSA`-to-contig conversion should be replaced with robust index metadata and reference identity verification before wider use. Terminal gaps, repetitive sequences, boundary behavior, alternative traces, wildcard behavior beyond these examples, and randomized exhaustive small-reference cases require broader tests. Input bounds and large-index overflow need hardening. Columba index memory, AGPL obligations, and divergence from production WFA selection need review. CFD remains outside the prototype.

## Bounded timings

After correctness checks, three single-thread runs on the 8 kb panel, existing index, output redirected, gave **0.70 / 0.70 / 0.70 s wall**, 0.70 s user each, peak RSS **4,424 / 4,176 / 4,316 KiB**, **8,313,485 DFS states** and 204 emitted anchors each. `states` counts every `walk` invocation including states pruned immediately; it is not directly comparable with Columba's node counter or WFA calls. Index construction is excluded. The retained single baseline run was **13.44 s wall, 13.39 s user, 35,520 KiB maximum RSS**, 204 anchors; its command, binary hash and measurement are in [search.time](../../docs/real_defaults_validation_49f20fc_2026-11-02_v1/search.time). Baseline includes Columba subprocess discovery, SAM production/reading, WFA verification, reporting and CFD, whereas the prototype emits only research TSV and does not score CFD. RSS for a process and separately spawned subprocesses is not necessarily combined. Runs were not interleaved and timing environments may differ. This is not an apples-to-apples speedup or a genome-scale forecast.

## Reproduction

Build against an existing clean Columba `edcc7ce` vanilla build:

```sh
COLUMBA_DIR=/path/to/columba COLUMBA_BUILD=/path/to/columba/build_Vanilla research/pam_anchored_poc/build.sh
P=docs/real_defaults_validation_49f20fc_2026-11-02_v1
# Build a fresh small index first: columba_build -r /tmp/panel-index -f "$P/ref.fa"
/tmp/crispr-pam-poc /tmp/panel-index "$P/ref.fa" "$P/guides.fa" 5 2 0.75 > /tmp/panel.tsv
python3 research/pam_anchored_poc/check_panel.py "$P" /tmp/panel.tsv
COLUMBA_BUILD_BIN=/path/to/columba_build POC_BIN=/tmp/crispr-pam-poc research/pam_anchored_poc/run_fixtures.sh /tmp/pam-poc-fixtures
```

On this host explicit Guix loader/library wrappers and GCC 12 paths were necessary; no build products are committed. `run_fixtures.sh` creates fresh small indexes and recreates fallback references from `src/verification.rs`. The test outputs and retained oracle are the acceptance criteria, not a production benchmark. Future work: independent code review, randomized exhaustive differential tests, reference/index binding, robust configuration and alternative-alignment selection, then bounded larger benchmarks. No whole-genome or pangenome experiment, merge, default change, production binary change, frozen-branch modification or CRISPR-Bulge integration was performed.
