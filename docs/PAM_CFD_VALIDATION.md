# PAM/CFD correction and validation

## Behavior contract

`-p GG` now requires a complete three-base Cas9 PAM whose last two bases
are GG, in guide orientation, on either strand. The wildcard first base must
be a resolved A/C/G/T base in the assembly. An incomplete or ambiguous PAM
is unavailable, not a canonical match. `-p` accepts two A/C/G/T bases,
case-insensitively; it does not accept an IUPAC pattern or the string NGG.

Use `--allow-noncanonical-pam` to keep candidates with other or unavailable
PAMs. This preserves permissive candidate eligibility while correcting score
inputs. CFD uses the actual observed final dinucleotide, never the requested
motif as a constant. With an unavailable PAM, a permissive row has no `cf` tag.
`--no-filter` also bypasses PAM eligibility and sequence thresholds; overlap
representative selection still applies. Neither flag reproduces buggy scores.

PAM eligibility is applied after interval verification and before choosing
representatives in overlapping groups. Candidate generation remains based on
sequence edit distance. Alternative-interval recovery continues to consider
nonmatching and unavailable PAMs. `-z` controls gap size, not PAM relaxation.

For zero-based, half-open target interval `[s,e)`, the full PAM is reference
`[e,e+3)` on plus and reverse-complement of `[s-3,s)` on minus. Only then is
the final two-base suffix used for CFD. Missing third bases do not qualify.

Output query coordinates cover the complete guide (`0..guide_length`).
Reference coordinates cover every reference-consuming `M/=/X/D` column;
`I` consumes query only. Terminal `D` bases remain in both the reference
interval and CIGAR; they do not advance query start or move target start.
Both scanning backends now align the original guide to an oriented target,
so CFD and `cg` use guide orientation on both strands. This local `cg`
convention for negative alignments is retained from the Columba backend;
consumers replaying against forward-reference sequence must orient the target
first. It is not a claim that an unmodified minimap2 CIGAR consumer will use
that same convention.

## Defects corrected

| Defect | Cause | Correction |
|---|---|---|
| Wrong CFD PAM factor | Two bases immediately beside protospacer were used; these are PAM positions 1–2 | Require three bases, orient once, pass positions 2–3 |
| Non-NGG rows survived `-p GG` | Requested motif only influenced overlap preference and fallback exploration | Strict motif eligibility before overlap selection; explicit permissive flag |
| Negative CC treated as requested GG | Matching reverse-complemented an already oriented pair a second time | Compare oriented pair directly on both strands |
| Shifted leading-D intervals | Reporting advanced start by leading D and also counted D in the full span | Preserve verified start; count all reference columns once |
| Query endpoints exceeded length | Reporting incremented query start for reference-only D | Query start stays zero; only query-consuming columns advance end |
| Negative window CFD differed by backend | Direct scanner aligned reverse-complement guide to forward target | Align original guide to reverse-complemented window, map coordinates back |
| Negative hits suppressed by plus hits | Direct scanner returned immediately on a successful plus alignment | Evaluate both orientations independently in each window |
| Compact CIGAR conversion lost operations | Formatter iterated characters and skipped `=`/counts | Parse operations and merge runs, accepting raw WFA and compact CIGARs |
| Compact Hit span/quality ignored counts | Hit helpers counted operation characters | Reuse parsed lengths for overlap endpoint and quality |
| Malformed CIGAR accepted | Zero length or unsupported operations could pass helper parser | Reject empty, zero, unsupported and overflowing operations |

The branch also carries the preexisting uncommitted CRISPRme-compatible CFD
DNA-bulge projection and its tests. It was already present in the production
binary; preserving it avoids comparing the PAM fix against an older scoring
algorithm. The original dirty checkout remains untouched. Its zero-byte
`main.rs` and `reporting.rs` were not copied: the isolated worktree starts from
their committed versions. The original work is archived separately.

## Automated validation

Final validation: **148 tests passed, zero failed, one existing live-oracle test ignored**.
`cargo fmt --all -- --check` and `git diff --check` pass.

The unit/integration suite covers positive/negative canonical and noncanonical
PAMs, observed CFD input, permissive behavior, missing/ambiguous PAM bases,
contig edges, lowercase, perfect/one/multiple mismatches, DNA/RNA bulges,
terminal deletions, raw/compact CIGARs and coordinate conservation. Direct
window tests place a hit on each strand in the same window. CLI tests use
production filtering rather than relying solely on the older `cfg!(test)`
scan-threshold override.

A separate existing live standalone CFD oracle test is ignored by default
because it requires external pinned Python source/pickles. Its stored numerical
fixtures and the new independent subset oracle run in the normal suite/review.
No CRISPR-Bulge implementation or environment was added.

## Real subset replay (2026-10-28)

A deterministic sample of 120 records from the previous HG00128_hap2 audit
includes both strands, perfect/mismatched alignments, RNA/DNA bulges and
leading/trailing reference deletions. Seed: 20261028. All available alignment
categories contribute up to five rows, then sampling fills to 120. This is a
representative debugging set, not a population-frequency sample.

`scripts/validate_pam_subset.py` fetches small windows around these records
from the existing BGZF/indexed assembly. For historical leading-D records it
reverses the known start shift before constructing candidate intervals. It
creates a new tiny FASTA/SAM and index, then runs old, corrected strict and
corrected permissive binaries on those same candidate intervals with one
thread, `-m 4 -b 1 -z 1`. No assembly scan or index construction is performed.
The old binary SHA256 is
`ee543e772ac6e956e8fc3eb88e6ea52517e5192eff9da255902c3bad6073faed`.

| Replay | Emitted | NGG | Non-NGG | CFD ≥0.20 |
|---|---:|---:|---:|---:|
| Old binary | 120 | 32 | 88 | 116 |
| Corrected strict | 88 | 88 | 0 | 50 |
| Corrected permissive | 120 | 88 | 32 | 50 |

By input candidate ID, strict mode removes 32 representatives. Of the 88
survivors, 20 scores increase, 37 decrease and 31 are unchanged; 57 intervals
change. Permissive mode retains all 120 IDs: 20 scores increase, 69 decrease
and 31 are unchanged; 63 intervals change. Corrected overlap preferences and
fallback intervals can select a different alignment for the same candidate,
so changes cannot be interpreted as a simple motif filter at fixed endpoints.

Every corrected output is independently fetched from the tiny indexed FASTA
using samtools, oriented, replayed against its guide using CIGAR, checked for
query/reference consumption and exact-match/mismatch assertions, and scored
by independent Python matrix multiplication (including the existing DNA-bulge
projection). All 208 corrected rows agree with CFD to output-rounding
precision; all 88 strict rows have NGG. The replay is not a reproduction of the
original production search: local window boundaries, imported intervals,
historical reporting repair and representative selection affect its outputs.
It does not replace the prior 6,123-record PAM audit.

Validation artifacts are in the analysis workspace:
`analysis/PAM_fix_validation_2026-10-28/`. The final replay is `real_subset_final/` and has per-row source
regions, FASTA/SAM/index, commands, binary hashes, old/new PAFs, independent
sequence/score verification and JSON counts. The validator refuses to reuse
an output directory and caps the requested sample at 1,000.

## Historical result impact and regeneration

All results produced with the identified old binary and its two-base extractor
are affected, including the Brunello CFD020 retained outputs and CFD-based
summaries/rankings. Canonical PAMs can gain scores while GGN contexts lose them.
Strict mode also changes eligibility and overlap representatives. Results
from other binaries must be classified by hash/source, not filename alone.

- For corrected permissive analysis, regenerate CFD and overlap selection from
  all available unfiltered candidates with exact guides and source assemblies.
- For canonical-only analysis, regenerate using corrected strict behavior;
  do not just remove non-NGG from old retained files.
- Current retained plus unscored sidecars omit scored-below-threshold rows.
  They cannot recover canonical candidates lost under the old scoring rule.
- Query/target-coordinate consumers and negative direct-window results need
  regeneration too; even perfect-match CFD is not universally unchanged.
- Recompute candidate counts, threshold summaries, downstream comparisons and
  biological conclusions that used these outputs. Preserve old artifacts with
  their original binary hashes for provenance.

No large regeneration was launched here. This branch does not resolve
completeness limitations of the direct scanner (one best alignment per strand
per window) or introduce a new biological model for multi-bulge CFD. Full
production reruns require separate scope and resources.

## Reproduction

Standard build (with a compatible existing WFA2 library):

```sh
WFA2LIB_PATH=/absolute/path/to/WFA2-lib cargo test --locked -j 2
cargo fmt --all -- --check
python3 scripts/validate_pam_subset.py --help
```

The validation host needed Rust 1.85.1 and explicit Guix glibc/GCC runtime
paths because its default runtime differs from the libraries selected by its
linker. These environment-specific settings are recorded only in the analysis
`cargo-local.sh`; no repository toolchain default or production runtime changed.
