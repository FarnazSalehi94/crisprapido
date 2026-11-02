# Bounded release validation — 2026-11-02

**Recommendation: C — needs additional bounded testing and independent review.**
The new held-out panels have 100% measured PAM-key recall and precision, and
all bounded regional searches completed. These tests do not establish that
Columba always supplies an overlapping raw candidate, or that the full
100-guide × 12-assembly experiment has an acceptable runtime/memory envelope.
No independent reviewer has examined this package. The frozen preflight gate
therefore **does not pass**; no full search was launched.

## Frozen implementation

Source commit `0aa9944e2eae4b7f8286415ee5c3fb3667cf6097`, tree
`45873bc889abc9366523f1089fbf7190d07c9261`; the isolated optimized executable
SHA256 is `c37f4ded6e7c6afb3a7d22eccc5566f291244731bf6872dbfbd0b4d47f53aa74`.
No application algorithm was changed in this validation turn. The manifest
`freeze_manifest.json` records Cargo.lock, score matrices, Columba 2.0.3,
WFA2, Rust 1.85.1, samtools 1.23.1, runtime dependencies and artifact hashes.
The test settings are strict `-p GG -m 4 -b 1 -z 1 -f 0.75`, Columba e=5 plus
separate e=0. The bounded pilot uses `-t 1`; the original full plan specifies
four threads. No relaxed eligibility or score cutoff was used.

## Independent held-out sensitivity

New synthetic fixtures use seed 20261103, two different frozen-guide sequences,
a different repetitive guide, fresh random flanks, both strands, exact sites,
mismatch positions 2/6/11/16/19, I and D at 0/3/8/13/19, four-mismatch and
combined mismatch/bulge cases, and adjacent NGG anchors. No complete synthetic
contig sequence is shared with the two recorded development panels. The real
fixtures contain four unchanged 1 kb windows from CHM13 and GRCh38 at
1,300,001–1,301,000 and 1,700,001–1,701,000 on their first two contigs,
different from prior correction windows. Sixteen deliberately designed real
queries include exact, four mismatch, I-oriented and D-oriented guides. These
are enriched fixtures, not unbiased population samples.

The independent direct-enumeration oracle scans both strands, every candidate
19/20/21-base span, and all one-base gap positions under frozen limits and
complete oriented NGG. Recall counts distinct `(guide, contig, strand, PAM
interval)` keys; precision counts emitted keys also in oracle truth. It does
not treat overlap as identity or claim two alignments at a key are physically
equivalent. Category membership can overlap. Raw Columba SAMs, oracle truth,
PAF, metadata and full independently audited summaries are retained in the
held-out directories.

| Panel | Oracle PAM keys | Recovered | Precision | Output strands |
|---|---:|---:|---:|---|
| Synthetic | 212 | 212/212 | 212/212 | 106 plus, 106 minus |
| Real windows | 18 | 18/18 | 18/18 | 18 plus, 0 minus |

There were **zero missed and zero unexpected** PAM keys in either panel.
Both-strand evidence comes from the synthetic panel; no natural minus-strand
site was observed in these four real windows. Among oracle keys, 96 synthetic
and five real pairs have overlapping genomic PAM intervals for the same guide,
contig and strand; all distinct keys were recovered. The independent emitted
row auditor passed all 230 rows, with no duplicate output PAM key. The strict
interval/category metrics remain separately in each `metrics.jsonl`; they
include alternate intervals and cannot be read as locus recall.

For the 413 synthetic and 38 real oracle-valid intervals, every interval had
at least one overlapping same-guide, same-contig, same-strand raw candidate
in e=5/e=0 SAM. `candidate_overlap.json` records this check. This is an
observation on these fixtures, **not a proof** that Columba cannot miss an
eligible site without an overlapping raw candidate. No such case was observed
here, but it remains an unresolved sensitivity risk.

## Bounded multi-assembly scaling

Four selected 2 Mb FASTA windows were indexed separately: chr1:100000001–
102000000 and chr9:34000001–36000000 in each of CHM13 and GRCh38. The first
four frozen guides were run **individually** for 16 search processes. Each
search and its diagnostic profile had a 90-second wall timeout and 4 GiB
per-process address-space limit; one thread was used. No candidate-count limit
was imposed. Every profiled PAF matched the frozen executable's PAF byte for
byte. All 43 emitted rows passed independent audit. The profile copy is
instrumented and separately hashed; WFA call counts are from the instrumented
runs, whereas timing/RSS are from the uninstrumented frozen binary.

The chr9 windows were selected as a higher-repeat **proxy**, not from a formal
repeat annotation. In nonoverlapping 20-mer samples, fractions belonging to
repeated 20-mers were 3.92% CHM13 and 3.72% GRCh38 versus 1.23% and 1.32%
in the chr1 windows. `repeat_proxy.json` includes sample sizes and frequencies;
this simple metric does not classify all biological repeats. The selected chr9
windows did not necessarily produce more candidates for these particular guides.

| Assembly / selected window | Mean seconds per guide | Unique raw candidates, four guides | WFA calls | Reported PAM keys | Max RSS KiB |
|---|---:|---:|---:|---:|---:|
| CHM13 chr1 | 0.583 | 1,010 | 111,810 | 9 | 55,632 |
| CHM13 chr9 | 0.503 | 685 | 75,595 | 12 | 57,308 |
| GRCh38 chr1 | 0.580 | 1,003 | 111,018 | 10 | 55,612 |
| GRCh38 chr9 | 0.503 | 682 | 75,170 | 12 | 57,044 |

Per guide/window, measured times ranged 0.45–0.72 s, e=5 candidate counts
126–378, WFA calls 13,994–42,285, WFA calls per unique raw candidate
109.3–111.9, reported keys 1–6, and peak RSS 55,292–57,308 KiB. There were
3,380 unique raw candidate keys, 373,593 WFA calls, and 43 output rows in
aggregate; no e=0 mapped row occurred in these regional inputs. Full per-guide
rows, timed stages, profiler counts, repeat proxy and all 16 auditor summaries
are in `scaling/`. One query per process repeatedly loads the small index;
these per-guide times are **not** equivalent to full-plan batch throughput.

This is a small non-random panel: two closely related assemblies, four guides,
four windows, one thread and small indexes. It shows candidate-density and
guide variation but cannot establish an acceptable genome-scale resource
bound or repeat worst case. The earlier 10 Mb pilot had 14,807 unique candidates
and 1.63 million WFA calls for ten guides. A linear 196-hour serial estimate
from that pilot remains only an illustrative, uncalibrated figure; full indexes,
reference loading, multi-guide batching, threading, candidate distribution,
I/O, auditing and scheduling make the true cost uncertain. Small-index RSS
here does not represent whole-assembly index or reference memory.

## Remaining risks and gate

- **Columba discovery:** its `-a all`/e=5 calls and these finite tests do not
  establish that every oracle-valid site yields an overlapping candidate.
  The regional runs do not exhaustively oracle-scan the full 2 Mb windows for
  imperfect matches. Non-ACGT/reference-index transformation and extreme repeat
  cases remain additional concerns.
- **Expansion:** for candidate span C, guide length L and B=b*z, the implemented
  start/length loops are bounded by overlapping geometry and configured limits;
  at L=20,b=z=1,C=20 there are at most 120 combinations before exclusions.
  This is conditional on an overlapping raw candidate and can be costly; large
  non-frozen parameter values are not characterized.
- **Alternative alignment:** for the frozen b=1,z=1 path, WFA2's failed choice
  is followed by finite ungapped or all one-base-gap placements, checked against
  m, f, gap count and size; the optimization was compared against the prior
  implementation across resolved/ambiguous cases. For other b/z settings this
  is not a general exhaustive constrained aligner, and direct-window scanning
  retains its own limitation. Strict PAM and CFD logic were unchanged.
- **Repeats/resources:** no timeout or memory limit was reached in these 16
  processes, but there is no validated acceptance threshold and the test did
  not approach worst-case repeats or whole-assembly memory. The earlier
  10 Mb pilot and the 110 WFA calls per candidate warrant further bounded,
  stratified profiling rather than a claim that growth is acceptable.
- **Review and completeness:** independent review is outstanding; the frozen
  all-guide/all-assembly exact oracle and imperfect sensitivity remain incomplete.
  `REVIEW_PACKAGE.md` is a concrete package for an independent reviewer, not
  an assertion that review occurred.

All 157 standard Rust tests, the normally ignored pinned CFD test, and all 51
Python tests passed. No application source or production binary was changed;
no full search, merge or CRISPR-Bulge integration occurred. Large FASTAs,
indexes and per-run raw files remain in the separate versioned directory
`/export/local/home/salehi/crisprapido-release-pilot-0aa9944-20261102/`.
