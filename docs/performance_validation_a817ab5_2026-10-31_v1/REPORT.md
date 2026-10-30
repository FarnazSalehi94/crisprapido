# Profile, optimization and held-out validation

**Decision: the frozen preflight does not pass.** The optimization preserved
all measured output and improved a bounded genomic pilot, but the pilot does
not certify genome-wide sensitivity or production-scale resource use. The
frozen 100-guide × 12-assembly validation was not launched, and independent
review remains outstanding.

## Profile method and scope

Baseline is commit `a817ab5`, optimized development binary is built separately
from the code accompanying this report. A disposable source copy was
instrumented with nanosecond wall-clock guards and counters; `make_profile_build.py`
records how. These are **instrumented**, single-run measurements, not the
unmodified benchmark timings. Their PAFs were compared byte-for-byte to the
corresponding uninstrumented outputs. Timings are inclusive only of the named
code scopes, can overlap, and do not sum to total wall time. In particular:

- Columba is wall time waiting for the two child searches, including startup.
- Expansion is construction of the interval vector, not the WFA calls or all
  filtering/iteration that consume it.
- WFA is the `align` call only; alternative enumeration is the eligibility-aware
  fallback called after the selected WFA CIGAR fails the limits.
- PAM time includes adjacent-PAM extraction and the timed primary/reporting
  eligibility checks; inline comparisons and unrelated overhead are not fully
  assigned. Selection measures the reporter's anchor-map construction and
  includes its PAM checks. CFD measures the scoring call, not output formatting.

The larger pilot contains two previously unused 5 Mb CHM13 regions,
`chr1:100000001-105000000` and `chr9:30000001-35000000`, indexed together,
and the first ten frozen guides. It includes repetitive and nonrepetitive
sequence; it is not a random sample. The following baseline profile is one
instrumented run. Counts are calls, not distinct biological loci:

| Stage | Calls | Time |
|---|---:|---:|
| Columba candidate discovery | 2 | 1.683 s |
| Candidate interval-vector expansion | 14,807 | 0.0138 s |
| WFA2 alignment | 1,632,892 | 14.183 s |
| Alternative alignment enumeration | 1,628,157 | 5.021 s |
| Timed PAM extraction/checks | 19,047 | 0.00116 s |
| CFD scoring | 185 | 0.00047 s |
| Representative selection | 10 | 0.00085 s |
| Total process wall | 1 | 24.003 s |

These measurements expose the dominant alignment work; the small vector-build
number is **not** the total cost attributable to expanded candidates. Columba
returned 14,807 e=5 mapped candidate rows, zero e=0 rows, and 14,807 unique
raw candidate keys. The profile counted 1,632,892 WFA calls, including primary
and expanded intervals, and the run reported 185 distinct PAM anchors. The
actual e=5/e=0 SAMs and 10 Mb FASTA/index remain outside Git. `scaling/manifest.json`
records paths, region identities, hashes and counts; timing/profile data and
the 185-row PAF/audit summary are retained here.

Separate held-out synthetic subpanels were profiled to expose sequence
composition. The repetitive subpanel had 45 interval expansions, 2,706 WFA
calls and 2,249 fallback calls: baseline WFA 17.93 ms, fallback 6.54 ms,
Columba 27.45 ms and total 56.38 ms. The nonrepetitive subpanel had 86
expansions, 3,965 WFA calls and 3,776 fallback calls: WFA 29.48 ms, fallback
11.74 ms, Columba 34.97 ms and total 81.24 ms. These are tiny and differ in
size, so they are not a direct per-base repeat/nonrepeat performance comparison.
The `repetitive/` and `nonrepetitive/` profile TSVs also retain PAM, CFD and
selection counts/times.

## Optimization and rejected cache

A per-guide cache of previously verified reference intervals was tested first
and then **removed**. On the pilot it produced identical output but hit only
3,082 of 1,632,892 lookups (0.19%); in one run baseline/cache elapsed times
were 22.39/24.05 s and peak RSS 89,312/175,592 KiB. On the small repetitive
subpanel it hit 172/2,706 calls and on the nonrepetitive subpanel 0/3,965.
The cache passed the then-existing Rust and Python regressions, five historical
panel byte comparisons, and independent audits; its memory/runtime tradeoff on
this pilot did not justify keeping it. Its profiling TSVs are preserved, but
no cache is present in the final implementation.

The retained optimization instead calculates match and mismatch accounting
for each candidate single-gap position before allocating and reconstructing a
CIGAR. It constructs a CIGAR only for an eligible alignment that can equal or
improve the current score. It preserves the same guide-N accounting, f32
fraction comparison, affine scoring and lexical CIGAR tie break; WFA2, PAM,
thresholds, expansion geometry, CFD and reporting are unchanged. A reference
implementation copied from `a817ab5` is used **only in a unit test** to compare
2,000 deterministic resolved/ambiguous random and near-match sequence pairs,
three mismatch limits and four fractions including NaN. Full returned values,
including CIGAR and score, were identical.

In the instrumented 10 Mb run, fallback enumeration fell from 5.021 s to
1.947 s with the same 1,628,157 calls; WFA remained 1,632,892 calls and took
14.401 s in the optimized run. Instrumented totals were 24.003 vs 19.201 s;
these include timer overhead and are single samples. Uninstrumented paired
runs are reported below.

## Correctness and held-out panels

After the cache experiment and after the retained optimization, the full
existing Rust and Python regressions were run. Each candidate implementation
was compared against `a817ab5` on the five historical indexed panels and
all resulting 297 PAF rows passed the independent emitted-record auditor.
The final retained implementation produced **byte-identical PAFs** to
`a817ab5` for those panels and for both new held-out panels. It also produced
byte-identical PAF to baseline on the 10 Mb pilot, whose 185 rows passed the
independent auditor. The pilot had 78 plus and 107 minus outputs, including
146 I-containing and 35 D-containing records; all 185 had mismatches. These
categories may overlap.

The held-out synthetic panel uses seed 20261102, three guides, new random
flanks and independently enumerated exact, mismatch, I, D and combined cases
at multiple positions on both strands; one guide is repetitive. It contains
126 small contigs. The held-out real panel consists of four unchanged 1 kb
windows from CHM13 and GRCh38 at `500001-501000` and `900001-901000` on the
first two FAI contigs; these windows were not used in the correction panels.
It uses 16 deliberately derived guide designs, so it is not an unbiased natural
sample. The same independent bounded oracle and exact PAM-position key used in
the earlier validation were applied.

| Held-out panel | Oracle PAM anchors | Recovered | Precision | Strand evidence |
|---|---:|---:|---:|---|
| Synthetic | 223 | 223/223 | 223/223 | 113 plus, 110 minus |
| Real windows | 18 | 18/18 | 18/18 | 17 plus, 1 minus |

There were no missing or unexpected PAM keys to investigate. The strict
interval/category comparison remains separate: synthetic exact 6/6, mismatch
46/72, I 41/48, D 30/45, mismatch+bulge 100/351; real exact 5/5, mismatch
0/5, I 0/5, D 6/11, mismatch+bulge 7/32. Category denominators overlap.
The nonidentical interval/category cases have a reported hit at the exact same
PAM key; this does **not** assert equivalent edit placement or CFD. Raw
SAM/PAF, oracle truth, metrics, small fixture FASTAs, provenance and audit
summaries are under `heldout_synthetic/` and `heldout_real/`.

## Runtime, memory and scaling estimate

Two alternating paired runs on the 10 Mb pilot produced median wall times of
22.51 s at `a817ab5` and 19.28 s optimized, about **14.4% lower** for this
sample. All four PAFs were byte-identical. Peak RSS ranged 89,348–89,632 KiB
baseline and 89,476–91,380 KiB optimized. Only two samples per binary were
collected, and concurrent host activity can affect them. Larger reference
loading and a whole-genome index will have very different memory behavior.

As a deliberately weak capacity illustration, scaling the optimized 19.28 s
linearly from 10 Mb × 10 guides to the frozen 36.55 Gb × 100 guides gives
about **196 serial wall-clock hours** (about 229 hours for the baseline) before
whole-assembly index loading, audits, preflight and scheduling. This is **not
a scheduling estimate**: the two selected regions are not representative,
repeat and candidate densities vary greatly, cache/index I/O and threading are
not linear, and the frozen task uses different guide/assembly distributions.
A larger stratified, resource-bounded pilot would be needed to establish a
usable range. No full validation job was submitted.

All 157 standard Rust tests and the normally ignored pinned standalone CFD
test passed. All 48 Python tests passed, including held-out oracle and
provenance checks. The baseline development binary SHA256 is
`b4cd93f57f512a4727a3561b6c13289f03646b58964117f13058c1884c8c8b45`; the
optimized development binary is
`c37f4ded6e7c6afb3a7d22eccc5566f291244731bf6872dbfbd0b4d47f53aa74`.
The original production-candidate binary remains
`f59b1bbf58e70b321359c7bd7b8687aad81bb874b22dc2f083319453726bc876`.
Large indexes, FASTA and raw SAMs remain outside Git at
`/export/local/home/salehi/crisprapido-scaling-pilot-20261031/`;
held-out full verbose reports remain under the corresponding heldout directory.
No independent human review, production-binary change, full search, merge or
CRISPR-Bulge integration occurred.
