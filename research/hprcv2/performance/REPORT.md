# CRISPRapido performance investigation — 2026-11-11

Bounded, guarded investigation. **No full Brunello or SBDS HPRCv2 rerun was
launched.** External root:
`/moosefs/salehi/HPRCv2_CRISPRapido/runs/performance_investigation_v1/`.
Machine-readable results: `04_reports/performance-report.json` (copy committed
under `research/hprcv2/performance/results/`).

## Task 1 — What the pilot jobs actually did

Slurm 2871689/2871690/2871691/2871692 each allocated **8 CPUs; the application
used exactly 1**. The experimental PAM-anchored backend is single-threaded by
design and rejects `-t != 1`. The native search is not multithreaded.

Stage walls (seconds; one fresh call each, quiet node, warm caches):

| Job | Task | Search wall | Native traversal | Independent audit | Neural annotation | Slurm elapsed | Peak tree RSS |
|---|---|---:|---:|---:|---:|---:|---:|
| 2871689 | SBDS m5/b1/z2 | 419.6 | 101.03 | 4.0 | 17.0 | 450 | 14.30 GiB |
| 2871691 | SBDS m2/b1/z1 | 317.5 | 0.196 | 1.0 | 10.0 | 334 | 14.26 GiB |
| 2871690 | Brunello A m4/b1/z1 | 335.5 | 22.97 | 2.0 | 19.0 | 363 | 14.31 GiB |
| 2871692 | Brunello B m5/b1/z2 | 575.9 | 262.80 | 7.0 | 25.0 | 614 | 14.32 GiB |

**The reported search wall DOES include setup, validation and CFD/output;
neural scoring and the independent audit are separate stages.** The fixed
overhead inside one search invocation was instrumented on this same assembly in
the attested optimized pilot (job 2871283, same code path):

| Stage inside one invocation | Seconds |
|---|---:|
| Private index+reference snapshot copy to scratch | 99.9 |
| Full snapshot verification before search | 77.7 |
| Full snapshot verification after search | 77.5 |
| Native index binding check + index load + reference load | 82.5+3.3+18.3 |
| Native traversal (m5/b1/z2, one guide) | 80.6 |
| Wrapper target extraction (re-reads whole reference) | 27.6 |
| Rust startup/parse/CFD/format/output residual (21,447 rows) | 3.5 |

So the "7 minutes 30 seconds for ONE guide" is **~310 s of fixed per-invocation
overhead + 101 s of m5/b1/z2 traversal + ~20 s of separate audit/scoring** —
not a per-guide search cost. At matched historical settings (m2/b1/z1) the same
job's traversal was **0.196 s**.

## Task 2 — Matched old-vs-new benchmark (HG00735_pat, identical guides)

Identical reference, identical first 1/100/1000 frozen Brunello guides,
identical `m2/b1/z1 -f 0.75 -p GG`, sequential isolated jobs on octopus02.
New engine: 8 CPUs allocated (1 used), strict NGG, full independent audit.
Old engine: preserved production binary `crisprapido_newCFD` (SHA256 verified
against the production pin), production `-t 16` and Guix runtime, index staged
to local scratch exactly as in the preserved full-library run; staged files
checksum-verified.

| Guides | New wall | New traversal | New rows (all audit-PASS) | Old wall | Old CPU% | Old rows | Old strict-eligible witnesses |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 309.5 | 0.163 | 20 | 11.1 | 99 | 722 | 15 |
| 100 | 328.5 | 17.10 | 3,197 | 37.5 | 85 | 63,506 | 2,839 |
| 1,000 | 489.8 | 173.40 | 33,094 | 383.2 | 79 | 575,350 | 28,750 |

**Detected hits vs validated eligible anchors:** the old pipeline emits ~575
rows per guide, but only ~5% are strict-eligible NGG/m2/b1/z1 witnesses; the
rest are the known legacy non-NGG/geometry rows. The new engine emits ~33 rows
per guide, every one independently audited eligible.

**Containment: all 28,750 old strict-eligible witnesses are present in the new
output; zero missing** (15/15, 2,839/2,839, 28,750/28,750). This is bounded
emitted-witness containment at matched settings on one assembly — not a
genome-wide recall proof — but it shows the corrected backend loses nothing the
old pipeline validly reported here.

## Task 3 — Scaling diagnosis

- **Indexes are NOT loaded per guide.** One invocation loads the index once and
  shares it across all guides: the 1000-guide call costs 489.8 s total, of which
  173.4 s is traversal — marginal wall 0.179 s/guide, of which traversal is
  0.174 s/guide. Fixed cost per invocation, not per guide.
- **Full reference preparation IS repeated per invocation** (snapshot copy
  ~100 s, pre+post verification ~155 s, wrapper re-read ~28 s). This is
  mandatory binding/integrity work, not waste; it amortizes to ~5 ms/guide at
  full-library scale. The neural stage additionally does its own pinned full
  reference pre/post hashing (also amortizable).
- **Guides are searched sequentially on one CPU** — yes, by design; the backend
  rejects `-t != 1`.
- **The 8 allocated CPUs are underutilized** — yes: 1 of 8 (12.5%).
- **Search-state expansion dominates only in the expanded domain:** m5/b1/z2
  traversal is ~101 s/guide (SBDS) vs 0.196 s at m2/b1/z1 — ~514× more states
  (1.00e9 vs 2.83e7). At matched settings traversal is 0.174 s/guide and fixed
  overhead dominates small batches.
- **CRISPR-Bulge is not a bottleneck:** separate stage, 10–25 s per 38k–70k
  rows, models load once per invocation; it scales with rows, not guides.
- **Batching works:** 309.5 s/guide at N=1 → 0.49 s/guide at N=1000 →
  marginal 0.179 s/guide. All guides of an assembly should be in ONE invocation.

Throughput after amortization: **~20,100 guides/hour per assembly on one core**
(new engine, matched settings). The old pipeline's preserved full-library run
achieved **48,275 guides/hour per assembly at ~1.9 cores average** (its `-t 16`
never exceeded ~2 effective cores on this workload: 99%/85%/79% of one CPU at
N=1/100/1000).

## Task 4 — Recommendations (no code changed without authorization)

1. **Batch all guides per assembly into one invocation** (already supported;
   no change needed). This alone turns the alarming 7.5-min/guide into
   ~0.18 s/guide at matched settings.
2. **Parallelize across assemblies, not within**: one 1-CPU Slurm job per
   assembly, as the historical array did. This is the intended scaling
   dimension; 464 independent jobs.
3. **Do not shard guides within an assembly** to use idle cores: each shard
   would re-pay ~310 s snapshot/verification and re-load the ~11.5 GB index.
   Only worth considering if single-assembly wall time matters, and requires
   explicit authorization.
4. **Keep all verification/hashing as is.** At full-library batch scale it is
   negligible; it is what makes results trustworthy.
5. **Keep neural scoring a separate stage** (it already is) and out of search
   timing.
6. **The real cost is the requested m5/b1/z2 domain, not the implementation.**
   At m5/b1/z2, traversal alone is ~87–101 s/guide (measured on 4 guides across
   pilots). For 77,441 guides that is ~1,900–2,200 core-hours **per assembly**
   (~78–90 core-days), i.e. ~100 core-years for 464 assemblies — infeasible as
   currently designed without authorized algorithmic change (parallelism,
   pruning, or domain reduction are all out of scope without discussion).
   At matched historical m2/b1/z1 the corrected backend needs only ~3.9 h per
   assembly on one core.
7. **Validate the extrapolation before any large launch:** one full-library
   single-assembly timing run (~4 h, one job) would replace the linear estimate.
   This requires explicit approval; it was not launched here.

## Task 5 — Comparison table

| Metric | Historical full library (old, preserved) | New engine, matched m2/b1/z1 (benchmark) | Old engine, matched (benchmark) |
|---|---:|---:|---:|
| Full-library wall per assembly (77,441 guides) | **5,775 s measured** (tux09, ~1.9 cores) | 14,188 s extrapolated (3.94 h, 1 core, unvalidated at scale) | 29,746 s extrapolated — **5.15× overestimate vs 5,775 s measured** |
| Throughput, guides/hour/assembly | 48,275 | ~20,100 (marginal, single core) | n/a at scale (measured value above) |
| CPU utilization | 190% of one CPU despite `-t 16` | 1 of 8 allocated (12.5%) | 99%/85%/79% of one CPU at N=1/100/1000 |
| Search-only per guide | not separable in preserved logs | 0.174 s traversal | 0.31–0.39 s CPU |
| Initialization overhead per invocation | ~14 s staging (scratch, warm) | ~310 s fixed (copy+verify+loads) | ~10 s (index from local scratch) |
| Rows per guide | ~696 (only ~5% strict-eligible) | 33 (100% independently audited) | 575 (5% strict-eligible) |
| Eligible-anchor containment | — | contains all 28,750 old strict-eligible witnesses, 0 missing | — |
| Estimated full Brunello per assembly, requested B domain (m5/b1/z2) | not run historically | **~1,883 core-hours (~78 core-days); infeasible without authorized change** | — |
| Estimated full-HPRCv2 (464 assemblies), matched settings | 741 h wall (measured basis) | **~1,829 CPU-h search** + ~62 h one-time index import/verify; e.g. ~2 days wall at 40-way parallelism | 3,834 h by (invalid) linear extrapolation |

**Bottom line:** the apparent 7.5-min/guide discrepancy is per-invocation
overhead (~310 s) that batching eliminates, plus a genuinely expensive requested
search domain (m5/b1/z2 ~500× matched m2/b1/z1 in traversal states). The
corrected backend at matched settings is within ~1.3× of the old pipeline per
core-hour while producing only valid strict-NGG rows and containing every valid
witness the old pipeline emitted. Full-library m5/b1/z2 Brunello across HPRCv2
remains computationally infeasible without authorized algorithmic or domain
changes.

## Provenance and limitations

Jobs 2871709–2871711 (new) and 2871718–2871720 (old) completed exit 0 under the
memory/disk supervisor; plans v1/v3 are checksum-pinned. Three earlier chain
attempts failed on tooling bugs (staged executable bit; validator flag
spelling) before producing benchmark data; their incomplete attempts
(2871712/2871715) are preserved, and the fixes are unit-tested.

- Guides are the first 1/100/1000 records of the frozen library in file order —
  deterministic but not a random sample; guide-mix effects on traversal cost
  are unmeasured.
- One assembly, one fresh call per point, warm caches, quiet node — not
  repeated medians.
- Old-engine small-N extrapolation is demonstrably invalid (5.15×); the
  new-engine estimate has the same structural risk until a full-library
  validation run is approved.
- Old-engine outputs are diagnostics of a preserved binary with known legacy
  PAM/geometry semantics; never valid search results.
- No search algorithm, defaults, representative policy, formulas, production
  binaries or historical artifacts were modified.
