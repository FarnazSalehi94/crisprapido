# HPRCv2: independent SBDS and Brunello projects

Development branch only. No search/scoring algorithm, defaults, representative
policy, production executable, historical result or main-branch changes.
Completed safety-pilot reports: [SBDS](SBDS_PILOT.md) and
[Brunello A/B, cost sensitivity and recommendation](BRUNELLO_PILOT.md).
Bounded old-vs-new performance investigation with matched-settings benchmark:
[performance/REPORT.md](performance/REPORT.md).
Both have external PDF/PNG pilot diagnostics. **The full 466-input SBDS rerun
and six population figures are not complete; full Brunello has not launched.**

Integrated CFD pruning is disabled. Native searches explicitly select
`--search-backend pam-anchored`; this route does **not** use WFA2 verification
(the development CLI still links its unchanged legacy backend).

## Initial audit — 2026-11-11

External independent roots:

- `/moosefs/salehi/HPRCv2_CRISPRapido/runs/sbds_pam_aware_v1/`
- `/moosefs/salehi/HPRCv2_CRISPRapido/runs/brunello_pam_aware_v1/`

The requested stage directories exist. Authoritative initial inventory is
`sbds_pam_aware_v1/00_manifest/audit-v6/`, completed in Slurm **2871688**.
Earlier audit development attempts remain incomplete, not authoritative.
Small summary: [results/input-audit.json](results/input-audit.json).

- Frozen manifest: **464 HPRCv2 haplotypes + GRCh38 + CHM13**. Brunello's
  requested cohort excludes both references. Alias mappings are explicit in
  `assemblies.tsv`: in particular HG002 PanSN names and long reference filenames
  differ from the human-readable analysis IDs. References are the frozen
  GRCh38 no-alt PanSN and CHM13v2 masked-Y/rCRS PanSN derivatives, not substitute
  assemblies.
- Brunello: **77,441 unique IDs, 77,441 unique resolved 20-nt sequences**;
  **19,114 gene symbols + the non-targeting-control category, 1,000 controls**.
  Frozen FASTA SHA256:
  `a4bac2dc4c8f5ad0c37713751835b416d1669f899de1cd05ca1a0e828be48971`.
  Raw/Unix source tables are identical after newline decoding; clean-table
  gene/sequence/PAM multisets match. The original FASTA's space-to-underscore
  control-label mapping is documented, not reapplied to change IDs. External
  release/version equivalence was not independently downloaded/certified.
  Provider reference: https://www.addgene.org/pooled-library/broadgpp-human-knockout-brunello/.
- Both historical search scripts use **m2/b1/z1/f0.75/GG**, not either new
  Brunello setting or the requested SBDS m5/b1/z2 domain. Historical lineage
  uses Columba candidate discovery/WFA2 verification and has known PAM/export
  issues; these are not equal-semantics baseline searches.
- Original SBDS: **464/466 PAFs fully reread**, 372,311 reported rows;
  459 assemblies with one `20=`, four with none, one with two. **HG03784_hap1
  and NA18983_hap1 PAF reads fail with EIO**. Their counts are unknown in this
  fresh read, never zero. Earlier complete analysis reported 373,901 rows and
  the 461/4/1 distribution. That historical report is not a new successful
  reread. Reported perfect-match representations are not yet duplication or
  homologous-locus claims.
- `Brunello_all_HPRCv2`: 466 directories, **no done markers**, not a completed
  run. `Brunello_all_HPRCv2_scored_unscored_v2`: **456/466** small-artifact
  checksum/arithmetic checks pass; **10 fail with EIO**. Large compressed PAFs
  were **not reread**. Historical completion QC previously reported 466
  successful tasks and 24,956,881,788 candidates, but only CFD>=0.20 plus
  1,308 unscored records were retained. Low-CFD rows cannot be recovered from
  that archive; it is not an unfiltered benchmark. New 456-set subtotals in
  the inventory must not be called full-cohort totals.
- Prior index recovery independently verified **12** sets and replaced damaged
  structures/text for nine. Use its resolved manifest, never those nine damaged
  original prefixes. Another **454** sets lack full-content verification in
  this evidence chain. Metadata presence is not a content/sequence-binding
  check. No automatic rebuild, source substitution or damaged-file repair was
  attempted in this continuation.

Sources: historical scripts under `$BASE/scripts/`; previous
`analysis/PAM_Index_Recovery_2026-10-29_v1/RECOVERY_REPORT.md`,
`analysis/Brunello_QC_PAM_audit_2026-10-20/FINAL_REPORT.md`, and
`analysis/SBDS/final_results/SBDS_results_summary.md`. Preserve the distinction
between prior successful QC and current storage readability.

## Guarded pilot plan, not full-run approval

Slurm capability job **2871687**, octopus02: `/export` is mounted from octopus01;
actual development CLI, embedded source helpers, bound HG00735 index and neural
interpreter are readable. Available memory was about 995 GiB; local scratch
about 504 GB. This observation is specific to this worker/time, not a promise
about all nodes. Initial free storage: about 53.0 TB MooseFS and 11.2 TB export.

Completed safety inputs are the already bound **HG00735_pat** assembly and:

| Project | Input | Configurations |
|---|---|---|
| SBDS | `ACTGAAATCTGTAAGCAGGC` | primary m5/b1/z2; separate m2/b1/z1 historical-domain control |
| Brunello | exact frozen QKI, STX12, TMEM145 IDs in inventory | A m4/b1/z1 versus B m5/b1/z2 |

All use f0.75/GG, one application thread, both scores, no CFD filtering and no
neural pruning. Three Brunello guides are deterministic GC q10/q50/q90
**safety/feasibility strata**, not a representative statistical sample of
77,441 guides or 464 assemblies. No full-cohort cost confidence interval or
final parameter recommendation can be justified by those three alone.

Resource envelope: eight CPUs reserved (one used by application), 64 GiB
requested, 48 GiB sampled process-tree stop, 64 GiB available-memory start
check, 8 GiB running available-memory floor, >=60 GB scratch before start,
5 GB running free-space floor, 10 GiB per-file limit, 6,600 s stage timeout,
two-hour Slurm tasks. After the dedicated worker capability check, two independent
serial pilot chains were submitted, at most two searches active on the same
otherwise-idle 1-TB worker (16/48 CPUs reserved; at most 96 GiB sampled RSS and
120 GB minimum scratch budget combined). A/B are single fresh-process calls,
not replicated benchmarks; concurrent neighbours and warm/uncontrolled caches
must not be hidden in comparisons.
The RSS sampler is **not** a kernel/cgroup hard limit and can miss spikes or
double-count shared pages. Slurm's historic CR_CPU scheduling did not reserve
physical memory; worker checks remain mandatory.

`make_plan.py` freezes separate tooling, pins code/binaries/matrices/input
manifests/models, and preserves the guide IDs. `worker.py` requires Slurm,
exclusive task locks, new attempt directories, native index binding, a verified
prepared FASTA offset manifest, and full existing backend verification. Search
and neural annotation are separate durable stages to avoid holding search
allocations while loading the ensemble. Each invocation loads five models once
across all three guides, rather than once per row. Source/HDF5 pins and prior
all-supported-row implementation-agreement evidence are retained; no training,
network downloads or inference pickle loading.

All emitted rows are independently checked using samtools/CIGAR/PAM replay and
CFD reconstruction. PAF identity, annotations, missing-score reasons and anchor
uniqueness are checked separately. A checksum-bound `SEARCH_COMPLETE.json`
allows **explicit** `--resume-search` after downstream failure without repeating
a completed genome search. It is not an eligibility PASS. Overall `SUCCESS.json`
is published only after every stage, checksum check and memory supervisor pass;
partial files alone are not success, missing counts are not zero. Re-entering a
completed task verifies hashes and skips it. No completed result is overwritten.

The full raw eligible PAF and both-score annotations remain authoritative.
Neural unavailability means absent numeric PAF tag and explicit TSV `N/A` plus
reason, never a zero or removed row. Analytical CFD views at none/0.05/0.20 are
**downstream only** and use independently reconstructed full precision rather
than rounded PAF scores. Unavailable CFD is counted/preserved separately;
neural unavailability never decides retention. Gzip-byte estimates are actual
stream compression at the same level, not row-count-only extrapolations.

## Expansion gates / work still pending

**Full SBDS is not yet launched:** finish the guarded pilot, resolve current
storage faults/copies and prepare verified bindings for the remaining assembly
inputs. Full verification alone requires reading multiple terabytes of index
components, plus frozen source FASTAs. Scope and concurrency must be budgeted;
stat checks must not silently authorize a 466-task search array.

Then compare old/new exact anchors (guide sequence, reference/contig, strand,
PAM-adjacent boundary), intervals, CIGARs and CFD. The historical-domain control
separates setting expansion from backend/reporting differences; valid alternative
representatives can have different scores. Assess all four no-perfect cases,
HG00321_hap2 and chr7/12/17 candidates using actual sequence and homologous
context. Do not infer a duplicated biological locus from two reported records.

The six final SBDS figures and group report require an explicit 466-input
completion/assessability ledger. Do not relabel an n=1 safety pilot as an updated
population figure. The historical five generated panels plus separately
hard-coded perfect-match plot must be regenerated from validated data, not copied
counts: CFD distribution, candidate reduction, per-assembly candidates, recurrent
high-CFD patterns, target-site variation, and perfect-match distribution. Add
explicit neural coverage/availability and old/new eligibility/locus supplements.

The prior independent SBDS PAM audit reports 26,245/373,901 canonical NGG rows,
including only 465/1,400 rows above CFD 0.20. Sequence follow-up reports observed
PAM **GGC** for the recurrent chr12 pattern and **AGT** for chr17; neither is NGG.
Those specific old anchors must not survive merely to preserve historical high
scores. This is an eligibility/reporting correction, not evidence that the
underlying homologous regions lack other valid NGG sites. Chr7 perfect sites
have GGG. These prior sequence/PAM audits are separate from our current 464-file
raw reread and do not establish Brunello's PAM proportions.

**Full Brunello remains prohibited until user review** of A/B pilot results,
expanded representative validation, cost/uncertainty, retention strategy and
parameter recommendation. Do not interpret lower A counts as an equivalent
faster B search: A and B have different biological search domains and expected
neural model coverage. No parameter default is changed. Before approval, plan
bounded guide batches (native process/index reuse), retention of every eligible
hit and score-status sidecars, and compression of large sequence/PAF outputs.

Later homologous-locus/allele/figure requirements are specified in
[ALLELE_PLAN.md](ALLELE_PLAN.md); IMPG 0.2.4 is installed but no mapping/full
library frequency analysis has been run.

## Tests

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s research/hprcv2 -p 'test_*.py' -v
```

**48 Python tests passed:** 15 new inventory/plan/accounting/collection/legacy
replay tests, four workflow auditors, nine CFD-model tests and 20 configured
scoring tests (none skipped). New tests include completed-search signature
identity, pending-not-zero accounting and unavailable-versus-actual-zero cases.
Independent live audits additionally checked all 121,959 emitted pilot rows and
95 stratified/tail official predictions. Slurm 2871695/2871696 check completed-task
idempotence without a new search. Existing core Rust evidence is reused because
core sources are unchanged; this continuation has not rebuilt production code.
