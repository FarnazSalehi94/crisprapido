# Brunello 10k pilot: scoring-only completion

## Submission record (2026-11-23)

**Scoring job 2880112 submitted; finalizer 2880113 requires afterok:2880112.**
At submission, the full 57,578,781-row scoring stage is not complete. No success
is claimed until full verification, successful supervisor and successful Slurm
scoring job. Authoritative eventual result: `SUCCESS.json` and `RESULTS.md` in:

`/moosefs/salehi/HPRCv2_CRISPRapido/brunello_10k_scoring_20261123_v1/`

This is only the existing **10,000-guide / HG00735_pat** pilot,
m4/b1/z1/f0.75/GG. **No SBDS, genome search, library expansion or cohort run.**
The full library remains 77,441 guides; the HPRC manifest is 464 haplotype
assemblies, not SBDS's 466 including references.

## Writable checkout and preservation

The login-host `/export` filesystem remains read-only. Existing MooseFS
checkouts were inspected, but were unrelated/stale or dirty and left untouched.
A new writable checkout was cloned from read-only HEAD `f1157eb` at:

`/moosefs/salehi/HPRCv2_CRISPRapido/brunello-writable-checkout`

Branch: `brunello/bounded-10k-scoring-20261123`; origin is the existing user's
`git@github.com:FarnazSalehi94/crisprapido.git`. Saved scoring/handoff patches
were checked and applied. The pre-existing 803-line handoff delta was separately
saved and carried into this checkout, not overwritten or silently committed as
new Brunello work. Original untracked outputs/SBDS scripts remain untouched in
the original checkout. Only explicit task source/documentation files are staged.

Previous attempts, especially 2876332, and old search/audit checkpoints are
never changed. The prior signed bounded-pilot evidence bundle remains intact:
`/moosefs/salehi/HPRCv2_CRISPRapido/brunello_resume_20261123_v1/REPORT.md`.
A continuation link, rather than rewriting that archived report, identifies
this new stage. Kernel OOM remains unconfirmed for 2876332: the recorded cause
was the sampled process-tree RSS guard, not an established unique-memory peak.

## Verification before launch

Reconfirmed the complete 101,000-row prior comparison: maximum numeric CB
absolute delta **1.51000000003787e-9**, within the predeclared **1e-6** tolerance.
Every non-CB PAF/annotation field, row coverage/order, CFD strings and CB
availability matched exactly. Both bundles' hashes and all three retained
checkpoint input hashes were verified. All real-data pilot rows were scored;
no N/A branch coverage is claimed for that sample.

Tests rerun in the writable checkout: **38 PASS** (14 bounded-wrapper,
12 completion-controller/validator, 11 adapter, 1 reference-index). Tests cover
cross-chunk duplicate rejection, changed-input hash rejection, corrupt
annotations, interrupted merge/resume, explicit N/A, row loss/extra rows,
CFD-format/availability drift, nonfinite CB, historical evidence hashes and
no-clobber publication. Shell syntax check passes.

The deployed bounded wrapper is byte-identical to the validated pilot wrapper;
unchanged pinned inner scorer, models, CFD, alignment/representative rules and
eligibility. Its global duplicate-anchor check uses disk-backed SQLite with
an 8-MiB page cache across the entire source, **before scoring**, not merely
within chunks. Checkpoint reuse requires exact chunk input SHA256, configuration
identity, COMPLETE status, output/annotation hashes and row counts. Input and
configuration content are checked again before merge publication.

## New stage identity and run

Frozen plan SHA256:
`d18707776e531285459c572a8a66a0a12c1938e153335c2c85ef47219cec1e6d`.

The migration explicitly retains old search/audit signature
`d5226c6c8002646581b135ebc4c1d13d9f5b9567dd26df6b5b5c4a57ca23fd9b`.
Worker preflight freshly verifies all original plan pins and all ten existing
search/audit output size/hash records; it never repeats either computation.
Source PAF: 9,473,582,518 bytes, SHA256
`fb09e1141f36e0f82561783dbacf6a302f92c2df547c2476a7ebe4edfb6da69c`.

Resources: **concurrency 1**, one CPU, 8G Slurm request, 4-GiB sampled summed-tree
RSS guard at configured 0.5 seconds, 36-hour cap, 50,000-row chunks (1,152 total;
last 28,781 rows), 200-GB running output free-space reserve. Per-PID RSS is
logged; summed RSS can double-count shared pages and is not unique memory or a
kernel-enforced peak. A prior approximately 51-second 50k chunk suggests roughly
16 hours of scoring alone, not a guarantee: full duplicate scan, I/O, repeated
reference hashing, merge and validation add work. No concurrency escalation.

Frozen scripts/plan live in the new evidence directory, independent of future
checkout edits. `research/hprcv2/brunello_expanded_benchmark/scoring_completion.py`
implements prepare/work/finalize and refuses scope/command drift.

Tracking:
- `attempt-2880112/PREFLIGHT.json`: verified checkpoint reuse and prior comparison.
- `attempt-2880112/progress.jsonl`: phase, elapsed time, completed chunks/rows,
  new per-chunk metrics, logical stage bytes and free bytes (10-second samples,
  storage rescanned at most once per minute). Logical bytes are not physical
  MooseFS allocation/replication.
- `scoring/metadata.json.incomplete`: current completed chunk records; safe
  readers must tolerate an in-progress JSON write.
- `scoring/chunks/*/chunk-config.json`: content-bound checkpoints. Failed or
  invalid chunks are retained with unique `.failed-*` names.
- `supervisor-2880112.log[.json]`, `score-2880112.log`, per-invocation logs and
  attempt `FAILED.json`: runtime/memory/failure evidence. No automatic retry
  loop or suppression of failures.

## Final validation and success gate

The controller independently streams all source PAF, merged PAF, annotations
and original 2876332 conversion spool in lockstep. The spool is complete
conversion/availability evidence, **not reused predictions**. Its archived
SHA256 is `b55abb94325d65af9eff565aea1d5ba2e626818e1115315e1579ef363f0b6615`
and size 29,881,742,568 bytes. Its full digest is rechecked during validation.

Require all 57,578,781 rows, exact raw-row prefix/order and CFD, all non-CB
annotation/context fields, exact original CB availability/N/A reasons, finite
numeric scores in [0,1], matching PAF/annotation CB strings, no extra/missing
rows, and matching source/merged/annotation/historical-spool hashes. Original
checkpoint metadata and all original input/tool/model pins are checked again.
Only then is an attempt CANDIDATE published.

The separate afterok finalizer requires successful scoring-job scheduler
state and supervisor child exit 0, verifies merged artifacts again, and only
then publishes **new-stage** SUCCESS/RESULTS. It does not overwrite old
`checkpoints/score.json` or old `SUCCESS.json`, or label a running wrapper as
complete. A failed/dependency-blocked finalizer is not success; query accounting
on login (compute-local sacct cannot reach its database). If scheduler state
has aged out before finalizer starts, finalization fails closed for manual
accounting reconciliation; it does not invent completion.

## Safe continuation

Inspect scheduler/progress and failures before any retry. To resume after a
failure, submit the frozen `tooling/submit.slurm ROOT work` under the same
resources; its new job creates a new controller attempt and reuses only
verified completed chunks. Submit a new afterok finalizer with that new scoring
job ID. Do not edit frozen plan/tooling, reuse partial predictions, delete old
attempts, or run search/SBDS/cohort. A completed wrapper after a controller
failure is independently revalidated rather than rescored.

**Publication job 2880114** requires afterok:2880113. Its hash-pinned publisher
checks scientific SUCCESS, expected local branch/HEAD, empty staged index and
expected origin before writing. It commits only generated actual results and
its own handoff section via a private index, leaving the pre-existing 803-line
handoff delta and unrelated working changes uncommitted. Push is fast-forward
only, with noninteractive SSH and normal strict host-key checking. Six Git
integration tests PASS (44 total relevant tests), including changed HEAD/branch,
remote or staged files and preservation of unrelated changes.

On success `PUBLICATION.json` records the pushed result commit; the new report
and handoff copies update to actual results. If Git/authentication/concurrency
checks fail, the publication log and FAILED marker require manual follow-up;
scientific SUCCESS is separate and must not be confused with successful Git
publication. The source repair/submission commit is already pushed, but no
full-10k completion/result commit is claimed while scoring runs. Reconcile
final detailed sacct accounting on login when available.
