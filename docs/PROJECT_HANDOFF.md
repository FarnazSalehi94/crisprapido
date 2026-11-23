# CRISPRapido project handoff

## 2026-11-23 Brunello 10k scoring-only completion submitted (NOT complete)

Writable checkout (original read-only checkout untouched):
`/moosefs/salehi/HPRCv2_CRISPRapido/brunello-writable-checkout`, branch
`brunello/bounded-10k-scoring-20261123`. Saved patches applied after preserving
the pre-existing 803-line handoff delta; that prior delta is not staged as new
Brunello work. Source/test/docs changes are confined to this task.

Full continuation: [BRUNELLO_10K_SCORING.md](BRUNELLO_10K_SCORING.md).
New evidence root:
`/moosefs/salehi/HPRCv2_CRISPRapido/brunello_10k_scoring_20261123_v1/`.
**2880112 scoring-only**, one CPU / concurrency 1 / 8G / 4-GiB sampled-tree
RSS guard / 50,000-row chunks / 36-hour cap. **2880113 finalizer** depends on
successful completion of 2880112. Both submitted; no full-10k success claimed.
Scope: existing 10,000-guide HG00735_pat pilot, 57,578,781 rows, m4/b1/z1/GG;
reuse original search/audit with fresh hashes, preserve all prior attempts.
Frozen plan SHA256 `d18707776e531285459c572a8a66a0a12c1938e153335c2c85ef47219cec1e6d`.

38 relevant tests PASS. Reconfirmed original 101k comparison and three input
checkpoint hashes: maximum CB delta 1.51e-9 <= predeclared 1e-6, exact coverage,
order, CFD and availability. Validated wrapper byte-identical to prior pilot;
full global disk-backed duplicate check precedes scoring, reuse binds input
SHA256 plus output/annotation hashes. Progress/memory/storage/failures recorded
under the new root. Independent final full-stream validation checks all source
rows, CFD, non-CB annotations and availability against the preserved complete
historical conversion spool (archived SHA256 required; predictions never reused).
New-stage SUCCESS requires validation + supervisor + successful Slurm scoring
job + afterok finalizer. Do not confuse wrapper COMPLETE or CANDIDATE with it.

No SBDS, new search, full library or full cohort launch. 2876332 remains an
RSS-guard failure, not confirmed kernel OOM. Original archived resume report
and signed evidence are unchanged; see its CONTINUATION.md for this stage.
**Publication 2880114** depends on afterok:2880113 and will commit/push actual
results and its own handoff section only after scientific SUCCESS. Hash-pinned
publisher uses a private index and branch/HEAD/staged-state/remote guards,
never commits the pre-existing 803-line delta, and never force pushes. Six Git
integration tests PASS (44 relevant tests total). Failure is visible and needs
manual follow-up; scientific completion and Git publication are separate.
Source/submission documentation already pushed; no actual full-10k result is
predeclared. Final detailed accounting remains queryable on login.

## 2026-11-23 Brunello reconciliation / bounded scoring pilot COMPLETE

User directs Brunello across HPRC; SBDS is complete enough and excluded.
No full-cohort launch is authorized. Evidence and fixed tooling:
`/moosefs/salehi/HPRCv2_CRISPRapido/brunello_resume_20261123_v1/`.
**Repository /export filesystem is read-only (EROFS)**. This is a reconciled
handoff COPY; original `docs/PROJECT_HANDOFF.md` could not be updated. Apply
`scoring-fix.patch` and the handoff patch only after home is writable again.
No remount, permission changes, relocation changes, or git commits were made.

Frozen full guide file: `guides/brunello_77441.fa` under the project root;
77,441 FASTA records / unique IDs / unique 20-nt ACGT sequences. Expected SHA256
`a4bac2dc4c8f5ad0c37713751835b416d1669f899de1cd05ca1a0e828be48971`.
The reusable pilot uses the first 10,000 records, HG00735_pat,
m4/b1/z1/f0.75/GG; it is NOT a 77,441-guide run.
Brunello `00_manifest/assemblies.tsv` is **464 HPRCv2 haplotype assemblies**,
not the SBDS 466 including two references. All 464 source FASTA paths exist.
Original manifest index-status fields are historical and must not be used as
current readiness evidence. `binding-paths.json` records all 464 committed
binding pointers and 5,568 current component paths/resolutions/sizes: 457
resolve under project `indices/`, seven under relocated index-recovery home.
Historical `binding.json:trusted=false` on imported inputs is not a current
readiness determination; never silently promote that flag. Cohort-wide path
checks are not a fresh multi-terabyte content hash or native-load validation.

2876332: Slurm FAILED 1:0, 06:35:49, octopus06, 128G; no MaxRSS supplied.
Exact recorded cause is sampled process-tree RSS guard, not confirmed kernel
OOM: 174,890,644 KiB summed RSS versus 115,343,360 KiB limit; child -9.
This can double-count shared pages; no unique-memory peak or cause of the
abrupt increase was established. Preserve all failed outputs and original
failure/evidence inventories. New small evidence copies/hashes and large-file
size/mtime inventory live here; large partial spools stay untouched in place.
No partial historical predictions are treated as completed scores.

Reconciliation job **2880107 COMPLETED 0:0 (0:35)** verified every frozen
plan pin and search/audit output size/hash; `reconciliation.json` is PASS.
The client-side compute mount reports read-write, but actual login edits fail
EROFS; the repository has not been altered by this task.

2874508 search and 57,578,781-row audit remain the reusable stages. Existing
search/audit signature must remain `d5226c6c8002646581b135ebc4c1d13d9f5b9567dd26df6b5b5c4a57ca23fd9b`.
No score checkpoint/SUCCESS exists. Earlier reconciliation checks 2880099 and
2880100 failed closed on incorrect assumptions (466 versus actual 464, and
historical trusted flag). 2880102 passed data checks but compute-local sacct
was unavailable; accounting is now captured on login. 2880106 reverified pins
and checkpoints but stopped on an existing evidence-copy directory; this
rerun/idempotence bug was fixed. All logs retained, original outputs untouched.

Fixed bounded wrapper is separate frozen tooling, not an environment override
of the pinned old timing plan. Default 50,000 rows, hard maximum 100,000;
content-bound chunk reuse, annotation hashes/completion checks, unique attempt
logs, preserved failed chunks, output-directory lock, pre/post content checks,
and disk-backed global duplicate-anchor rejection. The inner scorer, models,
CFD, representative choices, and eligibility formulas are unchanged. Fourteen
unit tests pass, including changed input, annotation corruption, incomplete
status, invocation collision, cross-chunk duplicates, and actual merge crash.

Pilot job 2880101 started after successful reconciliation. It reads
original search output only; selects first 1,025 rows + every 576th source row
+ first strand/edit stratum (100k–110k rows total), preserving source order.
Compares monolithic versus 50k-row scoring, deliberately interrupts after one
verified chunk, then verifies that resume reuses it. All PAF/annotation fields
must match except numeric CB within pre-existing 1e-6 absolute tolerance.
One CPU, 8G request, 4-GiB sampled-tree guard with per-PID RSS evidence, 1 hour.
**Pilot 2880101 COMPLETED 0:0 (6:20)**, validation/supervisor PASS:
101,000 rows / 9,995 guides, all 20 observed strand/edit strata, all scored;
1,660 ungapped, 77,832 target-gap/RNA-bulge, 21,508 guide-gap/DNA-bulge;
35,416 CFD-zero rows retained. No N/A examples occurred in this real-data
sample. Three chunks 50,000/50,000/1,000; interrupted after the first, which
was demonstrably reused on resume. All non-CB PAF/annotation fields match;
max numeric CB delta 1.51e-9 versus predeclared 1e-6 tolerance. Sampled summed
tree peak 890,124 KiB (0.848888 GiB), not a unique physical-memory measurement.
No full 10k scoring retry/full-cohort or SBDS launch. See REPORT.md and
pilot-validation.json. This is not exhaustive per-guide/cohort validation.
Full 10k scoring restart still requires a separate validated stage identity /
migration record; original search/audit need not be rerun. Do not launch full
library/cohort or reuse old unvalidated override paths.

## 2026-11-17 outcome reconciliation (supersedes earlier readiness claims)

Main model verified: `openai-codex/gpt-6-astra`; starting commit `c00c109`.
Actual accounting, not queue absence:
- Search 2875379: 48/48 COMPLETED, exit 0. Duplicate 2875384 was cancelled,
  not successful work.
- Wave 5 (2875393–95) and wave 6 (2875917–19): each build/recheck/promote
  stage 48/48 COMPLETED. All 96 promoted indexes and bindings are present.
- Wave-4 independent validation 2875903: 48/48 COMPLETED and 48 PASS JSONs.
- GRCh38 conversion 2875898 COMPLETED. Build 2875899 FAILED after 14:38:
  `source_fasta` was a Path and `json.dumps(manifest)` raised TypeError.
  Index construction and staged write-back had finished; scratch was retained.
- Brunello 2875915 FAILED after 1:41, before scoring: missing override digest
  in the submitted script. **Slurm snapshots scripts on submission**; editing
  the script file while pending did NOT update that submitted script. Also,
  the override code never assigned the replacement command. Prior claims
  that all reviewer issues were resolved and deployment was safe were wrong.

Snapshot `00_manifest/outcome-snapshot-20261117T151705Z.json`:
**335 unique COMPLETE / 466**, 131 remaining: 96 newly recovered + GRCh38 +
**34 EIO** (not 32). 432 READY_BOUND includes the then-old GRCh38 binding.
Counts: 64,118,953 hits/CFD-defined; 14,398,176 CB-scored;
49,720,777 explicit CB N/A. This snapshot checks SUCCESS-selected summaries,
not every output hash. Full output rehash job **2876348** writes
`00_manifest/outcome-review-20261117/verified-completions.json` on success.
That directory also preserves accounting and the two failure logs.
Initial rehash submission 2876340 failed immediately because Slurm's
`--wrap` used /bin/sh, which rejects `set -o pipefail`; retry 2876348
explicitly invokes bash. No output was produced by the failed attempt.
Current free space at snapshot: 27.78 TB decimal, not the older 34.7 TB.

### Work submitted and verified this session
- **2876318**: search/scoring for exactly 96 newly recovered assemblies,
  concurrency 4, octopus08; completed assemblies and GRCh38 excluded.
  Frozen plan/list: `sbds-plan-recovery56-20261117.json` and
  `sbds-assemblies-recovery56-20261117.txt` under `00_manifest/`.
- **2876323 → 2876324 → 2876325**: final 34 EIO rebuild/recheck/promote,
  concurrency 2, frozen `recovery-wave7-tooling-v1/`. No overlap with 96 searches.
- **2876334 / 2876335**: wave-5 / wave-6 independent indexed-text validation,
  48 tasks each, concurrency 2; separate from prior cross-node checksum checks.
- GRCh38 **2876328 → 2876329 → 2876330 all COMPLETED exit 0**.
  New `grch38-manifest-recovery-v1/` tooling rehashes the derived gzip,
  decompressed source, original scratch reference, and every scratch/staged
  component before recovering the manifest; no conversion or index rebuild
  repeated. Cross-node recheck on octopus07 reported TRUSTED_CROSS_NODE.
  Promotion/rebind completed; old references and failure evidence retained.
  Blocked jobs 2875900/01 cancelled. Rollback now renames partial bindings
  aside before restoring the old binding, rather than nesting directories.
  **2876345** submitted for GRCh38 search only, excluding active 2876318 list.
- Brunello **2876332** RUNNING: scoring-only retry using the **unchanged v3
  worker and pinned monolithic scorer**, 128G request / 110GiB sampled guard,
  octopus06. Frozen `timing-10000-execution-v5/`; plan SHA remains
  `d5226c6c8002646581b135ebc4c1d13d9f5b9567dd26df6b5b5c4a57ca23fd9b`.
  Search/audit checkpoints must exist and pass their hash validation; only
  `score-stage` exists in the new attempt (no repeated search/audit).
  This retains the scorer's global checks, but may still exceed its memory
  envelope. Do not cancel healthy work. No final scoring runtime, coverage,
  checksum PASS, or peak-memory result is available yet.

### Corrections and remaining gates
The chunked scorer is a **prototype**, not production-ready: config identity
uses paths rather than content hashes for several scientific inputs; global
cross-chunk duplicate validation, locking, failed-chunk invocation-log reuse,
annotation integrity, and crash-safe multi-file publication remain incomplete.
The repository timing worker now rejects environment scoring overrides until
an explicit validated checkpoint migration exists. The running v5 uses the
unchanged original worker and no override. The 500k-row probe is not a
measurement of 2M-row chunks, nor evidence of full-workflow CPU cost.
GRCh38's historical transformation JSON also overclaims “all hits preserved”
and its `bases` field counts only the final contig. Coordinate-preserving
IUPAC→N is a derived-reference policy, not proof of equivalence for hits
spanning ambiguous positions. Original evidence has not been rewritten.

Tests: SBDS suite **36 PASS**, timing-worker suite **11 PASS**. New recovery
verification test rejects changed source, scratch reference, or index bytes.
Next: monitor active searches; after wave-7 verified promotion, submit only
its 34 still-incomplete assemblies with active-list exclusions. Generate six
SBDS figures only after 466 COMPLETE. Full Brunello single-assembly search
still needs validated guide coverage/order/merge and restart safeguards;
IMPG allele-frequency/Figure 1C information must remain preserved. No
77,441-guide or 464-assembly cohort run launched.


## Session 2026-11-16 (night, addendum): all reviewer findings fixed pre-launch

Independent review of the chunked scorer (4 findings: 3×P1, 1×P2) — all
fixed, unit-tested (8 chunked + 10 worker tests green), and the **still-
PENDING** 2875915's frozen v4 tooling was refreshed before start (state
verified PENDING before and after; a mid-update start would fail closed on
the tooling checksum):
1. **Chunk reuse now binds a config identity** (`chunk-config.json` per
   chunk: hash of inner-scorer digest, scoring mode, reference/guides/model
   inputs, batch size, chunk size, flags). Stale chunks with any changed
   setting are rescored, never reused. Regression test covers config drift.
2. **Merge-crash recovery**: stale wrapper-owned `*.incomplete` merge
   artifacts are removed on resume; a crash mid-merge no longer wedges the
   out-dir.
3. **Worker override pinned**: `BRUNELLO_SCORING_SCRIPT_SHA256` must match
   the override's digest or the run refuses to start — the comment's
   "same pinned inner scorer" constraint is now enforced, not assumed.
4. **Annotations publish is conditional** (matches inner scorer contract).
Reviewer's memory note accepted: 2M-row chunks project ~6.4 GiB (linear
scaling), still under the 12 GiB stop.

Full-library readiness review (scout) confirms: not ready for cohort launch
(as planned). Arithmetic verified: search 163.8 h, audit 11.6 h, chunked
scoring 78.3 h ≈ **~254 h/assembly; ~118,000 CPU-h cohort; ~19.1–19.4 TB**
at measured 5,758 rows/guide (207.3B rows). Storage remainder vs current
free is tight (~15 TB after cohort data, before spools/staging) — a formal
reservation remains a launch prerequisite. The full-library worker still
needs guide-chunk search coverage/order/merge gates beyond the 10k worker.


## Session 2026-11-16 (night): chunked 10k scoring restart; waves 5–6 queued

**Brunello 10k scoring restart switched to the validated chunked scorer.**
The pending 128 G single-pass retry **2875902 was cancelled while still
PENDING** (never started, nothing lost): single-pass memory is unbounded by
design and already failed once at 48 GiB, while a 16 G job schedules far
sooner on the busy cluster. Measured scaling probe (1M rows, 2×500k chunks,
`chunked-scoring-probe-v1/`): **1.59 GiB peak per 500k-row chunk, 1,582 rows/s**
→ full 57.6 M-row scoring ≈ **10.2 h**. Restart **2875915**
(`timing-10000-execution-v4/`): byte-identical plan (search 76,267 s + audit
57,578,781 PASS checkpoints REUSED — only scoring runs), `BRUNELLO_SCORING_SCRIPT`
env override in timing_worker (narrow; actual command recorded in metrics),
`CHUNKED_SCORING_ROWS=2,000,000` (~29 chunks, ~5 GiB/chunk projected),
16 G request / 12 GiB sampled-RSS stop, octopus06. Inner scorer is the pinned
scorer-v2 `score_paf.py` byte-identical copy (SHA `51cde989…` verified).

**Recovery waves 5–6 + GRCh38 all queued** (cluster fully allocated by other
users; jobs advance as capacity frees): wave-5 builds running (4 done, 2
active, %2 on octopus10, array 2875393 → recheck 2875394 → promote 2875395);
**wave-6 chain 2875917 → 2875918 → 2875919** prepared and queued (48 disjoint
of the 80 EIO remaining after wave-5; exclusions wave1–6 all frozen);
GRCh38 chain **2875898 → 2875899 → 2875900 → 2875901**; wave-4 independent
validation **2875903** (48 tasks). After waves 5–6 promote: **32 EIO remain**
(wave-7) plus GRCh38 and HG06807 (in wave 5). Ledger: **SBDS 335/466**.

Chunked-scoring summary for the record: wrapper + equivalence validator
committed (`0393161`); real-data equivalence PASS (max cb Δ 2.2e-16,
batch-ULP only); global annotation row_index restored at merge; INCOMPLETE-
resume reuses verified chunks. Independent reviews of the chunked implementation
and full-library readiness are in flight.

Next: as promotions land — wave-7 for the last ~32 EIO, search waves with
frozen-list exclusions, GRCh38 search after rebind; 2875915 completes the 10k
benchmark (search/audit/scoring walls, CB availability, peak RSS, compressed
sizes) and finalizes the full-library single-assembly plan (~164 h search +
~11.6 h audit + ~79 h chunked scoring ≈ 255 h/assembly).


## Session 2026-11-16 (late): SBDS 335/466; chunked scoring implemented and validated

**Memory-safe chunked scoring is implemented, unit-tested, and validated
against the single-pass reference on real data**
(`research/crispr_bulge/score_paf_chunked.py` + `compare_scored_pafs.py`;
validation evidence: `sbds_ngg_nrr_v1/00_manifest/chunked-scoring-validation-v1/`).
Design: drop-in replacement with the exact score_paf.py CLI contract; splits
the input PAF into `--chunk-rows` chunks, scores each with the UNMODIFIED
pinned scorer (no model/CFD/representative/eligibility change), merges with
row-for-row input verification (stronger than the inner scorer's chunk-local
duplicate check), restores GLOBAL annotation row_index at merge, publishes
only after verification, and cleans derived chunk scratch only afterwards.
Restart-safe: `INCOMPLETE` out-dirs resume; verified chunks are reused.
**Validation on HG00706_mat NGG (38,085 rows, 4 chunks, pinned scorer-v2 +
models): PAF equivalence PASS — 0 differences beyond tolerance; 40 rows show
cb:f ULP differences from TF batch-composition effects, max |Δ| = 2.2e-16 at
scores < 1e-14; annotations identical modulo the same 40 ULP rows; peak
child RSS ~526 MB/chunk; wall 43 s.** Cause of the single-pass blowup
confirmed: the inner scorer holds a per-input duplicate-key set scaling with
rows (57.6 M rows → >48 GiB); chunking bounds it. This is the scoring design
for both the 10k restart (if 2875902's 110 GiB also trips) and the full
77,441-guide library (446 M rows → ~90 chunks).

**Brunello:** search (76,267 s) and audit (57,578,781 rows PASS) remain
checkpointed from 2874508; single-pass scoring restart **2875902** (128 G /
110 GiB stop) still queued behind cluster priority. If it trips again, the
restart uses chunked scoring with the same byte-identical plan. **Full-library
timing estimate (measured 10k basis): search ~164 h/assembly (~76,000 CPU-h
cohort), audit 11.6 h/assembly; scoring ~86 h/assembly extrapolated from 1k,
now chunk-safe; total ~262 h/assembly, ~121,000 CPU-h cohort.** Storage
~19.4 TB compressed unchanged.

**SBDS: 335/466 COMPLETE** (64.1 M hits, all CFD, 14.4 M CB-scored, 49.7 M
explicit N/A). Remaining 131 = 128 EIO + GRCh38 + HG06807 pair. Active:
wave-5 builds 2875393 (2 tasks running on octopus10, rest queued — cluster
busy with other users); wave-4 independent validation 2875903 queued; GRCh38
chain 2875898–78 requeued after the two diagnosed failures (cp-stdin bug,
then FileExistsError from the preserved-evidence dir; converter now takes
`--out-dir`, deterministic gzip, 35 tests). Brunello 2875902 queued.

Next: continue waves 5–6 to finish all indexes; GRCh38 search after rebind;
complete 10k scoring (2875902 or chunked) for the measured scoring time;
then the restartable full-library single-assembly validation run; SBDS final
figures at 466/466.


## Session 2026-11-16 (evening): SBDS 335/466; Brunello search+audit measured, scoring restart

**Brunello 10k job 2874508: search and audit COMPLETE and checkpointed;
scoring killed by the memory guard.** Measured (HG00735_pat, m4/b1/z1/GG,
10,000 first-in-file guides): search wall **76,266.6 s (21.2 h, 7.63 s/guide,
972.79e9 states, 76,071 s user CPU)**, raw PAF **9.47 GB**; independent audit
**PASS, 0 errors, 57,578,781 rows verified in 5,400.6 s** (5,758 rows/guide;
56.65 M bulged / 0.93 M ungapped; 1.47 M CFD≥0.02; strands ≈ balanced).
The both-score stage was killed at 94,015 s when sampled tree RSS hit
**50,333,920 KiB (>48 GiB stop)** — the scorer's memory grows with input rows,
so the 1k-derived 48 GiB envelope was insufficient at 57.6 M rows.
**Implication for full-library design: single-pass scoring is infeasible at
446 M rows/assembly; the full-library worker must score in verified chunks
(coverage/order/merge gates), exactly as the wrapper review required.**
Restart **2875902** (execution-v3) reuses the byte-identical plan (SHA
`d5226c6c…23fd9b`) so checksummed search+audit checkpoints are verified and
reused — **only scoring reruns**, under 128 GiB requested / 110 GiB sampled-RSS
stop on octopus10; deviation documented in `timing-10000-execution-v3/RESOURCE_DEVIATION.md`.
Measured-basis estimates (10k): search **164 h/assembly** (was 178.9 h
extrapolated), cohort search ≈ **76,000 CPU-h**; audit 11.6 h/assembly ≈ 5,400
CPU-h cohort; scoring time unmeasured (memory-bound). Storage scenario
(19.4 TB compressed) roughly confirmed at 5,758 rows/guide.

**GRCh38 chain failed twice more, both diagnosed and fixed:** 2874811 died at
45:28 because the plain-FASTA decompression called `cp` with stdin/stdout
(conversion itself completed; partial gz preserved as evidence); 2875375
died in 7 s on `FileExistsError` — the preserved-evidence output dir blocked
`mkdir(exist_ok=False)`. The converter now takes `--out-dir` (fresh v2 dir;
failed v1 dir preserved), writes a `fai_sha256` provenance record, and
produces byte-deterministic gzip. 35 unit tests + end-to-end mini-run pass.
Resubmitted chain: **2875898 convert → 2875899 build → 2875900 recheck →
2875901 promote+rebind** (dead 2875376–78 cancelled; STOPs preserved as
`STOP-2874811-evidence`/`STOP-2875375-evidence`).

**SBDS: 335/466 COMPLETE** (search array 2875379 all 48/48 exit 0, ~31 min
per assembly): 64,118,953 audited hits, all CFD-defined, 14,398,176
CB-scored, 49,720,777 explicit N/A. Remaining 131 = 128 INDEX_EIO (wave-5
2875393→94→95 pending on cluster occupancy — all octopus nodes currently
fully allocated by other users; this is normal, do not cancel) + GRCh38
(fix running) + HG06807_mat/pat (never-built; wave-5 includes them; source
FASTAs verified intact). Wave-4's 48 replacements promoted and searched;
**independent validation 2875903 queued** (same pinned comparator as
wave-2/wave-3 PASS gates).

Next: monitor 2875902/2875898–901/2875393–95/2875903; wave-6 for the final
~80 EIO; GRCh38 search after rebind; measured 10k report once scoring
completes; SBDS deliverables at 466/466.


## Session 2026-11-16 (afternoon): SBDS 287/466, waves 4 done, wave 5 queued

**GRCh38 convert 2874811 failure root-caused and fixed:** the conversion
itself completed (`GRCh38_acgtn.fa.gz`, 836 MB, written fully) but my
decompression step called `cp` with stdin/stdout redirection, which `cp`
rejects — the job died at 45:28 with a 0-byte `.fa` and no
`transformation.json`. Corrected to pure-Python streaming
(`shutil.copyfileobj`) and made the derived gzip byte-deterministic
(`gzip.GzipFile(filename='', mtime=0)` — the filename field and mtime are
the two nondeterminism sources; **5 conversion unit tests + end-to-end
mini-verification with the pinned samtools now pass**, including the
coordinate-preserving `.fai` equality check). Failed-run outputs preserved
(`*.partial-2874811`), STOP latch preserved (`STOP-2874811-evidence`), dead
chain 2874812–14 cancelled, and the reviewer's vacuous `.fai` check in
promote-rebind-run.sh replaced with a digest comparison against the new
conversion-written `fai_sha256` record. Corrected chain:
**2875375 convert → 2875376 build → 2875377 recheck → 2875378 promote+rebind**.

**Duplicate-submission incident (resolved):** piping the driver through
`tail -1` hid that the first invocation had already submitted **2875379**;
a second invocation submitted **2875384** on the same 48 assemblies. The
duplicate was cancelled within seconds (tasks ≤6 s old; restartable worker
cleans partial attempts); 2875379 retained. Lesson recorded: never pipe
`run_sbds_incremental.py` output through filters, and check `squeue` before
re-invoking.

**Overnight results (all exit 0):** SBDS search wave 2874795 **48/48**;
wave-4 chain 2874802/03/04 **48/48 each**; wave-3 independent validation
2874815 **48/48 PASS** (`wave3-independent-validation-v1/`, pinned
comparator). **Ledger v9: 336 READY_BOUND / 128 INDEX_EIO / 2 MISSING;
SBDS 287/466 COMPLETE** — 54,931,603 hits, all CFD-defined, 12,336,740
CB-scored, 42,594,863 explicit N/A.

Active: **2875379** search wave (48 wave-4 assemblies, octopus08 %4,
GRCh38 still excluded); **2875393 → 2875394 → 2875395** guarded wave-5
chain (48 disjoint, **includes both never-built HG06807 sets** whose
source FASTAs verified intact — 900 MB each, gzip-CRC clean);
**2875375–78** GRCh38 fix chain; **2874508** Brunello 10k at ~20.6 h,
~20.5 GiB peak tree RSS, still in buffered single-invocation search
(expected ~23.1 h traversal + audit + annotation).

Next: wave-6 for the final ~80 EIO + GRCh38 search once rebound; further
search waves as promotions complete; Brunello measured report on 2874508
completion; SBDS final deliverables at 466/466.


## Session 2026-11-16 (early): SBDS 239, wave-4 queued, GRCh38 root-caused

**Correction to the first GRCh38 chain (2874807–810):** the hand-written
slurm wrappers expanded `$D` to an empty string at generation time; job 2874807
failed in 1 s on `sha256sum` and the dependent jobs could never run. The
independent reviewer found the same blocker. All four wrappers were
regenerated with runtime expansion, and the reviewer's additional gaps were
fixed: superseded-binding target preflight (`test ! -e` prevents mv nesting),
an ERR trap restoring the binding rename if fresh binding fails, and a derived
`.fai` digest check before binding. Dead chain cancelled; corrected chain
resubmitted: **2874811 convert → 2874812 build → 2874813 cross-node recheck →
2874814 promote+fresh-bind** (frozen `00_manifest/grch38-fix-tooling-v1/`).
Also submitted **2874815**: wave-3 independent reference-to-indexed-text
validation for all 48 wave-3 replacements (`wave3-independent-validation-v1/`,
same pinned comparator as wave-2's 48/48 PASS), closing the audit gap the
SBDS scout flagged (wave-3 had only build verify + cross-node recheck +
promotion evidence until now).

Overnight the full wave-3 recovery chain completed (48/48 build, 48/48
cross-node recheck, 48/48 promote/bind; arrays 2874447/48/49 all exit 0) and
the wave-2-ready SBDS search wave (2874424) completed **48 of 49** tasks.
**Ledger v8: 288 READY_BOUND / 176 INDEX_EIO / 2 INDEX_FAILED_MISSING_INDEX;
SBDS 239/466 COMPLETE**, 45,739,759 audited hits, all CFD-defined,
10,273,454 CRISPR-Bulge-scored, 35,466,305 explicit N/A.

**GRCh38 root cause found (the one 2874424 failure):** its PanSN analysis
reference contains **94 IUPAC ambiguity positions** (Y/R/W/M/K/S/B) and the
native PAM search accepts only A/C/G/T/N, so every GRCh38 search fails fast
with `reference must be nonempty A/C/G/T/N` (attempts 2873071/2873638/2874424,
all preserved). CHM13 confirmed COMPLETE — lowercase soft-masking is accepted.
Fix implemented: audited IUPAC→N conversion (`grch38_ambiguity_fix.py`, full
per-position record, coordinate-preserving .fai check; 3 unit tests — one
caught a real write bug), `rebuild_sbds.py` gained a `--fasta-gz` override
recorded as `source_override` in completion manifests. Guarded chain submitted:
**2874807 convert → 2874808 build → 2874809 cross-node recheck → 2874810
promote+fresh-bind** (old binding renamed `GRCh38-ambiguity-superseded-v1`, never
deleted; promote uses the original_existed=False path — no /export original is
touched). Frozen tooling: `00_manifest/grch38-fix-tooling-v1/`.

Continuations launched this session: **2874795** SBDS search wave for the 48
newly promoted wave-3 assemblies (octopus08, %4; GRCh38 excluded until fixed
via `exclude-grch38-pendingfix.txt` — remove that exclusion after rebind);
**2874802 → 2874803 → 2874804** guarded recovery wave-4 for 48 disjoint
assemblies of the remaining 176 EIO (frozen `recovery-wave4-tooling-v1`, afterok
chain, %2 per stage). `prepare_recovery_wave.py` now accepts both wave-2 JSON
and wave-3 sha256sum tooling pin formats (tested).

**Brunello 10k benchmark (2874508) RUNNING HEALTHY at ~7 h**: steady ~16.5 GiB
sampled tree RSS (48 GiB stop), 5 processes, mid-search (single-invocation PAF
writes at traversal end; 1k-derived expectation ~23.1 h traversal + ~1.5 h
audit + ~11 h annotation ≈ 35.6 h). No guard trips, no stderr output.
Three read-only independent reviews in flight (GRCh38 fix soundness, SBDS/wave
progress audit, full-scale Brunello estimate framework).

Next: finalize the 10k measurement report when 2874508 completes (runtime,
guides, search/audit/score walls, hits, CB availability, peak RSS, compressed
sizes, audit PASS) and derive full-library estimates; after wave-4 completes,
launch wave-5 (128 EIO remain) and further search waves; rerun GRCh38 search
after rebind; SBDS final deliverables once 466/466.

## Session 2026-11-15 (later): wave-2 validated 48/48, SBDS at 199, Brunello 10k launched

Milestones since commit `7b93c4f`:

- **Wave-2 independent validation COMPLETE: 48/48 PASS, no STOP**
  (`00_manifest/wave2-independent-validation-v1/`, pinned wave-1 comparator).
  The `afterok` gate cleared and **wave-3 rebuilds (2874447) are running** on
  octopus10 — first builds already `STAGED_WRITTEN_AWAITING_CROSS_NODE_RECHECK`;
  recheck (2874448) and promote/bind (2874449) remain dependency-queued. All
  guarded; wave STOP absent.
- **SBDS search 2874424: 8 tasks COMPLETED** (0–7, ~30 min each), tasks 8–11
  running, rest pending at `%4`. **Ledger: 199 COMPLETE / 267 remaining**;
  38,099,500 hits, all CFD-defined, 8,556,941 CB-scored, 29,542,559 explicit
  N/A. No new FAILED markers.
- **Brunello bounded 10,000-guide timing benchmark LAUNCHED: job 2874508** on
  octopus06 (disjoint from SBDS octopus08 and recovery octopus10/07), one CPU,
  56 GiB, 72 h, 48 GiB sampled tree-RSS stop, guarded by the reviewed
  wave_guard → frozen supervise.py chain. Frozen tooling + checksummed
  execution plan: `runs/brunello_pam_aware_v1/00_manifest/timing-10000-execution-v2/`
  (plan SHA256 `d5226c6c8002646581b135ebc4c1d13d9f5b9567dd26df6b5b5c4a57ca23fd9b`).
  Output: `runs/brunello_pam_aware_v1/01_pilot/timing-10000-v1/HG00735_pat/`.
  Estimated ~35.6 h (traversal ~23.1 h + audit ~1.5 h + annotation ~11 h).

**Brunello wrapper** (`research/hprcv2/brunello_expanded_benchmark/timing_worker.py`):
bounded to the frozen 10k plan (rejects cohort scope/config drift); immutable
attempt dirs; checksummed search/audit/score checkpoints with atomic
no-clobber publication (crash before publication leaves no false checkpoint);
scratch reference alias with recomputed and pin-verified `.fai`; explicit
CB availability accounting (`scored` vs `unsupported_*`, never zero-filled);
scoring stdout/bundle byte-equality and `metadata.json` input-identity checks;
final `SUCCESS.json` only after the external supervisor reports success.
**10 unit tests pass**, including crash/tamper/pin-drift/lock/supervisor-gate
failures and restart reuse without duplicated work. Real on-node preflight
**2874502 PASSED** (all pins, FAI provenance, scratch 145 GB, output 34.8 TB
free). A first preflight (2874491) correctly failed closed on a node-dependent
system-interpreter pin; system python is now recorded, not pinned (v1 frozen
dir retained as evidence). Independent reviewer attempts were blocked by
provider usage limits and a scheduler admission timeout — three attempts,
all infrastructure failures — so launch proceeded on the user's stated gate
(wrapper tests + resources), with this handoff recording the deviation.

Next: monitor 2874508 (first full-scale m4 linearity evidence); continue SBDS
search wave to 49 completions; wave-3 recheck/promote then a new guarded
recovery wave for the ~178 still-EIO assemblies (exclude frozen wave-3 list);
schedule newly bound searches with active-list exclusions. Full Brunello
cohort remains gated on capacity review and the timing evidence.

## Astra resume 2026-11-15: SBDS restarted, guarded wave 3 queued

Verified main model **openai-codex/gpt-6-astra**, branch
`feature/pam-anchored-search-poc`, starting HEAD
`e2687245e972df18441aedef00f15e04aa6a487c`. Three independent **read-only**
`scout`/`gpt-6-luna` reviews completed successfully (SBDS, performance,
Brunello). The earlier `researcher` tool-allowlist failure was avoided by
using the repository scout profile; no main-model configuration was changed.

Fresh evidence-derived **manifest v7: 240 READY_BOUND / 224 INDEX_EIO /
2 INDEX_FAILED_MISSING_INDEX = 466**. SBDS baseline: **191 COMPLETE /
275 remaining**, 36,568,238 hits, all CFD-defined, 8,213,194 CB-scored and
28,355,044 explicit N/A. This is summary/completion evidence, not a fresh
full-output hash scan. All 48 tasks in each wave-2 build/recheck/promote array
2873752/53/54 completed exit 0.

Main-agent approved jobs (no completed searches included):

- **2874424**: 49 pending READY_BOUND assemblies (48 wave-2 + GRCh38),
  SBDS four-PAM m5/b1/z2/f0.75, `%4`, octopus08, 1 CPU / 56 GiB / 12 h.
  Four tasks running at checkpoint, with no new FAILED markers observed.
  Frozen plan/list/receipt: `00_manifest/sbds-{plan,assemblies}-wave2-ready-20261115.*`.
- **2874431**: wave-2 independent full-component/source/binding/indexed-text
  validation, `%2`, octopus07, 1 CPU / 8 GiB / 2 h. **13/48 PASS at the saved
  checkpoint**; still running. Uses preserved wave-1 comparator SHA256
  `2306463cb7a94f1da01842b554be4c1e5da83542a3d6905b6b698239ee54f89a`.
  Evidence: `00_manifest/wave2-independent-validation-v1/`.
- **2874447 → 2874448 → 2874449**: 48 disjoint wave-3 assemblies,
  build → cross-node recheck → rollback-safe promote/bind. **Queued, not yet
  running**: whole-array `afterok:2874431` before build, then whole-array
  `afterok` between stages. Build octopus10 `%2`, 1 CPU / 64 GiB / 4 h;
  recheck/promote octopus07 `%2`, 1 CPU / 8 GiB / 2 h. Frozen guarded tooling,
  assembly list, hashes and submission receipt:
  `00_manifest/recovery-wave3-tooling-v2/`. v1 was prepared but never submitted;
  retain it as superseded review evidence, do not launch it.

New `prepare_recovery_wave.py` only freezes a bounded disjoint wave; it never
submits. New `wave_guard.py` checks shared STOP, scratch/output free space and
available memory; any guarded failure blocks later tasks and terminates current
wave descendants across process sessions, even if STOP persistence fails.
Frozen tooling checks run inside the guard. Independent `reviewer`/`gpt-6-sol`
review found and resolved five fail-closed issues; final review accepted the
code with launch gates. **30 SBDS unit tests pass**, including 10 guard tests;
Slurm guard+supervisor smoke job **2874444** passed. Shell syntax and frozen
checksums verified, and `scontrol` confirms the exact whole-array dependencies.

Resource preflights 2874423/2874429/2874430: octopus08/10/07 each ~400 GB
available memory; scratch free ~411/623/341 GB respectively. MooseFS free
~34.7 TB. Early SBDS sampled tree peaks ~13.2–13.7 GiB. Slurm CPU/RSS
accounting is missing/zero for historical matched benchmark jobs; do not
interpret that as zero usage. Node SSH failed host-key validation; no SSH
security settings were weakened—resource checks used bounded Slurm tasks.
No running job was canceled; unrelated job 2873575 remains untouched.

Performance report corrected: 1k new wall 489.76 s versus old 383.25 s is a
1.28× wall ratio with unlike stage boundaries, **not a CPU-efficiency ratio**.
Old process CPU 305.95 s; new comparable CPU unavailable. Historical full-library
5,775 s measured versus new 14,188 s extrapolated is not a measured comparison.
All 28,750 eligible old witnesses contained, but 13,184 CFD comparisons disagree.
No formulas/algorithms/binaries changed. Brunello review and concrete wrapper
failure-test requirements are in `research/hprcv2/sbds_nrr/BRUNELLO_WRAPPER_REVIEW.md`;
**the 10k wrapper remains unfinished**. No performance/Brunello jobs submitted.

Next: monitor wave-2 validation and the dependency-gated wave-3 chain; halt
escalation on any new integrity error; refresh completion ledger and schedule
newly bound searches with active-list exclusions. Do not rerun matching completed
assemblies. Additional recovery waves must exclude the frozen wave-3 list and
any staged/bound sets. Rebuilt outputs remain under MooseFS `indices/` after
promotion. **No quarantine deletion this session**; all prior independent
validation, full-search-checksum, supervisor, and reviewed-deletion-manifest
gates remain mandatory. Full SBDS figures await all 466 completions. Session
receipt: `research/hprcv2/sbds_nrr/results/resume-20261115.json`.

## Resume checkpoint 2026-11-15: wave-2 finished, campaigns idle

Read-only recovery of live campaign state at 17:26 UTC:
`campaign_snapshot.py` reports **240 READY_BOUND**, 224 INDEX_EIO, and
2 INDEX_FAILED_MISSING_INDEX (466 total); **191 SBDS COMPLETE / 275 remaining**,
36,568,238 hits, all CFD-defined, 8,213,194 CB-scored and 28,355,044 explicit
CB N/A. These are completion/summary evidence counts, **not a fresh full-output
hash verification**. Snapshot: `/tmp/crispr-campaign-resume-20261115T172611.json`
(local temporary report, not authoritative persistent evidence).

All 48 tasks in each wave-2 array 2873752 (build), 2873753 (recheck), and
2873754 (promote/bind) are COMPLETED with exit 0:0 according to `sacct`.
No CRISPR campaign jobs remain queued/running; the user's unrelated
`bn_verkko_gt100kb` job is still running. Independent reference comparison of
wave-2 replacements has **not** been established by this resume check.
Only the two previously documented quarantine deletions are recorded;
MooseFS free space is approximately 34.7 TB decimal. No new jobs or deletions
were launched. Next: independent wave-2 validation, reconcile the pending
search cohort against a fresh authoritative manifest, then continue the
controlled rebuild/search campaign. Brunello remains launch-gated.

## Session 2026-11-14: wave-1 validated, wave-2 running, deletions gated

**Wave 1 complete and independently validated: 48/48 rebuilt indexes** passed
cross-worker recheck, promotion, binding, AND an independent reference-to-
indexed-text comparison (preserved comparator, pinned SHA256) — evidence in
`00_manifest/wave1-independent-validation-v1/`. Ledger v6: **192 READY_BOUND /
274 remaining** (44 wave-2 builds staged pending recheck, 2 still EIO-pending).
**SBDS: 89 assemblies COMPLETE** (17.03 M hits, 100 % CFD, 3.82 M
CRISPR-Bulge defined, rest explicit N/A); active: wave1c (121 tasks, %4,
octopus08) and wave1-new48 (48 tasks, %2, octopus09) — the driver now requires
explicit `--exclude-assemblies-file` lists whenever cohort jobs are queued
(no duplicate submissions) and writes checksummed submission records.

**First two quarantine deletions executed under full gates:** HG005_pat and
HG00099_hap2 corrupted originals (22.59 GB logical) removed only after
independent replacement validation + COMPLETE checksummed four-PAM search +
supervisor success, with pinned review manifests and per-file deletion logs
(`review_quarantine.py`; regular-file unlink only, inode/mtime re-checked,
symlink-swap detection, PARTIAL on failure; 5 unit tests). The old
`delete_original.py` is REFUSED/withdrawn (it "freed" storage by renaming and
hashed unreadable files). Wave 2 (48 disjoint assemblies) submitted with
monitored build→cross-check→promote chains (2873752/53/54, %2).

**Brunello (NOT launched):** timing claims reconciled — 178 h is
traversal-only; full audited both-score workflow ≈ **276 h (11.5 days)** per
assembly at m4; storage corrected: authoritative PAF+annotations ≈ **13 TB**,
audit records ≈ 6.4 TB, total ≈ 19.4 TB (previous 5.4/14 TB split was wrong —
see `TIMING_RECONCILIATION.md`). 10,000-guide timing inputs/commands frozen
(`brunello_pam_aware_v1/00_manifest/timing-10000-v1/`) but its
restart/monitor wrapper is NOT yet launch-ready; submission remains gated on
resource review. `campaign_snapshot.py` is the live evidence-derived ledger.

## Rebuild+SBDS campaigns RUNNING (2026-11-13, authorized)

User authorized proceeding without waiting for the admin: controlled rebuilding
AND concurrent SBDS searches, with ONE authoritative index directory
(`indices/` — no new /export datasets). Promotion pipeline implemented and
proven (`promote_index.py` quarantine+rollback renames, `repoint_binding.py`
keeps binding sidecars byte-identical so completed searches stay valid,
`delete_original.py` deletes quarantined originals only after promoted
replacement passes cross-node verification + binding + COMPLETE checksummed
SBDS search, preserving error reports and a deletion log). **All 11 previously
rebuilt sets promoted into `indices/`, originals quarantined; bindings
re-pointed with unchanged sidecars; HG03784_hap1's completed search remained
valid (no rerun).**

Active: SBDS four-PAM search wave (array 2873335, 141→129 remaining after 15
COMPLETE with 2.87 M audited hits, 640,646 CB-scored, rest explicit N/A;
first submission was correctly stopped 129× by the octopus05 scratch-full
guard and resubmitted on scratch-healthy nodes); rebuild wave 1 of the
remaining 322 sets (48 assemblies, chained arrays 2873212 rebuild → 2873213
reverify → 2873214 promote+bind, all restartable). Cluster is currently
heavily allocated by other users, so arrays advance as resources free —
all tooling is restart-safe; re-run `run_sbds_incremental.py --submit` and
`rebuild_wave` lists to continue. Deletion of quarantined originals runs via
`delete_original.py` as assemblies' searches complete. Brunello: full
readiness review prepared ([BRUNELLO_READINESS.md], includes storage math
(19.4 TB vs 41.7 TB free), the 128k CPU-h model, safe throughput-only
reductions, and the proposed 10,000-guide linearity timing experiment —
launch gated on final capacity review). `status_report.py` is the live
ledger. Manifest v5: 144 READY_BOUND / 318 INDEX_EIO / 4 MISSING_INDEX.

## Storage policy correction implemented (2026-11-13)

Per the user's correction, ALL new genomic data now lives under
`/moosefs/salehi/HPRCv2_CRISPRapido/` — rebuilt indexes in `indices_rebuilt_v1/`,
results in `runs/sbds_ngg_nrr_v1/`; no new persistent datasets under /export.
See [STORAGE_POLICY.md](../research/hprcv2/sbds_nrr/STORAGE_POLICY.md).
Rebuild pipeline rewritten: build on node scratch → full component verify →
MooseFS placement with fsync + write-back hash verification → **cross-worker
re-read gate** before a set is trusted; checksummed completion manifests;
originals never modified or deleted (deletion needs a reviewed manifest, not
implemented). Reliability evidence: 2 GB synthetic test + 3 migrated pilot
sets + 8-set rebuild batch — **~137 GB of new MooseFS writes verified from
independent workers, zero mismatches**; HG03784_hap1 search rerun after
index relocation gave **identical results** (192,353 hits / 42,415 CB /
149,938 N/A; old SUCCESS preserved as SUPERSEDED). All 12 trusted-recovery
assemblies bound. **Manifest v4: 144 READY_BOUND / 318 INDEX_EIO / 4
MISSING_INDEX = 466.** 124-assembly incremental search plan prepared but
NOT launched; mass rebuilding of the remaining ~322 sets is deliberately held
pending the administrators' answer on MooseFS health (updated
[ADMIN_MESSAGE.md](../research/hprcv2/sbds_nrr/ADMIN_MESSAGE.md)). /export
pilot copies preserved pending deletion review.

## Stage-1 verification campaign COMPLETE (2026-11-13, authorized)

[VERIFICATION_CAMPAIGN.md](../research/hprcv2/sbds_nrr/VERIFICATION_CAMPAIGN.md).
Authorized 438-assembly campaign (array 2871846, %12, 1 h 22 m wall, within
the 5.3 TB budget): **115 newly verified+bound → 126 READY_BOUND total;
326 INDEX_EIO; 4 MISSING_INDEX; 10 TRUSTED_EVIDENCE** (authoritative
manifest-v3, exactly 466, mutually exclusive). Sequential %1 rescan of a
10-assembly EIO sample reproduced identical persistent failures — real
chunk-level damage, not load artifacts. Actual usage: 1.41 TB measured
success reads, 360 GB committed bindings, 931 GB partial stages preserved
(archived, not deleted). Recovery plan: bind the 10 trusted (~26 min),
rebuild/build 330 sets onto /export (~116 one-CPU h, ~4.9 TB, RAM-bounded
~41 GB each) pending admin replica-recovery outcome (admin message ready,
updated scope ~333 damaged sets), then incremental SBDS search via the
driver (136 assemblies ≈ 68 CPU-h at first). No rebuilds or searches launched.

## Corrected estimates, admin message, campaign preparation (2026-11-13)

[CORRECTED_ESTIMATES.md](../research/hprcv2/sbds_nrr/CORRECTED_ESTIMATES.md):
the earlier "55–60 rebuild-days" figure was serial one-CPU HOURS mislabeled
as days. Corrected: rebuilds cost **0.34 one-CPU h per assembly end to end**;
~190 rebuilds ≈ **65 one-CPU hours (2.7 core-days)**, wall ≈17 h at %4, ≈8 h
at %8 (real limiter is ~41 GB RAM per build), ~2.8 TB on /export.
[ADMIN_MESSAGE.md](../research/hprcv2/sbds_nrr/ADMIN_MESSAGE.md) is ready to
send: affected store, reproducible errno-5 evidence, healthy source FASTAs,
request for chunk/replica recovery. Verification campaign for the 438
UNATTEMPTED is fully prepared (`verify-campaign-unattempted-438.txt`, sbatch
validated with --test-only; ≤5.3 TB reads, ~19 CPU-h, %12) — **awaiting
explicit I/O authorization**. SBDS incremental driver prepared and
dry-run (11 READY_BOUND, 8 pending, plan sha
`966f4cb9291646e74810ed396a36108b67330f62a5e251922d0128d6202cb5d9`): worker
now records per-PAM scientific search signatures so completed assemblies
are never rerun across plan expansions; 10 regression tests pass. Brunello
untouched.

## MooseFS incident, reconciliation and recovery pilot (2026-11-13)

[RECOVERY_PLAN.md](../research/hprcv2/sbds_nrr/RECOVERY_PLAN.md) and
[INCIDENT_REPORT.md](../research/hprcv2/sbds_nrr/INCIDENT_REPORT.md).
Ledger corrected: the previous "20+10+438=468" double-counted HG02055_mat
and NA18906_mat. Authoritative manifest v2 has **exactly 466 unique,
mutually exclusive** assemblies: 11 READY_BOUND + 10 TRUSTED_EVIDENCE +
7 INDEX_EIO + 438 UNATTEMPTED. Incident scan of all ten failed sets:
**persistent, reproducible errno-5 unreadable regions** (not transient, not
missing files, **not corrupted assemblies** — every source FASTA.gz verifies
with valid CRC/SHA256). Damage is distributed across the MooseFS index
store and grew since 2026-10-28; admin-facing report with paths/offsets/
reproducer is ready. Controlled rebuild pilot (3 assemblies, onto healthy
/export, jobs 2871822): **~16–17 min build, ~40–42 GB peak RSS, 11.1–11.55 GB
index each**, all components hashed; all three bound (2871828–30) and
**HG03784_hap1 passed a complete four-PAM verification search** (2871832,
30 min, 192,353 audited hits, full scoring). Campaign plans prepared but NOT
launched: verification ≤5.3 TB reads / ~19 CPU-h; rebuild scenario ~190 sets
≈ 55–60 rebuild-days at %4, ~2.9 TB on /export; do NOT rebuild onto MooseFS
before cluster-team diagnosis (replica recovery may avoid most rebuilds).

## SBDS deep-verification batch 1 (2026-11-13; campaign not launched)

[VERIFICATION_BATCH1.md](../research/hprcv2/sbds_nrr/VERIFICATION_BATCH1.md):
16-assembly controlled Stage-1 batch (array 2871790, %4 concurrency, one
CPU, octopus06/07/08; octopus05 excluded for 16 G scratch). **6 deep-verified
bindings committed; 10 EIO failures, originals untouched.** Measured per
successful assembly: 12.16 GB read, 147–172 s; committed bindings ~3.06 GB
each; reusability recheck passes. EIO prevalence: 6/7 historically notable
assemblies (all four historical no-perfect-match cases, HG00321_hap2,
HG03784_hap1) plus 4/9 of an ordinary spread ≈ **44% of ordinary assemblies
may need index rebuilds**. Ledger v2: **20 of 466 currently search-ready**
(12 trusted + 8 committed bindings), 10 EIO recorded, 438 unattempted.
Remaining-campaign scenario: ≤5.3 TB reads, ~19 CPU-hours, ~1.4 TB bindings —
**safe to proceed with verification in batches**; a complete search also
needs an authorized, separately piloted **index-rebuild campaign** and a
cluster-team escalation of the MooseFS degradation. Rebuild cost unmeasured.

## SBDS 466-genome readiness (2026-11-13; cohort NOT launched)

Infrastructure blockers resolved and a controlled multi-assembly test passed:
[COHORT_READINESS.md](../research/hprcv2/sbds_nrr/COHORT_READINESS.md).
Ledger audit: **12 assemblies** have trusted full-content index evidence
(incl. GRCh38/CHM13); **454 require deep verification** (~5.6 TB component
reads, ~20-35 CPU-hours at measured rates; EIO assemblies will stop with
explicit failures, never zero-filled). Expansion tooling rewritten and
unit-tested (restart-safe checkpoints, tamper detection, process-group kills,
per-assembly pins; 6 regression tests). Two new bindings (HG02055_mat,
NA18906_mat) committed via the repaired trusted binding helper. Cohort test
array 2871780 (one CPU, 56G): 192,378/192,148 total hits per assembly, all
CFD-defined, 42,465/42,389 CRISPR-Bulge-defined (NRR scorer byte-identical
to the officially validated bundle), 149,913/149,759 explicit N/A, all eight
(assembly,PAM) audits PASS. Searches ~28-30 min/assembly on one CPU; full
466 scenario ≈233 one-CPU job-hours, ~4.75 GB compressed outputs, bindings
~1.45 TB on /export. Remaining blockers: the 454-assembly deep-verification
campaign (needs I/O budget + EIO-rebuild authorization), noncanonical
biological-validity caveat, historical EIO/no-perfect input investigations.
Recommended submission settings are in the readiness report; nothing submitted
for the cohort.

## SBDS resumed: NGG/NRR pilot completed; NRR scoring validated; cohort still gated

Latest SBDS request supersedes the PAM pause below: standard SpCas9 NGG is
confirmed; NRR coverage is requested for alternative PAMs. See
[NGG/NRR report](../research/hprcv2/sbds_nrr/REPORT.md) and the
[saved-pilot NRR reannotation](../research/hprcv2/sbds_nrr/REANNOTATION.md).
Four unchanged literal-PAM searches on verified HG00735 completed under
one-CPU Slurm array 2871759: NGG 38,060; NAG 46,358; NGA 42,918; NAA 64,712.
All 192,048 hits independently audited; all CFD defined. NGG both-score output
is byte-identical to prior primary.

**Saved-pilot reannotation (searches NOT repeated):** the Python scorer now
supports NRR target suffixes within the unchanged mismatch/gap domain; guide
NGG / actual target PAM. Neural defined **42,425** total: NGG 8,104, NAG
10,452, NGA 9,421, NAA 14,448 (**+34,321 noncanonical**). Remaining 149,623
retain N/A with reasons (>24 columns 70,474; gapped mismatch domain 69,176;
multibase RNA 9,973). All **42,425** supported predictions independently match
official inference (max diff 2.2e-9), plus 8,272 synthetic domain fixtures
(max diff 6.1e-13); categorical features hash-equal. Default CFD pass-through
byte-identical on all four inputs. Reannotation job 2871764 (103.7 s, peak
0.944 GiB sampled); domain job 2871765. **Noncanonical biological validity is
NOT established**; scores are experimental. Full 466-genome rerun remains
BLOCKED and explicitly stopped pending further authorization.

Root: `/moosefs/salehi/HPRCv2_CRISPRapido/runs/sbds_ngg_nrr_v1/`.
Serial workflow 28.7 minutes, maximum sampled tree RSS 13.64 GiB; naive 466-fold
scenario 222.9 one-CPU job-hours and 4.75 GB compressed PAF+annotations, excluding
binding preparation and diagnostic copies. No 466-genome rerun or population
figures yet: complete trusted binding ledger and safe tested expansion tooling
are missing. Earlier claims below that draft full-run helpers are ready are
**withdrawn**: binding helper contains closed-file fsync/string-path/digest
errors; full-run worker has chunk stream/status/checkpoint defects. Do not use
them for this expansion until repaired and failure-tested. Pilot used only the
already verified HG00735 binding. Brunello remains separate and unlaunched.

## Current development state

- Repository: https://github.com/FarnazSalehi94/crisprapido
- Checkout: `/export/local/home/salehi/crisprapido-pam-validation`
- Development branch: `feature/pam-anchored-search-poc`.
- Whole-assembly experiment milestone: **`dfd1636`** (three-guide batch completed); handoff documentation baseline: **`7c4e77f`**; CFD feasibility milestone: **`252a6bb`**. The subsequent downstream filtering/isolated native pruning implementation is linked below; pruning is **not enabled in the integrated backend**. The earlier handoff ending at `bfbb8ac` is superseded; **do not submit the three-guide batch as an unperformed experiment**.
- Experimental backend is integrated but is neither the default nor a production release. The original **Columba candidate search → WFA2 verification → PAM filtering → CFD → PAF** backend remains the default and is preserved. Earlier original-versus-Columba/WFA2 benchmarking is separate from the experimental backend results below; do not repeat it without a scientific reason.
- `--search-backend pam-anchored` uses Columba FM-index primitives, seeds the PAM constraint, tracks mismatches/bulges, handles both strands and reference Ns, and reports one deterministic representative per PAM anchor without WFA2. CFD scoring remains downstream.
- Keep CLI-configurable settings **`-m 5 -b 1 -z 2 -f 0.75 -p GG`**; `-p GG` means strict NGG. Do not change defaults or representative selection without discussion. Representatives minimize `(mismatches + inserted bases + deleted bases, mismatches, compressed CIGAR, start)`, not CFD.

## Active HPRCv2 SBDS / Brunello request (2026-11-11)

See [independent project audit and guarded pilot plan](../research/hprcv2/README.md)
and [downstream allele/IMPG requirements](../research/hprcv2/ALLELE_PLAN.md).
Separate run trees now exist under `/moosefs/salehi/HPRCv2_CRISPRapido/runs/`.
Initial audit Slurm 2871688: frozen 464+2 cohort and 77,441-guide library verified;
two historical SBDS PAFs and ten Brunello small-artifact checks currently return
EIO. Earlier successful historical QC must not be described as current full
readability. Prior resolved index evidence covers 12 sets, not all 466.

The requested SBDS primary search is m5/b1/z2/f0.75/GG with both scores and no
CFD filter. Historical runs used m2/b1/z1, so a separate matched-domain control
is needed. Safety pilots on the already bound HG00735_pat assembly are now
complete: SBDS primary/control Slurm **2871689/2871691**, **38,060/61** eligible
hits, **8,104/61** neural defined. Brunello A=m4/b1/z1 and B=m5/b1/z2 share three
exact frozen-library GC-stratified guides: jobs **2871690/2871692**, **13,474/70,364**
hits, **13,474/14,924** neural defined. All **121,959** rows independently audited;
**95** sampled predictions match official inference within 1e-6 (jobs 2871693/4).
See [SBDS pilot](../research/hprcv2/SBDS_PILOT.md) and
[Brunello review report](../research/hprcv2/BRUNELLO_PILOT.md), including actual
resources/compression and explicitly unvalidated 76k/874k native-core-hour linear
full-library scenarios. Forty-eight Python tests passed; core Rust unchanged.

On HG00735, all **57/57** strict-eligible historical CIGAR witnesses recovered;
old raw had **756 non-NGG/813**. Prior chr12/17 patterns have **GGC/AGT** PAMs
and are correctly excluded. Historical gapped CFD discrepancies also remain:
25 original eligible CIGAR scores differ from independent current reconstruction;
24 shared unchanged representatives have changed CFD. Do not attribute every
score change to representatives or claim new biological calibration. Brunello B
contains every A anchor but changes 303 representatives/neural availabilities;
thresholded A/B views are not nested despite nested eligibility.

**No full Brunello launch is authorized before user review.** Full SBDS expansion
awaits remaining integrity/binding gates; no 466-task search array has been
submitted. Two historical SBDS PAFs and ten Brunello small-artifact checks return
EIO; 454 index sets still lack full-content evidence here. No completed
population rerun/updated six SBDS figures, IMPG mapping or allele-frequency
results are claimed yet. External figures are clearly labelled n=1/three-guide
pilot diagnostics. Completed prior HG00735 experiments below are unchanged and
must not be repeated as if new.

## Brunello expanded benchmark and full-run preparation (2026-11-12; SBDS paused)

**SBDS is PAUSED pending confirmation of the prime editor's PAM; no SBDS
search/parameter work is in progress.** The completed bounded benchmark is in
[brunello_expanded_benchmark/REPORT.md](../research/hprcv2/brunello_expanded_benchmark/REPORT.md):
16 timing-blind stratum-selected guides x 2 verified assemblies (HG01496_pat
safely bound this session, 371 s trusted import) x 4 configurations, plus a
1,000-guide library-mix scaling benchmark at m3 and m4; **11.7 M rows audited
with zero errors; nothing full-scale launched**. Key measured results: m3/b1/z1
≈ 1.45 s/guide, 523 rows/guide, 100 % both scores; m4/b1/z1 ≈ 8.3 s/guide,
5,855 rows/guide, 100 % both scores, and it adds ~100 CFD>=0.20 sites per guide
that m3 omits with *perfect nesting* (zero representative/CFD changes at shared
anchors); m5/b1/z2 ≈ 96 s/guide, 83 % neural N/A, 2,667 representative changes
and fewer additional high-CFD sites than m4 — rejected. Full-HPRCv2 campaign
estimates: m4 ≈ 128,000 CPU-hours/~19.4 TB; m3 ≈ 19,000 CPU-hours/~1.7 TB;
m2 rejected (misses high-CFD sites). Recommended scientific config: candidate-A
(m4); pragmatic staged alternative: m3 (forward-compatible because nesting is
perfect). All four configurations are fully compatible with -f 0.75 (the
match-fraction rule is matched/20 and clips nothing at m<=5).

Prepared but NOT launched (Task 6 of the request): prepare_binding.py (in-place
trusted binding for the 12 recovery-verified assemblies; 452 others need the
verification campaign first), make_full_plan.py (requires a recorded user
authorization token), full_run_worker.py (chunked, checkpointed, restart-safe,
1 CPU/assembly, full audits, explicit N/A), submit_full_run.slurm (464-task
array template), verify_full_run.py (ledger; never zero-fills). Storage/filtering
strategy: keep complete eligible both-score PAF + annotations compressed as
authoritative; CFD thresholds are derived views only and can never define
allele-frequency denominators. IMPG/Figure-1C needs remain specified in
research/hprcv2/ALLELE_PLAN.md; nothing mapped or plotted yet.

## Performance investigation (2026-11-11, bounded; no full reruns launched)

See [performance report](../research/hprcv2/performance/REPORT.md). The apparent
7.5-min-for-one-guide discrepancy is ~310 s of fixed per-invocation overhead
(snapshot copy, pre/post full verification, index/reference loads) that batching
eliminates, plus the genuinely expensive requested m5/b1/z2 domain (~101 s/guide
traversal vs 0.196 s at matched m2/b1/z1). Matched-settings benchmark
(jobs 2871709-2871711/2871718-2871720, one verified assembly, identical first
1/100/1000 frozen-library guides, m2/b1/z1/f0.75/GG): new backend marginal
0.179 s/guide wall, ~20,100 guides/h/assembly on one core, all rows
independently audited; old preserved binary 11.1/37.5/383.2 s with 99/85/79%
of one CPU. Old emits ~575 rows/guide with only ~5% strict-eligible; **all
28,750 old strict-eligible witnesses are contained in the new output, zero
missing**. Old-engine full-library wall is measured at 5,775 s (~1.9 cores);
linear extrapolation from N=1000 overestimates it 5.15x, so new-engine
extrapolations (3.94 h/assembly matched) remain unvalidated at scale.
Requested m5/b1/z2 full Brunello would be ~1,900 core-hours per assembly
(~100 core-years over 464 assemblies) — infeasible without authorized
algorithmic/domain change. Batch per assembly, parallelize across assemblies;
full-library/full-cohort launches remain gated on user review.

## Completed validation and whole-assembly experiments

Independent bounded oracle panels recovered **204/204** and **172/172** anchors, with no unexpected or duplicate anchors. Tests include both strands, mismatch/bulge cases, reference Ns, PAMs, alternative alignments, and supported options. Additional bounded large-coordinate probes exercised HG00735_pat positions around `2^31`, contig boundaries and near the assembly end; these are targeted checks, not an exhaustive recall oracle.

All three jobs below used the imported **HG00735_pat** assembly (3,029,527,752 bases, 82 contigs), pilot guide `GCCAGGTACAGGTTCCGTCC`, and the settings above. The batch additionally used exactly the two guides listed in the next table.

| Slurm job | Scope | Reported anchors | Native search including locate/replay | Total runtime |
|---|---|---:|---:|---|
| 2871278 | First successful single-guide pilot | 21,447 | 82.465 s | Supervisor 1,070.97 s |
| 2871283 | Optimized verification/reference handling, same search/scoring policy | 21,447 | 80.5854 s | Supervisor 645.6085467999801 s; Slurm 10 min 47 s |
| 2871284 | Completed three-guide batch, one snapshot/native process/index load | 78,582 | 257.493 s | Supervisor 837.874695419101 s; Slurm 13 min 59 s |

The optimized standalone run reduced observed total runtime by about 40% versus the first successful pilot. This is a historical same-workload observation on a shared system, not a general speedup guarantee. The two standalone PAFs are byte-identical. All 21,447 records passed independent eligibility and CFD audits in each run, with the known exact target recovered.

### Three-guide batch: job 2871284

| Guide ID | Exact guide sequence | Reported anchors / CFD validated | Search states | Native search s |
|---|---|---:|---:|---:|
| pilot | `GCCAGGTACAGGTTCCGTCC` | 21,447 / 21,447 | 856,229,607 | 80.4914 |
| beginning | `CATTATTCTCGGCTGCAAAG` | 28,835 / 28,835 | 926,942,591 | 87.5259 |
| before_2_31 | `TGAACTTCTGATAACTCCCC` | 28,300 / 28,300 | 931,599,680 | 89.4761 |

**All 78,582 reported sites passed independent eligibility and CFD audits. Genome-wide recall remains unproven.** The audit extracted reference targets/PAMs, replayed CIGARs, checked coordinates, strands, mismatch/bulge/match-fraction limits and independently recalculated CFD values (all were supported and matched at output precision). It found zero duplicate anchors and recovered all three independently known exact targets. Counts include those three targets; they are not counts excluding intended on-target sites.

The pilot subset matched job 2871283: 21,447 shared anchors, zero missing/unexpected anchors, zero representative-alignment or CFD differences, and byte-identical PAF records. This demonstrates output stability and within-batch reuse, not an independent genome-wide sensitivity estimate.

- Native index loading: **3.30994 s**, exactly one logged event; snapshot copy: **97.069908 s**, one event.
- CLI workflow: **654 s**, whole-second Bash timer, distinct from native search and total job time.
- Peak sampled summed process-tree RSS: **14,267,208 KiB (13.61 GiB)**.
- Ran on **octopus10**, one CPU, requested 64 GiB, two-hour limit, no requeue/restarts, sampled 48 GiB stop threshold at 100 ms intervals. **Slurm did not establish a hard 64 GiB memory limit.** Sampling can miss spikes and double-count shared pages.
- The other two guides were not run separately; there is **no measured whole-assembly three-single-job versus batch speedup**.

Authoritative batch report and small evidence: [batch_2871284/REPORT.md](../research/pam_anchored_poc/optimized_pilot/batch_2871284/REPORT.md), [audit.json](../research/pam_anchored_poc/optimized_pilot/batch_2871284/audit.json), [per-guide-metrics.json](../research/pam_anchored_poc/optimized_pilot/batch_2871284/per-guide-metrics.json), and [pilot-comparison.json](../research/pam_anchored_poc/optimized_pilot/batch_2871284/pilot-comparison.json).

## Reuse, artifacts and reproducibility

Do not rebuild existing HPRC indexes unnecessarily. Original indexes are under `/moosefs/salehi/HPRCv2_CRISPRapido/indices/`. The verified imported bundle is:

```
/export/local/home/salehi/crisprapido-hprc-import-hg00735-20261105-v1/
```

It contains `reference.fa`, the existing Columba index components and `index.crispr-binding.json`; it was imported **without rebuilding the FM-index**. Preserve mandatory reference/index verification, original N-containing reference bytes and fail-closed behavior.

Large outputs remain outside Git under `/export/local/home/salehi/crisprapido-resource-pilot-results/`:

- `hg00735-pilot-2871278/`
- `hg00735-optimized-2871283/`
- `hg00735-batch-2871284/`
- Corresponding `monitor-2871278/`, `monitor-optimized-2871283/`, `monitor-batch-2871284/` directories retain sampled resource logs/summaries.

Expected PAF SHA256:

- Both standalone PAFs and the batch pilot subset: `a28d0574a915a57d431851ab681b0602b8f5fb1740d1c8efeb6e3b4a378f23be`.
- Full batch PAF: `a038919786424e8312eebe5eedb64f5c9dd8bee2074615823da577d7987d36a4`.

The isolated CLI was built from `f7de425` at `/export/local/home/salehi/crisprapido-pilot-rebuild-f7de425-v1/target/release/crisprapido`. Measurement-only runner/native instrumentation is separately documented and checksum-pinned; production binaries were not modified. Preserve its archived source because the executable embeds the runner path.

Batch reproducibility files live in `research/pam_anchored_poc/optimized_pilot/`: `batch_guides.fa`, `batch_targets.json`, `batch.env`, `submit_batch.slurm`, `batch_workload.sh`, `runtime.sha256`, `audit_batch.py`, and `extract_batch_metrics.py`. **Use `batch.env`, not just `optimized.env`, for the per-guide instrumented runtime.** Submission is not needed now and requires new approval plus fresh worker/path/capacity/integrity checks for any future run. The historical preparation [report](../research/pam_anchored_poc/optimized_pilot/REPORT.md) is retained with a current-status notice.

## Important commits

| Commit | Milestone |
|---|---|
| `8dcc2ad` | Optional integrated PAM-anchored backend; existing backend retained |
| `671d0b6` | Verified import of existing Columba indexes without rebuilding |
| `53a3b88`, `02bacec` | Bounded large-coordinate and repetitive HG00735 probes |
| `1d3ed09` | Portable timing used for successful monitored pilot |
| `f7de425` | Verification/reference-loading overhead reduction, unchanged search/scoring |
| `bfbb8ac` | Optimized comparison and audited batch preparation |
| `dfd1636` | Completed and audited three-guide HG00735 batch |
| `252a6bb` | CFD bound feasibility study, independent bounded oracles and model tests |
| `178de00` | Optional downstream CFD filtering and isolated native pruning, not enabled in the integrated backend |
| `4756527` | Limited-domain isolated CRISPR-Bulge postprocessor and official-model validation |
| `90d958a` | Evaluated saved 21,447-row HG00735 PAF; initial CLI deferral, subsequently superseded below |
| `6fa96d3` | Optional live/saved-file CLI scoring, rich N/A annotations and pinned reference-offset reuse |
| `ce52c15` | Complete bounded workflow validation, fresh/offline installation evidence and research reproduction setup |

## Remaining limitations and proposed next step

- Emitted-site audits establish eligibility and score correctness for reported alignments, **not absence of missed sites**. Agreement between implementations is not an independent recall oracle. Bounded oracle recovery does not prove arbitrary genome-wide completeness.
- The experimental backend supports one thread, `b=0/1`, `z=0/1/2`, resolved A/C/G/T guides and A/C/G/T/N reference FASTA. Unsupported CLI modes fail explicitly; the existing backend remains available. The native format is the tested 32-bit vanilla, sparse-SA-4 layout, not an unrestricted index API.
- Different valid representatives can yield different CFD scores; neither backend's representative is assumed biologically superior. CFD is not a calibrated cleavage probability. Missing/unsupported scores must not silently become zero.
- One batch amortizes setup, but every invocation still copies the complete bundle and performs repeated full checks. There is no immutable verified cross-job cache. Unsigned manifests/read-only file modes do not protect against malicious same-user mutate-and-restore races; trusted inputs/process isolation are required. Output buffering and large/repetitive batches remain resource risks.
- **Integrated CFD-guided pruning remains disabled.** An isolated native pruning prototype exists. CRISPR-Bulge is now an experimental opt-in Rust CLI scoring mode using the limited-domain batched postprocessor; see the current integration section below. Wider guide/assembly performance and extreme-repeat behavior remain to be studied; three guides do not establish general scalability.

### Completed bounded CFD feasibility investigation

See [CFD feasibility report](../research/pam_anchored_poc/cfd_feasibility/REPORT.md) and its scripts, tests and small retained evidence. No production/search/scoring code or CLI changed; no whole-genome job was run. A PAM-proximal column-product upper bound is safe for defined scores under validated finite `[0,1]` factors, but branches that may yield unavailable scores must survive. Pure-D projection retains the last 20 columns; mixed I/D for a 20-base guide is unavailable. Reference Ns require conservative fallback or separately validated N-aware interval handling.

The independent forward path oracle matched unchanged native representatives in **25 tiny-index configurations** (both strands, all 16 PAM suffixes, mismatch limits 5/8/20). **125 model threshold comparisons** had no lost qualifying/unavailable paths or final anchor/representative differences after unpruned refinement. **15,478** distinct oracle scoring inputs agreed with Rust (4,861 defined, 10,617 unavailable). Nine study tests passed. These are finite synthetic checks, not whole-genome recall evidence. The forward oracle and reverse model share replay/scoring helpers; the native/probe comparisons add separate implementation checks.

At threshold `.1`, resolved-panel model states fell **70.156% without bulges**, but only **30.046% at m5/b1/z2/f.75** because mixed-gap possibilities must survive. These are Python occurrence-set model counts, not native FM-index speedups. The conservative N-containing-reference fallback gave **zero state reduction** and would do so on HG00735_pat. Naïve pruning changed representatives/membership; any future implementation must preserve the current comparator, for example with candidate discovery followed by unpruned anchor-local refinement.

### Implemented downstream filtering and isolated native pruning

The user approved canonical-representative threshold semantics. Fresh development CLIs now accept optional **`--min-cfd SCORE`**, finite in `[0,1]`, disabled by default. The threshold is applied after the existing representative is selected, using the formatter's full-precision CFD value once. Defined scores below the threshold are removed; unavailable scores are kept without promoting another alignment. Omitted `cf:f` still denotes unavailable scoring, not zero. Original CLI options/defaults and higher mismatch limits remain unchanged; integrated searching remains unpruned with or without this flag.

See the [implementation/validation report](../research/pam_anchored_poc/cfd_pruning/REPORT.md), code and small retained evidence. The separate native executable performs real FM-index pruning and unpruned anchor-local reconstruction. It reuses unchanged native search/emission for controls, rounds bounds conservatively, preserves mixed-gap possibilities, and falls back on references containing any N, non-20-base guides, zero threshold or unsupported table contracts. Its distinct protocol is rejected by the integrated adapter; **do not set it as a production `--pam-search-bin`**.

All **72 A/B/C comparisons** (18 configurations × thresholds `0,.1,.2,.5`) matched integrated unpruned+filter, isolated pruned, and exhaustive direct-string oracle final anchors/representatives. Default output was byte-identical to the preserved prior CLI in all 18 cases. Additional native/oracle checks covered the other 15 PAM suffixes at four thresholds, custom-table/non-20 fallbacks, and binding rejection. Both strands, mismatch+bulge cases, two-base mixed events and reference Ns were tested. All 32 unavailable canonical mixed-gap representatives in the resolved m0/z2 panel survived every threshold. **164 Rust tests, 5 native/CLI regressions and 9 prior model regressions passed**; one pre-existing external-Python CFD test remained ignored.

Actual native measurements (three repetitions per mode; threshold `.1`, m5/b1/z2/f.75): CHM13/GRCh38 **50 kb** windows explored about 51% fewer states and had **22.015% / 20.972% lower native algorithm times**, including reconstruction/filtering/TSV construction but excluding loading, verification and full CLI setup. Dense synthetic output explored 31.454% fewer states but was **5.766% slower** after reconstruction. N fallback reduced no states. These are small native benchmarks, **not modeled predictions, whole-workflow speedups or genome-wide guarantees**. Native peaks were about 7.3 MiB; process-tree sampling has the usual shared-page/spike limitations. The report gives exact times, states, memory and scopes.

Runtime root: `/export/local/home/salehi/crisprapido-cfd-pruning-dev-v1/`; final native `native-build-v4/native`, development CLI `target/release/crisprapido`, authoritative logs `validation-final-v4/`. Production binaries and existing indexes remain untouched. Reproduction uses a new output directory via `research/pam_anchored_poc/cfd_pruning/reproduce.sh`; no large genomic input/output is in Git.

**Stop point:** do not enable integrated or HPRC pruning without separate approval and further independent review. Next candidates are broader bounded/randomized oracle coverage, custom-table/numerical review, matched full-workflow profiling and completeness-preserving N-aware discovery. Preserve representative policy and higher mismatch options. Threshold-selected output is not exhaustive under mismatch/bulge settings alone; successful finite tests do not establish genome-wide recall.

### Initial limited CRISPR-Bulge postprocessor milestone (historical)

**The statements in this initial milestone about no core CLI option are superseded by the current CLI integration below.** Original results/limits are retained.

See [CRISPR-Bulge report](../research/crispr_bulge/REPORT.md), scripts, manifest and retained evidence. Official source is `OrensteinLab/CRISPR-Bulge`, pinned **`3eddcd5bfcaff00b2bdf29425116ec9756ace870`**; the selected upstream-recommended model is the five-component GRU-Embd (`c_2`) classification ensemble. Root license is MIT; attribution is retained, while weights, pickles, datasets and the dedicated dependency environment remain outside Git. Inference uses checked HDF5 models without unpickling; only original-reference validation uses the trusted pinned upstream pickle loader.

**Full compatibility with `m5/b1/z2` is not achieved.** Checkpoints accept 24 columns including PAM; a two-base DNA bulge requires 25 and the official encoder rejects it. Published GUIDE-seq preprocessing uses six mismatches without gaps, four with gaps, and one gapped nucleotide (not one arbitrary-length event). The adapter scores only resolved 20-base-guide NGG alignments within that conservative subset. Two-base RNA, mixed/multiple gaps, five-mismatch single-bulge cases, ambiguous reference/guide/PAM, non-NGG and other unsupported cases remain present with an unavailable reason. Some are numerically encodable but scientifically unvalidated, not hard model-shape failures. No cropping, CFD-style projection, realignment, retraining or other workaround was implemented; limitations were reported before implementation.

The isolated `research/crispr_bulge/score_paf.py` provides **`--scoring cfd|crispr-bulge|both`**, default **cfd**. CFD mode passes the existing PAF through byte-for-byte; `both` appends `cb:f` only when scored and `bu:Z` status, retaining existing CFD tags/rows/CIGARs/order. CB-only intentionally omits `cf` tags, although CFD was still calculated upstream. It never selects a different alignment or filters sites. Core Rust CLI/search/scorer code, defaults and optional `--min-cfd` remain unchanged; **there is no core Rust `--scoring` option yet**. Apply `--min-cfd` upstream if wanted, not to rounded PAF values or as a neural-score threshold.

Validation: **542 rows retained** (370 synthetic plus the existing 172-panel output), **385 scored / 157 unavailable**; retained-panel coverage was **43/172 scored**, with all 129 other rows retained. Default CFD output was byte-identical on both panels; all original fields/CFD tags in `both` were preserved. Independent reference/CIGAR extraction and exact categorical-feature parity were checked against unmodified official `ensemble_predict` on identical inputs. Maximum prediction difference was **0** at chunk 512, **2.0116567611694336e-7** at chunk 32 (385 records), within `1e-6`; actual annotated PAF scores matched within **4.990997060083657e-10**. Singleton and three repeated 1,024-row comparisons also passed. **10 stdlib-only regression tests passed**. This is implementation agreement, not biological accuracy or genome-wide recall evidence.

Measured median fresh-process CPU wall time for 1,024 supported repeated alignments: Rust CFD probe **0.008286481 s**, isolated CB worker **5.070746292 s**, official reference wrapper **6.322246979 s** (three repetitions; scopes include startup/loading/I/O, not matched kernels or full search). Maximum individual-child RSS was **7,052 / 562,532 / 650,904 KiB** respectively. Model import/loading dominate cold calls; five models load once per invocation and are reused across record chunks/guides. The report distinguishes exact worker peaks from sampled summed tree RSS and records numerical/batch variation.

Runtime paths: upstream/models `/export/local/home/salehi/crisprapido-crispr-bulge-review-v1/`; isolated Python 3.11.15/TF-CPU 2.12.0 environment `/export/local/home/salehi/crisprapido-bulge-env-v1/`; raw results and preserved small legacy inputs `/export/local/home/salehi/crisprapido-bulge-validation-v1/`. The 172-panel search was not rerun; no whole-genome scoring/query or Slurm job occurred.

**Stop point:** review limited coverage/availability semantics before any core integration or extension outside the documented domain. Complete two-base-DNA support needs an appropriate model/data strategy, not truncating alignments. Biological validation, larger bounded batches, licensing/packaging review and full-workflow profiling remain. Do not transfer CFD pruning bounds to this nonlinear scorer, claim either score is a calibrated cleavage probability, or promote pruning automatically. Manuscript updates and full HPRC analysis still require separate decisions.

### Evaluated CRISPR-Bulge on the saved whole-assembly output (historical milestone)

The deferral below records the earlier `90d958a` decision; the user subsequently approved the optional CLI work documented next.

See [saved HG00735_pat evaluation](../research/crispr_bulge/ASSEMBLY_EVALUATION.md), aggregate evidence and `evaluate_assembly.py`. Starting from `4756527`, applied the **unchanged** batched `both` postprocessor to the existing job **2871283** PAF, not a new search. All **21,447 original records and CFD scores** were preserved; input SHA256 remains `a28d0574a915a57d431851ab681b0602b8f5fb1740d1c8efeb6e3b4a378f23be`.

- **4,526 CB scores / both scores (21.10%)**, **16,921 unavailable**. Supported: 1,098 ungapped, 1,080 DNA1, 2,348 RNA1.
- Unavailable: **8,176 two-base DNA bulges** exceed 24 columns; **7,624 DNA1/five-mismatch** rows exceed conservative with-gap mismatch policy; **1,121 RNA2** rows exceed one-gapped-nucleotide policy. Zero malformed CIGAR/mixed/multiple-gap or ambiguity statuses. No unsupported prediction was invented.
- Independent samtools extraction and separate CIGAR reconstruction for **93 stratified/tail-enriched rows**; official feature hashes matched, official prediction maximum absolute difference **2.5598145825256324e-10**. **14 tests passed**. Implementation agreement, not biological accuracy or recall.
- Paired Spearman **0.4281**, Pearson **0.3447**; score distributions differ markedly, especially gapped versus ungapped. No direct calibration, common thresholds or interchangeability claim.
- One measured CPU annotation: **91.5757 s**, including **81.7841 s reference scan**, **6.0402 s neural worker**; peak individual RSS **561,608 KiB**, sampled summed tree **617,108 KiB**. Not search runtime or a repeated median; sampling caveats remain.
- External raw artifacts: `/export/local/home/salehi/crisprapido-bulge-hg00735-2871283-v1/`, annotated output `both/output.paf`, authoritative audit `evaluation-final.json`, independent check `independent-final/`. Models/reference/original results untouched; no large inputs/outputs committed.

**Readiness/stop point:** this is usable as an isolated research postprocessor and technically supports a future experimental opt-in wrapper, **not unrestricted production CLI promotion**. Core integration remains deferred pending agreement on partial availability, dependency/licensing packaging and CLI lifecycle/error tests. Full default-domain coverage remains unsupported. Findings are reported before any integration; no main CLI/search/default/binary changes, no new search or Slurm job, no main merge. Seek explicit user direction before proceeding.

### Current optional CLI scoring integration

See [CLI_INTEGRATION.md](../research/crispr_bulge/CLI_INTEGRATION.md), [small results](../research/crispr_bulge/results/cli-integration/summary.json), and `src/scoring.rs`. User approved partial-domain CLI integration starting from **`90d958a`**. Development CLI now supports **`--scoring cfd|crispr-bulge|both`**, default CFD unchanged. Neural modes spool original reporting once per invocation and call the existing five-model batched worker once, including multiple guides. They preserve all sites and existing CFD tags, append neural tags and write `annotations.tsv` with sequences/alignment/PAM/context plus explicit **N/A** values/reasons. CB-only intentionally retains existing CFD, like `both`; mode selection is recorded in metadata. No score fusion or implicit threshold is applied.

- **`--score-paf PATH`** annotates existing results without any search/CFD recalculation. Neural modes require a NEW `--scoring-output-dir`, absolute `--crispr-bulge-python` and existing `--crispr-bulge-model-root`; source-bundle/script requirements are documented. Original PAF is read only. Default saved-file CFD is byte pass-through with no model/reference scan. `--min-cfd` and search executable/index/candidate options are rejected with `--score-paf`; during search the existing threshold is still explicitly opt-in/disabled by default.
- Original default PAF bytes, both search algorithms, representative selection, CFD formulas/matrices, mismatch/bulge defaults, optional threshold semantics and integrated-unpruned policy are unchanged. Only reporting destination/lifecycle and a test fixture in verification changed. No production binary, existing index or original PAF was modified; only an external development target was built.
- Existing 21,447 HG00735 job 2871283 rows were checked in all three modes: all CFD tags/rows preserved, **4,526 CB / 16,921 N/A** in each neural mode. Same unavailable reasons: **8,176 DNA2 over-length, 7,624 DNA1/five-mismatch domain exclusions, 1,121 RNA2 unsupported**. Both neural PAFs equal the previously evaluated annotated PAF byte-for-byte.
- Independent samtools extraction/separate CIGAR replay plus **unmodified official ensemble inference for ALL 4,526 supported rows**, not just the prior sample. Feature hashes exact; max PAF-score difference **2.5598145825256324e-10** in each mode, within `1e-6`. This is implementation agreement, not biological accuracy/calibration/recall.
- Optional `--scoring-reference-index` reuses explicit trusted offset manifests prepared by `research/crispr_bulge/reference_index.py`. Requires **both trusted index and reference SHA256 pins**. Full reference hashes still occur before/after scoring and index digest is rechecked; no automatic unsigned cache trust, no weakening search/index verification. Fallback remains the full original Python FASTA scan. Plain regularly wrapped FASTA required in neural modes.
- One fresh-process saved-result measurement per mode on octopus01: **12.929960 s CB / 13.781680 s both**, worker ~6.1 s, initial pinned-index/hash check ~2.23 s; max individual RSS **573,596 / 577,588 KiB**, sampled summed tree **629,832 / 634,612 KiB**. One-time index preparation **75.4562 s** (separate, not included in annotation timing). Prior full-scan annotation was 91.576 s; this is amortized reference-validation improvement, not a search speedup or repeated median. CFD pass-through **0.007562 s** is not CFD calculation time.
- **164 Rust tests passed** (one pre-existing ignored), **20 configured Python tests passed**, including tiny live single/multi-guide reporting, prior default compatibility, all modes, unsupported DNA2/RNA2/mixed/N/A/zero, empty inputs, failure/timeout, model reuse and pinned-index identity checks.
- External root `/export/local/home/salehi/crisprapido-scoring-cli-v1/`: development binary `target/release/crisprapido`, trusted offsets `reference-index.json`, complete scoring bundles under `saved-results/{crispr-bulge,both}/scoring/`, full official checks `saved-results/official-all/`. No genomic/model outputs committed.

**Stop point/limitations:** optional experimental CLI scoring is implemented, not a production packaging release. Full default-domain model coverage remains unsupported; no fabricated scores, truncation, retraining or re-alignment. Source-bundle/TF dependency and licensing/packaging review remain; full pre/post reference hashing and disk spooling still cost resources. No hard CLI memory cap/global timeout, GPU/concurrent-writer protection, broad biological validation or genome-wide recall guarantee. Do not transfer CFD pruning bounds, turn on thresholds/pruning automatically, change defaults, launch a new whole-genome search or merge main without separate direction.

### Complete bounded workflow and reproducible research setup

Current continuation starts from **`6fa96d3`**. See [end-to-end/reproducibility report](../research/reproducibility/REPORT.md), [installation guide](RESEARCH_INSTALLATION.md), [workflow results](../research/reproducibility/results/workflow-summary.json) and [installation provenance](../research/reproducibility/results/installation.json). **No core algorithms/scoring/defaults/Cargo files changed.** Two stale compatibility tests now create proper bound-index fixtures instead of raw unbound indexes; native integrity enforcement was not weakened.

- Actual **PAM-anchored FM-index search → canonical reporting → CFD → batched neural scoring → final output**, not saved-PAF pass-through. All three scoring modes recovered **566/566 resolved**, **570/570 N-containing**, **80/80 mixed-only m0 control**, and **172/172 pinned retained-panel** expected anchors; **zero missing/unexpected/duplicate anchors**. Default CFD PAF was byte-identical to the preserved preintegration CLI in every panel. All CFD scores/availability were independently replayed against tables; all rows/CFD tags/coordinates/CIGAR/order survived neural modes.
- CB defined/unavailable by panel: **550/16**, **552/18**, **32/48**, **43/129**. N/A reasons explicitly retained in TSV (DNA2, RNA2, mixed, ambiguity and mismatch-domain exclusions). Default search limits remain m5/b1/z2/f.75/NGG, with explicit one-thread tests; only the mixed-only control explicitly requests m0. No score threshold/pruning flag is passed.
- **1,177 supported panel-row predictions** matched unmodified official ensemble inference within `1e-6`; maximum absolute difference **1.789291992349007e-7**, allowing known batch/float variation. Panels overlap intentionally; this is not a count of distinct biological sites. Finite oracle recovery/implementation agreement does not prove biological accuracy/calibration or arbitrary genome-wide recall.
- Identical **13,956-base / 606-contig N panel**, three fresh complete CLI calls per mode: median **0.214627 s CFD / 5.020807 s CB-requested / 4.970432 s both**. Maximum individual RSS **15,056 / 538,428 / 536,748 KiB**; sampled summed tree **51,008 / 589,128 / 587,724 KiB**. Includes search/snapshot verification/loading/reporting/scoring/publication, not one-time index construction or independent validation. CPU/shared-host/warm-cache/sampling caveats remain; tiny-run setup dominates, no genome-scale speedup claim. CB-requested mode still retains upstream CFD like both.
- Clean `git archive 6fa96d3` source and **empty Cargo target** built with Rust/Cargo **1.94.1**, coherent GCC/G++ **12.3.0**, fresh WFA2 v2.3 library at **df3dce4**, fresh pinned Columba **edcc7ce/2.0.3** libraries/builder/unpruned native. Trusted Cargo/dependency source caches reused, not target objects. Guix/system header/linker/-lrt pitfalls resolved and documented, not a portable container certification.
- New Python **3.11.15** inference venv (TF-CPU **2.12.0**, NumPy **1.23.5**, Keras **2.12.0**, h5py **3.16.0**) passed full execution. A second new venv installed the **46-wheel hash lock offline**, had identical frozen packages, and produced a byte-identical complete retained172 PAF. New lock: `research/reproducibility/requirements-inference-linux-py311.lock` (Linux x86_64/CPython 3.11 only). Standard-library Python is required for PAM backend integrity/snapshot helpers even in CFD mode; no neural environment is needed for CFD.
- **Installation blocker:** actual network model-fetch helper attempts failed (Git 2.41 internal partial-fetch abort; Git 2.20 HTTP 400/RPC failure). TLS was not disabled. Execution used the five existing reviewed weights copied into a new external directory and SHA256-verified, **not a successful fresh network model install**. The official-reference environment was reused. Reliable model distribution, rights/security review and portable packaging remain release gates.
- **164 Rust tests passed**, one pre-existing ignored; **20 scoring + 4 new validator + 9 prior offline CFD model tests passed**. Also **96 PAM/option oracle comparisons**, five invalid-option/two invalid-guide checks, binding/tamper/PATH-spoof checks and both mixed-gap strands passed. The final bounded reproduction wrapper itself completed successfully.
- External root **`/export/local/home/salehi/crisprapido-research-workflow-v1/`**. Authoritative all-pass run **`final/`**, development binary **`target-final/release/crisprapido`**, clean source **`source/`**, fresh native **`pam-search`**, builder **`columba-build-v2/columba_build`**, model root **`models-verified/`**, venvs **`inference-env/`** and **`offline-env/`**, offline check **`offline-smoke/`**. Retain sources because helpers are path-embedded. Earlier exploratory/failed installation/test attempts are preserved but are not the authoritative run. No genomic/model/wheel/environment outputs committed.

**Release stop point:** bounded research workflow passes; no production release, main merge, default change, new whole-genome experiment or CFD-pruning enablement is authorized. Resolve reliable model setup/rights, portable build/relocation/environment packaging, partial-domain scientific policy, broader independent recall/biological benchmarking and resource/concurrent-writer hardening before release consideration. Large live-search/neural memory overlap and repeat-heavy outputs remain unprofiled; hashes assume trusted preparation and cannot defeat malicious mutate-and-restore races. Existing default Columba/WFA2 completeness limitations remain unchanged.

### Systematic search-mode concordance and manuscript scaffold

Current continuation starts from **`ce52c15`**. See [search concordance](../research/search_concordance/REPORT.md), [implementation taxonomy](../research/search_concordance/IMPLEMENTATIONS.md), [complete cases](../research/search_concordance/results/cases.json), [all raw-record audits](../research/search_concordance/results/records.tsv) and [manuscript outline](../manuscript/OUTLINE.md). **No algorithms/defaults/production binaries changed; no large experiment or main merge.**

- The five requested names are **not five independent engines**. “Original” can mean historical standalone window/WFA2 (`46cf652`), the window branch in older `f75f82aa` figure benchmarks, or Sassy-labelled original on divergent main (`cfb1587`). Current window/WFA2, automatic Columba+WFA2 and PAM-aware Columba are three complete configurations of the same attested release executable (built 6fa96d3, core identical to ce52c15). WFA2 is an aligner/verifier shared by two families; native PAM-aware traversal performs no WFA2 query verification. Sassy uses a different branch/interface; no merge was attempted.
- New **700-base / 28-contig adversarial controls**, exact guide/ref/PAM/strand and equivalent modern limits: at m3/b1/z1/f.75/NGG, **window 20/26, Columba+WFA2 26/26, native 26/26**, shared by all current modes **20**. At m3/b1/z2: **24/30, 28/30, 30/30**, shared **24**. Separate 46-base m0/z2 mixed-only control: **0/2, 0/2, 2/2**, shared **0**. Zero unexpected eligible/duplicate standard anchors in current pipelines; every emitted current alignment/defined CFD independently valid. Exact shared keys are committed.
- **All 70 anchor cases / 254 first-replicate raw records** audited and linked; repeats also checked. Current window misses are localized by preferred-alignment probes: terminal-D/PAM boundary, invalid preferred gap versus valid constrained alternative, and one optimum hiding another PAM in a dense window. Four panel-anchor Columba+WFA2 omissions have valid oracle candidates explicitly injected and still fail: unfiltered anchored `4X` versus eligible mixed/three-X, or `1X` versus eligible mixed/m0. Raw candidates already contain the exact interval for both m0 mixed cases. This is verification/fallback constrained-alignment loss, not simply absent raw discovery. No fix was made under the no-algorithm-change scope.
- Shared anchor `D5_2+`, boundary22: two independently valid representatives `[1,22) 1X2=1D2X15=` vs `[3,22) 1I1=3X15=` give CFD **0.2700 vs 0.1736**. Existing separate 204-site audit reused: **54 valid differing pairs, 41 CFD differences, 108 independent/Rust scores agreed**. Site agreement is not score/alignment identity or biological preference.
- Historical standalone WFA2 and Sassy are **diagnostics, not equal Venn/ranking sets**. Historical eligible-witness counts **11/13/0**; Sassy **16/16/0**, under separate schema/validity labels. Sassy exports eleven mandatory-like fields (missing PAF block length) and fixed target spans for gaps; z2 radius is `m+b`, not `m+b*z`. No coordinate/CIGAR repair was hidden in normalization. Exact old source-to-build attestation unavailable; Sassy release rebuild blocked by missing pinned Git cache after using compatible Cargo. Do not call absent valid exports intrinsic engine recall failures. Historical/Sassy causal attribution remains TODO where reporting obscures it.
- Matched common controls, three fresh calls/current mode: medians **0.004152 s window / 0.016266 s Columba+WFA2 / 0.114388 s native**; max individual RSS **7,076 / 20,376 / 15,100 KiB**. Native is not fastest on this tiny input; setup/integrity work differs and window misses sites. Historical/debug-path timings retained but **no fair five-way speed rank**. Earlier 8 kb/64-guide native-vs-Columba timings are reused separately, not rerun/pooled into a universal speedup. No chromosome/genome search occurred.
- **Six new audit tests passed** (anchors/strand, labelled legacy schema, no geometry repair, mixed witness/N/A, valid injected SAM tags/orientation, duplicates); **4 existing auditor + 9 CFD-model + 20 scoring Python tests** also reran and passed. Prior 164-pass Rust evidence is reused because core is unchanged, not claimed as a fresh rebuild. Executable/source hashes checked before/after. External root **`/export/local/home/salehi/crisprapido-concordance-v1/`**, authoritative **`final-v4/`**, earlier development/failed-tag probes preserved. Existing current/historical/Sassy binaries were read only; no production build replacement.
- Manuscript proposals: [outline/abstract structure](../manuscript/OUTLINE.md), [six figure topics](../manuscript/FIGURES.md), [claim/evidence/TODO ledger](../manuscript/EVIDENCE.md), [journal fit](../manuscript/JOURNALS.md). Existing figure/legend materials reviewed and preserved. No full prose draft was located in inspected material paths. Old manuscript whole-genome/chr figures predate strict PAM fixes (including non-GG outputs) and need re-audit before current-primary reuse. No population-level HPRC, broad recall, neural biological accuracy or universal speed claim is invented.

**Next stop point:** choose the intended paper baseline; obtain source/profile-attested historical/Sassy comparators and separately review exporter/semantic reconciliation; expand held-out bounded validation and re-audit old raw figures before any new large experiment proposal. Current native finite-control recovery is not a general completeness proof. Manuscript figures/outline are proposals with TODOs. Journal assessment is conditional (Bioinformatics methods novelty/rigor, NAR GAB genomics application, Bioinformatics Advances research-software fit); exact current submission policies/citations/availability remain to verify. No release, main merge, algorithm/default change, pruning enablement or new chromosome/whole-genome run is authorized.

## Working rules for future sessions

Check branch, commit and status first; preserve existing experiments, frozen branches, deterministic selection, production binaries and configurable defaults. Never merge into main or submit jobs without approval. Add regression tests for code fixes and run relevant tests. Commit/push completed source and small documentation/evidence only; exclude credentials, FASTA, indexes and large outputs. The latest continuation adds bounded normalized search-mode concordance/forensics and a manuscript/evidence/figure/journal scaffold. It changes no algorithm/default/production binary, enables no pruning and runs no large experiment. Use the taxonomy, full case ledger and explicit comparability/TODO gates before claiming a five-engine comparison or planning larger work.
