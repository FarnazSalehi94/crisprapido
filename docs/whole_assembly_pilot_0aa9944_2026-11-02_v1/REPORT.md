# Controlled whole-assembly pilot — 0aa9944

**Status: complete for this bounded pilot.** All four jobs and checks passed.
The frozen preflight gate does not pass: independent review and unresolved
sensitivity and representativeness risks remain.

## Frozen implementation and sources

The tested algorithm is implementation commit
`0aa9944e2eae4b7f8286415ee5c3fb3667cf6097` with the optimized frozen binary
SHA256 `c37f4ded6e7c6afb3a7d22eccc5566f291244731bf6872dbfbd0b4d47f53aa74`.
The dependency and scoring matrix digests are in
[the freeze manifest](../release_validation_0aa9944_2026-11-02_v1/freeze_manifest.json).
All 12 resolved assemblies and their 144 index files, FASTAs, `.fai` and `.gzi`
files must pass a full current SHA256/size scan against
[the recorded assembly manifest](../index_recovery_2026-10-29_v1/verified_assemblies.tsv)
and [component manifest](../index_recovery_2026-10-29_v1/verified_component_checksums.tsv)
before any search can start. This is an integrity comparison with the recorded
manifests, not an independent reconstruction of each index.

## Design

The first two frozen 100-guide entries are ELMOD2
`TCTAGCGTTTAACTTCTTCG` and MED6 `TGCCACCAATACCCTTTGGA`. Their existing
whole-assembly independent exact-match truths include, respectively, a minus
and a plus strand hit in each assembly. They are a convenience selection,
not a representative random sample. Each is run separately against complete
CHM13 and GRCh38 with `-p GG -m 4 -b 1 -z 1 -f 0.75 -t 4`, strict NGG,
Columba e=5 plus e=0, and the frozen scoring matrices. No search of the full
100-guide by 12-assembly workload is launched.

Raw data and logs are retained outside Git in
`/moosefs/salehi/HPRCv2_CRISPRapido/results/whole_assembly_resource_pilot_0aa9944_20261102_v1`.
A rate-limited 128 MiB/s, low-priority verification job (2870005) precedes
four low-priority, serialized Slurm jobs (2870006–2870009) using `afterok`.
Each search/audit job requests four CPUs, 64 GiB and four hours. These limits
apply to the combined search, replay and auditing work; search wall time and
GNU time maximum resident set size refer to the frozen search invocation.
CPU/I/O niceness and shared cluster activity can affect elapsed time. The
peak RSS of a single process is not a measurement of the combined peak RSS of
all children.

For each run the raw mapped e=5 and e=0 SAM row counts and the union of the
application's candidate keys are recorded. A separate pre-existing
instrumented build replays that union to count WFA2 calls; its output must
match the frozen binary's output as a multiset of full PAF rows. WFA call
counts are diagnostic replay measurements, not instrumentation of the frozen
binary itself. Every frozen output row is independently checked against
indexed FASTA, guide, PAM, CIGAR, limits and CFD matrices, and duplicate PAM
keys are rejected. The previously produced independent genome-wide exact
scanner truth is compared against the e=0 SAM and against the PAM anchors in
final audited output. This exact check does not establish mismatch/bulge
completeness. Prior independent bounded coverage checks are described in the
[preceding validation report](../release_validation_0aa9944_2026-11-02_v1/REPORT.md).

## Results

Full hash verification passed: 181 files, 144 index components from 12
assemblies, 149,003,337,963 bytes; job 2870005 completed successfully.

All four searches, diagnostic replays and audits completed successfully.

| Assembly | Guide | Search seconds | Max RSS KiB | Unique raw | WFA calls (replay) | Audited rows | Exact truth recovered |
|---|---|---:|---:|---:|---:|---:|---:|
| CHM13 | ELMOD2 | 349.13 | 11,693,052 | 248,954 | 27,392,485 | 3,027/3,027 | 1/1 |
| CHM13 | MED6 | 320.70 | 11,671,328 | 219,251 | 24,366,556 | 3,540/3,540 | 1/1 |
| GRCh38 | ELMOD2 | 450.58 | 11,664,660 | 315,144 | 34,291,273 | 2,969/2,969 | 1/1 |
| GRCh38 | MED6 | 329.68 | 11,617,916 | 224,598 | 24,952,210 | 3,531/3,531 | 1/1 |

All four replay PAF row multisets matched their frozen search; the independent
auditor passed 13,067/13,067 rows, with no invalid or duplicate PAM keys.
This includes 6,652 plus, 6,415 minus and 12,794 bulged rows. That gives 100%
emitted-row audit precision for this pilot, **not** a whole-assembly estimate
of mismatch/bulge recall. The e=0 results recovered all four independent
whole-assembly exact sites, and all four corresponding PAM anchors appeared
in final output. Exact raw-locus recall and precision were 4/4 and 4/4; there
were no missed or unexpected exact loci. Exact truth is from the pre-existing
independent [CHM13](../broader_coverage_1eb2f40_2026-10-30/CHM13/oracle.jsonl)
and [GRCh38](../broader_coverage_1eb2f40_2026-10-30/GRCh38/oracle.jsonl)
scanner outputs, filtered to these frozen guide identities. A same-anchor
comparison does not assert that all possible alternative alignments were
reported. The auditor and scanner are independent implementations, not an
independent human review.

Across runs there were 1,007,947 unique raw candidate keys and 111,002,524
WFA2 calls in instrumented replay. The e=5 raw mapped counts equal those keys
in each run, while each e=0 SAM had one mapped row already present in e=5.
Search-time warnings for candidates rejected by WFA2 are retained in the
external logs; all four search processes and jobs exited successfully. Slurm
reported no timeout or out-of-memory event. Full job durations, including
replay and audit, were 11:17, 10:19, 14:33 and 10:33 respectively. Slurm
accounting did not provide reliable CPU or aggregate memory figures here.

The [previous held-out comparisons](../release_validation_0aa9944_2026-11-02_v1/REPORT.md)
found 212/212 synthetic and 18/18 real-window PAM anchors with 100% precision
on those finite oracle panels, including mismatches, both gap directions and
both strands across the fixtures. Those checks exercise the same frozen
implementation and settings; they do not exhaustively enumerate these four
whole assemblies, and their results must not be interpreted as whole-assembly
bulge sensitivity. The earlier [broader synthetic panel](../bounded_coverage_correction_2026-10-30_v1/REPORT.md)
also included these same two guide sequences (named `frozen1` and `frozen2`),
and recovered 237/237 anchors across the complete multi-guide panel; the
subsequent optimization retained the validated anchors as documented in the
[performance report](../performance_validation_a817ab5_2026-10-31_v1/REPORT.md). Original strict interval and category metrics remain in
those earlier reports and artifacts.

## Resource estimate and limitations

The four frozen search times average 362.52 seconds (range 320.70–450.58)
for an individual guide/assembly search with four CPUs requested. A direct,
illustrative 1,200-search multiplication gives **120.8 serial search hours**;
using the observed per-search range gives **106.9–150.2 hours**. At four
allocated CPUs per serial search, the mean extrapolation corresponds to about
483 allocated core-hours, not measured CPU consumption. The searches here
are separate processes; this is not a prediction for a batched execution
plan, a confidence interval, a worst-case bound, or an estimate including
index construction, diagnostic replays, independent audits, staging, queue
time or failures. The guides and assemblies are not a representative random
sample; repeats, high-hit guides, other assemblies, I/O, contention and
threading may change costs substantially. The low-priority I/O and CPU settings
also affect observed wall time.

The measured maximum process RSS ranged from 11,617,916 to 11,693,052 KiB
(about 11.1–11.2 GiB). This includes the frozen timed command's child-process
maximum, but is **not** the simultaneous aggregate peak across all processes.
Each pilot allocation requested 64 GiB; the pilot does not establish a safe
minimum allocation or whole-plan concurrency envelope. The e=5 candidates
varied by 44% between the lowest and highest runs, and the ELMOD2 GRCh38 search
took about 29% longer than its CHM13 counterpart. Output volume and verification
warnings grow with candidate density. The original index/FASTA paths remain
necessary to reproduce this test; no large indexes, references or results
were committed.

Whole-assembly exact testing and output auditing do not prove Columba discovers
every eligible mismatch or bulge locus. The two-guide panel cannot establish
worst-case memory or runtime, and no independent reviewer has yet examined the
evidence. The frozen preflight gate therefore **does not pass**. No release or
full-run authorization is implied by the successful pilot.

## Test and provenance checks

All 157 standard Rust tests and the separately executed pinned standalone CFD
test passed with Rust 1.85.1 and the pinned WFA library. The initial attempt on
the login node failed because its inherited library path selected an
incompatible glibc; rerunning with `LIBRARY_PATH` unset succeeded. All 52
Python tests passed, including the checksum-tampering test. The pilot runner
was also exercised end-to-end against a small existing independent indexed
fixture before the whole-assembly jobs. All Slurm job IDs, commands, stage
outputs and raw files are under the versioned external directory; small
summaries, exact checks, audit provenance, GNU time outputs and the full
[verification checks](verification_checks.jsonl) are retained here.
