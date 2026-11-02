# Controlled whole-assembly pilot of m5/b1/z2 defaults — stopped early

**Status: stopped for excessive candidate volume.** Only CHM13–ELMOD2 was
started. Columba's e=7 discovery completed, producing 17,653,389 reported
matches and a 1,944,321,145-byte SAM; the prior e=5 pilot for the same guide
and assembly produced 248,954 mapped candidates and a 26,705,890-byte SAM.
The increases are 70.9× and 72.8× respectively. Because this approached the
2 GiB file limit and implied a very large verification workload, the job was
cancelled at 5:59 while CRISPRapido was still processing candidates. No
other guide or assembly job was launched. The limit was not actually reached,
and there was no reported out-of-memory or program crash before cancellation.

This is **not a completed whole-genome validation**, and these candidate
numbers are not final eligible off-target counts. The e=7 search is broader
than e=5 because `m + b*z` rises from 5 to 7. Many extra raw candidates may
fail WFA2, match-fraction or PAM eligibility; without completed verification
we cannot establish the number or identity of additional final PAM anchors,
or whether any previous anchor is absent. Columba logs report the raw counts;
the 17,653,389 rows have not been independently deduplicated.

## Setup and observed measurements

The tested development code is `60f9d6c` (production source unchanged from
`49f20fc`), using the separately compiled release executable SHA256
`cb10e5d6be970bb6ef90aac4f2fe5ab035ae2b4d4ce643f28ec32f3895d7e42d`.
It used the same CHM13 reference, resolved Columba index, ELMOD2 guide, scoring
matrices and four-thread configuration as the [previous pilot](../whole_assembly_pilot_0aa9944_2026-11-02_v1/REPORT.md),
with `-p GG -m 5 -b 1 -z 2 -f 0.75 -t 4`. The e=7 and e=0 Columba phases
both completed. Slurm job 2870108 requested one hour, 64 GiB and four CPUs at
low priority, with a 2 GiB per-file ceiling. It was manually cancelled because
of the candidate volume; neither the ceiling nor the Slurm time limit was
hit. Large artifacts remain outside Git at
`/export/local/home/salehi/crisprapido-whole-defaults-60f9d6c-20261102-v1`.

| Measurement | Previous m4/b1/z1 | New m5/b1/z2 |
|---|---:|---:|
| Columba main edit bound | 5 | 7 |
| Columba logged matches | 248,954 | 17,653,389 |
| Main SAM size | 26,705,890 bytes | 1,944,321,145 bytes |
| Columba mapping phase | 11.24 s | 275.37 s |
| Entire Columba main phase | 15.46 s | 297.47 s |
| Complete search runtime | 349.13 s | unavailable: cancelled at 359 s |
| Exact raw e=0 site | 1/1 | 1/1 |
| Final PAM anchors | 3,027 | unavailable: no final output |
| Independently audited final anchors | 3,027/3,027 | unavailable: no final output |
| WFA2 calls | 27,392,485, instrumented replay | unavailable |

The independent exact scanner's one CHM13–ELMOD2 site was found in the new
e=0 SAM, with the reference and strict PAM rechecked (`exact_raw_chm13_elmod2.json`).
Its presence in a final alignment cannot be verified because the final PAF was
never emitted. The new SAM, e=0 result and logs were complete before
cancellation; the top-level `actual.paf` and GNU time output are empty and
must not be interpreted as a zero-site result. The new binary is
uninstrumented. Replaying millions of candidates merely to obtain a WFA2
counter would contradict the resource stop condition, so no such replay was
started.

The largest sampled Columba process RSS was 19,604,384 KiB (18.70 GiB); at
that sample the CRISPRapido parent used 3,047,976 KiB, for a sampled combined
21.60 GiB. Later a CRISPRapido sample reached 11,492,772 KiB while it was
processing the SAM. These are sampled RSS values, **not measured peaks** and
not necessarily the complete process-tree aggregate. Slurm accounting did not
provide reliable MaxRSS, and GNU time could not finish after cancellation.
The previous pilot's 11,693,052 KiB is a completed command's maximum-process
RSS; it is not directly comparable to the sampled combined figure. Raw
samples and Columba logs are included here.

## Audit and tests

The independent emitted-row auditor now accepts explicit mismatch, bulge and
fraction limits while retaining its original strict m4/b1/z1 defaults. Its
new-limits mode passed all 204 real-window fixture rows, including five
mismatches and two-base gaps, and tests verify those rows remain rejected
under the old defaults. It independently checks strict PAM, reference and
CIGAR, mismatch and gap limits, match fraction and CFD; the pilot wrapper also
checks duplicate guide/strand/PAM keys and compares old/new anchor sets.
Those checks were **not run on a completed new whole-assembly output**, because
there was none. Existing Rust and Python tests pass. No production binary,
validated branch, prior benchmark, main branch or full 100×12 search was
modified or launched.

## Remaining work

The other three guide–assembly pairs, full final-output auditing, final exact
anchor verification and old/new locus comparisons remain incomplete. The
observed 70.9× raw candidate growth is a resource warning for this single
case, not a validated estimate of all guides or assemblies. A separate
correctness-preserving performance strategy or an explicitly bounded
alternative study is needed before another whole-assembly attempt; no
candidate cap or eligibility rule should be silently relaxed. Independent
review and the frozen release preflight remain outstanding.

## Subsequent candidate-efficiency study

The [bounded efficiency investigation](../candidate_efficiency_61a4cfc_2026-11-02_v1/REPORT.md)
found no duplicate e=7 keys in the saved SAM and evaluated two distinct PAM
screens. A disposable interval-level precheck preserved bounded oracle and
regional records while substantially reducing WFA2 work; it has not been
implemented in production, and it does not eliminate Columba's e=7 cost.
