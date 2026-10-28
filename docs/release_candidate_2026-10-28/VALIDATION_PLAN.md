# Corrected CRISPRapido release candidate and complete-search validation plan

Prepared 2026-10-28. **PLAN ONLY: the 100-guide / 12-assembly search has not run.**
No merge, release tag, production-binary replacement, or CRISPR-Bulge integration
is authorized by this plan. An optimized candidate binary is prepared separately.

## 1. Git safety and release scope

Candidate source commit: `9c2d3f53673fd135f0b4379812a9a7a4a7f20961`, on existing
branch `fix/pam-cfd-validation`. Its parent is the deployed local production
commit `f75f82aa8a5301235bbd12501933530cb6b282e4` (`columba-wfa2-cfd`).
The remote production branch is at `115dbe2e1a8bc516cdef8aa86e57d9240c1481e4`;
the local batch-guide commit was already one commit ahead before the PAM fix.

The repository is `/moosefs/raid5/salehi/crispr2/crispr-progress/columba-crispr-benchmark/crisprapido`;
the isolated worktree is `/export/local/home/salehi/crisprapido-pam-validation`.
`origin` is `FarnazSalehi94/crisprapido`, with Pinellolab configured as `upstream`.
Documentation stays on the existing fix branch and is pushed to that same fork.
The original dirty working files still match the SHA256 manifest from the earlier
backup. They include zero-byte `main.rs`/`reporting.rs` and must not be overwritten
by checking out or resetting the production directory. Their archive remains in
`analysis/PAM_fix_validation_2026-10-28/` in the analysis workspace.

**Merge stopped:** `main` (`cfb15874fdaef082e1a9cf2f83d30e1d10d4252a`) is not an
ancestor of this branch. At the candidate commit, main/fix are 72/37 commits
apart, with common ancestor `89e9682cf1f8dda36b89cc03059475e85e929189`.
A read-only `git merge-tree` simulation produces 28 conflict blocks in
`README.md`, `src/main.rs`, `src/cfd_score.rs`, and `tests/cfd_tests.rs`.
Main includes Sassy/ropebwt3 backends and `--allow-indels` that are absent from
the older deployed Columba/WFA2 lineage. Replacing main with this branch would
discard those features. Neither an automatic merge nor an `ours`/`theirs` conflict
resolution is safe. A future reconciliation needs a separate reviewed scope.
See `merge_review.txt` and `git_review.txt`; the ancestry counts are recorded at
the code candidate, before these documentation commits.

Relative to the intended **deployed Columba/WFA2** baseline, no CLI option was
removed; `--allow-noncanonical-pam` is added. Cargo.toml, Cargo.lock,
`src/columba.rs`, batch-equivalence tests, scoring matrices and existing benchmark
assets are unchanged. Multi-guide input, exact-candidate recovery, imported SAM,
WFA2 bulge verification, gzip FASTA input and CFD remain present. The preexisting
production DNA-bulge projection is retained. This is evidence of preservation of
the production lineage, not certification of every feature on main.

## 2. Regression review and release gates

The preceding audit passed 148 tests and independently verified 208 corrected
emitted rows. Review checks for this preparation:

- Three-base PAM bounds are checked before orientation; only the final oriented
  dinucleotide enters CFD. Incomplete/ambiguous PAMs are rejected in strict mode.
- `-p GG` is a strict NGG suffix predicate on both strands, applied before overlap
  selection. `--allow-noncanonical-pam` and `--no-filter` are excluded from the run.
- Already oriented negative PAMs are never reverse complemented a second time.
- Reference-only D does not advance query start or shift target start twice;
  I/D and compact/raw CIGAR span accounting agree.
- Direct-window negatives now use forward guide / oriented target, matching the
  Columba path; plus success no longer suppresses scanning the other strand.
- CFD remains an empirical matrix score with the existing DNA-bulge projection,
  not a new learned bulge model. Its limitations do not disappear with this fix.

Known release limitations: overlap selection still collapses candidates; the
window scanner is not an exhaustive search oracle. Negative `cg` is stored in
**guide orientation**, so reverse-complement the target before replaying CIGAR.
Some older benchmark/PAM oracles assume the erroneous two-base adjacency or
allow the historical query-coordinate overrun; do not reuse them unchanged.
Existing index presence checks in CRISPRapido accept partial indexes, so the
external preflight below requires every component. No claim of genome-wide
sensitivity/completeness is supported by the earlier imported-interval replay.

Do not call this a production release until full-search results and their
independent checks are reviewed. Do not tag/merge/publish automatically after
successful execution; produce a reviewable report first.

## 3. Immutable inputs and versioned outputs

Use the prior seed-20261020 guide/assembly selection to permit comparison with the
prior PAM audit. `guides_100.tsv` contains exactly 100 unique IDs and unique
20-base A/C/G/T sequences; `guides_100.fa` preserves their order and sequences.
This avoids selecting guides based on corrected output. The source is the prior
`targeted_PAM_audit/sampled_guides.tsv`, not a new random draw.

`assemblies_12.tsv` freezes ten haplotypes (four mat, four pat, two other):
HG00128_hap2, HG00735_pat, HG01081_pat, HG01496_pat, HG02055_mat,
HG02630_mat, HG03369_hap2, HG04184_mat, NA18906_mat, NA21309_pat.
Their filenames specify HPRC r2 v1.0.1. The additional references are exactly:

- CHM13: `chm13v2.0_maskedY_rCRS.fa.PanSN.fa.gz` (masked Y / rCRS derivative,
  not an unspecified CHM13 build).
- GRCh38: `GCA_000001405.15_GRCh38_no_alt_analysis_set.PanSN.fa.gz`
  (no-alt analysis set, not a different patch/alternate-contig build).

The manifest records exact paths, versions, source file sizes/mtime, contig and
base totals, FAI/GZI hashes, index prefix and component sizes. All 12 FASTA/FAI/GZI
sets and all 12 required Columba components were checked for existence/nonzero
size during preparation. Total assembly bases: 36,550,242,121. Largest existing
index set: 11,787,263,350 bytes.

Full assembly and index-component hashes are deliberately marked
`DEFERRED_TO_RUN_PREFLIGHT`: this preparation did not read/have to hash hundreds
of gigabytes of indexes. Do not replace these with guessed hashes or claim that
FAI/GZI hashes authenticate every base. At run time hash each exact source FASTA
and every index component, then verify their recorded construction provenance
binds them to the same contig dictionary and sequence content. A three-line
Columba `.meta` alone is not that provenance. If provenance is missing, resolve
it or rebuild a new index in a new location before claiming search completeness.
Do not replace the existing indexes. This is a preflight gate, not a completed
validation claim.

Reserved, newly created run root:

```text
/moosefs/salehi/HPRCv2_CRISPRapido/results/PAM_CFD_complete_search_v1_9c2d3f5_2026-10-28
```

It initially contains `PLAN_ONLY_NOT_RUN`, `plan/` and `release/`; no assembly
search output. Put all future outputs beneath this root, never under the old
Brunello/SBDS output directories. Each assembly attempt must use `mkdir` without
`-p`, `set -C`, and a new attempt name. Never use `--force`, overwrite, append
partial results to a finished run, or remove prior artifacts to retry.

## 4. Release artifact provenance

`release_manifest.json` records the exact optimized binary path/hash, immutable
code commit, compiler, Cargo lock hash, WFA2 library hash, Columba executable/hash,
score matrices, library loader paths and frozen input hashes. A documentation
commit may follow the code commit without changing the source binary; record
both rather than falsely claiming the executable was built from a later commit.
The manifest is generated after a successful isolated `--release --locked
--offline -j 2` build with Rust 1.85.1. Build log and environment are copied into
`release/` in the run root. The production executable is not replaced.

This host needs the explicit Guix glibc/GCC runtime recorded in the manifest.
Use the recorded runtime loader to invoke the candidate, and an executable
Columba wrapper in the run root with the same loader/library path. CRISPRapido
launches Columba as a child, so configuring only the outer executable is
insufficient. Never change global LD_LIBRARY_PATH or shell startup files.

The optimized build succeeded in 1m 34s. All **five optimized CLI regression
tests passed**, including both-strand PAM scoring, ambiguity/edge handling,
indel coordinates and direct-window strand coverage. Formatting and whitespace
checks passed, and source/dependencies/tests remain unchanged from the candidate
commit. The existing unused `Hit.cfd_score` warning remains; it is not a new build
failure. The packaged binary and Columba wrapper both passed `--help` smoke checks.

Candidate executable SHA256:
`f59b1bbf58e70b321359c7bd7b8687aad81bb874b22dc2f083319453726bc876`.
Columba executable SHA256:
`7e44974c81dad4943292fa64a71e54073908b4909ce327480464d769322c40d1`.
These are isolated release artifacts; the previous debug/replay hashes are
preserved in their original reports.

## 5. Proposed commands — do not execute searches yet

First, review provenance and confirm the run root is the prepared one. The
following shell variables refer only to frozen inputs/artifacts in that root:

```bash
set -euo pipefail
set -C
RUN=/moosefs/salehi/HPRCv2_CRISPRapido/results/PAM_CFD_complete_search_v1_9c2d3f5_2026-10-28
PLAN="$RUN/plan"
BIN="$RUN/release/crisprapido-9c2d3f5"
COLUMBA="$RUN/release/columba-runtime"
SAMTOOLS=/moosefs/guarracino/condatools/samtools/1.22.1/bin/samtools
ZSTD=/home/salehi/miniconda3/bin/zstd
LOADER=/gnu/store/ln6hxqjvz6m9gdd9s97pivlqck7hzs99-glibc-2.35/lib/ld-linux-x86-64.so.2
LIBS=/gnu/store/ln6hxqjvz6m9gdd9s97pivlqck7hzs99-glibc-2.35/lib:/gnu/store/gmv6n5vy5qcsn71pkapg2hnknyn1p7g3-gcc-12.3.0-lib/lib:/gnu/store/ygm2rzg08w0w96h9587k3n2yvhndcqkn-zlib-1.2.13/lib
sha256sum -c "$RUN/release/artifacts.sha256"
```

Run-time preflight, separately logged for all 12 manifest rows:

1. Verify guide hashes/count/unique IDs/alphabet and all assembly metadata/hashes.
2. Compare source contig IDs/lengths with the index construction manifest. Verify
   all suffixes `bwt brt rev.brt sa.4 sa.bv.4 txt.bin cct pos meta fsid
   headerSN.bin sna`. Hash them with absolute paths, preserving stdout in a new
   `preflight/<assembly>.sha256` and source dictionary in a new TSV.
3. Confirm planned binary SHA, matrices, WFA/Columba build identity and loader
   resolution. Record `uname -a`, CPU model, hostname, Slurm allocation (if any),
   `/usr/bin/time --version`, compiler, samtools/zstd versions and free disk.
4. Confirm local scratch capacity, cap threading, and establish output/log size
   limits. Fail on mismatch or unavailable input; do not silently substitute
   another CHM13/GRCh38 build, guide set or binary.
5. Complete index provenance gate before allowing automatic candidate generation.

For a single assembly, read its exact row without globbing (example CHM13):

```bash
ASM=CHM13
mapfile -t INPUT < <(python3 - "$PLAN/assemblies_12.tsv" "$ASM" <<'PY'
import csv,sys
rows=[r for r in csv.DictReader(open(sys.argv[1]),delimiter='\t') if r['assembly_id']==sys.argv[2]]
assert len(rows)==1
print(rows[0]['fasta']); print(rows[0]['index_prefix'])
PY
)
FASTA=${INPUT[0]}
INDEX=${INPUT[1]}
TASK="$RUN/$ASM.attempt01"
mkdir "$TASK"                  # fail if an earlier attempt exists
mkdir "$TASK/tmp"
export TMPDIR="$TASK/tmp"      # preserved Columba SAM/log subdirectories
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
printf '%s\n' "$FASTA" "$INDEX" > "$TASK/input_paths.txt"
# Preflight hashes must already have passed; do not run until approved.
/usr/bin/time -v -o "$TASK/search.time" \
  "$LOADER" --library-path "$LIBS" "$BIN" \
  -r "$FASTA" -G "$PLAN/guides_100.fa" \
  -p GG -m 4 -b 1 -z 1 -f 0.75 -t 4 \
  --mismatch-scores "$RUN/release/mismatch_scores.txt" \
  --pam-scores "$RUN/release/pam_scores.txt" \
  --columba-bin "$COLUMBA" --columba-index "$INDEX" \
  --keep-columba-sam \
  > "$TASK/all_strict.paf.partial" 2> "$TASK/search.stderr"
# Only after a zero exit and basic structural checks:
mv -n "$TASK/all_strict.paf.partial" "$TASK/all_strict.paf"
sha256sum "$TASK/all_strict.paf" > "$TASK/output.sha256"
```

This is **fresh whole-assembly search**: no `--columba-sam`, retained-candidate
replay, subsetting of contigs or CFD cutoff. Columba is called with `-a all -m
edit -e 5`, plus a separate exact `-e 0` recovery call. WFA verification and
strict NGG selection then run. Preserve both raw SAMs, query FASTA and Columba
logs; the verifier sees their union. Search all contigs including the named
reference derivatives' nonautosomal content. Do not drop zero-CFD or missing-CFD
rows during validation. Record each count category explicitly.

After CHM13 passes, repeat this exact command by manifest assembly ID for the
remaining eleven rows, retaining 100 guides in every job. Proposed scheduling
only: six allocated CPUs, 64 GiB RAM, two-hour wall limit, maximum two assembly
jobs at once; start with one job to measure resources. No `sbatch` command or
array has been submitted in this preparation. If the scheduler needs a wrapper,
it must contain this frozen command, use new per-assembly attempts, and save its
own script hash. Do not invoke the older full-Brunello submission scripts.

## 6. Independent checks and acceptance criteria

Process **every emitted row**, not just cf>=0.20 or a random sample. Use separate
bounded batches of 500–1,000 rows and samtools faidx against the exact indexed
source assembly; do not use the Rust-extracted target/PAM as the oracle. Keep a
record-level verification TSV and machine-readable summary for each assembly.
The independent Python building blocks already exist in
`scripts/validate_pam_subset.py`: `fetch`, `rc`, and `score_row`. Reuse those pure
functions from the frozen code snapshot for a new streaming all-output audit;
**do not run its subset main() on a full production result**. The full-output
driver is a required pre-execution deliverable and must be reviewed/tested before
the planned search is approved. Its command contract should be:

```text
python3 audit_all_emitted.py \
  --paf TASK/all_strict.paf --fasta EXACT_MANIFEST_FASTA \
  --guides PLAN/guides_100.tsv --samtools PINNED_SAMTOOLS \
  --matrices RUN/release --output TASK/independent_audit
```

`audit_all_emitted.py` is a **specified future validation driver, not an existing
script or a completed check**. Do not attempt that command until implemented and
reviewed. The existing small subset validator cannot stand in for it. The driver
must enforce these assertions and preserve diagnostic failures:

| Check | Independent expectation |
|---|---|
| Guide coverage | Exactly 100 input IDs; every output ID belongs to the frozen list; zero-hit guides remain in the summary |
| Contig/coordinates | Exact contig ID, target length equals FAI length; `0 <= s < e <= contig_len`; guide coordinates `0..20` |
| Full target fetch | samtools `contig:s+1-e`; negative rows reverse-complemented once before replay |
| Full PAM fetch | plus `contig:e+1-e+3`; minus `contig:s-2-s`, reverse-complemented; exactly three A/C/G/T bases, suffix GG |
| CIGAR structure | Full parse, positive lengths, `=/X/I/D` only; sum(query operations)=20; sum(reference operations)=e-s; alignment/block/match counts agree |
| Reconstruction | Every `=` has identical guide/target bases; every `X` differs; I has guide base/opposite gap; D has genomic base/opposite gap |
| Thresholds | <=4 mismatches, <=1 gap run/event, gap size <=1, >=75% guide matches; derive from sequence/CIGAR, not only tags |
| CFD | Reconstruct full gapped strings, apply existing DNA-bulge projection, independent mismatch/PAM matrices; absolute error <=0.000051 for four-decimal output |
| Missing CFD | Expected unsupported/ambiguous sequence gets explicit unscored status; missing cf when independently scoreable is a failure; duplicate/invalid/nonfinite/out-of-range cf fails |
| Strand coverage | Count both; inspect all gap rows and 20 fixed-seed rows/strand when available; an absent stratum is recorded, never claimed as tested |
| Boundaries | No negative-index wrap or truncated PAM; synthetic tests supplement real data where no natural edge case exists |
| Summary accounting | emitted=scored+unscored; scored=above+below threshold; canonical=emitted; no hidden malformed/failed rows |

The score oracle must preserve the established DNA-bulge projection, never
truncate an arbitrary alignment to 20 nt. Save Python version and validator
source hash. Streaming must verify every row with bounded memory, close/check
subprocess exit status, and fail without publishing a DONE marker on any error.

Also retain candidate generation/rejection accounting from raw SAMs: candidate
row counts per `e=5`/`e=0`, exact duplicate counts, union query IDs, WFA-rejected
warnings and final emitted counts. Raw candidate counts differ from post-overlap
loci. Do not label the difference solely as PAM rejection: the current executable
does not expose separate complete stage counters.

### Search coverage, distinct from emitted-row correctness

A complete genome search is not a mathematical completeness proof. The release
checks must additionally include a separately implemented/reviewed genome-wide
multi-pattern exact-match oracle for all 100 guides on both strands (20-nt patterns
plus three-base PAM). Compare all exact NGG loci with raw Columba exact SAMs and
final output, allowing explicitly recorded overlap representative equivalence.
The direct-window scanner is not an acceptable oracle. If a shorter/bulged
representative replaces an exact locus, record the actual pair and boundary rule,
not just a count-based tolerance.

For mismatches and bulges, independent verification of emitted rows establishes
soundness, not sensitivity. On a fixed reproducible set of contigs/windows from
every assembly, use exhaustive alignment/PAM enumeration at the same bounds and
compare raw candidates and final representatives. Include synthetic inserted
positive/negative perfect, mismatch, DNA/RNA-bulge and contig-boundary fixtures
with independently enumerated truth. Synthetic fixtures must use separate new
indexes/references and cannot be mixed into natural-assembly biological counts.
These independent coverage oracles are additional pre-run deliverables, not
completed implementations. Global bulged-search completeness remains unproven
unless an independent exhaustive method is actually applied at that scope.

Report old/corrected comparisons on the same frozen panel where raw data exist.
Do not run another complete old-binary search automatically. The old audit is a
retained-only sample: it cannot provide lost/below-threshold candidate denominators.
Compare shared loci, shifted representatives, gained/lost loci, PAM distributions,
CFD changes and threshold categories; preserve both the strict and historical
eligibility definitions in the report.

## 7. Resource envelope and stop conditions

These are conservative planning requests, not measurements of this release.
The local HG00099 1,000-guide benchmark recorded 21:32 wall and ~10.9 GiB peak
RSS, but its binary, guide mix and execution mode differ. The older CHM13
20-guide benchmark recorded ~14 GiB peak RSS with repeated index loading;
linear scaling from it would be misleading for this batch workflow.

| Stage | CPU allocation | RAM request | Wall/storage planning |
|---|---:|---:|---|
| Release build (preparation only) | 2 workers | ~4 GiB | minutes; separate Cargo target directory |
| Provenance/hash preflight | 1–2 | <=2 GiB | I/O-bound; source+index hashes may read >150 GB across 12 assemblies |
| Search, per assembly | 6 allocated / 4 search threads | 64 GiB | target 15–60 min, limit 2 h; <=2 concurrent jobs after pilot |
| Independent emitted-row audit | 2 | 4–8 GiB | budget up to 2 h/assembly, BGZF/random-I/O dependent |
| Exact/exhaustive coverage oracles | bounded separately | 8–32 GiB | estimate after implementation/tiny-fixture benchmark; not included in search timing |

Reserve >=50 GiB local scratch per concurrent search when staging an index,
FASTA and raw logs/SAM; largest existing index is ~11.8 GB. Initially reserve
50 GiB persistent free space for the 12 outputs and audits; measure actual
SAM/stderr growth in the pilot before scheduling the rest. Keep compression to
one worker outside the search; preserve uncompressed raw files until checks pass.
At the nominal 15–60 minutes/search, twelve searches need roughly 3–12 job-hours
(18–72 allocated CPU-hours at six CPUs), or ~1.5–6 hours at concurrency two,
excluding audits/provenance. This is a budget estimate, not a promised duration.

Stop and report on provenance mismatch, missing index component, invalid emitted
row, exact-coverage failure, non-NGG strict output, unexpected unsupported CFD,
OOM/nonzero command exit, timeout or rapidly growing logs/disk use. Do not
silently relax `-m/-b/-z/-f`, discard rows, switch to permissive PAMs, overwrite an
attempt or increase resources across all assemblies. Review pilot measurements
and choose a new versioned attempt if changes are needed.

## 8. Completion and promotion checklist

Before executing the search, finish/review the full-output auditor and coverage
oracles, authenticate index provenance, inspect the release manifest, and obtain
approval for the concrete 12-assembly run plan. This preparation itself is not
that approval.

After execution: all commands must exit zero; all 12 assemblies must have audited
summaries, 100 guide entries including zero hits, sequence/PAM/coordinate/CFD
checks and coverage results, binary/input/validator hashes and retained raw logs.
Publish a new report distinguishing tested categories from absent strata and
complete-search execution from unproven global sensitivity. Keep the original
results and previous reports immutable. Do not promote or merge automatically.
