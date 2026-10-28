# CRISPRapido complete-search validation: STOPPED at prerequisite gate

Date: 2026-10-28. Repository: FarnazSalehi94/crisprapido, branch
`fix/pam-cfd-validation`. Starting documentation commit:
`db67867e4cc45def6f93f8a30a421e30f3a3fa16`. Candidate executable source remains
`9c2d3f53673fd135f0b4379812a9a7a4a7f20961`; no application code or binary was changed.

**Decision: do not launch the frozen 100-guide / 12-assembly search.** The full
indexed-text verification for HG00128_hap2 failed because the index file is not
fully readable. An independent Python binary read confirms operating-system
`OSError(5, 'Input/output error')`. This is a storage/readability prerequisite
failure, not evidence of a differing nucleotide or a newly identified PAM/CFD bug.

## Inputs and provenance checks

All checks are read-only. A new `preflight_attempt01/` holds prototype/diagnostic
logs; the structured preflight is in `preflight_attempt02/`. Both are under the
previously reserved versioned root:

```text
/moosefs/salehi/HPRCv2_CRISPRapido/results/PAM_CFD_complete_search_v1_9c2d3f5_2026-10-28
```

No historical result, frozen plan, production executable or original dirty
checkout was overwritten. The plan's `PLAN_ONLY_NOT_RUN` marker is retained.

| Check | Result |
|---|---|
| Frozen plan hashes (100-guide FASTA/TSV and 12-assembly manifest) | PASS |
| 100 unique guide IDs and sequences, each 20 A/C/G/T bases; FASTA/TSV agree | PASS |
| Packaged release-artifact and external-tool hashes | PASS |
| All twelve FASTA/FAI/GZI and required index components open/nonempty | PASS at initial-byte check; **not** proof of full readability |
| FASTA sizes/mtime and FAI/GZI hashes match frozen manifest | PASS, all twelve |
| Contig IDs, order and lengths match index SAM header | PASS, all twelve |
| Index position/name arrays and single-FASTA file-ID array match FAI | PASS, all twelve |
| Index format | All twelve report build tag 21, 32-bit offsets, VANILLA |
| Original index build logs | Present for all twelve |
| Complete FASTA-to-indexed-text comparison | **FAIL: HG00128_hap2 index read error** |
| Full FASTA/index-component SHA256 set | Not completed; deliberately stops before hashing after failed comparison |
| Full text verification for other eleven assemblies | Not attempted after failure |

Names checked: HG00128_hap2, HG00735_pat, HG01081_pat, HG01496_pat,
HG02055_mat, HG02630_mat, HG03369_hap2, HG04184_mat, NA18906_mat,
NA21309_pat, CHM13 and GRCh38. Exact source versions/paths remain in
[the frozen assembly manifest](../release_candidate_2026-10-28/assemblies_12.tsv).
Neither CHM13 nor GRCh38 is certified by the HG00128_hap2 check.

### Failed indexed file and independent diagnosis

```text
/moosefs/salehi/HPRCv2_CRISPRapido/indices/HG00128_hap2/HG00128_hap2.txt.bin
```

Its reported size is 3,025,377,917 bytes. The little-endian length prefix declares
3,025,377,913 text bytes (3,025,377,912 sequence bases plus the terminal `$`).
The streaming comparison stopped after 2,012,217,344 sequence bytes; the next
1,048,576-byte indexed read could not finish. The diagnostic did not find a
mismatching base before the failed read.

Independent Python `open(..., 'rb'); seek(offset); read(n)` results:

| Absolute file offset | Requested bytes | Result |
|---:|---:|---|
| 2,012,217,348 | 1,024 | 1,024 bytes read |
| 2,012,217,348 | 1,048,576 | OSError 5: Input/output error |
| 2,011,168,772 | 1,048,576 | 1,048,576 bytes read |
| 2,013,265,924 | 1,048,576 | OSError 5: Input/output error |

A successful initial read and correct file size therefore did not establish
full readability. The underlying MooseFS/chunk-storage cause has not been
determined; no repair/reconstruction was attempted. Do not label this as an
assembly-version mismatch, silently substitute an index, or run a search with
this index. Recovery must preserve the current evidence and use a separately
versioned index if rebuilding is necessary.

## Provenance verifier prepared and tested

New scripts:

- `scripts/release_validation/verify_inputs.py`: freezes/validates inputs and
  index dictionaries, invokes full transformed-text comparison, hashes complete
  source/index files only after the comparison passes, and records a failed
  summary without a PASS marker if any prerequisite fails.
- `scripts/release_validation/compare_index_text.cpp`: independent bounded-memory
  reconstruction of Columba's indexed text from the compressed FASTA.
- `scripts/release_validation/test_index_provenance.py`: four small-fixture tests.

The text format and transformation were checked against local Columba source
commit `edcc7ce2a1e28ded6784a44afbd8b9348a4805ad`, `buildindex.cpp` and
`definitions.h`: VANILLA build tag 21, seed length zero, uppercasing, replacing
non-ACGT bases using `minstd_rand(42)` with uniform A/C/G/T selection, concatenating
contigs, then a terminal `$`. `.pos`, `.sna`, `.fsid` and `.headerSN.bin` are
verified separately so identical concatenated text alone cannot hide different
contig boundaries. This reconstructs source binding; it does not alone prove
every BWT/suffix-array component is algorithmically correct.

**Four verifier tests passed.** They cover multiple records, lowercase and known
ambiguous-base normalization, crossing the 1 MiB buffer boundary, a wrong base,
and truncated/wrong-length/sentinel failures. A local compiler/include-path
mismatch was corrected by clearing inherited include/library environment variables
for this helper only; no installed toolchain was modified.

### Exact preparation and failing preflight commands

Executed in the isolated fix worktree:

```bash
env -u CPATH -u CPLUS_INCLUDE_PATH -u C_INCLUDE_PATH -u LIBRARY_PATH \
  /usr/bin/g++ -O2 -std=c++14 \
  scripts/release_validation/compare_index_text.cpp -lz \
  -o /export/local/home/salehi/compare_index_text_pam

python3 -m unittest discover -s scripts/release_validation -p 'test_*.py' -v

python3 scripts/release_validation/verify_inputs.py \
  --plan docs/release_candidate_2026-10-28 \
  --compare-bin /export/local/home/salehi/compare_index_text_pam \
  --output /moosefs/salehi/HPRCv2_CRISPRapido/results/PAM_CFD_complete_search_v1_9c2d3f5_2026-10-28/preflight_attempt02
```

The first text-comparison prototype did not print the failing offset; a second
build added it and an independent Python read confirmed the I/O error. The
committed helper adds errno context as well. Output-directory creation is
exclusive: any future retry must use `preflight_attempt03` or another new path,
not overwrite attempt02. No assembly search command was executed.

## Auditor readiness and deferred validation

The user's instruction to stop on a failed prerequisite takes precedence over
continuing later phases. Independent **full-search** PAM/CFD/CIGAR/coordinate/
strand auditors and coverage oracles were therefore not completed or certified
in this attempt. The prior 148 application tests, five optimized CLI tests and
208-row independent replay remain valid evidence at their documented scope;
they do not constitute the missing full-output or full-search audit.

Positive/negative emitted-row auditor tests were not run in this attempt. The
four tests above validate the new **index-provenance** checker, not CFD or strand
correctness. Do not conflate those test categories.

After index recovery, all twelve source/index bindings and full hashes must pass,
then the independently implemented streaming emitted-row auditor must be tested
on both strands and deliberately corrupted PAM/CIGAR/coordinate/score fixtures.
The exact-match and bounded exhaustive mismatch/bulge coverage oracles specified
in the existing plan are also still pre-run deliverables. Do not skip them to
launch the search.

## Parameters, resources and planned search command (NOT EXECUTED)

Frozen parameters remain `-p GG -m 4 -b 1 -z 1 -f 0.75 -t 4`, strict canonical NGG,
100 guides, ten haplotypes plus the exact frozen CHM13/GRCh38 derivatives.
Candidate code commit and binary SHA256 remain:

```text
9c2d3f53673fd135f0b4379812a9a7a4a7f20961
f59b1bbf58e70b321359c7bd7b8687aad81bb874b22dc2f083319453726bc876
```

Search resource requests remain six allocated CPUs, four search threads, 64 GiB
RAM, two-hour wall limit per assembly, first one pilot and at most two concurrent
after pilot review. Planning estimate: 15–60 minutes/search, 3–12 job-hours across
twelve, excluding provenance and independent auditing. These are unmeasured
estimates, not completed-run timings. The streaming provenance comparator uses
roughly a few MiB of memory and one CPU; full component hashing is I/O bound.
The failed large read can stall while storage retries, so wall time is not a
useful search-performance measurement. No full hashes/timings were fabricated.

The exact future loader/wrapper/runtime command is in the unchanged
[validation plan, section 5](../release_candidate_2026-10-28/VALIDATION_PLAN.md).
Its core invocation, valid only after all prerequisites pass, is:

```bash
"$LOADER" --library-path "$LIBS" "$BIN" \
  -r "$FASTA" -G "$PLAN/guides_100.fa" \
  -p GG -m 4 -b 1 -z 1 -f 0.75 -t 4 \
  --mismatch-scores "$RUN/release/mismatch_scores.txt" \
  --pam-scores "$RUN/release/pam_scores.txt" \
  --columba-bin "$RUN/release/columba-runtime" --columba-index "$INDEX" \
  --keep-columba-sam
```

## Outcome and next required action

No whole-assembly results, cross-assembly candidate counts, new CFD/CIGAR
validation totals or successful full-validation conclusion exist for this
attempt. Step 3 was not executed. Preserve/recover the unreadable index storage
or construct a new separately versioned index, then rerun prerequisites with
fresh output paths. Do not infer that all twelve indexes are good from their
metadata checks.

**Release and CRISPR-Bulge integration readiness cannot be certified.** The
corrected application has earlier unit/small-replay validation, but this release
validation is blocked by index I/O and the remaining auditor/full-search gates.
No merge into main, production binary modification, historical-result overwrite,
large-data push or CRISPR-Bulge implementation was performed.
