# Release validation tools

These tools do not launch the frozen complete search. Large outputs belong in
new directories outside Git. See
[the finalized recovery report](../../docs/index_recovery_2026-10-29_v1/RECOVERY_REPORT.md)
and its resolved manifest for the current index paths.

- `scan_index_integrity.py`: sequential, bounded full-file reads and SHA256s;
  continue to all assemblies after a failed file. Failure offsets are exact read
  request offsets/ranges, not claims about every corrupt byte. Incomplete reads
  never receive a full checksum. Resuming preserves prior per-file evidence.
- `verify_inputs.py`: frozen-input/index dictionary checks and a stop-on-first-fail
  preflight. `compare_index_text.cpp` performs the complete source/text comparison.
- `rebuild_index.py`: one affected set, verified source staging, new output only,
  pinned builder/runtime and resource logs. Refuses healthy-set rebuilds and
  unverified sources. Use the host's actual GNU time path.
- `rebuild_queue.py`: maximum two builders, no queries, failed-child propagation.
- `recover_text_component.py` / `regenerate_text.cpp`: when only `.txt.bin` failed,
  copy/checksum all other components and reconstruct text without rebuilding BWT/SA.
- `verify_recovery_set.py`: exactly twelve distinct assemblies, complete recovered
  component reads, source/text and dictionary checks, readable-original checksum
  comparisons; nonzero exit and no PASS marker on blocked verification.
- `audit_emitted.py`: independent streaming all-emitted-row auditor using indexed
  FASTA. This implements the prior plan's proposed full-output driver. It checks
  strict NGG, orientation, CIGAR, coordinates, thresholds/tags and CFD reconstruction.
  It is scoped to this release's guide-oriented negative-strand CIGAR convention.
  It does not establish search recall, and a zero-row soundness pass is not a
  completeness claim. Keep raw output until downstream checks finish.

Run small synthetic/real-fixture and failure-control tests:

```sh
python3 -m unittest discover -s scripts/release_validation -p 'test_*.py' -v
```

C++ helpers require g++ and zlib headers. The recorded host needed inherited
CPATH/CPLUS_INCLUDE_PATH/C_INCLUDE_PATH/LIBRARY_PATH cleared for its system compiler.
No Python packages beyond the standard library are needed. The end-to-end samtools
fixture uses the explicit tool path recorded in the release manifest; it reports
an explicit skip if that tool is absent elsewhere.

Full-content integrity and emitted-row soundness are distinct from candidate
sensitivity/completeness. Neither an index checksum nor a row auditor establishes
that every off-target was found. The full validation remains a separate task.
