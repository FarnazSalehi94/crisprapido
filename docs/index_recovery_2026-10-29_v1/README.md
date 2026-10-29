# Index recovery v1 — finalized 2026-10-29

[Recovery report](RECOVERY_REPORT.md) documents all twelve index sets, source
binding, failure offsets, recovery methods, tests and remaining validation scope.

- [Verified assembly/index paths](verified_assemblies.tsv): use these resolved
  paths rather than the nine damaged original prefixes.
- [Per-assembly status](index_status.tsv).
- [Original file checks](original_file_integrity.tsv) and [failure offsets](failure_offsets.tsv).
- [Verified component checksums](verified_component_checksums.tsv).
- [Recovery/auditor metadata](recovery_metadata.json) and [script hashes](script_hashes.json).

Six affected sets were rebuilt, three had only `.txt.bin` regenerated, and three
were retained intact. All twelve resolved sets passed complete verification.
The independent auditor passed 23 tests and all 88 prior real examples, including
checks at their original genomic coordinates. The complete search was not run.

Large artifacts are outside Git:
`/export/local/home/salehi/crisprapido-index-recovery-v1-2026-10-28`.
