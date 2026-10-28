# PAM/CFD release candidate preparation (2026-10-28)

Status: optimized candidate built; five optimized CLI tests passed; full-search
validation not launched; no merge/tag/promotion. A merge into main is stopped
because of divergent backends and conflicts.

- [Validation and release plan](VALIDATION_PLAN.md): safety review, exact commands,
  resource envelope, acceptance checks and pre-run deliverables.
- [Release manifest](release_manifest.json): source commit, binary/toolchain hashes,
  runtime paths and fixed parameters.
- [100-guide FASTA](guides_100.fa) and [guide table](guides_100.tsv).
- [12-assembly manifest](assemblies_12.tsv): ten HPRC haplotypes plus the exact
  CHM13/GRCh38 derivatives. Full large-input hashes remain an explicit preflight gate.
- [Merge review](merge_review.txt) and [Git state](git_review.txt): read-only evidence.
- [Optimized test log](optimized_cli_tests.txt).

The separate versioned run root is:
`/moosefs/salehi/HPRCv2_CRISPRapido/results/PAM_CFD_complete_search_v1_9c2d3f5_2026-10-28`.
It contains planning/release files and a `PLAN_ONLY_NOT_RUN` marker, not search
results. The original repository worktree and all prior analysis results remain
unchanged. CRISPR-Bulge integration has not started.
