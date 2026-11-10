# Optional CLI scoring: CFD and CRISPR-Bulge

## Implemented contract

Development baseline **`90d958a`**, branch `feature/pam-anchored-search-poc`. The Rust CLI now accepts **`--scoring cfd|crispr-bulge|both`**, default **`cfd`**. CRISPR-Bulge remains experimental and coverage-limited, not a calibrated cleavage probability or a replacement for CFD.

- **`cfd`**: original search/reporting/CFD path, original PAF bytes, no Python/model dependency, no added tags or sidecar by default.
- **`crispr-bulge`**: request the pinned neural scorer. Retain every original record **and any existing CFD tags**, append neural score/status tags and write full sequence/score annotations.
- **`both`**: explicitly request the combined reporting described above. The current compatibility-first implementation produces the same score columns/tags as `crispr-bulge`, with the selected mode recorded in metadata. It intentionally does **not** suppress existing CFD calculation or strip existing CFD values in neural-only mode. This preserves the user's requirement to retain both scores whenever available; the two scores are not fused.

For opted-in modes, reporting writes the ordinary PAF to a disk spool, then invokes the existing batched scorer **once for the entire invocation**, including multi-guide runs. All five pretrained `c_2` classification models load once and are reused across chunks/guides. No retraining, new model, alternate representative, alignment truncation, projection or domain relaxation was introduced. Inference is CPU-only in a user-supplied isolated interpreter with the existing TF-CPU 2.12.0 environment and checksum-pinned HDF5 files. No runtime downloads or inference-time pickle loading.

**Neither search algorithm, representative-selection policy, mismatch/bulge defaults nor CFD formulas changed.** Reporting now accepts a destination writer; formatting remains unchanged. The sole verification-module edit initializes the new CLI field in a test fixture. Integrated pruning remains disabled. Only isolated development executables were built; production binaries/indexes/original results were not modified. No new whole-genome search, Slurm job or main merge occurred. Tiny fixture searches were used for regression testing.

## Output and unavailable scores

Default CFD stdout remains exactly the original PAF. In opted-in modes, stdout and `SCORING_DIR/scoring/output.paf` retain all original PAF fields/tags/order and append:

- `cb:f:<value>` only for a defined neural score (nine significant digits);
- `bu:Z:scored` or the first unavailable reason for every record.

**Never put the string `N/A` in a numeric PAF tag.** Missing numeric tags retain PAF compatibility; explicit **`N/A`** values and reasons appear in `SCORING_DIR/scoring/annotations.tsv`, one row per original hit, in the same order:

`row_index, guide_id, contig, start_0based, end_0based, strand, guide_sequence, target_sequence, pam, cigar, aligned_guide, aligned_target, mismatches, bulge_events, max_bulge_size, cfd, cfd_reason, crispr_bulge, crispr_bulge_reason`.

Target/aligned sequences and the complete observed three-base PAM are guide-oriented, including reverse complementation for minus-strand hits. Coordinates are the original forward-reference zero-based half-open intervals. Gap characters represent the original CIGAR, not a new realignment. Ambiguous sequence bases remain literal ambiguous bases; an incomplete PAM is `N/A`. True score zero remains numeric zero, never unavailable.

Missing CFD is `N/A` / `upstream_cfd_unavailable`: a saved PAF's omitted CFD tag cannot retrospectively distinguish a missing matrix from an unsupported formula/input, so no more specific cause is fabricated. Default CFD PAF retains its historical omitted-tag representation for unavailable CFD, with no new default annotation file. Rich context and explicit reasons are provided in the opted-in sidecar.

No row is removed because neural scoring is unavailable. **No score-based filtering is enabled by selecting any scoring mode.** Existing optional `--min-cfd` remains disabled by default and applies only when explicitly requested during a search, after the original representative is selected. It is not a neural threshold. `--score-paf` rejects `--min-cfd` rather than re-filtering rounded stored scores.

The old standalone Python postprocessor retains its historical CB-only behavior unless passed `--preserve-cfd`; **the Rust integration always passes that flag** plus `--annotations`. Existing bounded tools/tests remain reproducible without silently changing their old mode semantics.

## CLI use

A normal search retains its existing inputs/backend options; add the following only when opting in:

```bash
# Existing search command/options + optional neural annotation:
crisprapido -r reference.fa -G guides.fa \
  --columba-sam candidates.sam \
  --scoring both --scoring-output-dir NEW_SCORING_DIR \
  --crispr-bulge-python /absolute/optional-env/bin/python \
  --crispr-bulge-model-root /trusted/pinned/CRISPR-Bulge > annotated.paf
```

The same scoring options apply to the PAM-anchored backend without changing its required index verification or unpruned search. The new output directory must not exist and its parent must exist. The source bundle is needed for Python inference; `--scoring-script /path/to/matching/score_paf.py` overrides the build-time source path when relocating. Keep its sibling adapter/annotation/reference-index/worker/manifest files together. No model packaging or production installer was added. See [requirements](requirements-inference.txt), [tested environment](requirements-validation.lock), [model setup](REPORT.md#7-reproducibility-artifacts-and-stop-point) and [licensing notice](UPSTREAM_NOTICE.md).

### Existing results: never launch a new search

```bash
crisprapido -r reference.fa -G exact-original-guides.fa \
  --score-paf saved-original.paf --scoring both \
  --scoring-output-dir NEW_SCORING_DIR \
  --crispr-bulge-python /absolute/optional-env/bin/python \
  --crispr-bulge-model-root /trusted/pinned/CRISPR-Bulge \
  --scoring-reference-sha256 TRUSTED_ORIGINAL_REFERENCE_SHA256 > annotated.paf
```

`--score-paf` dispatches before search, reference loading by the search engine, or CFD initialization. It uses existing CFD values, never recalculates them. Reference/guide CLI inputs remain required; for `-g`, legacy query ID `Guide` is expected, while `-G` must match all saved PAF IDs (HG00735 uses `pilot`). Search executables/candidate/index flags and `--min-cfd` are rejected with this mode. Other search limits are not reapplied to an existing PAF. Default/explicit CFD `--score-paf` is byte-preserving pass-through and requires no models or reference scan. Never direct stdout to the input file: shell redirection happens before the program can protect it.

### Failure behavior

Configuration/path errors fail before search where checked. Opted-in live-search PAF is retained as `SCORING_DIR/input.paf`; scoring-only mode only reads the supplied original. `SCORING_DIR/status` starts `INCOMPLETE`. Worker/protocol/checksum/timeouts are failures, **not biological unavailable predictions**. No PAF is published to stdout unless scoring succeeds and the child bundle reports complete. Complete bundles contain metadata, models/inputs/index digests, annotations, worker logs/metrics and status. On failure, preserve the incomplete directory for diagnosis; existing output directories are never reused/overwritten. Consumers should trust complete status, not the mere presence of a data file.

`--scoring-batch-size` defaults to 512, range 1–4096; `--scoring-worker-timeout` defaults to 300 seconds and covers the neural worker, not reference validation/search. There is no new hard memory limit or global invocation timeout in the CLI. Large disk spools/anchor-key sets still scale with record count. Default CFD retains its original streaming behavior.

## Reducing redundant reference validation

The previous 91.6 s annotation spent 81.8 s scanning Python FASTA lines to validate alphabet/wrapping/build offsets. This is now **amortizable**, without skipping reference identity checks:

```bash
# One-time explicit preparation into a NEW external file; no reference writes.
python3 research/crispr_bulge/reference_index.py \
  --reference reference.fa --output /trusted/new-reference-index.json
# Prints reference_sha256 and index_sha256. Retain these trusted pins.
```

On a subsequent neural invocation, add:

```text
--scoring-reference-index /trusted/new-reference-index.json
--scoring-reference-index-sha256 TRUSTED_INDEX_SHA256
--scoring-reference-sha256 TRUSTED_REFERENCE_SHA256
```

Preparation uses the original full wrapping/alphabet validator and a post-validation full hash. Reuse requires **both pins**; checks index digest/schema/reference binding/size, hashes the **entire reference before scoring**, verifies fetched intervals, and hashes the **entire reference again afterward**. Index digest is also rechecked. Without a pinned index, the original complete scan remains the fallback. No automatic sidecar or index is written beside the reference; no native FM-index is rebuilt. Existing mandatory search/index binding checks are untouched.

This is an explicit trusted local prepared index, not an unsigned automatic cross-job trust cache. Obtain pins from trusted preparation/search provenance, not from an untrusted supplied index. Hashes cannot defend against malicious same-user mutate-and-restore races; trusted inputs/process isolation remain necessary. Currently supports plain regularly wrapped FASTA, not gzip/irregular wrapping in opted-in neural modes; default search input support is unchanged. The CLI does not yet directly reuse in-memory reference data from search, so full pre/post hash I/O remains.

## Saved HG00735 validation and performance

All three CLI modes were run on the saved **21,447-record job 2871283** output, not a new genome search. Original input SHA256 remains **`a28d0574a915a57d431851ab681b0602b8f5fb1740d1c8efeb6e3b4a378f23be`**; exact reference SHA256 matches saved search provenance.

| Mode | Hits retained | Existing CFD retained | Neural defined | Neural unavailable |
|---|---:|---:|---:|---:|
| Default CFD | 21,447 | 21,447 | Not requested | Not evaluated |
| `crispr-bulge` | 21,447 | 21,447 | 4,526 | 16,921 |
| `both` | 21,447 | 21,447 | 4,526 | 16,921 |

Both neural sidecars have exactly 21,447 rows. All coordinates/CIGARs/fields/CFD tags/order were compared with the original; sidecar IDs, coordinates, guide sequences, CIGARs and scores were checked row-for-row. CFD output is byte-identical to the original. Neural outputs in both modes are byte-identical to each other **and the previous evaluated `both` PAF** (SHA256 `35a4c0441c90434c05334d5cb97221d83b4364c73ca324d35adedad0f09d15f1`).

Unchanged unavailability: **8,176 DNA2** exceed 24 columns, **7,624 DNA1/five-mismatch** cases exceed conservative with-gap mismatch policy, **1,121 RNA2** exceed the one-gapped-nucleotide policy. These are retained with `N/A` and the same reasons, not coerced to zero. Coverage remains **21.10%**; full `m5/b1/z2` model support is not claimed.

### Full official prediction validation

Unlike the earlier 93-record sample, this time **all 4,526 supported alignments** were independently reconstructed with samtools/htslib extraction and separate CIGAR replay, then passed to **unmodified official `ensemble_predict`** at upstream commit `3eddcd5bfcaff00b2bdf29425116ec9756ace870` with the pinned five-component pretrained ensemble. Source/HDF5/trusted reference-only pickle hashes were checked. Categorical features matched exactly (SHA256 `19481ec499dddfbcc0df92b9aab97d3cc13e1da8ca040c21d649f63f92c6fdc7`). All predictions in **both CLI modes** agreed within `1e-6`: maximum absolute difference **2.5598145825256324e-10** from actual PAF scores. This validates conversion/inference implementation agreement, not biological accuracy, calibration or genome-wide recall. Existing [score-association caveats](ASSEMBLY_EVALUATION.md#cfd-versus-crispr-bulge) remain unchanged.

### Measured costs

Single fresh-process CPU call per mode on `octopus01`, warm/uncontrolled storage caches; includes CLI/driver startup, full pre/post reference hashes, conversion, model loading/inference, annotations and PAF publication. Prepared-index creation is **separate**, not hidden in the annotation timings. No repeated median or whole-workflow search speedup is claimed.

| Mode | External child wall s | Maximum individual descendant RSS KiB | Sampled summed tree peak KiB |
|---|---:|---:|---:|
| CFD saved-file pass-through | **0.007562** | 12,812 | 24,492 |
| `crispr-bulge` | **12.929960** | 573,596 | 629,832 |
| `both` | **13.781680** | 577,588 | 634,612 |

CFD row is file pass-through, **not CFD calculation time**. Neural worker time was **6.115 / 6.133 s**, with exactly five loads and nine chunks per mode. Initial pinned-index/reference hash validation was **2.227 / 2.226 s**, replacing the previous 81.784 s full Python validation pass; the final reference hash is still performed. Compared with the previous single 91.576 s observation, annotation was substantially faster, but scopes now also include richer sidecars and storage/shared-host effects were uncontrolled.

Measured one-time full index preparation: **75.4562 s**, max individual RSS **19,588 KiB**. A separately prepared index had the exact same digest: `c1b5465af849d2d6e7e3e20827266304e77d64f0d3795dff891005de83b1761e`. This cost must be included for a first-use end-to-end comparison; the benefit is amortizing repeated wrapping/alphabet validation, not making verification free.

Validation used the existing 20 ms requested process-tree sampler with a 4 GiB sampled stop and 900 s outer timeout; neither fired. Actual intervals include sampling overhead. Summed RSS may double-count shared pages/miss spikes; `ru_maxrss` is a maximum individual descendant high-water mark, not aggregate unique memory. The full official-model check took 6.467 s separately and reached 759,784 KiB RSS; it is not included in annotation timing. Whole-assembly live-search-plus-scoring peak memory was not measured; the Rust process can still retain search allocations while the scoring child runs.

## Tests, artifacts and reproduction

- **164 Rust tests passed**, one pre-existing external-Python CFD test ignored. Existing unused-field/import warnings remain.
- **20 Python tests passed**, none skipped in the recorded configured run: previous tests plus CLI default/explicit CFD, both neural modes, tiny live single/multi-guide reporting, prior-CLI byte compatibility, unsupported DNA2/RNA2/mixed bulges, `N/A` versus actual zero, missing CFD, empty inputs, model load reuse, worker failure/timeout, exclusive directories, and pinned-index corruption/reference-mismatch rejection.
- No source changes to either search backend, CFD formulas/matrices, Cargo files or production binaries. Verification change is test setup only.

Small committed evidence: [summary](results/cli-integration/summary.json), [Rust tests](results/cli-integration/rust-tests.log), [Python tests](results/cli-integration/python-tests.log), [preparation measurement](results/cli-integration/index-preparation-measurement.json). Models, environments, prepared index, reference, original/annotated PAF, annotations and sequence-bearing official requests remain outside Git.

External root: `/export/local/home/salehi/crisprapido-scoring-cli-v1/`. Development executable: `target/release/crisprapido`; exact source/binary digests retained in the summary. `saved-results/{crispr-bulge,both}/scoring/` retain complete bundles; `saved-results/cfd.paf` is original-byte pass-through; `saved-results/official-all/` retains independent full-reference validation. `reference-index.json` is the trusted prepared index; `reference-index-reproduced.json` came from the measured preparation. No original results were overwritten.

Reproduce saved-result validation with a **new external output directory** and a separately built development executable:

```bash
PYTHONDONTWRITEBYTECODE=1 /absolute/optional-env/bin/python research/crispr_bulge/validate_cli.py \
  --cli /absolute/development/crisprapido \
  --input /export/local/home/salehi/crisprapido-resource-pilot-results/hg00735-optimized-2871283/output.paf \
  --reference /export/local/home/salehi/crisprapido-hprc-import-hg00735-20261105-v1/reference.fa \
  --reference-index /trusted/prepared-reference-index.json \
  --upstream /trusted/pinned/CRISPR-Bulge --model-python /absolute/optional-env/bin/python \
  --out /NEW/external-cli-validation

CRISPR_SCORING_CLI=/absolute/development/crisprapido \
CRISPR_MODEL_PYTHON=/absolute/optional-env/bin/python \
CRISPR_MODEL_ROOT=/trusted/pinned/CRISPR-Bulge \
CRISPR_OLD_CLI=/absolute/preserved-preintegration/crisprapido \
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover \
  -s research/crispr_bulge -p 'test_*.py' -v
```

Use an external `CARGO_TARGET_DIR` for `cargo test --release --locked --offline -- --test-threads=1`, preserving existing binaries; this site's validated Guix compiler/linker setup matches the earlier [bounded development recipe](../pam_anchored_poc/cfd_pruning/reproduce.sh). Only build/test commands were used, not that recipe's later experiment steps. Full official validation requires samtools and the reference environment's NumPy/pandas/auxiliary packages; ordinary CLI inference needs only the inference environment. CLI tests deliberately skip unless their optional environment variables are provided.

**Remaining limitations:** limited-domain model coverage, uncalibrated scores/representative sensitivity, optional source-bundle/TF dependency management and broader licensing/packaging review, plain-wrapped-reference restriction, pre/post full-reference hashing, disk spooling and anchor memory for huge outputs, no hard CLI memory cap, no GPU/concurrency/wide-guide biological validation or genome-wide recall guarantee. No neural pruning, automatic thresholds, retraining, domain extrapolation or further whole-assembly search is authorized by this integration.
