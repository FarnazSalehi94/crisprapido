# CRISPR-Bulge: limited optional scoring integration

**Subsequent saved whole-assembly evaluation:** see [HG00735_pat job 2871283](ASSEMBLY_EVALUATION.md). All 21,447 existing records were annotated without a new search: 4,526 CB-scored / 16,921 unavailable; 93 independently sampled predictions agreed with official inference. Measured annotation wall time was 91.576 s, dominated by reference scanning. This does not promote the scorer into the core CLI. The bounded-development results below are retained as their original milestone.

## Decision and current status

**Implemented an isolated, batched PAF postprocessor for a conservative supported subset. It is not complete coverage of CRISPRapido's defaults `-m 5 -b 1 -z 2`.** Unsupported alignments remain in the output with an explicit unavailable reason, never a fabricated zero score. No truncation, gap splitting, realignment, model resizing or retraining workaround was implemented. The limitations were reported before implementation.

**The Rust CLI and both search backends are unchanged.** CFD remains the default, `--min-cfd` remains optional/disabled by default, deterministic representative selection is unchanged, and integrated CFD pruning remains disabled. New scoring choices are available in `research/crispr_bulge/score_paf.py`, **not as a new Rust `--scoring` option yet**. No production executable, existing index, frozen branch or main was modified. No whole-genome search/scoring experiment or Slurm job was run during this initial bounded milestone; the subsequent saved-result scoring evaluation is linked above.

## 1. Official source, models, license and dependencies

Primary source: [OrensteinLab/CRISPR-Bulge](https://github.com/OrensteinLab/CRISPR-Bulge), pinned default-branch commit **`3eddcd5bfcaff00b2bdf29425116ec9756ace870`** (2024-07-04). Web search was unavailable due to quota; the review used the actual GitHub Git source, checkpoint bytes and primary documentation rather than search snippets. Relevant pinned sources:

- [README and requirements](https://github.com/OrensteinLab/CRISPR-Bulge/blob/3eddcd5bfcaff00b2bdf29425116ec9756ace870/README.md)
- [main_predict.py](https://github.com/OrensteinLab/CRISPR-Bulge/blob/3eddcd5bfcaff00b2bdf29425116ec9756ace870/main_predict.py), which recommends the five-component `c_2` classification ensemble
- [ensemble_predict](https://github.com/OrensteinLab/CRISPR-Bulge/blob/3eddcd5bfcaff00b2bdf29425116ec9756ace870/train_and_predict_scripts/utilities.py)
- [feature construction](https://github.com/OrensteinLab/CRISPR-Bulge/blob/3eddcd5bfcaff00b2bdf29425116ec9756ace870/OT_deep_score_src/data_processing_utilities.py)
- [published GUIDE-seq preprocessing](https://github.com/OrensteinLab/CRISPR-Bulge/blob/3eddcd5bfcaff00b2bdf29425116ec9756ace870/process_guide_seq/process_identified_files/process_with_SWOffinder.py)

Selected model: **GRU-Embd (`c_2`), classification, arithmetic mean of ensemble components 0–4**, transferred from CHANGE-seq to FullGUIDEseq, with RHAMPseq excluded according to the saved path. This is the upstream example, not an arbitrary substitute. Component-zero HDF5 files for the saved `c_3`/regression alternatives were also inspected: `c_3` takes `(None,24,25)` one-hot inputs, `c_2` takes `(None,24)` categorical inputs, and regression outputs are linear rather than sigmoid. All inspected variants still have 24 alignment positions; changing between these inspected variants does not solve the 25-column DNA-bulge problem. Only the `c_2` classification ensemble was prediction-validated/integrated; the other variants were not benchmarked. [Variant inspection](results/variant-inspection.json). Classification output is not a calibrated biological cleavage probability, read count, or a CFD-equivalent quantity.

All five selected HDF5 files were inspected: input shape **`(None,24)`**, categorical vocabulary **25**, embedding dimension **44**, output **one sigmoid unit**, saved Keras version **2.12.0**. Each file is **2,782,400 bytes**, totaling **13,912,000 bytes**. Checkpoint pickles were disassembled before reference execution: they specify 24 positions, 25 pair categories and one-hot-to-categorical conversion. The weights are regular Git blobs at this commit; `datasets.zip` is Git LFS and was **not downloaded or needed**. [Pinned paths/hashes](model_manifest.json), [inspection evidence](results/source-examination.json).

Upstream has an [MIT root license](UPSTREAM_LICENSE), copyright 2024 Orenstein Lab. No separate checkpoint-specific terms were found in the reviewed paths. This is not a legal certification of all underlying dataset/dependency rights. **No weights, pickles, upstream Python code or datasets are vendored into CRISPRapido.** Attribution and redistribution caveats are in [UPSTREAM_NOTICE.md](UPSTREAM_NOTICE.md); full dependency/data licensing still needs review before production packaging.

Upstream tested Python **3.10.10**, NumPy **1.23.5**, pandas **2.0.2**, scikit-learn **1.2.2**, SciPy **1.10.1**, TensorFlow **2.12.0**, XGBoost **2.0.1**, CatBoost **1.2.2**. This study used a new isolated **Python 3.11.15** venv, **tensorflow-cpu 2.12.0**, the same listed numeric/reference package versions, Keras **2.12.0**, and h5py **3.16.0**. `pip check` passed. This is an explicitly validated CPU runtime variation, not a claim of bitwise equivalence to the authors' Python 3.10/GPU environment. The worker needs only TF/NumPy/HDF5 dependencies; the unmodified reference imports the extra training/auxiliary packages even though they are not used for this ensemble. [Inference requirements](requirements-inference.txt), [reference requirements](requirements-reference.txt), [complete tested environment](requirements-validation.lock).

## 2. Compatibility with mismatches and bulges

The checkpoint has **24 alignment columns INCLUDING the three-base PAM**. A 20-base guide with D extra target bases has `23+D` aligned columns. This is not CFD's 20-column projection convention.

The upstream preprocessing call is `6 6 4 1 1 false 50 NGG true`. Its pinned SWOffinder submodule commit is **`25bbffb85a5a74bb54a97a469335b5ff17202497`**, not the newer current SWOffinder default branch. That version parses these as maximum edits 6, mismatches without gaps 6, mismatches with gaps 4, **one gapped nucleotide**, and one thread. Its recursion increments the bulge count per gapped base, unlike CRISPRapido's contiguous-event count. This documents a conservative preprocessing subset; we did not inventory all training rows or establish every pretrained example's domain.

| CRISPRapido representative | Encoder/checkpoint compatibility | Adapter policy |
|---|---|---|
| Resolved 20-base guide, no gap, strict NGG, ≤6 mismatches | 23 columns, left-padded to 24 | Score |
| One-base RNA/I bulge, ≤4 mismatches | 23 columns, left-padded to 24 | Score |
| One-base DNA/D bulge, ≤4 mismatches | 24 columns | Score |
| One-base bulge with 5 mismatches | Numerically encodable, outside documented with-gap preprocessing cap | Unavailable |
| Two-base DNA/D bulge | **25 columns: official encoder rejects it; checkpoint expects 24** | Unavailable; no cropping/projection |
| Two-base RNA/I bulge | Fits 23 columns but outside one-gapped-nucleotide policy | Unavailable, not a claim that the encoder rejects it |
| Adjacent mixed I/D, or multiple gap events | Some fit the tensor but exceed the documented gap budget | Unavailable |
| True reference ambiguity, guide ambiguity, incomplete/ambiguous PAM | Upstream maps some Ns to matches, while training loader excludes target Ns | Unavailable, not an inferred match |
| Non-NGG observed PAM, other guide lengths, partial guides | Not in this integration's validated subset | Unavailable |

Thus **all search settings remain configurable and unchanged**; the scoring layer reports availability per emitted representative. `-m 5 -b 1 -z 2` is **not** reduced to four mismatches/one-base bulges. Sites outside the scoring policy are retained. The neural network has no hard mismatch-count check in its input encoder; our mismatch caps are explicit conservative policy, not a claim that a high-mismatch tensor cannot produce a number. Wider data/PAM/gap-position generalization, and alignment-policy differences from upstream SWOffinder, remain scientifically unvalidated.

The official encoder's two-base DNA rejection was reproduced as **`ValueError: actual_seq_size > seq_len`**. A two-base RNA example produced feature shape **`(1,24,25)`**, but no model score was assigned to it. Supporting 25 columns would require an appropriate model/data validation, not merely changing the encoder's `seq_len` or applying CFD's trimming rule.

## 3. Alignment conversion and independent inference

For every input PAF representative:

1. Use the exact guide identified by PAF query ID, not a substituted guide.
2. Fetch the target and **complete observed PAM** from the same plain reference. For minus rows, reverse-complement the target and genomic PAM on its left. Existing CIGARs are already guide-oriented: **do not reverse the CIGAR**.
3. Replay the original CIGAR: `=`/`X`/`M` consume guide and target; `I` yields `(guide base, '-')`; `D` yields `('-', target base)`. Validate consumption and resolved-base assertions. Never re-align or choose another representative.
4. Append **`NGG`** to the aligned guide and the actual three-base PAM to the target. The structural guide PAM N is distinct from an unknown reference N. Upstream raw identification files also use guide `TargetSequence` values ending `NGG`.
5. Apply the explicit availability policy. For supported inputs, use upstream-equivalent left **gap** padding to 24 columns and categorical pair order `ACGT- × ACGT-`. Pair index is `5*guide_index+target_index`; a guide PAM N maps to the observed target base. A joint padding `('-', '-')` is category **24**, not zero/AA.

[`adapter.py`](adapter.py) implements this independently with an in-memory FAI dictionary and streamed reference hashing, not a whole-reference sequence copy. It requires plain, regularly wrapped FASTA (final lines may be short); gzip/irregular wrapping are not supported in this prototype. The original search backends retain their own input support. Pass a trusted `--expected-reference-sha256` from the search run where available. The supplied reference must be the actual search assembly: PAF metadata alone is not proof of identity. Pre/post hashes reject persistent input/reference/guide changes; they are not a defense against malicious same-user mutate-and-restore races.

[`worker.py`](worker.py) runs once per PAF batch, loads **five models once**, and processes up to 512 rows per chunk by default (`1..4096` configurable). Each Keras predict call uses the upstream internal batch-size setting 4096. It loads checksum-pinned HDF5 with `compile=False`; it does **not** load pickles or import upstream training/auxiliary models. CPU execution is explicit: GPU disabled, TF intra/inter-op, OpenMP and OpenBLAS compute-thread settings one. Auxiliary runtime threads can still exist. Models/manifest are checked around inference. Infrastructure/inference/protocol failures leave an `INCOMPLETE` run, not fabricated unavailable predictions.

## 4. Configurable scoring and output compatibility

This is a **post-search** tool; core Rust `--scoring` integration is deferred. Start from original, unannotated CRISPRapido PAF:

```bash
# Default: byte-preserving CFD-only pass-through; no TensorFlow or reference needed.
python3 research/crispr_bulge/score_paf.py --input original.paf --out-dir NEW_CFD_RESULT

# Optional second score, in an explicit isolated model environment.
python3 research/crispr_bulge/score_paf.py \
  --input original.paf --out-dir NEW_BOTH_RESULT --scoring both \
  --reference reference.fa --guides guides.fa \
  --expected-reference-sha256 TRUSTED_SEARCH_REFERENCE_SHA256 \
  --model-root /path/to/pinned/upstream \
  --model-python /absolute/path/to/isolated-env/bin/python
```

`--scoring` choices:

- **`cfd` (default):** preserve existing upstream CFD values/formatting byte-for-byte. This is not a new CFD implementation or a CFD calculation benchmark.
- **`crispr-bulge`:** report CB scores/statuses and intentionally omit `cf` tags. With this postprocessing design, the upstream CLI has still calculated CFD; suppressing its computation would require a later core reporting integration.
- **`both`:** preserve every original field/tag, including CFD, and append independent CB annotations. There is no fusion, reranking or score-dependent representative selection.

New auxiliary tags only in opted-in modes:

- `cb:f:<value>` — finite five-model mean, nine significant digits; **absent** if unavailable.
- `bu:Z:scored` or a reason such as `alignment_exceeds_24_columns`, `unsupported_multibase_bulge`, `unsupported_mixed_bulge`, `outside_documented_mismatch_domain`, `ambiguous_reference`, etc. Reasons are the first applicable failure, not an exhaustive list.

A true numeric zero is not unavailable. Every input row and its original coordinates/CIGAR/order are retained for unsupported scoring. Duplicate anchors or malformed/mispaired inputs fail rather than selecting an alternative. Already CB-annotated inputs are rejected to avoid mixing model versions. Results are a new exclusive directory containing `output.paf`, `metadata.json`, `status`, and worker logs/metrics; preserve the bundle for model/input provenance. Data publication occurs only after successful checks; failed runs remain incomplete.

`--min-cfd` is **not reimplemented against rounded PAF values** here. Apply it explicitly in the original CLI if desired. Its earlier filtering remains in the supplied PAF and cannot be undone by choosing CB-only output. It is never interpreted as a CB cutoff. Without it, original unfiltered behavior is unchanged.

## 5. Validation actually performed

No search was rerun. Inputs:

- **370 generated alignment fixtures:** both strands, all three resolved mismatch identities at every guide position, all single RNA-gap positions, all four bases at every single DNA-gap position, mismatch+bulge combinations, alternative gap placements, and deliberately unsupported two-base/mixed/multiple gaps, high mismatch counts, ambiguous bases/PAMs, non-NGG, incomplete PAM and non-20-base guides. These are supplied alignment tests, not new biological sites or search-recall results. Fixture CFD tags came from the existing unchanged Rust CFD probe; no new CFD formula was substituted.
- The retained **172-row independent N-panel PAF**, with its exact original guide/reference inputs: 64 guides on small genomic-derived/engineered-N windows. Its earlier eligibility/CFD audit was not replaced by this scoring experiment.

| Input panel | Rows preserved | CB scored | CB unavailable |
|---|---:|---:|---:|
| Synthetic | 370 | 342 | 28 |
| Retained independent panel | 172 | 43 | 129 |
| Total | **542** | **385** | **157** |

Retained-panel unavailability: **75 ambiguous-reference, 20 multi-base RNA bulges, 18 alignments exceeding 24 columns, 16 outside the conservative mismatch-with-gap domain**. All retained-panel original CFD tags (including absence) remain unchanged in `both`; default `cfd` mode is byte-identical on both panels. Availability on this enriched panel is not a genome-wide coverage estimate.

Independent checks:

- A separate bounded in-memory reference/CIGAR reconstruction matched adapter model inputs for all 385 supported records.
- The independently generated categorical encoding was **byte-identical** to official one-hot encoding plus argmax: SHA256 `21f5c5e24568daae2f971a1987e5207fc60b764c7143a77fe540df8f6f75a8c8`.
- [`official_reference.py`](official_reference.py) invoked **unmodified upstream `ensemble_predict`** on identical aligned inputs, with the pinned source/HDF5/pickles verified first. No mocked predictions or patched upstream inference were used. Its reference-only pickle loading is restricted to the pinned trusted files; arbitrary pickle loading is unsafe.
- All **385 predictions** matched: maximum absolute difference **0** at chunk size 512, **2.0116567611694336e-7** at chunk size 32 (tolerance `1e-6`). A singleton mismatch+DNA-bulge case matched exactly. The size-32 run used **13 chunks with only five model loads**.
- The actual annotated PAF tags in **both reporting modes** were independently compared to official predictions: **385 scores per mode**, maximum difference **4.990997060083657e-10** from PAF decimal formatting. All 542 rows, including 157 unavailable rows, were preserved.
- Three 1,024-row repeated-feature benchmark runs also matched official predictions within **1.564621925354004e-7** each. These repeated records measure inference overhead, not additional off-target sites.
- **10 stdlib-only regression tests passed**, including strand conversion, padding/N semantics, mismatch/gap availability, CIGAR failures, reference wrapping/identity, model checksum rejection, duplicate rejection, byte compatibility, unavailable-row preservation without loading a model, and worker failure leaving no published result.

This validates implementation/conversion agreement with the selected official model, **not biological accuracy, generalization, genome-wide recall, or superiority over CFD**. Different valid representative alignments can produce different neural scores, just as for CFD; no biological preference is inferred.

## 6. Runtime and memory overhead

CPU-only shared-host measurements, **three fresh processes per method**, alternating method order, **1,024 identical supported alignment records** per run. Storage caches were uncontrolled/warm; “fresh process” does not mean cold disk. Reference comparisons use the same model environment/weights and unmodified official function.

| Method | Median child wall seconds | Maximum individual child RSS KiB | Maximum sampled summed process-tree RSS KiB |
|---|---:|---:|---:|
| Existing Rust CFD probe | **0.008286481** | **7,052** | **27,396** |
| Isolated CB worker | **5.070746292** | **562,532** | **590,348** |
| Official CB reference wrapper | **6.322246979** | **650,904** | **675,412** |

Wall scopes include process startup, method-specific input/loading/scoring/output. The CB worker also verifies weights and encodes inputs; the official wrapper imports broader dependencies, validates source, builds official features, and probes encoder limits. These are **not isolated scoring-kernel or end-to-end CRISPRapido search timings**, and they should not be advertised as a general speed ratio. The CB method adds substantial cold-start/RAM cost versus CFD. Avoiding unused training imports/categorical one-hot materialization reduced overhead in these runs, not necessarily on every platform.

For the adapter, median TF import was about **1.746797 s**, model loading **0.864421 s**, and encoding all 1,024 rows **0.014337 s**. The first 512-row prediction chunk (including tracing) took a median **1.383432 s**; the second warm chunk **0.251480 s**. The five loaded model instances were reused. These warm phase observations are not a genome-scale throughput guarantee.

One end-to-end annotation per panel/mode took **5.637313 s / 4.627746 s** for synthetic/retained `both` runs, and **5.237102 s / 5.645714 s** for CB-only. Those single calls include driver/reference scan/hash/worker startup and were not repeated for a robust latency estimate. CFD pass-through file copying was approximately 1 ms and is **not** the time to calculate CFD.

RSS was captured by a fresh-process wrapper using Linux child `ru_maxrss`; it is a maximum individual descendant high-water mark, not aggregate unique physical memory. Tree RSS was additionally sampled every **20 ms**, including measurement/supervisor processes; shared pages may be double-counted and short spikes missed. A **4 GiB sampled process-tree stop** and 300 s child timeouts were used for benchmarks; no stop fired. This is risk reduction, not a hard aggregate memory limit.

Streaming keeps the model feature batch bounded; the driver still keeps a set of anchor keys and uses disk spools proportional to input rows. It scans/hashes the entire reference and rehashes it afterward; an unusually long unwrapped FASTA line can still require a large line buffer. Full-assembly I/O was unbenchmarked at this initial milestone; the subsequent report above measures one existing guide/assembly annotation. Huge outputs, irregular FASTA, concurrent writers, GPU inference, many-guide production lifetimes and packaging remain unbenchmarked. No general HPRC performance claim follows.

## 7. Reproducibility, artifacts and stop point

Source baseline: **`178de00`**. Actual external paths:

- upstream source/models: `/export/local/home/salehi/crisprapido-crispr-bulge-review-v1/`;
- isolated environment: `/export/local/home/salehi/crisprapido-bulge-env-v1/`;
- raw validation/benchmark artifacts: `/export/local/home/salehi/crisprapido-bulge-validation-v1/`;
- preserved copies of the small legacy inputs: `retained172-inputs/{reference.fa,guides.fa,input.paf}` under that validation root. The raw run read `/tmp/crispr-nfix-172/{ref.fa,guides.fa}` and `/tmp/crispr-f7-profile-p172-final.paf`; the preserved copies have identical digests.

Small committed evidence: [summary](results/summary.json), [annotated-output parity](results/output-parity.json), [official metrics](results/official-metrics.json), [adapter chunk-512](results/adapter-512-metrics.json), [chunk-32](results/adapter-32-metrics.json), [source/model inspection](results/source-examination.json), [tests](results/tests.log). The validation took **87.569 s** in total, including annotation, reference checks and benchmark processes; it was not a search timing. No large model, environment, FASTA/index or PAF output is committed.

Example reproduction (new external directories only):

```bash
# Explicit setup; scoring itself never downloads a model.
python3 research/crispr_bulge/fetch_models.py --out /NEW/upstream \
  --ca-file /etc/ssl/certs/ca-certificates.crt
python3.11 -m venv /NEW/env
/NEW/env/bin/python -m pip install -r research/crispr_bulge/requirements-reference.txt
# For the exact tested transitive versions, use requirements-validation.lock on Python 3.11.

PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover \
  -s research/crispr_bulge -p 'test_*.py' -v

PYTHONDONTWRITEBYTECODE=1 /NEW/env/bin/python research/crispr_bulge/validate.py \
  --out /NEW/validation --model-root /NEW/upstream --model-python /NEW/env/bin/python \
  --cfd-probe /export/local/home/salehi/crisprapido-cfd-feasibility-tools.XanPvlMD/cfd_probe \
  --panel-reference /export/local/home/salehi/crisprapido-bulge-validation-v1/retained172-inputs/reference.fa \
  --panel-guides /export/local/home/salehi/crisprapido-bulge-validation-v1/retained172-inputs/guides.fa \
  --panel-paf /export/local/home/salehi/crisprapido-bulge-validation-v1/retained172-inputs/input.paf

python3 research/crispr_bulge/check_output.py \
  --validation-dir /NEW/validation --output /NEW/validation/output-parity.json
```

`fetch_models.py` is an explicit sparse/checksummed setup helper, not a runtime dependency installer. It was syntax/help checked; this experiment used its equivalent component Git operations during source review, not a duplicate helper download. The saved Rust CFD probe was built from the unchanged existing CFD implementation; its source matches current `src/cfd_score.rs`. Compiler/source hashes and dependency versions are retained. The validation driver also needs NumPy (the shown model environment supplies it). Reference validation's upstream function writes a fixed scratch `predictions.csv` within the new external validation directory; per-call JSON predictions/metrics are retained separately. No old experiment/report was overwritten.

**Stop point:** optional scoring is isolated and coverage-limited, not promoted into the core Rust CLI. Before core integration, agree on this availability policy and whether complete default-setting coverage is essential. Two-base DNA bulges need an appropriate new model/data strategy; RNA-two-base/mixed/high-mismatch or noncanonical extensions need explicit scientific validation before reporting extrapolated scores. Broader biological benchmarking, licensing/packaging review and full-workflow profiling remain. Do not enable neural-score pruning or transfer CFD bounds to this nonlinear model, change search defaults, select a higher-scoring alternative alignment, or launch whole-genome experiments automatically.
