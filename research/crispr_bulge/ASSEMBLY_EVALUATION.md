# CRISPR-Bulge on the saved HG00735_pat result

**Subsequent user-approved CLI work:** see [CLI_INTEGRATION.md](CLI_INTEGRATION.md). Optional development CLI integration and pinned-offset reuse are now implemented, with all 4,526 predictions checked against official inference. The original evaluation/deferral below is retained as the `90d958a` milestone, not the current integration status.

## Decision

**The existing postprocessor works on this whole-assembly result, but provides only 21.10% neural-score coverage. It is not ready for unrestricted production CLI promotion or as a replacement for CFD.** It is technically suitable for a future explicitly experimental, opt-in CLI wrapper, provided partial availability is an accepted contract. Before integration, agree on that contract, review dependency/model packaging and licensing, and test CLI failure/timeout/output semantics. Avoid rescanning the full reference redundantly if integrating with already verified reference handling, without weakening identity checks. **No CLI integration was performed; findings are reported first.** Complete default-domain coverage requires a separately validated model/data strategy, not cropping or silently relaxing policy.

This evaluation used existing job **2871283**, not a new search. Both search algorithms, mismatch/bulge defaults, representative selection, CFD formulas, production binaries and original results were untouched. No Slurm submission or main merge occurred. CFD remains default and integrated pruning remains disabled.

## Input and preservation

Source baseline: **`4756527`**, branch `feature/pam-anchored-search-poc`. Evaluated 2026-11-10 on `octopus01`.

- Input: `/export/local/home/salehi/crisprapido-resource-pilot-results/hg00735-optimized-2871283/output.paf`
- Guide: `pilot`, `GCCAGGTACAGGTTCCGTCC`; original search settings `m5/b1/z2/f.75`, strict NGG.
- Exact imported reference: `/export/local/home/salehi/crisprapido-hprc-import-hg00735-20261105-v1/reference.fa` (82 contigs, 3,029,527,752 bases).
- Reference SHA256 checked against the saved search `preverify.json`: `bae8df3c5cf8ca96d3316b59eaa63fa37b7d3b58551c8a0700f2c530dc1e4f8b`.
- Original PAF SHA256 before/after: **`a28d0574a915a57d431851ab681b0602b8f5fb1740d1c8efeb6e3b4a378f23be`**.
- Annotated PAF SHA256: `35a4c0441c90434c05334d5cb97221d83b4364c73ca324d35adedad0f09d15f1`.

The unchanged batched postprocessor ran in **`both`** mode, batch size 512. Every original row, field, tag, CFD value, CIGAR and order was independently checked against the new annotated file. All 21,447 records were retained. The original PAF was only read, never rewritten. Reference faidx indexes for independent checks were written to the external evaluation directory, not beside the original reference. Full outputs, sample sequences and model files remain outside Git.

## Availability

| Score availability | Records |
|---|---:|
| CFD | **21,447** |
| CRISPR-Bulge | **4,526** |
| Both | **4,526** |
| CRISPR-Bulge unavailable | **16,921** |

Neural coverage: **21.1032%**; unavailable: **78.8968%**. Supported records comprise **1,098 ungapped, 1,080 single-base DNA bulges, 2,348 single-base RNA bulges**. These are one guide's reported representatives, not an estimate for arbitrary guides/assemblies or missed sites.

| First unavailable reason | Records | Independent CIGAR classification |
|---|---:|---|
| `alignment_exceeds_24_columns` | **8,176** | All pure **two-base DNA/D bulges**, 25 columns including PAM |
| `outside_documented_mismatch_domain` | **7,624** | All **one-base DNA/D bulges with five mismatches**; checkpoint encodable, outside conservative four-mismatch-with-gap policy |
| `unsupported_multibase_bulge` | **1,121** | All pure **two-base RNA/I bulges**; encodable, outside one-gapped-nucleotide policy |

Unavailable statuses were cross-checked from original CIGAR geometry and mismatch tags separately from the scorer's conversion routine. All alignment representations were parsed successfully; **zero malformed/unsupported CIGARs, mixed I/D bulges, multiple gap events, partial guides, ambiguous-reference/PAM or non-NGG statuses** occurred. Thus “unsupported alignment representations” here means the two-base bulge geometries above, not a parser failure. No unsupported prediction was invented, replaced with zero, truncated, realigned or assigned an alternative representative. Statuses are mutually exclusive first failures; a row can also violate a later policy check.

All original representatives consist of 1,098 ungapped, 8,704 DNA1, 8,176 DNA2, 2,348 RNA1 and 1,121 RNA2 alignments. The official encoder's DNA2 exception remains `actual_seq_size > seq_len`; RNA2 encoding succeeds but is not treated as validated scoring support.

## CFD versus CRISPR-Bulge

Comparison uses only the **4,526 paired supported representatives**, not all 21,447 sites. PAF CFD values have four decimal places; neural annotations have nine significant digits. Ties and zeros in stored CFD affect these descriptive statistics.

| Alignment class | Paired n | Spearman rho | Pearson r | Median CFD | Median CRISPR-Bulge |
|---|---:|---:|---:|---:|---:|
| All supported | 4,526 | **0.4281** | **0.3447** | 0.0052 | 1.2848e-13 |
| Ungapped | 1,098 | 0.5133 | 0.5765 | 0.0512 | 2.9633e-9 |
| One-base DNA | 1,080 | 0.1461 | 0.3299 | 0 | 1.8847e-19 |
| One-base RNA | 2,348 | 0.1571 | 0.2993 | 0.0017 | 1.0708e-13 |

Overall Kendall tau-b: **0.3122**. Paired CFD range: 0–1; neural range: **3.8038e-38–1**. There are **1,742 stored CFD zeros**, which may be actual zeros or rounding; no neural predictions are zero. Quantiles by alignment class are retained in [evaluation.json](results/hg00735-2871283/evaluation.json).

These models produce markedly different numerical distributions and only modest rank association, especially within gapped classes. **Do not transfer thresholds, interpret score differences as errors, fuse the values, or claim direct calibration/interchangeability.** Overall association partly reflects differences among alignment classes. Very small neural outputs are actual model outputs, not unavailable markers. The guide, representative alignment policy, conservative supported subset, rounding and repeated/homologous sites limit generalization. No cleavage labels were analyzed: this is not biological accuracy, superiority, calibration or recall validation.

## Independent official validation

A deterministic **93-record** sample includes up to four records per observed alignment-class/strand/mismatch stratum (seed 2871283), plus the bottom/top ten records under each score, deduplicated. This covers all observed supported strata, both strands, all four observed NGG PAM identities and the exact target. It is intentionally stratified/tail-enriched, not a simple random accuracy estimate.

- Reference targets and PAMs extracted independently with **samtools 1.14 / htslib faidx**, using an external index.
- A separate CIGAR replay implementation reconstructed aligned guide/target strings; checked consumption, strand, `=`/`X` assertions and mismatch counts without calling the adapter's conversion or indexed reference reader.
- Independently encoded categories matched official one-hot-plus-argmax features exactly: SHA256 `6008e1e59cee065d0259f4e9843ec1f4c098d8554bbb53ba266cac201a54caed`.
- Called **unmodified official `ensemble_predict`**, pinned source commit `3eddcd5bfcaff00b2bdf29425116ec9756ace870`, same five `c_2` classification checkpoints, TF-CPU 2.12.0. Source, checkpoint and trusted reference-only pickle hashes were checked.
- **93/93 predictions passed**: maximum absolute difference from actual annotated PAF scores **2.5598145825256324e-10**, below `1e-6` tolerance (decimal/batch precision differences allowed).

All score comparisons preserve record IDs/order. The initial development audit and final strengthened CIGAR-domain audit are retained externally; committed results use **`evaluation-final.json`** from the final audit only. Official predictions are sample-validated, not independently recomputed for every supported row. This establishes implementation agreement, not biological validity or genome-wide recall.

**14 regression tests passed**: ten prior scorer tests plus four saved-result audit/sampling/independent-reconstruction tests. [Test log](results/hg00735-2871283/tests.log).

## Runtime and memory

One measured fresh-process CPU invocation; storage cache state uncontrolled. This is **post-search annotation overhead only**, not a rerun of the original search, not a CFD recalculation benchmark, and not a repeated median.

| Measurement | Value |
|---|---:|
| Postprocessor wall, external wrapper | **91.5757 s** |
| Internal total | 91.5263 s |
| Reference index/validation/hash scan | **81.7841 s** |
| Alignment conversion | 0.7319 s |
| Neural worker including startup/loading/encoding/inference | **6.0402 s** |
| Sum of nine neural prediction calls | **3.3133 s** |
| TF import / five-model loading | 1.7481 / 0.8660 s |
| Maximum individual descendant RSS (worker) | **561,608 KiB (548.45 MiB)** |
| Sampled summed process-tree peak | **617,108 KiB (602.64 MiB)** |

The five models were loaded **once**, reused for eight 512-record chunks and one 430-record chunk. All 21,447 rows incur annotation/conversion work, but only 4,526 incur neural inference. Full-reference scanning validates wrapping/alphabet and computes hashes; its measured cost is not pure disk I/O. The reference is also rehashed after inference. This shows reference validation dominates this invocation; do not extrapolate a universal throughput or speedup.

Tree sampling requested 20 ms intervals (actual intervals include `/proc` sampling overhead); a 4 GiB sampled stop and 900 s driver timeout were used, neither fired. Summed RSS can double-count shared pages/miss spikes and includes supervisors; it is not unique physical memory or a hard limit. Linux `ru_maxrss` is a maximum individual descendant high-water mark, not concurrent summed memory. [Wall/RSS measurement](results/hg00735-2871283/measurement.json), [tree summary](results/hg00735-2871283/tree.json). Independent reference validation is separate from the annotation timing.

## Artifacts and reproduction

External root: `/export/local/home/salehi/crisprapido-bulge-hg00735-2871283-v1/`:

- `both/output.paf`, `both/metadata.json`, worker metrics/logs and complete status;
- `measurement.json`, `tree.json`, `tree.log`, supervisor/driver logs;
- `independent-final/` contains sample requests, region list, reference index, official predictions/metrics and logs;
- `evaluation-final.json` and `tests.log` are the authoritative final audit.

Committed evidence contains only small aggregate JSON/logs, not PAF/FASTA, sampled sequences, models or environments. Evaluation code: [evaluate_assembly.py](evaluate_assembly.py). Existing scoring/inference code was not changed.

To reproduce, choose a **new external ROOT**, create a `guides.fa` containing the exact `pilot` guide above, and retain all result directories:

```bash
ROOT=/NEW/external-evaluation
mkdir "$ROOT"
printf '>pilot\nGCCAGGTACAGGTTCCGTCC\n' > "$ROOT/guides.fa"
INPUT=/export/local/home/salehi/crisprapido-resource-pilot-results/hg00735-optimized-2871283/output.paf
REF=/export/local/home/salehi/crisprapido-hprc-import-hg00735-20261105-v1/reference.fa
UPSTREAM=/export/local/home/salehi/crisprapido-crispr-bulge-review-v1
PYTHON=/export/local/home/salehi/crisprapido-bulge-env-v1/bin/python

PYTHONDONTWRITEBYTECODE=1 python3 research/pam_anchored_poc/monitored_pilot/supervise.py \
  --limit-kib 4194304 --interval .02 --log "$ROOT/tree.log" --summary "$ROOT/tree.json" -- \
  python3 research/crispr_bulge/measure.py --metrics "$ROOT/measurement.json" \
  --stdout "$ROOT/driver.stdout" --stderr "$ROOT/driver.stderr" --timeout 900 -- \
  python3 research/crispr_bulge/score_paf.py --input "$INPUT" --out-dir "$ROOT/both" \
  --scoring both --reference "$REF" --guides "$ROOT/guides.fa" \
  --expected-reference-sha256 bae8df3c5cf8ca96d3316b59eaa63fa37b7d3b58551c8a0700f2c530dc1e4f8b \
  --model-root "$UPSTREAM" --model-python "$PYTHON"

PYTHONDONTWRITEBYTECODE=1 "$PYTHON" research/crispr_bulge/evaluate_assembly.py \
  --input "$INPUT" --scored-dir "$ROOT/both" --reference "$REF" \
  --work "$ROOT/independent-final" --upstream "$UPSTREAM" --model-python "$PYTHON" \
  --output "$ROOT/evaluation-final.json"

PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover \
  -s research/crispr_bulge -p 'test_*.py' -v
```

Requires existing trusted pinned models/environment and samtools; no model download or genome search occurs. Stop here for user review before core CLI work, any broader scoring experiment, model-domain expansion or new search.
