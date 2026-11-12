# NRR neural scoring: saved SBDS pilot reannotation

Baseline `ef99442`; **no searches repeated, no cohort launched, no production
binaries changed**. Only the Python scoring policy/encoder changed. Existing
attested development CLI invoked `--score-paf --scoring both` with an explicitly
selected frozen new scorer source directory. CFD tags are preserved from the
previous independently audited searches, not recomputed by `--score-paf`.

## Restriction and input policy

`adapter.convert()` and `encode_tokens()` previously required target suffix GG.
This was our conservative integration restriction, not a hard checkpoint or
official encoder restriction. Both checks now allow GG/AG/GA/AA, and still reject
other PAMs and ambiguous/incomplete PAMs. Guide-side aligned spacer always gets
**NGG**; target-side aligned sequence gets its **actual three-base PAM**. No PAM
replacement, cropping, realignment, new weights or retraining. Both strands are
reconstructed in guide orientation. Metadata explicitly states that NRR support
is computational, with noncanonical biological validity not established.

The remaining conservative domain is unchanged: resolved 20-base guide;
ungapped <=6 mismatches; one RNA/DNA gapped base <=4 mismatches; <=24 total
aligned columns. Higher mismatch counts, multibase/mixed/multiple gaps remain
unavailable even when an encoder could technically accept some of them.

## Results (all 192,048 original rows retained)

| PAM | Hits / CFD defined | Neural defined | N/A | >24 columns | Gapped mismatch-domain | Multibase RNA |
|---|---:|---:|---:|---:|---:|---:|
| NGG | 38,060 | 8,104 | 29,956 | 13,882 | 14,376 | 1,698 |
| NAG | 46,358 | 10,452 | 35,906 | 17,380 | 16,276 | 2,250 |
| NGA | 42,918 | 9,421 | 33,497 | 15,706 | 15,247 | 2,544 |
| NAA | 64,712 | 14,448 | 50,264 | 23,506 | 23,277 | 3,481 |
| **Total** | **192,048** | **42,425** | **149,623** | | | |

**34,321 additional noncanonical hits receive neural scores.** Unavailable
scores remain N/A with reasons, absent numeric PAF neural tags, never zero.
NAA's original defined CFD zero is unchanged and independent of neural scoring.

## Independent validation

Scoring-only Slurm **2871764**, one CPU, exit 0, supervisor completed:

- Compared **every supported pilot prediction (42,425)**, not just a sample,
  with unmodified official five-model `ensemble_predict` at upstream commit
  `3eddcd5bfcaff00b2bdf29425116ec9756ace870` with verified source/model hashes.
- Inputs reconstructed from the already independent samtools/CIGAR audit were
  compared row-for-row with new annotations: IDs, coordinates, strand, aligned
  sequences and actual PAMs. Official input uses those independent alignments.
- All categorical features exactly match official encoding (SHA256 equality).
- Maximum absolute prediction difference **2.2137804818467277e-9**, below 1e-6.
- Every original PAF field/tag and row retained. Default CFD `--score-paf`
  pass-through is byte-identical for all four files. NGG both-score output is
  byte-identical to the previous pilot; search and CFD sources are unchanged.

Synthetic scoring-only job **2871765**, exit 0: **8,272** fixtures cover all
16 NRR triplets, ungapped 0..6 mismatches, RNA/DNA single-base gaps with 0..4
mismatches, every feasible gap position and every inserted DNA base. Local
worker versus official predictions agree to **6.110667527536862e-13**; categorical
features match exactly. Separate adapter regression tests cover both strands,
actual PAM retention, mismatch limits, DNA2/RNA2 N/A and non-NRR rejection.
This is bounded domain/implementation validation, not every possible sequence.

**Biological limitation:** neither exact implementation agreement nor a model
producing a number establishes NAG/NGA/NAA calibration or experimental
prime-editing activity. These noncanonical predictions are experimental.

## Runtime and memory

Separate annotation walls (including reference identity checks/model startup):
NGG **18.0 s**, NAG **19.0 s**, NGA **19.0 s**, NAA **23.0 s**.
Full job including all four annotations, default-output checks, independent
all-supported reference comparison and input verification: **103.7 s**;
sampled summed-tree peak **990,340 KiB (0.944 GiB)**. One observation per PAM,
uncontrolled caches; RSS sampling can miss spikes/double-count shared pages.
Synthetic domain validation job elapsed 18 s.

## Artifacts and readiness

Root `/moosefs/salehi/HPRCv2_CRISPRapido/runs/sbds_ngg_nrr_v1/`:
- `03_scoring/nrr-enabled-v1/{NGG,NAG,NGA,NAA}/`: new both-score PAFs,
  annotations, metadata, timings; previous pilot outputs untouched.
- Same root: frozen scorer/plan, all independent requests/predictions and
  summary; `COMPLETE.json` is the post-supervisor checksum authority.
- `04_qc/nrr-domain-v1/`: synthetic features/predictions/agreement.

The original inner-stage summary inadvertently included a live Slurm log digest;
finalization verifies every scientific artifact and records the finished log's
hash separately. No scientific result or frozen tooling was overwritten. The
committed worker excludes live scheduler logs from future stage manifests.

The **NRR scoring component is computationally validated for this pilot and
ready for use with the pinned runtime/source and explicit validity caveat**.
This does not certify the complete 466-input workflow: input binding ledger and
safe expansion tooling remain separate blockers. Per user instruction, stop
here; do not launch the complete SBDS analysis. Brunello remains untouched.
