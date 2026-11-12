# SBDS NGG/NRR controlled pilot

**Scoring update:** [saved-pilot NRR reannotation](REANNOTATION.md) supersedes
this report's initial `unsupported_pam` results. No searches repeated. Neural
scores now defined for 42,425 hits (+34,321 noncanonical); all supported
predictions independently match official inference. Remaining 149,623 retain
N/A domain reasons. Noncanonical biological validity remains unestablished.
Full cohort launch is explicitly stopped pending further authorization.

SBDS resumed after confirmation of standard SpCas9 NGG and requested FDA NRR
coverage. Baseline HEAD: `6f2b53f`; core search/scoring binaries and formulas
unchanged. Brunello remains separate and was not launched. No subagent tool
was available.

External root:
`/moosefs/salehi/HPRCv2_CRISPRapido/runs/sbds_ngg_nrr_v1/`.
All eight requested directories exist. Frozen pilot plan/tooling are in
`00_manifest/`; four literal searches and scoring bundles in `01_pilot/`;
union/audits in `04_qc/`. Full-cohort stage directories remain empty.

## Search and complete emitted-hit validation

Guide ACTGAAATCTGTAAGCAGGC; HG00735_pat; m5/b1/z2/f0.75. Native and CLI accept
exactly two A/C/G/T PAM suffix bases, not IUPAC RR. Thus four unchanged searches
use GG, AG, GA, AA. There is no combined native NRR invocation. Unrestricted PAM
search was not needed or launched.

Array 2871759 ran serially on octopus07, **one allocated CPU per task**, 56 GiB
requested memory, 48 GiB sampled-tree stop. Task jobs 2871760/2871761/2871763/
2871759 completed exit 0. All four monitors report completed. Per-stage runtime,
CPU times, pins, checksums, actual target sequences/PAMs, original alignments,
independent CFD reconstruction and explicit neural reasons are retained.

| PAM | Hits / CFD defined | Neural defined | Search wall s | Workflow wall s |
|---|---:|---:|---:|---:|
| NGG | 38,060 | 8,104 | 405.2 | 426.8 |
| NAG | 46,358 | 0 | 416.2 | 429.6 |
| NGA | 42,918 | 0 | 406.2 | 419.1 |
| NAA | 64,712 | 0 | 430.2 | 446.7 |
| **NRR union** | **192,048** | **8,104** | **1,657.9** | **1,722.2** |

All 192,048 rows passed independent reference extraction, CIGAR/coordinate,
match-fraction, literal PAM and CFD checks. Audit defaults remain NGG; the new
optional suffix argument is explicit and recorded. Tests cover all 16 NRR
triplets on both strands and reject noncanonical PAMs under the old default.
The union enforces unique (guide, contig, strand, PAM boundary) identities,
checks all component output hashes and supervisor completion, and emits actual
PAM/classification alongside the complete both-score PAF. No hit was filtered.

All **38,060 NGG hits are preserved**. The NGG component is byte-identical to
the previous completed primary both-score output, SHA256
`595d95fef8146418373a6afcee2fffcd4a9fa803c6a462ae50cb14acdddbc836`.
Its extraction from the disjoint union is construction-preserving, not an
independent combined-engine recall experiment. Additional noncanonical hits:
**153,988**. These are emitted-site soundness checks, not exhaustive recall.

## Scorers: computation is not biological validation

CFD: all 16 local PAM penalties exactly match the separately retrieved CRISPOR
reference CFD table. NGG=1, NAG=0.259259259, NGA=0.069444444, NAA=0. The reference
calculator multiplies by the observed target's final two bases. NAA zero is a
**defined matrix penalty**, not an unavailable-score substitution and not proof
of no activity. Source URL/content digest and table agreement are saved under
`04_qc/reference-*`. The reference calculator describes ungapped CFD; its PAM
factor agreement is not independent biological validation of our bulge extension.
Every pilot alignment was also independently reconstructed with existing CFD
matrices/formula; scores were not silently repaired.

CRISPR-Bulge: the deployed adapter explicitly restricts target PAMs to NGG.
It appends NGG to the guide and retains the actual three-base target PAM.
Unmodified official five-model `ensemble_predict` at pinned upstream commit
3eddcd5bfcaff00b2bdf29425116ec9756ace870 successfully computed four synthetic
matched-protospacer inputs with AGG/AAG/AGA/AAA target PAMs (job 2871762).
This demonstrates **input acceptance/computability**, not biological validity
or full noncanonical adapter/reference agreement. Probe outputs are research
artifacts only, not substituted into pilot predictions.

Accordingly all 153,988 noncanonical rows retain **N/A / unsupported_pam** under
the unchanged conservative policy. NGG N/A counts: 13,882 length>24 columns,
14,376 outside gapped mismatch domain, 1,698 unsupported multibase bulge.
Total neural N/A: **183,944**. Existing NGG official-reference comparisons
remain applicable to the byte-identical reproduced scores; no new all-row
neural reference comparison is claimed. Broader noncanonical validation and an
explicitly labelled experimental scoring policy remain unfinished. Predicted
cleavage scores must never be described as confirmed prime-editing activity.

## Resources and expansion gate

Summed serial workflow: **28.7 minutes**. Maximum sampled tree RSS:
14,300,276 KiB = **13.64 GiB**, not 14.3 GiB. RSS sampling may miss spikes or
double-count shared pages. Caches were uncontrolled; one observation per PAM.
Compressed full PAF+annotations (gzip level 6): **10,190,698 bytes**.

Simple 466-fold scenarios: **222.9 one-CPU job-hours**, **4.75 GB** compressed
PAF+annotations, before index/reference preparation, full audits/extra artifacts,
retries and assembly variation. These are scenarios, not measured cohort costs
or guaranteed scheduling throughput. Four searches repeat binding/snapshot
checks intentionally; no integrity checks were bypassed to save time.

**Full 466-genome rerun NOT launched.** The user authorizes safe expansion, but
its prerequisites are not satisfied:

1. No complete 466-input trusted content/binding ledger exists. Prior recovery
   evidence covers 12 sets; ready pilot bindings cover HG00735 and HG01496.
   The remaining input preparation requires full reference/index comparison,
   not metadata-only approval. Historical EIO is never interpreted as zero.
2. Previously committed draft `prepare_binding.py` is unsafe to reuse without
   repair/testing: closed-file fsync, decompressed/compressed digest confusion,
   and `.open()` on a string binding path. The draft full-run worker also closes
   its merged stream inside chunk copies and uses exclusive status writes for
   repeat updates. Earlier claims of production readiness are withdrawn.
   No such helper was used for this pilot; it used the existing attested binding.
3. Noncanonical neural model acceptance is demonstrated, but the requested
   independent adapter/reference and biological-domain review is incomplete.

Next: repair and failure-test a dedicated SBDS binding/array workflow; freeze
per-input trusted source manifests; prepare the exceptional assemblies and two
references first; validate noncanonical neural feature/prediction agreement on
actual audited alignments without mislabelling extrapolations. Then expand
only passed bindings with per-assembly completion records and separate scoring.
Current pilot task restarts checksum-verify completed work; failed attempts are
retained and require a new Slurm attempt. The collector is the supervisor-gated
completion authority; a stage marker alone is not full experiment approval.

## Historical study and figures

The old historical m2/b1/z1 versus current m5/b1/z2 NGG comparisons remain in
`../SBDS_PILOT.md`: 57/57 eligible old witnesses retained, historical CFD
reconstruction discrepancies, four no-perfect cases, and possible duplicated
HG00321_hap2 representation are not newly resolved by this one-assembly pilot.
Chr7 GGG remains eligible; previously highlighted chr12 GGC and chr17 AGT are
**outside NRR as well as NGG**. NRR does not reinstate all historical non-NGG
hits. The 153,988 additions here isolate PAM expansion at fixed current edit
settings; do not attribute historical edit-domain differences to PAM alone.

No complete-population perfect-match distribution or updated six population
figures is claimed. Those require the 466-input completion/assessability ledger;
unreadable or unfinished inputs cannot be zero-filled. Brunello m4/b1/z1/GG
planning, full-library timing and IMPG remain separate, not launched here.
