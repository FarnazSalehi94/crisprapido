# Brunello downstream specification — not executed population results

**Full Brunello is still review-gated.** IMPG 0.2.4 availability is not proof of
usable alignment inputs, correct projections, homologous loci or allele counts.
Freeze and checksum its executable, alignment/index manifests, assembly aliases,
commands and mapping thresholds before a bounded mapping validation panel.

## Three distinct identities

1. **Reported site:** exact guide sequence + assembly/contig + strand +
   PAM-adjacent boundary (end on plus, start on minus). Preserve original span,
   CIGAR, guide-oriented target/PAM and both scores/reasons. Broad overlap and
   minus `start-3` are not site keys.
2. **Homologous locus/copy:** unique reciprocal mapping with adequate flanking
   context, coherent orientation and interval coverage. IMPG supplies mapping
   evidence, not an allele definition. Do not transitively merge overlapping
   repeats or paralogous copies. Distinct resolved copies retain distinct locus
   IDs. Ambiguous/split projections remain unresolved.
3. **Sequence allele:** actual complete guide-oriented sequence from frozen
   homologous protospacer/PAM interval boundaries; retain the sequence and a
   checksum-based ID. Homology defines *where* to compare; sequence defines
   *which* allele. Projected intervals can be variable length after indels.
   Preserve the original search-target/PAM sequence and alignment separately.
   A representative-dependent PAF start/CIGAR must not manufacture an allele.
   The interval rule and indel/copy handling require bounded validation before
   counting: do not silently crop/realign search or neural-model inputs.

## Assessability table, including haplotypes with no search hit

Minimum key: `(guide_sequence, locus_id, assembly_haplotype_id)`.

Fields include source/sample/haplotype IDs; query/mapping provenance; mapping
coverage/ambiguity and orientation; locus/copy ID; projected boundaries; complete
observed allele sequence; PAM; original hit IDs; eligibility status; CFD and
neural values with reasons; reference relation; assessability status/reason.

For every locus discovered in **any** haplotype, explicitly query all 464
haplotypes. A missing primary-search record is **not** a missing allele: it can
mean an ineligible PAM/mismatch count, another representative, incomplete
mapping, sequence uncertainty, a genuine deletion, or an execution failure.

Assessability states:

- `resolved_sequence`: complete uniquely assigned locus sequence, including
  alleles outside the search's m/b/z/PAM eligibility domain;
- `confirmed_deletion`: independently resolved flanking/path evidence of locus
  deletion, represented as a deletion state, never numeric score zero;
- `ambiguous_copy`, `incomplete_mapping`, `ambiguous_sequence`, `not_queried`,
  `input_or_run_failure`: unresolved, excluded from frequency denominators and
  reported explicitly with counts/reasons.

For a locus, denominator = distinct haplotypes with `resolved_sequence` or
`confirmed_deletion`; numerator = distinct assessable haplotypes carrying the
specific sequence allele. A confirmed deletion can have its own state frequency.
Each assessable haplotype contributes once per uniquely resolved locus/copy.
No division by 464 by default, no denominator restricted to hit-emitting or
CFD-filtered haplotypes, no duplicate-row counting. Report numerator, denominator,
missing fraction and confidence/uncertainty conventions, not frequency alone.
The 464 assemblies are not 464 independent people; paired haplotypes/sample
relationships and sampling ascertainment limit population interpretation.

## Reference absence versus failure to map

Keep separate:

- mapped reference locus with a search hit;
- mapped reference locus with a resolved but search-ineligible allele;
- **confirmed reference absence** with adequate context/deletion evidence;
- unresolved reference relationship (ambiguous/unmapped/insufficient coverage).

No reference hit does not mean off-reference. Unresolved mapping does not mean
confirmed reference absence. For confirmed off-reference loci, use non-reference
homology groups and actual sequences; never invent a GRCh38 coordinate/allele or
reference CFD of zero. GRCh38/CHM13 are comparators, not additional haplotypes in
the 464-denominator cohort.

## Extended CRISPRme-style Figure 1C

Proposed outputs: (i) allele/site CFD prioritization plot ordered within guide,
with frequency as a colour/size or companion panel, (ii) allele-frequency versus
CFD, and (iii) paired `CFD_alt - CFD_reference` when both values are defined and
comparably aligned/scored. Define the two comparison axes independently; do not
fuse CFD and neural scores or transfer a CFD cutoff to neural output.

Use distinguishable symbols/legend categories for reference-represented,
confirmed off-reference and unresolved relation. Unresolved rows cannot receive
a confirmed-off-reference colour. Missing CFD has a separate unavailable panel,
not a plotted zero. Neural coverage/distributions are a separate panel with
explicit N/A reasons and supported-hit denominators. If valid representatives
for the same allele yield different CFD values, retain/report that ambiguity;
do not select the maximum as a hidden biological-policy change.

Before final figures, validate bounded fixtures for an SNV, indel, PAM loss,
resolved non-hit allele, true deletion, duplicate copy, ambiguous projection,
confirmed off-reference locus and partial haplotype coverage. Quantify missing
mapping/score data and assert numerator/denominator conservation and one-count
per haplotype/locus. No such whole-library mapping/frequency result is claimed
by the initial audit or the three-guide search pilot.
