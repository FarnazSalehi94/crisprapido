# Search names, components and executable configurations

## Do not count five independent engines

The word **original** has at least three meanings in this project:

1. Historical standalone WFA2/window source, here inspected at **`46cf652bf878dd48df73de919d77ccb560ca6928`**.
2. The **window branch** of the Columba/WFA2 lineage used as “baseline/original” in the existing manuscript figures/benchmarks (production lineage **`f75f82aa...`**). It is not the same frozen artifact as (1).
3. The divergent **main** lineage at **`cfb15874fdaef082e1a9cf2f83d30e1d10d4252a`**, whose source calls Sassy the “original” backend. This is not the baseline in (2).

Authors must state the exact commit/command whenever using “original CRISPRapido.” `upstream/main` locally resolves to **`05a33a06887c010e07346071ed75dd45a0fe7a10`**, another WFA2 source revision, not an executable synonym for every original-mode result.

| Requested name | Actual role | Concrete configuration in this investigation |
|---|---|---|
| Original CRISPRapido | A versioned **complete pipeline**, not a separate search primitive | Historical `46cf652` window/WFA2 artifact, diagnostic only; old manuscript baseline is separately identified above |
| Sassy-based searching | Approximate-search engine **plus its own filtering/reporter** | Divergent main `cfb1587`, Sassy **0.2.4**, lock revision **`6edc9f325f888af34082b50d771eeb6473acf543`**, no ropebwt3 flag |
| WFA2-based searching | WFA2 is an **aligner/verifier**, not an indexed genome-discovery engine by itself | Current direct-window pipeline: generate overlapping windows → ends-free WFA2 → eligibility/PAM → per-anchor representative → CFD/PAF |
| Columba + WFA2 | Complete pipeline with **indexed candidate generation** and WFA2 verification | Same current CLI + `--columba-bin BIN --columba-index PREFIX`; edit candidates at `m+b*z` plus exact recovery, anchored WFA2/bounded fallback, eligibility/PAM/CFD/PAF |
| PAM-aware Columba without WFA2 | Complete pipeline with **PAM-seeded indexed alignment traversal**, no WFA2 query verification | Same CLI + `--search-backend pam-anchored --pam-search-bin NATIVE --columba-index PREFIX`; canonical native alignment → existing CFD/PAF |

In the current CLI, `columba-wfa2` is the default backend-selector value, **not a separate `window` flag**. Single-guide `-r/-g` without Columba/SAM arguments falls through to window scanning; adding Columba candidate inputs selects the indexed-candidate/WFA2 path. `pam-anchored` explicitly selects the native adapter. The JSON command arrays, not an invented mode name, define each execution.

There are **five executable configurations** in the exploratory records, but only **three current complete-pipeline configurations** with the same attested current executable, standard output and strict PAM contract. Historical and current window modes share an architectural family. Two pipelines use Columba differently; WFA2 is shared by the window/Columba-candidate families. Sassy and ropebwt3 are absent from this feature branch. No branch merge, new engine implementation or algorithm change was made to manufacture a five-mode CLI.

## Exact artifacts and source confidence

[implementations.json](results/implementations.json) records executable paths/SHA256s, associated source/lock object hashes and build limitations. [summary.json](results/summary.json) records every executed command, input parameters and pre/post executable/source checksums.

- **Three current configurations:** `/export/local/home/salehi/crisprapido-research-workflow-v1/target-final/release/crisprapido`, built from the clean **`6fa96d3`** source archive; core source verified identical to **`ce52c15`**. Release build, Rust/Cargo 1.94.1, GCC/G++ 12.3.0, WFA2 **df3dce4** / v2.3, Rust binding **c608c436**. Same executable and score tables in all three modes.
- **Native/builder/Columba:** `pam-search`, `columba-build-v2/columba_build`, `columba-build-v2/columba` beneath that same external root. Unmodified Columba **edcc7ce / 2.0.3**, vanilla/32-bit/alphabet-5/sparse-SA-4. Existing mandatory native/snapshot binding checks remain active.
- **Historical window artifact:** `/export/local/home/salehi/crisprapido/target/release/crisprapido`. Associated tracked source is clean at 46cf652; its checkout currently has WFA2 submodule **2ec2891**, with binding c608c436. **The binary's exact original compiler/library build attestation was not recovered**; do not infer it solely from today's checkout. Its measured bytes are hash-pinned, not overwritten/rebuilt.
- **Sassy diagnostic artifact:** `/moosefs/raid5/salehi/crispr2/crispr-progress/crispr-bwt/target/debug/crisprapido`. Associated tracked source is clean at cfb1587. This is a **debug-path existing artifact without a fresh source-to-executable rebuild attestation**. Source inspection and observed output agree on its exporter shape. Rebuild attempts were isolated: the default Cargo could not read lock format 4; Cargo 1.94.1 offline then lacked the pinned Sassy Git cache. Neither the checkout nor executable was modified. A release/profile-matched, attested build remains TODO before a fair Sassy timing/engine comparison.

Associated historical source is nominal provenance, not cryptographic proof of how an old binary was built. Exact executable hashes still identify precisely what was run. Raw historical outputs are diagnostics, not a claim of current official release quality.

## Equivalent requested limits versus actual semantics

Common resolved-reference stratum: one exact 20-base A/C/G/T guide, `m=3,b=1,z=1,f=.75`, strict observed NGG, one thread. Modern modes explicitly receive all parameters. Historical 46cf652 has no `-p`/`-f` options: strict complete NGG is imposed **only by independent post-output validation**, not by altering its search. For an eligible full 20-base alignment at m3/b1/z1, at least 16 guide positions match, so .75 adds no additional exclusion; this does not repair its query/span/orientation export defects or discovery choices.

Sassy source requires explicit **`--allow-indels`**; otherwise gaps are filtered regardless of b/z. Its search radius is **`m+b`**, not `m+b*z`. At z1 both candidate bounds are 4, but the engine optimizes unit edit cost rather than WFA affine cost and may select different representatives. At z2 the nominal Sassy radius remains 4 while the equivalent maximum edit budget is 5 (or 1 versus 2 in the m0 mixed control). **Z2 Sassy results are not a full-capability comparison under equivalent error budgets.** Increasing another parameter to hide this would change the requested semantics and was not done.

`-p GG` in the Sassy source supplies a CFD PAM argument; it does **not** enforce an observed genomic PAM predicate. Default ambiguity exclusion/IUPAC search semantics also differ from the current N-as-mismatch contract. New common controls therefore use only resolved sequence; no Sassy N-policy equivalence is claimed.

The Sassy reporter emits **eleven mandatory-like fields**, omitting the PAF block-length field, and fixes target end to `start+guide_length` even for gapped CIGARs. Its source also converts unavailable CFD to numeric zero via `unwrap_or(0.0)`; that is a version-specific reporting hazard, not a convention adopted here. The historical 46cf652 reporter lacks CFD and can overrun query coordinates on D operations / expose incorrectly oriented negative CIGARs. All raw records are retained, schema-labelled and independently audited; no missing field, repaired endpoint, reversed CIGAR or invented CFD is silently substituted.

Thus the main concordance uses current standards-valid eligible records only. A separate **historical diagnostic eligible-witness set** reads unambiguous known legacy field positions and validates the supplied original-guide alignment/reference/PAM; it is not standards-valid PAF and must not appear as an equal fifth pipeline in a publication Venn diagram.
