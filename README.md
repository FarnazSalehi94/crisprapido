# CRISPRapido

![CRISPRapido Logo](crisprapido.png)


CRISPRapido is a reference-free tool for comprehensive detection of CRISPR off-target sites using complete genome assemblies. Unlike traditional approaches that rely on reference genomes and variant files, CRISPRapido directly analyzes haplotype-resolved assemblies to identify potential off-targets arising from any form of genetic variation. By leveraging the efficient Wavefront Alignment (WFA) algorithm and parallel processing, CRISPRapido enables fast scanning of whole genomes while considering both mismatches and DNA/RNA bulges. The tool is particularly valuable for therapeutic applications, where comprehensive off-target analysis is critical for safety assessment. CRISPRapido can process both complete assemblies and raw sequencing data, providing flexibility for different analysis scenarios while maintaining high computational efficiency through its robust Rust implementation.

## Features

- Fast parallel scanning of genomic sequences
- Support for both gzipped and plain FASTA files
- Configurable mismatch and bulge tolerances
- Automatic reverse complement scanning
- PAF-format output compatible with downstream analysis tools
- Multi-threaded processing for improved performance
- CFD (Cutting Frequency Determination) scoring for off-targets
## Installation

You need to build `WFA2-lib` first, which is a submodule of this repository. To do so, run:

```bash
git clone --recursive https://github.com/pinellolab/crisprapido.git
cd crisprapido/WFA2-lib
make clean all
cd ..
```

Then, you can install CRISPRapido using Cargo:

```shell
# Point to your pre-built WFA2-lib directory
export WFA2LIB_PATH="./WFA2-lib"

# Install CRISPRapido
cargo install --git https://github.com/pinellolab/crisprapido.git
```

### For GUIX's users

```bash
git clone --recursive https://github.com/pinellolab/crisprapido.git
cd crisprapido/WFA2-lib
guix shell -C -D -f guix.scm
export CC=gcc; make clean all
exit
cd ..
env -i bash -c 'WFA2LIB_PATH="./WFA2-lib" PATH=/usr/local/bin:/usr/bin:/bin ~/.cargo/bin/cargo install --path .'
```

## Usage

```bash
crisprapido -r <reference.fa> -g <guide_sequence> -p <pam_sequence> [OPTIONS]

# Columba-backed batch mode
crisprapido -r <reference.fa> -G <guides.fa> -p <pam_sequence> \
  --columba-bin <columba> --columba-index <index_prefix> [OPTIONS]
```

The existing Columba + WFA2 search remains the default. An optional,
single-threaded PAM-anchored Columba backend can be selected with
`--search-backend pam-anchored --pam-search-bin <prototype>`. It requires a
verified reference-bound index (freshly built or independently imported) and rejects unsupported options. See the
[experimental backend report](research/pam_anchored_poc/integration/REPORT.md)
for build and usage instructions, index safety, supported settings, validation,
and limitations. Valid alternative alignments can yield different CFD scores.

### Required Arguments

- `-r, --reference <FILE>`: Input reference FASTA file (supports .fa and .fa.gz)
- `-g, --guide <SEQUENCE>`: One guide RNA sequence (without PAM)
- `-G, --guides <FASTA>`: Multi-guide FASTA with unique record IDs; mutually exclusive with `--guide` and requires Columba candidate input

In batch mode, each FASTA record ID is retained as the Columba QNAME and as
the PAF query name. The reference is loaded once, and the complete FASTA is
submitted once at the configured candidate edit bound (`e=3` for `-m 2 -b 1
-z 1`) plus one exact (`e=0`) recovery run when that bound is nonzero. An
existing multi-query SAM can instead be supplied with `--columba-sam`.
### Optional Arguments

- `-p, --pam <DINUCLEOTIDE>`: Required final two bases of a complete three-base PAM (default: `GG`, requiring NGG). Use `--allow-noncanonical-pam` for permissive observed-PAM scoring.

- `-m, --max-mismatches <NUM>`: Maximum number of mismatches allowed (default: 5)
- `-b, --max-bulges <NUM>`: Maximum number of bulges allowed (default: 1)
- `-z, --max-bulge-size <NUM>`: Maximum size of each bulge in bp (default: 2)
- `-f, --min-match-fraction <FRACTION>`: Minimum matching guide fraction (default: 0.75). For a 20-base guide, seven mismatches require `-f 0.65` or lower.
- `-w, --window-size <NUM>`: Size of sequence window to scan (default: 4x guide length)
- `-t, --threads <NUM>`: Number of threads to use (default: number of logical CPUs)
- `--allow-noncanonical-pam`: Keep nonmatching or unavailable PAMs; CFD still uses the actual observed PAM suffix.
- `--no-filter`: Bypass sequence thresholds and PAM eligibility; overlap representative selection still applies. An explicitly supplied `--min-cfd` still applies afterward.
- `--min-cfd <SCORE>`: Optional downstream minimum CFD, finite and in `[0,1]`; disabled by default. Applied **after** existing representative selection at each PAM anchor, using the full-precision score (not the four-decimal PAF tag). Lower-scoring representatives are excluded; unavailable scores are preserved without promoting a different alignment. This does not prune the integrated search or alter mismatch/bulge limits. Omit the flag for unfiltered output.

## Output Format

CRISPRapido outputs results in the Pairwise Alignment Format (PAF), which is widely used for representing genomic alignments. Each line represents a potential off-target site with the following tab-separated fields:

| Column | Field | Description |
|--------|-------|-------------|
| 1 | Query name | `Guide` in legacy single-guide mode; the FASTA record ID in batch mode |
| 2 | Query length | Length of the guide RNA |
| 3 | Query start | 0-based start position in the guide sequence |
| 4 | Query end | 0-based end position in the guide sequence |
| 5 | Strand | '+' (forward) or '-' (reverse complement) |
| 6 | Target name | Reference sequence name (e.g., chromosome) |
| 7 | Target length | Length of the target reference sequence |
| 8 | Target start | 0-based start position in reference |
| 9 | Target end | 0-based end position in reference |
| 10 | Matches | Number of matching bases |
| 11 | Block length | Total alignment block length |
| 12 | Mapping quality | Always 255 for CRISPRapido |

Additionally, CRISPRapido includes these custom tags:

| Tag | Description |
|-----|-------------|
| `as:i` | Alignment score (lower is better) |
| `nm:i` | Number of mismatches |
| `ng:i` | Number of gaps (indels) |
| `bs:i` | Biggest gap size in bases |
| `cg:Z` | CIGAR string representing alignment details |
| `cf:f` | CFD score


### CFD Score

The Cutting Frequency Determination (CFD) score is an empirical relative score, not a calibrated probability of cleavage.
With the supplied matrices it ranges from 0.0 to 1.0, taking into account:

- Position-specific mismatch penalties
- PAM sequence efficiency
- Bulge and gap effects

This implementation requires two data files:

- `mismatch_scores.txt` : Position-specific mismatch penalties
- `pam_scores.txt` : Efficiency scores for different PAM sequences

For example, add `--min-cfd 0.1` to an existing command to retain representatives with CFD at least 0.1 **or an unavailable score**. A missing `cf:f` tag denotes unavailable scoring under the existing output convention; it is not zero. Some pure bulges have supported CFD scores, while mixed I/D alignments or retained reference Ns can be unscorable. Threshold-selected output is not exhaustive under mismatch/bulge settings alone.

An isolated [native CFD-pruning study](research/pam_anchored_poc/cfd_pruning/REPORT.md) is research-only; it is **not enabled in either integrated search backend**.

### Optional CRISPR-Bulge postprocessing (experimental)

A separate [batched scoring postprocessor](research/crispr_bulge/REPORT.md) accepts `--scoring cfd|crispr-bulge|both` (default: `cfd`). These are **postprocessor options, not new Rust CLI flags**. It preserves the supplied representative alignments; `both` preserves CFD tags and appends `cb:f` scores with `bu:Z` availability reasons. Unsupported rows are retained without a fabricated score.

The pinned pretrained model does **not** cover all default `-m 5 -b 1 -z 2` alignments: two-base DNA bulges exceed its 24-column input, and other multi-base/mixed or out-of-domain cases remain unscored. No search setting is reduced to fit the model. CFD calculations, defaults and `--min-cfd` behavior remain unchanged. See the report for model licensing, a dedicated optional TensorFlow environment, exact input conversion, validation and overhead before use.

The [saved HG00735_pat evaluation](research/crispr_bulge/ASSEMBLY_EVALUATION.md) retained all 21,447 records and CFD scores; 4,526 received CRISPR-Bulge scores (21.10%). Annotation took 91.6 s, mostly reference validation; 93 sampled predictions matched official inference. Core CLI promotion remains deferred because of limited coverage and packaging/availability-policy review.

### Example Output

```
Guide   20      0       20      +       chr1    248956422       10050   10070   19      21      255     as:i:6  nm:i:1  ng:i:0  bs:i:0  cg:Z:19=1X  cf:f:0.0549
```

This indicates:
- A 20bp guide RNA aligned to chromosome 1
- Position 10050-10070 on the forward strand
- 19 bases match with 1 mismatch (nm:i:1)
- No gaps (ng:i:0)
- Alignment score of 6 (as:i:6)
- CIGAR string shows 19 matches followed by 1 mismatch

### PAF Format Specification

For more details on the PAF format, see the [official specification](https://github.com/lh3/miniasm/blob/master/PAF.md) from the developers of miniasm.

## Example

```bash
crisprapido -r genome.fa -g ATCGATCGATCG -p GG -m 3 -b 1 -z 2
```

## Testing

Run the test suite:

```bash
# Point to your pre-built WFA2-lib directory
export WFA2LIB_PATH="./WFA2-lib"

cargo test
```

Enable debug output during development:

```bash
cargo run --features debug
```

## License

See LICENSE file

## Citation

Stay tuned!

### PAM/CFD correction and compatibility

By default, `-p GG` now enforces NGG on both strands before overlap selection.
A full three-base PAM with resolved A/C/G/T bases is required. CFD receives its
last two bases after strand orientation. Use `--allow-noncanonical-pam` to keep
other PAMs with their actual observed score; a missing/ambiguous PAM has no CFD
score. `-p` accepts two DNA bases, not the string `NGG`. `-z` controls bulge size.

Earlier local binaries used the wrong adjacent two bases and did not enforce
this eligibility rule. Existing CFD thresholds, candidate counts and coordinate
based analyses may require regeneration from unfiltered candidates. The
permissive flag does not reproduce buggy historical scoring. See
[the fix and validation report](docs/PAM_CFD_VALIDATION.md) for corrected
coordinate conventions, tests, small-subset comparisons and result impact.
