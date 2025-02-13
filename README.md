# CRISPRapido

![CRISPRapido Logo](crisprapido.png)

A fast CRISPR guide RNA off-target scanner written in Rust. CRISPRapido uses the WFA2 (Wavefront Alignment) algorithm to efficiently identify potential off-target sites for CRISPR guide RNAs in genomic sequences.

## Features

- Fast parallel scanning of genomic sequences
- Support for both gzipped and plain FASTA files
- Configurable mismatch and bulge tolerances
- Automatic reverse complement scanning
- PAF-format output compatible with downstream analysis tools
- Multi-threaded processing for improved performance

## Installation

Requires Rust nightly toolchain due to WFA2 dependency. Install using:

```bash
# Install from GitHub repository using nightly toolchain
cargo +nightly install --git https://github.com/pinellolab/crisprapido.git
```

## Usage

```bash
crisprapido -r <reference.fa> -g <guide_sequence> [OPTIONS]
```

### Required Arguments

- `-r, --reference <FILE>`: Input reference FASTA file (supports .fa and .fa.gz)
- `-g, --guide <SEQUENCE>`: Guide RNA sequence (without PAM)

### Optional Arguments

- `-m, --max-mismatches <NUM>`: Maximum number of mismatches allowed (default: 4)
- `-b, --max-bulges <NUM>`: Maximum number of bulges allowed (default: 1)
- `-z, --max-bulge-size <NUM>`: Maximum size of each bulge in bp (default: 2)
- `-w, --window-size <NUM>`: Size of sequence window to scan (default: 4x guide length)
- `-t, --threads <NUM>`: Number of threads to use (default: number of logical CPUs)
- `--no-filter`: Disable all filtering (report every alignment)

## Output Format

Output is in PAF format with custom tags, including:
- Standard PAF columns for position and alignment information
- `as:i`: Alignment score
- `nm:i`: Number of mismatches
- `ng:i`: Number of gaps
- `bs:i`: Biggest gap size
- `cg:Z`: CIGAR string

## Example

```bash
crisprapido -r genome.fa -g ATCGATCGATCG -m 3 -b 1 -z 2
```

## Building from Source

```bash
git clone https://github.com/yourusername/crisprapido
cd crisprapido
cargo build --release
```

## Testing

Run the test suite:

```bash
cargo test
```

Enable debug output during development:

```bash
cargo run --features debug
```

## License

[Add your chosen license here]

## Citation

[Add citation information if applicable]
