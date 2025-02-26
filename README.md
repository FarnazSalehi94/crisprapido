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
cd trace_pocrisprapidoints/WFA2-lib
guix shell -C -D -f guix.scm
export CC=gcc; make clean all
exit
cd ..
env -i bash -c 'WFA2LIB_PATH="./WFA2-lib" PATH=/usr/local/bin:/usr/bin:/bin ~/.cargo/bin/cargo install'
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

See LICENSE file

## Citation

Stay tuned!
