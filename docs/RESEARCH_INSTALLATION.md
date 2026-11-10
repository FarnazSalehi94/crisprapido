# Reproducible research installation

This is the **development research branch**, not a production release. Use a new checkout and external build/output directories; never rebuild over installed production executables. Record the exact Git commit and retain the source bundle, lockfiles, trusted model files, index bindings and run metadata. The tested search/scoring code baseline is **`6fa96d3`**; the subsequent bounded validation adds tests/documentation, not search changes.

See [end-to-end findings and release blockers](../research/reproducibility/REPORT.md), [CLI scoring/output contract](../research/crispr_bulge/CLI_INTEGRATION.md), and [project handoff](PROJECT_HANDOFF.md).

## Tested stack, not a universal compatibility matrix

| Component | Tested version / contract |
|---|---|
| OS / architecture | Linux x86_64, shared `octopus01` host; not a clean OS/container |
| Rust / Cargo | **1.94.1 / 1.94.1**, edition 2021, repository `Cargo.lock`; minimum supported Rust version not established |
| C / C++ / build | **GCC/G++ 12.3.0**, GNU Make, **CMake 3.25.1**; use one coherent compiler/header/linker/libc stack |
| WFA2-lib | Submodule **`df3dce4b99c37ac4f34fdcef74b774c80b7fadc4`**, version file `v2.3`; static library, OpenMP runtime needed by Rust link |
| Rust WFA bindings | `lib_wfa2` **`c608c436a6753d2c21c97d9f5c338efae99d042b`**, pinned in Cargo files |
| Columba | **`edcc7ce2a1e28ded6784a44afbd8b9348a4805ad`**, vanilla **2.0.3**, non-RLC, **32-bit**, alphabet 5, sparse-SA-4 index layout |
| Columba dependencies | **libsais 2.8.5**, **fmt 10.2.1**, parallel-hashmap **`e5b892baed478513adcb6425773cae1eda033057`**, zlib **1.2.13**, pthreads |
| Python research driver | Python 3.9+ recommended; validation driver tested on **3.12.0** (see installation evidence for exact interpreter) |
| Optional inference | **CPython 3.11.15**, **tensorflow-cpu 2.12.0**, **Keras 2.12.0**, **NumPy 1.23.5**, **h5py 3.16.0**; CPU, one compute thread |
| Pretrained model | Five-component **GRU-Embd `c_2` classification ensemble**, upstream **`3eddcd5bfcaff00b2bdf29425116ec9756ace870`**, hashes/paths in [model_manifest.json](../research/crispr_bulge/model_manifest.json) |

Do not substitute a different Columba library/index layout just because its command-line version looks similar. The prototype uses internal Columba APIs. `-mpopcnt` is used by its build; TensorFlow wheels also require a compatible CPU. Other architectures, newer TensorFlow/Python versions, GPU execution and broader libc/compiler combinations are not certified here.

**CFD needs no neural Python environment or TensorFlow.** However, the **PAM-anchored backend requires standard-library `python3` on PATH in all scoring modes**, for snapshot/reference/index checks. The native build also pins an absolute Python interpreter for its binding verifier. Saved-file CFD pass-through and the ordinary Rust CFD scorer do not invoke neural Python.

## 1. Fresh source and Rust/WFA2 build

```bash
git clone --branch feature/pam-anchored-search-poc --recurse-submodules \
  https://github.com/FarnazSalehi94/crisprapido.git
cd crisprapido
# For exact repetition: checkout your recorded commit, then update its submodules.
git rev-parse HEAD
git submodule status
SRC=$(pwd)
OUT=/absolute/NEW/research-build
mkdir "$OUT"

# Build only the required WFA libraries, not optional benchmark/example executables.
# Run this only in the NEW source checkout.
make -C "$SRC/WFA2-lib" -j1 BUILD_TOOLS=0 BUILD_EXAMPLES=0 \
  CC="$(command -v gcc)" CXX="$(command -v g++)"
export WFA2LIB_PATH="$SRC/WFA2-lib"  # must be absolute
export CARGO_TARGET_DIR="$OUT/target"
cargo build --release --locked
cargo test --release --locked -- --test-threads=1
```

`--offline` can be added if the exact Cargo dependencies are already cached. The recorded clean build used an **empty Cargo target** and trusted existing Cargo download/source caches; it was not a network bootstrap of every dependency. Do not replace `Cargo.lock` or use an unrelated installed WFA library.

### Shared Guix/system installations

Do not mix `/usr/bin/cc`/G++ 8 headers with Guix GCC 12 headers/libraries. Clean inherited `CPATH`, `C_INCLUDE_PATH`, `CPLUS_INCLUDE_PATH`, `LIBRARY_PATH` and incompatible `RUSTFLAGS`, then use explicit coherent `CC`, `CXX` and CMake compiler paths. Sourcing a Guix profile may require temporarily disabling shell `nounset`.

On this host, the successful fresh build used Guix GCC 12.3 plus glibc 2.35 development/static compatibility libraries (including `librt.a`), explicit Rust linker/runtime paths, and a new target directory. The first attempts exposed mixed headers/linker paths and missing `-lrt`; the corrected build passed all Rust tests. Exact site commands remain under the external run root in the report. Those `/gnu/store` paths are **not portable installation instructions**. A versioned portable Guix/container environment remains a release blocker.

## 2. Columba libraries, unpruned PAM search and bound indexes

```bash
git clone --no-checkout https://github.com/biointec/columba.git "$OUT/columba"
git -C "$OUT/columba" checkout --detach edcc7ce2a1e28ded6784a44afbd8b9348a4805ad
git -C "$OUT/columba" submodule update --init --recursive

cmake -S "$OUT/columba" -B "$OUT/columba-build" \
  -DRUN_LENGTH_COMPRESSION=OFF -DTHIRTY_TWO=ON -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER="$(command -v gcc)" -DCMAKE_CXX_COMPILER="$(command -v g++)" \
  -DCMAKE_CXX_FLAGS='-include limits -include cstdint'
cmake --build "$OUT/columba-build" -j1

COLUMBA_DIR="$OUT/columba" COLUMBA_BUILD="$OUT/columba-build" \
CXX="$(command -v g++)" PYTHON="$(command -v python3)" OUTPUT="$OUT/pam-search" \
  bash research/pam_anchored_poc/build.sh

# Only a small bounded/new reference here; choose a NEW prefix.
python3 research/pam_anchored_poc/compatibility/index_manifest.py build \
  --builder "$OUT/columba-build/columba_build" --reference bounded.fa --index "$OUT/bounded-index"
```

CMake fetches the dependency versions above; network access/submodule initialization is required unless trusted source caches are supplied. The tested fresh C++ builds used clean local clones of pinned source caches via CMake `FETCHCONTENT_SOURCE_DIR_*` overrides, **not reused static libraries or object files**. The standard network CMake recipe above follows upstream configuration but was not independently network-bootstrapped in this run. Extra forced standard-library includes match the existing native build contract; no Columba algorithm was patched.

Use the reference **binding manifest**, not a raw unbound builder output. It verifies original N-containing reference bytes, component hashes and text/metadata. Existing HPRC indexes must use the documented verified import process, not unnecessary rebuilding. No whole-assembly example is authorized by these installation instructions. Never use the separate research pruning binary as `--pam-search-bin`.

## 3. Optional neural environment and model

Use a separate interpreter. Do not upgrade project/system Python packages in place.

```bash
python3.11 -m venv "$OUT/inference-env"
"$OUT/inference-env/bin/python" -m pip install --require-hashes \
  -r research/reproducibility/requirements-inference-linux-py311.lock
"$OUT/inference-env/bin/python" -m pip check
```

The hash lock contains the **46 tested wheels** for Linux x86_64/CPython 3.11, including transitive dependencies and installer packages. It is not a cross-platform lock or a security/license certification. A fresh venv installed from the original [top-level inference requirements](../research/crispr_bulge/requirements-inference.txt) passed; a second fresh venv installed the captured hash lock **offline** from a wheelhouse, had the identical frozen package set, and produced a byte-identical complete 172-row workflow result.

To prepare/reuse an external wheelhouse:

```bash
"$OUT/inference-env/bin/python" -m pip download --only-binary=:all: --require-hashes \
  -r research/reproducibility/requirements-inference-linux-py311.lock --dest "$OUT/wheelhouse"
# In another NEW venv, no package-index/network access is needed:
/absolute/new-env/bin/python -m pip install --no-index --find-links "$OUT/wheelhouse" \
  --require-hashes -r research/reproducibility/requirements-inference-linux-py311.lock
```

**Model setup is explicit and separate.** Runtime inference never downloads anything. The selected five HDF5 files total **13,912,000 bytes**; retain their paths/checksums from the manifest. Use only trusted upstream bytes and a compatible model root. Never load arbitrary pickle files. Ordinary inference loads hash-checked HDF5 with `compile=False`; only independent official-reference validation uses trusted checksum-pinned upstream pickles.

The setup helper is:

```bash
python3 research/crispr_bulge/fetch_models.py --out /NEW/model-checkout \
  --ca-file /path/to/trusted/ca-bundle.crt
```

**Known installation blocker:** fresh network execution of this helper failed on this host (Git 2.41.0 internal partial-fetch/config abort; a Git 2.20.1 retry failed HTTP/RPC negotiation). TLS checks were not disabled. These failures remain unresolved; do not claim that a full clean network installation is certified. This validation copied the already reviewed five checkpoints into a new external directory and verified all expected SHA256s before executing both clean environments. A trusted locally provisioned model root is currently the tested research path. A robust, tested model-distribution/setup path is required before release.

CRISPR-Bulge's root upstream license is MIT with attribution retained in [UPSTREAM_NOTICE.md](../research/crispr_bulge/UPSTREAM_NOTICE.md); WFA2-lib's checked-out license is also MIT. **Columba's checked-out root license is AGPL-3.0**, and this repository's [existing LICENSE](../LICENSE) declares an academic AGPL/commercial dual-license arrangement. Do not assume the model's MIT license covers the entire linked stack. Review linking/distribution/service obligations and checkpoint/dataset/dependency rights before packaging or deployment; this is not a legal certification or a new license grant. Weights and environments are not vendored into this repository.

Independent official validation additionally needs the pinned source checkout and [reference requirements](../research/crispr_bulge/requirements-reference.txt) (pandas/scikit-learn/SciPy/XGBoost/CatBoost, etc.). That larger environment was reused, not freshly reinstalled here. Inference does not need those training/auxiliary imports.

## 4. Run the complete bounded workflow

```bash
# CFD default: no neural environment, no implicit score threshold.
"$OUT/target/release/crisprapido" --search-backend pam-anchored \
  --pam-search-bin "$OUT/pam-search" --columba-index "$OUT/bounded-index" \
  -r bounded.fa -G guides.fa -t 1 \
  --mismatch-scores "$SRC/mismatch_scores.txt" --pam-scores "$SRC/pam_scores.txt" > cfd.paf

# Same inputs/search/settings + both scores where available.
"$OUT/target/release/crisprapido" --search-backend pam-anchored \
  --pam-search-bin "$OUT/pam-search" --columba-index "$OUT/bounded-index" \
  -r bounded.fa -G guides.fa -t 1 \
  --mismatch-scores "$SRC/mismatch_scores.txt" --pam-scores "$SRC/pam_scores.txt" \
  --scoring both --scoring-output-dir "$OUT/NEW-scores" \
  --crispr-bulge-python "$OUT/inference-env/bin/python" \
  --crispr-bulge-model-root /trusted/verified/model-root > both.paf
```

Use `--scoring crispr-bulge` for the neural-requested mode. For compatibility it **also retains upstream CFD calculation/tags**, like `both`; it is not a pure neural-only kernel benchmark. Default CFD PAF bytes remain unchanged. No scoring selection activates filtering or pruning. `--min-cfd` is a separate explicit search-time threshold, disabled by default; omit it to retain every reported representative.

Neural stdout is PAF with `cb:f` only when defined and `bu:Z` for availability. `NEW-scores/scoring/annotations.tsv` contains guide/target/aligned sequences, observed PAM, original coordinates/CIGAR/mismatch/bulge information, both scores and explicit **`N/A`** reasons. Missing is not zero. Every hit survives unsupported model configurations. Child/infrastructure failures leave `INCOMPLETE`, not fabricated N/A predictions; preserve bundles and check complete status. See the [full output/failure contract](../research/crispr_bulge/CLI_INTEGRATION.md).

Always supply the exact original reference/guides. Score tables default to paths relative to the current working directory; explicit absolute table paths prevent accidental missing-table scoring. The default script path embeds the source checkout: retain it or override `--scoring-script` with a matching complete script bundle. Pinned offset reuse can amortize reference validation but still requires trusted index/reference digests and full pre/post reference hashing; see the scoring guide.

## 5. Repeat the validation, not whole-genome experiments

The [bounded reproduction script](../research/reproducibility/reproduce.sh) creates new tiny bound indexes from the existing fixture generator, compares independent expected anchors/representatives/CFD, tests all three live scoring modes, and runs relevant Python/binding/PAM regressions. It accepts **no arbitrary genome reference**. Three fresh invocations per mode are benchmarked on the same N-containing bounded input.

```bash
export CRISPR_CLI="$OUT/target/release/crisprapido"
export CRISPR_OLD_CLI=/absolute/preserved/preintegration/crisprapido
export PAM_SEARCH="$OUT/pam-search"
export COLUMBA_BUILDER="$OUT/columba-build/columba_build"
export MODEL_PYTHON="$OUT/inference-env/bin/python"
export MODEL_ROOT=/trusted/verified/model-root
# Optional pinned retained172-inputs directory and official model comparison:
export RETAINED_PANEL=/absolute/retained172-inputs
export OFFICIAL_PYTHON=/absolute/reference-env/bin/python
export OFFICIAL_ROOT=/trusted/pinned/full-upstream-checkout
bash research/reproducibility/reproduce.sh /NEW/absolute/external-results
```

`CRISPR_OLD_CLI` is a preserved preintegration development binary (the recorded run used `178de00`). If unavailable, build that commit in another fresh checkout/worktree with a separate target directory and the same compatible WFA/toolchain; do not overwrite an installed binary or claim byte compatibility without a baseline.

The optional retained panel's three exact file hashes are enforced; omit that variable if unavailable. The three generated panels remain independently oracle-checkable and portable. Both official variables must be supplied together; omit both to skip official prediction validation (the summary records that omission). The recorded authoritative run included retained172 and official inference. Rust tests are a separate build step above; the script does not silently deploy/build arbitrary executables.

## Research scope and release gate

The bounded workflow recovered all expected anchors with no extras/duplicates and kept unsupported scores as N/A. That establishes **finite oracle recovery and implementation agreement**, not general genome-wide recall, biological accuracy, calibration or score interchangeability. Model domain remains resolved 20-base NGG guides with at most one gapped nucleotide under conservative mismatch caps; DNA2 exceeds 24 columns. Neural modes require plain regularly wrapped FASTA and add startup/RAM/disk costs. Experimental search remains single-threaded/limited-layout, and extreme repeats/huge outputs, concurrent writers, packaging/relocation, broad biological data and heterogeneous platforms remain unvalidated. Unsigned manifests/checksums assume trusted preparation and cannot defeat malicious mutate-and-restore races.

Do not consider a production release until model distribution/rights, a portable coherent environment, relocation/install tests, coverage/availability policy, resource/lifecycle hardening and independent broader scientific validation have been reviewed. Do not change defaults, enable integrated CFD-guided pruning, merge main or repeat whole-genome experiments automatically.
