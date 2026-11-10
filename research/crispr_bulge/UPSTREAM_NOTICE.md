# CRISPR-Bulge source/model attribution

Official repository: https://github.com/OrensteinLab/CRISPR-Bulge

Pinned commit: `3eddcd5bfcaff00b2bdf29425116ec9756ace870`.

The upstream root [MIT license](UPSTREAM_LICENSE), copyright (c) 2024 Orenstein Lab, is retained here for attribution. The checkpoint paths/hashes in `model_manifest.json` identify the five GRU-Embd (`c_2`) classification ensemble components recommended by upstream `main_predict.py`. Weights, upstream Python code, pickles and datasets are **not vendored into CRISPRapido**. Download them separately from the pinned official repository. No separate checkpoint-specific terms were found in the reviewed paths; this is not a legal certification of all underlying data/dependency rights.

The adapter implements input reconstruction/categorical encoding independently and validates it against the unmodified upstream prediction function. The reference-only validation imports upstream code and loads its checksum-pinned trusted pickles. The inference worker uses checksum-pinned HDF5 models without unpickling or importing upstream training/auxiliary model code.

Upstream cites: *Generating, modeling, and evaluating a large-scale set of CRISPR/Cas9 off-target sites with bulges*, preprint DOI [10.1101/2023.11.01.565099](https://doi.org/10.1101/2023.11.01.565099). Archive linked by upstream: [10.5281/zenodo.10902445](https://doi.org/10.5281/zenodo.10902445). Model agreement here is implementation validation, not independent confirmation of the paper's biological performance.

TensorFlow/Keras and other dependencies retain their own licenses. The full dependency/data redistribution audit required for packaging a production release is outside this bounded integration; no dependency environment is committed.
