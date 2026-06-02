# Entropy 0.1.0

🚀 Initial release of the Entropy package!

- Added `run_entropy()` for gene expression entropy estimation
- Supports Seurat objects with or without neighbor list

# Entropy 0.1.1
🛠 Minor feature updates:

- run_entropy() now supports:

- npca parameter: control number of PCA dimensions used for neighbor graph construction

- nfeatures parameter: control number of variable features selected for PCA

- nn.name parameter: optionally specify the neighbor graph name already stored in the Seurat object

Improved flexibility for entropy estimation in Seurat objects
Updated function documentation and vignettes


# Entropy 0.1.2
🚀 Speed improvements and feature enhancements:

- Major speed-up of entropy estimation functions (vectorization, precomputing constants)
- `nfeatures = 0` now skips the selection of highly variable genes (uses all features)
- Improved flexibility in `run_entropy()`
- Updated function documentation and vignettes


# Entropy 0.1.3
  Speed improvements

# Entropy 0.1.4
- Added binary entropy mode for accessibility data (e.g. ATAC-seq)
- Supports Bernoulli (0/1) entropy across cell neighborhoods
- New `mode` argument in `run_entropy()` to switch between count-based and binary entropy

# Entropy 0.1.5
- Optimized binary entropy mode for accessibility data (e.g. ATAC-seq)
- Binary entropy now uses sparse neighborhood 1-counts with a precomputed entropy lookup table
- Added internal binary_entropy_lookup resource loaded in binary mode, with safe fallback for larger neighborhood sizes


# Entropy 0.1.6
- Updated neighborhood input API: replaced legacy `nn_list` / `nn.name` with the new `neighbors` + `neighbors_method` interface
- Added support for computing neighborhoods from a different assay via the `neighbor_assay` parameter
- Added mean entropy per cell automatically to Seurat metadata (`mean_entropy_{assay}`)
- Added vignette explaining neighborhood definition (`auto`, `seurat`, and `custom` modes)
- Improved documentation and pkgdown website structure


# Entropy 0.2.0

## Performance improvements

* `run_entropy()` now supports a chunking strategy (`chunk_method = "auto"`) that groups
  cells into biologically similar chunks using Leiden clustering. This reduces memory
  pressure and computation time significantly for large datasets.
* Added parallelization support via the `n_cores` argument in `run_entropy()`.



# Entropy 0.2.1

### New features

* `run_entropy()` gains a `features` argument: pass a character vector of gene
  or feature names to restrict entropy computation to a subset. Defaults to
  `NULL` (all features, fully backward compatible). Invalid or missing names
  raise an informative error.

* `run_entropy()` gains `reduction` and `dims` arguments, replacing `npca`.
  `reduction` lets you point neighbor and chunk computation at any existing
  dimensional reduction (e.g. `reduction = "lsi"` for ATAC). `dims` accepts
  either a single integer (`dims = 20` → `1:20`) or a vector (`dims = 2:20`
  to skip the first LSI component). Both the auto-neighbors path and chunking
  now share the same reduction and dims, and PCA is computed at most once.
