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
