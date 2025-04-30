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
