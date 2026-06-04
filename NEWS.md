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
  
  
# Entropy 0.2.2

### New features

* `run_entropy()` gains a `restrict_by` argument: pass the name of a column in
  `seu@meta.data` (e.g. `"cell_type"`) to restrict each cell's neighbors to
  cells sharing the same group, so entropy is computed within groups rather
  than across them. Works in both normal and chunked mode, and is `NULL`
  (off) by default.

  - For `neighbors_method = "auto"`, neighbors are recomputed *within* each
    group, so each cell gets up to `k` neighbors drawn only from its own group.
    Groups smaller than `k + 1` are capped to the group size, with a warning.
  - For `neighbors_method = "custom"` or `"seurat"`, the provided neighbor
    graph is filtered to keep only within-group neighbors (it is never
    recomputed, preserving whatever method produced it). Neighborhoods shrink
    as a result, especially at group boundaries; a before/after summary of
    mean neighbors per cell is logged.
  - Cells with `NA` in the `restrict_by` column belong to no group, receive no
    neighbors, and return entropy 0. The number of such cells is reported.
    
    
# Entropy 0.2.3

### Bug fixes

* `run_entropy()` now works with both Seurat v4 and v5 objects.
* When no reduction exists and the assay has no scaled data, `run_entropy()`
  now stops with a clear message instead of failing inside `RunPCA()`.
* Silenced verbose PCA output and fixed a reduction-key warning for assay
  names containing underscores (e.g. `pca_sct`).
  
  
# Entropy 0.2.4

### Bug fixes

* Fixed `run_entropy()` creating the wrong assay class on Seurat v4 objects.
  The entropy assay now matches the input object's version (`Assay` for v4,
  `Assay5` for v5) instead of being chosen from the installed package version.
* Silenced remaining neighbor-computation console output during automatic
  neighbor finding.

### Internal

* Added a test covering Seurat v4 and v5 compatibility for `run_entropy()`.




# Entropy 0.2.5

### Bug fixes

* `run_entropy()` now requires `k >= 11`, since the counts entropy estimator is
  undefined below 11 neighbors (previously this could crash).
* In `restrict_by` auto mode, groups with fewer than `k + 1` cells are now
  skipped (their cells return entropy 0) instead of being computed with a
  reduced, incomparable neighbor count.
* For custom and Seurat neighbor inputs, cells with fewer than 11 neighbors are
  now orphaned (entropy 0) with a summary warning, preventing the same crash.

### Internal

* Added tests covering basic runs, the `k >= 11` check, within-group
  restriction, small-group orphaning, and custom-neighbor orphaning.
  
  
  
  
# Entropy 0.3.0

### New features

* New function `add_entropy_score()`: given a precomputed entropy assay and a
  named list of gene sets, it adds one metadata column per set
  (`entropy_<set name>`) holding the mean entropy of that set's genes. A fast
  way to score pathways or regulons without recomputing entropy. Missing genes
  raise an error unless `force_drop = TRUE`; a set with no present genes always
  errors.
* `run_entropy()` gains a `meta_name` argument: when set, the per-cell mean
  entropy is written to that metadata column regardless of `add_assay`. If
  either `add_assay = TRUE` or `meta_name` is set, the Seurat object is
  returned; otherwise a data frame of means is returned.
