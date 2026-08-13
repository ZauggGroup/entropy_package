# Entropy 0.6.0

* `run_entropy()` gains optional count-depth correction. Sequencing depth can
  inflate entropy because deeper libraries sample the expression distribution
  more completely. With `depth_correct = TRUE`, the function also computes
  entropy after thinning raw counts to a common depth. The default remains
  `FALSE`.

* Every `run_entropy()` call that returns a Seurat object now records its
  parameters in `seu@misc$entropy_provenance`. This provenance record is the
  only visible change when `depth_correct = FALSE`.

* Removed `entropy_density_diagnostic()` as a breaking change. Use
  `entropy_diagnose()` for computation and `entropy_diagnose_report()` for PDF
  reporting.

# Entropy 0.5.0

## Interventional entropy checks

* Added `entropy_recompute_check()` and `entropy_recompute_report()`.

* `entropy_diagnose()` remains an observational screen. It compares cells at
  matched covariate values, but it cannot distinguish a technical confounder
  from a biological mediator.

* The new recompute check is interventional. It can equalise count depth,
  equalise the number of cells sampled per cluster, or remove a gene programme,
  then call `run_entropy()` again. Its headline result asks whether the cluster
  ranking survives that physical change.

* Random interventions report repeat-to-repeat ranges and label inconsistent
  cluster movement as unstable. Gene removal reports ranking as the comparable
  result because removing features shifts the absolute entropy scale.

# Entropy 0.4.0

## New Class A entropy diagnostic

* Replaced the residual-based density diagnostic with `entropy_diagnose()` and
  `entropy_diagnose_report()`. Computation and PDF rendering are now separate.

* The old diagnostic had no overlap check, showed no raw entropy reference,
  used zero-sum regression residuals whose baseline depended on dataset
  composition, and exposed a log-density option that could change cluster
  conclusions.

* The new diagnostic compares each cluster with other cells in matched
  covariate bins. It reports the matched effect size, comparison support,
  within-cluster association, global association, raw cluster values, and an
  optional joint model.

* Matched bins must also pass a within-bin median balance check. This prevents
  wide quantile bins in a covariate tail from being treated as fair matches.

* Detailed reports show dropped bins in a separate balance strip instead of
  placing them on the zero-difference line. Loess fitting now uses a
  reproducible random subsample without changing the caller's random stream.

* Clusters without a fair comparison group now receive `cannot tell` instead
  of a confident result near zero.

* Added one-page summary PDFs and detailed five-part covariate reports using
  the base PDF device, ggplot2, and patchwork.

* `entropy_density_diagnostic()` remains as a deprecated wrapper. It now
  returns an `entropy_diagnostic` object rather than the old list shape.

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
  
  
  
# Entropy 0.3.1

* Fixed a bug where using `meta_name` together with `add_assay = TRUE` could error and leave the object unmodified when `meta_name` collided with an existing assay name. The metadata column is now automatically renamed with a `_meta` suffix and a warning is emitted instead.


# Entropy 0.3.5

## New functionality

* Added `run_entropy_h5ad()` for entropy estimation directly from AnnData
  `.h5ad` files, without conversion to a Seurat object.

* Added support for expression data in `adata.X` or a named AnnData layer,
  embeddings in `adata.obsm`, and grouping variables in `adata.obs`. Cell and
  feature identifiers are preserved in the returned results.

* The h5ad interface supports count-based and binary entropy, automatic or
  custom neighbors, within-group neighborhood restriction, selected features,
  chunked computation, optional parallelization, and saving the full
  gene-by-cell entropy matrix.

* Per-cell h5ad results are returned in the `Entropy_mean` column for direct
  alignment with `adata.obs`.

## Documentation and testing

* Added a Python and AnnData vignette covering Scanpy data preparation,
  count-based and binary analyses, within-group neighborhoods, chunking,
  selected features, result transfer, and visualization.

* Added tests using a synthetic AnnData object for automatic neighbors,
  count-based and binary entropy, within-cell-type restriction, and agreement
  between chunked and unchunked computation.

* Added validation for h5ad layers, embeddings, dimensions, feature names,
  count matrices, neighborhood size, and optional dependency versions.

## Compatibility

* Direct h5ad support requires R 4.5.0 or newer, `anndataR` 1.0.0 or newer,
  and `rhdf5`. These remain optional dependencies; the existing Seurat
  interface is unchanged.


# Entropy 0.3.6

## New density and covariate diagnostic

* Added `entropy_density_diagnostic()` to assess whether a per-cell
  entropy/plasticity signal is associated with local neighbourhood sparsity.
  The function measures the mean distance to each cell's `k` nearest
  neighbours, excluding the cell itself, and reports the global Spearman
  entropy-density correlation.

* Density is always tested in its own independent model
  (`entropy ~ mean_knn_dist`). Each optional numeric covariate is also tested
  independently (`entropy ~ covariate`); density and additional covariates are
  never combined in the same regression.

* Added per-cell diagnostic residuals, regression summaries, per-group mean
  and median residual summaries, one boxplot per independent model, optional
  group highlighting, log-transformed density, and a density/covariate
  Spearman correlation matrix.

* Plot saving is opt-in. When requested, plots are written as 300-dpi PNGs,
  with one file per independent model.

* Residuals are explicitly documented as diagnostics rather than corrected
  entropy scores. A sparse state may be genuinely plastic, so the diagnostic
  characterizes possible confounding but does not establish causality or
  remove biological signal.

## Documentation

* Added a self-contained vignette explaining the density confound, independent
  covariate checks, result interpretation, limitations, plotting, and saving.
* Added the new diagnostic to the pkgdown article and function-reference
  indexes.
