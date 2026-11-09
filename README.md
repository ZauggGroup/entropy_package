# Entropy

`Entropy` estimates cell-level transcriptional entropy from single-cell
expression data. Entropy is calculated across local cell neighborhoods and
measures how much neighbouring cells differ in their expression of each gene.
It is used to study transcriptional variability and changes in cell identity.
Interpreting a high score as plasticity is a claim about biology, not a
property of the measurement, so the package includes tools to test technical
alternatives before you make it.

[Package documentation](https://zaugggroup.github.io/entropy_package/)

## Overview

The package provides count-based entropy estimation using the
best upper bound (BUB) estimator and binary entropy estimation for
presence/absence data. Analyses can be performed directly from either Seurat
objects or AnnData `.h5ad` files.

Main capabilities include:

- automatic PCA-based nearest-neighbor construction and custom neighbor input;
- restriction of neighborhoods within biological or technical groups;
- count-based and binary entropy estimation;
- optional depth correction, which thins every cell to a common library depth
  before estimation because sequencing depth inflates entropy for a purely
  technical reason;
- scoring of selected feature sets;
- chunked and parallel computation for larger datasets;
- per-cell summaries and optional gene-by-cell entropy output.

## Choosing depth correction

A high entropy score can come from biology or from sequencing depth, gene
detection, neighborhood sparsity, or cell cycle. `entropy_depth_select()` tests
a grid of count-depth quantiles with one fixed neighbour graph and gene panel.
It compares removed depth-aligned variance with retained variation, includes
the uncorrected score as a candidate, and writes a compact PDF report.

```r
selection <- entropy_depth_select(
  object = seu,
  group = "celltype",
  outdir = "depth_selection",
  assay = "RNA",
  reduction = "pca_RNA"
)
selection
```

`entropy_recompute_check()` is interventional. It changes the input, by
equalizing count depth, equalizing the number of cells per cluster, or removing
a gene programme, then recomputes entropy and asks whether the cluster ranking
survives. That is stronger evidence than a statistical adjustment.

Neither approach proves plasticity. See
[Choosing depth correction](articles/choosing_depth_correction.html) for the
reasoning and two simulated frontiers.

## Installation

```r
devtools::install_github("ZauggGroup/entropy_package", force = TRUE)
```

## AnnData support

`run_entropy_h5ad()` provides direct support for AnnData `.h5ad` files without
requiring conversion to a Seurat object. It reads expression data, dimensional
reductions, and observation metadata while preserving AnnData cell and feature
identifiers. Results are returned as per-cell values that can be joined to
`adata.obs`.

Depth correction itself is available on the AnnData path. The automatic
depth-selection sweep is currently Seurat-only, so h5ad users should choose a
quantile in a Seurat analysis or run specified quantiles manually with
`run_entropy_h5ad()`.

Direct h5ad support requires R 4.5.0 or newer, `anndataR` 1.0.0 or newer, and
`rhdf5`. Installation, data preparation, analysis options, result validation,
and Scanpy visualization are described in the
[Python and AnnData workflow](articles/python_anndata_workflow.html).

## Documentation

- [Getting started](articles/getting_started.html)
- [Full Seurat workflow](articles/workflow.html)
- [Neighborhood definition](articles/neighborhood.html)
- [Choosing depth correction](articles/choosing_depth_correction.html)
- [Python and AnnData workflow](articles/python_anndata_workflow.html)
- [Function reference](reference/index.html)

## Support

For questions about the methodology, contact Judith B. Zaugg
(judith.zaugg@unibas.ch, judith.zaugg@embl.de) or Charles Dussiau
(charles.dussiau@embl.de).

For technical and implementation questions, contact Aryan Kamal
(aryan.kamal@unibas.ch).

Bug reports and feature requests can be submitted through the
[GitHub issue tracker](https://github.com/ZauggGroup/entropy_package/issues).

The package is also mirrored on the University of Basel sciCORE GitLab at
https://git.scicore.unibas.ch/zaugg-group/entropy_package, which is the institutional backup. Issues and releases live on GitHub.
