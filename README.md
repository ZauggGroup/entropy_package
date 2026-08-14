# Entropy

`Entropy` estimates cell-level transcriptional entropy from single-cell
expression data. Entropy is calculated across local cell neighborhoods and
measures how much neighbouring cells differ in their expression of each gene.
It is used to study transcriptional variability and changes in cell identity.
Interpreting a high score as plasticity is a claim about biology, not a
property of the measurement, so the package ships diagnostics to test that
interpretation before you make it.

[Package documentation](https://grp-zaugg.embl-community.io/entropy_package/)

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

## Diagnostics

A high entropy score can come from biology or from sequencing depth, gene
detection, neighborhood sparsity, or cell cycle. Two functions test which.

`entropy_diagnose()` is observational. For every cluster and covariate it asks
whether the cluster still has higher entropy once cells are compared at matched
covariate values, and reports how much of the elevation survives that matching.
`entropy_diagnose_report()` writes the result as PDF reports.

`entropy_recompute_check()` is interventional. It changes the input, by
equalizing count depth, equalizing the number of cells per cluster, or removing
a gene programme, then recomputes entropy and asks whether the cluster ranking
survives. That is stronger evidence than a statistical adjustment.

Neither proves plasticity. They rule out alternative explanations one at a
time. See
[Diagnosing confounds](articles/diagnosing_confounds.html) for the reasoning
and a worked example on simulated data with known causes.

## Installation

```r
devtools::install_gitlab(
  "grp-zaugg/entropy_package",
  host = "git.embl.de",
  force = TRUE
)
```

## AnnData support

`run_entropy_h5ad()` provides direct support for AnnData `.h5ad` files without
requiring conversion to a Seurat object. It reads expression data, dimensional
reductions, and observation metadata while preserving AnnData cell and feature
identifiers. Results are returned as per-cell values that can be joined to
`adata.obs`.

The diagnostics work from AnnData too. `run_entropy_h5ad(return_inputs = TRUE)`
returns the entropy vector, the embedding, and the `obs` metadata in one call,
and `entropy_diagnose()` accepts those directly with `seu = NULL`, so no Seurat
object is needed anywhere. Depth correction is available on this path as well.

Direct h5ad support requires R 4.5.0 or newer, `anndataR` 1.0.0 or newer, and
`rhdf5`. Installation, data preparation, analysis options, result validation,
and Scanpy visualization are described in the
[Python and AnnData workflow](articles/python_anndata_workflow.html).

## Documentation

- [Getting started](articles/getting_started.html)
- [Full Seurat workflow](articles/workflow.html)
- [Neighborhood definition](articles/neighborhood.html)
- [Diagnosing confounds](articles/diagnosing_confounds.html)
- [Python and AnnData workflow](articles/python_anndata_workflow.html)
- [Function reference](reference/index.html)

## Support

For questions about the methodology, contact Judith B. Zaugg
(judith.zaugg@embl.de) or Charles Dussiau
(charles.dussiau@embl.de).

For technical and implementation questions, contact Aryan Kamal
(aryan.kamal@embl.de).

Bug reports and feature requests can be submitted through the
[GitLab issue tracker](https://git.embl.de/grp-zaugg/entropy_package/-/issues).
