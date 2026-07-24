# Entropy

`Entropy` estimates cell-level transcriptional entropy from single-cell
expression data. Entropy is calculated across local cell neighborhoods and can
be used to investigate transcriptional variability, cellular plasticity, and
changes in cell identity.

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
- scoring of selected feature sets;
- chunked and parallel computation for larger datasets;
- per-cell summaries and optional gene-by-cell entropy output.

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

Direct h5ad support requires R 4.5.0 or newer, `anndataR` 1.0.0 or newer, and
`rhdf5`. Installation, data preparation, analysis options, result validation,
and Scanpy visualization are described in the
[Python and AnnData workflow](articles/python_anndata_workflow.html).

## Documentation

- [Getting started](articles/getting_started.html)
- [Full Seurat workflow](articles/workflow.html)
- [Neighborhood definition](articles/neighborhood.html)
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
