# Entropy: Estimate Cell-Level Entropy from Gene Expression Profiles

The *Entropy* package computes gene-expression entropy for individual cells
from either a Seurat object in R or an AnnData `.h5ad` file used by Python and
Scanpy. It provides a flexible and efficient way to measure transcriptional
variability and plasticity using cell neighborhoods.


### [Visit the Entropy package Webpage](https://grp-zaugg.embl-community.io/entropy_package/index.html)


---

## 🔍 Summary

Cellular entropy provides a quantitative measure of transcriptional uncertainty, useful in studies of differentiation, plasticity, or cellular identity. This package enables estimation of gene expression entropy across neighborhoods of cells, and adds the result back into the Seurat object as a new assay.

Key features:
- Uses BUB estimator with optional precomputed coefficients for speed
- Works directly on Seurat objects and AnnData `.h5ad` files
- Supports count-based and binary entropy
- Allows user-supplied neighbor lists or uses PCA-based kNN graph
- Restricts neighborhoods within cell types, samples, batches, or conditions
- Supports feature selection, chunking, and parallel processing
- Returns per-cell mean entropy and optionally saves the full gene-by-cell matrix

---

## 📦 Installation
```
library(devtools)
devtools::install_gitlab("grp-zaugg/entropy_package", host = "git.embl.de", force = TRUE)
```
---

## 🐍 Python and AnnData

Python users can calculate entropy directly from an `.h5ad` file without
converting it to a Seurat object. Direct h5ad support requires R 4.5.0 or newer
plus `anndataR` and `rhdf5`:

```r
install.packages("BiocManager")
BiocManager::install(c("anndataR", "rhdf5"))
```

Prepare the AnnData object in Python, preserving raw counts before
normalization:

```python
import scanpy as sc

adata = sc.read_h5ad("input.h5ad")
adata.layers["counts"] = adata.X.copy()

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.pca(adata, n_comps=50)
adata.write_h5ad("adata_for_entropy.h5ad")
```

Run entropy directly on the h5ad file in R:

```r
library(Entropy)

entropy_result <- run_entropy_h5ad(
  "adata_for_entropy.h5ad",
  layer = "counts",
  k = 15,
  mode = "counts"
)

head(entropy_result)
write.csv(entropy_result, "Entropy_mean.csv", quote = FALSE)
```

Join the result back into AnnData and visualize it:

```python
import pandas as pd

entropy = pd.read_csv(
    "Entropy_mean.csv",
    index_col=0,
)["Entropy_mean"]

adata.obs["Entropy_mean"] = entropy.reindex(adata.obs_names)

sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=15)
sc.tl.umap(adata, random_state=42)
sc.pl.umap(adata, color="Entropy_mean", color_map="viridis")
```

See the complete
[Python and AnnData workflow](articles/python_anndata_workflow.html) for binary
entropy, within-group neighborhoods, chunking, selected genes, output
validation, result comparisons, and UMAP examples.

---

### 📬 Get help and contact us

If you have questions or comments about the Entropy package, we’re happy to hear from you. We welcome inquiries related to both the scientific methodology and the technical implementation.

For method-related questions, please contact Judith B. Zaugg (judith.zaugg@embl.de) or Charles Dussiau (charles.dussiau@embl.de).

For technical or implementation-related questions, contact Aryan Kamal (aryan.kamal@embl.de).



If you have questions, doubts, ideas or problems, please use the [Gitlab Issue Tracker](https://git.embl.de/grp-zaugg/entropy_package/-/issues). We will respond in a timely manner.
