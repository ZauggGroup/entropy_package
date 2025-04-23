# Entropy: Estimate Cell-Level Entropy from Gene Expression Profiles

The *Entropy* package allows you to compute gene-level entropy for single cells based on a Seurat object. It provides a flexible and efficient way to measure transcriptional variability and plasticity, with support for neighbor-aware entropy computation and integration into Seurat workflows.

---

## 🔍 Summary

Cellular entropy provides a quantitative measure of transcriptional uncertainty, useful in studies of differentiation, plasticity, or cellular identity. This package enables estimation of gene expression entropy across neighborhoods of cells, and adds the result back into the Seurat object as a new assay.

Key features:
- Uses BUB estimator with optional precomputed coefficients for speed
- Works directly on any Seurat object
- Allows user-supplied neighbor lists or uses PCA-based kNN graph
- Returns either a summary matrix or adds a full assay to the Seurat object

---

## 📦 Installation
```
library(devtools)
devtools::install_gitlab("grp-zaugg/entropy_package", host = "git.embl.de", force = TRUE)
```
---

### 📬 Get help and contact us

If you have questions or comments about the Entropy package, we’re happy to hear from you. We welcome inquiries related to both the scientific methodology and the technical implementation.

For method-related questions, please contact Judith B. Zaugg (judith.zaugg@embl.de) or Charles Dussiau (charles.dussiau@embl.de).

For technical or implementation-related questions, contact Aryan Kamal (aryan.kamal@embl.de).



If you have questions, doubts, ideas or problems, please use the [Gitlab Issue Tracker](https://git.embl.de/grp-zaugg/entropy_package/-/issues). We will respond in a timely manner.

