test_that("run_entropy returns an assay matching the object's Seurat version", {

  skip_if_not_installed("Seurat")

  library(Seurat)
  data("pbmc_small")
  counts <- SeuratObject::GetAssayData(pbmc_small, layer = "counts")

  # build a v4-style object
  options(Seurat.object.assay.version = "v3")
  seu_v4 <- CreateSeuratObject(counts = counts)
  seu_v4 <- NormalizeData(seu_v4, verbose = FALSE)
  seu_v4 <- FindVariableFeatures(seu_v4, verbose = FALSE)
  seu_v4 <- ScaleData(seu_v4, verbose = FALSE)

  # build a v5-style object
  options(Seurat.object.assay.version = "v5")
  seu_v5 <- CreateSeuratObject(counts = counts)
  seu_v5 <- NormalizeData(seu_v5, verbose = FALSE)
  seu_v5 <- FindVariableFeatures(seu_v5, verbose = FALSE)
  seu_v5 <- ScaleData(seu_v5, verbose = FALSE)

  # sanity: the two objects really are different versions
  expect_true(inherits(seu_v4[["RNA"]], "Assay"))
  expect_true(inherits(seu_v5[["RNA"]], "Assay5"))

  feats <- head(rownames(seu_v4), 2)

  # run on both; should not error
  r4 <- run_entropy(seu_v4, "RNA", features = feats, add_assay = TRUE, mode = "counts")
  r5 <- run_entropy(seu_v5, "RNA", features = feats, add_assay = TRUE, mode = "counts")

  # the entropy assay class must match the input object's version
  expect_true(inherits(r4[["RNA_entropy"]], "Assay"))
  expect_true(inherits(r5[["RNA_entropy"]], "Assay5"))

  # and the matrix should be the right shape: 2 genes x all cells
  expect_equal(nrow(r4[["RNA_entropy"]]), 2)
  expect_equal(ncol(r4[["RNA_entropy"]]), ncol(seu_v4))
})
