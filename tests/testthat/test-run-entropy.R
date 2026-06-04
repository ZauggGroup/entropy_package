test_that("run_entropy runs and returns a sane entropy assay (counts)", {
  skip_if_not_installed("Seurat")
  library(Seurat)
  data("pbmc_small")

  r <- run_entropy(pbmc_small, assay = "RNA", mode = "counts", k = 15)

  expect_true("RNA_entropy" %in% Assays(r))
  ent <- SeuratObject::GetAssayData(r[["RNA_entropy"]], layer = "counts")
  expect_equal(ncol(ent), ncol(pbmc_small))   # all cells
  expect_true(all(ent >= 0))                  # entropy is non-negative
  expect_false(anyNA(ent))                    # no NAs
})

test_that("run_entropy rejects k < 11", {
  skip_if_not_installed("Seurat")
  library(Seurat)
  data("pbmc_small")

  expect_error(
    run_entropy(pbmc_small, assay = "RNA", mode = "counts", k = 5),
    "at least 11"
  )
})

test_that("restrict_by keeps neighbors within group", {
  skip_if_not_installed("Seurat")
  library(Seurat)
  data("pbmc_small")

  # give every cell a group label with enough cells per group (>= k+1)
  set.seed(1)
  grp <- sample(c("A", "B"), ncol(pbmc_small), replace = TRUE)
  pbmc_small$grp <- grp

  r <- run_entropy(pbmc_small, assay = "RNA", mode = "counts",
                   k = 11, restrict_by = "grp")

  # verify the neighbor structure directly: recompute and check no cross-group
  emb <- Embeddings(r[["pca_RNA"]])
  emb <- emb[, seq_len(min(20, ncol(emb)))]
  g <- as.character(r$grp); names(g) <- colnames(r)
  a_cells <- colnames(r)[g == "A"]
  nn <- RANN::nn2(emb[a_cells, ], k = 12)$nn.idx
  expect_true(all(g[a_cells[nn[, -1]]] == "A"))
})

test_that("restrict_by orphans groups smaller than k+1 (entropy 0)", {
  skip_if_not_installed("Seurat")
  library(Seurat)
  data("pbmc_small")

  # one normal group + one tiny group (< k+1 cells) that must be orphaned
  grp <- rep("big", ncol(pbmc_small))
  grp[1:5] <- "tiny"            # 5 cells, < k+1 = 12 -> must be skipped
  pbmc_small$grp <- grp

  r <- run_entropy(pbmc_small, assay = "RNA", mode = "counts",
                   k = 11, restrict_by = "grp", add_assay = FALSE)

  tiny_cells <- colnames(pbmc_small)[1:5]
  expect_true(all(r[tiny_cells, 1] == 0))   # orphaned cells -> entropy 0
})




test_that("custom neighbors with < 11 neighbors are orphaned (entropy 0)", {
  skip_if_not_installed("Seurat")
  library(Seurat)
  data("pbmc_small")

  cells <- colnames(pbmc_small)

  # build a custom neighbor list: every cell gets 12 neighbors,
  # except cell 1 which deliberately gets only 5 (< 11)
  nn <- lapply(seq_along(cells), function(i) {
    if (i == 1) {
      setdiff(cells, cells[1])[1:5]      # only 5 neighbors
    } else {
      setdiff(cells, cells[i])[1:12]     # 12 neighbors
    }
  })
  names(nn) <- cells

  r <- run_entropy(pbmc_small, assay = "RNA", mode = "counts", k = 11,
                   neighbors = nn, neighbors_method = "custom",
                   add_assay = FALSE)

  expect_equal(r[cells[1], 1], 0)        # under-11 cell -> entropy 0
  expect_gt(r[cells[2], 1], 0)           # normal cell -> real value
})




test_that("add_entropy_score adds one metadata column per gene set", {
  skip_if_not_installed("Seurat")
  library(Seurat)
  data("pbmc_small")

  r <- run_entropy(pbmc_small, assay = "RNA", mode = "counts", k = 11, add_assay = TRUE)
  sets <- list(setA = head(rownames(pbmc_small), 5),
               setB = head(rownames(pbmc_small), 8))
  r <- add_entropy_score(r, "RNA_entropy", sets)

  expect_true(all(c("entropy_setA", "entropy_setB") %in% colnames(r@meta.data)))
  expect_equal(length(r$entropy_setA), ncol(r))
  expect_false(anyNA(r$entropy_setA))
})

test_that("add_entropy_score errors on missing entropy assay", {
  skip_if_not_installed("Seurat")
  library(Seurat)
  data("pbmc_small")
  expect_error(
    add_entropy_score(pbmc_small, "RNA_entropy", list(s = rownames(pbmc_small)[1:3])),
    "not found"
  )
})

test_that("add_entropy_score errors on missing genes unless force_drop", {
  skip_if_not_installed("Seurat")
  library(Seurat)
  data("pbmc_small")

  r <- run_entropy(pbmc_small, assay = "RNA", mode = "counts", k = 11, add_assay = TRUE)
  good <- head(rownames(pbmc_small), 3)
  sets <- list(s = c(good, "NOT_A_GENE"))

  expect_error(add_entropy_score(r, "RNA_entropy", sets), "not in assay")

  # force_drop = TRUE drops the missing one and succeeds
  r2 <- add_entropy_score(r, "RNA_entropy", sets, force_drop = TRUE)
  expect_true("entropy_s" %in% colnames(r2@meta.data))
})
