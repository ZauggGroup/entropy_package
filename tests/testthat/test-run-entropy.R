test_that("depth correction disabled preserves released run_entropy behaviour", {
  skip_if_not_installed("Seurat")
  library(Seurat)
  data("pbmc_small")

  # Generated from git show HEAD:R/entropy.R before depth correction was added.
  # Regenerate it with that version and the parameters stored in the RDS only
  # when the entropy estimator is deliberately changed.
  reference <- readRDS(testthat::test_path(
    "fixtures", "run_entropy_reference.rds"
  ))
  cells <- colnames(pbmc_small)
  neighbors <- lapply(seq_along(cells), function(index) {
    setdiff(cells, cells[index])[seq_len(12L)]
  })
  names(neighbors) <- cells

  comparison <- run_entropy(
    pbmc_small,
    assay = reference$parameters$assay,
    mode = reference$parameters$mode,
    # k must match the list length now. The reference was generated with a
    # 12-neighbour list, and k was unused on the custom path back then, so
    # this reproduces the same numbers.
    k = length(neighbors[[1L]]),
    neighbors = neighbors,
    neighbors_method = reference$parameters$neighbors_method,
    add_assay = reference$parameters$add_assay,
    meta_name = reference$parameters$meta_name,
    depth_correct = FALSE
  )

  expect_identical(comparison$released_entropy, reference$means)
  expect_identical(
    colnames(comparison@meta.data), reference$metadata_columns
  )
  expect_false(any(grepl("depth_corrected", colnames(comparison@meta.data))))
})


make_depth_correction_fixture <- function() {
  n_cells <- 24L
  totals <- seq(240L, 470L, length.out = n_cells)
  totals <- as.integer(round(totals))
  counts <- vapply(totals, function(total) {
    c(total - 20L, 8L, 5L, 3L, 2L, 2L)
  }, integer(6L))
  dimnames(counts) <- list(
    paste0("gene", seq_len(nrow(counts))),
    paste0("cell", seq_len(ncol(counts)))
  )
  seu <- Seurat::CreateSeuratObject(counts = Matrix::Matrix(counts, sparse = TRUE))
  cells <- colnames(seu)
  neighbors <- lapply(seq_along(cells), function(index) {
    setdiff(cells, cells[index])[seq_len(12L)]
  })
  names(neighbors) <- cells
  list(seu = seu, counts = counts, totals = totals, neighbors = neighbors)
}


run_depth_correction_fixture <- function(fixture, ...) {
  arguments <- list(
    seu = fixture$seu,
    assay = "RNA",
    neighbors = fixture$neighbors,
    neighbors_method = "custom",
    k = 12,
    mode = "counts",
    add_assay = FALSE,
    depth_correct = TRUE,
    target_depth = 240L
  )
  arguments[names(list(...))] <- list(...)
  do.call(run_entropy, arguments)
}


test_that("depth correction adds expected metadata and keeps original means", {
  skip_if_not_installed("Seurat")
  fixture <- make_depth_correction_fixture()
  original <- run_entropy(
    fixture$seu, assay = "RNA", neighbors = fixture$neighbors,
    neighbors_method = "custom", k = 12, mode = "counts",
    add_assay = FALSE, meta_name = "mean_entropy_RNA"
  )
  set.seed(10)
  corrected <- run_depth_correction_fixture(fixture)

  added <- setdiff(colnames(corrected@meta.data), colnames(fixture$seu@meta.data))
  expect_setequal(added, c(
    "entropy_mean_RNA", "entropy_mean_RNA_depth_corrected",
    "depth_corrected_RNA"
  ))
  expect_identical(corrected$entropy_mean_RNA, original$mean_entropy_RNA)
  expect_type(corrected$depth_corrected_RNA, "logical")
})


test_that("downsampling caps deep cells and leaves shallow cells unchanged", {
  fixture <- make_depth_correction_fixture()
  target <- 300L
  set.seed(20)
  downsampled <- Entropy:::.entropy_downsample_matrix(
    Matrix::Matrix(fixture$counts, sparse = TRUE), target
  )
  after <- Matrix::colSums(downsampled)
  below <- fixture$totals <= target

  expect_true(all(after[!below] <= target))
  expect_equal(
    as.matrix(downsampled[, below, drop = FALSE]),
    fixture$counts[, below, drop = FALSE]
  )
})


test_that("depth correction records target counts and target override", {
  skip_if_not_installed("Seurat")
  fixture <- make_depth_correction_fixture()
  set.seed(30)
  result <- suppressWarnings(run_depth_correction_fixture(
    fixture, target_depth = 300L, depth_quantile = 0.9
  ))

  expect_identical(result@misc$entropy_provenance$target_depth, 300L)
  expect_identical(
    result@misc$entropy_provenance$n_cells_below_target,
    sum(fixture$totals < 300L)
  )
  expect_identical(
    result@misc$entropy_provenance$n_cells_corrected,
    sum(fixture$totals > 300L)
  )
})


test_that("depth correction validates counts and target arguments", {
  skip_if_not_installed("Seurat")
  fixture <- make_depth_correction_fixture()

  for (bad_quantile in list(0, 1, c(0.1, 0.2))) {
    expect_error(
      run_depth_correction_fixture(
        fixture, depth_quantile = bad_quantile, target_depth = NULL
      ),
      "strictly between 0 and 1"
    )
  }
  expect_error(
    run_depth_correction_fixture(fixture, target_depth = -1),
    "positive integer"
  )

  no_counts <- fixture$seu
  suppressWarnings(
    no_counts[["data_only"]] <- SeuratObject::CreateAssay5Object(
      data = log1p(Matrix::Matrix(fixture$counts, sparse = TRUE))
    )
  )
  expect_error(
    suppressWarnings(run_entropy(
      no_counts, assay = "data_only", neighbors = fixture$neighbors,
      neighbors_method = "custom", k = 11, add_assay = FALSE,
      depth_correct = TRUE, target_depth = 240L
    )),
    "requires raw counts"
  )

  non_integer <- fixture$seu
  non_integer_counts <- SeuratObject::GetAssayData(
    non_integer, assay = "RNA", layer = "counts"
  )
  non_integer_counts[1, 1] <- non_integer_counts[1, 1] + 0.5
  non_integer[["non_integer"]] <- SeuratObject::CreateAssay5Object(
    counts = non_integer_counts
  )
  expect_error(
    run_entropy(
      non_integer, assay = "non_integer", neighbors = fixture$neighbors,
      neighbors_method = "custom", k = 11, add_assay = FALSE,
      depth_correct = TRUE, target_depth = 240L
    ),
    "not integer-like"
  )
})


test_that("depth correction warns for shallow cells and binary mode", {
  skip_if_not_installed("Seurat")
  fixture <- make_depth_correction_fixture()

  expect_warning(
    run_depth_correction_fixture(fixture, target_depth = 300L),
    "below the target depth"
  )
  expect_warning(
    run_depth_correction_fixture(fixture, target_depth = 199L),
    "may be unreliable"
  )
  set.seed(40)
  binary <- NULL
  expect_warning(
    binary <- run_depth_correction_fixture(
      fixture, mode = "binary", target_depth = 240L
    ),
    "not been validated for binary mode"
  )
  expect_s4_class(binary, "Seurat")
})


test_that("depth corrected assay stores corrected entropy only", {
  skip_if_not_installed("Seurat")
  fixture <- make_depth_correction_fixture()
  set.seed(50)
  result <- NULL
  expect_message(
    result <- suppressWarnings(run_depth_correction_fixture(
      fixture, add_assay = TRUE, target_depth = 240L
    )),
    "DEPTH CORRECTED matrix only"
  )
  stored <- SeuratObject::GetAssayData(
    result, assay = "RNA_entropy", layer = "counts"
  )

  expect_identical(
    unname(Matrix::colMeans(stored)),
    unname(result$entropy_mean_RNA_depth_corrected)
  )
  expect_false(identical(
    unname(Matrix::colMeans(stored)), unname(result$entropy_mean_RNA)
  ))
})


test_that("metadata collision keeps corrected column names paired", {
  skip_if_not_installed("Seurat")
  fixture <- make_depth_correction_fixture()
  set.seed(55)
  result <- suppressWarnings(run_depth_correction_fixture(
    fixture, meta_name = "RNA"
  ))

  expect_true(all(c(
    "RNA_meta", "RNA_meta_depth_corrected"
  ) %in% colnames(result@meta.data)))
  expect_false("RNA_depth_corrected" %in% colnames(result@meta.data))
})


test_that("run_entropy writes provenance on corrected and uncorrected paths", {
  skip_if_not_installed("Seurat")
  fixture <- make_depth_correction_fixture()
  uncorrected <- run_entropy(
    fixture$seu, assay = "RNA", neighbors = fixture$neighbors,
    neighbors_method = "custom", k = 12, mode = "counts",
    add_assay = FALSE, meta_name = "entropy"
  )
  set.seed(60)
  corrected <- run_depth_correction_fixture(fixture)
  expected_names <- c(
    "k", "assay", "mode", "reduction", "dims", "restrict_by",
    "depth_correct", "depth_quantile", "target_depth",
    "n_cells_below_target", "n_cells_corrected"
  )

  expect_identical(names(uncorrected@misc$entropy_provenance), expected_names)
  expect_false(uncorrected@misc$entropy_provenance$depth_correct)
  expect_identical(names(corrected@misc$entropy_provenance), expected_names)
  expect_true(corrected@misc$entropy_provenance$depth_correct)
})


test_that("depth target uses the full library before feature subsetting", {
  skip_if_not_installed("Seurat")
  fixture <- make_depth_correction_fixture()
  expected <- floor(stats::quantile(
    fixture$totals, probs = 0.5, names = FALSE, type = 7
  ))
  set.seed(70)
  result <- suppressWarnings(run_depth_correction_fixture(
    fixture, target_depth = NULL, depth_quantile = 0.5,
    features = rownames(fixture$seu)[1:2]
  ))

  expect_identical(result@misc$entropy_provenance$target_depth, expected)
})


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
  expect_error(
    run_entropy(pbmc_small, assay = "RNA", mode = "counts", k = 11.5),
    "single integer"
  )
  expect_error(
    run_entropy(pbmc_small, assay = "RNA", mode = "counts", k = NA_real_),
    "single integer"
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




test_that("custom neighbors of unequal size are rejected", {
  skip_if_not_installed("Seurat")
  library(Seurat)
  data("pbmc_small")

  cells <- colnames(pbmc_small)
  nn <- lapply(seq_along(cells), function(i) {
    if (i == 1) setdiff(cells, cells[1])[1:5] else setdiff(cells, cells[i])[1:12]
  })
  names(nn) <- cells

  expect_error(
    run_entropy(pbmc_small, assay = "RNA", mode = "counts", k = 12,
                neighbors = nn, add_assay = FALSE),
    "different sizes"
  )
})


test_that("custom neighbors containing the cell itself are rejected", {
  skip_if_not_installed("Seurat")
  library(Seurat)
  data("pbmc_small")

  cells <- colnames(pbmc_small)
  nn <- lapply(seq_along(cells), function(i) c(cells[i], setdiff(cells, cells[i])[1:11]))
  names(nn) <- cells

  expect_error(
    run_entropy(pbmc_small, assay = "RNA", mode = "counts", k = 12,
                neighbors = nn, add_assay = FALSE),
    "their own neighbour"
  )
})


test_that("custom neighbors reject repeated and non-integer positions", {
  skip_if_not_installed("Seurat")
  library(Seurat)
  data("pbmc_small")

  cells <- colnames(pbmc_small)
  nn <- setNames(lapply(seq_along(cells), function(i) {
    setdiff(seq_along(cells), i)[seq_len(12)]
  }), cells)

  repeated <- nn
  repeated[[1L]][2L] <- repeated[[1L]][1L]
  expect_error(
    run_entropy(
      pbmc_small, assay = "RNA", mode = "counts", k = 12,
      neighbors = repeated, add_assay = FALSE
    ),
    "repeat a neighbour"
  )

  fractional <- nn
  fractional[[1L]][1L] <- 2.5
  expect_error(
    run_entropy(
      pbmc_small, assay = "RNA", mode = "counts", k = 12,
      neighbors = fractional, add_assay = FALSE
    ),
    "finite integers"
  )
})


test_that("k is read from the neighbour list when it is not supplied", {
  skip_if_not_installed("Seurat")
  library(Seurat)
  data("pbmc_small")

  cells <- colnames(pbmc_small)
  nn <- lapply(seq_along(cells), function(i) setdiff(cells, cells[i])[1:12])
  names(nn) <- cells

  r <- run_entropy(pbmc_small, assay = "RNA", mode = "counts",
                   neighbors = nn, add_assay = FALSE, meta_name = "e")
  expect_identical(r@misc$entropy_provenance$k, 12L)
})


test_that("k that disagrees with the neighbour list is rejected", {
  skip_if_not_installed("Seurat")
  library(Seurat)
  data("pbmc_small")

  cells <- colnames(pbmc_small)
  nn <- lapply(seq_along(cells), function(i) setdiff(cells, cells[i])[1:12])
  names(nn) <- cells

  expect_error(
    run_entropy(pbmc_small, assay = "RNA", mode = "counts", k = 15,
                neighbors = nn, add_assay = FALSE),
    "was requested but every cell has"
  )
})


test_that("an explicitly all-empty neighbour list is allowed", {
  skip_if_not_installed("Seurat")
  library(Seurat)
  data("pbmc_small")

  nn <- setNames(rep(list(integer(0)), ncol(pbmc_small)), colnames(pbmc_small))
  result <- run_entropy(
    pbmc_small, assay = "RNA", mode = "counts", k = 11,
    neighbors = nn, add_assay = FALSE
  )

  expect_true(all(result[[1L]] == 0))
})


test_that("restrict_by cannot create unequal supplied neighbourhoods", {
  skip_if_not_installed("Seurat")
  library(Seurat)
  data("pbmc_small")

  cells <- colnames(pbmc_small)
  pbmc_small$grp <- c("B", rep("A", length(cells) - 1L))
  nn <- setNames(lapply(seq_along(cells), function(i) {
    setdiff(cells, cells[[i]])[seq_len(12)]
  }), cells)

  expect_error(
    run_entropy(
      pbmc_small, assay = "RNA", mode = "counts", k = 12,
      neighbors = nn, restrict_by = "grp", add_assay = FALSE
    ),
    "Drop `neighbors`"
  )
})


test_that("auto neighbours exclude the cell itself and are cached", {
  skip_if_not_installed("Seurat")
  library(Seurat)
  data("pbmc_small")

  r <- run_entropy(pbmc_small, assay = "RNA", mode = "counts", k = 11)
  nn <- r@misc$entropy_neighbors

  expect_s3_class(nn, "entropy_neighbors")
  expect_identical(names(nn), colnames(pbmc_small))
  expect_true(all(lengths(nn) == 11L))
  expect_false(any(vapply(seq_along(nn), function(i) i %in% nn[[i]], logical(1))))

  # second call on the returned object reuses it
  r2 <- run_entropy(r, assay = "RNA", mode = "counts", k = 11, add_assay = FALSE,
                    meta_name = "again")
  expect_identical(
    attr(r2@misc$entropy_neighbors, "built"), attr(nn, "built")
  )
})


test_that("entropy_neighbors builds a valid non-self graph", {
  skip_if_not_installed("Seurat")
  library(Seurat)
  data("pbmc_small")

  pbmc_small <- suppressWarnings(RunPCA(pbmc_small, verbose = FALSE))
  nn <- entropy_neighbors(pbmc_small, reduction = "pca", dims = 10, k = 15)

  expect_s3_class(nn, "entropy_neighbors")
  expect_true(all(lengths(nn) == 15L))
  expect_false(any(vapply(seq_along(nn), function(i) i %in% nn[[i]], logical(1))))

  plain <- .entropy_nn_normalise(nn, colnames(pbmc_small), k = 15)
  expect_identical(attributes(plain), list(names = colnames(pbmc_small)))
})


test_that("self removal does not assume the first neighbour column", {
  nn_idx <- rbind(
    c(2L, 1L, 3L),
    c(2L, 1L, 3L),
    c(1L, 3L, 2L)
  )

  expect_identical(
    .entropy_nn_drop_self(nn_idx),
    list(c(2L, 3L), c(1L, 3L), c(1L, 2L))
  )
})


test_that("computed candidates are trimmed to k after self removal", {
  nn_idx <- rbind(
    c(2L, 3L, 4L),  # self absent: all three survive identity filtering
    c(2L, 1L, 3L),
    c(1L, 3L, 2L),
    c(4L, 2L, 1L)
  )

  trimmed <- .entropy_nn_trim(.entropy_nn_drop_self(nn_idx), k = 2L)
  expect_true(all(lengths(trimmed) == 2L))
  expect_identical(trimmed[[1L]], c(2L, 3L))
})


test_that("tie-heavy computed graphs contain exactly k neighbours", {
  skip_if_not_installed("Seurat")
  library(Seurat)
  data("pbmc_small")

  cells <- colnames(pbmc_small)
  tied <- matrix(
    0, nrow = length(cells), ncol = 2,
    dimnames = list(cells, c("TIE_1", "TIE_2"))
  )
  pbmc_small[["ties"]] <- SeuratObject::CreateDimReducObject(
    embeddings = tied, key = "TIE_", assay = "RNA"
  )

  fixed <- entropy_neighbors(
    pbmc_small, reduction = "ties", dims = 2, k = 11
  )
  expect_true(all(lengths(fixed) == 11L))

  global <- run_entropy(
    pbmc_small, assay = "RNA", reduction = "ties", dims = 2, k = 11,
    mode = "binary", add_assay = FALSE, meta_name = "e"
  )
  expect_true(all(lengths(global@misc$entropy_neighbors) == 11L))

  pbmc_small$grp <- rep(c("A", "B"), each = length(cells) / 2L)
  restricted <- run_entropy(
    pbmc_small, assay = "RNA", reduction = "ties", dims = 2, k = 11,
    mode = "binary", restrict_by = "grp",
    add_assay = FALSE, meta_name = "e"
  )
  expect_true(all(lengths(restricted@misc$entropy_neighbors) == 11L))
})


test_that("automatic builders reject fewer than k plus one cells", {
  skip_if_not_installed("Seurat")
  library(Seurat)
  data("pbmc_small")

  small <- pbmc_small[, seq_len(14L)]
  cells <- colnames(small)
  embedding <- matrix(
    0, nrow = length(cells), ncol = 2,
    dimnames = list(cells, c("SMALL_1", "SMALL_2"))
  )
  small[["small"]] <- SeuratObject::CreateDimReducObject(
    embeddings = embedding, key = "SMALL_", assay = "RNA"
  )

  expected <- "Cannot compute 15 neighbors from only 14 cells"
  expect_error(
    entropy_neighbors(small, reduction = "small", dims = 2, k = 15),
    expected, fixed = TRUE
  )
  expect_error(
    run_entropy(
      small, assay = "RNA", reduction = "small", dims = 2, k = 15,
      mode = "binary"
    ),
    expected, fixed = TRUE
  )
})


test_that("entropy_neighbors enforces the inferred k floor", {
  skip_if_not_installed("Seurat")
  library(Seurat)
  data("pbmc_small")

  cells <- colnames(pbmc_small)
  nn <- setNames(lapply(seq_along(cells), function(i) {
    setdiff(seq_along(cells), i)[seq_len(5)]
  }), cells)

  expect_error(entropy_neighbors(pbmc_small, neighbors = nn), "below 11")
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
