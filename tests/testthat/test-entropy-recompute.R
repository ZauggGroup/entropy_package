make_recompute_fixture <- function(n_per_group = 15L) {
  groups <- rep(LETTERS[1:4], each = n_per_group)
  n_cells <- length(groups)
  cells <- paste0("cell_", seq_len(n_cells))
  genes <- paste0("gene", seq_len(24))
  set.seed(501)
  counts <- matrix(
    stats::rpois(length(genes) * n_cells, lambda = 3),
    nrow = length(genes),
    dimnames = list(genes, cells)
  )
  counts[, seq_len(4)] <- 1
  counts <- Matrix::Matrix(counts, sparse = TRUE)
  seu <- Seurat::CreateSeuratObject(counts = counts)
  seu$cluster <- groups
  embedding <- matrix(
    stats::rnorm(n_cells * 5), nrow = n_cells,
    dimnames = list(cells, paste0("PC_", seq_len(5)))
  )
  seu[["pca_test"]] <- SeuratObject::CreateDimReducObject(
    embeddings = embedding, key = "PC_", assay = "RNA"
  )
  entropy <- rep(c(4, 3, 2, 1), each = n_per_group) +
    stats::rnorm(n_cells, sd = 0.01)
  seu$mock_after <- entropy
  list(
    seu = seu,
    entropy = entropy,
    args = list(
      assay = "RNA", k = 11, reduction = "pca_test", dims = 5,
      mode = "binary", chunk_method = "none", parallelize = FALSE
    )
  )
}


mock_recompute_from_metadata <- function() {
  testthat::local_mocked_bindings(
    .entropy_run_recompute = function(
        seu, entropy_args, neighbors, recompute_neighbors, features = NULL,
        rebuild_reduction = FALSE) {
      values <- as.numeric(seu$mock_after)
      names(values) <- colnames(seu)
      values
    },
    .package = "Entropy",
    .env = parent.frame()
  )
}


test_that("entropy_recompute_check validates its main arguments", {
  skip_if_not_installed("Seurat")
  fixture <- make_recompute_fixture()

  expect_error(
    entropy_recompute_check(
      fixture$seu, fixture$entropy, "cluster", mechanism = "unknown",
      entropy_args = fixture$args, verbose = FALSE
    ),
    "should be one of"
  )
  expect_error(
    entropy_recompute_check(
      fixture$seu, fixture$entropy, "cluster",
      entropy_args = list(), verbose = FALSE
    ),
    "include.*assay"
  )
  expect_error(
    entropy_recompute_check(
      fixture$seu, fixture$entropy, "cluster", mechanism = "remove_genes",
      entropy_args = fixture$args, genes = NULL, verbose = FALSE
    ),
    "genes"
  )
  expect_error(
    entropy_recompute_check(
      fixture$seu, fixture$entropy, "cluster", mechanism = "subsample_cells",
      entropy_args = fixture$args, recompute_neighbors = FALSE,
      verbose = FALSE
    ),
    "not valid"
  )
  expect_error(
    entropy_recompute_check(
      fixture$seu, fixture$entropy, "cluster",
      mechanism = "downsample_counts", entropy_args = fixture$args,
      target_depth = 50, drop_below_target = TRUE,
      recompute_neighbors = FALSE, verbose = FALSE
    ),
    "not valid"
  )
})


test_that("dropping shallow cells rebuilds neighbors but keeps the embedding", {
  skip_if_not_installed("Seurat")
  fixture <- make_recompute_fixture()
  seen <- NULL
  testthat::local_mocked_bindings(
    .entropy_run_recompute = function(
        seu, entropy_args, neighbors, recompute_neighbors, features = NULL,
        rebuild_reduction = FALSE) {
      seen <<- list(
        cells = colnames(seu), neighbors = neighbors,
        recompute_neighbors = recompute_neighbors,
        rebuild_reduction = rebuild_reduction
      )
      values <- as.numeric(seu$mock_after)
      names(values) <- colnames(seu)
      values
    },
    .package = "Entropy"
  )

  suppressWarnings(entropy_recompute_check(
    fixture$seu, fixture$entropy, "cluster",
    mechanism = "downsample_counts", entropy_args = fixture$args,
    target_depth = 50, drop_below_target = TRUE,
    n_repeats = 1, verbose = FALSE
  ))

  expect_true(seen$recompute_neighbors)
  expect_null(seen$neighbors)
  expect_false(seen$rebuild_reduction)
  expect_lt(length(seen$cells), ncol(fixture$seu))
})


test_that("downsampling reaches the target and records shallow cells", {
  skip_if_not_installed("Seurat")
  fixture <- make_recompute_fixture()
  mock_recompute_from_metadata()
  totals <- Matrix::colSums(SeuratObject::GetAssayData(
    fixture$seu, assay = "RNA", layer = "counts"
  ))
  target <- 50L
  result <- suppressWarnings(entropy_recompute_check(
    fixture$seu, fixture$entropy, "cluster",
    mechanism = "downsample_counts", entropy_args = fixture$args,
    target_depth = target, n_repeats = 2, verbose = FALSE
  ))

  expect_equal(result$provenance$n_below_target, sum(totals < target))
  expect_true(all(result$intervention_totals$total_after <= target))
  deep <- result$intervention_totals$total_before > target
  expect_true(all(result$intervention_totals$total_after[deep] == target))
})


test_that("gene removal is deterministic and excludes requested genes", {
  skip_if_not_installed("Seurat")
  fixture <- make_recompute_fixture()
  seen_features <- NULL
  testthat::local_mocked_bindings(
    .entropy_run_recompute = function(
        seu, entropy_args, neighbors, recompute_neighbors, features = NULL,
        rebuild_reduction = FALSE) {
      seen_features <<- features
      values <- as.numeric(seu$mock_after)
      names(values) <- colnames(seu)
      values
    },
    .package = "Entropy"
  )
  result <- entropy_recompute_check(
    fixture$seu, fixture$entropy, "cluster", mechanism = "remove_genes",
    entropy_args = fixture$args, genes = c("gene1", "gene2"),
    n_repeats = 4, verbose = FALSE
  )

  expect_equal(result$provenance$n_repeats, 1L)
  expect_false(any(c("gene1", "gene2") %in% seen_features))
  expect_true(all(is.na(result$cluster_summary$delta)))
  expect_false(result$provenance$absolute_delta_comparable)
})


test_that("a depth-driven leading cluster collapses after intervention", {
  skip_if_not_installed("Seurat")
  fixture <- make_recompute_fixture()
  fixture$seu$mock_after <- rep(c(0, 3, 2, 1), each = 15)
  mock_recompute_from_metadata()
  result <- suppressWarnings(entropy_recompute_check(
    fixture$seu, fixture$entropy, "cluster",
    mechanism = "downsample_counts", entropy_args = fixture$args,
    target_depth = 50, n_repeats = 2, rank_tolerance = 0,
    verbose = FALSE
  ))
  verdict <- result$cluster_verdict$verdict[
    result$cluster_verdict$cluster == "A"
  ]
  expect_equal(verdict, "collapses")
})


test_that("an intervention-independent cluster ranking survives", {
  skip_if_not_installed("Seurat")
  fixture <- make_recompute_fixture()
  mock_recompute_from_metadata()
  result <- suppressWarnings(entropy_recompute_check(
    fixture$seu, fixture$entropy, "cluster",
    mechanism = "downsample_counts", entropy_args = fixture$args,
    target_depth = 50, n_repeats = 2, rank_tolerance = 0,
    verbose = FALSE
  ))
  expect_true(all(result$cluster_verdict$verdict == "survives"))
})


test_that("a cluster that consistently rises is labelled improves", {
  skip_if_not_installed("Seurat")
  fixture <- make_recompute_fixture()
  fixture$seu$mock_after <- rep(c(4, 3, 2, 5), each = 15)
  mock_recompute_from_metadata()
  result <- suppressWarnings(entropy_recompute_check(
    fixture$seu, fixture$entropy, "cluster",
    mechanism = "downsample_counts", entropy_args = fixture$args,
    target_depth = 50, n_repeats = 2, rank_tolerance = 0,
    verbose = FALSE
  ))
  verdict <- result$cluster_verdict$verdict[
    result$cluster_verdict$cluster == "D"
  ]
  expect_equal(verdict, "improves")
  expect_false(verdict == "unstable")
  expect_output(print(result), "1 improves")
})


test_that("subsampling varies by seed and restores the caller RNG", {
  skip_if_not_installed("Seurat")
  fixture <- make_recompute_fixture()
  mock_recompute_from_metadata()
  set.seed(812)
  caller_seed <- .Random.seed
  result <- entropy_recompute_check(
    fixture$seu, fixture$entropy, "cluster",
    mechanism = "subsample_cells", entropy_args = fixture$args,
    cells_per_group = 10, n_repeats = 2, seed = 9,
    rank_tolerance = 0, verbose = FALSE
  )

  first <- result$cell_data$cell[result$cell_data$repeat_id == 1]
  second <- result$cell_data$cell[result$cell_data$repeat_id == 2]
  expect_false(setequal(first, second))
  expect_identical(.Random.seed, caller_seed)
  expect_false(is.null(result$density_check))
})


test_that("density checks are returned only for cell subsampling", {
  skip_if_not_installed("Seurat")
  fixture <- make_recompute_fixture()
  mock_recompute_from_metadata()
  downsample <- suppressWarnings(entropy_recompute_check(
    fixture$seu, fixture$entropy, "cluster",
    mechanism = "downsample_counts", entropy_args = fixture$args,
    target_depth = 50, n_repeats = 1, verbose = FALSE
  ))
  removed <- entropy_recompute_check(
    fixture$seu, fixture$entropy, "cluster", mechanism = "remove_genes",
    entropy_args = fixture$args, genes = "gene1", verbose = FALSE
  )
  expect_null(downsample$density_check)
  expect_null(removed$density_check)
})


test_that("a diagnostic connects to its declared intervention", {
  skip_if_not_installed("Seurat")
  fixture <- make_recompute_fixture()
  mock_recompute_from_metadata()
  diagnostic <- list(
    cell_data = data.frame(
      cell = colnames(fixture$seu), entropy = fixture$entropy,
      cluster = fixture$seu$cluster
    ),
    covariate_info = data.frame(
      covariate = "depth", intervention = "downsample_counts"
    ),
    provenance = list(group_column = "cluster")
  )
  class(diagnostic) <- "entropy_diagnostic"
  result <- suppressWarnings(entropy_recompute_check(
    diagnostic, entropy_args = fixture$args, covariate = "depth",
    expression_seu = fixture$seu, target_depth = 50,
    n_repeats = 1, verbose = FALSE
  ))
  expect_s3_class(result, "entropy_recompute")
  expect_equal(result$provenance$mechanism, "downsample_counts")

  diagnostic$covariate_info$intervention <- "none"
  expect_error(
    entropy_recompute_check(
      diagnostic, entropy_args = fixture$args, covariate = "depth",
      expression_seu = fixture$seu, verbose = FALSE
    ),
    "No recompute test exists"
  )
})


test_that("entropy_recompute_report writes files and protects existing output", {
  skip_if_not_installed("Seurat")
  fixture <- make_recompute_fixture()
  mock_recompute_from_metadata()
  result <- suppressWarnings(entropy_recompute_check(
    fixture$seu, fixture$entropy, "cluster",
    mechanism = "downsample_counts", entropy_args = fixture$args,
    target_depth = 50, n_repeats = 2, verbose = FALSE
  ))
  outdir <- tempfile("entropy-recompute-report-")
  paths <- entropy_recompute_report(
    result, outdir = outdir, write_tables = TRUE
  )
  expect_true(all(file.exists(paths)))
  expect_equal(length(paths), 5L)
  expect_error(
    entropy_recompute_report(result, outdir = outdir),
    "Refusing to overwrite"
  )
  expect_no_error(entropy_recompute_report(
    result, outdir = outdir, overwrite = TRUE
  ))
})


test_that("print method reports the headline result", {
  skip_if_not_installed("Seurat")
  fixture <- make_recompute_fixture()
  mock_recompute_from_metadata()
  result <- suppressWarnings(entropy_recompute_check(
    fixture$seu, fixture$entropy, "cluster",
    mechanism = "downsample_counts", entropy_args = fixture$args,
    target_depth = 50, n_repeats = 1, verbose = FALSE
  ))
  expect_output(print(result), "Headline rank correlation")
})
test_that("downsampled counts are returned in column-compressed storage", {
  counts <- Matrix::Matrix(
    matrix(
      as.integer(c(50, 8, 5, 3, 2, 2) * rep(1:4, each = 6)),
      nrow = 6L,
      dimnames = list(paste0("gene", 1:6), paste0("cell", 1:4))
    ),
    sparse = TRUE
  )
  set.seed(1)
  downsampled <- .entropy_downsample_matrix(counts, target = 40L)

  expect_true(methods::is(downsampled, "CsparseMatrix"))
  expect_identical(dim(downsampled), dim(counts))
  expect_identical(dimnames(downsampled), dimnames(counts))
  expect_true(all(Matrix::colSums(downsampled) <= Matrix::colSums(counts)))
})


test_that(".entropy_as_csparse preserves dimensions and dimnames", {
  dense <- matrix(
    c(1, 0, 2, 0, 3, 0),
    nrow = 2L,
    dimnames = list(c("a", "b"), c("x", "y", "z"))
  )
  converted <- .entropy_as_csparse(dense)

  expect_true(methods::is(converted, "CsparseMatrix"))
  expect_identical(dim(converted), c(2L, 3L))
  expect_identical(dimnames(converted), dimnames(dense))
  expect_equal(as.matrix(converted), dense)

  already <- Matrix::Matrix(dense, sparse = TRUE)
  expect_identical(.entropy_as_csparse(already), already)
})


test_that("rowCounts_CSC_cols rejects non-CSC input and bad columns", {
  sparse <- Matrix::Matrix(
    matrix(c(TRUE, FALSE, TRUE, TRUE), nrow = 2L),
    sparse = TRUE
  )
  expect_error(
    rowCounts_CSC_cols(as.matrix(sparse), 1L),
    "column-compressed sparse matrix"
  )
  expect_error(
    rowCounts_CSC_cols(sparse, 5L),
    "Column indices must lie between"
  )
  expect_identical(rowCounts_CSC_cols(sparse, 1:2), c(2L, 1L))
})
