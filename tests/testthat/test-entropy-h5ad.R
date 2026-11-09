expect_valid_h5ad_entropy <- function(result, n_cells = 1200L) {
  expect_s3_class(result, "data.frame")
  expect_identical(colnames(result), "Entropy_mean")
  expect_equal(nrow(result), n_cells)
  expect_equal(length(unique(rownames(result))), n_cells)
  expect_false(anyNA(result$Entropy_mean))
  expect_true(all(is.finite(result$Entropy_mean)))
  expect_true(all(result$Entropy_mean >= 0))
}


skip_if_h5ad_unavailable <- function() {
  skip_if(
    getRversion() < "4.5.0",
    "run_entropy_h5ad() requires R >= 4.5.0"
  )
  skip_if_not_installed("anndataR", minimum_version = "1.0.0")
  skip_if_not_installed("rhdf5")
}


test_that("run_entropy_h5ad explains its R version requirement", {
  skip_if(getRversion() >= "4.5.0")

  expect_error(
    run_entropy_h5ad("not-read-on-older-r.h5ad"),
    "requires R >= 4.5.0",
    fixed = TRUE
  )
})


test_that("run_entropy_h5ad works in counts mode with automatic neighbors", {
  skip_if_h5ad_unavailable()
  path <- test_path("fixtures", "pbmc_synth_small.h5ad")

  result <- run_entropy_h5ad(
    path,
    layer = "counts",
    k = 15,
    mode = "counts"
  )

  expect_valid_h5ad_entropy(result)
  expect_true(any(result$Entropy_mean > 0))
})


test_that("run_entropy_h5ad restricts automatic neighbors within celltype", {
  skip_if_h5ad_unavailable()
  path <- test_path("fixtures", "pbmc_synth_small.h5ad")

  result <- run_entropy_h5ad(
    path,
    layer = "counts",
    restrict_by = "celltype"
  )

  expect_valid_h5ad_entropy(result)
  expect_true(any(result$Entropy_mean > 0))
})


test_that("run_entropy_h5ad chunked computation matches unchunked computation", {
  skip_if_h5ad_unavailable()
  path <- test_path("fixtures", "pbmc_synth_small.h5ad")

  unchunked <- run_entropy_h5ad(
    path,
    layer = "counts",
    k = 15,
    mode = "counts"
  )
  chunked <- run_entropy_h5ad(
    path,
    layer = "counts",
    chunk_method = "auto",
    chunk_cells = 500
  )

  expect_valid_h5ad_entropy(chunked)
  expect_equal(chunked, unchunked, tolerance = 1e-12)
})


test_that("run_entropy_h5ad works in binary mode", {
  skip_if_h5ad_unavailable()
  path <- test_path("fixtures", "pbmc_synth_small.h5ad")

  result <- run_entropy_h5ad(
    path,
    layer = "counts",
    mode = "binary"
  )

  expect_valid_h5ad_entropy(result)
  expect_true(any(result$Entropy_mean > 0))
  expect_true(all(result$Entropy_mean <= 1))
})
test_that("run_entropy_h5ad returns its inputs on request", {
  skip_if_h5ad_unavailable()
  path <- testthat::test_path("fixtures", "pbmc_synth_small.h5ad")

  result <- run_entropy_h5ad(
    path, layer = "counts", k = 15, mode = "counts", return_inputs = TRUE
  )

  expect_type(result, "list")
  expect_named(
    result,
    c("entropy", "entropy_depth_corrected", "embedding", "obs", "results")
  )
  expect_null(result$entropy_depth_corrected)
  expect_true(is.numeric(result$entropy))
  expect_identical(names(result$entropy), rownames(result$results))
  expect_identical(rownames(result$embedding), rownames(result$results))
  expect_identical(rownames(result$obs), rownames(result$results))
  expect_equal(unname(result$entropy), result$results$Entropy_mean)
})


test_that("run_entropy_h5ad depth correction adds the expected columns", {
  skip_if_h5ad_unavailable()
  path <- testthat::test_path("fixtures", "pbmc_synth_small.h5ad")

  plain <- run_entropy_h5ad(path, layer = "counts", k = 15, mode = "counts")

  set.seed(11)
  corrected <- suppressWarnings(run_entropy_h5ad(
    path, layer = "counts", k = 15, mode = "counts",
    depth_correct = TRUE, depth_quantile = 0.05
  ))

  expect_identical(
    colnames(corrected),
    c("Entropy_mean", "Entropy_mean_depth_corrected", "depth_corrected")
  )
  expect_identical(corrected$Entropy_mean, plain$Entropy_mean)
  expect_identical(rownames(corrected), rownames(plain))
  expect_type(corrected$depth_corrected, "logical")
  expect_true(all(is.finite(corrected$Entropy_mean_depth_corrected)))
})


test_that("run_entropy_h5ad validates the depth correction arguments", {
  skip_if_h5ad_unavailable()
  path <- testthat::test_path("fixtures", "pbmc_synth_small.h5ad")

  expect_error(
    run_entropy_h5ad(path, layer = "counts", depth_correct = "yes"),
    "must be TRUE or FALSE"
  )
  expect_error(
    run_entropy_h5ad(path, layer = "counts", depth_correct = TRUE,
                     depth_quantile = 0),
    "strictly between 0 and 1"
  )
  expect_error(
    run_entropy_h5ad(path, layer = "counts", depth_correct = TRUE,
                     target_depth = -5),
    "positive integer"
  )
})


test_that("run_entropy_h5ad enforces custom-neighbour invariants", {
  skip_if_h5ad_unavailable()
  path <- test_path("fixtures", "pbmc_synth_small.h5ad")
  cells <- read_h5ad_inputs(
    path, layer = "counts", need_embedding = FALSE, read_obs = FALSE
  )$cell_names
  n <- length(cells)
  neighbors <- setNames(lapply(seq_len(n), function(i) {
    ((i - 1L + seq_len(11L)) %% n) + 1L
  }), cells)

  unequal <- neighbors
  unequal[[1L]] <- c(unequal[[1L]], 13L)
  expect_error(
    run_entropy_h5ad(
      path, layer = "counts", neighbors = unequal,
      neighbors_method = "custom", k = 11
    ),
    "different sizes"
  )

  self <- neighbors
  self[[1L]][1L] <- 1L
  expect_error(
    run_entropy_h5ad(
      path, layer = "counts", neighbors = self,
      neighbors_method = "custom", k = 11
    ),
    "their own neighbour"
  )

  repeated <- neighbors
  repeated[[1L]][2L] <- repeated[[1L]][1L]
  expect_error(
    run_entropy_h5ad(
      path, layer = "counts", neighbors = repeated,
      neighbors_method = "custom", k = 11
    ),
    "repeat a neighbour"
  )

  fractional <- neighbors
  fractional[[1L]][1L] <- 2.5
  expect_error(
    run_entropy_h5ad(
      path, layer = "counts", neighbors = fractional,
      neighbors_method = "custom", k = 11
    ),
    "finite integers"
  )

  too_small <- lapply(neighbors, utils::head, 5L)
  expect_error(
    run_entropy_h5ad(
      path, layer = "counts", neighbors = too_small,
      neighbors_method = "custom"
    ),
    "below 11"
  )
})


test_that("run_entropy_h5ad rejects partial restrict_by pruning", {
  skip_if_h5ad_unavailable()
  path <- test_path("fixtures", "pbmc_synth_small.h5ad")
  cells <- read_h5ad_inputs(
    path, layer = "counts", restrict_by = "celltype",
    need_embedding = FALSE, read_obs = FALSE
  )$cell_names
  n <- length(cells)
  neighbors <- setNames(lapply(seq_len(n), function(i) {
    ((i - 1L + seq_len(11L)) %% n) + 1L
  }), cells)

  expect_error(
    run_entropy_h5ad(
      path, layer = "counts", neighbors = neighbors,
      neighbors_method = "custom", k = 11, restrict_by = "celltype"
    ),
    "Drop `neighbors`"
  )
})
