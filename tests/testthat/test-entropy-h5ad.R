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
