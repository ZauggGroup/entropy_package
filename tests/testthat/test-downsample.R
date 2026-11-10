test_that("built-in thinning hits the target and leaves other cells alone", {
  set.seed(1)
  counts <- Matrix::rsparsematrix(
    200, 120, 0.15, rand.x = function(n) stats::rpois(n, 3) + 1
  )
  dimnames(counts) <- list(paste0("gene", 1:200), paste0("cell", 1:120))
  totals <- Matrix::colSums(counts)
  target <- floor(stats::quantile(totals, 0.8, names = FALSE))
  affected <- totals > target

  thinned <- Entropy:::.entropy_downsample_fallback(counts, target)

  expect_s4_class(thinned, "dgCMatrix")
  expect_identical(dim(thinned), dim(counts))
  expect_identical(dimnames(thinned), dimnames(counts))
  expect_true(all(Matrix::colSums(thinned)[affected] == target))
  expect_identical(
    as.matrix(thinned[, !affected]), as.matrix(counts[, !affected])
  )
  expect_true(all(as.matrix(thinned) <= as.matrix(counts)))
  expect_true(all(thinned@x > 0))
})

test_that("built-in thinning is unbiased for the expected counts", {
  set.seed(2)
  counts <- Matrix::Matrix(
    matrix(c(40, 25, 10, 5, 0, 20), ncol = 1), sparse = TRUE
  )
  target <- 50
  draws <- replicate(3000, as.numeric(
    Entropy:::.entropy_downsample_fallback(counts, target)
  ))
  expected <- as.numeric(counts) * target / sum(counts)
  expect_lt(max(abs(rowMeans(draws) - expected)), 0.3)
})

test_that("built-in thinning accepts dense input", {
  counts <- matrix(c(10, 0, 5, 3, 3, 3), nrow = 3)
  thinned <- Entropy:::.entropy_downsample_fallback(counts, 6)
  expect_s4_class(thinned, "dgCMatrix")
  expect_equal(Matrix::colSums(thinned), c(6, 6))
})
