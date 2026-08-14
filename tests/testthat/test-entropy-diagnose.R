make_diagnostic_fixture <- function(group, covariate, entropy, seed = 1L) {
  set.seed(seed)
  n <- length(group)
  counts <- Matrix::Matrix(matrix(
    stats::rpois(5L * n, lambda = 2),
    nrow = 5L,
    dimnames = list(paste0("gene", seq_len(5L)), paste0("cell", seq_len(n)))
  ), sparse = TRUE)
  seu <- Seurat::CreateSeuratObject(counts = counts)
  seu$cluster <- factor(group, levels = unique(group))
  embedding <- cbind(stats::rnorm(n), stats::rnorm(n))
  rownames(embedding) <- colnames(seu)
  names(entropy) <- colnames(seu)
  list(
    seu = seu,
    entropy = entropy,
    embedding = embedding,
    covariates = data.frame(test_covariate = covariate)
  )
}


run_fixture_diagnostic <- function(fixture, ...) {
  arguments <- list(
    seu = fixture$seu,
    entropy = fixture$entropy,
    embedding = fixture$embedding,
    k = 5,
    group = "cluster",
    covariates = fixture$covariates,
    n_bins = 10,
    min_bin_n = 5,
    min_cluster_n = 10,
    verbose = FALSE
  )
  arguments[names(list(...))] <- list(...)
  do.call(entropy_diagnose, arguments)
}


test_that("entropy_diagnose validates core inputs", {
  skip_if_not_installed("Seurat")
  fixture <- make_diagnostic_fixture(
    rep(c("A", "B"), each = 30),
    seq_len(60),
    seq_len(60) / 10
  )

  expect_error(
    entropy_diagnose(
      fixture$seu, fixture$entropy[-1], fixture$embedding,
      k = 5, group = "cluster", verbose = FALSE
    ),
    "one row per value"
  )
  expect_error(
    entropy_diagnose(
      fixture$seu, as.character(fixture$entropy), fixture$embedding,
      k = 5, group = "cluster", verbose = FALSE
    ),
    "numeric vector"
  )
  expect_error(
    entropy_diagnose(
      fixture$seu, fixture$entropy, fixture$embedding,
      k = 5, group = "missing", verbose = FALSE
    ),
    "was not found"
  )
  expect_error(
    entropy_diagnose(
      fixture$seu, fixture$entropy, fixture$embedding,
      k = 5, group = "cluster", covariates = list(bad = letters[1:60]),
      verbose = FALSE
    ),
    "must be numeric"
  )
  expect_error(
    entropy_diagnose(
      fixture$seu, fixture$entropy, fixture$embedding,
      k = 5, group = "cluster", covariates = fixture$covariates,
      covariate_type = c(test_covariate = "mechanical"), verbose = FALSE
    ),
    "must be one of"
  )
  expect_error(
    entropy_diagnose(
      fixture$seu, fixture$entropy, fixture$embedding,
      k = 5, group = "cluster", max_bin_imbalance = -0.1,
      verbose = FALSE
    ),
    "max_bin_imbalance"
  )
  expect_error(
    entropy_diagnose(
      fixture$seu, fixture$entropy, fixture$embedding,
      k = 5, group = "cluster", retained_min_delta = -0.1,
      verbose = FALSE
    ),
    "retained_min_delta"
  )
})


test_that("density is always present", {
  skip_if_not_installed("Seurat")
  fixture <- make_diagnostic_fixture(
    rep(c("A", "B"), each = 30),
    seq_len(60),
    stats::rnorm(60)
  )
  result <- entropy_diagnose(
    fixture$seu, fixture$entropy, fixture$embedding,
    k = 5, group = "cluster", min_bin_n = 2, verbose = FALSE
  )
  expect_identical(result$covariate_info$covariate, "mean_knn_dist")
  expect_identical(result$covariate_info$type, "unknown")
  expect_identical(result$covariate_info$intervention, "subsample_cells")
})


test_that("a cluster without overlap receives cannot tell", {
  skip_if_not_installed("Seurat")
  n <- 80L
  group <- rep(c("isolated", "other"), each = n / 2L)
  covariate <- c(seq(10, 11, length.out = n / 2L), seq(0, 1, length.out = n / 2L))
  entropy <- c(rep(2, n / 2L), rep(0, n / 2L)) + seq_len(n) / 10000
  result <- run_fixture_diagnostic(
    make_diagnostic_fixture(group, covariate, entropy),
    n_bins = 8,
    min_bin_n = 5
  )
  row <- result$cluster_table[
    result$cluster_table$cluster == "isolated" &
      result$cluster_table$covariate == "test_covariate",
    , drop = FALSE
  ]
  expect_identical(row$flag, "cannot tell")
  expect_lt(row$support, result$thresholds$support_threshold)
})


test_that("a tail cluster fails the within-bin balance check", {
  skip_if_not_installed("Seurat")
  n_cluster <- 50L
  n_other <- 950L
  group <- c(rep("tail", n_cluster), rep("other", n_other))
  covariate <- c(
    seq(20, 100, length.out = n_cluster),
    seq(8, 25, length.out = n_other)
  )
  entropy <- c(
    seq(2, 2.1, length.out = n_cluster),
    seq(0, 0.1, length.out = n_other)
  )
  fixture <- make_diagnostic_fixture(group, covariate, entropy)
  result <- entropy_diagnose(
    seu = fixture$seu,
    entropy = fixture$entropy,
    embedding = fixture$embedding,
    k = 5,
    group = "cluster",
    covariates = fixture$covariates,
    n_bins = 10,
    verbose = FALSE
  )
  row <- result$cluster_table[
    result$cluster_table$cluster == "tail" &
      result$cluster_table$covariate == "test_covariate",
    , drop = FALSE
  ]
  tail_bins <- result$bin_table[
    result$bin_table$cluster == "tail" &
      result$bin_table$covariate == "test_covariate",
    , drop = FALSE
  ]

  expect_identical(row$flag, "cannot tell")
  expect_lt(row$support, result$thresholds$support_threshold)
  expect_true(any(
    tail_bins$n_in >= result$thresholds$min_bin_n &
      tail_bins$n_out >= result$thresholds$min_bin_n &
      tail_bins$imbalance > result$thresholds$max_bin_imbalance
  ))
  expect_true(all(
    c("median_cov_in", "median_cov_out", "imbalance") %in% names(result$bin_table)
  ))
})


test_that("a cluster without elevated entropy receives not elevated", {
  skip_if_not_installed("Seurat")
  values <- seq(-1, 1, length.out = 60)
  fixture <- make_diagnostic_fixture(
    rep(c("A", "B"), each = 60),
    rep(values, 2),
    rep(values, 2)
  )
  result <- run_fixture_diagnostic(fixture)
  row <- result$cluster_table[
    result$cluster_table$cluster == "A" &
      result$cluster_table$covariate == "test_covariate",
    , drop = FALSE
  ]
  expect_identical(row$flag, "not elevated")
})


test_that("an independent cluster entropy shift receives not explained", {
  skip_if_not_installed("Seurat")
  set.seed(11)
  n_each <- 120L
  group <- rep(c("high", "other"), each = n_each)
  shared <- seq(-2, 2, length.out = n_each)
  covariate <- rep(shared, 2)
  entropy <- covariate + rep(c(0.8, 0), each = n_each) + stats::rnorm(2 * n_each, sd = 0.02)
  result <- run_fixture_diagnostic(make_diagnostic_fixture(group, covariate, entropy))
  row <- result$cluster_table[
    result$cluster_table$cluster == "high" &
      result$cluster_table$covariate == "test_covariate",
    , drop = FALSE
  ]
  expect_identical(row$flag, "not explained")
  expect_gte(row$support, result$thresholds$support_threshold)
})


test_that("covariate-generated cluster entropy receives explained", {
  skip_if_not_installed("Seurat")
  n_each <- 200L
  group <- rep(c("high", "other"), each = n_each)
  base <- seq(-2, 2, length.out = n_each)
  covariate <- c(base + 0.35, base)
  entropy <- covariate
  result <- run_fixture_diagnostic(
    make_diagnostic_fixture(group, covariate, entropy),
    n_bins = 20,
    min_bin_n = 3,
    effect_threshold = 0.1
  )
  row <- result$cluster_table[
    result$cluster_table$cluster == "high" &
      result$cluster_table$covariate == "test_covariate",
    , drop = FALSE
  ]
  expect_true(row$entropy_elevated)
  expect_identical(row$flag, "explained")
})


test_that("a reversed matched difference is explained and marked", {
  skip_if_not_installed("Seurat")
  n_each <- 400L
  group <- rep(c("shifted", "other"), each = n_each)
  base <- seq(-2, 2, length.out = n_each)
  covariate <- c(base + 0.8, base)
  entropy <- c(base + 0.8 - 0.4, base)
  result <- suppressWarnings(run_fixture_diagnostic(
    make_diagnostic_fixture(group, covariate, entropy),
    n_bins = 20,
    min_bin_n = 5,
    max_bin_imbalance = 0.1,
    retained_min_delta = 0.1
  ))
  row <- result$cluster_table[
    result$cluster_table$cluster == "shifted" &
      result$cluster_table$covariate == "test_covariate",
    , drop = FALSE
  ]

  expect_true(row$entropy_elevated)
  expect_lt(row$matched_diff_std, -result$thresholds$effect_threshold)
  expect_equal(row$retained, row$matched_diff_std / row$entropy_delta_std)
  expect_identical(row$flag, "explained")
  expect_true(row$reversed)
})


test_that("retained is missing when entropy elevation is zero", {
  skip_if_not_installed("Seurat")
  group <- rep(c("A", "B"), each = 60)
  entropy <- rep(c(-1, 1), 60)
  covariate <- rep(seq(-2, 2, length.out = 60), 2)
  result <- suppressWarnings(run_fixture_diagnostic(
    make_diagnostic_fixture(group, covariate, entropy),
    n_bins = 10,
    min_bin_n = 3
  ))

  expect_true(all(result$cluster_table$entropy_delta_std == 0))
  expect_true(all(is.na(result$cluster_table$retained)))
})


test_that("retained and reversed are missing near the global mean", {
  skip_if_not_installed("Seurat")
  n_each <- 120L
  group <- rep(c("slightly_high", "slightly_low"), each = n_each)
  base <- seq(-2, 2, length.out = n_each)
  entropy <- c(base + 0.05, base)
  covariate <- c(base, base)
  result <- suppressWarnings(run_fixture_diagnostic(
    make_diagnostic_fixture(group, covariate, entropy),
    n_bins = 12,
    min_bin_n = 3,
    elevated_threshold = 0.1
  ))
  rows <- result$cluster_table[
    result$cluster_table$covariate == "test_covariate", , drop = FALSE
  ]

  expect_true(all(abs(rows$entropy_delta_std) <= 0.1))
  expect_true(any(rows$entropy_delta_std != 0))
  expect_true(all(is.na(rows$retained)))
  expect_true(all(is.na(rows$reversed)))
  expect_identical(as.character(rows$direction), c("elevated", "depleted"))
})


test_that("retained is missing for a marginal standardized entropy difference", {
  skip_if_not_installed("Seurat")
  n_each <- 120L
  group <- rep(c("marginal_high", "marginal_low"), each = n_each)
  base <- seq(-2, 2, length.out = n_each)
  entropy <- c(base + 0.32, base)
  covariate <- c(base, base)
  result <- suppressWarnings(run_fixture_diagnostic(
    make_diagnostic_fixture(group, covariate, entropy),
    n_bins = 12,
    min_bin_n = 3,
    elevated_threshold = 0.1
  ))
  row <- result$cluster_table[
    result$cluster_table$cluster == "marginal_high" &
      result$cluster_table$covariate == "test_covariate",
    , drop = FALSE
  ]

  expect_equal(result$thresholds$retained_min_delta, 0.5)
  expect_equal(row$entropy_delta_std, 0.27, tolerance = 0.03)
  expect_gt(abs(row$entropy_delta_std), result$thresholds$elevated_threshold)
  expect_lt(abs(row$entropy_delta_std), result$thresholds$retained_min_delta)
  expect_true(is.na(row$retained))
})


test_that("depleted clusters report retained deficit but not reversed", {
  skip_if_not_installed("Seurat")
  n_each <- 160L
  group <- rep(c("high", "low"), each = n_each)
  base <- seq(-2, 2, length.out = n_each)
  result <- suppressWarnings(run_fixture_diagnostic(
    make_diagnostic_fixture(
      group, c(base, base), c(base + 0.5, base)
    ),
    n_bins = 16,
    min_bin_n = 3,
    retained_min_delta = 0.1
  ))
  row <- result$cluster_table[
    result$cluster_table$cluster == "low" &
      result$cluster_table$covariate == "test_covariate",
    , drop = FALSE
  ]

  expect_identical(row$direction, "depleted")
  expect_lt(row$entropy_delta_std, -result$thresholds$elevated_threshold)
  expect_true(is.finite(row$retained))
  expect_true(is.na(row$reversed))
})


test_that("strong covariates warn when bin balance is loose", {
  skip_if_not_installed("Seurat")
  values <- rep(seq(-2, 2, length.out = 100), 2)
  fixture <- make_diagnostic_fixture(
    rep(c("A", "B"), each = 100), values, values
  )
  expect_warning(
    run_fixture_diagnostic(fixture, max_bin_imbalance = 0.25),
    "Tighten max_bin_imbalance"
  )
})


test_that("matched difference is missing when no bins are usable", {
  skip_if_not_installed("Seurat")
  fixture <- make_diagnostic_fixture(
    rep(c("A", "B"), each = 20),
    c(seq(10, 11, length.out = 20), seq(0, 1, length.out = 20)),
    c(rep(2, 20), rep(0, 20))
  )
  result <- run_fixture_diagnostic(fixture, min_bin_n = 15)
  row <- result$cluster_table[
    result$cluster_table$cluster == "A" &
      result$cluster_table$covariate == "test_covariate",
    , drop = FALSE
  ]
  expect_true(is.na(row$matched_diff_std))
})


test_that("joint model reports shared explanatory power", {
  skip_if_not_installed("Seurat")
  set.seed(12)
  n <- 200L
  first <- stats::rnorm(n)
  second <- first + stats::rnorm(n, sd = 0.1)
  fixture <- make_diagnostic_fixture(
    rep(c("A", "B"), each = n / 2),
    first,
    first + stats::rnorm(n, sd = 0.1)
  )
  fixture$covariates$second <- second
  result <- run_fixture_diagnostic(fixture)
  individual_sum <- sum(result$covariate_info$r_squared)
  expect_lte(result$joint$joint_r_squared, individual_sum + 1e-12)
  expect_named(result$joint$incremental_r_squared, result$covariate_info$covariate)
  expect_named(result$joint$vif, result$covariate_info$covariate)
})


test_that("provenance k mismatch warns", {
  skip_if_not_installed("Seurat")
  fixture <- make_diagnostic_fixture(
    rep(c("A", "B"), each = 30),
    seq_len(60),
    stats::rnorm(60)
  )
  expect_warning(
    run_fixture_diagnostic(fixture, provenance = list(k = 10)),
    "records k = 10"
  )
})


test_that("loess subsampling is random, bounded, and preserves RNG state", {
  n <- 20000L
  predictor <- seq(-5, 5, length.out = n)
  entropy <- 0.5 * predictor + sin(predictor)
  set.seed(987)
  seed_before <- .Random.seed

  result <- Entropy:::.entropy_loess_curve(entropy, predictor)

  expect_identical(.Random.seed, seed_before)
  expected_seed <- .Random.seed
  set.seed(1L)
  keep <- sort(sample.int(n, 10000L))
  assign(".Random.seed", expected_seed, envir = .GlobalEnv)
  expect_equal(range(result$curve$x), range(predictor[keep]))
  expect_true(all(is.finite(result$curve$y[c(1L, nrow(result$curve))])))
  expect_true(result$subsampled)
})


test_that("print method runs and includes flag numbers", {
  skip_if_not_installed("Seurat")
  fixture <- make_diagnostic_fixture(
    rep(c("A", "B"), each = 30),
    seq_len(60),
    stats::rnorm(60)
  )
  result <- run_fixture_diagnostic(fixture)
  expect_output(print(result), "Entropy diagnostic")
  expect_output(print(result), "Use entropy_diagnose_report")
})


test_that("report writes expected files and protects existing files", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  fixture <- make_diagnostic_fixture(
    rep(c("A", "B"), each = 30),
    rep(seq(-1, 1, length.out = 30), 2),
    stats::rnorm(60)
  )
  result <- run_fixture_diagnostic(fixture)
  outdir <- tempfile("entropy-report-")
  paths <- entropy_diagnose_report(result, outdir = outdir, write_tables = TRUE)

  expect_length(paths, nrow(result$covariate_info) + 4L)
  expect_true(all(file.exists(paths)))
  expect_true(all(
    c("cluster_table", "bin_table", "covariate_info") %in% names(paths)
  ))
  expect_error(
    entropy_diagnose_report(result, outdir = outdir, overwrite = FALSE),
    "Refusing to overwrite"
  )
})


test_that("large reports paginate 20 clusters and 6 covariates", {
  skip_if_not_installed("Seurat")
  set.seed(71)
  n_clusters <- 20L
  cells_per_cluster <- 25L
  n <- n_clusters * cells_per_cluster
  group <- rep(sprintf("cluster_%02d", seq_len(n_clusters)), each = cells_per_cluster)
  entropy <- rep(rev(seq_len(n_clusters)) / 5, each = cells_per_cluster) +
    stats::rnorm(n, sd = 0.05)
  fixture <- make_diagnostic_fixture(group, stats::rnorm(n), entropy)
  fixture$covariates <- data.frame(
    covariate_1 = stats::rnorm(n),
    covariate_2 = stats::rnorm(n),
    covariate_3 = stats::rnorm(n),
    covariate_4 = stats::rnorm(n),
    covariate_5 = stats::rnorm(n)
  )
  result <- entropy_diagnose(
    seu = fixture$seu,
    entropy = fixture$entropy,
    embedding = fixture$embedding,
    k = 5,
    group = "cluster",
    covariates = fixture$covariates,
    n_bins = 5,
    min_bin_n = 2,
    min_cluster_n = 10,
    verbose = FALSE
  )

  expect_equal(nrow(result$covariate_info), 6L)
  expect_length(Entropy:::.entropy_summary_plots(result), 2L)
  expect_length(
    Entropy:::.entropy_covariate_pages(result, "covariate_1"),
    16L
  )
  expected_order <- sprintf("cluster_%02d", seq_len(n_clusters))
  expect_identical(Entropy:::.entropy_cluster_order(result), expected_order)
  expect_match(
    paste(Entropy:::.entropy_caveat_text("density", "subsample_cells"), collapse = " "),
    "neither purely technical nor purely biological"
  )
})
test_that("entropy_diagnose runs without a Seurat object", {
  skip_if_not_installed("Seurat")
  fixture <- make_diagnostic_fixture(
    rep(c("A", "B"), each = 60),
    rep(seq(-1, 1, length.out = 60), 2),
    c(seq(0, 1, length.out = 60), seq(0.5, 1.5, length.out = 60))
  )

  with_object <- entropy_diagnose(
    seu = fixture$seu, entropy = fixture$entropy,
    embedding = fixture$embedding, k = 5, group = "cluster",
    covariates = fixture$covariates, min_bin_n = 5, min_cluster_n = 10,
    verbose = FALSE
  )

  group_vector <- fixture$seu$cluster
  without_object <- entropy_diagnose(
    seu = NULL, entropy = fixture$entropy,
    embedding = fixture$embedding, k = 5, group = group_vector,
    covariates = fixture$covariates, min_bin_n = 5, min_cluster_n = 10,
    verbose = FALSE
  )

  expect_identical(
    without_object$cluster_table$flag,
    with_object$cluster_table$flag
  )
  expect_equal(
    without_object$cluster_table$matched_diff_std,
    with_object$cluster_table$matched_diff_std
  )
  expect_identical(without_object$provenance$group_column, "group")
})


test_that("entropy_diagnose rejects a column name when seu is NULL", {
  skip_if_not_installed("Seurat")
  fixture <- make_diagnostic_fixture(
    rep(c("A", "B"), each = 30),
    seq_len(60),
    stats::rnorm(60)
  )
  expect_error(
    entropy_diagnose(
      seu = NULL, entropy = fixture$entropy, embedding = fixture$embedding,
      k = 5, group = "cluster", verbose = FALSE
    ),
    "looks like a column name"
  )
  expect_error(
    entropy_diagnose(
      seu = NULL, entropy = fixture$entropy, embedding = fixture$embedding,
      k = 5, group = rep("A", 10L), verbose = FALSE
    ),
    "one value per cell"
  )
})
