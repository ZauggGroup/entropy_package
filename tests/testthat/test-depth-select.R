make_depth_select_fixture <- function() {
  n_cells <- 36L
  n_genes <- 30L
  cells <- paste0("cell_", seq_len(n_cells))
  genes <- paste0("gene", seq_len(n_genes))
  detected <- rep(8:25, length.out = n_cells)
  counts <- matrix(0, nrow = n_genes, ncol = n_cells,
                   dimnames = list(genes, cells))
  for (index in seq_len(n_cells)) {
    counts[seq_len(detected[[index]]), index] <- 1L
    counts[1L, index] <- counts[1L, index] + index
  }
  seu <- Seurat::CreateSeuratObject(
    counts = Matrix::Matrix(counts, sparse = TRUE)
  )
  seu$cluster <- rep(c("A", "B", "C"), each = 12L)
  set.seed(731)
  embedding <- matrix(
    stats::rnorm(n_cells * 4L), nrow = n_cells,
    dimnames = list(cells, paste0("PC_", 1:4))
  )
  seu[["pca_test"]] <- SeuratObject::CreateDimReducObject(
    embeddings = embedding, key = "PC_", assay = "RNA"
  )
  seu
}


depth_select_mock <- function(seu, profile = c("correct", "refuse")) {
  profile <- match.arg(profile)
  counts <- SeuratObject::GetAssayData(seu, assay = "RNA", layer = "counts")
  x <- as.numeric(scale(Matrix::colSums(counts > 0)))
  raw_other <- sin(seq_along(x)) + cos(seq_along(x) / 3)
  other <- as.numeric(scale(stats::residuals(stats::lm(raw_other ~ x))))
  calls <- new.env(parent = emptyenv())
  calls$values <- list()

  components <- if (profile == "correct") {
    data.frame(q = c(1, 0.8, 0.6), confound = c(1, 0.6, 0.25),
               leftover = c(1, 0.85, 0.45))
  } else {
    data.frame(q = c(1, 0.8, 0.6), confound = c(1, 0.95, 0.9),
               leftover = c(1, 0.85, 0.6))
  }

  mock <- function(seu, depth_correct, depth_quantile = NULL, neighbors,
                   features, meta_name = "Entropy_mean", ...) {
    q <- if (isTRUE(depth_correct)) depth_quantile else 1
    row <- components[match(q, components$q), , drop = FALSE]
    score <- sqrt(row$confound) * x + sqrt(row$leftover) * other
    names(score) <- colnames(seu)
    calls$values[[length(calls$values) + 1L]] <- list(
      q = q, neighbors = neighbors, features = features,
      depth_correct = depth_correct
    )
    if (!isTRUE(depth_correct)) {
      return(data.frame(Entropy_mean = score, row.names = colnames(seu)))
    }
    totals <- Matrix::colSums(SeuratObject::GetAssayData(
      seu, assay = "RNA", layer = "counts"
    ))
    target <- floor(stats::quantile(totals, q, names = FALSE, type = 7))
    seu[[meta_name]] <- score
    seu[[paste0(meta_name, "_depth_corrected")]] <- score
    seu@misc$entropy_provenance <- list(
      target_depth = as.integer(target),
      n_cells_corrected = as.integer(sum(totals > target))
    )
    seu
  }
  list(mock = mock, calls = calls)
}


run_mock_depth_select <- function(profile = "correct", outdir = tempfile(),
                                  per_cluster = FALSE) {
  seu <- make_depth_select_fixture()
  mocked <- depth_select_mock(seu, profile)
  testthat::local_mocked_bindings(
    run_entropy = mocked$mock, .package = "Entropy", .env = parent.frame()
  )
  result <- entropy_depth_select(
    seu, group = "cluster", q_grid = c(0.8, 0.6), w = 1.5,
    outdir = outdir, per_cluster = per_cluster, seed = 42,
    assay = "RNA", k = 11, reduction = "pca_test", dims = 1:4,
    mode = "counts", parallelize = FALSE
  )
  list(result = result, calls = mocked$calls, seu = seu, outdir = outdir)
}


pdf_page_count <- function(path) {
  executable <- Sys.which("pdfinfo")
  if (!nzchar(executable)) return(NA_integer_)
  information <- system2(executable, path, stdout = TRUE, stderr = TRUE)
  as.integer(sub("^Pages:[[:space:]]+", "",
                 information[grepl("^Pages:", information)][[1L]]))
}


test_that("a depth-driven fixture chooses an interior quantile", {
  skip_if_not_installed("Seurat")
  run <- run_mock_depth_select("correct")
  result <- run$result

  expect_s3_class(result, "entropy_depth_selection")
  expect_true(result$decision$correct)
  expect_equal(result$decision$q, 0.8)
  expect_equal(result$frontier$q, c(1, 0.8, 0.6))
  expect_equal(result$frontier$objective[result$frontier$q == 1], 1.5)
  expect_lt(
    result$frontier$price_cumulative[result$frontier$q == result$decision$q],
    1 / result$decision$w
  )
  expect_true(all(file.exists(file.path(
    run$outdir,
    c("depth_selection.rds", "depth_selection_report.pdf", "frontier.csv",
      "covariates.csv", "within_group.csv")
  ))))
})


test_that("a fixture without an advantageous trade refuses correction", {
  skip_if_not_installed("Seurat")
  run <- run_mock_depth_select("refuse")

  expect_false(run$result$decision$correct)
  expect_equal(run$result$decision$q, 1)
  if (nzchar(Sys.which("pdfinfo"))) {
    expect_equal(pdf_page_count(file.path(
      run$outdir, "depth_selection_report.pdf"
    )), 2L)
  }
})


test_that("larger retention weights cannot choose a lower quantile", {
  frontier <- data.frame(
    q = c(1, 0.8, 0.6, 0.4),
    confound_removed = c(0, 0.3, 0.6, 0.8),
    retained = c(1, 0.9, 0.65, 0.35)
  )
  weights <- c(0.5, 1, 1.5, 2)
  choices <- vapply(weights, function(weight) {
    Entropy:::.entropy_depth_choose(frontier, weight)$q
  }, numeric(1L))
  expect_true(all(diff(choices) >= 0))
})


test_that("the original x and fixed graph are unchanged across the sweep", {
  skip_if_not_installed("Seurat")
  run <- run_mock_depth_select("correct")
  counts <- SeuratObject::GetAssayData(
    run$seu, assay = "RNA", layer = "counts"
  )
  original_x <- as.numeric(Matrix::colSums(counts > 0))
  expect_equal(run$result$cell_data$detected_genes, original_x)

  calls <- run$calls$values
  expect_equal(length(calls), 3L)
  expect_true(all(vapply(calls, function(call) {
    identical(call$neighbors, calls[[1L]]$neighbors)
  }, logical(1L))))
  expect_true(all(vapply(calls, function(call) {
    identical(call$features, calls[[1L]]$features)
  }, logical(1L))))
  baseline <- calls[[which(vapply(calls, function(call) call$q == 1,
                                  logical(1L)))]]
  expect_false(baseline$depth_correct)
})


test_that("restart recomputes only the missing quantile", {
  skip_if_not_installed("Seurat")
  outdir <- tempfile()
  seu <- make_depth_select_fixture()
  mocked <- depth_select_mock(seu, "correct")
  testthat::local_mocked_bindings(
    run_entropy = mocked$mock, .package = "Entropy"
  )
  arguments <- list(
    object = seu, group = "cluster", q_grid = c(0.8, 0.6), w = 1.5,
    outdir = outdir, per_cluster = FALSE, seed = 42, assay = "RNA",
    k = 11, reduction = "pca_test", dims = 1:4, mode = "counts",
    parallelize = FALSE
  )
  first <- do.call(entropy_depth_select, arguments)
  expect_length(mocked$calls$values, 3L)
  unlink(file.path(outdir, "score_q0800.rds"))
  mocked$calls$values <- list()
  second <- do.call(entropy_depth_select, arguments)

  expect_length(mocked$calls$values, 1L)
  expect_equal(mocked$calls$values[[1L]]$q, 0.8)
  expect_equal(second$frontier, first$frontier)
})


test_that("reports redraw without entropy calls and have fixed page counts", {
  skip_if_not_installed("Seurat")
  skip_if_not(nzchar(Sys.which("pdfinfo")))
  run <- run_mock_depth_select("correct", per_cluster = FALSE)
  expect_equal(pdf_page_count(file.path(
    run$outdir, "depth_selection_report.pdf"
  )), 3L)

  testthat::local_mocked_bindings(
    run_entropy = function(...) stop("entropy was recomputed"),
    .package = "Entropy"
  )
  redraw <- tempfile()
  expect_no_error(entropy_depth_report(
    readRDS(file.path(run$outdir, "depth_selection.rds")),
    outdir = redraw, per_cluster = TRUE
  ))
  expect_equal(pdf_page_count(file.path(
    redraw, "depth_selection_report.pdf"
  )), 4L)
})


test_that("the real run_entropy return contract supports a two-point sweep", {
  skip_on_cran()
  skip_if_not_installed("Seurat")
  data("pbmc_small", package = "SeuratObject", envir = environment())
  if (!"pca" %in% names(pbmc_small@reductions)) {
    pbmc_small <- suppressWarnings(Seurat::RunPCA(
      pbmc_small, npcs = 10, verbose = FALSE
    ))
  }
  outdir <- tempfile("depth-select-integration-")
  result <- suppressWarnings(entropy_depth_select(
    pbmc_small,
    group = "groups",
    q_grid = 0.8,
    outdir = outdir,
    per_cluster = FALSE,
    assay = "RNA",
    k = 11,
    reduction = "pca",
    dims = 1:5,
    features = rownames(pbmc_small)[1:3],
    mode = "counts",
    parallelize = FALSE
  ))

  counts <- SeuratObject::GetAssayData(
    pbmc_small, assay = "RNA", layer = "counts"
  )
  totals <- Matrix::colSums(counts)
  expected_target <- floor(stats::quantile(
    totals, 0.8, names = FALSE, type = 7
  ))
  corrected <- result$frontier[result$frontier$q == 0.8, , drop = FALSE]
  baseline_cache <- readRDS(file.path(outdir, "score_q1000.rds"))
  corrected_cache <- readRDS(file.path(outdir, "score_q0800.rds"))

  expect_s3_class(result, "entropy_depth_selection")
  expect_equal(result$frontier$q, c(1, 0.8))
  expect_equal(corrected$target_depth, expected_target)
  expect_equal(corrected$n_cells_corrected, sum(totals > expected_target))
  expect_named(baseline_cache$score, colnames(pbmc_small))
  expect_named(corrected_cache$score, colnames(pbmc_small))
  expect_true(all(is.finite(baseline_cache$score)))
  expect_true(all(is.finite(corrected_cache$score)))
})
