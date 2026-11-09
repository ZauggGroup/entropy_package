make_gene_blocking_fixture <- function(k = 12L) {
  n_cells <- 24L
  totals <- as.integer(round(seq(240L, 470L, length.out = n_cells)))
  counts <- vapply(totals, function(total) {
    c(total - 30L, 10L, 7L, 5L, 4L, 3L, 1L, 0L)
  }, integer(8L))
  dimnames(counts) <- list(
    paste0("gene", seq_len(nrow(counts))),
    paste0("cell", seq_len(ncol(counts)))
  )
  seu <- Seurat::CreateSeuratObject(
    counts = Matrix::Matrix(counts, sparse = TRUE), min.cells = 0
  )
  cells <- colnames(seu)
  neighbors <- lapply(seq_along(cells), function(index) {
    setdiff(cells, cells[[index]])[seq_len(k)]
  })
  names(neighbors) <- cells
  list(seu = seu, neighbors = neighbors)
}


run_gene_blocking_fixture <- function(fixture, ...) {
  arguments <- list(
    seu = fixture$seu, assay = "RNA", neighbors = fixture$neighbors,
    neighbors_method = "custom", k = length(fixture$neighbors[[1L]]),
    mode = "counts", add_assay = FALSE
  )
  arguments[names(list(...))] <- list(...)
  do.call(run_entropy, arguments)
}


test_that("gene-blocked means equal unblocked means", {
  skip_if_not_installed("Seurat")
  fixture <- make_gene_blocking_fixture()
  block_dir <- tempfile("entropy-blocks-")

  unblocked <- run_gene_blocking_fixture(fixture)
  blocked <- run_gene_blocking_fixture(
    fixture, block_dir = block_dir, block_genes = 3L
  )

  expect_equal(blocked[[1L]], unblocked[[1L]])
  expect_length(list.files(block_dir, pattern = "^block_.*[.]rds$"), 3L)
  expect_true(file.exists(file.path(block_dir, "manifest.rds")))
  expect_identical(
    attr(blocked, "entropy_provenance")$block_dir,
    normalizePath(block_dir, winslash = "/")
  )
})


test_that("gene-blocked depth-corrected means equal unblocked means", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("scuttle")
  fixture <- make_gene_blocking_fixture()
  block_dir <- tempfile("entropy-depth-blocks-")

  set.seed(410)
  unblocked <- suppressWarnings(run_gene_blocking_fixture(
    fixture, depth_correct = TRUE, target_depth = 240L
  ))
  set.seed(410)
  blocked <- suppressWarnings(run_gene_blocking_fixture(
    fixture, depth_correct = TRUE, target_depth = 240L,
    block_dir = block_dir, block_genes = 3L
  ))

  expect_equal(blocked$entropy_mean_RNA, unblocked$entropy_mean_RNA)
  expect_equal(
    blocked$entropy_mean_RNA_depth_corrected,
    unblocked$entropy_mean_RNA_depth_corrected
  )
  expect_identical(
    blocked@misc$entropy_provenance$block_dir,
    normalizePath(block_dir, winslash = "/")
  )
})


test_that("a blocked run resumes only its missing block", {
  skip_if_not_installed("Seurat")
  fixture <- make_gene_blocking_fixture()
  block_dir <- tempfile("entropy-resume-")
  first <- run_gene_blocking_fixture(
    fixture, block_dir = block_dir, block_genes = 3L
  )
  files <- sort(list.files(
    block_dir, pattern = "^block_.*[.]rds$", full.names = TRUE
  ))
  old_time <- as.POSIXct("2020-01-01 00:00:00", tz = "UTC")
  Sys.setFileTime(files[c(1L, 3L)], old_time)
  unlink(files[[2L]])

  resumed <- run_gene_blocking_fixture(
    fixture, block_dir = block_dir, block_genes = 3L
  )

  expect_equal(resumed[[1L]], first[[1L]])
  expect_equal(as.numeric(file.info(files[c(1L, 3L)])$mtime),
               rep(as.numeric(old_time), 2L))
  expect_true(file.info(files[[2L]])$mtime > old_time)
})


test_that("a changed call cannot reuse a block manifest", {
  skip_if_not_installed("Seurat")
  fixture_12 <- make_gene_blocking_fixture(k = 12L)
  block_dir <- tempfile("entropy-mismatch-")
  run_gene_blocking_fixture(
    fixture_12, block_dir = block_dir, block_genes = 3L
  )
  fixture_11 <- make_gene_blocking_fixture(k = 11L)

  expect_error(
    run_gene_blocking_fixture(
      fixture_11, block_dir = block_dir, block_genes = 3L
    ),
    "manifest mismatch for field 'k'"
  )
})


test_that("blocking rejects operations that need the full matrix", {
  skip_if_not_installed("Seurat")
  fixture <- make_gene_blocking_fixture()
  block_dir <- tempfile("entropy-invalid-blocks-")

  expect_error(
    run_gene_blocking_fixture(
      fixture, block_dir = block_dir, block_genes = 3L,
      add_assay = TRUE
    ),
    "add_assay = TRUE"
  )
  expect_error(
    run_gene_blocking_fixture(
      fixture, block_dir = block_dir, block_genes = 3L,
      output_path = tempfile(fileext = ".rds")
    ),
    "output_path.*block_dir"
  )
})


test_that("large output warning uses the two GB threshold", {
  expect_warning(
    .entropy_warn_large_output(
      n_cells = 100000L, n_genes = 3000L,
      depth_correct = FALSE, block_dir = NULL
    ),
    "above the 2 GB.*set `block_dir`"
  )
  expect_no_warning(.entropy_warn_large_output(
    n_cells = 1000L, n_genes = 1000L,
    depth_correct = TRUE, block_dir = NULL
  ))
  expect_no_warning(.entropy_warn_large_output(
    n_cells = 100000L, n_genes = 3000L,
    depth_correct = FALSE, block_dir = tempfile()
  ))
})


test_that("entropy_read_genes reads only requested rows in requested order", {
  skip_if_not_installed("Seurat")
  fixture <- make_gene_blocking_fixture()
  block_dir <- tempfile("entropy-read-")
  full_path <- tempfile(fileext = ".rds")
  run_gene_blocking_fixture(fixture, output_path = full_path)
  run_gene_blocking_fixture(
    fixture, block_dir = block_dir, block_genes = 3L
  )
  full <- readRDS(full_path)
  requested <- c("gene7", "gene2", "gene7", "gene4")

  expect_equal(
    entropy_read_genes(block_dir, requested),
    full[requested, , drop = FALSE]
  )
  expect_error(
    entropy_read_genes(block_dir, c("gene2", "not_a_gene")),
    "not_a_gene"
  )
})


test_that("entropy_read_genes rejects missing or invalid block files", {
  skip_if_not_installed("Seurat")
  fixture <- make_gene_blocking_fixture()
  block_dir <- tempfile("entropy-read-errors-")
  run_gene_blocking_fixture(
    fixture, block_dir = block_dir, block_genes = 3L
  )
  second <- file.path(block_dir, "block_0002.rds")
  saved <- readRDS(second)
  unlink(second)
  expect_error(entropy_read_genes(block_dir, "gene4"), "missing")

  saveRDS(saved, second)
  first <- file.path(block_dir, "block_0001.rds")
  corrupt <- readRDS(first)
  corrupt$genes[[1L]] <- "wrong_gene"
  saveRDS(corrupt, first)
  expect_error(entropy_read_genes(block_dir, "gene1"), "invalid")
})
