# Gene-blocked entropy output -------------------------------------------------

.entropy_block_threshold_bytes <- 2 * 1024^3
.entropy_block_target_bytes <- 1024^3


#' @keywords internal
.entropy_output_bytes <- function(n_cells, n_genes, depth_correct) {
  as.double(n_cells) * as.double(n_genes) * 8 *
    if (isTRUE(depth_correct)) 2 else 1
}


#' @keywords internal
.entropy_format_gb <- function(bytes) {
  sprintf("%.2f GB", bytes / 1024^3)
}


#' @keywords internal
.entropy_warn_large_output <- function(n_cells, n_genes, depth_correct,
                                       block_dir) {
  bytes <- .entropy_output_bytes(n_cells, n_genes, depth_correct)
  if (is.null(block_dir) && bytes > .entropy_block_threshold_bytes) {
    warning(sprintf(
      paste0(
        "The estimated entropy output is %s (%d cells x %d genes%s), ",
        "above the 2 GB in-memory threshold. The calculation will continue ",
        "in memory; set `block_dir` to score genes in restartable blocks."
      ),
      .entropy_format_gb(bytes), n_cells, n_genes,
      if (isTRUE(depth_correct)) ", two entropy passes" else ""
    ), call. = FALSE)
  }
  invisible(bytes)
}


#' @keywords internal
.entropy_save_rds_atomic <- function(object, path) {
  temporary <- tempfile(
    pattern = paste0(basename(path), "."),
    tmpdir = dirname(path)
  )
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(object, temporary)
  if (!file.rename(temporary, path)) {
    stop(sprintf("Could not finalize entropy block file '%s'.", path),
         call. = FALSE)
  }
  invisible(path)
}


#' @keywords internal
.entropy_manifest_value <- function(value) {
  if (is.null(value)) "NULL" else paste(utils::head(value, 4L), collapse = ", ")
}


#' @keywords internal
.entropy_check_manifest <- function(stored, expected, path) {
  fields <- union(names(expected), names(stored))
  for (field in fields) {
    if (!field %in% names(stored) || !field %in% names(expected) ||
        !identical(stored[[field]], expected[[field]])) {
      stop(sprintf(
        paste0(
          "Block manifest mismatch for field '%s' (stored: %s; current: %s). ",
          "Use a different `block_dir`, or remove the old directory after ",
          "confirming it is no longer needed: %s"
        ),
        field,
        .entropy_manifest_value(stored[[field]]),
        .entropy_manifest_value(expected[[field]]),
        dirname(path)
      ), call. = FALSE)
    }
  }
  invisible(TRUE)
}


#' @keywords internal
.entropy_prepare_manifest <- function(block_dir, cells, genes, block_genes,
                                      k, assay, mode, reduction, dims,
                                      restrict_by, depth_correct, target_depth) {
  indices <- split(
    seq_along(genes),
    ceiling(seq_along(genes) / block_genes)
  )
  boundaries <- data.frame(
    block = seq_along(indices),
    start = vapply(indices, min, integer(1L)),
    end = vapply(indices, max, integer(1L)),
    file = sprintf("block_%04d.rds", seq_along(indices)),
    stringsAsFactors = FALSE
  )
  manifest <- list(
    cells = cells,
    genes = genes,
    blocks = boundaries,
    k = k,
    assay = assay,
    mode = mode,
    reduction = reduction,
    dims = dims,
    restrict_by = restrict_by,
    depth_correct = depth_correct,
    target_depth = target_depth,
    package_version = as.character(utils::packageVersion("Entropy"))
  )
  path <- file.path(block_dir, "manifest.rds")
  if (file.exists(path)) {
    stored <- tryCatch(
      readRDS(path),
      error = function(error) stop(sprintf(
        "Cannot read the existing block manifest '%s': %s",
        path, conditionMessage(error)
      ), call. = FALSE)
    )
    .entropy_check_manifest(stored, manifest, path)
  } else {
    existing_blocks <- list.files(
      block_dir, pattern = "^block_[0-9]+[.]rds$", full.names = TRUE
    )
    if (length(existing_blocks)) {
      stop(sprintf(
        paste0(
          "Block files exist but their manifest is missing in '%s'. ",
          "Use a different `block_dir`, or remove the orphaned block files ",
          "after confirming they are no longer needed."
        ),
        block_dir
      ), call. = FALSE)
    }
    .entropy_save_rds_atomic(manifest, path)
  }
  list(manifest = manifest, indices = indices)
}


#' @keywords internal
.entropy_validate_block <- function(block, manifest, block_index, path) {
  row <- manifest$blocks[block_index, , drop = FALSE]
  expected_genes <- manifest$genes[seq.int(row$start, row$end)]
  problem <- NULL

  if (!is.list(block) ||
      !identical(block$block, as.integer(block_index)) ||
      !identical(block$genes, expected_genes) ||
      !identical(block$cells, manifest$cells)) {
    problem <- "its block number, genes, or cells do not match the manifest"
  } else if (!is.matrix(block$entropy) ||
             !identical(dimnames(block$entropy),
                        list(expected_genes, manifest$cells))) {
    problem <- "its uncorrected entropy matrix does not match the manifest"
  } else if (any(!is.finite(block$entropy))) {
    problem <- "its uncorrected entropy matrix contains non-finite values"
  } else if (isTRUE(manifest$depth_correct) &&
             (!is.matrix(block$entropy_depth_corrected) ||
              !identical(dimnames(block$entropy_depth_corrected),
                         list(expected_genes, manifest$cells)))) {
    problem <- "its depth-corrected entropy matrix does not match the manifest"
  } else if (isTRUE(manifest$depth_correct) &&
             any(!is.finite(block$entropy_depth_corrected))) {
    problem <- "its depth-corrected entropy matrix contains non-finite values"
  } else if (!isTRUE(manifest$depth_correct) &&
             !is.null(block$entropy_depth_corrected)) {
    problem <- "it contains a depth-corrected matrix for an uncorrected run"
  }

  if (!is.null(problem)) {
    stop(sprintf(
      "Entropy block file '%s' is invalid: %s. Recompute it in a new `block_dir`.",
      path, problem
    ), call. = FALSE)
  }
  invisible(TRUE)
}


#' Read selected genes from blocked entropy output
#'
#' Reads gene-level entropy written by [run_entropy()] with `block_dir` set.
#' Only block files containing requested genes are loaded. For a run with
#' `depth_correct = TRUE`, the returned values are the depth-corrected entropy
#' values, matching the matrix that `run_entropy(add_assay = TRUE)` would store.
#'
#' @param block_dir Directory containing `manifest.rds` and the entropy block
#'   files written by [run_entropy()].
#' @param genes Character vector of gene names to read. Rows are returned in
#'   this order; repeated names are allowed.
#'
#' @return A numeric matrix with requested genes in rows and the cells from the
#'   blocked run in columns.
#' @export
entropy_read_genes <- function(block_dir, genes) {
  if (length(block_dir) != 1L || !is.character(block_dir) ||
      is.na(block_dir) || !nzchar(block_dir)) {
    stop("`block_dir` must be one non-empty directory path.", call. = FALSE)
  }
  if (!is.character(genes) || length(genes) < 1L || anyNA(genes)) {
    stop("`genes` must be a non-empty character vector without NA values.",
         call. = FALSE)
  }
  manifest_path <- file.path(block_dir, "manifest.rds")
  if (!file.exists(manifest_path)) {
    stop(sprintf("Block manifest is missing: '%s'.", manifest_path),
         call. = FALSE)
  }
  manifest <- tryCatch(
    readRDS(manifest_path),
    error = function(error) stop(sprintf(
      "Cannot read block manifest '%s': %s",
      manifest_path, conditionMessage(error)
    ), call. = FALSE)
  )
  required <- c(
    "cells", "genes", "blocks", "depth_correct", "package_version"
  )
  required_block_columns <- c("block", "start", "end", "file")
  if (!is.list(manifest) || !all(required %in% names(manifest)) ||
      !is.data.frame(manifest$blocks) ||
      !all(required_block_columns %in% names(manifest$blocks)) ||
      nrow(manifest$blocks) < 1L) {
    stop(sprintf("Block manifest '%s' is invalid.", manifest_path),
         call. = FALSE)
  }
  positions <- match(genes, manifest$genes)
  if (anyNA(positions)) {
    first <- genes[which(is.na(positions))[1L]]
    stop(sprintf("Requested gene '%s' is not in the block manifest.", first),
         call. = FALSE)
  }

  block_ids <- findInterval(positions, manifest$blocks$end) + 1L
  # A position equal to an end belongs to that block, not the next one.
  exact_end <- match(positions, manifest$blocks$end, nomatch = 0L)
  block_ids[exact_end > 0L] <- exact_end[exact_end > 0L]
  result <- matrix(
    NA_real_, nrow = length(genes), ncol = length(manifest$cells),
    dimnames = list(genes, manifest$cells)
  )

  for (block_index in unique(block_ids)) {
    block_path <- file.path(
      block_dir, manifest$blocks$file[[block_index]]
    )
    if (!file.exists(block_path)) {
      stop(sprintf("Entropy block file is missing: '%s'.", block_path),
           call. = FALSE)
    }
    block <- tryCatch(
      readRDS(block_path),
      error = function(error) stop(sprintf(
        "Cannot read entropy block file '%s': %s",
        block_path, conditionMessage(error)
      ), call. = FALSE)
    )
    .entropy_validate_block(block, manifest, block_index, block_path)
    matrix_name <- if (isTRUE(manifest$depth_correct)) {
      "entropy_depth_corrected"
    } else {
      "entropy"
    }
    requested <- which(block_ids == block_index)
    result[requested, ] <- block[[matrix_name]][genes[requested], , drop = FALSE]
    rm(block)
  }
  result
}


#' @keywords internal
.entropy_compute_matrix <- function(mat, mat_bin, mode, cell_names, neighbors,
                                    chunk_method, chunks, parallelize,
                                    n_cores) {
  n_cells <- length(cell_names)
  pb <- progress::progress_bar$new(
    format = "  Processing [:bar] :percent | :current/:total cells",
    total = n_cells, clear = FALSE, width = 60
  )
  entropy_matrix <- matrix(
    0, nrow = nrow(mat), ncol = n_cells,
    dimnames = list(rownames(mat), cell_names)
  )

  if (chunk_method == "none") {
    for (i in seq_along(cell_names)) {
      cell_name <- cell_names[[i]]
      if (i %% 100L == 0L) pb$tick(100L)
      nn_idx <- get_neighbors(cell_name, NULL, neighbors)
      if (mode == "binary") {
        entropy_matrix[, cell_name] <- calculate_entropy_for_cell_binary_fast(
          cell_name = cell_name, mat_bin = mat_bin, neighbors_idx = nn_idx
        )
      } else {
        entropy_matrix[, cell_name] <- calculate_entropy_for_cell(
          cell_name = cell_name, mat = mat, neighbors_idx = nn_idx
        )
      }
    }
    if (n_cells %% 100L != 0L) pb$tick(n_cells %% 100L)
    return(entropy_matrix)
  }

  run_chunk <- function(chunk_cell_names) {
    all_nn_idx <- unique(unlist(lapply(chunk_cell_names, function(cell) {
      get_neighbors(cell, NULL, neighbors)
    })))
    mat_chunk <- if (mode == "binary") {
      mat_bin[, all_nn_idx, drop = FALSE]
    } else {
      mat[, all_nn_idx, drop = FALSE]
    }
    idx_remap <- stats::setNames(seq_along(all_nn_idx), all_nn_idx)
    chunk_matrix <- matrix(
      0, nrow = nrow(mat), ncol = length(chunk_cell_names),
      dimnames = list(rownames(mat), chunk_cell_names)
    )
    for (cell_name in chunk_cell_names) {
      nn_idx_full <- get_neighbors(cell_name, NULL, neighbors)
      nn_idx_chunk <- as.integer(idx_remap[as.character(nn_idx_full)])
      if (mode == "binary") {
        chunk_matrix[, cell_name] <- calculate_entropy_for_cell_binary_fast(
          cell_name = cell_name, mat_bin = mat_chunk,
          neighbors_idx = nn_idx_chunk
        )
      } else {
        chunk_matrix[, cell_name] <- calculate_entropy_for_cell(
          cell_name = cell_name, mat = mat_chunk,
          neighbors_idx = nn_idx_chunk
        )
      }
    }
    chunk_matrix
  }

  if (parallelize) {
    futile.logger::flog.info("Running chunks in parallel with %d cores.", n_cores)
    chunk_results <- parallel::mclapply(chunks, run_chunk, mc.cores = n_cores)
  } else {
    futile.logger::flog.info("Running chunks sequentially.")
    chunk_results <- lapply(chunks, function(chunk_cell_names) {
      pb$tick(length(chunk_cell_names))
      run_chunk(chunk_cell_names)
    })
  }
  for (i in seq_along(chunk_results)) {
    names_i <- colnames(chunk_results[[i]])
    entropy_matrix[, names_i] <- chunk_results[[i]]
  }
  entropy_matrix
}
