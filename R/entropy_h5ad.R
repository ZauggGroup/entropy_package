#' Estimate gene expression entropy directly from an AnnData (.h5ad) file
#'
#' Seurat-free entry point for scanpy users. Reads the counts matrix and a PCA
#' embedding straight from an `.h5ad` file, computes neighbors internally (or
#' accepts a user-supplied neighbor list), and returns per-cell mean entropy.
#' The heavy lifting reuses the same engine as `run_entropy()`
#' (`calculate_entropy_for_cell*`, the BUB precomputed table), so results match.
#'
#' The default return is a data frame of per-cell mean entropy keyed by cell id,
#' meant to be joined back into `adata.obs` in Python:
#'   `adata.obs['Entropy_mean'] = pd.read_csv('out.csv',
#'   index_col=0)['Entropy_mean']`
#' The full gene x cell entropy matrix is opt-in via `output_path`.
#'
#' @details
#' Direct h5ad support is optional and requires R >= 4.5.0,
#' `anndataR` >= 1.0.0, and `rhdf5`. The other Entropy functions remain
#' available on older R versions and do not require these packages.
#'
#' @param path Character. Path to the `.h5ad` file.
#' @param layer Character or NULL. Name of the AnnData layer to score. `NULL`
#'   (default) uses `X`. For `mode = "counts"` this must hold integer counts.
#' @param neighbors Optional named list for `neighbors_method = "custom"`, where
#'   `neighbors[[cell_name]]` returns neighbor cell names or integer indices.
#' @param neighbors_method Character. One of `"auto"` or `"custom"`. `"auto"`
#'   computes neighbors from the PCA embedding with `RANN::nn2`. `"custom"` uses
#'   the user-provided `neighbors` list. (There is no `"seurat"` here: no object.)
#' @param features Optional character vector of feature names to restrict scoring
#'   to. `NULL` (default) scores all features.
#' @param restrict_by Optional character. Name of a column in `obs` to restrict
#'   neighbors within. Same semantics as `run_entropy()`: auto recomputes kNN
#'   within each group; groups smaller than `k + 1` return entropy 0; NA cells
#'   return entropy 0.
#' @param pca_key Character. Key in `obsm` holding the embedding. Default
#'   `"X_pca"`.
#' @param dims Integer vector or single integer. Columns of the embedding to use.
#'   A single `n` means `1:n`. A vector is used as given (e.g. `2:20` skips
#'   component 1 for LSI/ATAC). Default 20.
#' @param k Integer. Number of neighbors when computing internally. Must be at
#'   least 11 (counts estimator floor). Default 15.
#' @param mode Character. `"counts"` (BUB) or `"binary"` (Bernoulli Shannon).
#' @param chunk_method Character. `"none"` (default) or `"auto"`. `"auto"` splits
#'   cells into compact chunks by k-means on the embedding to cut sparse-subset
#'   cost. Recommended above ~10,000 cells.
#' @param chunk_cells Integer. Max cells per chunk when chunking. Default 2000.
#' @param parallelize Logical. Parallelize chunks across cores. Chunked mode
#'   only. Not supported on Windows. Default FALSE.
#' @param n_cores Integer. Cores when `parallelize = TRUE`.
#' @param depth_correct Logical. If `TRUE`, also compute entropy after thinning
#'   each cell's counts to a common library depth. Sequencing depth inflates
#'   entropy for a technical reason: entropy is concave, so subsampling biases
#'   its estimate downward and more deeply sequenced cells score higher. The
#'   returned data frame then carries `Entropy_mean_depth_corrected` alongside
#'   `Entropy_mean`, plus a logical `depth_corrected` marking which cells were
#'   actually thinned. Neighbourhoods are unchanged, only the counts fed to the
#'   estimator are. When `output_path` is also given, the saved matrix is the
#'   depth-corrected one. The default is `FALSE`.
#' @param depth_quantile Numeric. Quantile of the cell count totals used to
#'   choose the target depth when `target_depth` is `NULL`. Must be strictly
#'   between 0 and 1. It is effectively the fraction of cells left uncorrected,
#'   because cells already below the target keep their counts. The default is
#'   0.05.
#' @param target_depth Optional positive integer target for count thinning. If
#'   supplied, it overrides `depth_quantile`.
#' @param return_inputs Logical. If `TRUE`, return a list with the per-cell
#'   `entropy` vector, the `embedding` used for neighbours, the `obs` metadata
#'   read from the file, and the usual per-cell `results` data frame. This is
#'   the input [entropy_diagnose()] needs, so the diagnostic can be run without
#'   reading the file a second time and without a Seurat object. The default is
#'   `FALSE`, which returns only the data frame.
#' @param output_path Optional character. If set, saves the full gene x cell
#'   entropy matrix with `saveRDS()`. The per-cell mean data frame is still
#'   returned.
#'
#' @return A data frame with cell identifiers as row names and one column,
#'   `Entropy_mean`, containing the mean entropy across the scored genes.
#' @export
#'
#' @examples
#' \dontrun{
#' run_entropy_h5ad("adata.h5ad", layer = "counts", k = 15, mode = "counts")
#' }
run_entropy_h5ad <- function(
    path,
    layer            = NULL,
    neighbors        = NULL,
    neighbors_method = c("auto", "custom"),
    features         = NULL,
    restrict_by      = NULL,
    pca_key          = "X_pca",
    dims             = 20,
    k                = 15,
    mode             = c("counts", "binary"),
    chunk_method     = c("none", "auto"),
    chunk_cells      = 2000,
    parallelize      = FALSE,
    n_cores          = parallel::detectCores() - 1L,
    output_path      = NULL,
    return_inputs    = FALSE,
    depth_correct    = FALSE,
    depth_quantile   = 0.05,
    target_depth     = NULL
){

  neighbors_method <- match.arg(neighbors_method)
  mode <- match.arg(mode)
  chunk_method <- match.arg(chunk_method)

  if (length(return_inputs) != 1L || !is.logical(return_inputs) ||
      is.na(return_inputs)) {
    stop("`return_inputs` must be TRUE or FALSE.", call. = FALSE)
  }
  if (length(depth_correct) != 1L || !is.logical(depth_correct) ||
      is.na(depth_correct)) {
    stop("`depth_correct` must be TRUE or FALSE.", call. = FALSE)
  }
  if (length(depth_quantile) != 1L || !is.numeric(depth_quantile) ||
      !is.finite(depth_quantile) || depth_quantile <= 0 ||
      depth_quantile >= 1) {
    stop("`depth_quantile` must be one number strictly between 0 and 1.",
         call. = FALSE)
  }
  if (!is.null(target_depth) &&
      (length(target_depth) != 1L || !is.numeric(target_depth) ||
       !is.finite(target_depth) || target_depth != floor(target_depth) ||
       target_depth < 1)) {
    stop("`target_depth` must be one positive integer or `NULL`.", call. = FALSE)
  }

  check_h5ad_dependencies()

  if (!file.exists(path)) {
    stop(sprintf("File not found: %s", path))
  }

  ### check K
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) ||
      k != as.integer(k)) {
    stop("k must be a single integer.")
  }
  k <- as.integer(k)
  if (k < 11L) {
    stop(sprintf("k must be at least 11 (got %d). The counts entropy estimator is undefined below 11 neighbors.", k))
  }

  ### read the h5ad: counts (genes x cells), embedding (cells x dims), obs
  inputs <- read_h5ad_inputs(path, layer = layer, pca_key = pca_key,
                             restrict_by = restrict_by,
                             need_embedding = (neighbors_method == "auto" ||
                                                 chunk_method == "auto" ||
                                                 return_inputs),
                             read_obs = return_inputs)
  mat        <- inputs$mat          # dgCMatrix, genes x cells
  emb_full   <- inputs$emb          # cells x dims, or NULL
  grp_vec    <- inputs$groups       # named char vector or NULL
  cell_names <- inputs$cell_names
  n_cells    <- length(cell_names)

  validate_h5ad_matrix(mat, mode = mode, layer = layer)

  futile.logger::flog.info(
    "Loaded '%s': %d cells and %d features.", path, n_cells, nrow(mat)
  )

  ### check features
  if (!is.null(features)) {
    if (!is.character(features)) {
      stop("`features` must be a character vector of feature names, not ",
           class(features)[1], ". Indices are not accepted.")
    }
    features <- unique(features)
    missing <- setdiff(features, rownames(mat))
    if (length(missing) == length(features)) {
      stop(sprintf("None of the %d requested features are present. Check your feature names.", length(features)))
    }
    if (length(missing) > 0) {
      stop(sprintf(
        "%d of %d requested features are not present: %s%s",
        length(missing), length(features),
        paste(utils::head(missing, 10), collapse = ", "),
        if (length(missing) > 10) ", ..." else ""
      ))
    }
    futile.logger::flog.info("Restricting entropy computation to %d features.", length(features))
  }

  ### check parallelization
  if (parallelize && chunk_method == "none") {
    futile.logger::flog.warn("parallelize = TRUE has no effect when chunk_method = 'none'.")
  }
  if (parallelize && .Platform$OS.type == "windows") {
    futile.logger::flog.warn("parallelize = TRUE is not supported on Windows. Falling back to sequential.")
    parallelize <- FALSE
  }
  if (parallelize && n_cores > parallel::detectCores() - 1L) {
    futile.logger::flog.warn("n_cores = %d exceeds available cores. Setting to %d.",
                             n_cores, parallel::detectCores() - 1L)
    n_cores <- parallel::detectCores() - 1L
  }
  if (parallelize && n_cores < 1L) {
    stop("n_cores must be at least 1.")
  }

  ### chunk check
  if (chunk_method == "none" && n_cells > 50000) {
    futile.logger::flog.warn(
      "Dataset has %d cells but chunk_method = 'none'. This may be very slow. Consider chunk_method = 'auto'.",
      n_cells
    )
  }
  if (chunk_method == "auto") {
    if (!is.numeric(chunk_cells) || length(chunk_cells) != 1 || chunk_cells != as.integer(chunk_cells)) {
      stop("chunk_cells must be a single integer.")
    }
    if (chunk_cells < 500) {
      stop(sprintf("chunk_cells = %d is too low. Minimum recommended value is 500.", chunk_cells))
    }
    if (chunk_cells > 10000) {
      futile.logger::flog.warn(
        "chunk_cells = %d is very large. Values between 1000 and 5000 are recommended.", chunk_cells
      )
    }
    if (n_cells < chunk_cells * 2) {
      futile.logger::flog.warn(
        "Dataset has %d cells, too small to benefit from chunking with chunk_cells = %d. Falling back to 'none'.",
        n_cells, chunk_cells
      )
      chunk_method <- "none"
    }
  }

  ### restrict_by check
  if (!is.null(restrict_by)) {
    if (!is.character(restrict_by) || length(restrict_by) != 1) {
      stop("`restrict_by` must be a single character string naming an obs column.")
    }
    # presence already validated inside read_h5ad_inputs; here we just warn on NA
    n_na <- sum(is.na(grp_vec))
    if (n_na > 0) {
      futile.logger::flog.warn(
        "Column '%s' has %d cells with NA. These get no within-group neighbors (entropy 0).",
        restrict_by, n_na
      )
    }
    futile.logger::flog.info("Neighbors will be restricted within groups of '%s'.", restrict_by)
  }

  ### resolve dims against the embedding (only if we need it)
  need_embedding <- (neighbors_method == "auto") || (chunk_method == "auto")
  if (need_embedding) {
    if (is.null(emb_full)) {
      # DECISION POINT: we assume X_pca exists. If you want to support files
      # with no embedding, this is where PCA (and an nfeatures arg) would go.
      stop(sprintf(
        "obsm['%s'] not found in the file. Compute a PCA in scanpy first (sc.pp.pca), or pass a precomputed embedding.",
        pca_key
      ))
    }
    if (length(dims) == 1) {
      if (is.na(dims) || dims < 1) stop("dims must be a positive integer.")
      dims <- seq_len(dims)
    }
    dims <- as.integer(dims)
    if (any(is.na(dims)) || any(dims < 1)) stop("dims must be positive integers.")
    n_avail <- ncol(emb_full)
    if (max(dims) > n_avail) {
      futile.logger::flog.warn("Requested dims up to %d but embedding has %d. Capping.", max(dims), n_avail)
      dims <- dims[dims <= n_avail]
    }
    if (length(dims) < 1) stop("No valid dims after capping against the embedding.")
    emb <- emb_full[, dims, drop = FALSE]
    if (!is.numeric(emb) || any(!is.finite(emb))) {
      stop(sprintf("obsm['%s'] must contain only finite numeric values.", pca_key))
    }
  }

  ###################### NN creation -> normalize to nn_list (cell -> int idx) #

  restrict_method <- neighbors_method
  nn_list <- vector("list", n_cells)
  names(nn_list) <- cell_names

  if (neighbors_method == "auto") {

    if (is.null(restrict_by)) {
      if (n_cells < k + 1L) {
        stop(sprintf(
          paste0(
            "Cannot compute %d neighbors from only %d cells: global automatic ",
            "neighboring requires at least k + 1 = %d cells."
          ),
          k, n_cells, k + 1L
        ))
      }
      # unrestricted: one global nn2 (replaces Seurat::FindNeighbors)
      futile.logger::flog.info("Computing global neighbors with nn2 (k = %d).", k)
      nn <- RANN::nn2(emb, k = k + 1L)$nn.idx
      for (r in seq_len(n_cells)) {
        row <- nn[r, nn[r, ] != r]
        nn_list[[r]] <- utils::head(row, k)
      }

    } else {
      # restricted: recompute kNN within each group (lifted from run_entropy)
      futile.logger::flog.info("Recomputing neighbors within each group of '%s'.", restrict_by)
      grp <- grp_vec
      groups <- split(cell_names, grp)  # NA cells dropped by split()

      for (g in names(groups)) {
        g_cells <- groups[[g]]
        g_n     <- length(g_cells)
        if (g_n < k + 1L) {
          futile.logger::flog.warn(
            "Group '%s' has %d cells (< k+1 = %d); too small. These cells get entropy 0.",
            g, g_n, k + 1L
          )
          for (cell in g_cells) nn_list[[cell]] <- integer(0)
          next
        }
        nn_g       <- RANN::nn2(emb[g_cells, , drop = FALSE], k = k + 1L)$nn.idx
        global_idx <- match(g_cells, cell_names)
        for (r in seq_along(g_cells)) {
          nbr_row    <- nn_g[r, nn_g[r, ] != r]
          local_nbrs <- utils::head(nbr_row, k)
          nn_list[[g_cells[r]]] <- global_idx[local_nbrs]
        }
      }
      na_cells <- cell_names[is.na(grp)]
      for (cell in na_cells) nn_list[[cell]] <- integer(0)

      futile.logger::flog.info(
        "Within-group neighbors: mean %.1f per cell; %d cells with 0 neighbors.",
        mean(lengths(nn_list)), sum(lengths(nn_list) == 0)
      )
    }

  } else if (neighbors_method == "custom") {
    if (is.null(neighbors)) {
      stop("`neighbors` must be provided when neighbors_method = 'custom'.")
    }
    if (!is.list(neighbors) || is.null(names(neighbors))) {
      stop("`neighbors` must be a named list where neighbors[[cell_name]] returns neighbor names or indices.")
    }
    missing_cells <- setdiff(cell_names, names(neighbors))
    if (length(missing_cells) > 0) {
      stop(
        "The custom neighbors list is missing entries for ", length(missing_cells), " cells.\n",
        "Each cell must appear in names(neighbors), even if it has no neighbors."
      )
    }
    futile.logger::flog.info("Using user-provided custom neighbors.")
    # resolve names/indices to integer indices into cell_names, once
    for (cell in cell_names) {
      v <- neighbors[[cell]]
      if (length(v) == 0) {
        nn_list[[cell]] <- integer(0)
      } else if (is.character(v)) {
        idx <- match(v, cell_names)
        if (any(is.na(idx))) {
          stop(sprintf("Some neighbor names for cell '%s' were not found among cells.", cell))
        }
        nn_list[[cell]] <- idx
      } else {
        idx <- as.integer(v)
        if (any(is.na(idx)) || any(idx < 1) || any(idx > n_cells)) {
          stop(sprintf("Some neighbor indices for cell '%s' are out of range.", cell))
        }
        nn_list[[cell]] <- idx
      }
    }
  }

  #### Filter within groups (custom only; auto already handled itself)
  if (!is.null(restrict_by) && restrict_method == "custom") {
    futile.logger::flog.info("Filtering provided neighbors within groups of '%s'.", restrict_by)
    grp <- grp_vec
    n_before <- lengths(nn_list)
    for (cell in cell_names) {
      nn_idx   <- nn_list[[cell]]
      cell_grp <- grp[cell]
      if (is.na(cell_grp) || length(nn_idx) == 0) {
        nn_list[[cell]] <- integer(0)
      } else {
        nn_list[[cell]] <- nn_idx[grp[nn_idx] == cell_grp & !is.na(grp[nn_idx])]
      }
    }
    futile.logger::flog.info(
      "Within-group filtering: mean neighbors %.1f -> %.1f; %d cells now with 0 neighbors.",
      mean(n_before), mean(lengths(nn_list)), sum(lengths(nn_list) == 0)
    )
  }

  #### Orphan cells with < 11 neighbors (custom only; counts estimator floor)
  # NOTE: mirrors run_entropy, which applies this in both modes for custom.
  # Binary mode does not strictly need 11, so this is slightly conservative.
  if (restrict_method == "custom") {
    n_orphaned <- 0L
    for (cell in cell_names) {
      nn_idx <- nn_list[[cell]]
      if (length(nn_idx) < 11L) {
        if (length(nn_idx) > 0L) n_orphaned <- n_orphaned + 1L
        nn_list[[cell]] <- integer(0)
      }
    }
    if (n_orphaned > 0L) {
      futile.logger::flog.warn(
        "%d of %d cells have fewer than 11 neighbors and will return entropy 0.",
        n_orphaned, n_cells
      )
    }
  }

  ###################### NN finish ############################################

  depth_target <- NULL
  depth_corrected_mat <- NULL
  depth_corrected_cells <- NULL
  n_cells_below_target <- 0L
  n_cells_corrected <- 0L
  if (depth_correct) {
    depth_totals <- Matrix::colSums(mat)
    depth_target <- if (is.null(target_depth)) {
      floor(stats::quantile(
        depth_totals, probs = depth_quantile, names = FALSE, type = 7
      ))
    } else {
      as.integer(target_depth)
    }
    if (!is.finite(depth_target) || depth_target < 1) {
      stop("The depth correction target is below 1 count.", call. = FALSE)
    }
    n_cells_below_target <- sum(depth_totals < depth_target)
    n_cells_corrected <- sum(depth_totals > depth_target)
    depth_corrected_cells <- depth_totals > depth_target
    names(depth_corrected_cells) <- cell_names

    below_fraction <- n_cells_below_target / length(depth_totals)
    if (below_fraction > 0.1) {
      warning(sprintf(
        paste0(
          "%d cells (%.1f%%) are below the target depth and cannot be ",
          "downsampled; they will keep their original counts."
        ),
        n_cells_below_target, 100 * below_fraction
      ), call. = FALSE)
    }
    if (depth_target < 200) {
      warning(sprintf(
        paste0(
          "The target depth is %d counts. Entropy estimates at depths below ",
          "200 counts may be unreliable."
        ),
        depth_target
      ), call. = FALSE)
    }
    if (mode == "binary") {
      warning(
        paste0(
          "Depth correction has not been validated for binary mode. ",
          "The default depth_quantile was chosen using count based entropy."
        ),
        call. = FALSE
      )
    }
    futile.logger::flog.info(
      "Depth correction target: %d counts; %d of %d cells will be corrected.",
      depth_target, n_cells_corrected, length(depth_totals)
    )
    depth_corrected_mat <- .entropy_downsample_matrix(mat, depth_target)
  }

  if (!is.null(features)) {
    mat <- mat[features, , drop = FALSE]
    if (depth_correct) {
      depth_corrected_mat <- depth_corrected_mat[features, , drop = FALSE]
    }
  }

  if (mode == "binary") {
    futile.logger::flog.info("Binarizing matrix (counts > 0) for binary entropy.")
    mat_bin <- Matrix::Matrix(mat > 0, sparse = TRUE)
  }

  ## load the estimator table
  if (mode == "counts") {
    futile.logger::flog.info("Reading the precomputed table (counts mode)")
    .entropy_cache$precomputed_bub <- readRDS(system.file("extdata", "precomputed_bub_5000.rds", package = "Entropy"))
  } else {
    futile.logger::flog.info("Reading the binary entropy lookup (binary mode)")
    .entropy_cache$binary_entropy_lookup <- readRDS(system.file("extdata", "binary_entropy_lookup_100.rds", package = "Entropy"))
  }

  #### chunks (k-means on the embedding, no Leiden)
  if (chunk_method == "auto") {
    futile.logger::flog.info("Preparing chunks by k-means on the embedding.")
    chunks <- prepare_chunks_emb(emb, chunk_cells = chunk_cells, cell_names = cell_names)
  }

  ## One pass of the neighbourhood entropy loop. Called once normally, and a
  ## second time on the downsampled counts when depth_correct = TRUE.
  compute_entropy_matrix <- function(mat_use, mat_bin_use) {
    pb <- progress::progress_bar$new(
      format = "  Processing [:bar] :percent | :current/:total cells",
      total = n_cells, clear = FALSE, width = 60
    )


    entropy_matrix <- matrix(
      0, nrow = nrow(mat_use), ncol = n_cells,
      dimnames = list(rownames(mat_use), cell_names)
    )

    if (chunk_method == "none") {
      for (i in seq_len(n_cells)) {
        cell_name <- cell_names[i]
        if (i %% 100 == 0) pb$tick(100)
        nn_idx <- nn_list[[cell_name]]
        if (mode == "binary") {
          entropy_matrix[, cell_name] <- calculate_entropy_for_cell_binary_fast(
            cell_name = cell_name, mat_bin = mat_bin_use, neighbors_idx = nn_idx)
        } else {
          entropy_matrix[, cell_name] <- calculate_entropy_for_cell(
            cell_name = cell_name, mat = mat_use, neighbors_idx = nn_idx)
        }
      }
      if (n_cells %% 100 != 0) pb$tick(n_cells %% 100)
    } else {
      run_chunk <- function(chunk_cell_names) {
        all_nn_idx <- unique(unlist(nn_list[chunk_cell_names]))
        if (mode == "binary") {
          mat_chunk <- mat_bin_use[, all_nn_idx, drop = FALSE]
        } else {
          mat_chunk <- mat_use[, all_nn_idx, drop = FALSE]
        }
        idx_remap <- stats::setNames(seq_along(all_nn_idx), all_nn_idx)
        chunk_matrix <- matrix(
          0, nrow = nrow(mat_use), ncol = length(chunk_cell_names),
          dimnames = list(rownames(mat_use), chunk_cell_names)
        )
        for (cell_name in chunk_cell_names) {
          nn_idx_full  <- nn_list[[cell_name]]
          nn_idx_chunk <- as.integer(idx_remap[as.character(nn_idx_full)])
          if (mode == "binary") {
            chunk_matrix[, cell_name] <- calculate_entropy_for_cell_binary_fast(
              cell_name = cell_name, mat_bin = mat_chunk, neighbors_idx = nn_idx_chunk)
          } else {
            chunk_matrix[, cell_name] <- calculate_entropy_for_cell(
              cell_name = cell_name, mat = mat_chunk, neighbors_idx = nn_idx_chunk)
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
        cn <- colnames(chunk_results[[i]])
        entropy_matrix[, cn] <- chunk_results[[i]]
      }
    }

    entropy_matrix
  }

  entropy_matrix <- compute_entropy_matrix(
    mat, if (mode == "binary") mat_bin else NULL
  )

  means <- colMeans(entropy_matrix)

  if (!is.null(output_path) && !depth_correct) {
    saveRDS(entropy_matrix, output_path)
    futile.logger::flog.info("Full entropy matrix saved to: %s", output_path)
  }

  out <- data.frame(Entropy_mean = means)
  rownames(out) <- cell_names

  if (depth_correct) {
    # The uncorrected dense matrix is not needed after its means are taken.
    # Release it before allocating the corrected one.
    rm(entropy_matrix)
    gc()

    depth_corrected_mat_bin <- NULL
    if (mode == "binary") {
      depth_corrected_mat_bin <- Matrix::Matrix(
        depth_corrected_mat > 0, sparse = TRUE
      )
    }
    futile.logger::flog.info("Computing depth-corrected entropy.")
    depth_corrected_entropy_matrix <- compute_entropy_matrix(
      depth_corrected_mat, depth_corrected_mat_bin
    )

    if (!is.null(output_path)) {
      saveRDS(depth_corrected_entropy_matrix, output_path)
      futile.logger::flog.info(
        "Full DEPTH CORRECTED entropy matrix saved to: %s", output_path
      )
    }

    out$Entropy_mean_depth_corrected <- colMeans(depth_corrected_entropy_matrix)
    out$depth_corrected <- unname(depth_corrected_cells[cell_names])
  }

  if (return_inputs) {
    futile.logger::flog.info(
      "Returning per-cell mean entropy, embedding, and obs for %d cells.", n_cells
    )
    return(list(
      entropy   = stats::setNames(means, cell_names),
      entropy_depth_corrected = if (depth_correct) {
        stats::setNames(out$Entropy_mean_depth_corrected, cell_names)
      } else {
        NULL
      },
      embedding = emb_full,
      obs       = inputs$obs,
      results   = out
    ))
  }

  futile.logger::flog.info("Returning per-cell mean entropy for %d cells.", n_cells)
  out
}


#' Internal: read the pieces we need out of an .h5ad file
#'
#' Returns a genes x cells dgCMatrix, a cells x dims embedding (or NULL), an
#' optional group vector, and the cell/feature names. Reading is isolated here so
#' the AnnData/reader API lives in exactly one place.
#'
#' AnnData stores cells x genes; transposing its matrix gives the genes x cells
#' CSC representation expected by the entropy engine.
#'
#' @keywords internal
read_h5ad_inputs <- function(path, layer = NULL, pca_key = "X_pca",
                             restrict_by = NULL, need_embedding = TRUE,
                             read_obs = FALSE) {

  check_h5ad_dependencies()

  adata <- anndataR::read_h5ad(path)

  cell_names    <- adata$obs_names
  feature_names <- adata$var_names

  # counts: layer if named, else X. AnnData is cells x genes -> transpose.
  x <- if (is.null(layer)) adata$X else adata$layers[[layer]]
  if (is.null(x)) {
    if (is.null(layer)) {
      stop("X was not found in the file.")
    }
    stop(sprintf("Layer '%s' not found in the file.", layer))
  }
  mat <- Matrix::t(methods::as(x, "CsparseMatrix"))  # genes x cells, dgCMatrix
  rownames(mat) <- feature_names
  colnames(mat) <- cell_names

  # embedding
  emb <- NULL
  if (need_embedding && !is.null(adata$obsm) && pca_key %in% names(adata$obsm)) {
    emb <- as.matrix(adata$obsm[[pca_key]])
    rownames(emb) <- cell_names
  }

  # restrict_by group vector
  groups <- NULL
  if (!is.null(restrict_by)) {
    if (!restrict_by %in% colnames(adata$obs)) {
      stop(sprintf(
        "restrict_by column '%s' not found in obs. Available: %s",
        restrict_by, paste(colnames(adata$obs), collapse = ", ")
      ))
    }
    groups <- as.character(adata$obs[[restrict_by]])
    names(groups) <- cell_names
  }

  obs <- NULL
  if (read_obs) {
    obs <- as.data.frame(adata$obs)
    rownames(obs) <- cell_names
  }

  list(mat = mat, emb = emb, groups = groups, obs = obs,
       cell_names = cell_names, feature_names = feature_names)
}


#' Internal: check optional dependencies for direct h5ad support
#'
#' @noRd
check_h5ad_dependencies <- function() {
  if (getRversion() < "4.5.0") {
    stop(
      "run_entropy_h5ad() requires R >= 4.5.0 because it uses anndataR. ",
      "The other Entropy functions remain available on older R versions.",
      call. = FALSE
    )
  }

  if (!requireNamespace("anndataR", quietly = TRUE)) {
    stop(
      "run_entropy_h5ad() requires the optional package 'anndataR' >= 1.0.0. ",
      "Install it with BiocManager::install(\"anndataR\").",
      call. = FALSE
    )
  }

  if (utils::packageVersion("anndataR") < "1.0.0") {
    stop(
      sprintf(
        paste0(
          "run_entropy_h5ad() requires anndataR >= 1.0.0, but version %s is installed. ",
          "Update it with BiocManager::install(\"anndataR\")."
        ),
        as.character(utils::packageVersion("anndataR"))
      ),
      call. = FALSE
    )
  }

  if (!requireNamespace("rhdf5", quietly = TRUE)) {
    stop(
      "run_entropy_h5ad() requires the optional Bioconductor package 'rhdf5'. ",
      "Install it with BiocManager::install(\"rhdf5\").",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


#' Internal: validate the h5ad expression matrix for the requested mode
#'
#' @keywords internal
validate_h5ad_matrix <- function(mat, mode, layer = NULL) {
  values <- if ("x" %in% methods::slotNames(mat)) mat@x else rep.int(1, length(mat@i))
  source_name <- if (is.null(layer)) "X" else sprintf("layer '%s'", layer)

  if (any(!is.finite(values))) {
    stop(sprintf("%s contains NA, NaN, or infinite values.", source_name))
  }
  if (any(values < 0)) {
    stop(sprintf(
      "%s contains negative values, but entropy requires a non-negative expression matrix.",
      source_name
    ))
  }

  if (mode == "counts" &&
      any(abs(values - round(values)) > sqrt(.Machine$double.eps))) {
    stop(sprintf(
      paste0(
        "%s contains non-integer values. Counts mode requires raw integer ",
        "counts; select the AnnData counts layer (for example, layer = \"counts\") ",
        "or use mode = \"binary\"."
      ),
      source_name
    ))
  }

  invisible(TRUE)
}


#' Internal: split cells into compact chunks by k-means on the embedding
#'
#' Seurat-free replacement for prepare_chunks(). No Leiden, no igraph. k-means on
#' the PCA embedding gives spatially compact clusters, which is the property that
#' makes sparse column subsetting cheap.
#'
#' @keywords internal
prepare_chunks_emb <- function(emb, chunk_cells, cell_names) {
  n <- nrow(emb)
  n_chunks <- max(1L, ceiling(n / chunk_cells))
  if (n_chunks == 1L) {
    return(list(cell_names))
  }
  km <- stats::kmeans(emb, centers = n_chunks, iter.max = 50, nstart = 3)
  chunks <- split(cell_names, km$cluster)
  futile.logger::flog.info(
    "Prepared %d chunks with median size %d cells.",
    length(chunks), as.integer(stats::median(lengths(chunks)))
  )
  unname(chunks)
}
