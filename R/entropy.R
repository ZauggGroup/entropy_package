#' Estimate gene expression entropy on a Seurat object
#'
#' This function computes gene-level entropy across neighboring cells for a specified assay
#' in a Seurat object. It optionally saves the entropy matrix to file and/or adds it as a
#' new assay whose name is the assay name followed by `"_entropy"`.
#'
#' @param seu A Seurat object.
#' @param assay Character. Name of the assay used for entropy calculation.
#' @param neighbors Optional neighbor definition. If `neighbors_method = "custom"`,
#'   this should be a named list where `neighbors[[cell_name]]` returns a vector of
#'   neighbor cell names or indices. If `neighbors_method = "seurat"`, this should be
#'   a Seurat neighbor object. If `NULL` and `neighbors_method = "auto"`, neighbors
#'   are computed internally.
#' @param neighbors_method Character. How to obtain neighbors. One of `"auto"`,
#'   `"custom"`, or `"seurat"`. `"auto"` computes neighbors inside the function.
#'   `"custom"` uses the user-provided neighbor list in `neighbors`. `"seurat"`
#'   uses the user-provided Seurat neighbor object in `neighbors`.
#' @param neighbor_assay Character. Assay to use for automatic neighbor computation
#'   when `neighbors_method = "auto"`. Defaults to `assay`.
#' @param features Optional character vector of feature names (e.g. genes) to
#'   restrict the entropy computation to. If `NULL` (default), entropy is
#'   computed for all features in `assay`. Subsetting applies only to `assay`
#'   (the scored features), not to `neighbor_assay` (which defines neighborhoods).
#'   In `"counts"` mode this also reduces runtime proportionally.
#' @param restrict_by Optional character string naming a column in
#'   `seu@meta.data` (e.g. `"celltype"`). If provided, entropy is computed
#'   within groups rather than across them. For `neighbors_method = "auto"`,
#'   neighbors are recomputed within each group so every cell gets exactly `k`
#'   neighbors from its own group. Groups with fewer than `k + 1` cells are too
#'   small to provide `k` neighbors; their cells are treated as orphans and
#'   return entropy 0, with a warning. For `"custom"` or `"seurat"`, the provided
#'   neighbor graph is pruned to within-group neighbors instead. Cells with `NA`
#'   in the column, or with no within-group neighbors, return entropy 0. If
#'   `NULL` (default), neighbors are used as-is. Works in both normal and chunked
#'   mode.
#' @param output_path Optional character string. If provided, saves the full entropy
#'   matrix to this file using `saveRDS()`.
#' @param add_assay Logical. If `TRUE` (default), adds the entropy matrix to the
#'   Seurat object as a new assay whose name ends in `"_entropy"`.
#' @param meta_name Optional character string. If provided, the per-cell mean
#'   entropy is added to the Seurat object's metadata under this column name.
#'   Works independently of `add_assay`. If either `add_assay = TRUE` or
#'   `meta_name` is set, the function returns the Seurat object; if both are off,
#'   it returns a data frame of mean entropy per cell.
#' @param k Integer. Number of nearest neighbors to use when computing neighbors
#'   internally. Must be at least 11, because the counts entropy estimator is
#'   undefined below 11 neighbors. Default is 15.
#' @param reduction Character. Name of an existing dimensional reduction in
#'   `seu` to use for neighbor and chunk computation (e.g. `"lsi"` for ATAC).
#'   If `NULL` (default), a PCA named `"pca_<neighbor_assay>"` is computed on
#'   `neighbor_assay` and used.
#' @param dims Integer vector or single integer. Dimensions of `reduction` to
#'   use. A single value `n` is treated as `1:n` (so `dims = 20` means `1:20`).
#'   A vector is used as given, so `dims = 2:20` skips the first component
#'   (useful for LSI/ATAC, where component 1 mostly tracks sequencing depth).
#'   Default is 20.
#' @param nfeatures Integer. Number of variable features to select before PCA-based
#'   neighbor computation. Default is 0.
#' @param mode Character. Entropy mode. `"counts"` uses the count-based entropy estimator
#'   (BUB precomputed table). `"binary"` binarizes the assay counts (values > 0) and
#'   computes Shannon entropy for a Bernoulli distribution across neighbors.
#' @param chunk_method Character. Whether to use chunked computation.
#'   One of `"none"` or `"auto"`. `"none"` uses the standard approach (default).
#'   `"auto"` splits cells into biologically similar groups and processes
#'   each group independently, significantly reducing memory and runtime
#'   for large datasets. Recommended for datasets with more than 10,000 cells.
#' @param chunk_cells Integer. Maximum number of cells per chunk when
#'   `chunk_method = "auto"`. Defaults to 2000. Values between 1000 and 5000
#'   are recommended.
#' @param parallelize Logical. Whether to parallelize chunk processing across
#'   multiple cores. Only applies when \code{chunk_method = "auto"}.
#'   Defaults to \code{FALSE}. Not supported on Windows.
#' @param n_cores Integer. Number of cores to use when \code{parallelize = TRUE}.
#'   Defaults to \code{parallel::detectCores() - 1}.
#' @return If `add_assay = TRUE`, returns the modified Seurat object with a new assay
#'   whose name ends in `"_entropy"`. If `add_assay = FALSE`, returns a data frame of mean
#'   entropy values per cell.
#' @export
#'
#' @examples
#' \dontrun{
#' data("pbmc_small", package = "SeuratObject")
#' run_entropy(pbmc_small, assay = "RNA", k = 15, mode = "counts")
#' }
run_entropy <- function(
    seu,
    assay,
    neighbors = NULL,
    neighbors_method = c("auto", "custom", "seurat"),
    neighbor_assay = assay,
    features = NULL,
    restrict_by = NULL,
    output_path = NULL,
    add_assay = TRUE,
    meta_name = NULL,
    k = 15,
    reduction = NULL,
    dims = 20,
    nfeatures = 0,
    mode = c("counts", "binary"),
    chunk_method = c("none", "auto"),
    chunk_cells = 2000,
    parallelize = FALSE,
    n_cores = parallel::detectCores() - 1L
){


  neighbors_method <- match.arg(neighbors_method)
  mode <- match.arg(mode)
  chunk_method <- match.arg(chunk_method)

  if (!assay %in% Seurat::Assays(seu)) {
    stop(sprintf("Assay '%s' not found in Seurat object.", assay))
  }
  Seurat::DefaultAssay(seu) <- assay
  futile.logger::flog.info("Assay '%s' found and set as default assay.", assay)

  if (!neighbor_assay %in% Seurat::Assays(seu)) {
    stop(sprintf("neighbor_assay '%s' not found in Seurat object.", neighbor_assay))
  }

  futile.logger::flog.info("Starting entropy estimation for assay: %s", assay)
  futile.logger::flog.info(
    "Assay '%s' contains %d cells and %d features",
    assay, ncol(seu), nrow(seu)
  )
  ### cehck K

  if (k < 11) {
    stop(sprintf("k must be at least 11 (got %d). The counts entropy estimator is undefined below 11 neighbors.", k))
  }


  ### check  feature

  if (!is.null(features)) {
    if (!is.character(features)) {
      stop("`features` must be a character vector of feature names, not ",
           class(features)[1], ". Indices are not accepted.")
    }

    features <- unique(features)
    assay_features <- rownames(seu[[assay]])
    missing <- setdiff(features, assay_features)

    if (length(missing) == length(features)) {
      stop(sprintf(
        "None of the %d requested features are present in assay '%s'. Check your feature names.",
        length(features), assay
      ))
    }

    if (length(missing) > 0) {
      stop(sprintf(
        "%d of %d requested features are not in assay '%s': %s%s",
        length(missing), length(features), assay,
        paste(utils::head(missing, 10), collapse = ", "),
        if (length(missing) > 10) ", ..." else ""
      ))
    }

    futile.logger::flog.info(
      "Restricting entropy computation to %d of %d features in assay '%s'.",
      length(features), length(assay_features), assay
    )
  }


  ### check parallelization

  if (parallelize && chunk_method == "none") {
    futile.logger::flog.warn(
      "parallelize = TRUE has no effect when chunk_method = 'none'."
    )
  }

  if (parallelize && .Platform$OS.type == "windows") {
    futile.logger::flog.warn(
      "parallelize = TRUE is not supported on Windows. Falling back to sequential."
    )
    parallelize <- FALSE
  }

  if (parallelize && n_cores > parallel::detectCores() - 1L) {
    futile.logger::flog.warn(
      "n_cores = %d exceeds available cores. Setting to %d.",
      n_cores, parallel::detectCores() - 1L
    )
    n_cores <- parallel::detectCores() - 1L
  }

  if (parallelize && n_cores < 1L) {
    stop("n_cores must be at least 1.")
  }


  ### chunk check
  n_cells <- ncol(seu)

  if (chunk_method == "none" && n_cells > 50000) {
    futile.logger::flog.warn(
      "Dataset has %d cells but chunk_method = 'none'. This may be very slow. Consider setting chunk_method = 'auto'.",
      n_cells
    )
  }

  if (chunk_method == "auto") {

    if (!is.numeric(chunk_cells) || length(chunk_cells) != 1 || chunk_cells != as.integer(chunk_cells)) {
      stop("chunk_cells must be a single integer.")
    }

    if (chunk_cells < 500) {
      stop(
        sprintf("chunk_cells = %d is too low. Minimum recommended value is 500. Very small chunks increase overhead without benefit.", chunk_cells)
      )
    }

    if (chunk_cells > 10000) {
      futile.logger::flog.warn(
        "chunk_cells = %d is very large. Chunks may not be small enough to provide a speedup. Values between 1000 and 5000 are recommended.",
        chunk_cells
      )
    }

    if (n_cells < chunk_cells * 2) {
      futile.logger::flog.warn(
        "Dataset has %d cells which is too small to benefit from chunking with chunk_cells = %d. Falling back to chunk_method = 'none'.",
        n_cells, chunk_cells
      )
      chunk_method <- "none"
    }

  }

  #### restrict_by check
  if (!is.null(restrict_by)) {
    if (!is.character(restrict_by) || length(restrict_by) != 1) {
      stop("`restrict_by` must be a single character string naming a metadata column.")
    }

    if (!restrict_by %in% colnames(seu@meta.data)) {
      stop(sprintf(
        "restrict_by column '%s' not found in seu@meta.data. Available columns: %s",
        restrict_by,
        paste(colnames(seu@meta.data), collapse = ", ")
      ))
    }

    n_na <- sum(is.na(seu@meta.data[[restrict_by]]))
    if (n_na > 0) {
      futile.logger::flog.warn(
        "Column '%s' has %d cells with NA. These cells belong to no group and will have no within-group neighbors (entropy 0).",
        restrict_by, n_na
      )
    }

    futile.logger::flog.info(
      "Neighbors will be restricted within groups defined by metadata column '%s'.",
      restrict_by
    )
  }




  #### Chunk check finished

  need_reduction <- (neighbors_method == "auto") || (chunk_method == "auto")

  if (need_reduction) {

    if (is.null(reduction)) {
      reduction <- paste0("pca_", neighbor_assay)

      if (!reduction %in% names(seu@reductions)) {

        # cross-version check: does scaled data exist for this assay? (v4 + v5)
        sd <- tryCatch(
          SeuratObject::GetAssayData(seu, assay = neighbor_assay, layer = "scale.data"),
          error = function(e) NULL
        )
        has_scale <- !is.null(sd) && nrow(sd) > 0

        if (!has_scale) {
          stop(sprintf(
            paste0(
              "Assay '%s' has no scaled data, which is required to compute PCA. ",
              "Please normalize and scale the assay first (e.g. NormalizeData() then ScaleData(), ",
              "or the SCT/ATAC equivalent for your data), or pass a precomputed reduction via `reduction=`."
            ),
            neighbor_assay
          ))
        }

        futile.logger::flog.info("Reduction '%s' not found. Running PCA on assay '%s'.", reduction, neighbor_assay)
        old_default <- Seurat::DefaultAssay(seu)
        Seurat::DefaultAssay(seu) <- neighbor_assay
        if (nfeatures > 0) {
          futile.logger::flog.info("Selecting top %s variable features.", nfeatures)
          seu <- Seurat::FindVariableFeatures(seu, nfeatures = nfeatures)
        }
        seu <- Seurat::RunPCA(
          seu,
          reduction.name = reduction,
          reduction.key  = paste0(gsub("[^a-zA-Z0-9]", "", reduction), "_"),
          verbose        = FALSE
        )
        Seurat::DefaultAssay(seu) <- old_default
      } else {
        futile.logger::flog.info("Reusing existing reduction '%s'.", reduction)
      }

    } else {
      if (!reduction %in% names(seu@reductions)) {
        stop(sprintf("reduction '%s' not found. Available: %s", reduction, paste(names(seu@reductions), collapse = ", ")))
      }
      if (nfeatures > 0) {
        futile.logger::flog.warn("nfeatures = %d ignored because reduction '%s' was supplied.", nfeatures, reduction)
      }
      futile.logger::flog.info("Using user-provided reduction '%s'.", reduction)
    }

    # single n -> 1:n, vector -> as given
    if (length(dims) == 1) {
      if (is.na(dims) || dims < 1) stop("dims must be a positive integer.")
      dims <- seq_len(dims)
    }
    dims <- as.integer(dims)
    if (any(is.na(dims)) || any(dims < 1)) stop("dims must be positive integers.")

    n_avail <- ncol(Seurat::Embeddings(seu[[reduction]]))
    if (max(dims) > n_avail) {
      futile.logger::flog.warn("Requested dims up to %d but '%s' has %d. Capping.", max(dims), reduction, n_avail)
      dims <- dims[dims <= n_avail]
    }
    if (length(dims) < 1) stop(sprintf("No valid dims after capping against '%s'.", reduction))
  }

  ########### check finishd
  ## NN creation

  restrict_method <- neighbors_method


  if (neighbors_method == "auto") {

    futile.logger::flog.info(
      "No neighbors provided. Computing neighbors using assay '%s'.",
      neighbor_assay
    )

    if (is.null(restrict_by)) {
      # ---- unrestricted: global kNN via FindNeighbors (as before) ----
      seu_neighbor <- seu
      Seurat::DefaultAssay(seu_neighbor) <- neighbor_assay
      neighbor_name <- paste0("entropy_", neighbor_assay, "_nn")

      seu_neighbor <- suppressMessages(Seurat::FindNeighbors(
        seu_neighbor,
        reduction = reduction,
        dims = dims,
        return.neighbor = TRUE,
        k.param = k,
        graph.name = neighbor_name,
        verbose = FALSE
      ))

      if (!neighbor_name %in% names(seu_neighbor@neighbors)) {
        stop("Computed Neighbor object '", neighbor_name,
             "' was not found in the Seurat object.")
      }

      neighbors <- seu_neighbor@neighbors[[neighbor_name]]
      neighbors_method <- "seurat"

    } else {
      # restricted: recompute kNN WITHIN each group, build a named list
      futile.logger::flog.info(
        "Recomputing neighbors within each group of '%s'.", restrict_by
      )

      grp <- as.character(seu@meta.data[[restrict_by]])
      names(grp) <- colnames(seu)
      emb <- Seurat::Embeddings(seu[[reduction]])[, dims, drop = FALSE]

      if (!all(colnames(seu) %in% rownames(emb))) {
        stop(sprintf("Reduction '%s' does not contain all cells; cannot restrict by group.", reduction))
      }


      neighbors <- vector("list", ncol(seu))
      names(neighbors) <- colnames(seu)

      groups <- split(colnames(seu), grp)  # NA cells are dropped by split()

      for (g in names(groups)) {
        g_cells <- groups[[g]]
        g_n     <- length(g_cells)

        if (g_n < k + 1L) {
          futile.logger::flog.warn(
            "Group '%s' has %d cells (< k+1 = %d); too small for %d neighbors. These cells get entropy 0.",
            g, g_n, k + 1L, k
          )
          for (cell in g_cells) neighbors[[cell]] <- integer(0)
          next
        }

        nn_g <- RANN::nn2(emb[g_cells, , drop = FALSE], k = k + 1L)$nn.idx
        # column 1 is self; map local row/col indices back to global cell indices
        global_idx <- match(g_cells, colnames(seu))
        for (r in seq_along(g_cells)) {
          nbr_row    <- nn_g[r, ]
          local_nbrs <- nbr_row[nbr_row != r]
          neighbors[[g_cells[r]]] <- global_idx[local_nbrs]
        }
      }

      # NA cells never entered any group -> ensure they exist as empty
      na_cells <- colnames(seu)[is.na(grp)]
      for (cell in na_cells) neighbors[[cell]] <- integer(0)

      neighbors_method <- "custom"

      avg_nn <- mean(lengths(neighbors))
      n_zero <- sum(lengths(neighbors) == 0)
      futile.logger::flog.info(
        "Within-group neighbors: mean %.1f per cell; %d cells with 0 neighbors.",
        avg_nn, n_zero
      )
    }
  } else if (neighbors_method == "custom") {
    if (is.null(neighbors)) {
      stop("`neighbors` must be provided when neighbors_method = 'custom'.")
    }
    if (!is.list(neighbors) || is.null(names(neighbors))) {
      stop(
        "`neighbors` must be a named list where neighbors[[cell_name]] ",
        "returns a vector of neighbor cell names or indices."
      )
    }




    cell_names <- colnames(seu)

    missing_cells <- setdiff(cell_names, names(neighbors))

    if (length(missing_cells) > 0) {
      stop(
        "The custom neighbors list is missing entries for ",
        length(missing_cells), " cells.\n",
        "Each cell must appear in names(neighbors), even if it has no neighbors.\n",
        "Example structure:\n",
        "neighbors[['cell1']] -> c('cell2','cell3',...)"
      )
    }
    futile.logger::flog.info("Using user-provided custom neighbors.")

  } else if (neighbors_method == "seurat") {
    if (is.null(neighbors)) {
      stop("`neighbors` must be provided when neighbors_method = 'seurat'.")
    }

    if (!inherits(neighbors, "Neighbor")) {
      stop(
        "For neighbors_method = 'seurat', `neighbors` must be a Seurat Neighbor object ",
        "(e.g. the result of Seurat::FindNeighbors(..., return.neighbor = TRUE))."
      )
    }

    if (!all(colnames(seu) %in% neighbors@cell.names)) {
      stop("Cell names in Seurat object do not match the provided Seurat neighbor object.")
    }

    futile.logger::flog.info("Using user-provided Seurat neighbor object.")
  }


  #### Filter neighbors within groups (custom / seurat only; auto handled itself)
  if (!is.null(restrict_by) && restrict_method %in% c("custom", "seurat")) {
    futile.logger::flog.info(
      "Filtering provided neighbors within groups of '%s'.", restrict_by
    )

    grp <- as.character(seu@meta.data[[restrict_by]])
    names(grp) <- colnames(seu)
    all_cells <- colnames(seu)

    filtered <- vector("list", length(all_cells))
    names(filtered) <- all_cells

    n_before <- numeric(length(all_cells))
    n_after  <- numeric(length(all_cells))

    for (i in seq_along(all_cells)) {
      cell   <- all_cells[i]
      nn_idx <- get_neighbors(cell, seu, neighbors, restrict_method)
      n_before[i] <- length(nn_idx)

      cell_grp <- grp[cell]
      if (is.na(cell_grp) || length(nn_idx) == 0) {
        filtered[[cell]] <- integer(0)
      } else {
        filtered[[cell]] <- nn_idx[grp[nn_idx] == cell_grp & !is.na(grp[nn_idx])]
      }
      n_after[i] <- length(filtered[[cell]])
    }

    neighbors <- filtered
    neighbors_method <- "custom"

    futile.logger::flog.info(
      "Within-group filtering: mean neighbors %.1f -> %.1f; %d cells now with 0 neighbors.",
      mean(n_before), mean(n_after), sum(n_after == 0)
    )
  }


  #### Orphan cells with < 11 neighbors (custom / seurat only; counts estimator floor)
  if (restrict_method %in% c("custom", "seurat")) {

    all_cells <- colnames(seu)
    pruned    <- vector("list", length(all_cells))
    names(pruned) <- all_cells

    n_orphaned <- 0L
    for (i in seq_along(all_cells)) {
      cell   <- all_cells[i]
      nn_idx <- get_neighbors(cell, seu, neighbors, neighbors_method)
      if (length(nn_idx) < 11L) {
        pruned[[cell]] <- integer(0)
        if (length(nn_idx) > 0L) n_orphaned <- n_orphaned + 1L
      } else {
        pruned[[cell]] <- nn_idx
      }
    }

    neighbors <- pruned
    neighbors_method <- "custom"

    if (n_orphaned > 0L) {
      futile.logger::flog.warn(
        "%d of %d cells have fewer than 11 neighbors and will return entropy 0 (counts estimator requires >= 11).",
        n_orphaned, length(all_cells)
      )
    }
  }

  ###################### NN finish


  mat <- Seurat::GetAssayData(seu, assay = assay, layer = "counts")

  if (!is.null(features)) {
    mat <- mat[features, , drop = FALSE]
  }

  if (mode == "binary") {
    futile.logger::flog.info("Binarizing matrix (counts > 0) for binary entropy.")
    mat_bin <- Matrix::Matrix(mat > 0, sparse = TRUE)
  }

  # Get all cell names
  cell_names <- colnames(seu)
  n_cells <- length(cell_names)

  # Initialize the progress bar
  pb <- progress::progress_bar$new(
    format = "  Processing [:bar] :percent | :current/:total cells",
    total = n_cells,
    clear = FALSE,
    width = 60
  )

  #.data_loaded <<- FALSE

  ## prepare matrix
  if (mode == "counts") {
    futile.logger::flog.info("Reading the precomputed table (counts mode)")
    .entropy_cache$precomputed_bub <- readRDS(system.file(
      "extdata", "precomputed_bub_5000.rds",
      package = "Entropy"
    ))
  } else if (mode == "binary") {
    futile.logger::flog.info("Reading the binary entropy lookup (binary mode)")
    .entropy_cache$binary_entropy_lookup <- readRDS(system.file(
      "extdata", "binary_entropy_lookup_100.rds",
      package = "Entropy"
    ))
  } else {
    stop("Unknown mode: ", mode, ". Expected 'counts' or 'binary'.")
  }

  #### Prepare chunks
  if (chunk_method == "auto") {
    futile.logger::flog.info("Preparing chunks for chunked computation.")
    chunks <- prepare_chunks(seu, chunk_cells = chunk_cells, reduction = reduction, dims = dims, neighbor_assay = neighbor_assay)
  }

  # Apply with progress
  entropy_matrix <- matrix(
    0,
    nrow     = nrow(mat),
    ncol     = n_cells,
    dimnames = list(rownames(mat), cell_names)
  )

  if (chunk_method == "none") {
    for (i in seq_along(cell_names)) {
      cell_name <- cell_names[i]
      if (i %% 100 == 0) pb$tick(100)
      nn_idx <- get_neighbors(cell_name, seu, neighbors, neighbors_method)
      if (mode == "binary") {
        entropy_matrix[, cell_name] <- calculate_entropy_for_cell_binary_fast(
          cell_name = cell_name, mat_bin = mat_bin, neighbors_idx = nn_idx)
      } else {
        entropy_matrix[, cell_name] <- calculate_entropy_for_cell(
          cell_name = cell_name, mat = mat, neighbors_idx = nn_idx)
      }
    }
    if (n_cells %% 100 != 0) pb$tick(n_cells %% 100)
  } else {
    run_chunk <- function(chunk_cell_names) {
      all_nn_idx <- unique(unlist(lapply(chunk_cell_names, function(cell) {
        get_neighbors(cell, seu, neighbors, neighbors_method)
      })))
      if (mode == "binary") {
        mat_chunk <- mat_bin[, all_nn_idx, drop = FALSE]
      } else {
        mat_chunk <- mat[, all_nn_idx, drop = FALSE]
      }
      idx_remap <- stats::setNames(seq_along(all_nn_idx), all_nn_idx)
      chunk_matrix <- matrix(
        0,
        nrow     = nrow(mat),
        ncol     = length(chunk_cell_names),
        dimnames = list(rownames(mat), chunk_cell_names)
      )
      for (cell_name in chunk_cell_names) {
        nn_idx_full  <- get_neighbors(cell_name, seu, neighbors, neighbors_method)
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
      chunk_cell_names <- colnames(chunk_results[[i]])
      entropy_matrix[, chunk_cell_names] <- chunk_results[[i]]
    }
  }



  if (!is.null(output_path)) {
    saveRDS(entropy_matrix, output_path)
    futile.logger::flog.info("Entropy matrix saved to: %s", output_path)
  }


  means <- colMeans(entropy_matrix)

  if (add_assay) {
    entropy_assay <- if (inherits(seu[[assay]], "Assay5")) {
      SeuratObject::CreateAssay5Object(counts = entropy_matrix)
    } else {
      SeuratObject::CreateAssayObject(counts = entropy_matrix)
    }
    assay_name <- paste0(assay, "_entropy")
    seu[[assay_name]] <- entropy_assay
    futile.logger::flog.info("Entropy matrix added to Seurat object as '%s'.", assay_name)
  }

  if (!is.null(meta_name) || add_assay) {
    col <- if (!is.null(meta_name)) meta_name else paste0("entropy_mean_", assay)

    if (col %in% Seurat::Assays(seu)) {
      new_col <- paste0(col, "_meta")
      futile.logger::flog.warn(
        "meta_name '%s' collides with an assay name in the object. Using '%s' instead.",
        col, new_col
      )
      col <- new_col
    }

    seu[[col]] <- means
    futile.logger::flog.info("Mean entropy added to metadata column '%s'.", col)
  }

  # return rule: object if we modified it (assay or meta), else the data frame
  if (add_assay || !is.null(meta_name)) {
    return(seu)
  } else {
    return(as.data.frame(means))
  }

}







