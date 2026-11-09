#' Estimate gene expression entropy on a Seurat object
#'
#' This function computes gene-level entropy across neighboring cells for a specified assay
#' in a Seurat object. It optionally saves the entropy matrix to file and/or adds it as a
#' new assay whose name is the assay name followed by `"_entropy"`. When depth
#' correction is requested, neighbours still come from the original data. Only
#' the raw counts supplied to the entropy estimator are thinned.
#'
#' A cell is never its own neighbour. `k` always means `k` other cells,
#' whichever way the graph is supplied.
#'
#' When the graph is computed here it is cached at
#' `seu@misc$entropy_neighbors` and reused on the next call whose cells, `k`,
#' `reduction`, `dims`, `neighbor_assay` and `restrict_by` all match. The cache
#' survives only if the Seurat object is returned, so it is not kept when both
#' `add_assay` and `meta_name` are off and `depth_correct` is `FALSE`. A
#' `neighbors` argument supplied by the caller always wins over the cache.
#'
#' The dense gene-by-cell output is estimated before scoring. Above 2 GB the
#' function warns and continues normally unless `block_dir` is set. With a
#' block directory, genes are scored in restartable, ordered blocks and only
#' per-cell sums remain in memory; use [entropy_read_genes()] to retrieve rows.
#'
#' @param seu A Seurat object.
#' @param assay Character. Name of the assay used for entropy calculation.
#' @param neighbors Optional neighbor definition. Either a named list where
#'   `neighbors[[cell_name]]` gives that cell's neighbours as cell names or
#'   positions, or a Seurat `Neighbor` object, or `NULL` (default) to compute
#'   the graph internally. A supplied list must not contain the cell itself;
#'   it is rejected with a message if it does. A `Neighbor` object has its
#'   self entry removed automatically, so build it with `k.param = k + 1` if
#'   you want `k` neighbours left. Every cell must have the same number of
#'   distinct neighbours, empty lists aside. Numeric positions must be finite
#'   integers. [entropy_neighbors()] builds a valid list.
#' @param neighbors_method Character or `NULL` (default). Leave it unset and
#'   it is read from `neighbors`: `NULL` gives `"auto"`, a list gives
#'   `"custom"`, a `Neighbor` object gives `"seurat"`. Passing it explicitly
#'   still works, and errors if it contradicts what was supplied.
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
#'   return entropy 0, with a warning. For a supplied list or `Neighbor`
#'   object, the provided graph is pruned to within-group neighbors instead.
#'   Cells with `NA` in the column, or with no within-group neighbors, return
#'   entropy 0. If pruning leaves a non-empty entry with fewer than `k`
#'   neighbours, the call errors rather than comparing unequal neighbourhood
#'   sizes. For a global supplied graph, omit `neighbors` so the automatic path
#'   rebuilds `k` neighbours independently within each group. If `NULL`
#'   (default), neighbors are used as-is. Works in both normal and chunked mode.
#' @param output_path Optional character string. If provided, saves the full entropy
#'   matrix to this file using `saveRDS()`. It cannot be combined with
#'   `block_dir`, because doing so would reassemble the full matrix in memory.
#' @param block_dir Optional directory for restartable gene-blocked output. If
#'   set, gene-level entropy is written in ordered `block_*.rds` files plus a
#'   `manifest.rds`, while only per-cell sums are retained in memory. Existing
#'   valid blocks are reused; a manifest mismatch is an error. Blocking cannot
#'   be combined with `add_assay = TRUE`.
#' @param block_genes Optional positive integer giving genes per block. It is
#'   only valid when `block_dir` is set. If `NULL`, a size near 1 GB per block
#'   is chosen from the number of cells and entropy passes.
#' @param add_assay Logical. If `TRUE` (default), adds the entropy matrix to the
#'   Seurat object as a new assay whose name ends in `"_entropy"`. With
#'   `depth_correct = TRUE`, this assay contains only the depth-corrected entropy
#'   matrix.
#' @param meta_name Optional character string. If provided, the per-cell mean
#'   entropy is added to the Seurat object's metadata under this column name.
#'   Works independently of `add_assay`. If either `add_assay = TRUE` or
#'   `meta_name` is set, the function returns the Seurat object; if both are off,
#'   it returns a data frame of mean entropy per cell.
#' @param k Integer. Number of nearest neighbors, excluding the cell itself,
#'   used when computing neighbors internally. Must be at least 11, because the
#'   counts entropy estimator is undefined below 11 neighbors. Default is 15.
#'   With a supplied `neighbors` list or `Neighbor` object, `k` is checked
#'   rather than used: every cell must have exactly `k` non-self neighbours.
#'   Leave `k` unset and it is taken from the list instead.
#' @param force_neighbors Logical. If `TRUE`, ignore any cached neighbour graph
#'   on the object and recompute. Default `FALSE`.
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
#' @param depth_correct Logical. If `TRUE`, also compute entropy after thinning
#'   each cell's raw counts to a common library depth. The default is `FALSE`.
#' @param depth_quantile Numeric. Quantile of full-library cell totals used to
#'   choose the target depth when `target_depth` is `NULL`. Must be strictly
#'   between 0 and 1. The default is 0.05.
#' @param target_depth Optional positive integer target for count thinning. If
#'   supplied, it overrides `depth_quantile`.
#' @return If `add_assay = TRUE`, returns the modified Seurat object with a new
#'   assay whose name ends in `"_entropy"`. If `depth_correct = TRUE`, the object
#'   also contains uncorrected and corrected per-cell means, a logical correction
#'   indicator, and depth-correction provenance. If both `add_assay = FALSE` and
#'   `depth_correct = FALSE`, a data frame of mean entropy values is returned
#'   unless `meta_name` is supplied. In blocked mode the same return rules apply;
#'   `block_dir` is recorded in `seu@misc$entropy_provenance` when a Seurat object
#'   is returned, and in an `entropy_provenance` attribute on a returned data
#'   frame.
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
    neighbors_method = NULL,
    neighbor_assay = assay,
    features = NULL,
    restrict_by = NULL,
    output_path = NULL,
    block_dir = NULL,
    block_genes = NULL,
    add_assay = TRUE,
    meta_name = NULL,
    k = 15,
    reduction = NULL,
    dims = 20,
    nfeatures = 0,
    force_neighbors = FALSE,
    mode = c("counts", "binary"),
    chunk_method = c("none", "auto"),
    chunk_cells = 2000,
    parallelize = FALSE,
    n_cores = parallel::detectCores() - 1L,
    depth_correct = FALSE,
    depth_quantile = 0.05,
    target_depth = NULL
){


  k_given <- !missing(k)

  if (is.null(neighbors_method)) {
    neighbors_method <- if (is.null(neighbors)) {
      "auto"
    } else if (inherits(neighbors, "Neighbor")) {
      "seurat"
    } else {
      "custom"
    }
  } else {
    neighbors_method <- match.arg(neighbors_method,
                                  c("auto", "custom", "seurat"))
    if (neighbors_method == "auto" && !is.null(neighbors)) {
      stop(paste0(
        "`neighbors_method = \"auto\"` computes the graph internally, but ",
        "`neighbors` was also supplied. Drop one of the two."), call. = FALSE)
    }
    if (neighbors_method == "custom" && inherits(neighbors, "Neighbor")) {
      stop(paste0(
        "`neighbors` is a Seurat Neighbor object, not a list. Use ",
        "`neighbors_method = \"seurat\"`, or leave `neighbors_method` unset ",
        "and it is detected for you."), call. = FALSE)
    }
  }
  if (length(force_neighbors) != 1L || !is.logical(force_neighbors) ||
      is.na(force_neighbors)) {
    stop("`force_neighbors` must be TRUE or FALSE.", call. = FALSE)
  }
  mode <- match.arg(mode)
  chunk_method <- match.arg(chunk_method)

  if (!is.null(block_dir) &&
      (length(block_dir) != 1L || !is.character(block_dir) ||
       is.na(block_dir) || !nzchar(block_dir))) {
    stop("`block_dir` must be one non-empty directory path or `NULL`.",
         call. = FALSE)
  }
  if (!is.null(block_genes) &&
      (length(block_genes) != 1L || !is.numeric(block_genes) ||
       !is.finite(block_genes) || block_genes != floor(block_genes) ||
       block_genes < 1L)) {
    stop("`block_genes` must be one positive integer or `NULL`.",
         call. = FALSE)
  }
  if (is.null(block_dir) && !is.null(block_genes)) {
    stop("`block_genes` requires `block_dir`; set a directory for block output.",
         call. = FALSE)
  }
  if (!is.null(block_dir) && add_assay) {
    stop(paste0(
      "Gene blocking cannot be combined with `add_assay = TRUE`, because an ",
      "entropy assay requires the full matrix in memory. Set `add_assay = FALSE` ",
      "and read selected rows with `entropy_read_genes()`."
    ), call. = FALSE)
  }
  if (!is.null(block_dir) && !is.null(output_path)) {
    stop(paste0(
      "`output_path` cannot be combined with `block_dir`, because it would ",
      "reassemble the full entropy matrix. Use the block files and ",
      "`entropy_read_genes()` instead."
    ), call. = FALSE)
  }
  if (!is.null(block_genes)) block_genes <- as.integer(block_genes)

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
    stop("`target_depth` must be one positive integer or `NULL`.",
         call. = FALSE)
  }

  if (depth_correct && add_assay) {
    message(paste(
      "NOTE: depth_correct = TRUE and add_assay = TRUE.",
      "The stored entropy assay will contain the DEPTH CORRECTED matrix only.",
      "The uncorrected per-cell mean is still written to metadata.",
      sep = "\n"
    ))
  }

  if (!assay %in% Seurat::Assays(seu)) {
    stop(sprintf("Assay '%s' not found in Seurat object.", assay))
  }
  Seurat::DefaultAssay(seu) <- assay
  futile.logger::flog.info("Assay '%s' found and set as default assay.", assay)

  depth_target <- NULL
  n_cells_below_target <- 0L
  n_cells_corrected <- 0L
  depth_corrected_cells <- NULL
  full_depth_counts <- NULL
  if (depth_correct) {
    full_depth_counts <- suppressWarnings(tryCatch(
      SeuratObject::GetAssayData(seu, assay = assay, layer = "counts"),
      error = function(error) NULL
    ))
    if (is.null(full_depth_counts) || length(dim(full_depth_counts)) != 2L ||
        nrow(full_depth_counts) == 0L || ncol(full_depth_counts) == 0L) {
      stop(sprintf(
        paste0(
          "Depth correction requires raw counts in the counts layer of assay ",
          "'%s'."
        ),
        assay
      ), call. = FALSE)
    }
    count_values <- if (inherits(full_depth_counts, "sparseMatrix")) {
      full_depth_counts@x
    } else {
      as.vector(full_depth_counts)
    }
    integer_like <- all(
      is.finite(count_values) & count_values >= 0 &
        abs(count_values - round(count_values)) <= sqrt(.Machine$double.eps)
    )
    if (!integer_like) {
      stop(sprintf(
        paste0(
          "Depth correction requires raw counts in the counts layer of assay ",
          "'%s'; the stored values are not integer-like."
        ),
        assay
      ), call. = FALSE)
    }

    depth_totals <- Matrix::colSums(full_depth_counts)
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
    names(depth_corrected_cells) <- colnames(full_depth_counts)

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
    message(sprintf(
      "Depth correction target: %d counts; %d of %d cells will be corrected.",
      depth_target, n_cells_corrected, length(depth_totals)
    ))
  }

  if (!neighbor_assay %in% Seurat::Assays(seu)) {
    stop(sprintf("neighbor_assay '%s' not found in Seurat object.", neighbor_assay))
  }

  futile.logger::flog.info("Starting entropy estimation for assay: %s", assay)
  futile.logger::flog.info(
    "Assay '%s' contains %d cells and %d features",
    assay, ncol(seu), nrow(seu)
  )
  ### cehck K

  if (length(k) != 1L || !is.numeric(k) || !is.finite(k) ||
      k != floor(k)) {
    stop("k must be a single integer.", call. = FALSE)
  }
  k <- as.integer(k)
  if (k < 11L) {
    stop(sprintf("k must be at least 11 (got %d). The counts entropy estimator is undefined below 11 neighbors.", k))
  }


  ### check  feature

  assay_features <- rownames(seu[[assay]])
  if (!is.null(features)) {
    if (!is.character(features)) {
      stop("`features` must be a character vector of feature names, not ",
           class(features)[1], ". Indices are not accepted.")
    }

    features <- unique(features)
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
  scored_features <- if (is.null(features)) assay_features else features

  .entropy_warn_large_output(
    n_cells = ncol(seu), n_genes = length(scored_features),
    depth_correct = depth_correct, block_dir = block_dir
  )
  if (!is.null(block_dir)) {
    if (file.exists(block_dir) && !dir.exists(block_dir)) {
      stop(sprintf("`block_dir` exists but is not a directory: '%s'.", block_dir),
           call. = FALSE)
    }
    if (!dir.exists(block_dir) &&
        !dir.create(block_dir, recursive = TRUE, showWarnings = FALSE)) {
      stop(sprintf("Could not create `block_dir`: '%s'.", block_dir),
           call. = FALSE)
    }
    block_dir <- normalizePath(block_dir, winslash = "/", mustWork = TRUE)
    if (is.null(block_genes)) {
      passes <- if (depth_correct) 2 else 1
      block_genes <- max(1L, floor(
        .entropy_block_target_bytes / (ncol(seu) * 8 * passes)
      ))
      block_genes <- min(as.integer(block_genes), length(scored_features))
      futile.logger::flog.info(
        "Using %d genes per block (approximately 1 GB of entropy output).",
        block_genes
      )
    }
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


  # Cache lookup. Only consulted when the graph would be computed here; a
  # neighbour list handed in by the caller always wins.
  nn_key <- NULL
  if (neighbors_method == "auto") {
    nn_key <- list(
      cells = colnames(seu), k = k, reduction = reduction, dims = dims,
      neighbor_assay = neighbor_assay,
      restrict_by = if (is.null(restrict_by)) NA_character_ else restrict_by
    )
    cached <- seu@misc$entropy_neighbors
    if (!force_neighbors && !is.null(cached) &&
        identical(attr(cached, "signature"), nn_key)) {
      neighbors <- unclass(cached)
      attributes(neighbors) <- NULL
      names(neighbors) <- names(cached)
      neighbors_method <- "custom"
      futile.logger::flog.info(
        "Reusing cached neighbors from seu@misc$entropy_neighbors (built %s).",
        format(attr(cached, "built"))
      )
    } else if (!is.null(cached) && !force_neighbors) {
      futile.logger::flog.info(
        "Cached neighbors do not match this call (%s). Recomputing.",
        paste(.entropy_nn_diff(attr(cached, "signature"), nn_key),
              collapse = ", ")
      )
    }
  }

  if (neighbors_method == "auto") {

    futile.logger::flog.info(
      "No neighbors provided. Computing neighbors using assay '%s'.",
      neighbor_assay
    )

    if (is.null(restrict_by)) {
      # ---- unrestricted: global kNN on the embedding ----
      # k + 1 then drop self, so k means k OTHER cells. The restricted branch
      # below has always done this. The unrestricted branch used to include
      # the cell itself as one of the k observations.
      emb <- Seurat::Embeddings(seu[[reduction]])[, dims, drop = FALSE]
      if (nrow(emb) <= k) {
        stop(sprintf(paste0(
          "Cannot compute %d neighbors from only %d cells: global automatic ",
          "neighboring requires at least k + 1 = %d cells."),
          k, nrow(emb), k + 1L), call. = FALSE)
      }

      nn <- suppressMessages(Seurat::FindNeighbors(
        emb,
        k.param = k + 1L,
        return.neighbor = TRUE,
        verbose = FALSE
      ))

      neighbors <- .entropy_nn_trim(.entropy_nn_drop_self(nn@nn.idx), k)
      names(neighbors) <- colnames(seu)
      neighbors_method <- "custom"

      futile.logger::flog.info(
        "Computed %d non-self neighbors per cell on reduction '%s' (dims %d:%d).",
        k, reduction, min(dims), max(dims)
      )

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
        # Remove self by identity (ties mean it need not be column 1), then map
        # local row/column indices back to global cell indices.
        global_idx <- match(g_cells, colnames(seu))
        for (r in seq_along(g_cells)) {
          nbr_row    <- nn_g[r, ]
          local_nbrs <- utils::head(nbr_row[nbr_row != r], k)
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

    # Validate the graph before it becomes estimator input or enters the cache.
    .entropy_nn_check(neighbors, colnames(seu), k, k_given = TRUE)

    nn_store <- neighbors
    attributes(nn_store) <- NULL
    names(nn_store) <- colnames(seu)
    attr(nn_store, "signature") <- nn_key
    attr(nn_store, "source") <- "computed"
    attr(nn_store, "built") <- Sys.time()
    class(nn_store) <- "entropy_neighbors"
    seu@misc$entropy_neighbors <- nn_store

  } else if (restrict_method != "auto") {
    # Custom list or Seurat Neighbor object, normalised to one list form.
    # Self is never a neighbour: a Neighbor object has its self entry
    # dropped, a user list that contains self is rejected.
    if (is.null(neighbors)) {
      stop(sprintf("`neighbors` must be provided when neighbors_method = '%s'.",
                   neighbors_method), call. = FALSE)
    }
    if (neighbors_method == "seurat" && !inherits(neighbors, "Neighbor")) {
      stop(paste0(
        "For neighbors_method = 'seurat', `neighbors` must be a Seurat Neighbor ",
        "object (e.g. Seurat::FindNeighbors(..., return.neighbor = TRUE))."),
        call. = FALSE)
    }

    neighbors <- .entropy_nn_normalise(neighbors, colnames(seu), k)
    k_resolved <- .entropy_nn_check(neighbors, colnames(seu), k, k_given)
    if (!k_given) {
      futile.logger::flog.info(
        "k not supplied; set to %d from the neighbour list.", k_resolved)
      k <- k_resolved
      if (k < 11L) {
        stop(sprintf(paste0(
          "The supplied neighbour list gives k = %d. The counts entropy ",
          "estimator is undefined below 11 neighbors."), k), call. = FALSE)
      }
    }
    neighbors_method <- "custom"

    futile.logger::flog.info(
      "Using user-provided neighbors (%s): %d non-self per cell, %d cells with 0.",
      restrict_method, k, sum(lengths(neighbors) == 0L)
    )
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
      nn_idx <- get_neighbors(cell, seu, neighbors)
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

    # Filtering a caller-supplied graph must not silently create different
    # neighbourhood sizes. Empty entries remain valid for cells with no
    # within-group neighbours, but every non-empty entry must still contain k.
    .entropy_nn_check_restricted(neighbors, all_cells, k)

    futile.logger::flog.info(
      "Within-group filtering: mean neighbors %.1f -> %.1f; %d cells now with 0 neighbors.",
      mean(n_before), mean(n_after), sum(n_after == 0)
    )
  }

  ###################### NN finish


  if (depth_correct) {
    depth_corrected_full <- .entropy_downsample_matrix(
      full_depth_counts, depth_target
    )
  }


  full_mat <- Seurat::GetAssayData(seu, assay = assay, layer = "counts")
  cell_names <- colnames(seu)
  n_cells <- length(cell_names)

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
  } else {
    chunks <- NULL
  }

  if (is.null(block_dir)) {
    mat <- full_mat[scored_features, , drop = FALSE]
    mat_bin <- if (mode == "binary") {
      futile.logger::flog.info("Binarizing matrix (counts > 0) for binary entropy.")
      Matrix::Matrix(mat > 0, sparse = TRUE)
    } else {
      NULL
    }
    entropy_matrix <- .entropy_compute_matrix(
      mat, mat_bin, mode, cell_names, neighbors, chunk_method, chunks,
      parallelize, n_cores
    )
    means <- colMeans(entropy_matrix)
    if (!is.null(output_path)) {
      saveRDS(entropy_matrix, output_path)
      futile.logger::flog.info("Entropy matrix saved to: %s", output_path)
    }

    if (depth_correct) {
      rm(entropy_matrix)
      gc()
      depth_corrected_mat <- depth_corrected_full[
        scored_features, , drop = FALSE
      ]
      depth_corrected_mat_bin <- if (mode == "binary") {
        Matrix::Matrix(depth_corrected_mat > 0, sparse = TRUE)
      } else {
        NULL
      }
      depth_corrected_entropy_matrix <- .entropy_compute_matrix(
        depth_corrected_mat, depth_corrected_mat_bin, mode, cell_names,
        neighbors, chunk_method, chunks, parallelize, n_cores
      )
      depth_corrected_means <- colMeans(depth_corrected_entropy_matrix)
    }
  } else {
    prepared <- .entropy_prepare_manifest(
      block_dir = block_dir, cells = cell_names, genes = scored_features,
      block_genes = block_genes, k = k, assay = assay, mode = mode,
      reduction = reduction, dims = dims, restrict_by = restrict_by,
      depth_correct = depth_correct, target_depth = depth_target
    )
    entropy_sums <- stats::setNames(numeric(n_cells), cell_names)
    if (depth_correct) {
      depth_corrected_sums <- stats::setNames(numeric(n_cells), cell_names)
    }

    for (block_index in seq_along(prepared$indices)) {
      genes_i <- scored_features[prepared$indices[[block_index]]]
      block_path <- file.path(
        block_dir, prepared$manifest$blocks$file[[block_index]]
      )
      if (file.exists(block_path)) {
        block <- tryCatch(
          readRDS(block_path),
          error = function(error) stop(sprintf(
            "Cannot read entropy block file '%s': %s",
            block_path, conditionMessage(error)
          ), call. = FALSE)
        )
        .entropy_validate_block(
          block, prepared$manifest, block_index, block_path
        )
        futile.logger::flog.info(
          "Reusing gene block %d/%d (%d genes).",
          block_index, length(prepared$indices), length(genes_i)
        )
      } else {
        futile.logger::flog.info(
          "Computing gene block %d/%d (%d genes).",
          block_index, length(prepared$indices), length(genes_i)
        )
        mat_i <- full_mat[genes_i, , drop = FALSE]
        mat_bin_i <- if (mode == "binary") {
          Matrix::Matrix(mat_i > 0, sparse = TRUE)
        } else {
          NULL
        }
        entropy_i <- .entropy_compute_matrix(
          mat_i, mat_bin_i, mode, cell_names, neighbors, chunk_method, chunks,
          parallelize, n_cores
        )
        corrected_i <- NULL
        if (depth_correct) {
          corrected_mat_i <- depth_corrected_full[
            genes_i, , drop = FALSE
          ]
          corrected_bin_i <- if (mode == "binary") {
            Matrix::Matrix(corrected_mat_i > 0, sparse = TRUE)
          } else {
            NULL
          }
          corrected_i <- .entropy_compute_matrix(
            corrected_mat_i, corrected_bin_i, mode, cell_names, neighbors,
            chunk_method, chunks, parallelize, n_cores
          )
        }
        block <- list(
          block = as.integer(block_index), genes = genes_i,
          cells = cell_names, entropy = entropy_i,
          entropy_depth_corrected = corrected_i
        )
        .entropy_save_rds_atomic(block, block_path)
        rm(mat_i, mat_bin_i, entropy_i, corrected_i)
        if (exists("corrected_mat_i", inherits = FALSE)) {
          rm(corrected_mat_i, corrected_bin_i)
        }
      }
      entropy_sums <- entropy_sums + colSums(block$entropy)
      if (depth_correct) {
        depth_corrected_sums <- depth_corrected_sums +
          colSums(block$entropy_depth_corrected)
      }
      rm(block)
      gc(verbose = FALSE)
    }
    means <- entropy_sums / length(scored_features)
    if (depth_correct) {
      depth_corrected_means <- depth_corrected_sums / length(scored_features)
    }
  }

  if (add_assay) {
    assay_entropy_matrix <- if (depth_correct) {
      depth_corrected_entropy_matrix
    } else {
      entropy_matrix
    }
    entropy_assay <- if (inherits(seu[[assay]], "Assay5")) {
      SeuratObject::CreateAssay5Object(counts = assay_entropy_matrix)
    } else {
      SeuratObject::CreateAssayObject(counts = assay_entropy_matrix)
    }
    assay_name <- paste0(assay, "_entropy")
    seu[[assay_name]] <- entropy_assay
    futile.logger::flog.info("Entropy matrix added to Seurat object as '%s'.", assay_name)
  }

  if (!is.null(meta_name) || add_assay || depth_correct) {
    requested_col <- if (!is.null(meta_name)) {
      meta_name
    } else {
      paste0("entropy_mean_", assay)
    }
    col <- requested_col

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

  if (depth_correct) {
    corrected_col <- paste0(col, "_depth_corrected")
    if (corrected_col %in% Seurat::Assays(seu)) {
      new_corrected_col <- paste0(corrected_col, "_meta")
      futile.logger::flog.warn(
        "Corrected metadata name '%s' collides with an assay name in the object. Using '%s' instead.",
        corrected_col, new_corrected_col
      )
      corrected_col <- new_corrected_col
    }
    seu[[corrected_col]] <- depth_corrected_means
    seu[[paste0("depth_corrected_", assay)]] <-
      depth_corrected_cells[colnames(seu)]
    futile.logger::flog.info(
      "Depth-corrected mean entropy added to metadata column '%s'.",
      corrected_col
    )
  }

  if (add_assay || !is.null(meta_name) || depth_correct ||
      !is.null(block_dir)) {
    entropy_provenance <- list(
      k = k,
      assay = assay,
      mode = mode,
      reduction = reduction,
      dims = dims,
      restrict_by = restrict_by,
      depth_correct = depth_correct,
      depth_quantile = depth_quantile,
      target_depth = depth_target,
      n_cells_below_target = n_cells_below_target,
      n_cells_corrected = n_cells_corrected
    )
    if (!is.null(block_dir)) {
      entropy_provenance$block_dir <- block_dir
    }
    seu@misc$entropy_provenance <- entropy_provenance
  }

  # return rule: object if we modified it (assay or meta), else the data frame
  if (add_assay || !is.null(meta_name) || depth_correct) {
    return(seu)
  } else {
    result <- as.data.frame(means)
    if (!is.null(block_dir)) {
      attr(result, "entropy_provenance") <- entropy_provenance
    }
    return(result)
  }

}


