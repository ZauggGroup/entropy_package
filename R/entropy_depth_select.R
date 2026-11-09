#' Select a depth-correction quantile for entropy
#'
#' `entropy_depth_select()` runs a restartable sweep over count-depth
#' quantiles and compares two absolute variance components: variation aligned
#' with the original detected-gene count, and all remaining variation. The
#' detected-gene count, neighbour graph, and gene panel are fixed before the
#' sweep. The uncorrected score is included as `q = 1`, so correction is
#' recommended only when a tested quantile has the largest weighted objective.
#'
#' The retention weight `w` controls how much retained variation is valued
#' relative to removed depth-aligned variation. The default, `w = 1.5`, was
#' calibrated against three datasets: two bone marrow datasets where correction
#' was accepted and one DS Atlas dataset where it was refused. It is a
#' defensible default, not a derived constant.
#'
#' This function currently requires a Seurat object. AnnData and h5ad paths are
#' not supported by the depth-selection sweep.
#'
#' In the covariate tables and report, `thinnable = TRUE` is reserved for the
#' package-built `nFeature_<assay>` and `nCount_<assay>` covariates computed
#' from the original counts. User-supplied covariates are shown under
#' "Thinning does not act on this", even when their names or interpretation
#' are depth-related.
#'
#' @param object A Seurat object containing raw counts and an existing
#'   dimensional reduction.
#' @param group One character string naming a non-missing grouping column in
#'   `object` metadata.
#' @param q_grid Numeric depth quantiles to test. Values must be greater than
#'   zero and no greater than one, with at most three decimal places. `q = 1`
#'   is added automatically. The default tests 0.95 down to 0.05.
#' @param w Positive numeric retention weight. Larger values prefer milder or
#'   no correction. The default is 1.5.
#' @param covariates Optional additional covariates. Supply metadata column
#'   names, a numeric data frame, or a named list with one finite numeric value
#'   per cell. Original detected genes, original count depth, and mean neighbour
#'   distance are always included.
#' @param outdir Output directory for checkpoints, the saved selection object,
#'   three CSV tables, and the PDF report.
#' @param per_cluster Logical. Add the per-group before-and-after heatmap page
#'   when correction is recommended. The default is `TRUE`.
#' @param seed Integer seed used reproducibly for every thinning run. A separate
#'   deterministic seed derived from it selects the fixed gene panel.
#' @param ... Named arguments passed to [run_entropy()]. `assay` is required.
#'   Common arguments are `k`, `reduction`, `dims`, `mode`, `chunk_method`,
#'   `chunk_cells`, `parallelize`, and `n_cores`. A supplied `features` vector is
#'   used as the fixed panel; otherwise up to 5000 genes are sampled once.
#'
#' @return An object of class `entropy_depth_selection`. The expensive result
#'   is also saved as `depth_selection.rds` in `outdir`. Use
#'   [entropy_depth_report()] to redraw its report without rerunning entropy.
#' @export
entropy_depth_select <- function(
    object,
    group,
    q_grid = NULL,
    w = 1.5,
    covariates = NULL,
    outdir,
    per_cluster = TRUE,
    seed = 1,
    ...) {

  if (!inherits(object, "Seurat")) {
    stop("`object` must be a Seurat object. h5ad input is not supported here.",
         call. = FALSE)
  }
  if (length(group) != 1L || !is.character(group) || is.na(group) ||
      !nzchar(group) || !group %in% colnames(object[[]])) {
    stop("`group` must name one metadata column in `object`.", call. = FALSE)
  }
  group_values <- object[[]][[group]]
  if (anyNA(group_values) || length(unique(group_values)) < 2L) {
    stop("The grouping column must be non-missing and contain at least two groups.",
         call. = FALSE)
  }
  group_values <- as.character(group_values)
  names(group_values) <- colnames(object)

  q_grid <- .entropy_depth_q_grid(q_grid)
  w <- .entropy_depth_number(w, "w", lower = 0, lower_open = TRUE)
  seed <- .entropy_depth_integer(seed, "seed")
  .entropy_depth_logical(per_cluster, "per_cluster")
  outdir <- .entropy_depth_outdir(outdir)

  dots <- list(...)
  if (length(dots) &&
      (is.null(names(dots)) || any(!nzchar(names(dots))) ||
       anyDuplicated(names(dots)))) {
    stop("Every argument in `...` must have a unique name.", call. = FALSE)
  }
  if (!"assay" %in% names(dots) || length(dots$assay) != 1L ||
      !is.character(dots$assay) || is.na(dots$assay) || !nzchar(dots$assay)) {
    stop("`...` must include one non-empty `assay` value.", call. = FALSE)
  }
  forbidden <- intersect(
    names(dots),
    c("seu", "neighbors", "neighbors_method", "depth_correct",
      "depth_quantile", "target_depth", "add_assay", "meta_name",
      "output_path", "block_dir", "block_genes", "restrict_by",
      "force_neighbors")
  )
  if (length(forbidden)) {
    stop(sprintf(
      paste0(
        "The depth sweep controls `%s` so every quantile is comparable. ",
        "Remove this argument from `...`."
      ),
      forbidden[[1L]]
    ), call. = FALSE)
  }

  assay <- dots$assay
  if (!assay %in% Seurat::Assays(object)) {
    stop(sprintf("Assay '%s' was not found in `object`.", assay),
         call. = FALSE)
  }
  k <- if (is.null(dots$k)) 15L else dots$k
  k <- .entropy_depth_integer(k, "k", lower = 11L)
  mode <- if (is.null(dots$mode)) "counts" else match.arg(
    dots$mode, c("counts", "binary")
  )
  neighbor_assay <- if (is.null(dots$neighbor_assay)) assay else
    dots$neighbor_assay
  if (length(neighbor_assay) != 1L || !is.character(neighbor_assay) ||
      is.na(neighbor_assay) || !neighbor_assay %in% Seurat::Assays(object)) {
    stop("`neighbor_assay` must name one assay in `object`.", call. = FALSE)
  }
  reduction <- if (is.null(dots$reduction)) {
    paste0("pca_", neighbor_assay)
  } else {
    dots$reduction
  }
  if (length(reduction) != 1L || !is.character(reduction) ||
      is.na(reduction) || !reduction %in% names(object@reductions)) {
    stop(sprintf(
      paste0(
        "Reduction '%s' is not available. Build it before depth selection, ",
        "then pass its name with `reduction`."
      ),
      paste(reduction, collapse = ", ")
    ), call. = FALSE)
  }
  dims <- if (is.null(dots$dims)) 20L else dots$dims
  dims <- .entropy_nn_dims(dims, object, reduction)

  counts <- suppressWarnings(tryCatch(
    SeuratObject::GetAssayData(object, assay = assay, layer = "counts"),
    error = function(error) NULL
  ))
  if (is.null(counts) || length(dim(counts)) != 2L ||
      nrow(counts) < 1L || ncol(counts) < 2L) {
    stop(sprintf("Assay '%s' must contain a non-empty counts layer.", assay),
         call. = FALSE)
  }
  count_values <- if (inherits(counts, "sparseMatrix")) counts@x else
    as.vector(counts)
  integer_like <- all(
    is.finite(count_values) & count_values >= 0 &
      abs(count_values - round(count_values)) <= sqrt(.Machine$double.eps)
  )
  if (!integer_like) {
    stop(sprintf(
      "Assay '%s' contains non-integer counts. Supply the raw counts assay.",
      assay
    ), call. = FALSE)
  }

  cell_names <- colnames(object)
  assay_features <- rownames(counts)
  totals <- as.numeric(Matrix::colSums(counts))
  names(totals) <- cell_names
  # This is computed once from original counts because changing x during the
  # sweep changes the confound definition and invalidates the frontier.
  detected_genes <- as.numeric(Matrix::colSums(counts > 0))
  names(detected_genes) <- cell_names

  supplied_features <- dots$features
  panel_seed <- as.integer((as.double(seed) + 104729) %% 2147483647)
  if (panel_seed < 1L) panel_seed <- 1L
  if (is.null(supplied_features)) {
    panel_size <- min(5000L, length(assay_features))
    gene_panel <- .entropy_depth_with_seed(
      panel_seed, sort(sample(assay_features, panel_size))
    )
  } else {
    if (!is.character(supplied_features) || !length(supplied_features) ||
        anyNA(supplied_features) || anyDuplicated(supplied_features)) {
      stop("`features` in `...` must be a non-empty unique character vector.",
           call. = FALSE)
    }
    missing_features <- setdiff(supplied_features, assay_features)
    if (length(missing_features)) {
      stop(sprintf(
        "Feature '%s' is not present in assay '%s'. Fix the fixed gene panel.",
        missing_features[[1L]], assay
      ), call. = FALSE)
    }
    gene_panel <- supplied_features
  }

  fixed_neighbors <- entropy_neighbors(
    object, k = k, reduction = reduction, dims = dims,
    neighbor_assay = neighbor_assay
  )
  embedding <- Seurat::Embeddings(object[[reduction]])[, dims, drop = FALSE]
  embedding <- embedding[cell_names, , drop = FALSE]
  mean_knn_dist <- .entropy_depth_mean_distance(embedding, fixed_neighbors)

  covariate_data <- .entropy_depth_covariates(
    object = object, covariates = covariates, detected_genes = detected_genes,
    totals = totals, mean_knn_dist = mean_knn_dist, assay = assay
  )
  neighbor_signature_id <- .entropy_depth_signature_id(
    lapply(fixed_neighbors, as.integer)
  )
  analysis_signature_id <- .entropy_depth_signature_id(list(
    group = group_values,
    covariates = covariate_data
  ))

  run_dots <- dots
  run_dots$assay <- NULL
  run_dots$k <- NULL
  run_dots$mode <- NULL
  run_dots$neighbor_assay <- NULL
  run_dots$reduction <- NULL
  run_dots$dims <- NULL
  run_dots$features <- NULL
  count_signature <- .entropy_depth_count_signature(counts)
  signature <- list(
    cells = cell_names,
    genes = gene_panel,
    count_signature = count_signature,
    neighbor_signature_id = neighbor_signature_id,
    analysis_signature_id = analysis_signature_id,
    group = group,
    covariates = names(covariate_data),
    k = k,
    reduction = reduction,
    dims = dims,
    assay = assay,
    neighbor_assay = neighbor_assay,
    mode = mode,
    seed = seed,
    panel_seed = panel_seed,
    q_grid = q_grid,
    run_arguments = run_dots,
    package_version = as.character(utils::packageVersion("Entropy"))
  )
  signature_id <- .entropy_depth_signature_id(signature)
  .entropy_depth_prepare_manifest(outdir, signature)

  runs <- vector("list", length(q_grid))
  for (index in seq_along(q_grid)) {
    q <- q_grid[[index]]
    cache_path <- file.path(
      outdir, sprintf("score_q%04d.rds", as.integer(round(1000 * q)))
    )
    if (file.exists(cache_path)) {
      cached <- tryCatch(
        readRDS(cache_path),
        error = function(error) stop(sprintf(
          "Cannot read cached score '%s': %s",
          cache_path, conditionMessage(error)
        ), call. = FALSE)
      )
      .entropy_depth_validate_cache(
        cached, q, signature_id, cell_names, cache_path
      )
      futile.logger::flog.info("Reusing depth score for q = %.3f.", q)
      runs[[index]] <- cached
      next
    }

    futile.logger::flog.info(
      "Computing depth score %d/%d at q = %.3f.",
      index, length(q_grid), q
    )
    common <- c(list(
      seu = object,
      assay = assay,
      neighbors = fixed_neighbors,
      neighbor_assay = neighbor_assay,
      features = gene_panel,
      add_assay = FALSE,
      k = k,
      reduction = reduction,
      dims = dims,
      mode = mode
    ), run_dots)

    if (identical(q, 1)) {
      answer <- .entropy_depth_with_seed(
        seed,
        do.call(run_entropy, c(common, list(depth_correct = FALSE)))
      )
      if (!is.data.frame(answer) || nrow(answer) != length(cell_names)) {
        stop("The uncorrected entropy run returned an unexpected result.",
             call. = FALSE)
      }
      score <- as.numeric(answer[[1L]])
      target_depth <- max(totals)
      n_cells_corrected <- 0L
    } else {
      metadata_name <- "entropy_depth_selection"
      answer <- .entropy_depth_with_seed(
        seed,
        do.call(run_entropy, c(common, list(
          depth_correct = TRUE,
          depth_quantile = q,
          meta_name = metadata_name
        )))
      )
      corrected_name <- paste0(metadata_name, "_depth_corrected")
      if (!inherits(answer, "Seurat") ||
          !corrected_name %in% colnames(answer[[]])) {
        stop("The depth-corrected entropy run returned an unexpected result.",
             call. = FALSE)
      }
      score <- as.numeric(answer[[]][[corrected_name]])
      target_depth <- answer@misc$entropy_provenance$target_depth
      n_cells_corrected <- answer@misc$entropy_provenance$n_cells_corrected
    }
    names(score) <- cell_names
    expected_target <- if (identical(q, 1)) {
      max(totals)
    } else {
      floor(stats::quantile(totals, q, names = FALSE, type = 7))
    }
    if (length(score) != length(cell_names) || any(!is.finite(score)) ||
        !identical(as.numeric(target_depth), as.numeric(expected_target)) ||
        !identical(as.integer(n_cells_corrected),
                   as.integer(sum(totals > expected_target)))) {
      stop(sprintf(
        "The q = %.3f result does not match its cells or target depth. Use a new `outdir`.",
        q
      ), call. = FALSE)
    }
    cached <- list(
      q = q,
      target_depth = as.integer(target_depth),
      n_cells_corrected = as.integer(n_cells_corrected),
      score = score,
      signature_id = signature_id
    )
    .entropy_depth_save(cached, cache_path)
    runs[[index]] <- cached
    rm(answer, cached)
    gc(verbose = FALSE)
  }

  frontier <- .entropy_depth_frontier(runs, detected_genes, w)
  choice <- .entropy_depth_choose(frontier, w)
  selected_run <- runs[[match(choice$q, q_grid)]]
  baseline_run <- runs[[match(1, q_grid)]]

  baseline_screen <- .entropy_depth_screen(
    baseline_run$score, covariate_data, group_values,
    state = "uncorrected", q = 1
  )
  selected_screen <- if (identical(choice$q, 1)) {
    .entropy_depth_screen(
      baseline_run$score, covariate_data, group_values,
      state = "selected", q = 1
    )
  } else {
    .entropy_depth_screen(
      selected_run$score, covariate_data, group_values,
      state = "selected", q = choice$q
    )
  }

  decision_row <- frontier[frontier$q == choice$q, , drop = FALSE][1L, ]
  decision <- list(
    correct = !identical(choice$q, 1),
    q = choice$q,
    target_depth = as.integer(decision_row$target_depth),
    n_cells_thinned = as.integer(decision_row$n_cells_corrected),
    n_cells = length(cell_names),
    w = w,
    reason = choice$reason
  )
  cell_data <- data.frame(
    cell = cell_names,
    group = unname(group_values[cell_names]),
    detected_genes = unname(detected_genes[cell_names]),
    entropy_uncorrected = unname(baseline_run$score[cell_names]),
    entropy_selected = unname(selected_run$score[cell_names]),
    covariate_data,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  result <- list(
    decision = decision,
    frontier = frontier,
    covariates = baseline_screen$summary[, c(
      "covariate", "thinnable", "r_squared", "spearman_rho",
      "within_median", "joint_r_squared"
    ), drop = FALSE],
    within = baseline_screen$within[, c(
      "group", "covariate", "n_cells", "within_slope",
      "within_slope_std", "within_r_squared"
    ), drop = FALSE],
    comparison = rbind(baseline_screen$summary, selected_screen$summary),
    provenance = list(
      cells = cell_names,
      k = k,
      reduction = reduction,
      dims = dims,
      assay = assay,
      mode = mode,
      seed = seed,
      q_grid = q_grid,
      package_version = as.character(utils::packageVersion("Entropy")),
      date = as.character(Sys.Date()),
      group = group,
      gene_panel = gene_panel,
      panel_seed = panel_seed,
      neighbor_signature = attr(fixed_neighbors, "signature"),
      signature_id = signature_id,
      call = match.call()
    ),
    within_comparison = rbind(
      baseline_screen$within, selected_screen$within
    ),
    cell_data = cell_data
  )
  class(result) <- "entropy_depth_selection"

  .entropy_depth_save(result, file.path(outdir, "depth_selection.rds"))
  entropy_depth_report(result, outdir = outdir, per_cluster = per_cluster)
  result
}


#' Re-render a depth-selection report
#'
#' Writes the fixed-page PDF and the three CSV tables from an existing
#' `entropy_depth_selection` object. No entropy score is recomputed.
#'
#' @param x An `entropy_depth_selection` object returned by
#'   [entropy_depth_select()].
#' @param outdir Output directory. It is created recursively when needed.
#' @param per_cluster Logical. Add the per-group heatmap when correction was
#'   selected. The default is `TRUE`.
#'
#' @return A named character vector of normalized output paths, invisibly.
#' @export
entropy_depth_report <- function(x, outdir, per_cluster = TRUE) {
  if (!inherits(x, "entropy_depth_selection")) {
    stop("`x` must be an entropy_depth_selection object.", call. = FALSE)
  }
  .entropy_depth_logical(per_cluster, "per_cluster")
  outdir <- .entropy_depth_outdir(outdir)

  paths <- c(
    pdf = file.path(outdir, "depth_selection_report.pdf"),
    frontier = file.path(outdir, "frontier.csv"),
    covariates = file.path(outdir, "covariates.csv"),
    within = file.path(outdir, "within_group.csv")
  )
  utils::write.csv(x$frontier, paths[["frontier"]], row.names = FALSE)
  utils::write.csv(x$covariates, paths[["covariates"]], row.names = FALSE)
  utils::write.csv(x$within, paths[["within"]], row.names = FALSE)

  pages <- .entropy_depth_report_pages(x, per_cluster)
  .entropy_depth_write_pdf(paths[["pdf"]], pages)
  paths <- vapply(paths, normalizePath, character(1), mustWork = TRUE)
  invisible(paths)
}


#' Print a depth-selection recommendation
#'
#' @param x An `entropy_depth_selection` object.
#' @param ... Additional arguments, currently ignored.
#'
#' @return `x`, invisibly.
#' @export
print.entropy_depth_selection <- function(x, ...) {
  if (!inherits(x, "entropy_depth_selection")) {
    stop("`x` must be an entropy_depth_selection object.", call. = FALSE)
  }
  cat(paste(.entropy_depth_verdict_lines(x), collapse = "\n"), "\n")
  invisible(x)
}


.entropy_depth_q_grid <- function(q_grid) {
  if (is.null(q_grid)) {
    q_grid <- c(0.95, 0.90, 0.80, 0.70, 0.60, 0.50, 0.40,
                0.30, 0.25, 0.20, 0.15, 0.10, 0.05)
  }
  if (!is.numeric(q_grid) || !length(q_grid) || any(!is.finite(q_grid)) ||
      any(q_grid <= 0 | q_grid > 1)) {
    stop("`q_grid` must contain finite numbers greater than 0 and no greater than 1.",
         call. = FALSE)
  }
  q_grid <- sort(unique(c(as.numeric(q_grid), 1)), decreasing = TRUE)
  labels <- round(1000 * q_grid)
  if (any(abs(labels / 1000 - q_grid) > sqrt(.Machine$double.eps)) ||
      anyDuplicated(labels)) {
    stop("`q_grid` values must be distinct at three decimal places.",
         call. = FALSE)
  }
  q_grid
}


.entropy_depth_number <- function(value, argument, lower = -Inf,
                                  upper = Inf, lower_open = FALSE) {
  if (length(value) != 1L || !is.numeric(value) || !is.finite(value)) {
    stop(sprintf("`%s` must be one finite number in the allowed range.", argument),
         call. = FALSE)
  }
  valid_lower <- if (lower_open) value > lower else value >= lower
  if (!valid_lower || value > upper) {
    stop(sprintf("`%s` must be one finite number in the allowed range.", argument),
         call. = FALSE)
  }
  as.numeric(value)
}


.entropy_depth_integer <- function(value, argument, lower = -.Machine$integer.max) {
  if (length(value) != 1L || !is.numeric(value) || !is.finite(value) ||
      value != floor(value) || value < lower || value > .Machine$integer.max) {
    stop(sprintf("`%s` must be one integer of at least %d.", argument, lower),
         call. = FALSE)
  }
  as.integer(value)
}


.entropy_depth_logical <- function(value, argument) {
  if (length(value) != 1L || !is.logical(value) || is.na(value)) {
    stop(sprintf("`%s` must be TRUE or FALSE.", argument), call. = FALSE)
  }
  invisible(value)
}


.entropy_depth_outdir <- function(outdir) {
  if (missing(outdir) || length(outdir) != 1L || !is.character(outdir) ||
      is.na(outdir) || !nzchar(outdir)) {
    stop("`outdir` must be one non-empty directory path.", call. = FALSE)
  }
  if (file.exists(outdir) && !dir.exists(outdir)) {
    stop(sprintf("`outdir` exists but is not a directory: '%s'.", outdir),
         call. = FALSE)
  }
  if (!dir.exists(outdir) &&
      !dir.create(outdir, recursive = TRUE, showWarnings = FALSE)) {
    stop(sprintf("Could not create `outdir`: '%s'.", outdir), call. = FALSE)
  }
  normalizePath(outdir, winslash = "/", mustWork = TRUE)
}


.entropy_depth_with_seed <- function(seed, value) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv,
                                inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  force(value)
}


.entropy_depth_save <- function(value, path) {
  temporary <- tempfile(paste0(basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(value, temporary)
  if (!file.rename(temporary, path)) {
    stop(sprintf("Could not finalize depth-selection file '%s'.", path),
         call. = FALSE)
  }
  invisible(path)
}


.entropy_depth_signature_id <- function(signature) {
  path <- tempfile("entropy-depth-signature-", fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(signature, path, version = 2)
  unname(as.character(tools::md5sum(path)))
}


.entropy_depth_count_signature <- function(counts) {
  row_totals <- as.numeric(Matrix::rowSums(counts))
  col_totals <- as.numeric(Matrix::colSums(counts))
  values <- if (inherits(counts, "sparseMatrix")) counts@x else as.vector(counts)
  c(
    rows = nrow(counts),
    columns = ncol(counts),
    nonzero = sum(values != 0),
    total = sum(values),
    sum_squares = sum(values^2),
    row_weighted = sum(row_totals * seq_along(row_totals)),
    column_weighted = sum(col_totals * seq_along(col_totals))
  )
}


.entropy_depth_prepare_manifest <- function(outdir, signature) {
  path <- file.path(outdir, "depth_selection_manifest.rds")
  if (file.exists(path)) {
    stored <- tryCatch(
      readRDS(path),
      error = function(error) stop(sprintf(
        "Cannot read depth-selection manifest '%s': %s",
        path, conditionMessage(error)
      ), call. = FALSE)
    )
    fields <- union(names(stored), names(signature))
    for (field in fields) {
      if (!field %in% names(stored) || !field %in% names(signature) ||
          !identical(stored[[field]], signature[[field]])) {
        stop(sprintf(
          paste0(
            "Depth-selection manifest mismatch for field '%s'. ",
            "Use a new `outdir` for the changed run."
          ),
          field
        ), call. = FALSE)
      }
    }
    return(invisible(path))
  }
  existing <- list.files(outdir, pattern = "^score_q[0-9]+[.]rds$")
  if (length(existing)) {
    stop(paste0(
      "Cached depth scores exist without their manifest. ",
      "Use a new `outdir`, or remove the orphaned scores after checking them."
    ), call. = FALSE)
  }
  .entropy_depth_save(signature, path)
}


.entropy_depth_validate_cache <- function(cached, q, signature_id, cells, path) {
  valid <- is.list(cached) && identical(cached$q, q) &&
    identical(cached$signature_id, signature_id) &&
    identical(names(cached$score), cells) &&
    length(cached$score) == length(cells) &&
    all(is.finite(cached$score)) &&
    length(cached$target_depth) == 1L && is.finite(cached$target_depth) &&
    length(cached$n_cells_corrected) == 1L &&
    is.finite(cached$n_cells_corrected)
  if (!valid) {
    stop(sprintf(
      "Cached depth score '%s' does not match this run. Use a new `outdir`.",
      path
    ), call. = FALSE)
  }
  invisible(TRUE)
}


.entropy_depth_mean_distance <- function(embedding, neighbors) {
  values <- numeric(nrow(embedding))
  for (index in seq_len(nrow(embedding))) {
    neighbor_index <- neighbors[[index]]
    differences <- sweep(
      embedding[neighbor_index, , drop = FALSE], 2L,
      embedding[index, ], FUN = "-"
    )
    values[[index]] <- mean(sqrt(rowSums(differences^2)))
  }
  names(values) <- rownames(embedding)
  values
}


.entropy_depth_covariates <- function(object, covariates, detected_genes,
                                      totals, mean_knn_dist, assay) {
  cell_names <- colnames(object)
  detected_name <- paste0("nFeature_", assay)
  count_name <- paste0("nCount_", assay)
  result <- data.frame(
    detected_genes, totals, mean_knn_dist,
    check.names = FALSE, stringsAsFactors = FALSE
  )
  names(result)[1:2] <- c(detected_name, count_name)
  rownames(result) <- cell_names

  metadata <- object[[]]
  if (is.null(covariates)) {
    candidates <- c("percent_mt", "percent.mt", "S.Score", "G2M.Score")
    candidates <- candidates[candidates %in% colnames(metadata)]
    if (length(candidates)) {
      extra <- metadata[, candidates, drop = FALSE]
    } else {
      extra <- NULL
    }
  } else if (is.character(covariates)) {
    missing <- setdiff(covariates, colnames(metadata))
    if (length(missing)) {
      stop(sprintf(
        "Covariate column '%s' is not in `object` metadata.", missing[[1L]]
      ), call. = FALSE)
    }
    extra <- metadata[, covariates, drop = FALSE]
  } else {
    if (is.list(covariates) && !is.data.frame(covariates)) {
      extra <- tryCatch(
        as.data.frame(covariates, check.names = FALSE),
        error = function(error) NULL
      )
    } else {
      extra <- covariates
    }
    if (!is.data.frame(extra) || nrow(extra) != length(cell_names)) {
      stop("`covariates` must have one row or value per Seurat cell.",
           call. = FALSE)
    }
    extra_rows <- rownames(extra)
    automatic_rows <- identical(
      extra_rows, as.character(seq_len(nrow(extra)))
    )
    if (automatic_rows) {
      rownames(extra) <- cell_names
    } else if (!is.null(extra_rows)) {
      if (!setequal(extra_rows, cell_names)) {
        stop("Covariate row names must match the Seurat cells.", call. = FALSE)
      }
      extra <- extra[cell_names, , drop = FALSE]
    }
  }

  if (!is.null(extra) && ncol(extra)) {
    if (is.null(names(extra)) || any(!nzchar(names(extra))) ||
        anyDuplicated(names(extra))) {
      stop("Additional covariates must have unique non-empty names.",
           call. = FALSE)
    }
    reserved <- intersect(names(extra), names(result))
    if (length(reserved)) {
      stop(sprintf(
        "Covariate '%s' is computed by the sweep and must not be supplied.",
        reserved[[1L]]
      ), call. = FALSE)
    }
    numeric_columns <- vapply(extra, is.numeric, logical(1L))
    finite_columns <- vapply(
      extra, function(value) all(is.finite(value)), logical(1L)
    )
    if (!all(numeric_columns & finite_columns)) {
      stop("Every additional covariate must be finite and numeric.",
           call. = FALSE)
    }
    result <- cbind(result, extra)
  }
  result
}


.entropy_depth_split <- function(score, detected_genes) {
  total <- stats::var(score)
  if (!is.finite(total) || total <= 0) {
    stop("Entropy has no finite cell-to-cell variance. Depth selection is not possible.",
         call. = FALSE)
  }
  fit <- stats::lm(score ~ detected_genes)
  r_squared <- unname(summary(fit)$r.squared)
  if (!is.finite(r_squared)) r_squared <- 0
  c(
    V_total = total,
    R2 = r_squared,
    V_confound = r_squared * total,
    V_leftover = (1 - r_squared) * total
  )
}


.entropy_depth_frontier <- function(runs, detected_genes, w) {
  rows <- lapply(runs, function(run) {
    split <- .entropy_depth_split(run$score, detected_genes)
    data.frame(
      q = run$q,
      target_depth = run$target_depth,
      n_cells_corrected = run$n_cells_corrected,
      V_total = unname(split[["V_total"]]),
      R2 = unname(split[["R2"]]),
      V_confound = unname(split[["V_confound"]]),
      V_leftover = unname(split[["V_leftover"]]),
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  result <- result[order(result$q, decreasing = TRUE), , drop = FALSE]
  rownames(result) <- NULL
  reference <- result[result$q == 1, , drop = FALSE][1L, ]
  if (!is.finite(reference$V_leftover) || reference$V_leftover <= 0) {
    stop(paste0(
      "The uncorrected score has no variation outside its linear depth fit. ",
      "The retention frontier cannot be estimated."
    ), call. = FALSE)
  }
  if (reference$V_confound <= sqrt(.Machine$double.eps) * reference$V_total) {
    result$confound_removed <- 0
  } else {
    result$confound_removed <- 1 -
      result$V_confound / reference$V_confound
  }
  result$retained <- result$V_leftover / reference$V_leftover
  result$objective <- result$confound_removed + w * result$retained
  gain <- c(NA_real_, diff(result$confound_removed))
  result$price_step <- c(NA_real_, -diff(result$retained)) / gain
  result$price_step[!is.finite(result$price_step) | is.na(gain) | gain <= 0] <-
    NA_real_
  result$price_cumulative <- (1 - result$retained) / result$confound_removed
  result$price_cumulative[
    !is.finite(result$price_cumulative) | result$confound_removed <= 0
  ] <- NA_real_
  baseline <- which(result$q == 1)[[1L]]
  result$confound_removed[[baseline]] <- 0
  result$retained[[baseline]] <- 1
  result$objective[[baseline]] <- w
  result$price_step[[baseline]] <- NA_real_
  result$price_cumulative[[baseline]] <- NA_real_
  result
}


.entropy_depth_choose <- function(frontier, w) {
  objective <- frontier$confound_removed + w * frontier$retained
  objective[frontier$q == 1] <- w
  best <- max(objective, na.rm = TRUE)
  tolerance <- sqrt(.Machine$double.eps) * max(1, abs(best))
  candidates <- which(abs(objective - best) <= tolerance)
  index <- candidates[[which.max(frontier$q[candidates])]]
  q <- frontier$q[[index]]
  reason <- if (identical(q, 1)) {
    sprintf(
      "No tested quantile beats leaving the data uncorrected at w = %s.",
      format(w, trim = TRUE)
    )
  } else {
    sprintf(
      "q = %.3f has the highest objective over the grid including no correction.",
      q
    )
  }
  list(q = q, index = index, reason = reason)
}


.entropy_depth_model_r_squared <- function(response, predictor) {
  if (stats::sd(response) == 0 || stats::sd(predictor) == 0) return(NA_real_)
  value <- suppressWarnings(unname(summary(stats::lm(response ~ predictor))$r.squared))
  if (length(value) == 1L && is.finite(value)) value else NA_real_
}


.entropy_depth_screen <- function(score, covariates, groups, state, q) {
  varying <- vapply(covariates, function(value) stats::sd(value) > 0, logical(1L))
  joint_r_squared <- if (any(varying) && stats::sd(score) > 0) {
    model_data <- data.frame(.score = score, covariates[, varying, drop = FALSE],
                             check.names = FALSE)
    value <- suppressWarnings(unname(summary(
      stats::lm(.score ~ ., data = model_data)
    )$r.squared))
    if (is.finite(value)) value else NA_real_
  } else {
    NA_real_
  }
  detected_pattern <- "^nFeature_"
  count_pattern <- "^nCount_"
  summary_rows <- vector("list", ncol(covariates))
  within_rows <- vector("list", ncol(covariates) * length(unique(groups)))
  within_index <- 0L

  for (index in seq_along(covariates)) {
    covariate_name <- names(covariates)[[index]]
    predictor <- covariates[[index]]
    r_squared <- .entropy_depth_model_r_squared(score, predictor)
    rho <- if (stats::sd(score) == 0 || stats::sd(predictor) == 0) {
      NA_real_
    } else {
      suppressWarnings(stats::cor(score, predictor, method = "spearman"))
    }
    group_names <- unique(groups)
    this_within <- vector("list", length(group_names))
    for (group_index in seq_along(group_names)) {
      group_name <- group_names[[group_index]]
      use <- groups == group_name
      fit <- .entropy_depth_within_fit(score[use], predictor[use])
      within_index <- within_index + 1L
      this_within[[group_index]] <- data.frame(
        state = state,
        q = q,
        group = group_name,
        covariate = covariate_name,
        n_cells = sum(use),
        within_slope = fit$slope,
        within_slope_std = fit$slope_std,
        within_r_squared = fit$r_squared,
        stringsAsFactors = FALSE
      )
      within_rows[[within_index]] <- this_within[[group_index]]
    }
    this_within <- do.call(rbind, this_within)
    finite_within <- this_within$within_r_squared[
      is.finite(this_within$within_r_squared)
    ]
    within_median <- if (length(finite_within)) {
      stats::median(finite_within)
    } else {
      NA_real_
    }
    summary_rows[[index]] <- data.frame(
      state = state,
      q = q,
      covariate = covariate_name,
      thinnable = grepl(detected_pattern, covariate_name) ||
        grepl(count_pattern, covariate_name),
      r_squared = r_squared,
      spearman_rho = rho,
      within_median = within_median,
      joint_r_squared = joint_r_squared,
      stringsAsFactors = FALSE
    )
  }
  within_rows <- within_rows[seq_len(within_index)]
  list(
    summary = do.call(rbind, summary_rows),
    within = do.call(rbind, within_rows)
  )
}


.entropy_depth_within_fit <- function(response, predictor) {
  empty <- list(slope = NA_real_, slope_std = NA_real_, r_squared = NA_real_)
  if (length(response) < 3L || stats::sd(response) == 0 ||
      stats::sd(predictor) == 0) {
    return(empty)
  }
  fit <- stats::lm(response ~ predictor)
  standardized <- stats::lm(scale(response) ~ scale(predictor))
  list(
    slope = unname(stats::coef(fit)[[2L]]),
    slope_std = unname(stats::coef(standardized)[[2L]]),
    r_squared = unname(summary(fit)$r.squared)
  )
}


.entropy_depth_q_label <- function(value) {
  sub("[.]?0+$", "", formatC(value, format = "f", digits = 3))
}


.entropy_depth_percent <- function(value, digits = 0L) {
  if (!is.finite(value)) return("not available")
  paste0(formatC(100 * value, format = "f", digits = digits), "%")
}


.entropy_depth_broad_line <- function(x) {
  frontier <- x$frontier
  selected <- match(x$decision$q, frontier$q)
  if (is.na(selected)) return(character())
  neighbours <- c(selected - 1L, selected + 1L)
  neighbours <- neighbours[neighbours >= 1L & neighbours <= nrow(frontier)]
  if (!length(neighbours)) return(character())
  best <- frontier$objective[[selected]]
  if (!is.finite(best) || best == 0) return(character())
  close <- neighbours[
    is.finite(frontier$objective[neighbours]) &
      abs(frontier$objective[neighbours] - best) / abs(best) <= 0.02
  ]
  if (!length(close)) return(character())
  labels <- vapply(frontier$q[close], .entropy_depth_q_label, character(1L))
  compared <- if (length(labels) == 1L) {
    paste0("q = ", labels)
  } else {
    paste0("q = ", paste(labels[-length(labels)], collapse = ", q = "),
           " and q = ", labels[[length(labels)]])
  }
  sprintf(
    "  The optimum is broad: %s score within 2%% of q = %s.",
    compared, .entropy_depth_q_label(x$decision$q)
  )
}


.entropy_depth_verdict_lines <- function(x) {
  decision <- x$decision
  frontier <- x$frontier
  selected <- frontier[frontier$q == decision$q, , drop = FALSE][1L, ]
  budget <- 1 / decision$w
  broad <- .entropy_depth_broad_line(x)
  if (isTRUE(decision$correct)) {
    cells_percent <- decision$n_cells_thinned / decision$n_cells
    price <- selected$price_cumulative
    lines <- c(
      "Depth correction: RECOMMENDED",
      sprintf(
        "  retention weight      %s   (a step must cost less than %s per unit)",
        format(decision$w, trim = TRUE),
        formatC(budget, format = "f", digits = 3)
      ),
      sprintf("  recommended q         %s",
              .entropy_depth_q_label(decision$q)),
      sprintf("  target depth          %s counts",
              format(decision$target_depth, big.mark = ",", scientific = FALSE)),
      sprintf(
        "  cells thinned         %s of %s (%s)",
        format(decision$n_cells_thinned, big.mark = ",", scientific = FALSE),
        format(decision$n_cells, big.mark = ",", scientific = FALSE),
        .entropy_depth_percent(cells_percent)
      ),
      sprintf("  confound removed      %s",
              .entropy_depth_percent(selected$confound_removed)),
      sprintf("  variation retained    %s",
              .entropy_depth_percent(selected$retained)),
      sprintf("  average price paid    %s",
              formatC(price, format = "f", digits = 3)),
      broad
    )
    return(lines[nzchar(lines)])
  }

  finite_prices <- frontier$price_step[is.finite(frontier$price_step)]
  price_text <- if (length(finite_prices)) {
    sprintf(
      "The cheapest step costs %s per unit of confound removed, against a budget of %s.",
      formatC(min(finite_prices), format = "f", digits = 3),
      formatC(budget, format = "f", digits = 3)
    )
  } else {
    sprintf(
      "No step has a finite price, against a budget of %s per unit of confound removed.",
      formatC(budget, format = "f", digits = 3)
    )
  }
  lines <- c(
    "Depth correction: NOT RECOMMENDED",
    sprintf(
      "  No q beats leaving the data alone at w = %s. %s",
      format(decision$w, trim = TRUE), price_text
    ),
    "  The frontier is close to diagonal: thinning shrinks this score rather",
    "  than restructuring it. A common cause is depth and biological state",
    "  being collinear, for example when the deep cells are sorted progenitors.",
    "  Use the uncorrected score, or lower w if you accept the trade.",
    broad
  )
  lines[nzchar(lines)]
}


.entropy_depth_association_plot <- function(summary, title, worsened = character()) {
  summary <- summary[order(summary$thinnable, summary$r_squared,
                           decreasing = TRUE), , drop = FALSE]
  summary$section <- ifelse(
    summary$thinnable, "Thinning acts on this", "Thinning does not act on this"
  )
  labels <- summary$covariate
  labels[labels %in% worsened] <- paste0(labels[labels %in% worsened], "  * worse")
  summary$display <- factor(labels, levels = rev(labels))
  maximum <- max(c(summary$r_squared, summary$within_median),
                 na.rm = TRUE)
  if (!is.finite(maximum)) maximum <- 0
  text_x <- max(0.04, maximum * 1.12)
  limit <- max(0.08, maximum * 1.42)
  summary$rho_label <- ifelse(
    is.finite(summary$spearman_rho),
    sprintf("rho %.2f", summary$spearman_rho), "rho NA"
  )
  ggplot2::ggplot(summary, ggplot2::aes(y = display)) +
    ggplot2::geom_col(
      ggplot2::aes(x = r_squared), width = 0.62, fill = "#4C78A8",
      na.rm = TRUE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = within_median), shape = 124, size = 5,
      colour = "#E45756", na.rm = TRUE
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = text_x, label = rho_label), hjust = 0,
      size = 3, colour = "#333333"
    ) +
    ggplot2::facet_grid(
      section ~ ., scales = "free_y", space = "free_y", switch = "y"
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, limit), expand = ggplot2::expansion(mult = c(0, 0.01))
    ) +
    ggplot2::labs(
      title = title,
      subtitle = "Blue bar: global R2. Red tick: median within-group R2.",
      x = "Variance explained (R2)", y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      strip.placement = "outside",
      strip.text.y.left = ggplot2::element_text(angle = 0, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 8),
      plot.title = ggplot2::element_text(face = "bold")
    )
}


.entropy_depth_hex_data <- function(x, y, bins = 24L) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (!length(x) || diff(range(x)) == 0 || diff(range(y)) == 0) return(NULL)
  x_range <- range(x)
  y_range <- range(y)
  scaled_x <- (x - x_range[[1L]]) / diff(x_range) * bins
  scaled_y <- (y - y_range[[1L]]) / diff(y_range) * bins
  row_step <- sqrt(3) / 2
  base_row <- floor(scaled_y / row_step)
  candidates <- vapply(-1:2, function(offset) base_row + offset,
                       numeric(length(base_row)))
  columns <- round(scaled_x - 0.5 * (candidates %% 2))
  centers_x <- columns + 0.5 * (candidates %% 2)
  centers_y <- candidates * row_step
  distances <- (centers_x - scaled_x)^2 + (centers_y - scaled_y)^2
  chosen <- max.col(-distances, ties.method = "first")
  position <- cbind(seq_along(chosen), chosen)
  row_id <- candidates[position]
  column_id <- columns[position]
  keys <- paste(row_id, column_id, sep = ":")
  counts <- table(keys)
  unique_keys <- names(counts)
  pieces <- strsplit(unique_keys, ":", fixed = TRUE)
  rows <- as.numeric(vapply(pieces, `[[`, character(1L), 1L))
  columns <- as.numeric(vapply(pieces, `[[`, character(1L), 2L))
  center_x <- columns + 0.5 * (rows %% 2)
  center_y <- rows * row_step
  angles <- seq(30, 330, by = 60) * pi / 180
  radius <- 1 / sqrt(3)
  polygon <- do.call(rbind, lapply(seq_along(unique_keys), function(index) {
    data.frame(
      group = unique_keys[[index]],
      count = as.numeric(counts[[index]]),
      sx = center_x[[index]] + radius * cos(angles),
      sy = center_y[[index]] + radius * sin(angles),
      stringsAsFactors = FALSE
    )
  }))
  polygon$x <- x_range[[1L]] + polygon$sx / bins * diff(x_range)
  polygon$y <- y_range[[1L]] + polygon$sy / bins * diff(y_range)
  polygon
}


.entropy_depth_hex_plot <- function(x, y, x_name, title) {
  points <- data.frame(x = x, y = y)
  hexagons <- .entropy_depth_hex_data(x, y)
  if (nrow(points) > 10000L) {
    positions <- unique(round(seq(1, nrow(points), length.out = 10000L)))
    smooth_points <- points[positions, , drop = FALSE]
  } else {
    smooth_points <- points
  }
  plot <- ggplot2::ggplot()
  if (is.null(hexagons)) {
    plot <- plot + ggplot2::geom_point(
      data = points, ggplot2::aes(x = x, y = y), alpha = 0.3, size = 0.7
    )
  } else {
    plot <- plot +
      ggplot2::geom_polygon(
        data = hexagons,
        ggplot2::aes(x = x, y = y, group = group, fill = count),
        colour = NA
      ) +
      ggplot2::scale_fill_viridis_c(option = "C", trans = "sqrt",
                                   name = "Cells")
  }
  if (stats::sd(x) > 0 && stats::sd(y) > 0) {
    plot <- plot + ggplot2::geom_smooth(
      data = smooth_points, ggplot2::aes(x = x, y = y),
      method = "loess", formula = y ~ x, se = FALSE,
      colour = "#E45756", linewidth = 0.9
    )
  }
  plot +
    ggplot2::labs(title = title, x = x_name, y = "Uncorrected entropy") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold")
    )
}


.entropy_depth_report_page1 <- function(x) {
  summary <- x$comparison[x$comparison$state == "uncorrected", , drop = FALSE]
  association <- .entropy_depth_association_plot(
    summary, "A. Associations before correction"
  )
  finite_r2 <- which(is.finite(summary$r_squared))
  top <- if (length(finite_r2)) {
    finite_r2[[which.max(summary$r_squared[finite_r2])]]
  } else {
    1L
  }
  top_covariate <- summary$covariate[[top]]
  relationship <- .entropy_depth_hex_plot(
    x$cell_data[[top_covariate]], x$cell_data$entropy_uncorrected,
    top_covariate, paste0("B. Entropy against ", top_covariate)
  )
  joint <- summary$joint_r_squared[[1L]]
  subtitle <- sprintf(
    "Joint R2: %s. %s cells in %s groups.",
    if (is.finite(joint)) formatC(joint, format = "f", digits = 3) else "NA",
    format(nrow(x$cell_data), big.mark = ",", scientific = FALSE),
    length(unique(x$cell_data$group))
  )
  (association | relationship) +
    patchwork::plot_annotation(
      title = "Is there a confound, and can thinning touch it?",
      subtitle = subtitle,
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(size = 16, face = "bold"),
        plot.subtitle = ggplot2::element_text(size = 10)
      )
    )
}


.entropy_depth_report_page2 <- function(x) {
  frontier <- x$frontier
  chosen <- x$decision$q
  long <- rbind(
    data.frame(q = frontier$q, measure = "Confound removed",
               value = frontier$confound_removed),
    data.frame(q = frontier$q, measure = "Variation retained",
               value = frontier$retained)
  )
  colours <- c("Confound removed" = "#E45756",
               "Variation retained" = "#4C78A8")
  first <- ggplot2::ggplot(
    long, ggplot2::aes(x = q, y = value, colour = measure)
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::geom_vline(xintercept = chosen, linetype = 2) +
    ggplot2::scale_colour_manual(values = colours) +
    ggplot2::scale_x_reverse() +
    ggplot2::labs(title = "A. The two numbers", x = "Depth quantile q",
                  y = "Fraction", colour = NULL) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(legend.position = "bottom",
                   plot.title = ggplot2::element_text(face = "bold"))

  second <- ggplot2::ggplot(
    frontier, ggplot2::aes(x = q, y = price_step)
  ) +
    ggplot2::geom_hline(yintercept = 1 / x$decision$w,
                       linetype = 2, colour = "#E45756")
  if (sum(is.finite(frontier$price_step)) > 1L) {
    second <- second + ggplot2::geom_line(linewidth = 0.8, na.rm = TRUE)
  }
  second <- second +
    ggplot2::geom_point(size = 1.8, na.rm = TRUE) +
    ggplot2::geom_vline(xintercept = chosen, linetype = 2) +
    ggplot2::scale_x_reverse() +
    ggplot2::labs(title = "B. Step price", x = "Depth quantile q",
                  y = "Variation lost per confound removed") +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

  third <- ggplot2::ggplot(
    frontier,
    ggplot2::aes(x = confound_removed, y = retained, label =
                   vapply(q, .entropy_depth_q_label, character(1L)))
  ) +
    ggplot2::geom_abline(intercept = 1, slope = -1, linetype = 2,
                        colour = "#777777") +
    ggplot2::geom_path(linewidth = 0.8, colour = "#4C78A8") +
    ggplot2::geom_point(size = 2, colour = "#4C78A8") +
    ggplot2::geom_point(
      data = frontier[frontier$q == chosen, , drop = FALSE],
      size = 4, shape = 21, fill = "#E45756", colour = "black"
    ) +
    ggplot2::geom_text(nudge_y = 0.035, size = 2.6, check_overlap = TRUE) +
    ggplot2::labs(title = "C. Frontier and do-nothing diagonal",
                  x = "Confound removed", y = "Variation retained") +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

  verdict <- paste(.entropy_depth_verdict_lines(x), collapse = "\n")
  text_plot <- ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 1, label = verdict,
                      hjust = 0, vjust = 1, family = "mono", size = 3.2) +
    ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) + ggplot2::theme_void()
  ((first | second | third) / text_plot) +
    patchwork::plot_layout(heights = c(3, 1.25)) +
    patchwork::plot_annotation(
      title = "The price decision",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(size = 16, face = "bold")
      )
    )
}


.entropy_depth_report_page3 <- function(x) {
  baseline <- x$comparison[x$comparison$state == "uncorrected", , drop = FALSE]
  selected <- x$comparison[x$comparison$state == "selected", , drop = FALSE]
  selected <- selected[match(baseline$covariate, selected$covariate), , drop = FALSE]
  worse <- baseline$covariate[
    (is.finite(selected$r_squared) & is.finite(baseline$r_squared) &
       selected$r_squared > baseline$r_squared) |
      (is.finite(selected$within_median) & is.finite(baseline$within_median) &
         selected$within_median > baseline$within_median)
  ]
  left <- .entropy_depth_association_plot(baseline, "A. Before, q = 1")
  right <- .entropy_depth_association_plot(
    selected,
    paste0("B. After, q = ", .entropy_depth_q_label(x$decision$q)),
    worsened = worse
  )
  (left | right) + patchwork::plot_annotation(
    title = "Before and after",
    subtitle = "An asterisk marks a covariate whose global or within-group R2 increased.",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(size = 16, face = "bold")
    )
  )
}


.entropy_depth_report_page4 <- function(x) {
  within <- x$within_comparison
  baseline <- within[within$state == "uncorrected", , drop = FALSE]
  covariate_summary <- x$comparison[
    x$comparison$state == "uncorrected", , drop = FALSE
  ]
  finite <- which(is.finite(covariate_summary$r_squared))
  top_covariate <- if (length(finite)) {
    covariate_summary$covariate[[finite[[which.max(
      covariate_summary$r_squared[finite]
    )]]]]
  } else {
    covariate_summary$covariate[[1L]]
  }
  top <- baseline[baseline$covariate == top_covariate, , drop = FALSE]
  top <- top[order(top$within_r_squared, decreasing = TRUE, na.last = TRUE), ]
  group_levels <- unique(c(as.character(top$group), as.character(baseline$group)))
  within$group <- factor(within$group, levels = rev(group_levels))
  within$covariate <- factor(
    within$covariate, levels = unique(covariate_summary$covariate)
  )
  within$panel <- ifelse(
    within$state == "uncorrected", "q = 1",
    paste0("q = ", .entropy_depth_q_label(x$decision$q))
  )
  ggplot2::ggplot(
    within,
    ggplot2::aes(x = covariate, y = group, fill = within_r_squared)
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.2) +
    ggplot2::facet_grid(. ~ panel) +
    ggplot2::scale_fill_viridis_c(option = "C", na.value = "grey90",
                                 name = "Within-group R2") +
    ggplot2::labs(
      title = "Per-group associations before and after",
      subtitle = paste0("Groups are sorted by ", top_covariate,
                        ", highest within-group R2 first."),
      x = NULL, y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title = ggplot2::element_text(size = 16, face = "bold"),
      strip.text = ggplot2::element_text(face = "bold")
    )
}


.entropy_depth_report_pages <- function(x, per_cluster) {
  pages <- list(
    .entropy_depth_report_page1(x),
    .entropy_depth_report_page2(x)
  )
  if (isTRUE(x$decision$correct)) {
    pages[[length(pages) + 1L]] <- .entropy_depth_report_page3(x)
    if (isTRUE(per_cluster)) {
      pages[[length(pages) + 1L]] <- .entropy_depth_report_page4(x)
    }
  }
  pages
}


.entropy_depth_write_pdf <- function(path, pages) {
  grDevices::pdf(path, width = 11, height = 8.5, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (page in pages) print(page)
  invisible(path)
}


utils::globalVariables(c(
  "display", "r_squared", "within_median", "rho_label", "section",
  "measure", "value", "q", "price_step", "confound_removed", "retained",
  "count", "group", "within_r_squared", "covariate", "panel", "x", "y",
  "cluster", "label"
))
