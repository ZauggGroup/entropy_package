#' Recompute entropy after a targeted intervention
#'
#' `entropy_recompute_check()` changes one part of the input, calls
#' [run_entropy()] again, and asks whether the cluster ranking of entropy
#' survives. The available interventions equalise count depth, equalise the
#' number of sampled cells per group, or remove a specified gene programme.
#'
#' @details
#' This is an interventional check. It removes an effect physically rather than
#' adjusting for a covariate statistically, then recomputes entropy. A ranking
#' that survives is therefore stronger evidence, although it still does not
#' prove plasticity or distinguish every technical and biological cause.
#'
#' Fixed neighbours are the default for count downsampling and gene removal.
#' This isolates what the entropy estimator reads from a fixed local context.
#' For gene removal, `recompute_neighbors = TRUE` instead asks whether the
#' removed programme also shaped the embedding. Cell subsampling always
#' recomputes neighbours because removed cells cannot remain neighbours.
#'
#' When the first argument is an `entropy_diagnostic`, supply the expression
#' object in `expression_seu` and name the follow-up covariate with `covariate`.
#' The entropy vector, group name, and default mechanism are then read from the
#' diagnostic.
#'
#' @param seu A Seurat object, or an `entropy_diagnostic` object when
#'   `expression_seu` is supplied.
#' @param entropy Original finite per-cell entropy vector in Seurat cell order.
#'   It is read from the diagnostic when `seu` is an `entropy_diagnostic`.
#' @param group Metadata column defining clusters. It is read from the
#'   diagnostic when possible.
#' @param mechanism One of `"downsample_counts"`, `"subsample_cells"`, or
#'   `"remove_genes"`.
#' @param entropy_args Named list of arguments originally passed to
#'   [run_entropy()]. `assay` is required. Output arguments are ignored because
#'   the recomputed vector is captured directly.
#' @param genes Character vector of genes to remove.
#' @param depth_quantile Quantile used to choose the downsampling target.
#' @param target_depth Optional explicit integer target depth.
#' @param drop_below_target Logical. Drop cells below the target instead of
#'   retaining their original counts. Dropping cells forces
#'   `recompute_neighbors = TRUE`, because a dropped cell cannot remain a
#'   neighbour of the cells that stayed. The embedding is not rebuilt, so the
#'   graph is rebuilt from the existing reduction on the retained cells.
#' @param cells_per_group Number of cells retained per group. The default is
#'   the size of the smallest group.
#' @param recompute_neighbors Logical or `NULL`. `NULL` uses the mechanism's
#'   deliberate default.
#' @param n_repeats Positive number of random repeats. Gene removal is forced
#'   to one repeat because it is deterministic.
#' @param seed Integer base seed. Repeat `i` uses `seed + i - 1`.
#' @param top_fraction Fraction of clusters defining the high-entropy set.
#' @param rank_tolerance Non-negative number of rank positions allowed for a
#'   `"survives"` verdict. A cluster that moves down by more than this value in
#'   every repeat `"collapses"`; one that moves up by more than this value in
#'   every repeat `"improves"`.
#' @param covariate Covariate to follow up when `seu` is an
#'   `entropy_diagnostic`.
#' @param expression_seu Seurat expression object used with a diagnostic first
#'   argument.
#' @param verbose Logical. Print cost and progress messages.
#'
#' @return An object of class `entropy_recompute` containing per-repeat cluster,
#'   cell, rank, and optional density comparisons.
#' @export
entropy_recompute_check <- function(
    seu,
    entropy = NULL,
    group = NULL,
    mechanism = c("downsample_counts", "subsample_cells", "remove_genes"),
    entropy_args = list(),
    genes = NULL,
    depth_quantile = 0.1,
    target_depth = NULL,
    drop_below_target = FALSE,
    cells_per_group = NULL,
    recompute_neighbors = NULL,
    n_repeats = 3,
    seed = 1L,
    top_fraction = 1 / 3,
    rank_tolerance = 2,
    covariate = NULL,
    expression_seu = NULL,
    verbose = TRUE) {

  mechanism_missing <- missing(mechanism)
  diagnostic <- NULL
  if (inherits(seu, "entropy_diagnostic")) {
    diagnostic <- seu
    if (is.null(expression_seu) || !inherits(expression_seu, "Seurat")) {
      stop(
        "When the first argument is an entropy_diagnostic, `expression_seu` must be the Seurat expression object.",
        call. = FALSE
      )
    }
    if (length(covariate) != 1L || !is.character(covariate) ||
        is.na(covariate) || !nzchar(covariate)) {
      stop("`covariate` must name one diagnostic covariate.", call. = FALSE)
    }
    info <- diagnostic$covariate_info[
      diagnostic$covariate_info$covariate == covariate, , drop = FALSE
    ]
    if (nrow(info) != 1L) {
      stop(sprintf("Covariate '%s' is not present in the diagnostic.", covariate),
           call. = FALSE)
    }
    if (mechanism_missing) {
      mechanism <- info$intervention[[1L]]
      if (identical(mechanism, "none")) {
        stop(sprintf(
          paste0(
            "No recompute test exists for covariate '%s'. ",
            "The observational result is all the data can give."
          ),
          covariate
        ), call. = FALSE)
      }
    }
    seu <- expression_seu
    diagnostic_cells <- as.character(diagnostic$cell_data$cell)
    cell_index <- match(colnames(seu), diagnostic_cells)
    if (anyNA(cell_index) || length(diagnostic_cells) != ncol(seu)) {
      stop("The diagnostic cells do not match `expression_seu`.", call. = FALSE)
    }
    entropy <- diagnostic$cell_data$entropy[cell_index]
    if (is.null(group)) {
      group <- diagnostic$provenance$group_column
    }
  }

  if (!inherits(seu, "Seurat")) {
    stop("`seu` must be a Seurat object.", call. = FALSE)
  }
  mechanism <- match.arg(mechanism)
  if (!is.list(entropy_args) ||
      (length(entropy_args) &&
       (is.null(names(entropy_args)) || any(!nzchar(names(entropy_args))) ||
        anyDuplicated(names(entropy_args))))) {
    stop("`entropy_args` must be a named list.", call. = FALSE)
  }
  if (!"assay" %in% names(entropy_args) ||
      length(entropy_args$assay) != 1L || !is.character(entropy_args$assay) ||
      is.na(entropy_args$assay) || !nzchar(entropy_args$assay)) {
    stop("`entropy_args` must include one non-empty `assay` value.", call. = FALSE)
  }
  assay <- entropy_args$assay
  if (!assay %in% Seurat::Assays(seu)) {
    stop(sprintf("Assay '%s' from `entropy_args` is not in `seu`.", assay),
         call. = FALSE)
  }
  cell_names <- colnames(seu)
  if (!is.numeric(entropy) || length(entropy) != length(cell_names) ||
      any(!is.finite(entropy))) {
    stop("`entropy` must contain one finite numeric value per Seurat cell.",
         call. = FALSE)
  }
  entropy <- as.numeric(entropy)
  names(entropy) <- cell_names
  if (length(group) != 1L || !is.character(group) || is.na(group) ||
      !nzchar(group) || !group %in% colnames(seu[[]])) {
    stop("`group` must name one metadata column in `seu`.", call. = FALSE)
  }
  group_values <- seu[[]][[group]]
  if (anyNA(group_values)) {
    stop("The grouping column must not contain missing values.", call. = FALSE)
  }
  group_values <- as.character(group_values)
  names(group_values) <- cell_names

  n_repeats <- .entropy_recompute_positive_integer(n_repeats, "n_repeats")
  seed <- .entropy_recompute_integer(seed, "seed")
  top_fraction <- .entropy_recompute_number(
    top_fraction, "top_fraction", lower = 0, upper = 1, lower_open = TRUE
  )
  rank_tolerance <- .entropy_recompute_number(
    rank_tolerance, "rank_tolerance", lower = 0
  )
  depth_quantile <- .entropy_recompute_number(
    depth_quantile, "depth_quantile", lower = 0, upper = 1
  )
  .entropy_scalar_logical(drop_below_target, "drop_below_target")
  .entropy_scalar_logical(verbose, "verbose")
  if (!is.null(recompute_neighbors)) {
    .entropy_scalar_logical(recompute_neighbors, "recompute_neighbors")
  }

  if (mechanism == "subsample_cells") {
    if (identical(recompute_neighbors, FALSE)) {
      stop(
        "`recompute_neighbors = FALSE` is not valid for `subsample_cells`; removed cells cannot remain neighbours.",
        call. = FALSE
      )
    }
    recompute_neighbors <- TRUE
  } else if (mechanism == "downsample_counts" && isTRUE(drop_below_target)) {
    # Dropping cells below the target IS cell subsampling: the removed cells
    # cannot remain neighbours of the cells that stayed. Keeping the original
    # graph and intersecting it leaves each cell with a different number of
    # neighbours, which makes entropy values incomparable between cells.
    if (identical(recompute_neighbors, FALSE)) {
      stop(paste0(
        "`recompute_neighbors = FALSE` is not valid with ",
        "`drop_below_target = TRUE`; dropped cells cannot remain neighbours."),
        call. = FALSE)
    }
    recompute_neighbors <- TRUE
  } else if (is.null(recompute_neighbors)) {
    recompute_neighbors <- FALSE
  }
  if (mechanism == "remove_genes" && n_repeats != 1L) {
    if (verbose) message("Gene removal is deterministic, so `n_repeats` is set to 1.")
    n_repeats <- 1L
  }
  warning_messages <- character()
  record_warning <- function(message) {
    warning_messages <<- unique(c(warning_messages, message))
    warning(message, call. = FALSE)
  }
  if (ncol(seu) > 20000L && !identical(entropy_args$chunk_method, "auto")) {
    record_warning(
      "More than 20,000 cells are present and `chunk_method = \"auto\"` is not set. This recompute may be slow."
    )
  }
  if (verbose) {
    message(sprintf(
      "Entropy recomputation will process %d cells across %d repeat%s. This can take a long time.",
      ncol(seu), n_repeats, if (n_repeats == 1L) "" else "s"
    ))
  }

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  target <- NULL
  n_below_target <- NA_integer_
  used_genes <- NULL
  missing_genes <- character()
  n_genes_removed <- 0L
  intervention_totals <- NULL
  group_sizes <- table(group_values)
  if (mechanism == "downsample_counts") {
    counts <- SeuratObject::GetAssayData(seu, assay = assay, layer = "counts")
    totals <- Matrix::colSums(counts)
    if (is.null(target_depth)) {
      target <- floor(as.numeric(stats::quantile(
        totals, probs = depth_quantile, names = FALSE, type = 7
      )))
    } else {
      target <- .entropy_recompute_integer(target_depth, "target_depth", lower = 1)
    }
    if (!is.finite(target) || target < 1L) {
      stop("The downsampling target must be at least 1 count.", call. = FALSE)
    }
    n_below_target <- sum(totals < target)
    if (n_below_target > 0.1 * length(totals)) {
      record_warning(sprintf(
        "%d cells (%.1f%%) are below the target depth and cannot be downsampled.",
        n_below_target, 100 * n_below_target / length(totals)
      ))
    }
    if (drop_below_target && n_below_target > 0L) {
      record_warning(sprintf(
        "Dropping %d cells below the target depth changes the cell composition.",
        n_below_target
      ))
    }
  } else if (mechanism == "subsample_cells") {
    if (is.null(cells_per_group)) {
      cells_per_group <- min(group_sizes)
    } else {
      cells_per_group <- .entropy_recompute_positive_integer(
        cells_per_group, "cells_per_group"
      )
    }
    small <- names(group_sizes)[group_sizes < cells_per_group]
    if (length(small)) {
      record_warning(sprintf(
        "%d group%s contain fewer than `cells_per_group` cells and will be kept whole: %s.",
        length(small), if (length(small) == 1L) "" else "s",
        paste(small, collapse = ", ")
      ))
    }
  } else {
    if (!is.character(genes) || !length(genes) || anyNA(genes) ||
        any(!nzchar(genes))) {
      stop("`genes` must contain at least one non-empty gene for `remove_genes`.",
           call. = FALSE)
    }
    assay_genes <- rownames(seu[[assay]])
    genes <- unique(genes)
    present <- intersect(genes, assay_genes)
    missing_genes <- setdiff(genes, assay_genes)
    if (!length(present)) {
      stop(sprintf("None of the requested genes are present in assay '%s'.", assay),
           call. = FALSE)
    }
    if (length(missing_genes)) {
      record_warning(sprintf(
        "%d requested gene%s not found in assay '%s'.",
        length(missing_genes), if (length(missing_genes) == 1L) " was" else "s were",
        assay
      ))
    }
    used_genes <- setdiff(assay_genes, present)
    n_genes_removed <- length(present)
  }

  fixed_neighbors <- NULL
  if (!recompute_neighbors) {
    fixed_neighbors <- .entropy_recompute_neighbors(seu, entropy_args)
  }

  cluster_parts <- vector("list", n_repeats)
  repeat_parts <- vector("list", n_repeats)
  cell_parts <- vector("list", n_repeats)
  density_parts <- if (mechanism == "subsample_cells") vector("list", n_repeats) else NULL
  total_parts <- if (mechanism == "downsample_counts") vector("list", n_repeats) else NULL

  for (repeat_id in seq_len(n_repeats)) {
    set.seed(seed + repeat_id - 1L)
    if (verbose) message(sprintf("Recompute repeat %d of %d.", repeat_id, n_repeats))
    repeat_seu <- seu
    retained <- cell_names
    repeat_features <- NULL

    if (mechanism == "downsample_counts") {
      counts <- SeuratObject::GetAssayData(seu, assay = assay, layer = "counts")
      totals_before <- Matrix::colSums(counts)
      if (drop_below_target) {
        retained <- names(totals_before)[totals_before >= target]
        repeat_seu <- seu[, retained]
        counts <- counts[, retained, drop = FALSE]
        totals_before <- totals_before[retained]
      }
      downsampled <- .entropy_downsample_matrix(counts, target)
      repeat_seu[[assay]] <- SeuratObject::SetAssayData(
        repeat_seu[[assay]], layer = "counts", new.data = downsampled
      )
      total_parts[[repeat_id]] <- data.frame(
        cell = retained,
        repeat_id = repeat_id,
        total_before = as.numeric(totals_before),
        total_after = as.numeric(Matrix::colSums(downsampled)),
        stringsAsFactors = FALSE
      )
    } else if (mechanism == "subsample_cells") {
      retained <- unlist(lapply(names(group_sizes), function(level) {
        candidates <- names(group_values)[group_values == level]
        if (length(candidates) <= cells_per_group) candidates else
          sample(candidates, cells_per_group, replace = FALSE)
      }), use.names = FALSE)
      retained <- cell_names[cell_names %in% retained]
      repeat_seu <- seu[, retained]
    } else {
      repeat_features <- used_genes
    }

    repeat_neighbors <- NULL
    if (!recompute_neighbors) {
      repeat_neighbors <- .entropy_subset_neighbors(fixed_neighbors, retained)
    }
    recomputed <- .entropy_run_recompute(
      repeat_seu,
      entropy_args = entropy_args,
      neighbors = repeat_neighbors,
      recompute_neighbors = recompute_neighbors,
      features = repeat_features,
      rebuild_reduction = recompute_neighbors &&
        !(mechanism %in% c("subsample_cells", "downsample_counts"))
    )
    entropy_before <- entropy[retained]
    entropy_after <- recomputed[retained]
    if (anyNA(entropy_after)) {
      stop("The recomputed entropy vector did not contain every retained cell.",
           call. = FALSE)
    }
    repeat_groups <- group_values[retained]
    comparison <- .entropy_recompute_comparison(
      entropy_before, entropy_after, repeat_groups,
      repeat_id = repeat_id, top_fraction = top_fraction,
      delta_comparable = mechanism != "remove_genes"
    )
    cluster_parts[[repeat_id]] <- comparison$cluster_summary
    repeat_parts[[repeat_id]] <- comparison$repeat_summary
    cell_parts[[repeat_id]] <- data.frame(
      cell = retained,
      cluster = unname(repeat_groups),
      entropy_before = unname(entropy_before),
      entropy_after = unname(entropy_after),
      repeat_id = repeat_id,
      stringsAsFactors = FALSE
    )

    if (mechanism == "subsample_cells") {
      before_distance <- .entropy_recompute_mean_knn(seu, entropy_args)
      after_distance <- .entropy_recompute_mean_knn(repeat_seu, entropy_args)
      density_parts[[repeat_id]] <- do.call(rbind, lapply(
        unique(unname(repeat_groups)), function(cluster_name) {
          cells <- retained[repeat_groups == cluster_name]
          data.frame(
            cluster = cluster_name,
            repeat_id = repeat_id,
            mean_knn_before = mean(before_distance[cells]),
            mean_knn_after = mean(after_distance[cells]),
            stringsAsFactors = FALSE
          )
        }
      ))
    }
  }

  cluster_summary <- do.call(rbind, cluster_parts)
  rownames(cluster_summary) <- NULL
  repeat_summary <- do.call(rbind, repeat_parts)
  rownames(repeat_summary) <- NULL
  cell_data <- do.call(rbind, cell_parts)
  rownames(cell_data) <- NULL
  density_check <- if (is.null(density_parts)) NULL else do.call(rbind, density_parts)
  intervention_totals <- if (is.null(total_parts)) NULL else do.call(rbind, total_parts)

  verdict_parts <- lapply(split(cluster_summary$rank_change, cluster_summary$cluster),
                          function(change) {
    verdict <- if (all(abs(change) <= rank_tolerance)) {
      "survives"
    } else if (all(change > rank_tolerance)) {
      "collapses"
    } else if (all(change < -rank_tolerance)) {
      "improves"
    } else {
      "unstable"
    }
    data.frame(
      median_rank_change = stats::median(change),
      min_rank_change = min(change),
      max_rank_change = max(change),
      verdict = verdict,
      stringsAsFactors = FALSE
    )
  })
  cluster_verdict <- do.call(rbind, verdict_parts)
  cluster_verdict$cluster <- rownames(cluster_verdict)
  rownames(cluster_verdict) <- NULL
  cluster_verdict <- cluster_verdict[, c(
    "cluster", "median_rank_change", "min_rank_change", "max_rank_change", "verdict"
  )]
  rank_order <- stats::aggregate(
    cluster_summary$rank_before,
    list(cluster = cluster_summary$cluster), stats::median
  )
  cluster_verdict <- cluster_verdict[
    order(rank_order$x[match(cluster_verdict$cluster, rank_order$cluster)]), , drop = FALSE
  ]

  cells_after <- vapply(cell_parts, nrow, integer(1))
  result <- list(
    provenance = list(
      mechanism = mechanism,
      entropy_args = entropy_args,
      n_repeats = n_repeats,
      seed = seed,
      recompute_neighbors = recompute_neighbors,
      target_depth = target,
      n_below_target = n_below_target,
      drop_below_target = drop_below_target,
      cells_per_group = cells_per_group,
      n_genes_removed = n_genes_removed,
      n_genes_missing = length(missing_genes),
      genes_missing = missing_genes,
      features_used = used_genes,
      n_cells_before = ncol(seu),
      n_cells_after = if (length(unique(cells_after)) == 1L) cells_after[[1L]] else cells_after,
      group = group,
      covariate = covariate,
      absolute_delta_comparable = mechanism != "remove_genes",
      warnings = warning_messages,
      call = match.call()
    ),
    cluster_summary = cluster_summary,
    cluster_verdict = cluster_verdict,
    repeat_summary = repeat_summary,
    cell_data = cell_data,
    density_check = density_check,
    intervention_totals = intervention_totals,
    thresholds = list(
      rank_tolerance = rank_tolerance,
      top_fraction = top_fraction
    )
  )
  class(result) <- "entropy_recompute"
  result
}


#' @export
print.entropy_recompute <- function(x, ...) {
  if (!inherits(x, "entropy_recompute")) {
    stop("`x` must be an entropy_recompute object.", call. = FALSE)
  }
  rank_values <- x$repeat_summary$rank_correlation
  cat("Entropy recompute check\n")
  cat(sprintf("  Mechanism: %s\n", x$provenance$mechanism))
  cat(sprintf("  Repeats: %d\n", x$provenance$n_repeats))
  cat(sprintf(
    "  Cells: %d before, %s after per repeat\n",
    x$provenance$n_cells_before,
    paste(x$provenance$n_cells_after, collapse = ", ")
  ))
  cat(sprintf(
    "  Headline rank correlation: mean %.3f, range %.3f to %.3f\n",
    mean(rank_values), min(rank_values), max(rank_values)
  ))
  counts <- table(factor(
    x$cluster_verdict$verdict,
    levels = c("survives", "collapses", "improves", "unstable")
  ))
  cat(sprintf(
    "  Verdicts: %d survives, %d collapses, %d improves, %d unstable\n",
    counts[["survives"]], counts[["collapses"]], counts[["improves"]],
    counts[["unstable"]]
  ))
  if (length(x$provenance$warnings)) {
    cat("  Warnings recorded:\n")
    for (message in x$provenance$warnings) cat(sprintf("    - %s\n", message))
  }
  cat("  Use entropy_recompute_report() for the full report.\n")
  invisible(x)
}


.entropy_recompute_comparison <- function(
    entropy_before, entropy_after, groups, repeat_id, top_fraction,
    delta_comparable = TRUE) {
  levels <- unique(as.character(groups))
  before <- vapply(levels, function(level) {
    stats::median(entropy_before[groups == level])
  }, numeric(1))
  after <- vapply(levels, function(level) {
    stats::median(entropy_after[groups == level])
  }, numeric(1))
  rank_before <- rank(-before, ties.method = "average")
  rank_after <- rank(-after, ties.method = "average")
  top_n <- max(1L, ceiling(length(levels) * top_fraction))
  before_top <- levels[rank_before <= top_n]
  after_top <- levels[rank_after <= top_n]
  cluster_summary <- data.frame(
    cluster = levels,
    repeat_id = repeat_id,
    n_cells = vapply(levels, function(level) sum(groups == level), integer(1)),
    median_before = unname(before),
    median_after = unname(after),
    delta = if (delta_comparable) unname(after - before) else NA_real_,
    rank_before = unname(rank_before),
    rank_after = unname(rank_after),
    rank_change = unname(rank_after - rank_before),
    stringsAsFactors = FALSE
  )
  repeat_summary <- data.frame(
    repeat_id = repeat_id,
    rank_correlation = suppressWarnings(stats::cor(
      before, after, method = "spearman"
    )),
    cell_correlation = suppressWarnings(stats::cor(
      entropy_before, entropy_after, method = "spearman"
    )),
    top_k_retained = length(intersect(before_top, after_top)) / length(before_top),
    stringsAsFactors = FALSE
  )
  list(cluster_summary = cluster_summary, repeat_summary = repeat_summary)
}


.entropy_run_recompute <- function(
    seu, entropy_args, neighbors, recompute_neighbors, features = NULL,
    rebuild_reduction = FALSE) {
  if (rebuild_reduction) {
    seu <- .entropy_rebuild_reduction(seu, entropy_args, features)
  }
  ignored <- c(
    "seu", "neighbors", "neighbors_method", "add_assay", "meta_name",
    "output_path", "features"
  )
  call_args <- entropy_args[setdiff(names(entropy_args), ignored)]
  call_args <- c(list(seu = seu), call_args)
  if (!is.null(features)) call_args$features <- features
  if (recompute_neighbors) {
    call_args$neighbors_method <- "auto"
  } else {
    call_args$neighbors <- neighbors
    call_args$neighbors_method <- "custom"
  }
  call_args$add_assay <- FALSE
  call_args$meta_name <- NULL
  call_args$output_path <- NULL
  recomputed <- do.call(run_entropy, call_args)
  if (!is.data.frame(recomputed) || ncol(recomputed) < 1L) {
    stop("`run_entropy()` did not return a per-cell entropy table.", call. = FALSE)
  }
  values <- as.numeric(recomputed[[1L]])
  names(values) <- rownames(recomputed)
  values
}


.entropy_rebuild_reduction <- function(seu, entropy_args, features = NULL) {
  assay <- entropy_args$assay
  neighbor_assay <- if (is.null(entropy_args$neighbor_assay)) assay else
    entropy_args$neighbor_assay
  reduction <- if (is.null(entropy_args$reduction)) {
    paste0("pca_", neighbor_assay)
  } else {
    entropy_args$reduction
  }
  dims <- if (is.null(entropy_args$dims)) 20L else entropy_args$dims
  dims <- if (length(dims) == 1L) seq_len(as.integer(dims)) else as.integer(dims)
  rebuild_features <- if (is.null(features)) rownames(seu[[neighbor_assay]]) else
    intersect(features, rownames(seu[[neighbor_assay]]))
  old_assay <- Seurat::DefaultAssay(seu)
  Seurat::DefaultAssay(seu) <- neighbor_assay
  on.exit(Seurat::DefaultAssay(seu) <- old_assay, add = TRUE)
  seu <- Seurat::NormalizeData(seu, verbose = FALSE)
  seu <- Seurat::ScaleData(seu, features = rebuild_features, verbose = FALSE)
  seu <- Seurat::RunPCA(
    seu, features = rebuild_features, npcs = max(dims),
    reduction.name = reduction,
    reduction.key = paste0(gsub("[^A-Za-z0-9]", "", reduction), "_"),
    verbose = FALSE
  )
  seu
}


.entropy_recompute_neighbors <- function(seu, entropy_args) {
  embedding <- .entropy_recompute_embedding(seu, entropy_args)
  k <- if (is.null(entropy_args$k)) 15L else as.integer(entropy_args$k)
  if (nrow(embedding) <= k) {
    stop(sprintf("Cannot build %d fixed neighbours from %d cells.", k, nrow(embedding)),
         call. = FALSE)
  }
  fit <- RANN::nn2(embedding, query = embedding, k = k + 1L)
  cells <- rownames(embedding)
  result <- lapply(seq_len(nrow(embedding)), function(i) {
    index <- fit$nn.idx[i, fit$nn.idx[i, ] != i]
    cells[index[seq_len(k)]]
  })
  names(result) <- cells
  result
}


.entropy_subset_neighbors <- function(neighbors, cells) {
  result <- lapply(neighbors[cells], function(value) intersect(value, cells))
  names(result) <- cells
  result
}


.entropy_recompute_embedding <- function(seu, entropy_args) {
  assay <- entropy_args$assay
  neighbor_assay <- if (is.null(entropy_args$neighbor_assay)) assay else
    entropy_args$neighbor_assay
  reduction <- if (is.null(entropy_args$reduction)) {
    paste0("pca_", neighbor_assay)
  } else {
    entropy_args$reduction
  }
  if (!reduction %in% names(seu@reductions)) {
    stop(sprintf(
      paste0(
        "Reduction '%s' is required to build fixed neighbours. ",
        "Compute it first or set `recompute_neighbors = TRUE`."
      ),
      reduction
    ), call. = FALSE)
  }
  embedding <- Seurat::Embeddings(seu[[reduction]])
  dims <- if (is.null(entropy_args$dims)) 20L else entropy_args$dims
  dims <- if (length(dims) == 1L) seq_len(as.integer(dims)) else as.integer(dims)
  dims <- dims[dims >= 1L & dims <= ncol(embedding)]
  if (!length(dims)) {
    stop(sprintf("No requested dimensions exist in reduction '%s'.", reduction),
         call. = FALSE)
  }
  embedding[, dims, drop = FALSE]
}


.entropy_recompute_mean_knn <- function(seu, entropy_args) {
  embedding <- .entropy_recompute_embedding(seu, entropy_args)
  k <- if (is.null(entropy_args$k)) 15L else as.integer(entropy_args$k)
  k <- min(k, nrow(embedding) - 1L)
  if (k < 1L) return(stats::setNames(rep(NA_real_, nrow(embedding)), rownames(embedding)))
  fit <- RANN::nn2(embedding, query = embedding, k = k + 1L)
  result <- vapply(seq_len(nrow(embedding)), function(i) {
    distance <- fit$nn.dists[i, fit$nn.idx[i, ] != i]
    mean(distance[seq_len(k)])
  }, numeric(1))
  stats::setNames(result, rownames(embedding))
}


# Internal: thin each cell's counts to a target depth
#
# Cells whose total is above `target` are downsampled without replacement to
# exactly `target` counts. Cells at or below the target are returned
# unchanged. `scuttle::downsampleMatrix()` is used when it is installed.
# Otherwise a built-in sampler works directly on the column-compressed slots
# and only visits the cells above the target.
#
# @param counts A non-negative integer-like count matrix, genes by cells.
# @param target One positive count depth.
#
# @return The thinned counts as a column-compressed sparse matrix with the
#   same dimensions and dimnames.
.entropy_downsample_matrix <- function(counts, target) {
  totals <- Matrix::colSums(counts)
  proportions <- pmin(1, target / totals)
  proportions[!is.finite(proportions)] <- 1
  if (requireNamespace("scuttle", quietly = TRUE)) {
    downsampled <- scuttle::downsampleMatrix(
      counts, prop = proportions, bycol = TRUE
    )
    return(.entropy_as_csparse(downsampled))
  }
  futile.logger::flog.warn(paste(
    "Package 'scuttle' is not installed; using the built-in count thinning.",
    "Install it with BiocManager::install(\"scuttle\") for faster depth",
    "correction on large datasets."
  ))
  .entropy_downsample_fallback(counts, target, totals)
}

# Internal: built-in count thinning without scuttle
#
# Reads and writes the `@p` and `@x` slots of a `dgCMatrix` directly, so no
# per-column subsetting or full-matrix validation happens inside the loop.
# Each affected cell's counts are expanded to molecules, `target` molecules
# are drawn without replacement, and the draws are tabulated back onto the
# cell's nonzero genes. Totals of affected cells equal `target` exactly.
#
# @param counts A non-negative integer-like count matrix, genes by cells.
# @param target One positive count depth.
# @param totals Column totals of `counts`.
#
# @return A `dgCMatrix` with the same dimensions and dimnames.
.entropy_downsample_fallback <- function(counts, target,
                                         totals = Matrix::colSums(counts)) {
  counts <- .entropy_as_csparse(counts)
  if (!methods::is(counts, "dgCMatrix")) {
    counts <- methods::as(counts, "dgCMatrix")
  }
  affected <- which(totals > target)
  if (!length(affected)) {
    return(counts)
  }
  target <- as.integer(floor(target))
  pointers <- counts@p
  values <- counts@x
  for (column in affected) {
    from <- pointers[[column]] + 1L
    to <- pointers[[column + 1L]]
    index <- from:to
    cell_counts <- as.integer(round(values[index]))
    molecules <- rep.int(seq_along(cell_counts), cell_counts)
    kept <- molecules[sample.int(length(molecules), target)]
    values[index] <- tabulate(kept, nbins = length(cell_counts))
  }
  counts@x <- values
  Matrix::drop0(counts)
}


.entropy_recompute_positive_integer <- function(value, argument) {
  .entropy_recompute_integer(value, argument, lower = 1)
}


.entropy_recompute_integer <- function(value, argument, lower = -Inf) {
  if (length(value) != 1L || !is.numeric(value) || !is.finite(value) ||
      value != floor(value) || value < lower) {
    stop(sprintf("`%s` must be one integer of at least %s.", argument, lower),
         call. = FALSE)
  }
  as.integer(value)
}


.entropy_recompute_number <- function(
    value, argument, lower = -Inf, upper = Inf, lower_open = FALSE) {
  if (length(value) != 1L || !is.numeric(value) || !is.finite(value)) {
    stop(sprintf("`%s` is outside its allowed range.", argument), call. = FALSE)
  }
  valid_lower <- if (lower_open) value > lower else value >= lower
  if (!valid_lower || value > upper) {
    stop(sprintf("`%s` is outside its allowed range.", argument), call. = FALSE)
  }
  as.numeric(value)
}
#' Internal: coerce a matrix to column-compressed sparse storage
#'
#' Downsampling backends can return classes other than `dgCMatrix`. Newer
#' versions of `scuttle` return an `SVT_SparseMatrix` from the `SparseArray`
#' package. Downstream code indexes the `@p` and `@i` slots directly and
#' silently reshapes such objects, so every downsampled matrix is coerced to
#' column-compressed sparse storage here, at the single point where it is
#' produced.
#'
#' @param x A matrix-like object.
#'
#' @return The same values in `CsparseMatrix` storage, with dimensions and
#'   dimnames preserved.
#' @keywords internal
.entropy_as_csparse <- function(x) {
  if (inherits(x, "CsparseMatrix")) {
    return(x)
  }
  original_dim <- dim(x)
  original_dimnames <- dimnames(x)
  converted <- tryCatch(
    methods::as(x, "CsparseMatrix"),
    error = function(error) {
      Matrix::Matrix(as.matrix(x), sparse = TRUE)
    }
  )
  if (!identical(dim(converted), original_dim)) {
    stop(
      sprintf(
        paste0(
          "Coercing the downsampled counts to sparse storage changed the ",
          "dimensions from %d by %d to %d by %d."
        ),
        original_dim[[1L]], original_dim[[2L]],
        dim(converted)[[1L]], dim(converted)[[2L]]
      ),
      call. = FALSE
    )
  }
  dimnames(converted) <- original_dimnames
  converted
}
