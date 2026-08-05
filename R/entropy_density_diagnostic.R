#' Diagnose whether entropy is associated with neighbourhood density
#'
#' `entropy_density_diagnostic()` tests whether a per-cell entropy or
#' plasticity score is associated with local neighbourhood density in an
#' embedding. It is a diagnostic, not a correction: the returned residuals
#' must not be interpreted as corrected entropy because a genuinely plastic
#' biological state may itself occupy a diffuse region of the embedding.
#'
#' Local sparsity is measured as the mean distance to each cell's `k` nearest
#' neighbours (excluding the cell itself). Larger values therefore mean a
#' sparser neighbourhood. The function reports the global Spearman correlation
#' between entropy and this distance. It always fits density as its own model
#' (`entropy ~ density`). Every supplied covariate is then fitted in a completely
#' independent model (`entropy ~ covariate`): covariates are never combined with
#' density or with one another. Positive residuals indicate entropy above the
#' relevant fitted trend; negative residuals indicate entropy below it.
#'
#' Use the same embedding and `k` that were used to calculate entropy whenever
#' possible. The diagnostic is modality-agnostic and can be used with any
#' numeric cell embedding.
#'
#' @param seu A Seurat object whose cells correspond, in the same order, to
#'   `entropy` and the rows of `embedding`. When `group` is supplied, it is
#'   looked up directly in this object's metadata.
#' @param entropy Numeric vector with one entropy/plasticity value per cell.
#' @param embedding Numeric matrix (cells by dimensions), in the same cell order
#'   as `entropy`.
#' @param k Positive integer giving the number of nearest neighbours. Must be
#'   smaller than the number of cells. Default is 30.
#' @param group Optional single character string naming a categorical column in
#'   `seu[[]]`. The function checks that the column exists before fitting the
#'   model. When supplied, residuals are summarized and plotted by that column.
#' @param covariates Optional data frame or named list of numeric per-cell
#'   covariates. Each covariate is fitted alone and receives its own residual
#'   boxplot. Neighbourhood density is always included as a separate, independent
#'   model and is not added to the covariate models.
#' @param log_density Logical; if `TRUE`, regress on `log(mean_knn_dist)`.
#'   The returned `mean_knn_dist` values remain on the original distance scale.
#' @param cluster_of_interest Optional single group value to highlight in the
#'   boxplot. Ignored when `group` is `NULL`.
#' @param plot_file Optional base file path for 300-dpi PNG output when `group`
#'   is supplied. Default `NULL` does not save. With multiple covariates, the
#'   covariate name is appended to create one file per plot.
#'
#' @return A list with per-cell `mean_knn_dist`; a residual data frame with one
#'   column per separate model; the global Spearman `correlation`; named
#'   `regressions`; an optional combined `group_summary`; an optional Spearman
#'   `covariate_correlations` matrix; named `plots`; and optional named
#'   `plot_files`. Residuals are diagnostic values, not corrected scores.
#'
#' @examples
#' \dontrun{
#' embedding <- Seurat::Embeddings(seu[["pca"]])
#' result <- entropy_density_diagnostic(
#'   seu = seu, entropy = seu$entropy_mean_RNA, embedding = embedding,
#'   k = 15, group = "celltype",
#'   plot_file = NULL
#' )
#' result$correlation
#' result$group_summary
#' }
#'
#' @export
entropy_density_diagnostic <- function(
    seu,
    entropy,
    embedding,
    k = 30,
    group = NULL,
    covariates = NULL,
    log_density = FALSE,
    cluster_of_interest = NULL,
    plot_file = NULL) {

  if (!inherits(seu, "Seurat")) {
    stop("`seu` must be a Seurat object.", call. = FALSE)
  }

  if (!is.numeric(entropy) || !is.atomic(entropy) || !is.null(dim(entropy))) {
    stop("`entropy` must be a numeric vector.", call. = FALSE)
  }
  cell_names <- names(entropy)
  entropy <- as.numeric(entropy)
  if (length(entropy) < 2L || any(!is.finite(entropy))) {
    stop("`entropy` must contain at least two finite values.", call. = FALSE)
  }

  if (!is.matrix(embedding) || !is.numeric(embedding)) {
    stop("`embedding` must be a numeric matrix with cells in rows.", call. = FALSE)
  }
  if (nrow(embedding) != length(entropy)) {
    stop("`embedding` must have one row per value in `entropy`.", call. = FALSE)
  }
  if (ncol(seu) != length(entropy)) {
    stop("`seu`, `entropy`, and `embedding` must contain the same number of cells.",
         call. = FALSE)
  }
  seurat_cell_names <- colnames(seu)
  if (!is.null(cell_names) && !identical(cell_names, seurat_cell_names)) {
    stop("Names of `entropy` must exactly match the cell order in `seu`.", call. = FALSE)
  }
  if (!is.null(rownames(embedding)) &&
      !identical(rownames(embedding), seurat_cell_names)) {
    stop("Row names of `embedding` must exactly match the cell order in `seu`.",
         call. = FALSE)
  }
  cell_names <- seurat_cell_names
  if (ncol(embedding) < 1L || any(!is.finite(embedding))) {
    stop("`embedding` must have at least one column and only finite values.", call. = FALSE)
  }

  if (length(k) != 1L || !is.numeric(k) || !is.finite(k) ||
      k != as.integer(k) || k < 1L || k >= length(entropy)) {
    stop("`k` must be a single positive integer smaller than the number of cells.",
         call. = FALSE)
  }
  k <- as.integer(k)

  if (length(log_density) != 1L || !is.logical(log_density) || is.na(log_density)) {
    stop("`log_density` must be `TRUE` or `FALSE`.", call. = FALSE)
  }

  group_name <- group
  if (!is.null(group_name)) {
    if (!is.character(group_name) || length(group_name) != 1L ||
        is.na(group_name) || !nzchar(group_name)) {
      stop("`group` must be one non-empty metadata column name.", call. = FALSE)
    }
    metadata <- seu[[]]
    if (!group_name %in% names(metadata)) {
      stop(sprintf("Group column '%s' was not found in `seu` metadata.", group_name),
           call. = FALSE)
    }
    group_values <- metadata[[group_name]]
    if (!is.atomic(group_values) || is.list(group_values) || anyNA(group_values)) {
      stop("The selected group metadata column must contain non-missing categorical values.",
           call. = FALSE)
    }
    group <- if (is.factor(group_values)) {
      group_values
    } else {
      factor(as.character(group_values), levels = unique(as.character(group_values)))
    }
  }

  if (!is.null(cluster_of_interest)) {
    if (length(cluster_of_interest) != 1L || is.na(cluster_of_interest)) {
      stop("`cluster_of_interest` must be one non-missing group value.", call. = FALSE)
    }
    if (is.null(group)) {
      warning("`cluster_of_interest` is ignored because `group` is NULL.", call. = FALSE)
    } else if (!as.character(cluster_of_interest) %in% as.character(group)) {
      warning("`cluster_of_interest` is not present in `group`; no group will be highlighted.",
              call. = FALSE)
    }
  }

  covariate_data <- .entropy_diagnostic_covariates(covariates, length(entropy))

  knn <- RANN::nn2(data = embedding, query = embedding, k = k + 1L)
  mean_knn_dist <- vapply(seq_len(nrow(embedding)), function(i) {
    non_self <- knn$nn.dists[i, knn$nn.idx[i, ] != i]
    if (length(non_self) < k) {
      stop("Could not identify enough non-self neighbours for every cell.", call. = FALSE)
    }
    mean(non_self[seq_len(k)])
  }, numeric(1))
  names(mean_knn_dist) <- cell_names

  if (log_density && any(mean_knn_dist <= 0)) {
    stop(
      "Cannot use `log_density = TRUE` because at least one mean neighbour distance is zero.",
      call. = FALSE
    )
  }

  correlation_test <- suppressWarnings(stats::cor.test(
    entropy, mean_knn_dist, method = "spearman", exact = FALSE
  ))
  correlation <- c(
    rho = unname(correlation_test$estimate),
    p_value = correlation_test$p.value
  )

  density_term <- if (log_density) log(mean_knn_dist) else mean_knn_dist
  density_name <- if (log_density) "log_mean_knn_dist" else "mean_knn_dist"
  model_specs <- stats::setNames(list(density_term), density_name)
  if (!is.null(covariate_data)) {
    covariate_specs <- lapply(seq_len(ncol(covariate_data)), function(i) {
      covariate_data[[i]]
    })
    names(covariate_specs) <- names(covariate_data)
    model_specs <- c(model_specs, covariate_specs)
  }

  model_results <- lapply(seq_along(model_specs), function(i) {
    .entropy_diagnostic_fit(
      entropy = entropy,
      predictor = model_specs[[i]],
      predictor_name = names(model_specs)[i],
      cell_names = cell_names
    )
  })
  names(model_results) <- names(model_specs)

  residuals <- as.data.frame(
    do.call(cbind, lapply(model_results, `[[`, "residuals")),
    check.names = FALSE
  )
  names(residuals) <- names(model_results)
  rownames(residuals) <- cell_names
  regressions <- lapply(model_results, `[[`, "regression")

  group_summary <- NULL
  plots <- NULL
  saved_plot_files <- NULL
  if (!is.null(group)) {
    group_summary <- do.call(rbind, lapply(names(model_results), function(model_name) {
      split_residuals <- split(residuals[[model_name]], group, drop = TRUE)
      data.frame(
        model = model_name,
        group = names(split_residuals),
        n = vapply(split_residuals, length, integer(1)),
        mean_residual = vapply(split_residuals, mean, numeric(1)),
        median_residual = vapply(split_residuals, stats::median, numeric(1)),
        row.names = NULL,
        check.names = FALSE
      )
    }))
    message("Per-group residual summary (positive = entropy above density/covariate trend):")
    print(group_summary, row.names = FALSE)

    highlight_value <- as.character(cluster_of_interest)
    plots <- lapply(names(model_results), function(model_name) {
      plot_data <- data.frame(group = group, residual = residuals[[model_name]])
      plot_data$highlight <- if (!is.null(cluster_of_interest) &&
                                 highlight_value %in% as.character(group)) {
        ifelse(as.character(group) == highlight_value, "highlighted", "other")
      } else {
        "other"
      }
      model_subtitle <- paste("Independent model: entropy ~", model_name)
      model_title <- if (model_name %in% c("mean_knn_dist", "log_mean_knn_dist")) {
        "Entropy residual above/below density trend, per group"
      } else {
        "Entropy residual above/below covariate trend, per group"
      }

      ggplot2::ggplot(
        plot_data,
        ggplot2::aes(x = group, y = residual, fill = highlight)
      ) +
        ggplot2::geom_boxplot(outlier.alpha = 0.45) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
        ggplot2::scale_fill_manual(values = c(other = "grey75", highlighted = "#D55E00")) +
        ggplot2::labs(
          title = model_title,
          subtitle = model_subtitle,
          x = NULL,
          y = "residual (observed - predicted entropy)"
        ) +
        ggplot2::guides(fill = "none") +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    })
    names(plots) <- names(model_results)

    if (!is.null(plot_file)) {
      if (length(plot_file) != 1L || !is.character(plot_file) || is.na(plot_file) ||
          !nzchar(plot_file)) {
        stop("`plot_file` must be a non-empty file path or `NULL`.", call. = FALSE)
      }
      base_file <- if (grepl("\\.png$", plot_file, ignore.case = TRUE)) {
        plot_file
      } else {
        paste0(plot_file, ".png")
      }
      plot_paths <- if (length(plots) == 1L) {
        stats::setNames(base_file, names(plots))
      } else {
        suffixes <- make.unique(gsub("[^A-Za-z0-9_-]+", "_", names(plots)))
        stem <- sub("\\.png$", "", base_file, ignore.case = TRUE)
        stats::setNames(paste0(stem, "_", suffixes, ".png"), names(plots))
      }
      for (model_name in names(plots)) {
        ggplot2::ggsave(
          filename = plot_paths[[model_name]],
          plot = plots[[model_name]],
          width = 9,
          height = 6,
          units = "in",
          dpi = 300
        )
      }
      saved_plot_files <- vapply(plot_paths, normalizePath, character(1), mustWork = TRUE)
    }
  }

  covariate_correlations <- NULL
  if (!is.null(covariate_data)) {
    correlation_data <- data.frame(density_term, covariate_data, check.names = FALSE)
    names(correlation_data)[1L] <- density_name
    covariate_correlations <- stats::cor(correlation_data, method = "spearman")
  }

  list(
    mean_knn_dist = mean_knn_dist,
    residuals = residuals,
    correlation = correlation,
    regressions = regressions,
    group_summary = group_summary,
    covariate_correlations = covariate_correlations,
    plots = plots,
    plot_files = saved_plot_files
  )
}

.entropy_diagnostic_fit <- function(
    entropy,
    predictor,
    predictor_name,
    cell_names) {
  model_data <- data.frame(entropy = entropy, .predictor = predictor)
  fit <- stats::lm(entropy ~ .predictor, data = model_data)
  fit_summary <- summary(fit)

  coefficient_table <- as.data.frame(fit_summary$coefficients, check.names = FALSE)
  names(coefficient_table) <- c("estimate", "std_error", "t_value", "p_value")
  coefficient_table$term <- rownames(coefficient_table)
  coefficient_table$term[coefficient_table$term == ".predictor"] <- predictor_name
  rownames(coefficient_table) <- NULL
  coefficient_table <- coefficient_table[
    , c("term", "estimate", "std_error", "t_value", "p_value")
  ]

  f_statistic <- fit_summary$fstatistic
  model_p_value <- if (is.null(f_statistic)) {
    NA_real_
  } else {
    stats::pf(f_statistic[[1L]], f_statistic[[2L]], f_statistic[[3L]], lower.tail = FALSE)
  }

  residuals <- unname(stats::residuals(fit))
  names(residuals) <- cell_names
  list(
    residuals = residuals,
    regression = list(
      coefficients = coefficient_table,
      r_squared = unname(fit_summary$r.squared),
      adjusted_r_squared = unname(fit_summary$adj.r.squared),
      model_p_value = unname(model_p_value)
    )
  )
}

.entropy_diagnostic_covariates <- function(covariates, n_cells) {
  if (is.null(covariates)) {
    return(NULL)
  }
  if (is.list(covariates) && !is.data.frame(covariates)) {
    if (is.null(names(covariates)) || any(!nzchar(names(covariates)))) {
      stop("A covariate list must be named.", call. = FALSE)
    }
    covariates <- tryCatch(
      as.data.frame(covariates, check.names = FALSE),
      error = function(e) stop("All covariates must have the same per-cell length.", call. = FALSE)
    )
  }
  if (!is.data.frame(covariates) || ncol(covariates) < 1L) {
    stop("`covariates` must be a data frame or named list with at least one column.",
         call. = FALSE)
  }
  if (nrow(covariates) != n_cells) {
    stop("Every covariate must have one value per cell.", call. = FALSE)
  }
  if (is.null(names(covariates)) || any(!nzchar(names(covariates))) ||
      anyDuplicated(names(covariates))) {
    stop("Covariates must have unique, non-empty names.", call. = FALSE)
  }
  numeric_columns <- vapply(covariates, is.numeric, logical(1))
  if (!all(numeric_columns)) {
    stop("All covariates must be numeric.", call. = FALSE)
  }
  if (any(!vapply(covariates, function(x) all(is.finite(x)), logical(1)))) {
    stop("All covariates must contain only finite values.", call. = FALSE)
  }
  if (any(names(covariates) %in% c("entropy", "mean_knn_dist", "log_mean_knn_dist"))) {
    stop("Covariate names must not be `entropy`, `mean_knn_dist`, or `log_mean_knn_dist`.",
         call. = FALSE)
  }
  covariates
}
