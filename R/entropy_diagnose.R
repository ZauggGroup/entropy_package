#' Diagnose entropy against density and other covariates
#'
#' `entropy_diagnose()` asks whether entropy is elevated in each cluster and
#' whether that elevation remains when cells are compared at matched values of
#' neighbourhood density or another covariate. Neighbourhood density is always
#' included as `mean_knn_dist`. Larger values mean sparser neighbourhoods.
#'
#' The matched comparison reports effect sizes and support, not p values. With
#' tens of thousands of cells, very small differences can be statistically
#' significant, so a p value would make the flag misleading.
#'
#' The four flags mean:
#'
#' * `not elevated`: this cluster does not have high entropy, so the question of what explains its entropy does not apply.
#' * `cannot tell`: this cluster occupies a range of the covariate where too few other cells exist, so there is no fair comparison group. No conclusion can be drawn.
#' * `explained`: at matched covariate values, this cluster is not still higher than other cells. The covariate accounts for its elevation, and the direction may be reversed.
#' * `not explained`: at matched covariate values, this cluster is specifically still higher than other cells. The covariate does not account for it.
#'
#' If a joint model is requested, only joint and incremental R squared values
#' are returned. Joint-model coefficients are deliberately omitted because
#' correlated predictors make them unstable and easy to misread as independent
#' contributions.
#'
#' @param seu A Seurat object whose cells correspond, in the same order, to
#'   `entropy` and the rows of `embedding`.
#' @param entropy Numeric vector with one finite value per cell.
#' @param embedding Numeric matrix with cells in rows, in the same order as
#'   `entropy`.
#' @param k Positive integer number of nearest neighbours. It must be smaller
#'   than the number of cells.
#' @param group One non-empty character string naming a categorical column in
#'   `seu[[]]`. Missing group values are not allowed.
#' @param covariates Optional data frame or named list of finite numeric
#'   covariates, with one value per cell.
#' @param covariate_type Optional named character vector. Values must be
#'   `"technical"`, `"biological"`, or `"unknown"`. Unspecified user
#'   covariates default to `"unknown"`.
#' @param covariate_intervention Optional named character vector. Values must
#'   be `"downsample_counts"`, `"subsample_cells"`, `"remove_genes"`, or
#'   `"none"`. Unspecified user covariates default to `"none"`.
#' @param n_bins Positive integer number of quantile bins.
#' @param min_bin_n Positive integer minimum number of cells required inside
#'   and outside a cluster for a bin to be usable.
#' @param max_bin_imbalance Non-negative maximum allowed difference between the
#'   within-bin covariate medians for cluster and non-cluster cells, divided by
#'   the covariate IQR over all cells. Quantile bins alone do not guarantee
#'   comparable covariate values inside a bin, especially in a long tail.
#' @param min_cluster_n Positive integer minimum cluster size for the within
#'   cluster fit.
#' @param support_threshold Number from zero to one. Clusters below this
#'   matched-support fraction receive `"cannot tell"`.
#' @param effect_threshold Non-negative standardized matched difference below
#'   which an elevated cluster is called `"explained"`.
#' @param elevated_threshold Standardized entropy difference above which a
#'   cluster is considered elevated.
#' @param retained_min_delta Non-negative minimum absolute standardized entropy
#'   difference required to report `retained`. The ratio is unstable for
#'   marginal elevations or deficits. A missing value means the cluster is not
#'   far enough from the global mean for the fraction to be interpretable.
#' @param joint_model Logical. If `TRUE`, fit a joint model when more than one
#'   covariate is present.
#' @param provenance Optional named list describing how entropy was computed.
#' @param verbose Logical. If `TRUE`, show short progress messages.
#'
#' @return An object of class `entropy_diagnostic`. Its main result is
#'   `cluster_table`, with one row per cluster and covariate. Its `direction`
#'   column states whether entropy is elevated or depleted relative to other
#'   cells. `retained` is the matched standardized difference divided by the
#'   cluster's standardized entropy elevation or deficit. It is reported only
#'   when the absolute standardized difference is beyond `elevated_threshold`
#'   and at least `retained_min_delta`. A missing value means the cluster is not
#'   elevated or depleted enough for the fraction to be interpretable.
#'   `bin_table` contains every matched-bin comparison. No p values are computed
#'   for the matched comparison.
#' @export
entropy_diagnose <- function(
    seu,
    entropy,
    embedding,
    k = 30,
    group,
    covariates = NULL,
    covariate_type = NULL,
    covariate_intervention = NULL,
    n_bins = 10,
    min_bin_n = 20,
    max_bin_imbalance = 0.25,
    min_cluster_n = 30,
    support_threshold = 0.5,
    effect_threshold = 0.1,
    elevated_threshold = 0.1,
    retained_min_delta = 0.5,
    joint_model = TRUE,
    provenance = NULL,
    verbose = TRUE) {

  validated <- .entropy_validate_inputs(
    seu = seu,
    entropy = entropy,
    embedding = embedding,
    k = k,
    group = group,
    covariates = covariates
  )
  entropy <- validated$entropy
  embedding <- validated$embedding
  group_values <- validated$group
  group_name <- validated$group_name
  cell_names <- validated$cell_names
  covariate_data <- validated$covariates
  k <- validated$k

  n_bins <- .entropy_positive_integer(n_bins, "n_bins")
  min_bin_n <- .entropy_positive_integer(min_bin_n, "min_bin_n")
  max_bin_imbalance <- .entropy_bounded_number(
    max_bin_imbalance, "max_bin_imbalance", lower = 0
  )
  min_cluster_n <- .entropy_positive_integer(min_cluster_n, "min_cluster_n")
  support_threshold <- .entropy_bounded_number(
    support_threshold, "support_threshold", lower = 0, upper = 1
  )
  effect_threshold <- .entropy_bounded_number(
    effect_threshold, "effect_threshold", lower = 0
  )
  elevated_threshold <- .entropy_bounded_number(
    elevated_threshold, "elevated_threshold", lower = -Inf
  )
  retained_min_delta <- .entropy_bounded_number(
    retained_min_delta, "retained_min_delta", lower = 0
  )
  .entropy_scalar_logical(joint_model, "joint_model")
  .entropy_scalar_logical(verbose, "verbose")

  if (!is.null(provenance) &&
      (!is.list(provenance) || is.null(names(provenance)) ||
       any(!nzchar(names(provenance))) || anyDuplicated(names(provenance)))) {
    stop("`provenance` must be a named list or `NULL`.", call. = FALSE)
  }
  if (!is.null(provenance$k) &&
      (length(provenance$k) != 1L || is.na(provenance$k) || provenance$k != k)) {
    warning(
      sprintf(
        "Entropy provenance records k = %s, but this diagnostic uses k = %s.",
        paste(provenance$k, collapse = ", "), k
      ),
      call. = FALSE
    )
  }
  embedding_name <- .entropy_embedding_name(embedding)
  if (!is.null(provenance$reduction) && !is.null(embedding_name) &&
      (length(provenance$reduction) != 1L || is.na(provenance$reduction) ||
       provenance$reduction != embedding_name)) {
    warning(
      sprintf(
        "Entropy provenance records reduction '%s', but this embedding is named '%s'.",
        paste(provenance$reduction, collapse = ", "), embedding_name
      ),
      call. = FALSE
    )
  }

  user_covariate_names <- if (is.null(covariate_data)) character() else names(covariate_data)
  type_values <- .entropy_named_choices(
    covariate_type,
    user_covariate_names,
    allowed = c("technical", "biological", "unknown"),
    default = "unknown",
    argument = "covariate_type"
  )
  intervention_values <- .entropy_named_choices(
    covariate_intervention,
    user_covariate_names,
    allowed = c("downsample_counts", "subsample_cells", "remove_genes", "none"),
    default = "none",
    argument = "covariate_intervention"
  )

  if (verbose) {
    message("Computing mean nearest-neighbour distance.")
  }
  knn <- RANN::nn2(data = embedding, query = embedding, k = k + 1L)
  mean_knn_dist <- vapply(seq_len(nrow(embedding)), function(i) {
    non_self <- knn$nn.dists[i, knn$nn.idx[i, ] != i]
    if (length(non_self) < k) {
      stop("Could not identify enough non-self neighbours for every cell.", call. = FALSE)
    }
    mean(non_self[seq_len(k)])
  }, numeric(1))

  all_covariates <- data.frame(mean_knn_dist = mean_knn_dist, check.names = FALSE)
  if (!is.null(covariate_data)) {
    all_covariates <- cbind(all_covariates, covariate_data)
  }
  all_types <- c(mean_knn_dist = "unknown", type_values)
  all_interventions <- c(mean_knn_dist = "subsample_cells", intervention_values)

  entropy_sd <- stats::sd(entropy)
  if (!is.finite(entropy_sd) || entropy_sd == 0) {
    entropy_sd <- NA_real_
  }

  if (verbose) {
    message("Computing global covariate summaries and loess curves.")
  }
  q1 <- lapply(names(all_covariates), function(covariate_name) {
    predictor <- all_covariates[[covariate_name]]
    fit <- stats::lm(entropy ~ predictor)
    correlation <- suppressWarnings(stats::cor.test(
      entropy, predictor, method = "spearman", exact = FALSE
    ))
    loess_result <- .entropy_loess_curve(entropy, predictor)
    list(
      info = data.frame(
        covariate = covariate_name,
        type = unname(all_types[[covariate_name]]),
        intervention = unname(all_interventions[[covariate_name]]),
        r_squared = .entropy_model_r_squared(fit),
        spearman_rho = unname(correlation$estimate),
        spearman_p = unname(correlation$p.value),
        loess_subsampled = loess_result$subsampled,
        stringsAsFactors = FALSE,
        check.names = FALSE
      ),
      curve = loess_result$curve
    )
  })
  covariate_info <- do.call(rbind, lapply(q1, `[[`, "info"))
  rownames(covariate_info) <- NULL
  if (max_bin_imbalance > 0.1) {
    strong <- is.finite(covariate_info$spearman_rho) &
      abs(covariate_info$spearman_rho) > 0.7
    for (covariate_name in covariate_info$covariate[strong]) {
      rho <- covariate_info$spearman_rho[
        covariate_info$covariate == covariate_name
      ][[1L]]
      warning(sprintf(
        paste(
          "Covariate '%s' has a strong entropy association (Spearman rho = %.3f)",
          "and max_bin_imbalance = %.2f. Tighten max_bin_imbalance to 0.1 or",
          "below and increase n_bins to reduce residual mismatch."
        ),
        covariate_name, rho, max_bin_imbalance
      ), call. = FALSE)
    }
  }
  loess_curves <- lapply(q1, `[[`, "curve")
  names(loess_curves) <- names(all_covariates)

  if (verbose) {
    message("Computing cluster comparisons.")
  }
  clusters <- levels(group_values)
  cluster_rows <- vector("list", length(clusters) * ncol(all_covariates))
  bin_rows <- vector("list", length(clusters) * ncol(all_covariates))
  row_index <- 0L

  for (covariate_name in names(all_covariates)) {
    predictor <- all_covariates[[covariate_name]]
    predictor_sd <- stats::sd(predictor)
    if (!is.finite(predictor_sd) || predictor_sd == 0) {
      predictor_sd <- NA_real_
    }
    predictor_iqr <- stats::IQR(predictor)
    if (!is.finite(predictor_iqr) || predictor_iqr == 0) {
      predictor_iqr <- NA_real_
    }
    bins <- .entropy_quantile_bins(predictor, n_bins, covariate_name)

    for (cluster_name in clusters) {
      row_index <- row_index + 1L
      in_cluster <- group_values == cluster_name
      n_cluster <- sum(in_cluster)
      entropy_delta_std <- .entropy_standardized_difference(
        entropy[in_cluster], entropy[!in_cluster], entropy_sd
      )
      covariate_delta_std <- .entropy_standardized_difference(
        predictor[in_cluster], predictor[!in_cluster], predictor_sd
      )
      entropy_elevated <- is.finite(entropy_delta_std) &&
        entropy_delta_std > elevated_threshold

      this_bins <- lapply(seq_len(bins$n_bins), function(bin_index) {
        in_bin <- bins$index == bin_index
        n_in <- sum(in_cluster & in_bin)
        n_out <- sum(!in_cluster & in_bin)
        median_cov_in <- if (n_in > 0L) {
          stats::median(predictor[in_cluster & in_bin])
        } else {
          NA_real_
        }
        median_cov_out <- if (n_out > 0L) {
          stats::median(predictor[!in_cluster & in_bin])
        } else {
          NA_real_
        }
        imbalance <- if (is.finite(predictor_iqr) &&
                         is.finite(median_cov_in) &&
                         is.finite(median_cov_out)) {
          abs(median_cov_in - median_cov_out) / predictor_iqr
        } else {
          NA_real_
        }
        usable <- n_in >= min_bin_n && n_out >= min_bin_n &&
          is.finite(imbalance) && imbalance <= max_bin_imbalance
        difference <- if (usable) {
          stats::median(entropy[in_cluster & in_bin]) -
            stats::median(entropy[!in_cluster & in_bin])
        } else {
          NA_real_
        }
        data.frame(
          cluster = cluster_name,
          covariate = covariate_name,
          bin = bin_index,
          bin_lower = bins$lower[bin_index],
          bin_upper = bins$upper[bin_index],
          n_in = n_in,
          n_out = n_out,
          median_cov_in = median_cov_in,
          median_cov_out = median_cov_out,
          imbalance = imbalance,
          usable = usable,
          diff = difference,
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
      })
      this_bins <- do.call(rbind, this_bins)
      usable <- this_bins$usable
      support <- sum(this_bins$n_in[usable]) / n_cluster
      matched_diff <- if (any(usable)) {
        stats::weighted.mean(this_bins$diff[usable], this_bins$n_in[usable])
      } else {
        NA_real_
      }
      matched_diff_std <- matched_diff / entropy_sd
      direction <- if (is.finite(entropy_delta_std) && entropy_delta_std < 0) {
        "depleted"
      } else {
        "elevated"
      }
      retained <- if (is.finite(entropy_delta_std) &&
                      entropy_delta_std != 0 &&
                      abs(entropy_delta_std) > elevated_threshold &&
                      abs(entropy_delta_std) >= retained_min_delta) {
        matched_diff_std / entropy_delta_std
      } else {
        NA_real_
      }

      within <- .entropy_within_fit(
        entropy = entropy[in_cluster],
        predictor = predictor[in_cluster],
        minimum = min_cluster_n
      )
      flag <- if (!entropy_elevated) {
        "not elevated"
      } else if (support < support_threshold) {
        "cannot tell"
      } else if (is.finite(matched_diff_std) &&
                 matched_diff_std > effect_threshold) {
        "not explained"
      } else {
        "explained"
      }
      reversed <- if (entropy_elevated) {
        is.finite(matched_diff_std) && matched_diff_std < -effect_threshold
      } else {
        NA
      }

      cluster_rows[[row_index]] <- data.frame(
        cluster = cluster_name,
        covariate = covariate_name,
        n_cells = n_cluster,
        mean_entropy = mean(entropy[in_cluster]),
        median_entropy = stats::median(entropy[in_cluster]),
        entropy_delta_std = entropy_delta_std,
        entropy_elevated = entropy_elevated,
        direction = direction,
        mean_covariate = mean(predictor[in_cluster]),
        median_covariate = stats::median(predictor[in_cluster]),
        covariate_delta_std = covariate_delta_std,
        matched_diff = matched_diff,
        matched_diff_std = matched_diff_std,
        retained = retained,
        support = support,
        n_usable_bins = sum(usable),
        within_slope = within$slope,
        within_slope_std = within$slope_std,
        within_r_squared = within$r_squared,
        within_reason = within$reason,
        reversed = reversed,
        flag = flag,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      bin_rows[[row_index]] <- this_bins
    }
  }
  cluster_table <- do.call(rbind, cluster_rows)
  rownames(cluster_table) <- NULL
  bin_table <- do.call(rbind, bin_rows)
  rownames(bin_table) <- NULL

  joint <- NULL
  if (joint_model && ncol(all_covariates) > 1L) {
    if (verbose) {
      message("Computing the joint covariate model.")
    }
    joint <- .entropy_joint_model(entropy, all_covariates)
  }

  cell_data <- data.frame(
    cell = cell_names,
    cluster = as.character(group_values),
    entropy = entropy,
    all_covariates,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  result <- list(
    provenance = list(
      k = k,
      n_cells = length(entropy),
      n_covariates = ncol(all_covariates),
      group_column = group_name,
      embedding_dims = ncol(embedding),
      user_provenance = provenance,
      call = match.call()
    ),
    cell_data = cell_data,
    covariate_info = covariate_info,
    loess_curves = loess_curves,
    cluster_table = cluster_table,
    bin_table = bin_table,
    joint = joint,
    thresholds = list(
      n_bins = n_bins,
      min_bin_n = min_bin_n,
      max_bin_imbalance = max_bin_imbalance,
      min_cluster_n = min_cluster_n,
      support_threshold = support_threshold,
      effect_threshold = effect_threshold,
      elevated_threshold = elevated_threshold,
      retained_min_delta = retained_min_delta
    )
  )
  class(result) <- "entropy_diagnostic"
  result
}


#' Print an entropy diagnostic
#'
#' @param x An `entropy_diagnostic` object.
#' @param ... Additional arguments, currently ignored.
#'
#' @return `x`, invisibly.
#' @keywords internal
#' @export
print.entropy_diagnostic <- function(x, ...) {
  if (!inherits(x, "entropy_diagnostic")) {
    stop("`x` must be an entropy_diagnostic object.", call. = FALSE)
  }
  n_clusters <- length(unique(x$cluster_table$cluster))
  cat(sprintf(
    "Entropy diagnostic: %d cells, %d clusters, k = %d, %d covariates\n\n",
    x$provenance$n_cells,
    n_clusters,
    x$provenance$k,
    x$provenance$n_covariates
  ))

  q1 <- x$covariate_info[
    order(x$covariate_info$r_squared, decreasing = TRUE, na.last = TRUE),
    c("covariate", "type", "r_squared", "spearman_rho", "spearman_p"),
    drop = FALSE
  ]
  print(q1, row.names = FALSE, digits = 3)

  if (!is.null(x$joint)) {
    cat(sprintf("\nJoint R squared: %.3f\n", x$joint$joint_r_squared))
  }
  flag_counts <- table(factor(
    x$cluster_table$flag,
    levels = c("not elevated", "cannot tell", "explained", "not explained")
  ))
  cat("\nFlags:\n")
  print(flag_counts)
  cat("\nUse entropy_diagnose_report() to write the PDF report.\n")
  invisible(x)
}


.entropy_validate_inputs <- function(seu, entropy, embedding, k, group, covariates) {
  if (!inherits(seu, "Seurat")) {
    stop("`seu` must be a Seurat object.", call. = FALSE)
  }
  if (!is.numeric(entropy) || !is.atomic(entropy) || !is.null(dim(entropy))) {
    stop("`entropy` must be a numeric vector.", call. = FALSE)
  }
  entropy_names <- names(entropy)
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
  cell_names <- colnames(seu)
  if (!is.null(entropy_names) && !identical(entropy_names, cell_names)) {
    stop("Names of `entropy` must exactly match the cell order in `seu`.", call. = FALSE)
  }
  if (!is.null(rownames(embedding)) && !identical(rownames(embedding), cell_names)) {
    stop("Row names of `embedding` must exactly match the cell order in `seu`.",
         call. = FALSE)
  }
  if (ncol(embedding) < 1L || any(!is.finite(embedding))) {
    stop("`embedding` must have at least one column and only finite values.", call. = FALSE)
  }
  if (length(k) != 1L || !is.numeric(k) || !is.finite(k) ||
      k != as.integer(k) || k < 1L || k >= length(entropy)) {
    stop("`k` must be a single positive integer smaller than the number of cells.",
         call. = FALSE)
  }

  if (missing(group) || !is.character(group) || length(group) != 1L ||
      is.na(group) || !nzchar(group)) {
    stop("`group` must be one non-empty metadata column name.", call. = FALSE)
  }
  metadata <- seu[[]]
  if (!group %in% names(metadata)) {
    stop(sprintf("Group column '%s' was not found in `seu` metadata.", group),
         call. = FALSE)
  }
  group_values <- metadata[[group]]
  if ((!is.factor(group_values) && !is.character(group_values)) || anyNA(group_values)) {
    stop("The selected group metadata column must contain non-missing categorical values.",
         call. = FALSE)
  }
  group_values <- factor(
    as.character(group_values), levels = unique(as.character(group_values))
  )
  if (nlevels(group_values) < 2L) {
    stop("The selected group metadata column must contain at least two groups.",
         call. = FALSE)
  }

  covariate_data <- .entropy_diagnostic_covariates(covariates, length(entropy))
  list(
    entropy = entropy,
    embedding = embedding,
    k = as.integer(k),
    group = group_values,
    group_name = group,
    cell_names = cell_names,
    covariates = covariate_data
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
      error = function(e) {
        stop("All covariates must have the same per-cell length.", call. = FALSE)
      }
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
  if (any(!vapply(covariates, function(value) all(is.finite(value)), logical(1)))) {
    stop("All covariates must contain only finite values.", call. = FALSE)
  }
  if (any(names(covariates) %in% c("cell", "cluster", "entropy", "mean_knn_dist"))) {
    stop("Covariate names must not be `cell`, `cluster`, `entropy`, or `mean_knn_dist`.",
         call. = FALSE)
  }
  covariates
}


.entropy_positive_integer <- function(value, argument) {
  if (length(value) != 1L || !is.numeric(value) || !is.finite(value) ||
      value != as.integer(value) || value < 1L) {
    stop(sprintf("`%s` must be a single positive integer.", argument), call. = FALSE)
  }
  as.integer(value)
}


.entropy_bounded_number <- function(value, argument, lower, upper = Inf) {
  if (length(value) != 1L || !is.numeric(value) || !is.finite(value) ||
      value < lower || value > upper) {
    stop(
      sprintf("`%s` must be a single finite number between %s and %s.",
              argument, lower, upper),
      call. = FALSE
    )
  }
  as.numeric(value)
}


.entropy_scalar_logical <- function(value, argument) {
  if (length(value) != 1L || !is.logical(value) || is.na(value)) {
    stop(sprintf("`%s` must be `TRUE` or `FALSE`.", argument), call. = FALSE)
  }
  invisible(value)
}


.entropy_named_choices <- function(value, covariates, allowed, default, argument) {
  result <- stats::setNames(rep(default, length(covariates)), covariates)
  if (is.null(value)) {
    return(result)
  }
  if (!is.character(value) || is.null(names(value)) ||
      any(!nzchar(names(value))) || anyDuplicated(names(value))) {
    stop(sprintf("`%s` must be a named character vector.", argument), call. = FALSE)
  }
  unknown <- setdiff(names(value), covariates)
  if (length(unknown)) {
    stop(sprintf(
      "`%s` contains unknown covariate names: %s.",
      argument, paste(unknown, collapse = ", ")
    ), call. = FALSE)
  }
  if (anyNA(value) || any(!value %in% allowed)) {
    stop(sprintf(
      "Every `%s` value must be one of: %s.",
      argument, paste(allowed, collapse = ", ")
    ), call. = FALSE)
  }
  result[names(value)] <- value
  result
}


.entropy_embedding_name <- function(embedding) {
  possible <- c(attr(embedding, "reduction"), attr(embedding, "reduction_name"))
  possible <- possible[!is.na(possible) & nzchar(possible)]
  if (length(possible)) as.character(possible[[1L]]) else NULL
}


.entropy_model_r_squared <- function(fit) {
  value <- unname(suppressWarnings(summary(fit)$r.squared))
  if (length(value) != 1L || !is.finite(value)) NA_real_ else value
}


.entropy_loess_curve <- function(entropy, predictor) {
  n <- length(entropy)
  keep <- if (n > 10000L) {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    }
    on.exit({
      if (had_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(1L)
    sort(sample.int(n, 10000L))
  } else {
    seq_len(n)
  }
  x <- predictor[keep]
  y <- entropy[keep]
  grid <- seq(min(x), max(x), length.out = 200L)
  fitted <- rep(NA_real_, length(grid))
  if (length(unique(x)) >= 3L) {
    fit <- tryCatch(
      suppressWarnings(stats::loess(y ~ x, span = 0.75)),
      error = function(e) NULL
    )
    if (!is.null(fit)) {
      fitted <- tryCatch(
        as.numeric(stats::predict(fit, newdata = data.frame(x = grid))),
        error = function(e) rep(NA_real_, length(grid))
      )
    }
  }
  list(
    curve = data.frame(x = grid, y = fitted),
    subsampled = n > 10000L
  )
}


.entropy_quantile_bins <- function(predictor, n_bins, covariate_name) {
  requested_breaks <- stats::quantile(
    predictor,
    probs = seq(0, 1, length.out = n_bins + 1L),
    names = FALSE
  )
  breaks <- unique(requested_breaks)
  actual_n <- max(1L, length(breaks) - 1L)
  if (actual_n < n_bins) {
    warning(sprintf(
      "Covariate '%s' produced %d distinct quantile bins instead of %d because of ties.",
      covariate_name, actual_n, n_bins
    ), call. = FALSE)
  }
  if (length(breaks) == 1L) {
    index <- rep(1L, length(predictor))
    lower <- upper <- breaks
  } else {
    index <- as.integer(cut(
      predictor, breaks = breaks, labels = FALSE, include.lowest = TRUE
    ))
    lower <- breaks[-length(breaks)]
    upper <- breaks[-1L]
  }
  list(index = index, lower = lower, upper = upper, n_bins = actual_n)
}


.entropy_standardized_difference <- function(inside, outside, scale_value) {
  if (!is.finite(scale_value) || !length(outside)) {
    return(NA_real_)
  }
  (mean(inside) - mean(outside)) / scale_value
}


.entropy_within_fit <- function(entropy, predictor, minimum) {
  empty <- list(slope = NA_real_, slope_std = NA_real_, r_squared = NA_real_)
  if (length(entropy) < minimum) {
    return(c(empty, list(reason = sprintf("fewer than %d cells", minimum))))
  }
  if (stats::sd(entropy) == 0 || stats::sd(predictor) == 0) {
    return(c(empty, list(reason = "entropy or covariate has no within-cluster variation")))
  }
  fit <- stats::lm(entropy ~ predictor)
  standardized_fit <- stats::lm(scale(entropy) ~ scale(predictor))
  list(
    slope = unname(stats::coef(fit)[[2L]]),
    slope_std = unname(stats::coef(standardized_fit)[[2L]]),
    r_squared = .entropy_model_r_squared(fit),
    reason = NA_character_
  )
}


.entropy_joint_model <- function(entropy, covariates) {
  model_data <- data.frame(.entropy = entropy, covariates, check.names = FALSE)
  full_fit <- stats::lm(.entropy ~ ., data = model_data)
  joint_r_squared <- .entropy_model_r_squared(full_fit)
  incremental <- vapply(names(covariates), function(covariate_name) {
    reduced <- model_data[, c(".entropy", setdiff(names(covariates), covariate_name)),
                          drop = FALSE]
    reduced_fit <- stats::lm(.entropy ~ ., data = reduced)
    joint_r_squared - .entropy_model_r_squared(reduced_fit)
  }, numeric(1))
  correlations <- stats::cor(covariates, method = "spearman")
  vif <- vapply(names(covariates), function(covariate_name) {
    response <- covariates[[covariate_name]]
    others <- covariates[, setdiff(names(covariates), covariate_name), drop = FALSE]
    auxiliary <- stats::lm(response ~ ., data = others)
    r_squared <- .entropy_model_r_squared(auxiliary)
    if (!is.finite(r_squared) || r_squared >= 1) Inf else 1 / (1 - r_squared)
  }, numeric(1))
  list(
    joint_r_squared = joint_r_squared,
    incremental_r_squared = incremental,
    covariate_correlations = correlations,
    vif = vif,
    note = paste(
      "Individual R squared values can sum to more than the joint R squared",
      "when covariates are correlated. This is expected, not a bug."
    )
  )
}
