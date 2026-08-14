#' Write PDF reports for an entropy diagnostic
#'
#' `entropy_diagnose_report()` renders an existing `entropy_diagnostic` object.
#' It does not refit models or change any diagnostic result. It writes a
#' paginated summary and one detailed PDF per covariate.
#'
#' The four flags mean:
#'
#' * `not elevated`: this cluster does not have high entropy, so the question of what explains its entropy does not apply.
#' * `cannot tell`: this cluster occupies a range of the covariate where too few other cells exist, so there is no fair comparison group. No conclusion can be drawn.
#' * `explained`: at matched covariate values, this cluster is not still higher than other cells. The covariate accounts for its elevation, and the direction may be reversed.
#' * `not explained`: at matched covariate values, this cluster is specifically still higher than other cells. The covariate does not account for it.
#'
#' @param x An `entropy_diagnostic` object returned by
#'   [entropy_diagnose()].
#' @param outdir Output directory. It is created recursively when needed.
#' @param prefix Optional string prepended to every output file name.
#' @param width PDF width in inches.
#' @param height PDF height in inches.
#' @param overwrite Logical. If `FALSE`, existing target files cause an error.
#' @param write_tables Logical. If `TRUE`, also write the cluster, bin, and
#'   covariate tables as CSV files.
#'
#' @return A named character vector of normalized output paths, invisibly.
#' @export
entropy_diagnose_report <- function(
    x,
    outdir = "entropy_diagnostics",
    prefix = NULL,
    width = 9,
    height = 7,
    overwrite = FALSE,
    write_tables = FALSE) {

  if (!inherits(x, "entropy_diagnostic")) {
    stop("`x` must be an entropy_diagnostic object.", call. = FALSE)
  }
  if (length(outdir) != 1L || !is.character(outdir) || is.na(outdir) ||
      !nzchar(outdir)) {
    stop("`outdir` must be one non-empty path.", call. = FALSE)
  }
  if (!is.null(prefix) &&
      (length(prefix) != 1L || !is.character(prefix) || is.na(prefix))) {
    stop("`prefix` must be one character string or `NULL`.", call. = FALSE)
  }
  width <- .entropy_report_size(width, "width")
  height <- .entropy_report_size(height, "height")
  .entropy_scalar_logical(overwrite, "overwrite")
  .entropy_scalar_logical(write_tables, "write_tables")

  if (!dir.exists(outdir)) {
    if (!dir.create(outdir, recursive = TRUE, showWarnings = FALSE)) {
      stop(sprintf("Could not create output directory '%s'.", outdir), call. = FALSE)
    }
  }
  prefix <- if (is.null(prefix)) "" else prefix
  safe_names <- make.unique(
    gsub("[^A-Za-z0-9_-]+", "_", x$covariate_info$covariate),
    sep = "_"
  )
  file_names <- c(
    summary = paste0(prefix, "summary.pdf"),
    stats::setNames(
      paste0(prefix, safe_names, ".pdf"),
      x$covariate_info$covariate
    )
  )
  if (write_tables) {
    file_names <- c(
      file_names,
      cluster_table = paste0(prefix, "cluster_table.csv"),
      bin_table = paste0(prefix, "bin_table.csv"),
      covariate_info = paste0(prefix, "covariate_info.csv")
    )
  }
  paths <- file.path(outdir, unname(file_names))
  names(paths) <- names(file_names)
  existing <- paths[file.exists(paths)]
  if (!overwrite && length(existing)) {
    stop(sprintf(
      "Refusing to overwrite existing report file: %s.",
      existing[[1L]]
    ), call. = FALSE)
  }

  summary_plots <- .entropy_summary_plots(x)
  summary_height <- max(height, 7.8 + 0.42 * min(12L, length(.entropy_cluster_order(x))))
  .entropy_write_pdf(paths[["summary"]], width, summary_height, summary_plots)

  for (covariate_name in x$covariate_info$covariate) {
    pages <- .entropy_covariate_pages(x, covariate_name)
    .entropy_write_pdf(paths[[covariate_name]], width, max(height, 8), pages)
  }

  if (write_tables) {
    utils::write.csv(x$cluster_table, paths[["cluster_table"]], row.names = FALSE)
    utils::write.csv(x$bin_table, paths[["bin_table"]], row.names = FALSE)
    utils::write.csv(x$covariate_info, paths[["covariate_info"]], row.names = FALSE)
  }

  paths <- vapply(paths, normalizePath, character(1), mustWork = TRUE)
  invisible(paths)
}


.entropy_cluster_order <- function(x) {
  entropy_rows <- !duplicated(x$cluster_table$cluster)
  values <- x$cluster_table[entropy_rows, c("cluster", "median_entropy"), drop = FALSE]
  values$cluster[order(values$median_entropy, decreasing = TRUE, na.last = TRUE)]
}


.entropy_report_size <- function(value, argument) {
  if (length(value) != 1L || !is.numeric(value) || !is.finite(value) || value <= 0) {
    stop(sprintf("`%s` must be one positive finite number.", argument), call. = FALSE)
  }
  as.numeric(value)
}


.entropy_write_pdf <- function(path, width, height, plots) {
  local({
    grDevices::pdf(path, width = width, height = height, onefile = TRUE)
    on.exit(grDevices::dev.off(), add = TRUE)
    for (plot in plots) {
      print(plot)
    }
    grDevices::dev.off()
    on.exit(NULL)
  })
  invisible(path)
}


.entropy_summary_plots <- function(x) {
  base_theme <- ggplot2::theme_minimal(base_size = 9)
  n_clusters <- length(unique(x$cluster_table$cluster))
  ordering <- .entropy_cluster_order(x)
  user_provenance <- x$provenance$user_provenance
  provenance_text <- if (is.null(user_provenance) || !length(user_provenance)) {
    "Entropy provenance: not supplied"
  } else {
    values <- vapply(names(user_provenance), function(name) {
      paste0(name, " = ", paste(user_provenance[[name]], collapse = ", "))
    }, character(1))
    paste("Entropy provenance:", paste(values, collapse = "; "))
  }
  header_text <- paste0(
    x$provenance$n_cells, " cells, ", n_clusters, " clusters, k = ",
    x$provenance$k, ", group = ", x$provenance$group_column, "\n",
    provenance_text
  )
  header <- ggplot2::ggplot() +
    ggplot2::annotate(
      "text", x = 0, y = 1, label = "Entropy diagnostic summary",
      hjust = 0, vjust = 1, fontface = "bold", size = 5
    ) +
    ggplot2::annotate(
      "text", x = 0, y = 0.58, label = header_text,
      hjust = 0, vjust = 1, size = 3.1
    ) +
    ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) + ggplot2::theme_void()

  q1 <- x$covariate_info[
    order(x$covariate_info$r_squared, decreasing = TRUE, na.last = TRUE),
    , drop = FALSE
  ]
  q1$row <- rev(seq_len(nrow(q1)))
  table_columns <- c(1, 2.5, 3.8, 5)
  q1_labels <- rbind(
    data.frame(row = max(q1$row) + 1, column = table_columns,
               label = c("Covariate", "Type", "R squared", "Spearman rho"),
               header = TRUE),
    data.frame(
      row = rep(q1$row, each = 4),
      column = rep(table_columns, times = nrow(q1)),
      label = as.vector(t(cbind(
        q1$covariate,
        q1$type,
        .entropy_format_number(q1$r_squared),
        .entropy_format_number(q1$spearman_rho)
      ))),
      header = FALSE
    )
  )
  q1_table <- ggplot2::ggplot(q1_labels, ggplot2::aes(x = column, y = row)) +
    ggplot2::geom_text(
      ggplot2::aes(label = label, fontface = ifelse(header, "bold", "plain")),
      hjust = 0, size = 3
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0.8, 6.3), breaks = NULL, expand = c(0, 0)
    ) +
    ggplot2::labs(title = "Q1: global association") +
    ggplot2::theme_void() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 10),
                   legend.position = "none")

  if (is.null(x$joint)) {
    joint_panel <- ggplot2::ggplot() +
      ggplot2::annotate(
        "text", x = 0.5, y = 0.5,
        label = "Joint model not computed\nOnly one covariate was available.",
        size = 3.3
      ) +
      ggplot2::labs(title = "Q5: joint model") +
      ggplot2::theme_void() +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 10))
  } else {
    incremental <- data.frame(
      covariate = factor(
        names(x$joint$incremental_r_squared),
        levels = rev(names(sort(x$joint$incremental_r_squared)))
      ),
      value = unname(x$joint$incremental_r_squared)
    )
    joint_panel <- ggplot2::ggplot(
      incremental, ggplot2::aes(x = value, y = covariate)
    ) +
      ggplot2::geom_col(fill = "#009E73", width = 0.7) +
      ggplot2::labs(
        title = sprintf("Q5: joint R squared = %.3f", x$joint$joint_r_squared),
        x = "Incremental R squared", y = NULL,
        caption = paste(strwrap(x$joint$note, width = 58), collapse = "\n")
      ) +
      base_theme +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 10),
        plot.caption = ggplot2::element_text(size = 6.5, hjust = 0)
      )
  }

  flag_colours <- c(
    "not elevated" = "grey80",
    "explained" = "#D55E00",
    "cannot tell" = "#F0E442",
    "not explained" = "#0072B2"
  )
  chunks <- split(ordering, ceiling(seq_along(ordering) / 12L))
  lapply(seq_along(chunks), function(page_index) {
    chunk <- chunks[[page_index]]
    flag_data <- x$cluster_table[x$cluster_table$cluster %in% chunk, , drop = FALSE]
    flag_data$cluster <- factor(flag_data$cluster, levels = rev(chunk))
    flag_data$covariate <- factor(
      flag_data$covariate, levels = x$covariate_info$covariate
    )
    flag_data$flag <- factor(
      flag_data$flag,
      levels = c("not elevated", "explained", "cannot tell", "not explained")
    )
    flag_data$label_colour <- ifelse(
      flag_data$flag %in% c("explained", "not explained"), "white", "black"
    )
    flag_data$tile_label <- paste0(
      "m = ", .entropy_format_number(flag_data$matched_diff_std, digits = 2),
      ifelse(!is.na(flag_data$reversed) & flag_data$reversed, "*", "")
    )
    flag_matrix <- ggplot2::ggplot(
      flag_data,
      ggplot2::aes(x = covariate, y = cluster, fill = flag)
    ) +
      ggplot2::geom_tile(colour = "white", linewidth = 0.4, show.legend = TRUE) +
      ggplot2::geom_text(
        ggplot2::aes(label = tile_label, colour = label_colour), size = 2.6
      ) +
      ggplot2::scale_fill_manual(values = flag_colours, drop = FALSE) +
      ggplot2::scale_colour_identity() +
      ggplot2::coord_fixed(ratio = 0.62) +
      ggplot2::labs(
        title = sprintf(
          "Q3: matched comparison flags (page %d of %d)",
          page_index, length(chunks)
        ),
        subtitle = "Each tile shows matched difference (m); * marks reversed matching",
        x = NULL, y = NULL, fill = NULL
      ) +
      base_theme +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 10),
        axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
        axis.text.y = ggplot2::element_text(hjust = 1),
        panel.grid = ggplot2::element_blank(),
        legend.position = "bottom"
      )

    if (page_index == 1L) {
      patchwork::wrap_plots(
        header, q1_table, joint_panel, flag_matrix, .entropy_glossary_plot(),
        design = "AAAA\nBBCC\nDDDD\nEEEE",
        heights = c(0.65, 1.45, 3.4, 2.4)
      )
    } else {
      continuation_header <- ggplot2::ggplot() +
        ggplot2::annotate(
          "text", x = 0, y = 1,
          label = "Entropy diagnostic summary",
          hjust = 0, vjust = 1, fontface = "bold", size = 5
        ) +
        ggplot2::annotate(
          "text", x = 0, y = 0.45,
          label = sprintf("Matched comparison flags, continued (page %d of %d)",
                          page_index, length(chunks)),
          hjust = 0, vjust = 1, size = 3.1
        ) +
        ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) + ggplot2::theme_void()
      patchwork::wrap_plots(
        continuation_header, flag_matrix, .entropy_glossary_plot(),
        ncol = 1, heights = c(0.55, 3.4, 2.4)
      )
    }
  })
}


.entropy_covariate_pages <- function(x, covariate_name) {
  info <- x$covariate_info[x$covariate_info$covariate == covariate_name, , drop = FALSE]
  display_type <- if (identical(covariate_name, "mean_knn_dist")) {
    "density"
  } else {
    info$type
  }
  cluster_result <- x$cluster_table[
    x$cluster_table$covariate == covariate_name, , drop = FALSE
  ]
  bin_result <- x$bin_table[x$bin_table$covariate == covariate_name, , drop = FALSE]
  usable_imbalance <- bin_result$imbalance[
    bin_result$usable & is.finite(bin_result$imbalance)
  ]
  median_usable_imbalance <- if (length(usable_imbalance)) {
    stats::median(usable_imbalance)
  } else {
    NA_real_
  }
  retained_not_explained <- cluster_result$retained[
    cluster_result$flag == "not explained" & is.finite(cluster_result$retained)
  ]
  median_retained_not_explained <- if (length(retained_not_explained)) {
    stats::median(retained_not_explained)
  } else {
    NA_real_
  }
  plot_data <- data.frame(
    cluster = x$cell_data$cluster,
    entropy = x$cell_data$entropy,
    covariate_value = x$cell_data[[covariate_name]],
    stringsAsFactors = FALSE
  )
  ordering <- .entropy_cluster_order(x)
  plot_data$cluster <- factor(plot_data$cluster, levels = ordering)
  cluster_result$cluster <- factor(cluster_result$cluster, levels = ordering)
  bin_result$cluster <- factor(bin_result$cluster, levels = ordering)
  base_theme <- ggplot2::theme_bw(base_size = 10)
  global_entropy_median <- stats::median(plot_data$entropy)
  global_covariate_median <- stats::median(plot_data$covariate_value)

  page1 <- ggplot2::ggplot(
    plot_data, ggplot2::aes(x = covariate_value, y = entropy)
  ) +
    ggplot2::geom_point(alpha = 0.18, size = 0.55, colour = "grey35") +
    ggplot2::geom_line(
      data = x$loess_curves[[covariate_name]],
      ggplot2::aes(x = x, y = y),
      inherit.aes = FALSE,
      linewidth = 1,
      colour = "#D55E00",
      na.rm = TRUE
    ) +
    ggplot2::labs(
      title = paste("Q1: does", covariate_name, "matter?"),
      subtitle = sprintf(
        "R squared = %s, Spearman rho = %s",
        .entropy_format_number(info$r_squared),
        .entropy_format_number(info$spearman_rho)
      ),
      x = covariate_name, y = "Entropy"
    ) +
    ggplot2::annotate(
      "text", x = Inf, y = Inf,
      label = paste("Covariate type:", display_type),
      hjust = 1.05, vjust = 1.25, fontface = "bold",
      colour = "#C00000", size = 5
    ) +
    base_theme

  page2_entropy <- ggplot2::ggplot(
    plot_data, ggplot2::aes(x = cluster, y = entropy)
  ) +
    ggplot2::geom_hline(
      yintercept = global_entropy_median,
      linetype = "dashed", colour = "#C00000", linewidth = 0.5
    ) +
    ggplot2::geom_boxplot(fill = "grey80", outlier.alpha = 0.2) +
    ggplot2::annotate(
      "text", x = Inf, y = global_entropy_median,
      label = paste0(
        "Global median = ", .entropy_format_number(global_entropy_median, 2)
      ),
      hjust = 1.05, vjust = -0.45, colour = "#C00000", size = 2.7
    ) +
    ggplot2::labs(title = "Raw entropy", x = NULL, y = "Entropy") +
    base_theme +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  page2_covariate <- ggplot2::ggplot(
    plot_data, ggplot2::aes(x = cluster, y = covariate_value)
  ) +
    ggplot2::geom_hline(
      yintercept = global_covariate_median,
      linetype = "dashed", colour = "#C00000", linewidth = 0.5
    ) +
    ggplot2::geom_boxplot(fill = "#56B4E9", outlier.alpha = 0.2) +
    ggplot2::annotate(
      "text", x = Inf, y = global_covariate_median,
      label = paste0(
        "Global median = ", .entropy_format_number(global_covariate_median, 2)
      ),
      hjust = 1.05, vjust = -0.45, colour = "#C00000", size = 2.7
    ) +
    ggplot2::geom_text(
      data = cluster_result,
      ggplot2::aes(
        x = cluster, y = Inf,
        label = paste0("d = ", .entropy_format_number(covariate_delta_std, 2))
      ),
      inherit.aes = FALSE, vjust = 1.2, size = 2.5
    ) +
    ggplot2::labs(
      title = paste("Raw", covariate_name),
      subtitle = "d is the standardized covariate difference versus other cells",
      x = NULL, y = covariate_name
    ) +
    base_theme +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  page2 <- patchwork::wrap_plots(page2_entropy, page2_covariate, ncol = 2) +
    patchwork::plot_annotation(title = "Q2: is each cluster extreme?")

  cluster_chunks <- split(ordering, ceiling(seq_along(ordering) / 6L))
  page3 <- unlist(lapply(seq_along(cluster_chunks), function(page_index) {
    chunk <- cluster_chunks[[page_index]]
    chunk_bins <- bin_result[as.character(bin_result$cluster) %in% chunk, , drop = FALSE]
    chunk_bins$cluster <- factor(as.character(chunk_bins$cluster), levels = chunk)
    chunk_bins$plot_diff <- ifelse(chunk_bins$usable, chunk_bins$diff, NA_real_)
    too_few <- chunk_bins$n_in < x$thresholds$min_bin_n |
      chunk_bins$n_out < x$thresholds$min_bin_n
    top <- ggplot2::ggplot(
      chunk_bins,
      ggplot2::aes(x = bin, y = plot_diff, group = cluster)
    ) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey45") +
      ggplot2::geom_line(colour = "grey60", linewidth = 0.35, na.rm = TRUE) +
      ggplot2::geom_point(
        ggplot2::aes(size = n_in), colour = "#0072B2", na.rm = TRUE
      ) +
      ggplot2::facet_wrap(~cluster, scales = "free_y", ncol = 3, drop = FALSE) +
      ggplot2::labs(
        title = sprintf(
          "Q3: matched entropy difference for %s (page %d of %d)",
          covariate_name, page_index, length(cluster_chunks)
        ),
        subtitle = sprintf(
          paste0(
            "Only usable bins are drawn; median usable-bin imbalance = %s; ",
            "median retained among not explained clusters = %s"
          ),
          .entropy_format_number(median_usable_imbalance, 3),
          .entropy_format_number(median_retained_not_explained, 2)
        ),
        x = "Quantile bin", y = "Median entropy difference",
        size = "Cells in cluster"
      ) + base_theme +
      ggplot2::theme(legend.position = "bottom")

    status_data <- chunk_bins
    status_data$cluster <- factor(
      as.character(status_data$cluster), levels = rev(chunk)
    )
    status_data$bin_status <- ifelse(
      status_data$usable,
      "usable",
      ifelse(too_few, "too few cells", "imbalanced")
    )
    status_data$bin_status <- factor(
      status_data$bin_status,
      levels = c("usable", "too few cells", "imbalanced")
    )
    status_panel <- ggplot2::ggplot(
      status_data, ggplot2::aes(x = bin, y = cluster, fill = bin_status)
    ) +
      ggplot2::geom_tile(
        colour = "white", linewidth = 0.3, show.legend = TRUE
      ) +
      ggplot2::scale_fill_manual(
        values = c(
          "usable" = "#0072B2",
          "too few cells" = "grey75",
          "imbalanced" = "#E69F00"
        ),
        drop = FALSE
      ) +
      ggplot2::guides(fill = ggplot2::guide_legend(
        override.aes = list(fill = c("#0072B2", "grey75", "#E69F00"))
      )) +
      ggplot2::labs(
        title = "Bin matching status",
        subtitle = sprintf(
          "Imbalance limit = %.2f; median usable-bin imbalance = %s",
          x$thresholds$max_bin_imbalance,
          .entropy_format_number(median_usable_imbalance, 3)
        ),
        x = "Quantile bin", y = NULL, fill = "Bin status"
      ) +
      base_theme +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        legend.position = "bottom"
      )

    matched_page <- patchwork::wrap_plots(
      top, status_panel, ncol = 1, heights = c(1, 0.48)
    )

    background_index <- unique(round(seq(
      1, nrow(plot_data), length.out = min(5000L, nrow(plot_data))
    )))
    facet_data <- do.call(rbind, lapply(chunk, function(cluster_name) {
      use <- sort(unique(c(
        background_index,
        which(as.character(plot_data$cluster) == cluster_name)
      )))
      data.frame(
        facet_cluster = cluster_name,
        cluster = plot_data$cluster[use],
        entropy = plot_data$entropy[use],
        covariate_value = plot_data$covariate_value[use],
        highlighted = as.character(plot_data$cluster[use]) == cluster_name
      )
    }))
    facet_data$facet_cluster <- factor(facet_data$facet_cluster, levels = chunk)
    usable_ranges <- do.call(rbind, lapply(chunk, function(cluster_name) {
      rows <- chunk_bins[
        as.character(chunk_bins$cluster) == cluster_name & chunk_bins$usable,
        , drop = FALSE
      ]
      data.frame(
        facet_cluster = cluster_name,
        xmin = if (nrow(rows)) min(rows$bin_lower) else NA_real_,
        xmax = if (nrow(rows)) max(rows$bin_upper) else NA_real_
      )
    }))
    usable_ranges$facet_cluster <- factor(usable_ranges$facet_cluster, levels = chunk)
    bottom <- ggplot2::ggplot(
      facet_data, ggplot2::aes(x = covariate_value, y = entropy)
    ) +
      ggplot2::geom_rect(
        data = usable_ranges,
        ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
        inherit.aes = FALSE,
        fill = "#009E73", alpha = 0.12, na.rm = TRUE
      ) +
      ggplot2::geom_point(
        data = facet_data[!facet_data$highlighted, , drop = FALSE],
        colour = "grey75", alpha = 0.25, size = 0.35
      ) +
      ggplot2::geom_point(
        data = facet_data[facet_data$highlighted, , drop = FALSE],
        colour = "#D55E00", alpha = 0.65, size = 0.55
      ) +
      ggplot2::facet_wrap(~facet_cluster, ncol = 3, drop = FALSE) +
      ggplot2::labs(
        title = sprintf(
          "Q3: cluster overlap for %s (page %d of %d)",
          covariate_name, page_index, length(cluster_chunks)
        ),
        subtitle = "Green shading marks the usable covariate range",
        x = covariate_name, y = "Entropy"
      ) + base_theme
    list(matched_page, bottom)
  }), recursive = FALSE)

  facet_labels <- vapply(ordering, function(cluster_name) {
    row <- cluster_result[as.character(cluster_result$cluster) == cluster_name, , drop = FALSE]
    paste0(
      cluster_name, "\nslope = ", .entropy_format_number(row$within_slope_std, 2),
      ", R2 = ", .entropy_format_number(row$within_r_squared, 2)
    )
  }, character(1))
  names(facet_labels) <- ordering
  plot_data$facet_label <- factor(
    facet_labels[as.character(plot_data$cluster)], levels = facet_labels
  )
  line_data <- do.call(rbind, lapply(ordering, function(cluster_name) {
    row <- cluster_result[as.character(cluster_result$cluster) == cluster_name, , drop = FALSE]
    data.frame(
      facet_label = facet_labels[[cluster_name]],
      slope = row$within_slope,
      intercept = row$mean_entropy - row$within_slope * row$mean_covariate
    )
  }))
  line_data$facet_label <- factor(line_data$facet_label, levels = facet_labels)
  page4 <- lapply(seq_along(cluster_chunks), function(page_index) {
    chunk <- cluster_chunks[[page_index]]
    chunk_labels <- facet_labels[chunk]
    chunk_data <- plot_data[as.character(plot_data$cluster) %in% chunk, , drop = FALSE]
    chunk_data$facet_label <- factor(
      facet_labels[as.character(chunk_data$cluster)], levels = chunk_labels
    )
    chunk_lines <- line_data[as.character(line_data$facet_label) %in% chunk_labels, , drop = FALSE]
    chunk_lines$facet_label <- factor(
      as.character(chunk_lines$facet_label), levels = chunk_labels
    )
    ggplot2::ggplot(
      chunk_data, ggplot2::aes(x = covariate_value, y = entropy)
    ) +
      ggplot2::geom_point(alpha = 0.25, size = 0.45, colour = "grey35") +
      ggplot2::geom_abline(
        data = chunk_lines,
        ggplot2::aes(slope = slope, intercept = intercept),
        inherit.aes = FALSE,
        colour = "#D55E00", linewidth = 0.8, na.rm = TRUE
      ) +
      ggplot2::facet_wrap(~facet_label, scales = "free", ncol = 3, drop = FALSE) +
      ggplot2::labs(
        title = sprintf(
          "Q4: within-cluster association for %s (page %d of %d)",
          covariate_name, page_index, length(cluster_chunks)
        ),
        subtitle = "Facet labels show standardized slope and within-cluster R squared",
        x = covariate_name, y = "Entropy"
      ) + base_theme
  })

  verdict_pages <- .entropy_diagnostic_verdict_pages(cluster_result, covariate_name)

  caveat_text <- .entropy_caveat_text(display_type, info$intervention)
  caveat_lines <- unlist(lapply(caveat_text, function(paragraph) {
    c(strwrap(paragraph, width = 105), "")
  }))
  caveat_label <- paste(caveat_lines, collapse = "\n")
  page5 <- ggplot2::ggplot() +
    ggplot2::annotate(
      "text", x = 0, y = 1, label = paste("Caveats for", covariate_name),
      hjust = 0, vjust = 1, fontface = "bold", size = 5
    ) +
    ggplot2::annotate(
      "text", x = 0, y = 0.91, label = caveat_label,
      hjust = 0, vjust = 1, size = 3.05, lineheight = 1.08
    ) +
    ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) + ggplot2::theme_void()

  c(list(page1, page2), page3, page4, verdict_pages, list(page5))
}


.entropy_diagnostic_verdict_pages <- function(cluster_result, covariate_name) {
  flag_levels <- c("not elevated", "cannot tell", "explained", "not explained")
  order_index <- order(
    match(cluster_result$flag, flag_levels),
    -cluster_result$matched_diff_std,
    na.last = TRUE
  )
  table_data <- cluster_result[order_index, , drop = FALSE]
  chunks <- split(seq_len(nrow(table_data)), ceiling(seq_len(nrow(table_data)) / 20L))
  columns <- c(0.2, 0.95, 1.7, 2.5, 3.35, 4.2, 5.05, 5.85, 6.65, 7.4, 8.15, 9.35)
  headers <- c(
    "Cluster", "n", "Direction", "Entropy\ndelta", "Covariate\ndelta",
    "Matched\ndiff", "Retained\nelev./def.", "Support", "Usable\nbins",
    "Within\nR2", "Reversed", "Flag"
  )

  lapply(seq_along(chunks), function(page_index) {
    rows <- table_data[chunks[[page_index]], , drop = FALSE]
    rows$y <- rev(seq_len(nrow(rows)))
    labels <- cbind(
      as.character(rows$cluster),
      as.character(rows$n_cells),
      as.character(rows$direction),
      .entropy_format_number(rows$entropy_delta_std, 2),
      .entropy_format_number(rows$covariate_delta_std, 2),
      .entropy_format_number(rows$matched_diff_std, 2),
      .entropy_format_number(rows$retained, 2),
      .entropy_format_number(rows$support, 2),
      as.character(rows$n_usable_bins),
      .entropy_format_number(rows$within_r_squared, 2),
      ifelse(is.na(rows$reversed), "NA", ifelse(rows$reversed, "yes", "no")),
      as.character(rows$flag)
    )
    text_data <- data.frame(
      x = rep(columns, times = nrow(rows)),
      y = rep(rows$y, each = length(columns)),
      label = as.vector(t(labels)),
      stringsAsFactors = FALSE
    )
    header_data <- data.frame(
      x = columns, y = max(rows$y) + 1.2, label = headers,
      stringsAsFactors = FALSE
    )
    rows$retained_fill <- pmin(1, pmax(0, rows$retained))
    rows$retained_label_colour <- ifelse(
      is.finite(rows$retained_fill) & rows$retained_fill < 0.35,
      "white", "black"
    )
    rows$retained_label <- .entropy_format_number(rows$retained, 2)
    cannot_tell <- rows[rows$flag == "cannot tell", , drop = FALSE]
    table_plot <- ggplot2::ggplot() +
      ggplot2::geom_tile(
        data = rows,
        ggplot2::aes(x = 5.05, y = y, fill = retained_fill),
        width = 0.78, height = 0.82
      ) +
      ggplot2::geom_tile(
        data = rows,
        ggplot2::aes(x = 9.35, y = y),
        fill = "grey98", colour = "grey90", width = 1.25, height = 0.82
      ) +
      ggplot2::geom_tile(
        data = cannot_tell,
        ggplot2::aes(x = 9.35, y = y),
        fill = "#FFF4B3", colour = "grey90", width = 1.25, height = 0.82
      ) +
      ggplot2::geom_text(
        data = text_data[!text_data$x %in% c(5.05, 9.35), , drop = FALSE],
        ggplot2::aes(x = x, y = y, label = label),
        size = 2.35
      ) +
      ggplot2::geom_text(
        data = rows,
        ggplot2::aes(
          x = 5.05, y = y, label = retained_label,
          colour = retained_label_colour
        ),
        size = 2.5
      ) +
      ggplot2::geom_text(
        data = rows,
        ggplot2::aes(x = 9.35, y = y, label = flag),
        colour = "black", size = 2.15
      ) +
      ggplot2::geom_text(
        data = header_data,
        ggplot2::aes(x = x, y = y, label = label),
        fontface = "bold", size = 2.25
      ) +
      ggplot2::scale_fill_gradient(
        low = "#D55E00", high = "white", limits = c(0, 1),
        breaks = c(0, 0.5, 1), na.value = "grey95",
        name = "Retained fraction",
        guide = ggplot2::guide_colourbar(
          direction = "horizontal", title.position = "top"
        )
      ) +
      ggplot2::scale_colour_identity() +
      ggplot2::coord_cartesian(xlim = c(0, 10.05), clip = "off") +
      ggplot2::labs(
        title = sprintf(
          "Cluster verdicts for %s (page %d of %d)",
          covariate_name, page_index, length(chunks)
        )
      ) +
      ggplot2::theme_void() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 12),
        legend.position = "bottom",
        legend.justification = "left",
        plot.margin = ggplot2::margin(12, 10, 8, 10)
      )
    patchwork::wrap_plots(
      table_plot, .entropy_glossary_plot(compact = TRUE, include_global = FALSE),
      ncol = 1, heights = c(3.7, 2.2)
    )
  })
}


.entropy_glossary_plot <- function(compact = FALSE, include_global = TRUE) {
  definitions <- list(
    "Matched difference (m)" = paste(
      "The cluster's median entropy minus the median entropy of other cells with",
      "similar covariate values, averaged across bins and divided by the entropy",
      "standard deviation. Positive means the cluster is still higher. Around 0.1",
      "is small and around 0.5 is large."
    ),
    "Support (s)" = paste(
      "The fraction of cluster cells that could be compared. Cells are excluded",
      "where too few other cells exist or the two sides are not balanced in the",
      "covariate. Support of 1 means every cell was comparable."
    ),
    "Retained fraction" = paste(
      "The fraction of the cluster's entropy difference that survives matching.",
      "For elevated clusters this is retained elevation; for depleted clusters it",
      "is retained deficit. Near 1 means little is explained and near 0 means nearly",
      "all is explained. NA means the cluster is not elevated or depleted enough",
      "for the fraction to be interpretable."
    ),
    "R squared" = paste(
      "The fraction of cell-to-cell entropy variation accounted for by this covariate",
      "across the whole dataset. This differs from matching, and the two can disagree."
    ),
    "Incremental R squared" = paste(
      "How much this covariate adds when entered last, after all other covariates."
    )
  )
  if (!include_global) {
    definitions <- definitions[c(
      "Matched difference (m)", "Support (s)", "Retained fraction"
    )]
  }

  if (include_global) {
    positions <- data.frame(
      x = c(0, 0.52, 0, 0.52, 0.52),
      y = c(0.84, 0.84, 0.43, 0.43, 0.13),
      width = c(52, 52, 52, 52, 70)
    )
  } else {
    positions <- data.frame(
      x = rep(0, 3), y = c(0.82, 0.53, 0.24), width = rep(118, 3)
    )
  }

  plot <- ggplot2::ggplot() +
    ggplot2::annotate(
      "text", x = 0, y = 1, label = "Glossary",
      hjust = 0, vjust = 1, fontface = "bold",
      size = if (compact) 3.5 else 3.8
    ) +
    ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) + ggplot2::theme_void()

  for (index in seq_along(definitions)) {
    plot <- plot +
      ggplot2::annotate(
        "text", x = positions$x[index], y = positions$y[index],
        label = names(definitions)[index], hjust = 0, vjust = 1,
        fontface = "bold", size = if (compact) 2.75 else 2.9
      ) +
      ggplot2::annotate(
        "text", x = positions$x[index], y = positions$y[index] - 0.07,
        label = paste(strwrap(
          definitions[[index]], width = positions$width[index]
        ), collapse = "\n"),
        hjust = 0, vjust = 1, lineheight = 1.05,
        size = if (compact) 2.25 else 2.3
      )
  }
  plot
}


.entropy_format_number <- function(value, digits = 3) {
  ifelse(
    is.finite(value),
    formatC(value, format = "f", digits = digits),
    "NA"
  )
}


.entropy_caveat_text <- function(type, intervention) {
  type_text <- switch(
    type,
    technical = paste(
      "This covariate is declared technical. A flag of `explained` is bad news",
      "for a biological interpretation of entropy in this cluster."
    ),
    biological = paste(
      "This covariate is declared biological. A flag of `explained` is ambiguous,",
      "because adjusting for a biological covariate can remove real signal as easily",
      "as it removes an artifact."
    ),
    unknown = paste(
      "This covariate's type is not declared. Declare it as technical or biological",
      "to get a sharper interpretation."
    ),
    density = paste(
      "Neighbourhood density is neither purely technical nor purely biological.",
      "Sparse neighbourhoods can arise because few cells of that type were sampled,",
      "which is technical, or because the population is genuinely diffuse in expression",
      "space, which may be the biology you are trying to measure. Treat a flag of",
      "`explained` for density as ambiguous."
    )
  )
  intervention_text <- switch(
    intervention,
    downsample_counts = paste(
      "Follow up: downsample all cells to a common total UMI count, recompute entropy,",
      "and check whether the cluster ranking survives."
    ),
    subsample_cells = paste(
      "Follow up: subsample cells to equalise neighbourhood density, recompute entropy,",
      "and check whether the cluster ranking survives."
    ),
    remove_genes = paste(
      "Follow up: remove the genes that define this score, recompute entropy, and compare",
      "the rank correlation of cluster means before and after."
    ),
    none = paste(
      "No follow up test by recomputation is available for this covariate. The observational",
      "result above is all the data can give you."
    )
  )
  c(
    "`not elevated`: this cluster does not have high entropy, so the question of what explains its entropy does not apply.",
    "`cannot tell`: this cluster occupies a range of the covariate where too few other cells exist, so there is no fair comparison group. No conclusion can be drawn.",
    "`explained`: at matched covariate values, this cluster is not still higher than other cells. The covariate accounts for its elevation, and the direction may be reversed.",
    "`not explained`: at matched covariate values, this cluster is specifically still higher than other cells. The covariate does not account for it.",
    paste(
      "A covariate can be a confounder, where it mechanically raises entropy and should be",
      "treated as an artifact, or a mediator, where a genuinely plastic state causes both",
      "the covariate and the entropy. The mathematics is identical in both cases and this",
      "diagnostic cannot distinguish them. Only reasoning outside the data can."
    ),
    type_text,
    intervention_text,
    "This is a screening step. It does not correct entropy, establish causality, or prove plasticity."
  )
}


utils::globalVariables(c(
  "bin", "bin_status", "cluster", "column", "covariate", "covariate_delta_std",
  "covariate_value", "direction", "entropy", "flag", "imbalance_label",
  "intercept", "label", "label_colour", "n_in", "plot_diff", "retained_fill",
  "retained_label", "retained_label_colour", "slope", "status",
  "strip_label_colour", "tile_label", "value", "x", "xmax", "xmin", "y"
))
