#' Write a PDF report for an entropy recompute check
#'
#' `entropy_recompute_report()` renders an existing `entropy_recompute` object.
#' It does not rerun an intervention or change any result.
#'
#' @param x An `entropy_recompute` object returned by
#'   [entropy_recompute_check()].
#' @param outdir Output directory. It is created recursively when needed.
#' @param prefix Optional string prepended to every output file name.
#' @param width PDF width in inches.
#' @param height PDF height in inches.
#' @param overwrite Logical. If `FALSE`, existing target files cause an error.
#' @param write_tables Logical. If `TRUE`, also write the result tables as CSV
#'   files.
#'
#' @return A named character vector of normalized output paths, invisibly.
#' @export
entropy_recompute_report <- function(
    x,
    outdir = "entropy_recompute",
    prefix = NULL,
    width = 9,
    height = 7,
    overwrite = FALSE,
    write_tables = FALSE) {
  if (!inherits(x, "entropy_recompute")) {
    stop("`x` must be an entropy_recompute object.", call. = FALSE)
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
  if (!dir.exists(outdir) &&
      !dir.create(outdir, recursive = TRUE, showWarnings = FALSE)) {
    stop(sprintf("Could not create output directory '%s'.", outdir),
         call. = FALSE)
  }
  prefix <- if (is.null(prefix)) "" else prefix
  file_names <- c(report = paste0(prefix, "entropy_recompute.pdf"))
  if (write_tables) {
    file_names <- c(
      file_names,
      cluster_summary = paste0(prefix, "cluster_summary.csv"),
      cluster_verdict = paste0(prefix, "cluster_verdict.csv"),
      repeat_summary = paste0(prefix, "repeat_summary.csv"),
      cell_data = paste0(prefix, "cell_data.csv")
    )
    if (!is.null(x$density_check)) {
      file_names <- c(
        file_names,
        density_check = paste0(prefix, "density_check.csv")
      )
    }
  }
  paths <- file.path(outdir, unname(file_names))
  names(paths) <- names(file_names)
  existing <- paths[file.exists(paths)]
  if (!overwrite && length(existing)) {
    stop(sprintf("Refusing to overwrite existing report file: %s.",
                 existing[[1L]]), call. = FALSE)
  }

  pages <- .entropy_recompute_pages(x)
  .entropy_write_pdf(paths[["report"]], width, height, pages)
  if (write_tables) {
    for (table_name in intersect(names(paths), c(
      "cluster_summary", "cluster_verdict", "repeat_summary", "cell_data",
      "density_check"
    ))) {
      utils::write.csv(x[[table_name]], paths[[table_name]], row.names = FALSE)
    }
  }
  paths <- vapply(paths, normalizePath, character(1), mustWork = TRUE)
  invisible(paths)
}


.entropy_recompute_pages <- function(x) {
  colours <- c(
    survives = "#009E73",
    collapses = "#D55E00",
    improves = "#0072B2",
    unstable = "#E69F00"
  )
  rank_values <- x$repeat_summary$rank_correlation
  provenance_lines <- c(
    sprintf("Mechanism: %s", x$provenance$mechanism),
    sprintf("Repeats: %d, seed: %d", x$provenance$n_repeats, x$provenance$seed),
    sprintf("Cells: %d before, %s after per repeat",
            x$provenance$n_cells_before,
            paste(x$provenance$n_cells_after, collapse = ", ")),
    sprintf("Neighbours recomputed: %s",
            if (x$provenance$recompute_neighbors) "yes" else "no"),
    sprintf("Assay: %s, k: %s, reduction: %s, dims: %s",
            x$provenance$entropy_args$assay,
            .entropy_recompute_arg(x$provenance$entropy_args$k, 15),
            .entropy_recompute_arg(x$provenance$entropy_args$reduction, "automatic"),
            paste(.entropy_recompute_arg(x$provenance$entropy_args$dims, 20),
                  collapse = ", "))
  )
  if (x$provenance$mechanism == "downsample_counts") {
    provenance_lines <- c(
      provenance_lines,
      sprintf("Target depth: %d, cells below target: %d, dropped: %s",
              x$provenance$target_depth, x$provenance$n_below_target,
              if (x$provenance$drop_below_target) "yes" else "no")
    )
  } else if (x$provenance$mechanism == "subsample_cells") {
    provenance_lines <- c(
      provenance_lines,
      sprintf("Cells requested per group: %d", x$provenance$cells_per_group)
    )
  } else {
    provenance_lines <- c(
      provenance_lines,
      sprintf("Genes removed: %d, requested genes missing: %d",
              x$provenance$n_genes_removed, x$provenance$n_genes_missing),
      "Absolute entropy deltas are not comparable for gene removal."
    )
  }
  if (length(x$provenance$warnings)) {
    provenance_lines <- c(
      provenance_lines,
      paste("Recorded warnings:", paste(x$provenance$warnings, collapse = " | "))
    )
  }
  provenance_text <- paste(unlist(lapply(provenance_lines, strwrap, width = 105)),
                           collapse = "\n")
  page1 <- ggplot2::ggplot() +
    ggplot2::annotate(
      "text", x = 0, y = 1, label = "Entropy recompute check",
      hjust = 0, vjust = 1, fontface = "bold", size = 5.5
    ) +
    ggplot2::annotate(
      "text", x = 0.5, y = 0.76,
      label = sprintf("%.3f", mean(rank_values)),
      hjust = 0.5, vjust = 0.5, fontface = "bold", size = 15,
      colour = "#0072B2"
    ) +
    ggplot2::annotate(
      "text", x = 0.5, y = 0.63,
      label = sprintf(
        "Mean Spearman cluster-rank correlation\nrange %.3f to %.3f",
        min(rank_values), max(rank_values)
      ),
      hjust = 0.5, vjust = 0.5, size = 4
    ) +
    ggplot2::annotate(
      "text", x = 0, y = 0.48, label = provenance_text,
      hjust = 0, vjust = 1, size = 3, lineheight = 1.1
    ) +
    ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) + ggplot2::theme_void()

  verdict <- x$cluster_verdict[, c("cluster", "verdict"), drop = FALSE]
  cluster_medians <- stats::aggregate(
    x$cluster_summary[, c("median_before", "median_after")],
    list(cluster = x$cluster_summary$cluster), stats::median
  )
  cluster_medians <- merge(cluster_medians, verdict, by = "cluster", sort = FALSE)
  paired <- rbind(
    data.frame(
      cluster = cluster_medians$cluster, stage = "Before",
      entropy = cluster_medians$median_before, verdict = cluster_medians$verdict
    ),
    data.frame(
      cluster = cluster_medians$cluster, stage = "After",
      entropy = cluster_medians$median_after, verdict = cluster_medians$verdict
    )
  )
  paired$stage <- factor(paired$stage, levels = c("Before", "After"))
  page2 <- ggplot2::ggplot(
    paired, ggplot2::aes(x = stage, y = entropy, group = cluster, colour = verdict)
  ) +
    ggplot2::geom_line(linewidth = 0.75, alpha = 0.8) +
    ggplot2::geom_point(size = 2.2) +
    ggplot2::geom_text(
      data = paired[paired$stage == "After", , drop = FALSE],
      ggplot2::aes(label = cluster), hjust = -0.15, size = 2.6,
      show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(values = colours, drop = FALSE) +
    ggplot2::scale_x_discrete(expand = ggplot2::expansion(mult = c(0.15, 0.32))) +
    ggplot2::labs(
      title = "Cluster median entropy before and after",
      subtitle = "One line per cluster; values are medians across repeats",
      x = NULL, y = "Median entropy", colour = "Verdict"
    ) +
    ggplot2::theme_bw(base_size = 10)

  rank_data <- merge(
    x$cluster_summary, verdict, by = "cluster", sort = FALSE
  )
  page3 <- ggplot2::ggplot(
    rank_data,
    ggplot2::aes(x = rank_before, y = rank_after, colour = verdict, label = cluster)
  ) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    ggplot2::geom_point(size = 2.4) +
    ggplot2::geom_text(nudge_y = 0.18, size = 2.3, show.legend = FALSE) +
    ggplot2::facet_wrap(~repeat_id) +
    ggplot2::scale_colour_manual(values = colours, drop = FALSE) +
    ggplot2::scale_x_reverse(expand = ggplot2::expansion(mult = 0.08)) +
    ggplot2::scale_y_reverse(expand = ggplot2::expansion(mult = 0.08)) +
    ggplot2::labs(
      title = "Cluster rank before versus after",
      subtitle = "Rank 1 is the highest median entropy; the dashed line means no change",
      x = "Rank before", y = "Rank after", colour = "Verdict"
    ) +
    ggplot2::theme_bw(base_size = 10)

  verdict_pages <- .entropy_recompute_verdict_pages(x, colours)
  density_page <- if (is.null(x$density_check)) list() else {
    density <- merge(x$density_check, verdict, by = "cluster", sort = FALSE)
    density_long <- rbind(
      data.frame(
        cluster = density$cluster, repeat_id = density$repeat_id,
        stage = "Before", distance = density$mean_knn_before,
        verdict = density$verdict
      ),
      data.frame(
        cluster = density$cluster, repeat_id = density$repeat_id,
        stage = "After", distance = density$mean_knn_after,
        verdict = density$verdict
      )
    )
    density_long$stage <- factor(density_long$stage, levels = c("Before", "After"))
    list(ggplot2::ggplot(
      density_long,
      ggplot2::aes(x = stage, y = distance, group = cluster, colour = verdict)
    ) +
      ggplot2::geom_line(alpha = 0.75) +
      ggplot2::geom_point(size = 1.8) +
      ggplot2::facet_wrap(~repeat_id) +
      ggplot2::scale_colour_manual(values = colours, drop = FALSE) +
      ggplot2::labs(
        title = "Neighbourhood density check",
        subtitle = paste(
          "Mean neighbour distance should become more comparable across clusters.",
          "If it does not, the intervention did not achieve its aim."
        ),
        x = NULL, y = "Mean nearest-neighbour distance", colour = "Verdict"
      ) +
      ggplot2::theme_bw(base_size = 10))
  }
  caveat <- .entropy_recompute_caveat_page(x$provenance$mechanism)
  c(list(page1, page2, page3), verdict_pages, density_page, list(caveat))
}


.entropy_recompute_verdict_pages <- function(x, colours) {
  table_data <- x$cluster_verdict
  before_after <- stats::aggregate(
    x$cluster_summary[, c("median_before", "median_after")],
    list(cluster = x$cluster_summary$cluster), stats::median
  )
  if (x$provenance$absolute_delta_comparable) {
    delta <- stats::aggregate(
      x$cluster_summary$delta,
      list(cluster = x$cluster_summary$cluster), stats::median
    )
    before_after$delta <- delta$x[match(before_after$cluster, delta$cluster)]
  } else {
    before_after$delta <- NA_real_
  }
  table_data <- merge(table_data, before_after, by = "cluster", sort = FALSE)
  order_levels <- c("collapses", "unstable", "survives", "improves")
  table_data <- table_data[order(
    match(table_data$verdict, order_levels),
    -table_data$median_rank_change
  ), , drop = FALSE]
  chunks <- split(seq_len(nrow(table_data)), ceiling(seq_len(nrow(table_data)) / 20L))
  lapply(seq_along(chunks), function(page_index) {
    rows <- table_data[chunks[[page_index]], , drop = FALSE]
    rows$y <- rev(seq_len(nrow(rows)))
    columns <- c(0.3, 1.55, 2.75, 3.95, 5.15, 6.15, 7.15, 8.35)
    headers <- c(
      "Cluster", "Median before", "Median after", "Absolute\ndelta",
      "Median rank\nchange", "Minimum", "Maximum", "Verdict"
    )
    delta_label <- if (x$provenance$absolute_delta_comparable) {
      .entropy_format_number(rows$delta, 3)
    } else {
      rep("not comparable", nrow(rows))
    }
    labels <- cbind(
      rows$cluster,
      .entropy_format_number(rows$median_before, 3),
      .entropy_format_number(rows$median_after, 3),
      delta_label,
      .entropy_format_number(rows$median_rank_change, 1),
      .entropy_format_number(rows$min_rank_change, 1),
      .entropy_format_number(rows$max_rank_change, 1),
      rows$verdict
    )
    text_data <- data.frame(
      x = rep(columns, times = nrow(rows)),
      y = rep(rows$y, each = length(columns)),
      label = as.vector(t(labels)),
      stringsAsFactors = FALSE
    )
    headers_data <- data.frame(
      x = columns, y = max(rows$y) + 1.2, label = headers,
      stringsAsFactors = FALSE
    )
    rows$label_colour <- ifelse(rows$verdict == "unstable", "black", "white")
    ggplot2::ggplot() +
      ggplot2::geom_tile(
        data = rows,
        ggplot2::aes(x = 8.35, y = y, fill = verdict),
        width = 1.25, height = 0.82
      ) +
      ggplot2::geom_text(
        data = text_data[text_data$x != 8.35, , drop = FALSE],
        ggplot2::aes(x = x, y = y, label = label), hjust = 0, size = 2.8
      ) +
      ggplot2::geom_text(
        data = rows,
        ggplot2::aes(x = 8.35, y = y, label = verdict, colour = label_colour),
        size = 2.7
      ) +
      ggplot2::geom_text(
        data = headers_data,
        ggplot2::aes(x = x, y = y, label = label),
        hjust = 0, fontface = "bold", size = 2.8
      ) +
      ggplot2::scale_fill_manual(values = colours, guide = "none") +
      ggplot2::scale_colour_identity() +
      ggplot2::coord_cartesian(xlim = c(0, 9.05), clip = "off") +
      ggplot2::labs(title = sprintf(
        "Per-cluster verdicts (page %d of %d)", page_index, length(chunks)
      )) +
      ggplot2::theme_void() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 12),
        plot.margin = ggplot2::margin(15, 15, 15, 15)
      )
  })
}


.entropy_recompute_caveat_page <- function(mechanism) {
  paragraphs <- c(
    paste(
      "This test removes the effect physically rather than statistically, so",
      "surviving it is stronger evidence than surviving a regression adjustment."
    ),
    paste(
      "Surviving this test does not prove plasticity. It only rules out this one",
      "alternative explanation."
    ),
    paste(
      "Downsampling and subsampling are random. Results from a single repeat are",
      "not reliable. Check the range across repeats before drawing a conclusion."
    ),
    paste(
      "`improves` means the cluster rank rose by more than the tolerance in every",
      "repeat. This is a consistent positive result, not instability. `unstable` is",
      "reserved for repeat results that do not agree on a consistent direction."
    )
  )
  if (mechanism == "remove_genes") {
    paragraphs <- c(
      paragraphs,
      paste(
        "Removing genes shifts the entropy scale for every cell. Only the ranking",
        "is comparable before and after, not the absolute values."
      )
    )
  }
  label <- paste(unlist(lapply(paragraphs, function(value) {
    c(strwrap(value, width = 95), "")
  })), collapse = "\n")
  ggplot2::ggplot() +
    ggplot2::annotate(
      "text", x = 0, y = 1, label = "Caveats",
      hjust = 0, vjust = 1, fontface = "bold", size = 5
    ) +
    ggplot2::annotate(
      "text", x = 0, y = 0.86, label = label,
      hjust = 0, vjust = 1, size = 3.5, lineheight = 1.15
    ) +
    ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) + ggplot2::theme_void()
}


.entropy_recompute_arg <- function(value, default) {
  if (is.null(value)) default else value
}


utils::globalVariables(c(
  "distance", "entropy", "label_colour", "rank_after", "rank_before",
  "repeat_id", "stage", "verdict"
))
