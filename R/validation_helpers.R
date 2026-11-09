# Shared scalar validation --------------------------------------------------

.entropy_scalar_logical <- function(value, argument) {
  if (length(value) != 1L || !is.logical(value) || is.na(value)) {
    stop(sprintf("`%s` must be TRUE or FALSE.", argument), call. = FALSE)
  }
  invisible(value)
}


.entropy_report_size <- function(value, argument) {
  if (length(value) != 1L || !is.numeric(value) || !is.finite(value) ||
      value <= 0) {
    stop(sprintf("`%s` must be one positive finite number.", argument),
         call. = FALSE)
  }
  as.numeric(value)
}


.entropy_write_pdf <- function(path, width, height, plots) {
  local({
    grDevices::pdf(path, width = width, height = height, onefile = TRUE)
    on.exit(grDevices::dev.off(), add = TRUE)
    for (plot in plots) print(plot)
    grDevices::dev.off()
    on.exit(NULL)
  })
  invisible(path)
}


.entropy_format_number <- function(value, digits = 3L) {
  ifelse(
    is.finite(value),
    formatC(value, format = "f", digits = digits),
    "NA"
  )
}
