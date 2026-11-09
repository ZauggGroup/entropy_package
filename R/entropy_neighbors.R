#' Build or validate a fixed neighbour graph for entropy
#'
#' You do not need this for a normal run. [run_entropy()] builds the graph
#' itself and caches it on the object. Use `entropy_neighbors()` when the same
#' graph must exist in more than one place: saved to disk, reused by a second
#' script, or handed to a depth sweep so every run sees an identical graph.
#'
#' Neighbours never include the cell itself. `k` means `k` other cells.
#'
#' @param seu A Seurat object.
#' @param neighbors Optional. A named list of neighbour positions or cell names,
#'   or a Seurat `Neighbor` object. If `NULL`, the graph is computed from
#'   `reduction`. Positions must be distinct finite integers.
#' @param k Integer, at least 11. Number of non-self neighbours per cell.
#' @param reduction Name of an existing reduction. If `NULL`,
#'   `pca_<neighbor_assay>` is used. The reduction must already exist, this
#'   function never runs PCA.
#' @param dims One integer `n`, meaning `1:n`, or a vector of dimensions.
#' @param neighbor_assay Assay the reduction belongs to. Defaults to the
#'   object's default assay.
#' @param file Optional path. If given, the result is also written with
#'   `saveRDS()`.
#'
#' @return An object of class `entropy_neighbors`: a named list with one
#'   integer vector of neighbour positions per cell, in `colnames(seu)` order.
#' @export
entropy_neighbors <- function(seu,
                              neighbors = NULL,
                              k = 15,
                              reduction = NULL,
                              dims = 20,
                              neighbor_assay = NULL,
                              file = NULL) {

  k_given <- !missing(k)

  if (!inherits(seu, "Seurat")) {
    stop("`seu` must be a Seurat object.", call. = FALSE)
  }
  if (length(k) != 1L || !is.numeric(k) || is.na(k) || k != floor(k) || k < 11) {
    stop("`k` must be one integer of at least 11.", call. = FALSE)
  }
  k <- as.integer(k)
  cell_names <- colnames(seu)
  if (is.null(neighbor_assay)) neighbor_assay <- Seurat::DefaultAssay(seu)

  if (is.null(neighbors)) {
    if (is.null(reduction)) reduction <- paste0("pca_", neighbor_assay)
    if (!reduction %in% names(seu@reductions)) {
      stop(sprintf(paste0(
        "Reduction '%s' not found. Available: %s.\n",
        "entropy_neighbors() never runs PCA. Build the reduction first, ",
        "or let run_entropy() do it."),
        reduction, paste(names(seu@reductions), collapse = ", ")),
        call. = FALSE)
    }
    dims <- .entropy_nn_dims(dims, seu, reduction)
    emb <- Seurat::Embeddings(seu[[reduction]])[, dims, drop = FALSE]
    if (nrow(emb) <= k) {
      stop(sprintf(paste0(
        "Cannot compute %d neighbors from only %d cells: global automatic ",
        "neighboring requires at least k + 1 = %d cells."),
        k, nrow(emb), k + 1L), call. = FALSE)
    }
    nn <- suppressMessages(Seurat::FindNeighbors(
      emb, k.param = k + 1L, return.neighbor = TRUE, verbose = FALSE))
    lst <- .entropy_nn_trim(.entropy_nn_drop_self(nn@nn.idx), k)
    names(lst) <- cell_names
    source <- "computed"
  } else {
    lst <- .entropy_nn_normalise(neighbors, cell_names, k)
    source <- if (inherits(neighbors, "Neighbor")) "seurat" else "custom"
    reduction <- NA_character_
    dims <- NA_integer_
  }

  k <- .entropy_nn_check(
    lst, cell_names, k, is.null(neighbors) || k_given
  )
  if (k < 11L) {
    stop(sprintf(paste0(
      "The supplied neighbour list gives k = %d. The counts entropy ",
      "estimator is undefined below 11 neighbors."), k), call. = FALSE)
  }

  attr(lst, "signature") <- list(
    cells = cell_names, k = k, reduction = reduction, dims = dims,
    neighbor_assay = neighbor_assay, restrict_by = NA_character_
  )
  attr(lst, "source") <- source
  attr(lst, "built") <- Sys.time()
  class(lst) <- "entropy_neighbors"

  if (!is.null(file)) saveRDS(lst, file)
  lst
}


#' @export
print.entropy_neighbors <- function(x, ...) {
  s <- attr(x, "signature")
  cat("entropy_neighbors\n")
  cat("  cells         :", length(x), "\n")
  cat("  k (non-self)  :", s$k, "\n")
  cat("  reduction     :", s$reduction, "\n")
  cat("  dims          :", if (all(is.na(s$dims))) NA else
    paste0(min(s$dims), ":", max(s$dims)), "\n")
  cat("  neighbor_assay:", s$neighbor_assay, "\n")
  cat("  restrict_by   :", s$restrict_by, "\n")
  cat("  source        :", attr(x, "source"), "\n")
  cat("  built         :", format(attr(x, "built")), "\n")
  invisible(x)
}


# Resolve dims the same way run_entropy() does: one n means 1:n, cap to what
# the reduction actually has.
#' @keywords internal
.entropy_nn_dims <- function(dims, seu, reduction) {
  if (length(dims) == 1L) {
    if (is.na(dims) || dims < 1) stop("`dims` must be a positive integer.", call. = FALSE)
    dims <- seq_len(dims)
  }
  dims <- as.integer(dims)
  if (any(is.na(dims)) || any(dims < 1)) {
    stop("`dims` must be positive integers.", call. = FALSE)
  }
  n_avail <- ncol(Seurat::Embeddings(seu[[reduction]]))
  if (max(dims) > n_avail) {
    warning(sprintf("Requested dims up to %d but '%s' has %d. Capping.",
                    max(dims), reduction, n_avail), call. = FALSE)
    dims <- dims[dims <= n_avail]
  }
  if (length(dims) < 1L) {
    stop(sprintf("No valid dims after capping against '%s'.", reduction), call. = FALSE)
  }
  dims
}


# Drop the cell itself from each row.
#' @keywords internal
.entropy_nn_drop_self <- function(nn_idx) {
  n <- nrow(nn_idx)
  out <- vector("list", n)
  for (i in seq_len(n)) {
    v <- as.integer(nn_idx[i, ])
    out[[i]] <- v[v != i]
  }
  out
}


# Turn a user-supplied list or Neighbor object into the standard list form.
#' @keywords internal
.entropy_nn_normalise <- function(neighbors, cell_names, k) {

  if (inherits(neighbors, "entropy_neighbors")) {
    if (!identical(names(neighbors), cell_names)) {
      stop(paste0(
        "These neighbours were built on a different set of cells. ",
        "Neighbour values are positions in the full cell vector, so they ",
        "cannot be reused on a subset. Rebuild them for this object."),
        call. = FALSE)
    }
    lst <- unclass(neighbors)
    attributes(lst) <- NULL
    names(lst) <- cell_names
    return(lst)
  }

  if (inherits(neighbors, "Neighbor")) {
    nn_cells <- methods::slot(neighbors, "cell.names")
    if (!all(cell_names %in% nn_cells)) {
      stop("The Neighbor object does not cover every cell in the object.",
           call. = FALSE)
    }
    idx <- methods::slot(neighbors, "nn.idx")
    idx <- idx[match(cell_names, nn_cells), , drop = FALSE]
    remap <- match(nn_cells, cell_names)
    idx[] <- remap[idx]
    lst <- .entropy_nn_drop_self(idx)
    names(lst) <- cell_names
    return(lst)
  }

  if (!is.list(neighbors) || is.null(names(neighbors))) {
    stop(paste0(
      "`neighbors` must be a named list, a Seurat Neighbor object, or NULL. ",
      "neighbors[['cell1']] should give the neighbours of cell1."),
      call. = FALSE)
  }
  if (anyDuplicated(names(neighbors))) {
    stop("`neighbors` must not contain duplicate cell names.", call. = FALSE)
  }
  missing_cells <- setdiff(cell_names, names(neighbors))
  if (length(missing_cells)) {
    stop(sprintf(
      "The neighbours list is missing %d cells, for example '%s'. Every cell must appear in names(neighbors).",
      length(missing_cells), missing_cells[[1L]]), call. = FALSE)
  }

  lst <- vector("list", length(cell_names))
  names(lst) <- cell_names
  for (i in seq_along(cell_names)) {
    v <- neighbors[[cell_names[[i]]]]
    if (is.character(v)) {
      v <- match(v, cell_names)
      if (any(is.na(v))) {
        stop(sprintf("Some neighbour names for cell '%s' are not in the object.",
                     cell_names[[i]]), call. = FALSE)
      }
    } else if (length(v) &&
               (!is.numeric(v) || any(!is.finite(v)) ||
                any(v != floor(v)))) {
      stop(sprintf(
        "Neighbour positions for cell '%s' must be finite integers.",
        cell_names[[i]]), call. = FALSE)
    }
    lst[[i]] <- as.integer(v)
  }
  lst
}


# Validate a neighbour list and resolve k.
#  * positions must be inside 1:n
#  * no cell may be its own neighbour
#  * every non-empty list must have the SAME size
#  * if k was given explicitly, that size must equal k
# Empty lists are allowed: restrict_by produces them for groups under k + 1
# cells and for cells with NA in the grouping column.
# Returns the resolved k.
#' @keywords internal
.entropy_nn_check <- function(lst, cell_names, k = NULL, k_given = FALSE) {
  n <- length(cell_names)

  all_idx <- unlist(lst, use.names = FALSE)
  if (length(all_idx)) {
    if (anyNA(all_idx) || min(all_idx) < 1L || max(all_idx) > n) {
      stop("Some neighbour positions are outside 1:ncol(seu).", call. = FALSE)
    }
  }

  self_hit <- which(vapply(seq_len(n), function(i) any(lst[[i]] == i), logical(1)))
  if (length(self_hit)) {
    stop(sprintf(paste0(
      "%d cells list themselves as their own neighbour (for example '%s').\n",
      "Supply non-self neighbours only. If you built the list from a Seurat ",
      "Neighbor object, remove the row's own index wherever it occurs, or ",
      "pass the Neighbor object straight to `neighbors` and Entropy removes ",
      "it for you."),
      length(self_hit), cell_names[[self_hit[[1L]]]]), call. = FALSE)
  }

  duplicate_hit <- which(vapply(
    lst, function(x) anyDuplicated(x) > 0L, logical(1)
  ))
  if (length(duplicate_hit)) {
    stop(sprintf(
      "%d cells repeat a neighbour position (for example '%s'). Each of the k neighbours must be a different cell.",
      length(duplicate_hit), cell_names[[duplicate_hit[[1L]]]]),
      call. = FALSE)
  }

  lens <- lengths(lst)
  sizes <- sort(unique(lens[lens > 0L]))
  if (!length(sizes)) {
    if (k_given && !is.null(k)) return(as.integer(k))
    stop(paste0(
      "Every cell has an empty neighbour list, so `k` cannot be inferred. ",
      "Supply `k` explicitly."), call. = FALSE)
  }
  if (length(sizes) > 1L) {
    stop(sprintf(paste0(
      "Neighbour lists have different sizes (%s). Every cell must have the ",
      "same number of non-self neighbours. Empty lists are allowed, anything ",
      "between 1 and the common size is not."),
      paste(sizes, collapse = ", ")), call. = FALSE)
  }
  if (k_given && !identical(as.integer(sizes), as.integer(k))) {
    stop(sprintf(paste0(
      "k = %d was requested but every cell has %d non-self neighbours.\n",
      "Either drop `k` and it is taken from the list, or supply exactly k ",
      "neighbours. A Seurat Neighbor object loses its self entry here, so ",
      "build it with k.param = k + 1."),
      as.integer(k), as.integer(sizes)), call. = FALSE)
  }
  as.integer(sizes)
}


# Keep exactly k candidates after identity-based self removal. A row can still
# contain k + 1 entries when ties or approximate search omit the query cell.
#' @keywords internal
.entropy_nn_trim <- function(lst, k) {
  lapply(lst, utils::head, n = k)
}


# Add an actionable remedy to validation failures caused by pruning a supplied
# graph with restrict_by.
#' @keywords internal
.entropy_nn_check_restricted <- function(lst, cell_names, k) {
  tryCatch(
    .entropy_nn_check(lst, cell_names, k, k_given = TRUE),
    error = function(error) {
      stop(paste0(
        conditionMessage(error), "\n",
        "With `restrict_by`, a supplied graph must already give every cell ",
        "either exactly k within-group neighbours or none. Drop `neighbors` ",
        "(and leave `neighbors_method` unset) to rebuild k neighbours within ",
        "each group."), call. = FALSE)
    }
  )
}

# Name the signature fields that differ, for the "why did it recompute" message.
#' @keywords internal
.entropy_nn_diff <- function(old, new) {
  if (is.null(old)) return("no signature stored")
  changed <- character(0)
  for (f in union(names(old), names(new))) {
    if (!identical(old[[f]], new[[f]])) changed <- c(changed, f)
  }
  if (!length(changed)) return("unknown")
  changed[changed == "cells"] <- "cell set"
  paste(changed, "changed")
}
