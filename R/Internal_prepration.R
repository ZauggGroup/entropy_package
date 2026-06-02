#' Get neighbor indices for a cell from a custom list or Seurat Neighbor object
#'
#' @param cell_name Cell name (character).
#' @param seu Seurat object.
#' @param neighbors Neighbor definition. Either a named list where
#'   `neighbors[[cell_name]]` returns a vector of neighbor cell names or indices,
#'   or a Seurat Neighbor object.
#' @param neighbors_method Character. One of `"custom"` or `"seurat"`.
#' @return Integer vector of neighbor indices.
#' @keywords internal
get_neighbors <- function(cell_name, seu, neighbors, neighbors_method) {
  if (neighbors_method == "custom") {
    cell_neighbors <- neighbors[[cell_name]]

    if (length(cell_neighbors) == 0) {
      return(integer(0))
    }

    if (is.character(cell_neighbors)) {
      idx <- match(cell_neighbors, colnames(seu))

      if (any(is.na(idx))) {
        stop(
          "Some neighbor cell names for cell '", cell_name,
          "' were not found in the Seurat object."
        )
      }
      return(idx)
    }
    idx <- as.integer(cell_neighbors)

    if (any(is.na(idx)) || any(idx < 1) || any(idx > ncol(seu))) {
      stop(
        "Some neighbor indices for cell '", cell_name,
        "' are invalid. Indices must be between 1 and the number of cells in the Seurat object."
      )
    }

    return(idx)
  }

  if (neighbors_method == "seurat") {
    cell_idx <- match(cell_name, colnames(seu))
    return(neighbors@nn.idx[cell_idx, ])
  }

  stop("Unsupported neighbors_method in get_neighbors(): ", neighbors_method)
}


#' Calculate entropy for a single cell based on its neighbors
#'
#' @param cell_name Character. Name of the cell to compute entropy for.
#' @param mat Matrix. Expression matrix (genes x cells).
#' @param neighbors_idx Integer vector of neighbor indices into mat.
#'
#' @return A numeric vector of entropy values for each gene.
#' @keywords internal
calculate_entropy_for_cell <- function(cell_name, mat, neighbors_idx) {
  if (length(neighbors_idx) == 0) {
    ent <- rep(0, nrow(mat))
    names(ent) <- rownames(mat)
    return(ent)
  }
  neighbors_expr <- mat[, neighbors_idx, drop = FALSE]
  entropy_gene_expression <- apply(neighbors_expr, 1, estimateur_bub_entropie_discrete)
  return(entropy_gene_expression)
}


#' Calculate binary entropy for a single cell based on its neighbors (fast)
#'
#' @param cell_name Character. Name of the cell to compute entropy for.
#' @param mat_bin Matrix. Binary matrix (features x cells), values in {0,1}.
#' @param neighbors_idx Integer vector of neighbor indices into mat_bin.
#'
#' @return A numeric vector of entropy values for each feature.
#' @keywords internal
calculate_entropy_for_cell_binary_fast <- function(cell_name, mat_bin, neighbors_idx) {
  n <- length(neighbors_idx)
  if (n == 0) {
    ent <- rep(0, nrow(mat_bin))
    names(ent) <- rownames(mat_bin)
    return(ent)
  }
  counts <- rowCounts_CSC_cols(mat_bin, neighbors_idx)
  if (n <= length(.binary_entropy_lookup)) {
    ent <- .binary_entropy_lookup[[n]][counts + 1L]
  } else {
    p1 <- counts / n
    ent <- -(p1 * log2(p1) + (1 - p1) * log2(1 - p1))
    ent[p1 == 0 | p1 == 1] <- 0
  }
  names(ent) <- rownames(mat_bin)
  ent
}



#' Count row occurrences in selected columns of a sparse CSC matrix
#'
#' Efficiently counts how many times each row appears across a subset
#' of columns in a sparse matrix stored in CSC format.
#'
#' @param x Sparse matrix in CSC format.
#' @param cols Integer vector of column indices to inspect.
#'
#' @return Integer vector of counts per row.
#' @keywords internal
rowCounts_CSC_cols <- function(x, cols) {

  counts <- integer(nrow(x))
  p <- x@p
  i <- x@i

  for (j in cols) {
    start <- p[j] + 1L
    end   <- p[j + 1L]

    if (start <= end) {
      idx <- i[start:end] + 1L
      counts[idx] <- counts[idx] + 1L
    }
  }

  counts
}






############# EVERYTHING About Chunks

#' Internal: Prepare cell chunks for chunked entropy computation
#'
#' Splits cells into biologically similar chunks using Leiden clustering.
#' Clusters larger than chunk_cells are split into sub-chunks using k-means
#' on the reduction embeddings within the cluster.
#'
#' @param seu A Seurat object.
#' @param chunk_cells Integer. Maximum number of cells per chunk.
#' @param reduction Character. Name of an existing dimensional reduction in
#'   `seu` to use. Must already exist (resolved by the caller).
#' @param dims Integer vector. Dimensions of `reduction` to use for clustering
#'   and k-means sub-splitting.
#' @param neighbor_assay Character. Assay to set as default for clustering.
#'
#' @return A list where each element is a character vector of cell names
#'   belonging to that chunk.
#' @keywords internal
prepare_chunks <- function(seu, chunk_cells, reduction, dims, neighbor_assay) {
  # work on a copy - never touch the original seu
  seu_chunk <- seu
  Seurat::DefaultAssay(seu_chunk) <- neighbor_assay

  if (!reduction %in% names(seu_chunk@reductions)) {
    stop(sprintf("Reduction '%s' not found in object passed to prepare_chunks().", reduction))
  }

  futile.logger::flog.info(
    "Preparing chunks using Leiden clustering (resolution = 0.3) with %d dims of reduction '%s'.",
    length(dims), reduction
  )

  # FindNeighbors and FindClusters only on seu_chunk - original seu untouched
  seu_chunk <- Seurat::FindNeighbors(seu_chunk, reduction = reduction, dims = dims, verbose = FALSE)
  seu_chunk <- Seurat::FindClusters(seu_chunk, resolution = 0.3, verbose = FALSE)

  cluster_labels <- as.integer(Seurat::Idents(seu_chunk))
  names(cluster_labels) <- colnames(seu_chunk)
  clusters <- split(colnames(seu_chunk), cluster_labels)

  futile.logger::flog.info(
    "Found %d clusters. Splitting any cluster larger than %d cells.",
    length(clusters), chunk_cells
  )

  # split big clusters into sub-chunks
  chunks <- unlist(lapply(clusters, function(cluster_cell_names) {
    if (length(cluster_cell_names) <= chunk_cells) {
      return(list(cluster_cell_names))
    }
    n_subchunks <- ceiling(length(cluster_cell_names) / chunk_cells)
    pca_sub     <- Seurat::Embeddings(seu_chunk[[reduction]])[cluster_cell_names, dims]
    km          <- kmeans(pca_sub, centers = n_subchunks, iter.max = 50, nstart = 3)
    split(cluster_cell_names, km$cluster)
  }), recursive = FALSE)

  futile.logger::flog.info(
    "Prepared %d chunks with median size %d cells.",
    length(chunks), as.integer(median(lengths(chunks)))
  )

  return(chunks)
}



### OLD VERSIONs

### old versions
get_neighbors_old <- function(cell_name, seu, nn_list = NULL, nn.name = NULL) {
  if (!is.null(nn_list)) {
    neighbors <- nn_list[[cell_name]]

    if (is.character(neighbors)) {
      return(match(neighbors, colnames(seu)))
    } else {
      return(neighbors)
    }

  } else {
    assay <- Seurat::DefaultAssay(seu)

    if (is.null(nn.name)) {
      graph_name <- paste0(assay, ".nn")
    } else {
      graph_name <- paste0(nn.name, ".nn")
    }

    if (!graph_name %in% names(seu@neighbors)) {
      stop(sprintf("Graph '%s' not found in Seurat object's neighbors slot.", graph_name))
    }

    cell_idx <- match(cell_name, colnames(seu))
    return(seu@neighbors[[graph_name]]@nn.idx[cell_idx, ])
  }
}


calculate_entropy_for_cell_old <- function(cell_name, mat, seu, nn_list = NULL , nn.name = NULL) {
  neighbors_idx <- get_neighbors_old(cell_name, seu, nn_list, nn.name)
  neighbors_expr <- mat[, neighbors_idx, drop = FALSE]
  entropy_gene_expression <- apply(neighbors_expr, 1, estimateur_bub_entropie_discrete , force = T)
  return(entropy_gene_expression)
}


#' Calculate binary entropy for a single cell based on its neighbors
#'
#' This internal helper function extracts the binary accessibility matrix of a cell's
#' neighbors (0/1) and computes Shannon entropy per feature (Bernoulli p0/p1).
#'
#' @param cell_name Character. Name of the cell to compute entropy for.
#' @param mat_bin Matrix. Binary matrix (features x cells), values in {0,1}.
#' @param seu Seurat object.
#' @param neighbors Neighbor definition. Either a named list or a Seurat Neighbor object.
#' @param neighbors_method Character. One of `"custom"` or `"seurat"`.
#'
#' @return A numeric vector of entropy values for each feature.
#' @keywords internal
calculate_entropy_for_cell_binary <- function(cell_name, mat_bin, seu, neighbors, neighbors_method) {

  neighbors_idx <- get_neighbors(cell_name, seu, neighbors, neighbors_method)

  if (length(neighbors_idx) == 0) {
    ent <- rep(0, nrow(mat_bin))
    names(ent) <- rownames(mat_bin)
    return(ent)
  }

  neighbors_access <- mat_bin[, neighbors_idx, drop = FALSE]

  p1 <- Matrix::rowMeans(neighbors_access)
  p0 <- 1 - p1

  ent <- rep(0, length(p1))
  idx1 <- p1 > 0
  idx0 <- p0 > 0

  ent[idx1] <- ent[idx1] - p1[idx1] * log2(p1[idx1])
  ent[idx0] <- ent[idx0] - p0[idx0] * log2(p0[idx0])

  names(ent) <- rownames(mat_bin)
  return(ent)
}
