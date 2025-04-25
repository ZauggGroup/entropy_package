#' Get neighbor indices for a cell from either a neighbor list or Seurat object
#'
#' @param cell_name Cell name (character)
#' @param seu Seurat object
#' @param nn_list Optional list of neighbors by cell
#' @param nn.name Optional character. If neighbors are already present in the Seurat object, this specifies the name of the neighbor graph to use (e.g. `"SCT"` or `"RNA"`).
#' @return Integer vector of neighbor indices
#' @keywords internal
get_neighbors <- function(cell_name, seu, nn_list = NULL, nn.name = NULL) {
  if (!is.null(nn_list)) {
    neighbors <- nn_list[[cell_name]]
    if (is.character(neighbors)) {
      return(match(neighbors, colnames(seu)))
    } else {
      return(neighbors)
    }
  } else {
    # Try to infer the graph name from the default assay
    assay <- Seurat::DefaultAssay(seu)
    if(is.null(nn.name)) {
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






#' Calculate entropy for a single cell based on its neighbors
#'
#' This internal helper function extracts the expression matrix of a cell's neighbors
#' and computes entropy for each gene using the specified estimator.
#'
#' @param cell_name Character. Name of the cell to compute entropy for.
#' @param mat Matrix. Expression matrix (genes x cells).
#' @param seu Seurat object.
#' @param nn_list Optional list mapping cell names to neighbor indices or names.
#' @param nn.name Optional character. If neighbors are already present in the Seurat object, this specifies the name of the neighbor graph to use (e.g. `"SCT"` or `"RNA"`).
#'
#' @return A numeric vector of entropy values for each gene.
#' @keywords internal
calculate_entropy_for_cell <- function(cell_name, mat, seu, nn_list = NULL , nn.name = NULL) {
  neighbors_idx <- get_neighbors(cell_name, seu, nn_list, nn.name)
  neighbors_expr <- mat[, neighbors_idx, drop = FALSE]
  entropy_gene_expression <- apply(neighbors_expr, 1, estimateur_bub_entropie_discrete)
  return(entropy_gene_expression)
}












