#' Estimate gene expression entropy on a Seurat object
#'
#' This function computes gene-level entropy across neighboring cells for a specified assay
#' in a Seurat object. It optionally saves the entropy matrix to file and/or adds it as a
#' new assay named `"{assay}_entropy"`.
#'
#' @param seu A Seurat object.
#' @param assay Character. Name of the assay used for entropy calculation.
#' @param neighbors Optional neighbor definition. If `neighbors_method = "custom"`,
#'   this should be a named list where `neighbors[[cell_name]]` returns a vector of
#'   neighbor cell names or indices. If `neighbors_method = "seurat"`, this should be
#'   a Seurat neighbor object. If `NULL` and `neighbors_method = "auto"`, neighbors
#'   are computed internally.
#' @param neighbors_method Character. How to obtain neighbors. One of `"auto"`,
#'   `"custom"`, or `"seurat"`. `"auto"` computes neighbors inside the function.
#'   `"custom"` uses the user-provided neighbor list in `neighbors`. `"seurat"`
#'   uses the user-provided Seurat neighbor object in `neighbors`.
#' @param neighbor_assay Character. Assay to use for automatic neighbor computation
#'   when `neighbors_method = "auto"`. Defaults to `assay`.
#' @param output_path Optional character string. If provided, saves the full entropy
#'   matrix to this file using `saveRDS()`.
#' @param add_assay Logical. If `TRUE` (default), adds the entropy matrix to the
#'   Seurat object as a new assay called `"{assay}_entropy"`.
#' @param k Integer. Number of nearest neighbors to use when computing neighbors
#'   internally. Default is 15.
#' @param npca Integer. Number of PCA dimensions to use when computing neighbors
#'   internally. Default is 20.
#' @param nfeatures Integer. Number of variable features to select before PCA-based
#'   neighbor computation. Default is 0.
#' @param mode Character. Entropy mode. `"counts"` uses the count-based entropy estimator
#'   (BUB precomputed table). `"binary"` binarizes the assay counts (values > 0) and
#'   computes Shannon entropy for a Bernoulli distribution across neighbors.
#' @return If `add_assay = TRUE`, returns the modified Seurat object with a new assay
#'   named `"{assay}_entropy"`. If `add_assay = FALSE`, returns a data frame of mean
#'   entropy values per cell.
#' @export
#'
#' @examples
#' run_entropy(pbmc_small, assay = "RNA", k = 15, mode = "counts")
run_entropy <- function(
    seu,
    assay,
    neighbors = NULL,
    neighbors_method = c("auto", "custom", "seurat"),
    neighbor_assay = assay,
    output_path = NULL,
    add_assay = TRUE,
    k = 15,
    npca = 20,
    nfeatures = 0,
    mode = c("counts", "binary")
){


  neighbors_method <- match.arg(neighbors_method)
  mode <- match.arg(mode)

  if (!assay %in% Seurat::Assays(seu)) {
    stop(sprintf("Assay '%s' not found in Seurat object.", assay))
  }
  Seurat::DefaultAssay(seu) <- assay
  futile.logger::flog.info("Assay '%s' found and set as default assay.", assay)

  if (!neighbor_assay %in% Seurat::Assays(seu)) {
    stop(sprintf("neighbor_assay '%s' not found in Seurat object.", neighbor_assay))
  }

  futile.logger::flog.info("Starting entropy estimation for assay: %s", assay)
  futile.logger::flog.info(
    "Assay '%s' contains %d cells and %d features",
    assay, ncol(seu), nrow(seu)
  )

  if (neighbors_method == "auto") {

    futile.logger::flog.info(
      "No neighbors provided. Computing neighbors using assay '%s'.",
      neighbor_assay
    )

    seu_neighbor <- seu
    Seurat::DefaultAssay(seu_neighbor) <- neighbor_assay

    reduction_name <- paste0("pca_", neighbor_assay)

    if (!reduction_name %in% names(seu_neighbor@reductions)) {
      futile.logger::flog.info(
        "Reduction '%s' not found. Running PCA on assay '%s'.",
        reduction_name, neighbor_assay
      )

      if (nfeatures > 0) {
        futile.logger::flog.info("Selecting top %s variable features.", nfeatures)
        seu_neighbor <- Seurat::FindVariableFeatures(seu_neighbor, nfeatures = nfeatures)
      }

      seu_neighbor <- Seurat::RunPCA(
        seu_neighbor,
        reduction.name = reduction_name,
        reduction.key = paste0(tolower(reduction_name), "_")
      )
    }

    num_pcs <- ncol(Seurat::Embeddings(seu_neighbor[[reduction_name]]))
    npca <- min(npca, num_pcs)

    if (npca < 1) {
      stop("No PCA dimensions available for neighbor computation.")
    }

    neighbor_name <- paste0("entropy_", neighbor_assay, "_nn")

    seu_neighbor <- Seurat::FindNeighbors(
      seu_neighbor,
      reduction = reduction_name,
      dims = 1:npca,
      return.neighbor = TRUE,
      k.param = k,
      graph.name = neighbor_name
    )

    if (!neighbor_name %in% names(seu_neighbor@neighbors)) {
      stop(
        "Computed Neighbor object '", neighbor_name,
        "' was not found in the Seurat object."
      )
    }

    neighbors <- seu_neighbor@neighbors[[neighbor_name]]
    neighbors_method <- "seurat"

  } else if (neighbors_method == "custom") {
    if (is.null(neighbors)) {
      stop("`neighbors` must be provided when neighbors_method = 'custom'.")
    }
    if (!is.list(neighbors) || is.null(names(neighbors))) {
      stop(
        "`neighbors` must be a named list where neighbors[[cell_name]] ",
        "returns a vector of neighbor cell names or indices."
      )
    }

    cell_names <- colnames(seu)

    missing_cells <- setdiff(cell_names, names(neighbors))

    if (length(missing_cells) > 0) {
      stop(
        "The custom neighbors list is missing entries for ",
        length(missing_cells), " cells.\n",
        "Each cell must appear in names(neighbors), even if it has no neighbors.\n",
        "Example structure:\n",
        "neighbors[['cell1']] -> c('cell2','cell3',...)"
      )
    }
    futile.logger::flog.info("Using user-provided custom neighbors.")

  } else if (neighbors_method == "seurat") {
    if (is.null(neighbors)) {
      stop("`neighbors` must be provided when neighbors_method = 'seurat'.")
    }

    if (!inherits(neighbors, "Neighbor")) {
      stop(
        "For neighbors_method = 'seurat', `neighbors` must be a Seurat Neighbor object ",
        "(e.g. the result of Seurat::FindNeighbors(..., return.neighbor = TRUE))."
      )
    }

    if (!all(colnames(seu) %in% neighbors@cell.names)) {
      stop("Cell names in Seurat object do not match the provided Seurat neighbor object.")
    }

    futile.logger::flog.info("Using user-provided Seurat neighbor object.")
  }


  ######################


  mat <- Seurat::GetAssayData(seu, assay = assay, layer = "counts")

  if (mode == "binary") {
    futile.logger::flog.info("Binarizing matrix (counts > 0) for binary entropy.")
    mat_bin <- Matrix::Matrix(mat > 0, sparse = TRUE)
  }

  # Get all cell names
  cell_names <- colnames(seu)
  n_cells <- length(cell_names)

  # Initialize the progress bar
  pb <- progress::progress_bar$new(
    format = "  Processing [:bar] :percent | :current/:total cells",
    total = n_cells,
    clear = FALSE,
    width = 60
  )

  #.data_loaded <<- FALSE


  if (mode == "counts") {
    futile.logger::flog.info("Reading the precomputed table (counts mode)")
    .precomputed_bub <<- readRDS(system.file(
      "extdata", "precomputed_bub_5000.rds",
      package = "Entropy"
    ))
  } else if (mode == "binary") {
    futile.logger::flog.info("Reading the binary entropy lookup (binary mode)")
    .binary_entropy_lookup <<- readRDS(system.file(
      "extdata", "binary_entropy_lookup_100.rds",
      package = "Entropy"
    ))
  } else {
    stop("Unknown mode: ", mode, ". Expected 'counts' or 'binary'.")
  }


  # Apply with progress
  entropy_matrix_list <- lapply(cell_names, function(cell_name) {
    pb$tick()

    if (mode == "binary") {
      calculate_entropy_for_cell_binary_fast(
        cell_name = cell_name,
        mat_bin   = mat_bin,
        seu       = seu,
        neighbors = neighbors,
        neighbors_method = neighbors_method
      )
    } else {
      calculate_entropy_for_cell(
        cell_name = cell_name,
        mat       = mat,
        seu       = seu,
        neighbors = neighbors,
        neighbors_method = neighbors_method
      )
    }
  })

  names(entropy_matrix_list) <- cell_names

  entropy_matrix <- do.call(cbind, entropy_matrix_list)


  if (!is.null(output_path)) {
    saveRDS(entropy_matrix, output_path)
    futile.logger::flog.info("Entropy matrix saved to: %s", output_path)
  }


  if (add_assay) {
    entropy_assay <- Seurat::CreateAssayObject(counts = entropy_matrix)
    assay_name <- paste0(assay, "_entropy")
    seu[[assay_name]] <- entropy_assay
    futile.logger::flog.info("Entropy matrix added to Seurat object as '%s'.", assay_name)
    return(seu)
  } else {
    means <- colMeans(entropy_matrix)
    Ent_mean <- as.data.frame(means)
    return(Ent_mean)
  }

}









