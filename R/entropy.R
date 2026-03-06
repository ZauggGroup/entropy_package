#' Estimate gene expression entropy on a Seurat object
#'
#' This function computes gene-level entropy across neighboring cells for a specified assay
#' in a Seurat object. It optionally saves the entropy matrix to file and/or adds it as a
#' new assay named `"assay_entropy"`.
#'
#' @param seu A Seurat object.
#' @param assay Character. Name of the assay to use.
#' @param nn_list Optional named list of neighbor cell indices or names per cell. If NULL, neighbors are computed using PCA.
#' @param output_path Optional character string. If provided, saves the full entropy matrix to this file (using `qs::qsave()`).
#' @param add_assay Logical. If TRUE (default), adds the entropy matrix to the Seurat object as a new assay called `"{assay}_entropy"`.
#' @param k.param Integer. Number of nearest neighbors to use when computing the neighbor graph. Default is 15.
#' @param npca Integer. Number of PCA dimensions to use when computing neighbors. Default is 20.
#' @param nfeatures Integer. Number of variable features to select before PCA and neighbor finding. Default is 0.
#' @param nn.name Optional character. If neighbors are already present in the Seurat object, this specifies the name of the neighbor graph to use (e.g. `"SCT"` or `"RNA"`).
#' @param mode Character. Entropy mode. `"counts"` uses the count-based entropy estimator
#'   (BUB precomputed table). `"binary"` binarizes the assay counts (values > 0) and
#'   computes Shannon entropy for a Bernoulli distribution (p0/p1) across neighbors.
#'   Default is `"counts"`.
#' @return If `add_assay = TRUE`, returns the modified Seurat object with a new assay
#'         named `"{assay}_entropy"`. If `add_assay = FALSE`, returns a data frame
#'         of mean entropy values per cell.
#' @export
#'
#' @examples
#' # Assuming `pbmc` is a Seurat object with an "RNA" assay:
#' run_entropy(pbmc_small, assay = "RNA", k.param = 15, mode = "counts")
#' run_entropy(pbmc_small, assay = "ATAC", k.param = 15, mode = "binary")
run_entropy <- function(seu, assay, nn_list = NULL, nn.name = NULL, output_path = NULL, add_assay = TRUE, k.param = 15, npca = 20, nfeatures = 0, mode = c("counts", "binary")) {


  futile.logger::flog.info("Starting entropy estimation for assay: %s", assay)
  futile.logger::flog.info("Input Seurat object contains %d cells and %d features", ncol(seu), nrow(seu))


  if (assay %in% Seurat::Assays(seu)) {
    Seurat::DefaultAssay(seu) <- assay
    futile.logger::flog.info("Assay '%s' found and set as default assay.", assay)
  } else {
    stop(sprintf("Assay '%s' not found in Seurat object.", assay))
  }

  if( is.null(nn_list) ){
    futile.logger::flog.info("No neighbor list provided. Computing neighbors using PCA reduction.")
    num_pcs <- ncol(Seurat::Embeddings(seu[["pca"]]))
    npca = min(npca,num_pcs)
    seu <- Seurat::FindNeighbors(seu, reduction = "pca", dims = 1:npca, return.neighbor = TRUE, k.param = k.param)
    if( nfeatures > 0 ){
      futile.logger::flog.info("Selecting top %s variable features." , nfeatures)
      seu <- Seurat::FindVariableFeatures(seu, nfeatures = nfeatures)
    }

  }


  mat <- Seurat::GetAssayData(seu, layer = "counts")

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
    .precomputed_bub <<- qs::qread(system.file(
      "extdata", "precomputed_bub_5000.qs",
      package = "Entropy"
    ))
  } else if (mode == "binary") {
    futile.logger::flog.info("Reading the binary entropy lookup (binary mode)")
    .binary_entropy_lookup <<- qs::qread(system.file(
      "extdata", "binary_entropy_lookup_100.qs",
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
        nn_list   = nn_list,
        nn.name   = nn.name
      )
    } else {
      calculate_entropy_for_cell(
        cell_name = cell_name,
        mat       = mat,
        seu       = seu,
        nn_list   = nn_list,
        nn.name   = nn.name
      )
    }
  })

  names(entropy_matrix_list) <- cell_names

  entropy_matrix <- do.call(cbind, entropy_matrix_list)


  if (!is.null(output_path)) {
    qsave(entropy_matrix, output_path)
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









