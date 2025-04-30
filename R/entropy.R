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
#'
#' @return If `add_assay = TRUE`, returns the modified Seurat object with a new assay
#'         named `"{assay}_entropy"`. If `add_assay = FALSE`, returns a data frame
#'         of mean entropy values per cell.
#' @export
#'
#' @examples
#' # Assuming `pbmc` is a Seurat object with an "RNA" assay:
#' run_entropy(pbmc_small, assay = "RNA", k.param = 15)
run_entropy <- function(seu, assay, nn_list = NULL, nn.name = NULL, output_path = NULL, add_assay = TRUE, k.param = 15, npca = 20, nfeatures = 0) {


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


  futile.logger::flog.info("Reading the precomputed table")
  .precomputed_bub <<- qs::qread(system.file("extdata", "precomputed_bub_5000.qs", package = "Entropy"))


  # Apply with progress
  entropy_matrix_list <- lapply(cell_names, function(cell_name) {
    pb$tick()
    calculate_entropy_for_cell(cell_name, mat, seu, nn_list = nn_list , nn.name = nn.name)
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









