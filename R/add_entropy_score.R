#' Add mean entropy scores for gene sets to a Seurat object
#'
#' Computes the mean entropy of one or more gene sets from a precomputed entropy
#' assay (created by `run_entropy(..., add_assay = TRUE)`) and adds one metadata
#' column per set, named `entropy_<set name>`. This does not recompute entropy;
#' it averages already-stored per-gene entropy values, so it is fast.
#'
#' @param seu A Seurat object containing a precomputed entropy assay.
#' @param entropy_assay Character. Name of the entropy assay to read from
#'   (e.g. `"RNA_entropy"`). Must already exist; compute it first with
#'   `run_entropy(..., add_assay = TRUE)`.
#' @param features_list A named list of character vectors. Each element is a set
#'   of gene/feature names, and its name becomes the metadata column suffix
#'   (e.g. list element `pathwayA` -> column `entropy_pathwayA`).
#' @param force_drop Logical. If `FALSE` (default), any requested feature not
#'   present in the entropy assay raises an error. If `TRUE`, missing features
#'   are dropped with a warning. A set with no present features always errors.
#'
#' @return The Seurat object with one new metadata column per gene set.
#' @export
add_entropy_score <- function(seu, entropy_assay, features_list, force_drop = FALSE) {

  if (!entropy_assay %in% Seurat::Assays(seu)) {
    stop(sprintf(
      "Entropy assay '%s' not found. Compute it first with run_entropy(..., add_assay = TRUE).",
      entropy_assay
    ))
  }

  if (!is.list(features_list) || is.null(names(features_list)) || any(names(features_list) == "")) {
    stop("`features_list` must be a named list of character vectors, e.g. list(pathwayA = c('gene1','gene2')).")
  }

  ent_mat <- SeuratObject::GetAssayData(seu, assay = entropy_assay, layer = "counts")
  available <- rownames(ent_mat)

  for (set_name in names(features_list)) {
    genes <- unique(features_list[[set_name]])

    if (!is.character(genes)) {
      stop(sprintf("Set '%s' must be a character vector of feature names.", set_name))
    }

    missing <- setdiff(genes, available)
    present <- intersect(genes, available)

    if (length(present) == 0) {
      stop(sprintf(
        "None of the %d features in set '%s' are present in assay '%s'. Compute entropy for them first.",
        length(genes), set_name, entropy_assay
      ))
    }

    if (length(missing) > 0) {
      if (!force_drop) {
        stop(sprintf(
          "%d of %d features in set '%s' are not in assay '%s': %s%s. Compute entropy for them first, or set force_drop = TRUE.",
          length(missing), length(genes), set_name, entropy_assay,
          paste(utils::head(missing, 10), collapse = ", "),
          if (length(missing) > 10) ", ..." else ""
        ))
      }
      futile.logger::flog.warn(
        "Set '%s': dropping %d missing features, scoring on %d.",
        set_name, length(missing), length(present)
      )
    }

    col <- paste0("entropy_", set_name)
    seu[[col]] <- Matrix::colMeans(ent_mat[present, , drop = FALSE])
    futile.logger::flog.info("Added metadata column '%s' (mean entropy of %d features).", col, length(present))
  }

  seu
}
