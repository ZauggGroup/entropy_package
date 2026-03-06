args = commandArgs(TRUE)
n_parts =  as.integer(args[1]) # first argument of the bash script
part_number = as.integer(args[2]) # second argument of the bash script


n_parts
part_number

#library(renv)
#renv::activate(project = '/g/scb/zaugg/kurti/PBMCs/')

library(cli, lib.loc = '/g/scb/zaugg/kurti/PBMCs/renv/library/R-4.2/x86_64-pc-linux-gnu')
library(Seurat, lib.loc = '/g/scb/zaugg/kurti/PBMCs/renv/library/R-4.2/x86_64-pc-linux-gnu')
library(Signac, lib.loc = '/g/scb/zaugg/kurti/PBMCs/renv/library/R-4.2/x86_64-pc-linux-gnu')
library(qs, lib.loc = '/g/scb/zaugg/kurti/PBMCs/renv/library/R-4.2/x86_64-pc-linux-gnu')
library(Matrix, lib.loc = '/g/scb/zaugg/kurti/PBMCs/renv/library/R-4.2/x86_64-pc-linux-gnu')


# Set path to save the results
data_processed <- "/g/scb/zaugg/kurti/Hematopoiesis/Y1_T1_from_paper/ATAC/"

# Load data
seu <- qread('/g/scb/zaugg/kurti/Hematopoiesis/Y1_T1_from_paper/Y1_paper_set.qs')
seu <- FindNeighbors(seu, reduction ="pca", dims = 1:15, return.neighbor = T, k.param = 15)

DefaultAssay(seu) <- "ATAC"
#seu <- FindVariableFeatures(seu, nfeatures = 2000)

mat <- seu@assays$ATAC@counts
bin_mat <-  Matrix(apply(mat, 2, function(x) ifelse(x != 0, 1, 0)), sparse=T)
#bin_mat <- as.data.frame(mat) %>%
# mutate_all(~ ifelse(. != 0, 1, 0))


fun_ent = function (p) {
  p = p[p>0] #we keep only probability > 0 because log2(0) is not defined
  Ent = -sum(p*log2(p))
  return(Ent)
}

# function to compute the probability of having 0 and 1
cal_prob <- function(input_vector) {
  # Count the occurrences of 0 and 1
  count_0 <- sum(input_vector == 0)
  count_1 <- sum(input_vector == 1)

  # Calculate the probabilities
  probability_0 <- count_0 / length(input_vector)
  probability_1 <- count_1 / length(input_vector)

  result <- c(probability_0,
              probability_1)
  return(fun_ent(result))
}

library(parallel)
# separate cell Id in x parts
cell_idx_list <- splitIndices(ncol(seu), n_parts)

## separately for each part
calculate_entropy_for_cell <- function(cell_idx, bin_mat, seu) {
  neighbors_idx <- seu@neighbors$SCT.nn@nn.idx[cell_idx,] #or RNA.nn depending if you have or not multiome
  neighbors_access <- bin_mat[, neighbors_idx]
  entropy_access <- apply(neighbors_access, 1, cal_prob)
  return(entropy_access)
}

# Use mclapply to parallelize the loop
entropy_matrix_list <- mclapply(cell_idx_list[[part_number]], function(cell_idx) {
  calculate_entropy_for_cell(cell_idx, bin_mat, seu)
}, mc.cores = 16)

names(entropy_matrix_list) <- mat@Dimnames[[2]][cell_idx_list[[part_number]]]

# Convert the list of results into a matrix
entropy_matrix <- do.call(cbind, entropy_matrix_list)

qsave(entropy_matrix, paste0(data_processed,"entropy_matrix_all_15dim_15K_y1_best_" ,part_number,".qs"))
