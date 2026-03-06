



library(Rcpp)
sourceCpp("~/Desktop/EMBL/R_projects/Plasticity_Guard/Plasticity/src/Entropy/Rccp_bub_function.cpp")


precomputed_bub <- vector("list", 40)  # for N
for (N in 1:40) {
  precomputed_bub[[N]] <- vector("list", 5001)  # for m = 0:2000
  for (m in 0:5000) {
    precomputed_bub[[N]][[m + 1]] <- vector("list", 12)  # for k = 1:12
    for (k in 1:12) {
      if (k <= N) {
        precomputed_bub[[N]][[m + 1]][[k]] <- tryCatch(
          bub_opti_cpp(N, m, k),
          error = function(e) NA
        )
      }
    }
  }
}


saveRDS(precomputed_bub, "inst/extdata/precomputed_bub_5000.rds")






build_binary_entropy_lookup <- function(max_k = 30) {
  lookup <- vector("list", max_k)
  for (k in seq_len(max_k)) {
    p <- (0:k) / k
    H <- -(p * log2(p) + (1 - p) * log2(1 - p))
    H[!is.finite(H)] <- 0
    lookup[[k]] <- H
  }
  lookup
}

binary_entropy_lookup <- build_binary_entropy_lookup(100)
saveRDS(binary_entropy_lookup,
          "inst/extdata/binary_entropy_lookup_100.rds")


