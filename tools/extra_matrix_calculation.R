



library(Rcpp)
sourceCpp("src/Rccp_bub_function.cpp")


precomputed_bub <- vector("list", 40)  # for N
for (N in 1:40) {
  precomputed_bub[[N]] <- vector("list", 2001)  # for m = 0:2000
  for (m in 0:2000) {
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


qs::qsave(precomputed_bub, "inst/extdata/precomputed_bub.qs")









