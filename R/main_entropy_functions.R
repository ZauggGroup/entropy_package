#' Internal: Cached wrapper for the BUB discrete entropy estimator
#'
#' This function computes discrete entropy using the BUB estimator with caching.
#' It avoids re-running the C++ estimator for repeated values of N and m.
#'
#' @param X Integer vector. The discrete data to estimate entropy for.
#'
#' @return A single numeric entropy value.
#' @keywords internal
estimateur_bub_entropie_discrete <- function(X) {
  if (sum(X) == 0) return(0)

  N <- length(X)
  h <- tabulate(X + 1)
  m <- length(h)

  hist_eff <- numeric(N + 1)
  for (i in 1:(N + 1)) {
    hist_eff[i] <- sum(h == (i - 1))
  }

  A <- get_bub_precomputed(N, m, k_max = 11)
  H_BUB <- sum(A * hist_eff)

  return(H_BUB)
}



#' Internal: Get precomputed or fallback BUB entropy coefficients (R-only)
#'
#' @param N Integer. Total number of values.
#' @param m Integer. Number of unique values.
#' @param k_max Integer. Max number of basis terms.
#'
#' @return A numeric vector of coefficients from the BUB estimator.
#' @keywords internal
get_bub_precomputed <- function(N, m, k_max) {
  k_eff <- min(N, k_max)

  # Load precomputed table if needed
  if (!exists(".precomputed_bub", envir = .GlobalEnv)) {
    .precomputed_bub <<- qs::qread(system.file("extdata", "precomputed_bub.qs", package = "Entropy"))
  }

  # If within precomputed range
  if (N >= 1 && N <= 40 && m >= 0 && m <= 4000 && k_eff >= 1 && k_eff <= 12) {
    result <- .precomputed_bub[[N]][[m + 1]][[k_eff]]
    if (!is.null(result)) {
      return(result)
    }
  }

  # Fallback to pure R computation
  futile.logger::flog.info(
    "Using on-the-fly R computation for N = %d, m = %d, k = %d", N, m, k_eff
  )
  return(bub_opti(N, m, k_eff))
}



