#' Internal: BUB entropy estimator in pure R
#'
#' This is a full R implementation of the BUB entropy estimator adapted from Mathilde Gaillard’s version.
#' It avoids the need for C++ and can be used either for runtime fallback or to precompute values.
#'
#' @param N Integer. Number of total values (length of vector).
#' @param m Integer. Number of bins (unique values).
#' @param k_max Integer. Number of coefficients to use.
#'
#' @return A numeric vector of entropy coefficients.
#' @keywords internal
bub_opti <- function(N, m, k_max) {
  p <- (0:N) / N
  B <- mat_bernouilli_polynom(p, N)
  g <- sapply(p, f, m = m)
  g_col <- matrix(g, ncol = 1)
  h <- c(0, sapply(p[-1], H))
  Y_tot <- matrix(h, ncol = 1)

  a <- (0:N) / N
  a <- -a * log(a) + ((1 - a) / (2 * N))

  best_MM <- Inf

  for (i in 1:k_max) {
    h_exp <- a[(i + 1):length(a)] %*% B[(i + 1):nrow(B), ]
    Y_non_pond <- Y_tot - t(h_exp)

    G <- pracma::repmat(g, i, 1)
    X <- t(G * B[1:i, ])

    D <- pracma::Diag(rep(-1, i), 0) + pracma::Diag(rep(1, i - 1), 1)
    U <- crossprod(X) + (N / 4) * (crossprod(D))
    U[i, i] <- U[i, i] + (N / 4)

    Y <- g_col * Y_non_pond
    XY <- crossprod(X, Y)
    XY[i] <- XY[i] + (N / 4) * a[i + 1]

    a[1:i] <- corpcor::pseudoinverse(U) %*% XY

    biais <- 2 * g * abs((t(Y_tot) - a %*% B))
    maxbiais <- max(abs(biais))
    borne_var <- max((a - c(a[2:length(a)], 0))^2)
    MM <- sqrt(maxbiais^2 + N * borne_var)

    if (MM < best_MM) {
      best_MM <- MM
      best_a <- a
    }
  }

  return(best_a)
}

#' @keywords internal
mat_bernouilli_polynom <- function(p, N) {
  fa <- lgamma(1:(N + 1))
  Ni <- fa[N + 1] - fa[1:(N + 1)] - rev(fa)

  B <- pracma::zeros(n = N + 1, m = N + 1)

  p <- p[-1]
  p <- p[-length(p)]

  lp <- log(p)
  lq <- rev(lp)

  for (i in 0:N) {
    B[i + 1, 2:N] <- Ni[i + 1] + i * lp + (N - i) * lq
  }

  B <- exp(B)
  B[2:(N + 1), 1] <- 0
  B[1:N, N + 1] <- 0

  return(B)
}

#' @keywords internal
H <- function(x) {
  nz <- x > 0
  return(-sum(x[nz] * log(x[nz])))
}

#' @keywords internal
f <- function(x, m) {
  if (x < (1 / m)) {
    return(m)
  } else {
    return(1 / x)
  }
}
