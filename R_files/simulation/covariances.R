get_cor_mat <- function(kernel_type, distsqrd_mat, rate) {
  if (kernel_type == "rbf") {
    return(exp(-(rate^2 / 2) * distsqrd_mat))
  } else if (kernel_type == "matern_5_2") {
    dist_mat <- sqrt(distsqrd_mat)
    return(
      (1 + rate * sqrt(5) * dist_mat + (rate^2 * 5 / 3) * distsqrd_mat) *
        exp(-rate * sqrt(5) * dist_mat)
    )
  } else {
    stop("Invalid kernel type")
  }
}

fgn_cov <- function(t, hurst) {
  d <- t - t[1] # compute differences in time
  q <- 0.5 * (abs(d - 1)^(2 * hurst) + abs(d + 1)^(2 * hurst) - 2 * abs(d)^(2 * hurst))
  Q <- toeplitz(q)
  return(Q)
}

ar2_cov <- function(t, xi) {
  if (xi[1] == xi[2] || any(abs(xi) < 1)) {
    stop("Elements of xi must be distinct and larger than 1 in absolute value")
  }
  d <- t - t[1] # differences in time
  g <- xi[1]^2 * xi[2]^2 / ((xi[1] * xi[2] - 1) * (xi[2] - xi[1])) *
    (xi[1]^(1 - d) / (xi[1]^2 - 1) - xi[2]^(1 - d) / (xi[2]^2 - 1))
  Q <- toeplitz(g / g[0])
  return(Q)
}

# Non-separable space-time covariance kernel (Gneiting 2002 class, R^3)
# C(h; u | a, b, c) =
#   c / (a^2 u^2 + 1)^(3/2) * exp(-b^2 ||h||^2 / (a^2 u^2 + 1))
# a >= 0: non-separability; a = 0 gives a separable spatial RBF kernel
# b >= 0: spatial scale (analogous to phi_gamma)
# c >= 0: overall scale (typically 1; apply k_gamma externally)
# spatial_distsqrd: L x L matrix of squared spatial distances
# time_lags_sqrd:   M x M matrix of squared temporal lags
# Returns the LM x LM joint covariance matrix, ordering (l-1)*M + m
nonsep_cov <- function(spatial_distsqrd, time_lags_sqrd, a, b, c = 1) {
  L <- nrow(spatial_distsqrd)
  M <- nrow(time_lags_sqrd)

  denom <- a^2 * time_lags_sqrd + 1

  result <- matrix(0.0, nrow = L * M, ncol = L * M)
  row_base <- (seq_len(L) - 1L) * M
  for (m1 in seq_len(M)) {
    for (m2 in seq_len(M)) {
      d <- denom[m1, m2]
      block <- (c / d^1.5) * exp(-b^2 * spatial_distsqrd / d)
      result[row_base + m1, row_base + m2] <- block
    }
  }

  result
}