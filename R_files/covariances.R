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