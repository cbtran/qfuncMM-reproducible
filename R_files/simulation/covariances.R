get_cor_mat <- function(kernel_type, distsqrd_mat, rate) {
  if (kernel_type == "rbf") {
    return(exp(-(rate^2 / 2) * distsqrd_mat))
  } else if (kernel_type == "matern_1_2") {
    return(exp(-rate * sqrt(distsqrd_mat)))
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

  # time_lags_sqrd is symmetric Toeplitz: [m1,m2] = (m1-m2)^2.
  # Only M unique denominators exist (one per lag k = 0,...,M-1).
  # Precompute M spatial blocks; the placement loop then only copies — no exp.
  unique_denom <- a^2 * time_lags_sqrd[1, ] + 1  # denom for lag 0, 1, ..., M-1
  blocks <- lapply(unique_denom, function(d) {
    (c / d^1.5) * exp(-b^2 * spatial_distsqrd / d)
  })

  result <- matrix(0.0, nrow = L * M, ncol = L * M)
  row_base <- (seq_len(L) - 1L) * M
  for (m1 in seq_len(M)) {
    for (m2 in seq_len(M)) {
      result[row_base + m1, row_base + m2] <- blocks[[abs(m1 - m2) + 1]]
    }
  }
  result
}

# Nonseparable Matern space-time covariance (Ip & Li 2017, eq. 1.5, nu=9/2, d=3)
# Bessel order: nu - (d+1)/2 = 9/2 - 2 = 5/2.
# Closed form (K_{5/2}(r) = sqrt(pi/2r)*e^{-r}*(1+3/r+3/r^2), Gamma(5/2)=3*sqrt(pi)/4):
#   C(h; u | alpha, beta) = sigma2 * (1 + r + r^2/3) * exp(-r)
#   where r = sqrt(alpha^2 * ||h||^2 + beta^2 * u^2)
# Spatial marginal at u=0: sigma2 * (1 + alpha*h + (alpha*h)^2/3) * exp(-alpha*h),
#   which is Matern-5/2 with rate alpha/sqrt(5). Psi is calibrated via this marginal.
# alpha >= 0: spatial rate (psi calibrated using matern_5_2 at rate alpha/sqrt(5))
# beta  >= 0: temporal rate
# sigma2 >= 0: overall scale (typically 1; apply k_gamma externally)
# spatial_distsqrd: L x L matrix of squared spatial distances
# time_lags_sqrd:   M x M matrix of squared temporal lags
# Returns the LM x LM joint covariance matrix, ordering (l-1)*M + m
nonsep_matern_cov <- function(spatial_distsqrd, time_lags_sqrd, alpha, beta,
                              sigma2 = 1) {
  L <- nrow(spatial_distsqrd)
  M <- nrow(time_lags_sqrd)

  # r[l1,l2,m1,m2] = sqrt(alpha^2 * S[l1,l2] + beta^2 * T[m1,m2])
  # val = sigma2 * (1 + r + r^2/3) * exp(-r)
  # Build L^2 x M^2 via outer sum, then reshape and permute to (L*M) x (L*M).
  s_vec <- as.vector(spatial_distsqrd)  # L^2, column-major: l1 + L*(l2-1)
  t_vec <- as.vector(time_lags_sqrd)   # M^2, column-major: m1 + M*(m2-1)

  # r_sqrd_flat[i_s, i_t] = alpha^2 * S[i_s] + beta^2 * T[i_t]
  r_sqrd_flat <- outer(alpha^2 * s_vec, rep(1, M^2)) +
    outer(rep(1, L^2), beta^2 * t_vec)               # L^2 x M^2
  r_flat <- sqrt(r_sqrd_flat)
  val_flat <- sigma2 * (1 + r_flat + r_flat^2 / 3) * exp(-r_flat)

  # [l1,l2,m1,m2] -> permute to [m1,l1,m2,l2] -> (L*M) x (L*M)
  arr <- aperm(array(val_flat, dim = c(L, L, M, M)), c(3, 1, 4, 2))
  matrix(arr, nrow = L * M, ncol = L * M)
}