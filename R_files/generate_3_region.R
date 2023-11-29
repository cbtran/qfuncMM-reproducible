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

#' Simulate signals from 3 regions
#'
#' @param num_sim number of simulations
#' @param voxel_coords List of region coordinates. Each row is a voxel and
#'   the three columns are x, y, and z coordinates.
#' @param n_timept number of timepoints
#' @param true_corr vector of inter-regional correlations (r12, r13, r23)
#' @param shared_params vector of parameters shared across regions (tau_eta, k_eta, nugget)
#' @param region_params 3 x 5 dataframe of region-specific parameters.
#'    rows (region 1, region 2, region 3)
#'    columns (phi_gamma, tau_gamma, k_gamma, nugget_gamma, mean)
#' @param sigma2 noise variance
#' @param c_kernel_type Choice of spatial kernel. Defaul "matern_5_2".
#' @return obs_signal Simulated signal
generate_3_region <- function(
  num_sim, voxel_coords, n_timept, true_corr, region_params, shared_params,
  sigma2 = 1, seed = 1, c_kernel_type = "matern_5_2") {

  stopifnot(n_timept >= 1)
  stopifnot(length(true_corr) == 3)
  stopifnot(length(shared_params) == 3)
  stopifnot(nrow(region_params) == 3 && ncol(region_params) == 5)
  stopifnot(sigma2 >= 0)
  stopifnot(c_kernel_type %in% c("matern_5_2", "rbf"))

  n_voxel1 <- nrow(voxel_coords[[1]])
  n_voxel2 <- nrow(voxel_coords[[2]])
  n_voxel3 <- nrow(voxel_coords[[3]])
  # Get distance and time matrices
  dist_sqrd_mat_region1 <- as.matrix(dist(voxel_coords[[1]]))^2
  dist_sqrd_mat_region2 <- as.matrix(dist(voxel_coords[[2]]))^2
  dist_sqrd_mat_region3 <- as.matrix(dist(voxel_coords[[3]]))^2
  timesqrd_mat <- (outer(1:n_timept, 1:n_timept, `-`))^2

  # Create covariance matrices

  corr_mat <- matrix(nrow = 3, ncol = 3)
  corr_mat[upper.tri(corr_mat)] <- true_corr
  corr_mat[lower.tri(corr_mat)] <- true_corr
  diag(corr_mat) <- 1

  # Covariance of eta effect
  A <- shared_params["k_eta"] *
    get_cor_mat("rbf", timesqrd_mat, shared_params["tau_eta"]) +
      shared_params["nugget"] * diag(n_timept)
  eta_sigma <- kronecker(corr_mat, A)

  # Matrices of gamma effects. C is spatial correlation matrix, B is temporal covariance matrix
  phi_gamma <- region_params$phi_gamma
  tau_gamma <- region_params$tau_gamma
  k_gamma <- region_params$k_gamma
  nugget_gamma <- region_params$nugget_gamma

  ## Region 1
  C1 <- get_cor_mat(c_kernel_type, dist_sqrd_mat_region1, phi_gamma[1])
  B1 <- k_gamma[1] * get_cor_mat("rbf", timesqrd_mat, tau_gamma[1]) +
    nugget_gamma[1] * diag(n_timept)
  gamma_sigma1 <- kronecker(C1, B1)

  ## Region 2
  C2 <- get_cor_mat(c_kernel_type, dist_sqrd_mat_region2, phi_gamma[2])
  B2 <- k_gamma[2] * get_cor_mat("rbf", timesqrd_mat, tau_gamma[2]) +
    nugget_gamma[2] * diag(n_timept)
  gamma_sigma2 <- kronecker(C2, B2)

  ## Region 3
  C3 <- get_cor_mat(c_kernel_type, dist_sqrd_mat_region3, phi_gamma[3])
  B3 <- k_gamma[3] * get_cor_mat("rbf", timesqrd_mat, tau_gamma[3]) +
    nugget_gamma[3] * diag(n_timept)
  gamma_sigma3 <- kronecker(C3, B3)

  eta <- MASS::mvrnorm(
          num_sim, mu = rep(0, 3 * n_timept), Sigma = sigma2 * eta_sigma) |>
          matrix(nrow = num_sim)
  gamma_r1 <- MASS::mvrnorm(
               num_sim, mu = rep(0, n_voxel1 * n_timept), Sigma = sigma2 * gamma_sigma1) |>
               matrix(nrow = num_sim)
  gamma_r2 <- MASS::mvrnorm(
               num_sim, mu = rep(0, n_voxel2 * n_timept), Sigma = sigma2 * gamma_sigma2) |>
               matrix(nrow = num_sim)
  gamma_r3 <- MASS::mvrnorm(
               num_sim, mu = rep(0, n_voxel3 * n_timept), Sigma = sigma2 * gamma_sigma3) |>
               matrix(nrow = num_sim)


  mu_1 <- region_params$mean[1]
  mu_2 <- region_params$mean[2]
  mu_3 <- region_params$mean[3]
  n_voxel_total <- n_voxel1 + n_voxel2 + n_voxel3
  obs_signal <- vector(mode = "list", length = num_sim)
  for (i in seq_along(obs_signal)) {
      # Generate iid noise
      set.seed(seed + 1000 + i)
      noise <- diag(rep(sqrt(sigma2), n_voxel_total * n_timept)) %*% rnorm(n_voxel_total * n_timept)

      # Combine eta and gamma effects
      obs_signal[[i]] <- list(
        region1 = matrix(
          rep(mu_1, n_voxel1) + rep(eta[i, 1:n_timept], n_voxel1) + gamma_r1[i, ] + noise[1:(n_voxel1 * n_timept)],
          nrow = n_timept, ncol = n_voxel1),
        region2 = matrix(
          rep(mu_2, n_voxel2) +
            rep(eta[i, (n_timept + 1):(2 * n_timept)], n_voxel2) +
            gamma_r2[i, ] + noise[(n_voxel1 * n_timept + 1):((n_voxel1 + n_voxel2) * n_timept)],
          nrow = n_timept, ncol = n_voxel2),
        region3 = matrix(
          rep(mu_3, n_voxel3) +
            rep(eta[i, (2 * n_timept + 1):(3 * n_timept)], n_voxel3) +
            gamma_r3[i, ] + noise[((n_voxel1 + n_voxel2) * n_timept + 1):(n_voxel_total * n_timept)],
          nrow = n_timept, ncol = n_voxel3)
      )
  }

  if (num_sim == 1)
    return(obs_signal[[1]])
  return(obs_signal)
}


#' Simulate signals from 3 regions using the new model
#'
#' @param num_sim number of simulations
#' @param voxel_coords List of region coordinates. Each row is a voxel and
#'   the three columns are x, y, and z coordinates.
#' @param n_timept number of timepoints
#' @param true_corr vector of inter-regional correlations (r12, r13, r23)
#' @param shared_params vector of parameters shared across regions (tau_eta, nugget)
#' @param region_params 3 x 7 dataframe of region-specific parameters.
#'    rows (region 1, region 2, region 3)
#'    columns (k_eta, phi_gamma, tau_gamma, k_gamma, nugget_gamma, mean, sigma2)
#' @param c_kernel_type Choice of spatial kernel. Defaul "matern_5_2".
#' @return obs_signal Simulated signal
generate_3_region_new <- function(
  num_sim, voxel_coords, n_timept, true_corr, region_params, shared_params,
  seed = 1, c_kernel_type = "matern_5_2") {
  set.seed(seed + 1000)

  stopifnot(n_timept >= 1)
  stopifnot(length(true_corr) == 3)
  stopifnot(length(shared_params) == 2)
  stopifnot(nrow(region_params) == 3 && ncol(region_params) == 7)
  stopifnot(c_kernel_type %in% c("matern_5_2", "rbf"))

  n_voxel1 <- nrow(voxel_coords[[1]])
  n_voxel2 <- nrow(voxel_coords[[2]])
  n_voxel3 <- nrow(voxel_coords[[3]])
  # Get distance and time matrices
  dist_sqrd_mat_region1 <- as.matrix(dist(voxel_coords[[1]]))^2
  dist_sqrd_mat_region2 <- as.matrix(dist(voxel_coords[[2]]))^2
  dist_sqrd_mat_region3 <- as.matrix(dist(voxel_coords[[3]]))^2
  timesqrd_mat <- (outer(1:n_timept, 1:n_timept, `-`))^2

  # Create covariance matrices

  corr_mat <- matrix(nrow = 3, ncol = 3)
  corr_mat[upper.tri(corr_mat)] <- true_corr
  corr_mat[lower.tri(corr_mat)] <- true_corr
  diag(corr_mat) <- 1

  # Covariance of eta effect
  A <- get_cor_mat("rbf", timesqrd_mat, shared_params["tau_eta"]) +
        shared_params["nugget"] * diag(n_timept)

  k_eta_diag <- diag(sqrt(region_params$k_eta))
  var_epsilon_sqrt <- diag(sqrt(region_params$sigma2))
  eta_sigma <- kronecker(var_epsilon_sqrt %*% k_eta_diag %*% corr_mat %*% k_eta_diag %*% var_epsilon_sqrt, A)

  # Matrices of gamma effects. C is spatial correlation matrix, B is temporal covariance matrix
  phi_gamma <- region_params$phi_gamma
  tau_gamma <- region_params$tau_gamma
  k_gamma <- region_params$k_gamma
  nugget_gamma <- region_params$nugget_gamma

  ## Region 1
  C1 <- get_cor_mat(c_kernel_type, dist_sqrd_mat_region1, phi_gamma[1])
  B1 <- k_gamma[1] * get_cor_mat("rbf", timesqrd_mat, tau_gamma[1]) +
    nugget_gamma[1] * diag(n_timept)
  gamma_sigma1 <- kronecker(C1, B1)

  ## Region 2
  C2 <- get_cor_mat(c_kernel_type, dist_sqrd_mat_region2, phi_gamma[2])
  B2 <- k_gamma[2] * get_cor_mat("rbf", timesqrd_mat, tau_gamma[2]) +
    nugget_gamma[2] * diag(n_timept)
  gamma_sigma2 <- kronecker(C2, B2)

  ## Region 3
  C3 <- get_cor_mat(c_kernel_type, dist_sqrd_mat_region3, phi_gamma[3])
  B3 <- k_gamma[3] * get_cor_mat("rbf", timesqrd_mat, tau_gamma[3]) +
    nugget_gamma[3] * diag(n_timept)
  gamma_sigma3 <- kronecker(C3, B3)

  eta <- MASS::mvrnorm(
          num_sim, mu = rep(0, 3 * n_timept), Sigma = eta_sigma) |>
          matrix(nrow = num_sim)
  gamma_r1 <- MASS::mvrnorm(
                num_sim, mu = rep(0, n_voxel1 * n_timept),
                Sigma = region_params$sigma2[1] * gamma_sigma1) |>
              matrix(nrow = num_sim)
  gamma_r2 <- MASS::mvrnorm(
                num_sim, mu = rep(0, n_voxel2 * n_timept),
                Sigma = region_params$sigma2[2] * gamma_sigma2) |>
              matrix(nrow = num_sim)
  gamma_r3 <- MASS::mvrnorm(
                num_sim, mu = rep(0, n_voxel3 * n_timept),
                Sigma = region_params$sigma2[3] * gamma_sigma3) |>
              matrix(nrow = num_sim)


  mu_1 <- region_params$mean[1]
  mu_2 <- region_params$mean[2]
  mu_3 <- region_params$mean[3]
  n_voxel_total <- n_voxel1 + n_voxel2 + n_voxel3
  obs_signal <- vector(mode = "list", length = num_sim)
  for (i in seq_along(obs_signal)) {
      # Generate iid noise
    noise <- c(rnorm(n_voxel1 * n_timept, sd = sqrt(region_params$sigma2[1])),
               rnorm(n_voxel2 * n_timept, sd = sqrt(region_params$sigma2[2])),
               rnorm(n_voxel3 * n_timept, sd = sqrt(region_params$sigma2[3])))

      # Combine eta and gamma effects
      obs_signal[[i]] <- list(
        region1 = matrix(
          rep(mu_1, n_voxel1) +
            rep(eta[i, 1:n_timept], n_voxel1) +
            gamma_r1[i, ] + noise[1:(n_voxel1 * n_timept)],
          nrow = n_timept, ncol = n_voxel1),
        region2 = matrix(
          rep(mu_2, n_voxel2) +
            rep(eta[i, (n_timept + 1):(2 * n_timept)], n_voxel2) +
            gamma_r2[i, ] + noise[(n_voxel1 * n_timept + 1):((n_voxel1 + n_voxel2) * n_timept)],
          nrow = n_timept, ncol = n_voxel2),
        region3 = matrix(
          rep(mu_3, n_voxel3) +
            rep(eta[i, (2 * n_timept + 1):(3 * n_timept)], n_voxel3) +
            gamma_r3[i, ] + noise[((n_voxel1 + n_voxel2) * n_timept + 1):(n_voxel_total * n_timept)],
          nrow = n_timept, ncol = n_voxel3)
      )
  }

  if (num_sim == 1)
    return(obs_signal[[1]])
  return(obs_signal)
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


# Generate p independent AR(2) processes with autocorrelation cor1 and cor2
generate_ar2 <- function(n_sim, t, p, cor1, cor2, scale, spatial_lower) {
  autocor <- function(h, xi) {
    x <- xi[1]
    y <- xi[2]
    scale * (x^(1-h)/(x^2-1) - y^(1-h)/(y^2-1)) / (x/(x^2-1) - y/(y^2 - 1))
  }

  opt <- function(xi, cor1, cor2) {
    c(autocor(1, xi) - cor1, autocor(2, xi) - cor2)
  }

  xi <- nleqslv::nleqslv(c(2, -3), \(xi) opt(xi, cor1, cor2))$x

  ar_coef <- function(xi) {
    x <- xi[1]
    y <- xi[2]
    c(1/x + 1/y, -1/(x * y))
  }
  phi <- ar_coef(xi)

  result <- matrix(nrow = n_sim, ncol = t * p)
  # Length of time series is t and single is a single obs. as a p x t matrix.
  for (j in 1:n_sim) {
    single <- matrix(nrow = p, ncol = t)
    for (i in 1:p) {
      single[i, ] <- arima.sim(n = t, model = list(order = c(2, 0, 0), ar = phi))
    }
    # Vectorize the matrix for a single observation.
    # Need to transpose here since time comes first.
    result[j, ] <- as.numeric(t(spatial_lower %*% single))
  }
  result
}

dwt <- function(region) {
    dwt_ts <- lapply(seq_len(ncol(region)), \(i) {
      moddwtY <- waveslim::dwt.nondyadic(region[, i])
      vdwtY <- na.omit(waveslim::brick.wall(moddwtY, wf = "la8", method = "dwt")$d4)
      vdwtY
    })
    Reduce(cbind, dwt_ts)
}

generate_ar2_region <- function(num_sim, voxel_coords, n_timept,
                                true_corr, region_params, shared_params,
                                seed = 1, c_kernel_type = "matern_5_2") {
  set.seed(seed + 1000)

  stopifnot(n_timept >= 1)
  stopifnot(length(true_corr) == 3)
  stopifnot(length(shared_params) == 2)
  stopifnot(nrow(region_params) == 3 && ncol(region_params) == 7)
  stopifnot(c_kernel_type %in% c("matern_5_2", "rbf"))

  n_voxel1 <- nrow(voxel_coords[[1]])
  n_voxel2 <- nrow(voxel_coords[[2]])
  n_voxel3 <- nrow(voxel_coords[[3]])
  # Get distance and time matrices
  dist_sqrd_mat_region1 <- as.matrix(dist(voxel_coords[[1]]))^2
  dist_sqrd_mat_region2 <- as.matrix(dist(voxel_coords[[2]]))^2
  dist_sqrd_mat_region3 <- as.matrix(dist(voxel_coords[[3]]))^2

  # Create covariance matrices

  corr_mat <- matrix(nrow = 3, ncol = 3)
  corr_mat[upper.tri(corr_mat)] <- true_corr
  corr_mat[lower.tri(corr_mat)] <- true_corr
  diag(corr_mat) <- 1

  k_eta_diag <- diag(sqrt(region_params$k_eta))
  var_epsilon_sqrt <- diag(sqrt(region_params$sigma2))
  eta_spatial <- var_epsilon_sqrt %*% k_eta_diag %*% corr_mat %*% k_eta_diag %*% var_epsilon_sqrt
  eta_lower <- t(chol(eta_spatial))
  eta <- generate_ar2(num_sim, n_timept, 3, 0.4, 0.3, 1, eta_lower)

  phi_gamma <- region_params$phi_gamma
  k_gamma <- region_params$k_gamma # TODO is this a scale for ar2?

  ## Region 1
  C1 <- get_cor_mat(c_kernel_type, dist_sqrd_mat_region1, phi_gamma[1])
  gamma_r1 <- generate_ar2(num_sim, n_timept, n_voxel1, 0.4, 0.3, 1, t(chol(C1)))

  ## Region 2
  C2 <- get_cor_mat(c_kernel_type, dist_sqrd_mat_region2, phi_gamma[2])
  gamma_r2 <- generate_ar2(num_sim, n_timept, n_voxel2, 0.4, 0.3, 1, t(chol(C2)))

  ## Region 3
  C3 <- get_cor_mat(c_kernel_type, dist_sqrd_mat_region3, phi_gamma[3])
  gamma_r3 <- generate_ar2(num_sim, n_timept, n_voxel3, 0.4, 0.3, 1, t(chol(C3)))

  mu_1 <- region_params$mean[1]
  mu_2 <- region_params$mean[2]
  mu_3 <- region_params$mean[3]
  n_voxel_total <- n_voxel1 + n_voxel2 + n_voxel3
  obs_signal <- vector(mode = "list", length = num_sim)
  for (i in seq_along(obs_signal)) {
    # Generate iid noise
    noise <- c(rnorm(n_voxel1 * n_timept, sd = sqrt(region_params$sigma2[1])),
               rnorm(n_voxel2 * n_timept, sd = sqrt(region_params$sigma2[2])),
               rnorm(n_voxel3 * n_timept, sd = sqrt(region_params$sigma2[3])))

    region1 <-
      matrix(rep(mu_1, n_voxel1) + rep(eta[i, 1:n_timept], n_voxel1)
             + gamma_r1[i, ] + noise[1:(n_voxel1 * n_timept)],
             nrow = n_timept, ncol = n_voxel1)

    region2 <-
      matrix(rep(mu_2, n_voxel2) + rep(eta[i, (n_timept + 1):(2 * n_timept)], n_voxel2)
             + gamma_r2[i, ] + noise[(n_voxel1 * n_timept + 1):((n_voxel1 + n_voxel2) * n_timept)],
             nrow = n_timept, ncol = n_voxel2)

    region3 <-
      matrix(rep(mu_3, n_voxel3) + rep(eta[i, (2 * n_timept + 1):(3 * n_timept)], n_voxel3)
             + gamma_r3[i, ] + noise[((n_voxel1 + n_voxel2) * n_timept + 1):(n_voxel_total * n_timept)],
             nrow = n_timept, ncol = n_voxel3)

    # Combine eta and gamma effects
    obs_signal[[i]] <- list(
      region1 = dwt(region1),
      region2 = dwt(region2),
      region3 = dwt(region3))
  }

  if (num_sim == 1)
    return(obs_signal[[1]])
  return(obs_signal)
}