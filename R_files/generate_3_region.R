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
#' @param region_params 3 x 4 dataframe of region-specific parameters.
#'    rows (region 1, region 2, region 3)
#'    columns (phi_gamma, tau_gamma, k_gamma, mean)
#' @param sigma2 noise variance
#' @param c_kernel_type Choice of spatial kernel. Defaul "matern_5_2".
#' @return obs_signal Simulated signal
generate_3_region <- function(
  num_sim, voxel_coords, n_timept, true_corr, region_params, shared_params,
  sigma2 = 1, seed = 1, c_kernel_type = "matern_5_2") {

  stopifnot(n_timept >= 1)
  stopifnot(length(true_corr) == 3)
  stopifnot(length(shared_params) == 3)
  stopifnot(nrow(region_params) == 3 && ncol(region_params) == 4)
  stopifnot(sigma2 >= 0)
  stopifnot(c_kernel_type %in% c("matern_5_2", "rbf"))

  n_voxel1 <- nrow(voxel_coords[[1]])
  n_voxel2 <- nrow(voxel_coords[[2]])
  n_voxel3 <- nrow(voxel_coords[[3]])
  # Get distance and time matrices
  dist_sqrd_mat_region1 <- as.matrix(dist(voxel_coords[[1]]))
  dist_sqrd_mat_region2 <- as.matrix(dist(voxel_coords[[2]]))
  dist_sqrd_mat_region3 <- as.matrix(dist(voxel_coords[[3]]))
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

  ## Region 1
  C1 <- get_cor_mat(c_kernel_type, dist_sqrd_mat_region1, phi_gamma[1])
  B1 <- k_gamma[1] * get_cor_mat("rbf", timesqrd_mat, tau_gamma[1])
  gamma_sigma1 <- kronecker(C1, B1)

  ## Region 2
  C2 <- get_cor_mat(c_kernel_type, dist_sqrd_mat_region2, phi_gamma[2])
  B2 <- k_gamma[2] * get_cor_mat("rbf", timesqrd_mat, tau_gamma[2])
  gamma_sigma2 <- kronecker(C2, B2)

  ## Region 3
  C3 <- get_cor_mat(c_kernel_type, dist_sqrd_mat_region3, phi_gamma[3])
  B3 <- k_gamma[3] * get_cor_mat("rbf", timesqrd_mat, tau_gamma[3])
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
