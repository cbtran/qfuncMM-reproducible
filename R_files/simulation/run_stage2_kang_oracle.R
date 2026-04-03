args <- commandArgs(trailingOnly = TRUE)
# args <- c("std", "noisy", "data", "out", 100)
data_spec <- args[1]
cov_setting <- args[2]
data_dir <- args[3]
out_dir <- args[4]
if (length(args) < 4) {
  stop(
    "Usage: Rscript run_stage2_kang_oracle.R ",
    "<data_spec> <cov_setting> <data_dir> <out_dir>"
  )
}
if (!cov_setting %in% c("noiseless", "noisy")) {
  stop("cov_setting must be either 'noiseless' or 'noisy'")
}

# ---------------------------------------------------------------------------
# Helper: Matérn-5/2 and RBF kernels (mirrors qfuncMM/R/stage1_init.R:1-7)
# ---------------------------------------------------------------------------
rbf_fn <- function(tau, time_sqrd_mat) {
  exp(-tau^2 / 2 * time_sqrd_mat)
}

matern_fn <- function(phi, dist_mat) {
  (1 + phi * sqrt(5) * dist_mat + phi^2 * (5 / 3) * dist_mat^2) *
    exp(-phi * sqrt(5) * dist_mat)
}

# ---------------------------------------------------------------------------
# Helper: EBLUP of gamma_j via Kronecker eigendecomposition
#
# Computes gamma_hat = Lambda_j (Lambda_j + sigma2_ep I)^{-1} (X_j - U_j eblue)
# where Lambda_j = C_j (x) B_j, C_j is L×L spatial (Matern-5/2),
# B_j is M×M temporal (RBF + nugget).
#
# Returns M×L matrix.
# ---------------------------------------------------------------------------
compute_eblup_gamma <- function(data_std, coords, eblue, stage1) {
  M <- nrow(data_std)
  L <- ncol(data_std)

  dist_mat <- as.matrix(dist(coords))
  time_sqrd <- outer(seq_len(M), seq_len(M), function(i, j) (i - j)^2)

  C_j <- matern_fn(stage1$phi_gamma, dist_mat)
  B_j <- stage1$k_gamma * rbf_fn(stage1$tau_gamma, time_sqrd) +
    stage1$nugget_gamma * diag(M)

  eig_C <- eigen(C_j, symmetric = TRUE)
  eig_B <- eigen(B_j, symmetric = TRUE)

  # Eigenvalues of Lambda = C (x) B: outer(B-eigs, C-eigs) gives M×L matrix
  # consistent with column-major vec ordering (fast index = time/B, slow = voxel/C)
  lambda_eig <- outer(eig_B$values, eig_C$values)
  weights <- lambda_eig / (lambda_eig + stage1$sigma2_ep)

  # Residual: X - U * eta_hat  (U * eta = eblue repeated L times column-wise)
  R <- data_std - matrix(eblue, nrow = M, ncol = L)

  # Apply EBLUP via Kronecker identity: (C (x) B) vec(X) = vec(B X C^T)
  R_eig <- t(eig_B$vectors) %*% R %*% eig_C$vectors
  R_eig_weighted <- R_eig * weights
  gamma_hat <- eig_B$vectors %*% R_eig_weighted %*% t(eig_C$vectors)

  return(gamma_hat)
}

# ---------------------------------------------------------------------------
# Helper: Kang connectivity estimate (formula (1) in kang_estimator.md)
#
# r_star_1, r_star_2: M×L_j residual matrices (rows = timepoints, cols = voxels)
# Returns scalar rho_kang.
# ---------------------------------------------------------------------------
compute_kang_estimate <- function(r_star_1, r_star_2) {
  cov_cross <- cov(r_star_1, r_star_2) # L1×L2 cross-covariance
  cov_w1 <- cov(r_star_1) # L1×L1 within-region covariance
  cov_w2 <- cov(r_star_2) # L2×L2 within-region covariance

  numer <- mean(cov_cross)
  within_avg_1 <- mean(cov_w1[upper.tri(cov_w1)])
  within_avg_2 <- mean(cov_w2[upper.tri(cov_w2)])

  numer / sqrt(within_avg_1 * within_avg_2)
}

# ---------------------------------------------------------------------------
# Setup directories
# ---------------------------------------------------------------------------
data_spec_dir <- file.path(data_dir, data_spec)
if (!dir.exists(data_spec_dir)) {
  stop(paste0(
    "Data specification '", data_spec, "' does not exist in directory '", data_dir,
    "'. Generate the data first."
  ))
}

stage1_dir <- file.path(out_dir, data_spec, "stage1")
if (cov_setting == "noiseless") {
  stage1_dir <- paste0(stage1_dir, "_noiseless")
}
if (!dir.exists(stage1_dir)) {
  stop(paste("Stage 1 directory", normalizePath(stage1_dir, mustWork = FALSE), "does not exist. Run stage 1 first."))
}

out_dir_spec <- file.path(out_dir, data_spec, "stage2_kang_oracle")
if (cov_setting == "noiseless") {
  out_dir_spec <- paste0(out_dir_spec, "_noiseless")
}
dir.create(out_dir_spec, recursive = TRUE, showWarnings = FALSE)

coords_file <- if (data_spec == "std-hcp") "s111716_coords.rds" else "rat_coords.rds"
voxel_coords <- readRDS(file.path("R_files", "simulation", coords_file))

# ---------------------------------------------------------------------------
# Process all 9 delta/psi settings
# ---------------------------------------------------------------------------
all_settings <- expand.grid(
  delta = c("high", "mid", "low"),
  psi = c("high", "mid", "low"),
  stringsAsFactors = FALSE
)

for (s in seq_len(nrow(all_settings))) {
  delta_setting <- all_settings$delta[s]
  psi_setting <- all_settings$psi[s]
  setting_str <- paste0(delta_setting, "-", psi_setting)
  data_setting <- readRDS(file.path(data_spec_dir, paste0(setting_str, ".rds")))
  csv_file <- file.path(out_dir_spec, paste0("results_", setting_str, ".csv"))

  if (file.exists(csv_file)) {
    message(paste0("CSV file ", csv_file, " already exists. Skipping."))
    next
  }

  # True data-generation parameters for this delta/psi setting
  true_params <- data_setting$setting

  for (simid in seq_along(data_setting$data)) {
    d <- data_setting$data[[simid]]
    sim_name <- paste0(setting_str, "-", simid)

    s1_outfiles <- list.files(
      stage1_dir,
      pattern    = paste0("qfuncMM_stage1_intra_region_", sim_name, "_[0-9]+\\.json$"),
      full.names = TRUE
    )

    reg_pairs <- list(c(1, 2), c(1, 3), c(2, 3))
    for (rp in reg_pairs) {
      r1_id <- rp[1]
      r2_id <- rp[2]

      r1_data <- d[[r1_id]]
      r2_data <- d[[r2_id]]
      r1_data_std <- (r1_data - mean(r1_data)) / stats::sd(r1_data)
      r2_data_std <- (r2_data - mean(r2_data)) / stats::sd(r2_data)
      r1_coords <- voxel_coords[[r1_id]]
      r2_coords <- voxel_coords[[r2_id]]

      # Build oracle stage1 parameter lists from true generation parameters.
      # phi_gamma and tau_gamma are scale-free and need no adjustment.
      # k_gamma, nugget_gamma, and var_noise are all variance parameters and
      # must be divided by var(raw_data) to match the standardized data scale.
      r1_var <- stats::var(as.vector(r1_data))
      r1_oracle <- as.list(true_params$region_parameters[r1_id, ])
      r1_oracle$k_gamma      <- r1_oracle$k_gamma      / r1_var
      r1_oracle$nugget_gamma <- r1_oracle$nugget_gamma / r1_var
      r1_oracle$sigma2_ep    <- r1_oracle$var_noise    / r1_var

      r2_var <- stats::var(as.vector(r2_data))
      r2_oracle <- as.list(true_params$region_parameters[r2_id, ])
      r2_oracle$k_gamma      <- r2_oracle$k_gamma      / r2_var
      r2_oracle$nugget_gamma <- r2_oracle$nugget_gamma / r2_var
      r2_oracle$sigma2_ep    <- r2_oracle$var_noise    / r2_var

      # EBLUE still comes from the Stage 1 JSON (estimated, not oracle)
      j1 <- jsonlite::fromJSON(s1_outfiles[r1_id])
      j2 <- jsonlite::fromJSON(s1_outfiles[r2_id])

      t_start <- proc.time()[["elapsed"]]

      gamma_hat_1 <- compute_eblup_gamma(r1_data_std, r1_coords, j1$eblue, r1_oracle)
      gamma_hat_2 <- compute_eblup_gamma(r2_data_std, r2_coords, j2$eblue, r2_oracle)

      mu_hat_1 <- mean(j1$eblue)
      mu_hat_2 <- mean(j2$eblue)

      r_star_1 <- r1_data_std - mu_hat_1 - gamma_hat_1
      r_star_2 <- r2_data_std - mu_hat_2 - gamma_hat_2

      rho_kang <- compute_kang_estimate(r_star_1, r_star_2)
      wall_time <- proc.time()[["elapsed"]] - t_start

      result_row <- data.frame(
        simid = simid,
        region1_uniqid = r1_id,
        region2_uniqid = r2_id,
        rho_kang = rho_kang,
        wall_time_s = wall_time,
        stringsAsFactors = FALSE
      )

      first_write <- !file.exists(csv_file) || file.info(csv_file)$size == 0
      write.table(
        result_row,
        file      = csv_file,
        sep       = ",",
        row.names = FALSE,
        col.names = first_write,
        append    = !first_write
      )
    }
  }
} # end settings loop
