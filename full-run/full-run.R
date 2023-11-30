library(qfuncMM)
source("full-run/utils.R")

set.seed(100)

args <- commandArgs(trailingOnly = TRUE)
setting <- paste0(args[1], "-", args[2], "-M60-100-", args[3])

voxel_coords <- readRDS("full-run/rat_coords.rds")
allsignals <- readRDS(paste0("full-run/", setting, ".rds"))
num_timept <- nrow(allsignals$data[[1]][[1]])
num_voxel <- sapply(voxel_coords, nrow)
time_sqrd_mat <- outer(seq_len(num_timept), seq_len(num_timept), `-`)^2

error_na <- list(theta = rep(NA, 5), var_noise = rep(NA, 2))

run <- function(signal, runid) {
  region1_mx <- matrix(signal$region1, nrow = num_timept, ncol = num_voxel[1])
  region2_mx <- matrix(signal$region2, nrow = num_timept, ncol = num_voxel[2])
  region3_mx <- matrix(signal$region3, nrow = num_timept, ncol = num_voxel[3])
  coords1 <- voxel_coords$r1
  coords2 <- voxel_coords$r2
  coords3 <- voxel_coords$r3

  # Stage 1 param list: phi_gamma, tau_gamma, k_gamma, nugget_gamma, var_noise
  tictoc::tic("Finished intra all regions")
  stage1_region1 <- opt_intra_new(
    c(0, 0, 0, 0), matrix(region1_mx, ncol = 1),
    coords1, time_sqrd_mat, num_voxel[1], num_timept, kernel_dict("matern_5_2"))
  stage1params_r1 <- c(stage1_region1$theta, stage1_region1$var_noise)

  stage1_region2 <- opt_intra_new(
    c(0, 0, 0, 0), matrix(region2_mx, ncol = 1),
    coords2, time_sqrd_mat, num_voxel[2], num_timept, kernel_dict("matern_5_2"))
  stage1params_r2 <- c(stage1_region2$theta, stage1_region2$var_noise)

  stage1_region3 <- opt_intra_new(
    c(0, 0, 0, 0), matrix(region3_mx, ncol = 1),
    coords3, time_sqrd_mat, num_voxel[3], num_timept, kernel_dict("matern_5_2"))
  stage1params_r3 <- c(stage1_region3$theta, stage1_region3$var_noise)
  tictoc::toc()

  corr_avg <- computeCA(signal)
  theta_init <- c(softminus(1), softminus(1), 0, softminus(0.1))
  # phi_gamma, tau_gamma, k_gamma, nugget_gamma, noise_variance
  cat("Starting inter with corr_avg", sigmoid_inv(corr_avg), "\n")
  tictoc::tic("Finished region 12")
  region12 <- tryCatch(
    opt_inter_new(
      c(corr_avg[1], theta_init),
      region1_mx, region2_mx, coords1, coords2, time_sqrd_mat,
      stage1params_r1, stage1params_r2, kernel_dict("matern_5_2")
    ),
    warning = function(w) {
      warning(paste0("Warning in run ", runid, " region 12: ", w$message))
      return(error_na)
    },
    error = function(e) {
      warning(paste0("Error in run ", runid, " region 12: ", e$message))
      return(error_na)
    }
  )
  tictoc::toc()
  tictoc::tic("Finished region 13")
  region13 <- tryCatch(
    opt_inter_new(
      c(corr_avg[2], theta_init),
      region1_mx, region3_mx, coords1, coords3, time_sqrd_mat,
      stage1params_r1, stage1params_r3, kernel_dict("matern_5_2")
    ),
    warning = function(w) {
      warning(paste0("Warning in run ", runid, " region 13: ", w$message))
      return(error_na)
    },
    error = function(e) {
      warning(paste0("Error in run ", runid, " region 13: ", e$message))
      return(error_na)
    }
  )
  tictoc::toc()
  tictoc::tic("Finished region 23")
  region23 <- tryCatch(
    opt_inter_new(
      c(corr_avg[3], theta_init),
      region2_mx, region3_mx, coords2, coords3, time_sqrd_mat,
      stage1params_r2, stage1params_r3, kernel_dict("matern_5_2")
    ),
    warning = function(w) {
      warning(paste0("Warning in run ", runid, " region 23: ", w$message))
      return(error_na)
    },
    error = function(e) {
      warning(paste0("Error in run ", runid, " region 23: ", e$message))
      return(error_na)
    }
  )
  tictoc::toc()
  return(list(
    stage1 = c(stage1params_r1, stage1params_r2, stage1params_r3),
    stage2 = c(c(region12$theta, region12$var_noise),
               c(region13$theta, region13$var_noise),
               c(region23$theta, region23$var_noise)))
  )
}

nsim <- 100
results <- array(dim = c(7, 3, nsim))
dimnames(results) <- list(
  c("rho", "kEta1", "kEta2", "tauEta", "nugget", "var_noise_prof1", "var_noise_prof2"),
  c("r12", "r13", "r23"), NULL)
results_stage1 <- array(dim = c(5, 3, nsim))
dimnames(results_stage1) <- list(
  c("phi_gamma", "tau_gamma", "k_gamma", "nugget_gamma", "var_noise"),
  c("r1", "r2", "r3"), NULL)

warn_log <- file("full-run/warnings.log", open = "wt")
sink(warn_log, append = TRUE, type = "message")

runids <- seq_len(nsim)
for (i in runids) {
  signal <- allsignals$data[[i]]
  run_result <- run(signal, i)
  results_stage1[, , i] <- run_result$stage1
  results[, , i] <- run_result$stage2
  saveRDS(list(stage1 = results_stage1[, , 1:i], stage2 = results[, , 1:i]),
          paste0("full-run/", setting, "-result.rds"))
  cat("Finished sim", i, "\n")
}
saveRDS(list(stage1 = results_stage1, stage2 = results),
        paste0("full-run/", setting, "-result.rds"))

sink(type = "message")
close(warn_log)
