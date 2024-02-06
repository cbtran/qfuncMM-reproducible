# Run stage 2 for all pairs of regions in the correctly specified setting
# with true stage 1 parameters plugged in.

library(qfuncMM)
library(glue)

set.seed(1001)

args <- commandArgs(trailingOnly = TRUE)
setting <- paste0(args[1], "-", args[2], "-M60-100-rat")
startid <- args[3]
endid <- args[4]

voxel_coords <- readRDS("full-run/rat_coords.rds")
allsignals <- readRDS(paste0("full-run/data/", setting, ".rds"))
num_timept <- nrow(allsignals$data[[1]][[1]])
num_voxel <- sapply(voxel_coords, nrow)
time_sqrd_mat <- outer(seq_len(num_timept), seq_len(num_timept), `-`)^2

error_na <- list(theta = rep(NA, 5), var_noise = rep(NA, 2))

stage1_paramnames <- c("phi_gamma", "tau_gamma", "k_gamma", "nugget_gamma", "sigma2")
stage1_true <-
  allsignals$setting$region_parameters[, stage1_paramnames] |>
  as.matrix()

run <- function(signal, runid) {
  n_region <- length(signal)
  cor_mx <- matrix(1, nrow = n_region, ncol = n_region,
                   dimnames = list(paste0("r", seq_len(n_region)),
                                   paste0("r", seq_len(n_region))))
  stage2_inter <- array(dim = c(n_region, n_region, 6),
                        dimnames = list(paste0("r", seq_len(n_region)),
                                        paste0("r", seq_len(n_region)),
                                        c("k_eta1", "k_eta2",
                                          "tau_eta", "nugget_eta",
                                          "var_noise1", "var_noise2")))
  # Stage 1 param list: phi_gamma, tau_gamma, k_gamma, nugget_gamma, var_noise
  tictoc::tic(paste("Finished run", runid))
  for (reg1 in seq_len(n_region)) {
    for (reg2 in seq_len(reg1 - 1)) {
      stage2_result <- fit_inter_model(signal[[reg1]],
                                       voxel_coords[[reg1]],
                                       signal[[reg2]],
                                       voxel_coords[[reg2]],
                                       time_sqrd_mat,
                                       stage1_true[reg1, ],
                                       stage1_true[reg2, ],
                                       3L)
      cor_mx[reg1, reg2] <- stage2_result$rho
      cor_mx[reg2, reg1] <- stage2_result$rho
      stage2_inter[reg1, reg2, ] <- as.numeric(stage2_result[-1])
      stage2_inter[reg2, reg1, ] <- stage2_inter[reg1, reg2, ]
      cat("Finished region pair", reg1, "-", reg2, "\n")
    }
  }
  tictoc::toc()
  stage2 <- abind::abind(cor_mx, stage2_inter) |>
    apply(3, function(x) c(x[2, 1], x[3, 1], x[3, 2])) |>
    t()
  rownames(stage2) <- c("rho", rownames(stage2)[-1])
  colnames(stage2) <- c("r12", "r13", "r23")

  return(list(stage2 = stage2))
}

dataids <- startid:endid
nsim <- length(dataids)
results <- array(dim = c(7, 3, nsim))
dimnames(results) <- list(
  c("rho", "kEta1", "kEta2", "tauEta", "nugget", "var_noise_prof1", "var_noise_prof2"),
  c("r12", "r13", "r23"), NULL)

warn_log <- file("full-run/warnings.log", open = "wt")
sink(warn_log, append = TRUE, type = "message")

for (i in dataids) {
  signal <- allsignals$data[[i]]
  run_result <- run(signal, i)
  runid <- i - dataids[1] + 1
  results[, , runid] <- run_result$stage2
  saveRDS(list(stage2 = results[, , 1:runid]),
          paste0("full-run/out-oracle/", setting, glue("-result-{startid}-{endid}.rds")))
  cat("Finished sim", i, "\n")
}

sink(type = "message")
close(warn_log)