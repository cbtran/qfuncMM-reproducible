# Run this script in the terminal as
# >Rscript full-run/full-run.R <delta> <phi> <spec> <startid> <endid>
# where <delta> and <phi> is one of "high", "mid", "low",
# and <spec> is one of "std", "fgn", "ar2", "anisotropic",
# and <startid> and <endid> are between 1 and 100.

library(qfuncMM)
library(glue)

set.seed(100)

args <- commandArgs(trailingOnly = TRUE)
setting <- paste0(args[1], "-", args[2], "-M60-100-rat-", args[3])
if (args[3] == "std")
  setting <- paste0(args[1], "-", args[2], "-M60-100-rat")
startid <- args[4]
endid <- args[5]

voxel_coords <- readRDS("full-run/rat_coords.rds")
allsignals <- readRDS(paste0("full-run/data/", setting, ".rds"))
num_timept <- nrow(allsignals$data[[1]][[1]])
num_voxel <- sapply(voxel_coords, nrow)
time_sqrd_mat <- outer(seq_len(num_timept), seq_len(num_timept), `-`)^2

error_na <- list(theta = rep(NA, 5), var_noise = rep(NA, 2))

outpath <- paste0("full-run/out/", setting, sprintf("-result-%d-%d.rds", startid, endid))

run <- function(signal, runid) {
  # Stage 1 param list: phi_gamma, tau_gamma, k_gamma, nugget_gamma, var_noise
  # tictoc::tic(paste("Finished run", runid))
  result <- qfuncMM(signal, voxel_coords, verbose = TRUE)
  # tictoc::toc()
  stage2 <- abind::abind(result$rho, result$stage2) |>
    apply(3, function(x) c(x[2, 1], x[3, 1], x[3, 2])) |>
    t()
  rownames(stage2) <- c("rho", rownames(stage2)[-1])
  colnames(stage2) <- c("r12", "r13", "r23")

  return(list(
    stage1 = t(result$stage1),
    stage2 = stage2
  ))
}

dataids <- startid:endid
nsim <- length(dataids)
results <- array(dim = c(7, 3, nsim))
dimnames(results) <- list(
  c("rho", "kEta1", "kEta2", "tauEta", "nugget", "var_noise_prof1", "var_noise_prof2"),
  c("r12", "r13", "r23"), NULL)
results_stage1 <- array(dim = c(5, 3, nsim))
dimnames(results_stage1) <- list(
  c("phi_gamma", "tau_gamma", "k_gamma", "nugget_gamma", "var_noise"),
  c("r1", "r2", "r3"), NULL)

for (i in dataids) {
  signal <- allsignals$data[[i]]
  run_result <- run(signal, i)
  runid <- i - dataids[1] + 1
  results_stage1[, , runid] <- run_result$stage1
  results[, , runid] <- run_result$stage2
  saveRDS(list(stage1 = results_stage1[, , 1:runid],
               stage2 = results[, , 1:runid]),
          outpath)
  cat("Finished sim", i, "\n")
}
saveRDS(list(stage1 = results_stage1, stage2 = results), outpath)
