# Run this script in the terminal as
# >Rscript full-run/full-run.R <delta> <psi> <spec> <startid> <endid>
# where <delta> and <psi> is one of "high", "mid", "low",
# and <spec> is one of "std", "fgn", "ar2", "anisotropic",
# and <startid> and <endid> are between 1 and 100.

library(qfuncMM)

set.seed(100)

args <- commandArgs(trailingOnly = TRUE)
delta <- args[1]
psi <- args[2]
spec <- args[3]
startid <- args[4]
endid <- args[5]
dataids <- seq(startid, endid)
nsim <- length(dataids)

setting <- paste0(delta, "-", psi, "-M60-100-", spec)

voxel_coords <- readRDS("full-run/rat_coords.rds")
datapath <- paste0("full-run/data/", setting, ".rds")
if (!file.exists(datapath))
  stop(sprintf("%s not found. Generate the data first."))
all_data <- readRDS(datapath)
if (!dir.exists("full-run/out"))
  dir.create("full-run/out")
outpath <- paste0("full-run/out/", setting, "-result.rds")

message(sprintf("Running spec=%s, delta=%s, psi=%s, %d simulations\n",
                spec, delta, psi, nsim))

run <- function(signal, runid) {
  # Stage 1 param list: phi_gamma, tau_gamma, k_gamma, nugget_gamma, var_noise
  result <- qfuncMM(signal, voxel_coords, verbose = TRUE)
  stage2 <- result$stage2 |>
    apply(3, function(x) c(x[1, 2], x[1, 3], x[2, 3])) |>
    t()
  colnames(stage2) <- c("r12", "r13", "r23")

  rho_arr <- lapply(list(result$rho, result$rho_eblue, result$rho_ca),
      function(x) c(x[1, 2], x[1, 3], x[2, 3]))
  rho_arr <- Reduce(rbind, rho_arr)
  rownames(rho_arr) <- names(result)[1:3]
  colnames(rho_arr) <- c("r12", "r13", "r23")

  return(list(
    rho = rho_arr,
    stage1 = t(result$stage1),
    stage2 = stage2
  ))
}

results_rho <- array(dim = c(3, 3, nsim), dimnames = list(
  c("rho", "rho_eblue", "rho_ca"),
  c("r12", "r13", "r23"), NULL))
results_stage2 <- array(dim = c(4, 3, nsim))
dimnames(results_stage2) <- list(
  c("k_eta1", "k_eta2", "tau_eta", "nugget_eta"),
  c("r12", "r13", "r23"), NULL)
results_stage1 <- array(dim = c(5, 3, nsim))
dimnames(results_stage1) <- list(
  c("phi_gamma", "tau_gamma", "k_gamma", "nugget_gamma", "var_noise"),
  c("r1", "r2", "r3"), NULL)

for (i in dataids) {
  d <- all_data$data[[i]]
  run_result <- run(d, i)
  runid <- i - dataids[1] + 1
  results_rho[, , runid] <- run_result$rho
  results_stage1[, , runid] <- run_result$stage1
  results_stage2[, , runid] <- run_result$stage2
  saveRDS(list(rho = results_rho[, , 1:runid],
               stage1 = results_stage1[, , 1:runid],
               stage2 = results_stage2[, , 1:runid]),
          outpath)
  message("Finished sim ", i)
}
saveRDS(list(rho = results_rho, stage1 = results_stage1, stage2 = results_stage2), outpath)
message("Results saved to ", outpath)
