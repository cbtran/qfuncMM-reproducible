# Run stage 2 for all pairs of regions in the correctly specified setting
# with true stage 1 parameters plugged in.
# Run this script in the terminal as
# >Rscript full-run/full-run-oracle.R <delta> <psi> <startid> <endid>

library(qfuncMM)

set.seed(100)

args <- commandArgs(trailingOnly = TRUE)
setting <- paste0(args[1], "-", args[2], "-M60-100-std")
startid <- args[3]
endid <- args[4]

std_result_path <- paste0("full-run/out/", setting, "-result.rds")
if (!file.exists(std_result_path)) {
  stop(std_result_path, " does not exist. Run the correct specification before the oracle.")
}
std_result <- readRDS(std_result_path)

voxel_coords <- readRDS("full-run/rat_coords.rds")
allsignals <- readRDS(paste0("full-run/data/", setting, ".rds"))

num_timept <- nrow(allsignals$data[[1]]$region1)
time_sqrd_mat <- outer(seq_len(num_timept), seq_len(num_timept), `-`)^2

stage1_names <- c("phi_gamma", "tau_gamma", "k_gamma", "nugget_gamma", "var_noise")
stage1_true <-
  allsignals$setting$region_parameters[, stage1_names] |>
  as.matrix()

run <- function(signal, runid) {
  n_region <- length(signal)
  rho_mx <- matrix(1, nrow = n_region, ncol = n_region,
                   dimnames = list(paste0("r", seq_len(n_region)),
                                   paste0("r", seq_len(n_region))))
  stage2_inter <- array(dim = c(n_region, n_region, 4),
                        dimnames = list(paste0("r", seq_len(n_region)),
                                        paste0("r", seq_len(n_region)),
                                        c("k_eta1", "k_eta2",
                                          "tau_eta", "nugget_eta")))
  # Stage 1 param list: phi_gamma, tau_gamma, k_gamma, nugget_gamma, var_noise
  for (reg1 in seq_len(n_region - 1)) {
    for (reg2 in seq(reg1 + 1, n_region)) {
      pairname <- paste0("r", reg1, reg2)
      eblue_init <- std_result$rho["rho_eblue", pairname, runid]
      stage2_result <- qfuncMM:::fit_inter_model(signal[[reg1]],
                                                 voxel_coords[[reg1]],
                                                 signal[[reg2]],
                                                 voxel_coords[[reg2]],
                                                 time_sqrd_mat,
                                                 stage1_true[reg1, ],
                                                 stage1_true[reg2, ],
                                                 eblue_init,
                                                 3L)
      rho_mx[reg1, reg2] <- stage2_result["rho"]
      rho_mx[reg2, reg1] <- stage2_result["rho"]
      stage2_inter[reg1, reg2, ] <- stage2_result[-1]
      stage2_inter[reg2, reg1, ] <- stage2_inter[reg1, reg2, ]
      message("Finished region pair ", reg1, "-", reg2, "\n")
    }
  }
  stage2 <- stage2_inter |>
    apply(3, function(x) c(x[1, 2], x[1, 3], x[2, 3])) |>
    t()
  colnames(stage2) <- c("r12", "r13", "r23")
  return(list(rho = rho_mx, stage2 = stage2))
}

dataids <- seq(startid, endid)
nsim <- length(dataids)
rho_all <- matrix(nrow = nsim, ncol = 3, dimnames = list(NULL, c("r12", "r13", "r23")))
stage2_all <- array(dim = c(4, 3, nsim),
                    dimnames = list(c("k_eta1", "k_eta2", "tau_eta", "nugget_eta"),
                                    c("r12", "r13", "r23"), NULL))
for (i in dataids) {
  d <- allsignals$data[[i]]
  run_result <- run(d, i)
  runid <- i - dataids[1] + 1
  rho_all[i, ] <- c(run_result$rho[1, 2], run_result$rho[1, 3], run_result$rho[2, 3])
  stage2_all[, , runid] <- run_result$stage2
  saveRDS(list(rho = rho_all[1:runid, ], stage2 = stage2_all[, , 1:runid]),
          sprintf("full-run/out/%s-result-oracle.rds", setting))
  message("Finished sim ", i)
}
message("Results saved to ", sprintf("full-run/out/%s-result-oracle.rds", setting))
