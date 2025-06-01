library(qfuncMM)
library(tibble)
library(readr)
RhpcBLASctl::blas_set_num_threads(10)
RhpcBLASctl::omp_set_num_threads(1)

args <- commandArgs(trailingOnly = TRUE)
# args <- c("out", "std", "stage2_vecchia", "vecchia")
outdir <- args[1]
dataspec <- args[2]
s2_run_str <- args[3]
method <- args[4]
stopifnot(method %in% c("reml", "vecchia"))

s1_dir <- file.path(outdir, dataspec, "stage1")
s2_dir <- file.path(outdir, dataspec, s2_run_str)

coords <- readRDS("R_files/rat_coords.rds")
for (delta in c("high", "mid", "low")) {
  for (psi in c("high", "mid", "low")) {
    message(sprintf("Computing asymptotic confidence intervals for delta = %s, psi = %s", delta, psi))

    setting_str <- paste0(delta, "-", psi)
    outfile <- file.path(s2_dir, sprintf("asymp_ci_%s-%s.csv", delta, psi))
    if (file.exists(outfile)) {
      message(sprintf("Output file %s already exists, skipping...", outfile))
      next
    }
    stage2_out <- read_csv(file.path(s2_dir, sprintf("results_%s.csv", setting_str)), col_types = "iii")
    asymp_inf_mx <- matrix(nrow = nrow(stage2_out), ncol = 8)
    colnames(asymp_inf_mx) <- c(
      "simid", "region1_uniqid", "region2_uniqid",
      "rho_true", "rho_hat", "asymp_var_rho", "lower_95", "upper_95"
    )

    for (id in seq_len(nrow(stage2_out))) {
      cat(id, " ")
      sim_info <- stage2_out[id, ]
      stage1_r1_file <- file.path(
        s1_dir,
        sprintf("qfuncMM_stage1_intra_region_%s-%s-%d_%d.json", delta, psi, sim_info$simid, sim_info$region1_uniqid)
      )
      stage1_r2_file <- file.path(
        s1_dir,
        sprintf("qfuncMM_stage1_intra_region_%s-%s-%d_%d.json", delta, psi, sim_info$simid, sim_info$region2_uniqid)
      )
      r1_stage1_out <- jsonlite::fromJSON(stage1_r1_file)
      r1_stage1_out$coords <- coords[[sim_info$region1_uniqid]]
      r2_stage1_out <- jsonlite::fromJSON(stage1_r2_file)
      r2_stage1_out$coords <- coords[[sim_info$region2_uniqid]]

      theta <- stage2_out[id, c("rho", "k_eta1", "k_eta2", "tau_eta", "nugget_eta")]
      theta <- unlist(c(theta,
        sigma2_ep1 = r1_stage1_out$stage1$sigma2_ep,
        sigma2_ep2 = r2_stage1_out$stage1$sigma2_ep
      ))

      avrho <- qfuncMM::get_asymp_var_rho(theta, r1_stage1_out, r2_stage1_out, method)
      aci <- qfuncMM::get_asymp_ci_rho(theta, 0.95, asympvar_rho = avrho)
      rho_true <- 0.1
      if (sim_info$region1_uniqid == 1 && sim_info$region2_uniqid == 3) {
        rho_true <- 0.35
      } else if (sim_info$region1_uniqid == 2 && sim_info$region2_uniqid == 3) {
        rho_true <- 0.6
      }
      asymp_inf_mx[id, ] <- c(
        sim_info$simid, sim_info$region1_uniqid, sim_info$region2_uniqid,
        rho_true, theta["rho"], avrho, aci["lower"], aci["upper"]
      )

      write_csv(tibble::as_tibble(asymp_inf_mx), file = outfile, col_names = TRUE)
    }
    message(sprintf("\nResults saved to %s", outfile))
  }
}
