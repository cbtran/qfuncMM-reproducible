library(qfuncMM)
library(tibble)
library(readr)

# args <- commandArgs(trailingOnly = TRUE)
args <- c("out", "std", "stage2_reml")
outdir <- args[1]
dataspec <- args[2]
s2_run_str <- args[3]

s1_dir <- file.path(outdir, dataspec, "stage1")
s2_dir <- file.path(outdir, dataspec, s2_run_str)

coords <- readRDS("R_files/rat_coords.rds")
for (delta in c("high", "mid", "low")) {
  for (psi in c("high", "mid", "low")) {
    message(sprintf("Computing CA asymptotic confidence intervals for delta = %s, psi = %s", delta, psi))

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

    confidence_level <- 0.95
    alpha <- 1 - confidence_level
    z_critical <- qnorm(1 - alpha / 2)
    for (id in seq_len(nrow(stage2_out))) {
      cat(id, " ")
      sim_info <- stage2_out[id, ]
      rho_ca <- sim_info$rho_ca
      z_rho_ca <- atanh(rho_ca)
      se_z <- 1 / sqrt(60 - 3)
      z_lower <- z_rho_ca - z_critical * se_z
      z_upper <- z_rho_ca + z_critical * se_z
      rho_ca_lower <- tanh(z_lower)
      rho_ca_upper <- tanh(z_upper)
      cat(c(rho_ca, rho_ca_lower, rho_ca_upper), "\n")

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
