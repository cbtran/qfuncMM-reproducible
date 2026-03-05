library(tibble)
library(readr)
library(dplyr)
library(tidyr)

results_reml_dir <- file.path("out", "std", "stage2_reml")
results_vecchia_dir <- file.path("out", "std", "stage2_vecchia")

# Initialize data frame to store results
coverage_results <- tibble()
s1_dir <- file.path("out", "std", "stage1")
s2_dir_reml <- file.path("out", "std", "stage2_reml")
s2_dir_vecchia <- file.path("out", "std", "stage2_vecchia")
coords <- readRDS("R_files/simulation/rat_coords.rds")

get_asymp_ci_rho_ca <- function(rho_ca, level, adjustment = 1) {
  z_rho_ca <- atanh(rho_ca)
  intlen <- stats::qnorm(1 - (1 - level) / 2) / sqrt(60 - 3)
  ci <- tanh(c(z_rho_ca - intlen, z_rho_ca + intlen))
  names(ci) <- c("lower", "upper")
  ci * adjustment
}

coverage_results_list <- list()
for (delta in c("high", "mid", "low")) {
  for (psi in c("high", "mid", "low")) {
    setting_str <- paste0(delta, "-", psi)
    reml_file <- file.path(results_reml_dir, sprintf("asymp_ci_%s-%s.csv", delta, psi))
    reml_file_approx <- file.path(results_reml_dir, sprintf("asymp_ci_%s-%s_approx.csv", delta, psi))
    vecchia_file <- file.path(results_vecchia_dir, sprintf("asymp_ci_%s-%s.csv", delta, psi))
    vecchia_file_approx <- file.path(results_vecchia_dir, sprintf("asymp_ci_%s-%s_diag.csv", delta, psi))

    if (!file.exists(reml_file) || !file.exists(vecchia_file)) {
      next
    }

    reml_data <- read_csv(reml_file, col_types = "iii", show_col_types = FALSE)
    reml_data_approx <- read_csv(reml_file_approx, col_types = "iii", show_col_types = FALSE)
    vecchia_data <- read_csv(vecchia_file, col_types = "iii", show_col_types = FALSE)
    vecchia_data_approx <- read_csv(vecchia_file_approx, col_types = "iii", show_col_types = FALSE)
    stage2_out_reml <- read_csv(file.path(s2_dir_reml, sprintf("results_%s.csv", setting_str)), col_types = "iii")
    stage2_out_vecchia <- read_csv(file.path(s2_dir_vecchia, sprintf("results_%s.csv", setting_str)), col_types = "iii")

    nsim <- nrow(stage2_out_reml)
    levels <- seq(0.75, 0.95, by = 0.05)

    # Create detailed results dataframe instead of matrices
    detailed_results <- tibble()

    for (id in seq_len(nsim)) {
      sim_info <- stage2_out_reml[id, ]
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

      delta_value <- switch(delta,
        high = 0.7,
        mid = 0.5,
        low = 0.1
      )
      psi_value <- switch(psi,
        high = 0.8,
        mid = 0.5,
        low = 0.2
      )
      # delta1 <- sim_info$k_eta1 * (1 + sim_info$nugget_eta) / (sim_info$k_eta1 * (1 + sim_info$nugget_eta) + r1_stage1_out$stage1$k_gamma + r1_stage1_out$stage1$nugget_gamma)
      # delta2 <- sim_info$k_eta1 * (1 + sim_info$nugget_eta) / (sim_info$k_eta2 * (1 + sim_info$nugget_eta) + r2_stage1_out$stage1$k_gamma + r2_stage1_out$stage1$nugget_gamma)
      alpha1 <- delta_value + (1 - delta_value) * psi_value
      alpha2 <- delta_value + (1 - delta_value) * psi_value
      k_gamma <- 2
      nugget_gamma <- 0.1
      beta1 <- 1 / (nrow(r1_stage1_out$coords) * (sim_info$k_eta1 * (1 + sim_info$nugget_eta) + k_gamma + nugget_gamma))
      beta2 <- 1 / (nrow(r2_stage1_out$coords) * (sim_info$k_eta2 * (1 + sim_info$nugget_eta) + k_gamma + nugget_gamma))
      rho_adj_scale2 <- sqrt((alpha1 + beta1) * (alpha2 + beta2) / (delta_value * delta_value))

      theta_reml <- stage2_out_reml[id, c("rho", "k_eta1", "k_eta2", "tau_eta", "nugget_eta")]
      theta_reml <- unlist(c(theta_reml,
        sigma2_ep1 = r1_stage1_out$stage1$sigma2_ep,
        sigma2_ep2 = r2_stage1_out$stage1$sigma2_ep
      ))

      theta_vecchia <- stage2_out_vecchia[id, c("rho", "k_eta1", "k_eta2", "tau_eta", "nugget_eta")]
      theta_vecchia <- unlist(c(theta_vecchia,
        sigma2_ep1 = r1_stage1_out$stage1$sigma2_ep,
        sigma2_ep2 = r2_stage1_out$stage1$sigma2_ep
      ))

      rho_true <- reml_data[id, "rho_true"]

      # Calculate CIs for all levels and methods
      ci_reml <- sapply(levels, \(level) {
        asympvar_rho <- as.numeric(reml_data[id, "asymp_var_rho"])
        if (is.na(asympvar_rho)) {
          asympvar_rho <- -1
        }
        qfuncMM::get_asymp_ci_rho(
          theta_reml, level, asympvar_rho
        )
      })
      ci_reml_approx <- sapply(levels, \(level) {
        qfuncMM::get_asymp_ci_rho(
          theta_reml, level, as.numeric(reml_data_approx[id, "asymp_var_rho"])
        )
      })
      ci_vecchia <- sapply(levels, \(level) {
        qfuncMM::get_asymp_ci_rho(
          theta_vecchia, level, as.numeric(vecchia_data[id, "asymp_var_rho"])
        )
      })
      ci_vecchia_approx <- sapply(levels, \(level) {
        qfuncMM::get_asymp_ci_rho(
          theta_vecchia, level, as.numeric(vecchia_data_approx[id, "asymp_var_rho"])
        )
      })
      ci_ca <- sapply(levels, \(level) {
        get_asymp_ci_rho_ca(sim_info$rho_ca, level)
      })
      ci_ca_adjusted <- sapply(levels, \(level) {
        get_asymp_ci_rho_ca(sim_info$rho_ca, level, adjustment = rho_adj_scale2)
      })

      # Pre-calculate all data for this simulation
      sim_data <- list()
      idx <- 1

      for (i in seq_along(levels)) {
        level <- levels[i]
        rho_true_val <- as.numeric(rho_true)

        # REML
        sim_data[[idx]] <- list(
          delta = delta, psi = psi, simid = id, level = level, method = "reml",
          lower = ci_reml[1, i], upper = ci_reml[2, i], rho_true = rho_true_val,
          covered = ci_reml[1, i] <= rho_true_val && ci_reml[2, i] >= rho_true_val
        )
        idx <- idx + 1

        # REML Approx
        sim_data[[idx]] <- list(
          delta = delta, psi = psi, simid = id, level = level, method = "reml_approx",
          lower = ci_reml_approx[1, i], upper = ci_reml_approx[2, i], rho_true = rho_true_val,
          covered = ci_reml_approx[1, i] <= rho_true_val && ci_reml_approx[2, i] >= rho_true_val
        )
        idx <- idx + 1

        # Vecchia
        sim_data[[idx]] <- list(
          delta = delta, psi = psi, simid = id, level = level, method = "vecchia",
          lower = ci_vecchia[1, i], upper = ci_vecchia[2, i], rho_true = rho_true_val,
          covered = ci_vecchia[1, i] <= rho_true_val && ci_vecchia[2, i] >= rho_true_val
        )
        idx <- idx + 1

        # Vecchia Approx
        sim_data[[idx]] <- list(
          delta = delta, psi = psi, simid = id, level = level, method = "vecchia_approx",
          lower = ci_vecchia_approx[1, i], upper = ci_vecchia_approx[2, i], rho_true = rho_true_val,
          covered = ci_vecchia_approx[1, i] <= rho_true_val && ci_vecchia_approx[2, i] >= rho_true_val
        )
        idx <- idx + 1

        # CA
        sim_data[[idx]] <- list(
          delta = delta, psi = psi, simid = id, level = level, method = "ca",
          lower = ci_ca[1, i], upper = ci_ca[2, i], rho_true = rho_true_val,
          covered = ci_ca[1, i] <= rho_true_val && ci_ca[2, i] >= rho_true_val
        )
        idx <- idx + 1

        # CA Adjusted
        sim_data[[idx]] <- list(
          delta = delta, psi = psi, simid = id, level = level, method = "ca_adjusted",
          lower = ci_ca_adjusted[1, i], upper = ci_ca_adjusted[2, i], rho_true = rho_true_val,
          covered = ci_ca_adjusted[1, i] <= rho_true_val && ci_ca_adjusted[2, i] >= rho_true_val
        )
        idx <- idx + 1
      }

      # Add to main list
      if (id == 1) {
        all_sim_data <- sim_data
      } else {
        all_sim_data <- c(all_sim_data, sim_data)
      }
    }

    # Create dataframe once from all collected data
    detailed_results <- tibble(
      delta = sapply(all_sim_data, `[[`, "delta"),
      psi = sapply(all_sim_data, `[[`, "psi"),
      simid = sapply(all_sim_data, `[[`, "simid"),
      level = sapply(all_sim_data, `[[`, "level"),
      method = sapply(all_sim_data, `[[`, "method"),
      lower = sapply(all_sim_data, `[[`, "lower"),
      upper = sapply(all_sim_data, `[[`, "upper"),
      rho_true = sapply(all_sim_data, `[[`, "rho_true"),
      covered = sapply(all_sim_data, `[[`, "covered")
    )

    # Add alpha column
    setting_str <- paste0(delta, "-", psi)
    detailed_results$alpha <- switch(setting_str,
      "high-high" = 0.94,
      "high-mid" = 0.85,
      "high-low" = 0.76,
      "mid-high" = 0.9,
      "mid-mid" = 0.75,
      "mid-low" = 0.6,
      "low-high" = 0.82,
      "low-mid" = 0.55,
      "low-low" = 0.28
    )

    coverage_results_list[[setting_str]] <- detailed_results
  }
}
coverage_results_df <- bind_rows(coverage_results_list)
# write_csv(coverage_results_df, "out/std/asymp_ci_coverage_detailed.csv")
