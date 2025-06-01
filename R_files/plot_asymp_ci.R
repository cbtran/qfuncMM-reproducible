library(tibble)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggthemes)

results_reml_dir <- file.path("out", "std", "stage2_reml")
results_vecchia_dir <- file.path("out", "std", "stage2_vecchia")

# Initialize data frame to store results
coverage_results <- tibble()
s1_dir <- file.path("out", "std", "stage1")
s2_dir_reml <- file.path("out", "std", "stage2_reml")
s2_dir_vecchia <- file.path("out", "std", "stage2_vecchia")
coords <- readRDS("R_files/rat_coords.rds")

coverage_results <- tibble()
for (delta in c("high", "mid", "low")) {
  for (psi in c("high", "mid", "low")) {
    setting_str <- paste0(delta, "-", psi)
    reml_file <- file.path(results_reml_dir, sprintf("asymp_ci_%s-%s.csv", delta, psi))
    vecchia_file <- file.path(results_vecchia_dir, sprintf("asymp_ci_%s-%s.csv", delta, psi))

    if (!file.exists(reml_file) || !file.exists(vecchia_file)) {
      next
    }

    reml_data <- read_csv(reml_file, col_types = "iii", show_col_types = FALSE)
    vecchia_data <- read_csv(vecchia_file, col_types = "iii", show_col_types = FALSE)
    stage2_out_reml <- read_csv(file.path(s2_dir_reml, sprintf("results_%s.csv", setting_str)), col_types = "iii")
    stage2_out_vecchia <- read_csv(file.path(s2_dir_vecchia, sprintf("results_%s.csv", setting_str)), col_types = "iii")

    nsim <- nrow(stage2_out_reml)
    levels <- seq(0.75, 0.95, by = 0.05)
    coverage_reml <- matrix(nrow = nsim, ncol = length(levels), dimnames = list(NULL, levels))
    coverage_vecchia <- matrix(nrow = nsim, ncol = length(levels), dimnames = list(NULL, levels))
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
      ci_reml <- sapply(levels, \(level) {
        qfuncMM::get_asymp_ci_rho(
          theta_reml, level, as.numeric(reml_data[id, "asymp_var_rho"])
        )
      })
      ci_vecchia <- sapply(levels, \(level) {
        qfuncMM::get_asymp_ci_rho(
          theta_vecchia, level, as.numeric(vecchia_data[id, "asymp_var_rho"])
        )
      })
      coverage_reml[id, ] <- apply(ci_reml, 2, \(x) x[1] <= rho_true && x[2] >= rho_true)
      coverage_vecchia[id, ] <- apply(ci_vecchia, 2, \(x) x[1] <= rho_true && x[2] >= rho_true)
    }

    # print(sprintf("Delta: %s, Psi: %s", delta, psi))
    # print("REML Coverage:")
    # print(colMeans(coverage_reml))
    # print("Vecchia Coverage:")
    # print(colMeans(coverage_vecchia))

    coverage_results <- bind_rows(
      coverage_results,
      tibble(
        delta = as.factor(delta),
        psi = as.factor(psi),
        method = as.factor("reml"),
        coverage = colMeans(coverage_reml),
        level = as.factor(levels)
      ),
      tibble(
        delta = as.factor(delta),
        psi = as.factor(psi),
        method = as.factor("vecchia"),
        coverage = colMeans(coverage_vecchia),
        level = as.factor(levels)
      )
    )
  }
}

coverage_results$delta <- forcats::fct_relevel(coverage_results$delta, "low", "mid", "high")
coverage_results$psi <- forcats::fct_relevel(coverage_results$psi, "low", "mid", "high")

p <- coverage_results |>
  ggplot(aes(
    x = as.numeric(as.character(level)), y = coverage,
    color = method, group = method
  )) +
  geom_line() +
  geom_point(size = 2) +
  facet_grid(delta ~ psi,
    labeller = labeller(
      delta = function(x) paste0("\u03B4: ", x),
      psi = function(x) paste0("\u03C8: ", x)
    )
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
  labs(
    x = "1 - \u03B1",
    y = "Coverage"
  ) +
  theme_few() +
  scale_color_brewer(palette = "Set2") +
  theme(legend.position = "bottom", legend.title = element_blank())
p
