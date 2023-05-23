library(tidyverse)
library(ggpubr)
library(corrplot)
library(Matrix)
library(xtable)
library(ggridges)

########## Read data function ########################
read_results <- function(result_folder) {
  rho_pearson <- suppressMessages(
    read_csv(file.path(result_folder, "rho_pearson.csv"),
    show_col_types = FALSE))

  rho_inter <- suppressMessages(
    read_csv(file.path(result_folder, "rho_inter.csv"),
    show_col_types = FALSE))

  n_sim <- nrow(rho_pearson)
  corrtype <- factor(rep(c("REML", "Intra", "CA", "bCA"), each = n_sim))

  pair12 <- data.frame(
    corrtype = corrtype,
    corr = c(rho_inter$rho_12,
            rho_pearson$cor_intra_12,
            rho_pearson$cor_vanilla_12,
            rho_pearson$cor_bspline_12))

  pair13 <- data.frame(
    corrtype = corrtype,
    corr = c(rho_inter$rho_13,
            rho_pearson$cor_intra_13,
            rho_pearson$cor_vanilla_13,
            rho_pearson$cor_bspline_13))

  pair23 <- data.frame(
    corrtype = corrtype,
    corr = c(rho_inter$rho_23,
            rho_pearson$cor_intra_23,
            rho_pearson$cor_vanilla_23,
            rho_pearson$cor_bspline_23))

  list(pair12 = pair12,
       pair13 = pair13,
       pair23 = pair23)
}

boxplot_3_region <- function(results, true_corr) {

  plot_region  <- function(region, true_corr_region) {
    tidyr::pivot_longer(region, cols = c("Intra", "CA", "bCA", "qfunc"),
                names_to = "corrtype", values_to = "corr") |>
      ggplot(aes(x = fct_rev(corrtype), y = corr, fill = corrtype)) +
      geom_boxplot() +
      geom_hline(yintercept = true_corr_region,
                color = "red", linetype = "dashed") +
      ylim(-1, 1) +
      coord_flip() +
      labs(y = expression(rho),
          x = NULL) +
      theme_bw()
  }

  pair12 <- plot_region(
    select(as.data.frame(results$r12), -asympvar), true_corr$rho12)
  pair13 <- plot_region(
    select(as.data.frame(results$r13), -asympvar), true_corr$rho13)
  pair23 <- plot_region(
    select(as.data.frame(results$r23), -asympvar), true_corr$rho23)

  ggarrange(pair12, pair13, pair23, legend = "none", ncol = 1)
}


hist_1_pair <- function(pairdf, title, true_corr) {
  ggplot(pairdf,
    aes(x = corr, y = fct_rev(corrtype),
        fill = corrtype, height = after_stat(density))) +
    geom_density_ridges(
        stat = "binline", bins = 20, scale = 0.95, draw_baseline = FALSE) +
    geom_vline(xintercept = true_corr, linewidth = 0.75, color = "red") +
    xlim(-1, 1) +
    labs(title = title) +
    labs(y = NULL, x = expression(rho)) +
    theme_bw(base_size = 20) +
    theme(legend.position = "none")
}


summary_3_region <- function(summary_list, true_corr) {
  pair_12 <- summary_list$pair12
  pair_13 <- summary_list$pair13
  pair_23 <- summary_list$pair23

  pair_12_summary <- cbind(
    rho.true = true_corr$rho12,
    pair_12 %>%
    group_by(corrtype) %>%
    summarise(RMSE = sqrt(mean((corr - true_corr$rho12)^2)),
              BiasAbs = abs(mean(corr) - true_corr$rho12),
              SD = sd(corr)))

  pair_13_summary <- cbind(
    rho.true = true_corr$rho13,
    pair_13 %>%
    group_by(corrtype) %>%
    summarise(RMSE = sqrt(mean((corr - true_corr$rho13)^2)),
              BiasAbs = abs(mean(corr) - true_corr$rho13),
              SD = sd(corr)))

  pair_23_summary <- cbind(
    rho.true = true_corr$rho23,
    pair_23 %>%
    group_by(corrtype) %>%
    summarise(RMSE = sqrt(mean((corr - true_corr$rho23)^2)),
              BiasAbs = abs(mean(corr) - true_corr$rho23),
              SD = sd(corr)))
  df <- rbind(pair_12_summary, pair_13_summary, pair_23_summary)
  df




}
