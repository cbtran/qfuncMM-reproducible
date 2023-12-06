library(ggplot2)
library(dplyr)
library(patchwork)

fullrun_boxplot <- function(delta, psi, covar_setting) {
  stopifnot(covar_setting %in% c("std", "ar2", "anisotropic"))
  n_timept <- 60
  covar_prefix <- ""
  if (covar_setting != "std")
    covar_prefix <- paste0("-", covar_setting)
  setting <- paste0(delta, "-", psi, "-M", n_timept, "-100-rat", covar_prefix)

  resultpath <- paste0("full-run/", setting, "-result.rds")
  if (!file.exists(resultpath)) {
    cat("Setting", paste0(delta, "-", psi), "does not exist yet. Skipping.")
    empty_plot <- ggplot() + theme_minimal()
    return(empty_plot)
  }
  result <- readRDS(resultpath)
  stage2 <- result$stage2
  stage2rho <- stage2[1, , ]
  nsim <- dim(stage2)[3]

  # Convert matrix to data frame
  df <- data.frame(t(stage2rho))
  # Reshape data to long format
  df_long <- tidyr::gather(df, key = "pair", value = "value")
  df_long$method <- "ReML"

  nsim <- dim(result$stage2)[3]
  data <- readRDS(paste0("full-run/", setting, ".rds"))
  resultCA <- matrix(nrow = nsim, ncol = 3)
  computeCA <- function(r1avg, r2avg) {
    r1avgavg <- r1avg - mean(r1avg)
    r2avgavg <- r2avg - mean(r2avg)
    sum(r1avgavg * r2avgavg) / (sd(r1avg) * sd(r2avg) * length(r1avg))
  }

  for (i in seq_len(nsim)) {
    signal <- data$data[[i]]
    ca <- lapply(signal, \(regmat) apply(regmat, 1, mean))
    resultCA[i, 1] <- computeCA(ca$region1, ca$region2)
    resultCA[i, 2] <- computeCA(ca$region1, ca$region3)
    resultCA[i, 3] <- computeCA(ca$region2, ca$region3)
  }
  dfca <- data.frame(
    value = as.numeric(resultCA),
    pair = rep(c("r12", "r13", "r23"), each = nsim))
  dfca$method <- "CA"

  ggdf <- rbind(dfca, df_long)
  ggdf$yintercept <- ifelse(ggdf$pair == "r12", 0.1,
                          ifelse(ggdf$pair == "r13", 0.35, 0.6))

  param <- expression(rho)
  # Create boxplot
  p <- ggdf |>
    ggplot(aes(x = method, y = value, fill = method)) +
      geom_boxplot() +
      ylim(-1, 1) +
      labs(x = NULL, y = NULL) +
      theme_minimal() +
      guides(fill = "none") +
      facet_wrap(~pair, nrow = 1) +
      geom_hline(aes(yintercept = yintercept),
                 col = "black", lty = 2, lwd = 0.5, alpha = 0.7 )

  return(p)

}

args <- commandArgs(trailingOnly = TRUE)
covar_setting <- "std"

plots <- list()
for (delta in c("high", "mid", "low")) {
  for (psi in c("low", "mid", "high")) {
    plots <- c(plots, list(fullrun_boxplot(delta, psi, covar_setting)))
  }
}

grid_plot <- wrap_plots(plots, ncol = 3, byrow = TRUE) +
  plot_annotation(title = "Simulations with correctly specified covariances",
                  caption = "100 replications run for each setting")
print(grid_plot)
