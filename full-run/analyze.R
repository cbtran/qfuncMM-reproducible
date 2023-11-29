library(ggplot2)
library(dplyr)
library(patchwork)

n_timept <- 60
setting <- paste0("high-high-M", n_timept, "-100")

result <- readRDS(paste0("full-run/", setting, "-rat-result.rds"))
stage2 <- result$stage2
stage2rho <- stage2[1, , ]
nsim <- dim(stage2)[3]

# Convert matrix to data frame
df <- data.frame(t(stage2rho))
# Reshape data to long format
df_long <- tidyr::gather(df, key = "pair", value = "value")
df_long$method <- "ReML"

nsim <- dim(result$stage2)[3]
data <- readRDS(paste0("full-run/", setting, "-rat.rds"))
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

param <- expression(rho)
# Create boxplot
p12 <- ggdf |> filter(pair == "r12") |>
  ggplot(aes(x = method, y = value, fill = method)) +
    geom_boxplot() +
    ylim(-1, 1) +
    labs(x = "Method", y = param, title = "Correlation of regions 1 and 2") +
    theme_bw() +
    guides(fill = "none") +
    geom_hline(yintercept = c(0.1), col = "black", lty = 2, lwd = 2)

p13 <- ggdf |> filter(pair == "r13") |>
  ggplot(aes(x = method, y = value, fill = method)) +
    geom_boxplot() +
    ylim(-1, 1) +
    labs(x = "Method", y = param, title = "Correlation of regions 1 and 3") +
    theme_bw() +
    guides(fill = "none") +
    geom_hline(yintercept = c(0.35), col = "black", lty = 2, lwd = 2)

p23 <- ggdf |> filter(pair == "r23") |>
  ggplot(aes(x = method, y = value, fill = method)) +
    geom_boxplot() +
    ylim(-1, 1) +
    labs(x = "Method", y = param, title = "Correlation of regions 2 and 3") +
    theme_bw() +
    guides(fill = "none") +
    geom_hline(yintercept = c(0.6), col = "black", lty = 2, lwd = 2)

subtitle <- paste0(
  nsim, " replications with delta = ", data$setting$delta, ", psi = ", data$setting$psi, ", M = ", n_timept)

p12 + p13 + p23 + plot_annotation(
  title = "Correlation estimates",
  subtitle = subtitle,
  caption = "Full run using plug-in stage 1 noise variance estimate.")
ggsave(paste0("full-run/", setting, "-rat-plot.png"))
