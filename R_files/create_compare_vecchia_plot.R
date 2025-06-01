library(ggplot2)
library(reshape2)
library(dplyr)
library(tibble)
library(ggthemes)

# Example command to run the script:
# >Rscript create_plots.R out std noisy 1 FALSE plots
# args <- c("out", "std", "noisy", 1, FALSE, "plots")
args <- commandArgs(trailingOnly = TRUE)
args <- c("out", "plots")
results_dir <- args[1]
plots_dir <- args[2]

# List all CSV files in the directory that match the cov_setting
csv_pattern <- paste0("results_.*\\.csv$")
dir_path <- file.path(results_dir, file.path("std", "stage2_reml"))
csv_files <- list.files(path = dir_path, pattern = csv_pattern, full.names = TRUE)

sim_level <- sapply(csv_files, function(file) {
  filename <- basename(file)
  setting_part <- sub("results_(.*)\\.csv", "\\1", filename)
  setting_parts <- strsplit(setting_part, "-")[[1]]
  c(delta = setting_parts[1], psi = setting_parts[2])
})

dir_path <- file.path(results_dir, file.path("std", "stage2_vecchia"))
csv_files_vecchia <- list.files(path = dir_path, pattern = csv_pattern, full.names = TRUE)

sim_level_vecchia <- sapply(csv_files_vecchia, function(file) {
  filename <- basename(file)
  setting_part <- sub("results_(.*)\\.csv", "\\1", filename)
  setting_parts <- strsplit(setting_part, "-")[[1]]
  c(delta = setting_parts[1], psi = setting_parts[2])
})

# Read all CSV files and combine into one data frame
all_data <- do.call(rbind, lapply(
  seq_along(csv_files),
  \(i) {
    df <- read.csv(csv_files[i])
    df$delta <- sim_level["delta", i]
    df$psi <- sim_level["psi", i]
    df
  }
)) |>
  select(simid, region1_uniqid, region2_uniqid, delta, psi, rho, rho_ca, rho_eblue) |>
  as_tibble()

all_data_vecchia <- do.call(rbind, lapply(
  seq_along(csv_files_vecchia),
  \(i) {
    df <- read.csv(csv_files_vecchia[i])
    df$delta <- sim_level_vecchia["delta", i]
    df$psi <- sim_level_vecchia["psi", i]
    df
  }
)) |>
  rename(rho_vecchia = rho) |>
  select(simid, region1_uniqid, region2_uniqid, delta, psi, rho_vecchia) |>
  as_tibble()

all_data <- inner_join(
  all_data,
  all_data_vecchia,
  by = c("simid", "region1_uniqid", "region2_uniqid", "delta", "psi")
)

ggdf <- tidyr::pivot_longer(all_data,
  cols = c("rho", "rho_vecchia"),
  names_to = "method",
  values_to = "value"
)

ggdf$method <- factor(ggdf$method,
  levels = c("rho", "rho_vecchia"),
  labels = c("ReML", "Vecchia")
)

ggdf$pair <- factor(paste0("r", ggdf$region1_uniqid, ggdf$region2_uniqid))
ggdf$yintercept <- ifelse(ggdf$pair == "r12", 0.1, ifelse(ggdf$pair == "r13", 0.35, 0.6))
ggdf$delta <- forcats::fct_relevel(ggdf$delta, "low", "mid", "high")
ggdf$psi <- forcats::fct_relevel(ggdf$psi, "low", "mid", "high")

p <- ggplot(ggdf) +
  geom_boxplot(mapping = aes(x = pair, y = value, fill = method)) +
  geom_segment(
    aes(
      x = as.numeric(pair) - 0.5,
      xend = as.numeric(pair) + 0.5,
      y = yintercept,
      yend = yintercept
    ),
    lty = 2
  ) +
  ylim(-1, 1) +
  facet_grid(delta ~ psi,
    labeller = labeller(
      delta = function(x) paste0("\u03B4: ", x),
      psi = function(x) paste0("\u03C8: ", x)
    )
  ) +
  theme_few() +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Region pair", y = "\u03C1") +
  theme(legend.position = "bottom", legend.title = element_blank())
p

ggsave(file.path(plots_dir, "compare_vecchia.png"),
  p,
  width = 10, height = 7, dpi = 300
)
