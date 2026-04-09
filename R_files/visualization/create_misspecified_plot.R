library(ggplot2)
library(reshape2)
library(dplyr)
library(tibble)
# library(ggthemes)
ggthemr::ggthemr("fresh")

# Example command to run the script:
# >Rscript create_plots.R out std noisy 1 FALSE plots
args <- commandArgs(trailingOnly = TRUE)
# args <- c("out", "nonsep", "plots")
results_dir <- args[1]
data_spec <- args[2]
plots_dir <- args[3]

# Construct the directory path
stage2_dir <- "stage2_reml"
dir_path <- file.path(results_dir, data_spec, stage2_dir)

kang_dir <- file.path(results_dir, data_spec, "stage2_kang")

# List all CSV files in the directory that match the cov_setting
csv_pattern <- paste0("results_.*\\.csv$")
csv_files <- list.files(path = dir_path, pattern = csv_pattern, full.names = TRUE)
if (length(csv_files) == 0) {
  stop(sprintf("No simulation results found in %s", dir_path))
}

sim_level <- sapply(csv_files, function(file) {
  filename <- basename(file)
  setting_part <- sub("results_(.*)\\.csv", "\\1", filename)
  setting_parts <- strsplit(setting_part, "-")[[1]]
  c(delta = setting_parts[1], psi = setting_parts[2])
})

csv_files_kang <- list.files(path = kang_dir, pattern = csv_pattern, full.names = TRUE)
sim_level_kang <- sapply(csv_files_kang, function(file) {
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
))

all_data_kang <- do.call(rbind, lapply(
  seq_along(csv_files_kang),
  \(i) {
    df <- read.csv(csv_files_kang[i])
    df$delta <- sim_level_kang["delta", i]
    df$psi <- sim_level_kang["psi", i]
    df
  }
)) |>
  dplyr::select(simid, region1_uniqid, region2_uniqid, delta, psi, rho_kang)

all_data <- dplyr::inner_join(
  all_data,
  all_data_kang,
  by = c("simid", "region1_uniqid", "region2_uniqid", "delta", "psi")
)

all_data$delta <- forcats::fct_relevel(all_data$delta, "low", "mid", "high")
all_data$psi <- forcats::fct_relevel(all_data$psi, "low", "mid", "high")

# Convert data from wide to long format for plotting
id_vars <- setdiff(names(all_data), c("rho", "rho_eblue", "rho_ca", "rho_kang"))
ggdf <- tidyr::pivot_longer(all_data,
  cols = -all_of(id_vars),
  names_to = "method",
  values_to = "value"
)
ggdf$method <- factor(ggdf$method,
  levels = c("rho", "rho_eblue", "rho_ca", "rho_kang"),
  labels = c("ReML", "EBLUE", "CA", "ACE")
)
ggdf$pair <- factor(paste0("r", ggdf$region1_uniqid, ggdf$region2_uniqid))
ggdf$yintercept <- ifelse(ggdf$pair == "r12", 0, ifelse(ggdf$pair == "r13", 0.35, 0.6))

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
  # theme_few() +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Region pair", y = "\u03C1") +
  theme(
    legend.position = "bottom", legend.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )
# ggtitle("Nonseparable kernel with exponential and polynomial marginal")
# p


dir.create(plots_dir, showWarnings = FALSE)

outfile <- file.path(plots_dir, paste0(data_spec, ".pdf"))
cairo_pdf(outfile, width = 10, height = 7, onefile = TRUE)
print(p)
dev.off()
