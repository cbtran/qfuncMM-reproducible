library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(tibble)
library(ggthemes)

args <- c("out", "plots")
results_dir <- args[1]
plots_dir <- args[2]

stage2_dir <- c(file.path("std-hcp", "stage2_vecchia"), file.path("std-hcp", "stage2_kang"))

csv_pattern <- paste0("results_.*\\.csv$")
dir_path <- file.path(results_dir, stage2_dir[1])
csv_files <- list.files(path = dir_path, pattern = csv_pattern, full.names = TRUE)

sim_level <- sapply(csv_files, function(file) {
  filename <- basename(file)
  setting_part <- sub("results_(.*)\\.csv", "\\1", filename)
  setting_parts <- strsplit(setting_part, "-")[[1]]
  c(delta = setting_parts[1], psi = setting_parts[2])
})

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

dir_path_kang <- file.path(results_dir, stage2_dir[2])
csv_files_kang <- list.files(path = dir_path_kang, pattern = csv_pattern, full.names = TRUE)
sim_level_kang <- sapply(csv_files_kang, function(file) {
  filename <- basename(file)
  setting_part <- sub("results_(.*)\\.csv", "\\1", filename)
  setting_parts <- strsplit(setting_part, "-")[[1]]
  c(delta = setting_parts[1], psi = setting_parts[2])
})

all_data_kang <- do.call(rbind, lapply(
  seq_along(csv_files_kang),
  \(i) {
    df <- read.csv(csv_files_kang[i])
    df$delta <- sim_level_kang["delta", i]
    df$psi <- sim_level_kang["psi", i]
    df
  }
)) |>
  select(simid, region1_uniqid, region2_uniqid, delta, psi, rho_kang) |>
  as_tibble()

all_data <- inner_join(
  all_data,
  all_data_kang,
  by = c("simid", "region1_uniqid", "region2_uniqid", "delta", "psi")
)

ggdf <- tidyr::pivot_longer(all_data,
  cols = starts_with("rho"),
  names_to = "method",
  values_to = "value"
)

ggdf$method <- factor(ggdf$method,
  levels = c("rho", "rho_eblue", "rho_ca", "rho_kang"),
  labels = c("Vecchia", "EBLUE", "CA", "Kang")
)

ggdf$pair <- factor(paste0("r", ggdf$region1_uniqid, ggdf$region2_uniqid))
ggdf$yintercept <- ifelse(ggdf$pair == "r12", 0, ifelse(ggdf$pair == "r13", 0.35, 0.6))
ggdf$delta <- forcats::fct_relevel(ggdf$delta, "low", "mid", "high")
ggdf$psi <- forcats::fct_relevel(ggdf$psi, "low", "mid", "high")

ggthemr::ggthemr("fresh")
p <- ggdf |>
  ggplot() +
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
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Region pair", y = "\u03C1") +
  theme(
    legend.position = "bottom", legend.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

cairo_pdf(file.path(plots_dir, "hcp_sim.pdf"), width = 10, height = 7, onefile = TRUE)
print(p)
dev.off()
