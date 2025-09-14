library(ggplot2)
library(reshape2)
library(dplyr)
library(tibble)
library(ggthemes)
library(tidyr)

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

ggthemr::ggthemr("fresh")
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
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Region pair", y = "\u03C1") +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )
p

# ggsave(file.path(plots_dir, "compare_vecchia.pdf"),
#   p,
#   width = 10, height = 7, dpi = 300, device = cairo_pdf
# )

# Generate LaTeX table for ReML and Vecchia methods
ggdf_tbl <- ggdf |>
  group_by(delta, psi, method, pair) |>
  summarize(
    mse = mean((value - yintercept)^2),
    mse_sd = sd((value - yintercept)^2),
    mad = mean(abs(value - yintercept)),
    mad_sd = sd(abs(value - yintercept)),
  )

cat("\\begin{tabular}{|c|c|cc|cc|}\n")
cat("\\hline\n")
cat("\\multirow{2}{*}{$\\rho$} & \\multirow{2}{*}{Method} & \\multicolumn{2}{c|}{$\\psi = 0.2$} & \\multicolumn{2}{c|}{$\\psi = 0.8$} \\\\ \n")
cat("\\cline{3-6}\n")
cat(" & & MSE & MAD & MSE & MAD \\\\ \n")
cat("\\hline\n")

fmt <- function(mean, sd) sprintf("%.3f (%.3f)", mean, sd)

wide_data <- ggdf_tbl %>%
  filter(delta %in% c("low", "high"), psi %in% c("mid")) %>%
  mutate(
    pair = ifelse(pair == "r12", "0.1", ifelse(pair == "r13", "0.35", "0.6")),
  ) %>%
  select(pair, method, psi, mse, mse_sd, mad, mad_sd) %>%
  pivot_wider(
    id_cols = c(pair, method),
    names_from = delta, # psi
    values_from = c(mse, mse_sd, mad, mad_sd),
    names_sep = "_"
  )

unique_pairs <- unique(wide_data$pair)
for (i in seq_along(unique_pairs)) {
  current_pair <- unique_pairs[i]
  pair_rows <- wide_data %>% filter(pair == current_pair)

  min_mse_low <- min(pair_rows$mse_low, na.rm = TRUE)
  min_mad_low <- min(pair_rows$mad_low, na.rm = TRUE)
  min_mse_high <- min(pair_rows$mse_high, na.rm = TRUE)
  min_mad_high <- min(pair_rows$mad_high, na.rm = TRUE)

  first_row_in_pair <- TRUE

  for (j in seq_len(nrow(pair_rows))) {
    row <- pair_rows[j, ]

    pair_str <- if (first_row_in_pair) {
      paste0("\\multirow{", nrow(pair_rows), "}{*}{", current_pair, "}")
    } else {
      ""
    }

    method_str <- as.character(row$method)

    mse_low_str <- fmt(row$mse_low, row$mse_sd_low)
    if (row$mse_low == min_mse_low) {
      mse_low_str <- paste0("\\textbf{", mse_low_str, "}")
    }

    mad_low_str <- fmt(row$mad_low, row$mad_sd_low)
    if (row$mad_low == min_mad_low) {
      mad_low_str <- paste0("\\textbf{", mad_low_str, "}")
    }

    mse_high_str <- fmt(row$mse_high, row$mse_sd_high)
    if (row$mse_high == min_mse_high) {
      mse_high_str <- paste0("\\textbf{", mse_high_str, "}")
    }

    mad_high_str <- fmt(row$mad_high, row$mad_sd_high)
    if (row$mad_high == min_mad_high) {
      mad_high_str <- paste0("\\textbf{", mad_high_str, "}")
    }

    cat(paste(
      pair_str, method_str, mse_low_str, mad_low_str, mse_high_str, mad_high_str,
      sep = " & "
    ), "\\\\\n")

    first_row_in_pair <- FALSE
  }

  if (i < length(unique_pairs)) {
    cat("\\cline{2-6}\n")
  }
}

cat("\\hline\n")
cat("\\end{tabular}\n")
