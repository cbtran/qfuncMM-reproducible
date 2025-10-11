library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(tibble)
library(ggthemes)

args <- c("out", "plots")
results_dir <- args[1]
plots_dir <- args[2]
methods_list <- c("ReML", "oracle")
delta_list <- c("low", "high")
psi_list <- c("mid")

# Construct the directory path
stage2_dir <- c(file.path("std", "stage2_reml"), file.path("std", "stage2_reml_oracle"))

# List all CSV files in the directory that match the cov_setting
csv_pattern <- paste0("results_.*\\.csv$")
dir_path <- file.path(results_dir, stage2_dir[1])
csv_files <- list.files(path = dir_path, pattern = csv_pattern, full.names = TRUE)

sim_level <- sapply(csv_files, function(file) {
    filename <- basename(file)
    setting_part <- sub("results_(.*)\\.csv", "\\1", filename)
    setting_parts <- strsplit(setting_part, "-")[[1]]
    c(delta = setting_parts[1], psi = setting_parts[2])
})

dir_path <- file.path(results_dir, stage2_dir[2])
csv_files_oracle <- list.files(path = dir_path, pattern = csv_pattern, full.names = TRUE)

sim_level_oracle <- sapply(csv_files_oracle, function(file) {
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

all_data_oracle <- do.call(rbind, lapply(
    seq_along(csv_files_oracle),
    \(i) {
        df <- read.csv(csv_files_oracle[i])
        df$delta <- sim_level_oracle["delta", i]
        df$psi <- sim_level_oracle["psi", i]
        df
    }
)) |>
    rename(rho_oracle = rho) |>
    select(simid, region1_uniqid, region2_uniqid, delta, psi, rho_oracle) |>
    as_tibble()

all_data <- inner_join(
    all_data,
    all_data_oracle,
    by = c("simid", "region1_uniqid", "region2_uniqid", "delta", "psi")
)

ggdf <- tidyr::pivot_longer(all_data,
    cols = starts_with("rho"),
    names_to = "method",
    values_to = "value"
)

ggdf$method <- factor(ggdf$method,
    levels = c("rho", "rho_oracle", "rho_eblue", "rho_ca"),
    labels = c("ReML", "oracle", "EBLUE", "CA")
)

ggdf$pair <- factor(paste0("r", ggdf$region1_uniqid, ggdf$region2_uniqid))
ggdf$yintercept <- ifelse(ggdf$pair == "r12", 0.1, ifelse(ggdf$pair == "r13", 0.35, 0.6))
ggdf$delta <- forcats::fct_relevel(ggdf$delta, "low", "mid", "high")
ggdf$psi <- forcats::fct_relevel(ggdf$psi, "low", "mid", "high")

# Transform to ggdf_tbl
ggdf_tbl <- ggdf |>
    group_by(delta, psi, method, pair) |>
    summarize(
        mse = mean((value - yintercept)^2),
        mse_sd = sd((value - yintercept)^2),
        mad = mean(abs(value - yintercept)),
        mad_sd = sd(abs(value - yintercept)),
    )

# Generate LaTeX table
cat(r"(\begin{tabular}{|c|c|c|ccc|})", "\n")
cat(r"(\hline)", "\n")
cat(r"(Signal & Spatial covariance & $\rho_\text{True}$ & Type & MSE & MAD \\ )", "\n")
cat(r"(\hline)", "\n")

# 1. Filter and reshape the data
wide_data <- ggdf_tbl %>%
    filter(delta %in% delta_list, psi %in% psi_list, method %in% methods_list) %>%
    mutate(
        pair = ifelse(pair == "r12", "0.1", ifelse(pair == "r13", "0.35", "0.6")),
    ) %>%
    select(pair, method, delta, mse, mse_sd, mad, mad_sd) %>%
    pivot_wider(
        id_cols = c(pair, method),
        names_from = delta,
        values_from = c(mse, mse_sd, mad, mad_sd),
        names_sep = "_"
    )

# 2. Start LaTeX table and create the multi-level header
cat("\\begin{tabular}{|c|c|cc|cc|}\n")
cat("\\hline\n")
# cat("\\multirow{2}{*}{$\\rho$} & \\multirow{2}{*}{Method} & \\multicolumn{2}{c|}{$\\psi = 0.2$} & \\multicolumn{2}{c|}{$\\psi = 0.8$} \\\\\n")
cat("\\multirow{2}{*}{$\\rho$} & \\multirow{2}{*}{Method} & \\multicolumn{2}{c|}{$\\delta = 0.1$} & \\multicolumn{2}{c|}{$\\delta = 0.7$} \\\\\n")
cat("\\cline{3-6}\n")
cat(" & & MSE & MAD & MSE & MAD \\\\\n")
cat("\\hline\n")

# 3. Helper function to format the mean (sd) strings
fmt <- function(mean, sd) sprintf("%.3f (%.3f)", mean, sd)

# 4. Loop through the reshaped data to generate table rows
unique_pairs <- unique(wide_data$pair)
for (i in seq_along(unique_pairs)) {
    current_pair <- unique_pairs[i]
    pair_rows <- wide_data %>% filter(pair == current_pair)

    # Find the minimum MSE and MAD for this pair group
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

        # MSE Low
        mse_low_str <- fmt(row$mse_low, row$mse_sd_low)
        if (row$mse_low == min_mse_low) {
            mse_low_str <- paste0("\\textbf{", mse_low_str, "}")
        }

        # MAD Low
        mad_low_str <- fmt(row$mad_low, row$mad_sd_low)
        if (row$mad_low == min_mad_low) {
            mad_low_str <- paste0("\\textbf{", mad_low_str, "}")
        }

        # MSE High
        mse_high_str <- fmt(row$mse_high, row$mse_sd_high)
        if (row$mse_high == min_mse_high) {
            mse_high_str <- paste0("\\textbf{", mse_high_str, "}")
        }

        # MAD High
        mad_high_str <- fmt(row$mad_high, row$mad_sd_high)
        if (row$mad_high == min_mad_high) {
            mad_high_str <- paste0("\\textbf{", mad_high_str, "}")
        }

        # Combine all parts into a single LaTeX table row
        cat(paste(
            pair_str, method_str, mse_low_str, mad_low_str, mse_high_str, mad_high_str,
            sep = " & "
        ), "\\\\\n")

        first_row_in_pair <- FALSE
    }

    # Add a partial line to separate pair groups
    if (i < length(unique_pairs)) {
        cat("\\cline{2-6}\n")
    }
}

# 5. Add the final line and end the table
cat("\\hline\n")
cat("\\end{tabular}\n")
