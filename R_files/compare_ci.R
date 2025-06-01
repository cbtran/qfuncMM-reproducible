library(tibble)
library(readr)
library(dplyr)
library(gt)
library(tidyr)

results_reml_dir <- file.path("out", "std", "stage2_reml")
results_vecchia_dir <- file.path("out", "std", "stage2_vecchia")

# Initialize data frame to store results
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

    coverage_reml <- reml_data |>
      mutate(
        coverage = lower_95 <= rho_true & upper_95 >= rho_true,
        pair = as.factor(paste0("r", region1_uniqid, region2_uniqid)),
      ) |>
      select(pair, coverage) |>
      group_by(pair) |>
      summarize(coverage = mean(coverage, na.rm = TRUE)) |>
      mutate(
        delta = as.factor(delta), psi = as.factor(psi),
        method = as.factor("reml")
      )

    coverage_vecchia <- vecchia_data |>
      mutate(
        coverage = lower_95 <= rho_true & upper_95 >= rho_true,
        pair = as.factor(paste0("r", region1_uniqid, region2_uniqid))
      ) |>
      select(pair, coverage) |>
      group_by(pair) |>
      summarize(coverage = mean(coverage, na.rm = TRUE)) |>
      mutate(
        delta = as.factor(delta), psi = as.factor(psi),
        method = as.factor("vecchia")
      )

    # Store results
    coverage_results <- bind_rows(
      coverage_results,
      coverage_reml,
      coverage_vecchia
    )
  }
}

cat("Overall coverage across all pairs\n")
coverage_results |>
  group_by(delta, psi, method) |>
  summarize(
    coverage = mean(coverage, na.rm = TRUE),
    .groups = "drop"
  ) |>
  pivot_wider(
    names_from = method,
    values_from = coverage
  ) |>
  print()

# coverage_table <- coverage_results |>
#   arrange(delta, psi, pair, method) |>
#   gt(groupname_col = c("delta", "psi")) |>
#   cols_label(
#     pair = "Region Pair",
#     coverage = "Coverage Rate",
#     method = "Method"
#   ) |>
#   fmt_number(
#     columns = coverage,
#     decimals = 3
#   ) |>
#   tab_header(
#     title = "Confidence Interval Coverage Comparison",
#     subtitle = "REML vs Vecchia Methods Across Different Settings"
#   ) |>
#   tab_style(
#     style = cell_fill(color = "lightblue"),
#     locations = cells_body(
#       columns = everything(),
#       rows = method == "reml"
#     )
#   ) |>
#   tab_style(
#     style = cell_fill(color = "lightgreen"),
#     locations = cells_body(
#       columns = everything(),
#       rows = method == "vecchia"
#     )
#   )

# print(coverage_table)
