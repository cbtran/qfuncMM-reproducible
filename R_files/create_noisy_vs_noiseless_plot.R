suppressPackageStartupMessages(library(dplyr))
ggthemr::ggthemr("fresh")
library(readr)
library(stringr)
library(ggplot2)

# Define the directories to search
dirs <- c(
  "out/std/stage2_reml",
  "out/std/stage2_reml_noiseless",
  "out/std-noise-1e-2/stage2_reml",
  "out/std-noise-1e-2/stage2_reml_noiseless",
  "out/std-noise-1e-3/stage2_reml",
  "out/std-noise-1e-3/stage2_reml_noiseless",
  "out/std-noise-1e-7/stage2_reml",
  "out/std-noise-1e-7/stage2_reml_noiseless"
)

# Helper function to extract noise level and noisy/noiseless from path
get_labels <- function(path) {
  noise_match <- str_match(path, "out/(std(?:-noise-([\\deE.-]+))?)/")
  noise_level <- ifelse(is.na(noise_match[, 3]), "1", noise_match[, 3])
  noisy <- ifelse(str_detect(path, "noiseless"), "noiseless", "noisy")
  list(noise_level = noise_level, noisy = noisy)
}

# Helper function to extract delta and psi from filename
get_delta_psi <- function(filename) {
  # Assumes filename like results_mid-high.csv
  match <- str_match(basename(filename), "results_([^-]+)-([^.]+)\\.csv")
  delta <- match[, 2]
  psi <- match[, 3]
  list(delta = delta, psi = psi)
}

# Glob all csv files in the directories
csv_files <- unlist(lapply(dirs, function(dir) {
  list.files(path = dir, pattern = "results.*\\.csv$", full.names = TRUE)
}))

# Read and label all data
noisy_vs_noiseless_df <- lapply(csv_files, function(f) {
  df <- read_csv(f, show_col_types = FALSE)
  labels <- get_labels(f)
  dp <- get_delta_psi(f)
  df$noise_level <- labels$noise_level
  df$noisy <- labels$noisy
  df$delta <- dp$delta
  df$psi <- dp$psi
  df
}) %>% bind_rows()

# write_csv(noisy_vs_noiseless_df, "out/noisy_vs_noiseless.csv")

noisy_vs_noiseless_df <- read_csv("out/noisy_vs_noiseless.csv", col_types = "iiiddddddddddffff")

noisy_vs_noiseless_df <- noisy_vs_noiseless_df |>
  mutate(
    pair = factor(paste0("r", region1_uniqid, region2_uniqid)),
    rho_true = ifelse(pair == "r12", 0.1, ifelse(pair == "r13", 0.35, 0.6))
  ) |>
  select(delta, psi, pair, rho, rho_true, noise_level, noisy)

# Compute RMSE and standard deviation of error for each noise_level and noisy
rmse_df <- noisy_vs_noiseless_df |>
  group_by(noise_level, noisy, pair, delta, psi) |>
  summarise(
    rmse = sqrt(mean((rho - rho_true)^2, na.rm = TRUE)),
    mad = mean(abs(rho - rho_true), na.rm = TRUE),
    sd_error = sd(abs(rho - rho_true), na.rm = TRUE),
    .groups = "drop"
  )

p <- rmse_df |>
  filter(delta == "mid", psi == "mid") |>
  ggplot(aes(x = noise_level, y = mad, color = noisy, linetype = pair, group = interaction(noisy, pair))) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(x = "Noise Level", y = "Mean absolute deviation") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank())

outfile <- "plots/noisy_vs_noiseless.pdf"
cairo_pdf(outfile, width = 10, height = 7, onefile = TRUE)
print(p)
dev.off()
