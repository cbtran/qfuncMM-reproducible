# Run this script in the terminal as
# >Rscript R_files/generate.R <data_spec> <cov_setting> <data_dir> <out_dir> <seed>
# <data_spec> is one of "std", "fgn", "ar2", "anisotropic", "std-noise-1e-2", "std-noise-1e-3", "std-noise-1e-7",
# <cov_setting> is one of "noisy", "noiseless", "auto"
# <data_dir> is the directory containing the data files
# <out_dir> is the output directory
# <seed> is an integer seed for random number generation.

args <- commandArgs(trailingOnly = TRUE)
data_spec <- args[1]
cov_setting <- args[2]
data_dir <- args[3]
out_dir <- args[4]
seed <- as.numeric(args[5])

set.seed(seed)

if (!dir.exists(data_dir)) {
  stop(paste("Data directory", data_dir, "does not exist."))
}
if (!dir.exists(out_dir)) {
  message(paste("Creating output directory:", out_dir))
  dir.create(out_dir, recursive = TRUE)
} else {
  if (file.access(out_dir, 2) != 0) {
    stop(paste("Output directory", out_dir, "is not writable."))
  }
}
# Construct the full path for the specific data specification
data_spec_dir <- file.path(data_dir, data_spec)
out_spec_dir <- file.path(out_dir, data_spec, "stage1")
if (cov_setting == "noiseless") {
  out_spec_dir <- file.path(out_dir, paste0(data_spec, "-noiseless"), "stage1")
}

# Check if the directory exists
if (!dir.exists(data_spec_dir)) {
  stop(paste0(
    "Data specification '", data_spec, "' does not exist in directory '", normalizePath(data_dir),
    "'. Generate the data first."
  ))
}

dir.create(out_spec_dir, recursive = TRUE, showWarnings = FALSE)

voxel_coords <- readRDS(file.path("R_files", "rat_coords.rds"))

for (delta in c("high", "mid", "low")) {
  for (psi in c("high", "mid", "low")) {
    setting_str <- paste0(delta, "-", psi)
    data_setting <- readRDS(file.path(data_spec_dir, paste0(setting_str, ".rds")))

    for (simid in seq_along(data_setting$data)) {
      d <- data_setting$data[[simid]]
      sim_name <- paste0(setting_str, "-", simid)

      for (regid in 1:3) {
        region_data <- d[[regid]]
        region_coords <- voxel_coords[[regid]]
        qfuncMM::qfuncMM_stage1(sim_name, regid, sprintf("simulated region - sim %d, reg %d", simid, regid),
          region_data, region_coords,
          cov_setting = cov_setting,
          out_dir = out_spec_dir,
          save_data_and_coords = FALSE,
          overwrite = TRUE
        )
      }
    }
  }
}
