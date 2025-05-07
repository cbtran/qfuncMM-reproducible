args <- commandArgs(trailingOnly = TRUE)
# args <- c("std", "noisy", "data", "out", FALSE, 100)
# RhpcBLASctl::blas_set_num_threads(1)
# RhpcBLASctl::omp_set_num_threads(10)
data_spec <- args[1]
cov_setting <- args[2]
data_dir <- args[3]
out_dir <- args[4]
use_vecchia <- as.logical(args[5])
seed <- as.numeric(args[6])

if (!cov_setting %in% c("noiseless", "noisy")) {
  stop("cov_setting must be either 'noiseless' or 'noisy'")
}

set.seed(seed)
data_spec_dir <- file.path(data_dir, data_spec)
if (!dir.exists(data_spec_dir)) {
  stop(paste0(
    "Data specification '", data_spec, "' does not exist in directory '", data_dir,
    "'. Generate the data first."
  ))
}
stage1_dir <- file.path(out_dir, data_spec, "stage1")
if (cov_setting == "noiseless") {
  stage1_dir <- paste0(stage1_dir, "_noiseless")
}
if (!dir.exists(stage1_dir)) {
  stop(paste("Stage 1 directory", normalizePath(stage1_dir, mustWork = FALSE), "does not exist. Run stage 1 first."))
}
out_dir_spec <- file.path(out_dir, data_spec, "stage2")
if (use_vecchia) {
  out_dir_spec <- paste0(out_dir_spec, "_vecchia")
} else {
  out_dir_spec <- paste0(out_dir_spec, "_reml")
}
if (cov_setting == "noiseless") {
  out_dir_spec <- paste0(out_dir_spec, "_noiseless")
}
dir.create(out_dir_spec, recursive = TRUE, showWarnings = FALSE)

voxel_coords <- readRDS(file.path("R_files", "rat_coords.rds"))

for (delta in c("high", "mid", "low")) {
  for (psi in c("high", "mid", "low")) {
    setting_str <- paste0(delta, "-", psi)
    data_setting <- readRDS(file.path(data_spec_dir, paste0(setting_str, ".rds")))
    csv_file <- file.path(out_dir_spec, paste0("results_", setting_str, ".csv"))
    if (file.exists(csv_file)) {
      message(paste0("CSV file ", csv_file, " already exists. Skipping."))
      next
    }

    for (simid in seq_along(data_setting$data)) {
      d <- data_setting$data[[simid]]
      sim_name <- paste0(setting_str, "-", simid)
      s1_outfiles <- list.files(
        stage1_dir,
        pattern = paste0("qfuncMM_stage1_intra_region_", sim_name, "_[0-9]+\\.json$"),
        full.names = TRUE
      )

      reg_pairs <- list(c(1, 2), c(1, 3), c(2, 3))
      for (rp in reg_pairs) {
        r1_id <- rp[1]
        r2_id <- rp[2]
        r1_data <- d[[r1_id]]
        r2_data <- d[[r2_id]]
        r1_data_std <- (r1_data - mean(r1_data)) / stats::sd(r1_data)
        r2_data_std <- (r2_data - mean(r2_data)) / stats::sd(r2_data)
        r1_coords <- voxel_coords[[r1_id]]
        r2_coords <- voxel_coords[[r2_id]]

        data_and_coords <- list(
          data_std1 = r1_data_std,
          data_std2 = r2_data_std,
          coords1 = r1_coords,
          coords2 = r2_coords
        )

        result <- list()
        if (use_vecchia) {
          result <- qfuncMM::qfuncMM_stage2_vecchia(
            s1_outfiles[r1_id], s1_outfiles[r2_id],
            out_dir = NULL, data_and_coords = data_and_coords, overwrite = TRUE
          )
        } else {
          result <- qfuncMM::qfuncMM_stage2_reml(
            s1_outfiles[r1_id], s1_outfiles[r2_id],
            out_dir = NULL, data_and_coords = data_and_coords, overwrite = TRUE
          )
        }
        stage2_vec <- unlist(result$stage2, recursive = TRUE, use.names = TRUE)
        stage2_df <- as.data.frame(t(stage2_vec), check.names = FALSE)
        result_row <- cbind(
          data.frame(
            simid = simid,
            region1_uniqid = r1_id,
            region2_uniqid = r2_id,
            loglik = result$loglik,
            mu1 = result$mu[1],
            mu2 = result$mu[2],
            stringsAsFactors = FALSE
          ),
          stage2_df
        )
        first_write <- !file.exists(csv_file) || file.info(csv_file)$size == 0
        write.table(
          result_row,
          file = csv_file,
          sep = ",",
          row.names = FALSE,
          col.names = first_write, # only on first pass
          append = !first_write # after that, just tack on rows
        )
      }
    }
  }
}
