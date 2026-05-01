# Runs the supplemental result showing robustness of Vecchia neighbor selection to the number of neighbors m.

RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(10)
args <- commandArgs(trailingOnly = TRUE)
data_dir <- args[1]
out_dir <- args[2]
data_spec <- "std"
cov_setting <- "noisy"
delta_setting <- "mid"
psi_setting <- "mid"

if (length(args) < 2) {
  stop("Usage: Rscript run_vecchia_neighbors.R <data_dir> <out_dir>")
}

data_spec_dir <- file.path(data_dir, data_spec)
if (!dir.exists(data_spec_dir)) {
  stop(paste0(
    "Data specification '", data_spec, "' does not exist in directory '",
    data_dir, "'. Generate the data first."
  ))
}

stage1_dir <- file.path(out_dir, data_spec, "stage1")
if (!dir.exists(stage1_dir)) {
  stop(paste(
    "Stage 1 directory", normalizePath(stage1_dir, mustWork = FALSE),
    "does not exist. Run stage 1 first."
  ))
}

out_dir_spec <- file.path(out_dir, data_spec, "stage2_vecchia_neighbors")
dir.create(out_dir_spec, recursive = TRUE, showWarnings = FALSE)

coords_file <- "rat_coords.rds"
voxel_coords <- readRDS(file.path("R_files", "simulation", coords_file))

m_seq_values <- c(30, 50)

delta <- delta_setting
psi <- psi_setting
setting_str <- paste0(delta, "-", psi)
data_setting <- readRDS(file.path(data_spec_dir, paste0(setting_str, ".rds")))
csv_file <- file.path(out_dir_spec, paste0("results_", setting_str, ".csv"))

if (file.exists(csv_file)) {
  message(paste0("CSV file ", csv_file, " already exists. Skipping."))
} else {
  for (simid in seq_along(data_setting$data)) {
    d <- data_setting$data[[simid]]
    sim_name <- paste0(setting_str, "-", simid)
    s1_outfiles <- list.files(
      stage1_dir,
      pattern = paste0(
        "qfuncMM_stage1_intra_region_", sim_name, "_[0-9]+\\.json$"
      ),
      full.names = TRUE
    )

    reg_pairs <- list(c(1, 2), c(1, 3), c(2, 3))
    for (rp in reg_pairs) {
      r1_id <- rp[1]
      r2_id <- rp[2]
      r1_data <- d[[r1_id]]
      r2_data <- d[[r2_id]]
      r1_sd <- stats::sd(r1_data)
      r2_sd <- stats::sd(r2_data)
      r1_data_std <- (r1_data - mean(r1_data)) / r1_sd
      r2_data_std <- (r2_data - mean(r2_data)) / r2_sd
      r1_coords <- voxel_coords[[r1_id]]
      r2_coords <- voxel_coords[[r2_id]]

      j1 <- jsonlite::read_json(s1_outfiles[r1_id], simplifyVector = TRUE)
      j2 <- jsonlite::read_json(s1_outfiles[r2_id], simplifyVector = TRUE)
      if (identical(j1$stage1$sigma2_ep, "NA")) j1$stage1$sigma2_ep <- NA
      if (identical(j2$stage1$sigma2_ep, "NA")) j2$stage1$sigma2_ep <- NA
      j1$data_std <- r1_data_std
      j1$coords <- r1_coords
      j2$data_std <- r2_data_std
      j2$coords <- r2_coords

      for (m_seq in m_seq_values) {
        t_start <- proc.time()[["elapsed"]]
        result <- qfuncMM:::run_stage2(
          j1, j2,
          out_dir = NULL, method = "vecchia",
          m_seq = m_seq, overwrite = TRUE
        )
        wall_time <- proc.time()[["elapsed"]] - t_start

        stage2_vec <- unlist(result$stage2, recursive = TRUE, use.names = TRUE)
        stage2_df <- as.data.frame(t(stage2_vec), check.names = FALSE)
        result_row <- cbind(
          data.frame(
            simid = simid,
            m_seq = m_seq,
            region1_uniqid = r1_id,
            region2_uniqid = r2_id,
            loglik = result$loglik,
            mu1 = result$mu[1],
            mu2 = result$mu[2],
            stringsAsFactors = FALSE
          ),
          stage2_df,
          data.frame(wall_time_s = wall_time)
        )
        first_write <- !file.exists(csv_file) || file.info(csv_file)$size == 0
        write.table(
          result_row,
          file = csv_file,
          sep = ",",
          row.names = FALSE,
          col.names = first_write,
          append = !first_write
        )
      }
    }
  }
}
