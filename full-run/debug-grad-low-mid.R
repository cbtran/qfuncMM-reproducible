library(qfuncMM)
source("R/kernel_dict.R")

stage1 <- c(0.2203418, 0.23705384, 0.19398692,
            0.4490571, 0.51875970, 0.44994406,
            2.1344588, 2.06274039, 2.52525712,
            0.1314496, 0.04766074, 0.09164101,
            1.0240741, 0.98675431, 1.01683276)
stage1 <- matrix(stage1, nrow = 5, ncol = 3, byrow = TRUE)

setting <- "low-mid-M60-1"
voxel_coords <- readRDS("scratch/full-run/voxel_coords.rds")
allsignals <- readRDS(paste0("scratch/full-run/", setting, ".rds"))
num_timept <- nrow(allsignals$data[[1]][[1]])
num_voxel <- nrow(voxel_coords[[1]])
time_sqrd_mat <- outer(seq_len(num_timept), seq_len(num_timept), `-`)^2

signal <- allsignals$data[[43]]

region2_mx <- matrix(signal$region2, nrow = num_timept, ncol = num_voxel)
region3_mx <- matrix(signal$region3, nrow = num_timept, ncol = num_voxel)
coords2 <- voxel_coords$r2
coords3 <- voxel_coords$r3

softminus <- function(x) {
  log(exp(x) - 1)
}

# warn_log <- file("scratch/full-run/warnings.log", open = "wt")
# sink(warn_log, append = TRUE, type = "message")

badinit <- softminus(0.04889)
region23 <- tryCatch(
  opt_inter_new(
    c(badinit, 0, 0, 0, 0),
    region2_mx, region3_mx, coords2, coords3, time_sqrd_mat,
    stage1[, 2], stage1[, 3], kernel_dict("matern_5_2")),
  error = function(e) {
    message("error in node 13")
    message(e)
  }
)

# sink(type = "message")
# close(warn_log)