library(qfuncMM)
source("R/kernel_dict.R")

stage1 <- c(0.2311789, 0.2705925, 0.2526332,
            0.4449437, 0.5660230, 0.5213214,
            1.9680793, 1.6037095, 1.6141634,
            0.0937352, 0.1218195, 0.1201281,
            1.0625530, 0.9999726, 0.9993472)
stage1 <- matrix(stage1, nrow = 5, ncol = 3, byrow = TRUE)

setting <- "mid-mid-M60-1"
voxel_coords <- readRDS("scratch/full-run/voxel_coords.rds")
allsignals <- readRDS(paste0("scratch/full-run/", setting, ".rds"))
num_timept <- nrow(allsignals$data[[1]][[1]])
num_voxel <- nrow(voxel_coords[[1]])
time_sqrd_mat <- outer(seq_len(num_timept), seq_len(num_timept), `-`)^2

signal <- allsignals$data[[13]]

region1_mx <- matrix(signal$region1, nrow = num_timept, ncol = num_voxel)
region3_mx <- matrix(signal$region3, nrow = num_timept, ncol = num_voxel)
coords1 <- voxel_coords$r1
coords3 <- voxel_coords$r3

softminus <- function(x) {
  log(exp(x) - 1)
}

# warn_log <- file("scratch/full-run/warnings.log", open = "wt")
# sink(warn_log, append = TRUE, type = "message")

badinit <- 1.018142
region13 <- tryCatch(
  opt_inter_new(
    c(badinit, 0, 0, 0, 0),
    region1_mx, region3_mx, coords1, coords3, time_sqrd_mat,
    stage1[, 1], stage1[, 3], kernel_dict("matern_5_2")),
  error = function(e) {
    message("error in node 13")
    message(e)
  }
)

# sink(type = "message")
# close(warn_log)