library(qfuncMM)
source("full-run/utils.R")

setting <- "mid-low-M60"
result <- readRDS(paste0("full-run/", setting, "-100-result.rds"))
stage2 <- result$stage2

na_runs <- which(is.na(stage2[1, , ]), arr.ind = TRUE)[, 2]
large_runs <- unique(which(stage2[, , ] > 10, arr.ind = TRUE)[, 3])

runid <- 100
stage1 <- result$stage1[, , runid]
stage2run <- result$stage2[, , runid]

voxel_coords <- readRDS("full-run/voxel_coords.rds")
allsignals <- readRDS(paste0("full-run/", setting, "-100.rds"))
num_timept <- nrow(allsignals$data[[1]][[1]])
num_voxel <- c(nrow(voxel_coords[[1]]), nrow(voxel_coords[[2]]), nrow(voxel_coords[[3]]))
time_sqrd_mat <- outer(seq_len(num_timept), seq_len(num_timept), `-`)^2

signal <- allsignals$data[[runid]]

region1_mx <- matrix(signal$region1, nrow = num_timept, ncol = num_voxel[1])
region2_mx <- matrix(signal$region2, nrow = num_timept, ncol = num_voxel[2])
region3_mx <- matrix(signal$region3, nrow = num_timept, ncol = num_voxel[3])
coords1 <- voxel_coords$r1
coords2 <- voxel_coords$r2
coords3 <- voxel_coords$r3

corr_avg <- computeCA(signal)
badinit <- corr_avg[3]
region23 <-
  opt_inter_new(
    c(badinit, softminus(1), softminus(1), 0, softminus(0.1)),
    region2_mx, region3_mx, coords2, coords3, time_sqrd_mat,
    stage1[, 2], stage1[, 3], kernel_dict("matern_5_2"))
