library(qfuncMM)
here::i_am("full-run/full-run-eblue.R")
library(here)
source(here("R_files", "covariances.R"))

voxel_coords <- readRDS(here("full-run", "rat_coords.rds"))
dist_sqrd_mat_region1 <- as.matrix(dist(voxel_coords[[1]]))^2
dist_sqrd_mat_region2 <- as.matrix(dist(voxel_coords[[2]]))^2
dist_sqrd_mat_region3 <- as.matrix(dist(voxel_coords[[3]]))^2
region_dist_sqrd <- list(dist_sqrd_mat_region1,
                         dist_sqrd_mat_region2,
                         dist_sqrd_mat_region3)
num_voxel <- sapply(voxel_coords, nrow)

startid <- 1
endid <- 100

df_list <- vector("list", length = 9)
dfid <- 1

for (delta in c("low", "mid", "high")) {
  for (phi in c("low", "mid", "high")) {
    setting <- paste0(delta, "-", phi, "-M60-100-rat")
    outpath <- here("full-run", "out", paste0(setting, "-eblue.rds"))
    allsignals <- readRDS(here("full-run", "data", paste0(setting, ".rds")))
    eblues <- NULL
    if (file.exists(outpath)) {
      eblues <- readRDS(outpath)
    } else {
      num_timept <- nrow(allsignals$data[[1]][[1]])
      time_sqrd_mat <- outer(seq_len(num_timept), seq_len(num_timept), `-`)^2

      stage1_all <- readRDS(
        here("full-run", "out", paste0(setting, "-result.rds")))$stage1

      etahat <- function(regionid, simid) {
        s1 <- stage1_all[, regionid, simid]

        dist_sqrd <- region_dist_sqrd[[regionid]]
        cj <- get_cor_mat("matern_5_2", dist_sqrd, s1["phi_gamma"])
        bj <- s1["k_gamma"] * get_cor_mat("rbf", time_sqrd_mat, s1["tau_gamma"])
              + s1["nugget_gamma"] * diag(num_timept)
        vj <- kronecker(cj, bj)
        diag(vj) <- diag(vj) + 1
        vinv <- solve(vj)
        u <- kronecker(rep(1, num_voxel[regionid]), diag(num_timept))
        gt <- solve(t(u) %*% vj %*% u, t(u) %*% vinv)

        x <- allsignals$data[[simid]][[regionid]]
        dim(x) <- NULL
        result <- as.numeric(gt %*% x)
        return(result)
      }

      eblues <- matrix(nrow = endid - startid + 1, ncol = 3)
      colnames(eblues) <- c("r12", "r13", "r23")
      for (simid in startid:endid) {
        etahats <- lapply(1:3, \(rid) etahat(rid, simid))
        eblues[simid, 1] <- cor(etahats[[1]], etahats[[2]])
        eblues[simid, 2] <- cor(etahats[[1]], etahats[[3]])
        eblues[simid, 3] <- cor(etahats[[2]], etahats[[3]])
      }

      saveRDS(eblues, outpath)
      cat("Saved", outpath, "\n")
    }
    df <- tidyr::gather(data.frame(eblues), key = "pair", value = "value") |>
      dplyr::select(value, pair, everything())
    df$method <- factor("EBLUE")
    df$yintercept <- ifelse(df$pair == "r12", 0.1,
                        ifelse(df$pair == "r13", 0.35, 0.6))
    df$delta <- allsignals$setting$delta
    df$psi <- allsignals$setting$psi
    df_list[[dfid]] <- df
    dfid <- dfid + 1
  }
}

eblue_df <- Reduce(rbind, df_list)
outpath <- here("full-run", "out", "results_eblue.csv")
if (file.exists(outpath)) {
  cat(outpath, "already exists. Skipping.\n")
} else {
  cat("Writing", outpath, "\n")
  write.csv(eblue_df, outpath, row.names = FALSE, quote = FALSE)
}
