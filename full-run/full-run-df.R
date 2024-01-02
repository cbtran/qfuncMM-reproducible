library(tibble)

fullrun_df <- function(covar_setting) {
  stopifnot(covar_setting %in% c("std", "ar2", "anisotropic", "fgn"))
  n_timept <- 60
  covar_prefix <- ""
  if (covar_setting != "std")
    covar_prefix <- paste0("-", covar_setting)
  df_list <- vector("list", length = 9)

  dfid <- 1
  for (delta in c("high", "mid", "low")) {
    for (psi in c("low", "mid", "high")) {
      setting <- paste0(delta, "-", psi, "-M", n_timept, "-100-rat", covar_prefix)

      resultpath <- paste0("full-run/out/", setting, "-result.rds")
      cat(resultpath, '\n')
      if (!file.exists(resultpath)) {
        cat("Setting", paste0(delta, "-", psi), "does not exist yet. Skipping.\n")
        next
      }
      result <- readRDS(resultpath)
      stage2 <- result$stage2
      stage2rho <- stage2[1, , ]
      nsim <- dim(stage2)[3]

      # Convert matrix to data frame
      df <- data.frame(t(stage2rho))
      # Reshape data to long format
      df_long <- tidyr::gather(df, key = "pair", value = "value")
      df_long$method <- factor("ReML", levels = c("CA", "ReML"))

      nsim <- dim(result$stage2)[3]
      data <- readRDS(paste0("full-run/data/", setting, ".rds"))
      resultCA <- matrix(nrow = nsim, ncol = 3)
      computeCA <- function(r1avg, r2avg) {
        r1avgavg <- r1avg - mean(r1avg)
        r2avgavg <- r2avg - mean(r2avg)
        sum(r1avgavg * r2avgavg) / (sd(r1avg) * sd(r2avg) * length(r1avg))
      }

      for (i in seq_len(nsim)) {
        signal <- data$data[[i]]
        ca <- lapply(signal, \(regmat) apply(regmat, 1, mean))
        resultCA[i, 1] <- computeCA(ca$region1, ca$region2)
        resultCA[i, 2] <- computeCA(ca$region1, ca$region3)
        resultCA[i, 3] <- computeCA(ca$region2, ca$region3)
      }
      dfca <- data.frame(
        value = as.numeric(resultCA),
        pair = as.factor(rep(c("r12", "r13", "r23"), each = nsim)))
      dfca$method <- factor("CA", levels = c("CA", "ReML"))

      ggdf <- rbind(dfca, df_long)
      ggdf$yintercept <- ifelse(ggdf$pair == "r12", 0.1,
                            ifelse(ggdf$pair == "r13", 0.35, 0.6))

      ggdf$delta <- data$setting$delta
      ggdf$psi <- data$setting$psi
      df_list[[dfid]] <- ggdf
      dfid <- dfid + 1
    }
  }

  as_tibble(Reduce(rbind, df_list))
}
