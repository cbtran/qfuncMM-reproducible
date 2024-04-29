here::i_am("full-run/full-run-df.R")
library(here)

fullrun_df <- function(covar_setting) {
  stopifnot(covar_setting %in% c("std", "ar2", "anisotropic", "fgn", "eblue"))
  n_timept <- 60
  covar_prefix <- ""
  if (covar_setting != "std")
    covar_prefix <- paste0("-", covar_setting)
  df_list <- vector("list", length = 9)

  dfid <- 1
  for (delta in c("high", "mid", "low")) {
    for (psi in c("low", "mid", "high")) {
      setting <- paste0(delta, "-", psi, "-M", n_timept, "-100-rat", covar_prefix)
      resultpath <- here("full-run", "out", paste0(setting, "-result.rds"))
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
      data <- readRDS(here("full-run", "data", paste0(setting, ".rds")))
      resultCA <- matrix(nrow = nsim, ncol = 3)

      for (i in seq_len(nsim)) {
        signal <- data$data[[i]]
        ca <- lapply(signal, \(regmat) apply(regmat, 1, mean))
        resultCA[i, 1] <- cor(ca$region1, ca$region2)
        resultCA[i, 2] <- cor(ca$region1, ca$region3)
        resultCA[i, 3] <- cor(ca$region2, ca$region3)
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

  return(Reduce(rbind, df_list))
}

outpath <- here("full-run", "out", "results_std.csv")
if (file.exists(outpath)) {
  cat(outpath, "already exists. Skipping.\n")
} else {
  std_df <- fullrun_df("std")
  if (is.null(std_df)) {
    cat("No results found. Skipping.")
  } else {
    cat("Writing", outpath, "\n")
    write.csv(std_df, outpath, row.names = FALSE, quote = FALSE)
  }
}

outpath <- here("full-run", "out", "results_ar2.csv")
if (file.exists(outpath)) {
  cat(outpath, "already exists. Skipping.\n")
} else {
  ar2_df <- fullrun_df("ar2")
  if (is.null(ar2_df)) {
    cat("No AR2 results found. Skipping.")
  } else {
    cat("Writing", outpath, "\n")
    write.csv(ar2_df, outpath, row.names = FALSE, quote = FALSE)
  }
}

outphat <- here("full-run", "out", "results_anisotropic.csv")
if (file.exists(outpath)) {
  cat(outpath, "already exists. Skipping.\n")
} else {
  aniso_df <- fullrun_df("anisotropic")
  if (is.null(aniso_df)) {
    cat("No anisotropic results found. Skipping.")
  } else {
    cat("Writing", outpath, "\n")
    write.csv(aniso_df, outpath, row.names = FALSE, quote = FALSE)
  }
}

outpath <- here("full-run", "out", "results_fgn.csv")
if (file.exists(outpath)) {
  cat(outpath, "already exists. Skipping.\n")
} else {
  fgn_df <- fullrun_df("fgn")
  if (is.null(fgn_df)) {
    cat("No FGN results found. Skipping.")
  } else {
    cat("Writing", outpath, "\n")
    write.csv(fgn_df, outpath, row.names = FALSE, quote = FALSE)
  }
}
