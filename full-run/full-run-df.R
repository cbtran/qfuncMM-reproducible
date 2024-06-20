suppressMessages(here::i_am("full-run/full-run-df.R"))
library(here)
library(dplyr)

fullrun_df <- function(covar_setting) {
  stopifnot(covar_setting %in% c("std", "ar2", "anisotropic", "fgn"))
  n_timept <- 60
  covar_prefix <- paste0("-", covar_setting)
  df_list <- vector("list", length = 9)

  dfid <- 1
  for (delta in c("high", "mid", "low")) {
    for (psi in c("low", "mid", "high")) {
      setting <- paste0(delta, "-", psi, "-M", n_timept, "-100", covar_prefix)
      resultpath <- here("full-run", "out", paste0(setting, "-result.rds"))
      if (!file.exists(resultpath)) {
        # cat("Setting", paste0(delta, "-", psi), "does not exist yet. Skipping.\n")
        next
      }
      result <- readRDS(resultpath)
      df <- result$rho |>
        reshape2::melt(varnames = c("method", "pair", "sim"), value.name = "value") |>
        select(-sim) |>
        mutate(
          method =
            forcats::fct_recode(method,
                                "ReML" = "rho", "EBLUE" = "rho_eblue", "CA" = "rho_ca"),
          pair = as.factor(pair),
          yintercept = ifelse(pair == "r12", 0.1, ifelse(pair == "r13", 0.35, 0.6)),
          delta = factor(delta, levels = c("high", "mid", "low")),
          psi = factor(psi, levels = c("low", "mid", "high"))
        ) |>
        relocate(value, .before = method)
      df_list[[dfid]] <- df
      dfid <- dfid + 1
    }
  }

  return(Reduce(rbind, df_list))
}

outpath <- here("full-run", "out", "results_std.csv")
if (file.exists(outpath)) {
  message(outpath, " already exists. Skipping.")
} else {
  std_df <- fullrun_df("std")
  if (is.null(std_df)) {
    message("No results found. Skipping.")
  } else {
    message("Writing to ", outpath)
    write.csv(std_df, outpath, row.names = FALSE, quote = FALSE)
  }
}

outpath <- here("full-run", "out", "results_ar2.csv")
if (file.exists(outpath)) {
  message(outpath, " already exists. Skipping.")
} else {
  ar2_df <- fullrun_df("ar2")
  if (is.null(ar2_df)) {
    message("No AR2 results found. Skipping.")
  } else {
    message("Writing to ", outpath)
    write.csv(ar2_df, outpath, row.names = FALSE, quote = FALSE)
  }
}

outphat <- here("full-run", "out", "results_anisotropic.csv")
if (file.exists(outpath)) {
  message(outpath, " already exists. Skipping.")
} else {
  aniso_df <- fullrun_df("anisotropic")
  if (is.null(aniso_df)) {
    message("No anisotropic results found. Skipping.")
  } else {
    message("Writing to ", outpath)
    write.csv(aniso_df, outpath, row.names = FALSE, quote = FALSE)
  }
}

outpath <- here("full-run", "out", "results_fgn.csv")
if (file.exists(outpath)) {
  message(outpath, " already exists. Skipping.")
} else {
  fgn_df <- fullrun_df("fgn")
  if (is.null(fgn_df)) {
    message("No fGn results found. Skipping.")
  } else {
    message("Writing to ", outpath)
    write.csv(fgn_df, outpath, row.names = FALSE, quote = FALSE)
  }
}
