library(GpGp, lib.loc = "~/.R-gpgp-dev/library")
library(qfuncMM)
RhpcBLASctl::omp_set_num_threads(10)
RhpcBLASctl::blas_set_num_threads(1)

gpgp_opt_run <- function(region_list, voxel_coords,
                         stage1_cov_setting = "noisy",
                         num_init = 5L,
                         verbose = FALSE) {
  kernel_type_id <- 3L
  n_region <- length(region_list)
  n_timept <- nrow(region_list[[1]])
  n_voxel <- vector(length = n_region, mode = "integer")
  for (i in seq_along(region_list)) {
    if (!is.matrix(region_list[[i]])) {
      stop(sprintf("Region %d is not a matrix.", i))
    }
    if (!is.matrix(voxel_coords[[i]])) {
      stop(sprintf("Region coordinates %d is not a matrix.", i))
    }
    if (!is.numeric(voxel_coords[[i]]) || ncol(voxel_coords[[i]]) != 3) {
      stop(sprintf("Region %d: voxels must have three numeric coordinates.", i))
    }
    if (nrow(region_list[[i]]) != n_timept) {
      stop(sprintf("Region %d: inconsistent number of time points (rows)", i))
    }
    n_voxel[i] <- ncol(region_list[[i]])
    if (n_voxel[i] != nrow(voxel_coords[[i]])) {
      stop(sprintf("Region %d: Inconsistent number of voxels (columns)", i))
    }
  }

  message("Running QFunCMM with ", n_region, " regions and ", n_timept, " time points.")
  message("Stage 1: estimating intra-regional parameters...")

  # Standardize the data matrices
  # TODO: Should we keep the raw data matrices around?
  region_list_std <- lapply(region_list, \(reg) (reg - mean(reg)) / stats::sd(reg))
  stage1_info <- vector("list", length = n_region)

  for (regid in seq_along(region_list_std)) {
    inits <- qfuncMM:::stage1_init(region_list_std[[regid]], voxel_coords[[regid]], num_init, FALSE)
    intra_out <- qfuncMM:::fit_intra_model(
      region_list_std[[regid]], voxel_coords[[regid]], inits, kernel_type_id, stage1_cov_setting, verbose
    )

    stage1_region_info <- list()
    stage1_region_info$stage1 <- intra_out$intra_param
    stage1_region_info$eblue <- intra_out$eblue
    stage1_region_info$data_std <- region_list_std[[regid]]
    stage1_region_info$coords <- voxel_coords[[regid]]
    stage1_region_info$cov_setting <- stage1_cov_setting
    stage1_info[[regid]] <- stage1_region_info
  }

  message("Finished stage 1.\nStage 2: estimating inter-regional correlations...")

  region_dimnames <- list(
    paste0("r", seq_len(n_region)),
    paste0("r", seq_len(n_region))
  )
  rho <- matrix(1,
    nrow = n_region, ncol = n_region,
    dimnames = list(
      paste0("r", seq_len(n_region)),
      paste0("r", seq_len(n_region))
    )
  )
  stage2_inter <- array(
    dim = c(n_region, n_region, 4),
    dimnames = c(
      region_dimnames,
      list(c("k_eta1", "k_eta2", "tau_eta", "nugget_eta"))
    )
  )
  rho_eblue <- matrix(1,
    nrow = n_region, ncol = n_region,
    dimnames = region_dimnames
  )
  rho_ca <- matrix(1,
    nrow = n_region, ncol = n_region,
    dimnames = region_dimnames
  )

  # Run stage 2 for each pair of regions
  for (reg1 in seq_len(n_region - 1)) {
    for (reg2 in seq(reg1 + 1, n_region)) {
      region1_info <- stage1_info[[reg1]]
      region2_info <- stage1_info[[reg2]]
      eblue_r12 <- stats::cor(region1_info$eblue, region2_info$eblue)
      rho_eblue[reg1, reg2] <- eblue_r12
      rho_eblue[reg2, reg1] <- eblue_r12

      ca <- cor(rowMeans(region1_info$data), rowMeans(region2_info$data))
      rho_ca[reg1, reg2] <- ca
      rho_ca[reg2, reg1] <- ca

      start_parms <- c(eblue_r12, 0.5, 0.5, 0.5, 0.1)
      # start_parms <- c("rho" =  eblue_r12, start_parms_s2[paste0("r", reg1, reg2), ])
      # start_parms <- unlist(start_parms)
      # stage1_parms <- rbind(
      #     unlist(region1_info$stage1[c("k_gamma", "nugget_gamma", "tau_gamma", "phi_gamma", "sigma2_ep")]),
      #     unlist(region2_info$stage1[c("k_gamma", "nugget_gamma", "tau_gamma", "phi_gamma", "sigma2_ep")])
      # )
      # stage1_parms[, 5] <- sqrt(stage1_parms[, 5])
      qfit <- fit_qfuncmm(region1_info, region2_info, start_parms = start_parms, st_scale = c(10, 1), m_seq = 60)

      theta <- qfit$covparms

      rho[reg1, reg2] <- theta["rho"]
      rho[reg2, reg1] <- theta["rho"]
      stage2_inter[reg1, reg2, ] <- theta[-1]
      stage2_inter[reg2, reg1, ] <- stage2_inter[reg1, reg2, ]
      message("Finished region pair ", reg1, " - ", reg2, "\n")
    }
  }

  message("Finished stage 2.")
  stage1_regional <- do.call(rbind, lapply(stage1_info, \(x) x$stage1))
  list(rho = rho, rho_eblue = rho_eblue, rho_ca = rho_ca, stage1 = stage1_regional, stage2 = stage2_inter)
}

gpgp_run_m_seq <- function(region_list, voxel_coords, m_seq,
                           stage1_cov_setting = "noisy",
                           num_init = 5L,
                           verbose = FALSE) {
  kernel_type_id <- 3L
  n_region <- length(region_list)
  n_timept <- nrow(region_list[[1]])
  n_voxel <- vector(length = n_region, mode = "integer")
  for (i in seq_along(region_list)) {
    if (!is.matrix(region_list[[i]])) {
      stop(sprintf("Region %d is not a matrix.", i))
    }
    if (!is.matrix(voxel_coords[[i]])) {
      stop(sprintf("Region coordinates %d is not a matrix.", i))
    }
    if (!is.numeric(voxel_coords[[i]]) || ncol(voxel_coords[[i]]) != 3) {
      stop(sprintf("Region %d: voxels must have three numeric coordinates.", i))
    }
    if (nrow(region_list[[i]]) != n_timept) {
      stop(sprintf("Region %d: inconsistent number of time points (rows)", i))
    }
    n_voxel[i] <- ncol(region_list[[i]])
    if (n_voxel[i] != nrow(voxel_coords[[i]])) {
      stop(sprintf("Region %d: Inconsistent number of voxels (columns)", i))
    }
  }

  message("Running QFunCMM with ", n_region, " regions and ", n_timept, " time points.")
  message("Stage 1: estimating intra-regional parameters...")

  # Standardize the data matrices
  # TODO: Should we keep the raw data matrices around?
  region_list_std <- lapply(region_list, \(reg) (reg - mean(reg)) / stats::sd(reg))
  stage1_info <- vector("list", length = n_region)

  for (regid in seq_along(region_list_std)) {
    inits <- qfuncMM:::stage1_init(region_list_std[[regid]], voxel_coords[[regid]], num_init, FALSE)
    intra_out <- qfuncMM:::fit_intra_model(
      region_list_std[[regid]], voxel_coords[[regid]], inits, kernel_type_id, stage1_cov_setting, verbose
    )

    stage1_region_info <- list()
    stage1_region_info$stage1 <- intra_out$intra_param
    stage1_region_info$eblue <- intra_out$eblue
    stage1_region_info$data_std <- region_list_std[[regid]]
    stage1_region_info$coords <- voxel_coords[[regid]]
    stage1_region_info$cov_setting <- stage1_cov_setting
    stage1_info[[regid]] <- stage1_region_info
  }

  message("Finished stage 1.\nStage 2: estimating inter-regional correlations...")

  m_len <- length(m_seq)
  region_dimnames <- list(
    paste0("r", seq_len(n_region)),
    paste0("r", seq_len(n_region))
  )
  rho <- array(1,
    dim = c(n_region, n_region, m_len),
    dimnames = c(
      region_dimnames,
      list(paste0("m", m_seq))
    )
  )
  stage2_inter <- array(
    dim = c(n_region, n_region, 4, m_len),
    dimnames = c(
      region_dimnames,
      list(c("k_eta1", "k_eta2", "tau_eta", "nugget_eta"), paste0("m", m_seq))
    )
  )
  rho_eblue <- matrix(1,
    nrow = n_region, ncol = n_region,
    dimnames = region_dimnames
  )
  rho_ca <- matrix(1,
    nrow = n_region, ncol = n_region,
    dimnames = region_dimnames
  )

  # Run stage 2 for each pair of regions
  for (reg1 in seq_len(n_region - 1)) {
    for (reg2 in seq(reg1 + 1, n_region)) {
      region1_info <- stage1_info[[reg1]]
      region2_info <- stage1_info[[reg2]]
      eblue_r12 <- stats::cor(region1_info$eblue, region2_info$eblue)
      rho_eblue[reg1, reg2] <- eblue_r12
      rho_eblue[reg2, reg1] <- eblue_r12

      ca <- cor(rowMeans(region1_info$data), rowMeans(region2_info$data))
      rho_ca[reg1, reg2] <- ca
      rho_ca[reg2, reg1] <- ca

      start_parms <- c(eblue_r12, 0.5, 0.5, 0.5, 0.01)
      stage1_parms <- rbind(
          unlist(region1_info$stage1[c("k_gamma", "nugget_gamma", "tau_gamma", "phi_gamma", "sigma2_ep")]),
          unlist(region2_info$stage1[c("k_gamma", "nugget_gamma", "tau_gamma", "phi_gamma", "sigma2_ep")])
      )
      stage1_parms[, 5] <- sqrt(stage1_parms[, 5])

      for (mi in seq_along(m_seq)) {
        mnbh <- m_seq[mi]
        qfit <- fit_qfuncmm(region1_info$data_std, region2_info$data_std,
          region1_info$coords, region2_info$coords, start_parms, stage1_parms,
          st_scale = c(10, 1), m_seq = mnbh
          )

        theta <- qfit$covparms

        rho[reg1, reg2, mi] <- theta["rho"]
        rho[reg2, reg1, mi] <- theta["rho"]
        stage2_inter[reg1, reg2, , mi] <- theta[-1]
        stage2_inter[reg2, reg1, , mi] <- theta[-1]
      }
      message("Finished region pair ", reg1, " - ", reg2, "\n")
    }
  }

  message("Finished stage 2.")
  stage1_regional <- do.call(rbind, lapply(stage1_info, \(x) x$stage1))
  list(rho = rho)
}

# Run this script in the terminal as
# >Rscript full-run/full-run-gpgp.R <delta> <psi> <spec> <cov_setting> <startid> <endid>
# where <delta> and <psi> is one of "high", "mid", "low",
# and <spec> is one of "std", "fgn", "ar2", "anisotropic",
# and <cov_setting> is one of "standard", "diag_time", "noiseless", "noiseless_profiled",
# and <startid> and <endid> are between 1 and 100.

args <- commandArgs(trailingOnly = TRUE)
delta <- args[1]
psi <- args[2]
spec <- args[3]
cov_setting <- args[4]
startid <- args[5]
endid <- args[6]
dataids <- seq(startid, endid)
nsim <- length(dataids)

setting <- paste0(delta, "-", psi, "-M60-100-", spec)
voxel_coords <- readRDS("full-run/rat_coords.rds")

datapath <- paste0("full-run/data/", setting, ".rds")
if (!file.exists(datapath)) {
  stop(sprintf("%s not found. Generate the data first.", datapath))
}
all_data <- readRDS(datapath)
# ss <- all_data$setting
# start_parms_true <- cbind(ss$region_parameters[, c(1, 1)], ss$shared_parameters[1], ss$shared_parameters[2])
# colnames(start_parms_true) <- c("k_eta1", "k_eta2", "tau_eta", "nugget_eta")
# rownames(start_parms_true) <- c("r12", "r13", "r23")
# start_parms_true <- start_parms_true

if (!dir.exists("full-run/out/gpgp")) {
  dir.create("full-run/out/gpgp")
}
outpath <- sprintf("full-run/out/gpgp/%s-gpgp-earlystop60.rds", setting)

message(sprintf(
  "Running spec=%s, delta=%s, psi=%s, %d simulations\n",
  spec, delta, psi, nsim
))

run <- function(signal, runid) {
  # Stage 1 param list: phi_gamma, tau_gamma, k_gamma, nugget_gamma, var_noise
  result <- gpgp_opt_run(signal, voxel_coords)
  stage2 <- result$stage2 |>
    apply(3, function(x) c(x[1, 2], x[1, 3], x[2, 3])) |>
    t()
  colnames(stage2) <- c("r12", "r13", "r23")

  rho_arr <- lapply(
    list(result$rho, result$rho_eblue, result$rho_ca),
    function(x) c(x[1, 2], x[1, 3], x[2, 3])
  )
  rho_arr <- Reduce(rbind, rho_arr)
  rownames(rho_arr) <- names(result)[1:3]
  colnames(rho_arr) <- c("r12", "r13", "r23")

  return(list(
    rho = rho_arr,
    stage1 = t(result$stage1[, -5]),
    stage2 = stage2
  ))
}

run_vary_m <- function(signal, runid, m_seq) {
  # Stage 1 param list: phi_gamma, tau_gamma, k_gamma, nugget_gamma, var_noise
  result <- gpgp_run_m_seq(signal, voxel_coords, m_seq)

  rho_arr <- matrix(nrow = length(m_seq), ncol = 3)
  for (i in seq_along(m_seq)) {
    rhoi <- result$rho[, , i]
    rho_arr[i, ] <- c(rhoi[1, 2], rhoi[1, 3], rhoi[2, 3])
  }
  colnames(rho_arr) <- c("r12", "r13", "r23")
  rownames(rho_arr) <- paste0("m", m_seq)

  return(list(
    rho = rho_arr
  ))
}

results_rho <- array(dim = c(3, 3, nsim), dimnames = list(
  c("rho", "rho_eblue", "rho_ca"),
  c("r12", "r13", "r23"), NULL
))
results_stage2 <- array(dim = c(4, 3, nsim))
dimnames(results_stage2) <- list(
  c("k_eta1", "k_eta2", "tau_eta", "nugget_eta"),
  c("r12", "r13", "r23"), NULL
)
results_stage1 <- array(dim = c(5, 3, nsim))
dimnames(results_stage1) <- list(
  c("phi_gamma", "tau_gamma", "k_gamma", "nugget_gamma", "var_noise"),
  c("r1", "r2", "r3"), NULL
)

for (i in dataids) {
  d <- all_data$data[[i]]
  run_result <- run(d, i)
  runid <- i - dataids[1] + 1
  results_rho[, , runid] <- run_result$rho
  results_stage1[, , runid] <- run_result$stage1
  results_stage2[, , runid] <- run_result$stage2
  saveRDS(
    list(
      rho = results_rho[, , 1:runid],
      stage1 = results_stage1[, , 1:runid],
      stage2 = results_stage2[, , 1:runid]
    ),
    outpath
  )
  message("Finished sim ", i)
}
saveRDS(list(rho = results_rho, stage1 = results_stage1, stage2 = results_stage2), outpath)
message("Results saved to ", outpath)

# m_seq <- seq(10, 100, by = 10)
# results_rho <- array(dim = c(length(m_seq), 3, nsim), dimnames = list(
#   paste0("m", m_seq),
#   c("r12", "r13", "r23"), NULL
# ))
# for (i in dataids) {
#   d <- all_data$data[[i]]
#   run_result <- run_vary_m(d, i, m_seq)
#   runid <- i - dataids[1] + 1
#   results_rho[, , runid] <- run_result$rho
#   saveRDS(
#     list(
#       rho = results_rho[, , 1:runid]
#     ),
#     outpath
#   )
#   message("Finished sim ", i)
# }
# saveRDS(list(rho = results_rho), outpath)
# message("Results saved to ", outpath)
