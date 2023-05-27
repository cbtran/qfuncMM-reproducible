library(qfuncMM)
library(readr)

# Compute correlation of averages given two region matrices
corr_of_avgs <- function(region1, region2) {
  region1_avg <- apply(region1, 1, mean)
  region2_avg <- apply(region2, 1, mean)
  cor(region1_avg, region2_avg)
}

# Reshape data into a matrix with n_timept rows and n_voxel columns
reshape_signal <- function(signal) {
  n_timept <- 60
  n_voxel <- 50
  return(matrix(signal, nrow = n_timept, ncol = n_voxel))
}

args <- commandArgs(trailingOnly = TRUE)
if (!(args[1] %in% c("weak", "med", "strong"))) {
  stop("First argument must be 'weak', 'med', or 'strong'")
}
if (!(args[2] %in% c("weak", "strong"))) {
  stop("Second argument must be 'weak' or 'strong'")
}
setting <- paste0(args[1], "-signal-", args[2], "-intra")
if (!is.numeric(args[3])) {
  stop("Third argument must be an integer")
}
n_sim <- as.integer(args[3])
if (n_sim <= 0 || n_sim > 100) {
  stop("Number of simulations must be between 1 and 100")
}

voxel_coords <- readRDS(
  file.path("simulation-data", setting, "voxel_coords.rds"))

sim_all <- readRDS(
  file.path("simulation-data", setting, "sim-data.rds"))[1:n_sim]

simulate_single <- function(simid) {
  region_list <- lapply(sim_all[[simid]], reshape_signal)
  ca12 <- corr_of_avgs(region_list[[1]], region_list[[2]])
  ca13 <- corr_of_avgs(region_list[[1]], region_list[[3]])
  ca23 <- corr_of_avgs(region_list[[2]], region_list[[3]])

  cat("Simulation", simid, "\n")
  result <- qfuncmm(region_list, voxel_coords, n_basis = 45)
  cat("Finished simulation", simid, "\n")

  list(result = result,
       ca = c(ca12, ca13, ca23))
}

results <- parallel::mclapply(seq_len(n_sim), simulate_single)

pair12_result <- matrix(
  nrow = n_sim, ncol = 5,
  dimnames = list(NULL, c("qfunc", "Intra", "CA", "bCA", "asympvar")))
pair13_result <- matrix(
  nrow = n_sim, ncol = 5,
  dimnames = list(NULL, c("qfunc", "Intra", "CA", "bCA", "asympvar")))
pair23_result <- matrix(
  nrow = n_sim, ncol = 5,
  dimnames = list(NULL, c("qfunc", "Intra", "CA", "bCA", "asympvar")))


for (simid in seq_along(results)) {
  r <- results[[simid]]
  pair12_result[simid, "CA"] <- r$ca[1]
  pair12_result[simid, "bCA"] <- r$result$cor_fixed_bspline[1, 2]
  pair12_result[simid, "qfunc"] <- r$result$cor[1, 2]
  pair12_result[simid, "Intra"] <- r$result$cor_fixed_intra[1, 2]
  pair12_result[simid, "asympvar"] <- r$result$asymp_var[1, 2]

  pair13_result[simid, "CA"] <- r$ca[2]
  pair13_result[simid, "bCA"] <- r$result$cor_fixed_bspline[1, 3]
  pair13_result[simid, "qfunc"] <- r$result$cor[1, 3]
  pair13_result[simid, "Intra"] <- r$result$cor_fixed_intra[1, 3]
  pair13_result[simid, "asympvar"] <- r$result$asymp_var[1, 3]

  pair23_result[simid, "CA"] <- r$ca[3]
  pair23_result[simid, "bCA"] <- r$result$cor_fixed_bspline[2, 3]
  pair23_result[simid, "qfunc"] <- r$result$cor[2, 3]
  pair23_result[simid, "Intra"] <- r$result$cor_fixed_intra[2, 3]
  pair23_result[simid, "asympvar"] <- r$result$asymp_var[2, 3]
}

saveRDS(list(r12 = pair12_result, r13 = pair13_result, r23 = pair23_result),
         file.path("simulation-output", setting, "results.rds"))
