args <- commandArgs(trailingOnly = TRUE)

# check if args[1] is valid, throw exception if not

signal <- switch(args[1],
                 "0.5" = "weak",
                 "1" = "med",
                 "1.5" = "strong",
                 NULL)

if (is.null(signal)) {
  stop("Invalid signal argument")
}

intra <- switch(args[2],
                "1" = "weak",
                "1/4" = "strong",
                NULL)

if (is.null(intra)) {
    stop("Invalid intra argument")
}

setting <- paste0(signal, "-signal-", intra, "-intra")
source(paste0("simulation_setting/", setting, ".R"))

suppressPackageStartupMessages(library("qfuncMM"))
suppressPackageStartupMessages(library("splines"))
source("R_files/simulation_source.R")
source("R_files/summary.R")

n_voxel <- as.integer(args[3])
side_length <- as.integer(args[4])
# n_timept <- as.integer(args[5])
n_timept <- 60 # number of time points is fixed from simulation data
num_sim <- as.integer(args[6])
n_bspline <- as.integer(args[7])

# Generate voxels for each region
sampled_voxel_ids <- c()
suppressMessages({
    vxlID_1 <- read_csv(
        file.path("simulation-data", setting, "vxlID_1.csv"),
        show_col_types = FALSE, col_select = -1)
    vxlID_2 <- read_csv(
        file.path("simulation-data", setting, "vxlID_2.csv"),
        show_col_types = FALSE, col_select = -1)
    vxlID_3 <- read_csv(
        file.path("simulation-data", setting, "vxlID_3.csv"),
        show_col_types = FALSE, col_select = -1)
    sampled_voxel_ids <- cbind(vxlID_1, vxlID_2, vxlID_3)
    colnames(sampled_voxel_ids) <- c("x1", "x2", "x3")
})

# Simulate data
simulated_data <- readRDS(file.path("simulation-data", setting, "sim-data.rds"))
simulated_data <- simulated_data[seq_len(num_sim)]

simulate_single <- function(sim_id) {
    cat("Begin sim:", sim_id, "\n")
    result <- opt_simulation_3region(
        simulated_data[[sim_id]], sampled_voxel_ids,
        side_length, n_voxel, n_timept, n_bspline)
    cat("Finished sim:", sim_id, "\n")

    result
}

# TODO: make parallel optional
result_list <- parallel::mclapply(
    seq_len(num_sim), simulate_single, mc.cores = 20L)

# Rho Pearson's
rho_pearson <- as_tibble(reduce(lapply(result_list, \(x) x$rho_pearson), rbind))

# Rho REML
rho_inter <- as_tibble(reduce(lapply(result_list, \(x) x$rho_inter), rbind))

# Asymptotic variance of inter-regional correlation
asymp_var_inter <- as_tibble(
    reduce(lapply(result_list, \(x) x$asymp_var_inter), rbind))


#Save data if necessary
output_dir <- file.path("simulation-output", setting)
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}
filepath <- file.path(output_dir, "result_list.rds")
write_csv(rho_vec, file.path(output_dir, "rho_true.csv"))
saveRDS(result_list, file = filepath)

write_csv(rho_pearson, file.path(output_dir, "rho_pearson.csv"))
write_csv(rho_inter, file.path(output_dir, "rho_inter.csv"))
write_csv(asymp_var_inter, file.path(output_dir, "asymp_var_inter.csv"))


summary_list <- read_results(output_dir)
pdf(file.path(output_dir, paste0(setting, ".pdf")), width = 12, height = 12)
box_plot_3_region(summary_list, rho_vec)
invisible(dev.off())

write.csv(summary_3_region(summary_list, rho_vec),
    file.path(output_dir, paste0(setting, ".csv")), row.names = FALSE)
