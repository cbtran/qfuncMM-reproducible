args <- commandArgs(trailingOnly = TRUE)

signal <- switch(args[1],
                 "0.5" = "weak",
                 "1" = "med",
                 "1.5" = "strong")
intra <- switch(args[2],
                "1" = "weak",
                "1/2" = "med",
                "1/4" = "strong")

setting <- paste0(signal, "-signal-", intra, "-intra")
source(paste0("simulation_setting/", setting, ".R"))

suppressPackageStartupMessages(library("qfuncMM"))
suppressPackageStartupMessages(library("splines"))
source("R_files/simulation_source.R")
source("R_files/summary.R")

L <- as.integer(args[3])
side_length <- as.integer(args[4])
M <- as.integer(args[5])
num_sim <- as.integer(args[6])
K <- as.integer(args[7])
kernel_type = "matern_5_2"

# Generate voxels for each region
set.seed(1)
vxlID_1 <- sample(1:(side_length^3), L, replace = FALSE)
set.seed(2)
vxlID_2 <- sample(1:(side_length^3), L, replace = FALSE)
set.seed(4)
vxlID_3 <- sample(1:(side_length^3), L, replace = FALSE)

dist_sqrd_mat_region1 <- get_dist_sqrd_mat(L, side_length, vxlID_1)
dist_sqrd_mat_region2 <- get_dist_sqrd_mat(L, side_length, vxlID_2)
dist_sqrd_mat_region3 <- get_dist_sqrd_mat(L, side_length, vxlID_3)
time_sqrd_mat <- (outer(1:M, 1:M, "-"))^2


# Use different folder name to store results from different simulations.
# Folder name.
fdrName <-  paste0("simulation-data/", setting)
dir.create(file.path("simulation-data", setting), recursive = TRUE)
# Save data
write.csv(parameters_true, paste0(fdrName, "/parameters-true.csv")) 
write.csv(rho_vec, paste0(fdrName, "/rho_true.csv")) 
write.csv(vxlID_1, paste0(fdrName, "/vxlID_1.csv")) 
write.csv(vxlID_2, paste0(fdrName, "/vxlID_2.csv")) 
write.csv(vxlID_3, paste0(fdrName, "/vxlID_3.csv")) 


# Simulate data 
simulated_data = simulate_3_region(L1=L,
                                   L2=L,
                                   L3=L,
                                   side_length=side_length,
                                   M=M,
                                   theta=parameters_true,
                                   rho_vec=rho_vec,
                                   sigma_sqrd=1,
                                   mu_1=rep(parameters_true$mu_1, M),
                                   mu_2=rep(parameters_true$mu_2, M),
                                   mu_3=rep(parameters_true$mu_3, M),
                                   vxlID_1=vxlID_1,
                                   vxlID_2=vxlID_2,
                                   vxlID_3=vxlID_3,
                                   random_seed=1,
                                   num_sim=num_sim,
                                   out_folder_name=fdrName,
                                   plot_it=T,
                                   C_kernel_type="matern_5_2")
# Create result list
result_list <- list()

# Loop through each simulation
for (xx in 1:num_sim){
  result_list[[xx]] <- opt_simulation_3region (xx)
  cat(paste0("Done sim:", xx, "\n"))
}

# Rho Pearson's
rho_pearson = matrix(unlist(sapply(result_list, get, x="rho_pearson")), byrow = T, nrow=num_sim)
colnames(rho_pearson) = colnames(result_list[[1]]$rho_pearson)

# Rho REML
rho_inter = matrix(unlist(sapply(result_list, get, x="rho_inter")), byrow = T, nrow=num_sim)
colnames(rho_inter) = colnames(result_list[[1]]$rho_inter)

# Asymptotic variance of inter-regional correlation
asymp_var_inter = matrix(unlist(sapply(result_list, get, x="asymp_var_inter")), byrow = T, nrow=num_sim)
colnames(asymp_var_inter) = colnames(result_list[[1]]$asymp_var_inter)


#Save data if necessary
output_dir <- paste0("simulation-output/", setting)
dir.create(file.path('simulation-output', setting), recursive = TRUE)
filepath <- file.path(paste0(output_dir, "/result_list.rds"))
saveRDS(result_list, file=filepath)

write.csv(rho_pearson, paste0(output_dir, "/rho_pearson.csv"))
write.csv(rho_inter, paste0(output_dir, "/rho_inter.csv"))
write.csv(asymp_var_inter, paste0(output_dir, "/asymp_var_inter.csv"))


summary_list <- read_results(output_dir)
pdf(paste0(output_dir, "/", setting, ".pdf"), width = 12, height = 12)
box_plot_3_region(summary_list, rho_vec)
invisible(dev.off())







