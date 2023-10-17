here::i_am("R_files/generate.R")
library(here)
source(here("R_files/generate_3_region.R"))

nugget_gamma <- 0.1
nugget_eta <- 0.1
k_gamma <- 2
delta <- function(x) {
  x * (1 + nugget_eta) / (x * (1 + nugget_eta) + k_gamma + nugget_gamma)
}

kEta_seq <- sapply(c(0.1, 0.5, 0.7), function(y) {
  uniroot(function(x) delta(x) - y, interval = c(0, 20))$root
})

voxel_coords <- readRDS(here("simulation-output/new-model/voxel_coords.rds"))
sqrd_dist <- lapply(voxel_coords, \(coords) as.matrix(dist(coords))^2)

psi <- function(phi, dist_sqrd) {
  mean(get_cor_mat("matern_5_2", dist_sqrd, phi))
}

phi_seq <- lapply(sqrd_dist,
                  \(dist)
                    sapply(c(0.2, 0.5, 0.8), function(y) {
                      uniroot(\(x) psi(x, dist) - y, interval = c(0, 20))$root
                    })
                  )
phi_seq <- Reduce(rbind, phi_seq)

kEta <- kEta_seq[3] # mid
phi <- phi_seq[, 2] # mid

region_parameters <- data.frame(
  k_eta = rep(kEta, 3),
  phi_gamma = phi,
  tau_gamma = rep(0.5, 3),
  k_gamma = rep(2, 3),
  nugget_gamma = rep(nugget_gamma, 3),
  mean = c(1, 10, 20),
  sigma2 = c(1, 1, 1),
  row.names = c("region1", "region2", "region3")
)
shared_parameters <- c(tau_eta = 0.25, nugget = nugget_eta)
corr_true <- c(rho12 = 0.1, rho13 = 0.35, rho23 = 0.6)
n_timept <- 60
n_sim <- 1

three_region <- generate_3_region_new(
  n_sim, voxel_coords, n_timept, corr_true, region_parameters, shared_parameters, seed = 1234)

out <- list(data = three_region,
            setting = list(region_parameters = region_parameters,
                           shared_parameters = shared_parameters,
                           corr_true = corr_true,
                           signal = "high delta = 0.7, mid phi"))

saveRDS(out, here("simulation-output/new-model/high-mid-M60.rds"))
cat("Saved to", here("simulation-output/new-model/high-mid-M60.rds"), "\n")
