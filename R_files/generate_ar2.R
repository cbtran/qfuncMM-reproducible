here::i_am("R_files/generate_ar2.R")
library(here)
source(here("R_files/generate_3_region.R"))

nugget_gamma <- 0
nugget_eta <- 0
k_gamma <- 2
delta_fn <- function(x) {
  x * (1 + nugget_eta) / (x * (1 + nugget_eta) + k_gamma + nugget_gamma)
}
delta_seq <- c(0.1, 0.5, 0.7)
names(delta_seq) <- c("low", "mid", "high")

kEta_seq <- sapply(delta_seq, function(y) {
  uniroot(function(x) delta_fn(x) - y, interval = c(0, 20))$root
})

voxel_coords <- readRDS(here("full-run/rat_coords.rds"))
sqrd_dist <- lapply(voxel_coords, \(coords) as.matrix(dist(coords))^2)

psi_fn <- function(phi, dist_sqrd) {
  mean(get_cor_mat("matern_5_2", dist_sqrd, phi))
}

psi_seq <- c(0.2, 0.5, 0.8)
names(psi_seq) <- c("low", "mid", "high")
phi_seq <- lapply(sqrd_dist,
                  \(dist)
                    sapply(c(0.2, 0.5, 0.8), function(y) {
                      uniroot(\(x) psi_fn(x, dist) - y, interval = c(0, 20))$root
                    })
                  )
phi_seq <- Reduce(rbind, phi_seq)

delta <- 0.1
psi <- 0.2

kEta <- kEta_seq[which(delta_seq == delta)]
phi <- phi_seq[, which(psi_seq == psi)] # mid

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
n_timept <- 1070 # yields 60 wavelet coeffs
n_sim <- 100

three_region <- generate_ar2_region(
  n_sim, voxel_coords, n_timept, corr_true, region_parameters, shared_parameters, seed = 1234)

out <- list(data = three_region,
            setting = list(region_parameters = region_parameters,
                           shared_parameters = shared_parameters,
                           corr_true = corr_true,
                           delta = delta,
                           psi = psi))

outsetting <- paste0(names(delta_seq)[which(delta_seq == delta)],
                     "-",
                     names(psi_seq)[which(psi_seq == psi)],
                     "-",
                     "M",
                     n_timept,
                     "-",
                     n_sim, "-ar2")
saveRDS(out, here("full-run", paste0(outsetting, ".rds")))
cat("Saved to", here("full-run", paste0(outsetting, ".rds")), "\n")
