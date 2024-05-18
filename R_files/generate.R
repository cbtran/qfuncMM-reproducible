# Generate data for three-region regions using real voxel coordinates.
# Run this script in the terminal as
# >Rscript R_files/generate.R <delta> <psi> <spec> <nsim> <seed>
# where <delta> and <psi> is one of "high", "mid", "low",
# <spec> is one of "std", "fgn", "ar2", "anisotropic",
# <nsim> is the number of simulations to generate,
# <seed> is an integer seed for random number generation.


here::i_am("R_files/generate.R")
library(here)
source(here("R_files/generate_3_region.R"))
source(here("R_files/spatial-anisotropic.R"))

# Expect argument such as "mid mid std"
args <- commandArgs(trailingOnly = TRUE)
covar_setting <- args[3]
cat(sprintf("Generating data with %s-%s %s setting...\n", args[1], args[2], covar_setting))
stopifnot(covar_setting %in% c("std", "ar2", "fgn", "anisotropic"))
n_sim <- as.numeric(args[4])
seed <- as.numeric(args[5])

nugget_gamma <- 0.1
nugget_eta <- 0.1
if (covar_setting == "ar2") {
  nugget_gamma <- 0
  nugget_eta <- 0
}
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
phi_seq <- lapply(
  sqrd_dist,
  \(dist) {
    sapply(c(0.2, 0.5, 0.8), function(y) {
      uniroot(\(x) psi_fn(x, dist) - y, interval = c(0, 20))$root
    })
  }
)
phi_seq <- Reduce(rbind, phi_seq)

# Set delta and psi based on the input
delta <- delta_seq[args[1]]
psi <- psi_seq[args[2]]
print(paste0("delta = ", delta, ", psi = ", psi))

kEta <- kEta_seq[which(delta_seq == delta)]
phi <- phi_seq[, which(psi_seq == psi)] # mid

region_parameters <- data.frame(
  k_eta = rep(kEta, 3),
  phi_gamma = phi,
  tau_gamma = rep(0.5, 3),
  k_gamma = rep(2, 3),
  nugget_gamma = rep(nugget_gamma, 3),
  mean = c(1, 10, 20),
  var_noise = c(1, 1, 1),
  row.names = c("region1", "region2", "region3")
)
shared_parameters <- c(tau_eta = 0.25, nugget = nugget_eta)
corr_true <- c(rho12 = 0.1, rho13 = 0.35, rho23 = 0.6)

n_timept <- 60
if (covar_setting %in% c("fgn", "ar2")) {
  n_timept <- 1070
}

three_region <- switch(covar_setting,
  "std" = {
    spatial <- function(coords, phi_gamma) {
      dist_sqrd_mat <- as.matrix(dist(coords))^2
      get_cor_mat("matern_5_2", dist_sqrd_mat, phi_gamma)
    }
    generate_3_region_new(
      n_sim, voxel_coords, n_timept, corr_true,
      region_parameters, shared_parameters, spatial,
      seed = seed
    )
  },
  "ar2" = {
    temporal <- function(n, t, p, k, spatial_lower) {
      generate_ar2(n, t, p, 0.4, 0.3, k, spatial_lower)
    }
    generate_ar2_region(
      n_sim, voxel_coords, n_timept, corr_true,
      region_parameters, shared_parameters, temporal,
      seed = seed
    )
  },
  "fgn" = {
    temporal <- function(n, t, p, k, spatial_lower) {
      generate_fgn(n, t, p, 0.7, k, spatial_lower)
    }
    generate_ar2_region(
      n_sim, voxel_coords, n_timept, corr_true,
      region_parameters, shared_parameters, temporal,
      seed = seed
    )
  },
  "anisotropic" = {
    spatial <- function(coords, phi_gamma) {
      anisotropic(coords, phi_gamma, c(1, 1.2, 1.5), c(1, 1.2, 1.5))
    }
    generate_3_region_new(
      n_sim, voxel_coords, n_timept, corr_true,
      region_parameters, shared_parameters, spatial,
      seed = seed
    )
  }
)

out <- list(data = three_region,
            setting = list(region_parameters = region_parameters,
                           shared_parameters = shared_parameters,
                           corr_true = corr_true,
                           delta = delta,
                           psi = psi),
            seed = seed)

outsetting <- paste0(names(delta_seq)[which(delta_seq == delta)],
                     "-",
                     names(psi_seq)[which(psi_seq == psi)],
                     "-",
                     "M60",
                     "-",
                     n_sim, "-rat")
if (covar_setting != "std") {
  outsetting <- paste0(outsetting, "-", covar_setting)
}
outpath <- here("full-run/data", paste0(outsetting, ".rds"))
saveRDS(out, outpath)
cat("Saved to", outpath, "\n")
