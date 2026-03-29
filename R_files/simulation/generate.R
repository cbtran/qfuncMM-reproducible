# Generate data for three-region regions using real voxel coordinates.
# Run this script in the terminal as
# >Rscript R_files/generate.R <delta> <psi> <spec> <noise_level> <nsim> <seed> <outdir>
# where <delta> and <psi> is one of "high", "mid", "low",
# <spec> is one of "std", "fgn", "ar2", "anisotropic", "diag_time"
# <noise_level> specifies the overall noise variance
# <nsim> is the number of simulations to generate,
# <seed> is an integer seed for random number generation.

source("R_files/simulation/generate_3_region.R")
source("R_files/simulation/spatial-anisotropic.R")

# Expect argument such as "mid mid std"
args <- commandArgs(trailingOnly = TRUE)
delta_name <- args[1]
psi_name <- args[2]
covar_setting <- args[3]
message(sprintf("Generating data with %s-%s %s setting...\n", args[1], args[2], covar_setting))
stopifnot(covar_setting %in% c(
  "std", "ar2", "fgn", "anisotropic", "nonsep_matern"
))
noise_level_str <- args[4]
noise_level <- as.numeric(noise_level_str)
n_sim <- as.numeric(args[5])
seed <- as.numeric(args[6])
outdir <- args[7]

# Verify that outdir is a valid writeable directory
if (!dir.exists(outdir)) {
  stop(sprintf("Output directory '%s' does not exist", outdir))
}
if (!file.info(outdir)$isdir) {
  stop(sprintf("'%s' is not a directory", outdir))
}
if (file.access(outdir, 2) != 0) {
  stop(sprintf("Output directory '%s' is not writeable", outdir))
}

nugget_gamma <- 0.1
nugget_eta <- 0.1
if (covar_setting == "ar2") {
  nugget_gamma <- 0
  nugget_eta <- 0
}
k_gamma <- 2
if (covar_setting == "diag_time") {
  k_gamma <- 0
  nugget_eta <- 1
}
delta_fn <- function(x) {
  x * (1 + nugget_eta) / (x * (1 + nugget_eta) + k_gamma + nugget_gamma)
}
delta_seq <- c(0.1, 0.5, 0.7)
names(delta_seq) <- c("low", "mid", "high")

kEta_seq <- sapply(delta_seq, function(y) {
  uniroot(function(x) delta_fn(x) - y, interval = c(0, 20))$root
})

voxel_coords <- readRDS(file.path("R_files", "simulation", "rat_coords.rds"))
sqrd_dist <- lapply(voxel_coords, \(coords) as.matrix(dist(coords))^2)

# For "nonsep", psi is the average spatial correlation at zero temporal lag:
# C(h; 0 | a, b, 1) = exp(-b^2 ||h||^2), which is RBF with rate b.
# For all other settings, use the Matern-5/2 spatial kernel.
if (covar_setting == "nonsep") {
  # C(h; 0 | a, b, 1) = exp(-b^2 * ||h||^2): RBF with rate b
  psi_fn <- function(phi, dist_sqrd) {
    mean(get_cor_mat("rbf", dist_sqrd, phi))
  }
} else if (covar_setting == "nonsep_matern") {
  # Spatial marginal at u=0: Matern-5/2 with rate phi/sqrt(5)
  # (since C(h;0) = (1+phi*h+(phi*h)^2/3)*exp(-phi*h) = matern_5_2(h|phi/sqrt(5)))
  psi_fn <- function(phi, dist_sqrd) {
    mean(get_cor_mat("matern_5_2", dist_sqrd, phi / sqrt(5)))
  }
} else {
  psi_fn <- function(phi, dist_sqrd) {
    mean(get_cor_mat("matern_5_2", dist_sqrd, phi))
  }
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
delta <- delta_seq[delta_name]
psi <- psi_seq[psi_name]

kEta <- kEta_seq[which(delta_seq == delta)]
phi <- phi_seq[, which(psi_seq == psi)] # mid

tau_gamma <- 0.5
tau_eta <- 0.25
if (covar_setting == "diag_time") {
  tau_gamma <- 1
  tau_eta <- NA
}

region_parameters <- data.frame(
  k_eta = rep(kEta, 3),
  phi_gamma = phi,
  tau_gamma = rep(tau_gamma, 3),
  k_gamma = rep(k_gamma, 3),
  nugget_gamma = rep(nugget_gamma, 3),
  mean = c(1, 10, 20),
  var_noise = rep(noise_level, 3),
  row.names = c("region1", "region2", "region3")
)
shared_parameters <- c(tau_eta = tau_eta, nugget = nugget_eta)
corr_true <- c(rho12 = 0, rho13 = 0.35, rho23 = 0.6)

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
    generate_3_region(
      n_sim, voxel_coords, n_timept, corr_true,
      region_parameters, shared_parameters, spatial,
      seed = seed
    )
  },
  "std_diag_time" = {
    spatial <- function(coords, phi_gamma) {
      dist_sqrd_mat <- as.matrix(dist(coords))^2
      get_cor_mat("matern_5_2", dist_sqrd_mat, phi_gamma)
    }
    generate_3_region(
      n_sim, voxel_coords, n_timept, corr_true,
      region_parameters, shared_parameters, spatial,
      seed = seed, covar_spec = "diag_time"
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
    generate_3_region(
      n_sim, voxel_coords, n_timept, corr_true,
      region_parameters, shared_parameters, spatial,
      seed = seed
    )
  },
  "nonsep" = {
    # Non-separability parameter a = 0.5; b = phi_gamma (tuned via psi).
    # gamma covariance: k_gamma * C(h; u | a=0.5, b) + nugget_gamma * I_{LM}
    nonsep_a <- 0.5
    gamma_covar_fn <- function(coords, n_timept, params_row, timesqrd_mat) {
      distsqrd <- as.matrix(dist(coords))^2
      LM <- nrow(coords) * n_timept
      params_row$k_gamma *
        nonsep_cov(distsqrd, timesqrd_mat, nonsep_a, params_row$phi_gamma) +
        params_row$nugget_gamma * diag(LM)
    }
    generate_3_region(
      n_sim, voxel_coords, n_timept, corr_true,
      region_parameters, shared_parameters,
      gamma_covar_fn = gamma_covar_fn,
      seed = seed
    )
  },
  "nonsep_matern" = {
    # Ip & Li (2017) eq. 1.5, nu=5/2, d=3.
    # Closed form: exp(-sqrt(alpha^2 * ||h||^2 + beta^2 * u^2))
    # alpha = phi_gamma (tuned via psi using Matern-1/2 marginal at u=0).
    # beta = 0.5 (fixed temporal rate, analogous to tau_gamma in std).
    nonsep_matern_beta <- 0.5
    gamma_covar_fn <- function(coords, n_timept, params_row, timesqrd_mat) {
      distsqrd <- as.matrix(dist(coords))^2
      LM <- nrow(coords) * n_timept
      params_row$k_gamma *
        nonsep_matern_cov(
          distsqrd, timesqrd_mat,
          params_row$phi_gamma, nonsep_matern_beta
        ) +
        params_row$nugget_gamma * diag(LM)
    }
    generate_3_region(
      n_sim, voxel_coords, n_timept, corr_true,
      region_parameters, shared_parameters,
      gamma_covar_fn = gamma_covar_fn,
      seed = seed
    )
  }
)

out <- list(
  data = three_region,
  setting = list(
    region_parameters = region_parameters,
    shared_parameters = shared_parameters,
    corr_true = corr_true,
    delta = delta,
    psi = psi,
    noise_level = noise_level
  ),
  seed = seed
)

outname <- sprintf("%s-%s.rds", delta_name, psi_name)
outpath <- file.path(outdir, outname)
saveRDS(out, outpath)
message("Saved to ", normalizePath(outpath))
