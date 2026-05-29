

setwd("~/School/2025 - 2026/Thesis Statistics/Analyses/FINAL SCRIPTS/Simulation_260525")

# ==============================================================================
# FINAL GENERATIVE VALIDATION SCRIPT (STAN-CONSISTENT)
# ==============================================================================

rm(list = ls())
gc()

# Set your working directory
# setwd("...") 

library(sf)
library(Matrix)
library(data.table)
library(ggplot2)

set.seed(42)

# ==============================================================================
# 1. GRID & SPATIAL BASIS (PRE-FLIGHT)
# ==============================================================================
# (Replace grid_path with your actual file path)
grid_path <- "C:/Users/lieve/Documents/School/2025 - 2026/Thesis Statistics/Data/New/grid.gpkg"
my_grid <- st_read(grid_path, quiet = TRUE)
coords <- st_coordinates(st_centroid(my_grid))
S <- nrow(coords)

T_total <- 30
T_obs <- 25
N_spatial_bf <- 20

# Create knots
buffer_pct <- 0.10
xlim <- range(coords[,1]) + c(-1, 1) * diff(range(coords[,1])) * buffer_pct
ylim <- range(coords[,2]) + c(-1, 1) * diff(range(coords[,2])) * buffer_pct
knot_grid <- expand.grid(X = seq(xlim[1], xlim[2], length.out = 10),
                         Y = seq(ylim[1], ylim[2], length.out = 10))
km <- kmeans(knot_grid, centers = N_spatial_bf)
knots <- km$centers

# Basis functions
dist_full <- as.matrix(dist(rbind(coords, knots)))
dist_mat_raw <- dist_full[1:S, (S+1):(S+N_spatial_bf)]
dist_mat_anchors_sq <- as.matrix(dist(knots))^2
bandwidth <- quantile(as.matrix(dist(knots))[as.matrix(dist(knots)) > 0], 0.25)
spatial_bf <- exp(-(dist_mat_raw^2) / (2 * bandwidth^2))

# ==============================================================================
# 2. DISPERSAL (PRECISE INDEX MAPPING)
# ==============================================================================
dist_mat_disp <- as.matrix(dist(coords))
dist_mat_disp[dist_mat_disp > 200000] <- 0
diag(dist_mat_disp) <- 0

sp_mat <- drop0(Matrix(dist_mat_disp, sparse = TRUE))
triplet <- summary(sp_mat)

# triplet$i = row (destination), triplet$j = col (source)
to_idx <- triplet$i
from_idx <- triplet$j
dists_vals <- triplet$x / 1000

alpha <- exp(1)
w <- exp(-dists_vals / alpha)

# Row-stochastic normalization
out_sum <- numeric(S)
for (e in seq_along(to_idx)) {
  out_sum[from_idx[e]] <- out_sum[from_idx[e]] + w[e]
}
w_norm <- w / (out_sum[from_idx] + 1e-12)

# ==============================================================================
# 3. LATENT FIELDS & DYNAMICS (TRUE RECOVERY TEST)
# ==============================================================================
# "Logical" values (not the prior modes)
true_thermal_max <- 3.0    # Prior mean is 0
true_T_nodes <- c(2, 11)   # Keep these
true_slope_L <- 1.2
true_slope_R <- 1.3
true_r <- 0.1
true_logit_p <- -1.5

# TRUE SIGMAS (These must be recovered)
true_sigma_phi <- 0.12
true_sigma_gamma <- 0.08
# TRUE RHOs
true_rho_phi <- 80
true_rho_gamma <- 180

# 1. Generate Thermal Niche values
thermal_gradient <- seq(-5, 20, length.out = 400)
static_niche <- true_thermal_max - log1p(exp(-true_slope_L * (thermal_gradient - true_T_nodes[1]))) - log1p(exp(-true_slope_R * (true_T_nodes[2] - thermal_gradient)))

# 2. Covariance Matrices for Independent Lengthscales
Sigma_phi <- exp(-dist_mat_anchors_sq / (2 * true_rho_phi^2))
Sigma_gamma <- exp(-dist_mat_anchors_sq / (2 * true_rho_gamma^2))

L_phi <- t(chol(Sigma_phi + diag(1e-8, N_spatial_bf)))
L_gamma <- t(chol(Sigma_gamma + diag(1e-8, N_spatial_bf)))

# 3. Generate Fields (EXPLICITLY SCALED)
phi_eta <- rnorm(N_spatial_bf)
phi_true <- true_sigma_phi * as.vector(spatial_bf %*% (L_phi %*% phi_eta))
phi_true <- phi_true - mean(phi_true) # Gauge constraint

gamma_eta <- rnorm(N_spatial_bf)
gamma_space <- true_sigma_gamma * as.vector(spatial_bf %*% (L_gamma %*% gamma_eta))
gamma_space <- gamma_space - mean(gamma_space) # Gauge constraint

# 4. Dynamics (same as before)
eps_year <- rnorm(T_total, 0, 0.03)
t_scaled <- seq(-0.5, 0.5, length.out = T_total) 
temp_idx_mat <- matrix(sample(1:400, T_total * S, replace = TRUE), T_total, S)

N_true <- matrix(0, S, T_total)
                 
true_N0_prop <- 0.3 

# Generate N_curr at t=1 using the same thermal + spatial logic as Stan
logN_init <- static_niche[temp_idx_mat[1,]] + phi_true + (gamma_space * t_scaled[1])
N_curr <- exp(logN_init) * true_N0_prop

R0 <- exp(true_r)

for (t in 1:T_total) {
  logK <- static_niche[temp_idx_mat[t,]] + phi_true + (gamma_space * t_scaled[t])
  K <- pmin(pmax(exp(logK), 1e-12), 1e6)
  
  if (t > 1) {
    growth <- (R0 * N_curr) / (1 + ((R0 - 1) / (K + 1e-9)) * N_curr)
    N_new <- numeric(S)
    for (e in seq_along(to_idx)) {
      N_new[to_idx[e]] <- N_new[to_idx[e]] + w_norm[e] * growth[from_idx[e]]
    }
    N_curr <- N_new * exp(eps_year[t])
  }
  N_curr <- pmin(pmax(N_curr, 1e-12), 1e4)
  N_true[, t] <- N_curr
}

# ==============================================================================
# 4. OBSERVATION (MIXTURE MODEL GENERATIVE)
# ==============================================================================
K_visits <- 10
Y_obs <- matrix(0, S, T_total)
p_det <- plogis(true_logit_p)

for (s in 1:S) {
  for (t in (T_total - T_obs + 1):T_total) {
    psi <- 1 - exp(-N_true[s, t])
    # Occupancy mixture
    if (runif(1) < psi) {
      Y_obs[s, t] <- rbinom(1, K_visits, p_det)
    } else {
      Y_obs[s, t] <- 0
    }
  }
}

Y_vec <- as.vector(Y_obs[, (T_total - T_obs + 1):T_total])
map_state_idx <- 1:length(Y_vec)

# ==============================================================================
# 5. DATA EXPORT
# ==============================================================================
t_scaled <- seq(-0.5, 0.5, length.out = T_total)

stan_data <- list(
  S = S, T_total = T_total, T_obs = T_obs,
  N_obs_total = length(Y_vec), N_state_total = S * T_obs,
  Y_vec = as.integer(Y_vec), K_visits_vec = rep(K_visits, length(Y_vec)),
  map_state_idx = map_state_idx,
  E_disp = length(to_idx), dists = dists_vals,
  from_idx = from_idx, to_idx = to_idx,
  N_spatial_bf = N_spatial_bf, spatial_bf = spatial_bf,
  dist_mat_anchors_sq = dist_mat_anchors_sq,
  N_thermal_gradient = length(thermal_gradient),
  thermal_gradient = thermal_gradient,
  temp_idx = as.vector(t(temp_idx_mat)),
  t_scaled = t_scaled
)

saveRDS(stan_data, "stan_data_production.rds")
print("Production data saved.")

# ==============================================================================
# 10. FINAL CHECKS
# ==============================================================================
stopifnot(!anyNA(N_true))
stopifnot(all(is.finite(N_true)))
stopifnot(!anyNA(Y_vec))
stopifnot(all(Y_vec >= 0))

cat("\nSimulation complete: numerically stable + Stan-consistent\n")

my_grid$phi <- phi_true
ggplot(my_grid) +
  geom_sf(aes(fill = phi)) +
  scale_fill_viridis_c() +
  ggtitle("Simulated Spatial GP")