# ============================================================================== #
# SOFT-SCAM: HYBRID IDE MODEL (STATIC NICHE) - CSR SIMULATION ####
# ============================================================================== #
# NOTE: this script is designed to run on a virtual computer

library(cmdstanr)
library(posterior)
library(loo)
library(bayesplot)
library(ggplot2)
library(splines)
library(Matrix)
library(data.table)

# --- 0. DOWNLOAD FILES (OneDrive) ---
library(Microsoft365R)
od <- get_personal_onedrive(auth_type = "device_code")
# Make sure you upload the patched Stan script we just discussed as "Param_CSR_GP_ST_optimized.stan"
od$download_file("Param_CSR_GP_ST_optimized.stan", dest="Param_CSR_GP_ST_optimized.stan", overwrite=TRUE)


# ============================================================================== #
# COMPILATION WITH ABSOLUTE CORE ENFORCEMENT ####
# ============================================================================== #
cat("\n--- 3. Compiling Optimized Model Architecture ---\n")

# Turning stan_threads to FALSE guarantees Stan will 
# never spawn background threads. Each chain is strictly locked to 1 CPU core.
cpp_options <- list(
  stan_threads = FALSE, 
  CXXFLAGS = "-O3 -march=native -mtune=native"
)
cmdstan_make_local(cpp_options = cpp_options)

mod <- cmdstan_model("Param_CSR_GP_ST_optimized.stan", force_recompile = TRUE)


set.seed(42)

# --- 0. DIMENSIONS ---
S <- 1000   
T_total <- 30
T_obs <- 25
N_spatial_bf <- 20 

# ============================================================================== #
# 1. GEOGRAPHY & GP BASIS (20 Knots) ####
# ============================================================================== #
coords <- cbind(runif(S, 0, 10), runif(S, 0, 10))

# Generate knots with a 10% buffer
buffer_pct <- 0.10
xlim <- range(coords[,1]) + c(-1, 1) * diff(range(coords[,1])) * buffer_pct
ylim <- range(coords[,2]) + c(-1, 1) * diff(range(coords[,2])) * buffer_pct

knot_grid <- expand.grid(X = seq(xlim[1], xlim[2], length.out = 10),
                         Y = seq(ylim[1], ylim[2], length.out = 10))

km <- kmeans(knot_grid, centers = N_spatial_bf)
knots <- km$centers

# Build distances
full_dist <- as.matrix(dist(rbind(coords, knots)))
dist_mat_raw <- full_dist[1:S, (S+1):(S+N_spatial_bf)]
dist_mat_anchors <- as.matrix(dist(knots))

# Pre-calculate squared distances for the optimized Stan model
dist_mat_anchors_sq <- dist_mat_anchors^2

# Convert Linear Distances to Smooth Radial Basis Functions (RBF)
bandwidth <- median(dist_mat_anchors[dist_mat_anchors > 0])
spatial_bf <- exp(-(dist_mat_raw^2) / (2 * bandwidth^2))


# ============================================================================== #
# 2. CSR DISPERSAL BASIS (Sparse Network) ####
# ============================================================================== <"
# Calculate distances
dist_mat <- as.matrix(dist(coords)) 

# Create a cutoff threshold to make it "sparse" (4-cell radius)
threshold_distance <- 1.33 
dist_mat[dist_mat > threshold_distance] <- 0 
diag(dist_mat) <- 0 

sp_mat_gen <- as(drop0(dist_mat), "generalMatrix")
sp_mat_R <- as(sp_mat_gen, "RsparseMatrix")

# Extract the vectors for Row-Stochastic Normalization
dists <- sp_mat_R@x
to_idx <- sp_mat_R@j + 1                   # Column indices
row_ptr <- sp_mat_R@p + 1                  # Row pointers
row_ids <- rep(1:S, diff(sp_mat_R@p))      # Explicit Row IDs
E_disp <- length(dists)

cat("Generated CSR Network with", E_disp, "edges.\n")


# ============================================================================== #
# 3. TRUE PARAMETERS (Hybrid v2.1) ####
# ============================================================================== #
true_thermal_max <- 4.5
true_T_nodes <- c(10.0, 25.0)  
true_slope_L <- 1.2
true_slope_R <- 1.2
true_r <- 0.1

true_log_alpha <- 2.7  
true_alpha <- exp(true_log_alpha)

true_rho_phi <- 2.0
true_sigma_phi <- 0.15    
true_sigma_gamma <- 0.05  

# Generate True Habitat Field (phi) and Mean-Center it
L_phi <- exp(-(dist_mat_anchors^2) / (2 * true_rho_phi^2))
phi_true_raw <- spatial_bf %*% (true_sigma_phi * t(chol(L_phi + diag(1e-9, 20))) %*% rnorm(20))
phi_true <- as.vector(phi_true_raw - mean(phi_true_raw))


# ============================================================================== #
# 4. POPULATION DYNAMICS LOOP (CSR VERSION) ####
# ============================================================================== #
# A. Generate the actual temperatures for each site-year
temp_mat <- matrix(NA, T_total, S)

# Shift the temperature gradient to span the entire niche
# Y=0 (South) will be ~28°C. Y=10 (North) will be ~ -2°C.
raw_temp_base <- 28 - (coords[,2] * 3) 

# THIS LOOP IS CRITICAL: It actually fills the matrix with numbers!
for(t in 1:T_total) {
  temp_mat[t,] <- raw_temp_base + (t * 0.04) + rnorm(S, 0, 0.5)
}

# B. Create the Static Thermal Gradient & Map Indices
T_min <- min(temp_mat, na.rm = TRUE) - 2
T_max <- max(temp_mat, na.rm = TRUE) + 2
thermal_gradient <- seq(T_min, T_max, length.out = 400)

temp_idx_mat <- matrix(NA, T_total, S)
for(t in 1:T_total) {
  for(s in 1:S) {
    temp_idx_mat[t,s] <- which.min(abs(thermal_gradient - temp_mat[t,s]))
  }
}

# C. Calculate Sparse Dispersal Weights
W_sparse <- exp(-dist_mat / true_alpha)
W_sparse[dist_mat == 0] <- 0 
diag(W_sparse) <- 0

# Row-normalize so dragons don't spontaneously multiply
W_norm <- W_sparse / rowSums(W_sparse)
W_norm[is.nan(W_norm)] <- 0 

# D. Pre-Calculate the UNCENTERED True Thermal Niche
raw_niche_grad <- plogis(true_slope_L * (thermal_gradient - true_T_nodes[1]), log = TRUE) + 
  plogis(true_slope_R * (true_T_nodes[2] - thermal_gradient), log = TRUE)
static_niche_grad <- true_thermal_max + raw_niche_grad 


# hard cap to match the Stan model update 
N_true <- matrix(0, S, T_total)

# Calculate initial K using the hard limit
log_init_raw <- phi_true + static_niche_grad[temp_idx_mat[1, ]]
K_init <- pmin(pmax(exp(log_init_raw), 1e-6), 1e4)
N_curr <- K_init 

for (t in 1:T_total) {
  log_K_raw <- phi_true + static_niche_grad[temp_idx_mat[t, ]]
  K_cap <- pmin(pmax(exp(log_K_raw), 1e-6), 1e4)
  
  growth <- (exp(true_r) * N_curr) / (1 + ((exp(true_r) - 1) / K_cap) * N_curr)
  
  # Sparse Dispersal Step
  dispersed_raw <- as.numeric(W_norm %*% growth)
  
  # Positivity smoothing and hard cap
  dispersed_pos <- log1p(exp(5.0 * dispersed_raw)) / 5.0
  N_curr <- pmin(pmax(dispersed_pos, 1e-6), 1e4)
  
  N_true[, t] <- N_curr
}


# ============================================================================== #
# ASSEMBLE STAN DATA ####
# ============================================================================== #
Y_obs <- matrix(rbinom(S * T_total, 10, (1 - exp(-N_true)) * plogis(-1.5)), S, T_total)
Y_vec <- as.vector(Y_obs[, (T_total - T_obs + 1):T_total])

idx_seen <- which(Y_vec > 0)
idx_zero <- which(Y_vec == 0)

temp_idx_vec <- as.vector(t(temp_idx_mat))

stan_data_CSR <- list(
  S = S, 
  T_total = T_total, 
  T_obs = T_obs, 
  N_obs_total = S * T_obs,
  N_state_total = S * T_obs,   
  
  Y_vec = Y_vec, 
  K_visits_vec = rep(10, S * T_obs),
  map_state_idx = 1:(S * T_obs),
  
  N_seen = length(idx_seen), 
  N_zero = length(idx_zero),
  idx_seen = idx_seen, 
  idx_zero = idx_zero,
  
  N_spatial_bf = N_spatial_bf, 
  spatial_bf = spatial_bf, 
  dist_mat_anchors = dist_mat_anchors,
  
  E_disp = E_disp,
  dists = dists,
  to_idx = to_idx,
  row_ptr = row_ptr,
  row_ids = row_ids,
  
  dist_mat_anchors_sq = dist_mat_anchors_sq,
  
  N_thermal_gradient = length(thermal_gradient), 
  thermal_gradient = thermal_gradient,
  temp_idx = temp_idx_vec
)

make_init_hybrid <- function() {
  list(
    thermal_max = rnorm(1, -1.0, 0.2),  
    T_nodes = c(runif(1, 4, 6), runif(1, 14, 16)), 
    slope_L = runif(1, 0.8, 1.2), 
    slope_R = runif(1, 0.8, 1.2),
    r = runif(1, 0.05, 0.2),
    logit_p = rnorm(1, -2.0, 0.2),
    log_alpha = rnorm(1, mean = 2.7, sd = 0.5),
    N0_proportion = runif(1, 0.4, 0.6),
    
    rho_phi = runif(1, 1.5, 2.5),        
    rho_gamma = runif(1, 4.5, 5.5),      
    sigma_phi = runif(1, 0.01, 0.08),
    phi_eta = rnorm(20, 0, 0.01),
    
    sigma_gamma = runif(1, 0.001, 0.02),
    gamma_eta = rnorm(20, 0, 0.01),
    sigma_year = runif(1, 0.005, 0.03),
    eps_year_raw = rnorm(T_total, 0, 0.01) 
  )
}

# ============================================================================== #
# FIT THE MODEL ####
# ============================================================================== #
fit_param_ST <- mod$sample(
  data = stan_data_CSR, 
  chains = 2, parallel_chains = 2,
  iter_warmup = 200, iter_sampling = 200,
  adapt_delta = 0.9, max_treedepth = 12,
  init = make_init_hybrid, refresh = 50
)

print(fit_param_ST$summary(c("log_alpha", "thermal_max", "T_nodes", "slope_L", "slope_R")))
fit_param_ST$save_object("fit_parametric_sim.rds")
cat("Model successfully saved to disk!\n")

# ============================================================================== #
# PLOT THE STATIC POSTERIOR THERMAL NICHE (WITH TRUE CURVE) ####
# ============================================================================== #
cat("\n--- Generating Posterior Thermal Niche Plot ---\n")

draws_df <- as_draws_df(fit_param_ST$draws(variables = c(
  "thermal_max", "T_nodes[1]", "T_nodes[2]", "slope_L", "slope_R"
)))

plot_temps <- seq(-5, 30, length.out = 100)

# Match the exact hard-cap math fom Stan
calc_posterior_K <- function(temp) {
  log_K_raw <- draws_df$thermal_max + 
    plogis(draws_df$slope_L * (temp - draws_df$`T_nodes[1]`), log = TRUE) + 
    plogis(draws_df$slope_R * (draws_df$`T_nodes[2]` - temp), log = TRUE)
  
  # The new hard cap geometry (1e4)
  return(pmin(pmax(exp(log_K_raw), 1e-6), 1e4)) 
}

plot_data <- data.frame()

for (temp in plot_temps) {
  K_vals <- calc_posterior_K(temp)
  plot_data <- rbind(plot_data, data.frame(
    Temperature = temp,
    K_median = median(K_vals),
    K_lower = quantile(K_vals, 0.025),
    K_upper = quantile(K_vals, 0.975)
  ))
}

# Calculate the True Niche using the hard cap
calc_true_niche <- function(temp) {
  raw_temp <- plogis(true_slope_L * (temp - true_T_nodes[1]), log = TRUE) + 
    plogis(true_slope_R * (true_T_nodes[2] - temp), log = TRUE)
  
  log_K_raw <- true_thermal_max + raw_temp 
  return(pmin(pmax(exp(log_K_raw), 1e-6), 1e4))
}

plot_data$K_true <- calc_true_niche(plot_data$Temperature)

niche_plot <- ggplot(plot_data, aes(x = Temperature)) +
  geom_ribbon(aes(ymin = K_lower, ymax = K_upper, fill = "Estimated 95% CI"), alpha = 0.3) +
  geom_line(aes(y = K_median, color = "Estimated Median"), linewidth = 1.2) +
  geom_line(aes(y = K_true, color = "True Simulated Niche"), linewidth = 1.2, linetype = "dashed") +
  scale_color_manual(
    name = "Curve",
    values = c("Estimated Median" = "#440154FF", "True Simulated Niche" = "#D55E00")
  ) +
  scale_fill_manual(
    name = "Uncertainty",
    values = c("Estimated 95% CI" = "#440154FF")
  ) +
  labs(
    title = "Posterior Thermal Niche (CSR Simulation)",
    subtitle = "Comparing recovered boundaries to the true simulation (Hard Cap Geometry)",
    x = "Temperature (°C)",
    y = "Carrying Capacity (K)"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

print(niche_plot)


# ============================================================================== #
# TRACEPLOTS & DIAGNOSTICS ####
# ============================================================================== #
cat("\n--- Generating Traceplots ---\n")

# Use Bayesplot to check the health of the HMC chains
color_scheme_set("viridis")

core_params <- c(
  "r", "log_alpha", "logit_p", "N0_proportion",
  "thermal_max", "T_nodes[1]", "T_nodes[2]", "slope_L", "slope_R",
  "sigma_phi", "sigma_gamma", "sigma_year", "lp__"
)

p_trace <- mcmc_trace(fit_param_ST$draws(variables = core_params)) +
  labs(title = "Traceplots: Core Ecological Parameters (CSR Model)",
       subtitle = "Chains should be well-mixed, stationary, and rapidly overlapping") +
  theme_minimal()

print(p_trace)