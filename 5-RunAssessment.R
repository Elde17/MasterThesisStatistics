setwd("~/School/2025 - 2026/Thesis Statistics/Analyses/FINAL SCRIPTS/Real_Data_run")

## Read in data (from supercomputer)
library(cmdstanr)
library(posterior)

# 1. Look inside the folder where you manually extracted the files
# (Adjust "cmdstan_outputs" if you named the folder something else!)
csv_files <- list.files("Thesis_Data", full.names = TRUE, pattern = "ZwaHei_CSR-.*\\.csv$")

# 2. Sanity check: Make sure R actually sees all 6 chains!
print(csv_files)

# 3. Rebuild the model object
fit_final <- as_cmdstan_fit(csv_files)

# 4. Read your draws (ensure the path matches where you saved it)
# combined_draws <- readRDS("Thesis_Data/fit_final_ZwaHei_Parametric_ST.rds")
combined_draws <- fit_final$draws()

# ============================================================================== #
# QUICK CONVERGENCE CHECK SCRIPT ####
# ============================================================================== #
library(bayesplot)
library(ggplot2)
library(cmdstanr)

core_params <- c(
  # --- 1. Global Ecological Rates ---
  "r",                # Intrinsic growth rate
  "log_alpha",          # Dispersal (Spectral Diffusion coefficient)
  "logit_p",          # Detection probability (The "Anchor" for abundance)
  "N0_proportion",    # Initial state (Helps check if burn-in was successful)
  
  # --- 2. Thermal Niche (Double Logistic / Double Sigmoid) ---
  "thermal_max",      # Peak carrying capacity (log-scale)
  "T_nodes[1]",       # T_start: Cold-limit threshold
  "T_nodes[2]",       # T_end: Heat-limit threshold
  "slope_L",          # Sharpness of the cold-side boundary
  "slope_R",          # Sharpness of the heat-side boundary
  
  # --- 3. Hierarchical Variances (The "Noise" structure) ---
  "sigma_phi",        # Magnitude of the Spatial Baseline
  "sigma_gamma",      # Magnitude of the Space-Time Trends
  "sigma_year",       # Global year-to-year stochasticity
  
  "rho_phi",        # Magnitude of the Spatial Baseline
  "rho_gamma",      # Magnitude of the Space-Time Trends
  
  # --- 4. The "Health" Metric ---
  "lp__"              # Log-posterior (Crucial for detecting geometry issues)
)

# ============================================================================== #
# DYNAMICALLY SAMPLE LATENT PARAMETERS
# ============================================================================== #
library(posterior) # Ensure this is loaded

# 1. Extract the draws object (this is what posterior functions expect)
draws_obj <- fit_final$draws()

# 2. Now you can get the variable names successfully
all_vars <- variables(draws_obj)

# 3. Use grep to find the latent parameters
phi_vars   <- grep("^phi_eta\\[", all_vars, value = TRUE)
gamma_vars <- grep("^gamma_eta\\[", all_vars, value = TRUE)
eps_vars   <- grep("^eps_year_raw\\[", all_vars, value = TRUE)

# 4. Randomly sample 4 from each
set.seed(42) 
sampled_latent <- c(
  sample(phi_vars, min(length(phi_vars), 4)),
  sample(gamma_vars, min(length(gamma_vars), 4)),
  sample(eps_vars, min(length(eps_vars), 4))
)

# 5. Add them to your core parameters list
core_params <- c("r", "log_alpha", "logit_p", "N0_proportion", "thermal_max", 
                 "T_nodes[1]", "T_nodes[2]", "slope_L", "slope_R", 
                 "sigma_phi", "sigma_gamma", "sigma_year", "rho_phi", "rho_gamma")

extended_params <- c(core_params, sampled_latent)

# Verify
print(sampled_latent)

# Combine with your existing core_params
extended_params <- c(core_params, sampled_latent)

# Now use 'extended_params' for your summaries and plots
cat("\n--- Updated Summary with Latent Parameters ---\n")
print(fit_final$summary(extended_params),n=50)

cat("\n--- 1. Statistical Diagnostic Summary ---\n")

# 1. Print core ecological parameters
print(fit_final$summary(core_params))

# 2. Check for HMC-specific issues (Divergences)
# High divergences suggest the model geometry is too complex for the current step size
fit_final$diagnostic_summary()

cat("\n--- 2. Visual Traceplot Inspection ---\n")

# 3. Plot traceplots
color_scheme_set("viridis")
trace <- mcmc_trace(fit_final$draws(c(extended_params,"lp__"))) +
  theme_minimal()
trace
ggsave("Trace.pdf", plot = trace, 
       width = 10, height = 7, dpi = 300)

cat("\n--- 3. Numerical Convergence Thresholds ---\n")

# 4. Extract summary to check R-hat and ESS specifically
fit_summ <- fit_final$summary(extended_params)

# Check if any R-hat is above 1.05 (Common threshold for non-convergence)
bad_rhat <- fit_summ[fit_summ$rhat > 1.05, ]
if(nrow(bad_rhat) > 0) {
  cat("WARNING: The following parameters failed to converge (R-hat > 1.05):\n")
  print(bad_rhat$variable)
} else {
  cat("SUCCESS: All core R-hat values are < 1.05.\n")
}

# Check if Effective Sample Size is too low (below 400 is a common warning sign)
low_ess <- fit_summ[fit_summ$ess_bulk < 400, ]
if(nrow(low_ess) > 0) {
  cat("NOTE: The following parameters have low Effective Sample Size (ESS < 400):\n")
  print(low_ess$variable)
}


## ISSUE WITH BIMODALITY ##
# NOTE: the current model has unimodality, which is good
# Check log-likelihood
# Compare log-posterior across chains
library(bayesplot)
library(ggplot2)

lp <- fit_final$draws("lp__")
mcmc_trace(lp) + 
  labs(title = "Log-Posterior Traceplot", subtitle = "Higher is better")

library(posterior)

# Extract draws and subset to chains
combined_draws <- fit_final$draws()
final_draws <- subset_draws(combined_draws, chain = c(1,2,3,4,5,6))

# Re-check diagnostics on the subset
# summarise_draws(final_draws, "rhat", "ess_bulk")

# Plot the trace for the key biological parameters and the log-posterior
color_scheme_set("viridis")
arr <- as_draws_array(final_draws)
subset_arr <- arr[, , core_params]

mcmc_trace(subset_arr)

# Correlation between alpha and sigma_phi
corr <- mcmc_pairs(subset_arr)
corr
ggsave("Corr.pdf", plot = corr, 
       width = 15, height = 20)



## Parameter estimates + 95% CI
# Extract and summarize
# We use 0.025 and 0.975 for the 95% CI
summary_table <- summarise_draws(
  subset_draws(final_draws, variable = core_params),
  "mean", "median", "sd", 
  ~quantile(.x, probs = c(0.025, 0.975), na.rm = TRUE),
  "rhat", "ess_bulk"
) %>%
  dplyr::filter(variable %in% core_params)

# Clean up the names for the table
colnames(summary_table)[5:6] <- c("q2.5", "q97.5")

print(summary_table)


library(cmdstanr)
library(posterior)


library(posterior)
library(dplyr)

# 1. Extract all latent parameters
latent_vars <- all_vars[grep("^(phi_eta|gamma_eta|eps_year_raw)\\[", all_vars)]

# 2. Summarize (calculating the 95% CI)
latent_summ <- summarise_draws(
  subset_draws(final_draws, variable = latent_vars),
  "q2.5" = ~quantile(.x, 0.025),
  "q97.5" = ~quantile(.x, 0.975)
)

# 3. Assess significance correctly using backticks
latent_summ <- latent_summ %>%
  mutate(is_significant = `2.5%` > 0 | `97.5%` < 0)

# 4. Count and Report
summary_counts <- latent_summ %>%
  mutate(type = gsub("\\[.*", "", variable)) %>% 
  group_by(type) %>%
  summarise(
    n_total = n(),
    n_significant = sum(is_significant),
    pct_significant = round(n_significant / n_total * 100, 1)
  )

cat("\n--- Corrected Latent Parameter Significance Summary ---\n")
print(summary_counts)

# 5. View the parameters that ARE significant
significant_params <- latent_summ %>% filter(is_significant == TRUE)
print(significant_params)


## PLOT
## Eps
# 1. Prepare data and calculate significance in one pipe
eps_summ <- latent_summ %>%
  filter(grepl("eps_year_raw", variable)) %>%
  mutate(year = as.numeric(str_extract(variable, "\\d+"))) %>%
  # Ensure the significance calculation is done on THIS specific dataframe
  mutate(is_significant = `2.5%` > 0 | `97.5%` < 0)

# 2. Plot
ggplot(eps_summ, aes(x = year, y = median)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_pointrange(aes(ymin = `2.5%`, ymax = `97.5%`, 
                      color = is_significant), size = 0.5) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"))+
  theme_minimal() +
  labs(title = "Annual Stochastic Shocks (Epsilons)",
       x = "Year", y = "Deviation from mean", color = "Significant?")


## Phi & gamma
library(ggplot2)
library(patchwork) # For side-by-side plotting
library(dplyr)

# 1. Prepare data for Phi and Gamma
# We create 'sig_type' to hold: Positive, Negative, or Not Significant
prepare_spatial_data <- function(var_pattern) {
  latent_summ %>%
    filter(grepl(var_pattern, variable)) %>%
    mutate(id = as.numeric(str_extract(variable, "\\d+"))) %>%
    mutate(sig_type = case_when(
      `2.5%` > 0 ~ "Positive",
      `97.5%` < 0 ~ "Negative",
      TRUE ~ "Not Significant"
    )) %>%
    mutate(sig_type = factor(sig_type, levels = c("Positive", "Not Significant", "Negative")))
}

phi_plot_data <- anchors_df %>% mutate(id = 1:n()) %>% left_join(prepare_spatial_data("phi_eta"), by = "id")
gamma_plot_data <- anchors_df %>% mutate(id = 1:n()) %>% left_join(prepare_spatial_data("gamma_eta"), by = "id")

# 1. Combine data into one 'tidy' dataframe
# We add a 'Field' column to distinguish phi from gamma
plot_data_phi <- phi_plot_data %>% mutate(Field = "Static Habitat (phi)")
plot_data_gamma <- gamma_plot_data %>% mutate(Field = "Spatiotemporal Trend (gamma)")

combined_data <- bind_rows(plot_data_phi, plot_data_gamma)

# 2. Define Publication-Quality Colors
# Using a colorblind-friendly palette
pub_colors <- c("Positive" = "Green", "Negative" = "Red", "Not Significant" = "grey85")

# 3. Create the Publication Plot
ggplot(combined_data) +
  # Add the mainland background
  geom_sf(data = my_grid_mainland, color = "grey95", fill = "white") +
  # Plot the knots
  geom_point(aes(x = X, y = Y, color = sig_type, size = sig_type)) +
  # Force consistency
  scale_color_manual(values = pub_colors) +
  scale_size_manual(values = c("Positive" = 3.5, "Negative" = 3.5, "Not Significant" = 1.5)) +
  # Facet for clean side-by-side comparison
  facet_wrap(~Field) +
  # Remove clutter
  theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "bottom") +
  labs(color = "Significance:", size = "Significance:")



# ============================================================================== #
# POST-PROCESSING "BRIDGE" ####
# ============================================================================== #

# 1. Setup Data Environment
# Ensure stan_data_real is available (for the reconstruction branch)
# Thus, first run the StanDataPreparation script
stan_data <- stan_data_real

if (!exists("stan_data")) stan_data <- stan_data_real
list2env(stan_data, envir = .GlobalEnv)

# 2. Extract Variable List from current fit
all_vars <- variables(fit_final$draws())

# 3. Conditional Logic: Direct Extraction vs. Reconstruction
if ("N_state" %in% all_vars) {
  cat("\n--- [Bridge] Detected N_state in Stan output: Extracting directly ---\n")
  
  N_draws <- fit_final$draws("N_state")
  N_mat <- matrix(apply(as_draws_matrix(N_draws), 2, median), 
                  nrow = T_obs, ncol = S, byrow = TRUE)
  
} else {
  cat("\n--- [Bridge] N_state NOT found: Reconstructing in R (Numerical Safeguards active) ---\n")
  
  # --- STEP A: Extract posterior medians from fit_final ---
  phi_med      <- apply(as_draws_matrix(fit_final$draws("phi")), 2, median)
  gamma_med    <- apply(as_draws_matrix(fit_final$draws("gamma_smooth")), 2, median)
  eps_med      <- apply(as_draws_matrix(fit_final$draws("eps_year")), 2, median)
  r_med        <- median(as_draws_matrix(fit_final$draws("r")))
  N0_med       <- median(as_draws_matrix(fit_final$draws("N0_proportion")))
  alpha_med    <- median(as_draws_matrix(fit_final$draws("alpha")))
  
  # Thermal parameters
  thermal_max_med <- median(as_draws_matrix(fit_final$draws("thermal_max")))
  slope_L_med     <- median(as_draws_matrix(fit_final$draws("slope_L")))
  slope_R_med     <- median(as_draws_matrix(fit_final$draws("slope_R")))
  T_nodes_med     <- apply(as_draws_matrix(fit_final$draws("T_nodes")), 2, median)
  
  # Thermal Niche
  thermal_niche <- thermal_max_med - 
    log(1 + exp(-slope_L_med * (thermal_gradient - T_nodes_med[1]))) - 
    log(1 + exp(-slope_R_med * (T_nodes_med[2] - thermal_gradient)))
  
  # Dispersal weights
  w <- exp(-dists / alpha_med)
  out_sum <- aggregate(w, by=list(from_idx), FUN=sum)$x
  w_norm <- w / out_sum[from_idx]
  
  # --- STEP B: Hardened Simulation Loop ---
  N_state_total <- S * T_obs
  N_state_vec <- numeric(N_state_total)
  R0 <- exp(r_med)
  R0m1 <- R0 - 1.0
  
  logN <- thermal_niche[temp_idx[1:S]] + phi_med + (gamma_med * t_scaled[1]) + eps_med[1]
  N_curr <- exp(logN) * N0_med
  
  for (t in 2:T_total) {
    logK <- thermal_niche[temp_idx[((t - 1) * S + 1) : (t * S)]] + phi_med + (gamma_med * t_scaled[t])
    K <- pmax(exp(logK), 1e-6) 
    
    N_curr <- (R0 * N_curr) / (1.0 + (R0m1 / K) * N_curr)
    
    N_new <- numeric(S)
    for (e in 1:E_disp) {
      N_new[to_idx[e]] <- N_new[to_idx[e]] + w_norm[e] * N_curr[from_idx[e]]
    }
    N_curr <- N_new * exp(eps_med[t])
    
    N_curr[is.na(N_curr) | N_curr < 0] <- 0
    
    if (t > T_total - T_obs) {
      obs_t <- t - (T_total - T_obs)
      N_state_vec[((obs_t - 1) * S + 1):(obs_t * S)] <- N_curr
    }
  }
  N_mat <- matrix(N_state_vec, nrow = T_obs, ncol = S, byrow = TRUE)
}

# 4. Final Cleanup
gc()
cat("--- [Bridge] Success! N_mat is available for downstream plotting ---\n")

# 3. Final Verification
cat("NAs remaining:", sum(is.na(N_mat)), "\n")
summary(as.vector(N_mat))
cat("Success! N_mat is ready for plotting.\n")

library(ggplot2)
library(sf)
library(viridis)

# 1. Spatial Map: Mean Abundance
my_grid_final <- my_grid_mainland
my_grid_final$mean_N <- colMeans(N_mat)

p_spatial <- ggplot(my_grid_final) +
  geom_sf(aes(fill = mean_N), color = NA) +
  scale_fill_viridis_c(option = "magma", trans = "log1p", 
                       name = "Mean\nAbundance") +
  theme_minimal() +
  labs(title = "Mean Latent Abundance (2000-2024)",
       subtitle = "Reconstructed from latent parameters")

print(p_spatial)
ggsave("map_mean_abundance.pdf", width = 8, height = 8)

# 2. Temporal Plot: Total Population Trend
plot_df <- data.frame(Year = 2000:2024, Total_N = rowSums(N_mat))

p_temporal <- ggplot(plot_df, aes(x = Year, y = Total_N)) +
  geom_line(color = "#21918cFF", size = 1.2) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "Estimated Continental Population Trend",
       y = "Total Latent Abundance", x = "Year")

print(p_temporal)
ggsave("plot_population_trend.pdf", width = 8, height = 5)


# ============================================================================== #
# 9b. POSTERIOR PREDICTIVE CHECKS (Conditional Bridge: N_state or Reconstruction)
# ============================================================================== #
cat("\n--- 9. Simulating Posterior Predictive Distribution ---\n")

library(posterior)
library(bayesplot)
library(ggplot2)

# 1. Grab observation data
Y_obs <- stan_data_real$Y_vec
K_visits <- stan_data_real$K_visits_vec
map_idx <- stan_data_real$map_state_idx

# 2. Setup draws
set.seed(42)
n_ppc_draws <- 1000
draws_df <- as_draws_df(fit_final$draws())
sample_rows <- sample(1:nrow(draws_df), n_ppc_draws, replace=T)
y_rep_mat <- matrix(NA, nrow = n_ppc_draws, ncol = length(Y_obs))

# --- HELPER: Reconstruction Function --- 
reconstruct_N_for_ppc <- function(draw, S, T_total, T_obs, temp_idx, t_scaled, from_idx, to_idx, E_disp) {
  # Extract params for this draw
  phi <- as.numeric(draw[grep("phi\\[", names(draw))])
  gamma <- as.numeric(draw[grep("gamma_smooth\\[", names(draw))])
  eps_year <- as.numeric(draw[grep("eps_year\\[", names(draw))])
  r <- as.numeric(draw["r"])
  N0 <- as.numeric(draw["N0_proportion"])
  alpha <- as.numeric(draw["alpha"])
  
  # Calculate thermal niche and w_norm for this draw
  T_nodes <- c(as.numeric(draw["T_nodes[1]"]), as.numeric(draw["T_nodes[2]"]))
  thermal_niche <- as.numeric(draw["thermal_max"]) - 
    log(1 + exp(-as.numeric(draw["slope_L"]) * (thermal_gradient - T_nodes[1]))) - 
    log(1 + exp(-as.numeric(draw["slope_R"]) * (T_nodes[2] - thermal_gradient)))
  
  w <- exp(-dists / alpha)
  out_sum <- aggregate(w, by=list(from_idx), FUN=sum)$x
  w_norm <- w / out_sum[from_idx]
  
  # Simulation
  N_state_vec <- numeric(S * T_obs)
  R0 <- exp(r)
  R0m1 <- R0 - 1.0
  N_curr <- exp(thermal_niche[temp_idx[1:S]] + phi + (gamma * t_scaled[1]) + eps_year[1]) * N0
  
  for (t in 2:T_total) {
    logK <- thermal_niche[temp_idx[((t - 1) * S + 1) : (t * S)]] + phi + (gamma * t_scaled[t])
    K <- pmax(exp(logK), 1e-6)
    N_curr <- (R0 * N_curr) / (1.0 + (R0m1 / K) * N_curr)
    
    N_new <- numeric(S)
    for (e in 1:E_disp) N_new[to_idx[e]] <- N_new[to_idx[e]] + w_norm[e] * N_curr[from_idx[e]]
    N_curr <- N_new * exp(eps_year[t])
    N_curr[is.na(N_curr) | N_curr < 0] <- 0
    
    if (t > T_total - T_obs) {
      obs_t <- t - (T_total - T_obs)
      N_state_vec[((obs_t - 1) * S + 1):(obs_t * S)] <- N_curr
    }
  }
  return(N_state_vec)
}

# 3. Execution (Logic Branch)
all_vars <- variables(fit_final$draws())
cat("Simulating PPC and generating y_rep...\n")

for (i in 1:n_ppc_draws) {
  row_idx <- sample_rows[i]
  
  # A. Get N_state for this draw
  if ("N_state[1]" %in% names(draws_df)) {
    # Direct extraction
    N_sim <- as.numeric(draws_df[row_idx, grep("N_state\\[", names(draws_df))])
  } else {
    # Reconstruction
    N_sim <- reconstruct_N_for_ppc(draws_df[row_idx, ], S, T_total, T_obs, temp_idx, t_scaled, from_idx, to_idx, E_disp)
  }
  
  # B. Observation Process (Common to both branches)
  logit_p <- as.numeric(draws_df[row_idx, "logit_p"])
  p_draw <- plogis(logit_p)
  
  # Calculate psi and predict Y_rep
  N_at_obs <- N_sim[map_idx]
  psi <- 1 - exp(-pmax(N_at_obs, 0) - 1e-12)
  y_rep_mat[i, ] <- rbinom(length(Y_obs), size = K_visits, prob = psi * p_draw)
}

cat("y_rep generated successfully! Plotting...\n")

# 4. Plot the PPCs
color_scheme_set("brightblue")

# 1. Density Overlay
# The x-axis represents the counts of the observed variable
p1 <- ppc_dens_overlay(y = Y_obs, yrep = y_rep_mat) +
  coord_cartesian(xlim = c(0, max(Y_obs) + 5)) +
  labs(title = "PPC: Density Overlay",
       x = "Observed Counts", 
       y = "Density") +
  theme_minimal()

# 2. Proportion of Zeros
# The x-axis is the statistic calculated (proportion)
p2 <- ppc_stat(y = Y_obs, yrep = y_rep_mat, stat = function(y) mean(y == 0)) +
  labs(title = "PPC: Proportion of Zeros",
       x = "Proportion of Zero Counts", 
       y = "Frequency") +
  theme_minimal()

# 3. Maximum Count
# The x-axis is the maximum value observed in the dataset
p3 <- ppc_stat(y = Y_obs, yrep = y_rep_mat, stat = "max") +
  labs(title = "PPC: Maximum Count", 
       x = "Maximum Count Value", 
       y = "Frequency") +
  theme_minimal()

print(p1)
print(p2)
print(p3)

ggsave("ppc_density.pdf", p1, width = 8, height = 5)
ggsave("ppc_zeros.pdf", p2, width = 8, height = 5)
ggsave("ppc_max.pdf", p3, width = 8, height = 5)



# --- Calculate Pearson Chi-Square Statistic ---

# 1. Calculate the statistic for the real data (observed)
# Assuming Y_obs is your real data vector
chi2_obs <- sum(((Y_obs - mean(Y_obs))^2) / mean(Y_obs))

# 2. Calculate for every row in y_rep_mat (simulated data)
chi2_rep <- apply(y_rep_mat, 1, function(y_sim) {
  sum(((y_sim - mean(y_sim))^2) / mean(y_sim))
})

# 3. Calculate Bayesian p-value
ppp <- mean(chi2_rep > chi2_obs)

# 4. Plot the distribution
chi2_df <- data.frame(chi2 = chi2_rep)

p4 <- ggplot(chi2_df, aes(x = chi2)) +
  geom_histogram(fill = "skyblue", color = "white", bins = 50) +
  geom_vline(xintercept = chi2_obs, color = "red", linetype = "dashed", size = 1) +
  theme_minimal() +
  labs(title = "PPC: Chi-Squared Discrepancy",
       subtitle = paste("Bayesian p-value =", round(ppp, 3)),
       x = "Chi-squared value",
       y = "Frequency")
p4

library(patchwork)

# Arrange into a 2x2 layout
# Density (Top Left), Zeros (Top Right), Max (Bottom Left), Chi-Sq (Bottom Right)
combined_plot <- (p1 | p2) / (p3 | p4) + 
  plot_annotation(tag_levels = 'A', 
                  title = "Model Diagnostic Suite",
                  subtitle = "Posterior Predictive Checks and Goodness-of-Fit")

# Print to viewer
print(combined_plot)

# Save the final thesis-ready figure
ggsave("PPC_all.pdf", combined_plot, width = 10, height = 7)


# ============================================================================== # #
# 10. POSTERIOR DISTRIBUTIONS ####
# ============================================================================== # #
cat("\n--- 10. Posterior Distributions ---\n")
library(truncnorm)
# Convert subsetted draws to a data frame for ggplot
draws_df <- as_draws_df(subset_draws(final_draws, variable = core_params))
library(tidyr)
library(dplyr)
draws_long <- draws_df |>
  as_tibble() |>
  select(all_of(core_params)) |>
  pivot_longer(everything(), names_to = "Parameter", values_to = "Value")

estimates <- draws_long |>
  group_by(Parameter) |>
  summarize(Median = median(Value))

ggplot(draws_long, aes(x = Value)) +
  geom_density(fill = "skyblue", alpha = 0.6) +
  geom_text(data = estimates,
            aes(x = Inf, y = Inf,
                label = paste0("Median: ", round(Median, 3))),
            hjust = 1.1, vjust = 1.1,
            inherit.aes = FALSE,
            fontface = "bold") +
  facet_wrap(~Parameter, scales = "free") +
  theme_minimal() +
  labs(title = "Posterior Distributions")

# ============================================================================== #
# 10b. POSTERIOR VS. PRIOR DISTRIBUTIONS (TRUE PRIORS)
# ============================================================================== #

library(tidyr)
library(dplyr)
library(ggplot2)
library(posterior)

# 1. Define params
# Ensure these names are the "final" versions you want in your plots
core_params <- c("r", "log_alpha", "logit_p", "N0_proportion", "thermal_max", 
                 "T_nodes[1]", "T_nodes[2]", "slope_L", "slope_R", 
                 "sigma_phi", "sigma_gamma", "sigma_year", "rho_phi", "rho_gamma")

# 2. Extract and Clean Posterior Draws
# Since the colnames already match "T_nodes[1]", no renaming is required.
posterior_long <- as_draws_df(fit_final$draws()) %>%
  as_tibble() %>%
  # Just select the parameters directly. 
  # 'any_of' handles the bracketed names without needing backticks.
  select(any_of(core_params)) %>%
  pivot_longer(everything(), names_to = "Parameter", values_to = "Value") %>%
  mutate(Distribution = "Posterior")

# 3. Create and Clean Prior Data
# Use tibble() instead of data.frame() to prevent R from changing [ ] to .
n_draws <- 1000
T_nodes_prior <- matrix(rnorm(n_draws * 2, mean = 10, sd = 8), ncol = 2)

# Define the Half-t helper function
r_half_t <- function(n, df, mu, sigma) {
  # Generate from a t-distribution, then scale and shift, 
  # then take the absolute value for the 'half' aspect
  # Alternatively, use standard practice of folding a normal/t
  # This matches common usage for half-t priors:
  return(abs(rt(n, df = df) * sigma + mu))
}
prior_df <- tibble(
  r = rtruncnorm(n_draws, 0.3, 0.3),
  log_alpha = rnorm(n_draws, 2.5, 1),
  logit_p = rnorm(n_draws, -2, 1),
  N0_proportion = rbeta(n_draws, 2, 2),
  thermal_max = rnorm(n_draws, 0, 3),
  "T_nodes[1]" = T_nodes_prior[,1],
  "T_nodes[2]" = T_nodes_prior[,2],
  slope_L = rlnorm(n_draws, 0, 0.5),
  slope_R = rlnorm(n_draws, 0, 0.5),
  sigma_phi = abs(rnorm(n_draws, 0, 0.1)),
  sigma_gamma = abs(rnorm(n_draws, 0, 0.1)),
  sigma_year = r_half_t(n_draws, 3, 0, 0.2),
  rho_phi = rlnorm(n_draws, log(100), 0.5),
  rho_gamma = rlnorm(n_draws, log(150), 0.5)
)

# Continue with your existing code
prior_long <- prior_df %>%
  pivot_longer(everything(), names_to = "Parameter", values_to = "Value") %>%
  mutate(Distribution = "Prior")

# 4. Combine and Plot
combined_long <- bind_rows(prior_long, posterior_long)

# Plot
p_priors <- ggplot(combined_long, aes(x = Value, fill = Distribution, color = Distribution)) +
  geom_density(alpha = 0.5, linewidth = 0.6) +
  facet_wrap(~Parameter, scales = "free", ncol = 4) +
  scale_fill_manual(values = c("Posterior" = "skyblue", "Prior" = "grey70")) +
  scale_color_manual(values = c("Posterior" = "dodgerblue4", "Prior" = "grey40")) +
  theme_minimal() +
  labs(title = "Prior vs. Posterior Overlays",
       subtitle = "Assessing Prior Contraction and Model Learning",
       y = "Density", x = "Parameter Value")

print(p_priors)

ggsave("PostPrior.pdf", p_priors, width = 10, height = 7)


# ============================================================================== #
# 11. RECOVERING THE THERMAL NICHE ####
# ============================================================================== #
cat("\n--- 11. Reconstructing Thermal Performance Curve (TPC) ---\n")

library(ggplot2)
library(posterior)

# 1. Extract necessary parameters
# Ensure these exist in the draws_df
draws_mat <- as_draws_matrix(fit_final$draws(c("thermal_max", "slope_L", "slope_R", "T_nodes")))

# 2. Define the TPC reconstruction function (the inverse of your Stan logic)
# Note: Stan's log1p_exp(z) = log(1 + exp(z))
calc_tpc <- function(temp, t_max, s_L, s_R, t_node1, t_node2) {
  log_K <- t_max - log(1 + exp(-s_L * (temp - t_node1))) - log(1 + exp(-s_R * (t_node2 - temp)))
  return(exp(log_K)) # Exponentiate back to natural scale
}

# 3. Reconstruct across all draws (or a subset for speed if needed)
# thermal_gradient must be defined in your workspace (from stan_data)
# Let's take a sample of 500 draws to keep the ribbon calculation snappy
n_plot_draws <- 500
idx <- sample(1:nrow(draws_mat), n_plot_draws)
sub_draws <- draws_mat[idx, ]

# Create a matrix to store the curves (Rows = Draws, Cols = Temp Points)
tpc_matrix <- matrix(NA, nrow = n_plot_draws, ncol = length(thermal_gradient))

for(i in 1:n_plot_draws) {
  tpc_matrix[i, ] <- calc_tpc(thermal_gradient, 
                              sub_draws[i, "thermal_max"], 
                              sub_draws[i, "slope_L"], 
                              sub_draws[i, "slope_R"], 
                              sub_draws[i, "T_nodes[1]"], 
                              sub_draws[i, "T_nodes[2]"])
}

# 4. Summarize into Plotting Dataframe
df_tpc <- data.frame(
  Temp = thermal_gradient,
  K_Median = apply(tpc_matrix, 2, median),
  K_Lower  = apply(tpc_matrix, 2, quantile, 0.025),
  K_Upper  = apply(tpc_matrix, 2, quantile, 0.975)
)

# 5. Plot the Thermal Niche
# Update this section in your ggplot code
p_tpc <- ggplot(df_tpc, aes(x = Temp, y = K_Median)) +
  geom_ribbon(aes(ymin = K_Lower, ymax = K_Upper), fill = "firebrick", alpha = 0.3) +
  geom_line(color = "firebrick", linewidth = 1.5) +
  theme_minimal() +
  labs(title = "Thermal Scaling Effect on Carrying Capacity",
       subtitle = "Relative scaling factor of K as a function of temperature",
       y = "Thermal scaling factor (exp(niche))",
       x = "Mean annual temperature (°C)")

print(p_tpc)
ggsave("Niche.pdf", width = 7, height = 4)

# ============================================================================== #
# 12. RANGE DYNAMICS & VELOCITY ####
# ============================================================================== #
library(dplyr)
library(ggplot2)
library(posterior)
library(sf)

cat("\n--- 12. Range Velocity Analysis (Universal Bridge) ---\n")

# 0. Safety Check
if(!exists("grid_y_coords")) {
  grid_y_coords <- st_coordinates(st_centroid(my_grid_mainland))[,2]
}

# --- MARGIN FUNCTION ---
get_margin_draw <- function(abundance_vec, y_coords, prob, noise_threshold = 0.1) {
  clean_n <- ifelse(abundance_vec < noise_threshold, 0, abundance_vec)
  if(sum(clean_n) == 0) return(NA)
  ord <- order(y_coords)
  y_sorted <- y_coords[ord]
  n_sorted <- clean_n[ord]
  cum_n <- cumsum(n_sorted) / sum(n_sorted)
  return(y_sorted[which(cum_n >= prob)[1]])
}

# --- 1. THE BRIDGE: Extract or Reconstruct N_state ---
all_vars <- variables(fit_final$draws())
quantiles_to_track <- c(0.025, 0.25, 0.50, 0.75, 0.975)
quantile_names <- c("Trailing (2.5%)", "25th Percentile", "Core Median (50%)", "75th Percentile", "Leading (97.5%)")

if ("N_state[1]" %in% all_vars) {
  cat("--- [Bridge] Detected N_state in Stan output: Extracting from fit_final ---\n")
  N_draws_obj <- as_draws_matrix(fit_final$draws("N_state"))
  n_draws <- nrow(N_draws_obj)
  # Reshape: [draws, sites*years] -> [draws, sites, years]
  N_array <- array(as.matrix(N_draws_obj), dim = c(n_draws, S, T_obs))
  
} else {
  cat("--- [Bridge] N_state NOT found: Reconstructing in R (Subsetting to 200 draws for RAM safety) ---\n")
  n_draws <- 200 # RAM safeguard
  draws_df <- as_draws_df(fit_final$draws())
  sample_idx <- sample(1:nrow(draws_df), n_draws)
  
  N_array <- array(NA, dim = c(n_draws, S, T_obs))
  
  for(i in 1:n_draws) {
    N_array[i, , ] <- matrix(reconstruct_N_for_ppc(draws_df[sample_idx[i], ], S, T_total, T_obs, 
                                                   temp_idx, t_scaled, from_idx, to_idx, E_disp), 
                             nrow = S, ncol = T_obs)
  }
}

# --- 2. Calculate Margins ---
margin_list <- list()
for(q_name in quantile_names) margin_list[[q_name]] <- array(NA, dim = c(n_draws, T_obs))

cat("Calculating margins across", n_draws, "draws...\n")
for(d in 1:n_draws) {
  for(t in 1:T_obs) {
    for(i in seq_along(quantiles_to_track)) {
      margin_list[[quantile_names[i]]][d, t] <- get_margin_draw(N_array[d, , t], grid_y_coords, quantiles_to_track[i], 0.1)
    }
  }
}
rm(N_array); gc() # Clear memory

# --- 3. Summarize and Plot ---
summarize_margins <- function(margin_array) {
  data.frame(Year = 2000:2024, # Adjust to your actual years
             Median = apply(margin_array, 2, median, na.rm=TRUE) / 1000,
             Lower = apply(margin_array, 2, quantile, 0.025, na.rm=TRUE) / 1000,
             Upper = apply(margin_array, 2, quantile, 0.975, na.rm=TRUE) / 1000)
}

df_plot <- data.frame()
vel_stats <- list()

for(q_name in quantile_names) {
  df_temp <- summarize_margins(margin_list[[q_name]]) %>% mutate(Edge = q_name)
  df_plot <- bind_rows(df_plot, df_temp)
  
  vels <- apply(margin_list[[q_name]], 1, function(y) {
    if(any(is.na(y))) return(NA)
    coef(lm(y/1000 ~ c(2000:2024)))[2] # Update Year vector here
  })
  vel_stats[[q_name]] <- c(Median = median(vels, na.rm=TRUE), quantile(vels, c(0.025, 0.975), na.rm=TRUE))
}

print(vel_stats)

df_plot$Edge <- factor(df_plot$Edge, levels = rev(quantile_names))

# --- 4. Plot ---
edge_colors <- c("Leading (97.5%)" = "#0571b0", "75th Percentile" = "#92c5de", 
                 "Core Median (50%)" = "#999999", "25th Percentile" = "#f4a582", 
                 "Trailing (2.5%)" = "#ca0020")

p_margins <- ggplot(df_plot, aes(x = Year, y = Median, color = Edge, fill = Edge)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 1.5) +
  theme_minimal() +
  scale_color_manual(values = edge_colors) +
  scale_fill_manual(values = edge_colors) + # Fixed syntax here
  labs(title = "Spatio-Temporal Range Dynamics",
       subtitle = paste0("Leading Edge Velocity: ", round(vel_stats[["Leading (97.5%)"]][1], 2), " km/yr"),
       y = "Latitude (km)", 
       x = "Year")

print(p_margins)

ggsave("RangeDynamics.pdf", width = 8, height = 5)

# ============================================================================== #
# 13. LATENT ABUNDANCE & ANIMATIONS (UNIVERSAL BRIDGE VERSION) ####
# ============================================================================== #
library(gganimate)
library(gifski)

reconstruct_N_for_ppc <- function(draw, S, T_total, T_obs, temp_idx, t_scaled, from_idx, to_idx, E_disp, dists, thermal_gradient) {
  # Extract params for this draw
  phi <- as.numeric(draw[grep("phi\\[", names(draw))])
  gamma <- as.numeric(draw[grep("gamma_smooth\\[", names(draw))])
  eps_year <- as.numeric(draw[grep("eps_year\\[", names(draw))])
  r <- as.numeric(draw["r"])
  N0 <- as.numeric(draw["N0_proportion"])
  alpha <- as.numeric(draw["alpha"])
  
  # Calculate thermal niche and w_norm for this draw
  T_nodes <- c(as.numeric(draw["T_nodes[1]"]), as.numeric(draw["T_nodes[2]"]))
  thermal_niche <- as.numeric(draw["thermal_max"]) - 
    log(1 + exp(-as.numeric(draw["slope_L"]) * (thermal_gradient - T_nodes[1]))) - 
    log(1 + exp(-as.numeric(draw["slope_R"]) * (T_nodes[2] - thermal_gradient)))
  
  w <- exp(-dists / alpha)
  out_sum <- aggregate(w, by=list(from_idx), FUN=sum)$x
  w_norm <- w / out_sum[from_idx]
  
  # Simulation
  N_state_vec <- numeric(S * T_obs)
  R0 <- exp(r)
  R0m1 <- R0 - 1.0
  N_curr <- exp(thermal_niche[temp_idx[1:S]] + phi + (gamma * t_scaled[1]) + eps_year[1]) * N0
  
  for (t in 2:T_total) {
    logK <- thermal_niche[temp_idx[((t - 1) * S + 1) : (t * S)]] + phi + (gamma * t_scaled[t])
    K <- pmax(exp(logK), 1e-6)
    N_curr <- (R0 * N_curr) / (1.0 + (R0m1 / K) * N_curr)
    
    N_new <- numeric(S)
    for (e in 1:E_disp) N_new[to_idx[e]] <- N_new[to_idx[e]] + w_norm[e] * N_curr[from_idx[e]]
    N_curr <- N_new * exp(eps_year[t])
    N_curr[is.na(N_curr) | N_curr < 0] <- 0
    
    if (t > T_total - T_obs) {
      obs_t <- t - (T_total - T_obs)
      N_state_vec[((obs_t - 1) * S + 1):(obs_t * S)] <- N_curr
    }
  }
  return(N_state_vec)
}

cat("\n--- 13. Generating Latent Abundance Bridge ---\n")

# --- STEP 2: The Bridge Logic ---
cat("Detecting N_state or Reconstructing...\n")
draws_df <- as_draws_df(fit_final$draws())
all_vars <- colnames(draws_df)

if (any(grepl("N_state\\[", all_vars))) {
  cat("Detected N_state: Extracting median...\n")
  N_mat <- matrix(apply(as_draws_matrix(fit_final$draws("N_state")), 2, median), 
                  nrow = T_obs, ncol = S, byrow = TRUE)
} else {
  cat("N_state missing: Running R-reconstruction\n")
  n_reconst <- 200
  N_temp <- array(NA, dim = c(n_reconst, S * T_obs))
  for(i in 1:n_reconst) {
    N_temp[i,] <- reconstruct_N_for_ppc(draws_df[i,], S, T_total, T_obs, temp_idx, t_scaled, from_idx, to_idx, E_disp, dists, thermal_gradient)
  }
  N_mat <- matrix(apply(N_temp, 2, median), nrow = T_obs, ncol = S, byrow = TRUE)
}

cat("Bridge complete. N_mat is ready.\n")

# --- STEP 3: Generate Plot Data ---
coords <- as.data.frame(st_coordinates(st_centroid(my_grid_mainland)))
coords$grid_id <- valid_sites

df_n_mod <- data.frame(
  year = rep(obs_years, each = S),
  grid_id = rep(valid_sites, T_obs),
  abundance = as.vector(t(N_mat)) # Transpose if needed to match grid_id
) %>% left_join(coords, by = "grid_id")

# --- STEP 4: Animation ---
mod_cap <- quantile(df_n_mod$abundance, 0.99, na.rm=TRUE)

p_gif <- ggplot(df_n_mod, aes(X, Y, fill = pmin(abundance, mod_cap))) +
  geom_tile() +
  scale_fill_viridis_c(option = "magma", name = "Abundance (N)") +
  theme_void() +
  transition_states(year) +
  labs(title = "Latent Abundance Year: {closest_state}")

# Explicitly call the function using the namespace prefix
gganimate::animate(
  plot = p_gif, 
  renderer = gifski_renderer("latent_abundance.gif"), 
  width = 600, 
  height = 600
)
cat("GIF saved as 'latent_abundance.gif'\n")


# ============================================================================== #
# 14. RANGE EXPANSION MAP ####
# ============================================================================== #
cat("\n--- 14. Generating Latent Range Expansion Map ---\n")

library(ggplot2)
library(sf)
library(scales)
library(dplyr)

# 1. Calculate Change (Delta N)
# N_mat is [T_obs x S]. 
# N_End = Row T_obs, N_Start = Row 1
change_df <- data.frame(
  grid_id = valid_sites,
  N_Start = N_mat[1, ],
  N_End   = N_mat[nrow(N_mat), ],
  N_Change = N_mat[nrow(N_mat), ] - N_mat[1, ]
)

# 2. Join back to the spatial grid
# We filter the grid to match only the valid sites used in the model
grid_expansion <- my_grid_mainland %>%
  filter(grid_id %in% valid_sites) %>%
  left_join(change_df, by = "grid_id")

# 3. Define visualization limits (prevents outliers from washing out the map)
# We set limits to 95th percentile to keep the color scale sensitive
limit_val <- quantile(abs(grid_expansion$N_Change), 0.95, na.rm = TRUE)

p_expansion <- ggplot(grid_expansion) +
  geom_sf(aes(fill = N_Change), color = NA) +
  scale_fill_gradient2(
    low = "#d7191c",      # Red (Contraction/Loss)
    mid = "grey95",       # Neutral
    high = "#2c7bb6",     # Blue (Expansion/Gain)
    midpoint = 0, 
    limits = c(-limit_val, limit_val), 
    oob = scales::squish, # Squish values outside limits into the max color
    name = "Change in N\n(2000-2024)"
  ) +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(
    title = "Historical Latent Population Shift",
    subtitle = "Positive (Blue) indicates colonization/increase; Negative (Red) indicates decline",
    caption = "Values squished at 95th percentile to highlight spatial patterns"
  )

print(p_expansion)
ggsave("map_range_expansion.pdf", width = 6, height = 6)


# ============================================================================== #
# 14b. STATIC RANGE SNAPSHOTS ####
# ============================================================================== #
library(dplyr)
library(sf)

# 1. Ensure valid_sites and N_mat columns are aligned
# We assume N_mat columns correspond to 'valid_sites'
# If N_mat was created via the Bridge, it should be S columns wide.

target_years <- c(2000, 2008, 2016, 2024)
idx <- which(obs_years %in% target_years)

# 2. Build the snapshot data using an explicit loop (Alignment-Proof)
# This prevents the 'rep' recycling errors that cause striping
snapshot_list <- list()

for (i in seq_along(target_years)) {
  year_val <- target_years[i]
  mat_row  <- idx[i]
  
  snapshot_list[[i]] <- data.frame(
    grid_id = valid_sites, # This enforces the grid order
    year    = year_val,
    abundance = N_mat[mat_row, ] # This pulls the specific row for that year
  )
}

snapshot_data <- bind_rows(snapshot_list)

# 3. Explicit Spatial Join
# We join to the spatial grid using the grid_id, NOT position/index
grid_snapshots <- my_grid_mainland %>%
  filter(grid_id %in% valid_sites) %>%
  left_join(snapshot_data, by = "grid_id")

# 4. Check for NAs (If this is > 0, your alignment is still broken)
if(any(is.na(grid_snapshots$abundance))) {
  warning("Warning: Alignment issue detected! Some grid cells have no abundance data.")
}

# 5. Plot with consistent scaling
global_max <- quantile(snapshot_data$abundance, 0.99, na.rm = TRUE)

p_snapshots <- ggplot(grid_snapshots) +
  geom_sf(aes(fill = pmin(abundance, global_max)), color = NA) +
  facet_wrap(~year, ncol = 4) +
  scale_fill_viridis_c(
    option = "magma", 
    trans = "log1p", 
    breaks = c(0, 1, 2, 5, 10, 20, 50, 100),
    name = "Abundance (N)",
    guide = guide_colorbar(title.position = "top", 
                           title.hjust = 0.5, barwidth = unit(8, "cm"), 
                           barheight = unit(0.2, "cm"))
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "bottom"
  ) +
  labs(title = "Spatio-Temporal Range Expansion",
       subtitle = "Snapshot comparison of latent abundance")

print(p_snapshots)

ggsave("snapshots.pdf", width = 8, height = 5)

# ============================================================================== #
# 14c. PROPORTIONAL RANGE EXPANSION MAP ####
# ============================================================================== #
cat("\n--- 14. Generating Proportional Range Expansion Map (No limits/squish) ---\n")

# 1. Calculate Proportional Change
change_df <- data.frame(
  grid_id = valid_sites,
  Prop_Change = (N_mat[nrow(N_mat), ] - N_mat[1, ]) / (N_mat[1, ] + 1)
)

# 2. Join back to the spatial grid
grid_expansion <- my_grid_mainland %>%
  filter(grid_id %in% valid_sites) %>%
  left_join(change_df, by = "grid_id")

# 3. Plot without squishing (ggplot automatically uses the full data range)
p_prop_expansion <- ggplot(grid_expansion) +
  geom_sf(aes(fill = Prop_Change), color = NA) +
  scale_fill_gradient2(
    low = "#d7191c",      # Red (Relative Contraction)
    mid = "grey95",       # Neutral
    high = "#2c7bb6",     # Blue (Relative Expansion)
    midpoint = 0, 
    # By removing 'limits' and 'oob', ggplot will cover the full min to max
    labels = scales::percent, 
    name = "Rel. Change\n(2000-2024)"
  ) +
  theme_minimal() +
  labs(
    title = "Proportional Population Shift",
    subtitle = "Relative change in abundance (Full data range)"
  )

print(p_prop_expansion)
ggsave("map_prop_expansion_full.pdf", width = 6, height = 6)

# 3. Plot with fixed limits and outlier squishing
p_prop_expansion <- ggplot(grid_expansion) +
  geom_sf(aes(fill = Prop_Change), color = NA) +
  scale_fill_gradient2(
    low = "#d7191c",    # Red (Contraction)
    mid = "grey95",     # Neutral
    high = "#2c7bb6",   # Blue (Expansion)
    midpoint = 0, 
    limits = c(-0.7, 0.7),
    oob = scales::squish,        # Force values into the scale
    labels = scales::percent, 
    name = "Rel. Change\n(2000-2024)"
  ) +
  theme_minimal() +
  labs(
    title = "Proportional Population Shift",
    subtitle = "Relative change in abundance (Capped)" # Updated subtitle
  )

print(p_prop_expansion)
ggsave("map_prop_expansion_capped.pdf", width = 8, height = 8)


# ============================================================================== #
# 15. DISPERSAL KERNEL ####
# ============================================================================== #
library(ggplot2)
library(dplyr)
library(posterior)

# 1. Extract alpha draws (ensure you're using the final_draws object)
alpha_samples <- as_draws_df(subset_draws(final_draws, variable = "alpha"))$alpha

# 2. Define the distance range to plot (0 to 200km, matching your truncation)
dist_seq <- seq(0, 200, length.out = 100)

# 3. Calculate kernel density (Normalization for distance > 0)
kernel_matrix <- sapply(dist_seq, function(x) {
  (1 / alpha_samples) * exp(-x / alpha_samples)
})

# 4. Summarize (Median and 95% CI)
plot_df <- data.frame(
  Distance = dist_seq,
  Median = apply(kernel_matrix, 2, median),
  Lower = apply(kernel_matrix, 2, quantile, 0.025),
  Upper = apply(kernel_matrix, 2, quantile, 0.975)
)

# 5. Plot
ggplot(plot_df, aes(x = Distance)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "firebrick", alpha = 0.2) +
  geom_line(aes(y = Median), color = "firebrick", linewidth = 1.2) +
  geom_vline(xintercept = median(alpha_samples), linetype = "dashed") +
  annotate("text", x = median(alpha_samples) + 5, y = max(plot_df$Upper)*0.8, 
           label = paste("Mean (alpha) =", round(median(alpha_samples), 1), "km"), hjust = 0) +
  theme_minimal() +
  labs(title = "Estimated Dispersal Kernel",
       subtitle = "Probability density of movement distance (Laplace Kernel)",
       x = "Distance (km)",
       y = "Probability Density")
ggsave("Disp.pdf", width = 7, height = 4)



# ============================================================================== # #
# 17. FULL TREND LINES & UNCERTAINTY ANALYSIS ####
# ============================================================================== # #
library(sf); library(dplyr); library(tidyr); library(ggplot2); library(posterior); library(scales)

cat("\n--- Starting Trend Line Analysis (Bridge-Enabled) ---\n")

# --- 0. BRIDGE: Build Posterior Matrix ---
# This ensures we have the N_post_mat regardless of how the Stan model was saved
if(!exists("N_post_mat")) {
  all_vars <- variables(fit_final$draws())
  if ("N_state" %in% all_vars) {
    N_post_mat <- as_draws_matrix(fit_final$draws("N_state"))
  } else {
    cat("! N_state missing. Reconstructing draws from parameters (200 sample limit)...\n")
    n_draws_reconst <- 200 
    draws_df <- as_draws_df(fit_final$draws())
    N_post_mat <- matrix(NA, nrow = n_draws_reconst, ncol = S * T_obs)
    for(i in 1:n_draws_reconst) {
      N_post_mat[i, ] <- reconstruct_N_for_ppc(draws_df[i,], S, T_total, T_obs, temp_idx, t_scaled, from_idx, to_idx, E_disp, dists, thermal_gradient)
    }
  }
}

# Ensure Global Variables exist
N_sites <- S
T_obs   <- T_obs
obs_years <- 2000:(2000 + T_obs - 1)
N_mat   <- matrix(apply(N_post_mat, 2, median), nrow = T_obs, ncol = N_sites, byrow = TRUE)

# --- 1. PREP SPATIAL DATA ---
europe <- st_read("C:/Users/lieve/Documents/School/2025 - 2026/Thesis Statistics/Data/New/countries.gpkg", quiet = TRUE)
centroids <- st_centroid(my_grid_final)
europe <- st_transform(europe, st_crs(centroids))
nearest_idx <- st_nearest_feature(centroids, europe)

grid_country_map <- data.frame(
  grid_id = as.character(my_grid_final$grid_id),
  Country = europe$name[nearest_idx]
)

# --- 2. PREP ABUNDANCE DATA ---
df_n_long <- data.frame(
  year      = rep(obs_years, each = N_sites),
  grid_id   = as.character(rep(as.character(my_grid_final$grid_id), T_obs)),
  abundance = as.vector(t(N_mat))
) %>%
  left_join(grid_country_map, by = "grid_id") %>%
  filter(!is.na(Country))

# ============================================================================== # #
# 17a. AVERAGE N PER CELL (BY COUNTRY & EUROPE) ####
# ============================================================================== # #
country_trends <- df_n_long %>% 
  mutate(clean_abundance = ifelse(abundance < 0.1, 0, abundance)) %>%
  group_by(year, Country) %>%
  summarize(Mean_N = mean(clean_abundance, na.rm=TRUE), .groups = "drop") %>%
  group_by(Country) %>% filter(max(Mean_N) > 1) %>% ungroup()

# Sort legend
country_order <- country_trends %>% filter(year == max(year)) %>% arrange(desc(Mean_N)) %>% pull(Country)
country_trends <- country_trends %>% mutate(Country = factor(Country, levels = country_order))

overall_trend <- df_n_long %>%
  mutate(clean_abundance = ifelse(abundance < 0.1, 0, abundance)) %>%
  group_by(year) %>%
  summarize(Mean_N = mean(clean_abundance, na.rm=TRUE), .groups = "drop") %>%
  mutate(Country = "Total Europe (Average)")

p1 <- ggplot() +
  geom_line(data = country_trends, aes(x = year, y = Mean_N, color = Country), linewidth = 0.8) +
  geom_line(data = overall_trend, aes(x = year, y = Mean_N), linewidth = 1.2, color = "black", linetype = "dashed") +
  scale_color_viridis_d(option = "turbo") +
  theme_minimal() +
  labs(title = "Average Latent Abundance (N) per Grid Cell", x = "Year", y = "Average N")
print(p1)
ggsave("Trend_Avg_N_Country.pdf", width = 8, height = 5)

p1_faceted <- ggplot(country_trends, aes(x = year, y = Mean_N, color = Country)) +
  geom_line(linewidth = 1) + 
  # Facet by country: free_y allows each country to have its own vertical scale
  facet_wrap(~Country, scales = "free_y", drop = TRUE) + 
  scale_color_viridis_d(option = "turbo") +
  theme_minimal() + 
  theme(legend.position = "none", # Legend is redundant in facets
        strip.text = element_text(face = "bold", size = 10)) +
  labs(title = "Average Latent Abundance (N) per Grid Cell", 
       subtitle = "Faceted for independent trend analysis",
       x = "Year", y = "Average N")

print(p1_faceted)
ggsave("Trend_Avg_N_Faceted.pdf", width = 10, height = 8)

# ============================================================================== # #
# 17b. EUROPEAN AVERAGE TREND WITH 95% CI ####
# ============================================================================== # #
eu_avg_draws <- t(apply(N_post_mat, 1, function(row) {
  # Thresholding per draw
  row_mat <- matrix(ifelse(row < 0.1, 0, row), nrow = T_obs, ncol = N_sites, byrow = TRUE)
  rowMeans(row_mat)
}))

df_avg_ci <- data.frame(Year = obs_years, 
                        Median = apply(eu_avg_draws, 2, median),
                        Lower = apply(eu_avg_draws, 2, quantile, 0.025),
                        Upper = apply(eu_avg_draws, 2, quantile, 0.975))

p2 <- ggplot(df_avg_ci, aes(x = Year)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "dodgerblue4", alpha = 0.2) +
  geom_line(aes(y = Median), color = "dodgerblue4", linewidth = 1.5) +
  theme_minimal() +
  labs(title = "Overall European Expansion Trend (Avg N/cell)", x = "Year", y = "Average N")
print(p2)
ggsave("Trend_Europe_Avg_CI.pdf", width = 8, height = 5)

# ============================================================================== # #
# 17c. TOTAL ABUNDANCE PER COUNTRY ####
# ============================================================================== # #
# --- 1. CLEAN DATA & DROP LEVELS ---
country_totals <- df_n_long %>%
  mutate(clean_abundance = ifelse(abundance < 0.1, 0, abundance)) %>%
  # Filter NAs FIRST, before any other operations
  filter(!is.na(Country)) %>% 
  group_by(year, Country) %>%
  summarize(Total_N = sum(clean_abundance), .groups = "drop") %>%
  group_by(Country) %>% 
  filter(max(Total_N) > 10) %>% 
  ungroup() %>%
  # Crucial: This forgets that "NA" ever existed as a category
  droplevels() 

# --- 2. RE-CALCULATE ORDER ON CLEAN DATA ---
country_order <- country_totals %>%
  filter(year == max(year)) %>% 
  arrange(desc(Total_N)) %>% 
  pull(Country)

# --- 3. APPLY FACTOR ---
country_totals <- country_totals %>%
  mutate(Country = factor(Country, levels = country_order))

p3 <- ggplot(country_totals, aes(x = year, y = Total_N, color = Country)) +
  geom_line(linewidth = 1) + scale_color_viridis_d(option = "turbo") +
  theme_minimal() + scale_y_continuous(labels = scales::comma) +
  labs(title = "Total Latent Abundance (N) by Country", y = "Total Population (N)")
print(p3)
ggsave("Trend_Total_N_Country.pdf", width = 8, height = 5)

# Faceted Plot: Ideal for comparing trends across countries
p3_faceted <- ggplot(country_totals, aes(x = year, y = Total_N, color = Country)) +
  geom_line(linewidth = 1) + 
  # Added drop = TRUE to remove empty panels
  facet_wrap(~Country, scales = "free_y", drop = TRUE) + 
  scale_color_viridis_d(option = "turbo") +
  theme_minimal() +
  scale_y_continuous(labels = scales::comma) + 
  theme(legend.position = "none",  # Legend is redundant when faceted
        strip.text = element_text(face = "bold")) +
  labs(title = "Total Latent Abundance (N) by Country", 
       subtitle = "Faceted for independent scale comparison",
       x = "Year", y = "Total Population (N)")

print(p3_faceted)
ggsave("Trend_Total_N_Faceted.pdf", width = 10, height = 8)

# ============================================================================== # #
# 17c. TOTAL ABUNDANCE PER COUNTRY (SORTED + 95% CI) ####
# ============================================================================== # #
cat("\n--- Calculating CI for Total Abundance by Country ---\n")

# 1. Create a lookup table to map columns of N_post_mat to Countries
# This tells R which columns in the matrix belong to which country
col_lookup <- data.frame(
  col_idx = 1:(N_sites * T_obs),
  year    = rep(obs_years, each = N_sites),
  grid_id = rep(as.character(my_grid_final$grid_id), T_obs)
) %>% 
  left_join(grid_country_map, by = "grid_id") %>% 
  filter(!is.na(Country))

# 2. Iterate through posterior draws (Memory efficient)
# We calculate total N for each country/year for every single MCMC draw
n_draws <- nrow(N_post_mat)
all_draws_list <- list()

for(i in 1:n_draws) {
  # Apply threshold to this specific draw
  draw_values <- ifelse(N_post_mat[i, ] < 0.01, 0, N_post_mat[i, ])
  
  # Map draw to country/year
  temp_df <- data.frame(val = draw_values, col_idx = 1:(N_sites * T_obs)) %>%
    inner_join(col_lookup, by = "col_idx") %>%
    group_by(Country, year) %>%
    summarize(Total_N = sum(val), .groups = "drop")
  
  all_draws_list[[i]] <- temp_df
}

# 3. Summarize Quantiles (The CI values)
ci_data <- bind_rows(all_draws_list) %>%
  group_by(Country, year) %>%
  summarize(
    Median = median(Total_N),
    Lower = quantile(Total_N, 0.025),
    Upper = quantile(Total_N, 0.975),
    .groups = "drop"
  ) %>%
  # Keep only substantial populations for the plot
  group_by(Country) %>% filter(max(Median) > 1) %>% ungroup() %>%
  droplevels()

# 4. Apply the Sort Order (based on Median in the final year)
sort_order <- ci_data %>%
  filter(year == max(year)) %>% 
  arrange(desc(Median)) %>% 
  pull(Country)

ci_data <- ci_data %>% mutate(Country = factor(Country, levels = sort_order))

# 5. Faceted Plot with CI Ribbon
p3_faceted_ci <- ggplot(ci_data, aes(x = year)) +
  # Add the CI Ribbon
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "seagreen4", alpha = 0.2) +
  # Add the Median line
  geom_line(aes(y = Median), color = "seagreen4", linewidth = 1) +
  facet_wrap(~Country, scales = "free_y", drop = TRUE) +
  theme_minimal() +
  scale_y_continuous(labels = scales::comma, limits = c(0, NA)) +
  theme(strip.text = element_text(face = "bold", size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Total Latent Abundance (N) by Country (Sorted)",
       subtitle = "Ribbon indicates 95% Credible Interval from posterior draws",
       x = "Year", y = "Total Population (N)")

print(p3_faceted_ci)
ggsave("Trend_Total_N_Faceted_CI.pdf", width = 11, height = 9)

# ============================================================================== # #
# 17d. TOTAL EUROPEAN ABUNDANCE WITH 95% CI ####
# ============================================================================== # #
eu_tot_draws <- t(apply(N_post_mat, 1, function(row) {
  row_mat <- matrix(ifelse(row < 0.1, 0, row), nrow = T_obs, ncol = N_sites, byrow = TRUE)
  rowSums(row_mat)
}))

df_tot_ci <- data.frame(Year = obs_years, 
                        Median = apply(eu_tot_draws, 2, median),
                        Lower = apply(eu_tot_draws, 2, quantile, 0.025),
                        Upper = apply(eu_tot_draws, 2, quantile, 0.975))

p4 <- ggplot(df_tot_ci, aes(x = Year)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "seagreen4", alpha = 0.2) +
  geom_line(aes(y = Median), color = "seagreen4", linewidth = 1.5) +
  scale_y_continuous(labels = scales::comma) +
  theme_minimal() +
  labs(title = "Total European Population Trend", y = "Total Population (Sum of N)")
print(p4)
ggsave("Trend_Europe_Total_CI.pdf", width = 8, height = 5)


# ============================================================================== # #
# VISUALIZING THE SPATIAL FIELD (GAUSSIAN PROCESS) ####
# ============================================================================== # #
cat("\n--- Extracting and Visualizing the Spatial Field (phi) ---\n")
library(ggplot2)
library(sf)
library(posterior)

# 1. Extract the posterior draws for the spatial vector 'phi'
# 'phi' contains the spatial offset for each of the 1,835 sites
phi_draws <- as_draws_matrix(fit_final$draws("phi"))

# 2. Calculate the median spatial effect for each grid cell
phi_median <- apply(phi_draws, 2, median)

# 3. Attach the spatial effects to your pre-sorted sf object
# (my_grid_mainland is already perfectly aligned with the Stan indices)
my_grid_spatial <- my_grid_mainland
my_grid_spatial$spatial_effect <- phi_median

# 4. Plot the Median Spatial Field
p_spatial <- ggplot(my_grid_spatial) +
  # Using color = NA removes grid borders for a smooth, continuous look
  geom_sf(aes(fill = spatial_effect), color = NA) + 
  # 'magma' or 'viridis' are great for continuous spatial fields
  scale_fill_viridis_c(option = "magma", name = "Spatial\nEffect (φ)") + 
  theme_void() +
  labs(
    title = "Latent Spatial Field (Gaussian Process)",
    subtitle = "Intrinsic landscape suitability for Sympetrum danae (Independent of Climate)",
    caption = "Lighter colors = Intrinsic Hotspots | Darker colors = Intrinsic Coldspots"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, margin = margin(b = 15)),
    legend.position = "right"
  )

print(p_spatial)

cat("\n--- Visualizing Spatial Intercept (phi) and Spatial Trend (gamma) ---\n")

# 1. Extract posterior draws for both fields
phi_draws   <- as_draws_matrix(fit_final$draws("phi"))
gamma_draws <- as_draws_matrix(fit_final$draws("gamma_smooth"))

# 2. Calculate medians
my_grid_spatial$phi_median   <- apply(phi_draws, 2, median)
my_grid_spatial$gamma_median <- apply(gamma_draws, 2, median)

# 3. Plot Spatial Intercept (phi)
p_phi <- ggplot(my_grid_spatial) +
  geom_sf(aes(fill = phi_median), color = NA) + 
  scale_fill_viridis_c(option = "magma", name = "Intrinsic\nSuitability (phi)",
                       guide = guide_colorbar(title.position = "top", 
                                              title.hjust = 0.5, barwidth = unit(0.2, "cm"), 
                                              barheight = unit(4, "cm"))) + 
  theme_void()

# 4. Plot Spatial Trend (gamma)
# We use a divergent scale (RdBu) to show positive vs negative trends
p_gamma <- ggplot(my_grid_spatial) +
  geom_sf(aes(fill = gamma_median), color = NA) + 
  scale_fill_gradient2(low = "firebrick", mid = "white", high = "seagreen4", 
                       midpoint = 0, name = "Spatial trend\n (gamma)",
                       guide = guide_colorbar(title.position = "top", 
                                              title.hjust = 0.5, barwidth = unit(0.2, "cm"), 
                                              barheight = unit(4, "cm"))) + 
  theme_void()

# Display both
print(p_phi)
ggsave("phi.pdf", width = 5, height = 6)
print(p_gamma) 
ggsave("gamma.pdf", width = 5, height = 6)

library(patchwork)
# 1. Combine the plots
combined_plot <- p_phi + p_gamma +
  plot_layout(ncol = 2) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(face = "bold", size = 16)
    )
  )

# 2. Print to screen to verify
print(combined_plot)

# 3. Save as a single PDF
ggsave("phi_gamma.pdf", plot = combined_plot, width = 8, height = 5)

cat("\n--- Fixing Dimension Mismatch: Slicing to T_obs ---\n")

# We only care about the years we actually have observations for (25 years)
# We slice the eps_draws matrix to only keep the first T_obs columns
eps_draws_subset <- draws_df %>% 
  select(starts_with("eps_year")) %>% 
  select(1:all_of(T_obs))

# 1. Clean the subset by removing metadata columns (those starting with ".")
eps_draws_subset_clean <- eps_draws_subset %>% 
  select(-starts_with("."))

# 2. Now calculate the stats using only the clean subset
eps_summary <- data.frame(
  year = obs_years, 
  median = apply(eps_draws_subset_clean, 2, median),
  lower = apply(eps_draws_subset_clean, 2, quantile, probs = 0.025),
  upper = apply(eps_draws_subset_clean, 2, quantile, probs = 0.975)
)

# Verify the result
print(eps_summary)

# plot
p_year <- ggplot(eps_summary, aes(x = year, y = median)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "firebrick", alpha = 0.2) +
  geom_line(color = "firebrick", linewidth = 1) +
  geom_point(color = "firebrick") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  theme_minimal() +
  labs(
    title = "Global Inter-annual Shocks (ε_t)",
    subtitle = "Annual process deviations affecting population size across all sites",
    x = "Year", y = "Log-scale Deviation"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    panel.grid.minor = element_blank()
  )

print(p_year)
ggsave("epsilon.pdf", width = 7, height = 4)


# ============================================================================== # #
# 18. FORECASTING (SSP2-4.5) - UNIVERSAL BRIDGE & DGP ALIGNMENT ####
# ============================================================================== # #
library(geodata); library(terra); library(sf); library(dplyr); library(ggplot2); library(posterior)

# ------------------------------------------------------------------------------ #
# 1. BRIDGE: DETECT OR RECONSTRUCT 2024 STATE
# ------------------------------------------------------------------------------ #
cat("\n--- 1. Generating Hot-Start Bridge (Year 2024) ---\n")

draws_df <- as_draws_df(fit_final$draws())
all_vars <- colnames(draws_df)

if (any(grepl("N_state\\[", all_vars))) {
  cat("Detected N_state: Extracting median state for 2024...\n")
  N_post_mat <- as_draws_matrix(fit_final$draws("N_state"))
  # N_mat_full dimensions: [Samples, S*T_obs]
  N_mat_full <- matrix(apply(N_post_mat, 2, median), nrow = T_obs, ncol = S, byrow = TRUE)
  N_start_2024 <- N_mat_full[T_obs, ] # Last observed year
  
} else {
  cat("N_state missing: Running R-reconstruction for hot-start (this may take a moment)...\n")
  
  # 1. Define number of draws to use
  n_draws <- 200
  reconst_matrix <- matrix(NA, nrow = n_draws, ncol = S * T_obs)
  
  # 2. Iterate through rows one by one
  for (i in 1:n_draws) {
    # Extract a single row as a named vector (draws_df is a tibble, so row is a dataframe)
    # We unlist it to get the named vector the function expects
    single_draw_vec <- unlist(draws_df[i, ])
    
    reconst_matrix[i, ] <- reconstruct_N_for_ppc(
      draw = single_draw_vec, 
      S = S, T_total = T_total, T_obs = T_obs, 
      temp_idx = temp_idx, t_scaled = t_scaled, 
      from_idx = from_idx, to_idx = to_idx, E_disp = E_disp, 
      dists = dists, thermal_gradient = thermal_gradient
    )
  }
  
  # 3. Take median across the rows
  N_mat_reconstructed <- matrix(apply(reconst_matrix, 2, median), nrow = T_obs, ncol = S, byrow = TRUE)
  N_start_2024 <- N_mat_reconstructed[T_obs, ]
}

# ------------------------------------------------------------------------------ #
# 2. CLIMATE DATA LOADING & RAMP BUILDING
# ------------------------------------------------------------------------------ #
cat("\n--- 2. Loading Future Climate Data ---\n")
climate_dir <- "C:/Users/lieve/Documents/School/2025 - 2026/Thesis Statistics/Data/Climate_Rasters"
future_temp_raster <- cmip6_world(model = "MPI-ESM1-2-HR", ssp = "245", time = "2041-2060", var = "bioc", res = 5, path = climate_dir)[[1]]
current_temp_raster <- worldclim_global(var = "bio", res = 5, path = climate_dir)[[1]]

grid_centroids <- st_centroid(my_grid_final)
future_temp_vec <- terra::extract(future_temp_raster, vect(grid_centroids))[, 2]
current_temp_vec <- terra::extract(current_temp_raster, vect(grid_centroids))[, 2]
future_temp_vec[is.na(future_temp_vec)] <- mean(future_temp_vec, na.rm = TRUE)
current_temp_vec[is.na(current_temp_vec)] <- mean(current_temp_vec, na.rm = TRUE)

cat("\n--- Recovering missing historical_temp_matrix ---\n")

# Ensure the dependencies exist
if(!exists("template_obs") || !exists("valid_sites")) {
  stop("Missing 'template_obs' or 'valid_sites'. Please load your primary observation data first.")
}

# Rebuild the matrix
# 1. Cast LONG to WIDE format
wide_temps <- dcast(template_obs, grid_id ~ year, value.var = "mean_temp")

# 2. Match the row order to your valid_sites (this order is critical!)
wide_temps <- wide_temps[match(valid_sites, wide_temps$grid_id), ]

# 3. Convert to a pure numeric matrix
historical_temp_matrix <- as.matrix(wide_temps[, -1, with = FALSE])

cat("Successfully rebuilt 'historical_temp_matrix'. Dimensions:", 
    nrow(historical_temp_matrix), "sites x", ncol(historical_temp_matrix), "years.\n")

# Build Ramp (20 years)
warming_delta <- future_temp_vec - current_temp_vec
future_temp_matrix_smooth <- sapply(1:20, function(t) historical_temp_matrix[, 25] + ((t/20) * warming_delta))

# Indices for Stan Niche lookup
T_min <- min(stan_data$thermal_gradient)
N_thermal_gradient <- length(stan_data$thermal_gradient)
future_temp_idx <- apply(future_temp_matrix_smooth, 2, function(x) pmax(1, pmin(round((x - T_min) * 10) + 1, N_thermal_gradient)))


cat("\n--- Plotting Temperature Ramp ---\n")

# 1. Prepare Historical Mean Temps (Across all sites)
hist_temp_means <- colMeans(historical_temp_matrix, na.rm = TRUE)
hist_temp_df <- data.frame(
  Year = 2000:(2000 + ncol(historical_temp_matrix) - 1),
  Temp = hist_temp_means,
  Type = "Historical"
)

# 2. Prepare Forecast Mean Temps (Across all sites)
forecast_temp_means <- colMeans(future_temp_matrix_smooth, na.rm = TRUE)
forecast_temp_df <- data.frame(
  Year = 2025:2044,
  Temp = forecast_temp_means,
  Type = "Forecast (SSP2-4.5)"
)

# 3. Combine and Plot
temp_plot_data <- bind_rows(hist_temp_df, forecast_temp_df)

ggplot(temp_plot_data, aes(x = Year, y = Temp, color = Type)) +
  geom_line(linewidth = 1.2) +
  geom_point() +
  geom_vline(xintercept = 2024, linetype = "dashed", color = "grey50") +
  theme_minimal() +
  scale_color_manual(values = c("Historical" = "black", "Forecast (SSP2-4.5)" = "firebrick")) +
  labs(
    title = "Temperature Forcing: Historical Baseline & Future Projection",
    subtitle = "Aggregated mean annual temperature across all study sites",
    x = "Year", y = "Mean Temperature (°C)"
  ) +
  theme(legend.position = "bottom")

ggsave("TempForecast.pdf", width = 7, height = 5)

# ------------------------------------------------------------------------------ #
# 3. RUNNING STOCHASTIC FORECAST (50 iterations)
# ------------------------------------------------------------------------------ #
cat("\n--- 3. Running Stochastic Forecast ---\n")

# Define the simulation function
simulate_abundance_DGP <- function(S, T_total, thermal_niche, temp_idx, 
                                   phi, gamma_smooth, t_scaled, eps_year, 
                                   r, N0_proportion, w_norm, from_idx, to_idx, E_disp,
                                   N_init, max_K_cap = 5000) {
  
  R0 <- exp(r); R0m1 <- R0 - 1.0
  N_curr <- N_init 
  trajectory <- matrix(NA, nrow = S, ncol = T_total)
  
  for (t in 1:T_total) {
    # Extract indices for current year
    idx_range <- ((t - 1) * S + 1):(t * S)
    
    # Calculate Capacity K with safety cap
    logK <- thermal_niche[temp_idx[idx_range]] + phi + (gamma_smooth * t_scaled[t])
    K <- pmin(exp(logK), max_K_cap)
    
    # Beverton-Holt Growth
    N_curr <- (R0 * N_curr) / (1.0 + (R0m1 / (K + 1e-9)) * N_curr)
    
    # Dispersal
    N_new <- numeric(S)
    for (e in 1:E_disp) {
      N_new[to_idx[e]] <- N_new[to_idx[e]] + w_norm[e] * N_curr[from_idx[e]]
    }
    
    # Process Error
    N_curr <- N_new * exp(eps_year[t])
    trajectory[, t] <- N_curr
  }
  return(trajectory)
}

cat("\n--- Extracting posterior parameters to R workspace ---\n")

# Extract parameter draws
draws_df <- as_draws_df(fit_final$draws())

# Extract point estimates (medians) for the parameters
r_est            <- median(draws_df$r)
N0_prop_est      <- median(draws_df$N0_proportion)
alpha_est        <- median(exp(draws_df$log_alpha))

# Extract vectors/matrices
phi_est          <- colMeans(as_draws_matrix(fit_final$draws("phi")))
gamma_smooth_est <- colMeans(as_draws_matrix(fit_final$draws("gamma_smooth")))

# Re-calculate dispersal weights (w_norm) using the extracted alpha
w <- exp(-stan_data$dists / alpha_est)
out_sum <- rep(1e-12, S)
for (e in 1:stan_data$E_disp) {
  out_sum[stan_data$from_idx[e]] <- out_sum[stan_data$from_idx[e]] + w[e]
}
w_norm_est <- w / out_sum[stan_data$from_idx]

cat("\n--- Calculating thermal_niche_est ---\n")

# Reconstruct the niche using the same formula as your Stan 'transformed parameters'
# We access the T_nodes safely using the column names generated by posterior
thermal_niche_est <- median(draws_df$thermal_max) - 
  log(1 + exp(-median(draws_df$slope_L) * (stan_data$thermal_gradient - median(draws_df[["T_nodes[1]"]])))) - 
  log(1 + exp(-median(draws_df$slope_R) * (median(draws_df[["T_nodes[2]"]]) - stan_data$thermal_gradient)))

cat("thermal_niche_est defined. Length:", length(thermal_niche_est), "\n")

cat("All parameters extracted. Ready for simulation.\n")

n_sims <- 50
all_sims_list <- lapply(1:n_sims, function(i) {
  sim_res <- simulate_abundance_DGP(
    S = S, T_total = 20, thermal_niche = thermal_niche_est,
    temp_idx = as.vector(future_temp_idx),
    phi = phi_est, gamma_smooth = gamma_smooth_est,
    t_scaled = seq(max(stan_data$t_scaled)+1, length.out=20),
    eps_year = rep(0, 20), 
    r = r_est,
    N0_proportion = N0_prop_est, 
    w_norm = w_norm_est,
    from_idx = stan_data$from_idx, 
    to_idx = stan_data$to_idx, 
    E_disp = stan_data$E_disp,
    N_init = N_start_2024
  )
  data.frame(year = rep(2025:2044, each = S), 
             grid_id = rep(my_grid_final$grid_id, 20), 
             abundance = as.vector(sim_res), 
             sim_id = i)
})


# Extract sigma_year from your draws dataframe
sigma_year_est <- median(draws_df$sigma_year)

# Verification
cat("sigma_year_est is now defined as:", sigma_year_est, "\n")

n_sims <- 50
all_sims_list <- lapply(1:n_sims, function(i) {
  
  # Generate random process error using the extracted sigma_year_est
  sim_eps <- rnorm(20, mean = 0, sd = sigma_year_est)
  
  sim_res <- simulate_abundance_DGP(
    S = S, T_total = 20, thermal_niche = thermal_niche_est,
    temp_idx = as.vector(future_temp_idx),
    phi = phi_est, gamma_smooth = gamma_smooth_est,
    t_scaled = seq(max(stan_data$t_scaled)+1, length.out=20),
    eps_year = sim_eps, # <--- Now this will work!
    r = r_est,
    N0_proportion = N0_prop_est, 
    w_norm = w_norm_est,
    from_idx = stan_data$from_idx, 
    to_idx = stan_data$to_idx, 
    E_disp = stan_data$E_disp,
    N_init = N_start_2024
  )
  
  data.frame(
    year = rep(2025:2044, each = S), 
    grid_id = rep(my_grid_final$grid_id, 20), 
    abundance = as.vector(sim_res), 
    sim_id = i
  )
})



# 4. UNIFIED PLOTTING (Historical + Forecast - With Stitching)
cat("\n--- Final Alignment: Stitching Forecast to 2024 ---\n")

# 1. PREPARE HISTORICAL DATA (Up to 2024)
hist_data <- df_n_long %>% 
  mutate(Country = as.character(Country)) %>%
  filter(!is.na(Country), Country != "NA", Country != "Unassigned") %>%
  group_by(Country, year) %>%
  summarize(Median = sum(abundance), Lower = NA, Upper = NA, Type = "Historical", .groups="drop") %>%
  rename(Year = year)

# 2. PREPARE FORECAST DATA (Add 2024 as the 'Bridge' point)
# We extract the median 2024 state (N_start_2024) to ensure the line starts at the right height
forecast_bridge <- data.frame(
  Country = grid_country_map$Country,
  Year = 2024,
  Median = N_start_2024, # The exact state from 2024
  Lower = N_start_2024,
  Upper = N_start_2024,
  Type = "Forecast"
) %>%
  filter(!is.na(Country), Country != "NA", Country != "Unassigned") %>%
  group_by(Country, Year, Type) %>%
  summarize(Median = sum(Median), Lower = sum(Lower), Upper = sum(Upper), .groups = "drop")

# 3. COMBINE
forecast_main <- bind_rows(all_sims_list) %>%
  left_join(grid_country_map, by = "grid_id") %>%
  mutate(Country = as.character(Country)) %>%
  filter(!is.na(Country), Country != "NA", Country != "Unassigned") %>%
  group_by(Country, year, sim_id) %>% summarize(Total_N = sum(abundance), .groups="drop") %>%
  group_by(Country, year) %>%
  summarize(
    Median = median(Total_N),
    Lower = quantile(Total_N, 0.025),
    Upper = quantile(Total_N, 0.975),
    Type = "Forecast", .groups="drop"
  ) %>%
  rename(Year = year)

forecast_final <- bind_rows(forecast_bridge, forecast_main)

# 1. COMBINE AND CLEAN IN ONE PIPELINE
# We do not use the old 'country_order' vector because it is 'poisoned' with NA
full_plot_data <- bind_rows(hist_data, forecast_final) %>%
  mutate(Country = as.character(Country)) %>%
  # Aggressive filtering of anything that isn't a real country
  filter(!is.na(Country), 
         Country != "NA", 
         Country != "Unassigned", 
         Country != "") %>%
  # Drop empty levels immediately
  droplevels()

# 2. DYNAMIC SORTING
# Instead of using an external vector, we sort based on the remaining valid data
# This sorts by total population in the last year
last_year_data <- full_plot_data %>% 
  filter(Year == 2044) %>% 
  arrange(desc(Median))

full_plot_data$Country <- factor(full_plot_data$Country, levels = last_year_data$Country)

# 3. PLOT
ggplot(full_plot_data, aes(x = Year)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Type), alpha = 0.2) +
  geom_line(aes(y = Median, color = Type), linewidth = 1) +
  geom_vline(xintercept = 2024, linetype = "dashed", color = "grey50") +
  # Use your existing color scales
  scale_color_manual(values = c("Historical" = "grey40", "Forecast" = "seagreen4")) +
  scale_fill_manual(values = c("Historical" = "transparent", "Forecast" = "seagreen4")) +
  facet_wrap(~Country, scales = "free_y", drop = TRUE) +
  theme_minimal() +
  labs(
    title = "Historical & Forecasted Total Population (N) by Country",
    x = "Year", y = "Total Population (N)"
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


## TAKING INTO ACCOUNT PARAMETER UNCERTAINTY ####

cat("\n--- Running Full-Uncertainty Forecast (Propagating Parameters + Process Error) ---\n")

n_sims <- 200
all_sims_list <- lapply(1:n_sims, function(i) {
  
  # 1. SAMPLE PARAMETERS: Pick a random row from your posterior draws
  # This propagates the uncertainty of r, phi, thermal_niche, etc.
  draw_idx <- sample(1:nrow(draws_df), 1)
  draw <- draws_df[draw_idx, ]
  
  # Extract parameters for this specific draw
  curr_phi <- as.numeric(draw[grep("^phi\\[", names(draw))])
  curr_gamma <- as.numeric(draw[grep("^gamma_smooth\\[", names(draw))])
  curr_r <- as.numeric(draw$r)
  curr_N0 <- as.numeric(draw$N0_proportion)
  curr_alpha <- as.numeric(exp(draw$log_alpha))
  
  # Reconstruct Niche for THIS specific draw
  curr_niche <- as.numeric(draw$thermal_max) - 
    log(1 + exp(-as.numeric(draw$slope_L) * (stan_data$thermal_gradient - as.numeric(draw[["T_nodes[1]"]])))) - 
    log(1 + exp(-as.numeric(draw$slope_R) * (as.numeric(draw[["T_nodes[2]"]]) - stan_data$thermal_gradient)))
  
  # Recompute weights (alpha is different for every draw!)
  curr_w <- exp(-stan_data$dists / curr_alpha)
  curr_out_sum <- rep(1e-12, S)
  for (e in 1:stan_data$E_disp) curr_out_sum[stan_data$from_idx[e]] <- curr_out_sum[stan_data$from_idx[e]] + curr_w[e]
  curr_w_norm <- curr_w / curr_out_sum[stan_data$from_idx]
  
  # 2. GENERATE PROCESS ERROR
  # We use the sigma_year from this specific draw
  curr_sigma <- as.numeric(draw$sigma_year)
  sim_eps <- rnorm(20, mean = 0, sd = curr_sigma)
  
  # 3. RUN SIMULATION
  sim_res <- simulate_abundance_DGP(
    S = S, T_total = 20, thermal_niche = curr_niche,
    temp_idx = as.vector(future_temp_idx),
    phi = curr_phi, gamma_smooth = curr_gamma,
    t_scaled = seq(max(stan_data$t_scaled)+1, length.out=20),
    eps_year = sim_eps,
    r = curr_r,
    N0_proportion = curr_N0, 
    w_norm = curr_w_norm,
    from_idx = stan_data$from_idx, 
    to_idx = stan_data$to_idx, 
    E_disp = stan_data$E_disp,
    N_init = N_start_2024 # Note: You might technically want to draw N_init from the posterior too
  )
  
  data.frame(year = rep(2025:2044, each = S), 
             grid_id = rep(my_grid_final$grid_id, 20), 
             abundance = as.vector(sim_res), 
             sim_id = i)
})



# PLOTTING (Historical + Forecast - With Stitching)
cat("\n--- Final Alignment: Stitching Forecast to 2024 ---\n")

# 1. PREPARE HISTORICAL DATA (Up to 2024)
hist_data <- df_n_long %>% 
  mutate(Country = as.character(Country)) %>%
  filter(!is.na(Country), Country != "NA", Country != "Unassigned") %>%
  group_by(Country, year) %>%
  summarize(Median = sum(abundance), Lower = NA, Upper = NA, Type = "Historical", .groups="drop") %>%
  rename(Year = year)

# 2. PREPARE FORECAST DATA (Add 2024 as the 'Bridge' point)
# We extract the median 2024 state (N_start_2024) to ensure the line starts at the right height
forecast_bridge <- data.frame(
  Country = grid_country_map$Country,
  Year = 2024,
  Median = N_start_2024, # The exact state from 2024
  Lower = N_start_2024,
  Upper = N_start_2024,
  Type = "Forecast"
) %>%
  filter(!is.na(Country), Country != "NA", Country != "Unassigned") %>%
  group_by(Country, Year, Type) %>%
  summarize(Median = sum(Median), Lower = sum(Lower), Upper = sum(Upper), .groups = "drop")

# 3. COMBINE
forecast_main <- bind_rows(all_sims_list) %>%
  left_join(grid_country_map, by = "grid_id") %>%
  mutate(Country = as.character(Country)) %>%
  filter(!is.na(Country), Country != "NA", Country != "Unassigned") %>%
  group_by(Country, year, sim_id) %>% summarize(Total_N = sum(abundance), .groups="drop") %>%
  group_by(Country, year) %>%
  summarize(
    Median = median(Total_N),
    Lower = quantile(Total_N, 0.025),
    Upper = quantile(Total_N, 0.975),
    Type = "Forecast", .groups="drop"
  ) %>%
  rename(Year = year)

forecast_final <- bind_rows(forecast_bridge, forecast_main)

# 1. COMBINE AND CLEAN IN ONE PIPELINE
# We do not use the old 'country_order' vector because it is 'poisoned' with NA
full_plot_data <- bind_rows(hist_data, forecast_final) %>%
  mutate(Country = as.character(Country)) %>%
  # Aggressive filtering of anything that isn't a real country
  filter(!is.na(Country), 
         Country != "NA", 
         Country != "Unassigned", 
         Country != "") %>%
  # Drop empty levels immediately
  droplevels()

# 2. DYNAMIC SORTING
# Instead of using an external vector, we sort based on the remaining valid data
# This sorts by total population in the last year
last_year_data <- full_plot_data %>% 
  filter(Year == 2044) %>% 
  arrange(desc(Median))

full_plot_data$Country <- factor(full_plot_data$Country, levels = last_year_data$Country)

# 3. PLOT
ggplot(full_plot_data, aes(x = Year)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Type), alpha = 0.2) +
  geom_line(aes(y = Median, color = Type), linewidth = 1) +
  geom_vline(xintercept = 2024, linetype = "dashed", color = "grey50") +
  # Use your existing color scales
  scale_color_manual(values = c("Historical" = "grey40", "Forecast" = "seagreen4")) +
  scale_fill_manual(values = c("Historical" = "transparent", "Forecast" = "seagreen4")) +
  facet_wrap(~Country, scales = "free_y", drop = TRUE) +
  theme_minimal() +
  labs(
    title = "Historical & Forecasted Total Population (N) by Country",
    x = "Year", y = "Total Population (N)"
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )




cat("\n--- Running Counterfactual Forecast (No Gamma Trend) ---\n")

n_sims <- 200
no_trend_sims_list <- lapply(1:n_sims, function(i) {
  
  # 1. SAMPLE PARAMETERS
  draw_idx <- sample(1:nrow(draws_df), 1)
  draw <- draws_df[draw_idx, ]
  
  # Base parameters
  curr_phi <- as.numeric(draw[grep("^phi\\[", names(draw))])
  curr_r <- as.numeric(draw$r)
  curr_N0 <- as.numeric(draw$N0_proportion)
  curr_alpha <- as.numeric(exp(draw$log_alpha))
  
  # KEY CHANGE: Zero out the gamma trend (No climate trend)
  # This makes the environmental effect constant over time
  curr_gamma_no_trend <- rep(0, length(grep("^gamma_smooth\\[", names(draw)))) 
  
  # Reconstruct Niche
  curr_niche <- as.numeric(draw$thermal_max) - 
    log(1 + exp(-as.numeric(draw$slope_L) * (stan_data$thermal_gradient - as.numeric(draw[["T_nodes[1]"]])))) - 
    log(1 + exp(-as.numeric(draw$slope_R) * (as.numeric(draw[["T_nodes[2]"]]) - stan_data$thermal_gradient)))
  
  # Recompute weights
  curr_w <- exp(-stan_data$dists / curr_alpha)
  curr_out_sum <- rep(1e-12, S)
  for (e in 1:stan_data$E_disp) curr_out_sum[stan_data$from_idx[e]] <- curr_out_sum[stan_data$from_idx[e]] + curr_w[e]
  curr_w_norm <- curr_w / curr_out_sum[stan_data$from_idx]
  
  # 2. RUN SIMULATION
  sim_res <- simulate_abundance_DGP(
    S = S, T_total = 20, thermal_niche = curr_niche,
    temp_idx = as.vector(future_temp_idx),
    phi = curr_phi, 
    gamma_smooth = curr_gamma_no_trend, # Using zero-trend
    t_scaled = seq(max(stan_data$t_scaled)+1, length.out=20),
    eps_year = rnorm(20, 0, as.numeric(draw$sigma_year)),
    r = curr_r,
    N0_proportion = curr_N0, 
    w_norm = curr_w_norm,
    from_idx = stan_data$from_idx, 
    to_idx = stan_data$to_idx, 
    E_disp = stan_data$E_disp,
    N_init = N_start_2024
  )
  
  data.frame(year = rep(2025:2044, each = S), 
             grid_id = rep(my_grid_final$grid_id, 20), 
             abundance = as.vector(sim_res), 
             sim_id = i,
             Scenario = "No-Trend Forecast")
})

# Combine the two datasets
forecast_main_trend <- bind_rows(all_sims_list) %>% mutate(Scenario = "Climate-Trend Forecast")
forecast_main_notrend <- bind_rows(no_trend_sims_list)

cat("\n--- Standardizing Schemas for Comparison Plot ---\n")

# 1. PREPARE BRIDGE
# Already has Country and Median
bridge_clean <- forecast_bridge %>% 
  rename(Median = Median) %>% 
  mutate(Scenario = "Bridge")

# 2. PREPARE TREND FORECAST
trend_clean <- bind_rows(all_sims_list) %>%
  left_join(grid_country_map, by = "grid_id") %>%
  mutate(Country = as.character(Country), Scenario = "Climate-Trend") %>%
  group_by(Country, year, Scenario, sim_id) %>% summarize(Total_N = sum(abundance), .groups="drop") %>%
  group_by(Country, year, Scenario) %>%
  summarize(Median = median(Total_N), Lower = quantile(Total_N, 0.025), Upper = quantile(Total_N, 0.975), .groups = "drop") %>%
  rename(Year = year)

# 3. PREPARE NO-TREND FORECAST
notrend_clean <- bind_rows(no_trend_sims_list) %>%
  left_join(grid_country_map, by = "grid_id") %>%
  mutate(Country = as.character(Country), Scenario = "No-Trend") %>%
  group_by(Country, year, Scenario, sim_id) %>% summarize(Total_N = sum(abundance), .groups="drop") %>%
  group_by(Country, year, Scenario) %>%
  summarize(Median = median(Total_N), Lower = quantile(Total_N, 0.025), Upper = quantile(Total_N, 0.975), .groups = "drop") %>%
  rename(Year = year)

# 4. COMBINE EVERYTHING
# Now all 3 have the same columns: Country, Year, Scenario, Median, Lower, Upper
full_forecast_comparison <- bind_rows(bridge_clean, trend_clean, notrend_clean) %>%
  filter(!is.na(Country), Country != "NA", Country != "")

# 5. DYNAMIC SORTING (FIXED)
# We filter to ONE scenario to ensure each country appears exactly once in the order vector
country_order_vec <- full_forecast_comparison %>% 
  filter(Year == 2044, Scenario == "Climate-Trend") %>% 
  arrange(desc(Median)) %>% 
  pull(Country)

# Now apply this unique vector to the factors
full_forecast_comparison$Country <- factor(full_forecast_comparison$Country, levels = country_order_vec)

# 6. PLOT
ggplot(full_forecast_comparison, aes(x = Year, color = Scenario, fill = Scenario)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, color = NA) +
  geom_line(aes(y = Median), linewidth = 1) +
  geom_vline(xintercept = 2024, linetype = "dashed", color = "grey50") +
  facet_wrap(~Country, scales = "free_y") +
  theme_minimal() +
  labs(title = "Counterfactual Impact of Climate Trends",
       subtitle = "Comparison of Climate-Trend vs. No-Trend Scenarios",
       x = "Year", y = "Total Population (N)") +
  theme(strip.text = element_text(face = "bold"))





cat("\n--- Standardizing & Plotting Unified Comparison ---\n")

# 1. STANDARDIZE HISTORICAL DATA
# We ensure the schema matches the forecast scenarios exactly
hist_ready <- hist_data %>%
  mutate(Scenario = "Historical (Observed)",
         Lower = Median, # Set to Median so ribbon is invisible in history
         Upper = Median) %>%
  select(Country, Year, Scenario, Median, Lower, Upper)

# 2. STANDARDIZE FORECAST SCENARIOS
# Helper function to process the simulation lists
format_sims <- function(sim_list, scenario_name) {
  bind_rows(sim_list) %>%
    left_join(grid_country_map, by = "grid_id") %>%
    mutate(Country = as.character(Country), Scenario = scenario_name) %>%
    filter(!is.na(Country), Country != "NA", Country != "") %>%
    group_by(Country, year, Scenario, sim_id) %>% 
    summarize(Total_N = sum(abundance), .groups="drop") %>%
    group_by(Country, year, Scenario) %>%
    summarize(
      Median = median(Total_N), 
      Lower = quantile(Total_N, 0.025), 
      Upper = quantile(Total_N, 0.975), 
      .groups = "drop"
    ) %>%
    rename(Year = year)
}

trend_ready   <- format_sims(all_sims_list, "Projection (Climate + Trend)")
notrend_ready <- format_sims(no_trend_sims_list, "Projection (Climate Only)")

# 3. COMBINE & STITCH
# We bind the scenarios. Because 2024 exists in all, the lines will be continuous.
full_comparison <- bind_rows(hist_ready, trend_ready, notrend_ready) %>%
  filter(!is.na(Country))

# 4. DYNAMIC SORTING (Using the Climate-Trend median at 2044)
country_order <- full_comparison %>% 
  filter(Year == 2044, Scenario == "Projection (Climate + Trend)") %>% 
  arrange(desc(Median)) %>% 
  pull(Country)

full_comparison$Country <- factor(full_comparison$Country, levels = country_order)

# 5. FINAL PLOT
ggplot(full_comparison, aes(x = Year, color = Scenario, fill = Scenario)) +
  # Use alpha for the ribbons (Historical has alpha=0 effectively)
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, color = NA) +
  geom_line(aes(y = Median), linewidth = 1) +
  geom_vline(xintercept = 2024, linetype = "dashed", color = "grey50") +
  # Clean palette
  scale_color_manual(values = c("Historical (Observed)" = "black", 
                                "Projection (Climate + Trend)" = "seagreen4", 
                                "Projection (Climate Only)" = "firebrick")) +
  scale_fill_manual(values = c("Historical (Observed)" = "transparent", 
                               "Projection (Climate + Trend)" = "seagreen4", 
                               "Projection (Climate Only)" = "firebrick")) +
  facet_wrap(~Country, scales = "free_y", drop = TRUE) +
  theme_minimal() +
  labs(
    title = "Historical Baseline vs. Forecast Scenarios",
    subtitle = "Comparing climate-forced dynamics with vs. without intrinsic temporal trends",
    x = "Year", y = "Total Population (N)"
  ) +
  theme(strip.text = element_text(face = "bold", size = 11),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")


cat("\n--- Standardizing & Plotting Unified Comparison with Stitching & CI ---\n")

# 1. STANDARDIZE HISTORICAL DATA (From existing 'ci_data')
# Ensure column names match forecast: Year (uppercase)
hist_ready <- ci_data %>%
  rename(Year = year) %>%
  mutate(Scenario = "Historical (Observed)") %>%
  select(Country, Year, Scenario, Median, Lower, Upper)

# 2. STANDARDIZE FORECAST SCENARIOS (With 2024 Bridge)
# We add a 2024 data point to every forecast to ensure the stitch
format_sims_with_bridge <- function(sim_list, scenario_name) {
  # Calculate 2024 Bridge state
  bridge_df <- data.frame(
    Country = grid_country_map$Country,
    abundance = N_start_2024
  ) %>%
    group_by(Country) %>%
    summarize(Median = sum(abundance), .groups="drop") %>%
    mutate(Year = 2024, Scenario = scenario_name, Lower = Median, Upper = Median)
  
  # Process forecast sims (2025-2044)
  forecast_df <- bind_rows(sim_list) %>%
    left_join(grid_country_map, by = "grid_id") %>%
    mutate(Country = as.character(Country), Scenario = scenario_name) %>%
    filter(!is.na(Country), Country != "NA", Country != "") %>%
    group_by(Country, year, Scenario, sim_id) %>% 
    summarize(Total_N = sum(abundance), .groups="drop") %>%
    group_by(Country, year, Scenario) %>%
    summarize(
      Median = median(Total_N), 
      Lower = quantile(Total_N, 0.025), 
      Upper = quantile(Total_N, 0.975), 
      .groups = "drop"
    ) %>%
    rename(Year = year)
  
  # Bind bridge (2024) to forecast (2025+)
  bind_rows(bridge_df, forecast_df)
}

trend_ready   <- format_sims_with_bridge(all_sims_list, "Projection (Climate + Trend)")
notrend_ready <- format_sims_with_bridge(no_trend_sims_list, "Projection (Climate Only)")

# 3. COMBINE EVERYTHING
full_comparison <- bind_rows(hist_ready, trend_ready, notrend_ready) %>%
  filter(!is.na(Country), Country != "NA", Country != "")

# 4. DYNAMIC SORTING
# Use the median of the main climate projection in 2044 for consistent sorting
country_order <- full_comparison %>% 
  filter(Year == 2044, Scenario == "Projection (Climate + Trend)") %>% 
  arrange(desc(Median)) %>% 
  pull(Country)

full_comparison$Country <- factor(full_comparison$Country, levels = country_order)

# 5. FINAL PLOT
ggplot(full_comparison, aes(x = Year, color = Scenario, fill = Scenario)) +
  # Ribbons for both historical CI and forecast CI
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, color = NA) +
  geom_line(aes(y = Median), linewidth = 1) +
  # Visual stitch point
  geom_vline(xintercept = 2024, linetype = "dashed", color = "grey50", alpha=0.5) +
  # High-contrast colors
  scale_color_manual(values = c("Historical (Observed)" = "black", 
                                "Projection (Climate + Trend)" = "seagreen4", 
                                "Projection (Climate Only)" = "firebrick")) +
  scale_fill_manual(values = c("Historical (Observed)" = "grey70", 
                               "Projection (Climate + Trend)" = "seagreen4", 
                               "Projection (Climate Only)" = "firebrick")) +
  facet_wrap(~Country, scales = "free_y", drop = TRUE) +
  theme_minimal() +
  labs(
    title = "Historical Baseline vs. Forecast Scenarios",
    subtitle = "Continuous 95% Credible Intervals showing model fit and future uncertainty",
    x = "Year", y = "Total Population (N)"
  ) +
  theme(strip.text = element_text(face = "bold", size = 11),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")
ggsave("Forecast.pdf", width = 11, height = 9)



# LOCK 2024 GAMMA FIELD ####
# instead of zeroing out gamma field

# 1. Capture the trend state at 2024
t_2024 <- max(stan_data$t_scaled)
gamma_locked_vec <- gamma_smooth_est * t_2024

simulate_abundance_DGP_locked <- function(S, T_total, thermal_niche, temp_idx, 
                                          phi, gamma_locked, eps_year, 
                                          r, N0_proportion, w_norm, from_idx, to_idx, E_disp,
                                          N_init, max_K_cap = 5000) {
  
  R0 <- exp(r); R0m1 <- R0 - 1.0
  N_curr <- N_init 
  trajectory <- matrix(NA, nrow = S, ncol = T_total)
  
  for (t in 1:T_total) {
    # Extract indices for current year
    idx_range <- ((t - 1) * S + 1):(t * S)
    
    # Calculate Capacity K using the LOCKED gamma_locked_vec
    # Note: We removed (gamma_smooth * t_scaled[t]) and replaced it with gamma_locked
    logK <- thermal_niche[temp_idx[idx_range]] + phi + gamma_locked
    K <- pmin(exp(logK), max_K_cap)
    
    # Beverton-Holt Growth
    N_curr <- (R0 * N_curr) / (1.0 + (R0m1 / (K + 1e-9)) * N_curr)
    
    # Dispersal
    N_new <- numeric(S)
    for (e in 1:E_disp) {
      N_new[to_idx[e]] <- N_new[to_idx[e]] + w_norm[e] * N_curr[from_idx[e]]
    }
    
    # Process Error
    N_curr <- N_new * exp(eps_year[t])
    trajectory[, t] <- N_curr
  }
  return(trajectory)
}

# Run the stochastic forecast
all_sims_list_locked <- lapply(1:n_sims, function(i) {
  # 1. SAMPLE PARAMETERS
  draw_idx <- sample(1:nrow(draws_df), 1)
  draw <- draws_df[draw_idx, ]
  
  # 2. EXTRACT PARAMETERS
  curr_phi     <- as.numeric(draw[grep("^phi\\[", names(draw))])
  curr_gamma   <- as.numeric(draw[grep("^gamma_smooth\\[", names(draw))])
  curr_r       <- as.numeric(draw$r)
  curr_N0      <- as.numeric(draw$N0_proportion)
  curr_alpha   <- as.numeric(exp(draw$log_alpha))
  curr_sigma   <- as.numeric(draw$sigma_year)
  
  # 3. RECONSTRUCT NICHE (Must be inside to use draw parameters!)
  curr_niche <- as.numeric(draw$thermal_max) - 
    log(1 + exp(-as.numeric(draw$slope_L) * (stan_data$thermal_gradient - as.numeric(draw[["T_nodes[1]"]])))) - 
    log(1 + exp(-as.numeric(draw$slope_R) * (as.numeric(draw[["T_nodes[2]"]]) - stan_data$thermal_gradient)))
  
  # 4. RECOMPUTE DISPERSAL WEIGHTS (Must be inside because alpha changes!)
  curr_w <- exp(-stan_data$dists / curr_alpha)
  curr_out_sum <- rep(1e-12, S)
  for (e in 1:stan_data$E_disp) {
    curr_out_sum[stan_data$from_idx[e]] <- curr_out_sum[stan_data$from_idx[e]] + curr_w[e]
  }
  curr_w_norm <- curr_w / curr_out_sum[stan_data$from_idx]
  
  # 5. CALCULATE LOCKED TREND
  t_2024 <- max(stan_data$t_scaled)
  curr_gamma_locked <- curr_gamma * t_2024
  
  # 6. RUN SIMULATION
  sim_res <- simulate_abundance_DGP_locked(
    S = S, T_total = 20, 
    thermal_niche = curr_niche,
    temp_idx = as.vector(future_temp_idx),
    phi = curr_phi, 
    gamma_locked = curr_gamma_locked, 
    eps_year = rnorm(20, 0, curr_sigma),
    r = curr_r,
    N0_proportion = curr_N0, 
    w_norm = curr_w_norm,
    from_idx = stan_data$from_idx, 
    to_idx = stan_data$to_idx, 
    E_disp = stan_data$E_disp,
    N_init = N_start_2024
  )
  
  data.frame(year = rep(2025:2044, each = S), 
             grid_id = rep(my_grid_final$grid_id, 20), 
             abundance = as.vector(sim_res), 
             sim_id = i)
})


cat("\n--- Standardizing & Plotting Unified Comparison with Stitching & CI ---\n")

# 1. STANDARDIZE HISTORICAL DATA (From existing 'ci_data')
# Ensure column names match forecast: Year (uppercase)
hist_ready <- ci_data %>%
  rename(Year = year) %>%
  mutate(Scenario = "Historical (Observed)") %>%
  select(Country, Year, Scenario, Median, Lower, Upper)

# 2. STANDARDIZE FORECAST SCENARIOS (With 2024 Bridge)
# We add a 2024 data point to every forecast to ensure the stitch
format_sims_with_bridge <- function(sim_list, scenario_name) {
  # Calculate 2024 Bridge state
  bridge_df <- data.frame(
    Country = grid_country_map$Country,
    abundance = N_start_2024
  ) %>%
    group_by(Country) %>%
    summarize(Median = sum(abundance), .groups="drop") %>%
    mutate(Year = 2024, Scenario = scenario_name, Lower = Median, Upper = Median)
  
  # Process forecast sims (2025-2044)
  forecast_df <- bind_rows(sim_list) %>%
    left_join(grid_country_map, by = "grid_id") %>%
    mutate(Country = as.character(Country), Scenario = scenario_name) %>%
    filter(!is.na(Country), Country != "NA", Country != "") %>%
    group_by(Country, year, Scenario, sim_id) %>% 
    summarize(Total_N = sum(abundance), .groups="drop") %>%
    group_by(Country, year, Scenario) %>%
    summarize(
      Median = median(Total_N), 
      Lower = quantile(Total_N, 0.025), 
      Upper = quantile(Total_N, 0.975), 
      .groups = "drop"
    ) %>%
    rename(Year = year)
  
  # Bind bridge (2024) to forecast (2025+)
  bind_rows(bridge_df, forecast_df)
}

trend_ready   <- format_sims_with_bridge(all_sims_list, "Projection (Climate + Trend)")
notrend_ready <- format_sims_with_bridge(all_sims_list_locked, "Projection (Climate Only)")

# 3. COMBINE EVERYTHING
full_comparison <- bind_rows(hist_ready, trend_ready, notrend_ready) %>%
  filter(!is.na(Country), Country != "NA", Country != "")

# 4. DYNAMIC SORTING
# Use the median of the main climate projection in 2044 for consistent sorting
country_order <- full_comparison %>% 
  filter(Year == 2044, Scenario == "Projection (Climate + Trend)") %>% 
  arrange(desc(Median)) %>% 
  pull(Country)

full_comparison$Country <- factor(full_comparison$Country, levels = country_order)

# 5. FINAL PLOT
ggplot(full_comparison, aes(x = Year, color = Scenario, fill = Scenario)) +
  # Ribbons for both historical CI and forecast CI
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, color = NA) +
  geom_line(aes(y = Median), linewidth = 1) +
  # Visual stitch point
  geom_vline(xintercept = 2024, linetype = "dashed", color = "grey50", alpha=0.5) +
  # High-contrast colors
  scale_color_manual(values = c("Historical (Observed)" = "black", 
                                "Projection (Climate + Trend)" = "seagreen4", 
                                "Projection (Climate Only)" = "firebrick")) +
  scale_fill_manual(values = c("Historical (Observed)" = "grey70", 
                               "Projection (Climate + Trend)" = "seagreen4", 
                               "Projection (Climate Only)" = "firebrick")) +
  facet_wrap(~Country, scales = "free_y", drop = TRUE) +
  theme_minimal() +
  labs(
    title = "Historical Baseline vs. Forecast Scenarios",
    subtitle = "Continuous 95% Credible Intervals showing model fit and future uncertainty",
    x = "Year", y = "Total Population (N)"
  ) +
  theme(strip.text = element_text(face = "bold", size = 11),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")
ggsave("Forecast.pdf", width = 11, height = 9)



# ============================================================================== #
# 19. FORECAST PROPORTIONAL RANGE EXPANSION MAP (2024-2044) ####
# ============================================================================== #
cat("\n--- 19. Generating Forecast Proportional Range Expansion Map ---\n")

# 1. Extract 2044 Median Abundance per grid_id from simulation results
# We use the median of all 50 simulations (or 200, depending on your latest run)
df_2044 <- bind_rows(all_sims_list) %>%
  filter(year == 2044) %>%
  group_by(grid_id) %>%
  summarize(N_2044 = median(abundance), .groups = "drop")

# 2. Align with 2024 Bridge
# Create a dataframe for 2024 based on the existing N_start_2024 vector
# We need to ensure grid_id order matches my_grid_final
df_2024 <- data.frame(
  grid_id = as.character(my_grid_final$grid_id),
  N_2024 = N_start_2024
)

# 3. Join and calculate Prop_Change
forecast_change_df <- df_2044 %>%
  left_join(df_2024, by = "grid_id") %>%
  mutate(
    # Avoid division by zero with + 1, similar to your historical calculation
    Prop_Change = (N_2044 - N_2024) / (N_2024 + 1)
  )

# 4. Join back to spatial grid
grid_forecast_expansion <- my_grid_final %>%
  filter(grid_id %in% valid_sites) %>%
  left_join(forecast_change_df, by = "grid_id")

# 5. Plot with fixed limits (matching your preferred historical scale)
p_forecast_expansion <- ggplot(grid_forecast_expansion) +
  geom_sf(aes(fill = Prop_Change), color = NA) +
  scale_fill_gradient2(
    low = "#d7191c",    # Red (Contraction)
    mid = "grey95",     # Neutral
    high = "#2c7bb6",   # Blue (Expansion)
    midpoint = 0, 
    limits = c(-0.7, 0.7),       # Using your established scale
    oob = scales::squish,        # Force values into the scale
    labels = scales::percent, 
    name = "Rel. Change\n(2024-2044)"
  ) +
  theme_minimal() +
  labs(
    title = "Forecasted Proportional Population Shift",
    subtitle = "Relative change in abundance (2024 vs 2044)",
    caption = "Values capped at ±70% for visualization clarity"
  )

print(p_forecast_expansion)
ggsave("map_forecast_prop_expansion.pdf", width = 8, height = 8)


# ============================================================================== #
# 20. COMPARING PROPORTIONAL SHIFTS: TREND VS NO-TREND ####
# ============================================================================== #
cat("\n--- 20. Generating Comparative Proportional Maps ---\n")

# 1. Helper function to extract 2044 median and calculate Prop Change
get_prop_change <- function(sim_list, name) {
  df_2044 <- bind_rows(sim_list) %>%
    filter(year == 2044) %>%
    group_by(grid_id) %>%
    summarize(N_2044 = median(abundance), .groups = "drop")
  
  # Join with the N_2024 bridge (N_start_2024 is the 2024 state)
  df_2024 <- data.frame(grid_id = as.character(my_grid_final$grid_id), N_2024 = N_start_2024)
  
  df_2044 %>%
    left_join(df_2024, by = "grid_id") %>%
    mutate(
      Prop_Change = (N_2024 - N_2024) / (N_2024 + 1), # Placeholder
      Prop_Change = (N_2044 - N_2024) / (N_2024 + 1),
      Scenario = name
    )
}

# 2. Process both scenarios
change_trend   <- get_prop_change(all_sims_list, "With Spatiotemporal Trend")
change_notrend <- get_prop_change(all_sims_list_locked, "Climate Only (No Trend)")

# 3. Combine
full_comparison_map <- bind_rows(change_trend, change_notrend) %>%
  left_join(my_grid_final %>% select(grid_id), by = "grid_id") %>%
  st_as_sf()

# 4. Plot
p_comparison_maps <- ggplot(full_comparison_map) +
  geom_sf(aes(fill = Prop_Change), color = NA) +
  facet_wrap(~Scenario, ncol = 2) +
  scale_fill_gradient2(
    low = "#d7191c", mid = "grey95", high = "#2c7bb6",
    midpoint = 0, limits = c(-0.7, 0.7), oob = scales::squish,
    labels = scales::percent, name = "Rel. Change\n(2024-2044)"
  ) +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold", size = 12)) +
  labs(
    title = "Impact of Spatiotemporal Trends on Projected Abundance",
    subtitle = "Comparing scenarios with and without unmodeled temporal trends"
  )

print(p_comparison_maps)
ggsave("ForecastMap.pdf", width = 10, height = 5)




