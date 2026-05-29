setwd("~/School/2025 - 2026/Thesis Statistics/Analyses/FINAL SCRIPTS/Simulation_260525")

# ============================================================================== #
# LOCAL POST-PROCESSING & DIAGNOSTICS GENERATION (STRESS-TEST ALIGNED)
# ============================================================================== #
library(cmdstanr)
library(posterior)
library(bayesplot)
library(ggplot2)
library(dplyr)
library(tidyr)

theme_set(theme_minimal(base_size = 12))
color_scheme_set("viridis")

fit_local <- readRDS("fit_generative_validation.rds")

# ============================================================================== #
# DIAGNOSTIC PRINTS & TEXT OUTPUTS ####
# ============================================================================== #
sink("simulation_diagnostic_report.txt") 
cat("======================================================================\n")
cat("          MCMC CONVERGENCE & EFFICIENCY REPORT\n")
cat("======================================================================\n\n")

cat("--- Core Parameter Posterior Summaries ---\n")
param_summary <- fit_local$summary(c("log_alpha", "r", 
                                     "thermal_max", "T_nodes", 
                                     "slope_L", "slope_R",
                                     "logit_p", "N0_proportion", 
                                     "rho_phi", "rho_gamma", "sigma_gamma",
                                     "sigma_phi", "sigma_year"))
print(param_summary)
cat("\n")

cat("--- Latent High-Dimensional Layer Summary Statistics ---\n")
latent_summary <- fit_local$summary(c("phi_eta", "eps_year_raw"))

cat("Spatial & Temporal Latent Parameters Quantiles:\n")
cat("  R-hat range:    ", min(latent_summary$rhat, na.rm=TRUE), " to ", max(latent_summary$rhat, na.rm=TRUE), "\n")
cat("  Bulk ESS range: ", min(latent_summary$ess_bulk, na.rm=TRUE), " to ", max(latent_summary$ess_bulk, na.rm=TRUE), "\n")
cat("  Tail ESS range: ", min(latent_summary$ess_tail, na.rm=TRUE), " to ", max(latent_summary$ess_tail, na.rm=TRUE), "\n\n")

cat("--- HMC/NUTS Structural Diagnostics ---\n")
diagnostics <- fit_local$diagnostic_summary()
cat("Total Divergent Transitions: ", sum(diagnostics$num_divergent), "\n")
cat("Max Treedepth Exceedances:   ", sum(diagnostics$num_max_treedepth), "\n")
cat("Low E-BFMI Chains:           ", length(which(diagnostics$ebfmi < 0.2)), "\n\n")

cat("--- Execution Times per Chain (Seconds) ---\n")
print(fit_local$time()$chains)
sink() 

cat("Text report compiled as 'simulation_diagnostic_report.txt'\n")

param_summary
latent_summary
diagnostics


# 1. Update the Summary call to include gamma_eta
cat("--- Latent High-Dimensional Layer Summary Statistics ---\n")
latent_summary <- fit_local$summary(c("phi_eta", "gamma_eta", "eps_year_raw"))

# 2. Extract random subsets for a more granular look in the report
set.seed(123)
sampled_phi   <- sample(grep("phi_eta", latent_summary$variable, value = TRUE), 4)
sampled_gamma <- sample(grep("gamma_eta", latent_summary$variable, value = TRUE), 4)
sampled_eps   <- sample(grep("eps_year_raw", latent_summary$variable, value = TRUE), 4)

# 3. Print the subset details to the report
sink("simulation_diagnostic_report.txt", append = TRUE)
cat("\n--- Random Subsample of Hierarchical Parameters ---\n")
print(latent_summary[latent_summary$variable %in% c(sampled_phi, sampled_gamma, sampled_eps), ])
sink()

# 4. Print the general Latent Summary to console for quick review
print(latent_summary[latent_summary$variable %in% c(sampled_phi, sampled_gamma, sampled_eps), ])

# ============================================================================== #
# GRAPHICAL PLOTS FOR REPORT ####
# ============================================================================== #

# --- PLOT 1: PARALLEL MIXING TRACEPLOTS (CORE) ---
core_params <- c("r", "log_alpha", "logit_p", "N0_proportion", "thermal_max", 
                 "T_nodes", "slope_L", "slope_R", 
                 "rho_phi", "rho_gamma", "rho_gamma", "sigma_phi", "sigma_gamma", "sigma_year", "lp__")

p_trace <- mcmc_trace(fit_local$draws(variables = core_params)) +
  theme_minimal(base_size = 15)
p_trace
ggsave("plot_diagnostics_trace.pdf", plot = p_trace, width = 8, height = 4, dpi = 300)


# --- PLOT 1b: SUBSET TRACEPLOTS FOR LATENT SPATIAL & TEMPORAL PARAMETERS ---
set.seed(123) 
sampled_knots <- paste0("phi_eta[", sample(1:20, 4), "]")
sampled_trends <- paste0("gamma_eta[", sample(1:20, 4), "]")
sampled_years <- paste0("eps_year_raw[", sample(1:30, 4), "]")
latent_monitor_vars <- c(sampled_knots, sampled_trends, sampled_years)

p_latent_trace <- mcmc_trace(fit_local$draws(variables = latent_monitor_vars)) +
  theme_minimal(base_size = 15)
p_latent_trace
ggsave("plot_diagnostics_latent_trace.pdf", plot = p_latent_trace, width = 8, height = 4, dpi = 300)



library(bayesplot)
library(ggplot2)
library(posterior)

# 1. Combine all variables into one vector
all_vars <- c("r", "log_alpha", "logit_p", "N0_proportion", "thermal_max", 
              "T_nodes", "slope_L", "slope_R", "rho_phi", "rho_gamma", "rho_gamma", 
              "sigma_phi", "sigma_gamma", "sigma_year", "lp__",
              sampled_knots, sampled_trends, sampled_years)

# 2. Extract the unified draws
unified_draws <- fit_local$draws(variables = all_vars)

# 3. Create a single, massive trace plot
p_unified_trace <- mcmc_trace(unified_draws) +
  theme_minimal(base_size = 10) +
  theme(strip.text = element_text(size = 8)) # Shrink labels to fit 22 panels
p_unified_trace
# 4. Save with a very large height to keep 22 rows legible
ggsave("plot_diagnostics_unified_trace.pdf", plot = p_unified_trace, 
       width = 10, height = 7, dpi = 300)



# --- PLOT 2: THERMAL NICHE ENVELOPE (FIXED) ---

# 1. Extract variables correctly (Stan returns T_nodes as [1] and [2])
draws_df <- as_draws_df(fit_local$draws(variables = c("thermal_max", "T_nodes", "slope_L", "slope_R")))

plot_temps <- seq(-10, 25, length.out = 250)

calc_posterior_K <- function(temp) {
  # Use the Stan formula directly
  # logK = thermal_max - log1p_exp(-slope_L * (x - T1)) - log1p_exp(-slope_R * (T2 - x))
  
  # Accessing the nodes from the fit object (using backticks for index notation)
  T1 <- draws_df$`T_nodes[1]`
  T2 <- draws_df$`T_nodes[2]`
  
  log_K_raw <- draws_df$thermal_max - 
    log(1 + exp(-draws_df$slope_L * (temp - T1))) - 
    log(1 + exp(-draws_df$slope_R * (T2 - temp)))
  
  return(pmin(pmax(exp(log_K_raw), 1e-6), 1e5)) 
}

# 2. Update the "True Niche" calculation to match Stan's math exactly
# Using the values from your "Stress Test" script:
true_T_nodes <- c(2, 11) 
true_slope_L <- 1.2
true_slope_R <- 1.3
true_thermal_max <- 3.0

calc_true_niche <- function(temp) {
  log_K_raw <- true_thermal_max - 
    log(1 + exp(-true_slope_L * (temp - true_T_nodes[1]))) - 
    log(1 + exp(-true_slope_R * (true_T_nodes[2] - temp)))
  return(pmin(pmax(exp(log_K_raw), 1e-6), 1e4))
}

# 3. Rest of the loop stays the same
plot_data_95 <- data.frame()
for (temp in plot_temps) {
  K_vals <- calc_posterior_K(temp)
  plot_data_95 <- rbind(plot_data_95, data.frame(
    Temperature = temp,
    K_median    = median(K_vals),
    K_low_95    = quantile(K_vals, 0.025), # 95% Lower bound
    K_high_95   = quantile(K_vals, 0.975)  # 95% Upper bound
  ))
}
plot_data_95$K_true <- calc_true_niche(plot_data_95$Temperature)

p_niche <- ggplot(plot_data_95, aes(x = Temperature)) +
  # Single ribbon for 95%
  geom_ribbon(aes(ymin = K_low_95, ymax = K_high_95, fill = "95% credible interval"), alpha = 0.25) +
  geom_line(aes(y = K_median, color = "Posterior median"), linewidth = 1.2) +
  geom_line(aes(y = K_true, color = "True simulated niche"), linewidth = 1.2, linetype = "dashed") +
  scale_color_manual(name = "Curves", values = c("Posterior median" = "#440154FF", "True simulated niche" = "#D55E00")) +
  scale_fill_manual(name = "Uncertainty", values = c("95% credible interval" = "#440154FF")) +
  labs(title = "Posterior Thermal Niche (95% Credible Interval)",
       x = "Temperature (°C)", y = "Carrying Capacity (K)") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

print(p_niche)
ggsave("plot_niche.pdf", plot = p_niche, width = 8, height = 5.5, dpi = 300)


# --- PLOT 3: PRIOR-POSTERIOR UPDATING GRID WITH CORRECTED EXPANDED EVALUATION FIELDS ---

# 1. Define true values in a tidy format (One row per parameter)
# Make sure these names EXACTLY match your Stan output variables
true_values <- data.frame(
  variable = c("r", "log_alpha", "logit_p", "N0_proportion", "thermal_max", 
               "T_nodes[1]", "T_nodes[2]", "slope_L", "slope_R", 
               "sigma_phi", "sigma_gamma", "sigma_year", "rho_phi", "rho_gamma"),
  true_val = c(0.1, 1.0, -1.5, 0.3, 3.0, 2.0, 11.0, 1.2, 1.3, 
               0.12, 0.08, 0.03, 80.0, 180.0)
)

# 2. Extract draws (This will now grab the individual elements correctly)
draws_long <- as_draws_df(fit_local$draws(variables = true_values$variable)) %>%
  as.data.frame() %>%
  dplyr::select(dplyr::all_of(true_values$variable)) %>%
  pivot_longer(cols = everything(), names_to = "variable", values_to = "value")

draws_paneled <- merge(draws_long, true_values, by = "variable")

# 3. DYNAMICALLY GENERATE PRIOR CURVES
prior_curves <- draws_paneled %>%
  group_by(variable) %>%
  do({
    var_name <- unique(.$variable)
    
    # Define range for density evaluation
    eval_min <- min(.$value) - diff(range(.$value))*0.5
    eval_max <- max(.$value) + diff(range(.$value))*0.5
    grid <- seq(eval_min, eval_max, length.out = 400)
    
    # Define densities based on your specific priors
    dens <- case_when(
      var_name == "r"              ~ dnorm(grid, 0.3, 0.3),
      var_name == "log_alpha"      ~ dnorm(grid, 2.5, 1.0),
      var_name == "logit_p"        ~ dnorm(grid, -2.0, 1.0),
      var_name == "N0_proportion"  ~ dbeta(grid, 2, 2),
      var_name == "thermal_max"    ~ dnorm(grid, 0, 3.0),
      var_name == "T_nodes[1]"     ~ dnorm(grid, 10, 8.0),
      var_name == "T_nodes[2]"     ~ dnorm(grid, 10, 8.0),
      var_name == "slope_L"        ~ dlnorm(grid, 0, 0.5),
      var_name == "slope_R"        ~ dlnorm(grid, 0, 0.5),
      var_name == "sigma_phi"      ~ dt(grid/0.1, df=3)/0.1,
      var_name == "sigma_gamma"    ~ dnorm(grid, 0, 0.1),
      var_name == "sigma_year"     ~ dt(grid/0.2, df=3)/0.2,
      var_name == "rho_phi"        ~ dlnorm(grid, log(100), 0.5),
      var_name == "rho_gamma"      ~ dlnorm(grid, log(150), 0.5),
      TRUE                         ~ as.numeric(NA)
    )
    data.frame(value = grid, density = dens)
  }) %>% ungroup()

# 4. Plot
p_densities_expanded <- ggplot() +
  geom_line(data = prior_curves, aes(x = value, y = density, linetype = "Prior"), 
            color = "gray35", linewidth = 0.75) +
  geom_density(data = draws_paneled, aes(x = value, fill = "Posterior"), 
               alpha = 0.4, color = "#21918cFF", linewidth = 0.8) +
  geom_vline(data = true_values, aes(xintercept = true_val, color = "True Value"), 
             linewidth = 1.0) +
  facet_wrap(~variable, scales = "free", ncol = 4) + 
  scale_fill_manual(name = "Distributions", values = c("Posterior" = "#21918cFF")) +
  scale_color_manual(name = "Reference", values = c("True Value" = "#D55E00")) +
  scale_linetype_manual(name = "Priors", values = c("Prior" = "dashed")) +
  labs(title = "Prior-Posterior Updating Grid",
       x = "Parameter Space", y = "Density") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold"), legend.position = "bottom")

print(p_densities_expanded)
ggsave("plot_parameter_densities_grid.pdf", plot = p_densities_expanded, width = 11, height = 7.5, dpi = 300)


# --- PLOT 4: AUTOCORRELATION ---
p_acf <- mcmc_acf(fit_local$draws(variables = c("T_nodes", "log_alpha")))
p_acf
ggsave("appendix_autocorrelation.pdf", plot = p_acf, width = 8, height = 5)


# --- PLOT 5: POSTERIOR CREDIBLE INTERVAL WHISKERS ---
p_intervals <- mcmc_intervals(fit_local$draws(variables = c("T_nodes", "thermal_max", "log_alpha"))) +
  labs(title = "Posterior Credible Intervals")
p_intervals
ggsave("plot_parameter_intervals.pdf", plot = p_intervals, width = 7, height = 4)


# --- Correlations ---
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
combined_draws <- fit_local$draws()
final_draws <- subset_draws(combined_draws, chain = c(1,2,3,4,5,6))
# Plot the trace for the key biological parameters and the log-posterior
color_scheme_set("viridis")
arr <- as_draws_array(final_draws)
subset_arr <- arr[, , core_params]
corr <- mcmc_pairs(subset_arr)
corr
ggsave("Corr.pdf", plot = corr, 
       width = 15, height = 20)
