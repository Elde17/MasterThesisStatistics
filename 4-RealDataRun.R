# NOTE: this script is designed to run on a virtual computer

library(Microsoft365R)
library(cmdstanr)
library(posterior)

# "device_code" is required for headless servers
od <- get_personal_onedrive(auth_type = "device_code")
cat("Background Job is running in:", getwd(), "\n")

# ============================================================================== #
# 1. DOWNLOAD & PREPARE ####
# ============================================================================== #
cat("\n--- 1. Downloading Files from OneDrive ---\n")
# Download the specific data and the NEW optimized Stan script
od$download_file("Param_CSR_GP_ST_optimized.stan", dest="Param_CSR_GP_ST_optimized.stan", overwrite=TRUE)
od$download_file("stan_data_ZwaHei_HybridIDE.rds", dest="stan_data_ZwaHei_HybridIDE.rds", overwrite=TRUE)
stan_data <- readRDS("stan_data_ZwaHei_HybridIDE.rds")

# ==============================================================================
# 2. COMPILATION WITH ABSOLUTE CORE ENFORCEMENT
# ==============================================================================
cat("\n--- 2. Compiling Optimized Model Architecture ---\n")

# Turning stan_threads to FALSE guarantees Stan will 
# never spawn background threads. Each chain is strictly locked to 1 CPU core.
cpp_options <- list(
  stan_threads = FALSE, 
  CXXFLAGS = "-O3 -march=native -mtune=native"
)
cmdstan_make_local(cpp_options = cpp_options)

mod <- cmdstan_model("Param_CSR_GP_ST_optimized.stan", force_recompile = TRUE)

# ============================================================================== #
# 3. DATA PATCH & PREPARATION ####
# ============================================================================== #
cat("\n--- 3. Preparing Data with Safety Checks ---\n")

# 1. Add Needed Simulation/State Variables
stan_data$idx_seen <- which(stan_data$Y_vec > 0)
stan_data$idx_zero <- which(stan_data$Y_vec == 0)
stan_data$N_seen <- length(stan_data$idx_seen)
stan_data$N_zero <- length(stan_data$idx_zero)
stan_data$N_state_total <- length(stan_data$Y_vec)
stan_data$map_state_idx <- 1:length(stan_data$Y_vec)

# 2. Temperature Formatting
if(is.matrix(stan_data$temp_idx)) {
  stan_data$temp_idx <- as.vector(t(stan_data$temp_idx))
}
if(is.matrix(stan_data$thermal_gradient)) {
  stan_data$thermal_gradient <- as.vector(stan_data$thermal_gradient)
}

# 3. Final Dimensions Check
cat("Grid Cells (S):", stan_data$S, "\n")
cat("GP Knots (N_spatial_bf):", stan_data$N_spatial_bf, "\n")

cat("Data Patch Complete. Ready for Optimized CSR Run.\n")

# ============================================================================== #
# COMPILATION & UPDATED INITS ####
# ============================================================================== #
make_init_hybrid <- function() {
  list(
    thermal_max = rnorm(1, 0.0, 1.0),  # Weakly informative prior baseline matching chapter text
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
    eps_year_raw = rnorm(30, mean = 0, sd = 0.01)
  )
}

# ============================================================================== #
# MCMC RUN ####
# ============================================================================== #
cat("\n--- 5. Running Production MCMC ---\n")

run_name <- "ZwaHei_CSR" 

fit_real <- mod$sample(
  data = stan_data,
  chains = 6, 
  parallel_chains = 6,   
  iter_warmup = 200,     
  iter_sampling = 200,
  init = make_init_hybrid, 
  adapt_delta = 0.9,    
  max_treedepth = 12,    
  output_dir = ".",
  output_basename = run_name,
  refresh = 50,
  save_cmdstan_config = TRUE
)

# ============================================================================== #
# SAVE & UPLOAD ####
# ============================================================================== #
cat("\n--- 6. Saving and Uploading Results ---\n")

# Local saves
saveRDS(fit_real, paste0("fit_", run_name, ".rds"))
draws <- fit_real$draws()
saveRDS(draws, paste0("draws_", run_name, ".rds"))

# Upload to Thesis_Data folder
od$upload_file(paste0("fit_", run_name, ".rds"), 
               dest=paste0("Thesis_Data/fit_", run_name, ".rds"))
od$upload_file(paste0("draws_", run_name, ".rds"), 
               dest=paste0("Thesis_Data/draws_", run_name, ".rds"))

# Zip and upload raw CSV chains
csv_files <- list.files(pattern = paste0("^", run_name, ".*\\.csv$"))
zip_name <- paste0("raw_chains_", run_name, ".zip")
zip(zip_name, files = csv_files)
od$upload_file(zip_name, dest=paste0("Thesis_Data/", zip_name))

cat("\n--- Final Model Run Complete! ---\n")