# ============================================================================== #
# 0. SETUP & LIBRARIES ####
# ============================================================================== #
rm(list = ls())
gc()

library(readr)
library(sf)
library(data.table)
library(splines)
library(cmdstanr)

# Ensure cmdstan is ready
check_cmdstan_toolchain(fix = TRUE)

setwd("~/School/2025 - 2026/Thesis Statistics/Analyses/Param_CSR_260517/ZwaHei")

# File paths (Adjust these if your working directory changes)
data_path <- "C:/Users/lieve/Documents/School/2025 - 2026/Thesis Statistics/Data/New/data_processed.rds"
grid_path <- "C:/Users/lieve/Documents/School/2025 - 2026/Thesis Statistics/Data/New/grid.gpkg"
era5_monthly_path <- "C:/Users/lieve/Documents/School/2025 - 2026/Thesis Statistics/Data/era5_monthly_temps.csv"

# ==============================================================================
# 1. LOAD DATA & PHENOLOGY AUDIT (Using Precision DOY) ####
# ==============================================================================
cat("\n--- 1. Loading Data & Precision DOY Audit ---\n")
data <- read_rds(data_path)
setDT(data)

# Define your target flight window in DOY
start_doy <- 0 
end_doy   <- 366
target_species <- "Sympetrum danae"

# 1. Explicitly define all non-species columns from your new dataset
meta_cols <- c("visit", "year", "year_scaled", "doy_scaled", "doy", 
               "ll", "ll_raw", "spread_scaled", "grid_id",
               "parentDatasetID_2", "datasetID")

# 2. Dynamically grab all OTHER columns as species
species_cols <- setdiff(names(data), meta_cols)

# Check specifically for observations outside the DOY window
lost_tk <- data[(doy < start_doy | doy > end_doy) & get(target_species) > 0]
lost_tk_count <- sum(lost_tk[[target_species]], na.rm = TRUE)

if(lost_tk_count > 0) {
  cat(paste0("WARNING: ", lost_tk_count, " observations of ", target_species, 
             " found outside the DOY window (", start_doy, "-", end_doy, ")!\n"))
  cat("Earliest detection DOY:", min(lost_tk$doy), "\n")
  cat("Latest detection DOY:", max(lost_tk$doy), "\n")
} else {
  cat("Success: No observations found outside the DOY window.\n")
}

# Apply the filter to the raw data (Cleaning Effort and Detections)
data <- data[doy >= start_doy & doy <= end_doy]

# --- CORRECTED AGGREGATION ---
# This ensures that both Visits (K) and Detections (Y) are counted by UNIQUE DAYS.
# This guarantees that Y <= K, which keeps the Binomial distribution happy.
agg_dt <- data[, .(
  # 1. Total unique days the site was surveyed
  num_visits = uniqueN(doy), 
  
  # 2. Total unique days the species was detected at least once
  `Sympetrum danae` = uniqueN(doy[get(target_species) > 0])
  
), by = .(grid_id, year)]

# --- CRITICAL SAFETY CHECK ---
# This must return 0. If it's > 0, the model will crash again.
math_violation <- sum(agg_dt[[target_species]] > agg_dt$num_visits)
cat("Number of records violating Y <= K logic:", math_violation, "\n")

if(math_violation > 0) stop("Data alignment error: Detections exceed visits!")

cat("Aggregation complete. Unique grid-years available:", nrow(agg_dt), "\n")

# ============================================================================== #
# 2. PHENOLOGY FILTERING (95% Quantile via Spline) ####
# ============================================================================== #
cat("\n--- 2. Calculating 95% Phenology Window ---\n")

# Subset data where the species was actually detected to model the flight curve
pos_detections <- data[get(target_species) > 0, .(doy)]

if(nrow(pos_detections) < 50) {
  stop("Too few detections to fit a reliable phenology spline.")
}

# Fit a simple density spline to the Day of Year
# We use a density estimate to find where the bulk of the population is active
dens <- density(pos_detections$doy, bw = "nrd0", adjust = 1.5)

# Calculate the cumulative distribution to find the 95% interval
prob_dist <- data.table(x = dens$x, y = dens$y)
prob_dist[, cumulative := cumsum(y) / sum(y)]

# Extract the 2.5% and 97.5% quantiles (the 95% window)
start_95 <- round(prob_dist[cumulative >= 0.025, min(x)])
end_95   <- round(prob_dist[cumulative >= 0.975, min(x)])

cat(paste0("Calculated 95% flight window for ", target_species, ": DOY ", start_95, " to ", end_95, "\n"))

# Visual Check of the Spline/Density (Highly recommended for your thesis)
plot(dens, main = paste("Phenology Spline:", target_species), xlab = "Day of Year")
abline(v = c(start_95, end_95), col = "red", lty = 2)
mtext("95% Window", at = start_95, col = "red")

# Apply the 95% Filter
# This removes visits that happened outside the core flight season
n_before <- nrow(data)
data <- data[doy >= start_95 & doy <= end_95]
n_after <- nrow(data)

cat(paste("Filtered out", n_before - n_after, "visits outside the 95% window.\n"))

# ==============================================================================
# 2. DYNAMIC ERA5 CLIMATE INTEGRATION & MAINLAND CONNECTIVITY ####
# ==============================================================================
cat("\n--- 2. Integrating Full-Year ERA5 Temperatures & Mainland Connectivity ---\n")

library(sf)
library(data.table)
library(igraph)

# 1. Load the grid
my_grid_all <- st_read(grid_path, quiet = TRUE)
my_grid_all$grid_id <- as.character(my_grid_all$grid_id)

# 2. DEFINE THE MAINLAND (The "Geography-First" approach)
# Identify where surveys happened to find our "anchor" points
visited_ids <- unique(agg_dt$grid_id)
grid_visited <- my_grid_all[my_grid_all$grid_id %in% visited_ids, ]

# Buffer to capture gaps/coastline (200km ensures we bridge connectivity gaps)
cat("Identifying connectivity buffers and mainland extent...\n")
buffer_zone <- st_union(st_geometry(grid_visited)) %>% st_buffer(dist = 200000) 
my_grid_europe <- my_grid_all[st_intersects(my_grid_all, buffer_zone, sparse = FALSE), ]

# Connectivity Graph to remove islands (Cyprus/Ocean fragments)
centroids_eu <- st_centroid(my_grid_europe)
nb_geog <- st_is_within_distance(centroids_eu, dist = 115000) # Bridge 50km cells
edges_geog <- data.table(
  from = rep(1:length(nb_geog), times = sapply(nb_geog, length)),
  to = unlist(nb_geog)
)
g_geog <- graph_from_data_frame(edges_geog, directed = FALSE)
mainland_idx <- which.max(components(g_geog)$csize)
all_mainland_ids <- my_grid_europe$grid_id[which(components(g_geog)$membership == mainland_idx)]

# Final Grid Object for Model
my_grid_mainland <- my_grid_all[my_grid_all$grid_id %in% all_mainland_ids, ]

# 3. Load and Filter ERA5
era5_raw <- fread(era5_monthly_path)
era5_yearly <- era5_raw[grid_id %in% all_mainland_ids & year <= 2024, 
                        .(mean_temp = mean(temp_C, na.rm = TRUE)), 
                        by = .(grid_id, year)]

# 4. Patch Coastal Climate Gaps (The 92+ cells)
missing_climate_sites <- setdiff(all_mainland_ids, unique(era5_yearly$grid_id))

if(length(missing_climate_sites) > 0) {
  cat("Found", length(missing_climate_sites), "sites missing climate data. Patching...\n")
  grid_has_data <- my_grid_mainland[my_grid_mainland$grid_id %in% unique(era5_yearly$grid_id), ]
  grid_missing  <- my_grid_mainland[my_grid_mainland$grid_id %in% missing_climate_sites, ]
  
  nearest_idx <- st_nearest_feature(st_centroid(grid_missing), st_centroid(grid_has_data))
  impute_map <- data.table(missing_id = grid_missing$grid_id, donor_id = grid_has_data$grid_id[nearest_idx])
  
  patched_era5 <- merge(era5_yearly, impute_map, by.x = "grid_id", by.y = "donor_id", allow.cartesian = TRUE)
  patched_era5[, grid_id := missing_id][, missing_id := NULL]
  
  era5_yearly <- rbindlist(list(era5_yearly, patched_era5), use.names = TRUE)
  cat("Successfully patched climate for the full mainland grid.\n")
}

# 5. MERGE WITH BIOLOGY (Keeping the unvisited 'Grey Cells')
# all.y = TRUE forces the unvisited cells into the final dataset
final_dt <- merge(agg_dt, era5_yearly, by = c("grid_id", "year"), all.y = TRUE)

# 6. Fill zeros for the "Grey Cells"
# Since we only have Sympetrum danae in our final_dt now:
final_dt[is.na(num_visits), num_visits := 0]
final_dt[is.na(`Sympetrum danae`), `Sympetrum danae` := 0]

cat("Success: final_dt ready with", uniqueN(final_dt$grid_id), "connected cells.\n")

# ==============================================================================
## 2b. DYNAMIC ERA5 CLIMATE - visualization (SOCIALLY FILTERED) ####
# ==============================================================================
# This plot will now match your final timeline perfectly
eu_temp_trend <- era5_yearly[, .(avg_temp = mean(mean_temp, na.rm = TRUE)), by = year]
library(ggplot2)
ggplot(eu_temp_trend, aes(x = year, y = avg_temp)) +
  geom_line(color = "darkorange", linewidth = 1) +
  geom_point(color = "darkorange", size = 2) +
  geom_smooth(method = "lm", color = "black", linetype = "dashed", se = FALSE) +
  theme_minimal() +
  labs(
    title = "Thermal Evolution of Surveyed Sites",
    x = "Year", y = "Temperature (°C)"
  )

# 5. Extract the CORRECT Warming Rate for the Burn-in
temp_lm <- lm(avg_temp ~ year, data = eu_temp_trend)
warming_rate <- coef(temp_lm)["year"] 
cat("Calculated warming rate for study sites:", round(warming_rate, 4), "°C/year\n")


# Immediately after creating my_grid_mainland
my_grid_mainland <- my_grid_mainland[order(as.character(my_grid_mainland$grid_id)), ]


# ==============================================================================
# 3. BUILD TEMPLATE (Synchronized 1995-2024 window) ####
# ==============================================================================
cat("\n--- 3. Building Spatio-Temporal Template (1,835 Cells) ---\n")

# BRIDGE: Ensure all subsequent blocks use the 1,835 mainland IDs
valid_sites <- all_mainland_ids 
valid_sites <- sort(valid_sites)
n_sites <- length(valid_sites)

# 1. Setup Years (Explicitly 1995-2024)
obs_years <- 2000:2024
burn_years <- 1995:1999
all_years <- c(burn_years, obs_years)

# CJ is the data.table version of expand.grid - much faster
template_all <- CJ(grid_id = valid_sites, year = all_years)

# 2. Calculate Early Mean (2000-2004) for the hindcast baseline
early_mean_dt <- era5_yearly[year %between% c(2000, 2004), 
                             .(early_avg = mean(mean_temp, na.rm = TRUE)), 
                             by = grid_id]

# 3. Apply Hindcast Logic
template_all <- merge(template_all, early_mean_dt, by = "grid_id", all.x = TRUE)
template_all[year %in% burn_years, mean_temp := early_avg - (warming_rate * (2002 - year))]

# 4. Merge Observed ERA5 Data
obs_temps_to_merge <- era5_yearly[, .(grid_id, year, obs_t = mean_temp)]
template_all <- merge(template_all, obs_temps_to_merge, by = c("grid_id", "year"), all.x = TRUE)

# 5. Finalize and Cleanup
template_all[!is.na(obs_t), mean_temp := obs_t]
template_all[, c("early_avg", "obs_t") := NULL]
setorder(template_all, year, grid_id)

cat("Template created with", nrow(template_all), "total rows for", n_sites, "cells.\n")

# ------------------------------------------------------------------
## 3b. FINAL VISUALIZATION CHECK ####
# ------------------------------------------------------------------
full_trend <- template_all[, .(avg_temp = mean(mean_temp, na.rm=TRUE)), by = year]
full_trend[, Period := ifelse(year < 2000, "Burn-in (Ramp)", "Observed")]
# Add labels for the plot
full_trend[, Category := ifelse(year < 2000, "Burn-in (Hindcast)", 
                                ifelse(year <= 2004, "Baseline (2000-2004)", "Observed (2005-2024)"))]

ggplot(full_trend, aes(x = year, y = avg_temp)) +
  geom_smooth(data = full_trend[Period == "Observed"], method = "lm", 
              formula = y ~ x, color = "black", linetype = "dashed", se = FALSE) +
  geom_line(aes(color = Period), linewidth = 1) +
  geom_point(aes(color = Period), size = 2) +
  scale_color_manual(values = c("Burn-in (Ramp)" = "grey50", "Observed" = "darkorange")) +
  theme_minimal() +
  labs(title = "Final Thermal Timeline (Synchronized)",
       subtitle = paste0("Hindcast anchored to 2000-2024 trend: ", round(warming_rate, 4), "°C/yr"),
       y = "Mean Temp (°C)", x = "Year")

# ------------------------------------------------------------------
## 3b. VISUALIZE ANCHORED HINDCAST WITH ITS OWN REGRESSION ####
# ------------------------------------------------------------------
cat("\n--- Visualizing Hindcast Alignment ---\n")

# 1. Calculate the values for the "Hindcast Regression Line"
# This is the line: y = EarlyMean + (warming_rate * (Year - 2002))
early_grand_avg <- mean(early_mean_dt$early_avg) # Average across all sites
full_trend[, Hindcast_Line := early_grand_avg + (warming_rate * (year - 2002))]

# 2. Plot
ggplot(full_trend, aes(x = year, y = avg_temp)) +
  # A. The overall 2000-2024 trend (Dashed Black)
  geom_smooth(data = full_trend[year >= 2000], method = "lm", 
              formula = y ~ x, color = "black", linetype = "dashed", se = FALSE) +
  
  # B. The specific Hindcast Regression used for burn-in (Solid Purple)
  geom_line(aes(y = Hindcast_Line), color = "purple", linewidth = 0.8, alpha = 0.6) +
  
  # C. The Actual Data (Grey/Blue/Orange)
  geom_line(aes(color = Category), linewidth = 1) +
  geom_point(aes(color = Category), size = 2) +
  
  scale_color_manual(values = c(
    "Burn-in (Hindcast)" = "grey50", 
    "Baseline (2000-2004)" = "royalblue", 
    "Observed (2005-2024)" = "darkorange"
  )) +
  theme_minimal() +
  labs(
    title = "Thermal Timeline: Validation of Hindcast Slope",
    subtitle = "Purple line = The warming rate applied backward from the 2000-2004 mean",
    y = "Mean Temp (°C)", x = "Year"
  ) +
  theme(legend.position = "bottom")

# ==============================================================================
# 4. TEMPERATURE SPLINES (Monotonic Indexing) ####
# ==============================================================================
cat("\n--- 4. Generating Thermal Gradient and Indices ---\n")

# 1. Define the limits with a safe buffer (e.g., 5 degrees beyond observed data)
T_min <- floor(min(template_all$mean_temp, na.rm = TRUE)) - 5
T_max <- ceiling(max(template_all$mean_temp, na.rm = TRUE)) + 5

# 2. Create the clean gradient (by 0.1 °C)
thermal_gradient <- seq(T_min, T_max, by = 0.1)
N_thermal_gradient <- length(thermal_gradient)

# 3. Generate the Basis Functions for the gradient
N_thermal_bf <- 10
library(splines)
thermal_bf <- bs(thermal_gradient, df = N_thermal_bf, intercept = FALSE)

# 4. Map your actual site-year temperatures to the index
temp_indices <- round((template_all$mean_temp - T_min) * 10) + 1
temp_indices_vec <- as.integer(temp_indices)


# --- PRE-INTEGRATING THERMAL BASIS FUNCTIONS ---
cat("\n--- Pre-integrating Splines for Speed ---\n")

# Use the 0.1°C step from your thermal_gradient
dx <- 0.1 

# Turn the 'hills' (B-splines) into 'ramps' (Integrated B-splines)
# integrated_bf <- apply(thermal_bf, 2, function(x) cumsum(x) * dx)

# Quick check: The last value of the first integrated spline should be > 0
# print(tail(integrated_bf[,1]))

# ==============================================================================
## 4b. VISUALIZE THE B-SPLINE BASIS FUNCTIONS ####
# ==============================================================================
library(ggplot2)
library(tidyr)
library(dplyr)

cat("\n--- Visualizing the B-Spline Basis Functions ---\n")

# 1. Convert the basis function matrix to a dataframe
df_bf <- as.data.frame(as.matrix(thermal_bf))

# Name the columns so we know which spline is which
colnames(df_bf) <- paste0("Basis_", 1:N_thermal_bf)

# 2. Add the exact temperature gradient
df_bf$Temperature <- thermal_gradient

# 3. Pivot to long format for ggplot
df_bf_long <- df_bf %>%
  pivot_longer(
    cols = starts_with("Basis_"),
    names_to = "Basis_Function",
    values_to = "Value"
  ) %>%
  # Extract the number so it sorts correctly in the legend
  mutate(Basis_Index = factor(as.integer(gsub("Basis_", "", Basis_Function))))

# 4. Build the Plot
p_splines <- ggplot(df_bf_long, aes(x = Temperature, y = Value, color = Basis_Index)) +
  geom_line(linewidth = 1.2, alpha = 0.8) +
  scale_color_viridis_d(option = "turbo", name = "Spline\nIndex") +
  theme_minimal() +
  labs(
    title = "B-Spline Basis Functions for Thermal Niche",
    subtitle = paste0(N_thermal_bf, " overlapping functions generated from ", T_min, "°C to ", T_max, "°C"),
    x = "Mean Annual Temperature (°C)",
    y = "Basis Function Strength"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 14, face = "bold")
  )

print(p_splines)


# ==============================================================================
# 5. OBSERVATION VECTORS & INDEX MAPPING ####
# ==============================================================================
cat("\n--- 5. Preparing Observation Data and Index Mapping ---\n")

template_obs <- template_all[year %in% obs_years]
setorder(template_obs, year, grid_id)

obs_data <- merge(template_obs, 
                  final_dt[, .(grid_id, year, num_visits, `Sympetrum danae`)], 
                  by = c("grid_id", "year"), 
                  all.x = TRUE)

obs_data[is.na(num_visits), num_visits := 0]
obs_data[is.na(`Sympetrum danae`), `Sympetrum danae` := 0]

setorder(obs_data, year, grid_id) 

# Explicit indexing
obs_data[, obs_state_id := match(grid_id, valid_sites)]
obs_data[, obs_time_id := year - min(obs_years) + 1]

# state mapping
obs_data[, state_idx := (obs_time_id - 1) * n_sites + obs_state_id]

# ============================================================================== #
# 6. SPATIAL BASES (20 Habitat Splines + 40 Spectral Dispersal Modes) ####
# ============================================================================== #
cat("\n--- 6. Computing Hybrid Spatial Bases ---\n")

# --- 6a. HABITAT SPLINES (With Boundary Buffer) ---
# Force everything to match your 'valid_sites' order
my_grid_mainland <- my_grid_mainland[match(valid_sites, my_grid_mainland$grid_id), ]
centroids_raw <- st_coordinates(st_centroid(my_grid_mainland))

cat("Generating padded spatial basis functions...\n")
set.seed(42)
N_target_bf <- 20
coords_raw <- st_coordinates(st_centroid(my_grid_mainland))

# 1. Create a Bounding Box that is 10% larger than your study area
# This is the "Extra" that manages boundary effects
buffer_pct <- 0.10
xlim <- range(coords_raw[,1]) + c(-1, 1) * diff(range(coords_raw[,1])) * buffer_pct
ylim <- range(coords_raw[,2]) + c(-1, 1) * diff(range(coords_raw[,2])) * buffer_pct

# 2. Create a dense grid of potential knots over this expanded box
# We start with more than 20 and filter down
grid_size <- 10 # 10x10 = 100 potential knots
potential_knots <- expand.grid(
  X = seq(xlim[1], xlim[2], length.out = grid_size),
  Y = seq(ylim[1], ylim[2], length.out = grid_size)
)

# 3. Filter knots to only those near your study area (within 150km)
# This keeps the "Extra" but removes knots in the deep Atlantic
potential_sf <- st_as_sf(potential_knots, coords = c("X", "Y"), crs = st_crs(my_grid_mainland))
dist_to_mainland <- st_distance(potential_sf, st_union(my_grid_mainland))
keep_idx <- which(as.numeric(dist_to_mainland) < 150000) # 150km buffer

# 4. Use K-means to reduce these "relevant" knots to exactly 20
final_knots_sf <- potential_sf[keep_idx, ]
km_final <- kmeans(st_coordinates(final_knots_sf), centers = N_target_bf)
knots <- km_final$centers # Final 20 knots with boundary coverage

# 5. Build the Distance Matrices for Stan
dist_sk <- st_distance(st_centroid(my_grid_mainland), 
                       st_as_sf(as.data.frame(knots), coords = c("X", "Y"), crs = st_crs(my_grid_mainland)))

# STRIP UNITS and preserve matrix shape (Rows = Grid Cells, Cols = 20 Knots)
dist_mat_raw <- matrix(as.numeric(dist_sk), 
                       nrow = nrow(my_grid_mainland), 
                       ncol = nrow(knots)) / 1000 

dist_mat_anchors <- as.matrix(dist(knots)) / 1000

# --- Convert Linear Distances to Smooth Radial Basis Functions (RBF) ---
# Calculate a bandwidth (sigma) based on the distance between knots
bandwidth <- median(dist_mat_anchors[dist_mat_anchors > 0]) 

# Apply the Gaussian Kernel formula (No units = No errors!)
spatial_bf <- exp(-(dist_mat_raw^2) / (2 * bandwidth^2))
N_spatial_bf <- ncol(spatial_bf)

cat("GP Basis computed using RBF bandwidth of", round(bandwidth, 2), "km\n")

# --- 6b. CSR DISPERSAL (The "How they move") ---
library(Matrix)
cat("\n--- Generating CSR Network (200km cutoff) ---\n")

dist_mat_disp <- as.matrix(st_distance(st_centroid(my_grid_mainland)))
dist_mat_disp <- units::drop_units(dist_mat_disp) / 1000 

threshold_km <- 200.0
dist_mat_disp[dist_mat_disp > threshold_km] <- 0 
diag(dist_mat_disp) <- 0 

sp_mat_gen <- as(drop0(dist_mat_disp), "generalMatrix")
sp_mat_R <- as(sp_mat_gen, "RsparseMatrix")

# Extract the vectors for Row-Stochastic Normalization
dists <- sp_mat_R@x
to_idx <- sp_mat_R@j + 1                   # Column indices
row_ptr <- sp_mat_R@p + 1                  # Row pointers
row_ids <- rep(1:n_sites, diff(sp_mat_R@p)) # Explicit Row IDs
E_disp <- length(dists)

cat("Generated CSR Network with", E_disp, "edges.\n")

# --- Re-format the temperature vector into a T_total x S matrix ---
# This ensures that temp_idx_mat[t, s] is the temperature for year t, site s
temp_idx_mat <- matrix(temp_indices_vec, nrow = length(all_years), ncol = n_sites, byrow = TRUE)


# ============================================================================== #
# 7. ASSEMBLE STAN DATA LIST (Final Production Version) ####
# ============================================================================== #
cat("\n--- 7. Assembling Final Stan Data List ---\n")

# Calculate squared distances in R
dist_mat_anchors_sq <- dist_mat_anchors^2

stan_data_real <- list(
  S = n_sites,
  T_total = length(all_years),
  T_obs = length(obs_years),
  N_obs_total = nrow(obs_data),
  
  # Observations (Y <= K check passed earlier)
  N_state_total = n_sites * length(obs_years), # [ADDED]
  Y_vec = as.integer(obs_data$`Sympetrum danae`),
  K_visits_vec = as.integer(obs_data$num_visits),
  idx_seen = which(as.integer(obs_data$`Sympetrum danae`) > 0),
  idx_zero = which(as.integer(obs_data$`Sympetrum danae`) == 0),
  N_seen = sum(as.integer(obs_data$`Sympetrum danae`) > 0),
  N_zero = sum(as.integer(obs_data$`Sympetrum danae`) == 0),
  map_state_idx = obs_data$state_idx,
  
  # Habitat Splines
  N_spatial_bf = N_spatial_bf,
  spatial_bf = spatial_bf,
  dist_mat_anchors = dist_mat_anchors,
  
  # --- THE CSR REPLACEMENT ---
  E_disp = E_disp,
  dists = dists,
  to_idx = to_idx,        # Renamed
  row_ptr = row_ptr,
  row_ids = row_ids,      # New addition
  dist_mat_anchors_sq = dist_mat_anchors_sq,
  # ---------------------------
  
  # Thermal Niche (Matching the Parametric Double Logistic model)
  N_thermal_gradient = N_thermal_gradient,
  thermal_gradient = thermal_gradient,
  temp_idx = as.vector(t(temp_idx_mat))
)

saveRDS(stan_data_real, "stan_data_ZwaHei_HybridIDE.rds")
cat("Data list successfully saved for Hybrid IDE.\n")

# ==============================================================================
# 8. FINAL AUDIT ####
# ==============================================================================
# Audit 1: Smooth Surface Check
set.seed(123)
fake_weights <- rnorm(N_spatial_bf)
my_grid_mainland$Synthetic_Habitat <- as.vector(spatial_bf %*% fake_weights)

ggplot(my_grid_mainland) +
  geom_sf(aes(fill = Synthetic_Habitat), color = NA) +
  scale_fill_viridis_c(option = "magma") +
  theme_minimal() +
  labs(title = "Habitat Spline Audit", subtitle = "Must look smooth and cover all of Europe")

# Audit 2: Knot Check
anchors_df <- as.data.frame(knots)
colnames(anchors_df) <- c("X", "Y")

ggplot() +
  geom_sf(data = my_grid_mainland, color = "grey90", fill = NA) +
  geom_point(data = anchors_df, aes(x=X, y=Y), color="red", shape=4, size=3) +
  theme_void() + labs(title = "Knot Placement Audit (Exactly 20)")


cat("--- Temperature Data Audit ---\n")
cat("Absolute Min Temp:", min(template_all$mean_temp, na.rm=TRUE), "°C\n")
cat("Absolute Max Temp:", max(template_all$mean_temp, na.rm=TRUE), "°C\n")
cat("Median Temp:", median(template_all$mean_temp, na.rm=TRUE), "°C\n")


# ============================================================================== #
# 8. PREPARE STAN (DIRECT MCMC) ####
# ============================================================================== #
# IMPORTANT: Make sure this filename matches your newly saved NCP Stan script!
stan_GP_model <- cmdstan_model("SpatGaus_DynamicSCAM_SpatTempInt.stan")

make_init_stable <- function() {
  list(
    r = runif(1, 0.4, 0.6),
    log_alpha = log(runif(1, 10, 30)), 
    logit_p = runif(1, -2.5, -1.5),
    
    # --- NCP Spline Inits (thermal_sd removed!) ---
    thermal_intercept = rnorm(1, -2, 0.5), 
    weight_start_raw = rnorm(1, 0, 0.1),
    weight_steps_raw = rnorm(stan_data_real$N_thermal_bf - 1, 0, 0.1),
    
    # --- Temporal & Spatial ---
    sigma_year = runif(1, 0.005, 0.02),
    eps_year_raw = rnorm(stan_data_real$T_total, 0, 0.01),
    N0_proportion = runif(1, 0.4, 0.6),                        
    rho_phi = runif(1, 0.4, 0.6),
    sigma_phi = runif(1, 0.1, 0.3),
    phi_eta = rnorm(stan_data_real$N_spatial_bf, 0, 0.01),
    
    # --- Spatiotemporal Interaction ---
    sigma_gamma = runif(1, 0.01, 0.05), 
    gamma_eta = rnorm(stan_data_real$N_spatial_bf, 0, 0.01)
  )
}

cat("\n--- Alignment Check 3: Initialization Values ---\n")
test_inits <- make_init_stable()
cat("NCP Parameters present in Inits:", all(c("weight_start_raw", "weight_steps_raw") %in% names(test_inits)), "\n")

# Check for existence of GP parameters
gp_params <- c("rho_phi", "sigma_phi", "phi_eta")
gp_check <- all(gp_params %in% names(test_inits))
cat("GP Parameters present in Inits:", gp_check, "\n")

# Check for removal of ICAR parameters (should be FALSE)
icar_check <- "phi_raw" %in% names(test_inits)
cat("Old ICAR parameters correctly removed from Inits:", !icar_check, "\n")



# The max distance between anchors should be roughly 2 (since we scaled to -1 to 1)
anchor_mat <- do.call(rbind, stan_data_real$spatial_bf_range)
max_dist <- max(dist(anchor_mat))
cat("Maximum distance between GP anchors (should be ~2.0 to 2.8):", round(max_dist, 2), "\n")



fit_GP_Test <- stan_GP_model$sample(
  data = stan_data_real,
  chains = 1,              
  iter_warmup = 10,       # Extremely tiny
  iter_sampling = 5,      # Extremely tiny
  init = make_init_stable, 
  refresh = 1,            # Print every step
  max_treedepth = 12,
  adapt_delta = 0.85
)

# ------------------------------------------------------------------------------ #
# STEP B: RUN MCMC (Ready for your final run!)
# ------------------------------------------------------------------------------ #
cat("\n--- Running Final MCMC (Standard NUTS) ---\n")

output_folder <- "."

fit_spectral_final <- stan_GP_model$sample(
  data = stan_data_real,
  chains = 6,
  parallel_chains = 6,   
  iter_warmup = 300,
  iter_sampling = 300,
  init = make_init_stable, 
  refresh = 50,
  max_treedepth = 12,
  adapt_delta = 0.85,
  save_cmdstan_config = TRUE,
  output_dir = output_folder,
  output_basename = "bruinrode_run"
)


# 1. Safe entire model object
fit_spectral_final$draws()
fit_spectral_final$metadata()
saveRDS(fit_spectral_final, "fit_evo.rds")

# 2. safe draws seperately
library(posterior)
draw_base <- fit_spectral_final$draws()
saveRDS(draw_base, "draw_evo.rds")






## Read in data (from supercomputer)
library(cmdstanr)
library(posterior)

# 1. Look inside the folder where you manually extracted the files
# (Adjust "cmdstan_outputs" if you named the folder something else!)
csv_files <- list.files("Thesis_Data", full.names = TRUE, pattern = "ZwaHei_Hybrid_v2-.*\\.csv$")

# 2. Sanity check: Make sure R actually sees all 6 chains!
print(csv_files)

# 3. Rebuild the model object
fit_spectral_final <- as_cmdstan_fit(csv_files)

# 4. Read your draws (ensure the path matches where you saved it)
# combined_draws <- readRDS("Thesis_Data/fit_spectral_final_ZwaHei_Parametric_ST.rds")
combined_draws <- fit_spectral_final$draws()



## ##
library(cmdstanr)
library(posterior)
library(ggplot2)

# 2. Check basic parameter health (Look for Rhat < 1.05)
cat("\n--- Parameter Health Check ---\n")
print(fit_spectral_final$summary(
  variables = c("tau_raw", "thermal_max", "slope_L", "slope_R", "rho_phi", "sigma_phi"),
  measures = c("mean", "sd", "rhat", "ess_bulk")
))

# 3. Extract the Niche parameters to check the ceiling
draws_df <- as_draws_df(fit_spectral_final$draws(variables = c(
  "thermal_max", "T_nodes[1]", "T_nodes[2]", "slope_L", "slope_R"
)))

plot_temps <- seq(-10, 30, length.out = 100)
plot_data <- data.frame()

# Stan's soft limit math from cap_val = 8.0
cap_val <- 10.0
soft_limit <- exp(cap_val)

for (temp in plot_temps) {
  # Calculate raw K
  log_K_raw <- draws_df$thermal_max + 
    plogis(draws_df$slope_L * (temp - draws_df$`T_nodes[1]`), log = TRUE) + 
    plogis(draws_df$slope_R * (draws_df$`T_nodes[2]` - temp), log = TRUE)
  
  # Apply the soft limit EXACTLY as Stan does
  K_vals <- soft_limit * plogis(log_K_raw - cap_val)
  
  plot_data <- rbind(plot_data, data.frame(
    Temperature = temp,
    K_median = median(K_vals),
    K_lower = quantile(K_vals, 0.025),
    K_upper = quantile(K_vals, 0.975)
  ))
}

# 4. Plot the Niche against the Ceiling
ggplot(plot_data, aes(x = Temperature)) +
  geom_ribbon(aes(ymin = K_lower, ymax = K_upper), fill = "purple", alpha = 0.3) +
  geom_line(aes(y = K_median), color = "purple", linewidth = 1.2) +
  geom_hline(yintercept = soft_limit, linetype = "dashed", color = "red", linewidth = 1) +
  annotate("text", x = 0, y = soft_limit - 100, label = "cap_val = 8.0 Ceiling", color = "red") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Posterior Thermal Niche (Real Data)",
    subtitle = "Checking for 'Flat-Topping' against the soft limit",
    x = "Temperature (°C)", y = "Carrying Capacity (K)"
  )

## ##


library(cmdstanr)
library(ggplot2)
library(posterior)

# 1. Extract the latent abundance samples
# 'fit' is your cmdstanr fit object
N_samples <- fit$draws("N_state_vec", format = "df")

# 2. Reshape into a simple vector of all estimates across all space/time
# This removes the chain/iteration metadata to focus on the values
all_estimates <- as.numeric(unlist(N_samples[, -c(1:2)]))

# 3. Calculate the proportion hitting the bound
upper_bound <- 10000
prop_at_bound <- sum(all_estimates >= upper_bound) / length(all_estimates)

cat("Percentage of cells/years at the upper bound (10,000):", 
    round(prop_at_bound * 100, 4), "%\n")

# 4. Visualize the distribution
ggplot(data.frame(N = all_estimates), aes(x = N)) +
  geom_histogram(bins = 100, fill = "steelblue", color = "white") +
  geom_vline(xintercept = upper_bound, linetype = "dashed", color = "red") +
  labs(title = "Distribution of Latent Abundance Estimates",
       subtitle = paste("Vertical line marks the 10,000 cap. Percentage at cap:", 
                        round(prop_at_bound * 100, 3), "%"),
       x = "Estimated Abundance (N)",
       y = "Frequency") +
  theme_minimal()


# ==============================================================================
# QUICK CONVERGENCE CHECK SCRIPT ####
# ==============================================================================
library(bayesplot)
library(ggplot2)
library(cmdstanr)

cat("\n--- 1. Statistical Diagnostic Summary ---\n")

# 1. Print core ecological parameters
# Focus on Growth (r), Dispersal (alpha), Detection (theta), and Philopatry (rho) (now removed)
core_params <- c(
  # --- 1. Global Ecological Rates ---
  "r",                # Intrinsic growth rate
  "tau_raw",          # Dispersal (Spectral Diffusion coefficient)
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
  
  # --- 4. The "Health" Metric ---
  "lp__"              # Log-posterior (Crucial for detecting geometry issues)
)
print(fit_spectral_final$summary(core_params))

# 2. Check for HMC-specific issues (Divergences)
# High divergences suggest the model geometry is too complex for the current step size
fit_spectral_final$diagnostic_summary()

cat("\n--- 2. Visual Traceplot Inspection ---\n")

# 3. Plot traceplots
# Looking for "fuzzy caterpillars" where all 4 chains overlap and stay stationary
color_scheme_set("viridis")
mcmc_trace(fit_spectral_final$draws(c(core_params,"lp__"))) +
  labs(title = "Traceplots: Core Ecological Parameters",
       subtitle = "Chains should be well-mixed and stationary") +
  theme_minimal()

cat("\n--- 3. Numerical Convergence Thresholds ---\n")

# 4. Extract summary to check R-hat and ESS specifically
fit_summ <- fit_spectral_final$summary(core_params)

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
# Check log-likelihood
# Compare log-posterior across chains
library(bayesplot)
library(ggplot2)

lp <- fit_spectral_final$draws("lp__")
mcmc_trace(lp) + 
  labs(title = "Log-Posterior Traceplot", subtitle = "Higher is better")

library(posterior)

# Extract draws and subset to chains
combined_draws <- fit_spectral_final$draws()
final_draws <- subset_draws(combined_draws, chain = c(1,2,3,4,5,6))

# Re-check diagnostics on the subset
# summarise_draws(final_draws, "rhat", "ess_bulk")

# Plot the trace for the key biological parameters and the log-posterior
color_scheme_set("viridis")
arr <- as_draws_array(final_draws)
subset_arr <- arr[, , core_params]

mcmc_trace(subset_arr)

mcmc_pairs(subset_arr)
# Correlation between alpha and sigma_phi



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


# ==============================================================================
# THE POST-PROCESSING "BRIDGE" ####
# ==============================================================================
cat("\n--- Surgically Extracting Latent Abundance ---\n")

library(cmdstanr)

# 1. Define the missing variables from your dataset
T_obs <- 25      # 25 observation years (2000 to 2024)
N_sites <- 1835  # The number of connected grid cells

# 2. Re-locate the raw CSV files

# 3. Read ONLY the massive abundance vector to save RAM
N_fit <- read_cmdstan_csv(csv_files, variables = c("N_state"))
N_mat_clean <- N_fit$post_warmup_draws

# 4. Create the Median Vector (Collapsing 1800 draws down to 1 median per site-year)
# Note: we use apply over dimension 3 because post_warmup_draws is an array of [iterations, chains, parameters]
# Flattening it first makes apply easier:
N_mat_flat <- as_draws_matrix(N_mat_clean)
N_mat_vector <- apply(N_mat_flat, 2, median)

# 5. NOW we build the matrix! (T_obs rows, N_sites cols)
N_mat <- matrix(N_mat_vector, nrow = T_obs, ncol = N_sites, byrow = TRUE)

# 6. Delete the massive objects and run Garbage Collection to free RAM!
rm(N_fit, N_mat_clean, N_mat_flat)
gc()

# 7. Ensure grid variables are defined for the spatial plots
my_grid_final <- my_grid_mainland
grid_y_coords <- st_coordinates(st_centroid(my_grid_final))[,2]

cat("Success! N_mat built and RAM cleared.\n")

# ==============================================================================
# 9. POSTERIOR PREDICTIVE CHECKS (PPC) ####
# ==============================================================================
cat("\n--- 9. Posterior Predictive Check ---\n")

library(posterior)
library(bayesplot)
library(ggplot2)

# 1. Extract the generated y_rep directly from the Stan output
y_rep_draws <- fit_spectral_final$draws("y_rep")
y_rep_mat <- as_draws_matrix(y_rep_draws)
y_obs_vec <- stan_data_real$Y_vec

# 2. Subsample 100 draws for plotting (Standard practice)
set.seed(123)
y_rep_subset <- y_rep_mat[sample(nrow(y_rep_mat), 100), , drop = FALSE]

# Plot A: Rootogram (Best for count data)
p_ppc_root <- ppc_rootogram(y_obs_vec, y_rep_subset) +
  labs(title = "Posterior Predictive Check: Rootogram",
       subtitle = "Comparing observed counts to model-simulated counts") +
  theme_minimal()

# Plot B: Proportion of Zeros (Tests the Zero-Inflation psi parameter)
p_ppc_zeros <- ppc_stat(y_obs_vec, y_rep_subset, stat = function(y) mean(y == 0)) +
  labs(title = "Posterior Predictive Check: Proportion of Zeros",
       subtitle = "Did the model capture the zero-inflation accurately?") +
  theme_minimal()

print(p_ppc_root)
print(p_ppc_zeros)


# ==============================================================================
# 9b. POSTERIOR PREDICTIVE CHECKS (GENERATED IN R) ####
# ==============================================================================
cat("\n--- 9. Simulating Posterior Predictive Distribution ---\n")
library(posterior)
library(bayesplot)
library(ggplot2)

# 1. Grab the observation data from your original Stan Data List
Y_obs <- stan_data_real$Y_vec
K_visits <- stan_data_real$K_visits_vec
map_idx <- stan_data_real$map_state_idx

# 2. Extract Detection Probability (p) and convert from logit to probability
# We take a random sample of 100 draws to save RAM and plotting time
set.seed(42)
n_ppc_draws <- 100
total_draws <- fit_spectral_final$metadata()$iter_sampling * 6 # 6 chains
sample_rows <- sample(1:total_draws, n_ppc_draws)

logit_p_draws <- as_draws_matrix(fit_spectral_final$draws("logit_p"))[sample_rows, 1]
p_draws <- plogis(logit_p_draws)

# 3. Extract the Latent Abundance for those specific 100 draws
cat("Extracting N_state subset (This may take a few seconds)...\n")
N_draws <- as_draws_matrix(fit_spectral_final$draws("N_state"))[sample_rows, ]

# 4. Simulate the Observation Process
y_rep_mat <- matrix(NA, nrow = n_ppc_draws, ncol = length(Y_obs))

for(i in 1:n_ppc_draws) {
  # A. Map the site-year abundance to the specific field visits
  N_for_obs <- N_draws[i, map_idx]
  
  # B. Calculate Occupancy Probability: psi = 1 - exp(-N)
  psi <- 1 - exp(-N_for_obs)
  
  # C. Simulate True Presence/Absence (Z)
  Z <- rbinom(length(Y_obs), size = 1, prob = psi)
  
  # D. Simulate Observed Counts given effort (K) and detection (p)
  y_rep_mat[i, ] <- rbinom(length(Y_obs), size = K_visits, prob = p_draws[i] * Z)
}

# Clear the subset matrix to keep RAM happy
rm(N_draws)
gc()

cat("y_rep generated successfully! Plotting...\n")

# 5. Plot the PPCs
color_scheme_set("brightblue")

# A. Density Overlay (Does the simulated data match the overall shape of real data?)
p1 <- ppc_dens_overlay(y = Y_obs, yrep = y_rep_mat) +
  coord_cartesian(xlim = c(0, max(Y_obs) + 5)) + # Zoom in on the meaningful range
  labs(title = "PPC: Density Overlay", subtitle = "Dark line = Real Data, Light lines = Model Simulations")

# B. Proportion of Zeros (Did the model capture the zero-inflation correctly?)
p2 <- ppc_stat(y = Y_obs, yrep = y_rep_mat, stat = function(y) mean(y == 0)) +
  labs(title = "PPC: Proportion of Zeros", subtitle = "Does the model predict the right amount of non-detections?")

# C. Maximum Count (Can the model capture extreme swarms?)
p3 <- ppc_stat(y = Y_obs, yrep = y_rep_mat, stat = "max") +
  labs(title = "PPC: Maximum Count")

# Display them
print(p1)
print(p2)
print(p3)

# ============================================================================== #
# 10. POSTERIOR DISTRIBUTIONS ####
# ============================================================================== #
cat("\n--- 10. Posterior Distributions ---\n")

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

# ==============================================================================
# 10b. POSTERIOR VS. PRIOR DISTRIBUTIONS (TRUE PRIORS) ####
# ==============================================================================
library(tidyr)
library(dplyr)
library(ggplot2)

cat("\n--- 10. Posterior vs. Prior Distributions ---\n")

# 1. Extract Posterior Draws
draws_df <- as_draws_df(subset_draws(final_draws, variable = core_params))
n_draws <- nrow(draws_df)

posterior_long <- draws_df |>
  as_tibble() |>
  select(all_of(core_params)) |>
  pivot_longer(everything(), names_to = "Parameter", values_to = "Value") |>
  mutate(Distribution = "Posterior")

# 2. Simulate Prior Draws exactly as defined in Stan
set.seed(123)
# Properly simulate the truncated normal for 'r' (matching Stan's <lower=0>)
# We generate extra draws, throw away the negatives, and keep exactly n_draws
r_raw <- rnorm(n_draws * 3, mean = 0.5, sd = 0.5)
r_prior <- r_raw[r_raw > 0][1:n_draws]

# Simulate Prior Draws exactly as defined in Stan
set.seed(123)
prior_df <- data.frame(
  # Stan: r ~ normal(0.5, 0.5)
  r = r_prior,        
  
  # Stan: log_alpha ~ normal(2.7, 1); (Transformed: alpha = exp(log_alpha))
  # Updated mean from 4 to 2.7 to match your new Stan script
  alpha = exp(rnorm(n_draws, mean = 2.7, sd = 1)),
  
  # Stan: N0_proportion ~ beta(2, 2)
  N0_proportion = rbeta(n_draws, shape1 = 2, shape2 = 2),
  
  # Stan: rho_phi ~ inv_gamma(5, 5)
  # R equivalent: 1 / Gamma(shape, rate)
  rho_phi = 1 / rgamma(n_draws, shape = 5, rate = 5),  
  
  # --- NEW THERMAL PARAMETERS ---
  # Stan: thermal_intercept ~ normal(0, 5)
  thermal_intercept = rnorm(n_draws, mean = 0, sd = 5),
  
  # Stan: thermal_sd ~ normal(0, 3) <lower=0>
  thermal_sd = abs(rnorm(n_draws, mean = 0, sd = 3)),
  
  # --- SPATIAL & TEMPORAL ---
  # Stan: sigma_phi ~ normal(0, 0.5) <lower=0>
  sigma_phi = abs(rnorm(n_draws, mean = 0, sd = 0.5)),        
  
  # Stan: sigma_year ~ normal(0, 0.4) <lower=0>
  sigma_year = abs(rnorm(n_draws, mean = 0, sd = 0.4))           
)
prior_long <- prior_df |>
  pivot_longer(everything(), names_to = "Parameter", values_to = "Value") |>
  mutate(Distribution = "Prior")

# 3. Combine them into one dataset
combined_long <- bind_rows(prior_long, posterior_long)

# 4. Calculate Medians ONLY for the Posterior text labels
estimates <- posterior_long |>
  group_by(Parameter) |>
  summarize(Median = median(Value))

# 5. Plot the Overlay
ggplot(combined_long, aes(x = Value, fill = Distribution, color = Distribution)) +
  geom_density(alpha = 0.5, linewidth = 0.8) +
  geom_text(data = estimates,
            aes(x = Inf, y = Inf,
                label = paste0("Posterior Median: ", round(Median, 3))),
            hjust = 1.1, vjust = 1.1,
            inherit.aes = FALSE,
            color = "black",
            fontface = "bold") +
  facet_wrap(~Parameter, scales = "free") +
  scale_fill_manual(values = c("Posterior" = "skyblue", "Prior" = "grey70")) +
  scale_color_manual(values = c("Posterior" = "dodgerblue4", "Prior" = "grey40")) +
  theme_minimal() +
  labs(title = "Prior vs. Posterior Overlays",
       subtitle = "Checking for Prior Contraction (Did the data inform the model?)",
       y = "Density", x = "Parameter Value")

# ============================================================================== #
# 11. RECOVERING THE STATIC THERMAL NICHE ####
# ============================================================================== #
cat("\n--- 11. Recovering the Thermal Performance Curve (TPC) ---\n")

library(ggplot2)
library(posterior)

# 1. Extract the EXACT parameter Stan saved (the 100+ length thermal_niche vector)
vars_niche <- paste0("thermal_curve[", 1:N_thermal_gradient, "]")

# Extract the matrix (Rows = MCMC draws, Columns = temperature increments)
niche_mat <- as_draws_matrix(subset_draws(final_draws, variable = vars_niche))

# 2. Predict Carrying Capacity (K) 
# Math: Convert from log scale to natural scale
K_pred <- exp(niche_mat)

# 3. Summarize into a Plotting Dataframe
df_tpc <- data.frame(
  Temp = thermal_gradient,
  K_Median = apply(K_pred, 2, median),
  K_Lower  = apply(K_pred, 2, quantile, 0.025),
  K_Upper  = apply(K_pred, 2, quantile, 0.975)
)

# 4. Plot the Thermal Niche
p_tpc <- ggplot(df_tpc, aes(x = Temp, y = K_Median)) +
  geom_ribbon(aes(ymin = K_Lower, ymax = K_Upper), fill = "firebrick", alpha = 0.3) +
  geom_line(color = "firebrick", linewidth = 1.5) +
  coord_cartesian(ylim = c(0, quantile(df_tpc$K_Median, 0.99))) +
  theme_minimal() +
  labs(title = "Estimated Thermal Performance Curve",
       subtitle = "Shape-Constrained Monotonic Spline",
       y = "Relative Carrying Capacity (K)",
       x = "Mean Annual Temperature (°C)")

print(p_tpc)



# ============================================================================== #
# 12. RANGE DYNAMICS & VELOCITY (MEMORY-SAFE VERSION) ####
# ============================================================================== #
library(dplyr)
library(ggplot2)
library(cmdstanr)
library(posterior)

cat("\n--- 12. Range Velocity Analysis (2.5%, 25%, 50%, 75%, 97.5%) ---\n")

# --- ROBUST MARGIN FUNCTION (DRAWS VERSION) ---
get_margin_draw <- function(abundance_vec, y_coords, prob, noise_threshold = 0.1) {
  clean_n <- ifelse(abundance_vec < noise_threshold, 0, abundance_vec)
  if(sum(clean_n) == 0) return(NA)
  
  ord <- order(y_coords)
  y_sorted <- y_coords[ord]
  n_sorted <- clean_n[ord]
  
  cum_n <- cumsum(n_sorted) / sum(n_sorted)
  return(y_sorted[which(cum_n >= prob)[1]])
}

# 0. Safety Check: Ensure grid_y_coords exists in the environment
if(!exists("grid_y_coords")) {
  grid_y_coords <- st_coordinates(st_centroid(my_grid_mainland))[,2]
}

# 1. Temporarily bring the massive abundance matrix back into RAM
cat("Loading abundance draws for velocity calculation...\n")
# FIXED: Changed N_obs_vec to N_state to match the updated Stan model!
N_fit <- read_cmdstan_csv(csv_files, variables = c("N_state"))
N_mat_clean <- as_draws_matrix(N_fit$post_warmup_draws)

# 2. Setup Data and Target Quantiles
n_draws <- nrow(N_mat_clean)
N_array <- array(N_mat_clean, dim = c(n_draws, N_sites, T_obs))

quantiles_to_track <- c(0.025, 0.25, 0.50, 0.75, 0.975)
quantile_names <- c("Trailing (2.5%)", "25th Percentile", "Core Median (50%)", "75th Percentile", "Leading (97.5%)")

# Initialize a list to hold the margin arrays
margin_list <- list()
for(q_name in quantile_names) {
  margin_list[[q_name]] <- array(NA, dim = c(n_draws, T_obs))
}

# 3. Calculate Margins for EVERY DRAW across all 5 quantiles
cat("Calculating margins across all MCMC draws for 5 quantiles. This might take a minute...\n")
for(d in 1:n_draws) {
  for(t in 1:T_obs) {
    for(i in seq_along(quantiles_to_track)) {
      margin_list[[quantile_names[i]]][d, t] <- get_margin_draw(N_array[d, , t], grid_y_coords, quantiles_to_track[i], 0.1)
    }
  }
}

# 4. Delete the massive matrix to free RAM
rm(N_fit, N_mat_clean, N_array)
gc()
cat("RAM cleared! Proceeding to summarize and plot...\n")

# 5. Summarize into Medians and 95% CIs for Plotting
summarize_margins <- function(margin_array) {
  data.frame(
    Year = obs_years,
    Median = apply(margin_array, 2, median, na.rm=TRUE) / 1000,
    Lower = apply(margin_array, 2, quantile, 0.025, na.rm=TRUE) / 1000,
    Upper = apply(margin_array, 2, quantile, 0.975, na.rm=TRUE) / 1000
  )
}

# Combine all summarized quantiles into one dataframe
df_plot <- data.frame()
vel_stats <- list()

for(q_name in quantile_names) {
  df_temp <- summarize_margins(margin_list[[q_name]]) %>% mutate(Edge = q_name)
  df_plot <- bind_rows(df_plot, df_temp)
  
  # Calculate Observed Velocity Stats for each
  vels <- apply(margin_list[[q_name]], 1, function(y) {
    if(any(is.na(y))) return(NA)
    coef(lm(y/1000 ~ obs_years))[2] # Slope in km/year
  })
  vel_stats[[q_name]] <- c(Median = median(vels, na.rm=TRUE), quantile(vels, c(0.025, 0.975), na.rm=TRUE))
}

# Ensure correct factor ordering for the plot legend
df_plot$Edge <- factor(df_plot$Edge, levels = rev(quantile_names))

# ==============================================================================
# 12b. THEORETICAL VELOCITY (Fisher-KPP / IDE)
# ==============================================================================
r_draws <- as_draws_df(subset_draws(final_draws, variable = "r"))$r
alpha_draws <- as_draws_df(subset_draws(final_draws, variable = "alpha"))$alpha

# Approximated asymptotic spreading speed for Laplace Kernel
v_draws <- 2 * alpha_draws * sqrt(2 * r_draws)

v_theoretical <- c(Median = median(v_draws), quantile(v_draws, c(0.025, 0.975)))

# Print Combined Results
cat("\n--- FINAL VELOCITY SYNTHESIS (km/year) ---\n")
for(q_name in rev(quantile_names)) {
  cat(sprintf("%-20s: %5.2f [%5.2f, %5.2f]\n", 
              q_name, vel_stats[[q_name]][1], vel_stats[[q_name]][2], vel_stats[[q_name]][3]))
}
cat(sprintf("%-20s: %5.2f [%5.2f, %5.2f]\n", 
            "Theoretical (V)", v_theoretical[1], v_theoretical[2], v_theoretical[3]))

# ==============================================================================
# 12c. VISUALIZATION
# ==============================================================================
# Define a diverging color palette (Blue for leading, Red for trailing)
edge_colors <- c("Leading (97.5%)"   = "#0571b0",  # Dark Blue
                 "75th Percentile"   = "#92c5de",  # Light Blue
                 "Core Median (50%)" = "#999999",  # Grey
                 "25th Percentile"   = "#f4a582",  # Light Red
                 "Trailing (2.5%)"   = "#ca0020")  # Dark Red

p_margins <- ggplot(df_plot, aes(x = Year, y = Median, color = Edge, fill = Edge)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 1.5) +
  theme_minimal() +
  scale_color_manual(values = edge_colors) +
  scale_fill_manual(values = edge_colors) +
  labs(
    title = "Spatio-Temporal Range Dynamics",
    subtitle = paste0("Leading Edge Velocity: ", round(vel_stats[["Leading (97.5%)"]][1], 2), " km/yr"),
    y = "Latitude (km)",
    x = "Year",
    color = "Population Slice",
    fill = "Population Slice"
  ) +
  theme(legend.position = "right")

print(p_margins)

# ==============================================================================
# 13. LATENT ABUNDANCE GIF
# ==============================================================================
library(gganimate)

cat("\n--- 13. Generating Latent Abundance GIF ---\n")

coords <- as.data.frame(st_coordinates(st_centroid(my_grid_mainland)))
coords$grid_id <- valid_sites

df_n_mod <- data.frame(
  year = rep(obs_years, each = N_sites),
  grid_id = rep(valid_sites, T_obs),
  abundance = as.vector(t(N_mat))
) %>% left_join(coords, by = "grid_id")

mod_cap <- quantile(df_n_mod$abundance, 0.99)

p_gif <- ggplot(df_n_mod, aes(X, Y, fill = pmin(abundance, mod_cap))) +
  geom_tile() +
  scale_fill_viridis_c(option = "magma", name = "N") +
  theme_void() +
  transition_states(year) +
  labs(title = "Latent Abundance Year: {closest_state}")

animate(p_gif, renderer = gifski_renderer("latent_abundance.gif"))

# ============================================================================== #
## 13b. SIDE-BY-SIDE GIF: OBSERVED VS. ESTIMATED ABUNDANCE (ABSOLUTE SCALE) ####
# ============================================================================== #
library(dplyr)
library(gganimate)
library(gifski)
library(ggplot2)
library(sf)
library(data.table)

cat("\n--- 13b. Generating Side-by-Side Animation (Absolute) ---\n")

# 0. CRITICAL FIX: Aggregate the 54k observation rows back to 45k spatial states
obs_agg <- obs_data[, .(
  Total_Y = sum(`Sympetrum danae`, na.rm = TRUE),
  Total_K = sum(num_visits, na.rm = TRUE)
), by = .(grid_id, year)]
setorder(obs_agg, year, grid_id) # Ensure order matches the map perfectly

# 1. Prepare Spatial Coordinates
coords <- as.data.frame(st_coordinates(st_centroid(my_grid_mainland)))
coords$grid_id <- my_grid_mainland$grid_id

# 2. Extract Estimated Latent Abundance (N)
df_est <- data.frame(
  year = rep(obs_years, each = N_sites),
  grid_id = rep(valid_sites, T_obs),
  value = as.vector(t(N_mat)), 
  Metric = "Estimated Latent Abundance (N)"
)

# 3. Extract Aggregated Observed Data (Y)
df_obs <- data.frame(
  year = rep(obs_years, each = N_sites),
  grid_id = rep(valid_sites, T_obs),
  value = obs_agg$Total_Y, # <-- FIXED: Now exactly 45,875 long!
  Metric = "Observed Detections (Y)"
)

# 4. Combine and Cap for Outliers (No Relative Normalization)
cap_val <- max(df_est$value, na.rm = TRUE) 

df_combined <- bind_rows(df_est, df_obs) %>%
  left_join(coords, by = "grid_id") %>%
  mutate(
    value_cap = pmin(value, cap_val)
  )

# 5. Build the Plot 
p_side_by_side <- ggplot(df_combined, aes(x = X, y = Y, fill = value_cap)) +
  geom_tile() +
  scale_fill_viridis_c(option = "magma", 
                       labels = scales::comma, 
                       name = "Absolute\nCount", 
                       trans="log1p",
                       breaks = c(0, 1, 5, 10, 50, 100, 200, 300)) +
  facet_wrap(~ Metric, ncol = 2) +
  theme_void() +
  theme(
    strip.text = element_text(size = 14, face = "bold", margin = margin(b = 10, t = 10)),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(b = 15)),
    legend.position = "bottom",
    legend.key.width = unit(2, "cm")
  ) +
  labs(title = "Year: {closest_state}") +
  transition_states(year, transition_length = 2, state_length = 1)

# 6. Render the GIF 
gganimate::animate(
  p_side_by_side, 
  nframes = length(obs_years) * 4, 
  fps = 4, 
  width = 800, 
  height = 450, 
  renderer = gifski_renderer("observed_vs_estimated_absolute.gif")
)

cat("GIF saved as 'observed_vs_estimated_absolute.gif'\n")

# ==============================================================================
# 13c. TRIPLE GIF: OBSERVED vs. NAIVE RATE vs. ESTIMATED N (ABSOLUTE SCALE) ####
# ==============================================================================
cat("\n--- 13c. Generating Triple Comparison Animation (Absolute) ---\n")

library(ggnewscale)

# 1. Ensure Coordinates are synced
my_grid_final <- my_grid_mainland %>%
  filter(grid_id %in% valid_sites) %>%
  arrange(match(grid_id, valid_sites))

coords <- as.data.frame(st_coordinates(st_centroid(my_grid_final)))
coords$grid_id <- valid_sites

# 2. Extract Data Layers
df_est <- data.frame(
  year = rep(obs_years, each = N_sites),
  grid_id = rep(valid_sites, T_obs),
  value = as.vector(t(N_mat)),
  Metric = "3. Model Estimated Abundance (N)"
)

df_obs <- data.frame(
  year = rep(obs_years, each = N_sites),
  grid_id = rep(valid_sites, T_obs),
  value = obs_agg$Total_Y, # <-- FIXED
  Metric = "1. Raw Observed Detections (Y)"
)

df_naive <- data.frame(
  year = rep(obs_years, each = N_sites),
  grid_id = rep(valid_sites, T_obs),
  # <-- FIXED: Using aggregated Y and K
  value = obs_agg$Total_Y / pmax(obs_agg$Total_K, 1), 
  Metric = "2. Naive Detection Rate (Y/Effort)"
)

# 3. Combine and Cap
df_triple <- bind_rows(df_est, df_obs, df_naive) %>%
  left_join(coords, by = "grid_id") %>%
  mutate(
    value_cap = pmin(value, cap_val))

# 4. Build the Triple Facet Plot with INDEPENDENT SCALES
p_triple <- ggplot() +
  
  # --- LAYER 1: ABSOLUTE COUNTS (Y and N) ---
  geom_tile(
    data = df_triple %>% filter(Metric != "2. Naive Detection Rate (Y/Effort)"),
    aes(x = X, y = Y, fill = value_cap)
  ) +
  scale_fill_viridis_c(
    option = "magma", 
    labels = scales::comma, 
    name = "Absolute\nCount", 
    trans = "log1p",
    breaks = c(0, 1, 5, 20, 100, 300)
  ) +
  
  # --- THE MAGIC TRIGGER ---
  new_scale_fill() +
  
  # --- LAYER 2: NAIVE DETECTION RATE (0 to 1) ---
  geom_tile(
    data = df_triple %>% filter(Metric == "2. Naive Detection Rate (Y/Effort)"),
    aes(x = X, y = Y, fill = value) 
  ) +
  scale_fill_viridis_c(
    option = "mako", 
    labels = scales::percent, 
    name = "Naive Rate\n(Y / Visits)"
  ) +
  
  # --- FACET AND ANIMATE ---
  facet_wrap(~ Metric, ncol = 3) +
  theme_void() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 10)),
    legend.position = "bottom",
    legend.box = "horizontal" 
  ) +
  labs(title = "Species Expansion Dynamics | Year: {closest_state}") +
  transition_states(year, transition_length = 2, state_length = 1)

# 5. Render
cat("Rendering Triple GIF... this might take a minute...\n")
animate(
  p_triple, 
  nframes = length(obs_years) * 4,
  fps = 5,
  width = 1200, 
  height = 450,
  renderer = gifski_renderer("triple_comparison_expansion_absolute.gif")
)

cat("Triple GIF saved as 'triple_comparison_expansion_absolute.gif'\n")


# ==============================================================================
# 14. RANGE EXPANSION MAP
# ==============================================================================
cat("\n--- 14. Generating Latent Range Expansion Map ---\n")

# Summary table for change calculation
change_df <- data.table(
  grid_id = valid_sites,
  N_Start = N_mat[1, ],
  N_End = N_mat[T_obs, ],
  N_Change = N_mat[T_obs, ] - N_mat[1, ]
)

# Join back to spatial grid
grid_expansion <- my_grid_mainland %>%
  filter(grid_id %in% valid_sites) %>%
  left_join(change_df, by = "grid_id")

ggplot(grid_expansion) +
  geom_sf(aes(fill = N_Change), color = NA) +
  scale_fill_gradient2(low = "#d7191c", mid = "grey95", high = "#2c7bb6", 
                       midpoint = 0, name = "Delta N") +
  theme_minimal() +
  labs(title = "Latent Population Shift",
       subtitle = "Blue = Colonization/Increase | Red = Decline")


library(scales) # Make sure scales is loaded!

# ==============================================================================
# 14. RANGE EXPANSION MAP
# ==============================================================================
cat("\n--- 14. Generating Latent Range Expansion Map ---\n")

# Summary table for change calculation
change_df <- data.table(
  grid_id = valid_sites,
  N_Start = N_mat[1, ],
  N_End = N_mat[T_obs, ],
  N_Change = N_mat[T_obs, ] - N_mat[1, ]
)

# Join back to spatial grid
grid_expansion <- my_grid_mainland %>%
  filter(grid_id %in% valid_sites) %>%
  left_join(change_df, by = "grid_id")

ggplot(grid_expansion) +
  geom_sf(aes(fill = N_Change), color = NA) +
  scale_fill_gradient2(
    low = "#d7191c", mid = "grey95", high = "#2c7bb6", 
    midpoint = 0, 
    limits = c(-50, 50),            # Set the visual cap!
    oob = scales::squish,           # Squish the extreme outliers
    name = "Delta N\n(2000 to 2024)"
  ) +
  theme_minimal() +
  labs(
    title = "Historical Latent Population Shift (2000 - 2024)",
    subtitle = "Blue = Colonization/Increase | Red = Decline"
  )

# ==============================================================================
# 15. DISPERSAL KERNEL
# ==============================================================================
library(ggplot2)
library(dplyr)
library(posterior)

# 1. Extract alpha draws (ensure you're using the final_draws object)
alpha_samples <- as_draws_df(subset_draws(final_draws, variable = "alpha"))$alpha

# 2. Define the distance range to plot (0 to 200km, matching your truncation)
dist_seq <- seq(0, 200, length.out = 100)

# 3. Calculate kernel density for every posterior sample across the distances
# Result is a matrix: rows = samples, cols = distances
kernel_matrix <- sapply(dist_seq, function(x) {
  (1 / (2 * alpha_samples)) * exp(-x / alpha_samples)
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


# ============================================================================== #
# 17. TREND LINES ####
# ============================================================================== #
library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)

cat("\n--- Extracting Country-Level Abundance Trends (Custom GPKG) ---\n")

# 1. Load YOUR custom countries file 
# (Assuming it is saved in your main Data folder based on your previous paths)
countries_path <- "C:/Users/lieve/Documents/School/2025 - 2026/Thesis Statistics/Data/countries.gpkg"
europe <- st_read(countries_path, quiet = TRUE)

# 2. Extract Centroids of your valid grid cells
centroids <- st_centroid(my_grid_final)

# 3. Match CRS (Just to be absolutely safe!)
europe <- st_transform(europe, st_crs(centroids))

# 4. The Foolproof Spatial Join (Nearest Feature)
# This finds the closest country for EVERY single cell, eliminating ocean gaps.
nearest_idx <- st_nearest_feature(centroids, europe)

# IMPORTANT: I am assuming your country name column is called "name" or "NAME". 
# If you run `names(europe)` and it is called something else (like "CNTR_NAME"), 
# change `europe$name` below to match it!
grid_country_map <- data.frame(
  grid_id = my_grid_final$grid_id,
  Country = europe$name[nearest_idx] # <-- Adjust "name" if your column is different
)

# 5. Extract Abundance (N_mat) into a Long Dataframe
df_n_long <- data.frame(
  year = rep(obs_years, each = N_sites),
  grid_id = rep(valid_sites, T_obs),
  abundance = as.vector(t(N_mat)) 
)

# 6. Merge Abundance with Country Data
df_n_merged <- df_n_long %>%
  left_join(grid_country_map, by = "grid_id")

# ---------------------------------------------------------#
# 7. Aggregate: Mean abundance per cell (with biological threshold)
country_trends <- df_n_merged %>%
  mutate(clean_abundance = ifelse(abundance < 0.1, 0, abundance)) %>%
  group_by(year, Country) %>%
  summarize(Mean_N = mean(clean_abundance), .groups = "drop") %>%
  group_by(Country) %>%
  filter(max(Mean_N) > 1) %>% 
  ungroup()

# 7b. REORDER THE LEGEND (Based on the final year's abundance)
# Find the exact order of countries in the last year (2024)
country_order <- country_trends %>%
  filter(year == max(year)) %>%         # Look only at the end of the timeline
  arrange(desc(Mean_N)) %>%             # Sort Highest to Lowest
  pull(Country)                         # Extract the names in that order

# Apply that strict order to the dataframe as a factor
country_trends <- country_trends %>%
  mutate(Country = factor(Country, levels = country_order))

# 8. Aggregate for Europe
overall_trend <- df_n_merged %>%
  mutate(clean_abundance = ifelse(abundance < 0.1, 0, abundance)) %>%
  group_by(year) %>%
  summarize(Mean_N = mean(clean_abundance), .groups = "drop") %>%
  mutate(Country = "Total Europe (Average)")

# 9. Build the Plot
p_country_trends <- ggplot() +
  geom_line(data = country_trends, aes(x = year, y = Mean_N, color = Country), linewidth = 0.8) +
  geom_line(data = overall_trend, aes(x = year, y = Mean_N), linewidth = 1.2, color = "black", linetype = "dashed") +
  scale_color_viridis_d(option = "turbo") +
  theme_minimal() +
  labs(
    title = "Latent Abundance Over Time by Country",
    subtitle = "Average N per grid cell (Legend ordered by final abundance)",
    x = "Year", y = "Average Estimated Abundance (N)", color = "Country"
  ) +
  theme(
    legend.position = "right", 
    plot.title = element_text(size = 16, face = "bold")
  )

print(p_country_trends)


# ============================================================================== #
## 17b. CALCULATING OVERALL EUROPEAN TREND WITH 95% CI ####
# ============================================================================== #
library(ggplot2)
library(dplyr)
library(posterior)

cat("\n--- Calculating Overall European Trend with 95% CI ---\n")
cat("\n--- Re-extracting Latent Abundance Matrix ---\n")

# FIXED: Changed from "N_obs_vec" to "N_state" to match your updated Stan model!
N_draws_obj <- subset_draws(final_draws, variable = "N_state")

# Convert it to a usable matrix
N_mat_clean <- as_draws_matrix(N_draws_obj)

cat("Success! Matrix reconstructed with", nrow(N_mat_clean), "draws.\n")

# 1. Get the total number of posterior draws
n_draws <- nrow(N_mat_clean)

# 2. Create an empty matrix to hold the European average N for each draw, each year
# Rows = MCMC Draws, Columns = Years
Eu_N_draws <- matrix(NA, nrow = n_draws, ncol = T_obs)

# 3. Loop through each year, apply the threshold, and calculate the mean N
cat("Extracting posterior draws for each year...\n")
for(t in 1:T_obs) {
  
  # Find the exact columns in the Stan output for Year t
  # (Because the vector is ordered: all sites Year 1, all sites Year 2, etc.)
  start_col <- (t - 1) * N_sites + 1
  end_col   <- t * N_sites
  
  # Extract all draws for all 1,835 sites in Year t
  year_data <- N_mat_clean[, start_col:end_col]
  
  # Apply the biological threshold (kill the IDE ghost tails < 0.1)
  year_data_clean <- ifelse(year_data < 0.1, 0, year_data)
  
  # Calculate the average abundance per cell across all of Europe, for every draw
  Eu_N_draws[, t] <- rowMeans(year_data_clean)
}

# 4. Summarize the posterior distribution into Median and 95% CI bounds
df_overall_ci <- data.frame(
  Year = obs_years,
  Median = apply(Eu_N_draws, 2, median),
  Lower  = apply(Eu_N_draws, 2, quantile, probs = 0.025),
  Upper  = apply(Eu_N_draws, 2, quantile, probs = 0.975)
)

# 5. Build the Plot!
p_overall_ci <- ggplot(df_overall_ci, aes(x = Year)) +
  
  # A. Add the 95% Credible Interval ribbon
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "dodgerblue4", alpha = 0.2) +
  
  # B. Add the median trend line
  geom_line(aes(y = Median), color = "dodgerblue4", linewidth = 1.5) +
  
  # C. Add points for each year for clarity
  geom_point(aes(y = Median), color = "black", size = 2) +
  
  theme_minimal() +
  labs(
    title = "Overall European Expansion Trend (2000-2024)",
    subtitle = "Average Abundance (N) per grid cell with 95% Credible Interval",
    x = "Year",
    y = "Average Estimated Abundance (N)"
  ) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    panel.grid.minor = element_blank()
  )

print(p_overall_ci)

# ============================================================================== #
## 17c. TOTAL ABUNDANCE OVER TIME BY COUNTRY ####
# ============================================================================== #
cat("\n--- Extracting Total Country-Level Abundance Trends ---\n")

# 1. Aggregate: Total abundance per country (with biological threshold)
# We use df_n_merged which is already in your environment
country_totals <- df_n_merged %>%
  mutate(clean_abundance = ifelse(abundance < 0.1, 0, abundance)) %>%
  group_by(year, Country) %>%
  summarize(Total_N = sum(clean_abundance), .groups = "drop") %>%
  # Filter out countries with virtually no population to keep the plot clean
  group_by(Country) %>%
  filter(max(Total_N) > 50) %>% 
  ungroup()

# 2. REORDER THE LEGEND (Based on the final year's total abundance)
country_order_total <- country_totals %>%
  filter(year == max(year)) %>%         # Look only at the end of the timeline
  arrange(desc(Total_N)) %>%            # Sort Highest to Lowest
  pull(Country)                         # Extract the names in that order

# Apply that strict order to the dataframe as a factor
country_totals <- country_totals %>%
  mutate(Country = factor(Country, levels = country_order_total))

# Note: We omit the "Total Europe" line here because the sum of all Europe 
# would be so massive it would squish all the country lines to the bottom of the plot!

# 3. Build the Plot
p_country_totals <- ggplot() +
  geom_line(data = country_totals, aes(x = year, y = Total_N, color = Country), linewidth = 1) +
  scale_color_viridis_d(option = "turbo") +
  theme_minimal() +
  # Use commas for the massive numbers on the y-axis
  scale_y_continuous(labels = scales::comma) + 
  labs(
    title = "Total Latent Abundance Over Time by Country",
    subtitle = "Sum of N across all grid cells (Legend ordered by final abundance)",
    x = "Year", 
    y = "Total Estimated Population (Sum of N)", 
    color = "Country"
  ) +
  theme(
    legend.position = "right", 
    plot.title = element_text(size = 16, face = "bold")
  )

print(p_country_totals)

# ============================================================================== #
# 17d. TOTAL EUROPEAN ABUNDANCE OVER TIME (WITH 95% CI) ####
# ============================================================================== #
library(ggplot2)
library(dplyr)

cat("\n--- Calculating Total European Abundance with 95% CI ---\n")

# Assuming N_mat_clean is still in your environment from the previous block.
n_draws <- nrow(N_mat_clean)

# 1. Create an empty matrix to hold the European TOTAL N for each draw, each year
Eu_Total_draws <- matrix(NA, nrow = n_draws, ncol = T_obs)

# 2. Loop through each year, apply the threshold, and calculate the SUM N
cat("Calculating total posterior draws for each year...\n")
for(t in 1:T_obs) {
  
  # Find the exact columns in the Stan output for Year t
  start_col <- (t - 1) * N_sites + 1
  end_col   <- t * N_sites
  
  # Extract all draws for all 1,835 sites in Year t
  year_data <- N_mat_clean[, start_col:end_col]
  
  # Apply the biological threshold (kill the IDE ghost tails < 0.1)
  year_data_clean <- ifelse(year_data < 0.1, 0, year_data)
  
  # Calculate the TOTAL abundance across all of Europe, for every draw
  Eu_Total_draws[, t] <- rowSums(year_data_clean) # <-- CHANGED TO rowSums()
}

# 3. Summarize the posterior distribution into Median and 95% CI bounds
df_total_ci <- data.frame(
  Year = obs_years,
  Median = apply(Eu_Total_draws, 2, median),
  Lower  = apply(Eu_Total_draws, 2, quantile, probs = 0.025),
  Upper  = apply(Eu_Total_draws, 2, quantile, probs = 0.975)
)

# 4. Build the Plot!
p_total_ci <- ggplot(df_total_ci, aes(x = Year)) +
  
  # A. Add the 95% Credible Interval ribbon (Using a nice green to differentiate from the Average plot)
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "seagreen4", alpha = 0.2) +
  
  # B. Add the median trend line
  geom_line(aes(y = Median), color = "seagreen4", linewidth = 1.5) +
  
  # C. Add points for each year for clarity
  geom_point(aes(y = Median), color = "black", size = 2) +
  
  theme_minimal() +
  scale_y_continuous(labels = scales::comma) + # This formats millions cleanly (e.g., 1,500,000)
  labs(
    title = "Total Study Area Population Trend (2000-2024)",
    subtitle = "Total Abundance (Sum of N) across Europe with 95% Credible Interval",
    x = "Year",
    y = "Total Estimated Population (Sum of N)"
  ) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    panel.grid.minor = element_blank()
  )

print(p_total_ci)



# ============================================================================== #
# 18. FORECASTING (SSP2-4.5) ####
# ============================================================================== #
library(geodata)
library(terra)
library(sf)
library(units)
library(splines)
library(posterior)
library(ggplot2)
library(dplyr)
library(gganimate)
library(data.table)
library(scales) 

cat("\n--- Preparing Historical Temperature Matrix ---\n")

# 1. Cast the 'template_obs' data from LONG to WIDE format
wide_temps <- dcast(template_obs, grid_id ~ year, value.var = "mean_temp")

# 2. CRITICAL STEP: Ensure the rows perfectly match the Stan spatial order
wide_temps <- wide_temps[match(valid_sites, wide_temps$grid_id), ]

# 3. Convert to a pure numeric matrix (dropping the grid_id column)
historical_temp_matrix <- as.matrix(wide_temps[, -1, with = FALSE])

cat("Historical Matrix built! Dimensions:", nrow(historical_temp_matrix), "sites x", ncol(historical_temp_matrix), "years.\n")

# ------------------------------------------------------------------------------
cat("\n--- 1. Downloading Climate Data (Current Baseline & Future Projection) ---\n")

data_path_future <- "C:/Users/lieve/Documents/School/2025 - 2026/Thesis Statistics/Data"
climate_dir <- file.path(data_path_future, "Climate_Rasters")
if (!dir.exists(climate_dir)) dir.create(climate_dir, recursive = TRUE)

future_raster_file <- file.path(climate_dir, "future_temp_bio1_ssp245.tif")
current_raster_file <- file.path(climate_dir, "current_temp_bio1_worldclim.tif")

# A. Download FUTURE climate data (SSP2-4.5: 2041-2060)
if (file.exists(future_raster_file)) {
  cat("--- Future climate raster found locally! Loading instantly... ---\n")
  future_temp_raster <- terra::rast(future_raster_file)
} else {
  cat("--- Local file not found. Downloading CMIP6 climate data... ---\n")
  future_climate_ssp2 <- cmip6_world(model = "MPI-ESM1-2-HR", ssp = "245", time = "2041-2060", var = "bioc", res = 5, path = climate_dir)
  future_temp_raster <- future_climate_ssp2[[1]]
  terra::writeRaster(future_temp_raster, future_raster_file, overwrite = TRUE)
}

# B. Download CURRENT climate data (Baseline for calculating the Delta)
if (file.exists(current_raster_file)) {
  cat("--- Baseline climate raster found locally! Loading instantly... ---\n")
  current_temp_raster <- terra::rast(current_raster_file)
} else {
  cat("--- Local file not found. Downloading WorldClim baseline... ---\n")
  current_climate <- worldclim_global(var = "bio", res = 5, path = climate_dir)
  current_temp_raster <- current_climate[[1]]
  terra::writeRaster(current_temp_raster, current_raster_file, overwrite = TRUE)
}

# C. Extract temperatures for your specific grid cells
grid_centroids <- st_centroid(my_grid_final)
future_temp_vec <- terra::extract(future_temp_raster, vect(grid_centroids))[, 2]
current_temp_vec <- terra::extract(current_temp_raster, vect(grid_centroids))[, 2]

# Fill any NAs (coastal cells) with the mean
future_temp_vec[is.na(future_temp_vec)] <- mean(future_temp_vec, na.rm = TRUE)
current_temp_vec[is.na(current_temp_vec)] <- mean(current_temp_vec, na.rm = TRUE)

# ------------------------------------------------------------------------------
cat("\n--- 2. Building the Smooth Climate Ramp (The Delta Method) ---\n")

actual_2024_temps <- historical_temp_matrix[, 25]
warming_delta <- future_temp_vec - current_temp_vec

T_future <- 20
future_temp_matrix_smooth <- matrix(NA, nrow = N_sites, ncol = T_future)

for (t in 1:T_future) {
  fraction <- t / T_future 
  future_temp_matrix_smooth[, t] <- actual_2024_temps + (fraction * warming_delta)
}

# ------------------------------------------------------------------------------
cat("\n--- 3. Plotting the Overall Temperature Trend (2000 - 2044) ---\n")

hist_means <- colMeans(historical_temp_matrix)
future_means <- colMeans(future_temp_matrix_smooth)

df_temp_trend <- data.frame(
  Year = 2000:2044,
  Mean_Temp = c(hist_means, future_means),
  Era = c(rep("Historical (ERA5)", 25), rep("Projected (CMIP6 Delta)", 20))
)

p_temp_trend <- ggplot(df_temp_trend, aes(x = Year, y = Mean_Temp, color = Era)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  scale_color_manual(values = c("Historical (ERA5)" = "grey40", "Projected (CMIP6 Delta)" = "firebrick")) +
  theme_minimal() +
  labs(
    title = "European Mean Annual Temperature Trend",
    subtitle = "Bridging Historical ERA5 with CMIP6 SSP2-4.5 Projections",
    y = "Mean Annual Temperature (°C)", x = "Year"
  ) +
  theme(legend.position = "bottom")

print(p_temp_trend)

# ------------------------------------------------------------------------------
cat("\n--- 4. Building the Dispersal Kernel ---\n")

dist_m <- st_distance(grid_centroids)
distance_matrix <- matrix(as.numeric(dist_m) / 1000, nrow = N_sites, ncol = N_sites)
alpha_est <- median(as_draws_df(subset_draws(final_draws, variable = "alpha"))$alpha)

W_mat <- exp(-distance_matrix / alpha_est)
diag(W_mat) <- 1 # Self-retention
W_norm <- W_mat / rowSums(W_mat)

# ------------------------------------------------------------------------------
cat("\n--- 5. Running the IDE Forward Simulation ---\n")

r_est <- median(as_draws_df(subset_draws(final_draws, variable = "r"))$r)
phi_est <- apply(as_draws_matrix(subset_draws(final_draws, variable = "phi")), 2, median)

# --- SAFE EXTRACTION OF 2024 STARTING ABUNDANCE ---
start_idx <- (25 - 1) * N_sites + 1
end_idx   <- 25 * N_sites
N_start <- numeric(N_sites)
list_position <- 1

cat("Extracting exact starting abundance for Year 25...\n")
for (i in start_idx:end_idx) {
  param_name <- paste0("N_state\\[", i, "\\]")
  N_start[list_position] <- median(as_draws_df(subset_draws(final_draws, variable = param_name, regex = TRUE))[[1]])
  list_position <- list_position + 1
}

# --- SAFE EXTRACTION OF STATIC NICHE SPLINE ---
# We already calculated niche_mat in Section 11!
thermal_niche_est <- apply(niche_mat, 2, median)

# Define the Simulation Function
simulate_future_SDM <- function(N_init, future_temp_matrix, thermal_niche_vec, phi, r, W_kernel, 
                                T_min, N_thermal_gradient) {
  
  T_future <- ncol(future_temp_matrix)
  N_sites  <- nrow(future_temp_matrix)
  N_future_sim <- matrix(0, nrow = N_sites, ncol = T_future)
  N_current    <- N_init
  
  for (t in 1:T_future) {
    # 1. Convert future temperatures to integer indices
    idx <- round((future_temp_matrix[, t] - T_min) * 10) + 1
    
    # 2. Cap indices to ensure they don't exceed the bounds of the gradient
    idx <- pmax(1, pmin(idx, N_thermal_gradient))
    
    # 3. Look up the thermal capacity
    mu_future <- thermal_niche_vec[idx]
    K_future  <- pmax(exp(mu_future + phi), 0.01)
    
    # 4. Population Dynamics
    R0 <- exp(r)
    growth <- (R0 * N_current) / (1.0 + ((R0 - 1.0) / K_future) * N_current)
    N_next <- as.vector(W_kernel %*% growth)
    
    N_future_sim[, t] <- pmax(N_next, 1e-6)
    N_current         <- N_future_sim[, t]
  }
  
  return(N_future_sim)
}

future_abundance_smooth <- simulate_future_SDM(
  N_init = N_start, 
  future_temp_matrix = future_temp_matrix_smooth, 
  thermal_niche_vec = thermal_niche_est, 
  phi = phi_est, 
  r = r_est, 
  W_kernel = W_norm,
  T_min = T_min, 
  N_thermal_gradient = N_thermal_gradient
)

# ------------------------------------------------------------------------------
cat("\n--- 6. Plotting Future Range Shift Map (2024 vs 2044) ---\n")

future_change <- future_abundance_smooth[, T_future] - N_start
df_future <- data.frame(grid_id = my_grid_final$grid_id, Change = future_change)
map_future <- my_grid_final %>% left_join(df_future, by = "grid_id")

p_future_map <- ggplot(map_future) +
  geom_sf(aes(fill = Change), color = NA) +
  scale_fill_gradient2(
    low = "#d7191c", mid = "grey95", high = "#2c7bb6", 
    midpoint = 0, limits = c(-50, 50), oob = scales::squish, 
    name = "Change in N\n(2024 to 2044)"
  ) +
  theme_minimal() +
  labs(
    title = "Predicted Climate-Driven Range Shift",
    subtitle = "SSP2-4.5 Scenario (Projected 20 years into the future)",
    caption = "Blue = Projected Colonization/Increase | Red = Projected Decline"
  )

print(p_future_map)

# ------------------------------------------------------------------------------
cat("\n--- 7. Generating Future SDM Animation (Smooth Transition) ---\n")

coords <- as.data.frame(st_coordinates(grid_centroids))
coords$grid_id <- my_grid_final$grid_id

all_future_N <- cbind(N_start, future_abundance_smooth)
future_years <- 2024:(2024 + T_future)

df_future_anim <- data.frame(
  year = rep(future_years, each = nrow(all_future_N)),
  grid_id = rep(my_grid_final$grid_id, times = length(future_years)),
  abundance = as.vector(all_future_N)
) %>% left_join(coords, by = "grid_id")

p_future_gif <- ggplot(df_future_anim, aes(x = X, y = Y, fill = abundance)) +
  geom_tile() +
  scale_fill_viridis_c(
    option = "magma", name = "Predicted N", trans = "log1p", 
    breaks = c(0, 1, 5, 20, 100, 300, 1000) 
  ) +
  theme_void() +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(b = 10)),
    plot.subtitle = element_text(size = 14, hjust = 0.5, margin = margin(b = 15)),
    legend.position = "bottom"
  ) +
  transition_states(year, transition_length = 2, state_length = 1) +
  labs(title = "Projected Climate-Driven Expansion (SSP2-4.5)", subtitle = "Year: {closest_state}")

gganimate::animate(
  p_future_gif, nframes = length(future_years) * 4, fps = 6, 
  width = 600, height = 600, renderer = gifski_renderer("future_expansion_SSP245_smooth.gif")
)
cat("Success! Saved as 'future_expansion_SSP245_smooth.gif'\n")

# ------------------------------------------------------------------------------
cat("\n--- 8. Extracting Future Country-Level Trends (2024-2044) ---\n")

# Ensure grid_country_map exists from the historical trend blocks
future_years_seq <- 2024:(2024 + T_future)
df_future_long <- data.frame(
  year = rep(future_years_seq, each = N_sites),
  grid_id = rep(valid_sites, length(future_years_seq)),
  abundance = as.vector(all_future_N) 
)

df_future_merged <- df_future_long %>% left_join(grid_country_map, by = "grid_id")

future_country_trends <- df_future_merged %>%
  mutate(clean_abundance = ifelse(abundance < 0.1, 0, abundance)) %>%
  group_by(year, Country) %>%
  summarize(Mean_N = mean(clean_abundance), .groups = "drop") %>%
  group_by(Country) %>%
  filter(max(Mean_N) > 1) %>% 
  ungroup()

future_overall_trend <- df_future_merged %>%
  mutate(clean_abundance = ifelse(abundance < 0.1, 0, abundance)) %>%
  group_by(year) %>%
  summarize(Mean_N = mean(clean_abundance), .groups = "drop") %>%
  mutate(Country = "Total Europe (Average)")

# ORDER LEGEND
country_order <- future_country_trends %>%
  filter(year == max(year)) %>%         
  arrange(desc(Mean_N)) %>%             
  pull(Country)                         

future_country_trends <- future_country_trends %>%
  mutate(Country = factor(Country, levels = country_order))

p_future_trends <- ggplot() +
  geom_line(data = future_country_trends, aes(x = year, y = Mean_N, color = Country), linewidth = 1) +
  geom_line(data = future_overall_trend, aes(x = year, y = Mean_N), linewidth = 1.5, color = "black", linetype = "dashed") +
  scale_color_viridis_d(option = "turbo") +
  theme_minimal() +
  labs(
    title = "Forecasted Latent Abundance (SSP2-4.5)",
    subtitle = "Projected average N per grid cell (2024 - 2044)",
    x = "Year", y = "Average Estimated Abundance (N)", color = "Country"
  ) +
  theme(legend.position = "right", plot.title = element_text(size = 16, face = "bold"))

print(p_future_trends)




# ==============================================================================
# VISUALIZING THE SPATIAL FIELD (GAUSSIAN PROCESS) ####
# ==============================================================================
cat("\n--- Extracting and Visualizing the Spatial Field (phi) ---\n")
library(ggplot2)
library(sf)
library(posterior)

# 1. Extract the posterior draws for the spatial vector 'phi'
# 'phi' contains the spatial offset for each of the 1,835 sites
phi_draws <- as_draws_matrix(fit_spectral_final$draws("phi"))

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










#############
## Extra random stuf ####
#############

cat("\n--- Checking for Same-Day Multi-Dataset Overlaps ---\n")

# 1. Group by Cell, Year, and Day of Year (doy), and count unique datasets
overlap_summary <- data[, .(
  num_datasets = uniqueN(parentDatasetID_2),
  datasets_involved = paste(unique(parentDatasetID_2), collapse = " & "),
  total_visits_that_day = uniqueN(visit)
), by = .(grid_id, year, doy)]

# 2. Filter for only the days where more than 1 dataset was active
same_day_overlaps <- overlap_summary[num_datasets > 1]

# 3. Report the findings
if(nrow(same_day_overlaps) > 0) {
  cat("Found", nrow(same_day_overlaps), "instances where multiple datasets surveyed the exact same cell on the exact same day!\n\n")
  
  cat("Here is a peek at the top 5 overlapping events:\n")
  print(head(same_day_overlaps[order(-total_visits_that_day)], 5))
  
} else {
  cat("No overlapping datasets found! Every cell/day combination was only visited by a single dataset.\n")
}

cat("\n--- Extracting Datasets Involved in Overlaps ---\n")

# 1. Extract the column of concatenated strings
overlap_strings <- same_day_overlaps$datasets_involved

# 2. Split the strings by " & " and unlist them into one giant vector of individual names
all_overlapping_datasets <- unlist(strsplit(overlap_strings, " & "))

# 3. Get the unique list of dataset names
unique_overlapping_datasets <- unique(all_overlapping_datasets)

cat("There are", length(unique_overlapping_datasets), "unique datasets involved in same-day overlaps.\n\n")

# 4. Bonus: Let's see which datasets are the biggest culprits!
# We can tabulate how many times each dataset was involved in an overlapping event
dataset_counts <- table(all_overlapping_datasets)
dataset_counts_df <- as.data.frame(dataset_counts)
colnames(dataset_counts_df) <- c("Dataset", "Overlap_Occurrences")

# Sort from most overlaps to least
dataset_counts_df <- dataset_counts_df[order(-dataset_counts_df$Overlap_Occurrences), ]

cat("--- Top 15 Datasets by Overlap Frequency ---\n")
print(head(dataset_counts_df, 15))


cat("\n--- Analyzing Pairwise Dataset Overlaps ---\n")

# 1. Extract the list of overlapping datasets per event
split_datasets <- strsplit(same_day_overlaps$datasets_involved, " & ")

# 2. Create a function to generate all unique pairs from a single overlapping event
get_pairs <- function(x) {
  # If there's somehow less than 2 datasets, skip it
  if(length(x) < 2) return(NULL) 
  
  # Sort alphabetically so "FRAT & FRFF" is identical to "FRFF & FRAT"
  x <- sort(unique(x)) 
  
  # Generate all combinations of 2
  combos <- combn(x, 2) 
  
  # Return as a data.table
  data.table(Dataset_A = combos[1, ], Dataset_B = combos[2, ])
}

# 3. Apply the function to all 93,803 events and bind them into one massive table
# (This might take a few seconds!)
cat("Calculating all pairwise combinations...\n")
pair_list <- lapply(split_datasets, get_pairs)
pairwise_dt <- rbindlist(pair_list)

# 4. Count up the frequencies of each specific pair
pair_counts <- pairwise_dt[, .(Overlap_Count = .N), by = .(Dataset_A, Dataset_B)]

# 5. Sort from most frequent to least frequent
setorder(pair_counts, -Overlap_Count)

# 6. Display the top 20 overlapping pairs
cat("\n--- Top 20 Most Frequent Dataset Pairs ---\n")
print(head(pair_counts, 20))


cat("\n--- Visualizing the Dataset Overlap Network ---\n")

# Install required visualization packages if you don't have them:
# install.packages(c("tidygraph", "ggraph"))
library(tidygraph)
library(ggraph)
library(ggplot2)
library(dplyr)

# 1. Start with the 'pair_counts' data generated in the previous step.
# We apply a threshold to remove rare edges and keep the plot clean (max_N > 10).
# You can adjust this threshold later if you want more or less detail.
edge_list_filtered <- pair_counts[Overlap_Count > 10]

cat("Generating a network of", nrow(edge_list_filtered), "integrated dataset pairs.\n")

# 2. Convert the clean edge list into a graph object (Nodes and Edges)
graph_obj <- as_tbl_graph(edge_list_filtered, directed = FALSE)

# 3. Enhance the Graph with community detection and centrality (Node Size)
# This automatically calculates key metrics for the plot.
enhanced_graph <- graph_obj %>%
  # A. Calculate Node Size based on 'Centrality Degree' (number of connections)
  # Large nodes = Aggregator hubs. Small nodes = Specialists.
  activate(nodes) %>%
  mutate(centrality_degree = centrality_degree()) %>%
  
  # B. Calculate 'Louvain Community' to auto-detect clusters
  # This finds the tight national clusters (e.g., the full French group will all be Blue).
  mutate(community = factor(group_louvain()))

# 4. Design the Plot using ggraph
# We use the Fruchterman-Reingold ('fr') layout which is great for pushing clusters apart.
p_overlap_net <- ggraph(enhanced_graph, layout = "fr") +
  
  # A. The Edges (Connections)
  # The width of the line is proportional to the frequency (Overlap_Count)
  geom_edge_link(aes(edge_width = Overlap_Count), color = "dodgerblue4", alpha = 0.2) +
  
  # B. The Nodes (Datasets)
  # The color is set by the auto-detected community cluster. The size is set by centrality.
  geom_node_point(aes(color = community, size = centrality_degree), show.legend = c(color = TRUE, size = FALSE)) +
  
  # C. The Labels (Repel makes sure text doesn't overlap text)
  geom_node_text(aes(label = name), repel = TRUE, fontface = "bold", size = 4) +
  
  # D. Formatting
  scale_edge_width_continuous(range = c(0.1, 4), name = "Frequency") + # Edge thickness
  scale_size_continuous(range = c(5, 12)) + # Minimum/Maximum node size
  scale_color_viridis_d(option = "turbo", name = "National\nCluster") + # Nice coloring
  theme_graph() +
  labs(
    title = "European Dragonfly Database Co-occurrence Network",
    subtitle = paste0("Structure of dataset overlaps | (Filtered for pairs with >10 events)"),
    caption = "Large nodes are aggregators. Thick edges are intensely integrated pairs.\nColors are automatically detected communities (Likely countries)."
  ) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    legend.position = "right"
  )

print(p_overlap_net)




cat("\n--- Finding Isolated Datasets & Full Network Visualization ---\n")

library(tidygraph)
library(ggraph)
library(dplyr)

# 1. Get the master list of ALL unique datasets in your raw data
# (Excluding any NAs or the Unvisited Dummy)
all_datasets <- unique(data$parentDatasetID_2)
all_datasets <- all_datasets[!is.na(all_datasets)]

# 2. Get the list of datasets that DID overlap at least once
# (Using the string split from our earlier overlap check)
overlapping_datasets <- unique(unlist(strsplit(same_day_overlaps$datasets_involved, " & ")))

# 3. Find the isolated ones (Datasets in 'all' but NOT in 'overlapping')
isolated_datasets <- setdiff(all_datasets, overlapping_datasets)

cat("There are", length(isolated_datasets), "datasets that NEVER overlap with any other dataset on the same day:\n")
print(isolated_datasets)

# ------------------------------------------------------------------------------
# 4. Build the Full Graph Data
# ------------------------------------------------------------------------------

# To force ggraph to plot unconnected dots, we MUST explicitly define the nodes
nodes_df <- data.frame(name = all_datasets)

# We use the pair_counts table from earlier, keeping ALL > 0 connections this time
# so we see the true, completely unfiltered network
edges_df <- pair_counts[Overlap_Count > 0]

# Build the graph explicitly linking the full node list to the edge list
full_graph <- tbl_graph(nodes = nodes_df, edges = edges_df, directed = FALSE)

# Enhance the graph
enhanced_full_graph <- full_graph %>%
  activate(nodes) %>%
  mutate(centrality_degree = centrality_degree()) %>%
  mutate(community = factor(group_louvain())) %>%
  # Create a logical flag: Is this node completely isolated?
  mutate(is_isolated = centrality_degree == 0)

# ------------------------------------------------------------------------------
# 5. Plot the Full Network
# ------------------------------------------------------------------------------

p_full_net <- ggraph(enhanced_full_graph, layout = "fr") +
  
  # A. The Edges (Connections)
  geom_edge_link(aes(edge_width = Overlap_Count), color = "dodgerblue4", alpha = 0.2) +
  
  # B. The Nodes (We make isolated nodes a different shape so they stand out!)
  geom_node_point(
    aes(color = community, size = centrality_degree, shape = is_isolated), 
    show.legend = c(color = TRUE, size = FALSE, shape = FALSE)
  ) +
  
  # Shape 16 = Solid Circle (Connected), Shape 17 = Solid Triangle (Isolated)
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17)) + 
  
  # C. The Labels
  geom_node_text(aes(label = name), repel = TRUE, fontface = "bold", size = 3.5) +
  
  # D. Formatting
  scale_edge_width_continuous(range = c(0.1, 4), name = "Overlap Frequency") + 
  scale_size_continuous(range = c(4, 12)) + # Minimum size 4 so isolated nodes are visible
  scale_color_viridis_d(option = "turbo", name = "National\nCluster") + 
  theme_graph() +
  labs(
    title = "Complete European Dragonfly Database Network",
    subtitle = "Including completely isolated datasets (Triangles on the perimeter)",
    caption = "Triangles = Datasets with ZERO same-day overlaps.\nCircles = Datasets that act as cross-calibration anchors for the model."
  ) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    legend.position = "right"
  )

print(p_full_net)

