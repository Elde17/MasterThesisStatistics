# ============================================================================== #
# 0. SETUP & LIBRARIES ####
# ============================================================================== #
rm(list = ls())
gc()

## This script is for Sympetrum danae
# To adapt to another species: ctrl+f on Sympetrum danae and replace by target species.

library(readr)
library(sf)
library(data.table)
library(splines)
library(cmdstanr)

# Ensure cmdstan is ready
check_cmdstan_toolchain(fix = TRUE)

setwd("~/School/2025 - 2026/Thesis Statistics/Analyses/FINAL SCRIPTS")

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

# Define target flight window in DOY
start_doy <- 0 
end_doy   <- 366
target_species <- "Sympetrum danae"

# 1. Explicitly define all non-species columns from the dataset
meta_cols <- c("visit", "year", "year_scaled", "doy_scaled", "doy", 
               "ll", "ll_raw", "spread_scaled", "grid_id",
               "parentDatasetID_2", "datasetID")

# 2. Dynamically grab all other columns as species
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

# Aggregation
# This ensures that both Visits (K) and Detections (Y) are counted by unique days.
# This guarantees that Y <= K, which keeps the Binomial distribution happy.
agg_dt <- data[, .(
  # 1. Total unique days the site was surveyed
  num_visits = uniqueN(doy), 
  
  # 2. Total unique days the species was detected at least once
  `Sympetrum danae` = uniqueN(doy[get(target_species) > 0])
  
), by = .(grid_id, year)]

# SAFETY CHECK
# This must return 0. 
math_violation <- sum(agg_dt[[target_species]] > agg_dt$num_visits)
cat("Number of records violating Y <= K logic:", math_violation, "\n")

if(math_violation > 0) stop("Data alignment error: Detections exceed visits!")

cat("Aggregation complete. Unique grid-years available:", nrow(agg_dt), "\n")

# ============================================================================== #
## 1b. EXPLORATIONS: VISITS PER YEAR & GRID CELL ####
# ============================================================================== #
cat("\n--- 1b. Running Exploratory Analyses on Sampling Effort ---\n")

library(ggplot2)
library(viridis)

# 1. Temporal Exploration: Visits per Year
visits_per_year <- agg_dt[, .(
  total_visits = sum(num_visits),
  mean_visits_per_cell = mean(num_visits),
  active_cells = .N
), by = year][order(year)]

cat("\nSummary of Visits per Year:\n")
print(visits_per_year)

# Plot: Distribution of visits over the years
p_years <- ggplot(visits_per_year, aes(x = year, y = total_visits)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(color = "darkblue", size = 2) +
  theme_minimal() +
  labs(title = "Total Sampling Visits Per Year",
       x = "Year", y = "Total Number of Visits")
print(p_years)


# 2. Spatial Exploration: Load grid geography and merge data
cat("\nLoading spatial grid for map generation...\n")
grid_sf <- st_read(grid_path, quiet = TRUE)

# Calculate total visits per grid cell across ALL years
visits_per_cell <- agg_dt[, .(
  total_visits = sum(num_visits),
  years_active = uniqueN(year)
), by = grid_id]

# Merge data table back into the spatial sf object
# (Ensuring the merge keeps the sf structure intact)
grid_mapped <- merge(grid_sf, visits_per_cell, by = "grid_id", all.x = FALSE)

# Plot: Map of total sampling effort
p_map <- ggplot(data = grid_mapped) +
  geom_sf(aes(fill = total_visits), color = NA) +
  scale_fill_viridis_c(option = "viridis", trans = "log10", 
                       name = "Total Visits\n(Log10 Scale)") +
  theme_minimal() +
  labs(title = "Spatial Distribution of Sampling Effort",
       subtitle = "Aggregated visits across all years per grid cell")
print(p_map)


# 3. Spatiotemporal Exploration: Variation across both Year and Grid cell
# Spot checking the top 5% most heavily sampled grid-years
quantile_threshold <- quantile(agg_dt$num_visits, 0.95)
cat(paste0("\n95th percentile of visits per grid-year: ", quantile_threshold, "\n"))

# Distribution check of visits per cell-year combo
p_dist <- ggplot(agg_dt, aes(x = num_visits)) +
  # Using geom_bar because num_visits is discrete (integers)
  geom_bar(fill = "seagreen", color = "black", alpha = 0.7, width = 0.05) + 
  scale_x_log10(breaks = c(1, 2, 5, 10, 20, 50, 100, 300)) +
  theme_minimal() +
  labs(title = "Distribution of Visits per Grid-Year Combination",
       x = "Number of Visits (Log Scale)", 
       y = "Frequency of Grid-Years")

print(p_dist)


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

# ============================================================================== #
# 3. DYNAMIC ERA5 CLIMATE INTEGRATION & MAINLAND CONNECTIVITY ####
# ============================================================================== #
cat("\n--- 3. Integrating Full-Year ERA5 Temperatures & Mainland Connectivity ---\n")

library(sf)
library(data.table)
library(igraph)

# 1. Load the grid
my_grid_all <- st_read(grid_path, quiet = TRUE)
my_grid_all$grid_id <- as.character(my_grid_all$grid_id)

# 2. Define the mainland
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
my_grid_mainland <- my_grid_mainland[order(as.character(my_grid_mainland$grid_id)), ]

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

# ============================================================================== #
## 3b. DYNAMIC ERA5 CLIMATE - visualization ####
# ============================================================================== #
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


# ============================================================================== #
# 4. BUILD TEMPLATE (Synchronized 1995-2024 window) ####
# ============================================================================== #
cat("\n--- 4. Building Spatio-Temporal Template (1,835 Cells) ---\n")

# Ensure all subsequent blocks use the 1,835 mainland IDs
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

# ------------------------------------------------------------------ #
## 4b. VISUALIZATION CHECK ####
# ------------------------------------------------------------------ #
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

# ------------------------------------------------------------------ #
## 4c. VISUALIZE ANCHORED HINDCAST WITH ITS OWN REGRESSION ####
# ------------------------------------------------------------------ #
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

# ============================================================================== #
# 5. TEMPERATURE SPLINES (Monotonic Indexing) ####
# ============================================================================== #
# NOTE: These spline bases were created to model the thermal niche using spline-based
# methods. However, in the most recent version, a parametric double-logistic was 
# used. Thus, these spline bases are no longer used by the model. This part of the 
# code is left in to be able to move back to splines if wanted later on, since it
# does not interfere with the model anyway.

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

# ============================================================================== #
## 5b. VISUALIZE THE B-SPLINE BASIS FUNCTIONS ####
# ============================================================================== #
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


# ============================================================================== #
# 6. OBSERVATION VECTORS & INDEX MAPPING ####
# ============================================================================== #
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
# 7. SPATIAL BASES (20 Habitat Splines + CSR Dispersal) ####
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
# 8. ASSEMBLE STAN DATA LIST (Final Production Version) ####
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

# ============================================================================== #
# 9. FINAL AUDIT ####
# ============================================================================== #
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

