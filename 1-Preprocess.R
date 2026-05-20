setwd("~/School/2025 - 2026/Thesis Statistics/Analyses/FINAL SCRIPTS")

# NOTE: this preprocessing script is largely based on the preprocessing script by Lisa Nicvert

library(terra)
library(funique)
library(readr)
library(dplyr)
library(data.table)
library(sf)
library(ggplot2)

## --- PATHS ---
data_path <- "C:/Users/lieve/Documents/School/2025 - 2026/Thesis Statistics/Data/New/Raw/occ_all.rds"
taxa_path <- "C:/Users/lieve/Documents/School/2025 - 2026/Thesis Statistics/Data/New/Raw/taxo.rds"
grid_path <- "C:/Users/lieve/Documents/School/2025 - 2026/Thesis Statistics/Data/New/grid.gpkg"
countries_path <- "C:/Users/lieve/Documents/School/2025 - 2026/Thesis Statistics/Data/New/countries.gpkg"

## --- DATA PREPARATION ---
raw <- read_rds(data_path)
str(raw)

# View the unique combinations of Database ID and Name
unique(raw[, .(parentDatasetID, parentDatasetName)])[order(parentDatasetID)]

# Create the consolidated ID (eventually not used in the final model)
raw$parentDatasetID_2 <- ifelse(is.na(raw$parentDatasetID), 
                                as.character(raw$datasetID), 
                                as.character(raw$parentDatasetID))

# Clean LifeStage
raw$lifeStage <- ifelse(is.na(raw$lifeStage), "unknown", as.character(raw$lifeStage))

# Filter for adults & unknown (likely adults)
raw_2 <- raw %>% filter(lifeStage %in% c("adult", "exuvia", "teneral", "unknown"))
setDT(raw_2)

## --- TAXONOMY JOIN ---
taxa <- readRDS(taxa_path)
setDT(taxa)
taxa <- taxa[, .(taxonID, species, genus, family, gbifID)]

# Join - ensure we keep the ID and lifeStage columns from raw_2
occ <- taxa[raw_2, on = "taxonID"]

## --- SPATIAL FILTERS ---
occ <- occ[!is.na(eventDate), ]
occ <- occ[eventDate < as.IDate("2025-01-01") & eventDate > as.IDate("2000-01-01"), ]
occ <- occ[occurrenceStatus == "present", ]
occ[, doy := yday(eventDate)]
occ <- occ[doy != 1, ] # Remove Jan 1st noise
occ <- occ[taxonRank %in% c("SPECIES", "SUBSPECIES"), ]

## --- SPATIAL GRIDDING ---
coord <- funique(occ[, .(locationID, decimalLongitude, decimalLatitude)])

# Create SpatVector (raw coords are EPSG:3035 based)
vcoord <- vect(coord, geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:3035")

grid <- st_read(grid_path, quiet = TRUE)
names(grid)[names(grid) == "GRD_ID"] <- "grid_id"
grid_v <- vect(grid)

# Extract Grid IDs
id_grid <- terra::extract(grid_v, vcoord)
coord$grid_id <- id_grid$grid_id

# Map grid_id back to main table
occ <- merge(occ, coord[, .(locationID, grid_id)], by = "locationID", all.x = TRUE)
occ <- occ[!is.na(grid_id), ]

## --- COVARIATES ---
occ[, year := year(eventDate)]
year_scaled <- scale(occ$year)
occ[, year_scaled := as.numeric(year_scaled)]
doy_scaled <- scale(occ$doy)
occ[, doy_scaled := as.numeric(doy_scaled)]

ylong <- occ

# --- DEFINE VISIT BY METADATA ---
# This ensures that one "visit" row contains exactly one set of these IDs
ylong[, visit := .GRP, by = .(eventDate, grid_id, parentDatasetID_2, datasetID)]
ylong[, occ_score := 1]

## --- LIST LENGTH (per visit) ---
# calculated on a temporary wide table first
ywide_temp <- dcast(ylong, visit ~ species, value.var = "occ_score", fun.aggregate = length)
spcol <- grep(pattern = "[A-Z][a-z]+ [a-z]+", colnames(ywide_temp), value = TRUE)
ywide_temp[, ll_raw := rowSums(.SD > 0), .SDcols = spcol]
ywide_temp[, ll := fcase(ll_raw == 1, "1", ll_raw <= 3, "2-3", default = ">3")]

# Merge List Length back to ylong
ylong <- merge(ylong, ywide_temp[, .(visit, ll, ll_raw)], by = "visit")

## --- SPATIAL SPREAD ---
mean_dist <- function(lat, lon) {
  if (length(lat) <= 1) return(0)
  mx <- mean(lon); my <- mean(lat)
  sqrt((lon - mx)^2 + (lat - my)^2) |> mean()
}
ylong <- ylong[, spread_m := mean_dist(decimalLatitude, decimalLongitude), by = visit]
spread_scaled <- scale(ylong$spread_m)
ylong[, spread_scaled := as.numeric(spread_scaled)]

## --- FINAL GENERATE YWIDE ---
# 'species' must be on the Right Hand Side (RHS) of the ~ to become the columns
ywide <- dcast(ylong,
               visit + year + year_scaled + doy + doy_scaled + ll + 
                 spread_scaled + grid_id + parentDatasetID_2 + datasetID ~ species, 
               value.var = "occ_score",
               fun.aggregate = function(x) as.numeric(sum(x) > 0))

## --- VERIFICATION ---
ncol(ywide)
print(paste("Final ywide rows:", nrow(ywide)))
# Check that the columns exist
print(head(ywide[, .(visit, grid_id, parentDatasetID_2, datasetID, ll)]))

# Get sums of all species columns
species_counts <- colSums(ywide[, ..spcol])
# See if any species have 0 detections
sum(species_counts == 0)

# Count total detections for each species column
sp_counts <- colSums(ywide[, ..spcol])

# Identify species with zero detections
zero_sp <- names(sp_counts[sp_counts == 0])
print(paste("Species with zero detections to remove:", length(zero_sp)))

# Identify 'ultra-rare' species (e.g., < 10 detections)
# Most occupancy models will fail to converge with fewer than 10-30 detections
rare_sp <- names(sp_counts[sp_counts > 0 & sp_counts < 10])

# Plot number of visits per year, colored by the top 5 datasets
top_datasets <- ywide[, .N, by = parentDatasetID_2][order(-N)][1:5, parentDatasetID_2]

ggplot(ywide[parentDatasetID_2 %in% top_datasets], aes(x = year, fill = parentDatasetID_2)) +
  geom_bar() +
  theme_minimal() +
  labs(title = "Contribution of Top Datasets over Time")

# Check for artificial peaks on the 1st of every month
ggplot(ylong, aes(x = doy)) + 
  geom_histogram(bins = 365) + 
  geom_vline(xintercept = c(32, 60, 91), color = "red", linetype = "dashed") +
  labs(title = "Check for artificial date peaks (Red lines = 1st of Month)")

# Save result
saveRDS(ywide, "C:/Users/lieve/Documents/School/2025 - 2026/Thesis Statistics/Data/New/data_processed.rds")
