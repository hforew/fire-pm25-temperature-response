# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
########### Re-grid PM data to pop data resolution ##########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Remove all objects from the environment
rm(list = ls())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Packages #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(ncdf4)
library(here)
library(tidyverse)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Import #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Import annual average PM2.5 data and pop data (previously converted from NetCDF to dataframe)
annual_ave <- read_csv(here("output", "annual_ave_pm25.csv"))
pop <- read_csv(here("output", "pop_df_rev.csv")) #may use pop_df_rev.csv here for the latest version

# pm data dimension = (lon x lat) 288 x 192 = 55296
# pop data dimension = (lon x lat) 720 × 360 = 259200

nrow(pop)
nrow(annual_ave)
colnames(pop)
colnames(annual_ave)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Validate POP and PM2.5 Grid Coordinates: Centers vs. Edges #################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# NetCDF file data was previously extracted and converted to dataframe format. 
# check here if the lat and lon coordinates in netCDF represent grid cell center or edges

pm_2000 <- nc_open(here("input", "CESM_09x125_PM25_2000_Baseline.nc"))
print(pm_2000) # PM25 NetCDF files for other years follow identical format

pop_ssp1 <- nc_open(here("input", "SSP1_for_RCP45_2006-2100_population_density_c160701.nc"))
print(pop_ssp1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Validate PM2.5 Grid Coordinates: Centers vs. Edges
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Purpose: Determine whether lon/lat values in PM2.5 NetCDF represent
#          grid cell centers or grid cell edges. This affects regridding method.
#
# Expected for grid CENTERS at 1.25° resolution:
#   - First lon = -180° + 1.25°/2 = -179.375°
#   - Coordinates offset from edges by half the resolution
#
# Expected for grid EDGES:
#   - First lon = -180° (western boundary of globe)
#   - Last lon would be at or near 180° (eastern boundary)
# 

# Extract lon and lat coordinate vectors from NetCDF
lon_pm <- ncvar_get(pm_2000, "lon")
lat_pm <- ncvar_get(pm_2000, "lat")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Check Longitude Coordinates
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat("Longitude:\n")
cat("First few values:", head(lon_pm), "\n")
# Spacing should be consistent at 1.25° resolution
cat("Spacing:", unique(diff(lon_pm)), "\n")
cat("Expected resolution: 1.25°\n\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Check Latitude Coordinates
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat("Latitude:\n")
cat("First few values:", head(lat_pm), "\n")
# Spacing should be ~0.94° (192 cells covering 180° = 0.9375° average)
cat("Spacing:", unique(diff(lat_pm)), "\n")
cat("Expected resolution: ~0.9° to 1.0°\n\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Diagnostic Test: Grid Centers vs. Grid Edges
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# If coordinates are CENTERS of 1.25° cells starting at -180° edge:
#   First center = -180° + (1.25°/2) = -179.375°
#
# If coordinates are EDGES:
#   First edge = -180° exactly

cat("First lon value:", lon_pm[1], "\n")
cat("Expected if center of first cell starting at -180°:", -180 + 1.25/2, "\n\n")

cat("First lat value:", lat_pm[1], "\n")
cat("Last lat value:", lat_pm[length(lat_pm)], "\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RESULT: PM2.5 coordinates are GRID EDGES
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Evidence:
# 1. lon[1] = -180° (not -179.375°) --> starts at western edge of globe
# 2. lat[1] = -90°, lat[192] = 90° --> runs from pole to pole (edges)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Validate Population Grid Coordinates: Centers vs. Edges 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This file explicitly provides edge coordinates, making validation easier
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Extract lon and lat coordinate vectors
lon_pop <- ncvar_get(pop_ssp1, "lon")
lat_pop <- ncvar_get(pop_ssp1, "lat")

# Extract explicit edge coordinates (if wanting to verify)
edge_west <- ncvar_get(pop_ssp1, "EDGEW")
edge_east <- ncvar_get(pop_ssp1, "EDGEE")
edge_south <- ncvar_get(pop_ssp1, "EDGES")
edge_north <- ncvar_get(pop_ssp1, "EDGEN")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Check Longitude Coordinates
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("Population Longitude:\n")
cat("First few values:", head(lon_pop), "\n")
cat("Spacing:", unique(round(diff(lon_pop), 6)), "\n")
cat("Expected resolution: 0.5°\n\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Check Latitude Coordinates
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("Population Latitude:\n")
cat("First few values:", head(lat_pop), "\n")
cat("Spacing:", unique(round(diff(lat_pop), 6)), "\n")
cat("Expected resolution: 0.5°\n\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Diagnostic Test: Grid Centers vs. Grid Edges
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# For 0.5° resolution with CENTERS starting at -180° edge:
#   First center = -180° + (0.5°/2) = -179.75°
#
# For EDGES:
#   First edge = -180° exactly

cat("First lon value:", lon_pop[1], "\n")
cat("Expected if center of first cell starting at -180°:", -180 + 0.5/2, "\n\n")

cat("First lat value:", lat_pop[1], "\n")
cat("Expected if center of first cell starting at -90°:", -90 + 0.5/2, "\n")
cat("Last lat value:", lat_pop[length(lat_pop)], "\n")
cat("Expected if center of last cell ending at 90°:", 90 - 0.5/2, "\n\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Additional Check: Verify Against Explicit Edge Variables
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("Explicit edge coordinates:\n")
cat("Western edge:", edge_west, "\n")
cat("Eastern edge:", edge_east, "\n")
cat("Southern edge:", edge_south, "\n")
cat("Northern edge:", edge_north, "\n\n")

# Calculate what the first center should be based on western edge
cat("First lon calculated from western edge:", edge_west + 0.5/2, "\n")
cat("Actual first lon:", lon_pop[1], "\n\n")

nc_close(pm_2000)
nc_close(pop_ssp1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RESULT: Population coordinates are GRID CENTERS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Evidence:
# 1. lon[1] = -179.75° (not -180°) --> offset by 0.25° from edge
# 2. lat[1] = -89.75° (not -90°) --> offset by 0.25° from pole
# 3. File explicitly provides separate EDGE variables confirming this
#
# Implication:
# - Population grid uses cell centers 
# - PM2.5 grid uses cell edges 
# - PM2.5 must must be converted to cell centers before regridding


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Regrid PM2.5 from 1.25° to 0.5° Resolution ############ 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Method: Nearest neighbor assignment (uniform assignment within coarse cells)
# Each 1.25° PM2.5 cell maps to multiple 0.5° population cells

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 1: Convert PM2.5 coordinates from edges to centers
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Current PM2.5 coordinates are grid edges (lon starts at -180°, lat at -90°)
# Need to shift to cell centers for proper regridding

annual_ave <- annual_ave %>%
  mutate(
    lon_center = lon + 1.25/2,      # Shift east by half cell width (0.625°)
    lat_center = lat + 0.9424/2     # Shift north by half cell height (0.4712°)
  )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 2: Create distinct population coordinates for matching
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extract unique lon/lat combinations from population grid
# These are the target coordinates for regridding (0.5° resolution)

pop_coords <- pop %>%
  distinct(lon, lat) %>%           # Get unique coordinate pairs (removes duplicates if multiple years)
  arrange(lon, lat)                # Sort west to east, south to north for organized output

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 3: Find nearest PM2.5 cell for each population cell
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# For each 0.5° pop cell, calculate distance to all 1.25° PM2.5 cell centers
# Assign the PM2.5 value from the nearest cell (nearest neighbor interpolation)

cat("Regridding", nrow(pop_coords), "population cells...\n")

# Loop through each pop cell and find its nearest PM2.5 cell
nearest_pm_idx <- map2_int( # map2_int interates over 2 vectors
  pop_coords$lon,                  # Vector of pop cell longitudes
  pop_coords$lat,                  # Vector of pop cell latitudes
  function(plon, plat) {           # Function applied to each pop cell coordinate pair:
    # Calculate Euclidean distance to all PM2.5 cell centers
    dist <- sqrt((plon - annual_ave$lon_center)^2 +    # Distance in lon direction
                   (plat - annual_ave$lat_center)^2)     # Distance in lat direction
    # Return row index of closest PM2.5 cell
    which.min(dist)                # Find index where distance is minimum
  }
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 4: Create regridded PM2.5 dataset at 0.5° resolution
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Match each pop cell to its nearest PM2.5 cell and extract PM2.5 values

# Add nearest PM2.5 index to pop coordinates
pop_coords <- pop_coords %>%
  mutate(pm_idx = nearest_pm_idx)  # Store which PM2.5 row matches each pop cell

# Extract PM2.5 values for matched cells (all pm_ and fpm_ columns)
pm25_regridded <- pop_coords %>%
  left_join(
    annual_ave %>% 
      mutate(pm_idx = row_number()) %>%              # Create index column for PM2.5 rows
      select(pm_idx, starts_with("pm_"),            # Keep pm_2000, pm_2050_45, etc.
             starts_with("fpm_"), n_months),        # Keep fpm_2000, fpm_2050_45, etc.
    by = "pm_idx"                                    # Join based on nearest neighbor index
  ) %>%
  select(-pm_idx)                                    # Remove temporary index column

# Verify regridding succeeded
cat("\nRegridded PM2.5 dimensions:", nrow(pm25_regridded), "rows\n")
cat("Original pop dimensions:", nrow(pop_coords), "rows\n")
cat("Should match:", nrow(pm25_regridded) == nrow(pop_coords), "\n\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 5: Join regridded PM2.5 with population data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Combine regridded PM2.5 (0.5° resolution) with population data
# This creates a single dataset with both PM2.5 and population for each cell

pop_pm_combined <- pop %>%
  left_join(pm25_regridded,                          # Add PM2.5 columns to pop data
            by = c("lon", "lat"))                    # Match on exact coordinates

# Verify join succeeded (check if any PM2.5 values are missing)
cat("Combined dataset dimensions:", nrow(pop_pm_combined), "rows\n")
cat("All pop rows matched:", sum(is.na(pop_pm_combined$pm_2000)) == 0, "\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Save outputs
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

write_csv(pm25_regridded, here("output", "pm25_regridded_0.5deg.csv"))
write_csv(pop_pm_combined, here("output", "pop_pm_combined.csv"))

cat("\nRegridding complete!\n")
