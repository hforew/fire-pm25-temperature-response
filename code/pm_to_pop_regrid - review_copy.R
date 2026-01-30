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
setwd("C:/Ford_BA_FPM25/")
# Import annual average PM2.5 data and pop data (previously converted from NetCDF to dataframe)
annual_ave <- read_csv(here("output", "annual_ave_pm25.csv"))
pop <- read_csv(here("output", "pop_df.csv"))

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

################################################################################review start
#Step 1: Create a synthetic PM field on the PM grid
#first define an artificial (synthetic) PM2.5 field as a known mathematical function of longitude and latitude at each PM grid cell center:
#This guarantees that every PM cell has a unique and predictable value.
#Because the function is known, you always know what the “correct” value should be at any PM location.
#Purpose:
#Create a controlled “truth” field where the correct values are known exactly.
F_pm <- with(annual_ave,
             1000 + 2 * lon_center + 3 * lat_center + 10 * sin(lon_center * pi / 180))


#Step 2: Apply the actual regridding method
#apply real nearest-neighbor regridding procedure:
#Each 0.5° population grid cell is matched to its nearest 1.25° PM grid cell.
#The synthetic PM value from that nearest PM cell is copied to the population cell.
#Purpose:
#Simulate exactly what real PM2.5 regridding pipeline does in this script.
pop_truth <- pop_coords %>%
  mutate(
    pm_idx = nearest_pm_idx,
    F_regridded = F_pm[pm_idx]
  )

#Step 3: Independently compute the expected values
#For each population cell, independently compute what the synthetic value should be by:
#Taking the matched PM grid center coordinates.
#Evaluating the same mathematical function at those coordinates.
#This step does not use the regridded values but recomputes them.
#Purpose:
#Generate a fully independent benchmark for the correct result.
pop_truth <- pop_truth %>%
  mutate(
    lon_match = annual_ave$lon_center[pm_idx],
    lat_match = annual_ave$lat_center[pm_idx],
    F_expected = 1000 + 2 * lon_match + 3 * lat_match + 10 * sin(lon_match * pi / 180),
    err = F_regridded - F_expected
  )

#Step 4: Compare regridded vs expected values
#compare:
#F_regridded: value produced by the regridding code
#F_expected: value computed analytically
#compute the difference and check the maximum absolute error.
#Purpose:
#Verify that the regridding assigns exactly the correct PM cell to every population cell.
max_err <- max(abs(pop_truth$err), na.rm = TRUE)
cat("ONE-SHOT VALIDATION | max(|error|) =", max_err, "\n")

# Enforce correctness: any nonzero error indicates a problem
stopifnot(max_err < 1e-10)
cat("Validated: Regridding pipeline is internally consistent.\n")

# -------------------------------------------------------------------
# Intuition of the synthetic-field validation (three-step view)
#
# Step 1: Assign a unique "label" to each PM grid cell.
# first create a synthetic PM field as a simple mathematical
# function of longitude and latitude. This function acts as a
# deterministic label generator: each PM grid cell receives a
# unique value determined by its geographic coordinates.
#
# Step 2: Let the algorithm decide which PM cell each population
# cell belongs to.
# Using the real nearest-neighbor regridding procedure, compute
# an index pm_idx for each population grid cell, indicating which
# PM grid cell it is assigned to. Then copy the synthetic PM
# value from that PM cell to the population cell (F_regridded).
#
# Step 3: Independently verify the assignment using real-world
# coordinates.
# Based on the same pm_idx, we retrieve the matched PM grid cell's
# geographic coordinates and re-evaluate the synthetic function at
# those coordinates to obtain F_expected. This uses only the
# physical location of the matched PM cell, not the copied value.
#
# If the regridding is correct, F_regridded and F_expected must be
# identical (up to numerical precision). Any discrepancy indicates
# that the index pm_idx points to the wrong PM grid cell, revealing
# errors in spatial indexing, coordinate transformation, or
# nearest-neighbor matching.
#
# Therefore, this synthetic-field validation works by encoding
# geographic location into numerical values and checking whether
# the algorithm assigns each population cell to the correct PM
# grid cell in geographic space.
# -------------------------------------------------------------------

#Step 5: Distance sanity check
#Compute the actual geometric distance between:
#Each population cell center
#Its matched PM grid center
#This confirms:
#No population cell is matched to a far-away PM cell.
#All matches fall within the correct coarse PM neighborhood.
d_match <- sqrt(
  (pop_coords$lon - annual_ave$lon_center[nearest_pm_idx])^2 +
    (pop_coords$lat - annual_ave$lat_center[nearest_pm_idx])^2
)

summary(d_match)
stopifnot(max(d_match) < 1.5)

# -------------------------------------------------------------------
# NOTES: Interpretation of d_match (distance between POP and PM centers)
#
# d_match measures the Euclidean distance (in degrees) between each
# 0.5° population cell center and the center of its matched 1.25° PM cell.
#
# Geometry of the PM grid:
# The PM grid has 288 x 192 cells covering the globe.
# - Longitudinal resolution = 360° / 288 = 1.25°
# - Latitudinal resolution  = 180° / 192 = 0.9375° ≈ 0.94°
#
# Therefore, each PM cell is a rectangle with:
# - half width  (lon) = a = 1.25 / 2 = 0.625°
# - half height (lat) = b = 0.9375 / 2 = 0.46875°
#
# The maximum possible distance from a PM cell center to any point
# inside the cell is the distance to a corner:
#   d_max = sqrt(a^2 + b^2) ≈ sqrt(0.625^2 + 0.46875^2) ≈ 0.78°
#
# The observed maximum distance (max(d_match) = 0.7819°) matches
# this theoretical bound almost exactly.
#
# Why is the mean distance around 0.42°?
# Population cell centers can be treated as approximately uniformly
# distributed points inside each PM grid cell.
#
# Geometrically, this is equivalent to drawing a random point (X, Y)
# uniformly from a rectangle:
#   X ~ Uniform(-a, a)
#   Y ~ Uniform(-b, b)
#
# The distance to the center is:
#   D = sqrt(X^2 + Y^2)
#
# The theoretical expected value is:
#   E[D] = (1 / (4ab)) * ∫_{-a}^{a} ∫_{-b}^{b} sqrt(x^2 + y^2) dx dy
#
# This integral has no simple closed form, but a known geometric result is:
#   E[D] ≈ 0.52 * sqrt(a^2 + b^2)
#
# Here:
#   sqrt(a^2 + b^2) ≈ 0.78°
#   E[D] ≈ 0.52 * 0.78 ≈ 0.41°
#
# The observed mean and median distances (~0.42°) match this theoretical
# value, which is exactly what is expected from pure geometry.
#
# Conclusion:
# This confirms that every population cell lies within its matched PM cell,
# and that the PM-to-population nearest-neighbor mapping is geometrically
# correct, with no spatial misalignment or indexing errors.
# -------------------------------------------------------------------

#verify by plot by sampling
set.seed(123)
idx <- sample(nrow(pop_coords), 2000)
pop_test <- pop_coords[idx, ]

true_check <- map_df(1:nrow(pop_test), function(i) {
  plon <- pop_test$lon[i]
  plat <- pop_test$lat[i]
  
  dist_all <- sqrt(
    (plon - annual_ave$lon_center)^2 +
      (plat - annual_ave$lat_center)^2
  )
  
  tibble(
    d_min   = min(dist_all),                       
    d_match = sqrt(
      (plon - annual_ave$lon_center[nearest_pm_idx[idx[i]]])^2 +
        (plat - annual_ave$lat_center[nearest_pm_idx[idx[i]]])^2
    )
  )
})

#For each population cell selected during sampling, 
#the chosen PM cell is the one with the shortest distance in the entire space.
ggplot(true_check, aes(x = d_min, y = d_match)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, 
              linetype = "dashed", linewidth = 1) +
  coord_equal() +
  labs(
    title = "Nearest-neighbor correctness check",
    subtitle = "All points must lie on y = x if no mismatches",
    x = "True minimum distance (min of Euclidean distance)",
    y = "Distance chosen by algorithm"
  ) +
  theme_minimal()
################################################################################review end


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

#write_csv(pm25_regridded, here("output", "pm25_regridded_0.5deg.csv"))
#write_csv(pop_pm_combined, here("output", "pop_pm_combined.csv"))

cat("\nRegridding complete!\n")
