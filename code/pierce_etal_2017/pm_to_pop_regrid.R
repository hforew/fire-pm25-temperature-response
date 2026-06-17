# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
########### Re-grid PM data to pop data resolution ##########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Goal: Re-grid coarse PM2.5 (1.25°) onto the 0.5° population grid by nearest neighbor,
#       then merge PM2.5 + population into one per-cell table.
#   1. Read annual-average PM2.5 (288x192 = 55,296 cells) and the pop grid
#      (720x360 = 259,200 cells); shift PM2.5 edge coords to cell centers
#      (+1.25/2 lon, +0.9424/2 lat) so both grids are center-referenced.
#   2. For each 0.5° pop cell, find the nearest 1.25° PM2.5 center by Euclidean
#      distance (which.min) and copy that cell's pm_/fpm_ columns over.
#   3. Join the regridded PM2.5 back onto the full pop table on exact lon/lat;
#      verify every pop row matched, and write both outputs.
# Inputs : output/annual_ave_pm25.csv ; output/pop_df_rev.csv
# Outputs: output/pm25_regridded_0.5deg.csv ; output/pop_pm_combined.csv


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
############ Regrid PM2.5 from 1.25° to 0.5° Resolution ############ 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Method: Nearest neighbor assignment (uniform assignment within coarse cells)
# Each 1.25° PM2.5 cell maps to multiple 0.5° population cells

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 1: Convert PM2.5 coordinates from edges to centers
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Current PM2.5 coordinates are grid edges (lon starts at -180°, lat at -90°)
# The population data coordinates are grid centers
# Need to shift PM2.5 to cell centers for proper regridding

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


## THE END 

