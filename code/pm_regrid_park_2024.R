# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
########### Extract Year, Longitude, Latitude, and PM2.5 from NetCDF Files ##################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Remove all objects from the environment
rm(list = ls())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Packages #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(here)
library(tidyverse)
library(ncdf4)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ File Paths for All Scenarios ##################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

years <- c(1965, 1975, 1985, 1995, 2005, 2015)

# Scenario 1: withoutfire
file_paths_withoutfire <- c(
  here("input", "Park_etal_2024", "GEOSChem_output", "withoutfire", "05x05_CEDS_1965_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "withoutfire", "05x05_CEDS_1975_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "withoutfire", "05x05_CEDS_1985_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "withoutfire", "05x05_CEDS_1995_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "withoutfire", "05x05_CEDS_2005_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "withoutfire", "05x05_CEDS_2015_on_off_pm25_Surface_Re_yearavg.nc4")
)

# Scenario 2: classic
file_paths_classic <- c(
  here("input", "Park_etal_2024", "GEOSChem_output", "classic", "obsclim", "05x05_CEDS_1965_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "classic", "obsclim", "05x05_CEDS_1975_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "classic", "obsclim", "05x05_CEDS_1985_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "classic", "obsclim", "05x05_CEDS_1995_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "classic", "obsclim", "05x05_CEDS_2005_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "classic", "obsclim", "05x05_CEDS_2015_on_off_pm25_Surface_Re_yearavg.nc4")
)

# Scenario 3: jules
file_paths_jules <- c(
  here("input", "Park_etal_2024", "GEOSChem_output", "jules", "obsclim", "05x05_CEDS_1965_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "jules", "obsclim", "05x05_CEDS_1975_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "jules", "obsclim", "05x05_CEDS_1985_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "jules", "obsclim", "05x05_CEDS_1995_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "jules", "obsclim", "05x05_CEDS_2005_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "jules", "obsclim", "05x05_CEDS_2015_on_off_pm25_Surface_Re_yearavg.nc4")
)

# Scenario 4: ssib4
file_paths_ssib4 <- c(
  here("input", "Park_etal_2024", "GEOSChem_output", "ssib4", "obsclim", "05x05_CEDS_1965_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "ssib4", "obsclim", "05x05_CEDS_1975_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "ssib4", "obsclim", "05x05_CEDS_1985_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "ssib4", "obsclim", "05x05_CEDS_1995_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "ssib4", "obsclim", "05x05_CEDS_2005_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "ssib4", "obsclim", "05x05_CEDS_2015_on_off_pm25_Surface_Re_yearavg.nc4")
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Function to Convert Grid Edge to Grid Center ##################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

convert_edge_to_center <- function(coords) {
  # Check if coordinates are grid edges or centers
  # Grid edges: start at exact values like -180, 0, etc.
  # Grid centers: start at offset values like -179.75, 0.25, etc.
  
  # Get the first few values
  first_val <- coords[1]
  resolution <- mean(diff(coords))
  
  # Check if first value is a multiple of resolution (edge) or offset (center)
  # For 0.5 degree grid:
  # Edge starts at -180, -179.5, -179, ...
  # Center starts at -179.75, -179.25, -178.75, ...
  
  # If coordinates are already at center, return as is
  # If at edge, shift by half resolution
  
  # Simple check: if first value modulo resolution is close to 0, it's edge
  remainder <- abs(first_val %% resolution)
  
  if (remainder < 0.01 || abs(remainder - resolution) < 0.01) {
    # Grid edge detected - convert to center
    coords_center <- coords + (resolution / 2)
    return(coords_center)
  } else {
    # Already at grid center
    return(coords)
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Function to Extract Data #######################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

extract_pm25_data <- function(file_paths, years, scenario_name) {
  
  pm25_list <- list()
  
  for (i in 1:length(file_paths)) {
    
    # Open NetCDF file
    nc <- nc_open(file_paths[i])
    
    # Extract coordinates - use seq to ensure correct grid
    # For 0.5 degree resolution: -179.75 to 179.75
    lon <- seq(-179.75, 179.75, by = 0.5)  # 720 values
    lat <- seq(-89.75, 89.75, by = 0.5)    # 360 values
    
    # Extract PM2.5 data
    pm25_raw <- ncvar_get(nc, "PM25")
    
    # Check dimensions and extract surface data accordingly
    dims <- dim(pm25_raw)
    
    if (length(dims) == 4) {
      # 4D array: [lon, lat, lev, time]
      pm25_surface <- pm25_raw[, , 1, 1]
    } else if (length(dims) == 3) {
      # 3D array: could be [lon, lat, time] or [lon, lat, lev]
      pm25_surface <- pm25_raw[, , 1]
    } else if (length(dims) == 2) {
      # 2D array: [lon, lat]
      pm25_surface <- pm25_raw
    } else {
      stop("Unexpected PM25 dimensions: ", paste(dims, collapse = " x "))
    }
    
    # Close NetCDF file
    nc_close(nc)
    
    # Create data frame with lon varying first
    df <- expand.grid(
      lon = lon,
      lat = lat,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    
    # Add PM2.5 - transpose if needed to match grid order
    df$pm25 <- as.vector(pm25_surface)
    df$year <- years[i]
    
    # Reorder columns
    df <- df %>%
      select(year, lon, lat, pm25)
    
    # Store in list
    pm25_list[[i]] <- df
  }
  
  # Combine all years
  pm25_data <- bind_rows(pm25_list)
  
  # Rename pm25 column to include scenario name
  pm25_data <- pm25_data %>%
    rename(!!paste0("pm25_", scenario_name) := pm25)
  
  return(pm25_data)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Extract Data for All Scenarios #################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Extract data for each scenario
withoutfire <- extract_pm25_data(file_paths_withoutfire, years, "withoutfire")
classic <- extract_pm25_data(file_paths_classic, years, "classic")
jules <- extract_pm25_data(file_paths_jules, years, "jules")
ssib4 <- extract_pm25_data(file_paths_ssib4, years, "ssib4")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Merge All Scenarios by Year, Lon, Lat #########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Merge all four datasets by year, lon, lat
pm25_all <- withoutfire %>%
  left_join(classic, by = c("year", "lon", "lat")) %>%
  left_join(jules, by = c("year", "lon", "lat")) %>%
  left_join(ssib4, by = c("year", "lon", "lat"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Calculate Fire PM (Difference from Without Fire) ##############
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Calculate fire PM for each model by subtracting withoutfire PM
# Fire PM = Model PM - Without Fire PM
pm25_all <- pm25_all %>%
  mutate(
    classic_fpm = pm25_classic - pm25_withoutfire,
    jules_fpm = pm25_jules - pm25_withoutfire,
    ssib4_fpm = pm25_ssib4 - pm25_withoutfire
  )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Create Final Two Datasets ######################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Dataset 1: Complete dataset with all PM columns and fire PM columns
# Columns: year, lon, lat, pm25_withoutfire, pm25_classic, pm25_jules, pm25_ssib4, 
#          classic_fpm, jules_fpm, ssib4_fpm
pm25_all_complete <- pm25_all

# Dataset 2: Dataset with only fire PM columns (removed all pm25_* columns)
# Columns: year, lon, lat, classic_fpm, jules_fpm, ssib4_fpm
pm25_fpm_only <- pm25_all %>%
  select(year, lon, lat, classic_fpm, jules_fpm, ssib4_fpm)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Convert Fire PM to Wide Format by Year ########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Convert pm25_fpm_only from long to wide format
# Each model x year combination becomes a separate column
# Example: classic_1965_fpm, classic_1975_fpm, ..., jules_1965_fpm, etc.

pm25_fpm_wide <- pm25_fpm_only %>%
  pivot_longer(
    cols = c(classic_fpm, jules_fpm, ssib4_fpm),
    names_to = "model",
    values_to = "fpm"
  ) %>%
  mutate(
    # Create column name: model_year_fpm (e.g., classic_1965_fpm)
    model_year = paste0(gsub("_fpm", "", model), "_", year, "_fpm")
  ) %>%
  select(-model, -year) %>%
  pivot_wider(
    names_from = model_year,
    values_from = fpm
  ) %>%
  # Add "park_" prefix to all columns except lon and lat
  rename_with(
    ~ paste0("park_", .x),
    .cols = -c(lon, lat)
  )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ merge with pop_regrid_park_2024 and our master data ############
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pop_pm_combined <- read.csv(here("output", "pop_pm_combined.csv"))
pop_park <- read.csv(here("output", "pop_regrid_park_2024.csv"))

# Merge pop_wide with pm25_fpm_wide by lon and lat
pm25_pop_merged <- pm25_fpm_wide %>%
  left_join(pop_park, by = c("lon", "lat"))

# Merge pm25_pop_merged with pop_pm_combined by lon and lat
combined_data <- pop_pm_combined %>%
  left_join(pm25_pop_merged, by = c("lon", "lat"))

# Check the results
dim(combined_data)  # Check dimensions
head(combined_data)  # View first few rows
summary(combined_data)  # Summary statistics

# Check for missing values after merge (indicates non-matching rows)
sum(is.na(combined_data))  # Total NA count

# Check if coordinates match perfectly
anti_join(pop_pm_combined, pm25_pop_merged, by = c("lon", "lat"))  # Rows in pop_pm_combined not in pm25_pop_merged
anti_join(pm25_pop_merged, pop_pm_combined, by = c("lon", "lat"))  # Rows in pm25_pop_merged not in pop_pm_combined

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
######################## save final master data to local ####################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

write.csv(combined_data, here("output", "pop_pm_combined_with_park2024.csv"), row.names = FALSE)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ THE END ########################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# THE END