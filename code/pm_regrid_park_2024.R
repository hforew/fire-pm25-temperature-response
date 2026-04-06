# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
########### Import Park et al. 2024 PM2.5 Data by Model & Scenario ##########
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls())
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Packages #######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(ncdf4)
library(here)
library(tidyverse)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Import Population Data #########################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat("=== Importing population data ===\n\n")

pop_pm_combined <- read_csv(here("output", "pop_pm_combined.csv"))
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Function: Extract PM2.5 from NetCDF ############################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

extract_pm25 <- function(model, scenario, year) {
  
  file_path <- here("input", "Park_etal_2024", model, scenario,
                    paste0("05x05_CEDS_", year, "_on_off_pm25_Surface_Re_yearavg.nc4"))
  
  nc <- nc_open(file_path)
  lon <- ncvar_get(nc, "lon")
  lat <- ncvar_get(nc, "lat")
  pm25_var <- names(nc$var)[grep("pm25", names(nc$var), ignore.case = TRUE)][1]
  pm25 <- ncvar_get(nc, pm25_var)
  nc_close(nc)
  
  expand.grid(lon = lon, lat = lat) %>%
    mutate(pm25 = as.vector(pm25), year = year)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Import All 6 Datasets ##########################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

years <- c(1965, 1975, 1985, 1995, 2005, 2015)

# Dataset 1: classic_counterclim
classic_counterclim <- map_dfr(years, ~extract_pm25("classic", "counterclim", .x))

# Dataset 2: classic_obsclim
classic_obsclim <- map_dfr(years, ~extract_pm25("classic", "obsclim", .x))

# Dataset 3: jules_counterclim
jules_counterclim <- map_dfr(years, ~extract_pm25("jules", "counterclim", .x))

# Dataset 4: jules_obsclim
jules_obsclim <- map_dfr(years, ~extract_pm25("jules", "obsclim", .x))

# Dataset 5: ssib4_counterclim
ssib4_counterclim <- map_dfr(years, ~extract_pm25("ssib4", "counterclim", .x))

# Dataset 6: ssib4_obsclim
ssib4_obsclim <- map_dfr(years, ~extract_pm25("ssib4", "obsclim", .x))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Check if Coordinates are Cell Edges or Centers #################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat("=== Checking Coordinate Type ===\n\n")

# Extract unique lon and lat values
lon_vals <- sort(unique(classic_counterclim$lon))
lat_vals <- sort(unique(classic_counterclim$lat))

# Calculate grid spacing
lon_spacing <- diff(lon_vals)[1]
lat_spacing <- diff(lat_vals)[1]

cat("Longitude range:", min(lon_vals), "to", max(lon_vals), "\n")
cat("Latitude range:", min(lat_vals), "to", max(lat_vals), "\n")
cat("Longitude spacing:", lon_spacing, "degrees\n")
cat("Latitude spacing:", lat_spacing, "degrees\n\n")

# Check if coordinates align with cell edges or centers
# For 0.5° resolution:
# - Cell edges: -180, -179.5, -179, ... (starts at -180)
# - Cell centers: -179.75, -179.25, -178.75, ... (offset by 0.25)

cat("First longitude value:", lon_vals[1], "\n")
cat("Last longitude value:", tail(lon_vals, 1), "\n")
cat("First latitude value:", lat_vals[1], "\n")
cat("Last latitude value:", tail(lat_vals, 1), "\n\n")

# Determine if edges or centers
if (min(lon_vals) == -180 || max(lon_vals) == 180) {
  cat("DIAGNOSIS: Coordinates are CELL EDGES\n")
  cat("Converting to CELL CENTERS by adding half grid spacing...\n\n")
  is_edge <- TRUE
} else {
  cat("DIAGNOSIS: Coordinates are already CELL CENTERS\n")
  is_edge <- FALSE
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Convert Longitude from Edge to Center ##########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat("=== Checking Longitude Coordinates ===\n\n")

# Check current longitude values
lon_vals <- sort(unique(classic_counterclim$lon))
lon_spacing <- diff(lon_vals)[1]

cat("Current longitude range:", min(lon_vals), "to", max(lon_vals), "\n")
cat("Longitude spacing:", lon_spacing, "degrees\n")
cat("First few lon values:", head(lon_vals, 5), "\n\n")

# Target: -179.75, -179.25, -178.75, ... (cell centers)
# Current: -179.5, -179.0, -178.5, ... (cell edges)

cat("Converting longitude from cell edges to cell centers...\n")
cat("Shifting longitude by -0.25 degrees\n\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Apply Conversion to All 6 Datasets #############################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

classic_counterclim <- classic_counterclim %>%
  mutate(lon = lon - 0.25)

classic_obsclim <- classic_obsclim %>%
  mutate(lon = lon - 0.25)

jules_counterclim <- jules_counterclim %>%
  mutate(lon = lon - 0.25)

jules_obsclim <- jules_obsclim %>%
  mutate(lon = lon - 0.25)

ssib4_counterclim <- ssib4_counterclim %>%
  mutate(lon = lon - 0.25)

ssib4_obsclim <- ssib4_obsclim %>%
  mutate(lon = lon - 0.25)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Verify Conversion ##############################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat("=== Conversion Complete ===\n\n")

new_lon_vals <- sort(unique(classic_counterclim$lon))
cat("New longitude range:", min(new_lon_vals), "to", max(new_lon_vals), "\n")
cat("First few lon values:", head(new_lon_vals, 5), "\n\n")

cat("Sample data after conversion:\n")
print(head(classic_counterclim, 10))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Merge Counterclim and Obsclim to Calculate Fire PM2.5 ##########
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat("=== Merging scenarios and calculating fire PM2.5 ===\n\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Merge and Calculate fpm for Classic Model ######################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

classic <- classic_obsclim %>%
  rename(pm25_obsclim = pm25) %>%                      # Rename pm25 to pm25_obsclim
  left_join(
    classic_counterclim %>% 
      rename(pm25_counterclim = pm25),                 # Rename pm25 to pm25_counterclim
    by = c("lon", "lat", "year")                       # Join on coordinates and year
  ) %>%
  mutate(
    fpm25 = pm25_obsclim - pm25_counterclim            # Fire PM2.5 = observed - counterfactual
  ) %>%
  select(lon, lat, year, pm25_counterclim, pm25_obsclim, fpm25)  # Reorder columns

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Merge and Calculate fpm for JULES Model ########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

jules <- jules_obsclim %>%
  rename(pm25_obsclim = pm25) %>%
  left_join(
    jules_counterclim %>% 
      rename(pm25_counterclim = pm25),
    by = c("lon", "lat", "year")
  ) %>%
  mutate(
    fpm25 = pm25_obsclim - pm25_counterclim
  ) %>%
  select(lon, lat, year, pm25_counterclim, pm25_obsclim, fpm25)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Merge and Calculate fpm for SSIB4 Model ########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ssib4 <- ssib4_obsclim %>%
  rename(pm25_obsclim = pm25) %>%
  left_join(
    ssib4_counterclim %>% 
      rename(pm25_counterclim = pm25),
    by = c("lon", "lat", "year")
  ) %>%
  mutate(
    fpm25 = pm25_obsclim - pm25_counterclim
  ) %>%
  select(lon, lat, year, pm25_counterclim, pm25_obsclim, fpm25)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Verify Merged Datasets #########################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat("Classic dataset dimensions:", nrow(classic), "rows\n")
cat("JULES dataset dimensions:", nrow(jules), "rows\n")
cat("SSIB4 dataset dimensions:", nrow(ssib4), "rows\n\n")

cat("Sample of classic dataset:\n")
print(head(classic, 10))

cat("\n\nSummary of fire PM2.5 (fpm25) for each model:\n")
cat("\nClassic:\n")
print(summary(classic$fpm25))

cat("\nJULES:\n")
print(summary(jules$fpm25))

cat("\nSSIB4:\n")
print(summary(ssib4$fpm25))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Clean Up Individual Scenario Datasets ##########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Remove individual scenario datasets to save memory
rm(classic_counterclim, classic_obsclim, 
   jules_counterclim, jules_obsclim, 
   ssib4_counterclim, ssib4_obsclim)

cat("\n\n=== Merge complete! ===\n")
cat("3 final datasets: classic, jules, ssib4\n")
cat("Each contains: lon, lat, year, pm25_counterclim, pm25_obsclim, fpm25\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Reshape Data from Long to Wide Format by Year ##################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat("=== Reshaping datasets from long to wide format ===\n\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Reshape Classic Dataset ########################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

classic_wide <- classic %>%
  pivot_wider(
    id_cols = c(lon, lat),                             # Keep lon and lat as identifiers
    names_from = year,                                 # Create columns for each year
    values_from = c(pm25_counterclim, pm25_obsclim, fpm25),  # Spread these variables
    names_glue = "{.value}_{year}_classic"             # Column naming: variable_year_classic
  )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Reshape JULES Dataset ##########################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

jules_wide <- jules %>%
  pivot_wider(
    id_cols = c(lon, lat),
    names_from = year,
    values_from = c(pm25_counterclim, pm25_obsclim, fpm25),
    names_glue = "{.value}_{year}_jules"               # Column naming: variable_year_jules
  )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Reshape SSIB4 Dataset ##########################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ssib4_wide <- ssib4 %>%
  pivot_wider(
    id_cols = c(lon, lat),
    names_from = year,
    values_from = c(pm25_counterclim, pm25_obsclim, fpm25),
    names_glue = "{.value}_{year}_ssib4"               # Column naming: variable_year_ssib4
  )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Verify Wide Format #############################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat("Classic wide dimensions:", nrow(classic_wide), "rows x", ncol(classic_wide), "columns\n")
cat("JULES wide dimensions:", nrow(jules_wide), "rows x", ncol(jules_wide), "columns\n")
cat("SSIB4 wide dimensions:", nrow(ssib4_wide), "rows x", ncol(ssib4_wide), "columns\n\n")

cat("Expected rows per dataset:", 720 * 360, "(720 lon × 360 lat)\n")
cat("Expected columns per dataset:", 2 + 3*6, "(lon, lat + 18 variables)\n\n")

cat("Column names in classic_wide:\n")
print(colnames(classic_wide))

cat("\n\nSample of classic_wide:\n")
print(head(classic_wide))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Clean Up Long Format Datasets ##################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Remove long format datasets to save memory
rm(classic, jules, ssib4)

cat("\n\n=== Reshape complete! ===\n")
cat("3 wide datasets created: classic_wide, jules_wide, ssib4_wide\n")
cat("Each contains:\n")
cat("  - lon, lat (identifiers)\n")
cat("  - pm25_counterclim_YEAR_MODEL (6 columns)\n")
cat("  - pm25_obsclim_YEAR_MODEL (6 columns)\n")
cat("  - fpm25_YEAR_MODEL (6 columns)\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Merge Three Models and Rename Variables ########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat("=== Merging three models into single dataset ===\n\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Step 1: Merge Classic, JULES, and SSIB4 ########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pm_regrid_park <- classic_wide %>%
  left_join(jules_wide, by = c("lon", "lat")) %>%      # Add JULES columns
  left_join(ssib4_wide, by = c("lon", "lat"))          # Add SSIB4 columns

cat("Merged dimensions:", nrow(pm_regrid_park), "rows x", ncol(pm_regrid_park), "columns\n\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Step 2: Remove pm25_obsclim Variables ##########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pm_regrid_park <- pm_regrid_park %>%
  select(-starts_with("pm25_obsclim"))                 # Drop all obsclim columns

cat("After removing obsclim:", ncol(pm_regrid_park), "columns\n\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Step 3: Rename pm25_counterclim to pm25 ########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Rename pm25_counterclim_YEAR_MODEL to pm25_YEAR_MODEL
pm_regrid_park <- pm_regrid_park %>%
  rename_with(
    ~ str_replace(.x, "pm25_counterclim_", "pm25_"),   # Remove "counterclim" from name
    starts_with("pm25_counterclim")                    # Only for counterclim columns
  )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Step 4: Add _Park Suffix to All Variables ######################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pm_regrid_park <- pm_regrid_park %>%
  rename_with(
    ~ paste0(.x, "_Park"),                             # Add _Park suffix
    -c(lon, lat)                                       # Except lon and lat
  )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Verify Final Dataset ###########################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat("\n=== Final dataset structure ===\n\n")
cat("Dataset name: pm_regrid_park\n")
cat("Dimensions:", nrow(pm_regrid_park), "rows x", ncol(pm_regrid_park), "columns\n\n")

cat("Column names:\n")
print(colnames(pm_regrid_park))

cat("\n\nSample data:\n")
print(head(pm_regrid_park))

cat("\n\nVariable count:\n")
cat("  lon, lat: 2 columns\n")
cat("  pm25_YEAR_MODEL_Park: 18 columns (6 years × 3 models)\n")
cat("  fpm25_YEAR_MODEL_Park: 18 columns (6 years × 3 models)\n")
cat("  Total: 38 columns\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Clean Up Individual Model Datasets #############################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(classic_wide, jules_wide, ssib4_wide)

cat("\n=== Merge and rename complete! ===\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Join Park Data to Population Data ##############################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat("=== Joining Park et al. 2024 data to population data ===\n\n")

# Check dimensions before join
cat("Before join:\n")
cat("  pop_pm_combined:", nrow(pop_pm_combined), "rows\n")
cat("  pm_regrid_park:", nrow(pm_regrid_park), "rows\n\n")

# Perform left join
pop_pm_combined_withPark <- pop_pm_combined %>%
  left_join(pm_regrid_park, by = c("lon", "lat"))      # Join on coordinates

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Verify Join ####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat("After join:\n")
cat("  pop_pm_combined_withPark:", nrow(pop_pm_combined_withPark), "rows x", 
    ncol(pop_pm_combined_withPark), "columns\n\n")

# Check for missing values (indicates join issues)
park_vars <- colnames(pop_pm_combined_withPark)[grep("_Park$", colnames(pop_pm_combined_withPark))]
missing_count <- sum(is.na(pop_pm_combined_withPark[[park_vars[1]]]))

cat("Missing values in Park data:", missing_count, "\n")
cat("Join success rate:", 
    round((1 - missing_count/nrow(pop_pm_combined_withPark)) * 100, 2), "%\n\n")

cat("Sample of joined data:\n")
print(head(pop_pm_combined_withPark))

cat("\n\nColumn names:\n")
print(colnames(pop_pm_combined_withPark))

cat("\n=== Join complete! ===\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Save Final Dataset #############################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat("=== Saving pop_pm_combined_withPark to output folder ===\n\n")

# Save as CSV
write_csv(pop_pm_combined_withPark, 
          here("output", "pop_pm_combined_withPark.csv"))

cat("File saved to:", here("output", "pop_pm_combined_withPark.csv"), "\n")
cat("File size:", 
    round(file.size(here("output", "pop_pm_combined_withPark.csv")) / 1024^2, 2), 
    "MB\n\n")

cat("=== Save complete! ===\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
######################### THE END############################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~