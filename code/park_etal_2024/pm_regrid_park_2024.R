# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
########### Extract Year, Longitude, Latitude, and PM2.5 from NetCDF Files ##################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Process GEOS-Chem PM2.5 output from Park et al. (2024) to derive gridded
# fire-attributable PM2.5 (fPM) and merge it into the master dataset.
#
# Inputs : Decade-mean surface PM2.5 NetCDF files for 6 periods
#          (1965, 1975, 1985, 1995, 2005, 2015) under 4 scenarios:
#          withoutfire (no-fire baseline), classic, jules, ssib4 (fire models).
#
# Steps  : (1) Read NetCDFs, shift lon to [-180,180], regrid to 0.5° x 0.5°.
#          (2) Compute fPM = PM2.5(models) - PM2.5(withoutfire) for each model.
#          (3) Reshape to wide format (one column per model x year).
#          (4) Merge with population data and the master pop_pm dataset.
#          (5) Validate against authors' results1206.mat; render CONUS sanity-
#              check maps; diagnose the native 4° x 5° simulation resolution.

# Remove all objects from the environment

rm(list = ls())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Packages #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(here)
library(tidyverse)
library(ncdf4)
library(R.matlab)
library(fields)

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

# Scenario 2: classic factual
file_paths_classic <- c(
  here("input", "Park_etal_2024", "GEOSChem_output", "classic", "obsclim", "05x05_CEDS_1965_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "classic", "obsclim", "05x05_CEDS_1975_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "classic", "obsclim", "05x05_CEDS_1985_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "classic", "obsclim", "05x05_CEDS_1995_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "classic", "obsclim", "05x05_CEDS_2005_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "classic", "obsclim", "05x05_CEDS_2015_on_off_pm25_Surface_Re_yearavg.nc4")
)

# Scenario 3: jules factual
file_paths_jules <- c(
  here("input", "Park_etal_2024", "GEOSChem_output", "jules", "obsclim", "05x05_CEDS_1965_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "jules", "obsclim", "05x05_CEDS_1975_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "jules", "obsclim", "05x05_CEDS_1985_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "jules", "obsclim", "05x05_CEDS_1995_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "jules", "obsclim", "05x05_CEDS_2005_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "jules", "obsclim", "05x05_CEDS_2015_on_off_pm25_Surface_Re_yearavg.nc4")
)

# Scenario 4: ssib4 factual
file_paths_ssib4 <- c(
  here("input", "Park_etal_2024", "GEOSChem_output", "ssib4", "obsclim", "05x05_CEDS_1965_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "ssib4", "obsclim", "05x05_CEDS_1975_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "ssib4", "obsclim", "05x05_CEDS_1985_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "ssib4", "obsclim", "05x05_CEDS_1995_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "ssib4", "obsclim", "05x05_CEDS_2005_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "ssib4", "obsclim", "05x05_CEDS_2015_on_off_pm25_Surface_Re_yearavg.nc4")
)

# Scenario 5: classic counterclim
file_paths_classic_cf <- c(
  here("input", "Park_etal_2024", "GEOSChem_output", "classic", "counterclim", "05x05_CEDS_1965_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "classic", "counterclim", "05x05_CEDS_1975_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "classic", "counterclim", "05x05_CEDS_1985_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "classic", "counterclim", "05x05_CEDS_1995_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "classic", "counterclim", "05x05_CEDS_2005_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "classic", "counterclim", "05x05_CEDS_2015_on_off_pm25_Surface_Re_yearavg.nc4")
)

# Scenario 6: jules counterclim
file_paths_jules_cf <- c(
  here("input", "Park_etal_2024", "GEOSChem_output", "jules", "counterclim", "05x05_CEDS_1965_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "jules", "counterclim", "05x05_CEDS_1975_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "jules", "counterclim", "05x05_CEDS_1985_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "jules", "counterclim", "05x05_CEDS_1995_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "jules", "counterclim", "05x05_CEDS_2005_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "jules", "counterclim", "05x05_CEDS_2015_on_off_pm25_Surface_Re_yearavg.nc4")
)

# Scenario 7: ssib4 counterclim
file_paths_ssib4_cf <- c(
  here("input", "Park_etal_2024", "GEOSChem_output", "ssib4", "counterclim", "05x05_CEDS_1965_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "ssib4", "counterclim", "05x05_CEDS_1975_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "ssib4", "counterclim", "05x05_CEDS_1985_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "ssib4", "counterclim", "05x05_CEDS_1995_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "ssib4", "counterclim", "05x05_CEDS_2005_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "ssib4", "counterclim", "05x05_CEDS_2015_on_off_pm25_Surface_Re_yearavg.nc4")
)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Target Grid (0.5 deg x 0.5 deg centers) #######################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

target_lon <- seq(-179.75, 179.75, by = 0.5)  # 720
target_lat <- seq(-89.75,  89.75,  by = 0.5)  # 360

convert_edge_to_center <- function(coords) {
  return(coords)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ 0-360 lon -> -180-180 ##########################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

shift_lon_to_180 <- function(lon, mat) {
  if (max(lon) > 180) {
    lon_new <- ifelse(lon > 180, lon - 360, lon)
    ord     <- order(lon_new)
    lon_new <- lon_new[ord]
    mat     <- mat[ord, , drop = FALSE]
    return(list(lon = lon_new, mat = mat))
  }
  list(lon = lon, mat = mat)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Function to Extract + Regrid ##################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
extract_pm25_data <- function(file_paths, years, scenario_name,
                              target_lon, target_lat) {
  
  pm25_list <- list()
  
  for (i in seq_along(file_paths)) {
    
    nc <- nc_open(file_paths[i])
    
    # 1. Read actual lon/lat from the file
    lon_name <- if ("lon" %in% names(nc$dim)) "lon" else "longitude"
    lat_name <- if ("lat" %in% names(nc$dim)) "lat" else "latitude"
    lon_src  <- as.numeric(ncvar_get(nc, lon_name))
    lat_src  <- as.numeric(ncvar_get(nc, lat_name))
    
    cat("lon range:", range(lon_src), " | n =", length(lon_src), "\n")
    cat("lat range:", range(lat_src), " | n =", length(lat_src), "\n")
    cat("first 5 lon:", head(lon_src, 5), "\n")
    cat("first 5 lat:", head(lat_src, 5), "\n")
    
    # Edge -> center if needed (no-op for Park et al. files)
    lon_src <- convert_edge_to_center(lon_src)
    lat_src <- convert_edge_to_center(lat_src)
    
    # 2. Read PM2.5 and reduce to 2D
    pm25_raw <- ncvar_get(nc, "PM25")
    dims     <- dim(pm25_raw)
    pm25_surface <- switch(as.character(length(dims)),
                           "4" = pm25_raw[, , 1, 1],
                           "3" = pm25_raw[, , 1],
                           "2" = pm25_raw,
                           stop("Unexpected PM25 dimensions: ", paste(dims, collapse = " x "))
    )
    nc_close(nc)
    
    # 3. Ensure orientation is [lon, lat] (rows = lon, cols = lat)
    if (nrow(pm25_surface) == length(lat_src) &&
        ncol(pm25_surface) == length(lon_src)) {
      pm25_surface <- t(pm25_surface)
    }
    stopifnot(nrow(pm25_surface) == length(lon_src),
              ncol(pm25_surface) == length(lat_src))
    
    # 4. Convert lon to -180..180 if needed
    sh           <- shift_lon_to_180(lon_src, pm25_surface)
    lon_src      <- sh$lon
    pm25_surface <- sh$mat
    
    # 5. Make latitude ascending (required by interp.surface.grid)
    if (is.unsorted(lat_src)) {
      ord          <- order(lat_src)
      lat_src      <- lat_src[ord]
      pm25_surface <- pm25_surface[, ord, drop = FALSE]
    }
    
    # 5a. Pad latitude on both poles so target endpoints (-89.75, 89.75) fall
    # strictly inside the source lat range. Source lat endpoints equal target
    # endpoints, which can trigger out-of-bounds errors in interp.surface.grid
    # due to numerical precision. lat is not periodic so we use nearest-neighbor.
    # Matrix is [lon, lat], so lat is along columns
    lat_resolution <- mean(diff(lat_src))
    pm25_surface <- cbind(
      pm25_surface[, 1],
      pm25_surface,
      pm25_surface[, ncol(pm25_surface)]
    )
    lat_src <- c(min(lat_src) - lat_resolution,
                 lat_src,
                 max(lat_src) + lat_resolution)
    
    # 5b. Wrap longitude on both sides for Earth periodicity (lon = 180 == -180).
    # Target lon -179.75..179.75 must fall inside source lon -179.5..180.
    # Matrix is [lon, lat], so lon is along rows
    if (max(lon_src) >= 180 - 1e-6 && min(lon_src) > -180 + 1e-6) {
      lon_resolution <- mean(diff(lon_src))
      pm25_surface <- rbind(
        pm25_surface[nrow(pm25_surface), ],   # copy lon = 180  -> place at -180
        pm25_surface,
        pm25_surface[1, ]                      # copy lon = -179.5 -> place at +180.5
      )
      lon_src <- c(-180, lon_src, max(lon_src) + lon_resolution)
    }
    
    # 6. Regrid to target 0.5 deg grid 
    same_grid <-
      length(lon_src) == length(target_lon) &&
      length(lat_src) == length(target_lat) &&
      max(abs(lon_src - target_lon)) < 1e-6 &&
      max(abs(lat_src - target_lat)) < 1e-6
    
    if (same_grid) {
      pm25_regrid <- pm25_surface
    } else {
      obj <- list(x = lon_src, y = lat_src, z = pm25_surface)
      pm25_regrid <- interp.surface.grid(
        obj,
        grid.list = list(x = target_lon, y = target_lat)
      )$z
    }
    
    # 7. Long format (lon varies first)
    df <- expand.grid(
      lon = target_lon,
      lat = target_lat,
      KEEP.OUT.ATTRS  = FALSE,
      stringsAsFactors = FALSE
    )
    df$pm25 <- as.vector(pm25_regrid)
    df$year <- years[i]
    df <- df %>% select(year, lon, lat, pm25)
    pm25_list[[i]] <- df
    
    cat(sprintf("[%s] %d: src %d x %d -> target %d x %d %s\n",
                scenario_name, years[i],
                length(lon_src), length(lat_src),
                length(target_lon), length(target_lat),
                if (same_grid) "(no regrid)" else "(bilinear regridded)"))
  }
  
  bind_rows(pm25_list) %>%
    rename(!!paste0("pm25_", scenario_name) := pm25)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Extract Data for All Scenarios #################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

withoutfire    <- extract_pm25_data(file_paths_withoutfire, years, "withoutfire",
                                    target_lon, target_lat)
classic        <- extract_pm25_data(file_paths_classic,     years, "classic",
                                    target_lon, target_lat)
jules          <- extract_pm25_data(file_paths_jules,        years, "jules",
                                    target_lon, target_lat)
ssib4          <- extract_pm25_data(file_paths_ssib4,        years, "ssib4",
                                    target_lon, target_lat)
classic_counter <- extract_pm25_data(file_paths_classic_cf, years, "classic_counter",
                                     target_lon, target_lat)
jules_counter   <- extract_pm25_data(file_paths_jules_cf,   years, "jules_counter",
                                     target_lon, target_lat)
ssib4_counter   <- extract_pm25_data(file_paths_ssib4_cf,   years, "ssib4_counter",
                                     target_lon, target_lat)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ SANITY CHECK: PM2.5 by year for all 3 models ##################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Force 6 decimal places in console output
options(pillar.sigfig = 7, digits = 7)

# factual
classic %>%
  group_by(year) %>%
  summarise(
    min    = min   (pm25_classic, na.rm = TRUE),
    max    = max   (pm25_classic, na.rm = TRUE),
    mean   = mean  (pm25_classic, na.rm = TRUE),
    median = median(pm25_classic, na.rm = TRUE)
  ) %>%
  print()

# counterfactual
classic_counter %>%
  group_by(year) %>%
  summarise(
    min    = min   (pm25_classic_counter, na.rm = TRUE),
    max    = max   (pm25_classic_counter, na.rm = TRUE),
    mean   = mean  (pm25_classic_counter, na.rm = TRUE),
    median = median(pm25_classic_counter, na.rm = TRUE)
  ) %>%
  print()

# factual
jules %>%
  group_by(year) %>%
  summarise(
    min    = min   (pm25_jules, na.rm = TRUE),
    max    = max   (pm25_jules, na.rm = TRUE),
    mean   = mean  (pm25_jules, na.rm = TRUE),
    median = median(pm25_jules, na.rm = TRUE)
  ) %>%
  print()

# counterfactual
jules_counter %>%
  group_by(year) %>%
  summarise(
    min    = min   (pm25_jules_counter, na.rm = TRUE),
    max    = max   (pm25_jules_counter, na.rm = TRUE),
    mean   = mean  (pm25_jules_counter, na.rm = TRUE),
    median = median(pm25_jules_counter, na.rm = TRUE)
  ) %>%
  print()

# factual
ssib4 %>%
  group_by(year) %>%
  summarise(
    min    = min   (pm25_ssib4, na.rm = TRUE),
    max    = max   (pm25_ssib4, na.rm = TRUE),
    mean   = mean  (pm25_ssib4, na.rm = TRUE),
    median = median(pm25_ssib4, na.rm = TRUE)
  ) %>%
  print()

# counterfactual
ssib4_counter %>%
  group_by(year) %>%
  summarise(
    min    = min   (pm25_ssib4_counter, na.rm = TRUE),
    max    = max   (pm25_ssib4_counter, na.rm = TRUE),
    mean   = mean  (pm25_ssib4_counter, na.rm = TRUE),
    median = median(pm25_ssib4_counter, na.rm = TRUE)
  ) %>%
  print()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Merge All Scenarios by Year, Lon, Lat #########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Merge all scenarios by year, lon, lat
pm25_all <- withoutfire %>%
  left_join(classic,         by = c("year", "lon", "lat")) %>%
  left_join(jules,           by = c("year", "lon", "lat")) %>%
  left_join(ssib4,           by = c("year", "lon", "lat")) %>%
  left_join(classic_counter, by = c("year", "lon", "lat")) %>%
  left_join(jules_counter,   by = c("year", "lon", "lat")) %>%
  left_join(ssib4_counter,   by = c("year", "lon", "lat"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ PM2.5 Wide Format (separate dataset, PM not fPM) ##############
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pm25_wide <- pm25_all %>%
  select(year, lon, lat, pm25_withoutfire, pm25_classic, pm25_jules, pm25_ssib4) %>%
  pivot_longer(
    cols = c(pm25_withoutfire, pm25_classic, pm25_jules, pm25_ssib4),
    names_to = "model",
    values_to = "pm25"
  ) %>%
  mutate(
    model_year = paste0(gsub("pm25_", "", model), "_", year, "_pm25")
  ) %>%
  select(-model, -year) %>%
  pivot_wider(
    names_from = model_year,
    values_from = pm25
  ) %>%
  rename_with(
    ~ paste0("park_", .x),
    .cols = -c(lon, lat)
  ) %>%
  # Replace years like 1965 with their corresponding decades (e.g., 1960s) to ensure consistency with the FPM dataset.
  rename_with(
    .cols = starts_with("park_"),
    .fn = function(x) {
      for (y in years) {
        x <- str_replace_all(x, as.character(y), paste0(floor(y / 10) * 10, "s"))
      }
      x
    }
  )

#write.csv(pm25_wide, here("output", "pm25_park2024.csv"), row.names = FALSE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Calculate Fire PM (Difference from Without Fire) ##############
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Fire PM = Model PM - Without Fire PM
# Climate-attributable fire PM = factual - counterfactual
pm25_all <- pm25_all %>%
  mutate(
    classic_fpm     = pm25_classic - pm25_withoutfire,
    jules_fpm       = pm25_jules   - pm25_withoutfire,
    ssib4_fpm       = pm25_ssib4   - pm25_withoutfire,
    classic_fpm_cli = pm25_classic - pm25_classic_counter,
    jules_fpm_cli   = pm25_jules   - pm25_jules_counter,
    ssib4_fpm_cli   = pm25_ssib4   - pm25_ssib4_counter
  )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Create Final Two Datasets ######################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Dataset 1: Complete dataset with all PM columns and fire PM columns
pm25_all_complete <- pm25_all

# Dataset 2: Dataset with only fire PM columns (removed all pm25_* columns)
# Columns: year, lon, lat, classic_fpm, jules_fpm, ssib4_fpm,
#          classic_fpm_cli, jules_fpm_cli, ssib4_fpm_cli
pm25_fpm_only <- pm25_all %>%
  select(year, lon, lat,
         classic_fpm, jules_fpm, ssib4_fpm,
         classic_fpm_cli, jules_fpm_cli, ssib4_fpm_cli)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Convert Fire PM to Wide Format by Year ########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Convert pm25_fpm_only from long to wide format
# Each model x year combination becomes a separate column
# fpm:     classic_1965_fpm, jules_1965_fpm, ...
# fpm_cli: classic_1965_fpm_cli, jules_1965_fpm_cli, ...
# Two separate pivots to avoid gsub collision between _fpm and _fpm_cli

fpm_wide <- pm25_fpm_only %>%
  select(year, lon, lat, classic_fpm, jules_fpm, ssib4_fpm) %>%
  pivot_longer(
    cols = c(classic_fpm, jules_fpm, ssib4_fpm),
    names_to = "model",
    values_to = "fpm"
  ) %>%
  mutate(model_year = paste0(gsub("_fpm", "", model), "_", year, "_fpm")) %>%
  select(-model, -year) %>%
  pivot_wider(names_from = model_year, values_from = fpm)

fpm_cli_wide <- pm25_fpm_only %>%
  select(year, lon, lat, classic_fpm_cli, jules_fpm_cli, ssib4_fpm_cli) %>%
  pivot_longer(
    cols = c(classic_fpm_cli, jules_fpm_cli, ssib4_fpm_cli),
    names_to = "model",
    values_to = "fpm_cli"
  ) %>%
  mutate(model_year = paste0(gsub("_fpm_cli", "", model), "_", year, "_fpm_cli")) %>%
  select(-model, -year) %>%
  pivot_wider(names_from = model_year, values_from = fpm_cli)

pm25_fpm_wide <- fpm_wide %>%
  left_join(fpm_cli_wide, by = c("lon", "lat")) %>%
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
############ Validation: Convert Back to Matrix ############################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Select park_classic_2015_fpm for validation
fpm_vector <- combined_data$park_classic_2015_fpm

# Convert back to matrix format (720 lon x 360 lat)
fpm_matrix <- matrix(fpm_vector, nrow = 720, ncol = 360)

# Check results
cat("\n=== Validation: Converted to Matrix ===\n")
cat("Matrix dimensions:", dim(fpm_matrix), "\n")
cat("\nFirst 5x5 corner:\n")
print(fpm_matrix[1:5, 1:5])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Validation: Compare with Original MAT File ####################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Purpose: Verify that our processed Classic 2015 fire PM2.5 (NetCDF -> regrid
# -> classic minus withoutfire) matches the authors' pre-computed PMfire data
# in results1206.mat, as a sanity check on the full pipeline.

# Read the MAT file
mat_data <- readMat(here("input", "Park_etal_2024", "results1206.mat"))

# Extract Classic model data
PMfire_classic <- mat_data$PMfire[[1]][[1]]  # Classic model: 360(lat, north -> south) x 720(lon) x 6(years) x 2(scenarios)
cat("\n=== Original MAT File ===\n")
cat("PMfire_classic dimensions:", dim(PMfire_classic), "\n")

# Extract 2015 data (year 6, scenario 1)
PMfire_classic_2015 <- PMfire_classic[, , 6, 1]  # 360 x 720

# Convert our data to matrix and transpose to match
# expand.grid(lon, lat) -> lon varies first, so nrow = 720 (lon), ncol = 360 (lat).
# Transpose to get lat as rows, lon as cols (standard raster convention).
fpm_vector <- combined_data$park_classic_2015_fpm
fpm_matrix <- matrix(fpm_vector, nrow = 720, ncol = 360)
fpm_matrix_t <- t(fpm_matrix)  # Now 360 (lat: south -> north) x 720 (lon)

# Flip lat: our data is south -> north, but the MAT file stores lat north -> south.
# This is comparison-only; combined_data itself keeps the standard south -> north order.
fpm_matrix_t <- fpm_matrix_t[360:1, ]

# Round to some decimal places for comparison
PMfire_rounded <- round(PMfire_classic_2015, 6)
fpm_rounded <- round(fpm_matrix_t, 6)

# Compare
max_diff <- max(abs(fpm_rounded - PMfire_rounded))
mean_diff <- mean(abs(fpm_rounded - PMfire_rounded))
cat("Comparing: our park_classic_2015_fpm (lat-flipped) vs MAT PMfire[Classic, 2015, scenario 1]\n")
cat("Our matrix dimensions:", dim(fpm_matrix_t), "\n")
cat("MAT matrix dimensions:", dim(PMfire_classic_2015), "\n")
cat("Maximum difference:", max_diff, "\n")
cat("Mean difference:", mean_diff, "\n")
cat("\nFirst 5x5 from MAT file (rounded):\n")
print(PMfire_rounded[1:5, 1:5])
cat("\nFirst 5x5 from our data (rounded):\n")
print(fpm_rounded[1:5, 1:5])
if (max_diff == 0) {
  cat("\n SUCCESS: Matrices are identical at 6 decimal precision!\n")
} else {
  cat("\n WARNING: Matrices differ at 6 decimal precision!\n")
}

# Truncate to some decimal places (keep only some digits after decimal)
PMfire_truncated <- trunc(PMfire_classic_2015 * 1e6) / 1e6
fpm_truncated <- trunc(fpm_matrix_t * 1e6) / 1e6

# Compare
max_diff <- max(abs(fpm_truncated - PMfire_truncated))
mean_diff <- mean(abs(fpm_truncated - PMfire_truncated))
cat("Our matrix dimensions:", dim(fpm_matrix_t), "\n")
cat("MAT matrix dimensions:", dim(PMfire_classic_2015), "\n")
cat("Maximum difference:", max_diff, "\n")
cat("Mean difference:", mean_diff, "\n")
cat("\nFirst 5x5 from MAT file (truncated):\n")
print(PMfire_truncated[1:5, 1:5])
cat("\nFirst 5x5 from our data (truncated):\n")
print(fpm_truncated[1:5, 1:5])
if (max_diff == 0) {
  cat("\n SUCCESS: Matrices are identical after truncation!\n")
} else {
  cat("\n Difference remains:", max_diff, "\n")
}

# Find different cells
# the source NetCDF grid centers are at -179.5..180,
# while our target grid centers are at -179.75..179.75 (offset by 0.25 deg),
# so bilinear regridding introduces minor smoothing.
diff_matrix <- fpm_matrix_t - PMfire_classic_2015
diff_indices <- which(abs(diff_matrix) > 0, arr.ind = TRUE)
cat("\nNumber of different cells:", nrow(diff_indices), "\n")
if (nrow(diff_indices) > 0) {
  cat("\nFirst 10 different cells:\n")
  for (i in 1:min(20, nrow(diff_indices))) {
    row <- diff_indices[i, 1]
    col <- diff_indices[i, 2]
    cat(sprintf("Position[%d,%d]: MAT=%.10f, Ours=%.10f, Diff=%.10e\n", 
                row, col, PMfire_classic_2015[row, col], fpm_matrix_t[row, col], diff_matrix[row, col]))
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
######################## save final master data to local ####################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# rename columns for final save

combined_data <- combined_data %>%
  rename_with(
    .cols = starts_with("park_"),  # Target only "park_" columns
    .fn = function(x) {
      for (y in years) {
        # Replace each specific year with its corresponding decade
        x <- str_replace_all(
          x,
          as.character(y),
          paste0(floor(y / 10) * 10, "s")
        )
      }
      x
    }
  )

write.csv(combined_data, here("output", "pop_pm_combined_with_park2024.csv"), row.names = FALSE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ SANITY CHECK: Plot US fpm for ALL models x years ##############
############ Output as HTML to local folder                       ##########
############ does not affect final dataset           #######################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(ncdf4)
library(ggplot2)
library(dplyr)
library(maps)
library(patchwork)
library(htmltools)
library(base64enc)

# --- Output settings ---
out_dir <- here::here("images", "fpm_cell_usa")
out_name <- "fpm_cell_usa_park_sanitycheck"
out_html <- file.path(out_dir, paste0(out_name, ".html"))

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# --- helper: read PM25 surface field with original lon/lat ---
read_pm25 <- function(path) {
  nc <- nc_open(path)
  lon_name <- if ("lon" %in% names(nc$dim)) "lon" else "longitude"
  lat_name <- if ("lat" %in% names(nc$dim)) "lat" else "latitude"
  lon <- as.numeric(ncvar_get(nc, lon_name))
  lat <- as.numeric(ncvar_get(nc, lat_name))
  
  # ====== Print ORIGINAL lon/lat from the NetCDF file ======
  cat("\n----- ORIGINAL coordinates from file -----\n")
  cat("File: ", basename(path), "\n", sep = "")
  cat("Variable names:  lon = '", lon_name, "' | lat = '", lat_name, "'\n", sep = "")
  cat(sprintf("lon: n = %d | range = [%.4f, %.4f] | first 5 = %s | last 5 = %s\n",
              length(lon),
              min(lon), max(lon),
              paste(round(head(lon, 5), 4), collapse = ", "),
              paste(round(tail(lon, 5), 4), collapse = ", ")))
  cat(sprintf("lat: n = %d | range = [%.4f, %.4f] | first 5 = %s | last 5 = %s\n",
              length(lat),
              min(lat), max(lat),
              paste(round(head(lat, 5), 4), collapse = ", "),
              paste(round(tail(lat, 5), 4), collapse = ", ")))
  cat(sprintf("lon resolution = %.4f | lat resolution = %.4f\n",
              mean(diff(lon)), mean(diff(lat))))
  cat("------------------------------------------\n")
  # =========================================================
  
  pm <- ncvar_get(nc, "PM25")
  d  <- dim(pm)
  pm <- switch(as.character(length(d)),
               "4" = pm[, , 1, 1],
               "3" = pm[, , 1],
               "2" = pm,
               stop("Unexpected PM25 dimensions"))
  nc_close(nc)
  
  if (nrow(pm) == length(lat) && ncol(pm) == length(lon)) {
    pm <- t(pm)
  }
  
  # Ensure latitude is ascending (south -> north) for correct raster rendering
  if (is.unsorted(lat)) {
    ord <- order(lat)
    lat <- lat[ord]
    pm  <- pm[, ord, drop = FALSE]
  }
  
  # Ensure longitude is ascending too
  if (is.unsorted(lon)) {
    ord <- order(lon)
    lon <- lon[ord]
    pm  <- pm[ord, , drop = FALSE]
  }
  
  list(lon = lon, lat = lat, pm = pm)
}

# --- build fpm dataframe for one model x one year ---
build_fpm_df <- function(path_model, path_withoutfire) {
  a <- read_pm25(path_model)
  b <- read_pm25(path_withoutfire)
  stopifnot(length(a$lon) == length(b$lon),
            length(a$lat) == length(b$lat))
  
  fpm_mat <- a$pm - b$pm
  
  df <- expand.grid(lon = a$lon, lat = a$lat,
                    KEEP.OUT.ATTRS = FALSE) %>%
    mutate(fpm = as.vector(fpm_mat))
  
  if (max(df$lon) > 180) {
    df <- df %>% mutate(lon = ifelse(lon > 180, lon - 360, lon))
  }
  
  df %>% filter(lon >= -125, lon <= -66, lat >= 24, lat <= 50)
}

# --- helper: single plot ---
us_map <- map_data("state")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Color Palette: White -> Yellow -> Dark Coffee #################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Low fPM = white, High fPM = dark coffee
fpm_palette <- c("#3B1F0E", "#6B3E1B", "#A9651A", "#D9A441",
                 "#F2D04A", "#FFF2A8", "#FFFFFF")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ make_plot: matching plot_fpm_var styling ######################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

make_plot <- function(df, title, fill_max) {
  ggplot(df, aes(x = lon, y = lat, fill = fpm)) +
    geom_tile(width = 0.5, height = 0.5,
              color = "grey50", linewidth = 0.05) +     
    geom_polygon(data = us_map,
                 aes(x = long, y = lat, group = group),
                 fill = NA, color = "grey30", linewidth = 0.15,
                 inherit.aes = FALSE) +
    scale_fill_gradientn(
      colors   = rev(fpm_palette),
      name     = "fPM",
      limits   = c(0, fill_max),
      oob      = scales::squish,
      na.value = "white"
    ) +
    coord_quickmap(xlim = c(-125.5, -66.0), ylim = c(23.5, 50.0),
                   expand = FALSE) +
    labs(title = title) +
    theme_void(base_size = 10) +
    theme(
      plot.title        = element_text(hjust = 0.5, size = 10, face = "bold"),
      legend.key.height = unit(0.45, "cm"),
      legend.key.width  = unit(0.3, "cm"),
      legend.title      = element_text(size = 7),
      legend.text       = element_text(size = 6),
      plot.margin       = margin(3, 3, 3, 3)
    )
}

# --- build one combined image per model, save to temp PNG, then embed ---
models <- list(
  classic = file_paths_classic,
  jules   = file_paths_jules,
  ssib4   = file_paths_ssib4
)

png_to_base64 <- function(png_path) {
  paste0("data:image/png;base64,",
         base64enc::base64encode(png_path))
}

img_tags <- list()

for (m in names(models)) {
  
  cat("Building plots for model:", m, "\n")
  
  dfs <- lapply(seq_along(years), function(i) {
    build_fpm_df(models[[m]][i], file_paths_withoutfire[i])
  })
  
  fill_max <- quantile(unlist(lapply(dfs, `[[`, "fpm")),
                       0.99, na.rm = TRUE)
  
  year_imgs <- lapply(seq_along(years), function(i) {
    p <- make_plot(dfs[[i]],
                   title = paste0(m, " — ", years[i]),
                   fill_max = fill_max)
    tmp <- tempfile(fileext = ".png")
    ggsave(tmp, p, width = 9, height = 5.5, dpi = 220, bg = "white")
    
    tags$img(
      src   = png_to_base64(tmp),
      style = paste("width: 49%; margin: 0.5%; vertical-align: top;",
                    "border: 1px solid #ddd;",
                    "image-rendering: pixelated;",
                    "image-rendering: crisp-edges;")
    )
  })
  
  img_tags[[m]] <- tags$div(
    style = "margin-bottom: 40px;",
    tags$h2(paste0("Model: ", m),
            style = "font-family: sans-serif; color: #333;"),
    do.call(tags$div, year_imgs)
  )
}

# --- assemble HTML ---
page <- tags$html(
  tags$head(
    tags$title("Park et al. 2024 — Fire PM2.5 Sanity Check (CONUS)"),
    tags$meta(charset = "utf-8")
  ),
  tags$body(
    style = "font-family: sans-serif; max-width: 1400px; margin: 20px auto; padding: 0 20px;",
    tags$h1("Fire PM2.5 Sanity Check — Park et al. 2024"),
    tags$p(
      style = "color: #555;",
      paste0("Generated on ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
             ". Plots use original NetCDF coordinates (no regrid). ",
             "Color scale per model: 0 to 99th percentile of all years.")
    ),
    tags$hr(),
    img_tags
  )
)

save_html(page, out_html)

cat("\n HTML saved to:\n", out_html, "\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Sanity Check: Native Model Resolution ##########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Purpose: Verify the true (native) resolution of the GEOS-Chem PM2.5 fields
# distributed by Park et al. 2024. The files are stored on a 0.5° x 0.5° mesh,
# but inspection of block structure (consecutive identical values) reveals the
# underlying simulation grid is much coarser.
#
# Reference: GEOS-Chem Classic 14.7.1 horizontal grids documentation
#   https://geos-chem.readthedocs.io/en/latest/supplemental-guides/horizontal-grids.html
#
# Expected for a 4° x 5° global grid with half_size_polar_boxes = TRUE:
#   - Mid-latitude rows: 4° lat / 0.5° = 8 consecutive identical 0.5° cells
#   - Polar half-cells (lat = +/-89, width 2°): 4 consecutive identical cells
#   - Total lat runs: 44 mid-latitude + 2 polar = 46 runs spanning 360 cells

diagnose_resolution <- function(file_path, label, sample_lon = -100, sample_lat = 40) {
  cat("\n============================================================\n")
  cat("DIAGNOSTIC FOR:", label, "\n")
  cat("File:", basename(file_path), "\n")
  cat("============================================================\n")
  
  nc <- nc_open(file_path)
  lon_src <- as.numeric(ncvar_get(nc, "lon"))
  lat_src <- as.numeric(ncvar_get(nc, "lat"))
  pm      <- ncvar_get(nc, "PM25")
  nc_close(nc)
  
  # [1] Coordinate metadata
  cat("\n[1] Shape & coordinate metadata\n")
  cat(sprintf("    length(lon) = %d | range = [%g, %g] | res = %g\n",
              length(lon_src), min(lon_src), max(lon_src), mean(diff(lon_src))))
  cat(sprintf("    length(lat) = %d | range = [%g, %g] | res = %g\n",
              length(lat_src), min(lat_src), max(lat_src), mean(diff(lat_src))))
  cat(sprintf("    dim(pm)     = %d x %d\n", dim(pm)[1], dim(pm)[2]))
  
  # [2] Unique values along each axis at a sample point
  lon_idx <- which.min(abs(lon_src - sample_lon))
  lat_idx <- which.min(abs(lat_src - sample_lat))
  cat("\n[2] Unique values along each axis (at sample point)\n")
  cat(sprintf("    Sample column (lon = %g): %d unique values out of %d\n",
              lon_src[lon_idx], length(unique(pm[lon_idx, ])), length(lat_src)))
  cat(sprintf("    Sample row    (lat = %g): %d unique values out of %d\n",
              lat_src[lat_idx], length(unique(pm[, lat_idx])), length(lon_src)))
  
  # [3] Block-size structure (run-length encoding of identical values)
  # In a remapped coarse-to-fine field, identical 0.5° cells form blocks whose
  # length reveals the native cell width. Mode of run lengths = typical block.
  lat_runs <- rle(pm[lon_idx, ])
  lon_runs <- rle(pm[, lat_idx])
  most_common <- function(x) as.integer(names(sort(table(x), decreasing = TRUE)[1]))
  cat("\n[3] Detected block size (consecutive cells with identical value)\n")
  cat(sprintf("    Along lat axis: most common run length = %d   (samples: %s)\n",
              most_common(lat_runs$lengths),
              paste(head(lat_runs$lengths, 10), collapse = ",")))
  cat(sprintf("    Along lon axis: most common run length = %d   (samples: %s)\n",
              most_common(lon_runs$lengths),
              paste(head(lon_runs$lengths, 10), collapse = ",")))
  
  # [4] Inferred true resolution from block sizes
  eff_lat_res <- most_common(lat_runs$lengths) * abs(mean(diff(lat_src)))
  eff_lon_res <- most_common(lon_runs$lengths) * abs(mean(diff(lon_src)))
  cat("\n[4] Inferred TRUE resolution\n")
  cat(sprintf("    Stored lon res = %g deg | effective lon res ~ %g deg\n",
              abs(mean(diff(lon_src))), eff_lon_res))
  cat(sprintf("    Stored lat res = %g deg | effective lat res ~ %g deg\n",
              abs(mean(diff(lat_src))), eff_lat_res))
  
  # [5] 5x5 patch near sample point (visual confirmation of block structure)
  cat("\n[5] 5x5 PM25 patch near (lon =", sample_lon, ", lat =", sample_lat, ")\n")
  print(round(pm[lon_idx:(lon_idx + 4), lat_idx:(lat_idx + 4)], 4))
  
  # [6] Polar-cap check: half_size_polar_boxes leaves shorter runs at poles.
  # For a 4° x 5° grid, mid-latitude runs = 8 cells, polar runs = 4 cells.
  full_lat_runs <- rle(pm[lon_idx, ])
  cat("\n[6] Polar-cap structure (full lat axis run lengths)\n")
  cat(sprintf("    Number of lat runs: %d (expect 46 for 4 deg x 5 deg global grid)\n",
              length(full_lat_runs$lengths)))
  cat(sprintf("    Run length at far south: %d (expect 4 for half polar cell)\n",
              full_lat_runs$lengths[1]))
  cat(sprintf("    Run length at far north: %d (expect 4 for half polar cell)\n",
              full_lat_runs$lengths[length(full_lat_runs$lengths)]))
  
  # [7] Verdict
  cat("\n[7] Verdict\n")
  cat("    Data is stored on a 0.5 deg mesh BUT effective resolution is COARSER.\n")
  cat(sprintf("    Native model grid is approximately %g deg lon x %g deg lat.\n",
              eff_lon_res, eff_lat_res))
  cat("    Consistent with GEOS-Chem Classic 4 deg x 5 deg global grid\n")
  cat("    with half_size_polar_boxes = TRUE (mid-lat runs = 8, polar runs = 4,\n")
  cat("    46 total lat runs spanning 360 cells = 44 mid-lat + 2 polar).\n")
  cat("    The lon axis appears continuous because conservative regridding\n")
  cat("    (CDO remapcon) blends adjacent 5 deg native cells into smoothly\n")
  cat("    varying 0.5 deg sub-cells along the longitude direction.\n")
  cat("============================================================\n")
}

# Run diagnostics on representative files (one per scenario, plus a year check)
diagnose_resolution(file_paths_classic[6],     "classic 2015")
diagnose_resolution(file_paths_classic[1],     "classic 1965")
diagnose_resolution(file_paths_withoutfire[6], "withoutfire 2015")
diagnose_resolution(file_paths_jules[6],       "jules 2015")
diagnose_resolution(file_paths_ssib4[6],       "ssib4 2015")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Conclusion (consistent with GEOS-Chem Classic horizontal grids documentation):
#   - Native simulation grid: 4° lat × 5° lon, with half_size_polar_boxes = TRUE.
#   - Data are stored on a finer 0.5° × 0.5° grid after conservative remapping
#     (CDO remapcon), but this does NOT reflect the true model resolution.
#   - Latitude direction clearly shows the native resolution (~4°, ~2° at poles)
#     through repeated identical values (block structure).
#   - Longitude direction appears smooth because remapcon spreads each 5° grid
#     cell across multiple 0.5° cells, masking the original coarse resolution.
#   - Therefore, the true spatial resolution remains 4° × 5°, and any finer-scale
#     variation (e.g., <4°) should be interpreted as a regridding artifact,
#     not a physical signal.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#THE END

