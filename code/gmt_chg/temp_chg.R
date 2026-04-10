# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Global Mean Temperature Change — Area-Weighted, Three Periods
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Goal: Compute area-weighted GMT change for 2006-2010, 2041-2050, 2091-2100
#       using pre-computed grid cell areas from CESM area file
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# This code computes global mean temp change (GMTΔ) for RCP4.5 and RCP8.5.  

# Remove all objects from the environment
rm(list = ls())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Packages #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(ncdf4)
library(here)
library(tidyverse)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ File import #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

list.files(here("input/temperature"))


temp_45 <- nc_open(here("input/temperature", "af.tas.ccsm4.rcp45.2006-2300.nc"))
print(temp_45)

# Primary Data Variable temp_45
# ~~~~~~~~~~~~~
# tas: Reference height temperature anomaly
#   - Dimensions: 288 lon × 192 lat × 3540 time steps
#   - Units: "none" (anomaly = deviation from baseline, not absolute temperature)
#   - Time-dependent variable with monthly data
#   - Temporal span: ~295 years (3540 months / 12 = 295 years, 2006-2300)
#   - Fill value: 9.999e+35 for missing data
#   - Storage: 3D array [lon, lat, time]

# Spatial Coverage
# ~~~~~~~~~~~~~
# Grid structure: 288 × 192 grid cells
#   - Resolution: ~1.25° longitude × ~0.94° latitude (standard CCSM4 resolution)
#   - LONGXY/LATIXY: 2D arrays providing lon/lat coordinates for each grid cell
#   - Geographic boundaries defined by edge variables:
#     * EDGEE: eastern edge in atmospheric data (degrees_east)
#     * EDGEW: western edge in atmospheric data (degrees_east)
#     * EDGES: southern edge in atmospheric data (degrees_north)
#     * EDGEN: northern edge in atmospheric data (degrees_north)
#   - Likely global coverage given standard climate model grid dimensions

# Temporal Coverage
# ~~~~~~~~~~~~~
# Start date: January 1, 2006 (00:00:00)
# End date: December 2300 (2006 + 295 years - 1 = 2300)
# Calendar: "noleap" (365 days per year, no leap years in time calculations)
# Time steps: 3540 monthly observations
# Duration: 295 years of monthly data
# Time units: "days since 2006-01-01 00:00:00"
# Temporal resolution: Monthly (assumed from 3540 steps / 295 years = 12/year)
# Scenario: RCP4.5 (intermediate emissions pathway)
# Model: CCSM4 (Community Climate System Model version 4)


temp_85 <- nc_open(here("input/temperature", "af.tas.ccsm4.rcp85.2006-2300.nc"))
print(temp_85)

# temp_85 has identical file structure to temp_45, for RCP8.5.  

# Area will be computed directly from coordinates using Method 3 (WGS84) — see below




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ confirm if grid cell edges / centers #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Temperature file uses LATIXY/LONGXY as 2D coordinate variables
lat_temp <- ncvar_get(temp_45, "LATIXY")[1, ]   # one column suffices — lat varies along dim 2
lon_temp <- ncvar_get(temp_45, "LONGXY")[, 1]   # one row suffices — lon varies along dim 1

# Check first and last few values
cat("Temp file lat:", head(lat_temp, 3), "...", tail(lat_temp, 3), "\n")
cat("Temp file lon:", head(lon_temp, 3), "...", tail(lon_temp, 3), "\n")

# --- coordinates start at 0 / -90 and step by ~1.25° / ~0.94° — these are cell centers


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Extract data and close files ####################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Extract temperature anomaly arrays for both scenarios [lon, lat, time]
tas_45 <- ncvar_get(temp_45, "tas")    # RCP4.5 — dim: [288, 192, 3540]
tas_85 <- ncvar_get(temp_85, "tas")    # RCP8.5 — dim: [288, 192, 3540]

# Extract time vector — shared across both scenarios (same temporal structure)
# Units: days since 2006-01-01, noleap calendar
time_days <- ncvar_get(temp_45, "time")   # length: 3540

# Extract 2D coordinate arrays for area computation
# LATIXY/LONGXY are [lon, lat] arrays giving cell center coordinates
lat_2d <- ncvar_get(temp_45, "LATIXY")   # dim: [288, 192]
lon_2d <- ncvar_get(temp_45, "LONGXY")   # dim: [288, 192]

# Close all files — all needed data now in memory
nc_close(temp_45)
nc_close(temp_85)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Mask fill values ################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Replace fill values with NA so they are excluded from weighted mean
fill_val <- 9.99999961690316e+35

tas_45[tas_45 >= fill_val] <- NA
tas_85[tas_85 >= fill_val] <- NA

# Fraction of NAs after masking — confirms no unexpected data loss
cat("tas_45 NA fraction:", mean(is.na(tas_45)), "\n")
cat("tas_85 NA fraction:", mean(is.na(tas_85)), "\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Compute grid cell areas — Method 2 (spherical Earth) ############
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Method : area = R² × dlon_rad × (sin(lat_upper_rad) - sin(lat_lower_rad))
# Assumes spherical Earth with R = 6371 km. 
#
# Grid spacing in latitude is not perfectly uniform on the f09x125 grid,
# so dlat is derived per-row from neighbouring cell centers.
# dlon is fixed at 1.25° across all columns.

R    <- 6371       # Earth radius in km
dlon <- 1.25       # fixed longitude spacing (degrees)

# Extract 1D coordinate vectors from the 2D arrays
# (lat varies along dim 2, lon varies along dim 1)
lat_vec <- lat_2d[1, ]   # 192 latitude values — one per row
lon_vec <- lon_2d[, 1]   # 288 longitude values — one per column

# Goal: compute the latitude "thickness" of each grid cell (dlat) in degrees.
# Because the f09x125 grid spacing is not perfectly uniform in latitude,
# we cannot assume a fixed dlat — it must be derived from the actual
# coordinate values for each row.

# Pre-allocate a vector of zeros, one entry per latitude row (192 total)
dlat_vec <- numeric(length(lat_vec))

# --- South pole cell (index 1) ---
# Only has one neighbour (the cell above it), so dlat is simply the
# distance to that single neighbour
dlat_vec[1] <- abs(lat_vec[2] - lat_vec[1])

# --- North pole cell (last index) ---
# Only has one neighbour (the cell below it), same logic as south pole
dlat_vec[length(lat_vec)] <- abs(lat_vec[length(lat_vec)] - lat_vec[length(lat_vec) - 1])

# --- All interior cells (indices 2 to 191) ---
# Each interior cell has a neighbour above and below.
# The cell boundary sits halfway between cell centers, so:
#   - southern half-width = distance to cell below / 2
#   - northern half-width = distance to cell above / 2
#   - total dlat = (gap to cell below + gap to cell above) / 2
#
# diff(lat_vec) produces a vector of 191 consecutive differences:
#   diff(lat_vec)[k] = lat_vec[k+1] - lat_vec[k]
#
# For interior cell j (1-indexed):
#   gap to cell below = diff(lat_vec)[j-1]   → index [1:(n-2)] covers j = 2..191
#   gap to cell above = diff(lat_vec)[j]     → index [2:(n-1)] covers j = 2..191
dlat_vec[2:(length(lat_vec)-1)] <- (abs(diff(lat_vec)[1:(length(lat_vec)-2)]) +
                                      abs(diff(lat_vec)[2:(length(lat_vec)-1)])) / 2
# Compute lat bounds for each row (in radians)
lat_lower_rad <- (lat_vec - dlat_vec / 2) * pi / 180   # southern edge of each cell
lat_upper_rad <- (lat_vec + dlat_vec / 2) * pi / 180   # northern edge of each cell

# Longitude spacing in radians — same for all cells
dlon_rad <- dlon * pi / 180

# Cell area for each latitude row (km²): vectorised over 192 lat values
# Area depends only on latitude, not longitude — replicate across all 288 lon columns
area_by_lat <- R^2 * dlon_rad * (sin(lat_upper_rad) - sin(lat_lower_rad))

# Expand to full [lon, lat] matrix by repeating area_by_lat across all lon rows
# matrix() fills column-wise, so each column (= one lat) gets the same area value
grid_area <- matrix(
  rep(area_by_lat, each = length(lon_vec)),   # repeat each lat area 288 times
  nrow = length(lon_vec),                      # 288 rows (lon)
  ncol = length(lat_vec)                       # 192 cols (lat)
)

cat("Area range (km²):", round(min(grid_area), 2), "to", round(max(grid_area), 2), "\n")

# Sum all cell areas to get total global area in km²
total_area_km2 <- sum(grid_area)
cat("Total global area (km²):", round(total_area_km2, 0), "\n")

# Earth's true surface area is ~510,072,000 km² — compare as sanity check
cat("Expected (km²):           510,072,000\n")
cat("Difference (%):           ", round((total_area_km2 - 510072000) / 510072000 * 100, 3), "\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Convert time to years ###########################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Noleap calendar: 365 days/year exactly
# floor() assigns each month to its year; +2006 gives calendar year
time_years <- floor(time_days / 365) + 2006

cat("Year range in data:", min(time_years), "to", max(time_years), "\n")
cat("Total months:", length(time_years), "\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Compute annual mean at each grid cell ###########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

unique_years <- unique(time_years)   # 295 years: 2006-2300

# Goal: collapse monthly temperature data into annual means at each grid cell.
# tas_45 and tas_85 are currently [288 lon × 192 lat × 3540 months].
# We want [288 lon × 192 lat × 295 years] — one value per cell per year.

# Pre-allocate 3D arrays to hold annual means, filled with NA initially.
# dim = [288 lon, 192 lat, 295 years] — same spatial dims, time collapsed to years.
# Pre-allocating is faster than growing the array inside the loop.
tas_45_annual <- array(NA, dim = c(288, 192, length(unique_years)))
tas_85_annual <- array(NA, dim = c(288, 192, length(unique_years)))

# Loop over each of the 295 unique years (2006, 2007, ..., 2300)
for (i in seq_along(unique_years)) {
  
  # Get the calendar year corresponding to loop index i
  yr <- unique_years[i]
  
  # Find which of the 3540 time indices belong to this year.
  # time_years is a length-3540 vector where each entry is the year of that month.
  # which() returns the positions (e.g. 1:12 for 2006, 13:24 for 2007, etc.)
  month_idx <- which(time_years == yr)   # integer vector of length 12
  
  # tas_45[, , month_idx] subsets the full array to just the 12 monthly slices
  # for this year, giving a temporary [288, 192, 12] array.
  # apply(..., c(1, 2), mean) then collapses dimension 3 (the 12 months) by
  # taking the mean at each [lon, lat] position — result is [288, 192].
  # That 2D annual mean is assigned into the i-th year slot of the output array.
  tas_45_annual[, , i] <- apply(tas_45[, , month_idx], c(1, 2), mean, na.rm = TRUE)
  tas_85_annual[, , i] <- apply(tas_85[, , month_idx], c(1, 2), mean, na.rm = TRUE)
  
}

cat("Annual array dimensions:", dim(tas_45_annual), "\n")   # 288 192 295

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Compute area-weighted global mean for each year #################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Pre-allocate vectors for annual GMT values
gmt_45 <- numeric(length(unique_years))
gmt_85 <- numeric(length(unique_years))

for (i in seq_along(unique_years)) {
  
  # Extract annual temperature slice [lon, lat] for each scenario
  slice_45 <- tas_45_annual[, , i]
  slice_85 <- tas_85_annual[, , i]
  
  # tas has zero NAs — grid_area covers all cells — no masking needed
  # Area-weighted mean: sum(tas * area) / sum(area)
  gmt_45[i] <- sum(slice_45 * grid_area, na.rm = TRUE) / sum(grid_area)
  gmt_85[i] <- sum(slice_85 * grid_area, na.rm = TRUE) / sum(grid_area)
  
}

# Combine into tidy dataframe
gmt_df <- tibble(
  year   = unique_years,
  gmt_45 = gmt_45,
  gmt_85 = gmt_85
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Average GMT within each target period ###########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Define three target periods
periods <- list(
  "2006-2010" = 2006:2010,
  "2041-2050" = 2041:2050,
  "2091-2100" = 2091:2100
)

# Goal: for each of the 3 time periods, average the annual GMT values
# across all years in that period, for both RCP scenarios.
# Result: a 3-row dataframe, one row per period.

# map_dfr() iterates over a vector and row-binds the results into one dataframe.
# Here it iterates over the 3 period names: "2006-2010", "2041-2050", "2091-2100".
# For each period name, the anonymous function below is called once.
gmt_periods <- map_dfr(names(periods), function(period_name) {
  
  # Look up the integer vector of years for this period from the periods list.
  # e.g. for "2041-2050", yrs = c(2041, 2042, ..., 2050)
  yrs <- periods[[period_name]]
  
  gmt_df |>
    # Keep only rows where the year column falls within this period's year range.
    # %in% checks membership — equivalent to year == 2041 | year == 2042 | ...
    filter(year %in% yrs) |>
    
    # Collapse the filtered rows (e.g. 10 rows for a 10-year period) into
    # a single summary row by averaging gmt_45 and gmt_85 across those years.
    # period_name is carried through as a label column so we know which
    # period each row in the final output corresponds to.
    summarise(
      period      = period_name,
      mean_gmt_45 = mean(gmt_45, na.rm = TRUE),   # avg annual GMT across period, RCP4.5
      mean_gmt_85 = mean(gmt_85, na.rm = TRUE)    # avg annual GMT across period, RCP8.5
    )
  
})
# map_dfr() stacks the 3 single-row results into a 3-row × 3-column dataframe

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Print results ###################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# gmt_periods values are anomalies relative to 2005 baseline (per paper author)
cat("\nGlobal Mean Temperature Change by Period (relative to 2005):\n")
print(gmt_periods)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Adjust to pre-industrial baseline ###############################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Pre-industrial (PI) to 2005 warming = +0.91°C approximately  https://ourworldindata.org/grapher/temperature-anomaly
# Adding 0.91 to each anomaly shifts the reference from 2005 to pre-industrial
pi_offset <- 0.914 # 0.914 is average change 2001-2010

gmt_periods_pi <- gmt_periods |>
  mutate(
    mean_gmt_45 = mean_gmt_45 + pi_offset,   # RCP4.5 anomaly relative to PI
    mean_gmt_85 = mean_gmt_85 + pi_offset    # RCP8.5 anomaly relative to PI
  )

cat("\nGlobal Mean Temperature Change by Period (relative to pre-industrial):\n")
print(gmt_periods_pi)

# save output
write_csv(gmt_periods_pi, here("output", "gmt_periods_pi.csv"))


# THE END 