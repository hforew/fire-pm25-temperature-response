# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
########### Confirm if PM data grid cell centers or edges ##########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Remove all objects from the environment
rm(list = ls())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Packages #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(ncdf4)
library(here)
library(tidyverse)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
################################# read in data #################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pm_2000 <- nc_open(here("input", "CESM_09x125_PM25_2000_Baseline.nc"))
print(pm_2000) # PM25 NetCDF files for other years follow identical format

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Validate POP and PM2.5 Grid Coordinates: Centers vs. Edges #################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# NetCDF file data was previously extracted and converted to dataframe format. 
# check here if the lat and lon coordinates in netCDF represent grid cell center or edges


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

