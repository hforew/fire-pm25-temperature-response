# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Compute Global Mean Temp Change #####################################################
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
list.files(here("input/landmask_area"))


temp_45 <- nc_open(here("input/temperature", "af.tas.ccsm4.rcp45.2006-2300.nc"))
print(temp_45)

area_09x125 <- nc_open(here("input/landmask_area", "cesm130_clm5_firemodule_area_f09x125.nc"))
print(area_09x125)




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


