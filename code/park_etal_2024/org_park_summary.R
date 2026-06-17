# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
########### PM2.5 yearly statistics from NetCDF ###########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Goal: For each land-surface scheme, summarize global PM2.5 levels across six snapshot period.
#   1. Build file paths for three GEOS-Chem scenarios (classic / jules / ssib4),
#      each with one yearly-average NetCDF per year (1965–2015, every 10 yr).
#   2. Read the 'PM25' grid from each file, flattening whatever dims it has
#      (4-D, 3-D, or 2-D) down to a single vector of cell values.
#   3. Per scenario, print a year-by-year table of min / max / mean / median PM2.5.
# Inputs : Park et al. 2024, GEOS-Chem CEDS "on_off" surface PM2.5, 0.5° yearly avg (.nc4)

rm(list = ls())

library(here)
library(tidyverse)
library(ncdf4)

options(pillar.sigfig = 7, digits = 7)

years <- c(1965, 1975, 1985, 1995, 2005, 2015)

scenario_paths <- list(
  classic = vapply(years, function(y) here("input","Park_etal_2024","GEOSChem_output","classic","obsclim",
                                           sprintf("05x05_CEDS_%d_on_off_pm25_Surface_Re_yearavg.nc4", y)), character(1)),
  jules   = vapply(years, function(y) here("input","Park_etal_2024","GEOSChem_output","jules","obsclim",
                                           sprintf("05x05_CEDS_%d_on_off_pm25_Surface_Re_yearavg.nc4", y)), character(1)),
  ssib4   = vapply(years, function(y) here("input","Park_etal_2024","GEOSChem_output","ssib4","obsclim",
                                           sprintf("05x05_CEDS_%d_on_off_pm25_Surface_Re_yearavg.nc4", y)), character(1))
)

read_pm25_vec <- function(path) {
  nc <- nc_open(path)
  pm <- ncvar_get(nc, "PM25")
  nc_close(nc)
  d <- dim(pm)
  pm <- switch(as.character(length(d)),
               "4" = pm[, , 1, 1],
               "3" = pm[, , 1],
               "2" = pm)
  as.vector(pm)
}

for (sc in names(scenario_paths)) {
  cat("\n=== Scenario:", sc, "===\n")
  map2_dfr(scenario_paths[[sc]], years, function(p, y) {
    v <- read_pm25_vec(p)
    tibble(
      year   = y,
      min    = min   (v, na.rm = TRUE),
      max    = max   (v, na.rm = TRUE),
      mean   = mean  (v, na.rm = TRUE),
      median = median(v, na.rm = TRUE)
    )
  }) %>% print()
}