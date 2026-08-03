# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
########### Inspect a Single GEOS-Chem PM2.5 NetCDF (classic 2015) ##########
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read one decade-mean surface PM2.5 file from Park et al. (2024) and inspect
# its native grid before any regridding is applied.
#
# Input  : 05x05_CEDS_2015_on_off_pm25_Surface_Re_yearavg.nc4
#          (classic / obsclim scenario, 2015 decade mean)
# Output : none written to disk; objects pm_df and pm_mat left in the session.
#
# Steps  : (1) Read lon/lat from the file's own dimension values (never
#              reconstructed) and read the PM25 field.
#          (2) Diagnose whether stored coordinates are cell CENTERS or EDGES.
#          (3) Reduce PM25 to a 2D [lon, lat] surface field.
#          (4) Build two views of the same data: a long data frame and a
#              matrix carrying the file's coordinates as dimnames.
# Execution order:
#   files run before: none (standalone inspection)
#   files run after:  pm_regrid_park_2024.R

# Remove all objects from the environment

rm(list = ls())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Packages ######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(here)
library(tidyverse)
library(ncdf4)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ File Path #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

f <- here("input", "Park_etal_2024", "GEOSChem_output", "classic", "obsclim",
          "05x05_CEDS_2015_on_off_pm25_Surface_Re_yearavg.nc4")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Read Coordinates and PM2.5 ####################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# lon/lat are taken from the file's own dimension values. This file reports
# nvars = 2, so lon/lat are dimension variables rather than data variables.

nc  <- nc_open(f)
lon <- as.numeric(nc$dim$lon$vals)
lat <- as.numeric(nc$dim$lat$vals)
pm  <- ncvar_get(nc, "PM25")
nc_close(nc)

cat("dim(PM25):", dim(pm), "\n")
cat("lon: n =", length(lon), "| range =", range(lon), "\n")
cat("lat: n =", length(lat), "| range =", range(lat), "\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Check: Cell Centers or Cell Edges? #############################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# For a global axis of n points with spacing r:
#   CENTERS: values sit at half-offsets of r (e.g. -179.75, -179.25, ...)
#   EDGES:   values sit on integer multiples of r (e.g. -179.5, -179.0, ...)
# n * r == full span holds in both cases, so the offset pattern is what
# distinguishes them.

check_axis <- function(x, name) {
  r    <- mean(diff(x))
  frac <- (x / r) %% 1
  frac <- pmin(frac, 1 - frac)          # 0 = on a multiple, 0.5 = halfway
  
  verdict <- if (max(abs(frac - 0.5)) < 1e-6) "CENTERS" else
    if (max(abs(frac))       < 1e-6) "EDGES"   else "AMBIGUOUS"
  
  cat(sprintf("[%s] n = %d | res = %g | range = [%g, %g] -> %s\n",
              name, length(x), r, min(x), max(x), verdict))
  cat("     first 5:", head(x, 5), "| last 5:", tail(x, 5), "\n")
}

cat("\n=== Coordinate convention ===\n")
check_axis(lon, "lon")
check_axis(lat, "lat")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Reduce PM2.5 to a 2D [lon, lat] Surface Field ##################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pm <- switch(as.character(length(dim(pm))),
             "4" = pm[, , 1, 1],       # surface level, first time step
             "3" = pm[, , 1],
             "2" = pm,
             stop("Unexpected PM25 dimensions"))

if (nrow(pm) == length(lat) && ncol(pm) == length(lon)) pm <- t(pm)
stopifnot(nrow(pm) == length(lon), ncol(pm) == length(lat))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Output 1: Long Format Data Frame ##############################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# One row per grid cell; lon varies first, matching as.vector(pm).

pm_df <- expand.grid(lon = lon, lat = lat,
                     KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE) %>%
  mutate(pm25 = as.vector(pm)) %>%
  as_tibble()

print(pm_df)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Output 2: Matrix with File Coordinates as Dimnames ############
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pm_mat <- pm
dimnames(pm_mat) <- list(lon = lon, lat = lat)

cat("\nMatrix dimensions:", dim(pm_mat), "\n")
cat("Top-left corner (first 5 lon x first 5 lat):\n")
print(round(pm_mat[1:5, 1:5], 4))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Notes:
#   - pm_df and pm_mat hold identical values in two layouts; both carry the
#     file's native coordinates, with nothing reconstructed.
#   - The center-vs-edge verdict determines whether a 0.25 deg offset is
#     introduced when regridding onto the -179.75..179.75 target grid used
#     downstream in pm_regrid_park_2024.R.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# THE END