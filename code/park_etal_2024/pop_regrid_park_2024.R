# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Read and Process Population Grid Data #########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Remove all objects from the environment
rm(list = ls())
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Packages #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(R.matlab)
library(here)
library(tidyr)
library(dplyr)
library(terra)
library(sf)
library(rnaturalearth)

# Read the .mat file using here()
pop_data <- readMat(here("input", "Park_etal_2024", "pop_grid_6515.mat"))

# Extract the population array (360 lat x 720 lon x 6 years)
pop_array <- pop_data$pop.grid

# Define years (same as PM2.5 data)
years <- c(1965, 1975, 1985, 1995, 2005, 2015)

# Create latitude and longitude coordinates
# Source grid: matches PM2.5 data convention (lon -179.5 to 180, lat N-to-S)
lat <- seq(89.75, -89.75, by = -0.5)   # 360 values, N to S
lon <- seq(-179.5, 180.0, by = 0.5)    # 720 values

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Build Source Grid (long format) ################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create list to store data for each year
pop_list <- list()

for (i in 1:6) {
  # Extract population data for year i
  pop_matrix <- pop_array[, , i]
  
  # Create data frame with expand.grid (automatically creates lon-lat grid)
  # expand.grid creates combinations with first variable changing fastest
  df <- expand.grid(
    lon = lon,
    lat = lat,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  # Transpose and vectorize population matrix to match the grid order
  # pop_matrix is [lat, lon], we need to convert to vector in [lon, lat] order
  df$population <- as.vector(t(pop_matrix))
  
  # Add year
  df$year <- years[i]
  
  # Reorder columns: year, lon, lat, population
  df <- df %>%
    select(year, lon, lat, population)
  
  pop_list[[i]] <- df
}

# Combine all years into long format
pop_long <- bind_rows(pop_list)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Regrid to Target Grid (-179.75 to 179.75) ######################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Target grid: lon -179.75 to 179.75, lat -89.75 to 89.75 (step 0.5)
# Cell centers at .25/.75 — matches population/PM data convention in Park et al.
r_target <- rast(
  nrows = 360, ncols = 720,
  xmin = -180, xmax = 180,
  ymin = -90,  ymax = 90,
  crs = "EPSG:4326"
)

# Regrid each year separately using area-weighted sum, then recombine
pop_regrid_list <- list()

for (i in 1:6) {
  yr <- years[i]
  col_name <- paste0("park_", yr, "_pop")
  
  # Extract single year from long format
  df_yr <- pop_long %>%
    filter(year == yr) %>%
    select(lon, lat, population)
  
  # Convert to SpatRaster (source grid)
  r_source <- rast(df_yr, type = "xyz", crs = "EPSG:4326")
  
  # Resample to target grid using sum (area-weighted, preserves total population)
  r_resampled <- resample(r_source, r_target, method = "sum")
  
  # Convert back to data frame
  df_out <- as.data.frame(r_resampled, xy = TRUE) %>%
    rename(lon = x, lat = y) %>%
    rename(!!col_name := population)
  
  pop_regrid_list[[i]] <- df_out
}

# Join all years by lon/lat into wide format
pop_wide <- pop_regrid_list[[1]]
for (i in 2:6) {
  pop_wide <- left_join(pop_wide, pop_regrid_list[[i]], by = c("lon", "lat"))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Verify and Save ################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Verify coordinate range matches target grid
head(pop_wide, 30)

# Verify total population is preserved after regridding
pop_wide %>%
  select(starts_with("park_")) %>%
  colSums()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ US Population Sum (any cell overlap with US boundary) ##########
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load US boundary
usa <- ne_countries(country = "United States of America", returnclass = "sf") |>
  st_transform(4326)

# First filter to the CONUS bounding box
# This greatly reduces the number of cells before spatial filtering
pop_sub <- pop_wide %>%
  filter(
    lon >= -125, lon <= -66,
    lat >=   24, lat <=  50
  )

# Vectorized grid-based intersection — much faster than rowwise() mutate
# Instead of building a polygon for every cell, create a 0.5-degree raster grid once
r <- rast(
  xmin = min(pop_sub$lon) - 0.25,
  xmax = max(pop_sub$lon) + 0.25,
  ymin = min(pop_sub$lat) - 0.25,
  ymax = max(pop_sub$lat) + 0.25,
  resolution = 0.5,
  crs = "EPSG:4326"
)

# Rasterize the US boundary onto the grid
# touches = TRUE approximates the rule that any overlap counts
usa_r <- rasterize(
  vect(usa),
  r,
  field = 1,
  background = NA,
  touches = TRUE
)

# Match each lon/lat point to its grid cell
# This is much faster than st_buffer() + st_intersects() on many polygons
cell_id <- cellFromXY(usa_r, as.matrix(pop_sub[, c("lon", "lat")]))

# Filter to cells that intersect US boundary (any overlap counts)
keep <- !is.na(values(usa_r)[cell_id, 1])

# Sum population for each year within US-touching cells
us_pop_sums <- pop_sub[keep, ] %>%
  select(starts_with("park_")) %>%
  colSums(na.rm = TRUE)

cat("US population sums:\n")
print(us_pop_sums)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############  Save ##########################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Save data
write.csv(pop_wide, here("output", "pop_regrid_park_2024.csv"), row.names = FALSE)

##THE END