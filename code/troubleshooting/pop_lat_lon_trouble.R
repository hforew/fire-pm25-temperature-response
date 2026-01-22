

# Remove all objects from the environment
rm(list = ls())


# packages

library(ncdf4)
library(here)
library(tidyverse)

# troubleshooting:
    # - 1) cells not mapped to countries with high populations
    # - 2) higher than expected pop total for 2009


# Open the NetCDF file
pop_ssp1 <- nc_open(here("input", "SSP1_for_RCP45_2006-2100_population_density_c160701.nc"))

# Get coordinates
lon <- ncvar_get(pop_ssp1, "lon")
lat <- ncvar_get(pop_ssp1, "lat")

# Find the index for lon = 73.25, lat = 8.25
lon_idx <- which(lon == 73.25)
lat_idx <- which(lat == 8.25)

cat("Lon index:", lon_idx, "Lat index:", lat_idx, "\n")

# Get the 2009 year index
years <- ncvar_get(pop_ssp1, "year")
year_idx <- which(years == 2009)

# Extract the specific cell value
pop_density_cell <- ncvar_get(
  pop_ssp1,
  "hdm",
  start = c(lon_idx, lat_idx, year_idx),
  count = c(1, 1, 1)
)

cat("Population density at (73.25, 8.25) in 2009:", pop_density_cell, "people/km²\n")

# the population density for this location (middle of indian ocean) is too high

# Check neighboring cells to see if this is an isolated anomaly
library(tidyverse)

pop_df <- read_csv(here("output", "pop_df.csv"))

# Look at a 5x5 grid centered on the problematic cell
pop_df %>%
  filter(lon >= 71.25, lon <= 75.25,
         lat >= 6.25, lat <= 10.25) %>%
  select(lon, lat, pop_dens_2009, pop_tot_2009) %>%
  arrange(lat, lon) %>%
  print(n = 100)


# Find all cells with >10 million people
high_pop_cells <- pop_df %>%
  filter(pop_tot_2009 > 10e6) %>%
  arrange(desc(pop_tot_2009)) %>%
  select(lon, lat, pop_dens_2009, pop_tot_2009)

cat("\n=== CELLS WITH >10 MILLION PEOPLE ===\n")
cat("Total cells with >10M people:", nrow(high_pop_cells), "\n\n")

cat("Top 20 most populous cells:\n")
print(high_pop_cells %>% head(20), n = 20)

cat("\n=== ANALYSIS ===\n")
cat("Largest single cell:", format(max(high_pop_cells$pop_tot_2009), big.mark = ","), "people\n")
# is might be Male, a very densely populated island in Indian Ocean. 
# BUT this density would NOT apply over the entire grid cell area


