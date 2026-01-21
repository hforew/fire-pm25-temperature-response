#########################################
## Create grid cell to country mapping ##
#########################################

# Remove all objects from the environment
rm(list = ls())

############ Packages #####################################################

library(here)
library(ncdf4)
library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

############ Import data ##################################################

# Import the NetCDF file
pm25_data <- nc_open(here("input", "CESM_09x125_PM25_2000_Baseline.nc"))

# View the file structure (optional)
print(pm25_data)

# Load the list of countries to map (from GIVE)
countries_list <- read_csv(here("input", "GIVE_countries.csv"))


############ Map grid cell to country ##################################################

# Extract longitude and latitude
lon <- ncvar_get(pm25_data, "lon")
length(lon)

lat <- ncvar_get(pm25_data, "lat")
length(lat)

nc_close(pm25_data)

# Create grid dataframe with all lon/lat combinations
grid_df <- expand.grid(
  lon = lon,
  lat = lat
) %>%
  mutate(grid_id = row_number())

# Convert grid to sf points object
grid_sf <- st_as_sf(grid_df, 
                    coords = c("lon", "lat"), 
                    crs = 4326,
                    remove = FALSE)

# Get world country boundaries
world <- ne_countries(scale = "medium", returnclass = "sf")


# Filter to only countries in your list -- multiple matching strategies
countries_sf <- world %>%
  filter(
    iso_a3 %in% countries_list$ISO3 |           # Match by ISO3
      admin %in% countries_list$country |          # Match by admin name
      name %in% countries_list$country |           # Match by name
      sovereignt %in% countries_list$country       # Match by sovereignty
  )

# Spatial join: map each grid cell to a country
grid_country_mapping <- st_join(grid_sf, 
                                countries_sf %>% 
                                  select(name, iso_a3, iso_a2),
                                join = st_within)

# Convert back to regular dataframe and fix broken ISO codes
grid_lookup <- grid_country_mapping %>%
  st_drop_geometry() %>%
  mutate(
    lon = as.numeric(lon),
    lat = as.numeric(lat)
  ) %>%
  rename(country_name = name,
         country_code_iso3 = iso_a3,
         country_code_iso2 = iso_a2) %>%
  # Fix broken ISO codes by matching country names back to CSV
  left_join(
    countries_list %>% select(country, ISO3),
    by = c("country_name" = "country")
  ) %>%
  mutate(
    country_code_iso3 = ifelse(country_code_iso3 == "-99" & !is.na(ISO3), 
                               ISO3, 
                               country_code_iso3)
  ) %>%
  select(-ISO3)  # Remove the temporary join column

# Check if FRA and NOR codes are fixed
grid_lookup %>%
  filter(country_name %in% c("France", "Norway")) %>%
  distinct(country_name, country_code_iso3) %>%
  arrange(country_name)

# Check for any remaining "-99" codes
grid_lookup %>%
  filter(country_code_iso3 == "-99") %>%
  distinct(country_name, country_code_iso3)

# Verify the structure of grid_lookup
colnames(grid_lookup)
head(grid_lookup)

# Save the lookup table for reuse
saveRDS(grid_lookup, here("output", "grid_country_lookup.rds"))
write_csv(grid_lookup, here("output", "grid_country_lookup.csv"))


################## Analyze results of grid cell to country mapping ##################

# Summary statistics
summary_stats <- grid_lookup %>%
  group_by(country_name) %>%
  summarise(n_grid_cells = n()) %>%
  arrange(desc(n_grid_cells))

print(summary_stats)

# Check for unmapped grid cells
unmapped <- grid_lookup %>%
  filter(is.na(country_name)) %>%
  nrow()

total_cells <- nrow(grid_lookup)
fraction_unmapped <- unmapped / total_cells

print(paste("Number of unmapped grid cells:", unmapped))
print(paste("Total grid cells:", total_cells))
print(paste("Fraction unmapped:", round(fraction_unmapped, 3)))

# Count unique country codes
n_countries <- grid_lookup %>%
  filter(!is.na(country_code_iso3)) %>%
  distinct(country_code_iso3) %>%
  nrow()

print(paste("Number of unique countries mapped:", n_countries))

# Check which countries from CSV are not mapped in grid_lookup
countries_in_csv <- countries_list$ISO3
# Now recheck missing countries
countries_in_grid <- grid_lookup %>%
  filter(!is.na(country_code_iso3), country_code_iso3 != "-99") %>%
  distinct(country_code_iso3) %>%
  pull(country_code_iso3)

missing_countries <- countries_list$ISO3[!countries_list$ISO3 %in% countries_in_grid]

print(paste("Number of countries in CSV not mapped:", length(missing_countries)))
print(missing_countries)

# Show with country names for easier interpretation
missing_with_names <- countries_list %>%
  filter(ISO3 %in% missing_countries)

print(missing_with_names)
print(missing_with_names, n = Inf)

