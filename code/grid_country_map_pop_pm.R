# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
########### Map data frame with latitude and longitude coordinates to countries ##################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Remove all objects from the environment
rm(list = ls())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Packages #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(here)
library(tidyverse)
library(sf)              # Spatial data handling
library(rnaturalearth)   # Country boundary data

# Turn off spherical geometry (s2) to allow geometry repairs
# s2 is strict about geometry validity, which prevents repairs from working
sf_use_s2(FALSE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Import #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Import combined PM2.5 data and pop data (PM regridded, combined with pop, stored in data frame)
pop_pm_combined <- read_csv(here("output", "pop_pm_combined.csv"))

# pop data dimension = (lon x lat) 720 × 360 = 259200

nrow(pop_pm_combined)
colnames(pop_pm_combined)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ MAP GRID CELLS TO COUNTRIES ####################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# MAPPING APPROACH OVERVIEW:
#
# This script assigns each 0.5° grid cell (with lon/lat coordinates) to a 
# country using spatial point-in-polygon matching. The approach:
#
# 1. Convert dataframe coordinates --> spatial points (sf objects)
#    - Transforms lon/lat pairs into geometric point objects
#    - Assigns WGS84 coordinate reference system (EPSG:4326)
#
# 2. Load country boundary polygons from Natural Earth (scale = 10)
#    - High-resolution boundaries for accurate coastal/island coverage
#    - Includes maritime zones and territorial waters
#
# 3. Repair invalid geometries with st_buffer(dist = 0)
#    - Natural Earth data has topology errors (self-intersections, duplicate vertices)
#    - Zero-distance buffer cleans geometries without changing boundaries
#    - Requires s2 turned OFF (planar geometry) to allow repairs
#
# 4. Spatial join using st_intersects()
#    - Tests which country polygon contains each grid cell center point
#    - Returns country code/name for each matched cell
#    - Unmapped cells = ocean, disputed territories, or geometry issues
#
# 5. Visual validation
#    - Map shows countries with successfully matched cells (blue)
#    - Identifies missing countries or large unmapped regions
#
# IMPORTANT NOTES:
# - Grid cell centers near coastlines may fall in maritime zones
# - This can assign ocean cells to countries (intended behavior for territorial waters)
# - Planar geometry (s2 OFF) sacrifices some accuracy but enables mapping
# - For 0.5° cells and country-level analysis, accuracy loss is negligible
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat("\n=== Starting country mapping ===\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Convert coordinates to spatial points
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Convert lat/lon to sf points object with WGS84 coordinate system
# crs = 4326 specifies EPSG:4326 (WGS84), the standard geographic coordinate system
# remove = FALSE keeps original lon/lat columns in the dataframe

pop_pm_sf <- pop_pm_combined %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

cat("Created spatial points object:", nrow(pop_pm_sf), "points\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load country boundaries
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Use high-resolution Natural Earth boundaries (scale = 10)
# Scale 10 = 1:10m resolution, provides better coverage of coastal areas and small islands
# Scale 50 or 110 would be lower resolution with less accurate coastlines

world <- ne_countries(scale = 10, returnclass = "sf")

# Select only needed columns to keep data manageable
world_subset <- world %>%
  select(iso_a3, name, geometry)  # iso_a3 = ISO 3166-1 alpha-3 country codes

# Fix invalid geometries (self-intersections, duplicate vertices, etc.)
# st_buffer with distance 0 is a robust way to repair geometry issues
world_subset <- st_buffer(world_subset, dist = 0)
# warning message: issue with buffering lon/lat in planar mode, but not important since using dist = 0 

cat("Loaded country boundaries:", nrow(world_subset), "countries\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Perform spatial join
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Match each grid cell point to country using point-in-polygon intersection
# st_intersects() checks which country polygon each point falls within

pop_pm_with_countries <- st_join(pop_pm_sf, world_subset, join = st_intersects)

# Check how many cells were successfully mapped
mapped_count <- sum(!is.na(pop_pm_with_countries$iso_a3))
unmapped_count <- sum(is.na(pop_pm_with_countries$iso_a3))

cat("\nSpatial join results:\n")
cat("  Mapped cells:", mapped_count, "\n")
cat("  Unmapped cells (ocean/disputes/other):", unmapped_count, "\n")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Visual check: Plot mapped countries
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create map showing which countries were successfully matched to grid cells
# This helps identify any major gaps or mapping errors

# Get unique country codes from matched cells
matched_countries <- pop_pm_with_countries %>%
  filter(!is.na(iso_a3)) %>%
  pull(iso_a3) %>%
  unique()

cat("\nNumber of unique countries matched:", length(matched_countries), "\n")

# Filter world boundaries to only countries that were matched
world_matched <- world_subset %>%
  filter(iso_a3 %in% matched_countries)

# Create global map
ggplot() +
  # All countries in light gray as background
  geom_sf(data = world_subset, fill = "gray90", color = "white", linewidth = 0.1) +
  # Highlight matched countries in blue
  geom_sf(data = world_matched, fill = "steelblue", color = "white", linewidth = 0.1) +
  theme_minimal() +
  labs(
    title = "Countries with Mapped Grid Cells",
    subtitle = paste0(length(matched_countries), " countries matched"),
    caption = "Blue = countries with mapped cells, Gray = no cells"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank()
  )

# Save the plot
ggsave(
  here("images", "country_mapping_check.png"),
  width = 12, 
  height = 6, 
  dpi = 300
)

cat("Map saved to: images/country_mapping_check.png\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Convert back to regular dataframe and clean
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Remove geometry column (sf-specific, not needed for analysis)
# Rename country columns for consistency with previous workflow

pop_pm_final <- pop_pm_with_countries %>%
  st_drop_geometry() %>%
  rename(
    country_code_iso3 = iso_a3,
    country_name = name
  )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Save final output
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

write_csv(pop_pm_final, here("output", "pop_pm_with_countries.csv"))

cat("\nCountry mapping complete!\n")
cat("Output saved to: pop_pm_with_countries.csv\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Summary statistics
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat("\n=== Summary by Country (top 10 by cell count) ===\n")
pop_pm_final %>%
  filter(!is.na(country_code_iso3)) %>%
  count(country_name, country_code_iso3) %>%
  arrange(desc(n)) %>%
  head(10) %>%
  print()

cat("\n=== Unmapped cells by region ===\n")
pop_pm_final %>%
  filter(is.na(country_code_iso3)) %>%
  group_by(
    lat_band = cut(lat, breaks = seq(-90, 90, by = 30), include.lowest = TRUE)
  ) %>%
  summarize(
    n_cells = n(),
    .groups = "drop"
  ) %>%
  print()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ JOIN TO GIVE MODEL COUNTRIES  ####################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ THE END  ####################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# THE END 
