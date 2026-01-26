# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
########### Map data frame with latitude and longitude coordinates to countries ##################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Remove all objects from the environment
rm(list = ls())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Packages #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("C:/Ford_BA_FPM25/")
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
# ggsave(
#  here("images", "country_mapping_check.png"),
#  width = 12, 
#  height = 6, 
#  dpi = 300
#)

#cat("Map saved to: images/country_mapping_check.png\n")

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
# Country code check
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# View all entries with country code -99
pop_pm_final %>%
  filter(country_code_iso3 == "-99") %>%
  count(country_name, country_code_iso3) %>%
  arrange(country_name)

# Manually correct country codes for France and Norway
pop_pm_final <- pop_pm_final %>%
  mutate(
    country_code_iso3 = case_when(
      country_name == "France" & country_code_iso3 == "-99" ~ "FRA",
      country_name == "Norway" & country_code_iso3 == "-99" ~ "NOR",
      TRUE ~ country_code_iso3  # Keep all other values unchanged
    )
  )

# Verify the corrections
pop_pm_final %>%
  filter(country_name %in% c("France", "Norway")) %>%
  count(country_name, country_code_iso3)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Save final output
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# write_csv(pop_pm_final, here("output", "pop_pm_with_countries.csv"))

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
############ THE END  ####################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# THE END 

################################################################################review start
#Concern:
#Because the code assigns countries based only on whether the grid cell center point falls inside a country polygon, 
#coastal and border cells whose centers lie outside the polygon but whose areas (grid cell) contain land maybe misclassified as NA, 
#causing systematic undercounting of land and population.

# ------------------------------------------------------------------------------
#statistics may need
library(dplyr)

df <- pop_pm_final

# Define what counts as "unmapped"
# if ti is correct that a grid cell is considered unmapped if the country code is NA, "-99", or empty
is_unmapped <- function(x) is.na(x) | x == "-99" | x == ""

# 1) Fraction of ALL grid cells that are unmapped
n_all <- nrow(df)
n_unmapped_all <- sum(is_unmapped(df$country_code_iso3))

frac_unmapped_all <- n_unmapped_all / n_all

# 2) Fraction of grid cells with population > 0 that are unmapped
# These are the only cells that actually matter for economic story
df_pop <- df %>% filter(pop_tot_2009 > 0)

n_pop <- nrow(df_pop)
n_unmapped_pop <- sum(is_unmapped(df_pop$country_code_iso3))

frac_unmapped_pop <- n_unmapped_pop / n_pop

# 3) Population-weighted loss
# This measures how much real population is assigned to unmapped cells
total_pop <- sum(df$pop_tot_2009, na.rm = TRUE)
unmapped_pop <- sum(df$pop_tot_2009[is_unmapped(df$country_code_iso3)], na.rm = TRUE)

frac_pop_lost <- unmapped_pop / total_pop



# Print summary
cat("Total grid cells:", n_all, "\n")
cat("Unmapped cells:", n_unmapped_all, "\n")
cat("Fraction unmapped (all cells):", signif(frac_unmapped_all, 4), "\n\n")

cat("Cells with population > 0:", n_pop, "\n")
cat("Unmapped among populated:", n_unmapped_pop, "\n")
cat("Fraction unmapped (populated):", signif(frac_unmapped_pop, 4), "\n\n")

cat("Total population:", format(total_pop, scientific = FALSE), "\n")
cat("Population in unmapped cells:", format(unmapped_pop, scientific = FALSE), "\n")
cat("Fraction of population lost:", signif(frac_pop_lost, 5), "\n")

# ------------------------------------------------------------------------------
#land mark
# 1) Load a land mask (Earth land polygons)
land <- ne_download(scale = 10, type = "land", category = "physical", returnclass = "sf") %>%
  st_make_valid()

# 2) Take unmapped cells and convert to sf points (grid center points)
unmapped_pts <- df %>%
  filter(is_unmapped(country_code_iso3)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

# 3) Test whether each unmapped centroid falls on land
unmapped_pts$on_land <- lengths(st_intersects(unmapped_pts, land)) > 0

# 4) Summarize statistics
valid_summary <- unmapped_pts %>%
  st_drop_geometry() %>%
  summarize(
    n_total_cells          = nrow(df),
    n_unmapped_cells       = n(),
    frac_unmapped_all      = n() / nrow(df),
    
    n_unmapped_on_land     = sum(on_land),
    n_unmapped_ocean       = sum(!on_land),
    frac_unmapped_on_land  = ifelse(n() == 0, NA_real_, sum(on_land) / n()),
    
    pop_total              = sum(df$pop_tot_2009, na.rm = TRUE),
    pop_unmapped_total     = sum(pop_tot_2009, na.rm = TRUE),
    pop_unmapped_on_land   = sum(pop_tot_2009[on_land], na.rm = TRUE),
    pop_unmapped_ocean     = sum(pop_tot_2009[!on_land], na.rm = TRUE),
    
    frac_pop_lost_total    = ifelse(sum(df$pop_tot_2009, na.rm = TRUE) == 0, NA_real_,
                                    sum(pop_tot_2009, na.rm = TRUE) / sum(df$pop_tot_2009, na.rm = TRUE)),
    frac_pop_lost_on_land  = ifelse(sum(df$pop_tot_2009, na.rm = TRUE) == 0, NA_real_,
                                    sum(pop_tot_2009[on_land], na.rm = TRUE) / sum(df$pop_tot_2009, na.rm = TRUE))
  )

print(valid_summary)

# ------------------------------------------------------------------------------
library(leaflet)
#plot of unmapped cell
unmapped_pts <- df %>%
  filter(is_unmapped(country_code_iso3)) %>%
  st_as_sf(coords = c("lon","lat"), crs = 4326, remove = FALSE)

#convert to WGS84 for leaflet
unmapped_pts <- st_transform(unmapped_pts, 4326)
land <- st_transform(land, 4326)

#interactive map
leaflet() %>%
  addProviderTiles("CartoDB.Positron") %>%   # white/light basemap
  addPolylines(
    data = land,
    color = "red",
    weight = 1,
    opacity = 0.8,
    fill = FALSE,
    group = "Land boundary"
  ) %>%
  addCircleMarkers(
    data = unmapped_pts,
    radius = 2,
    color = "black",
    fillOpacity = 0.7,
    stroke = FALSE,
    group = "Unmapped cells",
    popup = ~paste0(
      "Lon: ", round(lon, 2), "<br>",
      "Lat: ", round(lat, 2), "<br>",
      "Pop 2009: ", format(pop_tot_2009, scientific = FALSE)
    )
  ) %>%
  addLayersControl(
    overlayGroups = c("Land boundary", "Unmapped cells"),
    options = layersControlOptions(collapsed = FALSE)
  )
################################################################################review end

