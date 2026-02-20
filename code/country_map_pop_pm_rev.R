# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
########### Map data frame with latitude and longitude coordinates to countries ##########
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Remove all objects from the environment
rm(list = ls())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Packages #######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(here)
library(tidyverse)
library(sf)              # Spatial data handling
library(rnaturalearth)   # Country boundary data

# (NEW) for 0.1° netCDF aggregation + fast joins
library(ncdf4)
library(data.table)
library(dplyr)
library(leaflet)

# Turn off spherical geometry (s2) to allow geometry repairs
# s2 is strict about geometry validity, which prevents repairs from working
sf_use_s2(FALSE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Import #########################################################
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
# 4. (NEW) Decision rules using 0.1° gridded country/land product
#    - For each 0.5° cell center, evaluate its 25 underlying 0.1° subcells (5x5)
#    - Rule 1 (majority): if a strict majority of the 25 subcells belong to one country,
#      assign that country to the 0.5° cell (in ISO2 if the netCDF provides ISO2)
#    - Rule 2 (land-any): if any of the 25 subcells overlaps land, treat the whole 0.5° cell
#      as land (not ocean), used as a guard against coastal/ocean misclassification
#
# 5. Spatial join using st_intersects()
#    - Tests which country polygon contains each grid cell center point
#    - Returns country code/name for each matched cell
#    - Unmapped cells = ocean, disputed territories, or geometry issues
#    - (NEW patch) if st_join returns NA but rules say land + majority ISO2 exists,
#      convert ISO2 -> ISO3 (via Natural Earth) and fill iso_a3/name
#
# 6. Visual validation
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
# Keep BOTH iso_a3 and iso_a2 so we can convert ISO2 -> ISO3 when needed
world_subset <- world %>%
  select(iso_a3, iso_a2, name, geometry)  # iso_a3 = ISO3, iso_a2 = ISO2

# Fix invalid geometries (self-intersections, duplicate vertices, etc.)
# st_buffer with distance 0 is a robust way to repair geometry issues
world_subset <- st_buffer(world_subset, dist = 0)
# warning message: issue with buffering lon/lat in planar mode, but not important since using dist = 0

cat("Loaded country boundaries:", nrow(world_subset), "countries\n")

# Build ISO2 -> ISO3 lookup (Natural Earth)
iso2_to_iso3 <- world_subset %>%
  st_drop_geometry() %>%
  distinct(iso_a2, iso_a3, name) %>%
  filter(!is.na(iso_a2), iso_a2 != "", !is.na(iso_a3), iso_a3 != "")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# (NEW) Decision rules: 0.1° -> 0.5° (inserted between Step 3 and Step 5)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n=== NEW: Decision rules from 0.1° grid to 0.5° cells ===\n")

# Read 0.1° gridded product (country/land)
nc_path <- here("input/landmask_area", "countries_gridded_0.1deg_v0.1.nc")
nc <- nc_open(nc_path)

cat("=== Method 2: Extent check ===\n")
lon01 <- ncvar_get(nc, "lon")
lat01 <- ncvar_get(nc, "lat")
cat("lat range:", min(lat01), "to", max(lat01), " (n =", length(lat01), ")\n")
cat("lon range:", min(lon01), "to", max(lon01), " (n =", length(lon01), ")\n")
cat("unique dlat:", paste(unique(round(diff(lat01), 10)), collapse = ", "), "\n")
cat("unique dlon:", paste(unique(round(diff(lon01), 10)), collapse = ", "), "\n")

# country is expected 2D: dims are lon,lat (or lat,lon)
country2d <- ncvar_get(nc, "country")
dm <- dim(country2d)
cat("country var dims:", paste(dm, collapse = " x "), "\n")

# If netCDF provides an ISO lookup in dimension "iso", capture it (optional)
iso_vals <- NULL
if ("iso" %in% names(nc$dim)) {
  iso_vals <- nc$dim$iso$vals
  if (!is.null(iso_vals)) {
    iso_vals <- as.character(iso_vals)
  }
}

nc_close(nc)

# Build 0.1° long table (lon, lat, country_id, land)
lon_len <- length(lon01)
lat_len <- length(lat01)

if (all(dm[1:2] == c(lon_len, lat_len))) {
  # [lon, lat]
  dt01 <- as.data.table(as.data.frame.table(country2d))
  setnames(dt01, c("lon_i", "lat_i", "country_id"))
  dt01[, lon := lon01[as.integer(lon_i)]]
  dt01[, lat := lat01[as.integer(lat_i)]]
} else if (all(dm[1:2] == c(lat_len, lon_len))) {
  # [lat, lon]
  dt01 <- as.data.table(as.data.frame.table(country2d))
  setnames(dt01, c("lat_i", "lon_i", "country_id"))
  dt01[, lon := lon01[as.integer(lon_i)]]
  dt01[, lat := lat01[as.integer(lat_i)]]
} else {
  stop("Unexpected dimensions for 'country' var vs lon/lat lengths.")
}

dt01[, lon := round(lon, 1)]
dt01[, lat := round(lat, 1)]

# country_id may be numeric IDs, ISO2 codes, or something else; keep as character
dt01[, country_id := as.character(country_id)]

# If country_id looks numeric AND we have iso_vals, map numeric index -> iso code (often ISO2)
# This step is optional and only runs when it is safe.
if (!is.null(iso_vals)) {
  # detect numeric-like ids (e.g., "0","1","2",...)
  is_num_like <- suppressWarnings(!is.na(as.integer(dt01$country_id)))
  if (any(is_num_like)) {
    # Map indices only for those >0 and within iso_vals length
    idx <- suppressWarnings(as.integer(dt01$country_id))
    ok <- !is.na(idx) & idx >= 1 & idx <= length(iso_vals)
    # Preserve "0" (ocean) as "0"
    dt01[ok, country_id := iso_vals[idx[ok]]]
  }
}

# Land rule: treat any non-missing, non-ocean code as land
# Common ocean codes: "0", "00", "" (keep all three guards)
dt01[, land := as.integer(!is.na(country_id) & country_id != "" & country_id != "NA" & country_id != "0" & country_id != "00")]

dt01 <- dt01[, .(lon, lat, country_id, land)]
setkey(dt01, lon, lat)

cat("0.1° cells in dt01:", nrow(dt01), "\n")
cat("Sample country_id values:", paste(head(unique(dt01$country_id), 20), collapse = ", "), "\n")

# Expand each 0.5° center to its 25 subcells (5x5 on the 0.1° grid, by index)
dt05 <- as.data.table(pop_pm_combined)
stopifnot(all(c("lon", "lat") %in% names(dt05)))

lon_vec <- lon01
lat_vec <- lat01
nlon <- length(lon_vec)
nlat <- length(lat_vec)

# nearest index on a sorted grid vector
nearest_idx <- function(x, v) {
  i <- findInterval(x, v)
  i[i < 1] <- 1
  i[i >= length(v)] <- length(v) - 1
  left <- v[i]
  right <- v[i + 1]
  i + (abs(right - x) < abs(x - left))
}

# force a 5-index window centered near i0, clipped to valid range
idx_window5 <- function(i0, n) {
  i0 <- pmax(1L, pmin(as.integer(i0), n))
  istart <- pmax(1L, i0 - 2L)
  iend <- pmin(n, istart + 4L)
  istart <- pmax(1L, iend - 4L)
  list(istart = istart, iend = iend)
}

# base table: 0.5 centers + nearest 0.1 indices
dt_base <- dt05[, .(lon05 = lon, lat05 = lat)]
dt_base[, lon_i0 := nearest_idx(lon05, lon_vec)]
dt_base[, lat_i0 := nearest_idx(lat05, lat_vec)]

wlon <- idx_window5(dt_base$lon_i0, nlon)
wlat <- idx_window5(dt_base$lat_i0, nlat)

dt_base[, lon_i_start := wlon$istart]
dt_base[, lat_i_start := wlat$istart]

# offsets in index space (exactly 5x5)
offs_i <- CJ(dlon_i = 0:4, dlat_i = 0:4)

# cartesian expand: each 0.5 cell -> 25 subcells (all exactly on lon01/lat01)
dt_exp <- data.table(
  lon05 = rep(dt_base$lon05, each = nrow(offs_i)),
  lat05 = rep(dt_base$lat05, each = nrow(offs_i)),
  lon_i = rep(dt_base$lon_i_start, each = nrow(offs_i)) + rep(offs_i$dlon_i, times = nrow(dt_base)),
  lat_i = rep(dt_base$lat_i_start, each = nrow(offs_i)) + rep(offs_i$dlat_i, times = nrow(dt_base))
)

# map indices -> exact 0.1 grid coordinates
dt_exp[, lon01 := lon_vec[lon_i]]
dt_exp[, lat01 := lat_vec[lat_i]]

# join subcells to 0.1 grid lookup dt01 (now exact matches, no rounding needed)
setkey(dt01, lon, lat)
dt_exp <- dt01[dt_exp, on = .(lon = lon01, lat = lat01)]

# Decision rules:
# - land_any: if any 0.1° subcell is land -> treat the 0.5° cell as land
# - country_id_major:
#     (1) if a strict majority (>=13 of 25) exists among LAND subcells, assign it
#     (2) otherwise, if land_any == TRUE, assign the plurality winner among LAND subcells
#         (with deterministic tie-break)

get_majority <- function(country_id, land_flag) {
  x <- country_id[
    land_flag == 1 &
      !is.na(country_id) & country_id != "" & country_id != "NA" &
      country_id != "0" & country_id != "00"
  ]
  if (length(x) == 0) return(NA_character_)
  tab <- sort(table(x), decreasing = TRUE)
  if (tab[1] >= 13) return(names(tab)[1])
  NA_character_
}

# UPDATED: deterministic tie-break for plurality (alphabetical among ties)
get_plurality <- function(country_id, land_flag) {
  x <- country_id[
    land_flag == 1 &
      !is.na(country_id) & country_id != "" & country_id != "NA" &
      country_id != "0" & country_id != "00"
  ]
  if (length(x) == 0) return(NA_character_)
  tab <- sort(table(x), decreasing = TRUE)
  
  max_n <- tab[1]
  winners <- names(tab)[tab == max_n]
  winners <- sort(winners)   # tie-break: alphabetical
  winners[1]
}

dt_rule <- dt_exp[, .(
  land_any = any(land == 1, na.rm = TRUE),
  country_id_major = {
    maj <- get_majority(country_id, land)
    if (!is.na(maj)) maj else if (any(land == 1, na.rm = TRUE)) get_plurality(country_id, land) else NA_character_
  }
), by = .(lon05, lat05)]

# Add rule columns back to original data (structure unchanged + 2 new columns)
pop_pm_combined <- pop_pm_combined %>%
  left_join(as.data.frame(dt_rule), by = c("lon" = "lon05", "lat" = "lat05"))

# Rebuild sf points so the new columns flow into the sf workflow (no other structure change)
pop_pm_sf <- pop_pm_combined %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

cat("Added columns: land_any, country_id_major\n")
cat("land_any TRUE count:", sum(pop_pm_combined$land_any, na.rm = TRUE), "\n")
cat("majority country assigned count:", sum(!is.na(pop_pm_combined$country_id_major)), "\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Perform spatial join
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Match each grid cell point to country using point-in-polygon intersection
# st_intersects() checks which country polygon each point falls within

pop_pm_with_countries <- st_join(pop_pm_sf, world_subset, join = st_intersects)

pop_pm_with_countries <- pop_pm_with_countries %>%
  mutate(
    iso_a3 = if_else(is.na(iso_a3) & land_any, "-99", iso_a3),
    name   = if_else(is.na(name)   & land_any, "Unmapped land (rule land_any=TRUE)", name)
  )

# (NEW) Patch: if st_join is NA but rules say land + majority ISO2 exists,
# convert ISO2 -> ISO3 and fill iso_a3/name
# This keeps your final output in ISO3.
pop_pm_with_countries <- pop_pm_with_countries %>%
  left_join(iso2_to_iso3, by = c("country_id_major" = "iso_a2"), suffix = c("", "_maj")) %>%
  mutate(
    iso_a3 = if_else(is.na(iso_a3) & land_any & !is.na(country_id_major) & !is.na(iso_a3_maj),
                     iso_a3_maj, iso_a3),
    name   = if_else(is.na(name)   & land_any & !is.na(country_id_major) & !is.na(name_maj),
                     name_maj, name)
  ) %>%
  select(-iso_a3_maj, -name_maj)

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
  here("images", "country_mapping_check_rev.png"),
  width = 12,
  height = 6,
  dpi = 300
)

cat("Map saved to: images/country_mapping_check_rev.png\n")

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

pop_pm_final <- pop_pm_final %>%
  select(-any_of(c("iso_a2", "country_id_major")))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Save final output
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

write_csv(pop_pm_final, here("output", "pop_pm_with_countries_rev.csv"))

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
############ interative map  ################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# function: center -> 0.5° grid polygon
make_cell_polygon <- function(lon, lat, half = 0.25) {
  coords <- matrix(c(
    lon - half, lat - half,
    lon + half, lat - half,
    lon + half, lat + half,
    lon - half, lat + half,
    lon - half, lat - half
  ), ncol = 2, byrow = TRUE)
  st_polygon(list(coords))
}

df <- pop_pm_final

# helpers
is_invalid <- function(x) {
  !is.na(x) & x == "-99"
}
is_unmapped <- function(x) {
  is.na(x) | x == ""
}

# ---- mapped: ONLY North America (ISO3) ----
na_iso3 <- c("USA", "CAN", "MEX")
mapped_df <- df %>%
  filter(!is_unmapped(country_code_iso3),
         !is_invalid(country_code_iso3),
         country_code_iso3 %in% na_iso3)

# ---- invalid: -99 (brown) ----
invalid_df <- df %>%
  filter(country_code_iso3 == "-99",
         !is.na(land_any), land_any == TRUE)

# ---- unmapped: NA or "" (yellow) ----
unmapped_df <- df %>%
  filter(is_unmapped(country_code_iso3),
         is.na(land_any) | land_any == FALSE)

# ---- quick diagnostics ----
cat("Total cells:", nrow(df), "\n")
cat("Mapped (all non-NA):", sum(!is.na(df$country_code_iso3)), "\n")
cat("Unmapped (NA or blank):", nrow(unmapped_df), "\n")
cat("Invalid (-99):", nrow(invalid_df), "\n")
cat("Mapped (NA only):", nrow(mapped_df), "\n")

# ---- sample ONLY mapped for performance (safe) ----
set.seed(1)
if (nrow(mapped_df) > 40000) {
  mapped_df <- mapped_df[sample(nrow(mapped_df), 100000), ]
}

# build mapped polygons (NA only)
mapped_polys <- mapply(make_cell_polygon, mapped_df$lon, mapped_df$lat, SIMPLIFY = FALSE)
mapped_cells <- st_sf(mapped_df, geometry = st_sfc(mapped_polys, crs = 4326))

# build invalid polygons (-99)
invalid_polys <- mapply(make_cell_polygon, invalid_df$lon, invalid_df$lat, SIMPLIFY = FALSE)
invalid_cells <- st_sf(invalid_df, geometry = st_sfc(invalid_polys, crs = 4326))

# build unmapped polygons (GLOBAL, full)
unmapped_polys <- mapply(make_cell_polygon, unmapped_df$lon, unmapped_df$lat, SIMPLIFY = FALSE)
unmapped_cells <- st_sf(unmapped_df, geometry = st_sfc(unmapped_polys, crs = 4326))

# land boundary (Natural Earth)
land <- ne_countries(scale = 10, returnclass = "sf") %>%
  st_transform(4326)

# palette for NA ISO3
pal_na <- colorFactor(
  palette = c("#1b7837", "#2166ac", "#b2182b"),
  domain = na_iso3
)

leaflet() %>%
  addProviderTiles("CartoDB.Positron") %>%
  
  # show the whole world
  fitBounds(lng1 = -180, lat1 = -60, lng2 = 180, lat2 = 85) %>%
  
  # land boundary
  addPolylines(
    data = land,
    color = "red",
    weight = 1,
    opacity = 0.7,
    fill = FALSE,
    group = "Land boundary"
  ) %>%
  
  # mapped NA cells
  addPolygons(
    data = mapped_cells,
    color = ~pal_na(country_code_iso3),
    weight = 0.4,
    fillColor = ~pal_na(country_code_iso3),
    fillOpacity = 0.40,
    stroke = TRUE,
    group = "Mapped cells (NA only)",
    popup = ~paste0(
      "<b>Mapped cell (NA only)</b><br>",
      "ISO3: ", country_code_iso3, "<br>",
      "Country: ", country_name, "<br>",
      "Lon: ", round(lon, 2), "<br>",
      "Lat: ", round(lat, 2)
    )
  ) %>%
  
  # invalid -99 cells (brown)
  addPolygons(
    data = invalid_cells,
    color = "#5c3b1e",
    weight = 0.8,
    fillColor = "#8b5a2b",
    fillOpacity = 0.70,
    stroke = TRUE,
    group = "Invalid cells (-99)",
    popup = ~paste0(
      "<b>Invalid cell (-99)</b><br>",
      "Lon: ", round(lon, 2), "<br>",
      "Lat: ", round(lat, 2)
    )
  ) %>%
  
  # unmapped NA/blank cells (yellow)
  addPolygons(
    data = unmapped_cells,
    color = "black",
    weight = 1,
    fillColor = "yellow",
    fillOpacity = 0.7,
    group = "Unmapped cells (NA/blank)",
    popup = ~paste0(
      "<b>Unmapped cell</b><br>",
      "Lon: ", round(lon, 2), "<br>",
      "Lat: ", round(lat, 2)
    )
  ) %>%
  
  addLegend(
    position = "bottomright",
    pal = pal_na,
    values = mapped_cells$country_code_iso3,
    title = "Mapped ISO3 (NA only)",
    opacity = 1
  ) %>%
  
  addLayersControl(
    overlayGroups = c("Land boundary", "Mapped cells (NA only)", "Invalid cells (-99)", "Unmapped cells (NA/blank)"),
    options = layersControlOptions(collapsed = FALSE)
  )
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ THE END  #######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# THE END