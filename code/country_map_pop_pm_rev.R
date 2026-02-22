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

# Import combined PM2.5 data and pop data (PM regridded, combined with pop)
pop_pm_combined <- read_csv(here("output", "pop_pm_combined.csv"))

# pop data dimension = (lon x lat) 720 × 360 = 259200
nrow(pop_pm_combined)
colnames(pop_pm_combined)


cat("\n=== Starting country mapping ===\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############Convert coordinates to spatial points############
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Convert lat/lon to sf points object with WGS84 coordinate system
# crs = 4326 specifies EPSG:4326 (WGS84), the standard geographic coordinate system
# remove = FALSE keeps original lon/lat columns in the dataframe

pop_pm_sf <- pop_pm_combined %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

cat("Created spatial points object:", nrow(pop_pm_sf), "points\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############Load country boundaries############
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

# ISO lookup that supports BOTH ISO2 and ISO3 matching
iso_lookup <- world_subset %>%
  st_drop_geometry() %>%
  distinct(iso_a2, iso_a3, name) %>%
  filter(!is.na(iso_a3), iso_a3 != "")

major_patch <- tibble::tribble(
  ~country_id_major, ~iso_a3_patch, ~name_patch,
  "FR", "FRA", "France",
  "NO", "NOR", "Norway",
  "TW", "TWN", "Taiwan",
  "GF", "GUF", "French Guiana",
  "GP", "GLP", "Guadeloupe",
  "MQ", "MTQ", "Martinique",
  "RE", "REU", "Réunion",
  "YT", "MYT", "Mayotte",
  "CX", "CXR", "Christmas Island",
  "BV", "BVT", "Bouvet Island",
  "SJ", "SJM", "Svalbard and Jan Mayen",
  # AN is deprecated; best modern mapping depends on netCDF meaning.
  # Often it's the old Netherlands Antilles; map to Curaçao (CUW) is a common practical choice,
  "AN", "CUW", "Netherlands Antilles (legacy → Curaçao)",
  # ZZ is "unknown/other"; keep it unmapped land
  "ZZ", NA_character_, "Unknown/Other (ZZ)"
)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############(NEW) Decision rules: 0.1° -> 0.5° (inserted between Step 3 and Step 5)############
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

# If netCDF provides an ISO lookup in dimension "iso", capture it
iso_vals <- NULL
if ("iso" %in% names(nc$dim)) {
  iso_vals <- nc$dim$iso$vals
  if (!is.null(iso_vals)) {
    iso_vals <- as.character(iso_vals)
  }
}

nc_close(nc)

lon_edge <- lon01
lat_edge <- lat01

lon_center <- (lon_edge[-1] + lon_edge[-length(lon_edge)]) / 2
lat_center <- (lat_edge[-1] + lat_edge[-length(lat_edge)]) / 2

# country2d should be cell-based (center grid): [nlon_center, nlat_center] or swapped
nlon_c <- length(lon_center)
nlat_c <- length(lat_center)

# If country2d is edge/node based (same dims as edges), trim to cell dims
if (all(dm[1:2] == c(length(lon_edge), length(lat_edge)))) {
  country2d <- country2d[1:nlon_c, 1:nlat_c]
  dm <- dim(country2d)
} else if (all(dm[1:2] == c(length(lat_edge), length(lon_edge)))) {
  country2d <- country2d[1:nlat_c, 1:nlon_c]
  dm <- dim(country2d)
}

# Build 0.1° lookup table on CENTER grid
if (all(dm[1:2] == c(nlon_c, nlat_c))) {
  dt01 <- as.data.table(as.data.frame.table(country2d))
  setnames(dt01, c("lon_i", "lat_i", "country_id"))
  dt01[, lon := lon_center[as.integer(lon_i)]]
  dt01[, lat := lat_center[as.integer(lat_i)]]
} else if (all(dm[1:2] == c(nlat_c, nlon_c))) {
  dt01 <- as.data.table(as.data.frame.table(country2d))
  setnames(dt01, c("lat_i", "lon_i", "country_id"))
  dt01[, lon := lon_center[as.integer(lon_i)]]
  dt01[, lat := lat_center[as.integer(lat_i)]]
} else {
  stop("Unexpected dimensions for 'country' var after trimming to center grid.")
}

# Keep exact centers rounded to 1 decimal (0.1°)
dt01[, lon := round(lon, 1)]
dt01[, lat := round(lat, 1)]
dt01[, country_id := as.character(country_id)]

# Map numeric ids -> ISO code if iso_vals exists
if (!is.null(iso_vals)) {
  is_num_like <- suppressWarnings(!is.na(as.integer(dt01$country_id)))
  if (any(is_num_like)) {
    idx <- suppressWarnings(as.integer(dt01$country_id))
    ok <- !is.na(idx) & idx >= 1 & idx <= length(iso_vals)
    dt01[ok, country_id := iso_vals[idx[ok]]]
  }
}

# land flag
dt01[, land := as.integer(!is.na(country_id) & country_id != "" & country_id != "NA" &
                            country_id != "0" & country_id != "00")]
dt01 <- dt01[, .(lon, lat, country_id, land)]
setkey(dt01, lon, lat)

cat("0.1° CENTER cells in dt01:", nrow(dt01), "\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ 0.5° -> 25 subcells via CENTER targets (no window/nearest)############
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dt05 <- as.data.table(pop_pm_combined)
dt_base <- dt05[, .(lon05 = lon, lat05 = lat)]

offsets <- c(-0.2, -0.1, 0.0, 0.1, 0.2)

# build 25 center points per 0.5° cell (exact 0.1° centers)
offs <- CJ(dlon = offsets, dlat = offsets)
dt_exp <- data.table(
  lon05 = rep(dt_base$lon05, each = nrow(offs)),
  lat05 = rep(dt_base$lat05, each = nrow(offs)),
  lon01 = round(rep(dt_base$lon05, each = nrow(offs)) + rep(offs$dlon, times = nrow(dt_base)), 1),
  lat01 = round(rep(dt_base$lat05, each = nrow(offs)) + rep(offs$dlat, times = nrow(dt_base)), 1)
)

# join to 0.1° center grid
dt_exp <- dt01[dt_exp, on = .(lon = lon01, lat = lat01)]

# ---- Decision rule helpers (MUST be defined before dt_rule) ----
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
# decision rules
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
############Perform spatial join############
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
# This keeps final output in ISO3.
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
############Visual check: Plot mapped countries############
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
############Convert back to regular dataframe and clean############
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
############Country code check############
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
############Save final output############
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

write_csv(pop_pm_final, here("output", "pop_pm_with_countries_rev.csv"))

cat("\nCountry mapping complete!\n")
cat("Output saved to: pop_pm_with_countries.csv\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############Summary statistics############
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