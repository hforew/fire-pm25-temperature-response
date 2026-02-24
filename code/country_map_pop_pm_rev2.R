#note: run time ~ 5 min
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
########### Map data frame with latitude and longitude coordinates to countries ##########
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Packages #######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(here)
library(tidyverse)
library(sf)           
library(rnaturalearth)   # Country boundary data
library(data.table)
library(dplyr)
library(leaflet)

sf_use_s2(FALSE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Import #########################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pop_pm_combined <- read_csv(here("output", "pop_pm_combined.csv"))
nrow(pop_pm_combined)
colnames(pop_pm_combined)

cat("\n=== Starting country mapping ===\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Load country boundaries (Natural Earth) #########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
world <- ne_countries(scale = 10, returnclass = "sf")

world_subset <- world %>%
  select(iso_a3, iso_a2, name, geometry)

world_subset <- st_buffer(world_subset, dist = 0)
cat("Loaded country boundaries:", nrow(world_subset), "countries\n")

# Keep ISO lookup
iso2_to_iso3 <- world_subset %>%
  st_drop_geometry() %>%
  distinct(iso_a2, iso_a3, name) %>%
  filter(!is.na(iso_a2), iso_a2 != "", !is.na(iso_a3), iso_a3 != "")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ NEW Decision rule (NO 0.1° netCDF) #############################
############ 0.5° cell assigned by MAX overlap area with country polygon ####
############ land_any = TRUE (land) if overlap exists, else FALSE (ocean) ###
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n=== NEW: Assign countries by 0.5° cell overlap area (Natural Earth) ===\n")

# 0.5° cell polygon from center
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

# Base table with row_id (for stable joins)
dt05 <- pop_pm_combined %>%
  mutate(row_id = row_number())

# Build cell polygons (WGS84)
cell_polys <- st_sfc(
  mapply(make_cell_polygon, dt05$lon, dt05$lat, SIMPLIFY = FALSE),
  crs = 4326
)
cells_sf <- st_sf(
  row_id = dt05$row_id,
  lon = dt05$lon,
  lat = dt05$lat,
  geometry = cell_polys
)

# Use an equal-area CRS for correct area comparison
cells_ea <- st_transform(cells_sf, 6933)
world_ea <- st_transform(world_subset %>% select(iso_a3, name), 6933)

# st_intersection on 259,200 polygons can be heavy; chunking by latitude bands
lat_breaks <- seq(-90, 90, by = 30)
cells_ea$lat_band <- cut(cells_ea$lat, breaks = lat_breaks, include.lowest = TRUE)

best_list <- list()

bands <- levels(cells_ea$lat_band)
cat("Processing bands:", paste(bands, collapse = " | "), "\n")

for (b in bands) {
  cat("  -> band:", b, "\n")
  cells_b <- cells_ea[cells_ea$lat_band == b, ]
  
  # intersect + area
  inter <- suppressWarnings(st_intersection(cells_b, world_ea))
  if (nrow(inter) == 0) next
  
  inter$overlap_area <- as.numeric(st_area(inter))
  
  best_b <- inter %>%
    st_drop_geometry() %>%
    group_by(row_id) %>%
    slice_max(overlap_area, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(
      row_id,
      land_any = TRUE,
      iso_a3 = iso_a3,
      name = name
    )
  
  best_list[[b]] <- best_b
}

best <- bind_rows(best_list)

# Merge back: no overlap -> ocean
assigned_df <- dt05 %>%
  left_join(best, by = "row_id") %>%
  mutate(
    land_any = if_else(is.na(land_any), FALSE, land_any),
    # keep a “major” column name for downstream compatibility (now ISO3 directly)
    country_id_major = iso_a3,
    # mark rare land cells still missing iso as -99 (should be near-zero)
    iso_a3 = if_else(is.na(iso_a3) & land_any, "-99", iso_a3),
    name   = if_else(is.na(name)   & land_any, "Unmapped land (area-overlap method)", name)
  )

cat("land_any TRUE count:", sum(assigned_df$land_any, na.rm = TRUE), "\n")
cat("assigned country count:", sum(!is.na(assigned_df$iso_a3) & assigned_df$iso_a3 != "-99"), "\n")
cat("invalid (-99) count:", sum(assigned_df$iso_a3 == "-99", na.rm = TRUE), "\n")

# Recreate sf points for consistency with downstream workflow
pop_pm_with_countries <- assigned_df %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

# NOTE: we assign ISO3 directly by polygons.

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Visual check: Plot matched countries ###########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
matched_countries <- pop_pm_with_countries %>%
  st_drop_geometry() %>%
  filter(!is.na(iso_a3)) %>%
  pull(iso_a3) %>%
  unique()

cat("\nNumber of unique countries matched:", length(matched_countries), "\n")

world_matched <- world_subset %>%
  filter(iso_a3 %in% matched_countries)

ggplot() +
  geom_sf(data = world_subset, fill = "gray90", color = "white", linewidth = 0.1) +
  geom_sf(data = world_matched, fill = "steelblue", color = "white", linewidth = 0.1) +
  theme_minimal() +
  labs(
    title = "Countries with Mapped Grid Cells (Area-Overlap Assignment)",
    subtitle = paste0(length(matched_countries), " countries matched"),
    caption = "Blue = countries with mapped cells, Gray = no cells"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank()
  )

ggsave(
  here("images", "country_mapping_check_rev.png"),
  width = 12,
  height = 6,
  dpi = 300
)

cat("Map saved to: images/country_mapping_check_rev.png\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Convert back to regular dataframe and clean ####################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pop_pm_final <- pop_pm_with_countries %>%
  st_drop_geometry() %>%
  rename(
    country_code_iso3 = iso_a3,
    country_name = name
  )

# Country code check (-99)
pop_pm_final %>%
  filter(country_code_iso3 == "-99") %>%
  count(country_name, country_code_iso3) %>%
  arrange(country_name) %>%
  print(n = 50)

# Keep manual fixes
pop_pm_final <- pop_pm_final %>%
  mutate(
    country_code_iso3 = case_when(
      country_name == "France" & country_code_iso3 == "-99" ~ "FRA",
      country_name == "Norway" & country_code_iso3 == "-99" ~ "NOR",
      TRUE ~ country_code_iso3
    )
  )

# Drop columns
pop_pm_final <- pop_pm_final %>%
  select(-any_of(c("iso_a2", "country_id_major")))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Save final output ##############################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write_csv(pop_pm_final, here("output", "pop_pm_with_countries_rev.csv"))

cat("\nCountry mapping complete!\n")
cat("Output saved to: output/pop_pm_with_countries_rev.csv\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Summary statistics #############################################
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
############ Interactive map ################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

is_invalid <- function(x) {
  !is.na(x) & x == "-99"
}
is_unmapped <- function(x) {
  is.na(x) | x == ""
}

na_iso3 <- c("USA", "CAN", "MEX")
mapped_df <- df %>%
  filter(!is_unmapped(country_code_iso3),
         !is_invalid(country_code_iso3),
         country_code_iso3 %in% na_iso3)

invalid_df <- df %>%
  filter(country_code_iso3 == "-99",
         !is.na(land_any), land_any == TRUE)

unmapped_df <- df %>%
  filter(is_unmapped(country_code_iso3),
         is.na(land_any) | land_any == FALSE)

cat("Total cells:", nrow(df), "\n")
cat("Mapped (all non-NA):", sum(!is.na(df$country_code_iso3)), "\n")
cat("Unmapped (NA or blank):", nrow(unmapped_df), "\n")
cat("Invalid (-99):", nrow(invalid_df), "\n")
cat("Mapped (NorthAmerica only):", nrow(mapped_df), "\n")


mapped_polys <- mapply(make_cell_polygon, mapped_df$lon, mapped_df$lat, SIMPLIFY = FALSE)
mapped_cells <- st_sf(mapped_df, geometry = st_sfc(mapped_polys, crs = 4326))

invalid_polys <- mapply(make_cell_polygon, invalid_df$lon, invalid_df$lat, SIMPLIFY = FALSE)
invalid_cells <- st_sf(invalid_df, geometry = st_sfc(invalid_polys, crs = 4326))

unmapped_polys <- mapply(make_cell_polygon, unmapped_df$lon, unmapped_df$lat, SIMPLIFY = FALSE)
unmapped_cells <- st_sf(unmapped_df, geometry = st_sfc(unmapped_polys, crs = 4326))

land <- ne_countries(scale = 10, returnclass = "sf") %>%
  st_transform(4326)

pal_na <- colorFactor(
  palette = c("#1b7837", "#2166ac", "#b2182b"),
  domain = na_iso3
)

leaflet() %>%
  addProviderTiles("CartoDB.Positron") %>%
  fitBounds(lng1 = -180, lat1 = -60, lng2 = 180, lat2 = 85) %>%
  addPolylines(
    data = land,
    color = "red",
    weight = 1,
    opacity = 0.7,
    fill = FALSE,
    group = "Land boundary"
  ) %>%
  addPolygons(
    data = mapped_cells,
    color = ~pal_na(country_code_iso3),
    weight = 0.4,
    fillColor = ~pal_na(country_code_iso3),
    fillOpacity = 0.40,
    stroke = TRUE,
    group = "Mapped cells (NorthAmerica only)",
    popup = ~paste0(
      "<b>Mapped cell (NorthAmerica only)</b><br>",
      "ISO3: ", country_code_iso3, "<br>",
      "Country: ", country_name, "<br>",
      "Lon: ", round(lon, 2), "<br>",
      "Lat: ", round(lat, 2)
    )
  ) %>%
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
    title = "Mapped ISO3 (NorthAmerica only)",
    opacity = 1
  ) %>%
  addLayersControl(
    overlayGroups = c("Land boundary", "Mapped cells (NorthAmerica only)", "Invalid cells (-99)", "Unmapped cells (NA/blank)"),
    options = layersControlOptions(collapsed = FALSE)
  )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ THE END ########################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ THE END