# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Combined pop_pm Dataset check ##############################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Purpose:
#   Sanity-check the population variables in the combined pop_pm dataset by
#   aggregating them at three nested geographic scales — global, United States,
#   and selected cities — and visualizing which 0.5-degree cells fall into each
#   city's catchment area.
#
# What this script does:
#   1. Reads the combined population and PM dataset (pop_pm_combined_final.csv)
#      and identifies every column whose name contains "pop".
#   2. Builds 0.5-degree polygon cells from each (lon, lat) center so cells can
#      be spatially intersected with country / city geometries.
#   3. Global table: sums each population column across all cells in the
#      dataset (one row, "Global").
#   4. United States table: sums each population column across cells that
#      intersect the US country polygon from rnaturalearth (one row,
#      "United States"). Note: this includes Alaska, Hawaii, and US territories.
#   5. City table: for 12 selected cities worldwide, builds a 50 km buffer
#      around each city center and sums each population column across all
#      0.5-degree cells intersecting the buffer (one row per city).
#   6. Combined table: row-binds the global, US, and city tables, rescales all
#      population columns to units of 10,000, and renders an interactive
#      DT::datatable for quick comparison across scales.
#   7. Interactive leaflet map: overlays the 0.5-degree cells assigned to each
#      city (color-coded by city, with toggleable layers) plus city-center
#      markers, so the spatial footprint of each city aggregation is verifiable.
#
# Notes / caveats:
#   - "Any intersecting cell" rule is used everywhere — a 0.5-degree cell is
#     fully attributed to a region if it overlaps the region at all, so cells
#     straddling boundaries can be double-counted across overlapping regions.
#   - City buffers are computed in EPSG:3857 (Web Mercator), which introduces
#     mild distortion at high latitudes (e.g. Anchorage, Moscow).
#
# Output:
#   - Inline HTML widgets only
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())

library(data.table)
library(sf)
library(knitr)
library(leaflet)
library(rnaturalearth)
library(rnaturalearthdata)
library(htmltools)
library(DT)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Read Final Combined pop_pm Dataset ###
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pop_pm_final <- fread("C:/Ford_BA_FPM25/output/pop_pm_combined_final.csv")

head(pop_pm_final)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Population Tables and Interactive Cell Map ##################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sf_use_s2(TRUE)

# All columns containing "pop"
pop_cols <- grep("pop", names(pop_pm_final), value = TRUE)
stopifnot(length(pop_cols) > 0)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Helper: Build 0.5-degree Cell Polygons ######################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
make_cell_polygon <- function(lon, lat, half_res = 0.25) {
  st_polygon(list(matrix(
    c(lon - half_res, lat - half_res,
      lon + half_res, lat - half_res,
      lon + half_res, lat + half_res,
      lon - half_res, lat + half_res,
      lon - half_res, lat - half_res),
    ncol = 2, byrow = TRUE
  )))
}

cell_sf <- st_sf(
  pop_pm_final[, .(cell_id = .I, lon, lat)],
  geometry = st_sfc(
    Map(make_cell_polygon, pop_pm_final$lon, pop_pm_final$lat),
    crs = 4326
  )
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ 1. Global Population Table ##################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
global_tbl <- data.table(region = "Global")
for (v in pop_cols) {
  global_tbl[, (v) := sum(pop_pm_final[[v]], na.rm = TRUE)]
}

global_tbl_html <- DT::datatable(
  global_tbl,
  rownames = FALSE,
  caption = htmltools::tags$caption(
    style = "caption-side: top; text-align: left;",
    "Global sums for all columns containing 'pop'."
  ),
  options = list(
    scrollX = TRUE,
    pageLength = 5,
    dom = "tip"
  )
)

global_tbl_html

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ 2. United States Population Table ###########################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
us_sf <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
us_sf <- us_sf[us_sf$admin == "United States of America", ]

us_idx <- lengths(st_intersects(cell_sf, us_sf)) > 0
cat("\nNumber of cells intersecting the United States:", sum(us_idx), "\n")

us_tbl <- data.table(region = "United States")
for (v in pop_cols) {
  us_tbl[, (v) := sum(pop_pm_final[[v]][us_idx], na.rm = TRUE)]
}

us_tbl_html <- DT::datatable(
  us_tbl,
  rownames = FALSE,
  caption = htmltools::tags$caption(
    style = "caption-side: top; text-align: left;",
    "United States sums for all columns containing 'pop', using any intersecting 0.5-degree cell."
  ),
  options = list(
    scrollX = TRUE,
    pageLength = 5,
    dom = "tip"
  )
)

us_tbl_html

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ 3. City Population Table ####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# City areas are approximated by a 50 km buffer around each city center.
# Any 0.5-degree cell intersecting the buffer is included.

city_buffer_km <- 50

city_dt <- data.table(
  city = c("New York City", "Los Angeles", "Denver", "Anchorage", "Honolulu",
           "Auckland", "Shanghai", "Delhi", "Paris", "Moscow",
           "Johannesburg", "Rio de Janeiro"),
  lon = c(-74.0060, -118.2437, -104.9903, -149.9003, -157.8583,
          174.7633, 121.4737, 77.1025, 2.3522, 37.6173,
          28.0473, -43.1729),
  lat = c(40.7128, 34.0522, 39.7392, 61.2181, 21.3069,
          -36.8485, 31.2304, 28.7041, 48.8566, 55.7558,
          -26.2041, -22.9068)
)

city_sf <- st_as_sf(city_dt, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
city_buffer_sf <- st_transform(city_sf, 3857)
city_buffer_sf <- st_buffer(city_buffer_sf, dist = city_buffer_km * 1000)
city_buffer_sf <- st_transform(city_buffer_sf, 4326)

city_tbl_list <- vector("list", nrow(city_buffer_sf))
city_cell_list <- vector("list", nrow(city_buffer_sf))

for (i in seq_len(nrow(city_buffer_sf))) {
  this_idx <- lengths(st_intersects(cell_sf, city_buffer_sf[i, ])) > 0
  
  tmp <- data.table(region = city_buffer_sf$city[i])
  for (v in pop_cols) {
    tmp[, (v) := sum(pop_pm_final[[v]][this_idx], na.rm = TRUE)]
  }
  city_tbl_list[[i]] <- tmp
  
  city_cell_list[[i]] <- st_sf(
    city = city_buffer_sf$city[i],
    lon = cell_sf$lon[this_idx],
    lat = cell_sf$lat[this_idx],
    geometry = st_geometry(cell_sf[this_idx, ])
  )
}

city_tbl <- rbindlist(city_tbl_list, fill = TRUE)

city_tbl_html <- DT::datatable(
  city_tbl,
  rownames = FALSE,
  caption = htmltools::tags$caption(
    style = "caption-side: top; text-align: left;",
    paste0(
      "City-area sums for all columns containing 'pop'. ",
      "Each city area is defined as a ", city_buffer_km,
      " km buffer around the city center, using any intersecting 0.5-degree cell."
    )
  ),
  options = list(
    scrollX = TRUE,
    pageLength = 12
  )
)

city_tbl_html

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Combined HTML View #######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
combined_tbl <- rbindlist(
  list(global_tbl, us_tbl, city_tbl),
  fill = TRUE,
  use.names = TRUE
)

# Round all pop columns to units of 10,000
pop_cols_combined <- grep("pop", names(combined_tbl), value = TRUE)
for (v in pop_cols_combined) {
  combined_tbl[, (v) := round(get(v) / 10000, 0)]
}

combined_tbl_html <- DT::datatable(
  combined_tbl,
  rownames = FALSE,
  caption = htmltools::tags$caption(
    style = "caption-side: top; text-align: left;",
    "Combined population summary table: Global, United States, and Cities. Population columns are shown in units of 10,000."
  ),
  options = list(
    scrollX = TRUE,
    pageLength = 20,
    autoWidth = TRUE
  )
)

combined_tbl_html

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ 4. Interactive Map of Selected Cells ########################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
city_cells_sf <- do.call(rbind, city_cell_list)
city_cells_sf$popup <- paste0(
  "<b>", city_cells_sf$city, "</b><br>",
  "Cell center: (", city_cells_sf$lon, ", ", city_cells_sf$lat, ")"
)

pal <- colorFactor(
  palette = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e",
              "#e6ab02", "#a6761d", "#666666", "#1f78b4", "#b2df8a",
              "#fb9a99", "#cab2d6"),
  domain = city_dt$city
)

city_map <- leaflet() |>
  addProviderTiles(providers$CartoDB.Positron)

for (ct in city_dt$city) {
  city_map <- city_map |>
    addPolygons(
      data = city_cells_sf[city_cells_sf$city == ct, ],
      fillColor = pal(ct),
      fillOpacity = 0.45,
      color = pal(ct),
      weight = 1,
      group = ct,
      popup = ~popup
    )
}

city_map <- city_map |>
  addCircleMarkers(
    data = city_dt,
    lng = ~lon, lat = ~lat,
    radius = 4,
    stroke = FALSE,
    fillOpacity = 1,
    group = "City centers",
    popup = ~city
  ) |>
  addLayersControl(
    overlayGroups = c(city_dt$city, "City centers"),
    options = layersControlOptions(collapsed = FALSE)
  )

city_map
