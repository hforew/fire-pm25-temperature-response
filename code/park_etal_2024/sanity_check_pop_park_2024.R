# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ City Population Extraction from Regridded Data ################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Goal: Pull 2015 gridded population for five US cities from the regridded Park grid
#       and show their 0.5° cells on an interactive map.
#   1. Load pop_regrid_park_2024.csv (259,200 cells, global 0.5° grid).
#   2. For each city bbox, keep cells whose 0.5° footprint (center ± 0.25°) touches
#      the box — any edge contact counts — and sum park_2015_pop; print cell count + total.
#   3. Rebuild each kept cell as a square polygon, stack all cities, and render a
#      color-by-city Leaflet map with per-cell popups (centroid lon/lat + pop_2015).
# Input  : output/pop_regrid_park_2024.csv
# Target : park_2015_pop ; Cities: NYC, Denver, LA, Honolulu, Anchorage


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Packages #######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(data.table)   # Fast CSV reading and data manipulation
library(here)         # Relative file paths anchored at project root
library(leaflet)      # Interactive map rendering
library(sf)           # Convert grid cells to polygon geometries

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Load Data ######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dt <- fread(here("output", "pop_regrid_park_2024.csv"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Helper Function ################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
get_city_pop <- function(dt, lon_min, lon_max, lat_min, lat_max,
                         col = "park_2015_pop") {
  subset_dt <- dt[
    lon - 0.25 < lon_max & lon + 0.25 > lon_min &
      lat - 0.25 < lat_max & lat + 0.25 > lat_min
  ]
  list(
    n_cells = nrow(subset_dt),
    pop_sum = sum(subset_dt[[col]], na.rm = TRUE),
    data    = subset_dt
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ City Bounding Boxes ############################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Approximate administrative extents (decimal degrees, WGS84):
#
#   NYC       lon: [-74.26, -73.70]    lat: [40.48, 40.92]
#   Denver    lon: [-105.11, -104.60]  lat: [39.61, 39.91]
#   LA        lon: [-118.67, -118.15]  lat: [33.70, 34.34]
#   Honolulu  lon: [-157.95, -157.65]  lat: [21.25, 21.40]
#   Anchorage lon: [-150.10, -149.60]  lat: [61.10, 61.30]

nyc       <- get_city_pop(dt, -74.26,  -73.70,  40.48, 40.92)
denver    <- get_city_pop(dt, -105.11, -104.60, 39.61, 39.91)
la        <- get_city_pop(dt, -118.67, -118.15, 33.70, 34.34)
honolulu  <- get_city_pop(dt, -157.95, -157.65, 21.25, 21.40)
anchorage <- get_city_pop(dt, -150.10, -149.60, 61.10, 61.30)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Results ########################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("NYC       —", nyc$n_cells,       "cells, park_2015_pop =", nyc$pop_sum,       "\n")
cat("Denver    —", denver$n_cells,    "cells, park_2015_pop =", denver$pop_sum,    "\n")
cat("LA        —", la$n_cells,        "cells, park_2015_pop =", la$pop_sum,        "\n")
cat("Honolulu  —", honolulu$n_cells,  "cells, park_2015_pop =", honolulu$pop_sum,  "\n")
cat("Anchorage —", anchorage$n_cells, "cells, park_2015_pop =", anchorage$pop_sum, "\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Build Cell Polygons ############################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
make_cell_sf <- function(result, city_name) {
  dt_sub <- result$data
  polys <- lapply(1:nrow(dt_sub), function(i) {
    x <- dt_sub$lon[i]
    y <- dt_sub$lat[i]
    st_polygon(list(matrix(c(
      x - 0.25, y - 0.25,
      x + 0.25, y - 0.25,
      x + 0.25, y + 0.25,
      x - 0.25, y + 0.25,
      x - 0.25, y - 0.25
    ), ncol = 2, byrow = TRUE)))
  })
  st_sf(
    city     = city_name,
    pop_2015 = dt_sub$park_2015_pop,
    geometry = st_sfc(polys, crs = 4326)
  )
}

sf_nyc       <- make_cell_sf(nyc,       "NYC")
sf_denver    <- make_cell_sf(denver,    "Denver")
sf_la        <- make_cell_sf(la,        "LA")
sf_honolulu  <- make_cell_sf(honolulu,  "Honolulu")
sf_anchorage <- make_cell_sf(anchorage, "Anchorage")

sf_all <- rbind(sf_nyc, sf_denver, sf_la, sf_honolulu, sf_anchorage)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Interactive Map ################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pal <- colorFactor(
  palette = c("firebrick", "steelblue", "forestgreen", "darkorchid", "darkorange"),
  domain  = sf_all$city
)

leaflet(sf_all) |>
  addProviderTiles(providers$CartoDB.Positron) |>
  addPolygons(
    fillColor   = ~pal(city),
    fillOpacity = 0.5,
    color       = "white",
    weight      = 1,
    popup       = ~paste0(
      "<b>", city, "</b><br>",
      "lon: ", round(st_coordinates(st_centroid(geometry))[, 1], 2), "<br>",
      "lat: ", round(st_coordinates(st_centroid(geometry))[, 2], 2), "<br>",
      "pop_2015: ", format(round(pop_2015), big.mark = ",")
    )
  ) |>
  addLegend(
    pal      = pal,
    values   = ~city,
    title    = "City",
    position = "bottomright"
  )

## THE END
