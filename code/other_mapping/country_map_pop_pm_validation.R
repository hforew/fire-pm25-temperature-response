# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
################ CLM Population vs Official Population Comparison ##############
#steps:
# Extract gridded population density from the CLM NetCDF dataset for the selected year (2010)
# and convert density to total population using grid-cell area.

# Assign each grid cell to a country using spatial intersection with national boundaries
# and aggregate grid-level population to obtain country-level totals.

# Compare the CLM-derived national populations with official World Bank population data
# for the same year and visualize over- and under-estimation on an interactive map.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Goal: Validate the CLM gridded population product against World Bank totals for 2010.
#   1. Read CLM 'hdm' (pop density) for compare_year; multiply by grid-cell area -> total pop.
#   2. Assign each 0.5° cell to a country by POINT-IN-POLYGON on the cell center
#      (st_intersects), not by area overlap; aggregate cell pop -> country totals.
#   3. Pull official pop (World Bank SP.POP.TOTL) live via API; merge, compute
#      diff & % diff, render a diverging green-white-red Leaflet choropleth.
# Inputs : CLM hdm NetCDF (Li 2018 SSP1 CMIP6, 0.5°) ; gridcell_area_0.5deg.nc ; World Bank API


rm(list = ls())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############################# Packages #########################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(ncdf4)
library(dplyr)
library(tidyr)
library(tibble)
library(sf)
library(rnaturalearth)
library(leaflet)
library(scales)
library(htmltools)
library(jsonlite)

sf_use_s2(FALSE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############################# INPUT DATA #######################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# pop file used clmforc.Li_2018_SSP1_CMIP6_hdm_0.5x0.5_AVHRR_simyr1850-2100_c181205.nc
# CLM population dataset
clm_nc_file <- "input/pierce_etal_2017/clmforc.Li_2018_SSP1_CMIP6_hdm_0.5x0.5_AVHRR_simyr1850-2100_c181205.nc"

# grid cell area dataset
area_nc_file <- "input/landmask_area/gridcell_area_0.5deg.nc"

# comparison year
compare_year <- 2010

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############################# READ CLM DATA ####################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

clm_nc <- nc_open(clm_nc_file)

lon   <- ncvar_get(clm_nc, "LONGXY")
lat   <- ncvar_get(clm_nc, "LATIXY")
years <- ncvar_get(clm_nc, "year")

year_index <- match(compare_year, years)

if (is.na(year_index)) {
  stop("compare_year not found in NetCDF year dimension.")
}

pop_slice <- ncvar_get(
  clm_nc,
  "hdm",
  start = c(1, 1, year_index),
  count = c(-1, -1, 1)
)

nc_close(clm_nc)

grid_df <- tibble(
  lon = as.vector(lon),
  lat = as.vector(lat),
  pop_density = as.vector(pop_slice)
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############################# READ AREA DATA ###################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

area_nc <- nc_open(area_nc_file)

area_m2 <- ncvar_get(area_nc, "area")

nc_close(area_nc)

if (!all(dim(area_m2) == dim(lon))) {
  stop("Area grid dimensions do not match population grid dimensions.")
}

area_km2 <- area_m2 / 1e6

area_df <- tibble(
  lon = as.vector(lon),
  lat = as.vector(lat),
  cell_area_km2 = as.vector(area_km2)
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
########################## COMPUTE GRID POPULATION #############################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pop_df <- grid_df %>%
  left_join(area_df, by = c("lon", "lat")) %>%
  mutate(
    pop_tot = pop_density * cell_area_km2
  )

cat("Global population from CLM in", compare_year, ":",
    round(sum(pop_df$pop_tot, na.rm = TRUE) / 1e6, 2), "million\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############################# COUNTRY SHAPES ###################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

world <- rnaturalearth::ne_countries(
  scale = 50,
  returnclass = "sf"
) %>%
  sf::st_make_valid() %>%
  dplyr::select(iso_a3, name_long, geometry) %>%
  dplyr::filter(!is.na(iso_a3), iso_a3 != "-99")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
######################## OFFICIAL POPULATION: API DIRECT #######################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

wb_url <- paste0(
  "https://api.worldbank.org/v2/country/all/indicator/SP.POP.TOTL",
  "?format=json&date=", compare_year, ":", compare_year,
  "&per_page=32500&page=1"
)

wb_raw <- jsonlite::fromJSON(wb_url)

if (length(wb_raw) < 2 || is.null(wb_raw[[2]])) {
  stop("World Bank API did not return usable data.")
}

official_pop <- as_tibble(wb_raw[[2]]) %>%
  filter(
    !is.na(value),
    !is.na(countryiso3code),
    countryiso3code != "",
    countryiso3code != "NA"
  ) %>%
  transmute(
    iso_a3 = countryiso3code,
    official_population = as.numeric(value)
  )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############################ GRID TO COUNTRY ###################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# convert grid to sf points
pop_sf <- pop_df %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

# spatial join: assign each grid-cell center to a country
idx <- st_intersects(pop_sf, world)

country_id <- sapply(idx, function(x) {
  if (length(x) >= 1) x[1] else NA_integer_
})

pop_sf$iso_a3 <- world$iso_a3[country_id]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
########################## AGGREGATE COUNTRY POP ###############################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

clm_country <- pop_sf %>%
  st_drop_geometry() %>%
  filter(!is.na(iso_a3)) %>%
  group_by(iso_a3) %>%
  summarise(
    clm_population = sum(pop_tot, na.rm = TRUE),
    .groups = "drop"
  )

# quick check
print(
  clm_country %>%
    arrange(desc(clm_population)) %>%
    head(10)
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############################ MERGE #############################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

compare_df <- world %>%
  left_join(clm_country, by = "iso_a3") %>%
  left_join(official_pop, by = "iso_a3") %>%
  mutate(
    diff_pop = clm_population - official_population,
    pct_diff = 100 * diff_pop / official_population
  )

# table check
compare_table <- compare_df %>%
  st_drop_geometry() %>%
  filter(!is.na(clm_population), !is.na(official_population)) %>%
  arrange(desc(pct_diff))

print(
  compare_table %>%
    dplyr::select(iso_a3, name_long, clm_population, official_population, diff_pop, pct_diff) %>%
    head(20)
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############################ INTERACTIVE MAP ###################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# cap extreme values for better color scaling
cap_val <- quantile(abs(compare_df$pct_diff), probs = 0.95, na.rm = TRUE)

compare_df <- compare_df %>%
  mutate(
    pct_diff_cap = pmax(pmin(pct_diff, cap_val), -cap_val)
  )

pal <- colorNumeric(
  palette = colorRampPalette(c("darkgreen", "white", "red"))(200),
  domain = compare_df$pct_diff_cap,
  na.color = "#d9d9d9"
)

labels <- sprintf(
  "<strong>%s</strong><br/>
   CLM population: %s<br/>
   Official population: %s<br/>
   Difference: %s<br/>
   Percent difference: %s",
  compare_df$name_long,
  ifelse(is.na(compare_df$clm_population), "NA", comma(round(compare_df$clm_population))),
  ifelse(is.na(compare_df$official_population), "NA", comma(round(compare_df$official_population))),
  ifelse(is.na(compare_df$diff_pop), "NA", comma(round(compare_df$diff_pop))),
  ifelse(is.na(compare_df$pct_diff), "NA", paste0(round(compare_df$pct_diff, 2), "%"))
) %>%
  lapply(htmltools::HTML)

leaflet(compare_df) %>%
  addProviderTiles("CartoDB.Positron") %>%
  addPolygons(
    fillColor = ~pal(pct_diff_cap),
    weight = 0.7,
    color = "#666666",
    fillOpacity = 0.8,
    popup = labels,
    highlightOptions = highlightOptions(
      weight = 1.2,
      color = "#000000",
      fillOpacity = 0.9,
      bringToFront = TRUE
    )
  ) %>%
  addLegend(
    pal = pal,
    values = ~pct_diff_cap,
    title = paste0("CLM vs Official Population (%) ", compare_year),
    position = "bottomright",
    labFormat = labelFormat(suffix = "%")
  )
