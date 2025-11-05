
# 1. Check how these countries are named in rnaturalearth
world %>%
  filter(iso_a3 %in% c("FRA", "NOR", "RWA", "LUX")) %>%
  select(name, iso_a3, admin, sovereignt) %>%
  st_drop_geometry()

# 2. Check all available name columns in world data
world %>%
  filter(iso_a3 %in% missing_countries) %>%
  select(name, name_long, admin, formal_en, iso_a3) %>%
  st_drop_geometry()

# 3. Compare with your CSV
countries_list %>%
  filter(ISO3 %in% c("FRA", "NOR", "RWA", "LUX"))


# Check if France and Norway exist with different ISO codes
world %>%
  filter(grepl("France|Norway", name, ignore.case = TRUE)) %>%
  select(name, iso_a3, admin, sovereignt) %>%
  st_drop_geometry()



# Check if these countries exist in the world data at all
world %>%
  filter(iso_a3 %in% missing_countries | 
           admin %in% missing_countries |
           name %in% countries_list$country[countries_list$ISO3 %in% missing_countries]) %>%
  select(name, iso_a3, admin, pop_est) %>%
  st_drop_geometry() %>%
  arrange(name)

# Check the actual size of these countries
world %>%
  filter(admin %in% c("Luxembourg", "Rwanda", "Bahrain", "Malta", "Cyprus")) %>%
  mutate(area_km2 = as.numeric(st_area(geometry)) / 1e6) %>%
  select(name, iso_a3, area_km2) %>%
  st_drop_geometry()


# Check if Rwanda exists in countries_sf after filtering
countries_sf %>%
  filter(admin == "Rwanda" | name == "Rwanda" | iso_a3 == "RWA") %>%
  select(name, iso_a3, admin) %>%
  st_drop_geometry()

# Check if Rwanda is in your original countries_list
countries_list %>%
  filter(ISO3 == "RWA")

# Check the exact spelling in both datasets
world %>%
  filter(iso_a3 == "RWA") %>%
  select(name, admin, formal_en) %>%
  st_drop_geometry()

# Look at what got mapped vs what's in countries_sf
nrow(countries_sf)
length(unique(countries_list$ISO3))


# Check if Rwanda grid cells exist before the ISO fix
grid_country_mapping %>%
  st_drop_geometry() %>%
  filter(name == "Rwanda") %>%
  select(name, iso_a3, lon, lat) %>%
  head(10)

# Check Rwanda's bounding box
rwanda_bbox <- world %>%
  filter(iso_a3 == "RWA") %>%
  st_bbox()

print(rwanda_bbox)

# Check if any grid cells fall within Rwanda's bounding box
grid_in_rwanda_area <- grid_sf %>%
  filter(lon >= rwanda_bbox["xmin"], lon <= rwanda_bbox["xmax"],
         lat >= rwanda_bbox["ymin"], lat <= rwanda_bbox["ymax"])

nrow(grid_in_rwanda_area)

# Visualize the issue - check if grid cells actually intersect Rwanda
rwanda_geom <- world %>%
  filter(iso_a3 == "RWA")

# Test the spatial join for just Rwanda
test_join <- st_join(grid_in_rwanda_area, 
                     rwanda_geom %>% select(name, iso_a3),
                     join = st_intersects)  # Try st_intersects instead of st_within

test_join %>%
  st_drop_geometry() %>%
  filter(!is.na(name))

# Reopen the NetCDF file
pm25_data <- nc_open(here("input", "CESM_09x125_PM25_2000_Baseline.nc"))

# Extract coordinates
lon <- ncvar_get(pm25_data, "lon")
lat <- ncvar_get(pm25_data, "lat")

# Rwanda bounding box: lon 28.86-30.88, lat -2.81 to -1.06

# Check which grid cells are near Rwanda
rwanda_area_lons <- lon[lon >= 28 & lon <= 31]
rwanda_area_lats <- lat[lat >= -3 & lat <= 0]

print("Longitudes near Rwanda:")
print(rwanda_area_lons)

print("Latitudes near Rwanda:")
print(rwanda_area_lats)

# Check actual PM2.5 values for Rwanda area (month 7 = July)
pm25_july <- ncvar_get(pm25_data, "pm25")[, , 7]

# Get values for Rwanda area
rwanda_pm25 <- pm25_july[lon >= 28 & lon <= 31, lat >= -3 & lat <= 0]
print("PM2.5 values in Rwanda area:")
print(rwanda_pm25)

nc_close(pm25_data)
