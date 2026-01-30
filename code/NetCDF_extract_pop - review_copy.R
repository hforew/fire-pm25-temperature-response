#########################################
## Extract Population grid cell data from NetCDF ##
#########################################


# Remove all objects from the environment
rm(list = ls())

############ Packages #####################################################
#############################################################################################reviewed start
library(tidyverse)
library(terra)
library(geosphere)
#############################################################################################review end
library(ncdf4)
library(here)

############ import #####################################################
#setwd("C:/Ford_BA_FPM25")
pop_ssp1 <- nc_open(here("input", "SSP1_for_RCP45_2006-2100_population_density_c160701.nc"))
print(pop_ssp1) # dimensions 720 x 360 x 95 --- longitude x latitude x time, for 95 years
  # pop density dimensions / attributes
    # 720 longitude cells (0.5° resolution: 360°/0.5° = 720)
    # 360 latitude cells (0.5° resolution: 180°/0.5° = 360)
    # 95 time steps (years 2006-2100)
    # Units: km^-2 = people per square kilometer (population density)
    # Data structure: 3D array where each [lon, lat] cell has 95 time values

# Extract coordinates
lon <- ncvar_get(pop_ssp1, "lon")
lat <- ncvar_get(pop_ssp1, "lat")

# Check global domain boundaries
ncvar_get(pop_ssp1, "EDGEE")  # Eastern edge: 180 degree
ncvar_get(pop_ssp1, "EDGEW")  # Western edge: -180
ncvar_get(pop_ssp1, "EDGEN")  # Northern edge: 90
ncvar_get(pop_ssp1, "EDGES")  # Southern edge: -90

# Verify uniform 0.5 degree spacing
diff(lon[1:5])  # Result: 0.5, 0.5, 0.5, 0.5 
diff(lat[1:5])  # Result: 0.5, 0.5, 0.5, 0.5 

# Confirm lon/lat are cell CENTERS (not edges)
c(lon[1], lon[720])    # Result: -179.75, 179.75 (values end in .75 = centers) 
c(lat[1], lat[360])    # Result: -89.75, 89.75 (values end in .75 = centers) 
# If these were edges, we'd see -180 to 180 and -90 to 90 instead
# Cell center at -179.75 degree means cell spans from -180 degree to -179.5 degree 

############ extract data from NetCDF #####################################################

# Extract the year variable to find 2009
years <- ncvar_get(pop_ssp1, "year")
year_2009_index <- which(years == 2009)

# Extract population density for 2009 only
pop_density_2009 <- ncvar_get(
  pop_ssp1,              # NetCDF file connection
  "hdm",                 # Variable name: population density
  start = c(1, 1, year_2009_index),  # Starting position: [lon_start, lat_start, time_start]
                                    # 1, 1 = begin at first lon and lat cell
                                    # year_2009_index = begin at 2009 time slice (4th year)
  count = c(-1, -1, 1)   # How many cells to read: [lon_count, lat_count, time_count]
                          # -1, -1 = read ALL longitude and latitude cells (720 × 360)
                            # 1 = read only ONE time step (just 2009)
)
# Result: 2D array [720, 360] containing population density for 2009 only

# Check for NAs in the extracted 2009 population density array
sum(is.na(pop_density_2009))      # Count of NA values
any(is.na(pop_density_2009))      # TRUE if any NAs exist, FALSE if none

# More detailed check
cat("Total cells in array:", length(pop_density_2009), "\n")
cat("Number of NAs:", sum(is.na(pop_density_2009)), "\n")
cat("Number of zeros:", sum(pop_density_2009 == 0, na.rm = TRUE), "\n")
cat("Number of non-zero values:", sum(pop_density_2009 > 0, na.rm = TRUE), "\n")

# Close the NetCDF file
nc_close(pop_ssp1)

############ convert to dataframe #####################################################

# Expand grid to get all lon/lat combinations
pop_df <- expand.grid(
  lon = lon,        # All 720 longitude values (-179.75 to 179.75)
  lat = lat         # All 360 latitude values (-89.75 to 89.75)
            # expand.grid creates all combinations: 720 × 360 = 259,200 rows
            # Each row represents one grid cell with its lon/lat coordinates
) %>%
  mutate(
    # Add population density values to each lon/lat combination
    pop_dens_2009 = as.vector(pop_density_2009),  
          # as.vector() flattens the 2D array [720, 360] into 1D vector [259,200]
          # Order matches expand.grid: lon varies fastest (fills column-wise)
          # Row 1 = lon[1], lat[1]; Row 2 = lon[2], lat[1]; etc. after all lon filled, next lat
    pop_year = 2009     # Add year column (same value for all rows)
  ) 

# Result: Long-format dataframe with 4 columns: lon, lat, pop_density, year
# Each row = one grid cell's population density at its coordinates

############ basic data check #####################################################

# Check the result
glimpse(pop_df)
head(pop_df)
# Check if there's non-zero population
sum(pop_df$pop_dens_2009 > 0)

# See populated cells (returns first 20 cells non zer0)
pop_df %>% 
  filter(pop_dens_2009 > 0) %>% 
  head(20)

# Summary statistics
summary(pop_df$pop_dens_2009)

# Check dimensions
cat("Total cells:", nrow(pop_df), "\n")
cat("Cells with population > 0:", sum(pop_df$pop_dens_2009 > 0), "\n")

############ convert pop. density to pop total #####################################################

# Calculate grid cell area and total population
pop_df1 <- pop_df %>%
  mutate(
    # Calculate grid cell area (km²) accounting for latitude
    # At 0.5 degree resolution:
    # - Longitude width: 0.5 degree × 111 km/degree = 55.5 km (constant)
    # - Latitude height: 0.5 degree × 111 km/degree = 55.5 km (constant)
    # - But longitude distance shrinks toward poles: multiply by cos(latitude)
    cell_area_km2 = (0.5 * 111) * (0.5 * 111 * cos(lat * pi/180)),
    
    # Calculate total population in grid cell
    # Units: (people/km^2) × (km^2) = people
    pop_tot_2009 = pop_dens_2009 * cell_area_km2
  )

#############################################################################################review start
#NOAA’s blog notes that 111km is a commonly used approximation rather than an exact value, 
#so using 111km introduces a small discrepancy.
#when using 111km, the area will smaller than actual area.
#https://oceanservice.noaa.gov/facts/latitude.html?utm_source=chatgpt.com

#When the Earth is approximated as a sphere: C = 2 * pi * r
#1 degree arc length = (2 * pi * r) / 360 = (pi * r) / 180
#Average radius of the Earth = R ≈ 6371 km
#1 degree arc length = (pi * 6371) / 180 = 111.138555

#method1
print(pop_df1)

#method2
#Advantages: More rigorous than 111*cos(lat)
#Disadvantages: Still assumes a spherical Earth (rely on average approximation radius of the Earth = R ≈ 6371 km), 
#but generally provides a much better approximation than 111km.
R <- 6371 # km
dlat <- 0.5
dlon <- 0.5

pop_df2 <- pop_df %>%
  mutate(
    phi1 = (lat - dlat/2) * pi/180,
    phi2 = (lat + dlat/2) * pi/180,
    dlam = dlon * pi/180,
    cell_area_km2 = (R^2) * dlam * (sin(phi2) - sin(phi1)),
    pop_tot_2009 = pop_dens_2009 * cell_area_km2
  ) %>%
  select(-phi1, -phi2, -dlam)

print(pop_df2)

#method3
#Advantages: Closest representation of the Earth's true shape, not a approximation.
#Disadvantages: Slow for millions of grid points.
dlat <- 0.5
dlon <- 0.5

cell_area_wgs84_km2 <- function(lon, lat) {
  poly <- cbind(
    c(lon - dlon/2, lon + dlon/2, lon + dlon/2, lon - dlon/2, lon - dlon/2),
    c(lat - dlat/2, lat - dlat/2, lat + dlat/2, lat + dlat/2, lat - dlat/2)
  )
  geosphere::areaPolygon(poly) / 1e6  # m² -> km²
}

pop_df3 <- pop_df %>%
  rowwise() %>%
  mutate(
    cell_area_km2 = cell_area_wgs84_km2(lon, lat),
    pop_tot_2009 = pop_dens_2009 * cell_area_km2
  ) %>%
  ungroup()
print(pop_df3)

#summary statistics for each method
summary(pop_df1[, c("cell_area_km2", "pop_dens_2009", "pop_tot_2009")])
summary(pop_df2[, c("cell_area_km2", "pop_dens_2009", "pop_tot_2009")])
summary(pop_df3[, c("cell_area_km2", "pop_dens_2009", "pop_tot_2009")])


# Combine the three area computation methods into a single long-format table.
areas_long <- bind_rows(
  pop_df1 %>% select(lon, lat, cell_area_km2) %>% mutate(method = "111km"),
  pop_df2 %>% select(lon, lat, cell_area_km2) %>% mutate(method = "sphere"),
  pop_df3 %>% select(lon, lat, cell_area_km2) %>% mutate(method = "wgs84")
)

# For a regular latitude–longitude grid, cell area is theoretically
# invariant with respect to longitude. Therefore, we average cell areas
# across all latitudes for each longitude and method.
# The resulting plot should display (approximately) flat horizontal lines,
# confirming that grid cell area does not vary systematically with longitude.
lon_area <- areas_long %>%
  group_by(method, lon) %>%
  summarise(area = mean(cell_area_km2), .groups = "drop")

ggplot(lon_area, aes(x = lon, y = area, color = method)) +
  geom_line(linewidth = 1) +
  labs(
    title = "Grid cell area vs longitude",
    x = "Longitude (degrees)",
    y = "Cell area (km²)"
  ) +
  theme_minimal()

#log
ggplot(lon_area, aes(x = lon, y = area, color = method)) +
  geom_line(linewidth = 1) +
  scale_y_log10() +
  labs(
    title = "Grid cell area vs longitude (log scale)",
    x = "Longitude (degrees)",
    y = "Cell area (km², log scale)"
  ) +
  theme_minimal()

# In contrast, grid cell area varies systematically with latitude.
# Due to the convergence of meridians toward the poles, the longitudinal
# width of each grid cell shrinks proportionally to cos(latitude).
# Cell area is approximately proportional to cos(latitude),
# reaching its maximum at the equator and approaching zero near the poles.
lat_area <- areas_long %>%
  group_by(method, lat) %>%
  summarise(area = mean(cell_area_km2), .groups = "drop")

ggplot(lat_area, aes(x = lat, y = area, color = method)) +
  geom_line(linewidth = 1) +
  labs(
    title = "Grid cell area vs latitude",
    x = "Latitude (degrees)",
    y = "Cell area (km²)"
  ) +
  theme_minimal()

#log
ggplot(lat_area, aes(x = lat, y = area, color = method)) +
  geom_line(linewidth = 1) +
  scale_y_log10() +
  labs(
    title = "Grid cell area vs latitude (log scale)",
    x = "Latitude (degrees)",
    y = "Cell area (km², log scale)"
  ) +
  theme_minimal()

#Area is independent of longitude and varies only with latitude according to cos(lat); 
#the 111km method yields slightly smaller values 
#overall, while the spherical method and WGS84 are almost identical numerically, 
#WGS84 should has the highest accuracy intuitively.
#############################################################################################review end

# Verification checks
cat("Total global population (millions):", sum(pop_df$pop_tot_2009) / 1e6, "\n")
cat("Expected ~6,800 million for 2009\n")

# Check area variation by latitude
pop_df %>%
  mutate(lat_band = cut(lat, breaks = seq(-90, 90, 30))) %>%
  group_by(lat_band) %>%
  summarise(
    mean_area = mean(cell_area_km2),
    n_cells = n()
  )
# Result confirms latitude-dependent area calculation is correct:
# 1. Symmetry: N and S hemispheres identical (788, 2153, 2941 km^2 at same latitudes)
# 2. Polar convergence: Cells shrink toward poles (2941 km^2 equator -> 788 km^2 poles)
# 3. cos(latitude) effect: 
      # Mid-latitude (30-60 degree) = 73% of equatorial area (2153/2941 = 0.73 approx cos(45 degree) = 0.71)
# 4. Equal counts: 43,200 cells per band (720 lon x 60 lat steps = correct)


############ save outputs #####################################################

# Convert array columns to plain numeric vectors
pop_df <- pop_df %>%
  mutate(
    lon = as.numeric(lon),
    lat = as.numeric(lat),
    cell_area_km2 = as.numeric(cell_area_km2),
    pop_tot_2009 = as.numeric(pop_tot_2009)
  )

# Now save
# write_csv(pop_df, here("output", "pop_df.csv"))


