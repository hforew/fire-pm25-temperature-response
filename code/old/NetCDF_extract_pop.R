#########################################
## Extract Population grid cell data from NetCDF ##
#########################################


# Remove all objects from the environment
rm(list = ls())

############ Packages #####################################################

library(ncdf4)
library(here)

############ import #####################################################

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
pop_df <- pop_df %>%
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
write_csv(pop_df, here("output", "pop_df.csv"))


