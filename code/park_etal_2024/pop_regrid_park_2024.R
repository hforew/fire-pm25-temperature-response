# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Read and Process Population Grid Data #########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Remove all objects from the environment
rm(list = ls())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Packages #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(R.matlab)
library(here)

# Read the .mat file using here()
pop_data <- readMat(here("input", "Park_etal_2024", "pop_grid_6515.mat"))

# Extract the population array (360 lat x 720 lon x 6 years)
pop_array <- pop_data$pop.grid

# Define years (same as PM2.5 data)
years <- c(1965, 1975, 1985, 1995, 2005, 2015)

# Create latitude and longitude coordinates
# 360 rows = latitudes from -90 to 90 (0.5 degree resolution)
# 720 columns = longitudes from -180 to 180 (0.5 degree resolution)
lat <- seq(-89.75, 89.75, by = 0.5)  # 360 values, centered
lon <- seq(-179.75, 179.75, by = 0.5)  # 720 values, centered

# Create list to store data for each year
pop_list <- list()

for (i in 1:6) {
  # Extract population data for year i
  pop_matrix <- pop_array[, , i]
  
  # Create data frame with expand.grid (automatically creates lon-lat grid)
  # expand.grid creates combinations with first variable changing fastest
  df <- expand.grid(
    lon = lon,
    lat = lat,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  # Transpose and vectorize population matrix to match the grid order
  # pop_matrix is [lat, lon], we need to convert to vector in [lon, lat] order
  df$population <- as.vector(t(pop_matrix))
  
  # Add year
  df$year <- years[i]
  
  # Reorder columns: year, lon, lat, population
  df <- df %>%
    select(year, lon, lat, population)
  
  pop_list[[i]] <- df
}

# Combine all years into long format
pop_long <- bind_rows(pop_list)

# Convert to wide format with park_year_pop naming
pop_wide <- pop_long %>%
  mutate(year_col = paste0("park_", year, "_pop")) %>%
  select(-year) %>%
  pivot_wider(
    names_from = year_col,
    values_from = population
  )

# Verify the coordinate order matches your PM2.5 data
head(pop_wide, 30)

# Column sums
pop_wide %>% 
  select(starts_with("park_")) %>% 
  colSums()

# Save data
write.csv(pop_wide, here("output", "pop_regrid_park_2024.csv"), row.names = FALSE)

##THE END
