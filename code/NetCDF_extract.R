#########################################
## Extract grid cell data from NetCDF ##
#########################################


# Remove all objects from the environment
rm(list = ls())

############ Packages #####################################################


library(here)
library(ncdf4)
library(tidyverse)

############ import #####################################################

# Import the grid lookup table
grid_lookup <- readRDS(here("output", "grid_country_lookup.rds"))

## BASELINE 
# PM2.5 all (fire plus other)
pm_2000 <- nc_open(here("input", "CESM_09x125_PM25_2000_Baseline.nc"))
print(pm_2000)
# PM2.5 no fire
pm_2000_nf <- nc_open(here("input", "CESM_09x125_PM25_2000_BaseLine_NoFire.nc"))
print(pm_2000_nf)

## 2050 RCP4.5 
# PM2.5 all (fire plus other)
pm_2050_45 <- nc_open(here("input", "CESM_09x125_PM25_2050_RCP45.nc"))
print(pm_2050_45)
# PM2.5 no fire
pm_2050_45_nf <- nc_open(here("input", "CESM_09x125_PM25_2050_RCP45_NoFire.nc"))
# PM2.5 human intervention
pm_2050_45_hi <- nc_open(here("input", "CESM_09x125_PM25_2050_RCP45_HI.nc"))

## 2050 RCP8.5
# PM2.5 all (fire plus other)
pm_2050_85 <- nc_open(here("input", "CESM_09x125_PM25_2050_RCP85.nc"))
# PM2.5 no fire
pm_2050_85_nf <- nc_open(here("input", "CESM_09x125_PM25_2050_RCP85_NoFire.nc"))
# PM2.5 human intervention
pm_2050_85_hi <- nc_open(here("input", "CESM_09x125_PM25_2050_RCP85_HI.nc"))

## 2100 RCP4.5 
# PM2.5 all (fire plus other)
pm_2100_45 <- nc_open(here("input", "CESM_09x125_PM25_2100_RCP45.nc"))
# PM2.5 no fire
pm_2100_45_nf <- nc_open(here("input", "CESM_09x125_PM25_2100_RCP45_NoFire.nc"))
# PM2.5 human intervention
pm_2100_45_hi <- nc_open(here("input", "CESM_09x125_PM25_2100_RCP45_HI.nc"))

## 2100 RCP8.5
# PM2.5 all (fire plus other)
pm_2100_85 <- nc_open(here("input", "CESM_09x125_PM25_2100_RCP85.nc"))
# PM2.5 no fire
pm_2100_85_nf <- nc_open(here("input", "CESM_09x125_PM25_2100_RCP85_NoFire.nc"))
# PM2.5 human intervention
pm_2100_85_hi <- nc_open(here("input", "CESM_09x125_PM25_2100_RCP85_HI.nc"))

############ extract data from NetCDFs #####################################################

### Extract all months of PM2.5 data 

## BASELINE 
pm_2000_all <- ncvar_get(pm_2000, "pm25")  # 3D array: [lon, lat, time]
pm_2000_nf_all <- ncvar_get(pm_2000_nf, "pm25")  # 3D array: [lon, lat, time]

## 2050 RCP4.5 
pm_2050_45_all <- ncvar_get(pm_2050_45, "pm25")  # 3D array: [lon, lat, time]
pm_2050_45_nf_all <- ncvar_get(pm_2050_45_nf, "pm25")  # 3D array: [lon, lat, time]
pm_2050_45_hi_all <- ncvar_get(pm_2050_45_hi, "pm25")  # 3D array: [lon, lat, time]

## 2050 RCP8.5 
pm_2050_85_all <- ncvar_get(pm_2050_85, "pm25")  # 3D array: [lon, lat, time]
pm_2050_85_nf_all <- ncvar_get(pm_2050_85_nf, "pm25")  # 3D array: [lon, lat, time]
pm_2050_85_hi_all <- ncvar_get(pm_2050_85_hi, "pm25")  # 3D array: [lon, lat, time]

## 2100 RCP4.5 
pm_2100_45_all <- ncvar_get(pm_2100_45, "pm25")  # 3D array: [lon, lat, time]
pm_2100_45_nf_all <- ncvar_get(pm_2100_45_nf, "pm25")  # 3D array: [lon, lat, time]
pm_2100_45_hi_all <- ncvar_get(pm_2100_45_hi, "pm25")  # 3D array: [lon, lat, time]

## 2100 RCP8.5
pm_2100_85_all <- ncvar_get(pm_2100_85, "pm25")  # 3D array: [lon, lat, time]
pm_2100_85_nf_all <- ncvar_get(pm_2100_85_nf, "pm25")  # 3D array: [lon, lat, time]
pm_2100_85_hi_all <- ncvar_get(pm_2100_85_hi, "pm25")  # 3D array: [lon, lat, time]

# Extract longitude. latitude, date (can extract from any NetCDF file)
lon <- ncvar_get(pm_2000, "lon")
length(lon)

lat <- ncvar_get(pm_2000, "lat")
length(lat)

dates <- ncvar_get(pm_2000, "date")  # Get dates

nc_close(pm_2000)
nc_close(pm_2000_nf)
nc_close(pm_2050_45)
nc_close(pm_2050_45_nf)
nc_close(pm_2050_45_hi)
nc_close(pm_2100_45)
nc_close(pm_2100_45_nf)
nc_close(pm_2100_45_hi)
nc_close(pm_2050_85)
nc_close(pm_2050_85_nf)
nc_close(pm_2050_85_hi)
nc_close(pm_2100_85)
nc_close(pm_2100_85_nf)
nc_close(pm_2100_85_hi)

############ convert to dataframe #####################################################

# Convert 3D array to long format dataframe

### example for 1 array 
## BASELINE 
# all PM2.5 (with fire)
pm_df_2000 <- expand.grid(  # create df with every combo of lat, lon and month = 192x288x12 
  lon_index = 1:dim(pm_2000_all)[1], 
  lat_index = 1:dim(pm_2000_all)[2],
  month = 1:12
) %>%
  mutate(
    lon = lon[lon_index], # Converts longitude indices to actual longitude coordinates
    lat = lat[lat_index], # Converts latitude indices to actual longitude coordinates
    pm_2000 = mapply(function(i, j, k) pm_2000_all[i, j, k], 
                  lon_index, lat_index, month) # Extracts the PM2.5 value for each combination
  ) %>%
  select(month, lon, lat, pm_2000) # select columns needed 

### loop to convert all arrays to df 

# Create a named list with the CORRECT array names (all should end in _all)
arrays_to_process <- list(
  # Year 2000
  pm_2000_all = "2000",
  pm_2000_nf_all = "2000_nf",
  
  # Year 2050 RCP 4.5
  pm_2050_45_all = "2050_45",
  pm_2050_45_nf_all = "2050_45_nf",
  pm_2050_45_hi_all = "2050_45_hi",
  
  # Year 2050 RCP 8.5
  pm_2050_85_all = "2050_85",
  pm_2050_85_nf_all = "2050_85_nf",
  pm_2050_85_hi_all = "2050_85_hi",
  
  # Year 2100 RCP 4.5
  pm_2100_45_all = "2100_45",
  pm_2100_45_nf_all = "2100_45_nf",
  pm_2100_45_hi_all = "2100_45_hi",
  
  # Year 2100 RCP 8.5
  pm_2100_85_all = "2100_85",
  pm_2100_85_nf_all = "2100_85_nf",
  pm_2100_85_hi_all = "2100_85_hi"
)

# Loop through each array in the list
for (array_name in names(arrays_to_process)) {
  
  # Check if the array exists before processing
  if (!exists(array_name)) {
    cat("WARNING: Array", array_name, "does not exist. Skipping.\n")
    next
  }
  
  # Get the actual array object from the environment
  pm_array <- get(array_name)
  
  # Check if array has 3 dimensions
  if (is.null(dim(pm_array)) || length(dim(pm_array)) != 3) {
    cat("WARNING: Array", array_name, "doesn't have 3 dimensions. Skipping.\n")
    next
  }
  
  # Get the desired column name suffix
  suffix <- arrays_to_process[[array_name]]
  
  # Create column name and data frame name
  col_name <- paste0("pm_", suffix)      # e.g., "pm_2000"
  df_name <- paste0("pm_df_", suffix)    # e.g., "pm_df_2000"
  
  # Process the 3D array into a tidy data frame
  result <- expand.grid(
    lon_index = 1:dim(pm_array)[1],   # 288 longitude points
    lat_index = 1:dim(pm_array)[2],   # 192 latitude points
    month = 1:12                       # 12 months
  ) %>%
    mutate(
      # Convert indices to actual coordinates
      lon = lon[lon_index],
      lat = lat[lat_index],
      # Extract PM2.5 value for each combination
      !!col_name := mapply(function(i, j, k) pm_array[i, j, k], 
                           lon_index, lat_index, month)
    ) %>%
    # Keep only final columns
    select(month, lon, lat, all_of(col_name))
  
  # Create the data frame in global environment
  assign(df_name, result)
  
  # Print success message
  cat("✓ Successfully processed:", df_name, 
      "(", nrow(result), "rows )\n")
}

# Join with grid_lookup to add country information
pm_country <- pm_df_2000 %>%
  left_join(grid_lookup %>% select(lon, lat, country_name, country_code_iso3),
            by = c("lon", "lat")) %>%
  left_join(pm_df_2000_nf %>% select(month, lon, lat, pm_2000_nf),
            by = c("month", "lon", "lat")) %>%
  left_join(pm_df_2050_45 %>% select(month, lon, lat, pm_2050_45), # 2050  RCP 4.5 
            by = c("month", "lon", "lat")) %>%
  left_join(pm_df_2050_45_nf %>% select(month, lon, lat, pm_2050_45_nf),
            by = c("month", "lon", "lat")) %>%
  left_join(pm_df_2050_45_hi %>% select(month, lon, lat, pm_2050_45_hi),
            by = c("month", "lon", "lat")) %>%
  left_join(pm_df_2050_85 %>% select(month, lon, lat, pm_2050_85), # 2050  RCP 8.5 
            by = c("month", "lon", "lat")) %>%
  left_join(pm_df_2050_85_nf %>% select(month, lon, lat, pm_2050_85_nf),
            by = c("month", "lon", "lat")) %>%
  left_join(pm_df_2050_85_hi %>% select(month, lon, lat, pm_2050_85_hi),
            by = c("month", "lon", "lat")) %>%
  left_join(pm_df_2100_45 %>% select(month, lon, lat, pm_2100_45), # 2100  RCP 4.5 
            by = c("month", "lon", "lat")) %>%
  left_join(pm_df_2100_45_nf %>% select(month, lon, lat, pm_2100_45_nf),
            by = c("month", "lon", "lat")) %>%
  left_join(pm_df_2100_45_hi %>% select(month, lon, lat, pm_2100_45_hi), 
            by = c("month", "lon", "lat")) %>%
  left_join(pm_df_2100_85 %>% select(month, lon, lat, pm_2100_85), # 2100  RCP 8.5 
            by = c("month", "lon", "lat")) %>%
  left_join(pm_df_2100_85_nf %>% select(month, lon, lat, pm_2100_85_nf),
            by = c("month", "lon", "lat")) %>%
  left_join(pm_df_2100_85_hi %>% select(month, lon, lat, pm_2100_85_hi),
            by = c("month", "lon", "lat"))

# Convert lon/lat to numeric
pm_country <- pm_country %>%
  mutate(
    lon = as.numeric(lon),
    lat = as.numeric(lat)
  )

# Create pm_fire column (difference between total PM2.5 and no-fire PM2.5)
pm_country <- pm_country %>%
  mutate(fpm_2000 = pm_2000 - pm_2000_nf)

pm_country <- pm_country %>%
  mutate(fpm_2050_45 = pm_2050_45 - pm_2050_45_nf)

pm_country <- pm_country %>%
  mutate(fpm_2050_85 = pm_2050_85 - pm_2050_85_nf)

pm_country <- pm_country %>%
  mutate(fpm_2100_45 = pm_2100_45 - pm_2100_45_nf)

pm_country <- pm_country %>%
  mutate(fpm_2100_85 = pm_2100_85 - pm_2100_85_nf)

# Preview the result
colnames(pm_country)
head(pm_country, 5)

# Summary by country and month -- averages across cells in country

### all PM2.5
pm_country_month_ave <- pm_country %>%
  filter(!is.na(country_code_iso3)) %>%
  group_by(country_code_iso3, country_name, month) %>%
  summarise(
    mean_pm_2000 = mean(pm_2000, na.rm = TRUE),
    mean_fpm_2000 = mean(fpm_2000, na.rm = TRUE),
    mean_pm_2050_45 = mean(pm_2050_45, na.rm = TRUE),
    mean_fpm_2050_45 = mean(fpm_2050_45, na.rm = TRUE),
    mean_pm_2050_85 = mean(pm_2050_85, na.rm = TRUE),
    mean_fpm_2050_85 = mean(fpm_2050_85, na.rm = TRUE),
    mean_pm_2100_45 = mean(pm_2100_45, na.rm = TRUE),
    mean_fpm_2100_45 = mean(fpm_2100_45, na.rm = TRUE),
    mean_pm_2100_85 = mean(pm_2100_85, na.rm = TRUE),
    mean_fpm_2100_85 = mean(fpm_2100_85, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  )

head(pm_country_month_ave, 5)

### Create annual average by averaging across all 12 months for each grid cell
pm_annual_ave <- pm_country %>%
  group_by(lon, lat, country_name, country_code_iso3) %>%
  summarise(
    pm_2000 = mean(pm_2000, na.rm = TRUE),
    fpm_2000 = mean(fpm_2000, na.rm = TRUE),
    pm_2050_45 = mean(pm_2050_45, na.rm = TRUE),
    fpm_2050_45 = mean(fpm_2050_45, na.rm = TRUE),
    pm_2050_85 = mean(pm_2050_85, na.rm = TRUE),
    fpm_2050_85 = mean(fpm_2050_85, na.rm = TRUE),
    pm_2100_45 = mean(pm_2100_45, na.rm = TRUE),
    fpm_2100_45 = mean(fpm_2100_45, na.rm = TRUE),
    pm_2100_85 = mean(pm_2100_85, na.rm = TRUE),
    fpm_2100_85 = mean(fpm_2100_85, na.rm = TRUE),
    n_months = n(),
    .groups = "drop"
  )


# identify FPM2.5 change within year between RCP
pm_annual_ave <- pm_annual_ave %>%
  mutate(fpm_2100_rcp_chg = fpm_2100_85 - fpm_2100_45)

pm_annual_ave <- pm_annual_ave %>%
  mutate(fpm_2050_rcp_chg = fpm_2050_85 - fpm_2050_45)

# identify FPM2.5 change w.r.t base year, for each RCP
pm_annual_ave <- pm_annual_ave %>%
  mutate(fpm_2050_45_base_chg = fpm_2050_45 - fpm_2000)

pm_annual_ave <- pm_annual_ave %>%
  mutate(fpm_2050_85_base_chg = fpm_2050_85 - fpm_2000)

pm_annual_ave <- pm_annual_ave %>%
  mutate(fpm_2100_45_base_chg = fpm_2100_45 - fpm_2000)

pm_annual_ave <- pm_annual_ave %>%
  mutate(fpm_2100_85_base_chg = fpm_2100_85 - fpm_2000)

### Create month 4-9 average (April to September)
pm_mon4to9_ave <- pm_country %>%
  filter(month >= 4 & month <= 9) %>%
  group_by(lon, lat, country_name, country_code_iso3) %>%
  summarise(
    pm_2000 = mean(pm_2000, na.rm = TRUE),
    fpm_2000 = mean(fpm_2000, na.rm = TRUE),
    pm_2050_45 = mean(pm_2050_45, na.rm = TRUE),
    fpm_2050_45 = mean(fpm_2050_45, na.rm = TRUE),
    pm_2050_85 = mean(pm_2050_85, na.rm = TRUE),
    fpm_2050_85 = mean(fpm_2050_85, na.rm = TRUE),
    pm_2100_45 = mean(pm_2100_45, na.rm = TRUE),
    fpm_2100_45 = mean(fpm_2100_45, na.rm = TRUE),
    pm_2100_85 = mean(pm_2100_85, na.rm = TRUE),
    fpm_2100_85 = mean(fpm_2100_85, na.rm = TRUE),
    n_months = n(),
    .groups = "drop"
  )

# Verify
head(pm_mon4to9_ave)
head(pm_annual_ave)
summary(pm_mon4to9_ave$n_months)  # Should be 6 for all rows


# Verify reduction
print(paste("Original observations:", nrow(pm_country)))
print(paste("Annual average observations:", nrow(pm_annual_ave)))
print(paste("Reduction factor:", nrow(pm_country) / nrow(pm_annual_ave)))

# Preview
head(pm_annual_ave, 6)

############ basic stats #####################################################

summary(pm_annual_ave$fpm_2050_rcp_chg) # In 2050, the difference in PM2.5 concentration between RCP8.5 and RCP4.5
summary(pm_annual_ave$fpm_2100_rcp_chg) # In 2100, the difference in PM2.5 concentration between RCP8.5 and RCP4.5
summary(pm_annual_ave$fpm_2050_45_base_chg) # For RCP4.5, the difference in PM2.5 concentration between 2050 and 2000
summary(pm_annual_ave$fpm_2050_85_base_chg) # For RCP8.5, the difference in PM2.5 concentration between 2050 and 2000
summary(pm_annual_ave$fpm_2100_45_base_chg) # For RCP4.5, the difference in PM2.5 concentration between 2100 and 2000
summary(pm_annual_ave$fpm_2100_85_base_chg) # For RCP8.5, the difference in PM2.5 concentration between 2100 and 2000

# Create summary table
summary_table <- tibble(
  Variable = c(
    "fpm_2050_rcp_chg",
    "fpm_2100_rcp_chg",
    "fpm_2050_45_base_chg",
    "fpm_2050_85_base_chg",
    "fpm_2100_45_base_chg",
    "fpm_2100_85_base_chg"
  ),
  Description = c(
    "2050: RCP8.5 - RCP4.5",
    "2100: RCP8.5 - RCP4.5",
    "RCP4.5: 2050 - 2000",
    "RCP8.5: 2050 - 2000",
    "RCP4.5: 2100 - 2000",
    "RCP8.5: 2100 - 2000"
  ),
  Min = c(
    min(pm_annual_ave$fpm_2050_rcp_chg, na.rm = TRUE),
    min(pm_annual_ave$fpm_2100_rcp_chg, na.rm = TRUE),
    min(pm_annual_ave$fpm_2050_45_base_chg, na.rm = TRUE),
    min(pm_annual_ave$fpm_2050_85_base_chg, na.rm = TRUE),
    min(pm_annual_ave$fpm_2100_45_base_chg, na.rm = TRUE),
    min(pm_annual_ave$fpm_2100_85_base_chg, na.rm = TRUE)
  ),
  Q1 = c(
    quantile(pm_annual_ave$fpm_2050_rcp_chg, 0.25, na.rm = TRUE),
    quantile(pm_annual_ave$fpm_2100_rcp_chg, 0.25, na.rm = TRUE),
    quantile(pm_annual_ave$fpm_2050_45_base_chg, 0.25, na.rm = TRUE),
    quantile(pm_annual_ave$fpm_2050_85_base_chg, 0.25, na.rm = TRUE),
    quantile(pm_annual_ave$fpm_2100_45_base_chg, 0.25, na.rm = TRUE),
    quantile(pm_annual_ave$fpm_2100_85_base_chg, 0.25, na.rm = TRUE)
  ),
  Median = c(
    median(pm_annual_ave$fpm_2050_rcp_chg, na.rm = TRUE),
    median(pm_annual_ave$fpm_2100_rcp_chg, na.rm = TRUE),
    median(pm_annual_ave$fpm_2050_45_base_chg, na.rm = TRUE),
    median(pm_annual_ave$fpm_2050_85_base_chg, na.rm = TRUE),
    median(pm_annual_ave$fpm_2100_45_base_chg, na.rm = TRUE),
    median(pm_annual_ave$fpm_2100_85_base_chg, na.rm = TRUE)
  ),
  Mean = c(
    mean(pm_annual_ave$fpm_2050_rcp_chg, na.rm = TRUE),
    mean(pm_annual_ave$fpm_2100_rcp_chg, na.rm = TRUE),
    mean(pm_annual_ave$fpm_2050_45_base_chg, na.rm = TRUE),
    mean(pm_annual_ave$fpm_2050_85_base_chg, na.rm = TRUE),
    mean(pm_annual_ave$fpm_2100_45_base_chg, na.rm = TRUE),
    mean(pm_annual_ave$fpm_2100_85_base_chg, na.rm = TRUE)
  ),
  Q3 = c(
    quantile(pm_annual_ave$fpm_2050_rcp_chg, 0.75, na.rm = TRUE),
    quantile(pm_annual_ave$fpm_2100_rcp_chg, 0.75, na.rm = TRUE),
    quantile(pm_annual_ave$fpm_2050_45_base_chg, 0.75, na.rm = TRUE),
    quantile(pm_annual_ave$fpm_2050_85_base_chg, 0.75, na.rm = TRUE),
    quantile(pm_annual_ave$fpm_2100_45_base_chg, 0.75, na.rm = TRUE),
    quantile(pm_annual_ave$fpm_2100_85_base_chg, 0.75, na.rm = TRUE)
  ),
  Max = c(
    max(pm_annual_ave$fpm_2050_rcp_chg, na.rm = TRUE),
    max(pm_annual_ave$fpm_2100_rcp_chg, na.rm = TRUE),
    max(pm_annual_ave$fpm_2050_45_base_chg, na.rm = TRUE),
    max(pm_annual_ave$fpm_2050_85_base_chg, na.rm = TRUE),
    max(pm_annual_ave$fpm_2100_45_base_chg, na.rm = TRUE),
    max(pm_annual_ave$fpm_2100_85_base_chg, na.rm = TRUE)
  )
)

# Print table
print(summary_table, n = Inf)



############ save outputs #####################################################

# Save annual average to output folder
write_csv(pm_annual_ave, here("output", "annual_ave_pm25.csv"))
write_csv(pm_mon4to9_ave, here("output", "month4to9_ave_pm25.csv"))

print(paste("File saved:", here("output", "pm_annual_ave_pm25.csv")))

