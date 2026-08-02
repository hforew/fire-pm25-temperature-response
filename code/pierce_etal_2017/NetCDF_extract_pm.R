# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## EXTRACT PM2.5 GRID CELL DATA FROM NETCDF: Pierce et al. (2017) CESM simulations
##
## Goal: Extract fire PM2.5 (fPM) at the grid-cell level from Pierce et al. NetCDF
##   files, convert to a tidy dataframe, and compute both monthly and annual-average
##   fPM (= total PM2.5 minus no-fire PM2.5) for each grid cell, year (2000, 2050,
##   2100), and RCP scenario (4.5, 8.5) — with and without human-influence adjustment.
##
## Inputs:
##   input/CESM_09x125_PM25_*.nc          (14 gridded NetCDF files: baseline, 2050/2100 x
##                                         RCP4.5/8.5, each with all-PM, no-fire, and
##                                         human-influence variants; var "pm25", dims lon x lat x month)
##
## Outputs:
##   output/pm_joined/pm_df_all_month_grid.csv  (monthly grid-cell PM/fPM, all variants,
##                                         pre-annual-averaging intermediate)
##   output/pm_joined/annual_ave_pm25.csv       (annual grid-cell average PM/fPM; primary
##                                         output — pm/fpm for 2000/2050/2100 x RCP4.5/8.5,
##                                         plus RCP and base-year change columns)
##
## Execution order:
##   files run before: NA
##   files run after: pm_to_pop_regrid.R --> reads annual_ave_pm25.csv
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Remove all objects from the environment
rm(list = ls())

############ Packages #####################################################


library(here)
library(ncdf4)
library(tidyverse)

############ import #####################################################

## BASELINE
# PM2.5 all (fire plus other)
pm_2000 <- nc_open(here("input", "CESM_09x125_PM25_2000_Baseline.nc"))
print(pm_2000)
# PM2.5 no fire
pm_2000_nf <- nc_open(here("input", "CESM_09x125_PM25_2000_BaseLine_NoFire.nc"))
print(pm_2000_nf) # dimensions 288x192x12 --- longitude x latitude x time, for 12 months

## 2050 RCP4.5 
# PM2.5 all (fire plus other)
pm_2050_45 <- nc_open(here("input", "CESM_09x125_PM25_2050_RCP45.nc"))
print(pm_2050_45)
# PM2.5 no fire
pm_2050_45_nf <- nc_open(here("input", "CESM_09x125_PM25_2050_RCP45_NoFire.nc"))
# PM2.5 human influence
pm_2050_45_hi <- nc_open(here("input", "CESM_09x125_PM25_2050_RCP45_HI.nc"))

## 2050 RCP8.5
# PM2.5 all (fire plus other)
pm_2050_85 <- nc_open(here("input", "CESM_09x125_PM25_2050_RCP85.nc"))
# PM2.5 no fire
pm_2050_85_nf <- nc_open(here("input", "CESM_09x125_PM25_2050_RCP85_NoFire.nc"))
# PM2.5 human influence
pm_2050_85_hi <- nc_open(here("input", "CESM_09x125_PM25_2050_RCP85_HI.nc"))

## 2100 RCP4.5 
# PM2.5 all (fire plus other)
pm_2100_45 <- nc_open(here("input", "CESM_09x125_PM25_2100_RCP45.nc"))
# PM2.5 no fire
pm_2100_45_nf <- nc_open(here("input", "CESM_09x125_PM25_2100_RCP45_NoFire.nc"))
# PM2.5 human influence
pm_2100_45_hi <- nc_open(here("input", "CESM_09x125_PM25_2100_RCP45_HI.nc"))

## 2100 RCP8.5
# PM2.5 all (fire plus other)
pm_2100_85 <- nc_open(here("input", "CESM_09x125_PM25_2100_RCP85.nc"))
# PM2.5 no fire
pm_2100_85_nf <- nc_open(here("input", "CESM_09x125_PM25_2100_RCP85_NoFire.nc"))
# PM2.5 human influence
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

pm_df_all <- pm_df_2000 %>%
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
pm_df_all <- pm_df_all %>%
  mutate(
    lon = as.numeric(lon),
    lat = as.numeric(lat)
  )

# Create fire pm column (difference between total PM2.5 and no-fire PM2.5)
# baseline 
pm_df_all <- pm_df_all %>%
  mutate(fpm_2000 = pm_2000 - pm_2000_nf)

# 2050
pm_df_all <- pm_df_all %>%
  mutate(fpm_2050_45 = pm_2050_45 - pm_2050_45_nf)

pm_df_all <- pm_df_all %>%
  mutate(fpm_2050_45_hi = pm_2050_45_hi - pm_2050_45_nf)

pm_df_all <- pm_df_all %>%
  mutate(fpm_2050_85 = pm_2050_85 - pm_2050_85_nf)

pm_df_all <- pm_df_all %>%
  mutate(fpm_2050_85_hi = pm_2050_85_hi - pm_2050_85_nf)

# 2100
pm_df_all <- pm_df_all %>%
  mutate(fpm_2100_45 = pm_2100_45 - pm_2100_45_nf)

pm_df_all <- pm_df_all %>%
  mutate(fpm_2100_45_hi = pm_2100_45_hi - pm_2100_45_nf)

pm_df_all <- pm_df_all %>%
  mutate(fpm_2100_85 = pm_2100_85 - pm_2100_85_nf)

pm_df_all <- pm_df_all %>%
  mutate(fpm_2100_85_hi = pm_2100_85_hi - pm_2100_85_nf)

# Preview the result
colnames(pm_df_all)
head(pm_df_all, 5)

# save pm_df_all as intermediate detailed data file
write_csv(pm_df_all, here("output", "pm_joined", "pm_df_all_month_grid.csv"))

### Create annual average by averaging across all 12 months for each grid cell
pm_annual_ave <- pm_df_all %>%
  group_by(lon, lat) %>%
  summarise(
    pm_2000 = mean(pm_2000, na.rm = TRUE), # baseline
    pm_2000_nf = mean(pm_2000_nf, na.rm = TRUE), # baseline no fire 
    fpm_2000 = mean(fpm_2000, na.rm = TRUE),  # baseline minus no fire
    pm_2050_45 = mean(pm_2050_45, na.rm = TRUE), # 2050 
    pm_2050_45_nf = mean(pm_2050_45_nf, na.rm = TRUE),  
    fpm_2050_45 = mean(fpm_2050_45, na.rm = TRUE),
    fpm_2050_45_hi = mean(fpm_2050_45_hi, na.rm = TRUE),
    pm_2050_85 = mean(pm_2050_85, na.rm = TRUE),
    pm_2050_85_nf = mean(pm_2050_85_nf, na.rm = TRUE),
    fpm_2050_85 = mean(fpm_2050_85, na.rm = TRUE),
    fpm_2050_85_hi = mean(fpm_2050_85_hi, na.rm = TRUE),
    pm_2100_45 = mean(pm_2100_45, na.rm = TRUE), # 2100
    pm_2100_45_nf = mean(pm_2100_45_nf, na.rm = TRUE), # 2100
    fpm_2100_45 = mean(fpm_2100_45, na.rm = TRUE),
    fpm_2100_45_hi = mean(fpm_2100_45_hi, na.rm = TRUE),
    pm_2100_85 = mean(pm_2100_85, na.rm = TRUE),
    pm_2100_85_nf = mean(pm_2100_85_nf, na.rm = TRUE),
    fpm_2100_85 = mean(fpm_2100_85, na.rm = TRUE),
    fpm_2100_85_hi = mean(fpm_2100_85_hi, na.rm = TRUE),
    n_months = n(),
    .groups = "drop"
  )


# identify FPM2.5 change within year between RCP
pm_annual_ave <- pm_annual_ave %>%
  mutate(fpm_2100_rcp_chg = fpm_2100_85 - fpm_2100_45)

pm_annual_ave <- pm_annual_ave %>% 
  mutate(fpm_2100_rcp_chg_hi = fpm_2100_85_hi - fpm_2100_45_hi) # human influence

pm_annual_ave <- pm_annual_ave %>%
  mutate(fpm_2050_rcp_chg = fpm_2050_85 - fpm_2050_45)

pm_annual_ave <- pm_annual_ave %>%
  mutate(fpm_2050_rcp_chg_hi = fpm_2050_85_hi - fpm_2050_45_hi) # human influence

# identify FPM2.5 change w.r.t base year, for each RCP

# without human influence
pm_annual_ave <- pm_annual_ave %>%
  mutate(fpm_2050_45_base_chg = fpm_2050_45 - fpm_2000)

pm_annual_ave <- pm_annual_ave %>%
  mutate(fpm_2050_85_base_chg = fpm_2050_85 - fpm_2000)

pm_annual_ave <- pm_annual_ave %>%
  mutate(fpm_2100_45_base_chg = fpm_2100_45 - fpm_2000)

pm_annual_ave <- pm_annual_ave %>%
  mutate(fpm_2100_85_base_chg = fpm_2100_85 - fpm_2000)


# with human influence --- allows USA comparison with Ford et al 2018
pm_annual_ave <- pm_annual_ave %>%
  mutate(fpm_2050_45_base_chg_hi = fpm_2050_45_hi - fpm_2000)

pm_annual_ave <- pm_annual_ave %>%
  mutate(fpm_2050_85_base_chg_hi = fpm_2050_85_hi - fpm_2000)

pm_annual_ave <- pm_annual_ave %>%
  mutate(fpm_2100_45_base_chg_hi = fpm_2100_45_hi - fpm_2000)

pm_annual_ave <- pm_annual_ave %>%
  mutate(fpm_2100_85_base_chg_hi = fpm_2100_85_hi - fpm_2000)

# Verify reduction
print(paste("Original observations:", nrow(pm_df_all)))
print(paste("Annual average observations:", nrow(pm_annual_ave)))
print(paste("Reduction factor:", nrow(pm_df_all) / nrow(pm_annual_ave)))

# Preview
head(pm_annual_ave, 6)

############ save outputs #####################################################

# Save annual average to output folder
write_csv(pm_annual_ave, here("output", "pm_joined", "annual_ave_pm25.csv"))

print(paste("File saved:", here("output", "pm_annual_ave_pm25.csv")))

