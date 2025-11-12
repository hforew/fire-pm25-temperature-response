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

# import and open NetCDF files

## BASELINE 
# PM2.5 all (fire plus other)
pm25_2000 <- nc_open(here("input", "CESM_09x125_PM25_2000_Baseline.nc"))
print(pm25_2000)
# PM2.5 no fire
pm25_2000_nf <- nc_open(here("input", "CESM_09x125_PM25_2000_BaseLine_NoFire.nc"))
print(pm25_2000_nf)

# 2. Extract all months of PM2.5 data
pm25_2000_all <- ncvar_get(pm25_2000, "pm25")  # 3D array: [lon, lat, time]
dates <- ncvar_get(pm25_2000, "date")  # Get dates

pm25_2000_nf_all <- ncvar_get(pm25_2000_nf, "pm25")  # 3D array: [lon, lat, time]

## 2050 RCP4.5 
# PM2.5 all (fire plus other)
pm25_2050_45 <- nc_open(here("input", "CESM_09x125_PM25_2050_RCP45.nc"))
print(pm25_2050_45)
# PM2.5 no fire
pm25_2050_45_nf <- nc_open(here("input", "CESM_09x125_PM25_2050_RCP45_NoFire.nc"))
print(pm25_2050_45_nf)

# 2. Extract all months of PM2.5 data
pm25_2000_all <- ncvar_get(pm25_2000, "pm25")  # 3D array: [lon, lat, time]
pm25_2000_nf_all <- ncvar_get(pm25_2000_nf, "pm25")  # 3D array: [lon, lat, time]

pm25_2050_45_all <- ncvar_get(pm25_2050_45, "pm25")  # 3D array: [lon, lat, time]
pm25_2050_45_nf_all <- ncvar_get(pm25_2050_45_nf, "pm25")  # 3D array: [lon, lat, time]


# Extract longitude. latitude, date (can extract from any NetCDF file)
lon <- ncvar_get(pm25_2000, "lon")
length(lon)

lat <- ncvar_get(pm25_2000, "lat")
length(lat)

dates <- ncvar_get(pm25_2000, "date")  # Get dates

nc_close(pm25_2000)
nc_close(pm25_2000_nf)
nc_close(pm25_2050_45)
nc_close(pm25_2050_45_nf)

############ convert to dataframe #####################################################

# 3. Convert 3D array to long format dataframe


## BASELINE 
# all PM2.5 (with fire)
pm25_df_2000 <- expand.grid(
  lon_index = 1:dim(pm25_2000_all)[1],
  lat_index = 1:dim(pm25_2000_all)[2],
  month = 1:12
) %>%
  mutate(
    lon = lon[lon_index],
    lat = lat[lat_index],
    pm25_2000 = mapply(function(i, j, k) pm25_2000_all[i, j, k], 
                  lon_index, lat_index, month)
  ) %>%
  select(month, lon, lat, pm25_2000)

# no fire 
pm25_df_2000_nf <- expand.grid(
  lon_index = 1:dim(pm25_2000_nf_all)[1],
  lat_index = 1:dim(pm25_2000_nf_all)[2],
  month = 1:12
) %>%
  mutate(
    lon = lon[lon_index],
    lat = lat[lat_index],
    pm25_2000_nf = mapply(function(i, j, k) pm25_2000_nf_all[i, j, k], 
                  lon_index, lat_index, month)
  ) %>%
  select(month, lon, lat, pm25_2000_nf)

## 2050 RCP4.5 

# all PM2.5 (with fire)
pm25_df_2050_45 <- expand.grid(
  lon_index = 1:dim(pm25_2050_45_all)[1],
  lat_index = 1:dim(pm25_2050_45_all)[2],
  month = 1:12
) %>%
  mutate(
    lon = lon[lon_index],
    lat = lat[lat_index],
    pm25_2050_45 = mapply(function(i, j, k) pm25_2050_45_all[i, j, k], 
                       lon_index, lat_index, month)
  ) %>%
  select(month, lon, lat, pm25_2050_45)

# no fire 
pm25_df_2050_45_nf <- expand.grid(
  lon_index = 1:dim(pm25_2050_45_nf_all)[1],
  lat_index = 1:dim(pm25_2050_45_nf_all)[2],
  month = 1:12
) %>%
  mutate(
    lon = lon[lon_index],
    lat = lat[lat_index],
    pm25_2050_45_nf = mapply(function(i, j, k) pm25_2050_45_nf_all[i, j, k], 
                          lon_index, lat_index, month)
  ) %>%
  select(month, lon, lat, pm25_2050_45_nf)


## 2050 RCP8.5 

## 2100 RCP4.5 

## 2100 RCP8.5 


# # Rename pm25 column in pm25_df_2000_nf
# pm25_df_2000_nf <- pm25_df_2000_nf %>%
#   rename(pm25_nf = pm25)

# 4. Join with grid_lookup to add country information
pm25_with_countries <- pm25_df_2000 %>%
  left_join(grid_lookup %>% select(lon, lat, country_name, country_code_iso3),
            by = c("lon", "lat")) %>%
  left_join(pm25_df_2000_nf %>% select(month, lon, lat, pm25_2000_nf),
            by = c("month", "lon", "lat")) %>%
  left_join(pm25_df_2050_45 %>% select(month, lon, lat, pm25_2050_45),
            by = c("month", "lon", "lat")) %>%
  left_join(pm25_df_2050_45_nf %>% select(month, lon, lat, pm25_2050_45_nf),
            by = c("month", "lon", "lat"))

# Create pm25_fire column (difference between total PM2.5 and no-fire PM2.5)
pm25_with_countries <- pm25_with_countries %>%
  mutate(fpm25_2000 = pm25_2000 - pm25_2000_nf)

pm25_with_countries <- pm25_with_countries %>%
  mutate(fpm25_2050_45 = pm25_2050_45 - pm25_2050_45_nf)

# Preview the result
colnames(pm25_with_countries)
head(pm25_with_countries, 20)

# Summary by country and month -- averages across cells in country

# all PM2.5
pm_country_month_ave <- pm25_with_countries %>%
  filter(!is.na(country_code_iso3)) %>%
  group_by(country_code_iso3, country_name, month) %>%
  summarise(
    mean_pm_2000 = mean(pm25_2000, na.rm = TRUE),
    mean_fpm_2000 = mean(fpm25_2000, na.rm = TRUE),
    mean_pm_2050_45 = mean(pm25_2050_45, na.rm = TRUE),
    mean_fpm_2050_45 = mean(fpm25_2050_45, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  )

head(pm_country_month_ave, 5)

# # fire PM2.5 - fpm
# country_month_fpm25 <- pm25_with_countries %>%
#   filter(!is.na(country_code_iso3)) %>%
#   group_by(country_code_iso3, country_name, month) %>%
#   summarise(
#     mean_pm25 = mean(fpm25, na.rm = TRUE),
#     median_pm25 = median(fpm25, na.rm = TRUE),
#     n_cells = n(),
#     .groups = "drop"
#   )

# head(country_month_fpm25, 20)



# Create annual average by averaging across all 12 months for each grid cell
pm_annual_ave <- pm25_with_countries %>%
  group_by(lon, lat, country_name, country_code_iso3) %>%
  summarise(
    pm25_annual = mean(pm25, na.rm = TRUE),
    fpm25_annual = mean(fpm25, na.rm = TRUE),
    n_months = n(),
    .groups = "drop"
  )

# Create month 4-9 average (April to September)
pm_month4to9_ave <- pm25_with_countries %>%
  filter(month >= 4 & month <= 9) %>%
  group_by(lon, lat, country_name, country_code_iso3) %>%
  summarise(
    pm25_m4to9 = mean(pm25, na.rm = TRUE),
    fpm25_m4to9 = mean(fpm25, na.rm = TRUE),
    n_months = n(),
    .groups = "drop"
  )

# Verify
head(pm_month4to9_ave)
summary(pm_month4to9_ave$n_months)  # Should be 6 for all rows


# Verify reduction
print(paste("Original observations:", nrow(pm25_with_countries)))
print(paste("Annual average observations:", nrow(pm_annual_ave)))
print(paste("Reduction factor:", nrow(pm25_with_countries) / nrow(pm_annual_ave)))

# Preview
head(pm_annual_ave, 20)

# Convert lon/lat to numeric before saving
pm_annual_ave <- pm_annual_ave %>%
  mutate(
    lon = as.numeric(lon),
    lat = as.numeric(lat)
  )

pm_month4to9_ave <- pm_month4to9_ave %>%
  mutate(
    lon = as.numeric(lon),
    lat = as.numeric(lat)
  )


############ save outputs #####################################################

# Save annual average to output folder
write_csv(pm_annual_ave, here("output", "annual_ave_pm25_2000.csv"))
write_csv(pm_month4to9_ave, here("output", "month4to9_ave_pm25_2000.csv"))


print(paste("File saved:", here("output", "pm_annual_ave_pm25_2000.csv")))
