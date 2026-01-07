#########################################
## Other Aggregation of extracted NetCDF data ##
#########################################


# Remove all objects from the environment
rm(list = ls())

############ Packages #####################################################


library(here)
library(ncdf4)
library(tidyverse)

############ Import #####################################################

# Import annual average PM2.5 data
pm_country <- read_csv(here("output", "pm_country_month_grid.csv"))

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
summary(pm_mon4to9_ave$n_months)  # Should be 6 for all rows

write_csv(pm_mon4to9_ave, here("output", "month4to9_ave_pm25.csv"))

