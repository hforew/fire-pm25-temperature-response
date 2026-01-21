###########################
## Basic data stats ##
###########################


# Remove all objects from the environment
rm(list = ls())

############ Packages #####################################################

library(here)
library(tidyverse)



############ Import #####################################################

# Import annual average PM2.5 data
annual_ave <- read_csv(here("output", "annual_ave_pm25.csv"))
month4to9_ave <- read_csv(here("output", "month4to9_ave_pm25.csv"))


# Preview the data
head(annual_ave)
str(annual_ave)

############ basic stats #####################################################

summary(annual_ave$fpm_2050_rcp_chg) # In 2050, the difference in PM2.5 concentration between RCP8.5 and RCP4.5
summary(annual_ave$fpm_2100_rcp_chg) # In 2100, the difference in PM2.5 concentration between RCP8.5 and RCP4.5
summary(annual_ave$fpm_2050_45_base_chg) # For RCP4.5, the difference in PM2.5 concentration between 2050 and 2000
summary(annual_ave$fpm_2050_85_base_chg) # For RCP8.5, the difference in PM2.5 concentration between 2050 and 2000
summary(annual_ave$fpm_2100_45_base_chg) # For RCP4.5, the difference in PM2.5 concentration between 2100 and 2000
summary(annual_ave$fpm_2100_85_base_chg) # For RCP8.5, the difference in PM2.5 concentration between 2100 and 2000

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
    min(annual_ave$fpm_2050_rcp_chg, na.rm = TRUE),
    min(annual_ave$fpm_2100_rcp_chg, na.rm = TRUE),
    min(annual_ave$fpm_2050_45_base_chg, na.rm = TRUE),
    min(annual_ave$fpm_2050_85_base_chg, na.rm = TRUE),
    min(annual_ave$fpm_2100_45_base_chg, na.rm = TRUE),
    min(annual_ave$fpm_2100_85_base_chg, na.rm = TRUE)
  ),
  Q1 = c(
    quantile(annual_ave$fpm_2050_rcp_chg, 0.25, na.rm = TRUE),
    quantile(annual_ave$fpm_2100_rcp_chg, 0.25, na.rm = TRUE),
    quantile(annual_ave$fpm_2050_45_base_chg, 0.25, na.rm = TRUE),
    quantile(annual_ave$fpm_2050_85_base_chg, 0.25, na.rm = TRUE),
    quantile(annual_ave$fpm_2100_45_base_chg, 0.25, na.rm = TRUE),
    quantile(annual_ave$fpm_2100_85_base_chg, 0.25, na.rm = TRUE)
  ),
  Median = c(
    median(annual_ave$fpm_2050_rcp_chg, na.rm = TRUE),
    median(annual_ave$fpm_2100_rcp_chg, na.rm = TRUE),
    median(annual_ave$fpm_2050_45_base_chg, na.rm = TRUE),
    median(annual_ave$fpm_2050_85_base_chg, na.rm = TRUE),
    median(annual_ave$fpm_2100_45_base_chg, na.rm = TRUE),
    median(annual_ave$fpm_2100_85_base_chg, na.rm = TRUE)
  ),
  Mean = c(
    mean(annual_ave$fpm_2050_rcp_chg, na.rm = TRUE),
    mean(annual_ave$fpm_2100_rcp_chg, na.rm = TRUE),
    mean(annual_ave$fpm_2050_45_base_chg, na.rm = TRUE),
    mean(annual_ave$fpm_2050_85_base_chg, na.rm = TRUE),
    mean(annual_ave$fpm_2100_45_base_chg, na.rm = TRUE),
    mean(annual_ave$fpm_2100_85_base_chg, na.rm = TRUE)
  ),
  Q3 = c(
    quantile(annual_ave$fpm_2050_rcp_chg, 0.75, na.rm = TRUE),
    quantile(annual_ave$fpm_2100_rcp_chg, 0.75, na.rm = TRUE),
    quantile(annual_ave$fpm_2050_45_base_chg, 0.75, na.rm = TRUE),
    quantile(annual_ave$fpm_2050_85_base_chg, 0.75, na.rm = TRUE),
    quantile(annual_ave$fpm_2100_45_base_chg, 0.75, na.rm = TRUE),
    quantile(annual_ave$fpm_2100_85_base_chg, 0.75, na.rm = TRUE)
  ),
  Max = c(
    max(annual_ave$fpm_2050_rcp_chg, na.rm = TRUE),
    max(annual_ave$fpm_2100_rcp_chg, na.rm = TRUE),
    max(annual_ave$fpm_2050_45_base_chg, na.rm = TRUE),
    max(annual_ave$fpm_2050_85_base_chg, na.rm = TRUE),
    max(annual_ave$fpm_2100_45_base_chg, na.rm = TRUE),
    max(annual_ave$fpm_2100_85_base_chg, na.rm = TRUE)
  )
)

# Print table
print(summary_table, n = Inf)

############ filter data select countries #####################################################

annual_ave_USA <- annual_ave %>%
  filter(country_code_iso3 == "USA")

month4to9_ave_USA <- month4to9_ave %>%
  filter(country_code_iso3 == "USA")


summary_table_usa <- tibble(
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
    min(annual_ave_USA$fpm_2050_rcp_chg, na.rm = TRUE),
    min(annual_ave_USA$fpm_2100_rcp_chg, na.rm = TRUE),
    min(annual_ave_USA$fpm_2050_45_base_chg, na.rm = TRUE),
    min(annual_ave_USA$fpm_2050_85_base_chg, na.rm = TRUE),
    min(annual_ave_USA$fpm_2100_45_base_chg, na.rm = TRUE),
    min(annual_ave_USA$fpm_2100_85_base_chg, na.rm = TRUE)
  ),
  Q1 = c(
    quantile(annual_ave_USA$fpm_2050_rcp_chg, 0.25, na.rm = TRUE),
    quantile(annual_ave_USA$fpm_2100_rcp_chg, 0.25, na.rm = TRUE),
    quantile(annual_ave_USA$fpm_2050_45_base_chg, 0.25, na.rm = TRUE),
    quantile(annual_ave_USA$fpm_2050_85_base_chg, 0.25, na.rm = TRUE),
    quantile(annual_ave_USA$fpm_2100_45_base_chg, 0.25, na.rm = TRUE),
    quantile(annual_ave_USA$fpm_2100_85_base_chg, 0.25, na.rm = TRUE)
  ),
  Median = c(
    median(annual_ave_USA$fpm_2050_rcp_chg, na.rm = TRUE),
    median(annual_ave_USA$fpm_2100_rcp_chg, na.rm = TRUE),
    median(annual_ave_USA$fpm_2050_45_base_chg, na.rm = TRUE),
    median(annual_ave_USA$fpm_2050_85_base_chg, na.rm = TRUE),
    median(annual_ave_USA$fpm_2100_45_base_chg, na.rm = TRUE),
    median(annual_ave_USA$fpm_2100_85_base_chg, na.rm = TRUE)
  ),
  Mean = c(
    mean(annual_ave_USA$fpm_2050_rcp_chg, na.rm = TRUE),
    mean(annual_ave_USA$fpm_2100_rcp_chg, na.rm = TRUE),
    mean(annual_ave_USA$fpm_2050_45_base_chg, na.rm = TRUE),
    mean(annual_ave_USA$fpm_2050_85_base_chg, na.rm = TRUE),
    mean(annual_ave_USA$fpm_2100_45_base_chg, na.rm = TRUE),
    mean(annual_ave_USA$fpm_2100_85_base_chg, na.rm = TRUE)
  ),
  Q3 = c(
    quantile(annual_ave_USA$fpm_2050_rcp_chg, 0.75, na.rm = TRUE),
    quantile(annual_ave_USA$fpm_2100_rcp_chg, 0.75, na.rm = TRUE),
    quantile(annual_ave_USA$fpm_2050_45_base_chg, 0.75, na.rm = TRUE),
    quantile(annual_ave_USA$fpm_2050_85_base_chg, 0.75, na.rm = TRUE),
    quantile(annual_ave_USA$fpm_2100_45_base_chg, 0.75, na.rm = TRUE),
    quantile(annual_ave_USA$fpm_2100_85_base_chg, 0.75, na.rm = TRUE)
  ),
  Max = c(
    max(annual_ave_USA$fpm_2050_rcp_chg, na.rm = TRUE),
    max(annual_ave_USA$fpm_2100_rcp_chg, na.rm = TRUE),
    max(annual_ave_USA$fpm_2050_45_base_chg, na.rm = TRUE),
    max(annual_ave_USA$fpm_2050_85_base_chg, na.rm = TRUE),
    max(annual_ave_USA$fpm_2100_45_base_chg, na.rm = TRUE),
    max(annual_ave_USA$fpm_2100_85_base_chg, na.rm = TRUE)
  )
)

# Print table
print(summary_table_usa, n = Inf)



