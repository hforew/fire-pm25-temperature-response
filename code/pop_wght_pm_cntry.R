# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## compute population weighted exposure by country ##
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Remove all objects from the environment
rm(list = ls())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Packages #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(here)
library(tidyverse)
library(tibble)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Import #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Import combined data
# pop_pm_country <- read_csv(here("output", "pop_pm_with_countries.csv"))
pop_pm_country <- read_csv(here("output", "pop_pm_country_death.csv"))
colnames(pop_pm_country)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Baseline population weights ########################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# pop_bar_ic = average population in grid cell i over baseline period 2001-2010
# Population is held fixed across all periods and scenarios 
pop_yr_cols <- paste0("pop_tot_", 2001:2010)

pop_pm_country <- pop_pm_country %>%
  mutate(pop_bar = rowMeans(across(all_of(pop_yr_cols)), na.rm = TRUE))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Fire PM2.5 exposure by country and scenario ########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Total exposure formula:
#   exposure_cps = sum_i( fPM_ics * pop_bar_ic )
#
# Per-capita exposure formula:
#   exposure_cps^per-capita = sum_i( fPM_ics * pop_bar_ic ) / pop_bar_c
#
# fpm_* columns are already decadal averages, so 1/10 * sum_t collapses to the column value

fpm_cols <- c(
  "fpm_2000",
  "fpm_2050_45",
  "fpm_2050_85",
  "fpm_2100_45",
  "fpm_2100_85"
)

pop_wght_pm_cntry <- pop_pm_country %>%
  group_by(country_code_iso3, country_name) %>%  # one group per country
  summarise(
    # Total baseline population in country c: pop_bar_c = sum_i( pop_bar_ic )
    pop_bar_c = sum(pop_bar, na.rm = TRUE),
    # Sum grid-cell population to country level for each year 2001-2010
    across(all_of(pop_yr_cols), ~ sum(.x, na.rm = TRUE)),
    # Total fire PM2.5 exposure in country c:
    # exposure_cps = sum_i( fPM_ics * pop_bar_ic )
    # i.e. weighted sum of grid-cell fire PM by baseline population
    across(
      all_of(fpm_cols),
      ~ sum(.x * pop_bar, na.rm = TRUE),
      .names = "exposure_{.col}"
    ),
    # Per-capita fire PM2.5 exposure in country c:
    # exposure_cps^per-capita = sum_i( fPM_ics * pop_bar_ic ) / pop_bar_c
    # denominator uses only cells where fPM is non-missing
    across(
      all_of(fpm_cols),
      ~ sum(.x * pop_bar, na.rm = TRUE) /
          sum(pop_bar[!is.na(.x)], na.rm = TRUE),
      .names = "exposure_percap_{.col}"
    ),
    .groups = "drop"
  )

# exposure_percap_fpm_2000 = X_c^baseline (used in Section 5 excess mortality formula)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Export #############################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

write_csv(pop_wght_pm_cntry, here("output", "pop_wght_pm_cntry.csv"))

glimpse(pop_wght_pm_cntry)



