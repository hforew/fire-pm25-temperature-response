
# Remove all objects from the environment
rm(list = ls())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Packages #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(here)
library(tidyverse)
library(tibble)
library(dplyr)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Import #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# pop_pm_country <- read_csv(here("output", "pop_pm_with_countries_rev.csv")) #may use pop_pm_with_countries_rev.csv for latest version
pop_pm_country <- read_csv(here("output", "pop_pm_with_countries_park.csv")) #may use pop_pm_with_countries_rev.csv for latest version


# check that fPM exposure in AUS is:
# LOWEST in 2100 under RCP85
# HIGHEST in 2100 under RCP45
# Yes, this is obs for AUS, matching pop-weighted results later
pop_pm_country %>%
  filter(country_code_iso3 == "AUS") %>%
  summarise(
    across(
      c(fpm_2000, fpm_2050_45, fpm_2050_85, fpm_2100_45, fpm_2100_85),
      ~ sum(.x, na.rm = TRUE)
    )
  )

# check that fPM exposure in CHN is:
# HIGHEST in 2100 under RCP 85
# LOWEST in 2000 under baseline
# Yes, this is obs for CHN, matching pop-weighted results later
pop_pm_country %>%
  filter(country_code_iso3 == "CHN") %>%
  summarise(
    across(
      c(fpm_2000, fpm_2050_45, fpm_2050_85, fpm_2100_45, fpm_2100_85),
      ~ sum(.x, na.rm = TRUE)
    )
  )


