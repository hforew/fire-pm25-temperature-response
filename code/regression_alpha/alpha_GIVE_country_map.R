# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## MAP ALPHA_C2 COEFFICIENTS TO GIVE COUNTRY LIST ##
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Maps country-level alpha_c2 regression coefficients (change in per-capita fire
# PM2.5 exposure per 1°C GMT increase) onto the country list used by the
# GIVE integrated assessment model. The goal is to produce a country-coverage-
# complete dataset suitable for use as a damage function input in GIVE.

rm(list = ls())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Packages #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(here)
library(tidyverse)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Import #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# GIVE country list: ISO3 code and country name for all countries recognised by GIVE.
# Used as the reference set to join alpha_c2 coefficients onto.
give_countries <- read_csv(here("input", "GIVE", "GIVE_countries.csv"))
colnames(give_countries)

alpha_coeff <- read_csv(here("output", "fpm_gmt_regression_coefs.csv"))
colnames(alpha_coeff)
head(alpha_coeff)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Join #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Left join: keep all GIVE countries, attach alpha coefficients where available.
# Join key: give_countries$ISO3 <-> alpha_coeff$country_code_iso3
give_alpha <- give_countries %>%
  left_join(alpha_coeff, by = c("ISO3" = "country_code_iso3"))

# Check which GIVE countries have no matching alpha coefficient.
missing_alpha <- give_alpha %>%
  filter(is.na(estimate_alpha_c2)) %>%
  select(ISO3, country)

print(missing_alpha)

# For countries with no regression match, populate coefficient columns with the
# global-average values from the "global" row in alpha_coeff.
global_row <- alpha_coeff %>% filter(country_code_iso3 == "global")

give_alpha <- give_alpha %>%
  mutate(across(
    starts_with("estimate_") | starts_with("std.error_") |
      starts_with("statistic_") | starts_with("p.value_"),
    ~ if_else(is.na(.x), global_row[[cur_column()]], .x)
  ))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Export #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Keep only GIVE country identifier and the key damage function parameter.
alpha_country_meta <- give_alpha %>%
  select(ISO3, estimate_alpha_c2)

write_csv(alpha_country_meta, here("output", "alpha_country_meta.csv"))
