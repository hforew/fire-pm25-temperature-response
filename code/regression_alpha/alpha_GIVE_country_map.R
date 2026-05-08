# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## MAP ALPHA_C2 COEFFICIENTS TO GIVE COUNTRY LIST ##
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Maps country-level alpha_c2 regression coefficients (change in per-capita fire
# PM2.5 exposure per 1°C GMT increase) onto the country list used by the
# GIVE integrated assessment model. The goal is to produce a country-coverage-
# complete dataset suitable for use as a damage function input in GIVE.
# outputs for regression with 1) pierce only, 2) pierce + park + zhao, 3) FE all

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

alpha_coeff_pierce <- read_csv(here("output", "fpm_gmt_regression_coefs_pierce.csv"))
colnames(alpha_coeff_pierce)
head(alpha_coeff_pierce)

alpha_coeff_final <- read_csv(here("output", "fpm_gmt_regression_coefs_final.csv"))
colnames(alpha_coeff_final)
head(alpha_coeff_final)

alpha_coeff_fe_all <- read_csv(here("output", "fpm_gmt_regression_coefs_FE_all.csv"))
colnames(alpha_coeff_fe_all)
head(alpha_coeff_fe_all)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Join #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# --- Pierce ---
give_alpha_pierce <- give_countries %>%
  left_join(alpha_coeff_pierce, by = c("ISO3" = "country_code_iso3"))

missing_alpha_pierce <- give_alpha_pierce %>%
  filter(is.na(estimate_alpha_c2)) %>%
  select(ISO3, country)

print(missing_alpha_pierce)

global_row_pierce <- alpha_coeff_pierce %>%
  filter(country_code_iso3 == "global")

give_alpha_pierce <- give_alpha_pierce %>%
  mutate(across(
    starts_with("estimate_") | starts_with("std.error_") |
      starts_with("statistic_") | starts_with("p.value_"),
    ~ if_else(is.na(.x), global_row_pierce[[cur_column()]], .x)
  ))

# --- Final (Pierce, Park, and Zhao) ---
give_alpha_final <- give_countries %>%
  left_join(alpha_coeff_final, by = c("ISO3" = "country_code_iso3"))

missing_alpha_final <- give_alpha_final %>%
  filter(is.na(estimate_alpha_c2)) %>%
  select(ISO3, country)

print(missing_alpha_final)

global_row_final <- alpha_coeff_final %>% filter(country_code_iso3 == "global")

give_alpha_final <- give_alpha_final %>%
  mutate(across(
    starts_with("estimate_") | starts_with("std.error_") |
      starts_with("statistic_") | starts_with("p.value_"),
    ~ if_else(is.na(.x), global_row_final[[cur_column()]], .x)
  ))

# --- FE All (Pierce, Park, and Zhao with fire-model fixed effects) ---
give_alpha_fe_all <- give_countries %>%
  left_join(alpha_coeff_fe_all, by = c("ISO3" = "country_code_iso3"))

missing_alpha_fe_all <- give_alpha_fe_all %>%
  filter(is.na(estimate_alpha_c2)) %>%
  select(ISO3, country)

print(missing_alpha_fe_all)

global_row_fe_all <- alpha_coeff_fe_all %>%
  filter(country_code_iso3 == "global")

give_alpha_fe_all <- give_alpha_fe_all %>%
  mutate(across(
    starts_with("estimate_") | starts_with("std.error_") |
      starts_with("statistic_") | starts_with("p.value_"),
    ~ if_else(is.na(.x), global_row_fe_all[[cur_column()]], .x)
  ))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Export #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Keep only GIVE country identifier and the key damage function parameter.
alpha_country_meta_pierce <- give_alpha_pierce %>%
  select(ISO3, estimate_alpha_c2, std.error_alpha_c2)

write_csv(
  alpha_country_meta_pierce,
  here("output", "alpha_country_meta_pierce.csv")
)

alpha_country_meta_final <- give_alpha_final %>%
  select(ISO3, estimate_alpha_c2, std.error_alpha_c2)

write_csv(
  alpha_country_meta_final,
  here("output", "alpha_country_meta_final.csv")
)

alpha_country_meta_fe_all <- give_alpha_fe_all %>%
  select(ISO3, estimate_alpha_c2, std.error_alpha_c2)

write_csv(
  alpha_country_meta_fe_all,
  here("output", "alpha_country_meta_fe_all.csv")
)

### THE END
