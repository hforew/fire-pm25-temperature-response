# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## MAP BETA_C COEFFICIENTS TO GIVE COUNTRY LIST ##
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Goal: Map country-level beta_c regression coefficients (change in per-capita
#   fire PM2.5 per 1°C GMT increase) onto the GIVE IAM country list, producing
#   a country-coverage-complete dataset for use as a damage function input.
# Input:
#   - input/GIVE/GIVE_countries.csv
#   - output/fpm_gmt_regression_coefs_pierce.csv     (pierce only)
#   - output/betas_fe_all/fpm_gmt_regression_coefs_FE_all.csv     (FE all)
#   - output/betas_fe_sensitivity/fpm_gmt_regression_coefs_FE_sensitivity_*.csv     (FE leave-one-fire-model-out sensitivity)
# Output:
#   - output/GIVE_betas/beta_country_meta_pierce.csv
#   - output/GIVE_betas/beta_country_meta_fe_all.csv
#   - output/GIVE_betas/beta_give_main.csv
#   - output/GIVE_betas/beta_give_excl_CESM.csv
#   - output/GIVE_betas/beta_give_excl_Zhao.csv
#   - output/GIVE_betas/beta_give_excl_classic.csv
#   - output/GIVE_betas/beta_give_excl_jules.csv
#   - output/GIVE_betas/beta_give_excl_ssib4.csv
# Execution order:
#   files run before: FE_sensitivity/fpm_gmt_regression_pierce.R   --> writes fpm_gmt_regression_coefs_pierce.csv
#                     FE/fpm_gmt_regression_FE_all.R             --> writes fpm_gmt_regression_coefs_FE_all.csv
#                     FE_sensitivity/fpm_gmt_regression_FE_cfact_sensitivity.R --> writes fpm_gmt_regression_coefs_FE_sensitivity_*.csv
#   files run after: NA

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
# Used as the reference set to join beta_c coefficients onto.
give_countries <- read_csv(here("input", "GIVE", "GIVE_countries.csv"))
colnames(give_countries)

# regression pierce only (beta for slope coef)
beta_coeff_pierce <- read_csv(here("output", "fpm_gmt_regression_coefs_pierce.csv"))
colnames(beta_coeff_pierce)
head(beta_coeff_pierce)

# regression w/ Pierce/Park/Zhao data WITH fixed effects (beta for slope coeff)
beta_coeff_fe_all <- read_csv(here("output", "betas_fe_all", "fpm_gmt_regression_coefs_FE_all.csv"))
colnames(beta_coeff_fe_all)
head(beta_coeff_fe_all)

# regression FE leave-one-fire-model-out sensitivity variants
# (full model + one variant per excluded fire-emissions product)
fe_sensitivity_files <- c(
  main         = "fpm_gmt_regression_coefs_FE_sensitivity_full.csv",
  excl_CESM    = "fpm_gmt_regression_coefs_FE_sensitivity_excl_CESM.csv",
  excl_Zhao    = "fpm_gmt_regression_coefs_FE_sensitivity_excl_Zhao.csv",
  excl_classic = "fpm_gmt_regression_coefs_FE_sensitivity_excl_classic.csv",
  excl_jules   = "fpm_gmt_regression_coefs_FE_sensitivity_excl_jules.csv",
  excl_ssib4   = "fpm_gmt_regression_coefs_FE_sensitivity_excl_ssib4.csv"
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Join #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# --- Pierce ---
give_beta_pierce <- give_countries %>%
  left_join(beta_coeff_pierce, by = c("ISO3" = "country_code_iso3"))

missing_beta_pierce <- give_beta_pierce %>%
  filter(is.na(estimate_beta_c)) %>%
  select(ISO3, country)

print(missing_beta_pierce)

global_row_pierce <- beta_coeff_pierce %>%
  filter(country_code_iso3 == "global")

give_beta_pierce <- give_beta_pierce %>%
  mutate(across(
    starts_with("estimate_") | starts_with("std.error_") |
      starts_with("statistic_") | starts_with("p.value_"),
    ~ if_else(is.na(.x), global_row_pierce[[cur_column()]], .x)
  ))

# --- FE All (Pierce, Park, and Zhao with fire-model fixed effects) ---
give_beta_fe_all <- give_countries %>%
  left_join(beta_coeff_fe_all, by = c("ISO3" = "country_code_iso3"))

missing_beta_fe_all <- give_beta_fe_all %>%
  filter(is.na(estimate_beta_c)) %>%
  select(ISO3, country)

print(missing_beta_fe_all)

global_row_fe_all <- beta_coeff_fe_all %>%
  filter(country_code_iso3 == "global")

give_beta_fe_all <- give_beta_fe_all %>%
  mutate(across(
    starts_with("estimate_") | starts_with("std.error_") |
      starts_with("statistic_") | starts_with("p.value_"),
    ~ if_else(is.na(.x), global_row_fe_all[[cur_column()]], .x)
  ))

# --- FE Sensitivity (leave-one-fire-model-out variants) ---
# Same join + global-fallback logic as above, applied once per variant.
give_beta_fe_sensitivity <- function(file_name) {
  # file_name <- fe_sensitivity_files[["main"]]  # DEBUG: uncomment to test line-by-line, re-comment before running 

  # Load one sensitivity variant's country-level coefficients.
  beta_coeff <- read_csv(here("output", "betas_fe_sensitivity", file_name))

  # Join onto the full GIVE country list; unmatched countries get NA.
  give_beta <- give_countries %>%
    left_join(beta_coeff, by = c("ISO3" = "country_code_iso3"))

  # Report GIVE countries with no country-specific estimate.
  missing_beta <- give_beta %>%
    filter(is.na(estimate_beta_c)) %>%
    select(ISO3, country)
  print(missing_beta)

  # Global-average row used to fill in for those missing countries.
  global_row <- beta_coeff %>%
    filter(country_code_iso3 == "global")

  # Fill all beta_c-related columns with the global estimate where NA.
  give_beta %>%
    mutate(across(
      starts_with("estimate_") | starts_with("std.error_") |
        starts_with("statistic_") | starts_with("p.value_"),
      ~ if_else(is.na(.x), global_row[[cur_column()]], .x)
    ))
}

give_beta_fe_sensitivity_list <- lapply(fe_sensitivity_files, give_beta_fe_sensitivity)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Export #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Keep only GIVE country identifier and the key damage function parameter.
beta_country_meta_pierce <- give_beta_pierce %>%
  select(ISO3, estimate_beta_c, std.error_beta_c)

write_csv(
  beta_country_meta_pierce,
  here("output", "GIVE_betas", "beta_country_meta_pierce.csv")
)

beta_country_meta_fe_all <- give_beta_fe_all %>%
  select(ISO3, estimate_beta_c, std.error_beta_c)

write_csv(
  beta_country_meta_fe_all,
  here("output", "GIVE_betas", "beta_country_meta_fe_all.csv")
)

# FE sensitivity variants: names(fe_sensitivity_files) --> beta_give_<variant>.csv
for (variant in names(give_beta_fe_sensitivity_list)) {
  beta_country_meta <- give_beta_fe_sensitivity_list[[variant]] %>%
    select(ISO3, estimate_beta_c, std.error_beta_c)

  write_csv(
    beta_country_meta,
    here("output", "GIVE_betas", paste0("beta_give_", variant, ".csv"))
  )
}

### THE END

