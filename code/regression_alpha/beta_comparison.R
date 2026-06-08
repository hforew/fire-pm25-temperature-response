# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## FPM–GMT RELATIONSHIP COMPARISON: Compare beta^(c) (per-capita fire PM2.5 change per 1°C GMT)
##                        from multiple regression specifications
##
## Goal: For each country c and regressions specification, compare betas from the following
##    pooled FE regression across Pierce/Park/Zhao data
##    regression with Pierce data
##    regression with Park data
## Inputs:
##   fpm_gmt_regression_coefs_classic.csv
##   fpm_gmt_regression_coefs_jules.csv
##   fpm_gmt_regression_coefs_ssib4.csv
##   fpm_gmt_regression_coefs_FE_all.csv
##   fpm_gmt_regression_coefs_pierce.csv
## Outputs:
##   
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Remove all objects from the environment to start fresh
rm(list = ls())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Packages #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(here)       # for portable file paths relative to project root
library(tidyverse)  # for data manipulation (dplyr, tidyr, purrr) and ggplot2
library(tibble)     # for tidy data frames (part of tidyverse, loaded explicitly for clarity)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Import #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# betas from Park et al models (CLASSIC, JULES, SSiB4)
beta_classic <- read_csv(here("output", "fpm_gmt_regression_coefs_classic.csv"))
beta_jules <- read_csv(here("output", "fpm_gmt_regression_coefs_jules.csv"))
beta_ssib4 <- read_csv(here("output", "fpm_gmt_regression_coefs_ssib4.csv"))

head(beta_classic)
head(beta_jules)
head(beta_ssib4)

# beta from FE regression (Park/Pierce/Zhao data)
beta_fe <- read_csv(here("output", "fpm_gmt_regression_coefs_FE_all.csv"))
head(beta_fe)

# beta from Pirce et al 
beta_pierce <- read_csv(here("output", "fpm_gmt_regression_coefs_pierce.csv"))
head(beta_pierce)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Join and Compare betas ##############################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

beta_comparison <- beta_fe |>
  select(country_code_iso3, country_name, estimate_beta_c, p.value_beta_c) |>
  rename(beta_fe = estimate_beta_c, p_fe = p.value_beta_c) |>
  left_join(
    beta_pierce |>
      select(country_code_iso3, estimate_beta_c, p.value_beta_c) |>
      rename(beta_pierce = estimate_beta_c, p_pierce = p.value_beta_c),
    by = "country_code_iso3"
  ) |>
  left_join(
    beta_classic |>
      select(country_code_iso3, estimate_beta_c, p.value_beta_c) |>
      rename(beta_classic = estimate_beta_c, p_classic = p.value_beta_c),
    by = "country_code_iso3"
  ) |>
  left_join(
    beta_jules |>
      select(country_code_iso3, estimate_beta_c, p.value_beta_c) |>
      rename(beta_jules = estimate_beta_c, p_jules = p.value_beta_c),
    by = "country_code_iso3"
  ) |>
  left_join(
    beta_ssib4 |>
      select(country_code_iso3, estimate_beta_c, p.value_beta_c) |>
      rename(beta_ssib4 = estimate_beta_c, p_ssib4 = p.value_beta_c),
    by = "country_code_iso3"
  ) |>
  mutate(across(where(is.numeric), \(x) round(x, 4)))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Summary statistics by model #############################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

beta_summary <- beta_comparison |>
  select(country_code_iso3, starts_with("beta_"), starts_with("p_")) |>
  pivot_longer(
    cols = -country_code_iso3,
    names_to  = c(".value", "model"),
    names_pattern = "^(beta|p)_(.*)"
  ) |>
  group_by(model) |>
  summarise(
    pct_beta_pos  = round(mean(beta > 0,      na.rm = TRUE) * 100, 4),
    pct_beta_sig5 = round(mean(p    < 0.05,   na.rm = TRUE) * 100, 4),
    beta_min      = round(min(beta,            na.rm = TRUE),       4),
    beta_p25      = round(quantile(beta, 0.25, na.rm = TRUE),       4),
    beta_median   = round(median(beta,         na.rm = TRUE),       4),
    beta_p75      = round(quantile(beta, 0.75, na.rm = TRUE),       4),
    beta_max      = round(max(beta,            na.rm = TRUE),       4)
  ) |>
  mutate(model = factor(model,
    levels = c("fe", "pierce", "classic", "jules", "ssib4"))) |>
  arrange(model)

beta_summary

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Save output ##############################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

write_csv(beta_comparison, here("output", "fpm_gmt_beta_comparison.csv"))


# THE END
