# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
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
library(broom)      # for tidy() and glance() to extract regression coefficients cleanly

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Import #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Import country-level per-capita fire PM2.5 exposure (µg/m^3/person/year).
# Each row is a country; columns include per-capita exposure for each period × scenario:
#   exposure_percap_fpm_2000         --> baseline period (~2001–2010)
#   exposure_percap_fpm_2050_45/85   --> 2041–2050 under RCP4.5 / RCP8.5
#   exposure_percap_fpm_2100_45/85   --> 2091–2100 under RCP4.5 / RCP8.5

#pop_wght <- read_csv(here("output", "pop_wght_pm_cntry.csv"))
pop_wght <- read_csv(here("output", "pop_wght_pm_cntry_final.csv"))

# Drop GIVE countries with no corresponding data in the Pierce or Park PM gridded data.
# These rows have pop_bar_c == 0, indicating no population-weighted exposure could be computed.
pop_wght <- pop_wght %>%
  filter(pop_bar_c != 0)

# Import decadal mean GMT anomaly (°C relative to 1850–1900 pre-industrial baseline)
# for each period × scenario combination.
# Rows: "2006-2010" (baseline), "2041-2050", "2091-2100"
# Columns: mean_gmt_45 (RCP4.5), mean_gmt_85 (RCP8.5)
gmt_chg <- read_csv(here("output", "gmt_periods_pi.csv"))

# Import Park decade GMT values (°C relative to pre-industrial baseline)
# Rows: one per Park snapshot decade (1960s–2010s)
# Columns: park_year, decade, mean_gmt_pi
gmt_park <- read_csv(here("output", "gmt_park_decades.csv"))

# Import Zhao et al. ~2095 GMT anomalies (°C relative to pre-industrial baseline)
# Rows: SSP245, SSP585 --- these temp changes assigned to Zhao but sourced from MimiSSPs 
# Columns: scenario, mean_gmt_pi
gmt_zhao <- read_csv(here("output", "gmt_zhao_mimi.csv"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Extract scalar GMT values for each period × scenario #######################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Baseline GMT: average of RCP4.5 and RCP8.5 for 2006–2010.
# At this early period the two scenarios have not yet diverged, so we average them
# to get a single baseline T value to pair with the observed baseline exposure.
gmt_baseline <- mean(c(gmt_chg$mean_gmt_45[gmt_chg$period == "2006-2010"],
                       gmt_chg$mean_gmt_85[gmt_chg$period == "2006-2010"]))

# Future GMT values for each period × scenario (used as regressors T_ps)
gmt_2040s_45 <- gmt_chg$mean_gmt_45[gmt_chg$period == "2041-2050"]
gmt_2040s_85 <- gmt_chg$mean_gmt_85[gmt_chg$period == "2041-2050"]
gmt_2090s_45 <- gmt_chg$mean_gmt_45[gmt_chg$period == "2091-2100"]
gmt_2090s_85 <- gmt_chg$mean_gmt_85[gmt_chg$period == "2091-2100"]

# Zhao et al. ~2095 GMT scalars from gmt_zhao_mimi.csv
gmt_2090s_245 <- gmt_zhao$mean_gmt_pi[gmt_zhao$scenario == "SSP245"]
gmt_2090s_585 <- gmt_zhao$mean_gmt_pi[gmt_zhao$scenario == "SSP585"]

# Park decade GMT scalars — one per decade, indexed by park_year (reliable numeric).
# Mirrors the future GMT scalar pattern above: one named object per observation type.
gmt_1960s <- gmt_park$mean_gmt_pi[gmt_park$park_year == 1965]   # 1960s decade mean GMT
gmt_1970s <- gmt_park$mean_gmt_pi[gmt_park$park_year == 1975]   # 1970s decade mean GMT
gmt_1980s <- gmt_park$mean_gmt_pi[gmt_park$park_year == 1985]   # 1980s decade mean GMT
gmt_1990s <- gmt_park$mean_gmt_pi[gmt_park$park_year == 1995]   # 1990s decade mean GMT
gmt_2000s <- gmt_park$mean_gmt_pi[gmt_park$park_year == 2005]   # 2000s decade mean GMT
gmt_2010s <- gmt_park$mean_gmt_pi[gmt_park$park_year == 2015]   # 2010s decade mean GMT


# Print GMT values
# pierce
cat("GMT baseline (2006-2010 avg):", gmt_baseline, "\n")
cat("GMT 2040s RCP4.5:", gmt_2040s_45, "  RCP8.5:", gmt_2040s_85, "\n")
cat("GMT 2090s RCP4.5:", gmt_2090s_45, "  RCP8.5:", gmt_2090s_85, "\n")

# park 
cat("GMT 1960s:", gmt_1960s, " 1970s:", gmt_1970s, " 1980s:", gmt_1980s, "\n")
cat("GMT 1990s:", gmt_1990s, " 2000s:", gmt_2000s, " 2010s:", gmt_2010s, "\n")

# Zhao (MimiSSPs)
cat("GMT SSP245 2095-99:", gmt_2090s_245, "GMT SSP585 2095-99:", gmt_2090s_585, "\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Reshape to long format for regression ######################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# We need one row per (country, period, scenario) with columns:
#   country_code_iso3, country_name, T_ps (GMT), exposure_percap (µg/m^3/person/yr)
#
# This long format lets us run a single lm() call per country across all 7 data points.

# Step 1: Select only the identifier columns and the seven per-capita exposure columns.
#         Drop rows where all exposure values are NA (e.g., uninhabited territories).
reg_data_wide <- pop_wght %>%
  select(country_code_iso3,
         country_name,
         exposure_percap_fpm_2000,             # baseline period exposure
         exposure_percap_fpm_2050_45,          # 2040s RCP4.5 exposure
         exposure_percap_fpm_2050_85,          # 2040s RCP8.5 exposure
         exposure_percap_fpm_2100_45,          # 2090s RCP4.5 exposure
         exposure_percap_fpm_2100_85,          # 2090s RCP8.5 exposure
         exposure_percap_fpm_2095_SSP245_Zhao, # Zhao ~2095 SSP2-4.5
         exposure_percap_fpm_2095_SSP585_Zhao) # Zhao ~2095 SSP5-8.5

# Step 2: Pivot to long format so each row is one (country, period×scenario) observation.
#         The column name encodes which GMT value applies to that row.
reg_data_long <- reg_data_wide %>%
  pivot_longer(
    cols      = starts_with("exposure_percap_fpm_"),  # the seven exposure columns
    names_to  = "period_scenario",                    # new column holding the column name
    values_to = "exposure_percap"                     # new column holding the exposure value
  ) %>%
  # Step 3: Map each period×scenario label to its corresponding GMT value (T_ps).
  #         This is the regressor in the country-level linear regression.
  mutate(T_ps = case_when(
    period_scenario == "exposure_percap_fpm_2000"    ~ gmt_baseline,   # baseline: avg of 45/85
    period_scenario == "exposure_percap_fpm_2050_45" ~ gmt_2040s_45,   # 2040s RCP4.5
    period_scenario == "exposure_percap_fpm_2050_85" ~ gmt_2040s_85,   # 2040s RCP8.5
    period_scenario == "exposure_percap_fpm_2100_45"          ~ gmt_2090s_45,   # 2090s RCP4.5
    period_scenario == "exposure_percap_fpm_2100_85"          ~ gmt_2090s_85,   # 2090s RCP8.5
    period_scenario == "exposure_percap_fpm_2095_SSP245_Zhao" ~ gmt_2090s_245,  # Zhao SSP2-4.5
    period_scenario == "exposure_percap_fpm_2095_SSP585_Zhao" ~ gmt_2090s_585   # Zhao SSP5-8.5
  )) %>%
  # Drop rows with missing exposure (uninhabited territories or missing data).
  # Countries with NA exposure cannot contribute to the regression.
  filter(!is.na(exposure_percap)) %>%
  # Label each Pierce/Zhao row with its source model for the fixed effect.
  mutate(fire_model = if_else(grepl("Zhao", period_scenario), "Zhao", "CESM"))

# Inspect the reshaped data to confirm structure (5 rows per country, fire_model visible)
print(head(reg_data_long, 15))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Reshape Park data to long format ###########################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Mirrors the future reshape above. One row per (country, fire_model, decade).
# 18 exposure columns total: 3 fire models × 6 decades.

# Step 1: Select identifier columns and all 18 Park per-capita exposure columns.
park_data_wide <- pop_wght %>%
  select(country_code_iso3,
         country_name,
         exposure_percap_park_classic_1960s_fpm,   # classic fire model, 1960s
         exposure_percap_park_classic_1970s_fpm,   # classic fire model, 1970s
         exposure_percap_park_classic_1980s_fpm,   # classic fire model, 1980s
         exposure_percap_park_classic_1990s_fpm,   # classic fire model, 1990s
         exposure_percap_park_classic_2000s_fpm,   # classic fire model, 2000s
         exposure_percap_park_classic_2010s_fpm,   # classic fire model, 2010s
         exposure_percap_park_jules_1960s_fpm,     # JULES fire model, 1960s
         exposure_percap_park_jules_1970s_fpm,     # JULES fire model, 1970s
         exposure_percap_park_jules_1980s_fpm,     # JULES fire model, 1980s
         exposure_percap_park_jules_1990s_fpm,     # JULES fire model, 1990s
         exposure_percap_park_jules_2000s_fpm,     # JULES fire model, 2000s
         exposure_percap_park_jules_2010s_fpm,     # JULES fire model, 2010s
         exposure_percap_park_ssib4_1960s_fpm,     # SSIB4 fire model, 1960s
         exposure_percap_park_ssib4_1970s_fpm,     # SSIB4 fire model, 1970s
         exposure_percap_park_ssib4_1980s_fpm,     # SSIB4 fire model, 1980s
         exposure_percap_park_ssib4_1990s_fpm,     # SSIB4 fire model, 1990s
         exposure_percap_park_ssib4_2000s_fpm,     # SSIB4 fire model, 2000s
         exposure_percap_park_ssib4_2010s_fpm)     # SSIB4 fire model, 2010s

# Step 2: Pivot to long format — one row per (country, fire_model × decade).
#         The column name is retained as period_scenario to match reg_data_long.
park_data_long <- park_data_wide %>%
  pivot_longer(
    cols      = starts_with("exposure_percap_park_"),   # the 18 Park exposure columns
    names_to  = "period_scenario",                      # column name becomes the label
    values_to = "exposure_percap"                       # exposure value
  ) %>%
  # Step 3: Map each column label to its GMT value (T_ps).
  #         Each decade maps to the same GMT regardless of fire model.
  mutate(T_ps = case_when(
    period_scenario == "exposure_percap_park_classic_1960s_fpm" ~ gmt_1960s,
    period_scenario == "exposure_percap_park_classic_1970s_fpm" ~ gmt_1970s,
    period_scenario == "exposure_percap_park_classic_1980s_fpm" ~ gmt_1980s,
    period_scenario == "exposure_percap_park_classic_1990s_fpm" ~ gmt_1990s,
    period_scenario == "exposure_percap_park_classic_2000s_fpm" ~ gmt_2000s,
    period_scenario == "exposure_percap_park_classic_2010s_fpm" ~ gmt_2010s,
    period_scenario == "exposure_percap_park_jules_1960s_fpm"   ~ gmt_1960s,
    period_scenario == "exposure_percap_park_jules_1970s_fpm"   ~ gmt_1970s,
    period_scenario == "exposure_percap_park_jules_1980s_fpm"   ~ gmt_1980s,
    period_scenario == "exposure_percap_park_jules_1990s_fpm"   ~ gmt_1990s,
    period_scenario == "exposure_percap_park_jules_2000s_fpm"   ~ gmt_2000s,
    period_scenario == "exposure_percap_park_jules_2010s_fpm"   ~ gmt_2010s,
    period_scenario == "exposure_percap_park_ssib4_1960s_fpm"   ~ gmt_1960s,
    period_scenario == "exposure_percap_park_ssib4_1970s_fpm"   ~ gmt_1970s,
    period_scenario == "exposure_percap_park_ssib4_1980s_fpm"   ~ gmt_1980s,
    period_scenario == "exposure_percap_park_ssib4_1990s_fpm"   ~ gmt_1990s,
    period_scenario == "exposure_percap_park_ssib4_2000s_fpm"   ~ gmt_2000s,
    period_scenario == "exposure_percap_park_ssib4_2010s_fpm"   ~ gmt_2010s
  )) %>%
  # Step 4: Label each row with its fire model, derived from the column name.
  mutate(fire_model = case_when(
    grepl("classic", period_scenario) ~ "classic",
    grepl("jules",   period_scenario) ~ "jules",
    grepl("ssib4",   period_scenario) ~ "ssib4"
  )) %>%
  filter(!is.na(exposure_percap))

# Inspect the reshaped Park data (18 rows per country — 3 models × 6 decades)
print(head(park_data_long, 18))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Combine historical (Park) and future projection data ########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# One regression per country across all observations:
#   18 Park historical obs (3 fire models × 6 decades) + 5 Pierce obs (1 base + 4 future)
#   + 2 Zhao obs (SSP2-4.5, SSP5-8.5) = 25 total.
reg_data_combined <- bind_rows(park_data_long, reg_data_long) %>%
  # Set fire_model as factor with "classic" as reference level so FE dummies
  # for jules, ssib4, CESM, Zhao are all interpreted relative to classic.
  mutate(fire_model = factor(fire_model,
                             levels = c("classic", "jules", "ssib4", "CESM", "Zhao")))

# Sanity check: total rows and per-country observation counts.
# Every country should have exactly 25 rows (18 Park + 5 Pierce + 2 Zhao).
cat("\nCombined rows:", nrow(reg_data_combined), "\n")
print(reg_data_combined %>% count(country_name) %>% summary())
# Expected: Combined rows: 4375 | n Min=25, Mean=25, Max=25 (175 countries * 25 obs) -- 25 per country

# Transparency: show observation count by fire_model across all countries.
# Per country: classic=6, jules=6, ssib4=6, CESM=5, Zhao=2 (total=25).
cat("\nTotal observations by fire_model (all countries combined):\n")
print(reg_data_combined %>% count(fire_model))
# Expected: classic=1050 = 175 x 6, jules=1050, ssib4=1050, CESM=875, Zhao=350 (175 countries each)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Regression ##############################################################
#
# Pooled OLS across all (country, period, scenario, fire model) observations:
#   PM_bar_ctsm = alpha_c + alpha_m + beta_c * T_tsm + epsilon_ctsm
#
# Implemented as:
#   exposure_percap ~ 0 + country_code_iso3 + fire_model + country_code_iso3:T_ps
#
# - 0 suppresses the global intercept; country_code_iso3 gives one alpha_c per country
# - fire_model adds shared fire-model FEs (classic = reference; jules/ssib4/CESM/Zhao as dummies)
# - country_code_iso3:T_ps gives one beta_c per country
#
# Key difference from fpm_gmt_regression_FE_all.R (per-country regressions):
#   fire-model FEs (alpha_m) are constrained equal across all countries here,
#   whereas the per-country approach allows country-specific fire-model offsets.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pooled_model <- lm(exposure_percap ~ 0 + country_code_iso3 + fire_model + country_code_iso3:T_ps,
                   data = reg_data_combined)

pooled_summary <- summary(pooled_model)
cat("R-squared:", round(pooled_summary$r.squared, 4), "\n")
cat("Adj. R-squared:", round(pooled_summary$adj.r.squared, 4), "\n")
cat("Residual df:", pooled_summary$df[2], "\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Regression Results ##############################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pooled_tidy <- tidy(pooled_model)

# beta_c: country_code_iso3XXX:T_ps interaction terms
beta_pooled <- pooled_tidy %>%
  filter(grepl(":T_ps$", term)) %>%
  mutate(country_code_iso3 = str_remove(term, "^country_code_iso3") %>%
                             str_remove(":T_ps$")) %>%
  rename(estimate_beta_c  = estimate,
         std.error_beta_c = std.error,
         statistic_beta_c = statistic,
         p.value_beta_c   = p.value) %>%
  select(country_code_iso3, estimate_beta_c, std.error_beta_c, statistic_beta_c, p.value_beta_c)

beta_pooled <- beta_pooled %>%
  left_join(reg_data_combined %>% distinct(country_code_iso3, country_name),
            by = "country_code_iso3")

# Degrees of freedom: N_obs - (N_country FEs + 4 model FEs + N_country slopes)
n_countries <- length(unique(reg_data_combined$country_code_iso3))
df_pooled   <- nrow(reg_data_combined) - (n_countries + 4 + n_countries)
t_crit      <- qt(0.975, df_pooled)
cat("\nPooled model df:", df_pooled, " | t_critical (95% CI):", round(t_crit, 4), "\n")

beta_pooled <- beta_pooled %>%
  mutate(
    lower_beta_c = estimate_beta_c - t_crit * std.error_beta_c,
    upper_beta_c = estimate_beta_c + t_crit * std.error_beta_c,
    gamma_beta_c = 0.008 * estimate_beta_c
  ) %>%
  select(country_code_iso3, country_name, gamma_beta_c,
         lower_beta_c, estimate_beta_c, upper_beta_c,
         std.error_beta_c, statistic_beta_c, p.value_beta_c)

cat("\n--- Summary of beta_c (pooled, shared fire-model FEs) ---\n")
print(summary(beta_pooled$estimate_beta_c))
cat("Countries with positive beta_c:", sum(beta_pooled$estimate_beta_c > 0, na.rm = TRUE), "\n")
cat("Countries with negative beta_c:", sum(beta_pooled$estimate_beta_c < 0, na.rm = TRUE), "\n")

print(head(beta_pooled %>% select(country_code_iso3, estimate_beta_c, lower_beta_c, upper_beta_c), 10))

# Shared fire-model FEs (relative to classic) — applies equally to all countries
fe_fire_model <- pooled_tidy %>%
  filter(grepl("^fire_model", term)) %>%
  mutate(fire_model = str_remove(term, "^fire_model"))

cat("\n--- Shared fire-model fixed effects (relative to classic) ---\n")
print(fe_fire_model %>% select(fire_model, estimate, std.error, statistic, p.value))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Save output ##############################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

write_csv(beta_pooled, here("output", "fpm_gmt_regression_coefs_pooled_fe.csv"))
cat("\nSaved pooled regression coefficients to output/fpm_gmt_regression_coefs_pooled_fe.csv\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ R implementation of regression ##############################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# The formula ~ 0 + country_code_iso3 + fire_model + country_code_iso3:T_ps expands to:
#
#   For observation (c, t, s, m):
#   PM_bar = alpha_c                    (country FE — one indicator per country, no global intercept)
#            + delta_jules  * D_jules   (shared fire-model FEs relative to classic)
#            + delta_ssib4  * D_ssib4
#            + delta_CESM   * D_CESM
#            + delta_Zhao   * D_Zhao
#            + beta_c       * T_ps      (country-specific slope via country_code_iso3:T_ps)
#            + epsilon
#
# delta_m is constrained equal across all countries — the key restriction vs. per-country regressions.
# beta_c is extracted from the country_code_iso3XXX:T_ps interaction coefficients.
# DF = 4375 - (175 country FEs + 4 model FEs + 175 slopes) = 4021.

# THE END
