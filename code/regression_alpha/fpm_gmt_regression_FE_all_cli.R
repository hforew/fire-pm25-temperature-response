# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## FPM–GMT RELATIONSHIP: Estimate beta_c (per-capita fire PM2.5 change per 1°C GMT)
##                        with Park et al. historical data integrated
##
## Goal: For each country c, estimate the linear regression with fire-model fixed effects:
##   PM_bar_ctsm = alpha_cm + beta_c * T_ts + epsilon_ctsm
##
## where PM_bar_ctsm is population-weighted per-capita fire PM2.5 exposure (µg/m^3/yr)
## for country c, time period t, scenario s, and fire model m;
## T_ts is GMT anomaly relative to the 1850-1900 pre-industrial baseline (°C);
## beta_c is the key damage function parameter (change in population-weighted
## per-capita fire PM2.5 per 1°C GMT), and alpha_cm are fire-model-specific
## intercepts (fixed effects absorbing model-level mean differences in exposure levels).
##
## Five fire-model fixed effect levels:
##   classic (reference), jules, ssib4  --> Park et al. historical models
##   CESM                               --> Pierce et al. projections
##   Zhao                               --> Zhao et al. projections
##
## Data spans three sources combined:
##   - Pierce et al. projections: baseline (~2001–2010), 2041–2050, and 2091–2100
##     under RCP4.5 and RCP8.5 — 5 (period × scenario) observations per country.
##   - Park et al. historical: 1960s–2010s across 3 fire models (classic, JULES, SSIB4)
##     - 18 (decade × fire model) observations per country. 
##     - Park **climate change-attributable** output
##   - Zhao et al. projections: ~2095 under SSP2-4.5 and SSP5-8.5 — 2 observations per country.
##   Combined: up to 25 observations per country for the regression.
##
## Inputs:
##   output/pop_wght_pm_cntry_final.csv   (country-level per-capita fPM exposure; Pierce,
##                                         Park, and Zhao columns)
##   output/gmt_periods_pi.csv            (decadal mean GMT anomaly for Pierce periods ×
##                                         scenarios, relative to 1850–1900 PI baseline)
##   output/gmt_park_decades.csv          (GMT anomaly for each Park snapshot decade,
##                                         relative to 1850–1900 PI baseline, sourced from OWID)
##   output/gmt_zhao_mimi.csv             (GMT anomaly for Zhao ~2095 SSP245 and SSP585,
##                                         sourced from MimiSSPs)
##
## Outputs:
##   output/fpm_gmt_regression_coefs_FE_all_cli.csv   (beta_c per country with SE, CI, p-value,
##                                                      and gamma_beta_c; primary downstream input)
##   output/reg_data_combined_fe_all_cli.csv           (combined long-format regression input:
##                                                      25 obs per country across all three sources)
##
## Execution order: 
##   files run before: pop_wght_pm_cntry.R      --> writes pop_wght_pm_cntry_final.csv
##   files run after: regression_FE_stats_all.R --> reads both outputs above
##                    alpha_GIVE_country_map.R   --> reads fpm_gmt_regression_coefs_FE_all.csv
##                    beta_comparison_latex.R    --> reads fpm_gmt_regression_coefs_FE_all.csv
##                    plots_beta_lobf_fe_all_cli.R   --> re-runs this file
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

# Future GMT values for each period × scenario (used as regressors T_ts; code variable: T_ps)
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

# Step 1: Select identifier columns and all 18 Park per-capita cli exposure columns.
park_data_wide <- pop_wght %>%
  select(country_code_iso3,
         country_name,
         exposure_percap_park_classic_1960s_fpm_cli,   # classic fire model, 1960s
         exposure_percap_park_classic_1970s_fpm_cli,   # classic fire model, 1970s
         exposure_percap_park_classic_1980s_fpm_cli,   # classic fire model, 1980s
         exposure_percap_park_classic_1990s_fpm_cli,   # classic fire model, 1990s
         exposure_percap_park_classic_2000s_fpm_cli,   # classic fire model, 2000s
         exposure_percap_park_classic_2010s_fpm_cli,   # classic fire model, 2010s
         exposure_percap_park_jules_1960s_fpm_cli,     # JULES fire model, 1960s
         exposure_percap_park_jules_1970s_fpm_cli,     # JULES fire model, 1970s
         exposure_percap_park_jules_1980s_fpm_cli,     # JULES fire model, 1980s
         exposure_percap_park_jules_1990s_fpm_cli,     # JULES fire model, 1990s
         exposure_percap_park_jules_2000s_fpm_cli,     # JULES fire model, 2000s
         exposure_percap_park_jules_2010s_fpm_cli,     # JULES fire model, 2010s
         exposure_percap_park_ssib4_1960s_fpm_cli,     # SSIB4 fire model, 1960s
         exposure_percap_park_ssib4_1970s_fpm_cli,     # SSIB4 fire model, 1970s
         exposure_percap_park_ssib4_1980s_fpm_cli,     # SSIB4 fire model, 1980s
         exposure_percap_park_ssib4_1990s_fpm_cli,     # SSIB4 fire model, 1990s
         exposure_percap_park_ssib4_2000s_fpm_cli,     # SSIB4 fire model, 2000s
         exposure_percap_park_ssib4_2010s_fpm_cli)     # SSIB4 fire model, 2010s

# Step 2: Pivot to long format — one row per (country, fire_model × decade).
#         The column name is retained as period_scenario to match reg_data_long.
park_data_long <- park_data_wide %>%
  pivot_longer(
    cols      = ends_with("_fpm_cli"),          # the 18 Park cli exposure columns
    names_to  = "period_scenario",              # column name becomes the label
    values_to = "exposure_percap"               # exposure value
  ) %>%
  # Step 3: Map each column label to its GMT value (T_ps).
  #         Each decade maps to the same GMT regardless of fire model.
  mutate(T_ps = case_when(
    period_scenario == "exposure_percap_park_classic_1960s_fpm_cli" ~ gmt_1960s,
    period_scenario == "exposure_percap_park_classic_1970s_fpm_cli" ~ gmt_1970s,
    period_scenario == "exposure_percap_park_classic_1980s_fpm_cli" ~ gmt_1980s,
    period_scenario == "exposure_percap_park_classic_1990s_fpm_cli" ~ gmt_1990s,
    period_scenario == "exposure_percap_park_classic_2000s_fpm_cli" ~ gmt_2000s,
    period_scenario == "exposure_percap_park_classic_2010s_fpm_cli" ~ gmt_2010s,
    period_scenario == "exposure_percap_park_jules_1960s_fpm_cli"   ~ gmt_1960s,
    period_scenario == "exposure_percap_park_jules_1970s_fpm_cli"   ~ gmt_1970s,
    period_scenario == "exposure_percap_park_jules_1980s_fpm_cli"   ~ gmt_1980s,
    period_scenario == "exposure_percap_park_jules_1990s_fpm_cli"   ~ gmt_1990s,
    period_scenario == "exposure_percap_park_jules_2000s_fpm_cli"   ~ gmt_2000s,
    period_scenario == "exposure_percap_park_jules_2010s_fpm_cli"   ~ gmt_2010s,
    period_scenario == "exposure_percap_park_ssib4_1960s_fpm_cli"   ~ gmt_1960s,
    period_scenario == "exposure_percap_park_ssib4_1970s_fpm_cli"   ~ gmt_1970s,
    period_scenario == "exposure_percap_park_ssib4_1980s_fpm_cli"   ~ gmt_1980s,
    period_scenario == "exposure_percap_park_ssib4_1990s_fpm_cli"   ~ gmt_1990s,
    period_scenario == "exposure_percap_park_ssib4_2000s_fpm_cli"   ~ gmt_2000s,
    period_scenario == "exposure_percap_park_ssib4_2010s_fpm_cli"   ~ gmt_2010s
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
############ Run country-level linear regression ########################################
##
## Model: PM_bar_ctsm = alpha_cm + beta_c * T_ts + epsilon_ctsm
##
## For each country c, we regress population-weighted per-capita fire PM2.5 exposure
## on T_ts with fire-model fixed effects across all 25 observations.
##
## beta_c (slope)        = change in per-capita fire PM2.5 (µg/m^3/yr) per 1°C GMT increase.
## alpha_cm (intercept)  = fire-model-specific intercept for country c. T_ts = 0 lies
##   outside the data range (all obs at ~1°C or higher), so alpha_cm is extrapolated
##   and should not be interpreted as a meaningful exposure estimate.
##
## Countries with fewer than 6 valid observations are dropped (cannot fit 6 parameters).
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ USA-only regression (diagnostic / inspection) ############################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Filter the long-format data to the United States only.
# ISO3 code for the USA is "USA".
usa_data <- reg_data_combined %>%
  filter(country_code_iso3 == "USA")

# Print the USA data to inspect all 25 observations (18 Park decade × fire model + 5 Pierce
# period × scenario + 2 Zhao) and their GMT values that will serve as inputs to the regression.
cat("\n--- USA regression input data ---\n")
print(usa_data %>% select(fire_model, period_scenario, T_ps, exposure_percap) %>% arrange(fire_model))

# Fit the OLS regression for the USA:
#   exposure_percap ~ T_ps + fire_model
#
# lm() fits an ordinary least squares (OLS) linear model.
#   - The first argument is a formula: response ~ predictor(s).
#     Here, exposure_percap is the dependent variable (y) and T_ps is the independent
#     variable (x). lm() automatically includes an intercept unless suppressed with -1.
#   - The second argument, data =, specifies the data frame containing those variables.
#   - lm() returns an object of class "lm" containing fitted coefficients, residuals,
#     the model matrix, and metadata needed for summary(), predict(), and other methods.
#   - The model is: PM_bar_ctsm = alpha_cm + beta_c * T_ts + epsilon_ctsm,
#     where alpha_cm is the fire-model-specific intercept for country c (CLASSIC
#     is the reference level; jules, ssib4, CESM, Zhao enter as dummies relative to it),
#     and beta_c is the GMT slope (key damage function parameter).
#   - OLS minimises the sum of squared residuals to find all 6 parameters.
#   - See last section for how R actually implements this specification. 
usa_model <- lm(exposure_percap ~ T_ps + fire_model, data = usa_data)

# summary() on an lm object prints:
#   - Residuals: min, Q1, median, Q3, max of the fitted residuals
#   - Coefficients table: estimate, standard error, t-statistic, p-value for each term
#   - Residual standard error (RSE): average magnitude of residuals on the response scale
#   - R-squared and adjusted R-squared: fraction of variance in y explained by the model
#   - F-statistic and its p-value: tests whether the model explains significant variance
#     relative to an intercept-only baseline
cat("\n--- USA lm() summary ---\n")
print(summary(usa_model))

usa_beta_c <- coef(usa_model)[["T_ps"]]  # beta_c: slope (µg/m^3/yr per 1°C GMT)

cat("USA beta^(c) (slope, fPM change per 1°C GMT):", round(usa_beta_c, 6), "µg/m^3/yr per °C\n")

# Extract R-squared to assess how well GMT explains USA per-capita fPM variation
# across all 25 observations (18 Park + 5 Pierce + 2 Zhao).
usa_r2 <- summary(usa_model)$r.squared
cat("USA R-squared:", round(usa_r2, 4), "\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Global regression (diagnostic / inspection) ##############################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

global_data <- reg_data_combined %>%
  filter(country_name == "global")

cat("\n--- Global regression input data ---\n")
print(global_data %>% select(fire_model, period_scenario, T_ps, exposure_percap) %>% arrange(fire_model))

global_model <- lm(exposure_percap ~ T_ps + fire_model, data = global_data)

cat("\n--- Global lm() summary ---\n")
print(summary(global_model))

global_beta_c <- coef(global_model)[["T_ps"]]  # beta_c: slope (µg/m^3/yr per 1°C GMT)

cat("Global beta^(c) (slope, fPM change per 1°C GMT):", round(global_beta_c, 6), "µg/m^3/yr per °C\n")

global_r2 <- summary(global_model)$r.squared
cat("Global R-squared:", round(global_r2, 4), "\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Group by country and run one lm() per group using purrr::map inside nest().
# The result is a list-column of tidy regression coefficient tables.
reg_results <- reg_data_combined %>%
  group_by(country_code_iso3, country_name) %>%
  # Keep only countries with at least 6 non-NA observations (minimum to fit 6 parameters;
  # in practice all countries have 25 observations after the pop_bar_c == 0 filter)
  filter(n() >= 6) %>%
  # Nest all observations for each country into a sub-dataframe
  nest() %>%
  mutate(
    # Fit OLS: PM_bar_ctsm = alpha_cm + beta_c * T_ts  (CLASSIC reference + 4 model FE dummies)
    model = purrr::map(data, ~ lm(exposure_percap ~ T_ps + fire_model, data = .x)),

    # Extract tidy coefficient table (term, estimate, std.error, statistic, p.value)
    tidied = purrr::map(model, tidy),

    # Extract model-level fit statistics (r.squared, adj.r.squared, p.value, etc.)
    glanced = purrr::map(model, glance)
  )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Extract beta_c ###########################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Unnest the tidy coefficient table and pivot so beta_c is a column.
# Only the T_ps slope term is retained — the CLASSIC intercept and FE dummies
# are not stored in the output.
reg_coefs <- reg_results %>%
  select(country_code_iso3, country_name, tidied) %>%
  unnest(tidied) %>%
  # Keep only the slope term. FE dummy terms and the intercept are dropped —
  # they absorb model-level mean differences but are not the output of interest.
  filter(term == "T_ps") %>%
  mutate(param = "beta_c") %>%
  select(country_code_iso3, country_name, param, estimate, std.error, statistic, p.value) %>%
  pivot_wider(
    id_cols     = c(country_code_iso3, country_name),
    names_from  = param,
    values_from = c(estimate, std.error, statistic, p.value)
  )

# Inspect coefficient table
print(head(reg_coefs, 10))
colnames(reg_coefs)
cat("\nDimensions of regression coefficient table:", nrow(reg_coefs), "countries\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Summarise beta_c distribution ############################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Print summary statistics for beta_c across all countries.
# beta_c > 0 means higher GMT --> more fire PM2.5 exposure (expected for most countries).
cat("\n--- Summary of beta^(c) (slope: per-capita fPM2.5 change per 1°C GMT) ---\n")
summary(reg_coefs$estimate_beta_c)

# Count countries with positive vs. negative beta^(c)
cat("\nCountries with positive beta^(c) (more fire PM with warming):",
    sum(reg_coefs$estimate_beta_c > 0, na.rm = TRUE), "\n")
cat("Countries with negative beta^(c) (less fire PM with warming):",
    sum(reg_coefs$estimate_beta_c < 0, na.rm = TRUE), "\n")
cat("Share with positive beta^(c):",
    round(mean(reg_coefs$estimate_beta_c > 0, na.rm = TRUE) * 100, 1), "%\n")

print(reg_coefs %>% filter(estimate_beta_c < 0) %>% select(country_code_iso3, estimate_beta_c), n = 50)

# results significant at 5% level
cat("Share with p-value < 0.05:",
    round(mean(reg_coefs$p.value_beta_c < 0.05, na.rm = TRUE) * 100, 1), "%\n")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Confidence intervals for beta_c ##########################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Compute the t critical value for a 95% confidence interval.
#
# Degrees of freedom (DF): each country regression has n = 25 observations and
#   estimates 6 parameters (intercept + slope + 4 model FE dummies), so DF = n - 6 = 19.
df <- nrow(usa_data) - 6
t_critical <- qt(0.975, df)   # ~2.093 for 95% CI with 19 degrees of freedom
cat("\nt_critical (95% CI, df = 19):", round(t_critical, 4), "\n")

t_critical_largesample <- 1.96

reg_coefs <- reg_coefs %>%
  mutate(
    # 95% CI bounds for beta_c using exact t_critical for df = 19
    lower_beta_c        = estimate_beta_c - t_critical * std.error_beta_c,
    upper_beta_c        = estimate_beta_c + t_critical * std.error_beta_c,

    # 95% CI bounds for beta_c using large-sample approximation (1.96)
    lower_beta_c_1.96   = estimate_beta_c - t_critical_largesample * std.error_beta_c,
    upper_beta_c_1.96   = estimate_beta_c + t_critical_largesample * std.error_beta_c,

    # Store both critical values as columns for reference
    t_critical_df19        = t_critical,
    t_critical_largesample = t_critical_largesample,

    # gamma: beta_c scaled by 0.008 RR (Pope et al.) 
    gamma_beta_c        = 0.008 * estimate_beta_c
  ) %>%
  select(
    country_code_iso3, country_name, gamma_beta_c,
    t_critical_df19, t_critical_largesample,
    lower_beta_c, lower_beta_c_1.96,
    estimate_beta_c,
    upper_beta_c_1.96, upper_beta_c,
    std.error_beta_c, statistic_beta_c, p.value_beta_c
  )

cat("\nSample of beta^(c) estimates with confidence bounds:\n")
print(head(reg_coefs %>% select(country_code_iso3, estimate_beta_c,
                                 lower_beta_c, upper_beta_c,
                                 lower_beta_c_1.96, upper_beta_c_1.96), 10))

# note, only country_code_iso3, estimate_beta_c, std.error_beta_c, p.value_beta_c used downstream

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Save output ##############################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Save the coefficient table (beta_c with SEs and p-values) to CSV.
# This is the primary output used downstream in the GIVE damage function.
write_csv(reg_coefs, here("output", "fpm_gmt_regression_coefs_FE_all_cli.csv"))

cat("\nSaved regression coefficients to output/fpm_gmt_regression_coefs_FE_all_cli.csv\n")

# Save combined long-format regression input data for use in descriptive statistics script.
write_csv(reg_data_combined, here("output", "reg_data_combined_fe_all_cli.csv"))
cat("Saved combined regression input data to output/reg_data_combined_fe_all_cli.csv\n")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ R implementation of regression ##############################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# R implementation of PM_bar_ctsm = alpha_cm + beta_c * T_ts + epsilon_ctsm:
#
# Because fire_model is a factor with classic as the reference level, R automatically
# creates four indicator (dummy) variables — one for each non-reference model:
#   D_jules = 1 if fire_model == "jules", 0 otherwise
#   D_ssib4 = 1 if fire_model == "ssib4", 0 otherwise
#   D_CESM  = 1 if fire_model == "CESM",  0 otherwise
#   D_Zhao  = 1 if fire_model == "Zhao",  0 otherwise
#
# The expanded model estimated by lm() is:
#   PM_bar = alpha_c_classic
#            + delta_jules  * D_jules
#            + delta_ssib4  * D_ssib4
#            + delta_CESM   * D_CESM
#            + delta_Zhao   * D_Zhao
#            + beta_c       * T_ts
#            + epsilon
#
# alpha_cm for each model is therefore:
#   classic --> (Intercept)
#   jules   --> (Intercept) + fire_modeljules
#   ssib4   --> (Intercept) + fire_modelssib4
#   CESM    --> (Intercept) + fire_modelCESM
#   Zhao    --> (Intercept) + fire_modelZhao
#
# beta_c is the T_ts coefficient — shared across all models for country c.


# THE END 
