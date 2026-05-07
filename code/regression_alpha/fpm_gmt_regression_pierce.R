# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## FPM–GMT RELATIONSHIP: Estimate alpha_c2 (per-capita fire PM2.5 change per 1°C GMT)
##                        Pierce et al. projections only
##
## Goal: For each country c, estimate the linear regression:
##   exposure_cps^per-capita = alpha_c1 + alpha_c2 * T_ps
##
## where T_ps is GMT change in period p under scenario s (°C relative to 1850–1900),
## and alpha_c2 is the key damage function parameter: the change in per-capita fire
## PM2.5 exposure (µg/m³/person/year) per 1°C increase in GMT.
##
## Data: Pierce et al. projections only:
##   - baseline (~2001–2010), 2041–2050, and 2091–2100 under RCP4.5 and RCP8.5
##   - 5 (period × scenario) observations per country.
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

# Import country-level per-capita fire PM2.5 exposure (µg/m³/person/year).
# Each row is a country; columns include per-capita exposure for each period × scenario:
#   exposure_percap_fpm_2000         --> baseline period (~2001–2010)
#   exposure_percap_fpm_2050_45/85   --> 2041–2050 under RCP4.5 / RCP8.5
#   exposure_percap_fpm_2100_45/85   --> 2091–2100 under RCP4.5 / RCP8.5

#pop_wght <- read_csv(here("output", "pop_wght_pm_cntry.csv"))
pop_wght <- read_csv(here("output", "pop_wght_pm_cntry_final.csv"))

# Drop GIVE countries with no corresponding data in the Pierce PM gridded data.
# These rows have pop_bar_c == 0, indicating no population-weighted exposure could be computed.
pop_wght <- pop_wght %>%
  filter(pop_bar_c != 0)

# Import decadal mean GMT anomaly (°C relative to 1850–1900 pre-industrial baseline)
# for each period × scenario combination.
# Rows: "2006-2010" (baseline), "2041-2050", "2091-2100"
# Columns: mean_gmt_45 (RCP4.5), mean_gmt_85 (RCP8.5)
gmt_chg <- read_csv(here("output", "gmt_periods_pi.csv"))


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

# Print GMT values
cat("GMT baseline (2006-2010 avg):", gmt_baseline, "\n")
cat("GMT 2040s RCP4.5:", gmt_2040s_45, "  RCP8.5:", gmt_2040s_85, "\n")
cat("GMT 2090s RCP4.5:", gmt_2090s_45, "  RCP8.5:", gmt_2090s_85, "\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Reshape to long format for regression ######################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# We need one row per (country, period, scenario) with columns:
#   country_code_iso3, country_name, T_ps (GMT), exposure_percap (µg/m³/person/yr)
#
# This long format lets us run a single lm() call per country across all 5 data points.

# Step 1: Select only the identifier columns and the five Pierce per-capita exposure columns.
#         Drop rows where all exposure values are NA (e.g., uninhabited territories).
reg_data_wide <- pop_wght %>%
  select(country_code_iso3,
         country_name,
         exposure_percap_fpm_2000,    # baseline period exposure
         exposure_percap_fpm_2050_45, # 2040s RCP4.5 exposure
         exposure_percap_fpm_2050_85, # 2040s RCP8.5 exposure
         exposure_percap_fpm_2100_45, # 2090s RCP4.5 exposure
         exposure_percap_fpm_2100_85) # 2090s RCP8.5 exposure

# Step 2: Pivot to long format so each row is one (country, period×scenario) observation.
#         The column name encodes which GMT value applies to that row.
reg_data_long <- reg_data_wide %>%
  pivot_longer(
    cols      = starts_with("exposure_percap_fpm_"), # the five Pierce exposure columns
    names_to  = "period_scenario",                   # new column holding the column name
    values_to = "exposure_percap"                    # new column holding the exposure value
  ) %>%
  # Step 3: Map each period×scenario label to its corresponding GMT value (T_ps).
  #         This is the regressor in the country-level linear regression.
  mutate(T_ps = case_when(
    period_scenario == "exposure_percap_fpm_2000"    ~ gmt_baseline, # baseline: avg of 45/85
    period_scenario == "exposure_percap_fpm_2050_45" ~ gmt_2040s_45, # 2040s RCP4.5
    period_scenario == "exposure_percap_fpm_2050_85" ~ gmt_2040s_85, # 2040s RCP8.5
    period_scenario == "exposure_percap_fpm_2100_45" ~ gmt_2090s_45, # 2090s RCP4.5
    period_scenario == "exposure_percap_fpm_2100_85" ~ gmt_2090s_85  # 2090s RCP8.5
  )) %>%
  # Drop rows with missing exposure (uninhabited territories or missing data).
  # Countries with NA exposure cannot contribute to the regression.
  filter(!is.na(exposure_percap))

# Inspect the reshaped data to confirm structure (5 rows per country)
print(head(reg_data_long, 15))

# Sanity check: every country should have exactly 5 rows (Pierce only).
cat("\nTotal rows:", nrow(reg_data_long), "\n")
print(reg_data_long %>% count(country_name) %>% summary())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Run country-level linear regression ########################################
##
## Model: exposure_cps^per-capita = alpha_c1 + alpha_c2 * T_ps
##
## For each country c, we regress per-capita fire PM2.5 exposure on GMT (T_ps)
## across the 5 (period × scenario) data points.
##
## alpha_c2 (slope) = change in per-capita fire PM2.5 (µg/m³/yr) per 1°C GMT increase.
## alpha_c1 (intercept) = OLS-fitted y-intercept; a mathematical artefact of the line fit,
##   not a meaningful estimate of exposure at pre-industrial temperatures. All observed
##   T_ps values are ~1°C or higher, so T=0 lies outside the data range.
##
## Countries with fewer than 2 valid observations are dropped (cannot fit a line).
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ USA-only regression (diagnostic / inspection) ############################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Filter the long-format data to the United States only.
# ISO3 code for the USA is "USA".
usa_data <- reg_data_long %>%
  filter(country_code_iso3 == "USA")

# Print the USA data to inspect all 5 Pierce observations and their GMT values.
cat("\n--- USA regression input data ---\n")
print(usa_data %>% select(period_scenario, T_ps, exposure_percap))

# Fit the OLS regression for the USA:
#   exposure_percap ~ T_ps
#
# lm() fits an ordinary least squares (OLS) linear model.
#   - The first argument is a formula: response ~ predictor(s).
#     Here, exposure_percap is the dependent variable (y) and T_ps is the independent
#     variable (x). lm() automatically includes an intercept unless suppressed with -1.
#   - The second argument, data =, specifies the data frame containing those variables.
#   - lm() returns an object of class "lm" containing fitted coefficients, residuals,
#     the model matrix, and metadata needed for summary(), predict(), and other methods.
#   - With one predictor, the model is: y = beta_0 + beta_1 * x + epsilon,
#     where beta_0 is the intercept (alpha_c1) and beta_1 is the slope (alpha_c2).
#   - OLS minimises the sum of squared residuals to find beta_0 and beta_1.
usa_model <- lm(exposure_percap ~ T_ps, data = usa_data)

# summary() on an lm object prints:
#   - Residuals: min, Q1, median, Q3, max of the fitted residuals
#   - Coefficients table: estimate, standard error, t-statistic, p-value for each term
#   - Residual standard error (RSE): average magnitude of residuals on the response scale
#   - R-squared and adjusted R-squared: fraction of variance in y explained by the model
#   - F-statistic and its p-value: tests whether the model explains significant variance
#     relative to an intercept-only baseline
cat("\n--- USA lm() summary ---\n")
print(summary(usa_model))

# Extract and label the two key coefficients explicitly for readability.
# Note: alpha_c1 is the OLS y-intercept — the value the fitted line takes at T_ps = 0.
# All observed T_ps values are ~1°C or above, so this is an extrapolation outside the
# data range and should not be interpreted as a meaningful exposure estimate.
usa_alpha_c1 <- coef(usa_model)[["(Intercept)"]]  # alpha_c1: OLS y-intercept (extrapolated)
usa_alpha_c2 <- coef(usa_model)[["T_ps"]]          # alpha_c2: slope (µg/m³/yr per 1°C GMT)

cat("\nUSA alpha_c1 (OLS intercept, extrapolated outside data range):", round(usa_alpha_c1, 6), "µg/m³/yr\n")
cat("USA alpha_c2 (slope, fPM change per 1°C GMT):", round(usa_alpha_c2, 6), "µg/m³/yr per °C\n")

# Extract R-squared to assess how well GMT explains USA per-capita fPM variation
# across the 5 Pierce observations.
usa_r2 <- summary(usa_model)$r.squared
cat("USA R-squared:", round(usa_r2, 4), "\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Global regression (diagnostic / inspection) ##############################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

global_data <- reg_data_long %>%
  filter(country_name == "global")

cat("\n--- Global regression input data ---\n")
print(global_data %>% select(period_scenario, T_ps, exposure_percap))

global_model <- lm(exposure_percap ~ T_ps, data = global_data)

cat("\n--- Global lm() summary ---\n")
print(summary(global_model))

global_alpha_c1 <- coef(global_model)[["(Intercept)"]]
global_alpha_c2 <- coef(global_model)[["T_ps"]]

cat("\nGlobal alpha_c1 (OLS intercept, extrapolated outside data range):", round(global_alpha_c1, 6), "µg/m³/yr\n")
cat("Global alpha_c2 (slope, fPM change per 1°C GMT):", round(global_alpha_c2, 6), "µg/m³/yr per °C\n")

global_r2 <- summary(global_model)$r.squared
cat("Global R-squared:", round(global_r2, 4), "\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Group by country and run one lm() per group using purrr::map inside nest().
# The result is a list-column of tidy regression coefficient tables.
reg_results <- reg_data_long %>%
  group_by(country_code_iso3, country_name) %>%
  # Keep only countries with at least 2 non-NA observations (minimum to fit a line;
  # in practice all countries have 5 observations after the pop_bar_c == 0 filter)
  filter(n() >= 2) %>%
  # Nest all observations for each country into a sub-dataframe
  nest() %>%
  mutate(
    # Fit OLS: exposure_percap ~ T_ps  (intercept + slope on GMT)
    model = purrr::map(data, ~ lm(exposure_percap ~ T_ps, data = .x)),

    # Extract tidy coefficient table (term, estimate, std.error, statistic, p.value)
    tidied = purrr::map(model, tidy),

    # Extract model-level fit statistics (r.squared, adj.r.squared, p.value, etc.)
    glanced = purrr::map(model, glance)
  )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Extract alpha_c1 and alpha_c2 ############################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Unnest the tidy coefficient table and pivot so alpha_c1 and alpha_c2 are columns.
# 'term' will be "(Intercept)" for alpha_c1 and "T_ps" for alpha_c2.
reg_coefs <- reg_results %>%
  select(country_code_iso3, country_name, tidied) %>%
  unnest(tidied) %>%
  # Rename the regression terms to the paper's notation for clarity
  mutate(param = case_when(
    term == "(Intercept)" ~ "alpha_c1",  # OLS y-intercept (extrapolated; T=0 outside data range)
    term == "T_ps"        ~ "alpha_c2"   # slope: exposure change per 1°C GMT (the key parameter)
  )) %>%
  select(country_code_iso3, country_name, param, estimate, std.error, statistic, p.value) %>%
  # Pivot wide so each country has one row with columns alpha_c1 and alpha_c2
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
############ Summarise alpha_c2 distribution ##########################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Print summary statistics for alpha_c2 across all countries.
# alpha_c2 > 0 means higher GMT → more fire PM2.5 exposure (expected for most countries).
cat("\n--- Summary of alpha_c2 (slope: per-capita fPM2.5 change per 1°C GMT) ---\n")
summary(reg_coefs$estimate_alpha_c2)

# Count countries with positive vs. negative alpha_c2
cat("\nCountries with positive alpha_c2 (more fire PM with warming):",
    sum(reg_coefs$estimate_alpha_c2 > 0, na.rm = TRUE), "\n")
cat("Countries with negative alpha_c2 (less fire PM with warming):",
    sum(reg_coefs$estimate_alpha_c2 < 0, na.rm = TRUE), "\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Confidence intervals for alpha_c2 ########################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Compute the t critical value for a 95% confidence interval.
#
# Degrees of freedom (DF): each country regression has n = 5 observations (Pierce only)
#   and estimates 2 parameters (intercept + slope), so DF = n - 2 = 3.
#
# With 3 DF, t_critical is ~3.18 — notably wider than the large-sample 1.96,
#   reflecting the limited number of observations.
df <- nrow(usa_data)-2
t_critical <- qt(0.975, df)   # ~3.182 for 95% CI with 3 degrees of freedom
cat("\nt_critical (95% CI, df = 3):", round(t_critical, 4), "\n")

t_critical_largesample <- 1.96

reg_coefs <- reg_coefs %>%
  mutate(
    # 95% CI bounds using exact t_critical for df = 3
    lower_alpha_c2        = estimate_alpha_c2 - t_critical * std.error_alpha_c2,
    upper_alpha_c2        = estimate_alpha_c2 + t_critical * std.error_alpha_c2,

    # 95% CI bounds using hardcoded large-sample approximation (1.96)
    lower_alpha_c2_1.96   = estimate_alpha_c2 - t_critical_largesample * std.error_alpha_c2,
    upper_alpha_c2_1.96   = estimate_alpha_c2 + t_critical_largesample * std.error_alpha_c2,

    # Store both critical values as columns for reference
    t_critical_df3         = t_critical,
    t_critical_largesample = t_critical_largesample,

    # gamma_alpha_c2: scaled slope used in the GIVE damage function. 0.0095 from Orellano et al 2024 RR
    gamma_alpha_c2        = 0.0095 * estimate_alpha_c2
  ) %>%
  select(
    country_code_iso3, country_name, gamma_alpha_c2,
    t_critical_df3, t_critical_largesample,
    lower_alpha_c2, lower_alpha_c2_1.96,
    estimate_alpha_c2,
    upper_alpha_c2_1.96, upper_alpha_c2,
    std.error_alpha_c2, statistic_alpha_c2, p.value_alpha_c2,
    estimate_alpha_c1, std.error_alpha_c1, statistic_alpha_c1, p.value_alpha_c1
  )

cat("\nSample of alpha_c2 estimates with confidence bounds:\n")
print(head(reg_coefs %>% select(country_code_iso3, estimate_alpha_c2,
                                 lower_alpha_c2, upper_alpha_c2,
                                 lower_alpha_c2_1.96, upper_alpha_c2_1.96), 10))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Save output ##############################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Save the full coefficient table (alpha_c1, alpha_c2 with SEs and p-values) to CSV.
write_csv(reg_coefs, here("output", "fpm_gmt_regression_coefs_pierce.csv"))

cat("\nSaved regression coefficients to output/fpm_gmt_regression_coefs_pierce.csv\n")


# THE END 
