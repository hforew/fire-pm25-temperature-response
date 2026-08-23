# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## FPM–GMT RELATIONSHIP: Main specification and leave-one-model-out sensitivity for beta_c
##                        (per-capita fire PM2.5 change per 1degC GMT)
##
## Goal: For each country c, estimate the linear regression with fire-model fixed effects:
##   PM_bar_ctsm = alpha_cm + beta_c * T_ts + epsilon_ctsm
##
## where PM_bar_ctsm is population-weighted per-capita fire PM2.5 exposure (µg/m^3/yr)
## for country c, time period t, scenario s, and fire model m;
## T_ts is GMT anomaly relative to the 1850-1900 pre-industrial baseline (degC);
## beta_c is the key damage function parameter (change in population-weighted
## per-capita fire PM2.5 per 1degC GMT), and alpha_cm are fire-model-specific
## intercepts (fixed effects absorbing model-level mean differences in exposure levels).
##
## This script refits the regression once per "exclusion group" -- the full 5-FE-level
## model, then 5 leave-one-model-out variants that each drop one fire-model source's
## observations entirely -- to test whether beta_c depends on any single data source.
## Park factual/counterfactual pairs of the same model
## (e.g. ssib4 factual + ssib4 counterfactual) are always dropped together, since they
## already share one FE group; CESM and Zhao are dropped alone.
##
## Five fire-model fixed effect levels (full model):
##   classic (reference), jules, ssib4 --> Park et al. factual (obsclim) AND counterfactual
##                                         (counterclim) runs of the same model pooled into
##                                         one intercept -- distinguished by T_ps (0 for
##                                         counterfactual), not by a separate FE level
##   CESM                               --> Pierce et al. projections
##   Zhao                               --> Zhao et al. projections
##
## Data spans three sources combined:
##   - Pierce et al. projections: baseline (~2001–2010), 2041–2050, and 2091–2100
##     under RCP4.5 and RCP8.5 — 5 (period × scenario) observations per country.
##   - Park et al. historical factual: 1960s–2010s across 3 fire models (classic, JULES, SSIB4)
##     — 18 (decade × fire model) observations per country.
##   - Park et al. historical counterfactual: 1960s–2010s across the same 3 fire models
##     — 18 observations per country, pooled into the same FE group as their factual runs.
##   - Zhao et al. projections: ~2095 under SSP2-4.5 and SSP5-8.5 — 2 observations per country.
##   Combined: up to 43 observations per country for the full model; fewer for each
##   leave-one-model-out exclusion group (see exclusion_groups below).
##
## Inputs:
##   output/pop_wght_pm_cntry_final.csv   (country-level per-capita fPM exposure; Pierce,
##                                         Park, and Zhao columns)
##   output/gmt/gmt_pierce_RCPs.csv       (decadal mean GMT anomaly for Pierce periods ×
##                                         scenarios, relative to 1850–1900 PI baseline)
##   output/gmt/gmt_park_hist.csv         (GMT anomaly for each Park snapshot decade,
##                                         relative to 1850–1900 PI baseline, sourced from OWID)
##   output/gmt/gmt_zhao_SSPs.csv         (GMT anomaly for Zhao ~2095 SSP245 and SSP585,
##                                         sourced from MimiSSPs)
##
## Outputs (one pair per exclusion group -- full, excl_classic, excl_jules, excl_ssib4,
## excl_CESM, excl_Zhao -- written to output/betas/):
##   fpm_gmt_betas_country_<group>.csv    (beta_c per country with SE, CI,
##                                        p-value, and gamma_beta_c for
##                                        that exclusion group)
##   regress_data_country_<group>.csv     (that group's combined long-format
##                                        regression input data)
##
## Execution order:
##   files run before: pop_wght_pm_cntry.R --> writes pop_wght_pm_cntry_final.csv
##   files run after: regress_desc_stats.R, beta_comparison_latex.R,
##                    beta_GIVE_country_map.R (all read this script's "full" group output)
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

pop_wght <- read_csv(here("output", "pop_wght_pm_cntry_final.csv"))

# Drop GIVE countries with no corresponding data in the Pierce or Park PM gridded data.
# These rows have pop_bar_c == 0, indicating no population-weighted exposure could be computed.
pop_wght <- pop_wght %>%
  filter(pop_bar_c != 0)

# Import decadal mean GMT anomaly (degC relative to 1850–1900 pre-industrial baseline)
# for each period × scenario combination.
# Rows: "2006-2010" (baseline), "2041-2050", "2091-2100"
# Columns: mean_gmt_45 (RCP4.5), mean_gmt_85 (RCP8.5)
gmt_chg <- read_csv(here("output", "gmt", "gmt_pierce_RCPs.csv"))

# Import Park decade GMT values (degC relative to pre-industrial baseline)
# Rows: one per Park snapshot decade (1960s–2010s)
# Columns: park_year, decade, mean_gmt_pi
gmt_park <- read_csv(here("output", "gmt", "gmt_park_hist.csv"))

# Import Zhao et al. ~2095 GMT anomalies (degC relative to pre-industrial baseline)
# Rows: SSP245, SSP585 --- these temp changes assigned to Zhao but sourced from MimiSSPs
# Columns: scenario, mean_gmt_pi
gmt_zhao <- read_csv(here("output", "gmt", "gmt_zhao_SSPs.csv"))

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

# Zhao et al. ~2095 GMT scalars from gmt_zhao_SSPs.csv
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

# Inspect the reshaped data to confirm structure (7 rows per country, fire_model visible)
nrow(reg_data_long[reg_data_long$country_code_iso3 == "AFG", ])
print(reg_data_long[reg_data_long$country_code_iso3 == "AFG", ])

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
############ Reshape Park counterfactual data to long format ############################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Mirrors the factual Park reshape above, using counterfactual (counterclim) exposure.
# One row per (country, fire_model, decade) — 18 obs per country. fire_model is labelled
# to match the factual reshape (e.g. "classic", not "classic_counter") so factual and
# counterfactual runs of the same Park model share one FE group / intercept below —
# they're distinguished by T_ps (0 for counterfactual) rather than by fire_model.

park_counter_data_wide <- pop_wght %>%
  select(country_code_iso3,
         country_name,
         exposure_percap_park_classic_1960s_fpm_counter,   # classic counterclim, 1960s
         exposure_percap_park_classic_1970s_fpm_counter,   # classic counterclim, 1970s
         exposure_percap_park_classic_1980s_fpm_counter,   # classic counterclim, 1980s
         exposure_percap_park_classic_1990s_fpm_counter,   # classic counterclim, 1990s
         exposure_percap_park_classic_2000s_fpm_counter,   # classic counterclim, 2000s
         exposure_percap_park_classic_2010s_fpm_counter,   # classic counterclim, 2010s
         exposure_percap_park_jules_1960s_fpm_counter,     # JULES counterclim, 1960s
         exposure_percap_park_jules_1970s_fpm_counter,     # JULES counterclim, 1970s
         exposure_percap_park_jules_1980s_fpm_counter,     # JULES counterclim, 1980s
         exposure_percap_park_jules_1990s_fpm_counter,     # JULES counterclim, 1990s
         exposure_percap_park_jules_2000s_fpm_counter,     # JULES counterclim, 2000s
         exposure_percap_park_jules_2010s_fpm_counter,     # JULES counterclim, 2010s
         exposure_percap_park_ssib4_1960s_fpm_counter,     # SSIB4 counterclim, 1960s
         exposure_percap_park_ssib4_1970s_fpm_counter,     # SSIB4 counterclim, 1970s
         exposure_percap_park_ssib4_1980s_fpm_counter,     # SSIB4 counterclim, 1980s
         exposure_percap_park_ssib4_1990s_fpm_counter,     # SSIB4 counterclim, 1990s
         exposure_percap_park_ssib4_2000s_fpm_counter,     # SSIB4 counterclim, 2000s
         exposure_percap_park_ssib4_2010s_fpm_counter)     # SSIB4 counterclim, 2010s

park_counter_data_long <- park_counter_data_wide %>%
  pivot_longer(
    cols      = starts_with("exposure_percap_park_"),
    names_to  = "period_scenario",
    values_to = "exposure_percap"
  ) %>%
  # Counterclim runs hold CO2 (and thus GMT) fixed at 1901 levels, so all
  # counterfactual decades share T_ps = 0 (degC anomaly rel. to 1850-1900 PI baseline).
  mutate(T_ps = 0) %>%
  mutate(fire_model = case_when(
    grepl("classic", period_scenario) ~ "classic",
    grepl("jules",   period_scenario) ~ "jules",
    grepl("ssib4",   period_scenario) ~ "ssib4"
  )) %>%
  filter(!is.na(exposure_percap))

# Inspect the reshaped counterfactual Park data (18 rows per country — 3 models × 6 decades)
print(head(park_counter_data_long, 18))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Combine historical (Park) and future projection data ########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# One regression per country across all observations:
#   18 Park factual obs (classic, jules, ssib4 × 6 decades)
#   + 18 Park counterfactual obs (classic, jules, ssib4 × 6 decades — same fire_model
#     labels as the factual obs, so each Park model's factual and counterfactual runs
#     pool into one FE group / intercept)
#   + 5 Pierce obs (1 base + 4 future)
#   + 2 Zhao obs (SSP2-4.5, SSP5-8.5) = 43 total.
reg_data_combined <- bind_rows(park_data_long, park_counter_data_long, reg_data_long) %>%
  # Set fire_model as factor with "classic" as reference level. FE dummies for
  # jules, ssib4, CESM, Zhao are all interpreted relative to classic. classic, jules,
  # and ssib4 each pool their Park factual + counterfactual observations into one
  # intercept (5 FE groups total, not 8).
  mutate(fire_model = factor(fire_model,
                             levels = c("classic", "jules", "ssib4",
                                        "CESM", "Zhao")))

# Sanity check: total rows and per-country observation counts.
# Every country should have exactly 43 rows (18 factual + 18 counter + 5 Pierce + 2 Zhao).
cat("\nCombined rows:", nrow(reg_data_combined), "\n")
print(reg_data_combined %>% count(country_name) %>% summary())
# Expected: Combined rows: 7525 | n Min=43, Mean=43, Max=43 (175 countries * 43 obs)

# Transparency: show observation count by fire_model across all countries.
# Per country: classic=12 (6 factual + 6 counterfactual), jules=12, ssib4=12,
#              CESM=5, Zhao=2 (total=43).
cat("\nTotal observations by fire_model (all countries combined):\n")
print(reg_data_combined %>% count(fire_model))
# Expected: classic=2100, jules=2100, ssib4=2100, CESM=875, Zhao=350 (175 countries each)

# Inspect the reshaped data to confirm structure ( fire_model visible)
afg <- reg_data_combined %>%
  filter(country_code_iso3 == "AFG")

print(head(afg, 23), n = 23)
print(tail(afg, 23), n = 23)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Leave-one-model-out exclusion groups #########################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Sensitivity test: does beta_c depend on any single fire-model data source? "full" is
## the baseline (all 5 FE levels). Each other entry
## drops one fire-model source entirely before refitting. Since Park factual and
## counterfactual runs of the same model already share one FE group (see fire_model
## above), excluding e.g. "classic" drops both its factual and counterfactual
## observations together; CESM (Pierce) and Zhao are each dropped alone.
exclusion_groups <- list(
  full         = character(0),
  excl_classic = c("classic"),
  excl_jules   = c("jules"),
  excl_ssib4   = c("ssib4"),
  excl_CESM    = c("CESM"),
  excl_Zhao    = c("Zhao")
)

# Create the sensitivity output directory if it doesn't already exist.
if (!dir.exists(here("output", "betas"))) {
  dir.create(here("output", "betas"), recursive = TRUE)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Run country-level linear regression, once per exclusion group ################
##
## Model: PM_bar_ctsm = alpha_cm + beta_c * T_ts + epsilon_ctsm
##
## For each country c, we regress population-weighted per-capita fire PM2.5 exposure
## on T_ts with fire-model fixed effects. The full model uses all 43 observations across
## 5 FE levels (1 intercept + 4 dummies + 1 slope = 6 parameters, df = 43 - 6 = 37); each
## exclusion group below drops one source's observations and re-derives its own parameter
## count and degrees of freedom from however many FE levels remain.
##
## beta_c (slope)        = change in per-capita fire PM2.5 (µg/m^3/yr) per 1degC GMT increase.
## alpha_cm (intercept)  = fire-model-specific intercept for country c. T_ts = 0 lies
##   outside the data range (all obs at ~1degC or higher), so alpha_cm is extrapolated
##   and should not be interpreted as a meaningful exposure estimate.
##
## Countries with fewer valid observations than parameters are dropped (cannot fit the model).
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Precompute each exclusion group's filtered data (excluded fire_model level(s) dropped,
# unused factor levels dropped so lm() doesn't try to estimate an all-zero dummy column)
# and parameter count once, up front. Every loop below reuses this instead of
# re-filtering, so the filtering logic lives in exactly one place.
group_data_list <- purrr::map(exclusion_groups, function(excluded_models) {
  reg_data_combined %>%
    filter(!fire_model %in% excluded_models) %>%
    droplevels()
})

# Parameters = (FE levels - 1) dummies + intercept + slope = FE levels + 1.
n_params_list <- purrr::map_dbl(group_data_list, ~ n_distinct(.x$fire_model) + 1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ USA-only regression (diagnostic / inspection) ############################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for (group_name in names(exclusion_groups)) {

  group_data <- group_data_list[[group_name]]

  cat("\n============================================================\n")
  cat("[", group_name, "] USA-only regression (diagnostic / inspection)\n", sep = "")
  cat("============================================================\n")

  # Filter this exclusion group's data to the United States only.
  usa_data <- group_data %>%
    filter(country_code_iso3 == "USA")

  # Print the USA data to inspect this group's observations and the GMT values that
  # will serve as inputs to the regression.
  cat("\n--- [", group_name, "] USA regression input data ---\n", sep = "")
  print(
    usa_data %>%
      select(fire_model, period_scenario, T_ps, exposure_percap) %>%
      arrange(fire_model),
    n = nrow(usa_data)
  )

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
  #     where alpha_cm is the fire-model-specific intercept for country c (CLASSIC is the
  #     reference level, unless dropped by this exclusion group, in which case the next
  #     remaining level becomes the reference; other remaining levels enter as dummies),
  #     and beta_c is the GMT slope (key damage function parameter).
  #   - OLS minimises the sum of squared residuals to find all n_params parameters.
  #   - See last section for how R actually implements this specification.
  usa_model <- lm(exposure_percap ~ T_ps + fire_model, data = usa_data)

  # summary() on an lm object prints:
  #   - Residuals: min, Q1, median, Q3, max of the fitted residuals
  #   - Coefficients table: estimate, standard error, t-statistic, p-value for each term
  #   - Residual standard error (RSE): average magnitude of residuals on the response scale
  #   - R-squared and adjusted R-squared: fraction of variance in y explained by the model
  #   - F-statistic and its p-value: tests whether the model explains significant variance
  #     relative to an intercept-only baseline
  cat("\n--- [", group_name, "] USA lm() summary ---\n", sep = "")
  print(summary(usa_model))

  usa_beta_c <- coef(usa_model)[["T_ps"]]  # beta_c: slope (µg/m^3/yr per 1degC GMT)

  cat("[", group_name, "] USA beta^(c) (slope, fPM change per 1°C GMT): ",
      round(usa_beta_c, 6), " µg/m^3/yr per °C\n", sep = "")

  # Extract R-squared to assess how well GMT explains USA per-capita fPM variation
  # across this exclusion group's observations.
  usa_r2 <- summary(usa_model)$r.squared
  cat("[", group_name, "] USA R-squared: ", round(usa_r2, 4), "\n", sep = "")
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Global regression (diagnostic / inspection) ##############################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for (group_name in names(exclusion_groups)) {

  group_data <- group_data_list[[group_name]]

  cat("\n============================================================\n")
  cat("[", group_name, "] Global regression (diagnostic / inspection)\n", sep = "")
  cat("============================================================\n")

  global_data <- group_data %>%
    filter(country_name == "global")

  cat("\n--- [", group_name, "] Global regression input data ---\n", sep = "")
  print(
    global_data %>%
      select(fire_model, period_scenario, T_ps, exposure_percap) %>%
      arrange(fire_model),
    n = nrow(global_data)
  )

  global_model <- lm(exposure_percap ~ T_ps + fire_model, data = global_data)

  cat("\n--- [", group_name, "] Global lm() summary ---\n", sep = "")
  print(summary(global_model))

  global_beta_c <- coef(global_model)[["T_ps"]]  # beta_c: slope (µg/m^3/yr per 1degC GMT)

  cat("[", group_name, "] Global beta^(c) (slope, fPM change per 1°C GMT): ",
      round(global_beta_c, 6), " µg/m^3/yr per °C\n", sep = "")

  global_r2 <- summary(global_model)$r.squared
  cat("[", group_name, "] Global R-squared: ", round(global_r2, 4), "\n", sep = "")
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Per-country regression, confidence intervals, and save ####################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# reg_coefs and reg_results are reassigned each loop iteration below, so anything
# sourcing this script (e.g. plots_beta_lobf.R) that wants a specific group's results
# after the loop finishes should pull them from these group-keyed lists rather than
# relying on the bare reg_coefs/reg_results objects, which only hold the last group run.
coefs_by_group   <- list()
results_by_group <- list()

for (group_name in names(exclusion_groups)) {

  group_data <- group_data_list[[group_name]]
  n_params   <- n_params_list[[group_name]]

  cat("\n============================================================\n")
  cat("[", group_name, "] Per-country regression, CI, and save\n", sep = "")
  cat("============================================================\n")

  # Needed only for the degrees-of-freedom calculation below (same USA row count used
  # in the USA-only diagnostic loop above).
  usa_data <- group_data %>%
    filter(country_code_iso3 == "USA")

  # Group by country and run one lm() per group using purrr::map inside nest().
  # The result is a list-column of tidy regression coefficient tables.
  reg_results <- group_data %>%
    group_by(country_code_iso3, country_name) %>%
    # Keep only countries with at least n_params non-NA observations (minimum to fit
    # the model; in practice all countries have the same observation count within a
    # given exclusion group after the pop_bar_c == 0 filter)
    filter(n() >= n_params) %>%
    # Nest all observations for each country into a sub-dataframe
    nest() %>%
    mutate(
      # Fit OLS: PM_bar_ctsm = alpha_cm + beta_c * T_ts (remaining FE dummies for this group)
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
  # Only the T_ps slope term is retained — the reference intercept and FE dummies
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

  # Store this group's results under its own name before the next loop iteration
  # overwrites reg_coefs/reg_results.
  coefs_by_group[[group_name]]   <- reg_coefs
  results_by_group[[group_name]] <- reg_results

  # Inspect coefficient table
  cat("\n--- [", group_name, "] Coefficient table (head) ---\n", sep = "")
  print(head(reg_coefs, 10))
  cat("[", group_name, "] Dimensions of regression coefficient table: ",
      nrow(reg_coefs), " countries\n", sep = "")

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ############ Summarise beta_c distribution ############################################
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Print summary statistics for beta_c across all countries.
  # beta_c > 0 means higher GMT --> more fire PM2.5 exposure (expected for most countries).
  cat("\n--- [", group_name, "] Summary of beta^(c) (slope: per-capita fPM2.5 change per 1°C GMT) ---\n", sep = "")
  print(summary(reg_coefs$estimate_beta_c))

  # Count countries with positive vs. negative beta^(c)
  cat("[", group_name, "] Countries with positive beta^(c) (more fire PM with warming): ",
      sum(reg_coefs$estimate_beta_c > 0, na.rm = TRUE), "\n", sep = "")
  cat("[", group_name, "] Countries with negative beta^(c) (less fire PM with warming): ",
      sum(reg_coefs$estimate_beta_c < 0, na.rm = TRUE), "\n", sep = "")

  cat("\n--- [", group_name, "] Countries with negative beta^(c) ---\n", sep = "")
  print(reg_coefs %>% filter(estimate_beta_c < 0) %>% select(country_code_iso3, estimate_beta_c), n = 50)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ############ Confidence intervals for beta_c ##########################################
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Compute the t critical value for a 95% confidence interval.
  #
  # Degrees of freedom (DF) depend on this exclusion group's observation count and
  # parameter count: DF = nrow(usa_data) - n_params (37 for the full model;
  # fewer for groups with fewer remaining FE levels).
  df <- nrow(usa_data) - n_params
  t_critical <- qt(0.975, df)
  cat("\n[", group_name, "] t_critical (95% CI, df = ", df, "): ", round(t_critical, 4), "\n", sep = "")

  t_critical_largesample <- 1.96

  reg_coefs <- reg_coefs %>%
    mutate(
      # 95% CI bounds for beta_c using this group's exact t_critical
      lower_beta_c        = estimate_beta_c - t_critical * std.error_beta_c,
      upper_beta_c        = estimate_beta_c + t_critical * std.error_beta_c,

      # 95% CI bounds for beta_c using large-sample approximation (1.96)
      lower_beta_c_1.96   = estimate_beta_c - t_critical_largesample * std.error_beta_c,
      upper_beta_c_1.96   = estimate_beta_c + t_critical_largesample * std.error_beta_c,

      # Store this group's critical value, df, and large-sample critical value for reference
      t_critical_df          = t_critical,
      df_beta_c               = df,
      t_critical_largesample = t_critical_largesample,

      # gamma: beta_c scaled by 0.008 RR (Pope et al.)
      gamma_beta_c        = 0.008 * estimate_beta_c
    ) %>%
    select(
      country_code_iso3, country_name, gamma_beta_c,
      t_critical_df, df_beta_c, t_critical_largesample,
      lower_beta_c, lower_beta_c_1.96,
      estimate_beta_c,
      upper_beta_c_1.96, upper_beta_c,
      std.error_beta_c, statistic_beta_c, p.value_beta_c
    )

  # note, only country_code_iso3, estimate_beta_c, std.error_beta_c, p.value_beta_c used downstream

  cat("\n--- [", group_name, "] Sample of beta^(c) estimates with confidence bounds ---\n", sep = "")
  print(head(reg_coefs %>% select(country_code_iso3, estimate_beta_c,
                                   lower_beta_c, upper_beta_c,
                                   lower_beta_c_1.96, upper_beta_c_1.96), 10))

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ############ Save output for this exclusion group ######################################
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  coefs_path <- here("output", "betas",
                      paste0("fpm_gmt_betas_country_", group_name, ".csv"))
  data_path  <- here("output", "betas",
                      paste0("regress_data_country_", group_name, ".csv"))

  # Save the coefficient table (beta_c with SEs and p-values) for this exclusion group.
  write_csv(reg_coefs, coefs_path)
  cat("[", group_name, "] Saved regression coefficients to output/betas/",
      "fpm_gmt_betas_country_", group_name, ".csv\n", sep = "")

  # Save this exclusion group's combined long-format regression input data.
  write_csv(group_data, data_path)
  cat("[", group_name, "] Saved combined regression input data to output/betas/",
      "regress_data_country_", group_name, ".csv\n", sep = "")
}

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
# Each Park model's factual and counterfactual runs are pooled into fire_model's
# classic/jules/ssib4 levels rather than getting separate intercepts, so within-group
# variation includes both
# T_ps ~= 0.25-1.1 (factual decades) and T_ps = 0 (counterfactual decades).
#
# The expanded model estimated by lm() is:
#   PM_bar = alpha_c_classic
#            + delta_jules * D_jules
#            + delta_ssib4 * D_ssib4
#            + delta_CESM  * D_CESM
#            + delta_Zhao  * D_Zhao
#            + beta_c      * T_ts
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
