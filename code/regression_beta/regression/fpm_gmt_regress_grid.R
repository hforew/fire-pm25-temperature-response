# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## FPM–GMT RELATIONSHIP: Main specification and grid-cell leave-one-model-out sensitivity for beta_i
##                        (raw fire PM2.5 concentration change per 1degC GMT)
##
## Goal: For each grid cell i, estimate the linear regression with fire-model fixed effects:
##   fPM_itsm = alpha_im + beta_i * T_ts + epsilon_itsm
##
## where fPM_itsm is fire PM2.5 concentration (µg/m^3) at cell i, time period t, scenario s,
## and fire model m; T_ts is GMT anomaly relative to the 1850-1900 pre-industrial baseline (degC);
## beta_i is the key damage function parameter (change in fPM concentration per 1degC GMT at
## cell i), and alpha_im are fire-model-specific intercepts (fixed effects absorbing
## model-level mean differences in concentration levels).
##
## This is the grid-cell counterpart to regression/fpm_gmt_regress_country.R
## (country-level): same leave-one-model-out design -- refit once for the full 5-FE-level
## model, then five exclusion groups that each drop one fire-model source's observations
## entirely -- but fit per grid cell rather than per country, and on raw fPM concentration
## rather than population-weighted per-capita exposure. Park factual/counterfactual pairs of the same model (e.g. ssib4
## factual + ssib4 counterfactual) are always dropped together, since they already share
## one FE group; CESM and Zhao are dropped alone.
##
## Five fire-model fixed effect levels (full model):
##   classic (reference), jules, ssib4 --> Park et al. factual (obsclim) AND counterfactual
##                                         (counterclim) runs of the same model pooled into
##                                         one intercept -- distinguished by T_ps (0 for
##                                         counterfactual), not by a separate FE level
##   CESM                               --> Pierce et al. projections
##   Zhao                               --> Zhao et al. projections
##
## Data spans three sources combined, up to 43 obs per cell:
##   - Pierce et al.: 5 obs (baseline + 2041-2050 RCP4.5/8.5 + 2091-2100 RCP4.5/8.5)
##   - Park et al. factual: 18 obs (3 fire models x 6 decades)
##   - Park et al. counterfactual: 18 obs (3 fire models x 6 decades, pooled into the same
##     FE group as their factual runs)
##   - Zhao et al.: 2 obs (SSP2-4.5 and SSP5-8.5 ~2095)
##   Fewer obs for each leave-one-model-out exclusion group (see exclusion_groups below).
##
## Inputs:
##   output/pm_joined/pop_pm_country_death_final.csv   (grid-cell fPM; Pierce, Zhao, and
##                                                       Park factual + counterfactual cols)
##   output/gmt/gmt_pierce_RCPs.csv          (GMT anomaly for Pierce periods × scenarios)
##   output/gmt/gmt_park_hist.csv            (GMT anomaly for each Park snapshot decade)
##   output/gmt/gmt_zhao_SSPs.csv            (GMT anomaly for Zhao ~2095 SSP245 and SSP585)
##
## Outputs (one pair per exclusion group -- full, excl_classic, excl_jules, excl_ssib4,
## excl_CESM, excl_Zhao -- written to output/betas/):
##   fpm_gmt_betas_grid_<group>.csv   (beta_i per grid cell with SE, CI, p-value, gamma_beta_i,
##                                     lon, lat, for that exclusion group)
##   regress_data_grid_<group>.csv    (that group's combined long-format regression input data)
##
## Execution order:
##   files run before: death_rate_processing.R         --> writes output/pm_joined/pop_pm_country_death_final.csv
##   files run after:  plots_beta_lobf_grid.R (reads fpm_gmt_betas_grid_<group>.csv for maps/histograms)
##
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Packages #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(here)
library(tidyverse)
library(broom)      # tidy() extracts lm() coefficients as a clean data frame

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Import #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Grid-cell fire PM2.5 data. One row per cell.
# fPM columns: Pierce (fpm_2000, fpm_2050_45/85, fpm_2100_45/85),
#              Zhao (fpm_2095_SSP245_Zhao, fpm_2095_SSP585_Zhao),
#              Park factual (park_{model}_{decade}_fpm: 18 columns),
#              Park counterfactual (park_{model}_{decade}_fpm_counter: 18 columns)
grid <- read_csv(here("output", "pm_joined", "pop_pm_country_death_final.csv"))

# GMT lookup tables — same as country-level regression
gmt_chg  <- read_csv(here("output", "gmt", "gmt_pierce_RCPs.csv"))  # Pierce periods
gmt_park <- read_csv(here("output", "gmt", "gmt_park_hist.csv"))    # Park decades
gmt_zhao <- read_csv(here("output", "gmt", "gmt_zhao_SSPs.csv"))    # Zhao ~2095

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Extract GMT scalars ##################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Baseline: average RCP4.5 and RCP8.5 for 2006-2010 — scenarios not yet diverged
gmt_baseline  <- mean(c(gmt_chg$mean_gmt_45[gmt_chg$period == "2006-2010"],
                        gmt_chg$mean_gmt_85[gmt_chg$period == "2006-2010"]))

# Pierce future scalars: one per period x scenario
gmt_2040s_45  <- gmt_chg$mean_gmt_45[gmt_chg$period == "2041-2050"]
gmt_2040s_85  <- gmt_chg$mean_gmt_85[gmt_chg$period == "2041-2050"]
gmt_2090s_45  <- gmt_chg$mean_gmt_45[gmt_chg$period == "2091-2100"]
gmt_2090s_85  <- gmt_chg$mean_gmt_85[gmt_chg$period == "2091-2100"]

# Zhao scalars: indexed by scenario string
gmt_2090s_245 <- gmt_zhao$mean_gmt_pi[gmt_zhao$scenario == "SSP245"]
gmt_2090s_585 <- gmt_zhao$mean_gmt_pi[gmt_zhao$scenario == "SSP585"]

# Park scalars: indexed by park_year (midpoint of each decade)
gmt_1960s <- gmt_park$mean_gmt_pi[gmt_park$park_year == 1965]
gmt_1970s <- gmt_park$mean_gmt_pi[gmt_park$park_year == 1975]
gmt_1980s <- gmt_park$mean_gmt_pi[gmt_park$park_year == 1985]
gmt_1990s <- gmt_park$mean_gmt_pi[gmt_park$park_year == 1995]
gmt_2000s <- gmt_park$mean_gmt_pi[gmt_park$park_year == 2005]
gmt_2010s <- gmt_park$mean_gmt_pi[gmt_park$park_year == 2015]

cat("GMT baseline (2006-2010 avg):", gmt_baseline, "\n")
cat("GMT 2040s RCP4.5:", gmt_2040s_45, "  RCP8.5:", gmt_2040s_85, "\n")
cat("GMT 2090s RCP4.5:", gmt_2090s_45, "  RCP8.5:", gmt_2090s_85, "\n")
cat("GMT 1960s:", gmt_1960s, " 1970s:", gmt_1970s, " 1980s:", gmt_1980s, "\n")
cat("GMT 1990s:", gmt_1990s, " 2000s:", gmt_2000s, " 2010s:", gmt_2010s, "\n")
cat("GMT Zhao SSP245:", gmt_2090s_245, "  SSP585:", gmt_2090s_585, "\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Reshape Pierce + Zhao fPM to long ####################################
##
## Convert from wide (one row per cell, one col per period/scenario) to long
## (one row per cell x period/scenario), then attach the GMT scalar for each row.
## Result: 7 rows per cell.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pierce_zhao_long <- grid %>%
  select(row_id, lon, lat, country_code_iso3, country_name,
         fpm_2000,
         fpm_2050_45, fpm_2050_85,
         fpm_2100_45, fpm_2100_85,
         fpm_2095_SSP245_Zhao,
         fpm_2095_SSP585_Zhao) %>%
  pivot_longer(
    cols      = starts_with("fpm_"),
    names_to  = "period_scenario",
    values_to = "fpm"
  ) %>%
  mutate(
    T_ps = case_when(
      period_scenario == "fpm_2000"             ~ gmt_baseline,
      period_scenario == "fpm_2050_45"          ~ gmt_2040s_45,
      period_scenario == "fpm_2050_85"          ~ gmt_2040s_85,
      period_scenario == "fpm_2100_45"          ~ gmt_2090s_45,
      period_scenario == "fpm_2100_85"          ~ gmt_2090s_85,
      period_scenario == "fpm_2095_SSP245_Zhao" ~ gmt_2090s_245,
      period_scenario == "fpm_2095_SSP585_Zhao" ~ gmt_2090s_585
    ),
    # fire_model: Zhao columns get "Zhao" FE; all Pierce columns get "CESM" FE
    fire_model = if_else(grepl("Zhao", period_scenario), "Zhao", "CESM")
  ) %>%
  filter(!is.na(fpm))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Reshape Park factual fPM to long ######################################
##
## Same structure as Pierce/Zhao reshape above, but for the 18 Park factual columns.
## Result: 18 rows per cell (3 fire models x 6 decades).
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

park_long <- grid %>%
  # ends_with("_fpm") excludes the "_fpm_cli" and "_fpm_counter" variants
  select(row_id, lon, lat, country_code_iso3, country_name,
         starts_with("park_") & ends_with("_fpm")) %>%
  pivot_longer(
    cols      = starts_with("park_"),
    names_to  = "period_scenario",
    values_to = "fpm"
  ) %>%
  mutate(
    T_ps = case_when(
      grepl("1960s", period_scenario) ~ gmt_1960s,
      grepl("1970s", period_scenario) ~ gmt_1970s,
      grepl("1980s", period_scenario) ~ gmt_1980s,
      grepl("1990s", period_scenario) ~ gmt_1990s,
      grepl("2000s", period_scenario) ~ gmt_2000s,
      grepl("2010s", period_scenario) ~ gmt_2010s
    ),
    fire_model = case_when(
      grepl("classic", period_scenario) ~ "classic",
      grepl("jules",   period_scenario) ~ "jules",
      grepl("ssib4",   period_scenario) ~ "ssib4"
    )
  ) %>%
  filter(!is.na(fpm))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Reshape Park counterfactual fPM to long ################################
##
## Mirrors the factual Park reshape above, using counterclim fPM. fire_model is labelled
## to match the factual reshape (e.g. "classic", not "classic_counter") so factual and
## counterfactual runs of the same Park model share one FE group / intercept below —
## they're distinguished by T_ps (0 for counterfactual) rather than by fire_model.
## Result: 18 rows per cell (3 fire models x 6 decades).
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

park_counter_long <- grid %>%
  select(row_id, lon, lat, country_code_iso3, country_name,
         starts_with("park_") & ends_with("_fpm_counter")) %>%
  pivot_longer(
    cols      = starts_with("park_"),
    names_to  = "period_scenario",
    values_to = "fpm"
  ) %>%
  mutate(
    # Counterclim runs hold CO2 (and thus GMT) fixed at 1901 levels, so all
    # counterfactual decades share T_ps = 0 (degC anomaly rel. to 1850-1900 PI baseline).
    T_ps = 0,
    fire_model = case_when(
      grepl("classic", period_scenario) ~ "classic",
      grepl("jules",   period_scenario) ~ "jules",
      grepl("ssib4",   period_scenario) ~ "ssib4"
    )
  ) %>%
  filter(!is.na(fpm))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Combine all sources #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# bind_rows stacks the three long data frames:
# 18 Park factual + 18 Park counterfactual (same fire_model labels, pooled FE group) + 7
# Pierce/Zhao = 43 rows per cell.
reg_data_combined <- bind_rows(park_long, park_counter_long, pierce_zhao_long) %>%
  mutate(fire_model = factor(fire_model,
                             levels = c("classic", "jules", "ssib4", "CESM", "Zhao")))

n_cells <- n_distinct(grid$row_id)

cat("\nCombined rows:", nrow(reg_data_combined), "\n")
cat("Expected (", n_cells, "cells x 43 obs):", n_cells * 43, "\n")

cat("\nTotal observations by fire_model (all cells combined):\n")
print(reg_data_combined %>% count(fire_model))
# Expected: classic=12*n_cells, jules=12*n_cells, ssib4=12*n_cells, CESM=5*n_cells, Zhao=2*n_cells

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Leave-one-model-out exclusion groups #########################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Sensitivity test: does beta_i depend on any single fire-model data source? "full" is
## the baseline (all 5 FE levels). Each other entry drops one fire-model source entirely
## before refitting. Since Park factual and counterfactual runs of the same model already
## share one FE group, excluding e.g. "classic" drops both its factual and counterfactual
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

# Precompute each exclusion group's filtered data (excluded fire_model level(s) dropped,
# unused factor levels dropped so lm() doesn't try to estimate an all-zero dummy column)
# and parameter count once, up front.
group_data_list <- purrr::map(exclusion_groups, function(excluded_models) {
  reg_data_combined %>%
    filter(!fire_model %in% excluded_models) %>%
    droplevels()
})

# Parameters = (FE levels - 1) dummies + intercept + slope = FE levels + 1.
n_params_list <- purrr::map_dbl(group_data_list, ~ n_distinct(.x$fire_model) + 1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Per-cell regression, confidence intervals, and save, once per exclusion group ##
##
## Model: fPM_itsm = alpha_im + beta_i * T_ts + epsilon_itsm
##
## For each grid cell i, we regress raw fPM concentration on T_ts with fire-model fixed
## effects. The full model uses all 43 observations across 5 FE levels (1 intercept + 4
## dummies + 1 slope = 6 parameters); each exclusion group below drops one source's
## observations and re-derives its own parameter count and degrees of freedom from
## however many FE levels remain.
##
## beta_i (slope)        = change in raw fire PM2.5 concentration (µg/m^3) per 1degC GMT increase.
## alpha_im (intercept)  = fire-model-specific intercept for cell i. T_ts = 0 lies outside
##   the data range (all obs at ~1degC or higher), so alpha_im is extrapolated and should not
##   be interpreted as a meaningful concentration estimate.
##
## Cells with fewer valid observations than parameters are dropped (cannot fit the model).
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for (group_name in names(exclusion_groups)) {

  # group_name <- "full"  # uncomment + set to run a single group manually, outside the loop

  group_data <- group_data_list[[group_name]]
  n_params   <- n_params_list[[group_name]]

  cat("\n============================================================\n")
  cat("[", group_name, "] Grid-cell regression, CI, and save\n", sep = "")
  cat("============================================================\n")

  # nest() collapses each cell's rows into a single list-column called "data";
  # purrr::map() applies lm() to each cell's sub-dataframe in turn; tidy() (broom)
  # extracts the coefficient table as a tidy data frame.
  reg_results <- group_data %>%
    group_by(row_id, lon, lat, country_code_iso3, country_name) %>%
    filter(n() >= n_params) %>%
    nest() %>%
    mutate(
      model  = purrr::map(data, ~ lm(fpm ~ T_ps + fire_model, data = .x)),
      tidied = purrr::map(model, tidy)
    )

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ############ Extract beta_i ###########################################################
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  reg_coefs <- reg_results %>%
    select(row_id, lon, lat, country_code_iso3, country_name, tidied) %>%
    unnest(tidied) %>%
    filter(term == "T_ps") %>%   # keep only the GMT slope; drop intercept and FE dummies
    mutate(param = "beta_i") %>%
    select(row_id, lon, lat, country_code_iso3, country_name,
           param, estimate, std.error, statistic, p.value) %>%
    pivot_wider(
      id_cols     = c(row_id, lon, lat, country_code_iso3, country_name),
      names_from  = param,
      values_from = c(estimate, std.error, statistic, p.value)
    )

  cat("[", group_name, "] Grid-cell regression complete: ", nrow(reg_coefs), " cells\n", sep = "")

  cat("\n--- [", group_name, "] Summary of beta^(i) (fPM concentration change per 1°C GMT, µg/m^3/°C) ---\n", sep = "")
  print(summary(reg_coefs$estimate_beta_i))

  cat("[", group_name, "] Cells with positive beta^(i) (more fPM with warming): ",
      sum(reg_coefs$estimate_beta_i > 0, na.rm = TRUE), "\n", sep = "")
  cat("[", group_name, "] Cells with negative beta^(i) (less fPM with warming): ",
      sum(reg_coefs$estimate_beta_i < 0, na.rm = TRUE), "\n", sep = "")

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ############ Confidence intervals for beta_i ##########################################
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Degrees of freedom depend on this exclusion group's observation count and parameter
  # count. Every retained cell has the same obs count within a group (constant 43 minus
  # whatever this group excludes), so the first fitted cell is a representative reference.
  ref_n_obs <- nrow(reg_results$data[[1]])
  df <- ref_n_obs - n_params
  t_critical <- qt(0.975, df)
  cat("\n[", group_name, "] t_critical (95% CI, df = ", df, "): ", round(t_critical, 4), "\n", sep = "")

  t_critical_largesample <- 1.96

  reg_coefs <- reg_coefs %>%
    mutate(
      lower_beta_i        = estimate_beta_i - t_critical * std.error_beta_i,
      upper_beta_i        = estimate_beta_i + t_critical * std.error_beta_i,

      lower_beta_i_1.96   = estimate_beta_i - t_critical_largesample * std.error_beta_i,
      upper_beta_i_1.96   = estimate_beta_i + t_critical_largesample * std.error_beta_i,

      t_critical_df           = t_critical,
      df_beta_i                = df,
      t_critical_largesample  = t_critical_largesample,

      # gamma: beta_i scaled by 0.008 RR (Pope et al.)
      gamma_beta_i        = 0.008 * estimate_beta_i
    ) %>%
    select(
      row_id, lon, lat, country_code_iso3, country_name, gamma_beta_i,
      t_critical_df, df_beta_i, t_critical_largesample,
      lower_beta_i, lower_beta_i_1.96,
      estimate_beta_i,
      upper_beta_i_1.96, upper_beta_i,
      std.error_beta_i, statistic_beta_i, p.value_beta_i
    )

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ############ Save output for this exclusion group ######################################
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  coefs_path <- here("output", "betas",
                      paste0("fpm_gmt_betas_grid_", group_name, ".csv"))
  data_path  <- here("output", "betas",
                      paste0("regress_data_grid_", group_name, ".csv"))

  write_csv(reg_coefs, coefs_path)
  cat("[", group_name, "] Saved regression coefficients to output/betas/",
      "fpm_gmt_betas_grid_", group_name, ".csv\n", sep = "")

  write_csv(group_data, data_path)
  cat("[", group_name, "] Saved combined regression input data to output/betas/",
      "regress_data_grid_", group_name, ".csv\n", sep = "")
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ R implementation of regression ###########################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# R implementation of fPM_itsm = alpha_im + beta_i * T_ts + epsilon_itsm:
#
# Because fire_model is a factor with classic as the reference level (unless dropped by
# an exclusion group, in which case the next remaining level becomes the reference), R
# automatically creates one indicator (dummy) variable per non-reference level, e.g. for
# the full model:
#   D_jules = 1 if fire_model == "jules", 0 otherwise
#   D_ssib4 = 1 if fire_model == "ssib4", 0 otherwise
#   D_CESM  = 1 if fire_model == "CESM",  0 otherwise
#   D_Zhao  = 1 if fire_model == "Zhao",  0 otherwise
#
# The expanded model estimated by lm() is:
#   fPM = alpha_i_classic
#         + delta_jules  * D_jules
#         + delta_ssib4  * D_ssib4
#         + delta_CESM   * D_CESM
#         + delta_Zhao   * D_Zhao
#         + beta_i       * T_ts
#         + epsilon
#
# Each Park model's factual and counterfactual runs are pooled into fire_model's
# classic/jules/ssib4 levels rather than getting separate intercepts, so within-group
# variation includes both T_ps ~= 0.25-1.1 (factual decades) and T_ps = 0 (counterfactual
# decades).
#
# alpha_im for each model is therefore:
#   classic --> (Intercept)
#   jules   --> (Intercept) + fire_modeljules
#   ssib4   --> (Intercept) + fire_modelssib4
#   CESM    --> (Intercept) + fire_modelCESM
#   Zhao    --> (Intercept) + fire_modelZhao
#
# beta_i is the T_ts coefficient — shared across all models for cell i.


### THE END
