# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## FPM–GMT RELATIONSHIP: Estimate beta_c (per-capita fire PM2.5 change per 1°C GMT)
##                        Park et al. historical data — one regression per fire model
##
## Goal: For each country c, estimate the linear regression:
##   PM_bar_ct = alpha_c + beta_c * T_t + epsilon_ct
##
## where T_t is GMT change in decade t (°C relative to 1850–1900),
## and beta_c is the key damage function parameter.
## Estimated separately for each Park fire model m (classic, jules, ssib4) —
## m is not a subscript in the equation above since each regression call
## uses only one model's data (no fire-model fixed effect, unlike fpm_gmt_regression_FE_all.R).
##
## Approach mirrors fpm_gmt_regression_pierce.R: simple OLS with no fixed effects,
## run separately for each fire model. 6 observations per country per model
## (one per decade: 1960s–2010s).
##
## Inputs:
##   output/pop_wght_pm_cntry_final.csv   (country-level per-capita fPM exposure;
##                                         Park columns: 3 fire models x 6 decades)
##   output/gmt/gmt_park_hist.csv         (GMT anomaly for each Park snapshot decade,
##                                         relative to 1850–1900 PI baseline)
##
## Outputs:
##   output/betas_single_study/fpm_gmt_regression_coefs_classic.csv  (beta_c for classic model)
##   output/betas_single_study/fpm_gmt_regression_coefs_jules.csv    (beta_c for JULES model)
##   output/betas_single_study/fpm_gmt_regression_coefs_ssib4.csv    (beta_c for SSIB4 model)
##
## Execution order:
##   files run before: pop_wght_pm_cntry.R --> writes pop_wght_pm_cntry_final.csv
##                     temp_chg.R          --> writes gmt_park_hist.csv
##   files run after: beta_comparison_latex.R --> reads all three coefficient outputs above
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
# Each row is a country; Park columns cover 3 fire models × 6 decades (1960s–2010s).
pop_wght <- read_csv(here("output", "pop_wght_pm_cntry_final.csv"))

# Drop GIVE countries with no corresponding data in the Park PM gridded data.
pop_wght <- pop_wght %>%
  filter(pop_bar_c != 0)

# Import Park decade GMT values (°C relative to pre-industrial baseline)
# Rows: one per Park snapshot decade (1960s–2010s)
# Columns: park_year, decade, mean_gmt_pi
gmt_park <- read_csv(here("output", "gmt", "gmt_park_hist.csv"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Extract scalar GMT values for each Park decade ##############################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Park decade GMT scalars — one per decade, indexed by park_year (reliable numeric).
gmt_1960s <- gmt_park$mean_gmt_pi[gmt_park$park_year == 1965]   # 1960s decade mean GMT
gmt_1970s <- gmt_park$mean_gmt_pi[gmt_park$park_year == 1975]   # 1970s decade mean GMT
gmt_1980s <- gmt_park$mean_gmt_pi[gmt_park$park_year == 1985]   # 1980s decade mean GMT
gmt_1990s <- gmt_park$mean_gmt_pi[gmt_park$park_year == 1995]   # 1990s decade mean GMT
gmt_2000s <- gmt_park$mean_gmt_pi[gmt_park$park_year == 2005]   # 2000s decade mean GMT
gmt_2010s <- gmt_park$mean_gmt_pi[gmt_park$park_year == 2015]   # 2010s decade mean GMT

cat("GMT 1960s:", gmt_1960s, " 1970s:", gmt_1970s, " 1980s:", gmt_1980s, "\n")
cat("GMT 1990s:", gmt_1990s, " 2000s:", gmt_2000s, " 2010s:", gmt_2010s, "\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Reshape Park data to long format ###########################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# One row per (country, fire_model × decade).
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
park_data_long <- park_data_wide %>%
  pivot_longer(
    cols      = starts_with("exposure_percap_park_"),
    names_to  = "period_scenario",
    values_to = "exposure_percap"
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
############ Helper: run Pierce-style regression for one Park fire model ################
##
## Runs simple OLS (exposure_percap ~ T_ps) per country across 6 historical decades.
## Includes USA and Global diagnostics, then the full nested regression.
## Returns a tibble matching the pierce output column layout (t_critical_df4 replaces
## t_critical_df3 to reflect df = 6 - 2 = 4).
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

run_park_regression <- function(model_data, model_label) {

  # Degrees of freedom: n = 6 observations, 2 parameters (intercept + slope) --> df = 4
  df         <- nrow(model_data %>% filter(country_code_iso3 == "USA")) - 2
  t_critical <- qt(0.975, df)   # ~2.776 for 95% CI with 4 degrees of freedom
  t_ls       <- 1.96

  cat("\n=======================================================\n")
  cat(" Fire model:", model_label, "\n")
  cat("=======================================================\n")
  cat("t_critical (95% CI, df =", df, "):", round(t_critical, 4), "\n")

  # --- USA diagnostic ---
  usa_data  <- model_data %>% filter(country_code_iso3 == "USA")
  cat("\n--- USA regression input data (", model_label, ") ---\n")
  print(usa_data %>% select(period_scenario, T_ps, exposure_percap))
  usa_model <- lm(exposure_percap ~ T_ps, data = usa_data)        # OLS: alpha_c + beta_c*T_t
  cat("\n--- USA lm() summary (", model_label, ") ---\n")
  print(summary(usa_model))                                        # coefs, SE, t-stat, p-value, R^2, F-stat
  cat("USA alpha_c:", round(coef(usa_model)[["(Intercept)"]], 6), "µg/m³/yr\n")
  cat("USA beta_c (slope, fPM change per 1°C GMT):", round(coef(usa_model)[["T_ps"]], 6), "µg/m³/yr per °C\n")
  cat("USA R-squared:", round(summary(usa_model)$r.squared, 4), "\n")

  # --- Global diagnostic ---
  global_data  <- model_data %>% filter(country_name == "global")
  cat("\n--- Global regression input data (", model_label, ") ---\n")
  print(global_data %>% select(period_scenario, T_ps, exposure_percap))
  global_model <- lm(exposure_percap ~ T_ps, data = global_data)  # OLS: alpha_c + beta_c*T_t
  cat("\n--- Global lm() summary (", model_label, ") ---\n")
  print(summary(global_model))                                     # coefs, SE, t-stat, p-value, R^2, F-stat
  cat("Global alpha_c:", round(coef(global_model)[["(Intercept)"]], 6), "µg/m³/yr\n")
  cat("Global beta_c (slope, fPM change per 1°C GMT):", round(coef(global_model)[["T_ps"]], 6), "µg/m³/yr per °C\n")
  cat("Global R-squared:", round(summary(global_model)$r.squared, 4), "\n")

  # --- Full country-level regression ---
  reg_results <- model_data %>%
    group_by(country_code_iso3, country_name) %>%
    filter(n() >= 2) %>%                                                    # skip countries with < 2 obs
    nest() %>%                                                               # one list-column of rows per country
    mutate(
      model   = purrr::map(data, ~ lm(exposure_percap ~ T_ps, data = .x)), # OLS per country
      tidied  = purrr::map(model, tidy),                                    # coefs as tibble (term, estimate, SE, t, p)
      glanced = purrr::map(model, glance)                                   # model-level stats (R^2, F, etc.)
    )

  reg_coefs <- reg_results %>%
    select(country_code_iso3, country_name, tidied) %>%
    unnest(tidied) %>%                                        # one row per (country, term)
    mutate(param = case_when(
      term == "(Intercept)" ~ "alpha_c",
      term == "T_ps"        ~ "beta_c"
    )) %>%
    select(country_code_iso3, country_name, param, estimate, std.error, statistic, p.value) %>%
    pivot_wider(                                              # one row per country, alpha_c and beta_c as columns
      id_cols     = c(country_code_iso3, country_name),
      names_from  = param,
      values_from = c(estimate, std.error, statistic, p.value)
    )

  print(head(reg_coefs, 10))
  cat("\nDimensions of regression coefficient table:", nrow(reg_coefs), "countries\n")

  cat("\n--- Summary of beta_c (slope: per-capita fPM2.5 change per 1°C GMT) ---\n")
  print(summary(reg_coefs$estimate_beta_c))
  n_total <- nrow(reg_coefs)
  n_pos   <- sum(reg_coefs$estimate_beta_c > 0, na.rm = TRUE)
  cat(sprintf("beta > 0 for %d of %d countries (%.0f%%)\n",
              n_pos, n_total, 100 * n_pos / n_total))

  # Add confidence intervals and gamma * beta_c
  reg_coefs <- reg_coefs %>%
    mutate(
      # 95% CI bounds using exact t_critical for df = 4
      lower_beta_c           = estimate_beta_c - t_critical * std.error_beta_c,
      upper_beta_c           = estimate_beta_c + t_critical * std.error_beta_c,

      # 95% CI bounds using large-sample approximation (1.96)
      lower_beta_c_1.96      = estimate_beta_c - t_ls * std.error_beta_c,
      upper_beta_c_1.96      = estimate_beta_c + t_ls * std.error_beta_c,

      # Store both critical values as columns for reference
      t_critical_df4         = t_critical,
      t_critical_largesample = t_ls,

      # gamma * beta_c: scaled slope used in the GIVE damage function. .008 from Pope HR
      gamma_beta_c           = 0.008 * estimate_beta_c
    ) %>%
    select(
      country_code_iso3, country_name, gamma_beta_c,
      t_critical_df4, t_critical_largesample,
      lower_beta_c, lower_beta_c_1.96,
      estimate_beta_c,
      upper_beta_c_1.96, upper_beta_c,
      std.error_beta_c, statistic_beta_c, p.value_beta_c,
      estimate_alpha_c, std.error_alpha_c, statistic_alpha_c, p.value_alpha_c
    )

  countries_to_plot <- c("global", "USA", "DEU", "RUS", "AUS", "CHN", "IND", "ARG", "BRA", "ETH", "NGA")

  cat("\nSample of beta_c estimates with confidence bounds:\n")
  print(reg_coefs %>%
          filter(country_code_iso3 %in% countries_to_plot) %>%
          select(country_code_iso3, estimate_beta_c,
                 lower_beta_c, upper_beta_c,
                 lower_beta_c_1.96, upper_beta_c_1.96))

  reg_coefs
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Run regression for each Park fire model ####################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

reg_coefs_classic <- run_park_regression(park_data_long %>% filter(fire_model == "classic"), "classic")
reg_coefs_jules   <- run_park_regression(park_data_long %>% filter(fire_model == "jules"),   "jules")
reg_coefs_ssib4   <- run_park_regression(park_data_long %>% filter(fire_model == "ssib4"),   "ssib4")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Save output ##############################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

write_csv(reg_coefs_classic, here("output", "betas_single_study", "fpm_gmt_regression_coefs_classic.csv"))
write_csv(reg_coefs_jules,   here("output", "betas_single_study", "fpm_gmt_regression_coefs_jules.csv"))
write_csv(reg_coefs_ssib4,   here("output", "betas_single_study", "fpm_gmt_regression_coefs_ssib4.csv"))

cat("\nSaved regression coefficients to:\n")
cat("  output/betas_single_study/fpm_gmt_regression_coefs_classic.csv\n")
cat("  output/betas_single_study/fpm_gmt_regression_coefs_jules.csv\n")
cat("  output/betas_single_study/fpm_gmt_regression_coefs_ssib4.csv\n")


# THE END
