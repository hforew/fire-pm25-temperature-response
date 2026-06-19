# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## FPM–GMT RELATIONSHIP: Grid-cell level fixed-effect regression
##
## Goal: For each grid cell i, estimate the linear regression with fire-model fixed effects:
##   fPM_itsm = alpha_im + beta_i * T_ts + epsilon_itsm
##
## where fPM_itsm is fire PM2.5 concentration (µg/m^3) at cell i, time period t, scenario s,
## and fire model m; T_ts is GMT anomaly relative to the 1850-1900 pre-industrial baseline (°C);
## beta_i is the slope (change in fPM concentration per 1°C GMT at cell i), and
## alpha_im are fire-model-specific intercepts (fixed effects absorbing model-level
## mean differences in concentration levels).
##
## Five fire-model fixed effect levels:
##   classic (reference), jules, ssib4  --> Park et al. historical models
##   CESM                               --> Pierce et al. projections
##   Zhao                               --> Zhao et al. projections
##
## Data: same three sources as country-level regression, 25 obs per cell:
##   - Pierce et al.: 5 obs (baseline + 2041-2050 RCP4.5/8.5 + 2091-2100 RCP4.5/8.5)
##   - Park et al.:  18 obs (3 fire models x 6 decades)
##   - Zhao et al.:   2 obs (SSP2-4.5 and SSP5-8.5 ~2095)
##
## Output: one row per grid cell with beta_i, SE, p-value, lon, lat.
##         Used for spatial mapping in a separate plotting script.
##
## Note: response is raw fPM concentration (µg/m^3), not per-capita exposure.
##       Grid-cell betas are NOT aggregated to country level here.
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
#              Park (park_{model}_{decade}_fpm: 18 columns)
grid <- read_csv(here("output", "pop_pm_country_death_final.csv"))

# GMT lookup tables — same as country-level regression
gmt_chg  <- read_csv(here("output", "gmt_periods_pi.csv"))   # Pierce periods
gmt_park <- read_csv(here("output", "gmt_park_decades.csv")) # Park decades
gmt_zhao <- read_csv(here("output", "gmt_zhao_mimi.csv"))    # Zhao ~2095

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
############ Descriptive Statistics: fPM distribution and zero-cell counts ###########
##
## Two tables for all 25 fPM source columns (7 Pierce/Zhao + 18 Park):
##   fpm_desc_stats: min, 25th pctile, mean, 75th pctile, max across grid cells
##   zero_fpm:       count and % of cells with fpm == 0 per column
##
## Run before reshaping to confirm data coverage and identify any degenerate columns
## (all-zero fPM) that would produce uninformative regression slopes.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Park fPM source columns: 18 total (3 fire models x 6 decades)
# grep() matches column names against a regex pattern and returns the matching names
park_fpm_cols <- grep("^park_(classic|jules|ssib4)_.+_fpm$", colnames(grid), value = TRUE)

# Pierce + Zhao fPM source columns: 7 total
pierce_zhao_fpm_cols <- c(
  "fpm_2000",
  "fpm_2050_45", "fpm_2050_85",
  "fpm_2100_45", "fpm_2100_85",
  "fpm_2095_SSP245_Zhao", "fpm_2095_SSP585_Zhao"
)

n_cells <- nrow(grid)

# all 25 fPM columns in order: Pierce/Zhao first, then Park
all_fpm_cols <- c(pierce_zhao_fpm_cols, park_fpm_cols)

# across(everything(), list(...)) applies each named function to every selected column,
# producing wide output with names like "fpm_2000_min", "fpm_2000_mean", etc.
# pivot_longer then collapses those wide stat columns back to long, using names_pattern
# to split "col_stat" into two pieces: the original column name and the stat name.
# .value in names_to tells pivot_longer to use the second capture group as column names,
# so each stat (min, p25, mean, p75, max) becomes its own output column.
fpm_desc_stats <- grid %>%
  select(all_of(all_fpm_cols)) %>%
  summarise(across(everything(), list(
    min  = ~ min(.x,            na.rm = TRUE),
    p25  = ~ quantile(.x, 0.25, na.rm = TRUE),
    mean = ~ mean(.x,           na.rm = TRUE),
    p75  = ~ quantile(.x, 0.75, na.rm = TRUE),
    max  = ~ max(.x,            na.rm = TRUE)
  ))) %>%
  pivot_longer(everything(),
               names_to      = c("column", ".value"),       # split name into col + stat
               names_pattern = "^(.+)_(min|p25|mean|p75|max)$") %>%  # regex captures both
  mutate(
    # label each row with its data source for readability
    source = case_when(
      grepl("^park_", column) ~ "Park",
      grepl("Zhao",   column) ~ "Zhao",
      TRUE                    ~ "Pierce (CESM)"
    ),
    across(c(min, p25, mean, p75, max), ~ round(.x, 3))     # round all stat cols at once
  ) %>%
  arrange(source, column) %>%
  select(source, column, min, p25, mean, p75, max)

cat("\n--- fPM descriptive statistics by model/scenario (n_cells =", n_cells, ") ---\n")
print(fpm_desc_stats, n = Inf)

# count cells where fpm == 0; pct_zero flags columns with widespread no-fire coverage
# across() with a single anonymous function applies it to every selected column at once
zero_fpm <- grid %>%
  select(all_of(all_fpm_cols)) %>%
  summarise(across(everything(), ~ sum(.x == 0, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "column", values_to = "n_zero") %>%
  mutate(
    source   = case_when(
      grepl("^park_", column) ~ "Park",
      grepl("Zhao",   column) ~ "Zhao",
      TRUE                    ~ "Pierce (CESM)"
    ),
    pct_zero = round(100 * n_zero / n_cells, 1)
  ) %>%
  arrange(source, column) %>%
  select(source, column, n_zero, pct_zero)

cat("\n--- Zero-fPM cell counts by model/scenario (n_cells =", n_cells, ") ---\n")
print(zero_fpm, n = Inf)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Reshape Pierce + Zhao fPM to long ####################################
##
## Convert from wide (one row per cell, one col per period/scenario) to long
## (one row per cell x period/scenario), then attach the GMT scalar for each row.
## Result: 7 rows per cell, columns: row_id, lon, lat, period_scenario, fpm, T_ts (code variable: T_ps),
##         fire_model.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pierce_zhao_long <- grid %>%
  select(row_id, lon, lat, country_code_iso3, country_name,
         fpm_2000,
         fpm_2050_45, fpm_2050_85,
         fpm_2100_45, fpm_2100_85,
         fpm_2095_SSP245_Zhao,
         fpm_2095_SSP585_Zhao) %>%
  # pivot_longer stacks the 7 fpm columns into two columns:
  # period_scenario (column name) and fpm (cell value)
  pivot_longer(
    cols      = starts_with("fpm_"),
    names_to  = "period_scenario",  # column name becomes a label
    values_to = "fpm"               # cell value becomes the response variable
  ) %>%
  mutate(
    # T_ps: assign the GMT scalar that corresponds to each period x scenario label
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
  filter(!is.na(fpm))  # drop cells with no fPM data for this period/scenario

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Reshape Park fPM to long #############################################
##
## Same structure as Pierce/Zhao reshape above, but for 18 Park columns.
## Result: 18 rows per cell (3 fire models x 6 decades).
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

park_long <- grid %>%
  # select only the 18 Park fPM columns (excludes park_*_pop and other non-fpm cols)
  select(row_id, lon, lat, country_code_iso3, country_name,
         starts_with("park_") & ends_with("_fpm")) %>%
  # pivot_longer stacks the 18 park fpm columns into period_scenario + fpm
  pivot_longer(
    cols      = starts_with("park_"),
    names_to  = "period_scenario",
    values_to = "fpm"
  ) %>%
  mutate(
    # T_ps: decade is encoded in the column name, so grepl() extracts it
    # without needing 18 explicit case_when entries
    T_ps = case_when(
      grepl("1960s", period_scenario) ~ gmt_1960s,
      grepl("1970s", period_scenario) ~ gmt_1970s,
      grepl("1980s", period_scenario) ~ gmt_1980s,
      grepl("1990s", period_scenario) ~ gmt_1990s,
      grepl("2000s", period_scenario) ~ gmt_2000s,
      grepl("2010s", period_scenario) ~ gmt_2010s
    ),
    # fire_model: land-surface model also encoded in column name
    fire_model = case_when(
      grepl("classic", period_scenario) ~ "classic",
      grepl("jules",   period_scenario) ~ "jules",
      grepl("ssib4",   period_scenario) ~ "ssib4"
    )
  ) %>%
  filter(!is.na(fpm))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Combine and run cell-level regression ################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# bind_rows stacks the two long data frames: 18 Park rows + 7 Pierce/Zhao rows = 25 per cell
reg_data <- bind_rows(park_long, pierce_zhao_long) %>%
  # factor() with explicit levels sets "classic" as the reference category;
  # lm() will create dummies for jules, ssib4, CESM, Zhao relative to classic
  mutate(fire_model = factor(fire_model,
                             levels = c("classic", "jules", "ssib4", "CESM", "Zhao")))

cat("\nCombined rows:", nrow(reg_data), "\n")
cat("Expected (", n_cells, "cells x 25 obs):", n_cells * 25, "\n")

# nest() collapses each cell's 25 rows into a single list-column called "data",
# giving one row per cell in reg_results.
# purrr::map() then applies lm() to each cell's sub-dataframe in turn.
# tidy() (from broom) extracts the coefficient table as a tidy data frame.
# Cells with fewer than 6 obs are dropped — cannot fit 6 parameters
# (intercept + T_ps + 4 fire-model dummies). In practice all cells have 25 obs.
# See last section for how R actually implements the regression.
reg_results <- reg_data %>%
  group_by(row_id, lon, lat, country_code_iso3, country_name) %>%
  filter(n() >= 6) %>%
  nest() %>%
  mutate(
    model  = purrr::map(data, ~ lm(fpm ~ T_ps + fire_model, data = .x)),
    tidied = purrr::map(model, tidy)  # one tidy coef table per cell
  )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Extract beta_i ########################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# unnest(tidied) expands the list-column so each coefficient term gets its own row.
# Only the T_ps slope term is retained — the CLASSIC intercept and FE dummies
# are not stored in the output.
# pivot_wider rotates beta_i from rows into columns, producing one output row
# per cell with estimate, SE, t-stat, and p-value.
reg_coefs <- reg_results %>%
  select(row_id, lon, lat, country_code_iso3, country_name, tidied) %>%
  unnest(tidied) %>%           # one row per term per cell
  filter(term == "T_ps") %>%   # keep only the GMT slope; drop intercept and FE dummies
  mutate(param = "beta_i") %>%
  select(row_id, lon, lat, country_code_iso3, country_name,
         param, estimate, std.error, statistic, p.value) %>%
  pivot_wider(
    id_cols     = c(row_id, lon, lat, country_code_iso3, country_name),
    names_from  = param,                                          # beta_i becomes col prefix
    values_from = c(estimate, std.error, statistic, p.value)     # four cols per param
  )

cat("\nGrid-cell regression complete:", nrow(reg_coefs), "cells\n")

cat("\n--- Summary of beta^(i) (fPM concentration change per 1°C GMT, µg/m^3/°C) ---\n")
print(summary(reg_coefs$estimate_beta_i))

cat("\nCells with positive beta^(i) (more fPM with warming):",
    sum(reg_coefs$estimate_beta_i > 0, na.rm = TRUE), "\n")
cat("Cells with negative beta^(i) (less fPM with warming):",
    sum(reg_coefs$estimate_beta_i < 0, na.rm = TRUE), "\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Export ##############################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

write_csv(reg_coefs, here("output", "fpm_gmt_regression_coefs_FE_grid.csv"))
cat("\nSaved to output/fpm_gmt_regression_coefs_FE_grid.csv\n")

glimpse(reg_coefs)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ R implementation of regression ###########################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# R implementation of fPM_itsm = alpha_im + beta_i * T_ts + epsilon_itsm:
#
# Because fire_model is a factor with classic as the reference level, R automatically
# creates four indicator (dummy) variables — one for each non-reference model:
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
# alpha_im for each model is therefore:
#   classic --> (Intercept)
#   jules   --> (Intercept) + fire_modeljules
#   ssib4   --> (Intercept) + fire_modelssib4
#   CESM    --> (Intercept) + fire_modelCESM
#   Zhao    --> (Intercept) + fire_modelZhao
#
# beta_i is the T_ts coefficient — shared across all models for cell i.


### THE END
