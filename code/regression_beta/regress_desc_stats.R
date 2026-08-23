# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## This code provides descriptive statistics for dataset used and points estimated obtained 
## from the below regression specification 

## FPM–GMT RELATIONSHIP: Estimate beta^(c) (per-capita fire PM2.5 change per 1degC GMT)
##
## specification: For each country c, we estimate the linear regression with fire-model fixed effects:
##   PM_bar^(c)_psm = alpha^(c)_m + beta^(c) * T_ps + epsilon^(c)_psm
##   (see regression code and manuscript for detail)
##
## Data spans three sources combined:
##   - Pierce et al. projections: baseline (~2001–2010), 2041–2050, and 2091–2100
##     under RCP4.5 and RCP8.5 — 5 (period × scenario) observations per country.
##   - Park et al. historical factual: 1960s–2010s across 3 fire models (CLASSIC, JULES, SSiB4)
##     — 18 (decade × fire model) observations per country.
##   - Park et al. historical counterfactual: 1960s–2010s across the same 3 fire models
##     — 18 observations per country, pooled into the same FE group as their factual runs.
##   - Zhao et al. projections: ~2095 under SSP2-4.5 and SSP5-8.5 — 2 observations per country.
##   Combined: up to 43 observations per country for the regression.
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
library(patchwork)  # for combining ggplot2 panels into multiplots

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

# Future GMT values for each period × scenario (used as regressors T_ps)
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

# Import country-level regression coefficients (beta^(c)) estimated from the FE regression.
# Each row is a country; columns include beta^(c) estimate, SE, CI bounds, p-value, and gamma.
reg_coefs <- read_csv(here("output", "betas", "fpm_gmt_betas_country_full.csv"))

# Import combined long-format regression input data (43 obs per country).
# One row per (country, fire_model, period x scenario); columns: country_code_iso3,
# country_name, period_scenario, exposure_percap, T_ps, fire_model.
reg_data <- read_csv(here("output", "betas", "regress_data_country_full.csv")) %>%
  mutate(
    fire_model = recode(fire_model,
                        "classic" = "CLASSIC",
                        "jules"   = "JULES",
                        "ssib4"   = "SSiB4"),
    fire_model = factor(fire_model,
                        levels = c("CLASSIC", "JULES", "SSiB4", "CESM", "Zhao"))
  )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Dataset / panel structure #######################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Classify each observation as historical (Park) or future projection (Pierce/Zhao).
# Park fire models: CLASSIC, JULES, SSiB4. Pierce: CESM. Zhao: Zhao.
n_countries       <- reg_data %>% distinct(country_code_iso3) %>%
                       filter(country_code_iso3 != "global") %>% nrow()
n_obs_total       <- nrow(reg_data %>% filter(country_code_iso3 != "global"))
obs_per_country   <- n_obs_total / n_countries
n_obs_historical  <- reg_data %>%
                       filter(country_code_iso3 != "global",
                              fire_model %in% c("CLASSIC", "JULES", "SSiB4")) %>%
                       nrow() / n_countries
n_obs_future      <- reg_data %>%
                       filter(country_code_iso3 != "global",
                              fire_model %in% c("CESM", "Zhao")) %>%
                       nrow() / n_countries
n_fire_models     <- nlevels(reg_data$fire_model)

panel_table <- tibble(
  Statistic = c(
    "Countries",
    "Observations per country regression",
    "Historical observations per country (Park et al.)",
    "Future projection observations per country (Pierce + Zhao)",
    "Fire models (fixed effect levels)"
  ),
  Value = c(
    n_countries,
    obs_per_country,
    n_obs_historical,
    n_obs_future,
    n_fire_models
  )
)

cat("\n--- Dataset / panel structure ---\n")
print(panel_table, n = Inf)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Summary statistics for regression variables #####################
##
## Reports mean, SD, min, median, and max for the two core regression variables:
##   exposure_percap  (µg/m^3/yr) -- the dependent variable PM_bar^(c)_psm
##   T_ps             (degC)        -- the GMT regressor
##
## Computed across all observations in reg_data (excluding the global pseudo-country),
## so the sample matches the one actually used in the regressions.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

reg_data_cntry <- reg_data %>% filter(country_code_iso3 != "global")

summary_stats <- tibble(
  Variable = c(
    "Fire PM2.5 exposure (µg/m³/yr)",
    "GMT anomaly (°C)"
  ),
  Mean   = c(mean(reg_data_cntry$exposure_percap, na.rm = TRUE),
             mean(reg_data_cntry$T_ps,            na.rm = TRUE)),
  SD     = c(sd(reg_data_cntry$exposure_percap,   na.rm = TRUE),
             sd(reg_data_cntry$T_ps,              na.rm = TRUE)),
  Min    = c(min(reg_data_cntry$exposure_percap,  na.rm = TRUE),
             min(reg_data_cntry$T_ps,             na.rm = TRUE)),
  Median = c(median(reg_data_cntry$exposure_percap, na.rm = TRUE),
             median(reg_data_cntry$T_ps,            na.rm = TRUE)),
  Max    = c(max(reg_data_cntry$exposure_percap,  na.rm = TRUE),
             max(reg_data_cntry$T_ps,             na.rm = TRUE))
) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

cat("\n--- Summary statistics for regression variables ---\n")
print(summary_stats, n = Inf)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ By-fire-model descriptive statistics ############################
##
## Shows mean exposure, SD, and observation count for each fire-model fixed effect level.
## Large differences in mean exposure across models justify the FE specification —
## the intercept shifts absorb model-level mean differences so beta^(c) is identified
## from within-country variation in T_ps alone.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

by_fire_model <- reg_data_cntry %>%
  group_by(fire_model) %>%
  summarise(
    Mean_exposure = round(mean(exposure_percap,             na.rm = TRUE), 4),
    SD            = round(sd(exposure_percap,               na.rm = TRUE), 4),
    p10           = round(quantile(exposure_percap, 0.10,   na.rm = TRUE), 4),
    p25           = round(quantile(exposure_percap, 0.25,   na.rm = TRUE), 4),
    Median        = round(median(exposure_percap,           na.rm = TRUE), 4),
    Mean          = round(mean(exposure_percap,             na.rm = TRUE), 4),
    p75           = round(quantile(exposure_percap, 0.75,   na.rm = TRUE), 4),
    p90           = round(quantile(exposure_percap, 0.90,   na.rm = TRUE), 4),
    p95           = round(quantile(exposure_percap, 0.95,   na.rm = TRUE), 4),
    Obs           = n(),
    .groups = "drop"
  ) %>%
  rename(`Fire model` = fire_model)

cat("\n--- By-fire-model descriptive statistics ---\n")
print(by_fire_model, n = Inf)

all_obs_stats <- reg_data_cntry %>%
  summarise(
    Mean_exposure = round(mean(exposure_percap,             na.rm = TRUE), 4),
    SD            = round(sd(exposure_percap,               na.rm = TRUE), 4),
    p10           = round(quantile(exposure_percap, 0.10,   na.rm = TRUE), 4),
    p25           = round(quantile(exposure_percap, 0.25,   na.rm = TRUE), 4),
    Median        = round(median(exposure_percap,           na.rm = TRUE), 4),
    Mean          = round(mean(exposure_percap,             na.rm = TRUE), 4),
    p75           = round(quantile(exposure_percap, 0.75,   na.rm = TRUE), 4),
    p90           = round(quantile(exposure_percap, 0.90,   na.rm = TRUE), 4),
    p95           = round(quantile(exposure_percap, 0.95,   na.rm = TRUE), 4),
    Obs           = n()
  )

cat("\n--- Descriptive statistics — all observations ---\n")
print(all_obs_stats)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Historical vs future comparison #################################
##
## Documents how the three data sources partition the GMT range used for identification.
## Park historical obs anchor the lower end of T_ps; Pierce and Zhao projections
## extend it to higher warming levels, enabling identification of beta^(c).
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

by_source <- reg_data_cntry %>%
  mutate(source = case_when(
    fire_model %in% c("CLASSIC", "JULES", "SSiB4") ~ "Park historical",
    fire_model == "CESM"                            ~ "Pierce projections",
    fire_model == "Zhao"                            ~ "Zhao projections"
  )) %>%
  group_by(source) %>%
  summarise(
    Mean_exposure = round(mean(exposure_percap, na.rm = TRUE), 4),
    Mean_GMT      = round(mean(T_ps,            na.rm = TRUE), 4),
    Obs           = n(),
    .groups = "drop"
  ) %>%
  arrange(Mean_GMT) %>%
  rename(`Data source` = source)

cat("\n--- Historical vs future comparison ---\n")
print(by_source, n = Inf)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Country-level distribution of estimated slopes ##################
##
## Summarises the cross-country distribution of beta^(c) — the main output of the
## regression. Reports central tendency, spread, direction, and significance.
## Excludes the global pseudo-country row if present.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

reg_coefs_cntry <- reg_coefs %>%
  filter(country_code_iso3 != "global")

n_total    <- nrow(reg_coefs_cntry)
n_positive <- sum(reg_coefs_cntry$estimate_beta_c > 0, na.rm = TRUE)
n_sig_05   <- sum(reg_coefs_cntry$p.value_beta_c < 0.05, na.rm = TRUE)

beta_table <- tibble(
  Statistic = c(
    "Mean beta^(c)",
    "Median beta^(c)",
    "SD beta^(c)",
    "Min beta^(c)",
    "Max beta^(c)",
    "Countries with beta^(c) > 0",
    "Countries with p < 0.05"
  ),
  Value = c(
    round(mean(reg_coefs_cntry$estimate_beta_c,   na.rm = TRUE), 6),
    round(median(reg_coefs_cntry$estimate_beta_c, na.rm = TRUE), 6),
    round(sd(reg_coefs_cntry$estimate_beta_c,     na.rm = TRUE), 6),
    round(min(reg_coefs_cntry$estimate_beta_c,    na.rm = TRUE), 6),
    round(max(reg_coefs_cntry$estimate_beta_c,    na.rm = TRUE), 6),
    paste0(n_positive, " of ", n_total,
           " (", round(100 * n_positive / n_total, 1), "%)"),
    paste0(n_sig_05, " of ", n_total,
           " (", round(100 * n_sig_05 / n_total, 1), "%)")
  )
)

cat("\n--- Country-level distribution of estimated slopes ---\n")
print(beta_table, n = Inf)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Manuscript Figure 4: Box-and-whisker plots for regression variables #################
##
## Two  panels:
##   left:    Per-capita fPM2.5 exposure by fire model (Overall + 5 models)
##   right: GMT anomaly by data source (Park et al. + Pierce et al. + Zhao et al.)
##
## Box: 25th-75th percentile; whiskers: 10th-90th percentile; dot: mean.
## Pre-computing quantiles and using stat = "identity" in geom_boxplot gives full
## control over whisker endpoints (default geom_boxplot uses 1.5 * IQR).
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~ Panel A: PM exposure by fire model ~~~~

exposure_bw <- bind_rows(
  reg_data_cntry %>% mutate(group = "Overall"),
  reg_data_cntry %>% mutate(group = recode(as.character(fire_model), "Zhao" = "Multi")) # relabel Zhao --> Multi
) %>%
  group_by(group) %>%
  summarise(
    ymin   = quantile(exposure_percap, 0.10, na.rm = TRUE), # bottom whisker
    lower  = quantile(exposure_percap, 0.25, na.rm = TRUE), # bottom of box
    middle = quantile(exposure_percap, 0.50, na.rm = TRUE), # median line
    upper  = quantile(exposure_percap, 0.75, na.rm = TRUE), # top of box
    ymax   = quantile(exposure_percap, 0.90, na.rm = TRUE), # top whisker
    mean   = mean(exposure_percap,     na.rm = TRUE),
    n      = n(),
    .groups = "drop"
  ) %>%
  mutate(group = factor(group,                              # fix display order
                        levels = c("Overall", "CLASSIC", "JULES", "SSiB4", "CESM", "Multi")))

p_exposure <- ggplot(exposure_bw, aes(x = group)) +
  geom_boxplot(aes(ymin = ymin, lower = lower, middle = middle,
                   upper = upper, ymax = ymax),
               stat = "identity",                          # use pre-computed quantiles, not raw data
               fill = "steelblue", alpha = 0.5, width = 0.5) +
  geom_point(aes(y = mean), shape = 19, size = 2.5, color = "darkred") +
  geom_text(aes(y = Inf, label = paste0("n=", n)),         # Inf pins label to top of panel
            vjust = 1.5, hjust = 0.5, size = 3) +
  labs(x = NULL, y = "Per-capita Fire PM2.5 (µg/m³/yr)",
       title = "A. Per-capita Fire PM2.5 exposure by fire model") +
  theme_bw(base_size = 11)

# ~~~~ Panel B: GMT anomaly by data source (Cleveland dot plot) ~~~~
# T_ps values are fixed scalars assigned by construction (7 Park, 5 Pierce, 2 Zhao),
# not a continuous distribution. A dot plot shows the actual discrete values honestly.

gmt_dots <- reg_data_cntry %>%
  mutate(source = case_when(
    fire_model %in% c("CLASSIC", "JULES", "SSiB4") ~ "Park et al.",
    fire_model == "CESM"                            ~ "Pierce et al.",
    fire_model == "Zhao"                            ~ "Zhao et al."
  )) %>%
  distinct(source, T_ps) %>%          # one row per unique GMT value per source
  mutate(source = factor(source,      # fix display order
                         levels = c("Park et al.", "Pierce et al.", "Zhao et al.")))

p_gmt <- ggplot(gmt_dots, aes(x = source, y = T_ps)) +
  geom_point(size = 3, color = "steelblue") +
  labs(x = NULL, y = "GMT anomaly (°C)",
       title = "B. GMT anomaly by data source") +
  theme_bw(base_size = 11)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Manuscript Figure 4 note (for paper) ################################################
##
## Manuscript Figure 4. Distribution of regression inputs across fire models and data sources.
## Panel A shows the cross-country distribution of population-weighted per-capita
## fire PM2.5 exposure (µg/m³/yr) for the full sample (Overall) and separately for
## each of the five fire-model fixed-effect levels: CLASSIC, JULES, and SSiB4
## (Park et al. historical), CESM (Pierce et al. projections), and Zhao
## (Zhao et al. projections).
## Panel B shows the discrete GMT anomaly values (degC relative to the 1850–1900
## pre-industrial baseline) used as regressors, by data source: 7 for Park et al.
## (6 historical decade means plus one shared T_ps = 0 point from counterfactual
## runs), 5 for Pierce et al. (baseline and four period x scenario combinations),
## and 2 for Zhao et al. (SSP2-4.5 and SSP5-8.5). Each dot is one unique T_ps value.
## In Panel A, the box spans the 25th–75th percentile, the center line is the
## median, whiskers extend to the 10th and 90th percentiles, and the red dot
## marks the mean.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
######  Combine and save ##############
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

caption_text <- "Panel A — box: 25th–75th pct; center line: median; whiskers: 10th–90th pct; dot: mean. Panel B — each dot is unique GMT value."

fig_wide <- (p_exposure | p_gmt) +
  plot_annotation(caption = caption_text)

print(fig_wide)

ggsave(here("images", "regression_beta", "beta_country",
            "fig4_desc_stats.png"),
       fig_wide, width = 8.5, height = 4, dpi = 300)
cat("Saved fig4_desc_stats to images/regression_beta/beta_country/fig4_desc_stats.png\n")


# THE END
