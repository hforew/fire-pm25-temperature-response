# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Descriptive statistics for the dataset and beta^(c) estimates from the FE regression.
##
## FPM–GMT RELATIONSHIP: Estimate beta^(c) (per-capita fire PM2.5 change per 1°C GMT)
##
## Goal: For each country c, we estimate the linear regression with fire-model fixed effects:
##   PM_bar^(c)_psm = alpha^(c)_m + beta^(c) * T_ps + epsilon^(c)_psm
##
## where PM_bar^(c)_psm is population-weighted per-capita fire PM2.5 exposure (µg/m^3/yr)
## for country c, period p, scenario s, and fire model m;
## T_ps is GMT anomaly relative to the 1850-1900 pre-industrial baseline (°C);
## beta^(c) is the key damage function parameter (change in population-weighted
## per-capita fire PM2.5 per 1°C GMT), and alpha^(c)_m are fire-model-specific
## intercepts (fixed effects absorbing model-level mean differences in exposure levels).
##
## Five fire-model fixed effect levels:
##   CLASSIC (reference), JULES, SSiB4  --> Park et al. historical models
##   CESM                               --> Pierce et al. projections
##   Zhao                               --> Zhao et al. projections
##
## Data spans three sources combined:
##   - Pierce et al. projections: baseline (~2001–2010), 2041–2050, and 2091–2100
##     under RCP4.5 and RCP8.5 — 5 (period × scenario) observations per country.
##   - Park et al. historical: 1960s–2010s across 3 fire models (CLASSIC, JULES, SSiB4)
##     — 18 (decade × fire model) observations per country; climate change-attributable output.
##   - Zhao et al. projections: ~2095 under SSP2-4.5 and SSP5-8.5 — 2 observations per country.
##   Combined: up to 25 observations per country for the regression.
##
## Inputs:
##   output/fpm_gmt_regression_coefs_FE_all_cli.csv  (beta^(c) per country; from fpm_gmt_regression_FE_all_cli.R)
##   output/reg_data_combined_fe_all_cli.csv          (combined long-format regression input; from same)
##
## Outputs:
##   images/regression_alpha/beta_FE_all_cli/bw_exposure_gmt_by_model_long.png
##   images/regression_alpha/beta_FE_all_cli/bw_exposure_gmt_by_model_wide.png
##
## Execution order:
##   files run before: fpm_gmt_regression_FE_all_cli.R   --> writes both input CSVs above
##   files run after: NA
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Remove all objects from the environment to start fresh
rm(list = ls())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Packages #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(here)       # for portable file paths relative to project root
library(tidyverse)  # for data manipulation (dplyr, tidyr, purrr) and ggplot2
library(tibble)     # for tidy data frames (part of tidyverse, loaded explicitly for clarity)
library(patchwork)  # for combining ggplot2 panels into multiplots

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Import #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Import country-level regression coefficients (beta^(c)) estimated from the FE regression.
# Each row is a country; columns include beta^(c) estimate, SE, CI bounds, p-value, and gamma.
reg_coefs <- read_csv(here("output", "fpm_gmt_regression_coefs_FE_all_cli.csv"))

# Import combined long-format regression input data (25 obs per country).
# One row per (country, fire_model, period x scenario); columns: country_code_iso3,
# country_name, period_scenario, exposure_percap, T_ps, fire_model.
reg_data <- read_csv(here("output", "reg_data_combined_fe_all_cli.csv")) %>%
  mutate(
    fire_model = recode(fire_model,
                        "classic" = "CLASSIC",
                        "jules"   = "JULES",
                        "ssib4"   = "SSiB4"),
    fire_model = factor(fire_model,
                        levels = c("CLASSIC", "JULES", "SSiB4", "CESM", "Zhao"))
  )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Table 1: Dataset / panel structure #######################################
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

cat("\n--- Table 1: Dataset / panel structure ---\n")
print(panel_table, n = Inf)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Table 2: Summary statistics for regression variables #####################
##
## Reports mean, SD, min, median, and max for the two core regression variables:
##   exposure_percap  (µg/m^3/yr) -- the dependent variable PM_bar^(c)_psm
##   T_ps             (°C)        -- the GMT regressor
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

cat("\n--- Table 2: Summary statistics for regression variables ---\n")
print(summary_stats, n = Inf)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Table 3: By-fire-model descriptive statistics ############################
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

cat("\n--- Table 3: By-fire-model descriptive statistics ---\n")
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

cat("\n--- Table 3b: Descriptive statistics — all observations ---\n")
print(all_obs_stats)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Table 4: Historical vs future comparison #################################
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

cat("\n--- Table 4: Historical vs future comparison ---\n")
print(by_source, n = Inf)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Table 5: Country-level distribution of estimated slopes ##################
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

cat("\n--- Table 5: Country-level distribution of estimated slopes ---\n")
print(beta_table, n = Inf)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Figure 1: Box-and-whisker plots for regression variables #################
##
## Two stacked panels:
##   Top:    Per-capita fPM2.5 exposure by fire model (Overall + 5 models)
##   Bottom: GMT anomaly by data source (Overall + Park grouped + CESM + Zhao)
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
# T_ps values are fixed scalars assigned by construction (6 Park, 5 Pierce, 2 Zhao),
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
############ Figure 1 note (for paper) ################################################
##
## Figure 1. Distribution of regression inputs across fire models and data sources.
## Panel A shows the cross-country distribution of population-weighted per-capita
## fire PM2.5 exposure (µg/m³/yr) for the full sample (Overall) and separately for
## each of the five fire-model fixed-effect levels: CLASSIC, JULES, and SSiB4
## (Park et al. historical), CESM (Pierce et al. projections), and Zhao
## (Zhao et al. projections).
## Panel B shows the discrete GMT anomaly values (°C relative to the 1850–1900
## pre-industrial baseline) used as regressors, by data source: 6 values for
## Park et al. (1960s–2010s decade means), 5 for Pierce et al. (baseline and
## four period x scenario combinations), and 2 for Zhao et al. (SSP2-4.5 and
## SSP5-8.5). Each dot is one unique T_ps value.
## In Panel A, the box spans the 25th–75th percentile, the center line is the
## median, whiskers extend to the 10th and 90th percentiles, and the red dot
## marks the mean.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
######  Combine and save ##############
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

caption_text <- "Panel A — box: 25th–75th pct; center line: median; whiskers: 10th–90th pct; dot: mean. Panel B — each dot is unique GMT value."

fig_long <- (p_exposure / p_gmt) +
  plot_annotation(caption = caption_text)

fig_wide <- (p_exposure | p_gmt) +
  plot_annotation(caption = caption_text)

print(fig_long)
print(fig_wide)

# ggsave(here("images", "regression_alpha", "beta_FE_all_cli",
#             "bw_exposure_gmt_by_model_long.png"),
#        fig_long, width = 8.5, height = 10, dpi = 300)
# cat("\nSaved _long to images/regression_alpha/beta_FE_all_cli/bw_exposure_gmt_by_model_long.png\n")

ggsave(here("images", "regression_alpha", "beta_FE_all_cli",
            "bw_exposure_gmt_by_model_wide.png"),
       fig_wide, width = 8.5, height = 4, dpi = 300)
cat("Saved _wide to images/regression_alpha/beta_FE_all_cli/bw_exposure_gmt_by_model_wide.png\n")


# THE END
