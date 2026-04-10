# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## COMPUTE EXCESS FIRE PM MORTALITY USING POPULATION-WEIGHTED EXPOSURE AND ALPHA ##
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#
# Mortality formula applied at country level:
#   ΔMort_c(yr, scenario) = (1 - exp(-β * Δfpm_c(scenario))) * (pop_c * Y_c(yr))

rm(list = ls())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Packages #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(here)       # project-relative file paths
library(tidyverse)  # dplyr, tidyr, ggplot2
library(tibble)     # tibble utilities

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Parameters #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Same parameters as pm_mortality_orellano_cell.R — no change
RR_orellano <- 1.095                  # relative risk per 10 μg/m³ increase in PM2.5 (Orellano et al Table 1)
beta_k      <- log(RR_orellano) / 10  # β = ln(RR) / ΔX: concentration-response coefficient
             # no threshold: Orellano state no evidence for a lower safe limit

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Import #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Country-level baseline mortality (base_death_ave_2001_10) and exposure
pop_wght <- read_csv(here("output", "pop_wght_pm_cntry.csv"))

# GMT change by period and scenario (T_{p,s})
gmt_periods <- read_csv(here("output", "gmt_periods_pi.csv"))

# Country-level fPM-GMT regression coefficients (alpha_c)
alpha_coefs <- read_csv(here("output", "fpm_gmt_regression_coefs.csv"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Prepare inputs #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Retain only columns needed: country identifiers, baseline mortality, alpha_c2
# Keep only columns needed for mortality computation
mort_base <- pop_wght %>%
  select(
    country_code_iso3, country_name,
    pop_bar_c,
    starts_with("exposure_percap_fpm_"),
    base_death_ave_2001_10
  )

# Diagnose alpha coverage: which countries in pop_wght are missing from alpha_coefs?
missing_alpha <- anti_join(
  mort_base,
  alpha_coefs,
  by = "country_name"
)
print(missing_alpha %>% select(country_code_iso3, country_name))

# Join alpha_c2 — countries absent from alpha_coefs will have NA for estimate_alpha_c2
mort_base <- mort_base %>%
  left_join(
    alpha_coefs %>% select(country_name, estimate_alpha_c2),
    by = "country_name"
  )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Excess mortality #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Extract T_ps as scalars from gmt_periods
# Baseline: average of ssp245 and ssp585 for 2006-2010
T_baseline <- gmt_periods %>%
  filter(period == "2006-2010") %>%
  summarise(T = mean(c(mean_gmt_45, mean_gmt_85))) %>%
  pull(T)

T_2050_45 <- gmt_periods %>% filter(period == "2041-2050") %>% pull(mean_gmt_45)
T_2050_85 <- gmt_periods %>% filter(period == "2041-2050") %>% pull(mean_gmt_85)
T_2100_45 <- gmt_periods %>% filter(period == "2091-2100") %>% pull(mean_gmt_45)
T_2100_85 <- gmt_periods %>% filter(period == "2091-2100") %>% pull(mean_gmt_85)

# Compute excess deaths — two formulations:
# Linear approximation: β * α_c2 * T_ps * base_death_ave_2001_10
# Exact health impact function (hif):   (1 - exp(-β * α_c2 * T_ps)) * base_death_ave_2001_10
mort_base <- mort_base %>%
  mutate(
    excess_deaths_2000    = beta_k * estimate_alpha_c2 * T_baseline * base_death_ave_2001_10,
    excess_deaths_2050_45 = beta_k * estimate_alpha_c2 * T_2050_45  * base_death_ave_2001_10,
    excess_deaths_2050_85 = beta_k * estimate_alpha_c2 * T_2050_85  * base_death_ave_2001_10,
    excess_deaths_2100_45 = beta_k * estimate_alpha_c2 * T_2100_45  * base_death_ave_2001_10,
    excess_deaths_2100_85 = beta_k * estimate_alpha_c2 * T_2100_85  * base_death_ave_2001_10,

    excess_deaths_2000_hif    = (1 - exp(-beta_k * estimate_alpha_c2 * T_baseline)) * base_death_ave_2001_10,
    excess_deaths_2050_45_hif = (1 - exp(-beta_k * estimate_alpha_c2 * T_2050_45))  * base_death_ave_2001_10,
    excess_deaths_2050_85_hif = (1 - exp(-beta_k * estimate_alpha_c2 * T_2050_85))  * base_death_ave_2001_10,
    excess_deaths_2100_45_hif = (1 - exp(-beta_k * estimate_alpha_c2 * T_2100_45))  * base_death_ave_2001_10,
    excess_deaths_2100_85_hif = (1 - exp(-beta_k * estimate_alpha_c2 * T_2100_85))  * base_death_ave_2001_10,

    diff_2000    = excess_deaths_2000    - excess_deaths_2000_hif,
    diff_2050_45 = excess_deaths_2050_45 - excess_deaths_2050_45_hif,
    diff_2050_85 = excess_deaths_2050_85 - excess_deaths_2050_85_hif,
    diff_2100_45 = excess_deaths_2100_45 - excess_deaths_2100_45_hif,
    diff_2100_85 = excess_deaths_2100_85 - excess_deaths_2100_85_hif
  )


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Export #############################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

write_csv(mort_base, here("output", "fpm_excess_mortality.csv"))

glimpse(mort_base)

