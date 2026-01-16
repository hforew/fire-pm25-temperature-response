# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## COMPUTE MORTALITY ##
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Remove all objects from the environment
rm(list = ls())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Packages #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(here)
library(tidyverse)
library(tibble)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Import #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Import all data
pop_pm_country <- read_csv(here("output", "pop_pm_with_countries.csv"))

colnames(pop_pm_country)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Mortality Equation Setup #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# SHOULD BE ANNUAL BUT CURRENTLY DECADE AVERAGE 

# this code computes mortality from total PM and fire PM
# ΔMortality = Pop (1 - exp^(-β*ΔX)) Y
# ΔMortality: increase in mortality from exposure to all PM or fire PM
# Pop: total population --> varies by cell and by year
# β = ln(RR)/ΔX (concentration response from literature)
# ΔX = PM2.5 concentration threshold above which deadly --> varies by cell and by year.
  # ΔX[lon, lat, year] = max(0, PM2.5 - threshold) --> if negative, then 0
# Y = baseline all-cause mortality --> is constant by cell but varies by year

# Parameters from Ford et al 2018 Table 1
RR_krewski <- 1.06
beta_k <- log(1.06)/10
threshold_crouse <- 1.9 # units μg/m^3
Y <- 0.0085 # USA 2000 mortality rate placeholder (must update to annual by country)

# fpm = pm_total - pm_noFire, where pm is simulated with and without fire
# Two approaches compute fpm mortality 
# 1: fpm_mortality = pm_total_mortality x (fpm/ pm_total)
# 2: fpm_mortality = pm_total_mortality - pm_noFire_mortality

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Calculate PM2000 Mortality #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Calculate ΔX = max(0, PM2.5 - threshold)
# ΔMortality = Pop × (1 - exp(-β×ΔX)) × Y

pop_pm_country <- pop_pm_country %>%
  mutate(
    # Step 1: Calculate ΔX (PM above threshold)
    delta_x_2000 = pmax(0, pm_2000 - threshold_crouse),
    
    # Step 2: Calculate mortality using Ford et al 2018 equation (5)
    # ΔMortality = Pop × (1 - exp(-β×ΔX)) × Y
    pm_2000_mort = pop_tot_2009 * (1 - exp(-beta_k * delta_x_2000)) * Y
  )

# Check results
summary(pop_pm_country$pm_2000_mort)

# View sample of results
pop_pm_country %>%
  select(lon, lat, country_name, pop_tot_2009, pm_2000, delta_x_2000, pm_2000_mort) %>%
  filter(!is.na(pm_2000_mort) & pm_2000_mort > 0) %>%
  arrange(desc(pm_2000_mort)) %>%
  head(10)

