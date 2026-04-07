# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## COMPUTE FIRE PM MORTALITY USING POPULATION-WEIGHTED EXPOSURE ##
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# METHODOLOGICAL DIFFERENCE from pm_mortality_orellano_cell.R:
#
# Notation (from writeup_GIVE_fPM_damage):
#   pop_ic  = average population in grid cell i of country c over baseline 2001-2010
#             held fixed across all periods and scenarios
#   pop_c   = total baseline population in country c = Σ pop_ic
#   Y_c     = baseline all-cause mortality rate, indexed by country c
#   β       = concentration-response coefficient = ln(RR) / ΔX
#   Δfpm_c  = pop-weighted per-capita fire PM2.5 exposure for country c
#
#   pm_mortality_orellano_cell.R (cell-level):
#     Step 1: compute total PM mortality per grid cell:
#               pop_tot_ic(yr) * (1 - exp(-β * pm_ic)) * Y_c(yr)
#     Step 2: scale by fire fraction per cell: mort * (fpm_ic / pm_ic)
#     Step 3: sum cell-level fire PM mortality to country
#     Step 4: average each year in period for 10-year annual ave
#   note: pop_tot_ic(yr) is year-varying, not held fixed at baseline
#
#   This script (country-level, pop-weighted):
#     Step 1: use pop-weighted per-capita fire PM per country (from pop_wght_pm_cntry.R)
#     Step 2: apply mortality equation directly to country-level fire PM:
#               (1 - exp(-β * Δfpm_c)) * (pop_c * Y_c)
#
#   KEY DIFFERENCES:
#     - Exposure: uses country pop-weighted Δfpm_c, not per-cell fire PM
#     - No total PM step: equation applied directly to fire PM, not via pm*(fpm/pm) ratio
#     - Because Δfpm_c is a fixed decadal average, (1-exp(-β*Δfpm_c)) is constant across
#       years and factors out of the decade average. The decade-average mortality reduces to:
#         (1 - exp(-β * Δfpm_c)) * (1/10) * Σ_t [ pop_ct(yr) * Y_c(yr) ]
#       This is equivalent to applying the mortality equation to the decade-average of
#       pop_ct * Y_c — i.e. mean(pop_ct * Y_ct), not mean(pop_ct) * mean(Y_ct).
#       pm_mortality_orellano_cell.R uses the same mean(pop_ct * Y_ct) approach.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
threshold_orellano <- 0               # no threshold: Orellano state no evidence for a lower safe limit

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Import #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Primary input: output from pop_wght_pm_cntry.R
# One row per country. Contains:
#   pop_bar_c             — total country baseline population: pop_c = Σ pop_ic
#   pop_tot_2001..2010    — total country population for each year: pop_ct(yr) = Σ pop_tot_ic(yr)
#   exposure_fpm_*        — total fire PM exposure: Σ(fPM_ic * pop_ic)
#   exposure_percap_fpm_* — per-capita fire PM: Δfpm_c = Σ(fPM_ic * pop_ic) / pop_c
# DIFFERENCE: pm_mortality_orellano_cell.R reads grid-cell data; here we read country aggregates
pop_wght <- read_csv(here("output", "pop_wght_pm_cntry.csv"))

# Cell-level data: used only to extract one Y_c per country per year
# Y_c is constant within a country across all grid cells, so one row per country suffices
# DIFFERENCE: pm_mortality_orellano_cell.R uses pop and Y_c at cell level throughout;
#             here Y_c is extracted at country level and multiplied by country-level pop_ct(yr)
pop_pm_country <- read_csv(here("output", "pop_pm_country_death.csv"))
colnames(pop_pm_country)  # inspect available columns

# GMT change scalars — same as pm_mortality_orellano_cell.R, included for context only
# GMT does not enter the mortality formula; used only when reporting results
gmt_chg <- read_csv(here("output", "gmt_periods_pi.csv"))
print(gmt_chg)

gmt_baseline  <- mean(c(gmt_chg$mean_gmt_45[gmt_chg$period == "2006-2010"],
                        gmt_chg$mean_gmt_85[gmt_chg$period == "2006-2010"]))
gmt_2040s_45  <- gmt_chg$mean_gmt_45[gmt_chg$period == "2041-2050"]
gmt_2040s_85  <- gmt_chg$mean_gmt_85[gmt_chg$period == "2041-2050"]
gmt_2090s_45  <- gmt_chg$mean_gmt_45[gmt_chg$period == "2091-2100"]
gmt_2090s_85  <- gmt_chg$mean_gmt_85[gmt_chg$period == "2091-2100"]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Country-Level Expected Deaths ######################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# In pm_mortality_orellano_cell.R, the mortality formula per cell is:
#   Step 1: ΔMort_ic(yr) = pop_tot_ic(yr) * (1 - exp(-β * pm_ic)) * Y_c(yr)
#   Step 2: fpm_mort_ic(yr) = ΔMort_ic(yr) * (fpm_ic / pm_ic)
# where pop_tot_ic(yr) is year-varying population (not held fixed at baseline)
#
# Y_c = baseline all-cause mortality rate, constant within country c
#
# In this script, the country-level formula is:
#   ΔMort_c(yr) = (1 - exp(-β * Δfpm_c)) * (pop_ct(yr) * Y_c(yr))
# where pop_ct(yr) * Y_c(yr) = total expected deaths in country c in year yr
#
# Y_c is constant within country c. We compute exp_deaths_c(yr) = pop_ct(yr) * Y_c(yr)
# directly at country level using pop_ct(yr) from pop_wght_pm_cntry.csv and Y_c(yr)
# extracted from one row per country in pop_pm_country_death.csv.
#
# DIFFERENCE: pm_mortality_orellano_cell.R computes pop_ic * Y_c at cell level then sums;
#             here we multiply at country level directly — no cell-level aggregation needed.

# Extract one Y_c per country per year from pop_pm_country_death.csv
# Y_c is constant within a country, so take the first non-NA row per country
# Result: one row per country, columns death_rate_2001 ... death_rate_2010
country_death_rate <- pop_pm_country %>%
  filter(!is.na(country_code_iso3) & country_code_iso3 != "-99") %>%  # drop unassigned ocean cells
  group_by(country_code_iso3, country_name) %>%
  slice(1) %>%                                                          # one row per country — Y_c is identical across cells
  ungroup() %>%
  select(country_code_iso3, country_name,
         all_of(paste0("death_rate_", 2001:2010)))                      # keep only Y_c columns

# Compute country-level expected deaths: exp_deaths_c(yr) = pop_ct(yr) * Y_c(yr)
# pop_ct(yr) comes from pop_wght_pm_cntry.csv; Y_c(yr) from country_death_rate above
# DIFFERENCE: pm_mortality_orellano_cell.R computes pop_ic * Y_c at cell level then sums;
#             here we multiply directly at country level — equivalent because Y_c is constant within country
country_exp_deaths <- pop_wght %>%
  left_join(country_death_rate, by = c("country_code_iso3", "country_name"))

for (yr in 2001:2010) {
  pop_col <- paste0("pop_tot_",    yr)   # country-level population for this year
  dr_col  <- paste0("death_rate_", yr)   # country-level Y_c for this year
  ed_col  <- paste0("exp_deaths_", yr)   # output column: pop_ct(yr) * Y_c(yr)

  country_exp_deaths <- country_exp_deaths %>%
    mutate(!!ed_col := .data[[pop_col]] * .data[[dr_col]])
}

country_exp_deaths <- country_exp_deaths %>%
  select(country_code_iso3, country_name,
         all_of(paste0("exp_deaths_", 2001:2010)))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Join Exposure and Expected Deaths ####################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Join country pop-weighted PM (from pop_wght_pm_cntry.R) with country expected deaths
# exp_deaths already computed at country level; join adds exposure columns to mortality inputs
pop_wght_cntry_mort <- pop_wght %>%
  left_join(
    country_exp_deaths,
    by = c("country_code_iso3", "country_name")
  )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ ΔX: Exposure Above Threshold ######################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Same threshold logic as pm_mortality_orellano_cell.R: ΔX = max(0, PM - threshold)
# threshold = 0, so ΔX = max(0, PM) = PM (all positive concentrations count)
# DIFFERENCE: ΔX here is the country pop-weighted per-capita fire PM (exposure_percap_fpm_*)
#             In pm_mortality_orellano_cell.R, ΔX is per-cell total PM (pm_*) — not fire PM

pop_wght_cntry_mort <- pop_wght_cntry_mort %>%
  mutate(
    delta_x_2000    = pmax(0, exposure_percap_fpm_2000    - threshold_orellano), # 2000s baseline fire PM
    delta_x_2050_45 = pmax(0, exposure_percap_fpm_2050_45 - threshold_orellano), # 2040s RCP 4.5 fire PM
    delta_x_2050_85 = pmax(0, exposure_percap_fpm_2050_85 - threshold_orellano), # 2040s RCP 8.5 fire PM
    delta_x_2100_45 = pmax(0, exposure_percap_fpm_2100_45 - threshold_orellano), # 2090s RCP 4.5 fire PM
    delta_x_2100_85 = pmax(0, exposure_percap_fpm_2100_85 - threshold_orellano)  # 2090s RCP 8.5 fire PM
  )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Compute Mortality — Year 2001 (Baseline Year) ############################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Mortality formula applied at country level:
#   ΔMort_c(yr, scenario) = (1 - exp(-β * Δfpm_c(scenario))) * (pop_c * Y_c(yr))
#
# pop_ct(yr) and Y_c(yr) are year-varying. exp_deaths_yr = pop_ct(yr) * Y_c(yr).
# 2001 exp_deaths are used for both baseline and future PM scenarios —
# matching pm_mortality_orellano_cell.R: 2001 pop/Y_c paired with 2041 and 2091 PM.
#
# DIFFERENCE: In pm_mortality_orellano_cell.R this first computes total PM mortality, then
#             scales by (fpm/pm). Here mortality is computed directly from fire PM exposure.

pop_wght_cntry_mort <- pop_wght_cntry_mort %>%
  mutate(

    # 2001 baseline: fire PM from 2000s decade, expected deaths from 2001
    fpm_2001_mort    = (1 - exp(-beta_k * delta_x_2000))    * exp_deaths_2001,

    # 2041 RCP 4.5: fire PM from 2040s decade under RCP4.5, expected deaths from 2001
    fpm_2041_45_mort = (1 - exp(-beta_k * delta_x_2050_45)) * exp_deaths_2001,

    # 2041 RCP 8.5: fire PM from 2040s decade under RCP8.5, expected deaths from 2001
    fpm_2041_85_mort = (1 - exp(-beta_k * delta_x_2050_85)) * exp_deaths_2001,

    # 2091 RCP 4.5: fire PM from 2090s decade under RCP4.5, expected deaths from 2001
    fpm_2091_45_mort = (1 - exp(-beta_k * delta_x_2100_45)) * exp_deaths_2001,

    # 2091 RCP 8.5: fire PM from 2090s decade under RCP8.5, expected deaths from 2001
    fpm_2091_85_mort = (1 - exp(-beta_k * delta_x_2100_85)) * exp_deaths_2001
  )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Loop Over Years 2002-2010 ############################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Mirror the loop in pm_mortality_orellano_cell.R: repeat mortality calculation for each
# base year 2002-2010, each time using that year's expected deaths with the same PM exposures.
# Population and death rate from year yr map to future years yr+40 (2040s) and yr+90 (2090s).
# Result: 45 new columns (9 years × 5 scenarios), e.g. fpm_2002_mort, fpm_2042_45_mort, ...

for (yr in 2002:2010) {

  yr_50  <- yr + 40   # future year in 2040s decade: 2042, 2043, ..., 2050
  yr_100 <- yr + 90   # future year in 2090s decade: 2092, 2093, ..., 2100
  ed_col <- paste0("exp_deaths_", yr)  # country expected deaths column for this year

  pop_wght_cntry_mort <- pop_wght_cntry_mort %>%
    mutate(

      # 200x baseline: fire PM from 2000s decade, expected deaths from year yr
      !!paste0("fpm_", yr, "_mort")        := (1 - exp(-beta_k * delta_x_2000))    * .data[[ed_col]],

      # 204x RCP 4.5: fire PM from 2040s RCP4.5, expected deaths from year yr
      !!paste0("fpm_", yr_50,  "_45_mort") := (1 - exp(-beta_k * delta_x_2050_45)) * .data[[ed_col]],

      # 204x RCP 8.5: fire PM from 2040s RCP8.5, expected deaths from year yr
      !!paste0("fpm_", yr_50,  "_85_mort") := (1 - exp(-beta_k * delta_x_2050_85)) * .data[[ed_col]],

      # 209x RCP 4.5: fire PM from 2090s RCP4.5, expected deaths from year yr
      !!paste0("fpm_", yr_100, "_45_mort") := (1 - exp(-beta_k * delta_x_2100_45)) * .data[[ed_col]],

      # 209x RCP 8.5: fire PM from 2090s RCP8.5, expected deaths from year yr
      !!paste0("fpm_", yr_100, "_85_mort") := (1 - exp(-beta_k * delta_x_2100_85)) * .data[[ed_col]]
    )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Decade Averages #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Average mortality across the 10 years within each decade, matching pm_mortality_orellano_cell.R
# Each country gets a single mortality value per scenario representing the decade mean

pop_wght_cntry_mort <- pop_wght_cntry_mort %>%
  mutate(

    # 2000s baseline: average of fpm_2001_mort through fpm_2010_mort
    fpm_2000s_mort    = rowMeans(pick(all_of(paste0("fpm_", 2001:2010, "_mort"))),       na.rm = TRUE),

    # 2040s RCP 4.5: average of fpm_2041_45_mort through fpm_2050_45_mort
    fpm_2040s_45_mort = rowMeans(pick(all_of(paste0("fpm_", 2041:2050, "_45_mort"))), na.rm = TRUE),

    # 2040s RCP 8.5: average of fpm_2041_85_mort through fpm_2050_85_mort
    fpm_2040s_85_mort = rowMeans(pick(all_of(paste0("fpm_", 2041:2050, "_85_mort"))), na.rm = TRUE),

    # 2090s RCP 4.5: average of fpm_2091_45_mort through fpm_2100_45_mort
    fpm_2090s_45_mort = rowMeans(pick(all_of(paste0("fpm_", 2091:2100, "_45_mort"))), na.rm = TRUE),

    # 2090s RCP 8.5: average of fpm_2091_85_mort through fpm_2100_85_mort
    fpm_2090s_85_mort = rowMeans(pick(all_of(paste0("fpm_", 2091:2100, "_85_mort"))), na.rm = TRUE)
  )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ USA Analysis ########################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Filter to USA for direct comparison with pm_mortality_orellano_cell.R results and literature
# DIFFERENCE: pm_mortality_orellano_cell.R filters grid cells to USA; here we filter country rows
usa_mort <- pop_wght_cntry_mort %>%
  filter(country_code_iso3 == "USA")

# Decade-average fire PM mortality for USA across all 5 scenarios
usa_fpm_mort_avg <- usa_mort %>%
  select(fpm_2000s_mort, fpm_2040s_45_mort, fpm_2040s_85_mort,
         fpm_2090s_45_mort, fpm_2090s_85_mort) %>%
  pivot_longer(
    cols = everything(),
    names_to  = "scenario",    # column names become scenario labels
    values_to = "fpm_deaths"   # mortality values
  ) %>%
  mutate(scenario = str_remove(scenario, "^fpm_")) %>%   # strip "fpm_" prefix: "fpm_2040s_45_mort" -> "2040s_45_mort"
  mutate(scenario = str_remove(scenario, "_mort$")) %>%  # strip "_mort" suffix: "2040s_45_mort" -> "2040s_45"
  separate(scenario, into = c("year", "rcp"), sep = "_", fill = "right") %>% # split "2040s_45" into year="2040s", rcp="45"
  select(year, rcp, fpm_deaths)

print("Fire PM Mortality Decade Averages — USA (pop-weighted PM method)")
print(usa_fpm_mort_avg)

# Literature comparison (same benchmarks as pm_mortality_orellano_cell.R)
cat("\nLiterature Comparison (2000s baseline):\n")
cat("Ford et al 2018: 17,000 deaths\n")
cat("Our estimate (pop-weighted PM, 2000s decade avg):",
    round(usa_fpm_mort_avg %>% filter(year == "2000s") %>% pull(fpm_deaths)), "\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Global Analysis #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Aggregate country-level mortality to global totals
# Same structure as global_totals_all in pm_mortality_orellano_cell.R

global_fpm_mort <- pop_wght_cntry_mort %>%
  summarise(
    fpm_2000s    = sum(fpm_2000s_mort,    na.rm = TRUE),  # global total, 2000s baseline
    fpm_2040s_45 = sum(fpm_2040s_45_mort, na.rm = TRUE),  # global total, 2040s RCP 4.5
    fpm_2040s_85 = sum(fpm_2040s_85_mort, na.rm = TRUE),  # global total, 2040s RCP 8.5
    fpm_2090s_45 = sum(fpm_2090s_45_mort, na.rm = TRUE),  # global total, 2090s RCP 4.5
    fpm_2090s_85 = sum(fpm_2090s_85_mort, na.rm = TRUE)   # global total, 2090s RCP 8.5
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to  = "scenario",
    values_to = "fpm_deaths"
  ) %>%
  mutate(scenario = str_remove(scenario, "^fpm_")) %>%
  separate(scenario, into = c("year", "rcp"), sep = "_", fill = "right") %>%
  mutate(
    gmt_chg = c(gmt_baseline, gmt_2040s_45, gmt_2040s_85, gmt_2090s_45, gmt_2090s_85),  # GMT for context
    all_death_2009 = 54100000,                                                            # total global deaths 2009 (same as pm_mortality_orellano_cell.R)
    death_pct_chg  = fpm_deaths / all_death_2009                                          # fire PM as share of all global deaths
  ) %>%
  select(year, rcp, fpm_deaths, gmt_chg, death_pct_chg)

print("\nGlobal Fire PM Mortality — Pop-Weighted PM Method")
print(global_fpm_mort)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Country Rankings ####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Top countries by fire PM mortality (2000s baseline), same format as pm_mortality_orellano_cell.R
pop_wght_cntry_mort_summary <- pop_wght_cntry_mort %>%
  select(country_code_iso3, country_name,
         fpm_2000s_mort, fpm_2040s_45_mort, fpm_2040s_85_mort,
         fpm_2090s_45_mort, fpm_2090s_85_mort) %>%
  arrange(desc(fpm_2000s_mort))  # rank by 2000s baseline mortality

print("\nTop 20 Countries by Fire PM Mortality (2000s baseline, pop-weighted PM method)")
print(pop_wght_cntry_mort_summary %>% head(20))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Export ##############################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Export country-level decade-average mortality
write_csv(
  pop_wght_cntry_mort_summary,
  here("output", "fpm_mort_pop_wght_summ.csv")
)

# Export country-level detailed data
write_csv(
  pop_wght_cntry_mort,
  here("output", "pop_wght_cntry_mort.csv")
)

glimpse(pop_wght_cntry_mort_summary)
