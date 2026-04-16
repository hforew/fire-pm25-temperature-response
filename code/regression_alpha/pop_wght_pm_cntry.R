# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## compute population weighted exposure by country ##
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

# Import combined data
# pop_pm_country <- read_csv(here("output", "pop_pm_with_countries.csv"))
# pop_pm_country <- read_csv(here("output", "pop_pm_country_death.csv"))
pop_pm_country <- read_csv(here("output", "pop_pm_country_death_park.csv"))

colnames(pop_pm_country)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Baseline population weights ########################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# pop_bar_ic = average population in grid cell i over baseline period 2001-2010
# Population is held fixed across all periods and scenarios 
pop_yr_cols <- paste0("pop_tot_", 2001:2010)

# compute ave pop 2001-2010 from the Pierce et al data 
pop_pm_country <- pop_pm_country %>%
  mutate(pop_bar = rowMeans(across(all_of(pop_yr_cols)), na.rm = TRUE))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Fire PM2.5 exposure by country and scenario ########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Total exposure formula:
#   exposure_cps = sum_i( fPM_ics * pop_bar_ic )
#
# Per-capita exposure formula:
#   exposure_cps^per-capita = sum_i( fPM_ics * pop_bar_ic ) / pop_bar_c
#
# fpm_* columns are already decadal averages, so 1/10 * sum_t collapses to the column value

# pierce data FPM 
fpm_cols <- c(
  "fpm_2000",
  "fpm_2050_45",
  "fpm_2050_85",
  "fpm_2100_45",
  "fpm_2100_85"
)

# Park et al. fPM column names: park_{model}_{decade}_fpm (18 columns: 3 models x 6 decades)
decades       <- c("1960s", "1970s", "1980s", "1990s", "2000s", "2010s")
models        <- c("classic", "jules", "ssib4")  # three land-surface model variants in Park et al.
park_fpm_cols <- as.vector(outer(models, decades, function(m, d) paste0("park_", m, "_", d, "_fpm")))

pop_wght_pm_cntry <- pop_pm_country %>%
  group_by(country_code_iso3, country_name) %>%  # one group per country
  summarise(
    # Total baseline population in country c: pop_bar_c = sum_i( pop_bar_ic )
    pop_bar_c = sum(pop_bar, na.rm = TRUE),
    # Sum grid-cell population to country level for each year 2001-2010
    across(all_of(pop_yr_cols), ~ sum(.x, na.rm = TRUE)),
    # Total fire PM2.5 exposure in country c:
    # exposure_cps = sum_i( fPM_ics * pop_bar_ic )
    # i.e. weighted sum of grid-cell fire PM by baseline population
    across(
      all_of(fpm_cols),
      ~ sum(.x * pop_bar, na.rm = TRUE),
      .names = "exposure_{.col}"
    ),
    # Per-capita fire PM2.5 exposure in country c:
    # exposure_cps^per-capita = sum_i( fPM_ics * pop_bar_ic ) / pop_bar_c
    # denominator uses only cells where fPM is non-missing
    across(
      all_of(fpm_cols),
      ~ sum(.x * pop_bar, na.rm = TRUE) /
          sum(pop_bar[!is.na(.x)], na.rm = TRUE),
      .names = "exposure_percap_{.col}"
    ),
    # Park et al. fPM exposure: all 18 columns weighted by pop_bar (same baseline as future scenarios).
    # Using pop_bar holds population geography fixed at 2001-2010, isolating PM concentration changes.
    across(
      all_of(park_fpm_cols),
      ~ sum(.x * pop_bar, na.rm = TRUE),
      .names = "exposure_{.col}"
    ),
    across(
      all_of(park_fpm_cols),
      ~ sum(.x * pop_bar, na.rm = TRUE) /
          sum(pop_bar[!is.na(.x)], na.rm = TRUE),
      .names = "exposure_percap_{.col}"
    ),
    .groups = "drop"
  )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Join Death Rate and Compute Baseline Death Country Level ########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Single baseline death rate per country: mean of annual rates over 2001-2010.
# Consistent with holding population fixed at pop_bar (2001-2010 average).
# death_rate is constant within a country across grid cells — slice(1) picks the representative row.
country_death_rate <- pop_pm_country %>%
  filter(!is.na(country_code_iso3)) %>%
  group_by(country_code_iso3, country_name) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(death_rate_base = rowMeans(
    across(all_of(paste0("death_rate_", 2001:2010))), na.rm = TRUE
  )) %>%
  select(country_code_iso3, country_name, death_rate_base)

# Join baseline death rate to country-level pop-weighted PM
pop_wght_pm_cntry <- pop_wght_pm_cntry %>%
  left_join(country_death_rate, by = c("country_code_iso3", "country_name"))

# Baseline expected deaths: base_death = pop_bar_c * death_rate_base
# Both population and death rate are held fixed at the 2001-2010 baseline
pop_wght_pm_cntry <- pop_wght_pm_cntry %>%
  mutate(base_death = pop_bar_c * death_rate_base)

# check: USA death_rate_base should be ~0.008 (WB crude death rate 2001-2010 avg)
pop_wght_pm_cntry %>%
  filter(country_code_iso3 == "USA") %>%
  select(country_name, pop_bar_c, death_rate_base, base_death)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Add Global Row #############################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# pierce data
fpm_suffixes <- c("2000", "2050_45", "2050_85", "2100_45", "2100_85")

# World average baseline death rate: mean of WLD crude death rate 2001-2010, from World Bank data
wb_death_rate <- read_csv(
  here("input", "WB_crude_death_rate", "API_SP.DYN.CDRT.IN_DS2_en_csv_v2_241.csv"),
  skip = 4
)

world_dr_base <- wb_death_rate %>%
  select(country_code_iso3 = `Country Code`, `2001`:`2010`) %>%
  filter(country_code_iso3 == "WLD") %>%
  mutate(across(`2001`:`2010`, ~ .x / 1000)) %>%       # convert per-1000 to proportion
  summarise(death_rate_base = rowMeans(across(`2001`:`2010`), na.rm = TRUE)) %>%
  pull(death_rate_base)

# Park exposure column names: exposure_park_{model}_{decade}_fpm (18 columns)
park_exp_cols <- paste0("exposure_park_", rep(models, times = length(decades)), "_",
                        rep(decades, each = length(models)), "_fpm")

# Step 1: sum population, exposure, and base_death
global_row <- pop_wght_pm_cntry %>%
  summarise(
    country_code_iso3 = "global",
    country_name      = "global",
    pop_bar_c         = sum(pop_bar_c,  na.rm = TRUE),
    across(all_of(paste0("pop_tot_",        2001:2010)),    ~ sum(.x, na.rm = TRUE)),
    across(all_of(paste0("exposure_fpm_",   fpm_suffixes)), ~ sum(.x, na.rm = TRUE)),
    base_death        = sum(base_death, na.rm = TRUE),
    # Park: sum total exposure to global level (per-capita computed in Step 2b below)
    across(all_of(park_exp_cols),                           ~ sum(.x, na.rm = TRUE))
  )

# Step 2: per-capita exposure = global total exposure / global pop_bar_c
for (s in fpm_suffixes) {
  global_row[[paste0("exposure_percap_fpm_", s)]] <-
    global_row[[paste0("exposure_fpm_", s)]] / global_row$pop_bar_c
}

# Step 2b: Park per-capita exposure = global total Park exposure / global pop_bar_c
for (exp_col in park_exp_cols) {
  percap_col <- sub("^exposure_", "exposure_percap_", exp_col)
  global_row[[percap_col]] <- global_row[[exp_col]] / global_row$pop_bar_c
}

# Step 3: assign world average baseline death rate to global row
global_row$death_rate_base <- world_dr_base

# Append global row
pop_wght_pm_cntry <- bind_rows(pop_wght_pm_cntry, global_row)

# check global baseline deaths: base_death = pop_bar_c * death_rate_base
# note: does NOT enter GIVE model — internal check only
pop_wght_pm_cntry %>%
  filter(country_code_iso3 == "global") %>%
  select(pop_bar_c, death_rate_base, base_death) %>%
  glimpse()

# compare global row vs sum of country rows — global pop_bar_c and base_death should match
bind_rows(
  pop_wght_pm_cntry %>%
    filter(country_code_iso3 == "global") %>%
    select(pop_bar_c, base_death) %>%
    mutate(source = "global_row"),
  pop_wght_pm_cntry %>%
    filter(country_code_iso3 != "global") %>%
    summarise(pop_bar_c = sum(pop_bar_c, na.rm = TRUE),
              base_death = sum(base_death, na.rm = TRUE)) %>%
    mutate(source = "country_sum")
) %>%
  select(source, pop_bar_c, base_death)

colnames(pop_wght_pm_cntry)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Export #############################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# write_csv(pop_wght_pm_cntry, here("output", "pop_wght_pm_cntry.csv"))
write_csv(pop_wght_pm_cntry, here("output", "pop_wght_pm_cntry_park.csv"))

glimpse(pop_wght_pm_cntry)



