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
pop_pm_country <- read_csv(here("output", "pop_pm_country_death.csv"))
colnames(pop_pm_country)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Baseline population weights ########################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# pop_bar_ic = average population in grid cell i over baseline period 2001-2010
# Population is held fixed across all periods and scenarios 
pop_yr_cols <- paste0("pop_tot_", 2001:2010)

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

fpm_cols <- c(
  "fpm_2000",
  "fpm_2050_45",
  "fpm_2050_85",
  "fpm_2100_45",
  "fpm_2100_85"
)

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
    .groups = "drop"
  )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Join Death Rate and Compute Baseline Death Country Level ########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Extract one death_rate per country per year
# death_rate is constant within country — take first non-NA row
country_death_rate <- pop_pm_country %>%
  filter(!is.na(country_code_iso3)) %>%
  group_by(country_code_iso3, country_name) %>%
  slice(1) %>%
  ungroup() %>%
  select(country_code_iso3, country_name,
         all_of(paste0("death_rate_", 2001:2010)))

# Join death rates to country-level pop-weighted PM
pop_wght_pm_cntry <- pop_wght_pm_cntry %>%
  left_join(country_death_rate, by = c("country_code_iso3", "country_name"))

# Compute baseline expected deaths: base_death_yr = pop_tot_yr * death_rate_yr
for (yr in 2001:2010) {
  pop_col <- paste0("pop_tot_",    yr)
  dr_col  <- paste0("death_rate_", yr)
  bd_col  <- paste0("base_death_", yr)

  pop_wght_pm_cntry <- pop_wght_pm_cntry %>%
    # .data refers to the current dataframe; [[col]] selects the column whose name is stored in the string variable col
    mutate(!!bd_col := .data[[pop_col]] * .data[[dr_col]])  # base_death_yr = pop_tot_yr * death_rate_yr
}

# Average baseline deaths across 2001-2010
pop_wght_pm_cntry <- pop_wght_pm_cntry %>%
  mutate(base_death_ave_2001_10 = rowMeans(
    across(all_of(paste0("base_death_", 2001:2010))), na.rm = TRUE
  ))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Add Global Row #############################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fpm_suffixes <- c("2000", "2050_45", "2050_85", "2100_45", "2100_85")

# Extract Somaliland death rates as global death rate (assigned earlier in death rate data)
somaliland_dr <- pop_pm_country %>%
  filter(country_code_iso3 == "-99" & country_name == "Somaliland") %>%
  slice(1) %>%
  select(all_of(paste0("death_rate_", 2001:2010)))

# Step 1: sum population, exposure, and base_death columns
global_row <- pop_wght_pm_cntry %>%
  summarise(
    country_code_iso3 = "global",
    country_name      = "global",
    pop_bar_c         = sum(pop_bar_c,  na.rm = TRUE),
    across(all_of(paste0("pop_tot_",      2001:2010)),    ~ sum(.x, na.rm = TRUE)),
    across(all_of(paste0("exposure_fpm_", fpm_suffixes)), ~ sum(.x, na.rm = TRUE)),
    across(all_of(paste0("base_death_",   2001:2010)),    ~ sum(.x, na.rm = TRUE))
  )

# Step 2: per-capita exposure = global total exposure / global pop_bar_c
for (s in fpm_suffixes) {
  global_row[[paste0("exposure_percap_fpm_", s)]] <-
    global_row[[paste0("exposure_fpm_", s)]] / global_row$pop_bar_c
}

# Step 3: assign Somaliland death rates to global row
for (yr in 2001:2010) {
  dr_col <- paste0("death_rate_", yr)
  global_row[[dr_col]] <- somaliland_dr[[dr_col]]
}

# Step 4: average baseline deaths across 2001-2010
global_row$base_death_ave_2001_10 <- mean(
  unlist(global_row[paste0("base_death_", 2001:2010)]), na.rm = TRUE
)

# Append global row
pop_wght_pm_cntry <- bind_rows(pop_wght_pm_cntry, global_row)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Export #############################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

write_csv(pop_wght_pm_cntry, here("output", "pop_wght_pm_cntry.csv"))

glimpse(pop_wght_pm_cntry)



