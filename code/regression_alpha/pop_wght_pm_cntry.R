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
############ Park et al. fPM exposure by country and decade ############################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Unlike the future-scenario fPM columns (which all share pop_bar as the weight),
# Park et al. decades each have their own population column: park_{decade}_pop.
# A single across() cannot vary the weight column per group of fPM columns,
# so we loop over decades, run one summarise() per decade, then join the results.

decades <- c("1960s", "1970s", "1980s", "1990s", "2000s", "2010s")
models  <- c("classic", "jules", "ssib4")  # three land-surface model variants in Park et al.

# lapply returns a list of 6 dataframes, one per decade.
# Each dataframe has columns: country_code_iso3, country_name,
#   exposure_park_{model}_{decade}_fpm        (total pop-weighted exposure)
#   exposure_percap_park_{model}_{decade}_fpm (per-capita pop-weighted exposure)
park_exposure_list <- lapply(decades, function(d) {
  
  # Name of the decade-specific population column used as the weight
  pop_col    <- paste0("park_", d, "_pop")           # e.g. "park_1960s_pop"
  
  # The three fPM columns for this decade (one per land-surface model)
  fpm_cols_d <- paste0("park_", models, "_", d, "_fpm")
  # e.g. c("park_classic_1960s_fpm", "park_jules_1960s_fpm", "park_ssib4_1960s_fpm")
  
  pop_pm_country %>%
    group_by(country_code_iso3, country_name) %>%
    summarise(

      # --- Country-level total population for this decade ---
      # park_{decade}_pop_c = sum_i( park_{decade}_pop_i )
      # Computed here alongside exposure since the same weight column is already in scope.
      # Used downstream to compute base_death_park_{decade}.
      !!paste0("park_", d, "_pop_c") := sum(.data[[pop_col]], na.rm = TRUE),

      # --- Total pop-weighted exposure ---
      # exposure_park_{model}_{decade}_fpm = sum_i( fPM_i * pop_i )
      # .data[[pop_col]] accesses the decade-specific pop column by string name at runtime.
      # na.rm = TRUE ignores grid cells where fPM or pop is missing.
      across(
        all_of(fpm_cols_d),
        ~ sum(.x * .data[[pop_col]], na.rm = TRUE), # ~ anonymous function; .x is column (fpm) iterating over
        .names = "exposure_{.col}"           # e.g. "exposure_park_classic_1960s_fpm"
      ),

      # --- Per-capita pop-weighted exposure ---
      # exposure_percap = sum_i( fPM_i * pop_i ) / sum_i( pop_i )  [only over cells where fPM is non-missing]
      # The denominator restricts to cells where .x (fPM) is non-NA so that
      # missing fPM cells do not inflate the population base and bias the per-capita figure.
      across(
        all_of(fpm_cols_d),
        ~ sum(.x * .data[[pop_col]], na.rm = TRUE) /
          sum(.data[[pop_col]][!is.na(.x)], na.rm = TRUE),
        .names = "exposure_percap_{.col}"    # e.g. "exposure_percap_park_classic_1960s_fpm"
      ),

      .groups = "drop"
    )
})

# reduce() (purrr, loaded via tidyverse) sequentially left_joins the 6 decade dataframes
# on country_code_iso3 + country_name, producing one wide dataframe with all 36 new columns.
park_exposure <- reduce(park_exposure_list, left_join, by = c("country_code_iso3", "country_name"))

# Join the 36 Park exposure columns onto the country-level summary built above.
# Left join preserves all rows in pop_wght_pm_cntry.
pop_wght_pm_cntry <- pop_wght_pm_cntry %>%
  left_join(park_exposure, by = c("country_code_iso3", "country_name"))

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
         all_of(paste0("death_rate_", 1960:2020)))

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
############ Park et al. Baseline Deaths by Decade ######################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# base_death_park_{decade} = park_{decade}_pop_c * mean( death_rate_{yr} for yr in decade )
#
# Unlike base_death_ave_2001_10 which averages the product of varying annual pop × death rate,
# here the population is fixed within each decade (one Park pop value per decade),
# so we average the death rates across the decade years first, then multiply by decade population.

# Map each decade label to its constituent years
decade_years <- list(
  "1960s" = 1960:1969,
  "1970s" = 1970:1979,
  "1980s" = 1980:1989,
  "1990s" = 1990:1999,
  "2000s" = 2000:2009,
  "2010s" = 2010:2019
)

for (d in names(decade_years)) {
  yrs     <- decade_years[[d]]
  pop_col <- paste0("park_", d, "_pop_c")     # e.g. "park_1960s_pop_c"
  dr_cols <- paste0("death_rate_", yrs)        # e.g. "death_rate_1960"..."death_rate_1969"
  bd_col  <- paste0("base_death_park_", d)     # e.g. "base_death_park_1960s"

  # Compute row-wise mean death rate across the decade, then multiply by decade population
  pop_wght_pm_cntry <- pop_wght_pm_cntry %>%
    mutate(!!bd_col := .data[[pop_col]] * rowMeans(across(all_of(dr_cols)), na.rm = TRUE))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Add Global Row #############################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fpm_suffixes <- c("2000", "2050_45", "2050_85", "2100_45", "2100_85")

# Extract Somaliland death rates as global death rate (assigned earlier in death rate data)
somaliland_dr <- pop_pm_country %>%
  filter(country_code_iso3 == "-99" & country_name == "Somaliland") %>%
  slice(1) %>%
  select(all_of(paste0("death_rate_", 1960:2020)))

# Park exposure column names: exposure_park_{model}_{decade}_fpm (18 columns)
park_exp_cols <- paste0("exposure_park_", rep(models, times = length(decades)), "_",
                        rep(decades, each = length(models)), "_fpm")

# Step 1: sum population, exposure, and base_death columns
global_row <- pop_wght_pm_cntry %>%
  summarise(
    country_code_iso3 = "global",
    country_name      = "global",
    pop_bar_c         = sum(pop_bar_c,  na.rm = TRUE),
    across(all_of(paste0("pop_tot_",        2001:2010)),    ~ sum(.x, na.rm = TRUE)),
    across(all_of(paste0("exposure_fpm_",   fpm_suffixes)), ~ sum(.x, na.rm = TRUE)),
    across(all_of(paste0("base_death_",     2001:2010)),    ~ sum(.x, na.rm = TRUE)),
    # Park: sum decade populations to global level
    across(all_of(paste0("park_", decades, "_pop_c")),      ~ sum(.x, na.rm = TRUE)),
    # Park: sum total exposure to global level (per-capita computed in Step 2b below)
    across(all_of(park_exp_cols),                           ~ sum(.x, na.rm = TRUE)),
    # Park: sum baseline deaths to global level
    across(all_of(paste0("base_death_park_", decades)),     ~ sum(.x, na.rm = TRUE))
  )

# Step 2: per-capita exposure = global total exposure / global pop_bar_c
for (s in fpm_suffixes) {
  global_row[[paste0("exposure_percap_fpm_", s)]] <-
    global_row[[paste0("exposure_fpm_", s)]] / global_row$pop_bar_c
}

# Step 2b: Park per-capita exposure = global total Park exposure / global park_{decade}_pop_c
for (d in decades) {
  pop_col <- paste0("park_", d, "_pop_c")
  for (m in models) {
    exp_col    <- paste0("exposure_park_", m, "_", d, "_fpm")
    percap_col <- paste0("exposure_percap_park_", m, "_", d, "_fpm")
    global_row[[percap_col]] <- global_row[[exp_col]] / global_row[[pop_col]]
  }
}

# Step 3: assign Somaliland death rates to global row (Somaliand previously given global rate)
for (yr in 1960:2020) {
  dr_col <- paste0("death_rate_", yr)
  global_row[[dr_col]] <- somaliland_dr[[dr_col]]
}

# Step 4: average baseline deaths across 2001-2010
global_row$base_death_ave_2001_10 <- mean(
  unlist(global_row[paste0("base_death_", 2001:2010)]), na.rm = TRUE
)

# Step 4b: Park baseline deaths for global row
# park_{decade}_pop_c (global) * mean Somaliland death rate across decade years
for (d in names(decade_years)) {
  yrs     <- decade_years[[d]]
  pop_col <- paste0("park_", d, "_pop_c")
  dr_cols <- paste0("death_rate_", yrs)
  bd_col  <- paste0("base_death_park_", d)
  global_row[[bd_col]] <- global_row[[pop_col]] *
    mean(unlist(somaliland_dr[dr_cols]), na.rm = TRUE)
}

# Append global row
pop_wght_pm_cntry <- bind_rows(pop_wght_pm_cntry, global_row)

# compare base death from the Park and Pierce data sets (different populuations, with same death rate)
# note that this base death data is only for internal checks but DOES NOT enter GIVE model
pop_wght_pm_cntry %>%
  filter(country_code_iso3 == "global") %>%
  select(base_death_ave_2001_10, starts_with("base_death_park_")) %>%
  glimpse()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Export #############################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# write_csv(pop_wght_pm_cntry, here("output", "pop_wght_pm_cntry.csv"))
write_csv(pop_wght_pm_cntry, here("output", "pop_wght_pm_cntry_park.csv"))

glimpse(pop_wght_pm_cntry)



