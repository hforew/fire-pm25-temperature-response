# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## POPULATION-WEIGHTED FIRE PM2.5 EXPOSURE BY COUNTRY ##
##
## Goal: Aggregate grid-cell fire PM2.5 to country-level population-weighted per-capita
##   exposure for Pierce, Park, and Zhao scenarios. Joins baseline death rates, appends
##   a global row, and exports for use in country-level FE regressions.
##
## Inputs:
##   output/pop_pm_country_death_final.csv         (grid-cell fPM with population and
##                                                  death rates; Pierce, Park, Zhao cols)
##   input/WB_crude_death_rate/API_SP.DYN.CDRT.IN_DS2_en_csv_v2_241.csv
##                                                 (World Bank crude death rates; used
##                                                  to assign global-row death_rate_base)
##
## Outputs:
##   output/pop_wght_pm_cntry_final.csv            (one row per country + global; columns:
##                                                  pop_bar_c, death_rate_base, base_death,
##                                                  exposure_percap_{scenario} for all
##                                                  Pierce, Park, and Zhao scenarios)
##
## Execution order:
##   files run before: mortality/death_rate_processing.R   --> writes pop_pm_country_death_final.csv
##   files run after:  FE/fpm_gmt_regression_FE_all.R
##                     FE_climate_attributable/fpm_gmt_regression_FE_all_cli.R
##                     FE_climate_attributable/fpm_gmt_regression_park_cli.R
##                     regression_FE_stats_all.R
##                     single_study/fpm_gmt_regression_pierce.R
##                     single_study/fpm_gmt_regression_park.R
##
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
# pop_pm_country <- read_csv(here("output", "pop_pm_country_death_park.csv"))
pop_pm_country <- read_csv(here("output", "pop_pm_country_death_final.csv"))

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

# Notation: i = grid cell, c = country, t = time period, s = scenario, m = fire model
#   PM_itsm   = fire PM2.5 in cell i; already a period avg
#   pop_bar_ic = baseline pop in cell i of country c; avg 2001-2010 from Pierce et al.
#   pop_bar_c  = sum_i( pop_bar_ic )
#
# Total exposure:   PM_ctsm    = sum_i( PM_itsm * pop_bar_ic )
# Per-capita:       PM_bar_ctsm = PM_ctsm / pop_bar_c
#
# Park et al.:   t in {1960s, 70s, 80s, 90s, 2000s, 2010s}, s in {historical}, m in {CLASSIC, JULES, SSiB4}
# Pierce et al.: t in {2000s, 2040s, 2090s}, s in {historical, RCP4.5, RCP8.5}, m in {CESM}
# Zhao et al.:   t in {2095-99}, s in {SSP2-4.5, SSP5-8.5}, m in {ML}

# pierce data FPM + Zhao et al. 2025
fpm_cols <- c(
  "fpm_2000",
  "fpm_2050_45",
  "fpm_2050_85",
  "fpm_2100_45",
  "fpm_2100_85",
  "fpm_2095_SSP245_Zhao",
  "fpm_2095_SSP585_Zhao"
)

# Park et al. PM column names: park_{model}_{decade}_fpm (18 columns: 3 models x 6 decades)
decades       <- c("1960s", "1970s", "1980s", "1990s", "2000s", "2010s")
models        <- c("classic", "jules", "ssib4")  # three land-surface model variants in Park et al.
park_fpm_cols     <- as.vector(outer(models, decades, function(m, d) paste0("park_", m, "_", d, "_fpm")))
park_fpm_cli_cols <- as.vector(outer(models, decades, function(m, d) paste0("park_", m, "_", d, "_fpm_cli")))

pop_wght_pm_cntry <- pop_pm_country %>%
  group_by(country_code_iso3, country_name) %>%  # one group per country
  summarise(
    # Total baseline population in country c: pop_bar_c = sum_i( pop_bar_ic )
    pop_bar_c = sum(pop_bar, na.rm = TRUE),
    # Sum grid-cell population to country level for each year 2001-2010
    across(all_of(pop_yr_cols), ~ sum(.x, na.rm = TRUE)),
    # PM_ctsm = sum_i( PM_itsm * pop_bar_ic )  -- Pierce + Zhao
    across(
      all_of(fpm_cols),
      ~ sum(.x * pop_bar, na.rm = TRUE),
      .names = "exposure_{.col}"
    ),
    # PM_bar_ctsm = PM_ctsm / pop_bar_c  -- Pierce + Zhao
    across(
      all_of(fpm_cols),
      ~ sum(.x * pop_bar, na.rm = TRUE) / pop_bar_c,
      .names = "exposure_percap_{.col}"
    ),
    # Park et al. PM exposure: all 18 columns weighted by pop_bar (same baseline as future scenarios).
    # Using pop_bar holds population geography fixed at 2001-2010, isolating PM concentration changes.
    across(
      all_of(park_fpm_cols),
      ~ sum(.x * pop_bar, na.rm = TRUE),
      .names = "exposure_{.col}"
    ),
    across(
      all_of(park_fpm_cols),
      ~ sum(.x * pop_bar, na.rm = TRUE) / pop_bar_c,
      .names = "exposure_percap_{.col}"
    ),
    # Park et al. climate-attributable fire PM exposure (factual minus counterfactual)
    across(
      all_of(park_fpm_cli_cols),
      ~ sum(.x * pop_bar, na.rm = TRUE),
      .names = "exposure_{.col}"
    ),
    across(
      all_of(park_fpm_cli_cols),
      ~ sum(.x * pop_bar, na.rm = TRUE) / pop_bar_c,
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

# pierce data + Zhao et al. 2025
fpm_suffixes <- c("2000", "2050_45", "2050_85", "2100_45", "2100_85", "2095_SSP245_Zhao", "2095_SSP585_Zhao")

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
park_exp_cols     <- paste0("exposure_park_", rep(models, times = length(decades)), "_",
                            rep(decades, each = length(models)), "_fpm")
park_exp_cli_cols <- paste0("exposure_park_", rep(models, times = length(decades)), "_",
                            rep(decades, each = length(models)), "_fpm_cli")

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
    across(all_of(park_exp_cols),                           ~ sum(.x, na.rm = TRUE)),
    across(all_of(park_exp_cli_cols),                       ~ sum(.x, na.rm = TRUE))
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

# Step 2c: Park cli per-capita exposure = global total Park cli exposure / global pop_bar_c
for (exp_col in park_exp_cli_cols) {
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
############ Descriptive Statistics: Per-Capita PM Exposure (excl. JULES) #############
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# select all per-capita exposure cols, drop JULES (excluded from analysis)
desc_cols <- colnames(pop_wght_pm_cntry) %>%
  .[grepl("^exposure_percap_", .)] %>%
  .[!grepl("jules", .)]

# pivot to long so each row is one (country, column, value) triple,
# then compute summary stats across countries for each exposure column
pop_wght_pm_cntry %>%
  select(all_of(desc_cols)) %>%
  pivot_longer(everything(), names_to = "column", values_to = "value") %>%
  group_by(column) %>%
  summarise(
    n      = sum(!is.na(value)),          # non-missing country count
    mean   = mean(value,             na.rm = TRUE),
    sd     = sd(value,               na.rm = TRUE),
    min    = min(value,              na.rm = TRUE),
    p25    = quantile(value, 0.25,   na.rm = TRUE),
    median = median(value,           na.rm = TRUE),
    p75    = quantile(value, 0.75,   na.rm = TRUE),
    max    = max(value,              na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print(n = Inf)


# subset to Park-only per-capita cols (non-jules) to investigate high max values
park_percap_cols <- desc_cols[grepl("park", desc_cols)]

# for each Park column, return the 3 highest-exposure countries —
# identifies which countries drive the distribution's upper tail
park_percap_cols %>%
  lapply(function(col) {
    pop_wght_pm_cntry %>%
      filter(!is.na(.data[[col]]), country_code_iso3 != "global") %>% # exclude synthetic global row
      slice_max(order_by = .data[[col]], n = 3) %>%
      select(country_code_iso3, country_name, value = all_of(col)) %>%
      mutate(column = col)
  }) %>%
  bind_rows() %>%
  select(column, country_code_iso3, country_name, value) %>%
  print(n = Inf)

# sanity check on CAF (Central African Republic), which has the highest per-capita
# PM across all Park decades: verify high values reflect real grid-cell PM,
# not a data/join artifact. mean_fpm ~42, pop ~3.5M → per-capita ~48 is expected.
pop_pm_country %>%
  filter(country_code_iso3 == "CAF") %>%
  summarise(
    n_cells        = n(),                                      # number of grid cells in CAF
    mean_fpm       = mean(park_classic_2010s_fpm, na.rm = TRUE), # avg fire PM across cells
    max_fpm        = max(park_classic_2010s_fpm,  na.rm = TRUE), # upper tail of cell distribution
    pop_bar_total  = sum(pop_bar,                 na.rm = TRUE)  # total baseline population
  )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ USA Sanity Check #############
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# choose USA only per capita exposure results
usa_percap <- pop_wght_pm_cntry %>%
  filter(country_code_iso3 == "USA") %>%
  select(country_code_iso3, pop_bar_c, starts_with("exposure_percap_"), -contains("jules"))

usa_percap

# print ssib4
usa_percap %>%
  select(country_code_iso3, pop_bar_c, contains("ssib4"))

# grid cell level desc stats for ssib4 for USA
pop_pm_country %>%
  filter(country_code_iso3 == "USA") %>%
  select(contains("ssib4")) %>%
  pivot_longer(everything(), names_to = "column", values_to = "value") %>%
  group_by(column) %>%
  summarise(
    n      = sum(!is.na(value)),
    mean   = mean(value,           na.rm = TRUE),
    sd     = sd(value,             na.rm = TRUE),
    min    = min(value,            na.rm = TRUE),
    p25    = quantile(value, 0.25, na.rm = TRUE),
    median = median(value,         na.rm = TRUE),
    p75    = quantile(value, 0.75, na.rm = TRUE),
    max    = max(value,            na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print(n = Inf)

# USA SSiB4 per-capita exposure: identify highest and lowest decade,
# then print grid-cell desc stats for those decades to substantiate
usa_ssib4_percap <- pop_wght_pm_cntry %>%
  filter(country_code_iso3 == "USA") %>%
  select(matches("^exposure_percap_park_ssib4.*_fpm$")) %>%
  pivot_longer(everything(), names_to = "column", values_to = "percap") %>%
  mutate(decade = gsub(".*ssib4_(\\w+)_fpm", "\\1", column))  # extract decade label

# highest and lowest per-capita SSiB4 decade for USA
usa_ssib4_max <- usa_ssib4_percap %>% slice_max(percap, n = 1)
usa_ssib4_min <- usa_ssib4_percap %>% slice_min(percap, n = 1)

# grid-cell desc stats for highest SSiB4 decade — substantiates the per-capita value
usa_ssib4_max_stats <- pop_pm_country %>%
  filter(country_code_iso3 == "USA") %>%
  summarise(
    n_cells = n(),
    mean    = mean(.data[[paste0("park_ssib4_", usa_ssib4_max$decade, "_fpm")]], na.rm = TRUE),
    sd      = sd(.data[[paste0("park_ssib4_", usa_ssib4_max$decade, "_fpm")]],   na.rm = TRUE),
    min     = min(.data[[paste0("park_ssib4_", usa_ssib4_max$decade, "_fpm")]],  na.rm = TRUE),
    median  = median(.data[[paste0("park_ssib4_", usa_ssib4_max$decade, "_fpm")]], na.rm = TRUE),
    max     = max(.data[[paste0("park_ssib4_", usa_ssib4_max$decade, "_fpm")]],  na.rm = TRUE)
  )

# grid-cell desc stats for lowest SSiB4 decade — substantiates the per-capita value
usa_ssib4_min_stats <- pop_pm_country %>%
  filter(country_code_iso3 == "USA") %>%
  summarise(
    n_cells = n(),
    mean    = mean(.data[[paste0("park_ssib4_", usa_ssib4_min$decade, "_fpm")]], na.rm = TRUE),
    sd      = sd(.data[[paste0("park_ssib4_", usa_ssib4_min$decade, "_fpm")]],   na.rm = TRUE),
    min     = min(.data[[paste0("park_ssib4_", usa_ssib4_min$decade, "_fpm")]],  na.rm = TRUE),
    median  = median(.data[[paste0("park_ssib4_", usa_ssib4_min$decade, "_fpm")]], na.rm = TRUE),
    max     = max(.data[[paste0("park_ssib4_", usa_ssib4_min$decade, "_fpm")]],  na.rm = TRUE)
  )

## comparison check 

cat("USA SSiB4 highest per-capita exposure:",
    usa_ssib4_max$decade, "—", round(usa_ssib4_max$percap, 3), "µg/m³\n")

cat("\nGrid-cell PM stats for USA SSiB4", usa_ssib4_max$decade, "(highest):\n")
print(usa_ssib4_max_stats)


cat("USA SSiB4 lowest  per-capita exposure:",
    usa_ssib4_min$decade, "—", round(usa_ssib4_min$percap, 3), "µg/m³\n")

cat("\nGrid-cell PM stats for USA SSiB4", usa_ssib4_min$decade, "(lowest):\n")
print(usa_ssib4_min_stats)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Export #############################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# write_csv(pop_wght_pm_cntry, here("output", "pop_wght_pm_cntry.csv"))
# write_csv(pop_wght_pm_cntry, here("output", "pop_wght_pm_cntry_park.csv"))
write_csv(pop_wght_pm_cntry, here("output", "pop_wght_pm_cntry_final.csv"))

glimpse(pop_wght_pm_cntry)


### THE END 
