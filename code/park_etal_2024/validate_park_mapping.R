# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## validate_park_mapping.R
## Diagnostic: compare Park et al. baseline deaths vs Pierce baseline deaths
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls())

library(here)
library(tidyverse)

# Run upstream script to build pop_wght_pm_cntry
source(here("code", "regression_alpha", "pop_wght_pm_cntry.R"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Diagnostic: Compare base_death_park_2000s vs base_death_ave_2001_10 ######
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 1. Country-level scatter: are the differences uniform or driven by specific countries?
death_compare <- pop_wght_pm_cntry %>%
  filter(country_code_iso3 != "global") %>%
  select(country_code_iso3, country_name,
         base_death_ave_2001_10, base_death_park_2000s,
         pop_bar_c,                        # Pierce baseline pop
         park_2000s_pop_c) %>%             # Park 2000s pop
  mutate(
    ratio_death = base_death_park_2000s / base_death_ave_2001_10,
    ratio_pop   = park_2000s_pop_c / pop_bar_c
  )

# Summary: how different are the populations driving the discrepancy?
death_compare %>%
  summarise(
    global_pierce_pop   = sum(pop_bar_c,        na.rm = TRUE),
    global_park_pop     = sum(park_2000s_pop_c, na.rm = TRUE),
    global_pierce_death = sum(base_death_ave_2001_10, na.rm = TRUE),
    global_park_death   = sum(base_death_park_2000s,  na.rm = TRUE),
    pop_ratio           = global_park_pop / global_pierce_pop,
    death_ratio         = global_park_death / global_pierce_death
  )

# 2. Countries where the ratio is most extreme (pop coverage differences)
death_compare %>%
  filter(!is.na(ratio_pop)) %>%
  arrange(desc(ratio_pop)) %>%
  select(country_name, pop_bar_c, park_2000s_pop_c, ratio_pop,
         base_death_ave_2001_10, base_death_park_2000s, ratio_death) %>%
  print(n = 20)

# 3. Countries present in Park but missing in Pierce (or vice versa)
death_compare %>%
  filter(is.na(pop_bar_c) | pop_bar_c == 0 | is.na(park_2000s_pop_c) | park_2000s_pop_c == 0) %>%
  select(country_name, pop_bar_c, park_2000s_pop_c)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Check: does Somaliland have WLD death rate? ##############################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# What death rate is actually stored for Somaliland in the source data?
somaliland_dr_check <- pop_pm_country %>%
  filter(country_code_iso3 == "-99" & country_name == "Somaliland") %>%
  slice(1) %>%
  select(country_code_iso3, country_name,
         death_rate_2000, death_rate_2001, death_rate_2005, death_rate_2009)

cat("Somaliland death rates in pop_pm_country:\n")
print(somaliland_dr_check)

# Expected WLD rate ~0.0085 for 2000s.
# If Somaliland shows ~0.016, the WLD assignment in death_rate_processing.R did not apply.

# Cross-check: global row implied death rate
global_check <- pop_wght_pm_cntry %>%
  filter(country_code_iso3 == "global") %>%
  mutate(
    implied_dr_park_2000s = base_death_park_2000s / park_2000s_pop_c,
    country_sum_park_2000s = sum(
      pop_wght_pm_cntry$base_death_park_2000s[pop_wght_pm_cntry$country_code_iso3 != "global"],
      na.rm = TRUE
    )
  ) %>%
  select(base_death_park_2000s, park_2000s_pop_c,
         implied_dr_park_2000s, country_sum_park_2000s)

cat("\nGlobal row implied death rate vs country sum:\n")
print(global_check)
