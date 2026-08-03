# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## PROCESS DEATH RATE DATA: Join World Bank crude death rates onto the pop/PM grid
##
## Goal: Clean World Bank crude death rate data (1960-2020, per capita) and join it
##   onto the pop/PM grid-cell table by country, filling grid cells whose country has
##   no WB match with the world-average death rate for that year.
##
## Inputs:
##   output/pm_joined/pop_pm_with_countries_final.csv           (grid-cell pop + PM2.5
##                                         + country assignment; Pierce, Park, Zhao data)
##   input/WB_crude_death_rate/API_SP.DYN.CDRT.IN_DS2_en_csv_v2_241.csv  (World Bank crude
##                                         death rate per 1,000 people, by country and year)
##
## Outputs:
##   output/pm_joined/pop_pm_country_death_final.csv  (grid-cell table with death_rate_YYYY
##                                         columns for 1960-2020 joined on, as a proportion)
##
## Execution order:
##   files run before: country_map_pop_pm_final.R --> writes pop_pm_with_countries_final.csv
##   files run after: pop_wght_pm_cntry.R, fpm_gmt_regress_FE_grid_all.R,
##                    plots_stats_fpm_all.R --> read pop_pm_country_death_final.csv
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

# pop_pm_country <- read_csv(here("output", "pop_pm_with_countries_rev.csv")) #may use pop_pm_with_countries_rev.csv for latest version
pop_pm_country <- read_csv(here("output", "pm_joined", "pop_pm_with_countries_final.csv")) # includes pierce, park, zhao data

colnames(pop_pm_country)

# Import world bank death rates (World Bank Death rate, crude (per 1,000 people)) skipping metadata rows
wb_death_rate <- read_csv(
  here("input", "WB_crude_death_rate", "API_SP.DYN.CDRT.IN_DS2_en_csv_v2_241.csv"),
  skip = 4  # Skip the first 4 rows (metadata)
)

colnames(wb_death_rate)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Process death rate data ################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Clean and prepare death rate data: keep 1960-2020, rename columns, convert to proportion
wb_death_rate_clean <- wb_death_rate %>%
  # Select country info and years 1960-2020
  select(
    country_name = `Country Name`,
    country_code_iso3 = `Country Code`,
    `1960`:`2020`
  ) %>%
  # Rename year columns to death_rate_YYYY format and divide by 1000 (convert to proportion)
  mutate(across( # applies the same operation to multiple columns at once.
    .cols = `1960`:`2020`,  # .cols = which columns to transform (years 1960-2020)
    .fns = ~ .x / 1000,  # .fns = function to apply (.x represents each column's values)
    .names = "death_rate_{.col}"  # .names = naming pattern for new columns ({.col} = original column name)
  )) %>%
  # Remove original year columns (keep only renamed ones)
  select(-(`1960`:`2020`))

# Verify structure
cat("Rows:", nrow(wb_death_rate_clean), "\n")
cat("Columns:", ncol(wb_death_rate_clean), "\n")
colnames(wb_death_rate_clean)

# Check USA death rate in 2001
usa_death_2001 <- wb_death_rate_clean %>%
  filter(country_code_iso3 == "USA") %>%  # Filter to United States
  select(country_name, death_rate_2001)  # Select country name and 2001 death rate

print(usa_death_2001) # should be about 0.0085

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Compare death rate data to pm_pop ################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Compare country codes between pop_pm_country and wb_death_rate_clean
# to identify matches and mismatches

# Get unique country codes from each dataset
codes_pop_pm <- pop_pm_country %>%
  distinct(country_code_iso3, country_name) %>%  # Get unique country codes and names from pop_pm data
  arrange(country_code_iso3)  # Sort alphabetically

codes_wb <- wb_death_rate_clean %>%
  distinct(country_code_iso3, country_name) %>%  # Get unique country codes and names from World Bank data
  arrange(country_code_iso3)  # Sort alphabetically

# Find countries in pop_pm_country but NOT in wb_death_rate_clean
codes_in_pop_not_wb <- codes_pop_pm %>%
  anti_join(codes_wb, by = "country_code_iso3")  # anti_join keeps rows from left table that have no match in right table

# Find countries in wb_death_rate_clean but NOT in pop_pm_country
codes_in_wb_not_pop <- codes_wb %>%
  anti_join(codes_pop_pm, by = "country_code_iso3")  # anti_join keeps rows from left table that have no match in right table

# Find countries present in BOTH datasets
codes_matched <- codes_pop_pm %>%
  inner_join(codes_wb, by = "country_code_iso3", suffix = c("_pop", "_wb"))  # inner_join keeps only rows that match in both tables

# Print summary
cat("\n=== Country Code Comparison ===\n")
cat("Countries in pop_pm_country:", nrow(codes_pop_pm), "\n")
cat("Countries in wb_death_rate_clean:", nrow(codes_wb), "\n")
cat("Countries matched:", nrow(codes_matched), "\n")
cat("Countries in pop_pm but NOT in WB:", nrow(codes_in_pop_not_wb), "\n")
cat("Countries in WB but NOT in pop_pm:", nrow(codes_in_wb_not_pop), "\n")

# View the mismatches
cat("\n=== Countries in pop_pm_country but NOT in wb_death_rate_clean ===\n")
print(codes_in_pop_not_wb)

cat("\n=== Countries in wb_death_rate_clean but NOT in pop_pm_country ===\n")
print(codes_in_wb_not_pop)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Join death rate data to pop_pm_country ################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Left join wb_death_rate_clean to pop_pm_country, retaining all grid cells
pop_pm_country_death <- pop_pm_country %>%
  left_join(
    wb_death_rate_clean,  # Join World Bank death rate data
    by = "country_code_iso3"  # Match on ISO3 country code
  )

pop_pm_country_death <- pop_pm_country_death %>%
  select(-country_name.y) %>%
  rename(country_name = country_name.x)

# Check for grid cells with missing death rates
missing_death_rates <- pop_pm_country_death %>%
  filter(is.na(death_rate_2001) & !is.na(country_name))  # Grid cells with country but no death rate

cat("Grid cells with missing death rates:", nrow(missing_death_rates), "\n")

# Calculate world average death rates for each year (1960-2020)
world_death_rates <- wb_death_rate_clean %>%
  filter(country_code_iso3 == "WLD") %>%  # Filter to world aggregate (WLD = World)
  select(starts_with("death_rate_"))  # Select only death rate columns

# For grid cells without country-specific death rates, assign world average
pop_pm_country_death <- pop_pm_country_death %>%
  mutate(across(
    .cols = starts_with("death_rate_"),  # Apply to all death_rate_YYYY columns
    .fns = ~ if_else(
      is.na(.) & !is.na(country_name),  # If death rate is NA AND country name exists
      world_death_rates[[cur_column()]],  # Assign world rate for this year (cur_column() gets current column name)
      .  # Otherwise keep existing value
    )
  ))

# Verify no missing death rates remain for grid cells with countries
remaining_missing <- pop_pm_country_death %>%
  filter(is.na(death_rate_2001) & !is.na(country_name))

cat("Grid cells with missing death rates after world assignment:", nrow(remaining_missing), "\n")

# Check 2001 death rates for Somaliland and Falkland Islands
# two countries missing from WB data

# World rate for comparison
world_wb <- wb_death_rate_clean %>%
  filter(country_code_iso3 == "WLD") %>%
  select(country_name, death_rate_2001)
cat("\nWorld (WLD):\n")
print(world_wb)

# In pop_pm_country_death (after join and world rate assignment)
cat("\n=== Death rates in pop_pm_country_death ===\n")

# countries dropped (not in GIVE)
# somaliland_final <- pop_pm_country_death %>%
#   filter(country_code_iso3 == "-99" & country_name == "Somaliland") %>%
#   distinct(country_code_iso3, country_name, death_rate_2001)
# cat("Somaliland:\n")
# print(somaliland_final)
# 
# flk_final <- pop_pm_country_death %>%
#   filter(country_code_iso3 == "FLK") %>%
#   distinct(country_code_iso3, country_name, death_rate_2001)
# cat("\nFalkland Islands:\n")
# print(flk_final)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ save output ################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# write_csv(pop_pm_country_death, here("output", "pop_pm_country_death.csv"))
# write_csv(pop_pm_country_death, here("output", "pop_pm_country_death_park.csv"))
write_csv(pop_pm_country_death, here("output", "pm_joined", "pop_pm_country_death_final.csv"))


## THE END 
