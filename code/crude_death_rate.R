# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Process Death Rate Data ##
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

# Import world bank death rates (World Bank Death rate, crude (per 1,000 people)) skipping metadata rows
wb_death_rate <- read_csv(
  here("input", "WB_crude_death_rate", "API_SP.DYN.CDRT.IN_DS2_en_csv_v2_241.csv"),
  skip = 4  # Skip the first 4 rows (metadata)
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Process death rate data and join to pm / pop ################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Reshape from wide to long format and keep only country name, code, year, and death rate
wb_death_rate <- wb_death_rate %>%
  select(`Country Name`, `Country Code`, starts_with("19"), starts_with("20")) %>%
  pivot_longer(
    cols = starts_with("19") | starts_with("20"),
    names_to = "year",
    values_to = "death_rate"
  ) %>%
  mutate(year = as.numeric(year)) %>%
  filter(!is.na(death_rate))  # Remove rows with missing death rates

# Check the result
head(wb_death_rate)
colnames(wb_death_rate)

# Rename columns
wb_death_rate <- wb_death_rate %>%
  rename(
    country_name = `Country Name`,
    country_code = `Country Code`
  ) %>%
  mutate(death_rate = death_rate / 1000)  # Convert to decimal form

head(wb_death_rate)


## OLD DROP BELOW 
# # Create new dataframe with 2011-2020 average
# wb_death_rate_2011_20 <- wb_death_rate %>%
#   filter(year >= 2011 & year <= 2020) %>%
#   group_by(country_name, country_code) %>%
#   summarise(
#     death_2011_20 = mean(death_rate, na.rm = TRUE),
#     .groups = "drop"
#   )

# Check the result
# head(wb_death_rate_2011_20)
## OLD DROP ABOVE 


# USA should now be around 0.0085
wb_death_rate %>%
  filter(country_code == "USA") %>%
  filter(year >= 2001 & year <= 2010)


# Count unique country codes
n_distinct(wb_death_rate_2011_20$country_code)
n_distinct(pop_pm_country$country_code_iso3)

# Create dataframe with unique country codes and names from pop_pm_country 
# left joined with death rates
country_death_rates <- pop_pm_country %>%
  distinct(country_code_iso3, country_name) %>%
  left_join(
    wb_death_rate_2011_20 %>% select(country_code, death_2011_20),
    by = c("country_code_iso3" = "country_code")
  )

# Check the result
head(country_death_rates)
nrow(country_death_rates) # 196 because -99 appears for certain countries 

# See which countries are missing death rate data
country_death_rates %>%
  filter(is.na(death_2011_20))

# Join pop_pm_country with death rate data

nrow(pop_pm_country)

pop_pm_country <- pop_pm_country %>%
  left_join(
    wb_death_rate_2011_20 %>% select(country_code, death_2011_20),
    by = c("country_code_iso3" = "country_code")
  )

# Check the result
nrow(pop_pm_country)
head(pop_pm_country)

# Check how many rows have death rate data
sum(!is.na(pop_pm_country$death_2011_20))
sum(is.na(pop_pm_country$death_2011_20))
sum(!is.na(pop_pm_country$death_2011_20)) /  sum(is.na(pop_pm_country$death_2011_20))

# Check if NAs correspond to low/zero population
pop_pm_country %>%
  summarise(
    # Grid cells WITH death rate
    pop_with_death = sum(pop_tot_2009[!is.na(death_2011_20)], na.rm = TRUE),
    # Grid cells WITHOUT death rate
    pop_without_death = sum(pop_tot_2009[is.na(death_2011_20)], na.rm = TRUE),
    # What % of population has death rate data?
    pct_pop_covered = pop_with_death / (pop_with_death + pop_without_death) * 100
  )



## THE END 
