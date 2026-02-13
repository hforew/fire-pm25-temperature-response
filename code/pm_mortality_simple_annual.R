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

# Import combined data
pop_pm_country <- read_csv(here("output", "pop_pm_with_countries.csv"))
colnames(pop_pm_country)

# Import world bank death rates (World Bank Death rate, crude (per 1,000 people)) skipping metadata rows
wb_death_rate <- read_csv(
  here("input", "WB_crude_death_rate", "API_SP.DYN.CDRT.IN_DS2_en_csv_v2_241.csv"),
  skip = 4  # Skip the first 4 rows (metadata)
)

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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Process death rate data and join to pm / pop ################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Create new dataframe with 2011-2020 average
wb_death_rate_2011_20 <- wb_death_rate %>%
  filter(year >= 2011 & year <= 2020) %>%
  group_by(country_name, country_code) %>%
  summarise(
    death_2011_20 = mean(death_rate, na.rm = TRUE),
    .groups = "drop"
  )

# Check the result
head(wb_death_rate_2011_20)

# USA should now be around 0.0085
wb_death_rate_2011_20 %>%
  filter(country_code == "USA")

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
beta_k <- log(RR_krewski)/10
threshold_crouse <- 1.9 # units μg/m^3

# fpm = pm_total - pm_noFire, where pm is simulated with and without fire
# Two approaches compute fpm mortality 
# 1: fpm_mortality = pm_total_mortality x (fpm/ pm_total)
# 2: fpm_mortality = pm_total_mortality - pm_noFire_mortality 
# Approach 2 is executed below

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Calculate Mortality for All PM Scenarios #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Calculate ΔX = max(0, PM2.5 - threshold)
# ΔMortality = Pop × (1 - exp(-β×ΔX)) × Y

# PM columns to process:
# pm_2000, pm_2000_nf, pm_2050_45, pm_2050_45_nf, pm_2050_85, pm_2050_85_nf,
# pm_2100_45, pm_2100_45_nf, pm_2100_85, pm_2100_85_nf
# population is held constant at the 2009 level.
# ΔX varies

pop_pm_country <- pop_pm_country %>%
  mutate(
    # 2000 baseline (with fire)
    delta_x_2000 = pmax(0, pm_2000 - threshold_crouse),
    pm_2000_mort = pop_tot_2009 * (1 - exp(-beta_k * delta_x_2000)) * death_2011_20,
    
    # 2000 baseline (no fire)
    delta_x_2000_nf = pmax(0, pm_2000_nf - threshold_crouse),
    pm_2000_nf_mort = pop_tot_2009 * (1 - exp(-beta_k * delta_x_2000_nf)) * death_2011_20,
    
    # 2050 RCP 4.5 (with fire)
    delta_x_2050_45 = pmax(0, pm_2050_45 - threshold_crouse),
    pm_2050_45_mort = pop_tot_2009 * (1 - exp(-beta_k * delta_x_2050_45)) * death_2011_20,
    
    # 2050 RCP 4.5 (no fire)
    delta_x_2050_45_nf = pmax(0, pm_2050_45_nf - threshold_crouse),
    pm_2050_45_nf_mort = pop_tot_2009 * (1 - exp(-beta_k * delta_x_2050_45_nf)) * death_2011_20,
    
    # 2050 RCP 8.5 (with fire)
    delta_x_2050_85 = pmax(0, pm_2050_85 - threshold_crouse),
    pm_2050_85_mort = pop_tot_2009 * (1 - exp(-beta_k * delta_x_2050_85)) * death_2011_20,
    
    # 2050 RCP 8.5 (no fire)
    delta_x_2050_85_nf = pmax(0, pm_2050_85_nf - threshold_crouse),
    pm_2050_85_nf_mort = pop_tot_2009 * (1 - exp(-beta_k * delta_x_2050_85_nf)) * death_2011_20,
    
    # 2100 RCP 4.5 (with fire)
    delta_x_2100_45 = pmax(0, pm_2100_45 - threshold_crouse),
    pm_2100_45_mort = pop_tot_2009 * (1 - exp(-beta_k * delta_x_2100_45)) * death_2011_20,
    
    # 2100 RCP 4.5 (no fire)
    delta_x_2100_45_nf = pmax(0, pm_2100_45_nf - threshold_crouse),
    pm_2100_45_nf_mort = pop_tot_2009 * (1 - exp(-beta_k * delta_x_2100_45_nf)) * death_2011_20,
    
    # 2100 RCP 8.5 (with fire)
    delta_x_2100_85 = pmax(0, pm_2100_85 - threshold_crouse),
    pm_2100_85_mort = pop_tot_2009 * (1 - exp(-beta_k * delta_x_2100_85)) * death_2011_20,
    
    # 2100 RCP 8.5 (no fire)
    delta_x_2100_85_nf = pmax(0, pm_2100_85_nf - threshold_crouse),
    pm_2100_85_nf_mort = pop_tot_2009 * (1 - exp(-beta_k * delta_x_2100_85_nf)) * death_2011_20
  )

# Check results for 2000 baseline
summary(pop_pm_country$pm_2000_mort)

# View sample of results
pop_pm_country %>%
  select(lon, lat, country_name, pop_tot_2009, pm_2000, delta_x_2000, pm_2000_mort) %>%
  filter(!is.na(pm_2000_mort) & pm_2000_mort > 0) %>%
  arrange(desc(pm_2000_mort)) %>%
  head(10)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Calculate Fire PM Mortality (Method 1: fraction) #################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Method 1: fpm_mortality = pm_total_mortality x (fpm/ pm_total)

pop_pm_country <- pop_pm_country %>%
  mutate(
    # 2000 baseline
    fpm_2000_mort_m1 = pm_2000_mort * (fpm_2000 / pm_2000),
    
    # 2050 RCP 4.5
    fpm_2050_45_mort_m1 = pm_2050_45_mort * (fpm_2050_45 / pm_2050_45),
    
    # 2050 RCP 8.5
    fpm_2050_85_mort_m1 = pm_2050_85_mort * (fpm_2050_85 / pm_2050_85),
    
    # 2100 RCP 4.5
    fpm_2100_45_mort_m1 = pm_2100_45_mort * (fpm_2100_45 / pm_2100_45),
    
    # 2100 RCP 8.5
    fpm_2100_85_mort_m1 = pm_2100_85_mort * (fpm_2100_85 / pm_2100_85)
  )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Calculate Fire PM Mortality (Method 2: Direct Difference) #################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Method 2: Fire PM mortality = Total PM mortality - No Fire PM mortality

pop_pm_country <- pop_pm_country %>%
  mutate(
    # 2000 baseline
    fpm_2000_mort_m2 = pm_2000_mort - pm_2000_nf_mort,
    
    # 2050 RCP 4.5
    fpm_2050_45_mort_m2 = pm_2050_45_mort - pm_2050_45_nf_mort,
    
    # 2050 RCP 8.5
    fpm_2050_85_mort_m2 = pm_2050_85_mort - pm_2050_85_nf_mort,
    
    # 2100 RCP 4.5
    fpm_2100_45_mort_m2 = pm_2100_45_mort - pm_2100_45_nf_mort,
    
    # 2100 RCP 8.5
    fpm_2100_85_mort_m2 = pm_2100_85_mort - pm_2100_85_nf_mort
  )

# Check fire PM mortality results
summary(pop_pm_country$fpm_2000_mort_m2)

# View top fire mortality grid cells for 2000
pop_pm_country %>%
  select(lon, lat, country_name, fpm_2000, pm_2000_mort, pm_2000_nf_mort, fpm_2000_mort_m2) %>%
  filter(!is.na(fpm_2000_mort_m2) & fpm_2000_mort_m2 > 0) %>%
  arrange(desc(fpm_2000_mort_m2)) %>%
  head(10)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Filter Data to USA ONLY #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Filter to USA only
pop_pm_country_usa <- pop_pm_country %>%
  filter(country_code_iso3 == "USA")

# Check how many USA grid cells
nrow(pop_pm_country_usa)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Total PM Mortality Analysis (USA) #########################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Calculate total deaths across all USA grid cells for each scenario
usa_total_pm_mort <- pop_pm_country_usa %>%
  summarise(
    total_2000 = sum(pm_2000_mort, na.rm = TRUE),
    total_2000_nf = sum(pm_2000_nf_mort, na.rm = TRUE),
    total_2050_45 = sum(pm_2050_45_mort, na.rm = TRUE),
    total_2050_45_nf = sum(pm_2050_45_nf_mort, na.rm = TRUE),
    total_2050_85 = sum(pm_2050_85_mort, na.rm = TRUE),
    total_2050_85_nf = sum(pm_2050_85_nf_mort, na.rm = TRUE),
    total_2100_45 = sum(pm_2100_45_mort, na.rm = TRUE),
    total_2100_45_nf = sum(pm_2100_45_nf_mort, na.rm = TRUE),
    total_2100_85 = sum(pm_2100_85_mort, na.rm = TRUE),
    total_2100_85_nf = sum(pm_2100_85_nf_mort, na.rm = TRUE)
  ) %>%
  # Reshape to long format immediately
  pivot_longer(
    cols = everything(),
    names_to = "scenario",
    values_to = "pm_deaths"
  ) %>%
  mutate(
    fire_scenario = if_else(str_detect(scenario, "_nf$"), "No Fire", "With Fire"),
    scenario_clean = str_remove(scenario, "^total_") %>% str_remove("_nf$")
  ) %>%
  separate(scenario_clean, into = c("year", "rcp"), sep = "_", fill = "right") %>%
  select(year, rcp, fire_scenario, pm_deaths)

print("Total PM Mortality (USA) - All Scenarios:")
print(usa_total_pm_mort)

# Literature comparison
cat("\nLiterature Comparison (2000 baseline):\n")
cat("Ford et al 2018 Table 3: 138,000 deaths (Method 1)\n")
cat("Pierce et al 2017 Figure 13: ~165,000 deaths (Method 2)\n")
cat("Our estimate (with fire):", 
    round(usa_total_pm_mort %>% filter(year == "2000", fire_scenario == "With Fire") %>% pull(pm_deaths)), "\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Fire PM Mortality Analysis (USA) ###########################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Calculate fire-attributable deaths for both methods
usa_fpm_mort <- pop_pm_country_usa %>%
  summarise(
    fpm_2000_m1 = sum(fpm_2000_mort_m1, na.rm = TRUE),
    fpm_2000_m2 = sum(fpm_2000_mort_m2, na.rm = TRUE),
    fpm_2050_45_m1 = sum(fpm_2050_45_mort_m1, na.rm = TRUE),
    fpm_2050_45_m2 = sum(fpm_2050_45_mort_m2, na.rm = TRUE),
    fpm_2050_85_m1 = sum(fpm_2050_85_mort_m1, na.rm = TRUE),
    fpm_2050_85_m2 = sum(fpm_2050_85_mort_m2, na.rm = TRUE),
    fpm_2100_45_m1 = sum(fpm_2100_45_mort_m1, na.rm = TRUE),
    fpm_2100_45_m2 = sum(fpm_2100_45_mort_m2, na.rm = TRUE),
    fpm_2100_85_m1 = sum(fpm_2100_85_mort_m1, na.rm = TRUE),
    fpm_2100_85_m2 = sum(fpm_2100_85_mort_m2, na.rm = TRUE)
  ) %>%
  # Reshape to long format immediately
  pivot_longer(
    cols = everything(),
    names_to = "scenario",
    values_to = "fpm_deaths"
  ) %>%
  mutate(
    method = if_else(str_detect(scenario, "_m2$"), "Method 2", "Method 1"),
    scenario_clean = str_remove(scenario, "^fpm_") %>% str_remove("_m[12]$")
  ) %>%
  separate(scenario_clean, into = c("year", "rcp"), sep = "_", fill = "right") %>%
  select(year, rcp, method, fpm_deaths)

print("\nFire PM Mortality (USA) - Both Methods:")
print(usa_fpm_mort)

# Side-by-side method comparison
usa_fpm_comparison <- usa_fpm_mort %>%
  pivot_wider(names_from = method, values_from = fpm_deaths) %>%
  mutate(
    difference = `Method 1` - `Method 2`,
    pct_difference = (difference / `Method 2`) * 100
  )

print("\nMethod 1 vs Method 2 Comparison:")
print(usa_fpm_comparison)

# Ford validation
cat("\n=== Ford et al 2018 Validation (2000 baseline) ===\n")
cat("Ford et al reported: 17,000 fire-attributable deaths\n")
cat("Our Method 1:", round(usa_fpm_mort %>% filter(year == "2000", method == "Method 1") %>% pull(fpm_deaths)), "\n")
cat("Our Method 2:", round(usa_fpm_mort %>% filter(year == "2000", method == "Method 2") %>% pull(fpm_deaths)), "\n")
cat("Ford et al used Method 1 for main results\n")

# Fire as percentage of total mortality
usa_fire_pct <- usa_fpm_mort %>%
  filter(method == "Method 1") %>%
  left_join(
    usa_total_pm_mort %>% filter(fire_scenario == "With Fire"),
    by = c("year", "rcp")
  ) %>%
  mutate(fpm_pm_pct = (fpm_deaths / pm_deaths) * 100)

print("Fire as % of Total PM Mortality (Method 1):")
print(usa_fire_pct)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Temporal & RCP Trends ######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Create summary tables for cleaner output
temporal_trends <- usa_fpm_mort %>%
  pivot_wider(names_from = method, values_from = fpm_deaths) %>%
  arrange(year, rcp)

cat("\n=== Fire PM Mortality Trends ===\n")
print(temporal_trends)

# Calculate percent changes from 2000 baseline
baseline_m1 <- usa_fpm_mort %>% filter(year == "2000", method == "Method 1") %>% pull(fpm_deaths)
baseline_m2 <- usa_fpm_mort %>% filter(year == "2000", method == "Method 2") %>% pull(fpm_deaths)

temporal_pct_change <- temporal_trends %>%
  mutate(
    m1_pct_change = (`Method 1` / baseline_m1 - 1) * 100,
    m2_pct_change = (`Method 2` / baseline_m2 - 1) * 100
  )

cat("\n=== Percent Change from 2000 Baseline ===\n")
print(temporal_pct_change)

# RCP comparisons
rcp_comparison <- temporal_trends %>%
  filter(year %in% c("2050", "2100")) %>%
  group_by(year) %>%
  summarise(
    m1_rcp45 = `Method 1`[rcp == "45"],
    m1_rcp85 = `Method 1`[rcp == "85"],
    m1_rcp_diff_pct = (m1_rcp85 / m1_rcp45 - 1) * 100,
    m2_rcp45 = `Method 2`[rcp == "45"],
    m2_rcp85 = `Method 2`[rcp == "85"],
    m2_rcp_diff_pct = (m2_rcp85 / m2_rcp45 - 1) * 100
  )

cat("\n=== RCP 8.5 vs RCP 4.5 Comparison ===\n")
print(rcp_comparison)

cat("\n=== Method Comparison Summary ===\n")
cat("Method 1 (Fractional): fpm_mort = pm_total_mort × (fpm / pm_total)\n")
cat("Method 2 (Difference): fpm_mort = pm_total_mort - pm_noFire_mort\n")
cat("\nMethod 2 gives HIGHER estimates\n")
cat("\nFord et al 2018 used Method 1\n")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Comparison with Pierce and Ford (method 1 ref only) ######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Create comparison object for Ford et al and Pierce et al benchmarking
usa_ford_pierce_comp <- usa_fire_pct %>% # method 1 data
  select(year, rcp, fpm_deaths, pm_deaths, fpm_pm_pct)

# Literature results from Ford et al 2018 and Pierce et al 2017
lit_results <- tibble(
  year = c(2000, 2050, 2050, 2100, 2100),
  rcp = c(NA, "45", "85", "45", "85"),
  # Ford et al 2018 results (Method 1)
  fpm_ford = c(17000, 42000, 32000, 32000, 44000),
  pm_ford = c(138000, 114000, 105000, 75000, 88000),
  fpm_pm_pct_ford = c(0.12, 0.37, 0.30, 0.43, 0.50),
  # Pierce et al 2017 results 
  fpm_pierce = c(20000, NA, 47000, NA, 65000),
  pm_pierce = c(165000, NA, 145000, NA, 125000),
  fpm_pm_pct_pierce = c(0.12, NA, 0.32, NA, 0.52)
)

# join as cbind (rows align perfectly):
usa_ford_pierce_comp <- usa_fire_pct %>%
  select(year, rcp, fpm_deaths, pm_deaths, fpm_pm_pct) %>%
  bind_cols(
    lit_results %>% select(fpm_ford, pm_ford, fpm_pm_pct_ford, 
                           fpm_pierce, pm_pierce, fpm_pm_pct_pierce)
  )

print("\n=== USA vs Literature Comparison ===")
print(usa_ford_pierce_comp)

# Prepare data for plotting
plot_data <- usa_ford_pierce_comp %>%
  mutate(
    # Create a combined x-axis label
    scenario = case_when(
      year == "2000" ~ "2000",
      TRUE ~ paste0(year, "\nRCP", rcp)
    ),
    # Order scenarios chronologically
    scenario = factor(scenario, levels = c("2000", "2050\nRCP45", "2050\nRCP85", 
                                           "2100\nRCP45", "2100\nRCP85"))
  ) %>%
  select(scenario, fpm_deaths, pm_deaths, fpm_ford, pm_ford, 
         fpm_pierce, pm_pierce) %>%
  # Reshape for plotting
  pivot_longer(
    cols = -scenario,
    names_to = "variable",
    values_to = "deaths"
  ) %>%
  mutate(
    # Separate metric (fire vs total) and source (ours, ford, pierce)
    metric = if_else(str_detect(variable, "^fpm"), "Fire PM", "Total PM"),
    source = case_when(
      str_detect(variable, "ford") ~ "Ford et al 2018",
      str_detect(variable, "pierce") ~ "Pierce et al 2017",
      TRUE ~ "Our Estimate"
    )
  )


# Plot 1: Total PM Deaths 
p1 <- plot_data %>%
  filter(metric == "Total PM") %>%
  ggplot(aes(x = scenario, y = deaths, color = source)) +
  geom_point(size = 3) +
  scale_color_manual(
    values = c("Our Estimate" = "#2E86AB", 
               "Ford et al 2018" = "#A23B72", 
               "Pierce et al 2017" = "#F18F01")
  ) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Total PM2.5 Mortality - USA",
    subtitle = "Comparison with Literature",
    x = "Scenario",
    y = "Annual Deaths",
    color = "Source"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# Plot 2: Fire PM Deaths 
p2 <- plot_data %>%
  filter(metric == "Fire PM") %>%
  ggplot(aes(x = scenario, y = deaths, color = source)) +
  geom_point(size = 3) +
  scale_color_manual(
    values = c("Our Estimate" = "#2E86AB", 
               "Ford et al 2018" = "#A23B72", 
               "Pierce et al 2017" = "#F18F01")
  ) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Fire PM2.5 Mortality - USA",
    subtitle = "Comparison with Literature",
    x = "Scenario",
    y = "Annual Deaths",
    color = "Source"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# Display plots
print(p1)
print(p2)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Mortality by Country Analysis ##############################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Aggregate mortality by country (Method 1 only, excluding no-fire scenarios)
country_mortality <- pop_pm_country %>%
  # Exclude disputed territories (-99) and unmapped (NA)
  filter(!is.na(country_code_iso3) & country_code_iso3 != "-99") %>%
  group_by(country_code_iso3, country_name) %>%
  summarise(
    # Population
    total_pop = sum(pop_tot_2009, na.rm = TRUE),
    
    # Total PM mortality (with fire)
    pm_2000_mort = sum(pm_2000_mort, na.rm = TRUE),
    pm_2050_45_mort = sum(pm_2050_45_mort, na.rm = TRUE),
    pm_2050_85_mort = sum(pm_2050_85_mort, na.rm = TRUE),
    pm_2100_45_mort = sum(pm_2100_45_mort, na.rm = TRUE),
    pm_2100_85_mort = sum(pm_2100_85_mort, na.rm = TRUE),
    
    # Fire PM mortality (Method 1)
    fpm_2000_mort_m1 = sum(fpm_2000_mort_m1, na.rm = TRUE),
    fpm_2050_45_mort_m1 = sum(fpm_2050_45_mort_m1, na.rm = TRUE),
    fpm_2050_85_mort_m1 = sum(fpm_2050_85_mort_m1, na.rm = TRUE),
    fpm_2100_45_mort_m1 = sum(fpm_2100_45_mort_m1, na.rm = TRUE),
    fpm_2100_85_mort_m1 = sum(fpm_2100_85_mort_m1, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  # Filter out countries with no population
  filter(total_pop > 0) %>%
  # Calculate fire as % of total for 2000 baseline
  mutate(
    fire_pct_2000 = (fpm_2000_mort_m1 / pm_2000_mort) * 100
  ) %>%
  # Sort by total PM mortality in 2000
  arrange(desc(pm_2000_mort))

print("\n=== Top 20 Countries by Total PM Mortality (2000) ===")
print(country_mortality %>% head(20))

# Summary statistics
cat("\n=== Global Totals (All Countries) ===\n")
country_mortality %>%
  summarise(
    total_population = sum(total_pop),
    pm_2000 = sum(pm_2000_mort),
    pm_2050_45 = sum(pm_2050_45_mort),
    pm_2050_85 = sum(pm_2050_85_mort),
    pm_2100_45 = sum(pm_2100_45_mort),
    pm_2100_85 = sum(pm_2100_85_mort),
    fpm_2000_m1 = sum(fpm_2000_mort_m1),
    fpm_2050_45_m1 = sum(fpm_2050_45_mort_m1),
    fpm_2050_85_m1 = sum(fpm_2050_85_mort_m1),
    fpm_2100_45_m1 = sum(fpm_2100_45_mort_m1),
    fpm_2100_85_m1 = sum(fpm_2100_85_mort_m1)
  ) %>%
  print()

# Top 10 countries by fire PM mortality
print("\n=== Top 10 Countries by Fire PM Mortality (2000, Method 1) ===")
country_mortality %>%
  select(country_name, country_code_iso3, fpm_2000_mort_m1, pm_2000_mort, fire_pct_2000) %>%
  arrange(desc(fpm_2000_mort_m1)) %>%
  head(10) %>%
  print()

# Countries with highest fire mortality as % of total
print("\n=== Top 10 Countries by Fire as % of Total PM Mortality (2000) ===")
country_mortality %>%
  filter(pm_2000_mort > 100) %>%  # Filter to countries with meaningful mortality
  select(country_name, country_code_iso3, fpm_2000_mort_m1, pm_2000_mort, fire_pct_2000) %>%
  arrange(desc(fire_pct_2000)) %>%
  head(10) %>%
  print()

# Temporal trends for top countries
print("\n=== Fire PM Mortality Trends for Top 5 Countries ===")
country_mortality %>%
  select(country_name, starts_with("fpm_")) %>%
  head(5) %>%
  print()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Global Mortality Comparison Plots ##########################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Get global totals from country_mortality
global_totals <- country_mortality %>%
  summarise(
    fpm_deaths = sum(fpm_2000_mort_m1),
    pm_deaths = sum(pm_2000_mort)
  ) %>%
  mutate(
    year = "2000",
    rcp = NA_character_,
    fpm_pm_pct = (fpm_deaths / pm_deaths) * 100
  )

# Add future scenarios
global_totals_all <- bind_rows(
  global_totals,
  country_mortality %>%
    summarise(
      year = "2050",
      rcp = "45",
      fpm_deaths = sum(fpm_2050_45_mort_m1),
      pm_deaths = sum(pm_2050_45_mort),
      fpm_pm_pct = (fpm_deaths / pm_deaths) * 100
    ),
  country_mortality %>%
    summarise(
      year = "2050",
      rcp = "85",
      fpm_deaths = sum(fpm_2050_85_mort_m1),
      pm_deaths = sum(pm_2050_85_mort),
      fpm_pm_pct = (fpm_deaths / pm_deaths) * 100
    ),
  country_mortality %>%
    summarise(
      year = "2100",
      rcp = "45",
      fpm_deaths = sum(fpm_2100_45_mort_m1),
      pm_deaths = sum(pm_2100_45_mort),
      fpm_pm_pct = (fpm_deaths / pm_deaths) * 100
    ),
  country_mortality %>%
    summarise(
      year = "2100",
      rcp = "85",
      fpm_deaths = sum(fpm_2100_85_mort_m1),
      pm_deaths = sum(pm_2100_85_mort),
      fpm_pm_pct = (fpm_deaths / pm_deaths) * 100
    )
)

print("\n=== Global Mortality Totals ===")
print(global_totals_all)

# Add additional columns for context
global_totals_all <- global_totals_all %>%
  mutate(
    all_death_2009 = 54100000,  # Total global deaths in 2009
    gmt_chg = c(0.91, 1.4, 2, 1.8, 3.7),  # Global mean temperature change (°C) for each scenario
    death_pct_chg = (fpm_deaths / all_death_2009)  # Fire PM deaths as % of all global deaths (decimal form)
  )


print(global_totals_all)

# Prepare data for plotting (matching USA format)
plot_data_global <- global_totals_all %>%
  mutate(
    # Create a combined x-axis label
    scenario = case_when(
      year == "2000" ~ "2000",
      TRUE ~ paste0(year, "\nRCP", rcp)
    ),
    # Order scenarios chronologically
    scenario = factor(scenario, levels = c("2000", "2050\nRCP45", "2050\nRCP85", 
                                           "2100\nRCP45", "2100\nRCP85"))
  ) %>%
  select(scenario, fpm_deaths, pm_deaths) %>%
  # Reshape for plotting
  pivot_longer(
    cols = c(fpm_deaths, pm_deaths),
    names_to = "metric",
    values_to = "deaths"
  ) %>%
  mutate(
    metric = if_else(metric == "fpm_deaths", "Fire PM", "Total PM")
  )

# Plot 1: Total PM Deaths (Global)
p1_global <- plot_data_global %>%
  filter(metric == "Total PM") %>%
  ggplot(aes(x = scenario, y = deaths)) +
  geom_point(size = 3, color = "#2E86AB") +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Total PM2.5 Mortality - Global",
    subtitle = "Our Estimates",
    x = "Scenario",
    y = "Annual Deaths"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# Plot 2: Fire PM Deaths (Global)
p2_global <- plot_data_global %>%
  filter(metric == "Fire PM") %>%
  ggplot(aes(x = scenario, y = deaths)) +
  geom_point(size = 3, color = "#A23B72") +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Fire PM2.5 Mortality - Global",
    subtitle = "Our Estimates",
    x = "Scenario",
    y = "Annual Deaths"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# Display plots
print(p1_global)
print(p2_global)

# Optional: Save plots
# ggsave("total_pm_mortality_global.png", p1_global, width = 8, height = 6, dpi = 300)
# ggsave("fire_pm_mortality_global.png", p2_global, width = 8, height = 6, dpi = 300)

# Print summary comparison
cat("\n=== Global vs USA Comparison (2000 baseline) ===\n")
cat("Global total PM deaths:", round(global_totals_all$pm_deaths[1]), "\n")
cat("USA total PM deaths:", round(usa_total_pm_mort %>% filter(year == "2000", fire_scenario == "With Fire") %>% pull(pm_deaths)), "\n")
cat("USA as % of global:", 
    round((usa_total_pm_mort %>% filter(year == "2000", fire_scenario == "With Fire") %>% pull(pm_deaths)) / 
            global_totals_all$pm_deaths[1] * 100, 1), "%\n")
cat("\nGlobal fire PM deaths:", round(global_totals_all$fpm_deaths[1]), "\n")
cat("USA fire PM deaths:", round(usa_fpm_mort %>% filter(year == "2000", method == "Method 1") %>% pull(fpm_deaths)), "\n")
cat("USA as % of global:", 
    round((usa_fpm_mort %>% filter(year == "2000", method == "Method 1") %>% pull(fpm_deaths)) / 
            global_totals_all$fpm_deaths[1] * 100, 1), "%\n")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Global Death and Temperature Change ########################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Create plot of fire PM deaths vs global mean temperature change
p_temp <- global_totals_all %>%
  ggplot(aes(x = gmt_chg, y = fpm_deaths)) +
  geom_smooth(method = "lm", se = TRUE, color = "#A23B72", fill = "#A23B72", alpha = 0.2) +
  geom_point(size = 3, color = "#A23B72") +
  scale_y_continuous(labels = scales::comma) +
  scale_x_continuous(breaks = c(0.91, 1.4, 1.8, 2, 3.7)) +
  labs(
    title = "Fire PM2.5 Mortality vs Global Mean Temperature Change",
    subtitle = "Global Estimates",
    x = "Global Mean Temperature Change (°C)",
    y = "Annual Fire PM Deaths"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# Display plot
print(p_temp)


# Calculate and display linear regression statistics
temp_model <- lm(fpm_deaths ~ gmt_chg, data = global_totals_all)
cat("\n=== Linear Regression: Fire PM Deaths ~ Temperature Change ===\n")
cat("Equation: Fire PM Deaths = ", round(coef(temp_model)[1]), " + ", 
    round(coef(temp_model)[2]), " × GMT Change\n", sep = "")
cat("R-squared:", round(summary(temp_model)$r.squared, 3), "\n")
cat("P-value:", format.pval(summary(temp_model)$coefficients[2,4], digits = 3), "\n")
cat("\nInterpretation: Each 1°C increase in global mean temperature is associated with\n")
cat("approximately", round(coef(temp_model)[2]), "additional fire PM deaths per year.\n")


# Optional: Save plot
# ggsave("fire_pm_mortality_vs_temperature.png", p_temp, width = 8, height = 6, dpi = 300)

# Summary statistics
cat("\n=== Fire PM Deaths by Temperature Scenario ===\n")
global_totals_all %>%
  select(year, rcp, gmt_chg, fpm_deaths, death_pct_chg) %>%
  mutate(
    fpm_deaths = round(fpm_deaths),
    death_pct_chg = round(death_pct_chg, 3)
  ) %>%
  print()

cat("\n=== Temperature-Mortality Relationship ===\n")
cat("At +0.91°C (2000 baseline):", round(global_totals_all$fpm_deaths[1]), "deaths\n")
cat("At +1.40°C (2050 RCP4.5):", round(global_totals_all$fpm_deaths[2]), "deaths\n")
cat("At +1.80°C (2100 RCP4.5):", round(global_totals_all$fpm_deaths[4]), "deaths\n")
cat("At +2.00°C (2050 RCP8.5):", round(global_totals_all$fpm_deaths[3]), "deaths\n")
cat("At +3.70°C (2100 RCP8.5):", round(global_totals_all$fpm_deaths[5]), "deaths\n")
cat("\nIncrease from +0.91°C to +3.70°C:", 
    round(global_totals_all$fpm_deaths[5] - global_totals_all$fpm_deaths[1]), 
    "deaths (", 
    round((global_totals_all$fpm_deaths[5] / global_totals_all$fpm_deaths[1] - 1) * 100, 1), 
    "% increase)\n")



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ All-cause death % chg due to FPM and Temperature Change ###############################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Calculate linear regression model first
temp_pct_model <- lm(death_pct_chg ~ gmt_chg, data = global_totals_all)

# Display regression statistics
cat("\n=== Linear Regression: Fire PM % of All Deaths ~ Temperature Change ===\n")
cat("Equation: Fire PM % = ", round(coef(temp_pct_model)[1], 5), " + ", 
    round(coef(temp_pct_model)[2], 5), " × GMT Change\n", sep = "")
cat("R-squared:", round(summary(temp_pct_model)$r.squared, 3), "\n")
cat("P-value:", format.pval(summary(temp_pct_model)$coefficients[2,4], digits = 3), "\n")
cat("\nInterpretation: Each 1°C increase in global mean temperature is associated with\n")
cat("an increase of", round(coef(temp_pct_model)[2], 5), "percentage points in the share of\n")
cat("all-cause deaths attributable to fire PM.\n")
cat("\nFor example, if fire PM deaths are currently 1% of all deaths, a 1°C warming would\n")
cat("increase this to approximately", round(1 + coef(temp_pct_model)[2], 5), "% of all deaths.\n")

# Create plot of fire PM deaths as % of all deaths vs global mean temperature change
# Add scenario labels to the data
global_totals_all <- global_totals_all %>%
  mutate(
    scenario_label = case_when(
      year == "2000" ~ "2000",
      year == "2050" & rcp == "45" ~ "RCP4.5 2050",
      year == "2050" & rcp == "85" ~ "RCP8.5 2050",
      year == "2100" & rcp == "45" ~ "RCP4.5 2100",
      year == "2100" & rcp == "85" ~ "RCP8.5 2100"
    )
  )


p_temp_pct <- global_totals_all %>%
  ggplot(aes(x = gmt_chg, y = death_pct_chg)) +
  geom_smooth(method = "lm", se = FALSE, color = "#A23B72") +
  geom_text(aes(label = scenario_label), # Add RCP labels for each point
            hjust = .5, vjust = -1, 
            size = 3.5, color = "#2E86AB") +
  geom_point(size = 3, color = "#A23B72") +
  scale_x_continuous(breaks = c(0.91, 1.4, 1.8, 2.0, 3.7)) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "All-cause death %Δ due to FPM vs Temperature Change",
    subtitle = paste0("All-cause death %Δ = ", round(coef(temp_pct_model)[1], 5), " + ", 
                      round(coef(temp_pct_model)[2], 5), " × GMT Change"),
    x = "Global Mean Temperature Change (°C)",
    y = "All-cause death %Δ due to FPM"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# Display plot
print(p_temp_pct)

# Optional: Save plot
# ggsave("fire_pm_pct_vs_temperature.png", p_temp_pct, width = 8, height = 6, dpi = 300)

# Summary statistics
cat("\n=== Fire PM Deaths as % of All Deaths by Temperature Scenario ===\n")
global_totals_all %>%
  select(year, rcp, gmt_chg, death_pct_chg) %>%
  mutate(
    death_pct_chg = round(death_pct_chg * 100, 3)
  ) %>%
  print()

cat("\n=== Temperature-Mortality Percentage Relationship ===\n")
cat("At +0.91°C (2000 baseline):", round(global_totals_all$death_pct_chg[1] * 100, 3), "% of all deaths\n")
cat("At +1.40°C (2050 RCP4.5):", round(global_totals_all$death_pct_chg[2] * 100, 3), "% of all deaths\n")
cat("At +2.00°C (2100 RCP4.5):", round(global_totals_all$death_pct_chg[3] * 100, 3), "% of all deaths\n")
cat("At +1.80°C (2050 RCP8.5):", round(global_totals_all$death_pct_chg[4] * 100, 3), "% of all deaths\n")
cat("At +3.70°C (2100 RCP8.5):", round(global_totals_all$death_pct_chg[5] * 100, 3), "% of all deaths\n")
cat("\nIncrease from +0.91°C to +3.70°C:", 
    round((global_totals_all$death_pct_chg[5] - global_totals_all$death_pct_chg[1]) * 100, 3), 
    "percentage points (", 
    round((global_totals_all$death_pct_chg[5] / global_totals_all$death_pct_chg[1] - 1) * 100, 1), 
    "% increase)\n")


#### THE END



