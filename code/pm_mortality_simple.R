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
    pm_2000_mort = pop_tot_2009 * (1 - exp(-beta_k * delta_x_2000)) * Y,
    
    # 2000 baseline (no fire)
    delta_x_2000_nf = pmax(0, pm_2000_nf - threshold_crouse),
    pm_2000_nf_mort = pop_tot_2009 * (1 - exp(-beta_k * delta_x_2000_nf)) * Y,
    
    # 2050 RCP 4.5 (with fire)
    delta_x_2050_45 = pmax(0, pm_2050_45 - threshold_crouse),
    pm_2050_45_mort = pop_tot_2009 * (1 - exp(-beta_k * delta_x_2050_45)) * Y,
    
    # 2050 RCP 4.5 (no fire)
    delta_x_2050_45_nf = pmax(0, pm_2050_45_nf - threshold_crouse),
    pm_2050_45_nf_mort = pop_tot_2009 * (1 - exp(-beta_k * delta_x_2050_45_nf)) * Y,
    
    # 2050 RCP 8.5 (with fire)
    delta_x_2050_85 = pmax(0, pm_2050_85 - threshold_crouse),
    pm_2050_85_mort = pop_tot_2009 * (1 - exp(-beta_k * delta_x_2050_85)) * Y,
    
    # 2050 RCP 8.5 (no fire)
    delta_x_2050_85_nf = pmax(0, pm_2050_85_nf - threshold_crouse),
    pm_2050_85_nf_mort = pop_tot_2009 * (1 - exp(-beta_k * delta_x_2050_85_nf)) * Y,
    
    # 2100 RCP 4.5 (with fire)
    delta_x_2100_45 = pmax(0, pm_2100_45 - threshold_crouse),
    pm_2100_45_mort = pop_tot_2009 * (1 - exp(-beta_k * delta_x_2100_45)) * Y,
    
    # 2100 RCP 4.5 (no fire)
    delta_x_2100_45_nf = pmax(0, pm_2100_45_nf - threshold_crouse),
    pm_2100_45_nf_mort = pop_tot_2009 * (1 - exp(-beta_k * delta_x_2100_45_nf)) * Y,
    
    # 2100 RCP 8.5 (with fire)
    delta_x_2100_85 = pmax(0, pm_2100_85 - threshold_crouse),
    pm_2100_85_mort = pop_tot_2009 * (1 - exp(-beta_k * delta_x_2100_85)) * Y,
    
    # 2100 RCP 8.5 (no fire)
    delta_x_2100_85_nf = pmax(0, pm_2100_85_nf - threshold_crouse),
    pm_2100_85_nf_mort = pop_tot_2009 * (1 - exp(-beta_k * delta_x_2100_85_nf)) * Y
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
############ Calculate Fire PM Mortality (Method 2: Direct Difference) #################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Fire PM mortality = Total PM mortality - No Fire PM mortality
# This is Method 2 from Ford et al.: fpm_mortality = pm_total_mortality - pm_noFire_mortality

pop_pm_country <- pop_pm_country %>%
  mutate(
    # 2000 baseline
    fpm_2000_mort = pm_2000_mort - pm_2000_nf_mort,
    
    # 2050 RCP 4.5
    fpm_2050_45_mort = pm_2050_45_mort - pm_2050_45_nf_mort,
    
    # 2050 RCP 8.5
    fpm_2050_85_mort = pm_2050_85_mort - pm_2050_85_nf_mort,
    
    # 2100 RCP 4.5
    fpm_2100_45_mort = pm_2100_45_mort - pm_2100_45_nf_mort,
    
    # 2100 RCP 8.5
    fpm_2100_85_mort = pm_2100_85_mort - pm_2100_85_nf_mort
  )

# Check fire PM mortality results
summary(pop_pm_country$fpm_2000_mort)

# View top fire mortality grid cells for 2000
pop_pm_country %>%
  select(lon, lat, country_name, fpm_2000, pm_2000_mort, pm_2000_nf_mort, fpm_2000_mort) %>%
  filter(!is.na(fpm_2000_mort) & fpm_2000_mort > 0) %>%
  arrange(desc(fpm_2000_mort)) %>%
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
    # 2000 baseline
    total_2000 = sum(pm_2000_mort, na.rm = TRUE),
    total_2000_nf = sum(pm_2000_nf_mort, na.rm = TRUE),
    
    # 2050 scenarios
    total_2050_45 = sum(pm_2050_45_mort, na.rm = TRUE),
    total_2050_45_nf = sum(pm_2050_45_nf_mort, na.rm = TRUE),
    total_2050_85 = sum(pm_2050_85_mort, na.rm = TRUE),
    total_2050_85_nf = sum(pm_2050_85_nf_mort, na.rm = TRUE),
    
    # 2100 scenarios
    total_2100_45 = sum(pm_2100_45_mort, na.rm = TRUE),
    total_2100_45_nf = sum(pm_2100_45_nf_mort, na.rm = TRUE),
    total_2100_85 = sum(pm_2100_85_mort, na.rm = TRUE),
    total_2100_85_nf = sum(pm_2100_85_nf_mort, na.rm = TRUE)
  )

# Display results
print("Total PM Mortality (USA) - All Scenarios:")
usa_total_pm_mort

# Reshape for easier comparison
usa_total_pm_mort_long <- usa_total_pm_mort %>%
  pivot_longer(
    cols = everything(),
    names_to = "scenario",
    values_to = "total_deaths"
  ) %>%
  # First, extract fire scenario before splitting
  mutate(
    fire_scenario = if_else(str_detect(scenario, "_nf$"), "No Fire", "With Fire"),
    scenario_clean = str_remove(scenario, "^total_"),
    scenario_clean = str_remove(scenario_clean, "_nf$")
  ) %>%
  separate(scenario_clean, into = c("year", "rcp"), sep = "_", fill = "right") %>%
  select(year, rcp, total_deaths, fire_scenario)

print(usa_total_pm_mort_long)

# Lit comparison
cat("\nFord et al 2018 Table 3 reported 138,000 deaths for 2000 baseline (APPROACH #1) \n")
cat("\nPierce et al 2017 Figure 13 reported ~ 165,000 deaths for 2000 baseline (APPROACH #2) \n")
cat("Using DECADAL average placeholder, our 2000 with fire estimate:", round(usa_total_pm_mort$total_2000), "\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Fire PM Mortality Analysis (USA) ###########################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Calculate total fire-attributable deaths
usa_fpm_mort <- pop_pm_country_usa %>%
  summarise(
    fpm_2000 = sum(fpm_2000_mort, na.rm = TRUE),
    fpm_2050_45 = sum(fpm_2050_45_mort, na.rm = TRUE),
    fpm_2050_85 = sum(fpm_2050_85_mort, na.rm = TRUE),
    fpm_2100_45 = sum(fpm_2100_45_mort, na.rm = TRUE),
    fpm_2100_85 = sum(fpm_2100_85_mort, na.rm = TRUE)
  )

print("\nFire PM Mortality (USA) - All Scenarios:")
usa_fpm_mort

# Reshape for comparison
usa_fpm_mort_long <- usa_fpm_mort %>%
  pivot_longer(
    cols = everything(),
    names_to = "scenario",
    values_to = "fire_deaths"
  ) %>%
  separate(scenario, into = c("fpm", "year", "rcp"), sep = "_", fill = "right") %>%
  select(-fpm)

print(usa_fpm_mort_long)

# Ford et al 2018 Table 3 comparison: 17,000 fire deaths for 2000 baseline
cat("\nFord et al 2018 reported 17,000 fire-attributable deaths for 2000 baseline\n")
cat("Our 2000 fire estimate:", round(usa_fpm_mort$fpm_2000), "\n")

# Calculate percentage of total mortality from fire
usa_fire_pct <- usa_fpm_mort_long %>%
  left_join(
    usa_total_pm_mort_long %>% 
      filter(fire_scenario == "With Fire") %>%
      select(year, rcp, total_deaths = deaths),
    by = c("year", "rcp")
  ) %>%
  mutate(fire_pct = (fire_deaths / total_deaths) * 100)

print("\nFire as % of Total PM Mortality:")
print(usa_fire_pct)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Temporal & RCP Comparison ##################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Compare changes over time
cat("\n=== Changes in Total PM Mortality (With Fire) ===\n")
cat("2000 baseline:", round(usa_total_pm_mort$total_2000), "\n")
cat("2050 RCP4.5:", round(usa_total_pm_mort$total_2050_45), 
    sprintf("(%.1f%% change)", (usa_total_pm_mort$total_2050_45/usa_total_pm_mort$total_2000 - 1)*100), "\n")
cat("2050 RCP8.5:", round(usa_total_pm_mort$total_2050_85), 
    sprintf("(%.1f%% change)", (usa_total_pm_mort$total_2050_85/usa_total_pm_mort$total_2000 - 1)*100), "\n")
cat("2100 RCP4.5:", round(usa_total_pm_mort$total_2100_45), 
    sprintf("(%.1f%% change)", (usa_total_pm_mort$total_2100_45/usa_total_pm_mort$total_2000 - 1)*100), "\n")
cat("2100 RCP8.5:", round(usa_total_pm_mort$total_2100_85), 
    sprintf("(%.1f%% change)", (usa_total_pm_mort$total_2100_85/usa_total_pm_mort$total_2000 - 1)*100), "\n")

cat("\n=== Changes in Fire PM Mortality ===\n")
cat("2000 baseline:", round(usa_fpm_mort$fpm_2000), "\n")
cat("2050 RCP4.5:", round(usa_fpm_mort$fpm_2050_45), 
    sprintf("(%.1f%% change)", (usa_fpm_mort$fpm_2050_45/usa_fpm_mort$fpm_2000 - 1)*100), "\n")
cat("2050 RCP8.5:", round(usa_fpm_mort$fpm_2050_85), 
    sprintf("(%.1f%% change)", (usa_fpm_mort$fpm_2050_85/usa_fpm_mort$fpm_2000 - 1)*100), "\n")
cat("2100 RCP4.5:", round(usa_fpm_mort$fpm_2100_45), 
    sprintf("(%.1f%% change)", (usa_fpm_mort$fpm_2100_45/usa_fpm_mort$fpm_2000 - 1)*100), "\n")
cat("2100 RCP8.5:", round(usa_fpm_mort$fpm_2100_85), 
    sprintf("(%.1f%% change)", (usa_fpm_mort$fpm_2100_85/usa_fpm_mort$fpm_2000 - 1)*100), "\n")

cat("\n=== RCP Comparison at 2050 ===\n")
cat("Total PM mortality - RCP4.5:", round(usa_total_pm_mort$total_2050_45), "\n")
cat("Total PM mortality - RCP8.5:", round(usa_total_pm_mort$total_2050_85), 
    sprintf("(%.1f%% higher)", (usa_total_pm_mort$total_2050_85/usa_total_pm_mort$total_2050_45 - 1)*100), "\n")
cat("Fire PM mortality - RCP4.5:", round(usa_fpm_mort$fpm_2050_45), "\n")
cat("Fire PM mortality - RCP8.5:", round(usa_fpm_mort$fpm_2050_85), 
    sprintf("(%.1f%% higher)", (usa_fpm_mort$fpm_2050_85/usa_fpm_mort$fpm_2050_45 - 1)*100), "\n")

cat("\n=== RCP Comparison at 2100 ===\n")
cat("Total PM mortality - RCP4.5:", round(usa_total_pm_mort$total_2100_45), "\n")
cat("Total PM mortality - RCP8.5:", round(usa_total_pm_mort$total_2100_85), 
    sprintf("(%.1f%% higher)", (usa_total_pm_mort$total_2100_85/usa_total_pm_mort$total_2100_45 - 1)*100), "\n")
cat("Fire PM mortality - RCP4.5:", round(usa_fpm_mort$fpm_2100_45), "\n")
cat("Fire PM mortality - RCP8.5:", round(usa_fpm_mort$fpm_2100_85), 
    sprintf("(%.1f%% higher)", (usa_fpm_mort$fpm_2100_85/usa_fpm_mort$fpm_2100_45 - 1)*100), "\n")