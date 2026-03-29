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
# pop_pm_country <- read_csv(here("output", "pop_pm_with_countries.csv"))
pop_pm_country <- read_csv(here("output", "pop_pm_country_death.csv"))
colnames(pop_pm_country)

# GMT change
gmt_chg <- read_csv(here("output", "gmt_periods_pi.csv"))
print(gmt_chg)

# Extract GMT scalars for each PM scenario period (pre-industrial baseline)
gmt_baseline <- mean(c(gmt_chg$mean_gmt_45[gmt_chg$period == "2006-2010"],
                       gmt_chg$mean_gmt_85[gmt_chg$period == "2006-2010"]))
gmt_2050_45  <- gmt_chg$mean_gmt_45[gmt_chg$period == "2041-2050"]
gmt_2050_85  <- gmt_chg$mean_gmt_85[gmt_chg$period == "2041-2050"]
gmt_2100_45  <- gmt_chg$mean_gmt_45[gmt_chg$period == "2091-2100"]
gmt_2100_85  <- gmt_chg$mean_gmt_85[gmt_chg$period == "2091-2100"]

gmt_baseline
gmt_2050_45
gmt_2050_85
gmt_2100_45
gmt_2100_85

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Mortality Equation Setup #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# SHOULD BE ANNUAL BUT CURRENTLY DECADE AVERAGE 

# this code computes mortality from total PM as:
# ΔMort_PM = Pop (1 - exp^(-β*ΔX)) Y
# ΔMortality: increase in mortality from exposure to all PM or fire PM
# Pop: total population --> varies by cell and by year
# β = ln(RR)/ΔX (concentration response from literature)
# ΔX = PM2.5 concentration threshold above which fatal --> varies by cell and by year.
# ΔX[lon, lat, year] = max(0, PM2.5 - threshold) --> if negative, then 0
# Y = baseline all-cause mortality --> is constant by cell but varies by year

# fire pm is:
# fpm = pm_total - pm_noFire, where pm is simulated with and without fire
# ΔMort_fPM = ΔMort_PM * (fpm / pm_total)

# Parameters from Orellano et al  Table 1
RR_orellano <- 1.095
beta_k <- log(RR_orellano)/10
# threshold_crouse <- 1.9 # units μg/m^3
threshold_orellano <- 0 # units μg/m^3. orellano state: "no evidence for a threshold below which no effects are observed"

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


# in computing future mortality, population and mortality rate is assumed constant
# e.g. 2001 pop and rate is applied to 2041 and 2091. 2002 to 2042 and 2092, etc. 
# the decade average is take after computing the yearly values
# only the PM exposure changes. this isolates the climate effect, without pop changes. 
# pm_2000 refers to the annual average pm2.5 exposure in the 2000s decade from 2001-2010
# pm_2050 refers to the annual average pm2.5 exposure in the 2040s decade from 2041-2050 (under rcp45 and 85)
# pm_2100 refers to the annual average pm2.5 exposure in the 2090s decade from 2091-2100 (under rcp45 and 85)

pop_pm_country <- pop_pm_country %>%
  mutate(
    # 2000 baseline (with fire)
    delta_x_2000 = pmax(0, pm_2000 - threshold_orellano),
    pm_2001_mort = pop_tot_2001 * (1 - exp(-beta_k * delta_x_2000)) * death_rate_2001,
    
    # 2050 RCP 4.5 (with fire)
    delta_x_2050_45 = pmax(0, pm_2050_45 - threshold_orellano),
    pm_2041_45_mort = pop_tot_2001 * (1 - exp(-beta_k * delta_x_2050_45)) * death_rate_2001,

    # 2050 RCP 8.5 (with fire)
    delta_x_2050_85 = pmax(0, pm_2050_85 - threshold_orellano),
    pm_2041_85_mort = pop_tot_2001 * (1 - exp(-beta_k * delta_x_2050_85)) * death_rate_2001,
    
    # 2100 RCP 4.5 (with fire)
    delta_x_2100_45 = pmax(0, pm_2100_45 - threshold_orellano),
    pm_2091_45_mort = pop_tot_2001 * (1 - exp(-beta_k * delta_x_2100_45)) * death_rate_2001,

    # 2100 RCP 8.5 (with fire)
    delta_x_2100_85 = pmax(0, pm_2100_85 - threshold_orellano),
    pm_2091_85_mort = pop_tot_2001 * (1 - exp(-beta_k * delta_x_2100_85)) * death_rate_2001,
  )

# Check results for 2000 baseline
summary(pop_pm_country$pm_2001_mort)

# View sample of results
pop_pm_country %>%
  select(lon, lat, country_name.x, pop_tot_2001, pm_2000, delta_x_2000, pm_2001_mort) %>%
  filter(!is.na(pm_2001_mort) & pm_2001_mort > 0) %>%
  arrange(desc(pm_2001_mort)) %>%
  head(10)


### loop thru other years 

# Loop over base years 2002-2010 (the 2000s decade, excluding 2001 which is already computed)
for (yr in 2002:2010) {
  
  # Build dynamic column name strings for pop and death rate for this base year
  # e.g. yr=2002 -> "pop_tot_2002", "death_rate_2002"
  pop_col   <- paste0("pop_tot_", yr)
  death_col <- paste0("death_rate_", yr)
  
  # Compute the corresponding future years for the 2050 and 2100 PM decades
  # 2042 corresponds to 2002 pop/death, 2043 to 2003, etc.
  yr_50  <- yr + 40  # 2042:2050 (2050s PM decade)
  yr_100 <- yr + 90  # 2092:2100 (2100s PM decade)
  
  pop_pm_country <- pop_pm_country %>%
    mutate(
      
      # --- 200x BASELINE ---
      # Uses this year's pop and death rate, with 2000s PM exposure (delta_x_2000)
      # delta_x_2000 already computed above: pmax(0, pm_2000 - threshold_orellano)
      !!paste0("pm_", yr, "_mort")        := .data[[pop_col]] * (1 - exp(-beta_k * delta_x_2000))    * .data[[death_col]],
      
      # --- 204x RCP 4.5 ---
      # Same pop/death as base year (e.g. 2002), but 2050s PM under RCP 4.5
      # delta_x_2050_45 already computed above: pmax(0, pm_2050_45 - threshold_orellano)
      !!paste0("pm_", yr_50,  "_45_mort") := .data[[pop_col]] * (1 - exp(-beta_k * delta_x_2050_45)) * .data[[death_col]],
      
      # --- 204x RCP 8.5 ---
      # Same pop/death as base year, but 2050s PM under RCP 8.5
      # delta_x_2050_85 already computed above: pmax(0, pm_2050_85 - threshold_orellano)
      !!paste0("pm_", yr_50,  "_85_mort") := .data[[pop_col]] * (1 - exp(-beta_k * delta_x_2050_85)) * .data[[death_col]],
      
      # --- 209x RCP 4.5 ---
      # Same pop/death as base year (e.g. 2002), but 2100s PM under RCP 4.5
      # delta_x_2100_45 already computed above: pmax(0, pm_2100_45 - threshold_orellano)
      !!paste0("pm_", yr_100, "_45_mort") := .data[[pop_col]] * (1 - exp(-beta_k * delta_x_2100_45)) * .data[[death_col]],
      
      # --- 209x RCP 8.5 ---
      # Same pop/death as base year, but 2100s PM under RCP 8.5
      # delta_x_2100_85 already computed above: pmax(0, pm_2100_85 - threshold_orellano)
      !!paste0("pm_", yr_100, "_85_mort") := .data[[pop_col]] * (1 - exp(-beta_k * delta_x_2100_85)) * .data[[death_col]]
    )
}
# Result: 45 new columns added (9 years x 5 scenarios)
# e.g. pm_2002_mort, pm_2042_45_mort, pm_2042_85_mort, pm_2092_45_mort, pm_2092_85_mort ... through yr=2010

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Calculate Fire PM Mortality  #################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# fpm_mortality = pm_total_mortality x (fpm/ pm_total)
# This isolates the fire-attributable share of total PM mortality
pop_pm_country <- pop_pm_country %>%
  mutate(
    # 2001 baseline
    fpm_2001_mort = pm_2001_mort * (fpm_2000 / pm_2000),
    
    # 2041 RCP 4.5
    fpm_2041_45_mort = pm_2041_45_mort * (fpm_2050_45 / pm_2050_45),
    
    # 2041 RCP 8.5
    fpm_2041_85_mort = pm_2041_85_mort * (fpm_2050_85 / pm_2050_85),
    
    # 2091 RCP 4.5
    fpm_2091_45_mort = pm_2091_45_mort * (fpm_2100_45 / pm_2100_45),
    
    # 2091 RCP 8.5
    fpm_2091_85_mort = pm_2091_85_mort * (fpm_2100_85 / pm_2100_85)
  )

# Loop over base years 2002-2010 to compute fire PM mortality
for (yr in 2002:2010) {
  
  # Corresponding future years (same offset logic as pm_ loop above)
  yr_50  <- yr + 40  # 2042:2050
  yr_100 <- yr + 90  # 2092:2100
  
  pop_pm_country <- pop_pm_country %>%
    mutate(
      
      # --- 200x BASELINE ---
      # Fire share of 2000s PM mortality for this base year
      # fpm_2000 / pm_2000 = fraction of total PM that is fire PM (from 2000s decade)
      !!paste0("fpm_", yr, "_mort")        := .data[[paste0("pm_", yr, "_mort")]]        * (fpm_2000    / pm_2000),
      
      # --- 204x RCP 4.5 ---
      # Fire share of 2050s RCP 4.5 PM mortality for this base year's pop/death
      !!paste0("fpm_", yr_50,  "_45_mort") := .data[[paste0("pm_", yr_50,  "_45_mort")]] * (fpm_2050_45 / pm_2050_45),
      
      # --- 204x RCP 8.5 ---
      !!paste0("fpm_", yr_50,  "_85_mort") := .data[[paste0("pm_", yr_50,  "_85_mort")]] * (fpm_2050_85 / pm_2050_85),
      
      # --- 209x RCP 4.5 ---
      # Fire share of 2100s RCP 4.5 PM mortality for this base year's pop/death
      !!paste0("fpm_", yr_100, "_45_mort") := .data[[paste0("pm_", yr_100, "_45_mort")]] * (fpm_2100_45 / pm_2100_45),
      
      # --- 209x RCP 8.5 ---
      !!paste0("fpm_", yr_100, "_85_mort") := .data[[paste0("pm_", yr_100, "_85_mort")]] * (fpm_2100_85 / pm_2100_85)
    )
}
# Result: 45 new columns (9 years x 5 scenarios)
# e.g. fpm_2002_mort, fpm_2042_45_mort, fpm_2042_85_mort, fpm_2092_45_mort, fpm_2092_85_mort ... through yr=2010


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Decade Average  #################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pop_pm_country <- pop_pm_country %>%
  mutate(
    
    # --- DECADE AVERAGE: Total PM Mortality ---
    # Each cell gets the mean mortality across all 10 years in the decade
    # pick() selects columns by name dynamically; paste0() builds the column name vector
    
    # 2000s baseline (average of pm_2001_mort through pm_2010_mort)
    pm_2000s_mort = rowMeans(pick(all_of(paste0("pm_", 2001:2010, "_mort"))),       na.rm = TRUE),
    
    # 2040s RCP 4.5 (average of pm_2041_45_mort through pm_2050_45_mort)
    pm_2040s_45_mort = rowMeans(pick(all_of(paste0("pm_", 2041:2050, "_45_mort"))), na.rm = TRUE),
    
    # 2040s RCP 8.5 (average of pm_2041_85_mort through pm_2050_85_mort)
    pm_2040s_85_mort = rowMeans(pick(all_of(paste0("pm_", 2041:2050, "_85_mort"))), na.rm = TRUE),
    
    # 2090s RCP 4.5 (average of pm_2091_45_mort through pm_2100_45_mort)
    pm_2090s_45_mort = rowMeans(pick(all_of(paste0("pm_", 2091:2100, "_45_mort"))), na.rm = TRUE),
    
    # 2090s RCP 8.5 (average of pm_2091_85_mort through pm_2100_85_mort)
    pm_2090s_85_mort = rowMeans(pick(all_of(paste0("pm_", 2091:2100, "_85_mort"))), na.rm = TRUE),
    
    # --- DECADE AVERAGE: Fire PM Mortality ---
    # Same logic as above but for fire-attributable mortality columns
    
    # 2000s baseline (average of fpm_2001_mort through fpm_2010_mort)
    fpm_2000s_mort = rowMeans(pick(all_of(paste0("fpm_", 2001:2010, "_mort"))),       na.rm = TRUE),
    
    # 2040s RCP 4.5 (average of fpm_2041_45_mort through fpm_2050_45_mort)
    fpm_2040s_45_mort = rowMeans(pick(all_of(paste0("fpm_", 2041:2050, "_45_mort"))), na.rm = TRUE),
    
    # 2040s RCP 8.5 (average of fpm_2041_85_mort through fpm_2050_85_mort)
    fpm_2040s_85_mort = rowMeans(pick(all_of(paste0("fpm_", 2041:2050, "_85_mort"))), na.rm = TRUE),
    
    # 2090s RCP 4.5 (average of fpm_2091_45_mort through fpm_2100_45_mort)
    fpm_2090s_45_mort = rowMeans(pick(all_of(paste0("fpm_", 2091:2100, "_45_mort"))), na.rm = TRUE),
    
    # 2090s RCP 8.5 (average of fpm_2091_85_mort through fpm_2100_85_mort)
    fpm_2090s_85_mort = rowMeans(pick(all_of(paste0("fpm_", 2091:2100, "_85_mort"))), na.rm = TRUE)
  )

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
    total_2001    = sum(pm_2001_mort,    na.rm = TRUE),
    total_2041_45 = sum(pm_2041_45_mort, na.rm = TRUE),
    total_2041_85 = sum(pm_2041_85_mort, na.rm = TRUE),
    total_2091_45 = sum(pm_2091_45_mort, na.rm = TRUE),
    total_2091_85 = sum(pm_2091_85_mort, na.rm = TRUE)
  ) %>%
  # Reshape from wide (one column per scenario) to long (one row per scenario)
  pivot_longer(
    cols = everything(),
    names_to = "scenario",   # column names (e.g. "total_2041_45") become values in a new "scenario" column
    values_to = "pm_deaths"  # summed deaths become values in "pm_deaths"
  ) %>%
  # Strip the "total_" prefix so "total_2041_45" becomes "2041_45"
  mutate(scenario = str_remove(scenario, "^total_")) %>%
  # Split "2041_45" into two columns: year = "2041", rcp = "45"
  # fill = "right" handles the 2001 baseline which has no RCP (rcp will be NA)
  separate(scenario, into = c("year", "rcp"), sep = "_", fill = "right") %>%
  select(year, rcp, pm_deaths)

print("Total PM Mortality (USA), 2001, 2041, 2091")
print(usa_total_pm_mort)

usa_total_pm_mort_avg <- pop_pm_country_usa %>%
  summarise(
    total_2000s    = sum(pm_2000s_mort,    na.rm = TRUE),
    total_2040s_45 = sum(pm_2040s_45_mort, na.rm = TRUE),
    total_2040s_85 = sum(pm_2040s_85_mort, na.rm = TRUE),
    total_2090s_45 = sum(pm_2090s_45_mort, na.rm = TRUE),
    total_2090s_85 = sum(pm_2090s_85_mort, na.rm = TRUE)
  ) %>%
  # Reshape from wide (one column per scenario) to long (one row per scenario)
  pivot_longer(
    cols = everything(),
    names_to = "scenario",   # column names (e.g. "total_2040s_45") become values in a new "scenario" column
    values_to = "pm_deaths"  # summed deaths become values in "pm_deaths"
  ) %>%
  # Strip the "total_" prefix so "total_2040s_45" becomes "2040s_45"
  mutate(scenario = str_remove(scenario, "^total_")) %>%
  # Split "2040s_45" into two columns: year = "2040s", rcp = "45"
  # fill = "right" handles the 2000s baseline which has no RCP (rcp will be NA)
  separate(scenario, into = c("year", "rcp"), sep = "_", fill = "right") %>%
  select(year, rcp, pm_deaths)

print("Total PM Mortality Decade Averages (USA)")
print(usa_total_pm_mort_avg)

# Literature comparison
cat("\nLiterature Comparison (2000 baseline):\n")
cat("Ford et al 2018 Table 3: 138,000 deaths \n")
cat("Pierce et al 2017 Figure 13: ~165,000 deaths \n")
cat("Our estimate (2000s decade avg):", 
    round(usa_total_pm_mort_avg %>% filter(year == "2000s") %>% pull(pm_deaths)), "\n")

# our result is higher because of 1) higher RR from Orellano and 2) not using a threshold, per Orellano

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Fire PM Mortality Analysis (USA) ###########################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Calculate fire-attributable deaths for both methods
usa_fpm_mort_avg <- pop_pm_country_usa %>%
  summarise(
    fpm_2000s    = sum(fpm_2000s_mort,    na.rm = TRUE),
    fpm_2040s_45 = sum(fpm_2040s_45_mort, na.rm = TRUE),
    fpm_2040s_85 = sum(fpm_2040s_85_mort, na.rm = TRUE),
    fpm_2090s_45 = sum(fpm_2090s_45_mort, na.rm = TRUE),
    fpm_2090s_85 = sum(fpm_2090s_85_mort, na.rm = TRUE)
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "scenario",    # column names (e.g. "fpm_2040s_45") become values in "scenario"
    values_to = "fpm_deaths"  # summed deaths become values in "fpm_deaths"
  ) %>%
  # Strip "fpm_" prefix so "fpm_2040s_45" becomes "2040s_45"
  mutate(scenario = str_remove(scenario, "^fpm_")) %>%
  # Split "2040s_45" into year = "2040s", rcp = "45"; 2000s baseline gets rcp = NA
  separate(scenario, into = c("year", "rcp"), sep = "_", fill = "right") %>%
  select(year, rcp, fpm_deaths)

print("\nFire PM Mortality (USA):")
print(usa_fpm_mort_avg)

# Ford validation
cat("\n=== Ford et al 2018 Validation (2001 baseline) ===\n")
cat("Ford et al reported: 17,000 fire-attributable deaths\n")
cat("Our estimate (2000s decade avg):", round(usa_fpm_mort_avg %>% filter(year == "2000s") %>% pull(fpm_deaths)), "\n")
# our result is higher because of 1) higher RR from Orellano and 2) not using a threshold, per Orellano

# Fire as percentage of total mortality
usa_fire_pct <- usa_fpm_mort_avg %>%
  left_join(usa_total_pm_mort_avg, by = c("year", "rcp")) %>%
  mutate(fpm_pm_pct = (fpm_deaths / pm_deaths) * 100)

print("Fire as % of Total PM Mortality (Method 1):")
print(usa_fire_pct)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Temporal & RCP Trends ######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Create summary tables for cleaner output
temporal_trends <- usa_fpm_mort_avg %>%
  arrange(year, rcp)

cat("\n=== Fire PM Mortality Trends ===\n")
print(temporal_trends)

# Calculate percent changes from 2000 baseline
baseline <- usa_fpm_mort_avg %>% filter(year == "2000s") %>% pull(fpm_deaths)

temporal_pct_change <- temporal_trends %>%
  mutate(
    pct_change = (fpm_deaths / baseline - 1) * 100
  )

cat("\n=== Percent Change from 2000 Baseline ===\n")
print(temporal_pct_change)

# RCP comparisons
rcp_comparison <- temporal_trends %>%
  filter(year %in% c("2040s", "2090s")) %>%
  group_by(year) %>%
  summarise(
    rcp45        = fpm_deaths[rcp == "45"],
    rcp85        = fpm_deaths[rcp == "85"],
    rcp_diff_pct = (rcp85 / rcp45 - 1) * 100
  )

cat("\n=== RCP 8.5 vs RCP 4.5 Comparison ===\n")
print(rcp_comparison)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Comparison with Pierce and Ford ######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Create comparison object for Ford et al and Pierce et al benchmarking
usa_lit_comp <- usa_fire_pct %>% 
  select(year, rcp, fpm_deaths, pm_deaths, fpm_pm_pct)

# Literature results from Ford et al 2018, Pierce et al 2017, Qiu
# Create comparison object for Ford et al and Pierce et al benchmarking
lit_results <- tibble(
  year = c("2000s", "2010s", "2040s", "2040s", "2046-2055", "2090s", "2090s"),
  rcp  = c(NA,      NA,      "45",    "85", "2-4.5",    "45",    "85"),
  # Ford et al 2018 results
  fpm_ford              = c(17000,  NA,     42000,  32000, NA,  32000,  44000),
  pm_ford               = c(138000, NA,     114000, 105000, NA, 75000,  88000),
  fpm_pm_pct_ford       = c(0.12,   NA,     0.37,   0.30, NA,   0.43,   0.50),
  # Pierce et al 2017 results
  fpm_pierce            = c(20000,  NA,     47000,  NA, NA,     65000,  NA),
  pm_pierce             = c(165000, NA,     145000, NA, NA,     125000, NA),
  fpm_pm_pct_pierce     = c(0.12,   NA,     0.32,   NA, NA,     0.52,   NA),
  # Qiu et al 2024 EarthArXiv results
  fpm_qiu2024_EarthArXiv = c(NA,    15800,  NA, NA, NA,     NA,     NA),
  # Qiu et al 2025 Nature results
  fpm_qiu2025_nature     = c(NA,    41380,  NA, NA, 57593,     NA,     NA)
)

# Join on year and rcp 
usa_lit_comp <- usa_fire_pct %>%
  select(year, rcp, fpm_deaths, pm_deaths, fpm_pm_pct) %>%
  full_join(lit_results, by = c("year", "rcp"))

usa_lit_comp <- usa_lit_comp %>% # custom row order
  slice(1, 6, 2, 3, 7, 4, 5)

print("\n=== USA vs Literature Comparison ===")
print(usa_lit_comp)

# Prepare data for plotting
plot_data <- usa_lit_comp %>%
  mutate(
    # Create x-axis label: baseline stays "2000s", future scenarios get "2040s\nRCP45" etc.
    scenario = case_when(
      year == "2000s" ~ "2000s",
      year == "2010s" ~ "2010s",
      TRUE ~ paste0(year, "\nRCP", rcp)
    ),
    # Set display order of scenarios on x-axis
    scenario = factor(scenario, levels = c("2000s", "2010s", "2040s\nRCP45", "2040s\nRCP85",
                                           "2046-2055\nRCP2-4.5",
                                           "2090s\nRCP45", "2090s\nRCP85"))
  ) %>%
  # Keep only columns needed for plotting (drop rcp, year, pct columns)
  select(scenario, fpm_deaths, pm_deaths, fpm_ford, pm_ford,
         fpm_pierce, pm_pierce, fpm_qiu2024_EarthArXiv, fpm_qiu2025_nature) %>%
  # Reshape from wide (one column per source/metric) to long (one row per source/metric)
  pivot_longer(
    cols = -scenario,       # pivot all columns except scenario
    names_to = "variable",  # column names (e.g. "fpm_ford") become values in "variable"
    values_to = "deaths"    # death counts become values in "deaths"
  ) %>%
  mutate(
    # Classify each row as Fire PM or Total PM based on column name prefix
    metric = if_else(str_detect(variable, "^fpm"), "Fire PM", "Total PM"),
    # Classify each row by data source based on column name suffix
    source = case_when(
      str_detect(variable, "ford")               ~ "Ford et al 2018",
      str_detect(variable, "pierce")             ~ "Pierce et al 2017",
      str_detect(variable, "qiu2024_EarthArXiv") ~ "Qiu et al 2024 EarthArXiv",
      str_detect(variable, "qiu2025_nature")     ~ "Qiu et al 2025 Nature",
      TRUE                                       ~ "Our Estimate"
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
    values = c("Our Estimate"              = "#2E86AB", 
               "Ford et al 2018"           = "#A23B72", 
               "Pierce et al 2017"         = "#F18F01",
               "Qiu et al 2024 EarthArXiv" = "#44BBA4",
               "Qiu et al 2025 Nature"     = "#E94F37")
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
    legend.text = element_text(size = 8),
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )


# Display plots
print(p1)
print(p2)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Mortality by Country Analysis ##############################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Aggregate mortality by country (excluding no-fire scenarios)
country_mortality <- pop_pm_country %>%
  filter(!is.na(country_code_iso3) & country_code_iso3 != "-99") %>%
  group_by(country_code_iso3, country_name.x) %>%
  summarise(
    total_pop_2001            = sum(pop_tot_2001,       na.rm = TRUE),
    pm_2000s_mort        = sum(pm_2000s_mort,      na.rm = TRUE),
    pm_2040s_45_mort     = sum(pm_2040s_45_mort,   na.rm = TRUE),
    pm_2040s_85_mort     = sum(pm_2040s_85_mort,   na.rm = TRUE),
    pm_2090s_45_mort     = sum(pm_2090s_45_mort,   na.rm = TRUE),
    pm_2090s_85_mort     = sum(pm_2090s_85_mort,   na.rm = TRUE),
    fpm_2000s_mort       = sum(fpm_2000s_mort,     na.rm = TRUE),
    fpm_2040s_45_mort    = sum(fpm_2040s_45_mort,  na.rm = TRUE),
    fpm_2040s_85_mort    = sum(fpm_2040s_85_mort,  na.rm = TRUE),
    fpm_2090s_45_mort    = sum(fpm_2090s_45_mort,  na.rm = TRUE),
    fpm_2090s_85_mort    = sum(fpm_2090s_85_mort,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(total_pop_2001 > 0) %>%
  mutate(fire_pct_2000s = (fpm_2000s_mort / pm_2000s_mort) * 100) %>%
  arrange(desc(pm_2000s_mort))

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
    gmt_chg = c(gmt_baseline, gmt_2050_45, gmt_2050_85, gmt_2100_45, gmt_2100_85),
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


print(global_totals_all)

cat("\n=== Temperature-Mortality Relationship ===\n")
cat("At +", round(gmt_baseline, 2), "°C (2000 baseline):", round(global_totals_all$fpm_deaths[1]), "deaths\n")
cat("At +", round(gmt_2050_45, 2), "°C (2050 RCP4.5):", round(global_totals_all$fpm_deaths[2]), "deaths\n")
cat("At +", round(gmt_2050_85, 2), "°C (2050 RCP8.5):", round(global_totals_all$fpm_deaths[3]), "deaths\n")
cat("At +", round(gmt_2100_45, 2), "°C (2100 RCP4.5):", round(global_totals_all$fpm_deaths[4]), "deaths\n")
cat("At +", round(gmt_2100_85, 2), "°C (2100 RCP8.5):", round(global_totals_all$fpm_deaths[5]), "deaths\n")
cat("\nIncrease from +", round(gmt_baseline, 2), "°C to +", round(gmt_2100_85, 2), "°C:",
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
  scale_x_continuous(breaks = c(gmt_baseline, gmt_2050_45, gmt_2050_85, gmt_2100_45, gmt_2100_85)) +
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


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ USA Death and Temperature Change ##########################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Add additional columns to existing usa_fire_pct object (matching global_totals_all structure)
usa_totals_all <- usa_fire_pct %>%
  mutate(
    all_death_usa_2009 = 2437163,  # Total USA deaths in 2009
    gmt_chg = c(gmt_baseline, gmt_2050_45, gmt_2050_85, gmt_2100_45, gmt_2100_85),
    death_pct_chg = (fpm_deaths / all_death_usa_2009)  # Fire PM deaths as % of all USA deaths (decimal form)
  )

# Create plot of fire PM deaths vs global mean temperature change
p_temp_usa <- usa_totals_all %>%
  ggplot(aes(x = gmt_chg, y = fpm_deaths)) +
  geom_smooth(method = "lm", se = TRUE, color = "#2E86AB", fill = "#2E86AB", alpha = 0.2) +
  geom_point(size = 3, color = "#2E86AB") +
  scale_y_continuous(labels = scales::comma) +
  scale_x_continuous(breaks = c(gmt_baseline, gmt_2050_45, gmt_2050_85, gmt_2100_45, gmt_2100_85)) +
  labs(
    title = "Fire PM2.5 Mortality vs Global Mean Temperature Change",
    subtitle = "USA Estimates",
    x = "Global Mean Temperature Change (°C)",
    y = "Annual Fire PM Deaths"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# Display plot
print(p_temp_usa)


# Calculate and display linear regression statistics
temp_model_usa <- lm(fpm_deaths ~ gmt_chg, data = usa_totals_all)
cat("\n=== USA Linear Regression: Fire PM Deaths ~ Temperature Change ===\n")
cat("Equation: Fire PM Deaths = ", round(coef(temp_model_usa)[1]), " + ", 
    round(coef(temp_model_usa)[2]), " × GMT Change\n", sep = "")
cat("R-squared:", round(summary(temp_model_usa)$r.squared, 3), "\n")
cat("P-value:", format.pval(summary(temp_model_usa)$coefficients[2,4], digits = 3), "\n")
cat("\nInterpretation: Each 1°C increase in global mean temperature is associated with\n")
cat("approximately", round(coef(temp_model_usa)[2]), "additional fire PM deaths per year in the USA.\n")


# Optional: Save plot
# ggsave("fire_pm_mortality_vs_temperature_usa.png", p_temp_usa, width = 8, height = 6, dpi = 300)

# Summary statistics
cat("\n=== USA Fire PM Deaths by Temperature Scenario ===\n")
usa_totals_all %>%
  select(year, rcp, gmt_chg, fpm_deaths, death_pct_chg) %>%
  mutate(
    fpm_deaths = round(fpm_deaths),
    death_pct_chg = round(death_pct_chg, 3)
  ) %>%
  print()

cat("\n=== USA Temperature-Mortality Relationship ===\n")
cat("At +", round(gmt_baseline, 2), "°C (2000 baseline):", round(usa_totals_all$fpm_deaths[1]), "deaths\n")
cat("At +", round(gmt_2050_45, 2), "°C (2050 RCP4.5):", round(usa_totals_all$fpm_deaths[2]), "deaths\n")
cat("At +", round(gmt_2050_85, 2), "°C (2050 RCP8.5):", round(usa_totals_all$fpm_deaths[3]), "deaths\n")
cat("At +", round(gmt_2100_45, 2), "°C (2100 RCP4.5):", round(usa_totals_all$fpm_deaths[4]), "deaths\n")
cat("At +", round(gmt_2100_85, 2), "°C (2100 RCP8.5):", round(usa_totals_all$fpm_deaths[5]), "deaths\n")
cat("\nIncrease from +", round(gmt_baseline, 2), "°C to +", round(gmt_2100_85, 2), "°C:",
    round(usa_totals_all$fpm_deaths[5] - usa_totals_all$fpm_deaths[1]), 
    "deaths (", 
    round((usa_totals_all$fpm_deaths[5] / usa_totals_all$fpm_deaths[1] - 1) * 100, 1), 
    "% increase)\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ USA All-cause death % chg due to FPM and Temperature Change #################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Calculate linear regression model first
temp_pct_model_usa <- lm(death_pct_chg ~ gmt_chg, data = usa_totals_all)

# Display regression statistics
cat("\n=== USA Linear Regression: Fire PM % of All Deaths ~ Temperature Change ===\n")
cat("Equation: Fire PM % = ", round(coef(temp_pct_model_usa)[1], 5), " + ", 
    round(coef(temp_pct_model_usa)[2], 5), " × GMT Change\n", sep = "")
cat("R-squared:", round(summary(temp_pct_model_usa)$r.squared, 3), "\n")
cat("P-value:", format.pval(summary(temp_pct_model_usa)$coefficients[2,4], digits = 3), "\n")
cat("\nInterpretation: Each 1°C increase in global mean temperature is associated with\n")
cat("an increase of", round(coef(temp_pct_model_usa)[2], 5), "percentage points in the share of\n")
cat("USA all-cause deaths attributable to fire PM.\n")
cat("\nFor example, if fire PM deaths are currently 1% of all USA deaths, a 1°C warming would\n")
cat("increase this to approximately", round(1 + coef(temp_pct_model_usa)[2], 5), "% of all deaths.\n")

# Create plot of fire PM deaths as % of all deaths vs global mean temperature change
# Add scenario labels to the data
usa_totals_all <- usa_totals_all %>%
  mutate(
    scenario_label = case_when(
      year == "2000" ~ "2000",
      year == "2050" & rcp == "45" ~ "RCP4.5 2050",
      year == "2050" & rcp == "85" ~ "RCP8.5 2050",
      year == "2100" & rcp == "45" ~ "RCP4.5 2100",
      year == "2100" & rcp == "85" ~ "RCP8.5 2100"
    )
  )


p_temp_pct_usa <- usa_totals_all %>%
  ggplot(aes(x = gmt_chg, y = death_pct_chg)) +
  geom_smooth(method = "lm", se = FALSE, color = "#2E86AB") +
  geom_text(aes(label = scenario_label), # Add RCP labels for each point
            hjust = .5, vjust = -1, 
            size = 3.5, color = "#A23B72") +
  geom_point(size = 3, color = "#2E86AB") +
  scale_x_continuous(breaks = c(gmt_baseline, gmt_2050_45, gmt_2050_85, gmt_2100_45, gmt_2100_85)) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "USA All-cause death %Δ due to FPM vs Temperature Change",
    subtitle = paste0("USA All-cause death %Δ = ", round(coef(temp_pct_model_usa)[1], 5), " + ", 
                      round(coef(temp_pct_model_usa)[2], 5), " × GMT Change"),
    x = "Global Mean Temperature Change (°C)",
    y = "USA All-cause death %Δ due to FPM"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# Display plot
print(p_temp_pct_usa)

# Optional: Save plot
# ggsave("fire_pm_pct_vs_temperature_usa.png", p_temp_pct_usa, width = 8, height = 6, dpi = 300)

# Summary statistics
cat("\n=== USA Fire PM Deaths as % of All Deaths by Temperature Scenario ===\n")
usa_totals_all %>%
  select(year, rcp, gmt_chg, death_pct_chg) %>%
  mutate(
    death_pct_chg = round(death_pct_chg * 100, 3)
  ) %>%
  print()

cat("\n=== Temperature-Mortality Percentage Relationship ===\n")
cat("At +", round(gmt_baseline, 2), "°C (2000 baseline):", round(global_totals_all$death_pct_chg[1] * 100, 3), "% of all deaths\n")
cat("At +", round(gmt_2050_45, 2), "°C (2050 RCP4.5):", round(global_totals_all$death_pct_chg[2] * 100, 3), "% of all deaths\n")
cat("At +", round(gmt_2050_85, 2), "°C (2050 RCP8.5):", round(global_totals_all$death_pct_chg[3] * 100, 3), "% of all deaths\n")
cat("At +", round(gmt_2100_45, 2), "°C (2100 RCP4.5):", round(global_totals_all$death_pct_chg[4] * 100, 3), "% of all deaths\n")
cat("At +", round(gmt_2100_85, 2), "°C (2100 RCP8.5):", round(global_totals_all$death_pct_chg[5] * 100, 3), "% of all deaths\n")
cat("\nIncrease from +", round(gmt_baseline, 2), "°C to +", round(gmt_2100_85, 2), "°C:",
    round((global_totals_all$death_pct_chg[5] - global_totals_all$death_pct_chg[1]) * 100, 3), 
    "percentage points (", 
    round((global_totals_all$death_pct_chg[5] / global_totals_all$death_pct_chg[1] - 1) * 100, 1), 
    "% increase)\n")



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ SAVE OUTPUTS #################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ggsave(here("images/death", "total_pm_mortality_usa.png"),  p1, width = 8, height = 6, dpi = 300)
ggsave(here("images/death", "fire_pm_mortality_usa.png"),   p2, width = 8, height = 6, dpi = 300)

write_csv(country_mortality, here("output", "country_mortality.csv")) # save death all countries / scenarios
write_csv(usa_fire_pct, here("output", "usa_fire_pct.csv")) # save death USA for scenarios 

#### THE END



