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


# pm2.5 concentration without fire under SSP5-8.5
pm_nf_SSP585 <- read.csv(
  here("input/Zhao_etal_2025/gridded_output", "SSP585_without fire.csv"),
  stringsAsFactors = FALSE
)

head(data)

# pm2.5 concentration without fire under SSP2-4.5
pm_nf_SSP245

# mortality without fire under SSP5-8.5
mort_SSP585 <- read.csv(
  here("input/Zhao_etal_2025/gridded_output", "SSP585_mortality.csv"),
  stringsAsFactors = FALSE
)

head(mort_SSP585)

# mortality without fire under SSP2-4.5
mort_SSP245 <- read.csv(
  here("input/Zhao_etal_2025/gridded_output", "SSP245_mortality.csv"),
  stringsAsFactors = FALSE
)

head(mort_SSP585)


mort_tot_585 <- sum(mort_SSP585$mortality_mean)
cat("Total deaths under SSP5-8.5:", mort_tot_585, "thousand", "\n")

mort_tot_245 <- sum(mort_SSP245$mortality_mean)
cat("Total deaths under SSP2-4.5:", mort_tot_245, "thousand", "\n")

