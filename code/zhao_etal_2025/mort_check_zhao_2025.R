# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Process Zhao et al 2025 PM data ##
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

## baseline simulation (varying pop)

# mortality without fire under SSP5-8.5
mort_SSP585 <- read.csv(
  here("input/Zhao_etal_2025/gridded_output", "SSP585_mortality.csv"),
  stringsAsFactors = FALSE
)

# mortality without fire under SSP2-4.5
mort_SSP245 <- read.csv(
  here("input/Zhao_etal_2025/gridded_output", "SSP245_mortality.csv"),
  stringsAsFactors = FALSE
)

head(mort_SSP245)
head(mort_SSP585)

## EM simulation (constant pop)

# mortality without fire under SSP5-8.5
mort_SSP585_EM <- read.csv(
  here("input/Zhao_etal_2025/gridded_output_EM", "SSP585-EM_mortality.csv"),
  stringsAsFactors = FALSE
)

# mortality without fire under SSP2-4.5
mort_SSP245_EM <- read.csv(
  here("input/Zhao_etal_2025/gridded_output_EM", "SSP245-EM_mortality.csv"),
  stringsAsFactors = FALSE
)

head(mort_SSP585_EM)
head(mort_SSP245_EM)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Global total mortality #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mort_tot_585 <- sum(mort_SSP585$mortality_mean)
cat("Total deaths under SSP5-8.5:", mort_tot_585, "thousand", "\n")

mort_tot_245 <- sum(mort_SSP245$mortality_mean)
cat("Total deaths under SSP2-4.5:", mort_tot_245, "thousand", "\n")

mort_tot_585_EM <- sum(mort_SSP585_EM$mortality_mean)
cat("Total deaths under SSP5-8.5 (EM Simulation):", mort_tot_585_EM, "thousand", "\n")

mort_tot_245 <- sum(mort_SSP245_EM$mortality_mean)
cat("Total deaths under SSP2-4.5 (EM Simulation):", mort_tot_245, "thousand", "\n")



### THE END 

