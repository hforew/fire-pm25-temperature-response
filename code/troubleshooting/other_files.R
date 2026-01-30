# Remove all objects from the environment
rm(list = ls())


# packages

library(ncdf4)
library(here)
library(tidyverse)


# Open the NetCDF file
area <- nc_open(here("input/landmask_area", "cesm130_clm5_firemodule_area_f09x125.nc"))
print(area)


HYDEv3.2 <- nc_open(here("input/other", "clmforc.Li_2017_HYDEv3.2_CMIP6_hdm_0.5x0.5_AVHRR_simyr1850-2016_c170828.nc"))
print(HYDEv3.2)

SSP1_CMIP6_hdm <- nc_open(here("input/other", "clmforc.Li_2018_SSP1_CMIP6_hdm_0.5x0.5_AVHRR_simyr1850-2100_c181205.nc"))
print(SSP1_CMIP6_hdm)

SSP3_CMIP6_hdm <- nc_open(here("input/other", "clmforc.Li_2018_SSP3_CMIP6_hdm_0.5x0.5_AVHRR_simyr1850-2100_c181205.nc"))
print(SSP3_CMIP6_hdm)





