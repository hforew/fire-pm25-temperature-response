# Remove all objects from the environment
rm(list = ls())


# packages

library(ncdf4)
library(here)
library(tidyverse)
library(haven) # for .sav file


list.files(here("input"))
list.files(here("input/landmask_area"))


# Open the NetCDF file
area <- nc_open(here("input/landmask_area", "cesm130_clm5_firemodule_area_f09x125.nc"))
print(area_09_125)

area_05_05 <- nc_open(here("input/landmask_area", "sfcarea4popgrid.nc"))
print(area_05_05)

land_frac <- nc_open(here("input/landmask_area", "cesm130_clm5_firemodule_landfrac_f09x125.nc"))
print(land_frac)

country_cell <- nc_open(here("input/landmask_area", "countries_gridded_0.1deg_v0.1.nc"))
print(country_cell)
iso_codes <- ncvar_get(country_cell, "iso")
print(iso_codes)

# pop data
pop_ssp1 <- nc_open(here("input", "SSP1_for_RCP45_2006-2100_population_density_c160701.nc"))
pop_ssp3 <- nc_open(here("input", "SSP3_for_RCP85_2006-2100_population_density_c160701.nc"))

years <- ncvar_get(pop_ssp1, "year")
min(years)
print(pop_ssp1)  # Full structure
years[1:10]      # First 10 year values

hdm_ssp1 <- ncvar_get(pop_ssp1, "hdm")
hdm_ssp3 <- ncvar_get(pop_ssp3, "hdm")  # Fixed: use ncvar_get, not nc_open

# Check 2006-2009 (indices 1-4)
identical(hdm_ssp1[,,1:4], hdm_ssp3[,,1:4])  # Should be TRUE (pop density is identical for historical period)

# Check 2010 onwards (index 5+)
identical(hdm_ssp1[,,5:95], hdm_ssp3[,,5:95])  # Should be FALSE (pop density is different for future period)


## other pop data
HYDEv3.2 <- nc_open(here("input/other", "clmforc.Li_2017_HYDEv3.2_CMIP6_hdm_0.5x0.5_AVHRR_simyr1850-2016_c170828.nc"))
print(HYDEv3.2)

SSP1_CMIP6_hdm <- nc_open(here("input/other", "clmforc.Li_2018_SSP1_CMIP6_hdm_0.5x0.5_AVHRR_simyr1850-2100_c181205.nc"))
print(SSP1_CMIP6_hdm)
### FILE CONTENTS -- SSP1_CMIP6_hdm
# A 0.5° × 0.5° global grid of population density (people per km² of land), for every year from 1850 to 2100.
# Historical (1850–2016): from HYDE v3.2
# Future (2017–2100): from SSP1 scenario
# Annual resolution
# Land-only (oceans zeroed out using a CLM land mask)

SSP3_CMIP6_hdm <- nc_open(here("input/other", "clmforc.Li_2018_SSP3_CMIP6_hdm_0.5x0.5_AVHRR_simyr1850-2100_c181205.nc"))
print(SSP3_CMIP6_hdm)


years <- ncvar_get(SSP1_CMIP6_hdm, "year")
min(years)
years[1:10]      # First 10 year values

hdm_HYDE <- ncvar_get(HYDEv3.2, "hdm")
hdm_SSP3_CMIP6 <- ncvar_get(SSP3_CMIP6_hdm, "hdm")  # Fixed: use ncvar_get, not nc_open

# Check 2006-2009 (indices 1-4)
identical(hdm_HYDE[,,1:4], hdm_SSP3_CMIP6[,,1:4])  # hyde not same of SSP3_CMIP6

# Check 2010 onwards (index 5+)
identical(hdm_HYDE[,,5:95], hdm_SSP3_CMIP6[,,5:95])  # Should be FALSE (pop density is different for future period)



