# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###################Population grid cell data ##################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())

library(ncdf4)
library(here)
library(dplyr)
library(readr)
library(data.table)
library(tibble)
library(rnaturalearth)
library(sf)
library(tidyr)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
################################# read in data #################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# tried:
#SSP3_for_RCP85_2006-2100_population_density_c160701.nc/different/8951715907(2006-2010)
#SSP1_for_RCP45_2006-2100_population_density_c160701.nc/different/8951715907 (2006-2010)
#clmforc.Li_2018_SSP3_CMIP6_hdm_0.5x0.5_AVHRR_simyr1850-2100_c181205.nc/different/7436984063
#clmforc.Li_2018_SSP1_CMIP6_hdm_0.5x0.5_AVHRR_simyr1850-2100_c181205.nc/different/7436984063
#clmforc.Li_2017_HYDEv3.2_CMIP6_hdm_0.5x0.5_AVHRR_simyr1850-2016_c170828.nc/same/1342428355
#clmforc.Li_2012_hdm_0.5x0.5_AVHRR_simyr1850-2010_c130401.nc/different/8951715907

# decide to use: clmforc.Li_2018_SSP1_CMIP6_hdm_0.5x0.5_AVHRR_simyr1850-2100_c181205.nc
clm_nc <- nc_open(here("input/google_drive/","clmforc.Li_2018_SSP1_CMIP6_hdm_0.5x0.5_AVHRR_simyr1850-2100_c181205.nc"))
print(clm_nc)
area_nc <- nc_open(here("input/google_drive/","gridcell_area_0.5deg.nc"))
print(area_nc)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############## Extract Population data from NetCDF ###################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
lon <- ncvar_get(clm_nc, "LONGXY")
lat <- ncvar_get(clm_nc, "LATIXY")

grid_df <- tibble(
  lon = as.vector(lon),
  lat = as.vector(lat)
)

years <- ncvar_get(clm_nc, "year")

# select years
#target_years <- 2001:2010
target_years <- 2001:2010
year_indices <- match(target_years, years) 
print(data.frame(target_years, year_indices, years_at_index = years[year_indices]))

stopifnot(!any(is.na(year_indices)))
stopifnot(all(years[year_indices] == target_years))


# ================== check if all years have same data =====================
# Use year 2009 to find a grid cell with population > 0
idx2009 <- year_indices[which(target_years == 2009)]

slice2009 <- ncvar_get(clm_nc, "hdm", start = c(1,1,idx2009), count = c(-1,-1,1))

# Identify all positions where population density > 0 (lon index, lat index)
pos <- which(slice2009 > 0, arr.ind = TRUE)

cat("Non-zero cells in 2009:", nrow(pos), "\n")

# Randomly select one non-zero grid cell (guaranteed land / populated)
set.seed(123)
pick <- pos[sample(nrow(pos), 1), ]
li <- pick[1]  # longitude index
la <- pick[2]  # latitude index

cat("Picked indices:", li, la, " lon/lat = ", lon[li], lat[la], "\n")

# Extract the time series (2001–2010) for this grid cell
vals <- sapply(year_indices, function(ti)
  ncvar_get(clm_nc, "hdm", start = c(li, la, ti), count = c(1,1,1))
)

print(data.frame(pop_year = target_years, hdm = as.numeric(vals)))


# =============== continue extract Population data from NetCDF =================
#loop to get pop density
pop_list <- vector("list", length(year_indices))

for (i in seq_along(year_indices)) {
  
  yr_index <- year_indices[i]
  yr_value <- years[yr_index]
  
  cat("Reading year:", yr_value, "\n")
  
  pop_slice <- ncvar_get(
    clm_nc,
    "hdm",
    start = c(1, 1, yr_index),
    count = c(-1, -1, 1)
  )
  
  stopifnot(length(pop_slice) == nrow(grid_df))
  
  pop_list[[i]] <- grid_df %>%
    mutate(
      pop_density = as.vector(pop_slice),
      pop_year    = as.integer(yr_value)
    ) %>%
    select(lon, lat, pop_density, pop_year)
}

#create dataframe
pop_df <- bind_rows(pop_list)
print(pop_df)

# Result: Long-format dataframe with 4 columns: lon, lat, pop_density, year
# Each row = one grid cell's population density at its coordinates

# basic checks
year_check <- pop_df %>%
  group_by(pop_year) %>%
  summarise(
    cells = n(),
    non_zero_cells = sum(pop_density > 0, na.rm = TRUE),
    na_cells = sum(is.na(pop_density)),
    zero_cells = sum(pop_density == 0, na.rm = TRUE),
    
    total_density = sum(pop_density, na.rm = TRUE),
    
    share_nonzero = mean(pop_density > 0, na.rm = TRUE),
    share_zero = mean(pop_density == 0, na.rm = TRUE),
    
    mean_density = mean(pop_density, na.rm = TRUE),
    median_density = median(pop_density, na.rm = TRUE),
    max_density = max(pop_density, na.rm = TRUE),
    .groups = "drop"
  )

print(year_check, n = Inf)

#close nc
nc_close(clm_nc)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############## Extract cell area data from NetCDF ###################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# area variable is "area" with dims [lon, lat], units m2
area_m2 <- ncvar_get(area_nc, "area")
nc_close(area_nc)

# make sure area grid matches population grid
# check by comparing dimension lengths and coordinate values if available
stopifnot(all(dim(area_m2) == dim(lon)))
stopifnot(all(dim(area_m2) == dim(lat)))

# Convert to km2
area_km2 <- area_m2 / 1e6

# Build area dataframe (one row per cell)
area_df <- tibble(
  lon = as.vector(lon),
  lat = as.vector(lat),
  cell_area_km2 = as.vector(area_km2)
) %>%
  select(lon, lat, cell_area_km2)

# join area onto pop_df and compute pop total
pop_df2 <- pop_df %>%
  left_join(area_df, by = c("lon", "lat")) %>%
  mutate(pop_tot = pop_density * cell_area_km2)

# sanity check: any missing area?
cat("Missing cell_area_km2 after join:", sum(is.na(pop_df2$cell_area_km2)), "\n")

# global population totals by year (millions)
global_pop_check <- pop_df2 %>%
  group_by(pop_year) %>%
  summarise(
    global_pop_million = sum(pop_tot, na.rm = TRUE) / 1e6,
    .groups = "drop"
  )

print(global_pop_check, n = Inf)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############## wide output: one row per cell, years as columns ##################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 1) pop density wide
pop_dens_wide <- pop_df2 %>%
  select(lon, lat, pop_year, pop_density) %>%
  pivot_wider(
    names_from = pop_year,
    values_from = pop_density,
    names_prefix = "pop_dens_"
  )

# 2) pop total wide
pop_tot_wide <- pop_df2 %>%
  select(lon, lat, pop_year, pop_tot) %>%
  pivot_wider(
    names_from = pop_year,
    values_from = pop_tot,
    names_prefix = "pop_tot_"
  )

# 3) final wide table: lon/lat + area + all year columns
pop_df_rev <- area_df %>%
  left_join(pop_dens_wide, by = c("lon", "lat")) %>%
  left_join(pop_tot_wide, by = c("lon", "lat"))

print(pop_df_rev)

#check total polulation and check total USA pop
sum(pop_df_rev$pop_tot_2009)

#USA pop
sf_use_s2(FALSE)

usa <- ne_countries(scale = 50, returnclass = "sf") %>%
  filter(iso_a3 == "USA") %>%
  st_make_valid()

pop_pts <- pop_df2 %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

usa_idx <- lengths(st_intersects(pop_pts, usa)) > 0
pop_usa <- pop_pts[usa_idx, ]

usa_pop_all_years <- pop_usa %>%
  st_drop_geometry() %>%
  group_by(pop_year) %>%
  summarise(
    usa_population = sum(pop_tot, na.rm = TRUE),
    usa_pop_million = usa_population / 1e6,
    .groups = "drop"
  )

print(usa_pop_all_years, n = Inf)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############## output final data set ##################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write_csv(pop_df_rev, here("output", "pop_df_rev.csv"))
