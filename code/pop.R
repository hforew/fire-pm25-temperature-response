
library(ncdf4)
library(here)


pop_ssp1 <- nc_open(here("input", "SSP1_for_RCP45_2006-2100_population_density_c160701.nc"))
print(pop_ssp1)

# Extract the population density variable for historical period (2006-2009)
pop_density <- ncvar_get(pop_ssp1, "hdm")  # dim: [720, 360, 95]
landmask <- ncvar_get(pop_ssp1, "LANDMASK")  # dim: [720, 360]

# Get first 4 time steps (2006-2009)
historical_years <- 1:4

# Calculate statistics for historical period
historical_pop <- pop_density[, , historical_years]

# Check how many land cells have population > 0
land_cells <- sum(landmask == 1)  # Total land cells
cells_with_pop <- apply(historical_pop, c(1, 2), function(x) {
  any(x > 0, na.rm = TRUE)  # TRUE if ANY year has pop > 0
})

# Mask by land only
cells_with_pop_on_land <- sum(cells_with_pop & landmask == 1)

# Calculate percentage
pct_coverage <- (cells_with_pop_on_land / land_cells) * 100

cat("Total land cells:", land_cells, "\n")
cat("Land cells with population (2006-2009):", cells_with_pop_on_land, "\n")
cat("Percentage coverage:", round(pct_coverage, 2), "%\n")

# Additional check: average across historical period
avg_historical_pop <- apply(historical_pop, c(1, 2), mean, na.rm = TRUE)

# Summary of non-zero population values
summary(avg_historical_pop[landmask == 1 & avg_historical_pop > 0])