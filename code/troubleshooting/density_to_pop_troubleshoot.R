

cat("Year index used:", year_2009_index, "\n")
cat("Year at that index:", years[year_2009_index], "\n")
cat("Mean pop density:", mean(pop_df$pop_density[pop_df$pop_density > 0]), "\n")

# Total land area being counted
total_land_area <- sum(pop_df$cell_area_km2[pop_df$pop_density > 0])
cat("Total land area (million km²):", total_land_area / 1e6, "\n")
cat("Expected ~150 million km² for Earth's land\n")

# Population should only be on land
pop_total_land_only <- sum(pop_df$pop_total[pop_df$pop_density > 0])
cat("Population (land cells only, millions):", pop_total_land_only / 1e6, "\n")


# The metadata showed a LANDMASK variable - extract it
landmask <- ncvar_get(pop_ssp1, "LANDMASK")
sum(landmask)  # Count of land cells


# Manual verification of the conversion
mean_density <- 51.62816  # people/km²
mean_area <- mean(pop_df$cell_area_km2[pop_df$pop_density > 0])
n_populated_cells <- sum(pop_df$pop_density > 0)

cat("Mean density (people/km²):", mean_density, "\n")
cat("Mean cell area (km²):", mean_area, "\n")
cat("Number of populated cells:", n_populated_cells, "\n")

# Expected total if formula is correct:
expected_total <- mean_density * mean_area * n_populated_cells
cat("Expected population (millions):", expected_total / 1e6, "\n")
cat("Actual population (millions):", pop_total_land_only / 1e6, "\n")

# Also check: does the mean of pop_total match mean_density × mean_area?
cat("Mean pop_total:", mean(pop_df$pop_total[pop_df$pop_density > 0]), "\n")
cat("Should equal:", mean_density * mean_area, "\n")



# Test the conversion on a few cells
test_cells <- pop_df %>%
  filter(pop_density > 0) %>%
  slice(1:5) %>%
  mutate(
    # Your current formula (lat in radians)
    area_current = (0.5 * 111) * (0.5 * 111 * cos(lat * pi/180)),
    
    # Alternative: maybe cos already expects degrees?
    area_alt = (0.5 * 111) * (0.5 * 111 * cos(lat)),
    
    # Or check if the data uses different Earth radius
    # Standard: 111 km/degree at equator
    # Some use 111.32 km
    area_precise = (0.5 * 111.32) * (0.5 * 110.57 * cos(lat * pi/180)),
    
    pop_current = pop_density * area_current,
    pop_alt = pop_density * area_alt
  ) %>%
  select(lat, pop_density, area_current, area_alt, area_precise, pop_current, pop_alt)

print(test_cells)



# The correct check: sum, not mean
sum_density_times_area <- sum(pop_df$pop_density[pop_df$pop_density > 0] * 
                                pop_df$cell_area_km2[pop_df$pop_density > 0])

sum_pop_total <- sum(pop_df$pop_total[pop_df$pop_density > 0])

cat("Sum(density × area):", sum_density_times_area / 1e6, "million\n")
cat("Sum(pop_total):", sum_pop_total / 1e6, "million\n")
cat("Should be identical:", sum_density_times_area == sum_pop_total, "\n")

# If those match, the formula is correct
# The mean mismatch is just because:
# mean(A × B) ≠ mean(A) × mean(B) when A and B vary together


# All cells (including ocean)
total_surface <- sum(pop_df$cell_area_km2)
cat("Total surface area (million km²):", total_surface / 1e6, "\n")
cat("Expected ~510 million km² for Earth's total surface\n")
