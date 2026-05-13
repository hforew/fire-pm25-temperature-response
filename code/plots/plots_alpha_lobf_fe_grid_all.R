# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## PLOTS OF GRID-CELL ALPHA COEFFICIENTS - FE regression (grid-level)
##
## Three outputs:
##   1. Histogram of alpha_i2 across all ~65k grid cells
##   2. World raster maps of alpha_i1 (intercept) and alpha_i2 (slope)
##   3. Line-of-best-fit scatter plots for two target locations:
##        Denver, CO  (39.74°N, 104.99°W)
##        Los Angeles, CA (34.05°N, 118.24°W)
##
## Coefficients read directly from fpm_gmt_regression_coefs_FE_grid.csv.
## Raw observations for LoBF are rebuilt from the grid input + GMT scalars
## (avoids re-running the full 65k-cell regression).
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Packages #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(here)
library(tidyverse)
library(ggplot2)
library(ggrepel)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Import: coefficients #############################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Grid-cell FE regression coefficients produced by fpm_gmt_regress_FE_grid_all.R.
# One row per 0.5° grid cell (~65k cells). Key columns:
#   row_id, lon, lat      -- cell identifier and centre coordinates
#   estimate_alpha_i1     -- OLS intercept (CLASSIC baseline; extrapolated, not interpretable)
#   estimate_alpha_i2     -- OLS slope: change in fire PM2.5 (µg/m^3) per 1°C GMT

reg_coefs <- read_csv(here("output", "fpm_gmt_regression_coefs_FE_grid.csv"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Histogram of alpha_i2 ####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

hist_alpha_i2 <- ggplot(reg_coefs, aes(x = estimate_alpha_i2)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(
    title    = "Distribution of alpha_i2 across grid cells",
    subtitle = "Change in fire PM2.5 concentration (µg/m³) per 1°C GMT increase",
    x        = expression(alpha[i2] ~ "(µg/m³ per °C)"),
    y        = "Number of grid cells"
  ) +
  theme_minimal()

hist_alpha_i2
ggsave(here("images/regression_alpha/alpha_FE_grid_all", "hist_alpha_i2_fe_grid.png"),
       plot = hist_alpha_i2, width = 8, height = 5, dpi = 300)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Descriptive stats: alpha_i2 distribution ################################
##
## Run before defining cut_alpha() breaks to see where the data actually falls.
## Key outputs: percentiles, share negative, share above candidate break thresholds.
## Use these to judge whether break values capture meaningful variation.
##
## Stats computed for: (1) all global land cells, (2) CONUS cells only.
## cut_alpha() applies to both maps, so check both distributions before setting breaks.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

thresholds <- c(0, 0.5, 1, 1.5, 2)

print_alpha_stats <- function(x, label) {
  cat("\n--- alpha_i2 descriptive stats:", label, "(n =", sum(!is.na(x)), "cells) ---\n")
  cat("Min:  ", round(min(x,  na.rm = TRUE), 3), "\n")
  cat("Max:  ", round(max(x,  na.rm = TRUE), 3), "\n")
  cat("Mean: ", round(mean(x, na.rm = TRUE), 3), "\n")
  cat("SD:   ", round(sd(x,   na.rm = TRUE), 3), "\n")
  pctiles <- quantile(x, probs = c(0.02, 0.05, 0.10, 0.25, 0.50,
                                    0.75, 0.90, 0.95, 0.98), na.rm = TRUE)
  cat("\nPercentiles:\n")
  print(round(pctiles, 3))
  cat("\nShare of cells below each candidate break (%):\n")
  for (thr in thresholds) {
    cat("  <", thr, ":", round(100 * mean(x < thr, na.rm = TRUE), 1), "%\n")
  }
  cat("  > 2 :", round(100 * mean(x > 2, na.rm = TRUE), 1), "%\n")
}

# Global: all land cells with a regression result
print_alpha_stats(reg_coefs$estimate_alpha_i2, "GLOBAL")

# CONUS: filtered here for stats and reused for the CONUS map below
reg_coefs_usa <- reg_coefs %>%
  filter(country_code_iso3 == "USA",  # keep only cells mapped to USA
         lon >= -125, lon <= -66,      # CONUS longitude bounds (excludes Alaska)
         lat >=   24, lat <=  50)      # CONUS latitude bounds (excludes Hawaii)

print_alpha_stats(reg_coefs_usa$estimate_alpha_i2, "CONUS")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Map color palette  #################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Exact color assigned to each band -- no gradient interpolation.
# cut() bins the continuous alpha value into a factor; scale_fill_manual()
# maps each factor level to its exact hex color. Both global and CONUS maps
# use the same band_labels and band_colors, so identical values always
# produce identical colors regardless of which cells are in the plot.
band_labels <- c("<0", "0-0.10", "0.10-0.25", "0.25-0.50", "0.50-1.00", ">1.00")

band_colors <- c(
  "<0"        = "#7FD67F",  # light green  -- warming reduces fire PM2.5
  "0-0.10"    = "#56C8FF",  # sky blue     -- low positive response
  "0.10-0.25" = "#C8A000",  # dark yellow  -- moderate-low
  "0.25-0.50" = "#CC5500",  # dark orange  -- moderate-high
  "0.50-1.00" = "#DD1111",  # red          -- high
  ">1.00"     = "#7A0000"   # dark red     -- extreme (Siberia, Congo, Amazon globally)
)

# Band boundaries: -Inf, 0, 0.1, 0.25, 0.5, 1.0, Inf (µg/m^3 per °C).
# Chosen to concentrate resolution in the 0-0.5 range where ~74% of global
# cells and ~79% of CONUS cells fall. Global max = 6.7 but CONUS max = 1.55,
# so upper breaks above 1.0 would be empty for CONUS and waste color bands.
# Boundaries are hardcoded inside cut_alpha() below and in band_labels/band_colors.
#   < 0          green      warming reduces fire PM2.5
#   0 - 0.10     blue       low positive response
#   0.10 - 0.25  dk.yellow  moderate-low
#   0.25 - 0.50  dk.orange  moderate-high
#   0.50 - 1.00  red        high
#   > 1.00       dk.red     extreme (global hotspots: Siberia, Congo, Amazon)

# Bins a continuous alpha vector into the 6 labelled factor bands.
# right = FALSE --> intervals are [lo, hi) so boundary values (0, 0.1, etc.)
# fall in the upper band. include.lowest = TRUE closes the final [1.0, Inf] bin.
cut_alpha <- function(x) {
  cut(x,
      breaks         = c(-Inf, 0, 0.1, 0.25, 0.5, 1.0, Inf),
      labels         = band_labels,
      right          = FALSE,   # [lo, hi) intervals
      include.lowest = TRUE)    # include Inf in last bin
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Map: alpha_i1 ####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# alpha_i1 is the OLS intercept from the grid-cell FE regression:
#   fPM_ips = alpha_i1 + alpha_i2 * T_ps + sum(delta_m * I[model=m])
# It is the CLASSIC-model fitted value at T_ps = 0 — outside the observed range,
# so it has no direct physical interpretation; plotted here for completeness.

map_alpha_i1 <- reg_coefs %>%
  mutate(fill_band = cut_alpha(estimate_alpha_i1)) %>%  # bin intercept into 6 factor bands
  ggplot(aes(x = lon, y = lat, fill = fill_band)) +
  geom_tile(width = 0.5, height = 0.5) +  # explicit 0.5deg tile size; avoids geom_raster() warning from ocean/unmapped gaps
  scale_fill_manual(
    values   = band_colors,            # exact color per band -- no gradient sampling
    name     = expression(alpha[i1]),  # legend title
    drop     = FALSE,                  # show all 6 bands in legend even if a band is empty
    na.value = "lightgray"             # cells with no regression data
  ) +
  coord_fixed(1.3) +
  theme_minimal() +
  labs(
    title    = "Distribution of alpha_i1 across grid cells",
    subtitle = "Fitted fire PM2.5 concentration (µg/m³) at GMT = 0 (CLASSIC baseline; extrapolated)"
  ) +
  theme(panel.grid  = element_blank(),
        axis.title  = element_blank(),
        axis.text   = element_blank(),
        axis.ticks  = element_blank())

map_alpha_i1
ggsave(here("images/regression_alpha/alpha_FE_grid_all", "map_alpha_i1_fe_grid.png"),
       map_alpha_i1, width = 12, height = 6)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Map: alpha_i2 ####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# alpha_i2 is the OLS slope: change in fire PM2.5 concentration (µg/m^3)
# per 1°C GMT increase at this grid cell.
#   alpha_i2 > 0 --> warming increases local fire PM2.5
#   alpha_i2 < 0 --> warming decreases local fire PM2.5

map_alpha_i2 <- reg_coefs %>%
  mutate(fill_band = cut_alpha(estimate_alpha_i2)) %>%  # bin slope into 6 factor bands
  ggplot(aes(x = lon, y = lat, fill = fill_band)) +
  geom_tile(width = 0.5, height = 0.5) +  # explicit 0.5deg tile size; avoids geom_raster() warning from ocean/unmapped gaps
  scale_fill_manual(
    values   = band_colors,            # exact color per band -- no gradient sampling
    name     = expression(alpha[i2]),  # legend title
    drop     = FALSE,                  # show all 6 bands in legend even if a band is empty
    na.value = "lightgray"             # cells with no regression data
  ) +
  coord_fixed(1.3) +
  theme_minimal() +
  labs(
    title    = "Distribution of alpha_i2 across grid cells",
    subtitle = "Change in fire PM2.5 concentration (µg/m³) per 1°C GMT increase"
  ) +
  theme(panel.grid  = element_blank(),
        axis.title  = element_blank(),
        axis.text   = element_blank(),
        axis.ticks  = element_blank())

map_alpha_i2
ggsave(here("images/regression_alpha/alpha_FE_grid_all", "map_alpha_i2_fe_grid.png"),
       map_alpha_i2, width = 12, height = 6)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Map: alpha_i2 — USA only ################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Same scale and color scheme as the global map; filtered to contiguous USA cells.
# reg_coefs_usa defined above in the descriptive stats section.

map_alpha_i2_usa <- reg_coefs_usa %>%
  mutate(fill_band = cut_alpha(estimate_alpha_i2)) %>%  # same cut_alpha() as global map
  ggplot(aes(x = lon, y = lat, fill = fill_band)) +
  geom_tile(width = 0.5, height = 0.5) +  # same explicit tile size as global map
  scale_fill_manual(
    values   = band_colors,            # same exact colors as global map -- directly comparable
    name     = expression(alpha[i2]),  # legend title
    drop     = FALSE,                  # show all 6 bands even if CONUS has none in ">1.00"
    na.value = "lightgray"             # cells with no regression data
  ) +
  coord_fixed(1.3) +
  theme_minimal() +
  labs(
    title    = "Distribution of alpha_i2 across grid cells (CONUS)",
    subtitle = "Change in fire PM2.5 concentration (µg/m³) per 1°C GMT increase"
  ) +
  theme(panel.grid  = element_blank(),
        axis.title  = element_blank(),
        axis.text   = element_blank(),
        axis.ticks  = element_blank())

map_alpha_i2_usa
ggsave(here("images/regression_alpha/alpha_FE_grid_all", "map_alpha_i2_fe_grid_usa.png"),
       map_alpha_i2_usa, width = 10, height = 6)

print("Alpha grid maps saved to images/regression_alpha/alpha_FE_grid_all/")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Line of Best Fit: selected grid cells ####################################
##
## For each target location, find the nearest 0.5° grid cell, rebuild its 25 raw
## observations (18 Park + 5 Pierce + 2 Zhao) from the grid input data and GMT
## scalars, then plot with the fitted line from reg_coefs.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ---- GMT scalars (mirrors fpm_gmt_regress_FE_grid_all.R) ----

gmt_chg  <- read_csv(here("output", "gmt_periods_pi.csv"))
gmt_park <- read_csv(here("output", "gmt_park_decades.csv"))
gmt_zhao <- read_csv(here("output", "gmt_zhao_mimi.csv"))

gmt_baseline  <- mean(c(gmt_chg$mean_gmt_45[gmt_chg$period == "2006-2010"],
                        gmt_chg$mean_gmt_85[gmt_chg$period == "2006-2010"]))
gmt_2040s_45  <- gmt_chg$mean_gmt_45[gmt_chg$period == "2041-2050"]
gmt_2040s_85  <- gmt_chg$mean_gmt_85[gmt_chg$period == "2041-2050"]
gmt_2090s_45  <- gmt_chg$mean_gmt_45[gmt_chg$period == "2091-2100"]
gmt_2090s_85  <- gmt_chg$mean_gmt_85[gmt_chg$period == "2091-2100"]
gmt_2090s_245 <- gmt_zhao$mean_gmt_pi[gmt_zhao$scenario == "SSP245"]
gmt_2090s_585 <- gmt_zhao$mean_gmt_pi[gmt_zhao$scenario == "SSP585"]
gmt_1960s     <- gmt_park$mean_gmt_pi[gmt_park$park_year == 1965]
gmt_1970s     <- gmt_park$mean_gmt_pi[gmt_park$park_year == 1975]
gmt_1980s     <- gmt_park$mean_gmt_pi[gmt_park$park_year == 1985]
gmt_1990s     <- gmt_park$mean_gmt_pi[gmt_park$park_year == 1995]
gmt_2000s     <- gmt_park$mean_gmt_pi[gmt_park$park_year == 2005]
gmt_2010s     <- gmt_park$mean_gmt_pi[gmt_park$park_year == 2015]

# ---- Target locations ----

locations <- tibble(
  city       = c("Denver, CO",  "Los Angeles, CA"),
  target_lon = c(-104.99,       -118.24),
  target_lat = c(  39.74,         34.05)
)

# Find the nearest grid cell to each target using Euclidean distance on lon/lat
locations <- locations %>%
  rowwise() %>%
  mutate(
    nearest     = list(reg_coefs %>%
                         mutate(dist = sqrt((lon - target_lon)^2 + (lat - target_lat)^2)) %>%
                         slice_min(dist, n = 1)),
    cell_row_id = nearest$row_id[1],
    cell_lon    = nearest$lon[1],
    cell_lat    = nearest$lat[1]
  ) %>%
  ungroup()

cat("\nTarget grid cells:\n")
print(locations %>% select(city, target_lon, target_lat, cell_lon, cell_lat, cell_row_id))

# ---- Read grid input data (for raw fPM observations) ----

grid <- read_csv(here("output", "pop_pm_country_death_final.csv"))

# ---- Scenario / model / trajectory lookup ----
# Maps each period_scenario label (from pivot_longer) to its fire model and
# warming trajectory — used to set color and shape aesthetics in the scatter plots.

scenario_lookup <- tibble(
  period_scenario = c(
    "park_classic_1960s_fpm", "park_classic_1970s_fpm",
    "park_classic_1980s_fpm", "park_classic_1990s_fpm",
    "park_classic_2000s_fpm", "park_classic_2010s_fpm",
    "park_jules_1960s_fpm",   "park_jules_1970s_fpm",
    "park_jules_1980s_fpm",   "park_jules_1990s_fpm",
    "park_jules_2000s_fpm",   "park_jules_2010s_fpm",
    "park_ssib4_1960s_fpm",   "park_ssib4_1970s_fpm",
    "park_ssib4_1980s_fpm",   "park_ssib4_1990s_fpm",
    "park_ssib4_2000s_fpm",   "park_ssib4_2010s_fpm",
    "fpm_2000",
    "fpm_2050_45",            "fpm_2050_85",
    "fpm_2100_45",            "fpm_2100_85",
    "fpm_2095_SSP245_Zhao",   "fpm_2095_SSP585_Zhao"
  ),
  model = c(
    rep("CLASSIC", 6), rep("JULES", 6), rep("SSiB4", 6),
    "CESM", "CESM", "CESM", "CESM", "CESM",
    "Multi-model", "Multi-model"
  ),
  trajectory = c(
    rep("Historical", 18),
    "Historical", "RCP4.5", "RCP8.5", "RCP4.5", "RCP8.5",
    "SSP2-4.5", "SSP5-8.5"
  )
)

model_colors <- c(
  "CESM"        = "#E69F00",
  "JULES"       = "#56B4E9",
  "SSiB4"       = "#009E73",
  "CLASSIC"     = "#CC79A7",
  "Multi-model" = "#000000"
)

trajectory_shapes <- c(
  "Historical" = 16,
  "RCP4.5"     = 17,
  "RCP8.5"     = 15,
  "SSP2-4.5"   = 18,
  "SSP5-8.5"   = 8
)

# ---- Build long data and plot for each target cell ----

for (i in seq_len(nrow(locations))) {

  city_name   <- locations$city[i]
  cell_rid    <- locations$cell_row_id[i]
  cell_lon_i  <- locations$cell_lon[i]
  cell_lat_i  <- locations$cell_lat[i]

  cell_grid <- grid %>% filter(row_id == cell_rid)

  # Pierce + Zhao: 7 obs
  pierce_zhao <- cell_grid %>%
    select(row_id, lon, lat,
           fpm_2000, fpm_2050_45, fpm_2050_85,
           fpm_2100_45, fpm_2100_85,
           fpm_2095_SSP245_Zhao, fpm_2095_SSP585_Zhao) %>%
    pivot_longer(cols = starts_with("fpm_"),
                 names_to = "period_scenario", values_to = "fpm") %>%
    mutate(
      T_ps = case_when(
        period_scenario == "fpm_2000"             ~ gmt_baseline,
        period_scenario == "fpm_2050_45"          ~ gmt_2040s_45,
        period_scenario == "fpm_2050_85"          ~ gmt_2040s_85,
        period_scenario == "fpm_2100_45"          ~ gmt_2090s_45,
        period_scenario == "fpm_2100_85"          ~ gmt_2090s_85,
        period_scenario == "fpm_2095_SSP245_Zhao" ~ gmt_2090s_245,
        period_scenario == "fpm_2095_SSP585_Zhao" ~ gmt_2090s_585
      )
    )

  # Park: 18 obs
  park <- cell_grid %>%
    select(row_id, lon, lat,
           starts_with("park_") & ends_with("_fpm")) %>%
    pivot_longer(cols = starts_with("park_"),
                 names_to = "period_scenario", values_to = "fpm") %>%
    mutate(
      T_ps = case_when(
        grepl("1960s", period_scenario) ~ gmt_1960s,
        grepl("1970s", period_scenario) ~ gmt_1970s,
        grepl("1980s", period_scenario) ~ gmt_1980s,
        grepl("1990s", period_scenario) ~ gmt_1990s,
        grepl("2000s", period_scenario) ~ gmt_2000s,
        grepl("2010s", period_scenario) ~ gmt_2010s
      )
    )

  df <- bind_rows(pierce_zhao, park) %>%
    filter(!is.na(fpm)) %>%
    left_join(scenario_lookup, by = "period_scenario") %>%
    mutate(
      model      = factor(model,      levels = names(model_colors)),
      trajectory = factor(trajectory, levels = names(trajectory_shapes)),
      label = case_when(
        period_scenario %in% c(
          "park_classic_2010s_fpm",
          "park_jules_2010s_fpm",
          "park_ssib4_2010s_fpm"
        ) ~ "Park et al.",
        period_scenario %in% c(
          "fpm_2000", "fpm_2050_45", "fpm_2050_85",
          "fpm_2100_45", "fpm_2100_85"
        ) ~ "Pierce et al.",
        period_scenario %in% c(
          "fpm_2095_SSP245_Zhao", "fpm_2095_SSP585_Zhao"
        ) ~ "Zhao et al.",
        TRUE ~ NA_character_
      )
    )

  # Pull coefficients for this cell from reg_coefs
  coef_row    <- reg_coefs %>% filter(row_id == cell_rid)
  alpha_i1    <- coef_row$estimate_alpha_i1
  alpha_i2    <- coef_row$estimate_alpha_i2
  se_alpha_i2 <- coef_row$std.error_alpha_i2
  n_obs       <- nrow(df)

  p <- ggplot(df, aes(x = T_ps, y = fpm, color = model, shape = trajectory)) +
    geom_point(size = 3) +
    geom_text_repel(data = ~ filter(.x, !is.na(label)),
                    aes(label = label),
                    size = 3, box.padding = 0.4, max.overlaps = Inf,
                    show.legend = FALSE) +
    # Fitted line using CLASSIC intercept (alpha_i1); slope alpha_i2 is shared across models
    geom_abline(intercept = alpha_i1, slope = alpha_i2,
                color = "darkred", linewidth = 0.8) +
    annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.5,
             label = "fitted line: CLASSIC intercept",
             color = "darkred", size = 2.5, fontface = "italic") +
    scale_color_manual(values = model_colors,      name = "Model") +
    scale_shape_manual(values = trajectory_shapes, name = "Trajectory") +
    labs(
      title    = paste0(city_name,
                        " (", cell_lon_i, "°, ", cell_lat_i, "°): Fire PM2.5 vs. GMT"),
      subtitle = bquote(alpha[i2] == .(round(alpha_i2, 3)) %+-% .(round(se_alpha_i2, 3)) ~
                          "µg/m^3 per °C  |  n =" ~ .(n_obs)),
      x        = "GMT Anomaly relative to 1850-1900 (°C)",
      y        = "Fire PM2.5 Concentration (µg/m^3)"
    ) +
    theme_minimal() +
    theme(panel.grid.minor = element_blank())

  print(p)

  fname <- paste0("lof_", gsub("[^a-z]", "", tolower(city_name)), "_fe_grid.png")
  ggsave(here("images/regression_alpha/alpha_FE_grid_all", fname),
         p, width = 7, height = 5, dpi = 300)
}

print("Line of best fit plots saved to images/regression_alpha/alpha_FE_grid_all/")

### THE END
