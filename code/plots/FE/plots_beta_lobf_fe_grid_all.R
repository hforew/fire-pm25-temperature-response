# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## PLOTS OF GRID-CELL BETA COEFFICIENTS - FE regression (grid-level)
##
## Goal: Visualise grid-cell beta_i estimates from fpm_gmt_regression_FE_grid_all.R.
##   Produces: (1) histogram of beta_i across all ~65k grid cells;
##   (2) world raster map of beta_i; (3) CONUS raster map of beta_i.
##
## Inputs:
##   output/fpm_gmt_regression_coefs_FE_grid.csv   (beta_i per grid cell with SE,
##                                                   p-value, lon, lat)
##
## Outputs:
##   images/regression_alpha/alpha_FE_grid_all/hist_beta_i_fe_grid.png
##   images/regression_alpha/alpha_FE_grid_all/map_beta_i_fe_grid.png
##   images/regression_alpha/alpha_FE_grid_all/map_beta_i_fe_grid_usa.png
##
## Execution order:
##   files run before: fpm_gmt_regression_FE_grid_all.R   --> writes fpm_gmt_regression_coefs_FE_grid.csv
##   files run after: NA
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Packages #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(here)
library(tidyverse)
library(maps)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Import: coefficients #############################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Grid-cell FE regression coefficients produced by fpm_gmt_regress_FE_grid_all.R.
# One row per 0.5° grid cell (~65k cells). Key columns:
#   row_id, lon, lat      -- cell identifier and centre coordinates
#   estimate_beta_i       -- OLS slope: change in fire PM2.5 (µg/m^3) per 1°C GMT

reg_coefs <- read_csv(here("output", "fpm_gmt_regression_coefs_FE_grid.csv"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Histogram of beta_i ####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

hist_beta_i <- ggplot(reg_coefs, aes(x = estimate_beta_i)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(
    title    = expression("Distribution of " ~ beta[i] ~ "across grid cells"),
    subtitle = "Change in fire PM2.5 concentration (µg/m³) per 1°C GMT increase",
    x        = expression(beta[i] ~ "(µg/m³ per °C)"),
    y        = "Number of grid cells"
  ) +
  theme_minimal()

hist_beta_i
ggsave(here("images/regression_alpha/alpha_FE_grid_all", "hist_beta_i_fe_grid.png"),
       plot = hist_beta_i, width = 8, height = 5, dpi = 300)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Descriptive stats: beta_i distribution ################################
##
## Run before defining cut_beta() breaks to see where the data actually falls.
## Key outputs: percentiles, share negative, share above candidate break thresholds.
## Use these to judge whether break values capture meaningful variation.
##
## Stats computed for: (1) all global land cells, (2) CONUS cells only.
## cut_beta() applies to both maps, so check both distributions before setting breaks.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

thresholds <- c(0, 0.5, 1, 1.5, 2)

print_beta_stats <- function(x, label) {
  cat("\n--- beta_i descriptive stats:", label, "(n =", sum(!is.na(x)), "cells) ---\n")
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

print_beta_stats(reg_coefs$estimate_beta_i, "GLOBAL")

# CONUS: filtered here for stats and reused for the CONUS map below
reg_coefs_usa <- reg_coefs %>%
  filter(country_code_iso3 == "USA",  # keep only cells mapped to USA
         lon >= -125, lon <= -66,      # CONUS longitude bounds (excludes Alaska)
         lat >=   24, lat <=  50)      # CONUS latitude bounds (excludes Hawaii)

print_beta_stats(reg_coefs_usa$estimate_beta_i, "CONUS")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Map color palette  #################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Exact color assigned to each band -- no gradient interpolation.
# cut() bins the continuous beta value into a factor; scale_fill_manual()
# maps each factor level to its exact hex color. Both global and CONUS maps
# use the same band_labels and band_colors, so identical values always
# produce identical colors regardless of which cells are in the plot.
band_labels <- c("<0", "0-0.10", "0.10-0.25", "0.25-0.50", "0.50-1.00", ">1.00")

band_colors <- c(
  # "<0"        = "#F5F0E8",  # off-white    -- warming reduces fire PM2.5
  "<0"        = "azure2",  # azure2    -- warming reduces fire PM2.5
  "0-0.10"    = "#FFE566",  # light yellow -- low positive response
  "0.10-0.25" = "#C8A000",  # dark yellow  -- moderate-low
  "0.25-0.50" = "#CC5500",  # dark orange  -- moderate-high
  "0.50-1.00" = "#DD1111",  # red          -- high
  ">1.00"     = "#7A0000"   # dark red     -- extreme (Siberia, Congo, Amazon globally)
)

# Bins a continuous beta vector into the 6 labelled factor bands.
# right = FALSE --> boundary values fall in the upper band; include.lowest closes [1.0, Inf].
cut_beta <- function(x) {
  cut(x,
      breaks         = c(-Inf, 0, 0.1, 0.25, 0.5, 1.0, Inf),
      labels         = band_labels,
      right          = FALSE,
      include.lowest = TRUE)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Map: beta_i ####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Country boundaries drawn as a gray polygon overlay on top of the raster tiles.
world_borders <- map_data("world")

# beta_i is the OLS slope: change in fire PM2.5 concentration (µg/m^3)
# per 1°C GMT increase at this grid cell.
#   beta_i > 0 --> warming increases local fire PM2.5
#   beta_i < 0 --> warming decreases local fire PM2.5

map_beta_i <- reg_coefs %>%
  mutate(fill_band = cut_beta(estimate_beta_i)) %>%
  ggplot(aes(x = lon, y = lat, fill = fill_band)) +
  geom_tile(width = 0.5, height = 0.5) +  # avoids geom_raster() warning from ocean/unmapped gaps
  geom_polygon(data = world_borders,
               aes(x = long, y = lat, group = group),
               fill = NA, color = "gray50", linewidth = 0.15) +
  scale_fill_manual(
    values   = band_colors,
    name     = NULL,
    drop     = FALSE,                 # show all 6 bands in legend even if a band is empty
    na.value = "lightgray"            # cells with no regression data
  ) +
  coord_fixed(0.90) +
  theme_minimal() +
  labs(
    title = "Grid cell estimates: Change in fire PM2.5 concentration (µg/m³) per 1°C GMT increase"
      ) +
  guides(fill = guide_legend(nrow = 1)) +
  theme(panel.grid       = element_blank(),
        panel.background = element_blank(),
        axis.title       = element_blank(),
        axis.text        = element_blank(),
        axis.ticks       = element_blank(),
        legend.position  = "bottom",
        plot.title       = element_text(size = 9.5))

map_beta_i
ggsave(here("images/regression_alpha/alpha_FE_grid_all", "map_beta_i_fe_grid.png"),
       map_beta_i, width = 6.5, height = 3.5, dpi = 300)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Map: beta_i — USA only ################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Same scale and color scheme as the global map; filtered to contiguous USA cells.
# reg_coefs_usa defined above in the descriptive stats section.

usa_state_borders <- map_data("state")

map_beta_i_usa <- reg_coefs_usa %>%
  mutate(fill_band = cut_beta(estimate_beta_i)) %>%
  ggplot(aes(x = lon, y = lat, fill = fill_band)) +
  geom_tile(width = 0.5, height = 0.5) +
  geom_polygon(data = usa_state_borders,
               aes(x = long, y = lat, group = group),
               fill = NA, color = "gray50", linewidth = 0.15) +
  scale_fill_manual(
    values   = band_colors,           # same exact colors as global map -- directly comparable
    name     = NULL,
    drop     = FALSE,                 # show all 6 bands even if CONUS has none in ">1.00"
    na.value = "lightgray"            # cells with no regression data
  ) +
  coord_fixed(0.97, ylim = c(24, 49)) +  # clip at 49°N (US-Canada border) to hide tiles overhanging state polygons
  theme_minimal() +
  labs(
    title = "Grid cell estimates: Change in fire PM2.5 concentration (µg/m³) per 1°C GMT increase"
  ) +
  guides(fill = guide_legend(nrow = 1)) +
  theme(panel.grid       = element_blank(),
        panel.background = element_blank(),
        axis.title       = element_blank(),
        axis.text        = element_blank(),
        axis.ticks       = element_blank(),
        legend.position  = "bottom",
        plot.title       = element_text(size = 9.5))

map_beta_i_usa
ggsave(here("images/regression_alpha/alpha_FE_grid_all", "map_beta_i_fe_grid_usa.png"),
       map_beta_i_usa, width = 6.5, height = 3.5, dpi = 300)

print("Beta grid maps saved to images/regression_alpha/alpha_FE_grid_all/")

### THE END ##########
