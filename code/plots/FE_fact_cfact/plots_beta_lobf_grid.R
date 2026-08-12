# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## PLOTS OF GRID-CELL BETA COEFFICIENTS - FE leave-one-model-out sensitivity (grid-level)
##
## Goal: Visualise grid-cell beta_i estimates from each exclusion group produced by
##   FE_sensitivity/fpm_gmt_regress_grid.R. For each of the 6 groups (full, excl_classic,
##   excl_jules, excl_ssib4, excl_CESM, excl_Zhao), produces:
##   (1) histogram of beta_i across all grid cells;
##   (2) world raster map of beta_i; (3) CONUS raster map of beta_i.
##   Same fixed color bands and breaks across all groups, so maps are directly comparable
##   to see how sensitive the spatial pattern is to dropping any one fire-model source.
##
## Inputs:
##   output/betas_fe_sensitivity/fpm_gmt_betas_grid_<group>.csv   (beta_i per grid cell with SE,
##                                                   p-value, lon, lat, per exclusion group)
##
## Outputs (one triplet per exclusion group):
##   images/regression_beta/beta_grid/<group>/hist_beta_i_grid.png
##   images/regression_beta/beta_grid/<group>/fig1_map_beta_i_grid.png
##   images/regression_beta/beta_grid/<group>/fig6_map_beta_i_grid_usa.png
##
## Execution order:
##   files run before: FE_sensitivity/fpm_gmt_regress_grid.R   --> writes fpm_gmt_betas_grid_<group>.csv
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
############ Exclusion groups #############################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Same 6 groups as FE_sensitivity/fpm_gmt_regress_grid.R -- one input file per group.
group_names <- c("full", "excl_classic", "excl_jules", "excl_ssib4", "excl_CESM", "excl_Zhao")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Map color palette (fixed across all groups) ###################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Exact color assigned to each band -- no gradient interpolation.
# cut() bins the continuous beta value into a factor; scale_fill_manual()
# maps each factor level to its exact hex color. Breaks are fixed (not recomputed
# per group) so identical beta values always produce identical colors across
# groups -- maps stay directly comparable.
band_labels <- c("<0", "0-0.10", "0.10-0.25", "0.25-0.50", "0.50-1.00", ">1.00")

band_colors <- c(
  "<0"        = "azure2",   # warming reduces fire PM2.5
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
############ Descriptive stats helper #######################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Key outputs: percentiles, share negative, share above candidate break thresholds.
## Run for GLOBAL and CONUS within every exclusion group, to see how sensitive the
## beta_i distribution is to dropping each fire-model source.
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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Borders (loaded once, reused across all groups) ################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

world_borders     <- map_data("world")
usa_state_borders <- map_data("state")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Plot loop: one histogram + 2 maps per exclusion group ##########
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for (group_name in group_names) {

  # group_name <- "full"  # uncomment + set to run a single group manually, outside the loop

  cat("\n============================================================\n")
  cat("[", group_name, "] Grid-cell beta_i plots\n", sep = "")
  cat("============================================================\n")

  # ---- Import: coefficients for this group ----
  # One row per grid cell. Key columns:
  #   row_id, lon, lat      -- cell identifier and centre coordinates
  #   estimate_beta_i       -- OLS slope: change in fire PM2.5 (µg/m^3) per 1°C GMT
  reg_coefs <- read_csv(here("output", "betas_fe_sensitivity",
                              paste0("fpm_gmt_betas_grid_", group_name, ".csv")))

  out_dir <- here("images", "regression_beta", "beta_grid", group_name)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  # ---- Histogram of beta_i ----
  hist_beta_i <- ggplot(reg_coefs, aes(x = estimate_beta_i)) +
    geom_histogram(bins = 50, fill = "steelblue", color = "white") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    labs(
      title    = expression("Distribution of " ~ beta[i] ~ "across grid cells"),
      subtitle = paste0("Change in fire PM2.5 concentration (µg/m³) per 1°C GMT increase  [",
                         group_name, "]"),
      x        = expression(beta[i] ~ "(µg/m³ per °C)"),
      y        = "Number of grid cells"
    ) +
    theme_minimal()

  ggsave(file.path(out_dir, "hist_beta_i_grid.png"),
         plot = hist_beta_i, width = 8, height = 5, dpi = 300)

  # ---- Descriptive stats: GLOBAL and CONUS ----
  print_beta_stats(reg_coefs$estimate_beta_i, paste0(group_name, " - GLOBAL"))

  reg_coefs_usa <- reg_coefs %>%
    filter(country_code_iso3 == "USA",  # keep only cells mapped to USA
           lon >= -125, lon <= -66,      # CONUS longitude bounds (excludes Alaska)
           lat >=   24, lat <=  50)      # CONUS latitude bounds (excludes Hawaii)

  print_beta_stats(reg_coefs_usa$estimate_beta_i, paste0(group_name, " - CONUS"))

  # ---- Map: beta_i, global ----
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
      title = paste0("Grid cell estimates: Change in fire PM2.5 concentration (µg/m³) per 1°C GMT increase  [",
                      group_name, "]")
    ) +
    guides(fill = guide_legend(nrow = 1)) +
    theme(panel.grid       = element_blank(),
          panel.background = element_blank(),
          axis.title       = element_blank(),
          axis.text        = element_blank(),
          axis.ticks       = element_blank(),
          legend.position  = "bottom",
          plot.title       = element_text(size = 9.5))

  ggsave(file.path(out_dir, "fig1_map_beta_i_grid.png"),
         plot = map_beta_i, width = 6.5, height = 3.5, dpi = 300)

  # ---- Map: beta_i, USA only ----
  # Same fixed scale and color scheme as the global map; filtered to contiguous USA cells.
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
      title = paste0("Grid cell estimates: Change in fire PM2.5 concentration (µg/m³) per 1°C GMT increase  [",
                      group_name, "]")
    ) +
    guides(fill = guide_legend(nrow = 1)) +
    theme(panel.grid       = element_blank(),
          panel.background = element_blank(),
          axis.title       = element_blank(),
          axis.text        = element_blank(),
          axis.ticks       = element_blank(),
          legend.position  = "bottom",
          plot.title       = element_text(size = 9.5))

  ggsave(file.path(out_dir, "fig6_map_beta_i_grid_usa.png"),
         plot = map_beta_i_usa, width = 6.5, height = 3.5, dpi = 300)

  cat("[", group_name, "] Plots saved to images/regression_beta/beta_grid/",
      group_name, "/\n", sep = "")
}

print("All grid sensitivity beta plots saved to images/regression_beta/beta_grid/")

### THE END ##########
