# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## PLOTS OF COUNTRY ALPHA COEFFICIENTS - ##
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls())

# Visualises country-level FE regression outputs from fpm_gmt_regression_FE_all.R.
# Produces: (1) histogram of alpha_c2 distribution; (2) world choropleth maps of
# alpha_c1 and alpha_c2; (3) individual country LOBF scatter plots; (4) 2x3
# multi-country LOBF grid. All outputs saved to images/regression_alpha/alpha_FE_all/.

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Packages #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(here)
library(tidyverse)
library(ggplot2)
library(maps)
library(countrycode)
library(ggtext)
library(ggrepel)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Import #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Import country-level regression coefficients produced by fpm_gmt_regression_FE_all.R.
# Each row is one country. Key columns:
#   country_code_iso3   -- ISO 3166-1 alpha-3 country code (used for map joining)
#   country_name        -- full country name
#   estimate_alpha_c1   -- OLS intercept: fitted exposure at T_ps = 0 (extrapolated
#                          outside the data range; not physically meaningful)
#   estimate_alpha_c2   -- OLS slope: change in per-capita fire PM2.5 exposure
#                          (µg/m^3/yr) per 1°C increase in GMT. This is the key
#                          damage function parameter used downstream in GIVE.
#                          Positive = more fire PM with warming; negative = less.

source(here("code/regression_alpha", "fpm_gmt_regression_FE_all.R")) # regression file and load all objects

head(reg_coefs)
str(reg_coefs)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Histogram of alpha_c2 ####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

hist_alpha_c2 <- ggplot(reg_coefs, aes(x = estimate_alpha_c2)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(
    title    = "Distribution of alpha_c2 across countries",
    subtitle = "Change in per-capita fire PM2.5 (µg/m³/yr) per 1°C GMT increase",
    x        = "alpha_c2 (µg/m³/yr per °C)",
    y        = "Number of countries"
  ) +
  theme_minimal()

hist_alpha_c2
ggsave(here("images/regression_alpha/alpha_FE_all", "hist_alpha_c2_fe.png"), plot = hist_alpha_c2, width = 8, height = 5, dpi = 300)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Prep map data #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load world polygon data from the maps package.
# Each row is one vertex of a country polygon; group identifies which polygon
# a vertex belongs to (needed by geom_polygon to draw borders correctly).
world_map <- map_data("world")

# Convert map region names (e.g., "United States", "Russia") to ISO3 codes.
# This avoids name-mismatch failures when joining to reg_coefs, which uses
# ISO3 codes (e.g., "USA", "RUS"). The countrycode package handles common
# discrepancies like "United States" vs "United States of America".
# warn = FALSE suppresses messages for territories that have no ISO3 code
# (e.g., Antarctica, disputed regions) — these will join as NA and show gray.
world_map$iso3 <- countrycode(world_map$region,
                               origin = "country.name",
                               destination = "iso3c",
                               warn = FALSE)

# Left join alpha coefficients onto the map polygon vertices by ISO3 code.
# Every polygon vertex row gets the country's alpha_c1 and alpha_c2 values.
# Countries not present in reg_coefs (no regression data) will have NA
# and will be rendered in the na.value color (lightgray) on the map.
world_alpha <- world_map %>%
  left_join(reg_coefs, by = c("iso3" = "country_code_iso3"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Descriptive stats: alpha_c2 distribution ################################
##
## Run before defining color palette breaks to see where the data actually falls.
## Key outputs: percentiles, share negative, share above candidate break thresholds.
## Use these to judge whether break values capture meaningful variation.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

thresholds <- c(0, 0.5, 1, 1.5, 2)

print_alpha_stats <- function(x, label) {
  cat("\n--- alpha_c2 descriptive stats:", label, "(n =", sum(!is.na(x)), "countries) ---\n")
  cat("Min:  ", round(min(x,  na.rm = TRUE), 3), "\n")
  cat("Max:  ", round(max(x,  na.rm = TRUE), 3), "\n")
  cat("Mean: ", round(mean(x, na.rm = TRUE), 3), "\n")
  cat("SD:   ", round(sd(x,   na.rm = TRUE), 3), "\n")
  pctiles <- quantile(x, probs = c(0.02, 0.05, 0.10, 0.25, 0.50,
                                    0.75, 0.90, 0.95, 0.98), na.rm = TRUE)
  cat("\nPercentiles:\n")
  print(round(pctiles, 3))
  cat("\nShare of countries below each candidate break (%):\n")
  for (thr in thresholds) {
    cat("  <", thr, ":", round(100 * mean(x < thr, na.rm = TRUE), 1), "%\n")
  }
  cat("  > 2 :", round(100 * mean(x > 2, na.rm = TRUE), 1), "%\n")
}

print_alpha_stats(reg_coefs$estimate_alpha_c2, "GLOBAL")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Map color palette  #################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

band_labels <- c("<0", "0-0.10", "0.10-0.25", "0.25-0.50", "0.50-1.00", ">1.00")

band_colors <- c(
  "<0"        = "#7FD67F",  # light green  -- warming reduces fire PM2.5
  "0-0.10"    = "#56C8FF",  # sky blue     -- low positive response
  "0.10-0.25" = "#C8A000",  # dark yellow  -- moderate-low
  "0.25-0.50" = "#CC5500",  # dark orange  -- moderate-high
  "0.50-1.00" = "#DD1111",  # red          -- high
  ">1.00"     = "#7A0000"   # dark red     -- extreme
)

# Bins a continuous alpha value into the 6 labelled factor bands.
# right = FALSE --> intervals are [lo, hi) so boundary values fall in the upper band.
# include.lowest = TRUE closes the final [1.0, Inf] bin.
cut_alpha <- function(x) {
  cut(x,
      breaks         = c(-Inf, 0, 0.1, 0.25, 0.5, 1.0, Inf),
      labels         = band_labels,
      right          = FALSE,
      include.lowest = TRUE)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Plot alpha_c1 #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# alpha_c1 is the OLS y-intercept from the country-level regression with fire-model FEs:
#   exposure_percap = alpha_c1 + alpha_c2 * T_ps + sum(delta_m * I[model=m])
# It represents the classic-model baseline intercept (T_ps = 0, model = classic).
# Because all observed T_ps values are ~1°C or above, T = 0 lies outside
# the data range — alpha_c1 is an extrapolation and has no direct physical
# interpretation. It is mapped here for completeness.

map_alpha_c1 <- world_alpha %>%
  mutate(fill_band = cut_alpha(estimate_alpha_c1)) %>%
  ggplot(aes(x = long, y = lat, group = group, fill = fill_band)) +
  geom_polygon(color = "white", linewidth = 0.1) +
  scale_fill_manual(
    values   = band_colors,
    name     = expression(alpha[c1]),
    drop     = FALSE,
    na.value = "lightgray"
  ) +
  # coord_fixed(1.3) preserves geographic aspect ratio (lat/lon are not
  # equal-distance units; 1.3 approximates a Mercator-like stretch).
  coord_fixed(1.3) +
  theme_minimal() +
  labs(title = expression("Country-Level Intercept Coefficient (" * alpha[c1] * "): FPM-GMT Regression"),
       x = NULL, y = NULL) +
  theme(panel.grid = element_blank())


map_alpha_c1
ggsave(here("images/regression_alpha/alpha_FE_all", "map_alpha_c1_fe.png"),
       map_alpha_c1, width = 12, height = 6)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Plot alpha_c2 #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# alpha_c2 is the OLS slope from the country-level regression with fire-model FEs:
#   exposure_percap = alpha_c1 + alpha_c2 * T_ps + sum(delta_m * I[model=m])
#
# It is the key damage function parameter: the change in per-capita fire PM2.5
# exposure (µg/m³/yr) per 1°C increase in global mean temperature (GMT).
#
# Interpretation:
#   alpha_c2 > 0  --> warming increases fire PM2.5 exposure in that country
#   alpha_c2 < 0  --> warming decreases fire PM2.5 exposure in that country
#   alpha_c2 = 0  --> no relationship between GMT and fire PM2.5 in that country
#
# This coefficient is used downstream in GIVE to translate GMT change scenarios
# into country-level fire PM2.5 exposure changes and associated mortality impacts.

map_alpha_c2 <- world_alpha %>%
  mutate(fill_band = cut_alpha(estimate_alpha_c2)) %>%
  ggplot(aes(x = long, y = lat, group = group, fill = fill_band)) +
  geom_polygon(color = "white", linewidth = 0.1) +
  scale_fill_manual(
    values   = band_colors,
    name     = expression(alpha[c2]),
    drop     = FALSE,
    na.value = "lightgray"
  ) +
  # coord_fixed(1.3) preserves geographic aspect ratio.
  coord_fixed(1.3) +
  theme_minimal() +
  labs(title    = "Distribution of alpha_c2 across countries",
       subtitle = "Change in per-capita fire PM2.5 (µg/m³/yr) per 1°C GMT increase",
       x = NULL, y = NULL) +
  theme(panel.grid = element_blank())


map_alpha_c2
ggsave(here("images/regression_alpha/alpha_FE_all", "map_alpha_c2_fe.png"),
       map_alpha_c2, width = 12, height = 6)

print("Alpha maps saved to images/regression_alpha/alpha_FE_all/")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Line of Best Fit #########################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# For each selected country, plot all 25 observations (18 Park + 5 Pierce + 2 Zhao) with the
# OLS fitted line: exposure_percap = alpha_c1 + alpha_c2 * T_ps + sum(delta_m * I[model=m])
#
# Data objects from source():
#   reg_data_combined -- 25 rows per country: T_ps (GMT regressor) and exposure_percap (response)
#   reg_coefs         -- alpha_c1 and alpha_c2 per country (define the fitted line)

# Lookup table: maps each period_scenario to its fire model and warming trajectory.
# Used to set color (model) and shape (trajectory) aesthetics in the scatter plots.
scenario_lookup <- tibble(
  period_scenario = c(
    "exposure_percap_park_classic_1960s_fpm", "exposure_percap_park_classic_1970s_fpm",
    "exposure_percap_park_classic_1980s_fpm", "exposure_percap_park_classic_1990s_fpm",
    "exposure_percap_park_classic_2000s_fpm", "exposure_percap_park_classic_2010s_fpm",
    "exposure_percap_park_jules_1960s_fpm",   "exposure_percap_park_jules_1970s_fpm",
    "exposure_percap_park_jules_1980s_fpm",   "exposure_percap_park_jules_1990s_fpm",
    "exposure_percap_park_jules_2000s_fpm",   "exposure_percap_park_jules_2010s_fpm",
    "exposure_percap_park_ssib4_1960s_fpm",   "exposure_percap_park_ssib4_1970s_fpm",
    "exposure_percap_park_ssib4_1980s_fpm",   "exposure_percap_park_ssib4_1990s_fpm",
    "exposure_percap_park_ssib4_2000s_fpm",   "exposure_percap_park_ssib4_2010s_fpm",
    "exposure_percap_fpm_2000",
    "exposure_percap_fpm_2050_45",            "exposure_percap_fpm_2050_85",
    "exposure_percap_fpm_2100_45",            "exposure_percap_fpm_2100_85",
    "exposure_percap_fpm_2095_SSP245_Zhao",   "exposure_percap_fpm_2095_SSP585_Zhao"
  ),
  model = c(
    rep("CLASSIC", 6), rep("JULES", 6), rep("SSiB4", 6),
    "CESM", "CESM", "CESM", "CESM", "CESM",
    "Multi", "Multi"
  ),
  trajectory = c(
    rep("Historical", 18),
    "Historical", "SSP1-4.5", "SSP3-8.5", "SSP1-4.5", "SSP3-8.5",
    "SSP2-4.5", "SSP5-8.5"
  )
)

model_colors <- c(
  "CESM"        = "#E69F00",
  "JULES"       = "#56B4E9",
  "SSiB4"       = "#009E73",
  "CLASSIC"     = "#CC79A7",
  "Multi"    = "#000000"
)

trajectory_shapes <- c(
  "Historical" = 16,
  "SSP1-4.5"     = 17,
  "SSP3-8.5"     = 15,
  "SSP2-4.5"   = 18,
  "SSP5-8.5"   = 8
)

countries_to_plot <- c("global", "USA","DEU", "RUS", "AUS", "CHN", "IND", "ARG", "BRA", "ETH", "NGA")

for (iso in countries_to_plot) {

  # Filter to this country's 25 observations, join model/trajectory lookup, and
  # create point labels. Labels shown for all Pierce/Zhao points and Park 2010s only;
  # all other Park decades are unlabelled to avoid overplotting.
  df <- reg_data_combined %>%
    filter(country_code_iso3 == iso) %>%
    left_join(scenario_lookup, by = "period_scenario") %>%
    mutate(
      model      = factor(model,      levels = names(model_colors)),
      trajectory = factor(trajectory, levels = names(trajectory_shapes)),
      label = case_when(
        period_scenario %in% c(
          "exposure_percap_park_classic_2010s_fpm",
          "exposure_percap_park_jules_2010s_fpm",
          "exposure_percap_park_ssib4_2010s_fpm"
        ) ~ "Park et al.",
        period_scenario %in% c(
          "exposure_percap_fpm_2000",
          "exposure_percap_fpm_2050_45",
          "exposure_percap_fpm_2050_85",
          "exposure_percap_fpm_2100_45",
          "exposure_percap_fpm_2100_85"
        ) ~ "Pierce et al.",
        period_scenario %in% c(
          "exposure_percap_fpm_2095_SSP245_Zhao",
          "exposure_percap_fpm_2095_SSP585_Zhao"
        ) ~ "Zhao et al.",
        TRUE ~ NA_character_
      )
    )

  # Pull alpha_c1 (intercept), alpha_c2 (slope), and its standard error from reg_coefs
  coef_row    <- reg_coefs %>% filter(country_code_iso3 == iso)
  alpha_c1    <- coef_row$estimate_alpha_c1
  alpha_c2    <- coef_row$estimate_alpha_c2
  se_alpha_c2 <- coef_row$std.error_alpha_c2

  # Pull R² from the glanced model-fit statistics and observation count from df
  r_squared <- reg_results %>%
    filter(country_code_iso3 == iso) %>%
    unnest(glanced) %>%
    pull(r.squared)
  n_obs <- nrow(df)

  p <- ggplot(df, aes(x = T_ps, y = exposure_percap, color = model, shape = trajectory)) +
    # Each point is one observation: 18 Park (decade × fire model) +
    # 5 Pierce (baseline, 2040s SSP1-4.5/SSP3-8.5, 2090s SSP1-4.5/SSP3-8.5) + 2 Zhao (~2095 SSP2-4.5/SSP5-8.5)
    geom_point(size = 3) +
    # Labels for all labelled points (Park 2010s, Pierce, Zhao).
    # ggrepel automatically nudges overlapping labels apart and draws
    # a line back to the point when a label is pushed away.
    geom_text_repel(data = ~ filter(.x, !is.na(label)),
                    aes(label = label),
                    size = 3, box.padding = 0.4, max.overlaps = Inf,
                    show.legend = FALSE) +
    # Classic-model reference line: intercept = alpha_c1, slope = alpha_c2.
    # With FEs, each model has its own intercept (alpha_c1 + delta_m); this line
    # shows the classic baseline. alpha_c2 (slope) is shared across all models.
    geom_abline(intercept = alpha_c1, slope = alpha_c2,  # fitted line w/ CLASSIC intercept
                color = "darkred", linewidth = 0.8) +
    annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.5,
             label = "fitted line: CLASSIC intercept",
             color = "darkred", size = 2.5, fontface = "italic") +
    scale_color_manual(values = model_colors,      name = "Model") +
    scale_shape_manual(values = trajectory_shapes, name = "Trajectory") +
    labs(
      title    = paste0(iso, ": Per-Capita Fire PM2.5 vs. GMT"),
      subtitle = bquote(alpha[c2] == .(round(alpha_c2, 2)) %+-% .(round(se_alpha_c2, 2)) ~
                          "µg/m³/yr per °C  |  n =" ~ .(n_obs) ~
                          "|  R² =" ~ .(round(r_squared, 2))),
      x        = "GMT Anomaly relative to 1850-1900 (°C)",
      y        = "Per-Capita Fire PM2.5 Exposure (µg/m³/yr)"
    ) +
    theme_minimal() +
    theme(panel.grid.minor  = element_blank(),
          plot.subtitle     = element_text(size = 12))

  print(p)
  ggsave(here("images/regression_alpha/alpha_FE_all", paste0("lof_", tolower(iso), "_fe.png")),
         p, width = 7, height = 5, dpi = 300)
}


print("Line of best fit plots saved to images/regression_alpha/alpha_FE_all/")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Line of Best Fit: multi-country grid (2 x 3) ############################
##
## Six-panel faceted plot: 2 columns x 3 rows, free y-axis scales, single legend.
## Strip label: country code + alpha_c2 + n.
## Fitted line uses CLASSIC model intercept (alpha_c1); slope alpha_c2 shared.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

countries_to_multiplot <- c("USA", "RUS", "CHN", "IND", "BRA", "NGA")

# Build combined long data for all 6 countries with strip_label column for faceting.
# map_dfr iterates over countries_to_multiplot; each element is passed as iso
# (the ISO3 code) to the anonymous function, which builds that country's
# 25-row data frame. Results are row-bound into one combined data frame.
df_multi <- purrr::map_dfr(countries_to_multiplot, function(iso) {
  coef_row    <- reg_coefs %>% filter(country_code_iso3 == iso)
  alpha_c2    <- coef_row$estimate_alpha_c2    # slope for strip label
  se_alpha_c2 <- coef_row$std.error_alpha_c2   # SE for strip label

  reg_data_combined %>%
    filter(country_code_iso3 == iso) %>%
    left_join(scenario_lookup, by = "period_scenario") %>%
    mutate(
      model       = factor(model,      levels = names(model_colors)),
      trajectory  = factor(trajectory, levels = names(trajectory_shapes)),
      # HTML subscript rendered by element_markdown(); Unicode ± for plus-minus
      strip_label = paste0(iso, "  |  α<sub>c2</sub> = ", round(alpha_c2, 2),
                           " ± ", round(se_alpha_c2, 2),
                           " µg/m³/yr per °C  |  n = ", n())
    )
})

# facet_wrap orders panels alphabetically by default; set factor levels to
# preserve the countries_to_multiplot order (USA, RUS, CHN, IND, BRA, NGA).
strip_levels <- df_multi %>%
  distinct(country_code_iso3, strip_label) %>%
  arrange(match(country_code_iso3, countries_to_multiplot)) %>%
  pull(strip_label)

df_multi <- df_multi %>%
  mutate(strip_label = factor(strip_label, levels = strip_levels))

# lines_multi must contain strip_label (the faceting variable) so geom_abline
# subsets it per panel and draws the correct country-specific fitted line.
lines_multi <- reg_coefs %>%
  filter(country_code_iso3 %in% countries_to_multiplot) %>%
  left_join(df_multi %>% distinct(country_code_iso3, strip_label),
            by = "country_code_iso3")

p_multi <- ggplot(df_multi,
                  aes(x = T_ps, y = exposure_percap, color = model, shape = trajectory)) +
  geom_point(size = 2) +
  # geom_abline with data = lines_multi draws a per-panel line using each
  # country's alpha_c1 (CLASSIC intercept) and shared slope alpha_c2.
  geom_abline(data  = lines_multi,
              aes(intercept = estimate_alpha_c1, slope = estimate_alpha_c2),
              color = "darkred", linewidth = 0.7) +  # fitted line w/ CLASSIC intercept
  # annotate() draws in every panel; Inf/-Inf anchors text to bottom-right corner.
  annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.5,
           label = "fitted line: CLASSIC intercept",
           color = "darkred", size = 2, fontface = "italic") +
  scale_color_manual(values = model_colors,      name = "Model") +
  scale_shape_manual(values = trajectory_shapes, name = "Trajectory") +
  guides(shape = guide_legend(order = 1),   # Trajectory row first
         color = guide_legend(order = 2)) + # Model row second
  # ncol = 2 --> 2x3 grid; scales = "free_y" --> each panel scaled to its own y range.
  facet_wrap(~ strip_label, ncol = 2, scales = "free_y") +
  labs(
    title = "Per-Capita Fire PM2.5 vs. GMT",
    x     = "GMT Anomaly relative to 1850-1900 (°C)",
    y     = "Per-Capita Fire PM2.5 Exposure (µg/m³/yr)"
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    strip.text       = element_markdown(size = 9),  # element_markdown renders HTML subscripts (ggtext)
    legend.position  = "bottom",                    # single shared legend below the grid
    legend.box       = "vertical",                  # stack Trajectory and Model rows vertically
    legend.spacing.y = unit(0.1, "cm"),             # tighten gap between the two legend rows
    plot.margin      = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm")
  )

print(p_multi)
ggsave(here("images/regression_alpha/alpha_FE_all", "lof_multiplot_fe.png"),
       p_multi, width = 8.5, height = 11, dpi = 300)

print("Multi-country grid plot saved to images/regression_alpha/alpha_FE_all/lof_multiplot_fe.png")

### THE END
