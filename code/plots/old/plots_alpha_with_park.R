# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## PLOTS OF COUNTRY ALPHA COEFFICIENTS - ##
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Packages #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(here)
library(tidyverse)
library(ggplot2)
library(maps)
library(countrycode)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Import #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Import country-level regression coefficients produced by fpm_gmt_relation_park.R.
# Each row is one country. Key columns:
#   country_code_iso3   -- ISO 3166-1 alpha-3 country code (used for map joining)
#   country_name        -- full country name
#   estimate_alpha_c1   -- OLS intercept: fitted exposure at T_ps = 0 (extrapolated
#                          outside the data range; not physically meaningful)
#   estimate_alpha_c2   -- OLS slope: change in per-capita fire PM2.5 exposure
#                          (µg/m³/yr) per 1°C increase in GMT. This is the key
#                          damage function parameter used downstream in GIVE.
#                          Positive = more fire PM with warming; negative = less.

source(here("code/regression_alpha", "fpm_gmt_relation_park.R")) # regression file and load all objects

# alpha_coefs <- read_csv(here("output", "fpm_gmt_regression_coefs.csv"))

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
ggsave(here("images/regression_alpha", "hist_alpha_c2_park.png"), plot = hist_alpha_c2, width = 8, height = 5, dpi = 300)


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
############ Plot alpha_c1 #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# alpha_c1 is the OLS y-intercept from the country-level regression:
#   exposure_percap = alpha_c1 + alpha_c2 * T_ps
# It represents the fitted line's value at T_ps = 0 (pre-industrial GMT).
# Because all observed T_ps values are ~1°C or above, T = 0 lies outside
# the data range — alpha_c1 is an extrapolation and has no direct physical
# interpretation. It is mapped here for completeness.

map_alpha_c1 <- ggplot() +
  # Draw country polygons filled by alpha_c1.
  # group = group ensures each country polygon is drawn as a closed shape
  # (without this, ggplot connects vertices across countries incorrectly).
  # color = "white" draws thin white country borders; linewidth = 0.1 keeps
  # borders subtle so the fill color is the visual focus.
  geom_polygon(data = world_alpha,
               aes(x = long, y = lat, group = group, fill = estimate_alpha_c1),
               color = "white", linewidth = 0.1) +
  # Diverging color scale centered at 0.
  # Blue = negative alpha_c1 (intercept below zero),
  # red  = positive alpha_c1 (intercept above zero),
  # gray = countries with no regression data (NA).
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "darkred",
    midpoint = 0,
    name = "Alpha C1",
    na.value = "lightgray"
  ) +
  # coord_fixed(1.3) preserves geographic aspect ratio (lat/lon are not
  # equal-distance units; 1.3 approximates a Mercator-like stretch).
  coord_fixed(1.3) +
  theme_minimal() +
  labs(title = "Country-Level Alpha Coefficients (C1): FPM-GMT Regression",
       x = "Longitude", y = "Latitude") +
  theme(panel.grid = element_blank())


map_alpha_c1
ggsave(here("images/regression_alpha", "map_alpha_c1_park.png"),
       map_alpha_c1, width = 12, height = 6)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Plot alpha_c2 #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# alpha_c2 is the OLS slope from the country-level regression:
#   exposure_percap = alpha_c1 + alpha_c2 * T_ps
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

map_alpha_c2 <- ggplot() +
  # Draw country polygons filled by alpha_c2.
  # group = group ensures each country polygon is drawn as a closed shape.
  # color = "white" draws thin white borders between countries.
  geom_polygon(data = world_alpha,
               aes(x = long, y = lat, group = group, fill = estimate_alpha_c2),
               color = "white", linewidth = 0.1) +
  # Diverging color scale centered at 0.
  # Blue  = negative alpha_c2 (warming reduces fire PM2.5 exposure),
  # red   = positive alpha_c2 (warming increases fire PM2.5 exposure),
  # gray  = countries with no regression data (NA).
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "darkred",
    midpoint = 0,
    name = "Alpha C2",
    na.value = "lightgray"
  ) +
  # coord_fixed(1.3) preserves geographic aspect ratio.
  coord_fixed(1.3) +
  theme_minimal() +
  labs(title = "Country-Level Alpha Coefficients (C2): FPM-GMT Regression",
       x = "Longitude", y = "Latitude") +
  theme(panel.grid = element_blank())


map_alpha_c2
ggsave(here("images/regression_alpha", "map_alpha_c2_park.png"),
       map_alpha_c2, width = 12, height = 6)

print("Alpha maps saved to images/regression_alpha/")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Line of Best Fit #########################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# For each selected country, plot all 23 observations (18 Park + 5 Pierce) with the
# OLS fitted line: exposure_percap = alpha_c1 + alpha_c2 * T_ps
#
# Data objects from source():
#   reg_data_combined -- 23 rows per country: T_ps (GMT regressor) and exposure_percap (response)
#   reg_coefs         -- alpha_c1 and alpha_c2 per country (define the fitted line)

countries_to_plot <- c("USA", "RUS", "AUS", "CHN", "IND", "ARG", "global")

for (iso in countries_to_plot) {

  # Filter to this country's 23 observations and create short label.
  # Labels are shown for all Pierce observations and only the Park 2010s points
  # (one per fire model). All other Park decades get NA and are unlabelled.
  df <- reg_data_combined %>%
    filter(country_code_iso3 == iso) %>%
    mutate(label = ifelse(
      !grepl("park", period_scenario) | grepl("2010s", period_scenario),
      str_remove(period_scenario, "exposure_percap_"),
      NA
    ))

  # Pull alpha_c1 (intercept) and alpha_c2 (slope) from reg_coefs
  coef_row  <- reg_coefs %>% filter(country_code_iso3 == iso)
  alpha_c1  <- coef_row$estimate_alpha_c1
  alpha_c2  <- coef_row$estimate_alpha_c2

  p <- ggplot(df, aes(x = T_ps, y = exposure_percap)) +
    # Each point is one observation: 18 Park (decade × fire model) +
    # 5 Pierce (baseline, 2040s RCP4.5/8.5, 2090s RCP4.5/8.5)
    geom_point(size = 3, color = "steelblue") +
    # Pierce labels: positioned above each point.
    geom_text(data = ~ filter(.x, !grepl("park", period_scenario)),
              aes(label = label), vjust = -0.8, size = 3) +
    # Park 2010s labels: positioned to the right of each point.
    geom_text(data = ~ filter(.x, grepl("2010s", period_scenario)),
              aes(label = label), hjust = -0.1, size = 3) +
    # OLS fitted line using alpha_c1 and alpha_c2 from reg_coefs.
    # geom_abline draws the exact same line as the stored regression coefficients,
    # ensuring consistency with the values used downstream in GIVE.
    geom_abline(intercept = alpha_c1, slope = alpha_c2,
                color = "darkred", linewidth = 0.8) +
    labs(
      title    = paste0(iso, ": Per-Capita Fire PM2.5 vs. GMT"),
      subtitle = paste0("alpha_c2 = ", round(alpha_c2, 4), " µg/m³/yr per °C"),
      x        = "GMT Anomaly relative to 1850-1900 (°C)",
      y        = "Per-Capita Fire PM2.5 Exposure (µg/m³/yr)"
    ) +
    theme_minimal() +
    theme(panel.grid.minor = element_blank())

  print(p)
  ggsave(here("images/regression_alpha", paste0("lof_", tolower(iso), "_park.png")),
         p, width = 7, height = 5, dpi = 300)
}


print("Line of best fit plots saved to images/regression_alpha/")

### THE END
