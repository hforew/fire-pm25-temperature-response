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
library(ggrepel)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Import #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Import country-level regression coefficients produced by fpm_gmt_relation_final.R.
# Each row is one country. Key columns:
#   country_code_iso3   -- ISO 3166-1 alpha-3 country code (used for map joining)
#   country_name        -- full country name
#   estimate_alpha_c1   -- OLS intercept: fitted exposure at T_ps = 0 (extrapolated
#                          outside the data range; not physically meaningful)
#   estimate_alpha_c2   -- OLS slope: change in per-capita fire PM2.5 exposure
#                          (µg/m^3/yr) per 1°C increase in GMT. This is the key
#                          damage function parameter used downstream in GIVE.
#                          Positive = more fire PM with warming; negative = less.

source(here("code/regression_alpha", "fpm_gmt_regression_zhao_pierce.R")) # regression file and load all objects

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
ggsave(here("images/regression_alpha/alpha_zhao_pierce", "hist_alpha_c2_zhao_pierce.png"), plot = hist_alpha_c2, width = 8, height = 5, dpi = 300)


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
############ Map color palette (Smokey Bear fire danger scale) #################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

plot_colors <- c(
  low       = "#7FD67F",  # light green   (Low)
  moderate  = "#56C8FF",  # bright sky blue (Moderate)
  high      = "#FFD700",  # gold yellow   (High)
  very_high = "#FF8C00",  # warm orange   (Very High)
  extreme   = "#CC1111"   # vivid red     (Extreme)
)

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
  scale_fill_stepsn(
    colors   = plot_colors,
    n.breaks = 4,
    name     = "Alpha C1",
    na.value = "lightgray"
  ) +
  # coord_fixed(1.3) preserves geographic aspect ratio (lat/lon are not
  # equal-distance units; 1.3 approximates a Mercator-like stretch).
  coord_fixed(1.3) +
  theme_minimal() +
  labs(title = expression("Country-Level Intercept Coefficient (" * alpha[c1] * "): FPM-GMT Regression"),
       x = "Longitude", y = "Latitude") +
  theme(panel.grid = element_blank())


map_alpha_c1
ggsave(here("images/regression_alpha/alpha_zhao_pierce", "map_alpha_c1_zhao_pierce.png"),
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
  scale_fill_stepsn(
    colors   = plot_colors,
    n.breaks = 4,
    name     = "Alpha C2",
    na.value = "lightgray"
  ) +
  # coord_fixed(1.3) preserves geographic aspect ratio.
  coord_fixed(1.3) +
  theme_minimal() +
  labs(title = expression("Country-Level Slope Coefficient (" * alpha[c2] * "): FPM-GMT Regression"),
       x = "Longitude", y = "Latitude") +
  theme(panel.grid = element_blank())


map_alpha_c2
ggsave(here("images/regression_alpha/alpha_zhao_pierce", "map_alpha_c2_zhao_pierce.png"),
       map_alpha_c2, width = 12, height = 6)

print("Alpha maps saved to images/regression_alpha/alpha_zhao_pierce/")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Line of Best Fit #########################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# For each selected country, plot all 7 observations (5 Pierce + 2 Zhao) with the
# OLS fitted line: exposure_percap = alpha_c1 + alpha_c2 * T_ps
#
# Data objects from source():
#   reg_data_long -- 7 rows per country: T_ps (GMT regressor) and exposure_percap (response)
#   reg_coefs     -- alpha_c1 and alpha_c2 per country (define the fitted line)

# Lookup table: maps each period_scenario to its fire model and warming trajectory.
# Used to set color (model) and shape (trajectory) aesthetics in the scatter plots.
scenario_lookup <- tibble(
  period_scenario = c(
    "exposure_percap_fpm_2000",
    "exposure_percap_fpm_2050_45",            "exposure_percap_fpm_2050_85",
    "exposure_percap_fpm_2100_45",            "exposure_percap_fpm_2100_85",
    "exposure_percap_fpm_2095_SSP245_Zhao",   "exposure_percap_fpm_2095_SSP585_Zhao"
  ),
  model = c(
    "CESM", "CESM", "CESM", "CESM", "CESM",
    "Multi-model", "Multi-model"
  ),
  trajectory = c(
    "Historical", "RCP4.5", "RCP8.5", "RCP4.5", "RCP8.5",
    "SSP2-4.5", "SSP5-8.5"
  )
)

model_colors <- c(
  "CESM"        = "#E69F00",
  "Multi-model" = "#000000"
)

trajectory_shapes <- c(
  "Historical" = 16,
  "RCP4.5"     = 17,
  "RCP8.5"     = 15,
  "SSP2-4.5"   = 18,
  "SSP5-8.5"   = 8
)

countries_to_plot <- c("global", "USA","DEU", "RUS", "AUS", "CHN", "IND", "ARG", "BRA", "ETH", "NGA")

for (iso in countries_to_plot) {

  # Filter to this country's 7 observations, join model/trajectory lookup, and
  # create point labels for all Pierce and Zhao points.
  df <- reg_data_long %>%
    filter(country_code_iso3 == iso) %>%
    left_join(scenario_lookup, by = "period_scenario") %>%
    mutate(
      model      = factor(model,      levels = names(model_colors)),
      trajectory = factor(trajectory, levels = names(trajectory_shapes)),
      label = case_when(
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
    # Each point is one observation: 5 Pierce (baseline, 2040s RCP4.5/8.5, 2090s RCP4.5/8.5)
    # + 2 Zhao (~2095 SSP2-4.5/SSP5-8.5)
    geom_point(size = 3) +
    # Labels for all labelled points (Pierce, Zhao).
    # ggrepel automatically nudges overlapping labels apart and draws
    # a line back to the point when a label is pushed away.
    geom_text_repel(data = ~ filter(.x, !is.na(label)),
                    aes(label = label),
                    size = 3, box.padding = 0.4, max.overlaps = Inf,
                    show.legend = FALSE) +
    # OLS fitted line using alpha_c1 and alpha_c2 from reg_coefs.
    # geom_abline draws the exact same line as the stored regression coefficients,
    # ensuring consistency with the values used downstream in GIVE.
    geom_abline(intercept = alpha_c1, slope = alpha_c2,
                color = "darkred", linewidth = 0.8) +
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
    theme(panel.grid.minor = element_blank())

  print(p)
  ggsave(here("images/regression_alpha/alpha_zhao_pierce", paste0("lof_", tolower(iso), "_zhao_pierce.png")),
         p, width = 7, height = 5, dpi = 300)
}


print("Line of best fit plots saved to images/regression_alpha/alpha_zhao_pierce/")

### THE END
