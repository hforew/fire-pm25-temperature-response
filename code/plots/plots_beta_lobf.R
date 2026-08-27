# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## PLOTS OF COUNTRY BETA COEFFICIENTS ##
##
## Goal: Visualise the main model's country-level FE regression outputs from
##   fpm_gmt_regress_country.R.
##   Produces: (1) histogram + world choropleth map of beta_c, combined into one
##   multiplot, plus the map (Panel B) saved standalone; (2) individual
##   country LOBF scatter plots; (3) 2x3 multi-country LOBF grid.
##
## Inputs:
##   code/regression_beta/regression/fpm_gmt_regress_country.R
##     (sourced directly; fits all 6 exclusion groups and stores each group's
##     results in coefs_by_group/results_by_group. This script pulls out the
##     "full" group as reg_coefs/reg_results; reg_data_combined is used as-is,
##     since it already holds the unfiltered full-model data before sourcing returns.)
##
## Outputs:
##   images/regression_beta/beta_country/fig3_multiplot_beta_c.png
##   images/regression_beta/beta_country/map_beta_c.png   (Panel B of fig3, standalone)
##   images/regression_beta/beta_country/country_lof/lof_{iso}_fe.png   (one per country in countries_to_plot)
##   images/regression_beta/beta_country/fig5_lof_multi.png
##
## Execution order:
##   files run before: fpm_gmt_regress_country.R   --> sourced directly; runs automatically
##   files run after: NA (no downstream files dependent on this one)
##
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
library(ggtext)
library(ggrepel)
library(patchwork)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Import #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# source() runs fpm_gmt_regress_country.R; reg_coefs/reg_results below are the "full"
# group pulled out of coefs_by_group/results_by_group (see extraction above). Key objects:
#   reg_coefs         -- one row per country; estimate_beta_c is the OLS slope:
#                        change in per-capita fire PM2.5 (µg/m^3/yr) per 1degC GMT.
#                        Positive = more fire PM with warming; negative = less.
#   reg_data_combined -- 43 rows per country: T_ps, exposure_percap, fire_model
#                        (18 Park factual + 18 Park counterfactual + 5 Pierce + 2 Zhao)
#   reg_results       -- nested list with tidied FE coefficients (intercept +
#                        4 model dummies) used to recover per-model fitted lines.

source(here("code/regression_beta/regression", "fpm_gmt_regress_country.R"))

# fpm_gmt_regress_country.R fits 6 exclusion groups in a loop and stores each group's
# results in coefs_by_group/results_by_group (keyed by group name) precisely because
# reg_coefs/reg_results are reassigned every iteration; pull out the "full" model here.
reg_coefs   <- coefs_by_group[["full"]]
reg_results <- results_by_group[["full"]]

head(reg_coefs)
str(reg_coefs)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Histogram of beta_c ####################################################
##
## Summarises the cross-country distribution of beta_c (OLS slope: change in per-capita
## fire PM2.5 per 1degC GMT increase). Three headline stats are computed and printed to
## console: total countries, share with positive beta_c (warming --> more fire PM),
## and share with p < .05 (statistically distinguishable from zero). These stats are
## also rendered as an annotation below the x-axis. The median is marked with a dotted
## vertical line and text label.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

reg_coefs_hist <- reg_coefs %>% filter(country_name != "global")  # exclude pooled global row from country-level histogram

n_total <- nrow(reg_coefs_hist)                                      # total countries in regression
n_pos   <- sum(reg_coefs_hist$estimate_beta_c > 0,    na.rm = TRUE)  # beta_c > 0: warming increases fire PM
n_sig   <- sum(reg_coefs_hist$p.value_beta_c  < 0.05, na.rm = TRUE)  # p < .05: slope sig. different from zero

# Format summary strings for console output; %d = integer, %% = literal percent sign
beta_c_pos <- sprintf("β(c) > 0: %d of %d (%d%%)", n_pos, n_total, round(100 * n_pos / n_total))
beta_c_sig <- sprintf("p < .05: %d of %d (%d%%)",  n_sig, n_total, round(100 * n_sig / n_total))

cat(beta_c_pos, "\n")
cat(beta_c_sig, "\n")

med_beta_c <- median(reg_coefs_hist$estimate_beta_c, na.rm = TRUE)  # median beta_c across countries

# Shared HTML label: used as subtitle on the console preview below, as annotate in multiplot
stats_label <- paste0("β<sub>c</sub> > 0: ", n_pos, " of ", n_total,
                      " (", round(100 * n_pos / n_total), "%) | ",
                      "p &lt; .05: ", n_sig, " of ", n_total,
                      " (", round(100 * n_sig / n_total), "%)")

# Base histogram: no annotation, no subtitle
hist_base <- ggplot(reg_coefs_hist, aes(x = estimate_beta_c)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white") +       # 50 bins; white borders separate bars
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +       # zero reference: no GMT effect
  geom_vline(xintercept = med_beta_c, linetype = "dotted", color = "black", linewidth = 0.8) +  # median marker
  annotate("text", x = med_beta_c, y = Inf,                              # y = Inf --> pin to top of panel
           label = paste0("Median = ", round(med_beta_c, 2)),
           vjust = 5, hjust = -0.1, size = 3.5, color = "black") +       # vjust pulls label down from top edge
  coord_cartesian(clip = "off") +                                         # allow drawing outside panel bounds
  labs(
    title = expression("Distribution of " ~ beta[c] ~ "across countries"),
    x     = expression(beta[c] ~ "(µg/m³/yr per °C)"),
    y     = "Number of countries"
  ) +
  theme_minimal() +
  theme(plot.margin = margin(b = 0.5, unit = "cm"))

# Console preview only (not saved): stats as subtitle; element_markdown renders the HTML superscript
hist_beta_c <- hist_base +
  labs(subtitle = stats_label) +
  theme(plot.subtitle = element_markdown(size = 10))

hist_beta_c


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

# Left join beta coefficients onto the map polygon vertices by ISO3 code.
# Every polygon vertex row gets the country's beta_c value.
# Countries not present in reg_coefs (no regression data) will have NA
# and will be rendered in the na.value color (lightgray) on the map.
world_beta <- world_map %>%
  left_join(reg_coefs, by = c("iso3" = "country_code_iso3"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Descriptive stats: beta_c distribution ################################
##
## Run before defining color palette breaks to see where the data actually falls.
## Key outputs: percentiles, share negative, share above candidate break thresholds.
## Use these to judge whether break values capture meaningful variation.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

thresholds <- c(0, 0.5, 1, 1.5, 2)

print_beta_stats <- function(x, label) {
  cat("\n--- beta_c descriptive stats:", label, "(n =", sum(!is.na(x)), "countries) ---\n")
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

print_beta_stats(reg_coefs$estimate_beta_c, "GLOBAL")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Map color palette  #################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

band_labels <- c("<0", "0-0.10", "0.10-0.25", "0.25-0.50", "0.50-1.00", ">1.00")

band_colors <- c(
  "<0"        = "azure2",  # azure2    -- warming reduces fire PM2.5
  "0-0.10"    = "#FFE566",  # light yellow -- low positive response
  "0.10-0.25" = "#C8A000",  # dark yellow  -- moderate-low
  "0.25-0.50" = "#CC5500",  # dark orange  -- moderate-high
  "0.50-1.00" = "#DD1111",  # red          -- high
  ">1.00"     = "#7A0000"   # dark red     -- extreme
)

# Bins a continuous beta value into the 6 labelled factor bands.
# right = FALSE --> intervals are [lo, hi) so boundary values fall in the upper band.
# include.lowest = TRUE closes the final [1.0, Inf] bin.
cut_beta <- function(x) {
  cut(x,
      breaks         = c(-Inf, 0, 0.1, 0.25, 0.5, 1.0, Inf),
      labels         = band_labels,
      right          = FALSE,
      include.lowest = TRUE)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Plot beta_c #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# beta_c is the OLS slope from the country-level regression with fire-model FEs:
#   exposure_percap = alpha_c_m + beta_c * T_ps + epsilon
#
# It is the key damage function parameter: the change in per-capita fire PM2.5
# exposure (µg/m³/yr) per 1degC increase in global mean temperature (GMT).
#
# Interpretation:
#   beta_c > 0  --> warming increases fire PM2.5 exposure in that country
#   beta_c < 0  --> warming decreases fire PM2.5 exposure in that country
#   beta_c = 0  --> no relationship between GMT and fire PM2.5 in that country
#
# This coefficient is used downstream in GIVE to translate GMT change scenarios
# into country-level fire PM2.5 exposure changes and associated mortality impacts.

map_beta_c <- world_beta %>%
  mutate(fill_band = cut_beta(estimate_beta_c)) %>%          # bin beta_c into 6 discrete colour bands
  ggplot(aes(x = long, y = lat, group = group, fill = fill_band)) +
  geom_polygon(color = "gray50", linewidth = 0.15) +          # gray50 borders separate country polygons
  scale_fill_manual(
    values       = band_colors,
    name         = NULL,                                      # no legend title
    drop         = FALSE,                                     # show all bands even if empty
    na.value     = "lightgray",                              # countries with no regression data --> gray
    na.translate = FALSE,                                     # exclude NA swatch from legend
    guide        = guide_legend(nrow = 1)                    # single-row legend
  ) +
  scale_x_continuous(expand = c(0, 0)) +                    # remove ggplot's default 5% x padding
  scale_y_continuous(expand = c(0, 0)) +                    # remove ggplot's default 5% y padding
  coord_fixed(1.3, xlim = c(-180, 180), ylim = c(-58, 83)) +  # fix aspect ratio; clip Antarctica
  theme_minimal() +
  labs(title = expression(beta[c] ~ "country estimates: Change in per-capita fire PM2.5 (µg/m³/yr) per 1°C GMT increase"),
       x = NULL, y = NULL) +
  theme(panel.grid       = element_blank(),                  # no graticule lines
        legend.position  = "bottom",                         # legend below map
        axis.text        = element_blank(),                  # hide lat/lon tick labels
        axis.ticks       = element_blank(),                  # hide tick marks
        plot.title       = element_text(size = 9.5),
        panel.background = element_blank())


map_beta_c

# theme_minimal() leaves plot.background transparent; add an opaque white
# background for the standalone save (patchwork gives the multiplot version
# its own background, so map_beta_c itself is left untouched for reuse there).
ggsave(here("images/regression_beta/beta_country", "map_beta_c.png"),
       map_beta_c + theme(plot.background = element_rect(fill = "white", colour = NA)),
       width = 8, height = 5, units = "in", dpi = 300)

print("Map (Panel B) saved to images/regression_beta/beta_country/")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Histogram + Map multiplot (Panel A / Panel B) ############################
##
## Combines hist_base (Panel A) and map_beta_c (Panel B) using patchwork.
## (hist_beta_c above is a console-only preview with its own subtitle; it is not
## reused here -- the multiplot restyles hist_base with its own stats annotation.)
## Layout: A left (40%), B right (60%); map given 2x the width.
## Individual panel titles are stripped and replaced with a single
## plot_annotation() title. A/B tags are added via tag_levels = "A".
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Strip individual titles; annotate carries stats text below x-axis in multiplot
hist_notitle <- hist_base +
  labs(title = NULL) +
  annotate("richtext", x = -Inf, y = -Inf,
           label = stats_label,
           hjust = .15, vjust = 3.6, size = 3.5,
           fill = NA, label.color = NA)                    # Panel A: stats as lower-left annotation
map_notitle  <- map_beta_c  + labs(title = NULL)

shared_annotation <- plot_annotation(
  title      = expression(beta[c] ~
                          "across countries: Change in per-capita fire PM2.5 (µg/m³/yr) per 1°C GMT increase"),
  tag_levels = "A"
)

# A left (40%), B right (60%); widths = c(2, 3) gives the 40/60 split
multiplot_beta_c <- (hist_notitle | map_notitle) +
  plot_layout(widths = c(2, 3)) +
  shared_annotation

print(multiplot_beta_c)
ggsave(here("images/regression_beta/beta_country", "fig3_multiplot_beta_c.png"),
       multiplot_beta_c, width = 8.5, height = 3.5, units = "in", dpi = 300)

print("Multiplot saved to images/regression_beta/beta_country/")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Line of Best Fit #########################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# For each selected country, plot all 43 observations (18 Park factual + 18 Park
# counterfactual + 5 Pierce + 2 Zhao) with five fitted lines -- one per fire model
# (CLASSIC, JULES, SSiB4, CESM, Zhao et al.) -- on raw exposure_percap.
#
# y-axis: exposure_percap (raw per-capita fire PM2.5 exposure).
# Fitted lines: one per fire model, each using the FE intercept for that model
# (alpha_classic + delta_m) and the shared slope beta_c. All five lines are
# parallel (common slope) but vertically offset by their model-level intercepts.
# Park factual and counterfactual runs of the same model share one FE group (see
# fpm_gmt_regress_country.R), so counterfactual points fall on the
# same line as their factual counterpart -- only the 5 base-model lines are drawn.
#
# Data objects from source():
#   reg_data_combined -- 43 rows per country: T_ps, exposure_percap, fire_model
#   reg_coefs         -- beta_c per country (shared slope across all five lines)
#   reg_results       -- tidied FE coefficients used to recover per-model intercepts

# Lookup table: maps each period_scenario to its fire model and warming trajectory.
# Used to set color (model) and shape (trajectory) aesthetics in the scatter plots.
scenario_lookup <- tibble(
  period_scenario = c(
    # Park factual (obsclim)
    "exposure_percap_park_classic_1960s_fpm", "exposure_percap_park_classic_1970s_fpm",
    "exposure_percap_park_classic_1980s_fpm", "exposure_percap_park_classic_1990s_fpm",
    "exposure_percap_park_classic_2000s_fpm", "exposure_percap_park_classic_2010s_fpm",
    "exposure_percap_park_jules_1960s_fpm",   "exposure_percap_park_jules_1970s_fpm",
    "exposure_percap_park_jules_1980s_fpm",   "exposure_percap_park_jules_1990s_fpm",
    "exposure_percap_park_jules_2000s_fpm",   "exposure_percap_park_jules_2010s_fpm",
    "exposure_percap_park_ssib4_1960s_fpm",   "exposure_percap_park_ssib4_1970s_fpm",
    "exposure_percap_park_ssib4_1980s_fpm",   "exposure_percap_park_ssib4_1990s_fpm",
    "exposure_percap_park_ssib4_2000s_fpm",   "exposure_percap_park_ssib4_2010s_fpm",
    # Park counterfactual (counterclim)
    "exposure_percap_park_classic_1960s_fpm_counter", "exposure_percap_park_classic_1970s_fpm_counter",
    "exposure_percap_park_classic_1980s_fpm_counter", "exposure_percap_park_classic_1990s_fpm_counter",
    "exposure_percap_park_classic_2000s_fpm_counter", "exposure_percap_park_classic_2010s_fpm_counter",
    "exposure_percap_park_jules_1960s_fpm_counter",   "exposure_percap_park_jules_1970s_fpm_counter",
    "exposure_percap_park_jules_1980s_fpm_counter",   "exposure_percap_park_jules_1990s_fpm_counter",
    "exposure_percap_park_jules_2000s_fpm_counter",   "exposure_percap_park_jules_2010s_fpm_counter",
    "exposure_percap_park_ssib4_1960s_fpm_counter",   "exposure_percap_park_ssib4_1970s_fpm_counter",
    "exposure_percap_park_ssib4_1980s_fpm_counter",   "exposure_percap_park_ssib4_1990s_fpm_counter",
    "exposure_percap_park_ssib4_2000s_fpm_counter",   "exposure_percap_park_ssib4_2010s_fpm_counter",
    # Pierce projections
    "exposure_percap_fpm_2000",
    "exposure_percap_fpm_2050_45",            "exposure_percap_fpm_2050_85",
    "exposure_percap_fpm_2100_45",            "exposure_percap_fpm_2100_85",
    # Zhao projections
    "exposure_percap_fpm_2095_SSP245_Zhao",   "exposure_percap_fpm_2095_SSP585_Zhao"
  ),
  model = c(
    rep("CLASSIC", 6), rep("JULES", 6), rep("SSiB4", 6),  # factual
    rep("CLASSIC", 6), rep("JULES", 6), rep("SSiB4", 6),  # counterfactual -- same model/color as factual
    "CESM", "CESM", "CESM", "CESM", "CESM",                # Pierce
    "Zhao et al.", "Zhao et al."                            # Zhao
  ),
  trajectory = c(
    rep("Historical", 18),                    # Park factual
    rep("Historical - counterfactual", 18),   # Park counterfactual
    "Historical", "RCP4.5", "RCP8.5", "RCP4.5", "RCP8.5",  # Pierce
    "SSP2-4.5", "SSP5-8.5"                                  # Zhao
  )
)

model_colors <- c(
  "CESM"    = "#E69F00",  # orange
  "JULES"   = "#56B4E9",  # blue
  "SSiB4"   = "#009E73",  # green
  "CLASSIC" = "#CC79A7",  # pink
  "Zhao et al." = "#000000"   # black
)

trajectory_shapes <- c(
  "Historical"                 = 16,  # filled circle
  "Historical - counterfactual" = 10, # circle with plus -- Park counterclim runs
  "RCP4.5"   = 17,
  "RCP8.5"   = 15,
  "SSP2-4.5"   = 18,
  "SSP5-8.5"   = 8
)

countries_to_plot <- c("global", "USA","DEU", "RUS", "AUS", "CHN", "IND", "ARG", "BRA", "ETH", "NGA")

country_lof_dir <- here("images/regression_beta/beta_country/country_lof")
if (!dir.exists(country_lof_dir)) dir.create(country_lof_dir, recursive = TRUE)

for (iso in countries_to_plot) {

  # Filter to this country's 43 observations, join model/trajectory lookup, and
  # create point labels. Labels shown for Pierce/Zhao points and Park 2010s only
  # (factual and counterfactual); all other Park decades are unlabelled to avoid overplotting.
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
          "exposure_percap_park_classic_2010s_fpm_counter",
          "exposure_percap_park_jules_2010s_fpm_counter",
          "exposure_percap_park_ssib4_2010s_fpm_counter"
        ) ~ "Park et al. (counter)",
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

  # Pull beta_c (slope) and its standard error from reg_coefs
  coef_row <- reg_coefs %>% filter(country_code_iso3 == iso)
  beta_c   <- coef_row$estimate_beta_c
  se_beta_c <- coef_row$std.error_beta_c

  # Pull R² from the glanced model-fit statistics
  r_squared <- reg_results %>%
    filter(country_code_iso3 == iso) %>%
    unnest(glanced) %>%
    pull(r.squared)

  # Extract per-model FE intercepts from the tidy coefficient table.
  # (Intercept) = classic baseline; other models add their delta offset.
  # NB: Park factual and counterfactual runs of the same model are pooled into one
  # FE group (see fpm_gmt_regress_country.R) -- there is no separate
  # *_counter dummy term, so only the 5 base-model lines are built below.
  tidied_c      <- reg_results %>% filter(country_code_iso3 == iso) %>% unnest(tidied)
  alpha_classic <- tidied_c %>% filter(term == "(Intercept)")     %>% pull(estimate)  # classic baseline intercept
  delta_jules   <- tidied_c %>% filter(term == "fire_modeljules") %>% pull(estimate)  # jules offset from classic
  delta_ssib4   <- tidied_c %>% filter(term == "fire_modelssib4") %>% pull(estimate)  # ssib4 offset from classic
  delta_CESM    <- tidied_c %>% filter(term == "fire_modelCESM")  %>% pull(estimate)  # CESM offset from classic
  delta_Zhao    <- tidied_c %>% filter(term == "fire_modelZhao")  %>% pull(estimate)  # Zhao offset from classic

  lines_df <- tibble(
    model       = factor(c("CLASSIC", "JULES", "SSiB4", "CESM", "Zhao et al."),
                         levels = names(model_colors)),
    intercept_m = c(alpha_classic,                             # classic: base intercept
                    alpha_classic + delta_jules,               # jules:   base + model offset
                    alpha_classic + delta_ssib4,               # ssib4:   base + model offset
                    alpha_classic + delta_CESM,                # CESM:    base + model offset
                    alpha_classic + delta_Zhao),               # Zhao:    base + model offset
    slope       = beta_c                                       # shared slope across all models
  )

  p <- ggplot(df, aes(x = T_ps, y = exposure_percap, color = model, shape = trajectory)) +
    geom_point(size = 3) +
    geom_text_repel(data = ~ filter(.x, !is.na(label)),   # only labelled points
                    aes(label = label),
                    size = 3, box.padding = 0.4, max.overlaps = Inf,
                    show.legend = FALSE) +
    # One fitted line per fire model: FE intercept for that model plus shared
    # slope beta_c. Parallel lines reflect the FE assumption of a common GMT
    # sensitivity with model-specific mean-level offsets.
    geom_abline(data = lines_df,
                aes(intercept = intercept_m, slope = slope, color = model),
                linetype = "dotted", linewidth = 0.5, show.legend = FALSE) +  # suppress from legend; points carry color
    scale_color_manual(values = model_colors,      name = "Model") +
    scale_shape_manual(values = trajectory_shapes, name = "Trajectory") +
    labs(
      title    = paste0(iso, ": Per-Capita Fire PM2.5 Exposure vs. GMT"),
      subtitle = bquote(beta[c] == .(round(beta_c, 2)) %+-% .(round(se_beta_c, 2)) ~
                          "µg/m³/yr per °C  |  R² =" ~ .(round(r_squared, 2))),
      x        = "GMT Anomaly relative to 1850-1900 (°C)",
      y        = "Per-Capita Fire PM2.5 Exposure (µg/m³/yr)"
    ) +
    theme_minimal() +
    theme(panel.grid.minor  = element_blank(),
          plot.subtitle     = element_text(size = 12))

  print(p)
  ggsave(file.path(country_lof_dir, paste0("lof_", tolower(iso), "_fe.png")),
         p, width = 7, height = 5, dpi = 300)
}


print("Line of best fit plots saved to images/regression_beta/beta_country/country_lof/")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Line of Best Fit: global + 5 most populous countries ######################
##
## Six-panel faceted LOBF plot: global aggregate plus CHN, IND, USA, IDN, PAK
## (the five most populous countries).
## Layout: 2 rows x 3 columns, single legend.
## Strip label: country code + beta_c + SE.
## y-axis: exposure_percap (raw per-capita fire PM2.5 exposure).
## Five parallel dotted fitted lines per panel: one per fire model (CLASSIC, JULES,
## SSiB4, CESM, Zhao et al.), each using its FE intercept (alpha_classic + delta_m) and the
## shared slope beta_c. Park counterfactual points share their factual model's FE
## group, so they fall on the same line rather than getting a separate one.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

multiplot_countries <- c("global", "CHN", "IND", "USA", "IDN", "PAK")

# Builds a faceted LOBF plot for a given vector of ISO3 codes.
# Returns a ggplot object; caller saves to disk.
build_lobf_multi <- function(countries) {

  df <- purrr::map_dfr(countries, function(iso) {   # one 43-row block per country, row-bound
    coef_row  <- reg_coefs %>% filter(country_code_iso3 == iso)
    beta_c    <- coef_row$estimate_beta_c    # slope for strip label text
    se_beta_c <- coef_row$std.error_beta_c   # SE for strip label text

    reg_data_combined %>%
      filter(country_code_iso3 == iso) %>%
      left_join(scenario_lookup, by = "period_scenario") %>%
      mutate(
        model       = factor(model,      levels = names(model_colors)),      # enforce color order
        trajectory  = factor(trajectory, levels = names(trajectory_shapes)), # enforce shape order
        strip_label = paste0(iso, "  |  β<sub>c</sub> = ", round(beta_c, 2),
                             " ± ", round(se_beta_c, 2),
                             " µg/m³/yr per °C")                            # HTML rendered by element_markdown
      )
  })

  # Preserve the supplied country order in facet panels
  strip_levels <- df %>%
    distinct(country_code_iso3, strip_label) %>%
    arrange(match(country_code_iso3, countries)) %>%
    pull(strip_label)

  df <- df %>%
    mutate(strip_label = factor(strip_label, levels = strip_levels))

  # One row per country x fire model: FE intercept + shared slope for geom_abline.
  # Park factual/counterfactual runs of a model are pooled into one FE group (see
  # fpm_gmt_regress_country.R), so only the 5 base-model lines exist.
  lines <- purrr::map_dfr(countries, function(iso) {   # 5 rows per country (one per model)
    tidied_c      <- reg_results %>% filter(country_code_iso3 == iso) %>% unnest(tidied)
    beta_c_val    <- tidied_c %>% filter(term == "T_ps")            %>% pull(estimate)
    alpha_classic <- tidied_c %>% filter(term == "(Intercept)")     %>% pull(estimate)
    delta_jules   <- tidied_c %>% filter(term == "fire_modeljules") %>% pull(estimate)
    delta_ssib4   <- tidied_c %>% filter(term == "fire_modelssib4") %>% pull(estimate)
    delta_CESM    <- tidied_c %>% filter(term == "fire_modelCESM")  %>% pull(estimate)
    delta_Zhao    <- tidied_c %>% filter(term == "fire_modelZhao")  %>% pull(estimate)

    strip_lbl <- df %>%
      filter(country_code_iso3 == iso) %>%
      pull(strip_label) %>%
      first()                                                                # all rows share the same label

    tibble(
      country_code_iso3 = iso,
      strip_label       = strip_lbl,                                         # must match df for facet subsetting
      model             = factor(c("CLASSIC", "JULES", "SSiB4", "CESM", "Zhao et al."),
                                 levels = names(model_colors)),
      intercept_m       = c(alpha_classic,
                            alpha_classic + delta_jules,
                            alpha_classic + delta_ssib4,
                            alpha_classic + delta_CESM,
                            alpha_classic + delta_Zhao),
      slope             = beta_c_val
    )
  }) %>%
    mutate(strip_label = factor(strip_label, levels = strip_levels))         # align with df for faceting

  ggplot(df, aes(x = T_ps, y = exposure_percap, color = model, shape = trajectory)) +
    geom_point(size = 2) +
    geom_abline(data = lines,
                aes(intercept = intercept_m, slope = slope, color = model),
                linetype = "dotted", linewidth = 0.5, show.legend = FALSE) +  # suppress from legend; points carry color
    scale_color_manual(values = model_colors,      name = "Model") +
    scale_shape_manual(values = trajectory_shapes, name = "Trajectory") +
    guides(shape = guide_legend(order = 1),   # trajectory first
           color = guide_legend(order = 2)) + # model second
    facet_wrap(~ strip_label, ncol = 3) +     # 3 cols --> 2 rows for 6 countries
    labs(
      title = "Per-Capita Fire PM2.5 Exposure vs. GMT",
      x     = "GMT Anomaly relative to 1850-1900 (°C)",
      y     = "Per-Capita Fire PM2.5 Exposure (µg/m³/yr)"
    ) +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      strip.text       = element_markdown(size = 9),
      legend.position  = "bottom",
      legend.box       = "vertical",
      legend.spacing.y = unit(-0.3, "cm"),
      plot.margin      = margin(t = 0.5, r = 0.5, b = 0.1, l = 0.5, unit = "cm")
    )
}

p_multi <- build_lobf_multi(multiplot_countries)

print(p_multi)
ggsave(here("images/regression_beta/beta_country", "fig5_lof_multi.png"),
       p_multi, width = 8.5, height = 5, dpi = 300)

print("Multi-country grid saved to images/regression_beta/beta_country/")

############ THE END  ############################

