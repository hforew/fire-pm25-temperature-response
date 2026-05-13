# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## PLOTS OF FIRE PM ALL DATASETS -- Park, Pierce, Zhao ##
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls())

## this code plots fire PM globally and for the USA for all datasets (park, pierce, zhao)
## desriptive stats are also provided. 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Fire PM column reference (25 cols; excludes _hi and _chg variants)
##
## Zhao (Multi-model)
##   fpm_2095_SSP245_Zhao      Multi-model   SSP2-4.5
##   fpm_2095_SSP585_Zhao      Multi-model   SSP5-8.5
##
## Pierce (CESM)
##   fpm_2000                  CESM          Historical
##   fpm_2050_45               CESM          SSP1-4.5
##   fpm_2050_85               CESM          SSP3-8.5
##   fpm_2100_45               CESM          SSP1-4.5
##   fpm_2100_85               CESM          SSP3-8.5
##
## Park (Historical)
##   park_classic_1960s_fpm    CLASSIC       Historical
##   park_classic_1970s_fpm    CLASSIC       Historical
##   park_classic_1980s_fpm    CLASSIC       Historical
##   park_classic_1990s_fpm    CLASSIC       Historical
##   park_classic_2000s_fpm    CLASSIC       Historical
##   park_classic_2010s_fpm    CLASSIC       Historical
##   park_jules_1960s_fpm      JULES         Historical
##   park_jules_1970s_fpm      JULES         Historical
##   park_jules_1980s_fpm      JULES         Historical
##   park_jules_1990s_fpm      JULES         Historical
##   park_jules_2000s_fpm      JULES         Historical
##   park_jules_2010s_fpm      JULES         Historical
##   park_ssib4_1960s_fpm      SSiB4         Historical
##   park_ssib4_1970s_fpm      SSiB4         Historical
##   park_ssib4_1980s_fpm      SSiB4         Historical
##   park_ssib4_1990s_fpm      SSiB4         Historical
##   park_ssib4_2000s_fpm      SSiB4         Historical
##   park_ssib4_2010s_fpm      SSiB4         Historical
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Packages ##########################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(here)
library(tidyverse)
library(ggplot2)
library(maps)
library(countrycode)
library(ggrepel)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Import ############################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

df <- read_csv(here("output", "pop_pm_country_death_final.csv"))

gmt_pierce  <- read_csv(here("output", "gmt_periods_pi.csv"))
gmt_park <- read_csv(here("output", "gmt_park_decades.csv"))
gmt_zhao <- read_csv(here("output", "gmt_zhao_mimi.csv"))

fpm_cols <- c(
  "fpm_2095_SSP245_Zhao", "fpm_2095_SSP585_Zhao",
  "fpm_2000", "fpm_2050_45", "fpm_2050_85", "fpm_2100_45", "fpm_2100_85",
  "park_classic_1960s_fpm", "park_classic_1970s_fpm", "park_classic_1980s_fpm",
  "park_classic_1990s_fpm", "park_classic_2000s_fpm", "park_classic_2010s_fpm",
  "park_jules_1960s_fpm",   "park_jules_1970s_fpm",   "park_jules_1980s_fpm",
  "park_jules_1990s_fpm",   "park_jules_2000s_fpm",   "park_jules_2010s_fpm",
  "park_ssib4_1960s_fpm",   "park_ssib4_1970s_fpm",   "park_ssib4_1980s_fpm",
  "park_ssib4_1990s_fpm",   "park_ssib4_2000s_fpm",   "park_ssib4_2010s_fpm"
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Descriptive Stats #################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

gmt_baseline  <- mean(c(gmt_pierce$mean_gmt_45[gmt_pierce$period == "2006-2010"],
                        gmt_pierce$mean_gmt_85[gmt_pierce$period == "2006-2010"]))
gmt_2040s_45  <- gmt_pierce$mean_gmt_45[gmt_pierce$period == "2041-2050"]
gmt_2040s_85  <- gmt_pierce$mean_gmt_85[gmt_pierce$period == "2041-2050"]
gmt_2090s_45  <- gmt_pierce$mean_gmt_45[gmt_pierce$period == "2091-2100"]
gmt_2090s_85  <- gmt_pierce$mean_gmt_85[gmt_pierce$period == "2091-2100"]
gmt_2090s_245 <- gmt_zhao$mean_gmt_pi[gmt_zhao$scenario == "SSP245"]
gmt_2090s_585 <- gmt_zhao$mean_gmt_pi[gmt_zhao$scenario == "SSP585"]
gmt_1960s     <- gmt_park$mean_gmt_pi[gmt_park$park_year == 1965]
gmt_1970s     <- gmt_park$mean_gmt_pi[gmt_park$park_year == 1975]
gmt_1980s     <- gmt_park$mean_gmt_pi[gmt_park$park_year == 1985]
gmt_1990s     <- gmt_park$mean_gmt_pi[gmt_park$park_year == 1995]
gmt_2000s     <- gmt_park$mean_gmt_pi[gmt_park$park_year == 2005]
gmt_2010s     <- gmt_park$mean_gmt_pi[gmt_park$park_year == 2015]

fpm_lookup <- tibble(
  col = fpm_cols,
  model = c(
    rep("Multi", 2), rep("CESM", 5),
    rep("CLASSIC", 6), rep("JULES", 6), rep("SSiB4", 6)
  ),
  trajectory = c(
    "SSP2-4.5", "SSP5-8.5",
    "Historical", "SSP1-4.5", "SSP3-8.5", "SSP1-4.5", "SSP3-8.5",
    rep("Historical", 18)
  ),
  year = c(
    "2095", "2095",
    "2000s", "2040s", "2040s", "2090s", "2090s",
    rep(c("1960s", "1970s", "1980s", "1990s", "2000s", "2010s"), 3)
  ),
  gmt = c(
    gmt_2090s_245, gmt_2090s_585,
    gmt_baseline, gmt_2040s_45, gmt_2040s_85, gmt_2090s_45, gmt_2090s_85,
    rep(c(gmt_1960s, gmt_1970s, gmt_1980s, gmt_1990s, gmt_2000s, gmt_2010s), 3)
  )
)

make_fpm_table <- function(data, cols) {
  bind_rows(lapply(cols, function(col) {
    x   <- data[[col]]
    pct <- quantile(x, probs = c(0.05, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99), na.rm = TRUE)
    tibble(
      col  = col,
      Min  = round(min(x,  na.rm = TRUE), 3),
      Max  = round(max(x,  na.rm = TRUE), 3),
      Mean = round(mean(x, na.rm = TRUE), 3),
      SD   = round(sd(x,   na.rm = TRUE), 3),
      `5%`  = round(pct[[1]], 3),
      `25%` = round(pct[[2]], 3),
      `50%` = round(pct[[3]], 3),
      `75%` = round(pct[[4]], 3),
      `90%` = round(pct[[5]], 3),
      `95%` = round(pct[[6]], 3),
      `99%` = round(pct[[7]], 3)
    )
  })) |>
    left_join(fpm_lookup, by = "col") |>
    rename(Model = model, Trajectory = trajectory, Year = year, GMT = gmt) |>
    select(Model, Trajectory, Year, GMT, everything(), -col)
}

#### Global ###

cat("=== GLOBAL ===\n")
print(make_fpm_table(df, fpm_cols), n = 25)

#### USA ###

df_usa   <- df |> filter(country_code_iso3 == "USA")
df_conus <- df |> filter(country_code_iso3 == "USA",
                          lon >= -125, lon <= -66,
                          lat >=   24, lat <=  50)

cat("\n=== USA ===\n")
print(make_fpm_table(df_usa, fpm_cols), n = 25)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Map color palette #################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

band_labels_fpm <- c("<0.1", "0.1-0.5", "0.5-1.5", "1.5-5", "5-10", ">10")

band_colors_fpm <- c(
  "<0.1"    = "#7FD67F",  # light green  -- very low
  "0.1-0.5" = "#56C8FF",  # sky blue     -- low
  "0.5-1.5" = "#C8A000",  # dark yellow  -- moderate-low
  "1.5-5"   = "#CC5500",  # dark orange  -- moderate-high
  "5-10"    = "#DD1111",  # red          -- high
  ">10"     = "#7A0000"   # dark red     -- extreme
)

cut_fpm <- function(x) {
  cut(x,
      breaks         = c(-Inf, 0.1, 0.5, 1.5, 5, 10, Inf),
      labels         = band_labels_fpm,
      right          = FALSE,
      include.lowest = TRUE)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Maps: Park multiplot (CLASSIC x JULES x SSiB4, 6 decades) ########
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Subset to the 18 Park columns (all cols starting with "park_").
park_cols <- fpm_cols[grepl("^park_", fpm_cols)]

# Pivot to long format: one row per grid cell x scenario combination.
# model and decade are extracted from the column name (e.g. "park_jules_1980s_fpm"
# --> model = "JULES", decade = "1980s"). Both are set as ordered factors so
# facet_grid() arranges panels in the correct column/row order.
# cut_fpm() bins the continuous fpm value into the 6 colour bands defined above.
park_decade_labels <- c(
  "1960s" = paste0("1960s: T=", round(gmt_1960s, 2)),
  "1970s" = paste0("1970s: T=", round(gmt_1970s, 2)),
  "1980s" = paste0("1980s: T=", round(gmt_1980s, 2)),
  "1990s" = paste0("1990s: T=", round(gmt_1990s, 2)),
  "2000s" = paste0("2000s: T=", round(gmt_2000s, 2)),
  "2010s" = paste0("2010s: T=", round(gmt_2010s, 2))
)

df_park_long <- df |>
  select(lon, lat, all_of(park_cols)) |>
  pivot_longer(cols = all_of(park_cols), names_to = "col", values_to = "fpm") |>
  mutate(
    model  = case_when(
      str_detect(col, "classic") ~ "CLASSIC",
      str_detect(col, "jules")   ~ "JULES",
      str_detect(col, "ssib4")   ~ "SSiB4"
    ),
    decade    = park_decade_labels[str_extract(col, "\\d{4}s")],
    model     = factor(model,  levels = c("CLASSIC", "JULES", "SSiB4")),
    decade    = factor(decade, levels = unname(park_decade_labels)),
    fill_band = cut_fpm(fpm)
  )

# 3 x 6 faceted raster map: columns = fire model, rows = decade.
# geom_tile() with explicit 0.5-degree width/height matches the grid resolution
# and avoids gaps that geom_raster() can produce near coastlines.
# drop = FALSE forces all 6 colour bands into the legend even if a band is
# empty in a given panel, keeping the legend identical across all 18 panels.
# coord_fixed(1.3) preserves the geographic aspect ratio.
p_park <- ggplot(df_park_long, aes(x = lon, y = lat, fill = fill_band)) +
  geom_tile(width = 0.5, height = 0.5) +
  scale_fill_manual(
    values   = band_colors_fpm,
    name     = "Fire PM2.5 (µg/m³)",
    drop     = FALSE,
    na.value = "lightgray"
  ) +
  coord_fixed(1.3) +
  facet_grid(decade ~ model) +
  labs(
    title = "Park et al. Historical Fire PM2.5 by Model and Decade",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid        = element_blank(),
    axis.text         = element_blank(),
    axis.ticks        = element_blank(),
    legend.position   = "bottom",
    legend.title      = element_text(size = 14),
    legend.text       = element_text(size = 13),
    legend.box.margin = margin(t = 0, r = 0, b = -8, l = 0),
    strip.text        = element_text(size = 14),
    plot.title        = element_text(size = 18, hjust = 0.5),
    plot.margin       = margin(t = 4, r = 4, b = 4, l = 4)
  ) +
  guides(fill = guide_legend(nrow = 1))  # single-row legend at bottom

p_park
ggsave(here("images", "pm", "map_park_multiplot.png"),
       p_park, width = 8.5, height = 11, dpi = 300)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Maps: Pierce + Zhao multiplot (CESM x Multi, 5 rows) #############
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Layout: 2 columns x 5 rows.
#   Col 1 (Pierce/CESM): all 5 rows -- Historical, SSP1-4.5/SSP3-8.5 2040s, SSP1-4.5/SSP3-8.5 2090s
#   Col 2 (Zhao/Multi):  rows 4-5 only -- SSP2-4.5 and SSP5-8.5 (~2090s)
#   Rows 1-3 of col 2 will render as empty panels; this is intentional.
# Rows 4-5 share labels across both columns because Pierce SSP1-4.5/SSP3-8.5 2090s
# and Zhao SSP2-4.5/SSP5-8.5 are roughly comparable future periods.

pierce_zhao_cols <- c(
  "fpm_2000", "fpm_2050_45", "fpm_2050_85", "fpm_2100_45", "fpm_2100_85",
  "fpm_2095_SSP245_Zhao", "fpm_2095_SSP585_Zhao"
)

# Row labels embed the GMT scalar values so the strip shows the actual
# temperature anomaly alongside the period. Rows 4-5 show both Pierce and
# Zhao GMT values separated by " - " since both datasets occupy that row.
pz_row_labels <- c(
  paste0("2000s: T=", round(gmt_baseline,  2)),
  paste0("2040s: T=", round(gmt_2040s_45, 2)),
  paste0("2040s: T=", round(gmt_2040s_85, 2)),
  paste0("2090s: T=", round(gmt_2090s_45, 2), "-", round(gmt_2090s_245, 2)),
  paste0("2090s: T=", round(gmt_2090s_85, 2), "-", round(gmt_2090s_585, 2))
)

pz_meta <- tibble(
  col_name  = pierce_zhao_cols,
  col_label = c(rep("Pierce (CESM)", 5), rep("Zhao (Multi)", 2)),
  row_label = c(
    pz_row_labels,
    pz_row_labels[4],  # Zhao SSP2-4.5 --> same row as Pierce SSP1-4.5 2090s
    pz_row_labels[5]   # Zhao SSP5-8.5 --> same row as Pierce SSP3-8.5 2090s
  )
)

# Pivot to long, join labels, bin into colour bands.
# facet_grid() creates all row x col combinations; the 3 missing Zhao panels
# (rows 1-3) have no data and render as empty lightgray tiles.
df_pz_long <- df |>
  select(lon, lat, all_of(pierce_zhao_cols)) |>
  pivot_longer(cols = all_of(pierce_zhao_cols), names_to = "col_name", values_to = "fpm") |>
  left_join(pz_meta, by = "col_name") |>
  mutate(
    col_label = factor(col_label, levels = c("Pierce (CESM)", "Zhao (Multi)")),
    row_label = factor(row_label, levels = pz_row_labels),
    fill_band = cut_fpm(fpm)
  )

p_pz <- ggplot(df_pz_long, aes(x = lon, y = lat, fill = fill_band)) +
  geom_tile(width = 0.5, height = 0.5) +
  scale_fill_manual(
    values   = band_colors_fpm,
    name     = "Fire PM2.5 (µg/m³)",
    drop     = FALSE,
    na.value = "lightgray"
  ) +
  coord_fixed(1.3) +
  facet_grid(row_label ~ col_label) +
  labs(
    title = "Pierce et al. (CESM) and Zhao et al. (Multi-model) Fire PM2.5",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid        = element_blank(),
    axis.text         = element_blank(),
    axis.ticks        = element_blank(),
    legend.position   = "bottom",
    legend.title      = element_text(size = 14),
    legend.text       = element_text(size = 13),
    legend.box.margin = margin(t = 0, r = 0, b = -8, l = 0),
    strip.text        = element_text(size = 14),
    plot.title        = element_text(size = 18, hjust = 0.5),
    plot.margin       = margin(t = 4, r = 4, b = 4, l = 4)
  ) +
  guides(fill = guide_legend(nrow = 1))

p_pz
ggsave(here("images", "pm", "map_pierce_zhao_multiplot.png"),
       p_pz, width = 8.5, height = 11, dpi = 300)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Maps: Park multiplot -- CONUS #####################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# CONUS version of the global Park multiplot. Data filtered to contiguous USA
# cells only (lon -125 to -66, lat 24 to 50) via df_conus, excluding Alaska
# and Hawaii. Identical structure, color palette, and factor ordering to the
# global plot so panels are directly comparable. Same 3 x 6 facet layout:
# columns = fire model (CLASSIC, JULES, SSiB4), rows = decade with GMT label.

df_park_long_usa <- df_conus |>
  select(lon, lat, all_of(park_cols)) |>
  pivot_longer(cols = all_of(park_cols), names_to = "col", values_to = "fpm") |>
  mutate(
    model  = case_when(
      str_detect(col, "classic") ~ "CLASSIC",
      str_detect(col, "jules")   ~ "JULES",
      str_detect(col, "ssib4")   ~ "SSiB4"
    ),
    decade    = park_decade_labels[str_extract(col, "\\d{4}s")],
    model     = factor(model,  levels = c("CLASSIC", "JULES", "SSiB4")),
    decade    = factor(decade, levels = unname(park_decade_labels)),
    fill_band = cut_fpm(fpm)
  )

p_park_usa <- ggplot(df_park_long_usa, aes(x = lon, y = lat, fill = fill_band)) +
  geom_tile(width = 0.5, height = 0.5) +
  scale_fill_manual(
    values   = band_colors_fpm,
    name     = "Fire PM2.5 (µg/m³)",
    drop     = FALSE,
    na.value = "lightgray"
  ) +
  coord_fixed(1.3) +
  facet_grid(decade ~ model) +
  labs(
    title = "Park et al. Historical Fire PM2.5 by Model and Decade — CONUS",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid        = element_blank(),
    axis.text         = element_blank(),
    axis.ticks        = element_blank(),
    legend.position   = "bottom",
    legend.title      = element_text(size = 14),
    legend.text       = element_text(size = 13),
    legend.box.margin = margin(t = 0, r = 0, b = -8, l = 0),
    strip.text        = element_text(size = 14),
    plot.title        = element_text(size = 18, hjust = 0.5),
    plot.margin       = margin(t = 4, r = 4, b = 4, l = 4)
  ) +
  guides(fill = guide_legend(nrow = 1))

p_park_usa
ggsave(here("images", "pm", "map_park_multiplot_conus.png"),
       p_park_usa, width = 8.5, height = 11, dpi = 300)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Maps: Pierce + Zhao multiplot -- CONUS ############################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# CONUS version of the global Pierce + Zhao multiplot. Uses df_conus for the
# same CONUS coordinate filter. Reuses pz_meta and pz_row_labels unchanged so
# column/row labels and factor ordering match the global plot exactly.
# Empty Zhao panels in rows 1-3 are retained for layout consistency.

df_pz_long_usa <- df_conus |>
  select(lon, lat, all_of(pierce_zhao_cols)) |>
  pivot_longer(cols = all_of(pierce_zhao_cols), names_to = "col_name", values_to = "fpm") |>
  left_join(pz_meta, by = "col_name") |>
  mutate(
    col_label = factor(col_label, levels = c("Pierce (CESM)", "Zhao (Multi)")),
    row_label = factor(row_label, levels = pz_row_labels),
    fill_band = cut_fpm(fpm)
  )

p_pz_usa <- ggplot(df_pz_long_usa, aes(x = lon, y = lat, fill = fill_band)) +
  geom_tile(width = 0.5, height = 0.5) +
  scale_fill_manual(
    values   = band_colors_fpm,
    name     = "Fire PM2.5 (µg/m³)",
    drop     = FALSE,
    na.value = "lightgray"
  ) +
  coord_fixed(1.3) +
  facet_grid(row_label ~ col_label) +
  labs(
    title = "Pierce et al. (CESM) and Zhao et al. (Multi-model) Fire PM2.5 — CONUS",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid        = element_blank(),
    axis.text         = element_blank(),
    axis.ticks        = element_blank(),
    legend.position   = "bottom",
    legend.title      = element_text(size = 14),
    legend.text       = element_text(size = 13),
    legend.box.margin = margin(t = 0, r = 0, b = -8, l = 0),
    strip.text        = element_text(size = 14),
    plot.title        = element_text(size = 18, hjust = 0.5),
    plot.margin       = margin(t = 4, r = 4, b = 4, l = 4)
  ) +
  guides(fill = guide_legend(nrow = 1))

p_pz_usa
ggsave(here("images", "pm", "map_pierce_zhao_multiplot_conus.png"),
       p_pz_usa, width = 8.5, height = 11, dpi = 300)

### THE END

