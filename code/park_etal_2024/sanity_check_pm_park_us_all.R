# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Plot PM2.5 maps for CONUS (0.5 deg cells, clipped to US) #####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read regridded PM2.5 (0.5 deg x 0.5 deg) from Park et al. 2024,
# keep cells overlapping CONUS, clip cell geometry to the US land boundary,
# and render one map per scenario x decade with visible 0.5 deg grid cells.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Packages ######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(here)
library(tidyverse)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(htmltools)
library(base64enc)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ 1. Read data ##################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pm <- read.csv(here("output", "pm25_park2024.csv"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ 2. Pre-filter cells that overlap CONUS bbox ###################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Cell centers are at *.25 / *.75 with cell size 0.5 deg, so a cell at center c
# spans [c - 0.25, c + 0.25]. CONUS bbox is roughly lon [-125, -66.5],
# lat [24, 49.5]. A cell overlaps CONUS iff its center lies within bbox
# expanded by 0.25 deg on each side. This step is just a coarse filter;
# the precise clipping to the US land boundary is done with sf below.

half <- 0.25
lon_min <- -125  - half
lon_max <- -66.5 + half
lat_min <-  24   - half
lat_max <-  49.5 + half

pm_us <- pm %>%
  filter(lon >= lon_min, lon <= lon_max,
         lat >= lat_min, lat <= lat_max)

# Identify all PM columns to plot
pm_cols <- grep("^park_.*_pm25$", names(pm_us), value = TRUE)
cat("PM columns found:", length(pm_cols), "\n")
print(pm_cols)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ 3. US land boundary (CONUS only) ##############################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Pull state-level polygons from Natural Earth and drop Alaska / Hawaii /
# territories so we have a clean CONUS multipolygon. We then dissolve to
# a single boundary used for clipping the 0.5 deg grid cells.

us_states <- rnaturalearth::ne_states(country = "United States of America",
                                      returnclass = "sf")

conus <- us_states %>%
  filter(!name %in% c("Alaska", "Hawaii")) %>%
  st_make_valid()

conus_union <- st_union(conus)   # single CONUS polygon, used for clipping

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ 4. Build 0.5 deg cell polygons and clip to CONUS ##############
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# For each row in pm_us, build a square polygon of side 0.5 deg centered on
# (lon, lat). Then intersect with CONUS so the cells are clipped exactly at
# the coastline / border. Cells that fall entirely outside CONUS disappear.

cell_half <- 0.25

make_cell <- function(lon, lat) {
  st_polygon(list(rbind(
    c(lon - cell_half, lat - cell_half),
    c(lon + cell_half, lat - cell_half),
    c(lon + cell_half, lat + cell_half),
    c(lon - cell_half, lat + cell_half),
    c(lon - cell_half, lat - cell_half)
  )))
}

cells_geom <- mapply(make_cell, pm_us$lon, pm_us$lat, SIMPLIFY = FALSE) %>%
  st_sfc(crs = 4326)

cells_sf <- st_sf(pm_us, geometry = cells_geom)

# Clip to CONUS land boundary
cells_clipped <- st_intersection(cells_sf, conus_union)
cat("Cells before clip:", nrow(cells_sf),
    " | after clip:", nrow(cells_clipped), "\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ 5. Color palette (high PM = dark coffee, low PM = white) #####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pm_palette <- c("#3B1F0E", "#6B3E1B", "#A9651A", "#D9A441",
                "#F2D04A", "#FFF2A8", "#FFFFFF")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ 6. Output dir #################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

out_dir <- here("images", "fpm_cell_usa")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ 7. Plot helper ################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# geom_sf draws each clipped cell as a filled polygon with a thin grey edge,
# so the 0.5 deg grid is clearly visible. State borders are overlaid for
# geographic reference; the outer CONUS boundary is drawn last in black.

make_pm_plot <- function(sf_data, col, fill_max) {
  ggplot() +
    geom_sf(data = sf_data,
            aes(fill = .data[[col]]),
            color = "grey35", linewidth = 0.12) +
    geom_sf(data = conus,
            fill = NA, color = "grey25", linewidth = 0.18) +
    geom_sf(data = conus_union,
            fill = NA, color = "black", linewidth = 0.4) +
    scale_fill_gradientn(
      colors   = rev(pm_palette),   # low -> white, high -> dark coffee
      name     = "PM2.5",
      limits   = c(0, fill_max),
      oob      = scales::squish,
      na.value = "white"
    ) +
    coord_sf(xlim = c(-125, -66.5), ylim = c(24, 49.5), expand = FALSE) +
    labs(title = col) +
    theme_void(base_size = 10) +
    theme(
      plot.title        = element_text(hjust = 0.5, size = 10, face = "bold"),
      legend.key.height = unit(0.6, "cm"),
      legend.key.width  = unit(0.45, "cm"),
      plot.margin       = margin(2, 2, 2, 2)
    )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ 8. Shared color scale (99.9th percentile) ######################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fill_max <- quantile(unlist(st_drop_geometry(cells_clipped)[pm_cols]),
                     0.999, na.rm = TRUE)
cat("Shared fill_max (99th percentile):", round(fill_max, 3), "\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ 9. Build one PNG per column, embed as base64 in HTML ##########
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

png_to_base64 <- function(p) {
  paste0("data:image/png;base64,", base64enc::base64encode(p))
}

img_tags <- lapply(pm_cols, function(col) {
  cat("Plotting:", col, "\n")
  p <- make_pm_plot(cells_clipped, col, fill_max)
  tmp <- tempfile(fileext = ".png")
  ggsave(tmp, p, width = 9, height = 5.5, dpi = 220, bg = "white")
  tags$img(
    src   = png_to_base64(tmp),
    style = paste("width: 49%; margin: 0.5%; vertical-align: top;",
                  "border: 1px solid #ddd;",
                  "image-rendering: pixelated;",
                  "image-rendering: crisp-edges;")
  )
})

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ 10. Assemble HTML and save ####################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

page <- tags$html(
  tags$head(
    tags$title("Park et al. 2024 - PM2.5 maps (CONUS, clipped to US land)"),
    tags$meta(charset = "utf-8")
  ),
  tags$body(
    style = "font-family: sans-serif; max-width: 1400px; margin: 20px auto; padding: 0 20px;",
    tags$h1("PM2.5 - Park et al. 2024 (CONUS, clipped to US boundary)"),
    tags$p(
      style = "color: #555;",
      paste0("Generated ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
             ". Each tile is a 0.5 deg x 0.5 deg cell, clipped to the CONUS ",
             "land boundary. Shared color scale: 0 to 99th percentile across ",
             "all maps (", round(fill_max, 2), ").")
    ),
    tags$hr(),
    do.call(tags$div, img_tags)
  )
)

out_html <- here("images", "fpm_cell_usa", "pm_cell_usa_park.html")
save_html(page, out_html)
cat("\nHTML saved to:\n", out_html, "\n")