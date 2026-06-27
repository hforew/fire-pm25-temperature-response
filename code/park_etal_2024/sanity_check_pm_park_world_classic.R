# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Plot PM2.5 maps - Classic model only, global #################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Goal: Render GLOBAL PM2.5 maps from Park et al. 2024 — Classic model only,
#       one map per decade, with visible 0.5° grid cells over country borders.
#   1. Read regridded PM2.5 CSV; keep only the Classic-model columns
#      (park_classic_<decade>_pm25); list how many were found.
#   2. Load Natural Earth country boundaries for overlay.
#   3. Plot each decade with ggplot/geom_tile (0.5° cells + country borders),
#      share one color scale (0 to 99.9th pctile), save each as base64 PNG, and
#      assemble them into a single self-contained HTML page.
# Input  : output/pm25_park2024.csv
# Output : images/fpm_cell_usa/pm_classic_global_park.html
# Execution order:
#  files run before: pm_regrid_park_2024.R


rm(list = ls())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Packages ######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(here)
library(tidyverse)
library(ggplot2)
library(rnaturalearth)
library(sf)
library(htmltools)
library(base64enc)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ 1. Read data ##################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pm <- read.csv(here("output", "pm25_park2024.csv"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ 2. Subset to Classic model columns only #######################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pm_cols <- grep("^park_classic_.*_pm25$", names(pm), value = TRUE)
cat("Classic PM columns found:", length(pm_cols), "\n")
print(pm_cols)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ 3. World country boundaries ###################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ 4. Color palette (high PM = dark coffee, low PM = white) #####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pm_palette <- c("#3B1F0E", "#6B3E1B", "#A9651A", "#D9A441",
                "#F2D04A", "#FFF2A8", "#FFFFFF")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ 5. Output dir #################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

out_dir <- here("images", "fpm_cell_usa")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ 6. Plot helper ################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# geom_tile draws each 0.5 deg cell with a thin grey edge so the global
# grid is clearly visible. Country borders are overlaid in black.

make_pm_plot <- function(df, col, fill_max) {
  ggplot() +
    geom_tile(data = df,
              aes(x = lon, y = lat, fill = .data[[col]]),
              width = 0.5, height = 0.5,
              color = "grey50", linewidth = 0.03) +
    geom_sf(data = world,
            fill = NA, color = "black", linewidth = 0.2,
            inherit.aes = FALSE) +
    scale_fill_gradientn(
      colors   = rev(pm_palette),   # low -> white, high -> dark coffee
      name     = "PM2.5",
      limits   = c(0, fill_max),
      oob      = scales::squish,
      na.value = "white"
    ) +
    coord_sf(xlim = c(-180, 180), ylim = c(-60, 85), expand = FALSE) +
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
############ 7. Shared color scale (99.9th percentile) ######################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fill_max <- quantile(unlist(pm[pm_cols]), 0.999, na.rm = TRUE)
cat("Shared fill_max (99th percentile):", round(fill_max, 3), "\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ 8. Build one PNG per decade, embed as base64 in HTML #########
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

png_to_base64 <- function(p) {
  paste0("data:image/png;base64,", base64enc::base64encode(p))
}

img_tags <- lapply(pm_cols, function(col) {
  cat("Plotting:", col, "\n")
  p <- make_pm_plot(pm, col, fill_max)
  tmp <- tempfile(fileext = ".png")
  ggsave(tmp, p, width = 12, height = 6.2, dpi = 200, bg = "white")
  tags$img(
    src   = png_to_base64(tmp),
    style = paste("width: 98%; margin: 1%; vertical-align: top;",
                  "border: 1px solid #ddd;",
                  "image-rendering: pixelated;",
                  "image-rendering: crisp-edges;")
  )
})

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ 9. Assemble HTML and save #####################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

page <- tags$html(
  tags$head(
    tags$title("Park et al. 2024 - Classic PM2.5 (Global)"),
    tags$meta(charset = "utf-8")
  ),
  tags$body(
    style = "font-family: sans-serif; max-width: 1500px; margin: 20px auto; padding: 0 20px;",
    tags$h1("PM2.5 - Classic model, Park et al. 2024 (Global)"),
    tags$p(
      style = "color: #555;",
      paste0("Generated ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
             ". Each tile is a 0.5 deg x 0.5 deg cell. Shared color scale: ",
             "0 to 99th percentile across all decades (",
             round(fill_max, 2), ").")
    ),
    tags$hr(),
    do.call(tags$div, img_tags)
  )
)

out_html <- here("images", "fpm_cell_usa", "pm_classic_global_park.html")
save_html(page, out_html)
cat("\nHTML saved to:\n", out_html, "\n")