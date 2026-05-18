# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ fPM data plot at cell level USA #########################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Purpose:
#   Visualize fire-related PM2.5 (fPM) estimates at the grid-cell level across the
#   United States for multiple data sources and time periods,
#   enabling side-by-side spatial comparison of fPM products.
#
# What this script does:
#   1. Reads the combined population and PM dataset (pop_pm_combined_final.csv) and
#      subsets all fPM variables together with their cell center coordinates
#      (lon, lat).
#   2. Groups fPM variables by source — Zhao, Park, and Other — and sorts each
#      group chronologically by the year embedded in the variable name.
#   3. Builds 0.5-degree polygon cells from each (lon, lat) center and keeps only
#      cells overlapping the US mainland, excluding Alaska, Hawaii, and
#      US territories via a bounding-box crop.
#   4. Plots each fPM variable as a choropleth map with state borders overlaid,
#      using a shared dark-coffee -> yellow -> white color palette. Each panel
#      uses its own color scale so within-panel spatial patterns are emphasized.
#   5. Arranges plots into one row per source group (oldest -> newest, left ->
#      right) and exports each row as a standalone.
#
# Output:
#   images/fpm_cell_usa/
#     - fpm_us_zhao.html
#     - fpm_us_park.html
#     - fpm_us_other.html
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(here)
library(data.table)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)
library(patchwork)
library(htmltools)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ read in fpm data from final dataset #########################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pop_pm_final <- fread(here("output", "pop_pm_combined_final.csv"))

head(pop_pm_final)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Subset fPM Variables with Coordinates #######################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fpm_cols <- grep("fpm", names(pop_pm_final), value = TRUE, ignore.case = TRUE)

fpm_data <- pop_pm_final[, c("lon", "lat", fpm_cols), with = FALSE]

head(fpm_data)
cat("\nNumber of fpm variables:", length(fpm_cols), "\n")
cat("Dimensions of fpm_data:", dim(fpm_data), "\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ fPM Maps over US Mainland by Variable Group #################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sf_use_s2(TRUE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Identify fpm Columns and Group Them #########################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fpm_cols  <- grep("fpm", names(fpm_data), value = TRUE, ignore.case = TRUE)
fpm_cols  <- setdiff(fpm_cols, c("lon", "lat"))

zhao_cols  <- grep("zhao", fpm_cols, value = TRUE, ignore.case = TRUE)
park_cols  <- grep("park", fpm_cols, value = TRUE, ignore.case = TRUE)
other_cols <- setdiff(fpm_cols, c(zhao_cols, park_cols))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Helper: Extract Earliest Year/Decade for Time Sorting #######################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
get_time_key <- function(x) {
  m <- regmatches(x, regexpr("[0-9]{4}", x))
  if (length(m) == 0) return(NA_integer_)
  as.integer(m)
}

zhao_cols  <- zhao_cols[order(sapply(zhao_cols, get_time_key))]
park_cols  <- park_cols[order(sapply(park_cols, get_time_key))]
other_cols <- other_cols[order(sapply(other_cols, get_time_key))]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Build 0.5-degree Cell Polygons ##############################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
make_cell_polygon <- function(lon, lat, half_res = 0.25) {
  st_polygon(list(matrix(
    c(lon - half_res, lat - half_res,
      lon + half_res, lat - half_res,
      lon + half_res, lat + half_res,
      lon - half_res, lat + half_res,
      lon - half_res, lat - half_res),
    ncol = 2, byrow = TRUE
  )))
}

cell_sf <- st_sf(
  fpm_data[, .(cell_id = .I, lon, lat)],
  geometry = st_sfc(
    Map(make_cell_polygon, fpm_data$lon, fpm_data$lat),
    crs = 4326
  )
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Filter Cells Overlapping the US Mainland (CONUS) ############################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
us_sf <- ne_countries(scale = "medium", returnclass = "sf")
us_sf <- us_sf[us_sf$admin == "United States of America", ]

# Crop to mainland (CONUS) bounding box: exclude Alaska, Hawaii, territories
conus_bbox <- st_as_sfc(st_bbox(c(
  xmin = -125, xmax = -66.5, ymin = 24, ymax = 49.5
), crs = 4326))
us_mainland <- st_intersection(us_sf, conus_bbox)

us_idx   <- lengths(st_intersects(cell_sf, us_mainland)) > 0
us_cells <- cell_sf[us_idx, ]

us_fpm <- cbind(
  us_cells[, c("lon", "lat")],
  fpm_data[us_idx, fpm_cols, with = FALSE]
)

cat("\nNumber of cells overlapping US mainland:", sum(us_idx), "\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ State Borders for Map Overlay ###############################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
states_sf <- ne_states(country = "United States of America", returnclass = "sf")
states_sf <- st_intersection(states_sf, conus_bbox)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Color Palette: Dark Coffee -> Yellow -> White ###############################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# High fPM = dark coffee, Low fPM = white. Order palette accordingly.
fpm_palette <- c("#3B1F0E", "#6B3E1B", "#A9651A", "#D9A441", "#F2D04A", "#FFF2A8", "#FFFFFF")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Plotting Function for One fpm Variable ######################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_fpm_var <- function(sf_obj, var) {
  ggplot(sf_obj) +
    geom_sf(aes(fill = .data[[var]]), color = NA) +
    geom_sf(data = states_sf, fill = NA, color = "grey30", linewidth = 0.2) +
    scale_fill_gradientn(
      colors = rev(fpm_palette),
      name = "fPM",
      na.value = "white"
    ) +
    coord_sf(xlim = c(-125, -66.5), ylim = c(24, 49.5), expand = FALSE) +
    labs(title = var) +
    theme_void(base_size = 9) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 9, face = "bold"),
      legend.key.height = unit(0.35, "cm"),
      legend.key.width  = unit(0.25, "cm"),
      legend.title = element_text(size = 7),
      legend.text  = element_text(size = 6),
      plot.margin = margin(2, 2, 2, 2)
    )
}

# ---- Compute shared fill_max per group (matches PM script's 99.9-pct rule) ----
shared_max <- function(cols) {
  vals <- unlist(st_drop_geometry(us_fpm)[, cols, drop = FALSE])
  quantile(vals, 0.999, na.rm = TRUE)
}

zhao_max  <- shared_max(zhao_cols)
park_max  <- shared_max(park_cols)
other_max <- if (length(other_cols)) shared_max(other_cols) else NA

cat("Shared fill_max — zhao:", round(zhao_max, 3),
    "| park:",  round(park_max, 3),
    "| other:", round(other_max, 3), "\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Build Plots: One Row per Group, Time Old -> New Left -> Right ###############
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
zhao_plots  <- lapply(zhao_cols,  function(v) plot_fpm_var(us_fpm, v))
park_plots  <- lapply(park_cols,  function(v) plot_fpm_var(us_fpm, v))
other_plots <- lapply(other_cols, function(v) plot_fpm_var(us_fpm, v))

zhao_row  <- wrap_plots(zhao_plots,  nrow = 1)
park_row  <- wrap_plots(park_plots,  nrow = 1)
other_row <- wrap_plots(other_plots, nrow = 1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Save Each Row to a Standalone HTML File #####################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
out_dir <- here::here("images", "fpm_cell_usa")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!requireNamespace("base64enc", quietly = TRUE)) install.packages("base64enc")
library(base64enc)

save_row_html <- function(plot_obj, n_cols, file_name, title) {
  png_path <- tempfile(fileext = ".png")
  
  # Cap per-panel width so total stays reasonable, but allow oversized output
  per_panel_w <- 2.4
  total_w     <- max(4, per_panel_w * n_cols)
  
  ggsave(
    filename  = png_path,
    plot      = plot_obj,
    width     = total_w,
    height    = 3.0,
    dpi       = 150,
    bg        = "white",
    limitsize = FALSE
  )
  
  png_b64 <- base64enc::base64encode(png_path)
  html <- tags$html(
    tags$head(tags$title(title)),
    tags$body(
      tags$h3(title),
      tags$div(
        style = "overflow-x: auto; white-space: nowrap;",
        tags$img(
          src   = paste0("data:image/png;base64,", png_b64),
          style = "height: auto; max-width: none;"
        )
      )
    )
  )
  save_html(html, file = file.path(out_dir, file_name))
}

save_row_html(zhao_row,  length(zhao_cols),  "fpm_us_zhao.html",  "fPM US Mainland - Zhao variables")
save_row_html(park_row,  length(park_cols),  "fpm_us_park.html",  "fPM US Mainland - Park variables")
save_row_html(other_row, length(other_cols), "fpm_us_other.html", "fPM US Mainland - Other fpm variables")

cat("\nHTML files saved to:", out_dir, "\n")

#THE END