# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ fPM Global Cell Maps by Variable Group ######################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())

library(here)
library(data.table)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)
library(patchwork)
library(htmltools)
library(base64enc)

sf_use_s2(TRUE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Read Combined pop_pm Dataset ################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pop_pm_final <- fread(here("output", "pop_pm_combined_final.csv"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Subset fpm Variables with Coordinates #######################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fpm_cols <- grep("fpm", names(pop_pm_final), value = TRUE, ignore.case = TRUE)
fpm_data <- pop_pm_final[, c("lon", "lat", fpm_cols), with = FALSE]

cat("Number of fpm variables:", length(fpm_cols), "\n")
cat("Dimensions of fpm_data:", dim(fpm_data), "\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Group fpm Columns: zhao / park / other ######################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
zhao_cols  <- grep("zhao", fpm_cols, value = TRUE, ignore.case = TRUE)
park_cols  <- grep("park", fpm_cols, value = TRUE, ignore.case = TRUE)
other_cols <- setdiff(fpm_cols, c(zhao_cols, park_cols))

# Sort each group by earliest 4-digit year embedded in the column name
get_time_key <- function(x) {
  m <- regmatches(x, regexpr("[0-9]{4}", x))
  if (length(m) == 0) return(NA_integer_)
  as.integer(m)
}

zhao_cols  <- zhao_cols[order(sapply(zhao_cols,  get_time_key))]
park_cols  <- park_cols[order(sapply(park_cols,  get_time_key))]
other_cols <- other_cols[order(sapply(other_cols, get_time_key))]

cat("\nzhao_cols  (n=", length(zhao_cols),  "):\n", sep = ""); print(zhao_cols)
cat("\npark_cols  (n=", length(park_cols),  "):\n", sep = ""); print(park_cols)
cat("\nother_cols (n=", length(other_cols), "):\n", sep = ""); print(other_cols)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Build Global 0.5-degree Cell Polygons #######################################
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

global_sf <- st_sf(
  fpm_data[, .(lon, lat)],
  fpm_data[, fpm_cols, with = FALSE],
  geometry = st_sfc(
    Map(make_cell_polygon, fpm_data$lon, fpm_data$lat),
    crs = 4326
  )
)

cat("\nGlobal sf object built. Total cells:", nrow(global_sf), "\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Country Borders for Map Overlay #############################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
world_sf <- ne_countries(scale = "medium", returnclass = "sf")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Color Palette: Dark Coffee -> Yellow -> White ###############################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# High fPM = dark coffee, low fPM = white. ggplot maps colors[1] to min, last to max.
fpm_palette <- c("#3B1F0E", "#6B3E1B", "#A9651A", "#D9A441", "#F2D04A", "#FFF2A8", "#FFFFFF")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Plotting Function: One Global Map per fpm Variable ##########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_global_fpm <- function(sf_obj, var) {
  ggplot(sf_obj) +
    geom_sf(aes(fill = .data[[var]]), color = NA) +
    geom_sf(data = world_sf, fill = NA, color = "grey30", linewidth = 0.15) +
    scale_fill_gradientn(
      colors   = rev(fpm_palette),
      name     = "fPM",
      na.value = "white"
    ) +
    coord_sf(xlim = c(-180, 180), ylim = c(-60, 85), expand = FALSE) +
    labs(title = var) +
    theme_void(base_size = 10) +
    theme(
      plot.title        = element_text(hjust = 0.5, size = 10, face = "bold"),
      legend.key.height = unit(0.45, "cm"),
      legend.key.width  = unit(0.3, "cm"),
      legend.title      = element_text(size = 7),
      legend.text       = element_text(size = 6),
      plot.margin       = margin(3, 3, 3, 3)
    )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Build One Row of Plots per Group ############################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
zhao_plots  <- lapply(zhao_cols,  function(v) plot_global_fpm(global_sf, v))
park_plots  <- lapply(park_cols,  function(v) plot_global_fpm(global_sf, v))
other_plots <- lapply(other_cols, function(v) plot_global_fpm(global_sf, v))

zhao_row  <- wrap_plots(zhao_plots,  nrow = 1)
park_row  <- wrap_plots(park_plots,  nrow = 1)
other_row <- wrap_plots(other_plots, nrow = 1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Save Each Group Row to Standalone HTML ######################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
out_dir <- "C:/Ford_BA_FPM25/images/fpm_cell_global"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

save_row_html <- function(plot_obj, n_cols, file_name, title) {
  png_path <- tempfile(fileext = ".png")
  per_panel_w <- 5            # global maps are wider (2:1 aspect)
  total_w     <- max(6, per_panel_w * n_cols)
  total_h     <- 3.2
  ggsave(
    filename  = png_path,
    plot      = plot_obj,
    width     = total_w,
    height    = total_h,
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

save_row_html(zhao_row,  length(zhao_cols),  "fpm_global_zhao.html",  "fPM Global - Zhao variables")
save_row_html(park_row,  length(park_cols),  "fpm_global_park.html",  "fPM Global - Park variables")
save_row_html(other_row, length(other_cols), "fpm_global_other.html", "fPM Global - Other fpm variables")

cat("\nHTML files saved to:", out_dir, "\n")

#THE END