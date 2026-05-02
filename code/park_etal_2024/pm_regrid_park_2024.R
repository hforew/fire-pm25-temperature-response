# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
########### Extract Year, Longitude, Latitude, and PM2.5 from NetCDF Files ##################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Remove all objects from the environment
rm(list = ls())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Packages #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(here)
library(tidyverse)
library(ncdf4)
library(R.matlab)
library(fields)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ File Paths for All Scenarios ##################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

years <- c(1965, 1975, 1985, 1995, 2005, 2015)

# Scenario 1: withoutfire
file_paths_withoutfire <- c(
  here("input", "Park_etal_2024", "GEOSChem_output", "withoutfire", "05x05_CEDS_1965_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "withoutfire", "05x05_CEDS_1975_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "withoutfire", "05x05_CEDS_1985_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "withoutfire", "05x05_CEDS_1995_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "withoutfire", "05x05_CEDS_2005_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "withoutfire", "05x05_CEDS_2015_on_off_pm25_Surface_Re_yearavg.nc4")
)

# Scenario 2: classic
file_paths_classic <- c(
  here("input", "Park_etal_2024", "GEOSChem_output", "classic", "obsclim", "05x05_CEDS_1965_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "classic", "obsclim", "05x05_CEDS_1975_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "classic", "obsclim", "05x05_CEDS_1985_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "classic", "obsclim", "05x05_CEDS_1995_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "classic", "obsclim", "05x05_CEDS_2005_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "classic", "obsclim", "05x05_CEDS_2015_on_off_pm25_Surface_Re_yearavg.nc4")
)

# Scenario 3: jules
file_paths_jules <- c(
  here("input", "Park_etal_2024", "GEOSChem_output", "jules", "obsclim", "05x05_CEDS_1965_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "jules", "obsclim", "05x05_CEDS_1975_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "jules", "obsclim", "05x05_CEDS_1985_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "jules", "obsclim", "05x05_CEDS_1995_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "jules", "obsclim", "05x05_CEDS_2005_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "jules", "obsclim", "05x05_CEDS_2015_on_off_pm25_Surface_Re_yearavg.nc4")
)

# Scenario 4: ssib4
file_paths_ssib4 <- c(
  here("input", "Park_etal_2024", "GEOSChem_output", "ssib4", "obsclim", "05x05_CEDS_1965_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "ssib4", "obsclim", "05x05_CEDS_1975_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "ssib4", "obsclim", "05x05_CEDS_1985_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "ssib4", "obsclim", "05x05_CEDS_1995_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "ssib4", "obsclim", "05x05_CEDS_2005_on_off_pm25_Surface_Re_yearavg.nc4"),
  here("input", "Park_etal_2024", "GEOSChem_output", "ssib4", "obsclim", "05x05_CEDS_2015_on_off_pm25_Surface_Re_yearavg.nc4")
)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Target Grid (0.5 deg x 0.5 deg centers) #######################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

target_lon <- seq(-179.75, 179.75, by = 0.5)  # 720
target_lat <- seq(-89.75,  89.75,  by = 0.5)  # 360

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Helper: edge -> center #########################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

convert_edge_to_center <- function(coords) {
  resolution <- mean(diff(coords))
  remainder  <- abs(coords[1] %% resolution)
  if (remainder < 1e-6 || abs(remainder - resolution) < 1e-6) {
    return(coords + resolution / 2)
  } else {
    return(coords)
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Helper: 0-360 lon -> -180-180 ##################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

shift_lon_to_180 <- function(lon, mat) {
  if (max(lon) > 180) {
    lon_new <- ifelse(lon > 180, lon - 360, lon)
    ord     <- order(lon_new)
    lon_new <- lon_new[ord]
    mat     <- mat[ord, , drop = FALSE]
    return(list(lon = lon_new, mat = mat))
  }
  list(lon = lon, mat = mat)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Function to Extract + Regrid ##################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

extract_pm25_data <- function(file_paths, years, scenario_name,
                              target_lon, target_lat) {
  
  pm25_list <- list()
  
  for (i in seq_along(file_paths)) {
    
    nc <- nc_open(file_paths[i])
    
    # 1. Read actual lon/lat from the file
    lon_name <- if ("lon" %in% names(nc$dim)) "lon" else "longitude"
    lat_name <- if ("lat" %in% names(nc$dim)) "lat" else "latitude"
    lon_src  <- as.numeric(ncvar_get(nc, lon_name))
    lat_src  <- as.numeric(ncvar_get(nc, lat_name))
    
    cat("lon range:", range(lon_src), " | n =", length(lon_src), "\n")
    cat("lat range:", range(lat_src), " | n =", length(lat_src), "\n")
    cat("first 5 lon:", head(lon_src, 5), "\n")
    cat("first 5 lat:", head(lat_src, 5), "\n")
    
    # Edge -> center if needed
    lon_src <- convert_edge_to_center(lon_src)
    lat_src <- convert_edge_to_center(lat_src)
    
    # 2. Read PM2.5 and reduce to 2D
    pm25_raw <- ncvar_get(nc, "PM25")
    dims     <- dim(pm25_raw)
    pm25_surface <- switch(as.character(length(dims)),
                           "4" = pm25_raw[, , 1, 1],
                           "3" = pm25_raw[, , 1],
                           "2" = pm25_raw,
                           stop("Unexpected PM25 dimensions: ", paste(dims, collapse = " x "))
    )
    nc_close(nc)
    
    # 3. Ensure orientation is [lon, lat]
    if (nrow(pm25_surface) == length(lat_src) &&
        ncol(pm25_surface) == length(lon_src)) {
      pm25_surface <- t(pm25_surface)
    }
    stopifnot(nrow(pm25_surface) == length(lon_src),
              ncol(pm25_surface) == length(lat_src))
    
    # 4. Convert lon to -180..180 if needed
    sh           <- shift_lon_to_180(lon_src, pm25_surface)
    lon_src      <- sh$lon
    pm25_surface <- sh$mat
    
    # 5. Make latitude ascending (required by interp.surface.grid)
    if (is.unsorted(lat_src)) {
      ord          <- order(lat_src)
      lat_src      <- lat_src[ord]
      pm25_surface <- pm25_surface[, ord, drop = FALSE]
    }
    
    # 6. Regrid to target 0.5 deg grid (skip if already identical)
    same_grid <-
      length(lon_src) == length(target_lon) &&
      length(lat_src) == length(target_lat) &&
      max(abs(lon_src - target_lon)) < 1e-6 &&
      max(abs(lat_src - target_lat)) < 1e-6
    
    if (same_grid) {
      pm25_regrid <- pm25_surface
    } else {
      obj <- list(x = lon_src, y = lat_src, z = pm25_surface)
      pm25_regrid <- interp.surface.grid(
        obj,
        grid.list = list(x = target_lon, y = target_lat)
      )$z
    }
    
    # 7. Long format (lon varies first)
    df <- expand.grid(
      lon = target_lon,
      lat = target_lat,
      KEEP.OUT.ATTRS  = FALSE,
      stringsAsFactors = FALSE
    )
    df$pm25 <- as.vector(pm25_regrid)
    df$year <- years[i]
    df <- df %>% select(year, lon, lat, pm25)
    pm25_list[[i]] <- df
    
    cat(sprintf("[%s] %d: src %d x %d -> target %d x %d %s\n",
                scenario_name, years[i],
                length(lon_src), length(lat_src),
                length(target_lon), length(target_lat),
                if (same_grid) "(no regrid)" else "(bilinear regridded)"))
  }
  
  bind_rows(pm25_list) %>%
    rename(!!paste0("pm25_", scenario_name) := pm25)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Extract Data for All Scenarios #################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

withoutfire <- extract_pm25_data(file_paths_withoutfire, years, "withoutfire",
                                 target_lon, target_lat)
classic     <- extract_pm25_data(file_paths_classic,     years, "classic",
                                 target_lon, target_lat)
jules       <- extract_pm25_data(file_paths_jules,       years, "jules",
                                 target_lon, target_lat)
ssib4       <- extract_pm25_data(file_paths_ssib4,       years, "ssib4",
                                 target_lon, target_lat)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Merge All Scenarios by Year, Lon, Lat #########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Merge all four datasets by year, lon, lat
pm25_all <- withoutfire %>%
  left_join(classic, by = c("year", "lon", "lat")) %>%
  left_join(jules, by = c("year", "lon", "lat")) %>%
  left_join(ssib4, by = c("year", "lon", "lat"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Calculate Fire PM (Difference from Without Fire) ##############
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Calculate fire PM for each model by subtracting withoutfire PM
# Fire PM = Model PM - Without Fire PM
pm25_all <- pm25_all %>%
  mutate(
    classic_fpm = pm25_classic - pm25_withoutfire,
    jules_fpm = pm25_jules - pm25_withoutfire,
    ssib4_fpm = pm25_ssib4 - pm25_withoutfire
  )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Create Final Two Datasets ######################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Dataset 1: Complete dataset with all PM columns and fire PM columns
# Columns: year, lon, lat, pm25_withoutfire, pm25_classic, pm25_jules, pm25_ssib4, 
#          classic_fpm, jules_fpm, ssib4_fpm
pm25_all_complete <- pm25_all

# Dataset 2: Dataset with only fire PM columns (removed all pm25_* columns)
# Columns: year, lon, lat, classic_fpm, jules_fpm, ssib4_fpm
pm25_fpm_only <- pm25_all %>%
  select(year, lon, lat, classic_fpm, jules_fpm, ssib4_fpm)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Convert Fire PM to Wide Format by Year ########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Convert pm25_fpm_only from long to wide format
# Each model x year combination becomes a separate column
# Example: classic_1965_fpm, classic_1975_fpm, ..., jules_1965_fpm, etc.

pm25_fpm_wide <- pm25_fpm_only %>%
  pivot_longer(
    cols = c(classic_fpm, jules_fpm, ssib4_fpm),
    names_to = "model",
    values_to = "fpm"
  ) %>%
  mutate(
    # Create column name: model_year_fpm (e.g., classic_1965_fpm)
    model_year = paste0(gsub("_fpm", "", model), "_", year, "_fpm")
  ) %>%
  select(-model, -year) %>%
  pivot_wider(
    names_from = model_year,
    values_from = fpm
  ) %>%
  # Add "park_" prefix to all columns except lon and lat
  rename_with(
    ~ paste0("park_", .x),
    .cols = -c(lon, lat)
  )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ merge with pop_regrid_park_2024 and our master data ############
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pop_pm_combined <- read.csv(here("output", "pop_pm_combined.csv"))
pop_park <- read.csv(here("output", "pop_regrid_park_2024.csv"))

# Merge pop_wide with pm25_fpm_wide by lon and lat
pm25_pop_merged <- pm25_fpm_wide %>%
  left_join(pop_park, by = c("lon", "lat"))

# Merge pm25_pop_merged with pop_pm_combined by lon and lat
combined_data <- pop_pm_combined %>%
  left_join(pm25_pop_merged, by = c("lon", "lat"))

# Check the results
dim(combined_data)  # Check dimensions
head(combined_data)  # View first few rows
summary(combined_data)  # Summary statistics

# Check for missing values after merge (indicates non-matching rows)
sum(is.na(combined_data))  # Total NA count

# Check if coordinates match perfectly
anti_join(pop_pm_combined, pm25_pop_merged, by = c("lon", "lat"))  # Rows in pop_pm_combined not in pm25_pop_merged
anti_join(pm25_pop_merged, pop_pm_combined, by = c("lon", "lat"))  # Rows in pm25_pop_merged not in pop_pm_combined

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Validation: Convert Back to Matrix ############################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Select park_classic_2015_fpm for validation
fpm_vector <- combined_data$park_classic_2015_fpm

# Convert back to matrix format (720 lon x 360 lat)
fpm_matrix <- matrix(fpm_vector, nrow = 720, ncol = 360)

# Check results
cat("\n=== Validation: Converted to Matrix ===\n")
cat("Matrix dimensions:", dim(fpm_matrix), "\n")
cat("\nFirst 5x5 corner:\n")
print(fpm_matrix[1:5, 1:5])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Validation: Compare with Original MAT File ####################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Read the MAT file
mat_data <- readMat(here("input", "Park_etal_2024", "results1206.mat"))

# Extract Classic model data
PMfire_classic <- mat_data$PMfire[[1]][[1]]  # Classic model: 360(lat) x 720(lon) x 6(years) x 2(scenarios)

cat("\n=== Original MAT File ===\n")
cat("PMfire_classic dimensions:", dim(PMfire_classic), "\n")

# Extract 2015 data (year 6, scenario 1)
PMfire_classic_2015 <- PMfire_classic[, , 6, 1]  # 360 x 720

# Convert our data to matrix and transpose to match
fpm_vector <- combined_data$park_classic_2015_fpm
fpm_matrix <- matrix(fpm_vector, nrow = 720, ncol = 360)
fpm_matrix_t <- t(fpm_matrix)  # Now 360 x 720

# Round to some decimal places for comparison
PMfire_rounded <- round(PMfire_classic_2015, 6)
fpm_rounded <- round(fpm_matrix_t, 6)

# Compare
max_diff <- max(abs(fpm_rounded - PMfire_rounded))
mean_diff <- mean(abs(fpm_rounded - PMfire_rounded))

cat("Our matrix dimensions:", dim(fpm_matrix_t), "\n")
cat("MAT matrix dimensions:", dim(PMfire_classic_2015), "\n")
cat("Maximum difference:", max_diff, "\n")
cat("Mean difference:", mean_diff, "\n")

cat("\nFirst 5x5 from MAT file (rounded):\n")
print(PMfire_rounded[1:5, 1:5])

cat("\nFirst 5x5 from our data (rounded):\n")
print(fpm_rounded[1:5, 1:5])

if (max_diff == 0) {
  cat("\n SUCCESS: Matrices are identical at 4 decimal precision!\n")
} else {
  cat("\n WARNING: Matrices differ at 4 decimal precision!\n")
}


# Truncate to some decimal places (keep only some digits after decimal)
PMfire_truncated <- trunc(PMfire_classic_2015 * 1e6) / 1e6
fpm_truncated <- trunc(fpm_matrix_t * 1e6) / 1e6

# Compare
max_diff <- max(abs(fpm_truncated - PMfire_truncated))
mean_diff <- mean(abs(fpm_truncated - PMfire_truncated))

cat("Our matrix dimensions:", dim(fpm_matrix_t), "\n")
cat("MAT matrix dimensions:", dim(PMfire_classic_2015), "\n")
cat("Maximum difference:", max_diff, "\n")
cat("Mean difference:", mean_diff, "\n")

cat("\nFirst 5x5 from MAT file (truncated):\n")
print(PMfire_truncated[1:5, 1:5])

cat("\nFirst 5x5 from our data (truncated):\n")
print(fpm_truncated[1:5, 1:5])

if (max_diff == 0) {
  cat("\n SUCCESS: Matrices are identical after truncation!\n")
} else {
  cat("\n Difference remains:", max_diff, "\n")
}

# Find different cells
diff_matrix <- fpm_matrix_t - PMfire_classic_2015
diff_indices <- which(abs(diff_matrix) > 0, arr.ind = TRUE)

cat("\nNumber of different cells:", nrow(diff_indices), "\n")
if (nrow(diff_indices) > 0) {
  cat("\nFirst 10 different cells:\n")
  for (i in 1:min(20, nrow(diff_indices))) {
    row <- diff_indices[i, 1]
    col <- diff_indices[i, 2]
    cat(sprintf("Position[%d,%d]: MAT=%.10f, Ours=%.10f, Diff=%.10e\n", 
                row, col, PMfire_classic_2015[row, col], fpm_matrix_t[row, col], diff_matrix[row, col]))
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
######################## save final master data to local ####################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# rename columns for final save

combined_data <- combined_data %>%
  rename_with(
    .cols = starts_with("park_"),  # Target only "park_" columns
    .fn = function(x) {
      for (y in years) {
        # Replace each specific year with its corresponding decade
        x <- str_replace_all(
          x,
          as.character(y),
          paste0(floor(y / 10) * 10, "s")
        )
      }
      x
    }
  )

write.csv(combined_data, here("output", "pop_pm_combined_with_park2024.csv"), row.names = FALSE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ SANITY CHECK: Plot US fpm for ALL models x years ##############
############ Output as HTML to local folder                       ##########
############ Standalone — does not affect final dataset           ##########
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(ncdf4)
library(ggplot2)
library(dplyr)
library(maps)
library(patchwork)
library(htmltools)
library(base64enc)

# --- Output settings ---
out_dir  <- "C:/Ford_BA_FPM25/images/fpm_cell_usa"
out_name <- "fpm_cell_usa_park_sanitycheck"
out_html <- file.path(out_dir, paste0(out_name, ".html"))

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# --- helper: read PM25 surface field with original lon/lat ---
read_pm25 <- function(path) {
  nc <- nc_open(path)
  lon_name <- if ("lon" %in% names(nc$dim)) "lon" else "longitude"
  lat_name <- if ("lat" %in% names(nc$dim)) "lat" else "latitude"
  lon <- as.numeric(ncvar_get(nc, lon_name))
  lat <- as.numeric(ncvar_get(nc, lat_name))
  
  # ====== Print ORIGINAL lon/lat from the NetCDF file ======
  cat("\n----- ORIGINAL coordinates from file -----\n")
  cat("File: ", basename(path), "\n", sep = "")
  cat("Variable names:  lon = '", lon_name, "' | lat = '", lat_name, "'\n", sep = "")
  cat(sprintf("lon: n = %d | range = [%.4f, %.4f] | first 5 = %s | last 5 = %s\n",
              length(lon),
              min(lon), max(lon),
              paste(round(head(lon, 5), 4), collapse = ", "),
              paste(round(tail(lon, 5), 4), collapse = ", ")))
  cat(sprintf("lat: n = %d | range = [%.4f, %.4f] | first 5 = %s | last 5 = %s\n",
              length(lat),
              min(lat), max(lat),
              paste(round(head(lat, 5), 4), collapse = ", "),
              paste(round(tail(lat, 5), 4), collapse = ", ")))
  cat(sprintf("lon resolution = %.4f | lat resolution = %.4f\n",
              mean(diff(lon)), mean(diff(lat))))
  cat("------------------------------------------\n")
  # =========================================================
  
  pm <- ncvar_get(nc, "PM25")
  d  <- dim(pm)
  pm <- switch(as.character(length(d)),
               "4" = pm[, , 1, 1],
               "3" = pm[, , 1],
               "2" = pm,
               stop("Unexpected PM25 dimensions"))
  nc_close(nc)
  
  if (nrow(pm) == length(lat) && ncol(pm) == length(lon)) {
    pm <- t(pm)
  }
  list(lon = lon, lat = lat, pm = pm)
}

# --- helper: build fpm dataframe for one model x one year ---
build_fpm_df <- function(path_model, path_withoutfire) {
  a <- read_pm25(path_model)
  b <- read_pm25(path_withoutfire)
  stopifnot(length(a$lon) == length(b$lon),
            length(a$lat) == length(b$lat))
  
  fpm_mat <- a$pm - b$pm
  
  df <- expand.grid(lon = a$lon, lat = a$lat,
                    KEEP.OUT.ATTRS = FALSE) %>%
    mutate(fpm = as.vector(fpm_mat))
  
  if (max(df$lon) > 180) {
    df <- df %>% mutate(lon = ifelse(lon > 180, lon - 360, lon))
  }
  
  df %>% filter(lon >= -125, lon <= -66, lat >= 24, lat <= 50)
}

# --- helper: single plot ---
us_map <- map_data("state")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Color Palette: White -> Yellow -> Dark Coffee #################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Low fPM = white, High fPM = dark coffee
fpm_palette <- c("#3B1F0E", "#6B3E1B", "#A9651A", "#D9A441",
                 "#F2D04A", "#FFF2A8", "#FFFFFF")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ make_plot: matching plot_fpm_var styling ######################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

make_plot <- function(df, title, fill_max) {
  ggplot(df, aes(x = lon, y = lat, fill = fpm)) +
    geom_raster() +
    geom_polygon(data = us_map,
                 aes(x = long, y = lat, group = group),
                 fill = NA, color = "grey30", linewidth = 0.2,
                 inherit.aes = FALSE) +
    scale_fill_gradientn(
      colors   = rev(fpm_palette),     # white -> ... -> dark coffee
      name     = "fPM",
      limits   = c(0, fill_max),
      oob      = scales::squish,
      na.value = "white"
    ) +
    coord_quickmap(xlim = c(-125, -66.5), ylim = c(24, 49.5),
                   expand = FALSE) +
    labs(title = title) +
    theme_void(base_size = 9) +
    theme(
      plot.title        = element_text(hjust = 0.5, size = 9, face = "bold"),
      legend.key.height = unit(0.35, "cm"),
      legend.key.width  = unit(0.25, "cm"),
      legend.title      = element_text(size = 7),
      legend.text       = element_text(size = 6),
      plot.margin       = margin(2, 2, 2, 2)
    )
}

# --- build one combined image per model, save to temp PNG, then embed ---
models <- list(
  classic = file_paths_classic,
  jules   = file_paths_jules,
  ssib4   = file_paths_ssib4
)

png_to_base64 <- function(png_path) {
  paste0("data:image/png;base64,",
         base64enc::base64encode(png_path))
}

img_tags <- list()

for (m in names(models)) {
  
  cat("Building plot for model:", m, "\n")
  
  dfs <- lapply(seq_along(years), function(i) {
    build_fpm_df(models[[m]][i], file_paths_withoutfire[i])
  })
  
  fill_max <- quantile(unlist(lapply(dfs, `[[`, "fpm")),
                       0.99, na.rm = TRUE)
  
  plots <- lapply(seq_along(years), function(i) {
    make_plot(dfs[[i]],
              title = paste0(m, " — ", years[i]),
              fill_max = fill_max)
  })
  
  combined <- wrap_plots(plots, nrow = 2, ncol = 3) +
    plot_annotation(
      title    = paste0("Sanity check: ", m, " fire PM2.5 (CONUS)"),
      subtitle = "Original NetCDF coordinates, common color scale (0 to 99th percentile)",
      theme    = theme(plot.title    = element_text(face = "bold", size = 13),
                       plot.subtitle = element_text(size = 10))
    ) +
    plot_layout(guides = "collect") & theme(legend.position = "right")
  
  # Save to a temp PNG first, then read back as base64 for embedding
  tmp_png <- tempfile(fileext = ".png")
  ggsave(tmp_png, combined, width = 14, height = 8, dpi = 130, bg = "white")
  
  img_tags[[m]] <- tags$div(
    style = "margin-bottom: 40px;",
    tags$h2(paste0("Model: ", m),
            style = "font-family: sans-serif; color: #333;"),
    tags$img(src   = png_to_base64(tmp_png),
             style = "max-width: 100%; height: auto; border: 1px solid #ddd;")
  )
}

# --- assemble HTML ---
page <- tags$html(
  tags$head(
    tags$title("Park et al. 2024 — Fire PM2.5 Sanity Check (CONUS)"),
    tags$meta(charset = "utf-8")
  ),
  tags$body(
    style = "font-family: sans-serif; max-width: 1400px; margin: 20px auto; padding: 0 20px;",
    tags$h1("Fire PM2.5 Sanity Check — Park et al. 2024"),
    tags$p(
      style = "color: #555;",
      paste0("Generated on ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
             ". Plots use original NetCDF coordinates (no regrid). ",
             "Color scale per model: 0 to 99th percentile of all years.")
    ),
    tags$hr(),
    img_tags
  )
)

save_html(page, out_html)

cat("\n HTML saved to:\n", out_html, "\n")

#THE END
